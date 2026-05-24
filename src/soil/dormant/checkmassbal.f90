!> @file src/soil/dormant/checkmassbal.f90
!! @brief DORMANT — per-period mass-balance audit for ANIMO/PEARL coupling.
!!
!! ## Status: DORMANT
!!
!! Extracted from `src/soil/waterbalance.f90` on 2026-05-24. Has had **no
!! callers** in the TOML-pipeline build since the migration: the legacy
!! dispatch site (called once per output period, gated by an
!! `output_csv.dwb` style switch when ANIMO/PEARL coupling was enabled)
!! was not ported. The body is preserved verbatim for future
!! reactivation because the ANIMO/PEARL coupling target is on the
!! long-term roadmap.
!!
!! Excluded from `meson.build` — this file is NOT compiled.
!!
!! ## What it does
!!
!! Walks the soil profile + ponding layer + (historical) macropore
!! domains at the end of every output period, computing the
!! mass-balance closure deviation for each subsystem. When the
!! deviation exceeds `CritDevMasBal`, the offending sub-balance is
!! written to `<outfil>.dwb.csv` for offline diagnosis. The output
!! file is the standard mass-balance audit channel for the
!! ANIMO/PEARL coupling — independent verification that SWAP's
!! water fluxes close against the same accounting the solute model
!! expects.
!!
!! ## Reactivation checklist
!!
!! 1. **Config gate** — add `output_csv.dwb : boolean` (or similar)
!!    to the TOML schema and surface it as `config%output%dwb` or
!!    similar. Dispatch this routine only when the gate is true.
!!
!! 2. **Globals to migrate** (currently still in `variables.f90`):
!!    - `NumNodNew`, `DZNew` — regridding temporaries (also needed
!!      by `regrid.f90`; see `src/soil/dormant/regrid.f90`).
!!      Migrate to `state%mesh` as regrid-temporary fields.
!!    - `outfil`, `pathwork` — output path globals. Either migrate
!!      to `state%cfg%output` (Pattern 2 snapshot) or pass through
!!      as subroutine arguments.
!!    - `CritDevMasBal` — mass-balance criterion. Migrate to
!!      `state%cfg%output` or `state%cfg%soil`.
!!    - `dev_cmb` — file-handle integer. Move to `state%output`
!!      (or a local file-handle helper) — file handles do not
!!      belong in `variables.f90`.
!!
!! 3. **Macropore branch was retired** (ADR 0040). The historical
!!    routine had Section 4 ("Macropore domains Dm1/Dm2") plus
!!    accumulation of `IQExcMtxDm1/2`, `IQInTopPreDm`,
!!    `IQInTopLatDm` from `qimmob`, `QExcMpMtx`, `QMaPo`, `qssdi`,
!!    `qssdisum`. Those globals were retired with the macropore
!!    arc. The dormant body below keeps the structure (variables
!!    left as locals initialised to zero) so reactivation does not
!!    need to re-add the macropore code paths — the per-domain
!!    write statements have been dropped from format labels 1-7
!!    that remain.
!!
!! 4. **Dispatch site** — legacy SWAP called `checkmassbal` once
!!    per output period from the output-coupling section in
!!    `swap_mod.f90`. Restore the call with the new gate.
!!
!! 5. **Add this file to `meson.build`** sources list, then run
!!    `pixi run check-fast` — any references to globals migrated
!!    in the meantime will surface as missing-symbol errors.
!!
!! ## Original docstring
!!
!! Date: 26-jun-2003
!!
!! "Purpose: Checking of mass balance per period OutPer for
!! ANIMO/PEARL output. File usage: outfil. SAVE removed - dev_cmb
!! now in variables.f90 module."
module checkmassbal_dormant_mod
   use swap_state_mod, only: swap_state_t
   use swap_array_dimensions, only: macp, madr
   use file_io_mod, only: file_open
   use variables, only: NumNodNew, outfil, pathwork, DZNew,    &
                        CritDevMasBal, dev_cmb
   implicit none
   private
   public :: checkmassbal

contains

   !> Per-period mass-balance audit.
   !!
   !! @param[inout] flopenfiledev      file-open flag (caller persists across calls)
   !! @param[in]    inqdranew          drainage flux per level/compartment
   !! @param[in]    iqexcmtxdm1cpnew   matrix->macropore exchange dm1 (always zero post-ADR-0040)
   !! @param[in]    iqexcmtxdm2cpnew   matrix->macropore exchange dm2 (always zero post-ADR-0040)
   !! @param[in]    inqnew             inter-compartment fluxes
   !! @param[in]    iqoutdrrapcpnew    outflow to drains/rapid drainage per compartment
   !! @param[in]    inqrotnew          root extraction per compartment
   !! @param[in]    ithetabegnew       initial water content per compartment
   !! @param[in]    thetanew           current water content per compartment
   !! @param[in]    state              swap state (sub-record associate)
   subroutine checkmassbal (flopenfiledev,inqdranew,iqexcmtxdm1cpnew,iqexcmtxdm2cpnew,inqnew,iqoutdrrapcpnew,inqrotnew,ithetabegnew,thetanew,state)
      implicit none

      real(8) IQExcMtxDm1CpNew(macp), IQExcMtxDm2CpNew(macp)
      real(8) inqdraNew(Madr,macp)
      real(8) inqNew(macp+1), IQOutDrRapCpNew(macp), inqrotNew(macp)
      real(8) IThetaBegNew(MaCp),thetaNew(macp)
      logical FlOpenFileDev
      type(swap_state_t), intent(in) :: state

      integer Level, ic
      real(8) DevMasBalDm1,DevMasBalDm2, DevMasBalCmp(MaCp)
      real(8) DevMasBalPnd, DevMasBalPrf, IQExcMtxDm1
      real(8) IQExcMtxDm2,IQInTopLatDm,  IQInTopPreDm, IQOutDrRap
      real(8) Qdra(MaCp), QdraPrf, QrotPrf, SrDif
      real(8) WaSr(MaCp), WaSrBeg(MaCp), WaSrPrf, WaSrPrfBeg
      character(len=300) filnam
      logical FlWriteDevCmp(MaCp), FlWriteDev, FlWriteDevDm1
      logical FlWriteDevDm2, FlWriteDevPnd, FlWriteDevPrf

      associate (soil => state%soilwater,    &
                 atmo => state%atmosphere,   &
                 drai => state%drainage,     &
                 time => state%timecontrol)

      FlWriteDev    = .false.
      FlWriteDevPnd = .false.
      FlWriteDevPrf = .false.
      do ic = 1, NumNodNew
         FlWriteDevCmp(ic) = .false.
      end do
      FlWriteDevDm1 = .false.
      FlWriteDevDm2 = .false.

      ! 1) Ponding layer
      SrDif        = soil%IPondBeg - soil%pond + atmo%ISsnowBeg - atmo%ssnow
      IQInTopPreDm = 0.d0
      IQInTopLatDm = 0.d0

      DevMasBalPnd = atmo%intr%igrai + atmo%intr%igsnow + soil%igird + soil%irunon + inqNew(1) + SrDif &
                   - (atmo%intr%igrai - atmo%intr%inrai - atmo%intr%isnrai + soil%igird - soil%inird &
                       + atmo%intr%isubl + atmo%intr%ievap + soil%iruno) &
                   - IQInTopPreDm - IQInTopLatDm

      if (abs(DevMasBalPnd) .gt. CritDevMasBal) then
         FlWriteDev    = .true.
         FlWriteDevPnd = .true.
      end if

      ! 2) Total Soil Profile
      WaSrPrfBeg  = 0.d0
      WaSrPrf     = 0.d0
      QrotPrf     = 0.d0
      QdraPrf     = 0.d0
      IQExcMtxDm1 = 0.d0
      IQExcMtxDm2 = 0.d0
      do ic = 1, numnodnew
         WaSrPrfBeg = WaSrPrfBeg + dzNew(ic) * IThetaBegNew(ic)
         WaSrPrf    = WaSrPrf    + dzNew(ic) * ThetaNew(ic)
         QrotPrf    = QrotPrf    + inqrotNew(ic)
         do level = 1, drai%nrlevs
            QdraPrf = QdraPrf + InqdraNew(level, ic)
         end do
      end do
      SrDif = WaSrPrfBeg - WaSrPrf

      DevMasBalPrf = inqNew(NumNodNew+1) + SrDif + IQExcMtxDm1 + IQExcMtxDm2 &
                   - (inqNew(1) + QrotPrf + QdraPrf)

      if (abs(DevMasBalPrf) .gt. CritDevMasBal) then
         FlWriteDev    = .true.
         FlWriteDevPrf = .true.
      end if

      ! 3) Individual Soil Compartments
      do ic = 1, numnodnew
         SrDif      = 0.d0
         Qdra(ic)   = 0.d0
         WaSrBeg(ic) = dzNew(ic) * IThetaBegNew(ic)
         WaSr(ic)    = dzNew(ic) * ThetaNew(ic)
         SrDif       = WaSrBeg(ic) - WaSr(ic)
         do level = 1, drai%nrlevs
            Qdra(ic) = Qdra(ic) + inqdraNew(level, ic)
         end do

         DevMasBalCmp(ic) = inqNew(ic+1) + SrDif &
                          - (inqNew(ic) + inqrotNew(ic) + Qdra(ic))

         if (abs(DevMasBalCmp(ic)) .gt. CritDevMasBal) then
            FlWriteDev        = .true.
            FlWriteDevCmp(ic) = .true.
         end if
      end do

      ! 4) Macropore domains Dm1/Dm2 — dropped (ADR 0040; FlMacropore always .false.)

      ! Open output file on first deviation
      if (FlWriteDev .and. .not. FlOpenFileDev) then
         filnam = trim(pathwork) // trim(outfil) // '.dwb'
         call file_open(dev_cmb, filnam, 'replace', 'write')
         write(dev_cmb, 1)
         FlOpenFileDev = .true.
      end if

      if (FlWriteDevPnd) write(dev_cmb, 3) time%daycum, DevMasBalPnd, &
         atmo%intr%igrai, atmo%intr%igsnow, soil%igird, soil%irunon, atmo%intr%isnrai, &
         atmo%intr%igrai - atmo%intr%inrai, soil%igird - soil%inird, &
         atmo%intr%isubl, atmo%intr%ievap, soil%iruno, inqNew(1), soil%pond, soil%IPondBeg, atmo%ssnow, &
         atmo%ISsnowBeg, IQInTopPreDm, IQInTopLatDm

      if (FlWriteDevPrf) write(dev_cmb, 4) time%daycum, DevMasBalPrf, &
         inqNew(1), inqNew(NumNodNew+1), QrotPrf, QdraPrf, WaSrPrf, &
         WaSrPrfBeg, IQExcMtxDm1, IQExcMtxDm2

      do ic = 1, numnodnew
         if (FlWriteDevCmp(ic)) write(dev_cmb, 5) time%daycum, ic, DevMasBalCmp(ic), &
            inqNew(ic), inqNew(ic+1), inqrotNew(ic), Qdra(ic), WaSr(ic), &
            WaSrBeg(ic), IQExcMtxDm1CpNew(ic), IQExcMtxDm2CpNew(ic)
      end do

 1    format(' DEVIATIONS WATERBALANCE for different subsystems: 1. Pon'&
           &'d.layer; 2. Whole profile; 3. Compartment',/,                  &
           &' Relevant terms of waterbalance per subsystem',                &
           &' (all terms in cm):',//,                                       &
           &' DayCum, 1. PONDLAY., DevMasBalAbs, IgRai, IgSnow, IgIrd, IRunon'&
           &', SnowFall,IntcpRai, IntcpIrd, ISubl, IEvap, IRuno, InQTop,   ',&
           &'Pond, IPondBeg, Ssnow, ISsnowBeg, IQInTopPreDm, IQInTopLatDm,',/,&
           &' , 2. PROFILE, DevMasBalPrf, InQTop, InQBot, QrotPrf, QdraPrf,',&
           &' WaSrPrf, WaSrPrfBeg, InQExcMtxDm1, InQExcMtxDm2',/,           &
           &' , 3. COMPno, DevMasBalCmp, InQNew(top), InQNew(bot),',        &
           &' InQrotNew, Qdra, WaSr, WaSrBeg, InQExcMtxDm1CpNew,',          &
           &' InQExcMtxDm2CpNew')
 3    format(i5,',',' Pondlay. : ',18(',',f12.8))
 4    format(i5,',',' Profile : ',9(',',f12.8))
 5    format(i5,',',' Comp',i3,': ',9(',',f12.8))

      end associate

      return
   end subroutine checkmassbal

end module checkmassbal_dormant_mod
