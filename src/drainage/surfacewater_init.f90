!> Runtime initialization for the surface-water management system.
!! Replaces the surviving math from the legacy `rddre` reader: build
!! the open-channel storage table sttab, set wls1/wlstar from the
!! adapter-computed wls1_init, and zero the running buffers
!! numadj/wlsbak. All inputs come from module globals already
!! populated by the TOML adapter (config_to_variables).
!!
!! Scope: swsrf=2, swsec=2, swqhr=1, swman=1, drainage.altcu=0 only.
!! Other branches are guarded with fatalerr_collected (defense in
!! depth — the config validator rejects them upstream too).
!!
!! ## L unit reconciliation (dramet=0 / swdra=2 only)
!! This module is only reachable for swdra=2 cases. The cross-field
!! validator in drainage_config_validate rejects swdra=2 + dramet/=0,
!! so at runtime dramet is guaranteed to be 0.
!!
!! Legacy rddre reads L from the .dra file in METERS (valid range
!! 1..100000 m), then immediately converts to CENTIMETRES at line 4399:
!!
!!     l(i) = l(i)*100.0d0
!!
!! before using l(ilev) in the sttab volume computation and before
!! drainage.f90 line 661 uses l(level) in the open-channel resistance
!! formula.  The TOML adapter pre-converts L to cm ONLY for
!! dramet=3 + swdivd=1 cases (config_to_variables.f90:404-408).  For
!! the dramet=0 / swdra=2 path the adapter does NOT pre-convert, so
!! the global L(:) array is still in metres when this module is called.
!! This module therefore applies the *100 conversion to the global L(:)
!! array (mirroring rddre) so that both the sttab geometry math and all
!! downstream resistance calculations (drainage.f90:661) operate in
!! centimetres.
!!
!! If a future case authors swdra=2 + dramet=3, the validator rejects
!! it before reaching this module, preventing a double-conversion of L.
module surfacewater_init_mod
   use iso_fortran_env, only: real64
   implicit none
   private

   public :: surfacewater_init

contains

   subroutine surfacewater_init(state)
      use swap_state_mod, only: swap_state_t
      use variables, only: nrlevs, numnod, swdtyp, zbotdr, widthr, taludr, l, &
                            wls1_init, wls, wlp, &
                            swst, &
                            swsrf, swsec, swqhr, swman, nmper
      use surfacewater_utils, only: swstlev
      use error_mod, only: fatalerr_collected
      type(swap_state_t), intent(inout) :: state

      integer      :: i, ilev
      real(real64) :: wdepth, wvolum, wbreadth
      integer      :: nrpri

      associate(sw => state%surfacewater)

      ! Defensive guards mirroring surface_water_config_validate.
      ! swman is a fixed-size array (dimensioned mamp); slice 1:nmper
      ! to compare authored periods only.
      if (swsrf == 3 .or. swsec == 1 .or. swqhr == 2) then
         call fatalerr_collected('surfacewater_init', &
            'swsrf=3, swsec=1, or swqhr=2 not supported on the TOML path')
         return
      end if
      if (any(swman(1:nmper) == 2)) then
         call fatalerr_collected('surfacewater_init', &
            'swman=2 (automatic weir) not supported on the TOML path')
         return
      end if

      ! For swsrf=2 (no primary system) nrpri = 0.
      nrpri = 0

      ! Convert L from metres -> centimetres (mirrors rddre line 4399:
      ! l(i) = l(i)*100.0d0). The adapter writes L in metres for the
      ! dramet=0/swdra=2 path; downstream drainage.f90 (line 661) and
      ! the sttab geometry math both expect L in centimetres.
      do i = 1, nrlevs
         l(i) = l(i) * 100.0_real64
      end do

      ! numadj and wlsbak written to state only; legacy globals are dead.
      sw%numadj = 0

      do i = 1, 4
         sw%wlsbak(i) = 0.0_real64
      end do

      ! Initial water level pre-computed by the adapter
      ! (adapter wrote wls1_init = wlact - altcu; altcu=0 is enforced by drainage_config_validate so this equals wlact).
      ! wls global kept: drainage.f90 (bocodre) reads it for previous-timestep surface water level.
      ! wlp global kept: drainage.f90 (bocodre) reads it for primary surface water level.
      ! wlstar global dropped: only output reads it, via state%surfacewater%wlstar.
      wls    = wls1_init
      wlp    = 0.0_real64    ! swsrf=2 has no primary system

      sw%wls    = wls1_init
      sw%wlstar = wls1_init

      ! sttab(:,1) — depths. Row 1 = +100cm above soil surface;
      ! row 2 = 0cm (soil surface); rows 3..22 divide
      ! [0, zbotdr(1+nrpri)] into 20 compartments.
      ! State is now authoritative; global sttab write dropped (ADR 0030 Phase 2 Task 2).
      sw%sttab(1, 1) = 100.0_real64
      sw%sttab(2, 1) =   0.0_real64
      do i = 3, 22
         sw%sttab(i, 1) = zbotdr(1 + nrpri) * (i - 2) / 20.0_real64
      end do

      ! sttab(:,2) — storage volume per unit area (cm), summed across
      ! open-channel levels (swdtyp=0). Verbatim port from
      ! readswap.f90:4878-4897.
      ! l(ilev) is now in centimetres (converted above), matching
      ! legacy rddre which converts l(i) = l(i)*100.0d0 before this loop.
      ! State is authoritative; global sttab write dropped (ADR 0030 Phase 2 Task 2).
      do i = 1, 22
         sw%sttab(i, 2) = 0.0_real64
         do ilev = 1 + nrpri, nrlevs
            if (swdtyp(ilev) == 0 .and. sw%sttab(i, 1) > zbotdr(ilev)) then
               if (sw%sttab(i, 1) <= 0.0_real64) then
                  ! Trapezium below soil surface
                  wdepth = sw%sttab(i, 1) - zbotdr(ilev)
                  wvolum = wdepth * (widthr(ilev) + wdepth / taludr(ilev))
               else
                  ! Trapezium up to surface, plus rectangle above
                  wdepth   = -zbotdr(ilev)
                  wvolum   = wdepth * (widthr(ilev) + wdepth / taludr(ilev))
                  wbreadth = widthr(ilev) + 2.0_real64 * wdepth / taludr(ilev)
                  wdepth   = sw%sttab(i, 1)
                  wvolum   = wvolum + wbreadth * wdepth
               end if
               sw%sttab(i, 2) = sw%sttab(i, 2) + wvolum / l(ilev)
            end if
         end do
      end do

      ! Initial storage state.
      ! swstini global dropped: only output reads it, via state%surfacewater%swstini.
      ! swst global kept: drainage.f90 (bocodre) reads it for previous-timestep storage.
      sw%swstini = swstlev(state, wls1_init)
      swst       = sw%swstini
      sw%swst    = sw%swstini

      ! Allocate per-level arrays in state (guard against repeated calls).
      if (.not. allocated(sw%cqdrain)) then
         allocate(sw%cqdrain(nrlevs))
         sw%cqdrain = 0.0_real64
      end if
      if (.not. allocated(sw%cqdrainin)) then
         allocate(sw%cqdrainin(nrlevs))
         sw%cqdrainin = 0.0_real64
      end if
      if (.not. allocated(sw%cqdrainout)) then
         allocate(sw%cqdrainout(nrlevs))
         sw%cqdrainout = 0.0_real64
      end if
      if (.not. allocated(sw%qdra)) then
         allocate(sw%qdra(nrlevs, numnod))
         sw%qdra = 0.0_real64
      end if
      if (.not. allocated(sw%inqdra)) then
         allocate(sw%inqdra(nrlevs, numnod))
         sw%inqdra = 0.0_real64
      end if
      if (.not. allocated(sw%inqdra_in)) then
         allocate(sw%inqdra_in(nrlevs, numnod))
         sw%inqdra_in = 0.0_real64
      end if
      if (.not. allocated(sw%inqdra_out)) then
         allocate(sw%inqdra_out(nrlevs, numnod))
         sw%inqdra_out = 0.0_real64
      end if

      end associate
   end subroutine surfacewater_init

end module surfacewater_init_mod
