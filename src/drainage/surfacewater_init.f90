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
!! ## L units (Phase 2 Task 9, spec D6)
!! L(:) is in centimetres when this module is called. The m→cm conversion
!! that legacy rddre performed at line 4399 is now done at TOML-read time
!! in read_drainage_toml.f90 (read_drainage_inner, per-level loop). The
!! typed config holds cm; the TOML adapter copies cm directly to the global
!! L(:) without further conversion. This module no longer mutates L.
module surfacewater_init_mod
   use iso_fortran_env, only: real64
   implicit none
   private

   public :: surfacewater_init

contains

   subroutine surfacewater_init(state)
      use swap_state_mod, only: swap_state_t
      ! SS-SWST Phase 2 Task 11: wls and swst global writes dropped; bocodre uses state aliases.
      use variables, only: nrlevs, numnod, swdtyp, zbotdr, widthr, taludr, l, &
                            wls1_init, wlp, &
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

      sw%numadj = 0

      do i = 1, 4
         sw%wlsbak(i) = 0.0_real64
      end do

      ! Initial water level pre-computed by the adapter
      ! (adapter wrote wls1_init = wlact - altcu; altcu=0 is enforced by drainage_config_validate so this equals wlact).
      ! SS-SWST Phase 2 Task 11: wls global write dropped (bocodre uses state%surfacewater%wls).
      ! wlp global kept: bocodre reads it for primary surface water level (not a SW-owned var).
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
      ! l(ilev) is in centimetres (converted at TOML-read time, spec D6),
      ! matching legacy rddre which converted l(i) = l(i)*100.0d0 before
      ! this loop. State is authoritative; global sttab write dropped
      ! (ADR 0030 Phase 2 Task 2).
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
      ! SS-SWST Phase 2 Task 11: swst global write dropped (bocodre uses state%surfacewater%swst).
      sw%swstini = swstlev(state, wls1_init)
      sw%swst    = sw%swstini

      ! Allocate per-level arrays in state cohort sub-records (guard against repeated calls).
      if (.not. allocated(sw%cumulative%cqdrain)) then
         allocate(sw%cumulative%cqdrain(nrlevs))
         sw%cumulative%cqdrain = 0.0_real64
      end if
      if (.not. allocated(sw%cumulative%cqdrainin)) then
         allocate(sw%cumulative%cqdrainin(nrlevs))
         sw%cumulative%cqdrainin = 0.0_real64
      end if
      if (.not. allocated(sw%cumulative%cqdrainout)) then
         allocate(sw%cumulative%cqdrainout(nrlevs))
         sw%cumulative%cqdrainout = 0.0_real64
      end if
      if (.not. allocated(sw%intermediate%inqdra)) then
         allocate(sw%intermediate%inqdra(nrlevs, numnod))
         sw%intermediate%inqdra = 0.0_real64
      end if
      if (.not. allocated(sw%intermediate%inqdra_in)) then
         allocate(sw%intermediate%inqdra_in(nrlevs, numnod))
         sw%intermediate%inqdra_in = 0.0_real64
      end if
      if (.not. allocated(sw%intermediate%inqdra_out)) then
         allocate(sw%intermediate%inqdra_out(nrlevs, numnod))
         sw%intermediate%inqdra_out = 0.0_real64
      end if

      end associate
   end subroutine surfacewater_init

end module surfacewater_init_mod
