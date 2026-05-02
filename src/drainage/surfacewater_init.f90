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
!! ## L unit reconciliation
!! Legacy rddre reads L from the .dra file in METERS (valid range
!! 1..100000 m), then immediately converts to CENTIMETRES at line 4399:
!!
!!     l(i) = l(i)*100.0d0
!!
!! before using l(ilev) in the sttab volume computation.  The TOML
!! adapter currently only multiplies L by 100 for dramet=3 cases
!! (config_to_variables.f90:404-408), so when this module is called
!! for the dramet=0 / swdra=2 path, the global L(:) array is still
!! in metres.  The multiplication is therefore applied locally inside
!! the sttab loop below so that the global stays consistent with what
!! the adapter writes (matching what downstream drainage.f90 expects
!! for resistance calculations), while the storage geometry math gets
!! the correct centimetre value.
module surfacewater_init_mod
   use iso_fortran_env, only: real64
   implicit none
   private

   public :: surfacewater_init

contains

   subroutine surfacewater_init(wls1, wlp1)
      use variables, only: nrlevs, swdtyp, zbotdr, widthr, taludr, l, &
                            wls1_init, wlstar, &
                            sttab, swstini, swst, wlsbak, numadj, &
                            swsrf, swsec, swqhr, swman, nmper
      use surfacewater_utils, only: swstlev
      use error_mod, only: fatalerr_collected
      real(real64), intent(out) :: wls1, wlp1

      integer      :: i, ilev
      real(real64) :: wdepth, wvolum, wbreadth
      integer      :: nrpri

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

      numadj = 0
      do i = 1, 4
         wlsbak(i) = 0.0_real64
      end do

      ! Initial water level pre-computed by the adapter
      ! (= wlact - altcu; altcu enforced = 0 by Task 7).
      wls1   = wls1_init
      wlp1   = 0.0_real64    ! swsrf=2 has no primary system
      wlstar = wls1

      ! sttab(:,1) — depths. Row 1 = +100cm above soil surface;
      ! row 2 = 0cm (soil surface); rows 3..22 divide
      ! [0, zbotdr(1+nrpri)] into 20 compartments.
      sttab(1, 1) = 100.0_real64
      sttab(2, 1) =   0.0_real64
      do i = 3, 22
         sttab(i, 1) = zbotdr(1 + nrpri) * (i - 2) / 20.0_real64
      end do

      ! sttab(:,2) — storage volume per unit area (cm), summed across
      ! open-channel levels (swdtyp=0). Verbatim port from
      ! readswap.f90:4878-4897.
      !
      ! L unit reconciliation: global l(ilev) is in metres for the
      ! dramet=0/swdra=2 path (the adapter does not multiply by 100
      ! for this case).  Legacy rddre converts l(i) = l(i)*100 before
      ! the sttab math, so we apply the factor locally here.  The
      ! global L(:) remains in metres for downstream resistance calcs.
      do i = 1, 22
         sttab(i, 2) = 0.0_real64
         do ilev = 1 + nrpri, nrlevs
            if (swdtyp(ilev) == 0 .and. sttab(i, 1) > zbotdr(ilev)) then
               if (sttab(i, 1) <= 0.0_real64) then
                  ! Trapezium below soil surface
                  wdepth = sttab(i, 1) - zbotdr(ilev)
                  wvolum = wdepth * (widthr(ilev) + wdepth / taludr(ilev))
               else
                  ! Trapezium up to surface, plus rectangle above
                  wdepth   = -zbotdr(ilev)
                  wvolum   = wdepth * (widthr(ilev) + wdepth / taludr(ilev))
                  wbreadth = widthr(ilev) + 2.0_real64 * wdepth / taludr(ilev)
                  wdepth   = sttab(i, 1)
                  wvolum   = wvolum + wbreadth * wdepth
               end if
               ! l(ilev) * 100 converts metres -> centimetres to match
               ! legacy rddre line 4399: l(i) = l(i)*100.0d0
               sttab(i, 2) = sttab(i, 2) + wvolum / (l(ilev) * 100.0_real64)
            end if
         end do
      end do

      ! Initial storage state.
      swstini = swstlev(wls1)
      swst    = swstini
   end subroutine surfacewater_init

end module surfacewater_init_mod
