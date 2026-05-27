!> Pure solute physics kernels — explicit args, no swap_state_t dependency.
!> Pilot extraction from solute_step (see spec 2026-05-27). Mirrors the
!> interception.f90 / et.f90 "pure function — explicit args" pattern.
module solute_kernels_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: bdenskf_coeff, bdenskfsatporos_coeff, ddiffwcs_coeff, decpotfdepth_coeff
   public :: solute_cml_from_cmsy
   public :: solute_ftemp, solute_ftheta, solute_decomp_ctrans

contains

   !> Freundlich bulk sorption coefficient: bdens * kf  [-].
   elemental function bdenskf_coeff(bdens, kf) result(v)
      real(real64), intent(in) :: bdens   ! dry soil bulk density (g/cm3)
      real(real64), intent(in) :: kf      ! Freundlich coefficient (cm3/g)
      real(real64)             :: v
      v = bdens*kf
   end function bdenskf_coeff

   !> Saturated-zone sorption + porosity: bdens * kfsat + poros  [-].
   elemental function bdenskfsatporos_coeff(bdens, kfsat, poros) result(v)
      real(real64), intent(in) :: bdens   ! dry soil bulk density (g/cm3)
      real(real64), intent(in) :: kfsat   ! saturated-zone Freundlich coeff (cm3/g)
      real(real64), intent(in) :: poros   ! aquifer porosity (-)
      real(real64)             :: v
      v = bdens*kfsat + poros
   end function bdenskfsatporos_coeff

   !> Tortuosity-scaled diffusion base: ddif / thetsl**2  (cm2/d).
   !> Precondition: thetsl > 0 (saturated water content; a positive soil property).
   elemental function ddiffwcs_coeff(ddif, thetsl) result(v)
      real(real64), intent(in) :: ddif    ! molecular diffusion coefficient (cm2/d)
      real(real64), intent(in) :: thetsl  ! saturated water content (-)
      real(real64)             :: v
      v = ddif / (thetsl**2)
   end function ddiffwcs_coeff

   !> Depth-weighted potential decomposition: decpot * fdepth  (1/d).
   elemental function decpotfdepth_coeff(decpot, fdepth) result(v)
      real(real64), intent(in) :: decpot  ! potential decomposition rate (1/d)
      real(real64), intent(in) :: fdepth  ! depth-decomposition factor (-)
      real(real64)             :: v
      v = decpot*fdepth
   end function decpotfdepth_coeff

   !> Recover mobile concentration cml from total cmsy via the Freundlich
   !> isotherm. Linear shortcut when frexp ~ 1; otherwise fixed-point iterate
   !> seeded from cml_guess. Caller handles the cmsy < vsmall zeroing.
   pure function solute_cml_from_cmsy(cmsy, theta, bdenskf, frexp, cref, cml_guess) result(cml)
      real(real64), intent(in) :: cmsy       ! total (dissolved+adsorbed) conc (M/L3 soil)
      real(real64), intent(in) :: theta      ! volumetric water content (-)
      real(real64), intent(in) :: bdenskf    ! bdens*kf (-)
      real(real64), intent(in) :: frexp      ! Freundlich exponent (-)
      real(real64), intent(in) :: cref       ! reference concentration (M/L3)
      real(real64), intent(in) :: cml_guess  ! previous cml, iteration seed (M/L3)
      real(real64)             :: cml

      real(real64), parameter :: rer    = 1.0d-3
      real(real64), parameter :: vsmall = 1.0d-15
      real(real64) :: old, dummy
      logical      :: differ

      if (abs(frexp - 1.0d0) .lt. 0.001d0) then
         cml = cmsy / (theta + bdenskf)
      else
         cml = cml_guess
         if (cml .lt. vsmall) cml = vsmall
         differ = .true.
         do while (differ)
            old   = cml
            dummy = bdenskf*(cml/cref)**(frexp - 1.0d0)
            cml   = cmsy/(theta + dummy)
            if (abs(cml - old) .lt. rer*cml) differ = .false.
         end do
      end if
   end function solute_cml_from_cmsy

   !> Temperature reduction factor for decomposition. Capped above 35 degC;
   !> zero when the temperature switch is off.
   elemental function solute_ftemp(tsoil, gampar, fl_temperature) result(ftemp)
      real(real64), intent(in) :: tsoil           ! soil temperature (degC)
      real(real64), intent(in) :: gampar          ! temperature coefficient (/C)
      logical,      intent(in) :: fl_temperature  ! temperature simulation on/off
      real(real64)             :: ftemp
      if (fl_temperature) then
         if (tsoil .lt. 35.0d0) then
            ftemp = exp(gampar*(tsoil - 20.0d0))
         else
            ftemp = exp(gampar*15.0d0)
         end if
      else
         ftemp = 0.0d0
      end if
   end function solute_ftemp

   !> Moisture reduction factor for decomposition, clamped to 1.
   !> Precondition: rtheta > 0 (reference moisture content).
   elemental function solute_ftheta(theta, rtheta, bexp) result(ftheta)
      real(real64), intent(in) :: theta   ! volumetric water content (-)
      real(real64), intent(in) :: rtheta  ! reference moisture content (-)
      real(real64), intent(in) :: bexp    ! moisture-decomposition exponent (-)
      real(real64)             :: ftheta
      ftheta = min(1.0d0, (theta/rtheta)**bexp)
   end function solute_ftheta

   !> Solute transformation (decomposition) rate per node.
   elemental function solute_decomp_ctrans(decact, theta, cml, bdenskfcref, cref, frexp) result(ctrans)
      real(real64), intent(in) :: decact       ! actual decomposition rate (1/d)
      real(real64), intent(in) :: theta        ! volumetric water content (-)
      real(real64), intent(in) :: cml          ! mobile concentration (M/L3)
      real(real64), intent(in) :: bdenskfcref  ! bdens*kf*cref (M/L3)
      real(real64), intent(in) :: cref         ! reference concentration (M/L3)
      real(real64), intent(in) :: frexp        ! Freundlich exponent (-)
      real(real64)             :: ctrans
      ctrans = decact*theta*cml + decact*bdenskfcref*((cml/cref)**frexp)
   end function solute_decomp_ctrans

end module solute_kernels_mod
