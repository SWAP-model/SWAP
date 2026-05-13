!> @file hydraulic_params_mod.f90
!! SS-GR-UTILS: typed van-Genuchten / MvG / PDI hydraulic parameter
!! record. Replaces the legacy `cofgen(j, node)` magic-row matrix.
!! One instance per soil node; aggregated as
!! `state%soilwater%vg_params(:)`.
!!
!! Row-name correspondence to legacy cofgen indices
!! (from src/soil/soilhydraulics.f90:830-880 + utils/WC_K_models):
!!   cofgen(1, n)  -> vg_params(n)%thetar       (residual water content)
!!   cofgen(2, n)  -> vg_params(n)%thetas       (saturated water content)
!!   cofgen(3, n)  -> vg_params(n)%ksat
!!   cofgen(4, n)  -> vg_params(n)%alpha
!!   cofgen(5, n)  -> vg_params(n)%lpar
!!   cofgen(6, n)  -> vg_params(n)%npar
!!   cofgen(7, n)  -> vg_params(n)%mpar
!!   cofgen(8, n)  -> vg_params(n)%alphaw_sentinel
!!   cofgen(9, n)  -> vg_params(n)%h_enpr
!!   cofgen(10, n) -> vg_params(n)%ksatexm
!!   cofgen(11, n) -> vg_params(n)%relsatthr
!!   cofgen(12, n) -> vg_params(n)%ksatthr
!!   cofgen(13, n) -> vg_params(n)%alpha_2
!!   cofgen(14, n) -> vg_params(n)%npar_2
!!   cofgen(15, n) -> vg_params(n)%mpar_2
!!   cofgen(16, n) -> vg_params(n)%omega_1
!!   cofgen(17, n) -> vg_params(n)%omega_2
!!   cofgen(18, n) -> vg_params(n)%h0
!!   cofgen(19, n) -> vg_params(n)%ha
!!   cofgen(20, n) -> vg_params(n)%apar
!!   cofgen(21, n) -> vg_params(n)%omega_k
module hydraulic_params_mod
   use iso_fortran_env, only: real64
   implicit none
   private
   public :: vanGenuchten_params_t

   type :: vanGenuchten_params_t
      ! Universal MvG (rows 1-7)
      real(real64) :: thetar           = 0.0_real64
      real(real64) :: thetas           = 0.0_real64
      real(real64) :: ksat             = 0.0_real64
      real(real64) :: alpha            = 0.0_real64
      real(real64) :: lpar             = 0.0_real64
      real(real64) :: npar             = 0.0_real64
      real(real64) :: mpar             = 0.0_real64
      ! Sentinel + extreme-K parameters (rows 8-12)
      real(real64) :: alphaw_sentinel  = 0.0_real64
      real(real64) :: h_enpr           = 0.0_real64
      real(real64) :: ksatexm          = 0.0_real64
      real(real64) :: relsatthr        = 0.0_real64
      real(real64) :: ksatthr          = 0.0_real64
      ! Bi-modal MvG (rows 13-17)
      real(real64) :: alpha_2          = 0.0_real64
      real(real64) :: npar_2           = 0.0_real64
      real(real64) :: mpar_2           = 0.0_real64
      real(real64) :: omega_1          = 0.0_real64
      real(real64) :: omega_2          = 0.0_real64
      ! PDI parameters (rows 18-21)
      real(real64) :: h0               = 0.0_real64
      real(real64) :: ha               = 0.0_real64
      real(real64) :: apar             = 0.0_real64
      real(real64) :: omega_k          = 0.0_real64
   end type vanGenuchten_params_t

end module hydraulic_params_mod
