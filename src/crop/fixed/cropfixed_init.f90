!> Runtime initialization for the type-1 (cropfixed) rotation, replacing
!! the per-rotation runtime side-effects of legacy `readcropfixed`.
!!
!! Operates on `cropfixed_config_t` (already validated and populated at
!! config-load time) and writes to the same `variables`-module globals
!! that legacy `readcropfixed` writes to. After the copy, runs the
!! tail-of-readcropfixed init math (cumdens build, swrdc default).
!!
!! Note: `lcc` is NOT a `variables` module global — it is a local SAVE
!! variable in `cropgrowth.f90` that is passed as an argument to the
!! legacy `readcropfixed`. This subroutine therefore returns it via the
!! `lcc` intent(out) argument so that the caller (cropgrowth task=1) can
!! assign it to its own local.
!!
!! Scope (Phase 1, ADR 0015): only the supported subset of switches.
!! Defense-in-depth runtime guards mirror the validator stub-errors.
module cropfixed_init_mod
   use iso_fortran_env, only: real64
   use cropfixed_config_mod, only: cropfixed_config_t
   implicit none
   private

   public :: cropfixed_init_from_config

contains

   subroutine cropfixed_init_from_config(cfg, icrop, lcc, state)
      ! [GR-CROP 2026-05-25] cropfixed_init writes max_resp_factor directly to
      ! state%crop%oxygen%max_resp_factor (the only reader); legacy global retired.
      use array_utils,  only: afgen
      use error_mod,    only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      type(cropfixed_config_t), intent(in)    :: cfg
      integer,                  intent(in)    :: icrop  ! reserved for swinco=3 inifil path (not yet ported)
      integer,                  intent(out)   :: lcc
      type(swap_state_t),       intent(inout) :: state  ! state writes for crop runtime init

      integer      :: i
      real(real64) :: depth, rootdis(202), sum

      ! ---- Defense-in-depth: stub-error gates mirror cropfixed_config_validate.
      if (cfg%swdrought == 2 .or. cfg%swoxygen == 2 .or.                    &
          cfg%swrd == 3     .or. cfg%swsalinity /= 0 .or. cfg%schedule_switch == 1) then
         call fatalerr_collected('cropfixed_init_from_config', &
            'Unsupported runtime branch reached on the TOML path. The validator ' // &
            'should have caught this earlier.')
         return
      end if

      ! ---- Copy scalars / switches to globals ------------------------
      state%crop%common%idev   = cfg%idev
      lcc    = cfg%lcc
      state%crop%common%tsumea = cfg%tsumea
      state%crop%common%tsumam = cfg%tsumam
      state%crop%common%tbase  = cfg%tbase

      state%crop%kdif = cfg%kdif
      state%crop%kdir = cfg%kdir

      state%crop%common%swgc = cfg%swgc
      state%crop%swcf = cfg%swcf
      state%crop%common%swrd     = cfg%swrd
      state%crop%common%swdmi2rd = cfg%swdmi2rd
      state%crop%common%swrdc    = cfg%swrdc
      state%crop%common%rdi = cfg%rdi
      state%crop%common%rri = cfg%rri
      state%crop%common%rdc = cfg%rdc

      state%crop%common%swoxygen = cfg%swoxygen
      state%crop%common%swWrtNonox = cfg%swwrtnonox
      state%crop%common%aeratecrit = cfg%aeratecrit
      state%crop%oxygen%max_resp_factor = cfg%max_resp_factor
      state%crop%common%hlim1  = cfg%hlim1
      state%crop%common%hlim2u = cfg%hlim2u
      state%crop%common%hlim2l = cfg%hlim2l

      state%crop%common%swdrought = cfg%swdrought
      state%crop%common%hlim3h = cfg%hlim3h
      state%crop%common%hlim3l = cfg%hlim3l
      state%crop%common%hlim4  = cfg%hlim4
      state%crop%common%adcrh  = cfg%adcrh
      state%crop%common%adcrl  = cfg%adcrl

      state%crop%common%swsalinity = cfg%swsalinity
      state%crop%common%saltmax   = cfg%saltmax
      state%crop%common%saltslope = cfg%saltslope
      state%crop%common%salthead  = cfg%salthead

      state%crop%common%swcompensate = cfg%swcompensate
      state%crop%common%swstressor = cfg%swstressor
      state%crop%common%alphacrit  = cfg%alphacrit
      state%crop%common%dcritrtz   = cfg%dcritrtz

      state%crop%common%swinter = cfg%swinter
      state%crop%cofab = cfg%cofab
      ! Adapted-Rutter storage interception (swinter=3): fimin lives on the
      ! atmosphere state; siccaplai feeds siccapact = siccaplai*lai in the runtime.
      state%atmosphere%fimin       = cfg%fimin
      state%crop%common%siccaplai  = cfg%siccaplai

      state%crop%common%schedule = cfg%schedule_switch

      state%crop%common%dvsend = cfg%dvsend
      state%crop%common%swharv = cfg%swharv

      ! ---- Reflection coefficients / crop resistance (legacy:2150-2160)
      ! Two branches in legacy:
      !   swcf=1 / swcf=3 (ETref crop factor / wet-crop factor): hardcoded
      !       ETref defaults (albedo=0.23, rsc=70, rsw=0)
      !   swcf=2 (crop height): user-authored values from .crp
      if (cfg%swcf == 1 .or. cfg%swcf == 3) then
         state%crop%common%albedo = 0.23_real64
         state%crop%common%rsc    = 70.0_real64
         state%crop%common%rsw    = 0.0_real64
      else if (cfg%swcf == 2) then
         state%crop%common%albedo = cfg%albedo
         state%crop%common%rsc    = cfg%rsc
         state%crop%common%rsw    = cfg%rsw
      end if

      ! ---- Copy tables ------------------------------------------------
      ! gctb (size up to 2*magrs in legacy; we copy what was authored).
      if (allocated(cfg%gctb))  call copy_pair_table(cfg%gctb, state%crop%fixed%gctb)
      if (allocated(cfg%cftb))    call copy_pair_table(cfg%cftb,    state%crop%fixed%cftb)
      if (allocated(cfg%chtb))    call copy_pair_table(cfg%chtb,    state%crop%fixed%chtb)
      if (allocated(cfg%cfeictb)) call copy_pair_table(cfg%cfeictb, state%crop%fixed%cfeictb)
      if (allocated(cfg%rdtb))  call copy_pair_table(cfg%rdtb, state%crop%common%rdtb)

      ! Gash forest-interception tables (swinter=2) live on the atmosphere
      ! state; legacy readcropfixed wrote the matching globals. Mirror that.
      if (allocated(cfg%pfreetb))   call copy_pair_table(cfg%pfreetb,   state%atmosphere%pfreetb)
      if (allocated(cfg%pstemtb))   call copy_pair_table(cfg%pstemtb,   state%atmosphere%pstemtb)
      if (allocated(cfg%scanopytb)) call copy_pair_table(cfg%scanopytb, state%atmosphere%scanopytb)
      if (allocated(cfg%avprectb))  call copy_pair_table(cfg%avprectb,  state%atmosphere%avprectb)
      if (allocated(cfg%avevaptb))  call copy_pair_table(cfg%avevaptb,  state%atmosphere%avevaptb)

      ! rdctb is sized 22 in legacy; we always copy.
      if (allocated(cfg%rdctb)) then
         do i = 1, min(size(cfg%rdctb), size(state%crop%common%rdctb))
            state%crop%common%rdctb(i) = cfg%rdctb(i)
         end do
      end if

      ! ---- Tail-of-readcropfixed RUNTIME init: cumdens (legacy:2449-2478)
      if (cfg%swdrought == 1) then
         do i = 0, 100
            depth = 0.01_real64 * real(i, real64)
            rootdis(i*2 + 1) = depth
            rootdis(i*2 + 2) = afgen(state%crop%common%rdctb, 22, depth)
         end do
         do i = 1, 202, 2
            state%crop%common%cumdens(i) = rootdis(i)
         end do
         sum = 0.0_real64
         state%crop%common%cumdens(2) = 0.0_real64
         do i = 4, 202, 2
            sum = sum + (rootdis(i-2) + rootdis(i)) * 0.5_real64 &
                      * (state%crop%common%cumdens(i-1) - state%crop%common%cumdens(i-3))
            state%crop%common%cumdens(i) = sum
         end do
         do i = 2, 202, 2
            state%crop%common%cumdens(i) = state%crop%common%cumdens(i) / sum
         end do
      end if
   end subroutine cropfixed_init_from_config

   ! Copy a flat (dvs, value) pair source into a destination global of the
   ! same flat-array shape. Destination is a fixed-size global; copy up to
   ! min(size).
   subroutine copy_pair_table(src, dst)
      real(real64), intent(in)    :: src(:)
      real(real64), intent(inout) :: dst(:)
      integer :: i
      do i = 1, min(size(src), size(dst))
         dst(i) = src(i)
      end do
   end subroutine copy_pair_table

end module cropfixed_init_mod
