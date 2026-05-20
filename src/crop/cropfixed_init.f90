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
      ! [SS-GR-FINAL B7] DEFERRED: all symbols are config→globals copy targets.
      !   Retirement requires Phase C3 adapter rewrite (config_to_variables.f90 dual-writes).
      use variables, only: &
                            ! DEFERRED: idev/tsumea/tsumam — crop development config; Phase C3; tbase retired
                            idev, tsumea, tsumam,                         &
                            ! DEFERRED: kdif/kdir/gctb/swgc — radiation/crop factor config; Phase C3
                            gctb, swgc,                       &  ! kdif/kdir retired
                            ! DEFERRED: cftb/chtb/cfeictb/swcf/albedo/rsc/rsw — crop factor config; Phase C3
                            cftb, chtb, cfeictb, swcf, rsw,  &  ! albedo/rsc retired
                            ! DEFERRED: rdtb/rdctb — root config; Phase C3; swrd/swdmi2rd/swrdc/rdi/rri/rdc retired
                            rdtb, rdctb,                                &
                            ! DEFERRED: swoxygen/swWrtNonox/aeratecrit/max_resp_factor — O2 stress config; Phase C3
                            swoxygen, swWrtNonox, aeratecrit, max_resp_factor, &
                            ! DEFERRED: hlim1/hlim2u/hlim2l — drought stress limits; Phase C3
                            hlim1, hlim2u, hlim2l,                        &
                            ! DEFERRED: swdrought/hlim3h/hlim3l/hlim4/adcrh/adcrl — drought config; Phase C3
                            swdrought, hlim3h, hlim3l, hlim4, adcrh, adcrl, &
                            ! DEFERRED: swsalinity/saltmax/saltslope/salthead — salinity config; Phase C3
                            swsalinity, saltmax, saltslope, salthead,     &
                            ! DEFERRED: swcompensate/swstressor/alphacrit/dcritrtz — compensation config; Phase C3
                            swcompensate, swstressor, alphacrit, dcritrtz, &
                            ! DEFERRED: swinter/cofab — interception config; Phase C3
                            swinter,                                      &  ! cofab retired
                            ! DEFERRED: schedule/dvsend/swharv — harvest schedule config; Phase C3
                            schedule, dvsend, swharv,                     &
                            ! DEFERRED: cumdens — root density cumulative; Phase C3
                            cumdens
      use array_utils,  only: afgen
      use error_mod,    only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      type(cropfixed_config_t), intent(in)    :: cfg
      integer,                  intent(in)    :: icrop  ! reserved for swinco=3 inifil path (not yet ported)
      integer,                  intent(out)   :: lcc
      type(swap_state_t),       intent(inout) :: state  ! [SS-GR-ATM A5.1] runtime dual-write target

      integer      :: i
      real(real64) :: depth, rootdis(202), sum

      ! ---- Defense-in-depth: stub-error gates mirror cropfixed_config_validate.
      if (cfg%swdrought == 2 .or. cfg%swoxygen == 2 .or. cfg%swcf == 3 .or. &
          cfg%swharv == 1   .or. cfg%swcompensate /= 0 .or.                 &
          cfg%swinter == 2  .or. cfg%swinter == 3 .or.                      &
          cfg%swrd /= 1     .or. cfg%swsalinity /= 0 .or. cfg%schedule_switch == 1) then
         call fatalerr_collected('cropfixed_init_from_config', &
            'Unsupported runtime branch reached on the TOML path. The validator ' // &
            'should have caught this earlier.')
         return
      end if

      ! ---- Copy scalars / switches to globals ------------------------
      idev   = cfg%idev
      lcc    = cfg%lcc
      tsumea = cfg%tsumea
      tsumam = cfg%tsumam
      state%crop%common%tbase = cfg%tbase

      state%crop%kdif = cfg%kdif
      state%crop%kdir = cfg%kdir

      swgc = cfg%swgc
      swcf = cfg%swcf
      state%crop%swcf = swcf   ! [SS-GR-ATM A5.1] runtime dual-write
      state%crop%common%swrd     = cfg%swrd
      state%crop%common%swdmi2rd = cfg%swdmi2rd
      state%crop%common%swrdc    = cfg%swrdc
      state%crop%common%rdi = cfg%rdi
      state%crop%common%rri = cfg%rri
      state%crop%common%rdc = cfg%rdc

      swoxygen   = cfg%swoxygen
      swWrtNonox = cfg%swwrtnonox
      aeratecrit = cfg%aeratecrit
      max_resp_factor = cfg%max_resp_factor
      hlim1      = cfg%hlim1
      hlim2u     = cfg%hlim2u
      hlim2l     = cfg%hlim2l

      swdrought = cfg%swdrought
      hlim3h    = cfg%hlim3h
      hlim3l    = cfg%hlim3l
      hlim4     = cfg%hlim4
      adcrh     = cfg%adcrh
      adcrl     = cfg%adcrl

      swsalinity = cfg%swsalinity
      saltmax    = cfg%saltmax
      saltslope  = cfg%saltslope
      salthead   = cfg%salthead

      swcompensate = cfg%swcompensate
      swstressor   = cfg%swstressor
      alphacrit    = cfg%alphacrit
      dcritrtz     = cfg%dcritrtz

      swinter = cfg%swinter
      state%crop%cofab = cfg%cofab

      schedule = cfg%schedule_switch

      dvsend = cfg%dvsend
      swharv = cfg%swharv

      ! ---- Reflection coefficients / crop resistance (legacy:2150-2160)
      ! Two branches in legacy:
      !   swcf=1 (ETref crop factor): hardcoded ETref defaults
      !       (albedo=0.23, rsc=70, rsw=0)
      !   swcf=2 (crop height): user-authored values from .crp
      ! swcf=3 is stub-errored upstream so it never reaches here.
      if (cfg%swcf == 1) then
         state%crop%common%albedo = 0.23_real64
         state%crop%common%rsc    = 70.0_real64
         rsw    = 0.0_real64
      else if (cfg%swcf == 2) then
         state%crop%common%albedo = cfg%albedo
         state%crop%common%rsc    = cfg%rsc
         rsw    = cfg%rsw
      end if

      ! ---- Copy tables ------------------------------------------------
      ! gctb (size up to 2*magrs in legacy; we copy what was authored).
      if (allocated(cfg%gctb))  call copy_pair_table(cfg%gctb,  gctb)
      if (allocated(cfg%cftb))  then
        call copy_pair_table(cfg%cftb,  cftb)
        state%crop%fixed%cftb = cftb   ! [SS-GR-CROP A5.2]
      endif
      if (allocated(cfg%chtb))  then
        call copy_pair_table(cfg%chtb,  chtb)
        state%crop%fixed%chtb = chtb   ! [SS-GR-CROP A5.2]
      endif
      if (allocated(cfg%cfeictb)) then
        call copy_pair_table(cfg%cfeictb, cfeictb)
        state%crop%fixed%cfeictb = cfeictb   ! [SS-GR-CROP A5.2]
      endif
      if (allocated(cfg%rdtb))  call copy_pair_table(cfg%rdtb,  rdtb)

      ! rdctb is sized 22 in legacy; we always copy.
      if (allocated(cfg%rdctb)) then
         do i = 1, min(size(cfg%rdctb), size(rdctb))
            rdctb(i) = cfg%rdctb(i)
         end do
      end if

      ! ---- Tail-of-readcropfixed RUNTIME init: cumdens (legacy:2449-2478)
      if (cfg%swdrought == 1) then
         do i = 0, 100
            depth = 0.01_real64 * real(i, real64)
            rootdis(i*2 + 1) = depth
            rootdis(i*2 + 2) = afgen(rdctb, 22, depth)
         end do
         do i = 1, 202, 2
            cumdens(i) = rootdis(i)
         end do
         sum = 0.0_real64
         cumdens(2) = 0.0_real64
         do i = 4, 202, 2
            sum = sum + (rootdis(i-2) + rootdis(i)) * 0.5_real64 &
                      * (cumdens(i-1) - cumdens(i-3))
            cumdens(i) = sum
         end do
         do i = 2, 202, 2
            cumdens(i) = cumdens(i) / sum
         end do
         state%crop%common%cumdens = cumdens   ! [SS-GR-CROPRT A5]
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
