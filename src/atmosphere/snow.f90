!> Module for snow accumulation and melt processes
!!
!! This module handles the simulation of snow accumulation, sublimation,
!! and melt processes in the SWAP model.
!!
!! @note
!! Originally developed: December 2004
!! @endnote
module snow_mod
   use error_mod, only: fatalerr_collected
   use swap_state_mod, only: swap_state_t

   implicit none

   private
   public :: snow

contains
!> Simulate snow accumulation and melt processes
!!
!! This subroutine handles snow pack dynamics including:
!! - Snow accumulation from precipitation
!! - Sublimation losses
!! - Temperature-driven snowmelt
!! - Rain-on-snow melt enhancement
!! - Liquid water storage and drainage
!!
!! The energy balance approach uses air temperature as a proxy
!! for available melt energy.
!!
   subroutine snow(task, state)

      use Variables
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none

      ! Arguments
      integer, intent(in) :: task
      ! SS-ATM Phase 1 Task A-1.3: non-optional intent(inout)
      ! SS-ATM Phase 1 Task A-1.7: dual-writes installed
      type(swap_state_t), intent(inout) :: state
    !! Task selector: 1=initialization, 2=calculation

      ! Local variables
      real(8) :: smelt
    !! Snowmelt by temperature [cm swe]
      real(8) :: smeltr
    !! Snowmelt by rain [cm swe]
      real(8) :: SnDefit
    !! Snow deficit when pack becomes negative
      real(8) :: SnLoss
    !! Total snow loss (melt + sublimation)
      real(8) :: tsoil_surf
    !! Surface soil temperature [deg C], read from state%heat or global tsoil(1)

      ! Constants
      real(8), parameter :: cwat = 4180.0d0
    !! Specific heat of water [J/kg/K]
      real(8), parameter :: lm = 333580d0
    !! Latent heat of melting [J/kg]
      real(8), parameter :: ts = 0.0d0
    !! Snow temperature [°C]
      real(8) :: slw_max
    !! Maximum storage of liquid water in snow [cm/d]
      real(8) :: qlw
    !! Drainage flux from snow pack [cm/d]

      ! ----------------------------------------------------------------------

      select case (task)
      case (1)

         ! === initialization ===================================================

         if (swinco .eq. 3) then
            snowinco = ssnow
            state%atmosphere%snowinco = snowinco  ! [SS-ATM] dual-write
         else
            ssnow = snowinco
            state%atmosphere%ssnow = ssnow  ! [SS-ATM] dual-write
         end if

         return

      case (2)

         ! SS-ATM Phase 1 Task A-1.7: ASSOCIATE aliases for dense dual-write body
         associate( &
            at_ssnow   => state%atmosphere%ssnow,       &
            at_snowinco => state%atmosphere%snowinco,    &
            at_melt    => state%atmosphere%melt,         &
            at_subl    => state%atmosphere%subl,         &
            at_slw     => state%atmosphere%slw,          &
            at_peva    => state%atmosphere%peva,         &
            at_empreva => state%atmosphere%empreva,      &
            at_igsnow  => state%atmosphere%intr%igsnow,  &
            at_isubl   => state%atmosphere%intr%isubl,   &
            at_isnrai  => state%atmosphere%intr%isnrai,  &
            at_cgsnow  => state%atmosphere%cumu%cgsnow,  &
            at_csubl   => state%atmosphere%cumu%csubl,   &
            at_csnrai  => state%atmosphere%cumu%csnrai,  &
            at_cmelt   => state%atmosphere%cumu%cmelt    &
         )

         ! === snow pack rate and state variables ===============================

         ! --- reset intermediate snow states
         if (flzerointr) then
            igsnow = 0.0d0
            at_igsnow = 0.0_real64  ! [SS-ATM] dual-write intr reset
            isubl = 0.0d0
            at_isubl = 0.0_real64   ! [SS-ATM] dual-write intr reset
            isnrai = 0.0d0
            at_isnrai = 0.0_real64  ! [SS-ATM] dual-write intr reset
            ISsnowBeg = Ssnow
         end if

         ! --- reset cumulative snow states
         if (flzerocumu) then
            cgsnow = 0.0d0
            at_cgsnow = 0.0_real64  ! [SS-ATM] dual-write cumu reset
            csubl = 0.0d0
            at_csubl = 0.0_real64   ! [SS-ATM] dual-write cumu reset
            csnrai = 0.0d0
            at_csnrai = 0.0_real64  ! [SS-ATM] dual-write cumu reset
            cmelt = 0.0d0
            at_cmelt = 0.0_real64   ! [SS-ATM] dual-write cumu reset
            snowinco = ssnow
            at_snowinco = snowinco  ! [SS-ATM] dual-write
         end if

         ! --- when there is snowpack calculate the amount of sublimation
         subl = 0.0d0
         at_subl = 0.0_real64  ! [SS-ATM] dual-write
         if (swsublim .eq. 0) then
            if (ssnow .gt. 0.0d0) then
               subl = peva
               at_subl = subl       ! [SS-ATM] dual-write
               if (swetsine .eq. 1) then
                  subl = pevaday
                  at_subl = subl    ! [SS-ATM] dual-write
               end if
               empreva = 0.0d0
               at_empreva = 0.0_real64  ! [SS-ATM] dual-write
               peva = 0.0d0
               at_peva = 0.0_real64     ! [SS-ATM] dual-write
            end if
         end if

         ! --- when the soil surface is above the freezing point there will be
         ! --- no accumulation of fresh snow.
         ! SS-ATM Phase 1 Task A-1.3: state is now mandatory — read tsoil(1) directly
         tsoil_surf = state%heat%tsoil(1)
         if (tsoil_surf .gt. 0.5d0 .and. ssnow .lt. 1.0d-6 .and. gsnow .gt. 0.0d0) then
            ssnow = 0.0d0
            at_ssnow = 0.0_real64  ! [SS-ATM] dual-write
            melt = gsnow
            at_melt = melt         ! [SS-ATM] dual-write
            subl = 0.d0
            at_subl = 0.0_real64   ! [SS-ATM] dual-write
         else

            ! --- amount of snowmelt [cm swe] negative values of smelt: see 'melt = '
            smelt = snowcoef*(tav - ts)

            ! --- extra snowmelt when there falls rain on the snowpack [cm swe]
            if (snrai .gt. 0.0d0) then
               smeltr = snrai*cwat*(tav - ts)/lm
            else
               smeltr = 0.0d0
            end if

            ! --- total snowmelt [cm swe]; negative values of smelt can partly compensate smeltr
            melt = max(0.0d0, (smelt + smeltr))
            at_melt = melt  ! [SS-ATM] dual-write

            ! --- amount of snow left [cm swe] without storage of liquid water slw
            ssnow = ssnow + gsnow - subl - melt - slw
            at_ssnow = ssnow  ! [SS-ATM] dual-write

            ! --- potential amount of liquid water storage
            slw = slw + snrai
            at_slw = slw  ! [SS-ATM] dual-write

            ! --- maximum retention of liquid water in snow is fraction 0.07 of total water storage
            slw_max = 0.07*(slw + ssnow)

            ! --- drainage of liquid water from snow
            qlw = max(0.0d0, slw - slw_max)

            ! --- remaining storage of liquid water in snow
            slw = slw - qlw
            at_slw = slw  ! [SS-ATM] dual-write

            ! --- reset total snow storage and total melt
            ssnow = ssnow + slw
            at_ssnow = ssnow  ! [SS-ATM] dual-write
            melt = melt + qlw
            at_melt = melt    ! [SS-ATM] dual-write

            ! --- in case of snow deficit: adapt snow loss terms melt and sublimation
            if (ssnow .lt. 0.0d0) then
               SnDefit = -Ssnow
               SnLoss = melt + subl
               melt = (1.d0 - SnDefit/SnLoss)*melt
               at_melt = melt  ! [SS-ATM] dual-write
               subl = (1.d0 - SnDefit/SnLoss)*subl
               at_subl = subl  ! [SS-ATM] dual-write
               Ssnow = 0.d0
               at_ssnow = 0.0_real64  ! [SS-ATM] dual-write
               slw = 0.d0
               at_slw = 0.0_real64    ! [SS-ATM] dual-write
            end if
         end if

         ! --- set cumulative amounts
         igsnow = igsnow + gsnow
         at_igsnow = igsnow  ! [SS-ATM] dual-write intr
         isubl = isubl + subl
         at_isubl = isubl    ! [SS-ATM] dual-write intr
         isnrai = isnrai + snrai
         at_isnrai = isnrai  ! [SS-ATM] dual-write intr
         cgsnow = cgsnow + gsnow
         at_cgsnow = cgsnow  ! [SS-ATM] dual-write cumu
         csubl = csubl + subl
         at_csubl = csubl    ! [SS-ATM] dual-write cumu
         cmelt = cmelt + melt
         at_cmelt = cmelt    ! [SS-ATM] dual-write cumu
         csnrai = csnrai + snrai
         at_csnrai = csnrai  ! [SS-ATM] dual-write cumu

         end associate

      case default
         call fatalerr_collected('Snow', 'Illegal value for TASK')
      end select

      return
   end subroutine snow

end module snow_mod
