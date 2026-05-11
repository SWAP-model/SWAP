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
      implicit none

      ! Arguments
      integer, intent(in) :: task
      ! SS-ATM Phase 1 Task A-1.3: non-optional intent(inout) — prepares for A-1.7 dual-writes
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
         else
            ssnow = snowinco
         end if

         return

      case (2)

         ! === snow pack rate and state variables ===============================

         ! --- reset intermediate snow states
         if (flzerointr) then
            igsnow = 0.0d0
            isubl = 0.0d0
            isnrai = 0.0d0
            ISsnowBeg = Ssnow
         end if

         ! --- reset cumulative snow states
         if (flzerocumu) then
            cgsnow = 0.0d0
            csubl = 0.0d0
            csnrai = 0.0d0
            cmelt = 0.0d0
            snowinco = ssnow
         end if

         ! --- when there is snowpack calculate the amount of sublimation
         subl = 0.0d0
         if (swsublim .eq. 0) then
            if (ssnow .gt. 0.0d0) then
               subl = peva
               if (swetsine .eq. 1) subl = pevaday
               empreva = 0.0d0
               peva = 0.0d0
            end if
         end if

         ! --- when the soil surface is above the freezing point there will be
         ! --- no accumulation of fresh snow.
         ! SS-ATM Phase 1 Task A-1.3: state is now mandatory — read tsoil(1) directly
         tsoil_surf = state%heat%tsoil(1)
         if (tsoil_surf .gt. 0.5d0 .and. ssnow .lt. 1.0d-6 .and. gsnow .gt. 0.0d0) then
            ssnow = 0.0d0
            melt = gsnow
            subl = 0.d0
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

            ! --- amount of snow left [cm swe] without storage of liquid water slw
            ssnow = ssnow + gsnow - subl - melt - slw

            ! --- potential amount of liquid water storage
            slw = slw + snrai

            ! --- maximum retention of liquid water in snow is fraction 0.07 of total water storage
            slw_max = 0.07*(slw + ssnow)

            ! --- drainage of liquid water from snow
            qlw = max(0.0d0, slw - slw_max)

            ! --- remaining storage of liquid water in snow
            slw = slw - qlw

            ! --- reset total snow storage and total melt
            ssnow = ssnow + slw
            melt = melt + qlw

            ! --- in case of snow deficit: adapt snow loss terms melt and sublimation
            if (ssnow .lt. 0.0d0) then
               SnDefit = -Ssnow
               SnLoss = melt + subl
               melt = (1.d0 - SnDefit/SnLoss)*melt
               subl = (1.d0 - SnDefit/SnLoss)*subl
               Ssnow = 0.d0
               slw = 0.d0
            end if
         end if

         ! --- set cumulative amounts
         igsnow = igsnow + gsnow
         isubl = isubl + subl
         isnrai = isnrai + snrai
         cgsnow = cgsnow + gsnow
         csubl = csubl + subl
         cmelt = cmelt + melt
         csnrai = csnrai + snrai

      case default
         call fatalerr_collected('Snow', 'Illegal value for TASK')
      end select

      return
   end subroutine snow

end module snow_mod
