!> Module for precipitation partitioning and snow/rain separation
!!
!! This module handles the partitioning of total precipitation into rain and snow
!! components based on temperature thresholds. It implements a simple temperature-index
!! approach with linear interpolation between transition temperatures.
!!
!! @author Original SWAP development team
!! @date Refactored February 2026
module precipitation_mod
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private

   public :: PartitionPrecipitation
contains
   !> Partition daily precipitation into rain and snow based on temperature.
   !!
   !! Daily mode (swmetdetail=0): apply temperature-index snow model to today's
   !! gross precipitation, partition into snowfall, rain-on-snowpack and
   !! soil-reaching fractions.
   !!
   !! Detailed mode (swmetdetail=1): aggregate sub-daily precipitation; snow
   !! calculations not supported (ssnow forced to 0).
   !!
   !! All inputs are pulled from state%atmosphere/state%timecontrol and
   !! config%meteo; all outputs are written to state%atmosphere.
   pure subroutine PartitionPrecipitation(state, config)
      implicit none

      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config

      integer :: i
      real(8) :: snow_fraction

      associate (atmo => state%atmosphere, time => state%timecontrol, meteo => config%meteo)

      if (meteo%swmetdetail == 0) then
         ! === DAILY METEOROLOGY ===
         ! Pull today's gross precipitation from the per-day array; convert mm → cm.
         atmo%grai = atmo%arai(time%daymeteo + 1 - atmo%daynrfirst) * 0.1_real64

         if (meteo%snow%swsnow == 1) then
            ! Temperature-based snow partitioning
            if (atmo%Tav > atmo%teprrain) then
               atmo%gsnow = 0.0_real64                             ! all rain (warm)
            else if (atmo%Tav < atmo%teprsnow) then
               atmo%gsnow = atmo%grai                              ! all snow (cold)
            else
               snow_fraction = (atmo%teprrain - atmo%Tav) / (atmo%teprrain - atmo%teprsnow)
               atmo%gsnow = atmo%grai * snow_fraction              ! linear blend
            end if

            ! Rainfall on existing snowpack (only when snowpack exists)
            if (atmo%ssnow > 1.0e-6_real64) then
               atmo%snrai = atmo%grai - atmo%gsnow
            else
               atmo%snrai = 0.0_real64
            end if

            ! Fraction of precipitation reaching soil surface
            if (atmo%grai > 0.0_real64) then
               atmo%fprecnosnow = 1.0_real64 - (atmo%gsnow + atmo%snrai) / atmo%grai
            else
               atmo%fprecnosnow = 0.0_real64
            end if
         else
            atmo%gsnow       = 0.0_real64
            atmo%snrai       = 0.0_real64
            atmo%fprecnosnow = 1.0_real64
         end if

      else if (meteo%swmetdetail == 1) then
         ! === DETAILED METEOROLOGY ===
         ! Aggregate sub-daily precipitation to daily total.
         atmo%grai = 0.0_real64
         do i = 1, meteo%nmetdetail
            atmo%grai = atmo%grai + atmo%arain_subdaily(i)
         end do

         ! Snow calculations not supported in detailed mode.
         atmo%restint     = 0.0_real64
         atmo%gsnow       = 0.0_real64
         atmo%ssnow       = 0.0_real64
         atmo%snrai       = 0.0_real64
         atmo%fprecnosnow = 1.0_real64
      end if

      end associate
   end subroutine PartitionPrecipitation

end module precipitation_mod
