!> Module for snow accumulation and melt processes
!!
!! This module handles the simulation of snow accumulation, sublimation,
!! and melt processes in the SWAP model.
!!
!! @note
!! Originally developed: December 2004
!! @endnote
module snow_mod
   use swap_state_mod, only: swap_state_t
   use atmosphere_constants_mod, only: SPECIFIC_HEAT_WATER, LATENT_HEAT_MELTING, &
                                       SNOW_TEMPERATURE_C, SNOW_LIQUID_WATER_FRACTION, &
                                       SOIL_SURFACE_FREEZE_THRESHOLD_C

   implicit none

   private
   public :: snow_init, snow_step

contains
!> Initialize snow pack state (formerly snow(task=1, ...))
!!
!! Seeds state%atmosphere%snowinco or state%atmosphere%ssnow based on
!! the initial-conditions switch swinco.
!!
   subroutine snow_init(state)

      implicit none

      type(swap_state_t), intent(inout) :: state

      associate (atmo => state%atmosphere)

         if (state%cfg%soil%swinco .eq. 3) then
            atmo%snowinco = atmo%ssnow
         else
            atmo%ssnow = atmo%snowinco
         end if

      end associate

   end subroutine snow_init

!> Simulate snow accumulation and melt processes (formerly snow(task=2, ...))
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
   subroutine snow_step(state)

      use Variables, only: ISsnowBeg   ! still a bare global; state-home migration pending
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none

      type(swap_state_t), intent(inout) :: state

      real(8) :: smelt       !! Snowmelt by temperature [cm swe]
      real(8) :: smeltr      !! Snowmelt by rain [cm swe]
      real(8) :: SnDefit     !! Snow deficit when pack becomes negative
      real(8) :: SnLoss      !! Total snow loss (melt + sublimation)
      real(8) :: tsoil_surf  !! Surface soil temperature [deg C]
      real(8) :: slw_max     !! Maximum storage of liquid water in snow [cm/d]
      real(8) :: qlw         !! Drainage flux from snow pack [cm/d]

      associate (atmo => state%atmosphere,      &
                 intr => state%atmosphere%intr, &
                 cumu => state%atmosphere%cumu, &
                 heat => state%heat,            &
                 time => state%timecontrol)

      ! --- reset intermediate snow states
      if (time%flZeroIntr) then
         intr%igsnow = 0.0_real64
         intr%isubl  = 0.0_real64
         intr%isnrai = 0.0_real64
         ISsnowBeg   = atmo%ssnow
      end if

      ! --- reset cumulative snow states
      if (time%flZeroCumu) then
         cumu%cgsnow   = 0.0_real64
         cumu%csubl    = 0.0_real64
         cumu%csnrai   = 0.0_real64
         cumu%cmelt    = 0.0_real64
         atmo%snowinco = atmo%ssnow
      end if

      ! --- when there is snowpack calculate the amount of sublimation
      atmo%subl = 0.0_real64
      if (state%cfg%soil%frost%swsublim .eq. 0) then
         if (atmo%ssnow .gt. 0.0d0) then
            atmo%subl = atmo%peva
            if (state%cfg%meteo%swetsine .eq. 1) then
               atmo%subl = atmo%pevaday
            end if
            atmo%empreva = 0.0_real64
            atmo%peva    = 0.0_real64
         end if
      end if

      ! --- when the soil surface is above the freezing point there will be
      ! --- no accumulation of fresh snow.
      tsoil_surf = heat%tsoil(1)
      if (tsoil_surf .gt. SOIL_SURFACE_FREEZE_THRESHOLD_C .and. &
          atmo%ssnow .lt. 1.0d-6 .and. atmo%gsnow .gt. 0.0d0) then
         atmo%ssnow = 0.0_real64
         atmo%melt  = atmo%gsnow
         atmo%subl  = 0.0_real64
      else

         ! --- amount of snowmelt [cm swe]; negative values can partly compensate smeltr
         smelt = state%cfg%meteo%snow%snowcoef*(atmo%Tav - SNOW_TEMPERATURE_C)

         ! --- extra snowmelt when rain falls on the snowpack [cm swe]
         if (atmo%snrai .gt. 0.0d0) then
            smeltr = atmo%snrai*SPECIFIC_HEAT_WATER*(atmo%Tav - SNOW_TEMPERATURE_C)/LATENT_HEAT_MELTING
         else
            smeltr = 0.0d0
         end if

         ! --- total snowmelt [cm swe]
         atmo%melt = max(0.0d0, (smelt + smeltr))

         ! --- amount of snow left [cm swe] without storage of liquid water slw
         atmo%ssnow = atmo%ssnow + atmo%gsnow - atmo%subl - atmo%melt - atmo%slw

         ! --- potential amount of liquid water storage
         atmo%slw = atmo%slw + atmo%snrai

         ! --- maximum retention of liquid water in snow
         slw_max = SNOW_LIQUID_WATER_FRACTION*(atmo%slw + atmo%ssnow)

         ! --- drainage of liquid water from snow
         qlw = max(0.0d0, atmo%slw - slw_max)

         ! --- remaining storage of liquid water in snow
         atmo%slw = atmo%slw - qlw

         ! --- reset total snow storage and total melt
         atmo%ssnow = atmo%ssnow + atmo%slw
         atmo%melt  = atmo%melt + qlw

         ! --- in case of snow deficit: adapt snow loss terms melt and sublimation
         if (atmo%ssnow .lt. 0.0d0) then
            SnDefit = -atmo%ssnow
            SnLoss  = atmo%melt + atmo%subl
            atmo%melt  = (1.d0 - SnDefit/SnLoss)*atmo%melt
            atmo%subl  = (1.d0 - SnDefit/SnLoss)*atmo%subl
            atmo%ssnow = 0.0_real64
            atmo%slw   = 0.0_real64
         end if
      end if

      ! --- set cumulative amounts
      intr%igsnow = intr%igsnow + atmo%gsnow
      intr%isubl  = intr%isubl  + atmo%subl
      intr%isnrai = intr%isnrai + atmo%snrai
      cumu%cgsnow = cumu%cgsnow + atmo%gsnow
      cumu%csubl  = cumu%csubl  + atmo%subl
      cumu%cmelt  = cumu%cmelt  + atmo%melt
      cumu%csnrai = cumu%csnrai + atmo%snrai

      end associate

   end subroutine snow_step

end module snow_mod
