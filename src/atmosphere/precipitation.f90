!> Module for precipitation partitioning and snow/rain separation
!!
!! This module handles the partitioning of total precipitation into rain and snow
!! components based on temperature thresholds. It implements a simple temperature-index
!! approach with linear interpolation between transition temperatures.
!!
!! @author Original SWAP development team
!! @date Refactored February 2026
module precipitation_mod
   use swap_state_mod, only: swap_state_t
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private

   public :: PartitionPrecipitation
contains
   !> Partition precipitation into rain and snow based on temperature
   !!
   !! Separates total precipitation (grai) into rainfall and snowfall components
   !! using temperature-based thresholds. For daily meteorology, applies a
   !! temperature-index snow model. For detailed meteorology, aggregates sub-daily
   !! precipitation (snow calculations not supported).
   !!
   !! ## Temperature-Index Snow Model
   !!
   !! The partitioning uses two threshold temperatures:
   !! - \( T_{rain} \): Above this temperature, all precipitation is rain
   !! - \( T_{snow} \): Below this temperature, all precipitation is snow
   !!
   !! ### Partitioning Logic:
   !! - \( T_{av} > T_{rain} \): \( f_{snow} = 0 \) (all rain)
   !! - \( T_{av} < T_{snow} \): \( f_{snow} = 1 \) (all snow)
   !! - Between thresholds: \( f_{snow} = \frac{T_{rain} - T_{av}}{T_{rain} - T_{snow}} \)
   !!
   !! ## Outputs
   !! - `gsnow`: Snowfall depth (cm)
   !! - `snrai`: Rainfall on existing snowpack (cm)
   !! - `fprecnosnow`: Fraction of precipitation that reaches soil surface (-)
   !!
   !! @warning For detailed meteorology (swmetdetail=1), snow calculations are disabled
   !!
   !! @note
   !! **Original documentation:**
   !! Section 2: Rain and Snow partitioning
   !! Converts precipitation from mm to cm and determines rain/snow split
   !! @endnote
   subroutine PartitionPrecipitation(swmetdetail, swsnow, tav, TePrRain, TePrSnow, &
                                     ssnow, nmetdetail, arain, grai, gsnow, snrai, &
                                     fprecnosnow, restint, state)
      implicit none

      ! Arguments
      integer, intent(in)    :: swmetdetail    !! Meteo detail switch: 0=daily, 1=detailed
      integer, intent(in)    :: swsnow         !! Snow calculation switch: 0=off, 1=on
      real(8), intent(in)    :: tav            !! Average air temperature (°C)
      real(8), intent(in)    :: TePrRain       !! Temperature above which all precip is rain (°C)
      real(8), intent(in)    :: TePrSnow       !! Temperature below which all precip is snow (°C)
      real(8), intent(inout) :: ssnow          !! Current snow water equivalent (cm)
      integer, intent(in)    :: nmetdetail     !! Number of detailed meteo records per day
      real(8), intent(in)    :: arain(:)       !! Sub-daily precipitation array (cm)
      real(8), intent(inout) :: grai           !! Gross precipitation (mm for daily, cm for detailed)
      real(8), intent(out)   :: gsnow          !! Snowfall depth (cm)
      real(8), intent(out)   :: snrai          !! Rainfall on snowpack (cm)
      real(8), intent(out)   :: fprecnosnow    !! Fraction of precip reaching soil surface (-)
      real(8), intent(out)   :: restint        !! Remaining interception for detailed mode (cm)
      type(swap_state_t), intent(inout) :: state  !! [SS-ATM] atmosphere dual-write target

      ! Local variables
      integer :: i                             !! Loop counter
      real(8) :: snow_fraction                 !! Fraction of precipitation falling as snow (-)

      if (swmetdetail == 0) then
         ! === DAILY METEOROLOGY ===

         ! Convert precipitation from mm to cm
         grai = grai*0.1d0
         state%atmosphere%grai = grai

         if (swsnow == 1) then
            ! Temperature-based snow partitioning
            if (tav > TePrRain) then
               ! All precipitation as rain (warm conditions)
               gsnow = 0.0d0
               state%atmosphere%gsnow = 0.0_real64

            elseif (tav < TePrSnow) then
               ! All precipitation as snow (cold conditions)
               gsnow = grai
               state%atmosphere%gsnow = grai

            else
               ! Linear interpolation between transition temperatures
               snow_fraction = (TePrRain - tav)/(TePrRain - TePrSnow)
               gsnow = grai*snow_fraction
               state%atmosphere%gsnow = gsnow
            end if

            ! Calculate rainfall on existing snowpack (only counted if snowpack exists)
            if (ssnow > 1.0d-6) then
               snrai = grai - gsnow
               state%atmosphere%snrai = snrai
            else
               snrai = 0.0d0
               state%atmosphere%snrai = 0.0_real64
            end if

            ! Fraction of precipitation reaching soil surface
            ! (excludes snow accumulation and rain on snowpack)
            if (grai > 0.0d0) then
               fprecnosnow = 1.0d0 - (gsnow + snrai)/grai
               state%atmosphere%fprecnosnow = fprecnosnow
            else
               fprecnosnow = 0.0d0
               state%atmosphere%fprecnosnow = 0.0_real64
            end if

         else
            ! Snow calculations disabled
            gsnow = 0.0d0
            state%atmosphere%gsnow = 0.0_real64
            snrai = 0.0d0
            state%atmosphere%snrai = 0.0_real64
            fprecnosnow = 1.0d0
            state%atmosphere%fprecnosnow = 1.0_real64
         end if

      elseif (swmetdetail == 1) then
         ! === DETAILED METEOROLOGY ===

         ! Aggregate sub-daily precipitation to daily total
         grai = 0.0d0
         do i = 1, nmetdetail
            grai = grai + arain(i)
         end do
         state%atmosphere%grai = grai

         ! Initialize remaining interception storage
         restint = 0.0d0

         ! Snow calculations not supported for detailed meteorology
         gsnow = 0.0d0
         state%atmosphere%gsnow = 0.0_real64
         ssnow = 0.0d0  ! Note: This modifies a state variable - consider refactoring
         state%atmosphere%ssnow = 0.0_real64        ! (D11: verbatim mutation preserved)
         snrai = 0.0d0
         state%atmosphere%snrai = 0.0_real64
         fprecnosnow = 1.0d0
         state%atmosphere%fprecnosnow = 1.0_real64
      end if

   end subroutine PartitionPrecipitation

end module precipitation_mod
