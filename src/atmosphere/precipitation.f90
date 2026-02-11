!> Module for precipitation partitioning and snow/rain separation
!!
!! This module handles the partitioning of total precipitation into rain and snow
!! components based on temperature thresholds. It implements a simple temperature-index
!! approach with linear interpolation between transition temperatures.
!!
!! @author Original SWAP development team
!! @date Refactored February 2026
module precipitation_mod
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
                                     fprecnosnow, restint)
      implicit none

      ! Arguments
      integer, intent(in)    :: swmetdetail  !! Meteo detail switch: 0=daily, 1=detailed
      integer, intent(in)    :: swsnow       !! Snow calculation switch: 0=off, 1=on
      real(8), intent(in)    :: tav          !! Average air temperature (°C)
      real(8), intent(in)    :: TePrRain     !! Temperature above which all precip is rain (°C)
      real(8), intent(in)    :: TePrSnow     !! Temperature below which all precip is snow (°C)
      real(8), intent(inout)    :: ssnow     !! Current snow water equivalent (cm)
      integer, intent(in)    :: nmetdetail   !! Number of detailed meteo records per day
      real(8), intent(in)    :: arain(:)     !! Sub-daily precipitation array (cm)
      real(8), intent(inout) :: grai         !! Gross precipitation (mm for daily, cm for detailed)
      real(8), intent(out)   :: gsnow        !! Snowfall depth (cm)
      real(8), intent(out)   :: snrai        !! Rainfall on snowpack (cm)
      real(8), intent(out)   :: fprecnosnow  !! Fraction of precip reaching soil surface (-)
      real(8), intent(out)   :: restint      !! Remaining interception for detailed mode (cm)

      ! Local variables
      integer :: i              !! Loop counter
      real(8) :: snow_fraction  !! Fraction of precipitation falling as snow (-)

      if (swmetdetail == 0) then
         ! === DAILY METEOROLOGY ===

         ! Convert precipitation from mm to cm
         grai = grai*0.1d0

         if (swsnow == 1) then
            ! Temperature-based snow partitioning
            if (tav > TePrRain) then
               ! All precipitation as rain (warm conditions)
               gsnow = 0.0d0

            elseif (tav < TePrSnow) then
               ! All precipitation as snow (cold conditions)
               gsnow = grai

            else
               ! Linear interpolation between transition temperatures
               snow_fraction = (TePrRain - tav)/(TePrRain - TePrSnow)
               gsnow = grai*snow_fraction
            end if

            ! Calculate rainfall on existing snowpack
            ! (only counted if snowpack exists)
            if (ssnow > 1.0d-6) then
               snrai = grai - gsnow
            else
               snrai = 0.0d0
            end if

            ! Fraction of precipitation reaching soil surface
            ! (excludes snow accumulation and rain on snowpack)
            if (grai > 0.0d0) then
               fprecnosnow = 1.0d0 - (gsnow + snrai)/grai
            else
               fprecnosnow = 0.0d0
            end if

         else
            ! Snow calculations disabled
            gsnow = 0.0d0
            snrai = 0.0d0
            fprecnosnow = 1.0d0
         end if

      elseif (swmetdetail == 1) then
         ! === DETAILED METEOROLOGY ===

         ! Aggregate sub-daily precipitation to daily total
         grai = 0.0d0
         do i = 1, nmetdetail
            grai = grai + arain(i)
         end do

         ! Initialize remaining interception storage
         restint = 0.0d0

         ! Snow calculations not supported for detailed meteorology
         gsnow = 0.0d0
         ssnow = 0.0d0  ! Note: This modifies a state variable - consider refactoring
         snrai = 0.0d0
         fprecnosnow = 1.0d0
      end if

   end subroutine PartitionPrecipitation
end module precipitation_mod
