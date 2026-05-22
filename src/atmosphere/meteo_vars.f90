!> Module containing shared meteorological variables
!!
!! This module provides shared variables for communication between meteorological
!! processing routines (ReadMeteoDay and ProcessMeteoDay). Variables are used to
!! store temporary meteorological data, intermediate calculations, and results
!! during the processing of daily or sub-daily meteorological inputs.
!!
!! ## Variable Categories
!! - Loop control: count, first, i, irecord, last, ndayparts
!! - Meteorological arrays: arain(96), awind(96)
!! - Interception: restint, interc, aintc, eintc
!! - Flux components: netrainflux, rainflux, Edirect, Tdirect, Tdirectwet, Edirectpond
!! - Meteorological scalars: etr, gctp, hum, svp, wfrac, win, dttp, sumtav
!!
!! @author Original SWAP development team
!! @date Refactored February 2026
module MeteoVars
   implicit none
   
   integer   count           !! Loop counter for meteorological records
   integer   first           !! First record index in processing range
   integer   i               !! General loop counter
   integer   irecord         !! Current record index
   integer   last            !! Last record index in processing range
   integer   ndayparts       !! Number of timesteps per day (1 for daily, >1 for sub-daily)
   
   real(8)   arain(96)       !! Array of precipitation values (cm)
   real(8)   awind(96)       !! Array of wind speed values (m/s)
   real(8)   restint         !! Remaining interception from previous timestep (cm)
   real(8)   interc          !! Current interception amount (cm)
   real(8)   Edirectpond     !! Direct evaporation from ponded water (mm/d)
   real(8)   aintc           !! Daily interception (cm)
   real(8)   dttp            !! Timestep duration for interception calculations (d)
   real(8)   eintc           !! Evaporated interception (cm)
   real(8)   etr             !! Reference evapotranspiration (mm/d)
   real(8)   gctp            !! Ground cover at current timestep (fraction)
   real(8)   hum             !! Humidity or vapor pressure (kPa)
   real(8)   netrainflux     !! Net rainfall flux after interception (cm)
   real(8)   rainflux        !! Gross rainfall flux (cm)
   real(8)   sumtav          !! Sum of temperatures for averaging
   real(8)   svp             !! Saturated vapor pressure (kPa)
   real(8)   wfrac           !! Wet fraction of crop canopy (fraction)
   real(8)   win             !! Wind speed (m/s)
   real(8)   Edirect         !! Direct soil evaporation (mm/d)
   real(8)   Tdirect         !! Direct transpiration (mm/d)
   real(8)   Tdirectwet      !! Direct transpiration from wet canopy (mm/d)
end module MeteoVars
