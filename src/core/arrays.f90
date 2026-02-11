!> Module containing global array dimension parameters for the SWAP model
!!
!! This module defines maximum array sizes used throughout the SWAP
!! (Soil-Water-Atmosphere-Plant) model. These parameters control memory
!! allocation and set upper limits for various model components.
!!
!! @note These parameters should be adjusted based on the specific
!!       simulation requirements and available system memory.
!!
!! @author Original SWAP development team
!! @date 2018-01-08
!! @version $Id: arrays.fi 362 2018-01-08 13:08:33Z kroes006 $
module swap_array_dimensions
  implicit none
  public

  !> Maximum number of years in the simulation period
  integer, parameter :: MAYRS = 200

  !> Maximum number of crops that can be simulated
  integer, parameter :: MACROP = 200

  !> Maximum number of soil compartments (discretization layers)
  integer, parameter :: MACP = 5000

  !> Maximum number of days in the simulation period (MAYRS × 366)
  integer, parameter :: MADAY = MAYRS*366

  !> Maximum number of drainage systems
  integer, parameter :: MADR = 5

  !> Maximum number of crop growth stages
  integer, parameter :: MAGRS = 366

  !> Maximum number of soil horizons
  integer, parameter :: MAHO = 1000

  !> Maximum number of time-dependent values for bottom boundary conditions
  integer, parameter :: MABBC = MADAY

  !> Maximum number of scaling factors for soil hydraulic properties
  integer, parameter :: MASCALE = 100

  !> Maximum number of rainfall records for detailed rainfall input
  integer, parameter :: MRAIN = 40000

  !> Maximum number of applied irrigation events
  integer, parameter :: MAIRG = 10000

  !> Maximum number of specified output dates
  integer, parameter :: MAOUT = 3000

  !> Maximum number of open water levels for basic drainage routine
  integer, parameter :: MAOWL = 10*366

  !> Maximum number of soil management events
  integer, parameter :: MASME = 1000

  !> Maximum number of data pairs in input tables for soil hydraulic relations
  integer, parameter :: MATAB = 1000

  !> Maximum number of total entries in soil hydraulic relation tables
  integer, parameter :: MATABENTRIES = 50005

  !> Maximum number of water levels in primary drainage system (extended drainage)
  integer, parameter :: MAWLP = 10*366

  !> Maximum number of water levels in secondary drainage system (extended drainage)
  integer, parameter :: MAWLS = 10*366

  !> Maximum number of surface water management periods (extended drainage)
  integer, parameter :: MAMP = 10*366

  !> Maximum number of surface water management table entries (extended drainage)
  integer, parameter :: MAMTE = 25

  !> Maximum number of macropore domains
  integer, parameter :: MADM = 20

  !> Maximum number of static equilibrium relations
  integer, parameter :: MASTEQ = 10000

  !> Maximum number of weather records in one year (48 records per day)
  integer, parameter :: NMETFILE = 17568

end module swap_array_dimensions
