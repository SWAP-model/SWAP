!> [simulation] section config: dates, output timing switches.
module simulation_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   use validation_mod, only: check_int_range, check_int_enum, check_ordered_pair
   implicit none
   private

   public :: simulation_config_t

   type :: simulation_config_t
      real(real64) :: tstart = 0.0_real64   !! Days since 1900
      real(real64) :: tend   = 0.0_real64
      integer :: nprintday = 1
      integer :: swmonth   = 0
      integer :: period    = 1
      integer :: swres     = 0
      integer :: swodat    = 0
      integer :: swyrvar   = 0
   contains
      procedure :: validate => simulation_config_validate
      procedure :: finalize => simulation_config_finalize
   end type simulation_config_t

contains

   subroutine simulation_config_validate(self, errors)
      class(simulation_config_t), intent(in)    :: self
      type(error_collection_t),   intent(inout) :: errors

      call check_int_range(self%nprintday, 1, 1440, "simulation.nprintday", errors)
      call check_int_enum(self%swmonth, [0, 1], "simulation.swmonth", errors)
      call check_int_range(self%period,  0, 366, "simulation.period",  errors)
      call check_int_enum(self%swres,   [0, 1], "simulation.swres",    errors)
      call check_int_enum(self%swodat,  [0, 1], "simulation.swodat",   errors)
      call check_int_enum(self%swyrvar, [0, 1], "simulation.swyrvar",  errors)
      call check_ordered_pair(self%tstart, self%tend, &
                              "tstart", "tend", "simulation", errors)
   end subroutine simulation_config_validate

   subroutine simulation_config_finalize(self, errors)
      class(simulation_config_t), intent(inout) :: self
      type(error_collection_t),   intent(inout) :: errors
   end subroutine simulation_config_finalize

end module simulation_config_mod
