!> [simulation] section config: dates, output timing switches.
module simulation_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   use validation_mod, only: check_int_range, check_real_range, &
                             check_int_enum, check_ordered_pair
   implicit none
   private

   public :: simulation_config_t
   public :: simulation_numerical_t

   !> Numerical solver controls (timestep + iteration limits).
   type :: simulation_numerical_t
      real(real64) :: dt        = 0.2_real64
      real(real64) :: dtmin     = 1.0e-6_real64
      real(real64) :: dtmax     = 0.2_real64
      integer      :: MaxIt     = 30
      integer      :: MaxBackTr = 3
      real(real64) :: taccur    = 1.0e-3_real64
   contains
      procedure :: validate => simulation_numerical_validate
   end type simulation_numerical_t

   type :: simulation_config_t
      real(real64) :: tstart = 0.0_real64   !! Days since 1900
      real(real64) :: tend   = 0.0_real64
      integer :: nprintday = 1
      integer :: swmonth   = 0
      integer :: period    = 1
      integer :: swres     = 0
      integer :: swodat    = 0
      integer :: swyrvar   = 0
      type(simulation_numerical_t) :: numerical
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
      call self%numerical%validate(errors)
   end subroutine simulation_config_validate

   subroutine simulation_numerical_validate(self, errors)
      class(simulation_numerical_t), intent(in)    :: self
      type(error_collection_t),      intent(inout) :: errors

      call check_real_range(self%dt,    1.0e-9_real64, 1.0_real64, &
                            "simulation.numerical.dt",    errors)
      call check_real_range(self%dtmin, 1.0e-12_real64, 1.0_real64, &
                            "simulation.numerical.dtmin", errors)
      call check_real_range(self%dtmax, 1.0e-9_real64, 1.0_real64, &
                            "simulation.numerical.dtmax", errors)
      call check_int_range(self%MaxIt, 1, 1000, &
                           "simulation.numerical.MaxIt", errors)
      call check_int_range(self%MaxBackTr, 1, 100, &
                           "simulation.numerical.MaxBackTr", errors)
      call check_real_range(self%taccur, 1.0e-9_real64, 1.0_real64, &
                            "simulation.numerical.taccur", errors)

      ! Cross-field invariants: dtmin <= dtmax, dtmin <= dt <= dtmax.
      if (self%dtmin > self%dtmax) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                            "dtmin > dtmax", "simulation.numerical")
      end if
      if (self%dt < self%dtmin) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                            "dt < dtmin", "simulation.numerical")
      end if
      if (self%dt > self%dtmax) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                            "dt > dtmax", "simulation.numerical")
      end if
   end subroutine simulation_numerical_validate

   subroutine simulation_config_finalize(self, errors)
      class(simulation_config_t), intent(inout) :: self
      type(error_collection_t),   intent(inout) :: errors
   end subroutine simulation_config_finalize

end module simulation_config_mod
