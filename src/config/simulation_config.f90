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

   !> Numerical solver controls (timestep + iteration limits + Richards
   !! convergence criteria + hydraulic-conductivity averaging).
   !! Defaults mirror the .swp template values used by the regression
   !! cases — so a TOML that omits this block runs identically to a
   !! .swp that authors all of these explicitly.
   type :: simulation_numerical_t
      real(real64) :: dt            = 0.2_real64
      real(real64) :: dtmin         = 1.0e-6_real64
      real(real64) :: dtmax         = 0.2_real64
      integer      :: MaxIt         = 30
      integer      :: MaxBackTr     = 3
      real(real64) :: taccur        = 1.0e-3_real64
      ! Richards-solver convergence criteria (legacy .swp Part 13).
      real(real64) :: gwlconv       = 100.0_real64    !! [1e-5..1000 cm]
      real(real64) :: critdevh1cp   = 0.01_real64     !! [1e-10..1e3 -]
      real(real64) :: critdevh2cp   = 0.1_real64      !! [1e-10..1e3 cm]
      real(real64) :: critdevponddt = 1.0e-4_real64   !! [1e-6..0.1 cm]
      ! Hydraulic conductivity mean type + implicitness flag.
      integer      :: swkmean       = 1               !! 1..6 (see readswap.f90:983-985)
      integer      :: swkimpl       = 0               !! 0=explicit, 1=implicit
      ! Hard-cap on iterations per day. The legacy default is 1e8 (effectively
      ! unlimited); we mirror that. The runtime aborts if a single day's
      ! iteration count exceeds this — the cap exists only as a circuit-breaker.
      integer      :: msteps        = 100000000
      ! Expert-only Richards switches (default off):
      logical      :: swcaprise           = .false.   !! Pin K at the deepest root node and below to 1e-10 to suppress capillary rise into the root zone (for experts).
      logical      :: dump_convergence_diagnostics = .false.   !! Emit per-step log_debug entries with the Richards convergence diagnostics (debug feature).
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
      integer :: swcrp     = 0     !! [IO-OUT/D] 0=no .crp output, 1=write legacy .crp crop output file
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
      call check_int_enum(self%swcrp,   [0, 1], "simulation.output.swcrp", errors)
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
      call check_real_range(self%gwlconv,       1.0e-5_real64,  1000.0_real64, &
                            "simulation.numerical.gwlconv", errors)
      call check_real_range(self%critdevh1cp,   1.0e-10_real64, 1.0e3_real64, &
                            "simulation.numerical.critdevh1cp", errors)
      call check_real_range(self%critdevh2cp,   1.0e-10_real64, 1.0e3_real64, &
                            "simulation.numerical.critdevh2cp", errors)
      call check_real_range(self%critdevponddt, 1.0e-6_real64,  0.1_real64, &
                            "simulation.numerical.critdevponddt", errors)
      call check_int_range(self%swkmean, 1, 6, &
                           "simulation.numerical.swkmean", errors)
      call check_int_enum(self%swkimpl, [0, 1], &
                          "simulation.numerical.swkimpl", errors)
      call check_int_range(self%msteps, 1, 1000000000, &
                           "simulation.numerical.msteps", errors)

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
