!> @file agetracer.f90
!! AgeTracer feature — extracted from solute.f90 during ADR 0032
!! to declutter the solute compute module.
!!
!! Status: DEAD. flAgeTracer is never set to .true. anywhere in the
!! codebase. The runtime body was removed during the solute Pattern 1
!! refactor (it was creating a perpetual dead-code dependency on
!! retiring bare globals). The Goode (1996) age-tracer implementation
!! is preserved in git history: see solute.f90:304-547 in commits
!! before ADR 0032, or agetracer.f90 before the 2026-05-24 deletion.
!!
!! Reactivation checklist (DO NOT skip — these decisions defer until
!! AgeTracer is genuinely needed):
!!
!!   1. Define agetracer_state_t. The 12 AgeTracer-specific globals
!!      formerly in variables.f90 (Ageirr, Agedrain, Agepre, Agepond,
!!      Agepondm1, icAgetopupw, icAgetopdwn, icAgeBot, icAgeDra,
!!      icAgeRot, icAgeSur, AgeGwl1m) are the target field set.
!!
!!   2. Resolve the cml dual-use (discovery hazard #1). AgeTracer
!!      formerly overwrote state%solute%cml(i) with age storage at
!!      the end of each task=2 call. Either:
!!      (a) Add a separate state%agetracer%cml_age(:) field and
!!          stop overwriting state%solute%cml; OR
!!      (b) Document the dual-use sequencing rules formally and
!!          retain the overwrite (faster but more fragile).
!!
!!   3. Wire flAgeTracer from typed config or external trigger
!!      (currently never assigned anywhere).
!!
!!   4. Update outage's flAgeTracer guard (in swapoutput.f90) to
!!      ungate output once the runtime is back online.
!!
!!   5. Restore the body from git history; rebuild with sub-record
!!      associate (Pattern 1) and state%agetracer access throughout.
!!
!! See docs/superpowers/specs/2026-05-10-state-migration-solute-discovery.md
!! Section 8 hazard #4 + the 2026-05-10 design spec D5.

module agetracer_mod
   use error_mod, only: fatalerr_collected
   implicit none
   private
   public :: AgeTracer

contains

   !> Stub: AgeTracer is currently inert. Reactivation requires the
   !! checklist in the file header above.
   subroutine AgeTracer(task, state)
      use swap_state_mod, only: swap_state_t
      implicit none

      integer,            intent(in) :: task
      type(swap_state_t), intent(in) :: state

      ! Silence unused-arg warnings (the dispatch signature is preserved
      ! for the day reactivation lands).
      associate (dummy_task => task, dummy_state => state)
      end associate

      call fatalerr_collected('AgeTracer', &
         'AgeTracer feature is currently inert. Reactivation requires '// &
         'agetracer_state_t and the checklist in src/solute/agetracer.f90.')
   end subroutine AgeTracer

end module agetracer_mod
