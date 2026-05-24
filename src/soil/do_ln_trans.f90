!> @file src/soil/do_ln_trans.f90
!! Single-parameter module extracted from `sptabulated.f90` on 2026-05-24.
!!
!! `do_ln_trans` is the global log-transform flag for hydraulic
!! conductivity values. It is consumed by the **live** Mualem-van
!! Genuchten compute path (`src/soil/soilhydraulics.f90`,
!! `src/utils/soilhydraulicsutils.f90`, `src/crop/oxygenstress.f90`)
!! independently of the tabulated soil-physics path.
!!
!! The rest of the original `sptabulated.f90` (TSPACK library +
!! `EvalTabulatedFunction` + `PreProcTabulatedFunction`) has been
!! moved to `src/soil/dormant/sptabulated.f90` — see that file's
!! header for the dormancy / reactivation notes.
module doln
   logical, parameter :: do_ln_trans = .true.
end module doln
