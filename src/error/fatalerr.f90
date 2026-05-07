!> Canonical FatalERR implementation.
!!
!! Provides a free `subroutine FatalERR(MODULE, MESSAG)` matching the
!! TTutil-era signature so legacy free-call sites in `src/` resolve.
!! Routes through `fatalerr_collected`, which appends the error to
!! the modern `error_collection_t` singleton and aborts via
!! `error stop "fatal error(s) in swap input pipeline"` (exit 1).
!!
!! Pre-2026-05-07 this was a "shim" overriding TTutil's own
!! `fatalerr.f90` at link time. With TTutil retired (ADR 0023), this
!! is the canonical implementation; no shimming required.
SUBROUTINE FatalERR(MODULE, MESSAG)
   use error_mod, only: fatalerr_collected
   implicit none
   character(len=*), intent(in) :: MODULE, MESSAG
   call fatalerr_collected(trim(MODULE), trim(MESSAG))
END SUBROUTINE FatalERR
