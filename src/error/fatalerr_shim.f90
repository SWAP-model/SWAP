!> TTutil `FatalERR` drop-in replacement.
!!
!! TTutil's upstream `FatalERR` (in `subprojects/ttutil/src/fatalerr.f90`)
!! ends with a bare `STOP`, which exits the process with status 0. That
!! masks failures from any parent that checks exit code (the pFUnit
!! harness, meson `test`, CI). It also blocks on `READ(*,*)` waiting for
!! `Press <Enter>`, which crashes with EOF when stdin is closed.
!!
!! This shim is a free subroutine with the same signature as TTutil's,
!! linked into every binary that depends on TTutil. Because object files
!! are processed before libraries by the linker, OUR `fatalerr_` symbol
!! wins over `libttutil`'s. All callers — TTutil's own internals
!! (`FOPENG`, `RDINIT`, `RDDATA`) AND legacy SWAP code in `src/` — route
!! through `fatalerr_collected`, which:
!!   - appends the error to the modern `error_collection_t` singleton
!!   - calls `error stop "fatal error(s) in swap input pipeline"` (exit 1)
!!
!! Net effect: legacy `fatalerr` calls produce the same error text as
!! before, but the process exits non-zero and the message goes through
!! `swap_log`.
!!
!! TO BE REMOVED when TTutil is phased out from the runtime
!! (legacy reader retirement umbrella spec; SS-11 closeout).
!! ADR 0018 documents the shim and its retirement gate.
SUBROUTINE FatalERR(MODULE, MESSAG)
   use error_mod, only: fatalerr_collected
   implicit none
   character(len=*), intent(in) :: MODULE, MESSAG
   call fatalerr_collected(trim(MODULE), trim(MESSAG))
END SUBROUTINE FatalERR
