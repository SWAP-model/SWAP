! File VersionID:
!   $Id: initialize.f90 374 2018-03-21 13:12:23Z heine003 $
! ----------------------------------------------------------------------
      Subroutine Initialize
! ----------------------------------------------------------------------
! [GR-IO 2026-05-25 Phase 6 Step 4] Bare-global zero-init body retired.
! Every legacy global previously zeroed here either:
!   - has been retired entirely (variables.f90 declarations are comments)
!   - has its canonical home on state%X (default-init via derived-type
!     component default) or config%X (default-init via type component
!     default + finalize defaults).
! The subroutine is preserved as a no-op stub because src/core/swap_mod.f90
! still calls `call Initialize` at startup. Once the swap_mod call retires
! (a separate cleanup), this file can be deleted.
! ----------------------------------------------------------------------
      implicit none
      return
      end subroutine

      Subroutine InitializeCrop
! ----------------------------------------------------------------------
! [GR-IO 2026-05-25 Phase 6 Step 4] Bare-global zero-init body retired.
! Crop-side state lives on state%crop%{common,wofost,grass,oxygen,irrigation};
! every legacy crop global has its derived-type default. Preserved as a
! no-op because src/crop/cropgrowth.f90 still calls InitializeCrop on
! new-rotation reset; can be deleted with that call site.
! ----------------------------------------------------------------------
      implicit none
      return
      end subroutine
