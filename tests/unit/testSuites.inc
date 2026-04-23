! pFUnit test-suite registry
!
! Empty at the rescue baseline (Phase 1 Task 2): every .pf file that
! existed at commit e256bc0 referenced at least one state/config module
! that was never actually committed, so the .pf files were removed and
! this registry was emptied to match.
!
! Add ADD_TEST_SUITE(<module_name>_suite) lines here as Phase 4
! reintroduces suites together with the state modules they exercise.
