! Add temporary test in your main code:
program test_reduceva
    use variables
    use et_mod, only: reduceva
    implicit none
    
    real(8) :: empreva_old, empreva_new
    
    ! Set up test conditions
    swredu = 1         ! Black model
    cofred = 3.5d0
    rsigni = 0.5d0
    peva = 0.4d0
    pond = 0.0d0
    ldwet = 5.0d0
    dt = 0.04167d0     ! 1 hour
    fldaystart = .false.
    
    ! Test daily (task=1)
    call reduceva(1, nrai=0.2d0)
    print *, 'Daily empreva:', empreva
    print *, 'Daily ldwet:', ldwet
    
    ! Test sub-daily (task=2)
    ldwet = 5.0d0      ! Reset
    call reduceva(2, nrai=0.2d0)
    print *, 'Sub-daily empreva:', empreva
    print *, 'Sub-daily ldwet:', ldwet
    
end program test_reduceva