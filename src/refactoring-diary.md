Here I'm writing about the progress in phase 2 of SWAP refactoring including some of the decisions I took.

# Modernization
SWAP has been developed for 50 years now and the coding practices changed. Even Fortran language itself changed a lot. Including the ways  

# Potential improvements
## Error messages
they are currently sent using the fatalerr (not sure where it comes from, probably ttutil). There is currently a way to use ISO_FORTRAN_ENV to make it cleaner:

```fortran
  ! --- drainage flux calculated according to hooghoudt or ernst
      if (dramet.eq.2) then
        if (shape .gt. small) dh = (gwldra-zbotdr(1)) / shape

  ! --- contributing layer below drains limited to 1/4 l
        zimp = max (basegw,zbotdr(1)-0.25*l(1))
        dbot = (zbotdr(1)-zimp)  
        if (dbot.lt.0.0d0) then
          messag = 'At the drainage section, the level of the'          &
     &       //' impervious layer is higher than the level of the'      &
     &       //' drain bottom. Adapt drain input!'
          call fatalerr ('Bocodrb',messag)
        endif
              
  ! --- no infiltration allowed
        if (dh.lt.1.0d-10) then
          qdrain(1) = 0.0d0
          return
        endif
```
After updete:

```fortran
  subroutine bocodrb(dh)
    real(8), intent(out) :: dh
    real(8) :: dbot, zimp
    
    ! ... calculations ...
    
    if (dbot < 0.0d0) then
      write(error_unit, '(A)') &
        'ERROR in bocodrb: At the drainage section, the impervious layer ' // &
        'is higher than the drain bottom. Adapt drain input!'
      error stop 1  ! Modern alternative to 'stop'
    endif
    
  end subroutine bocodrb

end module drainage_mod
```

2026-02-09 00:15AM
I just discovered that I could use a much simplier approach factoring in states into the codebase. It's partly intermediate, but will work and does not introduce overhead. I came up with an idea to use pointers instead of refactoring entire modules, and AI came with this answer that I should use the ASSOCIATE, which basically is like PARAMETER: at compile time the associated variables get switched with the state attributes I pass as aliases. Brilliant.

I should use ASSOCIATE construct

```fortran
subroutine bocodrb(dh, soil_state, drainage_config)
  use soil_state_type, only: soil_state_t
  use drainage_config_type, only: drainage_config_t
  
  type(soil_state_t), intent(inout) :: soil_state
  type(drainage_config_t), intent(in) :: drainage_config
  real(8), intent(out) :: dh
  
  ! Local variables
  real(8) :: zimp, dbot, pi, totres, x, fx, eqd
  integer :: i, lev
  parameter (pi=3.14159d0)
  
  ! === Use ASSOCIATE for clean aliasing ===
  associate( &
    gwl => soil_state%gwl, &
    qdrain => soil_state%qdrain, &
    zbotdr => drainage_config%zbotdr, &
    basegw => drainage_config%basegw, &
    dramet => drainage_config%method, &
    l => drainage_config%spacing, &
    nrlevs => drainage_config%nlevs &
  )
  
    ! === Code body exactly as before ===
    if (dramet.eq.2) then
      if (shape .gt. small) dh = (gwl-zbotdr(1)) / shape
      
      zimp = max(basegw, zbotdr(1)-0.25*l(1))
      dbot = (zbotdr(1)-zimp)  
      ! ... rest of code
    endif
  
  end associate
  
end subroutine bocodrb
```

Why ASSOCIATE is Better:

    ✅ No pointer overhead: Compiler can optimize better
    ✅ Scoped aliasing: Clear begin/end of aliased names
    ✅ Standard Fortran: Part of F2003/F2008 standard
    ✅ Intent clarity: Shows exactly which variables are used
    ✅ Safer: No memory aliasing concerns

Recommendation for For SWAP refactoring, go straight to Phase 3 (ASSOCIATE + direct):

    - Use ASSOCIATE for the 5-10 most frequently accessed variables per routine
    - Use state%variable for everything else
    - This keeps code readable while being modern and efficient
    - Avoids the pointer "trap" of appearing cleaner than it is

I also discovered a great linter. There are several things I need to get in order with the comments to the code, but ok. Everything at te right time. Although there are also other issues with that linter. Will need to look into it...

2026-02-10

Q:
I'm afraid a lot of SWAP code does just that, they update state instead of computing output
A:
You're describing a classic challenge with legacy scientific Fortran codes! SWAP (Soil Water Atmosphere Plant model) likely has decades of evolutionary development, and what you're encountering is extremely common in models from that era.
The Legacy Pattern

Many older Fortran models (pre-Fortran 90) used subroutines as "mini-programs" that directly manipulated global state through common blocks or argument mutation, rather than following modern functional decomposition principles. When these codes are modernized, the intent attributes expose these design patterns that were previously invisible.