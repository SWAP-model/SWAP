┌─────────────────────────────────────┐
│ Does function need 3+ fields from   │
│ the type?                           │
└──────────┬───────────────┬──────────┘
           │               │
          YES             NO
           │               │
           ▼               ▼
  ┌────────────────┐  ┌──────────────┐
  │ Pass the type  │  │ Extract the  │
  │ directly       │  │ field(s)     │
  └────────────────┘  └──────────────┘
           │               │
           ▼               ▼
  ┌────────────────┐  ┌──────────────┐
  │ Is function    │  │ Is function  │
  │ generic/       │  │ domain-      │
  │ reusable?      │  │ specific?    │
  └──┬─────────┬───┘  └──┬───────┬───┘
     │         │         │       │
    YES       NO        YES     NO
     │         │         │       │
     ▼         ▼         ▼       ▼
  Extract   Pass    Extract   Either
  values    type    values    works

Example:

```fortran
    !> @brief Public API - accepts high-level types
    subroutine DIVDRA(soil, drainage, time_info, qdrain, qdra)
        type(soil_state_t), intent(in) :: soil
        type(drainage_config_t), intent(in) :: drainage
        type(time_state_t), intent(in) :: time_info
        real(8), intent(inout) :: qdrain(:)
        real(8), intent(out) :: qdra(:,:)
        
        real(8) :: condsat_hor(MACP), condsat_ver(MACP)
        real(8) :: fac_aniso, wat_lev_av
        
        ! ✅ Pass whole soil to domain-specific helper
        call calc_conductivities(soil, condsat_hor, condsat_ver)
        
        ! ✅ Extract for simple conversions
        wat_lev_av = -1.0d0 * min(soil%gwl, 0.0d0)
        
        ! ✅ Pass extracted values to generic utility
        call Lev2Comp(soil%numnod, wat_lev_av, soil%dz, &
                     num_com_wat_lev, thick_cum, ...)
        
        ! ✅ Pass whole drainage to drainage-specific helper
        call distribute_discharge(soil, drainage, condsat_hor, &
                                 qdrain, qdra, ...)
    end subroutine DIVDRA

    !> @brief Domain helper - takes soil state (uses many fields)
    subroutine calc_conductivities(soil, condsat_hor, condsat_ver)
        type(soil_state_t), intent(in) :: soil  ! ✅ Pass whole state
        real(8), intent(out) :: condsat_hor(:), condsat_ver(:)
        
        ! Needs: numnod, layer, ksatfit, ksatexm, fluseksatexm, cofani
        ! That's 6 fields → definitely pass the type!
        
        where (soil%fluseksatexm(1:soil%numnod))
            condsat_hor = soil%ksatexm(soil%layer) * soil%cofani(soil%layer)
            condsat_ver = soil%ksatexm(soil%layer)
        elsewhere
            condsat_hor = soil%ksatfit(soil%layer) * soil%cofani(soil%layer)
            condsat_ver = soil%ksatfit(soil%layer)
        end where
    end subroutine
    
    !> @brief Generic utility - takes extracted values
    subroutine Lev2Comp(num_comp, level, thickness, &
                       num_com_2lev, thick_cum, ...)
        integer, intent(in) :: num_comp
        real(8), intent(in) :: level
        real(8), intent(in) :: thickness(:)  ! ✅ Extract, don't pass soil%
        
        ! This is GENERIC - could work on any array
        ! Should NOT depend on soil_state_t type
        ! Makes it reusable for other purposes
    end subroutine
    
    !> @brief Drainage helper - takes both types
    subroutine distribute_discharge(soil, drainage, condsat_hor, &
                                   qdrain, qdra, ...)
        type(soil_state_t), intent(in) :: soil       ! ✅ Pass type
        type(drainage_config_t), intent(in) :: drainage  ! ✅ Pass type
        real(8), intent(in) :: condsat_hor(:)        ! ✅ Derived data
        
        ! Needs soil%numnod, soil%dz, soil%gwl (3+ fields)
        ! Needs drainage%nrlevs, drainage%l, drainage%zbotdr (3+ fields)
        ! → Pass both types
    end subroutine

end module
```