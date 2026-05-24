module drainage_mod
   use error_mod, only: fatalerr_collected
!> Drainage Module - Lateral water flux calculations for subsurface drainage systems
!!
!! This module provides comprehensive routines for simulating lateral drainage and
!! infiltration fluxes in soil profiles with subsurface drainage systems. It supports
!! multiple drainage calculation methods and handles complex interactions between
!! groundwater, surface water, and drainage systems.
!!
!!### Calculation Methods
!!
!! The module supports three primary drainage calculation approaches:
!!
!! 1. **Table lookup (dramet=1)**: Direct interpolation from groundwater level vs flux tables
!! 2. **Hooghoudt/Ernst equations (dramet=2)**: Analytical solutions based on drainage theory
!! 3. **Resistance method (dramet=3)**: Flux calculated from drainage/infiltration resistances
!!
!!### Key Features
!!
!! - Multiple drainage levels (up to `nrlevs` systems)
!! - Layered soil profiles with anisotropic hydraulic conductivity
!! - Surface water - groundwater interaction
!! - Macropore rapid drainage pathways
!! - Interflow (surface runoff) simulation
!! - Dynamic wetted perimeter for open channels
!! - Storage capacity constraints for surface water systems
!! - Flux distribution over soil compartments using [[DIVDRA]]
!!
!!### Physical Concepts
!!
!! **Drainage resistance components:**
!!
!! - \(R_{ver}\): Vertical resistance from water table to drain level
!! - \(R_{hor}\): Horizontal resistance (radial flow to drains)
!! - \(R_{rad}\): Radial resistance near the drain
!! - \(R_{entry}\): Entry resistance into drain
!!
!! Total drainage flux: \(q = \frac{dh}{R_{tot}}\) where \(dh\) is head difference
!!
!! **Equivalent depth** for radial flow (Hooghoudt):
!!
!! \[ d_{eq} = \frac{\pi L}{8\left[\ln(L/r_0) + f(\pi d/L)\right]} \]
!!
!! where \(L\) is drain spacing, \(r_0\) is drain radius, and \(d\) is thickness
!! of saturated zone below drains.
!!
!!### Module Dependencies
!!
!! - **distribute_drainage**: Provides [[DIVDRA]] for flux distribution over compartments
!! - **variables**: Global state variables (gwl, qdrain, etc.)
!! - **params.fi**: Physical constants including `small`
!!
!!### Module History
!!
!!@history
!! - **2002-07**: Initial implementation
!! - **2004-08**: Added [[drainage]] main routine
!! - **2014-12**: Major update and refactoring
!! - **2017-09**: Added infiltration head limiting (P. van Walsum)
!! - **2026**: Converted to modern module structure
!!@endhistory
!!
!!### References
!!
!! - Hooghoudt, S.B. (1940): Bijdragen tot de kennis van eenige natuurkundige
!!   grootheden van den grond. Verslagen Landbouwkundige Onderzoekingen 46(14)
!! - Ernst, L.F. (1956): Calculation of the steady flow of groundwater in vertical
!!   cross sections. Netherlands Journal of Agricultural Science 4: 126-131
!! - van Dam, J.C. et al. (2008): SWAP version 3.2. Theory description and user manual.
!!   Wageningen University and Research Centre.
!!
   use distribute_drainage, only: DIVDRA
   use array_utils, only: afgen
   use swap_constants, only: small
   use swap_state_mod, only: swap_state_t
   implicit none

   public :: drainage
   public :: bocodrb, bocodre
   public :: drainage_init

contains

   !> Allocate and initialise all per-level arrays in state%drainage.
   !! Called from swap_main AFTER config_to_variables so that
   !! config-sourced geometry (e.g. wetper for dramet==2) can be
   !! seeded directly from config into state.
   !! ADR 0031 Phase 2 Task 5: drainl/wetper/ztopdislay/qdrd globals
   !! deleted; geometry is now seeded from config, flux arrays zeroed.
   subroutine drainage_init(state, config)
      use, intrinsic :: iso_fortran_env, only: real64
      use swap_state_mod, only: swap_state_t
      use swap_config_mod, only: swap_config_t
      ! [GR-BH Task 35] drainage_init is called AFTER CalcGrid (swap_mod.f90 line 108 vs 185),
      ! so state%mesh%numnod is already populated. numnod global deleted; use state%mesh%numnod.
      ! [GR-BH Task 37] nrlevs global deleted; use config%drain%nrlevs (authoritative single source).
      use swap_array_dimensions, only: MAOWL
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config
      integer :: nr

      nr = config%drain%nrlevs
      if (.not. allocated(state%drainage%qdrain))     allocate(state%drainage%qdrain(nr))
      if (.not. allocated(state%drainage%drainl))     allocate(state%drainage%drainl(nr))
      if (.not. allocated(state%drainage%wetper))     allocate(state%drainage%wetper(nr))
      if (.not. allocated(state%drainage%ztopdislay)) allocate(state%drainage%ztopdislay(nr))
      if (.not. allocated(state%drainage%qdra))       allocate(state%drainage%qdra(nr, state%mesh%numnod))
      if (.not. allocated(state%drainage%L))          allocate(state%drainage%L(nr))
      if (.not. allocated(state%drainage%zbotdr))     allocate(state%drainage%zbotdr(nr))
      if (.not. allocated(state%drainage%owltab))     allocate(state%drainage%owltab(nr, 2*MAOWL))

      ! Geometry arrays: start at zero; seed state%drainage%wetper(1) from config when
      ! dramet==2 (Hooghoudt/Ernst) — the only config-sourced geometry
      ! value.  drainl/ztopdislay are computed each timestep by bocodre;
      ! qdrd is computed by the secondary drainage block.
      state%drainage%drainl     = 0.0_real64
      state%drainage%wetper     = 0.0_real64
      state%drainage%ztopdislay = 0.0_real64
      state%drainage%qdrd       = 0.0_real64
      if (config%drain%dramet == 2) then
         state%drainage%wetper(1) = config%drain%wetper
      end if

      ! Flux arrays start at zero — computed each timestep.
      state%drainage%qdrain = 0.0_real64
      state%drainage%qdra   = 0.0_real64
   end subroutine drainage_init


   subroutine bocodrb(dh, state)
      !> Calculate drainage flux using Hooghoudt/Ernst or resistance methods
    !!
    !! This subroutine calculates the total drainage flux for a single drainage system
    !! or multiple levels based on the hydraulic head difference between groundwater
    !! level and drainage base level.
    !!
    !!### Calculation Methods
    !!
    !!#### Method 1: Table Lookup (dramet=1)
    !! Direct interpolation from pre-calculated gwl-flux tables.
    !!
    !!#### Method 2: Hooghoudt/Ernst Analytical (dramet=2)
    !!
    !! Uses classical drainage equations with five geometric configurations:
    !!
    !! - **ipos=1**: Homogeneous soil on impervious layer
    !! - **ipos=2**: Homogeneous soil, drain above impervious layer
    !! - **ipos=3**: Two-layer soil, drain at interface
    !! - **ipos=4**: Drain in bottom layer of two-layer system
    !! - **ipos=5**: Drain in top layer of two-layer system
    !!
    !! For each case, calculates equivalent depth and combines resistance components.
    !!
    !!#### Method 3: Resistance Method (dramet=3)
    !!
    !! Loops over all drainage levels. For each level:
    !!
    !! 1. Get surface water level from time series table
    !! 2. Calculate head difference \(dh = gwl - swl\)
    !! 3. Apply flux equation based on flow direction:
    !!    - Drainage: \(q = dh / R_{drain}\)
    !!    - Infiltration: \(q = dh / R_{infil}\)
    !!    - Interflow: \(q = C \cdot dh^E\) (power law)
    !!
    !!### Equivalent Depth Calculation
    !!
    !! For \(x = 2\pi d/L\):
    !!
    !! - If \(x > 0.5\): Use series expansion
    !!   \[ f(x) = \sum_{i=1,3,5} \frac{4e^{-2ix}}{i(1-e^{-2ix})} \]
    !!
    !! - If \(x < 0.5\): Use logarithmic approximation
    !!   \[ f(x) = \frac{\pi^2}{4x} + \ln\left(\frac{x}{2\pi}\right) \]
    !!
    !!### Special Features
    !!
    !! - Limits contributing layer below drains to maximum \(L/4\)
    !! - Handles dry channel conditions (sets drainage to zero)
    !! - Updates macropore drainage basis (`ZDraBas`) if enabled
    !! - Constrains infiltration head to water depth in channel (optional)
    !!
    !!@note
    !! For dramet=2, only infiltration is prevented (dh<0). For dramet=3,
    !! bidirectional flow is allowed based on resistance values.
    !!@endnote
    !!
      ! Residual `use variables`: drares/infres/shape are still bare globals; the rest of
      ! the symbols used in this routine come through the sub-record associate below.
      use variables, only: drares, infres, shape
      use array_utils, only: afgen

      real(8) dh
      type(swap_state_t), intent(inout) :: state

      integer :: i, lev
      real(8) :: zimp, dbot, pi, totres, x, fx, eqd, rver, rhor, rrad
      real(8) :: gwldra
      logical :: fldry
      character(len=200) :: messag

      parameter(pi=3.14159d0)

      associate (drai => state%drainage,    &
                 soil => state%soilwater,   &
                 time => state%timecontrol)

         gwldra = soil%gwl

         ! --- drainage flux according to hooghoudt or ernst
         if (drai%dramet .eq. 2) then
            if (shape .gt. small) dh = (gwldra - drai%zbotdr(1))/shape

            ! --- contributing layer below drains limited to 1/4 l
            zimp = max(drai%basegw, drai%zbotdr(1) - 0.25*drai%l(1))
            dbot = (drai%zbotdr(1) - zimp)
            if (dbot .lt. 0.0d0) then
               messag = 'At the drainage section, the level of the'          &
          &       //' impervious layer is higher than the level of the'      &
          &       //' drain bottom. Adapt drain input!'
               call fatalerr_collected('Bocodrb', messag)
            end if

            ! --- no infiltration allowed
            if (dh .lt. 1.0d-10) then
               drai%qdrain(1) = 0.0d0
               return
            end if

            ! --- case 1: homogeneous, on top of impervious layer
            if (drai%ipos .eq. 1) then
               totres = drai%l(1)*drai%l(1)/(4*drai%khtop*abs(dh)) + drai%entres
               drai%qdrain(1) = dh/totres

            ! --- case 2,3: in homogeneous profile or at interface of 2 layers
            elseif (drai%ipos .eq. 2 .or. drai%ipos .eq. 3) then

               ! --- equivalent depth
               x = 2*pi*dbot/drai%l(1)
               if (x .gt. 0.5d0) then
                  fx = 0.0d0
                  do 10 i = 1, 5, 2
                     fx = fx + (4*exp(-2*i*x))/(i*(1.0d0 - exp(-2*i*x)))
10                continue
                  eqd = pi*drai%l(1)/8/(log(drai%l(1)/drai%wetper(1)) + fx)
               else
                  if (x .lt. 1.0d-6) then
                     eqd = dbot
                  else
                     fx = pi**2/(4*x) + log(x/(2*pi))
                     eqd = pi*drai%l(1)/8/(log(drai%l(1)/drai%wetper(1)) + fx)
                  end if
               end if
               if (eqd .gt. dbot) eqd = dbot

               if (drai%ipos .eq. 2) then
                  totres = drai%l(1)*drai%l(1)/(8*drai%khtop*eqd + 4*drai%khtop*abs(dh)) + drai%entres
               elseif (drai%ipos .eq. 3) then
                  totres = drai%l(1)*drai%l(1)/(8*drai%khbot*eqd + 4*drai%khtop*abs(dh)) + drai%entres
               end if
               drai%qdrain(1) = dh/totres

            ! --- case 4: drain in bottom layer
            elseif (drai%ipos .eq. 4) then
               if (drai%zbotdr(1) .gt. drai%zintf) then
                  messag = 'At the drainage section, the level of the'        &
           &         //' impervious layer is higher than the level of the'    &
           &         //' drain bottom. Adapt drain input!'
                  call fatalerr_collected('bocodrb', messag)
               end if
               rver = max(gwldra - drai%zintf, 0.0d0)/drai%kvtop +                   &
        &             (min(drai%zintf, gwldra) - drai%zbotdr(1))/drai%kvbot
               rhor = drai%l(1)*drai%l(1)/(8*drai%khbot*dbot)
               rrad = drai%l(1)/(pi*dsqrt(drai%khbot*drai%kvbot))*log(dbot/drai%wetper(1))
               totres = rver + rhor + rrad + drai%entres
               drai%qdrain(1) = dh/totres

            ! --- case 5: drain in top layer
            elseif (drai%ipos .eq. 5) then
               if (drai%zbotdr(1) .lt. drai%zintf) then
                  messag = 'At the drainage section, the level of the'        &
           &         //' impervious layer is higher than the level of the'    &
           &         //' drain bottom. Adapt drain input!'
                  call fatalerr_collected('bocodrb', messag)
               end if
               rver = (gwldra - drai%zbotdr(1))/drai%kvtop
               rhor = drai%l(1)*drai%l(1)/(8*drai%khtop*(drai%zbotdr(1) - drai%zintf) +    &
        &                                  8*drai%khbot*(drai%zintf - zimp))
               rrad = drai%l(1)/(pi*dsqrt(drai%khtop*drai%kvtop))*log((drai%geofac*        &
        &              (drai%zbotdr(1) - drai%zintf))/drai%wetper(1))
               totres = rver + rhor + rrad + drai%entres
               drai%qdrain(1) = dh/totres
            end if

         ! --- drainage flux from drain/infiltration resistance (per level)
         elseif (drai%dramet .eq. 3) then
            do lev = 1, drai%nrlevs
               fldry = .false.

               x = afgen(drai%owltab(lev, 1:2*drai%nowltab(lev)), 2*drai%nowltab(lev), &
                         time%t1900 + time%dt - 1.d0)
               if ((x - drai%zbotdr(lev)) .lt. 1.0d-3) fldry = .true.
               dh = gwldra - x
               if (fldry) dh = gwldra - drai%zbotdr(lev)
               ! Macropore rapid-drainage block deleted (ADR 0040).

               ! --- drainage
               if (dh .ge. 0.0d0) then
                  ! interflow as power function
                  if ((lev .eq. drai%nrlevs) .and. (drai%swnrsrf .eq. 1)) then
                     drai%qdrain(lev) = drai%cofintfl*dh**drai%expintfl
                  else
                     drai%qdrain(lev) = dh/drares(lev)
                     if (drai%swallo(lev) .eq. 2) drai%qdrain(lev) = 0.0d0
                  end if

               ! --- infiltration
               else
                  ! Limit infiltration head (-dh) to channel waterdepth (KRO 2017-09-07)
                  if (drai%swdtyp(lev) .eq. 2) then
                     if (drai%swliminf .eq. 1) dh = max(dh, (drai%zbotdr(lev) - x))
                  end if
                  drai%qdrain(lev) = dh/infres(lev)
                  if (drai%swallo(lev) .eq. 3 .or. fldry) drai%qdrain(lev) = 0.0d0
               end if
            end do

         ! --- drainage flux from gwlevel-flux table
         elseif (drai%dramet .eq. 1) then
            drai%qdrain(1) = afgen(drai%qdrtab, 50, abs(gwldra))
         end if

      end associate
    end subroutine bocodrb

    subroutine drainage(state)
    !> Main drainage orchestration routine
    !!
    !! This subroutine coordinates all drainage-related calculations for each time step,
    !! including initialization, flux calculation, distribution, and flux accounting.
    !!
    !!### Algorithm Workflow
    !!
    !! **1. Initialization (first call only)**
    !! - Set macropore drainage basis level
    !! - Extract from drain bottom or surface water table
    !!
    !! **2. Reset flux accumulators**
    !! - Zero intermediate fluxes if `flzerointr = .true.`
    !! - Zero cumulative fluxes if `flzerocumu = .true.`
    !!
    !! **3. Handle deep groundwater**
    !! - Return with zero drainage if gwl > 998 cm (below profile)
    !!
    !! **4. Calculate total drainage rates**
    !! - Call [[bocodrb]] to compute `qdrain(level)` for all levels
    !! - Returns head difference `dh` for optional use
    !!
    !! **5. Distribute fluxes over compartments**
    !!
    !! Two distribution options controlled by `swdivd`:
    !!
    !! - **swdivd = 0**: All flux through bottom compartment only
    !! - **swdivd = 1**: Distributed using [[DIVDRA]] based on transmissivity
    !!
    !! **6. Adjust discharge layer tops (optional)**
    !!
    !! If `swdislay = 1` or `2`, recalculate top of discharge layers:
    !!
    !! \[ z_{top} = f \cdot gwl + (1-f) \cdot (gwl - dh) \]
    !!
    !! Then redistribute fluxes to exclude compartments above new top.
    !!
    !! **7. Sum total drainage**
    !! - Calculate `qdrtot` as sum of all level fluxes
    !!
    !!### Flux Distribution Options
    !!
    !! | swdivd | Description | Vertical profile |
    !! |--------|-------------|------------------|
    !! | 0 | Bottom only | All flux at `numnod` |
    !! | 1 | Transmissivity-weighted | Distributed over discharge layer |
    !!
    !!### Discharge Layer Adjustment
    !!
    !! | swdislay | swtopdislay | Behavior |
    !! |----------|-------------|----------|
    !! | 0 | - | Fixed discharge layer |
    !! | 1 | 1 | User-specified top level |
    !! | 2 | 1 | Dynamic top based on gwl and dh |
    !!
    !!@note
    !! The `flInitDraBas` flag ensures macropore initialization happens only once
    !! per simulation. The routine returns early after initialization.
    !!@endnote
    !!

      use swap_array_dimensions, only: madr
      use array_utils,           only: afgen

      type(swap_state_t), intent(inout) :: state

      integer :: node, level
      real(8) :: zCum, zTopDisLay(madr), difzTopDisLay(madr), ratio, ratiodz, sumqdr(madr), dh
      integer :: nodeTopDisLay(madr)
      character(len=33) :: messag

      associate (mesh => state%mesh,         &
                 drai => state%drainage,     &
                 soil => state%soilwater,    &
                 surf => state%surfacewater, &
                 time => state%timecontrol)

         ! Allocate per-level state arrays if not yet done (guard for
         ! fldrain path where surfacewater_init may not have been called).
         if (.not. allocated(surf%cqdrain)) then
            allocate(surf%cqdrain(drai%nrlevs));     surf%cqdrain    = 0.0d0
         end if
         if (.not. allocated(surf%cqdrainin)) then
            allocate(surf%cqdrainin(drai%nrlevs));   surf%cqdrainin  = 0.0d0
         end if
         if (.not. allocated(surf%cqdrainout)) then
            allocate(surf%cqdrainout(drai%nrlevs));  surf%cqdrainout = 0.0d0
         end if
         if (.not. allocated(surf%inqdra)) then
            allocate(surf%inqdra(drai%nrlevs, mesh%numnod));      surf%inqdra     = 0.0d0
         end if
         if (.not. allocated(surf%inqdra_in)) then
            allocate(surf%inqdra_in(drai%nrlevs, mesh%numnod));   surf%inqdra_in  = 0.0d0
         end if
         if (.not. allocated(surf%inqdra_out)) then
            allocate(surf%inqdra_out(drai%nrlevs, mesh%numnod));  surf%inqdra_out = 0.0d0
         end if

         ! Macropore initialisation: drainage basis for rapid drainage.
         if (surf%flInitDraBas) then
            if (drai%NumLevRapDra .gt. drai%nrlevs) then
               messag = ' NUMLEVRAPDRA greater then NRLEVS'
               call fatalerr_collected('MacroRead', messag)
            end if

            if (drai%dramet .lt. 3) then
               surf%ZDraBas = drai%zbotdr(1)
            else
               if (drai%swdtyp(drai%NumLevRapDra) .eq. 1) then
                  surf%ZDraBas = drai%zbotdr(drai%NumLevRapDra)
               else
                  surf%ZDraBas = afgen(drai%owltab(drai%NumLevRapDra, 1:2*drai%nowltab(drai%NumLevRapDra)), &
                                       2*drai%nowltab(drai%NumLevRapDra), time%t1900)
               end if
            end if

            surf%flInitDraBas = .false.

         else  ! .not. flInitDraBas — normal timestep path

            ! Reset intermediate surface-water and drainage fluxes.
            if (time%flZeroIntr) call surf%reset_intermediate()

            ! Reset cumulative drainage fluxes (cqdra/cqdrain*); reservoir-owned
            ! fields are reset by SurfaceWater(2) (ADR 0042).
            if (time%flZeroCumu) call surf%reset_cumulative_drainage()

            ! Skip when groundwater level is below the profile.
            if (soil%gwl .gt. 998.0d0) then
               do level = 1, drai%nrlevs
                  drai%qdrain(level) = 0.0d0
               end do
            else

               ! Total drainage rate per level.
               call bocodrb(dh, state)

               ! Partition the drainage flux over soil compartments.
               if (drai%swdivd .eq. 1) then
                  call divdra(mesh%numnod, drai%nrlevs, mesh%dz, soil%ksatfit, soil%ksatexm, &
                              soil%fluseksatexm, mesh%layer, soil%cofani, soil%gwl,          &
                              drai%L, drai%qdrain, drai%qdra,                                &
                              drai%swdivdinf, drai%swnrsrf, drai%swtopnrsrf, drai%zbotdr,    &
                              time%dt, drai%FacDpthInf, drai%owltab, drai%nowltab, time%t1900)

                  ! Redistribute qdrain with the new top boundary for discharge layers.
                  if (drai%swdislay .eq. 2) then
                     do level = 1, drai%nrlevs
                        if (drai%swtopdislay(level) .eq. 1) then
                           zTopDisLay(level) = drai%fTopDisLay(level)*soil%gwl +     &
                                               (1.0d0 - drai%fTopDisLay(level))*(soil%gwl - dh)
                        end if
                     end do
                  end if
                  if (drai%swdislay .eq. 1 .or. drai%swdislay .eq. 2) then
                     do level = 1, drai%nrlevs
                        if (drai%swtopdislay(level) .eq. 1) then
                           ! Find node number of the new top of the discharge layer.
                           nodeTopDisLay(level) = 1
                           zCum = -mesh%dz(1)
                           do while (zTopDisLay(level) .lt. zCum)
                              nodeTopDisLay(level) = nodeTopDisLay(level) + 1
                              zCum = zCum - mesh%dz(nodeTopDisLay(level))
                           end do
                           ! Saturated fraction of the partial compartment at the water level.
                           difzTopDisLay(level) = zTopDisLay(level) - zCum
                           ratiodz       = difzTopDisLay(level)/mesh%dz(nodeTopDisLay(level))
                           sumqdr(level) = ratiodz*drai%qdra(level, nodeTopDisLay(level))
                           do node = nodeTopDisLay(level) + 1, mesh%numnod
                              sumqdr(level) = sumqdr(level) + drai%qdra(level, node)
                           end do
                           if (dabs(sumqdr(level)) .lt. 1.0d-8) then
                              ratio = 1.0d0
                           else
                              ratio = drai%qdrain(level)/sumqdr(level)
                           end if
                           ! Redistribute drain-water fluxes.
                           do node = 1, nodeTopDisLay(level) - 1
                              drai%qdra(level, node) = 0.0d0
                           end do
                           drai%qdra(level, nodeTopDisLay(level)) =                    &
                              drai%qdra(level, nodeTopDisLay(level))*ratio*ratiodz
                           do node = nodeTopDisLay(level) + 1, mesh%numnod
                              drai%qdra(level, node) = drai%qdra(level, node)*ratio
                           end do
                        end if
                     end do
                  end if
               else
                  ! Drainage flux through lowest compartment only.
                  do level = 1, drai%nrlevs
                     do node = 1, mesh%numnod - 1
                        drai%qdra(level, node) = 0.0d0
                     end do
                     drai%qdra(level, mesh%numnod) = drai%qdrain(level)
                  end do
               end if

               surf%qdrtot = 0.0d0
               do level = 1, drai%nrlevs
                  surf%qdrtot = surf%qdrtot + drai%qdrain(level)
               end do

            end if  ! gwl > 998 skip

         end if  ! flInitDraBas vs. normal timestep

      end associate
    end subroutine drainage

            subroutine bocodre(dh, state)
               !> Calculate drainage/infiltration with surface water management
    !!
    !! This subroutine handles drainage calculations when surface water management
    !! is active (`swsrf >= 2`), including interactions with secondary drainage systems
    !! and storage capacity constraints.
    !!
    !!### Key Features
    !!
    !! **Dynamic wetted perimeter for open channels:**
    !!
    !! For trapezoidal channels:
    !! \[ P_{wet} = w + 2\sqrt{h^2 + (h/s)^2} \]
    !! where \(w\) is bottom width, \(h\) is water depth, \(s\) is side slope.
    !!
    !! **Drainage base level determination:**
    !! - If only gwl above bottom: drain bottom is base
    !! - If surface water above bottom: surface level is base
    !! - Updates macropore drainage basis
    !!
    !! **Resistance calculation:**
    !! - Drainage (dh>0): Uses `rdrain` and `rentry`
    !! - Infiltration (dh<0): Uses `rinfi` and `rexit`
    !! - Surface drainage: Dynamic resistance based on head
    !!
    !! **Interflow power law:**
    !! \[ q = C_{intfl} \cdot dh^{E_{intfl}} \]
    !!
    !!### Storage Constraint Handling
    !!
    !! For secondary systems (`swsec=2`, `swsrf>=2`):
    !!
    !! 1. Calculate potential storage change: \(\Delta V = (q_{tot} + q_{cap}) \cdot \Delta t\)
    !! 2. Check if storage becomes negative
    !! 3. If yes, reduce all secondary level fluxes proportionally:
    !!    \[ q_{adj} = q_{orig} \cdot \frac{q_{max}}{q_{tot}} \]
    !!    where \(q_{max} = -(Storage + q_{cap} \cdot \Delta t) / \Delta t\)
    !!
    !!### Surface Drainage (Vacuum Cleaner)
    !!
    !! When `swnrsrf=1` for top level:
    !! - Resistance decreases with increasing head
    !! - Minimum resistance constrained by `rsurfshallow`
    !! - \(R_d = R_{deep} - dh\)
    !!
    !!@note
    !! The `imper` variable tracks the current surface water management period.
    !! Periods are defined by `impend` array.
    !!@endnote
    !!
    !!@warning
    !! Uses GOTO statement for management period search (line 800).
    !! This is legacy code structure that should be refactored to DO/IF constructs.
    !!@endwarning
    !!
    !!@note
    !! ----------------------------------------------------------------------
    !!     Date               : 5/5/2001
    !!     Purpose            :
    !! --- Calculate drainage/infiltration fluxes for all levels: qdrain(level).
    !! --- Summate fluxes to secondary system: qdrd
    !! --- Given present storage (swst), qdrd and wscap(imper), check if the
    !! --- system falls dry. If so reduce drainage fluxes proportionally in
    !! --- such a way that total drainage flux equals available amount
    !! --- ( = swst + wscap).
    !!     Subroutines called : -
    !!     Functions called   : -
    !!     File usage         : -  Error handling
    !! ----------------------------------------------------------------------
    !!@endnote
  ! SS-DRST Phase 2 Task 4: qdrain removed from use-variables; written via state%drainage%qdrain.
  ! SS-SWC Phase 2 S-2.8: gwl and pond removed from use-list; read from state%soilwater.
  ! [SS-TC TC-14] t1900, dt read via state%timecontrol (ADR 0041)
  ! GR-BH Task 28: zbotdr, l, nrlevs, swnrsrf off variables → state%drainage%X.
  use variables, only: &   ! [SS-GR-FINAL B10] residuals — all DEFERRED
     ! swsec/swsrf retired (→state%cfg%surface_water%X); pilot
     ! [GR-DRA 2026-05-23] nrpri/nmper/impend/wscap retired — aliased from state%surfacewater below.
     ! DEFERRED: taludr/widthr/rdrain — drain geometry; config; Phase C3
     taludr, widthr, rdrain, &
     ! pondmx retired (→state%surfacewater%pondmx)

     ! [GR-DRA 2026-05-23] swdtyp retired — aliased from state%drainage below.
     ! DEFERRED: wlp/rinfi/rentry/rexit/gwlinf — surface water config; Phase C3
     wlp, rinfi, rentry, rexit, gwlinf, &
     ! DEFERRED: rsurfdeep/rsurfshallow — surface resistance config; Phase C3
     rsurfdeep, rsurfshallow
     ! [GR-DRA 2026-05-23] cofintfl/expintfl/NumLevRapdra retired — aliased from state%drainage below.
     ! [SS-GR-CROPRT A2] FlMacropore dropped — retired (ADR 0040)


! --- global
               real(8) dh
               type(swap_state_t), intent(inout) :: state

! --- local
               integer level, imper
               real(8) qdrdm, qdratio, swdepth, swexbrd, dvmax, swstmax, wl, rd, re
               character(len=200) messag
! Removed save statement for imper to avoid issues in parallel runs
! ----------------------------------------------------------------------

               ! SS-DRST Phase 2 Task 4: qdrain alias points directly to state; no legacy global written.
               ! SS-SWC Phase 2 S-2.8: gwl/pond read from state%soilwater via ASSOCIATE aliases.
               ! GR-BH Task 28: zbotdr/l/nrlevs/swnrsrf aliased from state%drainage.
               associate( &
                  wls     => state%surfacewater%wls,    &
                  swst    => state%surfacewater%swst,   &
                  ZDraBas => state%surfacewater%ZDraBas, &
                  qdrain  => state%drainage%qdrain,     &
                  zbotdr  => state%drainage%zbotdr,     &  ! GR-BH Task 28
                  l       => state%drainage%L,          &  ! GR-BH Task 28
                  nrlevs  => state%drainage%nrlevs,     &  ! GR-BH Task 28
                  swnrsrf => state%drainage%swnrsrf,    &  ! GR-BH Task 28
                  gwl     => state%soilwater%gwl,       &
                  pond    => state%soilwater%pond,      &
                  cofintfl => state%drainage%cofintfl,  &  ! [GR-DRA 2026-05-23]
                  expintfl => state%drainage%expintfl,  &  ! [GR-DRA 2026-05-23]
                  swdtyp   => state%drainage%swdtyp,    &  ! [GR-DRA 2026-05-23]
                  nrpri  => state%surfacewater%nrpri,   &
                  nmper  => state%surfacewater%nmper,   &
                  impend => state%surfacewater%impend,  &
                  wscap  => state%surfacewater%wscap    )

! --- Spec D7: zero drainage when groundwater is dry.
!     Was in SurfaceWater(2) in legacy code; relocated here per ADR 0030
!     since qdrain is drainage-owned.
               if (gwl .gt. 998.0d0) then
                  do level = 1, nrlevs
                     qdrain(level) = 0.0d0
                  end do
                  return
               end if

! --- summate fluxes for use by swballev and swlevbal
               state%drainage%qdrd = 0.0d0

               do 500 level = 1, nrlevs

! --- surface water level
                  if (state%cfg%surface_water%swsrf .ge. 2) then
                     if (level .gt. nrpri) then
                        wl = wls
                     else
                        wl = wlp
                     end if
                  end if

! --- drainage fluxes are set to zero if both groundwater level and surface
!     water level are above ponding sill (so the nonzero drainage flux
!     is only computed if either the gwl or the wl is below pondmx)
                  if (wl .lt. state%surfacewater%pondmx .or. gwl .lt. state%surfacewater%pondmx) then

! --- channel is active medium if either groundwater or surface water
!     level is above channel bottom
                     if (gwl .gt. (zbotdr(level) + 0.001d0) .or.                       &
                &        wl .gt. (zbotdr(level) + 0.001d0)) then
                        if (wl .le. (zbotdr(level) + 0.001d0) .or. state%cfg%surface_water%swsrf .eq. 1) then

! --- only groundw. level above channel bottom; bottom is dr. base
                           state%drainage%drainl(level) = zbotdr(level)

! --- wetted perimeter only computed for open channels
!     (for drains it is input)
                           if (swdtyp(level) .eq. 0) then
                              state%drainage%wetper(level) = widthr(level)
                           end if
                        else

! --- surface water level above channel bottom
                           state%drainage%drainl(level) = wl
                           if (swdtyp(level) .eq. 0) then
                              swdepth = wl - zbotdr(level)
                              swexbrd = (wl - zbotdr(level))/taludr(level)
                              state%drainage%wetper(level) = widthr(level) +           &
                   &            2*dsqrt(swdepth**2 + swexbrd**2)
                           end if
                        end if

! --- drainage flux (cm/d)
                        ! calculate head difference
                        dh = gwl - state%drainage%drainl(level)
                        if (gwl .gt. -0.1d0) dh = dh + pond
                        if (dh .lt. 0.0d0 .and. gwl .lt. gwlinf(level)) then
                           dh = gwlinf(level) - state%drainage%drainl(level)
                        end if
                        ! interflow flux calculated by a power function,
                        if ((level .eq. nrlevs) .and. (swnrsrf .eq. 2)) then
                           qdrain(level) = cofintfl*dh**expintfl
                        else
                           if (dh .gt. 0.0d0) then
                              rd = rdrain(level)
                              re = rentry(level)
! ---           surface drainage (vaccuum cleaner)
                              if ((level .eq. nrlevs) .and. (swnrsrf .eq. 1)) then
                                 rd = rsurfdeep - dh
                                 rd = max(rd, rsurfshallow)
                              end if
                           else
                              rd = rinfi(level)
                              re = rexit(level)
                           end if
                           if (swdtyp(level) .eq. 0) then
                              qdrain(level) = dh/((rd + re*l(level)/state%drainage%wetper(level)))
                           else
                              qdrain(level) = dh/rd
                           end if
                        end if
                     else
                        if (swdtyp(level) .eq. 0) then
                           state%drainage%wetper(level) = 0.0d0
                        end if
                        dh = 0.0d0
                        qdrain(level) = 0.0d0

!   - for determining drainage basis for rapid drainage through macropores
                        state%drainage%drainl(level) = zbotdr(level)

                     end if
!
                  else
                     qdrain(level) = 0.0d0

!   - for determining drainage basis for rapid drainage through macropores
                     state%drainage%drainl(level) = wl

                  end if
!
                  if (state%cfg%surface_water%swsrf .ge. 2 .and. level .gt. nrpri) then

! --- qdrd is total flux to or from secondary system
                     state%drainage%qdrd = state%drainage%qdrd + qdrain(level)
                  end if

500               continue

!   [SS-GR-CROPRT A2] rapid drainage macropore basis block dropped (ADR 0040; FlMacropore always .false.)

! ----------------------------------------------------------------------
! --- check for system falling dry (only for state%cfg%surface_water%swsec = 2):
                  if (state%cfg%surface_water%swsec .eq. 1) return
                  if (state%cfg%surface_water%swsrf .eq. 1) then
                     do 10 level = 1, nrlevs
                        if (qdrain(level) .lt. 0.0d0) then
                           qdrain(level) = 0.0d0
                        end if
10                      continue
                        elseif (state%cfg%surface_water%swsrf .ge. 2) then

! --- determine which management period the model is in:
                        imper = 0
800                     imper = imper + 1

! ---   Error handling
                        if (imper .gt. nmper) then
                           messag = 'sw-management periods(IMPER), more than defined'
                           call fatalerr_collected('Bocodre', messag)
                        end if
                        if (state%timecontrol%t1900 - 1.d0 + 0.1d-10 .gt. impend(imper)) goto 800

! ---   determine whether the system will become empty
                        dvmax = (state%drainage%qdrd + wscap(imper))*state%timecontrol%dt
                        swstmax = swst + dvmax

                        if (swstmax .lt. 0.0d0) then
! ---     storage decreases to below zero, then the surface water system
!         falls dry; make the total infiltration exactly equal to the
!         available amount:
                           qdrdm = -(swst + wscap(imper)*state%timecontrol%dt)/state%timecontrol%dt
                           qdratio = qdrdm/state%drainage%qdrd

! ---     Error handling
                           if (qdratio .gt. 1.0d0 .or. qdratio .lt. 0.0d0) then
                              messag = 'sw-management error with storage (qdratio)'
                              call fatalerr_collected('Bocodre', messag)
                           end if

                           do 820 level = 1 + NRPRI, nrlevs
                              qdrain(level) = qdrain(level)*qdratio
820                           continue
                              state%drainage%qdrd = qdrdm
                              end if
                           end if

                           end associate
                           end subroutine bocodre

                           end module drainage_mod
