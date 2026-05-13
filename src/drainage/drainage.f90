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
      use variables, only: nrlevs, numnod, MAOWL
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config

      if (.not. allocated(state%drainage%qdrain))     allocate(state%drainage%qdrain(nrlevs))
      if (.not. allocated(state%drainage%drainl))     allocate(state%drainage%drainl(nrlevs))
      if (.not. allocated(state%drainage%wetper))     allocate(state%drainage%wetper(nrlevs))
      if (.not. allocated(state%drainage%ztopdislay)) allocate(state%drainage%ztopdislay(nrlevs))
      if (.not. allocated(state%drainage%qdra))       allocate(state%drainage%qdra(nrlevs, numnod))
      if (.not. allocated(state%drainage%L))          allocate(state%drainage%L(nrlevs))
      if (.not. allocated(state%drainage%zbotdr))     allocate(state%drainage%zbotdr(nrlevs))
      if (.not. allocated(state%drainage%owltab))     allocate(state%drainage%owltab(nrlevs, 2*MAOWL))

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
      ! ADR 0031 Phase 2 Task 5: wetper removed from use-list; read from state%drainage%wetper(1).
      ! SS-SWC Phase 2 S-2.8: gwl removed from use-list; read from state%soilwater%gwl.
      use variables, only: dramet,zbotdr,basegw,l,ipos,khtop,khbot,kvtop,kvbot,entres,zintf,geofac,swdtyp,      &
owltab,swallo,drares,infres,qdrtab,nrlevs,swnrsrf,cofintfl,expintfl,shape,FlMacropore,NumLevRapDra,swliminf,nowltab
      ! [SS-TC TC-14] t1900, dt read via state%timecontrol (ADR 0041)
      use array_utils, only: afgen

      ! --- global
      real(8) dh
      type(swap_state_t), intent(inout) :: state

      ! ----------------------------------------------------------------------
      ! --- local
      integer i, lev

      real(8) zimp, dbot, pi, totres, x, fx, eqd, rver, rhor, rrad
      real(8) gwldra   !, temptab(2*maowl)

      logical fldry

      parameter(pi=3.14159d0)

      character(len=200) messag
      ! ----------------------------------------------------------------------

      ! SS-DRST Phase 2 Task 4: qdrain alias points directly to state; no legacy global written.
      associate(ZDraBas => state%surfacewater%ZDraBas, &
                qdrain  => state%drainage%qdrain)

      ! SS-SWC Phase 2 S-2.8: gwl read from state%soilwater
      gwldra = state%soilwater%gwl

      ! --- drainage flux calculated according to hooghoudt or ernst
      if (dramet .eq. 2) then
         if (shape .gt. small) dh = (gwldra - zbotdr(1))/shape

         ! --- contributing layer below drains limited to 1/4 l
         zimp = max(basegw, zbotdr(1) - 0.25*l(1))
         dbot = (zbotdr(1) - zimp)
         if (dbot .lt. 0.0d0) then
            messag = 'At the drainage section, the level of the'          &
       &       //' impervious layer is higher than the level of the'      &
       &       //' drain bottom. Adapt drain input!'
            call fatalerr_collected('Bocodrb', messag)
         end if

         ! --- no infiltration allowed
         if (dh .lt. 1.0d-10) then
            qdrain(1) = 0.0d0
            return
         end if

         ! --- case 1: homogeneous, on top of impervious layer
         if (ipos .eq. 1) then

            ! --- calculation of drainage resistance and drainage flux
            totres = l(1)*l(1)/(4*khtop*abs(dh)) + entres
            qdrain(1) = dh/totres

            ! --- case 2,3: in homogeneous profile or at interface of 2 layers
         elseif (ipos .eq. 2 .or. ipos .eq. 3) then

            ! --- calculation of equivalent depth
            x = 2*pi*dbot/l(1)
            if (x .gt. 0.5d0) then
               fx = 0.0d0
               do 10 i = 1, 5, 2
                  fx = fx + (4*exp(-2*i*x))/(i*(1.0d0 - exp(-2*i*x)))
10                continue
                  eqd = pi*l(1)/8/(log(l(1)/state%drainage%wetper(1)) + fx)
                  else
                  if (x .lt. 1.0d-6) then
                     eqd = dbot
                  else
                     fx = pi**2/(4*x) + log(x/(2*pi))
                     eqd = pi*l(1)/8/(log(l(1)/state%drainage%wetper(1)) + fx)
                  end if
                  end if
                  if (eqd .gt. dbot) eqd = dbot

                  ! --- calculation of drainage resistance & drainage flux
                  if (ipos .eq. 2) then
                     totres = l(1)*l(1)/(8*khtop*eqd + 4*khtop*abs(dh)) + entres
                  elseif (ipos .eq. 3) then
                     totres = l(1)*l(1)/(8*khbot*eqd + 4*khtop*abs(dh)) + entres
                  end if
                  qdrain(1) = dh/totres

                  ! --- case 4: drain in bottom layer
               elseif (ipos .eq. 4) then
                  if (zbotdr(1) .gt. zintf) then
                     messag = 'At the drainage section, the level of the'        &
              &         //' impervious layer is higher than the level of the'    &
              &         //' drain bottom. Adapt drain input!'
                     call fatalerr_collected('bocodrb', messag)
                  end if
                  rver = max(gwldra - zintf, 0.0d0)/kvtop +                   &
           &                (min(zintf, gwldra) - zbotdr(1))/kvbot
                  rhor = l(1)*l(1)/(8*khbot*dbot)
                  rrad = l(1)/(pi*dsqrt(khbot*kvbot))*log(dbot/state%drainage%wetper(1))
                  totres = rver + rhor + rrad + entres
                  qdrain(1) = dh/totres

                  ! --- case 5 : drain in top layer
               elseif (ipos .eq. 5) then
                  if (zbotdr(1) .lt. zintf) then
                     messag = 'At the drainage section, the level of the'        &
              &         //' impervious layer is higher than the level of the'    &
              &         //' drain bottom. Adapt drain input!'
                     call fatalerr_collected('bocodrb', messag)
                  end if
                  rver = (gwldra - zbotdr(1))/kvtop
                  rhor = l(1)*l(1)/(8*khtop*(zbotdr(1) - zintf) +               &
             &                             8*khbot*(zintf - zimp))
                  rrad = l(1)/(pi*dsqrt(khtop*kvtop))*log((geofac*               &
             &                (zbotdr(1) - zintf))/state%drainage%wetper(1))
                  totres = rver + rhor + rrad + entres
                  qdrain(1) = dh/totres
               end if

               ! --- drainage flux calc. using given drainage/infiltration resistance
            elseif (dramet .eq. 3) then

               do lev = 1, nrlevs
                  fldry = .false.

                  ! ---     first copy to 1-dimensional table temptab
                  !          do i = 1,2*nowltab(lev)
                  !          do i = 1,2*maowl
                  !            temptab(i) = owltab(lev,i)
                  !          end do
                  !          x = afgen(temptab,2*nowltab(lev),t1900+dt-1.d0)
                  !          x = afgen(temptab,2*maowl,t1900+dt-1.d0)
                  x = afgen(owltab(lev, 1:2*nowltab(lev)), 2*nowltab(lev), state%timecontrol%t1900 + state%timecontrol%dt - 1.d0)
                  !          x = afgen(owltab(lev,1:2*maowl),2*maowl,t1900+dt-1.d0)
                  if ((x - zbotdr(lev)) .lt. 1.0d-3) then
                     fldry = .true.
                  end if
                  dh = gwldra - x
                  if (fldry) dh = gwldra - zbotdr(lev)
                  ! ---     drainage basis for rapid drainage through macropores
                  if (FlMacropore .and. lev .eq. NumLevRapdra .and. swdtyp(lev) .eq. 2) then
                     ZDraBas = dmax1(x, zbotdr(NumLevRapDra))
                  end if

                  ! ---     drainage
                  if (dh .ge. 0.0d0) then

                     ! ---       interflow flux calculated by a power function
                     if ((lev .eq. nrlevs) .and. (swnrsrf .eq. 1)) then
                        qdrain(lev) = cofintfl*dh**expintfl
                     else
                        qdrain(lev) = dh/drares(lev)
                        if (swallo(lev) .eq. 2) qdrain(lev) = 0.0d0
                     end if

                     ! ---     infiltration
                  else
                     ! Pvw_begin , implemented by KRO_20170907
                     !           Limit the infiltration head (-dh, dh<0. here) to the waterdepth in the channel
                     if (swdtyp(lev) .eq. 2) then
                        if (swliminf .eq. 1) then
                           dh = max(dh, (zbotdr(lev) - x))
                        end if
                     end if
                     ! Pvw_end   , implemented by KRO_20170907
                     qdrain(lev) = dh/infres(lev)
                     if (swallo(lev) .eq. 3 .or. fldry) qdrain(lev) = 0.0d0
                  end if
               end do

               ! --- drainage flux from table with gwlevel - flux data pairs
            elseif (dramet .eq. 1) then
               qdrain(1) = afgen(qdrtab, 50, abs(gwldra))
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

               ! SS-SWST Phase 2 Task 11 B1: removed globals flInitDraBas,ZDraBas,inqdra*,iqdra,
               ! cqdra,cqdrain*,qdrtot — now written only via state%surfacewater.
               ! SS-DRST Phase 2 Task 3: qdra dropped — all qdra reads/writes use state%drainage%qdra.
               ! SS-DRST Phase 2 Task 4: qdrain dropped — bocodrb writes state%drainage%qdrain directly.
               ! ADR 0031 Phase 2 Task 5: zTopDisLay removed from use-list; declared local below.
               ! SS-SWC Phase 2 S-2.8: gwl removed from use-list; read from state%soilwater%gwl.
               ! [SS-SWC S-2.12B] fluseksatexm retired from variables — read via state%soilwater
               ! [SS-TC TC-14] t1900, dt read via state%timecontrol (ADR 0041)
               use variables, only: nrlevs,numnod,dramet,swdtyp,NumLevRapDra,owltab,nowltab, &
                  zbotdr,swdivd,swdislay,swtopdislay,fTopDisLay, &
                  dz,ksatfit,ksatexm,layer,cofani,l,Swdivdinf,Swnrsrf,    &
                  SwTopnrsrf,FacDpthInf,madr
               use array_utils, only: afgen

               type(swap_state_t), intent(inout) :: state

               !     local
               integer node, level
               real(8) zCum, zTopDisLay(madr), difzTopDisLay(madr), ratio, ratiodz, sumqdr(madr), dh
               !, temptab(2*maowl) ????
               integer nodeTopDisLay(madr)
               CHARACTER(len=33) messag

               ! Allocate per-level state arrays if not yet done (guard for
               ! fldrain path where surfacewater_init may not have been called).
               if (.not. allocated(state%surfacewater%cqdrain)) then
                  allocate(state%surfacewater%cqdrain(nrlevs))
                  state%surfacewater%cqdrain = 0.0d0
               end if
               if (.not. allocated(state%surfacewater%cqdrainin)) then
                  allocate(state%surfacewater%cqdrainin(nrlevs))
                  state%surfacewater%cqdrainin = 0.0d0
               end if
               if (.not. allocated(state%surfacewater%cqdrainout)) then
                  allocate(state%surfacewater%cqdrainout(nrlevs))
                  state%surfacewater%cqdrainout = 0.0d0
               end if
               if (.not. allocated(state%surfacewater%inqdra)) then
                  allocate(state%surfacewater%inqdra(nrlevs, numnod))
                  state%surfacewater%inqdra = 0.0d0
               end if
               if (.not. allocated(state%surfacewater%inqdra_in)) then
                  allocate(state%surfacewater%inqdra_in(nrlevs, numnod))
                  state%surfacewater%inqdra_in = 0.0d0
               end if
               if (.not. allocated(state%surfacewater%inqdra_out)) then
                  allocate(state%surfacewater%inqdra_out(nrlevs, numnod))
                  state%surfacewater%inqdra_out = 0.0d0
               end if

               !   - In case of macropores: initialise drainage basis for rapid drainage through macropores
               if (state%surfacewater%flInitDraBas) then
                  if (NumLevRapDra .gt. nrlevs) then
                     messag = ' NUMLEVRAPDRA greater then NRLEVS'
                     call fatalerr_collected('MacroRead', messag)
                  end if

                  ! SS-SWST Phase 2 Task 11: ZDraBas global dropped; write only to state.
                  if (dramet .lt. 3) then
                     state%surfacewater%ZDraBas = zbotdr(1)
                  else
                     if (swdtyp(NumLevRapDra) .eq. 1) then
                        state%surfacewater%ZDraBas = zbotdr(NumLevRapDra)
                     else
                        !do i = 1,2*maowl
                        !   temptab(i) = owltab(NumLevRapDra,i)
                        !end do
                        !state%surfacewater%ZDraBas = afgen (temptab,2*maowl,t1900)
                        state%surfacewater%ZDraBas = afgen(owltab(NumLevRapDra, 1:2*nowltab(NumLevRapDra)), 2*nowltab(NumLevRapDra), state%timecontrol%t1900)
                     end if
                  end if

                  ! SS-SWST Phase 2 Task 11 B1: flInitDraBas global write dropped.
                  state%surfacewater%flInitDraBas = .false.

                  Return

               end if

               ! --- reset intermediate surface-water and drainage fluxes
               ! Identical reset at both call sites (here and SurfaceWater(2));
               ! delegate to reset_intermediate(). See surfacewater_state_mod.
               if (state%timecontrol%flZeroIntr) call state%surfacewater%reset_intermediate()

               ! --- reset cumulative drainage fluxes
               ! reset_cumulative_drainage zeros only the drainage-owned fields
               ! (cqdra, cqdrain*) — gated by fldrain, active under swdra=1 OR
               ! swdra=2. The reservoir-owned fields (cqdrd, cwsupp, cwout) are
               ! reset by SurfaceWater(2) and never accumulate under swdra=1, so
               ! no reservoir reset site is needed here. See ADR 0042.
               if (state%timecontrol%flZeroCumu) call state%surfacewater%reset_cumulative_drainage()

               ! --- reset to zero if groundwater level under soil profile and return
               ! SS-SWC Phase 2 S-2.8: gwl read from state%soilwater
               if (state%soilwater%gwl .gt. 998.0d0) then
                  ! SS-DRST Phase 2 Task 3: zero state directly; legacy global qdrain
                  ! still zeroed so bocodrb's own use-variables path stays consistent.
                  do level = 1, nrlevs
                     state%drainage%qdrain(level) = 0.0d0
                  end do
                  return
               end if

               ! --- calculate total drainage rate and state variables
               ! SS-DRST Phase 2 Task 4: bocodrb writes state%drainage%qdrain directly; no bridge sync.
               call bocodrb(dh, state)

               ! --- partition drainage flux over compartments
               ! SS-DRST Phase 2 Task 3: divdra reads/writes state%drainage%qdrain and
               ! state%drainage%qdra directly — no legacy globals passed here.
               if (swdivd .eq. 1) then
                  ! SS-SWC Phase 2 S-2.8: gwl read from state%soilwater
                  call divdra(numnod, nrlevs, dz, ksatfit, ksatexm, state%soilwater%fluseksatexm,    &  ! [SS-SWC S-2.12B]
              &      layer, cofani, state%soilwater%gwl, l, state%drainage%qdrain, state%drainage%qdra, &
              &      Swdivdinf, Swnrsrf, SwTopnrsrf, Zbotdr, state%timecontrol%dt, FacDpthInf, owltab, state%timecontrol%t1900)  !  Divdra, infiltration
                  !       redistribute qdrain with new top boundary for discharge layers
                  if (swdislay .eq. 2) then
                     do level = 1, nrlevs
                        if (swtopdislay(level) .eq. 1) then
                           zTopDisLay(level) = fTopDisLay(level)*state%soilwater%gwl +        &
              &                       (1.0d0 - fTopDisLay(level))*(state%soilwater%gwl - dh)
                        end if
                     end do
                  end if
                  if (swdislay .eq. 1 .or. swdislay .eq. 2) then
                     do level = 1, nrlevs
                        if (swtopdislay(level) .eq. 1) then
                           !                 find node nr of new top of discharge layer
                           nodeTopDisLay(level) = 1
                           zCum = -dz(1)
                           do while (zTopDisLay(level) .lt. zCum)
                              nodeTopDisLay(level) = nodeTopDisLay(level) + 1
                              zCum = zCum - dz(nodeTopDisLay(level))
                           end do
                           !                 saturated part (difzTopDisLay(lev)) of compartment containing waterlevel
                           difzTopDisLay(level) = zTopDisLay(level) - zCum
                           ratiodz =                                             &
              &                     difzTopDisLay(level)/dz(nodeTopDisLay(level))
                           sumqdr(level) =                                       &
              &                        ratiodz*state%drainage%qdra(level, nodeTopDisLay(level))
                           do node = nodeTopDisLay(level) + 1, numnod
                              sumqdr(level) = sumqdr(level) + state%drainage%qdra(level, node)
                           end do
                           if (dabs(sumqdr(level)) .lt. 1.0d-8) then
                              ratio = 1.0d0
                           else
                              ratio = state%drainage%qdrain(level)/sumqdr(level)
                           end if
                           !                 redistribute drainwater fluxes
                           do node = 1, nodeTopDisLay(level) - 1
                              state%drainage%qdra(level, node) = 0.0d0
                           end do
                           state%drainage%qdra(level, nodeTopDisLay(level)) =       &
              &                state%drainage%qdra(level, nodeTopDisLay(level))*ratio*ratiodz
                           do node = nodeTopDisLay(level) + 1, numnod
                              state%drainage%qdra(level, node) = state%drainage%qdra(level, node)*ratio
                           end do
                        end if
                     end do
                  end if
               else
                  ! --- drainage flux through lowest compartment
                  do level = 1, nrlevs
                     do node = 1, numnod - 1
                        state%drainage%qdra(level, node) = 0.0d0
                     end do
                     state%drainage%qdra(level, numnod) = state%drainage%qdrain(level)
                  end do
               end if

               ! SS-SWST Phase 2 Task 11 B1: qdrtot global write dropped; only state written.
               state%surfacewater%qdrtot = 0.0d0
               do level = 1, nrlevs
                  state%surfacewater%qdrtot = state%surfacewater%qdrtot + state%drainage%qdrain(level)
               end do

               ! SS-DRST Phase 2 Task 4: state%drainage%qdrain is authoritative; legacy global
               ! qdrain no longer written here or by bocodrb.

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
  use variables, only: swsec,swsrf,nrlevs,nrpri,zbotdr,taludr,widthr,pondmx,swdtyp,wlp,l,rdrain,rinfi,          &
              rentry, rexit, gwlinf, impend, nmper, wscap, swnrsrf, rsurfdeep, rsurfshallow, cofintfl,              &
                                    expintfl, FlMacropore, NumLevRapdra

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
               associate( &
                  wls    => state%surfacewater%wls,    &
                  swst   => state%surfacewater%swst,   &
                  ZDraBas => state%surfacewater%ZDraBas, &
                  qdrain  => state%drainage%qdrain,    &
                  gwl    => state%soilwater%gwl,       &
                  pond   => state%soilwater%pond)

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
                  if (swsrf .ge. 2) then
                     if (level .gt. nrpri) then
                        wl = wls
                     else
                        wl = wlp
                     end if
                  end if

! --- drainage fluxes are set to zero if both groundwater level and surface
!     water level are above ponding sill (so the nonzero drainage flux
!     is only computed if either the gwl or the wl is below pondmx)
                  if (wl .lt. pondmx .or. gwl .lt. pondmx) then

! --- channel is active medium if either groundwater or surface water
!     level is above channel bottom
                     if (gwl .gt. (zbotdr(level) + 0.001d0) .or.                       &
                &        wl .gt. (zbotdr(level) + 0.001d0)) then
                        if (wl .le. (zbotdr(level) + 0.001d0) .or. swsrf .eq. 1) then

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
                  if (swsrf .ge. 2 .and. level .gt. nrpri) then

! --- qdrd is total flux to or from secondary system
                     state%drainage%qdrd = state%drainage%qdrd + qdrain(level)
                  end if

500               continue

!   - drainage basis for rapid drainage through macropores
                  if (FlMacropore) then
                     if (swdtyp(NumLevRapDra) .ne. 1) then
                        ZDraBas = state%drainage%drainl(NumLevRapDra)
                     end if
                  end if

! ----------------------------------------------------------------------
! --- check for system falling dry (only for swsec = 2):
                  if (swsec .eq. 1) return
                  if (swsrf .eq. 1) then
                     do 10 level = 1, nrlevs
                        if (qdrain(level) .lt. 0.0d0) then
                           qdrain(level) = 0.0d0
                        end if
10                      continue
                        elseif (swsrf .ge. 2) then

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
