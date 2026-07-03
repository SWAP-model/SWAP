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
   public :: redistribute_qdra_over_discharge_layers

contains

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
      use array_utils, only: afgen

      real(8) dh
      type(swap_state_t), intent(inout) :: state

      integer :: lev
      real(8) :: zimp, dbot, pi, totres, x, eqd, rver, rhor, rrad
      real(8) :: gwldra
      logical :: fldry
      character(len=200) :: messag

      parameter(pi=3.14159d0)
      character(len=*), parameter :: imp_above_drain_msg = &
         'At the drainage section, the level of the'       &
         //' impervious layer is higher than the level of the' &
         //' drain bottom. Adapt drain input!'

      associate (drai => state%drainage,    &
                 soil => state%soilwater,   &
                 exch => state%exchange,     &
                 time => state%timecontrol)

         gwldra = exch%drain_soil%gwl

         ! --- drainage flux according to hooghoudt or ernst
         if (drai%dramet .eq. 2) then
            if (drai%shape .gt. small) dh = (gwldra - drai%zbotdr(1))/drai%shape

            ! --- contributing layer below drains limited to 1/4 l
            zimp = max(drai%basegw, drai%zbotdr(1) - 0.25*drai%l(1))
            dbot = (drai%zbotdr(1) - zimp)
            if (dbot .lt. 0.0d0) then
               messag = imp_above_drain_msg
               call state%diag%fatal('Bocodrb', messag)
               return
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

               ! --- equivalent depth (Hooghoudt)
               eqd = equivalent_depth_hooghoudt(dbot, drai%l(1), drai%wetper(1))

               if (drai%ipos .eq. 2) then
                  totres = drai%l(1)*drai%l(1)/(8*drai%khtop*eqd + 4*drai%khtop*abs(dh)) + drai%entres
               elseif (drai%ipos .eq. 3) then
                  totres = drai%l(1)*drai%l(1)/(8*drai%khbot*eqd + 4*drai%khtop*abs(dh)) + drai%entres
               end if
               drai%qdrain(1) = dh/totres

            ! --- case 4: drain in bottom layer
            elseif (drai%ipos .eq. 4) then
               if (drai%zbotdr(1) .gt. drai%zintf) then
                  messag = imp_above_drain_msg
                  call state%diag%fatal('bocodrb', messag)
                  return
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
                  messag = imp_above_drain_msg
                  call state%diag%fatal('bocodrb', messag)
                  return
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
                     drai%qdrain(lev) = dh/drai%drares(lev)
                     if (drai%swallo(lev) .eq. 2) drai%qdrain(lev) = 0.0d0
                  end if

               ! --- infiltration
               else
                  ! Limit infiltration head (-dh) to channel waterdepth (KRO 2017-09-07)
                  if (drai%swdtyp(lev) .eq. 2) then
                     if (drai%swliminf .eq. 1) dh = max(dh, (drai%zbotdr(lev) - x))
                  end if
                  drai%qdrain(lev) = dh/drai%infres(lev)
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
      real(8) :: dh
      character(len=33) :: messag

      associate (mesh => state%mesh,         &
                 drai => state%drainage,     &
                 soil => state%soilwater,    &
                 exch => state%exchange,     &
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
         ! [FIX-ADAPTIVEDT 2026-06-11] One-time ZDraBas setup. In legacy SWAP
         ! 4.2.0 this was a separate init-phase call (Drainage task=1) made
         ! BEFORE the time loop; the first per-step call (task=2) then computed
         ! drainage. The strangler refactor collapsed both into this per-step
         ! routine but made the init and the rate computation mutually exclusive
         ! (init in the `if`, rates in the `else`). For basic-drainage cases
         ! (swdra=1) flInitDraBas is never cleared during init (the heavy
         ! surfacewater init that clears it is gated on swsec==2), so the FIRST
         ! timestep took the init-only branch and computed ZERO drainage — a
         ! wasted "no-op" step. With adaptive dt that no-op solve returns
         ! numbit=1, the controller doubles dt prematurely, and the modern run
         ! desyncs from 4.2.0 for the rest of the simulation (root cause of the
         ! hysteresis/winter/swinter3 known-divergences and the swdrought2 perf
         ! collapse). Fix: do the one-time init, then FALL THROUGH to compute
         ! drainage on the same call — matching legacy's init-then-step order.
         ! At this first call time%t1900 == tstart (timecontrol advances after
         ! drainage), identical to legacy's init-time value, so ZDraBas is
         ! computed bit-identically; for basic drainage ZDraBas is dead anyway
         ! (only macropore rapid drainage reads it, retired by ADR 0040).
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

         end if

         block  ! normal timestep path — runs on EVERY call, including the first

            ! Reset intermediate surface-water and drainage fluxes.
            if (time%flZeroIntr) call surf%reset_intermediate()

            ! Reset cumulative drainage fluxes (cqdra/cqdrain*); reservoir-owned
            ! fields are reset by SurfaceWater(2) (ADR 0042).
            if (time%flZeroCumu) call surf%reset_cumulative_drainage()

            ! Skip when groundwater level is below the profile.
            if (exch%drain_soil%gwl .gt. 998.0d0) then
               do level = 1, drai%nrlevs
                  drai%qdrain(level) = 0.0d0
               end do
            else

               ! Total drainage rate per level.
               call bocodrb(dh, state)

               ! Partition the drainage flux over soil compartments.
               if (drai%swdivd .eq. 1) then
                  call divdra(mesh%numnod, drai%nrlevs, mesh%dz, soil%ksatfit, soil%ksatexm, &
                              soil%fluseksatexm, mesh%layer, soil%cofani, exch%drain_soil%gwl,          &
                              drai%L, drai%qdrain, drai%qdra,                                &
                              drai%swdivdinf, drai%swnrsrf, drai%swtopnrsrf, drai%zbotdr,    &
                              time%dt, drai%FacDpthInf, drai%owltab, drai%nowltab, time%t1900)

                  ! Redistribute qdrain with the new top boundary for discharge layers.
                  call redistribute_qdra_over_discharge_layers(state, dh)
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

         end block  ! normal timestep path

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
      real(8) :: dh
      type(swap_state_t), intent(inout) :: state

      integer :: level

      associate (drai => state%drainage,    &
                 soil => state%soilwater,   &
                 exch => state%exchange,     &
                 surf => state%surfacewater)

         ! Spec D7: zero drainage when groundwater is dry.
         if (exch%drain_soil%gwl .gt. 998.0d0) then
            do level = 1, drai%nrlevs
               drai%qdrain(level) = 0.0d0
            end do
            return
         end if

         drai%qdrd = 0.0d0

         do level = 1, drai%nrlevs
            call compute_drainage_level_flux(state, level, dh)
            ! qdrd: total flux to/from secondary system.
            if (surf%swsrf .ge. 2 .and. level .gt. surf%nrpri) then
               drai%qdrd = drai%qdrd + drai%qdrain(level)
            end if
         end do

         ! [GR-SOIL …] Macropore rapid-drainage basis block deleted (ADR 0040).

         ! Check for system falling dry (only for swsec = 2):
         if (surf%swsec .eq. 1) return
         if (surf%swsrf .eq. 1) then
            do level = 1, drai%nrlevs
               if (drai%qdrain(level) .lt. 0.0d0) drai%qdrain(level) = 0.0d0
            end do
         elseif (surf%swsrf .ge. 2) then
            call apply_storage_dry_check(state)
         end if

      end associate
    end subroutine bocodre

   !> Compute drainage/infiltration flux for a single drainage level.
   !!
   !! Sets drai%qdrain(level), drai%drainl(level), optionally drai%wetper(level),
   !! and updates dh. dh is intent(inout): in the outermost else-branch (both wl
   !! and gwl above pondmx) neither dh nor qdrain is assigned, preserving the
   !! carry-forward value from the previous iteration.
   subroutine compute_drainage_level_flux(state, level, dh)
      use swap_state_mod, only: swap_state_t
      implicit none
      type(swap_state_t), intent(inout) :: state
      integer,            intent(in)    :: level
      real(8),            intent(inout) :: dh

      real(8) :: wl, swdepth, swexbrd, rd, re

      associate (drai => state%drainage,    &
                 soil => state%soilwater,   &
                 exch => state%exchange,     &
                 surf => state%surfacewater)

         ! Surface-water level for this drainage level.
         if (surf%swsrf .ge. 2) then
            if (level .gt. surf%nrpri) then
               wl = surf%wls
            else
               wl = surf%wlp
            end if
         end if

         ! Drainage fluxes set to zero if both gwl and surface water level are
         ! above the ponding sill (compute non-zero flux only when at least one
         ! is below pondmx).
         if (wl .lt. surf%pondmx .or. exch%drain_soil%gwl .lt. surf%pondmx) then

            ! Channel is active if either gwl or surface water is above bottom.
            if (exch%drain_soil%gwl .gt. (drai%zbotdr(level) + 0.001d0) .or.        &
        &       wl       .gt. (drai%zbotdr(level) + 0.001d0)) then
               if (wl .le. (drai%zbotdr(level) + 0.001d0) .or. surf%swsrf .eq. 1) then

                  ! Only groundwater above channel bottom: bottom is drainage base.
                  drai%drainl(level) = drai%zbotdr(level)
                  if (drai%swdtyp(level) .eq. 0) then
                     drai%wetper(level) = drai%widthr(level)
                  end if
               else
                  ! Surface-water level above channel bottom.
                  drai%drainl(level) = wl
                  if (drai%swdtyp(level) .eq. 0) then
                     swdepth = wl - drai%zbotdr(level)
                     swexbrd = (wl - drai%zbotdr(level))/drai%taludr(level)
                     drai%wetper(level) = drai%widthr(level) +           &
                                          2*dsqrt(swdepth**2 + swexbrd**2)
                  end if
               end if

               ! Drainage flux (cm/d): head difference.
               dh = exch%drain_soil%gwl - drai%drainl(level)
               if (exch%drain_soil%gwl .gt. -0.1d0) dh = dh + exch%drain_soil%pond
               if (dh .lt. 0.0d0 .and. exch%drain_soil%gwl .lt. drai%gwlinf(level)) then
                  dh = drai%gwlinf(level) - drai%drainl(level)
               end if
               ! Interflow as power function.
               if ((level .eq. drai%nrlevs) .and. (drai%swnrsrf .eq. 2)) then
                  drai%qdrain(level) = drai%cofintfl*dh**drai%expintfl
               else
                  if (dh .gt. 0.0d0) then
                     rd = drai%rdrain(level)
                     re = drai%rentry(level)
                     ! Surface drainage (vacuum cleaner).
                     if ((level .eq. drai%nrlevs) .and. (drai%swnrsrf .eq. 1)) then
                        rd = drai%rsurfdeep - dh
                        rd = max(rd, drai%rsurfshallow)
                     end if
                  else
                     rd = drai%rinfi(level)
                     re = drai%rexit(level)
                  end if
                  if (drai%swdtyp(level) .eq. 0) then
                     drai%qdrain(level) = dh/((rd + re*drai%L(level)/drai%wetper(level)))
                  else
                     drai%qdrain(level) = dh/rd
                  end if
               end if
            else
               if (drai%swdtyp(level) .eq. 0) drai%wetper(level) = 0.0d0
               dh = 0.0d0
               drai%qdrain(level) = 0.0d0
               ! Drainage basis for rapid drainage through macropores.
               drai%drainl(level) = drai%zbotdr(level)
            end if
         else
            drai%qdrain(level) = 0.0d0
            ! Drainage basis for rapid drainage through macropores.
            drai%drainl(level) = wl
         end if

      end associate
   end subroutine compute_drainage_level_flux

   !> Check if the secondary surface-water system falls dry and rescale fluxes.
   !!
   !! For swsec=2 + swsrf>=2: finds the current management period, computes
   !! potential storage change, and if storage would drop below zero rescales
   !! all secondary-level drainage fluxes proportionally so that total withdrawal
   !! exactly matches available storage (swst + wscap).
   subroutine apply_storage_dry_check(state)
      use swap_state_mod, only: swap_state_t
      implicit none
      type(swap_state_t), intent(inout) :: state

      integer :: imper, level
      real(8) :: qdrdm, qdratio, dvmax, swstmax
      character(len=200) :: messag

      associate (drai => state%drainage,    &
                 surf => state%surfacewater, &
                 time => state%timecontrol)

         ! Determine which management period the model is in.
         imper = 0
         do
            imper = imper + 1
            if (imper .gt. surf%nmper) then
               messag = 'sw-management periods(IMPER), more than defined'
               call state%diag%fatal('Bocodre', messag)
               return
            end if
            if (time%t1900 - 1.d0 + 0.1d-10 .le. surf%impend(imper)) exit
         end do

         ! Will the system become empty?
         dvmax   = (drai%qdrd + surf%wscap(imper))*time%dt
         swstmax = surf%swst + dvmax

         if (swstmax .lt. 0.0d0) then
            ! Storage would drop below zero: rescale infiltration to
            ! exactly match available (swst + wscap).
            qdrdm   = -(surf%swst + surf%wscap(imper)*time%dt)/time%dt
            qdratio = qdrdm/drai%qdrd

            if (qdratio .gt. 1.0d0 .or. qdratio .lt. 0.0d0) then
               messag = 'sw-management error with storage (qdratio)'
               call state%diag%fatal('Bocodre', messag)
               return
            end if

            do level = 1 + surf%nrpri, drai%nrlevs
               drai%qdrain(level) = drai%qdrain(level)*qdratio
            end do
            drai%qdrd = qdrdm
         end if

      end associate
   end subroutine apply_storage_dry_check

   !> Redistribute drai%qdra over compartments for moving discharge-layer tops.
   !! Mutates drai%qdra in place based on swdislay (1=user-specified top,
   !! 2=dynamic from gwl and dh). No-op when swtopdislay(level) /= 1.
   subroutine redistribute_qdra_over_discharge_layers(state, dh)
      use swap_array_dimensions, only: madr
      use swap_state_mod,        only: swap_state_t
      implicit none
      type(swap_state_t), intent(inout) :: state
      real(8),            intent(in)    :: dh

      integer :: level, node, nodeTopDisLay(madr)
      real(8) :: zCum, zTopDisLay(madr), difzTopDisLay(madr), ratio, ratiodz, sumqdr(madr)

      associate (mesh => state%mesh, drai => state%drainage, soil => state%soilwater, exch => state%exchange)
         if (drai%swdislay .eq. 2) then
            do level = 1, drai%nrlevs
               if (drai%swtopdislay(level) .eq. 1) then
                  zTopDisLay(level) = drai%fTopDisLay(level)*exch%drain_soil%gwl +     &
                                      (1.0d0 - drai%fTopDisLay(level))*(exch%drain_soil%gwl - dh)
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
      end associate
   end subroutine redistribute_qdra_over_discharge_layers

   !> Hooghoudt's equivalent depth for drainage radial-flow resistance.
   !!
   !! Given the thickness `dbot` of the saturated zone below the drains,
   !! drain spacing `L`, and wetted perimeter `wetper`, returns the
   !! equivalent depth `eqd` (clamped to `dbot`). Uses a series expansion
   !! for x = 2*pi*dbot/L > 0.5 and a logarithmic approximation otherwise;
   !! degenerate small-x case returns dbot directly.
   pure function equivalent_depth_hooghoudt(dbot, L, wetper) result(eqd)
      implicit none
      real(8), intent(in) :: dbot, L, wetper
      real(8)             :: eqd
      real(8) :: x, fx
      integer :: i
      real(8), parameter :: pi = 3.14159d0   ! matches bocodrb's pi for byte-identity
      x = 2*pi*dbot/L
      if (x .gt. 0.5d0) then
         fx = 0.0d0
         do i = 1, 5, 2
            fx = fx + (4*exp(-2*i*x))/(i*(1.0d0 - exp(-2*i*x)))
         end do
         eqd = pi*L/8/(log(L/wetper) + fx)
      else
         if (x .lt. 1.0d-6) then
            eqd = dbot
         else
            fx = pi**2/(4*x) + log(x/(2*pi))
            eqd = pi*L/8/(log(L/wetper) + fx)
         end if
      end if
      if (eqd .gt. dbot) eqd = dbot
   end function equivalent_depth_hooghoudt

end module drainage_mod
