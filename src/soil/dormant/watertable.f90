!> @file src/soil/dormant/watertable.f90
!! @brief DORMANT — perched-water-table search with volume-integral heuristic.
!!
!! ## Status: DORMANT
!!
!! Extracted from `src/soil/waterbalance.f90` on 2026-05-24. `watertable`
!! has had **no callers** in the TOML-pipeline build since the migration:
!! the legacy dispatch site that invoked `watertable` (gated by
!! `CritUndSatVol > 0`) was not ported. The body is preserved verbatim
!! for future reactivation.
!!
!! `level()` is also preserved here as a verbatim copy so this dormant
!! module compiles self-contained when reactivated. The **live** copy of
!! `level()` remains in `src/soil/waterbalance.f90` because `calcgwl`
!! calls it. On reactivation, drop the dormant `level` and import it
!! from `soilwaterbalance_mod` instead — there is no reason to maintain
!! two copies.
!!
!! Excluded from `meson.build` — this file is NOT compiled.
!!
!! ## What it does
!!
!! `watertable` is an alternative to `calcgwl`'s simple "first node with
!! `h < 0`" criterion. It walks the profile accumulating total unsaturated
!! air volume; a saturated zone embedded in unsaturated soil but with
!! less than `CritUndSatVol` cm of accumulated air above it is still
!! treated as part of the main saturated body. This suppresses spurious
!! perched-water-table detections caused by numerical noise.
!!
!! `level` is a helper used by `calcgwl` and `watertable` to derive the
!! groundwater level from the pressure-head profile, with two
!! interpolation modes (h=0 crossing, or h=±1 averaging).
!!
!! ## Reactivation checklist
!!
!! 1. Add a config gate: e.g., `soil.critundsatvol > 0` or a dedicated
!!    `swwatertable` switch in `config%soil`.
!! 2. Migrate `CritUndSatVol` from `variables.f90` to `state%soilwater`
!!    or `state%cfg%soil` (Pattern 2 snapshot).
!! 3. Restore the dispatch site: legacy SWAP called `watertable` from
!!    inside `calcgwl` when the criterion was active. The wrapper in
!!    `calcgwl` reads gwl with `level(state,2,node,nodheq1)`; the
!!    `watertable` variant should branch in based on the gate.
!! 4. Add this file to `meson.build` sources.
!! 5. Re-export `level` and `watertable` from `waterbalance_mod` (or
!!    a new `watertable_mod` if you prefer) and update the call sites.
!! 6. Run `pixi run check-fast` to surface any retired-global imports
!!    (`use variables` may reference symbols that have moved to state
!!    in the meantime; this body's references will need updating).
!!
!! ## Original docstring
!!
!! Date: April 2008
!!
!! "Searches for the watertable and perched watertable (if existing)
!! using a criterion based on total unsaturated volume. An unsaturated
!! zone embedded in saturated soil must contain at least CritUndSatVol
!! cm of air to be recognized as truly unsaturated."
module watertable_dormant_mod
   use error_mod, only: fatalerr_collected
   use swap_state_mod, only: swap_state_t
   implicit none
   private
   public :: level, watertable

contains

   !> Calculate water level (elevation head) from pressure head.
   !!
   !! @param[in] swoptlev  1 = elevation head where h=0; 2 = average of
   !!                      elevation heads at h=-1 and h=+1
   !! @param[in] node      reference compartment index
   !! @param[in] nodheq1   node where h equals +1 (for swoptlev=2)
   function level (state, swoptlev, node, nodheq1)
      implicit none

      type(swap_state_t), intent(in) :: state
      integer, intent(in) :: swoptlev, node, nodheq1
      integer :: i
      real(8) :: levm1, levp1
      real(8) :: level

      associate (mesh => state%mesh, soil => state%soilwater)

      if (swoptlev .eq. 1) then
         ! groundwater level equals elevation head where h = 0
         if (soil%h(node+1) .ge. 0.0d0) then
            level = mesh%z(node+1) + soil%h(node+1) / (soil%h(node+1) - soil%h(node)) * mesh%disnod(node+1)
         else
            level = mesh%zbotcp(node) - soil%h(node)
            level = min(mesh%z(node), max(mesh%zbotcp(node), level))
         end if

      elseif (swoptlev .eq. 2) then
         ! groundwater level equals average of elevation heads of h = -1 and h = +1
         ! elevation head of h = +1
         i = nodheq1
         if (nodheq1 .eq. mesh%numnod) then
            levp1 = mesh%z(i) - 0.5d0 * mesh%dz(i)
         else
            levp1 = mesh%z(i) - (mesh%z(i) - mesh%z(i+1)) * (1.d0 - soil%h(i)) / (soil%h(i+1) - soil%h(i))
         end if
         ! elevation head of h = -1
         i = node
         do while (soil%h(i) .gt. -1.d0 .and. i .gt. 1)
            i = i - 1
         end do
         if (i .eq. 1 .and. soil%h(1) .gt. -1.d0 .and. node .gt. 2) then
            ! no compartment with pressure head < -1 cm in top of profile:
            ! use elevation head of h = 0 as estimation for groundwater level
            levm1 = mesh%z(node+1) + soil%h(node+1) / (soil%h(node+1) - soil%h(node)) * mesh%disnod(node+1)
            levp1 = levm1
         else
            levm1 = mesh%z(i+1) + (mesh%z(i) - mesh%z(i+1)) * (1.d0 + soil%h(i+1)) / (soil%h(i+1) - soil%h(i))
         end if
         ! groundwater level = average of levp1 and levm1
         level = (levp1 + levm1) / 2.d0
      end if

      end associate

      return
   end function level

   !> Search for watertable and perched watertable using volume-integral heuristic.
   subroutine watertable (state, node, nodlev, nodhlp, nodheq1, CritUndSatVol, flsat, waterlevel)
      implicit none

      type(swap_state_t), intent(in)    :: state
      integer,            intent(inout) :: node
      integer,            intent(out)   :: nodlev, nodhlp
      integer,            intent(in)    :: nodheq1
      real(8),            intent(in)    :: CritUndSatVol
      logical,            intent(inout) :: flsat
      real(8),            intent(out)   :: waterlevel
      integer :: i
      real(8) :: TotUndSatVol
      logical :: flsat2

      associate (mesh => state%mesh, soil => state%soilwater)

      TotUndSatVol = 0.0d0
      flsat2 = .false.
      i = node
      do while (TotUndSatVol .lt. CritUndSatVol .and. .not. flsat2 .and. i .ge. 1)
         TotUndSatVol = TotUndSatVol + (soil%thetas(i) - soil%theta(i)) * mesh%dz(i)
         if (soil%h(i) .gt. -1.d-7) flsat2 = .true.
         i = i - 1
      end do

      if (i .eq. 0 .or. TotUndSatVol .gt. CritUndSatVol - 1.d-8) then
         flsat  = .false.
         ! NB: nodlev is the DEEPEST UNSATURATED NODE, not the GWL node.
         nodlev = node
         nodhlp = i
      elseif (flsat2) then
         node = i + 1
      end if

      if (.not. flsat) then
         ! find groundwater level containing node
         if (CritUndSatVol .gt. 0.d0) then
            waterlevel = level(state, 1, node, nodheq1)
         else
            waterlevel = level(state, 2, node, nodheq1)
         end if
         i = max(node - 2, 1)
         do while (mesh%z(i) - 0.5d0*mesh%dz(i) .gt. waterlevel .and. i .gt. 2 .and. i .lt. mesh%numnod)
            i = i + 1
         end do
         nodlev = min(max(i, 1), mesh%numnod)
      end if

      end associate

      return
   end subroutine watertable

end module watertable_dormant_mod
