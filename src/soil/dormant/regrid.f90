!> @file src/soil/dormant/regrid.f90
!! @brief DORMANT — mid-simulation vertical re-discretization (`ConvertDiscrVert`).
!!
!! ## Status: DORMANT
!!
!! Originally lived in `src/soil/soilgrid.f90`; deleted in commit
!! `f9d6dc6` ("refactor(soil): delete dead ConvertDiscrVert (no
!! callers)") because the `SwDiscrvert==1` dispatch site had been
!! lost in the TOML migration. Restored here as a dormant module on
!! 2026-05-24 — the regrid capability is on the long-term roadmap
!! (output coupling to coarser solute-transport grids).
!!
!! Excluded from `meson.build` — this file is NOT compiled.
!!
!! ## What it does
!!
!! Projects the soil-water state (h, theta, theta_beg, fluxes, soil
!! temperature) from the simulation compute grid onto a different
!! grid for output. Two modes:
!!
!! - `SwDiscrvert == 0`: simple node-for-node copy (the output grid
!!   has the same compartments as the compute grid).
!! - `SwDiscrvert == 1`: weighted projection — compute-grid
!!   compartments are aggregated into coarser output compartments
!!   using volume-weighted averaging (theta, h via the inverse VG
!!   curve) and integration (root extraction, inter-comp flux,
!!   drainage flux).
!!
!! Used internally by output coupling (originally for ANIMO/PEARL
!! which needed a fixed coarse grid independent of SWAP's compute
!! discretization).
!!
!! ## Differences from the pre-deletion code
!!
!! The original supported a `swop == 2` (macropore) branch that
!! aggregated macropore-specific quantities (`VlMpStDm1`,
!! `DiPoCp`, `IQExcMtxDm1Cp`, etc.). **ADR 0040** permanently
!! retired macropores from SWAP. The dormant copy below has:
!!
!! - Dropped all `if (swop == 2)` branches and macropore-output
!!   array arguments from the signature.
!! - Reduced the API to the core state arrays.
!!
!! The unused `swop` argument is retained for now in case
!! reactivation needs to gate on a soil-water option, but if
!! macropores are not reinstated it can be dropped.
!!
!! ## Reactivation checklist
!!
!! 1. Restore the dispatch site: legacy SWAP called
!!    `ConvertDiscrVert(1, swop, ...)` once per simulation at
!!    init (to build the projection mapping) and
!!    `ConvertDiscrVert(2, swop, ...)` once per output period.
!!    Likely belongs in the output-coupling section of
!!    `swap_mod.f90`.
!!
!! 2. Migrate the three globals it still reads
!!    (`SwDiscrvert`, `numnodNew`, `dzNew`, `numlay`, `botcom`) to
!!    state:
!!    - `SwDiscrvert` → `state%cfg%output` or `state%cfg%mesh`
!!      (one-time switch).
!!    - `numnodNew`, `dzNew` → `state%mesh` as
!!      regrid-output-temporary fields (allocated to `macp`).
!!    - `numlay`, `botcom` → already in `state%mesh` (the global
!!      readers in this file should be migrated to `mesh%numlay`,
!!      `mesh%botcom`).
!!
!! 3. Add this file to `meson.build` sources.
!!
!! 4. Add `use regrid_dormant_mod, only: ConvertDiscrVert` to the
!!    dispatch site.
!!
!! 5. Run `pixi run check-fast` — broken imports surface
!!    immediately. `prhead` signature in particular may have
!!    drifted; check `soilhydraulics_utils`.
module regrid_dormant_mod
   use error_mod,             only: fatalerr_collected
   use swap_state_mod,        only: swap_state_t
   use swap_array_dimensions, only: macp, maho, madr
   use hydraulic_params_mod,  only: vanGenuchten_params_t
   use soilhydraulics_utils,  only: prhead
   use variables, only: SwDiscrvert, numlay, botcom, numnodNew, dzNew
   implicit none
   private
   public :: ConvertDiscrVert

contains

   !> Project soil-water state onto a coarser output grid.
   !!
   !! @param[in]  part         1 = init (build NodeNew mapping); 2 = dynamic update
   !! @param[in]  swop         legacy soilwater option (macropore branches removed)
   !! @param[out] botcomNew    bottom-compartment-per-layer on the new grid
   !! @param[out] hNew         pressure head on the new grid
   !! @param[out] thetaNew     water content on the new grid
   !! @param[out] inqNew       inter-compartment flux on the new grid
   !! @param[out] inqrotNew    root extraction per compartment on the new grid
   !! @param[out] inqdraNew    drainage flux per level per compartment on the new grid
   !! @param[out] iThetaBegNew initial water content on the new grid
   !! @param[in]  Tsoil        soil temperature on the compute grid
   !! @param[out] TsoilNew     soil temperature on the new grid
   !! @param[in]  state        swap state (mesh, soilwater, drainage, surfacewater)
   subroutine ConvertDiscrVert(part, swop, botcomNew, hNew, thetaNew, inqNew, &
                               inqrotNew, inqdraNew, iThetaBegNew, Tsoil, TsoilNew, state)
      implicit none

      integer,            intent(in)    :: part, swop
      integer,            intent(out)   :: botcomNew(maho)
      real(8),            intent(out)   :: hNew(macp), thetaNew(macp), inqrotNew(macp)
      real(8),            intent(out)   :: inqNew(macp+1)
      real(8),            intent(out)   :: inqdraNew(madr, macp)
      real(8),            intent(in)    :: Tsoil(0:macp)
      real(8),            intent(out)   :: TsoilNew(0:macp)
      real(8),            intent(out)   :: iThetaBegNew(macp)
      type(swap_state_t), intent(in)    :: state

      integer :: lay, node, nodeN, nodeNew(macp, 2), i, level
      real(8) :: disnodNew(macp+1), total, zNew(macp)
      type(vanGenuchten_params_t) :: vg_tentative
      character(len=80) :: message
      character(len=*), parameter :: ModuleName = 'ConvertDiscrVert'

      ! `swop` retained for API parity with legacy callers but unused
      ! now that macropore branches are removed (ADR 0040). Suppress
      ! the unused-arg warning by touching it:
      if (swop < 0) continue

      associate (mesh => state%mesh,         &
                 soil => state%soilwater,    &
                 drai => state%drainage,     &
                 surf => state%surfacewater)

      if (part .lt. 1 .or. part .gt. 2) then
         write(message, *) 'fatal error in variable PART'
         call fatalerr_collected(ModuleName, message)
      end if

      if (SwDiscrVert .eq. 0) then
         ! Simple node-for-node copy
         numnodNew = mesh%numnod
         do lay = 1, numlay
            botcomNew(lay) = botcom(lay)
         end do
         do node = 1, numnodNew
            dzNew(node)        = mesh%dz(node)
            hNew(node)         = soil%h(node)
            thetaNew(node)     = soil%theta(node)
            IThetaBegNew(node) = soil%IThetaBeg(node)
            inqNew(node)       = soil%inq(node)
            inqrotNew(node)    = soil%inqrot(node)
            do level = 1, drai%nrlevs
               inqdraNew(level, node) = surf%inqdra(level, node)
            end do
         end do
         do node = 0, numnodNew
            TsoilNew(node) = Tsoil(node)
         end do
         inqNew(numnodNew+1) = soil%inq(mesh%numnod+1)

      else if (SwDiscrVert .eq. 1) then

         if (part .eq. 1) then
            ! zNew = depth at bottom of new compartment
            total = 0.0d0
            do node = 1, numnodNew
               total = total + dzNew(node)
               zNew(node) = total
            end do

            ! NodeNew(:,1) = top old comp, NodeNew(:,2) = bottom old comp in each new comp
            nodeN = 1
            node = 1
            total = 0.0d0
            do while (node .le. mesh%numnod)
               total = total + mesh%dz(node)
               if (abs(zNew(nodeN) - total) .lt. 1.0d-6) then
                  NodeNew(nodeN, 2) = node
                  nodeN = nodeN + 1
               end if
               node = node + 1
            end do
            NodeNew(1, 1) = 1
            do node = 2, numnodNew
               NodeNew(node, 1) = NodeNew(node-1, 2) + 1
            end do

            ! botcomNew
            do lay = 1, numlay
               do node = 1, numnodNew
                  if (NodeNew(node, 2) .eq. botcom(lay)) then
                     botcomNew(lay) = node
                  end if
               end do
            end do

            ! For hNew: position and distances between nodal points on new grid
            zNew(1)      = -0.5d0 * dzNew(1)
            disnodNew(1) = -zNew(1)
            do node = 2, numnodNew
               zNew(node)      = zNew(node-1) - 0.5d0*(dzNew(node-1) + dzNew(node))
               disnodNew(node) = zNew(node-1) - zNew(node)
            end do
         end if

         ! In both part 1 and 2: theta on new grid via volume-weighted average
         do node = 1, numnodNew
            total = 0.0d0
            do i = NodeNew(node, 1), NodeNew(node, 2)
               total = total + mesh%dz(i)
            end do
            thetaNew(node)     = 0.0d0
            IThetaBegNew(node) = 0.0d0
            do i = NodeNew(node, 1), NodeNew(node, 2)
               thetaNew(node)     = thetaNew(node)     + soil%theta(i)      * mesh%dz(i) / total
               IThetaBegNew(node) = IThetaBegNew(node) + soil%IThetaBeg(i)  * mesh%dz(i) / total
            end do
         end do

         ! Tsoil on new grid via volume-weighted average
         do node = 1, numnodNew
            total = 0.0d0
            do i = NodeNew(node, 1), NodeNew(node, 2)
               total = total + mesh%dz(i)
            end do
            TsoilNew(node) = 0.0d0
            do i = NodeNew(node, 1), NodeNew(node, 2)
               TsoilNew(node) = TsoilNew(node) + Tsoil(i) * mesh%dz(i) / total
            end do
         end do

         if (part .eq. 2) then
            ! hNew via inverse VG curve using tentative params from bottom old node
            do node = 1, numnodNew
               vg_tentative = soil%vg_params(NodeNew(node, 2))
               hNew(node) = prhead(disnodNew(node), thetaNew(node), hNew,    &
                                   soil%iHWCKmodel(soil%layer(node)),       &
                                   node, soil, vg_in=vg_tentative)
            end do

            ! inqrot on new grid via integration
            do node = 1, numnodNew
               inqrotNew(node) = 0.0d0
               do i = NodeNew(node, 1), NodeNew(node, 2)
                  inqrotNew(node) = inqrotNew(node) + soil%inqrot(i)
               end do
            end do

            ! inq on new grid
            do node = 1, numnodNew
               inqNew(node) = soil%inq(NodeNew(node, 1))
            end do
            inqNew(numnodNew+1) = soil%inq(mesh%numnod+1)

            ! inqdra on new grid via integration
            do level = 1, drai%nrlevs
               do node = 1, numnodNew
                  inqdraNew(level, node) = 0.0d0
                  do i = NodeNew(node, 1), NodeNew(node, 2)
                     inqdraNew(level, node) = inqdraNew(level, node) + surf%inqdra(level, i)
                  end do
               end do
            end do
         end if

         ! Macropore branches (swop == 2) — dropped (ADR 0040; macropores retired).

      else
         write(message, *) 'fatal error in variable SwDiscrVert'
         call fatalerr_collected(ModuleName, message)
      end if

      end associate

      return
   end subroutine ConvertDiscrVert

end module regrid_dormant_mod
