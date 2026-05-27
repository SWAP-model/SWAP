!> @file src/crop/dormant/jongvanlier.f90
!! @brief DORMANT — De Jong van Lier (2013) microscopic root-uptake path (`swdrought=2`).
!!
!! ## Status: DORMANT
!!
!! Moved here from `src/crop/rootextraction.f90` on 2026-05-25. The
!! De Jong van Lier microscopic-uptake compute path (`swdrought = 2`)
!! has no live dispatch in the TOML pipeline — all 6 regression cases
!! run with `swdrought = 1` (Feddes linear reduction). The dispatch
!! site at the top of `RootExtraction` (`if (swdrought .eq. 2) call
!! JongvanLier(state)`) has been replaced with a
!! `fatalerr_collected` stub pointing at this file.
!!
!! Excluded from `meson.build` and `tests/unit/meson.build` —
!! this file is NOT compiled.
!!
!! ## What it contains
!!
!! - `subroutine JongvanLier(state)` — the outer Newton-Raphson loop
!!   that searches for the leaf pressure head (`hleaf`) at which
!!   `qrosum = ptra`. Two convergence regimes: `flstress=false`
!!   (no drought stress; iterates on `hleaf`) and `flstress=true`
!!   (drought stress; iterates on `Tactual` with `hleaf` pinned at
!!   `wiltpoint`).
!! - `subroutine JongvanLierLoop(state)` — inner per-node loop that
!!   solves the soil-root pressure-head balance and computes
!!   `qrot(node)` from matric-flux-potential differences. Reads
!!   `criterhr`, `stephr`, `kroot`, `rxylem`, `rootcoefa`, `rooteff`,
!!   `flhydrlift`, `twilt`.
!!
!! Both bodies share live dependencies on:
!! - `state%soilwater%{h, hroot, hleaf, mflux, mroot, qrot, qrosum,
!!   rmax, rootrho, rootphi, alpJvLier, Hxylem, Tactual, ksatfit}`
!! - `state%crop%common%{rd, rdctb, noddrz, swsalinity, salthead}`
!! - `state%mesh%{layer, z, ztopcp, dz, disnod, numnod}`
!! - `state%atmosphere%ptra`
!! - `state%solute%cml(:)` (osmotic-head correction)
!! - `state%timecontrol%dt`
!! - `state%cfg%simulation%numerical%taccur`
!! - `MatricFlux(2, ...)` from `rootextraction_mod` (LIVE — still
!!   compiled because cropgrowth.f90 / cropfixed/grass/wofost runtime
!!   call `MatricFlux(1, ...)` to build the lookup table).
!!
!! ## Reactivation checklist
!!
!! 1. Restore the dispatch site in `RootExtraction(state)`:
!!    replace the `fatalerr_collected('RootExtraction', 'swdrought=2
!!    (de Jong van Lier) path is dormant ...')` with
!!    `call JongvanLier(state)`.
!! 2. Restore the alpdry assignment for `swdrought == 2`:
!!    `alpdry = soil%alpJvLier` inside the per-node combination loop
!!    (search for the same `fatalerr` stub in the drought reduction
!!    block).
!! 3. Wire `criterhr`, `stephr`, `kroot`, `kstem`, `rxylem`,
!!    `rootradius`, `rootcoefa`, `rooteff` through a new
!!    `crop_jvl_t` sub-record on `state%cfg%crop` (currently no TOML
!!    reader — `cropfixed_config_t` has slots but no `[crop.fixed]`
!!    section is read for them). Update `config_to_variables` or the
!!    crop init to write them into a runtime state field (e.g.
!!    `state%crop%common%jvl%{...}`).
!! 4. Drop the `use variables, only: ...` block in this file and
!!    cut the reads over to the new state fields.
!! 5. Add this file to `meson.build` AND `tests/unit/meson.build`
!!    source lists (under the Crop section).
!! 6. Add `use jongvanlier_dormant_mod, only: JongvanLier` to the
!!    head of `rootextraction.f90` (rename the module's `_dormant`
!!    suffix or keep as-is; the live caller will need the same `use`
!!    line).
!! 7. The retired-legacy declarations (`rootcoefa`, `rooteff`,
!!    `rootradius`, `kstem`) and orphan stubs (`CriterHr`, `Kroot`,
!!    `Rxylem`, `StepHr`) have been DROPPED in this commit — restore
!!    them or move directly to `state%cfg%crop` (preferred — they
!!    are scalar parameters, no per-node array).
!!
!! Original SVN revision:
!!   $Id: rootextraction.f90 374 2018-03-21 13:12:23Z heine003 $

module jongvanlier_dormant_mod
   use error_mod, only: fatalerr_collected
   use swap_state_mod, only: swap_state_t
   use rootextraction_mod, only: matric_flux
   implicit none
   private
   public :: JongvanLier, JongvanLierLoop

   contains

! ----------------------------------------------------------------------
      !> Calculate microscopic root extraction after De Jong van Lier et al. (2013).
      !!
      !! @note
      !! date: August 2016
      !! purpose: Calculate the root water extraction rate according to
      !!          De Jong van Lier et al. (2013)
      !! @endnote
        subroutine JongvanLier(state)
! ----------------------------------------------------------------------
!     date      : August 2016
!     purpose   : Calculate the root water extraction rate according to
!                 De Jong van Lier et al. (2013)
! ----------------------------------------------------------------------
      use swap_array_dimensions, only: macp
      use swap_log, only: log_warn
      ! [GR-CROP 2026-05-25] retained — dead-branch (swdrought=2 stub-errored in TOML)
      use variables, only: kroot, kstem, rootcoefa,        &
                           rootradius, rxylem, stephr,     &
                           wiltpoint
      use array_utils, only: afgen
      implicit none

      type(swap_state_t), intent(inout) :: state

! --- local variables
      integer node,counter
      real(8) Fy3
      real(8) y1,y2,y3,Fy1,Fy2,hwet,ratio
      real(8) hleafm1,hrootm1(macp)
      real(8) reldepth,rdensity,phi,meandepth
      real(8) maxstep,dHPLant,dHPlantMax
      logical flconverg,flstress
      character(len=200) messag

! --- [GR-CROP 2026-05-25] sub-record aliases (crop, soil, atmo, mesh, cfg_sim).
      associate( &
         crop    => state%crop,                       &
         soil    => state%soilwater,                  &
         atmo    => state%atmosphere,                 &
         mesh    => state%mesh,                       &
         cfg_sim => state%cfg%simulation,             &
         noddrz  => state%crop%common%noddrz          &  ! [GR-CROP 2026-05-25] alias for state read
      )

! --- initialization
      phi = 3.1415926d0
      counter = 0
      flstress = .false.
      hleafm1 = soil%hleaf
      flconverg = .false.

! --- reset values below root zone to zero
      do node = noddrz+1,mesh%numnod
        soil%mflux(node)   = 0.d0
        soil%mroot(node)   = 0.d0
        soil%hroot(node)   = 0.d0
        soil%rootrho(node) = 0.d0
        soil%rootphi(node) = 0.d0
      enddo

! --- give hroot a value when root zone becomes larger and store previous values
      do node = 1, noddrz
        if (abs(soil%hroot(node)) .lt. 1.0d-12) then
          soil%hroot(node) = soil%h(node)
        end if
        hrootm1(node) = soil%hroot(node)
      enddo

! --- initialization of rootrho and rootphi
      do node = 1,noddrz-1
        reldepth = -mesh%z(node)/crop%common%rd
        rdensity = afgen(crop%common%rdctb,22,reldepth)
        soil%rmax(node) = 1.d0/dsqrt(phi*rdensity)
        soil%rootrho(node) = 4.d0/(rootradius*rootradius-rootcoefa*soil%rmax(node)&
     &      *rootcoefa*soil%rmax(node) + 2.d0 * (soil%rmax(node)*soil%rmax(node) +     &
     &      rootradius*rootradius)*log(rootcoefa*soil%rmax(node)/rootradius))
        soil%rootphi(node) = soil%rootrho(node) * soil%rmax(node)*soil%rmax(node) *         &
     &     log(rootradius/rxylem) * 0.5d0 / kroot
      enddo
! --- last node, partly filled with roots
      node = noddrz
      meandepth = (mesh%ztopcp(node)-crop%common%rd)*0.5d0
      reldepth = meandepth/(-crop%common%rd)
      rdensity = afgen(crop%common%rdctb,22,reldepth)
      soil%rmax(node) = 1.d0/dsqrt(phi*rdensity)
      soil%rootrho(node) = 4.d0/(rootradius*rootradius-rootcoefa*soil%rmax(node)  &
     &      *rootcoefa*soil%rmax(node) + 2.d0 * (soil%rmax(node)*soil%rmax(node) +     &
     &      rootradius*rootradius)*log(rootcoefa*soil%rmax(node)/rootradius))
      soil%rootphi(node) = soil%rootrho(node) * soil%rmax(node)*soil%rmax(node) *           &
     &     log(rootradius/rxylem) * 0.5d0 / kroot

! --- calculate current matric flux potential in soil water
      do node = 1,noddrz
        soil%h(node) = min(1.0d3,max(soil%h(node),1.0d-8))
        call matric_flux(soil%h(node),node,soil%mflux(node),state)
      enddo

! --- interpolate matric flux potential in lowest compartment that is partly filled with roots
      node = noddrz
      if (node.gt.1) then
         soil%mflux(node) = soil%mflux(node) + (soil%mflux(node-1)-soil%mflux(node))/       &
     &                 mesh%disnod(node) * (meandepth-mesh%z(node))
      endif

! --- determine highest pressure head in root zone
      hwet = soil%h(1)
      do node = 2,noddrz
        if (soil%h(node) .gt. hwet) hwet = soil%h(node)
      enddo

! --- determine whether qrosum < ptra (flstress = true)
! --- calculate root water extraction qrosum at Hleaf = wiltpoint and Tactual = Ptra
      soil%hleaf = wiltpoint
      dHPlant = atmo%ptra/Kstem
      if (dHPlant .gt. (hwet - wiltpoint)) flstress = .true.

      soil%Hxylem = soil%hleaf + dHPLant
      call JongvanLierLoop(state)
      if (soil%qrosum .lt. atmo%ptra) flstress = .true.

! --- set hroot to previous values
      do node = 1,noddrz
        soil%hroot(node) = hrootm1(node)
      enddo

! --- FLSTRESS = FALSE: DETERMINE HLEAF WITH TACTUAL = PTRA
      if (.not. flstress) then

! --- calculate root extraction at previous value of hleaf
      soil%hleaf = Hleafm1
      soil%Hxylem = soil%hleaf + dHPLant
      call JongvanLierLoop(state)

! --- check convergence
      if (abs(soil%qrosum-atmo%ptra) .lt. cfg_sim%numerical%taccur) flconverg = .true.

! --- set initial values y3 and Fy3
      y3 = soil%hleaf
      Fy3 = atmo%ptra - soil%qrosum

! --- start search on correct value of pressure head in leaves
      do while (.not. flconverg)

! --- reset y1 and Fy1
      y1 = y3
      Fy1 = Fy3

      if (soil%qrosum .lt. atmo%ptra) then
! ---   root water extraction too small, decrease Hleaf
          y2 = y1 - StepHr*log10(max(abs(y1),1.d0))
      else
! ---   root water extraction too large, increase Hleaf
          y2 = y1 + StepHr*log10(max(abs(y1),1.d0))
      endif

! --- calculate root extraction at modified pressure head y2 in leaf
      soil%hleaf = y2
      soil%Hxylem = soil%hleaf + dHPLant
      call JongvanLierLoop(state)
      Fy2 = atmo%ptra - soil%qrosum

! --- determine new value of leaf pressure head with Newton Raphson algorithm
      if (abs(Fy1-Fy2).lt. 1.d-10) then
! ---   take maximum step
        if (soil%qrosum .gt. atmo%ptra) then
          maxstep = -0.05d0 * y1
          maxstep = max(100.d0,maxstep)
          y3 = y1 + maxstep
        else
          maxstep = 0.05d0 * (y1 - wiltpoint)
          maxstep = max(100.d0,maxstep)
          y3 = y1 - maxstep
        endif
      else
! ---   take Newton Raphson step
        y3 = y1 - (y2-y1)*Fy1 / (Fy2-Fy1)
        if ((y3-y1) .gt. 0.d0) then
          maxstep = -0.05d0 * y1
          maxstep = max(100.d0,maxstep)
          y3 = min((y1 + maxstep),y3)
        else
          maxstep = 0.05d0 * (y1 - wiltpoint)
          maxstep = max(100.d0,maxstep)
          y3 = max((y1 - maxstep),y3)
        endif
      endif
      y3 = min(y3,hwet)
      y3 = max(y3,wiltpoint)

! --- new estimated value of leaf pressure head is equal to y3!
 300  soil%hleaf = y3
      soil%Hxylem = soil%hleaf + dHPlant

! --- calculate root extraction at y3
      call JongvanLierLoop(state)
      Fy3 = atmo%ptra - soil%qrosum

! --- check convergence
      if (abs(soil%qrosum-atmo%ptra) .lt. cfg_sim%numerical%taccur .or.                             &
     &       abs(y3-y1) .lt. 1.0d0) flconverg = .true.

! --- fatal error if too many iterations
      counter = counter + 1
      if (counter .gt. 1000) then
         messag = '4 Too many iterations for microscopic root'          &
     &             //' water uptake. Please adapt input!'
         call log_warn('rootextraction', messag)
         call fatalerr_collected ('rootextraction',messag)
      endif

! --- apply linear interpolation when Fy1 and Fy3 have opposite sign
      if (.not. flconverg) then
        if (((Fy1 .gt. 0.d0 .and. Fy3 .lt. 0.d0) .or.                   &
     &                     (Fy1 .lt. 0.d0 .and. Fy3 .gt. 0.d0)) .and.   &
     &                      abs(y1-y3) .gt. 1.d0) then

          y3 = y1 + (y3 - y1) * Fy1 / (Fy1 - Fy3)
          y3 = min(y3,hwet)
          y3 = max(y3,wiltpoint)
          goto 300
        endif
      endif

! --- WHILE-DO LOOP TO DETERMINE HLEAF WITHOUT DROUGHT STRESS
      enddo

      soil%qrosum = atmo%ptra
      soil%hleaf = y3
      dHplant = soil%qrosum/kstem
      soil%Hxylem = soil%hleaf + dHPlant

      else
! --- FLSTRESS = TRUE: DETERMINE QROSUM WITH HLEAF = WILTPOINT
      counter = 0

! --- calculate root extraction at previous value of qrosum
      soil%hleaf = wiltpoint
      dHplant = soil%Tactual/kstem
      if (dHplant .gt. (hwet - wiltpoint)) then
        dHplant = 0.98 * (hwet - wiltpoint)
        soil%Tactual = dHplant * kstem
      endif
      soil%Hxylem = soil%hleaf + dHPLant
      call JongvanLierLoop(state)

! --- check convergence
      if (abs(soil%Tactual - soil%qrosum) .lt. cfg_sim%numerical%taccur) flconverg = .true.

! --- set initial values y3 and Fy3
      y3 = soil%Tactual
      Fy3 = soil%Tactual - soil%qrosum

! --- start search on correct value of qrosum
      do while (.not. flconverg)

! --- reset y1 and Fy1
      y1 = y3
      Fy1 = Fy3

      if (soil%Tactual .gt. soil%qrosum) then
! ---   root water extraction less than adopted, decrease tactual
          y2 = y1 - 0.001d0
      else
! ---   root water extraction larger than adopted, increase tactual
          y2 = y1 + 0.001d0
      endif

! --- calculate root extraction at modified pressure head y2 in leaf
      soil%Tactual = y2
      dHplant = soil%Tactual/kstem
      soil%Hxylem = soil%hleaf + dHPLant
      call JongvanLierLoop(state)
      Fy2 = soil%Tactual - soil%qrosum

! --- determine new value of tactual with Newton Raphson algorithm
        y3 = y1 - (y2-y1)*Fy1 / (Fy2-Fy1)
        maxstep = 0.05d0
        if ((y3-y1) .gt. 0.d0) then
          y3 = min((y1 + maxstep),y3)
        else
          y3 = max((y1 - maxstep),y3)
        endif
      dHPlantMax = hwet - wiltpoint
      y3 = min(y3,(dHPlantMax * Kstem))
      y3 = max(y3,0.d0)

! --- new estimated value of tactual is equal to y3!
 400  soil%Tactual = y3
      dHPlant = soil%Tactual/kstem
      soil%Hxylem = soil%hleaf + dHPlant

! --- calculate root extraction at y3
      call JongvanLierLoop(state)
      Fy3 = soil%Tactual - soil%qrosum

! --- check convergence
      if (abs(soil%Tactual - soil%qrosum) .lt. cfg_sim%numerical%taccur) flconverg = .true.

! --- fatal error if too many iterations
      counter = counter + 1
      if (counter .gt. 1000) then
         messag = '5 Too many iterations for microscopic root'          &
     &             //' water uptake. Please adapt input!'
         call log_warn('rootextraction', messag)
         call fatalerr_collected ('rootextraction',messag)
      endif

! --- apply linear interpolation when Fy1 and Fy3 have opposite sign
      if (.not. flconverg) then
        if ((Fy1 .gt. 0.d0 .and. Fy3 .lt. 0.d0) .or.                    &
     &                     (Fy1 .lt. 0.d0 .and. Fy3 .gt. 0.d0)) then

          y3 = y1 + (y3 - y1) * Fy1 / (Fy1 - Fy3)
          goto 400
        endif
      endif

! --- WHILE-DO LOOP TO DETERMINE HLEAF WITH DROUGHT STRESS
      enddo

      dHplant = soil%qrosum/kstem
      soil%Hxylem = soil%hleaf + dHPlant

      endif

! --- qrot(node) has been calculated based on Jong van Lier (2013)

      soil%qrosum = 0.0d0
      do node = 1,noddrz
        soil%qrosum = soil%qrot(node) + soil%qrosum
      enddo

      if ( (atmo%ptra - soil%qrosum) .lt. cfg_sim%numerical%taccur) then
! ---   compensate convergence error
        ratio = atmo%ptra / soil%qrosum
        do node = 1,noddrz
          soil%qrot(node) = soil%qrot(node) * ratio
        enddo
        soil%qrosum = atmo%ptra
        soil%alpJvLier = 1.d0
      else
! ---   drought reduction factor
        soil%alpJvLier = soil%qrosum / atmo%ptra
! ---   calculate potential qrot
        do node = 1,noddrz
          soil%qrot(node) = soil%qrot(node)/soil%alpJvLier
        enddo
      endif

      end associate
      return

  end subroutine JongvanLier

! ----------------------------------------------------------------------
!> Evaluate one microscopic uptake loop for a given xylem/leaf pressure state.
!!
!! @note
!! date: August 2016
!! purpose: Calculate microscopic root water uptake using hleaf
!! @endnote
  subroutine JongvanLierLoop(state)
! ----------------------------------------------------------------------
!     date      : August 2016
!     purpose   : Calculate microscopic root water uptake using hleaf
! ----------------------------------------------------------------------
      use swap_array_dimensions, only: macp
      use swap_log, only: log_warn
      ! [GR-CROP 2026-05-25] retained — dead-branch (swdrought=2 stub-errored in TOML)
      use variables, only: criterhr, flhydrlift, kroot,    &
                           rootcoefa, rooteff, rootradius, &
                           rxylem, stephr, twilt
      implicit none

      type(swap_state_t), intent(inout) :: state

! --- local variables
      integer counter,node,lay
      real(8) x1,x2,x3,Fx1,Fx2,qmax,conducsoil,conducroot
      real(8) step,mflux1,mflux2,depth
      character(len=200) messag

! --- [GR-CROP 2026-05-25] sub-record aliases (crop, soil, atmo, mesh, time).
      associate( &
         crop   => state%crop,                        &
         soil   => state%soilwater,                   &
         atmo   => state%atmosphere,                  &
         mesh   => state%mesh,                        &
         time   => state%timecontrol,                 &
         noddrz => state%crop%common%noddrz           &  ! [GR-CROP 2026-05-25] alias for state read
      )

! --  initialisatie
      soil%qrosum = 0.d0
      x1 = 999.d0
      counter = 0
      ConducRoot = KRoot / rootradius / log(rootradius/rxylem)

      do node = 1,noddrz
! ---   determine h and matricflux potential at root-soil interface
        if (soil%h(node) .gt. -1.d0) then
! ---      very wet conditions
           lay = mesh%layer(node)
           ConducSoil = soil%ksatfit(lay) / (rootcoefa*soil%rmax(node)) /  &
     &                  log(rootcoefa*soil%rmax(node)/rootradius)
           soil%hroot(node) = (ConducSoil*soil%h(node) + ConducRoot*soil%Hxylem) /  &
     &                   (ConducSoil + ConducRoot)
        else
! ---      common conditions
          x3 = soil%hroot(node)
          do while (abs(x3-x1) .gt. CriterHr*log10(max(abs(x3),1.d0)))
            Step = min(abs(x3-x1),StepHr*log10(max(abs(x3),1.d0)))
            x1 = x3
            x2 = x1 - Step
            call matric_flux(x1,node,mflux1,state)
            call matric_flux(x2,node,mflux2,state)
            Fx1 = soil%Hxylem -x1 + soil%rootphi(node)*(soil%mflux(node)-mflux1)
            Fx2 = soil%Hxylem -x2 + soil%rootphi(node)*(soil%mflux(node)-mflux2)
            if (abs(Fx2-Fx1).lt.1.d-10) then
              x3 = x2
            else
              x3 = x1 - (x2-x1)*Fx1 / (Fx2-Fx1)
            endif
            x3 = max(x3,min(soil%Hxylem,soil%h(node)))
            x3 = min(x3,max(soil%Hxylem,soil%h(node)))
! ---       fatal error if too many iterations
            counter = counter + 1
            if (counter .gt. 50) then
              messag = '1 Too many iterations for microscopic root'    &
     &               //' water uptake. Please adapt input!'
              call log_warn('rootextraction', messag)
              call fatalerr_collected ('rootextraction',messag)
            endif
          enddo
          x1 = 999.d0
          counter = 0
          soil%hroot(node) = x3
        endif
        call matric_flux(soil%hroot(node),node,soil%mroot(node),state)
! ---   calculate root water extraction flux
        if (node .lt. noddrz) then
          soil%qrot(node) = rooteff * soil%rootrho(node) *                        &
     &                 (soil%mflux(node)-soil%mroot(node)) * mesh%dz(node)
          if (soil%mflux(node) .gt. soil%mroot(node) ) then
! ---       water extraction, set maximum flux to 10% of available soil water
            qmax = (soil%theta(node) - twilt(node)) * mesh%dz(node) * 0.1d0 / time%dt
            qmax = min(qmax,atmo%ptra)
            soil%qrot(node) = min(soil%qrot(node),qmax)
          else
! ---       possible hydraulic lift, set maximum flux to 0.1% change water content
            if (flhydrlift) then
              qmax = -0.001d0 * mesh%dz(node) / time%dt
              soil%qrot(node) = max(soil%qrot(node),qmax)
            else
! ---         no hydraulic lift allowed
              soil%qrot(node) = 0.d0
              soil%hroot(node) = soil%h(node)
              soil%mroot(node) = soil%mflux(node)
            endif
          endif
        else
! ---     last node, partly filled with roots
          depth = mesh%ztopcp(node) + crop%common%rd
          soil%qrot(node) = rooteff * soil%rootrho(node) *                        &
     &                 (soil%mflux(node)-soil%mroot(node)) * depth
          if (soil%mflux(node) .gt. soil%mroot(node) ) then
! ---       water extraction, set maximum flux to 10% of available soil water
            qmax = (soil%theta(node) - twilt(node)) * depth * 0.1d0 / time%dt
            qmax = min(qmax,atmo%ptra)
            soil%qrot(node) = min(soil%qrot(node),qmax)
          else
! ---       possible hydraulic lift, set maximum flux to 0.1% change water content
            if (flhydrlift) then
              qmax = -0.001d0 * depth / time%dt
              soil%qrot(node) = max(soil%qrot(node),qmax)
            else
! ---         no hydraulic lift allowed
              soil%qrot(node) = 0.d0
              soil%hroot(node) = soil%h(node)
              soil%mroot(node) = soil%mflux(node)
            endif
          endif
        endif
        soil%qrosum = soil%qrot(node) + soil%qrosum
      enddo

      end associate
      return
  end subroutine JongvanLierLoop

end module jongvanlier_dormant_mod
