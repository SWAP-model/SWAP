! to do: default MvG parameters that are read: do they refer to BDENS or Rho_cons; or should BDENS and RhoCons be equal?
! to do: check if Rho_match differs from BDENS or Rho_cons (if not: division by zero possible)

module tillage_mod
   use error_mod, only: fatalerr_collected
   use swap_state_mod, only: swap_state_t

   ! [GR-CROP 2026-05-25] tillage_mod is `use variables`-free.
   ! - till_* Group AB → state%tillage (seeded by state%tillage%init in swap_mod).
   ! - ParamVG → state%soilwater%vg_params_layer(:) (typed per-layer VG
   !   store). Tillage mutates the layer-keyed store and rebuilds per-node
   !   vg_params(:) from it after each event. Legacy paramvg(21, maho)
   !   retired (this file was its last consumer).
   ! [state%cfg-retirement cluster 6] swtill/swhyst/swdiscrvert → state%soilwater%X;
   !   swsolu → state%solute%swsolu. All soil_cfg/solute_cfg associates retired.

   implicit none

!  by default: all in this module is private (local)
   private
!  except for these public routines/functions
   public :: tillage_seed, tillage_step, tillage_output

   contains

   subroutine tillage_seed (state)
   type(swap_state_t), intent(inout)            :: state
   ! local (not to be saved)
   integer                                   :: i
   logical                                   :: fine

   ! Sub-record aliases (canonical associate pattern).
   ! [state%cfg-retirement cluster 6] soil_cfg/solute_cfg aliases dropped;
   ! swtill/swhyst/swdiscrvert now on state%soilwater, swsolu now on state%solute.
   associate( &
      mesh       => state%mesh,             &
      soil       => state%soilwater,        &
      solu       => state%solute,           &
      time       => state%timecontrol,      &
      atmo       => state%atmosphere,       &
      tl         => state%tillage           )

      ! INITIALIZE

      ! some checks: some combinations not (yet) allowed
      if (soil%swtill == 1) then
         if (soil%swhyst == 1)                  call fatalerr_collected ('DoTillage', 'swhyst = 1 not allowed')
         if (solu%swsolu == 1)                  call fatalerr_collected ('DoTillage', 'swsolu = 1 not (yet) allowed')
         if (state%crop%common%swoxygen == 2)   call fatalerr_collected ('DoTillage', 'swoxygen = 2 not (yet) allowed')
         if (soil%flksatexm)                    call fatalerr_collected ('DoTillage', 'flksatexm not (yet) allowed')
         if (soil%swdiscrvert == 1) &
                                                call fatalerr_collected ('DoTillage', 'SwDiscrvert = 1 not (yet) allowed')
      end if

      ! currently: require all Z_tillage = Max_Z_tillage
      do i = 1, tl%Ntill
!         if (dabs(tl%Max_Z_tillage - tl%Z_tillage(i)) > 1.0d-2) call fatalerr ('DoTillage', 'For time being: all Z_tillage must equal Max_Z_tillage')
      end do

      ! determine entry point in tabulated tillage events based on start value t1900; check if input dates are sorted
      call set_iTill (state)

      ! determine number of horizon at depth Max_Z_tillage (MaxNumSoilHo)
      call det_MNSH (state)

      ! for special case i_n_model = 3: calculate slope per soil layer (remains constant over time)
      tl%Slope_match = 0.0d0

      ! TO ADD: CHECK THAT DEPTH OF EACH TILLAGE EVENT CORRESPONDS TO BOTTOM OF SOIL HORIZON; USER MAY NEED TO DEFINE MULTIPLE SUBS-HORIZONS WITHIN A SINGLE REAL SOIL HORIZON
      ! currently: require changes in horizon number (iSoilLayer) at depth Max_Z_tillage
      fine = .false.
      do i = 1, mesh%numlay
         if (mesh%botcom(i) == tl%MaxNumSoilCP) then
            fine = .true.
            exit
         end if
      end do
      if (.not. fine) call fatalerr_collected ('DoTillage', 'Bottom of soil horizon does not coincide with tillage depth(s)')

   end associate
   end subroutine tillage_seed

   subroutine tillage_step (state)
   type(swap_state_t), intent(inout)            :: state
   ! local (not to be saved)
   character(len=20)                         :: STRNG
   logical, parameter                        :: TEST = .false.

   ! Sub-record aliases (canonical associate pattern).
   ! [state%cfg-retirement cluster 6] soil_cfg/solute_cfg aliases dropped; swtill on state%soilwater.
   associate( &
      mesh       => state%mesh,             &
      soil       => state%soilwater,        &
      time       => state%timecontrol,      &
      atmo       => state%atmosphere,       &
      tl         => state%tillage           )

   if (soil%swtill == 0) return      ! no tillage to be considered: return immediately

      ! RATE/STATE EVENT
      tl%Rho_last(1:tl%MaxNumSoilHo) = soil%bdens(1:tl%MaxNumSoilHo)

      ! normal usage
      if (tl%iTill <= tl%Ntill .and. nint(time%t1900) == nint(tl%Date_tillage(tl%iTill))) then
         call DTDPST ("YEAR-MONTHST-DAY", time%t1900, STRNG)
         call Change_Tillage_Info (tl%iTill, state)
         call Change_Bdens(state)
         tl%iTill = tl%iTill + 1
      else
         call Consolidate_Bdens(state)
      end if

      call Change_MvGpars(state)

      Call Adapt_WC_H (TEST, state)

   end associate
   end subroutine tillage_step

   subroutine tillage_output (state)
   type(swap_state_t), intent(inout)            :: state

   ! Sub-record aliases (canonical associate pattern).
   ! [state%cfg-retirement cluster 6] soil_cfg/solute_cfg aliases dropped; swtill on state%soilwater.
   associate( &
      soil       => state%soilwater           )

   if (soil%swtill == 0) return      ! no tillage to be considered: return immediately

      ! OUTPUT
      ! (debug writes to units 222/224/226 removed — dead scaffolding)

   end associate
   end subroutine tillage_output

! **************************************************** Change_MvGpars *********************************************************
   subroutine Change_MvGpars (state)
   implicit none
   type(swap_state_t), intent(inout) :: state
   integer              :: i, node, lay
   integer, parameter   :: Delta = 4
   integer, parameter   :: DeltaMin7 = Delta - 7
   real(8), parameter   :: Omega = -3.97d0
   real(8), parameter   :: Rho_s = 2650d0      ! later as input?
   real(8)              :: wcs_last, Epsilon
   associate( &
      mesh => state%mesh,         &
      soil => state%soilwater,    &
      tl   => state%tillage)
   ! [GR-CROP 2026-05-25] mutate the per-layer VG store
   ! (soil%vg_params_layer) instead of the legacy paramvg(21, maho) global.
   ! Index correspondence (legacy → typed field):
   !   paramvg(1,lay)  → thetar           paramvg(7,lay)  → mpar
   !   paramvg(2,lay)  → thetas           paramvg(8,lay)  → alphaw_sentinel (not mutated)
   !   paramvg(3,lay)  → ksat             paramvg(9,lay)  → h_enpr          (not mutated)
   !   paramvg(4,lay)  → alpha            paramvg(10,lay) → ksatexm         (not mutated)
   !   paramvg(5,lay)  → lpar             (not mutated by tillage)
   !   paramvg(6,lay)  → npar
   do i = 1 , tl%MaxNumSoilHo
      wcs_last = soil%vg_params_layer(i)%thetas       ! help

      ! First, mutate per-layer VG params
      ! thetar (wcr)
      soil%vg_params_layer(i)%thetar = soil%vg_params_layer(i)%thetar * soil%bdens(i)/tl%Rho_last(i)
      ! thetas (wcs)
      soil%vg_params_layer(i)%thetas = soil%vg_params_layer(i)%thetas * (Rho_s - soil%bdens(i))/(Rho_s - tl%Rho_last(i))
      ! ksat (ks,fit)
      soil%vg_params_layer(i)%ksat   = soil%vg_params_layer(i)%ksat   * (soil%vg_params_layer(i)%thetas/wcs_last)**3 * &
                                       (soil%bdens(i)/tl%Rho_last(i))**DeltaMin7
      ! alpha
      soil%vg_params_layer(i)%alpha  = soil%vg_params_layer(i)%alpha  * (soil%bdens(i)/tl%Rho_last(i))**Omega
      ! lambda (lpar): do not change
      ! n (npar)
      select case (tl%i_n_model)
      case (1)
         ! npar: do not change
         continue
      case(2)
         Epsilon = -0.97d0 + 1.28d0 * soil%psilt(i) / soil%pclay(i)
         soil%vg_params_layer(i)%npar = 1.0d0 + (soil%vg_params_layer(i)%npar - 1.0d0) * &
                                        (soil%bdens(i)/tl%Rho_last(i))**Epsilon
      case(3)
         soil%vg_params_layer(i)%npar = dmax1(1.001d0, soil%vg_params_layer(i)%npar + &
                                              (soil%bdens(i) - tl%Rho_last(i)) * tl%Slope_match(i))
      end select

      ! m = 1-1/n
      soil%vg_params_layer(i)%mpar = 1.0d0 - 1.0d0/soil%vg_params_layer(i)%npar
      ! alpha_w / h_enpr / ksatexm: not modified by tillage
   end do

   ! Second: rebuild per-node vg_params(:) from the per-layer store after mutation.
   ! Indices 1-10 are mirrored; relsatthr (11) and ksatthr (12) are not touched
   ! by tillage (legacy paramvg path was the same — they remain at their init values).
   do node = 1, tl%MaxNumSoilCP
      lay = mesh%layer(node)
      soil%vg_params(node)%thetar          = soil%vg_params_layer(lay)%thetar
      soil%vg_params(node)%thetas          = soil%vg_params_layer(lay)%thetas
      soil%vg_params(node)%ksat            = soil%vg_params_layer(lay)%ksat
      soil%vg_params(node)%alpha           = soil%vg_params_layer(lay)%alpha
      soil%vg_params(node)%lpar            = soil%vg_params_layer(lay)%lpar
      soil%vg_params(node)%npar            = soil%vg_params_layer(lay)%npar
      soil%vg_params(node)%mpar            = soil%vg_params_layer(lay)%mpar
      ! alphaw_sentinel: per-node sentinel kept (not driven from per-layer store).
      soil%vg_params(node)%h_enpr          = soil%vg_params_layer(lay)%h_enpr
      soil%vg_params(node)%ksatexm         = soil%vg_params_layer(lay)%ksatexm
   end do
   end associate

   end subroutine Change_MvGpars

! **************************************************** Adapt_WC_H *********************************************************
   subroutine Adapt_WC_H (TEST, state)
   use soilhydraulics_utils, only: watcon, hconduc, prhead
   implicit none

   type(swap_state_t), intent(inout) :: state
   integer                          :: i
   real(8)                          :: sumWCtmin1, sumWCt, dwc, wcr, wcs, summ, dif
   real(8), dimension(state%tillage%MaxNumSoilCP) :: wc
   logical                          :: TEST

   associate( &
      mesh => state%mesh,         &
      soil => state%soilwater,    &
      heat => state%heat,         &
      time => state%timecontrol,  &
      tl   => state%tillage)

   select case (tl%iRedist)
   case (0)
      if (.not.TEST) call fatalerr_collected ('Adapt_WC_H', 'Option iRedist = 0 only allowed in combination with TEST option')
      continue

   case (1)
      ! keep current wc values (wc_new = wc_old) and only change h; exception: when wc_old > wcs_new: alternative redistribution required
      summ = 0.0d0
      do i = 1, tl%MaxNumSoilCP
         wcs = soil%vg_params_layer(mesh%layer(i))%thetas
         if (soil%theta(i) < wcs) then
            soil%h(i) = prhead(mesh%disnod(i), soil%theta(i), soil%h, &
                                          soil%iHWCKmodel(soil%layer(i)), &
                                          i, soil)
         else
            summ = summ + (wcs - soil%theta(i))*mesh%dz(i)
            soil%theta(i) = wcs
            soil%h(i) = 0.0d0
         end if
      end do
      if (summ > 0.0d0) then
         do i = tl%MaxNumSoilCP, 1, -1
            wcs = soil%vg_params_layer(mesh%layer(i))%thetas
            dif = wcs - soil%theta(i)
            if (dif > 0.0d0) then
               if (dif < summ) then
                  soil%theta(i) = wcs
                  summ = summ - dif
               else
                  soil%theta(i) = soil%theta(i) + dif
                  summ = 0.0d0
                  exit
               end if
            end if
         end do
      end if
      soil%pond = summ

   case (2)
      sumWCtmin1 = sum(soil%theta(1:tl%MaxNumSoilCP))
      sumWCt = 0.0d0
      tl%sumDWC = 0.0d0
      do i = 1, tl%MaxNumSoilCP
         wc(i) = watcon(soil%h(i), &
                         soil%vg_params(i), &
                         soil%iHWCKmodel(soil%layer(i)), &
                         i, soil)
         sumWCt = sumWCt + wc(i)
         dwc = soil%theta(i) - wc(i)
         tl%sumDWC = tl%sumDWC + dwc * mesh%dz(i)
      end do

      tl%sumAvail1 = 0.0d0
      tl%sumAvail2 = 0.0d0
      if (sumWCt < sumWCtmin1) then
         ! water to be added; same as sumDWC > 0.0
         do i = 1, tl%MaxNumSoilCP
            wcs = soil%vg_params_layer(mesh%layer(i))%thetas
            tl%sumAvail1 = tl%sumAvail1 + (wcs - wc(i))*mesh%dz(i)
         end do
         do i = 1, tl%MaxNumSoilCP
            wcs = soil%vg_params_layer(mesh%layer(i))%thetas
            if (tl%sumAvail1 > 0.0d0) then
               wc(i) = wc(i) + (wcs - wc(i)) * tl%sumDWC / tl%sumAvail1
               if (wc(i) > wcs) then
                  soil%pond = soil%pond + (wc(i) - wcs) * mesh%dz(i)
                  wc(i) = wcs
               end if
            end if
            soil%h(i)     = prhead(mesh%disnod(i), wc(i), soil%h, &
                                               soil%iHWCKmodel(soil%layer(i)), &
                                               i, soil)
            soil%theta(i) = wc(i)
         end do
      else if (sumWCt > sumWCtmin1) then
         ! water to be removed; same as sumDWC < 0.0
         do i = 1, tl%MaxNumSoilCP
            wcr = soil%vg_params_layer(mesh%layer(i))%thetar
            tl%sumAvail2 = tl%sumAvail2 + (wc(i) - wcr) * mesh%dz(i)
         end do
         do i = 1, tl%MaxNumSoilCP
            wcr = soil%vg_params_layer(mesh%layer(i))%thetar
            wc(i) = wc(i) + (wc(i) - wcr) * tl%sumDWC / tl%sumAvail2
            soil%h(i)     = prhead(mesh%disnod(i), wc(i), soil%h, &
                                               soil%iHWCKmodel(soil%layer(i)), &
                                               i, soil)
            soil%theta(i) = wc(i)
         end do

      endif

   end select

write(124,'(A,1P,12E12.5)') time%date, soil%bdens(1), soil%vg_params_layer(mesh%layer(1))%thetas, soil%theta(1), soil%h(1), &
   hconduc(soil%h(1),soil%theta(1),1.0d0,heat%tsoil(1), &
           soil%vg_params(1), &
           soil%iHWCKmodel(soil%layer(1)), &
           soil%fluseksatexm(1), 1, soil), soil%vg_params_layer(mesh%layer(1))%ksat,    &
   soil%bdens(2), soil%vg_params_layer(mesh%layer(2))%thetas, soil%theta(2), soil%h(2), &
   hconduc(soil%h(2),soil%theta(2),1.0d0,heat%tsoil(2), &
           soil%vg_params(2), &
           soil%iHWCKmodel(soil%layer(2)), &
           soil%fluseksatexm(2), 2, soil), soil%vg_params_layer(mesh%layer(2))%ksat

   end associate
   end subroutine Adapt_WC_H


! **************************************************** Consolidate_Bdens *********************************************************
   subroutine Consolidate_Bdens (state)
   implicit none
   type(swap_state_t), intent(inout) :: state
   integer :: i

   associate( &
      soil => state%soilwater,    &
      atmo => state%atmosphere,   &
      time => state%timecontrol,  &
      tl   => state%tillage)
   if (tl%iTill == 1) return        ! in the beginning before first tillage event: do nothing

   forall (i=1:tl%MaxNumSoilHo) soil%bdens(i) = tl%Rho_cons(i) - (tl%Rho_cons(i) - tl%Rho_last(i)) * dexp(-tl%K_R_cons(i)*atmo%nraida*10.0d0)    ! 10: to transform nraida from cm to mm
   write(123,'(A,1P,10E12.5)') time%date, atmo%nraida, soil%bdens(1:tl%MaxNumSoilHo)
   end associate
   end subroutine Consolidate_Bdens

! **************************************************** Change_Bdens *********************************************************
   subroutine Change_Bdens (state)
   implicit none
   type(swap_state_t), intent(inout) :: state
   integer :: i, NumSoilHo
   ! to check: why this loop to determine NumSoilHo?
   NumSoilHo = 1
   associate( &
      mesh => state%mesh,         &
      soil => state%soilwater,    &
      tl   => state%tillage)
   do i = 2, mesh%numnod
      if (tl%Z_tillage(tl%iTill) > -mesh%zbotcp(i-1) .and. tl%Z_tillage(tl%iTill) <= -mesh%zbotcp(i)) then
         NumSoilHo = mesh%layer(i)
         exit
      end if
   end do
   forall (i=1:NumSoilHo) soil%bdens(i) = tl%Rho_last(i) - tl%I_tillage(tl%iTill) * (tl%Rho_last(i) - tl%Rho_tillage(i))
   end associate

   end subroutine Change_Bdens

! **************************************************** set_iTill *********************************************************
   subroutine set_iTill (state)
   implicit none
   type(swap_state_t), intent(inout) :: state
   integer :: i
   associate( &
      time => state%timecontrol,  &
      tl   => state%tillage)
   tl%iTill = 1
   if (time%t1900 <= tl%Date_tillage(1)) tl%iTill = 1
   do i = 2, tl%Ntill
      if (tl%Date_tillage(i) < tl%Date_tillage(i-1)) call fatalerr_collected ('set_iTill', 'Dates in tabulated tillage events must be sorted')
      ! H-4 bug fix: second comparison was Date_tillage(i-1) (tautological); corrected to Date_tillage(i)
      if (time%t1900 >= tl%Date_tillage(i-1) .and. time%t1900 < tl%Date_tillage(i)) tl%iTill = i-1
   end do
   end associate
   end subroutine set_iTill

! **************************************************** det_MNSH *********************************************************
   subroutine det_MNSH (state)
   implicit none
   type(swap_state_t), intent(inout) :: state
   integer :: i
   associate( &
      mesh => state%mesh,         &
      tl   => state%tillage)
   do i = 2, mesh%numnod
      if (tl%Max_Z_tillage > -mesh%zbotcp(i-1) .and. tl%Max_Z_tillage <= -mesh%zbotcp(i)) then
         tl%MaxNumSoilHo = mesh%layer(i)
         tl%MaxNumSoilCP = i
         exit
      end if
   end do
   end associate
   end subroutine det_MNSH

   subroutine Change_Tillage_Info (iTill_in, state)
   implicit none
   ! global
   integer, intent(in) :: iTill_in
   type(swap_state_t), intent(inout) :: state
   ! local
   integer :: itype, nlay, k

   associate( &
      soil => state%soilwater,    &
      tl   => state%tillage)
   itype                  = tl%Type_Tillage(iTill_in)
   nlay                   = tl%iTT2(itype) - tl%iTT1(itype) + 1
   tl%Rho_tillage(1:nlay) = tl%TAB_Rho_tillage(tl%iTT1(itype):tl%iTT2(itype))
   tl%Rho_cons(1:nlay)    = tl%TAB_Rho_cons(tl%iTT1(itype):tl%iTT2(itype))
   tl%K_R_cons(1:nlay)    = tl%TAB_K_R_cons(tl%iTT1(itype):tl%iTT2(itype))
   tl%Rho_match(1:nlay)   = tl%TAB_Rho_match(tl%iTT1(itype):tl%iTT2(itype))
   tl%N_match(1:nlay)     = tl%TAB_N_match(tl%iTT1(itype):tl%iTT2(itype))
   do k = 1, nlay
      tl%Slope_match(k) = (soil%vg_params_layer(k)%npar - tl%N_match(k)) / (tl%Rho_cons(k) - tl%Rho_match(k))
   end do
   end associate
   end subroutine Change_Tillage_Info

end module tillage_mod
