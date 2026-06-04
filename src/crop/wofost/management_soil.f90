! File VersionID:
!   $Id: management_soil.f90 371 2018-02-19 09:25:14Z heine003 $
! ----------------------------------------------------------------------
!> Soil management module.
!!
!! Hosts `SoilManagement`, the legacy task-based routine for soil nutrient
!! management: per-step bookkeeping, mineralisation / nitrification /
!! denitrification, crop-residue amendments, nutrient-balance output,
!! and end-of-run dump.
!!
!! @note
!! The physics is unchanged from the legacy implementation (March 2015).
!! Cleanup so far:
!!   - Three write-only `state%nutrients%X = X` mirror blocks deleted
!!     (state%nutrients is write-only across the whole codebase; the
!!     mirrors were a no-op).
!!   - Sub-record ASSOCIATE pattern applied; `state%X%Y` chains aliased
!!     to short canonical names (time/soil/surf/crop/atmo/heat/mesh).
!!
!! Strangler-fig follow-up:
!!   The WSN compute path still reads/writes the bare-global pools in
!!   `Wofost_Soil_Declarations` (FOM_t, Bio_t, Hum_t, cNH4_t, cNO3_t, …).
!!   Retiring those globals onto `state%nutrients` is a separate arc and
!!   needs all of `wofost_soil_*.f90` ported alongside this driver.
module management_soil_mod
   use error_mod, only: fatalerr_collected

   implicit none
   private
   public :: SoilManagement

contains

!> Execute soil management task.
!!
!! Task-switching wrapper around the legacy nutrient implementation.
!!
!! @param[in] task Selector for init / per-step / amendment /
!!                 mineralisation / crop-residue / output / close.
   subroutine SoilManagement(task, state)
      use swap_array_dimensions, only: masme
      use Wofost_Soil_Declarations
      use Wofost_Soil_Interface
      use file_io_mod, only: file_open
      use swap_state_mod, only: swap_state_t

      implicit none

      type(swap_state_t), intent(inout) :: state

      character(len=300) :: filnam
      integer            :: task, sme, smm, nsmm, nsme
      integer            :: i, j, idum
      real(8)            :: dum, t4, t3, help
      real(8)            :: MatAmount(masme) ! amount of applied material (kg/ha)
      real(8)            :: smedate(masme)   ! soilmanagement event dates (-)

      integer :: le, fn
      real(8) :: ProdRate0, ProducPot, ProducAct, TCSF
      real(8) :: cNH4_t_Ndemand_rate_limited
      real(8) :: cNH4_av_Ndemand_rate_limited
      real(8) :: cNH4_t_Nsupply_rate_limited
      real(8) :: cNH4_av_Nsupply_rate_limited
      real(8) :: cNO3_t_Ndemand_rate_limited
      real(8) :: cNO3_av_Ndemand_rate_limited
      real(8) :: cNO3_t_Nsupply_rate_limited
      real(8) :: cNO3_av_Nsupply_rate_limited
      real(8) :: Nsupply_Ndemand_rate_limited
      real(8) :: Nsupply_Nsupply_rate_limited
      real(8) :: dum1, dum2, dum3, dum4
      character(len=1) :: comma

      real(8) :: xAmend, xAppAge, xOrgMatFrac, xOrgNFrac, xNH4NFrac, xNO3NFrac

      real(8) :: FactNuptJuvenil
      logical :: FlNuptJuvenil

      integer :: swexpertN
      real(8) :: AppAgeArableRt, AppAgeArableLv
      real(8) :: AppAgeArableSt, AppAgeArableSo
      real(8) :: AppAgeGrassRt, AppAgeGrassLv, AppAgeGrassSt

      associate (time => state%timecontrol,    &
                 soil => state%soilwater,      &
                 surf => state%surfacewater,   &
                 crop => state%crop,           &
                 atmo => state%atmosphere,     &
                 heat => state%heat,           &
                 mesh => state%mesh)

      select case (task)
      case (1)
         ! Legacy nutrient soil-management init (file-open + state
         ! initialisation) retired with the legacy readers. Initial
         ! soil-nutrient state now flows from the [nutrients] typed
         ! config via state%nutrients%init (ADR 0026 N2a).
         return

      case (2)
         ! Save start-of-timestep WSN snapshots.
         do fn = 1, nf
            FOM_t0(fn) = FOM_t(fn)
         end do
         Bio_t0    = Bio_t
         Hum_t0    = Hum_t
         cNH4_t0   = cNH4_t
         cNO3_t0   = cNO3_t
         WFrac_t0  = WFrac_t
         t1900Soil = time%t1900
         return

      case (3)
         ! Timed soil-management events (amendments).
         if (abs(state%nutrients%timeamend(state%nutrients%isme) + 1.0d0 - time%t1900) .lt. 1.d-3) then
            call Wofost_SoilAmendents(state)
            state%nutrients%isme = state%nutrients%isme + 1
         endif
         return

      case (4)
         ! Mineralisation (production of mineral nitrogen).

         dum1 = 0.0d0; dum2 = 0.0d0; dum3 = 0.0d0; dum4 = 0.0d0
         idum = 0
         do i = 1, mesh%numnod
            if (dum1 + 1.0d-2 * mesh%dz(i) .lt. dz_WSN) then
               idum = idum + 1
               dum1 = dum1 + 1.0d-2 * mesh%dz(i)
               dum2 = dum2 + heat%tsoil(i) * 1.0d-2 * mesh%dz(i)
               dum4 = dum4 + soil%theta(i)  * 1.0d-2 * mesh%dz(i)
            end if
         end do
         Temp     = dum2 / dum1
         WFrac_t  = dum4 / dum1
         dt_WSN   = time%t1900 - t_WSNold
         t_WSNold = time%t1900

         call Wofost_SoilRateConstants(1)
         call Wofost_SoilOrgMatN

         ! --- transformation and transport processes of soluble nitrogen ---

         ! Results from potential crop growth.
         Ndemand = 1.0d-4 * NdemandSoil

         ! Water-balance items of the WSN soil layer.
         dum1 = 0.0d0; dum2 = 0.0d0; dum3 = 0.0d0; dum4 = 0.0d0
         idum = 0
         do i = 1, mesh%numnod
            dum2 = dum2 + soil%inqrot(i) / time%outper
            if (dum1 + 1.0d-2 * mesh%dz(i) .lt. dz_WSN) then
               dum1 = dum1 + 1.0d-2 * mesh%dz(i)
               do le = 1, 5
                  dum3 = dum3 - min(0.0d0, (surf%inqdra(le, i) / time%outper))
               end do
            end if
         end do
         help = 1.0d-2 * (atmo%intr%igrai + atmo%intr%isnrai + atmo%intr%igsnow &
                          + soil%igird - soil%iintc + soil%irunon - soil%iruno) / time%outper
         SoilEvap     = 1.0d-2 * atmo%intr%ievap / time%outper
         Wflux_inTop  = help
         Wflux_inLat  = 1.0d-2 * dum3
         Wflux_transp = 1.0d-2 * dum2
         dum4         = help + 1.0d-2 * (dum3 - dum2) + &
                        (WFrac_t0 - WFrac_t) * dz_WSN / dt_WSN - SoilEvap
         Wflux_out    = max(0.0d0, dum4)
         Wflux_inBot  = -min(0.0d0, dum4)

         call Wofost_SoilRateConstants(2)

         ! --- ammonium ---

         ! Mineralisation; Nminer is negative under immobilisation.
         ProdRate0 = Nminer / dt_WSN

         ! Crop uptake.
         FactNuptJuvenil = 0.0d0
         FlNuptJuvenil   = .false.
         if (crop%common%flCropCalendar .and. crop%common%dvs .lt. 1.0d0 .and. &
             LaiCritNupt .gt. 1.0d-02 .and. crop%lai .lt. LaiCritNupt) then
            FactNuptJuvenil = (LaiCritNupt - crop%lai) / LaiCritNupt
            FlNuptJuvenil   = .true.
         end if

         if (FlNuptJuvenil) then
            ProducPot = ProdRate0 - FactNuptJuvenil * Ndemand / dz_WSN
            TCSF = TCSF_N * (0.5d0*(WFrac_t + WFrac_t0) + DryBD*SorpCoef) /  &
                            (0.5d0*(WFrac_t + WFrac_t0)) *                    &
                            (1.0d0 - FactNuptJuvenil)
            ! Rooting-depth limitation.
            TCSF = TCSF * max(1.0d0, (0.01d0 * crop%common%rd / dz_WSN))
            call Wofost_SoilWaterN(dz_WSN, dt_WSN, WFrac_t, WFrac_t0,           &
                                   Wflux_out, Wflux_transp, Wflux_inBot, Wflux_inTop, &
                                   Wflux_inLat, TCSF, RateConNitrif, ProducPot,      &
                                   ProducAct, DryBD, SorpCoef, cNH4N_seep, cNH4N_top,&
                                   cNH4N_lat, cNH4_t0, cNH4_t, cNH4_av)
            if (cNH4_t .ge. 0.0d0) then
               NsupplyNH4N = TCSF * Wflux_transp * cNH4_av / dz_WSN + ProdRate0 - ProducAct
            else
               ProducPot = ProdRate0
               TCSF = TCSF_N * (0.5d0*(WFrac_t + WFrac_t0) + DryBD*SorpCoef) /  &
                               (0.5d0*(WFrac_t + WFrac_t0))
               TCSF = TCSF * max(1.0d0, (0.01d0 * crop%common%rd / dz_WSN))
               call Wofost_SoilWaterN(dz_WSN, dt_WSN, WFrac_t, WFrac_t0,           &
                                      Wflux_out, Wflux_transp, Wflux_inBot, Wflux_inTop, &
                                      Wflux_inLat, TCSF, RateConNitrif, ProducPot,      &
                                      ProducAct, DryBD, SorpCoef, cNH4N_seep, cNH4N_top,&
                                      cNH4N_lat, cNH4_t0, cNH4_t, cNH4_av)
               NsupplyNH4N = TCSF * Wflux_transp * cNH4_av / dz_WSN
            end if

         else

            ! NH4N demand rate limiting.
            ProducPot = ProdRate0 - Ndemand / dz_WSN
            TCSF      = 0.0d0
            call Wofost_SoilWaterN(dz_WSN, dt_WSN, WFrac_t, WFrac_t0,           &
                                   Wflux_out, Wflux_transp, Wflux_inBot, Wflux_inTop, &
                                   Wflux_inLat, TCSF, RateConNitrif, ProducPot,      &
                                   ProducAct, DryBD, SorpCoef, cNH4N_seep, cNH4N_top,&
                                   cNH4N_lat, cNH4_t0, cNH4_t, cNH4_av)
            cNH4_t_Ndemand_rate_limited  = cNH4_t
            cNH4_av_Ndemand_rate_limited = cNH4_av
            Nsupply_Ndemand_rate_limited = ProdRate0 - ProducAct

            ! NH4N supply rate limiting.
            ProducPot = ProdRate0
            TCSF = TCSF_N * (0.5d0*(WFrac_t + WFrac_t0) + DryBD*SorpCoef) /  &
                            (0.5d0*(WFrac_t + WFrac_t0))
            ! Rooting-depth limitation.
            TCSF = TCSF * max(1.0d0, (0.01d0 * crop%common%rd / dz_WSN))
            call Wofost_SoilWaterN(dz_WSN, dt_WSN, WFrac_t, WFrac_t0,           &
                                   Wflux_out, Wflux_transp, Wflux_inBot, Wflux_inTop, &
                                   Wflux_inLat, TCSF, RateConNitrif, ProducPot,      &
                                   ProducAct, DryBD, SorpCoef, cNH4N_seep, cNH4N_top,&
                                   cNH4N_lat, cNH4_t0, cNH4_t, cNH4_av)
            cNH4_t_Nsupply_rate_limited  = cNH4_t
            cNH4_av_Nsupply_rate_limited = cNH4_av
            Nsupply_Nsupply_rate_limited = TCSF * Wflux_transp * cNH4_av / dz_WSN

            if (cNH4_t_Ndemand_rate_limited .lt. cNH4_t_Nsupply_rate_limited) then
               cNH4_t      = cNH4_t_Nsupply_rate_limited
               cNH4_av     = cNH4_av_Nsupply_rate_limited
               NsupplyNH4N = Nsupply_Nsupply_rate_limited
            else
               cNH4_t      = cNH4_t_Ndemand_rate_limited
               cNH4_av     = cNH4_av_Ndemand_rate_limited
               NsupplyNH4N = Nsupply_Ndemand_rate_limited
            end if

         end if

         ! --- nitrification ---
         ProdRate0 = 0.5d0 * (WFrac_t + WFrac_t0) * RateConNitrif * cNH4_av

         if (FlNuptJuvenil) then
            ProducPot = ProdRate0 - FactNuptJuvenil * (Ndemand / dz_WSN - NsupplyNH4N)
            TCSF      = TCSF_N * (1.0d0 - FactNuptJuvenil)
            ! Rooting-depth limitation.
            TCSF = TCSF * max(1.0d0, (0.01d0 * crop%common%rd / dz_WSN))
            call Wofost_SoilWaterN(dz_WSN, dt_WSN, WFrac_t, WFrac_t0,           &
                                   Wflux_out, Wflux_transp, Wflux_inBot, Wflux_inTop, &
                                   Wflux_inLat, TCSF, RateConDenitr, ProducPot,      &
                                   ProducAct, 0.0d0, 0.0d0, cNO3N_seep, cNO3N_top,   &
                                   cNO3N_lat, cNO3_t0, cNO3_t, cNO3_av)
            if (cNO3_t .ge. 0.0d0) then
               NsupplyNO3N = TCSF * Wflux_transp * cNO3_av / dz_WSN + ProdRate0 - ProducAct
            else
               ProducPot = ProdRate0
               TCSF      = TCSF_N
               TCSF = TCSF * max(1.0d0, (0.01d0 * crop%common%rd / dz_WSN))
               call Wofost_SoilWaterN(dz_WSN, dt_WSN, WFrac_t, WFrac_t0,           &
                                      Wflux_out, Wflux_transp, Wflux_inBot, Wflux_inTop, &
                                      Wflux_inLat, TCSF, RateConDenitr, ProducPot,      &
                                      ProducAct, 0.0d0, 0.0d0, cNO3N_seep, cNO3N_top,   &
                                      cNO3N_lat, cNO3_t0, cNO3_t, cNO3_av)
               NsupplyNO3N = TCSF * Wflux_transp * cNO3_av / dz_WSN
            end if

         else

            ! NO3N demand rate limiting.
            ProducPot = ProdRate0 - (Ndemand / dz_WSN - NsupplyNH4N)
            TCSF      = 0.0d0
            call Wofost_SoilWaterN(dz_WSN, dt_WSN, WFrac_t, WFrac_t0,           &
                                   Wflux_out, Wflux_transp, Wflux_inBot, Wflux_inTop, &
                                   Wflux_inLat, TCSF, RateConDenitr, ProducPot,      &
                                   ProducAct, 0.0d0, 0.0d0, cNO3N_seep, cNO3N_top,   &
                                   cNO3N_lat, cNO3_t0, cNO3_t, cNO3_av)
            cNO3_t_Ndemand_rate_limited  = cNO3_t
            cNO3_av_Ndemand_rate_limited = cNO3_av
            Nsupply_Ndemand_rate_limited = ProdRate0 - ProducAct

            ! NO3N supply rate limiting.
            ProducPot = ProdRate0
            TCSF      = TCSF_N
            ! Rooting-depth limitation.
            TCSF = TCSF * max(1.0d0, (0.01d0 * crop%common%rd / dz_WSN))
            call Wofost_SoilWaterN(dz_WSN, dt_WSN, WFrac_t, WFrac_t0,           &
                                   Wflux_out, Wflux_transp, Wflux_inBot, Wflux_inTop, &
                                   Wflux_inLat, TCSF, RateConDenitr, ProducPot,      &
                                   ProducAct, 0.0d0, 0.0d0, cNO3N_seep, cNO3N_top,   &
                                   cNO3N_lat, cNO3_t0, cNO3_t, cNO3_av)
            cNO3_t_Nsupply_rate_limited  = cNO3_t
            cNO3_av_Nsupply_rate_limited = cNO3_av
            Nsupply_Nsupply_rate_limited = TCSF * Wflux_transp * cNO3_av / dz_WSN

            if (cNO3_t_Ndemand_rate_limited .lt. cNO3_t_Nsupply_rate_limited) then
               cNO3_t      = cNO3_t_Nsupply_rate_limited
               cNO3_av     = cNO3_av_Nsupply_rate_limited
               NsupplyNO3N = Nsupply_Nsupply_rate_limited
            else
               cNO3_t      = cNO3_t_Ndemand_rate_limited
               cNO3_av     = cNO3_av_Ndemand_rate_limited
               NsupplyNO3N = Nsupply_Ndemand_rate_limited
            end if

         end if

         Nsupply = (NsupplyNH4N + NsupplyNO3N) * dz_WSN

         call Wofost_SoilBalanceCheck

         NsupplySoil = 1.0d+04 * Nsupply
         return

      case (5)
         ! Crop-residue amendments for arable crops.
         !
         ! Material codes (AppAge / OrgMatFrac / OrgNFrac / NH4NFrac / NO3NFrac):
         !   11 'Green leaves'              0.92  1.0  t.b.d.  0.0    0.0
         !   12 'Overground crop residues'  0.99  1.0  t.b.d.  0.0    0.0
         !   13 'Root and stubble residues' 1.57  1.0  t.b.d.  0.0    0.0
         !   14 'Grass shoots'              0.92  1.0  t.b.d.  0.0    t.b.d.
         !   15 'Grass roots'               1.20  1.0  t.b.d.  0.0    t.b.d.
         !   16 'Tree leaves'               2.25  1.0  t.b.d.  0.0    0.0
         !   17 'Spruce needles'            3.34  1.0  t.b.d.  0.0    0.0
         !
         ! The grassland branch below the arable block stays commented out
         ! until grassland nutrient coupling is implemented.

         if (idwrt .le. 1.0d-8 .and. idwst .le. 1.0d-8 .and.  &
             idwlv .le. 1.0d-8 .and. idwso .le. 1.0d-8) return

         xOrgMatFrac = 1.0d0
         xNH4NFrac   = 0.0d0
         xNO3NFrac   = 0.0d0
         if (idwrt .gt. 1.0d-8) then
            xAmend    = 1.0d-4 * idwrt
            xAppAge   = AppAgeArableRt
            xOrgNFrac = iNLOSSR / idwrt
            call Wofost_Soil_CropResidues(xAmend, xAppAge, &
                                          xOrgMatFrac, xOrgNFrac, xNH4NFrac, xNO3NFrac)
         end if
         if (idwlv .gt. 1.0d-8) then
            xAmend    = 1.0d-4 * idwlv
            xAppAge   = AppAgeArableLv
            xOrgNFrac = iNLOSSL / idwlv
            call Wofost_Soil_CropResidues(xAmend, xAppAge, &
                                          xOrgMatFrac, xOrgNFrac, xNH4NFrac, xNO3NFrac)
         end if
         if (idwst .gt. 1.0d-8) then
            xAmend    = 1.0d-4 * idwst
            xAppAge   = AppAgeArableSt
            xOrgNFrac = iNLOSSS / idwst
            call Wofost_Soil_CropResidues(xAmend, xAppAge, &
                                          xOrgMatFrac, xOrgNFrac, xNH4NFrac, xNO3NFrac)
         end if
         if (idwso .gt. 1.0d-8) then
            xAmend    = 1.0d-4 * idwso
            xAppAge   = AppAgeArableSo
            xOrgNFrac = iNLOSSO / idwso
            call Wofost_Soil_CropResidues(xAmend, xAppAge, &
                                          xOrgMatFrac, xOrgNFrac, xNH4NFrac, xNO3NFrac)
         end if

! ----------------------------------------------------------------------
! Grassland-detailed branch — kept commented until grassland nutrient
! coupling is implemented (CropType(icrop).eq.3).
!
!         xOrgMatFrac = 1.0d0
!         xNH4NFrac   = 0.0d0
!         if (idwrt .gt. 1.0d-8) then
!            xAmend    = 1.0d-4 * idwrt
!            xAppAge   = AppAgeGrassRt
!            xOrgNFrac = iNLOSSR / idwrt
!            xNO3NFrac = 0.0d0
!            call Wofost_Soil_CropResidues(xAmend, xAppAge, &
!                                          xOrgMatFrac, xOrgNFrac, xNH4NFrac, xNO3NFrac)
!         end if
!         if (idwlv .gt. 1.0d-8) then
!            xAmend    = 1.0d-4 * idwlv
!            xAppAge   = AppAgeGrassLv
!            xOrgNFrac = iNLOSSL / idwlv
!            xNO3NFrac = 0.0d0
!            call Wofost_Soil_CropResidues(xAmend, xAppAge, &
!                                          xOrgMatFrac, xOrgNFrac, xNH4NFrac, xNO3NFrac)
!         end if
!         if (idwst .gt. 1.0d-8) then
!            xAmend    = 1.0d-4 * idwst
!            xAppAge   = AppAgeGrassSt
!            xOrgNFrac = iNLOSSS / idwst
!            xNO3NFrac = 0.0d0
!            call Wofost_Soil_CropResidues(xAmend, xAppAge, &
!                                          xOrgMatFrac, xOrgNFrac, xNH4NFrac, xNO3NFrac)
!         end if
! ----------------------------------------------------------------------

         idwrt = 0.0d0
         idwlv = 0.0d0
         idwst = 0.0d0
         idwso = 0.0d0

         iNLOSSL = 0.0d0
         iNLOSSR = 0.0d0
         iNLOSSS = 0.0d0
         iNLOSSO = 0.0d0
         return

      case (6)
         ! Output and writing.

         comma   = ','
         t4      = 1.0d+4
         t3      = 1.0d+3
         pnratio = 0.1d0

         Ntotuptake = Ntotuptake + t4 * (NH4_upt + NO3_upt)
         Ptotuptake = Ptotuptake + pnratio * t4 * (NH4_upt + NO3_upt)
         DMcressur  = DMcressur  + idwlv_1 + idwst_1 + idwso_1
         Ncressurf  = Ncressurf  + iNLOSSL_1 + iNLOSSS_1 + iNLOSSO_1
         Pcressurf  = Pcressurf  + pnratio * (iNLOSSL_1 + iNLOSSS_1 + iNLOSSO_1)
         DMcresbott = DMcresbott + idwrt_1
         Ncresbott  = Ncresbott  + iNLOSSR_1
         Pcresbott  = Pcresbott  + pnratio * iNLOSSR_1

         write (nut, '(a11,a1,i3,a1,i6,99(a1,1pe11.4:))')                              &
              time%date, comma, time%daynr, comma, time%daycum, comma,                 &
              t4*FOM_old, comma, t4*FOM_end, comma, t4*(FOM_old-FOM_end), comma,       &
              t4*FOM_add, comma, t4*FOM_cres, comma,                                   &
              t4*FOM2Bio, comma, t4*FOM2Hum, comma, t4*FOM_dis, comma,                 &
              t4*Bio_old, comma, t4*Bio_end, comma, t4*(Bio_old-Bio_end), comma,       &
              t4*Bio2Bio, comma, t4*Bio2Hum, comma, t4*Bio_dis, comma,                 &
              t4*Hum_old, comma, t4*Hum_end, comma, t4*(Hum_old-Hum_end), comma,       &
              t4*Hum_add, comma, t4*Hum_cres, comma,                                   &
              t4*Hum2Bio, comma, t4*Hum2Hum, comma, t4*Hum_dis, comma,                 &
              t4*cDissi, comma, t4*NFOM_old, comma, t4*NFOM_end, comma,                &
              t4*(NFOM_old-NFOM_end), comma, t4*NFOM_add, comma,                       &
              t4*NFOM_cres, comma, t4*NFOM2Bio, comma, t4*NFOM2Hum, comma,             &
              t4*NFOM_min, comma, t4*NBio_old, comma, t4*NBio_end, comma,              &
              t4*(NBio_old-NBio_end), comma,                                           &
              t4*NBio2Bio, comma, t4*NBio2Hum, comma, t4*NBio_min, comma,              &
              t4*NHum_old, comma, t4*NHum_end, comma,                                  &
              t4*(NHum_old-NHum_end), comma, t4*NHum_add, comma,                       &
              t4*NHum_cres, comma, t4*NHum2Bio, comma, t4*NHum2Hum, comma,             &
              t4*NHum_min, comma, t4*NH4_miner, comma,                                 &
              t4*NH4_old, comma, t4*NH4_end, comma, t4*(NH4_old-NH4_end), comma,       &
              t4*NH4N_amend, comma, t4*NH4N_cres, comma, t4*NH4_intop, comma,          &
              t4*NH4_inlat, comma,                                                     &
              t4*NH4_inbot, comma, t4*NH4_upt, comma, t4*NH4_out, comma,               &
              t4*NH4N_volat, comma, t4*NH4_nitrif, comma,                              &
              t4*NO3_old, comma, t4*NO3_end, comma, t4*(NO3_old-NO3_end), comma,       &
              t4*NO3N_amend, comma, t4*NO3N_cres, comma, t4*NO3_intop, comma,          &
              t4*NO3_inlat, comma, t4*NO3_inbot, comma, t4*NO3_upt, comma,             &
              t4*NO3_out, comma, t4*NO3_denitr, comma, t4*Ndemand, comma,              &
              t4*Nsupply, comma,                                                       &
              t3*WFrac_t0*dz_WSN, comma, t3*WFrac_t*dz_WSN, comma,                     &
              t3*(WFrac_t0 - WFrac_t)*dz_WSN, comma,                                   &
              t3*Wflux_inTop*dt_WSN, comma,                                            &
              t3*Wflux_inLat*dt_WSN, comma,                                            &
              t3*Wflux_inBot*dt_WSN, comma, t3*SoilEvap*dt_WSN, comma,                 &
              t3*Wflux_transp*dt_WSN, comma, t3*Wflux_out*dt_WSN, comma,               &
              idwrt_1, comma, idwlv_1, comma, idwst_1, comma,                          &
              iNLOSSL_1, comma, iNLOSSR_1, comma, iNLOSSS_1, comma,                    &
              idwso_1, comma, iNLOSSO_1, comma,                                        &
              WFPS, comma, red_T, comma, red_W, comma, red_W_Nit, comma,               &
              red_W_Den, comma, red_Resp

         ! Reset OM/N balance accumulators for the next output period.
         FOM_add  = 0.0d0
         FOM_cres = 0.0d0
         FOM2Bio  = 0.0d0
         FOM2Hum  = 0.0d0
         FOM_dis  = 0.0d0

         Bio2Bio = 0.0d0
         Bio2Hum = 0.0d0
         Bio_dis = 0.0d0

         Hum_add  = 0.0d0
         Hum_cres = 0.0d0
         Hum2Bio  = 0.0d0
         Hum2Hum  = 0.0d0
         Hum_dis  = 0.0d0

         NFOM_add  = 0.0d0
         NFOM_cres = 0.0d0
         NFOM2Bio  = 0.0d0
         NFOM2Hum  = 0.0d0
         NFOM_min  = 0.0d0

         NBio2Bio = 0.0d0
         NBio2Hum = 0.0d0
         NBio_min = 0.0d0

         NHum_add  = 0.0d0
         NHum_cres = 0.0d0
         NHum2Bio  = 0.0d0
         NHum2Hum  = 0.0d0
         NHum_min  = 0.0d0

         NH4N_amend = 0.0d0
         NO3N_amend = 0.0d0
         NH4N_cres  = 0.0d0
         NO3N_cres  = 0.0d0
         NH4N_volat = 0.0d0

         ! Snapshot crop-residue DM / N losses for next-period reporting.
         idwrt_1   = idwrt
         idwlv_1   = idwlv
         idwst_1   = idwst
         idwso_1   = idwso
         iNLOSSL_1 = iNLOSSL
         iNLOSSR_1 = iNLOSSR
         iNLOSSS_1 = iNLOSSS
         iNLOSSO_1 = iNLOSSO

         ! CROP_EXT file output.
         if (flCropExt .and. time%floutput) then
            write (cropext, '(a10,8(a1,1pe11.4))')                                     &
                 time%date, comma, Ntotuptake, comma, Ptotuptake, comma, DMcressur,     &
                 comma, Ncressurf, comma, Pcressurf, comma, DMcresbott,                 &
                 comma, Ncresbott, comma, Pcresbott

            if (DMcressur .ge. 1.0d-06) then
               nishmi = min(nishmi, Ncressurf/DMcressur)
               nishma = max(nishma, Ncressurf/DMcressur)
               poshmi = min(poshmi, Pcressurf/DMcressur)
               poshma = max(poshma, Pcressurf/DMcressur)
            end if
            if (DMcresbott .ge. 1.0d-6) then
               niromi = min(niromi, Ncresbott/DMcresbott)
               niroma = max(niroma, Ncresbott/DMcresbott)
               poromi = min(poromi, Pcresbott/DMcresbott)
               poroma = max(poroma, Pcresbott/DMcresbott)
            end if
            Ntotuptake = 0.0d0
            Ptotuptake = 0.0d0
            DMcressur  = 0.0d0
            Ncressurf  = 0.0d0
            Pcressurf  = 0.0d0
            DMcresbott = 0.0d0
            Ncresbott  = 0.0d0
            Pcresbott  = 0.0d0
         end if
         return

      case (7)
         ! Legacy diagnostic dump retired (ADR 0026): the old path read
         ! <project>.snp as a template and emitted <project>_nut.end. The
         ! .snp files do not exist in the modern flow, so the read would
         ! fail at runtime if flCropNut were ever true via this path.
         ! A clean CSV-style pool dump can be added when someone asks.
         return

      case default
         call state%diag%fatal('SoilManagement', 'Illegal value for TASK')
      end select

      end associate
   end subroutine SoilManagement

end module management_soil_mod
