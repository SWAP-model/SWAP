! File VersionID:
!   $Id: soilwater.f90 369 2018-01-22 13:18:25Z heine003 $
! ----------------------------------------------------------------------
      subroutine soilwater(task) 
! ----------------------------------------------------------------------
!     Date               : Aug 2004   
!     Purpose            : calculate soil water state variables
! ----------------------------------------------------------------------
     use doln
      use Variables
      implicit none

! --- local variables
      integer task,lay,node,i,j

      real(8) tab(mabbc*2),afgen,watcon,moiscap,hconduc,hcomean
      character(len=200) messag

      select case (task)
      case (1)

! === initialize Soilwater rate/state variables ========================

! --- initialize miscellaneous
      hatm = -2.75d+05
      nraidt = 0.0d0
      nird = 0.0d0
      if (swinco.ne.3) then
        ldwet = 0.d0
        spev = 0.d0
        saev = 0.d0
      endif
      runon = 0.0d0
      qtop = 0.d0
      do i = 1,numnod+1
        q(i) = 0.0d0
      enddo
      do i = 1,macp
        evp(i) = 0.0d0
        rfcp(i) = 1.0d0
      enddo
      vtair = 0.0d0
      cQMpLatSs = 0.0d0

! --- Soil physics: tabulated or MualemVanGenuchten functions
      cofgen = 0.0d0
      if(swsophy.eq.1) then
!       tabulated functions (h,theta,k,dthetadh,dkdtheta) tabulated
        do node = 1,numnod
          numtab(node) = numtablay(layer(node))
          do i=0,matabentries
             ientrytab(node,i) = ientrytablay(layer(node),i) 
          end do
        end do
        do node = 1,numnod
          do i = 1,7
            do j = 1, numtab(node)
              sptab(i,node,j) = sptablay(i,layer(node),j)
            end do
          end do
!         assign values to cofgen
          cofgen(1,node) = 0.0d0                          ! thetar
          cofgen(2,node) = sptab(2,node,numtab(node))     ! thetas
          cofgen(3,node) = sptab(3,node,numtab(node))     ! ksat
          if (do_ln_trans) cofgen(3,node) = dexp(cofgen(3,node))
        end do
        do lay = 1,numlay
          ksatfit(lay) = cofgen(3,nod1lay(lay))
          thetsl(lay) = cofgen(2,nod1lay(lay))
        end do
      else
!       MvanG functions 
        do node = 1,numnod
          lay = layer(node)
          do i = 1, 10
            cofgen(i,node) = paramvg(i,lay)
          end do
!         assign dummy value to alphaw
          cofgen(8,node) = -9999.9d0
          if (cofgen(10,node) > 0.0d0) fluseksatexm(node) = .true.
          cofgen(11,node) = relsatthr(lay)
          cofgen(12,node) = ksatthr(lay)
          if (iHWCKmodel(lay) ==  3 .OR. iHWCKmodel(lay) ==  6 .OR. iHWCKmodel(lay) ==  7 .OR. &
              iHWCKmodel(lay) == 10 .OR. iHWCKmodel(lay) == 11) then
             cofgen(13:17,node) = paramvg(13:17,lay)
          end if
          if (iHWCKmodel(lay) ==  5 .OR. iHWCKmodel(lay) ==  7) then
             cofgen(18,node) = paramvg(18,lay)
          end if
          if (iHWCKmodel(lay) ==  8 .OR. iHWCKmodel(lay) ==  9 .OR. &
              iHWCKmodel(lay) == 10 .OR. iHWCKmodel(lay) == 11) then
             cofgen(18:21,node) = paramvg(18:21,lay)
          end if
        end do
        thetsl = 0.0d0
        do lay = 1, numlay
          thetsl(lay) = paramvg(2,lay)
        end do
      endif

! --- saturated and residual watercontent of each node; hysteresis parameters
      do node = 1,numnod
        lay = layer(node)
        thetar(node) = cofgen(1,node)
        thetas(node) = cofgen(2,node)
!!! Kroes: disable combi of swsophy=1 and swhyst=1
        if (swhyst.eq.1) then
! ---     wetting curve  
          indeks(node) = 1
          cofgen(4,node) = paramvg(8,lay)
        elseif (swhyst.eq.0.or.swhyst.eq.2) then
! ---     drying branch or simulation without hysteresis
          indeks(node) = -1
          cofgen(4,node) = paramvg(4,lay)
        endif
      end do

! --- additional Input checks
!   - ThetCrMP should be less then CofGen(2,Lay): requires additional inputcheck
      if (FlMacropore) then
         do lay = 1, NumLay
            if (SwSoilShr(lay).gt.0 .and.                               &
     &          ThetCrMp(lay).gt.CofGen(2,nod1lay(lay))) then
               messag = ' ThetCrMP.gt.ThetSat'
               call fatalerr('MacroRead',messag)
            endif
         enddo
      endif

      if (swinco.eq.1) then
! --- pressure head profile is input 
        do i = 1, nhead
          tab(i*2) = h(i)
          tab(i*2-1) = abs(zi(i))
        end do
        do i = 1, numnod
          h(i) = afgen(tab,macp*2,abs(z(i)))
        end do
      endif
      if (swinco.eq.2 .and. swbotb.ne.8) then
        if (abs(gwli-(z(numnod)-0.5d0*dz(numnod))) .lt.1.0d-4) then
          messag = 'Initial groundwaterlevel (SWINCO=2) is '//          &
     &    'too close to bottom of soil profile'//                       &
     &    ' must be corrected!'
          call fatalerr ('soilwater',messag)
        endif
      endif
      if (swinco.eq.3) then
        if (nhead.ne.numnod) then
          messag = 'Initial data are read from file (SWINCO=3) and '//  &
     &    'number of nodes/compartments is not consistent with NUMNOD'//&
     &    'must be corrected!'
          call fatalerr ('soilwater',messag)
        endif
      endif
      if (swinco.eq.1.or.swinco.eq.3) then
! ---   determine groundwater level
        if (h(numnod) .gt. -1.d-5) then
          i = numnod
          do while ((h(i) .gt. -1.d-5) .and. (i .gt. 1))
              i = i - 1
          end do
          if (h(i) .lt. -1.d-5) then 
            gwl = z(i+1) + h(i+1) / (h(i+1) - h(i)) * (z(i) - z(i+1))
! ---       assume hydrostatic equilibrium in saturated part
            do j = i+1, numnod
              h(j) = gwl - z(j)
            end do
          endif
        endif
      else
! --- pressure head profile is calculated from groundwater level
        if (swbotb.eq.1) then  
          gwl = afgen (gwltab,mabbc*2,t1900+dt-1.d0)

          if(abs(gwl-(z(numnod)-0.5d0*dz(numnod))) .lt.1.0d-4) then
          messag = 'Groundwaterlevel as bottom boundary (SWBOTB=1) is'//&
     &    'below or to close to bottom of soil profile'//               &
     &    ' must be corrected!'
            call fatalerr ('soilwater',messag)
          endif
        else
          gwl = gwli
        endif
        if (gwl.gt.0.0d0) then 
          pond = gwl
        else
          pond = 0.0d0
        endif
        do i = 1,numnod
          h(i) = gwl - z(i)
        end do
      endif

! --- in case of preferential flow, adjust Van Genuchten parameters
      do i = 1, numnod
        theta(i) = watcon(i,h(i))
      end do

! --- hydr. conductivities, differential moisture capacities
!     and mean hydraulic conductivities for each node
      do node = 1,numnod
        dimoca(node) = moiscap(node,h(node))

        FrArMtrx(node) = 1.d0
        k(node) = hconduc (node,h(node),theta(node),rfcp(node))
        if(FlMacropore)  k(node) = FrArMtrx(node) * k(node)

        if(node.gt.1) kmean(node) =  hcomean(swkmean,k(node-1),k(node),dz(node-1),dz(node))
      end do
      kmean(numnod+1) = k(numnod)

! --- initial soil water storage
      if (.not.flMacroPore) then
         do i = 1, NumNod
            FrArMtrx(i) = 1.d0
         enddo
         volact = 0.0d0
         call watstor ()
         volini = volact
      endif
      pondini = pond 

! --- initial groundwater level
      call calcgwl ()

      return

      case (2)

! === calculate Soilwater rate/state variables ========================

! --- reset intermediate soil water fluxes
      if (flDayStart) then
          iqredwet_day = 0.0d0
          iqreddry_day = 0.0d0
          iqredsol_day = 0.0d0
          iqredfrs_day = 0.0d0
          iptra_day    = 0.0d0
          do node = 1,numnod 
            qpotrot_day(node) = 0.d0
            qredtot_day(node) = 0.d0
          enddo
      end if
      
      if (flzerointr) then
        do node = 1,numnod
          inqrot(node) = 0.0d0
          inqssdi(node) = 0.0d0
          inq(node) = 0.0d0
        enddo
        inq(numnod+1) = 0.0d0
        iqrot = 0.0d0
        iqssdi = 0.0d0
        iqredwet = 0.0d0
        iqreddry = 0.0d0
        iqredsol = 0.0d0
        iqredfrs = 0.0d0
        ies0 = 0.0d0
        iet0 = 0.0d0
        iew0 = 0.0d0
        iintc = 0.0d0
        iptra = 0.0d0
        ipeva = 0.0d0
        ievap = 0.0d0
        iruno = 0.0d0
        irunoCN = 0.0d0
        iqbot = 0.0d0
        iqtdo = 0.0d0
        iqtup = 0.0d0
        irunon = 0.0d0
        iqdo(1:numnod+1) = 0.0d0
        iqup(1:numnod+1) = 0.0d0

        IPondBeg = Pond
        do node = 1, numnod
          IThetaBeg(node) = Theta(node)
        enddo

!   - macropore variables
        if (flMacroPore) call macropore(5)
      endif

! --- reset cumulative soil water fluxes
      if (flzerocumu) then
        cqssdi = 0.0d0
        cqrot = 0.0d0
        cqbot = 0.0d0
        cqbotdo = 0.0d0
        cqbotup = 0.0d0
        cptra = 0.0d0
        cpeva = 0.0d0
        cevap = 0.0d0
        cinund = 0.0d0
        crunon = 0.0d0
        crunoff = 0.0d0
        crunoffCN = 0.0d0
        cqtdo = 0.0d0
        cqtup = 0.0d0
        cqprai = 0.0d0

!   - macropore variables
        if (flMacroPore) call macropore(6)

! --- reset initial water storage and ponding
        volini = volact
        pondini = pond
      endif

! --- save state variables of time = t
      call SoilWaterStateVar(1)
     
! --- calculate new soil water state variables
      call headcalc

      return

      case (3)

! --- update hydraulic conductivities to time level t+1
      do i = 1,numnod
         k(i) = hconduc(i,h(i),theta(i),rfcp(i))
         if(FlMacropore)  k(i) = FrArMtrx(i) * k(i)
         if(i.gt.1)then
            kmean(i) = hcomean(swkmean,k(i-1),k(i),dz(i-1),dz(i))
         end if
      enddo
      kmean(numnod+1) = k(numnod)

! --- calculate actual water content of profile
      call watstor ()

! --- calculate water fluxes between soil compartments
      call fluxes ()

! --- calculation of states macropores and intermediate & cumulative values
      if (flMacroPore) call macropore(4)

! --- calculate cumulative fluxes
      call integral 

! --- update parameters for soil water hystereses
      if (swhyst.ne.0) call hysteresis ()

      case default
         call fatalerr ('SoilWater', 'Illegal value for TASK')
      end select

      return
      end

      subroutine SoilWaterStateVar(task) 
! ----------------------------------------------------------------------
!     Date               : January 2007
!     Purpose            : save and reset soil water state variables
! ----------------------------------------------------------------------

! --- global variables
      use Variables
      implicit none

! --- local variables
      integer task, i

      select case (task)
      case (1)

! --- save state variables of time = t
      do i = 1,numnod
        hm1(i) = h(i)
        thetm1(i) = theta(i)
      enddo
      gwlm1 = gwl
      pondm1 = pond

      return

      case (2)
! --- reset soil state variables
      do i = 1,numnod
        h(i) = hm1(i)
        theta(i) = thetm1(i)
      enddo
      kmean(numnod+1) = k(numnod)
      gwl    = gwlm1
      pond   = pondm1

      case default
         call fatalerr ('SoilWaterStateVar', 'Illegal value for TASK')
      end select

      return
      end


