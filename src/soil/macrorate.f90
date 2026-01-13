! File VersionID:
!   $Id: macrorate.f90 370 2018-02-09 13:29:10Z heine003 $
! ----------------------------------------------------------------------
      SUBROUTINE MACRORATE(ITask,ICpBtDm,ICpBtPerZon,ICpSatGWl,         &
     &            ICpSatPeGWl,ICpTpPerZon,ICpTpSatZon,ICpTpWaSrDm,      &
     &            ArMpTpDm,AwlCorFac,FrMpWalWet,KDCrRlRef,SorpDmCp,     &
     &            ThtSrpRefDmCp,TimAbsCumDmCp,VlMpDm,VlMpDmCp,WaSrMp,   &
     &            WaSrMpDm,ZBtDm,ZWaLevDm,         flDraTub,FlEndSrpEvt,&
     &            QExcMtxDmCp,QInIntSatDmCp,QInMtxSatDmCp,QInTopLatDm,  &
     &            QInTopVrtDm,QOutDrRapCp,QOutMtxSatDmCp,QOutMtxUnsDmCp)
! ----------------------------------------------------------------------
!     Date               : april 2008                                       
!     Purpose            : 
!     Subroutines called : ABSORPTION, RAPIDDRAIN, SATFLOW                               
!     Functions called   : VOLUNDR
! ----------------------------------------------------------------------
      use Variables
      implicit NONE

! --- global                                                          In
      integer ICpBtDm(MaDm), ICpTpWaSrDm(MaDm), ITask
      real(8) ArMpTpDm(MaDm), AwlCorFac(Macp), FrMpWalWet(MaDm,MaCp) 
      real(8) KDCrRlRef(MaDr), SorpDmCp(MaDm,MaCp)
      real(8) ThtSrpRefDmCp(MaDm,MaCp), TimAbsCumDmCp(MaDm,MaCp)
      real(8) VlMpDm(MaDm), VlMpDmCp(MaDm,MaCp), WaSrMp, WaSrMpDm(MaDm)
      real(8) ZBtDm(MaDm), ZWaLevDm(MaDm)
      logical flDraTub(Madr)
!     -                                                              Out
      real(8) QExcMtxDmCp(MaDm,MaCp), QInIntSatDmCp(MaDm,MaCp)
      real(8) QInMtxSatDmCp(MaDm,MaCp), QInTopLatDm(MaDm)
      real(8) QInTopVrtDm(MaDm), QOutDrRapCp(MaCp) 
      real(8) QOutMtxSatDmCp(MaDm,MaCp), QOutMtxUnsDmCp(MaDm,MaCp) 
      logical FlEndSrpEvt(MaDm,MaCp)
! ----------------------------------------------------------------------
! --- local
      integer ic, ICpBtPerZon, ICpSatGWl, ICpSatPeGWl, ICpTpPerZon  
      integer ICpTpSatZon, id, jd, Sq(MaDM), SqHlp
      real(8) CritHtop, FrFlwIn, FrFlwOut, FlwInDmExces 
      real(8) FlwInDmTot, FlwInIntSatDmCpPot(MaCp), FlwInIntSatDmPot
      real(8) FlwInMtxSatDmCpPot(MaCp), FlwInMtxSatDmPot
      real(8) FlwInTopExcesTot, FlwInTopLatDmPot(MaDm), FlwInTopPot
      real(8) FlwInTopXtr, FlwInTopVrtDmPot(MaDm), FlwOutDmExces(MaDm) 
      real(8) FlwOutDrRapCpPot(MaCp), FlwOutDrRapPot, FlwOutDmTot(MaDm)
      real(8) FlwOutMtxSatDmCpPot(MaDm,MaCp), FlwOutMtxSatDmPot(MaDm)
      real(8) FlwOutMtxUnsDmCpPot(MaDm,MaCp), FlwOutMtxUnsDmPot(MaDm)
      real(8) FrReduQ, FrXtr, Pp, PpTot, SatDefDm(MaDm)
      real(8) SatDefDmRel(MaDm), SatDefTot, VlMpUndrDrL, VlMpUndrGwL
      real(8) VOLUNDR, WaSrMpDmMax, WaSrMpDmMin, WaSrMpDmTmp
      real(8) Ld, pi, r0, w_geom, Henpr1  

      data    pi /3.14159d0/

! ----------------------------------------------------------------------
!   for FORCHECK
      t = t
      t1900 = t1900
!
!> I << ---------- INITIALISE AND CHECK MACROPORE STATUS ---------------

! --- Initialisation Fluxes and Derivatives
      QMaPo= 0.d0
      QRapDra= 0.d0
      do 10 id= 1,NumDm
         QInTopVrtDm(id)= 0.d0
         QInTopLatDm(id)= 0.d0
         do 9 ic= 1, NumNod
            QExcMtxDmCp(id,ic)   = 0.d0 
            QInIntSatDmCp(id,ic) = 0.d0
            QInMtxSatDmCp(id,ic) = 0.d0
            QOutMtxSatDmCp(id,ic)= 0.d0
            QOutMtxUnsDmCp(id,ic)= 0.d0
            if (id.eq.1) then
               QOutDrRapCp(ic) = 0.d0
               QExcMpMtx(ic)   = 0.d0
               dFdhMp(ic)      = 0.d0  
            endif
   9     continue
  10  continue

! --- Check whether the macropore status requires the working of this subroutine
      if (IcTopMp.eq.1) then                                        ! no covering layer on top of macropores
         FlwInTopPot= QMpLatSs + ArMpTp * (NRaiDt+NIrd+Melt) * DT
      else                                                          ! covering layer on top of macropores
         FlwInTopPot= 0.d0
      endif
      if (VlMp.lt.1.d-6 .or. (WaSrMp.lt.1.d-6 .and. FlwInTopPot.lt.1.d-6     &
     &    .and. ICpBtDm(1).le.NodGwL .and. NPeGwL.lt.IcTopMp) ) then ! not required
         do 20 id= 1, NumDm
            do 19 ic= IcTopMp, NumNod
               FlEndSrpEvt(id,ic)= .true.
  19        continue
  20     continue
         QMpLatSs= 0.d0
         return
      endif

!!!      FrReduQ= 0.1d0**0.d0 !(float(0.))   !IDecMpRat))     DIT HEEFT INVLOED: LANGERE RUN TIJD!!
      FrReduQ= 0.1d0**(float(IDecMpRat))
!      if (IDecMpRat.gt.0) write(94,*) t,IDecMpRat     
      

!> II < ----- CALCULATE RATES AND DERIVATIVES PER MACROPORE DOMAIN -----
!
! --- Subroutine is required; calculations for all macropore domains
!   - General initialisation
      FlwInTopExcesTot= 0.d0
      SatDefTot= 0.d0
!
      do 100 id= 1, NumDm
         FlwInDmTot= 0.d0  
         FlwOutDmTot(id)= 0.d0

!- A. CALCULATE POTENTIAL IN- AND OUTGOING WATER FLOWS

! - 1. CASE 0:   INITIALIZATION
         call SATFLOW(0,ICpBtDm(id),ICpBtPerZon,ICpSatPeGWl,            &
     &          ICpTpWaSrDm(id),ICpTpPerZon,id,FrMpWalWet,FrReduQ,Henpr1,&
     &          PeGWL,QInMtxSatDmCp,QOutMtxSatDmCp,ZWaLevDm(id),        &
     &          FlwInIntSatDmCpPot,FlwInIntSatDmPot,FlwOutMtxSatDmCpPot,&
     &          FlwOutMtxSatDmPot(id))

! - 1. CASE 1:   INCOMING WATER FLOWS
!   1.a. INFLoW into macropores at soil surface (FlwInTop...) by infiltration:
!        directly by vertical inflow (Vrt; precipitation Pre or seepage covering layer) 
!        and indirectly by lateral (Lat) overland flow of matrix infiltration excess (runoff)
         FlwInTopLatDmPot(id)= 0.d0

         if (IcTopMp.eq.1) then                                        ! no covering layer on top of macropores
            FlwInTopVrtDmPot(id)= ArMpTpDm(id) * (NRaiDt+NIrd+Melt) * DT
            if (ArMpTp.gt.0.d0) FlwInTopLatDmPot(id)= ArMpTpDm(id)/ArMpTp * QMpLatSs 
         else                                                          ! covering layer on top of macropores
            ic = IcTopMp-1
            Henpr1 = 0.d0 ! H_enpr(1) 
            CritHtop = Henpr1
            if (h(ic).gt.CritHtop) then  
               !Lay= Layer(ic)
               Ld   = dipomi
               r0   = 0.5d0 * dipomi * (1.d0 - dsqrt(1.d0-VlMpStCp(IcTopMp)-VlMpDyCp(IcTopMp)))              
               w_geom = 1.d0 / (1.d0 + Ld/(pi*dz(ic)/2.d0)*log((Ld/(pi*r0))) + Ld**2/(6.d0*dz(ic)**2))
!
               FlwInTopVrtDmPot(id)= w_geom * KsatCovLay * ((h(ic)-Henpr1) /(dz(ic)/2.d0) + 1.d0) * DT
               FlwInTopVrtDmPot(id)= ArMpTpDm(id)/ArMpTp * max(0.d0,FlwInTopVrtDmPot(id)) 
            else
               FlwInTopVrtDmPot(id)= 0.d0
            endif
         endif
!
!   1.b. INFLoW into macropores from SATurated top compartments by INTerflow: 
         call SATFLOW(1,ICpBtDm(id),ICpBtPerZon,ICpSatPeGWl,            &
     &          ICpTpWaSrDm(id),ICpTpPerZon,id,FrMpWalWet,FrReduQ,Henpr1,&
     &          PeGWL,QInMtxSatDmCp,QOutMtxSatDmCp,ZWaLevDm(id),        &
     &          FlwInIntSatDmCpPot,FlwInIntSatDmPot,FlwOutMtxSatDmCpPot,&
     &          FlwOutMtxSatDmPot(id))

!   1.c. INFLoW into macropores from SATurated MaTriX compartm. (exfiltration): 
         call SATFLOW(1,ICpBtDm(id),NumNod,ICpSatGWl,ICpTpWaSrDm(id),   &
     &          ICpTpSatZon,id,FrMpWalWet,FrReduQ,Henpr1,GWlFlCpZo,     &
     &          QInMtxSatDmCp,QOutMtxSatDmCp,ZWaLevDm(id),              &
     &          FlwInMtxSatDmCpPot,FlwInMtxSatDmPot,FlwOutMtxSatDmCpPot,&
     &          FlwOutMtxSatDmPot(id))

! - 2. CASE 2:   OUTGOING WATER FLOWS 
!   2.a. OUTFLoW out off macrop. into SATurated MaTriX compartm. (infiltrat.): 
         call SATFLOW(2,ICpBtDm(id),NumNod,ICpSatGWl,ICpTpWaSrDm(id),   &
     &          ICpTpPerZon,id,FrMpWalWet,FrReduQ,Henpr1,GWlFlCpZo,     &
     &          QInMtxSatDmCp,QOutMtxSatDmCp,ZWaLevDm(id),              &
     &          FlwInMtxSatDmCpPot,FlwInMtxSatDmPot,FlwOutMtxSatDmCpPot,&
     &          FlwOutMtxSatDmPot(id))

!   2.b. OUTFLoW out off macrop. into UNSaturated MaTriX compart. (absorption): 
         call ABSORPTION(1,ICpBtDm(id),ICpBtPerZon,ICpTpPerZon,         &
     &          ICpTpSatZon,ICpTpWaSrDm(id),id,AwlCorFac,FrMpWalWet,    &
     &          FrReduQ,QOutMtxUnsDmCp,SorpDmCp,ThtSrpRefDmCp,          &
     &          TimAbsCumDmCp,FlwOutMtxUnsDmCpPot,FlwOutMtxUnsDmPot(id),&
     &          FlEndSrpEvt,ZWaLevDm)

!   2.c. OUTFLoW out off macropores by RAPid DRainage. Only Main Domain (id= 1): 
         call RAPIDDRAIN(id,ICpBtDm(1),ICpTpWaSrDm(1),flDraTub,         &
     &          FrMpWalWet,FrReduQ,KDCrRlRef,VlMpDmCp,WaSrMpDm(1),      &
     &          ZBtDm(1),ZWaLevDm(1),FlwOutDrRapCpPot,FlwOutDrRapPot,   &
     &          VlMpUndrDrL)


!- B. CALCULATE NEW WATER STORAGE AND CHECK WITH AVAILABLE DOMAIN VOLUME

         FlwInDmTot= FlwInTopVrtDmPot(id) + FlwInTopLatDmPot(id) +      &
     &               FlwInIntSatDmPot + FlwInMtxSatDmPot
         FlwOutDmTot(id)= FlwOutMtxSatDmPot(id) + FlwOutMtxUnsDmPot(id)
         if (id.eq.1) FlwOutDmTot(1)= FlwOutDmTot(1) + FlwOutDrRapPot
!
!   - temporary storage in domain
         WaSrMpDmTmp= WaSrMpDm(id) + FlwInDmTot - FlwOutDmTot(id)
! 
!   - volume of domain below groundwater level in matrix
         VlMpUndrGwL= 0.d0
         if (ZBtDm(id).lt.GWlFlCpZo .and. GWlFlCpZo.lt.900.d0)          &  ! GWL = 999. indicates NO groundwater level
     &       VlMpUndrGwL= VOLUNDR(id,ICpBtDm(id),GWlFlCpZo,VlMpDmCp,    &
     &                            ZBtDm(id))
         VlMpUndrGwL= dmin1(VlMpUndrGwL,VlMpDm(id))
         
!   - if inflow is only from saturated matrix in real (not perched) groundwater 
!     then level in domain can not exceed ground waterlevel         
         WaSrMpDmMax= VlMpDm(id)      
         if((FlwInDmTot-FlwInMtxSatDmPot)/DT.lt.1.d-7 .and.             &
     &       FlwInMtxSatDmPot/DT.gt.1.d-7) WaSrMpDmMax= VlMpUndrGwL     
!
!   1. In case of inflow excess, decrease incoming fluxes
         FrFlwIn= 1.d0
         if (WaSrMpDmTmp.gt.WaSrMpDmMax+1.d-7) then
            FlwInDmExces= WaSrMpDmTmp - WaSrMpDmMax
            FrFlwIn= dmax1(0.d0,1.d0 - FlwInDmExces/FlwInDmTot)
!    - total inflow excess at top for decreasing top inflow in main SWAP
            FlwInTopExcesTot= FlwInTopExcesTot + (1.d0 - FrFlwIn) *     &
     &                       (FlwInTopVrtDmPot(id)+FlwInTopLatDmPot(id))

         if (frflwin.lt.0.1d0) then
         continue
         endif

         endif
!    - decrease inflow tems of domain
         QInTopVrtDm(id)= FrFlwIn * FlwInTopVrtDmPot(id) / DT
         QInTopLatDm(id)= FrFlwIn * FlwInTopLatDmPot(id) / DT
         do 90 ic= IcTopMp, NumNod
            QInIntSatDmCp(id,ic)= FrFlwIn * FlwInIntSatDmCpPot(ic) / DT
            QInMtxSatDmCp(id,ic)= FrFlwIn * FlwInMtxSatDmCpPot(ic) / DT
  90     continue
!
!   - determine minimum water storage in domain
         WaSrMpDmMin= dmin1(VlMpUndrGwL,WaSrMpDm(id))                   
!   - in case of MB domain drainage is possible                      
         if (id.eq.1 .and. VlMpUndrDrL.gt.1.d-6 .and.                   &
     &       VlMpUndrDrL.lt.WaSrMpDmMin) WaSrMpDmMin= VlMpUndrDrL
!
         FlwOutDmExces(id)= 0.d0
         if (WaSrMpDmTmp.lt.WaSrMpDmMin) FlwOutDmExces(id)=             &
     &                                   WaSrMpDmMin - WaSrMpDmTmp

!      saturation deficit of domain and of total of domains: for checkingavailable storage for excess from other domains
         SatDefDm(id)= dmax1(0.d0,(VlMpDm(id)-WaSrMpDmTmp))
         SatDefDmRel(id)= SatDefDm(id) / PpDmCp(id,IcTopMp)
         SatDefTot= SatDefTot + SatDefDm(id)
 100  continue
!
!-  End of loop domains
!

!    - Redistribute excess over available volume of domains that are not filled up yet
      if (FlwInTopExcesTot.gt.1.d-6 .and. SatDefTot.gt.1.d-6) then        ! action is needed and additional storage is available
!      put relative saturation deficit SatDefDmRel(id) in ascending order
         do 200 id= 1, NumDm
            Sq(id)= id
 200     continue
         do 210 id= 1, NumDm-1
            do 209 jd= id+1, NumDm
               if (SatDefDmRel(Sq(id)).gt.SatDefDmRel(Sq(jd))) then
                  SqHlp= Sq(id)
                  Sq(id)= SQ(jd)
                  Sq(jd)= SqHlp
               endif
 209        continue
 210     continue

!    - Redistribute excess 
         PpTot= 1.d0
         do 220 id= 1, NumDm
            jd= Sq(id)
            if (SatDefDm(jd).gt.1.d-7 .and.                             &
     &          FlwInTopVrtDmPot(jd)+FlwInTopLatDmPot(jd).gt.1.d-7) then
               Pp= PpDmCp(jd,IcTopMp) / PpTot
               FlwInTopXtr= dmin1(Pp*FlwInTopExcesTot,SatDefDm(jd))
               FrXtr= FlwInTopXtr /                                     &
     &                (FlwInTopVrtDmPot(jd) + FlwInTopLatDmPot(jd))
               QInTopVrtDm(jd)= (1.d0+FrXtr) * FlwInTopVrtDmPot(jd) / DT
               QInTopLatDm(jd)= (1.d0+FrXtr) * FlwInTopLatDmPot(jd) / DT
!
               FlwInTopExcesTot = FlwInTopExcesTot - FlwInTopXtr
               FlwOutDmExces(jd)= FlwOutDmExces(jd) - FlwInTopXtr
               FlwOutDmExces(jd)= dmax1(0.d0,FlwOutDmExces(jd))
               SatDefDm(jd)= SatDefDm(jd) - FlwInTopXtr
            else
               SatDefDm(jd)= 0.d0
            endif
            PpTot= PpTot - PpDmCp(jd,IcTopMp)
 220     continue
      endif

!   2. Calculate actual outgoing fluxes on basis of (possibly) outflow excess
      do 230 id= 1, NumDm
         FrFlwOut= 1.d0
         if (FlwOutDmExces(id)/DT.gt.1.d-7) then
            FrFlwOut= 1.d0 - FlwOutDmExces(id)/FlwOutDmTot(id)
            FrFlwOut= dmax1(0.d0,FrFlwOut)
         endif 
         do 228 ic= IcTopMp, NumNod
            QOutMtxSatDmCp(id,ic)= FrFlwOut *                           &
     &                             FlwOutMtxSatDmCpPot(id,ic) / DT
            QOutMtxUnsDmCp(id,ic)= FrFlwOut *                           &
     &                             FlwOutMtxUnsDmCpPot(id,ic) / DT
 228     continue
         if (id.eq.1) then
            do 229 ic= IcTopMp, NumNod
               QOutDrRapCp(ic)= FrFlwOut * FlwOutDrRapCpPot(ic) / DT
 229        continue
         endif
 230  continue

!   3. Remaining inflow excess is substracted from lateral inflow into macropores 
!      for use in main SWAP
      if (IcTopMp.eq.1) then
      QMpLatSs= QMpLatSs - FlwInTopExcesTot 
      if (dabs(QMpLatSs).lt.1.0d-7) QMpLatSs= 0.0d0
      else

      if(FlwInTopExcesTot.gt.0.d0) then
      continue
      endif

         do id = 1, numdm
!            QInTopVrtDm(id)= QInTopVrtDm(id) - ArMpTpDm(id)/ArMpTp * FlwInTopExcesTot  / DT 
!            QInTopVrtDm(id)= max(QInTopVrtDm(id),0.d0) 
            if (QInTopVrtDm(id).lt.0.d0) then
             continue
            endif
         enddo
      endif
!
!- C. CALCULATE FLUXES FOR USE IN SWAP
!
!   1. Aggregate fluxes for output
      do 240 id= 1, NumDm
         do 239 ic= IcTopMp, ICpBtDm(id)
            QExcMtxDmCp(id,ic)=                                         &
     &                   QOutMtxSatDmCp(id,ic) + QOutMtxUnsDmCp(id,ic) -&
     &                   (QInIntSatDmCp(id,ic) + QInMtxSatDmCp(id,ic))
 239     continue
 !
         if (IcTopMp.gt.1) then
            ic = IcTopMp - 1
            QInMtxSatDmCp(id,ic)= QInTopVrtDm(id)
            QExcMtxDmCp(id,ic)   = -QInMtxSatDmCp(id,ic)
         endif
 240  continue

!   2. Aggregate fluxes for main SWAP
      do 250 ic= 1, NumNod
         do 249 id= 1, NumDm
            QExcMpMtx(ic)   = QExcMpMtx(ic) + QExcMtxDmCp(id,ic)
 249     continue
         QMapo= QMapo + QExcMpMtx(ic)
         QRapDra= QRapDra + QOutDrRapCp(ic)
 250  continue
!
      if (ITask.eq.1) return
!              ------------- END of ITask 1 -------------
!
!- D. CALCULATE DERIVATIVES FOR USE IN SWAP (Headcalc)
!
!   1. exchange with unsaturated matrix
      do 260 id= 1, NumDm
         call ABSORPTION(2,ICpBtDm(id),ICpBtPerZon,ICpTpPerZon,         &
     &          ICpTpSatZon,ICpTpWaSrDm(id),id,AwlCorFac,FrMpWalWet,    &
     &          FrReduQ,QOutMtxUnsDmCp,SorpDmCp,ThtSrpRefDmCp,          &
     &          TimAbsCumDmCp,FlwOutMtxUnsDmCpPot,FlwOutMtxUnsDmPot(id),&
     &          FlEndSrpEvt,ZWaLevDm)
 260  continue
!
!   2. exchange with saturated matrix
      do 270 id= 1, NumDm
         call SATFLOW(4,ICpBtDm(id),NumNod,ICpSatGWl,ICpTpWaSrDm(id),   &
     &          ICpTpSatZon,id,FrMpWalWet,FrReduQ,Henpr1,GWL,           &
     &          QInMtxSatDmCp,                                          &
     &          QOutMtxSatDmCp,ZWaLevDm(id),FlwInMtxSatDmCpPot,         &
     &          FlwInMtxSatDmPot,FlwOutMtxSatDmCpPot,                   &
     &          FlwOutMtxSatDmPot(id))

 270  continue

      return
      END


!=====================================================================
      SUBROUTINE SATFLOW(ITask,ICpBtDm,ICpBtZon,ICpSatLev,ICpTpWaSrDm,  &
     &              ICpTpZon,id,FrMpWalWet,FrReduQ,Henpr1,Lev,          &
     &              QInMtxSatDmCp,QOutMtxSatDmCp,ZWaLevDm,              &
     &              FlwInDmCpPot,FlwInDmPot,FlwOutDmCpPot,FlwOutDmPot)
! ----------------------------------------------------------------------
!     Date               : April 2008                                        
!     Purpose            : To calculate water exchange between macropores
!                        : and the SATURATED soil matrix
!     Functions called   : -
! ----------------------------------------------------------------------
      use Variables
      implicit NONE

! --- global                                                          In
      integer id, ICpBtDm, ICpBtZon,ICpSatLev,ICpTpWaSrDm,ICpTpZon,ITask
      real(8) FrMpWalWet(MaDm,MaCp), FrReduQ, Henpr1, Lev, ZWaLevDm 
      real(8) QInMtxSatDmCp(MaDm,MaCp), QOutMtxSatDmCp(MaDm,MaCp)
!     -                                                              Out
      real(8) FlwInDmCpPot(MaCp), FlwInDmPot
      real(8) FlwOutDmCpPot(MaDm,MaCp), FlwOutDmPot
! ----------------------------------------------------------------------
! --- local
      integer ic, Lay
      real(8) DelH, DelHDmCp(MaDm,MaCp), FlwMtxSatDmCp(MaCp), HMp, HMa
      real(8) ksatHor, Pi, RatSeepFace, RecRes, ResHor, ResRad,ResVrt
      real(8) CritHtop

      data Pi, RatSeepFace / 3.14159d0, 10.d0 /
! ----------------------------------------------------------------------
      select case (itask)

! - 0. INITIALISATION  
      case (0)
! --- Initialisation global Flows and contribution Derivatives
      do 10 ic= 1, NumNod
         FlwMtxSatDmCp(ic)= 0.d0
         DelHDmCp(id,ic)  = 0.d0
  10  continue

! --- Calculate potential water EXCHANGE between macropores and saturated matrix
!
! - 1. Calculate potential flows INCOMING into macropores   
!
      case(1)
! --- Initialisation local Flows 
      FlwInDmPot= 0.d0
      do 15 ic= 1, NumNod
         FlwInDmCpPot(ic) = 0.d0
  15  continue

! --- Check whether this subroutine is relevant for this macropore domain
      if (ICpBtDm.ge.ICpTpZon .and. ICpBtZon.gt.0) then
         continue
      else
         return
      endif
!
!   - ICpTpZon = top compartment of relevant saturated zone
      do 100 ic= ICpTpZon, min0(ICpBtDm,ICpBtZon)
         Lay= Layer(ic)
         ksatHor=  ksatfit(Lay)   
!
!   - First calculation of Heads in Macropores and Matrix, for determing flow direction
         HMp = ZWaLevDm - Z(ic)  ! head in macropore domain, assuming static equilibrium
         HMp = dmax1(0.d0,HMp)   ! no negative head in macropores
         if (HMp.lt.1.d-8) HMp = 0.d0
         HMa = H(ic)             ! absolute comparison, no influence H_enpr
         DelH= HMp - HMa         ! differences between heads determines flow direction
         if (HMp.lt.1.d-8 .and. DelH.gt.0.d0) DelH= 0.d0
         if (abs(DelH).lt.1.d-8) DelH = 0.d0

! - 3 flow situations:
!    1. DelH = +  -> infiltration of macropore water into matrix; Darcy flow
!    2. and 3. DelH = -  -> exfiltration of matrix water into macropores  
!    2.        HMp > 0   -> Darcy flow
!    3         HMp = 0   -> seepage face; Hooghoudt above drain level
         if (Delh.gt.0.d0) then 
! - 1. infiltration out of macropores into matrix: = term 2.a. (Darcy)
!           RecRes = reciprocal resistance for calculation of derivative
            RecRes = ShapeFacMp * 8.0d0 * PpDmCp(id,ic) * DZ(ic) *      &
     &               ksatHor / DiPoCp(ic)**2
            if (ic.eq.ICpTpWaSrDm) RecRes = FrMpWalWet(id,ic) * RecRes
            FlwMtxSatDmCp(ic) = - FrReduQ * RecRes * DelH * DT
            DelHDmCp(id,ic)= DelH  ! save head difference for calculation of derivative
!
         elseif (DelH.lt.0.d0) then 
!   - exfiltration out of matrix into macropores
            if (HMp.gt.0.d0) then
! - 2. flow equation according to Darcy 
               RecRes = ShapeFacMp * 8.0d0 * PpDmCp(id,ic) * DZ(ic) *   &
     &                                       ksatHor / DiPoCp(ic)**2
               if (ic.eq.ICpSatLev) RecRes = RecRes *                   &
     &                       (GWlFlCpZo - (Z(ic)-0.5d0*DZ(ic))) / DZ(ic)
               DelHDmCp(id,ic)= DelH ! save head difference for calculation of derivative
            else
! - 3. flow equation according to seepage face: Ernst with radial resistance 
!               RecRes = -4.0d0 * PpDmCp(id,ic) * ksatHor * DelH / DiPoCp(ic)**2 
!               DelHDmCp(id,ic)= DelH / 2.d0 ! div. by 2 because of DelH**2 
               ResHor = DiPoCp(ic)**2 / (8.0d0 * DZ(ic) * ksatHor) 
               ResVrt = DZ(ic) / ksatHor
               ResRad = DiPoCp(ic) * DLog(RatSeepFace) / (Pi*ksatHor)
               RecRes = PpDmCp(id,ic) / (ResHor + ResVrt + ResRad) 
               if (ic.eq.ICpSatLev) RecRes = RecRes *                   &
     &                             (Lev - (Z(ic)-0.5d0*DZ(ic))) / DZ(ic)
               DelHDmCp(id,ic)= DelH ! save head difference for calculation of derivative
            endif
            FlwMtxSatDmCp(ic) = - FrReduQ * RecRes * DelH * DT
         endif
!   - save positive FlwMtxSatDmCp as INCOMING (into macropores) flux
         if (FlwMtxSatDmCp(ic).gt.0.d0) then
            FlwInDmCpPot(ic) = FlwMtxSatDmCp(ic)
            FlwInDmPot       = FlwInDmPot + FlwInDmCpPot(ic)
         endif
 100  continue 

      return

      case (2)
! - 2. Calculate potential flows OUTGOING out off macropores   
! --- Initialisation Flows
      FlwOutDmPot = 0.d0
      do 200 ic= 1, NumNod
         FlwOutDmCpPot(id,ic)= 0.d0
 200  continue

! --- Check whether this subroutine is relevant for this macropore domain
      if (ICpBtDm.ge.ICpTpZon) then
         continue
      else
         return
      endif

!   - save negative FlwMtxSatDmCp as OUTGOING (out off macropores) flux
!   - (calculation of FlwMtxSatDmCp(ic) done above, in ITask 1)
      do 250 ic= ICpTpZon, ICpBtDm
         if (FlwMtxSatDmCp(ic).lt.0.d0) then
            FlwOutDmCpPot(id,ic)= - FlwMtxSatDmCp(ic)
            FlwOutDmPot         = FlwOutDmPot + FlwOutDmCpPot(id,ic)
         endif
 250  continue 

      return

      case (4)
! - C. Determine derivatives 

      do 300 ic= ICpTpZon, ICpBtDm
         if (abs(DelHDmCp(id,ic)).gt.1.d-14) then
            dFdhMp(ic)= dFdhMp(ic) -                                    &
     &                  (QOutMtxSatDmCp(id,ic) - QInMtxSatDmCp(id,ic)) /& 
     &                  DelHDmCp(id,ic) 
         end if
 300  continue

!   - derivative for bottom compartiment of covering top layer in case of seepage out of this layer into macropores 
      if (IcTopMp.gt.1)  then
         ic = IcTopMp-1
         CritHtop = 0.d0   !        Henpr1                      ! Adapted for GEM: test for H_enpr
         if (h(ic).gt.CritHtop) then  
            DelHDmCp(id,ic) = 0.d0 - (h(ic)-Henpr1)  ! Adapted for GEM: test for H_enpr
            dFdhMp(ic)= dFdhMp(ic) + QInMtxSatDmCp(id,ic) / DelHDmCp(id,IcTopMp-1)   
         endif        
      endif
!
      case default
         call fatalerr ('SatFlow', 'Illegal value for TASK')
      end select

      Return
      End


!=====================================================================
      SUBROUTINE ABSORPTION(ITask,ICpBtDm,ICpBtPerZon,ICpTpPerZon,      &
     &              ICpTpSatZon,ICpTpWaSrDm,id,AwlCorFac,FrMpWalWet,    &
     &              FrReduQ,QOutMtxUnsDmCp,SorpDmCp,ThtSrpRefDmCp,      &
     &              TimAbsCumDmCp,FlwOutMtxUnsDmCpPot,FlwOutMtxUnsDmPot,&
     &              FlEndSrpEvt,ZWaLevDm)

! ----------------------------------------------------------------------
!     Date               : May 2005                                        
!     Purpose            : To calculate absorption of macropore water into the
!                          unsaturated matrix by sorptivity under dry conditions
!                          and by pressure head gradient with Darcy
!     Functions called   : moiscap
! ----------------------------------------------------------------------
      use Variables
      implicit NONE

! --- global                                                          In
      integer id, ICpBtDm, ICpBtPerZon, ICpTpPerZon, ICpTpSatZon
      integer ICpTpWaSrDm, ITask
      real(8) AwlCorFac(Macp), FrMpWalWet(MaDm,MaCp) 
      real(8) QOutMtxUnsDmCp(MaDm,MaCp), SorpDmCp(MaDm,MaCp)
      real(8) ThtSrpRefDmCp(MaDm,MaCp), TimAbsCumDmCp(MaDm,MaCp)
      real(8) ZWaLevDm(MaDm), FrReduQ
!     -                                                              Out
      real(8) FlwOutMtxUnsDmCpPot(MaDm,MaCp), FlwOutMtxUnsDmPot
      logical FlEndSrpEvt(MaDm,MaCp)
! ----------------------------------------------------------------------
! --- local
      integer ic, ICpBtUnsMtxDm, Lay
      real(8) AbsDarc, AbsSorp, CritSatDefMtx, DelH, DelHDmCp(MaDm,MaCp)
      real(8) Deriv, HMp, moiscap, QSorpMax, RecRes, SatDefMtx
      real(8) SorpAct, SorpFac, Time
      logical FlSorp(MaDm,MaCp)

!      integer tel
! ----------------------------------------------------------------------
! --- Initialization Flows and Flags
      FlwOutMtxUnsDmPot= 0.d0
      do 10 ic= 1, NumNod
         FlwOutMtxUnsDmCpPot(id,ic)= 0.d0
         FlEndSrpEvt(id,ic)        = .true.
  10  continue

! --- Initialization help variables
      ICpBtUnsMtxDm= min0(ICpBtDm,ICpTpSatZon-1) 
      CritSatDefMtx= 1.0d-8
      QSorpMax     = 1.0d3                                     
!
      select case (itask)
      case (1)
! --- Calculate potential Absorption
      do 100 ic= ICpTpWaSrDm, ICpBtUnsMtxDm
         if (ic.gt.ICpBtPerZon .or. ic.lt.ICpTpPerZon) then 
!
            Lay      = Layer(ic)
! - I By Sorptivity           
            Time     = TimAbsCumDmCp(id,ic)
            SatDefMtx= ThetaS(ic) - Theta(ic)
            if (SatDefMtx.lt.CritSatDefMtx) then
!    -  Matrix too wet: no sorptivity
               AbsSorp = 0.d0
            else
               if (Time.lt.1.d-8) then  
!    -  Start a new sorptivity event; set ThtSrpRefDmCp(id,ic) and SorpDmCp
!    - ThtSrpRefDmCp = Thet_sat + Delta_thet_theor. = Thet_sat + Thet_theor. - Thet_0
!      Delta_thet_theor. = theoretical increase of theta (calculated in MacroState) 
                  ThtSrpRefDmCp(id,ic)= ThetaS(ic)  ! initial, delta_theta = 0
                  SorpDmCp(id,ic)= SorpMax(Lay) *                       &
     &                          ((ThetaS(ic)-Theta(ic)) /               &
     &                          (ThetaS(ic)-ThetaR(ic)))** SorpAlfa(Lay)
                  SorpAct = SorpDmCp(id,ic)
               elseif((ThtSrpRefDmCp(id,ic)-Theta(ic)).gt.CritSatDefMtx) then 
!    -  Continue existing sorptivity event; set SorpAct
                  SorpAct= SorpMax(Lay) *                               &
     &                     ((ThtSrpRefDmCp(id,ic)-Theta(ic)) /          &
     &                     (ThetaS(ic)-ThetaR(ic)))** SorpAlfa(Lay)
               else
!    -  End existing sorptivity event
                  SorpAct= 0.d0
               endif
               AbsSorp = SorpAct * PpDmCp(id,ic) * (4.d0 *AwlCorFac(ic)*&
     &               DZ(ic) / DiPoCp(ic)) * (dsqrt(Time+DT)-dsqrt(Time))   
               if (ic.eq.ICpTpWaSrDm) AbsSorp= FrMpWalWet(id,ic)*AbsSorp
! --- Limit peak amount of absorbed water at beginning of sorptiviy event
               AbsSorp = dmin1(QSorpMax*DT*DZ(ic),AbsSorp)
            endif

!   II By Darcy 

            HMp = ZWaLevDm(id) - Z(ic)  ! head in macropore domain, assuming static equilibrium
            HMp = dmax1(0.d0,HMp)       ! no negative head in macropores
            if (H(ic).lt.H_enpr(lay)-1.d-8 .and. HMp.gt.1.d-8) then     ! Adapted for GEM: test for H_enpr        
               DelH= dmax1(0.d0,HMp - H(ic))
            else
               DelH= 0.d0
            endif
            DelHDmCp(id,ic)= DelH     ! save head difference for calculation of derivative
                                      ! reciporocal resistance for calculation of derivative
            RecRes  = ShapeFacMp * 8.0d0 * PpDmCp(id,ic) * DZ(ic) *     &
     &                K(ic) / DiPoCp(ic)**2
            if (ic.eq.ICpTpWaSrDm) RecRes = FrMpWalWet(id,ic) * RecRes
            AbsDarc = RecRes * DelH * DT
!
            SorpFac = SorpFacParl(lay) + (1.d0-SorpFacParl(lay)) *      &
     &                ( 1.d0 -(dmax1((ThetaS(ic)-Theta(ic)),0.d0) /     &
     &                         (ThetaS(ic)-ThetaR(ic)))**SorpAlfa(Lay) )

            if (SwDarcy.eq.0) AbsDarc = 0.d0

            if (AbsSorp.gt.SorpFac*AbsDarc) then
               FlwOutMtxUnsDmCpPot(id,ic) = FrReduQ * AbsSorp
               FlSorp(id,ic) = .true.
               AbsDarc = 0.d0
!            tel = tel + 1
!            write(99,99) t, absdarc/dt, abssorp/dt, theta(ic), time,ic,tel
            else
!            tel = tel + 1
!            write(99,99) t, absdarc/dt, abssorp/dt, theta(ic), time,ic,tel
               FlwOutMtxUnsDmCpPot(id,ic) = FrReduQ * SorpFac*AbsDarc
!            FlwOutMtxUnsDmCpPot(id,ic) = AbsDarc
               AbsSorp = 0.d0
               FlSorp(id,ic) = .false.
            endif
!    
            FlwOutMtxUnsDmPot= FlwOutMtxUnsDmPot +                      &
     &                         FlwOutMtxUnsDmCpPot(id,ic)
!
! --- Set Flag to END SoRPtivity EVenT in MACROSTATE
            if (AbsSorp/DT.gt.1.d-7) then
               FlEndSrpEvt(id,ic)= .false.
            endif
         endif
!
 100  continue
!
      return
!
      case (2)
! --- Determine derivatives 
      do ic= ICpTpWaSrDm, ICpBtUnsMtxDm
         if (ic.gt.ICpBtPerZon .or. ic.lt.ICpTpPerZon) then 
!
            if (QOutMtxUnsDmCp(id,ic).gt.1.d-7) then
               if (FlSorp(id,ic)) then
                  Lay   = Layer(ic)
                  Deriv = - QOutMtxUnsDmCp(id,ic) * SorpAlfa(Lay) /     &
     &                    (ThtSrpRefDmCp(id,ic) - Theta(ic))
                  Deriv = Deriv * moiscap(ic,H(ic))
               else
                  Deriv = 0.d0
                  if (abs(DelHDmCp(id,ic)).gt.1.d-14) then
                     Deriv = - QOutMtxUnsDmCp(id,ic) / DelHDmCp(id,ic)
                  endif
               endif
               dFdhMp(ic)= dFdhMp(ic) + Deriv   
            endif
         endif
      end do
!
      case default
         call fatalerr ('Absorption', 'Illegal value for TASK')
      end select

      return
      end


!=====================================================================
      SUBROUTINE RAPIDDRAIN(id,ICpBtDm,ICpTpWaSrDm,flDraTub,FrMpWalWet, &
     &              FrReduQ,KDCrRlRef,VlMpDmCp,WaSrMpDm,ZBtDm,ZWaLevDm, &
     &              FlwOutDrRapCpPot,FlwOutDrRapPot,VlMpUndrDrL)
! ----------------------------------------------------------------------
!     Date               : April 2008                                      
!     Purpose            : To calculate rapid drainage from Main Bypass flow domain
! ----------------------------------------------------------------------
      use Variables
      implicit NONE

! --- global                                                          In
      integer ICpBtDm, ICpTpWaSrDm, id
      real(8) FrMpWalWet(MaDm,MaCp), FrReduQ, KDCrRlRef(MaDr) 
      real(8) VlMpDmCp(MaDm,MaCp), WaSrMpDm, ZBtDm, ZWaLevDm
      logical flDraTub(Madr)
!     -                                                              Out
      real(8) FlwOutDrRapCpPot(MaCp), FlwOutDrRapPot, VlMpUndrDrL
      real(8) FacRes
! ----------------------------------------------------------------------
! --- local
      integer ic
      real(8) DelH, KDCrRl, KDCrRlCp(MaCp), RapDraRes
      real(8) VOLUNDR, WaSrDrainabl, WthCr
! ----------------------------------------------------------------------
!
! --- Check whether this subroutine is relevant for this macropore domain
!     Rapid drainage only in domain 1: Main Bypass flow domain
      if (id.eq.1) then
!   - Initialization
         FlwOutDrRapPot= 0.d0
         do 100 ic= 1, NumNod
            FlwOutDrRapCpPot(ic)= 0.d0
 100     continue
!
         if (SwDrRap.eq.1) then
            continue
         else
            return
         endif
      else
         return
      endif
!
! --- for drain tube rapid drainage only possible if ZBtDm .le. Zdrabas;
!   - for open drain, rapid drainage also possible when ZBtDm > Zdrabas
      if (ZBtDm.Lt.ZDraBas+1.0d-1 .or. .not.flDraTub(1)) then
! transmissivity kD
         KDCrRl= 0.d0
         do 200 ic= ICpTpWaSrDm, ICpBtDm
            WthCr= DiPoCp(ic) * (1.d0 - dsqrt(1.d0-VlMpDmCp(1,ic) /     &
     &             DZ(ic)))
            KDCrRlCp(ic)= ((WthCr**RapDraReaExp)/DiPoCp(ic)) * DZ(ic)     
            if (ic.eq.ICpTpWaSrDm) KDCrRlCp(ic)=                        &
     &                                FrMpWalWet(id,ic) * KDCrRlCp(ic)
            KDCrRl= KDCrRl + KDCrRlCp(ic)     
 200     continue
!
         VlMpUndrDrL= 0.d0
         if (ZBtDm.lt.ZDraBas) VlMpUndrDrL= VOLUNDR(                    &
     &                               1,ICpBtDm,ZDraBas,VlMpDmCp,ZBtDm)
         WaSrDrainabl= dmax1(0.d0,(WaSrMpDm - VlMpUndrDrL))
!
! for open drain: DelH =  zwalevdm - ditch water level or bottom of MB domain !
         DelH=  ZWaLevDm - dmax1(ZDraBas,ZBtDm)
         if (ZWaLevDm.gt.-1.d-7) DelH= DelH + Pond
! ! in this version: no infiltration
         DelH= dmax1(DelH,0.d0)
         if (KdCrRl.gt.1.d-10) then
            FacRes        = dmin1(KDCrRlRef(1)/KDCrRl,1.1d0)
            RapDraRes     = RapDraResRef(1) * FacRes
            FlwOutDrRapPot= FrReduQ * (DelH / RapDraRes) * DT
         else
            FlwOutDrRapPot= 0.d0
         endif
         FlwOutDrRapPot= dmin1(FlwOutDrRapPot,WaSrDrainabl)
! distribute over compartments according to kD
         do 300 ic= ICpTpWaSrDm, ICpBtDm
            if (KDCrRl.gt.1.d-15) then     
               FlwOutDrRapCpPot(ic)= FlwOutDrRapPot * KDCrRlCp(ic) /    & 
     &                               KDCrRl  
            else   
               FlwOutDrRapCpPot(ic)= 0.d0
               FlwOutDrRapPot      = 0.d0
            endif    
 300     continue
      endif

      return
      END


!=======================================================================
      real(8) FUNCTION VOLUNDR(id,ICpBtDm,Level,VlMpDmCp,ZBtDm)
! ----------------------------------------------------------------------
!     Date               : september 2002
!     Purpose            : calculates volume of macropores under drain or 
!                        : groundwater level
! ----------------------------------------------------------------------
      use Variables
      implicit NONE

! --- global
      integer ICpBtDm, id
      real(8) Level, VlMpDmCp(MaDm,MaCp), ZBtDm 
! ----------------------------------------------------------------------
! --- local
      integer ic
      real(8) ZHlp
! ----------------------------------------------------------------------
! --- macropore volume under drain or groundwater level
      VOLUNDR= 0.d0
      ic= ICpBtDm + 1
      ZHlp= ZBtDm
      do 10 while(ZHlp.lt.Level .and. ic.gt.1)
         ic= ic - 1
         ZHlp= ZHlp + DZ(ic)
         VOLUNDR= VOLUNDR + VlMpDmCp(id,ic)
  10  continue
      if (VOLUNDR.gt.0.d0) VOLUNDR =                                    &
     &                     VOLUNDR - (ZHlp-Level)*VlMpDmCp(id,ic)/DZ(ic)
      return
      END
