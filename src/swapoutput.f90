! File VersionID:
!   $Id: swapoutput.f90 377 2018-04-04 10:57:55Z heine003 $
! ----------------------------------------------------------------------
      subroutine swapoutput(task)
! ----------------------------------------------------------------------
!     Date               : Aug 2004
!     Purpose            : open and write general swap output files
! ----------------------------------------------------------------------

      use Variables
      implicit none

      integer task

      select case (task)
      case (1)

! === open output files ===============================

! --     bal file
         if (swbal .eq. 1) call outbal (task)

! --     blc file
         if (swblc .eq. 1) call outblc(task,                                                 &
                      (CQMpInIntSatDm1+CQMpInIntSatDm2+CQMpInMtxSatDm1+CQMpInMtxSatDm2),     &
                      (CQMpOutMtxSatDm1+CQMpOutMtxSatDm2+CQMpOutMtxUnsDm1+CQMpOutMtxUnsDm2))

         return

      case (2)

! ===    write actual data ===============================

         if (flbaloutput) then

! --        bal file
            if (swbal .eq. 1) call outbal (task)

! --        blc file
            if (swblc .eq. 1) call outblc(task,                                                 &
                         (CQMpInIntSatDm1+CQMpInIntSatDm2+CQMpInMtxSatDm1+CQMpInMtxSatDm2),     &
                         (CQMpOutMtxSatDm1+CQMpOutMtxSatDm2+CQMpOutMtxUnsDm1+CQMpOutMtxUnsDm2))

         endif

      return

      case (3)

! ===    close output files ================================

         close (bal)
         if (swblc .eq. 1) close (blc)

! ---    final message log file
         write(logf,'(/,a)') ' Swap simulation okay!'
         close (logf)

      case default
         call fatalerr ('SWAPoutput', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine soilwateroutput(task)
! ----------------------------------------------------------------------
!     Date               : Aug 2004
!     Purpose            : open and write soil water output files
! ----------------------------------------------------------------------

      use Variables
      use SWAP_csv_output
      use SWAP_csv_output_tz
      implicit none

      integer task

      select case (task)

      case (1)
! ===    open output files and write headers ===============================

! --     user-defined variables in CSV file
         if (swcsv == 1) call csv_out(1)                    ! call csv_write(1)
! --     user-defined variables in CSV file
         if (swcsv_tz == 1) call csv_out_tz(1)              ! call csv_write_tz(1)

! --     wba file
         if (swwba.eq.1) call outwba (1)

! --     inc file
         if (swinc.eq.1) call outinc (1)

! --     str file
         if (swstr.eq.1) call outstr (1)

! --     vap file
         if (swvap.eq.1) call outvap (1)

! --     rot file, only when drought stress according to De Jong van Lier
         if (swdrought.eq.2) call outrot(1)

! --     capillary rise output file
         if(swcapriseoutput) call capriseoutput(task)

! --     extensive formatted output file for solute studies
         if (swafo.ge.1) call outafo (task)

! --     generate file with soil physical parameters
         if (swini .eq. 1) call outsoilphys ()

! --     extensive unformatted output file for solute studies
         if (swaun.ge.1) call outaun (task)

! --     special output for RUME project
         if (swrum == 1) call outrume (task)

! --     output of recharge/storage for modflow
         if(swoutputmodflow.eq.1) call OutputModflow(1)

      case (2)
! ===    write actual data ===============================

! --     user-defined variables in CSV file
         if (swcsv == 1) call csv_out(2)                    ! csv_write(2)
! --     user-defined variables in CSV file
         if (swcsv_tz == 1) call csv_out_tz(2)              ! csv_write_tz(2)

! --     wba file
         if (swwba.eq.1) call outwba (2)

! --     inc file
         if (swinc.eq.1) call outinc (2)

! --     str file
         if (flCropCalendar .and. flCropOutput) then
            if (swstr.eq.1) call outstr (2)
         endif

! --     vap file
         if (swvap.eq.1) call outvap (2)

! --     rot file
         if (swdrought.eq.2 .and. ptra .gt. 1.0d-10) call outrot (2)

! --     capillary rise output file
         if(swcapriseoutput) call capriseoutput(task)

! --     extensive formatted output file for solute studies
         if (swafo.ge.1) call outafo (task)

! --     extensive unformatted output file for solute studies
         if (swaun.ge.1) call outaun (task)

! --     special output for RUME project
         if (swrum == 1) call outrume (task)

! --     output of recharge/storage for modflow
         if(swoutputmodflow.eq.1) call OutputModflow(2)

      case (3)
! ===    write final values end of a simulation day ===========================

! ---    write final pressure heads, solute concentrations and soil temperatures
         if (swend.ge.1) call outend ()

      case (4)
! ===    close output files ===========================

! --     user-defined variables in CSV file
         if (swcsv == 1) call csv_out(3)                    ! csv_write(3)
! --     user-defined variables in CSV file
         if (swcsv_tz == 1) call csv_out_tz(3)              ! csv_write_tz(3)

         close (wba)
         close (inc)
         close (str)
         if (swvap.eq.1) close (vap)
         if (swafo.ge.1) close (afo)
         if (swaun.ge.1) close (aun)
         if (swoutputmodflow.eq.1) call OutputModflow(3)
         if (swcapriseoutput) call capriseoutput(3)

! --     special output for RUME project
         if (swrum == 1) call outrume (3)

      case default
         call fatalerr ('SoilWaterOutput', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outwba (task)
! ----------------------------------------------------------------------
!     date               : July 2002
!     purpose            : write water balance data to outnam.wba file
! ---------------------------------------------------------------------
      use variables, only: wba,daynr,daycum,swscre,cevap,cgird,cgrai,csnrai,cnird,cnrai,cpeva,cptra,cqbot,cqdra,cqrot,crunon,                &
                           crunoff,gwl,cQMpOutDrRap,pond,t1900,date,volact,volini,wbalance,outfil,pathwork,project,flprintshort,floutput,       &
                           swsnow,cqprai,ssnow,snowinco,PondIni,flheader,flmacropore
      implicit none

! -   global
      integer   task

! -   local
      integer   getun
      real(8)   dstor
      character(len=80) filnam,filtext
      character(len=10) gwlout
      character(len=19) datexti
      character(len=1)  comma

! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)

! --- open output file
      filnam = trim(pathwork)//trim(outfil)//'.wba'
      wba = getun (20,90)
      call fopens(wba,filnam,'new','del')
      filtext = 'cumulative water balance components (cm)'
      call writehead (wba,1,filnam,filtext,project)
      if (swscre .eq. 1) then
! ---   header of screen
        filtext = 'quick view water balance components'
        call writehead (5,1,'screen',filtext,project)
      endif

! --- write header in wba file
      if (flprintshort) then
        write (wba,10)
      else
        if(FlMacropore) then
          write (wba,13)
        else
          write (wba,12)
        endif
      endif
 10     format ('*',/,                                                  &
     &  '       Date,    Time,Day,  Dcum,Rain_g,Rain_n,Irr_g,Irr_n,',   &
     &  '  RunOn, RunOff,  Tpot,  Tact,  Epot,  Eact,  Drain,    Bot,', &
     &  '  DStor,      Gwl ,  Pond, Wbal,      Date2')
 12     format ('*',/,                                                  &
     &  '       Date,Day,  Dcum,Rain_g,Rain_n,Irr_g,Irr_n,  RunOn,',    &
     &  ' RunOff,  Tpot,  Tact,  Epot,  Eact,   Drain,     Bot,  ',     &
     &  ' DStor,       Gwl,  Pond, Wbal,      Date2')
 13     format ('*',/,                                                  &
     &  '       Date,Day,  Dcum,Rain_g,Rain_n,Irr_g,Irr_n,  RunOn,',    &
     &  ' RunOff,  Tpot,  Tact,  Epot,  Eact,   Drain, RapDrain,     ', &
     &  'Bot,   DStor,       Gwl,  Pond, Wbal,      Date2')


! --- write header screen output
      if (swscre.eq.1) then
        write (*,20)
 20     format(/,t20,'cumulative water balance components (cm)')
        write (*,22)
 22     format (/,'        date   rain  irrig runoff transp evapor',    &
     &          '  drain   qbot    gwl   wbal'/,                        &
     &          '              gross  gross        actual actual',      &
     &          '    net    net           cum')
      endif
      return

      case (2)

! === write output water balance components ===============================

! --- write header in case of new balance period
      if (flheader) then
        if (flprintshort) then
          write (wba,10)
        else
          write (wba,12)
        endif
      endif

! --- determine date and date-time
      call dtdpst ('year-month-day,hour:minute:seconds',t1900,datexti)

! --- write output record wba file
      gwlout = "          "
      if (gwl.lt.998.0d0)  write(gwlout,'(f9.1)') gwl
      if (swsnow.eq.0) then
        dstor = (volact + pond) - (volini + PondIni)
        if (flprintshort) then
          write (wba,25) datexti,comma,daynr,comma,daycum,comma,        &
     &    cgrai+csnrai,comma,cnrai,comma,cgird,comma,cnird,comma,crunon,&
     &    comma,crunoff,comma,cptra,comma,cqrot,comma,cpeva,comma,cevap,&
     &    comma,(cqdra+cQMpOutDrRap),comma,cqbot,comma,dstor,comma,     &
     &    gwlout,comma,pond,comma,wbalance,comma,date
        else
          if(FlMacropore) then
            write (wba,30) date,comma,daynr,comma,daycum,comma,         &
     &    cgrai+csnrai,comma,cqprai,comma,cgird,comma,cnird,comma,      &
     &    crunon,comma,crunoff,comma,cptra,comma,cqrot,comma,cpeva,     &
     &    comma,cevap,comma,cqdra,comma,cQMpOutDrRap,comma,cqbot,comma, &   !!! aanpassing GEM
     &    dstor,comma,gwlout,comma,pond,comma,wbalance,comma,date
          else
            write (wba,31) date,comma,daynr,comma,daycum,comma,         &
     &    cgrai+csnrai,comma,cqprai,comma,cgird,comma,cnird,comma,      &
     &    crunon,comma,crunoff,comma,cptra,comma,cqrot,comma,cpeva,     &
     &    comma,cevap,comma,(cqdra+cQMpOutDrRap),comma,cqbot,comma,     &
     &    dstor,comma,gwlout,comma,pond,comma,wbalance,comma,date
          endif
        endif
      else
        dstor = (volact + pond + ssnow) - (volini + PondIni + snowinco)
        if (flprintshort) then
          write (wba,25) datexti,comma,daynr,comma,daycum,comma,        &
     &    cgrai+csnrai,comma,cnrai,comma,cgird,comma,cnird,comma,crunon,&
     &    comma,crunoff,comma,cptra,comma,cqrot,comma,cpeva,comma,cevap,&
     &    comma,(cqdra+cQMpOutDrRap),comma,cqbot,comma,dstor,comma,     &
     &    gwlout,comma,pond,comma,wbalance,comma,date
        else
          if(FlMacropore) then
            write (wba,30) date,comma,daynr,comma,daycum,comma,         &
     &    cgrai+csnrai,comma,cqprai,comma,cgird,comma,cnird,comma,      &
     &    crunon,comma,crunoff,comma,cptra,comma,cqrot,comma,cpeva,     &
     &    comma,cevap,comma,cqdra,comma,cQMpOutDrRap,comma,cqbot,comma, &   !!! aanpassing GEM
     &    dstor,comma,gwlout,comma,pond,comma,wbalance,comma,date
          else
            write (wba,31) date,comma,daynr,comma,daycum,comma,         &
     &    cgrai+csnrai,comma,cqprai,comma,cgird,comma,cnird,comma,      &
     &    crunon,comma,crunoff,comma,cptra,comma,cqrot,comma,cpeva,     &
     &    comma,cevap,comma,(cqdra+cQMpOutDrRap),comma,cqbot,comma,     &
     &    dstor,comma,gwlout,comma,pond,comma,wbalance,comma,date
          endif
        endif
      endif

 25   format (a19,a1,i3,a1,i6,2(a1,f7.2),2(a1,f5.1),a1,f7.2,a1,         &
     &     f7.2,a1,f7.2,3(a1,f8.2),3(a1,f7.2),2a,a1,f6.2,a1,f5.2,a1,a11)

 30   format (a11,a1,i3,a1,i6,2(a1,f7.2),2(a1,f5.1),a1,f7.2,a1,f7.2,    &
     &      4(a1,f6.2),4(a1,f8.3),2a,a1,f6.2,a1,f5.2,a1,a11)
 31   format (a11,a1,i3,a1,i6,2(a1,f7.2),2(a1,f5.1),a1,f7.2,a1,f7.2,    &
     &      4(a1,f6.2),3(a1,f8.3),2a,a1,f6.2,a1,f5.2,a1,a11)

! --- write output record screen
      if (swscre.eq.1 .and. floutput) then
        if (flheader) write (*,22)
        write (unit=*, fmt=40)                                          &
     &                date,cgrai,cgird,crunoff,cqrot,cevap,             &
     &                (cqdra+cQMpOutDrRap),cqbot,gwl,wbalance
 40     format(1x,a11,f8.1,6f7.2,f7.1,f7.2)
      endif

      case default
         call fatalerr ('outwba', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outinc (task)
! ----------------------------------------------------------------------
!     date               : july 2002
!     purpose            : write water balance increments to outnam.inc file
! ---------------------------------------------------------------------
      use variables, only: inc,daynr,daycum,igrai,isnrai,igsnow,igird,iintc,irunon,iruno,irunoCN,iptra,iqrot,ipeva,ievap,isubl,iqdra,  &
                           iQMpOutDrRap,iqbot,t1900,date,outfil,pathwork,project,flheader,gwl,volact,volini,pond,PondIni,ssnow,snowinco,flprintshort
      implicit none

! --- global
      integer   task

! --- local
      integer   getun
      real(8)  baldev,dstor
      character(len=80) filnam,filtext
      character(len=1)  comma
      character(len=10) gwlout
      character(len=19) datexti

      real(8), save ::   VolOld,PondOld,SnowOld


! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)

      t1900 = t1900

! --- open output file once
      filnam = trim(pathwork)//trim(outfil)//'.inc'
      inc = getun (20,90)
      call fopens(inc,filnam,'new','del')
      filtext = 'water balance increments (cm/day)'
      call writehead (inc,1,filnam,filtext,project)

! --- write header of inc file
      if (flprintshort) then
        write (inc,10)
      else
        write (inc,12)
      endif
  10  format('*',/,                                                     &
     & '      Date,    Time,Day,  Dcum,      Rain,     Snow,',          &
     & '      Irrig,    Interc,     Runon,    Runoff,      Tpot,',      &
     & '      Tact,      Epot,      Eact,  Drainage,   QBottom,',       &
     & '       Gwl,  dstorage,    baldev')
  12  format('*',/,                                                     &
     & '       Date,Day,  Dcum,      Rain,     Snow,',                  &
     & '      Irrig,    Interc,     Runon,    Runoff,      Tpot,',      &
     & '      Tact,      Epot,      Eact,  Drainage,   QBottom,',       &
     & '       Gwl,  dstorage,    baldev')


      VolOld  = volini
      PondOld = PondIni
      SnowOld = snowinco

      return

      case (2)

! --- write header in case of new balance period
      if (flheader) then
        if (flprintshort) then
          write (inc,10)
        else
          write (inc,12)
        endif
      endif

! --- determine date and date-time
      call dtdpst ('year-month-day,hour:minute:seconds',t1900,datexti)

! --- write output record .inc file
      gwlout = "          "
      if (gwl.lt.998.0d0)  write(gwlout,'(f9.1)') gwl
      dstor = (volact + pond + ssnow) - (VolOld + PondOld + SnowOld)
      baldev = (igrai+isnrai+igsnow+igird+irunon) - dstor -             &
     & (iintc+iruno+irunoCN+iqrot+ievap+isubl+iQMpOutDrRap+iqdra+(-1.0d0*iqbot))
      if (flprintshort) then
        write (inc,20) datexti,comma,daynr,comma,daycum,comma,          &
     &    igrai+isnrai,comma,igsnow,comma,igird,comma,iintc,comma,      &
     &    irunon,comma,iruno+irunoCN,comma,iptra,comma,iqrot,comma,ipeva,comma, &
     &    ievap,comma,(iQMpOutDrRap+iqdra),comma,iqbot,                 &
     &    comma,gwlout,comma,dstor,comma,baldev            !comma,storage
      else
        write (inc,22) date,comma,daynr,comma,daycum,comma,             &
     &    igrai+isnrai,comma,igsnow,comma,igird,comma,iintc,comma,      &
     &    irunon,comma,iruno+irunoCN,comma,iptra,comma,iqrot,comma,ipeva,comma, &
     &    ievap,comma,(iQMpOutDrRap+iqdra),comma,iqbot,                 &
     &    comma,gwlout,comma,dstor,comma,baldev           !comma,storage
      endif
 20   format (a19,a1,i3,a1,i6,12(a1,f10.5),2a,2(a1,f10.5))     !,(a1,e12.5)
 22   format (a11,a1,i3,a1,i6,12(a1,f10.5),2a,2(a1,f10.5))     !,(a1,e12.5)

      VolOld = volact
      PondOld = pond
      SnowOld = ssnow

      case default
         call fatalerr ('outinc', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outstr (task)
! ----------------------------------------------------------------------
!     date               : March 2008
!     purpose            : write ETpot and stress factors to outfil.str file
! ---------------------------------------------------------------------
      use variables, only: str,daynr,daycum,ies0,iet0,iew0,ipeva,iptra,iqrot,iqredwet,iqreddry,iqredsol,iqredfrs,t1900,   &
                           date,outfil,pathwork,project,flheader,flprintshort

      implicit none

! --- global
      integer   task

! --- local
      integer   getun
      character(len=80)  filnam
      character(len=132) filtext
      character(len=1)   comma
      character(len=19)  datexti
! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)

      t1900 = t1900

! --- open output file once
      filnam = trim(pathwork)//trim(outfil)//'.str'
      str = getun (20,90)
      call fopens(str,filnam,'new','del')
      filtext = 'ES0, ET0, EW0, Epot, Tpot, Tact and 4 Tact-stress '//  &
     & 'values for wetness, drought, salinity and frost (cm/period)'
      call writehead (str,1,filnam,filtext,project)

! --- write header of str file
      if (flprintshort) then
        write (str,10)
      else
        write (str,12)
      endif
  10  format('*',/,                                                     &
     & '       Date,    Time,Day,  Dcum,ESoilWet,Tcropdry,Tcropwet,',   &
     & '    Epot,    Tpot,    Tact, Tredwet, Treddry, Tredsol, Tredfrs')
  12  format('*',/,                                                     &
     & '       Date,Day,  Dcum,ESoilWet,Tcropdry,Tcropwet,',            &
     & '    Epot,    Tpot,    Tact, Tredwet, Treddry, Tredsol, Tredfrs')
      return

      case (2)

! --- write header in case of new balance period
      if (flheader) then
        if (flprintshort) then
          write (str,10)
        else
          write (str,12)
        endif
      endif

! --- determine date and date-time
      call dtdpst ('year-month-day,hour:minute:seconds',t1900,datexti)

! --- write output record .str file
      if (flprintshort) then
        write (str,20) datexti,comma,daynr,comma,daycum,                &
     &                 comma,ies0,comma,iet0,comma,iew0,comma,ipeva,    &
     &                 comma,iptra,comma,iqrot,comma,iqredwet,          &
     &                 comma,iqreddry,comma,iqredsol,comma,iqredfrs
      else
        write (str,22) date,comma,daynr,comma,daycum,                   &
     &                 comma,ies0,comma,iet0,comma,iew0,comma,ipeva,    &
     &                 comma,iptra,comma,iqrot,comma,iqredwet,          &
     &                 comma,iqreddry,comma,iqredsol,comma,iqredfrs
      endif
 20   format (a19,a1,i3,a1,i6,10(a1,f8.4))
 22   format (a11,a1,i3,a1,i6,10(a1,f8.4))

      case default
         call fatalerr ('outstr', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outvap (task)
! ----------------------------------------------------------------------
!     date               : December 2007
!     purpose            : write output of soil profile data
! ---------------------------------------------------------------------
      use variables, only: ztopcp, zbotcp, vap,daynr,numnod,daycum,z,cml,t1900,theta,h,k,tsoil,q,outfil,     &
                           pathwork,project,swheader,isqtop,isqbot,qdraincomp,qrot,cmsy,date,flprintshort
      implicit none

! --- global
      integer   task

! --- local
      integer   node,getun
      real(8)   sflux
      character(len=19) datexti
      character(len=11) inidate
      character(len=80) filnam,filtext
      character(len=1)  comma

! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)
! --- open output file
      filnam = trim(pathwork)//trim(outfil)//'.vap'
      vap = getun (20,90)
      call fopens(vap,filnam,'new','del')
      filtext = 'soil profile data'
      call writehead (vap,1,filnam,filtext,project)
      write (vap,100)

! --- write header in vap file
      if (flprintshort) then
        write (vap,200)
      else
        write (vap,210)
      endif

! --- write initial profile data to vap file

      if (flprintshort) then
! ---   determine date and date-time
        call dtdpst ('year-month-day,hour:minute:seconds',t1900,datexti)
        do node = 1,numnod
          if (node.eq.1) then
             sflux = isqtop
              else
             sflux = 0.5d0 * (cml(node) + cml(node-1)) * q(node)
           end if
           write (vap,300) datexti,comma,z(node),comma,theta(node),     &
     &       comma,h(node),comma,k(node),comma,qdraincomp(node),comma,  &
     &       qrot(node),comma,q(node),comma,tsoil(node),comma,cml(node),&
     &       comma,cmsy(node),comma,sflux,comma,ztopcp(node),           &
     &       comma,zbotcp(node),comma,daynr,comma,daycum
        end do
        sflux = isqbot
        write (vap,400) datexti,comma,zbotcp(numnod),                   &
     &    comma,comma,comma,comma,comma,                                &
     &    comma,q(numnod+1),comma,comma,comma,comma,sflux,comma,        &
     &    zbotcp(numnod),comma,comma,daynr,comma,daycum

      else
! ---   determine initial date
        call dtdpst ('year-month-day',t1900-0.1d0,inidate)
        do node = 1,numnod
          if (node.eq.1) then
             sflux = isqtop
              else
             sflux = 0.5d0 * (cml(node) + cml(node-1)) * q(node)
           end if
           write (vap,310) inidate,comma,z(node),comma,theta(node),     &
     &       comma,h(node),comma,k(node),comma,qdraincomp(node),comma,  &
     &       qrot(node),comma,q(node),comma,tsoil(node),comma,cml(node),&
     &       comma,cmsy(node),comma,sflux,comma,ztopcp(node),           &
     &       comma,zbotcp(node),comma,daynr,comma,daycum
        end do
        sflux = isqbot
        write (vap,410) inidate,comma,zbotcp(numnod),                   &
     &    comma,comma,comma,comma,comma,                                &
     &    comma,q(numnod+1),comma,comma,comma,comma,sflux,comma,        &
     &    zbotcp(numnod),comma,comma,daynr,comma,daycum

      endif

      return

      case (2)

! === write actual profile data ===========================================


! --- write header in vap file
      if (swheader .eq. 1) then
        if (flprintshort) then
          write (vap,200)
        else
          write (vap,210)
        endif
      endif

      if (flprintshort) then
! ---   determine date and date-time
        call dtdpst ('year-month-day,hour:minute:seconds',t1900,datexti)
        do node = 1,numnod
          if (node.eq.1) then
             sflux = isqtop
              else
             sflux = 0.5d0 * (cml(node) + cml(node-1)) * q(node)
           end if
           write (vap,300) datexti,comma,z(node),comma,theta(node),     &
     &       comma,h(node),comma,k(node),comma,qdraincomp(node),comma,  &
     &       qrot(node),comma,q(node),comma,tsoil(node),comma,cml(node),&
     &       comma,cmsy(node),comma,sflux,comma,ztopcp(node),           &
     &       comma,zbotcp(node),comma,daynr,comma,daycum
        end do
        sflux = isqbot
        write (vap,400) datexti,comma,zbotcp(numnod),                   &
     &    comma,comma,comma,comma,comma,                                &
     &    comma,q(numnod+1),comma,comma,comma,comma,sflux,comma,        &
     &    zbotcp(numnod),comma,comma,daynr,comma,daycum

      else
        do node = 1,numnod
          if (node.eq.1) then
             sflux = isqtop
              else
             sflux = 0.5d0 * (cml(node) + cml(node-1)) * q(node)
           end if
           write (vap,310) date,comma,z(node),comma,theta(node),        &
     &       comma,h(node),comma,k(node),comma,qdraincomp(node),comma,  &
     &       qrot(node),comma,q(node),comma,tsoil(node),comma,cml(node),&
     &       comma,cmsy(node),comma,sflux,comma,ztopcp(node),           &
     &       comma,zbotcp(node),comma,daynr,comma,daycum
        end do
        sflux = isqbot
        write (vap,410) date,comma,zbotcp(numnod),                      &
     &    comma,comma,comma,comma,comma,                                &
     &    comma,q(numnod+1),comma,comma,comma,comma,sflux,comma,        &
     &    zbotcp(numnod),comma,comma,daynr,comma,daycum

      endif

 100  format(                                                           &
     & '* Explanation:   instantaneous fluxes of drainage, ',           &
     & 'root extraction, water and solute ',/                           &
     & '*                fluxes of water and solute apply to ',         &
     & 'top of compartment',/                                           &
     & '*                solute1 = concentration in soil water',/       &
     & '*                solute2 = total concentration ',               &
     & '(dissolved + adsorbed)')

 200  format(/,                                                         &
     & '                         cm,  cm3/cm3,         cm,       cm/d,',&
     & '       cm/d,       cm/d,       cm/d,     ºC,     mg/cm3,     ', &
     &'mg/cm3,   mg/cm2/d,     cm,     cm,  nr,   nr',/                 &
     & '      date,    time,  depth, wcontent,      phead,    hconduc,',&
     & '   drainage,    rootext,  waterflux,   temp,    solute1,    ',  &
     & 'solute2, soluteflux,    top, bottom, day, dcum')

 210  format(/,                                                         &
     & '                 cm,  cm3/cm3,         cm,       cm/d,',        &
     & '       cm/d,       cm/d,       cm/d,     ºC,     mg/cm3,     ', &
     &'mg/cm3,   mg/cm2/d,     cm,     cm,  nr,   nr',/                 &
     & '       date,  depth, wcontent,      phead,    hconduc,',        &
     & '   drainage,    rootext,  waterflux,   temp,    solute1,    ',  &
     & 'solute2, soluteflux,    top, bottom, day, dcum')

 300  format(a19,a1,f7.1,a1,f9.3,1p,5(a1,e11.3),a1,0p,f7.2,1p,          &
     &       3(a1,e11.3),a1,0p,2(f7.1,a1),i4,a1,i5)

 310  format(a11,a1,f7.1,a1,f9.3,1p,5(a1,e11.3),a1,0p,f7.2,1p,          &
     &       3(a1,e11.3),a1,0p,2(f7.1,a1),i4,a1,i5)

 400  format(a19,a1,f7.1,a1,t38,a1,t50,a1,t62,a1,t74,a1,t86,a1,1p,e11.3,&
     &       a1,t106,a1,t118,a1,t130,a1,e11.3,a1,0p,f7.1,a1,t158,a1,i4, &
     &       a1,i5)

 410  format(a11,a1,f7.1,a1,9x,4(a1,11x),a1,1p,e11.3,a1,7x,2(a1,11x),a1,&
     &       e11.3,a1,0p,f7.1,a1,7x,a1,i4,a1,i5)

      case default
         call fatalerr ('outvap', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outrot (task)

! ----------------------------------------------------------------------
!     date               : February 2012
!     purpose            : write output of microscopic root water uptake
! ---------------------------------------------------------------------
      use variables, only: ztopcp, zbotcp,rot,daynr,noddrz,daycum,z,hxylem,t1900,theta,hm1,q,outfil,pathwork,project,swheader,qrot,      &
                           date,flprintshort,hroot,inq,inqrot,mroot,mflux,rootrho,rootphi,hleaf
      implicit none

! --- global
      integer   task

! --- local
      integer   node,getun
      character(len=19) datexti
!      character(len=11) inidate
      character(len=80) filnam,filtext
      character(len=1)  comma

! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)


! --- open output file
      filnam = trim(pathwork)//trim(outfil)//'.rot'
      rot = getun (20,90)
      call fopens(rot,filnam,'new','del')
      filtext = 'microscopic root water uptake'
      call writehead (rot,1,filnam,filtext,project)
      write (rot,100)

! --- write header in rot file
      if (swheader .eq. 0) then
        if (flprintshort) then
          write (rot,200)
        else
          write (rot,210)
        endif
      endif

      return

      case (2)

! === write actual profile data ===========================================

! --- write header in rot file
      if (swheader .eq. 1) then
        if (flprintshort) then
          write (rot,200)
        else
          write (rot,210)
        endif
      endif

      if (flprintshort) then
! ---   determine date and date-time
        call dtdpst ('year-month-day,hour:minute:seconds',t1900,datexti)
        do node = 1,noddrz
           write (rot,300) datexti,comma,z(node),comma,hleaf,comma,     &
     &       hxylem,comma,                                              &
     &       hroot(node),comma,hm1(node),comma,inqrot(node),comma,      &
     &       qrot(node),comma,inq(node),comma,q(node),comma,mroot(node),&
     &       comma,mflux(node),comma,rootrho(node),comma,rootphi(node), &
     &       comma,theta(node),comma,ztopcp(node),comma,                &
     &       zbotcp(node),comma,daynr,comma,daycum
        end do

      else
        do node = 1,noddrz
           write (rot,310) date,comma,z(node),comma,hleaf,comma,        &
     &       hxylem,comma,                                              &
     &       hroot(node),comma,hm1(node),comma,inqrot(node),comma,      &
     &       qrot(node),comma,inq(node),comma,q(node),comma,mroot(node),&
     &       comma,mflux(node),comma,rootrho(node),comma,rootphi(node), &
     &       comma,theta(node),comma,ztopcp(node),comma,                &
     &       zbotcp(node),comma,daynr,comma,daycum
        end do

      endif

 100  format(                                                           &
     & '* Explanation:   fluxes of soil water (qsoilw and iqsoilw)',    &
     & ' apply to top of compartment;',/                                &
     & '*                both instantaneous (qroot and qsoilw) and',    &
     & ' incremental fluxes (qroot and qsoilw) are listed')

 200  format(/,                                                         &
     & '                         cm,         cm,         cm,',          &
     & '         cm,',                                                  &
     & '         cm,         cm,       cm/d,         cm,       cm/d,',  &
     & '      cm2/d,      cm2/d,       /cm2,        d/m,  cm3/cm3,',    &
     & '     cm,     cm,  nr,   nr',/                                   &
     & '       date,  depth,      hleaf,     hxylem,      hroot,',      &
     & '      hsoil,     iqroot,      qroot,    iqsoilw,     qsoilw,',  &
     & '      Mroot,      Msoil,    RootRho,    RootPhi, wcontent,',    &
     & '    top, bottom, day, dcum')

 210  format(/,                                                         &
     & '                 cm,         cm,         cm,         cm,',      &
     & '         cm,         cm,       cm/d,         cm,       cm/d,',  &
     & '      cm2/d,      cm2/d,       /cm2,        d/m,  cm3/cm3,',    &
     & '     cm,     cm,  nr,   nr',/                                   &
     & '       date,  depth,      hleaf,     hxylem,      hroot,',      &
     & '      hsoil,     iqroot,      qroot,    iqsoilw,     qsoilw,',  &
     & '      Mroot,      Msoil,    RootRho,    RootPhi, wcontent,',    &
     & '    top, bottom, day, dcum')

 300  format(a19,a1,f7.1,4(a1,f11.0),8(a1,e11.3),a1,f9.3,a1,2(f7.1,a1), &
     &       i4,a1,i5)

 310  format(a11,a1,f7.1,4(a1,f11.0),8(a1,e11.3),a1,f9.3,a1,2(f7.1,a1), &
     &       i4,a1,i5)

      case default
         call fatalerr ('outrot', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outbal (task)
! ----------------------------------------------------------------------
!     date               : July 2002
!     purpose            : write overview balances to bal file
! ---------------------------------------------------------------------
      use variables, only: zbotcp,bal,logf,swscre,swdra,numnod,nrlevs,swsolu,ioutdat,cevap,cgird,cgrai,cqbot,tstart,cqrot,crunoff,crunoffCN,     &
                           crunon,cQMpOutDrRap,cqdra,cqdrain,samini,sampro,samcra,sqprec,sqirrig,sqbot,dectot,rottot,sqrap,sqdra,pond,volact,      &
                           volini,t1900,outdat,outfil,pathwork,project,caintc,csubl,PondIni,WaSrDm1,WaSrDm2,WaSrDm1Ini,WaSrDm2Ini,      &
                           swsnow,cgsnow,csnrai,snowinco,ssnow,cqssdi
      implicit none

! --- global
      integer   task

! --- local
      integer   i,getun
      real(8)   precip
      character(len=80) filnam,filtext
      character(len=11) datbegin,datend
! ----------------------------------------------------------------------

      select case (task)
      case (1)
      logf = logf ! for Forcheck
      swscre = swscre
! --- open output file
      filnam = trim(pathwork)//trim(outfil)//'.bal'
      bal = getun (20,90)
      call fopens(bal,filnam,'new','del')
      filtext='overview of actual water and solute balance components'
      call writehead (bal,1,filnam,filtext,project)

      return

      case (2)

! --- begin date of balance period
      if (ioutdat .eq. 2) then
        call dtdpst ('year-month-day',tstart+0.1d0,datbegin)
      else
        call dtdpst ('year-month-day',                                &
     &                                outdat(ioutdat-2)+1.1d0,datbegin)
      endif

! --- end date of balance period
      call dtdpst ('year-month-day',t1900-0.9d0,datend)

! --- write output record
      write (bal,20) datbegin,datend
      write (bal,22) -zbotcp(numnod)

      if (swsolu .eq. 1) then
          write (bal,24) (volact+pond+WaSrDm1+WaSrDm2+ssnow),           &
     &                                            (sampro+samcra),      &
     &                  (volini+PondIni+WaSrDm1Ini+WaSrDm2Ini+snowinco),&
     &                                             samini,              &
     &                  (volact+pond+WaSrDm1+WaSrDm2+ssnow-             &
     &                   volini-PondIni-WaSrDm1Ini-WaSrDm2Ini-snowinco),&
     &                                            (sampro+samcra-samini)
      else
          write (bal,25) (volact+pond+WaSrDm1+WaSrDm2+ssnow),           &
     &                  (volini+PondIni+WaSrDm1Ini+WaSrDm2Ini+snowinco),&
     &                  (volact+pond+WaSrDm1+WaSrDm2+ssnow-             &
     &                   volini-PondIni-WaSrDm1Ini-WaSrDm2Ini-snowinco)
      endif

      precip = cgrai
      if (swsnow.eq.1) then
         precip = precip + cgsnow + csnrai
      endif

      if (swsnow.eq.1) then
         write (bal,26) precip,caintc,crunon,crunoff,cqssdi,crunoffCN,cgird,cqrot,cqbot, &
     &               (cevap+csubl),cQMpOutDrRap
      else
         write (bal,27) precip,caintc,crunon,crunoff,cqssdi,crunoffCN,cgird,cqrot,cqbot, &
     &               (cevap+csubl),cQMpOutDrRap
      endif

      if (swdra .ne. 0) then
        do i = 1,nrlevs
          write (bal,28) i,cqdrain(i)
        end do
      endif
      write(bal,30) (precip+cgird+cqbot+crunon+cqssdi),                        &
     &   (caintc+crunoff+crunoffCN+cqrot+cevap+csubl+cQMpOutDrRap+cqdra)

      if (swsolu .eq. 1) then
        write (bal,34) sqprec,dectot,sqirrig,rottot,sqbot,sqrap,sqdra
        write (bal,36) (sqprec+sqirrig+sqbot),                          &
     &    (dectot+rottot+sqrap+sqdra)
      endif

 20   format(/'Period',t20,':',t23,a11,' until  ',a11)
 22   format('Depth soil profile',t20,':',f8.2,' cm')
 24   format(/,T13,'Water storage',t36,'Solute storage',/,              &
     &   'Final',t9,':',t15,f8.2,' cm',t31,e12.4,' mg/cm2',/,           &
     &   'Initial',t9,':',t15,f8.2,' cm',t31,e12.4,' mg/cm2',/,         &
     &    t13,13('='),t33,17('='),/,                                    &
     &   'Change',t15,f8.2,' cm',t31,e12.4,' mg/cm2')
 25   format(/,T13,'Water storage',/,'Final',t9,':',t15,f8.2,' cm',/,   &
     &   'Initial',t9,':',t15,f8.2,' cm',/,                             &
     &    t13,13('='),/,'Change',t15,f8.2,' cm')
     26   format(//,'Water balance components (cm)',//,                 &
     &   'In',t30,'Out',/,25('='),t30,28('='),/,                        &
     &   'Rain',t16,':',f9.2,t30,'Interception',t48,':',f9.2,/,         &
     &   'Runon',t16,':',f9.2,t30,'Runoff',t48,':',f9.2,/,              &
     &   'SSDI',t16,':',f9.2,t30,'Runoff_CN',t48,':',f9.2,/,            &
     &   'Irrigation',t16,':',f9.2,t30,'Transpiration',t48,':',f9.2,/,  &
     &   'Bottom flux',t16,':',f9.2,t30,'Soil evaporation',t48,':',f9.2,&
     &   /,t30,'Crack flux',t48,':',f9.2)
     27   format(//,'Water balance components (cm)',//,                 &
     &   'In',t30,'Out',/,25('='),t30,28('='),/,                        &
     &   'Rain + snow',t16,':',f9.2,t30,'Interception',t48,':',f9.2,/,  &
     &   'Runon',t16,':',f9.2,t30,'Runoff',t48,':',f9.2,/,              &
     &   'SSDI',t16,':',f9.2,t30,'Runoff_CN',t48,':',f9.2,/,            &
     &   'Irrigation',t16,':',f9.2,t30,'Transpiration',t48,':',f9.2,/,  &
     &   'Bottom flux',t16,':',f9.2,t30,'Soil evaporation',t48,':',f9.2,&
     &   /,t30,'Crack flux',t48,':',f9.2)
 28   format(t30,'Drainage level',i2,t48,':',f9.2)
 30   format(25('='),t30,28('='),/,                                     &
     &   'Sum',t16,':',f9.2,t30,'Sum',t48,':',f9.2)
 34   format(//,'Solute balance components (mg/cm2)',//,                &
     &   'In',t30,'Out',/,25('='),t30,28('='),/,                        &
     &   'Rain',t13,':',e12.4,t30,'Decomposition',t45,':',e12.4,/,      &
     &   'Irrigation',t13,':',e12.4,t30,'Root uptake',t45,':',e12.4,/,  &
     &   'Bottom flux',t13,':',e12.4,t30,'Cracks',t45,':',e12.4,/,      &
     &   t30,'Drainage',t45,':',e12.4)
 36   format(25('='),t30,28('='),/,                                     &
     &   'Sum',t13,':',e12.4,t30,'Sum',t45,':',e12.4,/)

      case default
         call fatalerr ('outbal', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outend ()
! ----------------------------------------------------------------------
!     Updates            : September 2015, June 2017
!     date               : July 2002
!     purpose            : write final result to .end file
! ---------------------------------------------------------------------
      use variables, only: t1900,swend,numnod,h,flSolute,flAgeTracer,cml,z,fltemperature,tsoil,                                                             &
                           ssnow,pond,dt,icrop,croptype,cropfil,flSurfaceWater,wls,swredu,ldwet,spev,saev,outfil,pathwork,project,                          &
                           rd,rdpot,dvs,flanthesis,tsum,ilvold,ilvoldpot,wrt,wrtpot,tadw,tadwpot,wst,wstpot,wso,wsopot,wlv,wlvpot,laiexp,lai,laipot,        &
                           dwrt,dwrtpot,dwlv,dwlvpot,dwst,dwstpot,dwlvSoil,dwlvCrop,gasst,gasstpot,mrest,mrestpot,                                          &
                           cwdm,cwdmpot,sla,slapot,lvage,lvagepot,lv,lvpot,daycrop,nofd,atmin7,tsumgerm,rid,flgrazingpot,idregr,                            &
                           idregrpot,laiexppot,laimax,tagp,tagppot,idaysgraz,idaysgrazpot,daygrowth,daygrowthpot,                                           &
                           tagpt,tagptpot,iharvest,iseqgm,iseqgmpot,flgrazing,slw,cuptgraz,cuptgrazpot,flIrrigate,dayfix,                                   &
                           flCropCalendar,flCropPrep,flCropSow,flCropGerm,flCropEmergence,flCropHarvest,PrepDelay,SowDelay,swinter,sicact,glaiex,glaiexpot
      implicit none

! --- global

! --- local
      integer   i,fin,getun,swanthesis,swgrazing,swgrazingpot,count
      integer   swIrrigate
      integer   swCropPrep,swCropSow,swCropGerm
      integer   swCropEmergence,swCropHarvest
      character(len=300) filnam
      character(len=80)  filtext
      character(len=10) date
! ----------------------------------------------------------------------

! --- open output file for final profile data and write heading
      if (swend .eq. 1) then
         filnam = trim(pathwork)//trim(outfil)//'.end'
      else
         call dtdpst ('year-month-day',t1900-0.5d0,date)
         filnam = trim(pathwork)//trim(outfil)//'_'//date(1:4)//date(6:7)//date(9:10)//'.end'
      end if
      fin = getun (20,90)
      call fopens(fin,filnam,'new','del')
      filtext = 'final state variables'
      call writehead (fin,1,filnam,filtext,project)

      ! write snow- and ponding-layer
      write(fin,'(/,"* Snow layer (Ssnow in cm)")')
      write(fin,'(" Ssnow = ", e12.5)') ssnow
      write(fin,'(/,"* Liquid water in snow layer (cm)")')
      write(fin,'(" Slw = ", e12.5)') slw
      write(fin,'(/,"* Ponding layer (Pond in cm)")')
      write(fin,'(" Pond = ", e12.5)') pond

      ! write soil water pressure heads
      write(fin,'(/,"* Soil water pressure heads  (z in cm; h in cm)")')
      write(fin,'("       z_h            h")')
      do i = 1, numnod
        write (fin,'(f10.1," ",1p,e12.5)') z(i), h(i)
      end do

      ! write solute concentrations
      if(flSolute .or. flAgeTracer) then
        write(fin,'(/,"* Solute concentrations (z in cm; Cml in mg/cm3)")')
        write(fin,'("     z_Cml          Cml")')
        do i = 1, numnod
          write (fin,'(f10.1," ",1p,e12.5)') z(i), cml(i)
        end do
      endif

      ! write soil temperatures
      if(fltemperature) then
        write(fin,'(/,"* Soil temperatures  (z in cm; Tsoil in oC)")')
        write(fin,'("   z_Tsoil        Tsoil")')
        do i = 1, numnod
          write (fin,'(f10.1," ",1p,e12.5)') z(i), tsoil(i)
        end do
      endif

      ! write surface water level
      if(flSurfaceWater) then
        write(fin,'("* Surface water")')
        write(fin,'("* Surface water level (cm)")')
        write(fin,'(" wls = ", e12.5/)') wls
      endif

      ! write soil evaporation reservoirs
      if(swredu.eq.1) then
        write(fin,'(/,"* Soil evaporation reservoir (Black)")')
        write(fin,'("* Time after significant rainfall (d)")')
        write(fin,'(" ldwet = ", e12.5)') ldwet
      endif
      if(swredu.eq.2) then
        write(fin,'(/,"* Soil evaporation reservoir (Boesten/Stroosnijder)")')
        write(fin,'("* Rainfall excess (cm)")')
        write(fin,'(" spev = ", e12.5)') spev
        write(fin,'(/,"* Cumulative actual evaporation (cm)")')
        write(fin,'(" saev = ", e12.5)') saev
      endif

      ! write length of final timestep (d)
      write(fin,'(/,"* Timing parameters")')
      write(fin,'("* Length of final time step (d)")')
      write(fin,'(" dt = ", e12.5)') dt
      
      ! write irrigation
      if (flIrrigate) then
        swIrrigate = 1
      else
        swIrrigate = 0
      endif
      write(fin,'(/,"* Irrigation")')
      write(fin,'(" swIrrigate = ",i4)') swIrrigate
      if (flIrrigate) then
        write(fin,'(" dayfix = ",i4)') dayfix
      endif

      ! write weather data
      write(fin,'(/,"* Weather data")')
      write(fin,'("* Minimum temperatures of last week")')
      write(fin,'(" atmin7 = ", 7e14.5)') atmin7(1:7)

      if (flCropCalendar) then

        write(fin,'(/,"* Crop growth data of last crop: ", a)') trim(cropfil(icrop))
        write(fin,'(" daycrop = ",i4)') daycrop

        ! write preparation
        if (croptype(icrop) .le. 2) then

          if (flCropPrep) then
            swCropPrep = 1
          else
            swCropPrep = 0
          endif
          write(fin,'(" swCropPrep = ",i4)') swCropPrep

          ! write sowing
          if (flCropSow) then
            swCropSow = 1
          else
            swCropSow = 0
          endif
          write(fin,'(" swCropSow = ",i4)') swCropSow

          ! write germination
          if (flCropGerm) then
            swCropGerm = 1
          else
            swCropGerm = 0
          endif
          write(fin,'(" swCropGerm = ",i4)') swCropGerm

          write(fin,'(" PrepDelay = ",i4)') PrepDelay
          write(fin,'(" SowDelay = ",i4)') SowDelay
          write(fin,'(" tsumgerm = ",e15.8)') tsumgerm

        endif

        ! write crop emergence
        if (flCropEmergence) then
          swCropEmergence = 1
        else
          swCropEmergence = 0
        endif
        write(fin,'(" swCropEmergence = ",i4)') swCropEmergence

        ! write crop harvest
        if (flCropHarvest) then
          swCropHarvest = 1
        else
          swCropHarvest = 0
        endif
        write(fin,'(" swCropHarvest = ",i4)') swCropHarvest

        ! write crop growth data of fixed crop growth
        if (flCropEmergence) then

          write(fin,'(" rd = ",e15.8)') rd
          write(fin,'(" rdpot = ",e15.8)') rdpot
          if (swinter.eq.3) then
            write(fin,'(" sicact = ",e15.8)') sicact
          endif
          write(fin,'(" nofd = ",i4)') nofd
          write(fin,'(" dvs = ",e15.8)') dvs
          write(fin,'(" tsum = ",e15.8)') tsum

          ! write crop growth data of dynamic crop growth (WOFOST and GRASS)
          if (croptype(icrop) .ge. 2) then

            ! write crop data of wofost or simulated grass
            if (flanthesis) then
              swanthesis = 1
            else
              swanthesis = 0
            endif
            write(fin,'(" swanthesis = ",i4)') swanthesis
            write(fin,'(" ilvold = ",i4)') ilvold
            write(fin,'(" ilvoldpot = ",i4)') ilvoldpot
            write(fin,'(" wrt = ",e15.8)') wrt
            write(fin,'(" wrtpot = ",e15.8)') wrtpot
            write(fin,'(" tadw = ",e15.8)') tadw
            write(fin,'(" tadwpot = ",e15.8)') tadwpot
            write(fin,'(" wst = ",e15.8)') wst
            write(fin,'(" wstpot = ",e15.8)') wstpot
            write(fin,'(" wso = ",e15.8)') wso
            write(fin,'(" wsopot = ",e15.8)') wsopot
            write(fin,'(" wlv = ",e15.8)') wlv
            write(fin,'(" wlvpot = ",e15.8)') wlvpot
            write(fin,'(" laiexp = ",e15.8)') laiexp
            write(fin,'(" laiexppot = ",e15.8)') laiexppot
            write(fin,'(" glaiex = ",e15.8)') glaiex
            write(fin,'(" glaiexpot = ",e15.8)') glaiexpot
            write(fin,'(" lai = ",e15.8)') lai
            write(fin,'(" laipot = ",e15.8)') laipot
            write(fin,'(" laimax = ",e15.8)') laimax
            write(fin,'(" dwrt = ",e15.8)') dwrt
            write(fin,'(" dwrtpot = ",e15.8)') dwrtpot
            write(fin,'(" dwlv = ",e15.8)') dwlv
            write(fin,'(" dwlvpot = ",e15.8)') dwlvpot
            write(fin,'(" dwst = ",e15.8)') dwst
            write(fin,'(" dwstpot = ",e15.8)') dwstpot
            write(fin,'(" dwlvSoil = ",e15.8)') dwlvSoil
            write(fin,'(" dwlvCrop = ",e15.8)') dwlvCrop
            write(fin,'(" gasst = ",e15.8)') gasst
            write(fin,'(" gasstpot = ",e15.8)') gasstpot
            write(fin,'(" mrest = ",e15.8)') mrest
            write(fin,'(" mrestpot = ",e15.8)') mrestpot
            write(fin,'(" cwdm = ",e15.8)') cwdm
            write(fin,'(" cwdmpot = ",e15.8)') cwdmpot
            write(fin,'(" iharvest = ",i4)') iharvest
            write(fin,'(" rid = ",e15.8)') rid
            write(fin,'(" idregr = ",i4)') idregr
            write(fin,'(" idregrpot = ",i4)') idregrpot
            write(fin,'(" tagp = ",e15.8)') tagp
            write(fin,'(" tagppot = ",e15.8)') tagppot
            write(fin,'(" tagpt = ",e15.8)') tagpt
            write(fin,'(" tagptpot = ",e15.8)') tagptpot
            write(fin,'(" cuptgraz = ",e15.8)') cuptgraz
            write(fin,'(" cuptgrazpot = ",e15.8)') cuptgrazpot
            write(fin,'(" daygrowth = ",i4)') daygrowth
            write(fin,'(" daygrowthpot = ",i4)') daygrowthpot
            write(fin,'(" iseqgm = ",i4)') iseqgm
            write(fin,'(" iseqgmpot = ",i4)') iseqgmpot
            if (flgrazing) then
              swgrazing = 1
            else
              swgrazing = 0
            endif
            write(fin,'(" swgrazing = ",i4)')    swgrazing
            if (flgrazingpot) then
              swgrazingpot = 1
            else
              swgrazingpot = 0
            endif
            write(fin,'(" swgrazingpot = ",i4)')    swgrazingpot
            write(fin,'(" idaysgraz = ",i4)')    idaysgraz
            write(fin,'(" idaysgrazpot = ",i4)')    idaysgrazpot

            write(fin,'("")')
            write(fin,'(" day            sla         slapot          lvage       lvagepot             lv          lvpot")')
            count = max(ilvold,ilvoldpot)
            do i = 1,count
              write (fin,'(i4,6e15.8)')  i,sla(i),slapot(i),lvage(i),lvagepot(i),lv(i),lvpot(i)
            enddo
          endif
        endif
      endif

      ! end write crop growth data
      close (fin)

      return
      end
! ----------------------------------------------------------------------
      subroutine OutCropFixed(task)
! ----------------------------------------------------------------------
!     Date               : Aug 2004
!     Purpose            : open and write fixed crop output files
! ----------------------------------------------------------------------
      use variables, only: date,t,daycrop,dvs,tsum,lai,cf,rd,crp,ch
      implicit none

! --- global variables ------------------
      integer task

! --- local
      character(len=1) comma
! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)

! --- write header of new crop ----------------------------------------------
      write (crp,100)
 100  format ('*',/,                                                    &
     & '*             day     day      -    grC       -      -      cm',&
     & '       -       cm    cm ',                                      &
     & '    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha',        &
     & '    kg/ha    kg/ha    kg/ha    kg/ha',                          &
     & '    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha',        &
     & '     kg/ha     kg/ha     kg/ha     kg/ha     kg/ha     kg/ha',/ &
     & '      Date, Daynr, Daycrp,   DVS,  TSUM, LAIpot,    LAI,',      &
     & ' Height,CrpFac,RootdPot, Rootd,    PWLV,     WLV,    ',         &
     & 'PWST,     WST,    PWRT,     WRT,   CPWDM,    CWDM,   CPWSO,',   &
     & '    CWSO,PGRASSDM, GRASSDM,  PMOWDM,   MOWDM, PGRAZDM,  GRAZDM',&
     & ', DWLVCROP, DWLVSOIL,     DWST,     DWRT,     DWSO,HarLosOrm')

      return

      case (2)

! --- write actual data ------------------------------------------------------

! --- write output record
      write (crp,200) date,comma,nint(t),comma,daycrop,comma,dvs,comma, &
     & tsum,comma,"       ",comma,lai,comma,ch,comma,cf,                &
     & comma,"       ",comma,nint(rd),                                  &
     & comma,comma,comma,comma,comma,comma,comma,comma,comma,comma,     &
     & comma,comma,comma,comma,comma,comma,                             &
     & comma,comma,comma,comma,comma,comma
 200  format (a11,a1,i5,a1,i7,a1,f6.2,a1,f6.0,a1,a,2(a1,f7.2),(a1,f6.2),&
     &      (a1,a),(a1,i7),16(a1,'        ') ,                          &
     &      6(a1,'         ') )

      case default
         call fatalerr ('OutCropFixed', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine OutWofost(task)
! ----------------------------------------------------------------------
!     UpDate             : Aug 2014
!     Date               : Oct 2004
!     Purpose            : Write detailed crop growth output files
! ----------------------------------------------------------------------
      use variables, only: date,daycrop,crp,t,dvs,tsum,laipot,lai,rdpot,rd,ch,cf,cwdmpot,cwdm,wsopot,wso,wlvpot,wlv,wstpot,   &
                           wst,wrtpot,wrt,dwlvCrop,dwlvSoil,dwst,dwrt,dwso,HarLosOrm_tot,swbulb,wblpot,wbl,dwblpot,dwbl,plwt
      implicit none

! --- global variables ------------------
      integer task

! --- local variables ------------------
      character(len=1) comma
! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)

! --- write header of new crop ----------------------------------------------

      if(swbulb.eq.0) then
         write (crp,100)
      else if(swbulb.eq.1) then
         write (crp,200)
      endif
 100  format ('*',/,                                                    &
     & '*             day     day      -    grC       -      -      cm',&
     & '       -       cm    cm ',                                      &
     & '    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha',        &
     & '    kg/ha    kg/ha    kg/ha    kg/ha',                          &
     & '    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha',        &
     & '     kg/ha     kg/ha     kg/ha     kg/ha     kg/ha     kg/ha',/ &
     & '       Date,Daynr, Daycrp,   DVS,  TSUM, LAIpot,    LAI,',      &
     & ' Height,CrpFac,RootdPot, Rootd,    PWLV,     WLV,    ',         &
     & 'PWST,     WST,    PWRT,     WRT,   CPWDM,    CWDM,   CPWSO,',   &
     & '    CWSO,PGRASSDM, GRASSDM,  PMOWDM,   MOWDM, PGRAZDM,  GRAZDM',&
     & ', DWLVCROP, DWLVSOIL,     DWST,     DWRT,     DWSO,HarLosOrm')
 200   format ('*',/,                                                   &
     & '*             day     day      -    grC       -      -      cm',&
     & '       -       cm    cm ',                                      &
     & '    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha',        &
     & '    kg/ha    kg/ha    kg/ha    kg/ha',                          &
     & '    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha',        &
     & '     kg/ha     kg/ha     kg/ha     kg/ha     kg/ha     kg/ha',  &
     &   '       -     kg/ha     kg/ha     kg/ha     kg/ha     kg/ha',/ &
     & '       Date,Daynr, Daycrp,   DVS,  TSUM, LAIpot,    LAI,',      &
     & ' Height,CrpFac,RootdPot, Rootd,    PWLV,     WLV,    ',         &
     & 'PWST,     WST,    PWRT,     WRT,   CPWDM,    CWDM,   CPWSO,',   &
     & '    CWSO,PGRASSDM, GRASSDM,  PMOWDM,   MOWDM, PGRAZDM,  GRAZDM',&
     & ', DWLVCROP, DWLVSOIL,     DWST,     DWRT,     DWSO,HarLosOrm',  &
     & ', swbulb,   wblpot,      wbl,  dwblpot,     dwbl,     plwt')

      return

      case (2)

! --- write actual data ------------------------------------------------------

      if(swbulb.eq.0) then
         write (crp,300) date,comma,nint(t),comma,daycrop,comma,dvs,    &
     &    comma,tsum,comma,laipot,comma,lai,comma,ch,comma,cf,comma,    &
     &    nint(rdpot),comma,nint(rd),comma,nint(wlvpot),comma,nint(wlv),&
     &    comma,nint(wstpot),comma,nint(wst),comma,nint(wrtpot),        &
     &    comma,nint(wrt),comma,nint(cwdmpot),comma,nint(cwdm),comma,   &
     &    nint(wsopot),comma,nint(wso),comma,comma,comma,comma,comma,   &
     &    comma,comma,dwlvCrop,comma,dwlvSoil,comma,dwst,comma,dwrt,    &
     &    comma,dwso,comma,HarLosOrm_tot
      elseif(swbulb.eq.1) then
         write (crp,400) date,comma,nint(t),comma,daycrop,comma,dvs,    &
     &    comma,tsum,comma,laipot,comma,lai,comma,ch,comma,cf,comma,    &
     &    nint(rdpot),comma,nint(rd),comma,nint(wlvpot),comma,nint(wlv),&
     &    comma,nint(wstpot),comma,nint(wst),comma,nint(wrtpot),        &
     &    comma,nint(wrt),comma,nint(cwdmpot),comma,nint(cwdm),comma,   &
     &    nint(wsopot),comma,nint(wso),comma,comma,comma,comma,comma,   &
     &    comma,comma,dwlvCrop,comma,dwlvSoil,comma,dwst,comma,dwrt,    &
     &    comma,dwso,comma,HarLosOrm_tot,comma,swbulb,comma,wblpot,     &
     &    comma,wbl,comma,dwblpot,comma,dwbl,comma,plwt
      endif
 300  format (a11,a1,i5,a1,i7,a1,f6.2,a1,f6.0,3(a1,f7.2),(a1,f6.2),     &
     &        2(a1,i7),  10(a1,i8), 6(a1,'        ') ,                  &
     &        6(a1,f9.2))
 400  format (a11,a1,i5,a1,i7,a1,f6.2,a1,f6.0,3(a1,f7.2),(a1,f6.2),     &
     &        2(a1,i7),  10(a1,i8), 6(a1,'        ') ,                  &
     &        6(a1,f9.2), a1,i7, 5(a1,f9.2) )

      case default
         call fatalerr ('OutWofost', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine OutGrass(task)
! ----------------------------------------------------------------------
!     UpDate             : Aug 2014
!     Date               : Oct 2004
!     Purpose            : Write detailed grass simulation output files
! ----------------------------------------------------------------------
      use variables, only: date,daycrop,crp,t,dvs,tsum,laipot,lai,rdpot,rd,ch,cf,tagppot,tagp,tagptpot,tagpt,          &
                           wlvpot,wlv,wstpot,wst,wrtpot,wrt,cuptgraz,cuptgrazpot
      implicit none

! --- global variables ------------------
      integer task

! --- local variables ------------------
      character(len=1) comma
! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)

! --- write header of new crop ----------------------------------------------

      write (crp,100)
 100  format ('*',/,                                                    &
     & '*             day     day      -    grC       -      -      cm',&
     & '       -       cm    cm ',                                      &
     & '    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha',        &
     & '    kg/ha    kg/ha    kg/ha    kg/ha',                          &
     & '    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha    kg/ha',        &
     & '     kg/ha     kg/ha     kg/ha     kg/ha     kg/ha     kg/ha',/ &
     & '      Date, Daynr, Daycrp,   DVS,  TSUM, LAIpot,    LAI,',      &
     & ' Height,CrpFac,RootdPot, Rootd,    PWLV,     WLV,    ',         &
     & 'PWST,     WST,    PWRT,     WRT,   CPWDM,    CWDM,   CPWSO,',   &
     & '    CWSO,PGRASSDM, GRASSDM,  PMOWDM,   MOWDM, PGRAZDM,  GRAZDM',&
     & ', DWLVCROP, DWLVSOIL,     DWST,     DWRT,     DWSO,HarLosOrm')

      return

      case (2)

! --- write actual data ------------------------------------------------------

! --- write output record
      write (crp,200) date,comma,nint(t),comma,daycrop,comma,dvs,comma, &
     & tsum,comma,laipot,comma,lai,comma,ch,comma,cf,                   &
     & comma,nint(rdpot),comma,nint(rd),                                &
     & comma,nint(wlvpot),comma,nint(wlv),comma,nint(wstpot),           &
     & comma,nint(wst),comma,nint(wrtpot),comma,nint(wrt),              &
     & comma,comma,comma,comma,comma,nint(tagppot),                     &
     & comma,nint(tagp),comma,nint(tagptpot),comma,nint(tagpt),         &
     & comma,nint(cuptgrazpot),comma,nint(cuptgraz),                    &
     & comma,comma,comma,comma,comma,comma
 200  format (a11,a1,i5,a1,i7,a1,f6.2,a1,f6.0,3(a1,f7.2),(a1,f6.2),     &
     &  2(a1,i7), 6(a1,i8), 4(a1,'        '), 6(a1,i8),                 &
     &  6(a1,'         '))

      case default
         call fatalerr ('OutGrass', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine SoluteOutput(task)
! ----------------------------------------------------------------------
!     Date               : November 2004
!     Purpose            : open and write solute output files
! ----------------------------------------------------------------------

      use Variables
      implicit none

      integer task

      select case (task)
      case (1)

! === open output files and write headers ===============================

! --  sba file
      if (swsba .eq. 1) call outsba (1)

      return

      case (2)

! === write actual data ===============================

! --  sba file
      if (swsba .eq. 1) call outsba (2)

      return

      case (3)

! === close output files ===========================

! --- close sba file
      if (swsba .eq. 1) close (sba)

      case default
         call fatalerr ('SoluteOutput', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outsba (task)
! ----------------------------------------------------------------------
!     date               : July 2002
!     purpose            : output of salt balance
! ---------------------------------------------------------------------
      use variables, only: sba,daynr,daycum,sampro,sqbot,project,sqdra,solbal,dectot,rottot,sqprec,                    &
                           date,sqirrig,outfil,pathwork,flheader
      implicit none

! --- global
      integer   task

! --- local variables ------------------
      integer   getun
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
! ----------------------------------------------------------------------
      comma = ','


      select case (task)
      case (1)

! --- open output file -------------------------------------------------
      filnam = trim(pathwork)//trim(outfil)//'.sba'
      sba = getun (20,90)
      call fopens(sba,filnam,'new','del')
      filtext = 'cumulative solute balance components (mg/cm2)'
      call writehead (sba,1,filnam,filtext,project)

! --- write header in sba file
      write (sba,10)
10    format ('*',/,                                                    &
     & '       Date, Day,  Dcum,      Flux top,   Root uptake, ',       &
     & 'Decomposition,      Drainage,   Flux bottom,       Storage,',   &
     & '       Balance')

      return

      case (2)

! === write actual data =================================================

! --- write header in case of new balance period
      if (flheader) write (sba,10)

! --- write output solute balance components ----------------------------

      write(sba,15) date,comma,daynr,comma,daycum,comma,(sqprec+sqirrig)&
     & ,comma,rottot,comma,dectot,comma,sqdra,comma,sqbot,comma,        &
     & sampro,comma,solbal

 15   format(a11,a1,i4,a1,i6,6(a1,e14.5),a1,e14.2)

      case default
         call fatalerr ('outsba', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine AgeTracerOutput(task)
! ----------------------------------------------------------------------
!     Date               : October 2010
!     Purpose            : open and write Groundwater Ageing output files
! ----------------------------------------------------------------------
      implicit none

! --- global variables ------------------
      integer task
! --- local variables ------------------
      integer agep,agee,ageq

      save   agee,agep,ageq

      select case (task)
      case (1)

! === open output files and write headers ===============================

! --  age files
      call outage (1,agep,agee,ageq)

      return

      case (2)

! === write actual data ===============================

! --  age files
      call outage (2,agep,agee,ageq)

      return

      case (3)

! === close output files ===========================

! --- close sba file
      close (agep)
      close (agee)

      case default
         call fatalerr ('AgeTracerOutput', 'Illegal value for Task')
      end select

      return
      end
! ----------------------------------------------------------------------
      subroutine outage(task,agep,agee,ageq)
! ----------------------------------------------------------------------
!     date               : October 2010
!     purpose            : output of groundwater age
! ---------------------------------------------------------------------
      use variables, only: daynr,daycum,date,outper,project,nrlevs,outfil,pathwork,numnod,z,cml,            &
                           AgeGwl1m,icAgeBot,icAgeDra,icAgeRot,icAgeSur,inqdra
      implicit none
      include 'arrays.fi'

! --- global
      integer   agep,agee,ageq,task

! --- local variables ------------------
      integer   getun,reclngth,node,level
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
      real(8)   iqdrainout(madr)   ! Cumulative (over 1 output timestep) drainage flux (L) for each drainage level
! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)

! --- open output files -------------------------------------------------
!     age of groundwater as profile
      filnam = trim(pathwork)//trim(outfil)//'.ageProfile.csv'
      agep = getun (20,90)
      reclngth = 50 + 12*numnod
      open(unit=agep,file=filnam,status='unknown',recl=reclngth)
      filtext = 'Groundwater age profiles (all age-values in days)'
      call writehead (agep,1,filnam,filtext,project)
!     age of groundwater in effluents: drains, transpiration, leaching, runoff
      filnam = trim(pathwork)//trim(outfil)//'.ageEffluent.csv'
      agee = getun (20,90)
      call fopens(agee,filnam,'new','del')
      filtext = 'Groundwater age effluent (all age-values in days)'
      call writehead (agee,1,filnam,filtext,project)
!     effluent drain water fluxes
      filnam = trim(pathwork)//trim(outfil)//'.ageEffluentqDrain.csv'
      ageq = getun (20,90)
      call fopens(ageq,filnam,'new','del')
      filtext = 'Drain water effluent (mm/day)'
      call writehead (ageq,1,filnam,filtext,project)

! --- write headers of files
      write (agep,9) (z(node),node=1,numnod)
      if (numnod.le.9) then
         write (agep,10) (node,node=1,numnod)
      else
         write (agep,11) (node,node=1,9), (node,node=10,numnod)
      endif
  9   format('*',t8,' NodeDepth (cm) =,,',18(',',f7.2) ,997(',',f8.2) )
 10   format('*',/, t8,'Date,Day,Daycum',   9(',Node',i3.3) )
 11   format('*',/, t8,'Date,Day,Daycum',   9(',Node',i3.3),            &
     &                                   1015(',Node',i4.4) )

      write (agee,'(3a)') ' Date,daynr,daycum,AgeGwl1m,AgeBottom,',     &
     &'AgeRootUpt,AgeRunoff,AgeDrainSys1,AgeDrainSys2,AgeDrainSys3,',   &
     &'AgeDrainSys4,AgeDrainSys5'

      write (ageq,'(3a)') ' Date,daynr,daycum,',                        &
     &'qDrainSys1,qDrainSys2,qDrainSys3,qDrainSys4,qDrainSys5'

      return

      case (2)

! === write actual data =================================================

!     age of groundwater as profile
      write(agep,15) date,comma,daynr,comma,daycum,                     &
     &               (comma,cml(node),node=1,numnod)
 15   format(a11,a1,i4,a1,i6,1p,1024(a1,e10.3))

!     age of groundwater in effluents: drains, transpiration, leaching, runoff
!     and age (d) of groundwater in upper 1 meter of saturated zone
      write(agee,16) date,comma,daynr,comma,daycum,comma,AgeGwl1m,comma,&
     &               icAgeBot/outper,comma,icAgeRot/outper,comma,       &
     &               icAgeSur/outper,                                   &
     &               (comma,icAgeDra(level)/outper,level=1,nrlevs)
 16   format(a11,a1,i4,a1,i6,1p,9(a1,e10.3))

!     qdrain discharge-effluent (without infiltration!)
      do level = 1,nrlevs
        iqdrainout(level) = 0.0d0
        do node = 1,numnod
          if (inqdra(level,node).gt.0.0d0) then
           iqdrainout(level) = iqdrainout(level) + inqdra(level,node)
          endif
        enddo
      enddo
      write(ageq,16) date,comma,daynr,comma,daycum,                     &
     &               (comma,iqdrainout(level),level=1,nrlevs)

      case default
         call fatalerr ('outage', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine TemperatureOutput(task)
! ----------------------------------------------------------------------
!     Date               : November 2004
!     Purpose            : open and write soil temperature output files
! ----------------------------------------------------------------------

! --- global variables ------------------
      use Variables

      implicit none
      integer task


      select case (task)
      case (1)

! === open output files and write headers ===============================

! --  heat params file
      if (swini .eq. 1 .and. swcalt.eq.2)  call outheapar()

! --  tem file
      if (swtem .eq. 1) call outtem (task)

      return

      case (2)

! === write actual data ===============================

! --  tem file
      if (swtem .eq. 1) call outtem (task)

      return

      case (3)

! === close output files ===========================

! --- close tem file
      if (swtem .eq. 1) close (tem)

      case default
         call fatalerr ('TemperatureOutput', 'Illegal value for Task')
      end select

      return
      end



! ----------------------------------------------------------------------
      subroutine outheapar ()
! ----------------------------------------------------------------------
!     date               : February 2005
!     purpose            : Output of soil heat conductivity and capacity
! ---------------------------------------------------------------------
! --- global
      use variables

! --- local variables ------------------
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
      integer   getun,hea,lay,node,j
      real(8)   heacnd(macp),thetadum(numnod)
! ---------------------------------------------------------------------
      comma = ','

! === open output file =================================================
      filnam = trim(pathwork)//'heatparam.csv'
      hea = getun (20,90)
      call fopens(hea,filnam,'new','del')
      filtext = 'soil heat conductivity and capacity'
      call writehead (hea,1,filnam,filtext,project)

! --- write header
      write(hea,'(a)') 'layer, theta,heacap(J/cm3/K),heacnd(J/cm/K/d)'

! --- write thermal properties: heat capacity and thermal conductivity
      do node = 1,numnod
        thetadum(node) = thetas(node)
      enddo
      do lay = 1,numlay
!         find first Node of the Layer
          Node= 1
          do while(Layer(Node).ne.lay)
            Node= Node + 1
          enddo
          do j = 1,21
            thetadum(node) = thetar(node) + dble(j-1) *               &
     &                    (thetas(node)-thetar(node)) / 20.0d0
            call devries (thetadum,heacap,heacnd)
            write(hea,22) lay, comma, thetadum(node), comma,            &
     &                 heacap(node), comma, heacnd(node)
22          format(i4,a1,f8.5,2(a1,e14.5))
          enddo
      enddo

! --- close heat params file
      close (hea)

      return
      end

! ----------------------------------------------------------------------
      subroutine outsoilphys ()
! ----------------------------------------------------------------------
!     date               : March 2008
!     purpose            : Output of soil physical parameters
!        cofgen(1,node) = ores
!        cofgen(2,node) = osat
!        cofgen(3,node) = ksatfit
!        cofgen(4,node) = alfa
!        cofgen(5,node) = lexp
!        cofgen(6,node) = npar
!        cofgen(7,node) = 1.d0 - (1.d0 / npar)
!        cofgen(8,node) = dummy
!        cofgen(9,node) = h_enpr
!        cofgen(10,node)= ksatexm
!        cofgen(11,node)= relsatthr
!        cofgen(12,node)= ksatthr
!        cofgen(13,node)= alfa_2
!        cofgen(14,node)= npar_2
!        cofgen(15,node)= mpar_2 (=1-1/npar_2)
!        cofgen(16,node)= omega_1
!        cofgen(17,node)= omega_2 (omega_2=1-omega_1)
!        cofgen(18,node)= h0
!        cofgen(19,node)= ha
!        cofgen(20,node)= apar
!        cofgen(21,node)= omega_K
! ---------------------------------------------------------------------
! --- global
      use variables
      use doln
      implicit none

! --- local variables ------------------
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
      integer   getun,soi,lay,hh,numhead,node
      real(8)   watcon,moiscap,hconduc,relsat, dhx!, dkx, dhconduc
      real(8)   dimocax,kx,thetax, FrArMtrx1,rfcpx
      parameter (numhead = 341)
      real(8)   hx(numhead)
! ---------------------------------------------------------------------
      comma = ','
      FrArMtrx1 = 1.d0
      rfcpx = 1.d0

! === open output file =================================================
      filnam = trim(pathwork)//'SoilPhysParam.csv'
      soi = getun (20,90)
      call fopens(soi,filnam,'new','del')
      filtext = 'soil physical parameters '
      call writehead (soi,1,filnam,filtext,project)

! --- write header
      write(soi,'(a)')                                                  &
     &    'layer,head(cm),theta(cm3.cm3),C(cm-1),RelSat(-),k(cm.d-1)'!,dK/dh(d-1)'

! --- determine range of heads for which soil physical params are calculated
      if(swsophy.eq.0)then
         hx(1) = 0.0d0
         hx(2) = -1.0d-7
         hx(numhead) = -1.0d07
         do hh = 3,numhead-1
           hx(hh) = hx(hh-1)*1.1d0
         enddo
      end if

! --- soil physical params: calculate and write
      do lay = 1,numlay
!        find first Node of the Layer
         Node = nod1lay(lay)
         if(swsophy.eq.1)then
            hx(1) = 0.0d0
            hx(2) = -1.0d-3 ! -1.0d-1

!## MH start
            if (do_ln_trans) then
               hx(numhead) = -(dexp(-sptab(1,node,1)) - 1.0d0)
            else
               hx(numhead) = sptab(1,node,1)
            end if
!            dhx = 10.d0**((log10(-sptab(1,node,1))+1.d0)/(numhead-2))
            dhx = 10.d0**((log10(-hx(numhead))+3.d0)/(numhead-2))
!## MH end

            do hh = 3,numhead
              hx(hh) = hx(hh-1)* dhx
            enddo
         end if

         do hh = 2, numhead
            thetax = watcon(node,hx(hh))
            dimocax = moiscap(node,hx(hh))
            kx = hconduc (node,hx(hh),thetax,rfcpx)
            !dkx = dhconduc (node,hx(hh),thetax,dimocax,rfcpx)
!           in case of static macropores FrArMtrx < 1
            if(FlMacropore)  kx = FrArMtrx1 * kx
!           write output
            relsat = (thetax-cofgen(1,Node)) /                          &
     &                 (cofgen(2,Node)-cofgen(1,Node))
            write (soi,22) lay, comma, hx(hh), comma, thetax, comma,    &
     &                     dimocax, comma, relsat, comma, kx!, comma, dkx
 22         format(i10,6(a,1pe15.7))
          enddo
      enddo

! --- close soil physical params file
      close (soi)

      return
      end

! ----------------------------------------------------------------------
      subroutine outtem (task)

! ----------------------------------------------------------------------
!     date               : November 2004
!     purpose            : Output of soil temperatures
! ---------------------------------------------------------------------
      use variables, only: numnod,date,daynr,tem,daycum,tav,tebot,tsoil,tetop,outfil,pathwork,flheader,project
      implicit none

! --- global
      integer   task

! --- local variables ------------------
      integer      getun,i, reclngth
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)

! === open output file =================================================
      filnam = trim(pathwork)//trim(outfil)//'.tem'
      tem = getun (20,90)
!      reclngth = 36 + 7*numnod
      reclngth = 50 + 7*numnod
      open(unit=tem,file=filnam,status='unknown',recl=reclngth)
      filtext = 'soil temperature profiles (oC)'
      call writehead (tem,1,filnam,filtext,project)

! --- write header
      if (numnod.le.9) then
         write (tem,10) (i,i=1,numnod)
      else
         write (tem,11) (i,i=1,9), (i,i=10,numnod)
      endif

 10   format('*',/,                                                     &
     & t8,'Date,Day,Daycum,   Tav, Tetop',9(',    T',i1),               &
     &', TeBot')

 11   format('*',/,                                                     &
     & t8,'Date,Day,Daycum,   Tav, Tetop',9(',    T',i1),               &
     &1024(',   T',i2),', TeBot')

      write (tem,'(a11,a1,i3,a1,i6,1024(a1,f6.1:))') '    Initial'      &
     &      ,comma,daynr,comma,daycum,comma,tav,comma,tetop,            &
     &      (comma,tsoil(i),i=1,numnod),comma,tebot

      return

      case (2)

! === write actual soil temperature data ================================

! --- write header in case of new balance period
      if (flheader) write (tem,10)

! --- write soil temperature profile
!     PWB: idem
!      write (tem,'(a11,a1,i3,a1,i6,<numnod+3>(a1,f6.1:))') date
!     &      ,comma,daynr,comma,daycum,comma,tav,comma,tetop,
!     &      (comma,tsoil(i),i=1,numnod),comma,tebot
      write (tem,'(a11,a1,i3,a1,i6,1024(a1,f6.1:))') date               &
     &      ,comma,daynr,comma,daycum,comma,tav,comma,tetop,            &
     &      (comma,tsoil(i),i=1,numnod),comma,tebot

      case default
         call fatalerr ('outtem', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine IrrigationOutput(task)
! ----------------------------------------------------------------------
!     Date               : November 2004
!     Purpose            : open and write irrigation output files
! ----------------------------------------------------------------------

      use variables
      implicit none

! --- global variables ------------------
      integer   task
! --- local variables ------------------
      integer   getun
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
! ----------------------------------------------------------------------
      comma = ','

      ! return in case of no output for irrigation
      if (swirg .eq. 0) return
      
      select case (task)
      case (1)

! === open output file and write header ===============================

      filnam = trim(pathwork)//trim(outfil)//'.irg'
      irg = getun (20,90)
      call fopens(irg,filnam,'new','del')
      filtext = 'irrigation data'
      call writehead (irg,1,filnam,filtext,project)

      write (irg,100)
 100  format ('*',/,                                                    &
     &  '*      Crop name        Date   Day   Day  Irrigation',         &
     &  '     Solute',/,                                                &
     &  '*                             Crop   Cum          cm',         &
     &  '     mg/cm3',/,                                                &
     &  '*<=============><==========><====><====><==========>           &
     &<=========>',/,                                                   &
     &  '       Crop name,       Date,  Day,  Day, Irrigation,',        &
     &  '    Solute')

      flIrg1Start = .false.

      return

      case (2)

! === write actual data ===============================

! --- write header in case of new balance period
      if (flheadirg) then
        write (irg,100)
       flheadirg = .false.
      endif

! --- write output record irrigation
      if (irrigevent.eq.1 .and. flCropCalendar                          &
     &                                   .and. .not. flCropHarvest) then
        write (irg,'(a16,a1,a11,a1,i4,a1,i6,a1,2x,f8.1,a1,e11.3)')      &
     &     cropfil(icrop),comma,date,comma,daycrop,                    &
     &     comma,daycum,comma,gird,comma,cirr
      elseif (irrigevent.eq.1) then
        write (irg,'(t1,a,t18,a11,a1,i4,a1,i6,a1,2x,f8.1,a1,e11.3)')    &
     &     'bare soil,',date,comma,daycrop,                             &
     &     comma,daycum,comma,gird,comma,cirr

      elseif (irrigevent.eq.2) then
        write (irg,'(a16,a1,a11,a1,i4,a1,i6,a1,2x,f8.1,a1,e11.3)')      &
     &     cropfil(icrop),comma,date,comma,daycrop,                    &
     &     comma,daycum,comma,gird,comma,cirr
      endif

      flIrrigationOutput = .false.

      return

      case (3)

! === close output files ===========================

! --- close irg file
      close (irg)

      case default
         call fatalerr ('IrrigationOutput', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine SnowOutput(task)
! ----------------------------------------------------------------------
!     Date               : December 2004
!     Purpose            : open and write snow pack data
! ----------------------------------------------------------------------

      use variables
      implicit none

! --- global variables ------------------
      integer task
! --- local variables ------------------
      integer   getun
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)

! === open output file =================================================
      filnam = trim(pathwork)//trim(outfil)//'.snw'
      snw = getun (20,90)
      call fopens(snw,filnam,'new','del')
      filtext = 'snow pack output data (cm/period)'
      call writehead (snw,1,filnam,filtext,project)

! --- write header
      write (snw,10)
 10   format ('*',/,                                                    &
     &   '    date,      dcum,  rainfall,  snowfall,snowstorage, ',     &
     &   'meltflux,sublimation')

      return

      case (2)

! === write actual soil temperature data ================================

! --- write header in case of new balance period
      if (flheader) write (snw,10)

! --- write actual data
      write (snw,20) date,comma,daycum,comma,snrai,comma,gsnow,         &
     &               comma,ssnow,comma,melt,comma,subl
20    format (a11,a1,i6,1x,5(a1,f10.4))

      return

      case (3)

! === close output file ===========================

! --- close snw file
      close (snw)

      case default
         call fatalerr ('SnowOutput', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outblc(task,CQMpExfMtx,CQMpInfMtx)
! ----------------------------------------------------------------------
!     date               : December 2004
!     purpose            : Write detailed overview of water balance
! ---------------------------------------------------------------------
      use variables, only: zbotcp,blc,outfil,pathwork,cgrai,cnrai,cgird,cnird,cqrot,cevap,volact,volini,FlMacropore,nrlevs,swirfix,       &
                           schedule,numnod,snowinco,ssnow,cgsnow,cmelt,caintc,csnrai,cqprai,ioutdat,t1900,outdat,tstart,project,pond,pondini,     &
                           cqdrainin,cqdrainout,cinund,crunoff,cqtdo,cqtup,cqbotdo,cqbotup,csubl,crunon,IcTopMp,CQMpInTopVrtDm1,CQMpInTopVrtDm2,                  &
                           CQMpInTopLatDm1,CQMpInTopLatDm2,cqssdi
      implicit none

!     global
      integer   task
      real(8)   CQMpInfMtx, CQMpExfMtx

!     local
      integer   getun,level

      real(8)   CQMpInTop,soilout
      real(8)   plantin,plantout,snowin,snowout,pondin,pondout,soilin

      character(len=300) filnam
      character(len=11)  datbegin,datend
      character(len=80)  filtext
! ----------------------------------------------------------------------

      select case (task)
      case (1)
! --- open output file
      filnam = trim(pathwork)//trim(outfil)//'.blc'
      blc = getun (50,90)
      call fopens(blc,filnam,'new','del')
      filtext='detailed overview of water balance components (cm)'
      call writehead (blc,1,filnam,filtext,project)

      return

      case (2)

! --- begin date of balance period
      if (ioutdat .eq. 2) then
        call dtdpst ('year-month-day',tstart+0.1d0,datbegin)
      else
        call dtdpst ('year-month-day',                                &
     &                                outdat(ioutdat-2)+1.1d0,datbegin)
      endif

! --- end date of balance period
      call dtdpst ('year-month-day',t1900-0.9d0,datend)

! --- initialize macropore variables
      if (.not.FlMacropore) then
         CQMpInfMtx = 0.0d0
         CQMpExfMtx = 0.0d0
         CQMpInTop = 0.0d0
      else
         if (IcTopMp.eq.1) then
            CQMpInTop = CQMpInTopVrtDm1 + CQMpInTopVrtDm2 +                &
     &               CQMpInTopLatDm1 + CQMpInTopLatDm2
         else
            CQMpInTop = 0.d0
         endif
      endif

! --- write output record
      write (blc,20) datbegin,datend
      write (blc,22) -zbotcp(numnod)
      write (blc,40)
      write (blc,41) snowinco,pondini,volini,ssnow,pond,volact
      write (blc,42) cgrai
      write (blc,44) csnrai,cqprai,cnrai
      if (swirfix.eq.1 .or. schedule.eq.1) then
        write (blc,46) cgird
        write (blc,67) cnird,cnird
      endif
      write (blc,45) caintc
      write (blc,43) cgsnow
      write (blc,47) cmelt,cmelt
      write (blc,68) csubl
      write (blc,48) cqssdi, cqrot
      write (blc,491) cevap
      write (blc,492) crunon,crunoff
      write (blc,50) cinund
      write (blc,52) cqtdo,cqtdo
      write (blc,53) cqtup,cqtup
      if (FlMacropore)  write(blc,54) CQMpInTop, CQMpInfMtx, CQMpExfMtx

      if (nrlevs .ge. 1) then
        write (blc,55)
        do level=1,nrlevs
          write (blc,56) level,cqdrainin(level),level, cqdrainout(level)
        enddo
      endif
      write (blc,61) cqbotup,cqbotdo

! --- sum
      plantin = cgrai+cgird
      snowin = snowinco+cgsnow+csnrai
      pondin = pondini+cqprai+cnird+cmelt+cinund+cqtup+crunon
      soilin = volini+cqtdo+cqbotup+CQMpInfMtx
      do level = 1,nrlevs
        soilin = soilin + cqdrainin(level)
      enddo
      plantout = cnrai+caintc+cnird
      snowout = cmelt+ssnow+csubl
      pondout = pond+crunoff+cqtdo+cevap+CQMpInTop
      soilout = volact+cqrot+cqtup+cqbotdo+CQMpExfMtx-cqssdi
      do level = 1,nrlevs
        soilout = soilout + cqdrainout(level)
      enddo

      write (blc,63) plantin,snowin,pondin,soilin,plantout,snowout,     &
     &       pondout,soilout
      write (blc,64) (ssnow-snowinco),(pond-pondini),(volact-volini)
      write (blc,65) (plantout-plantin),(snowout-snowin),               &
     &     (pondout-pondin),(soilout-soilin)

 20   format(/'Period',t20,':',t23,a11,' until  ',a11)
 22   format('Depth soil profile',t20,':',f8.2,' cm')
 40   FORMAT(49('='),'+',49('='),/,'INPUT',t50,'|',t52,'OUTPUT',/,t17,  &
     &'   PLANT    SNOW    POND    SOIL',t50,'|',t67,                   &
     &'   PLANT    SNOW    POND    SOIL',/,49('='),'+',49('='))
 41   FORMAT('Initially Present',t25,3f8.2,t50,'|',t52,                 &
     &       'Finally present',t75,3f8.2)
 42   FORMAT('Gross Rainfall',t17,f8.2,t50,'|')
 43   FORMAT('Snowfall',t25,f8.2,t50,'|')
 44   FORMAT('Nett Rainfall',t25,2f8.2,t50,'|',t52,'Nett Rainfall',     &
     &        t67,f8.2)
 45   FORMAT(t50,'|',t52,'Interception',t67,f8.2)
 46   FORMAT('Gross Irrigation',t17,f8.2,t50,'|')
 67   FORMAT('Nett Irrigation',t33,f8.2,t50,'|',t52,'Nett Irrigation',  &
     &        t67,f8.2)
 47   FORMAT('Snowmelt',t33,f8.2,t50,'|',t52,'Snowmelt',t75,f8.2)
 68   FORMAT(t50,'|',t52,'Sublimation',t75,f8.2)
 48   FORMAT('SSDI',t41,f8.2,t50,'|',t52,'Plant Transpiration',t91,f8.2)
 491  FORMAT(t50,'|',t52,'Soil Evaporation',t83,f8.2)
 492  FORMAT('Runon',t33,f8.2,t50,'|',t52,'Runoff',t83,f8.2)
 50   FORMAT('Inundation',t33,f8.2,t50,'|')
 52   FORMAT('Infiltr. Soil Surf.',t41,f8.2,t50,'|',t52,                &
     &       'Infiltr. Soil Surf.',t83,f8.2)
 53   FORMAT('Exfiltr. Soil Surf.',t33,f8.2,t50,'|',t52,                &
     &       'Exfiltr. Soil Surf.',t91,f8.2)
 54   FORMAT(t50,'|',t52,'Inflow top macroprs',t83,f8.2,/,              &
     &       'Infiltr. macropores',t41,f8.2,t50,'|',t52,                &
     &       'Exfiltr. macropores',t91,f8.2)
 55   FORMAT('Infiltr. subsurf.',t50,'|',t52,'Drainage')
 56   FORMAT('- system',i2,t41,f8.2,t50,'|',t52,'- system',i2,t91,f8.2)
 61   FORMAT('Upward seepage',t41,f8.2,t50,'|',t52,'Downward seepage',  &
     &        t91,f8.2,/,49('='),'+',49('='))
 63   FORMAT('Sum',t17,4f8.2,t50,'|',t52,'Sum',t67,4f8.2,/,49('='),'+', &
     &        49('='))
 64   FORMAT('Storage Change',t25,3f8.2)
 65   FORMAT('Balance Deviation',t18,f7.2,3f8.2,/,99('='),/)

      case default
         call fatalerr ('outbLC', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outafo (task)
! ----------------------------------------------------------------------
!     Date               : 29-jan-2003
!     Purpose            : ANIMO/PEARL output: formatted hydrological data
! ---------------------------------------------------------------------
      use variables, only: afo,outfil,pathwork,numnod,outper,period,pond,gwl,nrlevs,ievap,ipeva,iptra,    &
                           iruno,numlay,botcom,thetas,kdif,kdir,swafo,igrai,igird,inird,iintc,gc,lai,tav,tsoil,rd,cf,wbalance, &
                           project,tstart,tend,swdiscrvert,numnodnew,dznew,FlMacropore,Ssnow,igSnow,isnrai,iSubl,irunon,IcTopMp,     &
                           WalevDm1,VlMpDm1,WaSrDm1,VlMpDm2,WaSrDm2,IQInTopVrtDm1,IQInTopLatDm1,            &
                           IQInTopVrtDm2,IQInTopLatDm2,        &
                           CritDevMasBal,tcum,nod1lay,out_etr,out_hum,out_rad,out_tmn,out_tmx,out_wet,out_win, FlOpenFileDev
      implicit none
      include 'arrays.fi'

! -   global
      integer   task

! -   local
      integer   getun,ftype,node,bruny,eruny,brund,erund,lay,level
      integer   swop,botcomNew(maho)
      integer   datea(6)
      real(4)   fsec
      real(8)   timjan1,dval(maho),watcon
      real(8)   hNew(macp),inqNew(macp+1)
      real(8)   inqdraNew(Madr,macp),thetaNew(macp),inqrotNew(macp)
      real(8)   TsoilNew(0:macp),tsoili(0:macp)
      real(8)   DiPoCpNew(macp), IAvFrMpWlWtDm1New(macp)
      real(8)   IAvFrMpWlWtDm2New(macp), InQExcMtxDm1CpNew(macp)
      real(8)   InQExcMtxDm2CpNew(macp), InQOutDrRapCpNew(macp)
      real(8)   IThetaBegNew(MaCp)
      real(8)   VlMpStDm1New(macp), VlMpStDm2New(macp)
      real(8)   TopQprecDm1, TopQsoilDm1, TopQprecDm2, TopQsoilDm2

      character(len=80)  filtext
      character(len=300) filnam
      character(len=4)   afoext

!     save local variables
      save      brund, swop


! ----------------------------------------------------------------------
      select case (task)
      case (1)

! === initial output ====================================================

! -   verify reduced grid for output
      if(SwDiscrvert.eq.1) then
         call checkDiscrVert()
      endif

! -   set switch Swop
      Swop = 1
! -   BFO and CFO file for PEARL
      if (Swafo.ge.2 .and. FlMacropore) Swop = 2

! -   convert vertical discretization: initial values
      do node = 1,numnod
         tsoili(node) = tsoil(node)
      enddo
      call convertdiscrvert(1,swop,botcomnew,hnew,thetanew,inqnew,inqrotnew,inqdranew,ithetabegnew,tsoili,tsoilnew,   &
                            dipocpnew,iavfrmpwlwtdm1new,iavfrmpwlwtdm2new,inqexcmtxdm1cpnew,inqexcmtxdm2cpnew, &
                            inqoutdrrapcpnew,vlmpstdm1new,vlmpstdm2new)
! -   open file
      if (swafo.eq.1) then
         afoext = '.afo'
      elseif (swafo.eq.2) then
         afoext = '.bfo'
      elseif (swafo.eq.3) then
         afoext = '.cfo'
      endif
      filnam = trim(pathwork)//trim(outfil)//afoext
      afo = getun (20,90)
      call fopens(afo,filnam,'new','del')


! ---   write initial part
      if (swafo.ge.2) then
         ftype = 1
         if (swafo.eq.2) then
            filtext = 'formatted hydrological data; BFO version'
         else if (swafo.eq.3) then
            filtext = 'formatted hydrological data; CFO version'
         endif
         call writehead(afo,ftype,filnam,filtext,project)
      endif

      call dtdpar (tstart,datea,fsec)
      bruny = datea(1)
      datea(2) = 1
      datea(3) = 1
      fsec = 0.0
      call dtardp (datea,fsec,timjan1)
      brund = nint (tstart - timjan1 + 1.0d0)
      call dtdpar (tend,datea,fsec)
      eruny = datea(1)
      datea(2) = 1
      datea(3) = 1
      fsec = 0.0
      call dtardp (datea,fsec,timjan1)
      erund = nint (tend - timjan1 + 1.0d0)
      if (swafo.eq.1) then
         write (afo,4000) bruny,eruny,                                  &
     &                 real(brund-1),real(erund),real(period)
      else if (swafo.ge.2) then
         write (afo,4001) swop
         write (afo,4000) bruny,eruny,real(brund-1),real(erund)
      endif

      write (afo,4010) numnodNew,numlay,nrlevs
      write (afo,4010) (botcomNew(lay),lay=1,numlay)
!  -    ThetaS should be known for new soil layers; this works but is not very nice
      write (afo,4020) (real(thetas(botcom(lay))), lay=1,numlay)

      do lay = 1,numlay
!        find first Node of the Layer
         Node = nod1lay(lay)
         dval(lay) = watcon(node,-100.d0)
      enddo
      write (afo,4020) (real(dval(lay)),lay=1,numlay)

      do lay = 1,numlay
!        find first Node of the Layer
         Node = nod1lay(lay)
         dval(lay) = watcon(node,-15849.0d0)
      enddo
      write (afo,4020) (real(dval(lay)),lay=1,numlay)

      write (afo,4020) (0.01*real(dzNew(node)),node=1,numnodNew)

      if (swafo.ge.2 .and. swop.eq.2) then
         write (afo,4020)                                               &
     &         (0.01*real(VlMpStDm1New(node)),node=1,numnodNew)
         write (afo,4020)                                               &
     &         (0.01*real(VlMpStDm2New(node)),node=1,numnodNew)
         write (afo,4020) (0.01*real(DiPoCpNew(node)),node=1,numnodNew)
      endif

      write (afo,4020) (real(thetaNew(node)),node=1,numnodNew)
      write (afo,4020) -0.01*real(gwl),0.01*real(pond)

      if (swafo.ge.2) then
         write (afo,4020) 0.01*real(Ssnow)
         write (afo,4020) (real(tsoilNew(node)),node = 1,numnodNew)
      endif

      if (swafo.ge.2 .and. swop.eq.2) then
         write (afo,4020) 0.01*real(WaLevDm1), 0.01*real(VlMpDm1),      &
     &       0.01*real(WaSrDm1), 0.01*real(VlMpDm2), 0.01*real(WaSrDm2)
      endif

4000  format (2(1x,i8),3(1x,f8.0),1x,i8)
4001  format (1x,i8)
4010  format (5(1x,i8))
4020  format (8(1x,f10.6))

      FlOpenFileDev = .false.

      return

      case (2)

! === output during simulation ===============================================

! -     convert vertical discretization: dynamic part
        do node = 1,numnod
          tsoili(node) = tsoil(node)
        enddo
        call ConvertDiscrVert(2,swop,botcomnew,hnew,thetanew,inqnew,inqrotnew,inqdranew,ithetabegnew,tsoili,tsoilnew,    &
                              dipocpnew,iavfrmpwlwtdm1new,iavfrmpwlwtdm2new,inqexcmtxdm1cpnew,inqexcmtxdm2cpnew,  &
                              inqoutdrrapcpnew,vlmpstdm1new,vlmpstdm2new)

        if (swafo.eq.1) then
          write (afo,30) real(tcum),                                    &
     &    0.01*real((igrai+igird)/outper),                              &
     &    0.01*real(iintc/outper),                                      &
     &    0.01*real(ievap/outper),                                      &
     &    0.0,                                                          &
     &    0.01*real(ipeva/outper),0.01*real(iptra/outper),              &
     &    0.01*real(iruno/outper),                                      &
     &    -0.01*real(gwl),0.01*real(pond)
        elseif (swafo.ge.2) then
!          write (afo,*) real(t), real(outper),
          write (afo,33) (brund*1.0d0-1.0d0+tcum), real(outper),        &
     &        0.01*real((igrai+isnrai)/outper),                         &
     &        0.01*real(igsnow/outper),                                 &
     &        0.01*real(igird/outper),                                  &
     &        0.01*real( (iintc-(igird-inird))/outper ),                &
     &        0.01*real( (igird-inird)/outper ),                        &
     &        0.01*real(isubl/outper),                                  &
     &        0.01*real(ievap/outper),                                  &
     &        0.0,                                                      &
     &        0.01*real(ipeva/outper),                                  &
     &        0.01*real(iptra/outper),                                  &
     &        0.01*real(irunon/outper),                                 &
     &        0.01*real(iruno/outper),                                  &
     &        -0.01*real(gwl), 0.01*real(pond),                         &
     &        0.01*real(ssnow), 0.01*real(wbalance)
        endif
        write (afo,40) (real(hNew(node))           ,node=1,numnodNew)
        write (afo,50) (real(thetaNew(node))       ,node=1,numnodNew)
        write (afo,50) (0.01*real(inqrotNew(node)/outper),node=1,       &
     &   numnodNew)
        write (afo,50) (-0.01*real(inqNew(node)/outper),                &
     &   node=1,numnodNew+1)
        if (nrlevs.gt.0) then
          do level=1,nrlevs
            write (afo,50) (0.01*real(inqdraNew(level,node)/outper),    &
     &        node=1,numnodNew)
          enddo
        end if
        if (swafo.ge.2) then
!ckro-20141112   exponential relation between soil cover and lai
          gc = min( (1.0d0-exp(-1.0d0*kdif*kdir*lai)), 1.0d0)
          write (afo,'(4(1x,f10.4))') real(gc),real(lai),               &
     &     abs(real(0.01d0*rd)),real(cf)
! -   CFO file for PEARL
         if (swafo.eq.2) then
            write (afo,'(1x,f7.2)') real(tav)
         else if (swafo.eq.3) then
! -   if swafo = 3, CFO file with extra output concerning meteo:
!     RAD(kJ.m-2), Tmin(grC), Tmax(grC), HUM(kPA), WIND(m.s-1), ETref (m.d-1), WET (d)
!     and missing values voor WET: -1.0 if SWRAIN =/= 2
            write (afo,'(3x,f8.1,f9.2,2(2x,f9.2),f9.2,3x,2(3x,f8.5))')  &
     &     out_rad, out_tmn, out_tmx, out_hum, out_win, out_etr, out_wet
         endif
          write (afo,4020) (real(tsoilNew(node)),node = 1,numnodNew)
        endif

        if (swafo.ge.2 .and. swop.eq.2) then

          if (IcTopMp.ne.1) then              ! adaptation for GEM; covering lay on top of macropores:
             TopQsoilDm1 = IQInTopVrtDm1      ! request to put inflow into macropores from covering top layer
             TopQprecDm1 = 0.0d0              ! in variables IQInTopLatDm1 and IQInTopLatDm2 (the fluxes that contain pesticides in PEARL);
             TopQsoilDm2 = IQInTopVrtDm2      ! for calculations of macropore fluxes this inflow is in variables
             TopQprecDm2 = 0.0d0              ! IQInTopVrtDm1 and IQInTopVrtDm2
!            write(117,'(5f8.3)') tcum,IQInTopVrtDm1, IQInTopLatDm1, IQInTopVrtDm2, IQInTopLatDm2
          else
             TopQprecDm1 = IQInTopVrtDm1
             TopQsoilDm1 = IQInTopLatDm1
             TopQprecDm2 = IQInTopVrtDm2
             TopQsoilDm2 = IQInTopLatDm2
          endif
!
          write (AFO,4020) 0.01*real(WaLevDm1), 0.01*real(VlMpDm1),     &
     &             0.01*real(WaSrDm1), 0.01*real(TopQprecDm1/outper),   &
     &             0.01*real(TopQsoilDm1/outper)
          write (AFO,4020) (0.01*real(InQExcMtxDm1CpNew(node)/outper),  &
     &                      node=1,numnodNew)
          write (AFO,4020) (0.01*real(InQOutDrRapCpNew(node)/outper),   &
     &                      node=1,numnodNew)
          write (AFO,4020) (real(IAvFrMpWlWtDm1New(node)/outper),       &
     &                      node=1,numnodNew)
          write (AFO,4020) 0.01*real(VlMpDm2), 0.01*real(WaSrDm2),      &
     &                     0.01*real(TopQprecDm2/outper),               &
     &                     0.01*real(TopQsoilDm2/outper)
          write (AFO,4020) (0.01*real(InQExcMtxDm2CpNew(node)/outper),  &
     &                      node=1,numnodNew)
          write (AFO,4020) (real(IAvFrMpWlWtDm2New(node)/outper),       &
     &                      node=1,numnodNew)
        endif

        if (swafo.ge.2 .and. CritDevMasBal.gt.1.d-30) then
! --- Checking of mass balances of sub systems per period OutPer
           call checkmassbal (flopenfiledev,inqdranew,inqexcmtxdm1cpnew,inqexcmtxdm2cpnew,inqnew,inqoutdrrapcpnew,inqrotnew,ithetabegnew,thetanew)
        end if

30    format (1x,f6.0,7(1x,f9.6),3(1x,f7.4))
33    format (1(1x,f12.6),12(1x,f9.6),4(1x,f9.4))
40    format (8(1x,e10.3e2))
50    format (10(1x,f9.6))
!30    format (1x,f6.0,7(1x,f9.6),3(1x,f7.4))   ! 7
!33    format (1(1x,f12.6),12(1x,f9.6),4(1x,f9.4))
!40    format (8(1x,e10.3e2))  ! 9
!50    format (10(1x,f9.6))    ! 8

      case default
         call fatalerr ('outafo', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outaun (task)
! ----------------------------------------------------------------------
!     Date               : 29-jan-2003
!     Purpose            : ANIMO/PEARL output: unformatted hydrological data
! ---------------------------------------------------------------------
      use variables, only: aun,outfil,pathwork,numnod,outper,period,pond,gwl,nrlevs,ievap,ipeva,iptra,    &
                           iruno,numlay,botcom,thetas,kdif,kdir,swaun,igrai,igird,inird,iintc,gc,lai,tav,tsoil,rd,cf,wbalance, &
                           project,tstart,tend,swdiscrvert,numnodnew,dznew,SwAfo,FlMacropore,Ssnow,igSnow,isnrai,iSubl,irunon, &
                           WalevDm1,VlMpDm1,WaSrDm1,VlMpDm2,WaSrDm2,IQInTopVrtDm1,IQInTopLatDm1,   &
                           IQInTopVrtDm2,IQInTopLatDm2,CritDevMasBal,tcum,nod1lay, FlOpenFileDev
      implicit none
      include 'arrays.fi'

! -   global
      integer   task

! -   local
      integer   swop
      integer   getun,ftype,node,bruny,eruny,brund,erund,lay,level
      integer   botcomNew(maho)
      integer   datea(6)
      real(4)   fsec
      real(8)   timjan1,dval(maho),watcon
      real(8)   hNew(macp),inqNew(macp+1)
      real(8)   inqdraNew(Madr,macp),thetaNew(macp),inqrotNew(macp)
      real(8)   TsoilNew(0:macp),tsoili(0:macp)
      real(8)   DiPoCpNew(macp), IAvFrMpWlWtDm1New(macp)
      real(8)   IAvFrMpWlWtDm2New(macp), InQExcMtxDm1CpNew(macp)
      real(8)   InQExcMtxDm2CpNew(macp), InQOutDrRapCpNew(macp)
      real(8)   IThetaBegNew(MaCp),VlMpStDm1New(macp),VlMpStDm2New(macp)
      character(len=80)  filtext
      character(len=300) filnam
      character(len=4)   aunext

!     save local variables
      save      brund, swop

! ----------------------------------------------------------------------

      select case (task)
      case (1)

! === output initial condition ==========================================

! -   verify reduced grid for output
      if(SwDiscrvert.eq.1) then
         call checkDiscrVert()
      endif

! -   set switch Swop
      Swop = 1
      if (Swaun.eq.2 .and. FlMacropore) Swop = 2

! -     convert vertical discretization: initial values
      do node = 1,numnod
        tsoili(node) = tsoil(node)
      enddo
      call convertdiscrvert(1,swop,botcomnew,hnew,thetanew,inqnew,inqrotnew,inqdranew,ithetabegnew,tsoili,tsoilnew,      &
                            dipocpnew,iavfrmpwlwtdm1new,iavfrmpwlwtdm2new,inqexcmtxdm1cpnew,inqexcmtxdm2cpnew,    &
                            inqoutdrrapcpnew,vlmpstdm1new,vlmpstdm2new)

! -   open file
      if (swaun.eq.1) then
        aunext = '.aun'
      elseif (swaun.eq.2) then
        aunext = '.bun'
      endif
      filnam = trim(pathwork)//trim(outfil)//aunext
      aun = getun (20,90)
      open(unit=aun,file=filnam,status='unknown',form='unformatted')

! --- write initial part
      if (swaun.eq.2) then
          ftype = 2
          filtext = 'unformatted hydrological data'
          call writehead(aun,ftype,filnam,filtext,project)
      endif

      call dtdpar (tstart,datea,fsec)
      bruny = datea(1)
      datea(2) = 1
      datea(3) = 1
      fsec = 0.0
      call dtardp (datea,fsec,timjan1)
      brund = nint (tstart - timjan1 + 1.0d0)
      call dtdpar (tend,datea,fsec)
      eruny = datea(1)
      datea(2) = 1
      datea(3) = 1
      fsec = 0.0
      call dtardp (datea,fsec,timjan1)
      erund = nint (tend - timjan1 + 1.0d0)
      if (swaun.eq.1) then
        write (aun) bruny,eruny,                                        &
     &                   real(brund-1),real(erund),real(period)
      else if (swaun.eq.2) then
        write (aun) swop
        write (aun) bruny,eruny,real(brund-1),real(erund)
      endif

      write (aun) numnodNew,numlay,nrlevs
      write (aun) (botcomNew(lay),lay=1,numlay)
!  -  ThetaS should be known for new soil layers; this works but is not very nice
      write (aun) (real(thetas(botcom(lay))), lay=1,numlay)

      do lay = 1,numlay
!        find first Node of the Layer
         Node = nod1lay(lay)
         dval(lay) = watcon(node,-100.d0)
      enddo
      write (aun) (real(dval(lay)),lay=1,numlay)

      do lay = 1,numlay
!        find first Node of the Layer
         Node = nod1lay(lay)
         dval(lay) = watcon(node,-15849.0d0)
      enddo
      write (aun) (real(dval(lay)),lay=1,numlay)

      write (aun) (0.01*real(dzNew(node)),node=1,numnodNew)

      if (swaun.eq.2 .and. swop.eq.2) then
        write (aun) (0.01*real(VlMpStDm1New(node)),node=1,numnodNew)
        write (aun) (0.01*real(VlMpStDm2New(node)),node=1,numnodNew)
        write (aun) (0.01*real(DiPoCpNew(node)),node=1,numnodNew)
      endif

      write (aun) (real(thetaNew(node)),node=1,numnodNew)
      write (aun) -0.01*real(gwl),0.01*real(pond)

      if (swaun.eq.2) then
        write (aun) 0.01*real(Ssnow)
        write (aun) (real(tsoilNew(node)),node = 1,numnodNew)
      endif

      if (swaun.eq.2 .and. swop.eq.2) then
        write (aun) 0.01*real(WaLevDm1), 0.01*real(VlMpDm1),            &
     &        0.01*real(WaSrDm1), 0.01*real(VlMpDm2), 0.01*real(WaSrDm2)
      endif

      FlOpenFileDev = .false.

      return

      case (2)

! === output during simulation ===============================================

! -   convert vertical discretization: dynamic part
      do node = 1,numnod
        tsoili(node) = tsoil(node)
      enddo
      call convertdiscrvert(2,swop,botcomnew,hnew,thetanew,inqnew,inqrotnew,inqdranew,ithetabegnew,tsoili,tsoilnew,      &
                            dipocpnew,iavfrmpwlwtdm1new,iavfrmpwlwtdm2new,inqexcmtxdm1cpnew,inqexcmtxdm2cpnew,    &
                            inqoutdrrapcpnew,vlmpstdm1new,vlmpstdm2new)

      if (swaun.eq.1) then
        write (aun) real(tcum),                                         &
     &    0.01*real((igrai+igird)/outper),                              &
     &    0.01*real(iintc/outper),                                      &
     &    0.01*real(ievap/outper),                                      &
     &    0.0,                                                          &
     &    0.01*real(ipeva/outper),                                      &
     &    0.01*real(iptra/outper),                                      &
     &    0.01*real(iruno/outper),                                      &
     &    -0.01*real(gwl),0.01*real(pond)
      elseif (swaun.eq.2) then
        write (aun) real(brund*1.0d0-1.0d0+tcum), real(outper),         &
     &        0.01*real((igrai+isnrai)/outper),                         &
     &        0.01*real(igsnow/outper),                                 &
     &        0.01*real(igird/outper),                                  &
     &        0.01*real( (iintc-(igird-inird))/outper ),                &
     &        0.01*real( (igird-inird)/outper ),                        &
     &        0.01*real(isubl/outper),                                  &
     &        0.01*real(ievap/outper),                                  &
     &        0.0,                                                      &
     &        0.01*real(ipeva/outper),                                  &
     &        0.01*real(iptra/outper),                                  &
     &        0.01*real(irunon/outper),                                 &
     &        0.01*real(iruno/outper),                                  &
     &        -0.01*real(gwl), 0.01*real(pond),                         &
     &        0.01*real(ssnow), 0.01*real(wbalance)
      endif
      write (aun) (real(hNew(node))                 ,node=1,numnodNew)
      write (aun) (real(thetaNew(node))             ,node=1,numnodNew)
      write (aun) (0.01*real(inqrotNew(node)/outper),node=1,numnodNew)
      write (aun) (-0.01*real(inqNew(node)/outper),node=1,numnodNew+1)
      if (nrlevs.gt.0) then
        do level=1,nrlevs
          write (aun) (0.01*real(inqdraNew(level,node)/outper),         &
     &        node=1,numnodNew)
        enddo
      end if
      if (swaun.eq.2) then
!ckro-20141112   exponential relation between soil cover and lai
        gc = min( (1.0d0-exp(-1.0d0*kdif*kdir*lai)), 1.0d0)
        write (aun) real(gc),real(lai),abs(real(0.01d0*rd)),real(cf)
        write (aun) real(tav)
        write (aun) (real(tsoilNew(node)),node = 1,numnodNew)
      endif
      if (swaun.eq.2 .and. swop.eq.2) then
        write (AUN) 0.01*real(WaLevDm1), 0.01*real(VlMpDm1),            &
     &             0.01*real(WaSrDm1), 0.01*real(IQInTopVrtDm1/outper), &
     &             0.01*real(IQInTopLatDm1/outper)
        write (AUN) (0.01*real(InQExcMtxDm1CpNew(node)/outper),         &
     &                      node=1,numnodNew)
        write (AUN) (0.01*real(InQOutDrRapCpNew(node)/outper),          &
     &                      node=1,numnodNew)
        write (AUN) (real(IAvFrMpWlWtDm1New(node)/outper),              &
     &                      node=1,numnodNew)
        write (AUN) 0.01*real(VlMpDm2), 0.01*real(WaSrDm2),             &
     &                     0.01*real(IQInTopVrtDm2/outper),             &
     &                     0.01*real(IQInTopLatDm2/outper)
        write (AUN) (0.01*real(InQExcMtxDm2CpNew(node)/outper),         &
     &                      node=1,numnodNew)
        write (AUN) (real(IAvFrMpWlWtDm2New(node)/outper),              &
     &                      node=1,numnodNew)
      end if

      if ((swaun.ge.1.or.swafo.ge.1) .and. CritDevMasBal.gt.1.d-30) then

! --- Checking of mass balances of sub systems per period OutPer
        call checkmassbal (flopenfiledev,inqdranew,inqexcmtxdm1cpnew,inqexcmtxdm2cpnew,inqnew,inqoutdrrapcpnew,inqrotnew,ithetabegnew,thetanew)

      endif

      case default
         call fatalerr ('outaun', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine WriteSwapOk(Project)
!-----------------------------------------------------------------------
!      Date    : December 2004
!      Purpose : Create file Swap.ok at end of simulation
!-----------------------------------------------------------------------
      implicit none

      integer   cexf,getun
      character(len=80) project

!     create file Swap.ok to let environment programs verify termination
      cexf = getun (20,90)
      call fopens(cexf,'Swap.ok','new','del')
      call writehead (cexf,1,'Swap.ok',                                 &
     &  'this header only: simulation succesfully terminated',project)
      close (cexf)

      return
      end


! ----------------------------------------------------------------------
      subroutine SurfaceWaterOutput(task)
! ----------------------------------------------------------------------
!     Date               : Aug 2004
!     Purpose            : open and write surface water output files
! ----------------------------------------------------------------------

      use variables
      implicit none

      integer task

      select case (task)
      case (1)

! === open output files and write headers ===============================

! --  drf file
      if (swdrf.eq.1) call outdrf (task)

      if (swswb.eq.1) call outswb(task)

      return

      case (2)

! === write actual data ===============================

! --  drf file
      if (swdrf.eq.1) call outdrf (task)

      if (swswb.eq.1) call outswb(task)

      return

      case (3)

! === close output files ===========================

      close (drf)
      close (swb)

      case default
         call fatalerr ('SurfaceWaterOutput', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outdrf (task)
! ----------------------------------------------------------------------
!     Date               : 10/6/99
!     Purpose            : write drainage fluxes, surface runoff, rapid
!                          drainage to  OUTNAM.DRF file
! ---------------------------------------------------------------------
      use variables, only: outfil,drf,pathwork,daynr,date,nrpri,nrlevs,cqdrain,cqdrd,crunoff,cQMpOutDrRap,flheader
      implicit none
      include 'arrays.fi'

! --- global
      integer   task

! --- local
      integer   level,getun
      real(8)   c1qdrain(Madr),c1qdrd,c1runoff,c1qdrar
      character(len=80) filnam
      character(len=1)  comma

      save
      comma = ','
! ----------------------------------------------------------------------

      select case (task)
      case (1)
! === open output file and write header ================================

! --- open output file
      filnam = trim(pathwork)//trim(outfil)//'.drf'
      drf = getun (20,90)
      call fopens(drf,filnam,'new','del')

! --- write header
      if (nrpri .eq. 1) then
        write (DRF,100) nrlevs
      else
        write (DRF,110) nrlevs
      endif
      write (DRF,120)

! --- Output format
 100    format (                                                        &
     &'* CUMULATIVE DRAINAGE FLUXES, SURFACE RUNOFF AND RAPID',         &
     &' DRAINAGE:',/,'*',/,                                             &
     &'* Total number of levels for drainage fluxes :',I3,/,'*',/       &
     &'* First level (primary system) NOT included in sw-reservoir',/,  &
     &'*',/)

 110    format (                                                        &
     &'* CUMULATIVE DRAINAGE FLUXES, SURFACE RUNOFF AND RAPID',         &
     &' DRAINAGE:',/,'*',/,                                             &
     &'* Total number of levels for drainage fluxes :',I3,/,'*',/       &
     &'* First level (primary system) IS included in sw-reservoir',/,   &
     &'*')

 120    format (                                                        &
     &'* Meaning of symbols:',/,                                        &
     &'* - CQDinc(1-5): Drainage fluxes all level (1-5)',/,             &
     &'* - CQDRDinc   : Total of drainage fluxes into surface',         &
     &            ' water reservoir (SWSRF>=1)',/,                      &
     &'* - CRUNOFFinc : Surface runoff (increment can be <0.0)',/,      &
     &'* - CQDRARinc  : Rapid drainage (increment always > = 0.0)',/,   &
     &'* - CQDcum(1-5): Drainage fluxes all level (1-5)',/,             &
     &'* - CQDRDcum   : Total of drainage fluxes into surface',         &
     &            ' water reservoir (SWSRF>=1)',/,                      &
     &'* - CRUNOFFcum : Surface runoff (increment can be <0.0)',/,      &
     &'* - CQDRARcum  : Rapid drainage (increment always > = 0.0)',/,   &
     &'* - cum        : cumulative value',/,                            &
     &'* - inc        : incremental value',/,'*',/,                     &
     &'     Date, Daynr,CQDinc1,CQDinc2,CQDinc3,CQDinc4,CQDinc5,',      &
     &'CQDRDinc,CRUNOFFinc,CQDRARinc,',                                 &
     &'CQDcum1,CQDcum2,CQDcum3,CQDcum4,CQDcum5,',                       &
     &'CQDRDcum,CRUNOFFcum,CQDRARcum')

      do level=1,nrlevs
        c1qdrain(level)=0.0d0
      enddo
      c1qdrd = 0.0d0
      c1runoff = 0.0d0
      c1qdrar = 0.0d0

      return

      case (2)

! === output of data ===========================================================

! --- write header and reset cumulative fluxes in case of new balance period
      if (flheader) then
        if (nrpri .eq. 1) then
          write (DRF,100) nrlevs
        else
          write (DRF,110) nrlevs
        endif
        write (DRF,120)

        do level=1,nrlevs
          c1qdrain(level)=0.0d0
        enddo
        c1qdrd = 0.0d0
        c1runoff = 0.0d0
        c1qdrar = 0.0d0
      endif

! --- write output record
      write (DRF,20) date,comma,daynr,                                  &
     &  comma,(cqdrain(1)-c1qdrain(1)),                                 &
     &  comma,(cqdrain(2)-c1qdrain(2)),                                 &
     &  comma,(cqdrain(3)-c1qdrain(3)),                                 &
     &  comma,(cqdrain(4)-c1qdrain(4)),                                 &
     &  comma,(cqdrain(5)-c1qdrain(5)),                                 &
     &  comma,(cqdrd-c1qdrd),comma,(crunoff-c1runoff),                  &
     &  comma,(cQMpOutDrRap-c1qdrar),                                   &
     &  comma,cqdrain(1),comma,cqdrain(2),                              &
     &  comma,cqdrain(3),comma,cqdrain(4),comma,cqdrain(5),             &
     &  comma,cqdrd,comma,crunoff,comma,cQMpOutDrRap

 20   format (A11,a1,I4,8(a1,f8.2),8(a1,f8.1))

! --- store cumulative values
      do level=1,nrlevs
          c1qdrain(level)=cqdrain(level)
      enddo
      c1qdrd = cqdrd
      c1runoff = crunoff
      c1qdrar = cQMpOutDrRap

      case default
         call fatalerr ('outdrf', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine outswb(task)
! ----------------------------------------------------------------------
!     Date               : 21/08/99
!     Purpose            : write surface water balance data to
!                          OUTNAM.SWB file, surface water management
!                          data to OUTNAM.MAN. The files overlap
! ---------------------------------------------------------------------
      use variables, only: outfil,vtair,pathwork,daynr,daycum,hbweir,overfl,gwl,pond,wlstar,wls,swstini,swst,cqdrd,crunoff,           &
                           cQMpOutDrRap,cwsupp,swb,numadj,hwlman,cwout,swsec,swman,nmper,impend,imper,project,logf,swscre,date,t1900,t,outper,iyear
      implicit none

! --- Global
      integer   task

! --- Local
      integer   getun,nrOfDays,man
      real(8)   gwlev,cqdrf1,c1wsupp,c1wout,delbal,small,zero
      character(len=1)   spc,comma
      character(len=200) messag
      character(len=132) filnam
      character(len=80)  filtext
      character(len=19) datetime
      logical   dtleap

      save

      data      small     /0.0001d0/
      data      zero      /0.0d0/
      comma = ','
! ----------------------------------------------------------------------

      select case (task)
      case (1)
! === open output files and write headers ==============================

! --- open output file once
      filnam = trim(pathwork)//trim(outfil)//'.swb'
      swb = getun (20,90)
      call fopens(swb,filnam,'new','del')
      filtext = 'Surface water balance increments (cm/period)'
      call writehead (swb,1,filnam,filtext,project)

! --- write header
      write (swb,10)
      write (swb,11)
 10   format ('*',/,                                                    &
     &'* Surface water system:',/,                                      &
     &'* - GWL           : groundwater level + ponding (<0. = bel',     &
     &'ow s.s) (cm)',/,                                                 &
     &'* - WLST          : target level,  aut. weir (<0. = b.s.s.',     &
     &') (cm)',/,                                                       &
     &'* - WLST          : weir crest,   fixed weir (<0. = b.s.s.',     &
     &')  (cm)',/,                                                      &
     &'* - WLS           : surface water level (<0. = b.s.s.)  (c',     &
     &'m)',/,                                                           &
     &'* - swst           : storage in sw-reservoir, per unit area',    &
     &' of* subcatchment (cm)',/,                                       &
     &'* - DRORR         : drainage fluxes (>0. = into sw) + runo',     &
     &'ff +* rapid drainage (cm)',/,                                    &
     &'* - QSUPP         : external supply to sw-reservoir (cm)',/,     &
     &'* - QOUT          : outflow from sw-reservoir (cm)',/,           &
     &'* - DRORRcum      : drainage fluxes (>0. = into sw) + runo',     &
     &'ff +* rapid drainage (cumulative, cm)')
 11   format (                                                          &
     &'* - QSUPPcum      : external supply to sw-reservoir (cumul',     &
     &'ative, cm)',/,                                                   &
     &'* - QOUTcum       : outflow from sw-reservoir (cumulative,',     &
     &' cm)',/,'*')

! ---   determine which management period the model is in:
      if (swsec .eq. 2) then
        imper = 0
 100    imper = imper + 1
!       Error handling
        if (imper .gt. nmper) then
          messag = 'error sw-management(IMPER), more than defined'
          call fatalerr ('Outswb',messag)
        endif
        if (t1900-1.d0-0.1d-10 .gt. impend(imper)) goto 100
      endif

! --- add ponding to groundwater level, if gwl at soil surface:
      if (abs(gwl) .lt. 1.0d-7) then
          gwlev = pond
      else
          gwlev = gwl
      endif

! --- write record for initial state
      write (SWB,20) comma,gwlev,comma,wlstar,comma,wls,comma,swst,     &
     & comma,(0.),comma,(0.),comma,(0.),comma,(0.),comma,(0.),comma,(0.)
 20   format ('* initial  ,   - ,    - ',                               &
     & 3(a1,f7.1),a1,f6.1,6(a1,f7.2),/,'*',/,                           &
     &'      Date,Daynr,DayCum,    GWL,   WLST,    WLS,   swst,  ',     &
     &'DRORR,  QSUPP,  QOUT,DRORRcum,QSUPPcum,QOUTcum')

! --- compose filename for management data
      if (SWSEC .eq. 2) then

! --    open output file once
        filnam = trim(pathwork)//trim(outfil)//'.man'
        man = getun (20,90)
        call fopens(man,filnam,'new','del')
        filtext = 'Surface water management (cm/period)'
        call writehead (man,1,filnam,filtext,project)

! --    write header
        write (MAN,30)
 30     format ('*',/,                                                  &
     &'* Surface water management:',/,                                  &
     &'* - SWMAN         : type of weir (f = fixed; a = automatic)',/,  &
     &'* - GWL           : groundwater level + ponding (<0. = belo',    &
     &'w s.s) (cm)',/,                                                  &
     &'* - HWLMAN        : pressure head used for target level (cm)',/, &
     &'* - VTAIR         : total air volume in soil profile (cm)',/,    &
     &'* - WLST          : target level of autom. weir (<0. = b.s.',    &
     &'s.)  (cm)',/,                                                    &
     &'*              or : crest  level of fixed  weir (<0. = b.s.',    &
     &'s.)  (cm)',/,                                                    &
     &'* - WLS           : surface water level (<0. = b.s.s.) (cm)',/,  &
     &'* - QOUT          : surface water outflow (<0. = supply) (cm)',/,&
     &'* - NUMADJ        : number of adjustments of target level (',    &
     &'cm)',/,                                                          &
     &'* - OVERFL        : flag for overflow of aut. weir (o)',/,       &
     &'* - CREST         : Crest level (<0. = b.s.s.) (cm)',/,          &
     &'*',/,                                                            &
     &'     Date,DayNr,DayCum,SWMan,    GWL,HWLMAN, VTAIR,   WLST,',    &
     &'    WLS,   QOUT, NUMADJ,OVERFL,  CREST')
      endif

      return

      case (2)

! === actual output data =================================================

! --- reset cumulative fluxes each year
      if (abs(t-outper).lt.small) then
        cqdrf1 = zero
        c1wsupp = zero
        c1wout = zero
      endif

! --- water balance error
      delbal = (swst + cwout) - (swstini + cqdrd+crunoff + cwsupp +     &
     &                          cQMpOutDrRap)
      if (delbal .gt. 0.05d0) then
        call dtdpst                                                     &
     &        ('year-month-day,hour:minute:seconds',t1900,datetime)
        write(messag,1001) delbal,datetime
 1001   format('  Error in cumulative water balance of surfwater:',     &
     &         f5.2,'   at : ',a19)
        call warn ('Outswb',messag,logf,swscre)
      endif

! --- add ponding to groundwater level, if gwl at soil surface:
      if (abs(gwl) .lt. 1.0d-7) then
          gwlev = pond
      else
          gwlev = gwl
      endif

! --- write output record OUTNAM.SWB
      write (SWB,40) date,comma,daynr,comma,daycum,comma,gwlev,comma,   &
     &  wlstar,comma,wls,comma,swst,comma,                              &
     &  (cqdrd+crunoff+cQMpOutDrRap-cqdrf1),comma,(cwsupp-c1wsupp),     &
     &  comma,(cwout-c1wout),comma,(cqdrd+crunoff+cQMpOutDrRap),        &
     &  comma,cwsupp,comma,cwout
 40   format (a11,a1,I4,a1,I6,3(a1,f7.1),a1,f6.1,6(a1,f7.2))

! --- write output record OUTNAM.MAN
      if (swsec .eq. 2) then
        if (overfl) then
          spc = 'o'
        else
          spc = '-'
        endif
        if (swman(imper) .eq. 1) then
          write(man,50) date,comma,daynr,comma,daycum,comma,gwlev,comma,&
     &      wlstar,comma,wls,comma,((cwout-c1wout)-(cwsupp-c1wsupp)),   &
     &      comma,numadj,comma,spc,comma,hbweir(imper)
        else
          write(man,60) date,comma,daynr,comma,daycum,comma,gwlev,comma,&
     &      hwlman,comma,vtair,comma,wlstar,comma,wls,comma,            &
     &      ((cwout-c1wout)-(cwsupp-c1wsupp)),comma,numadj,comma,spc,   &
     &      comma,hbweir(imper)
        endif
 50     format (a11,a1,i4,a1,I6,',    f',a1,f7.1,2(',     -'),          &
     &          3(a1,f7.1),a1,i7,a1,5x,a1,a1,f7.1)

 60     format (a11,a1,i4,a1,I6,',    a',a1,f7.1,2(a1,f7.1),            &
     &          3(a1,f7.1),a1,i7,a1,5x,a1,a1,f7.1)
      endif

! --- store cumulative values
      cqdrf1 = (cqdrd+crunoff+cQMpOutDrRap)
      c1wsupp = cwsupp
      c1wout = cwout

! --- saves storage as initial values for next year
      nrOfDays = 365
      if (dtleap(iyear)) nrOfDays = nrOfDays + 1
      if (daynr.eq.nrOfDays) swstini = swst

      case default
         call fatalerr ('outswb', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine warn (modul,messag,logf,swscre)
      implicit none
! ----------------------------------------------------------------------
! --- global
      integer       logf,swscre
      character(len=*) modul,messag

! --- local
      character(len=400)  messages
! ----------------------------------------------------------------------

      messages = 'Warning from module '//modul//' : '//trim(messag)

      write (logf,'(a)') trim(messages)
      if (swscre .gt. 0) write (*,'(2x,a)') trim(messages)

      return
      end

      subroutine writehead(outf,ftype,filnam,filtext,project)
!-----------------------------------------------------------------------
!      Date    : January 2006
!      Purpose : Writes a header to output files of SWAP

!-----------------------------------------------------------------------
      implicit none

!     global
      integer       outf,ftype
      character(len=*) project,filnam,filtext

!     local
      integer       date_time(6)
      real(8)       dpactualtime
      real(4)       dum
      character(len=80)  model_id,dtstring,String
      character(len=132) Version
!-----------------------------------------------------------------------

! --- Version nr of model, to appear in all output files
!     together with revision nr
      include 'description.fi'

      model_id = 'Swap '//trim(Version)

! --- get actual time
      call dtnow (date_time)
      dum = 0.0
      call dtardp (date_time, dum, dpactualtime)
      call dtdpst ('year-month-day hour:min:sec', dpactualtime,       &
     &             dtstring)

      if (outf.eq.5) then
! ---   write to screen
        write (*,16)  trim(project)
        write (*,17)  trim(filtext)
        write (*,19)  trim(model_id)
        write (*,20)  trim(dtstring)
      else if (ftype.eq.1) then
!       write formatted output file
        write (outf,16)  trim(project)
        write (outf,17)  trim(filtext)
        write (outf,18)  trim(filnam)
        write (outf,19)  trim(model_id)
        write (outf,20)  trim(dtstring)
      else
!       write unformatted output file (header has fixed length of 80 characters)
        write (String,'(a80)')  '* Project:       '//trim(project)
        write (outf)  adjustl(String)
        write (String,'(a80)')  '* File content:  '//trim(filtext)
        write (outf)  adjustl(String)
        write (String,'(a80)')  '* File name:     '//trim(filnam)
        write (outf)  adjustl(String)
        write (String,'(a80)')  '* Model version: '//trim(model_id)
        write (outf)  adjustl(String)
        write (String,'(a80)')  '* Generated at:  '//trim(dtstring)
        write (outf)  adjustl(String)
      endif

 16   format('* Project:       ',a)
 17   format('* File content:  ',a)
 18   format('* File name:     ',a)
 19   format('* Model version: ',a)
 20   format('* Generated at:  ',a)

      return
      end

! ----------------------------------------------------------------------
      subroutine CloseTempFil
! ----------------------------------------------------------------------
!     date               : Aug 2004
!     purpose            : delete temporary files of TTUTIL
! ----------------------------------------------------------------------
      implicit  none

      integer cexf,getun
      logical fileopen

! --- delete temporary files
      cexf = getun (10,90)
      call rddtmp(cexf)
      inquire(unit=20,opened=fileopen)
      if(fileopen)close(20, STATUS = 'DELETE')

      return
      end


! ----------------------------------------------------------------------
      subroutine checkDiscrVert()
!     date               : 20081105
!     purpose            : verify reduced vertical discretizationface,
! global   formal parameters  : (i = input, o = output)
!     numnod       ! Number of nodes or compartments....................... i
!     numnodnew    ! Number of desired nodes for soil water quality models..i
!     dz(macp)     ! Compartment thickness (L) ............................ i
!     dznew(macp)  ! Desired dz for soil water quality models (L) ......... i
! local
! ----------------------------------------------------------------------
      use variables, only: numnod,dz,numnodnew,dznew
      implicit none
      include 'arrays.fi'

! global

! local
      integer   in,io,iotmp
      real(8)   cumdzN(macp),cumdzOld,cumdzNew
      character(len=200) messag
! local parameters
      real(8)   small
      data      small     /0.0001d0/
! ----------------------------------------------------------------------


! --  boundary of new compartments must equal boundaries of old compartmts
!     or: new size is sum of old sizes
      iotmp = 0
      do in = 1,numnodNew
        cumdzN(in) = 0.0d0
        do io = 1,numnod
          if(io.gt.iotmp .and. cumdzN(in).lt.dzNew(in)) then
             cumdzN(in) = cumdzN(in) + dz(io)
             iotmp = io
          endif
        enddo
      enddo
      do in = 1,numnodNew
        if (abs(dzNew(in)-cumdzN(in)).gt.small) then
          write(messag,'(a,i5)')                                        &
     &    'New discretization has error in size of new compartment ',in
          call fatalerr ('checkDiscrVert',messag)
        endif
      enddo

! --  cumulative thickness
      cumdzOld = 0.0d0
      do io = 1,numnod
         cumdzOld = cumdzOld + dz(io)
      enddo
      cumdzNew = 0.0d0
      do in = 1,numnodNew
         cumdzNew = cumdzNew + dznew(in)
      enddo
      if (abs(cumdzNew-cumdzOld).gt.small) then
        write(messag,'(a,f10.5)')                                       &
     &    'New discretization is wrong, total length differs >',small
        call fatalerr ('checkDiscrVert',messag)
      endif


      return
      end

! ----------------------------------------------------------------------
      subroutine OutputModflow(task)
! ----------------------------------------------------------------------
!     Date               : July 2009
!     Purpose            : open and write data to explore storage/recharge
! ----------------------------------------------------------------------

      use variables
      implicit none

! --- global variables ------------------
      integer task
! --- local variables ------------------
      integer   getun,sto,nod1m, nod
      real(8)   vsat,vt0,gwlt0,vt1,qre,qv1m,gwlt1,stocoav,stocot1
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
      integer   i,swBotbtmp
      real(8)   vair(2),dgwl(2),gwltmp,pondtmp,qbottmp, vt2
      real(8)   thetatmp(macp),htmp(macp),xd(6)


      save      sto,vsat,vt0,gwlt0, stocot1

! ----------------------------------------------------------------------
      comma = ','


      select case (task)
      case (1)

! === open output file =================================================
      filnam = trim(pathwork)//trim(outfil)//'.sto'
      sto = getun (20,90)
      call fopens(sto,filnam,'new','del')
      filtext = 'Storage and recharge output data '//                   &
     &'(q in cm/period; gwl in cm; sto in cm)'
      call writehead (sto,1,filnam,filtext,project)

! --- write header
      write (sto,10)
 10   format ('*',/,                                                    &
     &        'date,cumtime,qvss,qtr,qevs,qrun,qdr,qbot,qre,qv1m,'//    &
     &        'gwlt0,gwlt1,vt0,vt1,stocoav,stocot1')

      vsat = 0.0d0
      do nod = 1,numnod
         vsat = vsat + cofgen(2,nod) * dz(nod)
      end do

      return

      case (2)

! === write actual data ================================

! --- write header in case of new balance period
      if (flheader) write (sto,10)
! --- write actual data
      vt1 = 0.0d0
      do nod = 1,numnod
         if(nod.gt.1)then
           if(zbotcp(nod).le.-100.d0 .and. zbotcp(nod-1).gt.-100.d0) then
              nod1m = nod
           end if
         end if
         vt1 = vt1 + theta(nod) * dz(nod)
      end do
      vt1 = vsat - vt1
      qre = iqdra - iqbot
      qv1m = q(nod1m)
      gwlt1 = gwl
      if(dabs(gwlt1 - gwlt0) .gt.1.0d-6)then
         stocoav = - (vt1 - vt0) / (gwlt1 - gwlt0)
      else
         stocoav = 1.0d+6
      end if


      vt0 = vt1
      gwlt0 = gwlt1


! --- initilization
      xd(1) = q(1)
      xd(2) = iqrot
      xd(3) = ievap
      xd(4) = iruno
      xd(5) = iqdra
      xd(6) = iqbot

      swBotbtmp = swbotb
      qbottmp = qbot
      gwltmp = gwl
      pondtmp = pond
      do nod =1,numnod
         thetatmp(nod) = theta(nod)
         htmp(nod)     = h(nod)
      end do


      swbotb =1
      dgwl(1) = 1.0d0
      dgwl(2) = -1.0d0
! ---   determine SoilWater bottom boundary conditions

      do i=1,2

        gwlinp = gwltmp + dgwl(i)

!        call BoundBottom

        fldtreduce = .true.
        do while(fldtreduce)
           fldtreduce = .false.

! ---   calculate drainage fluxes
           if (fldrain)                           call Drainage
           if (.not.fldecdt .and. flSurfaceWater) call SurfaceWater(2)
           if (SwFrost.eq.1)                      call FrozenBounds

! ---   calculate SoilWater
           if (.not.fldecdt) then
! --- save state variables of time = t
              call SoilWaterStateVar(1)

! --- calculate new soil water state variables
              call headcalc

! ---   calculate surface water balace
              if (flSurfaceWater) call SurfaceWater(3)
           end if

! ---   update time variables and switches/flags
           if(fldecdt .or. (flMacroPore .and. FlDecMpRat))then
              call SoilWaterStateVar(2)
              call TimeControl(3)
              fldtreduce = .true.
           end if

        end do

        vt2 = 0.0d0
        do nod = 1,numnod
           vt2 = vt2 + theta(nod) * dz(nod)
        end do
        vair(i) = vsat - vt2

      end do


      stocot1 = - (vair(1) - vair(2))  / (dgwl(1) - dgwl(2))


      write (sto,20) date,comma,t1900,(comma,xd(i),i=1,6),comma,        &
     &               qre,comma,qv1m,comma,gwlt0,comma,gwlt1,comma,      &
     &               vt0,comma,vt1,comma,stocoav,comma,stocot1
20    format (a11,a1,f12.4,14(a1,f10.4))


      swbotb = swBotbtmp
      qbot = qbottmp
      gwl = gwltmp
      pond = pondtmp
      do nod =1,numnod
         theta(nod) = thetatmp(nod)
         h(nod)     = htmp(nod)
      end do


      return

      case (3)

! === close output file ===========================

! --- close snw file
      close (sto)

      case default
         call fatalerr ('OutputModFlow', 'Illegal value for Task')
      end select

      return
      end


! ----------------------------------------------------------------------
      subroutine capriseoutput(task)
! ----------------------------------------------------------------------
!     Date               : September 2016
!     Purpose            : open and write data to explore capillary rise to/from rootzone
! ----------------------------------------------------------------------
      use variables
      implicit none

! --- global variables ------------------
      integer task
! --- local variables ------------------
      integer   getun,sto
      character(len=300) filnam
      character(len=200) filtext
      character(len=1)   comma

      save      sto

! ----------------------------------------------------------------------
      comma = ','


      select case (task)
      case (1)

! === open output file =================================================
      filnam = trim(pathwork)//trim(outfil)//'.crz'
      sto = getun (20,90)
      call fopens(sto,filnam,'new','del')
      filtext = 'Capillary rise q to/from rootzone when roots '//      &
     &' are present (rd>0) (q in cm/period, pos=upward, neg=downward)'
      call writehead (sto,1,filnam,filtext,project)

! --- write header
      write (sto,10)
 10   format ('*',/,                                                    &
     &        'date,tcum,noddrz,rd,inq(noddrz)')
      return

      case (2)

! === write actual data ================================
      if(rd.gt.0.0d0) then
        write (sto,20) date,comma,tcum,comma,noddrz,comma,rd,comma,&
    &                  inq(noddrz+1)
 20     format (a11,a1,f12.4,a1,i4,a1,f10.4,a1,f20.12)
      endif

      return

      case (3)

! === close output file ===========================

! --- close snw file
      close (sto)

      case default
         call fatalerr ('capriseoutput', 'Illegal value for Task')
      end select

      return
      end

subroutine outrume (task)
use variables, only: tcum,outper,numnod,dz,theta,thetas,igrai,iruno,pond,gwl
implicit none

! global
integer, intent(in) :: task
! local
integer                       :: i, getun
integer,          save        :: LZnod             ! last node lying within LZ
integer,          save        :: LZnod2            ! last node lying within LZ2
real(8)                       :: sum, VT, WC
real(8),          parameter   :: LZ = 25.0d0       ! thickness of upper soil layer for air-filled pore volume will be calculated
real(8),          parameter   :: LZ2 = 10.0d0      ! thickness of upper soil layer for average water content
integer,          save        :: iunout2  !, iunout
!!!character(len=*), parameter   :: fout   = 'swap_rume.unf'
character(len=*), parameter   :: fout2  = 'swap_rume.csv'

logical, save :: Event
integer, save :: Nevent, Nrec, iDay, iDayOld, DayCum
real(8), save :: SumRain, SumRunOff, VTstart, GWLstart, hstart, RainRate, RunoffRate, DayOld, WCstart
real(8), save :: Tr, Ia, SS, hh, T, GG, P, Q, CN, S

select case (task)
case (1)
!!!   iunout  = getun (300, 900)
   iunout2 = getun (360, 900)
!!!   open (unit=iunout, file=fout, status='unknown', form='unformatted')
   open (unit=iunout2, file=fout2, status='unknown', form='formatted')
   ! write header
   write (iunout2,'(A)') 'DayCum (d), Day (d),Tr (d),P (mm),Q (mm),Ia (mm),SS (mm),CN (-),S (mm),h0 (mm),GWL (cm),WCini (cm3/cm3)'
   ! simple
   sum = 0.0d0
   do i = 1, numnod
      sum = sum + dz(i)
      if (sum > LZ) exit
   end do
   LZnod = i
   ! simple
   sum = 0.0d0
   do i = 1, numnod
      sum = sum + dz(i)
      if (sum > LZ2) exit
   end do
   LZnod2 = i-1

   Event     = .false.
   Nevent    = 0
   SumRain   = 0.0d0
   SumRunOff = 0.0d0
   VTstart   = 0.0d0
   GWLstart  = 0.0d0
   hstart    = 0.0d0
   DayOld    = 0.0d0
   iDayOld   = 0
   iDay      = 0
   Nrec      = 0

case (2)
   VT = 0.0d0
   do i = 1, LZnod
      VT = VT + dz(i)*(thetas(i)-theta(i))
   end do
   WC = 0.0d0
   do i = 1, LZnod2
      WC = WC + dz(i)*theta(i)
   end do
   WC = WC/LZ2

!  d, cm, cm/d, cm/d, cm, cm
!!!   write (iunout) real(tcum),real(VT),real(igrai/outper),real(iruno/outper),real(pond),real(gwl)
   RainRate   = igrai/outper
   RunoffRate = iruno/outper
   iDay       = int(tcum)
!!!   write (iunout,'(20F12.6)') real(tcum),real(VT)

   if (Nrec > 0) then
      if (.not. Event) then
         if (iDay > iDayOld) then  ! new day, no runoff occurring; restart summations
            SumRain   = RainRate*(tcum-iDay)
            SumRunOff = RunoffRate*(tcum-iDay) ! likely zero
            VTstart   = VT
            WCstart   = WC
            hstart    = Pond
            GWLstart  = GWL
            iDayOld   = iDay
         else if (RunoffRate > 0.0d0) then  ! start new runoff event!
            Event  = .true.
            Nevent = Nevent + 1
            Tr     = tcum-int(tcum)
            Ia     = SumRain*10.0d0 ! mm
            SS     = VTstart*10.0d0 ! mm
            hh     = hstart*10.0d0  ! mm
            T      = tcum
            GG     = GWLstart
            DayCum = int(tcum) + 1
         end if
         SumRain   = SumRain   + RainRate*(tcum-DayOld)
         SumRunoff = SumRunoff + RunoffRate*(tcum-DayOld)
      else
         if (RunoffRate < 1.0D-10) then  ! stop runoff event
            Event = .false.
            P     = SumRain*10.0 ! mm
            Q     = SumRunoff*10.0 ! mm
            CN    = 25400.0d0/((P-Ia)**2/Q - (P-Ia) + 254.0d0)
            S     = 25400.0d0/CN - 254.0d0
            write (iunout2,'(I12,11(A,F12.5))') DayCum, ',', T, ',', Tr, ',', P, ',', Q, ',', Ia, ',', SS, ',', CN, ',', S, ',', hh, ',', GG, ',', WCstart
         end if
         SumRain   = SumRain   + RainRate*(tcum-DayOld)
         SumRunoff = SumRunoff + RunoffRate*(tcum-DayOld)
      end if
      DayOld    = tcum

   else
      Nrec      = Nrec + 1
      SumRain   = RainRate*tcum
      SumRunOff = RunoffRate*tcum
      VTstart   = VT
      GWLstart  = GWL
      hstart    = Pond
      iDayOld   = iDay
      DayOld    = tcum
   end if

case (3)
!!!   close (iunout)
   close (iunout2)
case default
   call fatalerr ('outrume','Illegal task value')
end select

end subroutine outrume

subroutine csv_write (iTask)
! Routine designed for CSV output of user-selected vaiables (provided matchig defined variables in this routine).
! Currently only valid for time-dependent, single variables.
! Contains help routines: do_write_csv; check_list; remove_sqbr; det_node; Make_Header

! import global variables continaing possible output
use variables, only: pathwork, outfil, project, InList_csv,                   &
                     date, iptra, iqrot, iqreddry, iqredwet, iqredsol, gwl,   &
                     tsum, pgasspot, pgass, cwdmpot, cwdm, wsopot, wso,       &
                     laipot, lai, rdpot, rd,                                  &
                     tagppot, tagp, tagptpot, tagpt,                          &
                     cuptgrazpot, cuptgraz, plossdm, lossdm,                  &
                     h, theta, tsoil,                                         &
                     igrai, isnrai, igsnow, igird, iintc, irunon, iruno,      &
                     ipeva, ievap, iQMpOutDrRap, iqdra, iqbot, pond,          &
                     iqredfrs, wlvpot, wlv, wstpot, wst, wrtpot, wrt,         &
                     dwso,dwlv,dwlvpot,dwst,dwstpot,dwrt,dwrtpot,             &
                     dvs, ch, cf, K, cml, cmsy, c_top, flprintshort, t1900,   &
                     imsqprec, imsqirrig, imsqbot, imsqdra, imdectot,         &
                     imrottot, sampro, wc10, Runoff_CN, iqtdo, iqtup,         &
                     iqssdi, volact, ssnow, volini, pondini, irunocn, isubl

implicit none
! global
integer,          intent(in)        :: iTask

! local (some need to be saved)

! allowed variable names: number of items in Allowed must be exactly equal to Mlist
! Since Forcheck reports an error when elements have different number of characters, all now have same length
! Must be UPPERCASE
integer,                               parameter   :: Mlist = 73
character(len=*),    dimension(Mlist), parameter   :: Allowed = (/   &
'RAIN    ','SNOW    ','IRRIG   ','INTERC  ','RUNON   ','RUNOFF  ',   &
'EPOT    ','EACT    ','DRAINAGE','QBOTTOM ','GWL     ','POND    ',   &
'TPOT    ','TACT    ','TREDDRY ','TREDWET ','TREDSOL ','TREDFRS ',   &
'TSUM    ','DVS     ','PGASSPOT','PGASS   ','CPWDM   ','CWDM    ',   &
'CPWSO   ','CWSO    ','PWLV    ','WLV     ','PWST    ','WST     ',   &
'PWRT    ','WRT     ','DWSO    ','DWLV    ','DWLVPOT ',              &
'DWST    ','DWSTPOT ','DWRT    ','DWRTPOT ','HEIGHT  ','CRPFAC  ',   &
'LAIPOT  ','LAI     ','RDPOT   ','RD      ','PGRASSDM','GRASSDM ',   &
'PMOWDM  ','MOWDM   ','PGRAZDM ','GRAZDM  ','PLOSSDM ','LOSSDM  ',   &
'SQPREC  ','SQIRRIG ','SQBOT   ','SQDRA   ','DECTOT  ','ROTTOT  ',   &
'SAMPRO  ','WC10    ','RUNOFFCN','QTOPIN  ','QTOPOUT ','BALDEV  ',   &
'QSSDI   ','H       ','WC      ','TEMP    ','K       ','CONC    ',   &
'CONCADS ','O2TOP   '                                                &
/)
character(len=*),    dimension(Mlist), parameter   :: Units = (/           &
'(cm)     ','(cm)     ','(cm)     ','(cm)     ','(cm)     ','(cm)     ',   &
'(cm)     ','(cm)     ','(cm)     ','(cm)     ','(cm)     ','(cm)     ',   &
'(cm)     ','(cm)     ','(cm)     ','(cm)     ','(cm)     ','(cm)     ',   &
'(deg C)  ','(-)      ','(kgch/ha)','(kgch/ha)','(kg/ha)  ','(kg/ha)  ',   &
'(kg/ha)  ','(kg/ha)  ','(kg/ha)  ','(kg/ha)  ','(kg/ha)  ','(kg/ha)  ',   &
'(kg/ha)  ','(kg/ha)  ','(kg/ha)  ','(kg/ha)  ','(kg/ha)  ',               &
'(kg/ha)  ','(kg/ha)  ','(kg/ha)  ','(kg/ha)  ','(cm)     ','(-)      ',   &
'(m2/m2)  ','(m2/m2)  ','(cm)     ','(cm)     ','(kg/ha)  ','(kg/ha)  ',   &
'(kg/ha)  ','(kg/ha)  ','(kg/ha)  ','(kg/ha)  ','(kg/ha)  ','(kg/ha)  ',   &
'(g/cm2)  ','(g/cm2)  ','(g/cm2)  ','(g/cm2)  ','(g/cm2)  ','(g/cm2)  ',   &
'(g/cm2)  ','(cm3/cm3)','(cm)     ','(cm)     ','(cm)     ','(cm)     ',   &
'(cm)     ','(cm)     ','(cm3/cm3)','(deg C)  ','(cm/d)   ','(g/cm3 w)',   &
'(g/cm3)  ','(kg/m3)  '                                                    &
/)

! to add (.INC): dstorage?, baldev?
! to add (.WBA): wbal?
! to add (.STR):
! to add (.CRP): DWLVCROP, DWLVSOIL, DWST, DWRT, DWSO, HARLOSORM

integer,             dimension(Mlist), save        :: iCSV, iPOS

integer,             parameter                     :: ilw = Mlist
integer,             dimension(ilw)                :: iWbeg, iWend
integer,                               save        :: iuncsv, Nvars
integer                                            :: i, j
character(len=*),    parameter                     :: comma   = ','
character(len=160)                                 :: filtext
character(len=300)                                 :: filnam
character(len=2)                                   :: cval1, cval2
character(len=20),                     save        :: form_rea_F1, form_rea_E1
character(len=1024)                                :: InList_csv_x, Header, HeaderUnits
character(len=20),   dimension(ilw)                :: listVars
character(len=19)                                  :: datexti

! Specific nodal information
integer,             parameter                     :: Mnodes = 10
integer,                                save       :: NumNodes_H, NumNodes_WC, NumNodes_TEMP, NumNodes_K, NumNodes_CONC, NumNodes_CONCADS, NumNodes_O2TOP
integer,             dimension(Mnodes), save       :: Nodes_H, Nodes_WC, Nodes_TEMP, Nodes_K, Nodes_CONC, Nodes_CONCADS, Nodes_O2TOP
character(len=20),   dimension(Mnodes)             :: Str_nodes_H, Str_nodes_WC, Str_nodes_TEMP, Str_nodes_K, Str_nodes_CONC, Str_nodes_CONCADS, Str_nodes_O2TOP

! when to automatically swith from F to E formatting
integer,             parameter                     :: num_d = 5         ! # of decimals; later: user input?
integer,             parameter                     :: num_w = num_d+7   ! total width of format, for E-formatting: 7 positions are needed for "-x."at start and "E+00" at end
integer,             parameter                     :: expo = 4
real(8),             parameter                     :: t1 = 1.0d0/(10.0d0**expo)
real(8),             parameter                     :: t2 = 10.0d0**expo

! other
real(8), save :: dstor, baldev, volold, snowold, pondold

! functions
integer                                            :: getun

select case (iTask)
case (1)

   ! some basic info
   !  Note: field width of zero in I and F edit descriptors is allowed as of Fortran95 to ensure as little space usage
   !        in output file as possible (Metcalf et al., 2004, Fortran 95/2003 explained, Oxford Univ. Prrss, p. 199)
   write (cval1,'(I0)') num_d
   write (cval2,'(I0)') num_w
   form_rea_F1 = '(A1,F0.' // trim(cval1) // ')'
   form_rea_E1 = '(A1,1P,E' // trim(cval2) // '.' // trim(cval1) // ')'

   ! initialize
   NumNodes_H = 0; NumNodes_WC = 0; NumNodes_K = 0; NumNodes_TEMP = 0
   NumNodes_CONC = 0; NumNodes_CONCADS = 0; NumNodes_O2TOP = 0

   ! make user-supplied list UPPERCASE
   call upperc (InList_csv)

   ! reduce list: remove information between square brackets (if present) (square brackets don't appear in Allowed)
   call remove_sqbr (InList_csv, InList_csv_x)

   ! Nvars en ListVars
   call words (InList_csv_x, ilw, ', ', iWbeg, iWend, Nvars)
   do i = 1, Nvars
      ListVars(i) = InList_csv_x(iWbeg(i):iWend(i))
   end do

   ! check if ListVars contains valid data; determine iCSV and iPOS
   call check_list ()

   ! sort position in inList_csv in same way as in Allowed; ListVars is changed accordingly
   ! this needed since sequence in output columns are fixed by appearance in Allowed
   call sort_list (InList_csv_x, ListVars, Nvars)

   ! open file for output; existing file will be overwritten; formatted output
   !    to do: write some basic info at the top of the output file?
   iuncsv = getun (500, 900)
   filnam = trim(pathwork)//trim(outfil)//'_output.csv'
   open (unit=iuncsv, file=filnam, status='replace', form='formatted')

   ! column header
   filtext = 'specified output data of SWAP'
   call writehead (iuncsv,1,filnam,filtext,project)

   Header = 'DATE'
   call make_header (ListVars, Nvars, Header)
   HeaderUnits = '*'
   call make_headerunits (HeaderUnits)
   write (iuncsv,'(A)') trim(HeaderUnits)
   write (iuncsv,'(A)') trim(Header)

   volold  = volini
   snowold = 0.0d0
   pondold = pondini
   
case (2)

   ! date and time
   if (.not. flprintshort) then
      write (iuncsv,'(2A)',advance='no') trim(date)
   else
      ! determine date-time
      call dtdpst ('year-month-day hour:minute:seconds',t1900,datexti)
      write (iuncsv,'(2A)',advance='no') trim(datexti)
   end if

   ! all other variables
   ! programmer is responsible for correct correspondence between Names (and their position) in Allowed and
   ! actual SWAP variables as used below

   if (iCSV(1)  == 1) call do_write_csv (igrai+isnrai)         ! RAIN (PRECIPITATION)
   if (iCSV(2)  == 1) call do_write_csv (igsnow)               ! SNOW
   if (iCSV(3)  == 1) call do_write_csv (igird)                ! IRRIGATION
   if (iCSV(4)  == 1) call do_write_csv (iintc)                ! INTERCEPTION
   if (iCSV(5)  == 1) call do_write_csv (irunon)               ! RUNON
   if (iCSV(6)  == 1) call do_write_csv (iruno)                ! RUNOFF
   if (iCSV(7)  == 1) call do_write_csv (ipeva)                ! EPOT
   if (iCSV(8)  == 1) call do_write_csv (ievap)                ! EACT
   if (iCSV(9)  == 1) call do_write_csv (iQMpOutDrRap+iqdra)   ! QDRAIN
   if (iCSV(10) == 1) call do_write_csv (iqbot)                ! QBOTTOM
   if (iCSV(11) == 1) call do_write_csv (gwl)                  ! GWL
   if (iCSV(12) == 1) call do_write_csv (pond)                 ! POND

   if (iCSV(13) == 1) call do_write_csv (iptra)                ! TPOT
   if (iCSV(14) == 1) call do_write_csv (iqrot)                ! TACT
   if (iCSV(15) == 1) call do_write_csv (iqreddry)             ! TREDDRY
   if (iCSV(16) == 1) call do_write_csv (iqredwet)             ! TREDWET
   if (iCSV(17) == 1) call do_write_csv (iqredsol)             ! TREDSOL
   if (iCSV(18) == 1) call do_write_csv (iqredfrs)             ! TREDFRS

   ! CROP OUTPUT
   if (iCSV(19) == 1) call do_write_csv (tsum)                 ! TSUM
   if (iCSV(20) == 1) call do_write_csv (dvs)                  ! DVS

   if (iCSV(21) == 1) call do_write_csv (pgasspot)             ! PGASSPOT
   if (iCSV(22) == 1) call do_write_csv (pgass)                ! PGASS
   
   if (iCSV(23) == 1) call do_write_csv (cwdmpot)              ! CPWDM
   if (iCSV(24) == 1) call do_write_csv (cwdm)                 ! CWDM
   if (iCSV(25) == 1) call do_write_csv (wsopot)               ! CPWSO
   if (iCSV(26) == 1) call do_write_csv (wso)                  ! CWSO
   if (iCSV(27) == 1) call do_write_csv (wlvpot)               ! PWLV
   if (iCSV(28) == 1) call do_write_csv (wlv)                  ! WLV
   if (iCSV(29) == 1) call do_write_csv (wstpot)               ! PWST
   if (iCSV(30) == 1) call do_write_csv (wst)                  ! WST
   if (iCSV(31) == 1) call do_write_csv (wrtpot)               ! PWRT
   if (iCSV(32) == 1) call do_write_csv (wrt)                  ! WRT

   if (iCSV(33) == 1) call do_write_csv (dwso)                 ! DWSO
   if (iCSV(34) == 1) call do_write_csv (dwlv)                 ! DWLV
   if (iCSV(35) == 1) call do_write_csv (dwlvpot)              ! DWLVPOT
   if (iCSV(36) == 1) call do_write_csv (dwst)                 ! DWST
   if (iCSV(37) == 1) call do_write_csv (dwstpot)              ! DWSTPOT
   if (iCSV(38) == 1) call do_write_csv (dwrt)                 ! DWRT
   if (iCSV(39) == 1) call do_write_csv (dwrtpot)              ! DWRTPOT
   
   if (iCSV(40) == 1) call do_write_csv (ch)                   ! HEIGHT
   if (iCSV(41) == 1) call do_write_csv (cf)                   ! CRPFAC

   if (iCSV(42) == 1) call do_write_csv (laipot)               ! LAIPOT
   if (iCSV(43) == 1) call do_write_csv (lai)                  ! LAI
   if (iCSV(44) == 1) call do_write_csv (rdpot)                ! RDPOT
   if (iCSV(45) == 1) call do_write_csv (rd)                   ! RD

   ! for grass
   ! PGRASSDM   GRASSDM   PMOWDM   MOWDM   PGRAZDM     GRAZDM
   ! tagppot    tagp      tagptpot tagpt   cuptgrazpot cuptgraz
   if (iCSV(46) == 1) call do_write_csv (tagppot)              ! PGRASSDM
   if (iCSV(47) == 1) call do_write_csv (tagp)                 ! GRASSDM
   if (iCSV(48) == 1) call do_write_csv (tagptpot)             ! PMOWDM
   if (iCSV(49) == 1) call do_write_csv (tagpt)                ! MOWDM
   if (iCSV(50) == 1) call do_write_csv (cuptgrazpot)          ! PGRAZDM
   if (iCSV(51) == 1) call do_write_csv (cuptgraz)             ! GRAZDM
   if (iCSV(52) == 1) call do_write_csv (plossdm)              ! PLOSSDM
   if (iCSV(53) == 1) call do_write_csv (lossdm)               ! LOSSDM

   ! solute: imsqprec, imsqirrig, imsqbot, imsqdra, imdectot, imrottot, sampro
   if (iCSV(54) == 1) call do_write_csv (imsqprec)             ! SQPREC
   if (iCSV(55) == 1) call do_write_csv (imsqirrig)            ! SQIRRIG
   if (iCSV(56) == 1) call do_write_csv (imsqbot)              ! SQBOT
   if (iCSV(57) == 1) call do_write_csv (imsqdra)              ! SQDRA
   if (iCSV(58) == 1) call do_write_csv (imdectot)             ! DECTOT
   if (iCSV(59) == 1) call do_write_csv (imrottot)             ! ROTTOT
   if (iCSV(60) == 1) call do_write_csv (sampro)               ! SAMPRO

   ! surface: water content, runoff, net inflow, net outflow, max inf rate
   if (iCSV(61) == 1) call do_write_csv (wc10)                 ! WC10
   if (iCSV(62) == 1) call do_write_csv (Runoff_CN)            ! RUNOFF_CN
   if (iCSV(63) == 1) call do_write_csv (iqtdo)                ! QTOPIN
   if (iCSV(64) == 1) call do_write_csv (iqtup)                ! QTOPOUT

   dstor  = (volact + pond + ssnow) - (VolOld + PondOld + SnowOld)
   baldev = (igrai+isnrai+igsnow+igird+irunon) - dstor -             &
            (iintc+iruno+irunoCN+iqrot+ievap+isubl+iQMpOutDrRap+iqdra+(-1.0d0*iqbot) - iqssdi)
   if (iCSV(65) == 1) call do_write_csv (baldev)               ! BALDEV
      VolOld = volact
      PondOld = pond
      SnowOld = ssnow


   ! subsurface drip irrigation
   if (iCSV(66) == 1) call do_write_csv (iqssdi)               ! QSSDI

   ! others: H, WC, TEMP for specific nodes
   if (iCSV(67) == 1) then
      do j = 1, NumNodes_H
         call do_write_csv (h(Nodes_H(j)))
      end do
   end if
   if (iCSV(68) == 1) then
      do j = 1, NumNodes_WC
         call do_write_csv (theta(Nodes_WC(j)))
      end do
   end if
   if (iCSV(69) == 1) then
      do j = 1, NumNodes_TEMP
         call do_write_csv (tsoil(Nodes_TEMP(j)))
      end do
   end if
   if (iCSV(70) == 1) then
      do j = 1, NumNodes_K
         call do_write_csv (K(Nodes_K(j)))
      end do
   end if
   if (iCSV(71) == 1) then
      do j = 1, NumNodes_CONC
         call do_write_csv (cml(Nodes_CONC(j)))
      end do
   end if
   if (iCSV(72) == 1) then
      do j = 1, NumNodes_CONCADS
         call do_write_csv (cmsy(Nodes_CONCADS(j)))
      end do
   end if
   if (iCSV(73) == 1) then
      do j = 1, NumNodes_O2TOP
         call do_write_csv (c_top(Nodes_O2TOP(j)))
      end do
   end if

   ! finalize record (advance to next line)
   write (iuncsv,*)

case (3)

   close (unit=iuncsv)

case default

   call fatalerr ('csv_write','Illegal iTask value; range allowed: [1-3]')

end select

return

contains


subroutine do_write_csv (var)
implicit none
real(8),          intent(in)  :: var

if ( (var > 0.0d0 .AND. var < t1) .OR. (var < 0.0d0 .AND. var > -t1) ) then
   write (iuncsv,form_rea_E1,advance='no') comma, real(var)
else if ( (var > t2) .OR. (var < -t2) ) then
   write (iuncsv,form_rea_E1,advance='no') comma, real(var)
else
   write (iuncsv,form_rea_F1,advance='no') comma, real(var)
end if

end subroutine do_write_csv

subroutine check_list ()
implicit none

integer :: i, j
logical :: isthere

iCSV = 0
iPOS = 0
do i = 1, Nvars
   isthere = .false.
   do j = 1, Mlist
!      isthere = index(ListVars(i), trim(Allowed(j))) > 0  werkt niet bijv bij combinatie RAIN en DRAINAGE
      isthere = trim(ListVars(i)) == trim(Allowed(j))
      if (isthere) then
         iCSV(j) = 1    ! indicator for Yes/No for writing
         if (iPOS(j) > 0) call fatalerr ('check_list', 'Double entries in inList_csv are not allowed')
         iPOS(j) = i    ! help vector to be used later for sorting
         exit
      end if
   end do
   ! error if not there; or should we write a warning?
   if (.not. isthere) call fatalerr ('check_list','Illegal variable name in inList_csv for write_csv output')
end do

end subroutine check_list

subroutine sort_list (InList_csv, ListVars, Nvars)

implicit none
! global; NB: on output both variables are arranged in same sorted order as Allowed
integer,                            intent(in)     :: Nvars
character(len=*),                   intent(out)    :: InList_csv
character(len=*), dimension(Nvars), intent(inout)  :: ListVars
! local
integer              :: il, i, j
character(len=1024)  :: temp

il = 0
! iPOS contains for each of the Allowed names the position in the user-supplied list ListVars,
! and adds this to a temporary string
do j = 1, Mlist
   if (iPOS(j) > 0) then
      if (il == 0) then
         temp = trim(listVars(iPOS(j))); il = len_trim(temp)
      else
         temp = temp(1:il) // comma // trim(listVars(iPOS(j))); il = len_trim(temp)
      end if
   end if
end do

! temp now contains the user-supplied names in the same order of appearance is in Allowed
! decompose temp into ListVars (overwrite)
call words (temp, ilw, comma, iWbeg, iWend, Nvars)
do i = 1, Nvars
   ListVars(i) = temp(iWbeg(i):iWend(i))
end do

! on return: InList_csv becomes equal to sorted string temp
InList_csv = trim(temp)

end subroutine sort_list

subroutine remove_sqbr (Inlist_csv, InList_csv_x)
implicit none
! global
character(len=*), intent(in)  :: Inlist_csv
character(len=*), intent(out) :: Inlist_csv_x
! local

integer :: ipos1, ipos2
character(len=80) :: H_list, WC_list, TEMP_list, K_list, C_list, CA_list, O2_list

! start with output = input; later output string may be changed
Inlist_csv_x = trim(Inlist_csv)

1 continue
ipos1 = index(InList_csv_x,"[")
ipos2 = index(InList_csv_x,"]")
if (ipos1 == 0 .AND. ipos2 == 0) then
   ! no more square brackets present
   return
else
   if (ipos2 < ipos1 .OR. (ipos2 > 0 .AND. ipos1 == 0)) then
      call fatalerr ('remove_sqbr','Found ] before [')
   else if (index(InList_csv_x(ipos1:ipos2),"H") > 0 .OR. index(InList_csv_x(ipos1:ipos2),"WC") > 0 .OR. index(InList_csv_x(ipos1:ipos2),"TEMP") > 0 .OR.   &
            index(InList_csv_x(ipos1:ipos2),"K") > 0 .OR. index(InList_csv_x(ipos1:ipos2),"CONC") > 0 .OR. index(InList_csv_x(ipos1:ipos2),"CONCADS") > 0 .OR. &
            index(InList_csv_x(ipos1:ipos2),"O2TOP") > 0) then
      call fatalerr ('remove_sqbr','Possible mismatch in arranging information between [ and ]')
   else
      if (Inlist_csv_x(ipos1-1:ipos1-1) == 'H') then
         H_list = Inlist_csv_x(ipos1+1:ipos2-1)
         call det_node (H_list, Nodes_H, Mnodes, NumNodes_H, Str_nodes_H)
      else if (Inlist_csv_x(ipos1-2:ipos1-1) == 'WC') then
         WC_list = Inlist_csv_x(ipos1+1:ipos2-1)
         call det_node (WC_list, Nodes_WC, Mnodes, NumNodes_WC, Str_nodes_WC)
      else if (Inlist_csv_x(ipos1-4:ipos1-1) == 'TEMP') then
         TEMP_list = Inlist_csv_x(ipos1+1:ipos2-1)
         call det_node (TEMP_list, Nodes_TEMP, Mnodes, NumNodes_TEMP, Str_nodes_TEMP)
      else if (Inlist_csv_x(ipos1-1:ipos1-1) == 'K') then
         K_list = Inlist_csv_x(ipos1+1:ipos2-1)
         call det_node (K_list, Nodes_K, Mnodes, NumNodes_K, Str_nodes_K)
      else if (Inlist_csv_x(ipos1-4:ipos1-1) == 'CONC') then
         C_list = Inlist_csv_x(ipos1+1:ipos2-1)
         call det_node (C_list, Nodes_CONC, Mnodes, NumNodes_CONC, Str_nodes_CONC)
      else if (Inlist_csv_x(ipos1-7:ipos1-1) == 'CONCADS') then
         CA_list = Inlist_csv_x(ipos1+1:ipos2-1)
         call det_node (CA_list, Nodes_CONCADS, Mnodes, NumNodes_CONCADS, Str_nodes_CONCADS)
      else if (Inlist_csv_x(ipos1-5:ipos1-1) == 'O2TOP') then
         O2_list = Inlist_csv_x(ipos1+1:ipos2-1)
         call det_node (O2_list, Nodes_O2TOP, Mnodes, NumNodes_O2TOP, Str_nodes_O2TOP)
      else
         call fatalerr ('remove_sqbr','Found [] with illegal preceding variable name')
      end if

      Inlist_csv_x(ipos1:ipos2) = " "
   end if
   go to 1
end if

end subroutine remove_sqbr

subroutine det_node (String, Nodes, Mnodes, NumNodes, Str_nodes)
use variables, only: numnod, zbotcp
implicit none
! global
integer,                             intent(in)    :: Mnodes
integer,                             intent(out)   :: NumNodes
integer,          dimension(Mnodes), intent(out)   :: Nodes
character(len=*),                    intent(in)    :: String
character(len=*), dimension(Mnodes), intent(out)   :: Str_nodes
! local
integer, parameter :: ilw = 15  ! should be > Mnodes
integer, dimension(ilw) :: iw1, iw2
integer :: ifnd, i
real, dimension(ilw) :: dum

if (Mnodes > ilw) call fatalerr ('det_node', 'Mnodes > ilw')

call words (String, ilw, ', ', iw1, iw2, ifnd)
if (ifnd > Mnodes) call fatalerr ('det_node', 'No more than Mnodes allowed for nodal output')
NumNodes = ifnd
Nodes = 0
dum = 0.0
Str_nodes = ""
do i = 1, NumNodes
   read (String(iw1(i):iw2(i)),*) dum(i) ! Nodes(i)
   if (dum(i) < 0.0) then
      Str_nodes(i) = String(iw1(i):iw2(i))
   else
      write(Str_nodes(i),'(I0)') nint(dum(i))
   end if
end do

if (minval(dum(1:NumNodes))*maxval(dum(1:NumNodes)) < 0) call fatalerr ('det_node', 'All entries between [] must be either >0 or <0; not mixed')

do i = 1, NumNodes
   if (dum(i) > 0) then
      Nodes(i) = nint(dum(i))    ! any floating point positive number is rounded to nearest integer; no warning/error
      if (Nodes(i) > numnod) call fatalerr ('det_node', 'Entry between [] cannot exceed numnod')
   else
      do j = 1, numnod
         if (dum(i) > zbotcp(j)) then
            Nodes(i) = j   ! replace by correct node number
            exit
         end if
      end do
   end if
end do

end subroutine det_node

subroutine make_headerunits (HeaderUnits)
implicit none
! global
character(len=*),                   intent(inout)  :: HeaderUnits
! local
integer :: i, il

il = len_trim(HeaderUnits)
do i = 1, 66
   if (iCSV(i) == 1) then
      HeaderUnits = HeaderUnits(1:il) // "," // trim(Units(i)); il = len_trim(HeaderUnits)
   end if
end do

! others: H, WC, TEMP for specific nodes
if (iCSV(67) == 1) then
   i = 67
   do j = 1, NumNodes_H
      HeaderUnits = HeaderUnits(1:il) // "," // trim(Units(i)); il = len_trim(HeaderUnits)
   end do
end if
if (iCSV(68) == 1) then
   i = 68
   do j = 1, NumNodes_WC
      HeaderUnits = HeaderUnits(1:il) // "," // trim(Units(i)); il = len_trim(HeaderUnits)
   end do
end if
if (iCSV(69) == 1) then
   i = 69
   do j = 1, NumNodes_TEMP
      HeaderUnits = HeaderUnits(1:il) // "," // trim(Units(i)); il = len_trim(HeaderUnits)
   end do
end if
if (iCSV(70) == 1) then
   i = 70
   do j = 1, NumNodes_K
      HeaderUnits = HeaderUnits(1:il) // "," // trim(Units(i)); il = len_trim(HeaderUnits)
   end do
end if
if (iCSV(71) == 1) then
   i = 71
   do j = 1, NumNodes_CONC
      HeaderUnits = HeaderUnits(1:il) // "," // trim(Units(i)); il = len_trim(HeaderUnits)
   end do
end if
if (iCSV(72) == 1) then
   i = 72
   do j = 1, NumNodes_CONCADS
      HeaderUnits = HeaderUnits(1:il) // "," // trim(Units(i)); il = len_trim(HeaderUnits)
   end do
end if
if (iCSV(73) == 1) then
   i = 73
   do j = 1, NumNodes_O2TOP
      HeaderUnits = HeaderUnits(1:il) // "," // trim(Units(i)); il = len_trim(HeaderUnits)
   end do
end if

end subroutine make_headerunits

subroutine Make_Header (ListVars, Nvars, Header)
implicit none
! global
integer,                            intent(in)     :: Nvars
character(len=*), dimension(Nvars), intent(in)     :: ListVars
character(len=*),                   intent(inout)  :: Header
! local
integer :: i, il

il = len_trim(Header)
do i = 1, Nvars
   if (trim(ListVars(i)) == "H") then
      do j = 1, NumNodes_H
         Header = Header(1:il) // "," // trim(ListVars(i)) // "[" // trim(Str_nodes_H(j)) // "]"; il = len_trim(Header)
      end do
   else if (trim(ListVars(i)) == "WC") then
      do j = 1, NumNodes_WC
         Header = Header(1:il) // "," // trim(ListVars(i)) // "[" // trim(Str_nodes_WC(j)) // "]"; il = len_trim(Header)
      end do
   else if (trim(ListVars(i)) == "TEMP") then
      do j = 1, NumNodes_TEMP
         Header = Header(1:il) // "," // trim(ListVars(i)) // "[" // trim(Str_nodes_TEMP(j)) // "]"; il = len_trim(Header)
      end do
   else if (trim(ListVars(i)) == "K") then
      do j = 1, NumNodes_K
         Header = Header(1:il) // "," // trim(ListVars(i)) // "[" // trim(Str_nodes_K(j)) // "]"; il = len_trim(Header)
      end do
   else if (trim(ListVars(i)) == "CONC") then
      do j = 1, NumNodes_CONC
         Header = Header(1:il) // "," // trim(ListVars(i)) // "[" // trim(Str_nodes_CONC(j)) // "]"; il = len_trim(Header)
      end do
   else if (trim(ListVars(i)) == "CONCADS") then
      do j = 1, NumNodes_CONCADS
         Header = Header(1:il) // "," // trim(ListVars(i)) // "[" // trim(Str_nodes_CONCADS(j)) // "]"; il = len_trim(Header)
      end do
   else if (trim(ListVars(i)) == "O2TOP") then
      do j = 1, NumNodes_O2TOP
         Header = Header(1:il) // "," // trim(ListVars(i)) // "[" // trim(Str_nodes_O2TOP(j)) // "]"; il = len_trim(Header)
      end do
   else
      Header = Header(1:il) // "," // trim(ListVars(i)); il = len_trim(Header)
   end if
end do

end subroutine Make_Header

!subroutine det_jnodes(jnodes)
!implicit none
!integer, intent(out) :: jnodes
!jnodes = 0
!if (NumNodes_H       > 0) jnodes = jnodes + (NumNodes_H-1)
!if (NumNodes_WC      > 0) jnodes = jnodes + (NumNodes_WC-1)
!if (NumNodes_K       > 0) jnodes = jnodes + (NumNodes_K-1)
!if (NumNodes_TEMP    > 0) jnodes = jnodes + (NumNodes_TEMP-1)
!if (NumNodes_CONC    > 0) jnodes = jnodes + (NumNodes_CONC-1)
!if (NumNodes_CONCADS > 0) jnodes = jnodes + (NumNodes_CONCADS-1)
!if (NumNodes_O2TOP   > 0) jnodes = jnodes + (NumNodes_O2TOP-1)
!
!end subroutine det_jnodes


end subroutine csv_write

subroutine csv_write_tz (iTask)
! Routine designed for CSV output of user-selected vaiables (provided matching defined variables in this routine).
! Specifically for selected time-depth variables
! Contains help routines: do_write_csv; check_list; remove_sqbr; det_node; Make_Header

! import global variables contianing possible output
use variables, only: pathwork, outfil, project, InList_csv_tz, numnod, z, flprintshort, date, t1900, h, theta, tsoil, K, cml, cmsy, c_top, HEACAP, HEACON

implicit none
! global
integer,          intent(in)        :: iTask

! local (some need to be saved)

! allowed variable names: number of items in Allowed must be exactly equal to Mlist
! Since Forcheck reports an error when elements have different number of characters, all now have same length
! Must be UPPERCASE
integer,                               parameter   :: Mlist   = 9
character(len=*),    dimension(Mlist), parameter   :: Allowed =  (/'H      ',  'WC     ',  'TEMP   ',  'K      ',  'CONC   ',  'CONCADS',  'O2TOP  ',  'HEACAP ',  'HEACON '  /)
character(len=*),    dimension(Mlist), parameter   :: Units   =  (/'(cm)     ','(cm3/cm3)','(deg C)  ','(cm/d)   ','(g/cm3 w)','(g/cm3)  ','(kg/m3)  ','(J/m3/K) ','(W/m/K)  '/)

integer,             dimension(Mlist), save        :: iCSV, iPOS

integer,             parameter                     :: ilw = Mlist
integer,             dimension(ilw)                :: iWbeg, iWend
integer,                               save        :: iuncsv, Nvars
integer                                            :: i, j
character(len=*),    parameter                     :: comma   = ','
character(len=300)                                 :: filnam
character(len=160)                                 :: filtext
character(len=2)                                   :: cval1, cval2
character(len=20),                     save        :: formZ, form_rea_E
character(len=1024)                                :: Header, HeaderUnits
character(len=20),   dimension(ilw)                :: listVars
character(len=19)                                  :: datexti

! when to automatically swith from F to E formatting
integer,             parameter                     :: num_d = 5         ! # of decimals; later: user input?
integer,             parameter                     :: num_w = num_d+7   ! total width of format, for E-formatting: 7 positions are needed for "-x."at start and "E+00" at end

! functions
integer                                            :: getun

select case (iTask)
case (1)

   ! some basic info
   !  Note: field width of zero in I and F edit descriptors is allowed as of Fortran95 to ensure as little space usage
   !        in output file as possible (Metcalf et al., 2004, Fortran 95/2003 explained, Oxford Univ. Press, p. 199)
   write (cval1,'(I0)') num_d
   write (cval2,'(I0)') num_w
   form_rea_E = '(A1,1P,E' // trim(cval2) // '.' // trim(cval1) // ')'
   formZ = '(A1,F0.3)'

   ! counter for number of records
   !Nlines = 0

   ! make user-supplied list UPPERCASE
   call upperc (InList_csv_tz)

   ! Nvars en ListVars
   call words (InList_csv_tz, ilw, ', ', iWbeg, iWend, Nvars)
   do i = 1, Nvars
      ListVars(i) = InList_csv_tz(iWbeg(i):iWend(i))
   end do

   ! check if ListVars contains valid data; determine iCSV and iPOS
   call check_list_tz ()

   ! sort position in inList_csv in same way as in Allowed; ListVars is changed accordingly
   ! this needed since sequence in output columns are fixed by appearance in Allowed
   call sort_list_tz (InList_csv_tz, ListVars, Nvars)

   ! open file for output; existing file will be overwritten; formatted output
   !    to do: write some basic info at the top of the output file?
   iuncsv = getun (500, 900)
   filnam = trim(pathwork)//trim(outfil)//'_output_tz.csv'
   open (unit=iuncsv, file=filnam, status='replace', form='formatted')

   ! column header
   Header = 'DATE,DEPTH'
   HeaderUnits = '*,(cm)'
   call make_header_tz (ListVars, Nvars, Header)
   call make_headerunits_tz (HeaderUnits)
   filtext = 'specified output data of SWAP'
   call writehead (iuncsv,1,filnam,filtext,project)
   write (iuncsv,'(A)') trim(HeaderUnits)
   write (iuncsv,'(A)') trim(Header)

case (2)

  do j = 1, NumNod

    ! date and time
    if (.not. flprintshort) then
       write (iuncsv,'(2A)',advance='no') trim(date)
    else
       ! determine date-time
       call dtdpst ('year-month-day hour:minute:seconds',t1900,datexti)
       write (iuncsv,'(2A)',advance='no') trim(datexti)
    end if

    ! depth
    write (iuncsv,formZ,advance='no') comma, z(j)

    ! all other variables
    ! programmer is responsible for correct correspondence between Names (and their position) in Allowed and
    ! actual SWAP variables as used below

    if (iCSV(1)  == 1) call do_write_csv_tz (h(j))
    if (iCSV(2)  == 1) call do_write_csv_tz (theta(j))
    if (iCSV(3)  == 1) call do_write_csv_tz (tsoil(j))
    if (iCSV(4)  == 1) call do_write_csv_tz (k(j))
    if (iCSV(5)  == 1) call do_write_csv_tz (cml(j))
    if (iCSV(6)  == 1) call do_write_csv_tz (cmsy(j))
    if (iCSV(7)  == 1) call do_write_csv_tz (c_top(j))
    if (iCSV(8)  == 1) call do_write_csv_tz (HEACAP(j)/1.0d-6)    ! from J/cm3/K  to J/m3/K
    if (iCSV(9)  == 1) call do_write_csv_tz (HEACON(j)/864.0d0)   ! from J/cm/K/d to W/m/K

    ! finalize record (advance to next line)
    write (iuncsv,*)

  end do

case (3)

   close (unit=iuncsv)

case default

    call fatalerr ('csv_write_tz','Illegal iTask value; range allowed: [1-3]')

end select

return

contains

subroutine do_write_csv_tz (var)
implicit none
real(8),          intent(in)  :: var

write (iuncsv,form_rea_E,advance='no') comma, real(var)

end subroutine do_write_csv_tz

subroutine check_list_tz ()
implicit none

integer :: i, j
logical :: isthere

iCSV = 0
iPOS = 0
do i = 1, Nvars
   isthere = .false.
   do j = 1, Mlist
      isthere = trim(ListVars(i)) == trim(Allowed(j))
      if (isthere) then
         iCSV(j) = 1    ! indicator for Yes/No for writing
         if (iPOS(j) > 0) call fatalerr ('check_list_tz', 'Double entries in inList_csv_tz are not allowed')
         iPOS(j) = i    ! help vector to be used later for sorting
         exit
      end if
   end do
   ! error if not there; or should we write a warning?
   if (.not. isthere) call fatalerr ('check_list_tz','Illegal variable name in inList_csv_tz for write_csv output')
end do
end subroutine check_list_tz

subroutine sort_list_tz (InList_csv, ListVars, Nvars)

implicit none
! global; NB: on output both variables are arranged in same sorted order as Allowed
integer,                            intent(in)     :: Nvars
character(len=*),                   intent(out)    :: InList_csv
character(len=*), dimension(Nvars), intent(inout)  :: ListVars
! local
integer              :: il, i, j
character(len=1024)  :: temp

il = 0
! iPOS contains for each of the Allowed names the position in the user-supplied list ListVars,
! and adds this to a temporary string
do j = 1, Mlist
   if (iPOS(j) > 0) then
      if (il == 0) then
         temp = trim(listVars(iPOS(j))); il = len_trim(temp)
      else
         temp = temp(1:il) // comma // trim(listVars(iPOS(j))); il = len_trim(temp)
      end if
   end if
end do

! temp now contains the user-supplied names in the same order of appearance is in Allowed
! decompose temp into ListVars (overwrite)
call words (temp, ilw, comma, iWbeg, iWend, Nvars)
do i = 1, Nvars
   ListVars(i) = temp(iWbeg(i):iWend(i))
end do

! on return: InList_csv becomes equal to sorted string temp
InList_csv = trim(temp)

end subroutine sort_list_tz

subroutine Make_Header_tz (ListVars, Nvars, Header)
implicit none
! global
integer,                            intent(in)     :: Nvars
character(len=*), dimension(Nvars), intent(in)     :: ListVars
character(len=*),                   intent(inout)  :: Header
! local
integer :: i, il

il = len_trim(Header)
do i = 1, Nvars
   Header = Header(1:il) // "," // trim(ListVars(i)); il = len_trim(Header)
end do

end subroutine Make_Header_tz

subroutine make_headerunits_tz (HeaderUnits)
implicit none
! global
character(len=*),                   intent(inout)  :: HeaderUnits
! local
integer :: i, il

il = len_trim(HeaderUnits)
do i = 1, Mlist
   if (iCSV(i) == 1) then
      HeaderUnits = HeaderUnits(1:il) // "," // trim(Units(i)); il = len_trim(HeaderUnits)
   end if
end do

end subroutine make_headerunits_tz

end subroutine csv_write_tz
