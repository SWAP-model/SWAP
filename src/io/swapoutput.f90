! File VersionID:
!   $Id: swapoutput.f90 377 2018-04-04 10:57:55Z heine003 $
! ----------------------------------------------------------------------
      subroutine swapoutput(task, state)
      use error_mod, only: fatalerr_collected
! ----------------------------------------------------------------------
!     Date               : Aug 2004
!     Purpose            : open and write general swap output files
! ----------------------------------------------------------------------

      use Variables
      use swap_state_mod, only: swap_state_t
      implicit none

      integer task
      ! [SS-BMI2] inout: init/cleanup of water_balance_row buffer
      type(swap_state_t), intent(inout) :: state

      select case (task)
      case (1)

! === open output files ===============================
! ADR 0009 Phase 5+: outbal / outblc deleted (swbal=0, swblc=0).

         ! [SS-BMI2] allocate water balance row buffer (builder always runs)
         call init_water_balance_buffer(state)

         return

      case (2)

! ===    write actual data ===============================
! ADR 0009 Phase 5+: outbal / outblc deleted (swbal=0, swblc=0).
! [SS-BMI2] build_water_balance_row is called from outinc(2) which owns the data.

      return

      case (3)

! ===    close output files ================================
! ADR 0009 Phase 5+: bal / blc file units no longer opened.

! ---    final message log file
         write(logf,'(/,a)') ' Swap simulation okay!'
         close (logf)

         ! [SS-BMI2] deallocate water balance row buffer
         call cleanup_water_balance_buffer(state)

      case default
         call fatalerr_collected ('SWAPoutput', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine soilwateroutput(task, state)
      use error_mod, only: fatalerr_collected
! ----------------------------------------------------------------------
!     Date               : Aug 2004
!     Purpose            : open and write soil water output files
! ----------------------------------------------------------------------

      use Variables
      use SWAP_csv_output
      use SWAP_csv_output_tz
      use swap_state_mod, only: swap_state_t
      implicit none

      integer task
      ! B-2.7: inout for mini-sim writeback to state%soilwater%qbot / gwlinp
      type(swap_state_t), intent(inout) :: state

      select case (task)

      case (1)
! ===    open output files and write headers ===============================
! ADR 0009 Phase 5+: outwba/outstr/outvap/outafo/outaun/outsoilphys/
! OutputModflow/capriseoutput deleted (sw* switches forced to 0).

! --     user-defined variables in CSV file
         if (swcsv == 1) call csv_out(1, state)              ! call csv_write(1)
! --     user-defined variables in CSV file
         if (swcsv_tz == 1) call csv_out_tz(1, state)        ! call csv_write_tz(1)

! --     inc file
         if (swinc.eq.1) call outinc (1, state)

! --     rot file, only when drought stress according to De Jong van Lier
         if (swdrought.eq.2) call outrot(1, state)

! --     special output for RUME project
         if (swrum == 1) call outrume (task, state)

      case (2)
! ===    write actual data ===============================
! ADR 0009 Phase 5+: outwba/outstr/outvap/outafo/outaun/OutputModflow/
! capriseoutput deleted (sw* switches forced to 0).

! --     user-defined variables in CSV file
         if (swcsv == 1) call csv_out(2, state)              ! csv_write(2)
! --     user-defined variables in CSV file
         if (swcsv_tz == 1) call csv_out_tz(2, state)        ! csv_write_tz(2)

! --     inc file
         if (swinc.eq.1) call outinc (2, state)

! --     rot file
         ! SS-ATM A-2.5: ptra read from state%atmosphere (atmosphere home).
         if (swdrought.eq.2 .and. state%atmosphere%ptra .gt. 1.0d-10) call outrot (2, state)

! --     special output for RUME project
         if (swrum == 1) call outrume (task, state)

      case (3)
! ===    write final values end of a simulation day ===========================
! ADR 0009 Phase 5+: outend deleted (swend=0).

      case (4)
! ===    close output files ===========================
! ADR 0009 Phase 5+: wba/str/vap/afo/aun file units no longer opened;
! OutputModflow / capriseoutput close branches removed.

! --     user-defined variables in CSV file
         if (swcsv == 1) call csv_out(3, state)              ! csv_write(3)
! --     user-defined variables in CSV file
         if (swcsv_tz == 1) call csv_out_tz(3, state)        ! csv_write_tz(3)

         ! [SS-BMI2] headless guard: .inc file was only opened when not headless
         if (.not. state%timecontrol%headless) then
            if (swinc.eq.1) close (inc)
         end if

! --     special output for RUME project
         if (swrum == 1) call outrume (3, state)

      case default
         call fatalerr_collected ('SoilWaterOutput', 'Illegal value for Task')
      end select

      return
      end


! ----------------------------------------------------------------------
      subroutine outinc (task, state)
      use error_mod, only: fatalerr_collected
! ----------------------------------------------------------------------
!     date               : july 2002
!     purpose            : write water balance increments to outnam.inc file
! ---------------------------------------------------------------------
      ! SS-ATM A-2.5: igrai,isnrai,igsnow,iptra,ipeva,ievap,isubl removed from only-list; reads via state%atmosphere.
      ! SS-ATM A-2.5: ssnow,snowinco removed from only-list; reads via state%atmosphere.
      ! SS-SWC S-2.11: gwl,pond,volact,volini,PondIni,iqbot,iqrot,igird,iintc,irunon,iruno,irunoCN removed; reads via state%soilwater.
      ! SS-TC TC-7: daynr,daycum,t1900,date,flheader,flprintshort removed from only-list; reads via state%timecontrol.
      ! [SS-BMI2] inout: build_water_balance_row writes to state%water_balance_row
      use variables, only: inc,iQMpOutDrRap,outfil,pathwork,project
      use swap_state_mod, only: swap_state_t
      use file_io_mod, only: file_open
      implicit none

! --- global
      integer   task
      type(swap_state_t), intent(inout) :: state

! --- local
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

      ! SS-TC TC-7: t1900 = t1900 no-op removed; t1900 now via state%timecontrol%t1900.

! --- open output file once (headless: skip file I/O, buffer already init by SwapOutput(1))
      if (.not. state%timecontrol%headless) then
         filnam = trim(pathwork)//trim(outfil)//'.inc'
         call file_open(inc,filnam,'replace','write')
         filtext = 'water balance increments (cm/day)'
         call writehead (inc,1,filnam,filtext,project)

! --- write header of inc file
         if (state%timecontrol%flprintshort) then  ! TC-7
           write (inc,10)
         else
           write (inc,12)
         endif
      end if
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


      ! SS-ATM A-2.5: snowinco read from state%atmosphere (atmosphere home).
      ! SS-SWC S-2.11: volini,PondIni read from state%soilwater.
      VolOld  = state%soilwater%volini
      PondOld = state%soilwater%pondini
      SnowOld = state%atmosphere%snowinco

      return

      case (2)

! --- write header in case of new balance period
      ! SS-TC TC-7: flheader,flprintshort read via state%timecontrol (tc_* aliases below).
      ! SS-TC TC-7: t1900,daynr,daycum,date read via state%timecontrol (tc_* aliases below).
      ! SS-ATM A-2.5: ssnow,snowinco,igrai,isnrai,igsnow,iptra,ipeva,ievap,isubl from state%atmosphere.
      ! SS-SWC S-2.11: gwl,pond,volact,PondIni,igird,iintc,irunon,iruno,irunoCN,iqrot,iqbot via state%soilwater.
      associate( &
        tc_flheader     => state%timecontrol%flheader,      &  ! TC-7
        tc_flprintshort => state%timecontrol%flprintshort,  &  ! TC-7
        tc_t1900        => state%timecontrol%t1900,         &  ! TC-7
        tc_daynr        => state%timecontrol%daynr,         &  ! TC-7
        tc_daycum       => state%timecontrol%daycum,        &  ! TC-7
        tc_date         => state%timecontrol%date,          &  ! TC-7
        at_igrai  => state%atmosphere%intr%igrai,  &
        at_inrai  => state%atmosphere%intr%inrai,  &
        at_igsnow => state%atmosphere%intr%igsnow, &
        at_iptra  => state%atmosphere%intr%iptra,  &
        at_ipeva  => state%atmosphere%intr%ipeva,  &
        at_ievap  => state%atmosphere%intr%ievap,  &
        at_isubl  => state%atmosphere%intr%isubl,  &
        at_isnrai => state%atmosphere%intr%isnrai, &
        at_ssnow  => state%atmosphere%ssnow,       &
        sw_gwl    => state%soilwater%gwl,          &
        sw_pond   => state%soilwater%pond,         &
        sw_volact => state%soilwater%volact,       &
        sw_igird  => state%soilwater%igird,  &
        sw_iintc  => state%soilwater%iintc,  &
        sw_irunon => state%soilwater%irunon, &
        sw_iruno  => state%soilwater%iruno,  &
        sw_irunoCN => state%soilwater%irunoCN, &
        sw_iqrot  => state%soilwater%iqrot,  &
        sw_iqbot  => state%soilwater%iqbot   &
      )

! --- compute derived terms (always, for state buffer)
      dstor = (sw_volact + sw_pond + at_ssnow) - (VolOld + PondOld + SnowOld)
      baldev = (at_igrai+at_isnrai+at_igsnow+sw_igird+sw_irunon) - dstor -    &
     & (sw_iintc+sw_iruno+sw_irunoCN+sw_iqrot+at_ievap+at_isubl+iQMpOutDrRap+state%surfacewater%iqdra+(-1.0d0*sw_iqbot))

! --- [SS-BMI2] build water balance row into state buffer (always runs, headless-independent)
      call build_water_balance_row(state, dstor, baldev)

! --- write to .inc file (headless: skip)
      if (.not. state%timecontrol%headless) then
         if (tc_flheader) then
           if (tc_flprintshort) then
             write (inc,10)
           else
             write (inc,12)
           endif
         endif

! --- determine date and date-time
         call dtdpst ('year-month-day,hour:minute:seconds',tc_t1900,datexti)  ! TC-7

! --- write output record .inc file
         gwlout = "          "
         if (sw_gwl.lt.998.0d0)  write(gwlout,'(f9.1)') sw_gwl
         if (tc_flprintshort) then
           write (inc,20) datexti,comma,tc_daynr,comma,tc_daycum,comma,     &
     &       at_igrai+at_isnrai,comma,at_igsnow,comma,sw_igird,comma,sw_iintc, &
     &       comma,sw_irunon,comma,sw_iruno+sw_irunoCN,comma,at_iptra,comma,sw_iqrot, &
     &       comma,at_ipeva,comma,at_ievap,comma,                           &
     &       (iQMpOutDrRap+state%surfacewater%iqdra),comma,sw_iqbot,&
     &       comma,gwlout,comma,dstor,comma,baldev            !comma,storage
         else
           write (inc,22) tc_date,comma,tc_daynr,comma,tc_daycum,comma,     &
     &       at_igrai+at_isnrai,comma,at_igsnow,comma,sw_igird,comma,sw_iintc, &
     &       comma,sw_irunon,comma,sw_iruno+sw_irunoCN,comma,at_iptra,comma,sw_iqrot, &
     &       comma,at_ipeva,comma,at_ievap,comma,                           &
     &       (iQMpOutDrRap+state%surfacewater%iqdra),comma,sw_iqbot,&
     &       comma,gwlout,comma,dstor,comma,baldev           !comma,storage
         endif
      end if
 20   format (a19,a1,i3,a1,i6,12(a1,f10.5),2a,2(a1,f10.5))     !,(a1,e12.5)
 22   format (a11,a1,i3,a1,i6,12(a1,f10.5),2a,2(a1,f10.5))     !,(a1,e12.5)

      VolOld = sw_volact
      PondOld = sw_pond
      SnowOld = at_ssnow
      end associate

      case default
         call fatalerr_collected ('outinc', 'Illegal value for Task')
      end select

      return
      end


! ----------------------------------------------------------------------
! [SS-BMI2] Water balance buffer helpers (canonical output-sink pattern)
! These four subroutines are the canonical template for Tasks 9-15.
! ----------------------------------------------------------------------

      subroutine init_water_balance_buffer(state)
! ----------------------------------------------------------------------
!     Allocate state%water_balance_row and set column names.
!     Called from SwapOutput(1) — always runs, headless-independent.
!     N_COLS = 18: t1900(date), daynr, daycum, rain, snow, irrig,
!                  interc, runon, runoff, tpot, tact, epot, eact,
!                  drainage, qbottom, gwl, dstorage, baldev
! ----------------------------------------------------------------------
      use swap_state_mod, only: swap_state_t
      use iso_c_binding,  only: c_double
      implicit none
      type(swap_state_t), intent(inout) :: state
      integer, parameter :: N = 18

      state%water_balance_n_cols = N
      if (.not. allocated(state%water_balance_row))     allocate(state%water_balance_row(N))
      if (.not. allocated(state%water_balance_columns)) allocate(state%water_balance_columns(N))
      state%water_balance_row     = 0.0_c_double
      state%water_balance_columns(1)  = 'date'
      state%water_balance_columns(2)  = 'day'
      state%water_balance_columns(3)  = 'dcum'
      state%water_balance_columns(4)  = 'rain'
      state%water_balance_columns(5)  = 'snow'
      state%water_balance_columns(6)  = 'irrig'
      state%water_balance_columns(7)  = 'interc'
      state%water_balance_columns(8)  = 'runon'
      state%water_balance_columns(9)  = 'runoff'
      state%water_balance_columns(10) = 'tpot'
      state%water_balance_columns(11) = 'tact'
      state%water_balance_columns(12) = 'epot'
      state%water_balance_columns(13) = 'eact'
      state%water_balance_columns(14) = 'drainage'
      state%water_balance_columns(15) = 'qbottom'
      state%water_balance_columns(16) = 'gwl'
      state%water_balance_columns(17) = 'dstorage'
      state%water_balance_columns(18) = 'baldev'
      end subroutine init_water_balance_buffer


      subroutine build_water_balance_row(state, dstor, baldev)
! ----------------------------------------------------------------------
!     Fill state%water_balance_row(:) with the daily water balance values.
!     Called from outinc(2) — always runs, headless-independent.
!     Column order matches init_water_balance_buffer column names exactly.
! ----------------------------------------------------------------------
      use swap_state_mod, only: swap_state_t
      use variables,      only: iQMpOutDrRap
      use iso_c_binding,  only: c_double
      implicit none
      type(swap_state_t), intent(inout) :: state
      real(8), intent(in) :: dstor, baldev

      state%water_balance_row(1)  = real(state%timecontrol%t1900,       c_double)  ! date (days since 1900-01-01)
      state%water_balance_row(2)  = real(state%timecontrol%daynr,       c_double)  ! day of year
      state%water_balance_row(3)  = real(state%timecontrol%daycum,      c_double)  ! cumulative day
      state%water_balance_row(4)  = real(state%atmosphere%intr%igrai  + &
                                         state%atmosphere%intr%isnrai, c_double)  ! rain (gross + snow rain)
      state%water_balance_row(5)  = real(state%atmosphere%intr%igsnow,  c_double)  ! snow
      state%water_balance_row(6)  = real(state%soilwater%igird,         c_double)  ! irrig
      state%water_balance_row(7)  = real(state%soilwater%iintc,         c_double)  ! interc
      state%water_balance_row(8)  = real(state%soilwater%irunon,        c_double)  ! runon
      state%water_balance_row(9)  = real(state%soilwater%iruno        + &
                                         state%soilwater%irunoCN,      c_double)  ! runoff
      state%water_balance_row(10) = real(state%atmosphere%intr%iptra,   c_double)  ! tpot
      state%water_balance_row(11) = real(state%soilwater%iqrot,         c_double)  ! tact
      state%water_balance_row(12) = real(state%atmosphere%intr%ipeva,   c_double)  ! epot
      state%water_balance_row(13) = real(state%atmosphere%intr%ievap,   c_double)  ! eact
      state%water_balance_row(14) = real(iQMpOutDrRap                 + &
                                         state%surfacewater%iqdra,     c_double)  ! drainage
      state%water_balance_row(15) = real(state%soilwater%iqbot,         c_double)  ! qbottom
      state%water_balance_row(16) = real(state%soilwater%gwl,           c_double)  ! gwl (999.0 when not simulated)
      state%water_balance_row(17) = real(dstor,                          c_double)  ! dstorage
      state%water_balance_row(18) = real(baldev,                         c_double)  ! baldev
      end subroutine build_water_balance_row


      subroutine cleanup_water_balance_buffer(state)
! ----------------------------------------------------------------------
!     Deallocate state%water_balance_row and reset counter.
!     Called from SwapOutput(3) — always runs, headless-independent.
! ----------------------------------------------------------------------
      use swap_state_mod, only: swap_state_t
      implicit none
      type(swap_state_t), intent(inout) :: state

      if (allocated(state%water_balance_row))     deallocate(state%water_balance_row)
      if (allocated(state%water_balance_columns)) deallocate(state%water_balance_columns)
      state%water_balance_n_cols = 0
      end subroutine cleanup_water_balance_buffer


! ----------------------------------------------------------------------
      subroutine outrot (task, state)
      use error_mod, only: fatalerr_collected

! ----------------------------------------------------------------------
!     date               : February 2012
!     purpose            : write output of microscopic root water uptake
! ---------------------------------------------------------------------
      ! SS-CRP Phase 2 Task C-2.4: qrot, hroot, mroot, mflux, rootrho, rootphi, hleaf, hxylem
      !   removed from only-list; reads via state%soilwater.
      ! SS-SWC S-2.11: theta,hm1,q,inq,inqrot removed from only-list; reads via state%soilwater.
      ! SS-TC TC-7: daynr,daycum,t1900,date,flprintshort removed from only-list; reads via state%timecontrol.
      use variables, only: ztopcp, zbotcp,rot,noddrz,z,outfil,pathwork,project
      use swap_state_mod, only: swap_state_t
      use file_io_mod, only: file_open
      implicit none

! --- global
      integer   task
      type(swap_state_t), intent(in) :: state

! --- local
      integer   node
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
      call file_open(rot,filnam,'replace','write')
      filtext = 'microscopic root water uptake'
      call writehead (rot,1,filnam,filtext,project)
      write (rot,100)

! --- write header in rot file
      ! SS-TC TC-7: flprintshort,daynr,daycum,t1900,date read via state%timecontrol (tc_* below).
      if (state%timecontrol%swheader .eq. 0) then  ! [SS-BMI2 Task 4]
        if (state%timecontrol%flprintshort) then  ! TC-7
          write (rot,200)
        else
          write (rot,210)
        endif
      endif

      return

      case (2)

! === write actual profile data ===========================================

! --- write header in rot file and write actual data
      ! SS-TC TC-7: flprintshort,t1900,daynr,daycum,date -> state%timecontrol tc_* aliases.
      ! SS-SWC S-2.11: theta,hm1,q,inq,inqrot read from state%soilwater.
      associate( &
        tc_flprintshort => state%timecontrol%flprintshort,  &  ! TC-7
        tc_t1900        => state%timecontrol%t1900,         &  ! TC-7
        tc_daynr        => state%timecontrol%daynr,         &  ! TC-7
        tc_daycum       => state%timecontrol%daycum,        &  ! TC-7
        tc_date         => state%timecontrol%date           &  ! TC-7
      )
      if (state%timecontrol%swheader .eq. 1) then  ! [SS-BMI2 Task 4]
        if (tc_flprintshort) then
          write (rot,200)
        else
          write (rot,210)
        endif
      endif

      if (tc_flprintshort) then
! ---   determine date and date-time
        call dtdpst ('year-month-day,hour:minute:seconds',tc_t1900,datexti)  ! TC-7
        do node = 1,noddrz
           write (rot,300) datexti,comma,z(node),comma,state%soilwater%hleaf,comma,  &
     &       state%soilwater%Hxylem,comma,                                            &
     &       state%soilwater%hroot(node),comma,state%soilwater%hm1(node),comma,state%soilwater%inqrot(node),comma, &
     &       state%soilwater%qrot(node),comma,state%soilwater%inq(node),comma,state%soilwater%q(node), &
     &       comma,state%soilwater%mroot(node),                                       &
     &       comma,state%soilwater%mflux(node),comma,state%soilwater%rootrho(node),  &
     &       comma,state%soilwater%rootphi(node),                                     &
     &       comma,state%soilwater%theta(node),comma,ztopcp(node),comma,              &
     &       zbotcp(node),comma,tc_daynr,comma,tc_daycum
        end do

      else
        do node = 1,noddrz
           write (rot,310) tc_date,comma,z(node),comma,state%soilwater%hleaf,comma,  &
     &       state%soilwater%Hxylem,comma,                                            &
     &       state%soilwater%hroot(node),comma,state%soilwater%hm1(node),comma,state%soilwater%inqrot(node),comma, &
     &       state%soilwater%qrot(node),comma,state%soilwater%inq(node),comma,state%soilwater%q(node), &
     &       comma,state%soilwater%mroot(node),                                       &
     &       comma,state%soilwater%mflux(node),comma,state%soilwater%rootrho(node),  &
     &       comma,state%soilwater%rootphi(node),                                     &
     &       comma,state%soilwater%theta(node),comma,ztopcp(node),comma,              &
     &       zbotcp(node),comma,tc_daynr,comma,tc_daycum
        end do

      endif
      end associate  ! tc_flprintshort, tc_t1900, tc_daynr, tc_daycum, tc_date (TC-7)

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
         call fatalerr_collected ('outrot', 'Illegal value for Task')
      end select

      return
      end


! ----------------------------------------------------------------------
      subroutine OutCropFixed(task, state)
      use error_mod, only: fatalerr_collected
! ----------------------------------------------------------------------
!     Date               : Aug 2004
!     Purpose            : open and write fixed crop output files
! SS-TC TC-10: state added (intent in); date,t read via state%timecontrol
!   tc_* aliases. daycrop remains crop-domain variable (not in TC).
! ----------------------------------------------------------------------
      use variables, only: daycrop,dvs,tsum,lai,cf,rd,crp,ch
      use swap_state_mod, only: swap_state_t
      implicit none

! --- global variables ------------------
      integer task
      type(swap_state_t), intent(in) :: state

! --- local
      character(len=1) comma
! --- TC-10: date, t read via state%timecontrol tc_* aliases.
      associate( &
        tc_date => state%timecontrol%date,  &  ! TC-10
        tc_t    => state%timecontrol%t      &  ! TC-10
      )
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
      write (crp,200) tc_date,comma,nint(tc_t),comma,daycrop,comma,dvs,comma, &
     & tsum,comma,"       ",comma,lai,comma,ch,comma,cf,                &
     & comma,"       ",comma,nint(rd),                                  &
     & comma,comma,comma,comma,comma,comma,comma,comma,comma,comma,     &
     & comma,comma,comma,comma,comma,comma,                             &
     & comma,comma,comma,comma,comma,comma
 200  format (a11,a1,i5,a1,i7,a1,f6.2,a1,f6.0,a1,a,2(a1,f7.2),(a1,f6.2),&
     &      (a1,a),(a1,i7),16(a1,'        ') ,                          &
     &      6(a1,'         ') )

      case default
         call fatalerr_collected ('OutCropFixed', 'Illegal value for Task')
      end select

      end associate  ! tc_date, tc_t => state%timecontrol [TC-10]
      return
      end

! ----------------------------------------------------------------------
      subroutine OutWofost(task, state)
      use error_mod, only: fatalerr_collected
! ----------------------------------------------------------------------
!     UpDate             : Aug 2014
!     Date               : Oct 2004
!     Purpose            : Write detailed crop growth output files
! SS-TC TC-10: state added (intent in); date,t read via state%timecontrol
!   tc_* aliases. daycrop remains crop-domain variable (not in TC).
! ----------------------------------------------------------------------
      use variables, only: daycrop,crp,dvs,tsum,laipot,lai,rdpot,rd,ch,cf,cwdmpot,cwdm,wsopot,wso,wlvpot,wlv,wstpot,   &
                           wst,wrtpot,wrt,dwlvCrop,dwlvSoil,dwst,dwrt,dwso,HarLosOrm_tot,swbulb,wblpot,wbl,dwblpot,dwbl,plwt
      use swap_state_mod, only: swap_state_t
      implicit none

! --- global variables ------------------
      integer task
      type(swap_state_t), intent(in) :: state

! --- local variables ------------------
      character(len=1) comma
! --- TC-10: date, t read via state%timecontrol tc_* aliases.
      associate( &
        tc_date => state%timecontrol%date,  &  ! TC-10
        tc_t    => state%timecontrol%t      &  ! TC-10
      )
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
         write (crp,300) tc_date,comma,nint(tc_t),comma,daycrop,comma,dvs,    &
     &    comma,tsum,comma,laipot,comma,lai,comma,ch,comma,cf,comma,    &
     &    nint(rdpot),comma,nint(rd),comma,nint(wlvpot),comma,nint(wlv),&
     &    comma,nint(wstpot),comma,nint(wst),comma,nint(wrtpot),        &
     &    comma,nint(wrt),comma,nint(cwdmpot),comma,nint(cwdm),comma,   &
     &    nint(wsopot),comma,nint(wso),comma,comma,comma,comma,comma,   &
     &    comma,comma,dwlvCrop,comma,dwlvSoil,comma,dwst,comma,dwrt,    &
     &    comma,dwso,comma,HarLosOrm_tot
      elseif(swbulb.eq.1) then
         write (crp,400) tc_date,comma,nint(tc_t),comma,daycrop,comma,dvs,    &
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
         call fatalerr_collected ('OutWofost', 'Illegal value for Task')
      end select

      end associate  ! tc_date, tc_t => state%timecontrol [TC-10]
      return
      end

! ----------------------------------------------------------------------
      subroutine OutGrass(task, state)
      use error_mod, only: fatalerr_collected
! ----------------------------------------------------------------------
!     UpDate             : Aug 2014
!     Date               : Oct 2004
!     Purpose            : Write detailed grass simulation output files
! SS-TC TC-10: state added (intent in); date,t read via state%timecontrol
!   tc_* aliases. daycrop remains crop-domain variable (not in TC).
! ----------------------------------------------------------------------
      use variables, only: daycrop,crp,dvs,tsum,laipot,lai,rdpot,rd,ch,cf,tagppot,tagp,tagptpot,tagpt,          &
                           wlvpot,wlv,wstpot,wst,wrtpot,wrt,cuptgraz,cuptgrazpot
      use swap_state_mod, only: swap_state_t
      implicit none

! --- global variables ------------------
      integer task
      type(swap_state_t), intent(in) :: state

! --- local variables ------------------
      character(len=1) comma
! --- TC-10: date, t read via state%timecontrol tc_* aliases.
      associate( &
        tc_date => state%timecontrol%date,  &  ! TC-10
        tc_t    => state%timecontrol%t      &  ! TC-10
      )
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
      write (crp,200) tc_date,comma,nint(tc_t),comma,daycrop,comma,dvs,comma, &
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
         call fatalerr_collected ('OutGrass', 'Illegal value for Task')
      end select

      end associate  ! tc_date, tc_t => state%timecontrol [TC-10]
      return
      end

! ----------------------------------------------------------------------
      subroutine SoluteOutput(task, state)
      use error_mod, only: fatalerr_collected
! ----------------------------------------------------------------------
!     Date               : November 2004
!     Purpose            : open and write solute output files
! ----------------------------------------------------------------------

      use Variables
      use swap_state_mod, only: swap_state_t
      implicit none

      integer task
      type(swap_state_t), intent(in) :: state

      select case (task)
      case (1)

! === open output files and write headers ===============================
! ADR 0009 Phase 5+: outsba deleted (swsba=0).

      return

      case (2)

! === write actual data ===============================
! ADR 0009 Phase 5+: outsba deleted (swsba=0).

      return

      case (3)

! === close output files ===========================
! ADR 0009 Phase 5+: sba file unit no longer opened.

      case default
         call fatalerr_collected ('SoluteOutput', 'Illegal value for Task')
      end select

      return
      end


! ----------------------------------------------------------------------
      subroutine AgeTracerOutput(task, state)
      use error_mod, only: fatalerr_collected
! ----------------------------------------------------------------------
!     Date               : October 2010
!     Purpose            : open and write Groundwater Ageing output files
! ----------------------------------------------------------------------
      use swap_state_mod, only: swap_state_t
      implicit none

! --- global variables ------------------
      integer task
      type(swap_state_t), intent(in) :: state
! --- local variables ------------------
      integer agep,agee,ageq

      save   agee,agep,ageq

      select case (task)
      case (1)

! === open output files and write headers ===============================

! --  age files
      call outage (1,agep,agee,ageq,state)

      return

      case (2)

! === write actual data ===============================

! --  age files
      call outage (2,agep,agee,ageq,state)

      return

      case (3)

! === close output files ===========================

! --- close sba file
      close (agep)
      close (agee)

      case default
         call fatalerr_collected ('AgeTracerOutput', 'Illegal value for Task')
      end select

      return
      end
! ----------------------------------------------------------------------
      subroutine outage(task,agep,agee,ageq,state)
      use error_mod, only: fatalerr_collected
! ----------------------------------------------------------------------
!     date               : October 2010
!     purpose            : output of groundwater age
! ---------------------------------------------------------------------
! --- ADR 0032: AgeTracer is currently inert (flAgeTracer never .true.).
!     This routine is already unreachable via call-site guards in swap.f90,
!     but the explicit guard here makes the gating visible at the definition.
      ! SS-SWST Phase 2 Task 11: inqdra removed (now via state%surfacewater%inqdra).
      ! SS-SLST Phase 1 Task 5: cml migrated to state%solute (body is gated by flAgeTracer guard).
      ! SS-TC TC-7: daynr,daycum,date,outper removed from only-list; reads via state%timecontrol.
      use variables, only: project,nrlevs,outfil,pathwork,numnod,z,                                    &
                           AgeGwl1m,icAgeBot,icAgeDra,icAgeRot,icAgeSur,flAgeTracer
      use swap_state_mod, only: swap_state_t
      use swap_array_dimensions, only: madr
      use file_io_mod, only: file_open
                           implicit none

! --- global
      integer   agep,agee,ageq,task
      type(swap_state_t), intent(in) :: state

! --- local variables ------------------
      integer   reclngth,node,level
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
      real(8)   iqdrainout(madr)   ! Cumulative (over 1 output timestep) drainage flux (L) for each drainage level
! ----------------------------------------------------------------------
      if (.not. flAgeTracer) return   ! AgeTracer is currently inert (ADR 0032); output gated

      comma = ','

      select case (task)
      case (1)

! --- open output files -------------------------------------------------
!     age of groundwater as profile
      filnam = trim(pathwork)//trim(outfil)//'.ageProfile.csv'
      reclngth = 50 + 12*numnod
      open(newunit=agep,file=filnam,status='unknown',recl=reclngth)
      filtext = 'Groundwater age profiles (all age-values in days)'
      call writehead (agep,1,filnam,filtext,project)
!     age of groundwater in effluents: drains, transpiration, leaching, runoff
      filnam = trim(pathwork)//trim(outfil)//'.ageEffluent.csv'
      call file_open(agee,filnam,'replace','write')
      filtext = 'Groundwater age effluent (all age-values in days)'
      call writehead (agee,1,filnam,filtext,project)
!     effluent drain water fluxes
      filnam = trim(pathwork)//trim(outfil)//'.ageEffluentqDrain.csv'
      call file_open(ageq,filnam,'replace','write')
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

      ! SS-TC TC-7: date,daynr,daycum,outper read via state%timecontrol tc_* aliases.
      associate( &
        tc_date   => state%timecontrol%date,    &  ! TC-7
        tc_daynr  => state%timecontrol%daynr,   &  ! TC-7
        tc_daycum => state%timecontrol%daycum,  &  ! TC-7
        tc_outper => state%timecontrol%outper   &  ! TC-7
      )

!     age of groundwater as profile
      write(agep,15) tc_date,comma,tc_daynr,comma,tc_daycum,             &
     &               (comma,state%solute%cml(node),node=1,numnod)
 15   format(a11,a1,i4,a1,i6,1p,1024(a1,e10.3))

!     age of groundwater in effluents: drains, transpiration, leaching, runoff
!     and age (d) of groundwater in upper 1 meter of saturated zone
      write(agee,16) tc_date,comma,tc_daynr,comma,tc_daycum,comma,AgeGwl1m,comma, &
     &               icAgeBot/tc_outper,comma,icAgeRot/tc_outper,comma,  &
     &               icAgeSur/tc_outper,                                  &
     &               (comma,icAgeDra(level)/tc_outper,level=1,nrlevs)
 16   format(a11,a1,i4,a1,i6,1p,9(a1,e10.3))

!     qdrain discharge-effluent (without infiltration!)
      ! SS-SWST Phase 2 Task 11: inqdra global fallback removed; state is authoritative.
      do level = 1,nrlevs
        iqdrainout(level) = 0.0d0
        if (allocated(state%surfacewater%inqdra)) then
          do node = 1,numnod
            if (state%surfacewater%inqdra(level,node).gt.0.0d0) then
             iqdrainout(level) = iqdrainout(level) + state%surfacewater%inqdra(level,node)
            endif
          enddo
        end if
      enddo
      write(ageq,16) tc_date,comma,tc_daynr,comma,tc_daycum,             &
     &               (comma,iqdrainout(level),level=1,nrlevs)

      end associate  ! tc_date, tc_daynr, tc_daycum, tc_outper (TC-7)

      case default
         call fatalerr_collected ('outage', 'Illegal value for Task')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine TemperatureOutput(task, state)
      use error_mod, only: fatalerr_collected
! ----------------------------------------------------------------------
!     Date               : November 2004
!     Purpose            : open and write soil temperature output files
! ----------------------------------------------------------------------

! --- global variables ------------------
      use Variables
      use swap_state_mod, only: swap_state_t

      implicit none
      integer task
      type(swap_state_t), intent(in) :: state


      select case (task)
      case (1)

! === open output files and write headers ===============================

! ADR 0009 Phase 5+: outheapar deleted (swini=0).

! --  tem file
      if (swtem .eq. 1) call outtem (task, state)

      return

      case (2)

! === write actual data ===============================

! --  tem file
      if (swtem .eq. 1) call outtem (task, state)

      return

      case (3)

! === close output files ===========================

! --- close tem file
      if (swtem .eq. 1) close (tem)

      case default
         call fatalerr_collected ('TemperatureOutput', 'Illegal value for Task')
      end select

      return
      end





! ----------------------------------------------------------------------
      subroutine outtem (task, state)
      use error_mod, only: fatalerr_collected

! ----------------------------------------------------------------------
!     date               : November 2004
!     purpose            : Output of soil temperatures
! ---------------------------------------------------------------------
      ! SS-HEAT Phase 1 Task 5: tsoil, tebot, tetop migrated to state%heat.
      ! SS-TC TC-7: date,daynr,daycum,flheader removed from only-list; reads via state%timecontrol.
      use variables, only: numnod,tem,tav,outfil,pathwork,project
      use swap_state_mod, only: swap_state_t
      implicit none

! --- global
      integer   task
      type(swap_state_t), intent(in) :: state

! --- local variables ------------------
      integer      i, reclngth
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)

! === open output file =================================================
      filnam = trim(pathwork)//trim(outfil)//'.tem'
!      reclngth = 36 + 7*numnod
      reclngth = 50 + 7*numnod
      open(newunit=tem,file=filnam,status='unknown',recl=reclngth)
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

      ! SS-TC TC-7: daynr,daycum -> state%timecontrol (direct, case(1) only line).
      write (tem,'(a11,a1,i3,a1,i6,1024(a1,f6.1:))') '    Initial'      &
     &      ,comma,state%timecontrol%daynr,comma,state%timecontrol%daycum, &
     &      comma,tav,comma,state%heat%tetop,                              &
     &      (comma,state%heat%tsoil(i),i=1,numnod),comma,state%heat%tebot

      return

      case (2)

! === write actual soil temperature data ================================

      ! SS-TC TC-7: flheader,date,daynr,daycum -> state%timecontrol tc_* aliases.
      associate( &
        tc_flheader => state%timecontrol%flheader,  &  ! TC-7
        tc_date     => state%timecontrol%date,      &  ! TC-7
        tc_daynr    => state%timecontrol%daynr,     &  ! TC-7
        tc_daycum   => state%timecontrol%daycum     &  ! TC-7
      )
! --- write header in case of new balance period
      if (tc_flheader) write (tem,10)

! --- write soil temperature profile
!     PWB: idem
!      write (tem,'(a11,a1,i3,a1,i6,<numnod+3>(a1,f6.1:))') date
!     &      ,comma,daynr,comma,daycum,comma,tav,comma,state%heat%tetop,
!     &      (comma,state%heat%tsoil(i),i=1,numnod),comma,state%heat%tebot
      write (tem,'(a11,a1,i3,a1,i6,1024(a1,f6.1:))') tc_date             &
     &      ,comma,tc_daynr,comma,tc_daycum,comma,tav,comma,state%heat%tetop, &
     &      (comma,state%heat%tsoil(i),i=1,numnod),comma,state%heat%tebot
      end associate  ! tc_flheader, tc_date, tc_daynr, tc_daycum (TC-7)

      case default
         call fatalerr_collected ('outtem', 'Illegal value for Task')
      end select

      return
      end


! ----------------------------------------------------------------------
      subroutine SnowOutput(task, state)
      use error_mod, only: fatalerr_collected
! ----------------------------------------------------------------------
!     Date               : December 2004
!     Purpose            : open and write snow pack data
! ----------------------------------------------------------------------
      ! SS-ATM A-2.5: snrai,gsnow,ssnow,melt,subl reads migrated to state%atmosphere.
      ! SS-TC TC-7: date,daycum,flheader removed; reads via state%timecontrol.

      use variables, only: pathwork,outfil,project,snw
      use swap_state_mod, only: swap_state_t
      use file_io_mod, only: file_open
      implicit none

! --- global variables ------------------
      integer task
      type(swap_state_t), intent(in) :: state
! --- local variables ------------------
      character(len=300) filnam
      character(len=80)  filtext
      character(len=1)   comma
! ----------------------------------------------------------------------
      comma = ','

      select case (task)
      case (1)

! === open output file =================================================
      filnam = trim(pathwork)//trim(outfil)//'.snw'
      call file_open(snw,filnam,'replace','write')
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
      ! SS-TC TC-7: flheader,date,daycum -> state%timecontrol (direct refs, single-use).
      if (state%timecontrol%flheader) write (snw,10)  ! TC-7

! --- write actual data
      ! SS-ATM A-2.5: snrai,gsnow,ssnow,melt,subl read from state%atmosphere (atmosphere home).
      write (snw,20) state%timecontrol%date,comma,state%timecontrol%daycum, &  ! TC-7
     &               comma,state%atmosphere%snrai,comma,state%atmosphere%gsnow, &
     &               comma,state%atmosphere%ssnow,comma,state%atmosphere%melt,comma,state%atmosphere%subl
20    format (a11,a1,i6,1x,5(a1,f10.4))

      return

      case (3)

! === close output file ===========================

! --- close snw file
      close (snw)

      case default
         call fatalerr_collected ('SnowOutput', 'Illegal value for Task')
      end select

      return
      end




! ----------------------------------------------------------------------
      subroutine WriteSwapOk(Project)
!-----------------------------------------------------------------------
!      Date    : December 2004
!      Purpose : Create file Swap.ok at end of simulation
!-----------------------------------------------------------------------
      use file_io_mod, only: file_open
      implicit none

      integer   cexf
      character(len=80) project

!     create file Swap.ok to let environment programs verify termination
      call file_open(cexf,'Swap.ok','replace','write')
      call writehead (cexf,1,'Swap.ok',                                 &
     &  'this header only: simulation succesfully terminated',project)
      close (cexf)

      return
      end


! ----------------------------------------------------------------------
      subroutine SurfaceWaterOutput(task, state)
      use error_mod, only: fatalerr_collected
! ----------------------------------------------------------------------
!     Date               : Aug 2004
!     Purpose            : open and write surface water output files
! ----------------------------------------------------------------------

      use variables
      use swap_state_mod, only: swap_state_t
      implicit none

      integer task
      type(swap_state_t), intent(in) :: state

      select case (task)
      case (1)

! === open output files and write headers ===============================
! ADR 0009 Phase 5+: outdrf / outswb deleted (swdrf=0, swswb=0).

      return

      case (2)

! === write actual data ===============================
! ADR 0009 Phase 5+: outdrf / outswb deleted (swdrf=0, swswb=0).

      return

      case (3)

! === close output files ===========================
! ADR 0009 Phase 5+: drf / swb file units no longer opened.

      case default
         call fatalerr_collected ('SurfaceWaterOutput', 'Illegal value for Task')
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
!     purpose            : delete temporary files
! ----------------------------------------------------------------------
      implicit  none

      logical fileopen

! --- delete temporary files (TTutil scratch retired with ADR 0023).
! The unit-20 close-with-DELETE is the project's own cleanup, kept.
      inquire(unit=20,opened=fileopen)
      if(fileopen)close(20, STATUS = 'DELETE')

      return
      end


! ----------------------------------------------------------------------
      subroutine checkDiscrVert()
      use error_mod, only: fatalerr_collected
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
      use swap_array_dimensions, only: macp
      
      implicit none

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
          call fatalerr_collected ('checkDiscrVert',messag)
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
        call fatalerr_collected ('checkDiscrVert',messag)
      endif


      return
      end



subroutine outrume (task, state)
use error_mod, only: fatalerr_collected
! SS-ATM A-2.5: igrai removed from only-list; reads via state%atmosphere%intr.
! SS-SWC S-2.11: theta,thetas,iruno,pond,gwl removed from only-list; reads via state%soilwater.
! SS-TC TC-7: tcum,outper removed from only-list; reads via state%timecontrol.
use variables, only: numnod,dz
use swap_state_mod, only: swap_state_t
implicit none

! global
integer, intent(in) :: task
type(swap_state_t), intent(in) :: state
! local
integer                       :: i
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
!!!   open (unit=iunout, file=fout, status='unknown', form='unformatted')
   open (newunit=iunout2, file=fout2, status='unknown', form='formatted')
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
   ! SS-TC TC-7: tcum,outper -> state%timecontrol via tc_* aliases.
   associate( &
      tc_tcum  => state%timecontrol%tcum,   &  ! TC-7
      tc_outper => state%timecontrol%outper  &  ! TC-7
   )
   VT = 0.0d0
   do i = 1, LZnod
      VT = VT + dz(i)*(state%soilwater%thetas(i)-state%soilwater%theta(i))
   end do
   WC = 0.0d0
   do i = 1, LZnod2
      WC = WC + dz(i)*state%soilwater%theta(i)
   end do
   WC = WC/LZ2

!  d, cm, cm/d, cm/d, cm, cm
!!!   write (iunout) real(tc_tcum),real(VT),real(igrai/tc_outper),real(iruno/tc_outper),real(pond),real(gwl)
   ! SS-ATM A-2.5: igrai read from state%atmosphere%intr (atmosphere home).
   ! SS-SWC S-2.11: iruno read from state%soilwater.
   RainRate   = state%atmosphere%intr%igrai/tc_outper
   RunoffRate = state%soilwater%iruno/tc_outper
   iDay       = int(tc_tcum)
!!!   write (iunout,'(20F12.6)') real(tc_tcum),real(VT)

   if (Nrec > 0) then
      if (.not. Event) then
         if (iDay > iDayOld) then  ! new day, no runoff occurring; restart summations
            SumRain   = RainRate*(tc_tcum-iDay)
            SumRunOff = RunoffRate*(tc_tcum-iDay) ! likely zero
            VTstart   = VT
            WCstart   = WC
            hstart    = state%soilwater%pond
            GWLstart  = state%soilwater%gwl
            iDayOld   = iDay
         else if (RunoffRate > 0.0d0) then  ! start new runoff event!
            Event  = .true.
            Nevent = Nevent + 1
            Tr     = tc_tcum-int(tc_tcum)
            Ia     = SumRain*10.0d0 ! mm
            SS     = VTstart*10.0d0 ! mm
            hh     = hstart*10.0d0  ! mm
            T      = tc_tcum
            GG     = GWLstart
            DayCum = int(tc_tcum) + 1
         end if
         SumRain   = SumRain   + RainRate*(tc_tcum-DayOld)
         SumRunoff = SumRunoff + RunoffRate*(tc_tcum-DayOld)
      else
         if (RunoffRate < 1.0D-10) then  ! stop runoff event
            Event = .false.
            P     = SumRain*10.0 ! mm
            Q     = SumRunoff*10.0 ! mm
            CN    = 25400.0d0/((P-Ia)**2/Q - (P-Ia) + 254.0d0)
            S     = 25400.0d0/CN - 254.0d0
            write (iunout2,'(I12,11(A,F12.5))') DayCum, ',', T, ',', Tr, ',', P, ',', Q, ',', Ia, ',', SS, ',', CN, ',', S, ',', hh, ',', GG, ',', WCstart
         end if
         SumRain   = SumRain   + RainRate*(tc_tcum-DayOld)
         SumRunoff = SumRunoff + RunoffRate*(tc_tcum-DayOld)
      end if
      DayOld    = tc_tcum

   else
      Nrec      = Nrec + 1
      SumRain   = RainRate*tc_tcum
      SumRunOff = RunoffRate*tc_tcum
      VTstart   = VT
      GWLstart  = state%soilwater%gwl
      hstart    = state%soilwater%pond
      iDayOld   = iDay
      DayOld    = tc_tcum
   end if
   end associate  ! tc_tcum, tc_outper (TC-7)

case (3)
!!!   close (iunout)
   close (iunout2)
case default
   call fatalerr_collected ('outrume','Illegal task value')
end select

end subroutine outrume

