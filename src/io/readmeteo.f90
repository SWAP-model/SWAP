! File VersionID:
!   $Id: readmeteo.f90 372 2018-03-13 10:01:20Z heine003 $
!
!     This file contains the following subroutines, in order of calling:
!     1. ReadMeteo      : reads meteorological input data; called in SWAP
!     2. ReadRainEvents : reads input data on rain events; called in ReadMeteo (optional)

! SUBROUTINE 1.
! ----------------------------------------------------------------------
      subroutine ReadMeteoYear
      use error_mod, only: fatalerr_collected
! ----------------------------------------------------------------------
!     Last modified      : March 2014
!     Purpose            : read meteorological data of one calendar year
! ----------------------------------------------------------------------
      use variables
      use meteodt_mod, only: MeteoDT
      implicit none
    
! --- global
!   - general     

! --- local
      character(len=11)  datedum
      character(len=3)   ext
      character(len=200) filnam
      character(len=800) messag
      integer   datea(6),daynumber,i,ifnd
      integer   getun2,wth      
      real(4)   fsec
      real(8)   etrmax,etrmin,hummax,hummin,radmax,radmin,raimax,raimin
      real(8)   tmeteo,tmnmax,tmnmin,tmxmax,tmxmin,winmax,winmin
      
!========================= Read Meteo file =============================

! --- detection & handling of missing values ---------------------------
!       default ranges of parameter values
      radmin = -1.0d5
      radmax = 5.0d6
      tmnmin = -1.0d5
      tmnmax = 1.0d5
      tmxmin = -1.0d5
      tmxmax = 1.0d5
      hummin = -1.0d5
      hummax = 1.0d5
      winmin = -1.0d5
      winmax = 1.0d5
      raimin = 0.0d0
      raimax = 1.0d3
      etrmin = -1.0d5
      etrmax = 1.0d5
! --- no missing values allowed if penmon must be executed
      if (swetr.eq.0) then
        radmin = 0.0d0
        tmnmin = -50.0d0
        tmnmax = 35.0d0
        tmxmin = -50.0d0
        tmxmax = 60.0d0
        hummin = 0.0d0
        hummax = 10.0d0
        winmin = 0.0d0
        winmax = 150.0d0
      endif
! --- error in case etref missing
      if (swetr.eq.1) then
        etrmin = -0.00001d0
        etrmax = 1.0d2
      endif  
! --- no missing values for tmn and tmx allowed if crop development or
! ---   numerical soil temperatures must be simulated
      if (swmeteo.eq.2 .or. swcalt.eq.2) then
        tmnmin = -60.0d0
        tmnmax = 50.0d0
        tmxmin = -50.0d0
        tmxmax = 60.0d0
      endif
! --- no missing value for rad allowed in case the detailed crop model
! ---   or the grass routine is active 
      if (swmeteo .eq. 2) then
        radmin = 0.0d0
        radmax = 5.0d6
      endif
! --- end of handling missing values -----------------------------------

! --- compose filename meteorological file
      write (ext,'(i3.3)') mod(yearmeteo,1000)
      filnam = trim(pathatm)//trim(metfil)//'.'//trim(ext)


! --- get values from file; two options:
!   - 1 swmetdetail = 0; daily input
!   - 2 swmetdetail = 1; detailed input for nmetdetail time intervals per day

! --- CSV mode: extract year slice from the pre-loaded cache.
!     For daily mode (swmetdetail=0): MeteoCSVYear.
!     For sub-daily mode (swmetdetail=1): MeteoCSVDetYear.
      if (swMetCSV == 1) then
         if (swmetdetail == 0) then
            call MeteoCSVYear(ifnd)
         else
            call MeteoCSVDetYear(ifnd)
         end if
         goto 100
      end if

      if (swmetdetail.eq.0) then
         if (swMetFilAll == 1) then
            call MeteoInOneFile (2, ifnd)
         else
! ---       initialise and start reading
            wth = getun2 (10,90,2)
            call rdinit(wth,logf,filnam)
               call rdacha ('station',station,366,ifnd)
               call rdfinr ('dd',1,31,ad,366,ifnd)
               call rdfinr ('mm',1,12,am,366,ifnd)
               call rdfdor ('rad',radmin,radmax,arad,366,ifnd)
               call rdfdor ('tmin',tmnmin,tmnmax,atmn,366,ifnd)
               call rdfdor ('tmax',tmxmin,tmxmax,atmx,366,ifnd)
               call rdfdor ('hum',hummin,hummax,ahum,366,ifnd)
               call rdfdor ('wind',winmin,winmax,awin,366,ifnd)
               if (swrain.lt.3) call rdfdor('rain',raimin,raimax,arai,366,ifnd)
               call rdfdor ('etref',etrmin,etrmax,aetr,366,ifnd)
               if (swrain.eq.2) call rdfdor('wet',0.d0,1.d0,wet,366,ifnd)
! -            Convert radiation from kj/m2/d to j/m2/d
               arad = 1000.0d0 * arad
! ---       close meteorological file
            close (wth)
         end if
      elseif (swmetdetail.eq.1) then
! ---    initialise and start reading
         wth = getun2 (10,90,2)
         call rdinit(wth,logf,filnam)
            call rdatim ('date',dettime,nmetfile,ifnd)
            call rdfinr ('record',1,nmetdetail,detrecord,nmetfile,ifnd)
            call rdfdor ('rad',radmin,radmax,detrad,nmetfile,ifnd)   
            call rdfdor ('temp',tmxmin,tmxmax,dettav,nmetfile,ifnd)
            call rdfdor ('hum',hummin,hummax,dethum,nmetfile,ifnd)
            call rdfdor ('wind',winmin,winmax,detwind,nmetfile,ifnd)
            call rdfdor ('rain',raimin,raimax,detrain,nmetfile,ifnd)
! -         Convert radiation from kj/m2/d to j/m2/d
            detrad = 1000.0d0 * detrad
! ---    close meteorological file
         close (wth)
      endif

100   continue
!========================= tests and initialization ====================

! --- perform some reliability tests and some initialization
      if (swmetdetail.eq.0) then

! --- determine first and last day numbers        
         datea(1) = yearmeteo
         datea(2) = 1
         datea(3) = 1
         datea(4) = 0
         datea(5) = 0
         datea(6) = 0
         fsec = 0.0
         call dtardp (datea,fsec,timjan1)
         datea(2) = am(1)
         datea(3) = ad(1)
         call dtardp (datea,fsec,tmeteo)
         daynrfirst = nint ( tmeteo - timjan1 + 1.0d0 )
         datea(2) = am(ifnd)
         datea(3) = ad(ifnd)

         call dtardp (datea,fsec,tmeteo)
         daynrlast = nint ( tmeteo - timjan1 + 1.0d0 )

! --- check date 
         do i = 2, ifnd-1
            datea(2) = am(i)
            datea(3) = ad(i)
            call dtardp (datea,fsec,tmeteo)
            daynumber = nint ( tmeteo - timjan1 + 1.0d0 )
            if (daynumber .ne. daynrfirst+i-1) then 
!           wrong date after daynumber i-1
              datea(2) = am(i-1)
              datea(3) = ad(i-1)
              call dtardp (datea,fsec,tmeteo)
              call dtdpst ('year-month-day',tmeteo,datedum)           
              messag ='In meteo file '//trim(filnam)//' the date after '&
     &        //datedum//' is not correct! First adapt meteo file!'
              call fatalerr_collected ('meteo',messag)
            endif
         end do
         
! --- snow and frost calculation conditions require realistic air temperatures
         do i = 1, ifnd
            if ( (swsnow.eq.1.or.swfrost.eq.1) .and.                    &
     &         (atmn(i).lt.-98.9d0.or.atmx(i).lt.-98.9d0) ) then    
             messag ='In meteo file '//trim(filnam)//' temperatures'//  &
     &       ' must be input to calculate snow conditions (SWSNOW=1)'// &
     &       ' and/or frost conditions (SWFROST=1)! adapt meteo file.'
             call fatalerr_collected ('meteo',messag)
           endif
         enddo
! --- realistic air temperatures are required when soil temp. is simulated 
!     using the numerical model  (not if data from other file 
!! requires additional verification)
!      if (swhea.eq.1 .and. swcalt.eq.2) then
!        if (atmn(i).lt.-98.9d0 .or. atmx(i).lt.-98.9d0) then
!            messag ='In meteo file realistic temperatures are '//
!     &        ' required for simulation of temperature profiles (SWHEA=1) using'//
!     &        'num.model (SWCALT=2)!  Adapt meteo data! '
!            call fatalerr_collected ('meteo',messag)
!        endif
!      endif

! --- in case of swrain=2 then make sure that Rain and Wet correspond
         if (swrain.eq.2) then
            do i = 1, ifnd
              if ((arai(i).gt.1.d-10 .and. wet(i).lt.1.d-10) .or.       &
     &          (arai(i).lt.1.d-10 .and. wet(i).gt.1.d-10) ) then
                write(messag,1001) trim(filnam), am(i), ad(i)
 1001           format(' In meteo file ',a,';  month =',i3,';  day =',  &
     &             i3,'; SwRain= 2      ',                              &
     &             ' Rain and Wet donot correspond.  Adapt meteo data!')
               call fatalerr_collected ('ReadMeteoYear',messag)
              endif  
            enddo
         endif
!
! --- in case of swrain = 1 or 2 then store daily precipitation (rain) amount 
!     in rainamount and precipitation time in raintimearray
         if (swrain.eq.1 .or. swrain.eq.2) then
            nmrain = 0
            do i = 1, ifnd
               nmrain   = nmrain + 1
               datea(2) = am(i)
               datea(3) = ad(i)
               call dtardp (datea,fsec,tmeteo)
               raintimearray(i+1) = tmeteo
               rainamount(i)      = arai(i)
            enddo
         endif
! --- end of reliability tests and initialization of daily meteo   

      elseif (swmetdetail.eq.1) then
! --- initialization of detailed meteo
  
! ---   initialize total record number for new weather file
        irectotal = int(t1900-dettime(1)+0.1d0)*nmetdetail

! ---   initialize number of days for running average Tmin
        nofd = 0
      
      endif
! --- end of initialization of detailed meteo   
      
! --- end of reading meteo data file **********************************

! --- close present year for further reading      
      flYearStart = .false.

!========================= Read rain file =============================

      if (swrain .eq. 3) then
! ---   rainfall events are specified
        call ReadRainEvents()
      endif
! --- end of reading rain file ****************************************     

! --- reopen present year for processing rain intensity data at beginning MeteoDt      
      if (swrain .gt. 0) then
         flYearStart = .true.
         call MeteoDT
      endif

      return
      end subroutine ReadMeteoYear


! SUBROUTINE 2.
! ----------------------------------------------------------------------
      subroutine ReadRainEvents()
      use error_mod, only: fatalerr_collected
! ----------------------------------------------------------------------
!     Last modified      : February 2014
!     Purpose            : read rainfall data (events) of one calendar year
!     Interface:
!       I   - logf,yearmeteo,pathatm,rainfil
!       O   - nmrain,rainamount,raintimearray
! ----------------------------------------------------------------------
      use variables, only: logf,yearmeteo,pathatm,rainfil,nmrain,rainamount,raintimearray, &
                           swRainCSV, raincsv_dat, nraincsv
      use swap_array_dimensions, only: mrain

      implicit none

! --- local
      character(len=2)   chday,chmonth
      character(len=6)   chtime
      character(len=3)   ext
      character(len=200) filnam
      character(len=300) messag
      integer   ad(mrain),am(mrain),ay(mrain)
      integer   datea(6),getun2,i,ic,ifnd,pre
      real(4)   fsec,time4(mrain)
      real(8)   tday,tdayold,vsmall
      logical   flrnx
      ! CSV-path locals
      integer  :: jday
      external    jday
      integer, parameter :: jd1900 = 2415020
      real(8)  :: t_jan1, t_dec31, tfrac

      vsmall = 1.0d-8

!========================= CSV path ====================================
      if (swRainCSV == 1) then
         ! Extract current year's events from the pre-loaded cache.
         t_jan1  = real(jday(yearmeteo,  1,  1) - jd1900, 8)
         t_dec31 = real(jday(yearmeteo, 12, 31) - jd1900, 8) + 1.0d0

         ifnd = 0
         do i = 1, nraincsv
            if (raincsv_dat(i,1) >= t_jan1 - 0.5d0 .and. &
     &          raincsv_dat(i,1) <  t_dec31 + 0.5d0) then
               ifnd = ifnd + 1
               raintimearray(ifnd) = raincsv_dat(i,1)
               rainamount(ifnd)    = raincsv_dat(i,2)
            end if
         end do

         if (ifnd == 0) then
            call fatalerr_collected('ReadRainEvents', &
     &         'No rain events CSV records found for the requested year')
            return
         end if

         ! Zero-prepend: if first event is not at midnight, insert t=0 record.
         tfrac = raintimearray(1) - real(int(raintimearray(1)), 8)
         if (tfrac > vsmall) then
            do i = ifnd, 1, -1
               raintimearray(i+1) = raintimearray(i)
               rainamount(i+1)    = rainamount(i)
            end do
            ifnd = ifnd + 1
            raintimearray(1) = real(int(raintimearray(2)), 8)
            rainamount(1)    = 0.0d0
         end if

         ! Deduplication: drop midnight-crossover duplicates.
         ic = 1
         do i = 2, ifnd
            if ((raintimearray(i) - raintimearray(ic)) > vsmall) then
               ic = ic + 1
               raintimearray(ic) = raintimearray(i)
               rainamount(ic)    = rainamount(i)
            end if
         end do
         nmrain = ic

         ! Ascending order check.
         do i = 2, nmrain
            if ((raintimearray(i) - raintimearray(i-1)) .lt. vsmall) then
               messag = 'In rain events CSV file the time of a record ' //  &
     &            'is not greater than its predecessor. Adapt the file!'
               call fatalerr_collected('ReadRainEvents', messag)
            end if
         end do
         return
      end if

!========================= TTutil / legacy .YYY path ===================

      write (ext,'(i3.3)') mod(yearmeteo,1000)
      filnam = trim(pathatm)//trim(rainfil)//'.'//trim(ext)

      pre = getun2 (10,90,2)
      call rdinit(pre,logf,filnam)
      call rdainr ('day',1,31,ad,mrain,ifnd)
      call rdfinr ('month',1,12,am,mrain,ifnd)
      call rdfinr ('year',1,3000,ay,mrain,ifnd)
      call rdfdor ('time',0.0d0,1.0d0,raintimearray,mrain,ifnd)
      call rdfdor ('amount',0.0d0,10000.0d0,rainamount,mrain,ifnd)
      close (pre)

      if (raintimearray(1).gt.vsmall) then
         do i = ifnd, 1, -1
           ay(i+1) = ay(i)
           am(i+1) = am(i)
           ad(i+1) = ad(i)
           raintimearray(i+1) = raintimearray(i)
           rainamount(i+1)    = rainamount(i)
         enddo
         ifnd = ifnd + 1
         raintimearray(i+1) = 0.0d0
         rainamount(i+1)    = 0.0d0
      endif

      fsec = 0.0
      ic = 0
      datea = 0
      do i = 1,ifnd
         datea(1) = ay(i)
         datea(2) = am(i)
         datea(3) = ad(i)
         call dtardp (datea,fsec,tday)
         time4(i)         = real(raintimearray(i))
         raintimearray(i) = tday + raintimearray(i)
         if (ic.gt.0) then
            if ( (raintimearray(i)-raintimearray(ic)) .lt. vsmall .and. &
     &           tday.gt.tdayold) then
               flrnx = .false.
            else
               flrnx = .true.
            endif
         else
            flrnx = .true.
         endif
         if (flrnx) then
            ic = ic + 1
            raintimearray(ic) = raintimearray(i)
            rainamount(ic)    = rainamount(i)
         endif
         tdayold = tday
      end do

      nmrain = ic

      do i = 2,ic
         if ( (raintimearray(i)-raintimearray(i-1)) .lt. vsmall) then
            write(chmonth,'(i2)')  am(i)
            write(chday,'(i2)')    ad(i)
            write(chtime,'(f6.3)') time4(i)
           messag ='Error-message from module ReadRainEvents'//         &
     &        '. In detailed rain file '//trim(filnam)//                &
     &        ' the time of a rainrecord should be greater than the'//  &
     &        ' time of its preceding rainrecord! Adapt rain file.'//   &
     &        '      Month: '//chmonth//'; Day: '//chday//'; Time: '//  &
     &        chtime//'!'
            call fatalerr_collected ('ReadRain',messag)
         endif
      enddo

      return
      end subroutine ReadRainEvents
      
subroutine MeteoInOneFile (iTask, ifnd)
use error_mod, only: fatalerr_collected
use variables, only: arad, atmn, atmx, ahum, awin, arai, aetr, wet, station, ad, am, mayrs, maday, &
                     yearmeteo, tstart, tend, metfil, logf, pathatm
implicit none
! global
integer,           intent(in)             :: iTask
integer,           intent(out)            :: ifnd

! local
integer,           parameter              :: iReadType = 2     ! 1 = use TTutil; 2 = simple read (requires fixed sequence in columns; no value range checking)
integer,                             save :: Nall, nyrs, iyr
integer,           dimension(mayrs), save :: istart, ystart
character(len=10), dimension(maday), save :: all_station
integer,           dimension(maday), save :: all_dd, all_mm, all_yyyy
real(8),           dimension(maday), save :: all_rad, all_tmin, all_tmax, all_hum, all_wind, all_rain, all_etref, all_wet
integer                                   :: iunall, getun, i1, i2, i
character(len=80)                         :: fin
character(len=250)                        :: header

real(8)                                   :: radmin, radmax, tmnmin, tmnmax, tmxmin, tmxmax, hummin, hummax
real(8)                                   :: winmin, winmax, raimin, raimax, etrmin, etrmax
real(8)                                   :: wstart, wend
integer,           dimension(6)           :: datea
real                                      :: fsec

select case (iTask)
case (1)
   ! see ReadMeteoYear
   radmin =   0.0d0;    radmax =   5.0d6
   tmnmin = -50.0d0;    tmnmax =  60.0d0
   tmxmin = -50.0d0;    tmxmax =  60.0d0
   hummin =   0.0d0;    hummax =  10.0d0
   winmin =   0.0d0;    winmax = 150.0d0
   raimin =   0.0d0;    raimax =   1.0d3
   etrmin =  -1.0d5;    etrmax =   1.0d5
   
   ! read all data
   iunall = getun (500, 900)
   fin =  trim(pathatm) // trim(metfil)
   if (iReadType == 1) then
      call RDinit (iunall, logf, fin)
         ! station,dd,mm,yyyy,rad,tmin,tmax,hum,wind,rain,etref,wet
         call RDacha ('station',                 all_station, maday, Nall)
         call RDfinr ('dd',      1,      31,     all_dd,      maday, Nall)
         call RDfinr ('mm',      1,      12,     all_mm,      maday, Nall)
         call RDfint ('yyyy',                    all_yyyy,    maday, Nall)
         call RDfdor ('rad',     radmin, radmax, all_rad,     maday, Nall)
         call RDfdor ('tmin',    tmnmin, tmnmax, all_tmin,    maday, Nall)
         call RDfdor ('tmax',    tmxmin, tmxmax, all_tmax,    maday, Nall)
         call RDfdor ('hum',     hummin, hummax, all_hum,     maday, Nall)
         call RDfdor ('wind',    winmin, winmax, all_wind,    maday, Nall)
         call RDfdor ('rain',    raimin, raimax, all_rain,    maday, Nall)
         call RDfdor ('etref',   etrmin, etrmax, all_etref,   maday, Nall)
         call RDfdor ('wet',     0.d0,   1.d0,   all_wet,     maday, Nall)
      close (iunall)
   else
      open (unit = iunall, file = fin, status='old')
         ! skip comment lines at top of file; comment lines start with !, * or # as first character
         do
            read (iunall,'(A)',end=1) header
            if (header(1:1) /= '!' .AND. header(1:1) /= '*' .AND. header(1:1) /= '#') then
               exit
            end if
         end do
         call lowerc (header)
         if (trim(header) /= 'station,dd,mm,yyyy,rad,tmin,tmax,hum,wind,rain,etref,wet') call fatalerr_collected ('MeteoInOneFile', 'Wrong header line in metfil.met')
         Nall = 1
         do
            read (iunall,*,end=2,err=1) all_station(Nall), all_dd(Nall), all_mm(Nall), all_yyyy(Nall),   &
                                        all_rad(Nall), all_tmin(Nall), all_tmax(Nall), all_hum(Nall),    &
                                        all_wind(Nall), all_rain(Nall), all_etref(Nall), all_wet(Nall)
            Nall = Nall + 1
         end do
1        continue
         ! hopefully we never get here
         call fatalerr_collected ('MeteoInOneFile', 'reading error metfil.met')
2        continue
         Nall = Nall - 1
      close (iunall)
   end if
   
   ! check if start and end time are present in meteofile
   datea = 0; fsec = 0.0
   datea(1) = all_yyyy(1); datea(2) = all_mm(1); datea(3) = all_dd(1)
   call DTARDP (datea, fsec, wstart)
   datea(1) = all_yyyy(Nall); datea(2) = all_mm(Nall); datea(3) = all_dd(Nall)
   call DTARDP (datea, fsec, wend)
   if (wstart > tstart .or. wend < tend) call fatalerr_collected ('MeteoInOneFile', 'wstart > tstart .or. wend < tend')
   
   ! no check on daily step-size in time records nor on data in ascending times
   continue
   
   ! determine and store start positions for new years
   nyrs         = 1
   istart(nyrs) = 1
   ystart(nyrs) = all_yyyy(1)
   do i = 2, Nall
      if (all_yyyy(i) > all_yyyy(i-1)) then
         nyrs         = nyrs + 1
         istart(nyrs) = i
         ystart(nyrs) = all_yyyy(i)
      end if
   end do
   
   ! dummy
   ifnd = 0
   
   ! data read and stored preperly
   write (logf,'(A)') 'metfil.met: all meteo data read and stored.'
   
case (2)
   ! new year: determine where we are
   ! current year: yearmeteo
   do i = 1, nyrs
      if (ystart(i) == yearmeteo) then
         iyr = i
         exit
      end if
   end do
   
   i1 = istart(iyr)
   if (iyr < nyrs) then
      i2  = istart(iyr+1) - 1
   else
      i2  = Nall
   end if
   
   ! return data
   station  = all_station(i1:i2)
   ad       = all_dd(i1:i2)
   am       = all_mm(i1:i2)
   arad     = all_rad(i1:i2)
   atmn     = all_tmin(i1:i2)
   atmx     = all_tmax(i1:i2)
   ahum     = all_hum(i1:i2)
   awin     = all_wind(i1:i2)
   arai     = all_rain(i1:i2)
   aetr     = all_etref(i1:i2)
   wet      = all_wet(i1:i2)
   ! Convert radiation from kj/m2/d to j/m2/d
   arad = 1000.0d0 * arad

   ! number of days 
   ifnd = i2-i1+1

case default
   call fatalerr_collected ('MeteoInOneFile', 'Illegal value for iTask')
end select

return
end subroutine MeteoInOneFile


! SUBROUTINE: MeteoCSVYear
! Extract one year's daily meteo from the pre-loaded metcsv_dat cache.
! Called by ReadMeteoYear when swMetCSV == 1, replacing both the per-year
! file reader and MeteoInOneFile(2).  After return, arad/atmn/atmx/ahum/
! awin/arai/aetr/wet/ad/am are populated exactly as in the legacy path
! so that the validation and rain-array init code in ReadMeteoYear works
! unchanged.
subroutine MeteoCSVYear(ifnd)
use error_mod, only: fatalerr_collected
use variables, only: arad, atmn, atmx, ahum, awin, arai, aetr, wet, ad, am, &
                     yearmeteo, metcsv_dat, nmetcsv, &
                     daynrfirst, daynrlast, timjan1
implicit none
integer, intent(out) :: ifnd

integer  :: i, i1, i2, n
integer  :: datea(6)
real(4)  :: fsec
real(8)  :: t_jan1, t_dec31, tval
! Julian-day epoch used by csv_reader: JD - 2415020
integer, parameter :: jd1900 = 2415020
! Forward declaration: jday is defined later in this file.
integer :: jday
external jday

! Year boundaries in days-since-jd1900
t_jan1  = real(jday(yearmeteo,  1,  1) - jd1900, 8)
t_dec31 = real(jday(yearmeteo, 12, 31) - jd1900, 8)

! Find row range for this year (metcsv_dat is sorted by date).
i1 = 0; i2 = 0
do i = 1, nmetcsv
   if (metcsv_dat(i, 1) >= t_jan1 - 0.5d0 .and. &
       metcsv_dat(i, 1) <= t_dec31 + 0.5d0) then
      if (i1 == 0) i1 = i
      i2 = i
   end if
end do

if (i1 == 0) then
   call fatalerr_collected('MeteoCSVYear', &
      'No meteo CSV records found for the requested year')
   ifnd = 0; return
end if

n = i2 - i1 + 1
ifnd = n

! Populate per-day arrays.
arad(1:n) = metcsv_dat(i1:i2, 2) * 1000.0d0   ! kJ/m2/d → J/m2/d
atmn(1:n) = metcsv_dat(i1:i2, 3)
atmx(1:n) = metcsv_dat(i1:i2, 4)
ahum(1:n) = metcsv_dat(i1:i2, 5)
awin(1:n) = metcsv_dat(i1:i2, 6)
arai(1:n) = metcsv_dat(i1:i2, 7)
aetr(1:n) = metcsv_dat(i1:i2, 8)
wet(1:n)  = metcsv_dat(i1:i2, 9)

! Backfill ad/am from the date column (days-since-jd1900 → month/day).
! ReadMeteoYear's validation code uses am(i)/ad(i) to build raintimearray.
do i = 1, n
   call days1900_to_md(nint(metcsv_dat(i1+i-1, 1)), am(i), ad(i))
end do

! daynrfirst / daynrlast via existing DTARDP, matching ReadMeteoYear logic.
datea = 0; fsec = 0.0
datea(1) = yearmeteo; datea(2) = 1; datea(3) = 1
call dtardp(datea, fsec, t_jan1)
timjan1 = t_jan1

datea(2) = am(1); datea(3) = ad(1)
call dtardp(datea, fsec, tval)
daynrfirst = nint(tval - timjan1 + 1.0d0)

datea(2) = am(n); datea(3) = ad(n)
call dtardp(datea, fsec, tval)
daynrlast = nint(tval - timjan1 + 1.0d0)

end subroutine MeteoCSVYear


! SUBROUTINE: MeteoCSVDetYear
! Extract one year's sub-daily meteo from the pre-loaded metcsv_det cache.
! Called by ReadMeteoYear when swMetCSV==1 and swmetdetail==1.
! Populates dettime, detrecord, detrad, dettav, dethum, detwind, detrain.
! irectotal and nofd are set in ReadMeteoYear after goto 100.
subroutine MeteoCSVDetYear(ifnd)
use error_mod, only: fatalerr_collected
use variables, only: yearmeteo, metcsv_det, nmetcsv_det, &
                     dettime, detrecord, detrad, dettav, dethum, detwind, detrain
use swap_array_dimensions, only: NMETFILE
implicit none
integer, intent(out) :: ifnd

integer, parameter :: jd1900 = 2415020
integer :: jday
external jday

integer  :: i, i1, i2, n
real(8)  :: t_jan1, t_jan1_next

! Year boundaries in days-since-jd1900.
! All sub-daily timestamps for yearmeteo satisfy:
!   t_jan1 <= timestamp < t_jan1_next
t_jan1      = real(jday(yearmeteo,   1, 1) - jd1900, 8)
t_jan1_next = real(jday(yearmeteo+1, 1, 1) - jd1900, 8)

! Scan cache for this year (cache is sorted by datetime).
i1 = 0; i2 = 0
do i = 1, nmetcsv_det
   if (metcsv_det(i,1) >= t_jan1 - 0.5d0 .and. &
       metcsv_det(i,1) <  t_jan1_next - 0.5d0) then
      if (i1 == 0) i1 = i
      i2 = i
   end if
end do

if (i1 == 0) then
   call fatalerr_collected('MeteoCSVDetYear', &
      'No sub-daily meteo CSV records found for the requested year')
   ifnd = 0; return
end if

n = i2 - i1 + 1
if (n > NMETFILE) then
   call fatalerr_collected('MeteoCSVDetYear', &
      'Sub-daily meteo CSV record count exceeds NMETFILE (17568)')
   ifnd = 0; return
end if
ifnd = n

! Populate per-slot arrays.
! metcsv_det columns: 1=datetime, 2=record, 3=rad(kJ), 4=temp, 5=hum, 6=wind, 7=rain
dettime(1:n)   = metcsv_det(i1:i2, 1)
detrecord(1:n) = nint(metcsv_det(i1:i2, 2))
detrad(1:n)    = metcsv_det(i1:i2, 3) * 1000.0d0   ! kJ/m2 → J/m2
dettav(1:n)    = metcsv_det(i1:i2, 4)
dethum(1:n)    = metcsv_det(i1:i2, 5)
detwind(1:n)   = metcsv_det(i1:i2, 6)
detrain(1:n)   = metcsv_det(i1:i2, 7)

end subroutine MeteoCSVDetYear


! Helper: convert days-since-jd1900 to (month, day) via inverse Julian Day.
subroutine days1900_to_md(d1900, mm, dd)
implicit none
integer, intent(in)  :: d1900
integer, intent(out) :: mm, dd
integer :: jd, a, b, c, d, e, m
integer, parameter :: jd1900 = 2415020
jd = d1900 + jd1900
a = jd + 32044
b = (4*a + 3) / 146097
c = a - (146097*b) / 4
d = (4*c + 3) / 1461
e = c - (1461*d) / 4
m = (5*e + 2) / 153
dd = e - (153*m + 2)/5 + 1
mm = m + 3 - 12*(m/10)
end subroutine days1900_to_md


! Pure Julian Day Number for use inside MeteoCSVYear (avoids dependency on
! csv_reader_mod's private julian_day).
pure function jday(y, m, d) result(jd)
integer, intent(in) :: y, m, d
integer :: jd, a, yy, mm
a  = (14 - m) / 12
yy = y + 4800 - a
mm = m + 12*a - 3
jd = d + (153*mm + 2)/5 + 365*yy + yy/4 - yy/100 + yy/400 - 32045
end function jday
