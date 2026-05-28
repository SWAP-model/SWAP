! File VersionID:
!   $Id: readmeteo.f90 372 2018-03-13 10:01:20Z heine003 $
!     This file contains the following subroutines, in order of calling:
!     1. ReadMeteo      : reads meteorological input data; called in SWAP
!     2. ReadRainEvents : reads input data on rain events; called in ReadMeteo (optional)

! SUBROUTINE 1.
! ----------------------------------------------------------------------
      subroutine ReadMeteoYear(state, config)
      use error_mod, only: fatalerr_collected
! ----------------------------------------------------------------------
!     Last modified      : March 2014
!     Purpose            : read meteorological data of one calendar year
! ----------------------------------------------------------------------
      use meteodt_mod, only: MeteoDT
      use swap_state_mod, only: swap_state_t
      use swap_config_mod, only: swap_config_t
      use meteo_buffer_mod, only: get_meteo_mode, METEO_MODE_EXTERNAL_BUFFER
      implicit none

      type(swap_state_t), intent(inout) :: state
        !! Simulation state (passed through to MeteoDT/reduceva for atmosphere dual-writes)
      type(swap_config_t), intent(in) :: config
        !! Simulation configuration (switches: swetr/swrain/swmetdetail/swsnow/swfrost/swcalt)
    
! --- global
!   - general     

! --- local
      character(len=11)  datedum
      character(len=3)   ext
      character(len=200) filnam
      character(len=800) messag
      integer   datea(6),daynumber,i,ifnd
      real(4)   fsec
      real(8)   etrmax,etrmin,hummax,hummin,radmax,radmin,raimax,raimin
      real(8)   tmeteo,tmnmax,tmnmin,tmxmax,tmxmin,winmax,winmin

      associate (time => state%timecontrol,   &
                 atmo => state%atmosphere)

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
      if (config%meteo%swetr.eq.0) then
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
      if (config%meteo%swetr.eq.1) then
        etrmin = -0.00001d0
        etrmax = 1.0d2
      endif
! --- no missing values for tmn and tmx allowed if crop development or
! ---   numerical soil temperatures must be simulated
      if (time%swmeteo.eq.2 .or. config%heat%swcalt.eq.2) then
        tmnmin = -60.0d0
        tmnmax = 50.0d0
        tmxmin = -50.0d0
        tmxmax = 60.0d0
      endif
! --- no missing value for rad allowed in case the detailed crop model
! ---   or the grass routine is active
      if (time%swmeteo .eq. 2) then
        radmin = 0.0d0
        radmax = 5.0d6
      endif
! --- end of handling missing values -----------------------------------

! --- compose filename meteorological file
      write (ext,'(i3.3)') mod(time%yearmeteo,1000)
      filnam = trim(config%general%pathatm)//trim(config%meteo%metfile)//'.'//trim(ext)


! --- get values from file; two options:
!   - 1 swmetdetail = 0; daily input
!   - 2 swmetdetail = 1; detailed input for nmetdetail time intervals per day

! --- SS-BMI2: external buffer mode bypasses the CSV reader (daily path only).
!     Sub-daily / detailed meteo always stays on CSV in Phase 2.
      if (config%meteo%swmetdetail == 0 .and. get_meteo_mode() == METEO_MODE_EXTERNAL_BUFFER) then
         call read_meteo_from_external_buffer_year(ifnd, state)
! --- CSV mode is the only supported path otherwise (ADR 0014).
!     Daily mode → MeteoCSVYear. Sub-daily → MeteoCSVDetYear.
      else if (config%meteo%swmetdetail == 0) then
         call MeteoCSVYear(ifnd, state)
      else
         call MeteoCSVDetYear(ifnd, state)
      end if

!========================= tests and initialization ====================

! --- perform some reliability tests and some initialization
      if (config%meteo%swmetdetail.eq.0) then

! --- determine first and last day numbers
         datea(1) = time%yearmeteo
         datea(2) = 1
         datea(3) = 1
         datea(4) = 0
         datea(5) = 0
         datea(6) = 0
         fsec = 0.0
         call dtardp (datea,fsec,time%timjan1)
         datea(2) = atmo%am(1)
         datea(3) = atmo%ad(1)
         call dtardp (datea,fsec,tmeteo)
         atmo%daynrfirst = nint ( tmeteo - time%timjan1 + 1.0d0 )
         datea(2) = atmo%am(ifnd)
         datea(3) = atmo%ad(ifnd)

         call dtardp (datea,fsec,tmeteo)
         atmo%daynrlast = nint ( tmeteo - time%timjan1 + 1.0d0 )

! --- check date 
         do i = 2, ifnd-1
            datea(2) = atmo%am(i)
            datea(3) = atmo%ad(i)
            call dtardp (datea,fsec,tmeteo)
            daynumber = nint ( tmeteo - time%timjan1 + 1.0d0 )
            if (daynumber .ne. atmo%daynrfirst+i-1) then
!           wrong date after daynumber i-1
              datea(2) = atmo%am(i-1)
              datea(3) = atmo%ad(i-1)
              call dtardp (datea,fsec,tmeteo)
              call dtdpst ('year-month-day',tmeteo,datedum)           
              messag ='In meteo file '//trim(filnam)//' the date after '&
     &        //datedum//' is not correct! First adapt meteo file!'
              call fatalerr_collected ('meteo',messag)
            endif
         end do
         
! --- snow and frost calculation conditions require realistic air temperatures
         do i = 1, ifnd
            if ( (config%meteo%snow%swsnow.eq.1.or.config%soil%frost%swfrost.eq.1) .and.   &
     &         (atmo%atmn(i).lt.-98.9d0.or.atmo%atmx(i).lt.-98.9d0) ) then
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
         if (config%meteo%swrain.eq.2) then
            do i = 1, ifnd
              if ((atmo%arai(i).gt.1.d-10 .and. atmo%wet(i).lt.1.d-10) .or.  &
     &          (atmo%arai(i).lt.1.d-10 .and. atmo%wet(i).gt.1.d-10) ) then
                write(messag,1001) trim(filnam), atmo%am(i), atmo%ad(i)
 1001           format(' In meteo file ',a,';  month =',i3,';  day =',  &
     &             i3,'; SwRain= 2      ',                              &
     &             ' Rain and Wet donot correspond.  Adapt meteo data!')
               call fatalerr_collected ('ReadMeteoYear',messag)
              endif  
            enddo
         endif
! --- in case of swrain = 1 or 2 then store daily precipitation (rain) amount
!     in rainamount and precipitation time in raintimearray
! Store daily precipitation amount in rainamount and timing in raintimearray.
         if (config%meteo%swrain.eq.1 .or. config%meteo%swrain.eq.2) then
            atmo%nmrain = 0
            do i = 1, ifnd
               atmo%nmrain = atmo%nmrain + 1
               datea(2) = atmo%am(i)
               datea(3) = atmo%ad(i)
               call dtardp (datea,fsec,tmeteo)
               atmo%raintimearray(i+1) = tmeteo
               atmo%rainamount(i)      = atmo%arai(i)
            enddo
         endif
! --- end of reliability tests and initialization of daily meteo   

      elseif (config%meteo%swmetdetail.eq.1) then
! --- initialization of detailed meteo

! ---   initialize total record number for new weather file
        atmo%irectotal = int(time%t1900 - atmo%dettime(1) + 0.1d0) * config%meteo%nmetdetail

! ---   initialize number of days for running average Tmin
        atmo%nofd = 0

      endif
! --- end of initialization of detailed meteo   
      
! --- end of reading meteo data file **********************************

! --- close present year for further reading
      time%flYearStart = .false.

!========================= Read rain file =============================

      if (config%meteo%swrain .eq. 3) then
! ---   rainfall events are specified
        call ReadRainEvents(state, config)
      endif
! --- end of reading rain file ****************************************

! --- reopen present year for processing rain intensity data at beginning MeteoDt
      if (config%meteo%swrain .gt. 0) then
         time%flYearStart = .true.
         call MeteoDT(state)
      endif

      end associate

      return
      end subroutine ReadMeteoYear


! SUBROUTINE 2.
! ----------------------------------------------------------------------
      subroutine ReadRainEvents(state, config)
      use error_mod, only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      use swap_config_mod, only: swap_config_t
! ----------------------------------------------------------------------
!     Last modified      : February 2014
!     Purpose            : read rainfall data (events) of one calendar year
!     Interface:
!       I   - logf,yearmeteo (via state%timecontrol),pathatm,raincsv_dat,nraincsv
!       O   - nmrain,rainamount,raintimearray
! ----------------------------------------------------------------------
      implicit none
      ! changed to inout for rain timing dual-writes
      type(swap_state_t), intent(inout) :: state
      type(swap_config_t), intent(in) :: config
      integer :: yearmeteo   ! local copy

! --- local
      character(len=300) messag
      integer   i,ic,ifnd
      real(8)   vsmall
      real(8)   tfrac

      vsmall = 1.0d-8
      yearmeteo = state%timecontrol%yearmeteo

!========================= CSV path (only supported, ADR 0014) =========
      ! Extract current year's events from the pre-loaded cache.
      ! Store daily precipitation amount in rainamount and timing in raintimearray.
      ifnd = 0
      block
         integer :: i1, i2, k
         call state%atmosphere%rain_events%year_window(yearmeteo, i1, i2)
         if (i1 > 0) then
            do k = i1, i2
               ifnd = ifnd + 1
               state%atmosphere%raintimearray(ifnd) = state%atmosphere%rain_events%rows(k)%datetime
               state%atmosphere%rainamount(ifnd)    = state%atmosphere%rain_events%rows(k)%amount
            end do
         end if
      end block

      if (ifnd == 0) then
         call fatalerr_collected('ReadRainEvents', &
     &      'No rain events CSV records found for the requested year')
         return
      end if

      ! Zero-prepend: if first event is not at midnight, insert t=0 record.
      tfrac = state%atmosphere%raintimearray(1) - real(int(state%atmosphere%raintimearray(1)), 8)
      if (tfrac > vsmall) then
         do i = ifnd, 1, -1
            state%atmosphere%raintimearray(i+1) = state%atmosphere%raintimearray(i)
            state%atmosphere%rainamount(i+1)    = state%atmosphere%rainamount(i)
         end do
         ifnd = ifnd + 1
         state%atmosphere%raintimearray(1) = real(int(state%atmosphere%raintimearray(2)), 8)
         state%atmosphere%rainamount(1)    = 0.0d0
      end if

      ! Deduplication: drop midnight-crossover duplicates.
      ic = 1
      do i = 2, ifnd
         if ((state%atmosphere%raintimearray(i) - state%atmosphere%raintimearray(ic)) > vsmall) then
            ic = ic + 1
            state%atmosphere%raintimearray(ic) = state%atmosphere%raintimearray(i)
            state%atmosphere%rainamount(ic)    = state%atmosphere%rainamount(i)
         end if
      end do
      state%atmosphere%nmrain = ic

      ! Ascending order check.
      do i = 2, state%atmosphere%nmrain
         if ((state%atmosphere%raintimearray(i) - state%atmosphere%raintimearray(i-1)) .lt. vsmall) then
            messag = 'In rain events CSV file the time of a record ' //  &
     &         'is not greater than its predecessor. Adapt the file!'
            call fatalerr_collected('ReadRainEvents', messag)
         end if
      end do

      return
      end subroutine ReadRainEvents


! SUBROUTINE: read_meteo_from_external_buffer_year
! SS-BMI2 Task 7: populate per-day meteo arrays from the externally-supplied
! buffer when mode == METEO_MODE_EXTERNAL_BUFFER. Mirrors MeteoCSVYear output
! contract: arad/atmn/atmx/ahum/awin/arai/aetr/ad/am/daynrfirst/daynrlast/
! timjan1/ifnd are all set so ReadMeteoYear's post-call validation is unaffected.
! Canonical buffer column order (1-based):
!   1: date (days since JD 1900, csv_reader epoch — see meteo_daily_table_t)
!   2: rain        (mm/d  → arai)
!   3: tmin        (deg C → atmn)
!   4: tmax        (deg C → atmx)
!   5: et_ref      (mm/d  → aetr; -99.9 = compute internally)
!   6: radiation   (kJ/m2/d → arad, converted to J/m2/d by ×1000)
!   7: vapor       (kPa   → ahum)
!   8: wind        (m/s   → awin)
! wet(:) is set to -99.9 (missing) — caller can set swrain≠2 or supply a
! separate wet flag column in a future extension.
subroutine read_meteo_from_external_buffer_year(ifnd, state)
use meteo_buffer_mod, only: get_external_meteo_value, get_external_meteo_n_days, &
                            get_external_meteo_n_cols
use swap_state_mod,   only: swap_state_t
use swap_array_dimensions, only: NMETFILE
implicit none
integer,             intent(out)   :: ifnd
type(swap_state_t),  intent(inout) :: state

integer,  parameter :: jd1900 = 2415020
integer             :: i, i1, i2, n, n_buf
real(8)             :: t_jan1, t_dec31, tval, col1_val
integer             :: datea(6)
real(4)             :: fsec
integer             :: jday
external               jday

associate (atmo => state%atmosphere, &
           time => state%timecontrol)

! Year boundaries in days-since-jd1900 (csv_reader epoch — see meteo_daily_table_t)
t_jan1  = real(jday(time%yearmeteo,  1,  1) - jd1900, 8)
t_dec31 = real(jday(time%yearmeteo, 12, 31) - jd1900, 8)

n_buf = get_external_meteo_n_days()

! Find row range in the buffer for yearmeteo.
i1 = 0; i2 = 0
do i = 1, n_buf
   col1_val = get_external_meteo_value(i, 1)
   if (col1_val >= t_jan1 - 0.5d0 .and. col1_val <= t_dec31 + 0.5d0) then
      if (i1 == 0) i1 = i
      i2 = i
   end if
end do

if (i1 == 0) then
   ! No rows for this year — return 0; ReadMeteoYear will catch via validation.
   ifnd = 0
   return
end if

n = i2 - i1 + 1
ifnd = min(n, NMETFILE)

! Populate per-day arrays (canonical buffer column order).
do i = 1, ifnd
   atmo%arai(i) = get_external_meteo_value(i1+i-1, 2)           ! rain (mm/d)
   atmo%atmn(i) = get_external_meteo_value(i1+i-1, 3)           ! tmin (deg C)
   atmo%atmx(i) = get_external_meteo_value(i1+i-1, 4)           ! tmax (deg C)
   atmo%aetr(i) = get_external_meteo_value(i1+i-1, 5)           ! et_ref (mm/d)
   atmo%arad(i) = get_external_meteo_value(i1+i-1, 6) * 1000.0d0  ! kJ→J /m2/d
   atmo%ahum(i) = get_external_meteo_value(i1+i-1, 7)           ! vapor (kPa)
   atmo%awin(i) = get_external_meteo_value(i1+i-1, 8)           ! wind (m/s)
   atmo%wet(i)  = -99.9d0                                        ! missing — no wet-flag column yet
end do

! Backfill ad/am from the date column (days-since-jd1900 → month/day).
do i = 1, ifnd
   call days1900_to_md(nint(get_external_meteo_value(i1+i-1, 1)), atmo%am(i), atmo%ad(i))
end do

! daynrfirst / daynrlast and timjan1 — mirror MeteoCSVYear logic.
datea = 0; fsec = 0.0
datea(1) = time%yearmeteo; datea(2) = 1; datea(3) = 1
call dtardp(datea, fsec, t_jan1)
time%timjan1 = t_jan1

datea(2) = atmo%am(1); datea(3) = atmo%ad(1)
call dtardp(datea, fsec, tval)
atmo%daynrfirst = nint(tval - time%timjan1 + 1.0d0)

datea(2) = atmo%am(ifnd); datea(3) = atmo%ad(ifnd)
call dtardp(datea, fsec, tval)
atmo%daynrlast = nint(tval - time%timjan1 + 1.0d0)

end associate
end subroutine read_meteo_from_external_buffer_year


! SUBROUTINE: MeteoCSVYear
! Extract one year's daily meteo from the typed meteo table.
! Called unconditionally by ReadMeteoYear (the only supported daily path
! after ADR 0014). After return, arad/atmn/atmx/
! ahum/awin/arai/aetr/wet/ad/am are populated so that the validation and
! rain-array init code in ReadMeteoYear works unchanged.
subroutine MeteoCSVYear(ifnd, state)
use error_mod, only: fatalerr_collected
use swap_state_mod, only: swap_state_t
implicit none
integer, intent(out) :: ifnd
type(swap_state_t), intent(inout) :: state

integer  :: i, i1, i2, n
integer  :: datea(6)
real(4)  :: fsec
real(8)  :: t_jan1, tval

associate (atmo => state%atmosphere, &
           time => state%timecontrol)

call atmo%meteo%year_window(time%yearmeteo, i1, i2)

if (i1 == 0) then
   call fatalerr_collected('MeteoCSVYear', &
      'No meteo CSV records found for the requested year')
   ifnd = 0; return
end if

n = i2 - i1 + 1
ifnd = n

! Populate per-day arrays from typed rows.
do i = 1, n
   atmo%arad(i) = atmo%meteo%rows(i1+i-1)%rad * 1000.0d0   ! kJ/m2/d → J/m2/d
   atmo%atmn(i) = atmo%meteo%rows(i1+i-1)%tmin
   atmo%atmx(i) = atmo%meteo%rows(i1+i-1)%tmax
   atmo%ahum(i) = atmo%meteo%rows(i1+i-1)%hum
   atmo%awin(i) = atmo%meteo%rows(i1+i-1)%wind
   atmo%arai(i) = atmo%meteo%rows(i1+i-1)%rain
   atmo%aetr(i) = atmo%meteo%rows(i1+i-1)%etref
   atmo%wet(i)  = atmo%meteo%rows(i1+i-1)%wet
end do

! Backfill ad/am from the date column.
do i = 1, n
   call days1900_to_md(nint(atmo%meteo%rows(i1+i-1)%date), atmo%am(i), atmo%ad(i))
end do

! daynrfirst / daynrlast (unchanged logic).
datea = 0; fsec = 0.0
datea(1) = time%yearmeteo; datea(2) = 1; datea(3) = 1
call dtardp(datea, fsec, t_jan1)
time%timjan1 = t_jan1

datea(2) = atmo%am(1); datea(3) = atmo%ad(1)
call dtardp(datea, fsec, tval)
atmo%daynrfirst = nint(tval - time%timjan1 + 1.0d0)

datea(2) = atmo%am(n); datea(3) = atmo%ad(n)
call dtardp(datea, fsec, tval)
atmo%daynrlast = nint(tval - time%timjan1 + 1.0d0)

end associate
end subroutine MeteoCSVYear


! SUBROUTINE: MeteoCSVDetYear
! Extract one year's sub-daily meteo from the typed meteo_detail table.
! Called by ReadMeteoYear when swmetdetail==1 (the only supported sub-daily
! path after ADR 0014).
! Populates dettime, detrecord, detrad, dettav, dethum, detwind, detrain.
! irectotal and nofd are set by ReadMeteoYear after this returns.
subroutine MeteoCSVDetYear(ifnd, state)
use error_mod, only: fatalerr_collected
use swap_state_mod, only: swap_state_t
use swap_array_dimensions, only: NMETFILE
implicit none
integer, intent(out) :: ifnd
type(swap_state_t), intent(inout) :: state

integer  :: i, i1, i2, n

associate (atmo => state%atmosphere, &
           time => state%timecontrol)

call atmo%meteo_detail%year_window(time%yearmeteo, i1, i2)

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

! Allocate state-side detail arrays on first use (NMETFILE-sized for parity).
if (.not. allocated(atmo%dettime)) then
   allocate(atmo%dettime(NMETFILE))
   allocate(atmo%detrecord(NMETFILE))
   allocate(atmo%detrad(NMETFILE))
   allocate(atmo%dettav(NMETFILE))
   allocate(atmo%dethum(NMETFILE))
   allocate(atmo%detwind(NMETFILE))
   allocate(atmo%detrain(NMETFILE))
end if

do i = 1, n
   atmo%dettime(i)   = atmo%meteo_detail%rows(i1+i-1)%datetime
   atmo%detrecord(i) = atmo%meteo_detail%rows(i1+i-1)%record
   atmo%detrad(i)    = atmo%meteo_detail%rows(i1+i-1)%rad * 1000.0d0   ! kJ → J
   atmo%dettav(i)    = atmo%meteo_detail%rows(i1+i-1)%temp
   atmo%dethum(i)    = atmo%meteo_detail%rows(i1+i-1)%hum
   atmo%detwind(i)   = atmo%meteo_detail%rows(i1+i-1)%wind
   atmo%detrain(i)   = atmo%meteo_detail%rows(i1+i-1)%rain
end do

end associate
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
