!> Canonical implementations of TTutil utility functions.
!!
!! Provides free subroutines and functions matching the TTutil-era
!! signatures that are still referenced by src/ after Phases A-D of
!! the TTutil retirement arc. With TTutil retired (ADR 0023), this
!! file is the canonical source; no subproject required.
!!
!! Functions / subroutines implemented here (all free, no module):
!!
!!   LOGICAL FUNCTION DTLEAP(YEAR)
!!   SUBROUTINE DTDPAR(DPDTTM, DATEA, FSEC)
!!   SUBROUTINE DTARDP(DATEA, FSEC, DPDTTM)
!!   SUBROUTINE DTDPST(FORM, DPDTTM, STRNG)
!!   SUBROUTINE DTNOW(DATEA)
!!   SUBROUTINE LOWERC(STRING)
!!   SUBROUTINE UPPERC(STRING)
!!   SUBROUTINE ADDSTR(STRING, SIGLEN, TMP)
!!   SUBROUTINE WORDS(RECORD, ILW, SEPARS, IWBEG, IWEND, IFND)
!!   SUBROUTINE DECREA(IWAR, STRING, VALUE)
!!   INTEGER FUNCTION IFINDI(ILIS, ILDEC, IST, IEND, IINP)
!!
!! Date representation: DPDTTM counts days since 1900-01-01 00:00
!! (1900-01-01 = 1.0), using the TTUTIL OFFSET of 693594 absolute
!! days from 0001-01-01. DATEA(6) = [year, month, day, hour, min, sec].

! ============================================================================
! DTLEAP — leap-year predicate (mirrors TTutil DTLEAP.FOR)
! ============================================================================
LOGICAL FUNCTION DTLEAP(YEAR)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: YEAR
   DTLEAP = MOD(YEAR, 4) == 0
   IF (MOD(YEAR, 100) == 0) DTLEAP = .FALSE.
   IF (MOD(YEAR, 400) == 0) DTLEAP = .TRUE.
END FUNCTION DTLEAP

! ============================================================================
! DTARDP — date-array + fractional-seconds -> double-precision day number
! (mirrors TTutil DTARDP.FOR + DTSYS.FOR task 1)
! ============================================================================
SUBROUTINE DTARDP(DATEA, FSEC, DPDTTM)
   IMPLICIT NONE
   INTEGER,          INTENT(IN)  :: DATEA(6)
   REAL,             INTENT(IN)  :: FSEC
   DOUBLE PRECISION, INTENT(OUT) :: DPDTTM

   INTEGER, PARAMETER :: OFFSET = 693594
   INTEGER, PARAMETER :: IMOTB1(12) = [0,31,59,90,120,151,181,212,243,273,304,334]
   INTEGER, PARAMETER :: IMOTB2(12) = [31,28,31,30,31,30,31,31,30,31,30,31]

   LOGICAL :: DTLEAP
   INTEGER :: LASTY, INDAY, DAYMAX
   DOUBLE PRECISION :: DTMP

   ! Days up to end of previous year
   IF (DATEA(1) >= 1) THEN
      LASTY = DATEA(1) - 1
      INDAY = LASTY*365 + LASTY/4 - LASTY/100 + LASTY/400 - 1
   ELSE
      INDAY = OFFSET
   END IF

   ! Add months (cumulative day offset for month)
   IF (DATEA(2) >= 2 .AND. DATEA(2) <= 12) THEN
      INDAY = INDAY + IMOTB1(DATEA(2))
      IF (DTLEAP(DATEA(1)) .AND. DATEA(2) > 2) INDAY = INDAY + 1
   END IF

   ! Add day of month
   IF (DATEA(3) >= 1) THEN
      DAYMAX = IMOTB2(DATEA(2))
      IF (DTLEAP(DATEA(1)) .AND. DATEA(2) == 2) DAYMAX = DAYMAX + 1
      INDAY = INDAY + DATEA(3)
   END IF

   ! Subtract TTUTIL offset to get days-since-1900
   DTMP = DBLE(INDAY - OFFSET)

   ! Add time-of-day components
   DTMP = DTMP + DBLE(DATEA(4)) / 24.0D0
   DTMP = DTMP + DBLE(DATEA(5)) / 1440.0D0
   DTMP = DTMP + DBLE(DATEA(6)) / 86400.0D0
   IF (FSEC >= 0.0 .AND. FSEC < 1.0) &
      DTMP = DTMP + DBLE(FSEC) / 86400.0D0

   DPDTTM = DTMP
END SUBROUTINE DTARDP

! ============================================================================
! DTDPAR — double-precision day number -> date-array + fractional-seconds
! (mirrors TTutil DTDPAR.FOR + DTSYS.FOR task 2)
! ============================================================================
SUBROUTINE DTDPAR(DPDTTM, DATEA, FSEC)
   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(IN)  :: DPDTTM
   INTEGER,          INTENT(OUT) :: DATEA(6)
   REAL,             INTENT(OUT) :: FSEC

   INTEGER, PARAMETER :: OFFSET = 693594
   INTEGER, PARAMETER :: IMOTB2(12) = [31,28,31,30,31,30,31,31,30,31,30,31]

   LOGICAL :: DTLEAP
   INTEGER :: INDAY, ISECS, NDAYS
   INTEGER :: N400, N100, N4, N1, LEFT, LASTY, DLAST
   DOUBLE PRECISION :: FRDAY, DSECS, FSECS

   ! Split into integer days and fractional day
   IF (DPDTTM >= 0.0D0) THEN
      INDAY = INT(DPDTTM)
      FRDAY = DPDTTM - DBLE(INDAY)
   ELSE
      INDAY = -INT(ABS(DPDTTM))
      FRDAY = ABS(DPDTTM) - DBLE(ABS(INDAY))
      IF (FRDAY > 0.0D0) THEN
         FRDAY = 1.0D0 - FRDAY
         INDAY = INDAY - 1
      END IF
   END IF

   ! Fractional day -> integer seconds + sub-second remainder
   DSECS = 86400.0D0 * FRDAY
   ISECS = INT(DSECS)
   FSECS = DSECS - DBLE(ISECS)

   ! Round near-1 fractional seconds up
   IF (FSECS > 0.999999D0) THEN
      FSECS = 0.0D0
      ISECS = ISECS + 1
      IF (ISECS == 86400) THEN
         ISECS = 0
         INDAY = INDAY + 1
      END IF
   END IF

   ! Apply TTUTIL offset to get absolute day from 0001-01-01
   INDAY = INDAY + OFFSET

   ! Decompose year using 400/100/4/1-year period algorithm
   N400 = INDAY / 146097
   LEFT = INDAY - N400 * 146097

   N100 = MIN(3, LEFT / 36524)
   LEFT = LEFT - N100 * 36524

   N4   = LEFT / 1461
   LEFT = LEFT - N4 * 1461

   N1   = MIN(3, LEFT / 365)

   DATEA(1) = N400*400 + N100*100 + N4*4 + N1 + 1

   ! Find DOY from year base
   LASTY = DATEA(1) - 1
   DLAST = LASTY*365 + LASTY/4 - LASTY/100 + LASTY/400 - 1
   INDAY = INDAY - DLAST

   ! Iterate months
   DATEA(2) = 1
   NDAYS = IMOTB2(DATEA(2))
   IF (DTLEAP(DATEA(1)) .AND. DATEA(2) == 2) NDAYS = NDAYS + 1

   DO WHILE (INDAY > NDAYS)
      INDAY    = INDAY - NDAYS
      DATEA(2) = DATEA(2) + 1
      NDAYS    = IMOTB2(DATEA(2))
      IF (DTLEAP(DATEA(1)) .AND. DATEA(2) == 2) NDAYS = NDAYS + 1
   END DO

   DATEA(3) = INDAY
   DATEA(4) = MOD(ISECS, 86400) / 3600
   DATEA(5) = MOD(ISECS, 3600)  / 60
   DATEA(6) = MOD(ISECS, 60)
   FSEC     = REAL(FSECS)
END SUBROUTINE DTDPAR

! ============================================================================
! DTDPST — format double-precision day number as a date string
!!
!! Supported format tokens (case-insensitive, applied longest-first):
!!   YEAR (4-digit), MONTHLT (full month name), MONTHST (3-letter),
!!   MONTH (2-digit), DAY (2-digit), HOUR (2-digit),
!!   MINUTE / MIN (2-digit), SECONDS / SEC (2-digit).
!! All other characters in FORM are copied verbatim.
! ============================================================================
SUBROUTINE DTDPST(FORM, DPDTTM, STRNG)
   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(IN)  :: FORM
   DOUBLE PRECISION, INTENT(IN)  :: DPDTTM
   CHARACTER(LEN=*), INTENT(OUT) :: STRNG

   INTEGER  :: DATEA(6)
   REAL     :: FSEC
   CHARACTER(LEN=80) :: LFORM
   INTEGER :: I, IC, ILF

   CHARACTER(LEN=3), PARAMETER :: MONSHORT(12) = &
      ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
   CHARACTER(LEN=9), PARAMETER :: MONLONG(12) = &
      ['January  ','February ','March    ','April    ','May      ', &
       'June     ','July     ','August   ','September','October  ', &
       'November ','December ']

   ! Decompose day-number to date components
   CALL DTDPAR(DPDTTM, DATEA, FSEC)

   ! Make uppercase copy of format (max 80 chars)
   LFORM = ' '
   ILF = LEN_TRIM(FORM)
   IF (ILF == 0 .OR. ILF > LEN(LFORM)) THEN
      STRNG = '?bad-format?'
      RETURN
   END IF
   LFORM(1:ILF) = FORM(1:ILF)
   DO I = 1, ILF
      IC = ICHAR(LFORM(I:I))
      IF (IC >= 97 .AND. IC <= 122) LFORM(I:I) = CHAR(IC - 32)
   END DO

   ! Start with the uppercased format
   STRNG = LFORM(1:ILF)

   ! Replace tokens longest-first to avoid partial matches:
   ! MONTHLT before MONTHST before MONTH; SECONDS before SEC; MINUTE before MIN
   CALL DTDPST_REPLACE(STRNG, 'MONTHLT', TRIM(MONLONG(DATEA(2))))
   CALL DTDPST_REPLACE(STRNG, 'MONTHST', MONSHORT(DATEA(2)))
   CALL DTDPST_FMT4   (STRNG, 'YEAR',    DATEA(1))
   CALL DTDPST_FMT2   (STRNG, 'MONTH',   DATEA(2))
   CALL DTDPST_FMT2   (STRNG, 'DAY',     DATEA(3))
   CALL DTDPST_FMT2   (STRNG, 'HOUR',    DATEA(4))
   CALL DTDPST_FMT2   (STRNG, 'MINUTE',  DATEA(5))
   CALL DTDPST_FMT2   (STRNG, 'MIN',     DATEA(5))
   CALL DTDPST_FMT2   (STRNG, 'SECONDS', DATEA(6))
   CALL DTDPST_FMT2   (STRNG, 'SEC',     DATEA(6))

CONTAINS

   SUBROUTINE DTDPST_REPLACE(S, TOKEN, VAL)
      CHARACTER(LEN=*), INTENT(INOUT) :: S
      CHARACTER(LEN=*), INTENT(IN)    :: TOKEN, VAL
      INTEGER :: P, LTOK, LVAL
      LTOK = LEN(TOKEN)
      LVAL = LEN_TRIM(VAL)
      P = INDEX(S, TOKEN)
      IF (P > 0) THEN
         S = S(1:P-1) // VAL(1:LVAL) // S(P+LTOK:)
      END IF
   END SUBROUTINE DTDPST_REPLACE

   SUBROUTINE DTDPST_FMT2(S, TOKEN, IVAL)
      CHARACTER(LEN=*), INTENT(INOUT) :: S
      CHARACTER(LEN=*), INTENT(IN)    :: TOKEN
      INTEGER,          INTENT(IN)    :: IVAL
      CHARACTER(LEN=8) :: TMP
      WRITE(TMP, '(I2.2)') IVAL
      CALL DTDPST_REPLACE(S, TOKEN, TMP)
   END SUBROUTINE DTDPST_FMT2

   SUBROUTINE DTDPST_FMT4(S, TOKEN, IVAL)
      CHARACTER(LEN=*), INTENT(INOUT) :: S
      CHARACTER(LEN=*), INTENT(IN)    :: TOKEN
      INTEGER,          INTENT(IN)    :: IVAL
      CHARACTER(LEN=8) :: TMP
      WRITE(TMP, '(I4.4)') IVAL
      CALL DTDPST_REPLACE(S, TOKEN, TMP)
   END SUBROUTINE DTDPST_FMT4

END SUBROUTINE DTDPST

! ============================================================================
! DTNOW — fill DATEA(6) with the current wall-clock date and time.
! (mirrors TTutil dtnow.f90)
! ============================================================================
SUBROUTINE DTNOW(DATEA)
   IMPLICIT NONE
   INTEGER, INTENT(OUT) :: DATEA(6)
   INTEGER :: VALUES(8)
   CALL DATE_AND_TIME(VALUES=VALUES)
   DATEA(1) = VALUES(1)  ! year
   DATEA(2) = VALUES(2)  ! month
   DATEA(3) = VALUES(3)  ! day
   DATEA(4) = VALUES(5)  ! hour
   DATEA(5) = VALUES(6)  ! minute
   DATEA(6) = VALUES(7)  ! second
END SUBROUTINE DTNOW

! ============================================================================
! LOWERC — convert a character string to lowercase in place.
! (mirrors TTutil lowerc.for)
! ============================================================================
SUBROUTINE LOWERC(STRING)
   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(INOUT) :: STRING
   INTEGER :: I, IC, L
   L = LEN(STRING)
   DO I = 1, L
      IC = ICHAR(STRING(I:I))
      IF (IC >= 65 .AND. IC <= 90) STRING(I:I) = CHAR(IC + 32)
   END DO
END SUBROUTINE LOWERC

! ============================================================================
! UPPERC — convert a character string to uppercase in place.
! (mirrors TTutil upperc.f90)
! ============================================================================
SUBROUTINE UPPERC(STRING)
   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(INOUT) :: STRING
   INTEGER :: I, IC, L
   L = LEN_TRIM(STRING)
   DO I = 1, L
      IC = ICHAR(STRING(I:I))
      IF (IC >= 97 .AND. IC <= 122) STRING(I:I) = CHAR(IC - 32)
   END DO
END SUBROUTINE UPPERC

! ============================================================================
! ADDSTR — append TMP (stripped of leading/trailing spaces) to STRING,
!!         advancing SIGLEN to track the current significant length.
!! (mirrors TTutil addstr.for; ISTART inlined)
! ============================================================================
SUBROUTINE ADDSTR(STRING, SIGLEN, TMP)
   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(INOUT) :: STRING
   INTEGER,          INTENT(INOUT) :: SIGLEN
   CHARACTER(LEN=*), INTENT(IN)    :: TMP

   INTEGER :: IS, IL, L

   ! Find first non-blank character (ISTART)
   IS = 0
   DO L = 1, LEN(TMP)
      IF (TMP(L:L) /= ' ') THEN
         IS = L
         EXIT
      END IF
   END DO
   IL = LEN_TRIM(TMP)

   IF (IL > 0 .AND. IS > 0) THEN
      L = IL - IS + 1
      IF (SIGLEN + L > LEN(STRING)) THEN
         error stop 'ADDSTR: string buffer overflow'
      END IF
      STRING(SIGLEN+1:SIGLEN+L) = TMP(IS:IL)
      SIGLEN = SIGLEN + L
   END IF
END SUBROUTINE ADDSTR

! ============================================================================
! WORDS — split RECORD into words delimited by characters in SEPARS.
!!        Returns up to ILW words; IWBEG(i)/IWEND(i) are the begin/end
!!        positions of word i in RECORD; IFND is the number found.
!! (mirrors TTutil words.for)
! ============================================================================
SUBROUTINE WORDS(RECORD, ILW, SEPARS, IWBEG, IWEND, IFND)
   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(IN)  :: RECORD, SEPARS
   INTEGER,          INTENT(IN)  :: ILW
   INTEGER,          INTENT(OUT) :: IWBEG(ILW), IWEND(ILW), IFND

   INTEGER :: I, ICHR
   CHARACTER(LEN=1) :: CH
   LOGICAL :: SEPAR, GOING

   ICHR = LEN(RECORD)
   IFND = 0
   IWBEG = 0
   IWEND = 0
   GOING = .FALSE.

   DO I = 1, ICHR
      CH    = RECORD(I:I)
      SEPAR = (INDEX(SEPARS, CH) /= 0)

      IF (.NOT. SEPAR .AND. .NOT. GOING) THEN
         IFND = IFND + 1
         IWBEG(IFND) = I
         GOING = .TRUE.
      END IF

      IF (SEPAR .AND. GOING) THEN
         IWEND(IFND) = I - 1
         IF (IFND == ILW) RETURN
         GOING = .FALSE.
      END IF
   END DO

   IF (GOING) IWEND(IFND) = ICHR
END SUBROUTINE WORDS

! ============================================================================
! DECREA — parse a character string as a real number.
!!         IWAR=0 on success; IWAR=1 if STRING is not a valid number.
!!         VALUE is set to the parsed value, or 0.0 on failure.
!! (mirrors TTutil decrea.for; uses native Fortran internal READ)
! ============================================================================
SUBROUTINE DECREA(IWAR, STRING, VALUE)
   IMPLICIT NONE
   INTEGER,          INTENT(OUT) :: IWAR
   CHARACTER(LEN=*), INTENT(IN)  :: STRING
   REAL,             INTENT(OUT) :: VALUE

   INTEGER :: IOS
   CHARACTER(LEN=LEN(STRING)) :: STRTMP

   STRTMP = ADJUSTL(STRING)
   READ(STRTMP, *, IOSTAT=IOS) VALUE
   IF (IOS == 0) THEN
      IWAR = 0
   ELSE
      IWAR  = 1
      VALUE = 0.0
   END IF
END SUBROUTINE DECREA

! ============================================================================
! IFINDI — search integer array ILIS(1:ILDEC) for value IINP between
!!          indices IST and IEND (inclusive, signed step).
!!          Returns the index of the first match, or 0 if not found.
!! (mirrors TTutil ifindi.for)
! ============================================================================
INTEGER FUNCTION IFINDI(ILIS, ILDEC, IST, IEND, IINP)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: ILDEC, IST, IEND, IINP
   INTEGER, INTENT(IN) :: ILIS(ILDEC)

   INTEGER :: IN, STEP

   IF (IEND == 0) THEN
      IFINDI = 0
      RETURN
   END IF

   STEP = 1
   IF (IEND < IST) STEP = -1

   IFINDI = 0
   DO IN = IST, IEND, STEP
      IF (IINP == ILIS(IN)) THEN
         IFINDI = IN
         RETURN
      END IF
   END DO
END FUNCTION IFINDI
