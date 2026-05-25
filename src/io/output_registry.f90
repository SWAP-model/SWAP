!> Declarative catalogue of SWAP output variables. The module-level registry
!! is the single source of truth for the value-array layout (registry index ==
!! value index) and the emitted CSV header/units. Ported verbatim from the
!! former inline `data` block in swap_csv_output.f90:94-198.
module output_registry_mod
   use error_mod, only: error_collection_t, ERR_VALIDATION_ENUM
   implicit none
   private

   integer, parameter, public :: OUT_SCALAR    = 1
   integer, parameter, public :: OUT_NODE      = 2
   integer, parameter, public :: OUT_SUBREGION = 3

   type, public :: out_var_t
      character(len=32) :: name = ''
      character(len=16) :: unit = ''
      integer           :: kind = OUT_SCALAR
   end type out_var_t

   public :: var_count, var_name, var_unit, var_kind, resolve_inlist

   ! 103 entries — verbatim port of swap_csv_output.f90:94-198.
   ! Scalars (indices 1-77): RAIN … QINFMAX
   ! Node templates (indices 78-86): H[ … HEACON[
   ! Scalars (indices 87-88): TETOP, TEBOT
   ! Node templates (indices 89-92): DRAIN[ … SSDI[
   ! Subregion templates (indices 93-103): WTOT[ … QDRAINOUT[
   type(out_var_t), parameter :: REGISTRY(103) = [ &
      ! -- Scalars 1-29 --
      out_var_t('RAIN',       '(cm)',       OUT_SCALAR), &
      out_var_t('RAIN_NET',   '(cm)',       OUT_SCALAR), &
      out_var_t('SNOW',       '(cm)',       OUT_SCALAR), &
      out_var_t('IRRIG',      '(cm)',       OUT_SCALAR), &
      out_var_t('IRRIG_NET',  '(cm)',       OUT_SCALAR), &
      out_var_t('INTERC',     '(cm)',       OUT_SCALAR), &
      out_var_t('RUNON',      '(cm)',       OUT_SCALAR), &
      out_var_t('RUNOFF',     '(cm)',       OUT_SCALAR), &
      out_var_t('EPOT',       '(cm)',       OUT_SCALAR), &
      out_var_t('EACT',       '(cm)',       OUT_SCALAR), &
      out_var_t('SUBLIM',     '(cm)',       OUT_SCALAR), &
      out_var_t('DRAINAGE',   '(cm)',       OUT_SCALAR), &
      out_var_t('QBOTTOM',    '(cm)',       OUT_SCALAR), &
      out_var_t('GWL',        '(cm)',       OUT_SCALAR), &
      out_var_t('POND',       '(cm)',       OUT_SCALAR), &
      out_var_t('SSNOW',      '(cm)',       OUT_SCALAR), &
      out_var_t('TPOT',       '(cm)',       OUT_SCALAR), &
      out_var_t('TACT',       '(cm)',       OUT_SCALAR), &
      out_var_t('TREDDRY',    '(cm)',       OUT_SCALAR), &
      out_var_t('TREDWET',    '(cm)',       OUT_SCALAR), &
      out_var_t('TREDSOL',    '(cm)',       OUT_SCALAR), &
      out_var_t('TREDFRS',    '(cm)',       OUT_SCALAR), &
      out_var_t('ES0',        '(cm)',       OUT_SCALAR), &
      out_var_t('ET0',        '(cm)',       OUT_SCALAR), &
      out_var_t('EW0',        '(cm)',       OUT_SCALAR), &
      out_var_t('DSTOR',      '(cm)',       OUT_SCALAR), &
      out_var_t('BALDEV',     '(cm)',       OUT_SCALAR), &
      out_var_t('VOLACT',     '(cm)',       OUT_SCALAR), &
      out_var_t('QSSDI',      '(cm)',       OUT_SCALAR), &
      ! -- Scalars 30-50 --
      out_var_t('TSUM',       '(deg C)',    OUT_SCALAR), &
      out_var_t('DVS',        '(-)',        OUT_SCALAR), &
      out_var_t('PGASSPOT',   '(kgch/ha)', OUT_SCALAR), &
      out_var_t('PGASS',      '(kgch/ha)', OUT_SCALAR), &
      out_var_t('CPWDM',      '(kg/ha)',   OUT_SCALAR), &
      out_var_t('CWDM',       '(kg/ha)',   OUT_SCALAR), &
      out_var_t('CPWSO',      '(kg/ha)',   OUT_SCALAR), &
      out_var_t('CWSO',       '(kg/ha)',   OUT_SCALAR), &
      out_var_t('PWLV',       '(kg/ha)',   OUT_SCALAR), &
      out_var_t('WLV',        '(kg/ha)',   OUT_SCALAR), &
      out_var_t('PWST',       '(kg/ha)',   OUT_SCALAR), &
      out_var_t('WST',        '(kg/ha)',   OUT_SCALAR), &
      out_var_t('PWRT',       '(kg/ha)',   OUT_SCALAR), &
      out_var_t('WRT',        '(kg/ha)',   OUT_SCALAR), &
      out_var_t('DWSO',       '(kg/ha)',   OUT_SCALAR), &
      out_var_t('DWLV',       '(kg/ha)',   OUT_SCALAR), &
      out_var_t('DWLVPOT',    '(kg/ha)',   OUT_SCALAR), &
      out_var_t('DWST',       '(kg/ha)',   OUT_SCALAR), &
      out_var_t('DWSTPOT',    '(kg/ha)',   OUT_SCALAR), &
      out_var_t('DWRT',       '(kg/ha)',   OUT_SCALAR), &
      out_var_t('DWRTPOT',    '(kg/ha)',   OUT_SCALAR), &
      ! -- Scalars 51-77 --
      out_var_t('HEIGHT',     '(cm)',       OUT_SCALAR), &
      out_var_t('CRPFAC',     '(-)',        OUT_SCALAR), &
      out_var_t('LAIPOT',     '(m2/m2)',   OUT_SCALAR), &
      out_var_t('LAI',        '(m2/m2)',   OUT_SCALAR), &
      out_var_t('RDPOT',      '(cm)',       OUT_SCALAR), &
      out_var_t('RD',         '(cm)',       OUT_SCALAR), &
      out_var_t('PGRASSDM',   '(kg/ha)',   OUT_SCALAR), &
      out_var_t('GRASSDM',    '(kg/ha)',   OUT_SCALAR), &
      out_var_t('PMOWDM',     '(kg/ha)',   OUT_SCALAR), &
      out_var_t('MOWDM',      '(kg/ha)',   OUT_SCALAR), &
      out_var_t('PGRAZDM',    '(kg/ha)',   OUT_SCALAR), &
      out_var_t('GRAZDM',     '(kg/ha)',   OUT_SCALAR), &
      out_var_t('PLOSSDM',    '(kg/ha)',   OUT_SCALAR), &
      out_var_t('LOSSDM',     '(kg/ha)',   OUT_SCALAR), &
      out_var_t('SQPREC',     '(g/cm2)',   OUT_SCALAR), &
      out_var_t('SQIRRIG',    '(g/cm2)',   OUT_SCALAR), &
      out_var_t('SQBOT',      '(g/cm2)',   OUT_SCALAR), &
      out_var_t('SQDRA',      '(g/cm2)',   OUT_SCALAR), &
      out_var_t('DECTOT',     '(g/cm2)',   OUT_SCALAR), &
      out_var_t('ROTTOT',     '(g/cm2)',   OUT_SCALAR), &
      out_var_t('SAMPRO',     '(g/cm2)',   OUT_SCALAR), &
      out_var_t('SOLBAL',     '(g/cm2)',   OUT_SCALAR), &
      out_var_t('WC10',       '(cm3/cm3)', OUT_SCALAR), &
      out_var_t('RUNOFFCN',   '(cm) ',     OUT_SCALAR), &  ! trailing space preserved verbatim from legacy data block
      out_var_t('QTOPIN',     '(cm) ',     OUT_SCALAR), &  ! trailing space preserved verbatim from legacy data block
      out_var_t('QTOPOUT',    '(cm) ',     OUT_SCALAR), &  ! trailing space preserved verbatim from legacy data block
      out_var_t('QINFMAX',    '(cm)',       OUT_SCALAR), &
      ! -- Node templates 78-86 --
      out_var_t('H[',         '(cm)',       OUT_NODE), &
      out_var_t('WC[',        '(cm3/cm3)', OUT_NODE), &
      out_var_t('TEMP[',      '(deg C)',    OUT_NODE), &
      out_var_t('K[',         '(cm/d)',     OUT_NODE), &
      out_var_t('CONC[',      '(g/cm3 w)', OUT_NODE), &
      out_var_t('CONCADS[',   '(g/cm3)',   OUT_NODE), &
      out_var_t('O2TOP[',     '(kg/m3)',   OUT_NODE), &
      out_var_t('HEACAP[',    '(J/cm3/K)', OUT_NODE), &
      out_var_t('HEACON[',    '(J/cm/K/d)',OUT_NODE), &
      ! -- Scalars 87-88 --
      out_var_t('TETOP',      '(deg C)',    OUT_SCALAR), &
      out_var_t('TEBOT',      '(deg C)',    OUT_SCALAR), &
      ! -- Node templates 89-92 --
      out_var_t('DRAIN[',     '(cm)',       OUT_NODE), &
      out_var_t('RWU[',       '(cm)',       OUT_NODE), &
      out_var_t('FLUX[',      '(cm)',       OUT_NODE), &
      out_var_t('SSDI[',      '(cm)',       OUT_NODE), &
      ! -- Subregion templates 93-103 --
      out_var_t('WTOT[',      '(cm)',       OUT_SUBREGION), &
      out_var_t('QTRANS[',    '(cm)',       OUT_SUBREGION), &
      out_var_t('QTOP[',      '(cm)',       OUT_SUBREGION), &
      out_var_t('QBOT[',      '(cm)',       OUT_SUBREGION), &
      out_var_t('QDRA[',      '(cm)',       OUT_SUBREGION), &
      out_var_t('QTOPIN[',    '(cm)',       OUT_SUBREGION), &
      out_var_t('QTOPOUT[',   '(cm)',       OUT_SUBREGION), &
      out_var_t('QBOTIN[',    '(cm)',       OUT_SUBREGION), &
      out_var_t('QBOTOUT[',   '(cm)',       OUT_SUBREGION), &
      out_var_t('QDRAININ[',  '(cm)',       OUT_SUBREGION), &
      out_var_t('QDRAINOUT[', '(cm)',       OUT_SUBREGION)  &
      ]

contains

   pure integer function var_count()
      var_count = size(REGISTRY)
   end function var_count

   pure function var_name(i) result(s)
      integer, intent(in) :: i
      character(len=32)   :: s
      s = REGISTRY(i)%name
   end function var_name

   pure function var_unit(i) result(s)
      integer, intent(in) :: i
      character(len=16)   :: s
      s = REGISTRY(i)%unit
   end function var_unit

   pure integer function var_kind(i)
      integer, intent(in) :: i
      var_kind = REGISTRY(i)%kind
   end function var_kind

   !> Parse a comma-separated inlist, match each token against the registry,
   !! return matched indices in REGISTRY order. Bracket-aware: a token
   !! containing '[' matches a templated entry whose stored name is
   !! "<base>[" (e.g. 'QTOPIN[0:-30]' -> registry name 'QTOPIN['); a token
   !! with NO '[' matches a scalar entry "<base>" (e.g. 'QTOPIN' -> 'QTOPIN').
   !! This disambiguates dual names like QTOPIN (scalar) vs QTOPIN[ (subregion).
   !! Unknown tokens append a fatal error.
   subroutine resolve_inlist(inlist, sel, errors)
      character(len=*),         intent(in)    :: inlist
      integer, allocatable,     intent(out)   :: sel(:)
      type(error_collection_t), intent(inout) :: errors
      logical :: want(size(REGISTRY))
      character(len=:), allocatable :: rest, tok, key
      integer :: comma, br, i, n
      want = .false.
      rest = adjustl(inlist)
      do while (len_trim(rest) > 0)
         comma = index(rest, ',')
         if (comma > 0) then
            tok = trim(adjustl(rest(1:comma-1))); rest = adjustl(rest(comma+1:))
         else
            tok = trim(adjustl(rest)); rest = ''
         end if
         if (len_trim(tok) == 0) cycle
         br = index(tok, '[')
         if (br > 0) then
            key = tok(1:br-1) // '['        ! templated -> match "<base>["
         else
            key = tok                        ! scalar -> match "<base>"
         end if
         n = match(key)
         if (n == 0) then
            call errors%append(ERR_VALIDATION_ENUM, &
               "unknown output variable: '" // trim(tok) // "'", 'output_registry')
         else
            want(n) = .true.
         end if
      end do
      sel = pack([(i, i=1,size(REGISTRY))], want)
   end subroutine resolve_inlist

   pure integer function match(key)
      character(len=*), intent(in) :: key
      character(len=len(key))      :: ukey
      integer :: i
      ukey = to_upper(key)
      match = 0
      do i = 1, size(REGISTRY)
         if (trim(REGISTRY(i)%name) == trim(ukey)) then
            match = i; return
         end if
      end do
   end function match

   pure function to_upper(s) result(out)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: out
      integer :: i, c
      do i = 1, len(s)
         c = iachar(s(i:i))
         if (c >= iachar('a') .and. c <= iachar('z')) c = c - 32
         out(i:i) = achar(c)
      end do
   end function to_upper

end module output_registry_mod
