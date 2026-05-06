!> Reader for the [soil] section of a SWAP TOML.
!!
!! Phase 4a scope: top-level scalars + initial conditions. Per-layer
!! arrays (sublay, hcomp, orgmat, etc.) and the 2D cofgen array are
!! loader-skeleton only; the hupselbrook parity test drives completion
!! of per-layer paths.
!!
!! Phase 4f-prep extension: [soil.discretization], [soil.frost],
!! top-level cofani(:) and nrstaring.
module read_soil_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, get_value, len
   use soil_config_mod, only: soil_config_t
   use toml_field_helpers_mod, only: get_table, &
                                     get_optional_int_with_default, &
                                     get_optional_real_with_default, &
                                     get_optional_string_with_default
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   use read_soil_tillage_toml_mod, only: read_soil_tillage_toml
   implicit none
   private

   public :: read_soil_toml

contains

   subroutine read_soil_toml(doc, config, errors)
      type(toml_table), pointer, intent(in)    :: doc
      type(soil_config_t),       intent(inout) :: config
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), pointer :: sec, initial, discr, frost

      call get_table(doc, 'soil', sec, 'soil', errors)
      if (.not. associated(sec)) return

      call get_optional_int_with_default(sec, 'swsophy', config%swsophy, 0, 'soil.swsophy', errors)
      call get_optional_int_with_default(sec, 'swhyst',  config%swhyst,  0, 'soil.swhyst',  errors)
      call get_optional_int_with_default(sec, 'swinco',  config%swinco,  1, 'soil.swinco',  errors)
      call get_optional_int_with_default(sec, 'swmacro', config%swmacro, 0, 'soil.swmacro', errors)
      call get_optional_int_with_default(sec, 'swscal',  config%swscal,  0, 'soil.swscal',  errors)

      call get_optional_real_with_default(sec, 'ksatexm', config%ksatexm, 0.0_real64, 'soil.ksatexm', errors)
      call get_optional_real_with_default(sec, 'rsoil',   config%rsoil,   0.0_real64, 'soil.rsoil',   errors)
      call get_optional_real_with_default(sec, 'rsro',    config%rsro,    0.0_real64, 'soil.rsro',    errors)
      call get_optional_real_with_default(sec, 'rsroexp', config%rsroexp, 0.0_real64, 'soil.rsroexp', errors)
      call get_optional_int_with_default(sec, 'swrunon',  config%swrunon, 0,          'soil.swrunon', errors)
      call get_optional_int_with_default(sec, 'reva_top', config%reva_top, 0, 'soil.reva_top', errors)

      call get_optional_int_with_default(sec, 'nrstaring', config%nrstaring, 0, 'soil.nrstaring', errors)

      ! Top-level per-layer anisotropy ratios.
      call read_array_1d(sec, 'cofani', config%cofani, 'soil.cofani', errors)

      call get_table(sec, 'initial', initial, 'soil.initial', errors)
      if (associated(initial)) then
         call get_optional_real_with_default(initial, 'gwli',    config%gwli,    0.0_real64, 'soil.initial.gwli',    errors)
         call get_optional_real_with_default(initial, 'pondini', config%pondini, 0.0_real64, 'soil.initial.pondini', errors)
         call get_optional_real_with_default(initial, 'pondmx',  config%pondmx,  0.0_real64, 'soil.initial.pondmx',  errors)

         ! [soil.initial] sub-table — typed home for what swap.ini used to carry.
         call get_optional_int_with_default(initial,    'swirrigate', config%initial%swirrigate, 0, &
                                            'soil.initial.swirrigate', errors)
         call get_optional_real_with_default(initial,   'ssnow',  config%initial%ssnow,  0.0_real64, &
                                             'soil.initial.ssnow',  errors)
         call get_optional_real_with_default(initial,   'slw',    config%initial%slw,    0.0_real64, &
                                             'soil.initial.slw',    errors)
         call get_optional_real_with_default(initial,   'pond',   config%initial%pond,   0.0_real64, &
                                             'soil.initial.pond',   errors)
         call get_optional_real_with_default(initial,   'ldwet',  config%initial%ldwet,  0.0_real64, &
                                             'soil.initial.ldwet',  errors)
         call get_optional_real_with_default(initial,   'dt',     config%initial%dt,     0.0_real64, &
                                             'soil.initial.dt',     errors)
         call get_optional_string_with_default(initial, 'h_file',     config%initial%h_file,     '', &
                                               'soil.initial.h_file',     errors)
         call get_optional_string_with_default(initial, 'tsoil_file', config%initial%tsoil_file, '', &
                                               'soil.initial.tsoil_file', errors)
         call get_optional_string_with_default(initial, 'cml_file',   config%initial%cml_file,   '', &
                                               'soil.initial.cml_file',   errors)

         ! atmin7: fixed-length 7-element inline TOML array.
         block
            type(toml_array), pointer :: arr_ptr
            integer :: ilen, k, stat
            real(real64) :: tmp
            arr_ptr => null()
            call get_value(initial, 'atmin7', arr_ptr, requested=.false., stat=stat)
            if (associated(arr_ptr)) then
               ilen = min(len(arr_ptr), 7)
               do k = 1, ilen
                  call get_value(arr_ptr, k, tmp, stat=stat)
                  if (stat == 0) config%initial%atmin7(k) = tmp
               end do
            end if
         end block
      end if

      call get_table(sec, 'discretization', discr, 'soil.discretization', errors)
      if (associated(discr)) then
         call get_optional_int_with_default(discr, 'swdiscrvert', config%discretization%swdiscrvert, 0, &
                                            'soil.discretization.swdiscrvert', errors)
         call get_optional_int_with_default(discr, 'numnodnew',   config%discretization%numnodnew,   0, &
                                            'soil.discretization.numnodnew',   errors)
         call read_array_1d(discr, 'dznew', config%discretization%dznew, &
                            'soil.discretization.dznew', errors)
      end if

      call get_table(sec, 'frost', frost, 'soil.frost', errors)
      if (associated(frost)) then
         call get_optional_int_with_default(frost,  'swfrost',   config%frost%swfrost,   0,          'soil.frost.swfrost',   errors)
         call get_optional_real_with_default(frost, 'tfroststa', config%frost%tfroststa, 0.0_real64, 'soil.frost.tfroststa', errors)
         call get_optional_real_with_default(frost, 'tfrostend', config%frost%tfrostend, 0.0_real64, 'soil.frost.tfrostend', errors)
         call get_optional_int_with_default(frost,  'swsublim',  config%frost%swsublim,  0,          'soil.frost.swsublim',  errors)
      end if

      ! Vertical sub-layer discretization. Phase 4f Task B2: per-sub-layer
      ! arrays (sublay, isoillay, hsublay, ncomp) authored at top level
      ! of [soil]; calcgrid() consumes these to build numnod / dz / z /
      ! disnod / layer.
      call read_int_array_1d (sec, 'sublay',   config%sublay,   'soil.sublay',   errors)
      call read_int_array_1d (sec, 'isoillay', config%isoillay, 'soil.isoillay', errors)
      call read_array_1d     (sec, 'hsublay',  config%hsublay,  'soil.hsublay',  errors)
      call read_int_array_1d (sec, 'ncomp',    config%ncomp,    'soil.ncomp',    errors)
      ! hcomp is derived as hsublay/ncomp by the adapter; do not author
      ! it directly. (See readswap.f90:613-619 for the legacy derivation.)

      ! Per-soil-physical-layer Mualem-van Genuchten hydraulics. Authored
      ! as parallel 1D arrays (one cell per soil-physical layer) inside
      ! [soil.hydraulics]. The adapter copies these into the legacy
      ! variables%ores / variables%osat / ... arrays and then builds
      ! variables%paramvg(1..10, lay) just like readswap.f90:786-825.
      call read_hydraulics(sec, config, errors)

      ! Optional [soil.tillage] block — parsed by sibling module to keep
      ! this file focused. Block is absent when swtill = 0.
      call read_soil_tillage_toml(sec, config%tillage, errors)
   end subroutine read_soil_toml

   !> Read [soil.hydraulics] into the typed sub-config. Each key is an
   !! optional 1-D real array; missing arrays are left unallocated and
   !! the validator does not enforce presence (the runtime aborts later
   !! if needed). Section absence is silent.
   subroutine read_hydraulics(soil_sec, config, errors)
      type(toml_table), pointer, intent(in)    :: soil_sec
      type(soil_config_t),       intent(inout) :: config
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), pointer :: hyd

      call get_table(soil_sec, 'hydraulics', hyd, 'soil.hydraulics', errors)
      if (.not. associated(hyd)) return

      call read_array_1d(hyd, 'ores',    config%hydraulics%ores,    'soil.hydraulics.ores',    errors)
      call read_array_1d(hyd, 'osat',    config%hydraulics%osat,    'soil.hydraulics.osat',    errors)
      call read_array_1d(hyd, 'alfa',    config%hydraulics%alfa,    'soil.hydraulics.alfa',    errors)
      call read_array_1d(hyd, 'npar',    config%hydraulics%npar,    'soil.hydraulics.npar',    errors)
      call read_array_1d(hyd, 'ksatfit', config%hydraulics%ksatfit, 'soil.hydraulics.ksatfit', errors)
      call read_array_1d(hyd, 'lexp',    config%hydraulics%lexp,    'soil.hydraulics.lexp',    errors)
      call read_array_1d(hyd, 'alfaw',   config%hydraulics%alfaw,   'soil.hydraulics.alfaw',   errors)
      call read_array_1d(hyd, 'h_enpr',  config%hydraulics%h_enpr,  'soil.hydraulics.h_enpr',  errors)
      call read_array_1d(hyd, 'ksatexm', config%hydraulics%ksatexm, 'soil.hydraulics.ksatexm', errors)
      call read_array_1d(hyd, 'bdens',   config%hydraulics%bdens,   'soil.hydraulics.bdens',   errors)
   end subroutine read_hydraulics

   !> Decode a flat TOML int array at sec[key] into a 1-D integer
   !! allocatable. Mirrors `read_array_1d` semantics for ints.
   subroutine read_int_array_1d(sec, key, arr, context, errors)
      type(toml_table), pointer, intent(in)    :: sec
      character(len=*),          intent(in)    :: key
      integer, allocatable,      intent(out)   :: arr(:)
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors

      type(toml_array), pointer :: outer
      integer :: n, i, stat, val

      if (.not. associated(sec)) return

      outer => null()
      call get_value(sec, key, outer, requested=.false., stat=stat)
      if (.not. associated(outer)) return

      n = len(outer)
      if (n == 0) then
         allocate(arr(0))
         return
      end if

      allocate(arr(n))
      arr = 0

      do i = 1, n
         call get_value(outer, i, val, stat=stat)
         if (stat /= 0) then
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               "non-int cell", context)
            if (allocated(arr)) deallocate(arr)
            return
         end if
         arr(i) = val
      end do
   end subroutine read_int_array_1d

   !> Decode a flat TOML real array at sec[key] into a 1-D real(real64)
   !! allocatable. Absent key leaves arr unallocated. Empty array yields
   !! a 0-element allocation. Non-real cells append a parse-type-mismatch
   !! error and leave arr unallocated. Mirrors read_heat_toml's local helper.
   subroutine read_array_1d(sec, key, arr, context, errors)
      type(toml_table), pointer, intent(in)    :: sec
      character(len=*),          intent(in)    :: key
      real(real64), allocatable, intent(out)   :: arr(:)
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors

      type(toml_array), pointer :: outer
      integer :: n, i, stat
      real(real64) :: val

      if (.not. associated(sec)) return

      outer => null()
      call get_value(sec, key, outer, requested=.false., stat=stat)
      if (.not. associated(outer)) return

      n = len(outer)
      if (n == 0) then
         allocate(arr(0))
         return
      end if

      allocate(arr(n))
      arr = 0.0_real64

      do i = 1, n
         call get_value(outer, i, val, stat=stat)
         if (stat /= 0) then
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               "non-real cell", context)
            if (allocated(arr)) deallocate(arr)
            return
         end if
         arr(i) = val
      end do
   end subroutine read_array_1d

end module read_soil_toml_mod
