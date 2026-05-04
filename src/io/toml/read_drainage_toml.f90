!> TOML reader for the [drainage] section.
!! Supports inline drainage data OR external file via `[drainage].file = "..."`.
module read_drainage_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, toml_error, toml_load, get_value, len
   use drainage_config_mod, only: drainage_config_t
   use toml_field_helpers_mod, only: get_table, get_array_of_tables,   &
                                     get_optional_int_with_default,    &
                                     get_optional_real_with_default,   &
                                     get_optional_string_with_default
   use path_helpers_mod, only: resolve_relative_path
   use error_mod, only: error_collection_t, ERR_PARSE_MALFORMED_TOML, &
                        ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_drainage_toml

contains

   !> Decode `drainage.cofani = [...]` into a 1-D real array sized by
   !! the input length. Absent key leaves arr unallocated. Mirrors the
   !! local read_array_1d helpers in read_soil_toml / read_heat_toml.
   subroutine read_cofani(sec, arr, errors)
      type(toml_table), pointer, intent(in)    :: sec
      real(real64), allocatable, intent(out)   :: arr(:)
      type(error_collection_t),  intent(inout) :: errors

      type(toml_array), pointer :: outer
      integer :: n, i, stat
      real(real64) :: val

      if (.not. associated(sec)) return

      outer => null()
      call get_value(sec, 'cofani', outer, requested=.false., stat=stat)
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
                               "non-real cell", "drainage.cofani")
            if (allocated(arr)) deallocate(arr)
            return
         end if
         arr(i) = val
      end do
   end subroutine read_cofani

   subroutine read_drainage_surface_runoff(sec, config, errors)
      type(toml_table), pointer,  intent(in)    :: sec
      type(drainage_config_t),    intent(inout) :: config
      type(error_collection_t),   intent(inout) :: errors

      type(toml_table), pointer :: sr

      call get_table(sec, 'surface_runoff', sr, 'drainage.surface_runoff', errors)
      if (.not. associated(sr)) return

      call get_optional_int_with_default(sr, 'swnrsrf',      config%surface_runoff%swnrsrf,      0, &
                                         'drainage.surface_runoff.swnrsrf',      errors)
      call get_optional_int_with_default(sr, 'swtopnrsrf',   config%surface_runoff%swtopnrsrf,   0, &
                                         'drainage.surface_runoff.swtopnrsrf',   errors)
      call get_optional_int_with_default(sr, 'swdivdinf',    config%surface_runoff%swdivdinf,    0, &
                                         'drainage.surface_runoff.swdivdinf',    errors)
      call get_optional_int_with_default(sr, 'swtopdislay',  config%surface_runoff%swtopdislay,  0, &
                                         'drainage.surface_runoff.swtopdislay',  errors)
      call get_optional_int_with_default(sr, 'numlevrapdra', config%surface_runoff%numlevrapdra, 0, &
                                         'drainage.surface_runoff.numlevrapdra', errors)

      call get_optional_real_with_default(sr, 'facdpthinf',   config%surface_runoff%facdpthinf,   0.0_real64, &
                                          'drainage.surface_runoff.facdpthinf',   errors)
      call get_optional_real_with_default(sr, 'cofintfl',     config%surface_runoff%cofintfl,     0.0_real64, &
                                          'drainage.surface_runoff.cofintfl',     errors)
      call get_optional_real_with_default(sr, 'expintfl',     config%surface_runoff%expintfl,     0.0_real64, &
                                          'drainage.surface_runoff.expintfl',     errors)
      call get_optional_real_with_default(sr, 'geofac',       config%surface_runoff%geofac,       0.0_real64, &
                                          'drainage.surface_runoff.geofac',       errors)
      call get_optional_real_with_default(sr, 'gwlconv',      config%surface_runoff%gwlconv,      0.0_real64, &
                                          'drainage.surface_runoff.gwlconv',      errors)
      call get_optional_real_with_default(sr, 'ftopdislay',   config%surface_runoff%ftopdislay,   0.0_real64, &
                                          'drainage.surface_runoff.ftopdislay',   errors)
      call get_optional_real_with_default(sr, 'rsurfdeep',    config%surface_runoff%rsurfdeep,    0.0_real64, &
                                          'drainage.surface_runoff.rsurfdeep',    errors)
      call get_optional_real_with_default(sr, 'rsurfshallow', config%surface_runoff%rsurfshallow, 0.0_real64, &
                                          'drainage.surface_runoff.rsurfshallow', errors)
      call get_optional_real_with_default(sr, 'rapdrareaexp', config%surface_runoff%rapdrareaexp, 0.0_real64, &
                                          'drainage.surface_runoff.rapdrareaexp', errors)
      call get_optional_real_with_default(sr, 'rapdraresref', config%surface_runoff%rapdraresref, 0.0_real64, &
                                          'drainage.surface_runoff.rapdraresref', errors)
   end subroutine read_drainage_surface_runoff

   subroutine read_drainage_toml(doc, config, errors, base_path)
      type(toml_table), pointer,  intent(in)    :: doc
      type(drainage_config_t),    intent(inout) :: config
      type(error_collection_t),   intent(inout) :: errors
      character(len=*), optional, intent(in)    :: base_path

      type(toml_table), pointer             :: drain_tab, ext_root, ext_sec
      type(toml_table), allocatable, target :: ext_doc
      type(toml_error), allocatable         :: terr
      character(len=:), allocatable         :: file_rel, file_abs

      call get_table(doc, 'drainage', drain_tab, 'drainage', errors)
      if (.not. associated(drain_tab)) return

      call get_optional_string_with_default(drain_tab, 'file', file_rel, '', 'drainage.file', errors)
      if (len_trim(file_rel) > 0 .and. present(base_path)) then
         file_abs = resolve_relative_path(base_path, file_rel)
         call toml_load(ext_doc, trim(file_abs), error=terr)
         if (allocated(terr)) then
            call errors%append(ERR_PARSE_MALFORMED_TOML, trim(terr%message), trim(file_abs))
            return
         end if
         ext_root => ext_doc
         call get_table(ext_root, 'drainage', ext_sec, 'drainage', errors)
         if (.not. associated(ext_sec)) return
         call read_drainage_inner(ext_sec, config, errors)
      else
         call read_drainage_inner(drain_tab, config, errors)
      end if
   end subroutine read_drainage_toml

   subroutine read_drainage_inner(sec, config, errors)
      type(toml_table), pointer,  intent(in)    :: sec
      type(drainage_config_t),    intent(inout) :: config
      type(error_collection_t),   intent(inout) :: errors

      type(toml_table), pointer :: basic, item
      type(toml_array), pointer :: levels
      integer :: i, n, stat
      character(len=:), allocatable :: str_tmp

      call get_optional_int_with_default(sec, 'swdra',    config%swdra,    0, 'drainage.swdra',    errors)
      call get_optional_int_with_default(sec, 'dramet',   config%dramet,   0, 'drainage.dramet',   errors)
      call get_optional_string_with_default(sec, 'drfil', config%drfil, 'swap', 'drainage.drfil', errors)
      call get_optional_int_with_default(sec, 'swdivd',   config%swdivd,   0, 'drainage.swdivd',   errors)
      call get_optional_int_with_default(sec, 'swdislay', config%swdislay, 0, 'drainage.swdislay', errors)
      call get_optional_int_with_default(sec, 'nrlevs',   config%nrlevs,   0, 'drainage.nrlevs',   errors)
      call get_optional_int_with_default(sec, 'swliminf', config%swliminf, 0, 'drainage.swliminf', errors)
      call get_optional_real_with_default(sec, 'altcu',   config%altcu,    0.0_real64, 'drainage.altcu', errors)

      call get_table(sec, 'basic', basic, 'drainage.basic', errors)
      if (associated(basic)) then
         call get_optional_real_with_default(basic, 'basegw', config%basegw, 0.0_real64, 'drainage.basic.basegw', errors)
         call get_optional_real_with_default(basic, 'entres', config%entres, 0.0_real64, 'drainage.basic.entres', errors)
         call get_optional_real_with_default(basic, 'shape',  config%shape,  0.0_real64, 'drainage.basic.shape',  errors)
         ! DRAMET=2 (Hooghoudt/Ernst) scalars (legacy .dra Part 2).
         call get_optional_real_with_default(basic, 'lm',     config%lm,           0.0_real64, 'drainage.basic.lm',     errors)
         call get_optional_real_with_default(basic, 'wetper', config%wetper,       0.0_real64, 'drainage.basic.wetper', errors)
         call get_optional_real_with_default(basic, 'zbotdr', config%zbotdr_basic, 0.0_real64, 'drainage.basic.zbotdr', errors)
         call get_optional_int_with_default (basic, 'ipos',   config%ipos,         0,          'drainage.basic.ipos',   errors)
         call get_optional_real_with_default(basic, 'khtop',  config%khtop,        0.0_real64, 'drainage.basic.khtop',  errors)
         call get_optional_real_with_default(basic, 'khbot',  config%khbot,        0.0_real64, 'drainage.basic.khbot',  errors)
         call get_optional_real_with_default(basic, 'kvtop',  config%kvtop,        0.0_real64, 'drainage.basic.kvtop',  errors)
         call get_optional_real_with_default(basic, 'kvbot',  config%kvbot,        0.0_real64, 'drainage.basic.kvbot',  errors)
         call get_optional_real_with_default(basic, 'zintf',  config%zintf,        0.0_real64, 'drainage.basic.zintf',  errors)
         call get_optional_real_with_default(basic, 'geofac', config%geofac,       0.0_real64, 'drainage.basic.geofac', errors)
      end if

      ! Top-level `cofani` array (per soil-physical layer anisotropy).
      call read_cofani(sec, config%cofani, errors)

      call get_array_of_tables(sec, 'levels', levels, 'drainage.levels', errors)
      if (associated(levels)) then
         n = len(levels)
         if (n > 0) then
            allocate(config%swdtyp(n), config%zbotdr(n), config%drares(n), &
                     config%infres(n), config%L(n),      config%gwlinf(n), &
                     config%rdrain(n), config%rinfi(n),  config%rentry(n), &
                     config%rexit(n),  config%widthr(n), config%taludr(n), &
                     config%swallo(n), config%owltab_file(n))
            config%swdtyp      = 0
            config%zbotdr      = 0.0_real64
            config%drares      = 0.0_real64
            config%infres      = 0.0_real64
            config%L           = 0.0_real64
            config%gwlinf      = 0.0_real64
            config%rdrain      = 0.0_real64
            config%rinfi       = 0.0_real64
            config%rentry      = 0.0_real64
            config%rexit       = 0.0_real64
            config%widthr      = 0.0_real64
            config%taludr      = 0.0_real64
            config%swallo      = 0
            config%owltab_file = ''
            do i = 1, n
               call get_value(levels, i, item, stat=stat)
               if (stat /= 0 .or. .not. associated(item)) cycle
               call get_optional_int_with_default(item,  'swdtyp',     config%swdtyp(i),         0,           'drainage.levels.swdtyp',     errors)
               call get_optional_int_with_default(item,  'swallo',     config%swallo(i),         0,           'drainage.levels.swallo',     errors)
               call get_optional_real_with_default(item, 'zbotdr',     config%zbotdr(i),         0.0_real64,  'drainage.levels.zbotdr',     errors)
               call get_optional_real_with_default(item, 'drares',     config%drares(i),         0.0_real64,  'drainage.levels.drares',     errors)
               call get_optional_real_with_default(item, 'infres',     config%infres(i),         0.0_real64,  'drainage.levels.infres',     errors)
               call get_optional_real_with_default(item, 'L',          config%L(i),              0.0_real64,  'drainage.levels.L',          errors)
               call get_optional_real_with_default(item, 'gwlinf',     config%gwlinf(i),         0.0_real64,  'drainage.levels.gwlinf',     errors)
               call get_optional_real_with_default(item, 'rdrain',     config%rdrain(i),         0.0_real64,  'drainage.levels.rdrain',     errors)
               call get_optional_real_with_default(item, 'rinfi',      config%rinfi(i),          0.0_real64,  'drainage.levels.rinfi',      errors)
               call get_optional_real_with_default(item, 'rentry',     config%rentry(i),         0.0_real64,  'drainage.levels.rentry',     errors)
               call get_optional_real_with_default(item, 'rexit',      config%rexit(i),          0.0_real64,  'drainage.levels.rexit',      errors)
               call get_optional_real_with_default(item, 'widthr',     config%widthr(i),         0.0_real64,  'drainage.levels.widthr',     errors)
               call get_optional_real_with_default(item, 'taludr',     config%taludr(i),         0.0_real64,  'drainage.levels.taludr',     errors)
               call get_optional_string_with_default(item, 'owltab_file', str_tmp, '', 'drainage.levels.owltab_file', errors)
               if (allocated(str_tmp)) config%owltab_file(i) = str_tmp
            end do
         end if
      end if

      call read_drainage_surface_runoff(sec, config, errors)
   end subroutine read_drainage_inner

end module read_drainage_toml_mod
