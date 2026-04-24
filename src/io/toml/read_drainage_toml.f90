!> Reader for the [drainage] section of a SWAP TOML.
module read_drainage_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, get_value, len
   use drainage_config_mod, only: drainage_config_t
   use toml_field_helpers_mod, only: get_table, get_array_of_tables,   &
                                     get_optional_int_with_default,    &
                                     get_optional_real_with_default
   use error_mod, only: error_collection_t
   implicit none
   private

   public :: read_drainage_toml

contains

   subroutine read_drainage_toml(doc, config, errors)
      type(toml_table), pointer,  intent(in)    :: doc
      type(drainage_config_t),    intent(inout) :: config
      type(error_collection_t),   intent(inout) :: errors

      type(toml_table), pointer :: sec, basic, item
      type(toml_array), pointer :: levels
      integer :: i, n, stat

      call get_table(doc, 'drainage', sec, 'drainage', errors)
      if (.not. associated(sec)) return

      call get_optional_int_with_default(sec, 'swdra',    config%swdra,    0, 'drainage.swdra',    errors)
      call get_optional_int_with_default(sec, 'dramet',   config%dramet,   0, 'drainage.dramet',   errors)
      call get_optional_int_with_default(sec, 'swdivd',   config%swdivd,   0, 'drainage.swdivd',   errors)
      call get_optional_int_with_default(sec, 'swdislay', config%swdislay, 0, 'drainage.swdislay', errors)
      call get_optional_int_with_default(sec, 'nrlevs',   config%nrlevs,   0, 'drainage.nrlevs',   errors)
      call get_optional_real_with_default(sec, 'altcu',   config%altcu,    0.0_real64, 'drainage.altcu', errors)

      call get_table(sec, 'basic', basic, 'drainage.basic', errors)
      if (associated(basic)) then
         call get_optional_real_with_default(basic, 'basegw', config%basegw, 0.0_real64, 'drainage.basic.basegw', errors)
         call get_optional_real_with_default(basic, 'entres', config%entres, 0.0_real64, 'drainage.basic.entres', errors)
         call get_optional_real_with_default(basic, 'shape',  config%shape,  0.0_real64, 'drainage.basic.shape',  errors)
      end if

      call get_array_of_tables(sec, 'levels', levels, 'drainage.levels', errors)
      if (associated(levels)) then
         n = len(levels)
         if (n > 0) then
            allocate(config%swdtyp(n), config%zbotdr(n), config%drares(n), &
                     config%infres(n), config%L(n),      config%gwlinf(n), &
                     config%rdrain(n), config%rinfi(n),  config%rentry(n), &
                     config%rexit(n),  config%widthr(n), config%taludr(n), &
                     config%swallo(n))
            config%swdtyp  = 0
            config%zbotdr  = 0.0_real64
            config%drares  = 0.0_real64
            config%infres  = 0.0_real64
            config%L       = 0.0_real64
            config%gwlinf  = 0.0_real64
            config%rdrain  = 0.0_real64
            config%rinfi   = 0.0_real64
            config%rentry  = 0.0_real64
            config%rexit   = 0.0_real64
            config%widthr  = 0.0_real64
            config%taludr  = 0.0_real64
            config%swallo  = 0

            do i = 1, n
               call get_value(levels, i, item, stat=stat)
               if (stat /= 0 .or. .not. associated(item)) cycle
               call get_optional_int_with_default(item,  'swdtyp', config%swdtyp(i), 0,           'drainage.levels.swdtyp', errors)
               call get_optional_real_with_default(item, 'zbotdr', config%zbotdr(i), 0.0_real64,  'drainage.levels.zbotdr', errors)
               call get_optional_real_with_default(item, 'drares', config%drares(i), 0.0_real64,  'drainage.levels.drares', errors)
               call get_optional_real_with_default(item, 'infres', config%infres(i), 0.0_real64,  'drainage.levels.infres', errors)
               call get_optional_real_with_default(item, 'L',      config%L(i),      0.0_real64,  'drainage.levels.L',      errors)
            end do
         end if
      end if
   end subroutine read_drainage_toml

end module read_drainage_toml_mod
