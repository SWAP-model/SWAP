!> @file read_nutrients_toml.f90
!! Parses the optional top-level [nutrients] block into
!! nutrients_config_t. Block is optional; absent → leaves
!! defaults intact and sets present = .false..
!!
!! See ADR 0026 ([nutrients] N2a).
module read_nutrients_toml_mod

   use, intrinsic :: iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, get_value, len
   use toml_field_helpers_mod, only: get_optional_real_with_default, &
                                     get_optional_string_with_default
   use nutrients_config_mod, only: nutrients_config_t
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_nutrients_toml

contains

   subroutine read_nutrients_toml(doc, cfg, errors)
      type(toml_table), pointer, intent(in)    :: doc
      type(nutrients_config_t),  intent(inout) :: cfg
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), pointer :: nut_tbl, init_tbl
      type(toml_array), pointer :: fom_arr
      integer :: stat, n, i
      real(real64) :: v

      if (.not. associated(doc)) return

      nut_tbl => null()
      call get_value(doc, 'nutrients', nut_tbl, requested=.false., stat=stat)
      if (stat /= 0 .or. .not. associated(nut_tbl)) return

      cfg%present = .true.
      call get_optional_real_with_default(nut_tbl, 'sorp_coef', cfg%sorp_coef, &
                                          0.0_real64, 'nutrients.sorp_coef', errors)
      call get_optional_string_with_default(nut_tbl, 'events_file', &
                                            cfg%events_file, '', &
                                            'nutrients.events_file', errors)

      ! [nutrients.initial] sub-table is optional.
      init_tbl => null()
      call get_value(nut_tbl, 'initial', init_tbl, requested=.false., stat=stat)
      if (stat == 0 .and. associated(init_tbl)) then
         ! fom is a 1- to 8-element array; missing entries stay at 0.0.
         fom_arr => null()
         call get_value(init_tbl, 'fom', fom_arr, requested=.false., stat=stat)
         if (stat == 0 .and. associated(fom_arr)) then
            n = len(fom_arr)
            if (n > 8) then
               call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                  'nutrients.initial.fom: at most 8 entries (legacy FOM_t(1..8) cap)', &
                  'nutrients.initial.fom')
               n = 8
            end if
            do i = 1, n
               call get_value(fom_arr, i, v, stat=stat)
               if (stat == 0) cfg%initial%fom(i) = v
            end do
         end if

         call get_optional_real_with_default(init_tbl, 'bio',  cfg%initial%bio,  &
                                             0.0_real64, 'nutrients.initial.bio',  errors)
         call get_optional_real_with_default(init_tbl, 'hum',  cfg%initial%hum,  &
                                             0.0_real64, 'nutrients.initial.hum',  errors)
         call get_optional_real_with_default(init_tbl, 'cnh4', cfg%initial%cnh4, &
                                             0.0_real64, 'nutrients.initial.cnh4', errors)
         call get_optional_real_with_default(init_tbl, 'cno3', cfg%initial%cno3, &
                                             0.0_real64, 'nutrients.initial.cno3', errors)
      end if
   end subroutine read_nutrients_toml

end module read_nutrients_toml_mod
