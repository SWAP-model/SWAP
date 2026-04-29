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
                                     get_optional_real_with_default
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
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
      call get_optional_int_with_default(sec, 'reva_top', config%reva_top, 0, 'soil.reva_top', errors)

      call get_optional_int_with_default(sec, 'nrstaring', config%nrstaring, 0, 'soil.nrstaring', errors)

      ! Top-level per-layer anisotropy ratios.
      call read_array_1d(sec, 'cofani', config%cofani, 'soil.cofani', errors)

      call get_table(sec, 'initial', initial, 'soil.initial', errors)
      if (associated(initial)) then
         call get_optional_real_with_default(initial, 'gwli',    config%gwli,    0.0_real64, 'soil.initial.gwli',    errors)
         call get_optional_real_with_default(initial, 'pondini', config%pondini, 0.0_real64, 'soil.initial.pondini', errors)
         call get_optional_real_with_default(initial, 'pondmx',  config%pondmx,  0.0_real64, 'soil.initial.pondmx',  errors)
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
   end subroutine read_soil_toml

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
