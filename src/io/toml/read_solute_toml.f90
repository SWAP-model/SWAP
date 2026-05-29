!> Reader for the [solute] section of swap.toml.
!!
!! Populates a solute_config_t from a fully-formed swap.toml document.
!! Scalars are read one per line; the per-layer decomposition table
!! (`pertabsolu`) is encoded as a TOML array of arrays of two reals
!! (depth, factor) and decoded into a 2-D allocatable real(real64)
!! array.
!!
!! Mechanical pattern (mirrors read_heat_toml / read_bottom_boundary_toml):
!! no switch-conditional skipping — the validator handles cross-field
!! consistency. If [solute] is absent the reader returns silently and
!! leaves `config` at defaults (the validator's `swsolu=0` sentinel
!! short-circuits cleanly in that case).
module read_solute_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, get_value, len
   use solute_config_mod, only: solute_config_t
   use toml_field_helpers_mod, only: get_table,                         &
                                     get_optional_int_with_default,     &
                                     get_optional_real_with_default
   use toml_array_helpers_mod, only: read_table_2d
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_solute_toml

contains

   subroutine read_solute_toml(doc_root, config, errors)
      type(toml_table), pointer, intent(in)    :: doc_root
      type(solute_config_t),     intent(inout) :: config
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), pointer :: sec

      call get_table(doc_root, 'solute', sec, 'solute', errors)
      if (.not. associated(sec)) return

      ! Switches.
      call get_optional_int_with_default(sec, 'swsolu',   config%swsolu,   0, &
                                         'solute.swsolu',   errors)
      call get_optional_int_with_default(sec, 'swbotbc',  config%swbotbc,  0, &
                                         'solute.swbotbc',  errors)
      call get_optional_int_with_default(sec, 'swsoltyp', config%swsoltyp, 0, &
                                         'solute.swsoltyp', errors)
      call get_optional_int_with_default(sec, 'swdc',     config%swdc,     0, &
                                         'solute.swdc',     errors)

      ! Concentrations and dispersion.
      call get_optional_real_with_default(sec, 'cdrain', config%cdrain, 0.0_real64, &
                                          'solute.cdrain', errors)
      call get_optional_real_with_default(sec, 'cseep',  config%cseep,  0.0_real64, &
                                          'solute.cseep',  errors)
      call get_optional_real_with_default(sec, 'tscf',   config%tscf,   0.0_real64, &
                                          'solute.tscf',   errors)
      ! Phase 4f Task B5: ldis can be authored as a scalar OR as a
      ! per-layer array. Try the array form first; if absent fall back
      ! to the scalar (which the adapter broadcasts to all layers).
      call read_real_array(sec, 'ldis', config%ldis_array, errors)
      if (.not. allocated(config%ldis_array)) then
         call get_optional_real_with_default(sec, 'ldis',   config%ldis,   0.0_real64, &
                                             'solute.ldis',   errors)
      end if

      ! Root-uptake.
      call get_optional_real_with_default(sec, 'rtheta', config%rtheta, 0.0_real64, &
                                          'solute.rtheta', errors)
      call get_optional_real_with_default(sec, 'bexp',   config%bexp,   0.0_real64, &
                                          'solute.bexp',   errors)

      ! Salinity scalars.
      call get_optional_real_with_default(sec, 'ecmax',  config%ecmax,  0.0_real64, &
                                          'solute.ecmax',  errors)
      call get_optional_real_with_default(sec, 'ecslop', config%ecslop, 0.0_real64, &
                                          'solute.ecslop', errors)

      ! Per-layer decomposition table.
      call read_table_2d(sec, 'pertabsolu', config%pertabsolu, 2, &
                         'solute.pertabsolu', errors)

      ! Phase 0 (ADR 0032) — physics fields promoted from legacy globals.
      call get_optional_real_with_default(sec, 'cref',   config%cref,   0.0_real64, &
                                          'solute.cref',   errors)
      call get_optional_real_with_default(sec, 'cpre',   config%cpre,   0.0_real64, &
                                          'solute.cpre',   errors)
      call get_optional_real_with_default(sec, 'ddif',   config%ddif,   0.0_real64, &
                                          'solute.ddif',   errors)
      call get_optional_real_with_default(sec, 'frexp',  config%frexp,  0.0_real64, &
                                          'solute.frexp',  errors)
      call get_optional_real_with_default(sec, 'gampar', config%gampar, 0.0_real64, &
                                          'solute.gampar', errors)
      call get_optional_real_with_default(sec, 'daquif', config%daquif, 0.0_real64, &
                                          'solute.daquif', errors)
      call get_optional_real_with_default(sec, 'kfsat',  config%kfsat,  0.0_real64, &
                                          'solute.kfsat',  errors)
      call get_optional_real_with_default(sec, 'decsat', config%decsat, 0.0_real64, &
                                          'solute.decsat', errors)
      call get_optional_real_with_default(sec, 'poros',  config%poros,  0.0_real64, &
                                          'solute.poros',  errors)
      call get_optional_int_with_default (sec, 'swbr',   config%swbr,   0, &
                                          'solute.swbr',   errors)

      ! Per-layer arrays — 1D flat TOML arrays.
      call read_real_array(sec, 'kf',     config%kf,     errors)
      call read_real_array(sec, 'decpot', config%decpot, errors)
      call read_real_array(sec, 'fdepth', config%fdepth, errors)

      ! Seepage concentration table — 2D (rows × 2): col 1 = time, col 2 = concentration.
      ! Adapter flattens to the interleaved afgen layout used by cseeptab(mabbc*2).
      call read_table_2d(sec, 'cseeptab', config%cseeptab, 2, &
                         'solute.cseeptab', errors)
   end subroutine read_solute_toml

   !> Decode a flat TOML real array at sec[key] into a 1-D real(real64)
   !! allocatable. Absent key (or scalar value at the same key) leaves
   !! arr unallocated, letting the caller fall back to a scalar reader.
   !! Distinguished from toml_array_helpers_mod%read_real_array_1d by its
   !! absent context parameter — the error message embeds the key name
   !! directly (solute-specific scalar-fallback pattern).
   subroutine read_real_array(sec, key, arr, errors)
      type(toml_table), pointer, intent(in)    :: sec
      character(len=*),          intent(in)    :: key
      real(real64), allocatable, intent(out)   :: arr(:)
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
                               "non-real cell", "solute."//trim(key))
            if (allocated(arr)) deallocate(arr)
            return
         end if
         arr(i) = val
      end do
   end subroutine read_real_array

end module read_solute_toml_mod
