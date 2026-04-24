---
title: Validation
author: SWAP modernization team
---

# Validation

## Two layers

**Primitive checks** live in `validation_mod`
(`src/validation/validation.f90`): range, enum, ordered-pair,
not-empty, non-negative. Each is stateless and appends to a
supplied `error_collection_t`.

**Section and aggregate validators** live on the config types
themselves as type-bound procedures:

    call config%drain%validate(errors)
    call config%validate(errors)

`validate` is read-only on the config. `finalize` (a separate
type-bound procedure on the same types) is the only thing that
mutates the config post-parse.

## Primitive reference

| Primitive | Checks |
|---|---|
| `check_int_range(value, low, high, context, errors)` | `low <= value <= high` |
| `check_real_range(value, low, high, context, errors)` | same, real64 |
| `check_int_enum(value, allowed, context, errors)` | `value` in `allowed` |
| `check_not_empty(value, context, errors)` | `len_trim(value) > 0` |
| `check_nonnegative_real(value, context, errors)` | `value >= 0` |
| `check_ordered_pair(low, high, low_name, high_name, context, errors)` | `low <= high` |

Every primitive takes `context` (e.g. `"drainage.basic.general.dramet"`)
and appends the context into the resulting `error_t%context`, so error
messages point at the offending field path.

## Writing a section validator

    subroutine drainage_config_validate(self, errors)
       class(drainage_config_t),  intent(in)    :: self
       type(error_collection_t),  intent(inout) :: errors

       call check_int_enum(self%dramet, [1, 2, 3], "drainage.dramet", errors)

       if (self%dramet == 2 .and. self%swdivd == 0) then
          call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                             "swdivd must be 1 when dramet=2", &
                             "drainage", is_fatal=.true.)
       end if
    end subroutine

## Writing the aggregate validator

`swap_config_t%validate` first delegates to every section validator,
then runs cross-section invariants:

    subroutine swap_config_validate(self, errors)
       class(swap_config_t),      intent(in)    :: self
       type(error_collection_t),  intent(inout) :: errors

       call self%general%validate(errors)
       call self%simulation%validate(errors)
       call self%meteo%validate(errors)
       call self%drain%validate(errors)
       call self%soil%validate(errors)
       call self%crop%validate(errors)

       ! Cross-section rules go here.
    end subroutine

## Testing

One pFUnit test per documented rule. Each test builds a minimal
config with the rule's inputs, calls `validate`, asserts the
expected error code and context. See `tests/unit/config/` for
examples.
