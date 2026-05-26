---
title: Sub-state aliasing with ASSOCIATE
author: SWAP modernization team
---

# Sub-state aliasing with ASSOCIATE

## Context

After the state-rescue arcs, almost every physics routine receives the
aggregator `swap_state_t` and reaches into named sub-records: `state%soilwater`,
`state%drainage`, `state%timecontrol`, `state%cfg%drain`. The full chains read
poorly in tight numerical code — five-line equations turn into ten-line
equations once every variable is prefixed twice. Fortran's `ASSOCIATE` gives
scoped aliases without pointer overhead and without changing intent, so a
routine can open one `associate` block at the top of its body and write the
rest of the code in terms of short names that line up with the subsystem
labels we already use in commit messages and documentation.

The convention below is the pattern the rescue settled on once the
sub-record split was complete (drainage, soilwater, surfacewater, atmosphere,
crop, heat all use it). New code touching multiple sub-records of
`swap_state_t` should follow it. Short single-record blocks (one or two
`%` chains in the body) do not need an `associate`.

## Canonical example

`src/drainage/drainage.f90` is the reference implementation. The opening of
`bocodre` aliases four state sub-records and three config sub-records in one
block, then the entire 150-line body reads as if those names were locals:

```fortran
subroutine bocodre(dh, state)
   real(8) :: dh
   type(swap_state_t), intent(inout) :: state

   integer :: level, imper
   real(8) :: qdrdm, qdratio, swdepth, swexbrd, dvmax, swstmax, wl, rd, re

   associate (drai       => state%drainage,             &
              soil       => state%soilwater,            &
              surf       => state%surfacewater,         &
              time       => state%timecontrol,          &
              sw_cfg     => state%cfg%surface_water,    &
              drain_cfg  => state%cfg%drain,            &
              runoff_cfg => state%cfg%drain%surface_runoff)

      if (soil%gwl .gt. 998.0d0) then
         do level = 1, drai%nrlevs
            drai%qdrain(level) = 0.0d0
         end do
         return
      end if

      ! ... 150 lines reading drai%, soil%, surf%, sw_cfg%, drain_cfg%, ...

   end associate
end subroutine bocodre
```

Two things to notice. First, the block wraps the **entire** routine body,
not a tight inner loop — the routine reads from six places, and aliasing
once at the top is what makes the rest legible. Second, both the dynamic
state sub-records (`drai`, `soil`, `surf`, `time`) and the read-only config
sub-trees (`sw_cfg`, `drain_cfg`, `runoff_cfg`) are aliased side-by-side;
splitting them across two `associate` blocks adds nesting without buying
anything.

## Canonical alias names

Use the same short name for the same sub-record across every file. This is
not formally enforced, but it is consistent in the modernized code today
and matches the labels used in commit messages and ADRs. Reviewers should
push back on bespoke names.

| Sub-record                     | Alias    | Notes                              |
| ------------------------------ | -------- | ---------------------------------- |
| `state%mesh`                   | `mesh`   | Numerical grid                     |
| `state%soilwater`              | `soil`   | Soil water state                   |
| `state%atmosphere`             | `atmo`   | Meteorology + interception         |
| `state%crop`                   | `crop`   | Crop state                         |
| `state%drainage`               | `drai`   | Drainage-system state              |
| `state%surfacewater`           | `surf`   | Surface water state                |
| `state%heat`                   | `heat`   | Soil heat                          |
| `state%solute`                 | `solu`   | Solute transport (where present)   |
| `state%timecontrol`            | `time`   | Time loop + control flags          |

For config sub-trees, the established form is `<area>_cfg`, e.g.
`drain_cfg => state%cfg%drain`, `sw_cfg => state%cfg%surface_water`,
`runoff_cfg => state%cfg%drain%surface_runoff`. A few earlier files use
the inverse `cfg_<area>` order (`cfg_heat`, `cfg_crop`); both forms are
present in the tree today. Prefer `<area>_cfg` for new code and leave the
older `cfg_<area>` files for incidental cleanup.

When a routine reaches into a nested cohort of a sub-record (e.g. the
intermediate and cumulative cohorts of `state%atmosphere`), aliasing both
the parent and the children at once is acceptable and keeps the chains
short:

```fortran
associate (atmo => state%atmosphere,      &
           intr => state%atmosphere%intr, &
           cumu => state%atmosphere%cumu, &
           heat => state%heat)
```

## What `ASSOCIATE` does and does not change

The aliases exist only inside the `associate` block. Once control leaves
the block — including via a subroutine call inside the block — the alias
names are no longer visible to the callee, which sees the underlying
`state%X` chains again. They are a textual affordance, not a pointer or
reference type.

`ASSOCIATE` does not change argument intent. If the enclosing routine has
`type(swap_state_t), intent(inout) :: state`, then writes through the
aliases (`drai%qdrain(level) = …`) are mutations of `state%drainage%qdrain`
and require the `inout` on the dummy argument; the compiler will reject
the assignment otherwise. There is no special read-only `associate`
construct — read-only is enforced at the dummy argument, not at the
alias.

Aliasing an `allocatable` component is fine; the alias refers to the
component, and `allocated(drai%qdrain)` works as expected. Aliasing
something that is itself an `allocatable` array (rather than a record
containing one) is also fine but rarely needed — you can usually alias
the enclosing record and read fields off it.

## When to use it, when to skip

Use `ASSOCIATE` when the routine repeatedly references **two or more**
sub-records of `state` or `config`, and the body is long enough that the
`%` chains noticeably hurt readability. Most ported physics routines hit
this bar; many of them reach into four to six sub-records, which is where
the pattern pays off most.

Skip `ASSOCIATE` for very short routines where the chains appear once or
twice — `state%timecontrol%dt` once in a one-line guard does not need an
alias. `src/atmosphere/meteodt.f90` line 61 (`associate (time =>
state%timecontrol)`) is a borderline case that is allowed for consistency
with the rest of the file, but a standalone tiny routine is fine without
it.

`ASSOCIATE` is not a substitute for extracting a helper subroutine. If a
block is long enough that the alias list spans seven or eight entries
**and** the block has its own coherent algorithm, that algorithm should
become its own routine that takes the relevant slice of state by argument
— with the `associate` block moved inside the new routine if it still
helps. The drainage example above is at the upper end of what fits
comfortably in one block; anything significantly larger usually wants to
be decomposed first.

## Placement

The `associate` opens immediately after the local variable declarations
and the executable preamble (early-return guards, parameter unpacking)
and closes just before `end subroutine`. The body of the routine sits
inside it. Indentation follows the same rule as any other Fortran block:
the body is indented one level relative to `associate`, and `end
associate` lines up with the opening keyword.

The legacy `timecontrol.f90` ASSOCIATE binding `SAVE` locals to
`time_state_t` fields was the first appearance of the pattern in the
rescue and is still in the tree; it predates the sub-record convention
described here. New code should use sub-record aliasing as shown above
and not introduce new `SAVE`-binding blocks.

## Reference

- `src/drainage/drainage.f90` — `bocodrb`, `drainage`, `bocodre` all use
  the full sub-record aliasing pattern over six or seven entries.
- `src/soil/waterbalance.f90` — multiple `associate` blocks of varying
  width in one file; useful for seeing how the alias list scales with
  routine scope.
- `src/atmosphere/snow.f90` — example of aliasing nested cohorts
  (`intr`, `cumu`) alongside their parent sub-record.
- [`code-style.md`](code-style.html) — general ASSOCIATE policy in the
  context of other code-style conventions.
- [`state-management.md`](state-management.html) — `swap_state_t`
  composition and the per-domain state lifecycle that this pattern
  layers on top of.
