---
title: Code style
author: SWAP modernization team
---

# Code style

## Language target

Fortran 2008 is the target for new code. The build passes `-std=legacy`
(see `gfortran_flags` in `meson.build`) so legacy idioms in older modules
continue to compile during the rescue, but new code should not rely on
anything that requires that flag.

The rescue explicitly does not mandate a wholesale rewrite of legacy code.
Modernize modules incrementally as they are touched. A module that hasn't
been touched in Phase 2–3 is not a bug; it's waiting for Phase 4.

## Compiler flags

The flag set lives in `meson.build` under `gfortran_flags` and is applied
project-wide via `add_project_arguments(...)`. Each flag, in order:

- `-O2` — optimisation level; unchanged from the legacy build.
- `-ffree-line-length-none` — removes the 132-character line-length limit
  on free-form source. Needed by a handful of legacy files; new code
  should still aim for short lines.
- `-Wno-line-truncation`, `-Wno-compare-reals`, `-Wno-missing-include-dirs`,
  `-Wno-unused-variable`, `-Wno-unused-dummy-argument` — warnings
  suppressed during the rescue because the legacy code would otherwise
  produce hundreds of diagnostics and drown out anything actionable.
  Phase 4 gradually re-enables them as modules are cleaned up.
- `-std=legacy` — accept legacy Fortran idioms (Hollerith constants,
  non-standard DO termination, etc.) in legacy modules. New code should
  not need this; if a new file only compiles under `-std=legacy`, treat
  that as a code-review blocker.
- `-finit-local-zero` — zero-initialise local variables, matching Intel's
  `-init=zero`. Without this, `salinitystress` produces NaN values and
  `macropore` drifts substantially (see ADR 0001).

Compiler: gfortran only (ADR 0001). Any non-GCC Fortran compiler is
rejected at meson configure time.

## Module conventions

One module per file. The module name matches the filename stem plus the
`_mod` suffix: `src/atmosphere/atmosphere_state.f90` contains
`module atmosphere_state_mod`. `implicit none` goes immediately after the
`module` statement. Default visibility is `private`; every exported name
is re-declared with an explicit `public :: thing` line. For example:

```fortran
module atmosphere_state_mod
    implicit none
    private

    public :: atmosphere_state_t
    public :: atmosphere_state_init
    public :: atmosphere_state_finalize
    ...
```

No `COMMON` blocks and no `include` of legacy `.fi` headers in new code.
If a legacy module still uses `COMMON` or `include`, leave it until
Phase 4 touches it — do not refactor proactively. The goal is net
progress, not churn.

## Types and intent

Use derived types for config, initial, and state:

- `type :: <domain>_config_t` — values loaded from input that don't change
  during a run.
- `type :: <domain>_initial_t` — values required to set the initial
  condition.
- `type :: <domain>_state_t` — everything that evolves during a run.

Every procedure argument has an explicit `intent(in)`, `intent(out)`, or
`intent(inout)`. Choose the weakest one that fits:

- `intent(in)` for configuration-like values that the routine must not
  mutate.
- `intent(inout)` for state that the routine updates in place.
- `intent(out)` for purely-computed outputs.

For pure computations, prefer `intent(in)` arguments plus a function
return value.

Kind specifiers: `real(8)` and `real(kind=8)` are acceptable during the
rescue — that is what the legacy code uses and mass-renaming adds risk.
New code should prefer `real(kind=real64)` with
`use iso_fortran_env, only: real64` for clarity and portability.

## Naming

Names follow a simple, predictable scheme so that readers can locate the
defining module from a symbol alone.

- Modules: `<domain>_mod` (e.g. `atmosphere_state_mod`, `drainage_state_mod`).
- Types: `<domain>_state_t`, `<domain>_config_t`, `<domain>_initial_t`.
- Procedures on a type: `<domain>_state_init`, `<domain>_state_finalize`,
  `<domain>_state_reset_cumulative`, `<domain>_state_reset_intermediate`.
- Domain processes: `<domain>_<verb>` (e.g. `surfacewater_wballev`,
  `atmosphere_process_rain_events`).
- Variables: `snake_case`.

Prefer full words unless the abbreviation is a physics convention —
`theta` for volumetric water content, `psi` for pressure head, `K` for
hydraulic conductivity, `h` for head, `q` for flux, `ET` for
evapotranspiration components. Short notation that matches a textbook
equation is a feature; cryptic abbreviations that mean nothing outside
the original author's head are not.

Derived-type fields follow the same rule: use short physics-aligned
names when they match textbook notation (`theta`, `psi`, `tav`,
`rad`), and descriptive snake_case names otherwise
(`interception_capacity`, `swetr`).

## Control flow

Prefer `select case` over chained `if/elseif` when branching on a scalar.
Early returns are fine — the codebase does not enforce a single-exit-point
rule, and forcing one usually adds a flag variable that hurts readability
more than the extra `return` does.

`goto` is forbidden in new code. Legacy `goto` may stay until Phase 4
touches the module; do not refactor it proactively. Use named `do` loops
with `exit` and `cycle` when flow control inside nested loops becomes
non-obvious.

Always use `end subroutine foo` / `end function foo` / `end module foo`
rather than the bare `end` terminator — it makes mismatched blocks much
easier to spot.

## ASSOCIATE

Use `ASSOCIATE` for scoped aliasing when a block repeatedly references
fields of a state type. Example from `src/core/timecontrol.f90`:

```fortran
associate( datea     => tc_datea,     &
           nextyear  => tc_nextyear,  &
           tchange   => tc_tchange,   &
           dtEvent   => tc_dtEvent,   &
           tEvent    => tc_tEvent )
    ! ... use datea, nextyear, tchange, dtEvent, tEvent here ...
end associate
```

The alias names exist only inside the `associate` block. They are not
persistent aliases and they do not change the storage of the underlying
variables; they are a local readability affordance.

`ASSOCIATE` is not a substitute for extracting a helper subroutine. If a
block is long enough to need many aliases, it is usually long enough to
be its own routine that takes the state as an argument — pull the block
out first and use `ASSOCIATE` inside the new routine if it still helps.

## Comments and docstrings

Documentation is generated by FORD (the Fortran documenter). FORD
consumes `!>` docstring comments immediately above a procedure, type, or
module, and `!!` comments on declarations for per-field/per-argument
descriptions. Every new public procedure gets:

- a one-line `!>` summary immediately above its declaration, and
- per-argument `!!` descriptions on each dummy argument line.

See `src/atmosphere/atmosphere_state.f90` for the current template.

Comments should explain **why**, not **what**. If a reader can tell what
a line does by reading it, the comment adds nothing. Non-obvious
invariants, unit conventions (cm/d vs. m/s, °C vs. K), and references to
physics papers or textbooks are always worth writing down.

## Formatting

**`pixi run lint` (fprettify) is authoritative.** The rules are encoded
in `.fprettify.rc` at the repository root. Run `pixi run lint-check`
first to see the diffs the formatter would make without applying them;
`pixi run lint` rewrites files in place.

If the formatter disagrees with something written here, `.fprettify.rc`
wins — open an issue to reconcile rather than ignoring the formatter.
The linter is deliberately conservative during the rescue:
`enable-replacements` is off so that legacy operator forms (`.eq.`,
`.le.`, `.ne.`) are not auto-rewritten to `==`, `<=`, `/=`, and
`disable-indent-mod` is on so that hand-crafted module-scope alignment
is preserved. These choices will be revisited at Phase 4 exit, when the
legacy surface area is small enough that aggressive rewrites are safe.

Two house rules that exist independent of the formatter:

- Always indent with spaces; never tabs. Tab characters are not part of
  the Fortran character set and display inconsistently across editors.
- Case: lowercase for keywords, identifiers, and intrinsics. The
  Fortran language is case-insensitive, but consistent case makes
  grep-based navigation predictable.

## Tests

Every new public procedure comes with a pFUnit test in
`tests/unit/<domain>/`. Tests live close to their subject module: a test
for `src/atmosphere/foo.f90` goes under `tests/unit/atmosphere/`.

See `docs/build-and-test.md` for the pFUnit harness layout and the
`pixi run` gates (`check-fast`, `check-full`) that the CI pipeline
enforces, and `docs/contributing.md` for the fixture policy
(what to put under `tests/regression/` versus `tests/unit/`, and what
must **not** be regenerated casually).

## Legacy rule

`src/core/variables.f90` and `src/core/swap_state_sync.f90` are
rescue-era scaffolding. Do not extend them.

Concretely:

- If a new feature would need to add a variable to `variables.f90`, add
  the field to the appropriate `<domain>_state_t` and route access
  through `swap_state_t` instead.
- If a physics routine still reads from `variables` via the sync bridge,
  port it to take the relevant state as an argument **when you touch it**.
  Not as a separate cleanup task, and not as a condition for merging
  unrelated changes — incrementally, during normal work.

Phase 4's stated goal is to shrink `swap_state_sync.f90` to zero lines.
Every change that moves a field off the legacy module and onto a
`<domain>_state_t` is progress toward that; every change that adds a new
field to `variables.f90` is regress and will be rejected in review.
