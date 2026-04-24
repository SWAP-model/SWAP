---
title: "Phase 4a — Configuration Pipeline and Cross-Cutting Infrastructure"
author: Mateusz Zawadzki
date: 2026-04-24
status: draft
---

# Phase 4a: Infrastructure

A pre-Phase-4 diversion from the rescue spec that lands the cross-cutting infrastructure Phase 4 depends on — logger review, error handling, validation, composable TOML reader — before any physics refactor starts.

## Background

The rescue spec (`docs/superpowers/specs/2026-04-22-rescue-and-stabilize-design.md`) ends Phase 3 with full unit-test coverage of every `*_state_t`, both TOML readers, and a representative set of pure physics routines. Phase 4 as written starts by "purging `use variables` residue" and porting modules to the state-type / `intent(in)` pattern. That sequencing assumes the infrastructure for error reporting and input validation already exists. It does not.

Today:

- `src/error/` is an empty placeholder (README only).
- `src/validation/` does not exist.
- `src/config/` does not exist; `swap_config_t` does not exist. Input is read directly into `swap_state_t` by legacy `readswap.f90` (5161 lines) via ttutil, mixing deserialization, transformation, and validation.
- `swap_log.f90` exists and is complete but untested and undocumented.
- Two partial TOML readers exist (`readswaptoml.f90`, `readdrainagetoml.f90`), of which only the second is on a reusable pattern; both had null-pointer bugs surfaced and fixed during Phase 3.
- 275 `fatalerr` call sites across `src/` propagate errors by halting.

Phase 4a builds the infrastructure, tests it against the `hupselbrook` case end-to-end, authors TOML projects for all six regression cases, and leaves the legacy path untouched. Phase 4 then inherits a working, tested substrate and ports physics modules onto it incrementally.

## Scope

### Goal

A fully tested, runnable parallel path that loads a SWAP project from TOML, validates it, finalises it, and can be diffed against the legacy reader's output. For `hupselbrook` the parity diff is full (field-by-field). For the other five cases, authored TOML projects exist and smoke-parse cleanly.

### In scope

- New `src/error/` module with `error_t`, `error_collection_t`, typed error codes, logger integration.
- New `src/validation/` module with shared validator primitives (range, enum, ordered-pair, …).
- New `src/config/` module hierarchy with per-section `*_config_t` types and an aggregate `swap_config_t`, each with `validate` and `finalize` type-bound procedures.
- New `src/io/toml/` subtree with composite thematic TOML reader (field helpers → section readers → document dispatcher) and a matching emit writer.
- Temporary adapter `src/core/config_to_state.f90` that populates legacy `swap_state_t` from `swap_config_t`.
- Audit and documentation of existing `src/core/swap_log.f90`. No code change unless the audit finds a concrete bug.
- TOML project files for all six regression cases under `tests/swap-cases/toml/`.
- pFUnit unit tests for every new procedure.
- One field-by-field parity test for `hupselbrook`. Smoke tests (load + validate + finalize produces zero fatal errors) for the other five.
- One round-trip test for `hupselbrook` (`load → emit → load` produces identical config).
- Docs: `docs/logging.md`, `docs/error-handling.md`, `docs/validation.md`, `docs/configuration-schema.md` (update), two ADRs, per-subdir READMEs.

### Non-goals

- **No physics changes.** Not one line of non-infrastructure `src/` code is modified in Phase 4a. Physics-carrying modules (`src/atmosphere/*.f90` except `atmosphere_state.f90`, `src/soil/*.f90` except `soil_state.f90`, etc.) are read-only for the duration of this phase.
- **No legacy reader replacement.** `readswap.f90`, `readdra.f90`, `readcrop*.f90`, `readmeteo.f90` stay as they are. `readswap` is called as before during the regression suite.
- **No wiring the new path into main execution.** `swap_main.f90` is untouched. The new TOML path runs only inside pFUnit tests during Phase 4a.
- **No state-layout change.** `swap_state_t` and its sub-types stay as they are. The adapter is the only bridge.
- **No `fatalerr` migration in legacy modules.** The 275 existing sites stay. New infrastructure uses the new error types. Legacy sites migrate during Phase 4+ as each module is touched (matches `docs/code-style.md` §Legacy rule).
- **Full parity tests for cases 2–6.** Deferred to a pre-Phase-4 follow-on step after Phase 4a lands.
- **No Intel/ifx re-enablement, no multicore, no Python bindings.** Still gfortran-only.
- **No regression fixture regeneration.** `check-full` stays the canonical behavioural gate.
- **Retiring the two existing TOML readers.** `readswaptoml.f90` and `readdrainagetoml.f90` and their pFUnit suites disappear; the new architecture replaces them.

## Architecture

### Five new top-level subtrees

```
src/
├── error/                     NEW   error_t, error_collection_t, typed codes
├── validation/                NEW   shared validator primitives
├── config/                    NEW   swap_config_t + per-section *_config_t
├── io/toml/                   NEW   composite thematic reader + emit + field helpers
├── core/
│   ├── swap_log.f90                 EXISTING  reviewed, documented; no code change
│   └── config_to_state.f90    NEW   TEMPORARY adapter
```

### Data flow on program start

```
TOML bytes ──▶ Parse ──▶ config (raw) ──▶ Validate ──▶ Finalize ──▶ Adapter ──▶ state (legacy)
                  │                             │             │
                  └─▶ PARSE errors ─────▶ error_collection_t ◀─┘
                                                 │
                                                 ▼
                                        abort_if_fatal ──▶ log + error stop
```

Every stage writes to the same `error_collection_t`; the single `abort_if_fatal` checkpoint runs after finalize completes.

During Phase 4a this flow runs only inside pFUnit tests, not in `swap_main.f90`.

### Principles

- **Composability.** Each stage's input and output types are well-defined; any stage can be driven independently (e.g., eventual Python bindings construct `swap_config_t` in-memory, skip parse, run validate + finalize + adapter).
- **Error accumulation.** Parse, validate, and finalize all append to a shared `error_collection_t`. No early aborts inside those stages. Single abort point after finalize.
- **Logger is the universal sink.** Every `append` on `error_collection_t` also calls `swap_log%log_error` with the formatted message. No "remember to also log" pitfalls.
- **Legacy untouched.** Zero edits to any `src/**` file that isn't in the new infrastructure modules (plus `src/core/swap_log.f90`, which is reviewed but not changed unless the audit finds a concrete bug).

## Module design

### `src/error/error.f90`

```fortran
module error_mod
   use swap_log, only: log_error
   implicit none
   private

   integer, parameter, public :: ERR_NONE                     = 0
   integer, parameter, public :: ERR_IO_READ_FAILED           = 100
   integer, parameter, public :: ERR_IO_WRITE_FAILED          = 101
   integer, parameter, public :: ERR_PARSE_MALFORMED_TOML     = 200
   integer, parameter, public :: ERR_PARSE_TYPE_MISMATCH      = 201
   integer, parameter, public :: ERR_PARSE_MISSING_REQUIRED   = 202
   integer, parameter, public :: ERR_VALIDATION_OUT_OF_RANGE  = 300
   integer, parameter, public :: ERR_VALIDATION_ENUM          = 301
   integer, parameter, public :: ERR_VALIDATION_CROSS_FIELD   = 302
   integer, parameter, public :: ERR_VALIDATION_CROSS_SECTION = 303
   integer, parameter, public :: ERR_FINALIZE_DERIVATION      = 400
   integer, parameter, public :: ERR_ADAPTER_UNSUPPORTED      = 500

   type, public :: error_t
      integer                       :: code    = ERR_NONE
      character(len=:), allocatable :: message
      character(len=:), allocatable :: context  ! e.g. "drainage.basic.general.dramet"
      logical                       :: is_fatal = .false.
   end type

   type, public :: error_collection_t
      type(error_t), allocatable :: items(:)
   contains
      procedure :: append          ! appends + calls log_error
      procedure :: has_errors      ! any items?
      procedure :: has_fatals      ! any item with is_fatal?
      procedure :: count
      procedure :: summary         ! multi-line formatted report (string)
      procedure :: abort_if_fatal  ! if has_fatals(), write summary then error stop
      procedure :: clear           ! deallocate items
   end type
end module
```

**Conventions.** Every fallible procedure has signature `(..., errors)` with `type(error_collection_t), intent(inout) :: errors`. `append` auto-logs via `swap_log%log_error`. Callers never log manually. "Fatal" means execution cannot meaningfully continue past the next checkpoint. Phase 4a treats all parse errors and most validation errors as fatal; `is_fatal=.false.` exists so finalize or future code can report non-blocking warnings without upgrading them into halts.

### `src/validation/validation.f90`

Thin, reusable helpers that section validators call. No state, no side effects except appending to errors:

```fortran
module validation_mod
contains
   subroutine check_int_range(value, low, high, context, errors)
   subroutine check_real_range(value, low, high, context, errors)
   subroutine check_int_enum(value, allowed, context, errors)
   subroutine check_not_empty(value, context, errors)          ! strings
   subroutine check_nonnegative(value, context, errors)
   subroutine check_ordered_pair(low, high, low_name, high_name, context, errors)
end module
```

One pFUnit test per primitive, happy + one failing case each.

### `src/config/`

One file per section. Each file defines `<section>_config_t` with `validate` and `finalize` type-bound procedures:

```
src/config/
├── general_config.f90         general_config_t
├── simulation_config.f90      simulation_config_t
├── meteorology_config.f90     meteorology_config_t
├── drainage_config.f90        drainage_config_t
├── soil_config.f90            soil_config_t
├── crop_config.f90            crop_config_t
└── swap_config.f90            swap_config_t  (composes the six above)
```

`validate` does range / enum / cross-field checks, appends errors, mutates nothing. `finalize` runs only if validate produced no fatal errors, performs array expansion and derived-value computation, can itself append errors (typically `ERR_FINALIZE_DERIVATION`).

`swap_config_t%validate` delegates to section validators then runs cross-section invariants. `swap_config_t%finalize` mirrors that pattern: section finalize first, then cross-section derivation.

**Scope for hupselbrook:** every field the legacy `readswap` touches when loading `hupselbrook.swp`/`.dra`/`.crp` has a home on one of the six sections. Other cases may exercise more sections, but hupselbrook defines the minimum.

### `src/io/toml/`

```
src/io/toml/
├── toml_field_helpers.f90     get_required_int, get_optional_real_with_default,
│                               get_array_of_tables, parse_date, ...
│                               (wrap tomlf get_value; null-ptr safe;
│                                auto-append parse errors)
├── load_swap_config.f90       load_swap_config(path, config, errors)
├── write_swap_config.f90      write_swap_config(config, path, errors)
├── read_general_toml.f90      section readers (one subroutine per section)
├── read_simulation_toml.f90
├── read_meteorology_toml.f90
├── read_drainage_toml.f90
├── read_soil_toml.f90
└── read_crop_toml.f90
```

Each section reader has signature `read_<section>_toml(doc_table, config_section, errors)`:

- takes the root `toml_table` (or sub-table for nested sections)
- writes into its section's `*_config_t`
- appends PARSE errors to `errors`
- returns normally; no aborts, no direct halts

`load_swap_config` loads the file, extracts the root table, dispatches to each section reader in order. No section reader calls another; composition happens only at the top level.

`write_swap_config` mirrors the reader: one section-emit subroutine per section, walking a `swap_config_t` and writing a TOML file using toml-f's table-building API. Round-trip invariant (`load → emit → load` produces an identical `swap_config_t`) is enforced by a pFUnit test for hupselbrook.

### `src/core/config_to_state.f90`

Temporary adapter. Populates the subset of `state%*` that corresponds to config-loaded fields. State fields populated by other means (initial conditions, runtime accumulators) are left untouched. Header comment declares it scheduled for retirement when Phase 4 folds state types into config types. Contains no business logic: if a transformation needs logic, the logic lives in `finalize`; the adapter is a dumb field copy.

### `src/core/swap_log.f90`

No code change unless the audit finds a concrete bug. Audit checks:

- Log-level thresholds consistent with error codes.
- Thread-safety considerations flagged (none needed now, but surface anything that would break for future multicore work).
- Missing `to_str` overloads for any type the error module needs.

Audit output: `docs/logging.md` entry + mentions in `docs/index.md` and `src/core/README.md`.

## Test plan

### Unit coverage

| Module | Test scope |
|---|---|
| `error_mod` | `append` appends + auto-logs; `has_errors` / `has_fatals` semantics; `summary` format; `abort_if_fatal` halts when fatal, returns when not; `clear` resets |
| `validation_mod` | one test per primitive: happy path + failing case |
| `*_config_t` | per-section: `validate` catches each documented rule (one test per rule); `finalize` expands correctly; `finalize` idempotent; `finalize` reports errors when cross-field invariants broken |
| `toml_field_helpers` | one test per helper: happy, wrong-type input → `ERR_PARSE_TYPE_MISMATCH`, absent required key → `ERR_PARSE_MISSING_REQUIRED`, absent optional key → default, malformed nested table |
| section readers | per section: fixture with minimum valid content → populated config, zero errors; wrong-types fixture → parse errors; missing-optional fixture → defaults |
| `load_swap_config` | full hupselbrook TOML loads cleanly; intentionally malformed TOML surfaces correct error with file-path context |
| `write_swap_config` | round-trip (`load → emit → load → assert identical config`) for hupselbrook |
| `config_to_state` | one field-by-field adapter test per section for hupselbrook; state fields not touched by config remain at their init defaults |

Target: every new procedure has at least one direct test. Coverage tracked (not gated); Phase 4a re-runs the Phase 3 coverage baseline at exit to confirm no per-domain drop and to record the new numbers.

### Parity test — the headline deliverable

One pFUnit integration suite: `tests/unit/io/toml/test_hupselbrook_parity.pf`. Loads hupselbrook through both paths, runs the new path through validate / finalize / adapter, and diffs every adapter-written field with `@assertEqual` (tolerance `1.0d-12` on reals, exact on integers). A helper `assert_states_equal_for_config(expected, actual)` deduplicates the per-field assertions so the test body stays readable.

If it passes, the architecture works end-to-end for hupselbrook. This is the single load-bearing integration check for Phase 4a.

### Smoke tests for the other five cases

Per-case pFUnit suite: `load_swap_config(..., errors); call config%validate(errors); call config%finalize(errors); @assertFalse(errors%has_fatals())`. Confirms the reader accepts every real-world case, not just hupselbrook.

Full parity tests for cases 2–6 are deferred to a pre-Phase-4 follow-on step.

### TOML project authoring

For each case, author the TOML equivalent under `tests/swap-cases/toml/<name>/` (same naming as legacy case folders, e.g., `1.hupselbrook`, `2.grassgrowth`, …). Hand-authored, not auto-converted: every case needs one human pass with the legacy `.swp`/`.dra`/`.crp` open side-by-side. Six cases combined: ~2–3 days of mechanical authoring. Lives on the `main` branch of the `swap-cases` submodule; the outer repo bumps the submodule pointer once all six are written.

## Documentation

- `docs/logging.md` — logger review and usage guide.
- `docs/error-handling.md` — error infrastructure guide, code table, migration examples.
- `docs/validation.md` — per-section rules reference (what each section's `validate` checks).
- `docs/configuration-schema.md` — **updated** to reflect the new `*_config_t` layout (currently describes TOML keys only).
- `docs/adr/0007-config-validate-finalize-pipeline.md` — the four-stage pipeline decision.
- `docs/adr/0008-error-collection-over-fatalerr.md` — the error-model decision.
- New top-level subdirs (`src/error/`, `src/validation/`, `src/config/`, `src/io/toml/`) get `README.md` files matching the existing pattern (Responsibility / Public interface / Dependencies).
- `docs/index.md` updated to link all of the above.

## Exit criteria

1. All new modules compile; `pixi run -e test check-fast` stays under 90 seconds.
2. `pixi run -e test check-full` stays green.
3. Every new procedure has at least one pFUnit test.
4. `test_hupselbrook_toml_matches_legacy` passes — TOML path produces a `swap_state_t` identical (within tolerance) to legacy for hupselbrook.
5. Round-trip test passes for hupselbrook: `load → emit → load` → identical config.
6. Smoke test passes for each of the six cases: load + validate + finalize produces zero fatal errors.
7. All listed docs + ADRs committed; `pixi run -e docs docs-build` succeeds.
8. Coverage baseline re-run and recorded; no per-domain drop vs. Phase 3 baseline.
9. `main` fast-forwarded to `development` locally; `rescue/phase-4a-infrastructure` tag applied; no push.
10. `readswaptoml.f90` and `readdrainagetoml.f90` and their pFUnit suites removed; `tests/unit/io/test_readswaptoml.pf` and `tests/unit/io/test_readdrainagetoml.pf` removed; `tests/unit/io/fixtures/` refitted to the new architecture.

## Git workflow

- Local-only. Nothing pushed to `origin/main` for the duration of Phase 4a (continues the rescue-era policy).
- All commits on `development`. No per-task feature branches during Phase 4a (matches Phases 0–3 workflow; Phase 4's per-change branching starts with Phase 4 proper).
- `main` fast-forwards to `development` at phase exit, locally only.
- Tag `rescue/phase-4a-infrastructure` applied at exit.

## Risks and mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| Hand-writing six TOML files burns a week | Medium | Author hupselbrook first + write parity test; the test exposes every missing TOML field, so the other five follow the same reference once the schema is proven. |
| Legacy `readswap` has section-implicit defaults that aren't obvious in the TOML | Medium | Parity test surfaces them; field-by-field diff pinpoints exact slot; update `*_config_t` default or `finalize` as needed. One such quirk per section is realistic. |
| `config_to_state` adapter grows a special case that's really a state-layout bug | Low | Spec rule: adapter contains no business logic. If a transformation needs logic, the logic goes in `finalize`; adapter is dumb field copy. |
| Round-trip emit reveals toml-f quirks (array-of-tables ordering, scalar vs single-element array ambiguity) | Medium | The round-trip test is the surfacing mechanism; dedicate an extra ADR slot if quirks turn up and need documenting. |
| Six new subdirs + new `config_t` types + new tests = big PR | Always, for infrastructure phases | Same as Phases 2/3: commit per task, land on `development` only, tag at phase end. No per-task feature branches. |
| Phase 4a balloons as section schemas reveal surprises | Medium | Scope is hupselbrook-complete, not SWAP-complete. Fields not exercised by any of the six cases are out of scope; add them in Phase 4 when legacy paths port over. |

## Follow-on work

Deferred to a pre-Phase-4 step after Phase 4a tags:

- Full field-by-field parity tests for cases 2–6.
- Potential `write_swap_config` round-trip tests for cases 2–6 (lift the hupselbrook round-trip template).

Deferred to Phase 4 proper:

- Wiring `load_swap_config` → `config_to_state` into `swap_main.f90` replacing the legacy call.
- Retiring `readswap.f90` / `readdra.f90` / crop readers and the ttutil dependency.
- Folding `*_config_t` + `*_state_t` into the compartment-based state refactor.
- Migrating legacy `fatalerr` sites to the error collection as each module is touched.

## Definition of done for this phase

- `rescue/phase-4a-infrastructure` tag applied locally.
- Every exit criterion above verified.
- `docs/` complete and rendered by FORD.
- Zero changes in any existing `src/**` file outside the new infrastructure subtrees. Verified by `git diff rescue/phase-3-coverage..HEAD -- src/atmosphere/ src/soil/ src/crop/ src/boundary/ src/drainage/ src/macropore/ src/solute/ src/heat/ src/utils/` — output must be empty. If the adapter discovers a mismatch between `swap_config_t` and `swap_state_t` that would need a state-layout change, the mismatch is resolved by adjusting `swap_config_t` or `finalize`, not the state type. State-layout evolution happens in Phase 4.
