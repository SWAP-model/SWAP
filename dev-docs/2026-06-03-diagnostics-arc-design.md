# Diagnostics arc — design spec

**Date:** 2026-06-03
**Status:** design (approved direction; pending spec review → writing-plans)
**Branch:** `logging-review` (worktree off `development`)
**Companion:** `dev-docs/logging-review-2026-06-03.md` (the infrastructure review that
motivates this)

---

## 1. Goal

Make logging, warnings, and error-raising a single coherent, per-instance system so
that SWAP can run as a CLI **and** as many in-process columns under pyswap / BMI /
XMI / a coupler, where:

1. A critical physics violation **gracefully terminates the offending instance's run**
   without killing the host process or the sibling columns.
2. The terminated instance **returns a collection of errors** identifying *where* the
   problem was (subsystem/context, simulation date + node/compartment, the offending
   value and bound).
3. Diagnostics are **user-configurable** (level, file/sink, routing) without
   recompiling, from CLI, ops, and an embedding host.
4. Every change is **byte-identical** against `swap420gf` (logging touches only
   stderr/file/callback, never the compared output files or kernel numerics).

This **keeps** the deliberate "stop on a critical issue" policy. It changes the
*scope* (this instance, not the process), the *mechanism* (clean termination at a
step boundary, not `error stop` mid-kernel or continue-with-garbage), and makes the
behaviour *uniform across entry points*.

## 2. Decisions (from brainstorming)

| # | Decision | Choice |
|---|----------|--------|
| Scope | how much of the arc | **Full arc** (review items 1–6), decomposed into phased sub-specs A–E |
| Instance model | how diagnostics attach to an instance | **Threaded on `state`** (`state%diag`), reached via the already-threaded `swap_state_t` |
| Config surface | who sets level/file/routing, precedence | **All three**, precedence **C-API > env > TOML > built-in default** |
| Fatal mechanism | how a kernel stops the step | **Record on `state%diag` + early-return-on-detect + check at substep boundary**; entry point translates |
| Restructuring | logging/error/validation | **Authorized**; see §6 |

## 3. Problem recap (why, condensed from the review)

- The logger is initialised **only** in the CLI (`swap_main.f90`), hardcoded
  `INFO` + `swap_swap.log`; BMI/C-API/XMI never call `log_init`, so embedded runs get
  no log file but still spam the host's stdout. No user verbosity control anywhere.
- ~150 `fatalerr_collected` sites, **~105 inside compute kernels** (crop 61,
  soilwater 21, drainage 13, atmosphere 6, solute 3, core 4, heat 1). These are
  legacy **fire-and-abort** via the `global_errors` singleton, *not* the threaded
  accumulate-then-abort-at-the-edge discipline (which only reached `io/`/`config/`/
  `validation/`).
- Same fatal condition behaves three ways by entry point: CLI → clean `error stop`;
  XMI (`library_mode=.true.`) → **continue-with-garbage** (the kernel returns and keeps
  computing on the invalid state); direct BMI/C-API (never sets the flag) →
  **host-kill**. Two of three are wrong.
- Logger + error state are **process-global** (`current_level`, `log_unit`,
  `global_errors`, `fatal_was_raised`), so two columns share one file, one level, one
  fatal flag — incompatible with in-process multi-instance.
- ~17 stray `write(*,*)`/`stop`/magic-unit diagnostics bypass the logger; WARN/INFO/
  ERROR all go to stdout (not stderr); records carry wall-clock time but not
  simulation time; coverage is thin (heat/solute log nothing).

## 4. Architecture

### 4a. `diagnostics_t` on `state` (the spine)

A new per-instance record **type** defined in a foundation module
`src/error/diagnostics.f90` (it bundles `error_collection_t` from `error_mod` with the
sink primitives from `swap_log`, so it must sit *above* both and *below* `state`).
`swap_state_t` then **holds an instance** of it as `state%diag`. The
**process-default instance** also lives in this module, so the legacy free
`log_*`/`fatalerr` bridge can reach it without a dependency cycle
(`swap_log` ← `error_mod` ← `diagnostics` ← `state`):

```
type :: diagnostics_t
   integer                       :: level         = LOGLEVEL_INFO   ! min level emitted
   ! routing / sink
   logical                       :: to_stdout     = .false.         ! default OFF (CLI turns on)
   logical                       :: to_stderr     = .true.          ! WARN+ERROR destination
   integer                       :: file_unit     = -1              ! per-instance log file, or -1
   logical                       :: timestamps    = .false.         ! wall-clock; default OFF
   character(len=:), allocatable :: instance_id                     ! stamped on every record
   procedure(diag_sink_i), pointer, nopass :: sink => null()        ! optional host callback
   ! accumulation
   type(error_collection_t)      :: errors                          ! per-instance, accumulates
   logical                       :: fatal_raised   = .false.        ! sticky, per-instance
   ! sim-time context (updated each step by the driver; auto-stamped onto records)
   character(len=11)             :: sim_date       = ''
   integer                       :: daynr          = 0
   integer                       :: daycum         = 0
contains
   procedure :: debug, info, warn, error      ! leveled emit (context, message)
   procedure :: fatal                         ! record fatal + set fatal_raised (no abort here)
   procedure :: aborted                        ! .true. if fatal_raised
   procedure :: take_errors                    ! move-out the collection for the host
end type
```

Because compute code already receives `state`, kernels reach diagnostics with **no
new argument**: `call state%diag%warn('BoundBottom', msg)` /
`call state%diag%fatal('calcgwl', msg)`. The record auto-stamps simulation time and
`instance_id`; node/compartment location is included by the call site (in the
message/context), not auto-derived. State-less leaf kernels (the solver,
`WC_K_models`, `QROMBD`, …) take a small `diag` (or `errors`) argument, or return a
status their state-holding caller records. Until migrated, they use a thin
process-default instance (the bridge).

### 4b. `swap_log` demoted to the shared sink

`swap_log` stops owning policy. It keeps: formatting (`LEVEL context: message` with
optional sim-time/instance-id/wall-clock), `to_str`, and the low-level writers
(stdout / stderr / a file unit / a callback). It exposes a **process-default
`diagnostics_t`** used by the CLI and by the legacy free `log_*`/`fatalerr` bridge
during migration. The module-global `current_level`/`log_unit` are removed once
call sites read their level/unit from a `diagnostics_t`.

### 4c. `error_mod` role

Keep `error_t`, `error_collection_t`, the stable error-code parameters, `append`,
`summary`, `has_fatals`, `count`. **Retire at the end of the arc:** `global_errors`,
`fatalerr_collected`, the free `FatalERR`, `library_mode`/`fatal_was_raised` (their
job moves onto `state%diag`). `append` keeps auto-logging through the owning
instance's sink instead of the global logger.

### 4d. `validation_mod`

Unchanged in spirit — validators append to a passed `error_collection_t`. They now
naturally target `state%diag%errors` (or a load-phase collection) and so participate
in the same accumulation.

### 4e. Routing rules (terminal vs file)

| Level | stdout | stderr | file / sink |
|-------|:------:|:------:|:-----------:|
| ERROR | no | **yes** | yes |
| WARN  | no | **yes** | yes |
| INFO  | CLI only (`to_stdout`) | no | yes |
| DEBUG | no | no | yes (only when `level<=DEBUG`) |

- **CLI** default: `to_stdout=.true.` for INFO + a one-line start/end banner; WARN+
  ERROR to stderr; full record to `<project>.log`.
- **Embedded** default: `to_stdout=.false.` (never pollute the host's stdout); file
  optional (caller-supplied path) and/or `sink` callback the host drains.
- Records carry **simulation time** (date/daynr) and **instance_id**; wall-clock
  timestamps optional, default off (near-useless for fast batch runs).

### 4f. Config surface + precedence

A `resolve_diagnostics_config()` builds the initial `diagnostics_t`, applied in order
(later wins): **default → `[logging]` TOML block → env (`SWAP_LOG_LEVEL`,
`SWAP_LOG_FILE`, `SWAP_LOG_STDOUT`) → C-API/BMI setter**. Every entry point
(`swap_main`, `swap_bmi_mod`, `swap_capi_mod`, `swap_ensemble_mod`) calls it — logging
is initialised everywhere, not just the CLI. `[logging]` is **operational** config:
added to the schema but marked non-physics so it never affects regression.

### 4g. Fatal mechanism — graceful per-instance termination

1. A kernel detects a critical condition, calls `state%diag%fatal(context, msg)` which
   **appends a fatal `error_t`** (code, context, message, stamped with sim-date +
   node/compartment) to `state%diag%errors` and sets `state%diag%fatal_raised`.
2. The detecting routine **returns early from its own body** so it does *not* perform
   the dangerous operation it was guarding (no NaN/OOB). Multiple fatals may
   accumulate within a substep before the next checkpoint.
3. The step driver checks `state%diag%aborted()` at substep boundaries (top of
   `swap_run_step` and the main inner loops) and **stops that instance's step**.
4. The entry point translates the terminated state:

| Entry point | On `aborted()` |
|-------------|----------------|
| CLI (`swap_main`) | write `errors%summary()` to stderr, exit nonzero — **same user-visible result as today** |
| BMI / C-API | `update()` returns nonzero rc; instance marked failed; **host keeps stepping other instances**; `take_errors`/a view exposes the collection |
| XMI / ensemble | per-instance poll (replaces global `library_fatal_raised`); the column drops out, the run continues |

This satisfies the stated goal: the offending column terminates gracefully and hands
back a **collection of errors located by subsystem + sim-time + node**, while siblings
keep running. Once kernels record on `state%diag`, `global_errors` /
`fatalerr_collected` / `FatalERR` / `library_mode` have no callers and are deleted.

## 5. What gets retired (exit criteria for the arc)

- `swap_log` module-global mutable state (`current_level`, `log_unit`, `log_to_*`).
- `error_mod`: `global_errors`, `fatalerr_collected`, free `FatalERR`,
  `library_mode`/`fatal_was_raised`/`set_library_mode`/`library_fatal_raised`/
  `clear_library_fatal`.
- All stray diagnostics: the `TEMPORARY OUTPUT DELETE` writes to units 777/888, the
  `oxygenstress` `write(*,*) '… Too many steps.'` family, the bare `stop` in
  `interception.f90`, magic-unit `write(124/125,…)` in `dormant/sptabulated.f90`,
  and timecontrol progress `write(*,…)` (→ `info`).
- Unused `use … log_*` imports.

## 6. Module restructuring

Minimal, follows ADR 0048 (feature-first behaviour; layer-first shared data):

- **New:** `src/error/diagnostics.f90` — defines the `diagnostics_t` type + methods +
  the process-default instance (above `error_mod`/`swap_log`, below `state`).
  `src/state/swap_state.f90` holds an instance as `state%diag`.
- **`src/core/swap_log.f90`** — narrowed to sink primitives + the process-default
  instance + `to_str`. Stays a core leaf.
- **`src/error/error.f90`** — keeps types/codes/collection; sheds the global/legacy
  bridge by end of arc. Stays its own layer folder.
- **`src/validation/validation.f90`** — unchanged role; stays its own folder.

No new top-level folder; `error/` and `validation/` remain shared layers,
diagnostics-on-state lives in `state/`, the sink stays a `core/` leaf.

## 7. Phasing (each phase is its own implementation plan, byte-identical, check-fast green)

> The implementation plan produced next covers **Phase A only**. Phases B–E each get
> their own plan when their turn comes, so we stay in the loop and can re-scope as the
> physics-path work (D) is learned.

- **A — Foundation & wiring.** Add a `diagnostics_mod` config layer
  (`diagnostics_config_t` + `diag_overrides_t` + `resolve_diagnostics_config` with the
  precedence merge) and the **env + C-API** config surface. Add stderr routing to
  `swap_log` (WARN/ERROR → `error_unit`). Initialise logging in **all** entry points
  (CLI, BMI, C-API, ensemble); default `to_stdout` off when embedded. No kernel
  changes, no `state` schema change. *(Review items 1, 2, routing part of 3.)*
  *(Sequencing note 2026-06-03: the `[logging]` TOML block moved to the start of Phase
  B — it needs new config-schema surface and is lower-value than env/C-API, which
  already cover ops and pyswap.)*
- **B — Diagnostic cleanup & record content.** Add the `[logging]` TOML block (level,
  file, to_stdout, timestamps) feeding the resolver. Route the stray writes/`stop`/magic-unit
  diagnostics through `diag`; delete the `TEMPORARY DELETE` debug; add sim-time +
  instance-id stamping; normalise context strings; drop unused imports.
  *(Items 3, 4.)*
- **C — Per-instance isolation.** Move `swap_log` module-globals onto `state%diag`;
  make the process-default a real default instance; migrate state-holding compute code
  to `state%diag%{debug,info,warn,error}`. Multi-instance log isolation. *(Item 5
  foundation.)*
- **D — Kernel fatal migration (per subsystem).** Convert the ~105 kernel
  `fatalerr_collected` sites to `state%diag%fatal` + early-return + boundary checks,
  one subsystem per commit (order: soilwater → drainage → atmosphere → solute → heat →
  crop, crop last/largest). Add per-instance translation in BMI/C-API/XMI. Retire
  `global_errors`/`fatalerr_collected`/`FatalERR`/`library_mode`. *(Item 6, the core
  of the user goal.)*
- **E — Coverage uplift.** Solver convergence (DEBUG), boundary/switch fallbacks
  (WARN), per-run provenance summary (INFO) where it pays. Leave heat/solute thin
  unless a real need exists.

Dependency: A → B → C → D; E after C. A/B are pure infra and low-risk; D carries the
physics-path risk and is sliced per subsystem so a regression bisects to one subsystem.

## 8. Testing strategy

- **pFUnit unit tests** (register in `tests/unit/testSuites.inc`): config precedence
  resolution; level filtering; routing (WARN/ERROR→stderr, INFO→stdout only on CLI);
  sim-time/instance-id stamping; `diag%fatal` records + `aborted()` true + collection
  contents; **two independent `diagnostics_t` don't cross-talk** (multi-instance);
  callback `sink` receives formatted records (exercises the embedding path).
- **Multi-instance integration test:** two states; force a fatal in one substep of
  instance A; assert A `aborted()` and returns its error collection while B steps
  cleanly and logs to its own sink.
- **Byte-identical regression** (`check-fast` per task, `check-full` at phase tags) —
  logging never touches compared outputs; this is the guard that the migration is safe.
- **Capture-based CLI test:** redirect the default sink to a buffer, assert a known
  fatal still produces the same message + nonzero exit.

## 9. Non-negotiables / risks

- **gfortran only; byte-identical law; no new hidden state; clean rebuild after any
  `src/state/*` schema change** (`diagnostics_t` on `swap_state_t` is a schema change →
  `rm -rf builddir`). Per [[feedback_state_schema_clean_rebuild]] and
  [[feedback_state_component_storage]] keep `diagnostics_t` small / allocatable so it
  doesn't blow the stack.
- **Risk — code between raise and boundary.** Early-return-on-detect covers guard-style
  fatals (the majority). A fatal detected mid-computation where continuing crashes
  before the next checkpoint is the residual case; if any such site exists, that
  routine returns immediately after recording (local early-return) — we do **not** pay
  full error-return propagation across all callers.
- **Risk — `error_collection_append` currently auto-logs to the global logger.**
  Re-pointing it at the owning instance's sink must preserve ordering/output for the
  load-phase (`io`/`config`) collections so their byte output is unchanged.
- **Risk — XMI poll API change.** `library_fatal_raised()` becomes per-instance;
  `swap_ensemble_mod` and `swap_xmi_mod` callers update together in phase D.

## 10. Out of scope

- Rewriting physics, equations, FP ordering. Logging/error changes only.
- A structured (JSON) log format — plain text record with sim-time/instance-id is
  enough; revisit only if a downstream tool needs it.
- heat/solute coverage beyond what a real diagnostic need justifies.
