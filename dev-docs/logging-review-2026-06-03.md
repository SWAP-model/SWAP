# Logging, warnings & error-passing — infrastructure review and best-practice target

**Date:** 2026-06-03
**Branch:** `logging-review` (worktree off `development`)
**Scope:** review only — no code changed. Quantifies what exists, how consistently
it is used, and what a "good" target looks like for a scientific simulation model
that must run as a CLI *and* as an embedded library (pyswap / BMI / XMI).

---

## 1. What infrastructure exists

Three cooperating modules, all already sound in design:

### 1a. `src/core/swap_log.f90` — the leveled logger (290 L)
- Levels `DEBUG(10) / INFO(20) / WARN(30) / ERROR(40) / NONE(100)`.
- `log_init(level, file, to_stdout, timestamps)`, `log_set_level`, `log_close`.
- `log_debug/info/warn/error(context, message)` → `log_message`, which guards on
  `level < current_level`, formats `timestamp LEVEL context: message`, and writes
  to stdout and/or the log file (with `flush` after every file write).
- `to_str` generic for int/real/real64/logical so messages can interpolate values.
- `log_unit_handle()` escape hatch returning the raw file unit for legacy interop.
- **Module-level mutable state** (`current_level`, `log_unit`, `log_to_*`). This is
  acceptable for a logger but is hidden global state — worth noting against the
  repo's "no new hidden state" rule (it predates and is exempted, but it means the
  logger is process-global, not per-instance — see §3 multi-instance).

### 1b. `src/error/error.f90` — typed error accumulator (`error_mod`, 241 L)
- `error_t {code, message, context, is_fatal}` + stable integer code parameters
  (IO/PARSE/VALIDATION/FINALIZE/ADAPTER/LEGACY/DEPRECATED families).
- `error_collection_t` with `append / has_errors / has_fatals / count / summary /
  abort_if_fatal / clear`. `append` **auto-logs** through `log_error`.
- Threaded-`errors` pattern: fallible procedures take `errors` as `intent(inout)`
  and accumulate; a single edge call to `abort_if_fatal` turns fatals into
  `error stop`, writing `summary()` to `error_unit` first.
- `global_errors` singleton + `fatalerr_collected` / `FatalERR` drop-in for legacy
  physics paths that can't thread an argument.
- **`library_mode`**: when on, `abort_if_fatal` sets a sticky `fatal_was_raised`
  flag and *returns* instead of `error stop`, so an embedding host (Python/MODFLOW)
  isn't killed; the XMI layer polls `library_fatal_raised()`.

### 1c. `src/validation/validation.f90` — validator primitives (`validation_mod`, 101 L)
- Stateless `check_int_range / check_real_range / check_int_enum / check_not_empty /
  check_nonnegative_real / check_ordered_pair`, each appends to an `errors`
  collection on failure. Used by per-section config validators.

**Design verdict:** the architecture is good. Levels, typed errors, threaded
accumulation, an embedding-safe abort, and auto-logging on append are exactly the
right primitives. The problem is **adoption and wiring**, not design.

---

## 2. Consistency assessment (how much of the codebase actually uses it)

### 2a. The logger is barely turned on, and not at all when embedded
- `log_init` is called in **exactly one place**: `src/driver/swap_main.f90:17`,
  hardcoded `LOGLEVEL_INFO` + `swap_swap.log`.
- **BMI (`swap_bmi_mod`), C-API (`swap_capi_mod`), and ensemble/XMI
  (`swap_ensemble_mod`) never call `log_init`.** So an embedded run:
  - has **no log file** (`log_to_file` stays `.false.`),
  - still writes WARN/ERROR to **stdout** (module default `log_to_stdout=.true.`),
    polluting the Python/MODFLOW host's stdout,
  - runs at the hardcoded default `INFO`.
- **No user control of verbosity anywhere**: no CLI flag, no env var
  (`SWAP_LOG_LEVEL`), no TOML key. `log_set_level` exists but is never called.
  Changing the level requires editing source and recompiling.

### 2b. Logging coverage is thin and uneven
Call counts across all of `src/`:

| Level | Calls | Files | Notes |
|------:|------:|------:|-------|
| debug | ~31 | 1 real user | **All in `soilhydraulics.f90`**, guarded by `dump_convergence_diagnostics`. Unreachable at the default INFO level even when the flag is on (double gate). |
| info  | ~20 | 4 | init/lifecycle markers + itertime stats |
| warn  | ~28 | 9 | scattered |
| error | ~6  | 2 | almost only the indirect one inside `error%append` |

- Subsystems with **zero** logging: **heat, solute** (and effectively atmosphere,
  drainage, crop have 1–3 calls each).
- `simulation_config.f90` and `boundtop.f90` `use` a log routine they never call.
- Context strings are inconsistent: `swap`, `soilwater`, `Headcalc`, `BoundBottom`,
  `Calcgwl`, `Wlevbal`, `Astro`, `wofost` vs `CropGrowth_Wofost`, `file_io`,
  `itertime` — mixed case, mixed granularity (module vs subroutine vs phenomenon).

### 2c. Diagnostics still leak around the logger
- ~17 `write(*,...)` to stdout bypass the logger. Some legitimate (file headers in
  `io/file_headers.f90`, the `Swap normal completion!` banner, the log-open
  fallback). Several are not:
  - `src/timecontrol/timecontrol_mod.f90:146-149,515` — progress / day timeline to
    stdout (should be `log_info`).
  - `src/crop/oxygenstress.f90:1379,1418,1458,1566` — bare `write(*,*) 'QROMBD Too
    many steps.'`-style numeric-failure messages (should be `log_warn`/error).
  - `src/crop/cropgrowth_helpers.f90:417-425` — `write(777,...)`/`write(888,...)`
    explicitly marked `TEMPORARY OUTPUT DELETE`; orphaned debug to magic units.
  - `src/soilwater/dormant/sptabulated.f90` — `if (IER/=0) write(124,*) IER` to
    unopened magic units (silent no-op).
  - `src/atmosphere/interception.f90:259` — old-style bare `stop` after a
    `write(*,9199)`.

### 2d. Error/abort: modern at the edges, legacy in the core
- **Threaded `errors` pattern is essentially 100% confined to `io/`, `config/`,
  `validation/`** (~177 procedures, ~189 `append` sites). It has **not** penetrated
  any compute subsystem.
- **`fatalerr_collected` is the dominant abort path: ~150 call sites**, ~**99 of
  them (≈64%) inside physics/compute kernels** (crop 64, soilwater 21, drainage 13,
  atmosphere 6, core 6, solute 3, heat 1).
- `abort_if_fatal` is invoked at 16 edge sites, **all in init/load phase** (driver +
  state CSV loaders), none in the time-stepping loop. Good.

### 2e. The sharp consequence (the finding that matters most)
`fatalerr_collected` → `global_errors%abort_if_fatal()`, whose behaviour now depends
on a **process-global mode flag that the embedding paths set inconsistently**:

- **CLI**: a kernel `fatalerr` → `error stop`. Correct.
- **XMI/ensemble** (sets `library_mode=.true.`): a kernel `fatalerr` **returns
  normally** and the kernel **keeps computing with invalid state** until the host
  later polls `library_fatal_raised()`. Continue-with-garbage.
- **Direct BMI / C-API** (never set `library_mode`): a kernel `fatalerr` **`error
  stop`s and kills the host process** (the Python interpreter).

So the same fatal condition has three different behaviours depending on entry point,
two of them wrong. This directly undercuts the stated forward goal of "clean
BMI/Python binding + in-process multi-instance." It is the highest-value thing to
fix.

---

## 3. Multi-instance / embedding caveat (forward intent)

Both the logger (`current_level`, `log_unit`) and the error system (`global_errors`,
`library_mode`, `fatal_was_raised`) are **process-global**. The moment two SWAP
instances run in one process (the stated pyswap/coupling goal), they share one log
file, one level, and one fatal flag. Interleaved logs become unattributable and a
fatal in instance A trips the flag instance B polls. Any logging rework should
decide now whether diagnostics are process-global (simplest, tag every line with an
instance id) or threaded on `state` (cleanest, larger change).

---

## 4. Best practices — what good logging looks like for a model like this

A hydrological/scientific solver embedded as a library has specific needs that
differ from a web service. Target guidance:

### 4a. Level semantics (be disciplined about granularity)
- **ERROR** — the run (or this instance/timestep) cannot continue or produced an
  invalid result. Always emitted, always to file *and* stderr. Examples: input file
  missing/malformed, unsupported switch combination, solver produced NaN/Inf, mass
  balance closure violated beyond tolerance, groundwater below profile.
- **WARN** — physically/numerically suspicious but the run continues with a
  documented fallback. Examples: Richards solver hit max iterations but accepted the
  last iterate; switched to free drainage because the bottom compartment went oven
  dry; clipped a parameter to its valid range; a deprecated TOML key was ignored.
  Each WARN should say **what happened, where (date + node/compartment), and what
  was done about it**.
- **INFO** — lifecycle and provenance at human scale: run start with resolved
  config summary, project name, simulation window, key switches, output paths; phase
  boundaries (init complete, per-year rollovers if cheap); run end with status +
  wall time. INFO should be **O(1) per run / per phase, not per timestep** — a 50-yr
  daily run has ~18k steps; per-step INFO is noise.
- **DEBUG** — developer-facing, per-step / per-iteration detail, off by default:
  solver iteration counts and residuals, convergence dumps, intermediate fluxes,
  the actual numbers behind a WARN. Must be cheap to compile out / guard so it never
  costs anything in production runs.
- (Optional **TRACE** below DEBUG if you want per-node spew separated from
  per-step DEBUG.)

The repo's current taxonomy already matches this; the gap is that almost everything
is WARN-or-nothing, DEBUG is unreachable, and there's no per-step vs per-run
discipline.

### 4b. What every record should carry
- **Level, context, message** (already present).
- **Simulation time**, not just wall-clock timestamp — the model date / `daynr` /
  timestep index. For a solver, "diverged on 1998-07-14 step 3 node 12" is worth
  ten "diverged" lines. Wall-clock timestamps are near-useless for a fast batch run
  and should be optional/off by default; simulation time should be on.
- **Instance id** once multi-instance lands (see §3).
- For numeric warnings: the **offending value, the bound/threshold, and the
  action taken** — structured enough to grep.

### 4c. Terminal vs file — route by audience
- **Terminal (stdout/stderr):** for the CLI user, keep it quiet and meaningful —
  WARN+ERROR to **stderr**, a one-line start and one-line end to stdout, and a
  progress indicator only if interactive. INFO chatter and all DEBUG do **not**
  belong on the terminal by default. Errors to **stderr** (currently WARN/INFO/ERROR
  all go to stdout via one `write(*,...)`), so a pipeline can separate diagnostics
  from real output.
- **File (`*.log`):** the full record at the configured level (default INFO,
  DEBUG when asked). This is the forensic artifact for a 50-year batch run. It should
  always be written for a CLI run and **be opt-in with a caller-supplied path/level
  when embedded** (so pyswap can give each instance its own file or a callback).
- **Embedded:** never write to the host's stdout by default. Provide either a
  caller-set log file per instance or a registered callback/buffer the host drains —
  and route fatals to the existing sticky-flag mechanism, uniformly across BMI / C-API
  / XMI.

### 4d. Operational rules
- **One way to emit.** No `print *` / `write(*,*)` diagnostics; everything through
  `log_*`. Reserve raw writes for actual output files (`io/`) and the deliberate CLI
  banner.
- **User-settable verbosity** without recompiling: env var `SWAP_LOG_LEVEL` +
  optional `[logging]` TOML block (`level`, `file`, `to_stdout`, `timestamps`) +
  a C-API/BMI setter so pyswap can configure it. Initialize logging in **every**
  entry point, not just the CLI.
- **Never let logging change numerics.** Logging is stderr/file only and must not
  touch the compared output files, FP op-ordering, or kernel control flow — safe
  under the byte-identical regression law, but new DEBUG in hot loops must stay
  guard-gated and side-effect-free.
- **Fatal-in-kernel policy must be uniform.** Decide one behaviour for a compute-path
  fatal across all entry points (recommended: accumulate + set sticky flag + stop
  stepping that instance; never silently continue, never hard-kill the host).
  Long-term, migrate the highest-traffic kernel `fatalerr` sites (crop, soilwater) to
  return-an-error so a bad instance fails cleanly instead of aborting the world.

---

## 5. Suggested improvement arc (for discussion, not yet planned)

Ordered by value/See §2e, §2a first:

1. **Unify entry-point wiring + embedding safety** — `log_init` in BMI/C-API/XMI;
   route stdout default off when embedded; set `library_mode` consistently; one
   fatal-in-kernel policy. (Highest value: makes pyswap/coupling correct.)
2. **User-facing verbosity control** — `SWAP_LOG_LEVEL` env + `[logging]` TOML +
   C-API setter; stop hardcoding INFO/`swap_swap.log`.
3. **Route stray diagnostics through the logger** — kill the `TEMPORARY DELETE`
   debug, the bare `stop`, magic-unit writes; move timecontrol progress and
   oxygenstress numeric-failure messages to `log_*`; split WARN+ERROR to stderr.
4. **Add simulation-time + instance-id to records**; make wall-clock timestamps
   optional/off.
5. **Raise coverage where it pays** — solver convergence (DEBUG), boundary/switch
   fallbacks (WARN), per-run provenance summary (INFO). Leave heat/solute thin
   unless a real diagnostic need exists.
6. **(Long horizon) migrate hot kernel `fatalerr` sites to threaded errors** so a
   compute-path failure degrades one instance cleanly.

Each step is independently shippable and byte-identical (logging touches only
stderr/file, never the compared outputs).
