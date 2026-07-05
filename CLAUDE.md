# CLAUDE.md — SWAP modernization

SWAP (Soil–Water–Atmosphere–Plant) is a legacy Fortran hydrological model
being modernized in place. This file is the operating contract for anyone —
human or AI — working in this repo. Read it before making changes.

## Where the modernization stands

The **strangler-fig migration is complete**: input flows
`TOML → typed config (swap_config_t) → typed state (swap_state_t)` with **no
bare-globals layer**. The old `src/core/variables.f90` / `arrays.fi` bare-global
module and the `swapoutput.f90` C-API output machinery have been **deleted**.
Every runtime symbol now lives on a mutable `state%X` record or a read-only
`config%X` record.

The remaining work is *clarity and capability*: decompose oversized
subroutines, finish config-constant consolidation, and move toward
**in-process multi-instance execution + a clean BMI/Python (pyswap) binding**.
That forward intent — stateless kernels, no hidden state, I/O at the edges —
is what every change should push toward. See
`dev-docs/post-phase-4-modernization-summary.md` for the full arc.

## Architecture (current)

- **Config** — `src/config/swap_config.f90` defines `swap_config_t`, the
  read-only typed result of parsing TOML (+ CSV companion files). Compute code
  receives config as `intent(in)`; it never reads files.
- **State** — `src/state/` holds one typed record per subsystem
  (`*_state_t`: atmosphere, soilwater, solute, heat, drainage, crop*, …),
  aggregated under `swap_state_t` in `src/state/swap_state.f90`. State is
  threaded by argument and mutated via `associate`/direct field writes.
- **Two-phase init** — construction is a type-bound `state%X%init(config%Y)`
  plus a free `x_seed` procedure; there is no integer `case(task)` lifecycle
  dispatch (ADR 0043).
- **I/O at the edges** — readers/writers live in `src/io/`; numerically hot
  kernels are free of file I/O.
- **Source layout (ADR 0048)** — *feature-first for behavior, layer-first for
  the shared data model.* `src/core/` holds only foundation leaves (`constants`,
  `arrays`, `dtutil`, `swap_log`, `arrayutils`, `numericalsolvers`);
  orchestration lives in `src/driver/` (`swap_mod`, `swap_main`,
  `swap_ensemble_mod`) and the C-ABI facades in `src/bindings/` (`swap_capi/bmi/
  xmi`, `bmi_constants`). Compute subsystems are feature folders:
  `atmosphere/`, `soilwater/` (incl. boundary conditions), `heat/`, `solute/`,
  `drainage/`, `timecontrol/`, and `crop/` with `fixed/` `grass/` `wofost/`
  sub-packages (all WOFOST nutrient dynamics under `crop/wofost/`). The shared
  data + edge layers — `config/`, `state/`, `io/`, plus `error/`, `validation/`
  — stay as their own folders (they are depended on by every feature). There is
  no `utils/` or `boundary/` folder.

## Non-negotiables

- **gfortran only.** No Intel/ifx code paths (ADR 0001).
- **Byte-identical regression is the law.** The oracle is `swap420gf`
  (SWAP 4.2.0 recompiled with the modern toolchain); the harness compares
  against `tests/regression/*_expected_gfortran.json`. A refactor must
  reproduce reference output bit-for-bit.
- **Do not touch physics.** Don't change equations, coefficients, units, or
  floating-point operation ordering in compute kernels. Removing `SAVE`/state
  and improving interfaces is fine *only* if results are unchanged. Prefer
  wrapping legacy kernels over rewriting them.
- **No new hidden state.** No new `SAVE`, `COMMON`, `DATA`, or mutable
  module globals. `implicit none` everywhere; `iso_fortran_env` kinds
  (`real64`/`int32`), not `REAL*8`.
- **Single `builddir/`** with pFUnit enabled (ADR 0002).

## Build & test

```bash
pixi run build-linux          # meson compile (auto-configures builddir)
pixi run -e test check-fast   # pFUnit + 4 regression cases — GREEN before EVERY commit
pixi run -e test check-full   # pFUnit + all regression cases — at milestones / phase tags
pixi run -e test test-pfunit  # unit suite only
pixi run lint                 # fprettify (lint-check for a dry-run diff)
pixi run clean                # rm -rf builddir  (see clean-rebuild rule below)
```

- **Clean rebuild is mandatory after any `src/state/*.f90` or
  `src/config/*` schema change** — incremental Meson does not propagate `.mod`
  deps across the swap_modern↔swap_legacy boundary. Do `rm -rf builddir` (or
  `pixi run clean`) then rebuild, or you get stale-`.mod` SIGSEGVs.
- **The regression suite is pytest** (ADR 0054). `check-fast`/`check-full`/
  `regression` all invoke `pytest tests/regression` (xdist-parallel); a single
  case is `pytest tests/regression -k <case>`, the fast subset is `-m fast`.
  The case *inputs* live in the `swap-testcases` sibling repo (resolved via
  `SWAP_TESTCASES_PATH` / `TESTCASES_REF`); the expected-output fixtures stay
  in-repo.
- **Fixtures:** the live comparison is against `*_reference_gf.json` — the
  gfortran-compiled SWAP 4.2.0 oracle (`tests/reference/swap420gf`). Regenerate
  only for a *legitimate* physics change, with
  `python tests/regression/regen_reference.py [<case>]`, and document it in the
  commit. (`*_expected_gfortran.json` is a *diagnostic* snapshot of the modern
  build's own output via `regen_expected.py` — NOT the reference. Never delete
  the historical `*_expected.json` ifx reference.) Accepted divergences
  (winter frost-path; macropore case excluded since ADR 0011) are baked in as
  xfails — don't re-derive them.
- **New pFUnit `.pf` files must be registered in BOTH `tests/unit/meson.build`
  (the `pf_files` list) and `tests/unit/testSuites.inc`** or they compile but
  never run. Confirm the `OK (N tests)` count rises.

## Working conventions

- **Small, cohesive changes.** One subsystem / one concern per commit. Never
  mix refactoring with behavior changes. Decompose umbrella tasks.
- **Per-task regression gate.** Every change ends with `check-fast` green —
  not a deferred end-of-arc check (that makes regressions impossible to
  bisect). pFUnit alone misses global-default regressions; run `check-fast`.
- **Verify state, not reports.** Confirm the actual code/output changed; a
  subagent's summary is not evidence.
- **Decompose into narrative helpers**, not one-line micro-functions:
  orchestrator + comprehensive phase-helpers; keep control flags in the
  orchestrator.
- **State array fields default to `allocatable`** (fixed-size only if
  < ~64 KB) — large fixed arrays land on the stack via `swap_state_t` and
  SIGSEGV on smaller stack limits.
- **Big tabular inputs → CSV companion files**, not inlined TOML/arrays.
- For broad multi-file sweeps, prefer subagent-driven execution; grep the
  field shapes during pre-flight before migrating a symbol.

## Commits & branches

- Work lands on **`development`** (directly, or on a short-lived
  worktree branch merged into it). `main` advances only at milestones.
  **Never push to `origin/main`.**
- Conventional-commit prefixes (`refactor(state):`, `fix(nutrients):`,
  `docs:`, `chore:`). Explain the *why* in the body.
- **pyswap is a separate repository, not a submodule.** It lives standalone at
  `~/code/pyswap` (`github.com/zawadzkim/pySWAP`); SWAP depends on it only as the
  PyPI package (`pyswap = ">=0.3.9"` in `pixi.toml`). The former `pyswap/`
  submodule and its redundant `release/v1` branch (identical to `main`, frozen at
  the v0.3.9 release) were removed 2026-06-27 — nothing in SWAP's build or tests
  consumed the checkout. Do pyswap work in its own repo off `main`; never push to
  a branch you weren't asked to.

## Where to look

- `dev-docs/adr/` (start at `index.md`) — **why** each non-obvious choice was
  made. ADRs are the durable record; a reversed decision gets a *superseding*
  ADR, never an edit.
- `dev-docs/phase-4-modernization-summary.md` +
  `dev-docs/post-phase-4-modernization-summary.md` — **what** was done, arc by
  arc.
- `dev-docs/DEVELOPMENT_GUIDE.md` — the **deep** style/design reference: modern
  Fortran rules, in-memory/BMI/Python design, multicore/GPU patterns, and the
  full code-review checklist. CLAUDE.md is the summary; that guide is the depth.
- `tests/regression/INVESTIGATION_NOTES.md` — open fixture-divergence questions.

## When in doubt

Stop and ask. A partial, honest report beats a confident-but-wrong commit.
