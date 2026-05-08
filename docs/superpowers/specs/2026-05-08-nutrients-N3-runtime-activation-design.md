# [nutrients] N3 — Runtime Activation Design

**Date:** 2026-05-08
**Status:** accepted (pending implementation)
**Sub-arc of:** `[nutrients]` umbrella (N1 → N2a → N2b → **N3**)
**Predecessors:** ADR 0025 (N1 crop-side adapter), ADR 0026 (N2a soil-side initial state), ADR 0027 (N2b timed amendments)
**Successor ADR:** 0028

## Goal

Open the `flCropNut` runtime gate on the TOML build. After this arc, a rotation with `cropwofost.nutrient.flcropnut = true` runs the WOFOST nutrient subsystem end-to-end without hitting any stub-error, using **only** typed-config inputs — no TTutil reader, no `Read_Nutrient`, no legacy `.crp` parsing.

## Out of scope

- Authoring a regression test case. No legacy nutrient-enabled fixture exists; deferred to a separate review later.
- Validating nutrient-physics correctness. Without a parity baseline this cannot be done in this arc.
- Adding new nutrient features.

## Context

After N1+N2a+N2b, all nutrient typed-config flows to the legacy globals:

- N1 (ADR 0025): per-rotation crop-side block via `apply_cropwofost_nutrient`, called from `cropwofost_init_from_config` only when `cfg%nutrient%flcropnut = .true.`.
- N2a (ADR 0026): top-level `[nutrients]` block via `apply_nutrients`, called unconditionally from `config_to_variables`.
- N2b (ADR 0027): timed amendments from CSV companion via `apply_nutrients_events`.

The legacy global `flCropNut` (declared in `src/core/variables.f90:563`) gates nutrient physics inside the daily loop:

- `src/core/swap.f90` lines 320, 327, 330, 337, 374, 411 — `SoilManagement` calls
- `src/crop/cropgrowth.f90` ~10 sites — nutrient stress, root/leaf nutrient pools, NPK demand/supply
- `src/crop/tillage.f90:73` — runtime stub-error (the closed gate)
- `src/crop/management_soil.f90:89` — comment only

`flCropNut` is currently **never assigned anywhere** in `src/`. It defaults to `.false.` and stays `.false.`. That is the gap N3 closes.

Two stub-errors hold the gate shut today:

- **Validator stub-error** at `src/crop/cropwofost_init.f90:85-87` — fatal-errors when `cfg%nutrient%flcropnut = .true.`. Inputs cannot even reach the runtime.
- **Runtime stub-error** at `src/crop/tillage.f90:73` — fatal-errors when `flCropNut = .true.` inside `DoTillage`.

## Decision

Drive `flCropNut` per-rotation from the typed config, inside `cropwofost_init_from_config`, and remove both stub-errors.

### Rationale for per-rotation wiring

`flCropNut` is a single module-level logical, but `flcropnut` is per-rotation in typed config. `cropwofost_init_from_config` runs at every rotation start, before that rotation's first daily step. All `flCropNut` read sites execute during a rotation. Therefore assigning `flCropNut = cfg%nutrient%flcropnut` at rotation init correctly handles a sequence with mixed `flcropnut` values across rotations. This mirrors the existing per-rotation crop config cache pattern (ADR 0016).

An OR-reduction across all rotations was rejected — it would enable nutrient gates for rotations whose config explicitly disabled them.

### Rationale for deferring the regression fixture

No legacy nutrient-enabled `.crp` fixture exists in `tests/swap-cases/`. To establish a parity baseline we would need to hand-author a legacy fixture, run the legacy binary against it, and capture `swap.ok`. That is a separate workstream. N3 ships the runtime activation alone; the fixture arc can follow once we are ready to validate physics output.

## Code changes

### 1. `src/crop/cropwofost_init.f90`

**Add** `flCropNut` to the `use variables, only:` list at the top of the module.

**Replace** the existing N1 block at line ~464-468:

```fortran
   ! [nutrients] N1: apply nutrient block to legacy globals when
   ! flcropnut=true on this rotation. Validator and the runtime
   ! gate at tillage.f90:73 still control whether nutrient physics
   ! actually runs; this just makes the typed-config inputs
   ! available if/when the gate is lifted (N3).
   if (cfg%nutrient%flcropnut) call apply_cropwofost_nutrient(cfg%nutrient)
```

with:

```fortran
   ! [nutrients] N3: drive the legacy global flCropNut from the
   ! per-rotation typed config. cropwofost_init_from_config runs at
   ! every rotation start, so a sequence of rotations with mixed
   ! flcropnut values toggles the gate correctly.
   flCropNut = cfg%nutrient%flcropnut
   if (flCropNut) call apply_cropwofost_nutrient(cfg%nutrient)
```

**Delete** the validator stub-error at line 85-87:

```fortran
   if (cfg%nutrient%flcropnut) &
      call fatalerr_collected('cropwofost_init', &
         'flcropnut=.true. not supported on TOML path; validator should have rejected.')
```

### 2. `src/crop/tillage.f90`

**Delete** the runtime stub-error at line 73:

```fortran
   if (flCropNut)        call fatalerr_collected ('DoTillage', 'flCropNut = 1 not (yet) allowed')
```

### 3. Stale comment cleanup

Delete or rephrase `flCropNut=1 is stub-errored` references at:

- `src/core/swap.f90:192`
- `src/crop/management_soil.f90:89`
- `src/crop/cropwofost_init.f90:21, 414`
- `src/crop/cropgrowth.f90:1139`

These are stale narrative comments; they have no functional effect but mislead readers about current state.

## Non-regression: TTutil-free

The N3 change adds **zero** new dependencies. `cropwofost_init.f90` is already TTutil-clean (modernized in N1). All nutrient inputs flow from typed config through N1/N2a/N2b adapters. There is no reader, no `Read_Nutrient`, no `dtutil`/`ttutil` shape call introduced or revived.

## Testing

No regression case authored in this arc. Verification consists of:

1. **`pixi run -e test test-pfunit`** — full pFUnit suite stays green (607 tests + 1 disabled). The change is structural; no new pFUnit test exercises nutrient physics without a fixture.
2. **`pixi run check-full`** — all 5 existing TOML cases stay byte-identical against `swap.ok` baselines. None sets `flcropnut=true`, so `flCropNut` stays `.false.` and observed behaviour is unchanged.

### Hand-verification (one-off, not committed)

To prove the runtime gate is actually open, locally edit one TOML case to set `flcropnut = true` plus the required `[wofost.nutrient]` block, run the binary, and confirm:

- Build does **not** abort with `'flcropnut=.true. not supported on TOML path'` (validator stub gone).
- Build does **not** abort with `'flCropNut = 1 not (yet) allowed'` (runtime stub gone).
- Simulation completes. Output may differ numerically from the baseline; that is expected (no parity baseline yet).

Revert the local edit before commit. This is a smoke-test, not a fixture.

## ADR 0028

`docs/adr/0028-nutrients-N3-runtime-activation.md` captures:

- **Status:** accepted
- **Context:** N1+N2a+N2b populated all nutrient typed-config → legacy globals; `flCropNut` runtime gate held closed by two stub-errors.
- **Decision:** drive `flCropNut` per-rotation from `cfg%nutrient%flcropnut` inside `cropwofost_init_from_config`; remove validator and runtime stub-errors; defer the regression fixture.
- **Consequences:** `[nutrients]` umbrella complete on the TOML path; nutrient physics callable end-to-end without TTutil; correctness validation deferred until a fixture exists.
- **References:** ADRs 0025/0026/0027 (predecessors), 0023/0024 (TTutil retirement context).
