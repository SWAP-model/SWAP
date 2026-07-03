---
title: "Arc T1-E — crop-cluster reentrancy (spec + plan)"
date: 2026-07-03
status: in-progress
tags: [tier1, crop, save, multi-instance, thread-safety, byte-identical]
parent: 2026-07-03-modernization-review.md
---

# Arc T1-E — crop-cluster reentrancy

The keystone Tier-1 arc: remove the module-global mutable state / blanket `SAVE`
that makes the crop cluster unsafe to run concurrently across ensemble columns
(the blocker to `!$omp parallel do` in `swap_ensemble_mod`, and to per-instance
handles). Byte-identical throughout.

## Discovery (code-verified — corrects the review's scope)

1. **`crop_config_global` is already retired.** Every reference in `src/crop`
   and `src/driver` is a *comment/tombstone* (e.g. `swap_mod.f90:27`: "…
   crop_config_global retired"); there is no live declaration. **Drop from
   scope** — the review listed it as pending; it isn't.
2. **The crop *growth* path uses blanket `SAVE` on pure compute routines.**
   `cropgrowth.f90:595` (TOTASS) and `:672` (ASSIM) — WOFOST photosynthesis —
   plus `cropgrass_runtime.f90:93`. Their saved locals are (a) `data`-init Gauss
   constants (`xgauss`/`wgauss`/`scv`, implicitly `SAVE` regardless) and (b)
   working locals written-before-read each call. Removing the bare `SAVE` leaves
   the constants saved (via their `data` statements) and makes the working locals
   stack/thread-local — byte-identical (nothing is read before it is written) and
   reentrant. **This path is exercised by hupselbrook (maize/potato = WOFOST) and
   the grass cases, so `check-fast` verifies it.**
3. **The WOFOST *nutrient* cluster is a large module-global surface AND dark in
   the regression suite.** `Wofost_Soil_Interface` (13 `real(8)` module scalars)
   and `Wofost_Soil_Declarations` (~100 module vars: WSN organic-matter / N
   dynamics working state) are genuine mutable module globals. But the WSN
   compute is gated behind `flCropNut` (`cropgrowth.f90:424`), and **no
   regression case enables nutrients** (verified: no `[nutrients]`/`swcropnut` in
   any case TOML). So de-globaling them is byte-identical *by deadness* but
   **unverifiable** — and per the post-Phase-4 record there is no legacy
   nutrient-enabled `.crp` to build a fixture from.

## Scope split — T1-E is a multi-arc *program*, not one arc

Exhaustive discovery (7 `SAVE` sites, ~111 module-level mutable vars, ~95
carries-state) shows the crop cluster is far larger than "remove 3 saves". Only
one piece is a clean, verifiable de-save; the rest are genuine state migrations,
and one whole cluster is dark in the regression suite. **A cross-check caught the
discovery agent misclassifying `grass` as safe** — see the byte-identity note.

- **T1-E1 (this session — DONE): pure-compute de-save.** `TOTASS`
  (`cropgrowth.f90:595`) and `ASSIM` (`:672`) — WOFOST photosynthesis, pure
  functions of their arguments — keep only their `data`-init Gauss/scattering
  constants saved and drop the blanket `SAVE`. Byte-identical (`check-fast` 4/4);
  these routines are now reentrant.
- **T1-E2: WOFOST growth carries-state → `crop%wofost`.** `wofost()`
  (`cropwofost_runtime.f90:163`) blanket `SAVE` covers ≥10 verified cross-task
  carries-state locals: `vern`/`flvernalised` (vernalisation, read at `:434,437`),
  `gasst`/`gasstpot`, `mrest`/`mrestpot`, `tadw`/`tadwpot`, `fbl`/`drbl`/`drblpot`.
  Each must move to per-rotation `crop%wofost` state; the rest de-saved.
  Verifiable (WOFOST active in hupselbrook/swcf3). Careful, medium.
- **T1-E3: grass carries-state → `crop%grass`.** `grass()`
  (`cropgrass_runtime.f90:93`) — **discovery classified this write-then-read; it
  is not.** Verified carries-state: `flearlyhrvendact` set `.true.` at `:1156`,
  read at `:835`; `flearlyhrvendpot` similar. Move those to `crop%grass`; de-save
  the rest. Verifiable (grass active in grassgrowth/snow).
- **T1-E4: `cropfixed()` de-save** (`cropfixed_runtime.f90:47`). Agent says
  write-then-read — **must be re-verified per-local** (it made the same claim,
  wrongly, for grass) before de-saving.
- **T1-E5: `O2_pars%current_state` pointer** (`oxygenstress.f90:63`) — a
  module-level `state` pointer bound during `OxygenStress` for a ZBRENT callback.
  Thread-unsafe; thread the state through the callback instead. Verifiable
  (oxygenstress case active).
- **T1-E-b — SUPERSEDED (nutrients DELETED, not de-globalled).** Rather than
  migrate the ~85 nutrient module globals with no oracle, the whole WOFOST-N /
  ANIMO-derived soil-N subsystem was **removed** (ADR 0051 boundary + ADR 0052
  detachment). So `Wofost_Soil_Interface`, `Wofost_Soil_Declarations`, and the
  nutrient `cw_*` no longer exist. E-b is moot; the remaining E2–E5 (wofost/grass
  growth carries-state) stand as before. *(Original de-global analysis retained
  below for the record.)*

  **The "no coverage" blocker is resolvable — verified by construction
  (2026-07-03).** A nutrient-enabled case *can* be built from hupselbrook and it
  genuinely exercises the whole WSN cluster. Recipe + findings:

  - **Enable on a WOFOST rotation, not a fixed one.** The nutrient path is
    WOFOST-only: `flCropNut` is set solely at `cropwofost_init.f90:560` from
    `cfg%nutrient%flcropnut`. Crop `type` (`crop_config.f90`): 1 = fixed, 2 =
    WOFOST, 3 = grass. Hupselbrook's **maize is type 1 (fixed) — its `[nutrient]`
    block is ignored**; **potato is type 2 (WOFOST)** — put the block there. (First
    attempt on maize gave *zero* WSN activations; on potato, **579** — confirmed by
    instrumenting the `SoilManagement` dispatch at `swap_mod.f90:320`.)
  - **Inputs:** add `[nutrient] flcropnut=true` + the standard WOFOST-N params to
    `potatod.crp.toml` (`nmxlv` is an AFGEN table → DVS/value **pairs**, e.g.
    `[0.0,0.06,1.0,0.04,2.0,0.02]`, not bare values), and a top-level
    `[nutrients]` block with `sorp_coef` + initial pools (`fom(8)`, `bio`, `hum`,
    `cnh4`, `cno3`). Values from `tests/unit/io/toml/test_*nutrient*.pf`. Runs to
    completion, no NaN.
  - **Observability — the real subtlety.** The nutrient state does **not** change
    the water-balance output (transpiration follows the prescribed LAI/crop-factor
    tables, not N-limited biomass), and the output registry has **no soil-N
    variables at all**. But adding the crop-N residue columns `dwlvcrop`,
    `dwlvsoil` to the `[output.csv] inlist` **does** make it observable: nutrient-on
    vs -off differ in exactly the potato year in those columns (e.g. `293.09 /
    125.61` vs `0`). Those columns depend on N-uptake from the soil WSN, so they
    give an *indirect* characterization guard over the coupled crop↔soil-N state.

  **Two fixture routes, both now open:**
  1. **Characterization (modern golden-master) — unblocks E-b now.** Snapshot this
     case's output (default cols + `dwlvcrop`/`dwlvsoil`) as the reference; the
     de-globaling must reproduce it bit-for-bit. Sufficient for a *behavior-
     preserving* refactor. Needs harness support for a modern-vs-modern fixture
     (the suite currently only compares vs `swap420gf`), plus ideally dedicated
     soil-N output columns for a *direct* (not just indirect) guard.
  2. **`swap420gf` physics oracle — stronger, bigger.** Additionally build a
     matching legacy nutrient-enabled `.crp`/`.swp` so `swap420gf` produces a
     reference; validates the physics port, not just the refactor.

  **Recommendation:** do E-b as a characterization-guarded refactor (route 1),
  preceded by a small sub-task adding soil-N output variables to the registry so
  the ~70 `Wofost_Soil_Declarations` fields are *directly* observable.

The ~26 `cropwofost_init_mod` CONFIG-CONSTANT snapshot vars (`cw_rdrns`, …) are
read-only after init — safe to share across threads, no migration needed.

## Byte-identity analysis (T1-E)

`SAVE` + `-finit-local-zero` today: on the first call locals are zeroed, then
retained. Without `SAVE`: `-finit-local-zero` zeroes them every call. For a local
written before its first read in each call, the value at read time is identical
either way, so output is bit-for-bit unchanged. The only hazard is a
CARRIES-STATE local (read before written, relying on the previous call's value) —
discovery classifies every saved local for exactly this; any such local is kept
saved (or moved to `crop_state_t`) rather than dropped. `data`-init constants stay
saved via their `data` statements (implicit `SAVE`), so the Gauss weights are
untouched.

## Plan (each step ends `check-fast` green; verify each local from source first)

- **T1-E1 — DONE** (this session, one commit).
- **T1-E2..E5** — one focused arc each; per-local source verification of
  carries-state before touching a `SAVE`; migrate carries-state locals to the
  matching `crop%*` sub-record, de-save the rest; byte-identical `check-fast`
  after each. Order by value/confidence: E2 (wofost) → E3 (grass) → E5 (O2
  pointer) → E4 (cropfixed).
- **T1-E-b** — do not start until a nutrient regression fixture exists.

Only after E1–E5 land is the *growth* path reentrant; the threaded `!$omp` loop
(T2-C) additionally needs T1-E-b (nutrients) if coupled runs enable them, plus
per-instance handles (T2-B).

## Out of scope (this arc = T1-E1 only)
Everything above E1; the threaded loop (T2-C); config-constant consolidation
(T1-F).
