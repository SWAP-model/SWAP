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
- **T1-E-b (blocked): WOFOST nutrient de-globaling.** `Wofost_Soil_Interface`
  (7) + `Wofost_Soil_Declarations` (~70, blanket `save`) + `cropwofost_init_mod`
  carries-state subset (`cw_anlv`/`cw_anst`/`cw_nni`/`cw_fstr`/`cw_nmaxlv`/
  `cw_nmaxst`/`cw_nmaxrt`) → `state%nutrients`. **Gated on nutrient test
  coverage**, itself blocked (no oracle `.crp`; gated behind `flCropNut`, false in
  every regression case). Do NOT attempt this ~85-var migration with no oracle;
  the precursor is a physics-validated nutrient fixture (likely an ADR).

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
