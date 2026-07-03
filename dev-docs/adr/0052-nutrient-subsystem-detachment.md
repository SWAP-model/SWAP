# ADR 0052 — Nutrient subsystem detachment

**Status:** Accepted (2026-07-03)
**Implements:** ADR 0051 (SWAP-core vs coupled-component boundary), step 2 —
detach nutrients first.
**Supersedes the direction of:** ADR 0025 (N1 crop-side adapter), ADR 0026 (N2a
soil-side initial state), ADR 0027 (N2b timed amendments), ADR 0028 (N3 runtime
activation). Those arcs *reactivated* the nutrient subsystem on the modern
pipeline; this ADR reverses that reactivation and removes the subsystem.
**Precedent:** ADR 0040 (macropore retirement) — same profile (large, embedded,
unverified, blocking), same treatment (delete + preserve on the legacy branch +
re-implement-later path).

## Context

Per ADR 0051, SWAP's core is the vadose-zone water/heat/solute engine plus a
native table-driven "Plant"; detailed crop growth and nutrients are coupled
components, not embedded models. The nutrient subsystem is the first detachment
because it is the clearest case:

- **It is not SWAP's science.** The soil-N module (`Wofost_Soil_Declarations` /
  `Wofost_Soil_Interface` and the `wofost_soil_*` routines) is ANIMO's science
  embedded (the code cites *"according to ANIMO (Groenendijk et al, 2005)"*,
  *"Release of ANIMO4.0"*). The crop-side N is WOFOST-N / LINTUL4. Detailed
  nutrient fate belongs to the external **ANIMO** model, coupled via SWAP's
  hydrology (the `.afo` Aggregated Flux Output, parked by ADR 0009).
- **It is unverified and dark.** Gated behind `flCropNut`, which is **false in
  every regression case**; no case exercises it, and there is no byte-identical
  oracle (pristine SWAP 4.2.0 itself guards nutrients in the tillage path and
  ships no nutrient-enabled example — verified against `swap420gf`). The modern
  reactivation (ADRs 0025–0028) lifted that guard, going beyond stock 4.2.0.
- **It is the dominant remaining parallelism blocker.** ~85 module-global
  variables (`Wofost_Soil_Declarations` ~70 + `Wofost_Soil_Interface` 7 + the
  `cropwofost_init` `cw_*` nutrient subset ~7) — the largest obstacle to the
  threaded ensemble (T2-C), which a de-global (the abandoned T1-E-b) would have
  had to migrate with no oracle.

## Decision

**Remove the nutrient subsystem from the modern build**, in the ADR 0040 style:

1. Delete the nutrient-only modules (the `wofost_soil_*` cluster,
   `management_soil` / `SoilManagement`, `wofostnut`, `nutrients_config`,
   `nutrients_state`, `nutrients_csv`, `read_nutrients_toml`) and their build/test
   registrations. (~2,700 LoC in nutrient-only files, plus embedded crop-N.)
2. Excise the nutrient code embedded in shared files — the `flCropNut`-gated
   `SoilManagement` calls in `swap_mod`, the `wofost_apply_nstress` /
   `nutremrg` / `cw_*`-nutrient path in the WOFOST crop runtime and init, the
   `[nutrient]` reader block and `wofost_nutrient_t` config type, the top-level
   `[nutrients]` config wiring, and the `flCropNut` gate in `tillage`. Retain the
   non-nutrient crop config/state (the `cw_*` crop-growth constants that are not
   nutrient).
3. Remove `flCropNut` and the `[nutrients]` / `[cropwofost.nutrient]` TOML
   surface. A config that sets them becomes an error (unknown key) or is dropped —
   documented as a removed capability.

**Byte-identity basis.** Because every regression case has `flCropNut = .false.`,
the nutrient compute never runs on the tested path; excising the gated blocks and
removing the (never-true) gate is byte-identical. `check-fast` / `check-full`
must stay byte-identical at every step; `check-bindings` stays green.

**Preservation.** The modern reactivation work remains in git history (ADRs
0025–0028 and their commits); pristine 4.2.0's nutrient implementation remains on
`legacy/swap-4.2.0`. A future coupled-nutrient capability is re-enabling the
`.afo` hydrology output (ADR 0009's escape clause) and coupling to real ANIMO —
not resurrecting the embedded WSN.

## What this removes as a user-facing capability

- **WOFOST N-limited crop growth** (the `flcropnut` mode) and the in-SWAP soil-N
  simulation. Detailed nutrient work moves to external ANIMO (SWAP provides
  hydrology). Online N→crop-growth feedback is **not** preserved (offline ANIMO
  is one-way); reinstating it would be an online 3-way coupling — a separate,
  deliberate decision (ADR 0051 open question 1).

## Consequences

- SWAP sheds ~85 module globals and its largest embedded non-native subsystem;
  the threaded-ensemble blocker (T2-C) shrinks accordingly and T1-E-b (the
  nutrient de-global) is dropped.
- pFUnit count drops by the nutrient suites (`test_nutrients_config`,
  `test_apply_cropwofost_nutrient`, `test_nutrients_csv`); no regression case
  changes output.
- The WOFOST crop still runs — only its N-limited mode is gone; the simple crop
  and the water-side stress functions are untouched (ADR 0051 core).

## Manifest

From the exhaustive removal-surface discovery (verified against source).

**Files deleted (15, ~2,700 LoC of nutrient-only code):**
- `src/config/nutrients_config.f90`, `src/state/nutrients_state.f90`,
  `src/io/csv/nutrients_csv.f90`, `src/io/toml/read_nutrients_toml.f90`
- The WOFOST soil-N cluster (9): `wofost_soil_{declarations, interface,
  parameters, rateconstants, orgmatn, watern, amendments, cropresidues,
  balancecheck}.f90`
- `src/crop/wofost/wofostnut.f90` (legacy N subroutines), and
  `src/crop/wofost/management_soil.f90` (`SoilManagement`)

**Shared files surgically excised (11):** `swap_mod.f90` (7 `flCropNut`-gated
`SoilManagement` calls + 2 `use`), `cropgrowth.f90` (`wofost_apply_nstress`
dispatch), `cropwofost_runtime.f90` (the `wofost_apply_nstress` subroutine + all
`flCropNut`-gated N blocks + the `deaths()` N-stress branch and its signature),
`cropwofost_init.f90` (the `cw_*` nutrient vars + `apply_cropwofost_nutrient` +
the `flCropNut` write), `cropwofost_config.f90` (`wofost_nutrient_t` +
validate/finalize), `read_cropwofost_toml.f90` (`[nutrient]` block),
`crop_common_state.f90` (`flCropNut` field), `swap_config.f90` / `swap_state.f90`
(the `nutrients` sub-records), `load_swap_config.f90` (`[nutrients]` reader),
`tillage.f90` (the `flCropNut` gate).

**Build/test deregistration:** the nutrient sources from `meson.build` and
`tests/unit/meson.build`; 3 `ADD_TEST_SUITE` lines (`test_nutrients_config`,
`test_apply_cropwofost_nutrient`, `test_nutrients_csv`); the 3 `.pf` files + 5
nutrient CSV fixtures.

**Byte-identity:** all removed calls are `flCropNut`-gated and `flCropNut` is
`.false.` in every regression case; `deaths()`'s retained water-stress/LAI death
terms are unchanged. `check-fast`/`check-full` byte-identical; `check-bindings`
green.
</content>
