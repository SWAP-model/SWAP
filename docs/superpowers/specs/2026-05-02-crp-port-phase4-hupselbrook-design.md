# `.crp` port — Phase 4 (hupselbrook integration + fallback teardown) — design

**Status:** Draft, awaiting review
**Date:** 2026-05-02
**Predecessors:**
- `docs/superpowers/specs/2026-05-02-crp-port-phase1-cropfixed-design.md` (cropfixed via case 6)
- `docs/superpowers/specs/2026-05-02-crp-port-phase2-cropwofost-design.md` (cropwofost via case 5; authored in parallel)
- `docs/superpowers/specs/2026-05-02-crp-port-phase3-cropgrass-design.md` (cropgrass via cases 4+2; authored in parallel)
- `docs/adr/0015-strangler-narrow-scope-stub-errors.md`
- `docs/adr/0016-per-rotation-crop-config-cache.md`

---

## Goal

Land the case 1 (hupselbrook) integration test for the `.crp` port, then remove
the transitional `if/else` dispatch fallback introduced by Phases 1–3.

After Phase 4:

- `tests/swap-cases/toml/1.hupselbrook/` contains three fully-authored `.crp.toml`
  files (`maizes.crp.toml`, `potatod.crp.toml`, `grassd.crp.toml`) with 1:1 schema
  coverage of the legacy `.crp` sources.
- Case 1 regression runs green with no legacy `.crp` files present in the TOML case
  directory (all three are deleted from the TOML dir after the smoke test).
- `cropgrowth.f90` no longer contains `readcropfixed`, `readwofost`, or `readgrass`
  calls on the primary (TOML) execution path. The three `if/else` dispatch blocks are
  replaced with unconditional `*_init_from_config` calls.
- The legacy `readcropfixed` / `readwofost` / `readgrass` subroutines remain alive in
  `readswap.f90` as parity-test fixtures (per ADR 0015).
- All 5/5 regression cases remain green.

---

## Non-goals (Phase 4)

- Schema / parser / init-module work for any new fields. All schema, parser, and init
  work is completed by Phases 1–3. Phase 4 authors TOML, wires dispatch, and tears down
  the transitional fallback.
- Removing `readcropfixed` / `readwofost` / `readgrass` from `readswap.f90` — kept as
  parity-test fixtures. Their further removal is a future cleanup task.
- Removing the `crop_config_global` module-level pointer — deferred to the
  config-passing follow-on spec (ADR 0016 Part C). See Teardown plan below.
- Removing the `rotation_loaded(:)` sentinel array on `crop_config_t` — same trigger as
  `crop_config_global`.
- The config-passing direction (ADR 0016 Part C) — out of scope. Phase 4 ends with the
  legacy fallback removed and the per-rotation cache as the sole source of crop runtime
  data. The pointer-based access pattern remains until the follow-on spec.

---

## Status quo (depends on Phases 1–3 landing)

Phase 4 cannot begin until all three predecessor phases are merged and the regression
suite is 5/5 green. Specifically:

| Precondition | Provided by |
|---|---|
| `cropfixed_config_t` schema 1:1; `read_cropfixed_toml`; `cropfixed_init_from_config`; type-1 dispatch in `cropgrowth.f90:397` guarded by `rotation_loaded(icrop)` | Phase 1 |
| `cropwofost_config_t` schema 1:1; `read_cropwofost_toml`; `cropwofost_init_from_config`; type-2 dispatch in `cropgrowth.f90:984` guarded by `rotation_loaded(icrop)` | Phase 2 |
| `cropgrass_config_t` schema 1:1; `read_cropgrass_toml`; `cropgrass_init_from_config`; type-3 dispatch in `cropgrowth.f90:2068` guarded by `rotation_loaded(icrop)` | Phase 3 |
| `crop_config_global` module pointer set by `config_to_variables`; `rotation_loaded(:)` on `crop_config_t` | Phase 1 |
| Case 6 TOML dir has no `.crp` file; case 5 TOML dir has no `.crp` file; cases 4+2 TOML dirs have no `.crp` files | Phases 1–3 |

After Phases 1–3, the dispatch blocks in `cropgrowth.f90` look like this (one per crop
type, at the three task=1 callsites):

```fortran
! Example: type-1 (cropfixed) dispatch after Phase 1
if (associated(crop_config_global) .and. &
    allocated(crop_config_global%rotation_loaded) .and. &
    crop_config_global%rotation_loaded(icrop)) then
   call cropfixed_init_from_config(crop_config_global%rotation_fixed(icrop), icrop)
else
   call readcropfixed (icrop, cropfil(icrop), lcc, swhydrlift)   ! transitional
end if
```

Phase 4 removes the `else` branch from each of the three callsites.

---

## Case 1 rotation audit

`tests/swap-cases/toml/1.hupselbrook/swap.toml` defines exactly **3 rotations**
(not 5 — the brief's "5 rotations" was approximate):

| # | Period | File | Type | Crop | Active switches |
|---|---|---|---|---|---|
| 1 | 2002-05-01 → 2002-10-15 | `maizes.crp.toml` | 1 (cropfixed) | Maize | `SWCF=2`, `SWOXYGEN=1`, `SWRD=1`, `SWDROUGHT=1`, `SWINTER=1`, `SWGERM=0`, `SWHARV=0` |
| 2 | 2003-05-10 → 2003-09-29 | `potatod.crp.toml` | 2 (cropwofost) | Potato | `SWCF=2`, `SWGERM=2`, `IDSL=0`, `SWOXYGEN=1`, `SWWRTNONOX=1`, `SWRD=2`, `SWDROUGHT=1`, `SWCO2=0`, `SWINTER=1` |
| 3 | 2004-01-01 → 2004-12-31 | `grassd.crp.toml` | 3 (cropgrass) | Grass | `SWCF=2`, `SWRD=3` (`RLWTB`-biomass), `SWRDC=0`, `SWOXYGEN=1`, `SWJARVIS=4`, `SWDROUGHT=1`, `SWCO2=0`, `SWHARVEST=1+SWDMMOW=2` |

### Schema gap audit (Phase 4 scope check)

The audit's purpose is defensive: checking whether case 1 activates any switch
combination not exercised by cases 6 / 5 / 4+2, which would require retroactive schema
work. Findings:

**Type-1 (maizes.crp) vs. case-6 grass.crp:**
- `SWCF=2` (crop height + albedo/rsc/rsw) — case 6 uses `SWCF=1`. This is a material
  difference: the TOML file must author the `chtb` table and the `ALBEDO`/`RSC`/`RSW`
  fields in the `[crop_factor]` section. Phase 1's schema 1:1 coverage includes
  `swcf=2` fields; Phase 1's validator supports `swcf ∈ {1, 2}` (only `swcf=3` is
  stub-errored). **No schema gap** — Phase 1 covers it.
- `SWOXYGEN=1` with `HLIM1`/`HLIM2U`/`HLIM2L` — case 6 uses `SWOXYGEN=0`. Phase 1
  schema and init module cover `swoxygen=1` (Feddes). **No schema gap.**
- `GCTB` has 7 pairs (DVS 0.0–2.0, LAI 0.05–5.20) — larger than case 6 but same
  structure. **No gap.**
- `CHTB` is used (7 pairs, DVS/CH) alongside `CFTB` (7 pairs, DVS/CF) — both stored in
  case 6's schema. **No gap.**

**Type-2 (potatod.crp) vs. case-5:**
- `SWGERM=2` (temperature + hydrology) — Phase 2's scope; see Phase 2 spec for coverage.
  Assumption: Phase 2 covers `swgerm ∈ {0,1,2}` for case 5.
- `SWWRTNONOX=1` — Phase 2 scope; covered.
- Potato `potatod.crp.toml` is already 268 lines in the TOML dir with all sections
  present. It appears to be near-complete from a prior preparation task.

**Type-3 (grassd.crp) vs. cases 4+2:**
- `SWJARVIS=4` — this is the legacy `readgrass`-specific compensation switch (distinct
  from the `SWCOMPENSATE` field in cropfixed/cropwofost). If Phase 3 covers cases 4+2's
  grass crop, which use `SWJARVIS`, then this is not a gap. **Risk:** If Phase 3 only
  covers `SWJARVIS=0`, case 1's `SWJARVIS=4` would be a schema gap requiring retroactive
  Phase 3 work. Document as a risk; the implementing agent must verify Phase 3's
  `cropgrass_config_t` covers `SWJARVIS=4`.
- `SWRD=3` (biomass-based root growth via `RLWTB`) — Phase 3 scope. The grassd.crp
  already uses `SWRD=3` in both case 4+2 and case 1.
- `SWRDC=0` — covered.
- `SEQGRAZMOW` all-mowing (`2 2 2 ...` × 20) with `SWDMMOW=2` (flexible DM table
  `DMMOWTB`) — Phase 3 scope. `grassd.crp.toml` stub already extended in a prior
  task (Phase 4d Task 18) with a `[mowing]` section.
- `SWCO2=0` — Phase 3 scope.

**Summary:** No new schema gaps were found that Phase 4 must address itself. One
potential risk (SWJARVIS=4 coverage in Phase 3) is flagged in the Risk Register.

---

## Approach

Phase 4 has three workstreams executed sequentially:

### Workstream A — Author the three `.crp.toml` files

Translate the three legacy `.crp` files in
`tests/swap-cases/1.hupselbrook/` to fully-authored `.crp.toml` counterparts in
`tests/swap-cases/toml/1.hupselbrook/`. Each is a submodule pair commit.

1. **`maizes.crp.toml`** (type 1, cropfixed): expand the current skeletal 30-line file
   to full 1:1 coverage. The skeleton already has phenology, light, root (partial),
   water_stress, and interception. Missing: `[preparation]`, `[harvest]`, `[lai]` (with
   `swgc=1` + `gctb` 7-pair table), `[crop_factor]` (with `swcf=2` + `chtb` 7-pair
   table, `albedo`, `rsc`, `rsw`), `[root]` (with `swrd=1` + `rdtb` 6-pair table +
   `rdctb` 2-pair table), `[oxygen_stress]` (with `swoxygen=1` + `hlim1/hlim2u/hlim2l`),
   `[drought_stress]`, `[salinity]`, `[compensate]`, `[scheduling]`.

2. **`potatod.crp.toml`** (type 2, cropwofost): the existing 268-line file appears
   complete. The implementing agent must verify field-by-field against `potatod.crp`
   (704 lines) to confirm all sections are present. Likely only minor gaps or formatting
   issues remain.

3. **`grassd.crp.toml`** (type 3, cropgrass): the current 41-line stub has `[light]`,
   `[root]`, `[water_stress]`, `[interception]`, `[mowing]`, `[grazing]`. Missing from
   the full `grassd.crp`: `[initial]` (`tdwi`, `laiem`, `rgrlai`, `swtsum`), `[green_area]`
   (`ssa`, `span`, `tbase`, `slatb`), `[assimilation]` (`kdif`, `kdir`, `eff`, `amaxtb`,
   `tmpftb`, `tmnftb`), `[conversion]` (`cvl`, `cvr`, `cvs`), `[respiration]` (`q10`,
   `rml`, `rmr`, `rms`, `rfsetb`), `[partitioning]` (`frtb`, `fltb`, `fstb`), `[death]`
   (`perdl`, `rdrrtb`, `rdrstb`), `[root]` (expanded with `swrd=3`, `rlwtb`, `wrtmax`,
   `swrdc`, `rdctb`), `[oxygen_stress]`, `[drought_stress]`, `[salinity]`, `[compensate]`,
   `[co2]`, `[management]`.

### Workstream B — Regression test case 1 with full TOML

Run the regression with all three `.crp.toml` files authored and the `.crp` files still
present. Verify all 5 cases pass. If case 1 mismatches, diagnose and fix before
proceeding to the dispatch removal.

### Workstream C — Remove the legacy fallback dispatch (teardown)

Remove the `if/else` guard around each of the three legacy reader calls in
`cropgrowth.f90`. The unconditional replacement is described in the Teardown Plan
section below.

Then: smoke test by temporarily renaming all three `.crp` files in the TOML dir and
re-running the regression. On green: actually delete the three files and commit.

---

## Teardown plan

This section specifies precisely what gets removed, what the post-removal code looks
like, and what remains.

### The three dispatch removals

#### Type-1 dispatch (`cropgrowth.f90` — `subroutine cropfixed`, task=1)

**Before (after Phase 1):**
```fortran
      if (associated(crop_config_global) .and. &
          allocated(crop_config_global%rotation_loaded) .and. &
          crop_config_global%rotation_loaded(icrop)) then
         call cropfixed_init_from_config(crop_config_global%rotation_fixed(icrop), icrop)
      else
         call readcropfixed (icrop, cropfil(icrop), lcc, swhydrlift)   ! transitional
      end if
```

**After (Phase 4):**
```fortran
      call cropfixed_init_from_config(crop_config_global%rotation_fixed(icrop), icrop)
```

The `associated(...)` / `allocated(...)` checks become unnecessary because, after Phase
4, every type-1 rotation in every TOML case has a populated `rotation_fixed(icrop)` slot
(the loader is called unconditionally at config-load time). The call is therefore
unconditional. If `crop_config_global` is somehow null (e.g. a test harness that bypasses
`config_to_variables`), the null-dereference is a programming error, not a handled
runtime path; it surfaces immediately in testing.

#### Type-2 dispatch (`cropgrowth.f90` — `subroutine wofost`, task=1)

**Before (after Phase 2):**
```fortran
      if (associated(crop_config_global) .and. &
          allocated(crop_config_global%rotation_loaded) .and. &
          crop_config_global%rotation_loaded(icrop)) then
         call cropwofost_init_from_config(crop_config_global%rotation_wofost(icrop), icrop, &
                                          swhydrlift, swsoybean, mg, dvsi, dvrmax1, dvrmax2, &
                                          flrfphotoveg, tmaxdvr, tmindvr, toptdvr, popt, pcrt, &
                                          flphenodayl, FraDeceasedLvToSoil)
      else
         call readwofost (icrop, cropfil(icrop), swhydrlift, swsoybean, mg, dvsi, dvrmax1, dvrmax2, &
                          flrfphotoveg, tmaxdvr, tmindvr, toptdvr, popt, pcrt, flphenodayl, FraDeceasedLvToSoil)
      end if
```

**After (Phase 4):**
```fortran
      call cropwofost_init_from_config(crop_config_global%rotation_wofost(icrop), icrop, &
                                       swhydrlift, swsoybean, mg, dvsi, dvrmax1, dvrmax2, &
                                       flrfphotoveg, tmaxdvr, tmindvr, toptdvr, popt, pcrt, &
                                       flphenodayl, FraDeceasedLvToSoil)
```

**Note on wofost signature:** The exact argument list of `cropwofost_init_from_config`
may differ from the above — it mirrors whatever Phase 2 defined. The implementing agent
must read Phase 2's `src/crop/cropwofost_init.f90` to confirm the signature. The
`readwofost` legacy call (at `cropgrowth.f90:984`) provides the ground truth for which
local variables the call uses.

#### Type-3 dispatch (`cropgrowth.f90` — `subroutine grass`, task=1)

**Before (after Phase 3):**
```fortran
      if (associated(crop_config_global) .and. &
          allocated(crop_config_global%rotation_loaded) .and. &
          crop_config_global%rotation_loaded(icrop)) then
         call cropgrass_init_from_config(crop_config_global%rotation_grass(icrop), icrop, &
                                         swharvest, dmharvest, daylastharvest, dmlastharvest, &
                                         swdmmow, maxdaymow, swlossmow, swlossgrz, swdmgrz, &
                                         maxdaygrz, dmgrazing, LSDb, tagprest, swhydrlift)
      else
         call readgrass (icrop, cropfil(icrop), swharvest, dmharvest, daylastharvest, &
                         dmlastharvest, swdmmow, maxdaymow, swlossmow, swlossgrz, swdmgrz, &
                         maxdaygrz, dmgrazing, LSDb, tagprest, swhydrlift)
      end if
```

**After (Phase 4):**
```fortran
      call cropgrass_init_from_config(crop_config_global%rotation_grass(icrop), icrop, &
                                      swharvest, dmharvest, daylastharvest, dmlastharvest, &
                                      swdmmow, maxdaymow, swlossmow, swlossgrz, swdmgrz, &
                                      maxdaygrz, dmgrazing, LSDb, tagprest, swhydrlift)
```

**Note on grass signature:** Same caveat as wofost — the implementing agent must read
Phase 3's `src/crop/cropgrass_init.f90` to confirm the exact argument list from the
`readgrass` legacy call at `cropgrowth.f90:2068`.

### What stays after Phase 4

| Element | Status after Phase 4 | Teardown trigger |
|---|---|---|
| `readcropfixed` / `readwofost` / `readgrass` subs in `readswap.f90` | Remain alive as parity-test fixtures | Future cleanup task (rewrite parity tests against fixture values) |
| `crop_config_global` module-level pointer (`src/crop/crop_config_global.f90`) | Remains, still used by the three unconditional init calls | ADR 0016 Part C: config-passing follow-on spec |
| `rotation_loaded(:)` sentinel array on `crop_config_t` | Remains on the type | Same as `crop_config_global` |
| `rotation_fixed(:)` / `rotation_wofost(:)` / `rotation_grass(:)` parallel arrays | Remain — these are the long-term design per ADR 0016 | Not transitional |
| Legacy `.crp` files in `tests/swap-cases/1.hupselbrook/` (the legacy case dir) | Remain for the legacy executable | Legacy executable retirement (out of scope) |

### Tombstone: the three `.crp` files deleted from the TOML dir

Phase 4 deletes:
- `tests/swap-cases/toml/1.hupselbrook/maizes.crp`
- `tests/swap-cases/toml/1.hupselbrook/potatod.crp`
- `tests/swap-cases/toml/1.hupselbrook/grassd.crp`

(These are the TOML-dir copies that the legacy reader was finding; the originals in
`tests/swap-cases/1.hupselbrook/` are untouched.)

---

## Test plan

### Pre-flight (blocking on Phases 1–3)

```bash
pixi run test-pfunit 2>&1 | tail -5
pixi run regression 2>&1 | tail -10
```

Both must be green before any Phase 4 work begins.

In addition, the following pFUnit assertions (written in Phases 1–3) must pass:

- `test_load_case6_grass_crp_toml_populated` (Phase 1): all three rotation_loaded slots
  green for case 6.
- `test_load_case5_potatod_crp_toml_populated` (Phase 2): rotation_loaded for case 5.
- `test_load_case4_grassd_crp_toml_populated` and `test_load_case2_grassd_crp_toml_populated`
  (Phase 3): rotation_loaded for cases 4 and 2.

### New pFUnit assertion: case 1 rotation-loaded

Extend `tests/unit/io/toml/test_load_swap_config.pf` with:

```fortran
@test
subroutine test_load_case1_hupselbrook_all_rotations_loaded()
   use funit
   use error_mod, only: error_collection_t
   use swap_config_mod, only: swap_config_t
   use load_swap_config_mod, only: load_swap_config
   type(swap_config_t)      :: c
   type(error_collection_t) :: errors
   call load_swap_config('tests/swap-cases/toml/1.hupselbrook/swap.toml', c, errors)
   @assertFalse(errors%has_errors())
   @assertEqual(3, size(c%crop%rotation_loaded))
   ! Slot 1: type-1 (cropfixed) maize
   @assertTrue(c%crop%rotation_loaded(1))
   @assertEqual(1, c%crop%rotation_fixed(1)%idev)
   @assertEqual(168, c%crop%rotation_fixed(1)%lcc)
   @assertEqual(0.6d0, c%crop%rotation_fixed(1)%kdif, 1.0d-12)
   @assertEqual(2, c%crop%rotation_fixed(1)%swcf)
   @assertEqual(0.23d0, c%crop%rotation_fixed(1)%albedo, 1.0d-12)
   ! Slot 2: type-2 (cropwofost) potato
   @assertTrue(c%crop%rotation_loaded(2))
   @assertEqual(0, c%crop%rotation_wofost(2)%idsl)
   @assertEqual(2, c%crop%rotation_wofost(2)%swgerm)
   @assertEqual(150.0d0, c%crop%rotation_wofost(2)%tsumea, 1.0d-12)
   ! Slot 3: type-3 (cropgrass) grass
   @assertTrue(c%crop%rotation_loaded(3))
   @assertEqual(3, c%crop%rotation_grass(3)%swrd)
   @assertEqual(2, c%crop%rotation_grass(3)%swcf)
end subroutine
```

### Regression tests

- After each `.crp.toml` authoring task: run `pixi run regression` and verify 5/5 green.
- After dispatch-removal tasks: run `pixi run regression`, verify 5/5 green.
- Smoke test (before deletion): rename all three `.crp` files in the TOML dir, run
  regression, verify 5/5 green, restore.
- After actual deletion: run regression, verify 5/5 green as the acceptance gate.

### Acceptance gate

- `pixi run test-pfunit` → all green.
- `pixi run regression` → 5/5 green.
- `tests/swap-cases/toml/1.hupselbrook/maizes.crp` → does not exist.
- `tests/swap-cases/toml/1.hupselbrook/potatod.crp` → does not exist.
- `tests/swap-cases/toml/1.hupselbrook/grassd.crp` → does not exist.
- `git grep "call readcropfixed" src/` → zero results (all reachable calls removed).
- `git grep "call readwofost" src/` → zero results.
- `git grep "call readgrass" src/` → zero results.
- `git grep "readcropfixed\|readwofost\|readgrass" tests/` → parity-test references only.
- `crop_config_global` pointer still present; `rotation_loaded(:)` still present on
  `crop_config_t`.

---

## Implementation tasks (summary for the plan)

1. Pre-flight: confirm Phases 1–3 landed; all pFUnit assertions for cases 6/5/4/2 pass.
2. Author `maizes.crp.toml` 1:1 + regression smoke + submodule pair commit.
3. Author `potatod.crp.toml` verification + any gap-fill + submodule pair commit.
4. Author `grassd.crp.toml` 1:1 + regression smoke + submodule pair commit.
5. Regression test case 1 with all three `.crp.toml` files complete; fix any mismatches.
6. Write the pFUnit case-1 rotation-loaded assertion; verify it passes.
7. Remove the type-1 (cropfixed) `if/else` dispatch fallback from `cropgrowth.f90`.
8. Remove the type-2 (cropwofost) `if/else` dispatch fallback from `cropgrowth.f90`.
9. Remove the type-3 (cropgrass) `if/else` dispatch fallback from `cropgrowth.f90`.
10. Smoke test: rename all three `.crp` files in the TOML dir; run regression; confirm
    5/5 green; restore.
11. Delete all three `.crp` files from the TOML dir + submodule pair commit.
12. Docs update: add Phase 4 to `docs/csv-companion-files.md` migration history;
    confirm ADR 0016's teardown table is accurate.

---

## Risk register

| Risk | Severity | Likelihood | Mitigation |
|---|---|---|---|
| **Multi-type rotation integration** (all three init modules run side-by-side in case 1). If any two modules share module-global state that they each initialize independently, they may clobber each other's values mid-simulation (e.g. if `variables::swrd` is set by both `cropfixed_init` and `cropgrass_init`). | High | Low | Each rotation runs its own init at task=1 in sequence (one per rotation entry, not simultaneously). The rotation index `icrop` scopes the init. However, if any init module reads module globals left over from a prior rotation (e.g. `swrd` from the previous `cropfixed_init` being visible when `cropgrass_init` runs), it can cause subtle bugs. Mitigation: verify in the Phase 1/2/3 parity tests that each `*_init_from_config` writes **all** module globals its runtime path reads, so no stale values propagate. |
| **SWJARVIS=4 not covered by Phase 3** (grassd.crp uses `SWJARVIS=4`; Phase 3 may only support `SWJARVIS=0`). | Medium | Medium | If Phase 3 stub-errors `SWJARVIS /= 0`, case 1 rotation 3 will fail validation when loading grassd.crp.toml. Fix: retroactively extend Phase 3 to support `SWJARVIS ∈ {0,1,2,3,4}` (or the subset that case 1 uses). This is a Phase 3 retroactive concern, not Phase 4 work. Document and escalate to the Phase 3 implementer. |
| **potatod.crp.toml appears near-complete from a prior task** but may have minor gaps or stale field names. | Low | Medium | Phase 4 Task 3 includes a field-by-field diff of the TOML file against the 704-line legacy `.crp`. Any gap found there is minor authoring work. |
| **Dispatch argument lists for wofost/grass `_init_from_config`** may differ from the predicted signatures above. | Low | Low | The implementing agent reads the actual Phase 2/3 source before writing the post-removal code. The predicted signatures in this spec are derived from the legacy reader call signatures; adjust if Phase 2/3 changed them. |
| **Regression of cases 6/5/4/2** may regress after the dispatch removal if the unconditional calls expose a latent bug in the init modules that was masked by the fallback. | Low | Low | The smoke test (Task 10) catches this before deletion. The pFUnit assertions for all four predecessor cases also run in the pre-flight check. |
