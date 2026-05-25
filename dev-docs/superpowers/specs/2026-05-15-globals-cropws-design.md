---
title: "GR-CROPWS — Crop Write-Site Refactor (cropgrowth split + intent flip + write migration)"
date: 2026-05-15
status: approved
context: globals-retirement continuation; structural split + write-site cutover
---

# GR-CROPWS Design

## Goal

The cropgrowth retirement unlocker. Three interlocking efforts:

1. **Phase 0 (structural)**: split `src/crop/cropgrowth.f90` (4937L) into 5 focused files by subroutine boundary. Each file becomes 270-1400L. Public module preserved; byte-for-byte trivial.

2. **Phase A (architectural)**: flip `intent(in) :: state` → `intent(inout) :: state` in cropgrowth + irrigation + oxygenstress + tillage; migrate ~200-300 write sites to write directly into `state%crop%X` (dropping legacy global writes); retire Phase A.5 mirror lines (now redundant — state IS the write target).

3. **Phase B + C (retirement)**: drop adapter writes, dual-write seeds, zero-fills; tombstone 60-80 crop globals in `variables.f90`; iterative compile patches.

After this arc: variables.f90 ~1418L → <1000L (30% shrinkage); cropgrowth split into manageable files; future arcs operate per-file rather than per-mega-file.

## Architecture (4 phases)

- **Phase 0** — File split (5 commits, pure refactor)
- **Phase A** — Intent flip + write migration per file (7 commits)
- **Phase B** — Adapter cleanup (3-5 commits)
- **Phase C** — variables.f90 retirement (iterative, 3-5 commits)
- **Phase D** — Close + roadmap update

## Phase 0 — cropgrowth.f90 Split

5 new files (or 4 new + 1 rename), 5 commits. Public module name kept as `cropgrowth_mod` (umbrella) OR each file has its own module — implementer decides per file's natural interface.

| New file | Subroutines extracted | ~Lines |
|---|---|---|
| `src/crop/cropgrowth_helpers.f90` | `nocrop`, `ArableLandGerm`, `FacCO2`, `cropoutput`, `update_rootdistribution`, `sumttd` (plus any private helpers used only by these) | ~1300 |
| `src/crop/cropwofost_runtime.f90` | `wofost` | ~1300 |
| `src/crop/cropgrass_runtime.f90` | `grass` | ~1400 |
| `src/crop/cropfixed_runtime.f90` | `cropfixed` | ~270 |
| `src/crop/cropgrowth.f90` (renamed/kept) | `CropGrowth` (main dispatcher only) | ~600 |

**Approach for each split commit:**
1. Create new file with extracted subroutine
2. Adjust `meson.build` (add new sources before `cropgrowth.f90` in legacy list — preserving ordering)
3. Add `use ... only:` imports for the extracted subroutine's dependencies (state types, helper modules, etc.)
4. Update the caller (`CropGrowth`) to `use <new_module>, only: <subname>` if needed
5. Verify byte-for-byte regression (VG)
6. Commit

**Order of extraction (least-dependent first):**
- 0.1: `cropgrowth_helpers.f90` (helpers are leaf-level)
- 0.2: `cropwofost_runtime.f90` (calls helpers but not other runtimes)
- 0.3: `cropgrass_runtime.f90` (similar)
- 0.4: `cropfixed_runtime.f90` (similar)
- 0.5: rename `cropgrowth.f90` (keeps only CropGrowth dispatcher) + meson cleanup

VG byte-for-byte trivial across all 5 commits (no behavior change).

## Phase A — Intent Flip + Write-Site Migration

7 commits. Per-file: flip `intent(in)` → `intent(inout)` in subroutine signatures; migrate write sites; drop Phase A.5 mirror lines.

### Symbol → state%X table (target writes)

| Was (legacy global write) | Now (state-only write) |
|---|---|
| `daycrop = ...` | `state%crop%common%daycrop = ...` |
| `dvs = ...` | `state%crop%common%dvs = ...` |
| `tsum = ...` | `state%crop%common%tsum = ...` |
| `rd = ...`, `rdpot = ...`, etc. | `state%crop%common%X = ...` |
| `ch = ...`, `cf = ...`, `laipot = ...` | `state%crop%common%X = ...` |
| `flCrop*`/`swend` | `state%crop%common%X = ...` |
| `flCropEmergence` | `state%crop%flCropEmergence = ...` (top-level) |
| `lai/kdif/kdir/cofab/cfbs/swcf/swcfbs/gird/et0/ew0/es0` | `state%crop%X = ...` (top-level) |
| `wlv/wlvpot/wst/wstpot/wrt/wrtpot/wso/wsopot/...` | `state%crop%wofost%X = ...` |
| `dwlv/dwlvpot/dwst/dwstpot/dwrt/dwrtpot/dwso/dwlvCrop/dwlvSoil/plossdm/lossdm` | `state%crop%wofost%X = ...` |
| `tagp/tagppot/tagpt/tagptpot/cwdm/cwdmpot/pgass/pgasspot` | `state%crop%wofost%X = ...` |
| `swbulb/wbl/wblpot/dwbl/dwblpot/plwt/plwti` | `state%crop%wofost%X = ...` |
| `fco2amax/fco2eff/fco2tra` | `state%crop%wofost%X = ...` |
| `cropstart/cropend/PrepDelay/SowDelay/noddrz/cumdens/albedo/rsc` | `state%crop%common%X = ...` |
| `seqgrazmow/seqgrazmowpot/dateharvest/mowrest/cropstartpot/cropstartact/cropendpot/cropendact/swpotrelmf/relmf` | `state%crop%grass%X = ...` |
| `cftb/chtb/cfeic/cfeictb` | `state%crop%fixed%X = ...` |

### Commits

- **A.1**: `cropgrowth.f90` (CropGrowth dispatcher only after split) — flip intent + migrate writes
- **A.2**: `cropfixed_runtime.f90`
- **A.3**: `cropwofost_runtime.f90`
- **A.4**: `cropgrass_runtime.f90`
- **A.5**: `cropgrowth_helpers.f90`
- **A.6**: `irrigation.f90` + `oxygenstress.f90` + `tillage.f90` (any crop write sites in these)
- **A.7**: Phase A close marker + audit (all 4 original consumer files have inout state OR no crop writes)

### Per-file approach

1. **Identify intent change scope:** does this file's main subroutine WRITE crop fields? If yes, flip `intent(in)` → `intent(inout)`. If no, leave intent alone.
2. **Migrate write sites:** for each legacy `<sym> = <expr>` line, rewrite as `state%crop%X = <expr>`.
3. **Drop now-redundant mirror lines:** the Phase A.5 mirror lines (`state%crop%X = <sym>` immediately after a legacy write) become redundant once the write itself is state-only. Drop them.
4. **Callers don't break:** verify each caller of the affected subroutines passes a mutable state (most do; verify per `grep -rn "call <subname>"`).
5. VG byte-for-byte → commit.

**Critical zero-reset coverage:** the `nocrop()` subroutine resets many crop fields to zero. Each `<sym> = 0` must become `state%crop%X = 0`. Don't miss any (GR-ATM `8fb76cf` lesson).

## Phase B — Adapter Cleanup

After Phase A, the migrated crop globals are write-only-by-adapter. Drop the adapter machinery.

### Commits

- **B.1**: Drop swap_mod A14-A17 dual-write SEED lines for retire-able crop globals (the lines that read `legacy_sym` and write to `state%crop%X` at init — now redundant because the runtime writes state directly).

Note: the SEED lines remain ONLY where the legacy global carries non-trivial initial state from the adapter (e.g., crop config loaded once via config_to_variables that writes the legacy global). For most crop runtime fields, the seed lines were always redundant; this drops them cleanly.

- **B.2**: Drop adapter writes in `src/io/toml/config_to_variables.f90` for retire-able crop globals. Where state seeding from config is needed, source from `config%crop%X` directly inline in swap_mod (per Phase C of GR-FINAL precedent).

- **B.3**: Drop zero-fills in `src/core/initialize.f90` for retire-able crop globals.

- **B.4**: Phase B close marker.

## Phase C — variables.f90 Retirement

Iterative tombstoning of crop globals + compile-error patch loop.

### Approach

1. **Identify retire-able crop globals:** run audit grep:
```bash
for sym in $RETIRE_CANDIDATES; do
  hits=$(grep -rnE "\b${sym}\b" src/ --include="*.f90" \
    | grep -v "state%\|config%\|! \|variables.f90\|use variables, only:\|initialize.f90\|config_to_variables\|swap_mod.f90\|integer ::\|real(8) ::\|real(real64) ::\|character\|logical ::" | wc -l)
  if [ "$hits" -eq 0 ]; then echo "RETIRE: $sym"; fi
done
```

2. **Tombstone in batches:** 10-15 symbols per commit. Each commit:
   - Replace declaration with `! [SS-GR-CROPWS C] <sym> retired — state%crop%X is sole home`
   - Iterative `rm -rf builddir && pixi run build-linux 2>&1 | grep Error | head` — patch each compile error
   - VG → commit

### Expected retirement set (~60-80 symbols)

`daycrop, dvs, tsum, rd, rdpot, rdm, rri, rdi, rdc, ch, cf, laipot, cuptgraz, cuptgrazpot, HarLosOrm_tot, flCropCalendar, flCropOutput, flCropNut, flHarvestDay, flCropEmergence, flCropReadFile, flCropPrep, flCropSow, flCropGerm, flCropHarvest, cropstart, cropend, PrepDelay, SowDelay, noddrz, cumdens, albedo, rsc, wlv, wlvpot, wst, wstpot, wrt, wrtpot, wso, wsopot, tagp, tagppot, tagpt, tagptpot, cwdm, cwdmpot, pgass, pgasspot, dwlv, dwlvpot, dwst, dwstpot, dwrt, dwrtpot, dwso, dwlvCrop, dwlvSoil, plossdm, lossdm, swbulb, wbl, wblpot, dwbl, dwblpot, plwt, plwti, fco2amax, fco2eff, fco2tra, seqgrazmow, seqgrazmowpot, dateharvest, mowrest, cropstartpot, cropstartact, cropendpot, cropendact, swpotrelmf, relmf, cftb, chtb, cfeic, cfeictb, swcrp, icrop`

That's ~85 symbols. Final batch may be smaller depending on which still have non-write readers (mostly atmosphere/swap_mod orchestration).

### Final close marker

variables.f90 line count audit; check-full 5/5; commit close marker.

## Phase D — Close + Roadmap

- **D.1**: Update `docs/superpowers/plans/2026-05-13-globals-retirement-roadmap.md` with residuals
- **D.2**: ADR 0045 candidate (NEW) — "cropgrowth structural split + write-site refactor"
- **D.3**: Arc-complete marker

## Verification Discipline

Per `feedback_per_task_regression_gate.md` + `feedback_state_schema_clean_rebuild.md`:

```bash
rm -rf builddir && pixi run build-linux \
  && pixi run test-pfunit \
  && pixi run -e test python tests/regression/test_output_regression.py \
       hupselbrook surfacewater salinitystress grassgrowth
```

Expected: clean build, pFUnit 741, regression **4/4 byte-for-byte**. Per commit. Phase close: check-full 5/5.

## Controller Intervention Policy

Per user directive: subagents do NOT defer. If something blocks (caller chain extends out of crop, missing state field, etc.), controller intervenes inline.

Specific likely interventions:
- A caller in non-crop file passes `state` with `intent(in)`: thread inout
- Missing state field surfaces: extend schema inline
- Write site has subtle semantics (e.g., write inside macropore-retired branch): inspect carefully

## Risk & Mitigation

| Risk | Mitigation |
|---|---|
| Split breaks visible interface | All public symbols accessible via the umbrella module `cropgrowth_mod` OR explicit `use <new_module>, only: <subname>`. Verify after Phase 0 via grep on caller signatures. |
| Intent flip cascades to non-crop callers | Most are crop-internal; controller threads inout outward if needed |
| Write-site density ~300 sites → easy to miss | Per-file commits + VG per commit + audit grep at end of Phase A |
| Zero-resets in nocrop() critical (GR-ATM lesson) | Explicit nocrop audit step in A.5 |
| Phase C compile errors cascade | Batched retirement (10-15 symbols per commit) + iterative patch |

## What This Arc DOES NOT Change

- BMI / cffi surfaces unchanged
- Physics — byte-for-byte regression preserved
- Crop config schemas — already populated from prior arcs

## Effort

- **Phase 0**: 5 commits (~1-2 days)
- **Phase A**: 7 commits (~3-4 days — biggest work)
- **Phase B**: 4-5 commits (~1 day)
- **Phase C**: 3-5 commits (~1-2 days)
- **Phase D**: 2-3 commits (~half day)
- **Total**: ~25-35 commits, **6-9 days subagent-driven**

## Files & Artifacts

- **Spec**: `docs/superpowers/specs/2026-05-15-globals-cropws-design.md`
- **Plan**: `docs/superpowers/plans/2026-05-15-globals-cropws.md` (this arc's spec is detailed enough that the plan is mainly a task checklist mirroring the phases)
- **Roadmap context**: `docs/superpowers/plans/2026-05-13-globals-retirement-roadmap.md`
- **Prior arcs**: GR-CROPRT (`d2ed001`), GR-FINAL (`21e15dc`), GR-CROP (`476efba`), GR-ATM (`57b00b3`), GR-BH (`a65bdf3`), GR-UTILS (`26123ad`)
