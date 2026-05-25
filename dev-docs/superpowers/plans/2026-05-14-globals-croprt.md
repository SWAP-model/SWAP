# GR-CROPRT Implementation Plan — Crop Runtime Finalization + Cheap Wins

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Retire AgeTracer + macropore dead-code globals; migrate output control switches; extend state%crop with ~15-20 fields; migrate cropgrowth.f90's 10 use-variables sites (~738+ symbol substitutions); retire ~50-80 crop runtime globals from variables.f90.

**Architecture:** Five phases — Phase A (cheap wins + schema additions), Phase A.5 (runtime dual-writes), Phase B (cropgrowth migration site-by-site), Phase C (adapter shrinkage), Phase D (variables.f90 deletion), Phase E (close).

**Tech Stack:** Fortran 2018, gfortran, Meson, pFUnit, pixi, regression suite.

**Spec:** `docs/superpowers/specs/2026-05-14-globals-croprt-design.md`

**Builds on:** GR-FINAL partial close (`21e15dc`).

**Verification gate per task (VG):**

```bash
rm -rf builddir && pixi run build-linux \
  && pixi run test-pfunit \
  && pixi run -e test python tests/regression/test_output_regression.py \
       hupselbrook surfacewater salinitystress grassgrowth
```

Expected: clean build, pFUnit 741, regression **4/4 byte-for-byte**.

**Phase close gate:** `pixi run check-full` 5/5 byte-for-byte.

**Controller quality bar:** subagents do NOT defer due to missing infrastructure. Inline-fix per user directive.

---

## Phase A — Cheap Wins + Crop Schema Additions

### Task 1: Pre-flight baseline

- [ ] `rm -rf builddir && pixi run build-linux` — succeed
- [ ] `pixi run check-full` — 5/5
- [ ] Empty marker commit:
```bash
git commit --allow-empty -m "chore(gr-croprt): pre-flight baseline

check-full 5/5. Baseline locked before AgeTracer + macropore retirement,
output switch migration, crop schema expansion, and cropgrowth.f90
finalization (10 sites, ~738+ substitutions).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

### Task 2: A1 — AgeTracer dead-code retirement

**Files:** `src/core/variables.f90`, `src/io/toml/config_to_variables.f90`, `src/core/initialize.f90`, `src/solute/agetracer.f90`

ADR 0032: `flAgeTracer` is always `.false.`; body is unreachable.

- [ ] **Step 1: Verify guard.**
```bash
grep -rn "flAgeTracer" src/ --include="*.f90" | head -20
```
Confirm `flAgeTracer` is never set to `.true.` anywhere.

- [ ] **Step 2: Tombstone declarations** in `src/core/variables.f90` for: `flAgeTracer`, `AgeGwl1m`, `icAgeBot`, `icAgeDra`, `icAgeRot`, `icAgeSur`. Replace each declaration with:
```fortran
   ! [SS-GR-CROPRT A1] <sym> retired — ADR 0032 (AgeTracer dead-code)
```

- [ ] **Step 3: Drop adapter writes** in `src/io/toml/config_to_variables.f90`.

- [ ] **Step 4: Drop zero-fills** in `src/core/initialize.f90`.

- [ ] **Step 5: Audit + patch readers** in `src/solute/agetracer.f90` and any callers. Remove `if (flAgeTracer)` guards if the entire branch is dead, OR leave as `if (.false.)` (effectively dead-code but explicit).

- [ ] **Step 6: Iterative build.** `rm -rf builddir && pixi run build-linux 2>&1 | grep Error | head`. Patch each surfaced error.

- [ ] **Step 7: VG.**

- [ ] **Step 8: Commit:**
```bash
git add -A
git commit -m "retire(gr-croprt A1): AgeTracer dead-code globals (ADR 0032)

flAgeTracer/AgeGwl1m/icAge* retired. Body in agetracer.f90 is unreachable
per the always-false guard. <document any reader patches>.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

### Task 3: A2 — Macropore retired sentinel cleanup

**Files:** `src/core/variables.f90`, `src/io/toml/config_to_variables.f90`, `src/core/initialize.f90`, possibly `src/soil/*.f90`

ADR 0040: macropore retired.

- [ ] **Step 1: Inventory macropore sentinels.**
```bash
grep -rn "FlMacropore\|ArMpSs\|iQInTopLatDm\|IWaSrDm\|WaSrDm" src/ --include="*.f90" | head -30
```

- [ ] **Step 2: Verify no live readers.** Many `.BMA` writer globals are write-only-by-adapter — verify zero non-write consumers.

- [ ] **Step 3: Tombstone declarations + drop adapter writes + zero-fills** for each macropore sentinel.

- [ ] **Step 4: Remove `if (flmacropore)` guarded branches** if entire branches dead.

- [ ] **Step 5: Iterative build + patch.**

- [ ] **Step 6: VG.**

- [ ] **Step 7: Commit:**
```bash
git commit -m "retire(gr-croprt A2): macropore retired sentinels (ADR 0040)

FlMacropore + .BMA writer globals + ArMpSs residuals retired. Dead
branches dropped per ADR 0040.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

### Task 4: A3 — Output control switches migration

**Files:** `src/config/output_csv_config.f90` (may extend), `src/io/toml/config_to_variables.f90`, consumer files

- [ ] **Step 1: Inventory switches.**
```bash
grep -rnE "\b(swcsv|swcsv_tz|swinc|swdrought|swrum|swtem|swvap|swsnow_output)\b" src/ --include="*.f90" | head -30
```

- [ ] **Step 2: For each switch, determine homing.**
- Already in `config%output_csv` → migrate readers to use config
- Not yet in config → add field to `output_csv_config_t`, populator in adapter

- [ ] **Step 3: Migrate readers** in swapoutput.f90 + any other consumers. Thread config arg if needed.

- [ ] **Step 4: VG.**

- [ ] **Step 5: Commit:**
```bash
git commit -m "refactor(gr-croprt A3): output control switches → config%output_csv

<list switches migrated>. File unit globals (logf/inc/rot/crp/tem/snw)
retained as Arc 9.5 pragmatic edge.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

### Task 5: A4 — crop schema additions

**Files:** `src/state/crop_common_state.f90`, `src/state/crop_wofost_state.f90`, possibly `src/state/crop_grass_state.f90`, `src/state/nutrients_state.f90`

- [ ] **Step 1: Extend `crop_common_state_t`.**

Edit `src/state/crop_common_state.f90`. Add to the type body:

```fortran
      ! [SS-GR-CROPRT A4] lifecycle flags
      logical :: flCropReadFile = .false.
      logical :: flCropPrep     = .false.
      logical :: flCropSow      = .false.
      logical :: flCropGerm     = .false.
      logical :: flCropHarvest  = .false.

      ! Crop window dates
      real(real64) :: cropstart = 0.0_real64
      real(real64) :: cropend   = 0.0_real64

      ! Germination delay scratch
      integer :: PrepDelay = 0
      integer :: SowDelay  = 0

      ! Runtime root zone
      integer :: noddrz = 0
      real(real64), allocatable :: cumdens(:)   !! sized numnod in init

      ! Surface params (verify if config or state per audit)
      real(real64) :: albedo = 0.0_real64
      real(real64) :: rsc    = 0.0_real64
```

Update `crop_common_state_init` to allocate `cumdens(numnod)` (need to take numnod arg).

- [ ] **Step 2: Extend `crop_wofost_state_t`.**

Edit `src/state/crop_wofost_state.f90`:

```fortran
      ! [SS-GR-CROPRT A4] FCO2 derived values
      real(real64) :: fco2amax = 1.0_real64   !! verify legacy default
      real(real64) :: fco2eff  = 1.0_real64
      real(real64) :: fco2tra  = 1.0_real64
```

Verify legacy defaults via `grep -nE "fco2amax|fco2eff|fco2tra" src/core/initialize.f90`.

- [ ] **Step 3: Audit grass-specific needs.**
```bash
grep -n "use variables" src/crop/cropgrowth.f90 | head -3
sed -n '2598,2620p' src/crop/cropgrowth.f90
```
If grass-specific runtime state surfaces beyond what's already in `crop_grass_state_t`, add.

- [ ] **Step 4: Audit nutrient overflow.** If `anlv`, `anst`, `nni`, `nmaxlv`, etc. aren't in state%nutrients, decide whether to add (nutrients arc territory) OR retain narrow defer.

- [ ] **Step 5: Update swap_mod to pass numnod to crop_common_state%init.**

Find the call in `src/core/swap_mod.f90`:
```bash
grep -n "state%crop%common%init\|crop_common_state_init" src/core/swap_mod.f90
```

If init takes no args, change signature to accept numnod and update call site. Allocate cumdens inside.

- [ ] **Step 6: VG.**

- [ ] **Step 7: Commit:**
```bash
git add src/state/crop_common_state.f90 src/state/crop_wofost_state.f90 src/state/crop_grass_state.f90 src/state/nutrients_state.f90 src/core/swap_mod.f90
git commit -m "schema(gr-croprt A4): crop schema additions

crop_common: lifecycle flags, cropstart/cropend, delay scratch, noddrz,
cumdens, albedo, rsc. crop_wofost: fco2amax/eff/tra. Unpopulated until
A.5 + B migration.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

### Task 6: Phase A close marker

- [ ] Run VG (and check-full 5/5).
- [ ] Empty marker:
```bash
git commit --allow-empty -m "chore(gr-croprt): Phase A complete — cheap wins + schema landed

AgeTracer + macropore dead-code retired (~30-50 globals). Output switches
migrated to config. crop schema extended with ~15-20 fields.
check-full 5/5. Phase A.5 (runtime dual-writes) unblocked.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Phase A.5 — Runtime Dual-Write Coverage

### Task 7: A.5 — Dual-write mirrors for new crop fields

**Files:** `src/crop/cropgrowth.f90`, `src/crop/cropfixed_init.f90`, `src/crop/cropgrass_init.f90`, `src/crop/cropwofost_init.f90`, `src/crop/irrigation.f90`, `src/atmosphere/meteoday.f90`

- [ ] **Step 1: Per-symbol audit.** For each new field added in A4, grep all write sites:
```bash
for sym in flCropReadFile flCropPrep flCropSow flCropGerm flCropHarvest cropstart cropend PrepDelay SowDelay noddrz cumdens albedo rsc fco2amax fco2eff fco2tra; do
  echo "=== $sym ==="
  grep -rnE "^[[:space:]]*\b$sym\b\s*=" src/ --include="*.f90"
done
```

- [ ] **Step 2: For each write site, add mirror.**

Pattern (per Phase A.5 doctrine from prior arcs):
```fortran
<sym> = <expr>
state%crop%common%<sym> = <sym>   ! [SS-GR-CROPRT A5]
```

For zero-resets (e.g., `flCropPrep = .false.`), include mirrors in same conditional block.

**Allocate cumdens:** the `cumdens` allocatable should be sized to numnod. Allocate during crop_common%init OR at first write site if dynamic. Verify the legacy `cumdens` is declared as `cumdens(macp)` — match the size.

- [ ] **Step 3: VG.**

- [ ] **Step 4: Commit per file batch (3-5 commits total):**

Commit 1: cropgrowth.f90 dual-writes
Commit 2: crop init files (cropfixed/cropgrass/cropwofost) dual-writes
Commit 3: irrigation.f90 / meteoday.f90 if any

Example commit:
```bash
git add src/crop/cropgrowth.f90
git commit -m "schema(gr-croprt A5.1): cropgrowth.f90 dual-writes for new state fields

Mirror lifecycle flags + crop dates + noddrz/cumdens + albedo/rsc +
fco2amax/eff/tra writes to state%crop sub-records. Zero-resets covered.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

### Task 8: A.5 close marker

```bash
git commit --allow-empty -m "chore(gr-croprt): Phase A.5 complete — runtime dual-writes landed

State mirrors track legacy globals through simulation for all new crop
fields. Phase B reader migration unblocked.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Phase B — cropgrowth.f90 Reader Migration

### Tasks 9-18: Sites 1-10

One commit per site. Each task follows the same pattern.

**Pattern (apply to each site):**

1. **Inspect** the site:
```bash
sed -n '<line-5>,<line+50>p' src/crop/cropgrowth.f90
```

2. **Replace import** — drop `use variables, only: ...` (or narrow further to only ETSine/file-path deferrals):
```fortran
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   ! Narrow deferrals only: rad, daylp, difpp, dsinbe, atmtr (ETSine — deferred)
   ! And: outfil, pathwork, project, cropfil (file paths — deferred)
```

3. **Thread config arg** if not already present.

4. **Substitute body refs** per the symbol replacement table:

| Was | Now |
|---|---|
| daycrop/dvs/tsum/swcrp/icrop/swend/rd/rdpot/rdm/rri/rdi/rdc/ch/cf/laipot/cuptgraz/cuptgrazpot/HarLosOrm_tot | `state%crop%common%X` |
| flCropCalendar/flCropOutput/flCropNut/flHarvestDay/flCropEmergence/flCropReadFile/flCropPrep/flCropSow/flCropGerm/flCropHarvest | `state%crop%common%X` (per A4) or top-level `state%crop%flCropEmergence` |
| cropstart/cropend/PrepDelay/SowDelay/noddrz/cumdens/albedo/rsc | `state%crop%common%X` (per A4) |
| wlv/wlvpot/wst/wstpot/wrt/wrtpot/wso/wsopot/tagp/tagppot/tagpt/tagptpot/cwdm/cwdmpot/pgass/pgasspot/dwlv/dwlvpot/dwst/dwstpot/dwrt/dwrtpot/dwso/dwlvCrop/dwlvSoil/plossdm/lossdm/swbulb/wbl/wblpot/dwbl/dwblpot/plwt/plwti | `state%crop%wofost%X` |
| fco2amax/fco2eff/fco2tra | `state%crop%wofost%X` (per A4) |
| seqgrazmow/seqgrazmowpot/dateharvest/mowrest/cropstartpot/cropstartact/cropendpot/cropendact/swpotrelmf/relmf | `state%crop%grass%X` |
| cftb/chtb/cfeic/cfeictb | `state%crop%fixed%X` |
| lai/kdif/kdir/cofab/cfbs/swcf/swcfbs/gird/et0/ew0/es0 | `state%crop%X` (top-level) |
| Atmosphere refs (Tav/tavd/rh/arad/etc.) | `state%atmosphere%X` |
| Soil refs (orgmat/psand/etc.) | `state%soilwater%X` |
| Mesh refs (numnod/dz/z/etc.) | `state%mesh%X` |
| Config-loaded crop params (idev/tsumea/tbase/etc.) | `config%crop%X` (via threaded arg) |
| WOFOST nutrient params (lrnr/lsnr/nlue/etc.) | `config%crop%wofost%nutrient%X` |
| ETSine scratchpad (rad/daylp/difpp/dsinbe/atmtr/lat) | Retain narrow defer — `[SS-GR-CROPRT B] DEFERRED ETSine` |
| File paths (outfil/pathwork/project/cropfil) | Retain narrow defer |

5. **Update callers.**
```bash
grep -rn "call <subname>\b" src/
```

6. **VG → commit:**
```bash
git add src/crop/cropgrowth.f90 # + callers
git commit -m "refactor(gr-croprt B<N>): cropgrowth.f90 site <line> — <subroutine>

<symbol count migrated>. <narrow deferrals retained, if any>.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

**Sites in order (Task 9 = site 1, Task 10 = site 2, etc.):**

| Task | Line | Subroutine |
|---|---|---|
| 9 | 44 | CropGrowth (main dispatcher, ~80 symbols, biggest of small) |
| 10 | 621 | cropfixed (~18 symbols) |
| 11 | 894 | cropoutput (file paths — small) |
| 12 | 1005 | nocrop (**zero-reset CRITICAL**) |
| 13 | 1053 | ArableLandGerm (germination) |
| 14 | 1243 | FacCO2 (fco2amax/eff/tra) |
| 15 | 1313 | wofost (~60 symbols, **biggest site**) |
| 16 | 2598 | grass (~50 symbols) |
| 17 | 4024 | update_rootdistribution (noddrz/cumdens) |
| 18 | 5198 | sumttd (small) |

### Task 19: Phase B close marker

- [ ] Verify cropgrowth.f90 is clean:
```bash
grep -n "^[[:space:]]*use [Vv]ariables" src/crop/cropgrowth.f90
```
Expected: empty or only documented narrow deferrals (ETSine + file paths).

- [ ] check-full 5/5.

```bash
git commit --allow-empty -m "chore(gr-croprt): Phase B complete — cropgrowth.f90 migrated

10 sites migrated; ~738+ symbol substitutions. Remaining narrow imports:
ETSine astronomical scratchpad + file paths (documented deferrals).
check-full 5/5 byte-for-byte.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Phase C — Adapter Shrinkage

### Task 20: C1 — Drop adapter writes for retired crop globals

**Files:** `src/io/toml/config_to_variables.f90`

After Phase B migration, ~50-80 crop globals are write-only-by-adapter.

- [ ] **Step 1: Identify candidates.** For each crop runtime global that has a state mirror in swap_mod's A14-A17 blocks AND no remaining bare-global readers (post-Phase B), it's a candidate.

```bash
grep -nE "SS-GR-CROP A14|A15|A16|A17" src/core/swap_mod.f90 | head -30
```

For each `state%X = legacy_sym` line, audit `legacy_sym` for remaining readers:
```bash
sym=daycrop
hits=$(grep -rnE "\b${sym}\b" src/ --include="*.f90" \
  | grep -v "state%\|config%\|! \|variables.f90\|use variables, only:\|initialize.f90\|config_to_variables\|swap_mod.f90\|integer ::\|real(8) ::\|real(real64) ::\|character" | wc -l)
echo "$sym: $hits non-write hits"
```

- [ ] **Step 2: Drop adapter writes** for symbols with 0 non-write hits.

- [ ] **Step 3: VG → commit:**
```bash
git add src/io/toml/config_to_variables.f90
git commit -m "retire(gr-croprt C1): drop adapter writes for retired crop globals

<list symbols dropped>. config_to_variables shrinks by ~N lines.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

### Task 21: C2 — Drop swap_mod A14-A17 dual-writes for retired globals

**Files:** `src/core/swap_mod.f90`

- [ ] Drop dead `state%X = legacy_sym` lines from A14-A17 blocks.

- [ ] VG → commit.

### Task 22: C3 — Drop initialize.f90 zero-fills

**Files:** `src/core/initialize.f90`

- [ ] Drop zero-fills for retired crop globals.

- [ ] VG → commit.

### Task 23: Phase C close

```bash
git commit --allow-empty -m "chore(gr-croprt): Phase C complete — adapter shrinkage landed

config_to_variables + initialize + swap_mod dual-writes shrunk for
retired crop globals.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Phase D — variables.f90 Deletion

### Task 24: D1 — Retire AgeTracer + macropore from variables.f90

Already done in Tasks 2+3. This task is final verification:
```bash
grep -nE "flAgeTracer|AgeGwl1m|icAge|FlMacropore|iQInTopLatDm" src/core/variables.f90
```
Expected: only tombstones (no active declarations).

- [ ] If anything missed, tombstone + iterative patch.
- [ ] Empty marker commit if nothing to do.

### Task 25: D2 — Retire output switches from variables.f90

**Files:** `src/core/variables.f90`, plus compile-surfaced patches

- [ ] Tombstone retired output switches.
- [ ] Iterative compile + patch.
- [ ] VG → commit.

### Task 26: D3 — Retire crop runtime globals (~50-80)

**Files:** `src/core/variables.f90`, plus compile-surfaced patches

- [ ] **Step 1: Identify retire-able crop globals** (those with state mirrors and no remaining readers after Phase B):

```bash
for sym in daycrop dvs tsum swcrp icrop swend rd rdpot rdm rri rdi rdc ch cf laipot cuptgraz cuptgrazpot HarLosOrm_tot flCropCalendar flCropOutput flCropNut flHarvestDay flCropEmergence flCropReadFile flCropPrep flCropSow flCropGerm flCropHarvest cropstart cropend PrepDelay SowDelay noddrz cumdens albedo rsc wlv wlvpot wst wstpot wrt wrtpot wso wsopot tagp tagppot tagpt tagptpot cwdm cwdmpot pgass pgasspot dwlv dwlvpot dwst dwstpot dwrt dwrtpot dwso dwlvCrop dwlvSoil plossdm lossdm swbulb wbl wblpot dwbl dwblpot plwt plwti fco2amax fco2eff fco2tra seqgrazmow seqgrazmowpot dateharvest mowrest cropstartpot cropstartact cropendpot cropendact swpotrelmf relmf cftb chtb cfeic cfeictb; do
  hits=$(grep -rnE "\b${sym}\b" src/ --include="*.f90" \
    | grep -v "state%\|config%\|! \|variables.f90\|use variables, only:\|initialize.f90\|config_to_variables\|swap_mod.f90\|integer ::\|real(8) ::\|real(real64) ::\|character\|logical ::" | wc -l)
  if [ "$hits" -eq "0" ]; then echo "RETIRE: $sym"; else echo "KEEP: $sym ($hits hits)"; fi
done
```

- [ ] **Step 2: Tombstone each "RETIRE" symbol** in variables.f90.

- [ ] **Step 3: Iterative build + patch** any compile errors.

- [ ] **Step 4: VG → commit (may split across multiple commits if surface is large).**

### Task 27: D4 — Residual W-globals catch-all

Anything missed from prior tasks + GR-FINAL D1.

- [ ] Audit + tombstone + VG → commit.

### Task 28: D5 — initialize.f90 final zero-fill cleanup

```bash
grep -nE "^[[:space:]]*\b<sym>\b\s*=" src/core/initialize.f90 | head -30
```

For each retired symbol, drop the zero-fill.

- [ ] VG → commit.

### Task 29: Phase D close marker

- [ ] Sanity grep: `wc -l src/core/variables.f90` — measure shrinkage.

```bash
git commit --allow-empty -m "chore(gr-croprt): Phase D complete — variables.f90 shrunk significantly

variables.f90: 1418L → <new>L. ~<N> globals retired across D1-D5.
check-full 5/5.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Phase E — Close + Roadmap Update

### Task 30: E1 — Final verification

- [ ] `rm -rf builddir && pixi run build-linux`
- [ ] `pixi run test-pfunit`
- [ ] `pixi run check-full` — 5/5
- [ ] BMI + cffi-demo (optional)

### Task 31: E2 — Update roadmap

Edit `docs/superpowers/plans/2026-05-13-globals-retirement-roadmap.md`. Document residuals.

```bash
git add docs/superpowers/plans/2026-05-13-globals-retirement-roadmap.md
git commit -m "docs(roadmap): GR-CROPRT close — residual deferrals documented

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

### Task 32: E3 — Arc-complete marker

```bash
git commit --allow-empty -m "chore(gr-croprt): GR-CROPRT complete — crop runtime + cheap wins

check-full 5/5 byte-for-byte. <document final numbers>.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Plan Self-Review

**Spec coverage:** All 5 phases covered. Tasks 1-32 map to spec sections.

**Placeholder scan:** Audit-driven steps acknowledged where exact symbol lists depend on runtime inspection. Templates provided for substitution patterns.

**Type consistency:** state%crop%common new fields consistent across A4 + A.5 + B + D3 references.

**Known soft spots:**
- cropgrowth Task 9-18 are heavy lifting; some may need controller assistance
- Phase B Task 17 (update_rootdistribution) needs cumdens allocation to match legacy size
- A4 cumdens allocation needs numnod from mesh — verify pattern

---

## Execution

Subagent-driven; controller intervenes inline on infrastructure gaps.
