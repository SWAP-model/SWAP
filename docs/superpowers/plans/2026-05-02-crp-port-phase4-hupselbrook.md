# `.crp` port — Phase 4 (hupselbrook integration + fallback teardown) — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:executing-plans` to
> implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Author three `.crp.toml` files for case 1 (hupselbrook), run the integration
test across all three crop types in one rotation, remove the legacy `if/else` dispatch
fallback from `cropgrowth.f90` for all three types, smoke-test without `.crp` files,
and commit the final cleanup.

**Architecture:** All schema / parser / init-module work is complete (Phases 1–3). Phase
4 is authoring + wiring + teardown only. Per ADR 0016, `crop_config_global` and
`rotation_loaded(:)` remain after Phase 4 (teardown deferred to the config-passing
follow-on spec). Per ADR 0015, the legacy `read*` subs in `readswap.f90` stay alive as
parity-test fixtures.

**Tech stack:** Fortran 2008 (gfortran 13), pFUnit, Meson + Pixi, tomlf, git submodule.

**Spec:** `docs/superpowers/specs/2026-05-02-crp-port-phase4-hupselbrook-design.md`

**ADRs:** `docs/adr/0015-strangler-narrow-scope-stub-errors.md`,
`docs/adr/0016-per-rotation-crop-config-cache.md`

---

## Submodule discipline

`tests/swap-cases/` is a git submodule. Inner-commit + outer-bump pair is
non-negotiable; never one without the other. Affected tasks: Tasks 2, 3, 4 (one
submodule pair commit per .crp.toml authoring) and Task 11 (one submodule pair commit
for all three deletions).

When committing in the submodule, use file-scoped `git commit` (name the specific files)
to avoid bundling pre-existing dirty state from other case directories.

---

## File map

**Source — modify:**
- `src/crop/cropgrowth.f90` — remove the `if/else` dispatch guards at the three
  task=1 callsites (Tasks 7, 8, 9)

**Tests — modify:**
- `tests/unit/io/toml/test_load_swap_config.pf` — add case-1 rotation-loaded assertion
  (Task 6)

**Test fixture — modify (in submodule):**
- `tests/swap-cases/toml/1.hupselbrook/maizes.crp.toml` — expand skeletal 30-line file
  to 1:1 (Task 2)
- `tests/swap-cases/toml/1.hupselbrook/potatod.crp.toml` — verify + gap-fill (Task 3)
- `tests/swap-cases/toml/1.hupselbrook/grassd.crp.toml` — expand 41-line stub to 1:1
  (Task 4)
- `tests/swap-cases/toml/1.hupselbrook/maizes.crp` — DELETE (Task 11)
- `tests/swap-cases/toml/1.hupselbrook/potatod.crp` — DELETE (Task 11)
- `tests/swap-cases/toml/1.hupselbrook/grassd.crp` — DELETE (Task 11)

**Docs — modify:**
- `docs/csv-companion-files.md` — note Phase 4 completion + migration history (Task 12)

---

## Pre-flight

Confirm the starting state. **Do not proceed if any check fails.**

```bash
cd /home/zawadzkim/Code/swap
git log --oneline -5
pixi run test-pfunit 2>&1 | tail -5
pixi run regression 2>&1 | tail -10
```

Expected: `Ok: N, Fail: 0` for unit tests; `5 passed, 0 failed` for regression.

In addition, the following pFUnit assertions (authored by Phases 1–3) must exist and
pass:

```bash
grep -r "test_load_case6_grass_crp_toml_populated\|test_load_case5_potatod\|test_load_case4_grassd\|test_load_case2_grassd" tests/unit/io/toml/test_load_swap_config.pf
```

Expected: at least four matching subroutine names. If any are missing, Phases 1–3 are
not fully landed — stop and escalate.

Also confirm the three fallback dispatch blocks exist in `cropgrowth.f90`:

```bash
grep -n "readcropfixed\|readwofost\|readgrass" /home/zawadzkim/Code/swap/src/crop/cropgrowth.f90
```

Expected: exactly three callsites (one per crop type, at task=1 blocks), all inside
`else` branches. If these are already gone, Phases 1–3 may have done partial teardown
— read the current state before proceeding.

---

## Task 1: Read and understand the current state

Before authoring any TOML files, read the current state of the three `.crp.toml` stubs
and the three legacy `.crp` files to understand what's already present and what's
missing.

**Files:**
- Read: `tests/swap-cases/toml/1.hupselbrook/maizes.crp.toml` (current ~30 lines)
- Read: `tests/swap-cases/toml/1.hupselbrook/potatod.crp.toml` (current ~268 lines)
- Read: `tests/swap-cases/toml/1.hupselbrook/grassd.crp.toml` (current ~41 lines)
- Read: `tests/swap-cases/1.hupselbrook/maizes.crp` (legacy, type 1)
- Read: `tests/swap-cases/1.hupselbrook/potatod.crp` (legacy, type 2, 704 lines)
- Read: `tests/swap-cases/1.hupselbrook/grassd.crp` (legacy, type 3)
- Read: `src/crop/cropfixed_init.f90` (Phase 1 init module) — to understand which
  fields `cropfixed_init_from_config` reads
- Read: `src/crop/cropwofost_init.f90` (Phase 2 init module) — same purpose
- Read: `src/crop/cropgrass_init.f90` (Phase 3 init module) — same purpose
- Read: `src/crop/cropgrowth.f90:368-470` (cropfixed task=1 block — current dispatch)
- Read: `src/crop/cropgrowth.f90:979-1000` (wofost task=1 block — current dispatch)
- Read: `src/crop/cropgrowth.f90:2062-2090` (grass task=1 block — current dispatch)

- [ ] **Step 1: Read the three `.crp.toml` stubs and legacy `.crp` files**

For each pair, produce a gap list: fields present in the legacy `.crp` but absent from
the `.crp.toml` stub. Use the legacy `.crp` as the ground truth.

- [ ] **Step 2: Read the three init modules**

Note the exact signature of each `*_init_from_config` sub. You will need the exact
argument list for Task 7/8/9 when writing the post-removal calls.

```bash
grep -n "subroutine cropfixed_init_from_config\|subroutine cropwofost_init_from_config\|subroutine cropgrass_init_from_config" \
  /home/zawadzkim/Code/swap/src/crop/cropfixed_init.f90 \
  /home/zawadzkim/Code/swap/src/crop/cropwofost_init.f90 \
  /home/zawadzkim/Code/swap/src/crop/cropgrass_init.f90
```

- [ ] **Step 3: Read the current dispatch blocks**

In `cropgrowth.f90`, find and record the exact current form of the three `if/else`
dispatch blocks (the `if (associated(crop_config_global) ...` guards). These are your
BEFORE state for Tasks 7–9.

```bash
grep -n "associated(crop_config_global)\|readcropfixed\|readwofost\|readgrass" \
  /home/zawadzkim/Code/swap/src/crop/cropgrowth.f90
```

No commit for this task — it is reading only.

---

## Task 2: Author `maizes.crp.toml` 1:1

Translate `tests/swap-cases/1.hupselbrook/maizes.crp` (type 1, cropfixed, 463 lines)
into a fully-authored TOML file. This is a SUBMODULE PAIR COMMIT.

**Reference:** Phase 1's `grass.crp.toml` (case 6) and `cropfixed_config.f90` schema
as ground truth for section names and field names.

**Files:**
- Modify (in submodule): `tests/swap-cases/toml/1.hupselbrook/maizes.crp.toml`

Key differences from case-6 `grass.crp`:
- `lcc = 168` (not 366)
- `swcf = 2` (crop height, not crop factor) — include `chtb` table (7 pairs) + `albedo`,
  `rsc`, `rsw`
- `swgc = 1` — `gctb` table has 7 pairs (DVS 0.0–2.0, LAI values)
- `swoxygen = 1` — include `hlim1`, `hlim2u`, `hlim2l` fields
- `rdtb` has 6 pairs (DVS 0.0–2.0, root depth)
- Full TOML content to author (value-for-value from the legacy .crp):

```toml
# Hupselbrook maize crop (type 1, fixed) — 1:1 reproduction of
# tests/swap-cases/1.hupselbrook/maizes.crp

[preparation]
swprep = 0
swsow  = 0
swgerm = 0
dvsend = 3.0
swharv = 0

[phenology]
idev = 1
lcc  = 168
tsumea = 1050.0
tsumam = 1000.0
tbase  = 0.0

[light]
kdif = 0.6
kdir = 0.75

[lai]
swgc = 1
gctb = [
  [0.0,  0.05],
  [0.3,  0.14],
  [0.5,  0.61],
  [0.7,  4.10],
  [1.0,  5.00],
  [1.4,  5.80],
  [2.0,  5.20],
]

[crop_factor]
swcf   = 2
albedo = 0.23
rsc    = 61.0
rsw    = 0.0
cftb = [
  [0.0, 0.8],
  [0.3, 0.8],
  [0.5, 0.9],
  [0.7, 1.0],
  [1.0, 1.1],
  [1.4, 1.2],
  [2.0, 1.2],
]
chtb = [
  [0.0,   1.0],
  [0.3,  15.0],
  [0.5,  40.0],
  [0.7, 140.0],
  [1.0, 170.0],
  [1.4, 180.0],
  [2.0, 175.0],
]

[root]
swrd     = 1
swdmi2rd = 1
swrdc    = 0
rdi      = 10.0
rri      = 2.2
rdc      = 100.0
wrtmax   = 3000.0
rdtb = [
  [0.0,   5.0],
  [0.3,  20.0],
  [0.5,  50.0],
  [0.7,  80.0],
  [1.0,  90.0],
  [2.0, 100.0],
]
rdctb = [
  [0.0, 1.0],
  [1.0, 0.0],
]

[oxygen_stress]
swoxygen    = 1
swwrtnonox  = 0
aeratecrit  = 0.5
hlim1       = -15.0
hlim2u      = -30.0
hlim2l      = -30.0
q10_microbial           = 2.8
specific_resp_humus     = 0.0016
srl                     = 151375.0
swrootradius            = 2
dry_mat_cont_roots      = 0.075
air_filled_root_por     = 0.05
spec_weight_root_tissue = 1000.0
var_a                   = 0.000000000418
root_radiusO2           = 0.00015
q10_root                = 2.0
f_senes                 = 1.0
c_mroot                 = 0.016

[drought_stress]
swdrought = 1
hlim3h    = -325.0
hlim3l    = -600.0
hlim4     = -8000.0
adcrh     = 0.5
adcrl     = 0.1

[salinity]
swsalinity = 0
saltmax    = 3.0
saltslope  = 0.1
salthead   = 624.0

[compensate]
swcompensate = 0
swstressor   = 3
alphacrit    = 1.0
dcritrtz     = 5.0

[interception]
swinter = 1
cofab   = 0.25

[scheduling]
schedule = 0
```

- [ ] **Step 1: Write the full maizes.crp.toml**

Replace the skeletal file with the above. Verify every value matches the legacy
`maizes.crp` exactly. Pay special attention to:
- `lcc = 168` (not 366)
- `gctb`: 7 rows, LAI values are 0.05 / 0.14 / 0.61 / 4.10 / 5.00 / 5.80 / 5.20
- `chtb`: 7 rows, heights 1.0 / 15.0 / 40.0 / 140.0 / 170.0 / 180.0 / 175.0
- `rsc = 61.0` (from the SWCF=2 block)
- `swoxygen = 1` with hlim1=-15.0, hlim2u=-30.0, hlim2l=-30.0

- [ ] **Step 2: Run the regression to confirm no breakage**

```bash
cd /home/zawadzkim/Code/swap
pixi run regression 2>&1 | tail -10
```

Expected: 5/5 green. Case 1 may still fail (if the schema wasn't ready) or pass with
the legacy fallback. Either is acceptable at this stage.

- [ ] **Step 3: Submodule inner-commit**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
git commit toml/1.hupselbrook/maizes.crp.toml \
  -m "feat(toml/1.hupselbrook): maizes.crp.toml 1:1 with legacy maizes.crp

Replaces skeletal file with full 1:1 schema coverage of maizes.crp.
Active switches: idev=1, lcc=168, swgc=1, swcf=2 (chtb + albedo/rsc/rsw),
swrd=1 (rdtb 6-pairs), swoxygen=1 (hlim1/hlim2u/hlim2l), swdrought=1,
swsalinity=0, swcompensate=0, swinter=1, schedule=0."
```

- [ ] **Step 4: Outer-repo bump commit**

```bash
cd /home/zawadzkim/Code/swap
git add tests/swap-cases
git commit -m "feat(case1): author maizes.crp.toml 1:1; bump submodule"
```

---

## Task 3: Verify and complete `potatod.crp.toml`

The existing `potatod.crp.toml` (268 lines) was prepared in a prior task. Verify it is
a faithful 1:1 reproduction of `potatod.crp` (704 lines). This is a SUBMODULE PAIR
COMMIT only if gaps are found and filled.

**Files:**
- Read: `tests/swap-cases/1.hupselbrook/potatod.crp`
- Read: `tests/swap-cases/toml/1.hupselbrook/potatod.crp.toml`

- [ ] **Step 1: Field-by-field diff**

For each section in `potatod.crp`, verify the corresponding section exists in
`potatod.crp.toml` and all field values match.

Sections to check:
- `[preparation]` — swprep=0, swsow=0
- `[germination]` — **swgerm=2** (key difference from case 6): tsumemeopt=170.0,
  tbasem=3.0, teffmx=18.0, hdrygerm=-500.0, hwetgerm=-100.0, zgerm=-10.0, agerm=203.0
- `[harvest]` — dvsend=2.0, swharv=0
- `[crop_factor]` — swcf=2, cftb, chtb, albedo=0.19, rsc=207.0, rsw=0.0
- `[phenology]` — idsl=0, tsumea=150.0, tsumam=1550.0, dlo=14.0, dlc=8.0,
  vernsat=70.0, vernbase=14.0, verndvs=0.3, dtsmtb (4 rows), verntb (6 rows)
- `[initial]` — tdwi=75.0, laiem=0.0589, rgrlai=0.012
- `[green_area]` — spa=0.0, ssa=0.0, span=37.0, tbase=2.0, slatb (3 rows)
- `[assimilation]` — kdif=1.0, kdir=0.75, eff=0.45, amaxtb (3 rows), tmpftb (7 rows),
  tmnftb (2 rows)
- `[conversion]` — cvl=0.72, cvo=0.85, cvr=0.72, cvs=0.69
- `[respiration]` — q10=2.0, rml=0.03, rmo=0.0045, rmr=0.01, rms=0.015, rfsetb (2 rows)
- `[partitioning]` — frtb (4 rows), fltb (5 rows), fstb (5 rows), fotb (5 rows)
- `[death]` — perdl=0.03, rdrrtb (4 rows), rdrstb (4 rows)
- `[root]` — swrd=2, rdi=10.0, rri=1.2, rdc=50.0, swdmi2rd=1, wrtmax=3000.0,
  rdtb (3 rows), rlwtb (2 rows), rdctb (2 rows)
- `[oxygen_stress]` — swoxygen=1, swwrtnonox=1, aeratecrit=0.5, hlim1=-10.0,
  hlim2u=-25.0, hlim2l=-25.0, all swoxygen=2 subordinate fields
- `[drought_stress]` — swdrought=1, hlim3h=-300.0, hlim3l=-500.0, hlim4=-10000.0,
  adcrh=0.5, adcrl=0.1
- `[salinity]` — swsalinity=0, saltmax, saltslope, salthead
- `[compensate]` — swcompensate=0, swstressor, alphacrit, dcritrtz
- `[interception]` — swinter=1, cofab=0.25, gashtb (2 rows)
- `[co2]` — swco2=0, atmofil, co2amaxtb (5 rows), co2efftb (5 rows), co2tratb (5 rows)
- `[management]` — fraharlosorm_lv=0.2, fraharlosorm_st=0.1, fraharlosorm_so=0.01,
  fradeceasedlvtosoil=0.3, swpotrelmf=2, relmf=0.9
- `[scheduling]` — schedule=0

- [ ] **Step 2: If gaps found, fill them**

For any section or field found missing from `potatod.crp.toml`, add it. If the file is
verified complete as-is, no edit is needed.

- [ ] **Step 3: If any edits were made — run regression**

```bash
cd /home/zawadzkim/Code/swap
pixi run regression 2>&1 | tail -10
```

Expected: 5/5 green.

- [ ] **Step 4: If any edits were made — submodule pair commit**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
git commit toml/1.hupselbrook/potatod.crp.toml \
  -m "fix(toml/1.hupselbrook): potatod.crp.toml gap-fill for 1:1 coverage"
cd /home/zawadzkim/Code/swap
git add tests/swap-cases
git commit -m "fix(case1): potatod.crp.toml gap-fill; bump submodule"
```

If the file was already complete, skip the inner/outer commits. Note in the Task 3
completion log that no edits were required.

---

## Task 4: Author `grassd.crp.toml` 1:1

Expand the current 41-line `grassd.crp.toml` stub to full 1:1 coverage of
`tests/swap-cases/1.hupselbrook/grassd.crp` (type 3, cropgrass, 493 lines). This is a
SUBMODULE PAIR COMMIT.

**Reference:** Phase 3's cropgrass `.crp.toml` files (cases 4 and 2) as schema ground
truth for section and field names.

**Files:**
- Modify (in submodule): `tests/swap-cases/toml/1.hupselbrook/grassd.crp.toml`

The current stub has `[light]`, `[root]` (partial), `[water_stress]`, `[interception]`,
`[mowing]`, `[grazing]`. The following sections must be added:

```toml
# Hupselbrook grass crop (type 3, WOFOST grass) — 1:1 reproduction of
# tests/swap-cases/1.hupselbrook/grassd.crp

[crop_factor]
swcf   = 2
albedo = 0.23
rsc    = 100.0
rsw    = 0.0
# DNR / CH / CF table (SWCF=2: use CH column; CF column present but ignored by runtime)
chtb = [
  [  0.0, 12.0],
  [180.0, 12.0],
  [366.0, 12.0],
]
cftb = [
  [  0.0, 1.0],
  [180.0, 1.0],
  [366.0, 1.0],
]

[initial]
tdwi   = 1000.0
laiem  = 0.63
rgrlai = 0.007
swtsum = 1
tsumtemp  = 8.0
tsumdepth = 10.0
tsumtime  = 3

[green_area]
ssa   = 0.0004
span  = 30.0
tbase = 0.0
slatb = [
  [  1.0, 0.0015],
  [ 80.0, 0.0015],
  [300.0, 0.0020],
  [366.0, 0.0020],
]

[assimilation]
kdif = 0.60
kdir = 0.75
eff  = 0.50
amaxtb = [
  [  1.0, 40.0],
  [ 95.0, 40.0],
  [200.0, 35.0],
  [275.0, 25.0],
  [366.0, 25.0],
]
tmpftb = [
  [ 0.0, 0.00],
  [ 5.0, 0.70],
  [15.0, 1.00],
  [25.0, 1.00],
  [40.0, 0.00],
]
tmnftb = [
  [0.0, 0.0],
  [4.0, 1.0],
]

[conversion]
cvl = 0.685
cvr = 0.694
cvs = 0.662

[respiration]
q10 = 2.0
rml = 0.030
rmr = 0.015
rms = 0.015
rfsetb = [
  [  1.0, 1.0],
  [366.0, 1.0],
]

[partitioning]
frtb = [
  [  1.0, 0.30],
  [366.0, 0.30],
]
fltb = [
  [  1.0, 0.60],
  [366.0, 0.60],
]
fstb = [
  [  1.0, 0.40],
  [366.0, 0.40],
]

[death]
perdl = 0.050
rdrrtb = [
  [  1.0, 0.00],
  [180.0, 0.02],
  [366.0, 0.02],
]
rdrstb = [
  [  1.0, 0.00],
  [180.0, 0.02],
  [366.0, 0.02],
]

[root]
swrd   = 3
swrdc  = 0
rdi    = 20.0
rri    = 0.25
rdc    = 40.0
wrtmax = 3000.0
rdtb = [
  [  1.0, 20.0],
  [180.0, 40.0],
  [366.0, 40.0],
]
rlwtb = [
  [ 300.0, 20.0],
  [2500.0, 40.0],
]
rdctb = [
  [0.0, 1.0],
  [1.0, 0.0],
]

[oxygen_stress]
swoxygen = 1
hlim1    =   0.0
hlim2u   =   1.0
hlim2l   =  -1.0

[drought_stress]
swdrought = 1
# Legacy swjarvis=4 is dropped per user direction (2026-05-02).
# Translate to modern swcompensate. The closest equivalent of legacy
# swjarvis=4 ("compensate drought, wet, salt, frost") is swcompensate=1
# (Jarvis, all stressors per swstressor=1 default). If regression drifts,
# stub-error this rotation and run case 1 via the legacy executable.
swcompensate = 1
alphacrit    = 0.7
hlim3h    = -200.0
hlim3l    = -800.0
hlim4     = -8000.0
adcrh     =    0.5
adcrl     =    0.1

[salinity]
swsalinity = 0

[co2]
swco2   = 0
atmofil = "atmospheric"

[management]
swpotrelmf = 1
relmf      = 0.90

[mowing]
swharvest     = 1
swdmmow       = 2
dmharvest     = 4200.0
daylastharvest = 289
dmlastharvest = 2700.0
maxdaymow     = 42
swlossmow     = 0
mowrest       = 700.0
dmmowtb = [
  [120.0, 4700.0],
  [152.0, 3700.0],
  [182.0, 3200.0],
  [213.0, 2700.0],
  [366.0, 2700.0],
]
dmmowdelay = [
  [   0.0, 2],
  [2000.0, 3],
  [4000.0, 4],
]

[grazing]
swgraz   = 0
swdmgrz  = 2
maxdaygrz = 28
swlossgrz = 0
tagprest  = 700.0
dewrest   = 850.0
dmgrazing = 2400.0
dmgrztb = [
  [152.0, 2400.0],
  [244.0, 1800.0],
  [366.0, 1800.0],
]
lsdatb = [
  [21.25, 4.0, 16.0, 4.0],
]

[scheduling]
schedule = 0
```

- [ ] **Step 1: Read Phase 3's cropgrass `.crp.toml` for section name conventions**

```bash
cat /home/zawadzkim/Code/swap/tests/swap-cases/toml/4.oxygenstress/grassd.crp.toml 2>/dev/null | head -80
cat /home/zawadzkim/Code/swap/tests/swap-cases/toml/2.grassgrowth/grassd.crp.toml 2>/dev/null | head -80
```

Verify that the section names and field names used above match Phase 3's conventions.
If they differ, adjust the above TOML content to match Phase 3's actual schema.

- [ ] **Step 2: Write the full grassd.crp.toml**

Replace the current 41-line stub with the full content above (adjusted for any
conventions discovered in Step 1). Preserve the existing `[mowing]` and `[grazing]`
sections if they already have correct values.

**SWJARVIS dropped (2026-05-02 user direction).** The TOML pipeline does not
support the deprecated legacy `swjarvis` key. Use modern `swcompensate` instead.
Translate legacy `SWJARVIS=4` → `swcompensate = 1` with `swstressor = 1`
(Jarvis-compensate-all). If regression drifts unacceptably from the legacy
reference, stub-error this rotation and run case 1 via the legacy executable
rather than re-introducing `swjarvis` support.

- [ ] **Step 3: Run the regression**

```bash
cd /home/zawadzkim/Code/swap
pixi run regression 2>&1 | tail -10
```

Expected: 5/5 green. If case 1's grass rotation diverges from the legacy
reference because of the `swjarvis` → `swcompensate` translation, accept
the divergence (and document it) or stub-error the rotation type-3 in
case 1's swap.toml. Do NOT re-introduce `swjarvis` support to fix this.

- [ ] **Step 4: Submodule inner-commit**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
git commit toml/1.hupselbrook/grassd.crp.toml \
  -m "feat(toml/1.hupselbrook): grassd.crp.toml 1:1 with legacy grassd.crp

Expands 41-line stub to full 1:1 schema coverage of grassd.crp.
Active switches: swcf=2, tdwi=1000.0, swtsum=1, swrd=3 (rlwtb-biomass),
swrdc=0, swoxygen=1, swdrought=1, swcompensate=1 (translated from
legacy swjarvis=4 per 2026-05-02 user direction), swsalinity=0,
swco2=0, swharvest=1+swdmmow=2, swgraz=0, schedule=0."
```

- [ ] **Step 5: Outer-repo bump commit**

```bash
cd /home/zawadzkim/Code/swap
git add tests/swap-cases
git commit -m "feat(case1): author grassd.crp.toml 1:1; bump submodule"
```

---

## Task 5: Regression test case 1 with all three `.crp.toml` files complete

Run the full regression suite with all three `.crp.toml` files in place and the legacy
`.crp` files still present. This is the integration test.

- [ ] **Step 1: Run regression**

```bash
cd /home/zawadzkim/Code/swap
pixi run regression 2>&1 | tail -20
```

Expected: 5/5 green. If case 1 fails, proceed to Step 2.

- [ ] **Step 2 (if needed): Diagnose case 1 failure**

Case 1 failure at this point is a mismatch between what the init modules write to module
globals and what the simulation expects. Common causes:

a. **Value mismatch in a table:** e.g. `gctb(n)` value wrong — recheck the TOML against
   the legacy `.crp`.
b. **Missing field causes default-initialization mismatch:** e.g. a field that
   `*_init_from_config` reads but that isn't in the TOML, causing a 0.0 default to be
   used instead of the legacy value.
c. **Rotation-ordering issue:** if `icrop` index starts at 1 for the first active
   rotation but `rotation_loaded(1)` corresponds to the 2002 start date while the 2003
   rotation is at index 2.

Diagnose by running with verbose output, or by temporarily reverting to the legacy
fallback on one rotation type to narrow down which type is causing the mismatch.

- [ ] **Step 3 (if needed): Fix and re-run**

Fix the root cause — either in the TOML file (Task 2/3/4 submodule pair commit) or in
the init module (escalate to Phase 1/2/3 maintainer). Re-run until 5/5 green.

---

## Task 6: Write and verify the pFUnit case-1 rotation-loaded assertion

Add the integration test assertion to `test_load_swap_config.pf`.

**Files:**
- Modify: `tests/unit/io/toml/test_load_swap_config.pf`

- [ ] **Step 1: Append the test subroutine**

Append the following to `tests/unit/io/toml/test_load_swap_config.pf`:

```fortran
@test
subroutine test_load_case1_hupselbrook_all_rotations_loaded()
   use funit
   use error_mod,        only: error_collection_t
   use swap_config_mod,  only: swap_config_t
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

**Note:** Adjust the field names (`rotation_fixed`, `rotation_wofost`, `rotation_grass`,
`rotation_loaded`) to match whatever names Phase 1 actually implemented on
`crop_config_t`. The plan header notes that Phase 1 used `rotation_fixed(:)` and
`rotation_loaded(:)` (not `rotation_cropfixed(:)` and `populated`). Phase 2/3 used
analogous names.

- [ ] **Step 2: Run unit tests**

```bash
cd /home/zawadzkim/Code/swap
pixi run test-pfunit 2>&1 | tail -10
```

Expected: `Ok: N, Fail: 0` including the new test.

- [ ] **Step 3: Commit**

```bash
git add tests/unit/io/toml/test_load_swap_config.pf
git commit -m "test(toml): assert case-1 hupselbrook all three rotation slots loaded"
```

---

## Task 7: Remove the type-1 (cropfixed) fallback dispatch

Remove the `else call readcropfixed(...)` branch from `cropgrowth.f90:subroutine cropfixed:task=1`.

**Files:**
- Modify: `src/crop/cropgrowth.f90`

**Prerequisite:** Task 1's Step 3 recorded the exact current form of the type-1 dispatch
block. The BEFORE state is that block; the AFTER is the unconditional call shown below.

**BEFORE** (current state after Phase 1, found at the `! --- read crop data` comment
in `subroutine cropfixed, case(1)`, around line 397):

```fortran
! --- read crop data
      if (associated(crop_config_global) .and. &
          allocated(crop_config_global%rotation_loaded) .and. &
          crop_config_global%rotation_loaded(icrop)) then
         call cropfixed_init_from_config(crop_config_global%rotation_fixed(icrop), icrop)
      else
         call readcropfixed (icrop, cropfil(icrop), lcc, swhydrlift)   ! transitional
      end if
```

**AFTER** (replace the entire `if/else/end if` block with the single unconditional call):

```fortran
! --- read crop data
      call cropfixed_init_from_config(crop_config_global%rotation_fixed(icrop), icrop)
```

- [ ] **Step 1: Locate the exact block in the source**

```bash
grep -n "readcropfixed\|rotation_loaded\|associated(crop_config_global)" \
  /home/zawadzkim/Code/swap/src/crop/cropgrowth.f90 | head -20
```

Note the exact line numbers. Read those lines with the Read tool to confirm the exact
whitespace and continuation style before editing.

- [ ] **Step 2: Apply the edit**

Use the Edit tool to replace the `if/else/end if` block with the unconditional call.
Match the indentation of the surrounding code (legacy SWAP uses 6-space + column-7
indentation).

- [ ] **Step 3: Build**

```bash
cd /home/zawadzkim/Code/swap
pixi run swap 2>&1 | tail -20
```

Expected: clean build.

- [ ] **Step 4: Run unit tests + regression**

```bash
pixi run test-pfunit 2>&1 | tail -5
pixi run regression 2>&1 | tail -10
```

Expected: all green. If case 6 regresses, the init module has a bug (not a dispatch
issue) — escalate to Phase 1 maintainer.

- [ ] **Step 5: Commit**

```bash
git add src/crop/cropgrowth.f90
git commit -m "refactor(crop): remove cropfixed legacy fallback dispatch

Phase 4 of the .crp port (ADR 0016). After Phase 1 ported type-1 rotations
and Phase 4 authored the hupselbrook maizes.crp.toml, every type-1 rotation
in every TOML case has rotation_fixed(icrop) populated. The readcropfixed
fallback is no longer reachable on the TOML path. Replace the if/else guard
with an unconditional cropfixed_init_from_config call.

readcropfixed in readswap.f90 is kept alive as a parity-test fixture per
ADR 0015."
```

---

## Task 8: Remove the type-2 (cropwofost) fallback dispatch

Remove the `else call readwofost(...)` branch from `cropgrowth.f90:subroutine wofost:task=1`.

**Files:**
- Modify: `src/crop/cropgrowth.f90`

**Prerequisite:** Task 1's Step 3 recorded the exact current form of the type-2 dispatch
block (around line 984). The BEFORE state is that block.

**BEFORE** (current state after Phase 2, found at the `! --- read general crop data`
comment in `subroutine wofost, case(1)`, around line 984):

```fortran
! --- read general crop data
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

**AFTER:**

```fortran
! --- read general crop data
      call cropwofost_init_from_config(crop_config_global%rotation_wofost(icrop), icrop, &
                                       swhydrlift, swsoybean, mg, dvsi, dvrmax1, dvrmax2, &
                                       flrfphotoveg, tmaxdvr, tmindvr, toptdvr, popt, pcrt, &
                                       flphenodayl, FraDeceasedLvToSoil)
```

**IMPORTANT:** The argument list above is based on the legacy `readwofost` call at
`cropgrowth.f90:984`. The actual `cropwofost_init_from_config` signature may differ
slightly — Phase 2 may have added or reordered arguments. Use the argument list from
`src/crop/cropwofost_init.f90` (read in Task 1 Step 2) as the authoritative source.

- [ ] **Step 1: Locate the exact block**

```bash
grep -n "readwofost\|rotation_wofost\|cropwofost_init_from_config" \
  /home/zawadzkim/Code/swap/src/crop/cropgrowth.f90
```

Read those lines to confirm the exact form.

- [ ] **Step 2: Apply the edit**

Replace the `if/else/end if` block with the unconditional call. Use the exact argument
list from `src/crop/cropwofost_init.f90`.

- [ ] **Step 3: Build + test + regression**

```bash
pixi run swap 2>&1 | tail -20
pixi run test-pfunit 2>&1 | tail -5
pixi run regression 2>&1 | tail -10
```

Expected: all green. Case 5 (salinitystress) and case 1's potato rotation must both
pass.

- [ ] **Step 4: Commit**

```bash
git add src/crop/cropgrowth.f90
git commit -m "refactor(crop): remove cropwofost legacy fallback dispatch

Phase 4 of the .crp port. Every type-2 rotation in every TOML case has
rotation_wofost(icrop) populated after Phases 2+4. readwofost in
readswap.f90 kept alive as parity-test fixture (ADR 0015)."
```

---

## Task 9: Remove the type-3 (cropgrass) fallback dispatch

Remove the `else call readgrass(...)` branch from `cropgrowth.f90:subroutine grass:task=1`.

**Files:**
- Modify: `src/crop/cropgrowth.f90`

**BEFORE** (current state after Phase 3, found at the `! --- read grass input data`
comment in `subroutine grass, case(1)`, around line 2068):

```fortran
! --- read grass input data
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

**AFTER:**

```fortran
! --- read grass input data
      call cropgrass_init_from_config(crop_config_global%rotation_grass(icrop), icrop, &
                                      swharvest, dmharvest, daylastharvest, dmlastharvest, &
                                      swdmmow, maxdaymow, swlossmow, swlossgrz, swdmgrz, &
                                      maxdaygrz, dmgrazing, LSDb, tagprest, swhydrlift)
```

**IMPORTANT:** The `readgrass` call signature at `cropgrowth.f90:2068` is the ground
truth. Use the actual `cropgrass_init_from_config` signature from
`src/crop/cropgrass_init.f90` (read in Task 1 Step 2).

- [ ] **Step 1: Locate the exact block**

```bash
grep -n "readgrass\|rotation_grass\|cropgrass_init_from_config" \
  /home/zawadzkim/Code/swap/src/crop/cropgrowth.f90
```

- [ ] **Step 2: Apply the edit**

Replace the `if/else/end if` block with the unconditional call.

- [ ] **Step 3: Build + test + regression**

```bash
pixi run swap 2>&1 | tail -20
pixi run test-pfunit 2>&1 | tail -5
pixi run regression 2>&1 | tail -10
```

Expected: all green. Cases 4 (oxygenstress), 2 (grassgrowth), and 1 (hupselbrook
grass rotation) must all pass.

- [ ] **Step 4: Commit**

```bash
git add src/crop/cropgrowth.f90
git commit -m "refactor(crop): remove cropgrass legacy fallback dispatch

Phase 4 of the .crp port. Every type-3 rotation in every TOML case has
rotation_grass(icrop) populated after Phases 3+4. readgrass in
readswap.f90 kept alive as parity-test fixture (ADR 0015)."
```

---

## Task 10: Smoke test — run without `.crp` files in the TOML dir

Before permanently deleting the legacy `.crp` files, rename them to `.crp.disabled`
and run the regression to confirm the TOML path is fully decoupled.

- [ ] **Step 1: Rename the three `.crp` files**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases/toml/1.hupselbrook
mv maizes.crp  maizes.crp.disabled
mv potatod.crp potatod.crp.disabled
mv grassd.crp  grassd.crp.disabled
```

- [ ] **Step 2: Run regression**

```bash
cd /home/zawadzkim/Code/swap
pixi run regression 2>&1 | tail -10
```

Expected: 5/5 green. The runtime must not open the legacy `.crp` files at all.

If regression fails: the dispatch removal has a bug (a callsite that was missed, or an
init module that falls through to a file-open path). Diagnose and fix Tasks 7/8/9 before
proceeding.

- [ ] **Step 3: Restore the `.crp` files**

Do NOT commit the disabled state. Restore immediately:

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases/toml/1.hupselbrook
mv maizes.crp.disabled  maizes.crp
mv potatod.crp.disabled potatod.crp
mv grassd.crp.disabled  grassd.crp
```

No commit for this task.

---

## Task 11: Delete the three `.crp` files from the TOML dir

With the smoke test confirming the runtime is fully decoupled, permanently delete the
three legacy `.crp` files from the TOML case directory. This is a SUBMODULE PAIR COMMIT.

**Files to delete (in submodule):**
- `tests/swap-cases/toml/1.hupselbrook/maizes.crp`
- `tests/swap-cases/toml/1.hupselbrook/potatod.crp`
- `tests/swap-cases/toml/1.hupselbrook/grassd.crp`

**Files left intact (in legacy case dir — do not touch):**
- `tests/swap-cases/1.hupselbrook/maizes.crp`
- `tests/swap-cases/1.hupselbrook/potatod.crp`
- `tests/swap-cases/1.hupselbrook/grassd.crp`

- [ ] **Step 1: Delete the files in the submodule**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
git rm toml/1.hupselbrook/maizes.crp
git rm toml/1.hupselbrook/potatod.crp
git rm toml/1.hupselbrook/grassd.crp
```

- [ ] **Step 2: Run regression one final time**

```bash
cd /home/zawadzkim/Code/swap
pixi run regression 2>&1 | tail -10
```

Expected: 5/5 green (same as smoke test, but now permanent).

- [ ] **Step 3: Acceptance gate check**

```bash
# Confirm files are gone
ls tests/swap-cases/toml/1.hupselbrook/*.crp 2>&1
# Expected: "No such file or directory" or an empty listing

# Confirm legacy reader calls are absent from the runtime path
git grep "call readcropfixed\|call readwofost\|call readgrass" src/
# Expected: zero results

# Confirm legacy readers still exist as fixtures
git grep "subroutine readcropfixed\|subroutine readwofost\|subroutine readgrass" src/
# Expected: three matches (in readswap.f90)
```

- [ ] **Step 4: Submodule inner-commit**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
git commit toml/1.hupselbrook/maizes.crp toml/1.hupselbrook/potatod.crp toml/1.hupselbrook/grassd.crp \
  -m "feat(toml/1.hupselbrook): delete legacy .crp files; TOML path is standalone

Phase 4 of the .crp port. After authoring the three .crp.toml files and
removing the fallback dispatch in cropgrowth.f90, the runtime no longer
needs the legacy .crp files in the TOML case directory. The originals in
tests/swap-cases/1.hupselbrook/ are untouched (legacy executable support)."
```

- [ ] **Step 5: Outer-repo bump commit**

```bash
cd /home/zawadzkim/Code/swap
git add tests/swap-cases
git commit -m "feat(case1): delete legacy .crp files from TOML dir; bump submodule

Completes the .crp port for hupselbrook (case 1). All three crop types
(cropfixed/cropwofost/cropgrass) in the multi-type rotation are now
served entirely by the typed-config pipeline. The if/else legacy fallback
is removed; regression 5/5 green without any .crp files in the TOML dirs."
```

---

## Task 12: Docs update

Update the migration history and teardown documentation.

**Files:**
- Modify: `docs/csv-companion-files.md`
- Review (no modification): `docs/adr/0016-per-rotation-crop-config-cache.md`

- [ ] **Step 1: Update `docs/csv-companion-files.md`**

Read the existing migration history section of `docs/csv-companion-files.md`. Find the
`.crp` port entry (if it exists) and extend it with a Phase 4 completion note. If no
`.crp` entry exists, add one.

The note should record:
- Phase 4 completion date
- Three files removed from case 1 TOML dir
- The three fallback dispatch branches removed from `cropgrowth.f90`
- What remains deferred: `crop_config_global` pointer and `rotation_loaded(:)` sentinel
  (trigger: config-passing follow-on spec, ADR 0016 Part C)
- What remains as parity-test fixtures: `readcropfixed`, `readwofost`, `readgrass` in
  `readswap.f90`

- [ ] **Step 2: Verify ADR 0016 teardown table is accurate**

Read `docs/adr/0016-per-rotation-crop-config-cache.md` section "Teardown plan". Verify
that the table entries for:
- "Per-rotation legacy reader fallback" → Trigger: End of Phase 4 ✓ (now done)
- "`crop_config_global` module pointer" → Trigger: Config-passing follow-on spec (still
  pending — verify this entry still reads correctly given the Phase 4 outcome)
- "`populated :: logical` sentinel" → Trigger: Same as `crop_config_global`

If any entry needs a minor status update (e.g. changing "planned" to "done" on the
fallback removal), update the ADR. Do not change any entry's trigger or replacement
unless it is factually incorrect.

- [ ] **Step 3: Commit**

```bash
git add docs/csv-companion-files.md docs/adr/0016-per-rotation-crop-config-cache.md
git commit -m "docs(.crp port): Phase 4 completion — hupselbrook integration + teardown"
```

---

## Final state verification

After Task 12, run the complete acceptance gate:

```bash
cd /home/zawadzkim/Code/swap

# Unit tests
pixi run test-pfunit 2>&1 | tail -5
# Expected: Ok: N, Fail: 0

# Regression
pixi run regression 2>&1 | tail -10
# Expected: 5 passed, 0 failed

# No legacy .crp files in any TOML dir
ls tests/swap-cases/toml/*/**.crp 2>&1
# Expected: zero files

# No legacy reader runtime calls in source
git grep "call readcropfixed\|call readwofost\|call readgrass" src/
# Expected: zero results

# Legacy readers still exist as fixtures
git grep "subroutine readcropfixed\|subroutine readwofost\|subroutine readgrass" src/
# Expected: three matches

# crop_config_global still present (deferred teardown)
git grep "crop_config_global" src/crop/
# Expected: at least one match (the module declaration + its use in cropgrowth.f90)

# rotation_loaded still present on crop_config_t (deferred teardown)
git grep "rotation_loaded" src/
# Expected: at least one match (the field declaration + its use in the loader)
```

All checks green = Phase 4 complete.
