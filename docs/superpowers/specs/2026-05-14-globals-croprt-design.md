---
title: "GR-CROPRT — Crop Runtime Finalization + Cheap Wins"
date: 2026-05-14
status: approved
context: globals-retirement continuation after GR-FINAL partial completion
---

# GR-CROPRT Design

## Goal

Close the largest remaining retirement blockers in two coupled efforts:

1. **Cheap-wins cluster** — retire ADR-documented dead-code globals (AgeTracer, macropore sentinels) and small output control switches. ~30-50 globals retired with low-risk substitutions.

2. **cropgrowth.f90 finalization** — schema-add the ~15-20 missing fields to `state%crop` sub-records, then migrate cropgrowth's 10 use-variables sites (~738+ symbol substitutions). Unblocks ~50-80 crop runtime globals for variables.f90 retirement.

After this arc, variables.f90 shrinks ~50% from current state (1418L → <800L target).

## Architecture

Five phases:

- **Phase A** — Cheap-wins retirements + crop schema additions (additive only).
- **Phase A.5** — Runtime dual-write coverage for new crop fields.
- **Phase B** — cropgrowth.f90 reader migration (10 sites, one commit per site).
- **Phase C** — Drop adapter writes + Phase A.5 dual-write blocks for retired globals.
- **Phase D** — Iterative variables.f90 declaration deletion + initialize.f90 zero-fill cleanup.
- **Phase E** — Close marker + roadmap update for residual deferrals.

## Verification Discipline

Per `feedback_per_task_regression_gate.md` + `feedback_state_schema_clean_rebuild.md`:

- **Every implementer subagent ends with VG:**
  ```bash
  rm -rf builddir && pixi run build-linux \
    && pixi run test-pfunit \
    && pixi run -e test python tests/regression/test_output_regression.py \
         hupselbrook surfacewater salinitystress grassgrowth
  ```
  Expected: clean build, pFUnit 741, regression **4/4 byte-for-byte**.

- **Every phase ends with:** `pixi run check-full` 5/5 byte-for-byte.

- **Controller intervention policy:** subagents do NOT defer due to missing infrastructure. If schema needs extending mid-task, do it inline. If config threading chains through callers, do it inline. Per user directive after GR-ATM scaffolding pattern.

## Phase A — Cheap Wins + Crop Schema Additions

### A1: AgeTracer dead-code retirement

ADR 0032 — body unreachable per `flAgeTracer` guard (always `.false.`).

**Globals retired:** `flAgeTracer`, `AgeGwl1m`, `icAgeBot`, `icAgeDra`, `icAgeRot`, `icAgeSur`, plus any body-only locals in `src/solute/agetracer.f90`.

**Actions:**
- Tombstone declarations in `src/core/variables.f90`
- Drop adapter writes in `src/io/toml/config_to_variables.f90`
- Drop zero-fills in `src/core/initialize.f90`
- Remove `if (flAgeTracer)` guarded branches in `src/solute/agetracer.f90` and any callers (the body is no-op anyway)
- Verify no remaining bare-global readers via `grep -rnE "\b(flAgeTracer|AgeGwl1m|icAge.*)\b" src/`

### A2: Macropore retired sentinel cleanup

ADR 0040 — macropore is retired.

**Globals retired:** `FlMacropore`, `ArMpSs` (residual), `.BMA` writer globals (`iQInTopLatDm1`, `iQInTopLatDm2`, `IWaSrDm1Beg`, `IWaSrDm2Beg`, `WaSrDm1`, `WaSrDm2`, etc.). Per audit at `docs/superpowers/notes/2026-05-14-gr-final-retirement-inventory.md` — all W-category.

**Actions:**
- Tombstone declarations + adapter writes + zero-fills
- Remove `if (flmacropore)` dead branches if any remain
- Verify no remaining readers

### A3: Output control switches → config / state

**Switches to address:**
- `swcsv`, `swcsv_tz`, `swinc`, `swdrought`, `swrum`, `swtem`, `swvap`, `swsnow_output` (verify list)

**Decision per category:**
- Output-toggle switches → thread `config%output_csv` or `config%simulation` arg (verify per-symbol)
- File unit globals (`logf`, `inc`, `rot`, `crp`, `tem`, `snw`) → **retain via narrow `use variables, only:`** as Arc 9.5 / pragmatic edge. These are integer file handles, not architectural state.

**Actions:**
- For switches that have a config home, migrate readers to `config%X` via threaded arg
- Retire the bare global from variables.f90 (Phase D1/D2)
- File unit globals: document as deferred; no migration this arc

### A4: crop schema additions

Add ~15-20 fields to `state%crop` sub-records to give all cropgrowth-referenced runtime state a home.

#### `state%crop%common` additions

```fortran
! Lifecycle flags
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
integer :: noddrz = 0    !! number of nodes in root zone
real(real64), allocatable :: cumdens(:)   !! cumulative root density (allocated numnod)

! Surface params (verify if config-loaded or runtime)
real(real64) :: albedo = 0.0_real64
real(real64) :: rsc    = 0.0_real64
```

#### `state%crop%wofost` additions

```fortran
! FCO2 derived values (computed from atmospheric CO2 per timestep)
real(real64) :: fco2amax = 1.0_real64
real(real64) :: fco2eff  = 1.0_real64
real(real64) :: fco2tra  = 1.0_real64
```

(Defaults of 1.0 represent "no CO2 effect" — verify legacy default in `initialize.f90` or `FacCO2` subroutine.)

#### `state%crop%grass` additions

Audit-driven during implementation. If audit surfaces grass-specific scratch state not yet on grass subrecord, add.

#### Nutrient overflow

If `anlv`, `anst`, `nni`, `nmaxlv`, `nmaxst`, `nmaxrt` are cropgrowth-referenced and aren't on `state%nutrients`, add to nutrients (not crop).

**ETSine astronomical scratchpad** (rad/daylp/difpp/dsinbe/atmtr): **DEFERRED**. cropgrowth retains narrow `use variables, only: rad, daylp, difpp, dsinbe, atmtr` until a future mini-arc adds these to state%atmosphere. Document explicitly.

### Phase A close

- Cheap-wins retirements landed (~30-50 globals retired)
- crop schema extended with ~15-20 new fields (no readers yet — Phase B routes them)
- check-full 5/5 byte-for-byte

## Phase A.5 — Runtime Dual-Write Coverage

For each NEW field added in A4, audit all legacy-global write sites and add `state%crop%X = <sym>` mirror.

**Audit pattern:**
```bash
for sym in flCropReadFile flCropPrep flCropSow flCropGerm flCropHarvest cropstart cropend PrepDelay SowDelay noddrz cumdens albedo rsc fco2amax fco2eff fco2tra; do
  echo "=== $sym ==="
  grep -rnE "^[[:space:]]*\b$sym\b\s*=" src/ --include="*.f90" | head -10
done
```

Most write sites are in `cropgrowth.f90`, crop init files, `irrigation.f90`. One commit per file batch.

**Critical:** zero-resets count (GR-ATM `8fb76cf` lesson). Inside `nocrop()`, lifecycle transitions, `if (...) then` branches.

### A.5 close

check-full 5/5. State mirrors track legacy through simulation.

## Phase B — cropgrowth.f90 Reader Migration

10 sites, one commit per site.

| # | Line | Subroutine | Notes |
|---|---|---|---|
| B1 | 44 | CropGrowth | Main dispatcher, ~80 symbols. Thread `config` arg. Retain ETSine narrow defer. |
| B2 | 621 | cropfixed | ~18 symbols → state%crop%common/fixed + config%crop%fixed |
| B3 | 894 | cropoutput | File paths — retain narrow; lifecycle flags → state%crop%common |
| B4 | 1005 | nocrop | rd/rdpot/lai/laipot/cf/ch/dvs/tsum/wlv/wst/wrt/wso → state. **Zero-resets critical.** |
| B5 | 1053 | ArableLandGerm | Germination scratch (now in state%crop%common); lat → defer or config |
| B6 | 1243 | FacCO2 | fco2amax/eff/tra → state%crop%wofost; flco2 → config |
| B7 | 1313 | wofost | Biggest site, ~60 symbols. Most fields on state%crop%wofost or config%crop%wofost. |
| B8 | 2598 | grass | ~50 symbols. Most fields on state%crop%grass or config%crop%grass. |
| B9 | 4024 | update_rootdistribution | noddrz/cumdens/wrt/gwrt → state%crop%X |
| B10 | 5198 | sumttd | tsumdepth/tsumtemp/tsumtime — stub path (swtsum=2 errors). Retain or migrate per inspection. |

**Per-site:** inspect → substitute body refs → thread config → update callers → VG → commit.

### Phase B close

- `grep -n "^[[:space:]]*use [Vv]ariables" src/crop/cropgrowth.f90` returns empty OR only documented narrow deferrals (ETSine + file paths)
- check-full 5/5

## Phase C — Adapter Shrinkage

After Phase B, many crop globals become write-only-by-adapter.

| # | Task | Files |
|---|---|---|
| C1 | Drop adapter writes for retired crop globals | `src/io/toml/config_to_variables.f90` |
| C2 | Drop swap_mod A14-A17 dual-write block lines for retired globals | `src/core/swap_mod.f90` |
| C3 | Drop initialize.f90 zero-fills for retired globals | `src/core/initialize.f90` |

### Phase C close

Adapter shrinkage measurable; check-full 5/5.

## Phase D — Variables.f90 Deletion

Iterative tombstoning by category.

| # | Category | Globals |
|---|---|---|
| D1 | AgeTracer + macropore dead-code (from A1+A2) | ~20-30 globals |
| D2 | Output switches retired in A3 | ~5-10 globals |
| D3 | Crop runtime globals (post-B migration) | ~50-80 globals |
| D4 | Residual W-globals not yet covered | catch-all |
| D5 | initialize.f90 zero-fill cleanup | matching |

Each: tombstone + iterative compile-error patch + VG + commit.

### Phase D close

variables.f90 significantly reduced (target: 1418L → <800L). check-full 5/5.

## Phase E — Close + Roadmap

### E1: Final verification

- check-full 5/5
- BMI + cffi-demo passing
- pFUnit all-pass
- Sanity greps

### E2: Update roadmap

Edit `docs/superpowers/plans/2026-05-13-globals-retirement-roadmap.md` — document remaining work for full retirement:

```markdown
**Remaining after GR-CROPRT (as of 2026-05-14):**
- ETSine astronomical scratchpad (rad/daylp/difpp/atmtr/dsinbe/lat) — schema decision per-symbol; dedicated mini-arc
- File path globals (outfil/pathwork/project/cropfil) + log unit (logf) + file unit handles (inc/rot/crp/tem/snw) — state%io subrecord OR config%paths decision
- swapoutput.f90 remaining 13 narrow sites — extension of output schema
- WOFOST/grass lookup tables not in state — schema additions per crop type
- Crop config parameter threading residuals — narrow imports for read-only config still in `use variables, only:`
- variables.f90 + config_to_variables.f90 + initialize.f90 final deletion — when ALL consumers fully migrated
```

### E3: Arc-complete marker

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-croprt): GR-CROPRT complete — crop runtime finalization + cheap wins

check-full 5/5 byte-for-byte. BMI + cffi-demo passing. pFUnit all-pass.

Net retirement:
- AgeTracer dead-code globals (ADR 0032): N1 globals retired
- Macropore retired sentinels (ADR 0040): N2 globals retired
- Output control switches: N3 migrated to state/config
- Crop schema additions: ~15-20 fields in state%crop sub-records
- Phase A.5 runtime dual-write coverage for new fields
- cropgrowth.f90: 10 sites migrated; ~738+ symbol substitutions
- variables.f90: 1418L → <X>L (target <800L)
- config_to_variables.f90: 1742L → <X>L
- initialize.f90: 835L → <X>L

Inherited to future continuation arc (documented in roadmap):
- ETSine astronomical scratchpad
- File path + log unit globals
- swapoutput.f90 remaining narrow sites
- WOFOST/grass lookup tables (if not addressed in this arc)
- Full adapter deletion (variables.f90 + config_to_variables.f90 + initialize.f90)

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

## Risk & Mitigation

| Risk | Mitigation |
|---|---|
| cropgrowth.f90 substitution density (~738+ refs) | One site per commit (10 commits); VG after each |
| Phase A.5 misses zero-resets | Audit grep explicitly catches `<sym> = 0` / `.false.` patterns |
| Schema additions surface during Phase B (not in A4) | Controller extends schema inline; no deferral |
| FCO2 fields semantically tricky (multi-write paths) | Audit FacCO2 carefully; dual-write at every write site |
| AgeTracer guard removal might break ADR 0032 invariant | Verify `flAgeTracer == .false.` is unconditional; preserve the guard if uncertain |
| Phase D compile errors cascade | Iterative patch cycle; expect 5-15 patches per task |

## What This Arc DOES NOT Change

- BMI / cffi surfaces — continue working throughout
- Physics — every value bit-identical
- swapoutput.f90 remaining narrow sites (deferred to a continuation arc — not blocking cropgrowth)
- ETSine astronomical scratchpad (deferred)
- File path globals (pragmatic retain)

## Effort

- **Total commits:** ~30-40 estimated
- **Estimated effort:** 6-10 days subagent-driven
- **Phases gated by check-full 5/5 byte-for-byte**

## Memory & ADR Consequences

- After arc closes, update memory `project_state_rescue_complete_2026-05-12.md`.
- ADR candidate: if substantial schema additions surface, consider ADR 0045 (`state%crop` field set finalized).

## Files & Artifacts

- **This spec:** `docs/superpowers/specs/2026-05-14-globals-croprt-design.md`
- **Plan (next):** `docs/superpowers/plans/2026-05-14-globals-croprt.md`
- **Roadmap context:** `docs/superpowers/plans/2026-05-13-globals-retirement-roadmap.md`
- **Inventory:** `docs/superpowers/notes/2026-05-14-gr-final-retirement-inventory.md`
- **Prior arcs:** GR-UTILS (`26123ad`), GR-BH (`a65bdf3`), GR-ATM (`57b00b3`), GR-CROP (`476efba`), GR-FINAL (`21e15dc`)
