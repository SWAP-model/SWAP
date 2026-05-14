# GR-FINAL Implementation Plan — Final retirement (delete the 3 adapter files)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Migrate all 32 remaining `use variables` consumers, retire the swap_mod strangler-fig leftovers (transient buffers + swinco=3 inline block), then delete `src/core/variables.f90`, `src/io/toml/config_to_variables.f90`, `src/core/initialize.f90` — ~4040 lines of legacy plumbing gone. End state: all data flows through typed `state%X` + `config%X`, no bare globals.

**Architecture:** Five phases. Phase A audits retirement readiness. Phase B migrates major clusters (readmeteo, swapoutput, cropgrowth + crop runtime, timecontrol, soil/drainage). Phase C retires swap_mod strangler-fig + incrementally shrinks config_to_variables.f90. Phase D iteratively deletes globals from variables.f90 by category. Phase E deletes the 3 adapter files. Compile errors at deletion time are expected and patched in same task.

**Tech Stack:** Fortran 2018, gfortran, Meson, pFUnit 4.15, pixi, Python regression tests.

**Spec:** `docs/superpowers/specs/2026-05-14-globals-final-retirement-design.md`

**Builds on:** GR-CROP (`476efba`), GR-ATM (`57b00b3`), GR-BH (`a65bdf3`), GR-UTILS (`26123ad`).

**Verification gate per task (VG):**

```bash
rm -rf builddir && pixi run build-linux \
  && pixi run test-pfunit \
  && pixi run -e test python tests/regression/test_output_regression.py \
       hupselbrook surfacewater salinitystress grassgrowth
```

Expected: clean build, pFUnit 741, regression **4/4 byte-for-byte**. Any deviation = task NOT complete; fix or escalate.

**Phase close gate:** `pixi run check-full` 5/5 byte-for-byte.

**Controller quality bar:** subagents do NOT defer due to missing infrastructure. If a deferral surfaces, controller intervenes inline — thread state/config through caller chain, extend schema, migrate sibling reader in same commit.

---

## Phase A — Preparation (Tasks 1–4)

### Task 1: Pre-flight baseline

**Files:** None.

- [ ] **Step 1:** `rm -rf builddir && pixi run build-linux` — must succeed.
- [ ] **Step 2:** `pixi run test-pfunit` — record count.
- [ ] **Step 3:** `pixi run check-full` — 5/5 pass.
- [ ] **Step 4:** Empty marker commit:

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-final): pre-flight baseline — Arc 9 final retirement begins

check-full 5/5. Baseline locked before 32-file consumer migration,
swap_mod strangler-fig elimination, and deletion of variables.f90 +
config_to_variables.f90 + initialize.f90 (~4040 lines).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Audit retirement readiness

**Files:**
- Create: `docs/superpowers/notes/2026-05-14-gr-final-retirement-inventory.md`

- [ ] **Step 1:** Inventory declarations in `variables.f90`:

```bash
grep -nE "^[[:space:]]*\b(real\(8\)|integer|logical|character)\b" src/core/variables.f90 | grep -v "^[[:space:]]*!" | head -100
```

- [ ] **Step 2:** For each declaration, classify into one of 4 categories:
  - **W**: Write-only-by-adapter (no readers) → directly deletable
  - **R**: Read-only-by-readers (consumers still read; need migration)
  - **B**: Both (need migration + adapter-rewrite)
  - **C**: Config-loaded mirror (convert to direct config-sourcing then delete legacy)

Audit each symbol via:
```bash
sym=arad   # example
hits=$(grep -rnE "\b${sym}\b" src/ --include="*.f90" \
  | grep -v "state%\|config%\|! \|variables.f90\|use variables, only:\|integer ::\|real(8) ::\|real(real64) ::\|character\|^src/core/initialize\|^src/io/toml/config_to_variables\|^src/core/swap_mod" | wc -l)
echo "$sym: $hits non-write hits → category $cat"
```

- [ ] **Step 3:** Write inventory to `docs/superpowers/notes/2026-05-14-gr-final-retirement-inventory.md` with format:

```markdown
## Retirement inventory — variables.f90

| Symbol | Category | Notes |
|---|---|---|
| arad | R | Read by swapoutput, readmeteo writes → state%atmosphere — migrate in B1+B3 |
| ... | ... | ... |

## Summary
- W (deletable now): N1 symbols
- R (need reader migration): N2 symbols
- B (need both): N3 symbols
- C (config-loaded mirror): N4 symbols
```

- [ ] **Step 4:** Commit inventory:

```bash
git add docs/superpowers/notes/2026-05-14-gr-final-retirement-inventory.md
git commit -m "docs(gr-final A2): retirement readiness inventory

Per-symbol classification of variables.f90 declarations. Drives Phase B
migration ordering + Phase D retirement category sequencing.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 3: Add `max_resp_factor` to config (if missing)

**Files:** `src/config/crop_config.f90` or `src/config/cropfixed_config.f90`; `src/io/toml/config_to_variables.f90`

- [ ] **Step 1:** Verify `max_resp_factor` placement:

```bash
grep -rn "max_resp_factor" src/config/*.f90 src/io/toml/config_to_variables.f90 src/core/variables.f90 | head
```

- [ ] **Step 2:** If absent from config, add to most consistent location. Look at the legacy populator path to determine where (likely `cropfixed_config_t` or a sub-record).

```fortran
real(real64) :: max_resp_factor = 1.0_real64  !! oxygen stress max respiration factor (default per legacy)
```

(Verify legacy default in `initialize.f90` or `oxygenstress.f90`.)

- [ ] **Step 3:** Add adapter populator line in `config_to_variables.f90`:

```fortran
max_resp_factor = config%crop%fixed%max_resp_factor
```

- [ ] **Step 4:** Run VG.

- [ ] **Step 5:** Commit.

```bash
git add src/config/*.f90 src/io/toml/config_to_variables.f90
git commit -m "schema(gr-final A3): add max_resp_factor to config%crop

Closes one outstanding crop config gap. Required for Phase B6 oxygenstress
migration without retaining the narrow use-variables import.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

If already present in config, skip with explicit no-op commit:
```bash
git commit --allow-empty -m "chore(gr-final A3): max_resp_factor already in config — no-op"
```

---

### Task 4: Phase A close marker

**Files:** None.

- [ ] **Step 1:** check-full 5/5.

- [ ] **Step 2:**
```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-final): Phase A complete — retirement readiness audited

Inventory note + max_resp_factor schema gap closed. check-full 5/5.
Phase B (cluster migrations) unblocked.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Phase B — Major cluster migrations (Tasks 5–16)

### Task 5: Migrate readmeteo blanket (line 17)

**Files:** `src/io/readmeteo.f90`

- [ ] **Step 1:** Inspect the bare site:

```bash
sed -n '12,30p' src/io/readmeteo.f90
```

- [ ] **Step 2:** Audit ALL symbols read in `readmeteo`'s subroutines that ALSO appear in `variables.f90`. List candidates:

```bash
grep -nE "^[[:space:]]*[A-Za-z_].*=\|read\(" src/io/readmeteo.f90 | head -50
```

- [ ] **Step 3:** Replace `use variables` at line 17 with:

```fortran
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   ! Add narrow imports for symbols NOT yet on state/config (with [SS-GR-FINAL B1] DEFERRED comments)
```

- [ ] **Step 4:** For each affected subroutine, ensure `state` (and optionally `config`) is in scope. Add to signatures where needed. Update callers.

- [ ] **Step 5:** Substitute body references per the symbol replacement table from GR-CROP / GR-ATM specs. Atmosphere arrays → `state%atmosphere%X`. Time → `state%timecontrol%X`. etc.

- [ ] **Step 6:** Run VG.

- [ ] **Step 7:** Commit.

```bash
git add src/io/readmeteo.f90 # + any caller files
git commit -m "$(cat <<'EOF'
refactor(gr-final B1): readmeteo bare use variables — drop

Major unlock: closes the blocker for atmosphere multi-consumer global
retirement (D1). Bodies route to state%atmosphere / state%timecontrol /
state%mesh. Document any retained narrow deferrals.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: Migrate readmeteo CSV writers (lines 267, 366, 461, 554)

**Files:** `src/io/readmeteo.f90`

- [ ] **Step 1:** Inspect each narrow site:

```bash
for line in 267 366 461 554; do
  echo "=== Line $line ==="
  start=$((line - 5))
  end=$((line + 25))
  sed -n "${start},${end}p" src/io/readmeteo.f90
done
```

- [ ] **Step 2:** For each writer subroutine, migrate writes from legacy globals (`arad/atmn/atmx/ahum/awin/arai/aetr/wet/atav/ad/am/...`) to `state%atmosphere%X` directly.

Per-site:
- Line 267: `raincsv_dat, nraincsv` — rain CSV input; route to state%atmosphere or config if static
- Lines 366, 461: daily meteo arrays — route to state%atmosphere
- Line 554: `metcsv_det, nmetcsv_det` — meteo detail; route to state%atmosphere

- [ ] **Step 3:** After migration, drop the narrow `use variables, only:` line at each site (or retain with `[SS-GR-FINAL B2] DEFERRED` if symbol not yet on state).

- [ ] **Step 4:** Run VG.

- [ ] **Step 5:** Commit.

```bash
git add src/io/readmeteo.f90
git commit -m "$(cat <<'EOF'
refactor(gr-final B2): readmeteo CSV writers (4 sites) — drop use variables

Per-year reload routes meteo + rain CSV directly into state%atmosphere.
Closes atmosphere multi-consumer migration. legacy arad/atmn/atmx/... are
now write-only-by-adapter (retirement candidate for Phase D1).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 7: Migrate swapoutput.f90 — first half (sites 1-7)

**Files:** `src/io/swapoutput.f90`

- [ ] **Step 1:** Inventory all 15 `use variables` sites:

```bash
grep -nE "^[[:space:]]*use [Vv]ariables" src/io/swapoutput.f90
```

Identify which subroutine each site belongs to (search backward for `subroutine`).

- [ ] **Step 2:** Migrate sites 1-7 (typically the first half by line number). For each:
  - Replace `use variables, only: ...` with state/config imports
  - Substitute body references per spec table
  - Update callers if state/config arg added

Symbols in swapoutput typically include: numnod/dz/z/disnod/ztopcp/zbotcp/layer (→ state%mesh), cgrai/cnrai/caintc/cevap/etc. (→ state%atmosphere), crop fields (→ state%crop%X), output flag globals (→ state%atmosphere or config), water balance cumulatives.

- [ ] **Step 3:** Run VG.

- [ ] **Step 4:** Commit.

```bash
git add src/io/swapoutput.f90
git commit -m "$(cat <<'EOF'
refactor(gr-final B3): swapoutput.f90 sites 1-7 — drop use variables

Output writers route reads through state references. <document sites
migrated and any narrow deferrals retained>.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 8: Migrate swapoutput.f90 — second half (sites 8-15)

**Files:** `src/io/swapoutput.f90`

- [ ] **Step 1:** Migrate remaining sites.
- [ ] **Step 2:** After this task, `grep -n "^[[:space:]]*use [Vv]ariables" src/io/swapoutput.f90` should be empty (or only documented narrow deferrals).
- [ ] **Step 3:** Run VG.
- [ ] **Step 4:** Commit.

```bash
git commit -m "$(cat <<'EOF'
refactor(gr-final B4): swapoutput.f90 sites 8-15 — drop use variables

All 15 sites migrated. swapoutput.f90 clean of bare use variables
modulo documented deferrals.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 9: Migrate cropgrowth.f90 narrow imports (10 sites)

**Files:** `src/crop/cropgrowth.f90`

- [ ] **Step 1:** Inventory:

```bash
grep -nE "^[[:space:]]*use [Vv]ariables" src/crop/cropgrowth.f90
```

- [ ] **Step 2:** For each site, classify each `only:` symbol:
- Crop config (idev/tsumea/tsumam/tbase/cf/ch/etc.) → `config%crop%fixed/grass/wofost%X` (via threaded config arg)
- WOFOST nutrients (lrnr/lsnr/nlue/rnflv/rnfst/frnx/nmxlv) → `config%crop%wofost%nutrient%X`
- Oxygen stress params (c_mroot/f_senes/q10_root/q10_microbial/etc.) → `config%crop%fixed%X` per Phase A3 finding
- Array dim parameters (macp/magrs/maho) → `swap_array_dimensions`
- Crop runtime state already migrated → `state%crop%X`
- Any non-migrated → narrow defer with `[SS-GR-FINAL B5] DEFERRED` comment

- [ ] **Step 3:** Thread `config` through any subroutine that doesn't already have it. Controller intervenes inline if subagent reports difficulty.

- [ ] **Step 4:** Migrate body references per the symbol classification.

- [ ] **Step 5:** Run VG.

- [ ] **Step 6:** Commit.

```bash
git add src/crop/cropgrowth.f90 # + any callers
git commit -m "$(cat <<'EOF'
refactor(gr-final B5): cropgrowth.f90 (10 sites) — drop use variables

Config-loaded crop parameters route through config%crop%X. Array dim
parameters from swap_array_dimensions. Runtime crop state via state%crop.
Controller threaded config through <list subroutines>.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 10: Migrate other crop runtime files

**Files:**
- `src/crop/oxygenstress.f90` (4 sites)
- `src/crop/rootextraction.f90` (4 sites)
- `src/crop/irrigation.f90` (1 site at line 320)
- `src/crop/tillage.f90` (1 site)
- `src/crop/management_soil.f90` (1 narrow site)

- [ ] **Step 1:** Per-file inventory + migration. Same pattern as Task 9.

- [ ] **Step 2:** Special — `irrigation.f90` line 320: `mairg, irrigevent, qssdi, qssdisum, dt_SSDI_event` etc. Determine per-symbol whether on state%X, config%irrigation%X, or array dim. SSDI-event tracking may need state%X if not already.

- [ ] **Step 3:** Run VG after each file (or batch into one commit if commits stay small).

- [ ] **Step 4:** Commit.

```bash
git add src/crop/oxygenstress.f90 src/crop/rootextraction.f90 src/crop/irrigation.f90 src/crop/tillage.f90 src/crop/management_soil.f90
git commit -m "refactor(gr-final B6): crop runtime files — drop use variables

oxygenstress/rootextraction/irrigation/tillage/management_soil narrow
imports migrated to state%X + config%X per audit.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 11: Migrate crop init files

**Files:**
- `src/crop/cropfixed_init.f90` (1 site)
- `src/crop/cropgrass_init.f90` (2 sites)
- `src/crop/cropwofost_init.f90` (2 sites)
- `src/crop/wofost_soil_parameters.f90` (1 site)

- [ ] **Step 1-2:** Per-file migration. Mostly config threading for crop parameters.

- [ ] **Step 3:** Run VG.

- [ ] **Step 4:** Commit.

```bash
git commit -m "refactor(gr-final B7): crop init files — drop use variables

All crop init files clean. Crop config parameters thread via config arg.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 12: Migrate timecontrol_mod.f90 (5 sites)

**Files:** `src/core/timecontrol_mod.f90`

- [ ] **Step 1:** Inventory sites.
- [ ] **Step 2:** Migrate per-symbol. Most fields should already be on state%timecontrol.
- [ ] **Step 3:** Run VG.
- [ ] **Step 4:** Commit.

```bash
git commit -m "refactor(gr-final B8): timecontrol_mod.f90 — drop use variables

Residual narrow imports migrated to state%timecontrol / config%X.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 13: Migrate soil cluster

**Files:** `src/soil/soilhydraulics.f90` (4 sites), `src/soil/waterbalance.f90` (3 sites), `src/soil/soilgrid.f90` (residuals)

- [ ] **Step 1-2:** Per-file migration.

- [ ] **Step 3:** Close any boundtop-related deferrals if soilhydraulics calls boundtop (config threading via headcalc chain). Per user directive — inline fix.

- [ ] **Step 4:** Run VG.

- [ ] **Step 5:** Commit.

```bash
git commit -m "refactor(gr-final B9): soil cluster — drop use variables

soilhydraulics/waterbalance/soilgrid narrow imports migrated. Threading
config through headcalc → boundtop closes GR-BH residual deferral
(swkmean/swredu/flrunon/runonarr).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 14: Migrate drainage cluster

**Files:** `src/drainage/drainage.f90` (4 sites), `src/drainage/surfacewater.f90`, `divdra.f90` (residuals)

- [ ] **Step 1-2:** Per-file migration.
- [ ] **Step 3:** Run VG.
- [ ] **Step 4:** Commit.

```bash
git commit -m "refactor(gr-final B10): drainage cluster — drop use variables

drainage/surfacewater/divdra narrow imports migrated.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 15: Migrate small consumers (boundary/heat/solute/atmosphere residuals)

**Files:**
- `src/boundary/boundbottom.f90`, `boundtop.f90` (narrow deferrals)
- `src/heat/temperature.f90`, `frozencond.f90` (residuals)
- `src/solute/solute.f90`, `agetracer.f90` (residuals)
- `src/atmosphere/meteoday.f90`, `meteodt.f90`, `et.f90`, `interception.f90` (residual narrow imports)

- [ ] **Step 1:** Per-file inventory. For each remaining `use variables, only: ...` line:
  - If symbol is now state/config/array-dim accessible: migrate
  - If still genuinely deferred: document with comment + reason

- [ ] **Step 2:** Some narrow imports may stay if the symbol is a true edge case (e.g., `logf` global file unit number). Document these as Arc-9-final-deferrals if they survive.

- [ ] **Step 3:** Run VG.

- [ ] **Step 4:** Commit.

```bash
git commit -m "refactor(gr-final B11): small consumers — drop use variables

boundary/heat/solute/atmosphere residual narrow imports migrated.
Remaining narrow deferrals documented as final-arc edge cases.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 16: Phase B close marker

**Files:** None.

- [ ] **Step 1:** check-full 5/5.

- [ ] **Step 2:** Sanity grep:
```bash
grep -rn "^[[:space:]]*use [Vv]ariables" src/ --include="*.f90" | grep -v "config_to_variables\|src/core/initialize\|src/core/variables" | head -20
```

Expected: empty (or only minimal documented final deferrals). The remaining `use variables` lines should be in `config_to_variables.f90` and `initialize.f90` which will be deleted in Phase E.

- [ ] **Step 3:** Phase B close marker:

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-final): Phase B complete — all readers migrated

32-file consumer migration done. Only config_to_variables.f90 and
initialize.f90 still use variables (self-references; deleted in Phase E).
check-full 5/5 byte-for-byte.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Phase C — Strangler-fig elimination (Tasks 17–21)

### Task 17: Retire swap_mod transient buffer reads

**Files:** `src/core/swap_mod.f90`, `src/io/toml/config_to_variables.f90`

- [ ] **Step 1:** Inspect transient buffer reads:

```bash
grep -nE "_init_buf\b" src/core/swap_mod.f90 src/io/toml/config_to_variables.f90 | head -20
```

Identify each `tc_iyear_init_buf`, `tc_imonth_init_buf`, `tc_dt_init_buf`, `h_init_buf`, `pondini_init_buf`, `pond_init_buf` read site in swap_mod.

- [ ] **Step 2:** For each buffer, determine its source in config:

```bash
grep -nE "_init_buf\s*=" src/io/toml/config_to_variables.f90 | head
```

Each `<buf> = config%X%Y` line shows the canonical config path.

- [ ] **Step 3:** In swap_mod, replace each `state%X = <buf>` line with `state%X = config%path%to%source` directly.

Example:
```fortran
! Before:
state%timecontrol%iyear = tc_iyear_init_buf
! After:
state%timecontrol%iyear = config%simulation%iyear   ! per A2 audit + config path
```

For `h_init_buf` (array), the dual-write block at lines 297-310:
```fortran
! Before:
if (allocated(h_init_buf)) then
   do ki = 1, min(size(h_init_buf), size(state%soilwater%h))
      state%soilwater%h(ki) = h_init_buf(ki)
   end do
   deallocate(h_init_buf)
end if
! After:
if (allocated(config%soil%initial%h_init)) then
   do ki = 1, min(size(config%soil%initial%h_init), size(state%soilwater%h))
      state%soilwater%h(ki) = config%soil%initial%h_init(ki)
   end do
end if
```

Drop the `deallocate` since config arrays are owned by the config_t.

- [ ] **Step 4:** Drop the buffer declarations from `config_to_variables_mod` (or wherever they're declared).

- [ ] **Step 5:** Remove `use config_to_variables_mod, only: h_init_buf, ...` from swap_mod if it becomes empty.

- [ ] **Step 6:** Run VG.

- [ ] **Step 7:** Commit.

```bash
git add src/core/swap_mod.f90 src/io/toml/config_to_variables.f90
git commit -m "$(cat <<'EOF'
refactor(gr-final C1): retire swap_mod transient buffer reads

h_init_buf/pondini_init_buf/pond_init_buf/tc_iyear_init_buf/
tc_imonth_init_buf/tc_dt_init_buf reads in swap_mod replaced with
direct config%X reads. Buffer declarations dropped from
config_to_variables_mod. Strangler-fig leftover from SS-DRV/SS-TCM
era retired.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 18: Retire swap_mod swinco=3 inline block

**Files:** `src/core/swap_mod.f90`

- [ ] **Step 1:** Inspect block:

```bash
grep -n -A 30 "config%soil%swinco == 3" src/core/swap_mod.f90 | head -40
```

- [ ] **Step 2:** Determine if the block is redundant:
- If atmosphere state mirrors from prior arcs cover the warm-restart seeding → block is redundant → fully retire
- If the block does something unique (e.g., reads a specific config field that no other path reads) → simplify to use that config field directly

Likely outcome: the block reduces to a few `state%atmosphere%X = config%X%Y` lines or retires entirely.

- [ ] **Step 3:** Drop or simplify the block. Run VG (regression must hold for `swinco=3` test cases — verify).

- [ ] **Step 4:** Commit.

```bash
git commit -m "$(cat <<'EOF'
refactor(gr-final C2): retire swap_mod swinco=3 inline block

<document outcome: fully retired | simplified to N lines>. Warm-restart
seeding now uses state mirrors from prior arcs.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 19: Drop adapter dual-write blocks for write-only globals (incremental)

**Files:** `src/core/swap_mod.f90`

- [ ] **Step 1:** Identify `[SS-GR-ATM A...]` and `[SS-GR-CROP A...]` blocks:

```bash
grep -nE "SS-GR-(ATM|CROP) A" src/core/swap_mod.f90 | head -40
```

Each `state%X = legacy_sym` line is a candidate for retirement IF:
- `legacy_sym` is in the "write-only" category from A2 inventory
- `legacy_sym` will be deleted in Phase D
- OR `state%X` is now populated directly from `config%X` elsewhere

- [ ] **Step 2:** Drop each redundant line. Pace this — small commits per block.

- [ ] **Step 3:** Run VG after each commit.

- [ ] **Step 4:** Final commit (or multiple):

```bash
git commit -m "$(cat <<'EOF'
refactor(gr-final C3): drop adapter dual-write blocks (incremental)

state%X = legacy_sym lines retired where legacy_sym is write-only or
state is now sourced from config directly. <list blocks retired>.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

This task may produce multiple commits.

---

### Task 20: Shrink config_to_variables.f90 — drop legacy writes for retired globals

**Files:** `src/io/toml/config_to_variables.f90`

- [ ] **Step 1:** For each global being retired in Phase D, find the corresponding adapter write line.

```bash
grep -nE "^[[:space:]]*\b<sym>\b\s*=" src/io/toml/config_to_variables.f90
```

- [ ] **Step 2:** Drop the legacy write. The state mirror should now source from config%X directly (verified in Phase B/C migrations).

- [ ] **Step 3:** Run VG after each batch.

- [ ] **Step 4:** Commit incrementally:

```bash
git commit -m "refactor(gr-final C4): shrink config_to_variables — drop retired writes

<document symbols dropped>.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

May produce multiple commits.

---

### Task 21: Phase C close marker

**Files:** None.

- [ ] check-full 5/5.

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-final): Phase C complete — strangler-fig retired

swap_mod transient buffer reads + swinco=3 inline retired. Dual-write
blocks for write-only globals dropped. config_to_variables.f90 shrunk
significantly. check-full 5/5 byte-for-byte.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Phase D — Global retirement (Tasks 22–28)

Iterative; each task targets one category. Compile errors in any task → patch in same task.

### Task 22: Retire atmosphere multi-consumer globals

**Files:** `src/core/variables.f90`, `src/core/initialize.f90`, plus compile-surfaced patches.

- [ ] **Step 1:** Identify candidates (per A2 inventory + Phase B status):
- arad/atmn/atmx/ahum/awin/arai/aetr/wet (daily arrays)
- atav/epot/tpot/grain/nrain (sub-daily arrays)
- ad/am (date arrays)
- daynrfirst/daynrlast/atmin7/nofd/teprrain/teprsnow
- Tav/tav/tavd/rh

- [ ] **Step 2:** For each, verify no readers remain:
```bash
sym=arad
grep -rnE "\b${sym}\b" src/ --include="*.f90" | grep -v "state%\|config%\|! \|variables.f90\|use variables, only:\|initialize.f90\|config_to_variables\|swap_mod.f90\|integer ::\|real(8) ::\|real(real64) ::"
```

- [ ] **Step 3:** Tombstone declarations in `variables.f90`. Drop zero-fills in `initialize.f90`.

- [ ] **Step 4:** Iterative build — `rm -rf builddir && pixi run build-linux 2>&1 | grep Error | head` — patch each compile error in same task.

- [ ] **Step 5:** Run VG (regression must hold).

- [ ] **Step 6:** Commit:

```bash
git commit -m "$(cat <<'EOF'
retire(gr-final D1): delete atmosphere multi-consumer globals

<list deleted>. state%atmosphere is sole home. Compile-surfaced readers
patched in <list files>.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 23: Retire crop runtime globals

**Files:** Same pattern as Task 22 for crop fields (after Phase B5+B6+B7 migration).

Candidates: any crop fields still in variables.f90 after the GR-CROP retirement — these are the ones whose readers needed Phase B threading.

- [ ] Tombstone + zero-fill drop + iterative patch + VG + commit.

---

### Task 24: Retire nutrient globals (if any)

**Files:** Same pattern. Most nutrients are in `Wofost_Soil_Declarations` module, not `variables.f90`. Verify what (if anything) remains in variables.f90.

- [ ] Tombstone + VG + commit.

---

### Task 25: Retire timecontrol globals

**Files:** Same pattern for any remaining bare timecontrol fields.

- [ ] Tombstone + VG + commit.

---

### Task 26: Retire soil/drainage runtime globals

**Files:** Same pattern for soil/drainage residuals.

- [ ] Tombstone + VG + commit.

---

### Task 27: Retire residual scattered globals

**Files:** Same pattern for catch-all.

- [ ] Tombstone + VG + commit.

---

### Task 28: Phase D close marker — variables.f90 sanity check

**Files:** None.

- [ ] **Step 1:** Sanity grep:

```bash
grep -nE "^[[:space:]]*\b(real\(8\)|integer|logical|character)\b" src/core/variables.f90 | grep -v "^[[:space:]]*!" | head
```

Expected: empty (or only deprecated-comments).

- [ ] **Step 2:** check-full 5/5.

- [ ] **Step 3:** Phase D close marker:

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-final): Phase D complete — variables.f90 emptied

All ~150 bare global declarations retired. variables.f90 contains only
tombstone comments + module wrapper. check-full 5/5 byte-for-byte.
Ready for Phase E (file deletion).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Phase E — Final adapter deletion (Tasks 29–35)

### Task 29: Verify variables.f90 has no remaining declarations

**Files:** None.

- [ ] **Step 1:** Sanity grep (same as Phase D close).
- [ ] **Step 2:** If anything remains, return to Phase D. Otherwise advance.

---

### Task 30: Delete src/core/variables.f90

**Files:** `src/core/variables.f90` (delete), `meson.build`, `tests/unit/meson.build`

- [ ] **Step 1:** Delete the file:

```bash
git rm src/core/variables.f90
```

- [ ] **Step 2:** Remove `'src/core/variables.f90',` from `meson.build` legacy sources list.

- [ ] **Step 3:** Remove from `tests/unit/meson.build` if listed.

- [ ] **Step 4:** Iterative build:

```bash
rm -rf builddir && pixi run build-linux 2>&1 | grep Error | head -20
```

Expected: errors for any remaining `use variables` lines in `config_to_variables.f90` or `initialize.f90` (which are about to be deleted). Either:
- Drop the corresponding `use variables` line in those files (they'll be deleted in E3/E4 anyway)
- OR proceed with E3/E4 first, then come back to E2

If compile errors point to files OUTSIDE the 3 adapter files, those are missed readers — patch them.

- [ ] **Step 5:** Run VG.

- [ ] **Step 6:** Commit:

```bash
git add src/core/variables.f90 meson.build tests/unit/meson.build
git commit -m "$(cat <<'EOF'
retire(gr-final E2): delete src/core/variables.f90

~1414 lines of legacy global declarations deleted. End of the strangler-
fig demolition. <list compile-surfaced patches>.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 31: Delete src/io/toml/config_to_variables.f90

**Files:** `src/io/toml/config_to_variables.f90` (delete), `meson.build`, `src/core/swap_mod.f90`

- [ ] **Step 1:** Verify the file is now redundant:
```bash
wc -l src/io/toml/config_to_variables.f90
```

After Phase C+D, file should be mostly empty or only contain module wrapper + dead-code populator lines.

- [ ] **Step 2:** Delete:

```bash
git rm src/io/toml/config_to_variables.f90
```

- [ ] **Step 3:** Remove from `meson.build`.

- [ ] **Step 4:** Remove `use config_to_variables_mod` line(s) in `swap_mod.f90`.

- [ ] **Step 5:** Iterative build — patch any compile errors.

- [ ] **Step 6:** Run VG.

- [ ] **Step 7:** Commit.

```bash
git commit -m "$(cat <<'EOF'
retire(gr-final E3): delete src/io/toml/config_to_variables.f90

~1775 lines of legacy adapter deleted. State mirrors source directly
from config%X. swap_mod use-line dropped.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 32: Delete src/core/initialize.f90

**Files:** `src/core/initialize.f90` (delete), `meson.build`, `src/core/swap_mod.f90`

- [ ] **Step 1:** Verify the file is now redundant. After Phase D, only legacy zero-fills remain (and those were dropped per task).

- [ ] **Step 2:** Delete:

```bash
git rm src/core/initialize.f90
```

- [ ] **Step 3:** Remove from `meson.build`.

- [ ] **Step 4:** Remove `call initialize(...)` from `swap_mod.f90`:

```bash
grep -n "call initialize\b\|use initialize_mod" src/core/swap_mod.f90
```

- [ ] **Step 5:** Iterative build — patch any compile errors.

- [ ] **Step 6:** Run VG.

- [ ] **Step 7:** Commit.

```bash
git commit -m "$(cat <<'EOF'
retire(gr-final E4): delete src/core/initialize.f90

~851 lines of legacy zero-fill code deleted. State subrecords own their
own initialization (via their init type-bound procedures).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 33: meson.build + tests/unit/meson.build cleanup

**Files:** `meson.build`, `tests/unit/meson.build`

- [ ] **Step 1:** Inspect both files for any dangling references to deleted modules.

```bash
grep -nE "variables|config_to_variables|initialize" meson.build tests/unit/meson.build
```

- [ ] **Step 2:** Remove dangling entries.

- [ ] **Step 3:** Verify build still clean. VG.

- [ ] **Step 4:** Commit.

```bash
git commit -m "$(cat <<'EOF'
chore(gr-final E5): meson cleanup — drop deleted adapter file references

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 34: Final ADR — post-retirement state architecture

**Files:** `docs/superpowers/decisions/0044-state-architecture-final.md` (NEW)

- [ ] **Step 1:** Create ADR.

```markdown
# ADR 0044 — Post-Retirement State Architecture

## Status
Accepted — 2026-05-14

## Context
After ~600 commits across 5 arcs (GR-UTILS → GR-BH → GR-ATM → GR-CROP → GR-FINAL),
SWAP's legacy bare-global data flow has been retired. The 3 adapter files
(variables.f90, config_to_variables.f90, initialize.f90) — ~4040 lines of
legacy plumbing — are deleted.

## Decision
Two data-flow paradigms only:

1. **`state%X`** (runtime-mutable, typed):
   - state%mesh, state%soilwater, state%atmosphere, state%heat, state%drainage,
     state%surfacewater, state%solute, state%tillage, state%timecontrol,
     state%crop (sub-records: common/fixed/wofost/grass), state%nutrients.
   - Each subrecord owns its data and exposes a type-bound `init` procedure.

2. **`config%X`** (read-only after init, typed):
   - config%general, config%simulation, config%meteo, config%drain, config%soil,
     config%bottom_boundary, config%heat, config%irrigation, config%solute,
     config%surface_water, config%crop, config%output_csv, config%nutrients.

No bare globals. swap_mod orchestrates init via `state%X%init(config%X)` calls
in dependency order, then `swap_run_step(state, config)` per timestep.

## Consequences
- BMI / cffi surfaces depend only on state and config (already true).
- Multi-instance OpenMP parallelization is now structurally possible.
- pyswap rewrite proceeds against a clean architecture.
- variables.f90 / config_to_variables.f90 / initialize.f90 deleted.
- Meson incremental builds no longer suffer the swap_modern/swap_legacy
  cross-static-library .mod dep issue (memory `feedback_state_schema_clean_rebuild`
  becomes obsolete; evaluate retirement).

## Migration history
- 2026-05-12: GR-UTILS — cofgen retirement + utils cleanup (commit 26123ad)
- 2026-05-13: GR-BH — mesh extraction + boundary/heat migration (commit a65bdf3)
- 2026-05-14: GR-ATM — atmosphere globals + crop_state foundation (commit 57b00b3)
- 2026-05-14: GR-CROP — crop sub-records + nutrients_state + rain timing (commit 476efba)
- 2026-05-14: GR-FINAL — adapter deletion (this arc)

## Lessons codified
- Phase A.5 runtime dual-write coverage prevents simulator drift (zero-resets count)
- Controller inline-intervention for missing infrastructure prevents cascading deferrals
- Per-task byte-for-byte regression gate catches bugs at the source
- Audit-driven schema (let the code define the field list, not the spec)
```

- [ ] **Step 2:** Commit.

```bash
git add docs/superpowers/decisions/0044-state-architecture-final.md
git commit -m "$(cat <<'EOF'
docs(adr): ADR 0044 — post-retirement state architecture

Final ADR capping the GR-UTILS → GR-BH → GR-ATM → GR-CROP → GR-FINAL
arc series. Documents the two-paradigm data flow (state + config) and
the migration history.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 35: Final verification + arc-complete marker

**Files:** None.

- [ ] **Step 1:** Clean rebuild + check-full:

```bash
rm -rf builddir && pixi run build-linux && pixi run check-full
```

Expected: **5/5 byte-for-byte**.

- [ ] **Step 2:** BMI + cffi-demo:

```bash
pixi run -e test test-bmi
pixi run -e test test-cffi-demo
```

- [ ] **Step 3:** pFUnit:

```bash
pixi run test-pfunit
```

Expected: all-pass.

- [ ] **Step 4:** Final sanity:

```bash
# All 3 adapter files deleted
ls src/core/variables.f90 src/io/toml/config_to_variables.f90 src/core/initialize.f90 2>&1
# Expected: no such file or directory for all 3

# No use variables anywhere
grep -rn "^[[:space:]]*use [Vv]ariables" src/ --include="*.f90" | head
# Expected: empty

# state%X usage spread
grep -rl "state%mesh%\|state%crop%\|state%nutrients%\|state%atmosphere%" src/ | wc -l
# Expected: ~30-40 files
```

- [ ] **Step 5:** Arc-complete marker:

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-final): GR-FINAL complete — strangler-fig demolished

check-full 5/5 byte-for-byte. BMI + cffi-demo passing. pFUnit all-pass.

Final deletion summary:
- src/core/variables.f90 (~1414 lines) DELETED
- src/io/toml/config_to_variables.f90 (~1775 lines) DELETED
- src/core/initialize.f90 (~851 lines) DELETED
- Total: ~4040 lines of legacy plumbing eliminated

All simulation data flows through typed state%X (runtime) + config%X
(read-only after init). No bare globals. Architecture: clean.

End of the GR-UTILS → GR-BH → GR-ATM → GR-CROP → GR-FINAL arc series.
ADR 0044 caps the migration history.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Plan Self-Review

**Spec coverage:**
- ✅ Preparation + audit → Tasks 1-4
- ✅ readmeteo migration → Tasks 5-6
- ✅ swapoutput migration → Tasks 7-8
- ✅ Crop migration (cropgrowth + runtime + init) → Tasks 9-11
- ✅ timecontrol/soil/drainage → Tasks 12-14
- ✅ Small consumers → Task 15
- ✅ Phase B close → Task 16
- ✅ Strangler-fig retirement → Tasks 17-20
- ✅ Phase C close → Task 21
- ✅ Iterative global retirement → Tasks 22-27
- ✅ Phase D close → Task 28
- ✅ Adapter file deletion → Tasks 29-32
- ✅ Meson cleanup → Task 33
- ✅ Final ADR → Task 34
- ✅ Arc-complete → Task 35

**Placeholder scan:** Several "audit-driven" markers acknowledge that exact symbol lists and routing decisions are determined by the implementer per actual code state. The plan template + spec table provide direction; implementer applies per audit. No "TBD" / "implement later" / unfinished references.

**Type consistency:** State subrecord paths consistent across phases (state%crop%common, state%atmosphere, etc.). Config paths consistent (config%crop%fixed, config%meteo, etc.). swap_mod strangler-fig retirement (C1) precedes adapter deletion (E3) — correct ordering.

**Known soft spots:**
- swapoutput split (B3/B4) — task boundary is "first half / second half" by line number. Implementer determines exact split.
- Crop config threading scope (B5/B6/B7) — substantial work that may need 2-3 commits per file.
- Phase D order: D1 must happen AFTER B1+B2 (atmosphere readmeteo migration); D2 after B5-B7 (crop); etc. Explicit prerequisites stated.
- Compile-error iteration in Phase D + Phase E is expected; implementer patches each in same task.
- Phase E ordering: E2 (variables.f90) can happen before or after E3+E4 depending on residual use-variables in those adapter files. Plan flexes; implementer chooses.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-05-14-globals-final-retirement.md`. Two execution options:

**1. Subagent-Driven (recommended)** — Fresh subagent per task. Controller intervenes inline on infrastructure-gap deferrals per user directive.

**2. Inline Execution** — Batch execution with checkpoints. Higher main-context pressure given this is the capstone arc.

Which approach?
