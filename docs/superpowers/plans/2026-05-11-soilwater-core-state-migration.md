# Soil-Water Core State-Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Carve ~74 residual soil-water-owned globals into `state%soilwater` (extending the 35-field record from boundary + crop-uptake arcs); introduce TWO new cohort sub-records (`soilwater_intermediate_t` 22 fields + `soilwater_cumulative_t` 14 fields); resolve all deferred residuals from boundary/heat arcs (`pond` D5, `gwl` D6, `kmean` D7, `hconduc` `tsoil_node` 0.0 sentinel ADR 0034 residual); preserve byte-identical check-full at every commit. FINAL coupling-surface arc.

**Architecture:** Extension of `soilwater_state_t` with 30 flat instantaneous + 8 per-day + 22 intermediate cohort fields + 14 cumulative cohort fields. Mild ADR 0033 extension via `reset_per_day()` second method on `intr` cohort (`flDayStart` gate, distinct from `flzerointr`). Grid dimensions kept legacy (heat ADR 0034 precedent). Cohort reset consolidation collapses ~55 scattered reset lines into 3 type-bound `reset()` calls. `soilwater_init` signature `(sw, numnod, nlay)` unchanged.

**Tech Stack:** Fortran 2008+ derived types with type-bound procedures, ASSOCIATE, pixi+meson, pFUnit, check-full regression.

**Spec:** `docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md`
**Discovery:** `docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-discovery.md`
**Branch:** all commits go on `refactor/soilwater-core-state` (already created from `development`).

---

## Lessons-learned applied

- **Dual-write before cutover.** Every global write gets a paired `state%soilwater%X = …` mirror before readers are cut over. Drop globals only after all readers are migrated.
- **Verify-state-not-reports.** Run `check-full` output comparison, not just pFUnit, before every commit. pFUnit alone misses global-default regressions (`feedback_verify_before_committing.md`).
- **Pre-flight dual-write coverage check (atmosphere A-2.1 lesson).** Before cutting over any reader file in Phase 2, confirm each field has at least one non-zero dual-write site:
  ```bash
  grep -rn "state%soilwater%<field>\s*=" src/ | grep -v "0\.0\|0_real64"
  ```
  Must return at least one result. If only zero-init or `reset()` writes are found, dual-write coverage is missing — fix in Phase 1 first.
- **ASSOCIATE prefix `sw_`.** Shadow-safe aliasing in heavy compute bodies (convention: `sw_` for soilwater in files that use `use variables` for grid dims like `numnod`). In files where no name collision exists, bare-name aliases are acceptable.
- **Grep field shapes during pre-flight.** Before deleting a global, run both `use variables.*\b<var>\b` AND raw-symbol greps. Both must return zero hits. Also grep ASSOCIATE alias forms.
- **Compile-driven Phase 2.12 is expected — 20–25 hidden readers.** Largest Phase 2.6-style iteration of any arc. Plan for 10–15 iterative fixup commits. Each fixup is its own commit labeled `refactor(state): SS-SWC compile-fix N — <symbol> in <file>`.
- **Cohort sub-record paths are two levels deep.** `state%soilwater%intr%inq(:)`, `state%soilwater%cumu%cqrot`. ASSOCIATE blocks at call sites keep local names short.
- **reset_per_day() is distinct from reset().** `reset()` zeroes all 22 intermediate fields (flzerointr). `reset_per_day()` zeroes only the 8 per-day fields (flDayStart). Never conflate the two call sites.

---

## Cohort field allocation

### Flat instantaneous (~30 fields on `soilwater_state_t`)

**Per-node arrays (17):** `theta(:)`, `thetm1(:)`, `thetar(:)`, `thetas(:)`, `h(:)`, `hm1(:)`, `q(:)`, `k(:)`, `kmean(:)`, `dimoca(:)`, `cofgen(21,:)`, `FrArMtrx(:)`, `fluseksatexm(:)`, `indeks(:)`, `evp(:)`, `thetsl(:)` (per-layer, size numlay)

**Scalars (13):** `pond`, `pondm1`, `pondini`, `gwl`, `gwlm1`, `hatm`, `volact`, `volm1`, `volini`, `wbalance`, `runon`, `fllowgwl`, and the gwl-geometry helpers `nodgwl`, `pegwl`, `bpegwl`, `npegwl`, `gwlflcpzo`, `nodgwlflcpzo`

### `soilwater_intermediate_t` (22 fields, `flzerointr`-reset via `reset()` + 8 per-day fields via `reset_per_day()`)

**Per-node arrays (6):** `inq(:)`, `inqrot(:)`, `inqssdi(:)`, `iqdo(:)`, `iqup(:)`, `IThetaBeg(:)`

**Scalars (8 non-per-day):** `iqrot`, `iqssdi`, `iqredwet`, `iqreddry`, `iqredsol`, `iqredfrs`, `ies0`, `iet0`, `iew0`, `iintc`, `iruno`, `irunoCN`, `irunon`, `iqbot`, `iqtdo`, `iqtup`, `IPondBeg`, `iprec`, `igird`, `inird`

**Per-day scalars (6):** `tra`, `iqredwet_day`, `iqreddry_day`, `iqredsol_day`, `iqredfrs_day`, `iptra_day`

**Per-day arrays (2):** `qpotrot_day(:)`, `qredtot_day(:)`

### `soilwater_cumulative_t` (14 fields, `flzerocumu`-reset)

`cqssdi`, `cqrot`, `cqbot`, `cqbotdo`, `cqbotup`, `cinund`, `crunon`, `crunoff`, `crunoffCN`, `cqtdo`, `cqtup`, `cqprai`, `cgird`, `cnird`

**Total new fields: ~74 (30 F-INST + 8 F-DAY + 22 INTR + 14 CUMU).**

---

## File structure

**Modified (Phase 1 — state type + init + dual-write):**
- `src/state/soilwater_state.f90` — extend `soilwater_state_t`; add cohort types; extend `soilwater_init` (S-1.1, S-1.2)
- `tests/unit/state/test_soilwater_state.pf` — extend pFUnit tests for new cohorts (S-1.1)
- `src/soil/soilhydraulics.f90` — dual-write home (S-1.3, S-1.4a/b, S-1.5); SoilWaterStateVar + hysteresis plumbing (S-1.4)
- `src/soil/waterbalance.f90` — dual-write home (S-1.6); calcgwl intent-flip (S-1.6)
- `src/boundary/boundtop.f90` — pond + kmean(1) dual-write (S-1.7)
- `src/boundary/boundbottom.f90` — kmean(numnod+1) dual-write (S-1.8)
- `src/crop/tillage.f90` — cofgen + theta + h + pond co-write + internal sub plumbing (S-1.9)
- `src/heat/frozencond.f90` — q/theta verify + dual-write if needed (S-1.12)

**Modified (Phase 2 — reader cutover + retire):**
- `src/soil/soilhydraulics.f90` — reads cutover + cohort reset consolidation (S-2.1, S-2.3)
- `src/soil/waterbalance.f90` — reads cutover (S-2.4)
- `src/boundary/boundtop.f90` + `src/boundary/boundbottom.f90` — reads cutover (S-2.5)
- `src/crop/tillage.f90` — reads cutover (S-2.6)
- `src/crop/rootextraction.f90`, `src/crop/cropgrowth.f90`, `src/crop/oxygenstress.f90` — reads cutover (S-2.7)
- `src/solute/solute.f90`, `src/drainage/drainage.f90`, `src/drainage/surfacewater.f90` — reads cutover (S-2.8)
- `src/macropore/macropore.f90` — FrArMtrx reads cutover (S-2.9)
- `src/heat/temperature.f90`, `src/heat/frozencond.f90`, `src/utils/soilhydraulicsutils.f90` — reads cutover + hconduc sentinel (S-2.2, S-2.10)
- `src/io/swapoutput.f90`, `src/io/swap_csv_output.f90` — reads cutover + mini-sim writeback (S-2.11)
- `src/core/variables.f90` — retire ~74 globals (S-2.12)
- `docs/adr/0038-state-migration-soilwater-core.md` — ADR (S-2.13)
- `docs/superpowers/specs/state-migration-playbook.md` — playbook update (S-2.13)

---

# Phase 1 — State type extension + soilwater_init growth + dual-write home

### Task S-1.1: Create `soilwater_intermediate_t` + `soilwater_cumulative_t` cohort types; extend `soilwater_state_t`; pFUnit lifecycle tests

**Design decisions:** D2, D3
**Files:**
- Modify: `src/state/soilwater_state.f90`
- Modify: `tests/unit/state/test_soilwater_state.pf`

Define `soilwater_intermediate_t` (22 fields: 3 per-node arrays + 13 non-per-day scalars + 2 per-day arrays + 6 per-day scalars; two type-bound procedures: `reset()` and `reset_per_day()`). Define `soilwater_cumulative_t` (14 scalar fields, one `reset()` procedure). Extend `soilwater_state_t` with 30 flat-instantaneous fields (17 per-node arrays + 13 scalars) and the two cohort sub-records. **Do not extend `soilwater_init` body yet** (that is S-1.2). The new type components are declared with `= 0.0_real64` defaults for scalars and unallocated for arrays.

Field layout for the cohort types:

```fortran
type :: soilwater_intermediate_t
   ! per-node arrays (allocated in soilwater_init, S-1.2)
   real(real64), allocatable :: inq(:), inqrot(:), inqssdi(:)
   real(real64), allocatable :: iqdo(:), iqup(:), IThetaBeg(:)
   ! non-per-day scalars (16)
   real(real64) :: iqrot     = 0.0_real64
   real(real64) :: iqssdi    = 0.0_real64
   real(real64) :: iqredwet  = 0.0_real64
   real(real64) :: iqreddry  = 0.0_real64
   real(real64) :: iqredsol  = 0.0_real64
   real(real64) :: iqredfrs  = 0.0_real64
   real(real64) :: ies0      = 0.0_real64
   real(real64) :: iet0      = 0.0_real64
   real(real64) :: iew0      = 0.0_real64
   real(real64) :: iintc     = 0.0_real64
   real(real64) :: iruno     = 0.0_real64
   real(real64) :: irunoCN   = 0.0_real64
   real(real64) :: irunon    = 0.0_real64
   real(real64) :: iqbot     = 0.0_real64
   real(real64) :: iqtdo     = 0.0_real64
   real(real64) :: iqtup     = 0.0_real64
   real(real64) :: IPondBeg  = 0.0_real64
   real(real64) :: iprec     = 0.0_real64
   real(real64) :: igird     = 0.0_real64
   real(real64) :: inird     = 0.0_real64
   ! per-day cohort (flDayStart gate — reset_per_day() only, NOT reset())
   real(real64) :: tra           = 0.0_real64
   real(real64) :: iqredwet_day  = 0.0_real64
   real(real64) :: iqreddry_day  = 0.0_real64
   real(real64) :: iqredsol_day  = 0.0_real64
   real(real64) :: iqredfrs_day  = 0.0_real64
   real(real64) :: iptra_day     = 0.0_real64
   real(real64), allocatable :: qpotrot_day(:), qredtot_day(:)
contains
   procedure :: reset         => soilwater_intermediate_reset
   procedure :: reset_per_day => soilwater_per_day_reset
end type soilwater_intermediate_t

type :: soilwater_cumulative_t
   real(real64) :: cqssdi    = 0.0_real64
   real(real64) :: cqrot     = 0.0_real64
   real(real64) :: cqbot     = 0.0_real64
   real(real64) :: cqbotdo   = 0.0_real64
   real(real64) :: cqbotup   = 0.0_real64
   real(real64) :: cinund    = 0.0_real64
   real(real64) :: crunon    = 0.0_real64
   real(real64) :: crunoff   = 0.0_real64
   real(real64) :: crunoffCN = 0.0_real64
   real(real64) :: cqtdo     = 0.0_real64
   real(real64) :: cqtup     = 0.0_real64
   real(real64) :: cqprai    = 0.0_real64
   real(real64) :: cgird     = 0.0_real64
   real(real64) :: cnird     = 0.0_real64
contains
   procedure :: reset => soilwater_cumulative_reset
end type soilwater_cumulative_t
```

`reset()` zeroes all fields including per-day fields and uses `if (allocated(self%X)) self%X = 0.0_real64` guards for arrays. `reset_per_day()` zeroes only the 8 per-day fields (`tra`, `iqredwet_day`, `iqreddry_day`, `iqredsol_day`, `iqredfrs_day`, `iptra_day`, `qpotrot_day(:)`, `qredtot_day(:)`).

Extend `soilwater_state_t` declaration with the 30 flat-instantaneous fields (declare them; no allocate yet) and the two cohort components. Update the `public` export list in `soilwater_state_mod` to include `soilwater_intermediate_t`, `soilwater_cumulative_t`.

- [ ] **Add `soilwater_intermediate_t` type definition** (as above) to `soilwater_state.f90` above `soilwater_state_t`.
- [ ] **Add `soilwater_cumulative_t` type definition** (as above).
- [ ] **Extend `soilwater_state_t`** with all 30 flat-instantaneous declarations + `type(soilwater_intermediate_t) :: intr` + `type(soilwater_cumulative_t) :: cumu`. Use zero-defaults for scalars; `logical :: fllowgwl = .false.` for the flag.
- [ ] **Update `public` list** in the module to export the two new cohort types.
- [ ] **Implement `soilwater_intermediate_reset`** in the `contains` block — zeroes all 22 fields with `allocated()` guards for arrays.
- [ ] **Implement `soilwater_per_day_reset`** — zeroes only the 8 per-day fields with `allocated()` guards for the two arrays.
- [ ] **Implement `soilwater_cumulative_reset`** — zeroes all 14 fields (no arrays; no guards needed).
- [ ] **Extend `test_soilwater_state.pf`** with ~10 new tests:
  - Fresh `soilwater_state_t()` defaults: `pond`, `gwl`, `hatm` at their declared defaults; `fllowgwl = .false.`
  - `intr%reset()`: seed all non-per-day scalars with non-zero; call `reset()`; assert all 22 fields are zero
  - `intr%reset_per_day()`: seed ALL intr fields non-zero; call `reset_per_day()`; assert per-day fields zero, non-per-day unchanged
  - `cumu%reset()`: seed all 14 fields non-zero; call `reset()`; assert all zero
  - `allocated()` guards: allocate `intr%inq(3)`, set non-zero, call `intr%reset()`, assert zero; repeat for `cumu` (no arrays — N/A)
  - No aliasing between two independent `soilwater_state_t` instances
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
feat(state): SS-SWC Phase 1 S-1.1 — soilwater cohort types + state-type extension

New soilwater_intermediate_t: 22 fields (6 per-node arrays + 16 non-per-day
scalars), type-bound reset() (flzerointr) + reset_per_day() (flDayStart).
New soilwater_cumulative_t: 14 scalar fields, type-bound reset() (flzerocumu).
Mild ADR 0033 extension: two distinct reset procedures on intr cohort.
soilwater_state_t extended with 30 flat-instantaneous fields + intr + cumu.
pFUnit ~10 new tests green. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-1.2: Extend `soilwater_init` body — allocate 24 new per-node/per-layer arrays + cohort arrays + scalar inits

**Design decisions:** D4
**Files:**
- Modify: `src/state/soilwater_state.f90`

Extend the body of `soilwater_init(sw, numnod, nlay)` to allocate and zero all new fields. Signature unchanged.

New allocations (all sized by `numnod` unless noted):
- Per-node F-INST arrays (16): `theta`, `thetm1`, `thetar`, `thetas`, `h`, `hm1`, `q(numnod+1)`, `k(numnod+1)`, `kmean(numnod+1)`, `dimoca`, `cofgen(21, numnod)`, `FrArMtrx`, `fluseksatexm(logical)`, `indeks(integer)`, `evp`
- Per-layer F-INST: `thetsl(nlay)`
- Intermediate cohort arrays (6): `intr%inq(numnod+1)`, `intr%inqrot(numnod)`, `intr%inqssdi(numnod)`, `intr%iqdo(numnod+1)`, `intr%iqup(numnod+1)`, `intr%IThetaBeg(numnod)`
- Per-day cohort arrays (2): `intr%qpotrot_day(numnod)`, `intr%qredtot_day(numnod)`

New scalar inits: `pond = 0.0`, `pondm1 = 0.0`, `pondini = 0.0`, `gwl = 0.0`, `gwlm1 = 0.0`, `pegwl = 999.0` (legacy default), `gwlflcpzo = 999.0`, `hatm = -2.75d5` (legacy init from soilhydraulics:870), `volact = 0.0`, `volm1 = 0.0`, `volini = 0.0`, `wbalance = 0.0`, `runon = 0.0`, `fllowgwl = .false.`, `nodgwl = 0`, `bpegwl = -1`, `npegwl = -1`, `nodgwlflcpzo = 0`.

Call `call sw%intr%reset()` and `call sw%cumu%reset()` at end of init to ensure cohort scalar fields are also explicitly zeroed (forward-compat pattern from atmosphere_init).

- [ ] **Add F-INST per-node allocations and zero** (theta, thetm1, thetar, thetas, h, hm1; q/k/kmean as numnod+1; dimoca, cofgen(21,numnod), FrArMtrx, fluseksatexm, indeks, evp). Use `allocate(); sw%X = 0.0_real64` pattern.
- [ ] **Add F-INST per-layer:** `allocate(sw%thetsl(nlay)); sw%thetsl = 0.0_real64`.
- [ ] **Add scalar inits** for the 19 new scalars/integers/logicals.
- [ ] **Add cohort array allocations** (intr%inq, intr%inqrot, intr%inqssdi, intr%iqdo, intr%iqup, intr%IThetaBeg, intr%qpotrot_day, intr%qredtot_day).
- [ ] **Add cohort reset calls** at end of body.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
feat(state): SS-SWC Phase 1 S-1.2 — soilwater_init: allocate 24 new arrays + scalar inits

soilwater_init extended with 16 per-node F-INST arrays (including q/k/kmean
as numnod+1; cofgen(21,numnod); logical fluseksatexm; integer indeks),
1 per-layer thetsl(nlay), 8 intr cohort arrays, 19 new scalars/integers/
logicals. Signature (sw, numnod, nlay) unchanged. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-1.3: `soilhydraulics.f90` SoilWater(task=1) init dual-write — startup seed for ~30 fields

**Design decisions:** D2
**Files:**
- Modify: `src/soil/soilhydraulics.f90`

`SoilWater(task=1)` at `soilhydraulics.f90:842` initialises the Richards domain: sets `cofgen` from PdmVG tables (lines 891–937), initial h profile, hatm (line 870), FrArMtrx (1051, 1062), thetar/thetas/thetsl (948–952), indeks/fluseksatexm (955–960). Add dual-writes at all init assignment sites using an ASSOCIATE block: `associate(sw => state%soilwater)`.

- [ ] **Grep all init-write targets in SoilWater(1) body:**
  ```bash
  grep -n "cofgen\|hatm\|FrArMtrx\|thetar\|thetas\|thetsl\|indeks\|fluseksatexm\|theta\b\|h(\|runon\|pond\b" \
    src/soil/soilhydraulics.f90 | head -60
  ```
- [ ] **Add ASSOCIATE block** at top of task=1 body: `associate(sw => state%soilwater)`.
- [ ] **Dual-write each init assignment** (form: `cofgen(k,i) = X; sw%cofgen(k,i) = X`).
- [ ] **Dual-write `hatm = -2.75d5`** (line 870): `sw%hatm = -2.75e5_real64`.
- [ ] **Dual-write `runon = 0.d0`** (line 878): `sw%runon = 0.0_real64`.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 1 S-1.3 — soilhydraulics SoilWater(1) init dual-write

Dual-writes for ~30 init-path fields: cofgen (~28 sites in 891-937 block),
hatm, FrArMtrx (1051+1062), thetar/thetas/thetsl (948-952), indeks/
fluseksatexm (955-960), runon (878). ASSOCIATE sw => state%soilwater.
check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-1.4a: `soilhydraulics.f90` `headcalc` dual-write — theta + h (per-step Richards solve, ~40 sites)

**Design decisions:** D2
**Files:**
- Modify: `src/soil/soilhydraulics.f90`

`headcalc` (lines 25–840) is the Richards solver; `h` and `theta` are written at every iteration. This is split from S-1.4b to keep each commit reviewable. Also in this task: `SoilWaterStateVar(task, state)` signature addition (add `state` arg; retarget hm1/thetm1/pondm1/gwlm1 saves in task=1 and h/theta/pond/gwl/kmean(numnod+1) restores in task=2). `hysteresis(state)` signature addition (retarget thetar/thetas/cofgen/dimoca/h writes inside hysteresis body; ~10 sites).

- [ ] **Add `state` arg to `SoilWaterStateVar(task, state)`** — signature change + retarget task=1 writes (`hm1`, `thetm1`, `gwlm1`, `pondm1` → `sw%hm1`, etc.) and task=2 restores + the `kmean(numnod+1) = k(numnod)` retarget. Update 3 callers: `soilhydraulics.f90:soilwater:1161`, `swapoutput.f90:3903`, `swapoutput.f90:3915`. Add dual-writes (legacy side stays for now).
- [ ] **Add `state` arg to `hysteresis(state)`** — retarget ~10 write sites for `thetar`/`thetas`/`cofgen`/`dimoca`/`h` to `state%soilwater%X`. Add dual-writes.
- [ ] **Grep theta + h write sites in headcalc:**
  ```bash
  grep -n "\btheta(\|h(\|theta =\|h =" src/soil/soilhydraulics.f90 | grep -v "state%\|!" | head -50
  ```
- [ ] **Add ASSOCIATE in headcalc body** (or extend existing) to include `sw => state%soilwater`.
- [ ] **Dual-write `theta(:)` write sites** (~25 sites in headcalc: 1043, 237, 443, 519, 715, 1042, 1240 per discovery + iteration writes).
- [ ] **Dual-write `h(:)` write sites** (~15 sites: 123, 125, 435, 439, 724, 1037 + hysteresis path).
- [ ] **Dual-write `thetm1`/`hm1` (save-state)** via SoilWaterStateVar task=1 path (already handled above).
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 1 S-1.4a — headcalc theta/h dual-write; SoilWaterStateVar/hysteresis plumbing

SoilWaterStateVar gains state arg (3 callers updated). hysteresis gains state
arg; ~10 write sites dual-written. headcalc: theta (~25 sites) and h (~15 sites)
dual-written. Save-state path (hm1/thetm1/pondm1/gwlm1) dual-written via
SoilWaterStateVar task=1. ASSOCIATE sw => state%soilwater. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-1.4b: `soilhydraulics.f90` `headcalc` dual-write — q + k + kmean + dimoca + cohort reset dual-patch

**Design decisions:** D2, D5
**Files:**
- Modify: `src/soil/soilhydraulics.f90`

Continuation of headcalc dual-write for the remaining per-step fields.

- [ ] **Dual-write `q(:)` write sites** (~10 sites: 881 init-zero, soilwater(2/3) flux assignments).
- [ ] **Dual-write `k(:)` write sites** (~12 sites: 130, 170, 176, 181, 238, 453, 467, 471, 520, 1052, 1172, 1176).
- [ ] **Dual-write `kmean(:)` write sites** (~14 sites: 113, 133, 186, 189, 241, 262, 456, 474, 475, 523, 548, 1055, 1178, 1243).
- [ ] **Dual-write `dimoca(:)` write sites** (~3 sites: 312, 1049 + hysteresis:1339).
- [ ] **Dual-write the flzerointr reset block (soilhydraulics:1083–1132)** — for each field being zeroed (inq, inqrot, inqssdi, IThetaBeg, IPondBeg, iqrot, iqssdi, iqredwet/dry/sol/frs, ies0/iet0/iew0, iintc, iruno/irunon/irunoCN, iqbot/iqtdo/iqtup, iqdo/iqup), add `state%soilwater%intr%<field> = 0.0_real64` immediately after the legacy zero. Also dual-patch the per-day reset block (iqredwet_day etc., lines 1083–1092).
- [ ] **Dual-write the flzerocumu reset block (soilhydraulics:1134–1158)** — same dual-patch pattern for all 14 cumu fields + pondini/volini rebase writes.
- [ ] **Dual-write `IThetaBeg` assign on flzerointr** (line 1127): `sw%intr%IThetaBeg(:) = sw%theta(:)`.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 1 S-1.4b — headcalc q/k/kmean/dimoca dual-write; cohort reset dual-patch

headcalc: q (~10 sites), k (~12), kmean (~14), dimoca (~3) dual-written.
flzerointr reset block (1083-1132, ~30 fields) dual-patched to state%soilwater%intr.
flzerocumu reset block (1134-1158, 14 cumu fields + pondini/volini rebase)
dual-patched to state%soilwater%cumu. Inline resets remain (S-2.1 removes them).
check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-1.5: `waterbalance.f90` `integral` cumulative-cohort dual-write — ~50 accumulator write sites

**Design decisions:** D2
**Files:**
- Modify: `src/soil/waterbalance.f90`

`integral(state)` at `waterbalance.f90:379` is the primary accumulator for all intermediate + cumulative cohort fields. Already has `state` in scope. Dual-write every accumulate-write site.

- [ ] **Grep all accumulator write sites in `integral`:**
  ```bash
  grep -n "inq\|iqrot\|iqssdi\|iqred\|ies0\|iet0\|iew0\|iintc\|iruno\|iqbot\|iqtdo\|iqtup\|iqdo\|iqup\|iprec\|igird\|inird\|IPondBeg\|IThetaBeg\|tra\|iptra_day\|qpotrot_day\|qredtot_day\|cqssdi\|cqrot\|cqbot\|cqbotdo\|cqbotup\|cinund\|crunon\|crunoff\|crunoffCN\|cqtdo\|cqtup\|cqprai\|cgird\|cnird" \
    src/soil/waterbalance.f90 | grep -v "state%\|!" | head -80
  ```
- [ ] **Add ASSOCIATE block** at top of `integral` body: `associate(sw => state%soilwater)`.
- [ ] **Dual-write all ~50 accumulator write sites** — form: `cqrot = cqrot + X; sw%cumu%cqrot = sw%cumu%cqrot + X`. For the intermediate cohort: `iqrot = iqrot + X; sw%intr%iqrot = sw%intr%iqrot + X`. Per-day: `tra = tra + X; sw%intr%tra = sw%intr%tra + X`. Arrays: `inq(i) = inq(i) + X; sw%intr%inq(i) = sw%intr%inq(i) + X`.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 1 S-1.5 — waterbalance integral dual-write (~50 accumulator sites)

integral(state): ~50 accumulator write sites dual-written across intr cohort
(inq, inqrot, inqssdi, iqrot, iqssdi, iqredwet/dry/sol/frs, ies0/iet0/iew0,
iintc, iruno/irunon/irunoCN, iqbot/iqtdo/iqtup, iqdo/iqup, iprec/igird/inird,
per-day tra/iptra_day/qpotrot_day/qredtot_day, IPondBeg/IThetaBeg) and cumu
cohort (cqssdi/cqrot/cqbot/cqbotdo/up, cinund/crunon/crunoff/crunoffCN,
cqtdo/cqtup/cqprai/cgird/cnird). ASSOCIATE sw block. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-1.6: `waterbalance.f90` fluxes/calcgwl/watstor/checkmassbal — pond + gwl + volume dual-write; calcgwl intent flip

**Design decisions:** D7, D9 (watstor), D10 (mini-sim support)
**Files:**
- Modify: `src/soil/waterbalance.f90`

Four distinct sub-tasks in one file:
1. `calcgwl(state)` intent-flip: change `intent(in)` to `intent(inout)`; add dual-writes for ~12 write sites (gwl, nodgwl, pegwl, bpegwl, npegwl, gwlflcpzo, nodgwlflcpzo).
2. `watstor()` plumbing: add `state` arg; retarget volact + volm1 writes (3 sites). Update caller.
3. `fluxes(state)` scalar reads — verify q is available via state (already after S-1.4b).
4. `checkmassbal` — verify pond/gwl reads will be satisfied by dual-writes from S-1.3/S-1.4a.

- [ ] **Flip `calcgwl(state)` to `intent(inout)`** — add dual-writes for gwl, nodgwl, pegwl, bpegwl, npegwl, gwlflcpzo, nodgwlflcpzo at all ~12 write sites in calcgwl (lines 53–157). Form: `gwl = X; state%soilwater%gwl = X`.
- [ ] **Add `state` arg to `watstor(state)`** — retarget `volact` (lines 830–833) and `volm1` (line 830). Update caller.
- [ ] **Add `state` arg to `level(state, swoptlev, ...)` and `watertable(state, node, ...)`** — these are pure-read helpers; verify they have no writes needing dual-write; just add state arg for forward-compat.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 1 S-1.6 — waterbalance calcgwl intent-flip + watstor plumbing

calcgwl: intent(in) -> intent(inout); ~12 write sites dual-written (gwl,
nodgwl, pegwl, bpegwl, npegwl, gwlflcpzo, nodgwlflcpzo). watstor gains
state arg; volact/volm1 (3 sites) dual-written. level/watertable gain
state arg (pure reads). check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-1.7: `boundtop.f90` pond + kmean(1) dual-write (8 + 3 sites)

**Design decisions:** D6, D8
**Files:**
- Modify: `src/boundary/boundtop.f90`

`boundtop.f90` already takes `state intent(inout)` (ADR 0035). Dual-write `pond` at 8 sites (147, 166, 185, 258, 266, 274, 290, 308) and `kmean(1)` at 3 sites. Also dual-write `hatm` read site (boundtop reads `hatm` — this becomes `state%soilwater%hatm` after Phase 2; add dual-write source on the write side if any writes exist, else mark for reader cutover).

- [ ] **Add ASSOCIATE block** at relevant boundtop entry points: `associate(sw => state%soilwater)`.
- [ ] **Dual-write `pond`** at all 8 write sites: `pond = X; sw%pond = X`.
- [ ] **Dual-write `kmean(1)`** at 3 sites: `kmean(1) = X; sw%kmean(1) = X`.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 1 S-1.7 — boundtop pond + kmean(1) dual-write

boundtop.f90: pond dual-written at 8 sites (147,166,185,258,266,274,290,308);
kmean(1) dual-written at 3 sites. ASSOCIATE sw => state%soilwater.
Already state-plumbed (ADR 0035) — pure co-write retarget. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-1.8: `boundbottom.f90` kmean(numnod+1) dual-write (2 sites)

**Design decisions:** D8
**Files:**
- Modify: `src/boundary/boundbottom.f90`

`boundbottom.f90` already takes `state intent(inout)` (ADR 0035). Retarget `kmean(numnod+1)` writes at 2 sites.

- [ ] **Add dual-write for `kmean(numnod+1)`** at 2 sites (boundbottom:171 and one additional site per discovery). Form: `kmean(numnod+1) = X; state%soilwater%kmean(numnod+1) = X`.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 1 S-1.8 — boundbottom kmean(numnod+1) dual-write

boundbottom.f90: kmean(numnod+1) dual-written at 2 sites. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-1.9: `tillage.f90` pond + cofgen + theta + h dual-write — thread state into Change_MvGpars + Adapt_WC_H

**Design decisions:** D6, D2
**Files:**
- Modify: `src/crop/tillage.f90`

`DoTillage(iTask, state)` already takes `state` (A-2.6 windfall from ADR 0037). Internal subs `Change_MvGpars`, `Adapt_WC_H` read/write `cofgen`, `theta`, `h`, `pond` via bare `use Variables`. Thread state into these two subs.

Init-order note (discovery hazard #8): `DoTillage(1, state)` runs at `swap.f90:203` BEFORE `SoilWater(1, state)` at line 207. SoilWater(1) subsequently reinitialises cofgen. After migration, both write through `state%soilwater%cofgen`; the SoilWater(1) re-write preserves legacy semantics — verify the `cofgen = 0.0` init-blank in SoilWater(1) happens AFTER `DoTillage(1)`. The dual-write in tillage covers the DoTillage(1) side; the S-1.3 dual-write covers the SoilWater(1) side.

- [ ] **Add `state intent(inout)` arg to `Change_MvGpars`** and `Adapt_WC_H` (module-internal subs of tillage). Update the call sites inside `DoTillage`.
- [ ] **Dual-write `cofgen` writes in `Change_MvGpars`** — form: `cofgen(k,i) = X; state%soilwater%cofgen(k,i) = X`.
- [ ] **Dual-write `theta`, `h`, `pond` writes in `Adapt_WC_H`** — including `pond` write at tillage:304,330.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full` (including init-order regression check — run the full 5-case suite).

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 1 S-1.9 — tillage Change_MvGpars/Adapt_WC_H state threading + dual-write

Change_MvGpars and Adapt_WC_H gain state intent(inout) args. Dual-writes
for cofgen (Change_MvGpars), theta + h + pond (Adapt_WC_H, tillage:304+330).
Init-order DoTillage(1)->SoilWater(1) verified: SoilWater(1) re-init of cofgen
happens after DoTillage(1); legacy semantics preserved. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-1.10: `rootextraction.f90` h-clamp dual-write + verify cumu/intr links for qrot family

**Design decisions:** D2
**Files:**
- Modify: `src/crop/rootextraction.f90`

`rootextraction.f90` co-writes `h(node)` at a clamp site (~line 382). `qrot`/`qrosum` are already in state (ADR 0036). The `iqrot`/`cqrot`/`inqrot` accumulators are written in `integral` (covered by S-1.5). This task: dual-write the h-clamp.

- [ ] **Dual-write `h(node)` clamp** at rootextraction:382: `h(node) = X; state%soilwater%h(node) = X`.
- [ ] **Verify cumu/intr links** — grep `iqrot\|cqrot\|inqrot\|qredwet\|iqredwet` in rootextraction.f90 to confirm no direct writes (all writes go through integral; if any direct writes exist, dual-write them here).
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 1 S-1.10 — rootextraction h-clamp dual-write

rootextraction.f90: h(node) clamp at ~line 382 dual-written. iqrot/cqrot/
inqrot writes confirmed in waterbalance integral only — no direct writes
in rootextraction. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-1.11: `macropore.f90` FrArMtrx dual-write (2 sites) + soilhydraulics FrArMtrx non-macropore fallback

**Design decisions:** D2
**Files:**
- Modify: `src/macropore/macropore.f90`
- Modify (already touched): `src/soil/soilhydraulics.f90`

`macropore.f90:435,550` writes `FrArMtrx(:)` when macropore is active. `soilhydraulics.f90:1051,1062` writes `FrArMtrx(i) = 1.d0` as the non-macropore fallback (already covered in S-1.3 but verify). macropore.f90 already accepts state for tasks 1–6.

- [ ] **Dual-write `FrArMtrx(:)` at macropore:435,550**: `FrArMtrx(:) = X; state%soilwater%FrArMtrx(:) = X`.
- [ ] **Confirm soilhydraulics:1051,1062 dual-writes** were added in S-1.3 (verify grep). Add if missed.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 1 S-1.11 — macropore FrArMtrx dual-write (2 sites)

macropore.f90:435+550 FrArMtrx writes dual-written to state%soilwater%FrArMtrx.
soilhydraulics non-macropore fallback (1051+1062) confirmed dual-written (S-1.3).
check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-1.12: `frozencond.f90` q + theta verify + dual-write under frost path

**Design decisions:** D2
**Files:**
- Modify: `src/heat/frozencond.f90`

`frozencond.f90` reads `theta`, `thetas`, `gwl` and reads (but does not write) `qbot` (already in state per ADR 0035). Discovery Section 4 says "reads only — no new co-write." Verify no write sites for theta/q/gwl, then dual-write any surprise writes.

- [ ] **Grep frozencond.f90 for write sites:**
  ```bash
  grep -n "theta\|q(\|gwl" src/heat/frozencond.f90 | grep -v "state%\|!" | head -20
  ```
- [ ] **If write sites exist:** dual-write them following S-1.4 pattern.
- [ ] **If no write sites:** one documentation comment confirming read-only status + verify PASS.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 1 S-1.12 — frozencond.f90 write-site audit

frozencond.f90 confirmed read-only for theta/q/gwl (qbot already state per
ADR 0035). [OR: N write sites found and dual-written.] check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

# Phase 2 — Reader cutover + cohort consolidation + retire

### Task S-2.1: Cohort reset consolidation — replace ~55 scattered reset lines with 3 type-bound calls

**Design decisions:** D5
**Files:**
- Modify: `src/soil/soilhydraulics.f90`

Apply the atmosphere ADR 0037 D5 pattern. Replace the now-dual-patched `soilhydraulics.f90:1083–1132` (flzerointr block) and `soilhydraulics.f90:1134–1158` (flzerocumu block) with clean cohort reset() calls. The `pondini = pond` and `volini = volact` rebases immediately after the flzerocumu block must be retained as inline non-zero rebases (ADR 0033 pattern — not zeroing, so they stay).

**Pre-flight check before this task:** confirm all INTR fields have non-zero dual-write coverage:
```bash
grep -rn "state%soilwater%intr%\s*=" src/ | grep -v "0\.0\|0_real64\|reset"
# Must return results for inq, inqrot, iqrot, etc.
```

- [ ] **Replace `flzerointr` block (1083–1132)** with:
  ```fortran
  if (flzerointr) then
     call state%soilwater%intr%reset()
  end if
  ```
- [ ] **Replace per-day reset sub-block** with `if (flDayStart) call state%soilwater%intr%reset_per_day()` at the canonical per-day reset site.
- [ ] **Replace `flzerocumu` block (1134–1158)** with:
  ```fortran
  if (flzerocumu) then
     call state%soilwater%cumu%reset()
     ! non-zero rebases: volini = volact (state-side); pondini = pond (state-side)
     state%soilwater%volini  = state%soilwater%volact
     state%soilwater%pondini = state%soilwater%pond
  end if
  ```
- [ ] **Remove legacy global zeros** from both former reset blocks (the dual-write legacy side is now orphaned — remove).
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`. ~55 lines collapse to ~8.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 2 S-2.1 — cohort reset consolidation (~55 lines -> 3 calls)

soilhydraulics.f90 flzerointr block (1083-1132) replaced with
state%soilwater%intr%reset(). flDayStart per-day block replaced with
intr%reset_per_day(). flzerocumu block (1134-1158) replaced with
cumu%reset() + inline pondini/volini rebases. ~55 scattered legacy zeros
removed. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-2.2: `hconduc` `tsoil_node` sentinel fix — 14-site coordinated update

**Design decisions:** D9
**Files:**
- Modify: `src/soil/soilhydraulics.f90` (~10 sites)
- Modify: `src/boundary/boundbottom.f90` (1 site)
- Modify: `src/boundary/boundtop.f90` (1 site)
- Modify: `src/crop/rootextraction.f90` (2 sites)
- Modify: `src/crop/tillage.f90` (1 site)
- Modify: `src/io/swapoutput.f90` (1 site)
- Modify: `src/macropore/macropore.f90` (1 site)

Replace all 14+ `hconduc(...)` calls that omit `tsoil_node` or pass literals (`0.0d0`, `10.d0`, `1.0d0`) with `state%heat%tsoil(node)`. All callers have `state` in scope post-arc-7.

- [ ] **Grep all hconduc call sites:**
  ```bash
  grep -rn "hconduc(" src/ --include="*.f90" | grep -v "state%\|!"
  ```
- [ ] **For each call site:** add or replace the `tsoil_node` argument with `state%heat%tsoil(node)` (or `state%heat%tsoil(i)` as appropriate for the loop variable). At `rootextraction.f90:868,873` replace `10.d0` literal; at `tillage.f90:361` replace `1.0d0` literal; at `macropore.f90:1142` replace `Dum`.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
fix(physics): SS-SWC Phase 2 S-2.2 — hconduc tsoil_node sentinel: replace 0.0/literal with state%heat%tsoil(node)

14+ hconduc call sites across soilhydraulics (~10), boundbottom, boundtop,
rootextraction (10.d0 literal), tillage (1.0d0 literal), swapoutput,
macropore (Dum) updated to pass state%heat%tsoil(node). Resolves heat
ADR 0034 Task 9 residual (tsoil_loc=0.0 sentinel). check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-2.3: `soilhydraulics.f90` reads cutover — theta/h/q/k/kmean/dimoca/cofgen/FrArMtrx/hatm/pond/gwl/volact/wbalance

**Design decisions:** D2, D6, D7
**Files:**
- Modify: `src/soil/soilhydraulics.f90`

Pre-flight dual-write coverage check for all fields before starting. This is the largest single reader-cutover task; split into sub-passes if needed during execution.

- [ ] **Pre-flight check** for each field family (theta, h, q, k, kmean, dimoca, cofgen, FrArMtrx, hatm, pond, gwl, volact, wbalance, nodgwl, pegwl, hm1, thetm1, etc.):
  ```bash
  grep -rn "state%soilwater%<field>\s*=" src/ | grep -v "0\.0\|0_real64\|reset"
  # Confirm non-zero writes exist for all fields
  ```
- [ ] **Migrate all reads** in headcalc/soilwater/hysteresis/SoilWaterStateVar bodies from bare globals to `state%soilwater%X`. ASSOCIATE block: `associate(sw => state%soilwater)`.
- [ ] **Drop bare-global read sides** from `use variables, only:` clauses for migrated fields. Grid dims (`numnod`, `dz`, etc.) stay in `use variables`.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 2 S-2.3 — soilhydraulics reads cutover (~150 sites)

Reads migrated for theta/h/q/k/kmean/dimoca/cofgen/FrArMtrx/hatm/pond/
gwl/volact/wbalance/hm1/thetm1/pondm1/gwlm1/nodgwl/pegwl/fllowgwl etc.
Dual-writes dropped from home-tree write side. ASSOCIATE sw => state%soilwater.
Grid dims (numnod/dz/z/layer/etc.) retained in use variables. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-2.4: `waterbalance.f90` reads cutover — integral + fluxes + calcgwl + checkmassbal

**Design decisions:** D2, D7
**Files:**
- Modify: `src/soil/waterbalance.f90`

- [ ] **Pre-flight check** (per playbook atmosphere lesson) for all waterbalance-consumed fields.
- [ ] **Migrate reads in `integral`** — all cohort-field read-sides now point at state (they already point at state after S-1.5 dual-writes; drop legacy global side, update the read expressions to read from state exclusively).
- [ ] **Migrate reads in `fluxes`** — q reads from `state%soilwater%q(i)`.
- [ ] **Migrate reads in `calcgwl`** — gwl, h, theta, z (z stays legacy), nodgwl, etc.
- [ ] **Migrate reads in `checkmassbal`** — already reads atmosphere state; add soilwater reads.
- [ ] **Drop dual-write legacy global write side** from integral + calcgwl accumulator paths.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 2 S-2.4 — waterbalance reads cutover; dual-writes dropped

integral: all cohort accumulator reads migrated; dual-write legacy sides dropped.
fluxes: q reads migrated. calcgwl: gwl/nodgwl/pegwl write side already state
(S-1.6); read sides also migrated. checkmassbal: soilwater reads migrated.
check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-2.5: `boundtop.f90` + `boundbottom.f90` reads cutover (pond/gwl/kmean readbacks + hatm)

**Design decisions:** D6, D7, D8
**Files:**
- Modify: `src/boundary/boundtop.f90`
- Modify: `src/boundary/boundbottom.f90`

- [ ] **Migrate `hatm` read in boundtop** to `state%soilwater%hatm`.
- [ ] **Migrate `pond`/`pondm1` reads in boundtop** — all 8 write sites already point at state (S-1.7); drop dual-write legacy write side; migrate read side.
- [ ] **Migrate `kmean(1)` reads in boundtop** — drop dual-write; migrate reads.
- [ ] **Migrate `kmean(numnod+1)` reads in boundbottom** — drop dual-write; migrate reads.
- [ ] **Migrate `gwl`/`theta`/`h` reads in boundbottom** (for hbot calculation).
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 2 S-2.5 — boundtop/boundbottom reads cutover + dual-writes dropped

boundtop: hatm, pond/pondm1, kmean(1) reads migrated; pond/kmean dual-writes
dropped. boundbottom: kmean(numnod+1), gwl, theta, h reads migrated.
All boundary co-writes now state-exclusive. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-2.6: `tillage.f90` reads cutover

**Design decisions:** D2, D6
**Files:**
- Modify: `src/crop/tillage.f90`

- [ ] **Migrate `cofgen` reads in `Change_MvGpars`** to `state%soilwater%cofgen`.
- [ ] **Migrate `theta`/`h`/`pond` reads in `Adapt_WC_H`** to state; drop dual-write legacy sides from S-1.9.
- [ ] **Migrate any remaining bare-global reads** for dz/layer/z (these stay as legacy — grid dims).
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 2 S-2.6 — tillage reads cutover; cofgen/theta/h/pond dual-writes dropped

Change_MvGpars: cofgen reads migrated. Adapt_WC_H: theta/h/pond reads migrated.
Dual-write legacy sides dropped. Grid dims (dz/layer) retained as legacy.
check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-2.7: `rootextraction.f90` + `cropgrowth.f90` + `oxygenstress.f90` reads cutover

**Design decisions:** D2
**Files:**
- Modify: `src/crop/rootextraction.f90`
- Modify: `src/crop/cropgrowth.f90`
- Modify: `src/crop/oxygenstress.f90`

- [ ] **rootextraction.f90:** migrate theta, h, thetar, thetas, hm1, k, kmean reads (~30 sites); drop h-clamp dual-write (S-1.10 legacy side).
- [ ] **cropgrowth.f90:** migrate theta, h reads (~15 sites).
- [ ] **oxygenstress.f90:** migrate theta, thetas, h, dimoca, cofgen reads (~10 sites).
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 2 S-2.7 — rootextraction/cropgrowth/oxygenstress reads cutover

rootextraction: theta/h/thetar/thetas/hm1/k/kmean (~30 sites). cropgrowth:
theta/h (~15 sites). oxygenstress: theta/thetas/h/dimoca/cofgen (~10 sites).
All already state-plumbed (ADR 0036). check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-2.8: `solute.f90` + `drainage.f90` + `surfacewater.f90` + related reads cutover

**Design decisions:** D2, D6, D7
**Files:**
- Modify: `src/solute/solute.f90`
- Modify: `src/solute/agetracer.f90`
- Modify: `src/drainage/drainage.f90`
- Modify: `src/drainage/surfacewater.f90`
- Modify: `src/drainage/divdra.f90`
- Modify: `src/utils/surfacewaterutils.f90`

- [ ] **solute.f90:** theta, h, q, inq(→intr%inq), gwl, FrArMtrx reads.
- [ ] **agetracer.f90:** theta, q, h, gwl, pond reads.
- [ ] **drainage.f90:** gwl, pond, h, theta reads.
- [ ] **surfacewater.f90:** gwl, pond, theta reads.
- [ ] **divdra.f90:** gwl, h, theta reads.
- [ ] **surfacewaterutils.f90:** pond, rsro, rsroexp reads (rsro/rsroexp are config — verify only pond is a soil-water global here).
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 2 S-2.8 — solute/drainage/surfacewater/divdra/surfacewaterutils reads cutover

solute: theta/h/q/inq(intr)/gwl/FrArMtrx. agetracer: theta/q/h/gwl/pond.
drainage/surfacewater/divdra: gwl/pond/h/theta. surfacewaterutils: pond.
check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-2.9: `macropore.f90` + `macrorate.f90` reads cutover; FrArMtrx dual-writes dropped

**Design decisions:** D2
**Files:**
- Modify: `src/macropore/macropore.f90`
- Modify: `src/macropore/macrorate.f90`

- [ ] **macropore.f90:** h, theta, cofgen, kmean reads cutover; FrArMtrx writes drop dual-write legacy side (S-1.11). cQMpLatSs init-zero at soilhydraulics:889 remains as legacy (D11 — macropore arc territory).
- [ ] **macrorate.f90:** theta, h reads cutover.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 2 S-2.9 — macropore/macrorate reads cutover; FrArMtrx dual-write dropped

macropore: h/theta/cofgen/kmean reads cutover; FrArMtrx dual-write legacy side
dropped. cQMpLatSs init-zero (soilhydraulics:889) retained as legacy (macropore
arc, D11). macrorate: theta/h reads cutover. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-2.10: `temperature.f90` + `frozencond.f90` + `soilhydraulicsutils.f90` reads cutover

**Design decisions:** D2
**Files:**
- Modify: `src/heat/temperature.f90`
- Modify: `src/heat/frozencond.f90`
- Modify: `src/utils/soilhydraulicsutils.f90`

- [ ] **temperature.f90:** theta, thetm1, thetas reads cutover (~6 sites).
- [ ] **frozencond.f90:** theta, thetas, gwl reads cutover (~8 sites; qbot already state-side ADR 0035).
- [ ] **soilhydraulicsutils.f90:** cofgen, numtab, sptab reads — numtab/sptab are tabulation globals staying legacy (D13); cofgen reads migrate to `state%soilwater%cofgen`; `fluseksatexm` read migrates. The `hconduc` sentinel fix was already done in S-2.2.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 2 S-2.10 — temperature/frozencond/soilhydraulicsutils reads cutover

temperature: theta/thetm1/thetas (~6 sites). frozencond: theta/thetas/gwl
(~8 sites). soilhydraulicsutils: cofgen/fluseksatexm reads migrated; numtab/
sptab retained as legacy tabulation (D13). check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-2.11: `swapoutput.f90` + `swap_csv_output.f90` reads cutover + mini-sim writeback retarget

**Design decisions:** D10, D2
**Files:**
- Modify: `src/io/swapoutput.f90`
- Modify: `src/io/swap_csv_output.f90`

Heaviest output task in this arc (~50 sites in swapoutput alone). Also: the mini-sim writeback (`swapoutput.f90:3865–3948`) must retarget `gwl`, `pond`, `theta(:)`, `h(:)` snapshot and restore to state. `qbot` is already state-side (ADR 0035 B-2.6).

- [ ] **Grep all soil-water fields in swapoutput.f90:**
  ```bash
  grep -n "theta\b\|h(\|gwl\|pond\|kmean\|cofgen\|FrArMtrx\|dimoca\|hm1\|volact\|volini\|cqrot\|cqbot\|iqrot\|wbalance\|hatm\|inq\b\|iqdo\|iqup\|crunoff\|crunon\|cinund" \
    src/io/swapoutput.f90 | grep -v "state%\|!" | head -80
  ```
- [ ] **Retarget mini-sim writeback at lines 3865–3948:**
  - Save: `gwltmp = state%soilwater%gwl`, `pondtmp = state%soilwater%pond`, `thetatmp(:) = state%soilwater%theta(:)`, `htmp(:) = state%soilwater%h(:)`.
  - Restore: `state%soilwater%gwl = gwltmp`, etc. Match the ADR 0035 `qbot` pattern exactly.
- [ ] **Migrate all output reads** to `state%soilwater%X`, `state%soilwater%intr%X`, `state%soilwater%cumu%X`. ASSOCIATE: `associate(sw => state%soilwater)`.
- [ ] **swap_csv_output.f90:** migrate theta, h, gwl, pond, kmean, FrArMtrx, q, inq(→intr%inq), iqbot(→intr%iqbot), iqrot(→intr%iqrot), irunon(→intr%irunon), iruno(→intr%iruno) reads.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 2 S-2.11 — swapoutput/swap_csv_output reads cutover + mini-sim retarget

swapoutput: ~50 soil-water read sites migrated (theta/h/q/gwl/pond/kmean/
cofgen/FrArMtrx/dimoca/volact/cumu cohort/intr cohort). Mini-sim writeback
(3865-3948) retargets gwl/pond/theta(:)/h(:) snapshot+restore to state%soilwater
(completing ADR 0035 B-2.6 — qbot was already retargeted). swap_csv_output:
~15 reads migrated. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-2.12: Drop dual-writes + retire ~74 globals from `variables.f90`; compile-driven hidden-reader fixups (expect 20–25)

**Design decisions:** D15
**Files:**
- Modify: `src/soil/soilhydraulics.f90` (drop remaining dual-write legacy sides)
- Modify: `src/soil/waterbalance.f90`
- Modify: `src/boundary/boundtop.f90`, `src/boundary/boundbottom.f90`
- Modify: `src/crop/tillage.f90`
- Modify: `src/core/variables.f90` — retire ~74 globals
- Modify: `src/core/initialize.f90` — drop redundant zero-inits

Pre-flight — for each of the ~74 globals, run the two-grep test:
```bash
# For each field, e.g. theta, h, gwl, pond, q, k, kmean, etc.:
grep -rEn "use variables.*\b<var>\b" src/ --include="*.f90" | grep -v "variables.f90\|initialize.f90"
grep -rEn "\b<var>\b" src/ --include="*.f90" | grep -v "variables.f90\|initialize.f90\|state%soilwater\|config%" | grep -v "!"
# Also check ASSOCIATE alias forms:
grep -rEn "=> state%soilwater" src/ --include="*.f90"
```
Both greps must return zero hits before deletion. Expect 20–25 hidden readers surfacing after the drop. Each compile-driven fixup is its own commit.

- [ ] **Drop dual-write legacy-global write side** from all home-tree files (soilhydraulics, waterbalance, boundtop, boundbottom, tillage, macropore, rootextraction, frozencond).
- [ ] **Compile** — `pixi run -e build` — fix compiler errors one by one. Each fix is its own commit labeled `refactor(state): SS-SWC compile-fix N — <symbol> in <file>`.
- [ ] **Comment out ~74 globals in `variables.f90`** with provenance markers:
  ```fortran
  ! real(8) theta(macp)  ! [SS-SWC] retired 2026-05-11 — moved to state%soilwater%theta (ADR 0038)
  ```
  Full list: theta, thetm1, thetar, thetas, thetsl, h, hm1, q, k, kmean, dimoca, cofgen, FrArMtrx, fluseksatexm, indeks, evp, pond, pondm1, pondini, gwl, gwlm1, hatm, volact, volm1, volini, wbalance, runon, fllowgwl, nodgwl, pegwl, bpegwl, npegwl, gwlflcpzo, nodgwlflcpzo, tra, iqredwet_day, iqreddry_day, iqredsol_day, iqredfrs_day, iptra_day, qpotrot_day, qredtot_day, inq, inqrot, inqssdi, IThetaBeg, IPondBeg, iqrot, iqssdi, iqredwet, iqreddry, iqredsol, iqredfrs, ies0, iet0, iew0, iintc, iruno, irunoCN, irunon, iqbot, iqtdo, iqtup, iqdo, iqup, iprec, igird, inird, cqssdi, cqrot, cqbot, cqbotdo, cqbotup, cinund, crunon, crunoff, crunoffCN, cqtdo, cqtup, cqprai, cgird, cnird.
- [ ] **Drop zero-init lines from `initialize.f90`** for all retired globals.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full` (integration gate — must be byte-identical, 5/5).

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWC Phase 2 S-2.12 — retire ~74 soil-water globals; drop dual-writes

Drop dual-writes from all home-tree + co-writer files. state%soilwater is
now authoritative for all soil-water-owned fields. Legacy globals commented
out in variables.f90 with [SS-SWC] provenance markers. Redundant zero-inits
dropped from initialize.f90. Compile-driven fixups resolved (N hidden readers).
pFUnit green. check-full 5/5 byte-identical.

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

### Task S-2.13: ADR 0038 + playbook update + merge prep

**Design decisions:** all
**Files:**
- Create: `docs/adr/0038-state-migration-soilwater-core.md`
- Modify: `docs/superpowers/specs/state-migration-playbook.md`
- Modify: `docs/adr/index.md` (add ADR 0038 entry)

- [ ] **Author ADR 0038** — status accepted; migration #8 (FINAL coupling-surface arc); `soilwater_state_t` extended with 30 flat-instantaneous + 8 per-day + `soilwater_intermediate_t` (22 fields, two methods) + `soilwater_cumulative_t` (14 fields); mild ADR 0033 extension (reset_per_day()); boundary D5/D6/D7 resolved; heat ADR 0034 hconduc sentinel resolved; mini-sim writeback fully retargeted; grid dims retained legacy; ~74 globals retired; compile-driven fixup count (actual). Cross-reference discovery, design, plan, ADRs 0030–0037.
- [ ] **Update playbook** — append soil-water-core lessons. Key patterns: (1) ADR 0033 mild extension (two methods on one cohort type); (2) multi-arc deferral resolution cascade (D5/D6/D7/heat ADR 0034 all resolved in one FINAL arc); (3) scale of compile-driven Phase 2.12 vs prior arcs (20–25 vs 13 for atmosphere).
- [ ] **Update `docs/adr/index.md`** with ADR 0038 entry.
- [ ] **Final PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
docs(adr): ADR 0038 — soil-water core state-type migration; playbook update

ADR 0038: FINAL coupling-surface arc; soilwater_state_t extended with ~74
new fields across 4 cadences; soilwater_intermediate_t (22 fields) + mild
ADR 0033 extension (reset_per_day()); soilwater_cumulative_t (14 fields).
Resolves boundary D5 (pond), D6 (gwl), D7 (kmean), heat ADR 0034 hconduc
sentinel (14 sites), swapoutput mini-sim writeback. Grid dims stay legacy.

Playbook: soil-water-core lessons appended (ADR 0033 dual-method extension,
deferral cascade resolution, compile-driven Phase 2.12 scale).

Spec: docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)
EOF
)"
```

---

## Self-Review Notes

- **S-1.4a/S-1.4b split rationale.** headcalc is ~840 LoC with the Richards solver iterating theta/h/k/q over all nodes and all timestep sub-steps. The split keeps each commit reviewable (~40 sites in S-1.4a vs ~50 in S-1.4b). If 1.4a proves smaller than expected, 1.4b content can be merged. If either is still too large, split further into S-1.4a (theta only), S-1.4b (h only), S-1.4c (q/k/kmean/dimoca).

- **S-1.5 precondition.** S-1.5 must complete before S-2.1. S-2.1 drops the inline legacy zeros in soilhydraulics; if integral's accumulator dual-writes are not installed first, the cumu cohort will only ever be reset (never written non-zero), and the pre-flight check in S-2.1 will catch this. Apply the atmosphere A-2.1 pre-flight check rigorously.

- **S-2.3 is the densest Phase 2 task.** soilhydraulics.f90 is ~1344 LoC with theta and h in nearly every compute statement. Consider splitting into S-2.3a (theta/h/q/k/kmean) and S-2.3b (dimoca/cofgen/FrArMtrx/hatm/pond/gwl/vol/wbalance) if the review burden is too high.

- **S-2.12 compile-driven scope.** 20–25 hidden readers is the planning estimate (discovery Section 3, phase-2 estimate); actual count from atmosphere was 13. soil-water has 23 external reader files vs atmosphere's 9, so expect higher. Budgeting 10–15 fixup commits means S-2.12 may take longer than any single prior task. Do not rush; each fix must pass check-full before proceeding.

- **Mini-sim writeback in S-2.11.** The `SoilWaterStateVar(task, state)` plumbing added in S-1.4a covers the task=1/task=2 save/restore for hm1/thetm1 etc. The swapoutput mini-sim (3865–3948) directly snapshots gwl/pond/theta/h — DIFFERENT from SoilWaterStateVar. Verify both paths; they are independent. The mini-sim calls SoilWaterStateVar at lines 3903 + 3915 (already plumbed in S-1.4a); the gwl/pond/theta/h snapshots at 3865–3873 and restores at 3941–3948 are inline variable copies, not SoilWaterStateVar calls.

- **Tillage init-order hazard (discovery hazard #8).** After S-1.9, both `DoTillage(1)` (cofgen write via Change_MvGpars) and `SoilWater(1)` (cofgen init-blank via `cofgen = 0.0d0` then rebuild) write through `state%soilwater%cofgen`. This is physically correct: DoTillage(1) sets the bulk-density target; SoilWater(1) then builds cofgen from PdmVG using that target. Verify via check-full after S-1.9. If regression detected, inspect the ordering of the `cofgen = 0.0d0` init-blank line in SoilWater(1) — if it follows DoTillage(1), the dual-write is consistent; if it precedes, there is a masking issue.

- **cgird/cnird in soil-water CUMU cohort.** irrigation.f90:93–95 continues to reset these fields during irrigation flzerocumu. After S-2.12 retires the global, irrigation.f90 must write through `state%soilwater%cumu%cgird` and `state%soilwater%cumu%cnird` at its reset sites — these will surface as compile-driven fixups in S-2.12. Add them to the expected hidden-reader list.
