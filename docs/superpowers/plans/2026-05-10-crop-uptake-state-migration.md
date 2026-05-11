# Crop Water Uptake State-Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Carve 22 crop-uptake fields (`qrot`, `qrosum`, stress components, JvL machinery) into `state%soilwater` (extending the boundary subset); resolve `qrot` ownership ambiguity; preserve byte-identical check-full at every commit.

**Architecture:** Flat extension of `soilwater_state_t`. Per-node arrays (first time in this state type) require `soilwater_init` signature bump to `(sw, numnod, nlay)`. `cropgrowth.f90` Task-1 init plumbed via `CropGrowth(1, state)` extension. `mfluxtable` init relocated to `soilwater_init`. No new cohort sub-records.

**Tech Stack:** Fortran 2008, ASSOCIATE `cw_` prefix (crop-water), pixi+meson, pFUnit, check-full regression.

**Spec:** `docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md`
**Discovery:** `docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-discovery.md`
**Branch:** all commits go on `refactor/crop-uptake-state` (already created from `development`).

---

## Lessons-learned applied

- **Dual-write before cutover.** Every global write gets a paired `state%soilwater%X = …` mirror before readers are cut over. Drop globals only in Phase 2 after all readers are migrated.
- **Verify-state-not-reports.** Run `check-full` output comparison, not just pFUnit, before every commit. pFUnit alone misses global-default regressions.
- **ASSOCIATE prefix `cw_`.** Shadow-safe aliasing in `rootextraction.f90` compute body. Follows the `bw_` (boundary), `ht_` (heat), `dr_` (drainage), `sl_` (solute) convention.
- **Grep field shapes during pre-flight.** Before deleting a global, grep both `use variables.*\b<var>\b` AND raw-symbol forms (both must be empty).
- **Compile-driven Phase 2.7 is expected.** After dropping dual-writes the compiler will surface 2–4 hidden readers. Plan for 2–4 fix-up commits.
- **Cumulatives out of scope.** `cqrot`, `iqrot`, `inqrot`, `iqredXXX*` — soil-water-core arc territory. Do not touch even if they appear adjacent to in-scope writes.

---

## File structure

**Modified:**
- `src/state/soilwater_state.f90` — extend type + bump `soilwater_init` signature (C-1.1, C-1.2)
- `src/core/swap.f90` — update `soilwater_init` call at line 184 (C-1.2)
- `src/crop/cropgrowth.f90` — plumb `CropGrowth(1, state)`; update 3 Task-1 init paths (C-1.3)
- `src/crop/rootextraction.f90` — dual-write 22 fields (C-1.4); ASSOCIATE `cw_` (C-1.4)
- `src/soil/soilhydraulics.f90` — `qrot(:)` reads cutover (C-2.1)
- `src/soil/waterbalance.f90` — `qrot/qrosum/qpotrot/qredXXX/qredXXXsum` reads cutover (C-2.2)
- `src/crop/oxygenstress.f90` — any reads of in-scope fields (C-2.3)
- `src/solute/solute.f90` — `qrot(:)` reads cutover (C-2.3)
- `src/solute/agetracer.f90` — `qrot(:)` reads cutover (C-2.3)
- `src/crop/cropgrowth.f90` — `flWrtNonox` compute reads cutover (C-2.4)
- `src/io/swapoutput.f90` — output cutover `.afo/.aun/.vap` qrot + `.rot` JvL suite (C-2.5)
- `src/io/swap_csv_output.f90` — any direct reads cutover if compile-surfaced (C-2.5)
- `src/core/variables.f90` — retire 22 globals (C-2.6)
- `src/core/initialize.f90` — drop redundant zero-init lines for retired globals (C-2.6)
- `tests/unit/state/test_soilwater_state.pf` — extend with per-node alloc + scalar default tests (C-1.1)
- `docs/adr/0036-state-migration-crop-uptake.md` — create at end (C-2.6)

---

# Phase 0 — Regression coverage audit

### Task C-0.1: Audit `swdrought=2` / `swoxygen=2` / `swcompensate` / `swfrost` coverage

**Design decisions:** D11
**Files:** read-only audit (no code changes expected)

The JvL path (`rootextraction.f90` lines 313–760) is large. Zero integration coverage would be a notable risk. Check the five TOML regression cases.

- [ ] **Grep TOML fixtures for switch values:**
  ```bash
  grep -rn "swdrought\|swoxygen\|swcompensate\|swfrost" tests/ --include="*.toml" | sort
  ```
- [ ] **Note which values appear** — specifically `swdrought = 2`, `swoxygen = 2`, `swcompensate = 1` or `2`, `swfrost = 1`.
- [ ] **Document findings** — if any critical path is uncovered, add a comment to the ADR 0036 stub's known-issues section. Do not add new fixtures in this arc.
- [ ] **No commit needed** unless a trivially fixable gap is found. If a documentation artifact is warranted:
```bash
git commit -m "$(cat <<'EOF'
docs(test): SS-CRP Phase 0 C-0.1 — regression coverage audit (JvL + stress paths)

Audit of five TOML regression cases for swdrought=2 (JvL), swoxygen=2,
swcompensate, swfrost coverage. Findings documented as known-issues in
ADR 0036 stub. No new fixtures added this arc.

Spec: docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md
EOF
)"
```

---

# Phase 1 — State type extension, init plumbing, dual-write

### Task C-1.1: Extend `soilwater_state_t` with 22 new fields; add pFUnit tests

**Design decisions:** D2, D3, D7, D8, D10
**Files:**
- Modify: `src/state/soilwater_state.f90`
- Modify: `tests/unit/state/test_soilwater_state.pf`

Append the 22 crop-uptake fields to `soilwater_state_t` (after the existing 12 boundary fields) per the D2 layout. All 12 per-node arrays declared `allocatable` with no default initializers (allocation deferred to C-1.2). All scalars get explicit zero/false defaults. DO NOT change `soilwater_init` signature yet — that is C-1.2.

- [ ] **Read `src/state/soilwater_state.f90`** — note the closing `end type soilwater_state_t` line; insert the 22 new field declarations immediately before it.
- [ ] **Add the 22 fields** per D2 target layout: 6 primary per-node allocatables (`qrot`, `qpotrot`, `qredwet`, `qreddry`, `qredsol`, `qredfrs`), 5 primary scalars (`qrosum`, `qredwetsum`, `qreddrysum`, `qredsolsum`, `qredfrssum`), 1 flag (`flWrtNonox`), 6 JvL per-node allocatables (`mflux`, `mroot`, `hroot`, `rootrho`, `rootphi`, `rmax`), 1 layer×801 allocatable (`mfluxtable`), 4 JvL scalars (`Tactual`, `alpJvLier`, `hleaf`, `Hxylem`).
- [ ] **Write failing pFUnit tests** in `test_soilwater_state.pf`:
  - `@test` for each new scalar: assert default value = `0.0_real64` (logicals = `.false.`) on a fresh `soilwater_state_t()` literal.
  - `@test` that each allocatable is `.not. allocated(sw%qrot)` etc. before init.
  - `@test` independent-instance isolation: two `soilwater_state_t` variables share no state.
  - Confirm FAIL (allocatables unallocated until C-1.2).
- [ ] **Verify compile** — `pixi run -e test test-pfunit` (expect scalar tests pass; per-node tests fail until C-1.2).
- [ ] **Verify check-full byte-identical** — `pixi run check-full` (no code paths changed yet).
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
feat(state): SS-CRP Phase 1 C-1.1 — extend soilwater_state_t with 22 crop-uptake fields

Adds 22 owned fields (6 primary per-node allocatables, 5 primary scalars,
1 flag, 6 JvL per-node allocatables, 1 mfluxtable layer×801 allocatable,
4 JvL scalars) to soilwater_state_t. All scalars zero/false-defaulted;
allocatables unallocated until soilwater_init (C-1.2). pFUnit scalar
defaults green; per-node alloc tests pending init. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md (ADR 0036)
EOF
)"
```

---

### Task C-1.2: Bump `soilwater_init(sw, numnod, nlay)` + allocate arrays; update `swap.f90:184`

**Design decisions:** D4, D6
**Files:**
- Modify: `src/state/soilwater_state.f90`
- Modify: `src/core/swap.f90`

This is the most consequential task. It changes `soilwater_init`'s public signature, propagating to `swap.f90`. Also relocates the `mfluxtable` build from `cropgrowth.f90` Task-1 calls into `soilwater_init` (D6) — eliminating one co-write before the dual-write phase.

**Signature:** `subroutine soilwater_init(sw, numnod, nlay)`

`numnod` and `numlay` sources: `swap.f90:71` already imports both from `variables` via `use variables, only: …, numnod, numlay`. Both are set by `CalcGrid()` at line 183, immediately before the `soilwater_init` call at line 184.

`mfluxtable` relocation: `soilwater_init` calls the `MatricFlux(1, ...)` body equivalent (or delegates to `MatricFlux(1, h(1), 1, dummy)` with the state arg) to populate `sw%mfluxtable` when `swdrought == 2`. The existing `MatricFlux(1)` body in `rootextraction.f90` writes to legacy `mfluxtable(maho, 801)` — after relocation it writes to `sw%mfluxtable(nlay, 801)`. Verify that `h(1)`, `swdrought`, and the layer-config arrays (`iHWCKmodel`, `botcom`, etc.) required by `MatricFlux(1)` are available in the `swap.f90` context at line 184 (all are legacy globals set by `CalcGrid`). If they are not available at that call site, raise this as a blocker and consider a lazy-init fallback.

- [ ] **Read `soilwater_init` current body** — confirm the 12 scalar zero assignments; plan insertion point for allocations.
- [ ] **Update `soilwater_init` signature** to `subroutine soilwater_init(sw, numnod, nlay)`.
- [ ] **Add allocate statements** for all 12 per-node arrays (`qrot`, `qpotrot`, `qredwet`, `qreddry`, `qredsol`, `qredfrs`, `mflux`, `mroot`, `hroot`, `rootrho`, `rootphi`, `rmax`) to size `numnod`. Zero-init all after allocation.
- [ ] **Add allocate for `mfluxtable(nlay, 801)`** and invoke the lookup build (or call `MatricFlux(1, ...)` if that is cleaner). If `MatricFlux` cannot be called cleanly at this site, defer `mfluxtable` allocation to C-1.3 and document the decision.
- [ ] **Update `swap.f90:184`** — change `call soilwater_init(state%soilwater)` to `call soilwater_init(state%soilwater, numnod, numlay)`.
- [ ] **Complete the C-1.1 pFUnit tests** — the per-node alloc `@test`s that were pending should now pass.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full` (byte-identical; no compute paths touch new arrays yet).
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
feat(state): SS-CRP Phase 1 C-1.2 — soilwater_init(sw, numnod, nlay); per-node arrays allocated

soilwater_init gains numnod and nlay args; allocates and zeros all 12 per-node
arrays (qrot, qpotrot, qredwet/dry/sol/frs, mflux, mroot, hroot, rootrho,
rootphi, rmax) and mfluxtable(nlay,801). swap.f90:184 caller updated to pass
numnod/numlay. mfluxtable init relocated from cropgrowth Task-1 (D6).
pFUnit all-green. check-full 5/5 byte-identical.

Spec: docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md (ADR 0036)
EOF
)"
```

---

### Task C-1.3: Plumb `CropGrowth(1, state)` — update 3 Task-1 init paths in `cropgrowth.f90`

**Design decisions:** D5, D6
**Files:**
- Modify: `src/crop/cropgrowth.f90`
- Modify: `src/core/swap.f90` (Task-1 call site)

Current call: `swap.f90:272` — `call CropGrowth(1, state%heat%tsoil)`. Target: `call CropGrowth(1, state)`. The three Task-1 init paths in `cropgrowth.f90` (CropFixed: lines 595–612, CropWofost: 1325–1342, CropGrass: 2394–2411) each contain:

```fortran
if (swdrought .eq. 2) then
   call MatricFlux(1, h(1), 1, dummy)   ! build mfluxtable — NOW RELOCATED to soilwater_init (C-1.2)
   if (swhydrlift .eq. 1) then
      flhydrlift = .true.
   else
      flhydrlift = .false.
   endif
   do i = 1,numnod
      twilt(i) = watcon(i,wiltpoint)
      hroot(i) = h(i)                   ! seed → state%soilwater%hroot(i)
   enddo
   hleaf = -2000.d0                     ! seed → state%soilwater%hleaf
endif
```

After C-1.2, `MatricFlux(1)` calls at lines 601, 1331, 2400 are dead (mfluxtable already built in `soilwater_init`). This task: (a) removes those calls, (b) replaces the `hroot(i) = h(i)` and `hleaf = -2000.d0` writes with `state%soilwater%hroot(i) = h(i)` and `state%soilwater%hleaf = -2000.d0_real64`, (c) updates the `CropGrowth` signature to accept full `state`.

**Verify before implementing:** `grep -n "tsoil\|state" src/crop/cropgrowth.f90 | head -20` — confirm how `state%heat%tsoil` is currently passed in and that Tasks 2/3 already have `state` in scope from heat-arc Task 6 (per discovery Section 1 signature table).

- [ ] **Confirm `CropGrowth` signature and Task-2/3 state scope** — grep as above. Confirm Tasks 2/3 already receive `state` (or `state%heat%tsoil`) so extending to full `state` is incremental.
- [ ] **Update `CropGrowth` Task-1 dummy arg** — change from `tsoil(:)` slice (or equivalent) to `type(swap_state_t), intent(inout) :: state`. Add `use swap_state_mod` in `cropgrowth.f90` if not already present.
- [ ] **Remove `MatricFlux(1)` calls** at lines 601, 1331, 2400 (mfluxtable already built by C-1.2).
- [ ] **Replace `hroot(i) = h(i)` seeds** with `state%soilwater%hroot(i) = h(i)` at all three locations.
- [ ] **Replace `hleaf = -2000.d0` seeds** with `state%soilwater%hleaf = -2000.0_real64` at all three locations. Keep the legacy `hleaf = -2000.d0` global write as dual-write (drop in C-2.6).
- [ ] **Update `swap.f90:272`** — `call CropGrowth(1, state)`.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-CRP Phase 1 C-1.3 — CropGrowth(1, state) plumbing; hroot/hleaf seed to state

CropGrowth Task-1 call extended from tsoil slice to full state. Three
Task-1 init paths (CropFixed/Wofost/Grass) now write state%soilwater%hroot
and state%soilwater%hleaf alongside legacy globals (dual-write). MatricFlux(1)
calls removed from cropgrowth Task-1 (mfluxtable init relocated to
soilwater_init per D6). swap.f90:272 updated. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md (ADR 0036)
EOF
)"
```

---

### Task C-1.4: `rootextraction.f90` dual-write — all 22 owned fields

**Design decisions:** D1, D4, D7, D8, D9
**Files:**
- Modify: `src/crop/rootextraction.f90`

Install dual-writes for all 22 owned fields across `RootExtraction`, `JongvanLier`, `JongvanLierLoop`, and `MatricFlux`. Use ASSOCIATE with `cw_` prefix at the top of each subroutine body to alias `state%soilwater` fields. `RootExtraction` is medium-sized (lines 31–303); use explicit ASSOCIATE for the dense compute inside the Feddes path and the JvL path separately. The discovery identified all write sites (Section 2a); use those anchors.

Key write-site clusters:
- `RootExtraction`: `qrot` (lines 59, 88, 198, 199, 278), `qrosum` (61, 108, 200, 282), `qpotrot` (198), `qredwet/dry/sol/frs` (208–211, 219–222), `qredwetsum/…sum` (225–228, 286–295), `Tactual` (60), `flWrtNonox` (154, 156).
- `JongvanLier`: `mflux` (344, 383, 389), `mroot` (345), `hroot` (346, 353, 354, 411, 678, 709), `rootrho` (347, 362, 374), `rootphi` (348, 365, 377), `rmax` (362, 373), `alpJvLier` (628, 631), `Hxylem` (405, 419, 446, 480, 518, 531, 559, 578, 610), `hleaf` (339, 401, 416, 418, 445, 479, 516, 525), `qrosum` (616, 618, 627, 666).
- `JongvanLierLoop`: `hroot` (709, 729, 751), `mroot` (711, 731, 753), `qrot` (714, 720, 725, 728, 736, 742, 747, 750), `qrosum` (756).
- `MatricFlux(2)`: reads `state%solute%cml` (already state-plumbed). `mfluxtable` read at 836, 841, 856, 861 — after C-1.2 these can read `state%soilwater%mfluxtable` (or keep reading legacy global until C-2.6 if the dual-write approach is simpler).

Pattern recommendation: use explicit pair-writes for small routines; for `RootExtraction`'s dense body use ASSOCIATE block:

```fortran
associate( &
   cw_qrot     => state%soilwater%qrot, &
   cw_qrosum   => state%soilwater%qrosum, &
   cw_flWrtNon => state%soilwater%flWrtNonox &
   ! ... etc.
)
   ! body uses both bare names (legacy) and cw_ names (state)
end associate
```

- [ ] **List all write sites** using discovery Section 2a table as the authoritative source. Do not rely solely on grep — use line anchors from the table.
- [ ] **Add ASSOCIATE blocks** at top of `RootExtraction`, `JongvanLier`, `JongvanLierLoop`. For `MatricFlux(2)` which is smaller, use direct `state%soilwater%` writes.
- [ ] **Add dual-write at every write site** — legacy global write followed immediately by `cw_X = X` (or `state%soilwater%X = X`). Confirm `Tactual = qrosum` ordering: the read of prior `qrosum` at line 60 must precede the reset at line 61; preserve this order (`state%soilwater%Tactual = state%soilwater%qrosum` then `qrosum = 0` / `state%soilwater%qrosum = 0`).
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full` (byte-identical; dual-writes are additive).
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-CRP Phase 1 C-1.4 — dual-write 22 crop-uptake fields in rootextraction.f90

RootExtraction, JongvanLier, JongvanLierLoop, MatricFlux all write both
legacy globals and state%soilwater%X for all 22 owned fields. ASSOCIATE
cw_ prefix in primary routines. Legacy globals unchanged — check-full 5/5
byte-identical. State fields now mirror legacy values at every timestep.

Spec: docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md (ADR 0036)
EOF
)"
```

---

# Phase 2 — Reader cutover and global retirement

### Task C-2.1: `soilhydraulics.f90` reader cutover — 10 `qrot(:)` read sites

**Design decisions:** D1, D10
**Files:**
- Modify: `src/soil/soilhydraulics.f90`

10 `qrot(:)` read sites in `headcalc` (Richards F-vector) and `SoilWater(3)` integrand (lines 120, 202, 222, 244, 248, 481, 507, 527, 531, 719). `state` is already in scope throughout `soilhydraulics.f90`. Pure reader of `qrot` — no co-write, no cumulative touches. Drop `qrot` from `use variables, only:` after all 10 sites are migrated.

Note: do NOT touch `inqrot`, `iqrot`, `iqredXXX`, `cqrot` reads in this file — those are downstream cumulatives, out of scope (D10).

- [ ] **List all 10 read sites** — `grep -n "qrot" src/soil/soilhydraulics.f90 | grep -v "state%\|!" | head -20`. Cross-check with discovery Section 3.1.
- [ ] **Migrate reads** — replace `qrot(node)` / `qrot(i)` etc. with `state%soilwater%qrot(node)` etc. Use ASSOCIATE or direct access.
- [ ] **Drop `qrot` from `use variables, only:`**.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-CRP Phase 2 C-2.1 — soilhydraulics reads state%soilwater%qrot (10 sites)

All 10 qrot(:) read sites in headcalc/SoilWater(3) now read
state%soilwater%qrot. Dropped from use variables. Downstream cumulatives
(inqrot/iqrot/cqrot) untouched — soil-water-core territory (D10).
State authoritative for qrot reads; dual-write still active.

Spec: docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md (ADR 0036)
EOF
)"
```

---

### Task C-2.2: `waterbalance.f90` reader cutover — `qrot/qrosum/qpotrot/qredXXX/qredXXXsum` (~12 sites)

**Design decisions:** D1, D10
**Files:**
- Modify: `src/soil/waterbalance.f90`

Reads to migrate: `qrot` (lines 348, 420), `qrosum` (lines 339, 407, 415, 418, 490), `qpotrot` (line 421), `qredwet/dry/sol/frs` (line 422), `qredwetsum/…sum×4` (lines 428–435). `state` should already be in scope in `integral` and `checkmassbal`; verify before starting.

**Important constraint:** the 10-line accumulator block at lines 418–435 reads this arc's instantaneous scalars and writes into `iqrot`/`iqredXXX*` (downstream cumulatives). Migrate only the READ side (from `state%soilwater%X`); the WRITE side (into `iqrot` etc.) stays on legacy globals until soil-water-core arc (D10).

`qrosum` at line 339 enters the mass balance: `qbot = qtop + qrosum + …`. After this task: `state%soilwater%qbot = state%soilwater%qtop + state%soilwater%qrosum + …`. Verify existing line in context to avoid accidental double-state-field resolution.

- [ ] **Confirm state scope** in `integral` and `checkmassbal` — `grep -n "subroutine integral\|subroutine checkmassbal\|intent.*state" src/soil/waterbalance.f90 | head`. If `state` is not a dummy arg, trace the call chain and add it; update callers from `swap.f90`.
- [ ] **Migrate reads** of `qrot, qrosum, qpotrot, qredwet, qreddry, qredsol, qredfrs, qredwetsum, qreddrysum, qredsolsum, qredfrssum` to `state%soilwater%X`. Keep the write side of the accumulator block on legacy globals.
- [ ] **Drop migrated symbols from `use variables, only:`**.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-CRP Phase 2 C-2.2 — waterbalance reads state%soilwater (qrot + stress fields)

integral() and checkmassbal() in waterbalance.f90 read qrot/qrosum/qpotrot/
qredXXX/qredXXXsum from state%soilwater (~12 sites). Accumulator write side
(iqrot/iqredXXX) unchanged — soil-water-core territory (D10). qrosum in
mass-balance line confirmed consistent with state%soilwater%qtop/qbot.

Spec: docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md (ADR 0036)
EOF
)"
```

---

### Task C-2.3: `oxygenstress.f90`, `solute.f90`, `agetracer.f90` reader cutover — `qrot(:)` (6 sites)

**Design decisions:** D1
**Files:**
- Modify: `src/crop/oxygenstress.f90`
- Modify: `src/solute/solute.f90`
- Modify: `src/solute/agetracer.f90`

`solute.f90` lines 223–225 and `agetracer.f90` lines 224, 225, 265 — passive uptake and age-tracer uptake reads of `qrot(:)`. Both already take `state` from the solute migration arc.

`oxygenstress.f90` — verify whether it reads any of the 22 in-scope fields directly. Discovery Section 3 did not list it as a reader, but `OxygenStress` takes `state` as `optional intent(in)` from SS-HEAT Task 6. Quick grep to confirm before skipping.

- [ ] **Grep `oxygenstress.f90`** — `grep -n "qrot\|qrosum\|qredwet\|qreddry\|hroot\|hleaf\|mflux\|flWrtNonox" src/crop/oxygenstress.f90 | grep -v "!"`. If no hits, skip this file.
- [ ] **Migrate `qrot` reads in `solute.f90`** (lines 223–225) — `state%soilwater%qrot(node)` etc. Drop from `use variables, only:`.
- [ ] **Migrate `qrot` reads in `agetracer.f90`** (lines 224, 225, 265) — same pattern.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-CRP Phase 2 C-2.3 — solute/agetracer/oxygenstress reads state%soilwater

solute.f90 and agetracer.f90 read qrot(:) from state%soilwater (6 sites).
oxygenstress.f90 confirmed clean (no direct reads of in-scope fields).
Dropped from use variables in both files.

Spec: docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md (ADR 0036)
EOF
)"
```

---

### Task C-2.4: `cropgrowth.f90` `flWrtNonox` compute reads cutover (5 sites)

**Design decisions:** D1, D8
**Files:**
- Modify: `src/crop/cropgrowth.f90`

`flWrtNonox` is read at 5 sites in cropgrowth (lines 689, 1742, 2104, 3001, 3387) gating `rr = 0` / `grrt = 0` in WOFOST and grass dynamics. After C-1.3 cropgrowth Tasks 2/3 already have `state` in scope. Replace bare `flWrtNonox` reads with `state%soilwater%flWrtNonox`.

**Verify first:** confirm that Tasks 2 and 3 in `cropgrowth.f90` already receive full `state` (not just `state%heat%tsoil`) as a result of C-1.3. If Task-1's signature change did not propagate to Tasks 2/3, they may still use a `tsoil` slice — trace the `CropGrowth(task, ...)` dispatch and verify state is accessible at all 5 read sites.

- [ ] **Confirm state scope at all 5 `flWrtNonox` read sites** — `grep -n "flWrtNonox\|state" src/crop/cropgrowth.f90 | grep -v "!" | head -30`. Verify `state` (not just a slice) is in scope at each caller.
- [ ] **Migrate all 5 reads** — `state%soilwater%flWrtNonox` for each bare-name read. Drop `flWrtNonox` from `use variables, only:` clauses at those read sites.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-CRP Phase 2 C-2.4 — cropgrowth reads state%soilwater%flWrtNonox (5 sites)

CropFixed/Wofost/Grass dynamics (rr=0/grrt=0 gates) now read flWrtNonox
from state%soilwater (5 sites in cropgrowth.f90). Dropped from use variables.
flWrtNonox documented as soil-water-owned with future crop-state relocation
note in ADR 0036 (D8).

Spec: docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md (ADR 0036)
EOF
)"
```

---

### Task C-2.5: `swapoutput.f90` + `swap_csv_output.f90` output cutover (~14 sites)

**Design decisions:** D1
**Files:**
- Modify: `src/io/swapoutput.f90`
- Modify: `src/io/swap_csv_output.f90` (if compile-surfaced direct reads)

Output reads to migrate in `swapoutput.f90`:
- `.afo/.aun/.vap` qrot column output (~4 sites: lines 581, 630, 651, 690, 709)
- `.rot` rotation output — full JvL array suite: `hleaf`, `Hxylem`, `hroot`, `mflux`, `mroot`, `rootrho`, `rootphi` (~9 sites: lines 774–875 block)
- `.bal` output — read `qrosum` if present (confirm with grep)

`swap_csv_output.f90` reads only downstream cumulatives (`iqrot`, `iqredXXX`) — those are out of scope. Verify with grep; migrate only if direct reads of in-scope fields are found.

There is no mini-sim writeback touching qrot or the JvL arrays (boundary lesson #6 cross-check confirmed: `swapoutput.f90:3745–3820` snapshots `qbot`/`gwl`/`pond` only). No writeback retarget needed for this arc.

- [ ] **List all output read sites** — `grep -n "qrot\|qrosum\|qpotrot\|qredwet\|qreddry\|qredsol\|qredfrs\|hleaf\|Hxylem\|hroot\|mflux\|mroot\|rootrho\|rootphi\|rmax\|flWrtNonox\|Tactual\|alpJvLier" src/io/swapoutput.f90 src/io/swap_csv_output.f90 | grep -v "state%\|!" | head -40`. Note which fields appear.
- [ ] **Migrate output reads** to `state%soilwater%X` in `swapoutput.f90`.
- [ ] **Grep `swap_csv_output.f90`** for any direct reads — if present, migrate; if only downstream cumulatives, confirm and skip.
- [ ] **Drop migrated symbols from `use variables, only:`**.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-CRP Phase 2 C-2.5 — output reads state%soilwater (qrot + JvL suite)

swapoutput.f90 reads qrot(.afo/.aun/.vap columns) and JvL arrays
(hleaf/Hxylem/hroot/mflux/mroot/rootrho/rootphi in .rot output) from
state%soilwater. swap_csv_output.f90 confirmed indirect-only (no direct
in-scope field reads). No mini-sim writeback needed (qrot not snapshotted).

Spec: docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md (ADR 0036)
EOF
)"
```

---

### Task C-2.6: Drop dual-writes + retire 22 legacy globals; ADR 0036

**Design decisions:** D1, D8, D10, D12
**Files:**
- Modify: `src/crop/rootextraction.f90`
- Modify: `src/crop/cropgrowth.f90` (drop hroot/hleaf legacy dual-write side)
- Modify: `src/core/variables.f90`
- Modify: `src/core/initialize.f90`
- Create: `docs/adr/0036-state-migration-crop-uptake.md`

**This task has the highest compile-driven discovery risk.** After dropping dual-writes, expect 2–4 hidden readers to surface. Plan for iterative fix-up commits before finalizing the retirement commit. Do not bundle fix-ups into this commit — each one is a separate `refactor(state): SS-CRP compile-fix …` commit.

Pre-flight: for each of the 22 globals, run the standard two-grep test before any deletion:
```bash
# Run for each of: qrot qrosum qpotrot qredwet qreddry qredsol qredfrs
#   qredwetsum qreddrysum qredsolsum qredfrssum flWrtNonox
#   mflux mroot hroot rootrho rootphi rmax mfluxtable
#   Tactual alpJvLier hleaf Hxylem
grep -rEn "use variables.*\b${X}\b" src/ --include="*.f90" | grep -v "variables.f90\|initialize.f90"
grep -rEn "\b${X}\b" src/ --include="*.f90" | grep -v "variables.f90\|initialize.f90\|state%soilwater\|config%" | grep -v "!"
```
Both must return zero hits. Investigate any remaining hit before deletion.

**After pre-flight passes:**
- [ ] **Drop dual-write legacy-global write side** in `rootextraction.f90` (all 22 fields) and `cropgrowth.f90` (hroot, hleaf).
- [ ] **Compile** — `pixi run -e build` — fix any compiler errors one by one; each fix is its own commit.
- [ ] **Comment out 22 globals in `variables.f90`** with provenance markers:
  ```fortran
  ! real(8) qrot(macp)  ! [SS-CRP] retired 2026-05-11 — moved to state%soilwater%qrot (ADR 0036)
  ```
- [ ] **Drop zero-init lines from `initialize.f90`** at lines 353, 368, 401, 420–425, 804 (and nearby) that zero `qrot`, `qrosum`, `hroot`, `mflux`, `alphacrit`, etc. — redundant after state defaults.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full` (integration gate — must be byte-identical).
- [ ] **Author ADR 0036** — create `docs/adr/0036-state-migration-crop-uptake.md` with status accepted, migration #6, key decisions (D1–D12), consequences (22 globals retired, 1 signature change, ~47 reads migrated, compile-driven discoveries), known residuals (`flWrtNonox` semantic note, cumulatives deferred, swdrought=2 coverage gap if found in C-0.1). Reference discovery, design, plan.
- [ ] **Update playbook** — append crop-uptake lessons to `docs/superpowers/specs/state-migration-playbook.md` if new patterns emerged (mfluxtable relocation, per-node array allocation in state-init, CropGrowth slice-to-full-state progression).
- [ ] **Final commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-CRP Phase 2 C-2.6 — retire 22 crop-uptake globals; ADR 0036

Drop dual-writes; state%soilwater is now authoritative for all 22 crop-water-
uptake fields. Legacy globals commented out in variables.f90 with [SS-CRP]
provenance markers. Redundant zero-inits dropped from initialize.f90.
Compile-driven fixups resolved (N hidden readers surfaced and migrated).

pFUnit green. check-full 5/5 byte-identical.
ADR: docs/adr/0036-state-migration-crop-uptake.md

Spec: docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md
EOF
)"
```

---

## Self-Review Notes

- **C-1.2 mfluxtable relocation risk:** The `MatricFlux(1)` body reads layer-config globals (`iHWCKmodel`, `botcom`, `paramvg`, etc.) that are populated by CalcGrid or the soil-physics init chain. Verify these are all set by the time `soilwater_init` is called at `swap.f90:184`. If any required global is not yet set at that point, the relocation must be deferred — either to a later init call or to a lazy-init inside `RootExtraction` on first call. This is the key blocker risk for C-1.2; resolve by inspection during implementation.
- **C-1.3 CropGrowth signature propagation:** The current Task-1 signature takes `state%heat%tsoil` (a slice), not full `state`. Extending to full `state` requires that all `CropGrowth(task, ...)` dispatch paths accept `swap_state_t`. Confirm whether Tasks 2/3 already take full state or also use slices — if slices, the extension must cover all tasks simultaneously to avoid an inconsistent interface. Discovery Section 7 / Hazard #2 flagged this; the heat-arc `dummy_X_ => tsoil` rename trick may have left non-module routines out of the state-plumbing. Grep before implementing.
- **C-2.2 state scope in `waterbalance.f90`:** `integral()` and `checkmassbal()` may not yet have `state` as a dummy arg. If plumbing is needed, the callers in `swap.f90` must be updated — add ~30 min to C-2.2 estimate if so. Confirmed pattern from boundary B-2.2 self-review notes.
- **C-2.4 state scope at cropgrowth `flWrtNonox` sites:** If Tasks 2/3 in cropgrowth still use a `tsoil` slice rather than full state (see C-1.3 note above), the `flWrtNonox` cutover in C-2.4 may require resolving the slice-to-full-state issue first. Sequence: finish C-1.3 fully before starting C-2.4.
- **C-2.6 compile-driven discoveries:** boundary arc surfaced 4; crop-uptake's reader surface is smaller (47 vs 80 sites) but more spread across subsystem boundaries. Expect 2–4. The most likely surprise locations based on discovery Section 7 / Hazard #7: `swapoutput.f90` use clauses (JvL array reads), `waterbalance.f90` use clauses for `qpotrot`/`qredXXX`, and possibly `agetracer.f90` or `solute.f90` extra reads of auxiliary fields not captured in the per-field reader index.
- **D6 mfluxtable allocation decision:** If the C-1.2 risk materializes (layer-config not ready at `soilwater_init` call time), the fallback is: allocate `mfluxtable(nlay, 801)` in `soilwater_init` but keep the `MatricFlux(1)` calls in `cropgrowth.f90` Task-1 (they write into `state%soilwater%mfluxtable` instead of the legacy global). This is slightly less clean than full relocation but still correct. Document the deviation in ADR 0036 if taken.
