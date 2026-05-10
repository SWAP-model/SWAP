# Heat State Migration — Implementation Plan (all phases)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Migrate the heat (soil-temperature) subsystem off `variables.f90` globals into a typed `heat_state_t` aggregated under `swap_state_t`, plumb `Temperature(task)` with state, fix the `outheapar` output-writeback hazard, fix the `rfcp` co-write in soilhydraulics, and patch the `swcalt=1` config gap by promoting 6 missing fields.

**Architecture:** Three logical phases bundled into one plan (heat is small):
- **Phase 0** patches `heat_config_t` for 6 missing fields supporting `swcalt=1` (analytical method).
- **Phase 1** defines `heat_state_t` (flat — no cohorts; all 13 owned globals are instantaneous), threads `state` through `Temperature(task)`, dual-write transitional pattern, output reads switch, fixes `outheapar` writeback.
- **Phase 2** migrates cross-subsystem `tsoil` readers, fixes `rfcp` co-write in soilhydraulics, drops dual-write, removes globals.

**Tech Stack:** Fortran 2008, gfortran, meson + ninja + pixi, pFUnit, check-full byte-identical regression.

**Spec:** `docs/superpowers/specs/2026-05-10-state-migration-heat-design.md`
**Discovery:** `docs/superpowers/specs/2026-05-10-state-migration-heat-discovery.md`

**Branch:** all commits go on `refactor/heat-state` (new from development).

---

## Lessons-learned applied (from surfacewater + drainage + solute + cumulative-reset cohorts)

- **Section 2 categorization upfront.** Discovery already populated. Heat's owned globals are all instantaneous → flat `heat_state_t`, no cohort sub-records.
- **Activity-gate awareness.** No cumulative fields, so the Phase A correction lesson (partition cohorts by activity gate) doesn't apply here. But `flTemperature` gates the entire subsystem — note for design doc.
- **Verification grep BEFORE deletion.** Both `use variables.*\b<var>\b` AND raw-symbol AND alias-form (`<alias> => state%heat`) greps before global removal.
- **`ht_*` ASSOCIATE prefix** for shadow-safe aliasing in compute bodies that also import `use variables`.
- **Optional state args** for routines called from many call chains (the MatricFlux pattern from solute Phase 2).

---

## File Structure

**Created:**
- `src/state/heat_state.f90` — `heat_state_t` definition (Phase 1)
- `tests/unit/state/test_heat_state.pf` — pFUnit tests (Phase 1)
- `tests/unit/config/test_heat_phase0_fields.pf` — pFUnit tests for Phase 0 fields
- `docs/adr/0034-state-migration-heat.md` — finalized at end (Phase 2)

**Modified:**
- `src/config/heat_config.f90` — add 6 fields (Phase 0)
- `src/io/toml/read_heat_toml.f90` — read 6 fields (Phase 0)
- `src/io/toml/config_to_variables.f90` — `apply_heat` populates 6 legacy globals (Phase 0)
- `src/state/swap_state.f90` — add `heat` field (Phase 1)
- `src/heat/temperature.f90` — `Temperature(task)` gains `state` arg; ASSOCIATE in body; dual-write; drop dual-write (Phase 1, Phase 2)
- `src/heat/frozencond.f90` — already takes state; migrate `state%heat%*` reads/writes (Phase 1, Phase 2)
- `src/io/swapoutput.f90` — `outheapar` writeback fix (move devries to compute); `outtem`, `TemperatureOutput` switch reads to state (Phase 1)
- `src/io/swap_csv_output.f90` — `set_values` reads from state if applicable (Phase 1)
- `src/soil/soilhydraulics.f90` — `rfcp = 1.0d0` reset migrates to `state%heat%rfcp` (Phase 2)
- 9 other compute-reader files for `tsoil` (Phase 2): `crop/rootextraction.f90`, `crop/oxygenstress.f90`, `crop/management_soil.f90`, `crop/cropgrowth.f90`, `solute/solute.f90`, `atmosphere/snow.f90`, `utils/soilhydraulicsutils.f90`, `boundary/boundtop.f90`, `boundary/boundbottom.f90`
- `src/core/swap.f90` — pass state to `Temperature(task)` calls; pass to `outheapar` if needed (Phase 1)
- `src/core/variables.f90` — comment out 13 declarations with provenance (Phase 2)
- `src/core/initialize.f90` — drop matching zero-init lines (Phase 2)
- `tests/unit/testSuites.inc` — register new suites
- meson source list — add `src/state/heat_state.f90`

---

# Phase 0 — Config gap (1 task)

### Task 1: Promote 6 missing config fields for `swcalt=1`

Per discovery hazard #4: `swcalt=1` (analytical method) lacks 4 scalars (`ddamp`, `tmean`, `tampli`, `timref`) and 2 boundary-table arrays (`temtoptab(mabbc*2)`, `tembtab(mabbc*2)`) in `heat_config_t`. Same shape as solute Phase 0.

**Files:**
- Modify: `src/config/heat_config.f90` (add 6 fields + validators)
- Modify: `src/io/toml/read_heat_toml.f90` (read 6 fields)
- Modify: `src/io/toml/config_to_variables.f90` (`apply_heat` populates legacy globals)
- Create: `tests/unit/config/test_heat_phase0_fields.pf` (validator tests)
- Modify: `tests/unit/testSuites.inc`, `tests/unit/meson.build`

- [ ] **Step 1: Inspect existing heat_config_t**

```bash
sed -n '1,80p' src/config/heat_config.f90
grep -n "subroutine\|public\|type ::" src/config/heat_config.f90 | head
```

Note the existing field declarations and validator patterns.

- [ ] **Step 2: Inspect existing TOML reader**

```bash
sed -n '1,40p' src/io/toml/read_heat_toml.f90
```

- [ ] **Step 3: Inspect adapter**

```bash
grep -n "apply_heat\|swhea\b\|swcalt" src/io/toml/config_to_variables.f90 | head
```

- [ ] **Step 4: Write failing pFUnit tests**

Create `tests/unit/config/test_heat_phase0_fields.pf`. Mirror the pattern of `tests/unit/config/test_solute_phase0_fields.pf` (was created during solute Phase 0). For each scalar: default value test, in-range, out-of-range. For tables: default unallocated, allocation lifecycle.

```fortran
@test
subroutine test_heat_config_ddamp_default_zero()
   use funit
   use iso_fortran_env, only: real64
   use heat_config_mod, only: heat_config_t
   type(heat_config_t) :: cfg
   @assertEqual(0.0_real64, cfg%ddamp, 1.0e-12_real64)
end subroutine

! ... ddamp out-of-range, tmean default, tampli default, timref default
! ... temtoptab unallocated, allocation lifecycle, tembtab same
```

Add `ADD_TEST_SUITE(test_heat_phase0_fields_suite)` to `tests/unit/testSuites.inc`.

- [ ] **Step 5: Run tests — confirm FAIL**

```bash
pixi run -e test test-pfunit 2>&1 | tail -10
```
Expected: failures referencing missing fields on heat_config_t.

- [ ] **Step 6: Add the 6 fields to heat_config_t**

```fortran
! Phase 0 (ADR 0034): promoted from legacy globals to patch the
! silent-zero-defaults correctness gap on swcalt=1 (analytical method).

real(real64) :: ddamp  = 0.0_real64    !! Damping depth of temperature wave (L)
real(real64) :: tmean  = 0.0_real64    !! Prescribed mean annual surface temperature (degC)
real(real64) :: tampli = 0.0_real64    !! Amplitude of prescribed annual temperature wave (degC)
real(real64) :: timref = 0.0_real64    !! Time in year with top of prescribed sine wave (T)

! Boundary-condition tables — 2D (n_rows, 2): col 1 = time, col 2 = temperature
real(real64), allocatable :: temtoptab(:,:)   !! Soil-surface temperature vs time table
real(real64), allocatable :: tembtab(:,:)     !! Soil-bottom temperature vs time table
```

The flat 1D legacy form (`tembtab(mabbc*2)`) is interleaved per the existing pattern (see solute Phase 0 Task 1's `cseeptab` flattening). The 2D typed-config form is cleaner. The adapter flattens 2D → interleaved 1D.

- [ ] **Step 7: Add validators**

```fortran
call check_real_range(self%ddamp,  0.0_real64, 1.0e6_real64, "heat.ddamp",  errors)
call check_real_range(self%tmean, -100.0_real64, 100.0_real64, "heat.tmean", errors)
call check_real_range(self%tampli, 0.0_real64, 100.0_real64,  "heat.tampli", errors)
call check_real_range(self%timref, 0.0_real64, 366.0_real64,  "heat.timref", errors)
! Tables: when swcalt=1 and the top/bottom boundary uses table form, check allocated
```

- [ ] **Step 8: Add reads in read_heat_toml.f90**

```fortran
call get_optional_real_with_default(sec, "ddamp",  config%ddamp,  0.0_real64, "heat.ddamp",  errors)
call get_optional_real_with_default(sec, "tmean",  config%tmean,  0.0_real64, "heat.tmean",  errors)
call get_optional_real_with_default(sec, "tampli", config%tampli, 0.0_real64, "heat.tampli", errors)
call get_optional_real_with_default(sec, "timref", config%timref, 0.0_real64, "heat.timref", errors)
call get_optional_real_2d_array(sec, "temtoptab", config%temtoptab, "heat.temtoptab", errors)
call get_optional_real_2d_array(sec, "tembtab",   config%tembtab,   "heat.tembtab",   errors)
```

- [ ] **Step 9: Update apply_heat adapter**

```fortran
ddamp  = config%heat%ddamp
tmean  = config%heat%tmean
tampli = config%heat%tampli
timref = config%heat%timref

! Flatten 2D table → interleaved 1D legacy layout (verify pattern by reading solute Phase 0's cseeptab handling)
if (allocated(config%heat%temtoptab)) then
   ! ... interleaved flatten
end if
if (allocated(config%heat%tembtab)) then
   ! ... same
end if
```

- [ ] **Step 10: Run tests — confirm PASS**

```bash
pixi run -e test test-pfunit 2>&1 | tail -10
pixi run check-full
```

Expected: zero failures, count = baseline + however many tests added. check-full 5/5.

- [ ] **Step 11: Commit**

```bash
git add src/config/heat_config.f90 src/io/toml/read_heat_toml.f90 \
        src/io/toml/config_to_variables.f90 \
        tests/unit/config/test_heat_phase0_fields.pf tests/unit/testSuites.inc tests/unit/meson.build
git commit -m "$(cat <<'EOF'
feat(config): SS-HEAT Phase 0 — promote 6 missing heat physics fields

heat_config_t gains 4 scalars (ddamp, tmean, tampli, timref) and 2
2D tables (temtoptab, tembtab) supporting swcalt=1 (analytical
soil-temperature method). Each with typed default + validator +
TOML reader + adapter wiring.

Patches the silent-zero-defaults correctness gap: TOML runs with
swcalt=1 previously read zeroed legacy globals. None of the 5
regression cases activates this path, so check-full byte-identical
is preserved.

Spec: docs/superpowers/specs/2026-05-10-state-migration-heat-design.md
EOF
)"
```

---

# Phase 1 — State type, threading, dual-write, output reads, outheapar fix

### Task 2: Phase 1 — Create heat_state_t

**Files:**
- Create: `src/state/heat_state.f90`
- Create: `tests/unit/state/test_heat_state.pf`
- Modify: meson + testSuites.inc

- [ ] **Step 1: Identify exact field set from discovery Section 2**

Read `docs/superpowers/specs/2026-05-10-state-migration-heat-discovery.md` Section 2. The 13 owned globals include (verify the precise list):
- Per-node allocatable arrays: `tsoil(:)`, `heacap(:)`, `heacon(:)`, `rfcp(:)`
- Scalars: `tebot`, `tetop`
- Possibly others (check the discovery)

All 13 are instantaneous — no cohorts. Flat `heat_state_t`.

- [ ] **Step 2: Write failing tests**

```fortran
@test
subroutine test_heat_state_default_unallocated()
   use funit
   use heat_state_mod, only: heat_state_t
   type(heat_state_t) :: ht
   @assertFalse(allocated(ht%tsoil))
   @assertFalse(allocated(ht%heacap))
   @assertFalse(allocated(ht%heacon))
   @assertFalse(allocated(ht%rfcp))
end subroutine

@test
subroutine test_heat_state_scalar_defaults()
   use funit
   use iso_fortran_env, only: real64
   use heat_state_mod, only: heat_state_t
   type(heat_state_t) :: ht
   @assertEqual(0.0_real64, ht%tebot, 1.0e-12_real64)
   @assertEqual(0.0_real64, ht%tetop, 1.0e-12_real64)
end subroutine

@test
subroutine test_heat_state_array_allocation()
   use funit
   use iso_fortran_env, only: real64
   use heat_state_mod, only: heat_state_t
   type(heat_state_t) :: ht
   integer, parameter :: NUMNOD = 100

   allocate(ht%tsoil(NUMNOD));  ht%tsoil  = 0.0_real64
   allocate(ht%heacap(NUMNOD)); ht%heacap = 0.0_real64
   allocate(ht%heacon(NUMNOD)); ht%heacon = 0.0_real64
   allocate(ht%rfcp(NUMNOD));   ht%rfcp   = 1.0_real64

   @assertTrue(allocated(ht%tsoil))
   @assertEqual(NUMNOD, size(ht%tsoil))
end subroutine

@test
subroutine test_heat_state_independent_instances()
   use funit
   use iso_fortran_env, only: real64
   use heat_state_mod, only: heat_state_t
   type(heat_state_t) :: ht_a, ht_b
   ht_a%tebot = 5.0_real64
   ht_b%tebot = -3.0_real64
   @assertEqual(5.0_real64,  ht_a%tebot, 1.0e-12_real64)
   @assertEqual(-3.0_real64, ht_b%tebot, 1.0e-12_real64)
end subroutine
```

- [ ] **Step 3: Implement heat_state_t**

```fortran
!> @file heat_state.f90
!! Typed state record for the heat (soil-temperature) subsystem.
!! All 13 owned fields are instantaneous (no flzero* gating) — flat
!! layout with no cohort sub-records. See ADR 0034 / 2026-05-10
!! design spec.

module heat_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: heat_state_t

   type :: heat_state_t
      ! Per-node allocatable arrays — sized numnod, allocated by heat_init
      real(real64), allocatable :: tsoil(:)        !! soil temperature per node (degC)
      real(real64), allocatable :: heacap(:)       !! heat capacity per node (J/cm3/K)
      real(real64), allocatable :: heacon(:)       !! heat conductivity per node (J/cm/K/d)
      real(real64), allocatable :: rfcp(:)         !! frozen-conditions reduction factor per node

      ! Scalars
      real(real64) :: tebot = 0.0_real64           !! bottom-of-profile temperature (degC)
      real(real64) :: tetop = 0.0_real64           !! top-of-profile temperature (degC, under snow)

      ! Add others identified from discovery Section 2 to total 13 owned fields
   end type heat_state_t

end module heat_state_mod
```

The exact field list MUST come from discovery Section 2 — don't guess.

- [ ] **Step 4: Add to swap_state_t**

```fortran
type :: swap_state_t
   type(surfacewater_state_t) :: surfacewater
   type(drainage_state_t)     :: drainage
   type(solute_state_t)       :: solute
   type(heat_state_t)         :: heat
end type
```

- [ ] **Step 5: Build, test, commit**

```bash
pixi run -e test test-pfunit && pixi run check-full
git commit -m "feat(state): SS-HEAT Phase 1 Task 2 — heat_state_t typed state record"
```

### Task 3: Phase 1 — Plumb Temperature(task) with state + add heat_init

**Files:**
- Modify: `src/heat/temperature.f90`
- Modify: `src/core/swap.f90` — caller updates

- [ ] **Step 1: Inspect current Temperature signature**

```bash
grep -n "subroutine Temperature" src/heat/temperature.f90 | head
```

Currently `subroutine Temperature(task)`. Add `state`:
```fortran
subroutine Temperature(task, state)
   use swap_state_mod, only: swap_state_t
   type(swap_state_t), intent(inout) :: state
```

- [ ] **Step 2: Add heat_init(state) for allocation**

In `src/heat/temperature.f90` (or a new `src/heat/heat_init.f90` if cleaner — match pattern from drainage_init / solute_init), introduce:

```fortran
subroutine heat_init(state)
   use, intrinsic :: iso_fortran_env, only: real64
   use swap_state_mod, only: swap_state_t
   use Variables, only: numnod
   type(swap_state_t), intent(inout) :: state

   if (.not. allocated(state%heat%tsoil))  allocate(state%heat%tsoil(numnod))
   if (.not. allocated(state%heat%heacap)) allocate(state%heat%heacap(numnod))
   if (.not. allocated(state%heat%heacon)) allocate(state%heat%heacon(numnod))
   if (.not. allocated(state%heat%rfcp))   allocate(state%heat%rfcp(numnod))

   state%heat%tsoil  = 0.0_real64
   state%heat%heacap = 0.0_real64
   state%heat%heacon = 0.0_real64
   state%heat%rfcp   = 1.0_real64    ! NB: 1.0, not 0.0 — matches legacy initial value
end subroutine heat_init
```

- [ ] **Step 3: Wire heat_init from swap_main**

In `src/core/swap.f90`, after `drainage_init(state)` and `solute_init(state)`, add `call heat_init(state)`. Gate behind `if (flTemperature)` if appropriate (check existing init guards).

Update `Temperature(task)` callers in swap.f90 to pass `state`.

- [ ] **Step 4: Add use clause**

```fortran
use heat_state_mod  ! transitive via swap_state_mod, but if explicit access needed
```

- [ ] **Step 5: Build, verify**

```bash
pixi run check-full && pixi run -e test test-pfunit 2>&1 | tail
git commit -m "refactor(state): SS-HEAT Phase 1 Task 3 — Temperature(task) takes state, heat_init introduced"
```

### Task 4: Phase 1 — Heat compute dual-writes 13 owned fields

Heat compute (`Temperature(task)`, plus any related `FrozenCond` writes to heat-owned vars) writes the 13 owned globals. Add `state%heat%X = X` mirrors paralleling each global write.

**Files:**
- Modify: `src/heat/temperature.f90`
- Modify: `src/heat/frozencond.f90` (the `rfcp` write site, if heat-owned)

- [ ] **Step 1: Inventory write sites**

```bash
grep -nE "^\s+(tsoil|heacap|heacon|rfcp|tebot|tetop)\s*[=\(]" src/heat/temperature.f90 src/heat/frozencond.f90 | head -30
```

- [ ] **Step 2: Add dual-writes**

Use ASSOCIATE with `ht_*` prefix (per playbook):
```fortran
associate(ht_tsoil  => state%heat%tsoil, &
          ht_heacap => state%heat%heacap, &
          ht_heacon => state%heat%heacon, &
          ht_rfcp   => state%heat%rfcp,   &
          ht_tebot  => state%heat%tebot,  &
          ht_tetop  => state%heat%tetop)
   ! Wherever the body writes <global>, also write ht_<global>
end associate
```

- [ ] **Step 3: Verify, commit**

```bash
pixi run check-full && pixi run -e test test-pfunit 2>&1 | tail
git commit -m "refactor(state): SS-HEAT Phase 1 Task 4 — heat compute dual-writes state and globals"
```

### Task 5: Phase 1 — Output reads switch to state%heat + outheapar writeback fix

The plan's writeback hazard: `outheapar` calls `devries(thetadum, heacap, heacnd)` and writes `heacap`/`heacnd` from inside the OUTPUT path. Move that compute to where it belongs (likely `Temperature(task=1)` init), and `outheapar` becomes read-only on state.

**Files:**
- Modify: `src/io/swapoutput.f90` — `outtem`, `TemperatureOutput`, `outheapar` (writeback fix)
- Modify: `src/io/swap_csv_output.f90` — if `set_values` reads tsoil/heacap (verify)
- Modify: `src/heat/temperature.f90` — receive the moved devries call

- [ ] **Step 1: Find output reads of heat-owned fields**

```bash
grep -rEn "\b(tsoil|heacap|heacon|rfcp|tebot|tetop)\b" src/io/swapoutput.f90 src/io/swap_csv_output.f90 | head
```

- [ ] **Step 2: Migrate output reads**

Each `<heat-owned-global>` read becomes `state%heat%<global>`. Add `state` arg if not already present (most output routines already take state from prior migrations).

Drop migrated symbols from `use variables, only:` clauses.

- [ ] **Step 3: Fix outheapar writeback**

Find `outheapar` in swapoutput.f90. Locate the `call devries(thetadum, heacap, heacnd)` site. Move the `devries` call to the appropriate compute site (likely `Temperature(task=1)` or a heat-init helper). The `outheapar` body becomes read-only — reads `state%heat%heacap` / `state%heat%heacon`.

If the moved `devries` call needs different inputs at the compute site than at the output site (e.g., different `thetadum`), inspect carefully. The original output-side call may have used the most-recent `theta` from soilhydraulics; the compute-side equivalent uses the same. Verify by inspection.

- [ ] **Step 4: Verify, commit**

```bash
pixi run check-full && pixi run -e test test-pfunit 2>&1 | tail
git commit -m "refactor(state): SS-HEAT Phase 1 Task 5 — output reads state%heat; outheapar writeback fixed"
```

---

# Phase 2 — Cross-subsystem readers, soilhydraulics rfcp fix, global removal

### Task 6: Phase 2 — Migrate compute readers of tsoil

10 compute files read `tsoil` from globals. Migrate each. Per discovery:
- `src/crop/rootextraction.f90`
- `src/crop/oxygenstress.f90`
- `src/crop/management_soil.f90`
- `src/crop/cropgrowth.f90`
- `src/solute/solute.f90`
- `src/atmosphere/snow.f90`
- `src/utils/soilhydraulicsutils.f90`
- `src/soil/soilhydraulics.f90`
- `src/boundary/boundtop.f90`
- `src/boundary/boundbottom.f90`

Most already take `state` from prior migrations. The migration is mechanical: replace `tsoil(node)` with `state%heat%tsoil(node)` and drop `tsoil` from `use variables, only:` clauses.

**Caveat:** these 10 readers currently consume `tsoil=0.0` when `flTemperature=false` (pre-existing bug, out-of-scope for this arc). The migration preserves that behavior — `state%heat%tsoil` is still 0 when heat is disabled.

- [ ] **Step 1: Find tsoil reads in each file**

```bash
for f in src/crop/rootextraction.f90 src/crop/oxygenstress.f90 src/crop/management_soil.f90 src/crop/cropgrowth.f90 src/solute/solute.f90 src/atmosphere/snow.f90 src/utils/soilhydraulicsutils.f90 src/soil/soilhydraulics.f90 src/boundary/boundtop.f90 src/boundary/boundbottom.f90; do
  echo "=== $f ==="
  grep -nE "\btsoil\b" "$f" | head -5
done
```

- [ ] **Step 2: Migrate each file**

Per file: confirm `state` is in scope; replace `tsoil` reads with `state%heat%tsoil`; drop from `use variables, only:` if no other use.

If a file doesn't yet take `state`, add it (`type(swap_state_t), intent(in) :: state`) and update its caller. Match the optional-state-arg pattern from solute Phase 2 if the routine has many call chains.

- [ ] **Step 3: Verify, commit**

```bash
pixi run check-full && pixi run -e test test-pfunit 2>&1 | tail
git commit -m "refactor(state): SS-HEAT Phase 2 Task 6 — 10 cross-subsystem tsoil readers migrate to state"
```

### Task 7: Phase 2 — Fix rfcp co-write in soilhydraulics

`soilhydraulics.f90:884` resets `rfcp = 1.0d0` on every Richards timestep. Migrate to `state%heat%rfcp`.

**Files:**
- Modify: `src/soil/soilhydraulics.f90`

- [ ] **Step 1: Locate the rfcp write**

```bash
grep -nE "\brfcp\b" src/soil/soilhydraulics.f90 | head
```

- [ ] **Step 2: Migrate**

soilhydraulics already takes `state`. Replace:
```fortran
rfcp = 1.0d0
```
with:
```fortran
if (allocated(state%heat%rfcp)) state%heat%rfcp = 1.0_real64
```

The `allocated` guard handles the case where heat is disabled (`flTemperature=false`) — no allocation, no write.

Drop `rfcp` from `use variables, only:` clause in soilhydraulics.f90.

- [ ] **Step 3: Verify, commit**

```bash
pixi run check-full && pixi run -e test test-pfunit 2>&1 | tail
git commit -m "refactor(state): SS-HEAT Phase 2 Task 7 — soilhydraulics rfcp reset migrates to state%heat"
```

### Task 8: Phase 2 — Drop dual-write, state authoritative

After Task 6 + 7, all consumers of heat-owned globals read from state. Heat compute can stop writing legacy globals.

**Files:**
- Modify: `src/heat/temperature.f90` — remove legacy global writes

- [ ] **Step 1: Find dual-write sites**

```bash
grep -nE "^\s+(tsoil|heacap|heacon|rfcp|tebot|tetop)\s*[=\(]" src/heat/temperature.f90
grep -nE "^\s+state%heat%(tsoil|heacap|heacon|rfcp|tebot|tetop)\s*[=\(]" src/heat/temperature.f90
```

- [ ] **Step 2: Drop legacy global writes; keep state writes**

Each dual-write pair: delete the legacy global write. Keep the `state%heat%<global>` write.

Trim `use variables, only:` clauses for migrated symbols.

- [ ] **Step 3: Verify check-full byte-identical (integration gate)**

```bash
pixi run check-full
pixi run -e test test-pfunit 2>&1 | tail
```

If a case fails: a missed reader. Use the audit greps from the playbook to find it.

- [ ] **Step 4: Commit**

```bash
git commit -m "refactor(state): SS-HEAT Phase 2 Task 8 — drop dual-write; heat state authoritative"
```

### Task 9: Phase 2 — Comment out 13 globals

**Files:**
- Modify: `src/core/variables.f90` — comment out 13 declarations with provenance markers
- Modify: `src/core/initialize.f90` — drop matching zero-init lines

- [ ] **Step 1: Verify zero readers via 5-category grep**

```bash
for v in tsoil heacap heacon rfcp tebot tetop; do
  echo "=== $v ==="
  grep -rEn "use variables.*\b$v\b" src/ --include="*.f90" \
    | grep -v "src/core/variables.f90" \
    | grep -v "src/state/" \
    | grep -v "src/heat/" \
    | grep -v "src/core/initialize.f90" \
    | head
  grep -rEn "\b$v\b" src/ --include="*.f90" \
    | grep -v "src/core/variables.f90" \
    | grep -v "src/core/initialize.f90" \
    | grep -v "src/state/" \
    | grep -v "state%" \
    | grep -v "config%" \
    | head -3
done
```

Both must be empty before deletion. Investigate any remaining hits.

- [ ] **Step 2: Comment out declarations in variables.f90**

Match the surfacewater/drainage/solute pattern:
```fortran
! real(8) tsoil(macp)        ! Moved to heat_state_t%tsoil (ADR 0034)
! real(8) heacap(macp)       ! Moved to heat_state_t%heacap
! ... etc
```

- [ ] **Step 3: Drop initialize.f90 lines**

```bash
grep -nE "^\s+(tsoil|heacap|heacon|rfcp|tebot|tetop)\s*=" src/core/initialize.f90
```

Delete or comment out matching zero-init lines.

- [ ] **Step 4: Build, verify, commit**

```bash
pixi run check-full && pixi run -e test test-pfunit 2>&1 | tail
git commit -m "refactor(state): SS-HEAT Phase 2 Task 9 — comment out 13 heat globals"
```

### Task 10: Phase 2 — ADR 0034 finalize + index

**Files:**
- Create: `docs/adr/0034-state-migration-heat.md`
- Modify: `docs/adr/index.md`

- [ ] **Step 1: Author ADR 0034**

Match the structure of ADR 0030/0031/0032/0033. Required content:
- Status: accepted
- Date
- Migration #4
- Context: heat as the smallest migration so far; no cohorts
- Decision: flat heat_state_t; Phase 0 config promotion; outheapar writeback fix; rfcp co-write fix
- Phase 0 outcomes (6 fields promoted)
- Phase 1 outcomes (state type, threading, dual-write)
- Phase 2 outcomes (cross-subsystem readers, globals removed)
- Architectural learnings: first migration without cohorts; outheapar fix sets a precedent for moving compute out of output routines
- Known issue: 10 ungated tsoil readers in compute (pre-existing)
- References

- [ ] **Step 2: Add to index**

```markdown
- [ADR 0034 — Heat state-type migration](0034-state-migration-heat.html) — Fourth subsystem migrated. First with all-instantaneous fields (no cohort sub-records — flat heat_state_t). Temperature(task) plumbed with state. Phase 0 promoted 6 missing config fields for swcalt=1. outheapar writeback hazard fixed (devries call moved from output to compute). rfcp co-write in soilhydraulics migrated. 13 owned globals removed from variables.f90.
```

- [ ] **Step 3: Final verification**

```bash
pixi run check-full && pixi run -e test test-pfunit 2>&1 | tail
```

- [ ] **Step 4: Commit**

```bash
git commit -m "docs(adr): SS-HEAT Phase 2 Task 10 — ADR 0034 finalize"
```

---

## Self-Review Notes

- **Spec coverage:** D1 (scope) → all tasks. D2 (flat state shape) → Task 2. D3 (aggregator) → Task 2. D4 (Temperature(task) threading) → Task 3. D5 (rfcp co-write fix) → Task 7. D6 (outheapar writeback fix) → Task 5. D7 (Phase 0 config gap) → Task 1. D8 (10 ungated tsoil readers documented) → Task 10 ADR.
- **Lessons-learned applied:** Section 2 categorization upfront; activity-gate awareness (no cumulatives so trivially satisfied); ALL-alias-form grep before deletion; ht_* ASSOCIATE prefix; optional state args where appropriate.
- **Phase boundaries**: Phase 0 (Task 1) → Phase 1 (Tasks 2-5) → Phase 2 (Tasks 6-10). Each phase ships independently with check-full as the integration gate.
- **No physics changes**: ASSOCIATE preserves variable names; bodies of Temperature/FrozenCond/FrozenBounds change only in the writeback-fix routes.
- **Branch policy:** All commits on `refactor/heat-state`. After Task 10, ready to merge to development.
