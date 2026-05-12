# TimeControl State-Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Carve 61 runtime-state TimeControl globals into `state%timecontrol`; preserve byte-identical check-full at every commit.

**Architecture:** Flat `timecontrol_state_t` aggregated under `swap_state_t`. No `timecontrol_init` — init inlined in `TimeControl(case=1)`. Three signature touches: `TimeControl` intent(in) → intent(inout), `IterTime` gains state arg, DLL re-init path gets state write. Pre-init transient buffer for `iyear`/`imonth`/`dt`. Strategy B compile-driven retirement at the end.

**Tech Stack:** Fortran 2008+ derived types, ASSOCIATE, pixi+meson, pFUnit, check-full regression. Branch: `refactor/timecontrol-state` (branch off `development`).

---

## Standing instructions for every task

- **Branch.** All work on `refactor/timecontrol-state`. Verify before commit: `pwd && git rev-parse --abbrev-ref HEAD` → expected `refactor/timecontrol-state`.
- **Verification gate.** Every commit must pass `pixi run check-full` (5/5 byte-identical) AND `pixi run -e test test-pfunit` (727+N/0/0 where N = new tests).
- **Hupselbrook canary.** If any regression case hangs >5 seconds at 99% CPU with no output, kill and abort (`pkill -9 -f "test_output_regression"; pkill -9 -f "/builddir/swap"`); `git checkout -- .`; investigate.
- **Pre-flight dual-write coverage check** (atmosphere lesson, applied since A-2.1). Before each Phase 2 reader cutover task, verify each migrated field has non-zero state writes:
  ```bash
  grep -rn "state%timecontrol%<field>\s*=" src/ --include="*.f90" | grep -v "0\.0\|0_real64\|\.false\.\|0_int32\|0$"
  ```
  Must return at least one non-zero write per field. If gap, fix the dual-write in Phase 1 routines first.

---

## The 61 owned fields

Pulled from refreshed discovery Section 2 (post-purge). Grouped by cadence for state-type layout. Use these exact names — they match the legacy globals minus the `tc_` prefix where applicable.

**Clock state (16 fields)**
- Real: `t`, `t1900`, `tcum`, `dt`, `dtmin0`, `dtprevious`, `dtEvent`, `tEvent`, `tchange`, `tnext`, `tinit`, `tmptimestart`, `tmptimeend`, `dtcrit`
- Integer: `daynr`, `daycum`

**Schedule state (19 fields)**
- Integer: `iyear`, `imonth`, `daymonth`, `monthstart`, `monthend`, `yearstart`, `yearend`, `daycrop`, `yearcrop`, `isteps`, `MaxIterIterTime`, `ioutdat`, `cntper`, `nprintcount`, `nprintsecond`, `nprintdayfirst`, `outper`, `idstart`, `idend`
- Real: `cumirr` (cumulative irrigation — verify scope; may be irrigation-owned)

**Boolean state (26 fields)**
- `flDayStart`, `flDayEnd`, `flOutputShort`, `flprintdt`, `flprintshort`, `flprint`, `flprintsecond`, `floutputshort`, `flheader`, `flheadirg`, `flIrg1Start`, `flIrg2Start`, `flbaloutput`, `flintrout`, `flcumuout`, `flBegOfDay`, `flEndOfDay`, `flLastTimeStep`, `flFirstTimeStep`, `flCropCalendar`, `flCropEmergence`, `flCropEnd`, `flCropEvent`, `flMet`, `flMetDet`, `flUpdMetDet`

If any of these turns out to belong to another subsystem (e.g., `flCropCalendar` may be crop-owned), the implementer flags it during Task 1 and excludes it. The 61 count is approximate; final inventory comes from the implementer's confirmation grep.

---

## File Structure

**Create:**
- `src/state/timecontrol_state.f90` — module + type definition

**Modify:**
- `src/state/swap_state.f90` — add `type(timecontrol_state_t) :: timecontrol`
- `src/core/timecontrol.f90` — bump intent, plumb IterTime, dual-write all 61 fields
- `src/core/swap.f90` — update TimeControl/IterTime call sites; DLL re-init state write
- `src/io/toml/config_to_variables.f90` — keep iyear/imonth/dt pre-init writes (transient pattern)
- `src/core/variables.f90` — Strategy B comment-out of 61 globals (Phase 2 end)
- `src/core/initialize.f90` — drop dead zero-inits after retirement
- `tests/unit/state/test_timecontrol_state.pf` — new pFUnit suite
- `tests/unit/meson.build` and `tests/unit/testSuites.inc` — register the new test
- `meson.build` — add `src/state/timecontrol_state.f90` to state-module sources

**Reader-cutover modifies (Phase 2)** — ordered by site count:
- `src/soil/waterbalance.f90` (69 sites — heaviest)
- `src/io/swapoutput.f90` (53)
- `src/drainage/surfacewater.f90` (34)
- `src/soil/soilhydraulics.f90` (33)
- `src/io/readmeteo.f90` (32)
- `src/crop/cropgrowth.f90` (31)
- `src/atmosphere/meteodt.f90` (25)
- `src/atmosphere/meteoday.f90` (22)
- `src/core/initialize.f90` (20)
- `src/atmosphere/et.f90` (19)
- `src/io/swap_csv_output.f90` (~15)
- `src/atmosphere/precipitation.f90` (~10)
- `src/heat/temperature.f90` (~10)
- `src/heat/frozencond.f90` (~8)
- `src/atmosphere/snow.f90` (0 — already migrated per atmosphere arc; skip)
- `src/solute/solute.f90` (~6)
- `src/crop/rootextraction.f90` (~5)
- `src/utils/soilhydraulicsutils.f90` (~5)
- `src/utils/surfacewaterutils.f90` (~4)
- `src/crop/irrigation.f90` (~4)
- `src/solute/agetracer.f90` (~3)
- `src/crop/tillage.f90` (~2)

---

## Phase 1 — state type + signature plumbing + home-tree dual-write

### Task 1: Create `timecontrol_state_t`

**Files:**
- Create: `src/state/timecontrol_state.f90`
- Modify: `src/state/swap_state.f90`
- Modify: `meson.build`
- Create: `tests/unit/state/test_timecontrol_state.pf`
- Modify: `tests/unit/meson.build`, `tests/unit/testSuites.inc`

- [ ] **Step 1: Enumerate exact field set from variables.f90**

Run:
```bash
grep -nE "^[[:space:]]+(real|integer|logical)[^!]*\btc_" src/core/variables.f90 > /tmp/tc_inventory.txt
grep -nE "^[[:space:]]+(real|integer|logical)[^!]*\b(daynr|daycum|imonth|isteps|nprintday|ioutdat|cntper|outper|cumtime|t1900|tcum|dt|dtmin0|dtprevious|dtEvent|tEvent|tchange|tnext|tinit|dtcrit|iyear|daymonth|monthstart|monthend|yearstart|yearend|daycrop|yearcrop|MaxIterIterTime|nprintcount|nprintsecond|nprintdayfirst|idstart|idend|flDayStart|flDayEnd|flOutputShort|flprintdt|flprintshort|flprint|flprintsecond|floutputshort|flheader|flheadirg|flIrg1Start|flIrg2Start|flbaloutput|flintrout|flcumuout|flBegOfDay|flEndOfDay|flLastTimeStep|flFirstTimeStep|flCropCalendar|flCropEmergence|flCropEnd|flCropEvent|flMet|flMetDet|flUpdMetDet)\b" src/core/variables.f90 >> /tmp/tc_inventory.txt
cat /tmp/tc_inventory.txt | wc -l
```

Expected: ~61 lines (plus or minus a few — confirm and proceed with what's actually present).

- [ ] **Step 2: Write the failing test**

Create `tests/unit/state/test_timecontrol_state.pf`:
```fortran
@test
subroutine test_default_values_are_zero()
   use funit
   use timecontrol_state_mod, only: timecontrol_state_t
   use iso_fortran_env, only: real64
   type(timecontrol_state_t) :: tc
   @assertEqual(0.0_real64, tc%t, 1.0e-12_real64)
   @assertEqual(0.0_real64, tc%t1900, 1.0e-12_real64)
   @assertEqual(0, tc%daynr)
   @assertEqual(0, tc%daycum)
   @assertEqual(0, tc%iyear)
   @assertFalse(tc%flDayStart)
   @assertFalse(tc%flDayEnd)
end subroutine test_default_values_are_zero

@test
subroutine test_aggregator_access_via_swap_state()
   use funit
   use swap_state_mod, only: swap_state_t
   type(swap_state_t) :: state
   state%timecontrol%daynr = 42
   state%timecontrol%iyear = 2024
   state%timecontrol%flDayStart = .true.
   @assertEqual(42, state%timecontrol%daynr)
   @assertEqual(2024, state%timecontrol%iyear)
   @assertTrue(state%timecontrol%flDayStart)
end subroutine test_aggregator_access_via_swap_state

@test
subroutine test_two_instances_independent()
   use funit
   use timecontrol_state_mod, only: timecontrol_state_t
   use iso_fortran_env, only: real64
   type(timecontrol_state_t) :: a, b
   a%dt = 0.125_real64
   b%dt = 0.250_real64
   @assertEqual(0.125_real64, a%dt, 1.0e-12_real64)
   @assertEqual(0.250_real64, b%dt, 1.0e-12_real64)
end subroutine test_two_instances_independent

@test
subroutine test_real_scalar_round_trip()
   use funit
   use timecontrol_state_mod, only: timecontrol_state_t
   use iso_fortran_env, only: real64
   type(timecontrol_state_t) :: tc
   tc%t1900 = 47482.5_real64
   tc%tcum = 365.0_real64
   @assertEqual(47482.5_real64, tc%t1900, 1.0e-12_real64)
   @assertEqual(365.0_real64, tc%tcum, 1.0e-12_real64)
end subroutine test_real_scalar_round_trip

@test
subroutine test_integer_round_trip()
   use funit
   use timecontrol_state_mod, only: timecontrol_state_t
   type(timecontrol_state_t) :: tc
   tc%daynr = 200
   tc%imonth = 7
   tc%ioutdat = 5
   @assertEqual(200, tc%daynr)
   @assertEqual(7, tc%imonth)
   @assertEqual(5, tc%ioutdat)
end subroutine test_integer_round_trip

@test
subroutine test_logical_round_trip()
   use funit
   use timecontrol_state_mod, only: timecontrol_state_t
   type(timecontrol_state_t) :: tc
   tc%flDayStart = .true.
   tc%flDayEnd = .false.
   tc%flLastTimeStep = .true.
   @assertTrue(tc%flDayStart)
   @assertFalse(tc%flDayEnd)
   @assertTrue(tc%flLastTimeStep)
end subroutine test_logical_round_trip
```

Add `tests/unit/state/test_timecontrol_state.pf` to `tests/unit/meson.build` (`pf_files` array) and the suite name to `tests/unit/testSuites.inc`.

- [ ] **Step 3: Run the test (should fail — module doesn't exist)**

```bash
pixi run -e test test-pfunit 2>&1 | tail -10
```

Expected: compile error "Cannot open module file 'timecontrol_state_mod.mod'".

- [ ] **Step 4: Write the module**

Create `src/state/timecontrol_state.f90`. Define `timecontrol_state_t` with all confirmed fields (from Step 1's inventory), default-initialized:
- Real fields: `= 0.0_real64`
- Integer fields: `= 0`
- Logical fields: `= .false.`

Use the `use, intrinsic :: iso_fortran_env, only: real64` idiom matching `atmosphere_state.f90`.

In `src/state/swap_state.f90` add:
```fortran
use timecontrol_state_mod, only: timecontrol_state_t
```
and inside `type :: swap_state_t`:
```fortran
type(timecontrol_state_t) :: timecontrol
```
Position after `tillage` (the most recently added member).

In root `meson.build`, add `src/state/timecontrol_state.f90` to the state-module sources block (next to `tillage_state.f90`).

- [ ] **Step 5: Run the test (should pass)**

```bash
pixi run -e test test-pfunit 2>&1 | tail -10
```

Expected: 733 OK / 0 fail / 0 disabled (was 727; +6 new tests).

- [ ] **Step 6: Run check-full**

```bash
pixi run check-full 2>&1 | tail -5
```

Expected: `Results: 5 passed, 0 failed`. Byte-identical — no field migration yet, just type creation.

- [ ] **Step 7: Commit**

```bash
git add src/state/timecontrol_state.f90 src/state/swap_state.f90 meson.build \
        tests/unit/state/test_timecontrol_state.pf tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "$(cat <<'EOF'
feat(state): SS-TC Task 1 — timecontrol_state_t (61 runtime fields)

Introduces timecontrol_state_t with 16 clock + 19 schedule + 26
boolean fields (61 total). Flat layout — no cohorts (TC owns the
flzero reset gates rather than consuming them).

Aggregated under swap_state_t as state%timecontrol.

No field migration this commit — type only. Tasks 2-3 wire TC compute
to dual-write; Phase 2 cuts over readers.

pFUnit: 6 new tests covering defaults, aggregator access, round-trips.

Spec: docs/superpowers/specs/2026-05-12-state-migration-timecontrol-design.md
EOF
)"
```

---

### Task 2: Bump `TimeControl` intent + plumb `IterTime`

**Files:**
- Modify: `src/core/timecontrol.f90:1` (subroutine declaration)
- Modify: `src/core/swap.f90` (single TimeControl call site; IterTime call site)

- [ ] **Step 1: Locate current TimeControl signature**

```bash
grep -n "subroutine TimeControl\|subroutine IterTime" src/core/timecontrol.f90
grep -n "call TimeControl\|call IterTime" src/core/swap.f90
```

Expected output identifies the signatures and call sites.

- [ ] **Step 2: Bump TimeControl intent**

In `src/core/timecontrol.f90`, change:
```fortran
type(swap_state_t), intent(in) :: state
```
to:
```fortran
type(swap_state_t), intent(inout) :: state
```
in `TimeControl(task, state)`. No other change in this step.

- [ ] **Step 3: Plumb IterTime**

In `src/core/timecontrol.f90`, change `IterTime(task)` to `IterTime(task, state)`:
```fortran
subroutine IterTime(task, state)
   use swap_state_mod, only: swap_state_t
   integer, intent(in) :: task
   type(swap_state_t), intent(inout) :: state
   ! ... existing body unchanged
end subroutine
```

In `src/core/swap.f90` at the IterTime call site, change `call IterTime(N)` to `call IterTime(N, state)`. (Single call site — verify with grep above.)

- [ ] **Step 4: Run check-full**

```bash
pixi run check-full 2>&1 | tail -5
```

Expected: `Results: 5 passed, 0 failed`. No behavior change — only signatures.

- [ ] **Step 5: Run pFUnit**

```bash
pixi run -e test test-pfunit 2>&1 | tail -5
```

Expected: 733/0/0 unchanged.

- [ ] **Step 6: Commit**

```bash
git add src/core/timecontrol.f90 src/core/swap.f90
git commit -m "$(cat <<'EOF'
refactor(state): SS-TC Task 2 — TimeControl intent(inout) + IterTime state plumbing

TimeControl(task, state) intent promoted from intent(in) to
intent(inout) (was intent(in) since SS-SWC S-2.12B; now needs to write
state). IterTime(task) gains state arg as intent(inout) — single
call site in swap.f90 updated.

No state writes added yet; signature bumps only. Task 3 starts the
dual-write within TC.

Spec: docs/superpowers/specs/2026-05-12-state-migration-timecontrol-design.md
EOF
)"
```

---

### Task 3: Dual-write 61 fields in TimeControl

**Files:**
- Modify: `src/core/timecontrol.f90` (all TC bodies)

- [ ] **Step 1: Inventory write sites by field**

```bash
for v in t t1900 tcum dt dtmin0 dtprevious dtEvent tEvent tchange tnext tinit dtcrit daynr daycum iyear imonth daymonth monthstart monthend yearstart yearend daycrop yearcrop isteps MaxIterIterTime ioutdat cntper nprintcount nprintsecond nprintdayfirst outper idstart idend flDayStart flDayEnd flOutputShort flprintdt flprintshort flprint flprintsecond floutputshort flheader flheadirg flIrg1Start flIrg2Start flbaloutput flintrout flcumuout flBegOfDay flEndOfDay flLastTimeStep flFirstTimeStep flCropCalendar flCropEmergence flCropEnd flCropEvent flMet flMetDet flUpdMetDet; do
  hits=$(grep -nE "^[^!]*\b$v\b\s*=" src/core/timecontrol.f90 2>/dev/null | head -3)
  [ -n "$hits" ] && echo "=== $v ===" && echo "$hits"
done
```

Each write site needs a sibling `state%timecontrol%X = X` mirror.

- [ ] **Step 2: Add ASSOCIATE block in each TC case body**

Open the body of each `case (1)`, `case (2)`, `case (3)`, etc., and wrap with:
```fortran
associate( &
   tc_daynr     => state%timecontrol%daynr, &
   tc_daycum    => state%timecontrol%daycum, &
   tc_dt        => state%timecontrol%dt, &
   tc_t1900     => state%timecontrol%t1900, &
   ! ... etc — alias only the fields actually written in THIS case
)
   ! body unchanged at the syntactic level for legacy reads/writes;
   ! after each `daynr = X` write, add a sibling `tc_daynr = daynr`
end associate
```

For the most-written fields (dt, t1900, daynr, daycum, iyear, imonth, ioutdat, flDayStart), use the alias. For one-off fields, use direct `state%timecontrol%X = ...` after the legacy write.

- [ ] **Step 3: Add pair-write after each legacy write**

Pattern:
```fortran
! Legacy:
daynr = daynr + 1

! After (added mirror):
daynr = daynr + 1
tc_daynr = daynr     ! via ASSOCIATE alias
```

Or for fields outside an ASSOCIATE:
```fortran
flLastTimeStep = .true.
state%timecontrol%flLastTimeStep = flLastTimeStep
```

- [ ] **Step 4: Verify dual-write coverage**

```bash
for v in t t1900 tcum dt dtmin0 dtprevious daynr daycum iyear imonth ioutdat outper flDayStart flDayEnd flheader; do
  count=$(grep -cE "state%timecontrol%$v\s*=|tc_$v\s*=" src/core/timecontrol.f90 2>/dev/null)
  echo "$v: $count state-side writes"
done
```

Every field that the original code writes inside TC should have at least one state-side write here.

- [ ] **Step 5: Run check-full**

```bash
pixi run check-full 2>&1 | tail -5
```

Expected: `Results: 5 passed, 0 failed`. Byte-identical — legacy writes still happen; we only added mirrors.

- [ ] **Step 6: Run pFUnit**

```bash
pixi run -e test test-pfunit 2>&1 | tail -5
```

Expected: 733/0/0 unchanged.

- [ ] **Step 7: Commit**

```bash
git add src/core/timecontrol.f90
git commit -m "$(cat <<'EOF'
refactor(state): SS-TC Task 3 — TimeControl dual-write 61 fields

Every write to a TC-owned legacy global inside TimeControl(case=N)
and IterTime(task=N) now mirrors into state%timecontrol. ASSOCIATE
blocks (tc_* prefix) used in dense case bodies; explicit
state%timecontrol%X writes for sparse one-offs.

Legacy globals still written — dual-write era is live until Phase 2
retirement. <N> pair-writes added.

Spec: docs/superpowers/specs/2026-05-12-state-migration-timecontrol-design.md
EOF
)"
```

---

### Task 4: Pre-init transient pattern for `iyear` / `imonth` / `dt`

**Files:**
- Modify: `src/core/swap.f90` (DLL re-init path ~line 592-598; first-init path)
- Modify: `src/io/toml/config_to_variables.f90` (verify pre-init writes still happen for legacy)

- [ ] **Step 1: Locate pre-init writes**

```bash
grep -nE "iyear\s*=|imonth\s*=|^\s*dt\s*=" src/io/toml/config_to_variables.f90 | head
grep -nE "iyear\s*=|imonth\s*=" src/core/swap.f90 | head -5
```

Expected: `config_to_variables.f90` writes iyear/imonth/dt before state allocation. `swap.f90:592-598` (DLL re-init path) writes iyear directly.

- [ ] **Step 2: Add state seeds after `soilwater_init` returns**

The pattern from atmosphere A-2.6: legacy writes happen at config-time (before state allocation); right after `state%timecontrol` exists (no init routine — state struct allocated as part of `swap_state_t` declaration), copy legacy → state at the start of `TimeControl(case=1)`. Since the swap_state_t scalars default to zero, the COPY happens in TC(case=1) before any TC logic runs.

In `src/core/timecontrol.f90` at the very start of `case (1)`:
```fortran
case (1)
   ! Transient seed from legacy globals written by config_to_variables
   ! (pre-init pattern — atmosphere A-2.6 lesson).
   state%timecontrol%iyear = iyear
   state%timecontrol%imonth = imonth
   state%timecontrol%dt = dt
   ! ... continue with existing case(1) body
```

- [ ] **Step 3: DLL re-init state write**

In `src/core/swap.f90:592-598` (verify line numbers — DLL re-init path that mutates iyear externally), after the legacy `iyear = ...` write add:
```fortran
state%timecontrol%iyear = iyear
```

- [ ] **Step 4: Run check-full**

```bash
pixi run check-full 2>&1 | tail -5
```

Expected: 5/5 byte-identical. Pre-init seeds are no-ops at this stage (state still mirrors legacy from Task 3 dual-writes).

- [ ] **Step 5: Commit**

```bash
git add src/core/timecontrol.f90 src/core/swap.f90
git commit -m "$(cat <<'EOF'
refactor(state): SS-TC Task 4 — pre-init transient pattern for iyear/imonth/dt

config_to_variables writes iyear/imonth/dt before state%timecontrol
is populated (pre-init pattern). TC(case=1) now copies the three
legacy values into state at the start of init.

DLL re-init path (swap.f90:592-598) gets a state%timecontrol%iyear
mirror write alongside the existing legacy write.

Mirrors atmosphere A-2.6's transient-buffer pattern.

Spec: docs/superpowers/specs/2026-05-12-state-migration-timecontrol-design.md D7/D8
EOF
)"
```

---

### Task 5: swapoutput co-write `flheader` / `flheadirg` / `flIrg1Start`

**Files:**
- Modify: `src/io/swapoutput.f90` (write sites for these 3 flags)

- [ ] **Step 1: Locate writes**

```bash
grep -nE "^[^!]*\b(flheader|flheadirg|flIrg1Start)\s*=" src/io/swapoutput.f90
```

Expected: 3-6 write sites total.

- [ ] **Step 2: Add state mirrors**

For each write `flheader = .false.`, add `state%timecontrol%flheader = .false.` on the next line. Same for the other two flags.

- [ ] **Step 3: Run check-full**

```bash
pixi run check-full 2>&1 | tail -5
```

Expected: 5/5 byte-identical.

- [ ] **Step 4: Commit**

```bash
git add src/io/swapoutput.f90
git commit -m "$(cat <<'EOF'
refactor(state): SS-TC Task 5 — swapoutput co-write flheader/flheadirg/flIrg1Start

swapoutput.f90 dual-writes the three output-state flags it co-owns
with TC: flheader, flheadirg, flIrg1Start. Each legacy write now
has a sibling state%timecontrol%X mirror.

Spec: docs/superpowers/specs/2026-05-12-state-migration-timecontrol-design.md D9
EOF
)"
```

---

## Phase 2 — reader cutover (heaviest readers first)

### Task 6: waterbalance.f90 reader cutover (69 sites — heaviest)

**Files:**
- Modify: `src/soil/waterbalance.f90`

- [ ] **Step 1: Pre-flight dual-write coverage check**

```bash
for v in dt t1900 daynr daycum outper flDayStart flDayEnd flprintshort flbaloutput; do
  count=$(grep -cE "state%timecontrol%$v\s*=" src/ -r --include="*.f90" | awk -F: '{s+=$NF} END {print s}')
  echo "$v: $count state-side writes"
done
```

Expected: every field has ≥1 non-zero state write from Phase 1.

- [ ] **Step 2: Identify reader sites**

```bash
grep -nE "\b(t|t1900|tcum|dt|daynr|daycum|iyear|imonth|outper|cntper|ioutdat|flDayStart|flDayEnd|flbaloutput|flprintshort|flprint|flcumuout|flintrout)\b" src/soil/waterbalance.f90 | head -50
```

Filter out comments and dummy-arg lines. Expected ~69 read references.

- [ ] **Step 3: Add ASSOCIATE block**

In each subroutine in waterbalance.f90 that reads TC fields (integral, fluxes, calcgwl, watstor, checkmassbal), add an ASSOCIATE block aliasing the fields it reads:
```fortran
associate( &
   tc_dt        => state%timecontrol%dt, &
   tc_daynr     => state%timecontrol%daynr, &
   tc_outper    => state%timecontrol%outper, &
   ! ... etc
)
   ! body: replace bare reads with tc_*
end associate
```

- [ ] **Step 4: Replace reads**

For each read site, replace `dt` with `tc_dt` (or `state%timecontrol%dt` if outside an ASSOCIATE).

- [ ] **Step 5: Run check-full**

```bash
pixi run check-full 2>&1 | tail -5
```

Expected: 5/5 byte-identical.

- [ ] **Step 6: Commit**

```bash
git add src/soil/waterbalance.f90
git commit -m "$(cat <<'EOF'
refactor(state): SS-TC Task 6 — waterbalance reader cutover (69 sites)

Migrates ~69 TC field reads in waterbalance.f90 (integral, fluxes,
calcgwl, watstor, checkmassbal) to state%timecontrol via ASSOCIATE
tc_* aliases.

Pre-flight dual-write coverage verified for all migrated fields.

Spec: docs/superpowers/specs/2026-05-12-state-migration-timecontrol-design.md D12/D13
EOF
)"
```

---

### Task 7: swapoutput.f90 reader cutover (53 sites)

**Files:**
- Modify: `src/io/swapoutput.f90`

Same pattern as Task 6: pre-flight check, ASSOCIATE block, replace reads, byte-identical, commit. Routines: swapoutput, soilwateroutput, SoluteOutput, TemperatureOutput, SurfaceWaterOutput, SnowOutput, OutCropFixed, OutWofost, OutGrass, outinc, outrot, outtem, agetracer outputs, writehead, utilities.

Commit message:
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-TC Task 7 — swapoutput reader cutover (53 sites)

Migrates ~53 TC field reads in swapoutput.f90 (live output wrappers
+ utilities post Phase 5+ purge) to state%timecontrol via ASSOCIATE
tc_* aliases.

Spec: docs/superpowers/specs/2026-05-12-state-migration-timecontrol-design.md
EOF
)"
```

---

### Task 8: surfacewater.f90 + soilhydraulics.f90 reader cutover (34 + 33 sites)

**Files:**
- Modify: `src/drainage/surfacewater.f90`
- Modify: `src/soil/soilhydraulics.f90`

Same pattern. Two files in one commit (similar scope).

Commit message:
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-TC Task 8 — surfacewater + soilhydraulics reader cutover

surfacewater.f90 (~34) + soilhydraulics.f90 (~33) TC field reads
migrated to state%timecontrol via ASSOCIATE tc_* aliases.

Spec: docs/superpowers/specs/2026-05-12-state-migration-timecontrol-design.md
EOF
)"
```

---

### Task 9: readmeteo + meteoday + meteodt reader cutover (32 + 22 + 25 sites)

**Files:**
- Modify: `src/io/readmeteo.f90`
- Modify: `src/atmosphere/meteoday.f90`
- Modify: `src/atmosphere/meteodt.f90`

Atmosphere/meteo cluster. Note: meteodt has dtEventRain (cross-owned, out of scope per D10) — DON'T touch dtEventRain reads.

Commit message:
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-TC Task 9 — meteo cluster reader cutover

readmeteo (~32), meteoday (~22), meteodt (~25) TC field reads
migrated to state%timecontrol.

dtEventRain NOT touched (cross-owned with meteodt — out of scope
per D10).

Spec: docs/superpowers/specs/2026-05-12-state-migration-timecontrol-design.md D10
EOF
)"
```

---

### Task 10: cropgrowth.f90 reader cutover (31 sites)

**Files:**
- Modify: `src/crop/cropgrowth.f90`

The big crop file. Note: this file also reads `daycrop` and `yearcrop` (TC-owned per the inventory) — include those in the ASSOCIATE.

Commit message:
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-TC Task 10 — cropgrowth reader cutover (31 sites)

cropgrowth.f90 TC field reads (including daycrop/yearcrop) migrated
to state%timecontrol.

Spec: docs/superpowers/specs/2026-05-12-state-migration-timecontrol-design.md
EOF
)"
```

---

### Task 11: et + precipitation + interception + heat reader cutover

**Files:**
- Modify: `src/atmosphere/et.f90` (~19)
- Modify: `src/atmosphere/precipitation.f90` (~10)
- Modify: `src/atmosphere/interception.f90` (~5)
- Modify: `src/heat/temperature.f90` (~10)
- Modify: `src/heat/frozencond.f90` (~8)

Five small/medium files in one commit. Same pattern.

Commit message:
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-TC Task 11 — et/precipitation/interception/heat reader cutover

Atmosphere et/precipitation/interception + heat
temperature/frozencond TC field reads migrated to state%timecontrol.

Spec: docs/superpowers/specs/2026-05-12-state-migration-timecontrol-design.md
EOF
)"
```

---

### Task 12: solute + crop utilities + irrigation reader cutover

**Files:**
- Modify: `src/solute/solute.f90` (~6)
- Modify: `src/solute/agetracer.f90` (~3)
- Modify: `src/crop/rootextraction.f90` (~5)
- Modify: `src/crop/oxygenstress.f90` (~2)
- Modify: `src/crop/irrigation.f90` (~4)
- Modify: `src/crop/tillage.f90` (~2)
- Modify: `src/crop/management_soil.f90` (~3)
- Modify: `src/utils/soilhydraulicsutils.f90` (~5)
- Modify: `src/utils/surfacewaterutils.f90` (~4)

Long tail — all small files, group into one commit.

Commit message:
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-TC Task 12 — long-tail reader cutover

solute, agetracer, rootextraction, oxygenstress, irrigation, tillage,
management_soil, soilhydraulicsutils, surfacewaterutils TC field
reads migrated to state%timecontrol.

Spec: docs/superpowers/specs/2026-05-12-state-migration-timecontrol-design.md
EOF
)"
```

---

### Task 13: swap.f90 + initialize.f90 + swap_csv_output.f90 reader cutover

**Files:**
- Modify: `src/core/swap.f90`
- Modify: `src/core/initialize.f90` (~20)
- Modify: `src/io/swap_csv_output.f90` (~15)
- Modify: `src/io/toml/config_to_variables.f90`

Final reader files (the orchestrator + init + CSV output + adapter). swap.f90 has the main timestep loop and reads TC fields heavily. config_to_variables retains its pre-init writes for iyear/imonth/dt (Task 4 already handled).

Commit message:
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-TC Task 13 — core + init + CSV output reader cutover

swap.f90 (timestep loop), initialize.f90 (~20), swap_csv_output.f90
(~15), config_to_variables.f90 (remaining post-Task 4) TC field
reads migrated to state%timecontrol.

Spec: docs/superpowers/specs/2026-05-12-state-migration-timecontrol-design.md
EOF
)"
```

---

### Task 14: Strategy B compile-driven retirement

**Files:**
- Modify: `src/core/variables.f90` (comment out 61 globals)
- Modify: any files surfaced by compile errors (expected 5-10 hidden readers)
- Modify: `src/core/initialize.f90` (drop dead zero-inits)

- [ ] **Step 1: Baseline verify**

```bash
git status --short | head
pixi run check-full 2>&1 | tail -5
```

Expected: 5/5 byte-identical, clean working tree.

- [ ] **Step 2: Comment out 61 globals in variables.f90 with provenance**

For each TC-owned global declaration, prepend:
```fortran
! [SS-TC] retired 2026-05-12 — moved to state%timecontrol%X (ADR 0041)
```
and comment out the declaration line itself. Match boundary B-2.7 / soil-water S-2.12B style.

DO NOT touch:
- The 18 config-constants (tstart, tend, dtmin, dtmax, etc.)
- flZeroIntr / flZeroCumu
- dtEventRain
- Any non-TC global

- [ ] **Step 3: First compile pass**

```bash
pixi run check-full 2>&1 | tee /tmp/tc_p1.log | tail -50
```

Expected: build fails with errors only in TimeControl + reader files + initialize.f90. Errors are precise instructions: drop legacy write, migrate missed read, drop use-clause import.

- [ ] **Step 4: Iterate per error**

For each compile error:
- Legacy write line in TimeControl/IterTime → drop (state mirror is canonical).
- Missed reader → add `state%timecontrol%X` reference.
- `use variables, only: tc_*` → drop the retired names from only-list.
- initialize.f90 zero-init → drop.

Re-run `pixi run check-full` after each round. Expect 10-18 iterations.

- [ ] **Step 5: Final-state grep**

```bash
for v in t t1900 tcum dt dtmin0 dtprevious dtEvent tEvent tchange tnext tinit dtcrit daynr daycum iyear imonth daymonth monthstart monthend yearstart yearend daycrop yearcrop isteps MaxIterIterTime ioutdat cntper nprintcount nprintsecond nprintdayfirst outper idstart idend flDayStart flDayEnd flOutputShort flprintdt flprintshort flprint flprintsecond floutputshort flheader flheadirg flIrg1Start flIrg2Start flbaloutput flintrout flcumuout flBegOfDay flEndOfDay flLastTimeStep flFirstTimeStep flCropCalendar flCropEmergence flCropEnd flCropEvent flMet flMetDet flUpdMetDet; do
  hits=$(grep -rEn "\b$v\b" src/ --include="*.f90" \
    | grep -v "src/core/variables.f90" \
    | grep -v "src/state/" \
    | grep -v "%timecontrol%$v\b" \
    | grep -v "tc_$v\b\|_loc\|_tmp\|_old\|dummy_\|state_om%\|intent\|dimension" \
    | grep -v "^[^:]*:.*!" \
    | head -3)
  if [ -n "$hits" ]; then
    echo "=== $v ==="
    echo "$hits"
  fi
done
```

Expected: empty.

- [ ] **Step 6: Verify byte-identical**

```bash
pixi run check-full
pixi run -e test test-pfunit 2>&1 | tail -5
```

Expected: 5/5 byte-identical, 733/0/0.

- [ ] **Step 7: Commit**

```bash
git add -A
git commit -m "$(cat <<'EOF'
refactor(state): SS-TC Task 14 — Strategy B retire 61 timecontrol globals

Comments out 61 timecontrol-owned globals from variables.f90 with
[SS-TC] retired 2026-05-12 provenance markers. Compile-driven
discovery surfaced <N> remaining legacy references across <M>
files: <list>. Each fixed per the standard pattern.

Stays-legacy: 18 config-constants, flZeroIntr/flZeroCumu (future
reset-orchestration arc), dtEventRain (meteodt cross-ownership).

Strategy B compile passes: <N> (4th application of the pattern
established in soilwater-core S-2.12B).

check-full 5/5 byte-identical. Hupselbrook canary <T>s.

Spec: docs/superpowers/specs/2026-05-12-state-migration-timecontrol-design.md
EOF
)"
```

---

## Phase 3 — ADR + merge

### Task 15: ADR 0041 + playbook + merge prep

**Files:**
- Create: `docs/adr/0041-state-migration-timecontrol.md`
- Modify: `docs/adr/index.md`
- Modify: `docs/superpowers/specs/state-migration-playbook.md`

- [ ] **Step 1: Write ADR 0041**

Match ADR 0038/0039 structure (~120-150 lines). Sections: front-matter, Decision, Context, Decisions list (D1-D13 from design), Out of scope, Consequences (file counts, LoC, hupselbrook timing), Cross-references.

Emphasize:
- Largest reader inventory of any arc (~415 sites across 22 files post-purge; originally projected ~500/27 pre-purge).
- 4th application of Strategy B compile-driven retirement.
- Pre-init transient-buffer pattern (atmosphere A-2.6 lesson applied for iyear/imonth/dt).
- TC OWNS the flzero gates rather than consumes them — future "reset-orchestration arc" will consolidate the cross-subsystem consumers (deferred per D10).

- [ ] **Step 2: Update ADR index**

Append the ADR 0041 row to `docs/adr/index.md`.

- [ ] **Step 3: Update playbook**

Append 1-2 new lessons to `docs/superpowers/specs/state-migration-playbook.md` if any surfaced. Candidates:
- **Reader fanout dominates over field count.** TC has 61 fields but ~415 reads — the migration cost is determined by reader inventory, not field count. Arcs should be sized by reader sites, not owned globals.
- **"Reset gate owner" vs "reset consumer" boundary.** When a subsystem OWNS a reset gate (writes flzeroIntr/flzeroCumu), migrating the gate forces all consumers to be in scope. Defer to a dedicated cross-subsystem arc.

- [ ] **Step 4: Verify byte-identical**

```bash
pixi run check-full 2>&1 | tail -3
```

Expected: 5/5 byte-identical (docs-only changes shouldn't affect compile).

- [ ] **Step 5: Commit**

```bash
git add docs/adr/0041-state-migration-timecontrol.md docs/adr/index.md docs/superpowers/specs/state-migration-playbook.md
git commit -m "$(cat <<'EOF'
docs(adr): ADR 0041 — TimeControl state migration

Documents migration #10 (TimeControl) — largest reader inventory of
any arc to date (~415 sites across 22 files; was projected ~500/27
before the ADR 0009 Phase 5+ purge and macropore retirement).

61 owned globals retired (16 clock + 19 schedule + 26 boolean).
18 config-constants deferred. flZeroIntr/flZeroCumu reset gates
deferred to future reset-orchestration arc. dtEventRain
cross-ownership kept as legacy.

4th application of Strategy B compile-driven retirement (after
soilwater-core S-2.12B, atmosphere A-2.6, tillage T-5). TimeControl
intent bumped to inout. IterTime gained state arg.

Playbook: <N> new lessons captured.

Spec: docs/superpowers/specs/2026-05-12-state-migration-timecontrol-design.md
EOF
)"
```

- [ ] **Step 6: Merge to development**

```bash
git checkout development
git merge --ff-only refactor/timecontrol-state
git branch -d refactor/timecontrol-state
git log --oneline -3
```

---

## Self-review

**Spec coverage check:**
- D1 (scope, 61 fields, defer 18+2) → Tasks 1, 14
- D2 (flat layout) → Task 1
- D3 (aggregation) → Task 1
- D4 (no init routine) → Task 4 implements inline init via TC(case=1)
- D5 (TimeControl intent bump) → Task 2
- D6 (IterTime plumbing) → Task 2
- D7 (pre-init transient pattern) → Task 4
- D8 (DLL re-init state write) → Task 4
- D9 (swapoutput co-writes) → Task 5
- D10 (defer flZero* + dtEventRain) → covered in "Out of scope" + non-touches in Tasks 9, 14
- D11 (Strategy B retirement) → Task 14
- D12 (pre-flight dual-write check) → Standing instructions + Tasks 6-13 each
- D13 (waterbalance heaviest reader) → Task 6 dedicated
- D14 (monolithic arc) → 15 tasks, single branch

**Placeholder scan:** No TBDs, no "implement later" — every step has concrete code/commands.

**Type consistency:** Field names from refreshed discovery used throughout. ASSOCIATE prefix `tc_` consistent across all tasks. The 61-field count is approximate; Task 1 Step 1 verifies the exact set on current HEAD.
