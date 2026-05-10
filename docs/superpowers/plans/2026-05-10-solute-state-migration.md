# Solute State Migration — Implementation Plan (all phases)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Migrate the solute subsystem off `variables.f90` globals into a typed `solute_state_t` aggregated under `swap_state_t`, extract `AgeTracer` to its own module with a stub-error (dead-code preservation), and patch the physics config gap by promoting 14 currently-missing fields into `solute_config_t` (preventing silent zero-defaults on TOML runs with `swsolu=1`).

**Architecture:** Three logical phases bundled into one plan:
- **Phase 0** patches `solute_config_t` for the 14 missing physics fields. Mechanical schema extension; no regression coverage (none of the 5 cases activates `swsolu=1`).
- **Phase 1** defines `solute_state_t`, extracts `AgeTracer` to a stub-errored module, threads state, dual-write transitional pattern.
- **Phase 2** migrates cross-subsystem readers, drops dual-write, removes globals, finalizes ADR.

**Tech Stack:** Fortran 2008, gfortran, meson + ninja + pixi, pFUnit, check-full byte-identical regression.

**Spec:** `docs/superpowers/specs/2026-05-10-state-migration-solute-design.md`
**Discovery:** `docs/superpowers/specs/2026-05-10-state-migration-solute-discovery.md`

**Branch:** all commits go on `refactor/solute-state` (new branch from development; ready to merge after all tasks land green).

---

## Lessons-learned applied (from surfacewater + drainage)

- **Section 3.5 categorization** by intent (output / compute / working-buffer / init-seed / call-site argument). Discovery already populated.
- **Verification grep BEFORE deletion.** Run both `use variables.*\b<var>\b` and raw-symbol grep before removing declarations from `variables.f90`. Any non-zero count is an unmigrated reader.
- **`sw_*`/`dr_*` ASSOCIATE prefix** when bare names would shadow `use Variables`. Use `sl_*` (solute) here.
- **Integration-gate task usually exceeds planned scope.** Plan accordingly — solute Phase 2 final cleanup may catch hidden gaps (working-buffer reads, init-seed reads, etc.).

---

## File Structure

**Created:**
- `src/state/solute_state.f90` — `solute_state_t` definition (Phase 1)
- `src/solute/agetracer.f90` — extracted AgeTracer subroutine with stub-error (Phase 1)
- `tests/unit/state/test_solute_state.pf` — pFUnit tests (Phase 1)
- `tests/unit/config/test_solute_config_phase0_fields.pf` — pFUnit tests for the 14 promoted config fields (Phase 0)
- `docs/adr/0032-state-migration-solute.md` — finalized at end (Phase 2)

**Modified:**
- `src/config/solute_config.f90` — add 14 fields + validators (Phase 0)
- `src/io/toml/read_solute_toml.f90` — read 14 fields (Phase 0)
- `src/io/toml/config_to_variables.f90` — `apply_solute` populates the 14 legacy globals from typed config (Phase 0)
- `src/state/swap_state.f90` — add `solute` field (Phase 1)
- `src/solute/solute.f90` — remove AgeTracer subroutine (Phase 1); ASSOCIATE on `state%solute` in compute body; drop dual-write (Phase 2)
- `src/io/swapoutput.f90` — `outsba`, `outend`, `outvap`, `outbal`, `outage` switch to `state%solute`; `outage` gated on `flAgeTracer` (always false; effectively dead) (Phase 1, Phase 2)
- `src/io/swap_csv_output.f90` — `set_values` reads `state%solute%cml(:)` etc. (Phase 1, Phase 2)
- `src/crop/rootextraction.f90` — `cml` read switches to `state%solute%cml` (Phase 2)
- `src/crop/irrigation.f90` — same (Phase 2)
- `src/io/toml/config_to_variables.f90` — `cml(k)` init seed switches to `state%solute%cml(k)` (Phase 2)
- `src/core/variables.f90` — comment out 25 solute-owned declarations with provenance markers (Phase 2)
- `src/core/initialize.f90` — drop the corresponding zero-init lines (Phase 2)
- `tests/unit/testSuites.inc` — add new suites
- meson source list — add `src/state/solute_state.f90` and `src/solute/agetracer.f90`

**No changes (intentionally) to:**
- `ArMpSs` (shared working buffer co-written by 3 subsystems; macropore migration sorts ownership later)
- The 12 AgeTracer-specific globals in `variables.f90` (kept declared but unwritten; future agetracer migration removes them)

---

# Phase 0 — Physics config gap

The 14 missing fields per discovery hazard #6, all confirmed declared in `variables.f90` (lines 1042-1077 area):

| Field | variables.f90 | Type | Notes |
|---|---|---|---|
| `swbr` | 1042 | integer | breakthrough switch |
| `cpre` | 1053 | real | precipitation conc |
| `cref` | 1054 | real | Freundlich reference conc |
| `cseeptab(mabbc*2)` | 1056 | real array | seepage table — CSV companion candidate |
| `daquif` | 1058 | real | aquifer thickness |
| `ddif` | 1059 | real | molecular diffusion |
| `decpot(maho)` | 1060 | real per-layer | potential decomposition |
| `decsat` | 1061 | real | aquifer decomposition |
| `fdepth(maho)` | 1065 | real per-layer | depth correction |
| `frexp` | 1066 | real | Freundlich exponent |
| `gampar` | 1067 | real | low-temp reduction |
| `kf(maho)` | 1074 | real per-layer | Freundlich coefficient per layer |
| `kfsat` | 1075 | real | aquifer linear adsorption |
| `poros` | 1077 | real | aquifer porosity |

3 are per-layer (`decpot`, `fdepth`, `kf`) — sized `maho` (max horizons). 1 is a 2-column table (`cseeptab` — `mabbc*2` flat, 2 columns logically). The rest are scalars.

### Task 1: Phase 0 — Promote 14 physics fields into solute_config_t

**Files:**
- Modify: `src/config/solute_config.f90` — add 14 fields + validators
- Modify: `src/io/toml/read_solute_toml.f90` — read 14 fields
- Modify: `src/io/toml/config_to_variables.f90` — populate legacy globals in `apply_solute`
- Create: `tests/unit/config/test_solute_phase0_fields.pf` — validator tests
- Modify: `tests/unit/testSuites.inc`

- [ ] **Step 1: Inspect existing solute_config_t shape**

```bash
sed -n '1,50p' src/config/solute_config.f90
grep -n "subroutine\|public\|type ::" src/config/solute_config.f90 | head
```

Note the existing field declarations, the validator structure, and the type's `procedure :: validate` etc.

- [ ] **Step 2: Inspect the existing TOML reader**

```bash
sed -n '1,50p' src/io/toml/read_solute_toml.f90
```

Note the existing read pattern (`get_optional_real_with_default`, `get_optional_int_with_default`, etc.).

- [ ] **Step 3: Inspect the existing apply_solute adapter**

```bash
grep -n "apply_solute\|swsolu\|cdrain\|cseep" src/io/toml/config_to_variables.f90 | head -20
```

Find where the existing solute fields are populated into legacy globals. Add the 14 new fields nearby.

- [ ] **Step 4: Write failing pFUnit tests for the new validators**

Create `tests/unit/config/test_solute_phase0_fields.pf`. For each scalar field (cref, daquif, ddif, decsat, frexp, gampar, kfsat, poros, swbr, cpre): test default value, valid range, out-of-range rejection.

For per-layer arrays (`decpot`, `fdepth`, `kf`): test required-when-flag, length-matches-maho, range per element.

For `cseeptab` table: test allocated, expected dimensions, range.

Example:

```fortran
@test
subroutine test_solute_config_cref_default_zero()
   use funit
   use iso_fortran_env, only: real64
   use solute_config_mod, only: solute_config_t
   type(solute_config_t) :: cfg
   @assertEqual(0.0_real64, cfg%cref, 1.0e-12_real64)
end subroutine

@test
subroutine test_solute_config_cref_negative_rejected()
   use funit
   use iso_fortran_env, only: real64
   use solute_config_mod, only: solute_config_t
   use error_mod, only: error_collection_t
   type(solute_config_t) :: cfg
   type(error_collection_t) :: errors

   cfg%swsolu = 1
   cfg%cref   = -1.0_real64
   call cfg%validate(errors)
   @assertGreaterThan(0, errors%count())   ! at least one error
end subroutine
```

Add `ADD_TEST_SUITE(test_solute_phase0_fields_suite)` to `tests/unit/testSuites.inc`.

- [ ] **Step 5: Run tests — confirm FAIL**

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: failures referencing missing fields on solute_config_t.

- [ ] **Step 6: Add the 14 fields to solute_config_t**

For each scalar field, add a typed declaration with default:
```fortran
real(real64) :: cref    = 0.0_real64
real(real64) :: daquif  = 0.0_real64
real(real64) :: ddif    = 0.0_real64
real(real64) :: decsat  = 0.0_real64
real(real64) :: frexp   = 0.0_real64
real(real64) :: gampar  = 0.0_real64
real(real64) :: kfsat   = 0.0_real64
real(real64) :: poros   = 0.0_real64
real(real64) :: cpre    = 0.0_real64
integer      :: swbr    = 0
```

For per-layer arrays (sized by `maho`, the max horizons constant from the existing layered config — find via grep):
```fortran
real(real64), allocatable :: decpot(:)   !! per soil layer
real(real64), allocatable :: fdepth(:)   !! per soil layer
real(real64), allocatable :: kf(:)       !! per soil layer
```

For the seepage table `cseeptab(mabbc*2)` — flat in legacy. In typed config, prefer a 2-column 2D array:
```fortran
real(real64), allocatable :: cseeptab(:,:)   !! (n, 2) — col 1 = time, col 2 = concentration
```

Add validators in `solute_config_validate`. For scalars: `check_real_range(cfg%cref, 0.0_real64, 1.0e6_real64, "solute.cref", errors)`. For per-layer arrays: check allocated when needed (e.g., when `swsolu=1`), check size equals `maho_actual` (the active number of horizons).

- [ ] **Step 7: Add reads in read_solute_toml.f90**

For each scalar:
```fortran
call get_optional_real_with_default(sec, "cref", config%cref, 0.0_real64, "solute.cref", errors)
```

For per-layer arrays:
```fortran
call get_optional_real_array(sec, "kf", config%kf, "solute.kf", errors)
```

For `cseeptab`: use the established CSV-companion pattern (ADR 0012) if the table is large; inline TOML if small. Inspect what other tabular configs do. Recommend: inline 2D array `[[t1, c1], [t2, c2], …]` since seepage tables are typically <100 rows.

- [ ] **Step 8: Update apply_solute adapter**

In `src/io/toml/config_to_variables.f90`, find `apply_solute` (or wherever the existing solute fields like `cdrain`, `cseep` are written). Add lines populating the 14 legacy globals:

```fortran
cref   = config%solute%cref
daquif = config%solute%daquif
ddif   = config%solute%ddif
decsat = config%solute%decsat
frexp  = config%solute%frexp
gampar = config%solute%gampar
kfsat  = config%solute%kfsat
poros  = config%solute%poros
cpre   = config%solute%cpre
swbr   = config%solute%swbr

if (allocated(config%solute%decpot))  decpot(1:size(config%solute%decpot))  = config%solute%decpot
if (allocated(config%solute%fdepth))  fdepth(1:size(config%solute%fdepth))  = config%solute%fdepth
if (allocated(config%solute%kf))      kf(1:size(config%solute%kf))          = config%solute%kf

if (allocated(config%solute%cseeptab)) then
   ! flatten 2D table to 1D for legacy global storage (col1, col1, … col2, col2, …) — verify the legacy layout
   cseeptab(:) = ...
end if
```

The `cseeptab` flattening needs care — read the legacy reader (`readsolu` or wherever it parses) to confirm the 1D-vs-2D layout, OR read the consumers in `solute.f90` to see the indexing pattern. If unclear, raise.

- [ ] **Step 9: Run tests — confirm PASS**

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: zero failures, count = previous baseline + however many new tests Step 4 added.

Run: `pixi run check-full`
Expected: 5/5 byte-identical (none of the 5 cases activates swsolu, so the new fields default to zero — same as before).

- [ ] **Step 10: Commit**

```bash
git add src/config/solute_config.f90 src/io/toml/read_solute_toml.f90 \
        src/io/toml/config_to_variables.f90 \
        tests/unit/config/test_solute_phase0_fields.pf tests/unit/testSuites.inc
git commit -m "$(cat <<'EOF'
feat(config): SS-SLST Phase 0 — promote 14 missing solute physics fields

solute_config_t gains: cref, kf, frexp, gampar, decpot, fdepth, ddif,
kfsat, poros, daquif, decsat, cseeptab, swbr, cpre. Each with
typed default + validator + TOML reader + adapter wiring.

Patches a real correctness gap: TOML runs with swsolu=1 previously
silently used zeroed defaults for core Freundlich/decomposition/
breakthrough physics. None of the 5 regression cases activates
solute, so check-full byte-identical is preserved (defaults are
still zero), but future TOML runs with swsolu=1 will now correctly
honor the user-supplied values.

Spec: docs/superpowers/specs/2026-05-10-state-migration-solute-design.md
EOF
)"
```

---

# Phase 1 — State type, AgeTracer extraction, threading, dual-write + drop

### Task 2: Phase 1 — Create solute_state_t

**Files:**
- Create: `src/state/solute_state.f90`
- Create: `tests/unit/state/test_solute_state.pf`
- Modify: meson + testSuites.inc

- [ ] **Step 1: Identify exact state field set from discovery Section 2**

Read discovery doc `docs/superpowers/specs/2026-05-10-state-migration-solute-discovery.md` Section 2. The 37 owned globals minus the 12 AgeTracer-specific = ~25 fields for `solute_state_t`.

Likely set (verify against discovery):
- Per-node (allocatable, sized `numnod`): `cml(:)`, `cmsy(:)`, `cnh4(:)`, `cno3(:)` (and any other per-node solute concentration arrays)
- Per-soil-layer (allocatable, sized `numlay`/`maho`): possibly some decomposition state arrays
- Cumulative balance scalars: `samini`, `sampro`, `samcra`, `solbal`, `dectot`, `rottot`
- Cumulative source/sink scalars: `sqprec`, `sqirrig`, `sqdra`, `sqbot`, `sqrap` (if solute owns these — verify via discovery)
- Per-step intermediates: any `q*` scalars discovery flagged

Some scalars currently classified as "owned" by discovery may actually be water-balance-owned (e.g., `sqdra` is also referenced in waterbalance.f90 in surfacewater/drainage migrations). Reconfirm: a global is solute-owned if SOLUTE WRITES IT — even if other subsystems write it too, but in solute's case the typical pattern is SOLUTE_ONLY for solute-mass accumulators and SHARED for water-flux accumulators. Be discriminating.

- [ ] **Step 2: Write failing tests + implement module**

Same shape as drainage Phase 1 Task 1 / surfacewater Phase 1 Task 1. Default scalars at 0; allocatables initially unallocated; allocation lifecycle test.

- [ ] **Step 3: Add to swap_state_t**

In `src/state/swap_state.f90`:
```fortran
use solute_state_mod, only: solute_state_t
...
type(swap_state_t) ::
   type(surfacewater_state_t) :: surfacewater
   type(drainage_state_t)     :: drainage
   type(solute_state_t)       :: solute    ! NEW
end type
```

- [ ] **Step 4: Build, test, commit**

```bash
pixi run -e test test-pfunit
git commit -m "feat(state): SS-SLST Task 2 — solute_state_t typed state record"
```

### Task 3: Phase 1 — Extract AgeTracer to its own module with stub-error

**Files:**
- Create: `src/solute/agetracer.f90`
- Modify: `src/solute/solute.f90` — remove AgeTracer subroutine
- Modify: `src/core/swap.f90` — caller of AgeTracer (if any) imports from new module
- Modify: `src/io/swapoutput.f90` — outage routine gated on flAgeTracer (always false today)
- Modify: meson source list

- [ ] **Step 1: Locate AgeTracer**

```bash
grep -n "subroutine AgeTracer\|call AgeTracer" src/solute/solute.f90 src/core/swap.f90 src/io/swapoutput.f90
```

Confirm AgeTracer at line 304 of solute.f90 per discovery. Find every caller.

- [ ] **Step 2: Create the new module**

Create `src/solute/agetracer.f90`:

```fortran
!> @file agetracer.f90
!! AgeTracer feature — extracted from solute.f90 during ADR 0032
!! to clean the solute module's compute path.
!!
!! Status: DEAD. flAgeTracer is never set to .true. anywhere in
!! the codebase. The runtime body is preserved for source-completeness
!! but is gated behind a stub-error.
!!
!! To re-enable AgeTracer:
!!   1. Audit the 12 AgeTracer-specific globals currently declared
!!      in variables.f90 (cml is dual-used; the other 11 are dormant).
!!   2. Decide: separate agetracer_state_t, or share state%solute%cml
!!      with documented sequencing rules?
!!   3. Wire flAgeTracer from typed config or external trigger.
!!   4. Remove the stub-error guard in this file.
!!
!! See discovery doc Section 8 hazard #4 (cml dual-use).
module agetracer_mod
   use error_mod, only: fatalerr_collected
   use swap_state_mod, only: swap_state_t
   implicit none
   private
   public :: AgeTracer

contains

subroutine AgeTracer(task, state)
   use Variables, only: flAgeTracer
   integer,            intent(in)    :: task
   type(swap_state_t), intent(inout) :: state

   if (flAgeTracer) then
      call fatalerr_collected('AgeTracer', &
         'AgeTracer feature is currently inert. The runtime body was '// &
         'preserved during ADR 0032 (solute migration) for future '// &
         'reactivation. To re-enable: see agetracer.f90 file comment.')
      return
   end if

   ! Body kept for source completeness — never reached in current build.
   ! [Move the original AgeTracer body here verbatim from solute.f90:304+]

end subroutine AgeTracer

end module agetracer_mod
```

The body of AgeTracer (the original ~244 lines from solute.f90 lines 304-548) is moved verbatim. The `if (flAgeTracer)` guard is added at the top. The body keeps its existing `use Variables` clause for the 12 AgeTracer-specific globals; those stay declared in variables.f90 (the comment marker says "future agetracer_state_t target").

- [ ] **Step 3: Remove AgeTracer from solute.f90**

In `src/solute/solute.f90`, delete lines 304-548 (AgeTracer subroutine and any private helpers it owns). Update the module's `public ::` line to remove `AgeTracer`.

- [ ] **Step 4: Update callers**

Wherever `call AgeTracer(...)` appears, add `use agetracer_mod, only: AgeTracer` to the caller and keep the call site unchanged. Most likely caller: `src/core/swap.f90`.

- [ ] **Step 5: Gate outage in swapoutput.f90**

Find `subroutine outage` (or the AgeTracer output routine). Wrap the body in `if (flAgeTracer) then ... end if`. Since flAgeTracer is always false, this makes outage a no-op — but the code stays for the future reactivation.

- [ ] **Step 6: Add the new file to meson**

Add `src/solute/agetracer.f90` to the source lists in `meson.build` and `tests/unit/meson.build` (if it needs to be visible to tests).

- [ ] **Step 7: Build and test**

```bash
pixi run check-full
pixi run -e test test-pfunit 2>&1 | tail -10
```

Expected: 5/5 byte-identical (AgeTracer was never reached anyway), zero pFUnit failures.

- [ ] **Step 8: Commit**

```bash
git add src/solute/agetracer.f90 src/solute/solute.f90 src/core/swap.f90 \
        src/io/swapoutput.f90 meson.build tests/unit/meson.build
git commit -m "refactor(solute): SS-SLST Task 3 — extract AgeTracer to its own module"
```

### Task 4: Phase 1 — Solute compute writes state alongside globals (dual-write)

Solute's owned globals get parallel `state%solute%*` writes. Compute body still writes legacy globals — readers haven't migrated yet.

**Files:**
- Modify: `src/solute/solute.f90` — add `state%solute%X = ...` writes paralleling each `X = ...` global write

- [ ] **Step 1: Find write sites**

```bash
grep -nE "^\s+(cml|cmsy|cnh4|cno3|sampro|samcra|solbal|dectot|rottot|sqprec|sqirrig|sqdra|sqbot)\s*[=\(]" src/solute/solute.f90 | head -40
```

Per discovery Section 2, identify each write of a solute-owned global. Catalog file:line pairs.

- [ ] **Step 2: Add dual-writes**

Apply ASSOCIATE pattern with `sl_*` prefix where `use Variables` shadows:

```fortran
associate(sl_cml => state%solute%cml, sl_cmsy => state%solute%cmsy, &
          sl_sampro => state%solute%sampro, ...)
   ! body — wherever cml(i) = ... write happens, also write sl_cml(i) = cml(i)
end associate
```

Or for each global write, immediately after, add:
```fortran
state%solute%cml(i) = cml(i)
```

For per-node arrays, slice copies after a loop are cleaner:
```fortran
do i = 1, numnod
   cml(i) = ...
end do
state%solute%cml(:) = cml(:)
```

- [ ] **Step 3: Build, verify check-full + pFUnit, commit**

```bash
pixi run check-full && pixi run -e test test-pfunit 2>&1 | tail
git add src/solute/solute.f90
git commit -m "refactor(state): SS-SLST Task 4 — solute compute dual-writes state and globals"
```

### Task 5: Phase 1 — Output reads switch to state%solute

**Files:**
- Modify: `src/io/swapoutput.f90` (outsba, outend, outvap, outbal, possibly others)
- Modify: `src/io/swap_csv_output.f90` (set_values per-node solute reads)

- [ ] **Step 1: Find output reads of solute-owned globals**

```bash
grep -rEn "use variables.*\b(cml|cmsy|cnh4|cno3|sampro|samcra|solbal|dectot|rottot)\b" src/io/swapoutput.f90 src/io/swap_csv_output.f90 | head
```

- [ ] **Step 2: Migrate each reader**

For each output routine that reads solute globals, add `state` to its argument list (most output routines already take `state` from surfacewater Phase 2). Switch reads to `state%solute%*`. Drop migrated symbols from `use variables, only:` clauses.

- [ ] **Step 3: Build, verify, commit**

```bash
pixi run check-full && pixi run -e test test-pfunit 2>&1 | tail
git add src/io/swapoutput.f90 src/io/swap_csv_output.f90 src/core/swap.f90
git commit -m "refactor(state): SS-SLST Task 5 — output reads solute from state"
```

### Task 6: Phase 1 — Drop dual-write (state authoritative for solute home tree)

**Files:**
- Modify: `src/solute/solute.f90` — drop legacy global writes; keep state writes only

- [ ] **Step 1: Drop the legacy global writes**

For each dual-write pair from Task 4, delete the legacy global write line. Keep the `state%solute%*` write.

- [ ] **Step 2: Trim use variables clauses**

For each subroutine in solute.f90, drop the migrated symbol names from `use variables, only:` clauses where no longer referenced.

- [ ] **Step 3: Verify ALL solute-owned globals have zero remaining readers/writers in solute home tree**

```bash
grep -rEn "use variables.*\b(cml|cmsy|cnh4|cno3|sampro|samcra|solbal|dectot|rottot|sqprec|sqirrig|sqdra|sqbot)\b" src/solute/ --include="*.f90"
```

Expected: empty (or minimal — ASSOCIATE-bound names are not in the use clause anyway).

- [ ] **Step 4: Build, verify check-full + pFUnit, commit**

This is the Phase 1 integration gate. **Cross-subsystem readers (rootextraction, irrigation, output routines) still read globals.** Phase 2 migrates them. For Phase 1 to ship green, the dual-write pattern must keep the globals current via SOMETHING. Wait — at this commit, solute compute STOPS writing globals. The legacy globals will go stale. Cross-subsystem readers will see stale data. check-full could fail.

**Decision:** Task 6 (drop solute compute's global writes) is held until AFTER Phase 2 Task 7 (compute readers migrate). Re-order: Phase 1 ends after Task 5 (output reads migrated). Phase 2 starts with compute reader migration THEN drops dual-write THEN deletes globals.

**Revised ordering:** drop Task 6 from Phase 1; renumber. The dual-write stays through all of Phase 2's reader migrations. Phase 2's "drop dual-write" happens after the last reader migrates.

```bash
# (Task 6 from this plan is REMOVED. Phase 1 ends at Task 5. Phase 2 begins.)
```

---

# Phase 2 — Cross-subsystem reader migration + global cleanup + ADR

### Task 6 (was 7): Phase 2 — Migrate compute readers (rootextraction, irrigation)

**Files:**
- Modify: `src/crop/rootextraction.f90` — `cml` reads → `state%solute%cml`
- Modify: `src/crop/irrigation.f90` — `cml` reads → `state%solute%cml`

- [ ] **Step 1: Find reads**

```bash
grep -nE "\bcml\b" src/crop/rootextraction.f90 src/crop/irrigation.f90 | head
```

- [ ] **Step 2: Migrate**

Add `state` to the consuming subroutines if not already (likely already plumbed). Switch reads to `state%solute%cml`. Drop `cml` from `use variables, only:` if no longer referenced.

- [ ] **Step 3: Verify, commit**

```bash
pixi run check-full && pixi run -e test test-pfunit 2>&1 | tail
git commit -m "refactor(state): SS-SLST Task 6 — rootextraction/irrigation read solute from state"
```

### Task 7 (was 8): Phase 2 — Migrate config-time cml seed

**Files:**
- Modify: `src/io/toml/config_to_variables.f90` — the `cml(k)` init seed

- [ ] **Step 1: Find the seed**

```bash
grep -n "cml" src/io/toml/config_to_variables.f90
```

The seed is likely in `apply_solute` or similar — it copies a CSV-derived array into the legacy `cml(k)` global at config-load time.

- [ ] **Step 2: Decide: write to global, write to state, or both**

Two-stage `cml` seeding (discovery hazard #7): config-time seed → runtime `afgen` re-interpolation.

- If config_to_variables runs BEFORE drainage_init/solute_init (which allocate state), then writing to `state%solute%cml` here would fail (not allocated). Keep writing to global.
- If we add a `solute_init(state)` analogous to drainage_init, called after config_to_variables but before solute(task=1), then solute_init can copy from global to state.

Recommend the latter: introduce `solute_init(state)` that copies the cml seed from global to state%solute%cml, called after `config_to_variables` and before `solute(task=1)`.

- [ ] **Step 3: Add solute_init**

In `src/solute/solute.f90` (or a new init file), add:
```fortran
subroutine solute_init(state)
   use Variables, only: cml, numnod
   type(swap_state_t), intent(inout) :: state
   integer :: i

   if (.not. allocated(state%solute%cml))   allocate(state%solute%cml(numnod))
   ! ... allocate other per-node arrays ...
   state%solute%cml(:) = cml(:)   ! copy config-seeded values
end subroutine
```

- [ ] **Step 4: Wire from swap_main**

```fortran
call config_to_variables(config)
...
call drainage_init(state)
call solute_init(state)        ! NEW
call SurfaceWater(1, state, ...)
...
```

- [ ] **Step 5: Verify, commit**

### Task 8 (was 9): Phase 2 — Drop solute compute's dual-write

After Tasks 6 + 7 migrate all consumers, compute can stop writing legacy globals.

**Files:**
- Modify: `src/solute/solute.f90` — remove legacy global writes from compute

- [ ] Drop the legacy global writes (mirror of original Task 6).
- [ ] Verify check-full byte-identical (integration gate).
- [ ] Commit.

### Task 9 (was 10): Phase 2 — Comment out 25 solute-owned globals

**Files:**
- Modify: `src/core/variables.f90` — comment out 25 declarations with provenance markers
- Modify: `src/core/initialize.f90` — drop matching zero-init lines

- [ ] **Step 1: Verify zero readers via the 5-category grep**

Per drainage's Section 3.5 addendum (lessons-learned), run BOTH greps:

```bash
# Category 1+2 (use variables clauses):
for v in cml cmsy cnh4 cno3 sampro samcra solbal dectot rottot sqprec sqirrig sqdra sqbot dwatlay; do
  grep -rEn "use variables.*\b$v\b" src/ --include="*.f90" | grep -v "src/core/variables.f90 ; src/state/ ; src/solute/ ; src/io/toml/" | head
done

# Category 3+4+5 (raw symbol references after migration):
for v in cml cmsy cnh4 cno3 sampro samcra solbal dectot rottot sqprec sqirrig sqdra sqbot; do
  grep -rEn "\b$v\b" src/ --include="*.f90" \
    | grep -v "src/core/variables.f90" \
    | grep -v "src/core/initialize.f90" \
    | grep -v "src/state/" \
    | grep -v "state%" \
    | grep -v "config%" | head
done
```

Both must return zero hits before declarations are commented out. If any non-zero, those are unmigrated — handle in this task.

- [ ] **Step 2: Comment out declarations**

Match the surface-water/drainage pattern:
```fortran
! real(8) cml(macp)        ! Moved to solute_state_t%cml (ADR 0032)
! real(8) cmsy(macp)       ! Moved to solute_state_t%cmsy (ADR 0032)
! ...
```

- [ ] **Step 3: Drop initialize.f90 lines**

Same pattern.

- [ ] **Step 4: Build, verify, commit**

This is the final integration gate.

```bash
pixi run check-full && pixi run -e test test-pfunit 2>&1 | tail
git commit -m "refactor(state): SS-SLST Task 9 — comment out 25 solute globals"
```

### Task 10 (was 11): Phase 2 — ADR 0032 finalize + index + verification

**Files:**
- Create: `docs/adr/0032-state-migration-solute.md`
- Modify: `docs/adr/index.md`

- [ ] Author ADR 0032 matching the structure of 0030/0031.
- [ ] Add index entry in house style.
- [ ] Final verification: pFUnit zero failures, check-full 5/5.
- [ ] Commit.

---

## Self-Review Notes

- **Spec coverage:** D1 (scope) → all tasks. D2 (state shape) → Task 2. D3 (aggregator) → Task 2. D4 (threading) → already in place from surfacewater Phase 2. D5 (AgeTracer extract) → Task 3. D6 (Phase 0 physics gap) → Task 1. D7 (cross-subsystem readers) → Tasks 6, 7. D8 (outage gating) → Task 3. D9 (ArMpSs stays global) → no task; documented. D10 (cpre/cseeptab in config) → Task 1.
- **Lessons-learned applied:** Pre-flight inventory upfront via discovery Section 3.5; verification grep before deletion (Task 9 Step 1); ASSOCIATE `sl_*` prefix for shadow-safe aliasing; integration-gate task (Task 9) plans for surprise scope.
- **Phase 1 + Phase 2 boundary** is between Task 5 and Task 6. Plan-as-written re-ordered Task 6 (drop dual-write) to Task 8 (after readers migrate) — this is the cleanest invariant: dual-write stays until its last reader is gone.
- **AgeTracer body preservation:** Task 3's stub-errored module keeps the original 244 lines for future reactivation. The 12 AgeTracer-specific globals stay declared (commented as future targets) but are never written.
- **Branch policy:** All commits on `refactor/solute-state`. After Task 10 verifies green, branch is ready for review and merge to `development`.
