# SS-5 — readmeteo.f90 TTutil branch deletion (ADR 0014 Steps 2–3)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Retire the TTutil-based meteorology readers from `readmeteo.f90` (per-year `.YYY` daily reader, `.met` all-years reader, TTutil sub-daily detail reader, and TTutil tail of `ReadRainEvents`) and the matching `swMetFilAll` pre-load wiring in `config_to_variables.f90`. Reject any non-CSV `metfile` at the TOML validation boundary first to prevent the deleted code paths from being requested. Sweep the now-dead variables from `variables.f90`.

**Architecture:** Three sequential commits, each independently shippable with full check-full + pFUnit green:

1. **Stub-error commit.** Add a stub-error block to `meteorology_config_validate` that rejects `metfile` not ending in `.csv` with `ERR_VALIDATION_CROSS_FIELD`. Mirror the SS-4 / cropfixed_config pattern. Also create the audit doc (`docs/phase-4f-readmeteo-ttutil-audit.md`) and update ADR 0014 marking the prerequisite met.
2. **Deletion commit.** Delete the four TTutil branches in `readmeteo.f90` (lines ~106–144 daily/detail TTutil + ~442–593 `MeteoInOneFile` + ~367–end TTutil tail of `ReadRainEvents`) and the `swMetFilAll == 1` block in `config_to_variables.f90` (lines 176–200). Simplify control flow: drop `goto 100` / `100 continue`, drop the `if (swMetCSV == 1)` and `if (swRainCSV == 1)` guards (only one path remains). Mark ADR 0014 Step 2 done.
3. **Dead-variable sweep commit.** Delete `swMetFilAll`, `swMetCSV`, `swRainCSV`, `swMetDetCSV`, `rainfil`, `station(366)`, integer `ad`/`am` arrays from `variables.f90` per ADR 0014 Step 3. Flip umbrella-spec SS-5 row to DONE. Mark ADR 0014 Step 3 done.

**Tech Stack:** Fortran 2008 (gfortran 13), pFUnit, pixi, meson. `error_collection_t` + `ERR_VALIDATION_CROSS_FIELD` already imported in `meteorology_config.f90` (line 4).

**Sub-spec parent:** `docs/superpowers/specs/2026-05-04-legacy-reader-retirement-design.md` (SS-5).

---

## Pre-flight evidence (what makes deletion safe)

- All five regression cases use `.csv` metfiles — confirmed with `grep -n 'file = ".*\.csv"' tests/swap-cases/toml/*/swap.toml` (cases 1, 2, 4, 5, 6 all set `metfile/file = "<NNN>.csv"` in the `[meteorology.temporal]` section).
- Case 3 (3.macroporeflow) is excluded per ADR 0011 and is now stub-errored at the TOML boundary by SS-4.
- The TTutil branches (`if (swMetCSV /= 1)`) are therefore unreachable from working source after the stub-error in Commit 1 lands. Commit 2 deletes the dead code.
- `readswap.f90:419` sets `swMetFilAll = 0` in the legacy `.swp` pipeline, which is no longer invoked by working source code per the umbrella-spec retirement gate.

---

## File Structure

| File | Commit | Change |
|---|---|---|
| `src/config/meteorology_config.f90` | 1 | **Modify:** add a stub-error block to `meteorology_config_validate` rejecting `metfile` not ending in `.csv`. |
| `tests/unit/config/test_meteorology_config.pf` | 1 | **Modify:** append `test_meteo_metfile_must_be_csv` (rejects `.met`, `.YYY`, empty allocatable, etc.). |
| `docs/phase-4f-readmeteo-ttutil-audit.md` | 1 | **Create:** snapshot of pre-deletion state, evidence the branches are dead, deletion plan. |
| `docs/adr/0014-readmeteo-phaseout.md` | 1, 2, 3 | **Modify:** append progress notes after each commit (sequencing-constraint met → Step 2 done → Step 3 done). |
| `src/io/readmeteo.f90` | 2 | **Modify:** delete TTutil branches in `ReadMeteoYear` (~lines 106–144) and `ReadRainEvents` (~lines 367–end-of-subroutine), delete entire `MeteoInOneFile` subroutine (~lines 442–593), simplify control flow. |
| `src/io/toml/config_to_variables.f90` | 2 | **Modify:** delete the `swMetFilAll = 1` adapter block (~lines 176–200). The `swMetCSV` / `swMetDetCSV` / `swRainCSV` switches stay live in this commit but become trivially redundant — they are deleted in Commit 3. |
| `src/core/variables.f90` | 3 | **Modify:** delete `swMetFilAll`, `swMetCSV`, `swRainCSV`, `swMetDetCSV`, `rainfil`, `station(366)`, integer `ad(366)`/`am(366)` arrays per ADR 0014 Step 3. |
| `src/io/toml/config_to_variables.f90` | 3 | **Modify:** delete the now-redundant `swMetCSV = 0/1`, `swRainCSV = 0/1`, `swMetDetCSV = 0/1` assignments — the variables themselves are gone. |
| `src/io/readmeteo.f90` | 3 | **Modify:** delete `use variables, only: ..., swMetCSV, ..., swRainCSV, ...` references; the `if` guards are already gone from Commit 2. |
| `src/legacy/readswap.f90` | 3 | **Modify:** delete the `swMetFilAll = 0` assignment at line 419 (variable no longer exists). |
| `docs/superpowers/specs/2026-05-04-legacy-reader-retirement-design.md` | 3 | **Modify:** flip SS-5 row OPEN → DONE in the status snapshot table and the sub-spec roadmap table. |

No changes to schema, no new TOML keys, no adapter logic changes other than deletions.

---

## Conventions

1. **Test-first for the stub-error.** Write the pFUnit test before touching the validator. Verify it fails. Then add the validator change. Verify it passes.
2. **Stub-error message format** (mirrors SS-4): start with the offending key + value in parentheses, end with the ADR citation. Example: `'meteorology.metfile="weather.met" not supported in the TOML pipeline; only CSV metfiles are accepted (Phase 4f-extend SS-5; legacy .met/.YYY readers retired per ADR 0014).'`
3. **Verification gate at every commit.** Each of the three commits must individually pass `pixi run -e test test-pfunit` AND `pixi run -e test check-full` (5/5 cases). Do not roll commits together — the deletion is large enough that bisection guarantees matter.
4. **Single audit doc** at the start (Commit 1), then ADR 0014 progress notes appended at each subsequent commit.
5. **Case 3 (macroporeflow) expectation.** Already stub-errored by SS-4. After SS-5 it remains stub-errored. No change in behaviour.
6. **Commits are NOT amended after the fact.** Each commit's SHA is what it is — the umbrella spec status row in Commit 3 references the plan filename, not a SHA, matching the SS-4 hygiene fix.

---

## Task 1: Add the failing stub-error pFUnit test

**Commit:** 1
**Files:** `tests/unit/config/test_meteorology_config.pf`
**Verification:** `pixi run -e test test-pfunit -- --filter test_meteo_metfile`

- [ ] **Step 1: Append the new test to `tests/unit/config/test_meteorology_config.pf`**

Append at the bottom of the file (after the last `@test` block):

```fortran
@test
subroutine test_meteo_metfile_dot_met_stub_errors()
   use funit
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   use meteorology_config_mod, only: meteorology_config_t
   type(meteorology_config_t) :: m
   type(error_collection_t)   :: errors
   integer :: i
   logical :: found_metfile_msg

   m%lat = 52.0d0
   m%alt = 10.0d0
   m%metfile = "weather.met"
   call m%validate(errors)

   @assertTrue(errors%has_errors())

   found_metfile_msg = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'meteorology.metfile') > 0 .and. &
          index(errors%items(i)%message, 'ADR 0014') > 0) then
         found_metfile_msg = .true.
         exit
      end if
   end do
   @assertTrue(found_metfile_msg)
end subroutine

@test
subroutine test_meteo_metfile_per_year_extension_stub_errors()
   use funit
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   use meteorology_config_mod, only: meteorology_config_t
   type(meteorology_config_t) :: m
   type(error_collection_t)   :: errors
   integer :: i
   logical :: found

   m%lat = 52.0d0
   m%alt = 10.0d0
   m%metfile = "283"   ! base name, would historically grow .YYY suffix
   call m%validate(errors)

   @assertTrue(errors%has_errors())

   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'meteorology.metfile') > 0) then
         found = .true.
         exit
      end if
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_meteo_metfile_csv_passes()
   use funit
   use error_mod, only: error_collection_t
   use meteorology_config_mod, only: meteorology_config_t
   type(meteorology_config_t) :: m
   type(error_collection_t)   :: errors

   m%lat = 52.0d0
   m%alt = 10.0d0
   m%metfile = "283.csv"
   call m%validate(errors)

   @assertFalse(errors%has_errors())
end subroutine
```

- [ ] **Step 2: Build to make sure the tests compile**

Run: `pixi run -e test build-linux 2>&1 | tail -20`
Expected: build succeeds.

- [ ] **Step 3: Run the new tests to verify they FAIL appropriately**

Run: `pixi run -e test test-pfunit 2>&1 | tail -30`

Expected: 2 failures (the two stub-error tests). The `_csv_passes` test should already pass because the current validator doesn't reject `.csv`. If the per-year test PASSES (i.e. validator already rejects `"283"`), check whether some upstream `check_not_empty` already fires — that's fine, but make sure it produces the expected `ERR_VALIDATION_CROSS_FIELD` code, not `ERR_VALIDATION_RANGE` or similar. Adjust the test if needed.

---

## Task 2: Add the stub-error block to the validator

**Commit:** 1
**Files:** `src/config/meteorology_config.f90`
**Verification:** `pixi run -e test test-pfunit -- --filter test_meteo`

- [ ] **Step 1: Add the stub-error block at the top of `meteorology_config_validate`**

In `src/config/meteorology_config.f90`, locate `subroutine meteorology_config_validate` (line 66). Insert the following block immediately before the existing `call check_real_range(self%lat, ...)` at line 70:

```fortran
      ! ----- Phase 4f-extend SS-5 stub-error: only CSV metfiles supported -----
      ! ADR 0014 retires the TTutil-based readers in readmeteo.f90 (per-year
      ! .YYY daily reader, .met all-years reader, TTutil detail reader, and the
      ! TTutil tail of ReadRainEvents). After SS-5, swap_csv_dat is the only
      ! source of meteorology data. Reject any non-.csv metfile here so that
      ! configurations naming a legacy file produce an actionable error
      ! instead of falling through to deleted code paths.
      if (.not. allocated(self%metfile) .or. len_trim(self%metfile) == 0) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'meteorology.metfile is required and must reference a CSV file ' // &
            '(per ADR 0014; legacy .met / per-year .YYY readers retired ' // &
            'in Phase 4f-extend SS-5).', 'meteorology')
      else if (index(self%metfile, '.csv') == 0) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'meteorology.metfile="' // trim(self%metfile) // &
            '" not supported in the TOML pipeline; only CSV metfiles are ' // &
            'accepted (Phase 4f-extend SS-5; legacy .met/.YYY readers ' // &
            'retired per ADR 0014).', 'meteorology')
      end if

```

The trailing blank line keeps it visually separated from the existing range/enum checks below.

- [ ] **Step 2: Build**

Run: `pixi run -e test build-linux 2>&1 | tail -10`
Expected: build succeeds.

- [ ] **Step 3: Run the targeted tests to verify they PASS**

Run: `pixi run -e test test-pfunit 2>&1 | tail -30`
Expected: all `test_meteo*` tests pass; new tests included.

- [ ] **Step 4: Run check-full to confirm no regression**

Run: `pixi run -e test check-full 2>&1 | tail -15`
Expected: 5/5 cases pass. The validator change only rejects non-CSV metfiles; all working cases use `.csv`, so no behaviour change for them.

---

## Task 3: Write the audit doc

**Commit:** 1
**Files:** `docs/phase-4f-readmeteo-ttutil-audit.md` (new)

- [ ] **Step 1: Create `docs/phase-4f-readmeteo-ttutil-audit.md`**

Create the file with the following content:

```markdown
# Phase 4f-extend SS-5 — readmeteo.f90 TTutil audit

**Date:** 2026-05-05
**Sub-spec:** SS-5 (legacy reader retirement umbrella)
**Plan:** `docs/superpowers/plans/2026-05-05-ss5-readmeteo-ttutil-deletion.md`
**ADR:** 0014 (Steps 2 and 3)

## Pre-deletion state (HEAD before Commit 2)

`src/io/readmeteo.f90` (772 LoC total) contains the following TTutil-using
sections, all reachable only when `swMetCSV /= 1` (i.e. legacy `.met` or
per-year `.YYY` mode):

| Section | Lines | Description |
|---|---|---|
| `ReadMeteoYear` daily TTutil branch | 106–128 | `rdinit` / `rdacha` / `rdfinr` / `rdfdor` per-year `.YYY` reader |
| `ReadMeteoYear` `MeteoInOneFile` call | 107–108 | Legacy `.met` all-years branch |
| `ReadMeteoYear` detail TTutil branch | 129–144 | `rdinit` / `rdatim` / `rdfinr` / `rdfdor` per-year `.YYY` sub-daily reader |
| `MeteoInOneFile` subroutine | 442–593 | All-years `.met` reader |
| `ReadRainEvents` TTutil tail | 367–end of subroutine | `rdinit` / `rdainr` / `rdfinr` / `rdfdor` per-year `.YYY` rain events reader |
| Control-flow scaffolding | `goto 100`, `100 continue`, `if (swMetCSV == 1)`, `if (swRainCSV == 1)` guards | Strangler-fig harness routing CSV vs TTutil |

## Reachability evidence

Confirmed before SS-5 Commit 1:

- All 5 regression cases (1.hupselbrook, 2.grassgrowth, 4.oxygenstress,
  5.salinitystress, 6.surfacewater) set `metfile = "<NNN>.csv"` in
  `[meteorology.temporal]`. None use `.met` or per-year `.YYY` extensions.
  Verified with `grep -n 'file = ".*\.csv"' tests/swap-cases/toml/*/swap.toml`.
- Case 3 (3.macroporeflow) is excluded per ADR 0011 and is now also stub-
  errored at the TOML boundary by SS-4.
- `readswap.f90:419` sets `swMetFilAll = 0` in the legacy `.swp` pipeline,
  which is no longer invoked by working source code per the umbrella-spec
  retirement gate.

After SS-5 Commit 1 (the `.met`/non-CSV stub-error in
`meteorology_config_validate`), the TTutil branches above are unreachable
from any TOML configuration. Commit 2 deletes the dead code; Commit 3
sweeps the now-dead variables.

## Deletion order (Commit 2)

1. `src/io/toml/config_to_variables.f90`: delete the `swMetFilAll = 1` block
   (lines 176–200). After the stub-error, `metfile` always ends in `.csv`,
   so `swMetFilAll = 0` is the only outcome — the conditional and the
   `MeteoInOneFile(1, ...)` pre-load call are unreachable.
2. `src/io/readmeteo.f90` `ReadMeteoYear`:
   - Delete the `else` branch inside `if (swmetdetail == 0)` (the per-year
     `.YYY` daily TTutil reader, lines 109–128).
   - Delete the `if (swMetFilAll == 1)` line and the `else` keyword (line
     107–109); the remaining body becomes unreachable too — it is the
     `MeteoInOneFile(2, ifnd)` call.
   - Delete the `elseif (swmetdetail == 1)` branch (lines 129–144) — TTutil
     detail reader.
   - Delete the `if (swMetCSV == 1) ... goto 100 ... end if` guard
     (lines 97–104). The two `MeteoCSV*` calls become unconditional.
   - Delete the `100 continue` label (line 146).
3. `src/io/readmeteo.f90` `ReadRainEvents`:
   - Delete everything after the `if (swRainCSV == 1) ... return; end if`
     block (lines 367–439). The TTutil per-year `.YYY` rain reader is
     unreachable.
   - Delete the `if (swRainCSV == 1)` guard (line 312); the body becomes
     the unconditional CSV path.
4. `src/io/readmeteo.f90` end-of-file: delete the entire `MeteoInOneFile`
   subroutine (lines 442–593).

## Variable sweep (Commit 3)

Per ADR 0014 Step 3, the following are deleted from `src/core/variables.f90`:

| Variable | Why dead after Commit 2 |
|---|---|
| `swMetFilAll` | Adapter no longer sets it; reader no longer reads it |
| `swMetCSV` | Only one path remains (CSV) — guard is gone |
| `swRainCSV` | Same as `swMetCSV` |
| `swMetDetCSV` | Same |
| `rainfil` | Used only for `.YYY` rain filename construction |
| `station(366)` | Read by `rdacha` in the deleted per-year daily reader |
| `ad(366)`, `am(366)` integer arrays | Read by TTutil per-year daily reader; CSV path uses `days1900_to_md` |

Survivors (still consumed by `meteoday.f90` / `meteodt.f90`):
`metcsv_dat`, `nmetcsv`, `metcsv_det`, `nmetcsv_det`, `raincsv_dat`,
`nraincsv`, all `det*` simulation arrays, `raintimearray`, `rainamount`,
`nmrain`.

Also delete the `swMetFilAll = 0` assignment at `src/legacy/readswap.f90:419`
(variable no longer exists).
```

- [ ] **Step 2: Verify the audit doc renders**

Run: `head -20 docs/phase-4f-readmeteo-ttutil-audit.md`
Expected: front-matter and section headings render cleanly.

---

## Task 4: Update ADR 0014 (Commit 1 progress)

**Commit:** 1
**Files:** `docs/adr/0014-readmeteo-phaseout.md`

- [ ] **Step 1: Append a progress note section to ADR 0014**

At the bottom of the file, append:

```markdown

## Progress note 2026-05-05 — Sequencing constraint resolved (SS-5 Commit 1)

The sequencing constraint above is now closed by the second option:
`meteorology_config_validate` rejects any `metfile` not ending in `.csv`
with `ERR_VALIDATION_CROSS_FIELD`. The legacy `.met` and per-year `.YYY`
codepaths are therefore unreachable from any TOML configuration. The
legacy `.swp` pipeline is no longer invoked by working source code (umbrella
spec `2026-05-04-legacy-reader-retirement-design.md`), so no separate
`readmeteo_legacy.f90` is needed. Steps 2 and 3 may proceed.

Audit doc: `docs/phase-4f-readmeteo-ttutil-audit.md`.
Plan: `docs/superpowers/plans/2026-05-05-ss5-readmeteo-ttutil-deletion.md`.
```

---

## Task 5: Commit 1

**Commit:** 1
**Files:** all of the above + this plan file.

- [ ] **Step 1: Review the staged diff**

```bash
cd /home/zawadzkim/Code/swap
git status --short
git diff --stat
```

Expected files modified:
- `src/config/meteorology_config.f90`
- `tests/unit/config/test_meteorology_config.pf`
- `docs/adr/0014-readmeteo-phaseout.md`

Expected new files:
- `docs/phase-4f-readmeteo-ttutil-audit.md`
- `docs/superpowers/plans/2026-05-05-ss5-readmeteo-ttutil-deletion.md`

- [ ] **Step 2: Stage and commit**

```bash
cd /home/zawadzkim/Code/swap
git add src/config/meteorology_config.f90 \
        tests/unit/config/test_meteorology_config.pf \
        docs/adr/0014-readmeteo-phaseout.md \
        docs/phase-4f-readmeteo-ttutil-audit.md \
        docs/superpowers/plans/2026-05-05-ss5-readmeteo-ttutil-deletion.md
git commit -m "$(cat <<'EOF'
feat(config): stub-error non-CSV metfiles at TOML boundary (SS-5 Commit 1)

Phase 4f-extend SS-5 Commit 1. Reject any meteorology.metfile not ending
in .csv in meteorology_config_validate, mirroring the SS-4 stub-error
pattern. This closes the ADR 0014 sequencing constraint without needing a
separate readmeteo_legacy.f90: the TTutil branches in readmeteo.f90
become unreachable from working source.

- src/config/meteorology_config.f90 — new stub-error block at the top of
  meteorology_config_validate.
- tests/unit/config/test_meteorology_config.pf — three new tests
  (.met rejected, per-year extension rejected, .csv accepted).
- docs/phase-4f-readmeteo-ttutil-audit.md — pre-deletion audit and
  reachability evidence.
- docs/adr/0014-readmeteo-phaseout.md — progress note marking the
  sequencing constraint resolved; Steps 2 and 3 unblocked.
- docs/superpowers/plans/2026-05-05-ss5-readmeteo-ttutil-deletion.md —
  the SS-5 plan.

Verified: pFUnit suite green; check-full 5/5 green.
EOF
)"
```

---

## Task 6: Delete TTutil branches in `ReadMeteoYear`

**Commit:** 2
**Files:** `src/io/readmeteo.f90`
**Verification:** `pixi run -e test build-linux` and `pixi run -e test check-full`

Operate top-down (later line numbers will shift after each deletion — re-grep before each step).

- [ ] **Step 1: Delete the `if (swMetCSV == 1)` guard, keeping the body unconditional**

Locate the block at lines 94–104 (the guard around `MeteoCSVYear` / `MeteoCSVDetYear`). Replace it with an unconditional dispatch:

```fortran
! --- CSV mode is the only supported path (Phase 4f-extend SS-5,
!     ADR 0014). Daily mode → MeteoCSVYear. Sub-daily → MeteoCSVDetYear.
      if (swmetdetail == 0) then
         call MeteoCSVYear(ifnd)
      else
         call MeteoCSVDetYear(ifnd)
      end if
```

Delete the now-dangling `goto 100` (was at line 103).

- [ ] **Step 2: Delete the daily TTutil branch (`if (swmetdetail.eq.0)` block)**

Lines 106–128 (the `if (swmetdetail.eq.0) then ... endif` immediately after the deleted CSV guard). Delete the entire block, including the `if (swMetFilAll == 1) call MeteoInOneFile(2, ifnd)` line and the `else` branch with the `rdinit`/`rdacha`/`rdfinr`/`rdfdor` per-year reader.

- [ ] **Step 3: Delete the detail TTutil branch (`elseif (swmetdetail.eq.1)` block)**

Lines 129–144 (the entire `elseif` branch with the `rdinit` and detail-array reads). Delete the block. The remaining `endif` from the outer `if (swmetdetail.eq.0)` no longer has a partner — delete the `endif` too.

- [ ] **Step 4: Delete the `100 continue` label**

Line 146. Delete the line.

- [ ] **Step 5: Build and verify the file still compiles**

Run: `pixi run -e test build-linux 2>&1 | tail -15`
Expected: build succeeds. If a label-related error fires, double-check Step 4 caught all references to `100`.

---

## Task 7: Delete TTutil tail in `ReadRainEvents`

**Commit:** 2
**Files:** `src/io/readmeteo.f90`

- [ ] **Step 1: Delete the `if (swRainCSV == 1)` guard, making the CSV path unconditional**

Locate `subroutine ReadRainEvents` (currently line 277). Inside, find `if (swRainCSV == 1) then ... return; end if` (currently lines 312–365).

Replace with an unconditional CSV body — i.e. delete the `if (swRainCSV == 1) then` line, the matching `end if` (~line 365), and the trailing `return` that immediately follows the `end if`. The body of the former `if` becomes the body of the subroutine.

- [ ] **Step 2: Delete the TTutil tail**

Lines 367–end-of-subroutine (everything from the `!========================= TTutil / legacy .YYY path ===================` banner down to and including the final `end subroutine ReadRainEvents`). Re-add the `end subroutine ReadRainEvents` line so the subroutine is properly terminated.

Also delete the `swRainCSV` import in the `use variables, only: ...` clause near the top of the subroutine (line 286–287). Likewise drop unused locals: `chday`, `chmonth`, `chtime`, `ext`, `pre`, `flrnx`, `tdayold`, `time4`, integer arrays `ad`/`am`/`ay` of size `mrain`.

- [ ] **Step 3: Build**

Run: `pixi run -e test build-linux 2>&1 | tail -15`
Expected: build succeeds. Watch for unused-variable warnings — if gfortran complains, finish removing the unused locals listed in Step 2.

---

## Task 8: Delete `MeteoInOneFile` subroutine

**Commit:** 2
**Files:** `src/io/readmeteo.f90`

- [ ] **Step 1: Delete the entire `MeteoInOneFile` subroutine**

Lines 442–593 (from `subroutine MeteoInOneFile (iTask, ifnd)` to `end subroutine MeteoInOneFile`). Delete the whole block.

- [ ] **Step 2: Build**

Run: `pixi run -e test build-linux 2>&1 | tail -15`
Expected: build succeeds. The interface block in `config_to_variables.f90` will still be there from the next task — that's fine until Task 9 deletes it.

---

## Task 9: Delete `swMetFilAll` adapter block

**Commit:** 2
**Files:** `src/io/toml/config_to_variables.f90`

- [ ] **Step 1: Delete the `.met` extension branch**

Lines 176–186 (from `else if (index(trim(metfil), '.met') > 0) then` through the matching `end if`). Delete the block. The preceding `if (index(trim(metfil), '.csv') > 0) then ... end if` becomes a standalone `if`.

- [ ] **Step 2: Delete the `swMetFilAll == 1` pre-load block**

Lines 188–200 (from `! Pre-load all years into cache for legacy .met mode.` through the matching `end if`). Delete the block, including the `interface ... subroutine MeteoInOneFile ... end interface` declaration and the `call MeteoInOneFile(1, idum_meteo)`.

- [ ] **Step 3: Update the comment block above the dispatch**

Lines 139–142 (the `! Detect meteo file mode from the metfil extension.` comment block). Trim it to:

```fortran
      ! All metfile extensions other than .csv are rejected by
      ! meteorology_config_validate (Phase 4f-extend SS-5; ADR 0014).
      ! Pre-load CSV via read_csv_table.
      call lowerc(metfil)
      swMetCSV = 0
```

(`swMetFilAll = 0` line goes away — variable will be deleted in Commit 3, but for now leaving it set to 0 would still work; cleaner to drop the line now since we're touching the block.)

- [ ] **Step 4: Build**

Run: `pixi run -e test build-linux 2>&1 | tail -15`
Expected: build succeeds. `swMetFilAll` is still defined in `variables.f90` so no unresolved reference yet.

---

## Task 10: Verify Commit 2 candidate

**Commit:** 2
**Files:** none (verification only)

- [ ] **Step 1: Full pFUnit suite**

Run: `pixi run -e test test-pfunit 2>&1 | tail -15`
Expected: all green.

- [ ] **Step 2: check-full**

Run: `pixi run -e test check-full 2>&1 | tail -15`
Expected: 5/5 cases pass.

If anything regresses, do NOT commit. Investigate. Most likely failure modes:
- A test inadvertently constructed a `meteorology_config_t` with a non-CSV metfile (the new validator now rejects it — fix the test).
- A surviving reference to `MeteoInOneFile`, `swMetCSV`, or label `100` (re-grep and clean up).

- [ ] **Step 3: Update ADR 0014 (Commit 2 progress)**

Append to `docs/adr/0014-readmeteo-phaseout.md`:

```markdown

## Progress note 2026-05-05 — Step 2 complete (SS-5 Commit 2)

`readmeteo.f90` no longer contains TTutil calls. `MeteoInOneFile` deleted.
The `swMetFilAll = 1` adapter block in `config_to_variables.f90` deleted.
`ReadRainEvents` is CSV-only. Deleted ~150 LoC. Variables (`swMetCSV`,
`swRainCSV`, `swMetFilAll`, `swMetDetCSV`, `rainfil`, `station`, `ad`,
`am`) become trivially redundant — Step 3 sweeps them.
```

---

## Task 11: Commit 2

**Commit:** 2

- [ ] **Step 1: Review the staged diff**

```bash
cd /home/zawadzkim/Code/swap
git status --short
git diff --stat
```

Expected files modified:
- `src/io/readmeteo.f90` (large deletion)
- `src/io/toml/config_to_variables.f90`
- `docs/adr/0014-readmeteo-phaseout.md`

The diff stat for `readmeteo.f90` should show roughly `-150` lines.

- [ ] **Step 2: Stage and commit**

```bash
cd /home/zawadzkim/Code/swap
git add src/io/readmeteo.f90 \
        src/io/toml/config_to_variables.f90 \
        docs/adr/0014-readmeteo-phaseout.md
git commit -m "$(cat <<'EOF'
refactor(io): delete TTutil branches from readmeteo.f90 (SS-5 Commit 2)

Phase 4f-extend SS-5 Commit 2 / ADR 0014 Step 2. Delete the four TTutil
branches in readmeteo.f90 (per-year .YYY daily reader, .met all-years
reader via MeteoInOneFile, TTutil detail reader, TTutil tail of
ReadRainEvents) and the matching swMetFilAll pre-load block in the
adapter. Simplify control flow: remove the goto 100 / 100 continue label
dance and the swMetCSV / swRainCSV guards (only one path remains).

The non-CSV metfile stub-error from Commit 1 makes these branches
unreachable from any TOML configuration; the legacy .swp pipeline is
not invoked by working source per the umbrella spec retirement gate.

Files:
- src/io/readmeteo.f90 — ~150 LoC deleted (ReadMeteoYear TTutil daily +
  detail branches, MeteoInOneFile entire subroutine, ReadRainEvents
  TTutil tail, control-flow scaffolding).
- src/io/toml/config_to_variables.f90 — swMetFilAll = 1 block deleted;
  comment block trimmed.
- docs/adr/0014-readmeteo-phaseout.md — Step 2 progress note appended.

Variables still defined in variables.f90 (swMetCSV / swRainCSV /
swMetDetCSV / swMetFilAll / rainfil / station / ad / am) become
trivially redundant — Commit 3 sweeps them.

Verified: pFUnit suite green; check-full 5/5 green.
EOF
)"
```

---

## Task 12: Dead-variable sweep in `variables.f90`

**Commit:** 3
**Files:** `src/core/variables.f90`, `src/io/toml/config_to_variables.f90`, `src/io/readmeteo.f90`, `src/legacy/readswap.f90`
**Verification:** `pixi run -e test build-linux` and `pixi run -e test check-full`

- [ ] **Step 1: Find every reference to the doomed variables**

Run:

```bash
grep -rn "swMetCSV\|swRainCSV\|swMetFilAll\|swMetDetCSV\|rainfil\b\|\bstation(366)\|station\b" src/ --include="*.f90"
```

(For `station`, check carefully — there may be unrelated `station` identifiers in other modules; the deletion target is the `station(366)` array in `variables.f90`.)

Run this set of greps to enumerate references for `ad(366)` / `am(366)`:

```bash
grep -n "integer.*ad(366)\|integer.*am(366)\|\bad(\|\bam(" src/core/variables.f90 src/io/readmeteo.f90 src/legacy/readswap.f90 | head -40
```

Note the survivors in `readmeteo.f90` (the local `ad`/`am` scratch arrays inside subroutines should remain — only the module-level globals are deleted).

- [ ] **Step 2: Delete the variable declarations from `src/core/variables.f90`**

Lines 204, 205, 216, 218 (and any related lines for `rainfil`, `station(366)`, integer `ad(366)`/`am(366)` — exact line numbers may have shifted; re-grep). Delete:
- `integer   swMetFilAll`
- `integer   swMetCSV`
- `integer   swMetDetCSV`
- `integer   swRainCSV`
- `character(len=...) rainfil`
- `character(len=...) station(366)`
- `integer   ad(366)`
- `integer   am(366)`

- [ ] **Step 3: Delete adapter assignments**

In `src/io/toml/config_to_variables.f90`, delete the now-redundant assignments:
- `swMetCSV = 0` and `swMetCSV = 1` (around lines 144 / 148)
- `swMetDetCSV = 0` and `swMetDetCSV = 1` (around lines 203 / 235)
- `swRainCSV = 0` and `swRainCSV = 1` (around lines 240 / 243)

The CSV pre-load blocks themselves stay — only the switch assignments go.

- [ ] **Step 4: Delete `use variables` references**

In `src/io/readmeteo.f90`, find any `use variables, only: ..., swMetCSV, swRainCSV, swMetDetCSV, swMetFilAll, ...` clauses and remove the deleted names.

- [ ] **Step 5: Delete the legacy assignment**

In `src/legacy/readswap.f90:419`, delete the line `swMetFilAll = 0`.

- [ ] **Step 6: Build**

Run: `pixi run -e test build-linux 2>&1 | tail -20`
Expected: build succeeds. If gfortran flags an unresolved reference, grep for it across `src/` and delete that reference too. Common surviving references to look for:

```bash
grep -rn "swMetCSV\|swRainCSV\|swMetFilAll\|swMetDetCSV" src/ --include="*.f90"
```

(Should return zero results after this task.)

- [ ] **Step 7: Verify**

Run: `pixi run -e test test-pfunit 2>&1 | tail -15`
Expected: all green.

Run: `pixi run -e test check-full 2>&1 | tail -15`
Expected: 5/5 cases pass.

---

## Task 13: Update ADR 0014 (Step 3 done) and umbrella spec

**Commit:** 3
**Files:** `docs/adr/0014-readmeteo-phaseout.md`, `docs/superpowers/specs/2026-05-04-legacy-reader-retirement-design.md`

- [ ] **Step 1: Append Step 3 progress note to ADR 0014**

```markdown

## Progress note 2026-05-05 — Step 3 complete (SS-5 Commit 3)

Dead-variable sweep complete. Deleted from `variables.f90`: `swMetFilAll`,
`swMetCSV`, `swRainCSV`, `swMetDetCSV`, `rainfil`, `station(366)`, integer
`ad(366)`/`am(366)` arrays. Adapter assignments and the legacy
`swMetFilAll = 0` line in `readswap.f90:419` deleted alongside.

ADR 0014 phase-out is now complete. The TTutil library remains a
build-time dependency of the legacy `.swp` target only (if that target
is retained); the TOML-only build target no longer needs it.
```

- [ ] **Step 2: Flip umbrella-spec SS-5 row to DONE**

In `docs/superpowers/specs/2026-05-04-legacy-reader-retirement-design.md`, locate the status snapshot table row:

```
| CSV meteo Steps 2–3 (TTutil branch deletion + dead-var sweep) | OPEN | Daily-meteo TTutil fallback still alive |
```

Change to:

```
| CSV meteo Steps 2–3 (TTutil branch deletion + dead-var sweep) | DONE | SS-5 — TTutil branches retired (`2026-05-05-ss5-readmeteo-ttutil-deletion.md`) |
```

Also locate the sub-spec roadmap row:

```
| SS-5 | `readmeteo.f90` TTutil branch deletion (ADR 0014 Steps 2–3) | OPEN | new spec | — |
```

Change to:

```
| SS-5 | `readmeteo.f90` TTutil branch deletion (ADR 0014 Steps 2–3) | DONE | `2026-05-05-ss5-readmeteo-ttutil-deletion.md` | — |
```

(Convention matches SS-4: reference the plan filename, no SHA pin.)

---

## Task 14: Commit 3

**Commit:** 3

- [ ] **Step 1: Review the staged diff**

```bash
cd /home/zawadzkim/Code/swap
git status --short
git diff --stat
```

Expected files modified:
- `src/core/variables.f90`
- `src/io/toml/config_to_variables.f90`
- `src/io/readmeteo.f90`
- `src/legacy/readswap.f90`
- `docs/adr/0014-readmeteo-phaseout.md`
- `docs/superpowers/specs/2026-05-04-legacy-reader-retirement-design.md`

- [ ] **Step 2: Stage and commit**

```bash
cd /home/zawadzkim/Code/swap
git add src/core/variables.f90 \
        src/io/toml/config_to_variables.f90 \
        src/io/readmeteo.f90 \
        src/legacy/readswap.f90 \
        docs/adr/0014-readmeteo-phaseout.md \
        docs/superpowers/specs/2026-05-04-legacy-reader-retirement-design.md
git commit -m "$(cat <<'EOF'
refactor(core): sweep dead meteo variables; close SS-5 (Commit 3)

Phase 4f-extend SS-5 Commit 3 / ADR 0014 Step 3. Delete the variables
that became redundant after Commit 2 retired the TTutil branches in
readmeteo.f90: swMetFilAll, swMetCSV, swRainCSV, swMetDetCSV, rainfil,
station(366), integer ad(366)/am(366) arrays. Drop the adapter
assignments and the legacy swMetFilAll = 0 line in readswap.f90.

Files:
- src/core/variables.f90 — declarations deleted.
- src/io/toml/config_to_variables.f90 — switch assignments deleted; CSV
  pre-load blocks retained.
- src/io/readmeteo.f90 — `use variables, only:` clauses pruned.
- src/legacy/readswap.f90 — swMetFilAll = 0 (line 419) deleted.
- docs/adr/0014-readmeteo-phaseout.md — Step 3 progress note appended;
  ADR phase-out marked complete.
- docs/superpowers/specs/2026-05-04-legacy-reader-retirement-design.md —
  SS-5 row flipped OPEN → DONE in both tables.

Verified: pFUnit suite green; check-full 5/5 green.
EOF
)"
```

- [ ] **Step 3: Final verification**

```bash
cd /home/zawadzkim/Code/swap
git log --oneline -4
pixi run -e test check-full 2>&1 | tail -10
grep -rn "swMetCSV\|swRainCSV\|swMetFilAll\|swMetDetCSV" src/ --include="*.f90" | wc -l
```

Expected:
- Three new commits visible (`SS-5 Commit 1/2/3`).
- check-full: 5/5 passed.
- grep count: 0.

---

## Definition of done

- `meteorology_config_validate` rejects any non-CSV `metfile` with `ERR_VALIDATION_CROSS_FIELD` citing ADR 0014.
- `test_meteo_metfile_dot_met_stub_errors`, `test_meteo_metfile_per_year_extension_stub_errors`, and `test_meteo_metfile_csv_passes` all pass.
- `readmeteo.f90` contains zero TTutil calls (`grep -n "rdinit\|rdacha\|rdfinr\|rdfdor\|rdatim\|rdainr" src/io/readmeteo.f90` returns nothing).
- `MeteoInOneFile` deleted; no remaining references in any `.f90` under `src/`.
- `swMetFilAll`, `swMetCSV`, `swRainCSV`, `swMetDetCSV`, `rainfil`, `station(366)`, integer `ad(366)`/`am(366)` deleted from `src/core/variables.f90`; no remaining references anywhere under `src/`.
- pFUnit suite green at every commit.
- check-full 5/5 green at every commit.
- ADR 0014 records the phase-out as complete with three progress notes.
- Umbrella spec SS-5 status flipped DONE in both tables.
- Three commits, each independently shippable.
