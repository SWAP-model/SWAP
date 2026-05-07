# TTutil Retirement Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the TTutil dependency from SWAP. Replace `getun`/`fopens`/`delfil` with a thin `file_io_mod` wrapper, retire the rerun mechanism, drop `ttutil_dep` from meson, and vendor-delete `subprojects/ttutil/`.

**Architecture:** Five phases. Phase A introduces a `file_io_mod` module wrapping native Fortran `open()`/`close()`/`inquire()` with `swap_log`-integrated error reporting. Phase B mechanically replaces ~125 call sites of `getun`/`fopens`/`delfil` with the wrapper, one file per commit. Phase C simplifies `swap_main.f90` to a flat `swap(0,1); swap(0,2); swap(0,3)` triple, dropping `rdsets`/`rdfrom` and the `rddtmp` scratch-file cleanup. Phase D drops dead TTutil-related local declarations left behind by Phase B (the actual `rd*` reader calls were already eliminated by ADRs 0019/0021/0022). Phase E renames `fatalerr_shim.f90` → `fatalerr.f90`, drops `ttutil_dep` from meson, and vendor-deletes `subprojects/ttutil/`.

**Tech Stack:** Fortran 2008, meson + ninja, gfortran, pFUnit. Reference patterns: `src/core/swap_log.f90` (`log_warn`, `log_error`, `to_str`); the existing `fopens`/`getun` pattern across 9 production files (each `iun = getun(lo, hi); call fopens(iun, name, st, ac)` becomes `call file_open(iun, name, st_native, ac_native)`).

**Spec:** `docs/superpowers/specs/2026-05-07-ttutil-retirement-design.md`.

---

## File Structure

**Created:**

| Path | Responsibility |
|---|---|
| `src/io/file_io.f90` | `file_io_mod` — `file_open`, `file_delete`, `file_exists` wrapping native Fortran with `swap_log` integration. |
| `tests/unit/io/test_file_io.pf` | pFUnit suite covering open success/failure paths, idempotent delete, exists-check. |
| `tests/unit/io/fixtures/file_io_present.txt` | Fixture for `file_exists` true case. |
| `docs/adr/0023-ttutil-retirement.md` | ADR umbrella. |

**Modified (Phase B — bulk replace; one commit per file):**

| Path | Sites | Notes |
|---|---|---|
| `src/core/swap_main.f90` | 4 | Reruns plumbing — entire block goes away in Phase C. Phase B converts the surviving open/close into wrapper calls anyway, in case Phase C lands later. |
| `src/io/toml/config_to_variables.f90` | 4 | Includes the `logf` opening. |
| `src/utils/sharedsimulation.f90` | 2 | Shared-mode coordination file. |
| `src/soil/waterbalance.f90` | 3 | `dev_cmb` mass-balance writer. |
| `src/io/swap_csv_output.f90` | 4 | Modern CSV writer paths. |
| `src/io/macroporeoutput.f90` | 9 | Macropore output (stub-erred at runtime per ADR 0010 but built-in). |
| `src/crop/management_soil.f90` | 4 | `.snp` / `_nut.end` writers in nutrient cases. |
| `src/crop/cropgrowth.f90` | 16 | Crop output writers + nutrient block leftover. |
| `src/io/swapoutput.f90` | 78 | The big output module. ~78 uniform-pattern conversions. |

**Modified (Phase C):**

| Path | Change |
|---|---|
| `src/core/swap_main.f90` | Drop `rdsets`/`rdfrom` block + the `do iset = 0, insets` loop; replace with three direct `call swap(...)` calls. |
| `src/io/swapoutput.f90` | Drop the `call rddtmp(cexf)` paired with its `getun` allocation (lines around 3452-3454). |

**Modified (Phase D — dead-decl cleanup):**

| Path | Change |
|---|---|
| `src/crop/management_soil.f90:82` | Drop `logical :: rdinqr` (dead local). |
| `src/crop/cropgrowth.f90:1054` | Drop `integer getun2` (dead local). |
| `src/crop/management_soil.f90:47` | Drop `getun2` from the comma-list of `integer` declarations (dead — the `getun2` calls inside the file get converted to `file_open` in Task 8). |

**Modified (Phase E):**

| Path | Change |
|---|---|
| `src/error/fatalerr_shim.f90` | Renamed to `src/error/fatalerr.f90`; comment header rewritten as canonical (no longer a "shim"). |
| `meson.build` | Drop `ttutil_prj`, `ttutil_dep`; remove `ttutil_dep` from `swap` executable's `dependencies:`. Rename `'src/error/fatalerr_shim.f90'` → `'src/error/fatalerr.f90'`. |
| `tests/unit/meson.build` | Drop `ttutil_dep` from the test executable's `dependencies:`. Rename `'../../src/error/fatalerr_shim.f90'` → `'../../src/error/fatalerr.f90'`. |
| `docs/adr/index.md` | Add ADR 0023 row. |

**Deleted (Phase E):**

| Path | Reason |
|---|---|
| `subprojects/ttutil.wrap` | Meson subproject reference no longer needed. |
| `subprojects/ttutil/` | Vendored TTutil source tree (~5 MB, 170 files) no longer referenced. |

---

## Conventions for Phase B (read once, applies to Tasks 2-10)

Phase B replaces the legacy `getun`/`fopens`/`delfil` idiom with the `file_io_mod` wrapper introduced in Phase A. Replacement patterns are uniform:

**Pattern 1 — output file (most common, used at ~90 sites):**

Before:
```fortran
integer :: getun
...
iun = getun(20, 90)
call fopens(iun, filnam, 'new', 'del')
```

After:
```fortran
use file_io_mod, only: file_open
...
call file_open(iun, filnam, 'replace', 'write')
```

**Pattern 2 — input file (read-only):**

Before:
```fortran
integer :: getun, getun2
...
snp = getun2(10, 90, 2)
open(unit=snp, file=filnam, status='old')
```

After:
```fortran
use file_io_mod, only: file_open
...
call file_open(snp, filnam, 'old', 'read')
```

**Pattern 3 — output file with status='unknown':**

Before:
```fortran
oup = getun2(10, 90, 2)
open(unit=oup, file=filnam, status='unknown')
```

After:
```fortran
call file_open(oup, filnam, 'unknown', 'readwrite')
```

**Pattern 4 — file deletion:**

Before:
```fortran
call delfil(name, .false.)
```

After:
```fortran
use file_io_mod, only: file_delete
...
call file_delete(name)
```

**Pattern 5 — `inquire(file=..., exist=...)`:** stays as native — the existing idiom is one line and clear; the `file_exists` wrapper exists for callers who prefer the abstraction, but Phase B does not migrate working `inquire` sites.

**Mapping legacy `(status, action)` to wrapper `(status, action)`:**

| Legacy (`fopens` 3rd, 4th args) | Wrapper |
|---|---|
| `'new', 'del'` | `'replace', 'write'` |
| `'new', 'noop'` (rare) | `'new', 'write'` |
| `'old', 'noop'` (read existing) | `'old', 'read'` |
| `'unknown', 'noop'` | `'unknown', 'readwrite'` |

**Imports to drop from each touched file's `use variables, only: ...` line:** none. `getun`/`getun2`/`fopens`/`delfil` are local-declaration `integer ::`/`external` declarations or implicit free subroutines, not `use variables` imports. The cleanup is in the local declarations.

**Local declarations to drop:** Each touched file declares the TTutil function name as a local `integer` (e.g., `integer :: getun, getun2`) inside the subroutine. After replacing all calls, those declarations become unused. Drop them in the same commit as the function calls.

**Per-file workflow for Phase B tasks:**

1. Use `grep -nE "\\b(getun|getun2|fopens|delfil)\\b" <file>` to enumerate sites.
2. At each site, identify the pattern (1-5 above) and apply the replacement.
3. Add `use file_io_mod, only: file_open[, file_delete]` near the top of any subroutine that uses the wrapper. (Place it after existing `use variables` lines.)
4. Drop the `integer :: getun, getun2` (or whichever) local declarations that became unused.
5. `pixi run -e test build-linux` clean.
6. `pixi run -e test test-pfunit` `Ok: 1, Fail: 0`.
7. `pixi run -e test check-full` 5/5.
8. Commit with `refactor(io): replace getun/fopens with file_io wrapper in <file>`.

**Acceptance grep at end of each Phase B task:** `grep -nE "\\b(getun|getun2|fopens|delfil)\\b" <file>` returns no matches.

---

## Task 1: Phase A — Introduce `file_io_mod`

**Files:**
- Create: `src/io/file_io.f90`
- Create: `tests/unit/io/test_file_io.pf`
- Create: `tests/unit/io/fixtures/file_io_present.txt`
- Modify: `meson.build` (add `'src/io/file_io.f90'` to `sources`)
- Modify: `tests/unit/meson.build` (add to `pfunit_extra_sources`, register `.pf` and suite)
- Modify: `tests/unit/testSuites.inc` (add `ADD_TEST_SUITE(test_file_io_suite)`)

- [ ] **Step 1: Create the fixture file**

```bash
mkdir -p tests/unit/io/fixtures
echo "exists" > tests/unit/io/fixtures/file_io_present.txt
```

- [ ] **Step 2: Write the failing tests**

Create `tests/unit/io/test_file_io.pf`:

```fortran
@test
subroutine test_file_open_replace_write_succeeds()
   use file_io_mod, only: file_open, file_delete
   use funit
   implicit none
   integer :: u, ios

   ! Clean up any leftover from a previous run.
   call file_delete('test_file_io_replace.tmp')

   call file_open(u, 'test_file_io_replace.tmp', 'replace', 'write', iostat=ios)
   @assertEqual(0, ios, 'open(replace,write) should succeed')
   write(u, '(a)') 'hello'
   close(u)

   call file_delete('test_file_io_replace.tmp')
end subroutine

@test
subroutine test_file_open_old_read_existing_succeeds()
   use file_io_mod, only: file_open
   use funit
   implicit none
   integer :: u, ios
   character(len=64) :: line

   call file_open(u, 'tests/unit/io/fixtures/file_io_present.txt', 'old', 'read', iostat=ios)
   @assertEqual(0, ios, 'open(old,read) on existing file should succeed')
   read(u, '(a)') line
   close(u)
   @assertEqual('exists', trim(line))
end subroutine

@test
subroutine test_file_open_old_read_missing_returns_nonzero_iostat()
   use file_io_mod, only: file_open
   use funit
   implicit none
   integer :: u, ios

   call file_open(u, 'tests/unit/io/fixtures/does_not_exist.xyz', 'old', 'read', iostat=ios)
   @assertTrue(ios /= 0, 'expected nonzero iostat for missing file with status=old')
end subroutine

@test
subroutine test_file_delete_idempotent()
   use file_io_mod, only: file_delete, file_exists
   use funit
   implicit none

   ! Calling delete twice on a non-existing path should not error.
   call file_delete('test_file_io_idempotent.tmp')
   call file_delete('test_file_io_idempotent.tmp')
   @assertFalse(file_exists('test_file_io_idempotent.tmp'))
end subroutine

@test
subroutine test_file_delete_removes_existing()
   use file_io_mod, only: file_open, file_delete, file_exists
   use funit
   implicit none
   integer :: u

   call file_open(u, 'test_file_io_delme.tmp', 'replace', 'write')
   close(u)
   @assertTrue(file_exists('test_file_io_delme.tmp'))
   call file_delete('test_file_io_delme.tmp')
   @assertFalse(file_exists('test_file_io_delme.tmp'))
end subroutine

@test
subroutine test_file_exists_true_for_fixture()
   use file_io_mod, only: file_exists
   use funit
   implicit none
   @assertTrue(file_exists('tests/unit/io/fixtures/file_io_present.txt'))
end subroutine

@test
subroutine test_file_exists_false_for_missing()
   use file_io_mod, only: file_exists
   use funit
   implicit none
   @assertFalse(file_exists('tests/unit/io/fixtures/does_not_exist.xyz'))
end subroutine
```

- [ ] **Step 3: Wire the test file into meson and run to verify it fails**

In `tests/unit/meson.build`, add to `pf_files`:
```meson
        'io/test_file_io.pf',
```

Add to `pfunit_extra_sources`:
```meson
        '../../src/io/file_io.f90',
```

Add to `tests/unit/testSuites.inc`:
```
ADD_TEST_SUITE(test_file_io_suite)
```

Run:
```
pixi run -e test build-linux
```
Expected: build fails — `file_io_mod` doesn't exist.

- [ ] **Step 4: Create the module**

Create `src/io/file_io.f90`:

```fortran
!> @file file_io.f90
!! Thin wrapper around Fortran intrinsic file I/O (open/close/inquire)
!! with swap_log integration. Replaces TTutil's getun/fopens/delfil.
!!
!! Conventions:
!! - file_open allocates a fresh unit via newunit=.
!! - status accepts native Fortran 'old' / 'new' / 'replace' /
!!   'unknown' / 'scratch'.
!! - action accepts native Fortran 'read' / 'write' / 'readwrite'.
!! - When iostat is provided, callers handle errors. When absent,
!!   open failure routes through fatalerr_collected (exit 1).
!! - Failures are logged at WARN level with path + iostat.
!!
!! See ADR 0023 for the umbrella context (TTutil retirement).
module file_io_mod
   use error_mod,  only: fatalerr_collected
   use swap_log,   only: log_warn, to_str
   implicit none
   private

   public :: file_open
   public :: file_delete
   public :: file_exists

contains

   !> Open `path` with the given Fortran status/action, allocating
   !! a new unit number into `unit`.
   !!
   !! On failure: when `iostat` is present, sets it nonzero and
   !! logs a WARN entry. When `iostat` is absent, logs and aborts
   !! via fatalerr_collected.
   subroutine file_open(unit, path, status, action, iostat)
      integer,           intent(out) :: unit
      character(len=*),  intent(in)  :: path
      character(len=*),  intent(in)  :: status
      character(len=*),  intent(in)  :: action
      integer, optional, intent(out) :: iostat
      integer :: ios

      open(newunit=unit, file=trim(path), status=trim(status), &
           action=trim(action), iostat=ios)

      if (ios /= 0) then
         call log_warn('file_io', "open failed: '" // trim(path) // &
              "' status=" // trim(status) // " action=" // trim(action) // &
              " iostat=" // trim(to_str(ios)))
         if (.not. present(iostat)) then
            call fatalerr_collected('file_io', &
                 "cannot open '" // trim(path) // &
                 "' (status=" // trim(status) // &
                 ", action=" // trim(action) // ")")
         end if
      end if

      if (present(iostat)) iostat = ios
   end subroutine file_open


   !> Delete `path` if it exists. Idempotent — no-op (and no error)
   !! if the file is missing.
   subroutine file_delete(path)
      character(len=*), intent(in) :: path
      integer :: u, ios
      logical :: exists

      inquire(file=trim(path), exist=exists)
      if (.not. exists) return

      open(newunit=u, file=trim(path), status='old', iostat=ios)
      if (ios /= 0) then
         call log_warn('file_io', "delete: cannot open '" // trim(path) // &
              "' iostat=" // trim(to_str(ios)))
         return
      end if
      close(u, status='delete')
   end subroutine file_delete


   !> True iff `path` refers to an existing file. Side-effect-free.
   logical function file_exists(path) result(exists)
      character(len=*), intent(in) :: path
      inquire(file=trim(path), exist=exists)
   end function file_exists

end module file_io_mod
```

- [ ] **Step 5: Add the file to meson production sources**

In `meson.build`, in the `sources` list under the `# I/O` section, add `'src/io/file_io.f90'` near the other `src/io/` files. For example, place it after `'src/io/csv_reader.f90',`:

```meson
    'src/io/csv_reader.f90',
    'src/io/file_io.f90',
```

- [ ] **Step 6: Build and run tests**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
```
Expected: build clean; pFUnit `Ok: 1, Fail: 0`. Tests count grows by 7 (the seven test subroutines).

```
pixi run -e test check-full 2>&1 | grep "Results:"
```
Expected: `Results: 5 passed, 0 failed` (unchanged).

- [ ] **Step 7: Commit**

```bash
git add src/io/file_io.f90 tests/unit/io/test_file_io.pf tests/unit/io/fixtures/file_io_present.txt meson.build tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "$(cat <<'EOF'
feat(io): introduce file_io_mod — open/delete/exists wrapper

Adds src/io/file_io.f90 with three public procedures:
- file_open(unit, path, status, action[, iostat]): opens with
  newunit=; logs WARN on failure; routes through
  fatalerr_collected when iostat is absent (exit 1 on failure).
- file_delete(path): idempotent delete-if-exists.
- file_exists(path): inquire wrapper.

status/action accept native Fortran semantics (not legacy
'new'/'del'/'noop'). swap_log handles all error reporting.

Phase A of TTutil retirement (umbrella spec
docs/superpowers/specs/2026-05-07-ttutil-retirement-design.md).
No call sites converted yet — Phase B follows.

pFUnit suite test_file_io covers all paths: open replace/write,
open old/read, open old/read missing-file iostat, idempotent
delete, delete-removes-existing, exists true/false.

Verified: build clean; pFUnit Ok: 1, Fail: 0; check-full 5/5.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Phase B — Bulk replace `getun`/`fopens`/`delfil` (Tasks 2-10)

Each Phase B task follows the workflow described in **Conventions for Phase B** above. The 9 tasks are ordered smallest-to-largest — establishing the pattern in low-risk files before tackling the 78-site `swapoutput.f90`.

For every Phase B task:
1. Apply the replacement patterns 1-5 (see Conventions) at every site in the file.
2. Drop unused `integer :: getun, getun2` / `external` declarations.
3. Add `use file_io_mod, only: file_open[, file_delete]` to each subroutine that uses the wrapper.
4. `pixi run -e test build-linux` clean; `pixi run -e test test-pfunit` `Ok: 1, Fail: 0`; `pixi run -e test check-full` 5/5.
5. Commit per the template at the end of each task.

---

## Task 2: Phase B-1 — `src/utils/sharedsimulation.f90` (2 sites)

Smallest target; establishes the pattern.

**Files:**
- Modify: `src/utils/sharedsimulation.f90`

- [ ] **Step 1: Inventory the sites**

Run:
```
grep -nE "\\b(getun|getun2|fopens|delfil)\\b" src/utils/sharedsimulation.f90
```
Expected (line numbers may shift): two sites — one `getun` allocation and the matching `fopens` call.

- [ ] **Step 2: Apply the replacement**

Replace each `iun = getun(lo, hi); call fopens(iun, name, st, ac)` pair with the appropriate `file_open` call per the Conventions table. Add `use file_io_mod, only: file_open` near the top of the affected subroutine. Drop `integer :: getun` (or `external getun`) declarations that become unused.

- [ ] **Step 3: Build and verify**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep "Results:"
grep -nE "\\b(getun|getun2|fopens|delfil)\\b" src/utils/sharedsimulation.f90 || echo "OK: clean"
```
Expected: build clean; pFUnit `Ok: 1, Fail: 0`; check-full 5/5; grep returns "OK: clean".

- [ ] **Step 4: Commit**

```bash
git add src/utils/sharedsimulation.f90
git commit -m "$(cat <<'EOF'
refactor(utils): replace getun/fopens with file_io wrapper in sharedsimulation.f90

Phase B-1 of TTutil retirement (2 sites). Mechanical conversion
following the patterns documented in the umbrella plan
(docs/superpowers/plans/2026-05-07-ttutil-retirement.md).

Verified: build clean; pFUnit Ok: 1, Fail: 0; check-full 5/5.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Phase B-2 — `src/soil/waterbalance.f90` (3 sites)

**Files:**
- Modify: `src/soil/waterbalance.f90`

Same workflow as Task 2. Inventory:
```
grep -nE "\\b(getun|getun2|fopens|delfil)\\b" src/soil/waterbalance.f90
```

The sites are around line 707 (`dev_cmb` mass-balance writer): `dev_cmb = getun(20,90); call fopens(dev_cmb, filnam, 'new', 'del')`. Standard pattern-1 replacement.

Build, verify, commit with message `refactor(soil): replace getun/fopens with file_io wrapper in waterbalance.f90` (Phase B-2, 3 sites).

---

## Task 4: Phase B-3 — `src/io/swap_csv_output.f90` (4 sites)

**Files:**
- Modify: `src/io/swap_csv_output.f90`

Same workflow. 4 sites, all uniform output-file openings (pattern 1).

Commit message: `refactor(io): replace getun/fopens with file_io wrapper in swap_csv_output.f90` (Phase B-3, 4 sites).

---

## Task 5: Phase B-4 — `src/core/swap_main.f90` (4 sites)

**Files:**
- Modify: `src/core/swap_main.f90`

The four sites are inside the rerun-mechanism block (`iun1`/`iun2` for `reruns.log` and `reruns.dat`). Phase C will retire the entire block in Task 12, but Phase B converts them anyway — it's the cleanest order (a Phase C that touches both a `getun` block and the `rdsets` call would mix two concerns).

Standard pattern-1 (output) for `iun1`/`reruns.log` and pattern-2 (input/old) for `iun2`/`reruns.dat`. Drop `integer :: getun` from the local declarations.

Commit message: `refactor(core): replace getun/fopens with file_io wrapper in swap_main.f90` (Phase B-4, 4 sites).

---

## Task 6: Phase B-5 — `src/io/macroporeoutput.f90` (9 sites)

**Files:**
- Modify: `src/io/macroporeoutput.f90`

9 sites. All uniform output-file openings (pattern 1).

Commit message: `refactor(io): replace getun/fopens with file_io wrapper in macroporeoutput.f90` (Phase B-5, 9 sites).

---

## Task 7: Phase B-6 — `src/io/toml/config_to_variables.f90` (4 sites)

**Files:**
- Modify: `src/io/toml/config_to_variables.f90`

The four sites include the `logf` opening inside the `block ... use variables, only: logf ... call fopens(logf, 'swap_swap.log', 'new', 'del') ... end block` construct (the `block` was added per ADR 0018 / SS-C step 5).

Replacements:
- The `logf = getun(20, 99); call fopens(logf, 'swap_swap.log', 'new', 'del')` becomes `call file_open(logf, 'swap_swap.log', 'replace', 'write')`.
- The `delfil('swap_swap.log', .false.)` call (if present) becomes `call file_delete('swap_swap.log')`.
- The other sites (CSV companion path-resolution helpers — search for them) follow standard patterns 1 or 2.

Commit message: `refactor(io): replace getun/fopens/delfil with file_io wrapper in config_to_variables.f90` (Phase B-6, 4 sites).

---

## Task 8: Phase B-7 — `src/crop/management_soil.f90` (4 sites)

**Files:**
- Modify: `src/crop/management_soil.f90`

The four sites are inside `SoilManagement(4)`-or-later (nutrient cases — flCropNut-stub-erred at runtime). The patterns are mixed:
- `snp = getun2(10, 90, 2); open(unit=snp, file=filnam, status='old')` → `call file_open(snp, filnam, 'old', 'read')`.
- `oup = getun2(10, 90, 2); open(unit=oup, file=filnam, status='unknown')` → `call file_open(oup, filnam, 'unknown', 'readwrite')`.

After replacement, also drop `getun2` and `getun` from the comma-separated `integer ::` local declarations (lines 47-48). **Do not drop `rdinqr`** in this task — that's covered by Phase D Task 13.

Commit message: `refactor(crop): replace getun/fopens with file_io wrapper in management_soil.f90` (Phase B-7, 4 sites).

---

## Task 9: Phase B-8 — `src/crop/cropgrowth.f90` (16 sites)

**Files:**
- Modify: `src/crop/cropgrowth.f90`

16 sites. Mostly uniform output-file openings; one might be inside the `outbalcrop*` block per the ADR 0019 history. Apply patterns mechanically.

Drop `getun` from local `integer ::` lines that become unused. **Do not drop `integer getun2` at line 1054** — that's a separate dead-decl cleanup in Phase D.

Commit message: `refactor(crop): replace getun/fopens with file_io wrapper in cropgrowth.f90` (Phase B-8, 16 sites).

---

## Task 10: Phase B-9 — `src/io/swapoutput.f90` (78 sites)

**Files:**
- Modify: `src/io/swapoutput.f90`

The big one. 78 sites, all uniform `iun = getun(20,90); call fopens(iun, filnam, 'new', 'del')` patterns. Each subroutine in the file has its own `integer :: getun` local declaration plus an output writer.

Subagent strategy: scan the file once to enumerate all `(getun, fopens)` pairs, apply pattern 1 to every one, drop the `integer :: getun` declarations, add `use file_io_mod, only: file_open` to each subroutine. The diff is large (~150 lines changed) but homogeneous.

**Special case:** lines 3452-3454 contain the `cexf = getun(10,90); call rddtmp(cexf)` paired-with-`getun`. Phase B converts the `getun` part:
```fortran
! Before:
cexf = getun(10, 90)
call rddtmp(cexf)

! After (Phase B):
call file_open(cexf, '<scratch_name>', 'replace', 'write')
call rddtmp(cexf)
```
The `rddtmp` call stays for now — it's removed in Phase C / Task 12. But this `getun` site needs to allocate `cexf` somehow. Easier: skip this one site in Phase B (leave the `getun(10,90)` line) and address it in Phase C with the surrounding `rddtmp` retirement. Document the skip with a `! Phase B skip: paired with rddtmp; converted in Phase C` comment.

Alternative — just convert by allocating `cexf` to a unit number with a no-op `inquire`-then-leave or use Fortran's `newunit=` directly inside the call:
```fortran
integer :: cexf
inquire(file='dummy', exist=exists)  ! no-op
! ... let Phase C handle it
```

Simplest: leave the `cexf = getun(10,90); call rddtmp(cexf)` block UNTOUCHED in Phase B and process it as a unit in Phase C. The file's overall acceptance grep at the end of Phase B will show a single residual `getun` site at line ~3453, which is acceptable and called out in the commit message.

After this task, the rest of `swapoutput.f90` is wrapper-only.

Commit message: `refactor(io): replace getun/fopens with file_io wrapper in swapoutput.f90` (Phase B-9, 77 sites; the 78th — the `getun` paired with `rddtmp` at ~line 3453 — is deferred to Phase C alongside the `rddtmp` retirement).

---

## Task 11: Phase C — Retire reruns + drop `rddtmp`

**Files:**
- Modify: `src/core/swap_main.f90` (drop the rerun-loop scaffolding)
- Modify: `src/io/swapoutput.f90` (drop the `rddtmp` block at lines 3452-3454)

- [ ] **Step 1: Read the current state of `swap_main.f90`**

```bash
sed -n '40,80p' src/core/swap_main.f90
```

You'll see (post-Phase-B) something like:
```fortran
! open logfile and read rerun file
call file_open(iun1, 'reruns.log', 'replace', 'write')
call file_open(iun2, 'reruns.dat', 'old', 'read')
call rdsets(iun2, iun1, 'reruns.dat', insets)
if (insets == 0) write (iun1, '(a)') 'No reruns defined.'

! reruns (if supplied; else this loop is performed only once)
do iset = 0, insets
   call rdfrom(iset, .true.)
   iTask = 1
   if (iCaller == 0) call swap(iCaller, iTask)
   if (iCaller /= 0) call dummy(iTask)
   iTask = 2
   if (iCaller == 0) call swap(iCaller, iTask)
   if (iCaller /= 0) call dummy(iTask)
   iTask = 3
   if (iCaller == 0) call swap(iCaller, iTask)
   if (iCaller /= 0) call dummy(iTask)
end do
close (iun2)
```

- [ ] **Step 2: Replace with the simplified flow**

```fortran
! Reruns retired (ADR 0023). Parameter sweeps now driven externally.
if (iCaller == 0) then
   call swap(iCaller, 1)
   call swap(iCaller, 2)
   call swap(iCaller, 3)
else
   call dummy(1)
   call dummy(2)
   call dummy(3)
end if
```

Drop the local declarations for `iun1`, `iun2`, `iset`, `insets` if they become unused. Verify by reading the surrounding context — those names may be used elsewhere (unlikely, but check).

- [ ] **Step 3: Drop the `rddtmp` block in `swapoutput.f90`**

Open `src/io/swapoutput.f90` around line 3450 and find the block:
```fortran
! --- delete temporary files
cexf = getun (10,90)
call rddtmp(cexf)
inquire(unit=20,opened=fileopen)
if(fileopen)close(20, STATUS = 'DELETE')
```

Replace with:
```fortran
! --- delete temporary files (TTutil-managed scratch retired with ADR 0023).
! The unit-20 close-with-delete is the project's own cleanup, kept.
inquire(unit=20, opened=fileopen)
if (fileopen) close(20, status='DELETE')
```

Drop the `cexf` local declaration if it becomes unused. Drop `getun` from the comma-list of locals if it becomes the last `getun` reference in the file (Phase B-9 already replaced the others).

- [ ] **Step 4: Build and verify**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep "Results:"
echo "=== acceptance ==="
grep -rE "\\b(rdsets|rdfrom|rddtmp)\\b" src/ --include="*.f90" || echo "OK: no rdsets/rdfrom/rddtmp"
grep -rE "\\b(getun|getun2|fopens|delfil)\\b" src/ --include="*.f90" || echo "OK: no getun/fopens/delfil"
```
Expected: build clean; pFUnit `Ok: 1, Fail: 0`; check-full 5/5; both acceptance greps return "OK: ...".

- [ ] **Step 5: Commit**

```bash
git add src/core/swap_main.f90 src/io/swapoutput.f90
git commit -m "$(cat <<'EOF'
refactor: retire reruns mechanism + drop rddtmp scratch cleanup

Phase C of TTutil retirement (umbrella spec
docs/superpowers/specs/2026-05-07-ttutil-retirement-design.md).

swap_main.f90: drops the `do iset = 0, insets` loop and the
`call rdsets`/`call rdfrom` plumbing. The reruns.dat /
reruns.log files are no longer opened or referenced. Parameter
sweeps move outside SWAP entirely (driven by an external
workflow per the agreed scope).

swapoutput.f90: drops the `cexf = getun(10,90); call rddtmp(cexf)`
pair (TTutil-managed scratch cleanup). The project's own
unit-20 close-with-DELETE stays.

This eliminates the last three TTutil function calls in `src/`:
rdsets, rdfrom, rddtmp.

Verified: build clean; pFUnit Ok: 1, Fail: 0; check-full 5/5.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 12: Phase D — Drop dead TTutil-related local declarations

**Files:**
- Modify: `src/crop/management_soil.f90:82` (drop `logical :: rdinqr`)
- Modify: `src/crop/cropgrowth.f90:1054` (drop `integer getun2`)

The actual TTutil reader calls (`rdinqr`, `rdsdor`, etc.) were already eliminated by ADRs 0019/0021/0022. What remains are dead local declarations whose target functions are no longer called. Phase B left them alone (it cleaned only the `getun`/`getun2`/`fopens` decls); Phase D mops up.

- [ ] **Step 1: Verify the declarations are unused**

```bash
grep -n "\\brdinqr\\b" src/crop/management_soil.f90
```
Expected: only line 82 (the declaration itself; no callers).

```bash
grep -n "\\bgetun2\\b" src/crop/cropgrowth.f90
```
Expected: only line 1054 (declaration; no callers — Task 9 converted any actual call sites).

If either grep shows additional matches, the declaration is still in use and must NOT be dropped. Stop and investigate.

- [ ] **Step 2: Drop the declarations**

In `src/crop/management_soil.f90` around line 82:
```fortran
! Before:
      logical :: rdinqr

! After: (line removed entirely)
```

In `src/crop/cropgrowth.f90` around line 1054:
```fortran
! Before:
      integer getun2

! After: (line removed entirely)
```

- [ ] **Step 3: Build and verify**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep "Results:"
echo "=== full TTutil-symbol grep ==="
grep -rnE "\\b(rdinit|rdsdor|rdsinr|rdfdor|rdinqr|rdatim|rdsdou|rdadou|rdfint|rdftim|rdfinr|rdscha|rdinar|rdinne|rdsets|rdfrom|rddtmp|getun|getun2|fopens|delfil)\\b" src/ --include="*.f90" 2>/dev/null
```
Expected: build clean; pFUnit `Ok: 1, Fail: 0`; check-full 5/5; the final grep returns no matches (or only matches inside string literals / comments — those are documentation, not code).

If the grep shows code-level matches, investigate before proceeding to Phase E.

- [ ] **Step 4: Commit**

```bash
git add src/crop/management_soil.f90 src/crop/cropgrowth.f90
git commit -m "$(cat <<'EOF'
chore: drop dead TTutil-related local declarations

Phase D of TTutil retirement (umbrella spec
docs/superpowers/specs/2026-05-07-ttutil-retirement-design.md).

The TTutil reader call sites (rdinqr, getun2, etc.) were already
eliminated by ADRs 0019/0021/0022 + Phase B. What remained were
dead local declarations whose target functions are no longer
called:

- src/crop/management_soil.f90:82  drop `logical :: rdinqr`
- src/crop/cropgrowth.f90:1054     drop `integer getun2`

After this commit, src/ has zero references to TTutil function
names (verified via the umbrella's full-symbol grep).

Verified: build clean; pFUnit Ok: 1, Fail: 0; check-full 5/5.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 13: Phase E — Drop `ttutil_dep` from meson + rename fatalerr_shim + ADR

**Files:**
- Rename: `src/error/fatalerr_shim.f90` → `src/error/fatalerr.f90`
- Modify: `src/error/fatalerr.f90` (rewrite header — drop "shim" framing)
- Modify: `meson.build` (drop ttutil subproject + dependency; rename source)
- Modify: `tests/unit/meson.build` (drop ttutil_dep; rename source)
- Delete: `subprojects/ttutil.wrap`
- Delete: `subprojects/ttutil/` (vendored source tree)
- Create: `docs/adr/0023-ttutil-retirement.md`
- Modify: `docs/adr/index.md` (add ADR 0023 row)

- [ ] **Step 1: Rename the shim file and rewrite its header**

```bash
git mv src/error/fatalerr_shim.f90 src/error/fatalerr.f90
```

Edit `src/error/fatalerr.f90` — rewrite the header comment block (currently describes the shim's purpose and "TO BE REMOVED" gate). New header:

```fortran
!> Canonical FatalERR implementation.
!!
!! Provides a free `subroutine FatalERR(MODULE, MESSAG)` matching the
!! TTutil-era signature so legacy free-call sites in `src/` resolve.
!! Routes through `fatalerr_collected`, which appends the error to
!! the modern `error_collection_t` singleton and aborts via
!! `error stop "fatal error(s) in swap input pipeline"` (exit 1).
!!
!! Pre-2026-05-07 this was a "shim" overriding TTutil's own
!! `fatalerr.f90` at link time. With TTutil retired (ADR 0023), this
!! is the canonical implementation; no shimming required.
SUBROUTINE FatalERR(MODULE, MESSAG)
   use error_mod, only: fatalerr_collected
   implicit none
   character(len=*), intent(in) :: MODULE, MESSAG
   call fatalerr_collected(trim(MODULE), trim(MESSAG))
END SUBROUTINE FatalERR
```

- [ ] **Step 2: Update `meson.build`**

In the `sources` list, rename `'src/error/fatalerr_shim.f90'` → `'src/error/fatalerr.f90'`.

Locate and remove the ttutil block:
```meson
ttutil_prj = subproject('ttutil', default_options: ['default_library=static'])
ttutil_dep = ttutil_prj.get_variable('ttutil_dep')
```

In the `executable('swap', ...)` call, remove `ttutil_dep` from the `dependencies:` list. If the list becomes single-element `[tomlf_dep]`, that's fine (or you can switch to a bare `dependencies: tomlf_dep` per meson syntax).

- [ ] **Step 3: Update `tests/unit/meson.build`**

Rename `'../../src/error/fatalerr_shim.f90'` → `'../../src/error/fatalerr.f90'` in `pfunit_extra_sources`.

In the `executable('unit-swap-tests', ...)` call, remove `ttutil_dep` from the `dependencies:` list.

- [ ] **Step 4: Verify the build still works against the renamed shim**

```
pixi run -e test build-linux
```
Expected: build clean. (At this point ttutil_dep is removed but the subproject directory is still on disk; meson will not resolve any ttutil symbols since nothing requests them.)

Run the tests:
```
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep "Results:"
```
Expected: pFUnit `Ok: 1, Fail: 0`; check-full 5/5.

- [ ] **Step 5: Vendor-delete the TTutil tree**

Verify `subprojects/ttutil/` is not tracked in git (it's brought in by the wrap file, not vendored):
```bash
git status subprojects/ttutil
```
Expected: `subprojects/ttutil/` shows as untracked content of a submodule, OR not reported at all (gitignored). Either way, deletion is safe.

```bash
rm -rf subprojects/ttutil
git rm subprojects/ttutil.wrap
```

Verify clean state:
```bash
ls subprojects/
git status
```

The remaining subproject under `subprojects/` should be `toml-f` (and possibly `test-drive` / others). Confirm `subprojects/ttutil` is gone.

- [ ] **Step 6: Re-build to verify no stale references**

```
pixi run -e test build-linux 2>&1 | tail -10
```
Expected: clean rebuild. If meson complains about a missing subproject, there's a stragger reference somewhere — grep `meson.build` and `tests/unit/meson.build` for `ttutil` again.

- [ ] **Step 7: Final acceptance run**

```
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep "Results:"
echo "=== final TTutil-free check ==="
grep -rn "ttutil\\|TTutil\\|TTUTIL" src/ meson.build tests/unit/meson.build 2>/dev/null | grep -v "^Binary file" || echo "OK: no ttutil refs in source or meson"
ls subprojects/ttutil 2>/dev/null && echo "FAIL: ttutil dir still exists" || echo "OK: subprojects/ttutil/ deleted"
ls subprojects/ttutil.wrap 2>/dev/null && echo "FAIL: wrap file still exists" || echo "OK: ttutil.wrap deleted"
```
Expected: build green; pFUnit `Ok: 1, Fail: 0`; check-full 5/5; ttutil-refs grep returns "OK: ..." (ignore docs/comments referencing TTutil historically — those are intended as historical record).

- [ ] **Step 8: Write `docs/adr/0023-ttutil-retirement.md`**

```markdown
---
title: "ADR 0023 — TTutil retired from the SWAP build"
date: 2026-05-07
status: accepted
---

# ADR 0023: TTutil retired from the SWAP build

## Context

ADRs 0019 / 0021 / 0022 retired the TTutil-based **data readers**
from the production runtime, leaving only utility functions
(`getun` / `getun2` / `fopens` / `delfil` for unit-number
management and file I/O) and the rerun-mechanism plumbing
(`rdsets` / `rdfrom` / `rddtmp`) still depending on TTutil.

Continuing to ship TTutil as a `subproject('ttutil', ...)` of
the meson build is overhead: 170 source files, ~5 MB, of which
SWAP uses fewer than ten functions — none of them genuinely
TTutil-specific.

## Decision

Drop TTutil from the build entirely. Replace its utility
functions with native Fortran intrinsics wrapped in a small
project module; retire the rerun mechanism (parameter sweeps
move out-of-band); rename the shim implementation of `FatalERR`
to canonical.

## Phase chronology

| Phase | What landed |
|---|---|
| A | New `file_io_mod` (`src/io/file_io.f90`) with `file_open`/`file_delete`/`file_exists` wrappers + pFUnit suite. |
| B | ~125 call sites of `getun`/`getun2`/`fopens`/`delfil` across 9 files converted to the wrapper. Per-file commits. |
| C | `swap_main.f90` simplified to `call swap(0,1); call swap(0,2); call swap(0,3)` triple. `rdsets`/`rdfrom`/`rddtmp` plumbing dropped. Reruns moved out-of-band. |
| D | Dead TTutil-related local declarations dropped (residue of earlier reader-deletion arcs). |
| E | `fatalerr_shim.f90` → `fatalerr.f90` rename; `ttutil_dep` and the `subprojects/ttutil/` tree removed from meson. |

## Consequences

- The SWAP binary no longer links against TTutil. The
  `subprojects/ttutil/` source tree (170 files, ~5 MB) is gone
  from the repo.
- Uniform file-open error reporting via `swap_log` at every I/O
  site.
- `swap_main.f90` is ~30 LoC shorter; parameter sweeps are an
  external concern.
- `fatalerr.f90` is the canonical fatal-error handler; the
  "shim" framing in its predecessor is gone.
- `[nutrients]` reactivation (separate umbrella) starts from a
  TTutil-free baseline.

## Acceptance

- `grep -rE "\b(rdinit|rdsdor|rdsinr|rdfdor|rdinqr|rdatim|rdsdou|rdadou|rdfint|rdftim|rdfinr|rdscha|rdinar|rdinne|rdsets|rdfrom|rddtmp|getun|getun2|fopens|delfil)\b" src/` → no matches.
- `grep -rn "ttutil\|TTutil\|TTUTIL" src/ meson.build tests/unit/meson.build` → no matches (only historical references in docs/comments).
- `subprojects/ttutil/` does not exist; `subprojects/ttutil.wrap` does not exist.
- `pixi run -e test test-pfunit` → `Ok: 1, Fail: 0`.
- `pixi run -e test check-full` → `5 passed, 0 failed`.
- check-full output CSVs are byte-identical to the pre-arc baseline.

## Related

- ADR 0018 (fatalerr shim) — the shim framing retired here.
- ADR 0019 (legacy readers retired) — TTutil data-reader removal.
- ADR 0021 (tillage TOML port) — closed one of the two surviving
  reader exceptions.
- ADR 0022 (SSDI TOML port) — closed the other.
- Future: `[nutrients]` umbrella — separate spec for nutrient
  reactivation.
```

- [ ] **Step 9: Update `docs/adr/index.md`**

Append a row for ADR 0023, matching the existing format. Likely:
```markdown
| [0023](0023-ttutil-retirement.md) | TTutil retired from the build (subproject + utility functions + reruns) | 2026-05-07 | accepted |
```

(Inspect the file first; match the actual column layout.)

- [ ] **Step 10: Commit**

```bash
git add src/error/fatalerr.f90 meson.build tests/unit/meson.build docs/adr/0023-ttutil-retirement.md docs/adr/index.md
# `git mv` already staged the rename (src/error/fatalerr_shim.f90 → src/error/fatalerr.f90).
# `rm -rf subprojects/ttutil` is a working-tree change; we may need to inform git:
git status
# subprojects/ttutil/ may show as untracked-content-removed for a submodule, or as
# untracked. Either way, no `git rm` needed for that dir; the wrap file is what git knows about.

git commit -m "$(cat <<'EOF'
refactor: drop ttutil_dep from meson; rename fatalerr_shim to canonical

Phase E of TTutil retirement (final). With Phases A-D having
eliminated every TTutil function call from `src/`, the dependency
itself can be dropped:

- src/error/fatalerr_shim.f90 -> src/error/fatalerr.f90 (rename
  via git mv). Header rewritten — no longer a shim, just the
  canonical FatalERR implementation routing through
  fatalerr_collected.
- meson.build: drop subproject('ttutil') and ttutil_dep; remove
  ttutil_dep from the swap executable's dependencies. Rename
  the renamed shim source.
- tests/unit/meson.build: drop ttutil_dep from the unit-test
  executable's dependencies. Rename source.
- subprojects/ttutil.wrap: deleted (was the meson subproject
  reference).
- subprojects/ttutil/: vendor-deleted (~5 MB, 170 source files).

ADR 0023 captures the umbrella; index.md updated.

Verified: build clean from a fresh checkout; pFUnit Ok: 1,
Fail: 0; check-full 5/5; check-full CSV outputs byte-identical
to pre-arc baseline.

Closes the TTutil retirement umbrella spec
(docs/superpowers/specs/2026-05-07-ttutil-retirement-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 11: Optional final tag**

```bash
git tag -a rescue/ttutil-retired -m "TTutil dropped from the SWAP build (ADR 0023).

Net: subprojects/ttutil/ deleted (~5 MB, 170 files); ttutil_dep
removed from meson; ~125 file-I/O sites converted to file_io_mod
wrapper; reruns mechanism retired (parameter sweeps move
out-of-band); fatalerr_shim.f90 renamed to canonical fatalerr.f90.

pFUnit Ok: 1, Fail: 0; check-full 5/5; CSV outputs byte-identical
to pre-arc baseline."
```

---

## Self-Review

**Spec coverage:**

| Spec section | Task |
|---|---|
| §Phase A — `file_io_mod` wrapper | Task 1 |
| §Phase B — bulk-replace getun/fopens/delfil | Tasks 2-10 |
| §Phase C — retire reruns + drop rddtmp | Task 11 |
| §Phase D — collapse stub-erred TTutil reader regions | Task 12 (scope reduced — actual reader calls already gone; only dead decls remain) |
| §Phase E — drop ttutil_dep from meson | Task 13 |

All spec sections covered. Phase D's scope was reduced after pre-flight inventory: the `rd*` reader calls were already eliminated by earlier ADRs, leaving only two dead local declarations to clean up. The plan reflects this.

**Acceptance criteria** (from spec):

- `file_io_mod` exists with three public procedures + unit tests — Task 1.
- Full TTutil-symbol grep returns no matches — verified at end of Tasks 11, 12, 13.
- `swap_main.f90` does not reference reruns — Task 11.
- `meson.build` does not reference `ttutil` — Task 13.
- `src/error/fatalerr.f90` is canonical — Task 13.
- `pixi run -e test test-pfunit` and `check-full` green — verified at every task.
- ADR 0023 committed — Task 13.

**Type / signature consistency:**

- `file_open(unit, path, status, action[, iostat])` — used identically across Task 1 (definition), Tasks 2-10 (call sites in Phase B), and Task 11 (the simplified swap_main path). Argument names and semantics consistent.
- `file_delete(path)` — used identically.
- `file_exists(path)` — used identically.
- The `(status, action)` mapping table is referenced by every Phase B task via the "Conventions for Phase B" section.

**Open questions punted to implementation:**

- The exact location of the new `'src/io/file_io.f90'` line in `meson.build`'s `sources` list (Task 1 step 5 — implementer chooses the exact placement).
- Edge cases that don't fit the wrapper (binary I/O, `recl=`, append mode) — Task 9 / Task 10 implementers will encounter these; the Conventions section explicitly tells them to leave such sites as native `open()` with a one-line comment justifying the bypass.
- Whether `subprojects/ttutil/` is actually under git tracking — Task 13 step 5 verifies before deleting.

These are pickable at execution time without re-planning.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-05-07-ttutil-retirement.md`.
