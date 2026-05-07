---
title: "TTutil retirement (ADR 0023 umbrella)"
author: Mateusz Zawadzki
date: 2026-05-07
status: draft
---

# TTutil retirement — design spec

Sequel to the legacy-readers-deletion arc (ADRs 0019, 0021, 0022).
Those arcs eliminated TTutil-based **data readers** from the
production runtime. This umbrella eliminates the remaining TTutil
**utility** dependencies (`getun`/`getun2`/`fopens`/`delfil` for
unit-number management and file I/O; `rdsets`/`rdfrom`/`rddtmp`
for the rerun mechanism) and drops the TTutil subproject from
the build entirely.

After this lands, the SWAP binary builds without TTutil. The
`subprojects/ttutil/` tree (~5 MB, 170 source files) becomes
unreferenced.

## Goal

Reach a state where:

1. `meson.build` does not reference `subproject('ttutil', ...)`,
   `ttutil_dep`, or pass it to any `executable()`.
2. `subprojects/ttutil/` and `subprojects/ttutil.wrap` (if
   present) are deleted from the working tree.
3. `grep -rE "\b(rdinit|rdsdor|rdsinr|rdfdor|rdinqr|rdatim|rdsdou|rdadou|rdfint|rdftim|rdfinr|rdscha|rdinar|rdinne|rdsets|rdfrom|rddtmp|getun|getun2|fopens|delfil)\b" src/` returns no matches.
4. A new module `file_io_mod` (in `src/io/file_io.f90`) wraps
   the project's file-open / -delete / -inquire idioms with
   uniform `swap_log` integration and is unit-tested.
5. `swap_main.f90` invokes the simulation directly without a
   rerun loop. Parameter sweeps are an externally-driven concern.
6. The `fatalerr_shim.f90` module is renamed to canonical
   `fatalerr.f90` (no longer a shim — there is nothing to shim
   against once TTutil is gone).
7. `pixi run -e test test-pfunit` and `pixi run -e test check-full`
   stay green at every phase boundary; the binary passes
   regression byte-identical to pre-arc state.

End-state: TTutil-free production binary; uniform file-open
logging at all I/O sites; `swap_main.f90` simplified by ~30 LoC.

## Non-goals

- **Nutrient reactivation.** The `flCropNut`-stub-erred
  `SoilManagement(2..7)` bodies and the `cropgrowth.f90`
  nutrient leftovers are collapsed in Phase D (mirrors how
  `irrigation(1)` / `SoilManagement(1)` were collapsed in
  SS-C of the legacy-readers-deletion arc). The physics
  returns as a fresh TOML port in a separate `[nutrients]`
  umbrella, not part of this arc.
- **Reimplementing the rerun mechanism.** Parameter sweeps move
  outside SWAP entirely; you drive them from a wrapper script.
  No internal `reruns.dat` parser.
- **Refactoring the `close()` semantics.** `close(unit)` calls
  in `src/` stay native — they're trivial and vary too much by
  use case (binary vs text, deletion-on-close, etc.) for a
  uniform wrapper to add value.
- **Touching `subprojects/toml-f/`** — the TOML parser is and
  remains a build dependency.
- **Adding logging at every successful file open.** Too noisy
  at ~120 sites. WARN-on-failure only.

## Phase chronology

| Phase | Title | Effort |
|---|---|---|
| A | `file_io_mod` wrapper introduction | S |
| B | Bulk-replace `getun`/`fopens`/`delfil` in active code | M |
| C | Retire `reruns` mechanism + drop `rddtmp` | S |
| D | Collapse stub-erred TTutil reader regions | S-M |
| E | Drop `ttutil_dep` from meson; rename fatalerr shim; vendor-delete | S |

Each phase is independently shippable and leaves the repo green.

## Phase A — `file_io_mod` wrapper

### Goal

Introduce a uniform open / delete / inquire wrapper used by all
production I/O. No call sites are converted yet — that's Phase B.

### Module location and contents

`src/io/file_io.f90`:

```fortran
module file_io_mod
   implicit none
   private

   public :: file_open      ! open file with logging + iostat; aborts on failure if iostat absent
   public :: file_delete    ! delete-if-exists; never errors on missing file
   public :: file_exists    ! inquire wrapper

contains

   !> Open a file, allocating a fresh unit number.
   !!
   !! Mirrors Fortran's intrinsic open() with the conventions used
   !! across this codebase. On failure: when `iostat` is provided,
   !! sets it and returns; when absent, logs at WARN level and
   !! routes through fatalerr_collected (exit 1).
   !!
   !! Logged on failure: resolved path, iostat code, cwd.
   !!
   !! status accepts: 'old', 'new', 'replace', 'unknown', 'scratch'.
   !! action accepts: 'read', 'write', 'readwrite'.
   subroutine file_open(unit, path, status, action, iostat)
      integer,           intent(out) :: unit
      character(len=*),  intent(in)  :: path
      character(len=*),  intent(in)  :: status
      character(len=*),  intent(in)  :: action
      integer, optional, intent(out) :: iostat
   end subroutine file_open

   !> Delete a file if it exists; no-op if it doesn't. Idempotent.
   subroutine file_delete(path)
      character(len=*), intent(in) :: path
   end subroutine file_delete

   !> Wraps inquire(file=path, exist=...). Cheap, side-effect-free.
   logical function file_exists(path) result(exists)
      character(len=*), intent(in) :: path
   end function file_exists

end module file_io_mod
```

### Replacement equivalences

The Phase B bulk replacement uses these mappings (the wrapper's
`status`/`action` semantics are native-Fortran, not legacy):

| Legacy idiom | Wrapper call |
|---|---|
| `iun = getun(lo, hi); call fopens(iun, name, 'new', 'del')` | `call file_open(iun, name, 'replace', 'write')` |
| `iun = getun(lo, hi); call fopens(iun, name, 'new', 'noop')` | `call file_open(iun, name, 'new', 'write')` |
| `iun = getun(lo, hi); call fopens(iun, name, 'old', 'noop')` | `call file_open(iun, name, 'old', 'read')` |
| `iun = getun(lo, hi); call fopens(iun, name, 'unknown', 'noop')` | `call file_open(iun, name, 'unknown', 'readwrite')` |
| `call delfil(name, .false.)` | `call file_delete(name)` |
| `inquire(file=path, exist=fl)` | `fl = file_exists(path)` |

Special-case sites that don't fit the wrapper (binary I/O,
`recl=`, `position='append'`) stay as native `open()` with an
inline comment justifying the bypass.

### Tests

`tests/unit/io/test_file_io.pf` — pFUnit suite covering:
- `file_open` succeeds for `status='replace'`/`'old'`/`'new'`.
- `file_open` with `iostat` arg returns nonzero on missing-file
  with `status='old'` (not abort).
- `file_open` without `iostat` aborts via `fatalerr_collected`
  on missing-file (test catches the error in the collection).
- `file_delete` is idempotent (no error on missing path).
- `file_exists` returns true / false correctly.
- Path-with-cwd logging on failure includes the resolved path.

### Meson wiring

- `meson.build`: add `'src/io/file_io.f90'` to `sources` (place
  near other `src/io/` files, before `swap.f90` which depends
  on later subsystems).
- `tests/unit/meson.build`: add `'../../src/io/file_io.f90'` to
  `pfunit_extra_sources`; register `test_file_io.pf` in
  `pf_files` and `testSuites.inc`.

### Acceptance for Phase A

- New module + tests committed.
- pFUnit Ok: 1, Fail: 0; tests count grows.
- check-full unchanged.
- No call sites converted yet (that's Phase B).

## Phase B — Bulk-replace `getun`/`fopens`/`delfil`

### Goal

Convert all ~120 call sites to use `file_io_mod`.

### Strategy

Per-file commits, in dependency order (leaves first). Each commit
edits one Fortran file (or a small group of related files) and
must leave the build green + tests green at boundaries.

### Inventory

From a full-tree grep at spec authoring time:
- `getun` — 86 calls
- `getun2` — 4 calls (semantically identical to getun for our purposes)
- `fopens` — 33 calls
- `delfil` — 1 call

Distribution across `src/` (subject to drift; re-survey at
implementation):
- `src/io/macroporeoutput.f90` — 4 sites (file headers + output writers)
- `src/io/swapoutput.f90` — many sites (the big output module)
- `src/io/swap_csv_output.f90`, `src/io/readmeteo.f90` — output + meteo
- `src/utils/sharedsimulation.f90` — shared-mode coordination file
- `src/soil/waterbalance.f90` — `dev_cmb` mass-balance writer
- `src/crop/management_soil.f90` — output sites in nutrient cases (will be deleted in Phase D, but Phase B converts them anyway for hygiene)
- `src/io/toml/config_to_variables.f90` — the `logf` opening
- `src/core/swap_main.f90` — `iun1`/`iun2` for reruns (will be deleted in Phase C)

### Edge cases

Sites that **don't** fit the wrapper and stay as native `open()`:
- Binary unformatted I/O (if any — verify by grep for `form='unformatted'`).
- Direct-access I/O (`access='direct'`, `recl=`).
- Append mode (`position='append'`).

Each such site keeps a one-line `! native open(): <reason>` comment.

### Verification at each commit

- `pixi run -e test build-linux` clean.
- `pixi run -e test test-pfunit` `Ok: 1, Fail: 0`.
- `pixi run -e test check-full` 5/5.

### Acceptance for Phase B

- `grep -rE "\b(getun|getun2|fopens|delfil)\b" src/` returns no
  matches outside `src/io/file_io.f90` itself (which doesn't
  reference TTutil — it implements the wrapper).
- check-full unchanged (byte-identical CSV outputs across
  regression cases).

## Phase C — Retire reruns; drop `rddtmp`

### Goal

Eliminate the `rdsets`/`rdfrom`/`rddtmp` calls. Parameter sweeps
move out of SWAP entirely.

### `swap_main.f90` simplification

Current shape (post-Phase-B; `getun`/`fopens` already replaced):

```fortran
call file_open(iun1, 'reruns.log', 'replace', 'write')
call file_open(iun2, 'reruns.dat', 'old', 'read')
call rdsets(iun2, iun1, 'reruns.dat', insets)
if (insets == 0) write(iun1, '(a)') 'No reruns defined.'

do iset = 0, insets
   call rdfrom(iset, .true.)
   iTask = 1; if (iCaller == 0) call swap(iCaller, iTask)
   ! ... iTask = 2, 3 ...
end do
close(iun2)
```

Becomes:

```fortran
call swap(iCaller, 1)
call swap(iCaller, 2)
call swap(iCaller, 3)
```

The `iun1`, `iun2`, `iset`, `insets` locals go away. The
`reruns.log` and `reruns.dat` files are no longer opened or
referenced anywhere. The `dummy(iTask)` paths for `iCaller /= 0`
remain inline if other tooling uses them; otherwise they go too
(verify by grep).

### `rddtmp` retirement

`rddtmp` (called once in `swapoutput.f90`) deletes a TTutil
scratch file. Investigate at implementation time:
- If the scratch file is created by TTutil internals only (which
  Phase B + reader retirement already eliminated), the call is
  dead; delete it.
- If the scratch file has SWAP-side users, replace `call rddtmp(...)`
  with `call file_delete(<scratch_name>)`.

### Tests

No new tests. The simplified `swap_main` is exercised by
check-full (every regression case is one rerun-set worth of
work). pFUnit suites are unaffected.

### Acceptance for Phase C

- `grep -rE "\b(rdsets|rdfrom|rddtmp)\b" src/` returns no matches.
- `swap_main.f90` is ~30 LoC shorter.
- check-full 5/5; outputs byte-identical to pre-arc.

## Phase D — Collapse stub-erred TTutil reader regions

### Goal

Eliminate the remaining `rd*` reader calls. They live in
`flCropNut=1`-stub-erred regions of `SoilManagement` and any
straggler in `cropgrowth.f90`.

### Inventory

At spec-authoring time, post-Phase-B:
- `rdsdor` (8), `rdsinr` (7), `rdfdor` (4), `rdinqr` (2),
  `rdatim` (2), `rdsdou` (1), `rdinit` (1) — all confirmed
  inside `SoilManagement(4..7)` and `cropgrowth.f90` nutrient
  remnants.

`flCropNut=1` is stub-erred at `tillage.f90:79`
(`call fatalerr_collected('DoTillage', 'flCropNut = 1 not (yet)
allowed')`). All `SoilManagement(*)` callers in `swap.f90` are
gated by `if (flCropNut)`; case (1) was collapsed to `return` in
SS-C step 3 (commit `67f2545`). Cases 2-7 retain physics that
runs only when nutrients are active — currently never.

### Approach

Mirror the SS-C step 3 / step 4 pattern: replace each case body
with `return`, leaving a comment pointing at the future
`[nutrients]` umbrella. The original physics is recoverable from
git when nutrients are reactivated.

```fortran
case (4)
   ! Daily nutrient balance. Body collapsed in TTutil retirement
   ! Phase D; restored from git when [nutrients] umbrella
   ! reactivates the subsystem (typed [nutrients] TOML block +
   ! adapter, parallel to ADR 0021/0022 patterns).
   return
```

Same for cases 5, 6, 7 of `SoilManagement`. Case 2 (state copy)
and case 3 (timed events) are checked individually — if they
contain no `rd*` calls, leave them; if they do, collapse.

`cropgrowth.f90` nutrient leftovers (any `rd*` reads still
present after SS-C step 2) get the same treatment.

### Imports cleanup

After body collapse, drop `swpfile`, `logf`, TTutil-related
helpers from any `use variables, only: ...` lines that no longer
reference them.

### Tests

No new tests. pFUnit suites for the now-empty cases were never
exercised (flCropNut=1 path was always stub-erred). check-full
unchanged.

### Acceptance for Phase D

- `grep -rE "\b(rdinit|rdsdor|rdsinr|rdfdor|rdinqr|rdatim|rdsdou|rdadou|rdfint|rdftim|rdfinr|rdscha|rdinar|rdinne)\b" src/` returns no matches.
- pFUnit unchanged; check-full 5/5.

## Phase E — Drop `ttutil_dep` from meson

### Goal

Remove TTutil from the build entirely; rename `fatalerr_shim`
to canonical `fatalerr`.

### Steps

1. **Rename `src/error/fatalerr_shim.f90` → `src/error/fatalerr.f90`**
   (or merge its contents into `error.f90` — decide at
   implementation time based on file-size considerations).
   Update header comment: drop the "shim" framing, document as
   the canonical implementation.
2. **Update meson sources:**
   - `meson.build`: rename `'src/error/fatalerr_shim.f90'` →
     `'src/error/fatalerr.f90'` in the `sources` list. Drop
     `subproject('ttutil', ...)`, `ttutil_dep`, and remove
     `ttutil_dep` from the `dependencies:` list of the `swap`
     executable.
   - `tests/unit/meson.build`: same renames; drop `ttutil_dep`
     from the test executable's dependencies.
3. **Vendor-delete `subprojects/ttutil/`** and (if present)
   `subprojects/ttutil.wrap`. The directory is currently
   untracked (per `.gitignore` or because subproject artefacts
   aren't versioned); verify by `git status` before deletion.
4. **Final acceptance grep:**
   - `grep -rn "ttutil" src/ meson.build tests/unit/meson.build`
     → no matches.
   - `grep -rn "subproject.*ttutil" .` → no matches.

### ADR

`docs/adr/0023-ttutil-retirement.md` (new). Single ADR covering
all five phases — mirrors ADR 0019's umbrella shape (which
covered SS-1..SS-11 of the legacy-readers retirement). Sections:
Context, Decision, Phases (A-E with brief summary each),
Consequences, Acceptance.

`docs/adr/index.md` updated.

### Tests

No new tests. pFUnit + check-full unchanged across this commit.
The CI matrix now builds without TTutil.

### Acceptance for Phase E

- `grep -rn "ttutil\|TTutil\|TTUTIL" src/ meson.build tests/unit/meson.build` → no matches except possibly in comments/docs that reference TTutil historically (those stay as historical record).
- `subprojects/ttutil/` does not exist on disk.
- Build clean from a fresh checkout (no stale .mod files relying on TTutil symbols).
- pFUnit Ok: 1, Fail: 0; check-full 5/5.
- ADR 0023 committed.

## Acceptance criteria (umbrella)

- [ ] `src/io/file_io.f90` exists with three public procedures
      (`file_open`, `file_delete`, `file_exists`) and unit tests.
- [ ] `grep -rE "\b(rdinit|rdsdor|rdsinr|rdfdor|rdinqr|rdatim|rdsdou|rdadou|rdfint|rdftim|rdfinr|rdscha|rdinar|rdinne|rdsets|rdfrom|rddtmp|getun|getun2|fopens|delfil)\b" src/` returns no matches.
- [ ] `swap_main.f90` does not reference reruns; calls `swap()`
      thrice without a loop.
- [ ] `meson.build` does not reference `ttutil`. `subprojects/ttutil/`
      does not exist.
- [ ] `src/error/fatalerr.f90` is the canonical fatal-error
      handler; the "shim" framing is gone.
- [ ] `pixi run -e test test-pfunit` exits 0; `pixi run -e test
      check-full` exits 0 with `5 passed, 0 failed`.
- [ ] check-full output CSVs are byte-identical to the pre-arc
      baseline for all five regression cases.
- [ ] ADR 0023 committed.

## Commit cadence (umbrella)

Per phase, ~1-3 commits each:

| Phase | Commits |
|---|---|
| A | 1 (module + tests + meson wiring) |
| B | ~10-15 (one per file or related-file group) |
| C | 2 (rerun loop simplification + rddtmp removal) |
| D | 2-3 (SoilManagement case collapse + cropgrowth leftover + imports cleanup) |
| E | 2-3 (rename fatalerr_shim, drop ttutil_dep + vendor-delete, ADR) |

Total: ~20-25 commits across the umbrella.

## Risk

- **Phase B regression risk.** Replacing 120+ I/O sites is the
  widest-touch phase. Mitigation: per-file commits, check-full
  green at each, byte-identical CSV outputs verified at the
  umbrella tail.
- **Phase C scope creep.** Verify before deleting that no
  external tooling (pyswap, internal scripts) reads
  `reruns.log`. If they do, document the deprecation in the
  Phase C commit and ADR 0023.
- **Phase D boundary precision.** Cases 2 and 3 of
  `SoilManagement` may not contain `rd*` calls at all (they're
  state-copy + timed-event dispatch). Per-case inspection at
  implementation time; collapse only the cases that actually
  require it. Avoid over-collapsing — every line we delete in
  Phase D is a line `[nutrients]` umbrella has to restore.
- **Phase E vendor-delete.** Confirm `subprojects/ttutil/` is
  not tracked in git before `rm -rf`. The meson subproject
  system pulls these from a wrap file; deleting the wrap file
  + the cached source tree is the right way.
- **`fatalerr` rename.** Modules that `use error_mod, only:
  fatalerr_collected` are unaffected (the public name doesn't
  change). The free `subroutine FatalERR(MODULE, MESSAG)` in
  the shim file (called via TTutil's old uppercase convention)
  is no longer needed once TTutil-based callers are gone — it
  can be deleted in Phase E along with the rename.

## Out of scope (reaffirmed)

- Nutrient reactivation — separate `[nutrients]` umbrella;
  schema lands at `[nutrients]` (top-level TOML section), with
  its own ADR (24+).
- Reimplementing reruns inside SWAP — externally-driven.
- Wrapping `close()` — trivial, varies by use case.
- Touching the toml-f subproject.
