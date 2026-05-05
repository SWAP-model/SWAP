# Phase 4e — Unified Error Handling + Phase 4f Prep Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Unify `fatalerr` → `error_collection_t` across surviving files (skip readswap.f90 and crop sub-readers — those go in 4f/4g). Produce the Phase 4f adapter audit doc. Author macropore-case TOML so all 6 cases load via the new path before 4f rewires the main driver.

**Architecture:** No new modules; small extension to `src/error/error.f90` with a `fatalerr_collected` shim for deep physics call-chains. New helpers are test-only: `load_both_for_macroporeflow` in `parity_helpers.f90`, plus a `test_macroporeflow_parity.pf`.

**Tech Stack:** Unchanged from 4d.

Spec: `docs/superpowers/specs/2026-04-27-phase-4e-error-prep-design.md`

---

## Preamble: context every task needs

**Baseline:** Phase 4d complete at commit `b297116`, tag `rescue/phase-4d-remaining-configs`. 305 pFUnit tests, 6/6 regression green. main and development aligned. Submodule at `3ac2caf`.

**Branch discipline:**
- Work on `development`. No feature branches.
- One commit per task. Subject `<type>(<scope>): <what> (Phase 4e Task N)`.
- No pushes during 4e.
- Phase exit: fast-forward main to development locally; tag `rescue/phase-4e-error-prep`.

**Working directory:** `/home/zawadzkim/Code/swap` for all tasks unless noted (submodule work is inside `tests/swap-cases/`).

**Conventions** (unchanged from 4d).

**Verification after each task:**
```
pixi run -e test test-pfunit
./builddir/tests/unit/unit-swap-tests > /tmp/v.log 2>&1; \
   echo "exit=$?"; \
   grep -E "^[\.F ]+$|^\.+$" /tmp/v.log | tr -dc 'F' | wc -c
pixi run -e test check-fast
```
After each Part and at closeout:
```
pixi run -e test check-full
```

**Two concepts worth understanding before authoring tasks:**

1. **`fatalerr` is line-of-fire.** A `call fatalerr(routine, msg)` inside a deep physics subroutine writes the message to a log unit, prints to stderr, and triggers Fortran runtime abort via `read(5, *)` blocking on stdin. In a non-TTY context (pFUnit, CI), stdin returns EOF and the program dies with "Fortran runtime error: End of file." The replacement (`errors%append + abort_if_fatal`) preserves the abort but routes the message through `error_collection_t` which is testable, capturable, and (eventually) recoverable. The behavior change is small at the call site but big at the call hierarchy: the abort is now a checked control-flow point.

2. **Singleton vs threaded `errors`.** Threaded means: `subroutine X(args, errors); type(error_collection_t), intent(inout) :: errors`. Clean signature contract; clear caller responsibility. Singleton means: a module-level `type(error_collection_t), public :: global_errors` in `error_mod`; any code calls `global_errors%append(...)` without changing signatures. **D5 in the spec mandates singleton for physics paths and threaded for I/O / configuration paths.** The plan groups commits accordingly so each commit's reviewer sees a coherent slice of one strategy.

---

## File structure

### New

| File | Responsibility |
|---|---|
| `docs/phase-4f-config-to-variables-audit.md` | Part B output: `variables%` ↔ config map. |
| `tests/swap-cases/toml/3.macroporeflow/swap.toml` | Top-level macropore TOML (submodule). |
| `tests/swap-cases/toml/3.macroporeflow/swap.dra.toml` | Drainage TOML (submodule). |
| `tests/swap-cases/toml/3.macroporeflow/<crop>.crp.toml` | Crop TOML(s) (submodule). |
| `tests/unit/io/toml/test_macroporeflow_parity.pf` | Parity test mirroring 4d pattern. |
| `tests/regression/baselines/phase-4e-error-prep.log` | Closeout regression log. |

### Modify

| File | Change |
|---|---|
| `src/error/error.f90` | Add `fatalerr_collected` shim + module-level singleton (D5). |
| ~30 source files containing `call fatalerr` | Replace with `errors%append + abort_if_fatal()`. SKIP `src/io/readswap.f90`. SKIP crop sub-readers in `src/crop/cropgrowth.f90` (lines TBD per Task A2). |
| `tests/unit/error/test_error.pf` | Test for the new shim. |
| `tests/unit/io/toml/parity_helpers.f90` | Add `load_both_for_macroporeflow`. |
| `tests/unit/meson.build` | Register `test_macroporeflow_parity.pf`. |
| `tests/unit/testSuites.inc` | Register suite. |
| `docs/coverage-baseline.md` | Phase 4e section. |
| `docs/configuration-schema.md` | One-paragraph note pointing to audit doc. |

### Delete

None in 4e.

---

## Tasks

### Part A — Error replacement (singleton + threaded patterns)

The replacement is mechanical. We split it across multiple commits to keep diffs reviewable. **The order matters:** start with the `error_mod` extension (so callers have something to call), then the I/O paths (threaded), then physics paths (singleton via `fatalerr_collected` shim).

- [ ] **Task A1 — Extend `error_mod` with `fatalerr_collected` shim + global singleton.**
  Add to `src/error/error.f90`:
  - A module-level `type(error_collection_t), public, save :: global_errors`. Document its purpose: "physics-path drop-in for legacy `fatalerr` until the call sites are properly threaded with their own collections."
  - `subroutine fatalerr_collected(routine, message)` that calls `global_errors%append(ERR_FATAL, trim(message), trim(routine))` then `global_errors%abort_if_fatal()`.
  - Possibly add `ERR_FATAL` to the error code enum if it's not there.
  Add `tests/unit/error/test_error.pf` tests:
  - `test_fatalerr_collected_appends_and_aborts` — calls `fatalerr_collected('foo', 'bar')` after temporarily redirecting `global_errors`'s abort behavior (i.e., test the append step; abort can't be tested directly because it terminates the process — instead, test that the message lands in `global_errors` after a non-fatal append, demonstrating the channel works).
  - One pFUnit test verifying `global_errors%count() > 0` after a non-fatal-equivalent path.
  **Verify:** new tests pass; existing 305 tests still pass.
  **Commit:** `feat(error): add fatalerr_collected shim + module-level global_errors (Phase 4e Task A1)`.

- [ ] **Task A2 — Identify and document the SKIP regions in `cropgrowth.f90`.**
  This is a research/grep task with one tiny doc artifact. Run:
  ```
  grep -n "call fatalerr" src/crop/cropgrowth.f90
  ```
  Cross-reference each match against the file structure (`subroutine ...` headers around the match) to determine which calls are inside the crop sub-readers (`readwofost`, `readcropfixed`, `readgrass`-style sections — but those are actually in `readswap.f90`, NOT in `cropgrowth.f90`. Verify.). The 19 calls in `cropgrowth.f90` are likely all in the surviving physics path (CropGrowth task 1/2/3/4 dispatch), not in sub-readers. **Confirm and document.**
  Output: a short comment block at the top of `src/crop/cropgrowth.f90` listing which subroutine boundaries the sub-readers occupy (if any) and noting "all `call fatalerr` here is in surviving physics; safe to replace in 4e."
  No code change beyond the comment block. Commit alone or roll into A3.
  **Commit:** `docs(cropgrowth): document survive vs skip regions for fatalerr replacement (Phase 4e Task A2)`.

- [ ] **Task A3 — Replace fatalerr in I/O paths (threaded errors).**
  Files: `src/io/swapoutput.f90` (44), `src/io/swap_csv_output.f90` (11), `src/io/readmeteo.f90` (9), `src/io/macroporeoutput.f90` (1–4). Total ~70 calls.
  For each subroutine that calls `fatalerr`:
  1. Add `type(error_collection_t), intent(inout) :: errors` as the last argument.
  2. Replace `call fatalerr(routine, msg)` with `call errors%append(ERR_FATAL, msg, routine); call errors%abort_if_fatal()`.
  3. Update every caller of the modified subroutine to pass an `errors` from its scope. If the caller doesn't have one, propagate further up. The chain MUST end at a top-level entry point (e.g., `swap.f90`'s Initialize, or a public exported subroutine).
  4. Add `use error_mod, only: error_collection_t, ERR_FATAL` at the top of each modified file.
  **Risk:** the "chain ends at swap.f90" approach forces threading through MANY layers. If propagation is too painful for any I/O subroutine, fall back to the singleton (`call fatalerr_collected(routine, msg)` instead) and document the inconsistency.
  **Verify:** test-pfunit green; check-fast 4/4 green.
  **Commit:** `refactor(io): route fatalerr through error_collection_t in I/O paths (Phase 4e Task A3)`.

- [ ] **Task A4 — Replace fatalerr in physics paths (singleton via `fatalerr_collected`).**
  Files: everything under `src/atmosphere/`, `src/soil/`, `src/crop/`, `src/drainage/`, `src/heat/`, `src/macropore/`, `src/solute/`, `src/utils/` that contains `call fatalerr` AFTER excluding what A3 already did.
  Mechanical replace: `call fatalerr(routine, msg)` → `call fatalerr_collected(routine, msg)`. Add `use error_mod, only: fatalerr_collected` at the top of each modified file.
  No signature changes. No caller updates. The singleton handles the abort.
  Estimated ~80–100 calls across ~25 files.
  **Verify:** test-pfunit green; check-fast 4/4 green.
  **Commit:** `refactor(physics): route fatalerr through global_errors via fatalerr_collected (Phase 4e Task A4)`.

- [ ] **Task A5 — Replace fatalerr in `cropgrowth.f90` and `swap.f90`.**
  These are the "core" files with ~19 + a handful of calls. Use the threaded pattern (`errors` argument) since these are top-level dispatchers and threading is clean here. Per A2's audit, all `cropgrowth.f90` calls are in surviving physics — replace them all.
  In `swap.f90`, the few `call fatalerr` calls happen inside the DLL exchange handler — replace, threading `errors` from the top of `swap()`.
  **Verify:** test-pfunit + check-fast.
  **Commit:** `refactor(core,crop): route fatalerr through error_collection_t in core dispatchers (Phase 4e Task A5)`.

- [ ] **Task A6 — Part A check-full.**
  Run `pixi run -e test check-full`; expect 6/6 green. Save log to `tests/regression/baselines/phase-4e-task-a6-check-full.log`.
  **Commit:** `docs(baselines): record phase 4e Part A check-full output (Phase 4e Task A6)`.

### Part B — Adapter audit

- [ ] **Task B1 — Generate the variables list.**
  Research task. From `src/`, grep all `use variables, only: ...` lines and produce a deduplicated list of every field referenced. Also include bare `use variables` (no `, only:`) — those import the entire module; for those files, scan for actual usages of `variables%` or bare globals.
  Save the list to `/tmp/phase-4e-variables-list.txt`. This is the input for B2.
  **No commit.** This is interim data.

- [ ] **Task B2 — Author the audit doc.**
  Create `docs/phase-4f-config-to-variables-audit.md` with the markdown structure from spec D3:
  ```
  | variable | status (C/R/G) | source / target | notes |
  |---|---|---|---|
  ```
  For each variable from B1's list:
  - **C (covered):** find a `*_config_t` field that maps to it. Document the path (`config%simulation%tstart`).
  - **R (runtime):** field is set by simulation code, not input. Examples: `theta`, `h`, `t1900` (after TimeControl). Document briefly why it's runtime.
  - **G (gap):** field is read by execution paths and populated by legacy reader, but NO config path covers it. **Each G is a Phase 4f blocker** — add a "next step" note.
  **Important:** also produce a separate header section listing the legacy initialization order (which subroutines populate which fields, in `readswap.f90` and `cropgrowth.f90`). This is the order Phase 4f's `config_to_variables` adapter will need to preserve.
  Aim: 200–400 lines of audit doc, all variables classified.
  **Commit:** `docs(phase-4f): variables ↔ config audit (Phase 4e Task B2)`.

### Part C — Macropore TOML

- [ ] **Task C1 — Macropore audit.**
  Read `tests/swap-cases/3.macroporeflow/swap_linux.swp.template`, `swap.dra`, and any `.crp` files. Identify case-3-specific quirks. Particular interest:
  - Macropore physics section in the `.swp` (`SWMACRO=1` is the trigger; what other macropore-specific keys appear?).
  - Cross-reference against the schema (`src/config/swap_config.f90` and the typed configs) to find missing fields.
  Output: `docs/phase-4e-macroporeflow-audit.md` (~80–120 lines).
  **Commit:** `docs(phase-4e): macroporeflow legacy-file audit (Phase 4e Task C1)`.

- [ ] **Task C2 — Author macropore TOML files (in submodule).**
  Inside `tests/swap-cases/` submodule:
  - `toml/3.macroporeflow/swap.toml` — top-level
  - `toml/3.macroporeflow/swap.dra.toml` — drainage
  - `toml/3.macroporeflow/<crop>.crp.toml` for each rotation entry
  **If the audit (C1) found macropore-physics-specific keys not covered by the schema:** scope down — author what the schema covers, defer the rest with `# Phase 4f-prep:` comments. Don't author a `macropore_config_t` in this phase.
  Verify the smoke test `test_macroporeflow_loads_clean` (in `test_all_cases_smoke.pf`) passes.
  **DO NOT** commit the submodule from inside this task — leave the file additions in working tree; user will batch-commit at the end of Part C.

- [ ] **Task C3 — Macropore parity test.**
  Author `tests/unit/io/toml/test_macroporeflow_parity.pf` mirroring the 4d-extended pattern (7 sub-tests: general, simulation, meteorology, drainage, soil, crop_meta, plus 4d sections like bottom_boundary, heat, irrigation, solute, mowing/grazing, plus the macropore-specific assertions if any are schema-covered).
  Add `load_both_for_macroporeflow` to `parity_helpers.f90`. Register the test in `meson.build` + `testSuites.inc`.
  **Iterate:** test fails → fix TOML, validator, or reader. Each fix is a separate focused commit.
  **Verify:** test-pfunit green, no F's. test count grows by ~7.
  **Commit:** `test(parity): macroporeflow (case 3) full parity (Phase 4e Task C3)`.

- [ ] **User commits Task C2 submodule changes + outer-repo bump.**
  Standard pattern from prior phases:
  ```
  git -C /home/zawadzkim/Code/swap/tests/swap-cases add toml/3.macroporeflow/
  git -C /home/zawadzkim/Code/swap/tests/swap-cases commit -m "feat(toml/macroporeflow): swap.toml + swap.dra.toml + crp.toml (Phase 4e Task C2)"
  git -C /home/zawadzkim/Code/swap add tests/swap-cases
  git -C /home/zawadzkim/Code/swap commit -m "chore(swap-cases): bump submodule for macroporeflow TOML (Phase 4e Task C2)"
  ```

### Part D — Closeout

- [ ] **Task D1 — Schema doc note.**
  Add a 1-paragraph note at the top of `docs/configuration-schema.md` pointing to the new audit doc:
  > **Phase 4f preparation:** see `docs/phase-4f-config-to-variables-audit.md` for the field-by-field mapping between this schema and the legacy `variables` globals.
  **Commit:** `docs(schema): point to phase 4f audit doc (Phase 4e Task D1)`.

- [ ] **Task D2 — Coverage rebaseline.**
  Run `pixi run -e coverage coverage-report`. Update `docs/coverage-baseline.md` with a Phase 4e section. Numbers will likely be flat (the gcovr path-resolution issue still applies; error replacement doesn't move LoC by much).
  **Commit:** `docs(coverage): rebaseline at Phase 4e (Phase 4e Task D2)`.

- [ ] **Task D3 — Closeout: check-full + tag + main fast-forward.**
  - `pixi run -e test check-full` — confirm 6/6 green.
  - Save log to `tests/regression/baselines/phase-4e-error-prep.log`. Commit: `docs(baselines): record phase 4e check-full output (Phase 4e Task D3)`.
  - Tag: `git tag rescue/phase-4e-error-prep HEAD -m "Phase 4e closeout: error replacement + audit + macropore TOML"`.
  - Fast-forward main: `git checkout main && git merge --ff-only development && git checkout development`.
  - No push.

---

## Risks & mitigations

| Risk | Mitigation |
|---|---|
| Threading `errors` through I/O subroutines explodes the diff | Fall back to `fatalerr_collected` for any I/O subroutine where threading touches > 5 callers. Document the inconsistency in commit message. |
| `cropgrowth.f90` has hidden `fatalerr` calls inside crop sub-readers (not in `readswap.f90` as I claimed) | A2 verifies by structural grep before A3/A4 fire. If it does, A2's doc updates the SKIP list and A3/A4 work around it. |
| Audit (B2) surfaces large numbers of (G) gaps | Each G is a 4f blocker. Document them clearly; 4f's plan picks them up. If the count is large (>5), call a checkpoint and decide whether 4f scope grows or 4e adds another part to fix the gaps now. |
| Macropore TOML (C2) requires schema extensions | Spec D4 mandates the fallback: scope down to schema-covered fields, defer macropore-physics-specific keys. The goal is "case 3 loads via new path" not "case 3 has every field." |
| Singleton `global_errors` accumulates state across tests | pFUnit runs all tests in one process. Add a `call global_errors%clear()` to relevant test setup, or to `parity_helpers%reset_for_next_readswap`. |
| `abort_if_fatal()` semantics differ subtly from `fatalerr` (e.g., flushes log differently, doesn't prompt stdin) | This is intentional: the new behavior is correct. If a test relies on the old prompt, that test was wrong. |

---

## Verification gates

After Task A6 (Part A done): check-full 6/6 green. F count 0.
After Task B2 (Part B done): audit doc complete, every variable classified. Optionally: count of (G) entries reported.
After Task C3 (Part C done): macroporeflow parity test green. test-pfunit count grows by ~7. F count 0.
After Task D3 (closeout): tag exists; coverage rebaselined; main fast-forwarded.

Final audit-trail check:
```
git diff rescue/phase-4d-remaining-configs..HEAD -- src/ \
   | grep -v "fatalerr\|errors%append\|fatalerr_collected"
```
Output should show ONLY the `error_mod` extension (Task A1), nothing else under `src/`. Anything else means an accidental physics change — investigate.

## File count summary

- Modify: ~30 source files (error replacement) + 6 test/doc files.
- Create: ~6 files (audit doc, macropore audit, macropore parity test, baseline log, schema note, coverage update; macropore TOMLs are submodule-side).
- LoC delta: ~+1,000 / ~-500 (replacement is roughly 2:1 verbose; net +500).
