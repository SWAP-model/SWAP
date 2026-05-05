# SS-4 — Macropore TOML stub-error Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reject `soil.swmacro = 1` at the TOML validation boundary with a clear fatal error referencing ADR 0010 / ADR 0011, instead of silently letting the modern binary attempt to run macropore physics with unpopulated globals.

**Architecture:** Add one cross-field stub-error block to `soil_config_validate` (mirroring the established pattern in `cropfixed_config_validate`). No new schema, no new field, no adapter changes. Macropore reader code in `readswap.f90` and the orphan `macropore_config_t` stay untouched. Update ADR 0010 to record the new TOML-boundary behaviour.

**Tech Stack:** Fortran 2008 (gfortran 13), pFUnit, pixi, meson. `error_collection_t` + `ERR_VALIDATION_CROSS_FIELD` from `src/error/error.f90`.

**Sub-spec parent:** `docs/superpowers/specs/2026-05-04-legacy-reader-retirement-design.md` (SS-4).

---

## File Structure

| File | Change |
|---|---|
| `src/config/soil_config.f90` | **Modify:** add a stub-error block at the top of `soil_config_validate` (after line 138, before the existing `check_int_enum(self%swmacro, ...)` call) that appends `ERR_VALIDATION_CROSS_FIELD` when `self%swmacro == 1`. |
| `tests/unit/config/test_soil_config.pf` | **Modify:** add one pFUnit test `test_soil_swmacro_one_stub_errors` asserting that `swmacro=1` produces a validation error and that the error message mentions ADR 0010 or ADR 0011. |
| `docs/adr/0010-macropore-deferral.md` | **Modify:** append a "TOML-boundary stub-error (added 2026-05-04)" section documenting the new validator behaviour and citing ADR 0011. |
| `docs/superpowers/specs/2026-05-04-legacy-reader-retirement-design.md` | **Modify:** flip the SS-4 row in the status table from OPEN to DONE; flip the matching row in the "Sub-spec roadmap" table. |

No other files. No physics edits. No schema additions.

---

## Conventions

1. **Test-first.** Write the pFUnit test before touching the validator. Verify it fails. Then add the validator change. Verify it passes.
2. **Stub-error message format** (mirrors `cropfixed_config.f90:151–195`): start with the offending key + value in parentheses, end with the ADR citation. Example: `'soil.swmacro=1 (macropore physics) not yet supported in the TOML pipeline; case 3 is excluded from regression per ADR 0011 and the macropore module remains deferred per ADR 0010.'`
3. **Verification cadence:** build → targeted pFUnit → full pFUnit → check-full. Each step is its own command; do not chain.
4. **Single commit** at the end (the change is small enough that gap-by-gap commits are not warranted).
5. **Case 3 expectation.** After this change, invoking `pixi run -e test regression macroporeflow` will fatal-error during `config%validate(...)`. This is the intended behaviour: case 3 is excluded from check-full per ADR 0011 and should not run via the modern binary. Confirm the failure mode reads cleanly (the user sees the stub-error message, not a downstream segfault).

---

## Task 1: Add the failing pFUnit test

**Files:**
- Modify: `tests/unit/config/test_soil_config.pf`
- Verification: `pixi run -e test test-pfunit -- --filter test_soil_swmacro_one_stub_errors`

- [ ] **Step 1: Append the new test to `tests/unit/config/test_soil_config.pf`**

Append the following test at the bottom of the file (immediately before any trailing module statement; if the file ends with `@test` blocks, just append after the last one):

```fortran
@test
subroutine test_soil_swmacro_one_stub_errors()
   use funit
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   use soil_config_mod, only: soil_config_t
   type(soil_config_t)      :: s
   type(error_collection_t) :: errors
   integer :: i
   logical :: found_macro_msg

   ! Defaults are valid (swmacro = 0). Flip the switch to 1 and assert
   ! validation rejects it with a CROSS_FIELD error whose message
   ! references soil.swmacro=1 and at least one of ADR 0010 / 0011.
   s%swmacro = 1
   call s%validate(errors)

   @assertTrue(errors%has_errors())

   found_macro_msg = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'soil.swmacro=1') > 0 .and. &
          (index(errors%items(i)%message, 'ADR 0010') > 0 .or. &
           index(errors%items(i)%message, 'ADR 0011') > 0)) then
         found_macro_msg = .true.
         exit
      end if
   end do
   @assertTrue(found_macro_msg)
end subroutine
```

This mirrors the inspection pattern at `tests/unit/config/test_drainage_config.pf:219–223` (the `altcu` cross-field test): iterate `1..errors%count()`, read `errors%items(i)%code` and `errors%items(i)%message`.

- [ ] **Step 2: Build to make sure the test compiles**

Run: `pixi run -e test build-linux 2>&1 | tail -20`
Expected: build succeeds. If the test fails to compile, the most likely cause is a `use` statement: confirm `error_collection_t` and `ERR_VALIDATION_CROSS_FIELD` come from `error_mod`, and `soil_config_t` from `soil_config_mod` (matches the pattern in existing tests).

- [ ] **Step 3: Run the test to verify it FAILS**

Run: `pixi run -e test test-pfunit -- --filter test_soil_swmacro_one_stub_errors 2>&1 | tail -20`
Expected: 1 test, 1 failure. The failure should be on the `@assertTrue(errors%has_errors())` line, because the current `soil_config_validate` accepts `swmacro = 1` (it's a valid enum value).

If the test passes at this point, something is wrong — investigate before proceeding.

---

## Task 2: Add the stub-error block

**Files:**
- Modify: `src/config/soil_config.f90`
- Verification: `pixi run -e test test-pfunit -- --filter test_soil`

- [ ] **Step 1: Add the stub-error block at the top of `soil_config_validate`**

In `src/config/soil_config.f90`, locate `subroutine soil_config_validate` (currently at line 133). Find the existing first executable line (`call check_int_enum(self%swsophy, [0, 1], "soil.swsophy", errors)` near line 140). Insert the following block immediately before that first `check_int_enum` call:

```fortran
      ! ----- Phase 4f-extend SS-4 stub-error: macropore deferred -----
      ! ADR 0010 keeps macropore_config_t as orphan infrastructure and
      ! the macropore reader code in readswap.f90 unwired. ADR 0011
      ! excludes case 3 (3.macroporeflow) from regression. Authoring
      ! swmacro=1 in a TOML config would let the modern binary copy the
      ! switch into the legacy global without populating any macropore
      ! state, producing silent runtime corruption. Reject it here with
      ! a clear message instead.
      if (self%swmacro == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'soil.swmacro=1 (macropore physics) not yet supported in ' // &
            'the TOML pipeline; case 3 is excluded from regression per ' // &
            'ADR 0011 and the macropore module remains deferred per ' // &
            'ADR 0010.', 'soil')
      end if

```

The trailing blank line keeps it visually separated from the existing enum checks below.

- [ ] **Step 2: Verify the `use error_mod` import already includes `ERR_VALIDATION_CROSS_FIELD`**

Run: `grep -n "use error_mod" src/config/soil_config.f90 | head -5`
Expected: an existing `use error_mod, only: ...` line near the top of the module. Read it. If it does NOT already include `ERR_VALIDATION_CROSS_FIELD`, add it to the `only:` list:

```fortran
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD, ...
```

(Append to whatever the existing list is — do not replace it.) If the module uses bare `use error_mod` without `only`, no change needed.

- [ ] **Step 3: Build**

Run: `pixi run -e test build-linux 2>&1 | tail -10`
Expected: build succeeds (no compile errors).

- [ ] **Step 4: Run the targeted test to verify it PASSES**

Run: `pixi run -e test test-pfunit -- --filter test_soil_swmacro_one_stub_errors 2>&1 | tail -20`
Expected: 1 test, 0 failures.

- [ ] **Step 5: Run the full `test_soil` suite to ensure no other test regresses**

Run: `pixi run -e test test-pfunit -- --filter test_soil 2>&1 | tail -30`
Expected: all `test_soil_*` tests pass. (None of the existing tests in `test_soil_config.pf` set `swmacro = 1`, so none should regress.)

---

## Task 3: Full-suite verification

**Files:** none (verification only)

- [ ] **Step 1: Run the full pFUnit suite**

Run: `pixi run -e test test-pfunit 2>&1 | tail -15`
Expected: `Ok: 1` (the `unit-swap-tests` suite, all green).

If anything fails, investigate. The most likely failure mode would be an unrelated test that incidentally constructs a `soil_config_t` with `swmacro=1` for some other reason (e.g. an integration test fixture). Search for such constructions:

```bash
grep -rn "swmacro" tests/unit/ src/io/toml/
```

Fix any test that erroneously sets `swmacro = 1` by changing it to `0` (or by clearing the error collection where macropore is not the assertion target).

- [ ] **Step 2: Run check-full to confirm the 5 non-macropore cases stay green**

Run: `pixi run -e test check-full 2>&1 | tail -15`
Expected:
```
✓ hupselbrook: regression ok
✓ grassgrowth: regression ok
✓ oxygenstress: regression ok
✓ salinitystress: regression ok
✓ surfacewater: regression ok
Results: 5 passed, 0 failed
```

(Case 3 / macroporeflow is excluded by `tests/regression/test_output_regression.py` per ADR 0011; check-full does not invoke it.)

- [ ] **Step 3: Sanity-check the case-3 failure mode**

Run: `pixi run -e test regression macroporeflow 2>&1 | tail -25`
Expected: the run aborts with a fatal error mentioning `soil.swmacro=1` and citing ADR 0010 / 0011. The exact wording should match the validator message added in Task 2.

If the run instead segfaults, hangs, or produces an unrelated error, return to Task 2 to verify the validator wiring.

If the run succeeds (i.e. case 3 produces output despite swmacro=1), the stub-error did not fire — investigate whether `soil_config_validate` is actually being called from the validation pipeline (`grep -n "soil%validate" src/config/swap_config.f90`).

---

## Task 4: Update ADR 0010

**Files:**
- Modify: `docs/adr/0010-macropore-deferral.md`

- [ ] **Step 1: Append a new section to `docs/adr/0010-macropore-deferral.md`**

At the bottom of the file (after the existing content; do not modify the existing decision record), append:

```markdown

## Update 2026-05-04 — TOML-boundary stub-error (Phase 4f-extend SS-4)

Authoring `soil.swmacro = 1` in a TOML configuration is now rejected at
validation time by `soil_config_validate` (in `src/config/soil_config.f90`).
The error message reads:

> soil.swmacro=1 (macropore physics) not yet supported in the TOML
> pipeline; case 3 is excluded from regression per ADR 0011 and the
> macropore module remains deferred per ADR 0010.

Rationale: the modern binary's adapter copies `swmacro` into the legacy
global without populating any other macropore state. Without the
stub-error, a user who authored `swmacro = 1` would get silent runtime
corruption (or a crash deep in macropore physics that reads
unallocated arrays) instead of an immediate, actionable failure.

This update does not change the deferral itself: `macropore_config_t`
remains orphan infrastructure, no `read_macropore_toml` module exists,
and no `[macropore]` field is wired into `swap_config_t`. The
stub-error is the TOML-side counterpart to the regression exclusion
recorded in ADR 0011 — both close the macropore path cleanly without
removing the legacy code.

When future macropore work re-enables the module, this stub-error
must be removed in the same change that wires `[macropore]` into the
TOML pipeline.
```

- [ ] **Step 2: Verify the ADR file still renders cleanly**

Run: `head -20 docs/adr/0010-macropore-deferral.md && echo --- && tail -25 docs/adr/0010-macropore-deferral.md`
Expected: the original front-matter is untouched at the top; the new section appears at the bottom; no markdown syntax errors.

---

## Task 5: Update the umbrella spec status

**Files:**
- Modify: `docs/superpowers/specs/2026-05-04-legacy-reader-retirement-design.md`

- [ ] **Step 1: Mark SS-4 as DONE in the status snapshot table**

In `docs/superpowers/specs/2026-05-04-legacy-reader-retirement-design.md`, locate the status snapshot table (under "## Status snapshot (2026-05-04)"). Find the row:

```
| Macropore TOML stub-error | OPEN | `swmacro` currently passes validation silently |
```

Change it to:

```
| Macropore TOML stub-error | DONE | SS-4 (`6ab2b76`+); `swmacro=1` rejected by `soil_config_validate` |
```

(The commit SHA in parentheses will be the SS-4 commit hash from Task 6 — leave it as a placeholder during this step and update it in Task 6 after committing.)

- [ ] **Step 2: Mark SS-4 as DONE in the sub-spec roadmap table**

In the same file, locate the "## Sub-spec roadmap" section and the row:

```
| SS-4 | Macropore TOML stub-error | OPEN | new spec | — |
```

Change it to:

```
| SS-4 | Macropore TOML stub-error | DONE | `2026-05-04-ss4-macropore-toml-stub-error.md` | — |
```

---

## Task 6: Commit

**Files:** all of the above.

- [ ] **Step 1: Review the staged diff**

Run:

```bash
cd /home/zawadzkim/Code/swap
git status --short
git diff --stat
```

Expected files modified:
- `src/config/soil_config.f90`
- `tests/unit/config/test_soil_config.pf`
- `docs/adr/0010-macropore-deferral.md`
- `docs/superpowers/specs/2026-05-04-legacy-reader-retirement-design.md`

And one new file:
- `docs/superpowers/plans/2026-05-04-ss4-macropore-toml-stub-error.md` (this plan itself; commit it as part of the work)

If anything else shows up, investigate before committing.

- [ ] **Step 2: Stage and commit**

```bash
cd /home/zawadzkim/Code/swap
git add src/config/soil_config.f90 \
        tests/unit/config/test_soil_config.pf \
        docs/adr/0010-macropore-deferral.md \
        docs/superpowers/specs/2026-05-04-legacy-reader-retirement-design.md \
        docs/superpowers/plans/2026-05-04-ss4-macropore-toml-stub-error.md
git commit -m "$(cat <<'EOF'
feat(config): stub-error soil.swmacro=1 at TOML boundary (SS-4)

Phase 4f-extend SS-4. Reject swmacro=1 in soil_config_validate with a
clear message citing ADR 0010 (macropore deferral) and ADR 0011
(macropore regression exclusion), instead of letting the modern
binary copy the switch into the legacy global without populating any
macropore state.

- src/config/soil_config.f90 — new stub-error block at the top of
  soil_config_validate; mirrors the cropfixed_config pattern.
- tests/unit/config/test_soil_config.pf — new test asserting the
  stub-error fires and the message references the relevant ADRs.
- docs/adr/0010-macropore-deferral.md — appended a section
  documenting the new TOML-boundary behaviour.
- docs/superpowers/specs/2026-05-04-legacy-reader-retirement-design.md
  — flipped SS-4 status from OPEN to DONE.
- docs/superpowers/plans/2026-05-04-ss4-macropore-toml-stub-error.md
  — committed alongside the implementation.

Verified: pFUnit suite green; check-full 5/5 green; case-3 invocation
fails fast with the new validator message.
EOF
)"
```

- [ ] **Step 3: Update the commit SHA placeholder in the umbrella spec**

After committing, get the commit hash:

```bash
git log -1 --format=%h
```

Then update the umbrella spec's status snapshot row (changed in Task 5 Step 1) with the actual SHA. If the placeholder text reads `(SS-4 (`6ab2b76`+); ...)`, replace `6ab2b76` with the new hash from Task 6 Step 2. Then amend:

```bash
git add docs/superpowers/specs/2026-05-04-legacy-reader-retirement-design.md
git commit --amend --no-edit
```

(`--amend` is fine here because the commit has not been pushed and the only change is a doc-string SHA update.)

- [ ] **Step 4: Final verification**

Run:

```bash
git log -1 --stat
pixi run -e test check-full 2>&1 | tail -10
```

Expected: 5 cases passed, the commit shows the expected file list.

---

## Definition of done

- `soil_config_validate` rejects `swmacro = 1` with a fatal error citing ADR 0010 / 0011.
- `test_soil_swmacro_one_stub_errors` passes.
- Full pFUnit suite green.
- check-full green (5/5 non-macropore cases).
- `pixi run -e test regression macroporeflow` fails fast with the new validator message (intentional).
- ADR 0010 records the new TOML-boundary stub-error behaviour.
- Umbrella spec SS-4 status flipped to DONE with commit SHA.
- One commit with the message above.
