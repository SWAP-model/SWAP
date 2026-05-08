# [nutrients] N3 — Runtime Activation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Open the `flCropNut` runtime gate on the TOML build by wiring the legacy global from per-rotation typed config and removing both holdout stub-errors. Closes the `[nutrients]` umbrella.

**Architecture:** Per-rotation assignment inside `cropwofost_init_from_config` (Option A from the spec). All nutrient inputs already flow through N1/N2a/N2b adapters; this arc just opens the runtime gate. Zero new TTutil/dtutil dependencies introduced.

**Tech Stack:** Fortran 2008, gfortran, meson+ninja+pixi, pFUnit (test-pfunit), check-full byte-identical regression harness.

**Spec:** `docs/superpowers/specs/2026-05-08-nutrients-N3-runtime-activation-design.md`

---

## File Structure

**Modified:**
- `src/crop/cropwofost_init.f90` — remove validator stub-error, add `flCropNut` to `use variables, only:` list, assign `flCropNut = cfg%nutrient%flcropnut` per rotation
- `src/crop/tillage.f90` — remove runtime stub-error
- `src/core/swap.f90` — comment cleanup
- `src/crop/management_soil.f90` — comment cleanup
- `src/crop/cropgrowth.f90` — comment cleanup
- `tests/unit/crop/test_cropwofost_init.pf` — add two pFUnit tests for per-rotation wiring

**Created:**
- `docs/adr/0028-nutrients-N3-runtime-activation.md`

**Index updated:**
- `docs/adr/index.md` — add 0028 entry

**No changes to:**
- `testSuites.inc` (suite already registered)
- Any meson.build files
- `apply_cropwofost_nutrient`, `apply_nutrients`, `apply_nutrients_events` (already complete from N1/N2a/N2b)

---

### Task 1: Lift validator stub-error in cropwofost_init.f90

This unblocks `flcropnut=true` from reaching the runtime. After this task, calling `cropwofost_init_from_config` with `cfg%nutrient%flcropnut = .true.` no longer aborts — but `flCropNut` still stays `.false.` (wiring comes in Task 3). All existing tests must still pass because no current case sets `flcropnut=true`.

**Files:**
- Modify: `src/crop/cropwofost_init.f90:85-87`

- [ ] **Step 1: Delete the validator stub-error block**

Remove lines 85-87 of `src/crop/cropwofost_init.f90`:

```fortran
      if (cfg%nutrient%flcropnut) &
         call fatalerr_collected('cropwofost_init', &
            'flcropnut=.true. not supported on TOML path; validator should have rejected.')
```

The surrounding stub-error block (swsoybean, swbulb, swco2, etc.) stays intact; only this one triplet goes.

- [ ] **Step 2: Build and run check-full**

Run: `pixi run check-full`
Expected: `Results: 5 passed, 0 failed`. Existing cases don't set `flcropnut=true`, so this deletion is invisible to them.

- [ ] **Step 3: Run pFUnit suite**

Run: `pixi run -e test test-pfunit 2>&1 | tail -20`
Expected: `PASSED ... TESTS COMPLETE` line, suite count unchanged from current 607 passed / 1 disabled (or whatever the current count is — confirm no decrease and zero failures).

- [ ] **Step 4: Commit**

```bash
git add src/crop/cropwofost_init.f90
git commit -m "$(cat <<'EOF'
refactor(cropwofost): SS-N3 Task 1 — lift flcropnut validator stub-error

Drops the cropwofost_init_from_config defense-in-depth guard that
fatal-errored when cfg%nutrient%flcropnut = .true.. The runtime stub
at tillage.f90:73 still gates nutrient physics; that comes down in
Task 4. flCropNut wiring lands in Task 3.

Spec: docs/superpowers/specs/2026-05-08-nutrients-N3-runtime-activation-design.md
EOF
)"
```

---

### Task 2: Add failing pFUnit tests for per-rotation flCropNut wiring

Two test routines: one confirms `flcropnut=true` propagates to `flCropNut=.true.`, the other confirms `flcropnut=false` resets `flCropNut` to `.false.` (the per-rotation toggle case). Both will fail initially because Task 3 hasn't been done yet.

The `test_cropwofost_init_suite` is already registered in `tests/unit/testSuites.inc:60` — no suite-registration edit needed.

**Files:**
- Modify: `tests/unit/crop/test_cropwofost_init.pf` (append two new `@test` routines)

- [ ] **Step 1: Append the first test routine**

Append at the end of `tests/unit/crop/test_cropwofost_init.pf`:

```fortran

@test
subroutine test_flcropnut_true_propagates_to_global()
   ! N3: cfg%nutrient%flcropnut = .true. on a rotation must set the
   ! legacy global flCropNut = .true. so the daily loop's gates open.
   use funit
   use iso_fortran_env, only: real64
   use cropwofost_config_mod, only: cropwofost_config_t
   use cropwofost_init_mod,   only: cropwofost_init_from_config
   use variables, only: flCropNut
   type(cropwofost_config_t) :: cfg
   real(real64) :: FraDeceasedLvToSoil
   integer, parameter :: ICROP = 1

   ! Build the minimum viable fixture (same shape as test_cropwofost_init_basic_globals).
   cfg%crop_factor%swcf         = 2
   cfg%crop_factor%albedo       = 0.19_real64
   cfg%crop_factor%rsc          = 207.0_real64
   cfg%crop_factor%rsw          = 0.0_real64
   cfg%oxygen_stress%swoxygen   = 1
   cfg%oxygen_stress%swwrtnonox = 1
   cfg%oxygen_stress%aeratecrit = 0.5_real64
   cfg%oxygen_stress%hlim1      = -10.0_real64
   cfg%oxygen_stress%hlim2u     = -25.0_real64
   cfg%oxygen_stress%hlim2l     = -25.0_real64
   cfg%drought_stress%swdrought = 1
   cfg%drought_stress%hlim3h    = -300.0_real64
   cfg%drought_stress%hlim3l    = -500.0_real64
   cfg%drought_stress%hlim4     = -10000.0_real64
   cfg%drought_stress%adcrh     = 0.5_real64
   cfg%drought_stress%adcrl     = 0.1_real64
   cfg%salinity%swsalinity      = 1
   cfg%salinity%saltmax         = 0.732_real64
   cfg%salinity%saltslope       = 0.0868_real64
   cfg%compensate%swcompensate  = 0
   cfg%interception%swinter     = 1
   cfg%interception%cofab       = 0.25_real64
   cfg%co2%swco2                = 0
   cfg%root%swrd                = 2
   cfg%root%rdi                 = 10.0_real64
   cfg%root%rri                 = 1.2_real64
   cfg%root%rdc                 = 50.0_real64
   cfg%root%swdmi2rd            = 1
   cfg%root%swrdc               = 0
   cfg%management%swpotrelmf        = 2
   cfg%management%relmf             = 0.8_real64
   cfg%management%fradeceasedlvtosoil = 0.35_real64
   cfg%soybean%swsoybean            = 0
   cfg%bulb%swbulb                  = 0
   cfg%schedule%schedule            = 0
   allocate(cfg%root%rdctb(2, 2))
   cfg%root%rdctb(1,1) = 0.0_real64; cfg%root%rdctb(1,2) = 1.0_real64
   cfg%root%rdctb(2,1) = 1.0_real64; cfg%root%rdctb(2,2) = 0.0_real64

   ! Nutrient activation. nmxlv must be allocated and non-empty when
   ! flcropnut=true (validator rule, cropwofost_config.f90:1012-1019).
   cfg%nutrient%flcropnut = .true.
   allocate(cfg%nutrient%nmxlv(1))
   cfg%nutrient%nmxlv(1) = 60.0_real64

   ! Pre-set the global to .false. to prove the assignment is doing the work.
   flCropNut = .false.

   call cropwofost_init_from_config(cfg, ICROP, FraDeceasedLvToSoil)

   @assertTrue(flCropNut)
end subroutine test_flcropnut_true_propagates_to_global


@test
subroutine test_flcropnut_false_resets_global()
   ! N3 per-rotation toggle: a rotation with flcropnut=.false. must
   ! drive the global to .false. even if a previous rotation left it
   ! .true.. cropwofost_init_from_config runs at every rotation start.
   use funit
   use iso_fortran_env, only: real64
   use cropwofost_config_mod, only: cropwofost_config_t
   use cropwofost_init_mod,   only: cropwofost_init_from_config
   use variables, only: flCropNut
   type(cropwofost_config_t) :: cfg
   real(real64) :: FraDeceasedLvToSoil
   integer, parameter :: ICROP = 1

   cfg%crop_factor%swcf         = 2
   cfg%crop_factor%albedo       = 0.19_real64
   cfg%crop_factor%rsc          = 207.0_real64
   cfg%crop_factor%rsw          = 0.0_real64
   cfg%oxygen_stress%swoxygen   = 1
   cfg%oxygen_stress%swwrtnonox = 1
   cfg%oxygen_stress%aeratecrit = 0.5_real64
   cfg%oxygen_stress%hlim1      = -10.0_real64
   cfg%oxygen_stress%hlim2u     = -25.0_real64
   cfg%oxygen_stress%hlim2l     = -25.0_real64
   cfg%drought_stress%swdrought = 1
   cfg%drought_stress%hlim3h    = -300.0_real64
   cfg%drought_stress%hlim3l    = -500.0_real64
   cfg%drought_stress%hlim4     = -10000.0_real64
   cfg%drought_stress%adcrh     = 0.5_real64
   cfg%drought_stress%adcrl     = 0.1_real64
   cfg%salinity%swsalinity      = 1
   cfg%salinity%saltmax         = 0.732_real64
   cfg%salinity%saltslope       = 0.0868_real64
   cfg%compensate%swcompensate  = 0
   cfg%interception%swinter     = 1
   cfg%interception%cofab       = 0.25_real64
   cfg%co2%swco2                = 0
   cfg%root%swrd                = 2
   cfg%root%rdi                 = 10.0_real64
   cfg%root%rri                 = 1.2_real64
   cfg%root%rdc                 = 50.0_real64
   cfg%root%swdmi2rd            = 1
   cfg%root%swrdc               = 0
   cfg%management%swpotrelmf        = 2
   cfg%management%relmf             = 0.8_real64
   cfg%management%fradeceasedlvtosoil = 0.35_real64
   cfg%soybean%swsoybean            = 0
   cfg%bulb%swbulb                  = 0
   cfg%nutrient%flcropnut           = .false.
   cfg%schedule%schedule            = 0
   allocate(cfg%root%rdctb(2, 2))
   cfg%root%rdctb(1,1) = 0.0_real64; cfg%root%rdctb(1,2) = 1.0_real64
   cfg%root%rdctb(2,1) = 1.0_real64; cfg%root%rdctb(2,2) = 0.0_real64

   ! Pre-set the global to .true. — simulates a previous rotation that had
   ! nutrients on. The current rotation must reset it.
   flCropNut = .true.

   call cropwofost_init_from_config(cfg, ICROP, FraDeceasedLvToSoil)

   @assertFalse(flCropNut)
end subroutine test_flcropnut_false_resets_global
```

- [ ] **Step 2: Run the new tests to confirm they fail**

Run: `pixi run -e test test-pfunit 2>&1 | grep -E "flcropnut|FAILED|PASSED"`
Expected: both new tests fail. Reasoning: at this point Task 1 has lifted the validator stub-error, so the calls go through, but Task 3 (the wiring) has not landed yet — `flCropNut` is never assigned. The first test pre-sets `flCropNut=.false.` and asserts `@assertTrue(flCropNut)` → fails. The second pre-sets `flCropNut=.true.` and asserts `@assertFalse(flCropNut)` → fails.

Confirm: 2 new failing tests, all previously-passing tests still passing.

- [ ] **Step 3: Commit the failing tests**

```bash
git add tests/unit/crop/test_cropwofost_init.pf
git commit -m "$(cat <<'EOF'
test(cropwofost): SS-N3 Task 2 — failing pFUnit tests for flCropNut wiring

Two test routines covering the per-rotation wiring contract:
  - flcropnut=true propagates to flCropNut=.true.
  - flcropnut=false resets flCropNut=.false. (per-rotation toggle)

Both fail at this commit because cropwofost_init_from_config does not
yet assign flCropNut. Task 3 lands the wiring and turns these green.

Spec: docs/superpowers/specs/2026-05-08-nutrients-N3-runtime-activation-design.md
EOF
)"
```

---

### Task 3: Wire flCropNut from typed config

The actual wiring. Two sub-edits in the same file, same routine: add `flCropNut` to the `use variables, only:` list, then assign it from `cfg%nutrient%flcropnut` and use the local copy to gate the existing `apply_cropwofost_nutrient` call.

**Files:**
- Modify: `src/crop/cropwofost_init.f90:42-62` (use list) and `src/crop/cropwofost_init.f90` near line 464-468 (existing N1 block)

- [ ] **Step 1: Add flCropNut to the use list**

In `src/crop/cropwofost_init.f90`, locate the `use variables, only:` block inside `subroutine cropwofost_init_from_config` (starts around line 41). The block currently ends with `dvs, tsum, daycrop, nofd`. Append `flCropNut`:

```fortran
         dvs, tsum, daycrop, nofd, flCropNut
```

- [ ] **Step 2: Replace the N1 nutrient call block**

Find the existing block (around line 463-468):

```fortran
      ! [nutrients] N1: apply nutrient block to legacy globals when
      ! flcropnut=true on this rotation. Validator and the runtime
      ! gate at tillage.f90:73 still control whether nutrient physics
      ! actually runs; this just makes the typed-config inputs
      ! available if/when the gate is lifted (N3).
      if (cfg%nutrient%flcropnut) call apply_cropwofost_nutrient(cfg%nutrient)
```

Replace with:

```fortran
      ! [nutrients] N3: drive the legacy global flCropNut from the
      ! per-rotation typed config. cropwofost_init_from_config runs at
      ! every rotation start, so a sequence of rotations with mixed
      ! flcropnut values toggles the gate correctly.
      flCropNut = cfg%nutrient%flcropnut
      if (flCropNut) call apply_cropwofost_nutrient(cfg%nutrient)
```

- [ ] **Step 3: Run the two new pFUnit tests — they must pass**

Run: `pixi run -e test test-pfunit 2>&1 | grep -E "flcropnut|FAILED|PASSED"`
Expected: both `test_flcropnut_true_propagates_to_global` and `test_flcropnut_false_resets_global` pass. No other tests regress.

- [ ] **Step 4: Run check-full to confirm no regression on existing cases**

Run: `pixi run check-full`
Expected: `Results: 5 passed, 0 failed`. None of the 5 cases sets `flcropnut=true`, so `flCropNut` stays `.false.` and observed behaviour is unchanged.

- [ ] **Step 5: Commit**

```bash
git add src/crop/cropwofost_init.f90
git commit -m "$(cat <<'EOF'
feat(cropwofost): SS-N3 Task 3 — wire flCropNut from typed config (ADR 0028)

Per-rotation assignment of the legacy global flCropNut from
cfg%nutrient%flcropnut, inside cropwofost_init_from_config. Mirrors
the per-rotation crop config cache pattern (ADR 0016) — every rotation
init runs before the rotation's first daily step, so the global
toggles correctly across mixed-flcropnut rotation sequences.

The runtime stub-error at tillage.f90:73 still holds the gate closed;
that comes down in Task 4.

Spec: docs/superpowers/specs/2026-05-08-nutrients-N3-runtime-activation-design.md
EOF
)"
```

---

### Task 4: Lift runtime stub-error in tillage.f90

The actual gate. After this deletion, `DoTillage` runs through normally when `flCropNut=.true.`. No test change required — Task 3's pFUnit tests already lock the wiring; this task removes a dead-only-because-the-flag-was-always-false branch.

**Files:**
- Modify: `src/crop/tillage.f90:73`

- [ ] **Step 1: Delete the runtime stub-error**

In `src/crop/tillage.f90`, delete line 73:

```fortran
         if (flCropNut)        call fatalerr_collected ('DoTillage', 'flCropNut = 1 not (yet) allowed')
```

- [ ] **Step 2: Run check-full**

Run: `pixi run check-full`
Expected: `Results: 5 passed, 0 failed`. All 5 cases have `flCropNut=.false.` (no rotation sets `flcropnut=true`), so removing this guard is invisible to them.

- [ ] **Step 3: Run pFUnit suite**

Run: `pixi run -e test test-pfunit 2>&1 | tail -5`
Expected: zero failures. Test count unchanged from Task 3.

- [ ] **Step 4: Commit**

```bash
git add src/crop/tillage.f90
git commit -m "$(cat <<'EOF'
refactor(tillage): SS-N3 Task 4 — lift flCropNut runtime stub-error (ADR 0028)

Removes the DoTillage stub-error that fatal-errored on flCropNut=.true..
With Task 3 wiring the global from typed config, the path is now safe
to traverse — closes the [nutrients] umbrella's runtime gate.

Spec: docs/superpowers/specs/2026-05-08-nutrients-N3-runtime-activation-design.md
EOF
)"
```

---

### Task 5: Stale comment cleanup

Five comment blocks across the codebase still claim `flCropNut=1 is stub-errored upstream` or equivalent. After Task 4 those statements are no longer true. Update each.

**Files:**
- Modify: `src/core/swap.f90:191-194`
- Modify: `src/crop/management_soil.f90:87-93`
- Modify: `src/crop/cropwofost_init.f90:21-22` (docstring) and `:413-415`
- Modify: `src/crop/cropgrowth.f90:1135-1143`

- [ ] **Step 1: Update src/core/swap.f90 comment**

Locate lines 191-194:

```fortran
!  Soil Management init: SoilManagement(1) was the legacy reader entry
!  point and is now a no-op (SS-C step 3). flCropNut=1 is stub-errored
!  upstream in DoTillage; the SoilManagement(2..7) call sites below remain
!  for the eventual TOML-port reactivation (ADR 0021).
```

Replace with:

```fortran
!  Soil Management init: SoilManagement(1) was the legacy reader entry
!  point and is now a no-op (SS-C step 3). flCropNut is now driven by
!  the per-rotation typed config (ADR 0028); the SoilManagement(2..7)
!  call sites below run when a rotation has flcropnut=true.
```

- [ ] **Step 2: Update src/crop/management_soil.f90 comment**

Locate the comment in `case (1)` of `SoilManagement` (around line 87-93):

```fortran
      ! Legacy nutrient soil-management init (file-open + state
      ! initialization) deleted as part of legacy readers physical
      ! deletion. flCropNut=1 is stub-errored upstream so all
      ! SoilManagement(*) call sites are unreachable in the modern
      ! flow. A future TOML port (companion to ADR 0021) will
      ! reintroduce this with typed-config-driven init.
```

Replace with:

```fortran
      ! Legacy nutrient soil-management init (file-open + state
      ! initialization) deleted as part of legacy readers physical
      ! deletion. Initial soil-nutrient state now flows from the
      ! [nutrients] typed config via apply_nutrients (ADR 0026 N2a).
```

- [ ] **Step 3: Update src/crop/cropwofost_init.f90 docstring (line 21)**

Locate the docstring comment near line 21-22:

```fortran
!!     FraHarLosOrm_lv/st/so  (harvest/death losses)
!!       — set by the N-P-K block (readwofost lines 1086-1088), gated on
!!         flCropNut; stub-errored for TOML path so no action needed.
```

Replace with:

```fortran
!!     FraHarLosOrm_lv/st/so  (harvest/death losses)
!!       — set by apply_cropwofost_nutrient when flcropnut=true on the
!!         current rotation (ADR 0025 N1). For flcropnut=false rotations
!!         these stay at their wofost() local defaults.
```

- [ ] **Step 4: Update src/crop/cropwofost_init.f90 inline comment (line ~413)**

Locate the comment around line 413-415:

```fortran
      ! FraDeceasedLvToSoil — local SAVE in wofost(), returned via intent(out)
      ! so the dispatch block can assign it.  (FraHarLosOrm_* are set by the
      ! N-P-K block which is gated on flCropNut; stub-guarded above.)
```

Replace with:

```fortran
      ! FraDeceasedLvToSoil — local SAVE in wofost(), returned via intent(out)
      ! so the dispatch block can assign it.  (FraHarLosOrm_* are set by
      ! apply_cropwofost_nutrient below when flcropnut=true on this rotation.)
```

- [ ] **Step 5: Update src/crop/cropgrowth.f90 comment block**

Locate the block around line 1135-1143:

```fortran
!        Legacy nutrient parameters (LRNR, LSNR, NLAI, NLUE, NMAXSO,
!        NPART, NFIXF, NSLA, RNFLV/RT/ST, TCNT, DVSNLT, DVSNT, RDRNS,
!        FNTRT, FRNX, NMXLV, FraHarLosOrm_lv/st/so) used to be read here
!        from <cropfil>.crp via TTutil rdinit/rdsdou. Read block deleted
!        as part of legacy readers physical deletion. flCropNut=1 is
!        stub-errored upstream in DoTillage (tillage.f90), so this block
!        is unreachable in the modern flow. A future TOML port of the
!        nutrient sub-block (companion to ADR 0021) will reintroduce
!        these reads from typed config.
```

Replace with:

```fortran
!        Legacy nutrient parameters (LRNR, LSNR, NLAI, NLUE, NMAXSO,
!        NPART, NFIXF, NSLA, RNFLV/RT/ST, TCNT, DVSNLT, DVSNT, RDRNS,
!        FNTRT, FRNX, NMXLV, FraHarLosOrm_lv/st/so) used to be read here
!        from <cropfil>.crp via TTutil rdinit/rdsdou. Read block deleted
!        as part of legacy readers physical deletion. These globals are
!        now populated by apply_cropwofost_nutrient when flcropnut=true
!        on the active rotation (ADR 0025 N1, ADR 0028 N3).
```

- [ ] **Step 6: Build and run the full verification**

Run: `pixi run -e test test-pfunit 2>&1 | tail -5 && pixi run check-full`
Expected: zero pFUnit failures, `Results: 5 passed, 0 failed` from check-full.

- [ ] **Step 7: Commit**

```bash
git add src/core/swap.f90 src/crop/management_soil.f90 src/crop/cropwofost_init.f90 src/crop/cropgrowth.f90
git commit -m "$(cat <<'EOF'
docs(crop): SS-N3 Task 5 — refresh stale flCropNut stub-error comments

Five comment blocks across swap.f90, management_soil.f90,
cropwofost_init.f90 (×2), and cropgrowth.f90 still claimed flCropNut=1
was stub-errored upstream. After Tasks 1+3+4 that is no longer true —
the global is driven per-rotation from typed config (ADR 0028) and the
gate is open. Comments rewritten to reflect current state.

Spec: docs/superpowers/specs/2026-05-08-nutrients-N3-runtime-activation-design.md
EOF
)"
```

---

### Task 6: ADR 0028 + final umbrella verification

Author the ADR, update the index, and run the full umbrella verification one last time.

**Files:**
- Create: `docs/adr/0028-nutrients-N3-runtime-activation.md`
- Modify: `docs/adr/index.md`

- [ ] **Step 1: Read the existing ADR shape for consistency**

Run: `cat docs/adr/0027-nutrients-N2b-timed-amendments.md | head -40`
Expected: get the section headings, status line format, and reference style used by ADRs 0025/0026/0027.

- [ ] **Step 2: Author docs/adr/0028-nutrients-N3-runtime-activation.md**

Match the section structure of ADR 0027 (the immediate predecessor). Required sections:

```markdown
# ADR 0028 — [nutrients] N3 — Runtime Activation

**Status:** accepted
**Date:** 2026-05-08
**Sub-arc of:** [nutrients] umbrella (N1 → N2a → N2b → **N3**)

## Context

After ADRs 0025/0026/0027, all nutrient typed-config flows to the
legacy globals via three adapters: apply_cropwofost_nutrient (per
rotation, N1), apply_nutrients (top-level [nutrients], N2a), and
apply_nutrients_events (CSV companion, N2b). The legacy global
flCropNut, which gates nutrient physics inside the daily loop and
inside DoTillage, was never assigned anywhere on the modern path —
it defaulted to .false. and stayed there. Two stub-errors held the
gate closed: a validator at cropwofost_init.f90:85 and a runtime
guard at tillage.f90:73.

## Decision

Drive flCropNut per-rotation from cfg%nutrient%flcropnut, inside
cropwofost_init_from_config. Remove both stub-errors. Defer the
regression fixture authorship — no legacy nutrient-enabled .crp
exists to mirror, and the parity baseline question is a separate
workstream.

Per-rotation assignment is correct because cropwofost_init_from_config
runs at every rotation start, before that rotation's first daily step.
A sequence of rotations with mixed flcropnut values therefore toggles
the global correctly. This mirrors the per-rotation crop config cache
pattern (ADR 0016).

An OR-reduction across all rotations was rejected — it would enable
nutrient gates for rotations whose typed config explicitly disabled
them.

## Consequences

- The [nutrients] umbrella is complete on the TOML path. A rotation
  with cropwofost.nutrient.flcropnut = true now runs the WOFOST
  nutrient subsystem end-to-end.
- No new dependencies introduced — cropwofost_init.f90 is
  TTutil-clean (modernized in N1) and the change is one assignment
  plus the existing N1 adapter call.
- Correctness of nutrient outputs is not validated by this arc. No
  legacy parity fixture exists; first regression case is a separate
  workstream.
- Stale "flCropNut=1 is stub-errored upstream" comments across
  swap.f90, management_soil.f90, cropwofost_init.f90 and
  cropgrowth.f90 have been refreshed.

## References

- ADR 0025 — N1 crop-side adapter
- ADR 0026 — N2a soil-side initial state
- ADR 0027 — N2b timed amendments
- ADR 0023 / 0024 — TTutil retirement context
- ADR 0016 — per-rotation crop config cache (the precedent for the
  per-rotation assignment pattern)
- Spec: docs/superpowers/specs/2026-05-08-nutrients-N3-runtime-activation-design.md
```

- [ ] **Step 3: Add the entry to docs/adr/index.md**

Append (or insert in numerical order):

```markdown
- [0028 — [nutrients] N3 — Runtime Activation](0028-nutrients-N3-runtime-activation.md)
```

- [ ] **Step 4: Final full verification**

Run all in sequence:

```bash
pixi run -e test test-pfunit 2>&1 | tail -10
pixi run check-full
```

Expected:
- pFUnit: zero failures, test count = previous + 2 (the two new flcropnut tests).
- check-full: `Results: 5 passed, 0 failed`.

- [ ] **Step 5: Hand-verification (one-off, not committed)**

This is an out-of-band sanity check. Pick one TOML case (suggest `tests/swap-cases/toml/2.grassgrowth/grassd.crp.toml` because grass already has nutrient-friendly fields). Edit a copy or in-place, set `flcropnut = true` under `[wofost.nutrient]`, populate `nmxlv = [60.0]`, run:

```bash
cd tests/swap-cases/toml/2.grassgrowth && bash ../../run_case.sh
```

Confirm:
- Build does **not** abort with `'flcropnut=.true. not supported on TOML path'`.
- Build does **not** abort with `'flCropNut = 1 not (yet) allowed'`.
- Simulation completes (output may differ from baseline; that is expected).

**Revert any local file edit before commit.** Run `git status` and `git diff` to confirm clean.

- [ ] **Step 6: Commit**

```bash
git add docs/adr/0028-nutrients-N3-runtime-activation.md docs/adr/index.md
git commit -m "$(cat <<'EOF'
docs(adr): ADR 0028 — [nutrients] N3 — runtime activation (umbrella close-out)

Closes the [nutrients] umbrella (N1 → N2a → N2b → N3). Records the
decision to drive flCropNut per-rotation from typed config and to
defer the regression fixture authorship to a separate arc.

Spec: docs/superpowers/specs/2026-05-08-nutrients-N3-runtime-activation-design.md
EOF
)"
```

---

## Self-Review Notes

- **Spec coverage:** Section 1 (goal, out-of-scope) → covered by Tasks 1-6 collectively. Section 2 (stub-errors) → Tasks 1, 4. Section 3 (setting the global) → Task 3. Section 4 (testing) → Tasks 2+3 (pFUnit), Task 6 (check-full + hand-verify). Section 5 (ADR 0028) → Task 6.
- **Per-rotation rationale:** baked into Task 3 commit message and ADR 0028 text.
- **Non-regression TTutil-free claim:** no new module imports anywhere except `flCropNut` from `variables` (already in scope of `cropwofost_init.f90` via existing `use variables, only:`). Zero `dtutil` / `ttutil` references introduced.
- **Test fixture self-contained:** the two new pFUnit tests duplicate the fixture-setup boilerplate from `test_cropwofost_init_basic_globals`. Yes, that's a bit of repetition; pulling out a helper subroutine would touch more files and obscure the test logic. Acceptable for a 6-task arc.
