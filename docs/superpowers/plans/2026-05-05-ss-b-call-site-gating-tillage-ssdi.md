# SS-B — ADR 0020 + lift tillage/SSDI gates to call site

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bring `Read_Tillage` and `SSDI_irrigation` in line with the
existing call-site-gating convention used by every other optional
subsystem in the project. Codify the convention in ADR 0020. Remove
the SS-10.5 stub-error validators for `swtill=1` / `swssdi=1` —
those features are legitimate, not deprecated. The TOML port of the
tillage / SSDI parameter blocks themselves is deferred to a future
arc (ADR 0021 candidate); this sub-spec preserves the legacy TTutil
read path inside the subsystem bodies as the transition state.

**Sub-spec parent:**
`docs/superpowers/specs/2026-05-05-legacy-readers-physical-deletion-design.md`
(SS-B).

**Architecture.** Existing `swap.f90` pattern — `if (flX) call X(...)`
— is applied to `DoTillage(*)` and `SSDI_irrigation(*)`. Adapter
populates `flTillage` and `flSSDI` from `config%soil%swtill` and
`config%irrigation%swssdi`. Subsystem bodies drop the now-redundant
`if (swtill /= 1) return` self-check; entry implies the flag is
true. Validators stop rejecting `=1`.

**Tech stack:** Fortran 2008, pFUnit, pixi, meson. No new TTutil use;
no new schema fields beyond what SS-10.5 already added.

**Effort:** S-M (~6-7 commits).

---

## File structure

| File | Change |
|---|---|
| `docs/adr/0020-call-site-gating-convention.md` | **Create:** ADR documenting the convention. |
| `src/core/variables.f90` | **Modify:** add `logical :: flTillage = .false.` and `logical :: flSSDI = .false.`. |
| `src/io/toml/config_to_variables.f90` | **Modify:** set `flTillage = (config%soil%swtill == 1)`, `flSSDI = (config%irrigation%swssdi == 1)`. Restore `swpfile = 'swap.swp'` assignment (load-bearing if either flag is true). |
| `src/config/soil_config.f90` | **Modify:** delete the `swtill==1` stub-error block (keep enum check). |
| `src/config/irrigation_config.f90` | **Modify:** delete the `swssdi==1` stub-error block (keep enum check). |
| `src/crop/tillage.f90` | **Modify:** drop `if (swtill /= 1) return` self-check from `Read_Tillage`. Update header comment to point at ADR 0020. |
| `src/crop/irrigation.f90` (`SSDI_irrigation` `case (1)`) | **Modify:** drop the `if (swssdi == 0) return` self-check (the call-site gate replaces it). Move the swssdi=0 default-assignment up into the adapter or delete (the globals already initialise to those values). |
| `src/core/swap.f90` | **Modify:** wrap 3 `call DoTillage(*)` sites with `if (flTillage)`; wrap 2 `call SSDI_irrigation(*)` sites with `if (flSSDI)`. |
| `src/core/timecontrol.f90` | **Modify:** wrap the 1 `call SSDI_irrigation(9)` site with `if (flSSDI)`. |
| `tests/unit/config/test_soil_config.pf` | **Modify:** invert `test_soil_swtill_one_stub_errors` — assert `swtill=1` is now ACCEPTED. Or delete and add an enum-only test. |
| `tests/unit/config/test_irrigation_config.pf` | **Modify:** invert `test_irrigation_swssdi_one_stub_errors` similarly. |

---

## Conventions

1. **ADR first.** Write `docs/adr/0020-call-site-gating-convention.md`
   before touching any source. Reference it from updated subsystem
   header comments.
2. **One commit per change-type:**
   - Commit 1: ADR 0020
   - Commit 2: `flTillage` / `flSSDI` flags + adapter wiring + restore swpfile
   - Commit 3: validator un-stub + test inversion
   - Commit 4: subsystem self-check removal (tillage + SSDI)
   - Commit 5: call-site gate insertion in swap.f90 + timecontrol.f90
   - Commit 6: verify gates green; final merge-prep
3. **Verification at every commit:** `pixi run -e test build-linux` →
   `pixi run -e test test-pfunit` → `pixi run -e test check-full`.
   All three must be green before the next commit.
4. **No behaviour change for any current test case.** All cases have
   `swtill=0` and `swssdi=0`; with the new flags also defaulting to
   `.false.`, the call-site gates short-circuit and the subsystems
   never run. This must be visible in the gate output.

---

## Task 1: Write ADR 0020

**Commit:** 1
**Files:** `docs/adr/0020-call-site-gating-convention.md` (new),
`docs/adr/index.md` (append entry)

- [ ] **Step 1:** Create the ADR with the structure used by ADRs
  0017–0019. Required sections:
  - Context: existing call-site gating examples in `swap.f90` (~10
    sites with `if (flX) call X(*)`); the two outliers (`DoTillage`,
    `SSDI_irrigation`) currently using internal self-checks.
  - Decision: codify call-site gating as the convention; bring
    tillage and SSDI in line.
  - Consequences: each optional subsystem owns a `flX` flag in
    `variables.f90`; the adapter sets it from the typed config; the
    subsystem body assumes it's only invoked when needed.
  - Note that this ADR does NOT decide whether tillage and SSDI
    keep their TTutil-based reading. ADR 0021 (candidate, future)
    decides that.

- [ ] **Step 2:** Append to `docs/adr/index.md`:
  `- [ADR 0020 — Call-site gating convention](0020-call-site-gating-convention.html) — All optional subsystems are gated at the call site by an `flX` flag; subsystem implementations assume the flag is true. Brings tillage and SSDI in line with the existing pattern (flCropNut, flMacroPore, flSurfaceWater, …).`

- [ ] **Step 3:** Commit.

---

## Task 2: Add `flTillage` / `flSSDI` flags + adapter wiring

**Commit:** 2
**Files:** `src/core/variables.f90`, `src/io/toml/config_to_variables.f90`

- [ ] **Step 1:** In `src/core/variables.f90`, locate the section
  with other `flX` logical declarations (search for
  `logical\s*::\s*flCropNut` or similar). Add two new entries:

```fortran
logical :: flTillage = .false.   !! ADR 0020 call-site gate for DoTillage
logical :: flSSDI    = .false.   !! ADR 0020 call-site gate for SSDI_irrigation
```

- [ ] **Step 2:** In `src/io/toml/config_to_variables.f90`, immediately
  after the existing `till_swtill = config%soil%swtill` /
  `swssdi_irr = config%irrigation%swssdi` lines (added in SS-10.5),
  add:

```fortran
flTillage = (config%soil%swtill == 1)
flSSDI    = (config%irrigation%swssdi == 1)
```

- [ ] **Step 3:** Restore the `swpfile` assignment that SS-10.5
  dropped. In the same `config_to_variables` body, before the
  `logf`-opening block, add:

```fortran
! ADR 0020: Read_Tillage / SSDI_irrigation(1) read tillage / SSDI
! parameters from staged swap.swp via TTutil when their flag is set.
! Until ADR 0021 ports those blocks to TOML schema, swpfile must
! point at the staged template.
swpfile = 'swap.swp'
```

- [ ] **Step 4:** Build → pFUnit → check-full. All green. Commit.

---

## Task 3: Un-stub validators + invert tests

**Commit:** 3
**Files:** `src/config/soil_config.f90`,
`src/config/irrigation_config.f90`,
`tests/unit/config/test_soil_config.pf`,
`tests/unit/config/test_irrigation_config.pf`

- [ ] **Step 1:** In `src/config/soil_config.f90`, delete the
  Phase-4f-extend SS-10.5 stub-error block:

```fortran
! DELETE this block:
if (self%swtill == 1) then
   call errors%append(ERR_VALIDATION_CROSS_FIELD, &
      'soil.swtill=1 (tillage events) not yet supported in the ' // ...
end if
```

The `check_int_enum(self%swtill, [0, 1], "soil.swtill", errors)` call
stays — that's the legitimate value-range check.

- [ ] **Step 2:** Same in `src/config/irrigation_config.f90` for
  `self%swssdi == 1`. Keep the enum check.

- [ ] **Step 3:** In `tests/unit/config/test_soil_config.pf`, replace
  `test_soil_swtill_one_stub_errors` with:

```fortran
@test
subroutine test_soil_swtill_one_accepted()
   ! ADR 0020: swtill=1 is a legitimate value (tillage activated).
   ! Validation accepts it; the call-site gate (flTillage) controls
   ! whether the tillage subsystem actually runs at runtime.
   use funit
   use error_mod, only: error_collection_t
   use soil_config_mod, only: soil_config_t
   type(soil_config_t)      :: s
   type(error_collection_t) :: errors

   s%swtill = 1
   call s%validate(errors)
   @assertFalse(errors%has_errors())
end subroutine
```

- [ ] **Step 4:** Same shape for `test_irrigation_swssdi_one_accepted`.

- [ ] **Step 5:** Build → pFUnit → check-full. All green. Commit.

---

## Task 4: Drop subsystem self-checks

**Commit:** 4
**Files:** `src/crop/tillage.f90`, `src/crop/irrigation.f90`

- [ ] **Step 1:** In `src/crop/tillage.f90` `subroutine Read_Tillage`,
  delete the SS-10.5 early-return:

```fortran
! Delete:
if (swtill /= 1) return
```

Update the header comment to:

```fortran
! Phase 4f-extend SS-B (ADR 0020): entry to this routine implies the
! call-site gate flTillage was true → swtill=1 was set in the TOML
! config and copied to the global. Reads tillage parameters from
! staged swap.swp via TTutil; future ADR 0021 will replace this with
! a TOML schema port of the tillage block.
```

- [ ] **Step 2:** In `src/crop/irrigation.f90` `subroutine
  SSDI_irrigation` `case (1)`, delete the SS-10.5 early-return:

```fortran
! Delete:
if (swssdi == 0) then
   nod_ssdi = 0
   qssdi    = 0.0
   nirri    = 1
   dt_SSDI_event = 1.0d0
   return
end if
```

The defaults (`nod_ssdi = 0`, `qssdi = 0.0`, etc.) are now set by
the adapter or by Fortran's default initialisation; the
`SSDI_irrigation(1)` call no longer fires when `flSSDI = .false.`.

Update the case comment:

```fortran
case (1)
   ! ADR 0020 SS-B: entry implies flSSDI was true. Read SSDI
   ! parameters from staged swap.swp via TTutil; future ADR 0021
   ! will TOML-port the SSDI block.
   swp = getun2(10, 90, 2)
   call rdinit(swp, logf, swpfile)
      ! [no rdsinr('swssdi') needed — already known from TOML]
      if (rdinqr('ssdi_file')) call rdscha ('ssdi_file', ssdi_file)
   close(swp)
   call read_ssdi_input()
   ! ... rest unchanged ...
```

- [ ] **Step 3:** Build. The build SHOULD succeed because no other
  call-site invokes these routines yet (they're still
  unconditionally invoked from `swap.f90` / `timecontrol.f90`). With
  the self-check gone, the subsystem will now FIRE on every test run.
  pFUnit and check-full might break here.

- [ ] **Step 4:** Verify failure mode. Run `check-full`. If it fails
  with an SSDI/tillage runtime error, that's expected (and the
  reason Task 5 must follow before this commit lands). Do NOT commit
  yet; proceed directly to Task 5.

---

## Task 5: Add call-site gates in `swap.f90` and `timecontrol.f90`

**Commit:** 4 (combined with Task 4 — they must land together to
avoid breaking the build / regression)
**Files:** `src/core/swap.f90`, `src/core/timecontrol.f90`

- [ ] **Step 1:** In `src/core/swap.f90`, find each
  `call DoTillage(*)` line (3 sites: lines 165, 250, 355) and wrap
  with `if (flTillage)`:

```fortran
! Before:
call DoTillage(1)
! After:
if (flTillage) call DoTillage(1)
```

Same for `call DoTillage(2)` and `call DoTillage(3)`.

- [ ] **Step 2:** In `src/core/swap.f90`, find each
  `call SSDI_irrigation(*)` line (2 sites: lines 166, 345) and wrap
  with `if (flSSDI)`:

```fortran
if (flSSDI) call SSDI_irrigation(1)
if (flSSDI) call SSDI_irrigation(2)
```

- [ ] **Step 3:** In `src/core/timecontrol.f90`, find the 1
  `call SSDI_irrigation(9)` site (line 479) and wrap:

```fortran
if (flSSDI) call SSDI_irrigation(9)
```

- [ ] **Step 4:** Build → pFUnit → check-full. All gates green.
  Commit Tasks 4 and 5 together with message:

```
refactor(crop): SS-B — call-site gating for tillage and SSDI (ADR 0020)

Bring DoTillage and SSDI_irrigation in line with the existing
call-site gating convention used by ~10 other optional subsystems
in swap.f90. Per ADR 0020, optional subsystems are invoked as
`if (flX) call X(...)`; the subsystem assumes the flag is true and
no longer self-checks.

Changes:
- swap.f90: 3x `if (flTillage) call DoTillage(*)`, 2x `if (flSSDI)
  call SSDI_irrigation(*)`.
- timecontrol.f90: 1x `if (flSSDI) call SSDI_irrigation(9)`.
- tillage.f90:Read_Tillage: drop self-check; entry implies flTillage.
- irrigation.f90:SSDI_irrigation(1): drop self-check; entry implies
  flSSDI.

Behaviour preserved for all 5 regression cases (all default to
swtill=0, swssdi=0, so flTillage/flSSDI = .false. and the call-site
gates short-circuit).
```

---

## Task 6: Final verification + audit pass

**Commit:** 5
**Files:** none (verification only)

- [ ] **Step 1:** Confirm no production code path reaches a
  TTutil rdinit when flTillage = flSSDI = .false.:

```bash
# Search for any unconditional rdinit/RDinit in production source
grep -rn "call rdinit\|call RDinit" src/ --include='*.f90' \
  | grep -v "^src/io/readswap.f90:"
```

Expected:
- `src/crop/tillage.f90:NN call RDinit(IunIn, 0, swpfile)` —
  unreachable when `flTillage = .false.` (call-site gate)
- `src/crop/irrigation.f90:NN call rdinit(swp, logf, swpfile)` and
  `src/crop/irrigation.f90:NN call rdinit(swp, logf, ssdi_file)` —
  both unreachable when `flSSDI = .false.`
- `src/crop/irrigation.f90:79 call rdinit(irr,logf,filnam)` —
  legacy `case (1)` of `irrigation`, unreachable per SS-6 audit;
  scheduled for deletion in SS-C.

(`readswap.f90` calls — also dead in production, deleted in SS-C.)

- [ ] **Step 2:** Confirm test gates: `pixi run -e test test-pfunit`
  → `Ok: 1, Fail: 0`. `pixi run -e test check-full` → 5 passed.

- [ ] **Step 3:** No commit needed unless verification surfaced a
  fix. SS-B is closed when these checks pass.

---

## Definition of done (SS-B)

- ADR 0020 written and indexed.
- `flTillage` and `flSSDI` declared in `variables.f90`; adapter
  populates them from `config%soil%swtill` and
  `config%irrigation%swssdi`.
- `swpfile = 'swap.swp'` restored in adapter (load-bearing for
  Read_Tillage / SSDI when flags are true).
- `soil_config_validate` and `irrigation_config_validate` no longer
  reject `=1` for `swtill` / `swssdi` (enum check stays).
- `test_soil_swtill_one_accepted` and
  `test_irrigation_swssdi_one_accepted` (or equivalent) pass.
- `Read_Tillage` and `SSDI_irrigation(1)` no longer carry self-check
  early-returns; their headers reference ADR 0020 and ADR 0021
  (candidate).
- 6 call sites in `swap.f90` and `timecontrol.f90` wrapped with
  `if (flTillage)` / `if (flSSDI)`.
- `pixi run -e test test-pfunit` → `Ok: 1, Fail: 0`.
- `pixi run -e test check-full` → 5 passed.
- 4–5 commits, all individually bisectable.
