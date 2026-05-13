# GR-UTILS Implementation Plan — Utils Globals Retirement + cofgen Modernization

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Retire `use variables` and the `bind_*_target` pointer-binding pattern from `src/utils/soilhydraulicsutils.f90`, `src/utils/surfacewaterutils.f90`, and `src/soil/WC_K_models_04_11.f90` by (1) introducing a typed `vanGenuchten_params_t` record that replaces the 21-row `cofgen` magic-row matrix, and (2) migrating 17 additional bare-global symbols into the appropriate state subrecords. Utility functions get explicit, self-documenting signatures.

**Architecture:** Define `vanGenuchten_params_t` in a new small module. Replace `state%soilwater%cofgen(:,:)` with `state%soilwater%vg_params(:)` (array of derived types — same bit-level data, named fields). Update 5 functions in `soilhydraulicsutils` + `functionvalue_04_11` in `WC_K_models_04_11` to take per-node `vg` slices and other previously-implicit args (model index, dt, flags) as explicit dummy args. `surfacewaterutils` migrates by reading from `state%surfacewater%X` directly. All `bind_*_target` procedures and module-level state pointers retire.

**Tech Stack:** Fortran 2018, gfortran, Meson, pFUnit 4.15, pixi.

**Spec:** `docs/superpowers/specs/2026-05-13-globals-utils-design.md`

**Builds on:** SS-DRV Phase 1, SS-TCM, SS-BMI2.

**Verification per task** (per memory `feedback_state_schema_clean_rebuild.md` + `feedback_per_task_regression_gate.md`):
- **Clean rebuild** required after any state schema change: `rm -rf builddir && pixi run build-linux`
- `pixi run test-pfunit` — 741/741 expected (or current total)
- `pixi run -e test python tests/regression/test_output_regression.py hupselbrook surfacewater salinitystress grassgrowth` — 4/4 byte-for-byte

---

## Task 1: Pre-flight baseline

**Files:** No code changes; capture starting state.

- [ ] **Step 1: Build clean.** Run: `pixi run build-linux` — must succeed.
- [ ] **Step 2: pFUnit baseline.** Run: `pixi run test-pfunit` — record total pass count.
- [ ] **Step 3: check-full baseline.** Run: `pixi run check-full` — 5/5 pass.
- [ ] **Step 4: BMI + cffi-demo suites baseline.** Run: `pixi run -e test test-bmi && pixi run -e test test-cffi-demo` — both pass.
- [ ] **Step 5: Commit baseline marker.**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-utils): pre-flight baseline — globals-retirement Arc 1 begins

check-full: 5/5. pFUnit: passing. BMI + cffi-demo suites passing.
Baseline locked before cofgen modernization + utils globals retirement.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Define `vanGenuchten_params_t` + add `vg_params(:)` field to state

**Files:**
- Create: `src/state/hydraulic_params_mod.f90`
- Modify: `src/state/soilwater_state.f90` (add `vg_params` field + use of new module + allocation in init)
- Modify: `meson.build` (add new file to legacy sources list, before `soilwater_state.f90`)
- Modify: `tests/unit/meson.build` (add new file to `pfunit_extra_sources`, before `soilwater_state.f90`)

This task is additive — `vg_params(:)` is allocated but unused. `cofgen(:,:)` keeps working. No behavior change.

- [ ] **Step 1: Create the new module.**

Create `src/state/hydraulic_params_mod.f90`:

```fortran
!> @file hydraulic_params_mod.f90
!! SS-GR-UTILS: typed van-Genuchten / MvG / PDI hydraulic parameter
!! record. Replaces the legacy `cofgen(j, node)` magic-row matrix.
!! One instance per soil node; aggregated as
!! `state%soilwater%vg_params(:)`.
!!
!! Row-name correspondence to legacy cofgen indices
!! (from src/soil/soilhydraulics.f90:830-880 + utils/WC_K_models):
!!   cofgen(1, n)  -> vg_params(n)%thetar       (residual water content)
!!   cofgen(2, n)  -> vg_params(n)%thetas       (saturated water content)
!!   cofgen(3, n)  -> vg_params(n)%ksat
!!   cofgen(4, n)  -> vg_params(n)%alpha
!!   cofgen(5, n)  -> vg_params(n)%lpar
!!   cofgen(6, n)  -> vg_params(n)%npar
!!   cofgen(7, n)  -> vg_params(n)%mpar
!!   cofgen(8, n)  -> vg_params(n)%alphaw_sentinel  (=-9999.9 sentinel)
!!   cofgen(9, n)  -> vg_params(n)%h_enpr
!!   cofgen(10, n) -> vg_params(n)%ksatexm
!!   cofgen(11, n) -> vg_params(n)%relsatthr
!!   cofgen(12, n) -> vg_params(n)%ksatthr
!!   cofgen(13, n) -> vg_params(n)%alpha_2
!!   cofgen(14, n) -> vg_params(n)%npar_2
!!   cofgen(15, n) -> vg_params(n)%mpar_2
!!   cofgen(16, n) -> vg_params(n)%omega_1
!!   cofgen(17, n) -> vg_params(n)%omega_2
!!   cofgen(18, n) -> vg_params(n)%h0
!!   cofgen(19, n) -> vg_params(n)%ha
!!   cofgen(20, n) -> vg_params(n)%apar
!!   cofgen(21, n) -> vg_params(n)%omega_k
module hydraulic_params_mod
   use iso_fortran_env, only: real64
   implicit none
   private
   public :: vanGenuchten_params_t

   type :: vanGenuchten_params_t
      ! Universal MvG (rows 1-7)
      real(real64) :: thetar           = 0.0_real64
      real(real64) :: thetas           = 0.0_real64
      real(real64) :: ksat             = 0.0_real64
      real(real64) :: alpha            = 0.0_real64
      real(real64) :: lpar             = 0.0_real64
      real(real64) :: npar             = 0.0_real64
      real(real64) :: mpar             = 0.0_real64
      ! Sentinel + extreme-K parameters (rows 8-12)
      real(real64) :: alphaw_sentinel  = 0.0_real64   ! row 8: -9999.9 marker (set by SoilHydraulics(1))
      real(real64) :: h_enpr           = 0.0_real64   ! row 9
      real(real64) :: ksatexm          = 0.0_real64   ! row 10
      real(real64) :: relsatthr        = 0.0_real64   ! row 11
      real(real64) :: ksatthr          = 0.0_real64   ! row 12
      ! Bi-modal MvG (rows 13-17)
      real(real64) :: alpha_2          = 0.0_real64
      real(real64) :: npar_2           = 0.0_real64
      real(real64) :: mpar_2           = 0.0_real64
      real(real64) :: omega_1          = 0.0_real64
      real(real64) :: omega_2          = 0.0_real64
      ! PDI parameters (rows 18-21)
      real(real64) :: h0               = 0.0_real64
      real(real64) :: ha               = 0.0_real64
      real(real64) :: apar             = 0.0_real64
      real(real64) :: omega_k          = 0.0_real64
   end type vanGenuchten_params_t

end module hydraulic_params_mod
```

- [ ] **Step 2: Add `vg_params(:)` field to `soilwater_state_t`.**

Read `src/state/soilwater_state.f90:140` to confirm the current `cofgen` line. Just below (or above) it, add the new field declaration. Use Edit:

```fortran
      real(real64), allocatable :: cofgen(:,:)     !< Mualem-VG parameters (21 × numnod)
```
→
```fortran
      real(real64), allocatable :: cofgen(:,:)     !< Mualem-VG parameters (21 × numnod) — [SS-GR-UTILS] retiring in this arc
      type(vanGenuchten_params_t), allocatable :: vg_params(:)   !< [SS-GR-UTILS] typed VG parameters, one per node
```

Add `use hydraulic_params_mod, only: vanGenuchten_params_t` to the module's `use` list (top of `soilwater_state_mod`).

- [ ] **Step 3: Allocate `vg_params` in soilwater_init.**

In `src/state/soilwater_state.f90` around line 368 where `sw%cofgen` is allocated:

```fortran
      allocate(sw%cofgen(21, numnod));   sw%cofgen       = 0.0_real64
```
→
```fortran
      allocate(sw%cofgen(21, numnod));   sw%cofgen       = 0.0_real64
      allocate(sw%vg_params(numnod))    ! [SS-GR-UTILS] components default-init from type
```

The type's component-defaults zero each field; no explicit zeroing needed.

- [ ] **Step 4: Add `hydraulic_params_mod.f90` to `meson.build`.**

In `meson.build`, find the `sources = [...]` list, find where `'src/state/soilwater_state.f90'` is listed. Add immediately BEFORE it:

```
    'src/state/hydraulic_params_mod.f90',
```

- [ ] **Step 5: Add to `tests/unit/meson.build`.**

In `tests/unit/meson.build`'s `pfunit_extra_sources`, add `'../../src/state/hydraulic_params_mod.f90',` before the soilwater_state entry.

- [ ] **Step 6: Build + verify.**

Run: `rm -rf builddir && pixi run build-linux`. Build must succeed.

Run: `pixi run test-pfunit` — total pass count unchanged.

Run: `pixi run -e test python tests/regression/test_output_regression.py hupselbrook surfacewater salinitystress grassgrowth` — 4/4 byte-for-byte.

- [ ] **Step 7: Commit.**

```bash
git add src/state/hydraulic_params_mod.f90 src/state/soilwater_state.f90 \
        meson.build tests/unit/meson.build
git commit -m "$(cat <<'EOF'
schema(gr-utils): add vanGenuchten_params_t + vg_params(:) field

New small module src/state/hydraulic_params_mod.f90 defines the
typed VG/MvG/PDI parameter record (20 named fields covering all
21 legacy cofgen rows). state%soilwater%vg_params(:) added
alongside the existing cofgen(:,:) array. Allocated at
soilwater_init; populated in Task 3.

Additive change: vg_params is unused; cofgen remains the live data
source. No behavior change. Regression 4/4 byte-for-byte.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Populate `vg_params(:)` in `SoilHydraulics(1)` — transitional dual-write

**Files:**
- Modify: `src/soil/soilhydraulics.f90` (the `SoilHydraulics(task=1)` block — `cofgen` population at lines 830-880)

Add the parallel `vg_params(:)` writes alongside the existing `cofgen(:,:)` writes. Both data sources hold the same values; subsequent tasks migrate readers from `cofgen` to `vg_params`.

- [ ] **Step 1: Read the cofgen population block.**

Read `src/soil/soilhydraulics.f90:825-880`. Two branches:
- `swsophy == 1` (tabulated, lines 830-853): populates `sw%cofgen(1)`, `(2)`, `(3)` from `sptab`.
- `swsophy != 1` (analytical, lines 854-877): populates `sw%cofgen(1..10)` from `paramvg`, plus `(11)`, `(12)`, conditional `(13..17)`, `(18)`, `(18..21)`.

- [ ] **Step 2: Add parallel vg_params writes — tabulated branch (lines 845-848).**

Find the existing 4-line block:

```fortran
          sw%cofgen(1,node) = 0.0_real64                  ! thetar
          sw%cofgen(2,node) = sptab(2,node,numtab(node))  ! thetas
          sw%cofgen(3,node) = sptab(3,node,numtab(node))  ! ksat
          if (do_ln_trans) sw%cofgen(3,node) = dexp(sw%cofgen(3,node))
```

Append immediately after:

```fortran
          ! [SS-GR-UTILS] Mirror tabulated cofgen writes into typed vg_params
          sw%vg_params(node)%thetar = sw%cofgen(1,node)
          sw%vg_params(node)%thetas = sw%cofgen(2,node)
          sw%vg_params(node)%ksat   = sw%cofgen(3,node)
```

- [ ] **Step 3: Add parallel vg_params writes — analytical branch (lines 858-877).**

Find:

```fortran
          do i = 1, 10
            sw%cofgen(i,node) = paramvg(i,lay)            ! [SS-SWC S-1.3/S-2.12B]
          end do
          ! Assign dummy value to alphaw
          sw%cofgen(8,node) = -9999.9d0                   ! [SS-SWC S-1.3/S-2.12B]
          if (sw%cofgen(10,node) > 0.0d0) sw%fluseksatexm(node) = .true.  ! [SS-SWC S-2.12B]
          sw%cofgen(11,node) = relsatthr(lay)             ! [SS-SWC S-1.3/S-2.12B]
          sw%cofgen(12,node) = ksatthr(lay)               ! [SS-SWC S-1.3/S-2.12B]
          if (iHWCKmodel(lay) ==  3 .OR. iHWCKmodel(lay) ==  6 .OR. iHWCKmodel(lay) ==  7 .OR. &
              iHWCKmodel(lay) == 10 .OR. iHWCKmodel(lay) == 11) then
             sw%cofgen(13:17,node) = paramvg(13:17,lay)   ! [SS-SWC S-1.3/S-2.12B]
          end if
          if (iHWCKmodel(lay) ==  5 .OR. iHWCKmodel(lay) ==  7) then
             sw%cofgen(18,node) = paramvg(18,lay)         ! [SS-SWC S-1.3/S-2.12B]
          end if
          if (iHWCKmodel(lay) ==  8 .OR. iHWCKmodel(lay) ==  9 .OR. &
              iHWCKmodel(lay) == 10 .OR. iHWCKmodel(lay) == 11) then
             sw%cofgen(18:21,node) = paramvg(18:21,lay)   ! [SS-SWC S-1.3/S-2.12B]
          end if
```

Append immediately after the last `end if`:

```fortran
          ! [SS-GR-UTILS] Mirror analytical cofgen writes into typed vg_params
          sw%vg_params(node)%thetar          = sw%cofgen(1,node)
          sw%vg_params(node)%thetas          = sw%cofgen(2,node)
          sw%vg_params(node)%ksat            = sw%cofgen(3,node)
          sw%vg_params(node)%alpha           = sw%cofgen(4,node)
          sw%vg_params(node)%lpar            = sw%cofgen(5,node)
          sw%vg_params(node)%npar            = sw%cofgen(6,node)
          sw%vg_params(node)%mpar            = sw%cofgen(7,node)
          sw%vg_params(node)%alphaw_sentinel = sw%cofgen(8,node)
          sw%vg_params(node)%h_enpr          = sw%cofgen(9,node)
          sw%vg_params(node)%ksatexm         = sw%cofgen(10,node)
          sw%vg_params(node)%relsatthr       = sw%cofgen(11,node)
          sw%vg_params(node)%ksatthr         = sw%cofgen(12,node)
          sw%vg_params(node)%alpha_2         = sw%cofgen(13,node)
          sw%vg_params(node)%npar_2          = sw%cofgen(14,node)
          sw%vg_params(node)%mpar_2          = sw%cofgen(15,node)
          sw%vg_params(node)%omega_1         = sw%cofgen(16,node)
          sw%vg_params(node)%omega_2         = sw%cofgen(17,node)
          sw%vg_params(node)%h0              = sw%cofgen(18,node)
          sw%vg_params(node)%ha              = sw%cofgen(19,node)
          sw%vg_params(node)%apar            = sw%cofgen(20,node)
          sw%vg_params(node)%omega_k         = sw%cofgen(21,node)
```

The "mirror from cofgen" pattern is intentional — even though some fields are conditional (only set under specific iHWCKmodel branches), reading from `sw%cofgen` ensures both data sources hold identical bits.

- [ ] **Step 4: Build + verify.**

Run: `rm -rf builddir && pixi run build-linux && pixi run test-pfunit && pixi run -e test python tests/regression/test_output_regression.py hupselbrook surfacewater salinitystress grassgrowth`

All gates green; regression 4/4 byte-for-byte. No behavior change (vg_params still unused).

- [ ] **Step 5: Commit.**

```bash
git add src/soil/soilhydraulics.f90
git commit -m "$(cat <<'EOF'
refactor(gr-utils): transitional dual-write — vg_params populated alongside cofgen

SoilHydraulics(task=1) now populates state%soilwater%vg_params(:)
in parallel with state%soilwater%cofgen(:,:), via field-by-field
mirroring. Both data sources hold identical bits.

Subsequent tasks migrate utility-function readers from cofgen to
vg_params; once all readers are migrated, cofgen retires
(Task 15).

Regression: 4/4 byte-for-byte (vg_params still unread).

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Add the other 9 + 6 utils-imported state fields

**Files:**
- Modify: `src/state/soilwater_state.f90` (add 9 fields)
- Modify: `src/state/surfacewater_state.f90` (add 6 fields)
- Modify: `src/soil/soilhydraulics.f90` (populate the soilwater fields from existing globals as a transitional dual-write)
- Modify: `src/io/toml/config_to_variables.f90` (populate the surfacewater fields)

This task adds the additional bare globals consumed by the util files as state fields. Population uses transitional dual-writes from existing init paths. No behavior change.

- [ ] **Step 1: Add soilwater fields.**

In `src/state/soilwater_state.f90`, just after the new `vg_params` line added in Task 2:

```fortran
      type(vanGenuchten_params_t), allocatable :: vg_params(:)   !< [SS-GR-UTILS] typed VG parameters, one per node
```

Add this block:

```fortran
      ! [SS-GR-UTILS] Soil hydraulic property metadata (migrated from variables.f90)
      ! Populated by SoilHydraulics(1) alongside cofgen/vg_params.
      integer                       :: swsophy    = 0      !< 0=analytical, 1=tabulated
      integer,         allocatable  :: numtab(:)            !< per-node table entry count
      real(real64),    allocatable  :: sptab(:,:,:)         !< soil property table — shape verbatim from variables.f90 (7, macp, matab)
      integer,         allocatable  :: ientrytab(:,:)       !< per-node entry-table indices
      integer,         allocatable  :: iHWCKmodel(:)        !< per-layer hydraulic-K model selector
      integer,         allocatable  :: layer(:)             !< per-node soil-layer index
      integer                       :: swfrost    = 0      !< frost-reduction simulation switch
      logical,         allocatable  :: BiModal(:)           !< per-layer bi-modal flag (used by WC_K_models_04_11)
      logical,         allocatable  :: NoVap(:)             !< per-layer no-vapor flag (used by WC_K_models_04_11)
```

Confirm exact rank/kind by reading `variables.f90` for these symbols first:

```bash
grep -n "swsophy\|numtab\|sptab\|ientrytab\|iHWCKmodel\|layer\|swfrost\|BiModal\|NoVap" src/core/variables.f90
```

Adjust the new state declarations to match the legacy rank exactly. Particularly `sptab`'s rank-3 shape (likely `(7, macp, matab)`).

- [ ] **Step 2: Allocate soilwater fields in soilwater_init.**

In `src/state/soilwater_state.f90`'s `soilwater_init` subroutine (around line 360-380), add after the existing allocations:

```fortran
      allocate(sw%numtab(numnod));        sw%numtab     = 0
      allocate(sw%sptab(7, numnod, ...))  ! match the legacy 3rd-dimension upper bound
      sw%sptab = 0.0_real64
      ! ientrytab: similar — check legacy shape, allocate same.
      allocate(sw%iHWCKmodel(nlay));      sw%iHWCKmodel = 0
      allocate(sw%layer(numnod));         sw%layer      = 0
      allocate(sw%BiModal(nlay));         sw%BiModal    = .false.
      allocate(sw%NoVap(nlay));           sw%NoVap      = .false.
```

The `sptab` 3rd dimension upper bound (likely `matabentries` or `matab`) needs to be retrieved from the existing global — verify with grep against `variables.f90`. If `matab` is a module-level parameter, fine; if it's a runtime-set variable, the caller has to provide the value. Inspect `src/soil/soilhydraulics.f90:830-843` for the existing tabulated-branch allocation to see how the legacy code handles it.

- [ ] **Step 3: Populate soilwater fields in SoilHydraulics(1).**

In `src/soil/soilhydraulics.f90`, locate the `swsophy == 1` branch (lines 830-843). Add after the existing tabulated copy:

```fortran
            ! [SS-GR-UTILS] Mirror tabulated metadata into state
            sw%swsophy = swsophy
            sw%numtab(node) = numtab(node)
            ! ientrytab values mirror what was just written above
            ! sptab values mirror the loops 838-843
            ! Copy via:
            sw%ientrytab(node, :) = ientrytab(node, :)
            sw%sptab(:, node, :) = sptab(:, node, :)
            sw%iHWCKmodel(:) = iHWCKmodel(:)
            sw%layer(node) = layer(node)
```

The implementer adapts the exact assignment to match the legacy structure (e.g., assigning `iHWCKmodel(:)` once outside the `do node` loop is more efficient).

For `swfrost`, `BiModal`, `NoVap`: these are populated by `config_to_variables` (not `SoilHydraulics`). Add the mirror writes there. Inspect `src/io/toml/config_to_variables.f90` for the swfrost/BiModal/NoVap write sites and add `state%soilwater%swfrost = swfrost` etc. alongside.

Actually `config_to_variables` already takes state per SS-BMI2 Task 5. Look at how it writes `state%timecontrol%X` and follow the same pattern for these soilwater fields. Search for `state%timecontrol%tend` in `config_to_variables.f90` to find the established write style.

- [ ] **Step 4: Add surfacewater fields.**

In `src/state/surfacewater_state.f90`, find the `type :: surfacewater_state_t` block and add (near the bottom, before `contains` if present):

```fortran
      ! [SS-GR-UTILS] Surface-water utils config (migrated from variables.f90)
      real(real64),    allocatable  :: hqhtab(:)            !< Q-h table head entries
      real(real64),    allocatable  :: qqhtab(:)            !< Q-h table discharge entries
      integer                       :: swdra      = 0      !< drainage switch (0=none, 1=drainage, 2=surface water)
      real(real64)                  :: pondmx     = 0.0_real64 !< max ponding depth (cm)
      real(real64)                  :: rsro       = 0.0_real64 !< runoff resistance (d)
      real(real64)                  :: rsroexp    = 0.0_real64 !< runoff exponent
```

Confirm the rank of `hqhtab` and `qqhtab` against `variables.f90` first — likely allocatable rank-1 arrays.

- [ ] **Step 5: Allocate + populate surfacewater fields in state%surfacewater%init or config_to_variables.**

`state%surfacewater%init` already exists (the SS-* pilot). Check whether `hqhtab` / `qqhtab` are allocated there or in `config_to_variables`. Adopt the existing pattern. The 4 scalars (`swdra`, `pondmx`, `rsro`, `rsroexp`) are simple assignments from `config%surface_water%*` or similar config paths.

If the legacy global `swdra` is read from `config%drain%swdra` (verify via grep in `config_to_variables.f90`), mirror that into `state%surfacewater%swdra`.

- [ ] **Step 6: Build + verify all gates.**

```
rm -rf builddir
pixi run build-linux
pixi run test-pfunit
pixi run -e test python tests/regression/test_output_regression.py hupselbrook surfacewater salinitystress grassgrowth
```

All green. 4/4 regression byte-for-byte.

- [ ] **Step 7: Commit.**

```bash
git add src/state/soilwater_state.f90 src/state/surfacewater_state.f90 \
        src/soil/soilhydraulics.f90 src/io/toml/config_to_variables.f90
git commit -m "$(cat <<'EOF'
schema(gr-utils): add 9 soilwater + 6 surfacewater fields

state%soilwater gains: swsophy, numtab, sptab, ientrytab, iHWCKmodel,
layer, swfrost, BiModal, NoVap. state%surfacewater gains: hqhtab,
qqhtab, swdra, pondmx, rsro, rsroexp.

All populated as transitional dual-writes from existing
SoilHydraulics(1) / config_to_variables sources. Bare globals
still live in variables.f90 (other readers retire in later arcs).

Regression: 4/4 byte-for-byte. No behavior change.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Migrate `watcon` to typed signature

**Files:**
- Modify: `src/utils/soilhydraulicsutils.f90` (change `watcon` signature + body)
- Modify: 9 caller files (every site listed below)

`watcon`'s new signature: `watcon(head, vg, model)`. The `swsophy==1` tabulated branch needs `numtab` + `sptab` access — handled by adding a SECOND signature `watcon_tabulated(head, soilwater, node)` OR (chosen for simplicity here) by taking `soilwater_state` as an additional optional arg. Implementer picks; the simpler approach:

**Decision: `watcon` takes 3 dummy args** in the analytical case AND falls back to legacy behavior for `swsophy==1`. The tabulated branch is small (~10 lines inside `watcon`). Restructure as:

```fortran
function watcon(head, vg, model) result(theta)
   real(real64),                  intent(in) :: head
   type(vanGenuchten_params_t),   intent(in) :: vg
   integer,                       intent(in) :: model   ! iHWCKmodel value for this node
   real(real64) :: theta
   ...
end function
```

For the tabulated case (`swsophy == 1`), callers explicitly invoke a NEW separate function `watcon_tabulated`. This makes the dispatch explicit at the call site (matching the spec's preference for self-documenting signatures).

But for THIS task, focus on the analytical case only — the tabulated branch stays in `watcon` reading from module-level globals temporarily, OR the function takes a NEW optional `soilwater` arg for the tabulated case. To minimize blast radius, this task does the following:

**Simplified Task 5 scope:** `watcon` takes `(head, vg, model)` for the analytical path. For the tabulated path, it ALSO takes `soilwater` as an optional arg (read access to `swsophy`, `numtab`, `sptab`). Callers in the analytical path don't pass `soilwater`; callers in the tabulated path do.

- [ ] **Step 1: Read the current `watcon`.**

Read `src/utils/soilhydraulicsutils.f90:96-260` (the full function body).

- [ ] **Step 2: Rewrite the signature + body.**

Replace the current `watcon` function with the new version:

```fortran
   !> Calculate water content from pressure head
   function watcon(head, vg, model, soilwater) result(theta)
      use soilwater_state_mod, only: soilwater_state_t
      use hydraulic_params_mod, only: vanGenuchten_params_t
      implicit none

      real(real64),                  intent(in)           :: head
      type(vanGenuchten_params_t),   intent(in)           :: vg
      integer,                       intent(in)           :: model
      type(soilwater_state_t),       intent(in), optional :: soilwater
      real(real64) :: theta

      ! Local variables (unchanged from legacy)
      real(real64) :: h_enpr, help, m, n, s_enpr
      real(real64), parameter :: h_crit = -1.0d-2
      real(real64) :: h105, C105, a, b, dum, t105
      real(real64) :: alfamg, thetar, thetas
      real(real64) :: alfa_2, n_2, m_2, omega_1

      ! Branch on swsophy via the optional soilwater arg.
      ! Analytical branch (swsophy==0) uses vg directly; tabulated branch
      ! (swsophy==1) reads from soilwater%sptab/numtab. If soilwater is
      ! absent, default to analytical.
      if (present(soilwater)) then
         if (soilwater%swsophy == 1) then
            ! Tabulated branch — same as legacy body for swsophy==1.
            ! Read sptab/numtab from the soilwater dummy.
            ! ... port legacy lines 207-260 here, replacing:
            !       sptab(...,node,...) with soilwater%sptab(...,node,...)
            !       numtab(node)        with soilwater%numtab(node)
            ! Note: this branch still needs `node`. Pass it via the vg arg's
            ! position? No — vg has no `node` info. Instead, take an
            ! additional `node` arg, OR (cleaner) restructure callers to
            ! call a separate watcon_tabulated function.
            !
            ! For Phase 2 of this task: add an explicit `node` arg to
            ! watcon, so the tabulated branch can index sptab(...,node,...).
            !
            ! DECISION: add `node` as an OPTIONAL argument that defaults
            ! to 1 (sentinel — fails if used incorrectly). Callers in the
            ! tabulated path pass `node`. Callers in the analytical path
            ! omit it.
            return
         end if
      end if

      ! Analytical branch (swsophy == 0)
      thetar = vg%thetar
      thetas = vg%thetas
      alfamg = vg%alpha
      n      = vg%npar
      m      = vg%mpar
      h_enpr = vg%h_enpr

      if (model == 2) then
         ! Exponential
         theta = dmax1(1.0000001_real64*thetar, &
                       thetar + (thetas-thetar)*dexp(alfamg*head))
      else if (model == 3) then
         ! Bi-modal
         alfa_2  = vg%alpha_2
         n_2     = vg%npar_2
         m_2     = vg%mpar_2
         omega_1 = vg%omega_1
         if (head < 0.0_real64) then
            theta = omega_1 / (1.0_real64 + (dabs(alfamg*head))**n)**m
            theta = theta + (1.0_real64 - omega_1) / (1.0_real64 + (dabs(alfa_2*head))**n_2)**m_2
            theta = thetar + (thetas - thetar) * theta
         else
            theta = thetas
         end if
      else if (model > 3 .and. model < 12) then
         ! WC_K_models 4..11
         ! Note: functionvalue_04_11 is migrated in Task 9 — for now it
         ! still relies on the module-level cofgen pointer.
         theta = functionvalue_04_11(1, [node_placeholder], head)
         ! ↑ This call needs `node` — to avoid the placeholder, the
         ! implementer must thread `node` through or restructure.
         ! See implementer note below.
      else
         ! Default MvG with optional air-entry pressure modification
         ! ... port legacy lines 149-203 here, replacing every cofgen-derived
         ! local with the equivalent vg% field that was assigned above
         ! (thetar, thetas, alfamg, n, m, h_enpr).
      end if
   end function watcon
```

**Implementer reality check:** the cleanest split is to add `node` as an additional dummy arg to `watcon`, eliminating the optional pattern. Sketch:

```fortran
function watcon(head, vg, model, node, soilwater) result(theta)
   real(real64),                intent(in) :: head
   type(vanGenuchten_params_t), intent(in) :: vg
   integer,                     intent(in) :: model
   integer,                     intent(in) :: node       ! for tabulated branch + functionvalue_04_11
   type(soilwater_state_t),     intent(in) :: soilwater  ! for tabulated branch (swsophy/sptab/numtab)
```

The new signature is `(head, vg, model, node, soilwater)` — 5 args. Verbose but explicit. Use this form for the implementation.

- [ ] **Step 3: Update all callers.**

Search for `watcon(` callers via:

```bash
grep -rn "\\bwatcon\\s*(" src/ --include="*.f90" | grep -v "soilhydraulicsutils.f90"
```

The 8 caller files. At each site, transform:

```fortran
! before:
theta_val = watcon(node, head_val)

! after:
theta_val = watcon(head_val, &
                   state%soilwater%vg_params(node), &
                   state%soilwater%iHWCKmodel(state%soilwater%layer(node)), &
                   node, state%soilwater)
```

State is in scope at every site (callers were threaded with state in earlier SS-* arcs).

- [ ] **Step 4: Build + verify all gates.**

```
rm -rf builddir
pixi run build-linux
pixi run test-pfunit
pixi run -e test python tests/regression/test_output_regression.py hupselbrook surfacewater salinitystress grassgrowth
```

Regression must remain 4/4 byte-for-byte. If a case fails, the new `watcon` body either reads the wrong vg field or has a branch-condition mistake. Diff against the legacy carefully.

- [ ] **Step 5: Commit.**

```bash
git add src/utils/soilhydraulicsutils.f90 \
        src/soil/soilhydraulics.f90 src/soil/soilgrid.f90 \
        src/boundary/boundbottom.f90 src/boundary/boundtop.f90 \
        src/crop/tillage.f90 src/crop/irrigation.f90 \
        src/crop/cropgrowth.f90 src/crop/oxygenstress.f90 \
        src/crop/rootextraction.f90 \
        src/atmosphere/meteoday.f90
git commit -m "$(cat <<'EOF'
refactor(gr-utils): watcon — typed signature with vg_params

watcon now takes (head, vg, model, node, soilwater) — explicit
hydraulic parameter record (vg) replaces cofgen magic-row reads;
model = iHWCKmodel(layer(node)) is explicit; node + soilwater pass
through to the tabulated branch and (transitionally) to
functionvalue_04_11.

10+ call sites updated across crop/, boundary/, atmosphere/,
soil/. Module-level cofgen pointer still alive (used by Tasks 6-9
functions and by functionvalue_04_11 until Task 9). Bind retirement
happens in Task 13.

Regression: 4/4 byte-for-byte.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Migrate `hconduc` + `dhconduc` to typed signatures

**Files:**
- Modify: `src/utils/soilhydraulicsutils.f90`
- Modify: 4 caller files (`soilhydraulics.f90`, `boundtop.f90`, `boundbottom.f90`, `rootextraction.f90`)

`hconduc` and `dhconduc` both use cofgen / iHWCKmodel / layer / fluseksatexm. Migrate together since both are hydraulic-K functions.

- [ ] **Step 1: Rewrite `hconduc`.**

New signature:

```fortran
function hconduc(h, theta, rfcp, tsoil_node, vg, model, use_ksatexm, node, soilwater) result(k)
   real(real64),                intent(in) :: h, theta, rfcp, tsoil_node
   type(vanGenuchten_params_t), intent(in) :: vg
   integer,                     intent(in) :: model
   logical,                     intent(in) :: use_ksatexm
   integer,                     intent(in) :: node          ! for functionvalue_04_11
   type(soilwater_state_t),     intent(in) :: soilwater     ! for swsophy branch
   real(real64) :: k
```

Body replaces `cofgen(j, node)` with `vg%<field>`, `iHWCKmodel(layer(node))` with `model`, `fluseksatexm(node)` with `use_ksatexm`. The `swsophy` switch uses `soilwater%swsophy` for now (or is passed as an integer arg if cleaner).

- [ ] **Step 2: Rewrite `dhconduc` symmetrically.**

```fortran
function dhconduc(h, theta, dimoca, rfcp, vg, model, node, soilwater) result(dkdh)
   real(real64),                intent(in) :: h, theta, dimoca, rfcp
   type(vanGenuchten_params_t), intent(in) :: vg
   integer,                     intent(in) :: model
   integer,                     intent(in) :: node
   type(soilwater_state_t),     intent(in) :: soilwater
   real(real64) :: dkdh
```

- [ ] **Step 3: Update all `hconduc` callers.**

Search:

```bash
grep -rn "\\bhconduc\\s*(" src/ --include="*.f90" | grep -v "soilhydraulicsutils.f90"
```

Transform each site to pass the new args. State always in scope.

Example, `soilhydraulics.f90:151`:

```fortran
! before:
sw_k(i) = hconduc(i, sw_h(i), sw_theta(i), state%heat%rfcp(i), state%heat%tsoil(i))

! after:
sw_k(i) = hconduc(sw_h(i), sw_theta(i), state%heat%rfcp(i), state%heat%tsoil(i), &
                  state%soilwater%vg_params(i), &
                  state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                  state%soilwater%fluseksatexm(i), &
                  i, state%soilwater)
```

- [ ] **Step 4: Update `dhconduc` callers.**

Same pattern; only `soilhydraulics.f90` calls it.

- [ ] **Step 5: Build + verify all gates.** Same gate commands as Task 5 Step 4. 4/4 byte-for-byte.

- [ ] **Step 6: Commit.**

```bash
git add src/utils/soilhydraulicsutils.f90 \
        src/soil/soilhydraulics.f90 \
        src/boundary/boundbottom.f90 src/boundary/boundtop.f90 \
        src/crop/rootextraction.f90
git commit -m "$(cat <<'EOF'
refactor(gr-utils): hconduc + dhconduc — typed signatures

Both hydraulic-K functions now take (h, theta, rfcp, tsoil_node, vg,
model, use_ksatexm, node, soilwater) and (h, theta, dimoca, rfcp,
vg, model, node, soilwater) respectively. Explicit args replace
module-level cofgen pointer + bare-global iHWCKmodel/layer/
fluseksatexm reads.

Call sites updated in soilhydraulics, boundary, rootextraction.

Regression: 4/4 byte-for-byte.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Migrate `moiscap` to typed signature

**Files:**
- Modify: `src/utils/soilhydraulicsutils.f90` (rewrite `moiscap`)
- Modify: `src/soil/soilhydraulics.f90` (single caller site)

`moiscap` reads cofgen + iHWCKmodel + layer + tc_dt_ptr. New signature adds `dt` as a dummy arg.

- [ ] **Step 1: Rewrite `moiscap`.**

```fortran
function moiscap(h, vg, model, dt, node, soilwater) result(c)
   real(real64),                intent(in) :: h
   type(vanGenuchten_params_t), intent(in) :: vg
   integer,                     intent(in) :: model
   real(real64),                intent(in) :: dt
   integer,                     intent(in) :: node
   type(soilwater_state_t),     intent(in) :: soilwater
   real(real64) :: c
```

Body replaces every `cofgen(j, node)` with `vg%<field>` and `tc_dt_ptr` with `dt`.

- [ ] **Step 2: Update the caller in `soilhydraulics.f90`.**

```fortran
! before:
sw_dimoca(i) = moiscap(i, sw_h(i))

! after:
sw_dimoca(i) = moiscap(sw_h(i), &
                       state%soilwater%vg_params(i), &
                       state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                       state%timecontrol%dt, &
                       i, state%soilwater)
```

- [ ] **Step 3: Build + verify all gates.** 4/4.

- [ ] **Step 4: Commit.**

```bash
git add src/utils/soilhydraulicsutils.f90 src/soil/soilhydraulics.f90
git commit -m "$(cat <<'EOF'
refactor(gr-utils): moiscap — typed signature + explicit dt

moiscap(h, vg, model, dt, node, soilwater). The dt arg replaces the
tc_dt_ptr module-level pointer (which becomes dead code after
Task 13's bind cleanup).

Regression: 4/4 byte-for-byte.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Migrate `prhead` to typed signature

**Files:**
- Modify: `src/utils/soilhydraulicsutils.f90` (rewrite `prhead`)
- Modify: `src/soil/soilgrid.f90` (single non-trivial caller — passes a tentative `cofgenNew` array during grid redistribution)

`prhead` previously took `cofgen_in` as an explicit arg (the soilgrid caller passes `cofgenNew`, a temporary built during regridding). The new signature replaces `cofgen_in` with `vg_in`, an optional `vanGenuchten_params_t` argument that defaults to `soilwater%vg_params(node)`.

- [ ] **Step 1: Rewrite `prhead`.**

```fortran
function prhead(disnod, theta, h, model, node, soilwater, vg_in) result(h_out)
   real(real64),                intent(in)           :: disnod(:)
   real(real64),                intent(in)           :: theta
   real(real64),                intent(in)           :: h(:)
   integer,                     intent(in)           :: model
   integer,                     intent(in)           :: node
   type(soilwater_state_t),     intent(in)           :: soilwater
   type(vanGenuchten_params_t), intent(in), optional :: vg_in
   real(real64) :: h_out
   type(vanGenuchten_params_t) :: vg
   if (present(vg_in)) then
      vg = vg_in
   else
      vg = soilwater%vg_params(node)
   end if
   ! ... port legacy body, replacing cofgen_in(j, node) with vg%<field> ...
end function
```

- [ ] **Step 2: Update the soilgrid caller (`src/soil/soilgrid.f90:388`).**

The caller currently builds `cofgenNew(:,:)` as a temporary and passes it. The new pattern: build a single `type(vanGenuchten_params_t)` for the new node and pass via `vg_in`.

Read `src/soil/soilgrid.f90:370-400` to understand the redistribution loop. The transform:

```fortran
! before:
cofgenNew(:, node) = sw_cofgen(:, NodeNew(node, 2))      ! line 385
hNew(node) = prhead(node, disnodNew(node), thetaNew(node), cofgenNew, hNew)   ! line 388

! after:
type(vanGenuchten_params_t) :: vg_new
vg_new = state%soilwater%vg_params(NodeNew(node, 2))     ! copy from source node
hNew(node) = prhead(disnodNew(node), thetaNew(node), hNew, &
                    state%soilwater%iHWCKmodel(state%soilwater%layer(node)), &
                    node, state%soilwater, vg_in=vg_new)
```

The `cofgenNew` 2D temporary array can be retired entirely — replaced by a per-call `vg_new` scalar (of derived type). If the surrounding code uses `cofgenNew` for other purposes (verify with grep), keep it as a transient for those uses but stop passing it to `prhead`.

- [ ] **Step 3: Update the other `prhead` caller** in `src/crop/tillage.f90:291` (per the user's IDE selection — verify with grep).

Similar transform, but without the temporary `vg_new`:

```fortran
h_value = prhead(disnod_arr, theta_val, h_arr, &
                 state%soilwater%iHWCKmodel(state%soilwater%layer(node)), &
                 node, state%soilwater)
```

(omit `vg_in` to use the default).

- [ ] **Step 4: Build + verify all gates.** 4/4.

- [ ] **Step 5: Commit.**

```bash
git add src/utils/soilhydraulicsutils.f90 src/soil/soilgrid.f90 src/crop/tillage.f90
git commit -m "$(cat <<'EOF'
refactor(gr-utils): prhead — typed signature + optional vg_in

prhead(disnod, theta, h, model, node, soilwater, vg_in). The
optional vg_in handles soilgrid's grid-redistribution use case
(passes the vg from the source node directly). Otherwise defaults
to soilwater%vg_params(node).

cofgenNew 2D temporary in soilgrid retired (replaced by scalar
vg_new of derived type).

Regression: 4/4 byte-for-byte.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: Migrate `functionvalue_04_11` in `WC_K_models_04_11`

**Files:**
- Modify: `src/soil/WC_K_models_04_11.f90` (rewrite the function + delete bind_cofgen_target + delete `use variables`)
- Modify: `src/utils/soilhydraulicsutils.f90` (update call to `functionvalue_04_11`)

`functionvalue_04_11` uses module-level cofgen pointer + bare globals `iHWCKmodel`, `BiModal`, `NoVap`, `layer`. Migrate to take explicit args.

- [ ] **Step 1: Rewrite `functionvalue_04_11`.**

New signature:

```fortran
function functionvalue_04_11(iType, h, vg, model, is_bimodal, no_vap, wc, temp) result(val)
   integer,                     intent(in)           :: iType
   real(real64),                intent(in)           :: h
   type(vanGenuchten_params_t), intent(in)           :: vg
   integer,                     intent(in)           :: model
   logical,                     intent(in)           :: is_bimodal
   logical,                     intent(in)           :: no_vap
   real(real64),                intent(in), optional :: wc, temp
   real(real64) :: val
```

Body replaces `cofgen(j, iNode)` with `vg%<field>`, `iHWCKmodel(layer(iNode))` with `model`, `BiModal(layer(iNode))` with `is_bimodal`, `NoVap(layer(iNode))` with `no_vap`. `iNode` retires from the signature — callers pre-slice.

- [ ] **Step 2: Delete `bind_cofgen_target`.**

In `src/soil/WC_K_models_04_11.f90`, delete:

```fortran
! gone — module-level pointer
real(real64), pointer :: cofgen(:,:) => null()

! gone — bind setup
subroutine bind_cofgen_target(sw_cofgen_in)
   real(real64), target, intent(in) :: sw_cofgen_in(:,:)
   cofgen => sw_cofgen_in
end subroutine bind_cofgen_target

! also drop from public list:
public :: functionvalue_04_11, bind_cofgen_target
! → public :: functionvalue_04_11
```

- [ ] **Step 3: Delete `use variables` from `WC_K_models_04_11.f90`.**

Replace:

```fortran
use variables, only: iHWCKmodel, BiModal, NoVap, layer
```

with nothing (delete the line).

The fields are now read via the function args, not via module imports.

- [ ] **Step 4: Update the call inside `watcon` (in soilhydraulicsutils).**

Find `functionvalue_04_11` calls in `src/utils/soilhydraulicsutils.f90`. Transform:

```fortran
! before:
watcon = functionvalue_04_11(1, node, head)

! after:
theta = functionvalue_04_11(1, head, vg, model, &
                            soilwater%BiModal(soilwater%layer(node)), &
                            soilwater%NoVap(soilwater%layer(node)))
```

`watcon` already has `vg`, `model`, `node`, `soilwater` in scope (from Task 5).

- [ ] **Step 5: Build + verify all gates.** 4/4.

- [ ] **Step 6: Commit.**

```bash
git add src/soil/WC_K_models_04_11.f90 src/utils/soilhydraulicsutils.f90
git commit -m "$(cat <<'EOF'
refactor(gr-utils): functionvalue_04_11 — typed signature; bind_cofgen_target deleted

functionvalue_04_11(iType, h, vg, model, is_bimodal, no_vap, wc, temp).
Module-level cofgen pointer + bind_cofgen_target setup procedure
deleted. `use variables` dropped from WC_K_models_04_11 entirely.

watcon's internal call updated.

Regression: 4/4 byte-for-byte.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: Update inline cofgen readers in `soilhydraulics.f90`

**Files:**
- Modify: `src/soil/soilhydraulics.f90`

The file at lines 845-877 + 851-852 reads `sw%cofgen(j, node)` directly in non-utils code (it's the population block from Task 3 — those WRITES stay because cofgen is still the data source for any unmigrated reader, but the surrounding code also READS cofgen for thetas/ksat/etc.).

- [ ] **Step 1: Find inline reads.** Run:

```bash
grep -n "sw%cofgen\\b\\|state%soilwater%cofgen\\b" src/soil/soilhydraulics.f90
```

Identify each line that READS cofgen (vs. WRITES it). Examples from the earlier inventory:

- Line 851: `ksatfit(lay) = sw%cofgen(3,nod1lay(lay))` — reads ksat for first node of layer.
- Line 852: `sw%thetsl(lay) = sw%cofgen(2,nod1lay(lay))` — reads thetas.
- Line 863: `if (sw%cofgen(10,node) > 0.0d0) sw%fluseksatexm(node) = .true.` — reads ksatexm.

- [ ] **Step 2: Update each reader.**

Transform:

```fortran
! before:
ksatfit(lay) = sw%cofgen(3, nod1lay(lay))
sw%thetsl(lay) = sw%cofgen(2, nod1lay(lay))

! after:
ksatfit(lay) = sw%vg_params(nod1lay(lay))%ksat
sw%thetsl(lay) = sw%vg_params(nod1lay(lay))%thetas
```

```fortran
! before:
if (sw%cofgen(10,node) > 0.0d0) sw%fluseksatexm(node) = .true.

! after:
if (sw%vg_params(node)%ksatexm > 0.0d0) sw%fluseksatexm(node) = .true.
```

- [ ] **Step 3: Build + verify all gates.** 4/4.

- [ ] **Step 4: Commit.**

```bash
git add src/soil/soilhydraulics.f90
git commit -m "$(cat <<'EOF'
refactor(gr-utils): soilhydraulics inline cofgen readers → vg_params

Lines 845-877 in SoilHydraulics(1) still WRITE both cofgen and
vg_params (transitional). The inline READERS (ksatfit assignment,
thetsl assignment, fluseksatexm check) now read from vg_params.

Regression: 4/4 byte-for-byte.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 11: Update cofgenNew redistribution in `soilgrid.f90` if any reads remain

**Files:**
- Modify: `src/soil/soilgrid.f90`

Task 8 handled the `prhead` caller. This task cleans up any remaining `cofgenNew` references and the `sw_cofgen` associate alias.

- [ ] **Step 1: Find cofgen references in soilgrid.**

```bash
grep -n "cofgen" src/soil/soilgrid.f90
```

- [ ] **Step 2: Migrate each.**

If `sw_cofgen => state%soilwater%cofgen` associate alias is no longer used after Task 8 + this task, remove it. Replace any remaining `sw_cofgen(...)` reads with `state%soilwater%vg_params(...)%<field>`.

`cofgenNew` array can be deleted entirely if it's no longer referenced after Task 8.

- [ ] **Step 3: Build + verify all gates.** 4/4.

- [ ] **Step 4: Commit.**

```bash
git add src/soil/soilgrid.f90
git commit -m "$(cat <<'EOF'
refactor(gr-utils): soilgrid — finish cofgen → vg_params migration

Removes the sw_cofgen associate alias and any remaining inline
cofgenNew references. Grid redistribution now operates on
vg_params(:) directly (per Task 8's prhead caller update).

Regression: 4/4 byte-for-byte.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 12: Migrate `surfacewaterutils`

**Files:**
- Modify: `src/utils/surfacewaterutils.f90` (drop `use variables`, update bodies)
- Possibly modify: `src/drainage/surfacewater.f90` and `src/boundary/boundtop.f90` if `qhtab` needs a state arg

`swstlev` and `runoff` already take state. `wlevst(state, swstor)` takes state. `qhtab(wlev, imper_in)` does NOT — needs threading.

- [ ] **Step 1: Drop `use variables` from surfacewaterutils.**

In `src/utils/surfacewaterutils.f90`:

```fortran
! before:
use variables, only: hqhtab, qqhtab, swdra, pondmx, rsro, rsroexp

! after:
! (line deleted)
```

- [ ] **Step 2: Update each function body** to read from `state%surfacewater%<field>` instead of bare names.

For `runoff` (state in scope):
```fortran
! before:
! ... swdra ... pondmx ... rsro ... rsroexp ...

! after:
! ... state%surfacewater%swdra ... state%surfacewater%pondmx ... etc.
```

For `swstlev` (state in scope): no `use variables` reads to migrate (it uses `state%surfacewater%sttab` already).

For `wlevst` (state in scope): check the body for any `use variables` symbol reads, migrate to state%surfacewater%X.

For `qhtab(wlev, imper_in)` — currently doesn't take state. Two paths:
- **Option A:** add `state%surfacewater` as a new dummy arg. Update the caller in `surfacewater.f90:296` to pass `state%surfacewater`.
- **Option B:** add `hqhtab` and `qqhtab` as explicit dummy args. Caller passes `state%surfacewater%hqhtab, state%surfacewater%qqhtab`.

Option A is consistent with the rest of the module. Choose Option A.

- [ ] **Step 3: Update qhtab signature + caller.**

```fortran
! before:
function qhtab(wlev, imper_in)

! after:
function qhtab(state_sw, wlev, imper_in)
   use surfacewater_state_mod, only: surfacewater_state_t
   type(surfacewater_state_t), intent(in) :: state_sw
   real(real64),               intent(in) :: wlev
   integer,                    intent(in) :: imper_in
   ! ... body reads state_sw%hqhtab(...) etc.
```

Caller in `src/drainage/surfacewater.f90`:

```fortran
! before:
qh = qhtab(wlev, imper)

! after:
qh = qhtab(state%surfacewater, wlev, imper)
```

- [ ] **Step 4: Build + verify all gates.** 4/4.

- [ ] **Step 5: Commit.**

```bash
git add src/utils/surfacewaterutils.f90 src/drainage/surfacewater.f90
git commit -m "$(cat <<'EOF'
refactor(gr-utils): surfacewaterutils — drop use variables

All 6 imports (hqhtab, qqhtab, swdra, pondmx, rsro, rsroexp) replaced
by reads from state%surfacewater fields populated by Task 4.

qhtab signature gained a state_sw arg (was the only function not
already taking state). Caller in drainage/surfacewater.f90 updated.

Regression: 4/4 byte-for-byte.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 13: Drop `use variables` + delete module-level pointers + bind_*_target setup procedures in `soilhydraulicsutils`

**Files:**
- Modify: `src/utils/soilhydraulicsutils.f90`

After Tasks 5-9, no function in `soilhydraulicsutils` references the module-level pointers or the variables-imported globals anymore. Time to delete them.

- [ ] **Step 1: Delete the `use variables` line.**

```fortran
! gone:
use variables, only: swsophy, numtab, sptab, ientrytab, &
                     iHWCKmodel, layer, swfrost
```

- [ ] **Step 2: Delete the module-level pointers.**

```fortran
! gone:
real(real64), pointer :: cofgen(:,:) => null()
logical,      pointer :: fluseksatexm(:) => null()
real(real64), pointer :: tc_dt_ptr => null()
```

- [ ] **Step 3: Delete `bind_state_targets` and `bind_tc_target`.**

```fortran
! gone — whole subroutines:
subroutine bind_state_targets(sw_cofgen_in, sw_fluseksatexm_in)
   ...
end subroutine bind_state_targets

subroutine bind_tc_target(tc_dt_in)
   ...
end subroutine bind_tc_target
```

Update the module's `public ::` declaration to drop both names.

- [ ] **Step 4: Build + verify all gates.** 4/4.

If build fails with unresolved references, a Task 5-9 reader still uses a deleted pointer — find and migrate.

- [ ] **Step 5: Commit.**

```bash
git add src/utils/soilhydraulicsutils.f90
git commit -m "$(cat <<'EOF'
refactor(gr-utils): soilhydraulicsutils — drop use variables, retire bind pattern

Three module-level pointers (cofgen, fluseksatexm, tc_dt_ptr) and
the two bind setup procedures (bind_state_targets, bind_tc_target)
deleted. The seven `use variables` imports also deleted.

All functions now operate on explicit dummy args — no hidden state.
soilhydraulicsutils is a pure utility module.

Regression: 4/4 byte-for-byte.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 14: Delete bind_*_target calls + imports in `swap_mod.f90`

**Files:**
- Modify: `src/core/swap_mod.f90`

The bind setup at lines 73-74 + 111-113 are now dead code.

- [ ] **Step 1: Delete the 2 `use` imports.**

```fortran
! gone (line 73): use WC_K_models_04_11, only: bind_cofgen_target
! gone (line 74): use soilhydraulics_utils, only: bind_state_targets, bind_tc_target
```

- [ ] **Step 2: Delete the 3 `call bind_*` lines.**

```fortran
! gone (line 111): call bind_cofgen_target(state%soilwater%cofgen)
! gone (line 112): call bind_state_targets(state%soilwater%cofgen, state%soilwater%fluseksatexm)
! gone (line 113): call bind_tc_target(state%timecontrol%dt)
```

Also delete the surrounding `! [SS-SWC S-2.12B] bind module-level pointers...` comment block (3 lines).

- [ ] **Step 3: Build + verify all gates.** 4/4.

- [ ] **Step 4: Commit.**

```bash
git add src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
refactor(gr-utils): retire bind_*_target calls from swap_mod

3 `call bind_*_target(...)` lines + 2 `use ..., only: bind_*`
imports deleted (lines 73-74, 111-113 of swap_mod.f90, plus the
surrounding S-2.12B comment block).

The strangler-fig leftover for soilwater binding is fully gone.

Regression: 4/4 byte-for-byte.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 15: Delete `cofgen(:,:)` field from `soilwater_state_t`

**Files:**
- Modify: `src/state/soilwater_state.f90` (delete the field)
- Modify: `src/soil/soilhydraulics.f90` (drop the transitional cofgen writes added in Task 3)

All readers have migrated. The transitional dual-write is no longer needed.

- [ ] **Step 1: Verify no live cofgen readers remain.**

```bash
grep -rn "state%soilwater%cofgen\\|sw%cofgen" src/ --include="*.f90"
```

Expected hits:
- The DECLARATION in `src/state/soilwater_state.f90` (about to delete).
- The WRITES in `src/soil/soilhydraulics.f90` (about to delete in Step 3).

If ANY other line shows a live READ or WRITE, fix it first (migrate to vg_params).

- [ ] **Step 2: Delete the field from `soilwater_state_t`.**

In `src/state/soilwater_state.f90`:

```fortran
! gone:
real(real64), allocatable :: cofgen(:,:)     !< Mualem-VG parameters (21 × numnod) — [SS-GR-UTILS] retiring in this arc
```

Also delete the allocation line in `soilwater_init`:

```fortran
! gone:
allocate(sw%cofgen(21, numnod));   sw%cofgen       = 0.0_real64
```

- [ ] **Step 3: Drop transitional cofgen writes in `soilhydraulics.f90`.**

In `src/soil/soilhydraulics.f90:830-880`, delete every `sw%cofgen(...) = ...` line. The `paramvg` table is now only used to populate `sw%vg_params(node)%<field>` — which means the population needs to switch from "fill cofgen via paramvg(1..21, lay), then mirror to vg_params" to "fill vg_params directly from paramvg".

Restructure the population block:

```fortran
! before (analytical branch, after Tasks 3 and 10):
do i = 1, 10
   sw%cofgen(i,node) = paramvg(i,lay)
end do
sw%cofgen(8,node) = -9999.9d0
if (sw%cofgen(10,node) > 0.0d0) sw%fluseksatexm(node) = .true.
sw%cofgen(11,node) = relsatthr(lay)
sw%cofgen(12,node) = ksatthr(lay)
if (iHWCKmodel(lay) == 3 .OR. ...) then
   sw%cofgen(13:17,node) = paramvg(13:17,lay)
end if
! ... etc ...
sw%vg_params(node)%thetar = sw%cofgen(1,node)
sw%vg_params(node)%thetas = sw%cofgen(2,node)
! ... etc ...

! after (cofgen retired):
sw%vg_params(node)%thetar           = paramvg(1, lay)
sw%vg_params(node)%thetas           = paramvg(2, lay)
sw%vg_params(node)%ksat             = paramvg(3, lay)
sw%vg_params(node)%alpha            = paramvg(4, lay)
sw%vg_params(node)%lpar             = paramvg(5, lay)
sw%vg_params(node)%npar             = paramvg(6, lay)
sw%vg_params(node)%mpar             = paramvg(7, lay)
sw%vg_params(node)%alphaw_sentinel  = -9999.9_real64
sw%vg_params(node)%h_enpr           = paramvg(9, lay)
sw%vg_params(node)%ksatexm          = paramvg(10, lay)
if (sw%vg_params(node)%ksatexm > 0.0_real64) sw%fluseksatexm(node) = .true.
sw%vg_params(node)%relsatthr        = relsatthr(lay)
sw%vg_params(node)%ksatthr          = ksatthr(lay)
if (iHWCKmodel(lay) == 3 .OR. iHWCKmodel(lay) == 6 .OR. iHWCKmodel(lay) == 7 .OR. &
    iHWCKmodel(lay) == 10 .OR. iHWCKmodel(lay) == 11) then
   sw%vg_params(node)%alpha_2  = paramvg(13, lay)
   sw%vg_params(node)%npar_2   = paramvg(14, lay)
   sw%vg_params(node)%mpar_2   = paramvg(15, lay)
   sw%vg_params(node)%omega_1  = paramvg(16, lay)
   sw%vg_params(node)%omega_2  = paramvg(17, lay)
end if
if (iHWCKmodel(lay) == 5 .OR. iHWCKmodel(lay) == 7) then
   sw%vg_params(node)%h0 = paramvg(18, lay)
end if
if (iHWCKmodel(lay) == 8 .OR. iHWCKmodel(lay) == 9 .OR. &
    iHWCKmodel(lay) == 10 .OR. iHWCKmodel(lay) == 11) then
   sw%vg_params(node)%h0      = paramvg(18, lay)
   sw%vg_params(node)%ha      = paramvg(19, lay)
   sw%vg_params(node)%apar    = paramvg(20, lay)
   sw%vg_params(node)%omega_k = paramvg(21, lay)
end if
```

Do the analogous thing for the tabulated branch (lines 845-848).

- [ ] **Step 4: Clean rebuild + verify all gates.**

```
rm -rf builddir
pixi run build-linux
pixi run test-pfunit
pixi run -e test python tests/regression/test_output_regression.py hupselbrook surfacewater salinitystress grassgrowth
```

Regression must remain 4/4 byte-for-byte. If a case fails, a reader was missed — find and migrate.

- [ ] **Step 5: Commit.**

```bash
git add src/state/soilwater_state.f90 src/soil/soilhydraulics.f90
git commit -m "$(cat <<'EOF'
retire(gr-utils): delete state%soilwater%cofgen

All readers have migrated to state%soilwater%vg_params. The legacy
21-row magic-index matrix is retired. SoilHydraulics(1) populates
vg_params directly from paramvg + relsatthr + ksatthr — no
transitional duplication remains.

After this commit: state%soilwater is one cofgen-shaped record
smaller. No bind_*_target. No magic indices. Utility functions
are pure: signatures self-document dependencies.

Regression: 4/4 byte-for-byte.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 16: Final verification + arc-complete marker

**Files:** No code changes.

- [ ] **Step 1: check-full.** Run: `pixi run check-full` — 5/5 byte-for-byte.

- [ ] **Step 2: BMI + cffi-demo suites.** Run: `pixi run -e test test-bmi && pixi run -e test test-cffi-demo` — both pass.

- [ ] **Step 3: Retirement greps.**

```bash
echo "=== module-level state pointers in utils ==="
grep -rn ", pointer ::" src/utils/soilhydraulicsutils.f90 src/utils/surfacewaterutils.f90 \
                         src/soil/WC_K_models_04_11.f90
echo "=== bind_*_target ==="
grep -rn "bind_state_targets\\|bind_tc_target\\|bind_cofgen_target" src/
echo "=== use variables in utils + WC_K_models ==="
grep -n "use variables" src/utils/soilhydraulicsutils.f90 src/utils/surfacewaterutils.f90 \
                         src/soil/WC_K_models_04_11.f90
echo "=== state%soilwater%cofgen ==="
grep -rn "state%soilwater%cofgen\\|sw%cofgen" src/
echo "=== state%soilwater%vg_params readers ==="
grep -rn "vg_params" src/ --include="*.f90" | wc -l
```

Expected:
- Module-level pointers: zero hits (or only deletion-tombstone comments).
- `bind_*_target`: zero live-code hits.
- `use variables` in those 3 files: zero hits.
- `state%soilwater%cofgen`: zero live-code hits.
- `vg_params` readers: many (40+ across the codebase).

- [ ] **Step 4: Tag the arc complete.**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-utils): Arc 1 complete — utils globals retirement + cofgen modernization

Summary of the arc:

- New module src/state/hydraulic_params_mod.f90 defines
  vanGenuchten_params_t (20 named fields).
- state%soilwater%cofgen(21, numnod) retired in favor of
  state%soilwater%vg_params(:) (array of derived types).
- 9 additional soilwater + 6 additional surfacewater bare globals
  migrated to state subrecords.
- 5 soilhydraulicsutils functions + functionvalue_04_11 take
  explicit, self-documenting signatures: (head, vg, model, dt,
  use_ksatexm, ...) instead of relying on module-level pointers.
- bind_state_targets, bind_tc_target, bind_cofgen_target — all
  deleted. The 3 bind_*_target calls + 2 imports in swap_mod.f90
  retired.
- `use variables` dropped from soilhydraulicsutils, surfacewaterutils,
  WC_K_models_04_11 (3 files).
- 50+ call sites across crop/, boundary/, atmosphere/, drainage/,
  soil/ updated to pass the new typed args.

check-full: 5/5 byte-for-byte parity. BMI suite + cffi-demo:
passing. swap_mod.f90's strangler-fig leftover for soilwater
binding is gone.

Deferred to follow-on arcs (per globals-retirement roadmap):
- sptab(:,:,:) modernization (similar magic-row pattern;
  Arc 5 / soil cluster)
- Bare-global writes of cofgen in config_to_variables — none
  found; cofgen was populated in SoilHydraulics(1), not in
  config_to_variables. swsophy/numtab/etc. bare-globals continue
  to live in variables.f90 until Arc 9 retires the adapter.
- Other clusters: boundary, heat, atmosphere, soil, drainage, io,
  crop (Arcs 2-8).

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Plan self-review

**Spec coverage:**

| Spec requirement | Task |
|---|---|
| New vanGenuchten_params_t | Task 2 |
| state%soilwater%vg_params(:) field | Task 2 |
| Migrate cofgen writes → vg_params | Task 3 (transitional dual-write), Task 15 (cofgen retirement) |
| Add 9 soilwater + 6 surfacewater fields | Task 4 |
| watcon signature change | Task 5 |
| hconduc + dhconduc signature change | Task 6 |
| moiscap signature change | Task 7 |
| prhead signature change | Task 8 |
| functionvalue_04_11 signature change | Task 9 |
| Inline cofgen readers in soilhydraulics | Task 10 |
| cofgenNew in soilgrid | Task 8 (via prhead caller) + Task 11 (cleanup) |
| surfacewaterutils migration | Task 12 |
| Drop use variables + module pointers + bind setups | Task 13 |
| Drop bind calls from swap_mod | Task 14 |
| Delete state%soilwater%cofgen | Task 15 |

**Placeholder scan:** "Implementer chooses during the work" appears in Task 5 — but the decision is locked (the 5-arg form with `node` and `soilwater` is specified). Acceptable.

"Confirm exact rank/kind by reading variables.f90" appears in Task 4 — this is an explicit verification step, not a placeholder. The implementer literally runs the cited grep command.

The `sptab` 3rd-dimension upper bound is described as "verify the legacy declaration" — same; it's an explicit verification step.

**Type consistency check:** `vanGenuchten_params_t` is used consistently across Tasks 2, 5, 6, 7, 8, 9. Field names (`thetar`, `thetas`, `alpha`, `ksatexm`, etc.) match between the type definition and the function bodies.

`state%soilwater%vg_params` access path is consistent — singular field name `vg_params(:)` in all reads.

`state%soilwater%iHWCKmodel(state%soilwater%layer(i))` is the canonical pattern for fetching the model index at a node — used consistently in caller transforms.

**Ordering:**
- Task 2 (type + field) before Task 3 (population).
- Tasks 3 + 4 (populate state) before Tasks 5-9 (readers).
- Task 9 (functionvalue_04_11) after Task 5 (watcon, which calls it).
- Task 13 (drop bind setups) after Tasks 5-9 (all readers migrated).
- Task 14 (drop bind calls in swap_mod) after Task 13.
- Task 15 (delete cofgen) after all reads migrated (Tasks 5-12).

Dependency chain is correct.

---

Plan complete and saved to `docs/superpowers/plans/2026-05-13-globals-utils.md`. Sized at 16 tasks. Ready to execute.
