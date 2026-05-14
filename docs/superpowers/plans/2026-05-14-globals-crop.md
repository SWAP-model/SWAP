# GR-CROP Implementation Plan — Crop + Nutrients + Rain timing + Cross-arc closures

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Migrate the entire crop subsystem (12 files), the nutrients subsystem (management_soil + wofost_soil_*), and the rain timing arrays (in meteodt) off bare `use variables`. Expand `state%crop` with sub-records per crop type (common/fixed/wofost/grass). Introduce `state%nutrients`. Close inherited cross-arc deferrals (GET_MAX_RESP_FACTOR state threading, irrigation write-side, rain timing, atmosphere multi-consumer globals retirement). Retire crop + atmosphere multi-consumer globals from `variables.f90`.

**Architecture:** Four sub-phases. Phase A is additive schema (4 crop sub-records + nutrients_state_t + atmosphere rain fields + config%meteo raintab). Phase A.5 adds runtime dual-write coverage at every legacy mutation site (lesson from GR-ATM: zero-resets count as writes). Phase B migrates 12 crop files + meteodt's rain timing site. Phase C sweeps the codebase, audits, and retires legacy globals from `variables.f90`.

**Tech Stack:** Fortran 2018, gfortran, Meson, pFUnit 4.15, pixi, Python regression tests.

**Spec:** `docs/superpowers/specs/2026-05-14-globals-crop-design.md`

**Builds on:** GR-ATM (Arc 4, completed 2026-05-14 — commit `57b00b3`); GR-BH (Arc 2+3, completed 2026-05-13 — commit `a65bdf3`); GR-UTILS (Arc 1, completed 2026-05-12 — commit `26123ad`).

**Verification gate per task (VG)** — per memory `feedback_state_schema_clean_rebuild.md` + `feedback_per_task_regression_gate.md`:

```bash
rm -rf builddir && pixi run build-linux \
  && pixi run test-pfunit \
  && pixi run -e test python tests/regression/test_output_regression.py \
       hupselbrook surfacewater salinitystress grassgrowth
```

Expected: clean build, pFUnit all-pass, regression **4/4 byte-for-byte**. Any deviation = task NOT complete; fix or escalate.

**Phase close gate:** `pixi run check-full` 5/5 byte-for-byte.

**Controller quality bar (per user directive):** subagents do NOT defer due to missing infrastructure. If a deferral surfaces (state arg not threaded, sibling reader still on legacy, schema gap), the controller intervenes with an inline fix — thread state through the caller chain, extend Phase A.5 coverage, expand schema. NO cascading deferrals.

**Subagent dispatch convention:** every implementer prompt MUST include this VG verbatim. `rm -rf builddir` is mandatory because Meson's incremental build does not propagate `.mod` deps across the `swap_modern`/`swap_legacy` static-library boundary when state schema changes.

---

## File Structure

**New modules:**
- `src/state/crop_common_state.f90` — `crop_common_state_t` (shared across crop types)
- `src/state/crop_fixed_state.f90` — `crop_fixed_state_t`
- `src/state/crop_wofost_state.f90` — `crop_wofost_state_t`
- `src/state/crop_grass_state.f90` — `crop_grass_state_t`
- `src/state/nutrients_state.f90` — `nutrients_state_t`

**Extended modules:**
- `src/state/crop_state.f90` — adds 4 sub-record fields; existing 12 top-level fields preserved
- `src/state/swap_state.f90` — adds `type(nutrients_state_t) :: nutrients`
- `src/state/atmosphere_state.f90` — +4 rain timing fields
- `src/config/meteorology_config.f90` — +`raintab(60)`
- `src/io/toml/config_to_variables.f90` — adapter populator for `raintab` + any new config-sourced crop fields

**Migrated readers (full or near-full `use variables` retirement):**
- `src/crop/cropfixed_init.f90` (1 narrow site)
- `src/crop/cropgrass_init.f90` (2 narrow sites)
- `src/crop/cropwofost_init.f90` (2 narrow sites)
- `src/crop/wofost_soil_parameters.f90` (1 narrow site)
- `src/crop/cropgrowth.f90` (3 sites — bare lines 24, 549, 787)
- `src/crop/oxygenstress.f90` (2 sites incl. GET_MAX_RESP_FACTOR)
- `src/crop/rootextraction.f90` (3 bare sites)
- `src/crop/irrigation.f90` (2 sites — close gird/isua write-side)
- `src/crop/tillage.f90` (1 narrow site)
- `src/crop/management_soil.f90` (blanket use Variables)
- `src/crop/wofost_soil_*.f90` cluster (8 small files)
- `src/atmosphere/meteodt.f90` rain timing site (line 143)

**Codebase sweep readers (Phase C):**
- `src/core/`: `swap_mod.f90`, `initialize.f90`, `timecontrol_mod.f90`, `swap_bmi_mod.f90`, `swap_capi_mod.f90`
- `src/io/`: `swap_csv_output.f90`, `swapoutput.f90`, `readmeteo.f90`
- `src/atmosphere/`: `meteoday.f90`, `meteodt.f90`, `et.f90`, `interception.f90`
- `src/boundary/`, `src/heat/`, `src/soil/`, `src/drainage/`, `src/solute/` clusters

**Populator + adapter shrinks (Phase C):**
- `src/io/toml/config_to_variables.f90` — adds Phase A dual-writes; legacy global writes deleted in Phase C
- `src/io/readmeteo.f90` — rain timing populator routes to state in Phase C
- `src/core/initialize.f90` — loses zero-fills for deleted globals
- `src/core/variables.f90` — declarations of crop + atmosphere multi-consumer globals deleted

**Build:**
- `meson.build` — add new state modules **before** `src/state/crop_state.f90` (which `use`s them) and `src/state/swap_state.f90`. Ordering:
  1. `crop_common_state.f90`, `crop_fixed_state.f90`, `crop_wofost_state.f90`, `crop_grass_state.f90`
  2. `crop_state.f90` (existing, modified)
  3. `nutrients_state.f90`
  4. `swap_state.f90` (existing, modified)
- `tests/unit/meson.build` — mirror in `pfunit_extra_sources`

---

## Phase A — Schema + Dual-Write Seeding (Tasks 1–18)

Additive only. State mirrors populated alongside legacy globals; no reader migration.

### Task 1: Pre-flight baseline

**Files:** None.

- [ ] **Step 1: Clean rebuild.** `rm -rf builddir && pixi run build-linux` — must succeed.
- [ ] **Step 2: pFUnit baseline.** `pixi run test-pfunit` — record pass count.
- [ ] **Step 3: check-full baseline.** `pixi run check-full` — 5/5 pass.
- [ ] **Step 4: Empty marker commit.**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-crop): pre-flight baseline — crop+nutrients+rain timing arc begins

check-full 5/5. Baseline locked before crop_state sub-record expansion,
nutrients_state_t introduction, rain timing schema, 12-file reader
migration, and inherited cross-arc closures.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Create `crop_common_state_mod` (empty type body)

**Files:**
- Create: `src/state/crop_common_state.f90`
- Modify: `src/state/crop_state.f90` (add `use crop_common_state_mod, only: crop_common_state_t`; add `type(crop_common_state_t) :: common` field after the existing 12 fields)
- Modify: `meson.build` (add `'src/state/crop_common_state.f90',` **before** `'src/state/crop_state.f90'`)
- Modify: `tests/unit/meson.build` (mirror in `pfunit_extra_sources`)

- [ ] **Step 1: Create the new module with empty type body.**

```fortran
!> @file crop_common_state.f90
!! SS-GR-CROP: typed crop runtime state — shared across all crop types
!! (fixed/wofost/grass). Fields populated in Task A6 via audit of
!! cropgrowth + init files.
module crop_common_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: crop_common_state_t

   type :: crop_common_state_t
      ! Fields added in Task A6 — placeholder for build wiring.
      integer :: placeholder_ = 0
   contains
      procedure :: init => crop_common_state_init
   end type crop_common_state_t

contains

   subroutine crop_common_state_init(self)
      class(crop_common_state_t), intent(inout) :: self
      ! No-op for now; fields seeded in Task A14.
   end subroutine

end module crop_common_state_mod
```

- [ ] **Step 2: Wire into crop_state_mod.**

Edit `src/state/crop_state.f90`. Add to imports near the top (before `implicit none`):
```fortran
   use crop_common_state_mod, only: crop_common_state_t
```

Inside the type body (after the existing 12 atmospheric-exchange fields), add:
```fortran
      ! [SS-GR-CROP A2] sub-record for shared crop runtime fields
      type(crop_common_state_t) :: common
```

- [ ] **Step 3: Add to meson.build.**

Find the legacy sources list. Add `'src/state/crop_common_state.f90',` **immediately before** `'src/state/crop_state.f90'`.

- [ ] **Step 4: Add to tests/unit/meson.build.**

Find `pfunit_extra_sources`. Add `'../../src/state/crop_common_state.f90',` before `'../../src/state/crop_state.f90'`.

- [ ] **Step 5: Run VG.**

- [ ] **Step 6: Commit.**

```bash
git add src/state/crop_common_state.f90 src/state/crop_state.f90 meson.build tests/unit/meson.build
git commit -m "$(cat <<'EOF'
schema(gr-crop): add crop_common_state_mod (empty body)

Sub-record for crop fields shared across all crop types. Placeholder
field will be replaced with audited field list in Task A6.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: Create `crop_fixed_state_mod` (empty type body)

**Files:** Same pattern as Task 2 — create `src/state/crop_fixed_state.f90`, wire into `crop_state_mod` (add `type(crop_fixed_state_t) :: fixed` after `common`), add to meson and tests/unit/meson.

- [ ] **Step 1: Create the new module.**

```fortran
!> @file crop_fixed_state.f90
!! SS-GR-CROP: typed runtime state for fixed-crop simulations.
!! Fields populated in Task A9.
module crop_fixed_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: crop_fixed_state_t

   type :: crop_fixed_state_t
      integer :: placeholder_ = 0
   contains
      procedure :: init => crop_fixed_state_init
   end type crop_fixed_state_t

contains

   subroutine crop_fixed_state_init(self)
      class(crop_fixed_state_t), intent(inout) :: self
   end subroutine

end module crop_fixed_state_mod
```

- [ ] **Step 2: Wire into crop_state_mod.**

Add `use crop_fixed_state_mod, only: crop_fixed_state_t` to imports.
Add `type(crop_fixed_state_t) :: fixed` to the type body after `common`.

- [ ] **Step 3: Add to meson.build** before `crop_state.f90`.
- [ ] **Step 4: Add to tests/unit/meson.build** in same order.
- [ ] **Step 5: Run VG.**
- [ ] **Step 6: Commit.**

```bash
git add src/state/crop_fixed_state.f90 src/state/crop_state.f90 meson.build tests/unit/meson.build
git commit -m "schema(gr-crop): add crop_fixed_state_mod (empty body)

Sub-record for fixed-crop runtime. Fields populated in Task A9.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 4: Create `crop_wofost_state_mod` (empty type body)

**Files:** Same pattern. Create `src/state/crop_wofost_state.f90`, wire `type(crop_wofost_state_t) :: wofost` into `crop_state_t` after `fixed`, add to meson.

- [ ] **Step 1: Create module** (analogous to Task 3 but with type name `crop_wofost_state_t`).
- [ ] **Step 2: Wire into crop_state_mod.**
- [ ] **Step 3-5: Meson + tests + VG.**
- [ ] **Step 6: Commit.**

```bash
git commit -m "schema(gr-crop): add crop_wofost_state_mod (empty body)

Sub-record for WOFOST runtime (biomass pools + flows). Fields populated
in Task A7.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 5: Create `crop_grass_state_mod` (empty type body)

**Files:** Same pattern. Create `src/state/crop_grass_state.f90`, wire `type(crop_grass_state_t) :: grass` into `crop_state_t` after `wofost`.

- [ ] **Step 1: Create module** (analogous, type `crop_grass_state_t`).
- [ ] **Step 2: Wire into crop_state_mod.**
- [ ] **Step 3-5: Meson + tests + VG.**
- [ ] **Step 6: Commit.**

```bash
git commit -m "schema(gr-crop): add crop_grass_state_mod (empty body)

Sub-record for grass runtime. Fields populated in Task A8.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 6: Populate `crop_common_state_t` with shared crop fields

**Files:** `src/state/crop_common_state.f90`

This task adds the actual field list. Implementer audits `cropgrowth.f90` + init files to identify symbols read by ALL crop types.

- [ ] **Step 1: Audit shared symbols.**

```bash
grep -nE "use variables" src/crop/cropgrowth.f90 src/crop/cropfixed_init.f90 src/crop/cropgrass_init.f90 src/crop/cropwofost_init.f90 | head -20
```

Cross-reference each symbol with variables.f90 to determine if it's used by all crop types or only one. Symbols representing development, root depth, crop height, crop factor, calendar flags are typically shared.

- [ ] **Step 2: Replace placeholder with audited field list.**

Edit `src/state/crop_common_state.f90`. Replace `integer :: placeholder_ = 0` with:

```fortran
      ! Development & calendar
      integer      :: daycrop        = 0    !! days since crop start
      real(real64) :: dvs            = 0.0_real64    !! development stage
      real(real64) :: tsum           = 0.0_real64    !! temperature sum
      integer      :: daynrsta       = 0
      integer      :: daynrend       = 0
      integer      :: swcrp          = 0    !! 1=fixed, 2=wofost, 3=grass
      integer      :: icrop          = 0
      logical      :: flcropcalendar = .false.
      logical      :: flcropoutput   = .false.
      logical      :: flcropnut      = .false.
      logical      :: flharvestday   = .false.
      logical      :: swend          = .false.

      ! Root depth
      real(real64) :: rd             = 0.0_real64    !! actual rooting depth (cm)
      real(real64) :: rdpot          = 0.0_real64    !! potential rooting depth
      real(real64) :: rdm            = 0.0_real64    !! max rooting depth
      real(real64) :: rri            = 0.0_real64    !! max daily root depth increase
      real(real64) :: rdi            = 0.0_real64    !! initial rooting depth
      real(real64) :: rdc            = 0.0_real64    !! max crop rooting depth

      ! Crop physiology summary
      real(real64) :: ch             = 0.0_real64    !! crop height
      real(real64) :: cf             = 0.0_real64    !! crop factor
      real(real64) :: laipot         = 0.0_real64    !! potential LAI

      ! Grazing/harvest
      real(real64) :: cuptgraz       = 0.0_real64
      real(real64) :: cuptgrazpot    = 0.0_real64
      real(real64) :: HarLosOrm_tot  = 0.0_real64
```

**Audit gap:** if your audit finds additional shared fields (e.g., `Tbase`, `tsumea`, `tsumam`, etc.), add them. Use the same comment-block grouping.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/state/crop_common_state.f90
git commit -m "$(cat <<'EOF'
schema(gr-crop): populate crop_common_state_t with shared fields

Development (daycrop/dvs/tsum/calendar/flags) + root depth (rd/rdpot/
rdm/rri/rdi/rdc) + crop physiology summary (ch/cf/laipot) + grazing/
harvest (cuptgraz/HarLosOrm_tot). Audit-driven list from cropgrowth +
init files. Unpopulated until Task A14.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 7: Populate `crop_wofost_state_t` with WOFOST runtime fields

**Files:** `src/state/crop_wofost_state.f90`

- [ ] **Step 1: Audit WOFOST-specific symbols.**

```bash
grep -nE "wlv|wst|wrt|wso|tagp|pgass|dwlv|dwst|dwrt|dwso|swbulb|wbl|cwdm|plossdm|lossdm" src/crop/cropwofost_init.f90 src/crop/cropgrowth.f90 src/core/variables.f90 | head -30
```

- [ ] **Step 2: Replace placeholder with WOFOST fields.**

```fortran
      ! Biomass pools (actual + potential)
      real(real64) :: wlv    = 0.0_real64, wlvpot  = 0.0_real64
      real(real64) :: wst    = 0.0_real64, wstpot  = 0.0_real64
      real(real64) :: wrt    = 0.0_real64, wrtpot  = 0.0_real64
      real(real64) :: wso    = 0.0_real64, wsopot  = 0.0_real64
      real(real64) :: twlv   = 0.0_real64, twst    = 0.0_real64
      real(real64) :: tagp   = 0.0_real64, tagppot = 0.0_real64
      real(real64) :: tagpt  = 0.0_real64, tagptpot = 0.0_real64
      real(real64) :: cwdm   = 0.0_real64, cwdmpot = 0.0_real64
      real(real64) :: pgass  = 0.0_real64, pgasspot = 0.0_real64

      ! Death/decay flows
      real(real64) :: dwlv   = 0.0_real64, dwlvpot = 0.0_real64
      real(real64) :: dwst   = 0.0_real64, dwstpot = 0.0_real64
      real(real64) :: dwrt   = 0.0_real64, dwrtpot = 0.0_real64
      real(real64) :: dwso   = 0.0_real64, dwsopot = 0.0_real64
      real(real64) :: dwlvCrop = 0.0_real64, dwlvSoil = 0.0_real64

      ! Losses
      real(real64) :: plossdm = 0.0_real64, lossdm = 0.0_real64

      ! Bulb-crop fields
      logical      :: swbulb  = .false.
      real(real64) :: wbl    = 0.0_real64, wblpot  = 0.0_real64
      real(real64) :: dwbl   = 0.0_real64, dwblpot = 0.0_real64
```

If audit surfaces additional WOFOST runtime fields, add them.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git commit -m "schema(gr-crop): populate crop_wofost_state_t — biomass pools + flows

Audit-driven WOFOST runtime: actual+potential biomass pools (wlv/wst/wrt/
wso + twX + tagp/cwdm/pgass), death/decay flows (dwX), losses (plossdm),
bulb-crop fields. Unpopulated until Task A15.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 8: Populate `crop_grass_state_t` with grass-specific fields

**Files:** `src/state/crop_grass_state.f90`

- [ ] **Step 1: Audit grass-specific symbols.**

```bash
grep -nE "plwt|swardstate|mowing|dateharvest|grass" src/crop/cropgrass_init.f90 src/crop/cropgrowth.f90 src/core/variables.f90 | head -30
```

- [ ] **Step 2: Replace placeholder with grass fields.**

```fortran
      real(real64) :: plwt = 0.0_real64    !! plant weight (live)
      ! Grass-specific runtime — audit-driven additions:
      ! e.g. mowing schedule state, swardstate, grazing flags, etc.
```

Implementer adds the actual fields per audit. If the audit surfaces only `plwt` and the rest of grass state is shared with `crop_common_state_t`, document that.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git commit -m "schema(gr-crop): populate crop_grass_state_t — grass-specific fields

Audit-driven additions: plwt + <implementer documents what audit found>.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 9: Populate `crop_fixed_state_t` with fixed-crop fields

**Files:** `src/state/crop_fixed_state.f90`

- [ ] **Step 1: Audit fixed-crop-specific symbols.**

```bash
grep -nE "cftb|swcf|cfeic|cfeictb|fixed" src/crop/cropfixed_init.f90 src/crop/cropgrowth.f90 | head -20
```

- [ ] **Step 2: Replace placeholder.**

```fortran
      use swap_array_dimensions, only: magrs  ! add to module imports
      ! ...
      real(real64) :: cftb(2*magrs)    = 0.0_real64    !! crop factor or height vs DVS
      ! Other fixed-crop-only runtime — implementer documents audit findings
```

If `swcf/swcfbs/cfbs` are already on top-level `state%crop`, do NOT duplicate them here.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git commit -m "schema(gr-crop): populate crop_fixed_state_t — cftb + fixed-only fields

cftb table + <audit findings>.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 10: Create `nutrients_state_mod` (skeleton)

**Files:**
- Create: `src/state/nutrients_state.f90`
- Modify: `src/state/swap_state.f90` (add `use nutrients_state_mod, only: nutrients_state_t`; add `type(nutrients_state_t) :: nutrients` field after `type(crop_state_t) :: crop`)
- Modify: `meson.build` (add `'src/state/nutrients_state.f90',` between `crop_state.f90` and `swap_state.f90`)
- Modify: `tests/unit/meson.build` (mirror)

- [ ] **Step 1: Create the new module.**

```fortran
!> @file nutrients_state.f90
!! SS-GR-CROP: typed nutrients runtime state. Hosts management_soil +
!! wofost_soil_* nutrient pool and flow data. Fields populated in
!! Task A11 via audit.
module nutrients_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: nutrients_state_t

   type :: nutrients_state_t
      integer :: placeholder_ = 0
   contains
      procedure :: init => nutrients_state_init
   end type nutrients_state_t

contains

   subroutine nutrients_state_init(self, nlay)
      class(nutrients_state_t), intent(inout) :: self
      integer, intent(in) :: nlay
      ! Allocation/seeding deferred to Task A11.
   end subroutine

end module nutrients_state_mod
```

- [ ] **Step 2: Wire into swap_state_mod.**

Edit `src/state/swap_state.f90`. Add `use nutrients_state_mod, only: nutrients_state_t` to module imports. Add `type(nutrients_state_t) :: nutrients` after `type(crop_state_t) :: crop`.

- [ ] **Step 3-5: Meson + tests + VG.**

- [ ] **Step 6: Commit.**

```bash
git add src/state/nutrients_state.f90 src/state/swap_state.f90 meson.build tests/unit/meson.build
git commit -m "schema(gr-crop): add nutrients_state_mod (skeleton)

NEW nutrients_state_t for management_soil + wofost_soil_* nutrient
state. Wired into swap_state_t as state%nutrients. Fields populated
in Task A11.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 11: Populate `nutrients_state_t` from management_soil + wofost_soil_* audit

**Files:** `src/state/nutrients_state.f90`

- [ ] **Step 1: Audit nutrient symbols.**

```bash
grep -nE "use variables|use Variables" src/crop/management_soil.f90 src/crop/wofost_soil_*.f90 | head -10
grep -nE "nh4|no3|fom|bio|hum|cnflv|cnfst|nuptot|nflux" src/crop/management_soil.f90 src/crop/wofost_soil_*.f90 src/core/variables.f90 | head -30
```

Document every symbol read by management_soil + wofost_soil_*. Classify per-symbol:
- Per-layer nutrient pools (fom, bio, hum, nh4, no3, etc.) → allocatable arrays sized `nlay` or `numnod`
- Crop N concentration scalars (cnflv, cnfst, etc.) → real scalars
- Flux scratchpads → real scalars

- [ ] **Step 2: Replace placeholder with field list.**

Example skeleton (implementer expands per audit):

```fortran
      ! Soil nitrogen pools (per layer, allocatable)
      real(real64), allocatable :: fom(:)        !! fresh organic matter
      real(real64), allocatable :: bio(:)        !! microbial biomass
      real(real64), allocatable :: hum(:)        !! humus
      real(real64), allocatable :: nh4(:)        !! ammonium
      real(real64), allocatable :: no3(:)        !! nitrate

      ! Crop N demand/uptake (WOFOST)
      real(real64) :: cnflv  = 0.0_real64
      real(real64) :: cnfst  = 0.0_real64
      real(real64) :: nuptot = 0.0_real64

      ! Flux scratchpads (per-timestep)
      real(real64) :: nflux  = 0.0_real64

      ! Additional fields per audit
```

Update `nutrients_state_init` to allocate the per-layer arrays:

```fortran
   subroutine nutrients_state_init(self, nlay)
      class(nutrients_state_t), intent(inout) :: self
      integer, intent(in) :: nlay
      if (.not. allocated(self%fom)) allocate(self%fom(nlay))
      if (.not. allocated(self%bio)) allocate(self%bio(nlay))
      if (.not. allocated(self%hum)) allocate(self%hum(nlay))
      if (.not. allocated(self%nh4)) allocate(self%nh4(nlay))
      if (.not. allocated(self%no3)) allocate(self%no3(nlay))
      self%fom = 0.0_real64
      self%bio = 0.0_real64
      self%hum = 0.0_real64
      self%nh4 = 0.0_real64
      self%no3 = 0.0_real64
   end subroutine
```

- [ ] **Step 3: Add call to `state%nutrients%init(nlay)` in `swap_mod.f90`.**

Place after `soilwater_init` (which uses `nlay`). Find the relevant section via `grep -n "soilwater_init\|nlay\s*=" src/core/swap_mod.f90 | head`.

```fortran
   call state%nutrients%init(nlay)
```

- [ ] **Step 4: Run VG.**

- [ ] **Step 5: Commit.**

```bash
git add src/state/nutrients_state.f90 src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
schema(gr-crop): populate nutrients_state_t from management_soil audit

Soil nitrogen pools (fom/bio/hum/nh4/no3 per-layer allocatable) + crop
N (cnflv/cnfst/nuptot) + flux scratchpad (nflux). Field list audit-driven
from management_soil + wofost_soil_*. nutrients%init(nlay) called from
swap_mod after soilwater_init. Unpopulated by readers until Task A16.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 12: Extend `atmosphere_state_t` with rain timing (4 fields)

**Files:** `src/state/atmosphere_state.f90`

- [ ] **Step 1: Verify `mrain` is in scope.**

```bash
grep -n "swap_array_dimensions\|mrain" src/state/atmosphere_state.f90 | head
```

If not present, `use swap_array_dimensions, only: mrain` (or extend an existing `use swap_array_dimensions, only: ...` line).

- [ ] **Step 2: Add fields.**

In the atmosphere_state_t type body, place near other meteo input arrays:

```fortran
      ! [SS-GR-CROP A12] rain timing — per-year reload runtime state
      integer      :: nmrain                       = 0
      real(real64) :: rainamount(mrain)            = 0.0_real64
      real(real64) :: rainfluxarray(mrain)         = 0.0_real64
      real(real64) :: raintimearray(mrain)         = 0.0_real64
```

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/state/atmosphere_state.f90
git commit -m "schema(gr-crop): extend atmosphere_state_t with rain timing

nmrain + rainamount(mrain) + rainfluxarray(mrain) + raintimearray(mrain).
Per-year-reload runtime state. Closes GR-ATM/Arc 7 deferral. Unpopulated
until Task A17.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 13: Extend `config%meteo` with `raintab` + adapter populator

**Files:**
- Modify: `src/config/meteorology_config.f90` (add `raintab(60)` field)
- Modify: `src/io/toml/config_to_variables.f90` (route legacy `raintab` to read from `config%meteo%raintab`)

- [ ] **Step 1: Add field to meteo config.**

Edit `src/config/meteorology_config.f90`. Near other rain-related fields, add:

```fortran
      real(real64) :: raintab(60) = 0.0_real64    !! rain intensity (cm/d) vs time (T) — swrain==1 input
```

- [ ] **Step 2: Verify legacy default + populator path.**

```bash
grep -nE "raintab\s*=" src/io/toml/config_to_variables.f90 src/core/initialize.f90 | head
```

If a legacy `raintab = ...` line in config_to_variables exists, change to `raintab = config%meteo%raintab`. Else add the line.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/config/meteorology_config.f90 src/io/toml/config_to_variables.f90
git commit -m "schema(gr-crop): extend config%meteo with raintab + adapter populator

raintab(60) field added to meteorology_config_t. Legacy adapter routes
from config%meteo%raintab. Closes GR-ATM/Arc 7 deferral for the static-
input part of rain timing.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 14: Dual-write `state%crop%common` in swap_mod

**Files:** `src/core/swap_mod.f90`

- [ ] **Step 1: Locate placement.**

```bash
grep -n "atmosphere_init\|state%atmosphere%init\|state%crop%lai\s*=\|SS-GR-ATM" src/core/swap_mod.f90 | head -10
```

Find where state%crop top-level fields are seeded (GR-ATM A12 block). Add the common-block dual-write immediately after.

- [ ] **Step 2: Extend `use variables, only:` list.**

Add the legacy global names matching the Task A6 audit. Example:

```fortran
   use variables, only: ..., &
                        daycrop, dvs, tsum, daynrsta, daynrend, swcrp, icrop, &
                        flcropcalendar, flcropoutput, flcropnut, flharvestday, swend, &
                        rd, rdpot, rdm, rri, rdi, rdc, ch, cf, laipot, &
                        cuptgraz, cuptgrazpot, HarLosOrm_tot
```

(Add only names that are bare globals in variables.f90; verify each via `grep -nE "^[[:space:]]*(integer|real\(8\)|logical)\s+<sym>\b" src/core/variables.f90`.)

- [ ] **Step 3: Add dual-write block.**

```fortran
   ! [SS-GR-CROP A14] dual-write crop_common
   state%crop%common%daycrop        = daycrop
   state%crop%common%dvs            = dvs
   state%crop%common%tsum           = tsum
   state%crop%common%daynrsta       = daynrsta
   state%crop%common%daynrend       = daynrend
   state%crop%common%swcrp          = swcrp
   state%crop%common%icrop          = icrop
   state%crop%common%flcropcalendar = flcropcalendar
   state%crop%common%flcropoutput   = flcropoutput
   state%crop%common%flcropnut      = flcropnut
   state%crop%common%flharvestday   = flharvestday
   state%crop%common%swend          = swend
   state%crop%common%rd             = rd
   state%crop%common%rdpot          = rdpot
   state%crop%common%rdm            = rdm
   state%crop%common%rri            = rri
   state%crop%common%rdi            = rdi
   state%crop%common%rdc            = rdc
   state%crop%common%ch             = ch
   state%crop%common%cf             = cf
   state%crop%common%laipot         = laipot
   state%crop%common%cuptgraz       = cuptgraz
   state%crop%common%cuptgrazpot    = cuptgrazpot
   state%crop%common%HarLosOrm_tot  = HarLosOrm_tot
```

- [ ] **Step 4: Run VG.**

- [ ] **Step 5: Commit.**

```bash
git add src/core/swap_mod.f90
git commit -m "schema(gr-crop): dual-write crop_common in swap_mod

State sub-record state%crop%common seeded from legacy globals at swap
init. No reader migrated yet.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 15: Dual-write `state%crop%wofost`, `%fixed`, `%grass` in swap_mod

**Files:** `src/core/swap_mod.f90`

- [ ] **Step 1: Extend `use variables, only:` list.**

Add WOFOST + fixed + grass-specific legacy globals identified in Tasks A7, A8, A9.

- [ ] **Step 2: Add dual-write block.**

```fortran
   ! [SS-GR-CROP A15] dual-write crop_wofost
   state%crop%wofost%wlv     = wlv
   state%crop%wofost%wlvpot  = wlvpot
   state%crop%wofost%wst     = wst
   state%crop%wofost%wstpot  = wstpot
   state%crop%wofost%wrt     = wrt
   state%crop%wofost%wrtpot  = wrtpot
   state%crop%wofost%wso     = wso
   state%crop%wofost%wsopot  = wsopot
   state%crop%wofost%twlv    = twlv
   state%crop%wofost%twst    = twst
   state%crop%wofost%tagp    = tagp
   state%crop%wofost%tagppot = tagppot
   state%crop%wofost%tagpt   = tagpt
   state%crop%wofost%tagptpot = tagptpot
   state%crop%wofost%cwdm    = cwdm
   state%crop%wofost%cwdmpot = cwdmpot
   state%crop%wofost%pgass   = pgass
   state%crop%wofost%pgasspot = pgasspot
   state%crop%wofost%dwlv    = dwlv
   state%crop%wofost%dwlvpot = dwlvpot
   state%crop%wofost%dwst    = dwst
   state%crop%wofost%dwstpot = dwstpot
   state%crop%wofost%dwrt    = dwrt
   state%crop%wofost%dwrtpot = dwrtpot
   state%crop%wofost%dwso    = dwso
   state%crop%wofost%dwsopot = dwsopot
   state%crop%wofost%dwlvCrop = dwlvCrop
   state%crop%wofost%dwlvSoil = dwlvSoil
   state%crop%wofost%plossdm = plossdm
   state%crop%wofost%lossdm  = lossdm
   state%crop%wofost%swbulb  = swbulb
   state%crop%wofost%wbl     = wbl
   state%crop%wofost%wblpot  = wblpot
   state%crop%wofost%dwbl    = dwbl
   state%crop%wofost%dwblpot = dwblpot

   ! [SS-GR-CROP A15] dual-write crop_grass
   state%crop%grass%plwt = plwt
   ! Additional grass dual-writes per Task A8 audit

   ! [SS-GR-CROP A15] dual-write crop_fixed
   state%crop%fixed%cftb = cftb
   ! Additional fixed dual-writes per Task A9 audit
```

**Note:** seed all sub-records defensively. The fields are zero by default for inactive crop types — no harm in dual-writing.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/core/swap_mod.f90
git commit -m "schema(gr-crop): dual-write crop_wofost/fixed/grass in swap_mod

All three crop-type sub-records seeded from legacy globals at swap
init. Inactive crop types receive zero values (default).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 16: Dual-write nutrients in swap_mod

**Files:** `src/core/swap_mod.f90`

- [ ] **Step 1: Extend `use variables, only:` list** for nutrient symbols (per Task A11 audit).

- [ ] **Step 2: Add dual-write block** after `state%nutrients%init(nlay)` call:

```fortran
   ! [SS-GR-CROP A16] dual-write nutrients
   state%nutrients%fom(:) = fom(1:size(state%nutrients%fom))   ! slice if legacy is fixed-size
   state%nutrients%bio(:) = bio(1:size(state%nutrients%bio))
   state%nutrients%hum(:) = hum(1:size(state%nutrients%hum))
   state%nutrients%nh4(:) = nh4(1:size(state%nutrients%nh4))
   state%nutrients%no3(:) = no3(1:size(state%nutrients%no3))
   state%nutrients%cnflv  = cnflv
   state%nutrients%cnfst  = cnfst
   state%nutrients%nuptot = nuptot
   state%nutrients%nflux  = nflux
```

Verify legacy declarations in `src/core/variables.f90` to determine sizes and whether slicing is needed. Adjust accordingly.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git commit -m "schema(gr-crop): dual-write nutrients in swap_mod

state%nutrients seeded from legacy globals. Per-layer arrays slice-copied
from legacy maho-sized declarations.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 17: Dual-write rain timing in swap_mod

**Files:** `src/core/swap_mod.f90`

- [ ] **Step 1: Extend `use variables, only:`** with `nmrain, rainamount, rainfluxarray, raintimearray`.

- [ ] **Step 2: Add dual-write block** after atmosphere init / Phase A.10 block:

```fortran
   ! [SS-GR-CROP A17] dual-write rain timing
   state%atmosphere%nmrain        = nmrain
   state%atmosphere%rainamount    = rainamount
   state%atmosphere%rainfluxarray = rainfluxarray
   state%atmosphere%raintimearray = raintimearray
```

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git commit -m "schema(gr-crop): dual-write rain timing in swap_mod

state%atmosphere%(nmrain|rainamount|rainfluxarray|raintimearray) seeded
from legacy globals at swap init. Per-year reload routing happens in
Phase A.5.5 + B15.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 18: Phase A close marker

**Files:** None.

- [ ] **Step 1: Clean rebuild + check-full.** `rm -rf builddir && pixi run build-linux && pixi run check-full` — 5/5 byte-for-byte.

- [ ] **Step 2: Sanity greps.**

```bash
grep -n "state%crop%common%daycrop\s*=\|state%crop%wofost%wlv\s*=\|state%nutrients%fom\s*=\|state%atmosphere%nmrain\s*=" src/core/swap_mod.f90
```

Expected: ≥4 hits (the Task 14, 15, 16, 17 dual-writes).

- [ ] **Step 3: Phase A close marker.**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-crop): Phase A complete — schema + dual-write seeding landed

state%crop sub-records (common/fixed/wofost/grass), state%nutrients,
state%atmosphere rain timing, config%meteo raintab seeded alongside
legacy globals. check-full 5/5 byte-for-byte. No readers migrated yet.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Phase A.5 — Runtime Dual-Write Coverage (Tasks A5.1–A5.6)

Critical lesson from GR-ATM: runtime mutations of migrated symbols must mirror to state at every write site, including conditional zero-resets. Each task audits ALL writes for the targeted file(s), including ones inside `if/select case` branches.

### Task A5.1: cropgrowth.f90 dual-writes

**Files:** `src/crop/cropgrowth.f90`

- [ ] **Step 1: Audit ALL writes for the migrated symbols.**

For each symbol in `crop_common_state_t`, `crop_wofost_state_t`, `crop_grass_state_t`, `crop_fixed_state_t`, search:

```bash
grep -nE "^[[:space:]]*\b(daycrop|dvs|tsum|rd\b|rdpot|ch\b|cf\b|laipot|wlv\b|wlvpot|wst\b|wstpot|wrt\b|wrtpot|wso\b|wsopot|twlv|twst|tagp\b|tagppot|tagpt|tagptpot|cwdm|cwdmpot|pgass|pgasspot|dwlv|dwlvpot|dwst|dwstpot|dwrt|dwrtpot|dwso|dwsopot|dwlvCrop|dwlvSoil|plossdm|lossdm|swbulb|wbl|wblpot|dwbl|dwblpot|plwt|cftb|cuptgraz|cuptgrazpot|HarLosOrm_tot|flCropEmergence|flharvestday|flcropoutput|swend)\b\s*(\(|=)" src/crop/cropgrowth.f90 | head -60
```

Inspect each hit. If the write is to a bare global (not state%crop%...), add a mirror.

- [ ] **Step 2: Add mirror writes.**

For each legacy write `<sym> = <expr>`, add immediately after:

```fortran
state%crop%<subrecord>%<sym> = <sym>   ! [SS-GR-CROP A5.1] runtime dual-write
```

Where `<subrecord>` is `common` for shared fields, `wofost` for WOFOST-only, `grass` for grass-only, `fixed` for fixed-only.

**Critical:** if a write is inside a conditional block, place the mirror inside the same block. For zero-resets (e.g., `dvs = 0` on crop end), apply the same rule.

- [ ] **Step 3: Verify state arg is in scope.**

cropgrowth.f90 takes `state` in most subroutines (per prior arc work). For subroutines that don't, controller intervenes — thread `state` through the caller chain.

- [ ] **Step 4: Run VG.**

- [ ] **Step 5: Commit.**

```bash
git add src/crop/cropgrowth.f90
git commit -m "schema(gr-crop A5.1): cropgrowth.f90 dual-writes

Mirror all crop runtime writes to state%crop%(common|wofost|fixed|grass).
Includes conditional zero-resets per GR-ATM commit 8fb76cf lesson.
<implementer documents site count and any state-arg threading done>.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task A5.2: crop init files dual-writes

**Files:** `src/crop/cropfixed_init.f90`, `src/crop/cropgrass_init.f90`, `src/crop/cropwofost_init.f90`

- [ ] **Step 1: Per-file audit.**

```bash
for f in src/crop/cropfixed_init.f90 src/crop/cropgrass_init.f90 src/crop/cropwofost_init.f90; do
  echo "=== $f ==="
  grep -nE "^[[:space:]]*\b(daycrop|dvs|tsum|daynrsta|daynrend|swcrp|icrop|flcropcalendar|flcropoutput|flcropnut|flharvestday|swend|rd|rdpot|rdm|rri|rdi|rdc|ch|cf|laipot|cuptgraz|cuptgrazpot|HarLosOrm_tot|wlv|wlvpot|wst|wstpot|wrt|wrtpot|wso|wsopot|cwdm|cwdmpot|pgass|pgasspot|dwlv|dwlvpot|dwst|dwstpot|dwrt|dwrtpot|dwso|dwsopot|dwlvCrop|dwlvSoil|plossdm|lossdm|swbulb|wbl|wblpot|dwbl|dwblpot|plwt|cftb)\b\s*=" "$f" | head -30
done
```

- [ ] **Step 2: Add mirror writes per Task A5.1 pattern.**

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/crop/cropfixed_init.f90 src/crop/cropgrass_init.f90 src/crop/cropwofost_init.f90
git commit -m "schema(gr-crop A5.2): crop init files dual-writes

Mirror init-time crop assignments to state%crop sub-records. Tracks
state with legacy through crop init.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task A5.3: crop runtime cluster dual-writes (irrigation/oxygenstress/rootextraction/tillage)

**Files:** `src/crop/irrigation.f90`, `src/crop/oxygenstress.f90`, `src/crop/rootextraction.f90`, `src/crop/tillage.f90`

- [ ] **Step 1: Per-file audit** (same grep as Task A5.2).

- [ ] **Step 2: Add mirror writes.**

For irrigation.f90 — close the GR-ATM gird/isua write-side cleanup. The Phase A.5 mirrors from GR-ATM are at `state%crop%gird`, `state%atmosphere%isua` (top-level). Keep those + add any new common/sub-record writes.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/crop/irrigation.f90 src/crop/oxygenstress.f90 src/crop/rootextraction.f90 src/crop/tillage.f90
git commit -m "schema(gr-crop A5.3): crop runtime cluster dual-writes

irrigation/oxygenstress/rootextraction/tillage runtime writes mirrored
to state%crop sub-records.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task A5.4: nutrients dual-writes

**Files:** `src/crop/management_soil.f90`, `src/crop/wofost_soil_*.f90`

- [ ] **Step 1: Audit nutrient writes.**

```bash
for f in src/crop/management_soil.f90 src/crop/wofost_soil_*.f90; do
  echo "=== $f ==="
  grep -nE "^[[:space:]]*\b(fom|bio|hum|nh4|no3|cnflv|cnfst|nuptot|nflux)\b\s*(\(|=)" "$f" | head -20
done
```

- [ ] **Step 2: Add mirror writes to `state%nutrients%X`.**

```fortran
fom(i) = ...
state%nutrients%fom(i) = fom(i)   ! [SS-GR-CROP A5.4]
```

State arg threading: if a nutrient subroutine doesn't take state, thread it (controller intervenes).

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/crop/management_soil.f90 src/crop/wofost_soil_*.f90
git commit -m "schema(gr-crop A5.4): nutrients dual-writes

management_soil + wofost_soil_* mirror nutrient pool/flux writes to
state%nutrients. State arg threaded where needed.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task A5.5: rain timing dual-writes in meteodt

**Files:** `src/atmosphere/meteodt.f90`

- [ ] **Step 1: Audit rain timing writes.**

```bash
grep -nE "^[[:space:]]*\b(nmrain|rainamount|rainfluxarray|raintimearray)\b\s*(\(|=)" src/atmosphere/meteodt.f90 | head -20
```

- [ ] **Step 2: Add mirror writes.**

```fortran
nmrain = rainrec + 1
state%atmosphere%nmrain = nmrain   ! [SS-GR-CROP A5.5]

raintimearray(rainrec) = ...
state%atmosphere%raintimearray(rainrec) = raintimearray(rainrec)
```

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/atmosphere/meteodt.f90
git commit -m "schema(gr-crop A5.5): rain timing dual-writes in meteodt

Per-year rain timing recomputes mirror to state%atmosphere%(nmrain|
rainamount|rainfluxarray|raintimearray).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task A5.6: Phase A.5 close marker

**Files:** None.

- [ ] **Step 1: Clean rebuild + check-full** — 5/5 byte-for-byte.

- [ ] **Step 2: Phase A.5 close marker.**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-crop): Phase A.5 complete — runtime dual-writes landed

~120-160 dual-write sites added across crop runtime + init + nutrients
+ rain timing. State mirrors track legacy through simulation, not just
at init. Phase B reader migration unblocked.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Phase B — Reader Cutover (Tasks B1–B17)

### Symbol replacement reference (Phase B)

| Was (bare global) | Now |
|---|---|
| Common: `daycrop/dvs/tsum/rd/rdpot/rdm/rri/rdi/rdc/ch/cf/laipot/cuptgraz/cuptgrazpot/HarLosOrm_tot/daynrsta/daynrend/swcrp/icrop/flcropcalendar/flcropoutput/flcropnut/flharvestday/swend` | `state%crop%common%X` |
| WOFOST: `wlv/wlvpot/wst/wstpot/wrt/wrtpot/wso/wsopot/twlv/twst/tagp/tagppot/tagpt/tagptpot/cwdm/cwdmpot/pgass/pgasspot/dwlv/dwlvpot/dwst/dwstpot/dwrt/dwrtpot/dwso/dwsopot/dwlvCrop/dwlvSoil/plossdm/lossdm/swbulb/wbl/wblpot/dwbl/dwblpot` | `state%crop%wofost%X` |
| Grass: `plwt` + grass-specific | `state%crop%grass%X` |
| Fixed: `cftb` + fixed-specific | `state%crop%fixed%X` |
| Top-level (existing from GR-ATM): `lai/kdif/kdir/cofab/cfbs/swcf/swcfbs/gird/flCropEmergence/et0/ew0/es0` | `state%crop%X` |
| Nutrients: `fom/bio/hum/nh4/no3/cnflv/cnfst/nuptot/nflux` | `state%nutrients%X` |
| Rain timing: `nmrain/rainamount/rainfluxarray/raintimearray` | `state%atmosphere%X` |
| `raintab` | `config%meteo%raintab` |

### Task B1: Migrate `cropfixed_init.f90` (1 narrow site)

**Files:** `src/crop/cropfixed_init.f90`

- [ ] **Step 1: Inspect site at line 28.**

```bash
sed -n '25,40p' src/crop/cropfixed_init.f90
```

Symbols: `idev, tsumea, tsumam, tbase, ...`. Verify the full only-list.

- [ ] **Step 2: Replace narrow import.**

Drop `use variables, only: idev, tsumea, tsumam, tbase, ...`. Add:

```fortran
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
```

- [ ] **Step 3: Add config arg if not present.**

```bash
grep -n "subroutine.*cropfixed_init\|type(swap_state_t)\|type(swap_config_t)" src/crop/cropfixed_init.f90 | head
```

Signature should look like:
```fortran
   subroutine cropfixed_init(..., state, config)
      ...
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config
```

- [ ] **Step 4: Substitute references.**

`idev/tsumea/tsumam/tbase` are config-loaded crop parameters. Verify their paths in `src/config/cropfixed_config.f90`:

```bash
grep -nE "idev|tsumea|tsumam|tbase" src/config/cropfixed_config.f90 | head
```

Migrate per actual paths (e.g., `config%crop%fixed%idev`). For runtime symbols that should be on state%crop%common, use those.

- [ ] **Step 5: Update callers.**

```bash
grep -rn "call cropfixed_init\b\|call CropFixed\b" src/
```

Add `config` arg if needed.

- [ ] **Step 6: Run VG.**

- [ ] **Step 7: Commit.**

```bash
git add src/crop/cropfixed_init.f90 # + callers
git commit -m "refactor(gr-crop B1): cropfixed_init.f90 — drop use variables

Symbols routed to config%crop%fixed (config-loaded) or state%crop%common
(runtime). Signature gains config arg.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task B2: Migrate `cropgrass_init.f90` (2 narrow sites)

**Files:** `src/crop/cropgrass_init.f90`

Two sites: line 34 (broad use-only) + line 417 (`dateharvest`).

- [ ] **Step 1: Inspect both sites.**

```bash
sed -n '30,50p' src/crop/cropgrass_init.f90
sed -n '415,425p' src/crop/cropgrass_init.f90
```

- [ ] **Step 2: Migrate site 1 (line 34).**

Drop `use variables, only:`. Substitute body references per the spec replacement table.

- [ ] **Step 3: Migrate site 2 (line 417 — `dateharvest`).**

`dateharvest` is grass-specific harvest schedule. Verify if it's on `state%crop%grass` (Task A8 audit) or in `config%crop%grass`. Substitute accordingly. If it's a runtime-mutable array, it should be on state%crop%grass; if read-only schedule, config.

- [ ] **Step 4: Update callers + VG + commit.**

```bash
git add src/crop/cropgrass_init.f90 # + callers
git commit -m "refactor(gr-crop B2): cropgrass_init.f90 — drop use variables

Site 1 (line 34) + Site 2 (line 417 dateharvest) migrated to
state%crop%grass / config%crop%grass per audit.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task B3: Migrate `cropwofost_init.f90` (2 narrow sites)

**Files:** `src/crop/cropwofost_init.f90`

Two sites: line 42 + line 490.

- [ ] **Step 1: Inspect.**

```bash
sed -n '38,55p' src/crop/cropwofost_init.f90
sed -n '485,500p' src/crop/cropwofost_init.f90
```

Site 2 (line 490) symbols include `lrnr, lsnr, nlue, rnflv, rnfst, frnx, nmxlv, ...` — these are WOFOST nutrient parameters. Verify if they map to `config%crop%wofost%nutrient` or `state%nutrients`.

- [ ] **Step 2: Migrate both sites.**

WOFOST biomass init reads `wlv/wst/wrt/wso` etc. — those are NOW on state%crop%wofost. Use those refs. Config-loaded crop parameters go through `config%crop%wofost`.

- [ ] **Step 3: Update callers + VG + commit.**

```bash
git add src/crop/cropwofost_init.f90 # + callers
git commit -m "refactor(gr-crop B3): cropwofost_init.f90 — drop use variables

Site 1 (line 42) + Site 2 (line 490 nutrient params) migrated to
state%crop%wofost / config%crop%wofost / state%nutrients per audit.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task B4: Migrate `wofost_soil_parameters.f90` (1 narrow site)

**Files:** `src/crop/wofost_soil_parameters.f90`

Site at line 14: `BDENS` (bulk density).

- [ ] **Step 1: Locate the BDENS usage.**

```bash
grep -nE "BDENS|Bdens|bdens" src/crop/wofost_soil_parameters.f90 src/state/soilwater_state.f90 src/config/soil_config.f90 | head
```

`BDENS` is likely already on state%soilwater (from GR-BH layer flats) — verify. If yes, route there. If not, add to state%soilwater or use config%soil.

- [ ] **Step 2: Replace import and migrate reads.**

- [ ] **Step 3: VG + commit.**

```bash
git add src/crop/wofost_soil_parameters.f90 # + callers if state arg threaded
git commit -m "refactor(gr-crop B4): wofost_soil_parameters.f90 — drop use variables

BDENS → state%soilwater%bdens (or config path per audit).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task B5: Migrate `cropgrowth.f90` site 1 (line 24 bare)

**Files:** `src/crop/cropgrowth.f90`

Line 24 has bare `use variables, dummy_tsoil_cg_ => tsoil`. The rename is to avoid collision; preserve the intent.

- [ ] **Step 1: Audit subroutines reading this site.**

```bash
sed -n '20,40p' src/crop/cropgrowth.f90
```

The site is at module level (top of file). Subroutines inside the module inherit. Identify which subroutines read which symbols.

- [ ] **Step 2: Convert bare to explicit narrow imports.**

Drop bare `use variables`. Add explicit imports:

```fortran
      use swap_state_mod,         only: swap_state_t
      use swap_config_mod,        only: swap_config_t
      use variables,              only: dummy_tsoil_cg_ => tsoil  ! preserve rename if still needed; else drop
```

If `dummy_tsoil_cg_` is unused after migration, drop the rename entirely. If used inside the module (e.g., to avoid name collision with a local `tsoil`), retain the narrow import.

Verify usage:
```bash
grep -n "dummy_tsoil_cg_\|tsoil\b" src/crop/cropgrowth.f90 | head -20
```

- [ ] **Step 3: Migrate body reads per spec replacement table.**

This site covers MANY subroutines. Symbol substitutions per the table — crop fields → state%crop%X; nutrient fields → state%nutrients%X; etc.

**This is the largest single task.** Expect to substitute hundreds of references. Use sed-like or context-aware editing.

- [ ] **Step 4: Update callers.**

```bash
grep -rn "call CropGrowth\b\|call cropgrowth\b" src/
```

If signatures change to add state/config, update each call site.

- [ ] **Step 5: VG.**

- [ ] **Step 6: Commit.**

```bash
git add src/crop/cropgrowth.f90 # + callers
git commit -m "refactor(gr-crop B5): cropgrowth.f90 site 1 — drop bare use variables

Module-level bare use variables retired (modulo dummy_tsoil_cg_ rename
if still needed). Body reads route to state%crop / state%nutrients /
state%atmosphere / config%X per spec replacement table.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task B6: Migrate `cropgrowth.f90` site 2 (line 549 bare)

**Files:** `src/crop/cropgrowth.f90`

- [ ] **Step 1: Locate subroutine.**

```bash
sed -n '540,560p' src/crop/cropgrowth.f90
```

Find the `subroutine` declaration that contains line 549 (search backward).

- [ ] **Step 2: Replace import + migrate body** per spec replacement table.

- [ ] **Step 3: Update callers + VG + commit.**

```bash
git commit -m "refactor(gr-crop B6): cropgrowth.f90 site 2 (line 549) — drop bare use variables

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task B7: Migrate `cropgrowth.f90` site 3 (line 787 bare)

**Files:** `src/crop/cropgrowth.f90`

- [ ] **Step 1: Locate subroutine.**

```bash
sed -n '780,800p' src/crop/cropgrowth.f90
```

- [ ] **Step 2: Replace import + migrate body.**

- [ ] **Step 3: Update callers + VG + commit.**

```bash
git commit -m "refactor(gr-crop B7): cropgrowth.f90 site 3 (line 787) — drop bare use variables

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task B8: Verify cropgrowth.f90 clean

**Files:** None (verification only).

- [ ] **Step 1: Grep.**

```bash
grep -n "^[[:space:]]*use variables" src/crop/cropgrowth.f90
```

Expected: empty, OR only documented narrow deferrals (e.g., `dummy_tsoil_cg_` rename if still load-bearing).

- [ ] **Step 2: Run VG.**

- [ ] **Step 3: No commit** unless deferrals need documentation.

---

### Task B9: Migrate `oxygenstress.f90` (2 sites, close GET_MAX_RESP_FACTOR deferral)

**Files:** `src/crop/oxygenstress.f90`

Sites: line 67 narrow (`c_mroot, f_senes, max_resp_factor, q10_root, q10_microbial, shape_factor_rootr, specific_resp_humus`) + line 98 bare.

**Critical:** Close the GR-ATM deferral — `GET_MAX_RESP_FACTOR` doesn't take state today. Thread state through it (controller does this inline if subagent reports difficulty).

- [ ] **Step 1: Inspect both sites + GET_MAX_RESP_FACTOR.**

```bash
sed -n '60,75p' src/crop/oxygenstress.f90
sed -n '95,110p' src/crop/oxygenstress.f90
grep -n "GET_MAX_RESP_FACTOR\|get_max_resp_factor" src/crop/oxygenstress.f90 | head
```

- [ ] **Step 2: Thread state through GET_MAX_RESP_FACTOR.**

If the function signature is `function GET_MAX_RESP_FACTOR(...)` without state, change to `function GET_MAX_RESP_FACTOR(..., state)` with `type(swap_state_t), intent(in) :: state` arg. Update internal `tav` reads to `state%atmosphere%Tav`. Update callers within oxygenstress.f90.

- [ ] **Step 3: Migrate both sites.**

Replace narrow imports with state/config. Body reads migrate per spec replacement table.

`c_mroot/f_senes/max_resp_factor/q10_root/q10_microbial/shape_factor_rootr/specific_resp_humus` — likely config-loaded oxygen stress params. Verify path in `src/config/crop_config.f90`:

```bash
grep -nE "c_mroot|f_senes|max_resp_factor|q10_root|q10_microbial|shape_factor_rootr|specific_resp_humus" src/config/*.f90 | head
```

Route to config%crop or wherever they live.

- [ ] **Step 4: Update callers.**

```bash
grep -rn "call OxygenStress\b\|call oxygenstress\b\|call ReprodFunction\b\|call GET_MAX_RESP_FACTOR\b" src/
```

- [ ] **Step 5: VG + commit.**

```bash
git add src/crop/oxygenstress.f90 # + callers
git commit -m "refactor(gr-crop B9): oxygenstress.f90 — drop use variables; close GET_MAX_RESP_FACTOR

GET_MAX_RESP_FACTOR signature now takes state arg; tav reads route to
state%atmosphere%Tav. Closes GR-ATM deferral. Config-loaded oxygen stress
params → config%crop.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task B10: Migrate `rootextraction.f90` (3 bare sites)

**Files:** `src/crop/rootextraction.f90`

Sites at lines 41, 342, 701 (all bare).

- [ ] **Step 1: Inspect each site + subroutine.**

```bash
sed -n '35,55p' src/crop/rootextraction.f90
sed -n '335,355p' src/crop/rootextraction.f90
sed -n '695,710p' src/crop/rootextraction.f90
```

- [ ] **Step 2: Replace imports + migrate body per spec table.**

- [ ] **Step 3: Update callers + VG + commit.**

```bash
git commit -m "refactor(gr-crop B10): rootextraction.f90 (3 sites) — drop bare use variables

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task B11: Migrate `irrigation.f90` (2 sites, close gird/isua write-side)

**Files:** `src/crop/irrigation.f90`

Sites: line 40 bare + line 312 narrow (mostly already cleaned by GR-BH C7).

- [ ] **Step 1: Inspect sites.**

```bash
sed -n '35,55p' src/crop/irrigation.f90
sed -n '305,325p' src/crop/irrigation.f90
```

- [ ] **Step 2: Migrate.**

`gird` writes already mirror to `state%crop%gird` from GR-ATM A5.3. After this migration, the legacy global write can be dropped (state-only write).

`isua` writes mirror to `state%atmosphere%isua` from GR-ATM A5.3.

- [ ] **Step 3: Drop legacy mirrors at write sites.**

After migration, the bare `gird = ...` lines become `state%crop%gird = ...` (state-only). Verify there are no remaining bare-global readers of `gird`/`isua` outside irrigation.f90:

```bash
grep -rn "\bgird\b\|\bisua\b" src/ --include="*.f90" | grep -v "state%\|config%\|! \|variables.f90\|use variables\|integer ::\|real(8) ::"
```

If readers exist outside, retain dual-write until those migrate. Otherwise drop.

- [ ] **Step 4: VG + commit.**

```bash
git commit -m "refactor(gr-crop B11): irrigation.f90 — drop use variables; close gird/isua write-side

gird/isua writes route to state%crop/state%atmosphere only. Legacy
mirror dropped after Phase C audit confirms no remaining bare readers.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task B12: Migrate `tillage.f90` (1 narrow site)

**Files:** `src/crop/tillage.f90`

Site at line 10: `swhyst, swsolu, swoxygen, Bdens, ParamVG`.

- [ ] **Step 1: Inspect.**

```bash
sed -n '5,25p' src/crop/tillage.f90
```

These are mostly config-loaded soil parameters. `Bdens` and `ParamVG` may already be on state%soilwater (verify).

- [ ] **Step 2: Replace import + migrate per audit.**

- [ ] **Step 3: VG + commit.**

```bash
git commit -m "refactor(gr-crop B12): tillage.f90 — drop use variables

swhyst/swsolu/swoxygen → config%soil. Bdens/ParamVG → state%soilwater
(or config per audit).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task B13: Migrate `management_soil.f90` (blanket use Variables)

**Files:** `src/crop/management_soil.f90`

- [ ] **Step 1: Audit blanket use Variables.**

```bash
grep -n "^[[:space:]]*use [Vv]ariables" src/crop/management_soil.f90
```

Inspect every symbol read in the file body via:

```bash
grep -nE "\b[A-Za-z_][A-Za-z0-9_]*\b" src/crop/management_soil.f90 | head -50
```

(This is broad; narrow by intersecting with variables.f90 symbol list.)

- [ ] **Step 2: Convert blanket to explicit narrow imports.**

Replace blanket `use Variables` with explicit `use swap_state_mod, only: swap_state_t` + `use swap_config_mod, only: swap_config_t` + (if needed) narrow `use variables, only: <deferred-symbols>` with documentation.

- [ ] **Step 3: Migrate body reads** per spec replacement table.

Nutrient writes/reads use `state%nutrients%X`. Mesh/soil reads use `state%mesh%X` / `state%soilwater%X`. Crop runtime reads use `state%crop%X`.

- [ ] **Step 4: Add `state` arg to signatures.**

Most subroutines in management_soil should take state. Update callers (in cropgrowth + swap_mod).

- [ ] **Step 5: VG + commit.**

```bash
git commit -m "refactor(gr-crop B13): management_soil.f90 — drop blanket use Variables

Nutrient + soil + crop runtime reads routed to state%(nutrients|soilwater|
crop|mesh|atmosphere)%X. State threaded through all subroutines.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task B14: Migrate `wofost_soil_*.f90` cluster

**Files:** `src/crop/wofost_soil_amendments.f90`, `wofost_soil_balancecheck.f90`, `wofost_soil_cropresidues.f90`, `wofost_soil_declarations.f90`, `wofost_soil_interface.f90`, `wofost_soil_orgmatn.f90`, `wofost_soil_rateconstants.f90`, `wofost_soil_watern.f90`

Most are small; may or may not have `use variables` (audit per-file).

- [ ] **Step 1: Per-file audit.**

```bash
for f in src/crop/wofost_soil_*.f90; do
  echo "=== $f ==="
  grep -n "^[[:space:]]*use [Vv]ariables" "$f"
done
```

- [ ] **Step 2: Migrate each per spec.**

For files without `use variables`, skip (no work). For files with imports, replace with state-only.

- [ ] **Step 3: VG + commit (single commit for cluster).**

```bash
git add src/crop/wofost_soil_*.f90 # + callers
git commit -m "refactor(gr-crop B14): wofost_soil_* cluster — drop use variables

Per-file migration of WOFOST nutrient sub-modules to state%nutrients +
state%soilwater + state%crop refs.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task B15: Migrate `meteodt.f90` rain timing site (line 143)

**Files:** `src/atmosphere/meteodt.f90`

Site at line 143: `swrain, raintab, wet, nmrain, rainamount, arai, rainfluxarray, raintimearray`.

- [ ] **Step 1: Inspect site.**

```bash
sed -n '140,160p' src/atmosphere/meteodt.f90
```

- [ ] **Step 2: Replace import.**

```fortran
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
```

- [ ] **Step 3: Add config to signature.**

ProcessRainEvents likely needs config arg for `swrain` (config%meteo%swrain) and `raintab` (config%meteo%raintab). State for `nmrain`/`rainamount`/`rainfluxarray`/`raintimearray`/`wet`/`arai`.

- [ ] **Step 4: Migrate body references.**

| Was | Now |
|---|---|
| `swrain` | `config%meteo%swrain` |
| `raintab` | `config%meteo%raintab` |
| `wet`, `arai` | `state%atmosphere%X` |
| `nmrain`, `rainamount`, `rainfluxarray`, `raintimearray` | `state%atmosphere%X` |

- [ ] **Step 5: Update callers.**

```bash
grep -rn "call ProcessRainEvents\b\|call processrainevents\b" src/
```

- [ ] **Step 6: VG + commit.**

```bash
git commit -m "refactor(gr-crop B15): meteodt.f90 rain timing site — drop use variables

swrain/raintab → config%meteo. nmrain/rainamount/rainfluxarray/
raintimearray/wet/arai → state%atmosphere. Closes GR-ATM/Arc 7 rain
timing deferral.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task B16: Verify all crop files clean

**Files:** None (verification only).

- [ ] **Step 1: Grep all crop files.**

```bash
grep -rn "^[[:space:]]*use [Vv]ariables" src/crop/ | head -20
```

Expected: empty or only documented narrow deferrals.

- [ ] **Step 2: Grep rain timing site in meteodt.**

```bash
grep -n "^[[:space:]]*use variables" src/atmosphere/meteodt.f90
```

Should reduce by one site (line 143 retired by B15).

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: No commit** unless additional deferrals need documenting.

---

### Task B17: Phase B close marker

**Files:** None.

- [ ] **Step 1: Clean rebuild + check-full** — 5/5 byte-for-byte.

- [ ] **Step 2: Sanity greps:**
  - `grep -rn "^[[:space:]]*use [Vv]ariables" src/crop/` — empty (or only documented deferrals)
  - `grep -n "state%crop%common\|state%crop%wofost\|state%crop%grass\|state%crop%fixed" src/crop/cropgrowth.f90 | wc -l` — ≥20 (proves reads migrated)
  - `grep -n "state%nutrients" src/crop/management_soil.f90 | wc -l` — ≥5

- [ ] **Step 3: Phase B close marker.**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-crop): Phase B complete — 12-file reader cutover + rain timing migrated

All crop files migrated off bare use variables. meteodt's rain timing
site migrated. check-full 5/5 byte-for-byte. Closes GR-ATM deferrals:
GET_MAX_RESP_FACTOR state threading, irrigation gird/isua write-side,
rain timing schema.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Phase C — Codebase Sweep + Retirement (Tasks C1–C16)

**Scope discipline:** Phase C does NOT touch any other `use variables` symbols. Surgical replacement only — replace bare-name reads with state references; narrow `use variables, only:` lists by removing migrated names.

### Symbol replacement reference (Phase C)

Same as Phase B reference table, plus:
- Atmosphere multi-consumer globals (arad/atmn/atmx/ahum/awin/arai/aetr/wet/atav/tav) — `state%atmosphere%X`

### Task C1: Cluster sweep — `src/core/`

**Files:** `src/core/swap_mod.f90`, `initialize.f90`, `timecontrol_mod.f90`, `swap_bmi_mod.f90`, `swap_capi_mod.f90`

- [ ] **Step 1: Per-file inventory.**

```bash
for f in src/core/{swap_mod,initialize,timecontrol_mod,swap_bmi_mod,swap_capi_mod}.f90; do
  echo "=== $f ==="
  grep -nE "\b(daycrop|dvs|tsum|rd|rdpot|wlv|wst|wrt|wso|tagp|plwt|cftb|fom|bio|hum|nh4|no3|cnflv|cnfst|nuptot|nflux|nmrain|rainamount|rainfluxarray|raintimearray|arad|atmn|atmx|ahum|awin|arai|aetr|wet|atav|tav)\b" "$f" | grep -v "state%\|config%\|! \|integer ::\|real(8) ::\|real(real64) ::" | head -30
done
```

- [ ] **Step 2: Apply substitutions per reference table.**

**Caveat:** swap_mod's Phase A dual-write blocks still need legacy globals as source. Those imports STAY until retirement (Task C12+).

- [ ] **Step 3: VG + commit.**

```bash
git add src/core/*.f90
git commit -m "sweep(gr-crop C1): src/core/ — crop/nutrient/rain refs to state%X

Files retain use variables for unrelated globals (and dual-write imports
in swap_mod, retired in C12+).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task C2: Cluster sweep — `src/io/`

**Files:** `src/io/swap_csv_output.f90`, `src/io/swapoutput.f90`, `src/io/readmeteo.f90`

- [ ] **Step 1: Per-file inventory.**

Same grep pattern as C1.

- [ ] **Step 2: Apply substitutions.**

swap_csv_output + swapoutput likely read many crop/atmosphere fields for output. readmeteo writes meteo arrays — note write sites for C11.

- [ ] **Step 3: VG + commit.**

```bash
git commit -m "sweep(gr-crop C2): src/io/ — crop/nutrient/atmosphere refs to state%X

Output writers and readmeteo reads migrated.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task C3: Cluster sweep — `src/atmosphere/`

**Files:** `src/atmosphere/meteoday.f90`, `meteodt.f90`, `et.f90`, `interception.f90`

- [ ] **Step 1: Inventory.**

Same grep.

- [ ] **Step 2: Apply substitutions.**

meteoday + meteodt may have residual dual-write retentions (from GR-ATM Phase B.5). After C1+C2, more can be retired. Check `state%atmosphere%tav` write line — if no remaining legacy consumers, drop the dual-write.

- [ ] **Step 3: VG + commit.**

```bash
git commit -m "sweep(gr-crop C3): src/atmosphere/ — close residual dual-write retentions

After crop migration, tav/atav/arad/etc. have fewer legacy consumers.
Drop dual-write retentions where readers all migrated.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task C4: Cluster sweep — `src/boundary/`

**Files:** `src/boundary/boundbottom.f90`, `boundtop.f90`

- [ ] **Step 1: Inventory** any crop/nutrient bleed.

- [ ] **Step 2: Apply substitutions.**

- [ ] **Step 3: VG + commit.**

```bash
git commit -m "sweep(gr-crop C4): src/boundary/ — minimal crop/nutrient refs migrated

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task C5: Cluster sweep — `src/heat/`

**Files:** `src/heat/temperature.f90`, `src/heat/frozencond.f90`

- [ ] **Step 1: Inventory.**

- [ ] **Step 2: Apply substitutions** (likely minimal).

- [ ] **Step 3: VG + commit.**

---

### Task C6: Cluster sweep — `src/soil/`

**Files:** `src/soil/soilgrid.f90`, `soilhydraulics.f90`, `waterbalance.f90`

- [ ] **Step 1: Inventory.**

`soilhydraulics` may read `rd` for root depth feedback. Migrate to `state%crop%common%rd`.

- [ ] **Step 2: Apply substitutions.**

- [ ] **Step 3: VG + commit.**

```bash
git commit -m "sweep(gr-crop C6): src/soil/ — root depth + crop refs to state%crop

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task C7: Cluster sweep — `src/drainage/`

**Files:** `src/drainage/drainage.f90`, `surfacewater.f90`, `divdra.f90`

- [ ] **Step 1: Inventory** (likely minimal).

- [ ] **Step 2: VG + commit.**

---

### Task C8: Cluster sweep — `src/solute/`

**Files:** `src/solute/solute.f90`, `agetracer.f90`

- [ ] **Step 1: Inventory** — solute may read crop runtime for uptake feedback.

- [ ] **Step 2: Apply substitutions.**

- [ ] **Step 3: VG + commit.**

```bash
git commit -m "sweep(gr-crop C8): src/solute/ — crop refs to state%crop

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task C9: Audit pass — find stale reads

**Files:** Audit + targeted fixes. MANDATORY GATE before retirements.

- [ ] **Step 1: Crop common audit.**

```bash
grep -rnE "\b(daycrop|dvs|tsum|rd|rdpot|rdm|rri|rdi|rdc|ch|cf|laipot|cuptgraz|cuptgrazpot|HarLosOrm_tot|daynrsta|daynrend|swcrp|icrop|flcropcalendar|flcropoutput|flcropnut|flharvestday|swend)\b" src/ --include="*.f90" \
  | grep -v "state%\|config%\|! \|variables.f90\|use variables\|integer ::\|real(8) ::\|real(real64) ::\|logical :: \|character\|^src/core/initialize" | head -40
```

- [ ] **Step 2: Crop WOFOST/grass/fixed audit.**

```bash
grep -rnE "\b(wlv|wlvpot|wst|wstpot|wrt|wrtpot|wso|wsopot|twlv|twst|tagp|tagppot|tagpt|tagptpot|cwdm|cwdmpot|pgass|pgasspot|dwlv|dwlvpot|dwst|dwstpot|dwrt|dwrtpot|dwso|dwsopot|dwlvCrop|dwlvSoil|plossdm|lossdm|swbulb|wbl|wblpot|dwbl|dwblpot|plwt|cftb)\b" src/ --include="*.f90" \
  | grep -v "state%\|config%\|! \|variables.f90\|use variables\|integer ::\|real(8) ::\|real(real64) ::\|logical :: \|character\|^src/core/initialize" | head -40
```

- [ ] **Step 3: Nutrients audit.**

```bash
grep -rnE "\b(fom|bio|hum|nh4|no3|cnflv|cnfst|nuptot|nflux)\b" src/ --include="*.f90" \
  | grep -v "state%\|config%\|! \|variables.f90\|use variables\|integer ::\|real(8) ::\|real(real64) ::\|character\|^src/core/initialize" | head -40
```

- [ ] **Step 4: Rain timing audit.**

```bash
grep -rnE "\b(nmrain|rainamount|rainfluxarray|raintimearray|raintab)\b" src/ --include="*.f90" \
  | grep -v "state%\|config%\|! \|variables.f90\|use variables\|integer ::\|real(8) ::\|real(real64) ::\|^src/core/initialize" | head -20
```

- [ ] **Step 5: Atmosphere multi-consumer audit.**

```bash
grep -rnE "\b(arad|atmn|atmx|ahum|awin|arai|aetr|wet|atav|tav)\b" src/ --include="*.f90" \
  | grep -v "state%\|config%\|! \|variables.f90\|use variables\|integer ::\|real(8) ::\|real(real64) ::\|character\|^src/core/initialize" | head -40
```

- [ ] **Step 6: Patch each stale read in separate small commits.**

```bash
git add <file>
git commit -m "audit(gr-crop): <file> — patch stale <symbol> reference

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

- [ ] **Step 7: Re-run all 5 audits until clean.**

- [ ] **Step 8: Audit close marker.**

```bash
git commit --allow-empty -m "audit(gr-crop): codebase sweep clean — Phase C retirements unblocked

All 5 audit greps clean. Ready for global deletion.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

**MANDATORY GATE — do not advance to C10 until clean.**

---

### Task C10: Drop legacy writes from `config_to_variables.f90` for crop/nutrients/rain

**Files:** `src/io/toml/config_to_variables.f90`, `src/core/swap_mod.f90`

For each migrated symbol where state mirror reads from config:
- Rewrite the dual-write in swap_mod to source from `config%X` directly.
- Delete the legacy adapter write.

- [ ] **Step 1: Inventory adapter writes.**

```bash
grep -nE "^[[:space:]]*\b(daycrop|dvs|tsum|...|fom|...|raintab|nmrain|...)\b\s*(\(|=)" src/io/toml/config_to_variables.f90 | head -40
```

- [ ] **Step 2: Per-symbol rewrite + delete legacy.**

For each clean source path:
```fortran
! Before:
daycrop = config%crop%daycrop   ! legacy global write
! After:
state%crop%common%daycrop = config%crop%daycrop   ! direct state population
```

Delete the legacy line.

- [ ] **Step 3: VG + commit.**

```bash
git commit -m "retire(gr-crop C10): drop legacy crop/nutrient/rain writes from adapter

State mirrors sourced directly from config%X. Legacy adapter writes
removed. Sets up Tasks C12+ for declaration removal.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task C11: Drop legacy writes from `readmeteo.f90` for rain timing

**Files:** `src/io/readmeteo.f90`

Per-year reload writes rain timing arrays. Migrate to write state directly.

- [ ] **Step 1: Inspect rain timing writes.**

```bash
grep -nE "^[[:space:]]*\b(nmrain|rainamount|rainfluxarray|raintimearray)\b\s*(\(|=)" src/io/readmeteo.f90 | head
```

- [ ] **Step 2: Replace with state writes.**

```fortran
! Before:
nmrain = ...
! After:
state%atmosphere%nmrain = ...
```

Remove the GR-CROP A5.5 mirror retention in meteodt if readmeteo is the sole rain timing populator (verify).

- [ ] **Step 3: VG + commit.**

```bash
git commit -m "retire(gr-crop C11): drop legacy rain timing writes from readmeteo

State writes are sole. Per-year reloads now populate state%atmosphere
directly.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task C12: Retire write-only globals from `variables.f90`

**Files:** `src/core/variables.f90`, `src/core/initialize.f90`, plus compile-surfaced patches

For each migrated symbol where the bare global has NO remaining readers:
- Tombstone the declaration in `variables.f90`.
- Drop zero-fills in `initialize.f90`.

- [ ] **Step 1: Identify write-only candidates.**

For each candidate, run:
```bash
for sym in daycrop dvs tsum rd rdpot ch cf wlv wst wrt wso tagp pgass fom bio hum nh4 no3 cnflv cnfst nuptot nflux nmrain rainamount rainfluxarray raintimearray; do
  hits=$(grep -rnE "\b${sym}\b" src/ --include="*.f90" \
    | grep -v "state%\|config%\|! \|variables.f90\|use variables, only:\|initialize.f90\|config_to_variables\|swap_mod.f90\|integer ::\|real(8) ::\|real(real64) ::" | wc -l)
  echo "$sym: $hits non-write hits"
done
```

Symbols with 0 hits are candidates.

- [ ] **Step 2: Tombstone declarations.**

In `src/core/variables.f90`, replace each declaration with a tombstone comment:
```fortran
! [SS-GR-CROP C12] <sym> retired — state%X%Y is sole home
```

- [ ] **Step 3: Drop zero-fills in initialize.f90.**

- [ ] **Step 4: Iterative build cycles** — patch compile errors as they surface.

- [ ] **Step 5: VG + commit.**

```bash
git commit -m "retire(gr-crop C12): delete write-only crop/nutrient/rain globals

<list of retired symbols, ~30-50>. state%crop / state%nutrients /
state%atmosphere is sole home.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task C13: Retire atmosphere multi-consumer globals if last consumers migrated

**Files:** `src/core/variables.f90`, `src/core/initialize.f90`

Check whether arad/atmn/atmx/ahum/awin/arai/aetr/wet/atav/tav have any remaining readers:

```bash
for sym in arad atmn atmx ahum awin arai aetr wet atav tav; do
  hits=$(grep -rnE "\b${sym}\b" src/ --include="*.f90" \
    | grep -v "state%\|config%\|! \|variables.f90\|use variables, only:\|initialize.f90\|config_to_variables\|swap_mod.f90\|integer ::\|real(8) ::\|real(real64) ::" | wc -l)
  echo "$sym: $hits non-write hits"
done
```

- [ ] **Step 1: For symbols with 0 hits, retire.**

Tombstone in variables.f90; drop zero-fills.

- [ ] **Step 2: For symbols with remaining hits, document deferral.**

Note inheritance to Arc 9 with explicit consumer list.

- [ ] **Step 3: VG + commit.**

```bash
git commit -m "retire(gr-crop C13): atmosphere multi-consumer globals — partial retirement

<list retired>. Symbols with remaining readers deferred to Arc 9 with
explicit consumer documentation.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task C14: Drop swap_mod A12 init-time seeding for retired symbols

**Files:** `src/core/swap_mod.f90`

After C10+C12+C13, the [SS-GR-ATM A12] block + [SS-GR-CROP A14-A17] blocks become redundant for retired symbols. State mirrors now source directly from config (C10) or are written by primary writers (readmeteo, etc.).

- [ ] **Step 1: Inspect each dual-write line in swap_mod.**

```bash
grep -nE "SS-GR-ATM A|SS-GR-CROP A" src/core/swap_mod.f90 | head -30
```

For each line `state%X%Y = legacy_sym`, check if `legacy_sym` is still declared in variables.f90. If retired (per C12/C13), the line is dead — delete it. Also drop the corresponding entry from the `use variables, only:` list at top.

- [ ] **Step 2: VG + commit.**

```bash
git commit -m "retire(gr-crop C14): drop swap_mod A12/A14-A17 dual-writes for retired symbols

Phase A.5 dual-writes at write sites + state primary writers (readmeteo,
etc.) are now sole population paths. Init-time seeding for retired
globals is dead code.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task C15: Update Arc 9 plan with residual deferrals

**Files:** `docs/superpowers/plans/2026-05-13-globals-retirement-roadmap.md`

- [ ] **Step 1: Find the Arc 9 section.**

- [ ] **Step 2: Document residuals.**

```markdown
**Inherited from GR-CROP (as of 2026-05-14):**
- ETSine astronomical scratchpad (rad/daylp/difpp/atmtr/dsinbe/tsunrise_atm/tsunset_atm/lat) — schema decision per-symbol
- boundtop.f90 swkmean/swredu/flrunon/runonarr — config threading through soilhydraulics→headcalc (Arc 5 territory; not addressed by GR-CROP because not crop-related)
- swap_mod.f90 swinco=3 inline block + transient *_init_buf reads — strangler-fig leftovers in swap_init_from_loaded_config
- Remaining variables.f90 entries (those with deferred consumers; per audit)
- config_to_variables.f90 (1802 lines), initialize.f90 (845 lines) — deletion when all consumers migrate
```

- [ ] **Step 3: Commit.**

```bash
git add docs/superpowers/plans/2026-05-13-globals-retirement-roadmap.md
git commit -m "docs(roadmap): GR-CROP close — Arc 9 inherits final retirements

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task C16: Final verification + arc-complete marker

**Files:** None.

- [ ] **Step 1: Clean rebuild + check-full.**

```bash
rm -rf builddir && pixi run build-linux && pixi run check-full
```

Expected: **5/5 byte-for-byte**.

- [ ] **Step 2: BMI + cffi-demo suites.**

```bash
pixi run -e test test-bmi
pixi run -e test test-cffi-demo
```

- [ ] **Step 3: Final sanity greps.**

```bash
# Crop files clean
grep -rn "^[[:space:]]*use [Vv]ariables" src/crop/ | head

# State%crop usage spread
grep -rl "state%crop%" src/ | wc -l   # expect ≥10

# State%nutrients usage
grep -rl "state%nutrients%" src/ | wc -l   # expect ≥3

# variables.f90 shrinkage
wc -l src/core/variables.f90
```

- [ ] **Step 4: Arc-complete marker.**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-crop): GR-CROP complete — crop+nutrients+rain timing+cross-arc closures

check-full 5/5 byte-for-byte. BMI + cffi-demo passing. pFUnit all-pass.

Net retirement / scaffolding:
- NEW crop_state_t sub-records (common/fixed/wofost/grass) — ~50 fields
- NEW state%nutrients subrecord — ~15 fields
- state%atmosphere rain timing (4 fields)
- config%meteo raintab (1 field)
- ~120-160 Phase A.5 runtime dual-write sites (crop runtime + nutrients
  + rain timing)
- 12 crop files migrated off bare use variables
- meteodt.f90 rain timing site migrated
- Closed GR-ATM deferrals: GET_MAX_RESP_FACTOR state threading;
  irrigation gird/isua write-side; rain timing schema
- ~30-50 write-only globals retired from variables.f90
- Atmosphere multi-consumer globals (where consumers all migrated)

Inherited to Arc 9 (documented in roadmap): ETSine astronomical scratchpad;
boundtop config threading through soilhydraulics; swap_mod swinco=3 inline
block + transient *_init_buf reads; config_to_variables + initialize.f90
+ variables.f90 final deletion.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Plan Self-Review

**Spec coverage:**
- ✅ crop_state_t restructure (4 sub-records) → Tasks 2-9
- ✅ nutrients_state_t introduction → Tasks 10-11
- ✅ Rain timing schema → Tasks 12-13
- ✅ Dual-write seeding → Tasks 14-17
- ✅ Phase A close → Task 18
- ✅ Runtime dual-writes → Tasks A5.1-A5.6
- ✅ 12 crop files reader cutover → Tasks B1-B16
- ✅ Rain timing reader migration in meteodt → Task B15
- ✅ Phase B close → Task B17
- ✅ Codebase sweep clusters → Tasks C1-C8
- ✅ Audit pass → Task C9
- ✅ Adapter retirement → Tasks C10-C11
- ✅ Global declaration deletion → Tasks C12-C14
- ✅ Arc 9 roadmap update → Task C15
- ✅ Final verification → Task C16

**Placeholder scan:** Several "audit-driven" markers acknowledge that exact field lists for `crop_common`/`crop_wofost`/`crop_grass`/`crop_fixed`/`nutrients` will be determined by the implementer at Tasks A6-A11. This is a deliberate plan design — the spec lists representative fields; implementer fills in completely per audit. No "TBD" / "implement later" / "similar to Task N" patterns. Inline `<implementer documents audit findings>` markers are commit-message customization points.

**Type consistency:** crop_state_t sub-record naming consistent (common/fixed/wofost/grass). state%crop top-level fields preserved from GR-ATM. nutrients_state_t fields consistent across Tasks 11, 16, A5.4, B13. Rain timing fields consistent across Tasks 12, 17, A5.5, B15.

**Known soft spots flagged inline:**
- Field lists in Tasks A6-A9 + A11 are audit-driven; final list depends on grep results.
- Phase A.5 zero-reset coverage is critical (lesson from GR-ATM `8fb76cf`).
- Phase B cropgrowth split across 3 tasks (B5/B6/B7) due to file size (4937L).
- Controller intervenes inline if subagents defer due to missing infrastructure (per user directive).

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-05-14-globals-crop.md`. Two execution options:

**1. Subagent-Driven (recommended)** — Fresh subagent per task, two-stage review (spec compliance + code quality) between each. Same pattern as GR-UTILS + GR-BH + GR-ATM. Controller intervenes inline if subagents defer due to missing infrastructure (per user directive).

**2. Inline Execution** — Execute in this session via executing-plans, batch execution with checkpoints. Higher main-context pressure given this is the largest arc.

Which approach?
