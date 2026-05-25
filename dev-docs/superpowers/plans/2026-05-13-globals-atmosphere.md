# GR-ATM Implementation Plan — Atmosphere Globals Retirement + crop_state introduction

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Retire `use variables` from `et.f90`, `interception.f90`, `meteoday.f90`, `meteodt.f90` (4 atmosphere files) by introducing `state%crop` (NEW `crop_state_t`, 12 fields), extending `state%atmosphere` with ~50 fields in 6 blocks, extending `config%meteo` with `cofred`, and sweeping ~30 readers across the codebase. Closes prior-arc deferrals: `Tav/atav` (temperature.f90), `swkmean/swredu/flrunon/runonarr` (boundtop), `cofred` (et.f90).

**Architecture:** Three sub-phases inside one arc. Phase A is additive schema + dual-write (no reader changes). Phase B cuts the 4 target files over to state/config and introduces type-bound `state%atmosphere%init`. Phase C sweeps every remaining reader and deletes the bare globals from `variables.f90` (compile surfaces any missed reader). Each task ends with the full regression gate; each phase ends with `check-full` 5/5 byte-for-byte.

**Tech Stack:** Fortran 2018, gfortran, Meson, pFUnit 4.15, pixi, Python regression tests.

**Spec:** `docs/superpowers/specs/2026-05-13-globals-atmosphere-design.md`

**Builds on:** GR-BH (Arc 2+3 + mesh, completed 2026-05-13 — commit `a65bdf3`).

**Verification gate per task (VG)** — per memory `feedback_state_schema_clean_rebuild.md` + `feedback_per_task_regression_gate.md`:

```bash
rm -rf builddir && pixi run build-linux \
  && pixi run test-pfunit \
  && pixi run -e test python tests/regression/test_output_regression.py \
       hupselbrook surfacewater salinitystress grassgrowth
```

Expected: clean build, pFUnit all-pass, regression **4/4 byte-for-byte**. Any deviation = task NOT complete; fix or escalate. Do not commit on red.

**Phase close gate:** `pixi run check-full` 5/5 byte-for-byte.

**Subagent dispatch convention:** every implementer prompt MUST include this VG verbatim. `rm -rf builddir` is mandatory because Meson's incremental build does not propagate `.mod` deps across the `swap_modern`/`swap_legacy` static-library boundary when state schema changes.

---

## File Structure

**New module:**
- `src/state/crop_state.f90` — `crop_state_t` (12 fields + type-bound `init`)

**Extended modules:**
- `src/state/swap_state.f90` — adds `type(crop_state_t) :: crop` field
- `src/state/atmosphere_state.f90` — +~50 fields across 6 blocks + type-bound `init`
- `src/config/meteorology_config.f90` — +`cofred` field
- `src/config/simulation_config.f90` — possibly +`flrunon` + `runonarr` (if not yet present)

**Migrated readers (full `use variables` retirement):**
- `src/atmosphere/et.f90` (1 `use variables` site at line 614)
- `src/atmosphere/interception.f90` (4 sites: lines 40, 89, 171, 409)
- `src/atmosphere/meteoday.f90` (3 sites: lines 102, 289, 543 — long `only:` lists)
- `src/atmosphere/meteodt.f90` (4 sites: lines 58, 139, 358, 456 — mostly bare)

**Codebase sweep (atmosphere / crop globals swap to state):**
- `src/core/`: `swap_mod.f90`, `initialize.f90`, `timecontrol_mod.f90`, `swap_bmi_mod.f90`, `swap_capi_mod.f90`
- `src/io/`: `swap_csv_output.f90`, `swapoutput.f90`, `readmeteo.f90`
- `src/boundary/`: `boundbottom.f90`, `boundtop.f90` (close prior-arc deferrals: flrunon/runonarr/Tav)
- `src/heat/`: `temperature.f90` (close Tav/atav deferral)
- `src/soil/`: `soilgrid.f90`, `soilhydraulics.f90`, `waterbalance.f90`
- `src/drainage/`: `drainage.f90`, `surfacewater.f90`, `divdra.f90`
- `src/crop/`: all 9 crop files
- `src/solute/`: `solute.f90`, `agetracer.f90`

**Populator + adapter shrinks:**
- `src/io/toml/config_to_variables.f90` — adds Phase A dual-writes; legacy global writes deleted in Phase C
- `src/io/readmeteo.f90` — meteo arrays populator; routes into state in Phase C
- `src/core/initialize.f90` — loses zero-fills for deleted globals

**Build:**
- `meson.build` — add `src/state/crop_state.f90` to legacy sources, **before** `src/state/swap_state.f90`
- `tests/unit/meson.build` — add to `pfunit_extra_sources`, same ordering

---

## Phase A — Schema + Dual-Write (Tasks 1–13)

Additive only. After Phase A: every legacy global still load-bearing; state mirrors populated; no reader migrated.

### Task 1: Pre-flight baseline

**Files:** No code changes.

- [ ] **Step 1: Clean rebuild.** Run: `rm -rf builddir && pixi run build-linux` — must succeed.
- [ ] **Step 2: pFUnit baseline.** Run: `pixi run test-pfunit` — record count.
- [ ] **Step 3: check-full baseline.** Run: `pixi run check-full` — 5/5 pass.
- [ ] **Step 4: Empty marker commit.**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-atm): pre-flight baseline — atmosphere globals retirement begins

check-full 5/5. Baseline locked before crop_state introduction +
atmosphere_state extensions + 4-file reader migration + codebase sweep.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Create `crop_state_mod`

**Files:**
- Create: `src/state/crop_state.f90`
- Modify: `src/state/swap_state.f90` (add `type(crop_state_t) :: crop` field)
- Modify: `meson.build` (add new file to legacy sources, **before** `src/state/swap_state.f90`)
- Modify: `tests/unit/meson.build` (add to `pfunit_extra_sources`, same ordering)

- [ ] **Step 1: Create the new module.**

```fortran
!> @file crop_state.f90
!! SS-GR-ATM: typed crop runtime state. Foundation for Arc 8 (crop
!! cluster). Arc 4 (atmosphere) introduces this module to host the
!! ~12 crop-runtime symbols read by atmosphere code; Arc 8 later
!! migrates all remaining crop readers.
module crop_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: crop_state_t

   type :: crop_state_t
      real(real64) :: lai             = 0.0_real64    !! leaf area index (-)
      real(real64) :: kdif            = 0.0_real64    !! extinction for diffuse light (-)
      real(real64) :: kdir            = 0.0_real64    !! extinction for direct light (-)
      real(real64) :: cofab           = 0.0_real64    !! interception coefficient (cm)
      real(real64) :: cfbs            = 1.0_real64    !! bare-soil ET factor (-)
      integer      :: swcf            = 0             !! crop factor switch
      integer      :: swcfbs          = 0             !! bare-soil factor switch
      real(real64) :: gird            = 0.0_real64    !! gross irrigation depth (cm)
      logical      :: flCropEmergence = .false.       !! crop emerged flag
      real(real64) :: et0             = 0.0_real64    !! potential ET (cm/d)
      real(real64) :: ew0             = 0.0_real64    !! potential evap wet crop (cm/d)
      real(real64) :: es0             = 0.0_real64    !! potential evap bare soil (cm/d)
   contains
      procedure :: init => crop_state_init
   end type crop_state_t

contains

   subroutine crop_state_init(self, crop_cfg)
      use crop_config_mod, only: crop_config_t
      class(crop_state_t),  intent(inout) :: self
      type(crop_config_t),  intent(in)    :: crop_cfg
      ! Defaults retained from type initializers; runtime seeders update
      ! each timestep. crop_cfg arg reserved for future seed migration.
   end subroutine crop_state_init

end module crop_state_mod
```

- [ ] **Step 2: Wire into `swap_state_mod`.** Edit `src/state/swap_state.f90`:
  - Add `use crop_state_mod, only: crop_state_t` to module imports.
  - Add `type(crop_state_t) :: crop` after `type(mesh_state_t) :: mesh` (or wherever fits cleanly in the type body).

- [ ] **Step 3: Add to `meson.build`.** Find legacy sources list, add `'src/state/crop_state.f90',` **immediately before** `'src/state/swap_state.f90'`.

- [ ] **Step 4: Add to `tests/unit/meson.build`.** Find `pfunit_extra_sources`, add `'../../src/state/crop_state.f90',` before `'../../src/state/swap_state.f90'`.

- [ ] **Step 5: Run VG.** Expected: clean build, pFUnit pass, regression 4/4.

- [ ] **Step 6: Commit.**

```bash
git add src/state/crop_state.f90 src/state/swap_state.f90 meson.build tests/unit/meson.build
git commit -m "$(cat <<'EOF'
schema(gr-atm): add crop_state_mod + state%crop field

NEW crop_state_t with 12 runtime fields (lai/kdif/kdir/cofab/cfbs/swcf/
swcfbs/gird/flCropEmergence/et0/ew0/es0). Type-bound init(config%crop).
No reader uses it yet — schema is additive. Foundation for Arc 8.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: Extend `atmosphere_state_t` — Block 1 (8 daily meteo arrays)

**Files:** `src/state/atmosphere_state.f90`

- [ ] **Step 1: Add fields.** Place in a clearly-marked block near existing daily-flow fields:

```fortran
      ! [SS-GR-ATM A3] Block 1: daily meteo input arrays (366-sized fixed; matches legacy)
      real(real64) :: arad(366) = 0.0_real64    !! daily radiation input
      real(real64) :: atmn(366) = 0.0_real64    !! daily min temperature
      real(real64) :: atmx(366) = 0.0_real64    !! daily max temperature
      real(real64) :: ahum(366) = 0.0_real64    !! daily humidity
      real(real64) :: awin(366) = 0.0_real64    !! daily wind speed
      real(real64) :: arai(366) = 0.0_real64    !! daily rainfall
      real(real64) :: aetr(366) = 0.0_real64    !! daily ETref
      real(real64) :: wet(366)  = 0.0_real64    !! daily wet-fraction
```

- [ ] **Step 2: Run VG.**

- [ ] **Step 3: Commit.**

```bash
git add src/state/atmosphere_state.f90
git commit -m "$(cat <<'EOF'
schema(gr-atm): extend atmosphere_state_t with daily meteo arrays (Block 1)

arad/atmn/atmx/ahum/awin/arai/aetr/wet — 8 fixed-size(366) daily input
arrays. Unpopulated until Phase A dual-write tasks.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: Extend `atmosphere_state_t` — Block 2 (5 sub-daily arrays)

**Files:** `src/state/atmosphere_state.f90`

- [ ] **Step 1: Add fields.**

```fortran
      ! [SS-GR-ATM A4] Block 2: sub-daily detailed arrays (96-sized; swmetdetail==1)
      real(real64) :: atav(96)  = 0.0_real64    !! sub-daily air temp
      real(real64) :: epot(96)  = 0.0_real64    !! sub-daily potential evap
      real(real64) :: tpot(96)  = 0.0_real64    !! sub-daily potential transp
      real(real64) :: grain(96) = 0.0_real64    !! sub-daily gross rain
      real(real64) :: nrain(96) = 0.0_real64    !! sub-daily net rain
```

- [ ] **Step 2: Run VG.**

- [ ] **Step 3: Commit.**

```bash
git add src/state/atmosphere_state.f90
git commit -m "schema(gr-atm): extend atmosphere_state_t with sub-daily arrays (Block 2)

atav/epot/tpot/grain/nrain — 5 fixed-size(96) sub-daily arrays
populated when swmetdetail==1.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 5: Extend `atmosphere_state_t` — Block 3 (9 derived scalars)

**Files:** `src/state/atmosphere_state.f90`

- [ ] **Step 1: Add fields.**

```fortran
      ! [SS-GR-ATM A5] Block 3: derived meteo scalars
      real(real64) :: Tav        = 0.0_real64    !! daily mean air temp
      real(real64) :: tavd       = 0.0_real64    !! daytime mean air temp
      real(real64) :: rh         = 0.0_real64    !! relative humidity
      integer      :: daynrfirst = 0             !! first day in meteo year
      integer      :: daynrlast  = 0             !! last day in meteo year
      real(real64) :: atmin7(7)  = 0.0_real64    !! 7-day min-temp buffer
      integer      :: nofd       = 0             !! current day-of-running-avg
      real(real64) :: teprrain   = 0.0_real64    !! threshold rain temp
      real(real64) :: teprsnow   = 0.0_real64    !! threshold snow temp
```

- [ ] **Step 2: Run VG.**

- [ ] **Step 3: Commit.**

```bash
git add src/state/atmosphere_state.f90
git commit -m "schema(gr-atm): extend atmosphere_state_t with derived scalars (Block 3)

Tav/tavd/rh/daynrfirst/daynrlast/atmin7/nofd/teprrain/teprsnow — 9
derived meteo scalars. Closes Tav/atav deferral path for temperature.f90.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 6: Extend `atmosphere_state_t` — Block 4 (8 interception fields)

**Files:** `src/state/atmosphere_state.f90`

- [ ] **Step 1: Locate `magrs` import.** Check that `use swap_array_dimensions, only: magrs` is in atmosphere_state.f90 imports. If not, add it (needed for table sizes).

- [ ] **Step 2: Add fields.**

```fortran
      ! [SS-GR-ATM A6] Block 4: interception state/params
      real(real64) :: siccapact          = 0.0_real64
      real(real64) :: fimin              = 0.0_real64
      integer      :: isua               = 0
      real(real64) :: avevaptb(2*magrs)  = 0.0_real64  !! actual evap table
      real(real64) :: avprectb(2*magrs)  = 0.0_real64  !! actual precip table
      real(real64) :: pfreetb(2*magrs)   = 0.0_real64  !! free throughfall table
      real(real64) :: pstemtb(2*magrs)   = 0.0_real64  !! stemflow table
      real(real64) :: scanopytb(2*magrs) = 0.0_real64  !! canopy storage table
```

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/state/atmosphere_state.f90
git commit -m "schema(gr-atm): extend atmosphere_state_t with interception (Block 4)

siccapact/fimin/isua + avevaptb/avprectb/pfreetb/pstemtb/scanopytb
(5 tables sized 2*magrs). Unpopulated until dual-write tasks.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 7: Extend `atmosphere_state_t` — Block 5 (10 Runoff-CN fields)

**Files:** `src/state/atmosphere_state.f90`

- [ ] **Step 1: Add fields.**

```fortran
      ! [SS-GR-ATM A7] Block 5: Runoff-CN method state + tables
      real(real64) :: CNref             = 0.0_real64
      real(real64) :: CNdry             = 0.0_real64
      real(real64) :: CNwet             = 0.0_real64
      real(real64) :: ThetaRef          = 0.0_real64
      real(real64) :: Runoff_CN         = 0.0_real64
      real(real64) :: wc_cor            = 0.0_real64
      real(real64) :: wc10              = 0.0_real64
      integer      :: iCNtab            = 0
      real(real64) :: CNtimTAB(2*magrs) = 0.0_real64
      real(real64) :: CNrefTAB(2*magrs) = 0.0_real64
```

- [ ] **Step 2: Run VG.**

- [ ] **Step 3: Commit.**

```bash
git add src/state/atmosphere_state.f90
git commit -m "schema(gr-atm): extend atmosphere_state_t with Runoff-CN (Block 5)

CNref/CNdry/CNwet/ThetaRef/Runoff_CN/wc_cor/wc10/iCNtab + 2 tables
sized 2*magrs.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 8: Extend `atmosphere_state_t` — Block 6 (7 output flags)

**Files:** `src/state/atmosphere_state.f90`

- [ ] **Step 1: Add fields.**

```fortran
      ! [SS-GR-ATM A8] Block 6: output flag toggles
      logical :: out_tmn = .false.
      logical :: out_tmx = .false.
      logical :: out_hum = .false.
      logical :: out_win = .false.
      logical :: out_etr = .false.
      logical :: out_wet = .false.
      logical :: out_rad = .false.
```

- [ ] **Step 2: Run VG.**

- [ ] **Step 3: Commit.**

```bash
git add src/state/atmosphere_state.f90
git commit -m "schema(gr-atm): extend atmosphere_state_t with output flags (Block 6)

out_tmn/out_tmx/out_hum/out_win/out_etr/out_wet/out_rad — 7 logical
toggles set from config, read by meteoday to control output content.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 9: Extend `config%meteo` (+ verify flrunon/runonarr placement)

**Files:**
- Modify: `src/config/meteorology_config.f90` (add `cofred`)
- Possibly modify: `src/config/simulation_config.f90` or `src/config/general_config.f90` (add `flrunon`, `runonarr`)
- Possibly modify: `src/io/toml/config_to_variables.f90` (route legacy global from config%X — verify legacy default)

- [ ] **Step 1: Add `cofred` to meteo config.**

Edit `src/config/meteorology_config.f90`. Locate near `cofredbl`/`cofredbo` (around line 17-18 per spec):

```fortran
      real(real64) :: cofred = 0.35_real64  !! reduction coefficient for evap (legacy default 0.35)
```

- [ ] **Step 2: Verify legacy default for cofred.** Run:

```bash
grep -n "cofred\s*=" src/io/toml/config_to_variables.f90 src/core/initialize.f90 src/core/variables.f90 src/atmosphere/et.f90 | head
```

If the value differs from 0.35, adjust the default to match legacy.

- [ ] **Step 3: Inspect flrunon / runonarr.**

```bash
grep -n "flrunon\|runonarr" src/config/*.f90 src/core/variables.f90 src/io/toml/config_to_variables.f90 | head -20
```

If `flrunon` and `runonarr` are NOT in any config_t, add them to `simulation_config_t` (or wherever the legacy population path indicates). For example, in `src/config/simulation_config.f90`:

```fortran
      logical :: flrunon = .false.     !! runon flag
      ! runonarr: time-series — verify size from legacy declaration
      ! e.g.: real(real64), allocatable :: runonarr(:)
```

If `flrunon` and `runonarr` ARE already in a config (e.g., already moved to `config%general` or similar), no change needed — they were already retired from `variables.f90` in a prior arc.

- [ ] **Step 4: Add adapter route for `cofred`.** Inspect `config_to_variables.f90` for where `cofred` is set:

```bash
grep -n "^[[:space:]]*cofred\s*=" src/io/toml/config_to_variables.f90
```

If a legacy line `cofred = <value>` exists, change to `cofred = config%meteo%cofred`. If no such line exists, add:

```fortran
   cofred = config%meteo%cofred  ! [SS-GR-ATM A9] add cofred config path
```

- [ ] **Step 5: Run VG.**

- [ ] **Step 6: Commit.**

```bash
git add src/config/meteorology_config.f90 src/config/simulation_config.f90 src/io/toml/config_to_variables.f90
git commit -m "$(cat <<'EOF'
schema(gr-atm): extend config%meteo with cofred + verify flrunon/runonarr

cofred added to meteorology_config_t (default 0.35). flrunon/runonarr
placement verified — <document outcome in commit body>.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 10: Dual-write Block 1 + Block 2 in swap_mod

**Files:** `src/core/swap_mod.f90` (location pattern learned from GR-BH: dual-writes go in swap_init_from_loaded_config AFTER the relevant init).

- [ ] **Step 1: Locate placement.** Find where atmosphere_init is called:

```bash
grep -n "atmosphere_init\|state%atmosphere%init\|! [SS-GR-BH" src/core/swap_mod.f90 | head -20
```

Place the dual-write block **after** `atmosphere_init` is called (so state%atmosphere is allocated/initialized).

- [ ] **Step 2: Add `use variables` imports** at the top of `swap_init_from_loaded_config` for the legacy globals to be copied:

```fortran
   use variables, only: arad, atmn, atmx, ahum, awin, arai, aetr, wet, &
                        atav, epot, tpot, grain, nrain
```

(Add these to the existing `use variables` list; do NOT remove existing entries.)

- [ ] **Step 3: Add dual-write block.**

```fortran
   ! [SS-GR-ATM A10] dual-write Block 1 (daily meteo arrays) + Block 2 (sub-daily)
   state%atmosphere%arad  = arad
   state%atmosphere%atmn  = atmn
   state%atmosphere%atmx  = atmx
   state%atmosphere%ahum  = ahum
   state%atmosphere%awin  = awin
   state%atmosphere%arai  = arai
   state%atmosphere%aetr  = aetr
   state%atmosphere%wet   = wet
   state%atmosphere%atav  = atav
   state%atmosphere%epot  = epot
   state%atmosphere%tpot  = tpot
   state%atmosphere%grain = grain
   state%atmosphere%nrain = nrain
```

**Caveat — readmeteo population:** The legacy daily arrays are populated by `readmeteo` (per-year reload). The dual-write in swap_init_from_loaded_config seeds the initial values. For per-year reloads, Phase C C11 routes the readmeteo population directly into state%atmosphere. For now (Phase A), the legacy arrays remain authoritative.

- [ ] **Step 4: Run VG.**

- [ ] **Step 5: Commit.**

```bash
git add src/core/swap_mod.f90
git commit -m "schema(gr-atm): dual-write Block 1 + Block 2 meteo arrays

Daily (arad/atmn/atmx/ahum/awin/arai/aetr/wet) and sub-daily (atav/epot/
tpot/grain/nrain) arrays mirrored from legacy globals into state%atmosphere.
readmeteo populates legacy globals; state mirror is initial snapshot.
Per-year reload retirement deferred to Phase C C11.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 11: Dual-write Blocks 3 + 4 + 5 in swap_mod

**Files:** `src/core/swap_mod.f90`

- [ ] **Step 1: Add `use variables` entries** for legacy globals in Blocks 3, 4, 5.

```fortran
   use variables, only: ..., Tav, tavd, rh, daynrfirst, daynrlast, atmin7, nofd, teprrain, teprsnow, &
                        siccapact, fimin, isua, avevaptb, avprectb, pfreetb, pstemtb, scanopytb, &
                        CNref, CNdry, CNwet, ThetaRef, Runoff_CN, wc_cor, wc10, iCNtab, CNtimTAB, CNrefTAB
```

- [ ] **Step 2: Add dual-write block** after Block 1+2 dual-write:

```fortran
   ! [SS-GR-ATM A11] dual-write Block 3 (derived scalars) + Block 4 (interception) + Block 5 (CN)
   state%atmosphere%Tav        = Tav
   state%atmosphere%tavd       = tavd
   state%atmosphere%rh         = rh
   state%atmosphere%daynrfirst = daynrfirst
   state%atmosphere%daynrlast  = daynrlast
   state%atmosphere%atmin7     = atmin7
   state%atmosphere%nofd       = nofd
   state%atmosphere%teprrain   = teprrain
   state%atmosphere%teprsnow   = teprsnow

   state%atmosphere%siccapact = siccapact
   state%atmosphere%fimin     = fimin
   state%atmosphere%isua      = isua
   state%atmosphere%avevaptb  = avevaptb
   state%atmosphere%avprectb  = avprectb
   state%atmosphere%pfreetb   = pfreetb
   state%atmosphere%pstemtb   = pstemtb
   state%atmosphere%scanopytb = scanopytb

   state%atmosphere%CNref     = CNref
   state%atmosphere%CNdry     = CNdry
   state%atmosphere%CNwet     = CNwet
   state%atmosphere%ThetaRef  = ThetaRef
   state%atmosphere%Runoff_CN = Runoff_CN
   state%atmosphere%wc_cor    = wc_cor
   state%atmosphere%wc10      = wc10
   state%atmosphere%iCNtab    = iCNtab
   state%atmosphere%CNtimTAB  = CNtimTAB
   state%atmosphere%CNrefTAB  = CNrefTAB
```

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/core/swap_mod.f90
git commit -m "schema(gr-atm): dual-write Blocks 3+4+5

Derived meteo scalars (Tav/tavd/rh/...), interception (siccapact/...
+ 5 tables), Runoff-CN method (CNref/.../CNtabs) seeded from legacy
globals.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 12: Dual-write Block 6 + crop_state seeding

**Files:** `src/core/swap_mod.f90`

- [ ] **Step 1: Add `use variables` entries** for Block 6 + crop fields.

```fortran
   use variables, only: ..., out_tmn, out_tmx, out_hum, out_win, out_etr, out_wet, out_rad, &
                        lai, kdif, kdir, cofab, cfbs, swcf, swcfbs, gird, flCropEmergence, et0, ew0, es0
```

- [ ] **Step 2: Add dual-write + crop seeding block.**

```fortran
   ! [SS-GR-ATM A12] dual-write Block 6 (output flags)
   state%atmosphere%out_tmn = out_tmn
   state%atmosphere%out_tmx = out_tmx
   state%atmosphere%out_hum = out_hum
   state%atmosphere%out_win = out_win
   state%atmosphere%out_etr = out_etr
   state%atmosphere%out_wet = out_wet
   state%atmosphere%out_rad = out_rad

   ! [SS-GR-ATM A12] seed state%crop from legacy crop globals
   state%crop%lai             = lai
   state%crop%kdif            = kdif
   state%crop%kdir            = kdir
   state%crop%cofab           = cofab
   state%crop%cfbs            = cfbs
   state%crop%swcf            = swcf
   state%crop%swcfbs          = swcfbs
   state%crop%gird            = gird
   state%crop%flCropEmergence = flCropEmergence
   state%crop%et0             = et0
   state%crop%ew0             = ew0
   state%crop%es0             = es0
```

**Note:** `flCropEmergence` may be a different type (e.g., character flag in legacy). Verify and adjust state declaration if needed.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/core/swap_mod.f90
git commit -m "schema(gr-atm): dual-write Block 6 (output flags) + crop_state seeding

state%atmosphere%out_X (7 toggles) + state%crop (12 fields) seeded
from legacy globals.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 13: Phase A close — check-full 5/5

**Files:** None (verification only).

- [ ] **Step 1: Clean rebuild.** `rm -rf builddir && pixi run build-linux`.
- [ ] **Step 2: pFUnit.** `pixi run test-pfunit` — all-pass.
- [ ] **Step 3: check-full.** `pixi run check-full` — 5/5 byte-for-byte.
- [ ] **Step 4: Sanity greps.**
  - `grep -n "state%crop%lai\s*=" src/core/swap_mod.f90` ≥1
  - `grep -n "state%atmosphere%arad\s*=" src/core/swap_mod.f90` ≥1
  - `grep -n "state%atmosphere%Tav\s*=" src/core/swap_mod.f90` ≥1
  - `grep -n "state%atmosphere%CNref\s*=" src/core/swap_mod.f90` ≥1
  - `grep -n "state%atmosphere%out_tmn\s*=" src/core/swap_mod.f90` ≥1
- [ ] **Step 5: Phase A close marker.**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-atm): Phase A complete — schema + dual-write landed

state%crop (12 fields), state%atmosphere (+~50 fields across 6 blocks),
config%meteo (+cofred) populated alongside legacy globals.
check-full 5/5 byte-for-byte. No readers migrated yet.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Phase B — Reader Cutover + atmosphere init (Tasks 14–30)

The 4 target files drop `use variables` and read state/config exclusively. Type-bound `state%atmosphere%init` introduced.

### Task 14: Introduce type-bound `state%atmosphere%init`

**Files:** `src/state/atmosphere_state.f90`, `src/core/swap_mod.f90`

- [ ] **Step 1: Inspect existing init pattern.** Run:

```bash
grep -n "atmosphere_init\|procedure :: init" src/state/atmosphere_state.f90
```

If `atmosphere_init` exists as a standalone subroutine, relocate it to type-bound. If no init exists, create one.

- [ ] **Step 2: Add type-bound procedure declaration.** Edit `src/state/atmosphere_state.f90`, inside the type body before `end type`:

```fortran
   contains
      procedure :: init => atmosphere_state_init
```

- [ ] **Step 3: Add or relocate the subroutine.** After `end type atmosphere_state_t`, in a `contains` section at module scope:

```fortran
contains

   subroutine atmosphere_state_init(self, meteo_cfg)
      use meteorology_config_mod, only: meteorology_config_t
      class(atmosphere_state_t),  intent(inout) :: self
      type(meteorology_config_t), intent(in)    :: meteo_cfg

      ! Defaults retained from type initializers; runtime/per-year
      ! seeders update arrays. meteo_cfg arg reserved for future
      ! config-driven seed migration.
   end subroutine atmosphere_state_init

end module atmosphere_state_mod
```

- [ ] **Step 4: Update call site in swap_mod.f90.** Replace any existing `atmosphere_init(state%atmosphere, ...)` with `call state%atmosphere%init(config%meteo)`.

If no prior init existed, add the call AFTER any prior atmosphere setup (likely near the other typed init calls).

- [ ] **Step 5: Run VG.**

- [ ] **Step 6: Commit.**

```bash
git add src/state/atmosphere_state.f90 src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
refactor(gr-atm): introduce state%atmosphere%init type-bound

Mirrors heat_state / surfacewater_state pilot patterns. Signature:
state%atmosphere%init(config%meteo). config%meteo arg reserved for
future seed migration.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 15: Migrate `et.f90` body

**Files:** `src/atmosphere/et.f90` (line 614 site: `use variables, only: swredu, cofred, ...`)

- [ ] **Step 1: Identify the exact symbols imported at line 614.**

```bash
sed -n '610,625p' src/atmosphere/et.f90
```

Inspect the `use variables, only: swredu, cofred, ...` block.

- [ ] **Step 2: Replace import with explicit state/config.**

Drop the `use variables` line entirely. Replace with:

```fortran
      use swap_state_mod,  only: swap_state_t
      use swap_config_mod, only: swap_config_t
```

(If `swap_state_mod` is already imported elsewhere in the subroutine, merge.)

- [ ] **Step 3: Add `config` to signature.**

Find the subroutine declaration that contains line 614 (search backward for `subroutine` keyword). Add `config` as an `intent(in)` arg:

```fortran
   subroutine <name>(..., state, config)
      ...
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config
```

- [ ] **Step 4: Substitute symbol references in the subroutine body.**

| Was | Now |
|---|---|
| `swredu` | `config%meteo%swredu` |
| `cofred` | `config%meteo%cofred` |
| `...` (other symbols from the `only:` list) | verify per-symbol: config or state |

Run `grep -nE "\b(swredu|cofred)\b" src/atmosphere/et.f90 | head` to find all occurrences.

- [ ] **Step 5: Update callers.** Find all callers and add `config` to the call:

```bash
grep -rn "call <name>\b" src/
```

Update each `call X(...)` to `call X(..., config)`.

- [ ] **Step 6: Run VG.**

- [ ] **Step 7: Commit.**

```bash
git add src/atmosphere/et.f90 # + any updated callers
git commit -m "$(cat <<'EOF'
refactor(gr-atm): et.f90 — drop use variables

swredu/cofred → config%meteo. Signature gains config arg.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 16: Verify `et.f90` clean

**Files:** None (verification only).

- [ ] **Step 1: Grep.** Run: `grep -n "^[[:space:]]*use variables" src/atmosphere/et.f90`. **Expected:** empty (no `use variables` clause remaining).
- [ ] **Step 2: Run VG.**
- [ ] **Step 3: No commit needed.**

If any hit found, return to Task 15 to clean.

---

### Task 17: Migrate `interception.f90` — InterceptionAgric (line 40 site)

**Files:** `src/atmosphere/interception.f90`

Symbols at line 40: `gird, kdif, kdir, cofab, lai, isua`.

- [ ] **Step 1: Locate the subroutine containing line 40.**

```bash
sed -n '30,50p' src/atmosphere/interception.f90
```

- [ ] **Step 2: Replace import.**

Drop `use variables, only: gird, kdif, kdir, cofab, lai, isua`. Add:

```fortran
   use swap_state_mod,  only: swap_state_t
```

- [ ] **Step 3: Add `state` to signature** (if not already present).

```fortran
   subroutine <name>(..., state)
      ...
      type(swap_state_t), intent(inout) :: state
```

- [ ] **Step 4: Substitute references.**

| Was | Now |
|---|---|
| `gird` | `state%crop%gird` |
| `kdif` | `state%crop%kdif` |
| `kdir` | `state%crop%kdir` |
| `cofab` | `state%crop%cofab` |
| `lai` | `state%crop%lai` |
| `isua` | `state%atmosphere%isua` |

- [ ] **Step 5: Update callers.**

```bash
grep -rn "call <name>\b" src/
```

- [ ] **Step 6: Run VG.**

- [ ] **Step 7: Commit.**

```bash
git add src/atmosphere/interception.f90 # + any callers
git commit -m "refactor(gr-atm): interception.f90 InterceptionAgric — drop use variables

gird/kdif/kdir/cofab/lai → state%crop%X. isua → state%atmosphere%isua.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 18: Migrate `interception.f90` — VonHoyningenHune (line 89 site)

**Files:** `src/atmosphere/interception.f90`

Symbols at line 89: `gird, avevaptb, avprectb, pfreetb, pstemtb, scanopytb, isua`.

- [ ] **Step 1: Locate the subroutine.**

```bash
sed -n '85,100p' src/atmosphere/interception.f90
```

- [ ] **Step 2: Replace import.**

Drop `use variables, only: ...`. If state is not yet in scope, add:

```fortran
   use swap_state_mod, only: swap_state_t
```

- [ ] **Step 3: Add `state` to signature.**

- [ ] **Step 4: Substitute references.**

| Was | Now |
|---|---|
| `gird` | `state%crop%gird` |
| `avevaptb` | `state%atmosphere%avevaptb` |
| `avprectb` | `state%atmosphere%avprectb` |
| `pfreetb` | `state%atmosphere%pfreetb` |
| `pstemtb` | `state%atmosphere%pstemtb` |
| `scanopytb` | `state%atmosphere%scanopytb` |
| `isua` | `state%atmosphere%isua` |

- [ ] **Step 5: Update callers.**

- [ ] **Step 6: Run VG.**

- [ ] **Step 7: Commit.**

```bash
git add src/atmosphere/interception.f90
git commit -m "refactor(gr-atm): interception VonHoyningenHune — drop use variables

gird → state%crop. avevaptb/avprectb/pfreetb/pstemtb/scanopytb/isua →
state%atmosphere.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 19: Migrate `interception.f90` — UpdateInterception (line 171 site)

**Files:** `src/atmosphere/interception.f90`

Symbols at line 171: `logf, siccapact, fimin, ew0`.

- [ ] **Step 1: Locate the subroutine.**

```bash
sed -n '165,180p' src/atmosphere/interception.f90
```

- [ ] **Step 2: Replace import.** Retain `logf` via narrow `use variables, only: logf` (Arc 9 deferral):

```fortran
   use swap_state_mod, only: swap_state_t
   use variables,      only: logf   ! [SS-GR-ATM B6] DEFERRED to Arc 9
```

- [ ] **Step 3: Substitute references.**

| Was | Now |
|---|---|
| `siccapact` | `state%atmosphere%siccapact` |
| `fimin` | `state%atmosphere%fimin` |
| `ew0` | `state%crop%ew0` |
| `logf` | (unchanged — narrow import retained) |

- [ ] **Step 4: Run VG.**

- [ ] **Step 5: Commit.**

```bash
git add src/atmosphere/interception.f90
git commit -m "refactor(gr-atm): interception UpdateInterception — drop use variables (logf deferred)

siccapact/fimin → state%atmosphere. ew0 → state%crop. logf retained
via narrow use-only — Arc 9 deferral.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 20: Migrate `interception.f90` — InterceptionTrees (line 409 site)

**Files:** `src/atmosphere/interception.f90`

Symbols at line 409: `isua, gird, nird`.

- [ ] **Step 1: Locate the subroutine.**

```bash
sed -n '400,420p' src/atmosphere/interception.f90
```

- [ ] **Step 2: Replace import.** Retain `nird` via narrow import (Arc 8 deferral):

```fortran
   use swap_state_mod, only: swap_state_t
   use variables,      only: nird   ! [SS-GR-ATM B7] DEFERRED to Arc 8 (irrigation)
```

- [ ] **Step 3: Substitute references.**

| Was | Now |
|---|---|
| `isua` | `state%atmosphere%isua` |
| `gird` | `state%crop%gird` |
| `nird` | (unchanged — narrow import retained) |

- [ ] **Step 4: Run VG.**

- [ ] **Step 5: Commit.**

```bash
git add src/atmosphere/interception.f90
git commit -m "refactor(gr-atm): interception InterceptionTrees — drop use variables (nird deferred)

isua → state%atmosphere. gird → state%crop. nird retained via narrow
use-only — Arc 8 deferral.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 21: Verify `interception.f90` clean

**Files:** None.

- [ ] **Step 1: Grep.** Run: `grep -n "^[[:space:]]*use variables" src/atmosphere/interception.f90`.

**Expected:** at most TWO hits, both narrow imports:
- `use variables, only: logf` (line ~171, Arc 9 deferral)
- `use variables, only: nird` (line ~409, Arc 8 deferral)

- [ ] **Step 2: Run VG.**
- [ ] **Step 3: No commit.**

---

### Task 22: Migrate `meteoday.f90` — CN runoff (line 102 site)

**Files:** `src/atmosphere/meteoday.f90`

Symbols at line 102: `CNref, CNdry, CNwet, ThetaRef, Runoff_CN, wc_cor, iCNtab, CNtimTAB, CNrefTAB, wc10, ...`.

- [ ] **Step 1: Locate the subroutine.**

```bash
sed -n '95,115p' src/atmosphere/meteoday.f90
```

- [ ] **Step 2: Replace import.** Drop `use variables, only: ...`. Add:

```fortran
   use swap_state_mod, only: swap_state_t
```

- [ ] **Step 3: Add `state` to signature.**

- [ ] **Step 4: Substitute references** (all CN-method symbols → `state%atmosphere%X`):

| Was | Now |
|---|---|
| `CNref` | `state%atmosphere%CNref` |
| `CNdry` | `state%atmosphere%CNdry` |
| `CNwet` | `state%atmosphere%CNwet` |
| `ThetaRef` | `state%atmosphere%ThetaRef` |
| `Runoff_CN` | `state%atmosphere%Runoff_CN` |
| `wc_cor` | `state%atmosphere%wc_cor` |
| `wc10` | `state%atmosphere%wc10` |
| `iCNtab` | `state%atmosphere%iCNtab` |
| `CNtimTAB` | `state%atmosphere%CNtimTAB` |
| `CNrefTAB` | `state%atmosphere%CNrefTAB` |

If additional symbols appear in the `only:` list, classify each per the Phase C symbol replacement reference in the spec.

- [ ] **Step 5: Update callers.**

- [ ] **Step 6: Run VG.**

- [ ] **Step 7: Commit.**

```bash
git add src/atmosphere/meteoday.f90 # + callers
git commit -m "refactor(gr-atm): meteoday CN runoff — drop use variables

CN method symbols → state%atmosphere%X.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 23: Migrate `meteoday.f90` — meteo daily (line 289 site)

**Files:** `src/atmosphere/meteoday.f90`

Symbols at line 289 (very long `only:` list): `out_tmn, out_tmx, out_hum, out_win, out_etr, out_wet, out_rad, swrain, wet, rh, tav, tavd, arai, atmx, ahum, aetr, arad, teprrain, teprsnow, ...`.

- [ ] **Step 1: Locate the subroutine.**

```bash
sed -n '285,305p' src/atmosphere/meteoday.f90
```

- [ ] **Step 2: Replace import.**

```fortran
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
```

- [ ] **Step 3: Add `config` to signature** (state likely already there).

- [ ] **Step 4: Substitute references.**

| Was | Now |
|---|---|
| `out_tmn`, `out_tmx`, `out_hum`, `out_win`, `out_etr`, `out_wet`, `out_rad` | `state%atmosphere%out_X` |
| `arad`, `atmn`, `atmx`, `ahum`, `awin`, `arai`, `aetr` | `state%atmosphere%X` |
| `wet`, `rh`, `tav`, `tavd` | `state%atmosphere%X` (NOTE: spec uses `Tav` — verify the legacy global is `tav` or `Tav` and map to `state%atmosphere%Tav`) |
| `teprrain`, `teprsnow` | `state%atmosphere%X` |
| `swrain` | `config%meteo%swrain` |

- [ ] **Step 5: Update callers.**

- [ ] **Step 6: Run VG.**

- [ ] **Step 7: Commit.**

```bash
git add src/atmosphere/meteoday.f90
git commit -m "refactor(gr-atm): meteoday daily processing — drop use variables

Daily arrays + derived scalars + output flags → state%atmosphere.
swrain → config%meteo.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 24: Migrate `meteoday.f90` — ET/interception driver (line 543 site)

**Files:** `src/atmosphere/meteoday.f90`

Symbols at line 543: `lai, gird, swinter, swmetdetail, nmetdetail, swetr, flCropEmergence, et0, ew0, es0, swcf, swcfbs, cfbs, ...`.

- [ ] **Step 1: Locate the subroutine.**

```bash
sed -n '540,560p' src/atmosphere/meteoday.f90
```

- [ ] **Step 2: Replace import.**

```fortran
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
```

- [ ] **Step 3: Substitute references.**

| Was | Now |
|---|---|
| `lai`, `gird`, `swcf`, `swcfbs`, `cfbs`, `flCropEmergence`, `et0`, `ew0`, `es0` | `state%crop%X` |
| `swinter`, `swmetdetail`, `nmetdetail`, `swetr` | `config%meteo%X` |

- [ ] **Step 4: Update callers.**

- [ ] **Step 5: Run VG.**

- [ ] **Step 6: Commit.**

```bash
git add src/atmosphere/meteoday.f90
git commit -m "refactor(gr-atm): meteoday ET/interception driver — drop use variables

Crop refs → state%crop. Meteo config switches → config%meteo.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 25: Verify `meteoday.f90` clean

**Files:** None.

- [ ] **Step 1: Grep.** `grep -n "^[[:space:]]*use variables" src/atmosphere/meteoday.f90`. **Expected:** empty (or only documented deferrals).
- [ ] **Step 2: Run VG.**
- [ ] **Step 3: No commit.**

---

### Task 26: Migrate `meteodt.f90` — bare line 58

**Files:** `src/atmosphere/meteodt.f90`

Line 58 has bare `use variables`.

- [ ] **Step 1: Inspect what symbols are actually used.**

Find the subroutine starting around line 58, then grep all symbols read in its body:

```bash
sed -n '55,135p' src/atmosphere/meteodt.f90
```

Identify which globals are read in the body (could be 10–30 symbols).

- [ ] **Step 2: Convert bare to narrow imports.**

Replace `use variables` with explicit imports. For each read symbol, choose:
- **State**: `state%atmosphere%X` or `state%crop%X` — drop from imports.
- **Config**: `config%meteo%X` etc. — drop from imports.
- **Deferred**: retain via `use variables, only: X` with comment.

Example final state:

```fortran
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
```

- [ ] **Step 3: Add `state` and/or `config` to signature.**

- [ ] **Step 4: Substitute symbol references** per the spec's symbol replacement table.

- [ ] **Step 5: Update callers.**

- [ ] **Step 6: Run VG.**

- [ ] **Step 7: Commit.**

```bash
git add src/atmosphere/meteodt.f90 # + callers
git commit -m "refactor(gr-atm): meteodt line 58 — drop bare use variables

Symbols routed to state%atmosphere/state%crop or config%meteo per spec
replacement table. <document any deferrals>.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 27: Migrate `meteodt.f90` — rain timing (line 139 site)

**Files:** `src/atmosphere/meteodt.f90`

Symbols at line 139: `swrain, raintab, wet, nmrain, rainamount, arai, rainfluxarray, raintimearray`.

- [ ] **Step 1: Replace import.**

```fortran
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
```

- [ ] **Step 2: Verify config paths.** Run:

```bash
grep -nE "swrain|raintab|rainamount|rainfluxarray|raintimearray|nmrain" src/config/meteorology_config.f90
```

`raintab`/`rainamount`/`rainfluxarray`/`raintimearray`/`nmrain` may be in `config%meteo` (input data), or they may be runtime state. Classify per-symbol.

- [ ] **Step 3: Substitute references.**

| Was | Now |
|---|---|
| `swrain` | `config%meteo%swrain` |
| `raintab` | `config%meteo%raintab` (input table) |
| `wet` | `state%atmosphere%wet` |
| `arai` | `state%atmosphere%arai` |
| `nmrain` | `state%atmosphere%nmrain` OR `config%meteo%nmrain` — verify (rainfall record count; likely state if computed at runtime) |
| `rainamount`, `rainfluxarray`, `raintimearray` | `state%atmosphere%X` OR `config%meteo%X` — verify per-symbol; defer narrow if not yet in state schema |

**Note:** if `nmrain`/`rainamount`/etc. are NOT in atmosphere_state (they may have been deferred from Block 1/2 scope), retain via narrow `use variables, only: nmrain, rainamount, ...` with `[SS-GR-ATM B14] DEFERRED — rain timing arrays not yet on state`. Sub-arc cleanup possible later.

- [ ] **Step 4: Update callers.**

- [ ] **Step 5: Run VG.**

- [ ] **Step 6: Commit.**

```bash
git add src/atmosphere/meteodt.f90
git commit -m "refactor(gr-atm): meteodt rain timing — drop use variables

swrain/raintab → config%meteo. wet/arai → state%atmosphere. Rain timing
arrays → state or deferred (document).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 28: Migrate `meteodt.f90` — bare lines 358 + 456

**Files:** `src/atmosphere/meteodt.f90`

Both bare `use variables` sites.

- [ ] **Step 1: Per-site inventory.** For each site:

```bash
sed -n '355,400p' src/atmosphere/meteodt.f90  # for line 358 site
sed -n '450,480p' src/atmosphere/meteodt.f90  # for line 456 site
```

Identify symbols read in each subroutine body.

- [ ] **Step 2: Apply substitutions** per spec replacement table.

- [ ] **Step 3: Update callers.**

- [ ] **Step 4: Run VG.**

- [ ] **Step 5: Commit.**

```bash
git add src/atmosphere/meteodt.f90 # + callers
git commit -m "refactor(gr-atm): meteodt lines 358 + 456 — drop bare use variables

Substitutions per spec replacement table; deferrals documented.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 29: Verify `meteodt.f90` clean

**Files:** None.

- [ ] **Step 1: Grep.** `grep -n "^[[:space:]]*use variables" src/atmosphere/meteodt.f90`. **Expected:** empty or only documented deferrals.
- [ ] **Step 2: Run VG.**
- [ ] **Step 3: No commit.**

---

### Task 30: Phase B close — check-full 5/5

**Files:** None.

- [ ] **Step 1: Clean rebuild + check-full.** `rm -rf builddir && pixi run build-linux && pixi run check-full` — 5/5 byte-for-byte.
- [ ] **Step 2: Sanity greps.**
  - `grep -rn "^[[:space:]]*use variables" src/atmosphere/*.f90 | grep -v "logf\|nird" | head` — empty or only documented deferrals.
  - `grep -n "state%atmosphere%init\b" src/core/swap_mod.f90` ≥1.
- [ ] **Step 3: Phase B close marker.**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-atm): Phase B complete — 4-file reader cutover + atmosphere init

et/interception/meteoday/meteodt migrated off bare use variables.
state%atmosphere%init type-bound (per heat_state pilot). check-full 5/5
byte-for-byte.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Phase C — Codebase Sweep + Global Retirement (Tasks 31–46)

**Scope discipline:** Phase C does NOT touch any other `use variables` symbols. Files keep their `use variables` clause for unrelated globals. Surgical replacement only.

**Contingency:** if a file reads atmosphere/crop globals but doesn't already take `state`, add `type(swap_state_t), intent(inout) :: state` (or `intent(in)` if read-only).

### Symbol replacement reference (Phase C)

| Was | Now |
|---|---|
| `arad/atmn/atmx/ahum/awin/arai/aetr/wet` | `state%atmosphere%X` |
| `atav/epot/tpot/grain/nrain` | `state%atmosphere%X` |
| `Tav/tavd/rh/daynrfirst/daynrlast/atmin7/nofd/teprrain/teprsnow` | `state%atmosphere%X` |
| `siccapact/fimin/isua/avevaptb/avprectb/pfreetb/pstemtb/scanopytb` | `state%atmosphere%X` |
| `CNref/CNdry/CNwet/ThetaRef/Runoff_CN/wc_cor/wc10/iCNtab/CNtimTAB/CNrefTAB` | `state%atmosphere%X` |
| `out_tmn/out_tmx/out_hum/out_win/out_etr/out_wet/out_rad` | `state%atmosphere%X` |
| `lai/kdif/kdir/cofab/swcf/swcfbs/cfbs/gird/flCropEmergence/et0/ew0/es0` | `state%crop%X` |

---

### Task 31: Cluster `src/core/`

**Files:** `src/core/swap_mod.f90`, `initialize.f90`, `timecontrol_mod.f90`, `swap_bmi_mod.f90`, `swap_capi_mod.f90`.

- [ ] **Step 1: Per-file inventory.** For each, grep for atmosphere + crop globals:

```bash
for f in src/core/{swap_mod,initialize,timecontrol_mod,swap_bmi_mod,swap_capi_mod}.f90; do
  echo "=== $f ==="
  grep -nE "\b(arad|atmn|atmx|ahum|awin|arai|aetr|wet|atav|Tav|tavd|rh|daynrfirst|daynrlast|atmin7|nofd|teprrain|teprsnow|siccapact|fimin|isua|avevaptb|avprectb|pfreetb|pstemtb|scanopytb|CNref|CNdry|CNwet|ThetaRef|Runoff_CN|wc_cor|wc10|iCNtab|CNtimTAB|CNrefTAB|out_tmn|out_tmx|out_hum|out_win|out_etr|out_wet|out_rad|lai|kdif|kdir|cofab|swcf|swcfbs|cfbs|gird|flCropEmergence|et0|ew0|es0)\b" "$f" | grep -v "state%\|config%\|! \|character\|integer ::\|real(8) ::\|real(real64) ::" | head -40
done
```

- [ ] **Step 2: Apply substitutions** per the reference table.

**swap_mod.f90 caveat:** the Phase A dual-write blocks (`[SS-GR-ATM A10-A12]`) still need the legacy globals as their source. Those imports STAY until Phase C C11+C12 (when readers route directly into state). Leave them for now.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/core/swap_mod.f90 src/core/initialize.f90 src/core/timecontrol_mod.f90 src/core/swap_bmi_mod.f90 src/core/swap_capi_mod.f90
git commit -m "sweep(gr-atm): src/core/ — atmosphere + crop refs to state%X

Atmosphere data → state%atmosphere. Crop runtime → state%crop. Files
retain use variables for unrelated globals (and dual-write imports in
swap_mod, retired in Phase C C11-C13).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 32: Cluster `src/io/`

**Files:** `src/io/swap_csv_output.f90`, `src/io/swapoutput.f90`, `src/io/readmeteo.f90`.

- [ ] **Step 1: Per-file inventory** (same pattern as Task 31).

- [ ] **Step 2: Apply substitutions.** Special caveats:

- `swap_csv_output.f90` and `swapoutput.f90` read OUTPUT flags (out_X) heavily — these become `state%atmosphere%out_X`.
- Daily meteo arrays (arad/atmn/etc.) are also output — `state%atmosphere%X`.
- Tav/atav specifically: read for OUTTEM output buffer.

- `readmeteo.f90` is the per-year meteo loader. It WRITES the legacy globals. After this sweep, it writes to `state%atmosphere%X` directly (handled per readmeteo subroutine; treat as a dual-write retirement). Note: this is partial migration in Task 32; the full readmeteo retirement happens in Task 41 (C11).

For Task 32: migrate READS only in readmeteo. Leave WRITES (the populator paths) for Task 41 to ensure readmeteo continues populating legacy globals during the cluster sweep period.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/io/swap_csv_output.f90 src/io/swapoutput.f90 src/io/readmeteo.f90
git commit -m "sweep(gr-atm): src/io/ — meteo arrays + output flags to state%X

CSV output and swapoutput modules read state%atmosphere directly. readmeteo
READS migrated; writes deferred to Task 41 (C11) global-retirement bundle.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 33: Cluster `src/boundary/` (close prior-arc deferrals)

**Files:** `src/boundary/boundbottom.f90`, `src/boundary/boundtop.f90`.

These files have prior-arc deferrals (from GR-BH Tasks 18 + 20):
- `boundbottom.f90`: deferred `logf, SwBotb3ResVert, gwltab, qbotab, haqtab, hbotab` (none atmosphere-related)
- `boundtop.f90`: deferred `nird, swkmean, swredu, flrunon, runonarr`

This task retires the atmosphere-related deferrals: `swkmean` (already config%simulation), `swredu` (config%meteo), `flrunon` and `runonarr` (config — verified in Task 9), and `Tav` if read (closes GR-BH Tav/atav deferral target).

- [ ] **Step 1: Inventory deferred narrow imports.**

```bash
grep -n "use variables" src/boundary/boundbottom.f90 src/boundary/boundtop.f90
```

- [ ] **Step 2: For boundtop:** narrow the `use variables, only:` line. Remove migrated symbols (`swkmean, swredu, flrunon, runonarr`) — they're now imported from config. Body substitutions:

| Was | Now |
|---|---|
| `swkmean` | `config%simulation%swkmean` |
| `swredu` | `config%meteo%swredu` |
| `flrunon` | `config%simulation%flrunon` (or wherever Task 9 placed it) |
| `runonarr` | `config%simulation%runonarr` |

**If config is not in scope** at the boundtop call (per GR-BH Task 20 caveat — boundtop is called from soilhydraulics, not swap_mod), thread config through OR keep as deferred — preserve the GR-BH pragmatic approach if needed. The decision: if threading config requires touching soilhydraulics, defer further (Arc 5 territory). If it doesn't, migrate now.

- [ ] **Step 3: For boundbottom:** check if any atmosphere/crop refs exist. The deferred imports were all bottom-boundary-config-related — not atmosphere. Likely no changes here.

- [ ] **Step 4: Run VG.**

- [ ] **Step 5: Commit.**

```bash
git add src/boundary/boundbottom.f90 src/boundary/boundtop.f90
git commit -m "sweep(gr-atm): src/boundary/ — close GR-BH atmosphere deferrals

boundtop's swkmean/swredu/flrunon/runonarr (Arc 4 deferrals from GR-BH
Task 20) now read from config%simulation / config%meteo. boundbottom
unchanged (its deferrals are not atmosphere-related).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 34: Cluster `src/heat/` (close Tav/atav deferral)

**Files:** `src/heat/temperature.f90`

GR-BH Task 15 left `use variables, only: Tav, atav` as a deferred narrow import at line 83. This task closes that deferral.

- [ ] **Step 1: Inspect current state.**

```bash
grep -n "use variables\|Tav\|atav" src/heat/temperature.f90 | head
```

- [ ] **Step 2: Replace the narrow import.**

Drop the line `use variables, only: Tav, atav`. Substitute body references:

| Was | Now |
|---|---|
| `Tav` | `state%atmosphere%Tav` |
| `atav(...)` | `state%atmosphere%atav(...)` |

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/heat/temperature.f90
git commit -m "sweep(gr-atm): src/heat/temperature.f90 — close Tav/atav deferral

Closes GR-BH Task 15 deferral. Tav/atav → state%atmosphere%X.
temperature.f90 now fully clean of bare use variables.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 35: Cluster `src/soil/`

**Files:** `src/soil/soilgrid.f90`, `src/soil/soilhydraulics.f90`, `src/soil/waterbalance.f90`.

- [ ] **Step 1: Per-file inventory** for atmosphere + crop globals (same grep as Task 31).

- [ ] **Step 2: Apply substitutions.** Likely minimal exposure — soil reads vegetation/meteo only for special cases.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/soil/soilgrid.f90 src/soil/soilhydraulics.f90 src/soil/waterbalance.f90
git commit -m "sweep(gr-atm): src/soil/ — atmosphere + crop refs to state%X

Soil cluster reads atmosphere/crop data via state%X. Files retain use
variables for unrelated soil/drainage globals.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 36: Cluster `src/drainage/`

**Files:** `src/drainage/drainage.f90`, `src/drainage/surfacewater.f90`, `src/drainage/divdra.f90`.

- [ ] **Step 1: Per-file inventory.**
- [ ] **Step 2: Apply substitutions** — likely minimal exposure to atmosphere/crop.
- [ ] **Step 3: Run VG.**
- [ ] **Step 4: Commit.**

```bash
git add src/drainage/drainage.f90 src/drainage/surfacewater.f90 src/drainage/divdra.f90
git commit -m "sweep(gr-atm): src/drainage/ — atmosphere + crop refs to state%X

Drainage cluster reads atmosphere/crop via state%X.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 37: Cluster `src/crop/`

**Files:** All 9 crop files.

This cluster has the MOST atmosphere/crop refs since crop code BOTH reads and writes its own crop_state data, AND reads atmosphere data.

- [ ] **Step 1: Per-file inventory.**

```bash
for f in src/crop/*.f90; do
  echo "=== $f ==="
  grep -nE "\b(arad|atmn|atmx|ahum|awin|arai|aetr|wet|Tav|tavd|rh|lai|kdif|kdir|cofab|swcf|swcfbs|cfbs|gird|flCropEmergence|et0|ew0|es0)\b" "$f" | grep -v "state%\|config%\|! \|character\|integer ::\|real(8) ::" | head -30
done
```

- [ ] **Step 2: Apply substitutions.** Crop code references its OWN data — but crop code may need to UPDATE these fields. Watch for writes:

```bash
grep -nE "\b(lai|kdif|cofab|gird|et0|ew0|es0|swcf|swcfbs)\s*=" src/crop/*.f90 | head
```

For writes: update both legacy global AND `state%crop%X` (transitional dual-write) during the cluster sweep. Phase C global retirement (Task 45) drops the legacy writes.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/crop/*.f90
git commit -m "sweep(gr-atm): src/crop/ — atmosphere + crop self-refs to state%X

Crop cluster reads atmosphere via state%atmosphere; crop's own runtime
data routed through state%crop. Runtime writers dual-write to legacy
globals + state%crop until Task 45 (crop globals deletion).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 38: Cluster `src/solute/`

**Files:** `src/solute/solute.f90`, `src/solute/agetracer.f90`.

- [ ] **Step 1: Per-file inventory.**
- [ ] **Step 2: Apply substitutions** — solute reads rain/temp for partition + tracer.
- [ ] **Step 3: Run VG.**
- [ ] **Step 4: Commit.**

```bash
git add src/solute/solute.f90 src/solute/agetracer.f90
git commit -m "sweep(gr-atm): src/solute/ — atmosphere refs to state%X

solute/agetracer read rain + temperature via state%atmosphere.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 39: Self-audit `src/atmosphere/`

**Files:** `src/atmosphere/et.f90`, `src/atmosphere/interception.f90`, `src/atmosphere/meteoday.f90`, `src/atmosphere/meteodt.f90`.

These were migrated in Phase B but may still have stale refs that the cluster sweep needs to verify.

- [ ] **Step 1: Grep all atmosphere/crop bare refs.**

```bash
grep -rnE "\b(arad|atmn|atmx|ahum|awin|arai|aetr|wet|atav|Tav|tavd|rh|daynrfirst|daynrlast|atmin7|nofd|teprrain|teprsnow|siccapact|fimin|isua|avevaptb|avprectb|pfreetb|pstemtb|scanopytb|CNref|CNdry|CNwet|ThetaRef|Runoff_CN|wc_cor|wc10|iCNtab|CNtimTAB|CNrefTAB|out_tmn|out_tmx|out_hum|out_win|out_etr|out_wet|out_rad|lai|kdif|kdir|cofab|swcf|swcfbs|cfbs|gird|flCropEmergence|et0|ew0|es0)\b" src/atmosphere/ --include="*.f90" | grep -v "state%\|config%\|! \|character\|integer ::\|real(8) ::"
```

- [ ] **Step 2: Patch any stale refs found.**
- [ ] **Step 3: Run VG.**
- [ ] **Step 4: Commit (if any patches).**

```bash
git add src/atmosphere/*.f90
git commit -m "audit(gr-atm): atmosphere self-audit — patch any stale refs

Post-cluster-sweep verification. <document count> stale refs patched.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

If audit clean: no commit needed.

---

### Task 40: Audit pass — find stale reads

**Files:** Audit + targeted fixes. MANDATORY GATE before global deletion.

- [ ] **Step 1: Daily-meteo audit.**

```bash
grep -rnE "\b(arad|atmn|atmx|ahum|awin|arai|aetr|wet|atav|Tav|tavd|rh)\b" src/ --include="*.f90" \
  | grep -v "state%\|config%\|! \|variables.f90\|use variables\|integer ::\|real(8) ::\|real(real64) ::\|character\|^src/core/initialize\|^src/io/toml/config_to_variables\|^src/io/readmeteo\|deferred" | head -40
```

- [ ] **Step 2: Interception + CN audit.**

```bash
grep -rnE "\b(siccapact|fimin|isua|avevaptb|avprectb|pfreetb|pstemtb|scanopytb|CNref|CNdry|CNwet|ThetaRef|Runoff_CN|wc_cor|wc10|iCNtab|CNtimTAB|CNrefTAB)\b" src/ --include="*.f90" \
  | grep -v "state%\|config%\|! \|variables.f90\|use variables\|integer ::\|real(8) ::\|real(real64) ::\|character\|^src/core/initialize\|^src/io/toml/config_to_variables\|deferred" | head -40
```

- [ ] **Step 3: Output flags audit.**

```bash
grep -rnE "\b(out_tmn|out_tmx|out_hum|out_win|out_etr|out_wet|out_rad)\b" src/ --include="*.f90" \
  | grep -v "state%\|config%\|! \|variables.f90\|use variables\|logical :: " | head -20
```

- [ ] **Step 4: Crop refs audit.**

```bash
grep -rnE "\b(lai|kdif|kdir|cofab|swcf|swcfbs|cfbs|gird|flCropEmergence|et0|ew0|es0)\b" src/ --include="*.f90" \
  | grep -v "state%\|config%\|! \|variables.f90\|use variables\|integer ::\|real(8) ::\|real(real64) ::\|logical :: \|character\|^src/core/initialize\|^src/io/toml/config_to_variables\|deferred" | head -40
```

- [ ] **Step 5: Patch any stale refs in separate commits.** For each:

```bash
git add <file>
git commit -m "audit(gr-atm): <file> — patch stale <symbol> reference

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

- [ ] **Step 6: Re-run all 4 audits.** Each must return empty (or only documented false positives).

- [ ] **Step 7: Audit close marker.**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
audit(gr-atm): codebase sweep clean — no stale atmosphere/crop reads

All 4 audit greps clean. Ready for Phase C global deletion (Tasks 41-45).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

**This task is MANDATORY GATE.** Do NOT proceed to Task 41 if audits aren't clean.

---

### Task 41: Drop legacy meteo writes from adapter + readmeteo

**Files:** `src/io/toml/config_to_variables.f90`, `src/io/readmeteo.f90`, `src/core/swap_mod.f90`.

Atmosphere globals are populated by TWO paths:
1. `config_to_variables` — initial seeding (config tables, switches).
2. `readmeteo` — per-year meteo reload (daily/sub-daily arrays).

Both must route into `state%atmosphere` directly. The Phase A dual-writes can then be retired.

- [ ] **Step 1: Inspect adapter writes.**

```bash
grep -nE "^[[:space:]]*\b(arad|atmn|atmx|ahum|awin|arai|aetr|wet|atav|Tav|tavd|rh|daynrfirst|daynrlast|atmin7|nofd|teprrain|teprsnow|siccapact|fimin|isua|avevaptb|avprectb|pfreetb|pstemtb|scanopytb|CNref|CNdry|CNwet|ThetaRef|Runoff_CN|wc_cor|wc10|iCNtab|CNtimTAB|CNrefTAB|out_tmn|out_tmx|out_hum|out_win|out_etr|out_wet|out_rad)\b\s*(\(|=)" src/io/toml/config_to_variables.f90 | head -40
```

- [ ] **Step 2: For each adapter write:**
  - Either rewrite the dual-write in `swap_mod.f90` (Phase A blocks A10-A12) to source from `config%X` (where the value comes from config), then delete the legacy adapter write.
  - Or leave the adapter write if the value isn't directly from config (some interception tables are derived); delete only in Task 44 along with the declarations.

- [ ] **Step 3: Inspect readmeteo writes.**

```bash
grep -nE "^[[:space:]]*\b(arad|atmn|atmx|ahum|awin|arai|aetr|wet|atav|Tav|tavd|rh)\b\s*\(?[^=)]*\)?\s*=" src/io/readmeteo.f90 | head -30
```

`readmeteo` likely populates daily arrays per-year. Update to write to `state%atmosphere%X` instead of bare globals.

**Caveat:** if readmeteo doesn't have `state` in scope, add it (signature already updated in Task 32). After this, the dual-write in `swap_mod.f90 [SS-GR-ATM A10]` becomes redundant for arad/atmn/etc — readmeteo's direct writes are authoritative. Retire the A10 block.

- [ ] **Step 4: Run VG.**

- [ ] **Step 5: Commit.**

```bash
git add src/io/toml/config_to_variables.f90 src/io/readmeteo.f90 src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
retire(gr-atm): route adapter + readmeteo writes to state%atmosphere

config_to_variables and readmeteo write meteo data DIRECTLY into
state%atmosphere. Phase A dual-write blocks A10-A11 retired where
direct-write replaces them. Sets up Task 44 declaration deletion.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 42: Drop crop-state legacy writes

**Files:** `src/io/toml/config_to_variables.f90`, `src/core/swap_mod.f90`, possibly `src/crop/cropwofost_init.f90` or crop runtime files.

- [ ] **Step 1: Inspect crop writes in adapter.**

```bash
grep -nE "^[[:space:]]*\b(lai|kdif|kdir|cofab|swcf|swcfbs|cfbs|gird|flCropEmergence|et0|ew0|es0)\b\s*(\(|=)" src/io/toml/config_to_variables.f90 src/crop/*.f90 | head -30
```

- [ ] **Step 2: For each write:**
  - Adapter writes: rewrite the Phase A dual-write block A12 to source from `config%crop%X` where possible, then delete the legacy adapter write.
  - Crop-runtime writes: dual-write to legacy + state was added in Task 37. Drop the legacy write; keep state%crop write.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/io/toml/config_to_variables.f90 src/core/swap_mod.f90 src/crop/*.f90
git commit -m "retire(gr-atm): drop legacy crop writes from adapter + runtime

Adapter routes config-sourced crop values directly into state%crop.
Crop runtime writers drop legacy global writes; state%crop is sole
home. Phase A dual-write block A12 reduced.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 43: Drop output-flag + CN method legacy writes

**Files:** `src/io/toml/config_to_variables.f90`, `src/core/swap_mod.f90`.

- [ ] **Step 1: Inspect.**

```bash
grep -nE "^[[:space:]]*\b(out_tmn|out_tmx|out_hum|out_win|out_etr|out_wet|out_rad|CNref|CNdry|CNwet|ThetaRef|Runoff_CN|wc_cor|wc10|iCNtab|CNtimTAB|CNrefTAB)\b\s*(\(|=)" src/io/toml/config_to_variables.f90 | head
```

- [ ] **Step 2: Route output flags directly to state.**

For each `out_X = config%X` line in adapter, change to `state%atmosphere%out_X = config%output_csv%X` (or wherever output flags live in config). Delete legacy global write.

Same pattern for CN tables.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/io/toml/config_to_variables.f90 src/core/swap_mod.f90
git commit -m "retire(gr-atm): drop legacy output-flag + CN writes

Output flags and CN tables routed directly into state%atmosphere.
Phase A dual-write blocks A11-A12 reduced further.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 44: Delete meteo + interception + CN globals from `variables.f90`

**Files:** `src/core/variables.f90`, `src/core/initialize.f90`, plus any files surfaced by compile errors.

- [ ] **Step 1: Locate declarations.**

```bash
grep -nE "^[[:space:]]*(real\(8\)|integer|logical|character)\s+(arad|atmn|atmx|ahum|awin|arai|aetr|wet|atav|epot|tpot|grain|nrain|Tav|tavd|rh|daynrfirst|daynrlast|atmin7|nofd|teprrain|teprsnow|siccapact|fimin|isua|avevaptb|avprectb|pfreetb|pstemtb|scanopytb|CNref|CNdry|CNwet|ThetaRef|Runoff_CN|wc_cor|wc10|iCNtab|CNtimTAB|CNrefTAB|out_tmn|out_tmx|out_hum|out_win|out_etr|out_wet|out_rad|cofred)\b" src/core/variables.f90
```

- [ ] **Step 2: Replace each declaration line with a tombstone comment** (following GR-BH Task 35 pattern):

```fortran
   ! [SS-GR-ATM Task 44] arad/atmn/atmx/ahum/awin/arai/aetr/wet retired — state%atmosphere
   ! [SS-GR-ATM Task 44] atav/epot/tpot/grain/nrain retired — state%atmosphere
   ! ... etc per logical grouping ...
```

- [ ] **Step 3: Drop zero-fills in initialize.f90.**

```bash
grep -nE "^[[:space:]]*(arad|atmn|atmx|ahum|awin|arai|aetr|wet|atav|epot|tpot|grain|nrain|Tav|tavd|rh|atmin7|teprrain|teprsnow|siccapact|fimin|isua|CNref|CNdry|CNwet|ThetaRef|Runoff_CN|wc_cor|wc10|iCNtab|CNtimTAB|CNrefTAB|avevaptb|avprectb|pfreetb|pstemtb|scanopytb|out_tmn|out_tmx|out_hum|out_win|out_etr|out_wet|out_rad)\b\s*=" src/core/initialize.f90
```

Delete each zero-fill line.

- [ ] **Step 4: Iterative build cycles.**

```bash
rm -rf builddir && pixi run build-linux 2>&1 | grep -E "Error|undefined" | head -20
```

For each compile error, patch the file and retry. Common patterns:
- "Symbol 'arad' has no IMPLICIT type" — file still reads bare global; migrate to `state%atmosphere%arad`.

- [ ] **Step 5: Run VG when clean.**

- [ ] **Step 6: Commit.**

```bash
git add -A
git commit -m "$(cat <<'EOF'
retire(gr-atm): delete atmosphere globals from variables.f90

Mass deletion of meteo arrays + interception state + CN method + output
flags + derived scalars (~45 globals). state%atmosphere is the sole
home. Compile-surfaced readers patched in <list files>.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 45: Delete crop globals from `variables.f90`

**Files:** `src/core/variables.f90`, `src/core/initialize.f90`, plus compile-surfaced patches.

- [ ] **Step 1: Locate declarations.**

```bash
grep -nE "^[[:space:]]*(real\(8\)|integer|logical|character)\s+(lai|kdif|kdir|cofab|swcf|swcfbs|cfbs|gird|flCropEmergence|et0|ew0|es0)\b" src/core/variables.f90
```

- [ ] **Step 2: Replace each declaration with tombstone.**

- [ ] **Step 3: Drop zero-fills in initialize.f90.**

- [ ] **Step 4: Iterative build.** Common patches likely in `src/crop/` since crop code may still have stragglers.

- [ ] **Step 5: Run VG.**

- [ ] **Step 6: Commit.**

```bash
git add -A
git commit -m "retire(gr-atm): delete crop globals from variables.f90

12 crop runtime globals (lai/kdif/kdir/cofab/swcf/swcfbs/cfbs/gird/
flCropEmergence/et0/ew0/es0) deleted. state%crop is the sole home.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 46: Final verification + arc-complete marker

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

- [ ] **Step 3: pFUnit.**

```bash
pixi run test-pfunit
```

- [ ] **Step 4: Final sanity greps.**

```bash
# atmosphere globals retired
grep -nE "^[[:space:]]*real\(8\).*\b(arad|atmn|atmx|ahum|awin|arai|aetr|wet|atav|Tav|tavd|rh)\b" src/core/variables.f90
# Expected: empty.

# crop globals retired
grep -nE "^[[:space:]]*real\(8\).*\b(lai|kdif|gird|et0|ew0|es0)\b" src/core/variables.f90
# Expected: empty.

# atmosphere files clean (only documented deferrals)
grep -rn "^[[:space:]]*use variables" src/atmosphere/*.f90 | grep -v "logf\|nird"
# Expected: empty.

# state usage spread
grep -rl "state%atmosphere%" src/ | wc -l
# Expected: ≥ 25

grep -rl "state%crop%" src/ | wc -l
# Expected: ≥ 10
```

- [ ] **Step 5: Arc-complete marker.**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-atm): GR-ATM complete — atmosphere globals retired, crop_state introduced

check-full 5/5 byte-for-byte. BMI + cffi-demo passing. pFUnit all-pass.

Net retirement:
- 4 readers (et/interception/meteoday/meteodt) migrated off bare use
  variables (modulo documented narrow-only deferrals: logf → Arc 9,
  nird → Arc 8, plus any rain-timing arrays deferred sub-arc)
- NEW state%crop subrecord (12 fields) — foundation for Arc 8
- state%atmosphere extended with ~50 fields across 6 blocks
- state%atmosphere%init type-bound (per heat_state pilot)
- config%meteo extended with cofred
- ~45 atmosphere globals retired → state%atmosphere
- 12 crop globals retired → state%crop
- ~30 codebase readers swept to state references

Closes prior-arc deferrals:
- Tav/atav (GR-BH temperature.f90)
- swkmean/swredu/flrunon/runonarr (GR-BH boundtop)
- cofred (et.f90 long-standing)

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Plan Self-Review

**Spec coverage:**
- ✅ crop_state_mod introduction → Task 2
- ✅ atmosphere_state extensions (6 blocks) → Tasks 3-8
- ✅ config%meteo + cofred + flrunon/runonarr verification → Task 9
- ✅ Dual-writes → Tasks 10-12
- ✅ Phase A close → Task 13
- ✅ atmosphere init relocation → Task 14
- ✅ 4-file reader migration → Tasks 15-29
- ✅ Phase B close → Task 30
- ✅ Codebase cluster sweeps → Tasks 31-39
- ✅ Audit pass → Task 40
- ✅ Adapter retirement → Tasks 41-43
- ✅ Global declaration deletion → Tasks 44-45
- ✅ Final verification → Task 46

**Placeholder scan:** No "TBD" / "implement later" / "similar to Task N". Code blocks present in every code step. Inline `<document outcome>` markers are commit-message customization points, not plan placeholders — implementer fills them in based on actual findings.

**Type consistency:** crop_state_t field names consistent across Task 2 + Task 12 dual-write + Phase C tasks. atmosphere_state_t block names (Block 1-6) consistent in spec + plan tasks.

**Known soft spots flagged inline:**
- Several config field names need verification before substitution (e.g., flrunon/runonarr placement, swrain in config%meteo vs elsewhere).
- Rain timing arrays (`nmrain/rainamount/rainfluxarray/raintimearray`) — schema scope decision deferred to implementer at Task 27.
- readmeteo per-year reload retirement bundled in Task 41 (C11).
- Output flags read-by/written-by path verified in Task 9.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-05-13-globals-atmosphere.md`. Two execution options:

**1. Subagent-Driven (recommended)** — Fresh subagent per task, two-stage review (spec compliance + code quality) between each. Same pattern as GR-UTILS + GR-BH.

**2. Inline Execution** — Execute in this session via executing-plans, batch execution with checkpoints. Higher main-context pressure.

Which approach?
