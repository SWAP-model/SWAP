---
title: "GR-ATM — Globals Retirement Arc 4 (Atmosphere) with crop_state introduction"
date: 2026-05-13
status: approved
context: globals-retirement Arc 4 (atmosphere) — full scope, follows GR-BH
---

# GR-ATM Design

## Goal

Retire `use variables` entirely from `et.f90`, `interception.f90`, `meteoday.f90`, `meteodt.f90` (4 atmosphere files) by:

1. Introducing `state%crop` (NEW `crop_state_t`) with the ~12 crop-runtime symbols read by atmosphere code. Lays foundation for Arc 8 (crop).
2. Extending `state%atmosphere` with ~50 new fields across 6 blocks: daily meteo arrays, sub-daily meteo arrays, derived meteo scalars, interception state/params, Runoff-CN method state, output flag toggles.
3. Adding `cofred` (and verifying `flrunon`/`runonarr` placement) in `config%meteo` / `config%simulation`.
4. Type-bound `state%atmosphere%init(config%meteo)` per surfacewater_state pilot.

**Closes prior-arc deferrals:** `Tav, atav` (temperature.f90 from GR-BH), `swkmean, swredu, flrunon, runonarr` (boundtop from GR-BH), `cofred` (et.f90 read).

## Architecture

Three sub-phases inside one arc; each phase ends with `pixi run check-full` 5/5 byte-for-byte.

- **Phase A** — Schema + dual-write. Additive only. New `crop_state_mod`; atmosphere_state extended in 6 blocks (commits separated for review-ability); `config%meteo` extended with `cofred`. State mirrors populated alongside legacy globals.
- **Phase B** — Reader cutover + state%atmosphere%init introduction. The 4 target files drop `use variables` and read state/config exclusively.
- **Phase C** — Codebase mesh sweep + global retirement. Every other reader (~30 files) of the migrated globals swaps to `state%X`. Then bare globals are deleted from `variables.f90`.

## Verification Discipline

Per `feedback_per_task_regression_gate.md` + `feedback_state_schema_clean_rebuild.md`:

- **Every implementer subagent ends with VG:**
  ```bash
  rm -rf builddir && pixi run build-linux \
    && pixi run test-pfunit \
    && pixi run -e test python tests/regression/test_output_regression.py \
         hupselbrook surfacewater salinitystress grassgrowth
  ```
  Expected: clean build, pFUnit all-pass, regression **4/4 byte-for-byte**.
- **Every phase ends with:** `pixi run check-full` 5/5.
- **Each commit migrates one file / one block / one cohort.**

## Schema Additions

### NEW: `src/state/crop_state.f90`

```fortran
module crop_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: crop_state_t

   type :: crop_state_t
      ! Crop runtime — fields read by atmosphere (Arc 4 boot scope)
      real(real64) :: lai            = 0.0_real64    !! leaf area index (-)
      real(real64) :: kdif           = 0.0_real64    !! extinction for diffuse light (-)
      real(real64) :: kdir           = 0.0_real64    !! extinction for direct light (-)
      real(real64) :: cofab          = 0.0_real64    !! interception coefficient (cm)
      real(real64) :: cfbs           = 1.0_real64    !! bare-soil ET factor (-)
      integer      :: swcf           = 0             !! crop factor switch
      integer      :: swcfbs         = 0             !! bare-soil factor switch
      real(real64) :: gird           = 0.0_real64    !! gross irrigation depth (cm)
      logical      :: flCropEmergence = .false.      !! crop emerged flag
      real(real64) :: et0            = 0.0_real64    !! potential ET (cm/d)
      real(real64) :: ew0            = 0.0_real64    !! potential evap wet crop (cm/d)
      real(real64) :: es0            = 0.0_real64    !! potential evap bare soil (cm/d)
   contains
      procedure :: init => crop_state_init
   end type crop_state_t

contains
   subroutine crop_state_init(self, crop_cfg)
      use crop_config_mod, only: crop_config_t
      class(crop_state_t),  intent(inout) :: self
      type(crop_config_t),  intent(in)    :: crop_cfg
      ! Defaults; runtime seeders update each timestep.
      ! crop_cfg arg reserved for future crop_init expansion.
   end subroutine
end module
```

Wired into `swap_state_t` as `type(crop_state_t) :: crop`.

### `atmosphere_state_t` extensions

Six blocks, ~50 new fields total.

**Block 1 — Daily meteo input arrays (366-sized fixed):**

```fortran
real(real64) :: arad(366) = 0.0_real64    !! daily radiation
real(real64) :: atmn(366) = 0.0_real64    !! daily min temperature
real(real64) :: atmx(366) = 0.0_real64    !! daily max temperature
real(real64) :: ahum(366) = 0.0_real64    !! daily humidity
real(real64) :: awin(366) = 0.0_real64    !! daily wind speed
real(real64) :: arai(366) = 0.0_real64    !! daily rainfall
real(real64) :: aetr(366) = 0.0_real64    !! daily ETref
real(real64) :: wet(366)  = 0.0_real64    !! daily wet fraction
```

**Block 2 — Sub-daily detailed arrays (96-sized):**

```fortran
real(real64) :: atav(96)  = 0.0_real64    !! sub-daily air temp
real(real64) :: epot(96)  = 0.0_real64    !! sub-daily potential evap
real(real64) :: tpot(96)  = 0.0_real64    !! sub-daily potential transp
real(real64) :: grain(96) = 0.0_real64    !! sub-daily gross rain
real(real64) :: nrain(96) = 0.0_real64    !! sub-daily net rain
```

**Block 3 — Derived meteo scalars:**

```fortran
real(real64) :: Tav        = 0.0_real64
real(real64) :: tavd       = 0.0_real64
real(real64) :: rh         = 0.0_real64
integer      :: daynrfirst = 0
integer      :: daynrlast  = 0
real(real64) :: atmin7(7)  = 0.0_real64
integer      :: nofd       = 0
real(real64) :: teprrain   = 0.0_real64
real(real64) :: teprsnow   = 0.0_real64
```

**Block 4 — Interception state/params:**

```fortran
real(real64) :: siccapact          = 0.0_real64
real(real64) :: fimin              = 0.0_real64
integer      :: isua               = 0
real(real64) :: avevaptb(2*magrs)  = 0.0_real64
real(real64) :: avprectb(2*magrs)  = 0.0_real64
real(real64) :: pfreetb(2*magrs)   = 0.0_real64
real(real64) :: pstemtb(2*magrs)   = 0.0_real64
real(real64) :: scanopytb(2*magrs) = 0.0_real64
```

**Block 5 — Runoff-CN method:**

```fortran
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

**Block 6 — Output flag toggles:**

```fortran
logical :: out_tmn = .false.
logical :: out_tmx = .false.
logical :: out_hum = .false.
logical :: out_win = .false.
logical :: out_etr = .false.
logical :: out_wet = .false.
logical :: out_rad = .false.
```

Plus type-bound `procedure :: init => atmosphere_state_init` taking `config%meteo` arg.

### `config%meteo` extension

```fortran
real(real64) :: cofred = 0.35_real64  !! reduction coefficient (et.f90 evap reduction)
```

Verify legacy default of 0.35 in `config_to_variables` during A9.

### `flrunon` / `runonarr` placement

These are runon-related fields read by boundtop (deferred in GR-BH Task 20). Inspect `src/config/general_config.f90` and `simulation_config.f90`; if absent, add to `simulation_config_t` (or whichever section makes sense given the legacy population path). Verify in A9.

### `et0/ew0/es0` cross-scope

These appear in both interception (atmosphere code) and crop computations. **Decision:** they live on `state%crop` (crop-driven potential ET). Atmosphere code reads via `state%crop%X`. Revisit if Phase B/C surfaces ambiguity.

## Phase A — Tasks (13)

| # | Task | Files |
|---|---|---|
| A1 | Pre-flight baseline marker | (verify only) check-full 5/5 |
| A2 | Create `crop_state_mod` | NEW `src/state/crop_state.f90`; wire into `swap_state_t`; meson + tests/unit/meson |
| A3 | Extend `atmosphere_state_t` — Block 1 (8 daily meteo arrays) | `src/state/atmosphere_state.f90` |
| A4 | Extend `atmosphere_state_t` — Block 2 (5 sub-daily arrays) | `src/state/atmosphere_state.f90` |
| A5 | Extend `atmosphere_state_t` — Block 3 (9 derived scalars) | `src/state/atmosphere_state.f90` |
| A6 | Extend `atmosphere_state_t` — Block 4 (8 interception fields) | `src/state/atmosphere_state.f90` |
| A7 | Extend `atmosphere_state_t` — Block 5 (10 CN fields) | `src/state/atmosphere_state.f90` |
| A8 | Extend `atmosphere_state_t` — Block 6 (7 output flags) | `src/state/atmosphere_state.f90` |
| A9 | Extend `config%meteo` (+`cofred`) and verify `flrunon`/`runonarr` placement | `src/config/meteorology_config.f90`, possibly `simulation_config.f90`; populator in `config_to_variables.f90` |
| A10 | Dual-write Block 1 + Block 2 in swap_mod | `src/core/swap_mod.f90` |
| A11 | Dual-write Blocks 3 + 4 + 5 in swap_mod | `src/core/swap_mod.f90` |
| A12 | Dual-write Block 6 + crop_state seeding | `src/core/swap_mod.f90` |
| A13 | Phase A close marker — check-full 5/5 | (verify only) |

**Phase A close criteria:**
- `grep "state%crop%lai\|state%atmosphere%arad\|state%atmosphere%Tav\|state%atmosphere%CNref\|state%atmosphere%out_tmn" src/core/swap_mod.f90` returns at least 5 hits.
- No reader changes — every legacy global still load-bearing.
- check-full 5/5 byte-for-byte.

## Phase B — Tasks (17)

| # | Task | Files | Notes |
|---|---|---|---|
| B1 | Introduce type-bound `state%atmosphere%init` | `src/state/atmosphere_state.f90`, `src/core/swap_mod.f90` | Signature `init(config%meteo)`; mirrors heat pilot |
| B2 | Migrate `et.f90` body | `src/atmosphere/et.f90` | `swredu` → `config%meteo%swredu`; `cofred` → `config%meteo%cofred` |
| B3 | Verify `et.f90` clean | grep | Empty (or only logf deferral) |
| B4 | Migrate `interception.f90` — InterceptionAgric (line 40 site) | `src/atmosphere/interception.f90` | `gird/kdif/kdir/cofab/lai` → `state%crop%X`; `isua` → `state%atmosphere%isua` |
| B5 | Migrate `interception.f90` — VonHoyningenHune (line 89 site) | Same file | `gird` → `state%crop%gird`; `avevaptb/avprectb/pfreetb/pstemtb/scanopytb` → `state%atmosphere%X`; `isua` |
| B6 | Migrate `interception.f90` — UpdateInterception (line 171 site) | Same file | `logf` retained narrow (Arc 9); `siccapact/fimin/ew0` → state%atmosphere or state%crop (verify per-symbol) |
| B7 | Migrate `interception.f90` — InterceptionTrees (line 409 site) | Same file | `isua/gird` → state; `nird` retained narrow (Arc 8) |
| B8 | Verify `interception.f90` clean | grep | Empty modulo `logf` + `nird` |
| B9 | Migrate `meteoday.f90` — CN runoff (line 102 site) | `src/atmosphere/meteoday.f90` | CN symbols → state%atmosphere%X |
| B10 | Migrate `meteoday.f90` — meteo daily (line 289 site) | Same file | Daily arrays + derived scalars + output flags → state%atmosphere; `swrain` → `config%meteo%swrain` |
| B11 | Migrate `meteoday.f90` — ET/interception driver (line 543 site) | Same file | Crop refs → state%crop%X; `swinter/swmetdetail/nmetdetail/swetr` → `config%meteo%X` |
| B12 | Verify `meteoday.f90` clean | grep | Empty modulo documented deferrals |
| B13 | Migrate `meteodt.f90` — bare line 58 | `src/atmosphere/meteodt.f90` | Convert bare to explicit imports |
| B14 | Migrate `meteodt.f90` — rain timing (line 139 site) | Same file | `swrain/raintab` → `config%meteo`; `wet/arai/nmrain/rainamount/rainfluxarray/raintimearray` → state%atmosphere (verify input vs runtime per-symbol) |
| B15 | Migrate `meteodt.f90` — bare lines 358 + 456 | Same file | Same pattern |
| B16 | Verify `meteodt.f90` clean | grep | Empty modulo documented deferrals |
| B17 | Phase B close marker — check-full 5/5 | (verify only) |

**Phase B close criteria:**
- `grep -n "^[[:space:]]*use variables" src/atmosphere/*.f90 | grep -v "logf\|nird"` empty
- `state%atmosphere%init` is type-bound; init refactor complete
- check-full 5/5 byte-for-byte

## Phase C — Tasks (16)

Codebase sweep + global retirement. Phase C does NOT touch any other `use variables` symbols — surgical replacement only.

### Symbol replacement reference

| Was | Now |
|---|---|
| `arad/atmn/atmx/ahum/awin/arai/aetr/wet` | `state%atmosphere%X` |
| `atav/epot/tpot/grain/nrain` | `state%atmosphere%X` |
| `Tav/tavd/rh/daynrfirst/daynrlast/atmin7/nofd/teprrain/teprsnow` | `state%atmosphere%X` |
| `siccapact/fimin/isua/avevaptb/avprectb/pfreetb/pstemtb/scanopytb` | `state%atmosphere%X` |
| `CNref/CNdry/CNwet/ThetaRef/Runoff_CN/wc_cor/wc10/iCNtab/CNtimTAB/CNrefTAB` | `state%atmosphere%X` |
| `out_tmn/out_tmx/out_hum/out_win/out_etr/out_wet/out_rad` | `state%atmosphere%X` |
| `lai/kdif/kdir/cofab/swcf/swcfbs/cfbs/gird/flCropEmergence/et0/ew0/es0` | `state%crop%X` |

| # | Task | Files |
|---|---|---|
| C1 | Cluster `src/core/` | `swap_mod.f90`, `initialize.f90`, `timecontrol_mod.f90`, `swap_bmi_mod.f90`, `swap_capi_mod.f90` |
| C2 | Cluster `src/io/` | `swap_csv_output.f90`, `swapoutput.f90`, `readmeteo.f90` |
| C3 | Cluster `src/boundary/` (close prior deferrals) | `boundbottom.f90`, `boundtop.f90` |
| C4 | Cluster `src/heat/` (close Tav/atav deferral) | `temperature.f90` |
| C5 | Cluster `src/soil/` | `soilgrid.f90`, `soilhydraulics.f90`, `waterbalance.f90` |
| C6 | Cluster `src/drainage/` | `drainage.f90`, `surfacewater.f90`, `divdra.f90` |
| C7 | Cluster `src/crop/` | All 9 crop files |
| C8 | Cluster `src/solute/` | `solute.f90`, `agetracer.f90` |
| C9 | Self-audit src/atmosphere/ | The 4 atmosphere files — ensure no residual after schema-level dual-write cleanup |
| C10 | Audit pass — find stale reads (4 grep audits) | Mandatory gate before C11 |
| C11 | Drop legacy meteo writes from adapter | `src/io/toml/config_to_variables.f90`, `src/io/readmeteo.f90` |
| C12 | Drop crop-state legacy writes | `src/io/toml/config_to_variables.f90`, possibly readmeteo or cropwofost_init |
| C13 | Drop output-flag + CN method legacy writes | `src/io/toml/config_to_variables.f90` |
| C14 | Delete meteo + interception + CN globals from `variables.f90` | `src/core/variables.f90`, `src/core/initialize.f90` |
| C15 | Delete crop globals from `variables.f90` | `src/core/variables.f90`, `src/core/initialize.f90` |
| C16 | Final verification + arc-complete marker | (verify only) — check-full 5/5; BMI + cffi-demo; sanity greps; arc-complete commit |

**Phase C close criteria:**
- `grep -nE "^[[:space:]]*real\(8\).*\b(arad|atmn|atmx|ahum|awin|arai|aetr|wet|atav|Tav|tavd|rh|CNref|out_tmn|lai|gird)\b" src/core/variables.f90` empty
- `grep -nE "^[[:space:]]*use variables" src/atmosphere/*.f90 | grep -v "logf\|nird"` empty
- check-full 5/5 byte-for-byte
- BMI + cffi-demo passing

## Risk & Mitigation

| Risk | Mitigation |
|---|---|
| Per-symbol scope ambiguity in interception.f90 | Per-symbol verification by implementer in B6; spec table is starting point not gospel |
| Sub-daily array semantics (swmetdetail==1 only) | Fields exist regardless of mode; only populated in detail mode; legacy semantics preserved |
| Output flags read-by-atmosphere vs set-by-adapter | A9 verifies the config_to_variables → state%atmosphere%out_X seed path |
| 50+ field schema additions break incremental builds | `rm -rf builddir` mandatory per `feedback_state_schema_clean_rebuild.md` |
| Crop overlap symbols may be set by atmosphere code, not crop code | If A12 / B-tasks find writes that should originate from crop, dual-write transitionally — Arc 8 cleanup follows |
| ~50 schema additions explode atmosphere_state.f90 size | Acceptable — single subrecord, single responsibility (atmosphere meteo data). Alternative is splitting into atmosphere_meteo / atmosphere_interception / atmosphere_runoff sub-records, which would dilute the established pattern. |

## What this Arc DOES NOT Change

- BMI / cffi surfaces — continue working throughout
- Other `use variables` consumers in boundary/heat/soil/drainage/io/crop — keep their non-atmosphere globals
- Physics — every value computed by every formula stays bit-identical

## Effort

- **Total tasks:** 46 (A1–A13 + B1–B17 + C1–C16)
- **Estimated effort:** 8–11 days, subagent-driven
- **Phases gated by check-full 5/5 byte-for-byte**

## Memory & ADR Consequences

- ADR candidate: ADR 0043 — "crop_state_t introduced" — separate ADR documenting the new subrecord per ADR pattern.
- After arc closes, update memory `project_state_rescue_complete_2026-05-12.md` listing GR-ATM complete.
- Two strangler-fig deferrals close in this arc: Tav/atav (Arc 3 territory) and swkmean/swredu/flrunon/runonarr (Arc 3 territory).

## Files & Artifacts

- **This spec:** `docs/superpowers/specs/2026-05-13-globals-atmosphere-design.md`
- **Plan (next):** `docs/superpowers/plans/2026-05-13-globals-atmosphere.md`
- **Roadmap context:** `docs/superpowers/plans/2026-05-13-globals-retirement-roadmap.md` (Arc 4 implementation)
- **GR-BH precedent:** `docs/superpowers/specs/2026-05-13-globals-boundary-heat-design.md` (Arc 2+3, completed 2026-05-13)
