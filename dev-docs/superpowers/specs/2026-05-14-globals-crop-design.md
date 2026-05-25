---
title: "GR-CROP — Globals Retirement Arc 8 (Crop + Nutrients + Rain timing + Cross-arc closures)"
date: 2026-05-14
status: approved
context: globals-retirement Arc 8 — full scope, the one-large-arc finale before Arc 9 final retirement
---

# GR-CROP Design

## Goal

Migrate the entire crop subsystem off `use variables`, plus the nutrients subsystem (management_soil + wofost_soil_*), plus rain-timing arrays in meteodt — and close all inherited cross-arc deferrals from GR-BH and GR-ATM. This is the largest arc by far; after it, only Arc 9 (final retirement) remains.

Sub-objectives:

1. Expand `state%crop` with sub-records per crop type: `%common`, `%fixed`, `%wofost`, `%grass` plus the existing 12 atmospheric-exchange fields at top level (from GR-ATM Phase A).
2. NEW `state%nutrients` subrecord covering management_soil + wofost_soil_* nutrient state.
3. Rain timing: add `nmrain/rainamount/rainfluxarray/raintimearray` to `state%atmosphere`; add `raintab(60)` to `config%meteo`. Migrate meteodt readers.
4. Migrate 12 crop files off `use variables`:
   - 4 init files: `cropfixed_init`, `cropgrass_init`, `cropwofost_init`, `wofost_soil_parameters`
   - 5 crop runtime: `cropgrowth.f90` (4937L — the giant), `oxygenstress.f90` (1840L incl. GET_MAX_RESP_FACTOR), `rootextraction.f90` (966L), `tillage.f90`, `irrigation.f90`
   - `management_soil.f90` (632L) — nutrients subsystem
   - `wofost_soil_*.f90` cluster (8 files, mostly small)
   - `meteodt.f90` rain timing site (atmosphere file, addressed for cross-arc closure)
5. Close inherited deferrals: GET_MAX_RESP_FACTOR state threading (GR-ATM); irrigation `gird/isua` write-side cleanup (GR-ATM); rain timing (GR-ATM); `boundtop` config threading through `soilhydraulics→headcalc` (GR-BH) if it falls out naturally.
6. Retire crop + atmosphere multi-consumer globals from `variables.f90` once readers migrated.

**Closes inherited deferrals from prior arcs:**

| Deferral | Origin | Closure path |
|---|---|---|
| GET_MAX_RESP_FACTOR `tav` reads | GR-ATM | State arg threading in oxygenstress |
| `gird/isua` write-side cleanup | GR-ATM | Source state writes at write sites; drop legacy mirror |
| Rain timing (nmrain/rainamount/rainfluxarray/raintimearray) | GR-ATM/Arc 7 | New state%atmosphere fields + reader migration |
| `raintab` | GR-ATM/Arc 7 | New config%meteo field |
| swap_mod A12 init-time seeding cleanup | GR-ATM | After all consumers migrate, drop the dual-write |
| Atmosphere multi-consumer globals (arad/atmn/.../tav) | GR-ATM Phase C | Retire once cropgrowth/nutrients consumers migrate |
| Crop atmospheric-exchange field globals (lai/gird/kdif/.../es0) | GR-ATM Phase C | Retire once internal crop consumers migrate to state%crop |

## Architecture

Four sub-phases:

- **Phase A** — Schema additions: 4 crop sub-records, nutrients subrecord, rain timing fields, raintab config. Dual-write seeding at swap_mod init.
- **Phase A.5** — Runtime dual-write coverage. Every write site for the ~70 migrated symbols mirrors to state. Critical lesson from GR-ATM: zero-resets count as writes (see commit `8fb76cf` — nocrop() lai-reset).
- **Phase B** — 12-file reader cutover, ordered by file. cropgrowth.f90 split across B5/B6/B7 (one task per use-variables site) due to its 4937-line size.
- **Phase C** — Codebase sweep, audit pass, legacy retirement.

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

- **Controller quality bar (per user directive):** subagents do NOT defer due to missing infrastructure. If a deferral surfaces (state arg not threaded, sibling reader still on legacy, schema gap, etc.), the controller intervenes with an inline fix — thread state through the chain, extend Phase A.5, expand schema. No cascading deferrals.

## Schema Additions

### NEW: `crop_state_t` restructured (existing top-level fields preserved)

```fortran
module crop_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use crop_common_state_mod,  only: crop_common_state_t
   use crop_fixed_state_mod,   only: crop_fixed_state_t
   use crop_wofost_state_mod,  only: crop_wofost_state_t
   use crop_grass_state_mod,   only: crop_grass_state_t
   implicit none
   private
   public :: crop_state_t

   type :: crop_state_t
      ! Existing GR-ATM A12 atmospheric-exchange fields (preserved)
      real(real64) :: lai = 0, kdif = 0, kdir = 0, cofab = 0, cfbs = 1
      integer      :: swcf = 0, swcfbs = 0
      real(real64) :: gird = 0
      logical      :: flCropEmergence = .false.
      real(real64) :: et0 = 0, ew0 = 0, es0 = 0

      ! [SS-GR-CROP] runtime sub-records per crop type + shared common
      type(crop_common_state_t) :: common
      type(crop_fixed_state_t)  :: fixed
      type(crop_wofost_state_t) :: wofost
      type(crop_grass_state_t)  :: grass
   contains
      procedure :: init => crop_state_init
   end type
end module
```

### NEW: `crop_common_state_t`

Fields shared across all crop types. Representative list — implementer audits cropgrowth/init files to confirm full set:

```fortran
type :: crop_common_state_t
   integer      :: daycrop        = 0    !! days since crop start
   real(real64) :: dvs            = 0    !! development stage
   real(real64) :: tsum           = 0    !! temperature sum
   real(real64) :: rd             = 0    !! actual rooting depth (cm)
   real(real64) :: rdpot          = 0    !! potential rooting depth (cm)
   real(real64) :: rdm            = 0    !! max rooting depth from soil profile
   real(real64) :: rri            = 0    !! max daily root depth increase
   real(real64) :: rdi            = 0    !! initial rooting depth
   real(real64) :: rdc            = 0    !! max crop rooting depth
   real(real64) :: ch             = 0    !! crop height
   real(real64) :: cf             = 0    !! crop factor
   real(real64) :: laipot         = 0    !! potential LAI
   real(real64) :: cuptgraz       = 0
   real(real64) :: cuptgrazpot    = 0
   real(real64) :: HarLosOrm_tot  = 0    !! cumulative harvest loss organic matter
   integer      :: daynrsta       = 0
   integer      :: daynrend       = 0
   integer      :: swcrp          = 0    !! 1=fixed, 2=wofost, 3=grass
   integer      :: icrop          = 0
   logical      :: flcropcalendar = .false.
   logical      :: flcropoutput   = .false.
   logical      :: flcropnut      = .false.
   logical      :: flharvestday   = .false.
   logical      :: swend          = .false.
contains
   procedure :: init => crop_common_state_init
end type
```

### NEW: `crop_fixed_state_t`

```fortran
type :: crop_fixed_state_t
   real(real64) :: cftb(2*magrs) = 0    !! crop factor or height vs DVS
   ! Other fixed-crop-only runtime — audit during implementation
contains
   procedure :: init => crop_fixed_state_init
end type
```

### NEW: `crop_wofost_state_t`

WOFOST biomass pools, gross/death/decay flows:

```fortran
type :: crop_wofost_state_t
   ! Biomass pools (actual + potential)
   real(real64) :: wlv  = 0, wlvpot  = 0    !! leaf weight
   real(real64) :: wst  = 0, wstpot  = 0    !! stem weight
   real(real64) :: wrt  = 0, wrtpot  = 0    !! root weight
   real(real64) :: wso  = 0, wsopot  = 0    !! storage organ weight
   real(real64) :: twlv = 0, twst   = 0     !! cumulative
   real(real64) :: tagp = 0, tagppot = 0    !! total above-ground
   real(real64) :: tagpt = 0, tagptpot = 0  !! total above-ground transpiration
   real(real64) :: cwdm = 0, cwdmpot = 0    !! cumulative weight dry matter
   real(real64) :: pgass = 0, pgasspot = 0  !! gross assimilation

   ! Death/decay flows
   real(real64) :: dwlv = 0, dwlvpot = 0
   real(real64) :: dwst = 0, dwstpot = 0
   real(real64) :: dwrt = 0, dwrtpot = 0
   real(real64) :: dwso = 0, dwsopot = 0
   real(real64) :: dwlvCrop = 0, dwlvSoil = 0

   ! Losses
   real(real64) :: plossdm = 0, lossdm = 0

   ! Bulb-crop fields
   logical      :: swbulb  = .false.
   real(real64) :: wbl = 0, wblpot = 0, dwbl = 0, dwblpot = 0
contains
   procedure :: init => crop_wofost_state_init
end type
```

### NEW: `crop_grass_state_t`

```fortran
type :: crop_grass_state_t
   real(real64) :: plwt = 0     !! plant weight (live)
   ! mowing/grazing state + swardstate — audit during implementation
contains
   procedure :: init => crop_grass_state_init
end type
```

### NEW: `nutrients_state_t`

Hosts management_soil + wofost_soil_* nutrient state. Exact field list determined by Phase A11 audit. Representative skeleton:

```fortran
module nutrients_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: nutrients_state_t

   type :: nutrients_state_t
      ! Soil organic matter / nitrogen pools (per layer; allocatable to nlay or numnod)
      real(real64), allocatable :: fom(:)        !! fresh organic matter
      real(real64), allocatable :: bio(:)        !! microbial biomass
      real(real64), allocatable :: hum(:)        !! humus
      real(real64), allocatable :: nh4(:)        !! ammonium
      real(real64), allocatable :: no3(:)        !! nitrate
      ! Crop N demand/uptake (WOFOST)
      real(real64) :: cnflv  = 0    !! N concentration leaves
      real(real64) :: cnfst  = 0    !! N concentration stems
      real(real64) :: nuptot = 0    !! cumulative N uptake
      ! Process flux scratchpads (per-timestep)
      real(real64) :: nflux  = 0
      ! Additional fields surfaced during Phase A11 audit
   contains
      procedure :: init => nutrients_state_init
   end type
end module
```

### `state%atmosphere` rain timing additions

```fortran
integer      :: nmrain                       = 0
real(real64) :: rainamount(mrain)            = 0.0_real64
real(real64) :: rainfluxarray(mrain)         = 0.0_real64
real(real64) :: raintimearray(mrain)         = 0.0_real64
```

### `config%meteo` raintab addition

```fortran
real(real64) :: raintab(60) = 0.0_real64    !! intensity vs time (swrain==1 input)
```

### `swap_state_t` extensions

Add `type(nutrients_state_t) :: nutrients` field. `state%crop` is already wired; the 4 sub-records are inside it.

### Build wiring

- `meson.build`: add new sub-record sources **before** `src/state/swap_state.f90`. Dependency order:
  1. `crop_common_state.f90`, `crop_fixed_state.f90`, `crop_wofost_state.f90`, `crop_grass_state.f90`
  2. `crop_state.f90` (modify existing — add `use` lines for the 4 sub-records)
  3. `nutrients_state.f90`
  4. `swap_state.f90` (modify existing — add `nutrients` field)

- `tests/unit/meson.build`: mirror the same ordering in `pfunit_extra_sources`.

## Phase A — Tasks (18)

| # | Task | Files |
|---|---|---|
| A1 | Pre-flight baseline marker | (verify only) |
| A2 | Create `crop_common_state_mod` (empty type body) | NEW + meson + wire into crop_state |
| A3 | Create `crop_fixed_state_mod` | NEW + meson + wire |
| A4 | Create `crop_wofost_state_mod` | NEW + meson + wire |
| A5 | Create `crop_grass_state_mod` | NEW + meson + wire |
| A6 | Populate `crop_common_state_t` fields | Audit cropgrowth + init files for shared symbols |
| A7 | Populate `crop_wofost_state_t` fields | Audit cropwofost_init + cropgrowth wofost paths |
| A8 | Populate `crop_grass_state_t` fields | Audit cropgrass_init + cropgrowth grass paths |
| A9 | Populate `crop_fixed_state_t` fields | Audit cropfixed_init + cropgrowth fixed paths |
| A10 | Create `nutrients_state_mod` | NEW + meson + wire into swap_state |
| A11 | Populate `nutrients_state_t` from management_soil + wofost_soil_* audit | Inventory blanket-use Variables symbols, classify per-symbol |
| A12 | Extend `atmosphere_state_t` with rain timing (4 fields) | nmrain + 3 arrays |
| A13 | Extend `config%meteo` with raintab + adapter populator | meteorology_config.f90 + config_to_variables.f90 |
| A14 | Dual-write crop_common in swap_mod | After atmosphere init |
| A15 | Dual-write crop_wofost/fixed/grass in swap_mod | Branch by swcrp OR seed all defensively |
| A16 | Dual-write nutrients in swap_mod | After soilwater init |
| A17 | Dual-write rain timing in swap_mod | After atmosphere init |
| A18 | Phase A close marker — check-full 5/5 | (verify only) |

**Phase A close criteria:** all new state types declared and seeded; check-full 5/5; no reader changes.

## Phase A.5 — Tasks (6)

Runtime dual-write coverage. Critical lesson: GR-ATM commit `8fb76cf` — zero-reset writes count as writes. Audit ALL writes, including conditional zero-resets in state-machine code paths.

| # | Task | Files |
|---|---|---|
| A5.1 | Add dual-writes in `cropgrowth.f90` | Largest target — every `dvs/tsum/rd/...wlv/wst/wrt/wso/...` write |
| A5.2 | Add dual-writes in `cropfixed_init/cropgrass_init/cropwofost_init` | Per-file batch |
| A5.3 | Add dual-writes in `irrigation/oxygenstress/rootextraction/tillage` | Per-file batch |
| A5.4 | Add dual-writes in `management_soil` + `wofost_soil_*` for nutrients | Nutrients write sites |
| A5.5 | Add dual-writes in `meteodt.f90` for rain timing | nmrain/rainamount/... writes |
| A5.6 | Phase A.5 close marker — check-full 5/5 | (verify only) |

**Controller intervention triggered if:**
- Subagent reports "subroutine X doesn't take state" → controller threads state through the caller chain in a pre-task.
- Subagent reports "no obvious mirror target for symbol Y" → controller extends schema (add field to appropriate state sub-record) and re-dispatches.

## Phase B — Tasks (17)

| # | Task | Files | Notes |
|---|---|---|---|
| B1 | Migrate `cropfixed_init.f90` (1 site) | line 28 narrow | Reads `idev/tsumea/tsumam/tbase/...` — most are config-loaded; route to `config%crop` or migrate to state%crop%common |
| B2 | Migrate `cropgrass_init.f90` (2 sites) | line 34, 417 | |
| B3 | Migrate `cropwofost_init.f90` (2 sites) | line 42, 490 | Heavy nutrient interaction; nutrient assignments may need state%nutrients refs |
| B4 | Migrate `wofost_soil_parameters.f90` (1 narrow site) | line 14: `BDENS` | `BDENS` → `state%soilwater%X` (likely already there) or to nutrients |
| B5 | Migrate `cropgrowth.f90` site 1 (line 24 bare) | Carve subroutines reading via this site | The `dummy_tsoil_cg_ => tsoil` rename is interesting — preserve intent |
| B6 | Migrate `cropgrowth.f90` site 2 (line 549 bare) | |
| B7 | Migrate `cropgrowth.f90` site 3 (line 787 bare) | |
| B8 | Verify cropgrowth.f90 clean | grep |
| B9 | Migrate `oxygenstress.f90` (2 sites) | Close GR-ATM GET_MAX_RESP_FACTOR deferral via state arg threading |
| B10 | Migrate `rootextraction.f90` (3 sites) | |
| B11 | Migrate `irrigation.f90` (2 sites) | Close GR-ATM gird/isua write-side cleanup |
| B12 | Migrate `tillage.f90` (1 narrow site) | line 10: swhyst/swsolu/swoxygen/Bdens/ParamVG — config |
| B13 | Migrate `management_soil.f90` (blanket) | Largest non-cropgrowth file; nutrients heavy |
| B14 | Migrate `wofost_soil_*.f90` cluster | Small files, narrow imports |
| B15 | Migrate `meteodt.f90` rain timing site (line 143) | Close GR-ATM rain-timing deferral |
| B16 | Verify all crop files clean | grep |
| B17 | Phase B close marker — check-full 5/5 | (verify only) |

**Phase B close criteria:** all 12 target files migrated off bare `use variables`; only documented narrow deferrals remain (e.g., `logf`); check-full 5/5.

## Phase C — Tasks (16)

| # | Task | Files / Scope |
|---|---|---|
| C1 | Cluster sweep — `src/core/` | mesh/initialize/swap_mod for crop/nutrient reads (state references) |
| C2 | Cluster sweep — `src/io/` | swap_csv_output, swapoutput, readmeteo crop/nutrient/rain refs |
| C3 | Cluster sweep — `src/atmosphere/` | meteoday rain-timing reads; close meteodt-side leftovers |
| C4 | Cluster sweep — `src/boundary/` | Close any crop/nutrient bleed |
| C5 | Cluster sweep — `src/heat/` | minimal exposure expected |
| C6 | Cluster sweep — `src/soil/` | soilhydraulics root-depth feedback (reads `rd` from crop) |
| C7 | Cluster sweep — `src/drainage/` | minimal |
| C8 | Cluster sweep — `src/solute/` | solute+age tracer crop-runtime refs |
| C9 | Audit pass — find stale reads | MANDATORY gate; 5+ audit greps (crop fields, nutrient fields, atmosphere multi-consumer, rain timing, ETSine residuals) |
| C10 | Drop legacy writes from `config_to_variables.f90` for crop/nutrients/rain | Where state mirror sources directly from config |
| C11 | Drop legacy writes from `readmeteo.f90` for rain timing | State writes become sole |
| C12 | Retire write-only globals from `variables.f90` | Symbols with no remaining readers (most of the ~70 migrated set) |
| C13 | Retire atmosphere multi-consumer globals | arad/atmn/atmx/ahum/awin/arai/aetr/wet/atav/tav once all consumers migrated. Crop atmospheric-exchange globals (lai/gird/kdif/...) once internal consumers migrate |
| C14 | Drop swap_mod A12 init-time seeding for retired symbols | Per GR-ATM inheritance |
| C15 | Update Arc 9 (final retirement) plan with residual deferrals | Document boundtop config-threading, ETSine astronomical, swinco=3 inline block, transient buffers |
| C16 | Final verification + arc-complete marker | check-full 5/5; BMI + cffi-demo; sanity greps |

**Phase C close criteria:**
- `grep -nE "^[[:space:]]*real\(8\).*\b(lai|gird|kdif|kdir|cofab|swcf|swcfbs|cfbs|gird|flCropEmergence|et0|ew0|es0|wlv|wst|wrt|wso|tagp|dvs|tsum|rd|ch|cf|...)\b" src/core/variables.f90` empty
- `grep -nE "^[[:space:]]*use variables" src/crop/*.f90` empty modulo documented narrow deferrals
- check-full 5/5 byte-for-byte; BMI + cffi-demo passing

## Risk & Mitigation

| Risk | Mitigation |
|---|---|
| Schema explosion makes crop_state files large | Sub-records per crop type contain growth; each ~30-50 fields; manageable |
| cropgrowth.f90 (4937L) migration breaks regression mid-flight | Split into 3 tasks (B5/B6/B7), one per use-variables site; VG between each |
| Phase A.5 misses a write site (lesson from GR-ATM `8fb76cf`) | Audit greps in Phase A.5 explicitly include conditional zero-resets; controller cross-checks |
| Sub-records-per-crop-type contradicts ADR 0033 | Documented architectural deviation — crop types are disjoint algorithms, not redundant cohorts |
| Nutrients field inventory is unclear before A11 audit | A11 explicitly is the inventory step; A.5 dual-writes follow |
| GET_MAX_RESP_FACTOR state threading requires touching multiple callers | Controller does the threading inline if subagent defers; per user directive |
| Rain timing writes are tightly coupled to meteodt logic | A.5.5 audits all writes; B15 migrates reads |
| Phase C global deletion exposes missed readers | C9 mandatory audit before C10+; compile-error iterative patches |

## What this Arc DOES NOT Change

- BMI / cffi surfaces — continue working throughout
- Other arcs' deferred items not listed in cross-arc closures section (e.g., ETSine astronomical scratchpad)
- Physics — every value computed by every formula stays bit-identical

## Effort

- **Total tasks:** 57 (A1-A18 + A5.1-A5.6 + B1-B17 + C1-C16)
- **Estimated effort:** 12-18 days, subagent-driven
- **Phases gated by check-full 5/5 byte-for-byte**

## Memory & ADR Consequences

- ADR candidate: ADR 0044 — crop_state_t sub-records-per-type (deviation from ADR 0033 flatten principle, documented rationale).
- ADR candidate: ADR 0045 — nutrients_state_t introduction.
- After arc closes, update memory `project_state_rescue_complete_2026-05-12.md` listing GR-CROP complete.
- Phase A.5 lesson (zero-resets count as writes) — codify in `feedback_runtime_dualwrite_coverage.md` if recurs.

## Files & Artifacts

- **This spec:** `docs/superpowers/specs/2026-05-14-globals-crop-design.md`
- **Plan (next):** `docs/superpowers/plans/2026-05-14-globals-crop.md`
- **Roadmap context:** `docs/superpowers/plans/2026-05-13-globals-retirement-roadmap.md` (Arc 8 implementation)
- **GR-ATM precedent:** `docs/superpowers/specs/2026-05-13-globals-atmosphere-design.md` (Arc 4, completed 2026-05-14 — commit `57b00b3`)
