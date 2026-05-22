# `src/atmosphere` — Atmospheric boundary subsystem

## Overview

The `atmosphere` subsystem is the upper boundary of SWAP. It converts the raw
meteorological forcing — radiation, temperature, humidity, wind, precipitation
and (optionally) reference ET — into the per-timestep fluxes that the soil and
crop subsystems consume:

- daily and sub-daily meteo ingestion and bookkeeping,
- rain / snow partitioning,
- snow accumulation, sublimation and melt,
- canopy interception via four interchangeable schemes
  (Von Hoyningen-Hune & Braden, Gash 1995, adapted Rutter / Sparse Gash,
  MSW1 sparse canopy),
- reference and actual evapotranspiration with Penman–Monteith and two
  soil-evaporation reduction models (Black, Boesten–Stroosnijder),
- SCS Curve-Number surface runoff.

Outputs are written into `state%atmosphere` (aggregated `atmosphere_state_t`
defined in `src/state/atmosphere_state.f90`), and read by the rest of the
model through `state` aliases — no module-level globals are exported.

### Recent arcs

- **2026-05-22 — GR-ATM-CLEAN** (this README's current layout): split the
  bundled `meteoday.f90` into four single-module files, populated
  `atmosphere_constants`, replaced the three remaining task-selector
  subroutines with named init/step pairs, retired the `MeteoVars`
  shared-state bag, split `ProcessMeteoDay` into daily + sub-daily
  orchestrators with four private helpers, and swept the F77 intrinsics
  (`dexp`/`dlog`/`dmin1`/`dble`) to F90 generics.
  See `docs/superpowers/specs/2026-05-22-atmosphere-cleanup-arc-design.md`
  and `docs/superpowers/plans/2026-05-22-atmosphere-cleanup-arc.md`.
- **2026-05-22 — PenMon derived-type refactor** (commit `c29408b`): wrapped
  `PenMon` arguments in `pm_inputs_t` / `pm_outputs_t`.

### Module layout

| File | Module(s) | Role |
|------|-----------|------|
| `atmosphere_constants.f90` | `atmosphere_constants_mod` | Named physical constants (water/snow latent heats, pond threshold, CN ratios, FC/WP heads, …) |
| `meteo_io.f90` | `meteo_process_mod` | Daily meteo ingestion (`ReadMeteoDay`, `ResetMetFlx`); calls `PartitionPrecipitation` |
| `meteo_orchestrator.f90` | `meteo_mod` | Daily ET / interception coordinator; `ProcessMeteoDay` dispatches to `process_meteo_day_daily` and `process_meteo_day_subdaily`, with four private helpers (`apply_interception_step`, `compute_reference_et`, `compute_wet_fraction`, `partition_peva_ptra`) |
| `meteodt.f90` | `meteodt_mod` | Sub-daily meteo orchestration (rain events, time-step fluxes, sine-wave ET) |
| `et.f90` | `et_mod` | Penman–Monteith ET (`PenMon`, `PenMon_calc`, `pm_inputs_t`, `pm_outputs_t`) and soil-evaporation reduction (`reduceva_daily`, `reduceva_dt`) |
| `runoff.f90` | `runoff_mod` | SCS Curve-Number runoff (`cn_init`, `cn_step`) |
| `snow.f90` | `snow_mod` | Snow accumulation, sublimation, melt (`snow_init`, `snow_step`) |
| `precipitation.f90` | `precipitation_mod` | Rain / snow partitioning (`PartitionPrecipitation`) |
| `interception.f90` | `interception_mod` | Four interception schemes + interception splitter |

### Dependencies

Imports from
- `src/core/` — `swap_state_mod`, `swap_config_mod`, `swap_constants`, `swap_array_dimensions`, `variables` (legacy bridge — being retired),
- `src/utils/` — `array_utils` (`afgen`), `soilhydraulics_utils` (`watcon`), `error_mod`.

Has no downstream physics dependencies — every consumer reads atmosphere
outputs through `state%atmosphere`.

---

## `atmosphere_constants.f90`

Single-module file (`atmosphere_constants_mod`) hosting the named physical
and numerical constants used across the subsystem. Consumers: `snow.f90`,
`runoff.f90`, `et.f90` (Phase B of GR-ATM-CLEAN, commit `c7ee97a`).

### Public parameters (excerpt)
- `SPECIFIC_HEAT_WATER` (4180 J kg⁻¹ K⁻¹) — was `cwat` literal in `snow.f90`.
- `LATENT_HEAT_MELTING` (333580 J kg⁻¹) — was `lm`.
- `SNOW_TEMPERATURE_C` (0.0) — was `ts`.
- `SNOW_LIQUID_WATER_FRACTION` (0.07) — liquid-water cap in snowpack.
- `SOIL_SURFACE_FREEZE_THRESHOLD_C` (0.5).
- `POND_THRESHOLD_CM` (1.0e-10) — was duplicated in `reduceva` and
  `ProcessMeteoDay`.
- `DEPTH_10CM_CM` (10.0) — CN reference layer depth.
- `H_FIELD_CAPACITY_CM` (-100), `H_WILTING_POINT_CM` (-16000) — CN moisture
  correction reference heads.
- `INITIAL_ABSTRACTION_RATIO` (0.2) — `Iₐ / S` for SCS-CN.

### Modernisation notes
- Use `real64` from `iso_fortran_env` — already done.
- Additional constants (Karman, AFGEN table sizes, the 1e-6 snowpack-presence
  threshold in `precipitation.f90`) can join this module as their call sites
  migrate.

---

## `meteo_io.f90` — daily meteo I/O (`meteo_process_mod`)

### Public signatures
- **`ReadMeteoDay(state, config)`** — reads today's meteo record from the
  pre-loaded `state%atmosphere%a*` arrays. Two branches:
  - `swmetdetail = 0` (daily): copies `arad/atmn/atmx/ahum/awin/arai/aetr`
    into transient locals + state, computes 24 h and daytime mean
    temperature, saturation vapour pressure (Tetens), relative humidity,
    and mirrors the values into `state%atmosphere%out_*` for PEARL CFO
    output.
  - `swmetdetail = 1` (sub-daily): validates per-record indices and
    timestamps against the meteo file, fills `state%atmosphere%arad/ahum/atav`
    and `state%atmosphere%awind_subdaily`, `state%atmosphere%arain_subdaily`.
  - Always calls `ResetMetFlx` first, then `PartitionPrecipitation` to
    split today's `grai` into rain / snow / `fprecnosnow`.

- **`ResetMetFlx(state)`** — zeros the intermediate and cumulative
  atmosphere flux cohorts when the controller raises `flZeroIntr` /
  `flZeroCumu`. Delegates to `state%atmosphere%intr%reset()` and
  `state%atmosphere%cumu%reset()` (post ADR-0033 cohort-reset pattern).

### Modernisation notes
- F77 intrinsics swept to generics (`dexp`→`exp`, `dlog`→`log`, etc.).
- The sub-daily `arain`/`awind` arrays now live on the state object
  (`arain_subdaily`, `awind_subdaily`) instead of the retired `MeteoVars`
  bag — added in commit `5fc73ac`.
- Still imports `swmetdetail, swetr, swhydrlift, …` from `variables`
  (DEFERRED globals — config-side migration pending).

---

## `meteo_orchestrator.f90` — daily ET / interception coordinator (`meteo_mod`)

### Public signatures
- **`ProcessMeteoDay(state, config)`** — top-level dispatcher; branches on
  `swmetdetail` to either `process_meteo_day_daily` or
  `process_meteo_day_subdaily`.

### Private orchestrators
- **`process_meteo_day_daily(state, config)`** — daily-resolution path:
  interception → reference ET → wet-canopy fraction → peva/ptra
  partition → `cn_step` runoff (if `swuseCN = 1`) → `reduceva_daily`.
- **`process_meteo_day_subdaily(state, config)`** — sub-daily path:
  per-record interception → per-record ET → aggregated daily
  `tpot/epot/grain/nrain(i)` + recomputed daily means.

### Private helpers (extracted in Phase E)
- **`apply_interception_step(state, config, aintc)`** — selects between
  `VonHHBraden` / `Gash` / `ruttervw` per `swinter` and calls `DivIntercep`.
- **`compute_reference_et(state, config, irecord, etr, hum_in, win_in, rcs, pmo)`**
  — wraps `PenMon` / supplied reference ET into the three potentials
  `es0 / et0 / ew0`, with `swcf` / `swcfbs` crop-factor switches.
- **`compute_wet_fraction(state, config, aintc, eintc, Tdirectwet, interc, irecord, wfrac)`**
  — derives the wet-canopy fraction.
- **`partition_peva_ptra(state, config, wfrac, Edirect, Tdirect, Edirectpond)`**
  — splits potentials into `peva` (cover-corrected, ponding-corrected,
  PMdirect when `swdivide = 1`) and `ptra` (with optional CO₂
  correction `fco2tra`).

### Modernisation notes
- This file is the successor to the ~430-line bundled
  `ProcessMeteoDay` body. Each helper is a candidate for `pure` once the
  remaining `variables` reads are migrated.
- The intrinsics sweep (Phase F) cleaned up `dexp`/`dlog`/`dmin1`/`dble`
  call chains.

---

## `meteodt.f90`

Sub-daily meteorology coordinator. Activated when at least one of
`flrainintens`, `flmeteodt`, `fletsine` is set in `state%timecontrol`.

### Public signatures
- **`MeteoDT(state)`** — top-level dispatcher. On the first call of a
  calendar year (`flYearStart`) runs `ProcessRainEvents`. Each timestep,
  optionally runs `ProcessMeteoTsteps` and `ETSine`.

### Private signatures
- **`ProcessRainEvents(state)`** — builds the per-year
  `(raintimearray, rainfluxarray)` series for `swrain ∈ {1,2,3}`:
  - `swrain = 1`: daily depth × tabulated intensity (`raintab` via
    `afgen`) → duration `t = P / I` → two-point flux record per rainy
    day.
  - `swrain = 2`: daily depth × user-supplied wet duration (`wet(i)`).
  - `swrain = 3`: event-based input; bins event totals into daily sums
    (`arai`) handling events spanning midnight via weighted allocation,
    then builds the flux series directly.
- **`ProcessMeteoTsteps(state)`** — per timestep update of
  `graidt / nraidt / aintcdt` (and `dtEventRain`) for rain-intensity mode,
  or of `ptra / peva / graidt / nraidt` from the pre-loaded
  `tpot/epot/grain/nrain` arrays for detailed-meteo mode; ends with a
  call to `reduceva_dt(nraida, state)`.
- **`ETSine(state)`** — distributes the daily potentials `ptraday`,
  `pevaday` over a half-sine wave between sunrise and sunset. Computes
  the photoperiod-day length via `astro` at day start, the sunrise/sunset
  fractions `tsunrise_atm` / `tsunset_atm`, and integrates the wave
  fraction over `[t − dt, t]` (with careful handling of timesteps that
  span sunrise or sunset). Closes with `reduceva_dt(nraida, state)`.

### Modernisation notes
- `data pi/3.14159265d0/` should be a `real(real64), parameter :: pi`
  (or use `acos(-1.0_real64)` / the constant from
  `atmosphere_constants`).
- The four-branch `if … elseif …` chain in `ETSine` for sine-wave
  fraction is correct but unreadable — extract a pure
  `sine_wave_fraction(t0, t1, sunrise, sunset)` function.
- The `use variables` block still imports `rad, daylp, difpp, atmtr,
  dsinbe, tsunrise_atm, tsunset_atm, lat, finterception, dtEventRain` as
  DEFERRED globals — these are the remaining migration items.
- `int(t1900)` truncation for the day fraction is fragile near
  midnight; `floor` is more defensible.

---

## `et.f90`

Evapotranspiration and soil-evaporation reduction.

### Public signatures
- **`pm_inputs_t` / `pm_outputs_t`** — derived types bundling
  Penman–Monteith inputs / outputs (added in commit `4b2ced1`,
  consumed in `c29408b`).

- **`PenMon_calc(inputs, outputs)`** — `pure subroutine`, no I/O.
  Computes atmospheric properties (`λ`, `Δ`, `eₐ`, `eₐ−eₓ`, `γ`, `ρ`,
  `cₚ`, altitude-corrected pressure), aerodynamic resistances for
  crop / wet-crop / bare-soil using a log-wind profile, net short- and
  long-wave radiation plus soil heat flux, standard Penman–Monteith for
  the three potentials `es0 / et0 / ew0`, and (when `swdivide = 1`)
  PMdirect partitioning into `Edirect`, `Tdirect`, `Tdirectwet`,
  `Edirectpond`. Returns `warning_code` (0=ok, 1=polar 0 hrs,
  2=polar 24 hrs).

- **`PenMon(inputs, outputs, logf, swscre)`** — non-pure wrapper around
  `PenMon_calc` that calls `astro` for the daily branch and forwards
  `warning_code` to `warn` for screen / log output.

- **`reduceva_daily(nrai, state)`** — daily branch of soil-evaporation
  reduction (post-Phase-C; replaces `reduceva(1, …)`).
- **`reduceva_dt(nrai, state)`** — sub-daily branch (replaces
  `reduceva(2, …)`).

  Both call into `reduceva_apply`, which selects between Black
  (`swredu=1`) and Boesten–Stroosnijder (`swredu=2`) reduction. When
  ponding exceeds `POND_THRESHOLD_CM`, falls through with `empreva = peva`
  and resets `ldwet / spev / saev`. Otherwise dispatches to:
  - **`black_reduction(nrai, nird, peva, cofred, rsigni, ldwet, empreva, dt, fldaystart, task)`** —
    `Eₐ = β·(√t − √(t−Δt))/Δt`, with day-start reset of `ldwet`
    when input exceeds `rsigni`.
  - **`boesten_stroosnijder_reduction(nrai, nird, peva, cofred, spev, saev, empreva, dt)`** —
    accumulates the unmet evaporation deficit `spev` and the
    cumulative actual `saev`, using `Eₐ = E_p` (stage 1) below `β²`
    and `Eₐ = β·√spev` (stage 2) above.

### Modernisation notes
- `reduceva` still pulls `swredu, cofred, nird, rsigni` from
  `variables` (DEFERRED) — these are the last globals blocking a
  pure `et_mod` (config → `state%cfg%meteo%evaporation`, irrigation
  arc owns `nird`).
- Use `error stop` (or the project's `fatalerr_collected`) consistently
  for invalid tasks rather than mixing the two.

---

## `runoff.f90` — SCS Curve-Number runoff (`runoff_mod`)

### Public signatures
- **`cn_init(state)`** — initialisation. Validates `CNtimTAB`
  monotonicity, locates the soil-mesh node spanning the 0–10 cm
  reference layer (`nod10_cn`, `z10_cn`), and pre-computes `ThetaRef`
  (reference water content for moisture-correction modes 1 and 2).
- **`cn_step(state)`** — dynamic step. Advances the CN time-series
  pointer (`icn_atm`), derives `CNdry` / `CNwet` from `CNref`, applies
  moisture correction (`wc_cor`), computes maximum retention `S` and
  initial abstraction `Ia = INITIAL_ABSTRACTION_RATIO · S`, then sets
  `Runoff_CN` from `(P − Iₐ)² / (P − Iₐ + S)` using `nraidt + melt` as
  input depth.

### Modernisation notes
- Phase B consumed `DEPTH_10CM_CM`, `INITIAL_ABSTRACTION_RATIO`,
  `H_FIELD_CAPACITY_CM`, `H_WILTING_POINT_CM` from
  `atmosphere_constants`.
- Still imports CN switches from `variables` (DEFERRED).

---

## `snow.f90`

Snow water equivalent (SWE) accounting.

### Public signatures
- **`snow_init(state)`** — initialisation. If `swinco = 3` (warm restart)
  copy `state%atmosphere%ssnow` to `snowinco`; otherwise seed `ssnow`
  from `snowinco`.
- **`snow_step(state)`** — dynamic step:
  1. Optional sublimation: when there is a pack and `swsublim = 0`,
     drains `peva` into `subl` and zeroes `empreva`/`peva`.
  2. Above-freezing soil surface with `gsnow > 0` short-circuits to
     all fresh snow becoming melt.
  3. Otherwise: temperature-index melt
     `smelt = snowcoef · (Tav − SNOW_TEMPERATURE_C)`, augmented by
     rain-on-snow
     `smeltr = snrai · SPECIFIC_HEAT_WATER · (Tav − SNOW_TEMPERATURE_C) / LATENT_HEAT_MELTING`.
  4. Updates SWE: `ssnow = ssnow + gsnow − subl − melt − slw`, with
     liquid-water storage capped at
     `SNOW_LIQUID_WATER_FRACTION · (ssnow + slw)`.
  5. Handles the snow-deficit case (`ssnow < 0`) by scaling down
     `melt` and `subl` proportionally.
  6. Accumulates intermediate and cumulative totals.

### Modernisation notes
- Phase B consumed `SPECIFIC_HEAT_WATER`, `LATENT_HEAT_MELTING`,
  `SNOW_TEMPERATURE_C`, `SNOW_LIQUID_WATER_FRACTION`,
  `SOIL_SURFACE_FREEZE_THRESHOLD_C` from `atmosphere_constants`.
- Phase C replaced the `case (1)/case (2)` task selector with the
  named `snow_init` / `snow_step` pair.
- Still uses **bare `use Variables`** — `swinco`, `swsublim`,
  `swetsine`, `snowcoef`, `ISsnowBeg` are wide-import globals
  (config-side migration pending).

---

## `precipitation.f90`

### Public signatures
- **`PartitionPrecipitation(swmetdetail, swsnow, tav, TePrRain, TePrSnow, ssnow, nmetdetail, arain, grai, gsnow, snrai, fprecnosnow, restint, state)`**
  - Daily branch (`swmetdetail = 0`): converts `grai` from mm to cm and
    applies a three-region temperature index:
    - `Tav > TePrRain` → all rain,
    - `Tav < TePrSnow` → all snow,
    - in between → linear `f_snow = (TePrRain − Tav)/(TePrRain − TePrSnow)`.
    Counts `snrai` (rain falling on existing snowpack) only when
    `ssnow > 1e-6`, and writes `fprecnosnow = 1 − (gsnow+snrai)/grai`.
  - Detailed branch (`swmetdetail = 1`): sums `arain(i)` into `grai`,
    forces `gsnow = ssnow = snrai = 0`, `fprecnosnow = 1` (snow
    disabled).

### Modernisation notes
- Best-shaped routine in the subsystem: explicit interface, named
  intents, no globals, real64 via `iso_fortran_env`. **Pattern to
  copy.**
- The dual-write block (every local also assigning the
  `state%atmosphere%X` field) is a transitional artefact of the
  state-rescue arc — once all callers read from `state` directly,
  drop the locals and the duplicate writes.
- Move the `1.0d-6` snowpack-presence threshold into a named
  constant in `atmosphere_constants`.

---

## `interception.f90`

Four-method interception calculator plus the rain/irrigation splitter.
Intrinsics swept to generics in Phase F.

### Public signatures
- **`VonHHBraden(aintc, grai_in, state)`** — Von Hoyningen-Hune &
  Braden empirical exponential.

- **`Gash(aintc, grai_in, state)`** — analytical forest model (Gash
  1995).

- **`ruttervw(gctp, aintc, eintc, state)`** — wrapper that adapts the
  SWAP arguments (double precision) to the MetaSWAP **`msw1eic`**
  routine (single-precision arrays). Sole consumer of
  `state%atmosphere%sicact`.

- **`msw1eic(nuk, ibd, dc, dtsw, csk, vxick, fecmnk, ETw0, Pgdtsw, Sic, Sicolddtsw, Picdtsw, Eicdtsw, tcap, beta, zeta, fricdtsw, ib)`**
  — MetaSWAP "Sparse Gash" linear-ODE solver imported from Alterra's
  MSW1EIC.FOR (verbatim, OpenMP-parallel over SVATs). See open theme #1
  below.

- **`DivIntercep(aintc, state)`** — partitions today's interception
  between rain and sprinkler irrigation and writes the net depths
  `nraida` (to `state%atmosphere`) and `nird` (still legacy global).

### Modernisation notes
- See cross-cutting theme #1 (msw1eic / ruttervw quarantine).
- `DivIntercep` still writes the legacy `nird` global (irrigation-arc
  DEFERRED).
- `Gash` interpolates five separate AFGEN tables every call; consider
  a single `gash_tables_at(t)` helper returning a small derived type.

---

## Cross-cutting modernisation themes

**Resolved by GR-ATM-CLEAN (2026-05-22).** The following items from the
previous open list shipped and have been retired:

- ~~#1 Kill the shared-state scratch bag (MeteoVars)~~ — **DONE 2026-05-22
  (GR-ATM-CLEAN Phase D)**. Scratch scalars (`hum`, `etr`, `wfrac`, …)
  moved to locals/state%atmosphere; arrays `arain`/`awind` promoted to
  `state%atmosphere%arain_subdaily` / `awind_subdaily`. `meteo_vars.f90`
  deleted.
- ~~#2 Split meteoday.f90 and stop bundling modules~~ — **DONE 2026-05-22
  (GR-ATM-CLEAN Phase A)**. Three-modules-in-one file replaced with
  `meteo_vars.f90` (then deleted in Phase D), `runoff.f90`, `meteo_io.f90`,
  `meteo_orchestrator.f90`.
- ~~#3 Replace task-selector subroutines with named init/step pairs~~ —
  **DONE 2026-05-22 (GR-ATM-CLEAN Phase C)**. `snow(task)`,
  `CNmethod(Itask)` and `reduceva(task)` replaced with `snow_init`/
  `snow_step`, `cn_init`/`cn_step`, `reduceva_daily`/`reduceva_dt`.
- ~~#5 Promote `PenMon` to a derived-type interface~~ — **DONE 2026-05-22
  (commit `c29408b`)**. `pm_inputs_t` / `pm_outputs_t` added; the daily
  call site is no longer a 40-arg soup.
- ~~#6 Separate daily vs sub-daily orchestrators~~ — **DONE 2026-05-22
  (GR-ATM-CLEAN Phase E)**. `ProcessMeteoDay` now dispatches to
  `process_meteo_day_daily` / `process_meteo_day_subdaily`; four named
  private helpers extracted.
- ~~#7 Smaller efficiency wins worth bundling~~ — **DONE 2026-05-22
  (GR-ATM-CLEAN Phases B + F)**. `atmosphere_constants_mod` populated and
  consumed in `snow`/`runoff`/`et`; F77 intrinsics
  (`dexp`/`dlog`/`dmin1`/`dble`) swept to F90 generics across `et.f90`,
  `interception.f90`, `meteo_io.f90`, `meteodt.f90`.

**Still open:**

1. **Quarantine / replace `msw1eic` + `ruttervw`.** `msw1eic` is
   dead-on-arrival modern Fortran: single-precision fixed-size `(1)`
   arrays, `9199 format` numeric labels, `write(*,…) stop`, copyright
   header forbidding modification. Either replace the call with the
   SWAP-native double-precision linear-ODE solution (the algebra is in
   the doc comment), and delete the routine; or extract the kernel into
   a pure scalar function and keep the OpenMP loop only when there is a
   real SVAT array. `ruttervw` exists *only* to bridge SWAP↔MetaSWAP
   precision and dies the moment the replacement lands. Separate spec
   needed.
2. **Drop `real(8)`** in favour of `real(real64)` from `iso_fortran_env`
   (project-wide; not just atmosphere).
3. **Eliminate the `variables` import chain** — every routine still has
   DEFERRED globals (config switches `swredu/cofred/nird/rsigni`,
   irrigation `nird`, logging globals). Tracked in the wider
   globals-retirement arcs (GR-ATM Phase C3, irrigation arc, Arc 9
   logging).
4. **Pure procedures by default.** `PenMon_calc`,
   `PartitionPrecipitation`, and the two reduction kernels are already
   pure. Aim to make the interception schemes pure once `state` reads
   are replaced with explicit arguments or a small inputs derived type.
5. **Move legacy-API shims into a separate translation unit** (e.g. the
   `PenMon` wrapper) so the modern core can be audited and tested
   independently.
