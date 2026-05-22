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

### Module layout

| File | Module(s) | Role |
|------|-----------|------|
| `atmosphere_constants.f90` | *(none yet — header only)* | Placeholder for shared physical constants |
| `meteoday.f90` | `MeteoVars`, `runoff_mod`, `meteo_process_mod`, `meteo_mod` | Daily meteo ingestion, CN runoff, daily ET / interception orchestration |
| `meteodt.f90` | `meteodt_mod` | Sub-daily meteo orchestration (rain events, time-step fluxes, sine-wave ET) |
| `et.f90` | `et_mod` | Penman–Monteith ET and soil-evaporation reduction |
| `snow.f90` | `snow_mod` | Snow accumulation, sublimation, melt, liquid-water storage |
| `precipitation.f90` | `precipitation_mod` | Rain / snow partitioning |
| `interception.f90` | `interception_mod` | Four interception schemes + interception splitter |

### Dependencies

Imports from
- `src/core/` — `swap_state_mod`, `swap_config_mod`, `swap_constants`, `swap_array_dimensions`, `variables` (legacy bridge — being retired),
- `src/utils/` — `array_utils` (`afgen`), `soilhydraulics_utils` (`watcon`), `error_mod`.

Has no downstream physics dependencies — every consumer reads atmosphere
outputs through `state%atmosphere`.

---

## `atmosphere_constants.f90`

Header-only file. Declares the *intent* to host named constants for the
atmosphere subsystem (replacing magic numbers in the calculation code) but
contains no `module` block yet.

### Public signatures
- *(none)*

### Modernisation notes
- File ships with a doc comment but no actual `module … end module` block;
  Meson currently compiles it as an empty translation unit. Either create
  `module atmosphere_constants` with `real(real64), parameter :: …` declarations
  (e.g. `cwat`, `lm`, `ts`, `DEPTH_10CM`, `IA_RATIO`, `ALBEDO_PONDING`, the
  Karman constant, etc., currently scattered across `snow.f90`, `et.f90`,
  `meteoday.f90`) or delete the stub.
- Use `real64` from `iso_fortran_env` rather than the inherited `real(8)` kind.

---

## `meteoday.f90`

Bundles **three** modules in a single file (a candidate for splitting). It
covers everything that runs at the daily timestep: reading the meteo file,
computing reference ET and interception, partitioning into soil/canopy demand,
and applying the SCS Curve-Number runoff method.

### `module MeteoVars` — scratch / temporary state

Plain data module of module-level scalars and small fixed-size arrays
(`arain(96)`, `awind(96)`) used to pass intermediate meteo values between
`ReadMeteoDay`, `ProcessMeteoDay`, and `runoff_mod`. No procedures.

### `module runoff_mod` — SCS Curve-Number runoff

Public:

- **`CNmethod(Itask, state)`** — two-mode SCS-CN driver.
  - `Itask = 1`: initialisation. Validates `CNtimTAB` monotonicity, locates
    the soil-mesh node spanning the 0–10 cm reference layer
    (`nod10_cn`, `z10_cn`), and pre-computes `ThetaRef` (reference water
    content for the moisture-correction modes 1 and 2).
  - `Itask = 2`: dynamic step. Advances the CN time-series pointer
    (`icn_atm`), derives `CNdry` / `CNwet` from `CNref`, applies moisture
    correction (`wc_cor`), computes maximum retention `S` and initial
    abstraction `Ia = 0.2 · S`, then sets `Runoff_CN` from
    `(P − Iₐ)² / (P − Iₐ + S)` using `nraidt + melt` as input depth.

### `module meteo_process_mod` — daily meteo I/O

Public:

- **`ReadMeteoDay(state, config)`** — reads today's meteo record from the
  pre-loaded `state%atmosphere%a*` arrays. Two branches:
  - `swmetdetail = 0` (daily): copies `arad/atmn/atmx/ahum/awin/arai/aetr`
    into transient locals + state, computes 24 h and daytime mean
    temperature, saturation vapour pressure (Tetens), relative humidity,
    and mirrors the values into `state%atmosphere%out_*` for PEARL CFO
    output.
  - `swmetdetail = 1` (sub-daily): validates per-record indices and
    timestamps against the meteo file, fills `state%atmosphere%arad/ahum/atav`
    and the local `awind`, `arain` arrays.
  - Always calls `ResetMetFlx` first, then `PartitionPrecipitation` to
    split today's `grai` into rain / snow / `fprecnosnow`.

- **`ResetMetFlx(state)`** — zeros the intermediate and cumulative
  atmosphere flux cohorts when the controller raises `flZeroIntr` /
  `flZeroCumu`. Delegates to `state%atmosphere%intr%reset()` and
  `state%atmosphere%cumu%reset()` (post ADR-0033 cohort-reset pattern).

### `module meteo_mod` — daily ET / interception coordinator

Public:

- **`ProcessMeteoDay(state, config)`** — main daily orchestrator. Runs (in
  order):
  1. interception (`VonHHBraden` / `Gash` / `ruttervw` depending on
     `swinter`) and `DivIntercep`,
  2. ET₀ / EW₀ / ES₀ from `PenMon` or from supplied reference ET, with
     crop-factor switches (`swcf`, `swcfbs`),
  3. wet-canopy fraction `wfrac`,
  4. potential soil evaporation `peva` (cover-corrected, ponding-corrected,
     PMdirect when `swdivide = 1`) and potential transpiration `ptra`
     (with optional CO₂ correction `fco2tra`),
  5. for sub-daily: aggregates `tpot/epot/grain/nrain(i)` per record and
     recomputes daily mean / min / max temperature, RH and radiation,
  6. for daily: when `swrain = 0`, sets time-step rain fluxes
     (`graidt`, `nraidt`, `aintcdt`), applies `CNmethod(2, …)` if
     `swuseCN = 1`, then `reduceva(1, …)` for actual evaporation.

### Modernisation notes
- **Three modules per file** — split into `meteo_vars.f90`,
  `runoff.f90`, `meteo_io.f90`, `meteo_day.f90` so each module has its
  own translation unit and Meson dependency graph.
- `module MeteoVars` is a **shared mutable-state bag** — exactly the
  pattern the wider rescue arc is retiring. Migrate its scalars (`hum`,
  `etr`, `wfrac`, …) into `atmosphere_state_t` / local variables and
  pass `arain`/`awind` as arguments.
- Heavy use of inherited `real(8)`, `data` statements
  (e.g. `rcs/0.15d0/`), numeric loop labels (`1000 continue`), and bare
  `use variables` chains still flagged DEFERRED. Replace with
  `real64`, `parameter` constants, `do … end do`, and explicit imports.
- `ProcessMeteoDay` is ~430 lines with sentinel comment "Section N"
  blocks — each of the ten sections is a candidate `private` subroutine
  (`apply_interception`, `compute_ref_et`, `partition_peva_ptra`, etc.).
- The `dexp`/`dlog`/`dmin1`/`dble`/`int` chains are F66 holdovers; use
  the generic `exp`, `log`, `min`, `real`, `nint`.
- `1000 continue` should become `end do dayparts_loop` with a named loop.

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
  `graidt / nraidt / aintcdt` (and `dtEventRain`, the time-to-next-rain
  control) for rain-intensity mode, or of `ptra / peva / graidt / nraidt`
  from the pre-loaded `tpot/epot/grain/nrain` arrays for detailed-meteo
  mode; ends with a call to `reduceva(2, nraida, state)`.
- **`ETSine(state)`** — distributes the daily potentials `ptraday`,
  `pevaday` over a half-sine wave between sunrise and sunset. Computes
  the photoperiod-day length via `astro` at day start, the sunrise/sunset
  fractions `tsunrise_atm` / `tsunset_atm`, and integrates the wave
  fraction over `[t − dt, t]` (with careful handling of timesteps that
  span sunrise or sunset). Closes with `reduceva(2, nraida, state)`.

### Modernisation notes
- `data pi/3.14159265d0/` should be a `real(real64), parameter :: pi`
  (or use `acos(-1.0_real64)` / the constant from
  `atmosphere_constants`).
- The four-branch `if … elseif …` chain in `ETSine` for sine-wave
  fraction is correct but unreadable — extract a pure
  `sine_wave_fraction(t0, t1, sunrise, sunset)` function.
- The `use variables` block still imports `rad, daylp, difpp, atmtr,
  dsinbe, tsunrise_atm, tsunset_atm, lat, finterception, dtEventRain` as
  DEFERRED globals — these are the remaining migration items for the
  atmosphere arc (Phase C3 candidates).
- `int(t1900)` truncation for the day fraction is fragile near
  midnight; `floor` is more defensible.

---

## `et.f90`

Evapotranspiration and soil-evaporation reduction.

### Public signatures
- **`PenMon_calc(daynr, lat, alt, altw, a, b, rcs, rad, tav, hum, win, …, es0, et0, ew0, Edirect, Tdirect, Tdirectwet, Edirectpond, warning_code)`**
  — `pure subroutine`, no I/O. Computes:
  - atmospheric properties (`λ`, `Δ`, `eₐ`, `eₐ−eₓ`, `γ`, `ρ`, `cₚ`,
    altitude-corrected pressure),
  - aerodynamic resistances for crop / wet-crop / bare-soil using a
    log-wind profile, displacement `d = ⅔ h`, roughness `z₀ₘ = 0.123 h`
    and `z₀ₕ = 0.1 z₀ₘ`, with wind extrapolation to the crop level when
    measurement height < canopy height (`fmeas`, `fact` factors),
  - net short- and long-wave radiation plus soil heat flux (day/night
    split for sub-daily),
  - standard Penman–Monteith for the three potentials `es0 / et0 / ew0`,
  - optional PMdirect (`swdivide = 1`) partitioning into `Edirect`,
    `Tdirect`, `Tdirectwet`, `Edirectpond` using vegetation cover
    `vcover = 1 − exp(−kdif·kdir·LAI)` and effective LAI `LAI/(0.3·LAI+1.2)`.
  - Returns `warning_code` (0=ok, 1=polar 0 hrs, 2=polar 24 hrs).

- **`PenMon(logf, swscre, …)`** — non-pure wrapper around `PenMon_calc`
  that calls `astro` for the daily branch and forwards `warning_code`
  to `warn` for screen / log output. Preserved for the legacy call
  signature.

- **`reduceva(task, nrai, state)`** — selects between Black (`swredu=1`)
  and Boesten–Stroosnijder (`swredu=2`) reduction. When ponding
  exceeds 1e-10 cm, falls through with `empreva = peva` and resets
  `ldwet / spev / saev`. Otherwise dispatches to:
  - **`black_reduction(nrai, nird, peva, cofred, rsigni, ldwet, empreva, dt, fldaystart, task)`** —
    `Eₐ = β·(√t − √(t−Δt))/Δt`, with day-start reset of `ldwet`
    when input exceeds `rsigni`. Two task modes (daily / sub-daily).
  - **`boesten_stroosnijder_reduction(nrai, nird, peva, cofred, spev, saev, empreva, dt)`** —
    accumulates the unmet evaporation deficit `spev` and the
    cumulative actual `saev`, using `Eₐ = E_p` (stage 1) below `β²`
    and `Eₐ = β·√spev` (stage 2) above. Reduces `saev` when water
    input exceeds `E_p`.

### Modernisation notes
- `PenMon_calc` is the **template** the rest of the subsystem should
  follow: explicit-interface `pure` subroutine, `intent` on every
  argument, no globals, full doc strings. Use it as a reference when
  rewriting the other physics routines.
- The wrapper `PenMon` still owns the legacy ordering (~40 positional
  arguments). Promote it to a derived type (`pm_inputs_t`,
  `pm_outputs_t`) and have all callers pass that — the current
  call site in `ProcessMeteoDay` is a 9-line argument soup.
- `reduceva` still pulls `swredu, cofred, nird, rsigni` from
  `variables` (DEFERRED) — these are the last globals blocking a
  pure `et_mod` (config → `state%cfg%meteo%evaporation`, irrigation
  arc owns `nird`).
- `POND_THRESHOLD = 1.0d-10` is duplicated in `ProcessMeteoDay`;
  move to `atmosphere_constants`.
- Use `error stop` (or the project's `fatalerr_collected`) consistently
  for invalid tasks rather than mixing the two.

---

## `snow.f90`

Snow water equivalent (SWE) accounting.

### Public signatures
- **`snow(task, state)`** — two tasks:
  - `task = 1`: initialisation. If `swinco = 3` (warm restart) copy
    `state%atmosphere%ssnow` to `snowinco`; otherwise seed `ssnow`
    from `snowinco`.
  - `task = 2`: dynamic step.
    1. Optional sublimation: when there is a pack and `swsublim = 0`,
       drains `peva` into `subl` and zeroes `empreva`/`peva`.
    2. Above-freezing soil surface with `gsnow > 0` short-circuits to
       all fresh snow becoming melt.
    3. Otherwise: temperature-index melt
       `smelt = snowcoef · (Tav − 0)`, augmented by rain-on-snow
       `smeltr = snrai · cwat · (Tav − 0) / Lm`, with negative `smelt`
       allowed to partly cancel `smeltr` (`melt = max(0, smelt+smeltr)`).
    4. Updates SWE: `ssnow = ssnow + gsnow − subl − melt − slw`, with
       liquid-water storage capped at `0.07 · (ssnow + slw)` and excess
       routed to melt via drainage `qlw`.
    5. Handles the snow-deficit case (`ssnow < 0`) by scaling down both
       `melt` and `subl` proportionally.
    6. Accumulates intermediate (`igsnow / isubl / isnrai`) and
       cumulative (`cgsnow / csubl / csnrai / cmelt`) totals.

### Modernisation notes
- Physical constants `cwat`, `lm`, `ts`, the 0.07 liquid-water cap,
  and the 0.5 °C soil-surface threshold should live in
  `atmosphere_constants`.
- Still uses **bare `use Variables`** — `swinco`, `swsublim`,
  `swetsine`, `snowcoef`, `ISsnowBeg` are wide-import globals. These
  need to move to `state%cfg%meteo%snow` (config side) and
  `state%atmosphere` (runtime side).
- `tsoil_surf = state%heat%tsoil(1)` is the only mesh-coupled read;
  good — keep it explicit.
- The `case (1)/case (2)` task selector pattern is shared with
  `reduceva` and `CNmethod`; consider replacing with two named
  procedures (`snow_init`, `snow_step`) for clearer call sites.

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
    disabled). The forced `ssnow = 0` is flagged in-source as a
    refactor candidate.

### Modernisation notes
- Best-shaped routine in the subsystem: explicit interface, named
  intents, no globals, real64 via `iso_fortran_env`. **Pattern to
  copy.**
- The dual-write block (every local also assigning the
  `state%atmosphere%X` field) is a transitional artefact of the
  state-rescue arc — once all callers read from `state` directly,
  drop the locals and the duplicate writes.
- Move the `1.0d-6` snowpack-presence threshold into a named
  constant.

---

## `interception.f90`

Four-method interception calculator plus the rain/irrigation splitter.

### Public signatures
- **`VonHHBraden(aintc, grai_in, state)`** — Von Hoyningen-Hune &
  Braden empirical exponential. Computes a soil-cover coefficient
  `cofbb = 1 − exp(−kdif·kdir·LAI)` and
  `aintc = cofab·LAI·(1 − 1/(1 + rpd·cofbb/(cofab·LAI)))·0.1`.
  `rpd` includes irrigation when `isua = 0`.

- **`Gash(aintc, grai_in, state)`** — analytical forest model (Gash
  1995). Reads time-varying `pfree`, `pstem`, `scanopy`, `avprec`,
  `avevap` from the `state%atmosphere` AFGEN tables, computes the
  rainfall depth to saturate the canopy
  `Psat = −avprec·scanopy/avevap · ln(1 − avevap/avprec)`, then
  branches on whether today's `grai` exceeds `Psat`.

- **`ruttervw(gctp, aintc, eintc, state)`** — wrapper that adapts the
  SWAP arguments (double precision) to the MetaSWAP **`msw1eic`**
  routine (single-precision arrays). Sole consumer of `state%atmosphere%sicact`.

- **`msw1eic(nuk, ibd, dc, dtsw, csk, vxick, fecmnk, ETw0, Pgdtsw, Sic, Sicolddtsw, Picdtsw, Eicdtsw, tcap, beta, zeta, fricdtsw, ib)`**
  — MetaSWAP "Sparse Gash" linear-ODE solver. Imported from
  Alterra's MSW1EIC.FOR (verbatim, OpenMP-parallel over SVATs). Solves
  `dS/dt = csk·Pg − fecmn·ETw0 − (1−fecmn)·ETw0·S/vxick` for canopy
  storage `S`; computes `tcap` (time to fill the canopy) when the
  reservoir saturates within the timestep. Used only by `ruttervw`
  with `nuk = 1`.

- **`DivIntercep(aintc, state)`** — partitions today's interception
  between rain and sprinkler irrigation and writes the net depths
  `nraida` (to `state%atmosphere`) and `nird` (still legacy global).
  When `aintc < 1e-3` cm just removes snow / rain-on-snow from
  `grai`.

### Modernisation notes
- `msw1eic` is **dead-on-arrival modern Fortran**: single-precision
  fixed-size `(1)` arrays, `9199 format` numeric labels, `write(*,…)
  stop`, copyright header forbidding modification. Either:
  (a) replace the call with the SWAP-native double-precision linear-ODE
  solution (the algebra is in the doc comment), and delete this routine
  entirely; or (b) extract the kernel into a pure scalar function and
  keep the OpenMP loop only when there is a real SVAT array.
- `ruttervw` exists *only* to bridge SWAP↔MetaSWAP precision — kills
  itself the moment (a) above happens.
- `DivIntercep` still writes the legacy `nird` global (irrigation-arc
  DEFERRED).
- `Gash` interpolates five separate AFGEN tables every call; consider
  a single `gash_tables_at(t)` helper returning a small derived type.
- All four schemes accept `state` (only `DivIntercep`/`ruttervw`
  actually mutate it). Tag the others `intent(in)` for `state`
  (already done in `VonHHBraden`/`Gash`) — good. Make this consistent.

---

## Cross-cutting modernisation themes

These apply to every file in the subsystem and are worth listing once:

1. **Drop `real(8)`** in favour of `real(real64)` from
   `iso_fortran_env` (and the project's `real_kinds` module if/when it
   exists).
2. **Eliminate the `variables` import chain** — every routine still
   has DEFERRED globals. Track each in the existing globals-retirement
   arcs (GR-ATM Phase C3, irrigation arc, Arc 9 logging).
3. **Hoist physical constants** out of routine bodies into
   `atmosphere_constants` (Karman, latent heats, freezing point,
   pond threshold, the 0.07 liquid-water-storage fraction, AFGEN
   table sizes, …).
4. **Pure procedures by default.** `PenMon_calc`,
   `PartitionPrecipitation`, and the two reduction kernels are
   already pure. Aim to make the interception schemes pure too once
   `state` reads are replaced with explicit arguments or a small
   inputs derived type.
5. **Replace `data` statements and numeric loop labels** (`1000 continue`,
   `9199 format`) with `parameter` and named constructs.
6. **Stop bundling multiple modules per file.** `meteoday.f90`
   ships three; split for clearer dependencies and parallel builds.
7. **Switch task-selector subroutines to named procedures** where
   it improves clarity (`snow_init`/`snow_step`, `cnmethod_init`/
   `cnmethod_step`, `reduceva_daily`/`reduceva_dt`).
8. **Move legacy-API shims into a separate translation unit** (e.g.
   the `PenMon` wrapper, `ruttervw`) so the modern core can be
   audited and tested independently.
