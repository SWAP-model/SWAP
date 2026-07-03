---
title: "T2-A — compartment exchange records (design spec)"
date: 2026-07-03
status: design
tags: [tier2, exchange-records, component-isolation, coupling, wofost]
parent: dev-docs/adr/0051-swap-core-vs-coupled-component-boundary.md
---

# T2-A — compartment exchange records (design)

Implements ADR 0051's component boundary: turn each cross-compartment coupling
surface into an **explicit, typed, directional exchange record** — the single
channel a compartment talks through, instead of reaching into a sibling's
`state%X` fields. This is the foundation for swapping a compartment's
implementation (the near-term target: replace the in-house WOFOST with an
external crop model), and it removes the "everyone touches everything" residue
that the state migration left behind.

## Decisions (locked)

1. **Design the full set up front** (this doc), build incrementally.
2. **Records are fields on `swap_state_t`** — a new `state%exchange` aggregate of
   per-surface records, threaded like every other state record (consistent with
   today's arg-threading; can evolve to first-class MF6-style objects later).
3. **One combined "water-relations" exchange** for the crop↔atmosphere↔soil-water
   transpiration triangle (not strict bilateral) — because relative transpiration
   inherently combines an atmosphere quantity (`ptra`) and a soil quantity
   (`tra`).
4. **Copy semantics** — a record holds the crossing *values* (copied at the
   boundary), not pointers/aliases into sibling state. This is what makes the
   surface a real contract and makes an *external* component (which cannot share
   memory) a drop-in producer/consumer.
5. **Byte-identical** — each surface is cut over by routing its existing
   reach-across line(s) through the record; the value is a pure copy, so output
   is unchanged. `check-fast` + `check-bindings` green per surface.

## Principle: one record per coupling surface

Soil-water is the hub (SWAP *is* a soil-water model), so most surfaces are
soil-water ↔ X. One record per surface keeps each contract minimal — the
anti-pattern (rejected) is a single `swap_exchange_t` holding every crossing,
which would recreate the all-depend-on-all coupling.

The bottom-boundary surface **already exists** as a de-facto exchange record: the
XMI arrays `gwl` / `qbot_volume` / `storage_coef` (soil-water ↔ external
groundwater). The internal records apply the same shape to the couplings that are
currently implicit.

## Ownership & threading

```fortran
!> Aggregate of the cross-compartment exchange records (T2-A). One field per
!> coupling surface. A compartment reads/writes ONLY its side of the relevant
!> record; it never touches a sibling compartment's state directly.
type :: swap_exchange_t
   type(crop_water_exchange_t) :: crop_water   ! crop <-> atmosphere+soil-water (ET + sink + stress)
   type(heat_soil_exchange_t)  :: heat_soil    ! temperature <-> conductivity/water content
   type(drain_soil_exchange_t) :: drain_soil   ! drainage sink <-> heads/GWL
   type(atmos_soil_exchange_t) :: atmos_soil   ! top boundary: net flux <-> realized evaporation
   ! DEFERRED (do not fit the clean bilateral model — separate later pass):
   !   solute        — a many-to-one CONSUMER (reads water flux + temperature +
   !                   drainage flux + irrigation/rain conc from 4 compartments);
   !                   a 'transport-drivers' bundle, not a bilateral pair.
   !   surfacewater  — pond/runoff genuinely bidirectional with the top boundary.
end type
! swap_state_t gains: type(swap_exchange_t) :: exchange
```

**Excluded (not coupling — output accounting).** The `atmo%intr%*` / `atmo%cumu%*`
fields written by `waterbalance.f90:371-413` are per-period result accumulation,
not dynamic inter-compartment coupling. They stay in the results/water-balance
layer and get **no** exchange record.

## Update protocol (cadence)

Each record splits into `from_<producer>` and `to_<consumer>` groups. The
producer fills its group at the top of its step; the consumer reads after. The
crop↔water triangle is the trickiest because it spans three cadences:

- **crop day start:** crop fills `pot_et`, `lai`, `crop_height`, `rooting_depth`,
  `root_density` (stable across the day's water substeps).
- **atmosphere:** fills `pot_transp` (`ptra`) from `pot_et`.
- **soil-water substeps:** read `root_density` × `pot_transp` for the sink;
  accumulate `act_transp` over substeps.
- **crop day end / next day:** reads `rel_transp = act_transp/pot_transp`.

So `from_crop` fields are write-once-per-day; `to_crop` fields finalize at day
boundary. The record encodes this so the ordering is explicit, not conventional.

## The record set

> **Implementation finding (2026-07-04).** The crop↔*soil-water* coupling routes
> cleanly and is **done**: the root distribution (`rd`/`noddrz`/`cumdens`) and the
> transpiration feedback (`ptra`/`tra`→`reltr`) go through the record; the Feddes
> reduction thresholds stay native sink config. But the crop↔*atmosphere* ET
> coupling is **not** a clean handoff: `et0`/`es0`/`ew0` are computed *inside* the
> ET module (`meteo_orchestrator.f90:368-383`) from the crop factor × reference
> ET, stored on crop state — they are atmosphere outputs, not crop-provided. The
> genuine crop→atmosphere dynamic crossings are the canopy props (`lai`, `ch`);
> the potential-ET quantities are intra-atmosphere. Routing the ET interface
> therefore needs its own focused sub-pass (untangling the in-module ET
> computation), and `et0`/`es0`/`ew0` are **removed** from the from-crop list
> below. Deferred as `crop_water` sub-part 2.

### 1. `crop_water_exchange_t` — full crop↔SWAP interface (~18 dynamic fields)

The complete boundary an external crop model must satisfy: it feeds SWAP's ET +
interception *and* the root sink, and reads back the water-stress signal. Crop
config-*constants* (`albedo`, `rsc`, `rsw`, `kdir`, `kdif`, `cf`, the Feddes
`hlim*`/`adcr*` thresholds, `alphacrit`) stay **config-passed** — they don't
change per day, so they aren't in the dynamic exchange.

Seam lines routed through the record: `meteo_orchestrator.f90:368-472` (ET reads
`crop%et0/es0/ew0/lai/gird/ch/fco2tra`), `rootextraction.f90:71,95,122` (sink
reads `crop%cumdens/rd/noddrz`; writes `soil%qrot`), `cropgrowth.f90:292-293`
(crop reads `soil%hroot/hleaf`), `cropwofost_runtime.f90:656` (`reltr` from
`soil%tra`/`atmo%ptra`).

```fortran
type :: crop_water_exchange_t
   ! --- from crop (write-once per day) ---
   real(real64)              :: pot_et_crop  = 0.0  !! et0  (cm/d) potential transpiration driver
   real(real64)              :: pot_et_soil  = 0.0  !! es0  (cm/d) potential soil evaporation
   real(real64)              :: pot_et_wet   = 0.0  !! ew0  (cm/d) potential wet-canopy evaporation
   real(real64)              :: lai          = 0.0  !! -> interception + extinction
   real(real64)              :: crop_height  = 0.0  !! ch (cm)
   real(real64)              :: interc_demand = 0.0 !! gird
   real(real64)              :: co2_transp_fac = 1.0 !! fco2tra
   real(real64)              :: rooting_depth = 0.0 !! rd (cm)
   integer                   :: root_nodes   = 0    !! noddrz
   real(real64), allocatable :: root_density(:)     !! cumdens -> sink distribution
   ! --- from atmosphere (set by the ET step) ---
   real(real64)              :: pot_transp   = 0.0  !! ptra (cm/d)
   real(real64)              :: pot_evap     = 0.0  !! peva (cm/d)
   ! --- from soil-water (finalize at day end) ---
   real(real64), allocatable :: root_flux(:)        !! qrot per node (cm/d)
   real(real64)              :: act_transp   = 0.0  !! tra  (cm/d)
   real(real64)              :: root_head    = 0.0  !! hroot — root-zone matric potential
   real(real64)              :: leaf_head    = 0.0  !! hleaf — leaf water potential
   real(real64)              :: rel_transp   = 1.0  !! reltr = tra/ptra (growth-reduction signal)
end type
```

### 2. `heat_soil_exchange_t` — clean (~5 fields)

`temperature.f90`/`frozencond.f90` read soil water content; `boundtop`/
`soilhydraulics` read soil temperature + frozen fraction.

```fortran
type :: heat_soil_exchange_t
   ! from heat -> soil-water (K adjustment, frozen soil)
   real(real64), allocatable :: rfcp(:)     !! reduced-frozen-fraction on K
   real(real64), allocatable :: tsoil(:)    !! soil temperature (degC)
   ! from soil-water -> heat (de Vries thermal props / transport)
   real(real64), allocatable :: theta(:)    !! current water content
   real(real64), allocatable :: theta_prev(:) !! thetm1 (previous step)
end type
```
*(`thetas` (saturated) is init-only — excluded.)*

### 3. `drain_soil_exchange_t` — clean (~8 fields)

`drainage.f90`/`divdra.f90` read soil GWL / heads / K; soil-water reads the
drainage flux.

```fortran
type :: drain_soil_exchange_t
   ! from drainage -> soil-water
   integer                   :: n_levels = 0     !! nrlevs
   real(real64), allocatable :: qdra(:,:)        !! drainage flux per level/node
   real(real64), allocatable :: zbotdr(:)        !! level bottom depths
   ! from soil-water -> drainage
   real(real64)              :: gwl = 0.0         !! groundwater level
   real(real64), allocatable :: h(:)             !! pressure-head profile
   real(real64)              :: pond = 0.0        !! surface ponding
   ! K + anisotropy are crop/soil config-constants -> stay config-passed
end type
```

### 4. `atmos_soil_exchange_t` — top boundary (~10 dynamic fields)

`boundtop.f90` reads the atmospheric net flux; soil-water returns realized
evaporation + ponding. **Excludes** the `atmo%cumu`/`intr` accounting.

```fortran
type :: atmos_soil_exchange_t
   ! from atmosphere -> soil-water (top boundary flux)
   real(real64) :: net_rain   = 0.0  !! nraidt
   real(real64) :: irrig      = 0.0  !! nird
   real(real64) :: snowmelt   = 0.0  !! melt
   real(real64) :: pot_evap   = 0.0  !! peva
   real(real64) :: emp_evap   = 0.0  !! empreva
   integer      :: evap_reduce_method = 0 !! swredu
   ! from soil-water -> atmosphere (feedback)
   real(real64) :: act_soil_evap = 0.0 !! reva (realized soil evaporation)
   real(real64) :: pond          = 0.0 !! affects the peva adjustment
end type
```

### Deferred (separate later pass)
- **solute** — a transport-*drivers* bundle (water flux `q`, `theta`/`thetm1`,
  `qrot`, `qtop`, `pond`, `gwl` from soil; `tsoil` from heat; `nird`/`nraidt`
  from atmosphere). It is a downstream consumer, not a bilateral pair; a
  driver-bundle record is a clean design but lower value now.
- **surfacewater** — `pond`/`hsurf`/`runots`/`H0max`/`FlRunoff` are bidirectional
  with the top boundary within a substep; needs its own protocol design.

## Byte-identity strategy (per surface)

1. Populate the record's `from`/`to` fields at the boundary — a pure copy of the
   value the reach-across line currently reads.
2. Repoint the reach-across line(s) to read the record field instead of the
   sibling state.
3. `check-fast` byte-identical + `check-bindings` green. Because the record value
   is an exact copy inserted just before the existing read, output is bit-for-bit
   unchanged. No FP reordering (pure pass-through).
4. After cutover, the sibling-state field is no longer read cross-compartment
   (it may still be the record's *source* on the owner side).

## Implementation sequencing

Records are independent surfaces, so order by value/clarity. Each is a
byte-identical arc with its own commit; the set shares one architecture ADR
(0053) plus this design doc.

1. **crop_water** *(pattern-setter + WOFOST-enabler)* — introduces
   `swap_exchange_t` on `swap_state_t`, the per-day/per-substep protocol, and the
   ADR. Largest of the clean records (~18 fields), but the highest-value boundary.
   Route the ET / rootextraction / reltr seam lines through it.
2. **heat_soil** — small, clean (~5 fields); good confidence-builder right after
   the pattern is set.
3. **drain_soil** — small, clean (~8 fields).
4. **atmos_soil** — top-boundary flux (~10 fields), excluding the accounting.
5. *(deferred)* **solute** driver-bundle, **surfacewater** bidirectional loop —
   separate design pass after the clean four land.

**Rationale for crop-first despite its size:** it is the only surface that
*unlocks a capability* (external-WOFOST plug), and building the `swap_exchange_t`
machinery + protocol against the richest record de-risks the smaller ones. If a
smaller warm-up is preferred, `heat_soil` could go first purely to validate the
mechanics — flag if you want that.

## Payoff: the external-WOFOST plug

Once `crop_water` is the *only* crop↔rest channel, `crop_step()` is a swappable
producer of `from_crop`: native table crop, in-house WOFOST, or **external WOFOST
over BMI/XMI** (SWAP hands it `rel_transp`; it returns `lai`/`rooting_depth`/…).
Structurally identical to the SWAP↔MODFLOW 6 coupling already in place — same
transport (BMI/XMI), daily cadence, the record as the contract.

## Open design notes

- **Stress physics stays water-side.** Drought/salinity/oxygen reduction is
  computed in `rootextraction`; the crop only receives `act_transp`/`rel_transp`.
  Individual stress factors are a cheap future addition to `to_crop` if an
  external crop model wants them.
- **`root_density` is an allocatable array in the record** — the one non-scalar;
  sized to the mesh. Copy cost is per-day (not per-substep), negligible.
- **Evolution to first-class objects** (MF6-style, off `swap_state_t`) remains
  open; the `state%exchange` aggregate is the incremental first form.
</content>
