# ADR 0051 — SWAP-core vs coupled-component boundary

**Status:** Accepted (2026-07-03) — direction-setting. Concrete detachments are
deferred to follow-on ADRs/arcs; this ADR only draws the boundary and the
rationale. First detachment (nutrients): ADR 0052.
**Relates to:** ADR 0043 (state-threaded, no hidden state — the multi-instance
invariant), ADR 0050 (one `libswap.so` + variable registry — the coupling
substrate), ADR 0040 (macropore retirement — the precedent for detaching an
unverified embedded subsystem), ADR 0009 (non-CSV outputs discontinued — the
already-parked ANIMO `.afo` interface). Reframes the Tier-2 roadmap in
[`2026-07-03-modernization-review.md`](../2026-07-03-modernization-review.md).

## Context

SWAP is currently a monolith with copies of other Wageningen models compiled in.
Verified by measure: `src/crop/` is **~11,000 lines — the single largest
subsystem**, larger than the soil-water core it serves (`src/soilwater/` ≈ 3,700).
Of that crop code, `crop/wofost/` alone is ≈ 4,600, and it embeds:

- **WOFOST** — a full detailed crop-growth model (assimilation, phenology,
  biomass partitioning, dynamic LAI). Present as SWAP crop `type = 2`.
- **Detailed grass** — a grass-growth model. Crop `type = 3`.
- **WOFOST-N / LINTUL4** — nitrogen-limited crop growth (the `.crp` cites
  *"Data from: Linutl4"*).
- **An ANIMO-derived soil-N module** — `Wofost_Soil_Declarations` /
  `Wofost_Soil_Interface`: organic-matter pools, mineralization,
  nitrification/denitrification, NH₄/NO₃, with algorithms lifted from ANIMO (the
  code says *"according to ANIMO (Groenendijk et al, 2005)"*, *"Release of
  ANIMO4.0"*). ~85 module-global variables; the dominant remaining parallelism
  blocker; dark in the regression suite; reactivated beyond stock 4.2.0 (which
  guards it in the tillage path).

These embedded models create three concrete problems:

1. **Maintenance & drift.** SWAP carries forks of models that exist as canonical
   standalones (WOFOST, ANIMO); the forks drift from upstream and must be
   maintained inside SWAP.
2. **Architecture.** The embedded models are the largest source of module-global
   state and the main blocker to in-process multi-instance / threaded ensembles
   (ADR 0043's forward note; T2-C).
3. **Composition.** The field is moving to component-based modeling over standard
   interfaces (BMI/XMI, CSDMS, MODFLOW 6 exchanges, the iMOD coupler). SWAP↔MODFLOW 6
   (working, XMI) already proves the pattern for groundwater; ANIMO coupling
   already *is* external (offline, file-based via the `.afo` Aggregated Flux
   Output). The embedded crop/nutrient models are the outliers.

The coupling substrate now exists: one `libswap.so` exposing BMI + XMI over one
ensemble backing store, a variable registry (ADR 0050). This makes drawing —
and eventually enforcing — a component boundary feasible for the first time.

**Key enabling insight (the reason this is not a capability loss).** The crop→water
coupling is *model-agnostic at a small, already-existing seam.* Whatever computes
the crop feeds the same state fields — LAI (`lai`), rooting depth (`rd`), crop
height (`ch`), potential transpiration — and reads back actual transpiration and
water stress. SWAP's **simple crop** (`type = 1`, `cropfixed`) already fills those
from prescribed input tables (`gctb`, `cftb`, `chtb`, `rdtb`) and has **zero
dependency on any WOFOST/nutrient module** (verified: no such `use` statements).
So the simple crop *is* a native, standalone "Plant"; WOFOST is one detailed
*implementation* of the same interface.

## Decision

Adopt, as the project's architectural north star, a **SWAP-core vs
coupled-component boundary**:

### SWAP-native core (kept, owned, maintained)

The vadose-zone / soil-water process engine and its native forcing + plant
interface:

- **Soil water** — Richards flow, the numerical solver, boundary conditions
  (top and bottom), drainage / surface water.
- **Soil heat** and **basic solute transport** (convection–dispersion).
- **Atmosphere forcing** — meteo, ET, interception, snow, runoff.
- **The native "Plant" (water side)** — the **simple, table-driven crop**
  (`type = 1`: LAI/green-cover, crop factor, crop height, rooting depth from
  input tables), **root water uptake** (the sink term in Richards), and **all the
  stress-reduction functions** (Feddes drought, salinity, oxygen) and root
  distribution. These apply to *any* crop and are the socket a crop model plugs
  into — they are interface, not a borrowed model.

SWAP remains, semantically, Soil-Water-Atmosphere-**Plant**: it retains a native
plant component that closes the water balance.

### Coupled components (detached, external, optional)

- **Detailed crop growth** — WOFOST (arable, `type = 2`) and detailed grass
  (`type = 3`): a pluggable *high-fidelity crop-sink backend* that computes the
  same interface variables (LAI, rooting depth, height, transpiration demand)
  from simulated biomass, and consumes the returned actual transpiration / water
  stress. Coupled at a **daily** cadence (crop growth is daily; sub-daily water
  uptake stays in the core).
- **Nutrients** — the ANIMO-derived soil-N (WSN) and WOFOST-N. Detailed nutrient
  fate belongs to the external **ANIMO** model, which consumes SWAP's *hydrology*
  (the `.afo` output), not SWAP's internal nutrient state. The in-house WSN is
  removed rather than de-globalled.

### The interface

One model-agnostic crop↔water exchange record (a Tier-2 typed exchange object):
**out** = {LAI or soil-cover, rooting depth, crop height, potential
transpiration/interception}; **back** = {actual transpiration, water-uptake
stress}. The native simple crop is the default implementation; WOFOST/grass plug
into the same record via BMI/XMI. Nutrient coupling is SWAP-hydrology → ANIMO,
re-using the `.afo` concept (ADR 0009 parked it; re-enable if/when needed).

## Consequences

- **SWAP shrinks toward its core competency** — a maintainable vadose-zone engine
  — and stops carrying forks of WOFOST and ANIMO.
- **The dominant module-global / parallelism blocker leaves with the nutrient
  code** (unblocks T2-C threaded ensembles without a ~85-var de-global).
- **Detailed crop becomes *more* modular, not absent**: a clean interface with a
  native default and an optional detailed backend — strictly more composable than
  today's hard-wired WOFOST.
- **The crop↔water seam already exists as state** (`rd`/`ch`/`lai`/`ptra`), so
  formalizing it as an exchange record (T2-A) is low-friction and is the natural
  first structural step.
- **Preservation, not loss.** Detached code follows the ADR 0040 pattern:
  preserved on `legacy/swap-4.2.0` (and its extraction commit), re-implementable
  as a real coupled component later.

### Open questions this ADR deliberately does NOT settle

1. **Online N-limited growth.** Today the embedded WSN gives N-limited growth
   *within one SWAP run*. Offline ANIMO is one-way (SWAP→ANIMO) and does not feed
   N back to the crop. Preserving that feedback needs *online 3-way* coupling
   (SWAP water/heat ↔ ANIMO nutrients ↔ WOFOST crop). Whether to commit to that is
   a separate decision — do not assume the componentized world is a strict
   superset of today.
2. **WOFOST externalization target.** Assumes an external WOFOST with a live
   component interface (PCSE/Python WOFOST + BMI exists; maturity to be assessed).
   ANIMO's modern interface maturity is likewise unverified here.
3. **Product identity / users.** "Integrated soil-water-atmosphere-plant model"
   vs "vadose-zone engine in a coupled stack" is a community decision, not an
   engineering one. Many current users rely on the *integrated* WOFOST; a
   transition must keep an integrated mode during migration. This ADR sets
   engineering direction; the identity call belongs to the maintainers/stakeholders.

## Staged path (each stage its own ADR/arc)

1. **Boundary defined** — this ADR.
2. **Nutrients detached** (low-regret first step): remove the WSN / WOFOST-N,
   preserve on the legacy branch, document the external-ANIMO path. (Supersedes
   the "keep + de-global E-b" plan; ADR 0025–0028's reactivation is walked back.)
3. **Crop↔water interface formalized** as a typed exchange record (T2-A), with
   the embedded WOFOST refactored to sit *behind* it — the internal halfway house
   that preserves integrated mode while making the seam real.
4. **Detailed crop externalized** — swap the embedded WOFOST/grass for an external
   crop component over BMI/XMI, integrated mode retained during transition.

## Alternatives considered

- **Keep everything embedded, just de-global it (the original T1-E-b + T2-C
  plan).** Rejected as the primary direction: it commits to maintaining forks of
  WOFOST and ANIMO and to a ~85-variable nutrient migration with no oracle, for
  code that duplicates canonical standalone models. De-globaling remains the
  fallback if externalization stalls.
- **Delete crop entirely.** Rejected: the simple table-driven crop is native,
  cheap, and load-bearing for the water balance — SWAP needs a plant. Only the
  *detailed growth* models are components.
</content>
