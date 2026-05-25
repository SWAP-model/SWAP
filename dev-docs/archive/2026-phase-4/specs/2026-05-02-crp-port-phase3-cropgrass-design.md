# `.crp` port — Phase 3 (cropgrass via cases 4 + 2) — design

**Status:** Draft, awaiting review
**Date:** 2026-05-02
**Predecessors:**
- `docs/superpowers/specs/2026-05-01-swap-dra-port-design.md` (same strangler shape)
- `docs/adr/0015-strangler-narrow-scope-stub-errors.md` (stub-error pattern)
- `docs/adr/0016-per-rotation-crop-config-cache.md` (cache pattern)
- `docs/superpowers/specs/2026-05-02-crp-port-phase1-cropfixed-design.md` (Phase 1 — mirror this structure)
- **Sequencing dependency:** Phase 1 must land before Phase 3 is implemented.
  Phase 3 reuses `crop_config_global` (introduced in Phase 1) and the
  `rotation_loaded(:)` sentinel on `crop_config_t` without modification.

## Goal

Port case 4 (oxygenstress)'s `grassd.crp` (type 3, detailed grass) to TOML, with
case 2 (grassgrowth) as the secondary verification case. After this work:
- The new executable reads case 4's crop data only from `grassd.crp.toml` plus the
  typed pipeline; no `.crp` ASCII file is opened on the runtime path for type-3
  rotations.
- Case 4 and case 2 regression both remain 5/5 green throughout.
- The legacy `readgrass` reader stays alive in `readswap.f90` as a parity-test
  fixture.

This is **Phase 3 of a 4-phase `.crp` port**:
- **Phase 1** (cropfixed via case 6): see Phase 1 spec.
- **Phase 2** (cropwofost via case 5): see Phase 2 spec.
- **Phase 3** (this spec): cropgrass via case 4 (primary) + case 2 (secondary).
- **Phase 4**: case 1 (hupselbrook) integration; legacy fallback in
  `cropgrowth.f90` removed.

Each phase has its own spec, plan, and implementation cycle.

## Non-goals (Phase 3)

- Porting cropfixed or cropwofost — Phases 1 and 2.
- Porting case 1 (hupselbrook) — Phase 4.
- Removing legacy `readgrass` from `readswap.f90` — kept as parity fixture.
- Implementing `SWOXYGEN=2` Bartholomeus runtime semantics — see stub-error
  decision in "Branch scope" below.
- Implementing `SWCOMPENSATE=1` Jarvis runtime semantics — same.
- Implementing `SWDROUGHT=2` De Jong van Lier runtime — case 2 has
  `SWDROUGHT=1`; case 4 also has `SWDROUGHT=1`. Neither exercises it.
- Implementing CO2 assimilation correction (`SWCO2=1`) — both cases have
  `SWCO2=0`.
- Implementing `SWLOSSGRZ=1` / `SWLOSSMOW=1` treading-loss paths — both cases
  have these at 0.
- Implementing grazing-with-dewooling (`SEQGRAZMOW(i)=3`) — both cases use
  only value 2 (mowing); the grazing section in both files is defined but the
  `SEQGRAZMOW` sequence never activates it.
- Implementing `SWINTER=2` (Gash) or `SWINTER=3` (storage-cap) — both cases
  have `SWINTER=1`.

## Status quo

### Schema (partial)

`src/config/cropgrass_config.f90` (193 lines, Phase 4c-a + Phase 4d) covers:
- Phenology scalars: `idev`, `lcc`, `tbase`, `tsum1`, `tsum2`
- Light: `kdif`, `kdir`, `eff`, `amax`
- Tables (declared, partially parsed): `cftb`, `chtb`, `rdctb`
- Root: `rdi`, `rri`, `rdc`
- Water stress (Feddes): `hlim1`, `hlim2u`, `hlim2l`, `hlim3h`, `hlim3l`,
  `hlim4`, `adcrh`, `adcrl`, `rsc`
- Salinity: `ecmax`, `ecslop`
- Interception: `cofab`
- Mowing (Phase 4d): `swharv`, `nmow`, `dates_mowing`, `lai_after_mow`,
  `swdmmow`, `mowing_dates`, `mowing_heights`, `dmharvest`, `daylastharvest`,
  `dmlastharvest`, `maxdaymow`
- Grazing (Phase 4d): `swgraz`, `nstart_graz`, `nstop_graz`, `maxdaygrz`,
  `dmgrazing`, `swdmgrz`, `lsdb`, `tagprest`
- Per-crop irrigation: `schedule`

**Gap versus `readgrass`:** approximately 55 fields are read by `readgrass` that
have no counterpart on `cropgrass_config_t`. They fall into these groups:
- Initial crop state: `tdwi`, `laiem`, `rgrlai`, `swtsum`, `tsumtemp`,
  `tsumdepth`, `tsumtime`
- Green area: `slatb`, `ssa`, `span`
- Assimilation tables: `amaxtb`, `tmpftb`, `tmnftb`
- Biomass conversion: `cvl`, `cvr`, `cvs`
- Maintenance respiration: `q10`, `rml`, `rmr`, `rms`, `rfsetb`
- Partitioning: `frtb`, `fltb`, `fstb`
- Death rates: `perdl`, `rdrrtb`, `rdrstb`
- Crop factor / height switch: `swcf`, `albedo`, `rsw` (when `swcf=2`)
- Root extension switch and tables: `swrd`, `swdmi2rd`, `rdtb`, `rlwtb`,
  `wrtmax`, `swrdc`
- Oxygen stress: `swoxygen`, `swwrtnonox`, `aeratecrit`, `hlim2u`, `hlim2l`
  (Feddes params already partially in schema)
- Oxygen stress Bartholomeus params: `q10_microbial`, `specific_resp_humus`,
  `srl`, `swrootradius`, `dry_mat_cont_roots`, `air_filled_root_por`,
  `spec_weight_root_tissue`, `var_a`, `root_radiusO2`, `swoxygentype`,
  `swtopsub`, `nrstaring`
- Drought stress: `swdrought`
- Compensation: `swcompensate`, `swstressor`, `alphacrit`, `dcritrtz`
- Interception switch: `swinter`
- Management: `seqgrazmow`, `mowrest`, `dewrest`, `dmgrztb`, `dmmowtb`,
  `lsda`, `daysgrazing`, `uptgrazing`, `lossgrazing`, `dmmowdelay`, `daydelay`,
  `swpotrelmf`, `relmf`
- CO2: `swco2`, `co2amaxtb`, `co2efftb`, `co2tratb`, `co2year`, `co2ppm`
- Loss tables: `hlossmow`, `lossmow`, `hlossgrz`, `lossgrz`

### Parser (partial)

`src/io/toml/read_cropgrass_toml.f90` (157 lines) parses phenology, light, root,
water_stress, salinity, interception, mowing (7 scalars + 2 arrays), grazing
(7 scalars + `lsdb`), and irrigation_schedule. Does not yet handle:
- The ~55 missing schema fields above.
- The `amaxtb`, `tmpftb`, `tmnftb`, `rfsetb`, `frtb`, `fltb`, `fstb`,
  `rdrrtb`, `rdrstb`, `slatb`, `rdtb`, `rlwtb`, `dmgrztb`, `dmmowtb`,
  `dmmowdelay` tables.

### Loader (already dispatches type=3)

`src/io/toml/read_crop_toml.f90`'s rotation-loop dispatcher already routes
`rotation_type=3` entries to `read_cropgrass_toml`. In Phase 3, this dispatch
populates the `rotation_cropgrass(:)` slot on `crop_config_t`. The allocation
and slot-population logic mirrors Phase 1.

### Runtime callsite

`src/crop/cropgrowth.f90:2068` (`subroutine grass(task=1)`) calls:

```fortran
call readgrass (icrop, cropfil(icrop), swharvest, dmharvest, daylastharvest, &
                dmlastharvest, swdmmow, maxdaymow, swlossmow, swlossgrz,       &
                swdmgrz, maxdaygrz, dmgrazing, LSDb, tagprest, swhydrlift)
```

Phase 3 wraps this with an `if (rotation_loaded) ... else ... end if` guard
analogous to Phase 1's `cropfixed(task=1)` modification. The new
`cropgrass_init_from_config` subroutine replaces the legacy call on the
TOML path.

### Existing TOML case files

`tests/swap-cases/toml/4.oxygenstress/grassd.crp.toml` (53 lines) and
`tests/swap-cases/toml/2.grassgrowth/grassd.crp.toml` (59 lines) are both
skeletal: they contain only light, root, water_stress, interception, mowing,
and grazing sections. They are **separate files** because the two cases differ
in three key switches:
- `SWOXYGEN`: case 4 = 2 (Bartholomeus), case 2 = 1 (Feddes)
- `SWCOMPENSATE`: case 4 = 1 (Jarvis), case 2 = 0 (none)
- `SWHARVEST`: case 4 = 1 (DM-threshold mowing), case 2 = 2 (fixed-date mowing)

These three differences are fundamental to the crop physics and cannot be
collapsed into one shared TOML file. **Phase 3 authors and maintains two separate
`grassd.crp.toml` files.** Both reference the same underlying crop (same grass
species, same light/root/water parameters), but with different management and
stress-response settings.

## Case 4 vs case 2 branch-scope analysis

The table below summarises every switch read by `readgrass` and which value
each case exercises. The "Supported" column reflects whether Phase 3 provides
full runtime support or a stub-error.

| Switch | Case 4 value | Case 2 value | Phase 3 support |
|---|---|---|---|
| `swcf` | 2 (crop height) | 2 (crop height) | SUPPORTED |
| `swinter` | 1 (Von Hoyningen) | 1 (Von Hoyningen) | SUPPORTED |
| `swtsum` | 1 (air temp sum) | 1 (air temp sum) | SUPPORTED |
| `swrd` | 2 (max daily increase) | 2 (max daily increase) | SUPPORTED |
| `swdmi2rd` | 1 | 1 | SUPPORTED |
| `swrdc` | 0 | 0 | SUPPORTED |
| `swoxygen` | 2 (Bartholomeus) | 1 (Feddes) | SUPPORTED for value=1; STUB-ERROR for value=2 |
| `swwrtnonox` | 1 | 1 | SUPPORTED |
| `swdrought` | 1 (Feddes) | 1 (Feddes) | SUPPORTED |
| `swsalinity` | 0 (none) | 0 (none) | SUPPORTED (passthrough) |
| `swcompensate` | 1 (Jarvis) | 0 (none) | STUB-ERROR for value≠0 |
| `swdmmow` | 2 (flexible table) | 2 (flexible table) | SUPPORTED |
| `swharvest` | 1 (DM-threshold) | 2 (fixed dates) | SUPPORTED (both) |
| `swdmgrz` | 2 (flexible table) | 2 (flexible table) | SUPPORTED |
| `swlossgrz` | 0 | 0 | SUPPORTED (passthrough; stub-error for value=1) |
| `swlossmow` | 0 | 0 | SUPPORTED (passthrough; stub-error for value=1) |
| `SEQGRAZMOW` | all 2 (mow-only) | all 2 (mow-only) | SUPPORTED for value=2; STUB-ERROR for values=1 and =3 |
| `swco2` | 0 | 0 | SUPPORTED (passthrough; stub-error for value=1) |
| `swpotrelmf` | 1 | 1 | SUPPORTED |
| `schedule` (irrigation) | 0 | 0 | SUPPORTED (passthrough; stub-error for value=1) |

**SWOXYGEN=2 (Bartholomeus) decision:** case 4 is the PRIMARY validation case
for Phase 3 and it exercises `SWOXYGEN=2`. However, the Bartholomeus oxygen
stress module (`src/oxygenstress_mod`) reads its inputs through module globals
(`q10_microbial`, `srl`, `swrootradius`, `root_radiusO2`, etc.) that are
currently populated by `readgrass`. Providing full runtime support for
`SWOXYGEN=2` requires adding these ~12 fields to `cropgrass_config_t`, extending
the parser, writing them into module globals in `cropgrass_init_from_config`,
and verifying the oxygen stress output column (`treddry`/`tredwet`) in the
regression reference. This work is bounded and well-understood.

**Decision: SUPPORT `SWOXYGEN=2` in Phase 3, schema-1:1.** The Bartholomeus
parameters are purely inputs from the `.crp` file (no feedback loops at
init time); the runtime module itself is unchanged. Carrying the stub-error
would make case 4 unable to run on the TOML path at all, undermining the
primary validation goal of Phase 3.

**SWCOMPENSATE=1 (Jarvis) decision:** case 4 exercises `SWCOMPENSATE=1`.
Support requires threading `alphacrit` and `swstressor` into the compensation
module. This is straightforward schema-1:1 work. **Decision: SUPPORT
`SWCOMPENSATE=1` in Phase 3.** Case 2 has `SWCOMPENSATE=0` (passthrough, no
extra work). `SWCOMPENSATE=2` (Walsum) is stub-errored: neither case exercises
it and the Walsum-specific parameter `dcritrtz` introduces additional plumbing.

**SEQGRAZMOW with grazing (value=1) or dewooling (value=3):** both cases use
only mowing (value=2). The grazing block in `readgrass` is substantial (~100
lines) and requires `lsda`, `daysgrazingtab`, `uptgrazingtab`,
`lossgrazingtab` arrays that do not currently have counterparts on
`cropgrass_config_t`. **Decision: STUB-ERROR seqgrazmow values 1 and 3.**
The existing schema fields for grazing scalars (`dmgrazing`, `tagprest`, `lsdb`,
`swdmgrz`, `maxdaygrz`) stay in the schema (they are already there) but are
not written to module globals on the init path when `seqgrazmow` is all-2
(mowing-only). A validator guard rejects any configuration that includes
`seqgrazmow(i) ∈ {1, 3}` with a stub-error.

## Approach

Mirror Phase 1's structure: schema extension + parser extension + loader
dispatch + `populated` sentinel + new runtime init module + wiring change.

### Unit 1 — Schema extension 1:1 with `readgrass`

Extend `cropgrass_config_t` in `src/config/cropgrass_config.f90` so every field
that `readgrass` reads has a TOML home, regardless of whether the parent switch
is at its case-4 or case-2 value. Groups:

**New scalar fields:**
- Crop state init: `tdwi`, `laiem`, `rgrlai`
- Start-of-growth: `swtsum`, `tsumtemp`, `tsumdepth`, `tsumtime`
- Green area: `ssa`, `span` (note: `tbase` already on type)
- Assimilation: (already: `kdif`, `kdir`, `eff`, `amax` — note these are
  scalars on the existing type; `readgrass` reads these but also reads tables
  `amaxtb`, `tmpftb`, `tmnftb` — the tables are what the runtime uses)
- Biomass conversion: `cvl`, `cvr`, `cvs`
- Maintenance respiration: `q10`, `rml`, `rmr`, `rms`
- Death rates: `perdl`
- Crop factor/height switch: `swcf`, `albedo`, `rsw`
- Root depth/density switches: `swrd`, `swdmi2rd`, `swrdc`, `wrtmax`
- Oxygen stress: `swoxygen`, `swwrtnonox`, `aeratecrit`
- Oxygen stress Bartholomeus scalars: `swoxygentype`, `swrootradius`,
  `q10_microbial`, `specific_resp_humus`, `srl`, `dry_mat_cont_roots`,
  `air_filled_root_por`, `spec_weight_root_tissue`, `var_a`, `root_radiusO2`
- Drought stress switch: `swdrought`
- Compensation: `swcompensate`, `swstressor`, `alphacrit`, `dcritrtz`
- Interception switch: `swinter`
- Management: `mowrest`, `dewrest`, `swpotrelmf`, `relmf`, `nseqgrazmow`
- Grazing (already partially present; add): `dewrest`
- CO2: `swco2`

**New table fields** (`real(real64), allocatable :: name(:)` flat DVS/DNR-value
paired arrays):
- `slatb` — specific leaf area vs day-number (size ≤ 30 pairs)
- `amaxtb` — max CO2 assimilation rate vs day-number (size ≤ 15 pairs)
- `tmpftb` — AMAX reduction factor vs average day temperature (size ≤ 15)
- `tmnftb` — AMAX reduction factor vs minimum day temperature (size ≤ 15)
- `rfsetb` — senescence reduction factor vs day-number (size ≤ 15)
- `frtb` — fraction partitioned to roots vs day-number (size ≤ 15)
- `fltb` — fraction partitioned to leaves vs day-number (size ≤ 15)
- `fstb` — fraction partitioned to stems vs day-number (size ≤ 15)
- `rdrrtb` — relative death rate of roots vs day-number (size ≤ 15)
- `rdrstb` — relative death rate of stems vs day-number (size ≤ 15)
- `rdtb` — rooting depth vs day-number (when `swrd=1`; size ≤ 2*MAGRS)
- `rlwtb` — rooting depth vs root weight (when `swrd=3`; size ≤ 22)
- `dmmowtb` — flexible DM threshold for mowing vs day-number (when `swdmmow=2`)
- `dmgrztb` — flexible DM threshold for grazing vs day-number (when `swdmgrz=2`)
- `dmmowdelay_dm(:)` / `dmmowdelay_day(:)` — DM harvest and regrowth delay
  pairs (or a single flat array `dmmowdelay` mirroring legacy storage)

The `populated` sentinel is also added:

```fortran
logical :: populated = .false.   ! set .true. at end of read_cropgrass_toml
```

**Stub-error validators added to `cropgrass_config_validate`:**
- `swoxygen == 2 .and. swoxygentype == 2` — reproduction-function sub-branch
  of Bartholomeus (`swtopsub` / `nrstaring`); neither case exercises this.
  Stub-error only for that sub-branch; `swoxygentype=1` (physical) is supported.
- `swcompensate == 2` (Walsum) — neither case exercises it.
- `swinter == 2 .or. swinter == 3` (Gash / storage-cap) — neither case.
- `swdrought == 2` (De Jong van Lier) — neither case.
- `swsalinity /= 0` — neither case.
- `swco2 == 1` (CO2 correction) — neither case.
- `swlossgrz == 1` (treading losses during grazing) — neither case.
- `swlossmow == 1` (treading losses during mowing) — neither case.
- `any(seqgrazmow(1:nseqgrazmow) /= 2)` — grazing or dewooling periods — neither case.
- `swrd == 3` — root-biomass-based extension — neither case (both use `swrd=2`).
- `schedule%schedule == 1` (irrigation scheduling) — neither case.

### Unit 2 — Parser extension

`src/io/toml/read_cropgrass_toml.f90` extended to read:
- All new scalars under their natural TOML sections (phenology, light,
  water_stress, root, mowing, grazing, interception).
- New sections:
  - `[green_area]` — `ssa`, `span`; table `slatb` as flat array
  - `[assimilation]` — tables `amaxtb`, `tmpftb`, `tmnftb`
  - `[biomass_conversion]` — `cvl`, `cvr`, `cvs`
  - `[respiration]` — `q10`, `rml`, `rmr`, `rms`; table `rfsetb`
  - `[partitioning]` — tables `frtb`, `fltb`, `fstb`
  - `[death_rates]` — `perdl`; tables `rdrrtb`, `rdrstb`
  - `[crop_factor]` — `swcf`, `albedo`, `rsw`; tables `cftb` (already
    declared), `chtb` (already declared)
  - `[oxygen_stress]` — `swoxygen`, `swwrtnonox`, `aeratecrit`; Bartholomeus
    scalars under `[oxygen_stress.bartholomeus]`
  - `[drought_stress]` — `swdrought` (and subordinate hlim/adcr fields already
    under `[water_stress]`)
  - `[compensation]` — `swcompensate`, `swstressor`, `alphacrit`, `dcritrtz`
  - `[management]` — `swtsum`, `tsumtemp`, `tsumdepth`, `tsumtime`, `mowrest`,
    `dewrest`, `swpotrelmf`, `relmf`, `seqgrazmow` (integer array)
  - `[management.mowing_delay]` — flat pairs for `dmmowdelay` table
  - `[co2]` — `swco2` (scalar; remaining CO2 fields stub-errored by validator)

Table-reading uses the existing `read_array_1d` private helper (already in the
module).

`populated = .true.` is set at the end of `read_cropgrass_toml`.

A new file-based wrapper is added:

```fortran
subroutine read_cropgrass_file_toml(path, config, errors, base_path)
   ! Resolves <base_path>/<path>, opens via toml_load, calls
   ! read_cropgrass_toml(doc, config, errors).
end subroutine
```

### Unit 3 — Loader dispatch for type=3

`src/io/toml/read_crop_toml.f90` already dispatches `rotation_type=3` to
`read_cropgrass_toml`. In Phase 3 this dispatch is extended to:
1. Allocate `rotation_cropgrass(nrot)` on `crop_config_t` (using `rotation_grass`
   under the actual naming — see Phase 1 naming-reconciliation table).
2. For each `rotation_type=3` entry, call
   `read_cropgrass_file_toml(rotation_file(i), rotation_cropgrass(i), errors, base_path)`.
3. Leave slots for type=1 and type=2 entries at their default-initialized states
   (those are populated by Phases 1 and 2 respectively).

### Unit 4 — Runtime init module: `cropgrass_init`

New file: `src/crop/cropgrass_init.f90`. Public sub:

```fortran
subroutine cropgrass_init_from_config(cfg, icrop)
   class(cropgrass_config_t), intent(in) :: cfg
   integer,                   intent(in) :: icrop
end subroutine
```

Two halves:

**Half 1 — Config → globals copy.** One-for-one mirror of `readgrass`'s `rd*`
calls. Same module globals targeted (`tdwi`, `laiem`, `slatb`, `amaxtb`,
`kdif`, `kdir`, `eff`, `frtb`, `fltb`, `fstb`, `hlim1`, `hlim2u`, `hlim2l`,
`hlim3h`, `hlim3l`, `hlim4`, `adcrh`, `adcrl`, `swoxygen`, `swdrought`,
`swcompensate`, `alphacrit`, `mowrest`, `seqgrazmow`, `swharvest`, `swdmmow`,
`dmharvest`, `daylastharvest`, `dmlastharvest`, `maxdaymow`, `dmmowtb`,
`dmmowdelay` / `delayregrowthtab`, `cofab`, etc.).

When `swoxygen == 2` the Bartholomeus-path globals are also written
(`q10_microbial`, `srl`, `swrootradius`, `root_radiusO2`, etc.).

When `swharvest == 2` the `dateharvest` array is populated from
`cfg%mowing_dates` (already a DNR flat array in the schema).

**Half 2 — Runtime init math.** Verbatim port of `readgrass`'s tail:
- Build `cumdens(:)` from `rdctb` (identical to cropfixed; same 100-interval
  trapezium-sum when `swdrought=1`).
- Initialize `daycrop`, `nofd`, `rid`, `daygrowth`, `daygrowthpot` defaults.

**Defense-in-depth runtime guards** (matching the validator stub-errors):

```fortran
if (swdrought == 2 .or. swsalinity /= 0 .or. &
    swcompensate == 2 .or. swco2 == 1 .or. &
    swlossgrz == 1 .or. swlossmow == 1) then
   call fatalerr_collected('cropgrass_init', &
      'Unsupported runtime branch reached on the TOML path. ' // &
      'The validator should have caught this earlier.')
end if
if (swoxygentype == 2) then
   call fatalerr_collected('cropgrass_init', &
      'swoxygen=2 swoxygentype=2 (reproduction functions) not yet ' // &
      'supported on the TOML path. The validator should have caught this.')
end if
```

Note: the `readgrass` argument list passes `swharvest`, `dmharvest`,
`daylastharvest`, `dmlastharvest`, `swdmmow`, `maxdaymow`, `swlossmow`,
`swlossgrz`, `swdmgrz`, `maxdaygrz`, `dmgrazing`, `lsdb`, `tagprest`,
`swhydrlift` as OUT arguments. `cropgrass_init_from_config` writes the same
values to the same module globals but never needs to return them via
arguments — the globals themselves are the return channel.

### Unit 5 — Wiring change in `cropgrowth.f90:2068`

Before:
```fortran
call readgrass (icrop, cropfil(icrop), swharvest, dmharvest, daylastharvest, &
                dmlastharvest, swdmmow, maxdaymow, swlossmow, swlossgrz,       &
                swdmgrz, maxdaygrz, dmgrazing, LSDb, tagprest, swhydrlift)
```

After:
```fortran
if (associated(crop_config_global) .and. &
    allocated(crop_config_global%rotation_grass) .and. &
    crop_config_global%rotation_grass(icrop)%populated) then
   call cropgrass_init_from_config(crop_config_global%rotation_grass(icrop), icrop)
else
   call readgrass (icrop, cropfil(icrop), swharvest, dmharvest,              &
                   daylastharvest, dmlastharvest, swdmmow, maxdaymow,         &
                   swlossmow, swlossgrz, swdmgrz, maxdaygrz, dmgrazing,       &
                   LSDb, tagprest, swhydrlift)   ! transitional
end if
```

The `else` branch is **transitional** (Phase 4 removes it).

Note on naming: the actual name of the array on `crop_config_t` must be
confirmed against the current `src/config/crop_config.f90` before
implementation (it may be `rotation_grass` or `rotation_cropgrass`; Phase 1's
naming-reconciliation table shows the actual names diverged from the spec's
proposed names).

### Unit 6 — TOML content authoring

Two submodule files require full authoring from their respective legacy `.crp`
ASCII sources:

**`tests/swap-cases/toml/4.oxygenstress/grassd.crp.toml`** — extended from 53
lines to a full 1:1 reproduction of `tests/swap-cases/4.oxygenstress/grassd.crp`
(524 lines). Key entries:
- `swcf = 2`, `chtb` table (3 rows: DNR 0/180/366, CH=12.0)
- `albedo = 0.23`, `rsw = 0.0`
- `tdwi = 1000.0`, `laiem = 0.63`, `rgrlai = 0.007`
- `swtsum = 1`
- `slatb` table (4 rows: 1.0→0.0015, 80.0→0.0015, 300.0→0.0020, 366.0→0.0020)
- `amaxtb` (5 rows), `tmpftb` (5 rows), `tmnftb` (2 rows)
- `cvl=0.685`, `cvr=0.694`, `cvs=0.662`
- `q10=2.0`, `rml=0.03`, `rmr=0.015`, `rms=0.015`, `rfsetb` (2 rows)
- `frtb` (2 rows), `fltb` (2 rows), `fstb` (2 rows)
- `perdl=0.05`, `rdrrtb` (3 rows), `rdrstb` (3 rows)
- `swrd = 2`, `rdi=10.0`, `rri=1.0`, `rdc=40.0`, `swdmi2rd=1`, `swrdc=0`
- `rdtb` (3 rows), `rlwtb` (2 rows), `wrtmax=3000.0`
- `rdctb` (2 rows: 0.0→1.0, 1.0→0.0)
- `swoxygen = 2` (Bartholomeus), `swwrtnonox=1`, `aeratecrit=0.5`
- Bartholomeus params: `q10_microbial=2.8`, `specific_resp_humus=1.6e-3`,
  `srl=383571.0`, `swrootradius=2`, `root_radiusO2=0.000075`
- `swdrought = 1`, `hlim3h=-200.0`, ..., `adcrl=0.1`
- `swcompensate = 1`, `swstressor=1`, `alphacrit=0.7`
- `swinter = 1`, `cofab = 0.25`
- Management: `seqgrazmow = [2,2,...,2]` (20 values), `swharvest=1`,
  `swdmmow=2`, `dmmowtb` table (5 rows), `maxdaymow=42`, `mowrest=700.0`,
  `dmmowdelay` pairs (3 rows), `swdmgrz=2`, `dmgrztb` (3 rows),
  `maxdaygrz=28`, `tagprest=700.0`, `lsdb=[21.25]`,
  `daysgrazing=[4.0]`, `uptgrazing=[16.0]`, `lossgrazing=[4.00]`,
  `swpotrelmf=1`, `relmf=0.90`
- `swco2 = 0`
- `schedule = 0`

**`tests/swap-cases/toml/2.grassgrowth/grassd.crp.toml`** — extended from 59
lines to a full 1:1 reproduction of `tests/swap-cases/2.grassgrowth/grassd.crp`
(541 lines). Same structure as case 4 except:
- `swoxygen = 1` (Feddes), `hlim1=0.0`, `hlim2u=1.0`, `hlim2l=-1.0` (no
  Bartholomeus block)
- `swcompensate = 0` (no Jarvis block)
- `swharvest = 2` (fixed dates), `nmow=30`, `mowing_dates = [127.0, ...]`
  (30 DOY values derived from the legacy `dateharvest` calendar dates;
  already present in the existing TOML stub)

Both submodule updates each require an inner-repo commit + outer-repo pointer
bump (separate commits for the author step and the delete step — see
teardown).

### Unit 7 — Keep `readgrass` alive for parity tests

`readgrass` and its caller in `cropgrowth.f90` (the `else` branch) stay in
`readswap.f90` untouched. Used as the parity reference by
`tests/unit/io/toml/test_grasscrop_parity.pf`, which compares
`cropgrass_init_from_config(rotation_grass(i))` against
`readgrass`-populated globals starting from the same `grassd.crp` fixture.
Cleanup (deleting `readgrass` and the `else` fallback) is Phase 4's
responsibility.

### Unit 8 — Delete `grassd.crp` from TOML case directories

`tests/swap-cases/toml/4.oxygenstress/grassd.crp` is removed (submodule
inner commit + outer-repo bump). `tests/swap-cases/toml/2.grassgrowth/grassd.crp`
is also removed (separate submodule pair commit). The legacy ASCII copies
in `tests/swap-cases/4.oxygenstress/grassd.crp` and
`tests/swap-cases/2.grassgrowth/grassd.crp` stay for the legacy executable.

## Teardown plan (strangler-fig hygiene)

| Element | Introduced for | Teardown trigger | What gets removed | Replacement |
|---|---|---|---|---|
| `else call readgrass(...)` fallback in `cropgrowth.f90:2068` | Phase 3 ships before Phase 4; type-1/2 rotations still need legacy readers | Phase 4 (hupselbrook) lands all three types | The `if/else` dispatch around all three legacy reader calls | Unconditional `*_init_from_config` calls |
| Legacy `readgrass` in `readswap.f90` (~800 lines) | Parity-test fixture | Future cleanup: parity tests rewritten to compare against fixture values | The subroutine | Fixture-based parity tests |
| `crop_config_global` module-level pointer | Introduced in Phase 1; reused here | ADR 0016 config-passing direction lands (post-Phase 4) | `crop_config_global_mod` module | Explicit config+state argument threading |
| `populated :: logical` sentinel on `cropgrass_config_t` | Runtime dispatch guard | Same as `crop_config_global` | The `populated` field | Call-site dispatch on `rotation_type(icrop)` |

## Test plan

### New unit tests

`tests/unit/config/test_cropgrass_config.pf` (extend):
- `test_cropgrass_swdrought2_rejected` — ERR_VALIDATION_CROSS_FIELD
- `test_cropgrass_swcompensate2_rejected`
- `test_cropgrass_swinter2_rejected`
- `test_cropgrass_swsalinity1_rejected`
- `test_cropgrass_swco2_1_rejected`
- `test_cropgrass_swlossgrz1_rejected`
- `test_cropgrass_swlossmow1_rejected`
- `test_cropgrass_seqgrazmow_with_grazing_rejected` — seqgrazmow contains 1
- `test_cropgrass_swoxygentype2_rejected`
- `test_cropgrass_swcompensate1_with_swrd2_passes` — case-4 supported values
  produce no stub-errors
- `test_cropgrass_swharvest2_with_mowing_dates_passes` — case-2 supported values

`tests/unit/io/toml/test_read_cropgrass_toml.pf` (extend):
- `test_cropgrass_round_trip_case4_values` — round-trip all new fields from a
  fixture covering case-4 values (swoxygen=2, swcompensate=1, swharvest=1)
- `test_cropgrass_round_trip_case2_values` — round-trip case-2 values
  (swoxygen=1, swcompensate=0, swharvest=2 with mowing_dates)
- `test_cropgrass_populated_set_true` — `populated` flag is true after parsing

`tests/unit/io/toml/test_load_swap_config.pf` (extend):
- `test_swap_config_case4_rotation_grass_loaded` — loading case 4's swap.toml
  populates all 10 `rotation_grass(:)` slots with `populated=.true.`
- `test_swap_config_case2_rotation_grass_loaded` — loading case 2's swap.toml
  populates all 5 `rotation_grass(:)` slots with `populated=.true.`

`tests/unit/crop/test_cropgrass_init.pf` (new):
- `test_cropgrass_init_writes_kdif_kdir` — verifies light globals are set
- `test_cropgrass_init_cumdens_case4` — hand-computed `cumdens` for case-4's
  `rdctb = [(0.0, 1.0), (1.0, 0.0)]`; same trapezium-sum as cropfixed
- `test_cropgrass_init_swoxygen2_writes_srl` — verifies `srl` global is
  written when `swoxygen=2`
- `test_cropgrass_init_swcompensate1_writes_alphacrit` — verifies `alphacrit`
  global set when `swcompensate=1`
- `test_cropgrass_init_swharvest2_populates_dateharvest` — verifies
  `dateharvest` array is populated from `mowing_dates` when `swharvest=2`

### Parity tests

Add or extend `tests/unit/io/toml/test_grasscrop_parity.pf`:
- Assert that `cropgrass_init_from_config(rotation_grass(i))` produces the
  same module-global state as `readgrass` on the same `grassd.crp` fixture for
  case-4 configuration (swoxygen=2, swharvest=1).
- Repeat for case-2 configuration (swoxygen=1, swharvest=2).

### Regression (both cases)

Both case 4 (PRIMARY) and case 2 (SECONDARY) must remain 5/5 green after each
task. The smoke step (Task 11 in the plan) renames both `grassd.crp` files to
`.disabled` and re-runs regression before the delete commit.

### Acceptance gate

- `pixi run test-pfunit` → all green.
- `pixi run regression` → 5/5 cases green.
- `tests/swap-cases/toml/4.oxygenstress/grassd.crp` → does not exist.
- `tests/swap-cases/toml/2.grassgrowth/grassd.crp` → does not exist.
- `git grep "call readgrass" src/` → returns the legacy `else` fallback in
  `cropgrowth.f90:2068` only. Other reachability gone.
- `git grep "call readgrass" tests/` → returns only parity-test references.
- `crop_config_global` exists (introduced in Phase 1), `rotation_grass(:)`
  populated by `read_crop_toml` for type=3 entries.

## Implementation tasks (working draft)

Elaborated in `writing-plans`. 12 tasks total:

1. Audit `readgrass` (readswap.f90:3437-4270) line-by-line into READ/VALIDATE/
   NORMALIZE/RUNTIME/GUARDED buckets. Docs only.
2. Extend `cropgrass_config_t` schema 1:1 with grassd.crp. Add stub-error
   validators. Add `populated` sentinel.
3. Extend `read_cropgrass_toml.f90` parser for all new fields + tables. Set
   `populated = .true.` at end.
4. Add `read_cropgrass_file_toml` wrapper for file-based loading.
5. Author `tests/swap-cases/toml/4.oxygenstress/grassd.crp.toml` 1:1 with
   legacy case-4 values. Submodule pair commit.
6. Author `tests/swap-cases/toml/2.grassgrowth/grassd.crp.toml` 1:1 with
   legacy case-2 values. Submodule pair commit.
7. Create `src/crop/cropgrass_init.f90` with `cropgrass_init_from_config`
   (config-to-globals + runtime init math + defense-in-depth guards).
8. Extend `read_crop_toml.f90` loader dispatch to allocate and populate
   `rotation_cropgrass(:)` for type=3 entries.
9. Wire `cropgrowth.f90:2068` to dispatch on the `populated` sentinel; legacy
   `readgrass` in the `else` branch.
10. Write/extend unit tests (stub-error, round-trip, init globals, cumdens).
11. Smoke test: rename `grassd.crp` → `grassd.crp.disabled` in BOTH TOML case
    dirs, run regression, confirm 5/5, restore.
12. Delete `grassd.crp` from BOTH TOML case dirs (two separate submodule pair
    commits). Update `docs/csv-companion-files.md` path-resolution note.

## Risk register

- **`cumdens` math in `readgrass` tail:** the same trapezium-sum computation
  as cropfixed (lines 4091-4121). Verbatim port protects parity. Unit test
  with hand-computed expected value for case-4's `rdctb`.
- **`dateharvest` encoding in TOML for case-2 (`swharvest=2`):** the legacy
  file stores absolute calendar dates (e.g. `1980-05-06`). The existing TOML
  stub already converts them to day-of-year floats. Verify the conversion is
  correct before authoring the full TOML.
- **`seqgrazmow` integer array:** `readgrass` reads this as a variable-length
  sequence (up to 366 values). The TOML parser uses `read_array_1d` for
  `real64` — `seqgrazmow` is integer. Need an `read_array_1d_int` variant or
  inline the loop.
- **Phase 1 sequencing dependency:** `crop_config_global` and the
  `rotation_loaded(:)` sentinel are introduced in Phase 1. Phase 3 must not
  be implemented until Phase 1 is in the repository.
- **Bartholomeus init ordering:** `readgrass` validates `swoxygen=2` against
  `swhea` / `swcalt` / `bdens` (module globals set by the swap.ini reader).
  These cross-section checks must be replicated in `cropgrass_config_validate`
  or deferred to a `swap_config_validate` cross-section check, where both
  `heat.swcalt` and `cropgrass.swoxygen` are in scope. The validator approach
  is cleaner; note this as a potential scope creep point.
- **Grass `readgrass` argument list has 16 parameters:** the wiring change in
  `cropgrowth.f90` must not alter the signature of the legacy fallback call —
  only add the `if/else` guard around it.
