---
title: Phase 4c .crp file coverage audit
author: SWAP modernization team
date: 2026-04-25
---

# Phase 4c .crp file coverage audit

## Summary

SWAP legacy input uses three crop types in `.crp` (crop) files:
1. **Type 1**: Fixed crop (simple model) — e.g. `maizes.crp`
2. **Type 2**: WOFOST general (detailed physiological model) — e.g. `potatod.crp`
3. **Type 3**: WOFOST grass (detailed grass-specific model) — e.g. `grassd.crp`

The new TOML pipeline currently covers **rotation metadata only**: `swcrop` (enable/disable), plus parallel arrays `rotation_start(:)`, `rotation_end(:)`, `rotation_file(:)`, `rotation_type(:)`. The **`.crp` file contents are completely unparsed** by the new path; no data from the referenced `.crp` files is loaded.

Phase 4c must add:
- Three new config types: `cropfixed_config_t`, `cropgrass_config_t`, `cropwofost_config_t`
- Three new TOML readers (or unified reader with conditional branches)
- Cross-file TOML loading: `crop.rotation[N].file = "maizes.crp.toml"` → load that TOML file

## Cases × .crp files surveyed

| Case | .crp files present | Types used |
|---|---|---|
| 1.hupselbrook | maizes.crp, potatod.crp, grassd.crp | 1, 2, 3 |
| 2.grassgrowth | grassd.crp | 3 |
| 4.oxygenstress | grassd.crp | 3 |
| 5.salinitystress | potatod.crp | 2 |
| 6.surfacewater | grass.crp | 1 |
| 3.macroporeflow (excluded) | wintcer1.crp, wintcer2.crp | 1, 1 |

**Total: 8 unique `.crp` files across 6 cases; 5 non-macropore cases use types 1, 2, 3 all present.**

## What's already in the TOML pipeline

**File: `src/config/crop_config.f90` (75 lines)**
- `crop_config_t` holds:
  - `swcrop` (0=disable, 1=enable crops)
  - `rotation_start(:)`, `rotation_end(:)` (days since 1900)
  - `rotation_file(:)` (`.crp` filename without extension, e.g. "maizes")
  - `rotation_type(:)` (1, 2, or 3)
- Validation: checks array lengths match, `rotation_type` ∈ {1,2,3}, start < end
- **No support for per-crop data**

**File: `src/io/toml/read_crop_toml.f90` (77 lines)**
- Reads `[crop]` section and `[[crop.rotation]]` array
- Parses start/end dates, file name, type
- **Does not read per-crop `.crp.toml` files**

## What's missing per crop type

### Fixed crop (type 1) — ~80 parameters

Used in cases: 1.hupselbrook (maizes), 5.salinitystress (potatod), 6.surfacewater (grass).

**Development & phenology (6 fields)**
- `IDEV` (1=fixed, 2=temperature-sum-based)
- `LCC` (if IDEV=1; duration in days)
- `TSUMEA`, `TSUMAM`, `TBASE` (if IDEV=2; temperature sums and base)
- `SWHARV`, `DVSEND` (harvest timing: development stage threshold)

**Light interception (2 fields)**
- `KDIF`, `KDIR` (extinction coefficients for diffuse/direct light)

**Leaf area or soil cover (1 field + 1 table)**
- `SWGC` (1=LAI, 2=soil cover fraction)
- `GCTB` (table: DVS → LAI or SCF, up to 8 entries)

**Crop factor or height (1 field + 1–2 tables)**
- `SWCF` (1=crop factor, 2=crop height, 3=wet crop factor option)
- `CFTB` (table: DVS → crop factor, ~10 entries)
- `CHTB` (table: DVS → crop height cm, ~10 entries; if SWCF=2)

**Root development (5 fields + up to 3 tables)**
- `SWRD` (1=DVS-based, 2=daily increase, 3=biomass-based)
- `RDTB` (if SWRD=1; DVS → rooting depth)
- `RDI`, `RRI`, `RDC` (if SWRD=2; initial depth, daily rate, max)
- `RLWTB` (if SWRD=3; root weight → depth)
- `SWDMI2RD`, `WRTMAX`, `RDCTB` (root density profile)

**Oxygen stress (8–12 fields, depends on SWOXYGEN)**
- `SWOXYGEN` (0=none, 1=Feddes, 2=Bartholomeus)
- `SWWRTNONOX` (check aerobic conditions)
- `AERATECRIT` (threshold for root extension)
- If SWOXYGEN=1: `HLIM1`, `HLIM2U`, `HLIM2L` (3 pressure head limits)
- If SWOXYGEN=2: `Q10_MICROBIAL`, `SPECIFIC_RESP_HUMUS`, `SRL`, `SWROOTRADIUS`, `ROOT_RADIUSO2` (7 respiration params)

**Drought stress (5 fields, Feddes only)**
- `SWDROUGHT` (1=Feddes, 2=De Jong van Lier)
- If SWDROUGHT=1: `HLIM3H`, `HLIM3L`, `HLIM4`, `ADCRH`, `ADCRL` (5 pressure head / demand limits)

**Salinity stress (3 fields, optional)**
- `SWSALINITY` (0=none, 1=Maas & Hoffman, 2=osmotic head)
- If SWSALINITY≥1: `SALTMAX`, `SALTSLOPE` (threshold & slope)
- If SWSALINITY=2: `SALTHEAD` (conversion factor)

**Stress compensation (3 fields, optional)**
- `SWCOMPENSATE` (0=none, 1=Jarvis, 2=Walsum)
- `SWSTRESSOR` (which stress to compensate)
- `ALPHACRIT` or `DCRITRTZ` (depending on method)

**Interception (1 field + up to 5 tables, depending on SWINTER)**
- `SWINTER` (0=none, 1=Von Hoyningen-Hune & Braden, 2=Gash)
- If SWINTER=1: `COFAB` (interception coefficient)
- If SWINTER=2: up to 5 tables (PFREE, PSTEM, SCANOPY, AVPREC, AVEVAP vs. time)

**Pre-sowing, sowing, germination (10 fields, optional)**
- `SWPREP` (0=no, 1=yes)
- `SWSOW` (0=no, 1=yes)
- `SWGERM` (0=no, 1=temperature, 2=temperature+hydrology)
- If enabled: 3–6 fields per option (depths, temps, delays)

**Irrigation scheduling (20+ fields, optional)**
- `SCHEDULE` (0=no, 1=yes)
- If yes: `TCS` (timing criterion: 1–8), `DCS` (depth criterion: 1–2), plus ~20 subsidiary fields

**Estimated LoC for cropfixed_config_t:**
- Type definition + validators: ~200 lines (F90)
- TOML reader: ~300 lines
- **Total: ~500 lines**

### Grass (type 3) — ~110 parameters

Used in cases: 1.hupselbrook, 2.grassgrowth, 4.oxygenstress.

**Differences from fixed crop:**
- Use **day-of-year (DNR)** instead of development stage (DVS)
- Perennial crop: no explicit sowing/harvest, but mowing/grazing cycles
- Separate management section (Part 13–15): mowing, grazing, livestock density
- Grass-specific init: `SWTSUM` (temperature-based start)
- No storage organs: no CVO, RMO, FOTB

**All fixed-crop sections + grass-specific:**
- Grass init (4 fields): `TDWI`, `LAIEM`, `RGRLAI`, `SWTSUM` + `TSUMTEMP`, `TSUMDEPTH`, `TSUMTIME`
- Green surface area (4 fields + 1 table): `SSA`, `SPAN`, `TBASE`, `SLATB`
- Assimilation (3 fields + 3 tables): `KDIF`, `KDIR`, `EFF`, `AMAXTB`, `TMPFTB`, `TMNFTB`
- Respiration (4 fields + 1 table): `Q10`, `RML`, `RMS`, `RMR`, `RFSETB`
- Partitioning (3 tables): `FRTB`, `FLTB`, `FSTB` (no FOTB)
- Death rates (3 tables): `PERDL`, `RDRRTB`, `RDRSTB`
- Management (10+ fields + 4 tables): `SEQGRAZMOW`, `SWHARVEST`, `SWDMGRZ`, `DMGRAZING`, `DMGRZTB`, `MAXDAYGRZ`, `SWLOSSGRZ`, `TAGPREST`, `DEWREST`, `LSDA` table, `DMMOW`/`DMGRZ` tables
- CO2 impact (1 field): `SWCO2` + 3 tables if enabled
- Stress compensation, drought, salinity, oxygen: **same as fixed crop**

**Estimated LoC for cropgrass_config_t:**
- Type definition: ~300 lines (many fields + large tables)
- TOML reader: ~400 lines
- **Total: ~700 lines**

### WOFOST (type 2) — ~150+ parameters, mostly tables

Used in cases: 1.hupselbrook (potatod), 5.salinitystress (potatod).

WOFOST is a full physiological model with detailed tables. **Do NOT list every parameter**; instead group by category:

**Phenological development (7 fields + 2 tables)**
- `IDSL` (0=temperature, 1=temperature+daylength, 2=+vernalisation)
- `TSUMEA`, `TSUMAM` (temperature sums to anthesis, to maturity)
- `DTSMTB` (table: air temp → Tsum increment)
- `DLO`, `DLC` (daylength threshold)
- `VERNSAT`, `VERNBASE`, `VERNDVS`, `VERNTB` (if IDSL=2; vernalisation)

**Assimilation (3 fields + 3 tables)**
- `KDIF`, `KDIR`, `EFF` (light extinction & efficiency)
- `AMAXTB` (DVS → max CO2 assimilation)
- `TMPFTB`, `TMNFTB` (temp reduction factors)

**Conversion & respiration (5 fields + 1 table)**
- `CVL`, `CVO`, `CVR`, `CVS` (conversion efficiency for leaves, organs, roots, stems)
- `Q10`, `RML`, `RMO`, `RMR`, `RMS` (Q10 + maintenance rates for each organ)
- `RFSETB` (senescence reduction factor)

**Dry-matter partitioning (4 tables)**
- `FRTB`, `FLTB`, `FSTB`, `FOTB` (fractions to roots, leaves, stems, storage organs as f(DVS))
- ~20 entries per table

**Death rates & senescence (3 tables)**
- `PERDL` (leaf death rate)
- `RDRRTB`, `RDRSTB` (root, stem death rates vs. DVS)
- ~6 entries per table

**Root development (5 fields + 3 tables)**
- `SWRD`, `RDI`, `RRI`, `RDC`, `SWDMI2RD` (same as fixed)
- `RDTB`, `RLWTB`, `RDCTB` (same tables as fixed)

**Stress response (30 fields)**
- Oxygen stress (SWOXYGEN, HLIM1, HLIM2U/L, Q10_MICROBIAL, etc.; ~8 fields)
- Drought stress (SWDROUGHT, HLIM3H/L, HLIM4, ADCRH/L; ~5 fields)
- Salinity stress (SWSALINITY, SALTMAX, SALTSLOPE, SALTHEAD; ~4 fields)
- Compensation (SWCOMPENSATE, SWSTRESSOR, ALPHACRIT, DCRITRTZ; ~4 fields)

**CO2 impact (1 field + 3 tables)**
- `SWCO2` switch
- `CO2AMAXTB`, `CO2EFFTB`, `CO2TRATB` (CO2 correction factors)

**Nitrogen & management (20+ fields; optional in Phase 4c)**
- `RELMF` (relative management factor)
- Nitrogen parameters (optional; ~12 fields for detailed nitrogen model)
- Harvest losses (3 fields)
- Crop residues (1 field)

**Estimated LoC for cropwofost_config_t:**
- Type definition: ~400 lines (many tables, ~100 allocatable fields)
- TOML reader: ~500 lines
- **Total: ~900 lines**

**Note:** Phase 4c likely scopes **out WOFOST initially** and covers only types 1 and 3 (fixed & grass), as full WOFOST support requires additional infrastructure (table interpolation, nitrogen module). Revisit in Phase 4d.

## Current parity test gap

The parity test suite (`tests/swap-cases/` + `tests/unit/io/toml/test_read_*.pf`) currently:
1. ✓ Verifies rotation metadata reads correctly (swcrop, start/end dates, file, type)
2. ✗ **Does not verify per-crop data** — no assertions on crop-specific global variables (LAI tables, root density, stress thresholds, etc.)

To achieve full parity in Phase 4c:
- After loading `crop.rotation[N].file`, the reader must load the referenced `.crp.toml` file and populate the corresponding crop type config.
- The parity test must assert that globals populated by legacy `readcropfixed`, `readcropgrass`, `readcropwofost` subroutines match the fields in the new `cropfixed_config_t`, `cropgrass_config_t`, `cropwofost_config_t`.
- **Critical issue:** Legacy crop readers populate local variables (read once at init); the new path must also populate them to maintain in-memory behavior during simulation. Verify with `src/io/readswap.f90` line ~530 how globals are used downstream.

## Phase 4c parity test extension

For each case's crop rotation:

```
For each rotation entry (i = 1 .. num_rotations):
  - Load crop.rotation[i].file + type
  - If type = 1: load cropfixed_config_t, check vs. legacy readcropfixed globals
  - If type = 2: load cropwofost_config_t, check vs. legacy readcropwofost globals
  - If type = 3: load cropgrass_config_t, check vs. legacy readcropgrass globals
  - Example assertion:
    assert(config%gctb(:) == variables%gctb(:), "LAI/cover table mismatch")
    assert(config%rdtb(:) == variables%rdtb(:), "root density table mismatch")
```

## Recommendations

1. **Phase 4c-a (cropfixed + cropgrass only):**
   - Scopes out WOFOST to unblock cases 1, 2, 4, 5.
   - Estimated effort: ~1200 LoC (new types + readers) + ~400 LoC (test extensions).
   - Covers 5 non-macropore cases fully.

2. **Phase 4c-b (WOFOST; later):**
   - Adds cropwofost_config_t + reader.
   - Estimated effort: ~900 LoC + test coverage.
   - Covers case 6.surfacewater if it uses type 2; none of the test cases actually use type 2 in the current crop rotation (case 1 has type 2 but it's the `potatod` WOFOST variant, which is uncommon).

3. **Case 3.macroporeflow (excluded from Phase 4c):**
   - Uses `wintcer1.crp`, `wintcer2.crp` (both type 1).
   - Excluded from parity in Phase 4c per scope; its `.crp` files can be converted in Phase 4c or later if needed.

## Files to create/modify in Phase 4c

1. **`src/config/cropfixed_config.f90`** (new; ~200 lines)
2. **`src/io/toml/read_cropfixed_toml.f90`** (new; ~300 lines)
3. **`src/config/cropgrass_config.f90`** (new; ~300 lines)
4. **`src/io/toml/read_cropgrass_toml.f90`** (new; ~400 lines)
5. **`src/config/crop_config.f90`** (modify; add per-crop config refs)
6. **`src/io/toml/read_crop_toml.f90`** (modify; add cross-file loading)
7. **`tests/swap-cases/toml/*/crop-*.toml`** (new; one per `.crp` file → TOML conversion)
8. **`tests/unit/io/toml/test_read_cropfixed_toml.pf`** (new; unit tests)
9. **`tests/unit/io/toml/test_read_cropgrass_toml.pf`** (new; unit tests)
10. **`tests/parity/test_crop_parity.F90`** (new or extend; parity assertions)

## Estimated total LoC: Phase 4c-a (cropfixed + cropgrass)

- Config types: ~500 lines (cropfixed) + ~300 lines (cropgrass) = 800 lines
- TOML readers: ~300 lines (cropfixed) + ~400 lines (cropgrass) = 700 lines
- Test cases: ~500 lines (new TOML files for 5 cases, ~100 lines each) + ~400 lines (unit/parity tests)
- **Total: ~2400 lines** (mostly TOML file data, not complex logic)

**Timeline estimate: 2–3 weeks for Phase 4c-a (design, implementation, review, parity testing).**
