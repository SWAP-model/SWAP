# Phase 4e Task C1: Macroporeflow Regression Case Audit

## 1. Scope

Case 3 (`3.macroporeflow`) is a soil water simulation exercising **macropore flow physics**: preferential rapid flow through soil cracks and voids under saturated conditions. The case was deferred in Phase 4c-a because the existing `soil_config_t` schema (as of Phase 4d) lacks a `macropore_config_t` subtype — Phase 4f-prep will add this. This audit identifies case-3-specific inputs so Tasks C2 (TOML authoring) and C3 (parity test) can proceed with deferred macropore keys marked `# Phase 4f-prep:`.

## 2. .swp File Shape Comparison

| Section | Case 1 (hupselbrook) | Case 3 (macroporeflow) | Schema-Covered | Notes |
|---------|----------------------|------------------------|----------------|-------|
| **General** | PROJECT=hupsel | PROJECT=andelst | Y | Standard |
| **Simulation** | TSTART=2002-01-01, TEND=2004-12-31 | TSTART=1998-01-01, TEND=1999-04-26 | Y | Case 3: shorter (~1.3 yr vs 3 yr) |
| **Meteorology** | SWETR=0 (Penman-Monteith), SWRAIN=0 | SWETR=1 (ref ET), SWRAIN=3 (detailed rainfall) | Partial (G) | Case 3: requires separate rainfall file; SWRAIN=3 not yet in schema |
| **Crop rotation** | 3 entries (maizes, potatod, grassd) | 2 entries (wintcer1, wintcer2; both type 1) | Y | Case 3: simpler (winter cereals only; Type 1 = simple growth) |
| **SWCROP** | 1 | 1 | Y | Both cultivated |
| **Soil discretization** | 2 layers, 6 compartments (10cm each) | 8 layers, 94 compartments (1–5cm) | Y | Case 3: finer profile (macropore domain requires resolution) |
| **Soil hydraulics** | SWSOPHY=0 (van Genuchten) | SWSOPHY=0 (van Genuchten) | Y | Both analytical |
| **Hysteresis** | SWHYST=0 | SWHYST=0 | Y | None in either case |
| **SWMACRO** | 0 (no macropores) | **1 (macropores enabled)** | N (G) | **Key trigger for case 3** |
| **Bottom boundary** | SWBBCFILE=1 → swap.bbc | SWBBCFILE=1 → swap.bbc | Y | Standard file reference |
| **Heat** | SWHEA=1, SWCALT=2 (numerical) | SWHEA=1, SWCALT=2 (numerical) | Y | Both heat transport |
| **Solute** | SWSOLU=1 | SWSOLU=0 | Y | Case 1 includes solute; case 3 does not |

## 3. .dra File Shape Comparison

| Field | Case 1 | Case 3 | Schema-Covered | Notes |
|-------|--------|--------|----------------|-------|
| DRAMET | 2 (Hooghoudt/Ernst formula) | 3 (resistance multi-level) | Y | Case 3 uses more flexible resistance method |
| SWDIVD | 1 (vertical distribution) | 1 (vertical distribution) | Y | Standard |
| NRLEVS | 1 | 1 | Y | Both single drainage level |
| COFANI | 1.0 (× 2 layers) | 1.0 (× 8 layers) | Y | Anisotropy factor per soil layer |
| DRARES1 / INFRES1 | DRARES=500, INFRES=220 (read from separate levels) | DRARES=500, INFRES=220 | Y | Identical resistance values |
| L1 (drain spacing) | ~computed from formula | 10.0 m | Y | Case 3 explicit |
| ZBOTDR1 | -85.0 cm (inferred) | -85.0 cm | Y | Drainage level depth |

## 4. Crop File (.crp) Shape

**Case 3 rotation entries:**
1. `wintcer1.crp` (1998-01-01 to 1998-08-20, Type 1 = simple growth)
2. `wintcer2.crp` (1999-01-01 to 1999-04-26, Type 1 = simple growth)

Both are **Type 1 (simple growth)** — read via `readcropfixed()`. Comparison with case 1's same-type crop (maizes.crp, Type 1):
- Standard Type 1 fields: CROPSTART, CROPEND, LAI growth table, root depth, soil water stress factors.
- **No deviation expected** — winter cereal and maize are both Type 1; schema covers all fields.

## 5. Macropore-Specific Keys (All G — Phase 4f-prep Blockers)

The `.swp.template` lines **313–388** define 13 input keys with no current config slot:

| Legacy Name (variables.f90) | .swp Value | Units | Default | Trigger | Phase 4f Note |
|-----|-----|-----|-----|-----|-----|
| SWMACRO | 1 | — | 0 | Enables all below | needs `macropore_config_t.swmacro` |
| Z_AH | -26.0 | cm (absolute depth) | — | A-horizon bottom | needs `macropore_config_t.z_ah` |
| Z_IC | -100.0 | cm | — | IC domain bottom | needs `macropore_config_t.z_ic` |
| Z_ST | -220.0 | cm | — | Static macropore bottom | needs `macropore_config_t.z_st` |
| VLMPSTSS | 0.04 | cm³/cm³ | 0 | Static macropore volume @ surface | needs `macropore_config_t.vlmpstss` |
| PPICSS | 0.5 | — (fraction) | — | IC domain proportion @ surface | needs `macropore_config_t.ppicss` |
| NUMSBDM | 4 | — (count) | — | IC subdomains | needs `macropore_config_t.numsbdm` |
| POWM | 1.0 | — | 1.0 (optional) | Freq. dist. power | needs `macropore_config_t.powm` |
| RZAH | 0.0 | — (fraction) | 0.0 (optional) | Macropores ending @ A-horizon | needs `macropore_config_t.rzah` |
| SPOINT | 1.0 | — | 1.0 (optional) | Symmetry point | needs `macropore_config_t.spoint` |
| SWPOWM | 0 | Y/N | 0 (optional) | Double convex/concave switch | needs `macropore_config_t.swpowm` |
| DIPOMI | 10.0 | cm | — | Min polygon diameter | needs `macropore_config_t.dipomi` |
| DIPOMA | 50.0 | cm | — | Max polygon diameter | needs `macropore_config_t.dipoma` |

**Additional G keys (lines 330–388, soil profile modifiers for macropores):**

| Legacy Name | Rows/Scope | Schema Gap | Phase 4f Note |
|-----|-----|-----|-----|
| SWSOILSHR | Table (8 rows, per soil layer) | 0 | needs `macropore_config_t.swsoilshr` |
| SWSHRINP | Table (8 rows) | 0 | needs `macropore_config_t.swshrinp` |
| THETCRMP | Table (8 rows) | 0 | needs `macropore_config_t.thetcrmp` |
| GEOMFAC | Table (8 rows) | 0 | needs `macropore_config_t.geomfac` |
| SHRPARA...SHRPARE | Table (8 rows × 5 cols) | 0 | needs `macropore_config_t.shrpar_*` |
| SWSORP | Table (8 rows) | 0 | needs `macropore_config_t.swsorp` |
| SORPFACPARL | Table (8 rows) | 0 | needs `macropore_config_t.sorpfacparl` |
| SORPMAX | Table (8 rows) | 0 | needs `macropore_config_t.sorpmax` |
| SORPALFA | Table (8 rows) | 0 | needs `macropore_config_t.sorpalfa` |
| ZNCRAR | Scalar | 0 | needs `macropore_config_t.zncrar` |
| SHAPEFACMP | Scalar (1.5) | 0 | needs `macropore_config_t.shapefacmp` |
| CRITUNDSATVOL | Scalar (0.1) | 0 | needs `macropore_config_t.critundsatvol` |
| SWDARCY | Scalar (1) | 0 | needs `macropore_config_t.swdarcy` |
| SWDRRAP | Scalar (1) | 0 | needs `macropore_config_t.swdrrap` (rapid drainage) |
| RAPDRARESREF | Scalar (50.0) | 0 | needs `macropore_config_t.rapdraresref` |
| RAPDRAREAEXP | Scalar (1.0) | 0 | needs `macropore_config_t.rapdrareaexp` |
| NUMLEVRAPDRA | Scalar (1) | 0 | needs `macropore_config_t.numlevrapdra` |
| PNDMXMP | Scalar (0.0) | 0 | needs `macropore_config_t.pndmxmp` |

**Total unique G macropore keys: 22** (13 scalar + 9 table params per layer).
Per the variable audit, the **Macropore section has 10 confirmed G entries** (DiPoMa, DiPoMi, GeomFac, PowM, ShapeFacMp, swman, SwPowM, SwSoilShr, SwSorp, VlMpStSs); the additional 12 are likely reclassified as R or covered by runtime simulation (e.g., depth tables like SHRPAR are per-layer, not direct G inputs).

## 6. Schema-Covered Subset

**Fully or Partially Covered (% of case-3 .swp):**

- **General** (PROJECT, PATHWORK, PATHATM, PATHCROP, PATHDRAIN, SWSCRE, SWERROR) — 100% ✓
- **Simulation** (TSTART, TEND, NPRINTDAY, SWMONTH, PERIOD, SWRES, SWODAT, SWYRVAR, DATEFIX, OUTFIL) — 100% ✓
- **Meteorology** (METFIL, LAT, SWETR, SWETSINE) — 80% ✓; **SWRAIN=3 (G)**, RAINFIL (stubbed to SWRAIN=2)
- **Crop** (SWCROP, rotation table, RDS) — 100% ✓
- **Soil** (SWSOPHY, hydraulic parameters by layer, SWHYST, SWINCO, GWLI, PONDMX, vertical discretization) — 85% ✓; **SWMACRO=1 + macropore params (G)**
- **Drainage** (SWDRA, DRFIL, .dra.toml) — 100% ✓
- **Bottom Boundary** (SWBBCFILE, BBCFIL) — 100% ✓ (Phase 4d)
- **Heat** (SWHEA, SWCALT, soil texture, TSOIL init, SWTOPBHEA, SWBOTBHEA) — 100% ✓ (Phase 4d)
- **Solute** (SWSOLU) — 0% (SWSOLU=0 in case 3; no solute)

**Estimate: ~75–80% of the .swp is schema-covered** (excluding macropore section, SWRAIN=3, and inactive solute). The macropore block is ~50 lines out of 514; that's ~10% of the file but 100% of the G gaps.

## 7. Authoring Strategy for Tasks C2/C3

1. **TOML template:** Mirror case 1 / case 4 structure (already exists at `/home/zawadzkim/Code/swap/tests/swap-cases/toml/3.macroporeflow/swap.toml`), but **update crop.rotation entries to wintcer1/wintcer2 (Type 1)** and **update timing to TSTART=1998-01-01, TEND=1999-04-26**. Case-3-specific deviations: SWETR=1, SWRAIN=2 (stubbed from 3), DRAMET=3, SWMACRO=1.

2. **Macropore section:** **Omit entirely from .toml.** Add a single comment block at the file top:
   ```toml
   # Phase 4f-prep: SWMACRO=1 in legacy .swp; the [macropore] section is deferred.
   # The following parameters exist in swap_linux.swp.template (lines 313–388)
   # but have no schema slot until macropore_config_t lands:
   # Z_AH, Z_IC, Z_ST, VLMPSTSS, PPICSS, NUMSBDM, POWM, RZAH, SPOINT, SWPOWM,
   # DIPOMI, DIPOMA, shrinkage curves, sorptivity, rapid drainage params, etc.
   ```

3. **swap.dra.toml:** Standard pattern — DRAMET=3 (resistance multi-level), NRLEVS=1, single drainage level with DRARES=500, INFRES=220, L=10, ZBOTDR=-85.

4. **.crp.toml files:** `wintcer1.crp.toml` and `wintcer2.crp.toml` (both Type 1). Standard simple-growth template; no deviations from case 1's maizes.crp.

5. **Parity test scope:** Assert ~40–50 fields across [general], [simulation], [meteorology], [crop], [soil] (non-macropore), [drainage], [bottom_boundary], [heat]. **Do NOT assert macropore keys** — they have no config side yet. The test verifies round-trip (legacy .swp → TOML input → variables % state matches legacy readswap output), excluding the macropore block.

6. **check-full behavior:** Will continue using legacy `readswap()` for case 3 anyway (no Phase 4f entry point yet); TOML is surface for parity test only.

## 8. Files & Paths Summary

**Tasks C2/C3 will create or update:**

- `/home/zawadzkim/Code/swap/tests/swap-cases/toml/3.macroporeflow/swap.toml` (Task C2)
  - Update [simulation] dates, [meteorology] (SWETR=1, SWRAIN=2), [drainage] (DRAMET=3)
  - Keep [soil] with swmacro=1, comment macropore params defer
  - Update [crop.rotation] to wintcer1/wintcer2 Type 1
  
- `/home/zawadzkim/Code/swap/tests/swap-cases/toml/3.macroporeflow/swap.dra.toml` (Task C2)
  - DRAMET=3, NRLEVS=1, DRARES=500, INFRES=220, L=10, ZBOTDR=-85

- `/home/zawadzkim/Code/swap/tests/swap-cases/toml/3.macroporeflow/wintcer1.crp.toml` (Task C2)
- `/home/zawadzkim/Code/swap/tests/swap-cases/toml/3.macroporeflow/wintcer2.crp.toml` (Task C2)
  - Standard Type 1 templates; no case-3-specific overrides

- `/home/zawadzkim/Code/swap/tests/unit/io/toml/test_macroporeflow_parity.pf` (Task C3)
  - Assert TSTART, TEND, meteorology (swetr, lat, swrain=2), crop rotation (2 entries), soil (swmacro, hydraulics, discretization), drainage (dramet, nrlevs, drares, infres, L, zbotdr).
  - Skip all macropore-specific assertions (z_ah, vlmpstss, shrinkage, sorptivity, rapid drainage).

---

**Summary for Phase 4f-prep:**

- Case 3 regression enables **13–22 macropore input keys**, all of which are **G (gaps)** in the current schema.
- The macropore block spans ~10% of the .swp file but drives ~90% of case-3 specificity.
- Phase 4f must introduce `macropore_config_t` with slots for:
  - Depth zone markers (z_ah, z_ic, z_st)
  - Volume and proportion fractions (vlmpstss, ppicss)
  - IC domain parameters (numsbdm, powm, rzah, spoint, swpowm)
  - Polygon geometry (dipomi, dipoma)
  - Per-layer shrinkage and sorptivity curves (8 rows × 5–9 cols)
  - Rapid drainage switches and resistances
  - Darcy/macropore shape factor parameters

Once macropore_config_t lands, case 3's TOML can be completed with the deferred macropore section, and the parity test can extend to full coverage.
