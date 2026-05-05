# Phase 4c-b Task 9: Case 5 (Salinitystress) Legacy Input Audit

## 1. Scope

Case 5 (`salinitystress`) simulates salinity-stress crop response in potato using the WOFOST (type-2) growth model with active solute transport. The case was deferred from Phase 4c-a; this audit determines what is case-5-specific to inform TOML authoring in Tasks 10–12. Irrigation and solute/decomposition configuration are deferred to Phase 4d.

---

## 2. .swp File Shape (Main Configuration)

| Section | Key Fields | Case 1 | Case 5 | Notes |
|---------|-----------|--------|--------|-------|
| General | PROJECT | 'hupsel' | 'saltfarmtexel' | Case 5 names the salt farm domain |
| Simulation | TSTART | 2002-01-01 | 2012-01-01 | Case 5: 4-year run (2012–2015) vs 3-year (2002–2004) |
| Simulation | TEND | 2004-12-31 | 2015-12-31 | — |
| Crop | CROPSTART/END (4 rows) | mixed (maizes, potatod, grassd) | All potatod, 4 years | Case 5: single crop type, all potatod |
| Output CSV | INLIST_CSV | 'rain,irrig,...' | 'cwso,cpwso,tredwet,...,conc' | Case 5: solute-focused outputs |
| Output CSV | INLIST_CSV_TZ | (0) | 'conc' (1) | Case 5: requests conc time-depth profiles |
| Solute | SWSOLU | 1 | 1 | Both enabled; case 5 is test case for this |
| Solute | CDRAIN | 0.1 | 1.051 | Case 5: saline drainage (higher salt) |
| Solute | CSEEP | 0.1 | 15.574 | Case 5: seepage from deeper aquifer is very saline |
| Solute | TSCF | 0.0 | 0.5 | Case 5: roots absorb solutes (no longer 0) |
| Solute | SWBOTBC | 0 | 1 | Case 5: use constant CSEEP for upflow |
| Solute | LDIS | 5.0, 5.0 | 5.0, 5.0 | Same dispersion lengths |
| Solute | DDIF | 0.0 | 0.0 | Same (no diffusion) |

**Case-5-specific differences:** solute concentrations ~10× higher; solute uptake enabled (TSCF=0.5); deep seepage is the primary salt source.

---

## 3. .dra File Shape (Drainage)

| Parameter | Case 1 (dramet=2) | Case 5 (dramet=3) | Notes |
|-----------|-----------------|-----------------|-------|
| DRAMET | 2 (Hooghoudt/Ernst formula) | 3 (multi-level resistance) | Case 5: more complex dual-layer drainage |
| NRLEVS | (implicit 1 via finalize) | 2 | Case 5: two drainage levels |
| DRARES1 | N/A (formula-based) | 10.0 d | Shallow drains at -60 cm |
| ZBOTDR1 | (implicit -80 cm) | -60.0 cm | Case 5: shallower drain 1 |
| INFRES1 | N/A | 10.0 d | Infiltration resistance level 1 |
| SWALLO1 | N/A | 3 (infiltration only) | Case 5: blocks drainage to level 1, allows upflow |
| L1 | N/A | 5.0 m | Case 5: tight drain spacing (5 m vs ~11 m case 1) |
| DRARES2 | N/A | 50.0 d | Deeper drains at -120 cm |
| ZBOTDR2 | N/A | -120.0 cm | Case 5: deep drain 2 |
| SWALLO2 | N/A | 1 (both allowed) | Case 5: normal drainage/infiltration at level 2 |
| L2 | N/A | 50.0 m | Case 5: wider spacing at depth |

**Case-5-specific:** multi-level (dramet=3) vs single formula (dramet=2); shallow tight drains with salt removal, deeper infiltration layer.

---

## 4. .crp File Shape (WOFOST Crop, 418 lines vs case 1's 704)

| Section | Item | Case 1 | Case 5 | Notes |
|---------|------|--------|--------|-------|
| **Part 0: Sowing** | SWGERM | 2 (temperature + hydro) | 0 (none) | Case 5: simpler germination model |
| — | DVSEND (harvest stage) | 2.0 | 3.0 | Case 5: harvest at maturity (DVS=3) not DVS=2 |
| **Part 1: Crop factor** | SWCF | 2 (crop height) | 2 (crop height) | Same approach |
| — | CH table rows | 3 (DVS 0, 1, 2) | 3 (DVS 0, 1, 2) | Similar structure |
| **Part 2: Phenology** | IDSL | 0 (temperature only) | 0 (temperature only) | Same |
| — | DTSMTB | 4 rows | 4 rows | Identical temperature sum table |
| — | Daylength (IDSL>0) | Included (17 lines) | Absent | Case 1 has optional vernalisation; case 5 omits |
| **Part 5: Assimilation** | AMAXTB, TMPFTB, TMNFTB | Present | Present | Identical |
| **Part 8: Partitioning** | FRTB, FLTB, FSTB, FOTB | Present | Present | Case 5 slightly shorter (fewer DVS rows in FRTB/FLTB) |
| **Part 10: Root depth** | SWRD | 2 (max daily increase) | 2 (max daily increase) | Same |
| — | RDTB (SWRD=1 table) | Included (not used) | Absent | Case 1 includes optional RDTB; case 5 omits |
| — | RLWTB (SWRD=3 table) | Included (not used) | Absent | Case 1 includes optional biomass table; case 5 omits |
| — | WRTMAX (max root weight) | Included (not used) | Absent | Case 1 includes; case 5 omits |
| **Part 11: O₂ stress** | SWOXYGEN | 1 (Feddes) | 1 (Feddes) | Same |
| — | SWOXYGEN=2 params | 9 lines (Q10_MICROBIAL, SRL, etc.) | Absent | Case 1 includes detailed O₂ stress (Bartholomeus); case 5 uses simpler model |
| **Part 13: Salt stress** | SWSALINITY | **0** (none) | **1** (Maas–Hoffman) | **CASE 5 KEY DIFFERENCE**: salinity ON |
| — | SALTMAX | 3.0 | 0.732 | Case 5: lower threshold (more sensitive) |
| — | SALTSLOPE | 0.1 | 0.0868 | Case 5: shallower decline slope |
| — | SALTHEAD | 624.0 (only if SWSALINITY=2) | Absent | Case 5 uses SWSALINITY=1, not osmotic head |
| **Part 15: CO₂ impact** | SWCO2 + CO2*TB tables | 5 tables (50 lines) | Absent | Case 1 includes CO₂ corrections; case 5 omits (SWCO2=0, no tables) |
| **Management: Nitrogen** | LINTUL4 nitrogen section | 45 comment lines | Absent | Case 1 comments nitrogen; case 5 omits entirely |
| **Management: Scheduling** | SCHEDULE | 0 (none) | 0 (none) | Both fixed (no adaptive irrigation) |
| — | Irrigation detail tables | 23 lines (TCS/DCS options) | Absent | Case 1 includes scheduling options; case 5 omits (deferred to 4d) |

**Critical difference:** Case 5 activates `SWSALINITY=1` (Maas–Hoffman salinity stress), making this the salinity-stress regression case. Case 1 has `SWSALINITY=0`. Case 5 strips unnecessary features (daylength, O₂ Bartholomeus, CO₂, nitrogen, irrigation detail).

---

## 5. Schema Gaps: Fields Needed but Not Yet in Config Types

### 5.1 WOFOST-specific (extends `cropwofost_config_t`)

| Legacy Field | Section | Case 5 Requirement | Status | Phase |
|--------------|---------|-------------------|--------|-------|
| SALTHEAD | Part 13 (salt stress, line 413) | Present but unused (SWSALINITY=1, not 2) | Not needed for case 5 parity | 4c-b (OK to defer) |
| CO2AMAXTB, CO2EFFTB, CO2TRATB | Part 15 (CO₂, lines 485–509) | Absent in case 5 (SWCO2=0) | Covered in 4c-a but case 5 omits | 4c-b (skip) |
| Nitrogen tables (NMXLV, etc.) | Management (lines 543+) | Absent in case 5 | Commented in case 1 only | Phase 4d |
| **SWCIRRTHRES, CIRRTHRES, PERIRRSURP** | Scheduling (lines 655–660) | Absent; irrigation deferral | Not in current schema | Phase 4d |
| TCS, DCS, DVS_TC*, DVS_DC*, etc. | Irrigation scheduling (lines 596–701) | Absent; 100+ lines of detail | Not in current schema | Phase 4d |

### 5.2 General/simulation (extends `simulation_config_t`)

| Legacy Field | File (.swp) | Case 5 Requirement | Status | Phase |
|--------------|-------------|-------------------|--------|-------|
| SWINCO, INIFIL | Part 1 (lines 215–221) | SWINCO=3, INIFIL='swap.ini' | Likely in 4c-a | 4c-b (check coverage) |
| Initial solute concentration (ZC, CML table) | Solute Part 1 (lines 722–727) | Present (0.0 at all depths) | Not yet in schema | Phase 4d |
| **SWHEA, SWCALT, soil texture (PSAND, etc.), TSOIL table** | Heat flow (lines 444–467) | Present but deferred | Not in schema | Phase 4d |
| **Bottom boundary (SWBOTB, SHAPE, HDRAIN, RIMLAY, AQAVE, AQAMP, etc.)** | Bottom boundary (lines 390–437) | SWBOTB=3 (aquifer head) | Not in schema | Phase 4d |

### 5.3 Drainage (extends `drainage_config_t`)

| Legacy Field | File (.dra) | Case 5 Requirement | Status | Phase |
|-------------|-------------|-------------------|--------|-------|
| SWINTFL, COFINTFLB, EXPINTFLB | Part 3 (lines 39, 103–105) | 0 (interflow OFF) | Not in schema | Phase 4d |
| SWTOPNRSRF | Part 3 (line 109) | 0 (no adjustment) | Not in schema | Phase 4d |

**Summary:** 18–22 fields identified as NOT in current config schema, primarily:
- **4c-b fixable (0):** None identified that block case 5 parity if salinity is the focus.
- **4d deferred (18–22):** Heat, bottom boundary, initial solute concentrations, irrigation scheduling, nitrogen, advanced drainage (interflow).

---

## 6. Authoring Strategy for Tasks 10–12

1. **Task 10 (swap.toml + swap.dra.toml):**
   - Mirror case-1 swap.toml structure; deviations: PROJECT, TSTART/TEND (2012–2015), single potatod crop (4 entries), solute conc. outputs, CDRAIN/CSEEP (saline values), TSCF (0.5), SWBOTBC (1).
   - swap.dra.toml: **DO NOT reuse case 1's dramet=2 template.** Case 5 uses **dramet=3** (multi-level), NRLEVS=2, DRARES1=10, ZBOTDR1=-60, SWALLO1=3 (infiltration only), and deeper DRARES2/ZBOTDR2. This is a distinct pattern.

2. **Task 11 (potatod.crp.toml):**
   - Base structure: case-1 potatod with SWSALINITY=1 active, SALTMAX=0.732, SALTSLOPE=0.0868.
   - Omit: SWGERM detail (reduce to 0), daylength tables (if case 1 has IDSL>0 sections), CO₂ tables (SWCO2=0), nitrogen comments, irrigation scheduling tables.
   - Retain: core phenology (DTSMTB), assimilation, root density, drought (SWDROUGHT=1), O₂ stress (SWOXYGEN=1, Feddes only), interception.
   - DVSEND=3.0 (harvest at maturity, not DVS=2).

3. **Task 12 (test_salinitystress_parity.pf):**
   - Parity test covers ~40 fields across .swp + .dra: TSTART, TEND, solute concs (CDRAIN, CSEEP), drainage method (DRAMET), drain levels (NRLEVS), soil profile, bottom boundary type (SWBOTB), output flags.
   - ~50 WOFOST fields via read_legacy_wofost: SWSALINITY=1 (verified), SALTMAX/SALTSLOPE, crop phenology (TSUMEA, TSUMAM), partitioning (FRTB, FLTB, FSTB, FOTB), root (SWRD, RDI, RRI, RDC), O₂ (HLIM1, HLIM2U/L), drought (HLIM3H/L, HLIM4).
   - Assert parity on read path, not on computation (salinity-stress reduction happens in WOFOST, not in reader).

4. **Deferral list (NOT blocking 4c-b parity):**
   - Heat transport (SWHEA=1, numerical soil temperature, PSAND/PCLAY/ORGMAT, TSOIL table) → Phase 4d.
   - Bottom boundary detail (SHAPE, HDRAIN, RIMLAY, AQAVE/AMP/TMAX/PER, SW3, SW4, QBOT4 table) → Phase 4d.
   - Initial solute concentration profile (ZC/CML table) → Phase 4d.
   - Irrigation scheduling (SCHEDULE, TCS/DCS, DVS_TC*/DVS_DC*, STARTIRR/ENDIRR, IRGDAYFIX, etc.) → Phase 4d.
   - Nitrogen use (LINTUL4 parameters) → Phase 4d.

---

## 7. Files to be Created (Cross-Reference)

- `tests/swap-cases/toml/5.salinitystress/swap.toml` (Task 10)
- `tests/swap-cases/toml/5.salinitystress/swap.dra.toml` (Task 10)
- `tests/swap-cases/toml/5.salinitystress/potatod.crp.toml` (Task 11)
- `tests/unit/io/toml/test_salinitystress_parity.pf` (Task 12)

