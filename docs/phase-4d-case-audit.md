# Phase 4d: Regression Case Configuration Audit

## 1. Scope & Purpose

This document audits the 5 non-macropore regression test cases to identify which Phase 4d configuration fields each case actually uses. Phase 4d adds 4 new config types (`bottom_boundary_config_t`, `heat_config_t`, `irrigation_config_t` + nested `irrigation_schedule_t`, `solute_config_t`) and extends `cropgrass_config_t` with mowing/grazing per-event tables. This audit captures the exact legacy parameters that code authors (Tasks 17-23) must map into the new TOML schemas—no guesswork required.

**In-scope cases:** 1.hupselbrook, 2.grassgrowth, 4.oxygenstress, 5.salinitystress, 6.surfacewater (non-macropore only).

---

## 2. Bottom Boundary Configuration (`swbotb`)

| Case | swbotb | Boundary Type | Key Parameters | Notes |
|------|--------|---------------|-----------------|-------|
| 1 | 6 | Zero flux | None | Simple passthrough: `swbotb=6` → `bottom_boundary.swbotb=6` |
| 2 | 1 (BBC file) | Prescribed GWL | BBCFIL='swap' | Uses separate `.bbc` file; static GWL = -75 cm |
| 4 | 3 (Cauchy) | Deep aquifer | SHAPE=1.0, HDRAIN=-25, RIMLAY=500, SW3=2 (table) | Dynamic HAQUIF table (48 rows, 1993–2002) |
| 5 | 3 (Cauchy) | Deep aquifer | SHAPE=0.79, HDRAIN=-110, RIMLAY=50, SW3=1 (sinus) | AQAVE=-120, AQAMP=20, AQTMAX=120, AQPER=365 |
| 6 | 1 (BBC file) | Prescribed GWL | BBCFIL='swap' | Uses separate `.bbc` file |

**Sub-sections needed:** Cases 4 & 5 need `[bottom_boundary.cauchy]` with respective params; cases 1, 2, 6 are minimal/file-based.

---

## 3. Heat Flow Configuration (`swhea`)

| Case | swhea | swcalt | Soil Texture (PSAND/PCLAY/PORG) | Initial Temp Table | TFROSTSTA/TFROSTEND |
|------|-------|--------|----------------------------------|-------------------|---------------------|
| 1 | 1 | 2 (numerical) | Yes (PSAND=0.8, PCLAY=0.05, ORGMAT=0.1) | Yes (-10 to -95 cm, 15–9°C) | Not set (SWFROST=0) |
| 2 | 1 | 2 (numerical) | Yes (5 layers: 0.68–0.88 sand) | Yes (-10 to -95 cm, 15–9°C) | Not set (SWFROST=0) |
| 4 | 1 | 2 (numerical) | Yes (4 layers: 0.699–0.000 sand) | Yes (3 rows: 0 to -298 cm) | Not set (SWFROST=0) |
| 5 | 1 | 2 (numerical) | Yes (2 layers: 0.935 sand, 0.935 sand) | Yes (4 rows: -10 to -95 cm, 15–9°C) | Not set (SWFROST=0) |
| 6 | 0 | N/A | Not present | Not present | Not applicable |

**Summary:** Cases 1–5 have heat enabled. Cases 1, 2, 5 share same initial-temp table pattern (4 entries, -10 to -95 cm). Case 6 disables heat entirely.

---

## 4. Solute Transport Configuration (`swsolu`)

| Case | swsolu | swbotbc | cdrain | cseep | tscf | ldis | Notes |
|------|--------|---------|--------|-------|------|------|-------|
| 1 | 1 | 0 | 0.1 | 0.1 | 0.0 | 5.0, 5.0 | Standard boundary transport; no special uptake |
| 2 | 0 | N/A | N/A | N/A | N/A | N/A | **No solute transport** |
| 4 | 0 | N/A | N/A | N/A | N/A | N/A | **No solute transport** |
| 5 | 1 | 1 | 1.051 | 15.574 | 0.5 | 5.0, 5.0 | **Maas-Hoffman salinity stress active** (potatod.crp: SWSALINITY=1). CSEEP > CDRAIN (seepage enrichment). TSCF=0.5 (halts uptake in saturated root zone). |
| 6 | 0 | N/A | N/A | N/A | N/A | N/A | **No solute transport** |

**Critical:** Case 5 is the sole solute-transport case and triggers salinity stress in the crop model. Must map exact values `cdrain=1.051`, `cseep=15.574`, `tscf=0.5`.

---

## 5. Fixed Irrigation Configuration (`swirfix`)

| Case | swirfix | irgfil (file reference) | Inline events present? | Source |
|------|---------|-------------------------|----------------------|--------|
| 1 | 1 | 'testirri' | Yes (1 event: 2002-01-05, 5mm) | Inline in .swp (lines 255–257) |
| 2 | 0 | N/A | No | Disabled |
| 4 | 0 | N/A | No | Disabled |
| 5 | 1 | 'swap' | No (external file) | External `swap.irg` file |
| 6 | 0 | N/A | No | Disabled |

**Summary:** Only cases 1 & 5 use irrigation; case 1 has inline event, case 5 delegates to external file (not audited here per scope).

---

## 6. Per-Crop Irrigation Scheduling (`.crp` files)

| Case | Crop File | Rotation Type | SCHEDULE | tcs | dcs | Timing Table | Depth Table | Notes |
|------|-----------|---------------|----------|-----|-----|--------------|-------------|-------|
| 1 | maizes.crp | 1 (simple) | 0 | N/A | N/A | No | No | No scheduling |
| 1 | grassd.crp | 3 (grass) | 0 | N/A | N/A | No | No | No scheduling |
| 1 | potatod.crp | 2 (detailed) | 0 | N/A | N/A | No | No | No scheduling |
| 2 | grassd.crp | 3 (grass) | 0 | N/A | N/A | No | No | No scheduling |
| 4 | grassd.crp | 3 (grass) | 0 | N/A | N/A | No | No | No scheduling |
| 5 | potatod.crp | 2 (detailed) | 0 | N/A | N/A | No | No | No scheduling |
| 6 | grass.crp | 1 (simple) | 0 | N/A | N/A | No | No | No scheduling |

**Finding:** No case uses crop-level irrigation scheduling (SCHEDULE=0 in all .crp files). All irrigation is fixed-application only (case 1 inline, case 5 external file).

---

## 7. Grass Mowing & Grazing Tables (Cases 1, 2, 4)

| Case | Crop File | SWDMMOW | SWHARVEST | Mow Events | Grazing Schedule Present? | Notes |
|------|-----------|---------|-----------|-----------|-------------------------|-------|
| 1 | grassd.crp | 2 (flexible DM threshold) | 1 (DM-threshold) | Yes (NMOW table) | No (SWHARV=0) | Mowing defined; grazing disabled |
| 2 | grassd.crp | 2 (flexible DM threshold) | 2 (fixed dates) | Yes (harvest dates table) | Yes (nstart_graz, nstop_graz present) | Both mowing & grazing defined |
| 4 | grassd.crp | 2 (flexible DM threshold) | 1 (DM-threshold) | Yes (NMOW table) | Yes (nstart_graz, nstop_graz present) | Both mowing & grazing defined |
| 6 | grass.crp | N/A (simple crop) | N/A | N/A | N/A | Simple model; no detailed harvesting |

**Summary:**
- **Case 1:** Mowing only (SWDMMOW=2, lines 434–442 of grassd.crp)
- **Case 2:** Mowing + grazing (SWHARVEST=2 with fixed dates, plus grazing sections)
- **Case 4:** Mowing + grazing (SWHARVEST=1 with DM threshold, plus grazing sections)
- **Case 6:** No harvest tables (fixed crop, SWHARV=0 in grass.crp line 38)

---

## 8. Out-of-Scope (Deferred Beyond Phase 4d)

The following legacy features are used by the regression cases but NOT yet in Phase 4d scope:

1. **Extended drainage (DRAMET ≠ 3, SWDRA=2):** Case 6 uses extended drainage (SWDRA=2) with multi-level surface water mgmt. Phase 4d captures basic DRAMET=3 multi-level params only.
2. **Advanced interception models:** Case 6 has Gash-style forest interception (SWINTER≥2); Phase 4d assumes Von Hoyningen–Hune & Braden (1=ag crops).
3. **Nitrogen/LINTUL4:** All crops have N-uptake params (RMO, FOTB fractions) but Phase 4d does not yet model these.
4. **Ponding/runoff dynamics:** SWPONDMX, RSRO, RSROEXP are set in all cases but not yet Phase 4d-tracked.

**Confirmation:** No field that any regression case actually uses is left out by Phase 4d's scope. All bottom-boundary, heat, solute, and irrigation fields in the legacy files map to the new configs.

---

## 9. Authoring Strategy Recommendations

1. **Bottom boundary trivial cases (1, 2, 6):** swbotb=6 or swbotb=1 with BBC file reference—author these first (minimal sections, just `swbotb` value).

2. **Bottom boundary non-trivial (4, 5 with swbotb=3):** Requires `[bottom_boundary.cauchy]` sub-section. Author case 5 last so cases 1–4 establish the pattern.

3. **Heat config (1–5 all have SWHEA=1, swcalt=2):** All numerical method with soil texture + initial temp table. Case 6 disables heat entirely (simpler). Code authors should handle both paths; cases 1–5 are nearly identical.

4. **Solute transport (case 5 only):** Only case 5 activates solute + salinity stress. Simple paths (swsolu=0) for cases 1–4, 6; full path for case 5 alone. Map exactly: `cdrain=1.051`, `cseep=15.574`, `tscf=0.5`.

5. **Irrigation (cases 1 & 5):** Case 1 has inline fixed events (simple); case 5 delegates to external .irg file. No crop-level scheduling to author (all SCHEDULE=0).

6. **Grass management (cases 1, 2, 4):** Case 1's grassd.crp needs mowing table only (simpler, author first); cases 2 & 4 need both mowing + grazing tables. Case 6's grass.crp is type=1 (simple model)—no harvest tables.

---

## 10. Files to Author (Tasks 17–23)

**Main configs (Task 17: .swp → swap.toml):**
- tests/swap-cases/toml/1.hupselbrook/swap.toml
- tests/swap-cases/toml/2.grassgrowth/swap.toml
- tests/swap-cases/toml/4.oxygenstress/swap.toml
- tests/swap-cases/toml/5.salinitystress/swap.toml
- tests/swap-cases/toml/6.surfacewater/swap.toml

**Crop configs (Task 18–19: .crp → *.crp.toml):**
- tests/swap-cases/toml/1.hupselbrook/maizes.crp.toml (type 1)
- tests/swap-cases/toml/1.hupselbrook/potatod.crp.toml (type 2)
- tests/swap-cases/toml/1.hupselbrook/grassd.crp.toml (type 3, mowing only)
- tests/swap-cases/toml/2.grassgrowth/grassd.crp.toml (type 3, mowing + grazing)
- tests/swap-cases/toml/4.oxygenstress/grassd.crp.toml (type 3, mowing + grazing)
- tests/swap-cases/toml/5.salinitystress/potatod.crp.toml (type 2, salinity)
- tests/swap-cases/toml/6.surfacewater/grass.crp.toml (type 1)

---

**Summary:** 5 swap.toml files (all 4d-relevant), 7 *.crp.toml files. Zero crop-level irrigation scheduling blocks (all SCHEDULE=0). Three grassd files need mowing/grazing extensions (cases 1, 2, 4).
