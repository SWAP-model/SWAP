---
title: Phase 4b field coverage audit
author: SWAP modernization team
date: 2026-04-25
---

# Phase 4b field coverage audit

Per-case enumeration of which legacy `.swp` / `.dra` / `.crp` keys must be readable
by the new TOML pipeline. Macropore case (3.macroporeflow) is excluded from
Phase 4b; its schema fields are NOT in scope here.

## Coverage status legend

- ✅ Present in `*_config_t` and reader
- ⚠️ Present in `*_config_t` but reader doesn't populate it  
- ❌ Missing from `*_config_t` entirely (Phase 4b additions needed)
- — Not applicable to this case

## Case 1: hupselbrook

### `swap.swp` keys (General/Simulation/Meteorology)

| Key | Section | Status | Notes |
|---|---|---|---|
| PROJECT | [general] | ✅ | — |
| PATHWORK | [general] | ✅ | — |
| PATHATM | [general] | ✅ | — |
| PATHCROP | [general] | ✅ | — |
| PATHDRAIN | [general] | ✅ | — |
| SWSCRE | [general] | ✅ | Display switch |
| SWERROR | [general] | ✅ | Error switch |
| TSTART | [simulation] | ✅ | — |
| TEND | [simulation] | ✅ | — |
| NPRINTDAY | [simulation] | ✅ | — |
| SWMONTH | [simulation] | ✅ | — |
| PERIOD | [simulation] | ✅ | — |
| SWRES | [simulation] | ✅ | — |
| SWODAT | [simulation] | ✅ | — |
| SWYRVAR | [simulation] | ✅ | — |
| METFIL | [meteorology] | ✅ | — |
| LAT | [meteorology] | ✅ | — |
| SWETR | [meteorology] | ✅ | — |
| ALT | [meteorology] | ✅ | — |
| ALTW | [meteorology] | ✅ | — |
| ANGSTROMA | [meteorology] | ✅ | — |
| ANGSTROMB | [meteorology] | ✅ | — |
| SWDIVIDE | [meteorology] | ✅ | — |
| SWRAIN | [meteorology] | ✅ | Rainfall mode (0-3) |

### `swap.swp` keys (Output files & modelling options)

| Key | Section | Status | Notes |
|---|---|---|---|
| OUTFIL | [output] | ❌ | Generic output filename |
| SWHEADER | [output] | ❌ | Print header flag |
| SWCSV | [output] | ❌ | CSV output switch |
| INLIST_CSV | [output] | ❌ | CSV variable list |
| SWVAP | [output] | ❌ | Soil profile output flag |
| SWBLC | [output] | ❌ | Detailed balance output |
| SWSBA | [output] | ❌ | Cumulative solute balance |
| SWINC | [output] | ❌ | Water balance increments |

### `swap.swp` keys (Soil water section)

| Key | Section | Status | Notes |
|---|---|---|---|
| SWINCO | [soil] | ✅ | Initial condition mode (1-3) |
| GWLI | [soil] | ✅ | Initial groundwater level |
| INIFIL | [soil] | ✅ | Initial file (conditional on SWINCO=3) |
| PONDMX | [soil] | ✅ | Ponding threshold |
| SWPONDMX | [soil] | ✅ | Ponding variation switch |
| RSRO | [soil] | ❌ | Runoff drainage resistance |
| RSROEXP | [soil] | ❌ | Runoff exponent |
| SWRUNON | [soil] | ❌ | Runon switch |
| CFEVAPPOND | [soil] | ❌ | Evaporation coefficient |
| SWCFBS | [soil] | ❌ | Soil factor flag |
| CFBS | [soil] | ❌ | Soil evaporation factor |
| RSOIL | [soil] | ✅ | Soil resistance |
| SWREDU | [soil] | ❌ | Evaporation reduction method |
| COFREDBL | [soil] | ❌ | Black reduction coefficient |
| COFREDBO | [soil] | ❌ | Boesten/Stroosnijder coefficient |
| RSIGNI | [soil] | ❌ | Reset rainfall threshold |
| DTMIN | [soil] | ❌ | Minimum timestep |
| DTMAX | [soil] | ❌ | Maximum timestep |
| GWLCONV | [soil] | ❌ | GWL convergence criterion |
| CRITDEVH1CP | [soil] | ❌ | Relative pressure head tolerance |
| CRITDEVH2CP | [soil] | ❌ | Absolute pressure head tolerance |
| CRITDEVPONDDT | [soil] | ❌ | Ponding water balance error |
| MAXIT | [soil] | ❌ | Max iterations per timestep |
| MAXBACKTR | [soil] | ❌ | Max backtrack cycles |
| SWKMEAN | [soil] | ❌ | K averaging method |
| SWKIMPL | [soil] | ❌ | K update during iteration |

### `swap.swp` keys (Bottom boundary section - SWBOTB=3 path)

| Key | Section | Status | Notes |
|---|---|---|---|
| SWBOTB | [bottom_boundary] | ❌ | Boundary mode selector (1-8) |
| SWBBCFILE | [bottom_boundary] | ❌ | External file flag |
| SHAPE | [bottom_boundary] | ❌ | Shape factor for SWBOTB=3 |
| HDRAIN | [bottom_boundary] | ❌ | Mean drain base |
| RIMLAY | [bottom_boundary] | ❌ | Aquitard resistance |
| SWBOTB3RESVERT | [bottom_boundary] | ❌ | Vertical resistance switch (SWBOTB=3) |
| SWBOTB3IMPL | [bottom_boundary] | ❌ | Implicit solution switch |
| SW3 | [bottom_boundary] | ❌ | Table vs sinus for aquifer head |
| AQAVE | [bottom_boundary] | ❌ | Mean aquifer head |
| AQAMP | [bottom_boundary] | ❌ | Aquifer head amplitude |
| AQTMAX | [bottom_boundary] | ❌ | Time of max aquifer head |
| AQPER | [bottom_boundary] | ❌ | Aquifer head period |
| SW4 | [bottom_boundary] | ❌ | Extra flux flag |

### `swap.swp` keys (Heat flow section - SWHEA=1)

| Key | Section | Status | Notes |
|---|---|---|---|
| SWHEA | [heat] | ❌ | Heat transport switch |
| SWCALT | [heat] | ❌ | Calculation method (1=analytical, 2=numerical) |
| TAMPLI | [heat] | ❌ | Amplitude of annual temp wave |
| TMEAN | [heat] | ❌ | Mean annual temperature |
| TIMREF | [heat] | ❌ | Time of max temperature |
| DDAMP | [heat] | ❌ | Damping depth |
| SWTOPBHEA | [heat] | ❌ | Top BC (1=air temp, 2=measured) |
| TSOILFILE | [heat] | ❌ | Measured soil temp filename |
| SWBOTBHEA | [heat] | ❌ | Bottom BC |

### `swap.swp` keys (Solute section - SWSOLU=1)

| Key | Section | Status | Notes |
|---|---|---|---|
| SWSOLU | [solute] | ❌ | Solute transport switch |
| CPRE | [solute] | ❌ | Precip concentration |
| CDRAIN | [solute] | ❌ | Surface water concentration |
| SWBOTBC | [solute] | ❌ | Seepage concentration mode |
| CSEEP | [solute] | ❌ | Seepage concentration |
| DDIF | [solute] | ❌ | Molecular diffusion |
| TSCF | [solute] | ❌ | Transpiration stream concentration factor |
| SWSP | [solute] | ❌ | Adsorption switch |
| FREXP | [solute] | ❌ | Freundlich exponent |
| CREF | [solute] | ❌ | Reference concentration |
| SWDC | [solute] | ❌ | Decomposition switch |
| GAMPAR | [solute] | ❌ | Temperature reduction factor |
| RTHETA | [solute] | ❌ | Min water content for decomposition |
| BEXP | [solute] | ❌ | Dryness exponent |
| SWBR | [solute] | ❌ | Mixed reservoir switch |
| DAQUIF | [solute] | ❌ | Aquifer thickness |
| POROS | [solute] | ❌ | Aquifer porosity |
| KFSAT | [solute] | ❌ | Aquifer adsorption coefficient |
| DECSAT | [solute] | ❌ | Aquifer decomposition rate |
| CDRAINI | [solute] | ❌ | Initial GW concentration |

### `swap.dra` keys (DRAMET=2 path)

| Key | Section | Status | Notes |
|---|---|---|---|
| DRAMET | [drainage] | ✅ | Drainage method selector |
| SWDIVD | [drainage] | ✅ | Vertical distribution flag |
| SWDISLAY | [drainage] | ✅ | Discharge layer adjustment |
| NRLEVS | [drainage] | ✅ | Number of drainage levels |
| LM2 | [drainage] | ❌ | Drain spacing (Hooghoudt) |
| SHAPE | [drainage] | ❌ | Shape factor (conflicts with bottom_boundary.shape) |
| WETPER | [drainage] | ❌ | Wet perimeter |
| ZBOTDR | [drainage] | ✅ | Drain bottom level (array per level) |
| ENTRES | [drainage] | ✅ | Entry resistance (array) |
| IPOS | [drainage] | ❌ | Position code (1-5) |
| BASEGW | [drainage] | ✅ | Impervious layer depth |
| KHTOP | [drainage] | ❌ | Horizontal K top layer |
| KHBOT | [drainage] | ❌ | Horizontal K bottom layer |
| ZINTF | [drainage] | ❌ | Interface depth (for IPOS 3-5) |
| KVTOP | [drainage] | ❌ | Vertical K top layer |
| KVBOT | [drainage] | ❌ | Vertical K bottom layer |
| GEOFAC | [drainage] | ❌ | Ernst geometry factor |
| COFANI | [drainage] | ❌ | Anisotropy coefficient array |

### Crop files (maizes.crp, potatod.crp, grassd.crp)

| Key | Section | Status | Notes |
|---|---|---|---|
| SWPREP | [crop] | ❌ | Preparation switch |
| SWSOW | [crop] | ❌ | Sowing switch |
| SWGERM | [crop] | ❌ | Germination mode (0-2) |
| SWHARV | [crop] | ❌ | Harvest switch |
| SWGC | [crop] | ❌ | Green cover switch |
| IDEV | [crop] | ❌ | Development stage |
| TBASE | [crop] | ❌ | Base temperature |
| TSUMEA | [crop] | ❌ | Thermal sum emergence to anthesis |
| TSUMAM | [crop] | ❌ | Thermal sum anthesis to maturity |
| TBASEM | [crop] | ❌ | Base temp for germination |
| TEFFMX | [crop] | ❌ | Max eff temp for germination |
| TSUMEMEOPT | [crop] | ❌ | Temp sum for emergence |
| TEMPSOW | [crop] | ❌ | Sowing temperature threshold |
| AGERM | [crop] | ❌ | Additional germination param |
| DVSEND | [crop] | ❌ | DVS at end of season |
| GCTB | [crop] | ❌ | Green cover table |
| LCC | [crop] | ❌ | Leaf area development |
| KDIF | [crop] | ❌ | Extinction coeff diffuse radiation |
| KDIR | [crop] | ❌ | Extinction coeff direct radiation |

## Case 2: grassgrowth

Uses grassd.crp crop file. Config identical to case 1 in terms of top-level sections.
Additional specifics: SWMONTH=0, SWRAIN=2 (rainfall with duration).

Key differences:
- SWIRFIX=0 (no irrigation in this case)
- SWBBCFILE=1 (external BBC file, not covered inline)
- SWSOLU=0 (no solute transport)
- Single crop type (grass) vs rotation

## Case 4: oxygenstress

Uses grassd.crp. Config similar to case 2.
Key specifics:
- SWHEA=0 (no heat, unlike cases 1-2)
- SWSOLU=0 (no solute)
- SWDRA=2 (extended drainage with surface water)

## Case 5: salinitystress

Uses potatod.crp. MAJOR SOLUTE CASE.
Key specifics:
- SWSOLU=1 (active; solute section fully populated)
- SWCROP=1 with potato rotation
- SWIRFIX=1, SWIRGFIL=1 (irrigation from file)
- SWINCO=3 (read from init file)
- SWBOTB=3 (aquifer bottom boundary)
- Bottom section data: SWBOTB3IMPL=0, SW3=1 (sinus aquifer head)
- Solute: CPRE=0.0, CDRAIN=1.051, SWBOTBC=1, CSEEP=15.574
- TSCF=0.5 (root uptake), LDIS array, no adsorption/decomposition

## Case 6: surfacewater

Uses grass.crp. SURFACE WATER EXTENSION case.
Key specifics:
- SWDRA=2 (extended drainage = surface water management)
- SWBBCFILE=1 (external BBC file)
- SWHEA=0 (no heat)
- SWSOLU=0 (no solute)
- Drainage: Uses multi-level drainage (advanced DRA reader needed)

## Aggregate gaps

### `general_config_t` additions needed
None. All fields present.

### `simulation_config_t` additions needed
None. All fields present.

### `meteorology_config_t` additions needed
- SWRAIN 0-3 support (currently missing)
- SWETSINE (daily distribution sine switch)
- NMETDETAIL (for SWMETDETAIL=1)
- RAINFALL table support (if SWRAIN=1; or from separate file)

### `soil_config_t` additions needed (from SoilWaterSection)
- Initial conditions table (ZI, H arrays for SWINCO=1)
- Output switches: OUTFIL, SWHEADER, SWCSV, INLIST_CSV, SWCSV_TZ, SWVAP, SWBLC, SWSBA, SWINC
- Ponding/runoff: RSRO, RSROEXP, SWRUNON, RUFIL
- Evaporation: CFEVAPPOND, SWCFBS, CFBS, SWREDU, COFREDBL, COFREDBO, RSIGNI
- Numerical: DTMIN, DTMAX, GWLCONV, CRITDEVH1CP, CRITDEVH2CP, CRITDEVPONDDT, MAXIT, MAXBACKTR, SWKMEAN, SWKIMPL
- Soil water table (Z, HCOMP, NCOMP arrays for vertical discretization) — partially in place
- Hysteresis: TAU

### `drainage_config_t` additions needed
- IPOS (position code 1-5 for Hooghoudt method)
- LM1, LM2 (drain spacings for different methods)
- KHTOP, KHBOT, KVTOP, KVBOT (soil layer K values)
- ZINTF (interface depth for layered case)
- GEOFAC (Ernst geometry factor)
- COFANI array (anisotropy per soil layer)
- SWTOPDISLAY, ZTOPDISLAY, FTOPDISLAY arrays (discharge layer adjustment)
- COFINTFLB, EXPINTFLB, etc. (interflow params for DRAMET=3)
- SWSRF, SWNRSRF, NMPER, NRSRF, OSSWLM, RSURFDEEP, RSURFSHALLOW, SWQHR, SWSEC (surface water extension params)

### `crop_config_t` additions needed
All current crop rotation infrastructure is present (start, end, file, type). 
Missing are the crop-file-specific keys (parsed from .crp files, not from .swp):
- SWPREP, SWSOW, SWGERM, SWHARV, SWGC
- TBASE, TSUMEA, TSUMAM, TBASEM, TEFFMX, TSUMEMEOPT, TEMPSOW
- AGERM, DVSEND, GCTB, LCC, KDIF, KDIR, IDEV
- (These belong in a crop_variety or crop_parameters type, not in crop_config which is rotation only)

### Bottom boundary section (NEW)
Fully missing; ~70 lines of new type definition needed:
- bottom_boundary_config_t
- Fields: SWBOTB (1-8), SWBBCFILE, and switch-conditional branches:
  - SWBOTB=1: DATE1, GWLEVEL table
  - SWBOTB=2: SW2 (sinus/table), SINAVE, SINAMP, SINMAX; or DATE2, QBOT2 table
  - SWBOTB=3: SWBOTB3RESVERT, SWBOTB3IMPL, SHAPE, HDRAIN, RIMLAY, SW3, (AQAVE/AQAMP/AQTMAX/AQPER vs DATE3/HAQUIF), SW4, (DATE4/QBOT4)
  - SWBOTB=4: SWQHBOT, COFQHA, COFQHB, COFQHC; or HTAB/QTAB table
  - SWBOTB=5: DATE5, HBOT5 table
  - SWBOTB=6,7,8: No extra fields
- Also: BBCFIL (external file reference)

### Heat flow section (NEW)
Fully missing; ~60 lines of new type definition needed:
- heat_config_t
- Fields: SWHEA (0-1)
- If SWHEA=1:
  - SWCALT (1=analytical, 2=numerical)
  - If SWCALT=1: TAMPLI, TMEAN, TIMREF, DDAMP
  - If SWCALT=2: PSAND, PSILT, PCLAY, ORGMAT arrays (per soil layer)
  - SWTOPBHEA (1-2)
  - If SWTOPBHEA=2: TSOILFILE
  - SWBOTBHEA (1-2)
  - If SWBOTBHEA=2: DATET, TBOT table
  - Common: ZH, TSOIL initial condition table

### Solute transport section (NEW)
Fully missing; ~100 lines of new type definition needed:
- solute_config_t
- Fields: SWSOLU (0-1)
- If SWSOLU=1:
  - Part 1: CPRE, CDRAIN, SWBOTBC, CSEEP, (DATEC, CSEEPARR table if SWBOTBC=2), ZC, CML table
  - Part 2: LDIS, KF, DECPOT, FDEPTH arrays (per soil layer)
  - Part 3: DDIF, TSCF
  - Part 4: SWSP, (FREXP, CREF if SWSP=1)
  - Part 5: SWDC, (GAMPAR, RTHETA, BEXP if SWDC=1)
  - Part 6: SWBR, (DAQUIF, POROS, KFSAT, DECSAT, CDRAINI if SWBR=1)

### Irrigation section (NEW, needed for case 5)
Conditional structure:
- irrigation_config_t
- Fields: SWIRFIX (0-1), SWIRGFIL (0-1)
- If SWIRFIX=1 and SWIRGFIL=0: IRDATE, IRDEPTH, IRCONC, IRTYPE table
- If SWIRFIX=1 and SWIRGFIL=1: IRGFIL (external file)

### Output/Balance section (partial gap)
Some fields already in general_config (project, paths). 
Missing: OUTFIL, SWHEADER, SWCSV, INLIST_CSV, SWCSV_TZ, INLIST_CSV_TZ, SWVAP, SWBLC, SWSBA, SWINC.
These may belong in a dedicated output_config_t (10-15 fields).

## Recommendations for Task 7-9

### Priority 1 (Foundation, blocks 4 of 5 cases)
1. **bottom_boundary_config_t + reader** (~70 LoC)
   - Cases: 1, 2, 4, 5, 6 (all non-macropore)
   - Blocking: Bottom BC data reading entirely absent
   - Complexity: High (SWBOTB=1-8 branching logic)
   
2. **heat_config_t + reader** (~60 LoC)
   - Cases: 1, 2, 5, 6 (SWHEA=1)
   - Blocking: Heat transport simulation setup
   - Complexity: Medium (SWCALT=1-2 branching, soil texture arrays)

### Priority 2 (Extensions, blocks 1-2 cases)
3. **solute_config_t + reader** (~100 LoC)
   - Cases: 1, 5 (SWSOLU=1)
   - Blocking: Solute transport model initialization
   - Complexity: High (6-part structure, conditional fields, tables)

4. **drainage_config_t extensions** (~40 LoC additions)
   - Cases: 1, 6 (extended params)
   - Blocking: IPOS-dependent Hooghoudt method, surface water extension
   - Complexity: Medium (IPOS=1-5 branching, DRAMET=3 advanced features)

### Priority 3 (Enhancements)
5. **irrigation_config_t + reader** (~25 LoC)
   - Cases: 5 only
   - Blocking: Fixed irrigation applications
   - Complexity: Low (simple table or file reference)

6. **output_config_t** (~15 LoC)
   - All cases (but optional, may defer)
   - Blocking: CSV output specification, balance file flags
   - Complexity: Very low (mostly boolean switches + CSV list)

7. **meteorology_config_t extensions** (~10 LoC)
   - Cases: 1, 2, 5 (SWRAIN, rainfall table)
   - Blocking: Advanced rainfall intensity/duration handling
   - Complexity: Low (SWRAIN=0-3 branching, optional table)

8. **soil_config_t numerical solver params** (~20 LoC)
   - All cases
   - Blocking: Richards solver tuning
   - Complexity: Very low (direct assignment, no branching)

### Estimated effort (in LoC, combined type + reader + tests)
- Task 7 (bottom_boundary + heat + drainage ext): ~200-250 LoC
- Task 8 (solute + irrigation): ~150-200 LoC
- Task 9 (output + meteorology + misc): ~80-100 LoC

### Dependency order
1. Define all new config types first (bottom_boundary, heat, solute, irrigation, output)
2. Update main_config_t to include new sections
3. Write readers for each in readers/ module
4. Integrate with TOML parser
5. Test on cases in order: 1 (hupselbrook, baseline), 2 (grassgrowth), 4 (oxygen, no heat), 5 (salinity, all features), 6 (surface water, extended drainage)

