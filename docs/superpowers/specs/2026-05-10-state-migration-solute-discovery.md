# Subsystem Migration Discovery: Solute

**Date:** 2026-05-10
**Status:** discovery (read-only inventory)
**Migration #:** 3 of N (surface-water + drainage already shipped)
**Branch:** `refactor/solute-state`
**Predecessor playbooks:**
- `docs/superpowers/specs/2026-05-08-state-migration-surfacewater-discovery.md`
- `docs/superpowers/specs/2026-05-09-state-migration-drainage-discovery.md` (with the Phase 1+2 lessons-learned addendum on Section 3.5 categorization)

## 1. Source files

| File | LoC | Role |
|------|-----|------|
| `src/solute/solute.f90` | 548 | Main compute: `solute(task, state)` (task 1 = init, task 2 = rate) and `AgeTracer(task, state)` (task 1 = init, task 2 = rate). Both already accept `state` (plumbed in SS-SWST Phase 2 Task 5). Both use bare `use Variables` (no `only:` guard). |
| `src/config/solute_config.f90` | 132 | Typed config struct `solute_config_t`. Fields: `swsolu`, `swbotbc`, `cdrain`, `cseep`, `tscf`, `ldis`/`ldis_array`, `rtheta`, `bexp`, `ecmax`, `ecslop`, `swsoltyp`, `swdc`, `pertabsolu`. `validate` + `finalize` procedures. |
| `src/io/toml/read_solute_toml.f90` | 187 | TOML reader: populates `solute_config_t` from `[solute]` section. Handles scalar vs. per-layer `ldis` array; decodes `pertabsolu` as array-of-arrays. Returns silently when section is absent. |
| `tests/unit/config/test_solute_config.pf` | 93 | pFUnit tests for `solute_config_validate`: sentinel (`swsolu=0`), happy path (case 5 values), and invalid-enum rejection. |
| `tests/unit/io/toml/test_read_solute_toml.pf` | 136 | pFUnit tests for `read_solute_toml`: absent section, minimal `swsolu=1`, ldis as scalar, ldis as array, malformed `pertabsolu`. |
| `tests/unit/io/toml/fixtures/solute_swsolu0.toml` | — | Fixture: empty/absent sentinel case. |
| `tests/unit/io/toml/fixtures/solute_swsolu1_minimal.toml` | — | Fixture: minimal swsolu=1 case. |
| `tests/unit/io/toml/fixtures/solute_swsolu1.toml` | — | Fixture: full swsolu=1 with pertabsolu. |
| `tests/unit/io/toml/fixtures/solute_malformed_pertabsolu.toml` | — | Fixture: bad pertabsolu for error-path test. |
| `tests/unit/io/toml/fixtures/swap_with_solute.toml` | — | Fixture: full swap TOML with [solute] section. |
| `tests/swap-cases/toml/5.salinitystress/` | — | End-to-end regression case with `swsolu=1` (TOML path); covers cml profile output and CSV concentration columns. |

## 2. Owned globals (this subsystem writes)

Both `solute` and `AgeTracer` use a bare `use Variables` clause (no `only:` list) so the owned set is determined by grepping for `<symbol>\s*=` assignments inside `src/solute/solute.f90`.

| Variable | Type | `variables.f90` line | Description | Write sites (file:line) |
|----------|------|----------------------|-------------|--------------------------|
| `cml(macp)` | `real(8)` | 1050 | Solute concentration (M/L3 water) in mobile region | `solute.f90:54,219,222,224,228` (init profile, iteration solve); `solute.f90:539` (AgeTracer overwrites with age data — see Hazard #1) |
| `cmsy(macp)` | `real(8)` | 1051 | Dissolved+adsorbed solute concentration (M/L3 soil volume) | `solute.f90:65,213,218` |
| `samini` | `real(8)` | 1085 | Total solute (M/L2) in profile at start of balance period | `solute.f90:59,66,70,94,370,373` |
| `sampro` | `real(8)` | 1086 | Total solute (M/L2) in soil column | `solute.f90:70,281,283,285` |
| `cdrain` | `real(8)` | 1048 | Mean solute conc in aquifer/drainage system | `solute.f90:244,247,250` (AgeTracer init: line 365 sets `Agedrain=cdrain`, reads cdrain) |
| `cseep` | `real(8)` | 1055 | Mean solute conc in upward seepage at bottom | `solute.f90:107,250` |
| `csurf` | `real(8)` | 1057 | Total solutes (M/L2) in ponding layer | `solute.f90:86,93,141,145` |
| `cpond` | `real(8)` | 1052 | Mean solute conc (M/L3) in ponding layer | `solute.f90:143,146,148` |
| `isqbot` | `real(8)` | 1072 | Instantaneous solute flux at profile bottom (M/L2/T) | `solute.f90:97,273,275,392` |
| `isqtop` | `real(8)` | 1073 | Instantaneous solute flux through soil surface (M/L2/T) | `solute.f90:98,146,393,428` |
| `sqdra` | `real(8)` | 1090 | Cumulative solute to drainage canals (M/L2) | `solute.f90:89,209,481` |
| `imsqdra` | `real(8)` | 1091 | Intermediate solute to drainage (M/L2) | `solute.f90:81,210` |
| `sqprec` | `real(8)` | 1094 | Cumulative solute in precipitation (M/L2) | `solute.f90:86,289` |
| `imsqprec` | `real(8)` | 1095 | Intermediate solute in precipitation | `solute.f90:78,290` |
| `sqirrig` | `real(8)` | 1092 | Cumulative solute in irrigation water (M/L2) | `solute.f90:87,291` |
| `imsqirrig` | `real(8)` | 1093 | Intermediate solute in irrigation | `solute.f90:79,292` |
| `sqbot` | `real(8)` | 1088 | Cumulative solute through profile bottom (M/L2) | `solute.f90:88,260,261,263,264` |
| `imsqbot` | `real(8)` | 1089 | Intermediate solute through bottom | `solute.f90:80,261,264` |
| `sqsur` | `real(8)` | 1097 | Cumulative solute to surface water (M/L2) | `solute.f90:90,255` |
| `dectot` | `real(8)` | 1062 | Cumulative solute decomposition (M/L2) | `solute.f90:91,186` |
| `imdectot` | `real(8)` | 1063 | Intermediate decomposition | `solute.f90:82,187` |
| `rottot` | `real(8)` | 1078 | Cumulative solute extracted by roots (M/L2) | `solute.f90:92,191,461` (both `solute` and `AgeTracer` accumulate into this) |
| `imrottot` | `real(8)` | 1079 | Intermediate root extraction | `solute.f90:83,192` |
| `solbal` | `real(8)` | 1087 | Cumulative solute balance (M/L2) | `solute.f90:295` |
| `dtsolu` | `real(8)` | 1064 | Solute sub-timestep (T) | `solute.f90:111,126,136,137,401,408,417,418` |
| `ArMpSs` | `real(8)` | 1111 | Area fraction of macropores at soil surface (-) | `solute.f90:102,103,397,398` (also written by `soilhydraulics.f90:101,102` and `boundtop.f90:129,130` — see Hazard #3) |
| `Ageirr` | `real(8)` | 1104 | Age of irrigation water (d) — AgeTracer | `solute.f90:364` |
| `Agedrain` | `real(8)` | 1105 | Age of drainage water (d) — AgeTracer | `solute.f90:365` |
| `Agepre` | `real(8)` | 1106 | Age of precipitation (d) — AgeTracer | `solute.f90:363` |
| `Agepond` | `real(8)` | 1107 | Age of ponding water (d) — AgeTracer | `solute.f90:366,425,432,535` |
| `Agepondm1` | `real(8)` | 1108 | Age of ponding water previous timestep (d) — AgeTracer | `solute.f90:367,535` |
| `icAgetopdwn` | `real(8)` | 1110 | Incremental age entering top compartment downward (d) | `solute.f90:383,429` |
| `icAgetopupw` | `real(8)` | 1109 | Incremental age leaving top compartment upward (d) | `solute.f90:384,498` |
| `icAgeSur` | `real(8)` | 1071 | Incremental age leaving by surface runoff (d) | `solute.f90:382,435` |
| `icAgeRot` | `real(8)` | 1070 | Incremental age by root uptake (d) | `solute.f90:501` |
| `icAgeBot` | `real(8)` | 1068 | Incremental age leaving bottom (d) | `solute.f90:389,508` |
| `icAgeDra(madr)` | `real(8)` | 1069 | Incremental age per drainage level (d) | `solute.f90:387,504` |
| `AgeGwl1m` | `real(8)` | 1046 | Age of groundwater in upper 1 m of saturated zone (d) | `solute.f90:530,531` |

**Note on AgeTracer local arrays:** `Ageml(macp)` and `Agemsy(macp)` are declared as subroutine-local automatic arrays at `solute.f90:329,330`. They are NOT module-level globals. Their values are NOT preserved between timesteps — the comment at line 335–337 notes that the age-tracer state variables that need persistence between timesteps have been moved to `variables.f90` (`Ageirr`, `Agedrain`, `Agepre`, `Agepond`, `Agepondm1`, `icAgetopdwn`, `icAgetopupw`, and `ArMpSs`).

## 3. Borrowed globals (this subsystem reads, owned elsewhere)

Variables read in `src/solute/solute.f90` but not written by it. `state%surfacewater%*` and `state%drainage%*` reads are excluded (already typed state).

| Variable | Best-guess owner subsystem | Read sites (file:line) |
|----------|----------------------------|--------------------------|
| `numnod` | core/grid | `solute.f90:53,61,112,154,198,281,354,371,402,439,500,515,539` |
| `dt` | core/timecontrol | `solute.f90:111,136,289,290,291,292` |
| `dtmin` | core/timecontrol | `solute.f90:137,418` |
| `theta(macp)` | soilhydraulics | `solute.f90:65,67,113,183,185,372,403,489` |
| `thetm1(macp)` | soilhydraulics | `solute.f90:484` |
| `q(macp+1)` | soilhydraulics | `solute.f90:116,121,164,166,169,259,263,272,405,445,449,452,455,498,508` |
| `qbot` | soilhydraulics | `solute.f90:259,260,263,264` |
| `qtop` | soilhydraulics/boundtop | `solute.f90:142,144,146,424,426,428` |
| `bdens(maho)` | soil physics (config) | `solute.f90:62,64` |
| `kf(maho)` | soil physics (config) | `solute.f90:62` |
| `kfsat` | soil physics (config) | `solute.f90:64` |
| `poros` | soil physics (config) | `solute.f90:64` |
| `cref` | solute config (legacy global, not yet in typed config) | `solute.f90:63,65,185,227` |
| `frexp` | solute config (legacy global) | `solute.f90:65,66,185,220,227` |
| `rtheta` | solute config (set by adapter) | `solute.f90:183` |
| `bexp` | solute config (set by adapter) | `solute.f90:183` |
| `tscf` | solute config (set by adapter) | `solute.f90:190,191` |
| `thetsl(maho)` | soil physics (config) | `solute.f90:67,404,446` |
| `tsoil(macp)` | heat subsystem | `solute.f90:175,176` |
| `fltemperature` | heat subsystem | `solute.f90:174` |
| `gampar` | solute config (legacy global, not yet in typed config) | `solute.f90:176,178` |
| `decpot(maho)` | solute config (legacy global) | `solute.f90:68` |
| `fdepth(maho)` | solute config (legacy global) | `solute.f90:68` |
| `ddif` | solute config (legacy global) | `solute.f90:67,404` |
| `ldis(maho)` | solute config (set by adapter) | `solute.f90:117,121,405,447` |
| `layer(macp)` | core/grid | `solute.f90:62,63,64,67,68,117,121,404,405` |
| `z(macp)` | core/grid | `solute.f90:54,355` |
| `dz(macp)` | core/grid | `solute.f90:66,124,186,190,200,202,208,213,283,407,480,488` |
| `inpola(macp)` | core/grid | `solute.f90:113,158,355,403,443` |
| `inpolb(macp)` | core/grid | `solute.f90:113,158,355,403,443` |
| `disnod(macp+1)` | core/grid | `solute.f90:164,450` |
| `nraidt` | meteo | `solute.f90:141,289,290` |
| `nird` | irrigation | `solute.f90:141,291,292` |
| `cirr` | irrigation (set from legacy reader) | `solute.f90:141,291,364` |
| `cpre` | unowned (only initialized to 0; see Hazard #5) | `solute.f90:141,289,290,363,366` |
| `pond` | surface/boundary | `solute.f90:143` |
| `pondm1` | surface/boundary | `solute.f90:423` |
| `t1900` | core/timecontrol | `solute.f90:107` |
| `nrlevs` | drainage | `solute.f90:198,386,386,503` |
| `runots` | meteo/surface | `solute.f90:435` |
| `nodgwl` | soilhydraulics | `solute.f90:515` |
| `gwl` | soilhydraulics | `solute.f90:520,522,531` |
| `zbotcp(macp)` | core/grid | `solute.f90:516` |
| `ztopcp(macp)` | core/grid | `solute.f90:517,527` |
| `thetas(macp)` | soilhydraulics | `solute.f90:524` |
| `qrot(macp)` | crop/rootextraction | `solute.f90:190,191,460,461,501` |
| `FlMacropore` | macropore | `solute.f90:103,398` |
| `ArMpTp` | macropore | `solute.f90:103,398` |
| `Z_Tp` | macropore | `solute.f90:103,398` |
| `flzerointr` | core/timecontrol | `solute.f90:77` |
| `flzerocumu` | core/timecontrol | `solute.f90:85,381` |
| `swbotbc` | solute config (set by adapter) | `solute.f90:106` |
| `cseeptab(mabbc*2)` | solute config (legacy global, not in typed config) | `solute.f90:107` |
| `swbr` | solute config (legacy global, not in typed config) | `solute.f90:242,254` |
| `nconc` | solute config (set by adapter from CSV) | `solute.f90:49,52` |
| `zc(macp)` | solute config (set by adapter from CSV) | `solute.f90:51` |
| `daquif` | solute config (legacy global, not in typed config) | `solute.f90:245,248` |
| `decsat` | solute config (legacy global, not in typed config) | `solute.f90:245,248` |
| `swinco` | core/config | `solute.f90:48,349` |
| `swsolu` | solute config (set by adapter) | used via `flSolute` guard in `swap.f90:200,325` |

## 3.5 External readers of owned globals

The two grep templates from the drainage Section 3.5 addendum were applied for each owned global. Results are classified by the 5-category framework.

### Grep commands used

```bash
# Categories 1+2 (use variables clauses):
grep -rEn "use variables.*\b<owned-var>\b" src/ --include="*.f90" \
  | grep -v "src/core/variables.f90" \
  | grep -v "src/state/" \
  | grep -v "src/solute/" \
  | grep -v "src/config/solute_config" \
  | grep -v "src/io/toml/read_solute_toml"

# Categories 3+4+5 (raw symbol references):
grep -rEn "\b<owned-var>\b" src/ --include="*.f90" \
  | grep -v "src/core/variables.f90" \
  | grep -v "src/core/initialize.f90" \
  | grep -v "src/state/" \
  | grep -v "state%" \
  | grep -v "config%" \
  | grep -v "src/solute/" \
  | grep -v "src/config/solute_config" \
  | grep -v "src/io/toml/read_solute_toml"
```

### Results per owned global

**`cml(macp)`** — most heavily read global in the codebase:

- **Output (Cat 1):** `src/io/swapoutput.f90:577` (`outvap` — `use variables, only: ...,cml,...`; writes each node's cml to `.vap`); `src/io/swapoutput.f90:1041` (`outend` — writes final cml profile to `.end` file); `src/io/swapoutput.f90:1696` (`outage` — uses cml as AgeTracer age profile in `.ageProfile.csv` output)
- **Output (Cat 1, inline CSV):** `src/io/swapoutput.f90:4636` (module-level `use variables, only: ...,cml,...`; `4739` writes `cml(j)` to per-node time-series CSV)
- **Output (Cat 1, swap_csv_output):** `src/io/swap_csv_output.f90:10` (module-level `use variables, only: ...,cml,...`); `481` (`lp_C%vals = cml(nodes)` — per-node concentration column); `952,1088` (CSV time-series writes)
- **Compute (Cat 2):** `src/crop/rootextraction.f90:176,177` (reads `cml(node)` for salinity stress check: `if (cml(node) .gt. saltmax)`); `src/crop/rootextraction.f90:838` (reads `cml(node)` for osmotic head: `hosm = salthead * cml(node)`)
- **Compute (Cat 2):** `src/crop/irrigation.f90:269` (reads `cml(nodsen)` to check if irrigation concentration exceeds threshold)
- **Init-seed (Cat 4):** `src/io/toml/config_to_variables.f90:628,629,635` (populates `cml(k)` from initial concentration CSV file at config time)

**`cmsy(macp)`:**

- **Output (Cat 1):** `src/io/swapoutput.f90:578` (`outvap` reads `cmsy`); `src/io/swapoutput.f90:4017` (CSV module use clause); `src/io/swapoutput.f90:4291` (writes `cmsy(Nodes_CONCADS(j))` per node to CSV)
- **Output (Cat 1, swap_csv_output):** `src/io/swap_csv_output.f90:10` (module use); `src/io/swap_csv_output.f90:952,1089` (time-series CSV writes)

**`sampro`:**

- **Output (Cat 1):** `src/io/swapoutput.f90:897,946,948,951` (`outbal` — reads for .bal water+solute balance summary); `src/io/swapoutput.f90:1581` (`outsba` — reads for `.sba` cumulative solute balance output); `src/io/swapoutput.f90:1626` (writes to `.sba`)
- **Output (Cat 1, swap_csv_output):** `src/io/swap_csv_output.f90:8` (module use); `src/io/swap_csv_output.f90:296` (`set_values` maps `sampro` to SAMPRO column)

**`samini`:**

- **Output (Cat 1):** `src/io/swapoutput.f90:897,948,951` (`outbal` reads `samini` for balance)

**`sqdra`, `sqprec`, `sqirrig`, `sqbot`, `dectot`, `rottot`, `solbal`:**

- **Output (Cat 1):** `src/io/swapoutput.f90:897,984,985,986` (`outbal` reads all for `.bal` balance summary); `src/io/swapoutput.f90:1581,1624,1625,1626` (`outsba` reads `sampro,sqbot,sqdra,solbal,dectot,rottot,sqprec,sqirrig` for `.sba` output)

**`imsqprec`, `imsqirrig`, `imsqbot`, `imsqdra`, `imdectot`, `imrottot`:**

- **Output (Cat 1):** `src/io/swapoutput.f90:4018,4019,4236–4241` (time-series CSV writes of intermediate fluxes for each output interval)
- **Output (Cat 1, swap_csv_output):** `src/io/swap_csv_output.f90:8` (module use); `src/io/swap_csv_output.f90:290–295` (`set_values` maps each to named CSV columns SQPREC, SQIRRIG, SQBOT, SQDRA, DECTOT, ROTTOT)

**`isqtop`, `isqbot`:**

- **Output (Cat 1):** `src/io/swapoutput.f90:577` (`outvap` use clause includes `isqtop,isqbot`); `swapoutput.f90:619,621,629,640,642,650,679,681,689,698,700,708` (`outvap` uses both as local `sflux` for node-level flux output in `.vap`)

**`icAgeBot`, `icAgeDra`, `icAgeRot`, `icAgeSur`, `AgeGwl1m`:**

- **Output (Cat 1):** `src/io/swapoutput.f90:1696,1697` (`outage` — `use variables, only: ...,AgeGwl1m,icAgeBot,icAgeDra,icAgeRot,icAgeSur`); `swapoutput.f90:1749,1769,1770,1771,1772` (writes to `.ageEffluent.csv`)

**`ArMpSs`:**

- **Compute (Cat 2):** `src/soil/soilhydraulics.f90:101,102` (resets to 0 then computes from `ArMpTp` — ALSO WRITES this variable); `soilhydraulics.f90:109,651` (reads after own write for boundary flux)
- **Compute (Cat 2):** `src/boundary/boundtop.f90:129,130` (resets to 0 then computes from `ArMpTp` — ALSO WRITES); `boundtop.f90:131,174` (reads after own write)

Note: `ArMpSs` is a working-buffer variable computed fresh each timestep by three subsystems (soilhydraulics, boundtop, and solute) from the same macropore inputs. It is not truly "owned" by solute — it is a shared re-derivation (see Hazard #3).

### Section 3.5 tally

| Owned global | Output readers (Cat 1) | Compute readers (Cat 2) | Working-buffer (Cat 3) | Init-seed (Cat 4) | Call-site arg (Cat 5) |
|---|---|---|---|---|---|
| `cml` | `swapoutput.f90` (outvap, outend, outage, CSV-tz), `swap_csv_output.f90` (set_values, per-node) | `rootextraction.f90` (salt stress), `irrigation.f90` (conc threshold) | — | `config_to_variables.f90` (cml_file CSV seed) | — |
| `cmsy` | `swapoutput.f90` (outvap, CSV), `swap_csv_output.f90` | — | — | — | — |
| `sampro` | `swapoutput.f90` (outbal, outsba), `swap_csv_output.f90` (set_values) | — | — | — | — |
| `samini` | `swapoutput.f90` (outbal) | — | — | — | — |
| `sqdra`, `sqprec`, `sqirrig`, `sqbot`, `dectot`, `rottot`, `solbal` | `swapoutput.f90` (outbal, outsba) | — | — | — | — |
| `imsqprec`…`imrottot` (6 vars) | `swapoutput.f90` (CSV), `swap_csv_output.f90` (set_values) | — | — | — | — |
| `isqtop`, `isqbot` | `swapoutput.f90` (outvap) | — | — | — | — |
| `icAgeBot/Dra/Rot/Sur`, `AgeGwl1m` | `swapoutput.f90` (outage) | — | — | — | — |
| `ArMpSs` | — | `soilhydraulics.f90`, `boundtop.f90` (also write it) | working-buffer: recomputed fresh each call | — | — |
| `cdrain`, `cseep`, `csurf`, `cpond`, `sqsur`, `imdectot`, `imrottot`, `dtsolu`, `Ageirr`, `Agedrain`, `Agepre`, `Agepond`, `Agepondm1`, `icAgetopdwn/upw` | — | — | — | — | — |

**External reader files (union across all 5 categories):** 5 files:
1. `src/io/swapoutput.f90` (Cat 1: output — outvap, outend, outage, outbal, outsba, CSV inline)
2. `src/io/swap_csv_output.f90` (Cat 1: output — set_values, per-node time-series)
3. `src/crop/rootextraction.f90` (Cat 2: compute — salt/osmotic stress)
4. `src/crop/irrigation.f90` (Cat 2: compute — concentration threshold check)
5. `src/io/toml/config_to_variables.f90` (Cat 4: init-seed — cml CSV seed)

Additionally: `src/soil/soilhydraulics.f90` and `src/boundary/boundtop.f90` both write AND read `ArMpSs`, making them working-buffer co-writers (Cat 3).

## 4. Entry points (subroutines called from outside this subsystem)

| Subroutine | Called from (file:line) | Purpose / lifecycle stage |
|-------------|--------------------------|----------------------------|
| `solute(task=1, state)` | `src/core/swap.f90:200` | Initialize: set initial cml profile, compute cmsy and samini |
| `solute(task=2, state)` | `src/core/swap.f90:325` | Rate calculation: sub-timestepped transport, decomposition, drainage losses, balance update |
| `AgeTracer(task=1, state)` | `src/core/swap.f90:203` | Initialize: set age boundary conditions from solute BC globals |
| `AgeTracer(task=2, state)` | `src/core/swap.f90:328` | Rate calculation: age transport, writes age values into `cml(i)` at end |

Both are guarded at call sites: `solute` by `if (flSolute)` and `AgeTracer` by `if (flAgeTracer)`. `flSolute` is set in `src/core/timecontrol.f90:120–121` (`flSolute = (swsolu == 1)`). `flAgeTracer` is currently **only set to `.false.`** in `initialize.f90:584` and **never set to `.true.`** anywhere in the codebase (see Hazard #4).

Additionally, output entry points called from `swap.f90`:
| Subroutine | Called from | Purpose |
|-------------|-------------|---------|
| `SoluteOutput(task=1,2,3)` | `swap.f90:216,377,426` | Open/write/close `.sba` file (delegates to `outsba`) |
| `AgeTracerOutput(task=1,2,3, state)` | `swap.f90:217,378,427` | Open/write/close `.ageProfile.csv`, `.ageEffluent.csv` (delegates to `outage`) |

## 5. Internal call graph

```
solute_mod
├── subroutine solute(task, state)
│   ├── task=1 (init)
│   │   ├── afgen(tab, macp*2, abs(z(i)))     [array_utils, pure function]
│   │   └── [no further calls; fatalerr_collected on bad task]
│   └── task=2 (rate)
│       ├── afgen(cseeptab, mabbc*2, t1900+dt)  [array_utils, swbotbc=2 only]
│       ├── [associate qdra => state%drainage%qdra]
│       ├── [associate qdrtot => state%surfacewater%qdrtot]
│       └── fatalerr_collected('Solute', ...)   [error_mod, default case]
│
└── subroutine AgeTracer(task, state)
    ├── task=1 (init)
    │   └── afgen(tab, macp*2, abs(z(i)))     [array_utils]
    └── task=2 (rate)
        ├── [associate qdra => state%drainage%qdra]
        └── fatalerr_collected('AgeTracer', ...) [error_mod, default case]
```

Neither `solute` nor `AgeTracer` calls the other. Both are flat loops with no private helper subroutines.

## 6. Output coupling

| Output writer (file:subroutine) | Reads which solute fields | Destination file |
|-----------------------------------|--------------------------|------------------|
| `swapoutput.f90:outvap` | `cml(node)`, `cmsy(node)`, `isqtop`, `isqbot` | `*.vap` (soil profile time-series) |
| `swapoutput.f90:outend` | `cml(i)` (conditional on `flSolute .or. flAgeTracer`) | `*.end` (final state snapshot) |
| `swapoutput.f90:outbal` | `samini`, `sampro`, `sqprec`, `sqirrig`, `sqbot`, `dectot`, `rottot`, `sqdra` | `*.bal` (water + solute balance overview) |
| `swapoutput.f90:outsba` | `sampro`, `sqbot`, `sqdra`, `solbal`, `dectot`, `rottot`, `sqprec`, `sqirrig` | `*.sba` (cumulative solute balance) |
| `swapoutput.f90:outage` | `cml(node)` (age profile), `AgeGwl1m`, `icAgeBot`, `icAgeRot`, `icAgeSur`, `icAgeDra(:)` | `*.ageProfile.csv`, `*.ageEffluent.csv`, `*.ageEffluentqDrain.csv` |
| `swapoutput.f90` (inline CSV block ~4017–4291) | `imsqprec`, `imsqirrig`, `imsqbot`, `imsqdra`, `imdectot`, `imrottot`, `sampro`, `cml(nodes)`, `cmsy(nodes)` | `swap.csv` (per-timestep CSV) |
| `swap_csv_output.f90:set_values` | `imsqprec`, `imsqirrig`, `imsqbot`, `imsqdra`, `imdectot`, `imrottot`, `sampro`, `solbal`, `cml(nodes)`, `cmsy(nodes)` | `swap.csv` columns SQPREC…SAMPRO, C[node], CA[node] |

## 7. Config inputs

The TOML adapter (`src/io/toml/config_to_variables.f90`, Solute section ~999–1020) writes the following legacy globals from `config%solute%*`:

| TOML key | Typed-config field | Legacy global written | Used at (file:line) |
|---|---|---|---|
| `solute.swsolu` | `config%solute%swsolu` | `swsolu` (→ `flSolute` in timecontrol) | `config_to_variables.f90:1002`; `timecontrol.f90:120–121` |
| `solute.swbotbc` | `config%solute%swbotbc` | `swbotbc` | `config_to_variables.f90:1003`; `solute.f90:106` |
| `solute.cdrain` | `config%solute%cdrain` | `cdrain` | `config_to_variables.f90:1004`; `solute.f90:200–248` |
| `solute.cseep` | `config%solute%cseep` | `cseep` | `config_to_variables.f90:1005`; `solute.f90:107,167,250,260,263,273` |
| `solute.tscf` | `config%solute%tscf` | `tscf` | `config_to_variables.f90:1006`; `solute.f90:190,191` |
| `solute.rtheta` | `config%solute%rtheta` | `rtheta` | `config_to_variables.f90:1007`; `solute.f90:183` |
| `solute.bexp` | `config%solute%bexp` | `bexp` | `config_to_variables.f90:1008`; `solute.f90:183` |
| `solute.ldis` (scalar) or `solute.ldis` (array) | `config%solute%ldis` / `ldis_array` | `ldis(maho)` (broadcast or element-wise) | `config_to_variables.f90:1014–1020`; `solute.f90:117,121,405,447` |

**Solute physics globals NOT yet in the TOML adapter** (read from legacy `.swp` input only; absent from `[solute]` section and `solute_config_t`):

| Legacy global | `variables.f90` line | Role in solute.f90 | Comment |
|---|---|---|---|
| `cpre` | 1053 | Solute conc in precipitation | Only initialized to 0; never set via TOML or any source beyond initialize.f90 — see Hazard #5 |
| `cref` | 1054 | Freundlich reference conc | `solute.f90:63,65,185,227` |
| `kf(maho)` | 1074 | Freundlich adsorption coeff per layer | `solute.f90:62` |
| `kfsat` | 1075 | Linear adsorption coeff in aquifer | `solute.f90:64` |
| `poros` | 1077 | Aquifer porosity | `solute.f90:64` |
| `daquif` | 1058 | Aquifer thickness | `solute.f90:245,248` |
| `decsat` | 1061 | Decomposition rate in aquifer | `solute.f90:245,248` |
| `gampar` | 1067 | Temp reduction factor for decomp | `solute.f90:176,178` |
| `decpot(maho)` | 1060 | Potential decomposition rate per layer | `solute.f90:68` |
| `fdepth(maho)` | 1065 | Depth reduction factor for decomp per layer | `solute.f90:68` |
| `ddif` | 1059 | Molecular diffusion coefficient | `solute.f90:67,404` |
| `cseeptab(mabbc*2)` | 1056 | Seepage conc as function of time (tabulated) | `solute.f90:107` — only used when `swbotbc=2` |
| `swbr` | 1042 | Aquifer breakthrough switch (0/1) | `solute.f90:242,254` |
| `frexp` | 1066 | Freundlich exponent | `solute.f90:65,185,220,221,227` |
| `cirr` | 1049 | Solute conc in irrigation water | `solute.f90:141,291,364` — set by `irrigation.f90` from per-crop schedule |
| `nconc` | 1040 | Number of initial concentration values | `solute.f90:49,52` — set by adapter from CSV |
| `zc(macp)` | 1099 | Depths for initial conc profile | `solute.f90:51` — set by adapter from CSV |

The adapter section also reads `config%soil%initial%cml_file` to seed `cml(k)` and `zc(k)` via CSV (`config_to_variables.f90:619–637`).

## 8. Open questions / cross-subsystem hazards

1. **Hazard: AgeTracer overwrites `cml` at end of task=2** (`solute.f90:537–540`). At the end of every `AgeTracer(2, state)` call, the local `Ageml(i)` array (which holds water age in days, not a solute concentration) is copied into `cml(i)`. This means that after `AgeTracer` runs, `cml` no longer holds the solute concentration — it holds groundwater age. The `outage` subroutine relies on this deliberately (it writes `cml(node)` to `.ageProfile.csv` as "age"). However, it is an architectural smell: `cml` is dual-purposed as solute concentration (via `solute`) and age storage (via `AgeTracer`). A solute state-type must accommodate this dual use: either separate `cml_age(macp)` or a runtime discriminator flag. Since `flSolute` and `flAgeTracer` are mutually guarded, in practice only one ever runs, but the current implementation assumes this remains true.

2. **Hazard: `rottot` and `sqdra` accumulated by both `solute` and `AgeTracer`** (`solute.f90:191,461` for rottot; `209,481` for sqdra). If both were enabled (they cannot currently be simultaneously), these cumulative balances would be double-counted. The design for solute state-type migration must ensure the two routines have clearly separated output fields.

3. **Hazard: `ArMpSs` is co-written by three subsystems** — `solute.f90:102–103`, `soilhydraulics.f90:101–102`, `boundtop.f90:129–130`. All three reset it to 0 and then re-derive it from `ArMpTp` and `FlMacropore` via the same identical pattern. It is a working-buffer variable, not true subsystem state. It should NOT be migrated into `solute_state_t`; instead it should remain a module-level scratch variable or become a local in each routine that computes it.

4. **Hazard: `flAgeTracer` is never set to `true`**. `initialize.f90:584` sets it to `.false.` and no other code ever sets it to `.true.`. The AgeTracer feature is effectively dead in the current codebase. The TOML adapter has no `swAgeTracer` key and `timecontrol.f90` has no `if (swAgeTracer)` block. Before solute state migration, a decision is needed: retire AgeTracer or explicitly re-enable it. The migration must account for the AgeTracer-owned globals (`Ageirr`, `Agedrain`, `Agepre`, `Agepond`, `Agepondm1`, `icAgetopdwn`, `icAgetopupw`, `icAgeBot`, `icAgeDra`, `icAgeRot`, `icAgeSur`, `AgeGwl1m`) regardless of whether the feature is active.

5. **Hazard: `cpre` is stuck at 0**. The solute concentration in precipitation (`cpre`, `variables.f90:1053`) is initialized to 0 in `initialize.f90:537` and has no write site in any other source file. It is not a field in `solute_config_t` and not in the TOML adapter. The TOML pipeline silently drops this parameter. Case 5 may not notice because the legacy `.swp` path sets it. This needs to be resolved as part of the solute config completeness audit.

6. **Large swath of solute physics globals not in the TOML adapter**: `cref`, `kf`, `frexp`, `gampar`, `decpot`, `fdepth`, `ddif`, `kfsat`, `poros`, `daquif`, `decsat`, `cseeptab`, `swbr` are all consumed by `solute.f90` but absent from `solute_config_t` and `config_to_variables.f90`. The TOML case 5 currently relies on these being set by the legacy `.swp` reader path. Before state migration can complete, these fields must be added to `solute_config_t` and the adapter, or the migration will leave the TOML pipeline with silently zeroed physics parameters.

7. **`cml` seeding at config time vs. runtime**: `config_to_variables.f90:629–635` populates `cml(k)` from a CSV file when `swsolu=1`. This same array is then used as input to `solute(task=1, state)` which re-interpolates via `afgen` to the actual node depths (`solute.f90:54`). This double-write is intentional (CSV gives the depth-vs-conc table; `solute` init turns it into per-node values). Migration must preserve this two-stage seed: config writes initial points, `solute(1)` interpolates to all nodes. The `nconc` and `zc` arrays support this but are also not in `solute_config_t`.

8. **`solute` reads `state%drainage%qdra` and `state%surfacewater%qdrtot`** — these are already typed-state reads (not legacy globals). This is correct per Phase 2 work and should be carried forward unchanged in solute state migration.

9. **`swinco` flag gate in `solute(task=1)` and `AgeTracer(task=1)`**: `swinco` controls whether initial concentrations are read from the CSV (`swinco != 3` triggers the afgen interpolation; `swinco == 3` uses cml directly). `swinco` is a general soil-initial-condition switch, not a solute-specific config. It is borrowed from the core config.

10. **`swbr` gate for aquifer breakthrough**: when `swbr=1`, `solute` computes `cdrain` evolution using a mixed reservoir model (`solute.f90:242–255`). This uses `daquif`, `decsat`, `kfsat`, `poros` — all legacy globals. The breakthrough model is a significant physics block not covered by the TOML adapter at all.

## 9. Test surface

| Test file | What it covers |
|-----------|----------------|
| `tests/unit/config/test_solute_config.pf` | `solute_config_validate`: sentinel (`swsolu=0`), happy path (case 5 values), invalid `swsolu` enum |
| `tests/unit/io/toml/test_read_solute_toml.pf` | `read_solute_toml`: absent section, `swsolu=0`, minimal `swsolu=1`, scalar `ldis`, array `ldis`, malformed `pertabsolu` error path |
| `tests/swap-cases/toml/5.salinitystress/` | End-to-end: TOML path with `swsolu=1`; covers cml profile, CSV concentration node output at `-5/-25/-55 cm`, `SQPREC/SQIRRIG/SQBOT/SQDRA/DECTOT/ROTTOT/SAMPRO` CSV columns |
| `tests/swap-cases/5.salinitystress/` | Legacy `.swp` path with full solute physics (case 5 reference) — important parity baseline |

There is **no pFUnit unit test for the compute logic** of `solute.f90` or `AgeTracer`. Unit test coverage is limited to config parsing and validation. The main runtime exercise is the regression case 5.

## 10. Summary statistics

- Total LoC across home files: **867** (548 + 132 + 187)
- Owned globals (Section 2): **37** distinct variables (including 9 AgeTracer-specific: `Ageirr`, `Agedrain`, `Agepre`, `Agepond`, `Agepondm1`, `icAgetopdwn`, `icAgetopupw`, plus `AgeGwl1m`, `icAgeBot`, `icAgeDra`, `icAgeRot`, `icAgeSur`)
- Borrowed globals (Section 3): **47** variables from 10+ external subsystems
- External reader files (Section 3.5 union of all 5 categories): **7** files (`swapoutput.f90`, `swap_csv_output.f90`, `rootextraction.f90`, `irrigation.f90`, `config_to_variables.f90`, `soilhydraulics.f90` [ArMpSs co-writer], `boundtop.f90` [ArMpSs co-writer])
- Entry points (Section 4): **2** public subroutines (`solute`, `AgeTracer`), each with task=1 and task=2
- Output coupling sites (Section 6): **7** output writers reading solute fields
- Cross-subsystem hazards (Section 8): **10**
