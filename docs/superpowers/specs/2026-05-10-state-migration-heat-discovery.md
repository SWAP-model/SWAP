# Subsystem Migration Discovery: Heat

**Date:** 2026-05-10
**Status:** discovery (read-only inventory)
**Migration #:** 4 of N (surfacewater + drainage + solute already shipped; cumulative-reset cohort refactor merged)
**Branch:** `refactor/meteo-state` (created from development; note: branch name differs from the convention `refactor/heat-state` — heat migration will proceed on this branch or a rename)
**Predecessor playbooks:**
- `docs/superpowers/specs/state-migration-playbook.md` (living document)
- ADRs 0030, 0031, 0032, 0033

---

## 1. Source files

| File | LoC | Role |
|------|-----|------|
| `src/heat/temperature.f90` | 497 | Main compute module: `temperature(task)` (task 1 = init, task 2 = rate) and `devries(theta, HeaCap, HeaCon)` (internal thermal-property calculator). Uses bare `use variables` in `temperature` subroutine (no `only:` guard); `devries` uses `use variables, only: NumNod, THETAS, FQUARTZ, FCLAY, FORG`. Does NOT yet take `state` — signature is `subroutine temperature(task)`. |
| `src/heat/frozencond.f90` | 307 | Frozen-soil module: `FrozenCond()` (computes rfcp, zfrostbot, zfrosttop) and `FrozenBounds(state)` (modifies drainage and bottom fluxes under frost). `FrozenBounds` already takes `state intent(inout)` per SS-DRST Phase 2 Task 1; `FrozenCond` still uses bare `use variables` with no state arg. |
| `src/config/heat_config.f90` | 142 | Typed config struct `heat_config_t`. Fields: `swhea`, `swcalt`, `swtopbhea`, `swbotbhea`, `tfroststa`, `tfrostend`, `psand(:)`, `psilt(:)`, `pclay(:)`, `porg(:)`, `tsoil_init(:,:)`. Has `validate` + `finalize`. |
| `src/io/toml/read_heat_toml.f90` | 173 | TOML reader: populates `heat_config_t` from `[heat]` section. Handles scalar scalars and per-layer arrays; decodes `tsoil_init` as array-of-arrays. Returns silently when section absent. |
| `tests/unit/config/test_heat_config.pf` | 158 | pFUnit tests for `heat_config_validate`: sentinel (`swhea=0`), invalid `swhea`, invalid `swcalt`, happy path (numerical), mismatched texture arrays. |
| `tests/unit/io/toml/test_read_heat_toml.pf` | 154 | pFUnit tests for `read_heat_toml`: absent section, `swhea=0`, `swhea=1` minimal, mismatched arrays, malformed `tsoil_init`. |
| `tests/unit/io/toml/fixtures/heat_swhea0.toml` | — | Fixture: absent/sentinel case. |
| `tests/unit/io/toml/fixtures/heat_swhea1.toml` | — | Fixture: full `swhea=1` case. |
| `tests/unit/io/toml/fixtures/heat_mismatched_arrays.toml` | — | Fixture: mismatched psand/pclay/porg lengths. |
| `tests/unit/io/toml/fixtures/heat_malformed_tsoil.toml` | — | Fixture: bad `tsoil_init` for error-path test. |
| `tests/unit/io/toml/fixtures/swap_with_heat.toml` | — | Fixture: full swap TOML with `[heat]` section. |
| `tests/unit/heat/` | — | Directory exists (`.gitkeep` only); no heat compute unit tests yet. |
| `tests/swap-cases/toml/1.hupselbrook/` | — | TOML regression case: `swhea=1`, `swcalt=2` (numerical). Covers `tsoil` profile output via `.tem` file. |
| `tests/swap-cases/toml/2.grassgrowth/` | — | TOML regression case: `swhea=1`, `swcalt=2`. |
| `tests/swap-cases/toml/3.macroporeflow/` | — | TOML regression case: `swhea=1`, `swcalt=2`. |
| `tests/swap-cases/toml/4.oxygenstress/` | — | TOML regression case: `swhea=1`, `swcalt=2`. Heat required for oxygen-stress physiology. |
| `tests/swap-cases/toml/5.salinitystress/` | — | TOML regression case: `swhea=1`, `swcalt=2`. `tsoil` used for temperature-dependent decomposition factor in solute. |

---

## 2. Owned globals (this subsystem writes)

`temperature.f90:subroutine temperature` uses bare `use variables`. The owned set is determined by assignment sites in both source files.

**No flzerointr or flzerocumu blocks exist anywhere in the heat home tree.** All heat globals are recomputed fresh each timestep (tsoil via the numerical/analytical solver, rfcp from tsoil, heacap/heacon via devries). There are NO intermediate or cumulative accumulators.

| Variable | Type | `variables.f90` line | Description | Reset cadence | Cohort | Activity gate | Write sites |
|----------|------|----------------------|-------------|---------------|--------|---------------|-------------|
| `tsoil(macp)` | `real(8)` | 1009 | Soil temperature (°C) per compartment | Instantaneous — recomputed each timestep | flat | `flTemperature` (swhea=1) | `temperature.f90:111,122,234` (via tridag), `243`; also written at init-time from TOML adapter by `config_to_variables.f90:613,937` |
| `heacap(macp)` | `real(8)` | 1012 | Heat capacity (J/cm³/K) per compartment | Instantaneous — recomputed each step (numerical only) | flat | `flTemperature` + `swcalt=2` | `temperature.f90:devries:431,486` (via devries called at 193) |
| `heacon(macp)` | `real(8)` | 1013 | Heat conductivity (J/cm/K/d) per compartment | Instantaneous — recomputed each step | flat | `flTemperature` + `swcalt=2` | `temperature.f90:158,194,196` (also reads heacon at 203,204,210,211,221,227 in tridag setup); `devries:439,450,480,489` |
| `fquartz(macp)` | `real(8)` | 996 | Gravimetric sand+silt fraction per node | Instantaneous — written once at init (task=1), read-only thereafter | flat | `flTemperature` + `swcalt=2` | `temperature.f90:133` (task=1 only) |
| `fclay(macp)` | `real(8)` | 994 | Gravimetric clay fraction per node | Instantaneous — written once at init | flat | `flTemperature` + `swcalt=2` | `temperature.f90:134` (task=1 only) |
| `forg(macp)` | `real(8)` | 995 | Gravimetric organic matter fraction per node | Instantaneous — written once at init | flat | `flTemperature` + `swcalt=2` | `temperature.f90:135` (task=1 only) |
| `tetop` | `real(8)` | 1010 | Temperature (°C) at top of soil profile (under snow cover) | Instantaneous — recomputed each step (numerical only) | flat | `flTemperature` + `swcalt=2` | `temperature.f90:151,161,163,167,169` |
| `tebot` | `real(8)` | 1002 | Temperature (°C) at bottom of soil profile | Instantaneous — recomputed each step (numerical only) | flat | `flTemperature` + `swcalt=2` + `swbotbhea=2` | `temperature.f90:176,179` |
| `rfcp(macp)` | `real(8)` | 924 | Reduction factor for frozen conditions per node (0–1) | Instantaneous — recomputed each step | flat | `swfrost=1` (ungated on flTemperature directly; FrozenCond called unconditionally in the swfrost=1 branch) | `frozencond.f90:80,83,85,88` (FrozenCond) |
| `zfrostbot` | `real(8)` | 971 | Depth of bottom of frozen layer (L) | Instantaneous — recomputed each step | flat | `swfrost=1` | `frozencond.f90:96,104` (FrozenCond) |
| `zfrosttop` | `real(8)` | 972 | Depth of top of frozen layer (L) | Instantaneous — recomputed each step | flat | `swfrost=1` | `frozencond.f90:97,121,123,129,133` (FrozenCond) |
| `nodfrostbot` | `integer` | 665 | Node number of deepest frozen node | Instantaneous — recomputed each step | flat | `swfrost=1` | `frozencond.f90:95,108` (FrozenCond) |
| `fltemperature` | `logical` | 1014 | Flag: soil heat flow active | Config-time set, never resets | flat (logical flag) | self (set if `swhea=1`) | `timecontrol.f90:116,117` (set at run init, not heat-file write) |

**Summary: all 13 owned globals are instantaneous (cadence: step). Zero intermediate. Zero cumulative.** No cohort sub-records needed — all fields go flat on `heat_state_t`. The cohort refactor from ADR 0033 has no bearing on this subsystem.

**Note on `fltemperature`:** written by `timecontrol.f90`, not by the heat home tree. It is a derived flag (from `swhea`), analogous to `flSolute` from `timecontrol`. It is listed here because it gates all heat compute, but ownership is `timecontrol` / `core`.

**Note on `rfcp` vs `swfrost` gating:** `FrozenCond()` is called unconditionally (not behind `if (flTemperature)`) at `swap.f90:284`. It is protected only by `swfrost=1`. However, `rfcp` is always initialized to 1.0 at the start of each `soilwater()` call (`soilhydraulics.f90:884`), so when `swfrost=0`, FrozenCond still runs and resets rfcp to 1.0 — no net effect.

---

## 3. Borrowed globals (this subsystem reads, owned elsewhere)

Variables read by the heat home files but not written by them. `state%drainage%*` and `state%surfacewater%*` reads in `frozencond.f90` are typed-state reads already — not listed here.

| Variable | Best-guess owner subsystem | Read sites (file:line) |
|----------|----------------------------|------------------------|
| `numnod` | core/grid | `temperature.f90:110,121,182,195,209,217,241`; `frozencond.f90:79,100,201,240` |
| `dt` | core/timecontrol | `temperature.f90:203,204,210,211,221,222,223,227,228,229,230` |
| `daynr` | core/timecontrol | `temperature.f90:111,243` |
| `t1900` | core/timecontrol | `temperature.f90:151,179` |
| `z(macp)` | core/grid | `temperature.f90:112,120,123,244`; `frozencond.f90:104,123,129` |
| `dz(macp)` | core/grid | `temperature.f90:159,203,204,210,211,221,222,227,228` |
| `disnod(macp+1)` | core/grid | `temperature.f90:203,204,210,211,221,227` |
| `layer(macp)` | core/grid | `temperature.f90:130` |
| `theta(macp)` | soilhydraulics | `temperature.f90:189` |
| `thetm1(macp)` | soilhydraulics | `temperature.f90:189` |
| `thetas(macp)` | soilhydraulics | `frozencond.f90:207` (via devries: `temperature.f90:415`); |
| `swnco` / `swinco` | core/config | `temperature.f90:116` |
| `swcalt` | heat config (legacy global, set by adapter) | `temperature.f90:108,116,127,145`; `swapoutput.f90:1840` |
| `swtopbhea` | heat config (legacy global) | `temperature.f90:149` |
| `swbotbhea` / `SwBotbHea` | heat config (legacy global) | `temperature.f90:174,218,224` |
| `swhea` | heat config (legacy global, read by timecontrol) | `timecontrol.f90:117` |
| `swfrost` | soil/frost config (legacy global) | `frozencond.f90:81`; `soilhydraulics.f90:884` region |
| `nheat` | heat config (legacy global, set by adapter) | `temperature.f90:117` |
| `zh(macp)` | heat config (legacy global, set by adapter) | `temperature.f90:119` |
| `tmean` | heat config (legacy global, NOT in typed config or adapter) | `temperature.f90:111,243` |
| `tampli` | heat config (legacy global, NOT in typed config or adapter) | `temperature.f90:111,243` |
| `timref` | heat config (legacy global, NOT in typed config or adapter) | `temperature.f90:111,243` |
| `ddamp` | heat config (legacy global, NOT in typed config or adapter) | `temperature.f90:112,244` |
| `temtoptab(mabbc*2)` | heat config (legacy global, NOT in typed config or adapter) | `temperature.f90:151` |
| `tembtab(mabbc*2)` | heat config (legacy global, NOT in typed config or adapter) | `temperature.f90:179` |
| `tfroststa` | heat config (set by adapter) | `frozencond.f90:82,84,88` |
| `tfrostend` | heat config (set by adapter) | `frozencond.f90:84,86,88,103,118` |
| `psand(maho)` | heat config (set by adapter) | `temperature.f90:133` |
| `psilt(maho)` | heat config (set by adapter) | `temperature.f90:133` |
| `pclay(maho)` | heat config (set by adapter) | `temperature.f90:134` |
| `orgmat(maho)` | soil config / heat adapter (dual-set — see Hazard #3) | `temperature.f90:131` |
| `ssnow` | snow subsystem | `temperature.f90:152` |
| `flmetdetail` | meteo subsystem | `temperature.f90:160,166` |
| `atav(wrecord)` | meteo subsystem | `temperature.f90:161,167` |
| `wrecord` | meteo subsystem | `temperature.f90:161,167` |
| `Tav` | meteo subsystem | `temperature.f90:163,169` |
| `macp` | core/dimensions | `temperature.f90:95,96,97,98,99` (dimension constant) |
| `mabbc` | core/dimensions | `temperature.f90:95,151,179` |
| `swtem` | heat config (legacy output switch) | `swapoutput.f90:1843,1852` |
| `swini` | core/config | `swapoutput.f90:1840` |
| `swdra` | drainage config | `frozencond.f90:220` |
| `swdivd` | drainage config | `frozencond.f90:276` |
| `nrlevs` | drainage | `frozencond.f90:233,251,260,286` |
| `zbotdr(madr)` | drainage | `frozencond.f90:234,252,265` |
| `gwl` | soilhydraulics | `frozencond.f90:277` |
| `L` | drainage | `frozencond.f90:280` |
| `cofani(maho)` | drainage | `frozencond.f90:249` |
| `ksatexm(maho)` | soilhydraulics | `frozencond.f90:242` |
| `ksatfit(maho)` | soilhydraulics | `frozencond.f90:245` |
| `fluseksatexm(macp)` | soilhydraulics | `frozencond.f90:241` |
| `qbot_nonfrozen` | soilhydraulics | `frozencond.f90:200` |
| `qbot` | soilhydraulics | `frozencond.f90:200,222,266,268` (FrozenBounds also WRITES qbot) |
| `owltab` | drainage | `frozencond.f90:281` |
| `Swdivdinf` | drainage | `frozencond.f90:280` |
| `Swnrsrf` | drainage | `frozencond.f90:280` |
| `SwTopnrsrf` | drainage | `frozencond.f90:280` |
| `FacDpthInf` | drainage | `frozencond.f90:281` |

**Note:** `forg(i) = dummy*gmineral/1.4d0` at `temperature.f90:135` writes `forg` (heat-owned). `orgmat(lay)` at `temperature.f90:131` is read-only in temperature.f90; `orgmat` is owned by soil/heat config adapter.

**Note on `qbot`:** `FrozenBounds` reads `qbot` (via `qbot = qbot_nonfrozen` at line 200) and may also zero it (lines 222, 266). `qbot` is primarily a soilhydraulics variable. This is the same consumer-writes-borrower pattern as `FrozenBounds` writing `qdrain` / `qdra` for drainage.

---

## 3.5 External readers of owned globals — 5-category framework

Grep templates applied for each owned global. Alias-form grep (`=> state%heat`) returns empty — no state%heat exists yet.

```bash
# Categories 1+2 (use variables clauses):
grep -rEn "use variables.*\b<owned-var>\b" src/ --include="*.f90" \
  | grep -v "src/core/variables.f90" \
  | grep -v "src/state/" \
  | grep -v "src/heat/" \
  | grep -v "src/config/heat_config" \
  | grep -v "src/io/toml/read_heat_toml"

# Categories 3+4+5 (raw symbol references):
grep -rEn "\b<owned-var>\b" src/ --include="*.f90" \
  | grep -v "src/core/variables.f90" \
  | grep -v "src/core/initialize.f90" \
  | grep -v "src/state/" \
  | grep -v "state%" \
  | grep -v "config%" \
  | grep -v "src/heat/"

# Alias form:
grep -rEn "=> state%heat" src/ --include="*.f90"
```

### Results per owned global

**`tsoil(macp)`** — most heavily read heat global:

- **Output (Cat 1):** `swapoutput.f90:578` (`outvap`: `use variables, only: ...,tsoil,...`; writes `tsoil(node)` per-node); `swapoutput.f90:1049` (`outend`: reads `tsoil(i)` for final profile snapshot); `swapoutput.f90:2048` (`outtem`: reads `tsoil(i)` for `.tem` file); `swapoutput.f90:4661` (CSV tz module: `use variables, only: ...,tsoil,...`; `4762` writes per node)
- **Output (Cat 1, CSV tz):** `swapoutput.f90:4037,4302` (CSV inline block — writes `tsoil(Nodes_TEMP(j))`)
- **Output (Cat 1, afo/bfo):** `swapoutput.f90:2458,2519,2626,2751,2807,2901` (formatted/unformatted output routines pass `tsoil` to binary output blocks)
- **Output (Cat 1, swap_csv_output):** `swap_csv_output.f90:954,1090` (`use variables, only: ...,tsoil,...`; writes `tsoil(j)` per node)
- **Compute (Cat 2):** `solute/solute.f90:206,207` (temperature-dependent decomposition factor: `if (tsoil(i) .lt. 35.0)…ftemp = exp(gampar*(tsoil(i)-20.0))`); `crop/rootextraction.f90:144,190` (oxygen stress ODE + frost root cutoff: `if (swfrost .eq.1 .and. tsoil(node) .lt. 0.0)`); `crop/oxygenstress.f90:214,1550` (soil temperature for Arrhenius: `soil_temp = tsoil(node)+273.d0`); `crop/management_soil.f90:131` (soil temperature depth-integral for tillage decision); `crop/cropgrowth.f90:876,4609,4640` (sowing temperature sum, crop phenology); `atmosphere/snow.f90:110` (snow melt threshold: `if (tsoil(1) .gt. 0.5)`); `utils/soilhydraulicsutils.f90:457` (temperature-dependent hydraulic conductivity model 04_11)
- **Call-site arg (Cat 5):** `soil/soilgrid.f90:179` (passed as explicit argument to `ConvertDiscrVert` for grid-change remapping of temperature profile); `crop/oxygenstress.f90:1529,1541` (passed as argument to `OxygenReproFunction`)
- **Init-seed (Cat 4):** `io/toml/config_to_variables.f90:613,937` (seeds `tsoil(k)` from CSV file and from `config%heat%tsoil_init` at config time)

**`heacap(macp)`:**

- **Output (Cat 1):** `swapoutput.f90:4661` (CSV tz module: `use variables, only: ..., HEACAP,...`; `4767` writes `HEACAP(j)/1.0d-6`); `swap_csv_output.f90:12,487,1095` (module use and writes `HEACAP(nodes)` to CSV)
- **Working-buffer / Compute-writeback (Cat 3 hazard):** `swapoutput.f90:1914` — `outheapar()` calls `devries(thetadum, heacap, heacnd)` WITH the GLOBAL `heacap` as output argument. `outheapar` uses bare `use variables` (line 1879) and declares `heacnd` locally but NOT `heacap`. This means `outheapar` OVERWRITES the global `heacap` array each time it is called. Called only once at init (`TemperatureOutput(1)` when `swini=1 .and. swcalt=2`), but this is still an output routine mutating heat state — matches playbook gotcha #6. See Hazard #1.

**`heacon(macp)`:**

- **Output (Cat 1):** `swapoutput.f90:4661` (CSV tz module: `use variables, only: ..., HEACON,...`; `4768` writes `HEACON(j)/864.0d0`); `swap_csv_output.f90:12,954,1095` (module use and writes)

**`fquartz(macp)`, `fclay(macp)`, `forg(macp)`:**

- No external readers found outside the heat home tree. Only read within `devries` (which is in the heat home tree). These composition arrays are internal to the heat numerical method.

**`tetop`:**

- **Output (Cat 1):** `swapoutput.f90:2048` (`outtem`: `use variables, only: ...,tetop,...`; `2089,2107` write `tetop` to `.tem`); `swap_csv_output.f90:10,168,304` (`use variables` import, header `TETOP`, `set_values` maps `TeTop`)

**`tebot`:**

- **Output (Cat 1):** `swapoutput.f90:2048` (`outtem`; `2090,2108` write `tebot`); `swap_csv_output.f90:10,169,305` (same module, `TEBOT` column)

**`rfcp(macp)`:**

- **Compute (Cat 2):** `soil/soilhydraulics.f90:113,130,170,238,262,318,453,520,548,1050,1170` (reads `rfcp(i)` as argument to `hconduc`/`dhconduc` — scales hydraulic conductivity); `boundary/boundtop.f90:105,146,148` (reads `rfcp(1)` for surface K and ks at top boundary); `boundary/boundbottom.f90:159` (reads `rfcp(numnod)` for bottom mean K)
- **Co-write (Cat 3):** `soil/soilhydraulics.f90:884` (inside `soilwater` init block: `rfcp(i) = 1.0d0` for all nodes — resets rfcp to unfrozen default before each Richards-step). This is a cross-subsystem co-write; soilhydraulics resets rfcp, then FrozenCond overwrites with the actual frozen-state values. The ordering is enforced by the call sequence in `swap.f90` (soilwater → FrozenCond → FrozenBounds).
- **Call-site arg (Cat 5):** `utils/soilhydraulicsutils.f90:327,332,399,405,410,526` (`rfcp` is a function argument to `hconduc` and `dhconduc` — these are pure utility functions receiving rfcp by value)

**`zfrostbot`, `zfrosttop`, `nodfrostbot`:**

- No external readers outside the heat home tree found. `FrozenBounds` reads them internally to determine which drains fall inside the frozen zone (via `zfrostbot.lt.zbotdr(level)` checks, but `FrozenBounds` is in the heat home tree).

**`fltemperature`:**

- **Call-site gate (Cat 2):** `core/swap.f90:66,196,217,284` (used to gate `Temperature(1)`, `TemperatureOutput(1)`, and `Temperature(2)` calls); `solute/solute.f90:205` (gates temperature-dependent decomposition branch: `if (fltemperature) then`); `swapoutput.f90:1049,1110` (`outend` reads `fltemperature` to decide whether to output the tsoil profile in `.end` file)

### Section 3.5 tally

| Owned global | Output readers (Cat 1) | Compute readers (Cat 2) | Working-buffer / co-write (Cat 3) | Init-seed (Cat 4) | Call-site arg (Cat 5) |
|---|---|---|---|---|---|
| `tsoil` | `swapoutput.f90` (outvap, outend, outtem, CSV-tz, afo/bfo), `swap_csv_output.f90` | `solute.f90`, `rootextraction.f90`, `oxygenstress.f90`, `management_soil.f90`, `cropgrowth.f90`, `snow.f90`, `soilhydraulicsutils.f90` | — | `config_to_variables.f90` | `soilgrid.f90` (ConvertDiscrVert), `oxygenstress.f90` (OxygenReproFunction) |
| `heacap` | `swapoutput.f90` (CSV-tz), `swap_csv_output.f90` | — | `swapoutput.f90:outheapar` (WRITES global heacap — hazard!) | — | — |
| `heacon` | `swapoutput.f90` (CSV-tz), `swap_csv_output.f90` | — | — | — | — |
| `fquartz`, `fclay`, `forg` | — | — | — | — | — |
| `tetop` | `swapoutput.f90` (outtem), `swap_csv_output.f90` | — | — | — | — |
| `tebot` | `swapoutput.f90` (outtem), `swap_csv_output.f90` | — | — | — | — |
| `rfcp` | — | `soilhydraulics.f90`, `boundtop.f90`, `boundbottom.f90` | `soilhydraulics.f90:884` (resets rfcp to 1.0 each step) | — | `soilhydraulicsutils.f90` (hconduc/dhconduc args) |
| `zfrostbot`, `zfrosttop`, `nodfrostbot` | — | — | — | — | — |
| `fltemperature` | — | `swap.f90` (call-site gate), `solute.f90` (branch gate), `swapoutput.f90` (outend gate) | — | — | — |

**External reader files (union across all 5 categories):** 12 distinct files:
1. `src/io/swapoutput.f90` (Cat 1: outvap, outend, outtem, outheapar, CSV-tz, afo/bfo)
2. `src/io/swap_csv_output.f90` (Cat 1: set_values, per-node time-series)
3. `src/solute/solute.f90` (Cat 2: decomposition temperature gate)
4. `src/crop/rootextraction.f90` (Cat 2: frost root cutoff, oxygen ODE)
5. `src/crop/oxygenstress.f90` (Cat 2: Arrhenius temperature; also Cat 5)
6. `src/crop/management_soil.f90` (Cat 2: tillage temperature sum)
7. `src/crop/cropgrowth.f90` (Cat 2: sowing temp sum, phenology)
8. `src/atmosphere/snow.f90` (Cat 2: snow melt threshold)
9. `src/utils/soilhydraulicsutils.f90` (Cat 2: temp-dependent K model 04_11; Cat 5: hconduc/dhconduc args)
10. `src/soil/soilhydraulics.f90` (Cat 2: hconduc calls; Cat 3: rfcp co-write)
11. `src/boundary/boundtop.f90` (Cat 2: surface K)
12. `src/boundary/boundbottom.f90` (Cat 2: bottom K)
13. `src/soil/soilgrid.f90` (Cat 5: ConvertDiscrVert arg)
14. `src/core/swap.f90` (Cat 2: call-site gate on fltemperature)
15. `src/io/toml/config_to_variables.f90` (Cat 4: tsoil init-seed)

---

## 4. Entry points (subroutines called from outside this subsystem)

| Subroutine (this subsystem) | Called from (file:line) | Purpose / lifecycle stage | Takes `state`? |
|-----------------------------|-------------------------|--------------------------|----------------|
| `Temperature(1)` | `swap.f90:196` | Initialization: set initial tsoil profile, initialize fquartz/fclay/forg (numerical only) | NO — bare task-int signature |
| `Temperature(2)` | `swap.f90:324` | Rate: update tsoil (analytical or numerical solver), compute heacap/heacon, set tetop/tebot | NO — bare task-int signature |
| `FrozenCond()` | `swap.f90:284` | Rate: compute rfcp(node), zfrostbot, zfrosttop, nodfrostbot from tsoil | NO — bare no-arg signature |
| `FrozenBounds(state)` | `swap.f90:302` (gated `SwFrost.eq.1`); `swapoutput.f90:3754` (mini-simulation arm) | Rate: modify drainage/bottom fluxes under frost using rfcp | YES — already plumbed with `state intent(inout)` |
| `TemperatureOutput(1)` | `swap.f90:217` (gated `flTemperature`) | Init output: open `.tem` file; if `swini=1 .and. swcalt=2` call `outheapar()` | NO — bare task-int signature |
| `TemperatureOutput(2)` | `swap.f90:378` | Periodic output: write `tsoil`, `tetop`, `tebot`, `tav` to `.tem` | NO |
| `TemperatureOutput(3)` | `swap.f90:427` | Final output: close `.tem` | NO |

**Key asymmetry:** `FrozenBounds` already takes `state intent(inout)` (migrated in SS-DRST Task 3). `Temperature(task)` and `FrozenCond()` do NOT — they still write all heat globals directly via `use variables`. This is the primary migration target.

---

## 5. Internal call graph

```
temperature_mod
├── subroutine temperature(task)
│   ├── use variables          [bare wildcard import]
│   ├── use array_utils, only: afgen
│   ├── use numericalsolvers_mod, only: tridag
│   │
│   ├── task=1 (init)
│   │   ├── if (swcalt=1): tsoil(i) = analytical formula     [writes tsoil]
│   │   ├── if (swcalt=2 .and. swinco≠3):
│   │   │   └── afgen(tab, macp*2, z(i)) → tsoil(i)         [writes tsoil]
│   │   └── if (swcalt=2):
│   │       └── fquartz(i), fclay(i), forg(i) from psand/psilt/pclay/orgmat
│   │                                                          [writes fquartz,fclay,forg]
│   │
│   └── task=2 (rate)
│       ├── if (swcalt=2):
│       │   ├── afgen(temtoptab,...) or Tav/atav → TeTop     [writes tetop]
│       │   ├── Tsoil(Numnod) or afgen(tembtab,...) → TeBot  [writes tebot]
│       │   ├── call devries(theave, heacap, heacnd)          [writes heacap via devries]
│       │   ├── heacon(1) = heacnd(1)                        [writes heacon]
│       │   ├── heacon(i) = 0.5*(heacnd(i)+heacnd(i-1))     [writes heacon]
│       │   └── call tridag(..., tsoil, ierror)               [writes tsoil]
│       └── if (swcalt=1):
│           └── tsoil(i) = analytical formula                 [writes tsoil]
│
└── subroutine devries(theta, HeaCap, HeaCon)
    ├── use variables, only: NumNod, THETAS, FQUARTZ, FCLAY, FORG
    ├── [pure computation — no external calls]
    ├── HeaCap(Node) = volume-weighted heat capacity          [writes heacap]
    └── HeaCon(Node) = de Vries thermal conductivity         [writes heacon]
        (with unit conversions: J/m³/K → J/cm³/K, W/m/K → J/cm/K/d)

frozencond_mod
├── subroutine FrozenCond()
│   ├── use variables          [bare wildcard import]
│   ├── rfcp(node) = 1.0 to 0.0 based on tsoil vs tfroststa/tfrostend
│   ├── nodfrostbot = …        [writes nodfrostbot]
│   ├── zfrostbot = …          [writes zfrostbot]
│   └── zfrosttop = …          [writes zfrosttop]
│
└── subroutine FrozenBounds(state)
    ├── use variables          [bare wildcard import]
    ├── use distribute_drainage, only: DIVDRA
    ├── use swap_state_mod, only: swap_state_t
    ├── reads: rfcp, thetas, theta, nodfrostbot, nrlevs, zbotdr, swdra, swdivd, gwl
    ├── writes: qbot                                         [borrowed global write]
    ├── associate(qdrain => state%drainage%qdrain, qdra => state%drainage%qdra)
    │   └── modifies state%drainage%qdrain, state%drainage%qdra
    ├── writes: state%surfacewater%qdrtot
    └── call divdra(…) [if swdivd=1]                        [redistribute drainage]
```

---

## 6. Output coupling

| Output writer (file:subroutine) | Reads which heat fields | Destination file |
|---------------------------------|------------------------|------------------|
| `swapoutput.f90:outvap` | `tsoil(node)` | `*.vap` (soil profile time-series) |
| `swapoutput.f90:outend` | `tsoil(i)` (gated `fltemperature`) | `*.end` (final state snapshot) |
| `swapoutput.f90:outtem` | `tsoil(i)`, `tetop`, `tebot` | `*.tem` (daily soil temperature) |
| `swapoutput.f90:outheapar` | WRITES `heacap` (global!) via `devries(thetadum, heacap, heacnd)`, then reads `heacap(node)` to write | `heatparam.csv` |
| `swapoutput.f90` (CSV tz block) | `tsoil(j)`, `HEACAP(j)`, `HEACON(j)` | per-timestep CSV (tz) |
| `swapoutput.f90` (afo/bfo blocks, ~2458+) | `tsoil` (passed as argument to binary output blocks) | `*.afo`, `*.bfo` |
| `swap_csv_output.f90:set_values` | `tsoil(nodes)`, `HEACAP(nodes)`, `HEACON(nodes)`, `tetop`, `tebot` | `swap.csv` columns TEMP[n], HEACAP[n], HEACON[n], TETOP, TEBOT |

**Critical: `outheapar` anti-pattern (Hazard #1).** The `outheapar` routine calls `devries(thetadum, heacap, heacnd)` where `heacap` is the MODULE-GLOBAL array (not a local). `heacnd` is declared local but `heacap` is not — it comes from `use variables`. This means the first call to `outheapar` (at init) overwrites the global `heacap` array with at-saturation values across a theta sweep. Since `outheapar` is called only once (`swini=1 .and. swcalt=2`), the damage is limited to that moment, but `Temperature(2)` immediately follows and recomputes `heacap` correctly — so in practice there is no physics error. However, this is a textbook playbook gotcha #6 (output routine mutates state). Migration must resolve this by either (a) using a local `heacap_local(macp)` in `outheapar`, or (b) having `outheapar` call `devries` into a local array and not touch `state%heat%heacap`.

---

## 7. Config inputs

The TOML adapter section for heat in `src/io/toml/config_to_variables.f90` (lines 885–937):

| TOML key | Typed-config field path | Legacy global written | Used at (file:line) |
|---|---|---|---|
| `heat.swhea` | `config%heat%swhea` | `swhea` | `config_to_variables.f90:885`; `timecontrol.f90:117` (→ `flTemperature`) |
| `heat.swcalt` | `config%heat%swcalt` | `swcalt` | `config_to_variables.f90:886`; `temperature.f90:108,116,127,145` |
| `heat.swtopbhea` | `config%heat%swtopbhea` | `swtopbhea` | `config_to_variables.f90:887`; `temperature.f90:149` |
| `heat.swbotbhea` | `config%heat%swbotbhea` | `swbotbhea` | `config_to_variables.f90:888`; `temperature.f90:174,218,224` |
| `heat.tfroststa` | `config%heat%tfroststa` | `tfroststa` | `config_to_variables.f90:889`; `frozencond.f90:82,84,88` |
| `heat.tfrostend` | `config%heat%tfrostend` | `tfrostend` | `config_to_variables.f90:890`; `frozencond.f90:84,86,88,103,118` |
| `heat.psand` | `config%heat%psand(:)` | `psand(i)` | `config_to_variables.f90:892–895`; `temperature.f90:133` |
| `heat.psilt` | `config%heat%psilt(:)` | `psilt(i)` | `config_to_variables.f90:897–900`; `temperature.f90:133` |
| `heat.pclay` | `config%heat%pclay(:)` | `pclay(i)` | `config_to_variables.f90:902–905`; `temperature.f90:134` |
| `heat.porg` | `config%heat%porg(:)` | `forg(i)` via `orgmat` (dual path) | `config_to_variables.f90:907–918`; `temperature.f90:131,135` |
| `heat.tsoil_init` | `config%heat%tsoil_init(:,:)` | `nheat`, `zh(i)`, `tsoil(i)` | `config_to_variables.f90:932–937`; `temperature.f90:117–123` |

**Heat physics globals NOT in the TOML adapter** (used by `temperature.f90` but absent from `heat_config_t` and `config_to_variables.f90`):

| Legacy global | `variables.f90` line | Role | Comment |
|---|---|---|---|
| `ddamp` | 993 | Damping depth (L) for analytical temperature wave | Required for `swcalt=1`. Only set by legacy `.swp` reader. TOML path for `swcalt=1` is silently broken. |
| `tmean` | 1008 | Mean annual temperature (°C) at soil surface | Required for `swcalt=1`. Same issue. |
| `tampli` | 1001 | Amplitude of annual temperature wave (°C) | Required for `swcalt=1`. Same issue. |
| `timref` | 1007 | Time of year (T) with peak temperature wave | Required for `swcalt=1`. Same issue. |
| `temtoptab(mabbc*2)` | 1004 | Prescribed surface temperature table (time → T°C) | Required for `swtopbhea=2`. Not in typed config. |
| `tembtab(mabbc*2)` | 1003 | Prescribed bottom temperature table (time → T°C) | Required for `swbotbhea=2`. Not in typed config. |

All 5 TOML regression cases use `swcalt=2` (numerical method). The `swcalt=1` (analytical) path is not covered by the TOML regression suite and has missing adapter fields.

---

## 8. Open questions / cross-subsystem hazards

1. **Hazard: `outheapar` overwrites global `heacap`.** `swapoutput.f90:1914` calls `devries(thetadum, heacap, heacnd)` with the global `heacap` as the output array. Since `outheapar` uses bare `use variables` and `heacap` is not declared locally, this mutates the global during output. The physics consequence is mild (only at run init, `Temperature(2)` immediately follows), but it is a textbook playbook gotcha #6 anti-pattern. Resolution: use a local `heacap_local(numnod)` in `outheapar` and write only from it. This must happen before or simultaneous with migrating `heacap` to `state%heat`.

2. **Hazard: `FrozenBounds` called from within `swapoutput.f90:3754` (mini-simulation).** `swapoutput.f90` runs an independent mini-simulation arm for output at non-standard groundwater levels. This arm calls `FrozenBounds(state_om)` where `state_om` is a SAVE-local clone state. This is the same pattern as `Drainage(state_om)` in the output mini-simulation (established in SS-DRST). No new wrinkle for heat migration, but must be tracked: when `FrozenBounds` signature changes (if it needs `state%heat` reads), the `state_om` call must also receive the correct heat-state snapshot.

3. **Hazard: `orgmat` is dual-written by soil adapter and heat adapter.** `config_to_variables.f90:538–540` writes `orgmat(i)` from `config%soil%orgmat`. Later at `config_to_variables.f90:916–918`, the heat adapter writes `orgmat(i) = config%heat%porg(i)` if and only if `config%soil%orgmat` is not allocated. This is intentional (heat provides a fallback), but it means `orgmat` has two potential sources. Heat state migration must not migrate `orgmat` into `state%heat` — it is shared soil configuration.

4. **Hazard: `rfcp` co-written by `soilhydraulics.f90:884`.** Before each `soilwater` call, `soilhydraulics` resets `rfcp(i) = 1.0` for all nodes. This is the correct design (initialize to unfrozen, then FrozenCond overwrites). But it means `rfcp` has a cross-subsystem reset writer that is not in the heat home tree. If `rfcp` is migrated to `state%heat%rfcp`, then `soilhydraulics.f90` must be given read access to `state%heat` (or `rfcp` must remain a legacy global for the soilhydraulics reset path until a full heat migration phase completes).

5. **Hazard: `Temperature(task)` does NOT take `state` — it is the only compute entry point in the first three migrations not yet plumbed.** `FrozenBounds(state)` was plumbed in SS-DRST Task 3 (already done). `FrozenCond()` has no state coupling. `Temperature(task)` still uses bare `use variables` for all reads and writes. This subsystem has more work to do in Phase 1 than previous subsystems because the main compute entry point must be threaded with state before any dual-write can start.

6. **Hazard: `swcalt=1` (analytical) path is uncovered and adapter-incomplete.** `ddamp`, `tmean`, `tampli`, `timref` are absent from `heat_config_t` and the TOML adapter. All five TOML regression cases use `swcalt=2`. The analytical path cannot be tested via TOML inputs as currently implemented. This is a pre-existing gap, not created by state migration, but migration must not worsen it.

7. **Hazard: `FrozenCond()` directly reads `tsoil` (owned by heat) without a state argument.** When `tsoil` is migrated to `state%heat%tsoil`, `FrozenCond()` must either receive state or be refactored to accept `tsoil` as an explicit argument. The simplest fix: add `state intent(in)` to `FrozenCond`, matching the pattern already established by `FrozenBounds`.

8. **Hazard: `soilhydraulicsutils.f90:457` uses `tsoil(node)` for temperature-dependent K model (04_11).** This is a compute reader (Cat 2) currently importing `tsoil` via `use variables`. After migration, `soilhydraulicsutils` would need either a `state%heat%tsoil(node)` read path or `tsoil` passed as an explicit argument. Same pattern as `OxygenReproFunction` in oxygenstress.f90 which already receives `tsoil` by argument.

9. **`swfrost` vs `flTemperature` independence.** `FrozenCond()` is called at `swap.f90:284` unconditionally (inside the Richards timestep loop, not gated by `flTemperature`). This means frost reduction of K can operate even when `swhea=0`, as long as `swfrost=1`. In practice this would be incoherent (rfcp depends on tsoil, which would be all-zero if heat is off), but it is a latent design tension. The heat state migration should document the intended lifecycle: `flTemperature` gates temperature computation; `swfrost` gates frost-K reduction independently.

10. **Many crop/solute compute readers of `tsoil` (Cat 2) are not gated by `flTemperature`.** `cropgrowth.f90`, `rootextraction.f90`, `management_soil.f90`, `oxygenstress.f90`, `solute.f90` all read `tsoil` without checking `flTemperature`. When `swhea=0`, `tsoil` is all-zero (initialize.f90 zeros it). These callers silently use zero temperatures, which may produce incorrect results (e.g., solute decomposition at 20°C when `gampar*(0-20)` gives a very small factor, not the intended default). This is a pre-existing issue. After migration, `state%heat%tsoil` will carry the same behavior — but an `is_active` flag on the heat state could make the zero-temperature case more explicit.

---

## 9. Test surface

| Test file | What it covers |
|-----------|----------------|
| `tests/unit/config/test_heat_config.pf` (158 LoC) | `heat_config_validate`: sentinel (`swhea=0`), invalid `swhea`/`swcalt`/`swtopbhea`/`swbotbhea` enum, happy path (numerical method), mismatched texture array lengths |
| `tests/unit/io/toml/test_read_heat_toml.pf` (154 LoC) | `read_heat_toml`: absent section, `swhea=0`, `swhea=1` minimal, mismatched arrays error, malformed `tsoil_init` error |
| `tests/swap-cases/toml/1.hupselbrook/` | End-to-end: `swhea=1`, `swcalt=2`, `swfrost=0`. Full numerical heat transport with soil texture per layer. Covers `tsoil` profile output. |
| `tests/swap-cases/toml/2.grassgrowth/` | End-to-end: `swhea=1`, `swcalt=2`. Same as above with grass growth enabled. |
| `tests/swap-cases/toml/3.macroporeflow/` | End-to-end: `swhea=1`, `swcalt=2`. Same heat path with macropore flow active. |
| `tests/swap-cases/toml/4.oxygenstress/` | End-to-end: `swhea=1`, `swcalt=2`. `tsoil` used in Bartholomeus oxygen-stress model — most heat-sensitive of all regression cases. |
| `tests/swap-cases/toml/5.salinitystress/` | End-to-end: `swhea=1`, `swcalt=2`. `tsoil` used in solute temperature-factor `gampar*(tsoil-20)`. |
| `tests/unit/heat/` | **Empty** (`.gitkeep` only). No unit tests for `temperature.f90` or `FrozenCond`. |

**No pFUnit unit tests for the compute logic of `temperature.f90` or `frozencond.f90`.** Unit test coverage is limited to config parsing and validation. The five regression cases are the only runtime exercise of the heat solver.

**Analytical method (`swcalt=1`) has zero regression coverage** in the TOML path — none of the five TOML cases use it. The legacy `.swp` cases (e.g., `tests/swap-cases/1.hupselbrook/swap_linux.swp.template`) use `swcalt=2` as well. Investigation of `.swp` templates showed all use `SWHEA = 1` + `SWCALT = 2`.

---

## 10. Summary statistics

- **Total LoC across home files:** 1119 (497 + 307 + 142 + 173)
- **Owned globals (Section 2):** 13 variables (including `fltemperature` as a logical flag owned by timecontrol but listed for gating clarity)
  - Instantaneous: 12 (`tsoil`, `heacap`, `heacon`, `fquartz`, `fclay`, `forg`, `tetop`, `tebot`, `rfcp`, `zfrostbot`, `zfrosttop`, `nodfrostbot`)
  - Intermediate: 0
  - Cumulative: 0
- **Activity gates governing cumulative accumulation:** N/A — no cumulative fields. The heat subsystem is cadence-clean: all state is either step-computed or init-time set.
- **Cohort design consequence:** no `heat_intermediate_t` or `heat_cumulative_t` cohort types needed. A single flat `heat_state_t` with all fields covers the full owned set. ADR 0033 partitioning rules do not apply.
- **Borrowed globals (Section 3):** ~50 variables from 10+ subsystems (core/grid, core/timecontrol, soilhydraulics, meteo, snow, drainage, heat-config legacy globals)
- **External reader files (Section 3.5 distinct files):** 15 files across output (2), compute (10), co-write (1), init-seed (1), call-site arg (1) categories
- **Entry points (Section 4):** 3 public subroutines: `Temperature(task)`, `FrozenCond()`, `FrozenBounds(state)`. Plus `TemperatureOutput(task)` in swapoutput.
- **Output coupling sites (Section 6):** 7 output writers reading heat fields
- **Cross-subsystem hazards (Section 8):** 10 items
  - Critical: `outheapar` writes global `heacap` (output mutates state)
  - Critical: `Temperature(task)` takes no state arg — most migration work lives here
  - Critical: `rfcp` co-reset by `soilhydraulics.f90` — cross-subsystem init dependency
  - Notable: `swcalt=1` TOML adapter incomplete (4 missing physics globals)
  - Notable: crop/solute tsoil readers (10 files) ungated on `flTemperature`
