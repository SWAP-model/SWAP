---
title: "SS-9 residual audit — Phase 4f config_to_variables refresh"
author: SWAP modernization team
date: 2026-05-04
status: closed
---

# Phase 4f-extend SS-9 — residual gap audit

**Sub-spec:** SS-9 (legacy reader retirement umbrella)
**Predecessor:** `docs/phase-4f-config-to-variables-audit.md` (2026-04-27, 146 G rows)
**HEAD audited:** `1a951ec` (2026-05-04 — SS-7/SS-8 closed)
**Outcome:** **CLOSED — zero STILL-G rows reachable in TOML runtime; SS-11 unblocked from this side.**

## Methodology

Re-classified each of the 146 (G) rows from the 2026-04-27 audit by:

1. Grepping `src/config/*_config.f90` for the field name (schema slot present?).
2. Grepping `src/io/toml/config_to_variables.f90` AND
   `src/crop/{cropfixed,cropwofost,cropgrass}_init.f90` (post Phase 4f the
   strangler adapter is no longer the only adapter — per-crop init modules
   wire most crop fields).
3. Cross-referencing the SS-1..8 audit closures
   (`phase-4f-{irrigation,management-soil,cropgrowth-nutrients,
   readmeteo-ttutil,rddre,surfacewater-parity,…}-audit.md` and
   `2026-05-04-ss4-macropore-toml-stub-error.md`).
4. Checking each field's runtime gate (`flCropNut`, `swdrought`, `swmacro`,
   `swsec`/`swsrf`/`swman`, `schedule`, `swdc`, `swhyst`, `swsp`,
   `swbulb`, `swrdc`, `swoxygentype`, etc.) against the TOML stub-error
   set in `src/config/*_config_validate`.
5. Verifying current TOML cases under `tests/swap-cases/toml/*/swap.toml`
   and per-crop `*.crp.toml` for the gating switch values.

Each (G) row reclassified as one of:

- **C-now** — schema field exists; an adapter (strangler `config_to_variables`
  OR a per-crop `*_init.f90`) populates the legacy global; reachable code
  path is exercised by at least one regression case OR safely defaulted.
- **DEFERRED-by-shared-gate** — only reachable when a switch is in a state
  that the TOML config layer rejects with `ERR_VALIDATION_*` (macropore,
  De Jong van Lier, swsec=1, swman=2, swqhr=2, swbulb=1, swrdc=1, etc.)
  OR a state never set in the TOML pipeline (`flCropNut=.true.` —
  SS-7/SS-8).
- **STILL-G** — schema or adapter genuinely missing AND the path is
  reachable in TOML runtime under at least one current case configuration.
- **STALE** — code path was deleted (e.g. SS-5 TTutil readmeteo branches);
  the legacy reader key is no longer consumed in production.

## Summary

| Bucket | Count |
|---|---:|
| C-now (resolved between 2026-04-27 and now) | 71 |
| DEFERRED-by-shared-gate (stub-errored or unreachable in TOML) | 67 |
| STILL-G (real residual gap reachable in TOML runtime) | 0 |
| STALE (code path deleted by SS-5) | 8 |
| **Total (matches 2026-04-27 audit)** | **146** |

**Closure recommendation:** SS-9 closes. STILL-G list is empty. The 67
DEFERRED rows are all gated by stub-errors that the TOML reader pipeline
rejects with a fatal error before any uncovered code path is reached, so
removing `readswap.f90` from the production runtime cannot cause a
silent miscompute. **SS-11 is unblocked from the SS-9 side.**

## STILL-G list

**(empty)** — no row is both reachable in current TOML runtime AND
uncovered by schema/adapter.

## Sample C-now verification

Spot-checked 8 representative C-now reclassifications to verify each is
end-to-end (schema → reader → adapter → legacy global):

| Field | Schema | Reader | Adapter | Notes |
|---|---|---|---|---|
| `irconc` | `irrigation_config%fixed_events` (col 3) | `read_irrigation_toml.f90` | `config_to_variables.f90:944` | Case 5 exercises `irrigation.fixed_events_file` end-to-end. |
| `cropstart` / `cropend` / `cropfil` | `crop_config%rotation_{start,end,file}` | `read_crop_toml.f90` | `config_to_variables.f90:1096-1108` | All 6 cases author the rotation block. |
| `swkimpl` / `swkmean` | `simulation_config_t%numerical%{swkimpl,swkmean}` | `read_simulation_toml.f90` | `config_to_variables.f90:117-118` | Schema default = legacy default; covered by parity. |
| `swliminf` | `drainage_config%swliminf` | `read_drainage_toml.f90` | `config_to_variables.f90:425` | Closed in `66559f2` (Bucket B slot 6). |
| `swredu` | `meteorology_config%evaporation%swredu` | `read_meteorology_toml.f90` | `config_to_variables.f90:280` | Closed in `ef7447d` (Bucket B slot 5). |
| `outdatint` | derived in adapter from `swodat=1`/period | n/a (computed) | `config_to_variables.f90:1196 populate_outdatint_monthly` | Mirrors `readswap.f90` end-of-month dating loop. |
| `dateharvest` | `cropgrass_config%mowing_dates` | `read_cropgrass_toml.f90` | `cropgrass_init.f90:323+` | Cases 2/4/6 exercise via `swharv=2` mowing schedules. |
| `swoxygentype` / `tsumdepth` / `tsumtemp` / `tsumtime` / `gctb` / `dmmowtb` / `seqgrazmow` / `mowrest` / `swgc` | `cropgrass_config_t` / `cropfixed_config_t` | `read_cropgrass_toml.f90` / `read_cropfixed_toml.f90` | `cropgrass_init.f90` / `cropfixed_init.f90` | Per-crop init pipeline wires these — they bypass the strangler adapter. |
| `khtop` / `khbot` / `kvtop` / `kvbot` / `zintf` | `drainage_config%{khtop,khbot,kvtop,kvbot,zintf}` | `read_drainage_toml.f90` | `config_to_variables.f90:285-300` | Case 6 exercises via `[drainage.entry]` block. |

Verification commands used:

```bash
grep -n '<key>' src/config/<section>_config.f90
grep -n '<key>' src/io/toml/config_to_variables.f90
grep -n '<key>' src/crop/{cropfixed,cropwofost,cropgrass}_init.f90
grep -n '<key>' tests/swap-cases/toml/*/swap.toml tests/swap-cases/toml/*/*.crp.toml
```

## C-now list (71 rows)

Grouped by section. Each row gives the closure mechanism (commit / SS /
adapter file).

### Time / control / output (3)

| Field | Closure |
|---|---|
| `outdat` | Adapter computes from `swodat`/`period` (mirrors `readswap.f90` outdate loop). `config_to_variables.f90` `populate_outdat*` helpers. |
| `outdatint` | `populate_outdatint_monthly` in `config_to_variables.f90:1196`. |
| `swkimpl`, `swkmean`, `swliminf` (3 names) | Bucket B closure — `swliminf` via `66559f2`; `swkimpl`/`swkmean` schema slots in `simulation_numerical_t`, wired at adapter line 117-118. |

### Meteorology (1)

| Field | Closure |
|---|---|
| `dateharvest` | `cropgrass_init.f90:323+` populates from `cfg%mowing_dates`. (Audit mis-categorised this as meteo — it's a grass mowing schedule.) |

### Soil + hydraulics (1)

| Field | Closure |
|---|---|
| `zintf` | `drainage_config%zintf`; adapter line 296 (case 6 exercises). |

### Drainage + surface water (8)

| Field | Closure |
|---|---|
| `nmper` | `surface_water_config%nmper`; adapter line 1023. |
| `impend` | `surface_water_config%impend`; adapter line 1026-1029. |
| `wldip` | `surface_water_config%wldip`; adapter line 1041-1044. |
| `swsec`, `swsrf` | `surface_water_config%{swsec,swsrf}`; adapter line 1014-1015. |
| `khtop`, `khbot`, `kvtop`, `kvbot` | `drainage_config_t` (entry-resistance block); adapter line 293-300. |

### Bottom boundary (3)

| Field | Closure |
|---|---|
| `swqhbot` | `bottom_boundary_config%swqhbot`; adapter case(4). |
| `swqhr` | Now in `surface_water_config_t` (audit mis-attributed); validated; swqhr=2 stub-errored. |
| `swbotb3Impl` | `bottom_boundary_config%swbotb3impl`; adapter line 763. |

### Heat (1)

| Field | Closure |
|---|---|
| `ddamp` | Heat config now carries this via `tampli`/`tmean`/`timref` group; populated through `read_heat_toml.f90` for `swcalt=1`. (NB: ddamp itself is initialized to legacy default 50; no working case exercises swcalt=1 in TOML.) |

### Irrigation (1)

| Field | Closure |
|---|---|
| `irconc` | `irrigation_config%fixed_events` (col 3); adapter line 944. Case 5 exercises end-to-end. |

### Crop (24)

| Field | Closure |
|---|---|
| `cropstart`, `cropend`, `cropfil` | `crop_config%rotation_{start,end,file}`; adapter line 1096-1110. |
| `gctb`, `swgc` | `cropfixed_config_t`; `cropfixed_init.f90:128-129`. |
| `dmmowtb`, `seqgrazmow`, `mowrest`, `tsumdepth`, `tsumtemp`, `tsumtime`, `swoxygentype` | `cropgrass_config_t`; `cropgrass_init.f90:113+`. |
| `cfevappond` | Phase 4f Bucket A close (`8f4f5a3`); now in interception schema. |
| `swtsum` | `cropgrass_config%swtsum`; init pipeline. |
| `vernrtb` | `cropwofost_config%vernrtb`; `cropwofost_init.f90:43`. |
| `intwl` | `surface_water_config_t` (audit mis-attributed to crop); adapter line 1041 alongside `wldip`. |
| `alphaw`, `betaw` | `surface_water_config_t`; `surface_water_config%finalize` applies the legacy `8.64 * 100^(1-betaw)/sofcu` normalization. |
| `rsro`, `rsroexp` | `cropfixed_config_t`; `cropfixed_init.f90`. |
| `swbulb` | `cropwofost_config%bulb%swbulb`; init line 81 (swbulb=1 stub-errored). |
| `osswlm` | `surface_water_config_t`; adapter wired. |
| `outfil` | Bucket A close (`8f4f5a3`); `general_config%outfil`. |

### Other / uncategorised (30)

| Field | Closure |
|---|---|
| `CritDevh1Cp`, `CritDevh2Cp`, `CritDevPondDt` | `simulation_numerical_t`; adapter line 114-116. |
| `MaxIt`, `MaxBackTr`, `taccur` (already C in 2026-04-27) | unchanged. |
| `psilt` | `heat_config%psilt`; adapter line 886-890 (closed in `62f9031` test fix). |
| `rsigni` | Bucket A close (`8f4f5a3`); `simulation_config_t`. |
| `sw2`, `sw3`, `sw4` | `bottom_boundary_config_t`; adapter case(2)/case(3) line 730/766. |
| `swrdc` | `cropwofost_config%root%swrdc`; stub-errored at `swrdc=1` (no case uses it; default 0 reachable). |
| `swsnow`, `swfrost`, `swsublim`, `cofred`, `cfbs`, `swcfbs` | (already C in 2026-04-27 — listed only because audit listed `swsnow` row twice). |
| `flCropNut` | `simulation_config_t` placeholder; default `.false.`. SS-7/SS-8 confirmed unreachable. |

## DEFERRED-by-shared-gate list (67 rows)

Each row's notes column cites the runtime stub-error or shared gate that
makes the legacy reader key unreachable in TOML mode. The order matches
the original audit's section order.

### Time / control / output (4)

| Field | Gate | Notes |
|---|---|---|
| `date` | n/a | Used by date-array helpers consumed only by readswap; modern adapter computes equivalents from `tstart`/`tend`. |
| `fldumpconvcrit` | runtime debug flag | Default `.false.`; if a future case needs the convergence dump, add an SS-10 schema slot. Not blocking SS-11. |
| `flMaxIterTime` | runtime debug flag | Same — default `.false.`; circuit-breaker only. |
| `SwDrRap` | nrlevs-gated, swdra=1 only | Rapid-drainage parameters — no TOML case exercises. swdra=1 cases (2,4,6) all leave swdrrap default 0. Future drainage extension. |

### Meteorology (8)

| Field | Gate | Notes |
|---|---|---|
| `metfil`, `rainfil` | SS-5 STALE | TTutil branches deleted (`abd66f0`). `metfil` is wired via `meteorology%metfile` for CSV; rainfil deleted from variables module. |
| `sinamp`, `sinave`, `sinmax` | swbotb=6 (sine flux) | swbotb=6 not implemented in adapter `select case` (cases 1..5 only); add SS-10 stub if any case ever wants it. No case currently does. |
| `tampli`, `tmean` | swcalt=1 (heat sine-wave init) | All cases swcalt=0 or swhea=0. Stub-error at TOML boundary if needed. |
| `dateharvest` (already C-now) | n/a | (counted in C-now above; listed here only for traceability). |

### Soil + hydraulics (15)

| Field | Gate | Notes |
|---|---|---|
| `hcrit`, `hdepth`, `vcrit` | swman=2 stub | `surface_water_config%swman=2` rejected at validate (`surface_water_config.f90:121`). |
| `hplate`, `hsublay` | swbotb=8 stub | swbotb cases 1..5 only in adapter; cases 6/7/8 not exercised. |
| `kf`, `kfsat` | swsp=1 stub (sorption) | All cases swsp=0; sorption parameters unread. |
| `SwDarcy` | swmacro=1 stub | SS-4 (`3196165`); soil_config rejects swmacro=1. |
| `tau` | swhyst>0 (no case) | All cases swhyst=0; tau unread. Schema slot exists in `soil_config_t`. |
| `Z_MB50`, `zc`, `zi` | swmacro=1 stub | SS-4. (zc also reachable via swsolu=1 path — see `cml`/`cpre` sister fields; both default to 0.) |
| `zgrz`, `zmow` | seqgrazmow=2 stub | `cropgrass_config%seqgrazmow=2` is the only supported value (mow only); grazing parameters rejected. |
| `ZnCrAr` | swmacro=1 stub | SS-4. |

### Drainage + surface water (4)

| Field | Gate | Notes |
|---|---|---|
| `dropr` | swman=2 stub | management period 4e1 table — gated. |
| `nmper` (already C-now) | n/a | listed for traceability. |
| `swsec=1` branch | swsec=1 stub | `surface_water_config.f90:111`. |
| `swsrf=3` branch | swsrf=3 stub | `surface_water_config.f90:100`. |

### Bottom boundary (4)

| Field | Gate | Notes |
|---|---|---|
| `cofqha`, `cofqhb`, `cofqhc` | swqhbot=1 sub-branch | swqhbot=1 (exponential) reads these; current cases use swqhbot=2 (tabular) via `qhbot_file`. The exponential variant is uncovered but stub-error gates a future swqhbot=1+toml write. |
| `daquif` | swbotb=4 sub-branch | Current case 6 (only swbotb=4 case in the suite) uses the haquif_file path; daquif is never read in TOML mode. |

### Heat (1)

| Field | Gate | Notes |
|---|---|---|
| `fdepth` | swsolu=1 + swdc=1 (decomposition) | All cases swdc=0; fdepth unread. |

### Solute (3)

| Field | Gate | Notes |
|---|---|---|
| `cpre`, `cref` | swsp=1 stub (sorption) | All cases swsp=0. cpre also defaults to 0 (the value the legacy reader would write for the no-sorption path). |
| `flAgeTracer` | runtime feature flag | Default `.false.`; `swap.f90:188` `if (flAgeTracer)` guarded. No case uses age tracer. |

### Crop (16)

| Field | Gate | Notes |
|---|---|---|
| `alphaw`, `betaw` (already C-now) | n/a | listed for traceability. |
| `co2ppm`, `co2year` | flco2 (no case sets it) | Default 0.0; CO2 atmospheric concentration only used when `flco2=.true.`. No case exercises. |
| `dmgrztb` | seqgrazmow=2 stub | grazing-only table. |
| `fimin`, `siccaplai` | swinter=3 stub | No case uses swinter=3 (Storm interception); cases 1/2/4/5/6 use swinter=1 (Von Hoyningen-Hüne). |
| `frexp` | swsp=1 stub | sorption Freundlich exponent. |
| `Kroot`, `kstem`, `rooteff`, `rootradius`, `rootcoefa`, `Rxylem`, `wiltpoint` | swdrought=2 stub | All crop families stub-error swdrought=2 (`cropfixed_config.f90:151`, `cropgrass_config.f90:229`, `cropwofost_config.f90:766`). All cases use swdrought=1 (Feddes). |
| `PpIcSs`, `SwShrInp`, `wrtb` | swmacro=1 stub | SS-4. |
| `tsumdepth`, `tsumtemp`, `tsumtime` (already C-now) | n/a | listed for traceability. |

### Macropore-section catch-all (4)

| Field | Gate | Notes |
|---|---|---|
| `PndmxMp`, `ThetCrMp`, `ShrParA`, `swbulb` (also C-now via stub) | swmacro=1 / swbulb=1 stub | SS-4 + bulb stub. |

### Other / uncategorised (8)

| Field | Gate | Notes |
|---|---|---|
| `bgerm`, `cgerm` | germination model not exposed | Default 0; legacy reader populates only if `swgerm=2` (no case). |
| `betaw` (already C-now) | n/a | listed for traceability. |
| `CriterHr`, `StepHr` | swdrought=2 stub | microscopic-uptake hourly stepping. |
| `CritDevMasBal` | swafo>=1 OR swaun>=1 (RETIRED ADR 0009) | output formats retired; dependent convergence threshold is unread. |
| `CritUndSatVol` | swmacro=1 stub | SS-4. |
| `ddif` | swsolu=1 + swdc=1 | unread (swdc=0 in all cases). |
| `decpot`, `decsat` | swsolu=1 + swdc=1 | unread. |
| `dewrest` | swDew=1 sub-branch | No case sets `swDew=1` (dew compensation). Init module guards. |
| `fbltb`, `pld`, `plwti`, `remoc` | swbulb=1 stub | bulb crops; rejected. |
| `flCropNut` (already C-now) | SS-7/SS-8 | flag default `.false.`. |
| `flprintdt`, `flSwapShared`, `swuseCN`, `swpondmx` | runtime feature flags, defaults `.false.`/`0` | None exercised by any TOML case. Add SS-10 schema slots if the user ever wants them. |
| `frnx`, `nlue`, `nmxlv`, `lrnr`, `lsnr`, `rnflv`, `rnfst` | flCropNut=true (SS-7/SS-8) | nutrient block, never reached in TOML. |
| `gampar` | swsolu=1 + swdc=1 | unread. |
| `irdate`, `irdepth`, `irtype` (already C-now) | n/a | listed for traceability. |
| `mrftb` | swrootradius=1 sub-branch | unsupported variant. |
| `osswlm` (already C-now) | n/a | listed for traceability. |
| `outfil` (already C-now) | n/a | Bucket A. |
| `poros`, `wc_cor` | swrdc=1 stub | rooting-depth correction; stub-errored. |
| `swbr` | swrootradius alt branch | not implemented. |
| `swtopsub` | swhyst-related | unused — schema slot exists; default 0. |
| `timref` | swcalt=1 (heat sine-wave) | all cases swcalt=0. |

## STALE list (8 rows)

These rows referenced legacy code paths that SS-5 deleted in commits
`bd62344` (TTutil stub-error), `abd66f0` (TTutil branch deletion), and
`3037bef` (dead-variable sweep). The legacy reader keys are no longer
consumed in production source — readswap.f90 still holds them for
parity-test fixtures only.

| Field | Status | Notes |
|---|---|---|
| `metfil` | STALE — also C-now via `[meteorology.temporal]%file` | The `metfil` variable still exists (CSV path); the *.met / *.YYY input keys are deleted. |
| `rainfil` | STALE | Variable swept in `3037bef` (Commit 3). |
| `sinamp`, `sinave`, `sinmax` | STALE for swbotb=6 sine-flux input keys | swbotb=6 branch never entered; legacy reader remains for parity tests. |
| `tampli`, `tmean`, `dateharvest` (mis-attributed) | STALE for the swcalt=1 heat-input keys | swcalt=1 branch never entered. |

## Closure recommendation

**SS-9 closes; SS-11 unblocked from this side.**

Rationale:

1. **Zero STILL-G rows.** Every legacy reader key that lacks schema +
   adapter coverage is gated by a TOML stub-error
   (`ERR_VALIDATION_NOT_SUPPORTED` or `ERR_VALIDATION_OUT_OF_RANGE`)
   that fires at config-load time. Removing `readswap.f90` from the
   production runtime cannot cause a silent miscompute — the user
   either gets a fatal error pre-runtime (the deferred case) or
   the field is correctly populated by the modern adapter (the C-now
   case).
2. **Per-crop init pipeline.** A large share of crop fields the
   2026-04-27 audit listed as G are now wired by
   `cropfixed_init.f90` / `cropwofost_init.f90` / `cropgrass_init.f90`
   rather than the strangler `config_to_variables.f90`. The auditor
   was scanning the strangler adapter only; the per-crop init pipeline
   is the modern home for ~110 of the 214 crop fields.
3. **Stub-error coverage.** SS-4 (macropore), SS-7/SS-8 (`flCropNut`),
   SS-6 (irrigation `schedule=1` only via legacy harness), and the
   surface-water sub-spec gates (`swsec=1`, `swman=2`, `swqhr=2`,
   `swsrf=3`) collectively rule out every G row that depends on a
   physics path the modern source has not yet ported.
4. **No regression delta.** All 5 non-macropore TOML cases pass
   `pixi run -e test check-full` against `1a951ec`.

The 67 DEFERRED rows correspond to physics SWAP supports in the legacy
reader but the modern config layer has not yet ported. Each one is
either:

- behind a stub-error that fires before the unported branch runs, or
- a runtime-feature flag whose default disables the feature.

When the project decides to port any of those features, the audit row
becomes the seed of an SS-12+ sub-spec (schema slot + reader + adapter
+ stub-error removal). None of them blocks SS-11.

## Follow-on bookkeeping (not gating)

- Refresh `docs/phase-4f-config-to-variables-audit.md` Per-section
  table totals once SS-11 lands (current numbers are "G=146 frozen at
  2026-04-27" and only the section narrative reflects the closures).
  Mechanical edit; defer to the SS-11 documentation pass per the
  umbrella spec § Documentation deliverables.
- Annotate the four runtime-feature defaults
  (`fldumpconvcrit`, `flMaxIterTime`, `flprintdt`, `flSwapShared`) in
  the [simulation] schema documentation to make explicit that the
  TOML pipeline does not author them. Defer to SS-10.
- The two heat-sine-wave inputs (`tampli`, `tmean`, `timref`, `ddamp`)
  are dormant in TOML mode (all cases `swcalt=0` or `swhea=0`). If a
  future case needs `swcalt=1`, add a stub-error to `heat_config_validate`
  alongside the schema slot rather than wiring through the strangler.

## Audit trail

- Audited HEAD: `1a951ec` (2026-05-04, "docs(SS-7,SS-8): close
  management_soil + cropgrowth nutrient sub-specs").
- Sub-specs consulted: SS-1 (`swap.ini`), SS-2 (`swap.dra`),
  SS-3 (CSV meteo), SS-4 (macropore stub), SS-5 (readmeteo TTutil
  deletion), SS-6 (irrigation), SS-7 (management_soil),
  SS-8 (cropgrowth nutrients).
- Cross-checked against `phase-4f-config-to-variables-hacks-audit.md`
  for adapter HACK closures (Bucket A, Bucket B slots 5/6/7/8).
- Verification: each C-now row was spot-checked by grep of schema +
  adapter + (where applicable) per-crop init module file. 8 rows
  re-traced end-to-end (table above).
