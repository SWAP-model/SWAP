# Subsystem Migration Discovery: TimeControl

**Date:** 2026-05-12
**Status:** discovery (read-only inventory)
**Migration #:** 10 of N — cross-cutting infrastructure arc following tillage
**Branch:** `development`

**Scope note:** This arc carves the **TimeControl runtime state** (clock, calendar
counters, output-schedule counters, runtime-evaluated flags) out of legacy
`variables.f90` (plus the pre-existing `tc_*` synchronisable SAVE-local migration block)
into a new `timecontrol_state_t`. Static simulation-clock configuration (`tstart`,
`tend`, `dtmin`, `dtmax`, `nprintday`, `period`, `swres`, `swodat`, `outdat(:)`,
`outdatint(:)`, etc.) is already in `simulation_config_t` and **stays as legacy**
(config-constants — tillage lesson #1). The cross-subsystem `flzerointr` / `flzerocumu`
gate flags are read from ~25 files and are intentionally **deferred**: they are the
gates that drive `*_state_t%intr%reset()` / `cumu%reset()` and should keep a single
global truth source until a future "subsystem-reset orchestration" arc converts them
into return values from `state%timecontrol%advance()`. See Section 7 H-1 and Section 10
OQ-1.

**Predecessor docs:**
- `docs/superpowers/specs/state-migration-playbook.md` (29 lessons — config-constant
  vs runtime-state distinction is the most relevant; Strategy B is master)
- `docs/superpowers/specs/2026-05-12-state-migration-tillage-discovery.md` (style + 6-task precedent)
- `docs/superpowers/specs/2026-05-11-state-migration-atmosphere-discovery.md` (cross-cutting precedent)
- `docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-discovery.md` (largest-scope precedent)
- `docs/adr/0039-state-migration-tillage.md`
- `docs/adr/0037-state-migration-atmosphere.md`

> Read-only discovery: no code changes made. All file:line references are anchors for
> the design phase.

---

## 1. Big picture

### Subsystem role

`timecontrol.f90` is the central time-stepping driver: it advances the simulation
clock (`t1900`, `tcum`, `t`), maintains the day/year/month/cumulative counters
(`daynr`, `daycum`, `iyear`, `imonth`), drives the variable-step Newton-Raphson size
(`dt`, `dtprevious`, `dtEvent`, `tEvent`) and emits the cross-subsystem gate flags
(`flDayStart`, `flDayEnd`, `flRunEnd`, `flYearStart`, `flZeroIntr`, `flZeroCumu`,
`flOutput`, `flBalOutput`, …). Every other subsystem reads these flags and counters.

The companion routine `IterTime(task)` (lines 639–697) is a CPU-time-budget watchdog
backed by `tc_tmptimestart` / `tc_tmptimeend`.

### Lines of code

```
698  src/core/timecontrol.f90   (TimeControl 1..636 + IterTime 639..697)
```

Single home file. Mid-sized — slightly larger than tillage (446), smaller than
boundtop+boundbottom (494) and atmosphere (1 398).

### Subroutine structure

| Subroutine | Lines | Role |
|---|---|---|
| `TimeControl(task, state)` | 4–636 | Dispatcher — case 1 init, 2 advance, 3 dt-reduce, 9 end-of-day SSDI hook |
| `IterTime(task)` | 639–697 | CPU-time watchdog (1=start, 2=intermediate, 3=stats) |

The `associate` block at lines 53–64 aliases the legacy `tc_*` "synchronisable
SAVE-locals" back to their pre-migration names (`datea => tc_datea`,
`dtEvent => tc_dtEvent`, …). The migration in this arc REPLACES that associate
block with `state%timecontrol%*` references.

### Entry points from outside the home file

| # | Call site | File:line | Phase | Task arg |
|---|---|---|---|---|
| 1 | `TimeControl(1, state)` | `swap.f90:191` | init | task=1 |
| 2 | `TimeControl(3, state)` | `swap.f90:376` | per-step (dt-reduction) | task=3 |
| 3 | `TimeControl(2, state)` | `swap.f90:395` | per-step (advance) | task=2 |
| 4 | `TimeControl(9, state)` | `swap.f90:430` | end-of-day (SSDI hook) | task=9 |
| 5 | `TimeControl(1, state)` | `swap.f90:598` | DLL exchange re-init | task=1 |
| 6 | `TimeControl(3, state_om)` | `swapoutput.f90:4018` | mini-sim (orgmat balance) | task=3 |

Call site #6 (mini-sim) uses a private `state_om` (not the main `state`) — same
pattern as SS-SWC S-2.12B mini-sim split. The new `state%timecontrol` participates
automatically because `state_om` is a `swap_state_t` clone.

### Signature status

```bash
grep -n "type(swap_state_t)" src/core/timecontrol.f90
# 14:      use swap_state_mod, only: swap_state_t            ! [SS-SWC S-2.12B]
# 18:      type(swap_state_t), intent(in) :: state
```

`TimeControl(task, state)` ALREADY takes `state` (added by SS-SWC S-2.12B for the
SSDI hook on case 9). Today the arg is declared `intent(in)` because TimeControl
itself does not write to state — it reads state only to forward to `SSDI_irrigation(9)`.

**This arc needs:** promote `intent(in)` → `intent(inout)` so that TimeControl can
write into `state%timecontrol` (formerly assigning legacy globals). `IterTime` will
also need a `state` argument added — it currently has no state plumbing.

---

## 2. Owned-globals inventory

The TimeControl-owned legacy globals split into FOUR groups by mutability and
cadence. Following the tillage-lesson distinction: only the runtime-mutated fields
migrate (Groups C+D); the static config copies (Group A) stay as legacy globals,
and the cross-subsystem reset gate flags (Group B) are deferred for a future
subsystem-reset-orchestration arc.

### Group A — Static simulation config (read-only after `config_to_variables`)

These have an existing typed home in `simulation_config_t` (`src/config/simulation_config.f90`)
and are populated once by `config_to_variables.f90:79–135`. They are NOT mutated at
runtime. Per ADR 0039 lesson #1 they stay as legacy globals until a future
config-consolidation arc.

| Variable | Type | Shape | Source | Notes |
|---|---|---|---|---|
| `tstart` | real(8) | scalar | `config%simulation%tstart` | sim start (days since 1900) |
| `tend` | real(8) | scalar | `config%simulation%tend` | sim end |
| `dtmin` | real(8) | scalar | `config%simulation%numerical%dtmin` | min step length |
| `dtmax` | real(8) | scalar | `config%simulation%numerical%dtmax` | max step length |
| `nprintday` | integer | scalar | `config%simulation%nprintday` | outputs per day |
| `period` | integer | scalar | `config%simulation%period` | interm. period length |
| `swres` | integer | scalar | `config%simulation%swres` | reset counter at year |
| `swodat` | integer | scalar | `config%simulation%swodat` | extra output dates? |
| `swheader` | integer | scalar | `config%output…` | print header switch |
| `swscre` | integer | scalar | `config%general%swscre` | screen verbosity |
| `outdat(maout)` | real(8) | array | adapter | balance output dates |
| `outdatint(maout)` | real(8) | array | adapter (incl. monthly auto-fill) | intermediate output dates |
| `flprintdt` | logical | scalar | derived from `nprintday`/period config | (Initialized once in `config_to_variables`; never mutated.) |
| `msteps` | integer | scalar | `config%simulation%numerical%msteps` | iteration cap |
| `MaxIt` | integer | scalar | `config%simulation%numerical%MaxIt` | NR iter cap |
| `MaxIterTime` | integer | scalar | (config TBD) | CPU-time cap |
| `nmetdetail`, `swmetdetail`, `swrain`, `swetsine` | mixed | scalars | `config%meteo.*` | gates that drive `flmetdetail`/`flrainintens`/`flmeteodt`/`fletsine` derived flags |
| `swdra`, `swhea`, `swsnow`, `swsolu`, `swirfix` | integer | scalars | per-subsystem configs | drive `flDrain`/`flSurfaceWater`/`flTemperature`/`flSnow`/`flSolute`/`flIrrigate` |

Subtotal: ~18 config-constants kept as legacy globals (out of scope for THIS arc).

### Group B — Cross-subsystem reset-gate flags (DEFERRED)

| Variable | Type | Shape | Cadence | Readers (file count) |
|---|---|---|---|---|
| `flZeroIntr` | logical | scalar | flips on/off per interm-period boundary | 12 files (drainage, surfacewater, atmosphere, soilwater, solute, snow, meteoday, soilhydraulics, waterbalance, initialize) |
| `flZeroCumu` | logical | scalar | flips on/off per cumu-output boundary | 11 files (similar set) |

These two flags are the canonical "reset gate" signal read by every subsystem's
`*_state_t%intr%reset()` / `cumu%reset()` call site. Migrating them into
`state%timecontrol` would force every subsystem to read `state%timecontrol%flZeroIntr`
instead of a bare global. That is a mechanical retarget BUT it cascades into:
- 14 read sites for `flzerointr`
- 14 read sites for `flzerocumu`

across 25+ files (including the state modules themselves). **Recommendation:** defer
to a follow-on arc (or fold into the future "subsystem-reset orchestration" arc
where `flzero*` become return values from `state%timecontrol%advance()`). See
Section 10 OQ-1.

### Group C — Runtime clock + calendar state (MIGRATE)

These are mutated each timestep / each day and have no static config home. They
are the canonical runtime-state to migrate.

| Variable | Type | Shape | Cadence | Activity gate | Representative write |
|---|---|---|---|---|---|
| `t1900` | real(8) | scalar | per-step | always | `timecontrol.f90:138, 280` |
| `t` | real(8) | scalar | per-step | always | `timecontrol.f90:147, 278, 338` |
| `tcum` | real(8) | scalar | per-step | always | `timecontrol.f90:149, 279` |
| `daynr` | integer | scalar | per-day | always | `timecontrol.f90:150, 329, 337` |
| `daycum` | integer | scalar | per-day | always | `timecontrol.f90:151, 330` |
| `iyear` | integer | scalar | per-day (sometimes per-year) | always | `timecontrol.f90:156, 322`; also adapter `config_to_variables.f90:97` and DLL `swap.f90:597` |
| `iyearm1` | integer | scalar | per-day | always | `timecontrol.f90:321` |
| `imonth` | integer | scalar | per-day | always | `timecontrol.f90:157, 323`; adapter `config_to_variables.f90:98` |
| `daymeteo` | integer | scalar | per-day | always | `timecontrol.f90:152, 559, 565` |
| `yearmeteo` | integer | scalar | per-year | always | `timecontrol.f90:158, 561` |
| `date` | character(11) | scalar | per-day | always | `timecontrol.f90:161, 167, 169, 326` (via `dtdpst`) |
| `timjan1` | real(8) | scalar | per-year | always | `timecontrol.f90:146` |
| `tc_datea(6)` | integer | array | per-call | always | `timecontrol.f90:141, 155, 320, 555` — already in tc_ namespace |
| `tc_fsec` | real(4) | scalar | per-call | always | `timecontrol.f90:141, 145, 155, 320, 555` |
| `tc_nextyear` | integer | scalar | per-day | always | `timecontrol.f90:557` |
| `swmeteo` | integer | scalar | per-year (mutated when meteo year flips) | always | `timecontrol.f90:196, 200, 567, 571` |

Subtotal: 16 runtime clock+calendar fields.

### Group D — Timestep / output-schedule machinery (MIGRATE)

| Variable | Type | Shape | Cadence | Activity gate | Representative write |
|---|---|---|---|---|---|
| `dt` | real(8) | scalar | per-step (heaviest mutator) | always | `timecontrol.f90:234, 237, 242, 251, 405, 414, 415, 422, 429, 435, 441, 442, 591, 608, 618, 626` (also config seed `config_to_variables.f90:121, 576`) |
| `dtold` | real(8) | scalar | per-step | always | `timecontrol.f90:256` |
| `tc_dtEvent` | real(8) | scalar | per-event | always | `timecontrol.f90:206, 211, 217, 223, 226, 365, 370, 372, 380, 392, 393, 397` |
| `tc_tEvent` | real(8) | scalar | per-event | always | `timecontrol.f90:148, 226, 397` |
| `tc_dtprevious` | real(8) | scalar | per-step | always | `timecontrol.f90:239, 405, 409, 412, 416, 579, 600, 612, 619` |
| `tc_tchange` | real(8) | scalar | per-meteo-period | flmetdetail | `timecontrol.f90:379, 474` |
| `tc_tcumold` | real(8) | scalar | per-output-event | flprintshort | `timecontrol.f90:134, 453, 463` |
| `tc_flprevious` | integer | scalar | per-step | always | `timecontrol.f90:241, 245, 404, 411, 419, 423, 599, 611` |
| `tc_flTnext` | logical | scalar | per-step | always | `timecontrol.f90:243, 246, 352, 398, 424, 601` |
| `isteps` | integer | scalar | per-step (advances + resets per-day) | always | `timecontrol.f90:129, 265, 552` |
| `daycrop` | integer | scalar | per-day | flCropCalendar | (write commented in TC; daycrop is crop-owned — moved to crop arc already) |
| `nprintcount` | integer | scalar | per-day | flprintshort | `timecontrol.f90:135, 467` |
| `cntper` | integer | scalar | per-day | always | `timecontrol.f90:132, 331, 341, 520` |
| `outper` | real(8) | scalar | per-output-event | always | `timecontrol.f90:133, 302, 452, 462, 501` |
| `ioutdat` | integer | scalar | per-bal-output | always | `timecontrol.f90:130, 514, 544` |
| `ioutdatint` | integer | scalar | per-interm-output | flprintdt | `timecontrol.f90:131, 529, 535` |
| `rainrec` | integer | scalar | per-rain-event | swrain>0 | `timecontrol.f90:389` |
| `wrecord` | integer | scalar | per-meteo-record | flmetdetail | `timecontrol.f90:99, 353` |
| `dtEventRain` | real(8) | scalar | per-rain-event | swrain>0 | `meteodt.f90:345` (atmosphere writes; tc reads at 392) — **co-owned with atmosphere/meteodt** |
| `metperiod` | real(8) | scalar | per-year-init | swmetdetail | `timecontrol.f90:98` (derived from `nmetdetail`) |

Subtotal: ~19 timestep / schedule fields (counting `metperiod` as runtime-derived).

### Group E — Runtime evaluated boolean flags (MIGRATE)

These are mutated within TimeControl (or by close partners) per-step or per-day.
They differ from Group B (flzero*) which are cross-subsystem reset gates.

| Variable | Type | Cadence | Reset by | Notes |
|---|---|---|---|---|
| `flDayStart` | logical | per-step | TimeControl | flips at first step of each day |
| `flDayEnd` | logical | per-step | TimeControl | flips at last step of each day |
| `flRunEnd` | logical | per-step | TimeControl | flips when t1900 reaches tend |
| `flYearStart` | logical | per-day | TimeControl | flips at iyear change |
| `flprintshort` | logical | per-day | TimeControl init | derived from `nprintday>1` or `flprintdt` |
| `floutputshort` | logical | per-step | TimeControl | derived from `flprintshort` schedule |
| `floutput` | logical | per-day | TimeControl | output-day signal |
| `flbaloutput` | logical | per-day | TimeControl | bal-output-day signal |
| `flheader` | logical | per-day | TimeControl + swapoutput.f90 | swapoutput.f90:2299 also resets (co-writer) |
| `flheadirg` | logical | per-day | TimeControl + swapoutput.f90 | swapoutput.f90:2310 also resets (co-writer) |
| `flIrg1Start` | logical | per-day | TimeControl + swapoutput.f90 | swapoutput.f90:2299 also resets |
| `flUpdMetDet` | logical | per-step | TimeControl + meteodt.f90 | meteodt.f90:361 sets `.false.` after consuming |
| `fldecdt` | logical | per-step | timestep_control_mod (already module) | **already migrated** out of variables.f90 |
| `fldecdtmin` | logical | per-step | TimeControl + boundtop | boundtop sets `.true.`; TC resets `.false.` |
| `fldtmin` | logical | per-step | TimeControl | derived state from `dt` vs `dtmin` |
| `fldtreduce` | logical | per-step | swap.f90 only | swap.f90 main-loop local-pattern flag |
| `flmetdetail` | logical | init-once | TimeControl | derived from `swmetdetail` |
| `flmeteodt` | logical | init-once | TimeControl | derived from `flmetdetail` or `flrainintens` |
| `flrainintens` | logical | init-once | TimeControl | derived from `swrain>0` |
| `fletsine` | logical | init-once | TimeControl | derived from `swetsine` |
| `flSnow` | logical | init-once | TimeControl | derived from `swsnow` |
| `flDrain` | logical | init-once | TimeControl | derived from `swdra==1` |
| `flSurfaceWater` | logical | init-once | TimeControl | derived from `swdra==2` |
| `flTemperature` | logical | init-once | TimeControl | derived from `swhea` |
| `flSolute` | logical | init-once | TimeControl | derived from `swsolu` |
| `flIrrigate` | logical | init-once | TimeControl | derived from `swirfix` |

Subtotal: ~26 runtime-evaluated booleans. The 11 "init-once" booleans
(`flmetdetail`, `flmeteodt`, `flrainintens`, `fletsine`, `flSnow`, `flDrain`,
`flSurfaceWater`, `flTemperature`, `flSolute`, `flIrrigate`, plus `flprintshort`)
are derived from config switches inside TimeControl(task=1) — they are
runtime-mutated under TimeControl ownership but never change after init. They are
on the boundary between config-constant and runtime-state; the migration could
either include them (own them in state%timecontrol since they ARE set by TC) or
leave them as legacy globals (since they could equivalently be config-derived).
Recommendation: **include** — they fit naturally with the per-day flags, and
they are written by TC which is the migration target.

### Roll-up totals (RUNTIME state to migrate)

| Group | Count | Notes |
|---|---|---|
| C — clock + calendar | 16 | t1900, t, tcum, daynr, daycum, iyear, iyearm1, imonth, daymeteo, yearmeteo, date, timjan1, tc_datea(6), tc_fsec, tc_nextyear, swmeteo |
| D — timestep / schedule | 19 | dt, dtold, tc_dtEvent, tc_tEvent, tc_dtprevious, tc_tchange, tc_tcumold, tc_flprevious, tc_flTnext, isteps, nprintcount, cntper, outper, ioutdat, ioutdatint, rainrec, wrecord, metperiod, [dtEventRain co-owned] |
| E — runtime booleans | 26 | flDayStart/End/RunEnd/YearStart, flprintshort, floutputshort, floutput, flbaloutput, flheader, flheadirg, flIrg1Start, flUpdMetDet, fldecdtmin, fldtmin, fldtreduce, flmetdetail, flmeteodt, flrainintens, fletsine, flSnow, flDrain, flSurfaceWater, flTemperature, flSolute, flIrrigate |
| **Total runtime state** | **~61** | (Group B `flzero*` deferred adds 2 more later) |

### Reset cadence summary

No `flzerocumu`/`flzerointr` cohort fits here — TimeControl OWNS those gates (it
emits them); it doesn't consume them. All Group C+D fields are either:
- per-step / per-day / per-event mutated (no shared reset gate)
- init-once (set in case 1, then constant)
- monotonic counter (advances forever, never resets — e.g. `daycum`, `t1900`,
  `daynr` resets per-year not per-period)

**Flat layout is appropriate** (no cohort sub-records). See Section 8.

---

## 3. External readers (5-category framework)

For each Group C+D+E owned field, files outside `src/core/timecontrol.f90` and
`src/state/` that **read** the field. Categories:
(1) Output, (2) Compute, (3) Working buffer, (4) Init-seed, (5) Call-site arg.

### 3.1 Pre-flight reader count by field

```
dt          : 264 read sites  (heaviest — read from every solver loop)
t1900       : 143 read sites  (time-of-day calculations, balance dates)
period      :  92 read sites  (period-aware logic)
outper      :  77 read sites  (per-period accumulation denominators)
daynr       :  62 read sites  (day-of-year logic — meteo, output, crop)
daycum      :  55 read sites  (sim-day logic — output, banner, mass-balance)
tstart      :  53 read sites  (config-constant; outside this arc scope)
tcum        :  38 read sites  (cumulative time — meteodt, surfacewater)
tend        :  36 read sites  (config-constant; outside scope)
flprintshort:  27 read sites  (output gating)
dtmin       :  21 read sites  (config-constant; outside scope)
flmetdetail :  16 read sites
flzerointr  :  14 read sites  (DEFERRED — Group B)
flzerocumu  :  14 read sites  (DEFERRED — Group B)
dtmax       :  14 read sites  (config-constant)
ioutdat     :  13 read sites
fldaystart  :  11 read sites
fldecdt     :  10 read sites  (already in timestep_control_mod)
iyear       :   9 read sites
imonth      :   5 read sites
floutput    :   4 read sites
flbaloutput :   3 read sites
```

### 3.2 Distinct external reader files (Groups C+D+E only)

Excluding `src/core/timecontrol.f90`, `src/core/variables.f90`, `src/core/swap.f90`,
and `src/state/` (which already takes state):

```
src/atmosphere/et.f90
src/atmosphere/meteoday.f90
src/atmosphere/meteodt.f90
src/atmosphere/snow.f90
src/boundary/boundbottom.f90
src/boundary/boundtop.f90
src/config/simulation_config.f90    (init-seed cat 4)
src/core/initialize.f90             (init-seed cat 4 — zero-init)
src/crop/cropgrass_init.f90
src/crop/cropgrowth.f90
src/crop/irrigation.f90
src/crop/management_soil.f90
src/crop/oxygenstress.f90
src/crop/tillage.f90
src/drainage/divdra.f90
src/drainage/drainage.f90
src/drainage/surfacewater.f90
src/heat/frozencond.f90
src/heat/temperature.f90
src/io/macroporeoutput.f90
src/io/readmeteo.f90
src/io/swap_csv_output.f90
src/io/swapoutput.f90
src/io/toml/read_simulation_toml.f90  (config — out of arc)
src/io/toml/write_swap_config.f90     (config — out of arc)
src/macropore/macropore.f90
src/macropore/macrorate.f90
src/soil/soilhydraulics.f90
src/soil/waterbalance.f90
src/solute/agetracer.f90
src/solute/solute.f90
src/state/atmosphere_state.f90        (state — already plumbed)
src/state/soilwater_state.f90         (state — already plumbed)
src/state/solute_state.f90            (state — already plumbed)
src/state/surfacewater_state.f90      (state — already plumbed)
```

**Distinct in-scope external reader files: ~27** after excluding state/config files.
This is by far the largest reader inventory of any arc (heat 11, boundary 14,
crop-uptake 7, atmosphere 9, soil-water-core 22, tillage 3).

### 3.3 Per-file read-site estimate (in-scope)

Cumulative across `daynr daycum t1900 imonth iyear isteps cntper ioutdat outper
nprintcount flprintshort floutputshort flmetdetail floutput flbaloutput dt
fldaystart fldecdt`:

| File | Estimated read sites | Category mix |
|---|---|---|
| `swapoutput.f90` | 60+ (heaviest — output banner, dates, headers, bal section) | Cat 1 |
| `swap_csv_output.f90` | 15 | Cat 1 |
| `waterbalance.f90` | 20 | Cat 1 + Cat 2 |
| `soilhydraulics.f90` | 18 (dt heavy, fldaystart) | Cat 2 |
| `meteoday.f90` | 15 (daynr/iyear/daymeteo/yearmeteo) | Cat 2 + 4 |
| `meteodt.f90` | 15 (t1900/dt/tcum/dtEventRain) | Cat 2 |
| `cropgrowth.f90` | 12 (daynr/iyear) | Cat 2 |
| `management_soil.f90` | 12 (daynr/daycum/t1900) | Cat 2 |
| `drainage.f90` + `divdra.f90` | 14 (period/dt/t1900) | Cat 2 |
| `surfacewater.f90` | 10 (tcum/t1900/period) | Cat 2 |
| `temperature.f90` | 8 (daynr/dt/t1900) | Cat 2 |
| `frozencond.f90` | 5 | Cat 2 |
| `solute.f90` + `agetracer.f90` | 12 (dt) | Cat 2 |
| `macropore.f90` + `macrorate.f90` | 10 (t1900/dt) | Cat 2 |
| `boundtop.f90` + `boundbottom.f90` | 8 (daycum/t1900/dt) | Cat 2 |
| `oxygenstress.f90` | 3 | Cat 2 |
| `tillage.f90` | 4 (t1900) | Cat 2 |
| `cropgrass_init.f90` | 3 (t1900) | Cat 4 |
| `irrigation.f90` | 6 (period/t1900) | Cat 2 |
| `readmeteo.f90` | 5 (daymeteo/yearmeteo) | Cat 2 |
| `et.f90` | 4 (dt/fldaystart) | Cat 2 |
| `snow.f90` | 5 (dt) | Cat 2 |
| `macroporeoutput.f90` | 4 (period/t1900) | Cat 1 |

**Total external read sites: ~270–290** for the runtime-state subset (Groups C+D+E
counter/clock/flag reads only — excluding `dt` reads, which adds another 200+).
The `dt` field alone contributes ~264 reads. **Order-of-magnitude estimate:
~500 external read sites** (heaviest of any arc by a factor of 3+).

### 3.4 Cat 3 (working-buffer / mini-sim) hazard

The `swapoutput.f90` mini-sim writeback at line 4018 already passes a cloned
`state_om` containing `state_om%timecontrol`, so the SS-SWC mini-sim plumbing
is reusable. No new Cat 3 surface beyond the existing one. Confirmed clean.

---

## 4. Co-writers

Files outside `src/core/timecontrol.f90` that WRITE any Group C+D+E owned field.

| Co-writer file | Fields written | Status |
|---|---|---|
| `src/core/initialize.f90` | Zero-init for daynr, daycum, imonth, ioutdat, ioutdatint, isteps, iyear, iyearm1, nprintcount, cntper, outper, t, t1900, tcum, dt, daymeteo, yearmeteo, rainrec, wrecord + 11 flag fields | bare `use variables`; redundant after `timecontrol_init` |
| `src/core/swap.f90:597` | `iyear = datea(1)` — DLL exchange path re-derives iyear from Tstart before re-init | bare `use variables`; rewrite via state |
| `src/io/toml/config_to_variables.f90:97, 98` | `iyear = datea_init(1)`, `imonth = datea_init(2)` — pre-init seed of iyear/imonth from `tstart` | pre-init pattern: keep legacy until `timecontrol_init` reads from a transient buffer |
| `src/io/toml/config_to_variables.f90:121, 576` | `dt = config%simulation%numerical%dt`; `dt = config%soil%initial%dt` (swinco=3) | pre-init pattern: transient `dt_init_buf` then drain into state at init time, or rely on Group A `dt_config` field |
| `src/atmosphere/meteodt.f90:361` | `flUpdMetDet = .false.` (after consuming meteodt update) | dual-write needed; `meteodt` already takes `state` (SS-ATM A-1.6) |
| `src/atmosphere/meteodt.f90:345` | `dtEventRain = …` | co-write to a TimeControl-read field; `dtEventRain` is technically atmosphere/meteo-owned (atmosphere arc deferred it). **Cross-ownership** — leave as legacy global for THIS arc; TimeControl reads it through legacy until atmosphere meteo refactor migrates it. |
| `src/soil/soilhydraulics.f90:810` | `fldecdt = .true.` | `fldecdt` is already in `timestep_control_mod` — out of arc |
| `src/io/swapoutput.f90:2299` | `flIrg1Start = .false.` | dual-write needed; SwapOutput already takes state |
| `src/io/swapoutput.f90:2310` | `flheadirg = .false.` | dual-write needed; SwapOutput already takes state |
| `src/io/swapoutput.f90:3999, 4012` | `fldecdt = .true.` (mini-sim) | already in timestep_control_mod — out of arc |
| `src/io/swapoutput.f90:442, 579` | `t1900 = t1900` | self-assignment (OpenMP private hint); no write — confirmed safe |
| `src/macropore/macropore.f90:1716`, `src/macropore/macrorate.f90:82` | `t1900 = t1900` | self-assignment idiom; no real write |

**Co-writer summary:**
- `initialize.f90` — zero-init duplicates of state init (will be dropped).
- `config_to_variables.f90` — pre-init seed of `iyear`, `imonth`, `dt`. Use the
  transient-buffer pattern (soil-water-core lesson #4).
- `swap.f90:597` — DLL re-init path; mechanical state retarget.
- `meteodt.f90:361` — `flUpdMetDet = .false.` after consume; mechanical dual-write.
- `swapoutput.f90:2299, 2310` — output-side flag flips; mechanical dual-write.
- `dtEventRain` — cross-ownership (meteodt writes, TC reads); leave as legacy.

**Total real co-writers: 5 files** (initialize, swap, config_to_variables, meteodt,
swapoutput). All except `initialize.f90` are already state-plumbed.

---

## 5. Init-order analysis

```
swap.f90:184  call config_to_variables(config)   ! seeds tstart/tend/iyear/imonth/dt + 18 config-constants
swap.f90:188  if (flSwapShared) call SharedSimulation(1)
swap.f90:191  call TimeControl(1, state)         ! ← canonical TC init; sets all Group C+D+E
swap.f90:194  call CalcGrid()
swap.f90:195  call soilwater_init(state%soilwater, numnod, numlay)
swap.f90:215  call atmosphere_init(state%atmosphere)
swap.f90:228  call tillage_init(state%tillage, numlay)
swap.f90:235  call heat_init(state)
swap.f90:238  call SoilWater(1, state)
swap.f90:246  call drainage_init(state, config)
...
swap.f90:295   per-step loop
swap.f90:395    call TimeControl(2, state)
```

**Init order observations:**

1. `TimeControl(1, state)` runs at line 191 — BEFORE every other subsystem's init.
   This is where the new `state%timecontrol` is initialized. **The natural place
   for `timecontrol_init`** is to either:
   - (a) add a new `call timecontrol_init(state%timecontrol)` call BEFORE
     `TimeControl(1, state)` — pure zero-init similar to `atmosphere_init`; OR
   - (b) inline the zero-init inside `TimeControl(case=1)` itself — since case 1
     IS the init, adding a `state%timecontrol = timecontrol_state_t()` line at
     the top of case 1 is the cleanest pattern.

   **Recommended (b)** because case 1 already does the canonical init work; adding
   a separate `timecontrol_init` would split init logic across two routines.
   Subordinate atmosphere/snow precedent: `atmosphere_init` is a SEPARATE routine
   only because atmosphere has no init-task dispatcher. TimeControl HAS one
   (case=1), which already plays the init role.

2. `config_to_variables.f90:97, 98, 121` writes `iyear, imonth, dt` BEFORE
   `TimeControl(1)` runs. The transient buffer pattern (`iyear_init_buf`,
   `imonth_init_buf`, `dt_init_buf`) lets the adapter write into module-level
   buffers; `TimeControl(case=1)` then drains them into `state%timecontrol` at
   line 142, 156–157, 234/237. Identical to soil-water-core's `pondini_init_buf`
   / `h_init_buf` pattern.

3. `swap.f90:191` passes `state` to `TimeControl(1)`. Today the arg is `intent(in)`
   — promote to `intent(inout)` (mechanical signature change).

4. The DLL re-init path at `swap.f90:598` calls `TimeControl(1, state)` again with
   `iyear` re-derived. After migration, the re-init writes `state%timecontrol%iyear`.

### Pre-write reader gap

None observed. Every reader of `state%timecontrol%*` consumes values produced by
either `TimeControl(case=1)` (init), `config_to_variables.f90` (pre-init seed
that is drained at TC init time), or the per-step `TimeControl(case=2)` advance.

### Allocation needs

All Group C+D+E fields are **scalars** (or fixed-size `tc_datea(6)` which is a
6-element integer array — still treatable as inline). **No allocatables.** No
per-node array allocation. The `outdat(maout)` / `outdatint(maout)` arrays are
Group A config-constants (stay legacy).

**Recommendation:** initialize inline inside `TimeControl(case=1)` — no separate
`timecontrol_init` subroutine needed. Optionally introduce a tiny
`timecontrol_init(state%timecontrol)` for symmetry with other subsystems (one
line of zero-assignment + drain buffer). The design phase picks one.

---

## 6. Config / Phase 0 candidates

```bash
grep -nE "tstart|tend|nprint|outper|dt(min|max)?|swres|swodat|swheader|nmetdetail|swetsine|swetr|swrain|swmetdetail|swsnow|swhea|swsolu|swdra|swirfix|swinco|swscre|MaxIterTime|msteps|MaxIt|period\b" src/config/simulation_config.f90 src/config/general_config.f90 src/config/meteorology_config.f90 2>/dev/null | head
```

Confirmed covered by typed config:
- `simulation_config_t`: `tstart`, `tend`, `nprintday`, `swmonth`, `period`, `swres`,
  `swodat`, `swyrvar`, `numerical{dt, dtmin, dtmax, MaxIt, MaxBackTr, taccur,
  gwlconv, swkmean, swkimpl, msteps}`
- `meteorology_config_t`: `swetr`, `swetsine`, `swmetdetail`, `swrain`,
  `nmetdetail`, snow{`swsnow`}
- `general_config_t`: `swscre`
- per-subsystem configs: `swdra`, `swhea`, `swsolu`, `swirfix`

Phase 0 candidates:
- `MaxIterTime` — verify covered. (Search needed.)
- `swheader` — verify covered. (Initialized to 0 in adapter:1247.)

```bash
grep -nE "MaxIterTime|swheader" src/io/toml/config_to_variables.f90
```

**Phase 0 candidates: ~0–2 fields.** Likely zero — the static-config side of TC is
the best-covered subsystem (since it was the first migrated to TOML).

---

## 7. Coupling hazards

### H-1: `flzerointr` / `flzerocumu` are cross-subsystem reset gates (DEFER)

`flzerointr` and `flzerocumu` are read by 11–12 distinct files. They are the
canonical gates that drive every `*_state_t%intr%reset()` / `cumu%reset()`
invocation. Migrating them into `state%timecontrol` would force every subsystem
to read `state%timecontrol%flZeroIntr` instead of a bare global.

**Resolution:** **defer both flags to a follow-on arc** ("subsystem-reset
orchestration"). Keep them as legacy globals owned by TimeControl-side writes but
not migrated to state%timecontrol in THIS arc. Justification:
- They form a 14+14 read-site fanout across 25 files (HEAVY).
- The natural follow-on arc would convert them into a `state%timecontrol%advance()`
  return value (`subsystem_reset_gates_t`) rather than retargeting bare reads —
  that is a structurally different refactor.
- Same precedent as tillage's `Bdens`/`ParamVG`: cross-subsystem side-effects stay
  legacy until a dedicated ownership arc (playbook lesson 2026-05-12 #2).

### H-2: `dt` is the most-read time-control field (264 sites)

`dt` is read from every solver/compute path. Today TimeControl writes it 16 times,
soilhydraulics writes it (via `fldecdt → dt = dt/3` in TC case=3, not directly),
and `config_to_variables.f90:121, 576` seeds it pre-init.

**Resolution:** standard reader cutover. Promote `dt → state%timecontrol%dt`.
Strategy B (compile-driven) will surface every hidden import. Estimated 8–12
compile passes.

**Sub-hazard H-2a (`dt` mutated inside Newton-Raphson loop):** `headcalc.f90` and
the macropore iteration cause TimeControl(case=3) to be called WITH `fldecdt =
.true.`, which then writes `dt`. The mutation is single-owner (TC writes) but the
SIGNAL (`fldecdt`) is written by 4 different sites. Already handled by
`timestep_control_mod` since SS-SWST. No extra hazard for THIS arc.

### H-3: `dtEventRain` is meteodt-owned but TimeControl reads it

`meteodt.f90:345` writes `dtEventRain`; `timecontrol.f90:392` reads it. This is a
cross-ownership field where meteodt is the conceptual owner (it's a meteo-derived
quantity). Since the atmosphere arc deferred meteodt (excluded), `dtEventRain`
stays as a legacy global. After meteodt is refactored under a future "meteo arc",
it can move into `state%atmosphere%dtEventRain` or `state%meteodt%dtEventRain`.

**Resolution:** **leave `dtEventRain` as legacy global.** TimeControl's read at
line 392 continues to use the bare global. No design impact.

### H-4: Pre-init writes from `config_to_variables.f90`

`iyear`, `imonth`, and `dt` are written by `config_to_variables.f90:97, 98, 121,
576` BEFORE `TimeControl(1, state)` runs. The state is allocated but
`state%timecontrol` defaults to zero (Fortran default-init of derived type
components). If `config_to_variables` wrote to `state%timecontrol%iyear` directly,
the subsequent `TimeControl(1)` would either overwrite or seed from it.

**Resolution:** **transient-buffer pattern** (soil-water-core lesson #4). Define
module-level buffers `iyear_init_buf`, `imonth_init_buf`, `dt_init_buf` in the
adapter; populate them in `config_to_variables`; let `TimeControl(case=1)` read
them at lines 142, 156–157, 234/237 and write into `state%timecontrol`. Same
pattern as `pondini_init_buf`, `pond_init_buf`, `h_init_buf` (ADR 0038).

Alternative: do the inverse — let `TimeControl(case=1)` initialize state from
`tstart`/config (which it ALREADY does — line 138 `t1900 = tstart` etc.), and
delete the redundant `config_to_variables` writes. This is cleaner but touches
the adapter. Recommendation in design phase.

### H-5: `IterTime(task)` needs state plumbing

`IterTime(task)` at line 639 takes only `task`. It reads/writes `tc_tmptimestart`
and `tc_tmptimeend` (Group D-adjacent globals; CPU-watchdog only). Three call
sites in swap.f90 (`IterTime(1)`, `IterTime(2)`, `IterTime(3)`).

**Resolution:** add `state` arg to `IterTime` (mechanical signature change). Move
`tc_tmptimestart` / `tc_tmptimeend` into `state%timecontrol`. Minor.

### H-6: `flheader`/`flheadirg`/`flIrg1Start` co-written by swapoutput.f90

`swapoutput.f90:2299, 2310` flips these flags as part of output header logic.
After migration these are `state%timecontrol%*` writes. SwapOutput already takes
`state`. **Mechanical dual-write.**

### H-7: `swap.f90:191` `TimeControl(1, state)` signature promotion

Currently `intent(in)`. Promote to `intent(inout)` so TC can write into state.
Single call-site signature change at line 191 + 376 + 395 + 430 + 598 +
swapoutput:4018. All callers pass full state, so the change is mechanical.

### H-8: DLL re-init path at swap.f90:592–598

The DLL exchange path re-derives `iyear` from `Tstart` (line 596–597) and then
calls `TimeControl(1, state)` to re-init. After migration: `state%timecontrol%iyear`
must be set to the re-derived value. The adapter write goes via the transient
buffer or directly into state (since state is available at this call site).

**Resolution:** in the DLL re-init block, write `state%timecontrol%iyear = datea(1)`
directly (state is available). Same for `state%timecontrol%flrunend = .false.` and
`state%timecontrol%fldaystart = .true.` (lines 592–593).

### H-9: Output banner reads dozens of TC fields

`swapoutput.f90` reads `daynr`, `daycum`, `iyear`, `imonth`, `date`, `t1900`,
`period`, `outper`, `nprintcount`, `cntper`, `flprintshort`, `floutput`,
`flbaloutput`, plus many more — for header lines, balance summaries, banner
display, CSV columns. Estimated 60+ read sites in `swapoutput.f90` alone.

**Resolution:** the heaviest reader cutover task in the arc (likely a dedicated
subtask). SwapOutput already takes state.

### H-10: Compile-driven Phase 2 fanout

Expected: **15–25 compile passes** to fully retire the 60+ runtime-state globals,
because every subsystem subprogram reads at least one of `dt`, `daynr`, `t1900`,
`fldaystart`, `period`. This is BY FAR the largest Strategy-B fanout — comparable
to or larger than soil-water-core's 81-call-site `cofgen` pointer-pattern surface,
but here distributed across many small reader updates rather than one binding.

**Resolution:** Strategy B (commenting out globals, iterative compile fixes) is
appropriate. Budget 3–5 compile-driven commits (each one batch of ~40 fixes).

### H-11: `daycrop` is crop-owned (already migrated)

`daycrop` is now owned by the crop subsystem (since SS-CRP), and the commented-out
write at TC lines 296–297 / 495–497 reflects this. Not in scope for THIS arc.

### H-12: `swmeteo` is set by TimeControl but read by ReadMeteoDay/ProcessMeteoDay

TC writes `swmeteo = 1` or `2` depending on crop type (lines 196, 200, 567, 571).
`readmeteo.f90` reads it. **Resolution:** included in Group C migration —
`state%timecontrol%swmeteo`.

---

## 8. Cohort vs flat decision

**Recommendation: FLAT layout.** No `flzerointr` or `flzerocumu` reset clusters
inside TC (TC OWNS those gates; it doesn't consume them). The 61 fields are:

- Per-step / per-day mutated (clock, dt, counters) — no reset gate.
- Monotonic counters (`daycum`, `t1900`, `tcum`) — never reset (except DLL re-init).
- Init-once flags (`flSnow`, `flDrain`, etc.) — set in case 1, then constant.

**Per-day reset cluster?** Examined: `isteps` is reset to 0 at line 552 inside
`flDayEnd` branch. `cntper` increments per day but resets at line 520 inside
`cntper == period` branch. `nprintcount` increments inside per-step branch and
indirectly resets per day. These are scattered per-day "rollovers", not a unified
reset gate. They do NOT justify a `reset_per_day()` cohort method.

**Per-event clusters?** `outper`/`tcumold` reset together inside the print-event
branch. `dtEvent`/`tEvent`/`flTnext` reset together inside the event-completion
branch. Again, scattered single-purpose, not a generalisable cohort.

Boundary precedent (ADR 0035) is flat. Tillage (ADR 0039) is flat. TimeControl
also fits flat.

**Recommended type shape:**

```fortran
type :: timecontrol_state_t
   ! Clock / calendar (Group C — 16 fields)
   real(real64) :: t1900   = 0.0_real64
   real(real64) :: t       = 0.0_real64
   real(real64) :: tcum    = 0.0_real64
   integer      :: daynr   = 0
   integer      :: daycum  = 0
   integer      :: iyear   = 0
   integer      :: iyearm1 = 0
   integer      :: imonth  = 0
   integer      :: daymeteo  = 0
   integer      :: yearmeteo = 0
   character(len=11) :: date = ''
   real(real64) :: timjan1 = 0.0_real64
   integer      :: datea(6) = 0
   real(real32) :: fsec    = 0.0_real32
   integer      :: nextyear = 0
   integer      :: swmeteo = 1

   ! Timestep + schedule (Group D — 19 fields)
   real(real64) :: dt          = 0.0_real64
   real(real64) :: dtold       = 0.0_real64
   real(real64) :: dtEvent     = 0.0_real64
   real(real64) :: tEvent      = 0.0_real64
   real(real64) :: dtprevious  = 0.0_real64
   real(real64) :: tchange     = 0.0_real64
   real(real64) :: tcumold     = 0.0_real64
   integer      :: flprevious  = 1
   logical      :: flTnext     = .false.
   integer      :: isteps      = 0
   integer      :: nprintcount = 1
   integer      :: cntper      = 0
   real(real64) :: outper      = 0.0_real64
   integer      :: ioutdat     = 1
   integer      :: ioutdatint  = 1
   integer      :: rainrec     = 0
   integer      :: wrecord     = 0
   real(real64) :: metperiod   = 0.0_real64

   ! Runtime evaluated booleans (Group E — 26 fields)
   logical      :: flDayStart      = .true.
   logical      :: flDayEnd        = .false.
   logical      :: flRunEnd        = .false.
   logical      :: flYearStart     = .true.
   logical      :: flprintshort    = .false.
   logical      :: floutputshort   = .false.
   logical      :: floutput        = .false.
   logical      :: flbaloutput     = .false.
   logical      :: flheader        = .false.
   logical      :: flheadirg       = .false.
   logical      :: flIrg1Start     = .true.
   logical      :: flUpdMetDet     = .true.
   logical      :: fldecdtmin      = .false.
   logical      :: fldtmin         = .false.
   logical      :: fldtreduce      = .false.
   logical      :: flmetdetail     = .false.
   logical      :: flmeteodt       = .false.
   logical      :: flrainintens    = .false.
   logical      :: fletsine        = .false.
   logical      :: flSnow          = .false.
   logical      :: flDrain         = .false.
   logical      :: flSurfaceWater  = .false.
   logical      :: flTemperature   = .false.
   logical      :: flSolute        = .false.
   logical      :: flIrrigate      = .false.
end type timecontrol_state_t
```

`soilwater_init` signature does NOT change — TimeControl state is independent.

---

## 9. Scope estimate

| Metric | TimeControl (THIS) | Atmosphere (ADR 0037) | Soil-water-core (ADR 0038) | Tillage (ADR 0039) |
|---|---|---|---|---|
| Owned runtime-state globals | ~61 (16 C + 19 D + 26 E) | 40 | 74 | 13 |
| Config-constants (deferred) | ~18 | — | — | 18 |
| External reader files (in-scope) | ~27 | 9 | 22 | 3 |
| External read sites | **~500** (dt heavy) | ~155 | ~370 | ~5 |
| Co-writers (real) | 5 (init, swap, config_to_var, meteodt, swapoutput) | 5 | 8 | 1 |
| Phase 0 candidates | 0–2 | 1–2 | 0 | 0 |
| Cohorts | 0 (flat) | 2 | 2 (twin) | 0 (flat) |
| New init routine | inline in `TimeControl(case=1)` | atmosphere_init | (in soilwater_init) | tillage_init |
| Strategy B compile passes | 15–25 | 8 | 30 | 2 |
| Suggested task count | **~14–18** | 17 | 35+ subtasks | 6 |

### Suggested task decomposition

**Phase 0 — config gaps (optional)**
- TC-0.1: Audit `MaxIterTime` config home; add to `simulation_config_t` if missing.
- TC-0.2: Audit `swheader` config home; verify covered.

**Phase 1 — state-type + dual-write**
- TC-1.1: Create `src/state/timecontrol_state.f90` with `timecontrol_state_t`
  (flat record; 61 fields with explicit defaults). pFUnit init/zero tests.
- TC-1.2: Add `type(timecontrol_state_t) :: timecontrol` to `swap_state_t` in
  `swap_state.f90`. Promote `TimeControl` arg `intent(in)` → `intent(inout)`.
- TC-1.3: Inside `TimeControl(case=1)`: dual-write — every legacy assignment
  (`daynr = nint(t)`, `daycum = 0`, `flZeroIntr = .true.`, etc.) ALSO writes
  `state%timecontrol%daynr`, etc. Drain transient buffers `iyear_init_buf`,
  `imonth_init_buf`, `dt_init_buf` into state.
- TC-1.4: Inside `TimeControl(case=2)`: dual-write every clock/counter/flag
  mutation. Heaviest body.
- TC-1.5: Inside `TimeControl(case=3)` and `case=9`: dual-write `dt`, `dtprevious`,
  `fldecdt`, `fldtmin`, `flTnext`.
- TC-1.6: Co-writer dual-writes — `swap.f90:597` (iyear), `meteodt.f90:361`
  (flUpdMetDet), `swapoutput.f90:2299, 2310` (flheader, flheadirg, flIrg1Start),
  `config_to_variables.f90` adapter writes (transient-buffer pattern).
- TC-1.7: `IterTime(task, state)` signature change + dual-write
  `tc_tmptimestart`/`tc_tmptimeend`.
- TC-1.8: pre-init buffer plumbing (`iyear_init_buf`, `imonth_init_buf`,
  `dt_init_buf` in `config_to_variables.f90`; drain in TC case 1).

**Phase 2 — reader cutover** (heaviest fanout — recommend grouping by subsystem)
- TC-2.1: `swapoutput.f90` + `swap_csv_output.f90` + `macroporeoutput.f90` cutover
  (60+ output-side reads).
- TC-2.2: `soilhydraulics.f90` + `waterbalance.f90` cutover (40+ reads).
- TC-2.3: `meteoday.f90` + `meteodt.f90` + `et.f90` + `snow.f90` cutover (40+).
- TC-2.4: `cropgrowth.f90` + `management_soil.f90` + `oxygenstress.f90` +
  `cropgrass_init.f90` + `tillage.f90` cutover (35+).
- TC-2.5: `drainage.f90` + `divdra.f90` + `surfacewater.f90` + `irrigation.f90`
  cutover (35+).
- TC-2.6: `temperature.f90` + `frozencond.f90` cutover (15+).
- TC-2.7: `boundtop.f90` + `boundbottom.f90` cutover (10+).
- TC-2.8: `macropore.f90` + `macrorate.f90` cutover (10+).
- TC-2.9: `solute.f90` + `agetracer.f90` cutover (15+).
- TC-2.10: `readmeteo.f90` + `swap.f90` cutover (state-machine guards). swap.f90
  main loop reads `flDayStart`, `flDayEnd`, `flRunEnd`, `flYearStart`,
  `flOutput`, `fldecdt`, `fldtreduce` etc. (15+).
- TC-2.11: Strategy B — comment out the ~61 runtime-state globals in
  `variables.f90`; iteratively fix compile errors until clean. 3–5 commits.
- TC-2.12: Drop redundant zero-inits from `initialize.f90` (Group C+D+E lines).

**Phase 3 — ADR + merge**
- TC-3.1: ADR 0040 + playbook lessons + cross-references + merge.

**Total task count estimate: ~17 tasks.** Similar shape to atmosphere (17–19) but
with heavier Phase 2 (10 reader-cutover subtasks vs atmosphere's 10) and slightly
lighter Phase 1 (no cohort design work).

---

## 10. Open questions

- **OQ-1 (Group B `flZeroIntr` / `flZeroCumu` — defer or migrate?):** **Recommendation:
  DEFER.** They have a 14+14 read-site fanout across 25+ files (state modules,
  drainage, surfacewater, meteoday, snow, soilhydraulics, waterbalance, solute,
  agetracer, initialize). They are the canonical "reset gate" signals consumed by
  every `*_state_t%intr%reset()` / `%cumu%reset()` call. The natural future arc is
  a "subsystem-reset orchestration" arc that converts them into a return value
  from `state%timecontrol%advance()`. Migrating them in THIS arc as bare reads
  (without changing semantics) would be mechanical but would add ~28 sites of work
  with no architectural payoff.

- **OQ-2 (Inline init in `TimeControl(case=1)` vs separate `timecontrol_init`):**
  **Recommendation: inline in case=1.** TimeControl ALREADY has a case-1 init
  body (lines 73–258). Adding a separate `timecontrol_init` would split init
  logic across two routines. The atmosphere precedent of a separate `_init`
  routine exists only because atmosphere has no case-1 dispatcher.

- **OQ-3 (Transient-buffer pattern for `iyear`/`imonth`/`dt` pre-init writes):**
  **Recommendation: use the buffer pattern** (soil-water-core lesson #4). Three
  buffers in `config_to_variables.f90` module: `iyear_init_buf`, `imonth_init_buf`,
  `dt_init_buf`. Drained at the top of `TimeControl(case=1)` via `state%timecontrol%* = *_init_buf`.
  Alternative: delete the adapter writes entirely and rely on TimeControl's own
  derivation from `tstart`/config. Cleaner but riskier (changes DLL re-init path).

- **OQ-4 (`IterTime` state plumbing):** **Recommendation: add `state` arg
  unconditionally.** Three callers (swap.f90 lines for `IterTime(1)`,
  `IterTime(2)`, `IterTime(3)`). Move `tc_tmptimestart`/`tc_tmptimeend` into
  state%timecontrol. Trivial cost.

- **OQ-5 (DLL re-init path at `swap.f90:592–598` — direct state write vs second
  TimeControl(1) call):** **Recommendation: keep `TimeControl(1, state)` call;
  write the re-derived `iyear` to `state%timecontrol%iyear = datea(1)` BEFORE
  the call (line 597), and let TC's case-1 init seed everything else from
  `tstart`.** This preserves the existing behavior exactly.

- **OQ-6 (Init-once boolean derived flags — Group E "config-derived" subset):**
  Eleven booleans (`flmetdetail`, `flmeteodt`, `flrainintens`, `fletsine`,
  `flSnow`, `flDrain`, `flSurfaceWater`, `flTemperature`, `flSolute`,
  `flIrrigate`, `flprintshort`) are derived once from config switches at TC
  case=1 and never change. They are on the boundary between config-constant and
  runtime-state. **Recommendation: INCLUDE in `state%timecontrol`** — they are
  written by TC, fit the same record, and avoid splitting the flag set across
  legacy + state.

- **OQ-7 (`dtEventRain` cross-ownership):** **Recommendation: leave as legacy
  global.** meteodt writes it, TC reads it. Will be migrated when the meteo
  refactor delivers `state%atmosphere%dtEventRain` or equivalent. Same pattern
  as atmosphere's deferred `pond` read in `reduceva`.

- **OQ-8 (Group A static config — true config-constants stay legacy):** **Recommendation:
  DEFER all 18 config-constant fields** (tstart, tend, dtmin, dtmax, nprintday,
  period, swres, swodat, swheader, swscre, outdat, outdatint, msteps, MaxIt,
  MaxIterTime, nmetdetail, swmetdetail, swrain, swetsine, swdra, swhea, swsnow,
  swsolu, swirfix). They have typed config homes already; migrating them to
  `state%timecontrol` would duplicate config. Future "config consolidation" arc
  can retire them by routing reads through `config%simulation%*`.

- **OQ-9 (Strategy B fanout estimate):** With ~500 read sites across ~27 files,
  Strategy B compile passes are estimated at **15–25 compile cycles** (each
  cycle batches ~30–50 fixes). Plan for 3–5 separate compile-driven commits to
  keep PR review manageable, similar to soil-water-core's 30-pass profile.

- **OQ-10 (Pre-flight verification — does the arc justify decomposition?):**
  Possible decomposition cuts:
  - **By task type:** clock (Group C) + schedule (Group D) + flags (Group E)
    as three separate sub-arcs.
  - **By reader subsystem:** output side first (TC-2.1) since it's the heaviest;
    physics side second; meteo side third.

  **Recommendation: monolithic arc** — TC's runtime state is a tightly-coupled
  unit. Splitting Groups C/D/E would require maintaining intermediate dual-write
  states across multiple PRs and would amplify Strategy B work. The full
  ~17-task plan should ship as ONE arc.

- **OQ-11 (Daycrop migration status):** Verify that `daycrop` is already owned by
  `state%crop` (or similar) and that the commented-out writes at TC lines 296–297
  / 495–497 are safe to delete during this arc. Not in scope for THIS arc, but
  worth a one-line confirmation in design.

- **OQ-12 (`fldecdt` already in module — verify it's stable in this arc):**
  `fldecdt` lives in `timestep_control_mod` (since SS-SWST). Verify reader
  cutover for `fldecdt` is NOT in scope for THIS arc (it's already migrated).
  The other timestep-decrement flag `fldecdtmin` is still in `variables.f90` and
  IS in scope (Group E).

---

## Discovery summary

- **Owned runtime state: ~61 fields** across 3 groups (16 clock + 19 schedule + 26
  booleans). New record: `timecontrol_state_t`, flat layout, no cohorts.
- **Config-constants deferred: 18 fields** (`tstart`, `tend`, `dt(min/max)`,
  `nprintday`, `period`, etc. — already in `simulation_config_t`).
- **`flZeroIntr` / `flZeroCumu` deferred: 2 fields** with 28 combined read sites
  — natural target for a follow-on "subsystem-reset orchestration" arc.
- **External reader files: ~27 distinct** (output, every solver subsystem, every
  meteo subsystem).
- **External read sites: ~500** (heaviest: `dt` 264 + `t1900` 143 + `period` 92 +
  `outper` 77 + `daynr` 62 + `daycum` 55). 3× the next-heaviest arc.
- **Co-writers: 5 real** (`initialize.f90`, `swap.f90`, `config_to_variables.f90`,
  `meteodt.f90`, `swapoutput.f90`). All except `initialize.f90` already
  state-plumbed.
- **Signature status: `TimeControl(task, state)` already takes state**
  (intent(in)). Promote to intent(inout). `IterTime(task)` needs `state` arg
  added.
- **Init pattern: inline in `TimeControl(case=1)`** — no separate
  `timecontrol_init` routine recommended (TC has its own dispatcher).
- **Hazards: 12 items.** Critical: H-1 deferred `flzero*` gates; H-2 `dt`-as-most-read;
  H-4 pre-init writes via transient buffers; H-9 swapoutput.f90 60+ reads;
  H-10 Strategy B 15–25 compile passes. Manageable: H-5 IterTime, H-6
  swapoutput dual-writes, H-7/8 signature/DLL paths.
- **Cohort: FLAT** — no shared reset gate clusters within the runtime state. ADR
  0035/0039 precedent.
- **Phase 0: 0–2 fields** (verify `MaxIterTime` and `swheader` coverage).
- **Task count estimate: ~17 tasks** across Phase 0/1/2/3. Comparable to
  atmosphere (17–19); larger Phase 2 fanout but no cohort design overhead.
- **soilwater_init signature change needed: NO.** TimeControl is independent of
  per-node arrays.
- **Strategy B applicability: YES** — recommended as the master method for
  retiring the 61 globals after dual-write coverage is verified.

End of discovery.
