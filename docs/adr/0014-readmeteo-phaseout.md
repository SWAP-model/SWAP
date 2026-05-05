---
title: "ADR 0014 — readmeteo.f90 TTutil phase-out"
date: 2026-05-01
status: accepted
---

# ADR 0014: readmeteo.f90 TTutil phase-out

## Context

After ADRs 0012 and 0013, the TOML meteo pathway reads daily meteorology and
rain events from CSV companion files.  The TTutil code that reads per-year
`.YYY` files and the all-years `.met` file still exists in `readmeteo.f90`,
kept as the fallback for the legacy `.swp` pipeline.

The result is a strangler-fig: each reader is doubled — a CSV branch at the
top and the old TTutil branch below it, guarded by `swMetCSV` / `swRainCSV`
flags.  Once the legacy `.swp` pipeline is discontinued these branches become
permanently dead.

Two gaps remain before the dead code can actually be deleted:

1. **Sub-daily meteorology (`swmetdetail = 1`)** — the `goto 100` bypass in
   `ReadMeteoYear` already skips the TTutil detail reader, but no CSV
   replacement populates `dettime`/`detrad`/`detrecord`/… .  The simulation
   would run with uninitialised arrays.  The stub allocatables `metcsv_det`
   and `nmetcsv_det` were added in variables.f90 but never wired.

2. **Flow cleanup** — after all TTutil branches are gone the `goto 100` /
   label dance and the CSV flags (`swMetCSV`, `swRainCSV`, `swMetFilAll`) can
   be removed.

## Decision

Phase out the TTutil reader from `readmeteo.f90` in three sequential steps.
Each step is a shippable unit; the regression suite must pass after every step.

---

### Step 1 — Implement CSV sub-daily meteorology

**Prerequisite for deleting the TTutil detail reader.**

**Schema** (`[meteorology.temporal].detail_file`):

```
datetime,record,rad,temp,hum,wind,rain
```

| Column     | Unit         | Description                                              |
|------------|--------------|----------------------------------------------------------|
| `datetime` | ISO datetime | `YYYY-MM-DD HH:MM:SS`; parsed to fractional days since JD 2415020 |
| `record`   | –            | Intra-day slot index 1 … `nmetdetail`                    |
| `rad`      | kJ m⁻² d⁻¹  | Global radiation for the slot                            |
| `temp`     | °C           | Air temperature (single value, not min/max)              |
| `hum`      | kPa          | Actual vapour pressure                                   |
| `wind`     | m s⁻¹        | Wind speed                                               |
| `rain`     | mm           | Rainfall for the slot                                    |

`nmetfile = 17568` (48 slots × 366 days) is the existing fixed-array ceiling
for one year; the same limit applies to the CSV path.

**Changes required:**

| File | Change |
|---|---|
| `meteorology_config.f90` | Add `detail_file` to `meteorology_config_t` |
| `read_meteorology_toml.f90` | Read `detail_file` from `[meteorology.temporal]` |
| `config_to_variables.f90` | Pre-load block: `read_csv_table(['datetime','record','rad','temp','hum','wind','rain'], …)` → `metcsv_det`, `nmetcsv_det`; set `swMetDetCSV = 1` |
| `variables.f90` | Promote `nmetcsv_det`/`metcsv_det` from stubs to live; add `swMetDetCSV` flag |
| `readmeteo.f90` | Add `MeteoCSVDetYear` subroutine: scan `metcsv_det` for current year, populate `dettime`, `detrecord`, `detrad`, `dettav`, `dethum`, `detwind`, `detrain`; initialise `irectotal` and `nofd` |
| `readmeteo.f90` | In `ReadMeteoYear`, replace the bare `goto 100` with `if (swMetCSV==1 .and. swmetdetail==0) … else if (swMetCSV==1 .and. swmetdetail==1) call MeteoCSVDetYear … end if; goto 100` |
| `docs/meteorology.md` | Document `detail_file` key and CSV schema |

**Notes on `irectotal` initialisation:** `meteoday.f90` uses `irectotal` as a
running index into `dettime(:)`.  The existing TTutil path sets
`irectotal = int(t1900 - dettime(1) + 0.1) * nmetdetail`, where `t1900` is
the simulation start time.  `MeteoCSVDetYear` must replicate this exactly
using `dettime(1)` populated from the CSV.

---

### Step 2 — Delete all TTutil branches from `readmeteo.f90`

**Prerequisite:** Step 1 complete; legacy `.swp` pipeline discontinued.

Deletions, in order:

1. **Per-year `.YYY` daily reader** (current lines ~104–122):
   the `else` branch inside `if (swmetdetail == 0)` that calls
   `rdinit` / `rdacha` / `rdfinr` / `rdfdor`.

2. **`MeteoInOneFile`** (current lines ~384–535):
   the entire subroutine; also delete its call site in `config_to_variables.f90`
   and the `swMetFilAll == 1` pre-load block there.

3. **TTutil detail reader** (current lines ~123–138):
   the `elseif (swmetdetail == 1)` block inside `ReadMeteoYear`.

4. **TTutil path in `ReadRainEvents`** (second half of the subroutine):
   everything after the `if (swRainCSV == 1)` block.

After these deletions, simplify the flow:

- Remove `goto 100` and `100 continue` label; `MeteoCSVYear` /
  `MeteoCSVDetYear` are now unconditional calls directly followed by the
  validation block.
- Remove the `if (swMetCSV == 1)` guard — there is only one path.
- `ReadRainEvents` becomes the CSV-only logic with no guard flag needed.

---

### Step 3 — Delete dead variables from `variables.f90`

**Prerequisite:** Step 2 complete.

| Variable | Why dead after Step 2 |
|---|---|
| `swMetFilAll` | Gated the `.met` all-years path; that path is gone |
| `swMetCSV` | No branching needed when only one path exists |
| `swRainCSV` | Same |
| `swMetDetCSV` | Same (added in Step 1, deleted here) |
| `rainfil` | Only used to construct the `.YYY` rain filename |
| `station(366)` | Read by `rdacha` in the old per-year daily reader |
| `ad(mrain)`, `am(mrain)` | Integer day/month arrays read by TTutil; `MeteoCSVYear` uses `days1900_to_md` instead |

Variables that survive (remain load-bearing after Step 2):

- `metcsv_dat`, `nmetcsv` — pre-loaded cache for daily CSV
- `metcsv_det`, `nmetcsv_det` — pre-loaded cache for detail CSV (wired in Step 1)
- `raincsv_dat`, `nraincsv` — pre-loaded cache for rain events CSV
- All `det*` simulation arrays (`dettime`, `detrad`, `detrecord`, …) — consumed by `meteoday.f90`
- `raintimearray`, `rainamount`, `nmrain` — consumed by `meteodt.f90`

---

## Sequencing constraint

Steps 2 and 3 are blocked until:

- The legacy `.swp` pipeline is formally discontinued (its test cases removed
  or ported to TOML), **or**
- The legacy code is moved into a separate `readmeteo_legacy.f90` that the
  `.swp` pathway calls directly, leaving `readmeteo.f90` as TOML-only from
  Step 1 onward.

The second option avoids coupling the meteo phase-out to the full `.swp`
retirement and can proceed earlier.

## Consequences

After all three steps:

- `readmeteo.f90` contains only `ReadMeteoYear`, `ReadRainEvents`,
  `MeteoCSVYear`, `MeteoCSVDetYear`, and the small helpers
  (`days1900_to_md`, `jday`).  No TTutil calls remain.
- The TTutil library can be dropped as a dependency of the TOML build
  target (it may still be needed for the legacy `.swp` target, if that is
  retained as a separate artifact).
- ~150 LoC deleted from `readmeteo.f90`, ~8 variables deleted from
  `variables.f90`.

## Revisit trigger

If `swmetdetail = 1` is determined to be unused by any real-world case and
not worth porting, Step 1 can be replaced by dropping `swmetdetail` as a
supported option in the TOML schema (add it to the deprecated-key list per
ADR 0009's pattern) and deleting the detail arrays from `variables.f90`
outright.

## Progress note 2026-05-05 — Sequencing constraint resolved (SS-5 Commit 1)

The sequencing constraint above is now closed by the second option:
`meteorology_config_validate` rejects any `metfile` not ending in `.csv`
with `ERR_VALIDATION_CROSS_FIELD`. The legacy `.met` and per-year `.YYY`
codepaths are therefore unreachable from any TOML configuration. The
legacy `.swp` pipeline is no longer invoked by working source code (umbrella
spec `2026-05-04-legacy-reader-retirement-design.md`), so no separate
`readmeteo_legacy.f90` is needed. Steps 2 and 3 may proceed.

Audit doc: `docs/archive/2026-phase-4/audits/phase-4f-readmeteo-ttutil-audit.md`.
Plan: `docs/archive/2026-phase-4/plans/2026-05-05-ss5-readmeteo-ttutil-deletion.md`.

## Progress note 2026-05-05 — Step 2 complete (SS-5 Commit 2)

`readmeteo.f90` no longer contains TTutil calls. `MeteoInOneFile` deleted.
The `swMetFilAll = 1` adapter block in `config_to_variables.f90` deleted.
`ReadRainEvents` is CSV-only. Deleted ~280 LoC across `readmeteo.f90`
(file shrank from 772 to 491 lines). The dangling `MeteoInOneFile`
production-side caller in `src/io/readswap.f90` (the `if (swMetFilAll == 1)
call MeteoInOneFile (1, idum)` block, ~line 1747) was also removed to keep
the linker happy — this code path is unreachable in working source per
the umbrella spec retirement gate. Variables (`swMetCSV`, `swRainCSV`,
`swMetFilAll`, `swMetDetCSV`, `rainfil`, `station`, `ad`, `am`) become
trivially redundant — Step 3 sweeps them.

## Progress note 2026-05-05 — Step 3 complete (SS-5 Commit 3)

Dead-variable sweep complete. Deleted from `src/core/variables.f90`:
`swMetFilAll`, `swMetCSV`, `swRainCSV`, `swMetDetCSV`, `rainfil`,
`station(366)`. The `ad(mrain)` and `am(mrain)` arrays are *retained* —
contrary to the audit's initial plan they are still read by `MeteoCSVYear`
(it backfills them from the CSV date column) and by `ReadMeteoYear`'s
date-validation and rain-array init code, so they are alive, not dead.
Adapter assignments deleted from `src/io/toml/config_to_variables.f90`
(`swMetCSV`, `swMetDetCSV`, `swRainCSV`, plus the `rainfil = …` copy).
Legacy `.swp` references deleted from `src/io/readswap.f90`: the
`call rdscha ('rainfil',rainfil)` line and the entire 9-line
`swMetFilAll = 0 / if (.met) … end if` block. Stale comments in
`readmeteo.f90` (`MeteoCSVYear` / `MeteoCSVDetYear` docstrings) updated
to reflect the post-deletion reality. Unused locals `wth` and `getun2`
in `ReadMeteoYear` pruned.

ADR 0014 phase-out is now complete. The TTutil library remains a
build-time dependency of the legacy `.swp` target only (if that target is
retained); the TOML-only build target no longer needs it.

Plan: `docs/archive/2026-phase-4/plans/2026-05-05-ss5-readmeteo-ttutil-deletion.md`.
