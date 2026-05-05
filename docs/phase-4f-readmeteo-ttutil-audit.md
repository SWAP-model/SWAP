# Phase 4f-extend SS-5 — readmeteo.f90 TTutil audit

**Date:** 2026-05-05
**Sub-spec:** SS-5 (legacy reader retirement umbrella)
**Plan:** `docs/superpowers/plans/2026-05-05-ss5-readmeteo-ttutil-deletion.md`
**ADR:** 0014 (Steps 2 and 3)

## Pre-deletion state (HEAD before Commit 2)

`src/io/readmeteo.f90` (772 LoC total) contains the following TTutil-using
sections, all reachable only when `swMetCSV /= 1` (i.e. legacy `.met` or
per-year `.YYY` mode):

| Section | Lines | Description |
|---|---|---|
| `ReadMeteoYear` daily TTutil branch | 106–128 | `rdinit` / `rdacha` / `rdfinr` / `rdfdor` per-year `.YYY` reader |
| `ReadMeteoYear` `MeteoInOneFile` call | 107–108 | Legacy `.met` all-years branch |
| `ReadMeteoYear` detail TTutil branch | 129–144 | `rdinit` / `rdatim` / `rdfinr` / `rdfdor` per-year `.YYY` sub-daily reader |
| `MeteoInOneFile` subroutine | 442–593 | All-years `.met` reader |
| `ReadRainEvents` TTutil tail | 367–end of subroutine | `rdinit` / `rdainr` / `rdfinr` / `rdfdor` per-year `.YYY` rain events reader |
| Control-flow scaffolding | `goto 100`, `100 continue`, `if (swMetCSV == 1)`, `if (swRainCSV == 1)` guards | Strangler-fig harness routing CSV vs TTutil |

## Reachability evidence

Confirmed before SS-5 Commit 1:

- All 5 regression cases (1.hupselbrook, 2.grassgrowth, 4.oxygenstress,
  5.salinitystress, 6.surfacewater) set `metfile = "<NNN>.csv"` in
  `[meteorology.temporal]`. None use `.met` or per-year `.YYY` extensions.
  Verified with `grep -n 'file = ".*\.csv"' tests/swap-cases/toml/*/swap.toml`.
- Case 3 (3.macroporeflow) is excluded per ADR 0011 and is now also stub-
  errored at the TOML boundary by SS-4.
- `readswap.f90:419` sets `swMetFilAll = 0` in the legacy `.swp` pipeline,
  which is no longer invoked by working source code per the umbrella-spec
  retirement gate.

After SS-5 Commit 1 (the `.met`/non-CSV stub-error in
`meteorology_config_validate`), the TTutil branches above are unreachable
from any TOML configuration. Commit 2 deletes the dead code; Commit 3
sweeps the now-dead variables.

## Deletion order (Commit 2)

1. `src/io/toml/config_to_variables.f90`: delete the `swMetFilAll = 1` block
   (lines 176–200). After the stub-error, `metfile` always ends in `.csv`,
   so `swMetFilAll = 0` is the only outcome — the conditional and the
   `MeteoInOneFile(1, ...)` pre-load call are unreachable.
2. `src/io/readmeteo.f90` `ReadMeteoYear`:
   - Delete the `else` branch inside `if (swmetdetail == 0)` (the per-year
     `.YYY` daily TTutil reader, lines 109–128).
   - Delete the `if (swMetFilAll == 1)` line and the `else` keyword (line
     107–109); the remaining body becomes unreachable too — it is the
     `MeteoInOneFile(2, ifnd)` call.
   - Delete the `elseif (swmetdetail == 1)` branch (lines 129–144) — TTutil
     detail reader.
   - Delete the `if (swMetCSV == 1) ... goto 100 ... end if` guard
     (lines 97–104). The two `MeteoCSV*` calls become unconditional.
   - Delete the `100 continue` label (line 146).
3. `src/io/readmeteo.f90` `ReadRainEvents`:
   - Delete everything after the `if (swRainCSV == 1) ... return; end if`
     block (lines 367–439). The TTutil per-year `.YYY` rain reader is
     unreachable.
   - Delete the `if (swRainCSV == 1)` guard (line 312); the body becomes
     the unconditional CSV path.
4. `src/io/readmeteo.f90` end-of-file: delete the entire `MeteoInOneFile`
   subroutine (lines 442–593).

## Variable sweep (Commit 3)

Per ADR 0014 Step 3, the following are deleted from `src/core/variables.f90`:

| Variable | Why dead after Commit 2 |
|---|---|
| `swMetFilAll` | Adapter no longer sets it; reader no longer reads it |
| `swMetCSV` | Only one path remains (CSV) — guard is gone |
| `swRainCSV` | Same as `swMetCSV` |
| `swMetDetCSV` | Same |
| `rainfil` | Used only for `.YYY` rain filename construction |
| `station(366)` | Read by `rdacha` in the deleted per-year daily reader |
| `ad(366)`, `am(366)` integer arrays | Read by TTutil per-year daily reader; CSV path uses `days1900_to_md` |

Survivors (still consumed by `meteoday.f90` / `meteodt.f90`):
`metcsv_dat`, `nmetcsv`, `metcsv_det`, `nmetcsv_det`, `raincsv_dat`,
`nraincsv`, all `det*` simulation arrays, `raintimearray`, `rainamount`,
`nmrain`.

Also delete the `swMetFilAll = 0` assignment at `src/legacy/readswap.f90:419`
(variable no longer exists).
