---
title: "Phase 4f — Legacy reader analysis"
date: 2026-05-01
---

# Phase 4f: Legacy reader analysis

Survey of what code became dead after the meteo CSV migration, and which
legacy file formats are still consumed in the TOML pathway.  Written to
inform the next brainstorm on reader migration scope.

---

## Dead code confirmed after meteo CSV migration

### 1. `src/utils/ioutils.f90` — fully dead

`io_utils_mod` defines one public routine (`parse_output_extensions`) that
maps 3-letter output-format extension strings to the legacy `sw*` output
switches.  It is **never imported** by any other source file and is not
listed in any `meson.build`.  The output switches it maps to were retired
by ADR 0009.  The module also contains a latent bug (`swcsv = 1` instead
of `csv = 1` in the `'csv'` branch).

**Verdict: safe to delete.**

### 2. `MeteoInOneFile` (readmeteo.f90:384–535) — unreachable from active TOML tests

`MeteoInOneFile` is called when `swMetFilAll == 1`, which is set in
`config_to_variables.f90` only when `metfil` ends in `.met`.  All five
active TOML regression cases now use `.csv`.  Case 3 (macroporeflow) is
excluded from regression and still references a legacy stem without
extension.

The function is still reachable from the legacy `.swp` path
(`readswap.f90:1747–1748`).

**Verdict: unreachable from any active test; not safe to delete yet (still
needed by `.swp` path).**

### 3. Per-year TTutil reader block (readmeteo.f90:104–122) — unreachable from TOML

The `else` branch at line 104 (entered when `swMetCSV==0 &&
swMetFilAll==0`) reads per-year `.YYY` files via TTutil `rdinit` /
`rdacha` / `rdfinr` / `rdfdor`.  In the TOML path this is only reachable
if `metfil` has neither `.csv` nor `.met` extension — an untested scenario.

Still used by `.swp` path.

**Verdict: unreachable from active TOML tests; not safe to delete yet.**

---

## Legacy file formats still consumed in the TOML pathway

These formats have *not* been migrated to CSV companions.  They are called
either from `config_to_variables.f90` directly or from sub-readers invoked
during the simulation loop.

### A. Crop files — `.crp` (TTutil keyword format)

**Where:** `src/crop/cropgrowth.f90`, `src/crop/irrigation.f90`,
`src/crop/management_soil.f90` — all called from the simulation loop.

**Mechanism:** `config_to_variables` strips the `.crp.toml` suffix from
`config%crop%rotation_file(i)` and writes the stem into `cropfil(i)`. The
legacy sub-readers reconstruct `<stem>.crp` and open it via TTutil.

**Files in test cases:** `maizes.crp`, `potatod.crp`, `grassd.crp`,
`grass.crp`.

**Note:** The `.crp.toml` files in the TOML case directories are *not* read
at runtime — they are authoring sources that happen to share the same stem.
The actual runtime file is the TTutil-format `.crp` file.

**Migration complexity:** High — `.crp` files contain nested sections
(crop growth tables, root parameters, WOFOST coefficients, etc.).

### B. Initial soil conditions — `.ini` (TTutil keyword format)

**Where:** `config_to_variables.f90:556–584`, `swinco == 3` branch.

**Mechanism:** When `[soil].inifil` is set, the adapter opens the named
file via `rdinit` and reads `h`, `Tsoil`, `Cml` arrays directly into the
legacy globals.

**Files in test cases:** not exercised by any current TOML case (all
use `swinco != 3`).

**Migration complexity:** Low — scalar arrays, straightforward CSV schema.

### C. Detail meteorology files — per-year `.YYY` (TTutil keyword format)

**Where:** `readmeteo.f90:123–138`, `swmetdetail == 1` branch.

**Mechanism:** When `swmetdetail=1`, the per-year TTutil reader is used for
sub-daily (`nmetdetail` time steps) meteorology regardless of `swMetCSV`.
No CSV path exists for detail files yet; the `datetime_keyed` branch in
`read_csv_table` and the `metcsv_det` cache stubs in `variables.f90` are
infrastructure placeholders.

**Files in test cases:** not exercised by any current TOML case.

**Migration complexity:** Low — same column set as daily; only the
`date` column becomes a `datetime` column.

---

## Summary table

| Format | Reader | Still needed by TOML path | Deletable |
|---|---|---|---|
| Meteo `.met` (all-years) | `MeteoInOneFile` | Only if `.met` specified; no active case | No (`.swp` path) |
| Meteo `.YYY` (per-year) | TTutil `rdinit` stack | Only if no ext; no active case | No (`.swp` path) |
| Meteo `.csv` | `read_csv_table` + `MeteoCSVYear` | Yes (all 5 active cases) | — |
| Crop `.crp` | TTutil `rdinit` stack | Yes (all 5 active cases) | No |
| Initial `.ini` | TTutil `rdinit` stack | Conditionally (swinco=3) | No |
| Detail meteo `.YYY` | TTutil `rdinit` stack | Conditionally (swmetdetail=1) | No |
| `io_utils_mod` | (never called) | No | **Yes** |
