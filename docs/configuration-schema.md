---
title: Configuration schema
author: SWAP modernization team
---

# Configuration schema

This document is the reference for the SWAP TOML input format as it exists at
the rescue baseline. It is extracted directly from the two modernised readers,
which are the schema source of truth:

- `src/io/readswaptoml.f90` — main `.swp` reader (`ReadSwapToml_state`).
- `src/io/readdrainagetoml.f90` — drainage `.dra` reader
  (`ReadDrainageToml_state`).

If a key is not listed here, the TOML reader does not currently parse it,
regardless of whether it appears in example files under `tests/swap-cases/`.
Sample files in the tree contain aspirational keys for future phases; only the
subset documented below actually reaches `swap_state_t`.

## Input file layout

A SWAP simulation at the rescue baseline consumes the following input files:

- `*.swp` — main simulation configuration (TOML). Parsed by
  `ReadSwapToml_state`.
- `*.dra` or `*.dra.toml` — drainage configuration (TOML). Loaded by
  `read_drainage_toml_if_present` when `[drainage].drainage_file` resolves to a
  file that exists. Parsed by `ReadDrainageToml_state`.
- `*.crp` — crop configuration (fixed-format at the baseline). A TOML crop
  reader is a Phase 4 deliverable; current crop inputs are still parsed by the
  legacy fixed-format readers in `src/crop/`.
- `*.met` — meteorology data (fixed-format, day or sub-daily records). The
  path is composed from `[general.paths].atmosphere` and
  `[meteorology].file`.
- `*.bbc` (optional) — bottom boundary data (fixed-format), referenced only
  when the corresponding switch selects an external file.

The `.swp` file is the single entry point: every other input path is resolved
relative to the `.swp` file via the `[general.paths]` subsection.

## TOML conventions

SWAP uses [toml-f](https://github.com/toml-f/toml-f) to parse inputs. toml-f
implements the TOML 1.0 specification strictly, which has consequences that
matter in practice:

- **Tables** are written as `[section]` and nested with dotted headers:
  `[general.paths]`, `[meteorology.evapotranspiration]`,
  `[drainage.basic.resistance]`. Order of sections does not matter.
- **Arrays of tables** use `[[double.brackets]]` and are used in the drainage
  reader for per-level configuration (`[[drainage.basic.resistance.levels]]`
  and `[[drainage.extended.systems]]`). Each occurrence of the doubled header
  starts a new element.
- **Scalars** are typed: integers (`42`), floats (`0.25`), strings (`"hupsel"`),
  booleans (`true` / `false`), and TOML datetimes / local dates
  (`2002-01-01` or `2002-01-01T00:00:00`). toml-f will reject an integer literal
  where a float is expected and vice versa.
- **Arrays** must be homogeneous. `[1, 2, 3]` is legal; `[1, 2.0, 3]` is not.
- **Comments** use `#`.
- Unknown keys are **silently ignored** by the current readers. A Phase 4
  strict-validation pass will reject them; new code should not rely on keys
  that the readers silently discard.
- **Missing required keys** (for example `project` or `start_date`) cause
  `toml-f`'s `get_value` to return a non-zero `stat`; the downstream field
  simply keeps its default value. Missing files (the `.swp` path itself, or a
  referenced drainage TOML) trigger `fatalerr` at read time with the file path
  and toml-f's error message.

## `.swp` reference

The `.swp` reader walks a small set of top-level sections and their
subsections. Every key that reaches `swap_state_t` is enumerated below,
grouped by the TOML section the reader descends into. "Required" means the
reader will only populate the target state field if the key is present;
"optional" means the key is read with a default equal to the current state
value, so omitting it leaves the corresponding default untouched.

### `[general]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `project` | string | optional | `state%time%project` | Project name used in output headers. |
| `swscre` | integer | optional | `state%time%swscre` | Screen output mode. |

### `[general.paths]`

Each path is read as a string; the reader applies `trim` but does not append a
trailing slash. If unset, the default empty string is used and paths resolve
relative to the process working directory.

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `work` | string | optional | `state%time%pathwork` | Directory for output files. |
| `atmosphere` | string | optional | `state%atm%pathatm` | Directory containing `.met` and rainfall files. |
| `crop` | string | optional | `state%crop%pathcrop` | Directory containing crop files. |
| `drain` | string | optional | `state%drain%pathdrain` | Directory containing the drainage TOML. |

### `[simulation]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `start_date` | TOML datetime | optional | `state%time%tstart` | Simulation start. Converted to days since 1900-01-01. |
| `end_date` | TOML datetime | optional | `state%time%tend` | Simulation end. Converted to days since 1900-01-01. |

Both are read as full TOML datetimes; plain local dates (`2002-01-01`) are
accepted and the time component defaults to midnight.

### `[output]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `file_prefix` | string | optional | `state%time%outfil` | Prefix for generated output files. |
| `nprintday` | integer | optional | `state%time%nprintday` | Number of output samples per day. |
| `swheader` | integer | optional | `state%time%swheader` | Emit balance-period headers (0/1). |

### `[output.timing]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `swmonth` | integer | optional | `state%time%swmonth` | Monthly output mode (0=use `period`, 1=end of each month). |
| `period`  | integer | optional | `state%time%period`  | Fixed-interval output in days (when `swmonth=0`). |
| `swres`   | integer | optional | `state%time%swres`   | Reset the interval counter each calendar year. |
| `swodat`  | integer | optional | `state%time%swodat`  | Enable the extra-dates output list. |

### `[meteorology]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `file` | string | optional | `state%atm%metfil` | Name of the `.met` file (resolved against `[general.paths].atmosphere`). |
| `lat`  | real   | optional | `state%atm%lat`   | Station latitude in degrees. |

### `[meteorology.evapotranspiration]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `swetr`       | integer | optional | `state%atm%swetr`       | 0 = Penman–Monteith, 1 = reference ET with crop factors. |
| `alt`         | real    | optional | `state%atm%alt`         | Station altitude (m). |
| `altw`        | real    | optional | `state%atm%altw`        | Anemometer height (m). |
| `angstrom_a`  | real    | optional | `state%atm%angstroma`   | Angstrom A (overcast fraction of extraterrestrial radiation). |
| `angstrom_b`  | real    | optional | `state%atm%angstromb`   | Angstrom B (additional clear-sky fraction). |
| `swdivide`    | integer | optional | `state%atm%swdivide`    | ET/E partitioning: 0 = crop/soil factors, 1 = direct Penman–Monteith. |

### `[meteorology.temporal]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `swmetdetail` | integer | optional | `state%atm%swmetdetail` | 0 = daily records, 1 = sub-daily. |
| `nmetdetail`  | integer | optional | `state%atm%nmetdetail`  | Records per day when `swmetdetail=1`. |
| `swetsine`    | integer | optional | `state%atm%swetsine`    | Sine-wave distribution of daily Tp/Ep. |

### `[meteorology.rainfall]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `swrain`        | integer | optional | `state%atm%swrain`  | Rain input mode (0=daily, 1=daily+intensity, 2=daily+duration, 3=detailed file). |
| `rainfall_file` | string  | optional | `state%atm%rainfil` | Detailed rainfall file (used when `swrain=3`). |

### `[crop.rotation]`

Three parallel arrays define the crop rotation. The reader allocates crop
arrays sized to the length of `crop_file`; shorter `start_date` / `end_date`
arrays are honoured (`min(len(arr), size(state%crop%cropstart))`), longer ones
are truncated.

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `crop_file`  | array of string | optional | `state%crop%cropfil`   | Per-crop file names (without `.CRP` extension). |
| `start_date` | array of string | optional | `state%crop%cropstart` | ISO dates `YYYY-MM-DD` for each crop start. |
| `end_date`   | array of string | optional | `state%crop%cropend`   | ISO dates for each crop end. |

Note that these are **string** arrays, not TOML date arrays — the reader calls
`iso_date_to_t1900` on the trimmed string. Bare TOML local-date literals would
also round-trip through `to_string`, but strings are the documented form.

### `[soil.initial]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `gwli` | real | optional | `state%soil%gwli` (and `state%soil%gwl`) | Initial groundwater level (cm). |

### `[soil.surface]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `pondmx` | real | optional | `state%soil%pondmx` | Minimum ponding depth before runoff (cm). |

### `[soil.evaporation]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `cofredbl` | real | optional | `state%soil%cofred` | Black evaporation coefficient. |
| `cofredbo` | real | optional | `state%soil%cofred` | Boesten/Stroosnijder coefficient. |

The reader writes both keys to the same `cofred` field; whichever appears last
wins. This is a known quirk — see the Discoveries note at the end.

### `[soil.hysteresis]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `swhyst` | integer | optional | `state%soil%swhyst` | Hysteresis switch (0=off). |
| `tau`    | real    | optional | `state%soil%tau`    | Minimum pressure-head difference for wetting/drying transition. |

### `[soil.numerical]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `gwlconv`     | real    | optional | `state%soil%gwlconv`     | Max allowed groundwater-level difference (cm). |
| `critdevh1cp` | real    | optional | `state%soil%CritDevh1Cp` | Relative pressure-head convergence tolerance. |
| `critdevh2cp` | real    | optional | `state%soil%CritDevh2Cp` | Absolute pressure-head convergence tolerance (cm). |
| `maxit`       | integer | optional | `state%soil%msteps`      | Maximum Picard iterations per time step. |

### `[drainage]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `drainage_file` | string | optional | `state%drain%drfil` | Name of the drainage TOML (without `.toml` if you want the `.dra.toml` suffix auto-appended). |

If `drainage_file` resolves to a file that exists, the reader then calls
`ReadDrainageToml_state` on it (see next section). The filename is resolved
against `[general.paths].drain`; an absolute path (leading `/`) is used
verbatim. If the file is missing, the reader logs a warning and continues.

### `[boundary.bottom]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `swbotb`  | integer | optional | `state%boundary%swbotb`  | Bottom boundary type (1..8). |
| `swqhbot` | integer | optional | `state%boundary%swqhbot` | q(h) relationship type for `swbotb=4`. |
| `deepgw`  | real    | optional | `state%boundary%deepgw`  | Hydraulic head in deep aquifer (cm). |
| `hbot`    | real    | optional | `state%boundary%hbot`    | Bottom pressure head (cm). |
| `qbot`    | real    | optional | `state%boundary%qbot`    | Bottom flux (cm/d). |

### `[boundary.top]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `swpondmx` | integer | optional | `state%boundary%swpondmx` | Time-varying ponding threshold switch. |
| `pondmx`   | real    | optional | `state%boundary%pondmx`   | Maximum ponding depth (cm). |
| `rsro`     | real    | optional | `state%boundary%rsro`     | Runoff resistance (d). |
| `runon`    | real    | optional | `state%boundary%runon`    | Runon flux (cm/d). |
| `qtop`     | real    | optional | `state%boundary%qtop`     | Surface flux (cm/d). |

## `.dra` reference

The drainage reader descends from the root `[drainage]` table into two
subsections (`basic`, `extended`) and then parses per-level data from arrays
of tables inside `basic.resistance` and `extended`. The layout below is the
one the reader expects; note that the sample `tests/swap-cases/1.1.hupselbrook-toml/dra.toml`
uses a **different** (older, aspirational) layout and is not accepted by the
current reader. See the Discoveries note at the end.

### `[drainage.basic.general]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `dramet`   | integer | optional | `drain%dramet`   | Drainage method (1, 2, or 3). |
| `swdivd`   | integer | optional | `drain%swdivd`   | Distribute drainage over soil profile. |
| `swdislay` | integer | optional | `drain%swdislay` | Discharge-layer option (0/1/2). |

### `[drainage.basic.general.dislay_table]`

Three parallel arrays, read with `min(len(arr), size(target))`. Only populated
when the corresponding state arrays are already allocated.

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `swtopdislay` | array of integer | optional | `drain%swtopdislay` | Per-level top-of-discharge option. |
| `ztopdislay`  | array of real    | optional | `drain%zTopDisLay`  | Per-level top-of-discharge depth (cm). |
| `ftopdislay`  | array of real    | optional | `drain%fTopDisLay`  | Per-level top-of-discharge fraction. |

### `[drainage.basic.hooghoudt_ernst]`

Used by the Hooghoudt / Ernst drainage method (`dramet=2`). The reader only
writes into index 1 of the per-level arrays for this method.

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `shape`  | real | optional | `drain%shape`  | Shape factor. |
| `entres` | real | optional | `drain%entres` | Drain entry resistance (d). |
| `basegw` | real | optional | `drain%basegw` | Depth of impervious layer (cm). |
| `zbotdr` | real | optional | `drain%zbotdr(1)` | Drain bottom depth, level 1 (cm). |
| `wetper` | real | optional | `drain%wetper(1)` | Wet perimeter of drain, level 1 (cm). |

### `[drainage.basic.resistance]`

Scalar parameters governing the multi-level resistance method (`dramet=3`).

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `nrlevs`      | integer | optional | `drain%nrlevs`      | Number of drainage levels. |
| `swintfl`     | integer | optional | `drain%swnrsrf`     | Interflow switch (note: writes to `swnrsrf`). |
| `cofintflb`   | real    | optional | `drain%cofintfl`    | Interflow coefficient. |
| `expintflb`   | real    | optional | `drain%expintfl`    | Interflow exponent. |
| `swtopnrsrf`  | integer | optional | `drain%SwTopnrsrf`  | Top-of-interflow switch. |

### `[[drainage.basic.resistance.levels]]` (array of tables)

Each table entry is read up to `min(len(arr), drain%nrlevs)`. The reader
supports an explicit `level` key to place a table at a specific index;
otherwise the 1-based array position is used. Entries that fall outside the
allocated array range are skipped silently.

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `level`  | integer | optional | loop index    | Destination index (defaults to array position). |
| `drares` | real    | optional | `drain%drares(idx)` | Drainage resistance (d). |
| `infres` | real    | optional | `drain%infres(idx)` | Infiltration resistance (d). |
| `swallo` | integer | optional | `drain%swallo(idx)` | 1=drainage+infiltration, 2=drainage only, 3=infiltration only. |
| `l`      | real    | optional | `drain%L(idx)`      | Drain spacing (m). |
| `zbotdr` | real    | optional | `drain%zbotdr(idx)` | Drain bottom depth (cm). |
| `swdtyp` | integer | optional | `drain%swdtyp(idx)` | Drain type (0=channel, 1=tube, 2=interflow). |

### `[drainage.extended.characteristics]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `nrsrf`        | integer | optional | `drain%nrlevs`       | Number of levels (also written by `resistance.nrlevs`). |
| `swnrsrf`      | integer | optional | `drain%swnrsrf`      | Interflow switch. |
| `rsurfdeep`    | real    | optional | `drain%rsurfdeep`    | Deep interflow resistance. |
| `rsurfshallow` | real    | optional | `drain%rsurfshallow` | Shallow interflow resistance. |

### `[[drainage.extended.systems]]` (array of tables)

Same indexing rule as `resistance.levels`, bounded by `drain%nrlevs` and
`size(drain%rdrain)`.

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `lev`     | integer | optional | loop index        | Destination index. |
| `swdtyp`  | integer | optional | `drain%swdtyp(idx)` | Drain type. |
| `l`       | real    | optional | `drain%L(idx)`      | Drain spacing (m). |
| `zbotdre` | real    | optional | `drain%zbotdr(idx)` | Drain bottom depth (cm). Note: `zbotdre` here, `zbotdr` in basic. |
| `gwlinf`  | real    | optional | `drain%gwlinf(idx)` | GWL below which no infiltration (cm). |
| `rdrain`  | real    | optional | `drain%rdrain(idx)` | Drainage resistance per level (d). |
| `rinfi`   | real    | optional | `drain%rinfi(idx)`  | Infiltration resistance per level (d). |
| `rentry`  | real    | optional | `drain%rentry(idx)` | Entry resistance (d). |
| `rexit`   | real    | optional | `drain%rexit(idx)`  | Exit resistance (d). |
| `widthr`  | real    | optional | `drain%widthr(idx)` | Drain / channel width (cm). |
| `taludr`  | real    | optional | `drain%taludr(idx)` | Talus slope of channel. |

## `.crp` reference (baseline)

At the rescue baseline, crop configuration uses fixed-format `.crp` files, not
TOML. The TOML crop reader is a Phase 4 deliverable. Current fixed-format crop
inputs are read by `readcropfixed`, `readgrass`, and `readwofost` in
`src/crop/`. See the upstream SWAP 4.2.0 manual (preserved at
`legacy/swap-4.2.0/doc/`) for the fixed-format crop-file layout. This section
will be filled in when Phase 4 adds a `readcrop_toml.f90` module with a
defined schema.

## Example

The smallest `.swp` that exercises the supported keys looks like this. Every
line is optional; sections may be omitted entirely and the state field keeps
its default value.

```toml
# Minimal SWAP TOML configuration

[general]
project = "hupsel"
swscre  = 0

[general.paths]
work       = "./"
atmosphere = "./"
crop       = "./"
drain      = "./"

[simulation]
start_date = 2002-01-01
end_date   = 2002-12-31

[output]
file_prefix = "result"
nprintday   = 1
swheader    = 0

[output.timing]
swmonth = 1
period  = 1
swres   = 0
swodat  = 0

[meteorology]
file = "283.met"
lat  = 52.0

[meteorology.evapotranspiration]
swetr      = 0
alt        = 10.0
altw       = 10.0
angstrom_a = 0.25
angstrom_b = 0.5
swdivide   = 1

[meteorology.temporal]
swmetdetail = 0
nmetdetail  = 24
swetsine    = 0

[meteorology.rainfall]
swrain = 0

[crop.rotation]
crop_file  = ["maizes"]
start_date = ["2002-05-01"]
end_date   = ["2002-10-15"]

[soil.initial]
gwli = -75.0

[soil.surface]
pondmx = 0.2

[soil.numerical]
gwlconv     = 100.0
critdevh1cp = 0.01
critdevh2cp = 0.1
maxit       = 30

[drainage]
drainage_file = "swap"

[boundary.bottom]
swbotb = 6

[boundary.top]
rsro = 0.5
```

A fuller example lives at `tests/swap-cases/1.1.hupselbrook-toml/swap.toml`.
That file contains many additional keys not yet read by `ReadSwapToml_state`
(irrigation, heat, solute, soil hydraulic tables, macropores, snow, frost);
those keys are ignored today and are targets for later rescue phases.

## Validation and errors

Validation in the current readers is deliberately thin; the rescue plan
promotes it in Phase 4.

- **Missing `.swp` file**, or any `.swp` that fails to parse as TOML, triggers
  `fatalerr` from `ReadSwapToml_state` with the file path and the toml-f error
  message (typically a line and column).
- **Missing `[drainage].drainage_file`** is tolerated — the reader simply
  skips the detailed drainage step. If `drainage_file` is set but the
  resolved path does not exist, the reader logs an info message and continues
  with defaults instead of aborting.
- **Missing required keys** within a present section return a non-zero `stat`
  from `get_value`; the downstream state field keeps whatever value it held
  before the call (normally a type default from `swap_state_init`). There is
  currently no "this key was required" check — every read uses `default=...`
  or is guarded by the `istat == 0` conditional.
- **Wrong types** (for example a string where an integer is expected) are
  rejected by toml-f itself with a type-mismatch error; because the readers
  pass `stat=istat`, the single mis-typed key is silently skipped rather than
  aborting the whole run.
- **Unknown keys** are ignored. toml-f does not raise on them, and the
  readers do not iterate over the full table to look for extras. Phase 4 will
  add a strict-mode pass that diffs the observed keys against a declared
  schema; new code should not rely on silently-ignored keys as a free-form
  extension point.
- **Duplicate table headers**, **unterminated strings**, and other
  syntactically invalid TOML are rejected by `toml_load` before any reader
  logic runs; the error bubbles up through `err%message`.

## Discoveries

Worth noting for readers who will extend these schemas:

1. `[soil.evaporation].cofredbl` and `cofredbo` are both read into the same
   `state%soil%cofred` field. Whichever key appears later in the file wins.
   The underlying physics uses a single reduction coefficient but the
   intent of the two keys differs (Black vs. Boesten/Stroosnijder); this is
   a latent bug or at minimum an ambiguous schema.
2. `[drainage.basic.resistance].swintfl` writes to `drain%swnrsrf`, the same
   target as `[drainage.extended.characteristics].swnrsrf`. The names suggest
   independent switches but they share state.
3. `[drainage.basic.resistance].nrlevs` and
   `[drainage.extended.characteristics].nrsrf` also share `drain%nrlevs`;
   setting both to different values produces last-write-wins.
4. The sample file `tests/swap-cases/1.1.hupselbrook-toml/dra.toml` predates
   the current reader and uses top-level sections (`[general]`,
   `[drainage_resistance]`, `[drainage_extended]`) that the reader does not
   recognise. A real drainage TOML must live under the `[drainage.basic]` /
   `[drainage.extended]` hierarchy shown above.
5. Crop rotation dates in `[crop.rotation]` are read as **string arrays**
   (parsed by a custom `iso_date_to_t1900`), not as TOML date arrays. The
   simulation start/end in `[simulation]` are TOML datetimes. The two
   conventions exist side-by-side in the same file.

## Typed config hierarchy

The TOML key layout in this document is mirrored by a typed config
hierarchy in `src/config/`:

| TOML section | Fortran type | Module |
|---|---|---|
| `[general]` | `general_config_t` | `general_config_mod` |
| `[simulation]` | `simulation_config_t` | `simulation_config_mod` |
| `[meteorology]` | `meteorology_config_t` | `meteorology_config_mod` |
| `[drainage]` | `drainage_config_t` | `drainage_config_mod` |
| `[soil]` | `soil_config_t` | `soil_config_mod` |
| `[crop]` | `crop_config_t` | `crop_config_mod` |
| (top-level) | `swap_config_t` | `swap_config_mod` |

Every type exposes `validate(errors)` and `finalize(errors)` as
type-bound procedures. See `docs/validation.md` for the rules and
`docs/error-handling.md` for how errors propagate.
