---
title: Configuration schema
author: SWAP modernization team
---

# Configuration schema

This document is the reference for the SWAP TOML input format as it exists at
the Phase 4c-b rescue baseline. It is extracted directly from the modernised
readers, which are the schema source of truth:

- `src/io/toml/load_swap_config.f90` — top-level dispatcher.
- `src/io/toml/read_general_toml.f90` — `[general]` section.
- `src/io/toml/read_simulation_toml.f90` — `[simulation]` section (incl. `[simulation.output]`).
- `src/io/toml/read_meteorology_toml.f90` — `[meteorology]` section.
- `src/io/toml/read_drainage_toml.f90` — `[drainage]` section + external file reference.
- `src/io/toml/read_soil_toml.f90` — `[soil]` section.
- `src/io/toml/read_crop_toml.f90` — `[[crop.rotation]]` array-of-tables + cross-file dispatch.
- `src/io/toml/read_cropfixed_toml.f90` — type-1 `.crp.toml` schema.
- `src/io/toml/read_cropwofost_toml.f90` — type-2 `.crp.toml` schema (Phase 4c-b).
- `src/io/toml/read_cropgrass_toml.f90` — type-3 `.crp.toml` schema.

If a key is not listed here, the TOML reader does not currently parse it,
regardless of whether it appears in example files under `tests/swap-cases/`.
Sample files in the tree contain aspirational keys for future phases; only the
subset documented below actually reaches `swap_config_t`.

**See also:**
- `docs/toml-format-guide.md` — general TOML conventions used in SWAP.
- Each module's source (listed above) for authoritative field and validator definitions.

## Input file layout

A SWAP simulation at the rescue baseline consumes the following input files:

- `*.swp` — main simulation configuration (TOML). Parsed by
  `ReadSwapToml_state`.
- `*.dra.toml` — drainage configuration (TOML). Loaded when `[drainage].file`
  is present and the referenced file exists. Parsed by `read_drainage_toml`.
- `*.crp.toml` — crop configuration (TOML). Loaded for each `[[crop.rotation]]`
  entry that specifies `file = "..."` and the file exists. Parsed by
  `read_cropfixed_toml` (type 1) or `read_cropgrass_toml` (type 3).
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
| `start_date` | TOML date | **required** | `config%simulation%tstart` | Simulation start. Converted to days since 1900-01-01. |
| `end_date`   | TOML date | **required** | `config%simulation%tend`   | Simulation end. Converted to days since 1900-01-01. |
| `nprintday`  | integer   | optional     | `config%simulation%nprintday` | Number of output samples per day (default 1). |

Both dates are read as TOML local dates (`2002-01-01`) or datetimes
(`2002-01-01T00:00:00`). `start_date` and `end_date` are the only
**required** keys in the entire `.swp` file.

### `[simulation.output]`

Output timing is nested under `[simulation.output]`, **not** at top level.

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `swmonth` | integer | optional | `config%simulation%swmonth` | Monthly output mode (0=use `period`, 1=end of each month). |
| `period`  | integer | optional | `config%simulation%period`  | Fixed-interval output in days (when `swmonth=0`; default 1). |
| `swres`   | integer | optional | `config%simulation%swres`   | Reset the interval counter each calendar year (0/1). |
| `swodat`  | integer | optional | `config%simulation%swodat`  | Enable the extra-dates output list (0/1). |
| `swyrvar` | integer | optional | `config%simulation%swyrvar` | Year-varying output switch (0/1). |

> **Phase 4b/4c change:** Keys `file_prefix`, `swheader` from the Phase-2
> aspirational schema were never implemented in the typed TOML pipeline. The
> legacy `readswap.f90` still populates `state%time%outfil` from a separate
> fixed-format section; the TOML path does not yet map these.

### `[meteorology]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `file` | string  | optional | `config%meteo%metfil` | Name of the `.met` file (resolved against `[general.paths].atmosphere`). |
| `lat`  | real    | optional | `config%meteo%lat`    | Station latitude in degrees. |
| `alt`  | real    | optional | `config%meteo%alt`    | Station altitude (m). **At section root, not under evapotranspiration.** |
| `altw` | real    | optional | `config%meteo%altw`   | Anemometer height (m; default 2.0). **At section root, not under evapotranspiration.** |

> **Phase 4b change:** `alt` and `altw` were previously documented under
> `[meteorology.evapotranspiration]`. The reader (`read_meteorology_toml.f90`)
> reads them directly from the `[meteorology]` table root.

### `[meteorology.evapotranspiration]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `swetr`       | integer | optional | `config%meteo%swetr`      | 0 = Penman–Monteith, 1 = reference ET with crop factors. |
| `angstrom_a`  | real    | optional | `config%meteo%angstroma`  | Angstrom A (overcast fraction of extraterrestrial radiation; default 0.25). |
| `angstrom_b`  | real    | optional | `config%meteo%angstromb`  | Angstrom B (additional clear-sky fraction; default 0.50). |
| `swdivide`    | integer | optional | `config%meteo%swdivide`   | ET/E partitioning: 0 = crop/soil factors, 1 = direct Penman–Monteith. |

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

### `[crop]`

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `swcrop` | integer | optional | `config%crop%swcrop` | Crop simulation switch (0 = no crops; default 0). |

### `[[crop.rotation]]` (array of tables)

The crop rotation is an **array of tables** (double-bracket `[[crop.rotation]]`),
not parallel string arrays. Each entry in the array represents one crop period.

> **Phase 4b change:** The old Phase-2 aspirational schema used three parallel
> arrays (`crop_file`, `start_date`, `end_date`) under `[crop.rotation]`.
> The Phase 4c-a reader uses `[[crop.rotation]]` array-of-tables with
> TOML date values.

| Key | Type | Required | State target | Description |
|---|---|---|---|---|
| `start` | TOML date | optional | `config%crop%rotation_start(i)` | Crop period start (days since 1900). |
| `end`   | TOML date | optional | `config%crop%rotation_end(i)`   | Crop period end (days since 1900). |
| `file`  | string    | optional | `config%crop%rotation_file(i)`  | Path to the `.crp.toml` file for this period (relative to `.swp` file). |
| `type`  | integer   | optional | `config%crop%rotation_type(i)`  | Crop type: 1=fixed, 2=WOFOST general (Phase 4c-b), 3=WOFOST grass. Default 0. |

When `file` is present and the referenced file exists, the loader reads the
per-crop TOML and dispatches by `type` (see "Cross-file references" below).
If the file is absent, the entry is silently skipped — this allows test cases
not yet converted to `.crp.toml` to continue loading.

Example:

```toml
[crop]
swcrop = 1

[[crop.rotation]]
start = 2002-05-01
end   = 2002-09-30
file  = "maizes.crp.toml"
type  = 1

[[crop.rotation]]
start = 2003-05-01
end   = 2003-09-30
file  = "maizes.crp.toml"
type  = 1
```

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
| `file` | string | optional | _(resolved by loader)_ | Path to an external `*.dra.toml` drainage file, relative to the `.swp` file. |

When `file` is present, the loader reads the referenced file and parses its
`[drainage]` table (see "Cross-file references" below and the `.dra` reference
section). If `file` is absent, drainage keys may appear inline under
`[drainage]` in the same `.swp` file.

> **Phase 4b→4c change:** The old key `drainage_file` under `[general.paths]`
> is superseded. The Phase 4c-a reader uses `[drainage].file` for the
> external file reference. If `file` is absent and no inline drainage keys are
> present, the drainage config struct retains its defaults.

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

## `.crp.toml` reference (Phase 4c-a)

Phase 4c-a introduces TOML crop files. Each `[[crop.rotation]]` entry with
`file = "..."` points to a `.crp.toml` file. The schema of that file depends
on the `type` key.

### `*.crp.toml` — type 1 (fixed crop)

Parsed by `read_cropfixed_toml.f90`; config type `cropfixed_config_t`
(`src/config/cropfixed_config.f90`).

Sections: `[phenology]`, `[light]`, `[root]`, `[water_stress]`,
`[salinity]`, `[interception]`.

Example:

```toml
[phenology]
idev = 1
lcc  = 168

[light]
kdif = 0.6
kdir = 0.6

[root]
rdi = 5.0
rri = 1.2
rdc = 100.0

[water_stress]
hlim1  = -10.0
hlim2u = -25.0
hlim2l = -200.0
hlim3h = -400.0
hlim3l = -600.0
hlim4  = -8000.0
adcrh  = 0.5
adcrl  = 0.1
rsc    = 70.0

[salinity]
ecmax  = 1.7
ecslop = 12.0

[interception]
cofab = 0.25
```

#### `[phenology]` (type 1)

| Key | Type | Default | Description |
|---|---|---|---|
| `idev` | integer | 1 | Development mode: 1=fixed period, 2=temperature-sum-based. |
| `lcc`  | integer | 0 | Length of crop cycle (days); used when `idev=1`. |

#### `[light]` (type 1)

| Key | Type | Default | Description |
|---|---|---|---|
| `kdif` | real | 0.0 | Diffuse light extinction coefficient. |
| `kdir` | real | 0.0 | Direct light extinction coefficient. |

#### `[root]` (type 1)

| Key | Type | Default | Description |
|---|---|---|---|
| `rdi` | real | 0.0 | Initial rooting depth (cm). |
| `rri` | real | 0.0 | Daily root extension rate (cm/d). |
| `rdc` | real | 0.0 | Maximum rooting depth (cm). |

#### `[water_stress]` (type 1 and 3)

Feddes pressure-head thresholds (all in cm, negative values):

| Key | Type | Default | Description |
|---|---|---|---|
| `hlim1`  | real | 0.0 | Saturation threshold (near 0, least negative). |
| `hlim2u` | real | 0.0 | Upper anaerobiosis threshold. |
| `hlim2l` | real | 0.0 | Lower anaerobiosis threshold. |
| `hlim3h` | real | 0.0 | High-transpiration wilting start. |
| `hlim3l` | real | 0.0 | Low-transpiration wilting start. |
| `hlim4`  | real | 0.0 | Wilting point (most negative). |
| `adcrh`  | real | 0.0 | Critical fraction reduction at hlim3h. |
| `adcrl`  | real | 0.0 | Critical fraction reduction at hlim3l. |
| `rsc`    | real | 0.0 | Crop resistance for Penman-Monteith ET (s/m). |

#### `[salinity]` (type 1 and 3)

| Key | Type | Default | Description |
|---|---|---|---|
| `ecmax`  | real | 0.0 | Threshold EC above which yield reduction starts (dS/m). |
| `ecslop` | real | 0.0 | Slope of yield reduction per unit EC above `ecmax` (%/dS/m). |

#### `[interception]` (type 1 and 3)

| Key | Type | Default | Description |
|---|---|---|---|
| `cofab` | real | 0.0 | Interception coefficient (cm/LAI per event). |

### `*.crp.toml` — type 2 (WOFOST general)

Phase 4c-b adds the WOFOST type-2 schema. Parsed by
`read_cropwofost_toml.f90`; config type `cropwofost_config_t`
(`src/config/cropwofost_config.f90`).

The schema mirrors the legacy `readcropwofost`/`crpgrowth` parameter set,
broken into 21 nested sub-tables. All sections are optional: missing
sections leave their fields at type defaults and the validator only flags
inconsistencies when the relevant switch is enabled.

Sections (read in order):

`[preparation]`, `[sowing]`, `[germination]`, `[harvest]`,
`[crop_factor]`, `[phenology]`, `[initial]`, `[green_area]`,
`[assimilation]`, `[conversion]`, `[respiration]`, `[partitioning]`,
`[death]`, `[root]`, `[oxygen_stress]`, `[drought_stress]`, `[salinity]`,
`[compensate]`, `[interception]`, `[co2]`, `[management]`.

#### Table encoding

WOFOST uses many AFGEN tables (development-stage → coefficient curves).
TOML encodes them as **arrays of arrays**, one inner array per row, with a
fixed column count (2 for AFGEN x/y pairs, 6 for `gashtb` Gash interception
parameters). The reader (`read_table_2d`) leaves the field unallocated when
the key is absent and emits `ERR_PARSE_TYPE_MISMATCH` for ragged or
wrong-width rows. Validators then enforce 2 ≤ nrows ≤ 15 for AFGEN tables.

```toml
[phenology]
dtsmtb = [
  [ 0.0,  0.0],
  [30.0, 30.0],
  [45.0, 30.0],
]
```

#### `[preparation]`

| Key | Type | Default | Description |
|---|---|---|---|
| `swprep`       | integer | 0   | 0=no land preparation, 1=delay sowing until conditions are met. |
| `zprep`        | real    | 0.0 | Soil-depth at which `hprep` is evaluated (cm, ≤0). |
| `hprep`        | real    | 0.0 | Pressure-head threshold for preparation (cm). |
| `maxprepdelay` | integer | 0   | Maximum allowed days to delay preparation (1..366 when active). |

#### `[sowing]`

| Key | Type | Default | Description |
|---|---|---|---|
| `swsow`       | integer | 0   | 0=fixed sowing date, 1=delayed sowing on conditions. |
| `zsow`        | real    | 0.0 | Depth at which `hsow`/`tempsow` are evaluated (cm). |
| `hsow`        | real    | 0.0 | Pressure-head threshold for sowing (cm). |
| `ztempsow`    | real    | 0.0 | Depth at which soil temperature is evaluated (cm). |
| `tempsow`     | real    | 0.0 | Soil-temperature threshold for sowing (°C, 0..30). |
| `maxsowdelay` | integer | 0   | Maximum allowed sowing delay in days (1..366). |

#### `[germination]`

| Key | Type | Default | Description |
|---|---|---|---|
| `swgerm`     | integer | 0   | 0=immediate emergence, 1=tsum-driven, 2=tsum + soil-water. |
| `tsumemeopt` | real    | 0.0 | Optimum temperature sum for emergence (°C·d, 0..1000). |
| `tbasem`     | real    | 0.0 | Base temperature for emergence (°C, 0..40). |
| `teffmx`     | real    | 0.0 | Maximum effective temperature (°C, 0..40). |
| `hdrygerm`   | real    | 0.0 | Pressure head at which germination stalls (dry side, cm). |
| `hwetgerm`   | real    | 0.0 | Pressure head at which germination stalls (wet side, cm). |
| `zgerm`      | real    | 0.0 | Depth at which germination is evaluated (cm). |
| `agerm`      | real    | 0.0 | Empirical germination coefficient. |

#### `[harvest]`

| Key | Type | Default | Description |
|---|---|---|---|
| `dvsend` | real    | 0.0 | Development stage at harvest (typically 2.0 for cereals). |
| `swharv` | integer | 0   | 0=harvest at `dvsend`, 1=harvest at `crop.rotation.end`. |

#### `[crop_factor]`

| Key | Type | Default | Description |
|---|---|---|---|
| `swcf`   | integer            | 0   | 1=crop factor for ETref, 2=Penman–Monteith with `rsc`/`rsw`. |
| `albedo` | real               | 0.0 | Crop albedo (0..1). |
| `rsc`    | real               | 0.0 | Crop resistance (s/m). |
| `rsw`    | real               | 0.0 | Wet canopy resistance (s/m). |
| `cftb`   | array of [dvs, cf] | —   | Crop factor vs. development stage (or vs. day for type=1). |
| `chtb`   | array of [dvs, h]  | —   | Crop height vs. development stage (cm). |

#### `[phenology]`

| Key | Type | Default | Description |
|---|---|---|---|
| `idsl`     | integer                | 0   | 0=temperature-only, 1=temperature+daylength, 2=add vernalisation. |
| `tsumea`   | real                   | 0.0 | Temperature sum from emergence to anthesis (°C·d). |
| `tsumam`   | real                   | 0.0 | Temperature sum from anthesis to maturity (°C·d). |
| `dlo`      | real                   | 0.0 | Optimum daylength (h, used when `idsl≥1`). |
| `dlc`      | real                   | 0.0 | Critical daylength (h). |
| `vernsat`  | real                   | 0.0 | Saturated vernalisation requirement (d, used when `idsl=2`). |
| `vernbase` | real                   | 0.0 | Base vernalisation requirement (d). |
| `verndvs`  | real                   | 0.0 | Development stage at which vernalisation completes. |
| `dtsmtb`   | array of [t, dtsm]     | —   | Daily temperature-sum response curve. |
| `verntb`   | array of [t, vern]    | —   | Vernalisation rate vs. temperature (used when `idsl=2`). |

#### `[initial]`

| Key | Type | Default | Description |
|---|---|---|---|
| `tdwi`   | real | 0.0 | Initial total dry-matter weight (kg/ha). |
| `laiem`  | real | 0.0 | Leaf area index at emergence (m²/m²). |
| `rgrlai` | real | 0.0 | Maximum daily LAI growth rate (m²/m²/d). |

#### `[green_area]`

| Key | Type | Default | Description |
|---|---|---|---|
| `spa`   | real                | 0.0 | Specific pod area (ha/kg). |
| `ssa`   | real                | 0.0 | Specific stem area (ha/kg). |
| `span`  | real                | 0.0 | Maximum life-span of leaves under optimum conditions (d). |
| `tbase` | real                | 0.0 | Lower threshold temperature for ageing of leaves (°C). |
| `slatb` | array of [dvs, sla] | —   | Specific leaf area vs. development stage (ha/kg). |

#### `[assimilation]`

| Key | Type | Default | Description |
|---|---|---|---|
| `kdif`   | real                  | 0.0 | Diffuse-light extinction coefficient. |
| `kdir`   | real                  | 0.0 | Direct-light extinction coefficient. |
| `eff`    | real                  | 0.0 | Light-use efficiency (kg·ha⁻¹·h⁻¹/(J·m⁻²·s⁻¹)). |
| `amaxtb` | array of [dvs, amax]  | —   | Maximum assimilation rate vs. development stage (kg·ha⁻¹·h⁻¹). |
| `tmpftb` | array of [t, factor]  | —   | Day-temperature reduction factor for `amax`. |
| `tmnftb` | array of [t, factor]  | —   | Night-temperature reduction factor for `amax`. |

#### `[conversion]`

| Key | Type | Default | Description |
|---|---|---|---|
| `cvl` | real | 0.0 | Conversion efficiency assimilates → leaves (kg/kg). |
| `cvo` | real | 0.0 | Conversion efficiency assimilates → storage organs. |
| `cvr` | real | 0.0 | Conversion efficiency assimilates → roots. |
| `cvs` | real | 0.0 | Conversion efficiency assimilates → stems. |

#### `[respiration]`

| Key | Type | Default | Description |
|---|---|---|---|
| `q10`    | real                | 0.0 | Q10 of maintenance respiration. |
| `rml`    | real                | 0.0 | Maintenance respiration rate of leaves (kg·CH₂O·kg⁻¹·d⁻¹). |
| `rmo`    | real                | 0.0 | Maintenance respiration of storage organs. |
| `rmr`    | real                | 0.0 | Maintenance respiration of roots. |
| `rms`    | real                | 0.0 | Maintenance respiration of stems. |
| `rfsetb` | array of [dvs, frac]| —   | Fraction of maintenance-respiration reduction with senescence. |

#### `[partitioning]`

Carbohydrate partitioning fractions vs. development stage. All four
tables must sum to 1 at each row.

| Key | Type | Default | Description |
|---|---|---|---|
| `frtb` | array of [dvs, fr] | — | Fraction of total dry matter to roots. |
| `fltb` | array of [dvs, fl] | — | Fraction of above-ground DM to leaves. |
| `fstb` | array of [dvs, fs] | — | Fraction of above-ground DM to stems. |
| `fotb` | array of [dvs, fo] | — | Fraction of above-ground DM to storage organs. |

#### `[death]`

| Key | Type | Default | Description |
|---|---|---|---|
| `perdl`  | real                  | 0.0 | Maximum relative leaf death rate due to water stress (d⁻¹). |
| `rdrrtb` | array of [dvs, rdrr]  | —   | Relative root death rate vs. development stage. |
| `rdrstb` | array of [dvs, rdrs]  | —   | Relative stem death rate vs. development stage. |

#### `[root]`

| Key | Type | Default | Description |
|---|---|---|---|
| `swrd`     | integer             | 0   | 0=fixed `rdc`, 1=tabulated `rdtb`, 2=`rri` rate-driven. |
| `rdi`      | real                | 0.0 | Initial rooting depth (cm). |
| `rri`      | real                | 0.0 | Daily increase of rooting depth (cm/d). |
| `rdc`      | real                | 0.0 | Maximum rooting depth (cm). |
| `swdmi2rd` | integer             | 0   | 1=root growth tied to dry-matter increase. |
| `wrtmax`   | real                | 0.0 | Maximum root weight (kg/ha). |
| `rdtb`     | array of [dvs, rd]  | —   | Tabulated rooting depth (used when `swrd=1`). |
| `rlwtb`    | array of [rd, rlw] | —   | Relative root-water-uptake distribution. |
| `rdctb`    | array of [rd, frac]| —   | Root-density correction vs. depth. |

#### `[oxygen_stress]`

| Key | Type | Default | Description |
|---|---|---|---|
| `swoxygen`                | integer | 0   | 0=off, 1=Feddes (`hlim1/2u/2l`), 2=Bartholomeus root-respiration model. |
| `swwrtnonox`              | integer | 0   | 1=allow root growth into anoxic layers. |
| `aeratecrit`              | real    | 0.0 | Critical aeration fraction (when `swoxygen=2`). |
| `hlim1` / `hlim2u` / `hlim2l` | real | 0.0 | Feddes thresholds (cm, used when `swoxygen=1`). |
| `q10_microbial`           | real    | 0.0 | Q10 of microbial respiration. |
| `specific_resp_humus`     | real    | 0.0 | Specific respiration of humus (mg·O₂·kg⁻¹·d⁻¹). |
| `srl`                     | real    | 0.0 | Specific root length (m/g). |
| `swrootradius`            | integer | 0   | 0=fixed `root_radiusO2`, 1=derive from `srl` & root density. |
| `dry_mat_cont_roots`      | real    | 0.0 | Dry-matter content of roots (g/g). |
| `air_filled_root_por`     | real    | 0.0 | Air-filled root porosity (cm³/cm³). |
| `spec_weight_root_tissue` | real    | 0.0 | Specific weight of dry root tissue (g/cm³). |
| `var_a`                   | real    | 0.0 | Empirical Bartholomeus coefficient. |
| `root_radiusO2`           | real    | 0.0 | Root radius for O₂ diffusion (cm). |

#### `[drought_stress]`

| Key | Type | Default | Description |
|---|---|---|---|
| `swdrought` | integer | 0   | 0=Feddes pressure-head, 1=Jarvis stress factor. |
| `hlim3h`    | real    | 0.0 | Pressure head where stress starts (high transpiration, cm). |
| `hlim3l`    | real    | 0.0 | Pressure head where stress starts (low transpiration, cm). |
| `hlim4`     | real    | 0.0 | Wilting-point pressure head (cm). |
| `adcrh`     | real    | 0.0 | Critical transpiration rate at `hlim3h` (cm/d). |
| `adcrl`     | real    | 0.0 | Critical transpiration rate at `hlim3l` (cm/d). |

#### `[salinity]`

| Key | Type | Default | Description |
|---|---|---|---|
| `swsalinity` | integer | 0   | 0=off, 1=Maas–Hoffman, 2=osmotic-head Feddes extension. |
| `saltmax`    | real    | 0.0 | EC threshold (dS/m) above which yield reduction starts. |
| `saltslope`  | real    | 0.0 | Slope of yield reduction (%/dS/m above `saltmax`). |
| `salthead`   | real    | 0.0 | Osmotic-head conversion factor (cm/(dS/m)) when `swsalinity=2`. |

#### `[compensate]`

Stress-compensation parameters governing root water uptake redistribution.

| Key | Type | Default | Description |
|---|---|---|---|
| `swcompensate` | integer | 0   | 0=no compensation, 1=Jarvis 1989, 2=de Jong van Lier. |
| `swstressor`   | integer | 0   | Bitfield selecting which stressors to compensate. |
| `alphacrit`    | real    | 0.0 | Critical stress factor for compensation onset. |
| `dcritrtz`     | real    | 0.0 | Critical depth in root zone (cm). |

#### `[interception]`

| Key | Type | Default | Description |
|---|---|---|---|
| `swinter` | integer                       | 0   | 0=off, 1=`cofab` × LAI, 2=Gash table model. |
| `cofab`   | real                          | 0.0 | Interception coefficient (cm/LAI; required when `swinter=1`). |
| `gashtb`  | array of [m1,m2,m3,m4,m5,m6] | —   | Gash six-parameter table (required when `swinter=2`). |

> **Reading note.** `gashtb` is only consumed when `swinter=2`; when
> `swinter=1` the field is left unallocated and the validator does not
> require it. The legacy reader populates `gashtb` only on the `swinter=2`
> branch; the TOML reader follows the same convention.

#### `[co2]`

| Key | Type | Default | Description |
|---|---|---|---|
| `swco2`     | integer                | 0   | 0=ambient, 1=read CO₂ history file. |
| `atmofil`   | string                 | ""  | Path to CO₂ atmospheric history file (used when `swco2=1`). |
| `co2amaxtb` | array of [co2, factor] | —   | CO₂ → `amax` correction. |
| `co2efftb`  | array of [co2, factor] | —   | CO₂ → light-use-efficiency correction. |
| `co2tratb`  | array of [co2, factor] | —   | CO₂ → transpiration-efficiency correction. |

#### `[management]`

| Key | Type | Default | Description |
|---|---|---|---|
| `fraharlosorm_lv`     | real    | 0.0 | Harvest loss fraction of leaves left on field. |
| `fraharlosorm_st`     | real    | 0.0 | Harvest loss fraction of stems left on field. |
| `fraharlosorm_so`     | real    | 0.0 | Harvest loss fraction of storage organs left on field. |
| `fradeceasedlvtosoil` | real    | 0.0 | Fraction of dead leaves returned to soil. |
| `swpotrelmf`          | integer | 0   | 0=potential growth, 1=management factor `relmf` applied. |
| `relmf`               | real    | 0.0 | Relative management factor (0..1, used when `swpotrelmf=1`). |

#### Validator constraints (selected)

The `cropwofost_config_t%validate` walks each sub-table; common rules:

- Enum keys (`swprep`, `swsow`, `swgerm`, `swharv`, `swcf`, `idsl`, `swrd`,
  `swoxygen`, `swdrought`, `swsalinity`, `swcompensate`, `swinter`, `swco2`,
  `swpotrelmf`) accept only the values listed above.
- AFGEN tables (`dtsmtb`, `slatb`, `amaxtb`, `frtb`, `fltb`, `fstb`,
  `fotb`, `cftb`, `chtb`, `rdtb`, `rlwtb`, `rdctb`, `rdrrtb`, `rdrstb`,
  `tmpftb`, `tmnftb`, `co2amaxtb`, `co2efftb`, `co2tratb`, `verntb`,
  `rfsetb`) must satisfy `2 ≤ nrows ≤ 15` and the column count noted in
  the encoding section. `gashtb` requires 6 columns.
- Range checks apply only when the relevant switch is on (e.g. `tempsow`
  is only validated when `swsow=1`).
- Cross-field checks: `tsumea > 0` when `idsl ≥ 0`; `dlo ≥ dlc` when
  `idsl ≥ 1`; `vernsat ≥ vernbase` when `idsl = 2`.

#### Cross-file loading

A type-2 entry is selected via the rotation array:

```toml
[[crop.rotation]]
start = 1980-01-01
end   = 1980-12-31
file  = "potatod.crp.toml"
type  = 2
```

`load_swap_config` resolves `file` relative to the `.swp` directory and
calls `read_cropwofost_toml` on the loaded TOML root. The result is stored
under `swap_config%crop%rotation_wofost(i)`. As with type-1 / type-3,
missing files are silently skipped to keep partially-converted cases
loadable.

### `*.crp.toml` — type 3 (WOFOST grass)

Parsed by `read_cropgrass_toml.f90`; config type `cropgrass_config_t`
(`src/config/cropgrass_config.f90`).

Sections: same as type 1 (`[phenology]`, `[light]`, `[root]`, `[water_stress]`,
`[salinity]`, `[interception]`) plus `[mowing]` and `[grazing]`.

#### `[phenology]` (type 3, additional keys)

| Key | Type | Default | Description |
|---|---|---|---|
| `idev`  | integer | 2   | Development mode (2=temperature-sum-based typical for grass). |
| `lcc`   | integer | 0   | Crop cycle length in days (if `idev=1`). |
| `tbase` | real    | 0.0 | Base temperature for development sum (°C). |
| `tsum1` | real    | 0.0 | Temperature sum for vegetative stage (°C·d). |
| `tsum2` | real    | 0.0 | Temperature sum for generative stage (°C·d). |

#### `[light]` (type 3, additional keys)

| Key | Type | Default | Description |
|---|---|---|---|
| `kdif` | real | 0.0 | Diffuse light extinction coefficient. |
| `kdir` | real | 0.0 | Direct light extinction coefficient. |
| `eff`  | real | 0.0 | Light use efficiency (kg·ha⁻¹·h⁻¹/(J·m⁻²·s⁻¹)). |
| `amax` | real | 0.0 | Maximum assimilation rate (kg·ha⁻¹·h⁻¹). |

#### `[mowing]` (type 3 only)

| Key | Type | Default | Description |
|---|---|---|---|
| `swharv` | integer | 0 | 0=no scheduled mowing, 1=scheduled. |
| `nmow`   | integer | 0 | Number of mowing events (required when `swharv=1`). |

#### `[grazing]` (type 3 only)

| Key | Type | Default | Description |
|---|---|---|---|
| `swgraz`      | integer | 0   | 0=no grazing, 1=scheduled. |
| `nstart_graz` | real    | 0.0 | Start day-of-year for grazing. |
| `nstop_graz`  | real    | 0.0 | Stop day-of-year for grazing. |

## Example

The smallest `.swp` that exercises the supported keys looks like this. Every
line except `start_date` and `end_date` is optional; sections may be omitted
entirely and the config field keeps its default value.

```toml
# Minimal SWAP TOML configuration (Phase 4c-a)

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
nprintday  = 1

[simulation.output]
swmonth = 1
period  = 1
swres   = 0
swodat  = 0

[meteorology]
file = "283.met"
lat  = 52.0
alt  = 10.0
altw = 10.0

[meteorology.evapotranspiration]
swetr      = 0
angstrom_a = 0.25
angstrom_b = 0.5
swdivide   = 1

[meteorology.temporal]
swmetdetail = 0
nmetdetail  = 24
swetsine    = 0

[meteorology.rainfall]
swrain = 0

[crop]
swcrop = 1

[[crop.rotation]]
start = 2002-05-01
end   = 2002-10-15
file  = "maizes.crp.toml"
type  = 1

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
file = "swap.dra.toml"

[boundary.bottom]
swbotb = 6

[boundary.top]
rsro = 0.5
```

A fuller example lives at `tests/swap-cases/1.1.hupselbrook-toml/swap.toml`.
That file contains many additional keys not yet read by the TOML pipeline
(irrigation, heat, solute, soil hydraulic tables, macropores, snow, frost);
those keys are silently ignored and are targets for later rescue phases.

## Cross-file references

Phase 4c-a introduces explicit `file = "..."` references that point at
external TOML files for drainage and per-crop configuration. The loader
follows these references relative to the directory of the referencing file
(via `path_helpers_mod%resolve_relative_path`).

### Drainage

```toml
# in swap.toml
[drainage]
file = "swap.dra.toml"   # loader follows this; swap.dra.toml has its own
                         # [drainage] table with all the drainage keys
```

`swap.dra.toml` itself uses the same `[drainage]` schema as the inline form:

```toml
# swap.dra.toml
[drainage]
swdra  = 1
dramet = 3
# ... etc.
```

The loader calls `read_drainage_toml` with `base_path=` set to the directory
of the `.swp` file. If `file` is absent the inline `[drainage]` table is
parsed instead.

### Per-crop files

Each `[[crop.rotation]]` entry references its `.crp.toml` file via the `file`
key:

```toml
[[crop.rotation]]
start = 2002-05-01
end   = 2002-09-30
file  = "maizes.crp.toml"
type  = 1
```

The `.crp.toml` schema depends on `type`:
- `type = 1` → fixed crop schema (see "type 1 (fixed crop)" above)
- `type = 2` → WOFOST general (Phase 4c-b; placeholder accepted in 4c-a)
- `type = 3` → WOFOST grass (see "type 3 (grass)" above)

The loader uses `[[crop.rotation]].type` to dispatch to the appropriate
section reader. Phase 4c-a handles types 1 and 3; type 2 entries are
populated only as placeholders.

If `file` is absent for a rotation entry, or if the referenced file does not
exist on disk, the entry is skipped silently — this allows incremental
migration of cases without breaking cases that have not yet been converted.

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
5. The Phase-2 aspirational schema used string arrays under `[crop.rotation]`.
   Phase 4c-a uses `[[crop.rotation]]` array-of-tables with actual TOML date
   values. The `parse_date_to_days1900` helper is shared with the simulation
   reader so both paths now use real TOML date parsing, eliminating the prior
   inconsistency.
6. `[[crop.rotation]].file` is silently skipped when the referenced file does
   not exist on disk. This is deliberate: it allows test cases that have not
   yet been converted to `.crp.toml` to continue loading without error. The
   consequence is that missing files produce no diagnostic — a strict mode
   that warns on absent references is deferred to Phase 4d.

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
| `*.crp.toml` type 1 | `cropfixed_config_t` | `cropfixed_config_mod` |
| `*.crp.toml` type 3 | `cropgrass_config_t` | `cropgrass_config_mod` |
| (top-level) | `swap_config_t` | `swap_config_mod` |

Every type exposes `validate(errors)` and `finalize(errors)` as
type-bound procedures. See `docs/validation.md` for the rules and
`docs/error-handling.md` for how errors propagate.
