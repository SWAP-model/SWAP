---
title: SWAP TOML configuration format guide
author: SWAP modernization team
---

# SWAP TOML configuration format guide

A practical guide to writing SWAP project configurations in TOML format.
For a full key-by-key reference, see `configuration-schema.md`.

## Overview

SWAP project configuration is moving from the legacy ttutil-based
key=value `.swp` / `.dra` / `.crp` files to a structured TOML format.
The new format is split into thematic sections, validated through a
typed pipeline (`load_swap_config` → `validate` → `finalize`), and
loaded into per-section `*_config_t` types.

The pipeline lives in `src/io/toml/`: a composite reader at
`load_swap_config.f90` dispatches to per-section readers
(`read_general_toml.f90`, `read_meteorology_toml.f90`, etc.), and a
round-trip writer is provided at `write_swap_config.f90`.

This guide explains the format conventions and walks through a real
example. It does not enumerate every key — see `configuration-schema.md`
for that.

## File organization

A SWAP simulation reads a primary file named `swap.toml` (or any
`.swp`-style name resolved by the runner). Every other path is resolved
relative to that file via the `[general.paths]` subsection.

The format supports splitting configuration across multiple files using
the `file = "..."` convention inside a section. Current status:

- **`swap.toml`** — main file, all sections supported as of Phase 4b.
- **`swap.dra.toml`** — drainage detail file; cross-file loading is
  Phase 4c work. For now, drainage fields are inlined directly in
  `swap.toml`.
- **`*.crp.toml`** — per-crop detail file; TOML crop reader is Phase 4c.
  The `[[crop.rotation]]` entries store `.crp` filenames as strings;
  crop detail is still read by the legacy fixed-format readers in
  `src/crop/`.

In practice, every currently-shipped `swap.toml` is self-contained: all
drainage configuration lives inline under `[drainage]` and
`[[drainage.levels]]`, not in a separate `.dra.toml`.

## Sections

Six top-level sections cover the full Phase 4b schema. Each maps to a
`*_config_t` type in `src/config/`:

- **`[general]`** — project name, screen-output and error switches, and
  the `[general.paths]` subsection that sets directories for output,
  atmosphere, crop, and drainage files.
- **`[simulation]`** — simulation period (`start_date` / `end_date`),
  output frequency, and the `[simulation.output]` subsection for
  monthly / fixed-period output control.
- **`[meteorology]`** — met-file name, station latitude and altitude, and
  three subsections: `[meteorology.evapotranspiration]` (ET method,
  Angstrom coefficients), `[meteorology.rain]` (rainfall input mode),
  and `[meteorology.interception]` (interception switch).
- **`[drainage]`** — drainage mode switch (`swdra`), drainage method
  (`dramet`), profile-distribution flag (`swdivd`), and per-level data
  in `[[drainage.levels]]` array-of-tables.
- **`[soil]`** — soil physics switches (`swsophy`, `swhyst`, `swinco`,
  `swmacro`) and the `[soil.initial]` subsection for initial groundwater
  level and ponding.
- **`[crop]`** — crop simulation switch (`swcrop`) and per-season entries
  in `[[crop.rotation]]`.

See `configuration-schema.md` for the complete key listing per section.

## Type and value conventions

### Scalars

**Integer keys (switches)** are small enums, typically 0–3. TOML
requires a bare integer literal — no decimal point:

```toml
swetr    = 0
swrain   = 2
dramet   = 3
```

Section validators enforce valid membership (e.g. `dramet ∈ {1,2,3}`).
A float literal where an integer is expected (`swetr = 0.0`) will be
rejected by toml-f at parse time.

**Real keys** must always carry a decimal point, even when the value is
whole-number-like. The Fortran reader expects a real type:

```toml
lat        = 52.0
angstrom_a = 0.25
basegw     = -200.0
```

**String keys** are always double-quoted:

```toml
project = "hupsel"
file    = "283.met"
```

**Date keys** use TOML local-date literals — no quotes, no time
component. The reader converts them internally to days since 1900:

```toml
start_date = 2002-01-01
end_date   = 2004-12-31
```

This applies to simulation start/end and to `[[crop.rotation]]` entry
dates. The two conventions (TOML date vs. quoted ISO string) existed
side-by-side in earlier readers; the Phase 4b readers use TOML-native
dates throughout.

### Arrays

**Inline arrays** are used for fixed-length scalar sequences. All
elements must be the same type:

```toml
swtopdislay = [0, 0, 1]
ztopdislay  = [-30.0, -60.0, -90.0]
```

An array mixing integers and floats (`[1, 2.0, 3]`) is a TOML error.

**Array-of-tables** use `[[double.brackets]]` and are the preferred form
when each entry has its own named keys. SWAP uses this for drainage
levels and crop rotation seasons:

```toml
[[drainage.levels]]
swdtyp = 2
zbotdr = -55.0
drares = 750.0
L      = 500.0

[[drainage.levels]]
swdtyp = 2
zbotdr = -20.0
drares = 50.0
L      = 15.0
```

Each `[[drainage.levels]]` block appends a new element; the array is
built in file order.

### Tables (sections and subsections)

Top-level sections use a single bracket:

```toml
[general]
[simulation]
[drainage]
```

Nested subsections use dotted headers:

```toml
[meteorology.evapotranspiration]
[soil.initial]
[simulation.output]
```

Section order within the file does not matter to TOML. By convention,
SWAP files list sections top-to-bottom in the order they appear in this
guide (general → simulation → meteorology → drainage → soil → crop).

## Validation

The pipeline applies two layers of validation. See `validation.md` for
the full rule catalogue and `error-handling.md` for how errors propagate.

**Per-section validation** runs via `*_config_t%validate(errors)` after
the section is populated. Examples:

- `meteorology_config_t` checks `swetr ∈ {0,1}`, `swrain ∈ {0,1,2}`,
  `lat ∈ [-90, 90]`.
- `drainage_config_t` checks `dramet ∈ {1,2,3}` and enforces the
  cross-field rule that `dramet=2` (Hooghoudt/Ernst) requires `swdivd=1`.

**Aggregate validation** runs via `swap_config_t%validate(errors)`,
which delegates to every section validator and then checks cross-section
invariants (e.g. that `start_date < end_date`).

An example of a cross-field validation error: if you write

```toml
[drainage]
dramet = 2
swdivd = 0   # wrong — Hooghoudt/Ernst requires swdivd = 1
```

`drainage_config_t%validate` appends a fatal `ERR_VALIDATION_CROSS_FIELD`
error with context `"drainage"` and the message
`"swdivd must be 1 when dramet=2"`.

## Error reporting

Errors accumulate in an `error_collection_t`. The caller checks for
fatal errors after the full pipeline:

```fortran
call load_swap_config(path, config, errors)
call config%validate(errors)
call config%finalize(errors)
call errors%abort_if_fatal()
```

`abort_if_fatal` writes `errors%summary()` to stderr and calls
`error stop` if any appended error has `is_fatal=.true.`. Non-fatal
warnings pass through without aborting. See `error-handling.md` for the
full contract, error codes, and testing patterns.

Unknown keys in a TOML file are currently **silently ignored** by the
readers — toml-f does not raise on them. A strict-validation pass
(Phase 4c) will reject unknown keys; do not rely on silently-ignored
keys as an extension point.

Wrong types (integer literal where a real is expected, or vice versa)
are rejected by toml-f itself, but because the readers pass `stat=istat`,
a single mis-typed key is silently skipped rather than aborting the run.
Pay attention to the scalar type conventions above.

## Worked example: hupselbrook swap.toml

The most complete reference file is
`tests/swap-cases/toml/1.hupselbrook/swap.toml`. The annotated version
below covers the sections most likely to need explanation.

```toml
# SWAP hupselbrook project

[general]
project = "hupsel"
swscre  = 0       # 0 = minimal screen output
swerror = 0       # 0 = non-fatal errors do not abort

[general.paths]
work       = "./"  # directory for output files
atmosphere = "./"  # directory for .met file
crop       = "./"  # directory for .crp files
drain      = "./"  # directory for drainage files

[simulation]
start_date = 2002-01-01   # TOML local-date — no quotes
end_date   = 2004-12-31
nprintday  = 1

[simulation.output]
swmonth = 1   # 1 = monthly output (ignores period below)
period  = 1   # output interval in days when swmonth=0
swres   = 0
swodat  = 0
swyrvar = 0

[meteorology]
file = "283.met"   # resolved against [general.paths].atmosphere
lat  = 52.0
alt  = 10.0
altw = 10.0        # anemometer height (m)

[meteorology.evapotranspiration]
swetr      = 0     # 0 = Penman-Monteith; 1 = reference ET with crop factors
swdivide   = 1     # 1 = direct Penman-Monteith partitioning
angstrom_a = 0.25
angstrom_b = 0.50

[meteorology.rain]
swrain   = 0   # 0 = daily totals
swetsine = 0

[meteorology.interception]
swinter = 0

# dramet=2 (Hooghoudt/Ernst) requires swdivd=1.
# The validator enforces this; mismatching values produce a fatal error.
[drainage]
swdra    = 1   # 1 = basic drainage
dramet   = 2   # 2 = Hooghoudt/Ernst formula
swdivd   = 1   # must be 1 when dramet=2
swdislay = 0
nrlevs   = 0   # 0 = single implicit level for dramet=2
altcu    = 0.0

[drainage.basic]
basegw = -200.0   # depth of impervious layer (cm)
entres = 20.0     # drain entry resistance (d)
shape  = 0.8      # shape factor

[soil]
swsophy = 0   # 0 = use built-in hydraulic functions
swhyst  = 0   # 0 = no hysteresis
swinco  = 2   # 2 = init from gwli below
swmacro = 0   # 0 = no macropore flow

[soil.initial]
gwli    = -75.0   # initial groundwater level (cm, negative = below surface)
pondini = 0.0
pondmx  = 0.2     # minimum ponding depth before runoff (cm)

[crop]
swcrop = 1   # 1 = simulate crop

# Each [[crop.rotation]] block is one growing season.
# Blocks are ordered by start date; they must not overlap.
[[crop.rotation]]
start = 2002-05-01
end   = 2002-10-15
file  = "maizes"   # .crp filename without extension
type  = 1          # crop type code (1=simple, 2=detailed, 3=grass)

[[crop.rotation]]
start = 2003-05-10
end   = 2003-09-29
file  = "potatod"
type  = 2

[[crop.rotation]]
start = 2004-01-01
end   = 2004-12-31
file  = "grassd"
type  = 3
```

A few things worth noting in this example:

- `dramet=2` (Hooghoudt/Ernst) requires `swdivd=1`. This is a cross-field
  rule enforced by `drainage_config_t%validate`. Omitting `swdivd` or
  setting it to 0 produces a fatal validation error.
- `start_date` and `end_date` (in both `[simulation]` and
  `[[crop.rotation]]`) are TOML local-date literals. No quotes, no time
  component. The reader converts them to days since 1900.
- `[[crop.rotation]]` uses array-of-tables syntax. Each double-bracketed
  block appends a new season. The alternative (parallel arrays of
  strings) was used in an earlier reader version and is shown in the
  Discoveries section of `configuration-schema.md`; the Phase 4a/4b
  readers use the array-of-tables form.
- `nrlevs = 0` under `[drainage]` is correct for `dramet=2`; the
  per-level data lives in `[drainage.basic]` instead. Cases using
  `dramet=3` (resistance formula) set `nrlevs ≥ 1` and provide one
  `[[drainage.levels]]` block per level.

### Multi-level drainage example (dramet=3)

Cases 2–5 use the resistance-formula drainage method with one or two
levels. The pattern is:

```toml
[drainage]
swdra    = 1
dramet   = 3
swdivd   = 1
swdislay = 0
nrlevs   = 2    # two drainage levels
altcu    = 0.0

[[drainage.levels]]
swdtyp = 2       # 2 = open channel
zbotdr = -140.0
drares = 130.0
infres = 150.0
L      = 40.0

[[drainage.levels]]
swdtyp = 2
zbotdr = -20.0
drares = 50.0
infres = 40.0
L      = 15.0
```

`nrlevs` declares how many levels to expect; the reader processes that
many `[[drainage.levels]]` blocks. Blocks beyond `nrlevs` are silently
ignored.

## Naming conventions

- **Section names** are lowercase, no separators. Nested sections use
  a dot: `[meteorology.evapotranspiration]`, `[soil.initial]`.
- **Key names** are lowercase with underscores:
  `angstrom_a`, `start_date`, `nprintday`.
- **Legacy name mapping**: keys map as closely as possible to the
  corresponding legacy SWAP variable names, with underscores inserted
  for readability. For example:
  - `ANGSTROMA` → `angstrom_a`
  - `TSTART` / `TEND` → `start_date` / `end_date`
  - `METFIL` → `file` (inside `[meteorology]`)
  - `PATHWORK` / `PATHATM` / ... → `work` / `atmosphere` / ... (inside
    `[general.paths]`)
- **Switch keys** keep the legacy `sw` prefix: `swscre`, `swetr`,
  `swrain`, `swcrop`, `swdra`, `swmacro`.

## Known gaps — Phase 4c

The following features from the legacy `.swp` format are not yet covered
by the TOML schema:

- **Bottom boundary** (`SWBOTB` and related keys) — not yet in any
  `*_config_t`. Cases 4 and 5 stub this.
- **Heat transport** (`SWHEA`) — not yet in the schema.
- **Solute transport** (`SWSOLU`) — not yet in the schema. Case 5 stubs
  this.
- **Extended drainage / surface-water management** — `swdra=2` is
  accepted by `drainage_config_t` but the extended-drainage fields
  (management periods, weir tables, `NRSRF`, etc.) are not yet
  populated. Case 6 stubs this.
- **Irrigation** (`SWIRFIX` and schedule) — not yet in the schema.
- **`.crp.toml` cross-file loading** — rotation entries store `.crp`
  filenames as strings; the Phase 4c crop reader will consume them.
- **`.dra.toml` cross-file loading** — drainage section is inlined in
  `swap.toml` for now.
- **Detailed rainfall file** (`swrain=3`, `RAINFIL`) — `swrain` range
  is validated as 0–2 in the current readers; `swrain=3` is a Phase 4c
  extension.

See `archive/2026-phase-4/audits/phase-4b-field-coverage.md` for a detailed per-case gap inventory.

## Related docs

- `configuration-schema.md` — full key reference per section
- `error-handling.md` — error model, abort checkpoint, testing patterns
- `validation.md` — primitive checks, how to write section validators
- `logging.md` — logger usage
- `architecture.md` — overall pipeline and module boundaries
