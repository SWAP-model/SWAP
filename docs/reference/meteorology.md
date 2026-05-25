---
title: Meteorology
---

# Meteorology

The `[meteorology]` section of `swap.toml` controls weather input and the
reference evapotranspiration method. This page covers the available options,
their TOML keys, and the CSV file formats expected for each.

## Quick reference

```toml
[meteorology]
file = "283.csv"        # daily meteo CSV (ISO date)
lat  = 52.0             # latitude (degrees N)
alt  = 30.0             # altitude (m)
altw = 2.0              # anemometer height (m)

[meteorology.evapotranspiration]
swetr    = 1            # 0 = Penman-Monteith (requires full meteo); 1 = ETref from file
swdivide = 0            # 0 = use Penman; 1 = divide ETref into Ep/Et by crop coefficient

[meteorology.temporal]
swmetdetail = 0         # 0 = daily; 1 = sub-daily (not yet wired to CSV path)
swmetfilall = 0         # legacy flag; set automatically from file extension

[meteorology.rain]
swrain       = 0        # 0 = from meteo file; 1 = uniform; 2 = ETref-based; 3 = events CSV
swetsine     = 0        # 0 = uniform within day; 1 = sine distribution
events_file  = ""       # required when swrain = 3

[meteorology.interception]
swinter = 0             # 0 = no interception; 1 = von Hoyningen-Huene; 2 = Gash

[meteorology.evaporation]
swcfbs   = 0            # 0 = fixed reduction; 1 = Black/Boesten-Stroosnijder model
cofredbl = 0.35         # Black reduction coefficient (-)
cofredbo = 0.35         # Boesten-Stroosnijder reduction coefficient (-)

[meteorology.snow]
swsnow   = 0            # 0 = no snow; 1 = degree-day snow model
snowcoef = 0.0          # degree-day melt coefficient (cm °C⁻¹ d⁻¹)
teprrain = 1.0          # temperature threshold: rain vs snow (°C)
teprsnow = -1.0         # temperature threshold: snow accumulation (°C)
```

---

## Daily meteorology CSV (`file`)

Used when `[meteorology].file` ends in `.csv`. The file must cover every
year of the simulation period.

**Schema:**

```
date,rad,tmin,tmax,hum,wind,rain,etref,wet
```

| Column  | Unit         | Description                                          |
|---------|--------------|------------------------------------------------------|
| `date`  | ISO date     | Calendar date `YYYY-MM-DD`                           |
| `rad`   | kJ m⁻² d⁻¹  | Global radiation                                     |
| `tmin`  | °C           | Daily minimum temperature                            |
| `tmax`  | °C           | Daily maximum temperature                            |
| `hum`   | kPa          | Actual vapour pressure                               |
| `wind`  | m s⁻¹        | Mean wind speed at anemometer height (`altw`)        |
| `rain`  | mm d⁻¹       | Daily precipitation (ignored when `swrain = 3`)      |
| `etref` | mm d⁻¹       | Reference evapotranspiration (required if `swetr=1`) |
| `wet`   | d d⁻¹        | Fraction of day with precipitation (required if `swrain=2`) |

**Example:**

```csv
# Station: Hupsel (283), source: KNMI
# Units: rad kJ/m2/d, temp °C, hum kPa, wind m/s, rain mm, etref mm, wet d/d
date,rad,tmin,tmax,hum,wind,rain,etref,wet
2002-01-01,3810.0,-3.2,-0.1,0.524,4.90,0.000,0.400,0.000
2002-01-02,1200.0, 0.1, 4.3,0.710,5.20,3.500,0.180,0.320
2002-01-03, 580.0, 2.5, 6.1,0.780,8.40,8.200,0.100,0.620
```

**Notes:**

- The reader converts `rad` from kJ m⁻² d⁻¹ to J m⁻² d⁻¹ internally; author
  the file in kJ (matching the legacy `.met` format).
- Missing values are not supported. All rows must contain all 9 columns.
- Rows are sorted by date ascending. Gaps between years are not allowed.

---

## Rain events CSV (`rain.events_file`)

Used when `swrain = 3`. Each row is a sub-daily rainfall event.  The adapter
pre-loads all events at startup; `ReadRainEvents` extracts the current year's
slice at the start of each simulation year.

**Schema:**

```
datetime,amount
```

| Column     | Unit         | Description                                               |
|------------|--------------|-----------------------------------------------------------|
| `datetime` | ISO datetime | Event timestamp `YYYY-MM-DD HH:MM:SS`                    |
| `amount`   | mm           | Rainfall amount for this interval                         |

**Example:**

```csv
# Sub-daily rainfall events, station Andelst
# datetime in local time; amount in mm per interval
datetime,amount
1998-01-01 00:00:00,0.0
1998-01-01 06:00:00,2.4
1998-01-01 12:00:00,5.1
1998-01-01 18:00:00,0.0
1998-01-02 00:00:00,0.0
1998-01-02 08:30:00,1.2
```

**Notes:**

- The `datetime` column is parsed to fractional days since Julian Day 2415020
  (1899-12-31), the same epoch used throughout SWAP.
- If the first event of a year is not at midnight (`HH:MM:SS = 00:00:00`), a
  zero-amount record at midnight is automatically prepended.
- Events must be in strictly ascending time order across the whole file.
- The file must cover every year in the simulation range.
- The `rain` column in the daily meteo file is ignored when `swrain = 3`.

**TOML configuration:**

```toml
[meteorology.rain]
swrain      = 3
events_file = "mysite.rain.csv"
```

---

## Sub-daily detail meteorology CSV (`temporal.detail_file`)

Used when `[meteorology.temporal].swmetdetail = 1`. The file must cover every
year of the simulation period. Each row is one sub-daily time slot.

**Schema:**

```
datetime,record,rad,temp,hum,wind,rain
```

| Column     | Unit         | Description                                                    |
|------------|--------------|----------------------------------------------------------------|
| `datetime` | ISO datetime | Slot timestamp `YYYY-MM-DD HH:MM:SS`; fractional days since JD 2415020 |
| `record`   | –            | Intra-day slot index 1 … `nmetdetail`                          |
| `rad`      | kJ m⁻² d⁻¹  | Global radiation for the slot (converted to J m⁻² d⁻¹ internally) |
| `temp`     | °C           | Air temperature (single value per slot)                        |
| `hum`      | kPa          | Actual vapour pressure                                         |
| `wind`     | m s⁻¹        | Wind speed                                                     |
| `rain`     | mm           | Rainfall for the slot                                          |

**Example (48 slots per day):**

```csv
# Sub-daily meteo, station Hupsel, 48 half-hourly slots per day
datetime,record,rad,temp,hum,wind,rain
2002-01-01 00:00:00,1,0.0,-3.2,0.524,4.90,0.000
2002-01-01 00:30:00,2,0.0,-3.3,0.525,4.85,0.000
2002-01-01 23:30:00,48,0.0,-2.9,0.520,5.10,0.000
```

**TOML configuration:**

```toml
[meteorology]
file = "hupsel.csv"

[meteorology.temporal]
swmetdetail = 1
nmetdetail  = 48
detail_file = "hupsel.det.csv"
```

**Notes:**

- The adapter (`config_to_variables.f90`) pre-loads all rows into `metcsv_det`
  at startup. `MeteoCSVDetYear` extracts the current year's slice on each
  `ReadMeteoYear` call, exactly mirroring the daily `MeteoCSVYear` pattern.
- The maximum record count per year is `NMETFILE = 17568` (48 slots × 366 days).
- `irectotal` is initialised in `ReadMeteoYear` using `dettime(1)` after
  `MeteoCSVDetYear` returns — no changes needed in `meteoday.f90`.
- Year boundaries are half-open `[t_jan1, t_jan1_next)` to correctly bucket
  `YYYY-12-31 23:59:59` into year YYYY (not YYYY+1).
- If no rows are found for the requested year, SWAP aborts via
  `fatalerr_collected`.

---

## Legacy formats (ASCII `.swp` pathway only)

The following formats are used by the legacy `swap420` binary and the ASCII
`.swp` pathway. They are **not** read by the TOML pipeline and should not be
placed in the TOML case directory.

### Per-year files (`<stem>.<YYY>`)

One file per calendar year; extension is the last three digits of the year
(e.g. `.002` for 2002, `.998` for 1998).

**Daily meteo schema** (keyword-value, TTutil format):

```
 Station  DD MM YYYY  RAD    Tmin   Tmax   HUM    WIND  RAIN  ETref  wet
'ST01'    1  1  2002  3810.0 -3.2  -0.1   0.524  4.90  0.000 0.400  0.000
```

**Rain events schema** (keyword-value, TTutil format):

```
 day  month  year  time   amount
  1     1    1998  0.000   0.0
  1     1    1998  0.250   2.4
```

Where `time` is the fraction of the day (0.0 = midnight, 0.5 = noon).

### All-years file (`<stem>.met`)

Same column layout as the per-year format but all years concatenated into one
file. Activated automatically when `file` ends in `.met` and `swmetdetail = 0`.

---

## Choosing the right option

| Scenario                                  | Recommended setting         |
|-------------------------------------------|-----------------------------|
| TOML pipeline, daily rain from meteo file | `swrain = 0`, `file = "*.csv"` |
| TOML pipeline, sub-daily rain events      | `swrain = 3`, `events_file = "*.rain.csv"` |
| Legacy `.swp` pipeline, per-year files    | `swrain = 0`, legacy `.YYY` files |
| Legacy `.swp` pipeline, all-years file    | `swmetfilall = 1`, `*.met` file |
