---
title: TOML config
author: Mateusz Zawadzki
---

# New TOML Configuration Files

In the refactored SWAP, configuration is read from TOML files, phasing out the legacy ASCII 
.swp, .crp, .dra formats. TOML (Tom's Obvious Minimal Language) is a widely recognized syntax 
for configuration files used in software development and modeling applications (e.g., MODFLOW 6). 
SWAP leverages the toml-f library to parse these files, providing type-safe, validated input data.

## Why TOML?

- **Human-readable**: Clear, intuitive syntax
- **Type-safe**: Explicit data types (integers, floats, strings, booleans, dates)
- **Hierarchical**: Natural nesting of related parameters
- **Standardized**: Well-defined specification with mature parsers
- **Validation-friendly**: Schema validation tools available
- **Comment support**: Extensive documentation within config files
## Novelties
@note
The list of switches for the output files has been removed. Now the overarching parameter is output and it is a list of files (extensions) that user wants to get:

INSTEAD OF:
```toml
# Control which output files are generated
swwba = 0     # Cumulative water balance
swend = 0     # End conditions
swvap = 1     # Soil profiles (moisture, solute, temperature)
swbal = 0     # Yearly water balance
swblc = 1     # Detailed yearly water balance
swsba = 1     # Cumulative solute balance
swate = 0     # Soil temperature profiles
swbma = 0     # Water fluxes (macropore flow only)
swdrf = 0     # Drainage fluxes (extended drainage only)
swswb = 0     # Surface water reservoir (extended drainage only)
swini = 0     # Initial soil physical and heat parameters
swinc = 1     # Water balance increments
swcrp = 0     # Crop growth output
swstr = 0     # Stress values (wetness, drought, salinity, frost)
swirg = 0     # Irrigation gifts
swcsv = 1     # Specific CSV output file
swcsv_tz = 0  # Specific CSV output file with depth information
```

USER PASSES:
```toml
output = ['vap', 'blc', 'sba', 'inc', 'csv']
```
For now parsing happens in the fortran code, but all those files will be slated for removal in favour of csv output.
@endnote
## File Organization

Configuration is split across multiple files for modularity:

- `swp.toml` - Main simulation configuration (meteorology, soil, boundaries, output)
- `crp.toml` - Crop-specific parameters (separate file per crop type)
- `dra.toml` - Drainage system configuration (if extended drainage is used)
- `bbc.toml` - Bottom boundary conditions (DEPRECATED; should be slated for removal in favour of boundary conditions in the `swp.toml` and a separate `csv` file with any additional data)

This guide describes the structural conventions and best practices for organizing TOML 
configuration files in SWAP.

---

## Sections

Sections organize related parameters into logical groups using square brackets `[section_name]`.
Sections can be nested using dot notation to create hierarchies.

### Top-Level Sections

Top-level sections typically correspond to major model domains or functional areas:

```toml
[general]
project = "hupsel"

[simulation]
start_date = 2002-01-01T00:00:00
end_date = 2004-12-31T23:59:59

[atmosphere]
file = "meteo.met"
lat = 52.0

[soil]
# Soil configuration...

[boundary]
# Boundary conditions...
```

**Convention**: Top-level sections align with SWAP's refactored state modules:
- `[atmosphere]` → `atmosphere_state_t`
- `[soil]` → `soil_state_t`
- `[boundary]` → `boundary_state_t`
- `[crop]` → `crop_state_t`

### Nested Sections (Subsections)

Related parameters are grouped under parent sections using dotted keys:

```toml
[atmosphere]
file = "meteo.met"
lat = 52.0

[atmosphere.rainfall]
swrain = 1  # Rainfall mode: 0=daily, 1=+intensity, 2=+duration, 3=file
file = "wagrain"

[soil.initial]
swinco = 2  # Initial condition mode: 1=pressure profile, 2=hydrostatic, 3=file
gwli = -75.0  # Initial groundwater level [cm]

[soil.surface]
pondmx = 0.2  # Maximum ponding depth [cm]
rsro = 0.5    # Runoff resistance [d]
```

**Nesting Depth**: Typically 2-3 levels deep. Avoid excessive nesting (>3 levels) for readability.

### Tables Subsection Pattern

Tabular data (arrays of related parameters) are grouped under a `.tables` subsection:

```toml
[atmosphere]
file = "meteo.met"
lat = 52.0

# Tabular data grouped under .tables
[atmosphere.tables.rainfall_intensity]
day = [1.0, 90.0, 180.0, 360.0]    # Day of year
rate = [15.0, 25.0, 20.0, 15.0]    # Intensity [mm/d]

[soil.hydraulic]
swsophy = 0  # 0=analytical, 1=tables

[soil.tables.hydraulic_parameters]
ores = [0.01, 0.02]              # Residual water content per layer
osat = [0.42, 0.38]              # Saturated water content per layer
alfa = [0.0276, 0.0213]          # Van Genuchten alpha [1/cm]
npar = [1.491, 1.951]            # Van Genuchten n [-]
```

**Benefits**:
- Clear visual separation of scalar parameters vs. tabular data
- Easy navigation: all tables in one subsection
- Scalable: easy to add more tables without cluttering parent section

### Section Access in Fortran

```fortran
use toml_f, only: toml_table, get_value

type(toml_table), pointer :: doc, atm_tab, rain_tab, tables_tab, intensity_tab

! Navigate to nested section
call get_value(doc, 'atmosphere', atm_tab, stat=istat)
call get_value(atm_tab, 'rainfall', rain_tab, stat=istat)
call get_value(rain_tab, 'swrain', state%atm%swrain)

! Access tables subsection
call get_value(atm_tab, 'tables', tables_tab, stat=istat)
call get_value(tables_tab, 'rainfall_intensity', intensity_tab, stat=istat)
call get_value(intensity_tab, 'day', day_array)
call get_value(intensity_tab, 'rate', rate_array)
```

---

## Scalar Values

Scalar values are single parameters specified as key-value pairs. TOML supports multiple data types:

### Integers

```toml
[output]
nprintday = 1       # Number of output times per day
swheader = 0        # Print header: 0=no, 1=yes

[soil.numerical]
maxit = 30          # Maximum iteration cycles
maxbacktr = 3       # Maximum backtrack cycles
```

**Fortran reading**:
```fortran
integer :: nprintday, swheader
call get_value(tab, 'nprintday', nprintday, default=1)
call get_value(tab, 'swheader', swheader, default=0)
```

### Floating-Point Numbers

```toml
[atmosphere]
lat = 52.0            # Latitude [degrees]
alt = 10.0            # Altitude [m]
angstrom_a = 0.25     # Angstrom coefficient [-]

[soil.numerical]
dtmin = 1.0e-6        # Minimum time step [d]
dtmax = 0.04          # Maximum time step [d]
gwlconv = 100.0       # GWL convergence criterion [cm]
```

**Fortran reading**:
```fortran
real(8) :: lat, dtmin
call get_value(tab, 'lat', lat, default=52.0d0)
call get_value(tab, 'dtmin', dtmin, default=1.0e-6d0)
```

### Strings

```toml
[general]
project = "hupsel"

[atmosphere]
file = "meteo.met"

[general.paths]
work = "./"
atmosphere = "/data/meteo/"
crop = "/data/crops/"
```

**Fortran reading**:
```fortran
character(len=:), allocatable :: s
call get_value(tab, 'project', s, stat=istat)
if (istat == 0) state%project = trim(s)
```

### Booleans

```toml
[crop]
flcropnut = false     # Enable nutrient simulation

[soil.numerical]
fldumpconvcrit = false   # Output convergence warnings
swcaprise = true         # Reduce capillary rise below root zone
```

**Fortran reading**:
```fortran
logical :: flcropnut
call get_value(tab, 'flcropnut', flcropnut, default=.false.)
```

### Dates and Times

TOML has native date/time support using ISO 8601 format:

```toml
[simulation]
start_date = 2002-01-01T00:00:00
end_date = 2004-12-31T23:59:59

[[irrigation.schedule]]
date = 2002-06-01        # Date only (implicit 00:00:00)
depth = 10.0
```

**Fortran reading**:
```fortran
type(toml_datetime) :: dt
integer :: year, month, day
call get_value(tab, 'start_date', dt, stat=istat)
year = dt%date%year
month = dt%date%month
day = dt%date%day
! Convert to SWAP's internal time representation (e.g., day since 1900)
```

### Comments and Documentation

Use `#` for inline documentation:

```toml
[soil.evaporation]
swredu = 1         # Reduction method: 0=Darcy, 1=Darcy+Black, 2=Darcy+Boesten
cofredbl = 0.35    # Black coefficient [cm/d^0.5, range: 0..1]
rsigni = 0.5       # Min rainfall to reset reduction [cm/d]
```

**Best practices**:
- Include units in comments: `[cm]`, `[m/d]`, `[-]` (dimensionless)
- Document valid ranges: `[0..1]`, `[-90..90]`
- Explain switch values: `# 0=option_a, 1=option_b, 2=option_c`

---

## Tables

Tables represent multi-column tabular data (e.g., time series, depth profiles, layer parameters).
In SWAP, tables use **column-wise arrays** where each column is a separate key.

### Single Table with Multiple Columns

Use single square brackets `[section.tables.table_name]` for one table per type:

```toml
[atmosphere.tables.rainfall_intensity]
day = [1.0, 90.0, 180.0, 360.0]       # Day of year
rate = [15.0, 25.0, 20.0, 15.0]       # Rainfall intensity [mm/d]

[boundary.tables.prescribed_flux]
date = [2002-01-01, 2002-06-30, 2002-12-31]
flux = [0.10, 0.20, 0.15]             # Bottom flux [cm/d, positive=upward]

[soil.tables.initial_pressure]
depth = [-0.5, -10.0, -50.0, -195.0]  # Depth [cm, negative=below surface]
pressure_head = [-93.0, -80.0, -50.0, 120.0]  # Pressure head [cm]

[soil.tables.hydraulic_parameters]
ores = [0.01, 0.02, 0.015]            # Residual water content per layer
osat = [0.42, 0.38, 0.40]             # Saturated water content per layer
alfa = [0.0276, 0.0213, 0.0250]       # Van Genuchten alpha [1/cm]
npar = [1.491, 1.951, 1.750]          # Van Genuchten n [-]
ksatfit = [12.52, 12.68, 10.30]       # Saturated hydraulic conductivity [cm/d]
bdens = [1315.0, 1315.0, 1400.0]      # Bulk density [mg/cm³]
```

### Fortran Reading (Column-Wise Arrays)

```fortran
real(8), allocatable :: day(:), rate(:)
real(8), allocatable :: ores(:), osat(:), alfa(:), npar(:)
type(toml_table), pointer :: tables_tab, intensity_tab, params_tab
integer :: istat

! Navigate to tables subsection
call get_value(atm_tab, 'tables', tables_tab, stat=istat)

! Read rainfall intensity table
call get_value(tables_tab, 'rainfall_intensity', intensity_tab, stat=istat)
if (istat == 0) then
    call get_value(intensity_tab, 'day', day, stat=istat)
    call get_value(intensity_tab, 'rate', rate, stat=istat)
    
    ! Validate: columns must have same length
    if (size(day) /= size(rate)) then
        call log_error("rainfall_intensity: 'day' and 'rate' must have same length")
    end if
    
    ! Store in state
    state%atm%rain_intensity_day = day
    state%atm%rain_intensity_rate = rate
end if

! Read soil hydraulic parameters table
call get_value(soil_tab, 'tables', tables_tab, stat=istat)
call get_value(tables_tab, 'hydraulic_parameters', params_tab, stat=istat)
if (istat == 0) then
    call get_value(params_tab, 'ores', ores)
    call get_value(params_tab, 'osat', osat)
    call get_value(params_tab, 'alfa', alfa)
    call get_value(params_tab, 'npar', npar)
    
    ! Validate: all columns same length (one entry per layer)
    if (.not. (size(ores) == size(osat) .and. &
            size(osat) == size(alfa) .and. &
            size(alfa) == size(npar))) then
        call log_error("hydraulic_parameters: all columns must have same length")
    end if
end if
```

### Helper Function for 2-Column Tables

```fortran
!> Read a 2-column table (e.g., depth vs. value profile)
subroutine read_table_2col(parent_tab, table_name, col1_name, col2_name, &
                            col1_array, col2_array, required)
    type(toml_table), intent(in) :: parent_tab
    character(len=*), intent(in) :: table_name, col1_name, col2_name
    real(8), allocatable, intent(out) :: col1_array(:), col2_array(:)
    logical, intent(in), optional :: required
    
    type(toml_table), pointer :: tables_tab, table_tab
    integer :: istat
    logical :: req
    
    req = .false.
    if (present(required)) req = required
    
    ! Navigate to tables subsection
    call get_value(parent_tab, 'tables', tables_tab, stat=istat)
    if (istat /= 0) then
        if (req) call log_error("Missing [*.tables] section")
        return
    end if
    
    ! Get specific table
    call get_value(tables_tab, table_name, table_tab, stat=istat)
    if (istat /= 0) then
        if (req) call log_error("Missing table: " // trim(table_name))
        return
    end if
    
    ! Read columns
    call get_value(table_tab, col1_name, col1_array, stat=istat)
    if (istat /= 0) call log_error("Missing column: " // trim(col1_name))
    
    call get_value(table_tab, col2_name, col2_array, stat=istat)
    if (istat /= 0) call log_error("Missing column: " // trim(col2_name))
    
    ! Validate same length
    if (size(col1_array) /= size(col2_array)) then
        call log_error("Table columns must have same length: " // trim(table_name))
    end if
end subroutine

! Usage:
call read_table_2col(soil_tab, 'initial_pressure', 'depth', 'pressure_head', &
                    depth_array, phead_array, required=.true.)
```

### When NOT to Use Tables

**Don't use tables for unrelated arrays**:

```toml
# BAD: These are separate parameters, not a table
[soil.tables.misc]
dz = [10.0, 20.0, 30.0]         # Compartment thickness
bdens = [1315.0, 1315.0]        # Bulk density (different length!)
```

**Use separate keys instead**:

```toml
[soil.discretization]
dz = [10.0, 20.0, 30.0]         # Compartment thickness [cm]

[soil.tables.hydraulic_parameters]
bdens = [1315.0, 1315.0]        # Bulk density per layer [mg/cm³]
```

---

## Arrays

Arrays represent lists of values. TOML supports inline arrays and array-of-tables.

### Inline Arrays (Homogeneous Data)

Simple lists of values:

```toml
[output]
# Array of strings
variables = ["rain", "irrig", "drainage", "gwl", "runoff"]

[soil.discretization]
# Compartment thicknesses per sublayer
hsublay = [10.0, 20.0, 30.0, 140.0]    # Heights [cm]
hcomp = [1.0, 5.0, 5.0, 10.0]          # Compartment size [cm]
ncomp =                  # Number of compartments per sublayer[1][2][3][4]

[output.yearly_balance]
datefix =      # [day, month] for fixed yearly output date[5][6]
```

**Fortran reading**:
```fortran
real(8), allocatable :: hsublay(:), hcomp(:)
integer, allocatable :: ncomp(:)
character(len=:), allocatable :: variables(:)

call get_value(tab, 'hsublay', hsublay, stat=istat)
call get_value(tab, 'ncomp', ncomp, stat=istat)
call get_value(tab, 'variables', variables, stat=istat)
```

### Array of Tables (Heterogeneous Records)

Use double square brackets `[[array_name]]` for multiple records of the same type:

```toml
# Crop rotation: multiple crop periods
[[crop.rotation]]
start_date = 2002-05-01
end_date = 2002-10-15
file = "maize"
type = 1                # 1=simple, 2=WOFOST general, 3=WOFOST grass

[[crop.rotation]]
start_date = 2003-05-10
end_date = 2003-09-29
file = "potato"
type = 2

[[crop.rotation]]
start_date = 2004-01-01
end_date = 2004-12-31
file = "grass"
type = 3

# Irrigation schedule: multiple events
[[irrigation.schedule]]
date = 2002-06-01
depth = 10.0            # Irrigation depth [mm]
concentration = 0.0     # Solute concentration [mg/cm³]
type = 1                # 0=sprinkling, 1=surface

[[irrigation.schedule]]
date = 2002-07-15
depth = 15.0
concentration = 0.0
type = 0
```

**Fortran reading**:
```fortran
type(toml_array), pointer :: rotation_array, irrigation_array
type(toml_table), pointer :: item_tab
type(toml_datetime) :: date_val
character(len=:), allocatable :: file_str
integer :: i, crop_type, irrig_type
real(8) :: depth, concentration

! Read crop rotation array
call get_value(crop_tab, 'rotation', rotation_array, stat=istat)
if (istat == 0) then
    do i = 1, size(rotation_array)
        call get_value(rotation_array, i, item_tab, stat=istat)
        
        call get_value(item_tab, 'start_date', date_val)
        call get_value(item_tab, 'end_date', date_val)
        call get_value(item_tab, 'file', file_str)
        call get_value(item_tab, 'type', crop_type)
        
        ! Store in crop rotation schedule
        state%crop%rotation(i)%start_date = convert_to_swap_time(date_val)
        state%crop%rotation(i)%file = trim(file_str)
        state%crop%rotation(i)%type = crop_type
    end do
end if

! Read irrigation schedule array
call get_value(irrig_tab, 'schedule', irrigation_array, stat=istat)
if (istat == 0) then
    do i = 1, size(irrigation_array)
        call get_value(irrigation_array, i, item_tab)
        
        call get_value(item_tab, 'date', date_val)
        call get_value(item_tab, 'depth', depth)
        call get_value(item_tab, 'concentration', concentration)
        call get_value(item_tab, 'type', irrig_type)
        
        ! Store in irrigation event array
        state%irrig%events(i)%date = convert_to_swap_time(date_val)
        state%irrig%events(i)%depth = depth
        state%irrig%events(i)%concentration = concentration
        state%irrig%events(i)%type = irrig_type
    end do
end if
```

### When to Use Array-of-Tables vs. Table-with-Columns

**Use array-of-tables `[[...]]` when:**
- Each entry represents a distinct **event** or **record** (irrigation, crop periods)
- Entries may have different structures or optional fields
- Natural to process sequentially (loop over events)

**Use table-with-columns `[....tables.xxx]` when:**
- Data forms a **lookup table** or **profile** (depth vs. property)
- All rows have the same structure (same columns)
- Need to access columns independently (e.g., interpolation)

**Example comparison**:

```toml
# Array-of-tables: Irrigation events (distinct records)
[[irrigation.schedule]]
date = 2002-06-01
depth = 10.0
type = 1

[[irrigation.schedule]]
date = 2002-07-15
depth = 15.0
type = 0

# Table-with-columns: Bottom flux time series (lookup table)
[boundary.tables.prescribed_flux]
date = [2002-01-01, 2002-06-30, 2002-12-31]
flux = [0.10, 0.20, 0.15]
```

---

## Validation Best Practices

### Check Section Existence

```fortran
call get_value(doc, 'atmosphere', atm_tab, requested=.true., stat=istat)
if (istat /= 0) then
    call log_error("Missing required section: [atmosphere]")
    stop
end if
```

### Use Defaults for Optional Parameters

```fortran
! Optional parameter with default
call get_value(tab, 'swmetdetail', state%atm%swmetdetail, default=0)

! Required parameter (no default)
call get_value(tab, 'lat', state%atm%lat, stat=istat)
if (istat /= 0) call log_error("Missing required parameter: lat")
```

### Validate Switch-Dependent Sections

```fortran
! Read switch first
call get_value(tab, 'swrain', swrain, default=0)

! Then validate required subsection based on switch
if (swrain == 1) then
    ! swrain=1 requires intensity table
    call get_value(tab, 'tables', tables_tab, stat=istat)
    if (istat /= 0) then
        call log_error("swrain=1 requires [atmosphere.tables] section")
    end if
    
    call get_value(tables_tab, 'rainfall_intensity', intensity_tab, stat=istat)
    if (istat /= 0) then
        call log_error("swrain=1 requires [atmosphere.tables.rainfall_intensity]")
    end if
end if
```

### Validate Array Lengths

```fortran
! For tables: all columns must have same length
if (size(day_array) /= size(rate_array)) then
    call log_error("rainfall_intensity: columns 'day' and 'rate' must have same length")
end if

! For layer parameters: must match number of layers
if (size(ores) /= state%soil%numlay) then
    write(msg, '(A,I0,A,I0)') "Expected ", state%soil%numlay, &
        " layer parameters, got ", size(ores)
    call log_error(trim(msg))
end if
```

### Validate Value Ranges

```fortran
call get_value(tab, 'lat', lat)
if (lat < -90.0d0 .or. lat > 90.0d0) then
    write(msg, '(A,F0.2,A)') "Invalid latitude: ", lat, " (must be -90..90)"
    call log_error(trim(msg))
end if
```

---

## Complete Example

### TOML Configuration File

```toml
# =============================================================================
# SWAP Model Configuration
# =============================================================================

[general]
project = "hupsel_test"

[simulation]
start_date = 2002-01-01T00:00:00
end_date = 2004-12-31T23:59:59

[atmosphere]
file = "meteo.met"
lat = 52.0         # Latitude [degrees, -90..90]
alt = 10.0         # Altitude [m]
angstrom_a = 0.25  # Angstrom coefficient [-]
angstrom_b = 0.5

[atmosphere.rainfall]
swrain = 1         # 0=daily, 1=+intensity table, 2=+duration, 3=file

[atmosphere.tables.rainfall_intensity]
day = [1.0, 90.0, 180.0, 360.0]
rate = [15.0, 25.0, 20.0, 15.0]  # mm/d

[soil.initial]
swinco = 2         # 1=pressure profile, 2=hydrostatic, 3=file
gwli = -75.0       # Initial GWL [cm]

[soil.tables.hydraulic_parameters]
ores = [0.01, 0.02]
osat = [0.42, 0.38]
alfa = [0.0276, 0.0213]
npar = [1.491, 1.951]

[[crop.rotation]]
start_date = 2002-05-01
end_date = 2002-10-15
file = "maize"
type = 1

[[crop.rotation]]
start_date = 2003-05-10
end_date = 2003-09-29
file = "potato"
type = 2
```

### Fortran Reading Code

```fortran
subroutine read_config(filename, state)
    character(len=*), intent(in) :: filename
    type(swap_state_t), intent(inout) :: state
    
    type(toml_table), pointer :: doc, atm_tab, tables_tab, rain_tab
    type(toml_array), pointer :: rotation_array
    integer :: istat
    
    ! Parse TOML file
    call toml_parse(doc, filename, error, stat=istat)
    if (istat /= 0) then
        call log_error("Failed to parse " // trim(filename))
        stop
    end if
    
    ! Read general section
    call get_value(doc, 'general', gen_tab, stat=istat)
    call get_value(gen_tab, 'project', s)
    state%project = trim(s)
    
    ! Read atmosphere section
    call get_value(doc, 'atmosphere', atm_tab, requested=.true., stat=istat)
    if (istat /= 0) call log_error("Missing [atmosphere] section")
    
    call get_value(atm_tab, 'file', s)
    state%atm%metfil = trim(s)
    
    call get_value(atm_tab, 'lat', state%atm%lat)
    if (state%atm%lat < -90.0d0 .or. state%atm%lat > 90.0d0) then
        call log_error("Invalid latitude")
    end if
    
    ! Read rainfall configuration
    call get_value(atm_tab, 'rainfall', rain_tab, stat=istat)
    call get_value(rain_tab, 'swrain', state%atm%swrain, default=0)
    
    ! Conditional: read intensity table if swrain=1
    if (state%atm%swrain == 1) then
        call get_value(atm_tab, 'tables', tables_tab, stat=istat)
        call read_table_2col(atm_tab, 'rainfall_intensity', 'day', 'rate', &
                            day_array, rate_array, required=.true.)
        state%atm%rain_day = day_array
        state%atm%rain_rate = rate_array
    end if
    
    ! Read crop rotation (array of tables)
    call get_value(doc, 'crop', crop_tab, stat=istat)
    call get_value(crop_tab, 'rotation', rotation_array, stat=istat)
    if (istat == 0) then
        allocate(state%crop%rotation(size(rotation_array)))
        do i = 1, size(rotation_array)
            call get_value(rotation_array, i, item_tab)
            call get_value(item_tab, 'start_date', date_val)
            call get_value(item_tab, 'file', file_str)
            state%crop%rotation(i)%file = trim(file_str)
        end do
    end if
end subroutine
```

---

## Migration from Legacy .swp Format

### Key Differences

| Feature | Legacy .swp | Modern TOML |
|---------|-------------|-------------|
| **Syntax** | Free-form, whitespace-delimited | Structured key-value pairs |
| **Types** | All strings (manual conversion) | Native types (integer, float, string, bool, date) |
| **Comments** | `*` or `!` | `#` |
| **Sections** | Comment-based, no enforcement | Explicit `[section]` headers |
| **Tables** | Inline or separate files | Inline arrays with column names |
| **Validation** | Manual in Fortran code | Schema-based (external tools) |

### Conversion Guidelines

1. **Group related parameters** under logical sections (`[atmosphere]`, `[soil]`, etc.)
2. **Use descriptive keys** instead of positional arguments
3. **Leverage native types** (dates as `YYYY-MM-DD`, not day-since-1900)
4. **Inline tables** where appropriate (avoid separate `.MET`, `.IRG` files for simple data)
5. **Document units** in comments for every parameter
6. **Use array-of-tables** for repeating structures (crop rotation, irrigation events)

---

## Additional Resources

- **TOML Specification**: https://toml.io/
- **toml-f Library**: https://github.com/toml-f/toml-f
- **SWAP TOML Schema** (for validation): `docs/schema/swp_schema.toml`
- **Example Configurations**: `examples/toml/`
