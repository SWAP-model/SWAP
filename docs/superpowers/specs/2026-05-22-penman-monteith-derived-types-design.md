# Penman–Monteith derived-type refactor — design

**Date:** 2026-05-22
**Scope:** [src/atmosphere/et.f90](../../../src/atmosphere/et.f90) — Penman–Monteith subroutines only
**Goal:** Replace the ~40-positional-argument signatures of `PenMon_calc` and `PenMon` with two derived types — `pm_inputs_t` (intent in) and `pm_outputs_t` (intent out) — so the single call site and the characterization test become readable and adding/removing a parameter no longer touches every caller.

## Background

`PenMon_calc` is the pure Penman–Monteith kernel in `et_mod`; `PenMon` is its non-pure wrapper that adds `astro()` and warning logging. Both expose ~40 positional arguments, which produces a 9-line call site in [meteoday.f90:685](../../../src/atmosphere/meteoday.f90#L685) and a fragile positional test fixture in [tests/unit/atmosphere/test_et.pf:36](../../../tests/unit/atmosphere/test_et.pf#L36).

This refactor is the first execution of theme #6 ("Collapse the PenMon argument soup into derived types") from the [architecture review](../../../src/atmosphere/README.md#cross-cutting-modernisation-themes) drafted on 2026-05-22.

## Inventory of impacted code

| Site | Role | File |
|------|------|------|
| `PenMon_calc` | pure kernel, public | `src/atmosphere/et.f90` |
| `PenMon` | I/O wrapper, public | `src/atmosphere/et.f90` |
| `ProcessMeteoDay` | sole production caller (calls `PenMon`) | `src/atmosphere/meteoday.f90` |
| `test_et.pf` | characterization test (calls `PenMon_calc` directly) | `tests/unit/atmosphere/test_et.pf` |
| `et_test_stubs.f90` | link-only stubs for the wrapper's I/O dependencies | `tests/unit/atmosphere/et_test_stubs.f90` |

No other callers, no other tests. The stub file is unaffected by the signature change.

## Design

### New types (in `et_mod`, top of file, both `public`)

```fortran
type, public :: pm_inputs_t
    ! Time / control
    integer :: daynr           = 0       ! day of year
    integer :: irecord         = 1       ! current sub-daily record
    integer :: nmetdetail      = 1       ! records/day (sub-daily)
    logical :: flmetdetail     = .false. ! sub-daily mode
    logical :: flCropEmergence = .false. ! crop present
    integer :: swcf            = 1       ! crop-factor switch (1 cf, 2 ch, 3 cfeic)
    integer :: swdivide        = 0       ! 0 = standard PM, 1 = PMdirect

    ! Site
    real(real64) :: lat  = 0.0_real64
    real(real64) :: alt  = 0.0_real64
    real(real64) :: altw = 0.0_real64
    real(real64) :: a    = 0.0_real64    ! Angstrom a
    real(real64) :: b    = 0.0_real64    ! Angstrom b
    real(real64) :: rcs  = 0.0_real64    ! soil reflection coefficient

    ! Atmospheric forcing
    real(real64) :: rad   = 0.0_real64
    real(real64) :: tav   = 0.0_real64
    real(real64) :: tmn   = 0.0_real64
    real(real64) :: tmx   = 0.0_real64
    real(real64) :: hum   = 0.0_real64
    real(real64) :: win   = 0.0_real64
    real(real64) :: atmtr = 0.0_real64
    real(real64) :: difpp = 0.0_real64
    real(real64) :: dsinbe = 0.0_real64
    real(real64) :: daylp = 0.0_real64

    ! Crop / canopy
    real(real64) :: rsc    = 0.0_real64
    real(real64) :: rsw    = 0.0_real64
    real(real64) :: ch     = 0.0_real64
    real(real64) :: albedo = 0.0_real64
    real(real64) :: kdif   = 0.0_real64
    real(real64) :: kdir   = 0.0_real64
    real(real64) :: lai    = 0.0_real64

    ! PMdirect-only
    real(real64) :: rsoil  = 0.0_real64
end type pm_inputs_t

type, public :: pm_outputs_t
    real(real64) :: es0         = 0.0_real64
    real(real64) :: et0         = 0.0_real64
    real(real64) :: ew0         = 0.0_real64
    real(real64) :: Edirect     = 0.0_real64
    real(real64) :: Tdirect     = 0.0_real64
    real(real64) :: Tdirectwet  = 0.0_real64
    real(real64) :: Edirectpond = 0.0_real64
    integer      :: warning_code = 0  ! 0=ok, 1=polar 0hrs, 2=polar 24hrs
end type pm_outputs_t
```

Field names preserve the existing SWAP abbreviations (`rad`, `tav`, `rsc`, `kdif`, …). This minimises diff inside the kernel body — only the variable references change, not the computation.

Default initialisers (`= 0.0_real64`) guard against accidentally reading an unset field at the call site.

### New procedure signatures

```fortran
pure subroutine PenMon_calc(inputs, outputs)
    type(pm_inputs_t),  intent(in)  :: inputs
    type(pm_outputs_t), intent(out) :: outputs
    ! body unchanged except identifiers like `lat`  ->  `inputs%lat`
    !                                        and `es0` -> `outputs%es0`.

subroutine PenMon(inputs, outputs, logf, swscre)
    type(pm_inputs_t),  intent(in)    :: inputs
    type(pm_outputs_t), intent(inout) :: outputs   ! inout so warning_code survives
    integer,            intent(in)    :: logf, swscre
    ! still calls astro() for daily branch, then PenMon_calc, then forwards
    ! outputs%warning_code to warn().
```

### Call-site change ([meteoday.f90:685](../../../src/atmosphere/meteoday.f90#L685))

Replace:

```fortran
call PenMon (logf,swscre,tc_daynr,state%cfg%meteo%lat,state%cfg%meteo%alt, &
             state%cfg%meteo%altw,angstroma,angstromb,rcs,rad,state%atmosphere%Tav, &
             hum,win,state%crop%common%rsc,state%crop%es0,state%crop%et0, &
             state%crop%ew0,state%crop%swcf,state%crop%common%ch, &
             state%crop%flCropEmergence,daylp,tc_flmetdetail,irecord, &
             config%meteo%nmetdetail,state%crop%common%albedo,tmn,tmx, &
             state%crop%common%rsw,difpp,dsinbe,atmtr,Edirect,Tdirect, &
             Tdirectwet,rsoil,state%cfg%meteo%swdivide,state%crop%kdif, &
             state%crop%kdir,state%crop%lai,Edirectpond)
```

With (grouped by source for readability):

```fortran
! Pack PM inputs.
pmi%daynr       = tc_daynr
pmi%irecord     = irecord
pmi%nmetdetail  = config%meteo%nmetdetail
pmi%flmetdetail = tc_flmetdetail
pmi%flCropEmergence = state%crop%flCropEmergence
pmi%swcf        = state%crop%swcf
pmi%swdivide    = state%cfg%meteo%swdivide

pmi%lat  = state%cfg%meteo%lat
pmi%alt  = state%cfg%meteo%alt
pmi%altw = state%cfg%meteo%altw
pmi%a    = angstroma
pmi%b    = angstromb
pmi%rcs  = rcs

pmi%rad    = rad
pmi%tav    = state%atmosphere%Tav
pmi%tmn    = tmn
pmi%tmx    = tmx
pmi%hum    = hum
pmi%win    = win
pmi%atmtr  = atmtr
pmi%difpp  = difpp
pmi%dsinbe = dsinbe
pmi%daylp  = daylp

pmi%rsc    = state%crop%common%rsc
pmi%rsw    = state%crop%common%rsw
pmi%ch     = state%crop%common%ch
pmi%albedo = state%crop%common%albedo
pmi%kdif   = state%crop%kdif
pmi%kdir   = state%crop%kdir
pmi%lai    = state%crop%lai

pmi%rsoil  = rsoil

call PenMon(pmi, pmo, logf, swscre)

! Unpack PM outputs to existing state / locals.
state%crop%es0 = pmo%es0
state%crop%et0 = pmo%et0
state%crop%ew0 = pmo%ew0
Edirect     = pmo%Edirect
Tdirect     = pmo%Tdirect
Tdirectwet  = pmo%Tdirectwet
Edirectpond = pmo%Edirectpond
```

`pmi` and `pmo` declared as locals at the top of `ProcessMeteoDay` (alongside `rcs`), before the `associate` block.

**Position in control flow.** Pack/call/unpack sits **inside the existing `do 1000 irecord = 1, ndayparts` loop, inside the `elseif (config%meteo%swmetdetail.eq.1 .or. config%meteo%swetr.eq.0)` branch, immediately before the call** — replacing the current `call PenMon(...)` block in place. The other branch (daily + `swetr=1`) bypasses PenMon entirely and is unchanged.

**Why repopulate every iteration:** `ndayparts ≤ ~96`, ~30 scalar assignments per call is invisible overhead next to the kernel's transcendentals, and keeping the full pack at one site means future field additions cannot be silently forgotten. See [Call frequency](#call-frequency) below for numbers.

**Why `pm_outputs_t` is safe to reuse across iterations:** Fortran's `intent(out)` derived-type semantics reset every component to its default-init value on entry to `PenMon_calc`. The `= 0.0_real64` / `= 0` defaults in `pm_outputs_t` are what make this guarantee load-bearing.

### Call frequency

| Mode | `ndayparts` | Calls / day | Typical 10-yr run |
|------|-------------|-------------|-------------------|
| Daily, `swetr=0` (PM computes ET) | 1 | 1 | ~3,650 |
| Daily, `swetr=1` (PM bypassed) | 1 | 0 | 0 |
| Sub-daily, hourly (`nmetdetail=24`) | 24 | 24 | ~87,600 |
| Sub-daily, 15-min (`nmetdetail=96`) | 96 | 96 | ~350,000 |

Per-call cost split: kernel body has dozens of `exp`/`log`/`sin`/`cos`/`asin`/`sqrt` (microseconds); pack of 30 scalar fields is 30 memcpys (nanoseconds). Three+ orders of magnitude apart — pack overhead is invisible even at the high-end case.

### Test change ([test_et.pf:36](../../../tests/unit/atmosphere/test_et.pf#L36))

Same shape as the production call-site change — field assignments to a `pm_inputs_t` local, single `call PenMon_calc(pmi, pmo)`, then `@assertEqual` against `pmo%es0`, `pmo%et0`, etc.

The 27 expected values in the characterization fixture are **byte-identical** to today (no physics change), so the assertion values are copied verbatim.

`et_test_stubs.f90` does not need changes — it only stubs procedures the wrapper depends on (`astro`, `warn`), not anything in the signature.

## Out of scope

Deliberately not in this refactor:

- Physics changes of any kind.
- `real(8)` → `real(real64)` migration of the rest of `et_mod`.
- `reduceva` / `black_reduction` / `boesten_stroosnijder_reduction` refactor.
- Hoisting constants to `atmosphere_constants`.
- Touching the legacy `variables` imports inside `et_mod`.
- Migrating `MeteoVars` or `Vars`.
- Splitting `meteoday.f90`'s four modules.

These are catalogued in the [atmosphere README](../../../src/atmosphere/README.md#cross-cutting-modernisation-themes) for later passes.

## Why not precompute PenMon up front for the whole simulation?

A question that comes up reading the code: meteo is read year-by-year — why not also precompute the daily ET₀ / EW₀ / ES₀ at meteo-load time and cache them?

The answer is that **PenMon depends on internal model state, not just on meteo**. Categorising the inputs:

| Category | Inputs | Source |
|----------|--------|--------|
| Site-static | `lat`, `alt`, `altw`, `a`, `b`, `rcs`, `nmetdetail`, `swdivide` | config — set once |
| Daily meteo | `rad`, `tmn`, `tmx`, `hum`, `win`, `Tav` | meteo file |
| Astronomical | `atmtr`, `difpp`, `dsinbe`, `daylp` | pure function of `daynr`, `lat`, `rad` |
| **Crop state — produced by the simulation** | `lai`, `ch`, `flCropEmergence`, `albedo`, `kdif`, `kdir`, `rsc`, `rsw` | crop growth model output, changes daily |
| **Soil state — produced by the simulation** | `rsoil` (PMdirect only) | soil resistance model output |

The bottom two rows block whole-simulation precomputation. LAI and crop height are not inputs — they are computed by the crop growth model during the same timestep, in an explicit forward-coupled cycle:

```
crop growth (day N)  ->  today's LAI, height, albedo
                          |
                          v
                       PenMon  ->  today's ET₀ / EW₀ / ES₀
                          |
                          v
                  reduceva + water uptake  ->  actual T, E (today)
                          |
                          v
              water stress = T_actual / T_potential
                          |
                          v
              crop growth (day N+1)  ->  tomorrow's LAI ...
```

The `swetr=1` branch in `ProcessMeteoDay` *is* the "precompute externally" escape hatch — when the user supplies reference ET in the meteo file (e.g. FAO-56 ET₀ from a separate tool), SWAP skips PenMon entirely and uses crop factors. PenMon is only invoked when the simulation needs to respond dynamically to the modelled crop's LAI/height/resistances.

## Forward-compatibility with partial precomputation

The categorisation above also suggests a real efficiency optimisation **for a future arc** (out of scope here): the air-only thermodynamics, astronomy, soil aerodynamics, and net longwave terms depend only on (meteo + site) and could be precomputed once per day into a `pm_atmos_baseline_t` cache. The kernel would then only compute the crop-coupled bits.

If/when that arc happens, the natural shape is:

```fortran
type :: pm_atmos_baseline_t   ! computed once per day from meteo + site
   real(real64) :: lambda, delta, ea, ed, vpd, gamma, palt, rho, cp
   real(real64) :: rns, rnp, rnl, sinld, cosld, dayl, sunrise, sunset
   real(real64) :: atmtr, difpp, dsinbe, daylp, ras, gs
end type

type :: pm_crop_state_t       ! produced per day by the crop growth model
   real(real64) :: lai, ch, albedo, kdif, kdir, rsc, rsw
   logical      :: flCropEmergence
   integer      :: swcf
end type

pure subroutine PenMon_calc(baseline, crop, outputs)
```

The field groupings in today's flat `pm_inputs_t` (`! Site`, `! Atmospheric forcing`, `! Crop / canopy`) map cleanly to that split — when the time comes, the migration is field-relocation, not redesign. The deciding signal for the future arc is a profiler showing PenMon is a hotspot, not a design preference.

Splitting today would add boilerplate without payoff: it makes the call site more complex, not less, and locks in answers (cache location, baseline ownership, recomputation triggers) we have not yet had to ask. So we keep one flat `pm_inputs_t` and document the migration path.

## Implementation order (TDD)

1. **Update `test_et.pf`** to expect the new signature: build `pm_inputs_t`, call `PenMon_calc(pmi, pmo)`, assert on `pmo%X`. Compile — fails (types don't exist, signature wrong).
2. **Add `pm_inputs_t` + `pm_outputs_t`** at the top of `et_mod`. Compile — test still fails (kernel signature still positional).
3. **Rewrite `PenMon_calc`** to take the new signature; rewrite the body's variable references (`tav` → `inputs%tav`, `es0` → `outputs%es0`, etc.). Compile — test should pass with byte-identical values. **Full build still fails** because `PenMon` (the wrapper) still calls the old `PenMon_calc` positionally and `meteoday.f90` still calls the old `PenMon` positionally. This is the TDD red-green moment for the pure kernel — verify the test before continuing.
4. **Combined commit:** rewrite `PenMon` wrapper to take the new signature *and* migrate the `ProcessMeteoDay` call site in `meteoday.f90`. These must land together because changing the wrapper signature breaks the only caller. Compile — full build passes.
5. **Run check-fast** (the per-task regression gate): full pFUnit + check-fast must be green.

Steps 1-4 are each one commit (4 total). Step 5 is the verification gate, not a commit.

## Risks and mitigations

| Risk | Mitigation |
|------|------------|
| `intent(out)` on derived type with default initialisers zeros every field on entry — risks silently dropping a forgotten input | Standard Fortran behaviour for `intent(out)` derived types; this is the *intent*, since `PenMon_calc` is supposed to populate all of `outputs`. Test will catch any output regression. |
| 27 byte-identical assertions in `test_et.pf` rely on the kernel not changing — but step 3 touches every line of the body | Only identifier substitution (`tav` → `inputs%tav`), not arithmetic. Run the test after step 3 before doing anything else. If even one value drifts, revert and investigate. |
| Default initialisers (`= 0.0_real64`) require Fortran 2003; need to confirm the SWAP toolchain supports it | SWAP already uses derived types with default initialisers in `swap_state_mod` and `swap_config_mod`. Confirmed safe. |
| Two new public exports (`pm_inputs_t`, `pm_outputs_t`) widen `et_mod`'s API surface | Acceptable — these are the new public contract. No downstream module imports them yet; only `meteoday.f90` will. |

## Verification

- **Per-task regression gate:** every commit ends with `check-fast` (per the project's `feedback_per_task_regression_gate` playbook entry — pFUnit alone misses global-default regressions).
- **Specifically:** `test_et.pf` must remain green after step 3; full check-fast must be green after step 4.
- **No clean rebuild needed:** this refactor does not touch the `state` schema, so the `feedback_state_schema_clean_rebuild` constraint does not apply. Incremental Meson build is fine.

## Definition of done

- `PenMon_calc` and `PenMon` take `(pm_inputs_t, pm_outputs_t [, logf, swscre])`.
- `pm_inputs_t` and `pm_outputs_t` are public, exported from `et_mod`.
- `meteoday.f90` call site uses field assignments + one `call PenMon(pmi, pmo, logf, swscre)`.
- `test_et.pf` uses the new signature with byte-identical expected values.
- Full `check-fast` is green after step 4.
- No legacy positional-argument `PenMon` / `PenMon_calc` signature remains in the codebase.
