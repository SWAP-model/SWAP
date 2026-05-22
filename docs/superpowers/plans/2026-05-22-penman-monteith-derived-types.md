# Penman–Monteith Derived-Type Refactor — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the ~40-positional-argument signatures of `PenMon_calc` and `PenMon` in `src/atmosphere/et.f90` with two derived types — `pm_inputs_t` (intent in) and `pm_outputs_t` (intent out) — so the single call site and the characterization test become readable and signature changes no longer touch every caller.

**Architecture:** Two new public types declared at the top of `et_mod`. `PenMon_calc` stays `pure`; its body's identifier references are prefixed with `inputs%` / `outputs%`. `PenMon` wrapper keeps `(logf, swscre)` and threads the two types through. The sole production caller in `ProcessMeteoDay` (`meteoday.f90`) packs the inputs from existing state/config/locals just before the call and unpacks the outputs immediately after. The characterization test in `test_et.pf` is rewritten in the same pack-call-unpack shape, with byte-identical expected values.

**Tech Stack:** Fortran 2003+ (derived types with default initialisers); Meson; pFUnit; pixi task runner.

**Spec:** [docs/superpowers/specs/2026-05-22-penman-monteith-derived-types-design.md](../specs/2026-05-22-penman-monteith-derived-types-design.md)

---

## File structure

- **Modify:** `src/atmosphere/et.f90` — add two derived types at the top of `et_mod`; rewrite signatures and bodies of `PenMon_calc` and `PenMon`. ~410 lines of `et.f90` change; the `reduceva` / `black_reduction` / `boesten_stroosnijder_reduction` block is untouched.
- **Modify:** `src/atmosphere/meteoday.f90` — replace the 9-line positional `call PenMon(...)` in `ProcessMeteoDay` (around line 685) with a pack-call-unpack block. Add two local declarations (`pmi`, `pmo`) at the top of `ProcessMeteoDay`.
- **Modify:** `tests/unit/atmosphere/test_et.pf` — rewrite the single test to use the new signature. Expected values stay byte-identical.

No new files. No file deletions. `tests/unit/atmosphere/et_test_stubs.f90` is untouched (stubs `astro` and `warn`, neither of which is in the new signature).

## Commit plan

Three commits, executed in this order. The build is **broken between commit 2 and commit 3** because the `unit-swap-tests` executable links production sources — there is no isolated build-just-the-test target. Commit 3 must be created in the same session as commit 2; do not push commit 2 alone.

| Commit | What it does | Build state after | Test state after |
|--------|--------------|-------------------|------------------|
| 1 | Update `test_et.pf` to use new types/signature | broken (test only) | won't compile |
| 2 | Add `pm_inputs_t` + `pm_outputs_t` to `et_mod` | still broken | still won't compile |
| 3 | Rewrite `PenMon_calc`, `PenMon`, and `meteoday.f90` caller | green | green |

Plus a final verification gate (not a commit): `check-fast`.

(The spec proposed 4 commits with a clean intermediate state between kernel and wrapper changes; that's unworkable because the wrapper's body calls the kernel, so they must change together. Consolidating to 3 commits — the spec's intent of "test-first" is preserved.)

---

## Task 1: Update `test_et.pf` to use the new signature

**Files:**
- Modify: `tests/unit/atmosphere/test_et.pf`

Test-first commit. After this commit, the test won't compile — types and signature don't exist yet. The next two commits fix it.

- [ ] **Step 1.1: Read the current test to confirm expected values**

Read `tests/unit/atmosphere/test_et.pf` and note the three `@assertEqual` values. They must reappear byte-identical:
- `et0 == 4.6108539555674373d0` (tol `1.0d-4`)
- `ew0 == 6.2483978923071479d0`
- `es0 == 5.0880237214078035d0`

- [ ] **Step 1.2: Rewrite `tests/unit/atmosphere/test_et.pf`**

Replace the entire file contents with:

```fortran
! Characterization test for PenMon_calc (Penman-Monteith reference ET).
!
! Locks in the current gfortran build's output for a fixed set of inputs.
! If the physics intentionally changes, update the expected values in a
! separate commit from the physics change.
!
! PenMon_calc is a pure subroutine; it requires no I/O and has no side
! effects. Inputs represent a mid-summer, mid-latitude day (day 180,
! latitude 52 N) with standard meteorological conditions.
!
! Expected values are the observed outputs of the current gfortran build
! and are characterization-only (not hand-validated against FAO-56).

@test
subroutine test_penmon_midsummer_midlatitude_characterization()
   use funit
   use et_mod, only: PenMon_calc, pm_inputs_t, pm_outputs_t
   implicit none

   type(pm_inputs_t)  :: pmi
   type(pm_outputs_t) :: pmo

   ! Inputs: late-June midlatitude day (daily mode, flmetdetail=.false.)
   pmi%daynr           = 180
   pmi%irecord         = 1
   pmi%nmetdetail      = 1
   pmi%flmetdetail     = .false.
   pmi%flCropEmergence = .true.
   pmi%swcf            = 1
   pmi%swdivide        = 0

   pmi%lat   = 52.0d0
   pmi%alt   = 10.0d0
   pmi%altw  = 2.0d0
   pmi%a     = 0.25d0
   pmi%b     = 0.50d0
   pmi%rcs   = 0.15d0

   pmi%rad    = 25.0d6
   pmi%tav    = 16.0d0
   pmi%tmn    = 10.0d0
   pmi%tmx    = 22.0d0
   pmi%hum    = 1.2d0
   pmi%win    = 3.0d0
   pmi%atmtr  = 0.6d0
   pmi%difpp  = 0.0d0
   pmi%dsinbe = 0.5d0
   pmi%daylp  = 0.0d0

   pmi%rsc    = 70.0d0
   pmi%rsw    = 0.0d0
   pmi%ch     = 12.0d0
   pmi%albedo = 0.23d0
   pmi%kdif   = 0.5d0
   pmi%kdir   = 0.5d0
   pmi%lai    = 3.0d0

   pmi%rsoil  = 0.0d0

   call PenMon_calc(pmi, pmo)

   @assertEqual(4.6108539555674373d0, pmo%et0, 1.0d-4)
   @assertEqual(6.2483978923071479d0, pmo%ew0, 1.0d-4)
   @assertEqual(5.0880237214078035d0, pmo%es0, 1.0d-4)
end subroutine
```

- [ ] **Step 1.3: Confirm the test won't compile yet (expected)**

Run: `pixi run build-linux`
Expected: compile failure — `et_mod` does not export `pm_inputs_t`, `pm_outputs_t`. This is intentional; the next commits fix it.

Do **not** run `pixi run test-pfunit` — it will also fail.

- [ ] **Step 1.4: Commit**

```bash
git add tests/unit/atmosphere/test_et.pf
git commit -m "test(et): switch PenMon_calc test to derived-type signature

Prepares the characterization test for the upcoming derived-type
refactor. Will not compile until pm_inputs_t / pm_outputs_t and the
new PenMon_calc signature land in et_mod.

Spec: docs/superpowers/specs/2026-05-22-penman-monteith-derived-types-design.md"
```

---

## Task 2: Add `pm_inputs_t` and `pm_outputs_t` to `et_mod`

**Files:**
- Modify: `src/atmosphere/et.f90` (declarations only; subroutine bodies untouched)

After this commit the types exist but the kernel signature is still positional, so the test still won't compile. The build is still broken.

- [ ] **Step 2.1: Add `real64` import and type declarations**

In `src/atmosphere/et.f90`, modify the module header. Find:

```fortran
module et_mod
   use error_mod, only: fatalerr_collected
   use swap_state_mod, only: swap_state_t
    implicit none
    private

  public :: PenMon, PenMon_calc, reduceva
contains
```

Replace with:

```fortran
module et_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use error_mod, only: fatalerr_collected
   use swap_state_mod, only: swap_state_t
    implicit none
    private

  public :: PenMon, PenMon_calc, reduceva
  public :: pm_inputs_t, pm_outputs_t

  type :: pm_inputs_t
      ! Time / control
      integer :: daynr           = 0
      integer :: irecord         = 1
      integer :: nmetdetail      = 1
      logical :: flmetdetail     = .false.
      logical :: flCropEmergence = .false.
      integer :: swcf            = 1
      integer :: swdivide        = 0

      ! Site
      real(real64) :: lat  = 0.0_real64
      real(real64) :: alt  = 0.0_real64
      real(real64) :: altw = 0.0_real64
      real(real64) :: a    = 0.0_real64
      real(real64) :: b    = 0.0_real64
      real(real64) :: rcs  = 0.0_real64

      ! Atmospheric forcing
      real(real64) :: rad    = 0.0_real64
      real(real64) :: tav    = 0.0_real64
      real(real64) :: tmn    = 0.0_real64
      real(real64) :: tmx    = 0.0_real64
      real(real64) :: hum    = 0.0_real64
      real(real64) :: win    = 0.0_real64
      real(real64) :: atmtr  = 0.0_real64
      real(real64) :: difpp  = 0.0_real64
      real(real64) :: dsinbe = 0.0_real64
      real(real64) :: daylp  = 0.0_real64

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

  type :: pm_outputs_t
      real(real64) :: es0         = 0.0_real64
      real(real64) :: et0         = 0.0_real64
      real(real64) :: ew0         = 0.0_real64
      real(real64) :: Edirect     = 0.0_real64
      real(real64) :: Tdirect     = 0.0_real64
      real(real64) :: Tdirectwet  = 0.0_real64
      real(real64) :: Edirectpond = 0.0_real64
      integer      :: warning_code = 0   ! 0=ok, 1=polar 0hrs, 2=polar 24hrs
  end type pm_outputs_t

contains
```

- [ ] **Step 2.2: Confirm the types compile in isolation**

Build only the `et_mod` object to verify the types are syntactically valid:

Run: `pixi run build-linux 2>&1 | head -60`
Expected: `et.f90` compiles successfully on its own. Linking still fails (test references the new signature; `meteoday.f90` still calls old signature internally is fine — but the test references `pm_inputs_t`). Note: this step is just a sanity check that the type declarations are syntactically correct; failure to link is expected.

If `et.f90` itself fails to compile, fix the type declarations before continuing.

- [ ] **Step 2.3: Commit**

```bash
git add src/atmosphere/et.f90
git commit -m "feat(et): add pm_inputs_t and pm_outputs_t derived types

Declares the new derived types for the PenMon refactor; the kernel
and wrapper still use the old positional signatures. Test_et.pf and
the full build remain broken until Task 3 lands.

Spec: docs/superpowers/specs/2026-05-22-penman-monteith-derived-types-design.md"
```

---

## Task 3: Rewrite `PenMon_calc`, `PenMon`, and the `meteoday.f90` caller

**Files:**
- Modify: `src/atmosphere/et.f90` (rewrite `PenMon_calc` and `PenMon` signatures + bodies)
- Modify: `src/atmosphere/meteoday.f90` (rewrite the `call PenMon(...)` block in `ProcessMeteoDay`)

All three changes land in one commit because the wrapper's body calls the kernel and the caller calls the wrapper — they must move in lockstep. After this commit the build is green and the test passes.

- [ ] **Step 3.1: Rewrite `PenMon_calc` signature and body**

In `src/atmosphere/et.f90`, find the existing `PenMon_calc` declaration (currently starting `pure subroutine PenMon_calc(daynr, lat, alt, altw, ...`).

Replace the entire `pure subroutine PenMon_calc(...) end subroutine PenMon_calc` block with the version below. The body's computation is identical — only identifier references are prefixed with `inputs%` (for the 31 input fields) or `outputs%` (for the 7 outputs + `warning_code`). The argument-doc comment block is replaced with a one-line reference to the types.

```fortran
  !> Pure Penman-Monteith calculation (no I/O, fully deterministic)
  !!
  !! Performs all ET calculations without side effects. This pure version
  !! enables compiler optimizations, parallelization, and easier testing.
  !!
  !! @note Warnings are returned via outputs%warning_code; see pm_outputs_t.
  pure subroutine PenMon_calc(inputs, outputs)
      use swap_constants, only: vlarge, small, KARMAN_CONSTANT, &
                                GRASS_HEIGHT_CM, MEASUREMENT_HEIGHT_CM, &
                                BARE_SOIL_HEIGHT_CM, PI, ALBEDO_PONDING
      implicit none

      type(pm_inputs_t),  intent(in)  :: inputs
      type(pm_outputs_t), intent(out) :: outputs

      ! Local variables - atmospheric
      real(8) :: lambda, delta, ea, ed, vpd, gamma, palt, rho, cp
      real(8) :: tavk, tmnk, tmxk, tkv

      ! Local variables - radiation
      real(8) :: rns, rnc, rnw, rnp, rnl, relssd
      real(8) :: gs, gc, gw
      real(8) :: sinld, cosld, dayl, radial, dec, aob
      real(8) :: sunrise, sunset, startrec, endrec

      ! Local variables - aerodynamic
      real(8) :: chplant, zmeasw
      real(8) :: ud, rac, ras, raw, rss
      real(8) :: zm, zh, d, zom, zoh
      real(8) :: dgrass, zomgrass, zact, dact, zomact, fact, fmeas

      ! Local variables - Penman-Monteith terms
      real(8) :: gammos, gammoc, gammow, gammop
      real(8) :: etaers, etaerc, etaerw, etaerp
      real(8) :: etrads, etradc, etradw, etradp

      ! Local variables - PMdirect
      real(8) :: vcover, laieff

      ! Outputs are auto-initialised by intent(out) on the derived type
      ! (defaults in pm_outputs_t).

      ! ========================================================================
      ! 1. PREPROCESSING: Unit conversions and crop height determination
      ! ========================================================================

      zmeasw = 100.0d0 * inputs%altw  ! Convert wind measurement height to cm

      ! Determine effective crop height
      if (.not. inputs%flCropEmergence) then
          chplant = GRASS_HEIGHT_CM
      else
          if (inputs%swcf == 1 .or. inputs%swcf == 3) then
              chplant = GRASS_HEIGHT_CM
          else
              chplant = max(inputs%ch, 0.1d0)
          endif
      endif

      ! ========================================================================
      ! 2. ATMOSPHERIC PROPERTIES
      ! ========================================================================

      ! Temperature conversions [K]
      tavk = inputs%tav + 273.15d0
      tmnk = inputs%tmn + 273.15d0
      tmxk = inputs%tmx + 273.15d0

      ! Atmospheric pressure at elevation [kPa]
      palt = 101.3d0 * ((tavk - 0.0065d0*inputs%alt) / tavk)**5.26d0

      ! Latent heat of vaporization [MJ/kg]
      lambda = 2.501d0 - 0.002361d0*inputs%tav

      ! Saturation vapour pressure [kPa]
      if (inputs%flmetdetail) then
          ea = 0.611d0 * exp(17.27d0*inputs%tav / (inputs%tav + 237.3d0))
      else
          ea = 0.3055d0 * (exp(17.27d0*inputs%tmn / (inputs%tmn + 237.3d0)) + &
                          exp(17.27d0*inputs%tmx / (inputs%tmx + 237.3d0)))
      endif

      ! Measured vapour pressure (capped at saturation)
      ed = min(inputs%hum, ea)

      ! Vapour pressure deficit [kPa]
      vpd = ea - ed

      ! Slope of vapour pressure curve [kPa/C]
      delta = 4098.0d0 * ea / (inputs%tav + 237.3d0)**2

      ! Psychrometric constant [kPa/C]
      gamma = 0.00163d0 * palt / lambda

      ! Atmospheric density [kg/m3]
      tkv = tavk / (1.0d0 - 0.378d0*ed/palt)
      rho = 3.486d0 * palt / tkv

      ! Specific heat of moist air [kJ/kg/C]
      cp = 622.0d0 * gamma * lambda / palt

      ! ========================================================================
      ! 3. WIND SPEED AND AERODYNAMIC RESISTANCE
      ! ========================================================================

      ! Day wind speed [m/s], avoid zero
      ud = max(inputs%win, 0.0001d0)

      ! Adjust wind speed for height differences
      if (chplant > MEASUREMENT_HEIGHT_CM .or. zmeasw > MEASUREMENT_HEIGHT_CM) then
          dgrass = 2.0d0/3.0d0 * GRASS_HEIGHT_CM
          zomgrass = 0.123d0 * GRASS_HEIGHT_CM
          fmeas = log((1.0d4 - dgrass) / zomgrass) / &
                  log((zmeasw - dgrass) / zomgrass)

          zact = max(chplant, 200.0d0)
          dact = 2.0d0/3.0d0 * chplant
          zomact = 0.123d0 * chplant
          fact = log((zact - dact) / zomact) / log((1.0d4 - dact) / zomact)

          ud = ud * fact * fmeas
      endif

      ! Aerodynamic parameters for crop
      zm = max(chplant, MEASUREMENT_HEIGHT_CM)
      zh = zm
      d = 2.0d0/3.0d0 * chplant
      zom = 0.123d0 * chplant
      zoh = 0.1d0 * zom

      ! Aerodynamic resistance for crop (dry and wet) [s/m]
      rac = log((zm - d)/zom) * log((zh - d)/zoh) / KARMAN_CONSTANT**2 / ud
      raw = rac

      ! Aerodynamic resistance for bare soil [s/m]
      d = 2.0d0/3.0d0 * BARE_SOIL_HEIGHT_CM
      zom = 0.123d0 * BARE_SOIL_HEIGHT_CM
      zoh = 0.1d0 * zom
      ras = log((zm - d)/zom) * log((zh - d)/zoh) / KARMAN_CONSTANT**2 / ud

      ! Surface resistance of soil [s/m]
      if (inputs%swdivide == 1) then
          rss = inputs%rsoil  ! PMdirect partitioning
      else
          rss = 0.0d0
      endif

      ! Modified psychrometric constants [kPa/C]
      gammos = gamma * (1.0d0 + rss/ras)
      gammoc = gamma * (1.0d0 + inputs%rsc/rac)
      gammow = gamma * (1.0d0 + inputs%rsw/raw)

      ! ========================================================================
      ! 4. RADIATION CALCULATIONS
      ! ========================================================================

      ! Net shortwave radiation [MJ/m2/d]
      rns = (1.0d0 - inputs%rcs) * inputs%rad / 1.0d6
      rnc = (1.0d0 - inputs%albedo) * inputs%rad / 1.0d6
      rnw = (1.0d0 - inputs%albedo) * inputs%rad / 1.0d6
      rnp = (1.0d0 - ALBEDO_PONDING) * inputs%rad / 1.0d6

      ! Extraterrestrial radiation and daylength
      if (inputs%flmetdetail) then
          ! Sub-daily: calculate astronomical parameters
          radial = PI / 180.0d0
          dec = -asin(sin(23.45d0*radial) * &
                      cos(2.0d0*PI*dble(inputs%daynr+10) / 365.0d0))

          sinld = sin(radial*inputs%lat) * sin(dec)
          cosld = cos(radial*inputs%lat) * cos(dec)
          aob = sinld / cosld

          ! Daylength calculation with polar circle handling
          if (aob < -1.0d0) then
              dayl = 0.0d0
              outputs%warning_code = 1   ! Polar circle, 0 hours
          else if (aob > 1.0d0) then
              dayl = 24.0d0
              outputs%warning_code = 2   ! Polar circle, 24 hours
          else
              dayl = 12.0d0 * (1.0d0 + 2.0d0*asin(aob)/PI)
          endif

          sunrise = 0.5d0 - dayl / 48.0d0
          sunset = 0.5d0 + dayl / 48.0d0
          startrec = dble(real(inputs%irecord-1) / real(inputs%nmetdetail))
          endrec = dble(real(inputs%irecord) / real(inputs%nmetdetail))
      else
          ! Daily: use provided parameters
          ! Note: sinld, cosld would come from astro() call in wrapper
          sinld = inputs%dsinbe
          cosld = sqrt(max(0.0d0, 1.0d0 - sinld**2))
          startrec = 0.0d0
          endrec = 1.0d0
          sunrise = 0.0d0
          sunset = 1.0d0
      endif

      ! Net longwave radiation [MJ/m2/d]
      relssd = max(min((inputs%atmtr - inputs%a) / inputs%b, 1.0d0), 0.0d0)
      rnl = 4.9d-9 * 0.5d0 * (tmxk**4 + tmnk**4) * &
            (0.34d0 - 0.14d0*sqrt(ed)) * (0.1d0 + 0.9d0*relssd)

      ! Soil heat flux [MJ/m2/d]
      if (inputs%flmetdetail) then
          if ((startrec + endrec)/2.0d0 > sunrise .and. &
              (startrec + endrec)/2.0d0 < sunset) then
              ! Daytime
              gs = 0.1d0 * (rns - rnl)
              gc = 0.1d0 * (rnc - rnl)
              gw = 0.1d0 * (rnw - rnl)
          else
              ! Nighttime
              gs = 0.5d0 * (rns - rnl)
              gc = 0.5d0 * (rnc - rnl)
              gw = 0.5d0 * (rnw - rnl)
          endif
      else
          ! Daily: negligible net flux
          gs = 0.0d0
          gc = 0.0d0
          gw = 0.0d0
      endif

      ! ========================================================================
      ! 5. PENMAN-MONTEITH EQUATION - STANDARD METHOD
      ! ========================================================================

      ! Aerodynamic term [mm/d]
      etaers = (86.4d0/lambda) * (1.0d0/(delta + gammos)) * (rho*cp*vpd/ras)
      etaerc = (86.4d0/lambda) * (1.0d0/(delta + gammoc)) * (rho*cp*vpd/rac)
      etaerw = (86.4d0/lambda) * (1.0d0/(delta + gammow)) * (rho*cp*vpd/raw)

      ! Radiation term [mm/d]
      etrads = delta/(delta + gammos) * (rns - rnl - gs) / lambda
      etradc = delta/(delta + gammoc) * (rnc - rnl - gc) / lambda
      etradw = delta/(delta + gammow) * (rnw - rnl - gw) / lambda

      ! Total potential rates [mm/d]
      outputs%es0 = max(0.0d0, etaers + etrads)
      outputs%et0 = max(0.0d0, etaerc + etradc)
      outputs%ew0 = max(0.0d0, etaerw + etradw)

      ! ========================================================================
      ! 6. PENMAN-MONTEITH DIRECT PARTITIONING (PMdirect)
      ! ========================================================================

      if (inputs%swdivide == 1) then
          ! Vegetation cover fraction
          vcover = 1.0d0 - exp(-inputs%kdif * inputs%kdir * inputs%lai)

          ! Adjust aerodynamic resistances for partial cover
          if (vcover > 1.0d-6) then
              rac = rac / vcover
          else
              rac = 1.0d12
          endif
          raw = rac

          if ((1.0d0 - vcover) > 1.0d-6) then
              ras = ras / (1.0d0 - vcover)
          else
              ras = 1.0d12
          endif

          ! Effective LAI for resistance scaling
          laieff = inputs%lai / (0.3d0*inputs%lai + 1.2d0)

          ! Modified psychrometric constants with resistance scaling
          gammos = vlarge
          if (ras > small) gammos = gamma * (1.0d0 + rss/ras)

          gammoc = vlarge
          if ((rac*laieff) > small) gammoc = gamma * (1.0d0 + inputs%rsc/(rac*laieff))

          gammow = vlarge
          if ((raw*laieff) > small) gammow = gamma * (1.0d0 + inputs%rsw/(raw*laieff))

          gammop = vlarge
          if (ras > small) gammop = gamma

          ! Aerodynamic terms [mm/d]
          etaers = (86.4d0/lambda) * (1.0d0/(delta + gammos)) * (rho*cp*vpd/ras)
          etaerc = (86.4d0/lambda) * (1.0d0/(delta + gammoc)) * (rho*cp*vpd/rac)
          etaerw = (86.4d0/lambda) * (1.0d0/(delta + gammow)) * (rho*cp*vpd/raw)
          etaerp = (86.4d0/lambda) * (1.0d0/(delta + gammop)) * (rho*cp*vpd/ras)

          ! Radiation terms [mm/d] weighted by cover fraction
          etrads = delta/(delta + gammos) * (rns - rnl - gs) * (1.0d0 - vcover) / lambda
          etradc = delta/(delta + gammoc) * (rnc - rnl - gc) * vcover / lambda
          etradw = delta/(delta + gammow) * (rnw - rnl - gw) * vcover / lambda
          etradp = delta/(delta + gammop) * (rnp - rnl - gs) * (1.0d0 - vcover) / lambda

          ! Direct partitioned rates [mm/d]
          outputs%Edirect     = max(0.0d0, etaers + etrads)
          outputs%Tdirect     = max(0.0d0, etaerc + etradc)
          outputs%Tdirectwet  = max(0.0d0, etaerw + etradw)
          outputs%Edirectpond = max(0.0d0, etaerp + etradp)
      endif

  end subroutine PenMon_calc
```

Notes for the agent:
- The body local declarations (`lambda`, `delta`, …) are unchanged — they're scratch variables, not arguments.
- The eight `if (present(warning_code)) warning_code = N` lines from the old code become unconditional `outputs%warning_code = N` (no `optional`, no `present` check — the field always exists).
- The explicit zero-init of the seven output reals at the start of the old body is **removed** because `intent(out)` on `pm_outputs_t` resets them to the default initialisers (`= 0.0_real64`) automatically. The default for `warning_code` (`= 0`) covers that one too.

- [ ] **Step 3.2: Rewrite `PenMon` wrapper**

In `src/atmosphere/et.f90`, find the existing `PenMon` wrapper (currently starting `subroutine PenMon(logf, swscre, daynr, lat, alt, altw, ...`).

Replace the entire `subroutine PenMon(...) end subroutine PenMon` block with:

```fortran
  !> Penman-Monteith evapotranspiration calculation (wrapper with I/O)
  !!
  !! Thin wrapper around PenMon_calc that calls astro() for the daily
  !! branch and forwards outputs%warning_code to warn().
  subroutine PenMon(inputs, outputs, logf, swscre)
      implicit none

      type(pm_inputs_t),  intent(in)    :: inputs
      type(pm_outputs_t), intent(inout) :: outputs
      integer,            intent(in)    :: logf
        !! Internal number of logbook output file
      integer,            intent(in)    :: swscre
        !! Switch of screen display: 0=none, 1=summary, 2=daynumber

      ! Local variables
      real(8) :: dayl, sinld, cosld
      character(len=200) :: messag
      type(pm_inputs_t)  :: inputs_local

      ! Make a mutable working copy so astro() can backfill the
      ! daily-branch astronomical fields (dayl/sinld/cosld are local;
      ! daylp is already in inputs).
      inputs_local = inputs

      ! Call astro() for daily radiation if needed
      if (.not. inputs_local%flmetdetail) then
          call astro(inputs_local%daynr, inputs_local%lat, inputs_local%rad, &
                     dayl, inputs_local%daylp, sinld, cosld, &
                     inputs_local%difpp, inputs_local%atmtr, inputs_local%dsinbe)
      endif

      ! Call pure calculation core
      call PenMon_calc(inputs_local, outputs)

      ! Handle warnings
      if (outputs%warning_code == 1) then
          messag = 'Warning: latitude above polar circle, daylength = 0hrs'
          call warn('Astro', messag, logf, swscre)
      else if (outputs%warning_code == 2) then
          messag = 'Warning: latitude within polar circle, daylength = 24hrs'
          call warn('Astro', messag, logf, swscre)
      endif

  end subroutine PenMon
```

Notes for the agent:
- The wrapper takes a mutable working copy `inputs_local` so it can update the astronomical fields from `astro()` before passing them to the pure kernel. This preserves the existing behaviour where `astro()` was called *before* `PenMon_calc` and its outputs were threaded in as positional arguments.
- `intent(inout)` on `outputs` is deliberate (spec section "New procedure signatures"): the kernel uses `intent(out)`, which reinitialises every field, so the wrapper preserves nothing across the call — but `inout` keeps the contract honest in case future warning-emission logic ever needs to read a prior state.

- [ ] **Step 3.3: Migrate the `meteoday.f90` caller**

In `src/atmosphere/meteoday.f90`, find the existing `call PenMon(...)` block in `ProcessMeteoDay` (currently around line 685, inside the `do 1000 irecord = 1, ndayparts` loop, inside the `elseif (config%meteo%swmetdetail.eq.1 .or. config%meteo%swetr.eq.0)` branch).

First, add the type imports. Locate the existing `use et_mod, only: PenMon, reduceva` near the top of `ProcessMeteoDay` (currently around line 519) and change it to:

```fortran
    use et_mod, only: PenMon, reduceva, pm_inputs_t, pm_outputs_t
```

Then, declare `pmi` and `pmo` as locals at the top of `ProcessMeteoDay`, alongside `rcs`. Locate:

```fortran
    real(8)  rcs
    data     rcs/0.15d0/
```

and add immediately after it:

```fortran
    type(pm_inputs_t)  :: pmi
    type(pm_outputs_t) :: pmo
```

Now replace the `call PenMon(...)` block. Find the block (it spans nine continuation lines):

```fortran
        call PenMon (logf,swscre,tc_daynr,state%cfg%meteo%lat,state%cfg%meteo%alt,state%cfg%meteo%altw,angstroma, &
                     angstromb,rcs,rad,state%atmosphere%Tav,hum,win,state%crop%common%rsc, &
                     state%crop%es0,state%crop%et0,state%crop%ew0, &
                     state%crop%swcf,state%crop%common%ch, &
                     state%crop%flCropEmergence,daylp,tc_flmetdetail,irecord, &
                     config%meteo%nmetdetail,state%crop%common%albedo,tmn,tmx,state%crop%common%rsw,difpp,dsinbe,atmtr, &
                     Edirect,Tdirect,Tdirectwet,rsoil,state%cfg%meteo%swdivide, &
                     state%crop%kdif,state%crop%kdir, &
                     state%crop%lai,Edirectpond)
```

Replace with:

```fortran
        ! Pack PM inputs.
        pmi%daynr           = tc_daynr
        pmi%irecord         = irecord
        pmi%nmetdetail      = config%meteo%nmetdetail
        pmi%flmetdetail     = tc_flmetdetail
        pmi%flCropEmergence = state%crop%flCropEmergence
        pmi%swcf            = state%crop%swcf
        pmi%swdivide        = state%cfg%meteo%swdivide

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
        Edirect        = pmo%Edirect
        Tdirect        = pmo%Tdirect
        Tdirectwet     = pmo%Tdirectwet
        Edirectpond    = pmo%Edirectpond
```

Notes for the agent:
- `Edirect`, `Tdirect`, `Tdirectwet`, `Edirectpond` are module-level locals from `MeteoVars` (declared at the top of `meteoday.f90`) — they remain locals consumed later in `ProcessMeteoDay`. The unpack keeps them populated.
- The pack/call/unpack sits **inside the `do 1000 irecord = 1, ndayparts` loop, inside the `elseif (config%meteo%swmetdetail.eq.1 .or. config%meteo%swetr.eq.0)` branch**, replacing the existing call in place. The other branch (daily + `swetr=1`) does not call PenMon and is unchanged.

- [ ] **Step 3.4: Build the full project**

Run: `pixi run build-linux`
Expected: clean build with no errors. If errors appear, fix them before continuing.

- [ ] **Step 3.5: Run the pFUnit suite**

Run: `pixi run test-pfunit 2>&1 | tail -40`
Expected: all tests pass, including `test_penmon_midsummer_midlatitude_characterization`. The three `@assertEqual` checks must produce byte-identical-within-1e-4 values to the previous build.

If `test_penmon_midsummer_midlatitude_characterization` fails on a value mismatch:
- The kernel body has drifted from identifier-only substitution. Diff the new `PenMon_calc` body against git HEAD and confirm only `inputs%X` / `outputs%X` substitutions, no arithmetic changes.
- Do not "fix" the expected values. The whole point is byte-identical preservation.

- [ ] **Step 3.6: Commit**

```bash
git add src/atmosphere/et.f90 src/atmosphere/meteoday.f90
git commit -m "refactor(et): wrap PenMon args in pm_inputs_t / pm_outputs_t

Collapses the 40-positional-arg signatures of PenMon_calc and PenMon
into derived types. The sole production caller in ProcessMeteoDay
packs inputs from existing state/config/locals immediately before the
call and unpacks outputs immediately after. No physics changes;
characterization test values are byte-identical.

Spec: docs/superpowers/specs/2026-05-22-penman-monteith-derived-types-design.md"
```

---

## Task 4: Verification gate (check-fast)

**Not a commit.** This is the project's standard per-task regression gate (per the `feedback_per_task_regression_gate` playbook entry — pFUnit alone misses global-default regressions).

- [ ] **Step 4.1: Run check-fast**

Run: `pixi run check-fast 2>&1 | tail -30`
Expected: pFUnit green + the four fast regression cases (hupselbrook, surfacewater, salinitystress, grassgrowth) all match their golden outputs.

If any regression case fails:
- Inspect the regression diff (the test runner prints the failing variable and cell).
- The most likely cause is a typo in the identifier substitution inside `PenMon_calc` (e.g. `inputs%tav` vs `inputs%Tav` — Fortran is case-insensitive so this shouldn't bite, but `inputs%hum` vs `inputs%win` could). Diff `et.f90` against HEAD~1 and audit each `inputs%X` reference against the old positional name.
- Do not push commits until check-fast is green.

- [ ] **Step 4.2: Report completion to the user**

Once check-fast is green, the refactor is done. Report:
- Three commits landed.
- check-fast green.
- `pm_inputs_t` and `pm_outputs_t` are the new public contract for PenMon.
- No further work needed for this spec.
