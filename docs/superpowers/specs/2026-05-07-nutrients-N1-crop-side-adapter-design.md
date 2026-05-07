---
title: "[nutrients] N1 — Crop-side nutrient adapter (ADR 0025 candidate)"
author: Mateusz Zawadzki
date: 2026-05-07
status: draft
---

# [nutrients] N1 — Crop-side nutrient adapter

First sub-arc of the `[nutrients]` umbrella. Wires the existing
`wofost_nutrient_t` typed config (already defined in
`src/config/cropwofost_config.f90`) to the legacy `variables`
globals via a new per-rotation adapter. Lifts the
`wofost_nutrient_validate` stub-error so configurations with
`flcropnut=true` can be loaded and validated.

**Does NOT lift the runtime gate at `tillage.f90:73`.** The crop-side
schema flows through, but the simulation still aborts at runtime if
`flCropNut=true` reaches `DoTillage`. That gate is lifted in N3
after N2 has wired the soil-side initial state.

## Goal

Reach a state where:

1. A TOML config that sets `crop.rotation.<n>.cropwofost.nutrient.flcropnut = true`
   plus all required nutrient parameters loads without error.
2. The legacy globals consumed by the WOFOST nutrient physics
   (`LRNR`, `LSNR`, `NLAI`, `NLUE`, `NMAXSO`, `NPART`, `NFIXF`,
   `NSLA`, `RNFLV`, `RNFRT`, `RNFST`, `TCNT`, `DVSNLT`, `DVSNT`,
   `RDRNS`, `FNTRT`, `FRNX`, `NMXLV(:)`, `ILNMXL`,
   `FraHarLosOrm_lv/st/so`) are populated from the typed config at
   per-rotation init time.
3. `wofost_nutrient_validate` enforces the same input contract the
   legacy reader implicitly assumed: `nmxlv` non-empty when
   `flcropnut=true`; harvest fractions in `[0, 1]`.
4. `pixi run -e test test-pfunit` passes; new pFUnit suite covers
   adapter + validator behaviour. `pixi run -e test check-full` is
   unchanged (no regression case enables nutrients yet — that's N3).
5. ADR 0025 captures the decision and the unfinished business N2/N3
   close.

End-state for this sub-arc: the typed-config layer fully
round-trips nutrient inputs to the legacy globals. The runtime
gate at `tillage.f90:73` still rejects `flCropNut=true` because
soil-side state isn't initialized yet (N2's job).

## Non-goals

- **Soil-side `[nutrients]` TOML block.** Initial state for the
  FOM/Bio/Hum/NH4/NO3 pools, `DryBD`, `SorpCoef`, etc. — N2's
  scope. Hardcoded scientific defaults in `Wofost_SoilParameters`
  remain untouched.
- **Lifting `tillage.f90:73`.** The DoTillage runtime stub-error
  for `flCropNut=true` stays. Lifting it would let a partially-wired
  config slip into runtime where SoilManagement(2..7) reads
  uninitialised soil-pool state. N3's job, after N2.
- **Adding a regression case** with `flcropnut=true`. Cannot
  meaningfully do this until the runtime path is unblocked. N3.
- **Touching `wofost_nutrient_t` shape.** The type already has all
  the right fields (`lrnr`, `lsnr`, `nlai`, `nlue`, `nmaxso`,
  `npart`, `nfixf`, `nsla`, `rnflv`, `rnfrt`, `rnfst`, `tcnt`,
  `dvsnlt`, `dvsnt`, `rdrns`, `fntrt`, `frnx`, `nmxlv` allocatable,
  `frahar_los_orm_lv/st/so`). N1 only adds validation rules; no
  schema changes.
- **CSV companions for `nmxlv`.** Per-DVS table is up to 30
  rows in the legacy reader (`rdadou(... 30, ILNMXL)`). Inline
  TOML array is fine at that size. (Mirror of tillage tables
  decision.)

## Schema

No new types. Existing `wofost_nutrient_t` in
`src/config/cropwofost_config.f90`:

```fortran
type :: wofost_nutrient_t
   logical      :: flcropnut = .false.
   real(real64) :: lrnr      = 0.0_real64
   real(real64) :: lsnr      = 0.0_real64
   real(real64) :: nlai      = 0.0_real64
   real(real64) :: nlue      = 0.0_real64
   real(real64) :: nmaxso    = 0.0_real64
   real(real64) :: npart     = 0.0_real64
   real(real64) :: nfixf     = 0.0_real64
   real(real64) :: nsla      = 0.0_real64
   real(real64) :: rnflv     = 0.0_real64
   real(real64) :: rnfrt     = 0.0_real64
   real(real64) :: rnfst     = 0.0_real64
   real(real64) :: tcnt      = 0.0_real64
   real(real64) :: dvsnlt    = 0.0_real64
   real(real64) :: dvsnt     = 0.0_real64
   real(real64) :: rdrns     = 0.0_real64
   real(real64) :: fntrt     = 0.0_real64
   real(real64) :: frnx      = 0.0_real64
   real(real64), allocatable :: nmxlv(:)
   real(real64) :: frahar_los_orm_lv = 0.0_real64
   real(real64) :: frahar_los_orm_st = 0.0_real64
   real(real64) :: frahar_los_orm_so = 0.0_real64
contains
   procedure :: validate => wofost_nutrient_validate
   procedure :: finalize => wofost_nutrient_finalize
end type
```

## TOML

```toml
[[crop.rotation]]
type    = "wofost"
file    = "maize.crp.toml"

# In maize.crp.toml:
[wofost.nutrient]
flcropnut = true
lrnr      = 0.5
lsnr      = 0.5
nlai      = 1.0
nlue      = 1.1
nmaxso    = 0.0176
npart     = 1.0
nfixf     = 0.0
nsla      = 0.5
rnflv     = 0.004
rnfrt     = 0.0048
rnfst     = 0.0015
tcnt      = 10.0
dvsnlt    = 1.0
dvsnt     = 0.8
rdrns     = 0.05
fntrt     = 0.15
frnx      = 0.5
nmxlv     = [0.06, 0.06, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.0176]
frahar_los_orm_lv = 0.5
frahar_los_orm_st = 0.5
frahar_los_orm_so = 0.0
```

When `flcropnut = false` (default), the entire `[wofost.nutrient]`
block can be omitted; field defaults are zero, validator silent.

## Validation rules

In `wofost_nutrient_validate`, fired only when `self%flcropnut == .true.`:

| Rule | Error context |
|---|---|
| `nmxlv` allocated and non-empty | `cropwofost.nutrient.nmxlv` |
| `size(nmxlv) <= 30` (legacy `rdadou` cap) | `cropwofost.nutrient.nmxlv` |
| `frahar_los_orm_lv` in `[0.0, 1.0]` | `cropwofost.nutrient.frahar_los_orm_lv` |
| `frahar_los_orm_st` in `[0.0, 1.0]` | `cropwofost.nutrient.frahar_los_orm_st` |
| `frahar_los_orm_so` in `[0.0, 1.0]` | `cropwofost.nutrient.frahar_los_orm_so` |

Removed: the existing stub-error block:
```fortran
if (self%flcropnut) then
   call errors%append(ERR_VALIDATION_CROSS_FIELD, &
      'cropwofost.nutrient.flcropnut=.true. (N-P-K nutrient model) ' // &
      'not yet supported in the TOML pipeline; use the legacy executable.', &
      'cropwofost.nutrient')
end if
```

The 17 scalars (`lrnr`, `lsnr`, ..., `frnx`) are NOT range-checked
in N1 — the legacy reader (`rdsdou`) didn't enforce ranges either,
so the legacy contract is "any real". Adding sensible bounds is
out of scope; can be a future tightening when nutrient regression
fixtures exist (post-N3).

When `flcropnut == .false.`, the validator stays silent (defaults
are tolerated unconditionally).

## Adapter

`src/crop/cropwofost_init.f90` — extend with one new public
procedure and one wire-up call:

```fortran
public :: apply_cropwofost_nutrient

subroutine apply_cropwofost_nutrient(cfg)
   use variables, only:                                   &
        LRNR, LSNR, NLAI, NLUE, NMAXSO, NPART, NFIXF,    &
        NSLA, RNFLV, RNFRT, RNFST, TCNT,                 &
        DVSNLT, DVSNT, RDRNS, FNTRT, FRNX,               &
        NMXLV, ILNMXL,                                     &
        FraHarLosOrm_lv, FraHarLosOrm_st, FraHarLosOrm_so
   type(wofost_nutrient_t), intent(in) :: cfg

   integer :: n

   LRNR   = cfg%lrnr
   LSNR   = cfg%lsnr
   NLAI   = cfg%nlai
   NLUE   = cfg%nlue
   NMAXSO = cfg%nmaxso
   NPART  = cfg%npart
   NFIXF  = cfg%nfixf
   NSLA   = cfg%nsla
   RNFLV  = cfg%rnflv
   RNFRT  = cfg%rnfrt
   RNFST  = cfg%rnfst
   TCNT   = cfg%tcnt
   DVSNLT = cfg%dvsnlt
   DVSNT  = cfg%dvsnt
   RDRNS  = cfg%rdrns
   FNTRT  = cfg%fntrt
   FRNX   = cfg%frnx

   n = 0
   if (allocated(cfg%nmxlv)) n = size(cfg%nmxlv)
   ILNMXL = n
   if (n > 0) NMXLV(1:n) = cfg%nmxlv(1:n)

   FraHarLosOrm_lv = cfg%frahar_los_orm_lv
   FraHarLosOrm_st = cfg%frahar_los_orm_st
   FraHarLosOrm_so = cfg%frahar_los_orm_so
end subroutine apply_cropwofost_nutrient
```

`NMXLV` in `module variables` is a fixed-size array (likely
`real(8) :: NMXLV(30)`) per the legacy reader's `rdadou(... 30,
ILNMXL)` signature. Verify by `grep -n "NMXLV" src/core/variables.f90`
during implementation; adapt the slice assignment if the actual
declaration is allocatable.

### Wire-up

Inside `cropwofost_init_from_config(cfg, icrop, FraDeceasedLvToSoil)`,
after the existing wofost-param copies, add:

```fortran
if (cfg%nutrient%flcropnut) call apply_cropwofost_nutrient(cfg%nutrient)
```

Place it near the end of the subroutine, after the per-rotation
output-file initialization (`outbalcropOM1`, `outbalcropOM2`,
`outbalcropN`) — the order doesn't matter functionally (those
calls don't read the nutrient globals), but keeps the file
narratively coherent.

## Tests

New pFUnit suite `tests/unit/io/toml/test_apply_cropwofost_nutrient.pf`:

| Test | Asserts |
|---|---|
| `test_apply_populates_legacy_globals` | Build minimal `wofost_nutrient_t` with all 17 scalars + non-empty `nmxlv` + 3 harvest fractions; call `apply_cropwofost_nutrient`; assert each global matches expected value, `ILNMXL == size(nmxlv)`, `NMXLV(1:n)` matches input. |
| `test_validate_passes_clean_flcropnut_true` | Construct `wofost_nutrient_t` with `flcropnut=true`, valid `nmxlv`, valid harvest fractions; `validate(errors)`; assert `errors%count() == 0`. |
| `test_validate_rejects_flcropnut_true_empty_nmxlv` | Same but `.not. allocated(nmxlv)`; expect at least one error. |
| `test_validate_rejects_flcropnut_true_oversized_nmxlv` | `size(nmxlv) == 31`; expect error. |
| `test_validate_rejects_flcropnut_true_out_of_range_harvest` | `frahar_los_orm_lv = 1.5`; expect error. |
| `test_validate_silent_flcropnut_false` | `flcropnut=false`, all other fields at defaults; `validate(errors)`; assert `errors%count() == 0`. |

Suite registered in `tests/unit/meson.build` `pf_files` and
`tests/unit/testSuites.inc`.

`pfunit_extra_sources` already includes `cropwofost_init.f90` and
`cropwofost_config.f90`; no new entries needed.

`check-full`: 5/5 unchanged. No existing regression case sets
`flcropnut=true`, so the new adapter is unreachable from
production at the end of N1.

## ADR

`docs/adr/0025-nutrients-N1-crop-side-adapter.md` — sub-arc N1 of
the `[nutrients]` umbrella. Sections:

- Context (umbrella overview, three sub-arcs N1/N2/N3, why split).
- Decision (lift the validator stub; add the adapter; runtime gate
  stays at `tillage.f90:73` until N3).
- Consequences (typed-config now flows through; soil-side and
  runtime-gate work outstanding).
- Related (ADR 0016 — per-rotation cache; ADR 0021/0022 — same
  shape; future ADR 0026 — soil-side; future ADR 0027 — runtime
  reactivation).

`docs/adr/index.md` updated.

## Acceptance criteria

- [ ] `pixi run -e test build-linux` clean.
- [ ] `pixi run -e test test-pfunit` exits 0; suite count grows by
      6 (the six tests above).
- [ ] `pixi run -e test check-full` exits 0 with `5 passed,
      0 failed`.
- [ ] `grep -n "not yet supported in the TOML pipeline.*nutrient" src/config/cropwofost_config.f90`
      → no match (the stub-error string is gone).
- [ ] A hand-authored fixture in the test suite exercises every
      validator rule.
- [ ] ADR 0025 committed.
- [ ] `tillage.f90:73` UNCHANGED — the runtime stub-error stays.
      `grep -n "flCropNut.*not.*allowed" src/crop/tillage.f90` →
      one match.

## Commit cadence

~5-6 commits, each independently buildable + green:

1. `feat(config): drop wofost_nutrient stub-error; add validation rules`
2. `test: pFUnit suite for wofost_nutrient validation`
3. `feat(crop): apply_cropwofost_nutrient adapter (cropwofost_init)`
4. `test: pFUnit adapter tests for apply_cropwofost_nutrient`
5. `feat(crop): wire apply_cropwofost_nutrient into cropwofost_init_from_config`
6. `docs: ADR 0025 — [nutrients] N1 crop-side adapter`

(Commits 1+2 may merge into one; same for 3+4. Plan-time decision.)

## Risk

- **`NMXLV` storage class.** If `module variables` declares it
  allocatable rather than fixed-size, the adapter must allocate
  on first use. Verified at implementation time.
- **`ILNMXL` global.** Some legacy code may treat `ILNMXL` as a
  global "active length"; the adapter sets it to `size(cfg%nmxlv)`,
  which is the legacy `rdadou` semantic.
- **No regression-test signal until N3.** The adapter populates
  legacy globals that nothing reads until `tillage.f90:73` is
  lifted. The pFUnit suite covers the adapter contract; runtime
  byte-identity vs the legacy binary is verified in N3.
- **Validator changes don't break existing cases.** All five
  current regression cases default to `flcropnut=false` (or omit
  the block entirely). Validator stays silent at `flcropnut=false`,
  so existing cases are unaffected.

## Out of scope (reaffirmed)

- Soil-side `[nutrients]` TOML block (N2).
- Lifting `tillage.f90:73` runtime gate (N3).
- Adding a `flcropnut=true` regression case (N3).
- Range-checking the 17 scalars (legacy didn't; future tightening
  when fixtures exist).
- Touching `cropgrowth.f90`'s deleted `if (flCropNut) outbalcrop*`
  output-file initialization (those stay; they're unrelated to the
  parameter-load path).
