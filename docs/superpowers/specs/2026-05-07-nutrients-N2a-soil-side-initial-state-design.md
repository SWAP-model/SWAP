---
title: "[nutrients] N2a — Soil-side initial state + SorpCoef (ADR 0026 candidate)"
author: Mateusz Zawadzki
date: 2026-05-07
status: draft
---

# [nutrients] N2a — Soil-side initial state + SorpCoef

Second sub-arc of the `[nutrients]` umbrella. Adds a top-level
`[nutrients]` TOML block carrying the soil-side initial pool
concentrations and the `SorpCoef` sorption coefficient. Wires
the values to the existing legacy globals via a new adapter.
Retires the broken `<project>_nut.end` diagnostic dump in
`SoilManagement(7)`.

**Does NOT lift the runtime gate at `tillage.f90:73`.** The
soil-side initial state flows through, but the simulation
still aborts at runtime if `flCropNut=true` reaches
`DoTillage`. That gate is lifted in N3 once N2a + N2b
(amendments) are both in place.

## Goal

Reach a state where:

1. A TOML config that sets `[nutrients]` (top-level) loads
   without error and the typed config flows through to the
   legacy `variables`-module globals: `FOM_t(1..8)`, `Bio_t`,
   `Hum_t`, `cNH4_t`, `cNO3_t`, `SorpCoef`.
2. The `[nutrients]` block is **optional**. When absent, the
   adapter falls back to: pools at zero (current behaviour;
   natural-mineralization-only simulation), `SorpCoef = 0.0`
   (no sorption — physically conservative; matches
   uninitialised-variable behaviour seen pre-arc).
3. `SoilManagement(7)`'s `_nut.end` diagnostic dump is gone
   (the case body collapses to `return`, mirroring SS-C step 3
   of the legacy-readers-deletion arc). The dump was
   template-driven by reading `<project>.snp` line by line —
   without `.snp` files in the modern flow, the read would
   fail at runtime if `flCropNut` were ever true.
4. `pixi run -e test test-pfunit` passes; new pFUnit suite
   covers parser + validator + adapter behaviour.
   `pixi run -e test check-full` is unchanged (no regression
   case enables `flCropNut` yet — that's N3).
5. ADR 0026 captures the decision and the unfinished business
   N2b/N3 close.

End-state: `[nutrients]` initial state flows through to legacy
globals; `SorpCoef` (genuine uninitialised-variable bug
discovered during N2 brainstorm) gets a deterministic default;
broken `_nut.end` write-path retired.

## Non-goals

- **Timed soil management events (amendments).** The legacy
  `<project>.sme` mechanism (smedate, MatNum, MatAmount,
  VolaFrac for fertilizer applications) is N2b's scope. Without
  N2b, simulations run with **zero amendments** — natural
  mineralization only. That's a useful baseline behaviour for
  testing.
- **Material-properties overrides** (`<project>.smm` legacy).
  The 17 hardcoded materials in `Wofost_SoilParameters` (cattle
  manure → spruce needles, with their AppAge / OrgMatFrac /
  OrgNFrac / NH4NFrac / NO3NFrac defaults) stay hardcoded.
  TOML port deferred indefinitely; user-tunable on demand.
- **Lifting `tillage.f90:73`.** N3's job after both N2a and N2b
  land. With no amendments wired (N2b pending), a flCropNut=true
  simulation would run from zero pools with no fertilizer
  inputs — useful but not the full picture.
- **Adding a regression case** with `flcropnut=true`. Cannot
  meaningfully do this until the runtime gate is lifted (N3).
- **Replacing the `_nut.end` dump with a CSV-style writer.**
  Future enhancement; not blocked on N2a.

## Schema

### New types in `src/config/nutrients_config.f90`

```fortran
module nutrients_config_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use error_mod, only: error_collection_t
   implicit none
   private

   public :: nutrients_initial_t
   public :: nutrients_config_t

   !> Initial concentrations of soil organic-matter and N pools.
   !! Units: kg/m^2 for FOM/Bio/Hum (depth-integrated); kg/m^3 for
   !! cNH4/cNO3 (per-water-volume concentrations). Legacy ranges
   !! (rdsdor [0, 1000]) preserved.
   type :: nutrients_initial_t
      real(real64) :: fom(8) = 0.0_real64    !! FOM_t(1..8) — eight FOM pools
      real(real64) :: bio    = 0.0_real64    !! Bio_t — microbial biomass
      real(real64) :: hum    = 0.0_real64    !! Hum_t — humus
      real(real64) :: cnh4   = 0.0_real64    !! cNH4_t — ammonium concentration
      real(real64) :: cno3   = 0.0_real64    !! cNO3_t — nitrate concentration
   end type nutrients_initial_t

   !> [nutrients] top-level block.
   type :: nutrients_config_t
      logical      :: present   = .false.        !! true iff [nutrients] was supplied
      real(real64) :: sorp_coef = 0.0_real64     !! SorpCoef — sorption coefficient (m^3/kg). Default 0 = no sorption.
      type(nutrients_initial_t) :: initial
   contains
      procedure :: validate => nutrients_config_validate
      procedure :: finalize => nutrients_config_finalize
   end type nutrients_config_t

contains

   subroutine nutrients_config_validate(self, errors)
      ! ... see §Validation rules ...
   end subroutine

   subroutine nutrients_config_finalize(self, errors)
      class(nutrients_config_t), intent(inout) :: self
      type(error_collection_t),  intent(inout) :: errors
      ! No finalization needed; included for interface uniformity.
      return
   end subroutine

end module nutrients_config_mod
```

`swap_config_t` (in `src/config/swap_config.f90`) gains:

```fortran
type(nutrients_config_t) :: nutrients
```

### TOML shape

```toml
[nutrients]
sorp_coef = 0.005

[nutrients.initial]
fom  = [0.5, 0.3, 0.2, 0.1, 0.5, 0.3, 0.2, 0.1]   # FOM_t(1..8) kg/m^2
bio  = 0.4                                          # Bio_t kg/m^2
hum  = 8.0                                          # Hum_t kg/m^2
cnh4 = 0.001                                        # cNH4_t kg/m^3
cno3 = 0.005                                        # cNO3_t kg/m^3
```

The entire `[nutrients]` block is optional. Absent block →
`present = .false.`, defaults intact (zero pools, zero
sorption).

`[nutrients.initial]` is optional within `[nutrients]`. A user
might want to set `sorp_coef` while leaving pools at zero.

`fom` accepts a 1- to 8-element array. Missing entries default
to 0.0 (matches behaviour of unused legacy pool slots —
`Wofost_SoilParameters` sets `nf = 8` but the higher-indexed
pools are reference-rate-constant variants that may not all
carry initial mass in every site).

## Validation rules

In `nutrients_config_validate`, fired only when `self%present`
is `.true.`:

| Rule | Error context |
|---|---|
| `sorp_coef >= 0.0` | `nutrients.sorp_coef` |
| `initial.fom(i)` in `[0.0, 1000.0]` for each populated `i` | `nutrients.initial.fom` |
| `initial.bio` in `[0.0, 1000.0]` | `nutrients.initial.bio` |
| `initial.hum` in `[0.0, 1000.0]` | `nutrients.initial.hum` |
| `initial.cnh4` in `[0.0, 1000.0]` | `nutrients.initial.cnh4` |
| `initial.cno3` in `[0.0, 1000.0]` | `nutrients.initial.cno3` |

The legacy `rdsdor` enforced `[0, 1000]` for all 12 values —
preserve that contract. No upper bound on `sorp_coef`; values
above ~0.1 are unphysical but the legacy code didn't enforce it
and unusual scenarios (synthetic test cases, etc.) might use
high values intentionally.

When `self%present == .false.`, validator stays silent —
defaults are tolerated unconditionally.

No cross-section rules at the typed-config level. A future
"`flCropNut=true` requires `[nutrients]` present" cross-validator
can land in N3 alongside the runtime-gate lift.

## TOML reader

New module `src/io/toml/read_nutrients_toml.f90`:

```fortran
module read_nutrients_toml_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, get_value, len
   use toml_field_helpers_mod, only: get_optional_real_with_default
   use nutrients_config_mod, only: nutrients_config_t
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_nutrients_toml

contains

   !> Populate `cfg` from the optional top-level [nutrients] table on `doc`.
   !! Missing block is benign — leaves defaults; sets cfg%present = .false..
   subroutine read_nutrients_toml(doc, cfg, errors)
      type(toml_table), pointer, intent(in)    :: doc
      type(nutrients_config_t),  intent(inout) :: cfg
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), pointer :: nut_tbl, init_tbl
      type(toml_array), pointer :: fom_arr
      integer :: stat, n, i
      real(real64) :: v

      if (.not. associated(doc)) return

      nut_tbl => null()
      call get_value(doc, 'nutrients', nut_tbl, requested=.false., stat=stat)
      if (stat /= 0 .or. .not. associated(nut_tbl)) return

      cfg%present = .true.
      call get_optional_real_with_default(nut_tbl, 'sorp_coef', cfg%sorp_coef, &
                                          0.0_real64, 'nutrients.sorp_coef', errors)

      ! [nutrients.initial] sub-table
      init_tbl => null()
      call get_value(nut_tbl, 'initial', init_tbl, requested=.false., stat=stat)
      if (stat == 0 .and. associated(init_tbl)) then
         ! fom is a 1- to 8-element array; missing entries stay at 0.0.
         fom_arr => null()
         call get_value(init_tbl, 'fom', fom_arr, requested=.false., stat=stat)
         if (stat == 0 .and. associated(fom_arr)) then
            n = len(fom_arr)
            if (n > 8) then
               call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                  'nutrients.initial.fom: at most 8 entries (legacy FOM_t(1..8) cap)', &
                  'nutrients.initial.fom')
               n = 8
            end if
            do i = 1, n
               call get_value(fom_arr, i, v, stat=stat)
               if (stat == 0) cfg%initial%fom(i) = v
            end do
         end if

         call get_optional_real_with_default(init_tbl, 'bio',  cfg%initial%bio,  &
                                             0.0_real64, 'nutrients.initial.bio',  errors)
         call get_optional_real_with_default(init_tbl, 'hum',  cfg%initial%hum,  &
                                             0.0_real64, 'nutrients.initial.hum',  errors)
         call get_optional_real_with_default(init_tbl, 'cnh4', cfg%initial%cnh4, &
                                             0.0_real64, 'nutrients.initial.cnh4', errors)
         call get_optional_real_with_default(init_tbl, 'cno3', cfg%initial%cno3, &
                                             0.0_real64, 'nutrients.initial.cno3', errors)
      end if
   end subroutine read_nutrients_toml

end module read_nutrients_toml_mod
```

Called from `load_swap_config` after the existing top-level
section parses (general / simulation / soil / drainage /
crop / etc.):

```fortran
call read_nutrients_toml(doc, config%nutrients, errors)
```

## Adapter

`src/io/toml/config_to_variables.f90` gets a new helper:

```fortran
subroutine apply_nutrients(cfg)
   use, intrinsic :: iso_fortran_env, only: real64
   use nutrients_config_mod, only: nutrients_config_t
   use wofost_soil_declarations, only: FOM_t, Bio_t, Hum_t, &
                                        cNH4_t, cNO3_t, SorpCoef
   type(nutrients_config_t), intent(in) :: cfg

   integer :: i

   ! SorpCoef is set unconditionally — fixes the genuine
   ! uninitialised-variable bug discovered during the N2
   ! brainstorm. Default 0.0 = no sorption (matches legacy
   ! behaviour for any case that didn't supply a .snp file).
   SorpCoef = cfg%sorp_coef

   ! Initial pool concentrations. When cfg%present = .false.,
   ! cfg%initial fields are at their default 0.0; the writes
   ! become no-ops (variables already at zero from module init).
   do i = 1, 8
      FOM_t(i) = cfg%initial%fom(i)
   end do
   Bio_t  = cfg%initial%bio
   Hum_t  = cfg%initial%hum
   cNH4_t = cfg%initial%cnh4
   cNO3_t = cfg%initial%cno3
end subroutine apply_nutrients
```

Called unconditionally from `config_to_variables` (no `if
(flSomething)` guard) — the "absent block" case is handled by
defaults.

## Runtime change in `src/crop/management_soil.f90`

`SoilManagement(7)`'s body reads `<project>.snp` line-by-line as
a template and writes a partially-substituted copy to
`<project>_nut.end`. The case body is several dozen lines (see
`management_soil.f90:587-657`). Replace with:

```fortran
case (7)
   ! [nutrients] N2a (ADR 0026): the legacy diagnostic dump used
   ! <project>.snp as a template (read each line; if a known key
   ! match, write the current pool value to <project>_nut.end;
   ! else echo). Without .snp files in the modern flow, the read
   ! would fail. The dump is retired here. A future arc can add
   ! a clean CSV-style nutrient-pool dump if anyone asks.
   !
   ! The flCropExt write-back at the top of the legacy case (7)
   ! is also retired — that branch was driven by a separate
   ! mechanism (cropext) tied to the N-P-K case-1 init that's
   ! now also gone.
   return
```

Drop unused locals from the `SoilManagement` subroutine that
were specific to case (7) (e.g., `snp`, `oup`, `cropext`,
`filnam` if no longer referenced after the case-7 collapse).
The `flCropExt` global stays as-is — referenced elsewhere.

`SoilManagement` cases 2-6 retain their bodies (they're the
nutrient runtime physics; not touched by N2a).

## Tests

New pFUnit suite `tests/unit/io/toml/test_nutrients_config.pf`:

| Test | Asserts |
|---|---|
| `test_parse_full_block` | TOML with `[nutrients]` + `[nutrients.initial]` populated; assert all 12 values + `sorp_coef` round-trip; `present == .true.` |
| `test_parse_missing_block_keeps_defaults` | TOML without `[nutrients]`; `present == .false.`, all initial pools at 0.0, `sorp_coef == 0.0` |
| `test_parse_partial_fom_array` | `fom = [0.5, 0.3, 0.2]` (3 elements); assert `fom(1..3)` populated, `fom(4..8) == 0.0` |
| `test_parse_oversized_fom_array` | `fom = [...9 entries...]`; expect parse error |
| `test_validate_silent_when_absent` | `present = .false.`, all defaults; `validate(errs)` → `errs%count() == 0` |
| `test_validate_passes_clean` | `present = .true.`, all values in range; `errs%count() == 0` |
| `test_validate_rejects_negative_sorp_coef` | `sorp_coef = -0.1`; expect error |
| `test_validate_rejects_fom_above_1000` | `fom(3) = 1500.0`; expect error |
| `test_validate_rejects_negative_bio` | `bio = -0.5`; expect error |
| `test_apply_populates_globals` | Build config with all values; call `apply_nutrients`; assert `FOM_t(1..8)`, `Bio_t`, `Hum_t`, `cNH4_t`, `cNO3_t`, `SorpCoef` match |
| `test_apply_defaults_zero_pools_zero_sorp` | `present = .false.`; call `apply_nutrients`; assert all globals are 0.0 |

Suite registered in `tests/unit/meson.build` `pf_files` and
`tests/unit/testSuites.inc`. New module added to
`pfunit_extra_sources`.

`check-full`: 5/5 unchanged. No regression case sets
`flCropNut=true`, so even with the adapter populating soil-side
state, the nutrient runtime path stays gated.

## ADR

`docs/adr/0026-nutrients-N2a-soil-side-initial-state.md` — sub-arc
N2a of the `[nutrients]` umbrella. Sections:

- Context (umbrella overview, three sub-arcs N1/N2a/N2b/N3,
  why N2 was split into N2a + N2b).
- Decision (top-level `[nutrients]` block; optional with
  zero-defaults; SorpCoef default 0.0; `_nut.end` dump retired).
- What N2a does NOT do (amendments → N2b; runtime gate → N3;
  material-overrides → indefinite).
- Consequences (one new top-level config section; fixes the
  SorpCoef uninitialised-variable bug; case-7 dump gone).
- Acceptance criteria.
- Related (ADR 0008, ADR 0024 architectural-direction note,
  ADR 0025 N1, future ADR 0027 N2b, future ADR 0028 N3).

`docs/adr/index.md` updated.

## Acceptance criteria

- [ ] `pixi run -e test build-linux` clean.
- [ ] `pixi run -e test test-pfunit` exits 0; suite count grows
      by 11 (the eleven tests above).
- [ ] `pixi run -e test check-full` exits 0 with `5 passed,
      0 failed` (byte-identical CSV outputs).
- [ ] `grep -nE "snp\b" src/crop/management_soil.f90` → no
      matches (case-7 body retired).
- [ ] `grep -nE "_nut\.end" src/crop/management_soil.f90` → no
      matches.
- [ ] A hand-authored fixture in the test suite exercises every
      validator rule.
- [ ] ADR 0026 committed.
- [ ] `tillage.f90:73` UNCHANGED — runtime stub-error stays.
      `grep -n "flCropNut.*not.*allowed" src/crop/tillage.f90` →
      one match.
- [ ] `apply_nutrients` is called from `config_to_variables`
      unconditionally (not gated on `flCropNut`); a future N3
      cross-validator can require `[nutrients]` when
      `flCropNut=true`.

## Commit cadence

~7 commits, each independently buildable + green:

1. `feat(config): add nutrients_config_t types + validator`
2. `test: pFUnit suite for nutrients_config validate`
3. `feat(io): parse [nutrients] in TOML reader; wire into load_swap_config`
4. `test: pFUnit reader tests + adapter placeholder`
5. `feat(io): apply_nutrients adapter populates legacy globals`
6. `refactor(crop): retire SoilManagement(7) _nut.end dump`
7. `docs: ADR 0026 + index update + configuration-schema.md`

Commit (5) replaces the adapter placeholder from (4); commits
(1)+(2) and (3)+(4) may be merged at plan time. Plan-time
decision.

## Risk

- **`SorpCoef = 0.0` default.** This makes a real behavioural
  choice: no sorption when the user doesn't set it. The legacy
  code left `SorpCoef` uninitialised — if any case in legacy
  use ran with `flCropNut=true` and uninitialised SorpCoef, its
  output was already non-deterministic. Setting 0.0
  deterministically may produce different numbers vs that
  legacy non-deterministic baseline. Verified at N3 against a
  deliberately-constructed test case (or accepted as a
  documented behavioural change in N3's ADR).
- **Pool-zero default.** Same shape: legacy startup with no
  `.snp` file would fail at runtime; legacy startup with `.snp`
  populated the pools. Modern with no `[nutrients]` block →
  zero pools (defined behaviour). N3's regression case will
  use a fully-populated `[nutrients]` block.
- **Case-7 retirement collateral.** `flCropExt` write-back
  branch at the top of the legacy case (7) is also dropped
  (it's gated by `flCropExt`, which is set elsewhere; verify
  no orphan globals dangle after removal). If `flCropExt`
  is unreachable in the modern flow (likely — it's tied to the
  cropext file mechanism that's no longer wired), drop the
  associated globals at plan time.
- **`apply_nutrients` always called.** Setting `SorpCoef = 0.0`
  on a case where `flCropNut = false` means the pre-existing
  behaviour (`SorpCoef` uninitialised → unspecified value) is
  replaced by deterministic 0.0. None of the five regression
  cases reach a code path that reads `SorpCoef` (all use
  `flCropNut=false`), so check-full stays byte-identical. If
  any test case is added later that reads `SorpCoef` outside
  the nutrient runtime path, this change is the floor.

## Out of scope (reaffirmed)

- N2b (timed soil management events / amendments).
- N3 (lift `tillage.f90:73`; first nutrient regression case).
- Material-properties override (legacy `.smm`).
- Replacing `_nut.end` dump with a CSV-style writer.
- Range-validation tightening on `SorpCoef`.
- Cross-section validator (`flCropNut=true` requires
  `[nutrients]` populated). Future N3 work.
