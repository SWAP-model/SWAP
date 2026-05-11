# SS-ATM Phase 0 Audit — A-0.1 Findings

**Date:** 2026-05-11
**Branch:** `refactor/atmosphere-state`
**Plan:** `docs/superpowers/plans/2026-05-11-atmosphere-state-migration.md`
**Design:** `docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md` (D12)

---

## Item 1 — `swsublim` typed-config coverage

**Status: COVERED — no gap.**

`swsublim` is fully covered by the typed-config stack:

| Layer | Location | Detail |
|---|---|---|
| Config struct | `src/config/soil_config.f90:59` | `soil_frost_config_t :: swsublim = 0` (default 0) |
| Validation | `src/config/soil_config.f90:402` | `check_int_enum(self%swsublim, [0,1], ...)` |
| TOML reader | `src/io/toml/read_soil_toml.f90:113` | `get_optional_int_with_default(frost, 'swsublim', ..., 0, ...)` |
| Adapter seed | `src/io/toml/config_to_variables.f90:708` | `swsublim = config%soil%frost%swsublim` |
| Global declaration | `src/core/variables.f90:1065` | `integer swsublim` |
| Usage site | `src/atmosphere/snow.f90:104` | `if (swsublim .eq. 0) then` |

The field lives in `soil_frost_config_t` (not `meteorology_config_t`), which is the correct home (snow sublimation suppression was introduced as Adaptation 3 for PEARL-MACRO and is gated alongside other frost/snow soil switches). No corrective action needed.

**pFUnit coverage:** `tests/unit/io/toml/test_read_soil_toml.pf:108` asserts `swsublim == 1` when `soil_frost_full.toml` sets it. Round-trip is tested.

**check-full regression coverage:** None of the 6 check-full cases (`1.hupselbrook` through `6.surfacewater`) sets `swsublim = 1`. The `swsublim = 1` (sublimation active) code path is untested at the integration level. This is a **known gap** but is explicitly excluded from A-0.1 scope — the design doc (line 97) lists "Adding regression coverage for `swsnow=1`" (which gates the entire snow path, including sublimation) as a non-goal for this arc. No action required here.

---

## Item 2 — `spev`/`saev` adapter seeding

**Status: COVERED via `initialize.f90` zero-init — no TOML adapter path needed.**

`spev` and `saev` are per-event accumulators for the Boesten-Stroosnijder soil evaporation reduction method:

- `saev` — cumulative actual evaporation (`src/core/variables.f90:976`)
- `spev` — cumulative potential evaporation (`src/core/variables.f90:981`)

**Seeding chain:**

| Layer | Location | Detail |
|---|---|---|
| Global zero-init | `src/core/initialize.f90:452` | `saev = 0.0d0` |
| Global zero-init | `src/core/initialize.f90:457` | `spev = 0.0d0` |
| Per-event reset | `src/atmosphere/et.f90:639-641` | reset to 0 on tillage/new rain event inside `reduceva` |

Neither `spev` nor `saev` appears anywhere in `src/config/` or `src/io/toml/` — they are not config-driven. Zero-init at startup is correct: the Boesten-Stroosnijder accumulators are built up by `reduceva` within each dry-period event and are reset when a new wetting event begins (see `et.f90:639–641`). There is no warm-restart path for these fields (they are also not in the `swinco=3` warm-restart block at `config_to_variables.f90:558–570`). The design plan confirms these are atmosphere-arc fields (`src/atmosphere/et.f90` is the owner).

**Note for A-1.4:** The `state` argument added to `et.f90` in A-1.4 will expose `spev` and `saev` as `state%atmosphere%spev` / `state%atmosphere%saev`. The dual-write for these fields is via `reduceva`'s new signature (`subroutine reduceva(task, nrai, state)`; see design D7). The `initialize.f90` zero-inits will remain correct until A-2.6 retires the globals; at that point `atmosphere_init` zeroes the state fields.

---

## Item 3 — `config_to_variables.f90` ldwet/ssnow init-seeds (A-2.6 hazard)

**Status: DOCUMENTED — compile-error hazard confirmed for A-2.6.**

The `swinco=3` warm-restart block in `src/io/toml/config_to_variables.f90` seeds two atmosphere-owned globals:

| Variable | Line | Expression |
|---|---|---|
| `ssnow` | 561 | `ssnow = config%soil%initial%ssnow` |
| `ldwet` | 565 | `ldwet = config%soil%initial%ldwet` |

Additional conditional on line 569:

```fortran
if (config%meteo%snow%swsnow /= 1) ssnow = 0.0d0
```

**Hazard for A-2.6:** When A-2.6 retires `ssnow` and `ldwet` from `variables.f90` (removing `use variables, only: ..., ssnow, ldwet, ...` in `config_to_variables.f90`), these three assignment lines will become compile errors. They must be replaced with dual writes:

```fortran
! A-2.6 migration target:
ssnow           = config%soil%initial%ssnow          ! line 561 — retire + write state%atmosphere%ssnow
state%atmosphere%ssnow = config%soil%initial%ssnow
ldwet           = config%soil%initial%ldwet           ! line 565 — retire + write state%atmosphere%ldwet
state%atmosphere%ldwet = config%soil%initial%ldwet

if (config%meteo%snow%swsnow /= 1) ssnow = 0.0d0    ! line 569 — retire + write state%atmosphere%ssnow
if (config%meteo%snow%swsnow /= 1) state%atmosphere%ssnow = 0.0_real64
```

**Context:** `ssnow` is a flat top-level field of `atmosphere_state_t` (design plan, "Instantaneous (11)" group). `ldwet` is in the "Per-event / tillage-reset (3)" group. Both will be initialized by `atmosphere_init` to 0 at startup; the `swinco=3` path then overwrites with warm-restart values. The `config_to_variables` call happens inside `read_swap_toml` before the main time loop, so the assignment order is: `atmosphere_init` (zeros) → `config_to_variables` (swinco=3 overwrites) → time loop.

**Source field locations:**

- `config%soil%initial%ssnow` — `src/config/soil_config.f90:71`
- `config%soil%initial%ldwet` — `src/config/soil_config.f90:74`

Both fields have validation at `soil_config.f90:422–429`.

---

## Summary

| Audit item | Gap found? | Action taken |
|---|---|---|
| `swsublim` typed-config + adapter | No gap | Audit-only (doc) |
| `spev`/`saev` adapter seeding | No gap — zero-init correct | Audit-only (doc) |
| `config_to_variables` ldwet/ssnow seeds | Hazard documented | A-2.6 must migrate lines 561, 565, 569 |
| `swsublim=1` check-full coverage | Gap (known non-goal) | Excluded per design doc line 97 |

**Phase 0 gap closures: zero.** All three items passed the audit. The A-2.6 hazard was already flagged in the design self-review; this document records the exact line numbers.
