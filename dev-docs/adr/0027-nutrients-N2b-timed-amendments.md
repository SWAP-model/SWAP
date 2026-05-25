---
title: "ADR 0027 — [nutrients] N2b: timed amendments via CSV companion"
date: 2026-05-08
status: accepted
---

# ADR 0027: [nutrients] N2b — timed amendments via CSV companion

## Context

ADR 0025 (N1) wired the crop-side nutrient parameters; ADR 0026
(N2a) wired the soil-side initial state and `SorpCoef`. The
remaining gap before the runtime gate at `tillage.f90:73` can
be lifted is the timed soil management events — fertilizer
applications and manure spreading driven by the legacy
`<project>.sme` file.

Pre-flight inventory found:
- `<project>.sme` legacy schema: `smedate` (date), `MatNum`
  (1..20), `Dosagekgha` (0..500000), `VolatFraction` (legacy
  range 0..500000 — typo; physically a fraction in [0, 1]).
  Up to `maxamn = 1000` rows.
- The pre-collapse `SoilManagement(1)` block (deleted in SS-C
  step 3, commit `67f2545`) sorted events by date, converted
  kg/ha → kg/m², and grouped same-day dosages into
  `TimeAmend(:)`, `NuAmend(:)`, `iamend(:,:)`.
- The runtime path (`SoilManagement(3)` → `Wofost_SoilAmendents`)
  consumes those grouped globals; no code change needed there.

User feedback during brainstorming: "the toml would explode" —
inline `[[nutrients.events]]` arrays were rejected in favour
of a CSV companion. Same shape as SSDI events (ADR 0022).

## Decision

`nutrients_config_t` gains an optional
`character(len=:), allocatable :: events_file` field. The TOML
reader parses `[nutrients].events_file = "amendments.csv"`.
When non-empty, the adapter resolves the path against
`pathwork` and stages a 4-column CSV
(`date,material,amount_kgha,volat_fraction`) via
`read_csv_table` (column 1's ISO date auto-converts to
days-since-1900 per ADR 0012).

`apply_nutrients_events` (new public helper inside
`config_to_variables_mod`, called from the existing
`apply_nutrients`) does:

1. Default to `namend = 0`, `isme = 1` (no amendments).
2. If `events_file` allocated and `len_trim > 0`: load CSV,
   validate per-row ranges, sort in-place by date (bubble
   sort), populate `MatNum(:)`, `Amend(:)` (kg/ha → kg/m²),
   `VolaFrac(:)`, then group same-day dosages into
   `TimeAmend(:)`, `NuAmend(:)`, `iamend(:,:)`, set `namend`,
   `isme = 1`.
3. Sort + group algorithm is a verbatim port of the deleted
   pre-collapse `SoilManagement(1)` legacy block.

The `volat_fraction` validator range tightens from legacy's
typo'd `[0, 500000]` to physically meaningful `[0, 1]`.

Note: `get_optional_string_with_default` always allocates
`cfg%events_file` (defaults to empty string when the TOML
key is missing), so the no-events check uses
`len_trim(events_file) == 0`, not `.not. allocated`.

## Schema

```toml
[nutrients]
sorp_coef   = 0.005
events_file = "amendments.csv"

[nutrients.initial]
fom  = [0.5, 0.3, 0.2, 0.1, 0.5, 0.3, 0.2, 0.1]
bio  = 0.4
hum  = 8.0
cnh4 = 0.001
cno3 = 0.005
```

```
date,material,amount_kgha,volat_fraction
2003-04-15,10,100.0,0.05
2003-06-20,1,25000.0,0.10
2003-06-20,3,500.0,0.02
```

## Consequences

- With N1 + N2a + N2b in place, every soil-side and crop-side
  nutrient input flows through TOML + CSV. No legacy file
  reads remain on the nutrient-init path.
- `apply_nutrients` is called unconditionally per N2a; the
  `events_file` branch is also unconditional but takes the
  no-op path when the field is empty. check-full
  byte-identical.
- N3 (next sub-arc) lifts `tillage.f90:73`'s runtime
  stub-error and adds a regression case with populated
  `[nutrients]`. Verifying byte-identical against the legacy
  binary is N3's correctness gate.
- `volat_fraction` range tightening is a behavioural change.
  Legacy values outside `[0, 1]` were physically meaningless;
  ADR 0027 documents the intentional tightening.

## What N2b does NOT do

- **Lift `tillage.f90:73`.** N3.
- **Add a regression case.** N3.
- **Material-property overrides** (legacy `<project>.smm`).
  Out of scope across the whole `[nutrients]` umbrella.
- **`error_collection_t` consolidation.** Adapter validation
  routes through `fatalerr_collected` — matches SSDI mode-0
  precedent. Consolidation is a separate architectural arc.
- **Range-validate `events_file` shape at validate-time.**
  Per-row CSV validation lives in the adapter (file isn't
  parsed at config-validate time).

## Acceptance

- `pixi run -e test test-pfunit` → `Ok: 1, Fail: 0`; suite
  count +6 (607 from 601 pre-N2b).
- `pixi run -e test check-full` → `5 passed, 0 failed`
  (byte-identical CSVs).
- `grep -n "subroutine apply_nutrients_events" src/io/toml/config_to_variables.f90`
  → one match.
- CSV fixture committed at
  `tests/unit/io/toml/fixtures/nutrients_events_small.csv`.
- `tillage.f90:73` UNCHANGED — `grep -n "flCropNut.*not.*allowed" src/crop/tillage.f90` returns one match.
- ADR 0027 committed.

## Related

- ADR 0008 — error collection over fatalerr.
- ADR 0012 — CSV companion input files.
- ADR 0022 — SSDI events CSV companion (mirror shape).
- ADR 0024 — `dtutil.f90` shim — same architectural thread.
- ADR 0025 — [nutrients] N1 (crop-side adapter).
- ADR 0026 — [nutrients] N2a (soil-side initial state).
- Future: ADR 0028 (N3: lift runtime gate + regression case).
