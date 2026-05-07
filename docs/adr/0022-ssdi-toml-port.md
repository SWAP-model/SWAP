---
title: "ADR 0022 — SSDI parameters ported to [irrigation.ssdi] TOML block"
date: 2026-05-07
status: accepted
---

# ADR 0022: SSDI parameters ported to `[irrigation.ssdi]` TOML block

## Context

ADR 0019's "Update 2026-05-06" closed the umbrella physical
deletion of the legacy fixed-format reader, but two TTutil-based
reads survived as foot-noted exceptions: `Read_Tillage` (closed
by ADR 0021) and `read_ssdi_input` (closed by this ADR).

This ADR ports the SSDI subsystem and retires the standalone
`src/io/checkdate.f90` along with it.

## Decision

SSDI parameters move to a new `[irrigation.ssdi]` sub-table with
an explicit `schedule = 0|1` discriminator. Mode-specific
parameters live in disjoint sub-tables:

- `[irrigation.ssdi.fixed]` (`schedule = 0`) — the date-rate-amount
  table is held in a CSV companion (`events_file = "..."`); CSV-only,
  no inline path. Tables can hit 366 rows and inline TOML becomes
  unreadable at that scale.
- `[irrigation.ssdi.scheduled]` (`schedule = 1`) — trigger-based
  parameters (`sched_type`, thresholds, amount, application rate,
  interval gating) inline.

`swssdi` stays at `[irrigation].swssdi` (least disruption to
existing TOML cases; parallel to the soil.tillage / soil.swtill
pattern from ADR 0021).

`read_ssdi_input` is deleted. A new helper
`apply_irrigation_ssdi` in `config_to_variables.f90` populates
the legacy `_irr`-suffixed globals (`ssdi_date_irr`,
`ssdi_rate_f_irr`, `ssdi_amount_f_irr`, `ssdi_sched_type_irr`,
`ssdi_threshold_irr`, `ssdi_threshold_z_irr`, `ssdi_amount_irr`,
`ssdi_appl_rate_irr`, `sw_interval_irr`, `days_interval_irr`,
`days_counter_irr`, `nod_ssdi_irr`, `nirri_ssdi_irr`, plus
`dt_SSDI_event`, `qssdi`) at config-load time.

`src/io/checkdate.f90` is deleted; the date-window check lives
in `apply_irrigation_ssdi`'s deferred-validation tail, expressed
against `fatalerr_collected` instead of `fatalerr` + `STOP`.

## Schema

```toml
[irrigation]
swssdi = 1

[irrigation.ssdi]
schedule = 0          # 0 fixed-date, 1 scheduled-trigger
ssdi_z   = [-30.0, -50.0]   # cm; both equal for single-depth

# When schedule = 0:
[irrigation.ssdi.fixed]
events_file = "ssdi_events.csv"

# When schedule = 1:
[irrigation.ssdi.scheduled]
sched_type      = 1   # 1=Tred, 2=presh, 3=watc
threshold       = 0.7
threshold_depth = -50.0   # cm; required when sched_type > 1
ssdi_amount     = 10.0    # mm
ssdi_appl_rate  = 5.0     # mm/h
sw_interval     = 0       # 0|1
# days_interval required only when sw_interval = 1
```

CSV companion (`ssdi_events.csv`):
```
date,rate_f,amount_f
2003-04-15,5.0,10.0
2003-05-02,5.0,8.0
```

When `swssdi = 0`, the entire `[irrigation.ssdi]` block can be
omitted.

## Consequences

- `swpfile` global pointer no longer read by `irrigation.f90`. The
  remaining `swpfile`/`logf` plumbing in
  `src/io/toml/config_to_variables.f90` becomes orphan; cleanup
  is the next step (separate, smaller follow-on commit).
- `read_ssdi_input` deleted; `apply_irrigation_ssdi` is the sole
  SSDI-init path.
- `src/io/checkdate.f90` deleted; `fatalerr_collected` is the
  sole validation channel for the date-window check.
- With ADR 0021 + ADR 0022 landed, the only TTutil runtime calls
  remaining in `src/` are utility (`rdsets`/`rdfrom` in
  `swap_main.f90` + `rddtmp` in `swapoutput.f90`) — not data
  readers.

## Tests

pFUnit suites under `tests/unit/io/toml/`:
- `test_read_irrigation_ssdi_toml`
- `test_irrigation_ssdi_validate`
- `test_apply_irrigation_ssdi`

CSV fixture: `tests/unit/io/toml/fixtures/ssdi_events_small.csv`.

No `swssdi = 1` regression case is added; deferred until a
known-good legacy comparator is available.

## Acceptance

- `grep -n "subroutine read_ssdi_input" src/` -> no matches.
- `grep -n "swpfile" src/crop/irrigation.f90` -> no matches.
- `grep -rn "checkdate" src/ --include="*.f90"` -> no matches
  (excluding the historical-note comment in
  `config_to_variables.f90`).
- `src/io/checkdate.f90` does not exist.
- `pixi run -e test test-pfunit` -> `Ok: 1, Fail: 0` (576 tests,
  1 disabled).
- `pixi run -e test check-full` -> `5 passed, 0 failed`.
