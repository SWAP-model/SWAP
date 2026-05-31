# src/timecontrol

## Responsibility

The simulation clock and schedule. `timecontrol_mod.f90` (`timecontrol_mod`)
owns the time loop bookkeeping: calendar advance, day/year boundaries, output
schedule, and the adaptive `dt` control (reduce-on-non-convergence). It pairs
with the data record `state/timecontrol_state.f90` (`timecontrol_state_t`).

Extracted from `src/core/` per ADR 0048 — it is a compute subsystem, not a
foundation leaf.

## Dependencies

Reads/writes `state%timecontrol`; consumed by `driver/` (the lifecycle
procedures call `timecontrol_init`/`_advance`/`_day_end`/…) and read across the
physics folders via `dt`/`t1900`/`daynr`/`daycum`.
