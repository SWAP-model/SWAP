# dtutil De-Shim from Physics Layer — Design

**Date:** 2026-05-08
**Status:** accepted (queued for execution)
**Predecessor ADR:** 0024 (`dtutil.f90` as TTutil-API compatibility shim)
**Successor ADR:** 0029

## Goal

Remove all `dtutil.f90` references from the physics layer (`src/crop/`,
`src/soil/`, `src/drainage/`, `src/atmosphere/`, `src/macropore/`),
replacing them with native helpers in non-physics modules. After this
arc, physics subroutines no longer call `dtdpst`, `ifindi`, or any
other TTutil-shape utility.

`dtutil.f90` itself stays — the I/O layer (`src/io/`,
`src/core/timecontrol.f90`, `src/core/swap.f90`,
`src/utils/surfacewaterutils.f90`) continues to use it. Removing
dtutil entirely is out of scope.

## Out of scope

- Restructuring error reporting (the deeper "physics emits codes, I/O
  formats" refactor). For this arc, physics keeps building log/error
  strings inline; we just sever the *dependency* on dtutil-shape
  helpers, not the responsibility of message construction.
- I/O-layer dtutil calls. Those are appropriate per ADR 0024 and are
  not touched.
- Deleting `dtutil.f90`. It stays as the I/O-layer helper.

## Context

ADR 0024 introduced `dtutil.f90` as a same-signature shim for ~11
TTutil utility functions, deferring ~75 call-site edits. Inspection
shows the actual physics-layer footprint is much smaller than that
overall count suggests — the bulk of dtutil callers are in I/O code
where dtutil is the right home for the helper.

The physics-layer call sites:

| File | Calls | Functions used |
|------|-------|----------------|
| `src/soil/soilhydraulics.f90` | 4 | `dtdpst` |
| `src/crop/irrigation.f90` | 4 | `dtdpar`, `dtardp` (all commented-out, dead) |
| `src/crop/cropgrowth.f90` | 3 | `dtdpst` (1), `ifindi` (1 + integer-function decl) |
| `src/atmosphere/meteoday.f90` | 2 | `dtdpst` |
| `src/soil/waterbalance.f90` | 1 | `dtdpst` |
| `src/drainage/surfacewater.f90` | 1 | `dtdpst` |

Net: 9 live `dtdpst` calls + 1 live `ifindi` call + 4 commented-out
dead lines in `irrigation.f90`.

All `dtdpst` calls in the physics layer build a date string for
inclusion in a log or error message (`messag` strings concatenated
into `fatalerr_collected` payloads). They have one of two formats:
`'year-month-day'` or `'year-month-day,hour:minute:seconds'`.

The single `ifindi` call (cropgrowth.f90:993) searches a sorted
integer array (`mayrs`) for `iyear`, returning an index. Trivial to
replace with a native helper.

## Decision

Three changes:

1. **Introduce a native ISO date formatter** in `src/utils/` (a new
   module, e.g. `date_format_mod`). Public surface: a single function
   `format_iso_date(t1900, with_time)` returning a `character(len=19)`
   string. Two output forms based on the boolean:
   - `with_time = .false.` → `'YYYY-MM-DD'`
   - `with_time = .true.` → `'YYYY-MM-DD,HH:MM:SS'`

   Implementation reuses the date-arithmetic core already in
   `dtutil.f90` (the `DPDTTM → DATEA` conversion), but exposes it as
   a clean function rather than a TTutil-shape subroutine.

2. **Move `ifindi` to `src/utils/array_utils.f90`.** That module
   already exists and already hosts `afgen`. Add a pure
   `index_in_sorted_int(values, target, lo, hi)` function with the
   same semantics as `ifindi`. Update the cropgrowth.f90 call site.

3. **Delete the commented-out `dtdpar`/`dtardp` lines** at
   `src/crop/irrigation.f90:113-127`.

Physics-layer files then `use date_format_mod, only: format_iso_date`
and `use array_utils, only: index_in_sorted_int`. No physics file
references dtutil.

## Non-goals

- The output strings are **byte-identical** to what `dtdpst` produced
  for the same input (same format codes, same date arithmetic). This
  is required so check-full passes byte-identical comparisons against
  the existing baselines, since some of these strings flow into log
  outputs that are part of regression fixtures.

## Testing

- pFUnit tests for `format_iso_date` covering both formats and a few
  edge cases (leap day, year boundary, end-of-day rollover).
- pFUnit tests for `index_in_sorted_int` covering hit, miss, boundary
  positions.
- check-full byte-identical regression as the integration gate.

## ADR 0029

`docs/adr/0029-dtutil-deshim-physics-layer.md` captures the decision,
the boundary (physics vs I/O), and the residual `dtutil.f90` purpose
post-arc.
