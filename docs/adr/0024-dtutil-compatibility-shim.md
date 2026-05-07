---
title: "ADR 0024 — `dtutil.f90` as TTutil-API compatibility shim; hoist utility calls out of physics"
date: 2026-05-07
status: accepted
---

# ADR 0024: `dtutil.f90` as TTutil-API compatibility shim; hoist utility calls out of physics

## Context

ADR 0023 retired the TTutil subproject from the SWAP build. The
umbrella spec covered the TTutil **reader** family (`rdinit`,
`rdsdor`, `rdsinr`, `rdfdor`, `rdinqr`, `rdatim`, …), the
**file I/O** primitives (`getun`, `getun2`, `fopens`, `delfil`),
the **rerun mechanism** (`rdsets`, `rdfrom`, `rddtmp`), and
renamed the existing `fatalerr` shim.

The spec **missed** the TTutil date/string utility family:
`DTLEAP`, `DTARDP`, `DTDPAR`, `DTDPST`, `DTNOW`, `LOWERC`,
`UPPERC`, `ADDSTR`, `WORDS`, `DECREA`, `IFINDI`. After Phase D,
roughly 75 sites across 12 source files still called these
functions. With `ttutil_dep` dropped from meson, every one of
those call sites would fail to link.

Two ways forward were available at the start of Phase E:

**(a) Touch every call site.** Replace each `call dtdpst('year-month-day', t, str)` (and friends) with native Fortran or a project-native helper. ~75 scattered edits across 12 files; harder to verify byte-identical CSV outputs because each replacement is locally-shaped.

**(b) Drop-in compatibility shim.** Create a new project-owned file `src/core/dtutil.f90` providing native-Fortran reimplementations with the **exact same signatures** as the TTutil functions. Zero call-site changes. Mirrors the pre-existing `fatalerr_shim.f90` strategy from ADR 0018 — except this shim is the canonical implementation, since TTutil is gone.

The Phase E implementer chose **(b)**. This ADR documents the
decision and sets the future direction.

## Decision

`src/core/dtutil.f90` is the canonical project-owned home for
the 11 TTutil-shaped utility functions still referenced by
`src/`. The file is ~444 LoC of native Fortran 2008. Each
function preserves the upstream signature and behaviour for
inputs that valid SWAP callers actually produce (verified:
check-full byte-identical across all five regression cases).

**Defensive validation is deliberately softened** relative to
upstream: where TTutil's `DTSYS` would have called `FATALERR` on
malformed input (year=0, day > DAYMAX, hour < 0, search bounds
out of array, unknown format descriptor in `DTDPST`, etc.), the
new code silently returns a default. Internal callers in `src/`
always pass valid inputs; treating a defensive `FATALERR` as
load-bearing for those callers would mean propagating that
strict input contract through 75 sites — exactly the work this
ADR set aside.

**`DTDPST`'s tokenizer is replaced with a sequential
first-occurrence `INDEX`-based replacer.** Behaviourally
equivalent for every format string used in `src/` today
(documented in the source review). Edge cases — formats with
alphabetic content beyond the supported tokens — would diverge;
no such caller exists.

`fatalerr_shim.f90` was renamed to canonical `fatalerr.f90` in
the same Phase E commit. Same shape: project-owned, TTutil-API
shape, drop-in.

## Future direction — utility calls don't belong in physics

The compatibility shim solves the build problem but leaves an
architectural smell: physics-focused subroutines (soil
hydraulics, drainage, crop growth, water balance, time control)
make direct calls to date-formatting and string-manipulation
helpers. By layer, the residual call counts are:

| Layer | DTDPST sites | Other dtutil sites |
|---|---|---|
| `src/io/` | 24 | 55 |
| `src/soil/` | 5 | 0 |
| `src/core/` | 4 | 7 |
| `src/atmosphere/` | 2 | 0 |
| `src/drainage/` | 1 | 0 |
| `src/crop/` | 1 | 6 |

`src/io/` calls are appropriate — that's the I/O layer's job.
The 26 physics/core sites are the architectural concern. They
are mostly date-string formatting for log/error messages
(`call dtdpst(...)` → `call warn(messag, ...)`) — the message
construction has crept into the physics layer when it belongs
at the boundary.

**The intended end-state:**

1. Physics subroutines receive parsed-and-validated inputs from
   the adapter layer (typed `swap_config_t`, `swap_state_t`,
   …) and emit results without ever opening a file, formatting
   a date, or splitting a string.
2. Logging from physics goes through `swap_log` (`log_warn`,
   `log_error`), which already accepts pre-formatted strings —
   the date-formatting work happens in the I/O layer or in a
   small message-building helper, not in the soil-water solver.
3. Once every physics-layer caller of `dtutil` is gone,
   `dtutil.f90` shrinks to whatever the I/O layer still needs.
   The TTutil-shaped signatures can be retired in favour of
   project-native APIs (`format_date_iso(t, str)` etc.).
4. Long-term, `src/core/dtutil.f90` may be deletable entirely.

This is captured as a follow-on direction, not a blocking arc.
The `[nutrients]` umbrella is the next priority; this cleanup
can run in parallel or wait.

## Consequences

- **TTutil-the-dependency is gone** (no subproject, no
  `ttutil_dep` in meson, no vendored 5 MB tree). That goal of
  ADR 0023 is achieved.
- **TTutil-shaped APIs persist** in `src/core/dtutil.f90` until
  the physics-layer cleanup happens. Future maintainers see
  legacy-style names (`DTDPST`, `LOWERC`, `IFINDI`) in the
  codebase and may assume TTutil is still in play.
- **Build is fully self-contained** for the date/string utility
  surface. No upstream dependency to track, version-pin, or
  patch.
- **Net code accounting:** dropped ~5 MB of vendored TTutil;
  added ~550 LoC of project-owned code (`dtutil.f90` 444 +
  `file_io.f90` ~80 + `fatalerr.f90` ~30).
- **Validation softening** could mask bad input from new
  callers. Mitigation: future call sites should be in the I/O
  layer (where validation is the explicit job) or use
  project-native APIs (which can keep strict semantics).

## Acceptance for the future direction

When the physics-layer hoist is done:

- `grep -rE "\b(dtdpst|dtardp|dtdpar|dtnow|dtleap|lowerc|upperc|addstr|words|decrea|ifindi)\b" src/soil src/atmosphere src/drainage src/crop src/heat src/macropore src/solute src/boundary` → no matches.
- `grep -rn "\bcall warn\b" src/soil ...` → uses `swap_log` or pre-formatted strings, not dtutil-via-warn.
- `dtutil.f90` either shrinks substantially (only I/O-layer functions remain) or is replaced by project-native helpers and deleted.
- check-full remains 5/5 with byte-identical CSV outputs.

This work is not gated. It can be done in small per-subsystem
arcs (e.g., one ADR per layer: `src/soil/`, `src/crop/`, …) or
opportunistically as physics modules are touched for other
reasons (a "boy-scout" cleanup). When fully complete, supersede
this ADR with one that captures the deletion of `dtutil.f90`.

## Related

- ADR 0008 — error collection over fatalerr (one abort
  checkpoint; physics shouldn't call fatalerr directly either).
- ADR 0016 — per-rotation crop config cache (future direction:
  physics receives explicit arguments). Same architectural
  thread.
- ADR 0018 — fatalerr shim (predecessor; same drop-in pattern).
- ADR 0023 — TTutil retired from the build (the umbrella that
  triggered this ADR).
- Future: `[nutrients]` umbrella (orthogonal; takes priority).
- Future: per-layer "hoist utilities to I/O boundary" arcs (not
  yet specced).
