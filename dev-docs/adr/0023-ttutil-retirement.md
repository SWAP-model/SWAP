---
title: "ADR 0023 — TTutil retired from the SWAP build"
date: 2026-05-07
status: accepted
---

# ADR 0023: TTutil retired from the SWAP build

## Context

ADRs 0019 / 0021 / 0022 retired the TTutil-based **data readers**
from the production runtime, leaving only utility functions
(`getun` / `getun2` / `fopens` / `delfil` for unit-number
management and file I/O) and the rerun-mechanism plumbing
(`rdsets` / `rdfrom` / `rddtmp`) still depending on TTutil.

Continuing to ship TTutil as a `subproject('ttutil', ...)` of
the meson build is overhead: 170 source files, ~5 MB, of which
SWAP uses fewer than ten functions — none of them genuinely
TTutil-specific.

## Decision

Drop TTutil from the build entirely. Replace its utility
functions with native Fortran intrinsics wrapped in a small
project module; retire the rerun mechanism (parameter sweeps
move out-of-band); rename the shim implementation of `FatalERR`
to canonical.

## Phase chronology

| Phase | What landed |
|---|---|
| A | New `file_io_mod` (`src/io/file_io.f90`) with `file_open`/`file_delete`/`file_exists` wrappers + pFUnit suite. |
| B | ~125 call sites of `getun`/`getun2`/`fopens`/`delfil` across 9 files converted to the wrapper. Per-file commits. |
| C | `swap_main.f90` simplified to `call swap(0,1); call swap(0,2); call swap(0,3)` triple. `rdsets`/`rdfrom`/`rddtmp` plumbing dropped. Reruns moved out-of-band. |
| D | Dead TTutil-related local declarations dropped (residue of earlier reader-deletion arcs). |
| E | `fatalerr_shim.f90` → `fatalerr.f90` rename; `ttutil_dep` and the `subprojects/ttutil/` tree removed from meson. TTutil date utilities (`dtdpst`, `dtdpar`, `dtardp`, `dtleap`, `dtnow`) and string utilities (`lowerc`, `upperc`, `addstr`, `words`, `decrea`, `ifindi`) ported to `src/core/dtutil.f90` as canonical native-Fortran implementations. |

## Consequences

- The SWAP binary no longer links against TTutil. The
  `subprojects/ttutil/` source tree (170 files, ~5 MB) is gone
  from the repo.
- Uniform file-open error reporting via `swap_log` at every I/O
  site.
- `swap_main.f90` is ~30 LoC shorter; parameter sweeps are an
  external concern.
- `fatalerr.f90` is the canonical fatal-error handler; the
  "shim" framing in its predecessor is gone.
- `src/core/dtutil.f90` provides canonical native-Fortran
  implementations of the TTutil date and string utility
  routines that remain in active use throughout `src/`, without
  touching the ~75 call sites.
- `[nutrients]` reactivation (separate umbrella) starts from a
  TTutil-free baseline.

## Acceptance

- `grep -rE "\b(rdinit|rdsdor|rdsinr|rdfdor|rdinqr|rdatim|rdsdou|rdadou|rdfint|rdftim|rdfinr|rdscha|rdinar|rdinne|rdsets|rdfrom|rddtmp|getun|getun2|fopens|delfil)\b" src/` → no active-code matches (only `!!` doc comments referencing legacy ranges).
- `grep -rn "ttutil\|TTutil\|TTUTIL" src/ meson.build tests/unit/meson.build` → comment/documentation references only; no active code or build references.
- `subprojects/ttutil/` does not exist; `subprojects/ttutil.wrap` does not exist.
- `pixi run -e test test-pfunit` → `Ok: 1, Fail: 0`.
- `pixi run -e test check-full` → `5 passed, 0 failed`.
- check-full output CSVs are byte-identical to the pre-arc baseline.

## Related

- ADR 0018 (fatalerr shim) — the shim framing retired here.
- ADR 0019 (legacy readers retired) — TTutil data-reader removal.
- ADR 0021 (tillage TOML port) — closed one of the two surviving
  reader exceptions.
- ADR 0022 (SSDI TOML port) — closed the other.
- Future: `[nutrients]` umbrella — separate spec for nutrient
  reactivation.
