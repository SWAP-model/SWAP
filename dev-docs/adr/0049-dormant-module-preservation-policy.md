# ADR 0049 — Dormant-module preservation policy (+ T1-B dead-code sweep)

**Status:** Accepted (2026-07-03)
**Arc:** T1-B (first Tier-1 arc of the 2026-07-03 modernization review)
**Relates to:** ADR 0032 (AgeTracer extraction), ADR 0040 (macropore retirement),
ADR 0047 (jongvanlier + legacy_state deletion). Does not supersede any ADR.

## Context

The 2026-07-03 modernization review's §4 estimated "~9 K lines of dead code to
remove," naming the `src/**/dormant/` modules and a list of "dead" state fields.
Discovery for the T1-B arc revised that sharply, and the revision is the point of
this ADR:

1. The `dormant/` modules are **not dead code** — they are a deliberate,
   documented preservation of legacy capabilities the TOML migration did not
   re-wire (alternative drought/oxygen formulations, tabulated hydraulics,
   AgeTracer, water-table/regrid/mass-balance utilities). Each `dormant/` dir
   carries a `README.md` with a per-module reactivation checklist; the files are
   excluded from the meson build (zero binary cost).
2. Ten of the twelve "dead state fields" named in the review are in fact **live**
   (read in compute or output). Only two were genuinely dead.

Without a written policy, this ambiguity recurs every review: someone greps,
finds an uncompiled module or an unread-looking field, and proposes deleting it.

## Decisions

### 1. Dormant-module preservation policy (formalized)

A **dormant module** is legacy behavior preserved as source but not compiled:

- It lives under a feature folder's `dormant/` subdirectory.
- It is **excluded from `meson.build`** (and `tests/unit/meson.build`) — it does
  not build and cannot run.
- Its `dormant/README.md` records, per file: the feature, the switch that would
  dispatch it, the `fatalerr_collected` stub that currently guards that switch,
  and the reactivation checklist.
- The **raw** legacy body is independently preserved on the `legacy/swap-4.2.0`
  branch; the `dormant/` copy adds the *modern-integration* restoration notes.
- **Deleting a dormant module requires an ADR** (as macropore did, ADR 0040, and
  jongvanlier, ADR 0047). Grep-driven "looks unused" is not sufficient grounds.
- Dormant bodies are frozen at extraction time and are **not** kept in lockstep
  with active-module refactors; harmonizing against the current data model
  (`state%…` / `config%…`, post-`variables.f90`, post-`state%cfg`) is part of the
  reactivation work, not ongoing maintenance.

Current dormant inventory (2026-07-03): `soilwater/dormant/{sptabulated,
watertable, checkmassbal, regrid}.f90`, `solute/dormant/agetracer.f90`,
`crop/dormant/oxygenrepro.f90`.

### 2. Dead-code removed in this arc (byte-identical)

- **`stepnr()`** (`arrayutils.f90`) — a complete public function, self-documented
  `@warning not called from any SWAP code`, with zero callers in `src/` or
  `tests/`. Removed (definition + public export).
- **`pegwl`, `npegwl`** (`soilwater_state.f90`) — written in `waterbalance.f90`
  but read nowhere (compute, output, or bindings). Removed, along with their
  write-only assignments. Byte-identity holds because every removed RHS is
  side-effect-free: the `level()` call that fed `pegwl` takes `state` as
  `intent(in)` (pure) and is still called for the live `gwl` at
  `waterbalance.f90:87`. The sibling `bpegwl` (node index) is the only live
  output of that block and is retained. Removing the dead `level()` call also
  drops a ~59 K-calls/run pure evaluation (minor bonus).

### 3. Correction of record — fields that are LIVE (do not re-delete)

A full-tree read/write grep (including `associate` aliases, output writers, and
BMI/CAPI bindings) confirmed the following review-§4 candidates are **live** and
must be kept: `evp` (`soilhydraulics.f90:101`), `runonarr`/`flrunon`
(`boundtop.f90:53-54`), `alfaw_layer` (`soilhydraulics.f90:1349`), `ksatexm`
(~8 sites: drainage/surfacewater/boundtop/frozencond/soilhydraulics(utils)/
tillage), `cofani` (drainage/surfacewater/frozencond), `flksatexm`
(`tillage.f90:50`), `bpegwl` (waterbalance perched-gwl loop), `psilt`/`pclay`
(tillage + temperature).

## Consequences

- The "is this live?" ambiguity has a durable answer: `dormant/` = intentional
  (see its README + this ADR); anything else unread is a genuine-dead candidate
  that still needs read/write verification before removal.
- The three `dormant/README.md` files were refreshed to match reality (stale
  `jongvanlier.f90` row removed; `variables.f90`/`state%cfg` restoration steps
  annotated as targeting the typed records).
- Net `src/` change is small (~1 function + 2 fields + their writes); the arc's
  value is correctness (no false-dead deletions) and clarity, not LoC.
- `check-fast` byte-identical (4/4) and pFUnit green (833) after the removals.
</content>
