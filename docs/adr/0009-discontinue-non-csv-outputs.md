---
title: "ADR 0009 — Discontinue non-CSV output formats"
date: 2026-04-27
status: accepted
---

# ADR 0009: Discontinue non-CSV output formats

## Context

SWAP today writes a fan-out of legacy output formats: `.afo` and `.aun`
binary files, `.vap` vapor profile dumps, `.bal` / `.wba` water balance
text reports, `.sba` solute balance, `.blc` balance check, `.drf`
formatted drainage, `.str` stress, `.IRG` irrigation, `.RUM` macropore
runoff, `.SWB` surface water balance, plus initial-state (`swini`) and
end-state (`swend`) dumps and a MODFLOW exchange writer
(`swoutputmodflow`). Each is gated by an `sw*` switch read from the
`.swp` file — about 18 switches in total.

The Phase 4f variables audit flagged these as gaps because the new
`swap_config_t` schema has no slot for them. Adding `output_config_t`
fields for each would re-implement the legacy fan-out in TOML, then
require us to maintain ~20 distinct output writer subroutines indefin-
itely.

The CSV writers (`swcsv`, `swcsv_tz`) cover what these legacy outputs
were used for in modern usage. Test cases in this repo exercise CSV
output exclusively. Maintaining the legacy fan-out for input
compatibility — when no test exercises any of those formats —
amounts to dragging dead code through every refactor.

## Decision

The new TOML pipeline supports **only** CSV output (`swcsv`,
`swcsv_tz`). The 18 legacy output switches are retired.

The retired switches are:

`swafo`, `swaun`, `swvap`, `swbal`, `swwba`, `swsba`, `swblc`,
`swdrf`, `swstr`, `swirg`, `swini`, `swend`, `swheader`, `swcaprise`,
`swcapriseoutput`, `swrum`, `swswb`, `swoutputmodflow`.

Retirement happens in two stages:

**Phase 4e** (this phase) — purely declarative:
- This ADR documents the policy.
- `docs/configuration-schema.md` lists the retired keys explicitly so
  authors know they will be ignored.
- `docs/phase-4f-config-to-variables-audit.md` reclassifies these
  18 entries from `G` (gap) to `RETIRED`.
- A `warn_deprecated_key(routine, key)` helper is added to
  `error_mod` for use by the new TOML readers when they encounter a
  deprecated key.
- **No legacy code is touched.** Existing `.swp` files keep working
  with the legacy reader; existing output writer subroutines still
  exist and run.

**Phase 4f** (next phase) — enforced on the new path:
- The new TOML reader has no schema slot for these keys; if a TOML
  file lists one, the reader appends a non-fatal deprecation warning
  via `warn_deprecated_key` and ignores it.
- `config_to_variables(config)` explicitly forces all 18 globals to
  0 before the legacy execution path runs, so the legacy output
  writers exist as dead code (unreachable but not yet deleted).

**Phase 5+** (later, not gated on this ADR):
- The dead legacy writer subroutines and their state arrays are
  deleted. Estimated ~800–1200 LoC of `src/io/swapoutput.f90`,
  `src/io/macroporeoutput.f90`, and related state in
  `src/core/variables.f90` becomes deletable.

## Consequences

Positive:

- Cuts ~17 entries from the Phase 4f variables-audit gap list.
- Removes a class of "what does this format actually mean?" questions
  from the schema design.
- Sets up ~1000 LoC of legacy code for deletion in a follow-up phase.
- One supported output format means one path to test, one path to
  document, one path to migrate when output requirements change.

Negative:

- Users with workflows that consume `.afo` / `.aun` / `.bal` / etc.
  will need to migrate to CSV. The CSV stream covers the same
  underlying data; only the file format differs.
- `swsublim` (snow sublimation toggle) and `swirg` look like output
  switches by name but `swsublim` is physics-affecting (`snow.f90:99`
  gates the sublimation simulation) and `swirg` is pure output. The
  retirement list distinguishes them — only `swirg` is retired.
  Verified by spot-checks on each switch's call sites.

## Revisit trigger

If a user workflow surfaces that genuinely depends on a non-CSV
legacy format and cannot be served by extending the CSV writer, this
ADR is open for amendment. Until then, the policy stands.
