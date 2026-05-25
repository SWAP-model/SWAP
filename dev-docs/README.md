# Dev documentation (internal)

This directory holds development-only material that used to live under `docs/`
but is no longer published on the FORD-generated public site. It's kept here
so contributors can still find it via repository search and editor file
navigation.

The public site source is at `docs/`. The FORD config (`docs.md`) does not
scan this directory.

## Contents

### Architecture decisions

- `adr/` — 42+ Architecture Decision Records. Each ADR captures a non-obvious
  choice made during the modernization arcs and the reasoning behind it. The
  ADR carries custom metadata keys (`status:`, `supersedes:`,
  `supersedes-portion-of:`) that FORD doesn't recognise — that's why these
  files don't belong on the public site.

### Plans, specs, and notes

- `superpowers/specs/` — design documents brainstormed before implementation.
- `superpowers/plans/` — implementation plans corresponding to each spec.
- `superpowers/notes/` — ad-hoc per-arc notes.

### Phase 4 archive

- `archive/2026-phase-4/` — frozen specs, plans, and per-reader audits from
  the Phase 4 / Phase 4f-extend modernization (completed 2026-05-05; tags
  `rescue/phase-4f-extend-complete` and `rescue/phase-4f-extend-followups`).

### Historical reports

- `phase-4-modernization-summary.md` — capstone summary of the Phase 4 arc.
- `coverage-baseline.md` — Phase 3 line/branch coverage snapshot. Coverage
  is tracked, not gated (see `adr/0006-coverage-tracked-not-gated.md`).
