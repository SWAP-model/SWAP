# Dev documentation (internal)

This directory holds development-only material that used to live under `docs/`
but is no longer published on the FORD-generated public site. It's kept here
so contributors can still find it via repository search and editor file
navigation.

The public site source is at `docs/`. The FORD config (`docs.md`) does not
scan this directory.

## Contents

### Architecture decisions

- `adr/` — 47 Architecture Decision Records (start at `adr/index.md`). Each ADR
  captures a non-obvious choice made during the modernization arcs and the
  reasoning behind it. The ADR carries custom metadata keys (`status:`,
  `supersedes:`, `supersedes-portion-of:`) that FORD doesn't recognise — that's
  why these files don't belong on the public site. ADRs are the durable record:
  a reversed decision gets a new superseding ADR rather than an edit.

### Capstone summaries

- `phase-4-modernization-summary.md` — the Phase 4 arc (input-pipeline
  modernization: TTutil readers → typed TOML, 2026-04-22 → 2026-05-05).
- `post-phase-4-modernization-summary.md` — everything after Phase 4
  (typed-state migrations, globals retirement, `variables.f90` deletion,
  orchestrator/`state%cfg` retirement, polish; 2026-05-12 → 2026-05-31).

### Historical reports

- `coverage-baseline.md` — Phase 3 line/branch coverage snapshot. Coverage
  is tracked, not gated (see `adr/0006-coverage-tracked-not-gated.md`).

> The per-arc working docs (brainstorm specs, implementation plans, per-reader
> audits) that used to live under `superpowers/` and `archive/2026-phase-4/`
> were retired on 2026-05-31. Their substance is consolidated into the two
> capstone summaries above; the originals remain recoverable from git history.
