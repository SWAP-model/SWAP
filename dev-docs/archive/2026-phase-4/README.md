---
title: Archive — Phase 4 modernization (2026-04-22 to 2026-05-05)
---

# Phase 4 modernization archive

Frozen historical record of the Phase 4 / Phase 4f-extend modernization
work that ran from 2026-04-22 (rescue baseline) to 2026-05-05 (legacy
reader retirement). The capstone summary is at
[`docs/PHASE-4-MODERNIZATION-SUMMARY.md`](../../PHASE-4-MODERNIZATION-SUMMARY.md).

## What's in here

```
2026-phase-4/
├── README.md                ← you are here
├── specs/                   ← 18 design specs (one per phase + sub-phases)
├── plans/                   ← 22 implementation plans (matched to specs)
└── audits/                  ← 23 per-reader / per-case audit docs
```

## How to navigate

- **Want to see why a code path is the way it is?** Look first at the
  matching ADR under `docs/adr/`. ADRs were intentionally kept in place
  because they remain authoritative decision records, even when the
  decisions are now frozen-in-code.
- **Want the full chronology?** See the capstone summary
  ([`PHASE-4-MODERNIZATION-SUMMARY.md`](../../PHASE-4-MODERNIZATION-SUMMARY.md))
  — it lists every sub-spec with its closing commit and audit pointer.
- **Want a specific phase's design?** Search `specs/` by date prefix
  (chronological).
- **Want a specific phase's tasks?** The matching plan in `plans/` has
  the same date prefix.
- **Want to see how a particular legacy reader was retired?** Look in
  `audits/` for `phase-4f-<reader>-audit.md`.

## Tags marking the boundaries

| Tag | What it marks |
|---|---|
| `rescue/phase-0-baseline` | Pre-modernization snapshot |
| `rescue/phase-1-infra` | Build + test infrastructure live |
| `rescue/phase-2-docs` | Architecture + ADR backbone written |
| `rescue/phase-3-coverage` | Coverage baseline established |
| `rescue/phase-4a-infrastructure` | TOML config skeleton |
| `rescue/phase-4e-error-prep` | Error-collection pipeline live |
| `rescue/phase-4f-prep-gap-closure` | Strangler swap prerequisite gates |
| `rescue/phase-csv-meteo-complete` | CSV meteorology done |
| `rescue/phase-swap-ini-port` | swap.ini → TOML done |
| `rescue/phase-4f-extend-complete` | **Legacy reader retirement DONE** |
| `rescue/phase-4f-extend-followups` | SS-10.5 + SS-5 M1-M5 polish |

## Why archive instead of delete

Two reasons:

1. **Audit trail.** Each spec / plan / audit captures *why* a decision
   was made and what alternatives were considered. The git log alone
   doesn't carry that depth. Future debate ("why did we DEFER this?",
   "what was the original spec?") is settled by re-reading the
   contemporaneous doc.
2. **Pattern reuse.** The Phase 4 sub-spec staging
   (audit → schema → adapter → reader-deletion → close-doc) became the
   modernization template. Future modernization arcs (compartment-state
   refactor, performance work, ifx port, Python bindings) can re-use
   the same shape.

## What's NOT here (still active in main `docs/`)

- All 19 ADRs (decision records — permanent)
- `architecture.md` — current architecture
- `configuration-schema.md` — current TOML schema reference
- `state-management.md`, `error-handling.md`, `logging.md`,
  `validation.md` — current domain reference
- `csv-companion-files.md`, `meteorology.md`,
  `toml-format-guide.md` — input-format reference
- `build-and-test.md`, `code-style.md`, `contributing.md`,
  `dependency-management.md` — process docs
- `coverage-baseline.md` — current coverage status
- `PHASE-4-MODERNIZATION-SUMMARY.md` — capstone summary
