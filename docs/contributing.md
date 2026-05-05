---
title: Contributing
author: SWAP modernization team
---

# Contributing

## You are here

This repository was rescued through Phase 4 / Phase 4f-extend (completed
2026-05-05); see `docs/PHASE-4-MODERNIZATION-SUMMARY.md` for the
capstone overview and `docs/archive/2026-phase-4/` for the frozen
specs, plans, and audits. Trivial typo fixes and documentation edits
do not need the full context; anything touching source, build, tests,
or fixtures does.

## Branch model during the rescue

The rescue uses a small set of long-lived branches rather than the usual
per-change feature-branch flow:

- `main` advances only at phase completion, marked by a
  `rescue/phase-N-<shortname>` tag. It does not accept direct commits
  mid-phase.
- `development` is where all in-phase commits land. Phases 1–3 (including
  Phase 2, the current phase) work directly on `development`.
- `archive/*` branches preserve the pre-rescue state and must not be lost:
  `archive/main-pre-rescue`, `archive/swaplib`, `archive/swaplib-simple`,
  and `archive/wip-drifted`.
- `legacy/swap-4.2.0` is an orphan branch carrying the upstream
  Intel-compiled SWAP 4.2.0 tree for reference.
- Starting in Phase 4, per-change feature branches branch off
  `development` and merge back by fast-forward once green.
- `origin/main` is **not** advanced during the rescue. Nothing is pushed
  to it until the rescue completes.

## Commit conventions

Commits follow a conventional-commits prefix scheme:

- `feat:` — new functionality.
- `fix:` — bug fix.
- `docs:` — documentation only.
- `refactor:` — code reshaping without behavior change.
- `chore:` — infrastructure or tooling.
- `test:` — test additions or modifications.
- `style:` — formatting or lint config.

An optional scope goes in parentheses: `chore(build):`, `fix(tests):`,
`docs(adr):`.

Subject line: imperative mood, present tense, no trailing period. Keep it
under 72 characters. The body, if present, is wrapped at 72 columns and
separated from the subject by a blank line. Rescue commits often carry a
substantial body explaining rationale — err on the side of more context,
especially when the change is non-obvious from the diff alone.

## Verification gates

Every commit on `development` must pass the fast gate:

    pixi run -e test check-fast

This is the everyday iteration gate. It runs in roughly 90 seconds and
covers build, pFUnit unit tests (currently empty at the rescue baseline),
and the four fast regression cases.

Every phase-end tag must additionally pass the full gate:

    pixi run -e test check-full

`check-full` extends `check-fast` with the two slow regression cases
(`oxygenstress` and `macropore`). It runs in roughly 10 minutes. This is
the gate for fast-forwarding `main` to `development` at phase boundaries.

When you write new Fortran source, run `pixi run lint-check` against the
files you changed to inspect what the formatter would change. New code
is expected to be lint-clean at commit time; legacy code is not. See
`docs/code-style.md` Section 8 for the advisory-only rescue lint policy.

Tag names follow the convention `rescue/phase-N-<shortname>`:
`rescue/phase-0-baseline`, `rescue/phase-1-infra`, `rescue/phase-2-docs`,
and so on.

## Fixture policy

The regression harness compares against
`tests/regression/*_expected_gfortran.json`. The companion
`*_expected.json` files are historical ifx-compiled references —
preserved in git as a diff target but never compared against at runtime.
Do not delete them.

If a deliberate physics change legitimately moves the expected output,
regenerate fixtures:

    python tests/regression/test_output_regression.py --regenerate-fixtures

Document the regeneration in the commit message. Include at least:

1. What physics change caused the output to move.
2. Before and after numeric values for at least one affected variable
   (e.g., "macropore DRAINAGE 1998: 114.92 → 115.30").
3. A rationale for why the new values are correct.

**Never regenerate silently to make a failing test pass without
understanding why.** The fixture is the truth about what the code should
produce; if the code produces something else, that is a physics
question, not a fixture question.

Known accepted tolerances already baked into the `_gfortran` fixtures:

- `macropore` DRAINAGE: ~21-unit divergence at year 1998 versus the ifx
  reference.
- `oxygenstress` MOWDM: ~85-unit deviation at year 1995 (pre-existing
  under both compilers).

See `tests/regression/INVESTIGATION_NOTES.md` for the open
compiler-vs-physics questions on these.

## Style

Run `pixi run lint-check` on the files you changed. New code must be
lint-clean; legacy code that happens to sit in the same commit is not
expected to clean up as a side effect. `.fprettify.rc` at the repo root
is the config source of truth.

Indent and whitespace conventions that apply independent of the
formatter: spaces not tabs, lowercase keywords, snake_case identifiers.
See `docs/code-style.md` for the full style reference.

## Adding new code

- Adding a new domain state type: see `docs/state-management.md`
  Section 7 ("How to add a new domain state").
- Adding a new TOML configuration key or section: see
  `docs/configuration-schema.md`.
- Adding a new test suite: see `docs/build-and-test.md`.

When a change touches multiple concerns, prefer multiple small commits
(one per concern) over one large commit. The conventional-commits prefix
tells readers what kind of change each commit contains — keeping commits
focused lets readers skip categories they do not care about.

## Running agents

The `.github/agents/swap-fortran.agent.md` file specifies the context an
automated agent should load before making changes: the rescue spec, the
baseline record, the current phase plan, the ADR directory, and the
build/test guide. If you invoke an agent that does not follow these
conventions (for example, an agent that starts modifying code without
reading the rescue spec), stop it and point it at the agent definition
file.

## When in doubt

In order: look at adjacent similar code, read the relevant ADR, check
the relevant spec or plan under `docs/archive/2026-phase-4/` (or for
new modernization arcs, the active spec/plan), and open a
discussion (GitHub issue or direct conversation) before guessing in a
commit. The rescue has enough moving pieces that undocumented decisions
are expensive — explicit disagreement is cheaper than silent drift.
