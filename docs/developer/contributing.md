---
title: Contributing
author: SWAP modernization team
---

# Contributing

## You are here

This repository was rescued through Phase 4 / Phase 4f-extend (completed
2026-05-05) and then modernized further; see
`dev-docs/phase-4-modernization-summary.md` and
`dev-docs/post-phase-4-modernization-summary.md` for the capstone
overviews, and `dev-docs/adr/` for the architecture decision records.
Trivial typo fixes and documentation edits
do not need the full context; anything touching source, build, tests,
or fixtures does.

## Branch model

| Branch | Role | Push policy |
|---|---|---|
| `main` | Modernization release line. Each merge into `main` corresponds to a tagged release. | **No direct commits, no untagged merges.** Merge from `development` only when shipping a release. |
| `development` | Active day-to-day work. Topic branches feed into it. | Push freely; CI runs on tag pushes only (see below). |
| topic branches (optional) | Short-lived feature work. | Merge into `development` once green. |
| `v4.2.0` tag | Legacy SWAP 4.2.0 reference snapshot. Downstream consumers (pyswap) pin to this. | Frozen. |

### Release flow

```bash
# On development (or a topic branch merged into it), confirm gates green:
pixi run -e test test-pfunit       # → Ok: 1, Fail: 0
pixi run -e test check-full        # → 5 passed, 0 failed

# Merge into main and tag the release. The tag triggers GitHub Actions
# (build + test + release artifacts).
git checkout main
git merge --no-ff development -m "Release vX.Y.Z"
git tag -a vX.Y.Z -m "<release notes>"
git push origin main vX.Y.Z
```

**Never push to `main` without a corresponding tag.** Untagged pushes
to `main` will not trigger CI under the current workflow trigger
(`tags: ['v*']`), so an artifactless push silently weakens the
release line.

### Other branches

- `archive/*` branches preserve pre-rescue state and must not be
  deleted: `archive/main-pre-rescue`, `archive/swaplib`,
  `archive/swaplib-simple`, `archive/wip-drifted`.
- `legacy/swap-4.2.0` is an orphan branch carrying the upstream
  Intel-compiled SWAP 4.2.0 tree for reference.

See [`branches.md`](branches.html) for the full convention narrative.

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

`CLAUDE.md` (repo root) is the operating contract for any automated agent:
non-negotiables, build/test commands, the current architecture, commit
conventions, and where to look (ADRs, the modernization summaries, and the
deep `.github/DEVELOPMENT_GUIDE.md`). If you invoke an agent that ignores
these conventions — for example, one that starts changing code without
running `check-fast` or that edits physics formulas — stop it and point it
at `CLAUDE.md`.

## When in doubt

In order: look at adjacent similar code, read the relevant ADR, check
the capstone modernization summaries under `dev-docs/`, and open a
discussion (GitHub issue or direct conversation) before guessing in a
commit. The rescue has enough moving pieces that undocumented decisions
are expensive — explicit disagreement is cheaper than silent drift.
