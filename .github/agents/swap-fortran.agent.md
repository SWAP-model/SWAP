---
name: swap-fortran
description: Fortran agent for the SWAP modernization repo — operates under the rescue spec.
tools: [read, edit, search, bash, todo]
---

# swap-fortran agent

Use this agent to work on the SWAP modernization. It operates under the rescue-and-stabilize workflow; its decisions must remain consistent with that workflow until the workflow completes.

## Required context before acting

Read these before making any change:

1. `docs/superpowers/specs/2026-04-22-rescue-and-stabilize-design.md` — the rescue spec.
2. `docs/superpowers/specs/2026-04-22-baseline-record.md` — what was observed at the green baseline `e256bc0`.
3. The current Phase's plan under `docs/superpowers/plans/`.
4. `docs/adr/` — every ADR (architecture decision record). They record the non-obvious choices.
5. `tests/regression/INVESTIGATION_NOTES.md` — open questions about macropore / oxygenstress fixture divergences.
6. `docs/pfunit-vendoring.md` — the pFUnit gitlink peculiarity and how to rebuild the install.

## Non-negotiables

- **Compiler**: gfortran only. Do not re-introduce Intel/ifx code paths. See `docs/adr/0001-gfortran-first.md`.
- **Build**: a single `builddir/` with `enable_pfunit=true` by default. See `docs/adr/0002-single-builddir.md`.
- **Verification**: every commit must pass `pixi run -e test check-fast`. Phase tags also require `pixi run -e test check-full` to be green at documented tolerances.
- **Accepted tolerances**: the `oxygenstress` MOWDM deviation and the `macropore` DRAINAGE divergence are captured in `tests/regression/*_expected_gfortran.json` fixtures — they are baked in as the current truth, NOT to be re-derived case by case.
- **Branch model**: `main` only moves forward on phase completion; `development` is where in-phase commits land; Phase 4 onward uses per-change feature branches. Nothing is pushed to `origin/main` during the rescue.
- **Fixture policy**: do not delete the `*_expected.json` files (historical ifx reference). The harness compares against `*_expected_gfortran.json`. If a physics change legitimately updates expected output, use `python tests/regression/test_output_regression.py --regenerate-fixtures` and document the regeneration in the commit message.

## When in doubt

Stop and ask. A partial, honest report beats a confident-but-wrong commit.
