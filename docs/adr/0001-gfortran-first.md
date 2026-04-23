# ADR 0001 — gfortran-first (rescue phase)

Status: accepted (2026-04-22, during rescue Phase 1)

## Context

At the rescue baseline commit `e256bc0`, `meson.build` branched between Intel Fortran (`ifx` / `ifort`) and GCC (`gfortran`), with a separate prebuilt-binary download path for ttutil under Intel and a source-subproject path under GCC. Separately, `pixi.toml` hardcoded `FC=ifx` in the production configure tasks, so the Phase 0 baseline was inadvertently recorded under Intel even though a later survey had assumed gfortran-only. Two supported compilers means two sets of flags, two build behaviours, two failure modes to investigate when a regression surfaces.

The rescue's explicit ordering is: *working car first, then fine-tuned bolid*. Reducing compiler variance removes one axis of uncertainty while physics and test-coverage debt are being paid down.

## Decision

During the rescue (Phases 1–4), the only officially supported Fortran compiler is **GCC / gfortran** as provided by pixi. `meson.build` rejects any other compiler at configure time with an explicit error message pointing at this ADR. `pixi.toml` configure tasks set `FC=gfortran` explicitly.

The gfortran flag set includes `-finit-local-zero` to match Intel's `-init=zero` behaviour. SWAP's legacy Fortran code relies on zero-initialisation of local variables; without this flag, `salinitystress` produces NaN values and `macropore` drifts substantially.

## Consequences

- **Positive**: deterministic builds. One flag set. ttutil always resolves via the meson subproject. No Intel-specific download-retry / corporate-proxy path.
- **Positive**: simpler meson.build — the Intel branch (~30 lines including the ttutil download logic) is gone.
- **Positive**: simpler pixi.toml — no `oneapi setvars` wrappers around build tasks.
- **Negative**: the Intel binary has historically been faster at some physics hot paths. A performance regression follow-on spec (between Phase 4 and the compartment refactor) will re-evaluate.
- **Negative**: users who habitually source Intel oneAPI in their shell will see a clear configure error. The error message points at this ADR.
- **Negative**: the swap from ifx to gfortran produced non-trivial numeric divergence in `macropore` (~21-unit DRAINAGE drift at year 1998). Rather than pick a side, the rescue keeps both reference sets — historical ifx-produced fixtures at `tests/regression/*_expected.json`, and current gfortran+finit-local-zero fixtures at `tests/regression/*_expected_gfortran.json`. The harness compares against the latter. See `tests/regression/INVESTIGATION_NOTES.md` for open questions.

## Re-enablement

Intel/ifx support may return after Phase 4 exit, once the codebase is stabilised. The criteria for re-enablement:
1. All pre-compartment gates from the rescue spec satisfied.
2. A dedicated follow-on spec (Phase 4.x or separate) scoped to restore Intel with regression green on both compilers.
3. The regression harness extended to exercise both compilers and select between `_expected.json` / `_expected_gfortran.json` fixtures per compiler.

Reverting this ADR is a deliberate, test-gated act — not an environment-variable flip.
