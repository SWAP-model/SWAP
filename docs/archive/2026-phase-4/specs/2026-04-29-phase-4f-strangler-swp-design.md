---
title: "Phase 4f — Strangler-fig: disconnect readswap, wire TOML into swap.f90"
author: Mateusz Zawadzki
date: 2026-04-29
status: draft
---

# Phase 4f: Strangler-fig — disconnect readswap, wire TOML into swap.f90

The big leap. After this phase, `src/core/swap.f90` no longer calls `readswap()` — it loads `swap.toml` via the new TOML pipeline and copies fields to legacy `variables%` globals via a new `config_to_variables` adapter. The five non-macropore regression cases run end-to-end through the new path and produce output matching their fixtures. Macropore (case 3) is excluded from `check-full` per ADR 0011.

## Background

Phases 4a–4e built the typed-config + TOML reader pipeline; Phase 4f-prep closed enough schema gaps that the strangler-fig is feasible. The Phase 4e variables audit identified 236 (C) fields across 11 typed configs; the adapter copies those fields to `variables%` globals so legacy physics code keeps reading what it expects, but from a TOML-driven source rather than the legacy `readswap()`.

Macropore physics (case 3) was deferred per ADR 0010 — its 22-key schema lives at `src/config/macropore_config.f90` as orphan infrastructure. With case 3 now excluded from `check-full`, Phase 4f doesn't need a runtime fallback to the legacy reader: cases 1, 2, 4, 5, 6 take the new path; case 3 simply isn't run in regression.

`readswap()` itself **stays intact**. Phase 4f only stops *calling* it from `swap.f90` — the file remains in the tree as the reference for Phase 4f-extend, which will walk it line by line and add schema slots for every legacy input key not already covered (the 5 test cases exercise only a subset).

## Scope

### In

- **`src/io/toml/config_to_variables.f90`** — new module with one public subroutine `config_to_variables(config)`. Walks the typed config in audit order (general → simulation → meteorology → drainage → soil → bottom_boundary → heat → irrigation → solute → surface_water → crop), assigning each (C)-classified field to its `variables%` global. Plus zero-forcing of the 18 RETIRED switches per ADR 0009. ~250–350 LoC.
- **`src/core/swap.f90`** — replace `call ReadSwap()` in the Initialize block with the TOML pipeline:
  ```fortran
  call load_swap_config('swap.toml', config, errors)
  call config%validate(errors)
  call config%finalize(errors)
  call errors%abort_if_fatal()
  call config_to_variables(config)
  ```
  Binary always reads `swap.toml` from cwd; no CLI arg support in this phase.
- **`tests/swap-cases/run_case.sh`** — add `--toml` flag. When set, copy `tests/swap-cases/toml/<case>/swap.toml` (plus cross-file siblings `swap.dra.toml` and `*.crp.toml`) into the case dir. Without `--toml`, script behaves as today (legacy `.swp` mode).
- **`tests/regression/test_output_regression.py`** — drop `macroporeflow` from the case iteration list.
- **ADR 0011** — "Macropore case exclusion from regression suite" — supersedes the relevant clause of ADR 0010 (which had said "case 3 stays in regression via legacy fallback").
- **Per-test-case TOML completeness** — verify each of cases 1, 2, 4, 5, 6 has a complete `swap.toml` (+ cross-file siblings) that loads, validates, and finalizes without errors. Iteratively fix gaps as the binary surfaces them.
- **Tag** `rescue/phase-4f-strangler-swp` once 5/5 regression green via the new path.

### Out

- **`readswap.f90` deletion or trimming** — file stays intact as Phase 4f-extend's reference.
- **Phase 4g** — disconnect crop sub-readers (`readwofost`, `readcropfixed`, `readgrass`) from `cropgrowth.f90`. Separate phase.
- **Phase 4f-extend** — walk `readswap()` to add schema slots for ALL legacy input keys (not just test-case subset). Separate phase, opens after Phase 4f.
- **Macropore reactivation** — future macropore physics phase.
- **CLI argument for project name** — defer (binary always opens `swap.toml` in cwd).
- **Cleanup of `--toml` flag** — keep both modes through Phase 4g; cleanup post-readswap-deletion.
- **Coverage rebaseline / schema doc bump beyond what's already there** — Phase 4f doesn't add new schema; it just connects the existing one.

## Critical decisions

### D1 — `readswap.f90` stays intact, untouched (per user direction)

Phase 4f does NOT delete `readswap()`. The strangler-fig disconnects it from `swap.f90`'s execution path, but the subroutine and its sub-readers (`readwofost`, `readcropfixed`, `readgrass`) keep compiling. Two reasons:

1. `legacy_crop_helper.f90` (test-only) calls the sub-readers from parity tests. Keeping them callable preserves Phase 4c-b's validation infrastructure.
2. Phase 4f-extend will walk `readswap()` line by line. Having it in the tree as a reference is essential.

Once Phase 4f-extend is done and the schema covers every legacy key, a later phase deletes `readswap()` and orphans the sub-readers (or moves them to a focused `legacy_crop_readers.f90`).

### D2 — `config_to_variables` shape: single flat subroutine

One file, one public subroutine, no per-domain helpers. Body is a long sequence of `variables%foo = config%<section>%foo` lines plus allocatable-array copies. Reasons:

- Mirrors the audit doc's layout 1:1, which makes review trivial.
- Single point to grep when a global isn't being populated (failure mode is "I added a field to the schema but forgot to wire it through" — easy to spot in one file).
- No abstraction tax; the file is mechanical, not algorithmic.

If the file grows past ~500 lines or readability degrades, Phase 4f-extend can split it. For now, flat.

### D3 — Failure mode for missing `swap.toml`

Binary at startup calls `load_swap_config('swap.toml', config, errors)`. If the file doesn't exist, `load_swap_config` appends `ERR_IO_READ_FAILED` (fatal); `abort_if_fatal()` prints summary and `error stop`s. Clear failure message: "Could not open swap.toml in current directory."

No fallback to legacy `swap.swp`. The whole point of Phase 4f is to make the new path the only path.

### D4 — Macropore case (3) exclusion from `check-full`

Drop the case from `tests/regression/test_output_regression.py`'s case list. Case dir stays on disk for archival; `swap_linux.swp.template`, the legacy `.crp` files, and the partial TOML files (which omit the `[macropore]` block per ADR 0010) all stay. The case is preserved but not exercised in CI / regression.

ADR 0011 documents the decision and supersedes the `swap.f90` runtime-branch clause from ADR 0010 (which said "case 3 stays in regression via legacy fallback" — it doesn't anymore).

### D5 — `run_case.sh` `--toml` flag semantics

```bash
./run_case.sh --case hupselbrook --toml --exec ../../builddir/swap
```
1. `cd` into `tests/swap-cases/1.hupselbrook/`.
2. Copy `tests/swap-cases/toml/1.hupselbrook/swap.toml` → cwd.
3. Copy `tests/swap-cases/toml/1.hupselbrook/swap.dra.toml` → cwd (if exists).
4. Copy `tests/swap-cases/toml/1.hupselbrook/*.crp.toml` → cwd (glob).
5. Run binary.
6. Cleanup: `rm -f swap.toml swap.dra.toml *.crp.toml swap_swap.log Swap.ok *rd$*.tmp` plus the existing cleanup.

Without `--toml`: legacy mode (no change to today's behavior).

### D6 — Iteration order

Per user guidance: get one case working at a time. Order:
1. **hupselbrook** (case 1) — has all 3 crop types (1=fixed maizes, 2=WOFOST potatod, 3=grass grassd); good coverage of the typed-config crop dispatch.
2. **grassgrowth** (case 2) — type-3 grass, dramet=3 multi-level drainage.
3. **oxygenstress** (case 4) — similar shape to grassgrowth, exercises oxygen-stress crop section.
4. **salinitystress** (case 5) — type-2 WOFOST with `swsalinity=1`.
5. **surfacewater** (case 6) — `swsrf=2` with full surface-water management; the most complex case.

Each case is one or more focused commits: typically (a) attempt to run; (b) fix the assignment(s) in `config_to_variables` that surface; (c) re-run; repeat until output matches fixture.

After all 5 individually pass, run `pixi run -e test check-full` (with macropore already excluded). Should be 5/5 green.

## Critical files

### Create
- `src/io/toml/config_to_variables.f90` — the adapter.
- `tests/unit/io/toml/test_config_to_variables.pf` — adapter unit tests (assignment correctness).
- `docs/adr/0011-macropore-exclusion-from-regression.md` — supersedes ADR 0010's runtime-fallback clause.
- `tests/regression/baselines/phase-4f-strangler-swp.log` — closeout regression log.

### Modify
- `src/core/swap.f90` — Initialize block: replace `call ReadSwap()` with TOML pipeline + adapter.
- `tests/swap-cases/run_case.sh` — add `--toml` flag; matching cleanup.
- `tests/regression/test_output_regression.py` — drop macroporeflow.
- `meson.build` — register `config_to_variables.f90` in `sources`.
- `tests/unit/meson.build`, `tests/unit/testSuites.inc` — register the new pFUnit suite.
- `docs/adr/0010-macropore-deferral.md` — note that ADR 0011 supersedes the runtime-fallback clause.
- `docs/adr/index.md` — index ADR 0011.
- `docs/configuration-schema.md` — add a one-line note that the binary now reads `swap.toml` directly.

### Delete
- None.

## Verification gates

After each iteration step:
- `pixi run -e test build-linux` — binary builds.
- `pixi run swap --case <name> -- --toml` — case runs to completion (no segfault, no FOPENG abort).
- `pixi run -e test regression <name>` — output matches fixture.

After all 5 cases pass individually:
- `pixi run -e test check-full` — 5/5 green (macropore excluded).
- F count from dot stream stays at 0 (pFUnit suite still green).
- Submit as the canonical phase-end gate.

After tag:
- `git diff rescue/phase-4f-prep-gap-closure..HEAD -- src/` — should show ONLY `swap.f90` (one block change), `config_to_variables.f90` (new), and `meson.build` registrations. No other source changes. (Adapter authoring may surface bugs in the typed configs / readers / finalizers — those get separate focused commits as they appear.)

## Risks & mitigations

| Risk | Mitigation |
|---|---|
| Adapter omits a `variables%` field that physics depends on; segfault or wrong output | Iterate one case at a time. The audit doc enumerates every C-classified field; the adapter is mechanical. Unit tests for the adapter check assignment per section. |
| TOML schema is missing a key that case N uses (gap not caught by Phase 4f-prep audit) | Iterate surfaces it. Fix surfaces as schema-extension commits in Phase 4f. The 76 residual real-domain G entries in the audit may include candidates; Phase 4f-extend addresses systematically. |
| Validator/finalize fails on a real regression-case TOML that worked at parity-test time | Loosen validators only when the legacy reader has the same permissive behavior. Document each loosening with a `! READING NOTE`. |
| Cross-file TOML reader (`drainage.file = "..."`, `[[crop.rotation]].file = "..."`) doesn't resolve paths correctly when running from case dir | Path resolution already works in parity tests; the run_case.sh `cd` puts cwd in the case dir before the binary launches, so relative paths resolve as expected. Smoke-tested by case 1's first run. |
| `bbcfil` external file references stop working | Phase 4d's `read_bottom_boundary_toml` already handles the bbc reference; preserve that handling. If a case's bottom boundary breaks, surface in the iteration. |
| Output writer code paths still depend on RETIRED legacy switches being non-zero | Adapter explicitly forces them to 0 per ADR 0009. The check-full output comparisons use only the CSV stream which is governed by `swcsv`/`swcsv_tz` (kept). |

## File count summary

- Create: ~4 files (adapter, adapter tests, ADR 0011, baseline log).
- Modify: ~5 files (swap.f90, run_case.sh, test_output_regression.py, meson configs, schema doc note).
- Per-iteration fix commits: variable count (0 to ~10 depending on what surfaces).
- LoC delta: ~+400 / ~-15 (adapter is ~300 LoC; swap.f90 net ~+5 LoC; run_case.sh ~+25 LoC).
