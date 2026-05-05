---
title: "Legacy readers physical deletion (umbrella)"
author: Mateusz Zawadzki
date: 2026-05-05
status: draft
---

# Legacy readers physical deletion — umbrella spec

Sequel to the Phase 4f-extend modernization
(`docs/PHASE-4-MODERNIZATION-SUMMARY.md`). Phase 4f-extend retired the
legacy fixed-format readers from the production runtime. ADR 0019
deferred their physical deletion to a follow-on spec because the
parity-test infrastructure (`tests/unit/io/toml/test_*_parity.pf`)
depends on `readswap('swap')` to produce ground-truth `variables`
globals for legacy-vs-TOML comparison.

This spec executes the deletion: convert the parity tests to
self-contained literal-value assertions, then physically remove the
legacy reader code from `src/`.

## Goal

Reach a state where:

1. `src/io/readswap.f90` is gone.
2. The dead `case(1)` init blocks in `src/crop/irrigation.f90`,
   `src/crop/management_soil.f90`, and the nutrient block in
   `src/crop/cropgrowth.f90` are gone (along with
   `readarablelandgerm`).
3. The parity-test infrastructure files
   (`tests/unit/io/toml/parity_helpers.f90`,
   `legacy_crop_helper.f90`, `readswap_stubs.f90`) are gone.
4. The per-suite isolation wrapper `tests/unit/run_pfunit.sh` is gone
   (it existed only because legacy-reader state leaked between parity
   suites; with no legacy reader invocations, the wrapper has no job).
5. `meson test()` invokes the pFUnit binary directly again.
6. **All 6 parity test suites still pass**, asserting against
   captured-once literal values instead of live legacy-reader output.
7. ~5500 LoC of production source is removed; the binary is smaller
   and the `src/` tree no longer contains code that doesn't run in
   production.

End-state tag: `rescue/legacy-readers-deleted`.

## Non-goals

- **Tillage and SSDI TOML port** — `Read_Tillage` and
  `SSDI_irrigation(1)` continue to read tillage / SSDI parameters
  from staged `swap.swp` via TTutil when activated. Schema-porting
  those blocks (so the modern binary doesn't need TTutil for tillage
  / SSDI either) is a future arc tracked as ADR 0021 candidate. Out
  of scope here.
- **TTutil utility-call retirement** (`rdsets`/`rdfrom` in
  `swap_main.f90`, `rddtmp` in `swapoutput.f90`) — these are not
  data readers (rerun mechanism + temp-file cleanup). Separate,
  smaller arc; out of scope.
- **TTutil utility calls** (`rdsets`, `rdfrom` in
  `swap_main.f90:47,54`; `rddtmp` in `swapoutput.f90:3454`) — these
  are not data readers. Out of scope.
- **The `archive/2026-phase-4/audits/phase-4f-config-to-variables-audit.md`
  per-section table refresh** — mechanical doc cleanup; bundle into
  SS-C closeout if convenient, otherwise defer.
- **Macropore deletion** — `macroporeflow_parity` already disabled
  (test moved to `@disable` in SS-11 closeout). No additional work
  here.

## Sub-spec roadmap

| SS | Title | Effort | Depends on |
|---|---|---|---|
| SS-A | Per-parity-suite literal-value capture + helper-file drop | M | — |
| SS-B | ADR 0020: call-site gating convention; lift tillage/SSDI gates; un-stub validators | S-M | — |
| SS-C | Code deletion (readswap.f90 + dead case(1) blocks + wrapper retirement) | M | SS-A |
| SS-D | Closeout (ADR 0019 update, architecture.md update, summary update, tag) | S | SS-A, SS-B, SS-C |

**Note on SS-B framing.** Tillage and SSDI are legitimate SWAP
features whose TOML port hasn't been written yet — not deprecated /
deferred features. The `swtill=1` / `swssdi=1` stub-errors added in
SS-10.5 sent the wrong signal ("this feature is gone"); SS-B removes
those stub-errors and instead gates the subsystems at the call site
(matching the existing convention used by `flCropNut`,
`flMacroPore`, `flSurfaceWater`, etc.). The subsystems remain
runnable via TTutil reading of staged `swap.swp` until a future arc
ports them to TOML schema. See "Out of scope" for the follow-on.

## SS-A — Per-parity-suite literal-value capture + helper-file drop

**Goal.** Each parity test currently asserts shape `legacy_global ==
config%X`. Convert each assertion to `<literal> == config%X` where
`<literal>` is the value the legacy reader would have written. After
SS-A, parity tests are self-contained: they don't invoke `readswap`,
they don't read from `legacy globals`, they don't share state with
each other.

**Approach.** For each parity test:

1. Add a temporary instrumentation pass that prints the current
   `legacy_global` value when the test runs (e.g., `write(*,*)
   'EXPECT', '<config_field>', '=', <legacy_global>`).
2. Run the test; capture stdout.
3. Replace each `@assertEqual(legacy_global, config%X)` with
   `@assertEqual(<captured_literal>, config%X)`.
4. Drop the `use variables, only: <legacy globals>` lines.
5. Drop the `call load_both_for_*` helper invocation; replace with a
   plain `call load_swap_config(...)` + `validate` + `finalize`.
6. Verify the test still passes (same literal in both branches).

**Files affected.**

| File | Change |
|---|---|
| `tests/unit/io/toml/test_hupselbrook_parity.pf` | Convert assertions; drop `parity_helpers_mod` import |
| `tests/unit/io/toml/test_grassgrowth_parity.pf` | Same |
| `tests/unit/io/toml/test_oxygenstress_parity.pf` | Same |
| `tests/unit/io/toml/test_salinitystress_parity.pf` | Same |
| `tests/unit/io/toml/test_surfacewater_parity.pf` | Same |
| `tests/unit/io/toml/test_legacy_crop_helper.pf` | This suite tests the legacy helper itself; either delete the suite or convert to TOML-only assertions. Decide during SS-A. |
| `tests/unit/io/toml/parity_helpers.f90` | DELETE (no callers after the conversion) |
| `tests/unit/io/toml/legacy_crop_helper.f90` | DELETE |
| `tests/unit/io/toml/readswap_stubs.f90` | DELETE |
| `tests/unit/meson.build` | Drop the deleted helper files from `pfunit_extra_sources` |

**Commit cadence.** One commit per parity suite (5-6 commits), then a
final helper-deletion commit. Verify pFUnit + check-full green at
every commit.

**Risk.** A parity test asserts against a `variables` field that the
TOML adapter doesn't yet populate correctly → would silently start
asserting `0 == 0` (default-vs-default). Mitigation: spot-check that
the captured literal is non-zero and matches the documented
expectation in each test's docstring.

## SS-B — Code deletion

**Prerequisite.** SS-A complete; no parity test invokes
`readswap()` or any legacy reader.

**Deletions.**

| File / region | Action |
|---|---|
| `src/io/readswap.f90` | DELETE entire file (5080 LoC) |
| `src/crop/irrigation.f90` lines 70-303 (`case(1)` block) | DELETE |
| `src/crop/management_soil.f90` lines 88-343 (`SoilManagement(1)` body) | DELETE; replace with `case (1); return` |
| `src/crop/cropgrowth.f90` lines 1151-1194 (nutrient `rdinit` block) | DELETE |
| `src/crop/cropgrowth.f90` `subroutine readarablelandgerm` and the cache-miss `else` branch in `ArableLandGerm(1)` | DELETE (both unreachable per ADR 0017 cache-hit dispatch) |
| `src/io/toml/config_to_variables.f90:1175-1196` (`logf` opening hack) | KEEP — `logf` is still consumed by `Read_Tillage` and `SSDI_irrigation(1)` for log writes. SS-10.5 already trimmed the `swpfile` line. |
| `tests/unit/run_pfunit.sh` | DELETE — its only purpose was per-suite isolation from legacy reader state leakage; with no readers running, no leakage |
| `tests/unit/meson.build` `test()` invocation | Restore `test('unit-swap-tests', unit_tests, …)` directly (drop the wrapper invocation) |
| `meson.build` (top-level production sources) | Drop `'src/io/readswap.f90'` |
| `tests/unit/meson.build` `pfunit_extra_sources` | Drop `'../../src/io/readswap.f90'` |

**Verification at each step.**

- After each file/block deletion: `pixi run -e test build-linux` → no
  unresolved symbols.
- After all deletions: `pixi run -e test test-pfunit` → 548 passed
  (or close — minus any tests that were specifically testing legacy
  reader behaviour and got removed in SS-A).
- After all deletions: `pixi run -e test check-full` → 5/5 green.

**Commit cadence.** One commit per deletion target (file or block);
~6-8 commits total. Each independently bisectable.

## SS-C — Closeout

**Goal.** Capture the new state in docs and tag.

**Updates.**

- `docs/adr/0019-legacy-readers-retired.md` — append "Update
  2026-MM-DD: physical deletion executed" section. The "retirement
  gate" forward-look becomes a closing note.
- `docs/architecture.md` — I/O layer section: drop the "remain
  compiled but unreachable" paragraph; describe `src/` as
  modernization-only.
- `docs/PHASE-4-MODERNIZATION-SUMMARY.md` — add a final-state row to
  the "What was achieved" table; remove "physical deletion of
  readswap.f90" from "What's next".
- `docs/contributing.md` — drop any remaining mentions of legacy
  readers being in tree.

**Tag.** `rescue/legacy-readers-deleted` placed at the closing
commit.

**Effort.** S (mostly mechanical doc edits + tag).

## Acceptance criteria (umbrella)

- [ ] All 6 parity test suites pass against literal-value assertions
  (no `readswap` invocation in any test).
- [ ] `src/io/readswap.f90` does not exist.
- [ ] `grep -rn "subroutine readswap\b" src/` returns no matches.
- [ ] `grep -rn "call irrigation(1)\|call SoilManagement(1)" src/`
  returns no matches.
- [ ] `tests/unit/run_pfunit.sh` does not exist.
- [ ] `pixi run -e test test-pfunit` exits 0 with `Ok: 1, Fail: 0`.
- [ ] `pixi run -e test check-full` exits 0 with 5 passed.
- [ ] ADR 0019, `architecture.md`, and
  `PHASE-4-MODERNIZATION-SUMMARY.md` reflect the deletion.
- [ ] Tag `rescue/legacy-readers-deleted` placed.
