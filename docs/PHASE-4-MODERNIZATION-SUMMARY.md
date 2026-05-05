---
title: Phase 4 modernization — capstone summary
date: 2026-05-05
status: complete
tags: [rescue/phase-4f-extend-complete, rescue/phase-4f-extend-followups]
---

# Phase 4 modernization — capstone summary

**Window:** 2026-04-22 → 2026-05-05 (15 days)
**Outcome:** SWAP's runtime input pipeline modernized from a TTutil-based
fixed-format reader chain (.swp / .dra / .crp / .met / .YYY) to a single
typed TOML pipeline. All six regression cases pass against the modern
binary running TOML-only inputs.

## What was achieved

| Goal | State at start | State at end |
|---|---|---|
| **Build & test infrastructure** | Out-of-date drift; tests not running | meson + pixi unified; pFUnit harness with per-suite isolation; check-full regression suite green |
| **Input format** | TTutil-based fixed-format `.swp` / `.dra` / `.crp` / `.met` / `.YYY` | Typed TOML pipeline (`load_swap_config` → `swap_config_t` → `config_to_variables` → legacy globals) |
| **Crop config** | Three legacy readers (`readcropfixed` / `readwofost` / `readgrass`) read from per-rotation `.crp` | Per-rotation TOML cache (ADR 0016) + per-crop init modules under `src/crop/{cropfixed,cropwofost,cropgrass}_init.f90` |
| **Meteorology** | TTutil per-year `.YYY` and all-years `.met` | CSV via `[meteorology.temporal].file` + optional `detail_file`/`events_file` (ADR 0014) |
| **Output** | 18 toggleable output formats (`swafo`, `swaun`, `swvap`, …) | CSV-only via `[output.csv]` (ADR 0009) |
| **Error handling** | TTutil `FatalERR` → bare `STOP` (exit 0) — silent test failures | `error_collection_t` accumulates errors; `error stop "..."` (exit 1) at one abort checkpoint (ADR 0008, 0018) |
| **Macropore physics** | Legacy reader populated globals; case 3 ran via legacy path | Stub-errored at validation; case 3 excluded from regression (ADR 0010, 0011) |
| **Legacy reader runtime calls** | Whole-stack: `readswap()` → `irrigation(1)` → `SoilManagement(*)` → `cropgrowth.f90` nutrient block, all reading from `.swp` / `.crp` | Zero TTutil **data-reader** calls reachable from the production main loop. Source code retained as parity-test fixtures (ADR 0019). |

## How it was done — the staging template

The work used a repeatable **audit → schema → adapter → reader-deletion
→ close-doc** loop, applied per legacy reader / per phase. The shape:

1. **Audit doc.** Enumerate every legacy key the reader touches, where
   each value lands in `variables` module globals, and which TOML
   schema field (existing or new) should mirror it.
2. **Schema additions.** Add typed fields to the relevant `*_config_t`
   in `src/config/`. Validators stub-error switches no case exercises
   (ADR 0015 — narrow scope, defer broadly).
3. **Reader.** Wire `read_<section>_toml.f90` to populate the schema.
4. **Adapter.** Extend `config_to_variables.f90` (or per-crop init
   modules from ADR 0016) to copy schema → legacy globals at the
   correct point (init time / per-rotation activation).
5. **Reader retirement.** Remove the runtime call to the legacy reader
   from the production call graph, or stub-error the runtime gate
   (e.g. `flCropNut`, `swmacro`).
6. **Close-doc.** Mark the sub-spec DONE in the umbrella roadmap with
   pointers to the closing commit + audit.

Verification at each step: `pixi run -e test test-pfunit` (548 tests
green) + `pixi run -e test check-full` (5 of 6 regression cases pass;
case 3 = macropore = excluded per ADR 0011).

## Phase chronology

The work split into **rescue phases** (Phase 1-3, infrastructure) and
**modernization phases** (Phase 4 with sub-phases 4a-4f-extend).

| Phase | Date | What landed | Tag |
|---|---|---|---|
| 0 | 2026-04-22 | Baseline snapshot | `rescue/phase-0-baseline` |
| 1 | 2026-04-22 | meson + pixi + CI infrastructure | `rescue/phase-1-infra` |
| 2 | 2026-04-23 | Architecture + ADR backbone (`docs/architecture.md`, ADRs 0001-0006) | `rescue/phase-2-docs` |
| 3 | 2026-04-24 | Coverage baseline established | `rescue/phase-3-coverage` |
| 4a | 2026-04-24 | TOML config skeleton (`load_swap_config`, `swap_config_t`, ADR 0007) | `rescue/phase-4a-infrastructure` |
| 4b | 2026-04-25 | Per-case field coverage; gap inventory | (no tag) |
| 4c-a/b | 2026-04-25/26 | Crop config readers (cropfixed/cropgrass + cropwofost) | (no tag) |
| 4d | 2026-04-27 | Remaining configs (drainage entries, bottom_boundary, heat, irrigation, output_csv, solute, surface_water) | (no tag) |
| 4e | 2026-04-27 | Error-collection pipeline live (ADR 0008) | `rescue/phase-4e-error-prep` |
| 4f-prep | 2026-04-28 | Strangler swap prerequisites (gap closure) | `rescue/phase-4f-prep-gap-closure` |
| 4f-strangler | 2026-04-29 | `readswap()` removed from `src/core/swap.f90` call graph | (no tag) |
| 4f .swp.ini/.dra port | 2026-04-30 / 05-01 | swap.ini → `[soil.initial]`; swap.dra → `surfacewater_init` | `rescue/phase-swap-ini-port` |
| 4f CSV meteo | 2026-04-30 / 05-01 | CSV daily + sub-daily detail; ADR 0012, 0013 | `rescue/phase-csv-meteo-complete` |
| 4f .crp port (Phases 1-4) | 2026-05-02 | Per-rotation crop cache (ADR 0016, 0017); cropfixed → cropwofost → cropgrass; hupselbrook full | (no tag) |
| **4f-extend** | **2026-05-04 → 05-05** | **Legacy reader retirement umbrella (SS-1..SS-11); ADR 0014, 0015, 0018, 0019** | **`rescue/phase-4f-extend-complete`** |
| Follow-ups | 2026-05-05 | SS-10.5 (swtill/swssdi schema port) + SS-5 M1-M5 polish | `rescue/phase-4f-extend-followups` |

## Why these decisions — ADR landing zone

ADRs 0001-0006 are infrastructure decisions (compiler, build dir,
state design, test framework, licensing, coverage). ADRs 0007-0019
shape the modernized runtime:

| ADR | Decision | Why |
|---|---|---|
| [0007](adr/0007-config-validate-finalize-pipeline.md) | Four-stage pipeline: parse → validate → finalize → adapter | Single canonical input flow with each stage independently testable |
| [0008](adr/0008-error-collection-over-fatalerr.md) | `error_collection_t` over scattered `fatalerr` calls | Full-pass validation; users see all errors at once, not one-at-a-time |
| [0009](adr/0009-discontinue-non-csv-outputs.md) | Retire 18 legacy output switches | CSV is the only output format the modern pipeline supports; legacy switches become deprecation warnings |
| [0010](adr/0010-macropore-deferral.md) | Macropore stays in legacy code; new pipeline does not wire it | Avoid blocking modernization on a complex deferred feature |
| [0011](adr/0011-macropore-exclusion-from-regression.md) | Drop case 3 (macroporeflow) from check-full | The only case exercising the deferred path; cleaner cut than runtime fallback |
| [0012](adr/0012-csv-companion-input-files.md) | `read_csv_table` is the canonical CSV reader | One reader for all tabular inputs (meteo, irrigation events, soil-initial, surface-water management) with ISO date support and typed errors |
| [0013](adr/0013-csv-meteorology-input.md) | TOML pipeline reads daily meteo from CSV | Drop the per-year `.YYY` file proliferation |
| [0014](adr/0014-readmeteo-phaseout.md) | Three-step phase-out of `readmeteo.f90` TTutil branches | Sub-daily CSV first, then delete TTutil branches, then sweep dead variables |
| [0015](adr/0015-strangler-narrow-scope-stub-errors.md) | Validator-level stub-errors for branches no case exercises | Defer breadth; cover depth; loud rejection beats silent miscompute |
| [0016](adr/0016-per-rotation-crop-config-cache.md) | Per-rotation crop content loaded into `crop_config_global` at config-load time | Avoid the per-rotation legacy reader cycle; future direction is explicit-argument passing |
| [0017](adr/0017-sibling-reader-dispatch-from-cache.md) | Sibling-reader dispatch around legacy `ArableLandGerm` | Cache-hit path uses typed config; cache-miss is fatal — no silent legacy fallback |
| [0018](adr/0018-fatalerr-shim.md) | Shim TTutil's `FatalERR` through `fatalerr_collected` | Fix the silent-exit-0 bug at the root; one-file linker shadow |
| [0019](adr/0019-legacy-readers-retired.md) | Leave legacy reader code in `src/`; document unreachability | Parity tests still need the readers; deletion deferred until parity tests retire |

## Test infrastructure (worth knowing)

The journey surfaced several pFUnit harness quirks that future work
inherits:

- **`tests/unit/run_pfunit.sh`** wraps the pFUnit binary and invokes
  each test suite in a separate process. Without this, state leakage
  in legacy `variables` globals between parity suites causes RDDATA
  fatal errors mid-suite. See ADR 0018 + `phase-4f-ss11-closeout.md`.
- **`call readswap('swap')`** — readswap accepts an optional
  `project_name` argument (added in SS-5). Parity tests pass `'swap'`
  explicitly so pFUnit's `--tap` / `-f` flags don't get interpreted as
  the project name by `Get_Command_Argument(1)`.
- **6 known parity divergences** between legacy and TOML are
  intentional (case 5 swinco; case 6 swdra=2 dramet/swdivd/swdislay;
  case-4 swoxygen=2 hlim1/2u/2l; case-5 ldis array; metfile basename
  comparison). Documented in the per-case parity audits.

## What's next (out of scope for Phase 4f-extend)

Pulled forward from the rescue plan and from sub-spec follow-ups:

| Topic | Origin | Effort |
|---|---|---|
| Compartment-state refactor | rescue plan §"Out of scope item 2" | L |
| Performance work (state-sync 344s macropore, vectorization, GPU, ifx re-enable) | rescue plan | L each |
| Python bindings (pyswap) | rescue plan | M-L |
| Iso-iso physical deletion of `readswap.f90` + dead `case(1)` blocks | ADR 0019 retirement gate | L (requires parity-test refactor first) |
| `iHWCKmodel` + `paramvg(10,:)` ksatexm threshold (Buckets B/C) | SS-10 audit | M each, on-demand |
| Refresh `phase-4f-config-to-variables-audit.md` per-section table totals | SS-9 audit | S |

## Pointers for new contributors

If you've just landed in this codebase:

1. **Run the gates** to make sure your environment is healthy:
   ```
   pixi run -e test build-linux
   pixi run -e test test-pfunit       # → Ok: 1, Fail: 0  (548 tests)
   pixi run -e test check-full        # → 5 passed, 0 failed
   ```
2. **Read [`architecture.md`](architecture.html)** for the current
   call-graph and I/O layer picture.
3. **Read [`configuration-schema.md`](configuration-schema.html)** for
   the TOML schema reference.
4. **Browse [`adr/`](adr/)** if you wonder *why* the code looks the
   way it does.
5. **Pop into [`archive/2026-phase-4/`](archive/2026-phase-4/README.html)**
   only when investigating a specific historical decision.

The legacy fixed-format readers in `src/io/readswap.f90`,
`src/crop/irrigation.f90` (case 1), `src/crop/management_soil.f90`
(SoilManagement(1)), and `src/crop/cropgrowth.f90` (nutrient block
+ readarablelandgerm) are **dead in the production runtime** but still
compiled. They serve as parity-test fixtures
(`tests/unit/io/toml/test_*_parity.pf`) that cross-check
legacy-vs-TOML output equivalence on the six regression cases.
