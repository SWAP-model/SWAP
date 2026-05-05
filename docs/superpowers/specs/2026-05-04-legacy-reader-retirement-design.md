---
title: "Legacy reader retirement — Phase 4f-extend complete (umbrella)"
author: Mateusz Zawadzki
date: 2026-05-04
status: draft
---

# Legacy reader retirement — Phase 4f-extend complete

Umbrella roadmap that finalizes the input-format port stage. After this work
the modern SWAP runtime never calls a TTutil-based legacy reader, every legacy
input key has a TOML schema mirror with verified finalized-data parity, and
the macropore switch is explicitly stub-errored at the TOML boundary. Each
sub-spec is independently shippable; the umbrella defines the methodology, the
sub-spec roadmap, and the retirement gate.

## Goal

Reach a state where:

1. The runtime call graph for the modern binary is free of TTutil reader calls
   (`rdinit`, `rdsdor`, `rdador`, `rdfdor`, `rdsror`, `rdaror`, …) outside of
   parity-test fixtures.
2. Every legacy `rd*` input key in `src/io/readswap.f90`,
   `src/io/readmeteo.f90` (TTutil branches), `src/crop/irrigation.f90`,
   `src/crop/management_soil.f90`, and `src/crop/cropgrowth.f90` has either:
   - a TOML schema mirror in `src/config/*_config.f90` populated by the typed
     reader pipeline, with the corresponding adapter wiring in
     `src/io/toml/config_to_variables.f90`, **or**
   - an explicit stub-error at the TOML boundary (the deferral case), **or**
   - an entry in `docs/adr/0009-discontinue-non-csv-outputs.md` (RETIRED) /
     `docs/adr/0010-macropore-deferral.md` (DEFERRED).
3. Every finalized-data transformation that the legacy reader performs is
   reproduced by the TOML adapter and verified — by regression for code paths
   the 6 test cases exercise, by **unit tests on pure-helper extractions** for
   code paths no test case exercises.
4. `swmacro = 1` in TOML triggers an explicit `ERR_VALIDATION_*` fatal error;
   the macropore reader code stays in tree (per ADR 0010) but cannot be
   reached from the TOML pipeline.
5. Adapter HACK marker count drops to zero, or each remaining marker is
   pinned to a deferral ADR.

End-state tag: `rescue/phase-4f-extend-complete`.

## Non-goals

- **Macropore physics re-port.** Case 3 stays excluded per ADR 0011. We add a
  TOML stub-error so authoring `swmacro=1` fails clearly, but we do not
  populate macropore globals from TOML. The legacy macropore reader code in
  `readswap.f90` lines 844+ stays in tree as orphan infrastructure.
- **Deletion of `readswap.f90` itself.** Sub-readers stay alive in the file as
  parity-test fixtures (per ADR 0015). A follow-up cleanup converts those
  parity tests to fixture-data comparisons and then deletes the file. That
  cleanup is out of scope here.
- **Compartment-state refactor.** Follow-on spec (per
  `docs/superpowers/specs/2026-04-22-rescue-and-stabilize-design.md` § Out of
  scope item 2). Phase 4f-extend completion unblocks it but does not start
  it.
- **Performance work** (state-sync 344s macropore, vectorization, GPU,
  Python bindings, ifx re-enable). All deferred to follow-on specs per the
  rescue & stabilize plan.

## Status snapshot (2026-05-04)

| Item | Status | Notes |
|---|---|---|
| Phase 4f strangler-swap | DONE | `readswap()` no longer called from `swap.f90` |
| `.crp` port phases 1–4 | DONE | Three crop families load from TOML init modules |
| `swap.ini` port | DONE | `[soil.initial]` + 3 CSV companions; `inifil` slot dropped (`7cea63d`) |
| `swap.dra` port | DONE | `surfacewater_init` replaces `rddre` (`27f8e65`) |
| CSV meteo Step 1 (sub-daily detail) | DONE | `MeteoCSVDetYear` + `detail_file` schema |
| CSV meteo Steps 2–3 (TTutil branch deletion + dead-var sweep) | DONE | SS-5 — TTutil branches retired (`2026-05-05-ss5-readmeteo-ttutil-deletion.md`) |
| Phase 0 test-dir reorg | DONE | `run_case.sh` defaults to in-place TOML |
| Adapter HACK Bucket A + B | DONE | 7 of 12 slots closed |
| Adapter HACK residuals (4) | OPEN | `iHWCKmodel`, `ksatexm`, `cropfil` suffix-strip, `swpfile`/`logf` |
| Macropore TOML stub-error | DONE | SS-4 — `swmacro=1` rejected by `soil_config_validate` (see `2026-05-04-ss4-macropore-toml-stub-error.md`) |
| `irrigation.f90` TTutil reader | DONE | SS-6 — closed by Phase 4f strangler swap; `irrigation(1)` reachable only via parity-test harness (see `phase-4f-irrigation-audit.md`) |
| `management_soil.f90` TTutil readers | DONE | SS-7 — unreachable via `flCropNut=.false.` default (see `phase-4f-management-soil-audit.md`) |
| `cropgrowth.f90` nutrient TTutil reader | DONE | SS-8 — same `flCropNut` gate as SS-7 (see `phase-4f-cropgrowth-nutrients-audit.md`) |
| `readswap.f90` full audit (unauthored keys) | OPEN | Existing audit only covered keys exercised by 6 cases |

## Methodology

### Audit format

Every sub-spec produces (or extends) a per-reader audit document. The
existing per-case audits and the HACKs audit set the format precedent:
`docs/phase-4f-hupselbrook-parity-audit.md`,
`docs/phase-4f-config-to-variables-hacks-audit.md`,
`docs/phase-4f-rddre-audit.md`,
`docs/phase-4f-readwofost-audit.md`. New audits follow the same column
layout:

| Legacy key | Status | Notes |
|---|---|---|

with the status legend extended:

- **OK** — TOML authors it; reader copies it; adapter writes the legacy
  global; finalize matches.
- **TOML** — schema slot exists; case TOML missing the value.
- **READER** — schema field exists; no TOML reader populates it.
- **ADAPTER** — schema + reader fine; adapter doesn't copy it.
- **FINALIZE** — adapter copies the value but skips a unit-conversion or
  derived transformation the legacy reader performs.
- **STUB_ERROR** — switch is rejected at the TOML boundary by an explicit
  fatal-error validator (the deferral case; macropore is the principal
  example).
- **RETIRED** — ADR 0009 (output switches collapsed into [output.csv] /
  forced-zero).
- **DEFERRED** — ADR 0010 (macropore) or ADR 0011 (regression-exclusion).
- **N/A** — switch is read but its value is irrelevant on the path taken.

The audit closes when every legacy key is OK, STUB_ERROR, RETIRED,
DEFERRED, or N/A.

### Parity-test mechanism for unexercised paths

For code paths no regression case exercises, parity is verified by **pure-helper extraction**:

1. Identify the finalize transformation in the legacy reader (e.g. a
   non-trivial conversion, a per-layer derivation, a unit change).
2. Port the transformation to a pure helper subroutine or function — no
   side effects, no I/O, no global access. Place it in the relevant
   `src/io/toml/` file or a new `src/io/toml/<domain>_finalize.f90`.
3. Wire the adapter to call the helper.
4. Add a pFUnit test under `tests/unit/io/toml/` that calls the helper with
   hand-computed inputs and asserts the outputs to documented numerical
   precision (typically 1e-12 relative for real64 arithmetic, exact for
   integer/string).
5. Reference the legacy reader line numbers in the test docstring so a
   future reviewer can verify the math against `readswap.f90` directly.

This replaces the alternative — synthesizing a minimal `.swp`+TOML pair and
running both pipelines — which was rejected because the modern source tree
no longer expects to open `.swp` files.

### Retirement gate

A legacy reader runtime call may be deleted only when **all** of these hold:

1. The reader's audit document is closed (every key OK / STUB_ERROR /
   RETIRED / DEFERRED / N/A).
2. Every `FINALIZE` row that was closed via pure-helper extraction has a
   passing unit test.
3. Every regression case affected by the reader passes
   `pixi run -e test regression <case>` at current tolerance.
4. The HACK marker count attributable to the reader drops to zero (or each
   remaining one points to a deferral ADR with rationale).

When the gate passes for a reader, its runtime call site is removed in a
single commit (mechanical change), and the reader's sub-functions stay alive
only if a parity test still depends on them.

## Sub-spec roadmap

Each sub-spec gets its own design + plan + execution cycle, following the
project's normal `superpowers:writing-plans` → `superpowers:executing-plans`
flow. Sub-specs marked DONE are listed for completeness only — no further
work expected.

| ID | Title | Status | Driver | Dependencies |
|---|---|---|---|---|
| SS-1 | `swap.ini` port | DONE | `2026-05-01-swap-ini-port-design.md` | — |
| SS-2 | `swap.dra` port | DONE | `2026-05-01-swap-dra-port-design.md` | — |
| SS-3 | CSV meteo finalize Step 1 (sub-daily) | DONE | `2026-05-01-csv-meteo-finalize.md` Phase 0 + Phase 1 | — |
| SS-4 | Macropore TOML stub-error | DONE | `2026-05-04-ss4-macropore-toml-stub-error.md` | — |
| SS-5 | `readmeteo.f90` TTutil branch deletion (ADR 0014 Steps 2–3) | DONE | `2026-05-05-ss5-readmeteo-ttutil-deletion.md` | — |
| SS-6 | `src/crop/irrigation.f90` TOML port | DONE | `phase-4f-irrigation-audit.md` (no code work — closed by Phase 4f strangler swap) | — |
| SS-7 | `src/crop/management_soil.f90` TOML port (SMM/SME/SNP) | DONE | `phase-4f-management-soil-audit.md` (no code work — `flCropNut=.false.` default) | — |
| SS-8 | `src/crop/cropgrowth.f90` nutrient block port | DONE | `phase-4f-cropgrowth-nutrients-audit.md` (same `flCropNut` gate as SS-7) | — |
| SS-9 | `readswap.f90` full audit + unauthored-key closure | OPEN | new spec | SS-4..8 close first (their keys overlap `readswap.f90`'s switch tree) |
| SS-10 | Adapter HACK residual cleanup | OPEN | new spec | — |
| SS-11 | Closeout — runtime-call deletion + retirement tag | OPEN | new spec | SS-4..10 all closed |

### Dependency graph

```
SS-1 ─┐
SS-2 ─┤  (already complete)
SS-3 ─┘

SS-4 ─┐
SS-5 ─┤
SS-6 ─┼──► SS-9 ──► SS-11
SS-7 ─┤
SS-8 ─┘

SS-10 ───────────────► SS-11
```

SS-4..8 are mutually independent. SS-9 walks `readswap.f90` end-to-end,
including any keys not picked up by SS-4..8; it must come after them so the
audit captures the actual residual surface. SS-10 (adapter HACK cleanup)
runs in parallel — it touches the adapter, not the legacy readers. SS-11 is
the closeout: deletes runtime calls, tags `rescue/phase-4f-extend-complete`.

### Per-sub-spec scope summary

**SS-4 — Macropore TOML stub-error.** Add a validator in
`soil_config_validate` (or a top-level `swap_config_validate` cross-section
check) that fatal-errors when `soil.swmacro = 1`. Update
`docs/adr/0010-macropore-deferral.md` with the new TOML behaviour.
Confirm that case 3 (excluded per ADR 0011) is unaffected and that no other
case authors `swmacro = 1`. Effort: S.

**SS-5 — `readmeteo.f90` TTutil branch deletion (ADR 0014 Steps 2–3).**
Steps 2 and 3 of ADR 0014 (TTutil reader removal + dead-variable sweep).
After SS-3 the runtime always takes the CSV path; the TTutil branches at
lines 112–132 and the daily-meteo `.met` fallback are unreachable. Delete
them, then sweep for now-dead variables (`swMetCSV`, `swRainCSV`,
`swMetFilAll`, `goto 100` label dance). Audit doc:
`docs/phase-4f-readmeteo-ttutil-audit.md`. Effort: M.

**SS-6 — `src/crop/irrigation.f90` TOML port.** Port the embedded
irrigation block (lines 79–326) currently read via `rdinit` from each
`.crp` file. The block contains `irrigevent`, `irgthreshold`, `irgdayfix`,
`irgdepmin`/`irgdepmax`, `dvs_tc1..tc8` + `trel/raw/taw/dwa/hcri/tcri`
arrays, `phFieldCapacity`, `dcrit`. Schema slot:
`crop.<rotation>.irrigation` sub-section in each crop config (cropfixed,
cropwofost, cropgrass). Wire the existing `salinitystress.irg.csv`
fixed-events file path through the new schema. Audit doc:
`docs/phase-4f-irrigation-audit.md`. Effort: M-L.

**SS-7 — `src/crop/management_soil.f90` TOML port.** Port the SMM
(soil-management measures), SME (soil-management events), and SNP
(soil-nutrient parameters) `rdinit` blocks at lines 129, 142, 197. These
read auxiliary `.crp`-companion files for organic matter / fertilization
management. Schema location TBD during the SS-7 design pass — likely a
`crop.management` sub-section or a new `[soil_management]` top-level
section. Audit doc: `docs/phase-4f-management-soil-audit.md`. Effort: M.

**SS-8 — `src/crop/cropgrowth.f90` nutrient block port.** Port the `rdinit`
at line 1155 reading nutrient inputs from a `.crp` file. Schema location
TBD during SS-8 design pass — likely `crop.<rotation>.nutrients` mirroring
the `FraHarLosOrm_*` keys. Audit doc:
`docs/phase-4f-cropgrowth-nutrients-audit.md`. Effort: M.

**SS-9 — `readswap.f90` full audit + unauthored-key closure.** Walk every
`call rd*` line in `readswap.f90` (16 distinct procedures). Cross-reference
against the existing audit docs and against the keys closed by SS-4..8.
For each unauthored key, decide: schema slot + reader + adapter, or
STUB_ERROR. Pure-helper extraction + unit tests for FINALIZE rows on
unexercised paths. Audit doc:
`docs/phase-4f-readswap-residual-audit.md`. Effort: L (the file is 5,080
lines, but most of it is already covered).

**SS-10 — Adapter HACK residual cleanup.** Close the four remaining
markers in `src/io/toml/config_to_variables.f90`:
- `iHWCKmodel` (line 716, Bucket C) — pin to a new ADR for deferred
  multi-model hydraulics, since no test case exercises models 2–11. The
  marker stays but cites the ADR.
- `ksatexm` threshold branch (line 743, Bucket B) — port the `flksatexm`
  branch from `readswap.f90:802–815` via pure-helper extraction. Unit-test
  the helper against hand-computed `relsatthr`/`ksatthr` values for a
  case with `ksatexm > ksatfit`.
- `cropfil` suffix-strip (line 1144, Bucket D) — closes naturally as
  SS-6/SS-8 complete the per-rotation crop adapter; remove the suffix
  strip, let `cropfil(i)` be empty, and drop the unused legacy global if
  the legacy crop sub-readers no longer need it.
- `swpfile` + `logf` (line 1197, Bucket D) — closes once SS-5 + SS-6 + SS-7
  + SS-8 retire all per-crop / drainage legacy readers that depend on
  `RDinit(unit, logf, swpfile)`. Effort: S–M (mostly mechanical after
  SS-4..9).

**SS-11 — Closeout.** When SS-4..10 are all closed:
1. Confirm zero TTutil calls in runtime src/ via grep.
2. Run full check-full + pFUnit; both green.
3. Update `docs/architecture.md` with the post-retirement I/O picture.
4. Write closing ADR (0018) — "TOML is the only runtime input path;
   legacy readers retained as parity-test fixtures only."
5. Tag `rescue/phase-4f-extend-complete`.

## Documentation deliverables

In addition to per-sub-spec audit docs:

- **Refresh `phase-4f-config-to-variables-hacks-audit.md`** — mark SS-1
  (CSV slot 11), SS-2 (rddre / surfacewater), and SS-3 (CSV meteo Step 1)
  as closed; document the SS-10 plan for the four remaining HACK slots.
- **Refresh `phase-4f-config-to-variables-audit.md`** — recompute the C/R/G
  classification after SS-4..9 close.
- **New ADR 0018** — "TOML-only runtime input path." Records the
  retirement decision and cites the parity-test fixture role of the
  surviving legacy reader code.
- **Update `docs/adr/0010-macropore-deferral.md`** — add the SS-4 stub-error
  at the TOML boundary as part of the deferral mechanism.

## Risks and mitigations

| Risk | Mitigation |
|---|---|
| SS-9 surfaces a switch combination that requires significant new schema work | Treat as a sub-sub-spec under SS-9; do not expand the umbrella. The audit drives the gap list; the gap list drives the work. |
| Pure-helper extraction breaks an unrelated path that secretly depended on the legacy reader's side effects | Run full check-full after every helper-extraction commit; the 5/5 regression suite is the primary backstop. |
| `cropfil` and `swpfile` HACKs (Bucket D) prove load-bearing for an unported reader missed in SS-4..8 | SS-9 catches the residual; SS-10 closes the HACKs only after SS-9 confirms zero callers. The dependency graph enforces ordering. |
| ADR 0014 Step 2–3 reveals more dead variables than expected | Sweep aggressively in SS-5; commit per file; don't batch. |
| Sub-spec count creates merge friction | SS-4..8 are mutually independent and small; each ships in its own session. SS-9..11 land sequentially. |

## Verification

End-to-end gate for the umbrella spec:

```bash
# 1. No TTutil calls in runtime sources outside parity-test fixtures.
grep -rn "call rd\(init\|sdor\|ador\|fdor\|sror\|aror\|scha\|acha\|ftim\|atim\|slog\|sinr\|finr\|ainr\)" src/ \
  | grep -v "readswap.f90\|readmeteo.f90:.*! parity-fixture"
# Expected: empty.

# 2. No HACK Phase 4f-extend markers without an ADR citation.
grep -n "HACK Phase 4f-extend" src/io/toml/config_to_variables.f90
# Expected: each remaining marker on the next line cites an ADR (e.g. ADR 0010, 0019).

# 3. All regression cases green via the modern binary.
pixi run -e test check-full
# Expected: 5 passed, 0 failed (case 3 excluded per ADR 0011).

# 4. pFUnit suite green, including new finalize-helper unit tests.
pixi run -e test test-pfunit
# Expected: OK.

# 5. Tag exists.
git tag --list rescue/phase-4f-extend-complete
# Expected: tag present.
```

## Definition of done

- All sub-specs (SS-4..11) closed and tagged in their own commits.
- Verification commands above all pass.
- ADR 0018 committed.
- `docs/phase-4f-config-to-variables-hacks-audit.md` and
  `docs/phase-4f-config-to-variables-audit.md` reflect post-retirement
  state.
- `rescue/phase-4f-extend-complete` tag present locally.

After this spec closes, the project is ready to start the
**compartment-based state refactor** follow-on (per
`docs/superpowers/specs/2026-04-22-rescue-and-stabilize-design.md` § Out of
scope item 2).
