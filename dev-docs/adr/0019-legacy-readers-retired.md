---
title: "ADR 0019 — Legacy readers retired from runtime; retained as parity fixtures"
date: 2026-05-05
status: accepted
---

# ADR 0019: Legacy readers retired from production runtime; retained as parity fixtures

## Context

Phase 4f-extend's umbrella spec
defined eleven sub-specs (SS-1 … SS-11) to retire the legacy
fixed-format `.swp` / `.dra` / `.crp` / `.YYY` / `.met` input readers
in favour of a single typed TOML pipeline. SS-1 through SS-10 are now
closed (commits `7b338c4` … `7879d71`). SS-11 is the closeout.

The remaining design question is what to do with the legacy reader
**source code**. Three options were considered:

1. **Delete entirely** — drop `src/io/readswap.f90`, `irrigation.f90`
   case(1), `management_soil.f90` SoilManagement(1..7),
   `cropgrowth.f90` nutrient block + `readarablelandgerm`. ~~3000 LoC
   removed~~.
2. **Move to a separate `tests/legacy/` tree** — keep the code
   compileable for parity tests but drop it from the production binary.
   Restructures meson build.
3. **Leave in place; document unreachability** — production call graph
   is verified clean (SS-1..9 audits); parity tests still invoke the
   legacy paths intentionally.

## Decision

**Option 3.** Leave the legacy reader code in place; document its
production-unreachability in `docs/architecture.md` and this ADR. The
parity-test infrastructure (`tests/unit/io/toml/test_*_parity.pf`,
`legacy_crop_helper.pf`, `parity_helpers.f90`) calls `readswap('swap')`
to drive the legacy reader against the six regression cases and
cross-check the resulting `variables` module globals against the typed
TOML config. This is a load-bearing safety net — it would be lost by
options 1 and 2 unless the parity tests are simultaneously refactored
to fixture-data comparisons (a non-trivial ~1000 LoC effort).

Verified by the SS-1..9 audits:

- **`readswap()` is not called from production source.** Confirmed by
  `grep -rn "call readswap\b" src/ --include='*.f90'` returning a
  comment-only match. Phase 4f strangler swap (`6abea3` and successors)
  removed the call from `src/core/swap.f90`.
- **`irrigation(1)`** is reachable only inside `readswap.f90` (gated on
  `SCHEDULE = 1` from the legacy `.crp`). Production-runtime
  unreachable. (SS-6 audit, `phase-4f-irrigation-audit.md`.)
- **`SoilManagement(1..7)`** is gated on `flCropNut`, which has no
  TOML adapter writer and defaults to `.false.`. The
  `tillage.f90:79` rejection (`fatalerr_collected('flCropNut = 1 not
  (yet) allowed')`) acts as belt-and-suspenders. (SS-7 audit,
  `phase-4f-management-soil-audit.md`.)
- **`cropgrowth.f90`'s nutrient init block (line 1151)** shares the
  `flCropNut` gate. (SS-8 audit, `phase-4f-cropgrowth-nutrients-audit.md`.)
- **`ArableLandGerm(1)`** is unreachable: ADR 0017's sibling-reader
  dispatch (`cropgrowth.f90:92-195`) treats cache-miss as a fatal
  error, not a fallback to `readarablelandgerm`.
- **All 146 G (gap) rows** from the 2026-04-27
  `phase-4f-config-to-variables-audit.md` were re-classified in SS-9:
  zero rows are both reachable in TOML runtime AND uncovered by
  schema/adapter (`phase-4f-ss9-residual-audit.md`).
- **All 12 HACK slots** in the strangler adapter
  (`config_to_variables.f90`) are either resolved (8/12), deferred
  with no-case-triggers rationale (2/12), reclassified as
  not-actually-a-hack (1/12), or logged as the SS-10.5 follow-up
  (1/12: `swpfile`/`logf`). (SS-10 closure,
  `phase-4f-config-to-variables-hacks-audit.md` 2026-05-05 update.)

## Consequences

**Production runtime:**
- Single canonical input path: `swap.toml` + companions →
  `load_swap_config` → `swap_config_t` → `config_to_variables` →
  legacy `variables` globals → physics solvers.
- Zero TTutil-based **data readers** in the production call graph.
- Three TTutil **utility** calls survive — `swap_main.f90:47,54`
  (`rdsets`, `rdfrom` for reruns) and `swapoutput.f90:3454` (`rddtmp`
  for cleanup). These are not data readers; they are TTutil's own
  bookkeeping primitives, load-bearing for the rerun mechanism.
- Two stub-readers survive — `Read_Tillage` (`tillage.f90:433`) and
  `SSDI_irrigation(1)` (`irrigation.f90:580`) both call
  `RDinit(unit, 0, swpfile)` to look up `swtill` / `swssdi` switches.
  Both default to 0 in every regression case. Retiring requires schema
  slots and adapter wiring; tracked as the SS-10.5 follow-up.

**Test infrastructure:**
- `readswap('swap')` is invoked by parity tests under
  `tests/unit/io/toml/`. The optional `project_name` argument added in
  SS-5 (commit `996476b`) lets these tests bypass the
  `Get_Command_Argument(1)` / pFUnit-CLI-flag collision.
- The `tests/unit/run_pfunit.sh` wrapper (ADR 0018) runs each test
  suite in a separate process to prevent legacy-reader state leakage
  between parity suites.
- TTutil's `FatalERR` is shimmed to `error_collection_t` /
  `fatalerr_collected` (ADR 0018) so legacy-reader fatal paths exit
  with status 1 instead of 0 — pFUnit harness now catches them.

**Code-shape:**
- `src/io/readswap.f90` (5080 LoC, 1017 TTutil calls), `case(1)`
  blocks in `irrigation.f90` / `management_soil.f90`, the nutrient
  block in `cropgrowth.f90`, and `readarablelandgerm` all remain in
  `src/`. They compile into the production binary but are statically
  unreachable.
- Future option to delete (Option 1 above) is preserved; it requires
  retiring the parity-test infrastructure first. Out of scope for
  Phase 4f-extend.

## Retirement gate (forward-looking)

When/if the parity tests are converted to fixture-data comparisons (no
longer needing the legacy reader binary), the legacy reader files can
be physically deleted from `src/`:

1. Convert all `tests/unit/io/toml/test_*_parity.pf` to read pre-cached
   `variables` snapshots from disk instead of invoking `readswap`.
2. Delete `tests/unit/io/toml/parity_helpers.f90`,
   `legacy_crop_helper.f90`, `readswap_stubs.f90`.
3. Delete `src/io/readswap.f90`. Delete `case(1)` from
   `irrigation.f90`, `management_soil.f90`. Delete the nutrient block
   from `cropgrowth.f90`. Delete `readarablelandgerm`.
4. Drop the corresponding entries from the meson `pfunit_extra_sources`
   list and the production `sources` list.
5. Keep TTutil as a build dependency only if the SS-10.5 follow-up
   isn't done yet (the rerun mechanism + the swtill/swssdi RDinit
   stubs both still depend on TTutil).

## Tag

`rescue/phase-4f-extend-complete` placed at the closing commit.

## Update 2026-05-06: physical deletion executed

The "retirement gate" forward-look above was executed in a follow-on
spec, landed across SS-A → SS-D on the `development` branch:

- **SS-A** — All five parity suites (hupselbrook, grassgrowth,
  oxygenstress, salinitystress, surfacewater) converted to literal-value
  assertions; helper files (`parity_helpers.f90`, `legacy_crop_helper.f90`,
  test_legacy_crop_helper suite) deleted.
- **SS-B** — Tillage and SSDI gating moved from internal self-checks to
  call-site flags (ADR 0020); SS-10.5 stub-errors retired so the
  subsystems are reachable when their flags are set.
- **SS-C** — Code deletion (this ADR's "Option 1"):
  - `src/io/readswap.f90` deleted (-5080 LoC).
  - `src/crop/cropgrowth.f90`: `ArableLandGerm(1)` case removed; legacy
    nutrient `rdinit` block deleted.
  - `src/crop/management_soil.f90`: `SoilManagement(1)` body collapsed
    to `return` (file-open + state-init was readswap's responsibility).
  - `src/crop/irrigation.f90`: `irrigation(1)` body collapsed to
    `return`.
  - `src/core/swap.f90`: dead `call SoilManagement(1)` removed.
  - `tests/unit/io/toml/readswap_stubs.f90` deleted (existed only to
    satisfy linker references from readswap).
  - `tests/unit/run_pfunit.sh` retired; meson `test()` invokes the
    pFUnit binary directly (no global state to leak between suites
    after readswap is gone).
  - Net: ~5,500 LoC removed from `src/`; one ~30-line `checkdate`
    helper extracted to `src/io/checkdate.f90` because
    `read_ssdi_input` (the swssdi=1 path) still calls it pending ADR
    0021 (TOML port of the SSDI block).
- **SS-D** — This update + sibling docs refresh.

End-state acceptance, all green:
- `src/io/readswap.f90` does not exist.
- `grep -rnE "subroutine readswap\b|call irrigation\(1\)|call SoilManagement\(1\)" src/` → no matches.
- `tests/unit/run_pfunit.sh` does not exist.
- `pixi run -e test test-pfunit` → `Ok: 1, Fail: 0` (540 tests, 1
  disabled = macroporeflow_parity per ADR 0011).
- `pixi run -e test check-full` → `5 passed, 0 failed` in ~42s.

Tag at the closing commit: `rescue/legacy-readers-deleted` (on
`development`; not yet merged to `main` per the post-Phase-4 release
gate).

