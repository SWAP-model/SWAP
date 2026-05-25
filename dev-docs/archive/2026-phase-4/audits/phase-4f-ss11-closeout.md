# Phase 4f-extend SS-11 — closeout

**Date:** 2026-05-05
**Sub-spec:** SS-11 (legacy reader retirement umbrella, final)
**Tag:** `rescue/phase-4f-extend-complete`
**ADR:** 0019 — Legacy readers retired from runtime; retained as parity fixtures

## Summary

Phase 4f-extend is complete. The production SWAP binary uses TOML as
its only input format. All eleven sub-specs of the umbrella spec
(`docs/superpowers/specs/2026-05-04-legacy-reader-retirement-design.md`)
have closed.

## Sub-spec roll-up

| SS | Title | Closure |
|---|---|---|
| SS-1 | `swap.ini` port | DONE pre-umbrella (`7cea63d`) |
| SS-2 | `swap.dra` port | DONE pre-umbrella (`27f8e65`) |
| SS-3 Step 1 | CSV meteo Step 1 (sub-daily detail) | DONE pre-umbrella (`MeteoCSVDetYear`) |
| SS-3 Steps 2–3 | CSV meteo TTutil branch deletion + dead-var sweep | DONE — SS-5 (`bd62344`/`abd66f0`/`3037bef`) |
| SS-4 | Macropore TOML stub-error | DONE (`3196165`) |
| SS-5 | `readmeteo.f90` TTutil branch deletion (ADR 0014 Steps 2–3) | DONE (3 commits) |
| SS-6 | `irrigation.f90` TOML port | DONE — closed by Phase 4f strangler swap (audit only, `9ed4a00`) |
| SS-7 | `management_soil.f90` TOML port | DONE — `flCropNut=.false.` default makes it unreachable (audit only, `1a951ec`) |
| SS-8 | `cropgrowth.f90` nutrient block port | DONE — same `flCropNut` gate (audit only, `1a951ec`) |
| SS-9 | `readswap.f90` full audit + unauthored-key closure | DONE — zero STILL-G rows (`94e1bdc`) |
| SS-10 | Adapter HACK residual cleanup | DONE — zero actionable slots (`7879d71`) |
| SS-11 | Closeout — runtime-call deletion + retirement tag | **DONE** (this doc) |

## Acceptance criteria

| Criterion | Result |
|---|---|
| Zero TTutil **data-reader** calls reachable from production main loop | ✅ — verified via SS-1..9 audits + per-sub-spec gate checks |
| Two TTutil **utility** calls survive | Documented (rerun mechanism + temp-file cleanup); not data readers |
| Two stub-readers survive (`Read_Tillage`, `SSDI_irrigation(1)`) | Documented as SS-10.5 follow-up; both default to 0 in every case |
| `pixi run -e test test-pfunit` | ✅ 546 passed, 0 failed |
| `pixi run -e test check-full` | ✅ 5 passed, 0 failed |
| `docs/architecture.md` updated | ✅ I/O layer section rewritten |
| Closing ADR | ✅ ADR 0019 |
| Retirement tag | ✅ `rescue/phase-4f-extend-complete` |

## Test gate snapshot

```
$ pixi run -e test test-pfunit | grep -E "^Ok:|^Fail:"
Ok:                 1
Fail:               0

$ pixi run -e test check-full | grep -E "passed|failed"
pFUnit summary: 546 passed, 0 failed
Results: 5 passed, 0 failed
```

## Outstanding follow-ups (NOT gating Phase 4f-extend)

| Task | Origin | Effort |
|---|---|---|
| SS-10.5: schema-port `swtill`/`swssdi`; retire `RDinit(swpfile)` in `tillage.f90:433` and `irrigation.f90:580`; drop the swpfile/logf hack | SS-10 audit | M |
| SS-5 follow-ups M1-M5 (test helper refactor, comment tightening) | SS-5 code review | S |
| Iso-iso physical deletion of `readswap.f90` + dead `case(1)` blocks | ADR 0019 retirement gate | L (requires parity-test refactor first) |
| Refresh `phase-4f-config-to-variables-audit.md` per-section table totals | SS-9 audit | S (mechanical) |
| Add SS-10 schema slots for runtime feature flags (`fldumpconvcrit`, `flMaxIterTime`, `flprintdt`, `flSwapShared`) if needed | SS-9 audit | S each, on-demand |

## Pointers for the post-Phase-4f-extend reader

- **Modern input pipeline:** `src/io/toml/load_swap_config.f90` →
  `src/config/swap_config.f90` → `src/io/toml/config_to_variables.f90`.
- **Per-rotation crop init:** `src/crop/{cropfixed,cropwofost,cropgrass}_init.f90`.
- **Adapter audit:** `docs/phase-4f-config-to-variables-audit.md`
  (553 rows; 2026-04-27, refreshed by SS-9).
- **Hack tracker:** `docs/phase-4f-config-to-variables-hacks-audit.md`
  (12 original slots; 2026-05-05 update).
- **Per-sub-spec audits:** `docs/phase-4f-{irrigation,management-soil,
  cropgrowth-nutrients,readmeteo-ttutil,rddre,…}-audit.md`.
- **Test harness:** `tests/unit/run_pfunit.sh` (per-suite isolation,
  ADR 0018) → `unit-swap-tests` binary → `--tap` per suite.
- **Closing ADR:** `docs/adr/0019-legacy-readers-retired.md`.
