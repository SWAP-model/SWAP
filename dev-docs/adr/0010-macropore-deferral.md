---
title: "ADR 0010 — Macropore module deferral"
date: 2026-04-29
status: accepted
---

# ADR 0010: Macropore module deferral

## Context

Phase 4f-prep was scheduled to deliver `macropore_config_t` (22-key
schema), `read_macropore_toml`, wiring into `swap_config_t`, plus
per-case TOML and parity assertions for case 3 (3.macroporeflow). The
new TOML path would then drive macropore physics for case 3 once
Phase 4f's strangler-fig replaced `readswap()`.

While Phase 4f-prep was in flight, the upstream macropore developers
indicated that the current implementation is being phased out for
reliability reasons. Future macropore work will rebuild the physics
on the existing legacy code as a foundation, but no Phase-4-era
investment in TOML parity for the present implementation is
warranted.

## Decision

**The current macropore implementation stays in place but is
"switched off" in the new TOML pipeline.** Concretely:

1. **`macropore_config_t` is kept** as orphan infrastructure under
   `src/config/macropore_config.f90` (committed `87bea98`). The
   22-key schema, validators, and 9 unit tests remain in the
   codebase, ready for future macropore work to wire up.
2. **No `read_macropore_toml` module is authored** in Phase 4f-prep.
3. **No `[macropore]` field is added to `swap_config_t`.** Nothing
   in the new TOML pipeline references `macropore_config_t`.
4. **No per-case TOML extension** for case 3 (3.macroporeflow). Its
   `tests/swap-cases/toml/3.macroporeflow/swap.toml` carries no
   `[macropore]` block.
5. **No parity assertions** for macropore fields in case 3's
   parity test. The existing parity test already scopes to the
   schema-covered subset (per Phase 4e Task C3); no change needed.
6. **Phase 4f's strangler-fig keeps `readswap.f90` as a fallback**
   for macropore-using cases. The runtime branch is:
   ```
   if (toml_path_available) then
      call load_swap_config(...); call config_to_variables(config)
   else
      call readswap()   ! legacy fallback (case 3 with SWMACRO=1)
   end if
   ```
   Case 3 retains its `swap_linux.swp.template` and runs on the
   legacy reader. Cases 1, 2, 4, 5, 6 take the new TOML path.

   **SUPERSEDED by [ADR 0011](0011-macropore-exclusion-from-regression.md):**
   case 3 is now excluded from `check-full` entirely; the runtime
   branch is dropped and Phase 4f makes a clean cut to the new
   TOML pipeline. `readswap.f90` still stays in the tree (per the
   broader Phase 4f design) as Phase 4f-extend's reference.
7. **Audit-doc reclassification.** The 10 G entries in the macropore
   section of the Phase 4f config-to-variables audit (since retired)
   reclassify from G (gap) to **DEFERRED** — a new fourth status
   distinct from RETIRED. Semantics: "schema exists, wiring deferred
   to a future phase," vs RETIRED's "going away forever per ADR
   0009."

## Consequences

Positive:

- Phase 4f-prep skips three sub-tasks (A2, A3, D1) that would have
  produced throw-away wiring once the macropore implementation
  changes.
- Phase 4f's strangler-fig has a clean cutover criterion: cases
  with `swmacro=0` use the new path; the one case with `swmacro=1`
  keeps using the legacy reader. No partial-coverage edge cases.
- Future macropore work has a fully-typed schema waiting in
  `src/config/macropore_config.f90` to extend or replace.

Negative:

- `src/io/readswap.f90` is NOT deletable in Phase 4f. It stays
  ~5,000 LoC of legacy Fortran in the tree until the future
  macropore phase resolves.
- `macropore_config_t` is orphan code (compiled, tested, never
  consumed) until that phase. Some readers may find the asymmetry
  confusing without this ADR as context.
- Case 3's parity test asserts the schema-covered subset only
  (general/simulation/meteorology/drainage/soil/bottom_boundary/
  crop_meta + wintcer1 crop). Macropore physics regressions caught
  only by `check-full`, not by unit-level parity.

## Revisit trigger

When upstream commits to the next-generation macropore
implementation. At that point: extend or replace
`macropore_config_t`, author the reader, wire it into
`swap_config_t`, populate case 3's TOML with `[macropore]`, extend
the parity test, and finally remove `readswap.f90`.

## Update 2026-05-04 — TOML-boundary stub-error (Phase 4f-extend SS-4)

Authoring `soil.swmacro = 1` in a TOML configuration is now rejected at
validation time by `soil_config_validate` (in `src/config/soil_config.f90`).
The error message reads:

> soil.swmacro=1 (macropore physics) not yet supported in the TOML
> pipeline; case 3 is excluded from regression per ADR 0011 and the
> macropore module remains deferred per ADR 0010.

Rationale: the modern binary's adapter copies `swmacro` into the legacy
global without populating any other macropore state. Without the
stub-error, a user who authored `swmacro = 1` would get silent runtime
corruption (or a crash deep in macropore physics that reads
unallocated arrays) instead of an immediate, actionable failure.

This update does not change the deferral itself: `macropore_config_t`
remains orphan infrastructure, no `read_macropore_toml` module exists,
and no `[macropore]` field is wired into `swap_config_t`. The
stub-error is the TOML-side counterpart to the regression exclusion
recorded in ADR 0011 — both close the macropore path cleanly without
removing the legacy code.

When future macropore work re-enables the module, this stub-error
must be removed in the same change that wires `[macropore]` into the
TOML pipeline.
