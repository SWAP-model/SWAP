# ADR 0045 — Orchestrator dissolution (`seed_state_from_config` deleted)

**Status:** Accepted (2026-05-28)
**Closes:** W1 + W2 from the 2026-05-28 model-setup-layer analysis.
**Predecessor:** ADR 0043 (two-phase init: type-bound init + free `x_seed`). ADR 0044 (typed CSV record tables) — the W3/W4 follow-on.

## Context

The strangler-pattern adapter `src/io/toml/seed_state_from_config.f90` was created to bridge the legacy `case(task)` lifecycle dispatcher and the modern typed-state initialization. It held cross-subsystem writes that no single `state%X%init` could own, plus a graveyard of historical comments documenting what had already been folded out (each `[GR-SEED Task N]` tag marked a prior migration).

Post-meteo-pilot and csv-families arcs, the module held only 8 live writes amid ~290 lines of historical comments. The W1 + W2 items from the model-setup-layer analysis identified two related symptoms:

- **W1**: manual post-init seeding blocks in `swap_mod.f90`'s `swap_init_body` (soilwater layer-flats, drainage `L`/`zbotdr`, etc.) — bypassed the `state%X%init` contract because of multi-source resolution rules.
- **W2**: `seed_state_from_config.f90` was 90% comment-graveyard.

Both were the same architectural debt: cross-subsystem writes living outside the type-bound init they belonged to.

## Decision

Fold every remaining write into the appropriate `state%X%init()` method, then delete `seed_state_from_config.f90` entirely.

The arc proceeded in 9 small steps, each one byte-identical, each one removing a write or block:

1. `mesh%numlay` from `config%soil%isoillay` → `state%mesh%init`
2. `surfacewater pondmx/rsro/rsroexp` → `state%surfacewater%init` (signature gains `config_soil`; init becomes unconditional with `swdra==2` gating heavy work internally)
3. `swdra` dependency relocation → `state%timecontrol%init` gains `config_drain` parameter; reads `swdra` directly from config rather than through `state%surfacewater` pre-write
4. Drainage post-init manual block (`swtopnrsrf`/`swdivdinf`/`FacDpthInf`/`L`/`zbotdr` with `dramet==2` branch) → `state%drainage%init`
5. Soilwater layer-flats STRANGLER (`ksatexm`/`ksatfit`/`cofani`/`flksatexm`/`orgmat`/`psand`/`psilt`/`pclay`/`swbotb_runtime` + `q0`/`k1max`/`H0max` zeroing) → `state%soilwater%init` (signature gains `config_drain` + `config_heat`; multi-source rules absorbed: `cofani` drain-then-soil-override, `orgmat` soil-then-heat.porg-backfill)
6. Snow handshake STRANGLER (`snowinco ↔ ssnow` for `swinco==3`) → `state%atmosphere%init`
7. `swinco==3` atmosphere warm-restart (`ssnow`/`ldwet`/`slw`/`atmin7` from `config%soil%initial`) → `state%atmosphere%init`
8. `swinco==3` solute cml profile (typed-table load + `cml_init`/`zc_init`/`nconc` populate) → `state%solute%init` (signature gains `config_soil`)
9. **Delete `seed_state_from_config.f90`** (module file + test file + build registrations). The single remaining `call state%timecontrol%init(...)` is inlined directly into `swap_init_body`.

## Consequences

+ **`swap_init_body` is the orchestrator.** No more separate adapter file; the init order is explicit in one place and reads top-to-bottom.
+ **Every state field has a single owner.** Each `state%X` field is written by exactly one `state%X%init` method (modulo intentional cross-cutting writes documented in code).
+ **Init signatures encode dependencies.** `state%timecontrol%init(config%simulation, config%general, config%drain)`, `state%soilwater%init(config%soil, config%drain, config%heat, config%bottom_boundary, ...)`, `state%solute%init(config%solute, config%soil, ...)` — the config sub-trees a subsystem reads are explicit at the call site.
+ **One fewer indirection layer for new developers.** No "where does this write come from?" → "look at the X subsystem's init" works for every field.
+ **All `state%cfg%X` reads inside init paths are gone or accounted for.** The remaining `state%cfg` reads are in compute paths (next arc).

- Init signatures grew. `state%soilwater%init` now takes 7 parameters; `state%solute%init` takes 3. The compute-path equivalent (`state%cfg` reach-through) was 2 parameters. Trade: explicit dependencies vs. less typing at call site. The former wins for maintainability.

- Some init methods are now unconditional even when their feature flag is off (`state%surfacewater%init` always runs; `swdra==2` gates the heavy work internally). This is correct architecturally — the lightweight scalars are always needed — but it's a small departure from the old "init only if you'll use it" pattern.

## Verification

- **812 → 806 pFUnit tests** (the 6 tests in `test_seed_state_from_config_suite` were deleted along with the module they tested).
- **check-fast PASS**: 4/4 byte-identical (hupselbrook, grassgrowth, surfacewater, salinitystress).
- **check-full PASS**: 5/5 non-xfail byte-identical (adds oxygenstress). Pre-existing xfails (winter, soilhysteresis) unchanged.
- 9 commits on branch `orchestrator-dissolution`.

## Next arc: `state%cfg` retirement

`state%soilwater%init`, `state%solute%init`, `state%surfacewater%init`, and `state%timecontrol%init` now all take their config sub-trees as parameters and read directly. The init layer is `state%cfg`-free.

The remaining `state%cfg%X` reads are in compute paths (`soilhydraulics.f90`, `temperature.f90`, `cropgrowth.f90`, etc.). The retirement strategy is per-subsystem: snapshot what each compute path needs into corresponding `state%X` fields, migrate read sites, drop the `state%cfg` pointer. ADR 0044 + ADR 0045 set up the conditions for this work — every CSV-derived value and every init-time config value already lives on state; the remaining work is compute-time config reads.
