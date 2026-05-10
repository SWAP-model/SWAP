# Solute State-Type Migration — Design

**Date:** 2026-05-10
**Status:** accepted (pending implementation)
**Migration #:** 3 of N
**Branch:** `refactor/solute-state` (new branch from development; subsequent subsystem migrations may continue on this branch or fork further)
**Discovery:** `docs/superpowers/specs/2026-05-10-state-migration-solute-discovery.md`
**Predecessor playbook:** ADR 0030 (surfacewater) + ADR 0031 (drainage)
**Successor ADR:** 0032

## Goal

Migrate the solute subsystem off `variables.f90` globals into a typed `solute_state_t` aggregated under `swap_state_t`. Extract `AgeTracer` (currently dead — `flAgeTracer` never set) into a separate module with a stub-error so it can be reactivated cleanly later. Promote ~14 physics fields currently missing from `solute_config_t` into the typed config (Phase 0), patching a real correctness gap where TOML runs silently use zeroed defaults for core Freundlich/decomposition/breakthrough physics.

## Out of scope

- Adding a regression test case that exercises `swsolu=1` with non-trivial physics. None of the 5 current regression cases activates solute transport. Phase 0 fixes the schema gap; coverage is a separate workstream.
- Re-enabling the AgeTracer feature (currently dead). Stub-error preserves the code; reactivation is a future arc.
- Full solute coverage of `swsoltyp=2` (pesticide) branches if they touch globals not listed in Section 2 of discovery.
- Resolving working-buffer ownership of `ArMpSs` (co-written by soilhydraulics/boundtop/solute). It stays as a shared global until macropore migration sorts it.

## Context

The solute subsystem (`solute.f90`, 548 LoC, plus `AgeTracer` co-located there) computes solute transport, sorption, decomposition, and tracer aging. The discovery doc cataloged 37 owned globals (12 of which are AgeTracer-specific), 47 borrowed globals, 7 external reader files, 10 cross-subsystem hazards.

Distinctive concerns vs surfacewater/drainage:

- **AgeTracer is dead code.** `flAgeTracer` is never set to `.true.` anywhere in `src/`. The feature has been dormant.
- **Major physics config gap.** 14 fields consumed by `solute.f90` are absent from `solute_config_t` and the adapter. TOML-driven runs silently use zero defaults for `cref`, `kf`, `frexp`, `gampar`, `decpot`, `fdepth`, `ddif`, `kfsat`, `poros`, `daquif`, `decsat`, `cseeptab`, `swbr`, plus `cpre` (concentration in precipitation). This is a real correctness bug for any TOML run that activates `swsolu=1` and would silently produce wrong outputs.
- **`cml` dual-use** (solute concentration / age storage). Resolved by extracting AgeTracer.
- **Two-stage `cml` seeding** at config-time (CSV) + runtime `afgen` re-interpolation. Migration must preserve.

## Decisions

### D1. Pilot scope: solute subsystem home tree

The three home files: `src/solute/solute.f90`, `src/config/solute_config.f90`, `src/io/toml/read_solute_toml.f90`. Plus call-site adjustments wherever the subsystem is invoked.

`AgeTracer` extracts to a NEW home file `src/solute/agetracer.f90` (Phase 1).

### D2. State-type shape

`solute_state_t` in `src/state/solute_state.f90` holds the solute-owned variables. AgeTracer-specific globals are NOT in this state — they go with the AgeTracer module (Phase 1, Decision D5).

Field set (preliminary — implementer refines from discovery Section 2):

```fortran
type :: solute_state_t
   real(real64), allocatable :: cml(:)        ! soil solute concentration (per node)
   real(real64), allocatable :: cmsy(:)       ! solute mass storage (per node, soil + solid)
   real(real64), allocatable :: cnh4(:), cno3(:)  ! ammonia, nitrate per node
   ! cumulative balance fields
   real(real64) :: samini = 0.0_real64        ! initial total mass
   real(real64) :: sampro = 0.0_real64        ! cumulative production
   real(real64) :: samcra = 0.0_real64        ! cumulative crop uptake
   real(real64) :: solbal = 0.0_real64        ! balance error
   ! per-step intermediates
   real(real64) :: rottot = 0.0_real64        ! root uptake total
   real(real64) :: sqprec = 0.0_real64, sqirrig = 0.0_real64
   real(real64) :: sqdra  = 0.0_real64, sqbot = 0.0_real64
   real(real64) :: dectot = 0.0_real64
   ! AgeTracer fields are NOT here — they go with the AgeTracer module
end type solute_state_t
```

Implementer derives the precise list from discovery doc Section 2 (37 owned globals minus the 12 AgeTracer-specific).

### D3. State-type aggregation

`src/state/swap_state.f90` adds `type(solute_state_t) :: solute`.

### D4. Argument-threading

`solute(task, state)` already takes state (added during surfacewater Phase 2 Task 5). The signature stays. `AgeTracer(task, state)` will move to its new module (Phase 1).

ASSOCIATE in compute body. Use `sl_*` prefix where bare names would shadow `use Variables` imports (drainage's lesson — gfortran's host-association rules tolerate the shadow but the explicit prefix is clearer).

### D5. AgeTracer extraction + stub-error

Move the `AgeTracer` subroutine from `src/solute/solute.f90` to a new file `src/solute/agetracer.f90` declaring a `module agetracer_mod`. Move the AgeTracer-specific globals (12 fields per discovery) into the new module's documentation as a "would-be `agetracer_state_t`" placeholder, but **don't migrate them off variables.f90 yet**. The runtime path stub-errors if `flAgeTracer = .true.`.

Pattern (mirroring tillage/SSDI before they had typed configs):

```fortran
subroutine AgeTracer(task, state)
   if (flAgeTracer) then
      call fatalerr_collected('AgeTracer', &
         'AgeTracer feature is not currently supported on the TOML path. ' // &
         'It was extracted from solute.f90 during ADR 0032 with the design ' // &
         'decision deferred. To re-enable, define agetracer_state_t, migrate ' // &
         'the 12 owned globals (currently in variables.f90), and remove this ' // &
         'guard. See discovery doc Section 8 hazard #4 for context.')
      return
   end if
   ! ...stub body kept for source-completeness...
end subroutine
```

The 12 AgeTracer-specific globals stay declared in `variables.f90` (commented as "AgeTracer-specific; future agetracer_state_t target"). They're not written by anything because the gate is closed. After this arc:
- `outage` (the AgeTracer output routine in `swapoutput.f90`) gets the same stub-error treatment OR is gated behind `flAgeTracer` (which is always false).
- The dual-use of `cml` disappears — solute owns it cleanly.

### D6. Phase 0: promote 14+ missing physics fields into `solute_config_t`

Currently absent fields (per discovery hazard #6):

| Field | Type | Physical meaning |
|---|---|---|
| `cref` | real(real64) | Reference concentration (Freundlich) |
| `kf` | real(real64) | Freundlich coefficient |
| `frexp` | real(real64) | Freundlich exponent |
| `gampar` | real(real64) | Decomposition rate parameter |
| `decpot` | real(real64) | Potential decomposition rate |
| `fdepth` | real(real64) (allocatable per layer?) | Depth correction for decomposition |
| `ddif` | real(real64) | Molecular diffusion coefficient |
| `kfsat` | real(real64) | Saturated-zone Freundlich coefficient |
| `poros` | real(real64) | Porosity (solute-specific) |
| `daquif` | real(real64) | Aquifer thickness |
| `decsat` | real(real64) | Saturated-zone decomposition rate |
| `cseeptab` | real(real64) (allocatable, table?) | Seepage concentration vs time |
| `swbr` | integer | Breakthrough curve switch |
| `cpre` | real(real64) | Solute concentration in precipitation |

Plus any siblings the implementer finds during the schema audit. Each gets:
- A field in `solute_config_t` with a documented default and validation range
- A line in the TOML reader to populate it
- A line in the `apply_solute` adapter to write the legacy global at config-load time

Phase 0 doesn't add regression coverage — none of the 5 cases activates `swsolu=1`. The fix prevents future TOML runs from silently using zeros.

For tabular fields (`fdepth`, `cseeptab`), use the established CSV-companion pattern (ADR 0012) if they're large; inline TOML arrays if small. Implementer's judgment.

### D7. Cross-subsystem reader migration (Phase 2)

7 external reader files per discovery Section 3.5:

| Reader file | Fields read | Category |
|---|---|---|
| `src/io/swapoutput.f90` (outvap, outend, outage, outbal, outsba) | `cml`, `cmsy`, `samcra`, `sampro`, `solbal`, etc. | Output |
| `src/io/swap_csv_output.f90` (set_values, per-node CSV) | `cml(:)`, others | Output |
| `src/crop/rootextraction.f90` (salt/osmotic stress) | `cml` | Compute |
| `src/crop/irrigation.f90` (concentration threshold) | `cml` | Compute |
| `src/soil/soilhydraulics.f90`, `src/boundary/boundtop.f90` | `ArMpSs` (NOT solute-owned — shared working buffer) | (skip — not migrated) |
| `src/io/toml/config_to_variables.f90` | `cml(k)` init seed | Init-seed |

Phase 2 migrates each compute/output reader to `state%solute%*`. The init-seed pattern (CSV → `cml(k)` at config time) preserves; the runtime `afgen` re-interpolation in `solute(task=1)` switches to writing `state%solute%cml`.

### D8. AgeTracer-related output routines

`outage` in swapoutput.f90 reads AgeTracer-specific globals. Two options:

- **(A) Gate behind flAgeTracer** (which is always false today). Don't migrate; let it stay as dead code parallel to AgeTracer.
- **(B) Stub-error inside outage** like the AgeTracer compute routine.

Recommend **(A)**: gate behind `if (flAgeTracer) then ... end if`. Cleaner than stub-error since `outage` was never reachable from the TOML path anyway. The `outage` body keeps reading legacy globals — they're never written so the read returns zeros, which is fine because the gate prevents write-out.

### D9. `ArMpSs` stays as a shared global

Discovery hazard #3: `ArMpSs` is co-written by `soilhydraulics.f90`, `boundtop.f90`, and `solute.f90`. Document in solute_state_t's file comment that `ArMpSs` is intentionally NOT in the state — it's a shared working buffer. A future macropore migration may sort the ownership.

### D10. `cpre` and `cseeptab` are config-time (Phase 0)

`cpre` is a constant solute concentration in precipitation; `cseeptab` is a time-table of seepage concentration. Both belong in config, not state. Phase 0 promotes them.

## Phasing

Three phases, each shipping with check-full byte-identical (5/5) and pFUnit green:

- **Phase 0 — Physics config gap.** Promote 14+ missing fields into `solute_config_t` + TOML reader + adapter. Each TOML field gets a default + validation. No regression test (none of the 5 cases activates swsolu=1) but pFUnit covers the validators. Closes the silent-zero-defaults bug for any future TOML run with `swsolu=1`.
- **Phase 1 — State type + AgeTracer extraction + threading + dual-write + drop.** Define `solute_state_t`. Add to `swap_state_t`. Extract AgeTracer to `src/solute/agetracer.f90` with stub-error. Gate `outage` behind `flAgeTracer`. Solute compute writes state alongside globals (dual-write). Output reads switch to typed state. Drop dual-write — state authoritative.
- **Phase 2 — Cross-subsystem reader migration + global removal + ADR finalize.** Migrate `rootextraction`, `irrigation` (compute readers of `cml`), `swapoutput`, `swap_csv_output` (output readers) to `state%solute`. Drop legacy global writes in solute compute. Comment out (with provenance markers) the 25 solute-owned globals in `variables.f90` — keep the 12 AgeTracer-specific declared (gated; future agetracer migration removes them).

## Testing

- **pFUnit:** new `test_solute_state` (default values + array allocation lifecycle). New tests for the 14 promoted config fields' validators. Migrate any existing solute tests that read globals.
- **check-full:** 5/5 byte-identical at every commit. The 5 cases don't activate solute, so this gate proves the migration doesn't break the inactive path. Active-path testing is a separate future workstream.

## Non-goals (explicit)

- Rewriting solute physics. ASSOCIATE preserves variable names; bodies of `solute(task=1/2)` should not change semantically.
- Adding regression coverage for `swsolu=1`. Phase 0 fixes the schema gap; activation testing is a future arc.
- Migrating the AgeTracer-specific globals. They stay in variables.f90 (commented as future targets) until the AgeTracer feature is reactivated.

## ADR 0032

`docs/adr/0032-state-migration-solute.md` records:
- Decision: solute is migration #3; subsystems migrate one at a time.
- AgeTracer extraction with stub-error rationale (dead code preservation).
- Phase 0 physics config patch — the discovery's correctness finding.
- Cross-references to discovery doc, ADRs 0030/0031, plan files.
