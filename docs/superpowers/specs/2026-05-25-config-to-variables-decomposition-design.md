# config_to_variables.f90 Decomposition — Design Spec

**Date:** 2026-05-25
**Author:** zawadzkim + Claude
**Status:** Draft — awaiting review

## Background

After the multi-week bare-globals retirement arc (terminal commits Phase 6
Steps 1–6, 2026-05-25), `src/io/toml/config_to_variables.f90` is the last
piece of the strangler pattern. It started life as a dual-write shim
(TOML → typed config → bare globals), but every bare-global write has now
been retired. Today it is a **legitimate TOML → typed state adapter**: it
reads `config%X` (typed, TOML-driven) and writes `state%X` (typed runtime).

However, the file is still a 1631-line god-object — every subsystem's
config→state seeding logic lives in one place, including non-trivial work
(CSV reads, derived computations, allocations, sub-record dispatch). The
right home for each block is the corresponding `state%X%init(config%Y)`
method on the typed state. This spec captures the decomposition plan.

## Goals

1. Move every computation out of `config_to_variables.f90` into the
   typed init method of the subsystem it belongs to.
2. Where a subsystem has no `init` method yet, create one with the
   pattern `state%X%init(config%Y, …)` (type-bound procedure).
3. Where a subsystem has a free `<subsystem>_init` subroutine (older
   pattern), promote it to type-bound + extend its signature so it
   absorbs the work currently in the adapter.
4. Replace `config_to_variables(config, state)` with a thin orchestrator
   (~20–30 lines) that just dispatches to subsystem inits in the correct
   order.
5. **Rename** the orchestrator's module/file to reflect the new role
   (e.g. `seed_state_from_config_mod` → `src/io/toml/seed_state_from_config.f90`).
6. Preserve byte-identical regression at every commit (the standard
   check-fast 5/5 gate).

## Non-goals

- Changing the TOML schema or the `*_config_t` types.
- Reworking the readers in `src/io/toml/read_*_toml.f90` (those are fine).
- Eliminating the `apply_*` helpers' algorithms — only their location and
  signatures change.
- Cleanup of `src/state/legacy_state.f90` (separate small sweep).

## Current state — inventory

### Main adapter sections

| Section | Lines | Destination subsystem | Computations (non-trivial work) |
|---|---|---|---|
| General + simulation | 70–124 | `state%timecontrol` | • Derive `iyear`/`imonth` from `tstart` via `dtdpar`. • Allocate `outdat`/`outdatint` to MAOUT. • Call `populate_outdatint_monthly` when `swmonth==1`. |
| Simulation.numerical | 126–136 | `state%timecontrol` | Plain assignments (dt, dtmin, dtmax, MaxIt, msteps, MaxIterTime/flMaxIterTime defaults). |
| Meteorology | 137–263 | `state%atmosphere` | • Daily metfile CSV → `metcsv_dat`/`nmetcsv` (9 cols, lowercase-extension guard). • Detail metfile CSV → `metcsv_det`/`nmetcsv_det` (7 cols, `swmetdetail==1`). • Rain events CSV → `raincsv_dat`/`nraincsv` (2 cols, `swrain==3`). • Snow params scalars (`TePrRain`, `TePrSnow`). |
| Drainage | 265–418 | `state%drainage` (existing `drainage_init`) + `state%surfacewater` | • DRAMET=2 dispatch (`ipos` chain → `khtop/khbot/zintf/kvtop/kvbot/geofac`). • Per-level arrays (`swdtyp`, `swallo`). • owltab CSV pre-load (per drainage level → `state%drainage%owltab(nrlevs, 2*MAOWL)` interleaved). • `surface_runoff` sub-section: broadcast `swtopdislay`/`ftopdislay` to MADR. |
| Soil | 420–595 | `state%soilwater` (existing `soilwater_init`) + `state%mesh%numlay` | • Reads `numlay` from `config%soil%isoillay`. • `bdens` element-wise copy. • `swinco==3` warm-restart: `h_file` CSV (z,h) → `config%soil%initial%z_init`; `atmin7` scalar; optional Cml CSV → `state%solute%{cml_init,zc_init}`. |
| Bottom boundary | 597–795 | TBD — case dispatch | Per-case CSV loads: case(1) gwl_file→`gwltab`; case(2) `qbot2_file`→`qbotab`; case(3) `haquif_file`→`haqtab` + optional `qbot4_file`→`qbotab`; case(4) `qhbot_file`→`qbotab` (with `abs(htab)`); case(5) `hbot5_file`→`hbotab`; case(6,7,8) no-op. All write into `state%soilwater%X` (existing arrays). |
| Heat | 797–816 | `state%heat` (existing `heat_state_init`) | No active work — all legacy mirrors retired. Already done. |
| Irrigation (fixed) | 818–881 | `state%crop%irrigation` | • Fixed events from inline array OR CSV (`swirfix==1`). • mm→cm conversion on `irdepth`. • `nirri_fixed=1` reset. |
| Solute | 883–957 | `state%solute` (existing `solute_init`) | Per-layer arrays (`ldis`, `kf`, `decpot`, `fdepth`) — broadcast scalar to ldis(1) when no array. Flatten 2D `cseeptab` to interleaved afgen layout. |
| Surface water (management periods) | 959–1046 | `state%surfacewater` (existing `surfacewater_state_init`) | Per-period arrays (`impend`, `swman`, `wscap`, `wldip`, `intwl`, `hbweir`, `alphaw`, `betaw`) — guarded allocate-then-copy. |
| Crop | 1048–1087 | `state%crop%common` + `state%crop` (existing `crop_state_init`) | • If `swcrop==1` set `flCropReadFile`/`flCropOpenFile`. • Allocate + copy `state%crop%common%croptype` from `config%crop%rotation_type`. • Point `crop_config_global` to `config%crop` (legacy ABI). |
| Output switches | 1089–1112 | `state%timecontrol%swheader` | Plain assignment. |

### Helper subroutines

| Helper | Lines | Destination |
|---|---|---|
| `populate_outdatint_monthly` | 1122–1155 | Private inside `timecontrol_state_mod` — called by `state%timecontrol%init`. |
| `parse_iso_date_to_days1900`, `strip_crp_toml_suffix` | small | Either keep as module helpers, or relocate to `src/io/toml/toml_field_helpers_mod` (date helper) and `src/io/toml/read_crop_toml` (stem stripper). |
| `apply_soil_tillage` | 1196–1288 | `state%tillage%init(config%soil%tillage, tend)` — promote `tillage_init` to type-bound + absorb. |
| `apply_irrigation_ssdi` | 1291–1354 | Private inside `crop_irrigation_state_mod` — called by `state%crop%irrigation%init`. |
| `apply_ssdi_mode0`/`apply_ssdi_mode1` | 1502–1629 | Same as above — private helpers. |
| `apply_nutrients` | 1357–1386 | `state%nutrients%init(config%nutrients, pathwork)` — extend existing init signature. |
| `apply_nutrients_events` | 1389–1499 | Same — private helper. |

### Init landscape

| Subsystem | Has init? | Current signature | Refactor verdict |
|---|---|---|---|
| `atmosphere` | ✓ type-bound | `init(self, config)` | **EXTEND** — absorb meteo CSV pre-loads + snow params. |
| `crop_common` | ✓ type-bound | `init(self)` | **EXTEND** — accept config; absorb crop flags + croptype copy. |
| `crop_irrigation` | ✓ type-bound (stub) | `init(self)` | **EXTEND** — accept config; absorb fixed-irrigation + SSDI logic. |
| `crop` (umbrella) | ✓ type-bound | `init(self, crop_cfg)` | **EXTEND** — call child inits in order; handle crop_config_global. |
| `crop_fixed`/`grass`/`wofost`/`oxygen` | ✓ type-bound | various | KEEP — already encapsulated, no work in adapter. |
| `heat` | ✓ type-bound | `init(self, heat_cfg, numnod)` | KEEP — already fully encapsulated. |
| `mesh` | ✓ type-bound | `init(self, nn, dz, z, …)` | KEEP — explicit-args pattern; called from CalcGrid in swap_mod. |
| `nutrients` | ✓ type-bound | `init(self, nlay)` | **EXTEND** — add `config_nut` + `pathwork` args; absorb `apply_nutrients` logic. |
| `surfacewater` | ✓ type-bound | `init(self, config_sw, config_drain, numnod)` | **EXTEND** — absorb management-period arrays. |
| `drainage` | ✗ free sub | `drainage_init(state, config)` | **PROMOTE** to `state%drainage%init(config_drain, mesh_numnod)`; absorb owltab CSV + surface_runoff broadcasting. |
| `soilwater` | ✗ free sub | `soilwater_init(self, numnod, numlay)` | **PROMOTE** to type-bound; extend signature with `config_soil` to absorb bdens/swfrost/swinco=3 logic. |
| `solute` | ✗ free sub | `solute_init(state)` | **PROMOTE** to `state%solute%init(config_solute, numlay)`; absorb per-layer broadcast + cseeptab flatten. |
| `tillage` | ✗ free sub | `tillage_init(self, numlay)` | **PROMOTE** to `state%tillage%init(config_soil_tillage, tend, numlay)`; absorb `apply_soil_tillage`. |
| `timecontrol` | ✗ free sub | `timecontrol_init(state)` (already exists in `src/core/timecontrol_mod.f90`!) | **CREATE** new `state%timecontrol%init(config_simulation, config_general)` for *seeding* (distinct from the existing per-run `timecontrol_init` which is a daily-step initializer). |
| bottom_boundary | ✗ | — | **DECIDE** — either fold into `state%soilwater%init` (since tables write to `state%soilwater%X`) or create a small dedicated helper. |

## Proposed target architecture

### Top-level orchestrator (replaces `config_to_variables`)

`src/io/toml/seed_state_from_config.f90`:

```fortran
module seed_state_from_config_mod
   use swap_config_mod, only: swap_config_t
   use swap_state_mod,  only: swap_state_t
   implicit none
   private
   public :: seed_state_from_config

contains

   subroutine seed_state_from_config(config, state)
      type(swap_config_t), target, intent(inout) :: config
      type(swap_state_t),          intent(inout) :: state

      ! Top-level wiring — must run first.
      state%cfg => config

      ! Order matters: mesh and timecontrol seed the size/time context
      ! that downstream inits need.
      call state%timecontrol%init(config%simulation, config%general)
      call state%mesh%init_from_config(config%soil)            ! deferred — currently CalcGrid
      call state%atmosphere%init(config)                       ! extended: absorbs CSV pre-loads
      call state%heat%init(config%heat, state%mesh%numnod)
      call state%soilwater%init(config%soil, state%mesh%numnod, state%mesh%numlay)
      call state%solute%init(config%solute, state%mesh%numlay)
      call state%drainage%init(config%drain, state%mesh%numnod)
      call state%surfacewater%init(config%surface_water, config%drain, state%mesh%numnod)
      call state%tillage%init(config%soil%tillage, state%timecontrol%tend, state%mesh%numlay)
      call state%nutrients%init(config%nutrients, config%general%pathwork, state%mesh%numlay)
      call state%crop%init(config%crop, config%general%pathwork)
   end subroutine seed_state_from_config

end module seed_state_from_config_mod
```

The exact arg list per init is finalized during implementation. Some calls
may take `config` whole instead of a sub-record where they need cross-section
data (e.g. `state%atmosphere%init` already does this).

### Call-chain flattening (Task 12)

The current chain is **three deep** with two layers conflated:

```
swap_main → swap_init → swap_init_from_loaded_config → config_to_variables
```

`swap_init_from_loaded_config` (~330 lines) is two responsibilities mashed
together:

- **Layer A — config → state seeding.** Everything seed_state_from_config
  absorbs: the big `config_to_variables` call, plus the explicit
  `state%X%init(...)` calls already split out (atmosphere, heat, nutrients,
  soilwater, tillage, drainage).
- **Layer B — first-day compute init.** Work that must run *after* state is
  fully seeded: `timecontrol_init` (derives day-1 flags from seeded
  `tstart/iyear`), `CalcGrid` (builds mesh from layer config),
  `SoilWater(1, state)`, `SwapOutput(1, state)`,
  `SoilWaterOutput(1, state, config)`.

After Tasks 1–11 collapse Layer A into `seed_state_from_config`,
Layer B is ~5 inline calls. The `swap_init_from_loaded_config` wrapper
becomes trivially small and can be deleted, with its remaining body
folded directly into `swap_init`:

```fortran
subroutine swap_init(config_file, state, config)
   ! Phase 1: load + validate + finalize TOML
   call load_swap_config(config_file, config, errors)
   call config%validate(errors)
   call config%finalize(errors)
   call errors%abort_if_fatal()

   ! Phase 2: seed state from config (pure config→state, no compute)
   call seed_state_from_config(config, state)

   ! Phase 3: first-day compute init (needs state fully seeded)
   call timecontrol_init(state)
   call CalcGrid(state, config)
   call SoilWater(1, state)
   call SwapOutput(1, state)
   call SoilWaterOutput(1, state, config)

   call log_info('swap', &
        'Initialization complete for project: ' // trim(config%general%project))
end subroutine swap_init
```

End-state call chain:

```
swap_main → swap_init(config_file, state, config)
              ├─ load + validate + finalize config
              ├─ seed_state_from_config(config, state)   ← Layer A
              └─ [Layer B: 5 inline calls]
```

Two routines instead of three; phases are explicit. `swap_init_from_loaded_config`
retires (with a deprecation note if any external consumer of the public API
needs warning).

### Phasing

The decomposition is naturally per-subsystem. Each task is:

1. Identify the section of code in `config_to_variables` that belongs to subsystem X.
2. Create or extend `state%X%init` to accept the right args + absorb that logic.
3. Move private helpers into the state's module (e.g. `apply_ssdi_mode0/1` → `crop_irrigation_state_mod`).
4. Replace the section in `config_to_variables` with a single `call state%X%init(...)`.
5. check-fast 5/5 byte-identical.
6. Commit.

After all sections are migrated, `config_to_variables` is reduced to the thin
orchestrator above. Then a final commit:
7. Rename file/module to `seed_state_from_config`.
8. Update the one caller in `swap_mod.f90`.
9. Update test files.

### Suggested task order (low-risk first)

| Task | Subsystem | Risk | Notes |
|---|---|---|---|
| 1 | `state%timecontrol` (new init) | Low | Self-contained section; depends on nothing else. |
| 2 | `state%atmosphere` (extend init) | Low | Existing init already takes `config`; just absorb CSV pre-loads. |
| 3 | `state%solute` (promote + extend) | Low | Self-contained per-layer arrays + table flatten. |
| 4 | `state%nutrients` (extend) | Medium | Absorbs `apply_nutrients` + `apply_nutrients_events` (CSV + sort + group). |
| 5 | `state%crop%irrigation` (extend) | Medium | Absorbs `apply_irrigation_ssdi` + mode0/mode1 helpers + fixed events. |
| 6 | `state%surfacewater` (extend) | Medium | Per-period arrays; large slot. |
| 7 | `state%drainage` (promote + extend) | Medium-high | Owltab CSV + surface_runoff broadcasting. |
| 8 | `state%soilwater` (promote + extend) | High | Largest subsystem, plus bottom-boundary case dispatch. |
| 9 | `state%tillage` (promote + extend) | Medium | Absorbs `apply_soil_tillage`. |
| 10 | `state%crop` (extend) | Medium | Umbrella + croptype + global pointer. |
| 11 | Final orchestrator rename + cleanup | Low | `config_to_variables` → `seed_state_from_config`, file/module rename, test update. |
| 12 | Flatten swap_init: retire `swap_init_from_loaded_config` | Low | Fold Layer B into `swap_init` directly. Call chain becomes 2-deep (`swap_main → swap_init → seed_state_from_config`). |

Risk is mostly proportional to amount of moved code; high-risk items have
non-trivial computations (CSVs, allocations, sort, dispatch).

## Open questions

- **`bottom_boundary` home.** Currently writes into `state%soilwater%X`. Two
  choices: fold into `state%soilwater%init`, or create `state%soilwater%init_bottom_boundary`
  as a separate method called sequentially. Recommendation: **fold in**, with
  a `select case (swbotb)` private helper for clarity.

- **`crop_config_global` pointer.** Set at line ~1083 of the current adapter
  (`crop_config_global => config%crop`). This is a legacy ABI hook for code
  that doesn't yet take `config` as arg. Should the rename also retire this
  pointer? Likely yes, but as a separate follow-up — out of this spec's scope.

- **`config%soil%initial%z_init` write-back.** Current adapter mutates
  `config` (in-place write of `z_init` from h_file CSV). After decomposition,
  this needs an explicit `intent(inout)` on whichever init absorbs it (likely
  `state%soilwater%init`). Naming: keep `config_soil` as `intent(inout)` for
  that init, document the mutation.

- **`swap_config_t%target` attribute on the adapter input.** Current adapter
  takes `config` with `target` attribute (for `state%cfg => config` pointer).
  After decomposition, the orchestrator still needs `target`. Internal subsystem
  inits don't, but should accept `intent(in)` (or `inout` for the soilwater
  exception above). Each init's intent is decided per-task.

- **Helpers' new homes.**
  - `populate_outdatint_monthly` → private helper inside `timecontrol_state_mod`.
  - `parse_iso_date_to_days1900` → relocate to `toml_field_helpers_mod` (shared utility).
  - `strip_crp_toml_suffix` → relocate to `read_crop_toml_mod` (only user there).

## Acceptance criteria

- After full decomposition: `src/io/toml/seed_state_from_config.f90` is
  ≤ 50 lines (orchestrator + module boilerplate).
- Every subsystem with active work in the old adapter has a type-bound
  `init` method on its state type, called from the orchestrator.
- check-fast 5/5 byte-identical at every commit throughout the
  decomposition.
- One commit per task; spec-compliant `[GR-SEED YYYY-MM-DD]` tag on each
  commit.
- File renamed in a single dedicated commit at the end (no semantic
  changes in the rename commit).

## Risks

- **Init-order dependencies.** Some inits need `state%mesh%numnod` /
  `state%timecontrol%tend` before they can run. The orchestrator order
  in the architecture sketch above is the validated order. Any subsystem
  init that needs a value from another subsystem's init must take that
  value as an explicit arg, not reach into `state`.
- **Test fallout.** Several pFUnit tests construct `swap_state_t` and call
  inits directly. Their signatures will change; test calls must update too.
  Most tests are small (one or two inits each) so this is manageable.
- **`crop_config_global` legacy pointer.** Some code paths read this
  module-level pointer rather than taking `config` as arg. The orchestrator
  must continue to set the pointer until those readers are migrated.

## Out of scope (next-arc follow-ups)

- Retire `crop_config_global` legacy pointer (separate refactor).
- Clean up `src/state/legacy_state.f90` retirement comments.
- Delete `nheat`/`zh`/`tsoil` legacy decls in heat_state_mod once
  temperature.f90's `tab` local replaces them entirely.
- Migrate the bottom_boundary case-dispatch into typed `bottom_boundary_state_t`
  if a clearer subsystem boundary is wanted.
