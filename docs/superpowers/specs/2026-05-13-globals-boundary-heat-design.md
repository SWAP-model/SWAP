---
title: "GR-BH — Globals Retirement Arc 2+3 (Boundary + Heat) with Mesh Extraction"
date: 2026-05-13
status: approved
context: globals-retirement Arc 2 + Arc 3 + mesh extraction (combined), follows GR-UTILS
---

# GR-BH Design

## Goal

Retire `use variables` from `boundbottom`, `boundtop`, `frozencond`, `temperature` (4 files) by:

1. Extracting a `state%mesh` subrecord (new module) for the 7 mesh globals.
2. Extending `state%soilwater` with 8 soil-layer flats + 5 runtime scalars.
3. Extending `state%drainage` with 10 drainage geometry / switch fields.
4. Sweeping every other reader of those globals across the codebase.
5. Deleting the migrated globals from `variables.f90` and their populators from `config_to_variables.f90` / `initialize.f90`.
6. Relocating `heat_init` to a type-bound `state%heat%init(config%heat, numnod)` per the surfacewater_state pilot pattern.

Net deletion target: ~30 declarations from `variables.f90`, ~30 populator lines from `config_to_variables.f90`, zero-fills from `initialize.f90`, all 4 target files' bare `use variables` clauses (boundtop retains one narrow `use variables, only: nird` deferred to Arc 8).

## Architecture

Three sub-phases inside one arc; each phase ends with `pixi run check-full` 5/5 byte-for-byte.

- **Phase A — Schema + dual-write.** Additive. New module `mesh_state_mod`; schema extensions on `soilwater_state_t` and `drainage_state_t`; `config_to_variables.f90` populates both legacy globals and state mirrors. No reader change — byte-for-byte trivial.
- **Phase B — Reader cutover + heat_init relocation.** The 4 target files drop `use variables` and read state/config exclusively. `heat_init` moves to type-bound on `heat_state_mod`.
- **Phase C — Codebase mesh sweep + global retirement.** Every other reader (~30 files) of the migrated globals swaps to `state%X`. Then bare globals are deleted from `variables.f90`; legacy populator lines removed from `config_to_variables.f90` and `initialize.f90`. Compile-time errors surface any missed reader.

## Verification Discipline

Per `feedback_per_task_regression_gate.md` + `feedback_state_schema_clean_rebuild.md`:

- **Every implementer subagent ends with:** `rm -rf builddir && pixi run build-linux && pixi run test-pfunit && pixi run -e test python tests/regression/test_output_regression.py hupselbrook surfacewater salinitystress grassgrowth` — 4/4 byte-for-byte required.
- **Every phase ends with:** `pixi run check-full` 5/5.
- **Each commit migrates one file or one tight cohort.** Clean bisection lane preserved.
- **No skipped tests, no commented-out checks, no `--no-verify` commits.**

## Schema Additions

### NEW: `src/state/mesh_state.f90`

```fortran
module mesh_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use swap_array_dimensions, only: macp
   implicit none
   private
   public :: mesh_state_t

   type :: mesh_state_t
      integer :: numnod = 0                       !! number of soil compartments
      integer,      allocatable :: layer(:)       !! soil-layer index per compartment
      real(real64), allocatable :: dz(:)          !! compartment thickness [cm]
      real(real64), allocatable :: z(:)           !! mid-depth (negative downward) [cm]
      real(real64), allocatable :: disnod(:)      !! inter-node distance [cm]
      real(real64), allocatable :: ztopcp(:)      !! top of compartment [cm]
      real(real64), allocatable :: zbotcp(:)      !! bottom of compartment [cm]
   contains
      procedure :: init => mesh_init
   end type mesh_state_t

contains
   subroutine mesh_init(self, numnod_in, dz_in, z_in, disnod_in, ztopcp_in, zbotcp_in, layer_in)
      class(mesh_state_t), intent(inout) :: self
      integer, intent(in) :: numnod_in
      real(real64), intent(in) :: dz_in(:), z_in(:), disnod_in(:), ztopcp_in(:), zbotcp_in(:)
      integer,      intent(in) :: layer_in(:)
      ! allocate self%dz/z/disnod/ztopcp/zbotcp/layer to numnod_in, copy
      self%numnod = numnod_in
      ! ... allocate + assign each ...
   end subroutine
end module
```

Wired into `swap_state_t` as `type(mesh_state_t) :: mesh`.

### `soilwater_state_t` extensions

| Field | Type | Semantics |
|---|---|---|
| `ksatexm(:)` | real(real64), allocatable | Layer Ksat (examined extension) |
| `ksatfit(:)` | real(real64), allocatable | Layer fitted Ksat |
| `cofani(:)` | real(real64), allocatable | Layer anisotropy |
| `flksatexm` | logical | Global flag — Ksatexm present in input |
| `orgmat(:)` | real(real64), allocatable | Layer organic matter (gravimetric) |
| `psand(:)` | real(real64), allocatable | Layer sand fraction |
| `psilt(:)` | real(real64), allocatable | Layer silt fraction |
| `pclay(:)` | real(real64), allocatable | Layer clay fraction |
| `q0` | real(real64) | Runtime surface flux (boundtop ↔ pondrunoff) |
| `k1max` | real(real64) | Runtime max conductivity at z=0 |
| `H0max` | real(real64) | Runtime max ponding pre-runoff |
| `swbotb_runtime` | integer | Runtime-overridable bottom-boundary switch (legacy mutates to -2) |

`ArMpSs` is NOT migrated — it's a macropore retirement artifact, always 0. Phase B simply deletes its remaining assignment and global declaration.

### `drainage_state_t` extensions

| Field | Type | Semantics |
|---|---|---|
| `nrlevs` | integer | Number of drainage levels |
| `swdra` | integer | Drainage method switch |
| `swdivd` | integer | Drainage subdivision switch |
| `swnrsrf` | integer | NRSRF table switch |
| `swtopnrsrf` | integer | NRSRF top switch |
| `swdivdinf` | integer | Infiltration subdivision switch |
| `FacDpthInf` | real(real64) | Infiltration depth factor |
| `L(:)` | real(real64), allocatable | Drainage spacing per level |
| `zbotdr(:)` | real(real64), allocatable | Drainage depth per level |
| `owltab(:)` | real(real64), allocatable | Open-water level table |

### `heat_state_t` — no new fields

Adds `procedure :: init => heat_state_init` only; relocates the existing `heat_init` body. Signature becomes `state%heat%init(config%heat, numnod)`. The `config%heat` arg is unused for now (no fields read), but lands the signature for future heat-init seed migration.

### swbotb runtime mutation

Legacy code mutates the bare global `swbotb` to `-2` at runtime when bottom node goes oven-dry (boundbottom.f90:102). This persists across timesteps — genuine state. **Decision: `swbotb_runtime` on `state%soilwater`, seeded from `config%bottom_boundary%swbotb` at init.** Readers use `state%soilwater%swbotb_runtime`. Config value stays immutable.

## Phase A — Tasks (10)

Schema + dual-write. Each task ends with full regression gate.

| # | Task | Files |
|---|---|---|
| A1 | Pre-flight baseline | (verify only) — record baseline SHA, `pixi run check-full` 5/5 |
| A2 | Create `mesh_state_mod` | `src/state/mesh_state.f90` (NEW), `src/state/swap_state.f90` (add field), `meson.build` |
| A3 | `state%mesh` dual-write populator | `src/io/toml/config_to_variables.f90` — call `state%mesh%init(...)` after mesh globals are populated |
| A4 | Extend `soilwater_state_t` — layer flats (8 fields) | `src/state/soilwater_state.f90` |
| A5 | Extend `soilwater_state_t` — runtime scalars (4 fields) | `src/state/soilwater_state.f90` (q0, k1max, H0max, swbotb_runtime) |
| A6 | Dual-write soilwater layer flats | `src/io/toml/config_to_variables.f90` |
| A7 | Dual-write soilwater runtime scalars | `src/io/toml/config_to_variables.f90` (seed `swbotb_runtime = config%bottom_boundary%swbotb`) |
| A8 | Extend `drainage_state_t` (10 fields) | `src/state/drainage_state.f90` |
| A9 | Dual-write drainage geometry | `src/io/toml/config_to_variables.f90` |
| A10 | Phase A close | (verify only) — `pixi run check-full` 5/5; commit `arc-phase-A-complete` marker |

**Phase A close criteria:**
- All schema fields populated alongside legacy globals.
- `grep "state%mesh%numnod" src/` returns ≥1 hit (the populator).
- No reader changes — every legacy global still load-bearing.
- check-full 5/5 byte-for-byte.

## Phase B — Tasks (13)

Reader cutover + heat_init relocation.

| # | Task | Files | Notes |
|---|---|---|---|
| B1 | Relocate `heat_init` to type-bound | `src/state/heat_state.f90`, `src/heat/temperature.f90` (delete), `src/core/swap_mod.f90` (call site) | New call: `state%heat%init(config%heat, state%mesh%numnod)` |
| B2 | Migrate `FrozenCond` | `src/heat/frozencond.f90` | Drop `use variables`. Mesh→`state%mesh%`. `swfrost`/`tfroststa`/`tfrostend`→ `config` arg + `config%soil%frost%X` / `config%heat%X` |
| B3 | Migrate `FrozenBounds` | `src/heat/frozencond.f90` | Drop second `use variables`. Mesh→`state%mesh%`; drainage geo→`state%drainage%`; layer flats→`state%soilwater%` |
| B4 | Verify `frozencond.f90` clean | grep | `grep "^[[:space:]]*use variables" src/heat/frozencond.f90` empty |
| B5 | Migrate `temperature` | `src/heat/temperature.f90` | Drop `use variables`. Mesh→`state%mesh%`; config switches→`config%heat%X` + `config%simulation%swinco`; soil composition→`state%soilwater%`; `Tav`/`atav` already on `state%atmosphere` |
| B6 | Migrate `Devries` helper | `src/heat/temperature.f90` | Replace `use variables, only: NumNod` with explicit `numnod` arg passed by caller |
| B7 | Verify `temperature.f90` clean | grep | empty |
| B8 | Migrate `BoundBottom` | `src/boundary/boundbottom.f90` | Drop `use variables`. Mesh→`state%mesh%`; config refs (`gwltab`, `qbotab`, `haqtab`, `hbotab`, `sw2/3/4`, sine params, etc.)→ `config` arg + `config%bottom_boundary%X`. `swbotb` reads/writes→`state%soilwater%swbotb_runtime` |
| B9 | Verify `boundbottom.f90` clean | grep | empty |
| B10 | Migrate `boundtop` | `src/boundary/boundtop.f90` | Drop `use variables`. Mesh→`state%mesh%`; runtime (`q0`,`k1max`,`H0max`)→`state%soilwater%`; `ArMpSs` deleted; config switches→`config%X`; layer flats→`state%soilwater%`. **Retains:** `use variables, only: nird` (irrigation, Arc 8 territory) |
| B11 | Migrate `PONDRUNOFF` | `src/boundary/boundtop.f90` | Drop the `only:` block. Same replacement set |
| B12 | Verify `boundtop.f90` clean | grep | Empty EXCEPT the documented `nird` deferred line |
| B13 | Phase B close | (verify only) — check-full 5/5; commit `arc-phase-B-complete` marker |

**Phase B close criteria:**
- `grep -n "^[[:space:]]*use variables" src/boundary/*.f90 src/heat/*.f90` returns at most one expected hit (`nird` in boundtop).
- heat_init lives on `heat_state_mod`; legacy symbol deleted from `temperature.f90`.
- check-full 5/5.

## Phase C — Tasks (16)

Codebase mesh sweep + global retirement. Phase C does NOT touch any other `use variables` symbols. Files keep their `use variables` clause for unrelated globals (those belong to Arcs 4–8). Surgical replacement only.

**Contingency:** if a file reads `numnod` but doesn't already take `state`, add `type(swap_state_t), intent(in) :: state` to its signature. All existing callers hold state.

| # | Task | Files | Replace |
|---|---|---|---|
| C1 | Cluster: `src/core/` | `initialize.f90`, `timecontrol_mod.f90`, `swap_bmi_mod.f90`, `swap_mod.f90`, `swap_capi_mod.f90` | Mesh→`state%mesh%X` |
| C2 | Cluster: `src/atmosphere/` | `et.f90`, `interception.f90`, `meteoday.f90`, `meteodt.f90` | Mesh→`state%mesh%X` |
| C3 | Cluster: `src/io/` | `readmeteo.f90`, `swap_csv_output.f90`, `swapoutput.f90` | Mesh→`state%mesh%X`; `nrlevs` where present→`state%drainage%X` |
| C4 | Cluster: `src/soil/` | `soilgrid.f90`, `soilhydraulics.f90`, `waterbalance.f90` | Mesh→`state%mesh%X`; layer flats (`orgmat`/`psand`/`psilt`/`pclay`/`ksatexm`/`ksatfit`/`cofani`/`flksatexm`)→`state%soilwater%X` |
| C5 | Cluster: `src/drainage/` | `drainage.f90`, `surfacewater.f90`, `divdra.f90` | Mesh→`state%mesh%X`; drainage geo (`nrlevs`/`swdra`/`swdivd`/`swnrsrf`/`swtopnrsrf`/`swdivdinf`/`FacDpthInf`/`L`/`zbotdr`/`owltab`)→`state%drainage%X`; layer flats→`state%soilwater%X` |
| C6 | Cluster: `src/crop/` init files | `cropfixed_init.f90`, `cropgrass_init.f90`, `cropwofost_init.f90`, `wofost_soil_parameters.f90` | Mesh + relevant layer flats only |
| C7 | Cluster: `src/crop/` runtime | `cropgrowth.f90`, `irrigation.f90`, `oxygenstress.f90`, `rootextraction.f90`, `tillage.f90` | Mesh + relevant layer flats |
| C8 | Audit pass — find stale reads | (grep-only audit) | `grep -rnw "numnod\|dz\|disnod\|ztopcp\|zbotcp" src/ --include="*.f90" \| grep -v "state%\|config%\|variables.f90\|mesh_state\|! "` — any bare reference missed is fixed in this task. **Mandatory gate before C9.** |
| C9 | Drop mesh writes from `config_to_variables.f90` | `src/io/toml/config_to_variables.f90` | Delete legacy `numnod`/`dz`/`z`/`disnod`/`ztopcp`/`zbotcp`/`layer` mirrors. Keep `state%mesh%init(...)` call |
| C10 | Drop layer-flats + runtime-state writes | `src/io/toml/config_to_variables.f90` | Delete legacy `ksatexm`/`ksatfit`/`cofani`/`flksatexm`/`orgmat`/`psand`/`psilt`/`pclay`/`q0`/`k1max`/`H0max`/`swbotb` mirrors |
| C11 | Drop drainage-geo writes | `src/io/toml/config_to_variables.f90` | Delete legacy `nrlevs`/`swdra`/`swdivd`/`swnrsrf`/`swtopnrsrf`/`swdivdinf`/`FacDpthInf`/`L`/`zbotdr`/`owltab` mirrors |
| C12 | Drop zero-fills from `initialize.f90` | `src/core/initialize.f90` | Delete zero-fills for now-deleted globals |
| C13 | Delete mesh globals from `variables.f90` | `src/core/variables.f90` | Delete declarations of `numnod`, `dz(:)`, `z(:)`, `disnod(:)`, `ztopcp(:)`, `zbotcp(:)`, `layer(:)`. Compile surfaces any missed reader |
| C14 | Delete layer-flat + runtime globals | `src/core/variables.f90` | Delete `ksatexm`/`ksatfit`/`cofani`/`flksatexm`/`orgmat`/`psand`/`psilt`/`pclay`/`q0`/`k1max`/`H0max`/`ArMpSs`/`swbotb` |
| C15 | Delete drainage-geo globals | `src/core/variables.f90` | Delete `nrlevs`/`swdra`/`swdivd`/`swnrsrf`/`swtopnrsrf`/`swdivdinf`/`FacDpthInf`/`L(:)`/`zbotdr(:)`/`owltab(:)` |
| C16 | Final verification + arc-complete marker | (verify only) | `pixi run check-full` 5/5; greps: `state%mesh%` ≥30 hits, `state%soilwater%ksatexm` ≥5 hits, `state%drainage%nrlevs` ≥5 hits; `grep -nw "numnod" src/core/variables.f90` empty; commit `arc-bh-complete` |

**Phase C close criteria:**
- `grep -nw "numnod\|dz\|z\|disnod\|ztopcp\|zbotcp\|layer" src/core/variables.f90` empty for mesh
- Layer-flat + runtime + drainage-geo globals deleted from `variables.f90`
- `pixi run check-full` 5/5 byte-for-byte

## Risk & Mitigation

| Risk | Mitigation |
|---|---|
| Phase C C13–C15 deletion exposes missed readers | C8 audit pass is mandatory gate; deletions surface remaining gaps as compile errors |
| Mesh blast radius (~30 files) creates merge complexity | Cluster commits (C1–C7); each is independent and bisectable |
| Signature changes for files not taking state | Documented contingency: add `intent(in)` arg; all callers already hold state |
| State schema changes break incremental builds | `rm -rf builddir` mandatory per `feedback_state_schema_clean_rebuild.md` — every implementer prompt includes this |
| Subtle byte-for-byte regression in Richards solver path | Per-task regression gate catches immediately; bisectable to single-file commit |

## What this Arc DOES NOT Change

- BMI / cffi surfaces (`swap_bmi_mod`, `swap_capi_mod`) — they continue working throughout. They depend only on `state` and `config`, which is exactly what this arc extends.
- The Python demo `tests/cffi-demo/run_ensemble.py` — unchanged.
- Physics — every value computed by every formula stays bit-identical. Only data-access paths change.
- Other `use variables` consumers (atmosphere/io/soil/drainage/crop) keep their bare-globals readers for OTHER symbols not in scope here. Those are Arcs 4–8 territory.

## Effort

- **Total tasks:** 39 (A1–A10 + B1–B13 + C1–C16)
- **Estimated effort:** 6–8 days, subagent-driven
- **Phases gated by check-full 5/5 byte-for-byte**

## Memory & ADR Consequences

- ADR candidate: ADR 0042 supersedes — new ADR documenting the mesh subrecord extraction and the soilwater/drainage extensions.
- After arc closes, update memory `project_state_rescue_complete_2026-05-12.md` listing GR-BH as completed.
- The mesh extraction completes the largest single architectural debt remaining toward Arc 9 — after GR-BH, the remaining `use variables` readers are concentrated in atmosphere/soil/drainage/crop runtime, and the bare globals in `variables.f90` are dominated by atmosphere/crop data.

## Files & Artifacts

- **This spec:** `docs/superpowers/specs/2026-05-13-globals-boundary-heat-design.md`
- **Plan (next):** `docs/superpowers/plans/2026-05-13-globals-boundary-heat.md`
- **Roadmap context:** `docs/superpowers/plans/2026-05-13-globals-retirement-roadmap.md` (this arc combines roadmap Arcs 2 + 3 + mesh extraction)
- **GR-UTILS precedent:** `docs/superpowers/specs/2026-05-13-globals-utils-design.md` (Arc 1, completed 2026-05-12)
