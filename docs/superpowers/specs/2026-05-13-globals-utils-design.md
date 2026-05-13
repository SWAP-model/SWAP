---
title: "Globals Retirement — Arc 1: utils — Design Spec"
date: 2026-05-13
status: draft
arc: GR-UTILS (globals retirement, utils cluster)
roadmap: docs/superpowers/plans/2026-05-13-globals-retirement-roadmap.md
---

# Globals Retirement Arc 1 — `src/utils/`

## Context

Arc 1 of the Globals Retirement Roadmap. Two utility files (`src/utils/soilhydraulicsutils.f90`, `src/utils/surfacewaterutils.f90`) still import legacy globals via `use variables, only: ...`. Additionally, `soilhydraulicsutils` carries the `bind_state_targets` / `bind_tc_target` module-level-pointer pattern (binding `cofgen`, `fluseksatexm`, `tc_dt_ptr` from state at init time). Both are bridges over the same chasm — the underlying functions don't take `state`, so they reach for it via module-level state instead of explicit args.

This arc retires both. Functions that need state are refactored to take the relevant subrecord (`state%soilwater`, `state%surfacewater`) as an explicit dummy argument. All module-level state pointers and their `bind_*_target` setup procedures disappear. The 13 legacy globals (7 in soilhydraulicsutils, 6 in surfacewaterutils) become fields on the appropriate state subrecord, populated from config at init.

The pilot that established this pattern: `docs/superpowers/specs/2026-05-13-state-init-pilot-surfacewater-design.md` (surfacewater_state — type-bound `init` consumes config slices; no transient buffers; no bind_*_target reliance). Arc 1 applies the same discipline to the utility functions called by physics modules.

## Decision

1. **Add the 13 fields** to `state%soilwater` (7 fields) and `state%surfacewater` (6 fields). Populated from `config_to_variables` at init alongside the existing `cofgen` / `fluseksatexm` writes.
2. **Refactor 5 soilhydraulicsutils functions** (`watcon`, `moiscap`, `hconduc`, `dhconduc`, `prhead`) to take `soilwater_state_t` as a new first argument. `hcomean` and `dkmean` stay unchanged — they use no module-level state.
3. **Update `moiscap`** to also take `dt` as a scalar argument (replaces `tc_dt_ptr` module pointer).
4. **Refactor surfacewaterutils functions** to read from `state%surfacewater%X` in their bodies (they already take `state` as an argument — no signature changes here, just data-source changes).
5. **Delete** module-level pointers (`cofgen`, `fluseksatexm`, `tc_dt_ptr`) and `bind_*_target` setup procedures from soilhydraulicsutils. Delete the 3 `bind_*_target` calls and 2 `bind_*_target` imports from `swap_mod.f90`.
6. **Drop `use variables` entirely** from both utility files. `WC_K_models_04_11.f90` also drops its `bind_cofgen_target` exit point (the public `bind_cofgen_target` setup procedure retires).

## Architecture

### What goes

- **Module-level pointers in `soilhydraulicsutils.f90`:** `cofgen(:,:)`, `fluseksatexm(:)`, `tc_dt_ptr` — deleted.
- **`bind_state_targets(sw_cofgen_in, sw_fluseksatexm_in)`** — deleted.
- **`bind_tc_target(tc_dt_in)`** — deleted.
- **`bind_cofgen_target(sw_cofgen_in)`** in `src/soil/WC_K_models_04_11.f90` — deleted.
- **`use variables, only: ...`** lines in both utility files — deleted (drops 7 + 6 symbol imports total).
- **`use WC_K_models_04_11, only: bind_cofgen_target`** in `swap_mod.f90:73` — deleted.
- **`use soilhydraulics_utils, only: bind_state_targets, bind_tc_target`** in `swap_mod.f90:74` — deleted.
- **Three `call bind_*_target(...)` lines** in `swap_mod.f90:111-113` — deleted.

### What stays

- The 7 utility-function algorithms in `soilhydraulicsutils` and the 4 in `surfacewaterutils`. **The compute logic is unchanged.** Only data-access paths shift from module-level globals/pointers to explicit dummy arguments.
- `state%soilwater%cofgen` and `state%soilwater%fluseksatexm` — already present. Become accessed directly via `state%soilwater%cofgen` instead of via module-level pointer.
- `config_to_variables.f90` — still writes the bare globals (other readers still consume them; this arc only retires utils). Will also write the new state fields. The bare-global writes retire in later arcs.

### State schema additions

**`src/state/soilwater_state.f90` — `soilwater_state_t` gains 7 fields:**

```fortran
   ! [GR-UTILS] Soil hydraulic property metadata — populated by config_to_variables
   ! from the soil config section. Consumed by soilhydraulicsutils functions
   ! (watcon, hconduc, moiscap, dhconduc, prhead). Static throughout the run.
   integer                       :: swsophy    = 0      !< 0 = analytical (VG/MvG), 1 = tabulated
   integer,         allocatable  :: numtab(:)            !< per-layer table entry count (size nlay)
   real(real64),    allocatable  :: sptab(:,:,:)         !< soil property table (4, numnod, max_entries)
   integer                       :: ientrytab  = 0      !< entry-table pointer
   integer,         allocatable  :: iHWCKmodel(:)        !< per-layer hydraulic-K model selector (size nlay)
   integer,         allocatable  :: layer(:)             !< per-node soil-layer index (size numnod)
   integer                       :: swfrost    = 0      !< frost-reduction simulation switch
```

(The exact rank/shape of `sptab`, `numtab`, `iHWCKmodel`, `layer` is taken VERBATIM from `variables.f90` — copy the type, kind, and rank of each declaration. The arrays are populated by `config_to_variables` at init; the rank/shape doesn't change in this arc.)

`swfrost` strictly belongs on `state%heat` or `state%timecontrol`, but it's only read by `soilhydraulicsutils` in this arc. Lands on `state%soilwater` for cohabitation; a follow-on arc can relocate when `heat/` migrates.

**`src/state/surfacewater_state.f90` — `surfacewater_state_t` gains 6 fields:**

```fortran
   ! [GR-UTILS] Surface-water runoff / discharge metadata — populated by
   ! config_to_variables from the surface_water + soil config sections.
   ! Consumed by surfacewaterutils functions (runoff, qhtab).
   real(real64),    allocatable  :: hqhtab(:)            !< Q-h table heads
   real(real64),    allocatable  :: qqhtab(:)            !< Q-h table discharges
   integer                       :: swdra      = 0      !< drainage switch
   real(real64)                  :: pondmx     = 0.0_real64 !< max ponding depth (cm)
   real(real64)                  :: rsro       = 0.0_real64 !< runoff resistance (d)
   real(real64)                  :: rsroexp    = 0.0_real64 !< runoff exponent
```

### Function signature changes (soilhydraulicsutils)

5 signatures change. The new dummy argument is the SOILWATER subrecord (not the whole `swap_state_t`):

| Function | Old signature | New signature |
|---|---|---|
| `watcon` | `(node, head)` | `(soilwater, node, head)` |
| `moiscap` | `(node, head)` | `(soilwater, dt, node, head)` |
| `hconduc` | `(node, h, theta, rfcp, tsoil_node)` | `(soilwater, node, h, theta, rfcp, tsoil_node)` |
| `dhconduc` | `(node, h, theta, dimoca, rfcp)` | `(soilwater, node, h, theta, dimoca, rfcp)` |
| `prhead` | `(node, disnod, theta, cofgen_in, h)` | `(soilwater, node, disnod, theta, h)` (the explicit `cofgen_in` argument retires — pulled from `soilwater%cofgen` instead) |

| Function | Unchanged |
|---|---|
| `hcomean` | `(swkmean, kup, klow, dzup, dzlow)` |
| `dkmean` | `(swkmean, kup, klow, dzup, dzlow)` |

The `intent(in)` of `soilwater` is sufficient — these are read-only functions. Compiler can pass by reference, no copy.

### Function body changes (soilhydraulicsutils)

Every bare reference to a migrated global becomes prefixed by `soilwater%`. Examples from `watcon`:

| Before | After |
|---|---|
| `if (swsophy == 0)` | `if (soilwater%swsophy == 0)` |
| `cofgen(1,node)` | `soilwater%cofgen(1,node)` |
| `iHWCKmodel(layer(node))` | `soilwater%iHWCKmodel(soilwater%layer(node))` |
| `numtab(layer(node))` | `soilwater%numtab(soilwater%layer(node))` |
| `sptab(:, node, j)` | `soilwater%sptab(:, node, j)` |

`moiscap` body additionally replaces `tc_dt_ptr` with the new `dt` dummy argument.

`fluseksatexm` was bound via `bind_state_targets` along with `cofgen`. It's not in our 7-field migration list because `state%soilwater%fluseksatexm` already exists. The body change is the same pattern: `fluseksatexm(node)` → `soilwater%fluseksatexm(node)`.

`prhead` currently takes `cofgen_in` as an explicit argument (one of its callers — `soilgrid.f90:388` — passes `cofgenNew`, a temporary built during grid-redistribution). **Decision: keep `cofgen_in` as an optional argument** that defaults to `soilwater%cofgen` when absent. New signature: `prhead(soilwater, node, disnod, theta, h, cofgen_in)` with `real(real64), optional, intent(in) :: cofgen_in(:,:)`. The soilgrid caller continues passing its tentative `cofgenNew`; everyone else omits the arg and gets the state-backed default. No temporary state mutation needed.

### Surfacewaterutils — different shape

`surfacewaterutils.f90`'s public functions are: `wlevst`, `swstlev`, `swstlev_from_table`, `qhtab`, `runoff`. Three of them already take state:

- `swstlev(state, wlev)` — already takes state.
- `runoff(state, ...)` — already takes state.

`wlevst`, `qhtab`, and `swstlev_from_table` — verify their existing signatures via `grep -n "function wlevst\|function qhtab\|function swstlev_from_table" src/utils/surfacewaterutils.f90 src/state/surfacewater_state.f90`. If they don't take state today, thread `state%surfacewater` in as a new first arg (same pattern as soilhydraulicsutils). If they do, just change the body.

Migration here is simpler than soilhydraulicsutils because the signatures already accommodate state for the call sites that need it:

1. **Add 6 fields to `state%surfacewater`** (schema only).
2. **Inside function bodies**, replace `hqhtab(...)` etc. with `state%surfacewater%hqhtab(...)` (or via the existing subrecord arg if already in scope).
3. **Drop `use variables, only: ...`** from the module header.
4. **`config_to_variables` writes the 6 fields** alongside its existing surfacewater writes.

No function signatures change here.

### Caller updates

**40 call sites for soilhydraulicsutils functions**, listed in 8 caller files:

- `src/soil/soilhydraulics.f90` — heaviest user (most calls)
- `src/soil/soilgrid.f90` — `prhead`
- `src/boundary/boundtop.f90`, `src/boundary/boundbottom.f90` — `watcon`, `hconduc`, `hcomean`
- `src/crop/tillage.f90`, `src/crop/irrigation.f90`, `src/crop/cropgrowth.f90` (3 sites), `src/crop/oxygenstress.f90`, `src/crop/rootextraction.f90` — `watcon`
- `src/atmosphere/meteoday.f90` — `watcon`

Each call site changes from `watcon(node, head)` → `watcon(state%soilwater, node, head)`. **State is in scope at almost every call site** (the SS-* arcs threaded state through all the physics modules). The implementer must verify per call site; if a site doesn't have state in scope, the surrounding subroutine signature is threaded with state as a precondition — same pattern used throughout the SS-DRV / SS-TCM / SS-BMI2 arcs.

For `moiscap` callers (in `soilhydraulics.f90`): add `state%timecontrol%dt` as a second new arg after `state%soilwater`.

**13 call sites for surfacewaterutils functions**, listed in 2 caller files:

- `src/drainage/surfacewater.f90` — `wlevst`, `swstlev`, `qhtab`
- `src/boundary/boundtop.f90` — `runoff`

These functions already take state in their existing signatures. Callers DO NOT change — the migration is purely internal to the utility functions.

### `swap_mod.f90` cleanup (consequence of this arc)

The bind setup at `swap_mod.f90:73-74, 111-113` (5 lines) deletes:

```fortran
! gone (line 73): use WC_K_models_04_11, only: bind_cofgen_target
! gone (line 74): use soilhydraulics_utils, only: bind_state_targets, bind_tc_target
...
! gone (line 111): call bind_cofgen_target(state%soilwater%cofgen)
! gone (line 112): call bind_state_targets(state%soilwater%cofgen, state%soilwater%fluseksatexm)
! gone (line 113): call bind_tc_target(state%timecontrol%dt)
```

Note `state%soilwater%cofgen` and `state%soilwater%fluseksatexm` remain populated by their existing init paths (S-2.12B) — only the binding step retires.

The comment `! [SS-SWC S-2.12B] bind module-level pointers...` block (3 lines) also goes.

## Verification

Per memory `feedback_per_task_regression_gate.md` + `feedback_state_schema_clean_rebuild.md`:

1. **Clean rebuild after schema changes** (`rm -rf builddir && pixi run build-linux`).
2. **`pixi run test-pfunit`** — 741/741.
3. **`pixi run -e test python tests/regression/test_output_regression.py hupselbrook surfacewater salinitystress grassgrowth`** — 4/4 byte-for-byte.
4. **End of arc:** `pixi run check-full` — 5/5 byte-for-byte.

Final greps to confirm retirement:

```bash
grep -rn "bind_state_targets\|bind_tc_target\|bind_cofgen_target" src/
# → only the cleanup tombstone comments, no live code

grep -n "use variables" src/utils/soilhydraulicsutils.f90 src/utils/surfacewaterutils.f90
# → no hits

grep -n ", pointer ::" src/utils/soilhydraulicsutils.f90
# → no module-level pointers remain
```

## Scope

| Component | Estimated lines |
|---|---|
| `soilwater_state_t` schema additions (7 fields + populate path) | ~50 |
| `surfacewater_state_t` schema additions (6 fields + populate path) | ~30 |
| `config_to_variables.f90` — extend to write the 13 new state fields | ~60 |
| `soilhydraulicsutils.f90` — 5 signature changes, body rewrites, delete pointers + bind subs | ~120 (mostly mechanical) |
| `surfacewaterutils.f90` — drop `use variables`, body rewrites | ~30 |
| `WC_K_models_04_11.f90` — delete `bind_cofgen_target` | ~10 |
| Caller updates (40 sites in 8 files) | ~50 |
| `swap_mod.f90` — delete bind imports + bind calls | ~10 |
| **Total** | **~360 lines** |

Sized like a normal SS-* arc (≈ SS-DRV Phase 1 / SS-TCM). Expected duration: 1–2 days under subagent-driven execution.

## Out of scope (deferred to later arcs)

- **All other reader clusters** (`boundary/`, `heat/`, `atmosphere/`, `soil/`, `drainage/`, `io/`, `crop/`) — separate arcs per the roadmap. After this arc the `use variables` count in `src/` drops from 33 to 31 files; further reductions come arc by arc.
- **Bare-global readers of the 13 fields** (e.g. `soilhydraulics.f90` reads `swsophy` directly from variables; it will continue doing so until its own arc migrates it to `state%soilwater%swsophy`). The bare-global declarations in `variables.f90` STAY for this arc; they retire in Arc 9.
- **`config_to_variables.f90` retirement** — Arc 9.
- **Other transient buffers** (`pondini_init_buf`, `pond_init_buf`, `h_init_buf`, `tc_*_init_buf`) — Arc 9 / the state-init pilot rollout.
- **swfrost relocation** to `state%heat` or `state%timecontrol` — deferred to the heat-cluster arc.

## Tests

- Existing tests verify byte-for-byte regression. No new tests required for this arc (the algorithmic behavior is unchanged; the migration is a pure data-access refactor).
- The existing `tests/unit/state/test_soilwater_state.pf` may need extension to verify the 7 new fields default-zero correctly — implementer decides during the work.

## References

- Globals retirement roadmap: `docs/superpowers/plans/2026-05-13-globals-retirement-roadmap.md`.
- Surfacewater state-init pilot (precedent pattern): `docs/superpowers/specs/2026-05-13-state-init-pilot-surfacewater-design.md`.
- Memory: `feedback_state_schema_clean_rebuild.md` (clean rebuild on state schema changes).
- Existing `bind_*_target` setup (will be deleted): `src/utils/soilhydraulicsutils.f90:42-55`, `src/soil/WC_K_models_04_11.f90:34-37`.
- Existing caller catalogue: see "Caller updates" section above. 40 + 13 sites.
