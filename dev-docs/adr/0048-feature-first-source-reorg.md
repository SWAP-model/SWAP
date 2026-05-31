# ADR 0048 — Feature-first source-tree reorganization

**Status:** Accepted (2026-05-31)
**Plan:** `dev-docs/superpowers/plans/2026-05-31-feature-first-reorg.md` (executed inline, 6 commits).

## Context

The strangler-fig migration is complete (TOML → typed config → typed state, no
bare-globals layer; ADRs 0030–0047). What remained was *clarity*: the directory
tree had grown organically across ~20 migration arcs and no longer told a clean
story. A structural review (2026-05-31) found four concrete smells, each backed
by the inter-folder `use`-dependency graph:

1. **`src/core/` is apex + floor + facades at once.** It holds the leaves every
   other folder depends on (`constants`, `arrays`, `dtutil`, `swap_log`) *and*
   the top-level orchestrator (`swap_mod`, `swap_main`) *and* the three C-ABI
   facades (`swap_capi/bmi/xmi`) *and* a whole compute subsystem
   (`timecontrol_mod`). `core` both depends on everything and is depended on by
   everything — there is no readable "bottom" to the tree.

2. **`src/utils/` mixes true leaves with domain compute.** `arrayutils` and
   `numericalsolvers` are pure leaf math, but `soilhydraulicsutils` and
   `surfacewaterutils` operate on `soilwater_state_t` / `surfacewater_state_t` —
   domain logic mislabeled as utilities (the `utils → state` ×12 edge is a
   layering inversion).

3. **`src/crop/` is four subsystems in one folder.** Crop growth proper
   (dispatcher + fixed/grass/wofost models), plus irrigation, tillage, and the
   entire WOFOST soil-nutrient cycling (`wofost_soil_*`, `wofostnut`,
   `management_soil`) — none of which is crop *biology*.

4. **`src/boundary/` (447 L) is soil-water boundary conditions**, not its own
   domain (the `soil → boundary` edge confirms it).

The forward intent — stateless kernels, no hidden state, I/O at the edges,
in-process multi-instance + BMI/Python binding — is a *data-oriented* statement:
explicit shared state threaded through pure transforms. The tree should make
that legible.

## Decision

Reorganize on a single principle: **feature-first for behavior, layer-first for
the shared data model, contracts, and edges.**

### Target structure

```
src/
  core/          foundation leaves ONLY (flat): constants, arrays(dims),
                 dtutil, swap_log + numericalsolvers, arrayutils (from utils/)
  error/  validation/     standalone leaves (unchanged)
  config/  state/  io/     the shared data model, input contract, and edges

  atmosphere/
  soilwater/     = old soil/ + boundary/ + utils/soilhydraulicsutils
  heat/  solute/
  drainage/      + utils/surfacewaterutils
  timecontrol/   pulled out of core/ (pairs with state/timecontrol_state)
  crop/
    cropgrowth, cropgrowth_helpers, rootextraction, oxygenstress,
    irrigation, tillage                  (irrigation + tillage kept here for now)
    fixed/    cropfixed_*
    grass/    cropgrass_*
    wofost/   cropwofost_*, wofostnut, wofost_soil_*×9, management_soil

  driver/        swap_mod (lifecycle), swap_main (exe), swap_ensemble_mod
  bindings/      swap_capi, swap_bmi, swap_xmi, bmi_constants
```

### Why `state/` and `config/` stay as layers (not folded into features)

A deliberate exception to feature-first. State and config are the
*depended-upon foundation*, not feature-private:

- The dependency graph shows every compute folder pointing at them
  (`crop → state` ×25, `soil → state` ×21, `state → config` ×34); they sit
  *below* the features.
- `swap_state_t` / `swap_config_t` aggregate **all** subsystem records into one
  type — the data model is inherently cross-cutting. At runtime the whole
  `swap_state_t` is threaded through `swap_run_step` and kernels read sibling
  state (soilwater reads `heat%tsoil`, crop reads soil, …).
- The BMI/multi-instance goal wants one coherent "complete runtime state" to
  expose/serialize and one "complete TOML schema" — not a model fragmented
  across a dozen folders.
- It keeps the project's "clean rebuild after any `src/state/*` or `src/config/*`
  change" rule (the `.mod` cross-boundary hazard) a one-glance answer.

This is the standard functional-core / data-oriented shape: data definitions
separate from the transforms over them.

### Why WOFOST nutrients live under `crop/wofost/`

All nutrient dynamics (`wofost_soil_*`, `wofostnut`, `management_soil`) are part
of the WOFOST model and move with it, rather than spinning out a separate
`nutrients/` compute folder. The pre-existing `nutrients_state` / `nutrients_config`
/ `nutrients_csv` *layer* records keep their layer homes (per the rule above).

### Mechanics

The move is **`git mv` + build-path edits only** — no `use` statements change,
because Fortran module names are independent of file paths and ninja computes
Fortran module dependencies automatically (intra-target source order is
irrelevant). Two build files are edited in lockstep: the root `meson.build`
(source lists + `src_inc` include dirs) and `tests/unit/meson.build` (134
explicit `../../src/...` paths). Each commit moves one cohesive group and ends
`check-fast` green and byte-identical.

## Consequences

+ `core/` becomes purely the foundation floor; orchestration (`driver/`) and the
  C-ABI facades (`bindings/`) are separated. The "facades cannot co-link"
  constraint (currently a long `meson.build` comment) becomes legible from the
  tree.
+ The three `*utils → state` layering inversions disappear (`utils/` dissolved).
+ `crop/wofost/` is the single self-contained home for the whole WOFOST model.
+ `crop/` stops being four subsystems in a trenchcoat; the three growth models
  split into `fixed/ grass/ wofost/`.
+ `state/` + `config/` remain the coherent, complete data model and input schema
  the BMI binding will expose.
- More top-level directories. Accepted: each now has one clear responsibility.

## Out of scope (separate follow-ups)

- **`state → io` init inversion** — `state%X%init` still triggers file reads;
  pushing those to the io layer is behavioral, deferred.
- **`wofost_soil_interface.f90`'s implicit-`SAVE` exchange state** — moving it
  into `crop/wofost/` does not fix that it is hidden cross-instance state; a real
  multi-instance blocker (sibling to `crop_config_global`), flagged on its own.
- **`irrigation` / `tillage`** kept in `crop/` for now; promoting them to
  top-level features is a later call.

## Verification

Executed inline in 6 commits, each byte-identical at the commit boundary:

1. `core/` → `core/` + `driver/` + `bindings/`
2. `soil/` → `soilwater/` (+ `boundary/` folded in)
3. `utils/` dissolved (→ `core/`, `soilwater/`, `drainage/`)
4. `timecontrol/` extracted from `core/`
5. `crop/` → `fixed/` `grass/` `wofost/` (WOFOST nutrients under `wofost/`)
6. docs + full-suite gate

Final state:
- `check-fast` 4/4 byte-identical (held at every commit).
- `check-full` 5/5 non-xfail byte-identical; the two pre-existing xfails
  (`winter`, `soilhysteresis`) unchanged.
- pFUnit **810 tests** OK — unchanged from baseline at every commit (no test
  ran "dark").
- `swap` executable + `libswap_bmi.so` + `libswap_xmi.so` all link.

No `.f90` content changed — pure `git mv` + build-path edits (root `meson.build`
source lists + `src_inc`; `tests/unit/meson.build` paths). Module names are
path-independent, so zero `use`-statement edits were needed.
