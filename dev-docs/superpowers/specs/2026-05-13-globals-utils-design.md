---
title: "Globals Retirement Arc 1 — utils + cofgen modernization — Design Spec"
date: 2026-05-13
status: draft
arc: GR-UTILS (globals retirement + structural reform, utils cluster)
roadmap: docs/superpowers/plans/2026-05-13-globals-retirement-roadmap.md
supersedes-portion-of: (this file replaces the earlier draft of GR-UTILS that only covered globals retirement)
---

# Globals Retirement Arc 1 — utils cluster + cofgen modernization

## Context

Arc 1 of the Globals Retirement Roadmap. Two motivations land in one arc:

**Motivation 1 — retire `use variables` from `src/utils/`.** `soilhydraulicsutils.f90` imports 7 symbols; `surfacewaterutils.f90` imports 6. Both files also carry the `bind_*_target` pattern (`bind_state_targets`, `bind_tc_target` in soilhydraulicsutils; `bind_cofgen_target` in WC_K_models_04_11.f90 — a sibling utility). The bind pattern is a coupling crutch because the underlying functions don't take state.

**Motivation 2 — modernize the `cofgen` magic-row array.** `state%soilwater%cofgen(j, node)` is a 21-row `(parameter_index, node)` matrix where every row is a *named* hydraulic parameter. Every access site has to look up "what does row 13 mean" in a comment:

```fortran
! From soilhydraulicsutils.f90:121-126
thetar = cofgen(1,node)        ! row 1 = residual water content
thetas = cofgen(2,node)        ! row 2 = saturated water content
alfamg = cofgen(4,node)        ! row 4 = MvG alpha
n      = cofgen(6,node)        ! row 6 = MvG n
m      = cofgen(7,node)        ! row 7 = MvG m
h_enpr = cofgen(9,node)        ! row 9 = air-entry pressure
```

Replacing `cofgen(:,:)` with `vg_params(:)` — an array of typed records with named fields — makes every access site self-documenting and removes the magic indices.

**Why combine the two motivations:** the cofgen reform and the utils retirement touch *exactly the same files and exactly the same call sites*. Splitting them into two arcs would force the utils retirement to use either (a) the verbose explicit-fields signature pattern that gets thrown away when cofgen modernizes, or (b) the pass-subrecord pattern that obscures what each function reads. Combining them gives the utility functions clean, self-documenting signatures from day one with no redundant churn.

## Decision

1. **Define a new typed record `vanGenuchten_params_t`** with named fields replacing the 21-row `cofgen` schema. Lives in a new small module `src/state/hydraulic_params_mod.f90` (kept separate from `soilwater_state_mod` to avoid the module's growing further).
2. **Migrate `state%soilwater%cofgen(:,:)` → `state%soilwater%vg_params(:)`** — an allocatable array of `vanGenuchten_params_t`, one entry per node.
3. **Refactor 5 soilhydraulicsutils functions** + `functionvalue_04_11` (in `WC_K_models_04_11.f90`) to take the relevant `vg_params(node)` slice as an explicit argument. Other static metadata (`iHWCKmodel`, `layer`, `swsophy`, etc.) becomes additional explicit dummy args.
4. **Refactor surfacewaterutils functions** to read from `state%surfacewater%X` directly (they already take state as an argument — no signature changes needed, just data-source changes).
5. **Add the remaining utils-imported globals** (`numtab`, `sptab`, `ientrytab`, `iHWCKmodel`, `layer`, `swsophy`, `swfrost` for soilhydraulicsutils; `hqhtab`, `qqhtab`, `swdra`, `pondmx`, `rsro`, `rsroexp` for surfacewaterutils) as fields on `state%soilwater` / `state%surfacewater`. Populated from `config_to_variables` at init.
6. **Delete** all module-level pointers (`cofgen`, `fluseksatexm`, `tc_dt_ptr`) and `bind_*_target` setup procedures from soilhydraulicsutils and WC_K_models_04_11. **Delete** the corresponding calls + imports from `swap_mod.f90`.
7. **Drop `use variables`** entirely from soilhydraulicsutils, WC_K_models_04_11, and surfacewaterutils.

**Explicit non-scope:** `sptab(:,:,:)` (the tabulated-property table) has the same magic-row pattern but is consumed by different code paths (`soilhydraulicsutils` only via the `swsophy=1` branch; `sptabulated.f90`'s separate spline math). It is migrated in a later arc (probably the `soil/` arc that touches `soilhydraulics.f90` and `sptabulated.f90`). For this arc, `sptab` moves into `state%soilwater` as a 3D array unchanged in shape; its modernization to a derived type is deferred.

## Architecture

### What goes

- **`state%soilwater%cofgen(:,:)`** — the legacy `(21, numnod)` matrix. Deleted from `soilwater_state_t`.
- **Module-level pointers** in `src/utils/soilhydraulicsutils.f90`: `cofgen(:,:)`, `fluseksatexm(:)`, `tc_dt_ptr`. Deleted.
- **Module-level pointer** in `src/soil/WC_K_models_04_11.f90`: `cofgen(:,:)`. Deleted.
- **`bind_state_targets`** (in soilhydraulicsutils) — deleted.
- **`bind_tc_target`** (in soilhydraulicsutils) — deleted.
- **`bind_cofgen_target`** (in WC_K_models_04_11) — deleted.
- **`use variables, only: ...`** lines in soilhydraulicsutils, WC_K_models_04_11, and surfacewaterutils. Deleted (drops 7 + 4 + 6 = 17 symbol imports total).
- **5 `bind_*_target` calls + 2 `use ..., only: bind_*` imports** in `swap_mod.f90:73-74, 111-113`. Deleted.

### What stays

- **The numerical algorithms** in all 7 affected functions (`watcon`, `moiscap`, `hconduc`, `dhconduc`, `prhead`, `functionvalue_04_11` and its helpers). Unchanged. Only data-access paths shift.
- **`state%soilwater%fluseksatexm`** — already on state, stays. Now read directly as a dummy arg into functions that need it.
- **`config_to_variables.f90`** — extended to populate the new `vg_params` array and the other migrated state fields, alongside its existing writes. Retirement of the file itself is Arc 9.
- **The bare globals in `variables.f90`** — still alive (other reader clusters not yet migrated). Their retirement is Arc 9.

### New type: `vanGenuchten_params_t`

New module: `src/state/hydraulic_params_mod.f90`

```fortran
!> @file hydraulic_params_mod.f90
!! SS-GR-UTILS: typed van-Genuchten / MvG / PDI hydraulic parameter
!! record. Replaces the legacy `cofgen(j, node)` magic-row matrix.
!! One instance per soil node; aggregated as `state%soilwater%vg_params(:)`.
!!
!! Row-name correspondence to legacy cofgen indices:
!!   cofgen(1, n)  -> vg_params(n)%thetar       (residual water content)
!!   cofgen(2, n)  -> vg_params(n)%thetas       (saturated water content)
!!   cofgen(3, n)  -> vg_params(n)%ksat
!!   cofgen(4, n)  -> vg_params(n)%alpha
!!   cofgen(5, n)  -> vg_params(n)%lpar         (pore-size index L)
!!   cofgen(6, n)  -> vg_params(n)%npar
!!   cofgen(7, n)  -> vg_params(n)%mpar
!!   cofgen(9, n)  -> vg_params(n)%h_enpr       (air-entry pressure modification)
!!   cofgen(13, n) -> vg_params(n)%alpha_2      (bi-modal 2nd alpha)
!!   cofgen(14, n) -> vg_params(n)%npar_2       (bi-modal 2nd n)
!!   cofgen(15, n) -> vg_params(n)%mpar_2       (bi-modal 2nd m)
!!   cofgen(16, n) -> vg_params(n)%omega_1
!!   cofgen(17, n) -> vg_params(n)%omega_2
!!   cofgen(18, n) -> vg_params(n)%h0
!!   cofgen(19, n) -> vg_params(n)%ha
!!   cofgen(20, n) -> vg_params(n)%apar
!!   cofgen(21, n) -> vg_params(n)%omega_k
!!
!! Rows 8, 10, 11, 12 of legacy cofgen are unused in the codebase
!! (verified via `grep -n "cofgen( ?\(8\|10\|11\|12\)" src/`). The
!! corresponding fields are omitted; if a future iHWCKmodel branch
!! needs them, add named fields then.
module hydraulic_params_mod
   use iso_fortran_env, only: real64
   implicit none
   private
   public :: vanGenuchten_params_t

   type :: vanGenuchten_params_t
      ! Universal MvG (rows 1-7 in legacy cofgen)
      real(real64) :: thetar  = 0.0_real64  !< residual water content (m³/m³)
      real(real64) :: thetas  = 0.0_real64  !< saturated water content (m³/m³)
      real(real64) :: ksat    = 0.0_real64  !< saturated hydraulic K (cm/d)
      real(real64) :: alpha   = 0.0_real64  !< MvG α (1/cm)
      real(real64) :: lpar    = 0.0_real64  !< pore-size index L
      real(real64) :: npar    = 0.0_real64  !< MvG n
      real(real64) :: mpar    = 0.0_real64  !< MvG m (typically 1 - 1/n)
      ! MvG with air-entry modification (row 9)
      real(real64) :: h_enpr  = 0.0_real64  !< air-entry pressure (cm)
      ! Bi-modal MvG (rows 13-17 — used by iHWCKmodel ∈ {3, 6, 7, 10, 11})
      real(real64) :: alpha_2 = 0.0_real64  !< 2nd α
      real(real64) :: npar_2  = 0.0_real64  !< 2nd n
      real(real64) :: mpar_2  = 0.0_real64  !< 2nd m
      real(real64) :: omega_1 = 0.0_real64  !< 1st mode weight
      real(real64) :: omega_2 = 0.0_real64  !< 2nd mode weight
      ! PDI parameters (rows 18-21 — used by iHWCKmodel ∈ {5, 7, 8..11})
      real(real64) :: h0      = 0.0_real64  !< PDI h0
      real(real64) :: ha      = 0.0_real64  !< PDI ha
      real(real64) :: apar    = 0.0_real64  !< PDI A parameter
      real(real64) :: omega_k = 0.0_real64  !< PDI Ω_K
   end type vanGenuchten_params_t

end module hydraulic_params_mod
```

16 named fields total. Implementer verifies during execution that no cofgen accesses to unmapped row indices exist (the `grep` cited in the doc comment).

### State schema changes

**`soilwater_state_t`** in `src/state/soilwater_state.f90`:

- **Delete:** `real(real64), allocatable :: cofgen(:,:)`
- **Add:** `type(vanGenuchten_params_t), allocatable :: vg_params(:)` (sized `numnod`)
- **Add (utils-imported globals migrated):**
  ```fortran
  integer                       :: swsophy    = 0
  integer,         allocatable  :: numtab(:)
  real(real64),    allocatable  :: sptab(:,:,:)    ! shape preserved verbatim from variables.f90
  integer                       :: ientrytab  = 0
  integer,         allocatable  :: iHWCKmodel(:)
  integer,         allocatable  :: layer(:)
  integer                       :: swfrost    = 0
  logical,         allocatable  :: BiModal(:)      ! used by WC_K_models_04_11
  logical,         allocatable  :: NoVap(:)        ! used by WC_K_models_04_11
  ```

`fluseksatexm` is already on `state%soilwater` — unchanged.

**`surfacewater_state_t`** in `src/state/surfacewater_state.f90`:

- **Add (utils-imported globals migrated):**
  ```fortran
  real(real64),    allocatable  :: hqhtab(:)
  real(real64),    allocatable  :: qqhtab(:)
  integer                       :: swdra      = 0
  real(real64)                  :: pondmx     = 0.0_real64
  real(real64)                  :: rsro       = 0.0_real64
  real(real64)                  :: rsroexp    = 0.0_real64
  ```

(`swdra` is a high-level switch; reasonable to live here since it gates surfacewater activation. Could also live on `state%timecontrol` or `state%simulation` later — out of scope for this arc.)

### `config_to_variables.f90` updates

After this arc, `config_to_variables`:
- **Writes the legacy `cofgen` matrix** as a transient step (so other readers still consuming `cofgen` keep working).
- **Also writes `state%soilwater%vg_params(:)`** with the same values, indexed by row meaning. Example: `state%soilwater%vg_params(node)%thetar = config%soil%cofgen_residual_water_content(layer)` (or whatever the config-side path is).
- **Also writes** the additional migrated state fields (swsophy, numtab, sptab, ientrytab, iHWCKmodel, layer, swfrost, BiModal, NoVap on state%soilwater; hqhtab, qqhtab, swdra, pondmx, rsro, rsroexp on state%surfacewater).

Both `cofgen` AND `vg_params` are populated during the transition. The `cofgen` write retires when the LAST non-utils reader migrates (Arc 5 — soil cluster). Tracked as a deferred item in the roadmap.

### Function signatures — soilhydraulicsutils

Five signatures change to take explicit, self-documenting arguments. The numerical bodies are essentially unchanged — only the source of each value shifts from "magic row index" or "module-level pointer" to "named dummy argument."

| Function | Old signature | New signature |
|---|---|---|
| `watcon` | `(node, head)` | `(head, vg, model)` — `node` no longer needed since `vg` is the slice for that node |
| `moiscap` | `(node, head)` | `(head, vg, model, dt)` |
| `hconduc` | `(node, h, theta, rfcp, tsoil_node)` | `(h, theta, rfcp, tsoil_node, vg, model, use_ksatexm)` |
| `dhconduc` | `(node, h, theta, dimoca, rfcp)` | `(h, theta, dimoca, rfcp, vg, model)` |
| `prhead` | `(node, disnod, theta, cofgen_in, h)` | `(disnod, theta, h, vg, model)` — `cofgen_in` arg retires; for the soilgrid-redistribution callsite, the caller passes a tentative `vg_params_new(NodeNew(node,2))` instead of the old `cofgenNew(:,node)` |

Where:
- `vg` is `type(vanGenuchten_params_t), intent(in)` — the slice for this node.
- `model` is `integer, intent(in)` — the `iHWCKmodel(layer(node))` value, computed by the caller.
- `dt` is `real(real64), intent(in)` — replaces `tc_dt_ptr`.
- `use_ksatexm` is `logical, intent(in)` — replaces the `fluseksatexm(node)` lookup.

Functions that don't read any module-level state stay unchanged: `hcomean(swkmean, kup, klow, dzup, dzlow)`, `dkmean(swkmean, kup, klow, dzup, dzlow)`. These are pure mean-computation utilities.

### Function bodies — soilhydraulicsutils

Every bare-name access to a migrated parameter becomes the dummy-arg field access. Example diff for `watcon`:

**Before:**
```fortran
function watcon(node, head)
   ! ... declarations ...
   if (swsophy == 0) then
      thetar = cofgen(1,node)
      thetas = cofgen(2,node)
      alfamg = cofgen(4,node)
      n      = cofgen(6,node)
      m      = cofgen(7,node)
      h_enpr = cofgen(9,node)
      if (iHWCKmodel(layer(node)) == 2) then
         watcon = dmax1(1.0000001_real64*thetar, thetar + (thetas-thetar)*dexp(alfamg*head))
      else if (iHWCKmodel(layer(node)) == 3) then
         ! ...
      end if
   end if
```

**After (this arc):**
```fortran
function watcon(head, vg, model) result(theta)
   real(real64),                 intent(in) :: head
   type(vanGenuchten_params_t),  intent(in) :: vg
   integer,                      intent(in) :: model
   real(real64) :: theta
   ! ... local declarations (swsophy moved to a dummy arg of its own
   !     if/when it's needed — for the analytical path, model is the
   !     dispatcher; swsophy is handled at the caller via dispatching
   !     to watcon or to a tabulated lookup function). See "swsophy
   !     dispatch" below.
   if (model == 2) then
      theta = dmax1(1.0000001_real64*vg%thetar, &
                    vg%thetar + (vg%thetas - vg%thetar) * dexp(vg%alpha*head))
   else if (model == 3) then
      ! ... bi-modal branch using vg%alpha_2, vg%npar_2, etc.
   end if
```

The body is essentially the same — `vg%thetar` / `vg%alpha` replace `cofgen(1,node)` / `cofgen(4,node)`. The intent is self-documenting at every site.

#### swsophy dispatch

Currently `watcon` has an outer `if (swsophy == 0) then ... else if (swsophy == 1) then ...` switch — analytical vs. tabulated. The tabulated branch uses `sptab` and `numtab`. **Decision:** keep the `swsophy` dispatch inside `watcon`. The function signature gains a `tabulated_data` argument that bundles `(swsophy, numtab, sptab)` — or the caller passes them as additional dummy args. **For this arc, `swsophy` is sufficient as an additional integer dummy arg**, with `numtab` and `sptab` accessible from a passed `state%soilwater` reference *only on the tabulated branch*. The tabulated branch is the messier corner of the migration; implementer may decide during execution whether to split `watcon` into `watcon_analytical(head, vg, model)` and `watcon_tabulated(head, soilwater, node)` for full clarity.

Per the YAGNI principle: do the cleanest analytical-side migration first; revisit the tabulated dispatch when it actively blocks something.

### Function signature — `functionvalue_04_11`

In `src/soil/WC_K_models_04_11.f90`, the function currently uses module-level pointer `cofgen` and bare globals `iHWCKmodel`, `BiModal`, `NoVap`, `layer`.

**Old:**
```fortran
function functionvalue_04_11(iType, iNode, h, wc, temp)
```

**New:**
```fortran
function functionvalue_04_11(iType, h, vg, model, is_bimodal, no_vap, wc, temp) result(val)
   integer,                     intent(in)           :: iType
   real(real64),                intent(in)           :: h
   type(vanGenuchten_params_t), intent(in)           :: vg
   integer,                     intent(in)           :: model     ! iHWCKmodel for this layer
   logical,                     intent(in)           :: is_bimodal
   logical,                     intent(in)           :: no_vap
   real(real64),                intent(in), optional :: wc, temp
   real(real64) :: val
```

`iNode` retires from the signature — callers pre-slice. The function body assigns `WCr = vg%thetar`, `WCs = vg%thetas`, etc., instead of `WCr = cofgen(1,iNode)`. The case-by-case branches (`case (5): h0 = vg%h0`, etc.) follow the same pattern.

### `surfacewaterutils.f90` migration

Three of `surfacewaterutils`'s public functions already take state (`swstlev(state, wlev)`, `runoff(state, ...)`). The migration here is simpler:

1. **No signature changes.** The functions already have state access.
2. **Function bodies** replace `hqhtab(...)` / `qqhtab(...)` / `swdra` / `pondmx` / `rsro` / `rsroexp` reads with `state%surfacewater%X` access via the existing dummy.
3. **Drop `use variables`** from the module header.
4. **`config_to_variables` writes the 6 fields** into `state%surfacewater` at init.

`wlevst`, `qhtab`, and `swstlev_from_table` — implementer checks whether each currently takes state. If not, the surrounding caller's state is in scope and a state arg gets threaded the same way.

### Caller updates — soilhydraulicsutils

40 call sites across 8 files (listed in [Arc 1 spec inventory grep, run before drafting](#references)). At each:

- Caller computes per-call slices: `vg = state%soilwater%vg_params(node)`, `model = state%soilwater%iHWCKmodel(state%soilwater%layer(node))`.
- Caller passes the slices to the function.

**Example transformation (the hottest call site, `soilhydraulics.f90:151`):**

Before:
```fortran
sw_k(i) = hconduc(i, sw_h(i), sw_theta(i), state%heat%rfcp(i), state%heat%tsoil(i))
```

After:
```fortran
sw_k(i) = hconduc(sw_h(i), sw_theta(i), state%heat%rfcp(i), state%heat%tsoil(i), &
                  state%soilwater%vg_params(i), &
                  state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                  state%soilwater%fluseksatexm(i))
```

3 extra arg lines. Verbose, but each arg is self-explanatory.

**Where it might be worth introducing a local alias** to clean up repeated boilerplate (e.g. inside a tight node-loop):

```fortran
associate( vg => state%soilwater%vg_params, &
           model => state%soilwater%iHWCKmodel, &
           layer => state%soilwater%layer, &
           kexm => state%soilwater%fluseksatexm )
   do i = 1, numnod
      sw_k(i) = hconduc(sw_h(i), sw_theta(i), state%heat%rfcp(i), state%heat%tsoil(i), &
                        vg(i), model(layer(i)), kexm(i))
   end do
end associate
```

The implementer chooses per-site whether the associate-block aliasing is worth it; either form is acceptable.

### `swap_mod.f90` cleanup (consequence)

The 5 lines at `swap_mod.f90:73-74, 111-113` delete:

```fortran
! gone (line 73): use WC_K_models_04_11, only: bind_cofgen_target
! gone (line 74): use soilhydraulics_utils, only: bind_state_targets, bind_tc_target
...
! gone (line 111): call bind_cofgen_target(state%soilwater%cofgen)
! gone (line 112): call bind_state_targets(state%soilwater%cofgen, state%soilwater%fluseksatexm)
! gone (line 113): call bind_tc_target(state%timecontrol%dt)
```

After this arc, `swap_mod`'s strangler-fig leftover for soilwater binding is GONE entirely.

## Verification

Per memory `feedback_per_task_regression_gate.md` + `feedback_state_schema_clean_rebuild.md`:

- **Clean rebuild** required (state schema changes): `rm -rf builddir && pixi run build-linux`.
- **pFUnit:** 741/741 expected.
- **Regression byte-for-byte:** 4/4 (then 5/5 at end of arc).

Because the bit-level data is unchanged (a `real(real64)` field at offset `k` in `vanGenuchten_params_t` holds the same bits as `cofgen(k, node)` would have), regression must pass byte-for-byte. **Any deviation = a bug in the migration.**

Final greps to confirm retirement:

```bash
grep -rn "bind_state_targets\|bind_tc_target\|bind_cofgen_target" src/
# → no live code; only deletion-tombstone comments allowed

grep -n "use variables" src/utils/soilhydraulicsutils.f90 \
                       src/utils/surfacewaterutils.f90 \
                       src/soil/WC_K_models_04_11.f90
# → no hits

grep -n ", pointer ::" src/utils/soilhydraulicsutils.f90 src/soil/WC_K_models_04_11.f90
# → no module-level pointers

grep -rn "state%soilwater%cofgen\b" src/
# → either none, or only the temporary writer in config_to_variables
#   (kept for non-utils readers, retires in Arc 5 / Arc 9)
```

## Scope

| Component | Estimated lines |
|---|---|
| New `hydraulic_params_mod.f90` module | ~50 |
| `soilwater_state_t` — add `vg_params` + 9 utils-related fields; delete `cofgen` | ~50 |
| `surfacewater_state_t` — add 6 fields | ~25 |
| `config_to_variables.f90` — populate vg_params + 9+6 = 15 new state fields; keep legacy cofgen write transitional | ~200 |
| `soilhydraulicsutils.f90` — 5 signature changes, body rewrites, delete pointers + bind subs | ~250 |
| `WC_K_models_04_11.f90` — signature change for `functionvalue_04_11`, body rewrite, delete pointer + bind | ~150 |
| `surfacewaterutils.f90` — drop `use variables`, body rewrites | ~30 |
| `soil/soilhydraulics.f90` — update inline `sw%cofgen(j,node) = ...` writes to write vg_params fields too; update internal cofgen reads | ~50 |
| `soil/soilgrid.f90` — update cofgenNew redistribution to redistribute vg_params | ~30 |
| Caller updates (40 sites in 8 files for utils functions; ~10 sites for surfacewaterutils) | ~150 |
| `swap_mod.f90` — delete bind imports + bind calls | ~10 |
| **Total** | **~995 lines** |

Sized as a **medium-large arc** — between SS-DRV Phase 1 (~600) and SS-BMI2 (~960). Estimated duration: **3–5 days under subagent-driven execution**.

## Out of scope (deferred to later arcs)

- **`sptab(:,:,:)` modernization.** The tabulated soil-property table has the same magic-row pattern but consumes `swsophy=1` branch and `sptabulated.f90`'s spline math. Deferred to the soil cluster arc (Arc 5).
- **`config_to_variables.f90` retirement.** Arc 9. This arc keeps both `cofgen` and `vg_params` populated; the `cofgen` writes retire as readers migrate.
- **Bare-global readers of cofgen / sptab outside utils.** `soilhydraulics.f90` reads `cofgen` inline. These migrate in Arc 5.
- **Other clusters** — boundary, heat, atmosphere, soil, drainage, io, crop — separate arcs per roadmap.
- **`swfrost` relocation to `state%heat` or `state%timecontrol`.** Lands on `state%soilwater` here for cohabitation with hydraulic data; future arc relocates.

## Tests

Existing pFUnit + regression suites cover the migration:
- Byte-for-byte regression on `check-full` is the strong correctness signal.
- `tests/unit/state/test_soilwater_state.pf` should be extended with assertions verifying `vg_params(:)` allocates correctly post-init and that a known vg_params field round-trips through `state%soilwater%vg_params(node)%thetar`. Implementer adds during execution.

A bonus test: `tests/unit/state/test_hydraulic_params.pf` (new) verifies the new derived type defaults are all zero, allocation works, and a populated value round-trips correctly. ~30 lines, mechanical.

## References

- Globals retirement roadmap: `docs/superpowers/plans/2026-05-13-globals-retirement-roadmap.md`.
- Surfacewater state-init pilot (precedent for state%X%init pattern): `docs/superpowers/specs/2026-05-13-state-init-pilot-surfacewater-design.md`.
- Memory: `feedback_state_schema_clean_rebuild.md` (clean rebuild on state schema changes).
- Cofgen row catalogue (source for the vg_params field list): `src/utils/soilhydraulicsutils.f90:111-126` (analytical-MvG path) + `src/soil/WC_K_models_04_11.f90:54-103` (extended models).
- Existing `bind_*_target` setup (will be deleted): `src/utils/soilhydraulicsutils.f90:42-55`, `src/soil/WC_K_models_04_11.f90:34-37`.
- 40 utils caller sites: confirmed via `grep -rn "use soilhydraulics_utils" src/ --include="*.f90"` — see inventory in earlier brainstorm.
