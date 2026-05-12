# Subsystem Migration Discovery: Tillage

**Date:** 2026-05-12
**Status:** discovery (read-only inventory)
**Migration #:** 9 of N — small finishing arc following soil-water 4-arc decomposition
**Branch:** `development`
**Scope note:** This arc carves the 31 `till_*` globals in `variables.f90` out into a new
`tillage_state_t` and adds `state%tillage` to `swap_state_t`. DoTillage/Change_MvGpars/
Adapt_WC_H already take `state` (SS-ATM A-2.6, SS-SWC S-1.9, S-2.6). Phase 0 is zero —
tillage is already 100% TOML (ADR 0021).
**Predecessor docs:**
- `docs/superpowers/specs/state-migration-playbook.md` (26 lessons)
- `docs/superpowers/specs/2026-05-10-state-migration-boundary-discovery.md` (style template)
- `docs/adr/0035-state-migration-boundary.md`
- `docs/adr/0036-state-migration-crop-uptake.md`
- `docs/adr/0038-state-migration-soilwater-core.md`

> Read-only discovery: no code changes made. All file:line references are anchors for
> the design phase.

---

## 1. Big picture

### Subsystem role

`tillage.f90` models tillage-event-driven changes to soil bulk density and van Genuchten
hydraulic parameters. On each day a scheduled tillage event fires, it:
1. Updates `Bdens` (per-horizon bulk density) via `Change_Bdens`.
2. Recalculates `ParamVG` (per-layer VG parameters) via `Change_MvGpars` — which also
   writes `state%soilwater%cofgen` (S-2.6 cutover complete).
3. Redistributes water/head to account for pore-space change via `Adapt_WC_H` — which
   reads and writes `state%soilwater%theta`, `state%soilwater%h`, `state%soilwater%pond`
   (S-2.6 cutover complete).
4. On non-event days, `Consolidate_Bdens` drives exponential re-consolidation using
   `state%atmosphere%nraida` (A-2.6 cutover complete).

### Lines of code

```
446  src/crop/tillage.f90
```

Single home file. Smallest home-file LoC in the arc series (boundary: 494 across 2 files).

### Subroutine structure

| Subroutine | Lines | Role |
|---|---|---|
| `DoTillage(iTask, state)` | 40–196 | Dispatcher — iTask=1 init, 2 rate/event, 3 output, 4 closure |
| `Change_MvGpars(state)` | 199–254 | Recalculate ParamVG + cofgen on Bdens change |
| `Adapt_WC_H(TEST, state)` | 257–372 | Redistribute theta/h/pond after pore-space change |
| `Consolidate_Bdens(state)` | 376–386 | Exponential re-consolidation each non-event step |
| `Change_Bdens` | 389–402 | Apply tillage-event intensity to Bdens (no state) |
| `set_iTill` | 405–414 | Find entry point in event table from simulation start |
| `det_MNSH` | 417–427 | Determine MaxNumSoilHo/MaxNumSoilCP at Max_Z_tillage |
| `Change_Tillage_Info(iTill)` | 429–444 | Copy type-table row into per-layer working arrays |

### Call sites in swap.f90

```
swap.f90:194  call soilwater_init(state%soilwater, numnod, numlay)   ! SS-CRP C-1.2
swap.f90:227  if (flTillage) call DoTillage(1, state)                ! iTask=1 init
swap.f90:328  if (flTillage) call DoTillage(2, state)                ! iTask=2 rate
swap.f90:437  if (flTillage) call DoTillage(3, state)                ! iTask=3 output
```

`DoTillage(1)` fires at line 227 — **after** `soilwater_init` at line 194. The order
concern from S-1.9 is confirmed resolved. `DoTillage(4)` (closure) has a `continue` body
and is never explicitly called from swap.f90.

### Signature status

All three routines that matter already accept `state`:
- `DoTillage(iTask, state)` — `state` is `intent(inout)` (A-2.6 + S-1.9)
- `Change_MvGpars(state)` — `state` is `intent(inout)` (S-1.9; reads/writes `state%soilwater%cofgen`)
- `Adapt_WC_H(TEST, state)` — `state` is `intent(inout)` (S-1.9; reads/writes `state%soilwater%{theta,h,pond}`)
- `Consolidate_Bdens(state)` — `state` is `intent(in)` (A-2.6; reads `state%atmosphere%nraida`)

No signature changes needed in Phase 1 of this arc. The new `state%tillage` component will
be threaded through by extending the existing `state` parameter.

---

## 2. Owned-globals inventory (31 till_* fields)

All 31 `till_*` fields are declared in `variables.f90:593–623`. They divide into three
natural groups:

### Group A — Config parameters (scalars, written once by config_to_variables)

| Variable | Type | Shape | Cadence | Notes |
|---|---|---|---|---|
| `till_swtill` | integer | scalar | config | switch 0/1; also sets `flTillage` (ADR 0020) |
| `till_i_n_model` | integer | scalar | config | n-parameter treatment model (1-3) |
| `till_iRedist` | integer | scalar | config | redistribution type after MvG change |
| `till_Max_Z_tillage` | real(8) | scalar | config | max tillage depth (cm) |
| `till_Ntill` | integer | scalar | config | number of tabulated tillage events |
| `till_Ntypes` | integer | scalar | config | number of tillage types |

These six are read from TOML by `config_to_variables.f90` and then never mutated. They are
closest to config-constant parameters rather than runtime state. Candidate for exclusion
(stay as config-only or move to typed config) — see Section 8 open question.

### Group B — Event-table arrays (allocated + populated once from config)

| Variable | Type | Shape | Cadence | Notes |
|---|---|---|---|---|
| `till_Date_tillage` | real(8)(:) | per-event allocatable | config | event dates (t1900) |
| `till_Z_tillage` | real(8)(:) | per-event allocatable | config | tillage depths |
| `till_I_tillage` | real(8)(:) | per-event allocatable | config | intensities (0-1) |
| `till_Type_Tillage` | integer(:) | per-event allocatable | config | type index |
| `till_iType_Tillage` | integer(:) | per-event allocatable | config | type identifier |
| `till_iTT1` | integer(:) | per-type allocatable | config | first row in type table |
| `till_iTT2` | integer(:) | per-type allocatable | config | last row in type table |
| `till_TAB_Rho_tillage` | real(8)(:) | per-type-entry allocatable | config | post-tillage bulk density table |
| `till_TAB_Rho_cons` | real(8)(:) | per-type-entry allocatable | config | consolidated bulk density table |
| `till_TAB_K_R_cons` | real(8)(:) | per-type-entry allocatable | config | consolidation rate table |
| `till_TAB_Rho_match` | real(8)(:) | per-type-entry allocatable | config | matching density table |
| `till_TAB_N_match` | real(8)(:) | per-type-entry allocatable | config | matching n-value table |

These 12 are read-only-after-init (populated by `config_to_variables.f90` via the typed
`soil_tillage_t` config record). Never mutated at runtime. Strong candidate for exclusion
from `tillage_state_t` — could be moved to a typed `tillage_config_t` parameter bag, or
simply left as legacy globals if they are read-only after init.

### Group C — Per-layer working state (allocated at DoTillage(1), written each step)

| Variable | Type | Shape | Cadence | Notes |
|---|---|---|---|---|
| `till_Rho_tillage` | real(8)(:) | per-layer (NumLay) | step | post-event target density per layer |
| `till_Rho_cons` | real(8)(:) | per-layer | step | current consolidated density target |
| `till_Rho_last` | real(8)(:) | per-layer | step | density at start of current step |
| `till_K_R_cons` | real(8)(:) | per-layer | step | consolidation rate per layer |
| `till_Rho_match` | real(8)(:) | per-layer | step | matching-point density |
| `till_N_match` | real(8)(:) | per-layer | step | matching-point n-value |
| `till_Slope_match` | real(8)(:) | per-layer | step | slope at matching point |

### Group D — Per-step working scalars (written by Adapt_WC_H each call)

| Variable | Type | Shape | Cadence | Notes |
|---|---|---|---|---|
| `till_sumDWC` | real(8) | scalar | step | sum of water content changes |
| `till_sumAvail1` | real(8) | scalar | step | available pore space (wetting) |
| `till_sumAvail2` | real(8) | scalar | step | available water (draining) |

### Group E — Computed geometry scalars (computed once by det_MNSH / DoTillage(1))

| Variable | Type | Shape | Cadence | Notes |
|---|---|---|---|---|
| `till_MaxNumSoilHo` | integer | scalar | init-once | horizon count up to Max_Z_tillage |
| `till_MaxNumSoilCP` | integer | scalar | init-once | node count up to Max_Z_tillage |
| `till_iTill` | integer | scalar | step | current event-table index (advances as events fire) |

**Summary:** Groups A+B (18 fields) are read-only-after-init constants loaded from config.
Groups C+D+E (13 fields) are genuine runtime state that change during simulation.

**Reset cadence:** No `flzerointr` or `flzerocumu` anywhere in tillage. All 31 fields are
either init-once or instantaneous per-step. **Flat layout, no cohorts.**

---

## 3. External readers (5-category)

The grep survey (`src/` minus home file, variables.f90, initialize.f90, config_to_variables.f90)
found **8 external read/write sites** across **3 files**:

| File | Lines | Fields accessed | Category |
|---|---|---|---|
| `src/core/swap.f90` | 227, 328, 437 | (flTillage gate, not till_* directly) | Call-site arg |
| `src/io/toml/config_to_variables.f90` | 492, 1366–1382 | `till_swtill` (direct write); 17 fields via renaming use (pre-init write) | Co-writer (config path) |
| `src/io/toml/read_soil_tillage_toml.f90` | 29, 36–66 | `till_tbl` is a local TOML pointer, not a `till_*` global; writes to typed `tillage%` config record | Config reader (no till_* globals) |

**Key finding:** `read_soil_tillage_toml.f90` writes to the typed `soil_tillage_t` config
record (not to `till_*` globals directly). `config_to_variables.f90` is the bridge that
copies from the typed config into `till_*` globals. This is the standard pre-init pattern.

**Output files:** `grep -rn "till_"` on `swapoutput.f90` and `swap_csv_output.f90` returned
empty — till_* fields have no output readers in the I/O layer. DoTillage(3) writes debug
output directly via Fortran `write` statements to hardcoded unit numbers (222, 224, 226),
not via the output subsystem.

**No external compute or working-buffer uses** of `till_*` fields outside tillage.f90 itself.
The subsystem is unusually self-contained.

**Note on Bdens/ParamVG:** DoTillage writes legacy globals `Bdens` (per-horizon bulk density)
and `ParamVG` (per-layer VG parameters) as a coupling side-effect. These are NOT `till_*`
prefixed and are NOT owned by tillage in the naming convention. They are shared legacy
globals read by `soilhydraulics.f90`, `solute.f90`, and `oxygenstress.f90`. Tillage writes
them but does not own them. This coupling hazard is the arc's main design tension — see
Section 7.

---

## 4. Co-writers

Files that write `till_*` fields from outside `tillage.f90`:

| File | What it writes | Pre-init? | State-in-scope? |
|---|---|---|---|
| `src/io/toml/config_to_variables.f90:492` | `till_swtill = config%soil%swtill` (direct) | Yes — runs before `swap_main` loop | No `state` param; pre-init pattern |
| `src/io/toml/config_to_variables.f90:1366–1382` | 17 event/type-table fields via renaming use-association | Yes — bulk copy from typed config | No `state` param; pre-init pattern |
| `src/core/initialize.f90` | Zero-init (grep returned empty — no till_* in initialize.f90) | n/a | n/a |

**Finding:** `initialize.f90` does NOT zero-init `till_*` fields. This means the 31
globals rely entirely on Fortran default initialization (variables declared without `= 0`)
or on `config_to_variables.f90` to set them before first use. The allocatable arrays are
set by `config_to_variables.f90` before `DoTillage(1)` runs. The per-layer Group C arrays
(Rho_tillage etc.) are explicitly re-allocated inside `DoTillage(1, state)` at lines 62–68.

**config_to_variables pre-init pattern:** same as the soilwater-core arc's `cofgen` table
and the atmosphere arc's `nraida` seed. The adapter runs before any state allocator, so
`state%tillage` is not yet usable at that point. The design will follow the same pre-init
pattern: `config_to_variables.f90` continues to write the legacy globals; `DoTillage(1)`
(which runs after `state` is initialized) copies them into `state%tillage`.

---

## 5. Init-order analysis

```
swap.f90:194  call soilwater_init(state%soilwater, numnod, numlay)
swap.f90:227  if (flTillage) call DoTillage(1, state)
```

`soilwater_init` runs at line 194, `DoTillage(1)` at line 227. Confirmed: DoTillage does
not depend on soilwater_init for allocation (it allocates its own per-layer arrays), but
`Change_MvGpars` (which may be called transitively if a test-mode flag were set in iTask=1,
but is not called in the normal iTask=1 path) writes into `state%soilwater%cofgen` — which
IS allocated by soilwater_init. Order is safe.

**Tillage init pattern:** there is no `tillage_init` subroutine today. `DoTillage(1)` plays
that role: it allocates per-layer arrays (Group C), calls `set_iTill`, calls `det_MNSH`
(which writes `MaxNumSoilHo` / `MaxNumSoilCP`), and validates compatibility flags.

**Proposed `tillage_init` for this arc:** add a `tillage_init(state%tillage, numlay)` call
after `soilwater_init` (around line 195-200), allocating the per-layer Group C arrays into
`state%tillage`. Then DoTillage(1) calls `tillage_init` internally or the call site in
swap.f90 calls it. This follows the `soilwater_init` pattern exactly.

---

## 6. Config / Phase 0 candidates

```
src/config/soil_config.f90   # soil_tillage_t, soil_tillage_event_t, soil_tillage_type_t
```

ADR 0021 ported tillage input to TOML. The typed config `soil_tillage_t` (in
`soil_config.f90`) already covers:
- `i_n_model`, `iRedist` (scalars with defaults)
- `events(:)` (array of `soil_tillage_event_t`: date, z_tillage, i_tillage, type)
- `types(:)` (array of `soil_tillage_type_t`: rho_tillage, rho_cons, k_r_cons, rho_match, n_match per entry)

`till_swtill` is set from `config%soil%swtill` at `config_to_variables.f90:492`.

**Phase 0 is zero.** Every input parameter already has a typed config home. No stub-error
or legacy-only parameter work needed before the state migration can begin.

---

## 7. Coupling hazards

### H-1: Bdens and ParamVG are shared legacy globals, not till_* prefixed

`Change_Bdens` and `Consolidate_Bdens` mutate `Bdens(1:maho)` in-place. `Change_MvGpars`
mutates `ParamVG(1:21, 1:maho)` in-place. Both are declared in `variables.f90` (lines 761,
918) without `till_` prefix. They are read by:
- `soilhydraulics.f90:976-1014` — builds `cofgen` from `paramvg` at SoilWater(1) init
- `solute.f90:75-77` — reads `bdens(layer(i))` for sorption calculations
- `oxygenstress.f90:198` — reads `bdens(lay)` for soil density

These shared globals cannot be moved into `tillage_state_t` without also migrating their
readers (soilhydraulics, solute, oxygenstress). The safe option for this arc is:
**keep Bdens and ParamVG as legacy globals**; `state%tillage` will own only the `till_*`
fields. The arc does NOT attempt to move Bdens/ParamVG. This matches the playbook Lesson 6:
"don't expand scope mid-arc."

### H-2: DoTillage(3) output writes to hardcoded Fortran unit numbers

Lines 180–185: `write(222,...)`, `write(224,...)`, `write(226,...)`. These are active
unconditionally (the `if (TEST)` guard at line 178 is followed by un-guarded writes at
lines 184-185). This is a latent bug independent of this arc but worth noting: the output
section will serialize debug data to unit 222 and 226 even when `TEST=.false.`. Not a
migration hazard but flagged for cleanup.

### H-3: Pre-init pattern for Group A+B config fields

`config_to_variables.f90` populates 18 of the 31 `till_*` globals before `swap_state_t`
is even allocated. The migration strategy must keep these writes pointing at legacy
`till_*` globals during the pre-init window, then copy into `state%tillage` at
`tillage_init` / `DoTillage(1)` time. Same pattern as soilwater-core's `cofgen` handling.

### H-4: `till_iTill` is a mutable event-table cursor

`till_iTill` (the current event index) is mutated at runtime (`iTill = iTill + 1` at
tillage.f90:149). It is NOT a config constant. It belongs in the instantaneous flat layer
of `tillage_state_t`. The dual-write for this field must be established before the legacy
global is retired.

### H-5: Adapt_WC_H test outputs (units 124, 333, 444)

Lines 335, 363–364, 367–370: additional `write` statements to units 124, 333, 444 —
similarly not guarded by the `TEST` flag at the function level. These will produce output
files in every simulation run. Pre-existing issue, not introduced by this arc.

---

## 8. Cohort vs flat decision

**Recommendation: flat layout.** No `flzerointr` or `flzerocumu` accumulators exist in
tillage. All 31 `till_*` fields are either:
- Config-constant (Groups A+B): initialized once from config, never reset.
- Per-step instantaneous (Groups C+D): updated each step that tillage fires, no
  accumulation semantics.
- Init-once geometry (Group E): set by `det_MNSH` / `set_iTill` at DoTillage(1), never reset.

This matches the boundary subsystem precedent (ADR 0035: flat layout, 12 instantaneous
scalars, no cohorts).

**Allocation note:** the per-layer arrays (Group C, 7 fields) need `tillage_init` to
allocate them against `NumLay` — exactly as `soilwater_init` allocates per-node arrays.

---

## 9. Scope estimate

| Metric | Tillage | Boundary (ADR 0035) |
|---|---|---|
| Owned globals | 31 | 12 |
| External readers | 3 files / 8 sites | 14 files / 80 sites |
| Co-writers (outside home) | 1 file (config_to_variables) | 5 files |
| Phase 0 candidates | 0 | 7–8 |
| Allocatable arrays | 19 | 0 |
| Suggested tasks | 6–8 | 16 |

**Suggested task decomposition (6 tasks):**

| Task | Description |
|---|---|
| T-1 | `tillage_state_t` type + `tillage_init` subroutine in new `src/state/tillage_state.f90` (flat; 13 instantaneous fields from Groups C+D+E; allocatable per-layer; no cohorts) |
| T-2 | Add `type(tillage_state_t) :: tillage` to `swap_state_t`; add `call tillage_init(state%tillage, numlay)` in swap.f90 after soilwater_init |
| T-3 | Dual-write phase: DoTillage(1) copy config-seeded values → `state%tillage`; dual-write Group C+D+E writes in DoTillage(2), Consolidate_Bdens, Change_Tillage_Info, set_iTill, det_MNSH |
| T-4 | Reader cutover: replace all `till_*` local-alias reads in tillage.f90 with `state%tillage%*` reads (Groups C+D+E only; Groups A+B remain as aliases into legacy globals) |
| T-5 | Drop dual-write for Groups C+D+E; verify parity; retire 13 `till_*` globals from variables.f90 |
| T-6 | ADR 0039 + playbook update + merge |

Groups A+B (18 config-constant fields) are NOT migrated to `state%tillage` in this arc.
They remain as legacy globals populated by `config_to_variables.f90`. Rationale: they are
config parameters, not runtime state; their correct home is the typed `soil_tillage_t`
config record (already exists in `soil_config.f90`). A future "config consolidation" arc
can retire them. This keeps the arc tight.

---

## 10. Open questions

- **OQ-1 (Groups A+B exclusion):** Should Groups A+B (18 config-constant `till_*` fields)
  be migrated to `state%tillage` in this arc, or deferred? Recommendation: **defer.**
  They are already in typed `soil_tillage_t` config; migrating them to `state%tillage`
  would duplicate them. A config-consolidation arc (or the "retire all legacy globals"
  final arc) is the right venue. Migrating 13 vs 31 fields keeps this arc at 6 tasks.

- **OQ-2 (Bdens/ParamVG ownership):** Tillage writes `Bdens` and `ParamVG` as coupling
  side-effects. They cannot be owned by `state%tillage` without also migrating soilhydraulics,
  solute, oxygenstress readers. Recommendation: **keep as legacy globals for this arc.**
  Document as a future "soil parameters" arc candidate.

- **OQ-3 (Hardcoded unit numbers in DoTillage(3) and Adapt_WC_H):** Lines 180–185, 335,
  363–370, 367 write to units 122–124, 222, 224, 226, 333, 444 without conditional guards.
  These produce output files in every run. Recommendation: **fix the guard logic as a
  non-arc housekeeping task** before or during T-3 dual-write, since the output section
  is touched anyway.

- **OQ-4 (set_iTill logic bug):** `set_iTill` at line 412 has a tautological condition:
  `if (t1900 >= Date_tillage(i-1) .and. t1900 < Date_tillage(i-1))` — the same array
  element on both sides means this branch is always false. The loop never sets `iTill` for
  any event after the first. This appears to be a pre-existing latent bug. Recommendation:
  **note in arc design doc, fix in T-3 as a co-located cleanup** (condition should likely
  be `Date_tillage(i-1)` vs `Date_tillage(i)`).

- **OQ-5 (DoTillage(4) not called):** The `case (4)` (closure) has a `continue` body and
  is never called from `swap.f90`. Recommendation: **leave as-is** for this arc; it is a
  placeholder for future cleanup logic.

- **OQ-6 (swap_state_t component name):** Use `state%tillage` (matches existing naming
  convention: `state%soilwater`, `state%atmosphere`, `state%heat`, etc.).

- **OQ-7 (tillage_init signature):** Recommend `tillage_init(state, numlay)` where `state`
  is `type(tillage_state_t), intent(out)`. Per-layer arrays (7 fields) allocated to `numlay`.
  Geometry scalars (MaxNumSoilHo, MaxNumSoilCP) are computed by `det_MNSH` and set later
  during `DoTillage(1)`, not at `tillage_init` time.
