# Crop Subsystem `use variables` Sweep Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Each Task below is one sequential subagent dispatch. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate every `use variables` import from active files in `src/crop/`. Replace per-field associate aliases with sub-record aliases (`soil`, `time`, `atmo`, …) matching the soil-sweep convention. Resolve cross-file symbols inline; no deferrals.

**Architecture:** Bottom-up file ordering (leaves → runtimes → dispatcher). Each file is one self-contained subagent sub-arc with multiple internal commits (associate refactor + per-symbol-cluster retirements + final cleanup). The orchestrator maintains a cross-file symbol ledger and dispatches subagents sequentially.

**Tech Stack:** Fortran 2003+, gfortran, Meson/Ninja build, pixi task runner, pFUnit test framework. State lives on `state%X%Y` typed records; config lives on `state%cfg%X` from TOML-loaded `swap_config_t`.

**Spec:** `docs/superpowers/specs/2026-05-25-crop-use-variables-sweep-design.md`

---

## Conventions (apply to every task)

**Associate pattern** (target form after refactor):

```fortran
associate( &
  crop => state%crop,            &  ! crop sub-state (common/wofost/grass/fixed/oxygen)
  soil => state%soilwater,       &  ! soil-water runtime
  mesh => state%mesh,            &  ! soil-mesh discretization
  atmo => state%atmosphere,      &  ! atmospheric runtime
  heat => state%heat,            &  ! heat runtime
  drai => state%drainage,        &  ! drainage runtime
  surf => state%surfacewater,    &  ! surface-water runtime
  time => state%timecontrol,     &  ! time control
  cfg_crop => state%cfg%crop,    &  ! crop config (when consumed)
  cfg_soil => state%cfg%soil     &  ! soil config (when consumed)
)
```

Inside the body, reads become `soil%qrot`, `time%t1900`, `cfg_soil%swhyst`, etc. — **no per-field prefix aliases** like `cw_qrot` / `tc_t1900`.

**Symbol classification** (every imported symbol must land in one bucket):

- **A — state-rebind**: state field exists; switch read to `state%X%Y`.
- **B — config-direct**: lives in `state%cfg%X`; read via `cfg_X` alias.
- **C — state-migrate**: no state field today; add to most natural sub-record in this commit.
- **D — config-migrate**: no config field today; add field + TOML reader wiring + validator in this commit.
- **E — retired-zero / orphan**: zero live readers OR always-zero; drop with `! [GR-CROP 2026-05-25] retired — <reason>` breadcrumb.
- **F — dormant**: legacy capability with no live dispatch; move body to `src/crop/dormant/<name>.f90` + fatalerr stub at unreachable call sites.

**Per-commit verification**:

```bash
pixi run check-fast
```

Must report `748 pFUnit + 4 regression cases passed`. State-schema or config-schema changes require:

```bash
rm -rf builddir && pixi run check-fast
```

**Allow-list for subagent edits** (each task's prompt names the exact files):

- The crop file(s) for this sub-arc.
- `src/core/variables.f90`, `src/state/legacy_state.f90`, `src/core/initialize.f90`, `src/io/toml/config_to_variables.f90` (for legacy global retirement).
- `src/state/crop_*.f90` (for Class C migrations).
- `src/config/*.f90`, `src/io/toml/read_*_toml.f90` (for Class D migrations).
- Test files in `tests/unit/` only when an import retirement breaks a test (e.g., `use variables, only: X` in a test).

**Cross-file ledger** — the orchestrator maintains a table:

| Symbol | Last consumer(s) | Status |
|---|---|---|
| `noddrz` | rootextraction, irrigation, helpers | reads-only; declaration retires in rootextraction sub-arc (largest) |

Each subagent's output must include any new ledger entries.

**Commit-message convention**:

```
refactor(crop): <file> — <one-line summary>

<body explaining what was retired, classification breakdown,
new state/config fields if any, cross-file leftovers if any>

check-fast: 748 pFUnit + 4/4 regression — byte-identical.
```

---

## Task 0: Cross-file symbol pre-flight discovery

**Goal:** Build the full cross-file symbol inventory and ledger before dispatching any retirement subagent. No code changes; output is a markdown table.

**Files:**
- Read: every file under `src/crop/*.f90`.
- Write: `docs/superpowers/plans/2026-05-25-crop-symbol-ledger.md` (new) — the cross-file symbol ledger.

- [ ] **Step 1: Subagent dispatch — discovery only**

Dispatch a single read-only subagent (Explore agent type) with this prompt:

> Read every `.f90` file in `src/crop/` that imports `use variables`. For each file, list every imported symbol with the file/line of the `use variables, only:` block it appears in. Then produce a consolidated cross-file ledger:
>
> | Symbol | Imported by (files) | Likely classification (A-F) | Best-guess target home | Already on state? |
> |---|---|---|---|---|
>
> For each symbol, run `grep -rn "\\b<symbol>\\b" src/` to find all read sites and writers. Flag any symbol that has writers in multiple files (potential dual-write trap).
>
> Also build a per-file summary:
>
> | File | use-block count | Total imports | Class A est. | Class B est. | Class C est. | Class D est. | Class E est. | Class F est. |
>
> Save the result to `docs/superpowers/plans/2026-05-25-crop-symbol-ledger.md`. Do not modify any source files.
>
> Suggested target homes to consult: `src/state/crop_common_state.f90`, `crop_wofost_state.f90`, `crop_grass_state.f90`, `crop_fixed_state.f90`, `crop_state.f90`, and `src/config/soil_config.f90` / `simulation_config.f90` for cross-domain symbols. The soil sweep recently retired ~21 globals — check `git log --oneline -30 -- src/core/variables.f90` to avoid re-classifying retired symbols.
>
> Output a final summary: total unique symbols across all 12 files; symbols imported by >1 file (cross-file); count of likely Class A/B (cheap) vs C/D (expensive) per file.

- [ ] **Step 2: Verify ledger written and committed**

```bash
ls -la docs/superpowers/plans/2026-05-25-crop-symbol-ledger.md
```

Expected: file exists with > 100 lines.

```bash
git add docs/superpowers/plans/2026-05-25-crop-symbol-ledger.md
git commit -m "docs(plans): crop subsystem symbol ledger (pre-flight discovery)"
```

- [ ] **Step 3: Read ledger and decide ordering refinements**

The orchestrator (parent) reads the ledger, confirms the planned ordering is still appropriate, and adjusts the per-file scope notes in subsequent tasks if surprises were found (e.g., a "leaf" file turns out to share many symbols with a runtime file).

---

## Task 1: `tillage.f90` sub-arc

**Goal:** Make `src/crop/tillage.f90` `use variables`-free. Resolve the `ParamVG` mutator (deferred from soil sweep).

**Files (allow-list for subagent edits):**
- `src/crop/tillage.f90`
- `src/state/soilwater_state.f90` (if (a) chosen: add `vg_params_layer(:)`)
- `src/state/legacy_state.f90`, `src/core/variables.f90`, `src/core/initialize.f90`, `src/io/toml/config_to_variables.f90` (for legacy global retirement)
- `tests/unit/state/test_soilwater_state.pf` (only if a new state field changes a test signature)

- [ ] **Step 1: Dispatch subagent for tillage sub-arc**

Subagent prompt (dispatch via `Agent` tool with `subagent_type: general-purpose`):

> Execute the `src/crop/tillage.f90` sub-arc of the crop `use variables` sweep. See spec at `docs/superpowers/specs/2026-05-25-crop-use-variables-sweep-design.md` and ledger at `docs/superpowers/plans/2026-05-25-crop-symbol-ledger.md` for context.
>
> Allow-list for edits: `src/crop/tillage.f90`, `src/state/soilwater_state.f90`, `src/state/legacy_state.f90`, `src/core/variables.f90`, `src/core/initialize.f90`, `src/io/toml/config_to_variables.f90`, and `tests/unit/state/test_soilwater_state.pf` only if a test breaks.
>
> Symbol inventory: read tillage.f90's `use variables, only:` block (lines ~14-33). It imports `swsolu`, `ParamVG`, `SwDiscrvert`, plus the `till_X` aliases for tillage SAVE state (already migrated to state%tillage but kept as legacy aliases for body compat).
>
> **Sub-arc commits**:
>
> 1. **Decision commit (no code)**: post a short message describing the ParamVG mutator situation and ask the user to pick (a) add `state%soilwater%vg_params_layer(:)` as a mutable layer-keyed VG store, or (b) keep `paramvg` alive permanently as tillage scratch space. Default (a). Wait for the user's response before proceeding.
>
> 2. **Associate refactor**: rewrite tillage's existing structure to use sub-record aliases (`mesh`, `soil`, `crop`, `time`, `cfg_soil`). Drop any `cw_X`/`tc_X`/`at_X` per-field alias prefixes. Byte-identical.
>
> 3. **`swsolu` cluster**: read `swsolu` via `cfg_soil%swsolu` if it has a config home; otherwise retire as orphan if `state%cfg%soil%swsolu` doesn't exist (check `src/config/soil_config.f90`). Drop legacy global if last consumer.
>
> 4. **`SwDiscrvert` cluster**: same pattern. The soil sweep marked it dormant via regrid.f90 — likely already retired or config-direct.
>
> 5. **`ParamVG` cluster** (the big one): if (a) is chosen — add `vg_params_layer(:)` to `soilwater_state_t`, populate from `state%cfg%soil%hydraulics` in `soilwater_init`, mutate it inside tillage instead of the bare global, and rebuild per-node `vg_params(:)` from it after each event. Retire `paramvg(21,maho)` from `variables.f90`, `legacy_state.f90`, `initialize.f90`, and `config_to_variables.f90`. If (b) — leave paramvg alive with a clear breadcrumb explaining it's a tillage-only mutable store.
>
> 6. **Final cleanup**: drop tillage.f90's `use variables` line entirely (or keep only a deprecation comment). Strip stale retirement-round tags (`[SS-GR-FINAL B6]`, `[SS-TIL T-5]`, etc.) from comments in the file body. Verify `grep -n "use variables" src/crop/tillage.f90` returns zero non-comment hits.
>
> **Every commit**:
>
> ```bash
> pixi run check-fast
> ```
>
> Must report `748 pFUnit + 4/4 regression — byte-identical`. State-schema commits (the ParamVG migration) need `rm -rf builddir && pixi run check-fast`.
>
> **Output contract**: a final message listing:
> - Each commit (hash + one-line summary).
> - Per-symbol disposition (Class A-F, target home, retirement status).
> - Any cross-file leftover (e.g., "swsolu still consumed by X — retire after X's sub-arc").
>
> **Do not** edit files outside the allow-list. **Do not** touch other crop files. **Do not** retire a legacy global until you verify (`grep -rn`) it has no remaining bare-name reader.

- [ ] **Step 2: Review subagent output**

Parent reviews:
- All commits passed `check-fast` byte-identical.
- ParamVG decision matches the user's pick.
- Cross-file ledger updated (likely no new entries since tillage is mostly self-contained).

If issues: dispatch a fix subagent. Else proceed.

- [ ] **Step 3: Update ledger and confirm next dispatch**

Update `docs/superpowers/plans/2026-05-25-crop-symbol-ledger.md` with the tillage sub-arc's leftovers (if any). Commit the ledger update.

---

## Task 2: `irrigation.f90` sub-arc

**Goal:** Make `src/crop/irrigation.f90` `use variables`-free.

**Files (allow-list):**
- `src/crop/irrigation.f90`
- `src/state/crop_common_state.f90` or new `src/state/irrigation_state.f90` (for Class C)
- `src/config/soil_config.f90` (likely; new `irrigation_t` sub-record) + `src/io/toml/read_soil_toml.f90`
- `src/state/legacy_state.f90`, `src/core/variables.f90`, `src/core/initialize.f90`, `src/io/toml/config_to_variables.f90`

- [ ] **Step 1: Dispatch subagent**

Subagent prompt:

> Execute the `src/crop/irrigation.f90` sub-arc. See spec and ledger.
>
> Allow-list: `src/crop/irrigation.f90`, `src/state/crop_common_state.f90`, `src/config/soil_config.f90`, `src/io/toml/read_soil_toml.f90`, `src/io/toml/config_to_variables.f90`, `src/state/legacy_state.f90`, `src/core/variables.f90`, `src/core/initialize.f90`. If you need a new state/config sub-record for irrigation runtime/config, propose it in your first message and wait for the parent's confirmation before creating it.
>
> Symbol inventory: irrigation.f90 has 2 `use variables` blocks (one at module level ~line 28 for `wclos`/`wcmes`/`wchis`/`isua`/`swcirrthres`/`noddrz`/`gird`; one in `SSDI_irrigation` at line ~330 for SSDI scheduling state). The SSDI block's persistent state needs `state%cfg%crop%irrigation` (new) or similar — propose the structure.
>
> **Sub-arc commits**:
>
> 1. **Associate refactor**: introduce sub-record aliases in both subroutines. Byte-identical.
>
> 2. **`noddrz` cluster** (reads-only this file): cut `noddrz` reads at lines 171-179 over to `state%crop%common%noddrz`. Drop `noddrz` from the use-list. **Do NOT** retire the `noddrz` declaration in `variables.f90` — rootextraction (sub-arc 3) is the last consumer. Update the cross-file ledger entry.
>
> 3. **`wclos`/`wcmes`/`wchis` cluster** (per-layer water content thresholds): these are crop-irrigation config staging. Class B if `state%cfg%X%wclos(:)` exists; Class D if not (add to a new `state%cfg%crop%irrigation` sub-record). Propose the structure in your first message.
>
> 4. **`isua`/`swcirrthres`/`gird` cluster**: irrigation runtime state. Likely Class C (add fields on `state%crop%common%irrigation_runtime` or a new sub-record).
>
> 5. **SSDI persistent state cluster** (`irrigevent`, `dt_SSDI_event`, `swssdi_irr`, `nod_ssdi_irr`, `ssdi_schedule_irr`, `ssdi_sched_type_irr`, `nod_ssdi_sensor_irr`, `ssdi_threshold_irr`, `ssdi_threshold_z_irr`, `ssdi_amount_irr`, `ssdi_appl_rate_irr`, `sw_interval_irr`, `days_interval_irr`, `days_counter_irr`, `nirri_ssdi_irr`, `ssdi_date_irr`, `ssdi_rate_f_irr`, `ssdi_amount_f_irr`): all currently legacy module variables. Already partially migrated (qssdi/qssdisum moved during soil sweep). The remaining 18 fields need a home — propose `state%crop%common%irrigation_ssdi` (or similar) as a new sub-record.
>
> 6. **Final cleanup**: drop irrigation.f90's `use variables` lines entirely. Strip stale `[SS-GR-FINAL B6]`, `[GR-BH C7]`, `[GR-SOIL 2026-05-24]`, `[SS-SWC S-2.12B]`, `[SS-TC TC-12]` tags from comments. Verify `grep -n "use variables" src/crop/irrigation.f90` returns zero non-comment hits.
>
> **Verification**: same gates as Task 1.
>
> **Output contract**: same format as Task 1.

- [ ] **Step 2: Review and update ledger**

- [ ] **Step 3: Confirm `noddrz` ledger entry is intact for Task 3**

---

## Task 3: `rootextraction.f90` sub-arc

**Goal:** Make `src/crop/rootextraction.f90` `use variables`-free. **Last-consumer retirement of `noddrz`**.

**Files (allow-list):**
- `src/crop/rootextraction.f90`
- `src/state/crop_common_state.f90`, `src/state/crop_state.f90` (Class C for any rootextraction-specific runtime state)
- `src/config/soil_config.f90` (Class D for `criterhr`, `wiltpoint`, `twilt`, etc. if they go to config)
- `src/io/toml/read_soil_toml.f90`
- `src/state/legacy_state.f90`, `src/core/variables.f90`, `src/core/initialize.f90`, `src/io/toml/config_to_variables.f90`
- After last-consumer retirement of `noddrz`: `src/crop/cropgrowth.f90` (drop the legacy writes at lines 171, 399-400 — keep only the state writes)

- [ ] **Step 1: Dispatch subagent**

Subagent prompt:

> Execute the `src/crop/rootextraction.f90` sub-arc. 3 subroutines (`RootExtraction`, `JongvanLier`, `JongvanLierLoop`) each have their own `use variables` block (~lines 44, 363, 736) plus a wofostnut-adjacent one at ~189.
>
> Symbols: `criterhr`, `flhydrlift`, `kroot`, `kstem`, `noddrz`, `oxygenintercept`, `oxygenslope`, `rootcoefa`, `rooteff`, `rootradius`, `rxylem`, `stephr`, `swfrost`, `twilt`, `wiltpoint`. Some are config (Class B/D), some are crop-runtime state (Class A/C), some may be retired-zero (Class E).
>
> **Sub-arc commits**:
>
> 1. **Associate refactor** (3 subroutines): replace per-field aliases (e.g., `cw_qrot` etc.) with sub-record aliases. Byte-identical.
>
> 2. **Config-side cluster** (Class B/D): `criterhr`, `twilt`, `wiltpoint`, `stephr`, etc. — most likely config snapshots. Check `state%cfg%crop` and `state%cfg%soil` for existing homes; create new fields on `state%cfg%crop%rootextraction` or similar if needed.
>
> 3. **Runtime cluster** (Class A/C): `kroot`/`kstem` (root depth densities — likely already on `state%crop%common`?), `oxygenintercept`/`oxygenslope`/`rootcoefa`/`rooteff` (already might be on `state%crop%common`).
>
> 4. **Retired-zero cluster** (Class E): `flhydrlift`, `swfrost` — check if they guard dead branches.
>
> 5. **`noddrz` last-consumer retirement**: after all 3 subroutines cut their ~20 read sites over to `state%crop%common%noddrz`, drop the legacy declaration from `variables.f90`, `legacy_state.f90`, `initialize.f90`. Also drop the dual writes in `src/crop/cropgrowth.f90` lines 171 and 399-400 (keep only the `state%crop%common%noddrz = ...` writes). Verify `grep -rn "\\bnoddrz\\b" src/ --include="*.f90"` returns no non-comment, non-state-field hits anywhere outside src/state/.
>
> 6. **Final cleanup**: drop rootextraction.f90's `use variables` blocks entirely. Strip stale tags. Verify zero non-comment hits in this file.
>
> **Verification**: same gates. The `noddrz` retirement is a cross-file change — verify check-fast passes on all 4 regression cases byte-identical.
>
> **Output contract**: same format as prior tasks; explicitly call out the noddrz retirement commit hash.

- [ ] **Step 2: Review and update ledger (noddrz row removed)**

---

## Task 4: `oxygenstress.f90` sub-arc

**Goal:** Make `src/crop/oxygenstress.f90` `use variables`-free. Create `state%crop%oxygen` sub-record.

**Files (allow-list):**
- `src/crop/oxygenstress.f90`
- `src/state/crop_common_state.f90` (or new `src/state/crop_oxygen_state.f90`) — add `crop_oxygen_state_t`
- `src/state/crop_state.f90` — add `oxygen :: crop_oxygen_state_t` field on `crop_state_t`
- `src/config/soil_config.f90` (Class D for oxygenstress config params)
- `src/io/toml/read_soil_toml.f90`
- `src/state/legacy_state.f90`, `src/core/variables.f90`, `src/core/initialize.f90`, `src/io/toml/config_to_variables.f90`

- [ ] **Step 1: Dispatch subagent**

Subagent prompt:

> Execute the `src/crop/oxygenstress.f90` sub-arc. 3 `use variables` blocks at the module level + inside subroutines. Major symbols: `c_mroot`, `f_senes`, `q10_root`, `q10_microbial`, `shape_factor_rootr`, `specific_resp_humus`, `max_resp_factor`, `c_top`, `rid`, `w_root_ss`, `tsoil` (staging), plus the persistent `module O2_pars` SAVE state at line 65: `o2_d_soil_term1`, `o2_d_soil_term2`, `o2_gfp100`, `o2_capac_term`, `o2_nmin1`, `o2_mplus1`, `o2_ini_stress`.
>
> **Sub-arc commits**:
>
> 1. **Add `state%crop%oxygen` sub-record**: create `src/state/crop_oxygen_state.f90` defining `crop_oxygen_state_t` with the O2 SAVE state fields (logical + 6 allocatable arrays sized to macp/numnod). Add an `init`/`aggregate` procedure following the same pattern as `crop_common_state_t`. Add the field `oxygen :: crop_oxygen_state_t` to `crop_state_t` in `src/state/crop_state.f90`. Wire allocation in the crop state aggregator. Byte-identical because no readers yet.
>
> 2. **Associate refactor**: replace any per-field aliases with sub-record style. Byte-identical.
>
> 3. **O2 SAVE-state migration**: cut over the `module O2_pars` reads/writes inside oxygenstress.f90 to `state%crop%oxygen%X`. Drop `module O2_pars` entirely if it has no other consumer. Verify check-fast byte-identical.
>
> 4. **Config snapshot cluster**: `c_mroot`, `f_senes`, `q10_root`, `q10_microbial`, `shape_factor_rootr`, `specific_resp_humus`, `max_resp_factor`, `c_top` — most are config-side. Add a `crop_oxygen_config_t` sub-record on `state%cfg%crop%oxygen` (new) and wire TOML reader. Cut over reads.
>
> 5. **Misc cluster**: `rid`, `w_root_ss`, `tsoil`. `rid` is in cropgrowth_helpers' use-list too — coordinate (likely cropgrowth_helpers needs it longer; leave declaration). `tsoil` staging buffer — drop the `=> tsoil` renaming in this file only.
>
> 6. **Final cleanup**: drop oxygenstress.f90's `use variables` blocks. Strip stale tags. Verify zero non-comment hits.
>
> **Verification**: state-schema changes require `rm -rf builddir && pixi run check-fast`. byte-identical.
>
> **Output contract**: same format. New sub-record fields list with their types and default values.

- [ ] **Step 2: Review and update ledger**

---

## Task 5: `cropgrowth_helpers.f90` sub-arc

**Goal:** Make `src/crop/cropgrowth_helpers.f90` `use variables`-free. 5 use-blocks (one per subroutine).

**Files (allow-list):**
- `src/crop/cropgrowth_helpers.f90`
- `src/state/crop_common_state.f90`, `src/state/crop_wofost_state.f90` (Class A/C as needed)
- `src/config/soil_config.f90` or new crop-specific config (`state%cfg%crop%germination`, `state%cfg%crop%co2`, `state%cfg%crop%wofost%rootgrowth`)
- `src/io/toml/read_soil_toml.f90` or new readers
- Standard legacy-cleanup files

- [ ] **Step 1: Dispatch subagent**

Subagent prompt:

> Execute the `src/crop/cropgrowth_helpers.f90` sub-arc. 5 subroutines each with their own `use variables`:
>
> - `cropoutput` (~line 45): `flCropOpenFile`, `outfil`, `cropfil`, `pathwork`, `project`, `crp`.
> - `ArableLandGerm` (~line 206): germination params (`flCropPrep`, `flCropSow`, `flCropGerm`, `dhPrep`, `hPrep`, `zPrep`, `dhSow`, `hSow`, `zSow`, `zTempSow`, `dtempSow`, `TempSow`, `MaxPrepDelay`, `MaxSowDelay`, `PrepDelay`, `SowDelay`, `tsumemeopt`, `tsumgerm`, `hdrygerm`, `hwetgerm`, `zgerm`, `TBASEM`, `TEFFMX`, `agerm`, `bgerm`, `cgerm`, `dummy_tsoil_alg_ => tsoil`).
> - `FacCO2` (~line 402): CO2 atmospheric tables (`co2year`, `mayrs`, `co2ppm`, `co2amaxtb`, `co2efftb`, `co2tratb`).
> - `update_rootdistribution` (~line 459): `gwrt`, `wrtmin`.
> - `sumttd` (~line 592): `tsumdepth`, `tsumtemp`, `tsumtime`, `pathwork`, `outfil`, `project`, `dummy_tsoil_sumttd_ => tsoil`.
>
> **Sub-arc commits**:
>
> 1. **Associate refactor** (all 5 subroutines): introduce sub-record aliases per subroutine. Byte-identical.
>
> 2. **Germination cluster** (`ArableLandGerm`): all 26 germination params likely Class D (config). Add a `state%cfg%crop%germination` sub-record on `soil_config_t` (or on `crop_config_global`'s rotation if it's per-rotation). Wire TOML reader. Cut over.
>
> 3. **CO2 cluster** (`FacCO2`): config tables, likely Class D. Add `state%cfg%crop%co2` or `state%cfg%atmosphere%co2` (CO2 is atmospheric — coordinate with atmosphere state).
>
> 4. **Root distribution cluster** (`update_rootdistribution`): `gwrt`, `wrtmin` likely Class D (config). Add to `state%cfg%crop%wofost%rootgrowth` or similar.
>
> 5. **Sumttd cluster**: `tsumdepth`/`tsumtemp`/`tsumtime` — stub-path Class D. The user previously deferred this; resolve now.
>
> 6. **File-path cluster** (`pathwork`, `outfil`, `project`, `cropfil`, `crp`, `flCropOpenFile`): file-path globals are out-of-scope for crop migration in principle. Check if they're already in a config home; if used by other subsystems too, leave the declarations alive but drop the imports from this file.
>
> 7. **Final cleanup**: drop all 5 `use variables` blocks. Strip stale tags. Verify zero non-comment hits.
>
> **Verification**: state-schema + config-schema changes require clean rebuild. byte-identical.
>
> **Output contract**: same format. New config sub-records listed.

- [ ] **Step 2: Review and update ledger**

---

## Task 6: `cropfixed_init.f90` + `cropfixed_runtime.f90` sub-arc

**Goal:** Make both `cropfixed_*.f90` files `use variables`-free.

**Files (allow-list):**
- `src/crop/cropfixed_init.f90`, `src/crop/cropfixed_runtime.f90`
- `src/state/crop_fixed_state.f90`, `src/state/crop_common_state.f90` (Class A/C)
- Standard legacy-cleanup files

- [ ] **Step 1: Dispatch subagent**

Subagent prompt:

> Execute the `src/crop/cropfixed_init.f90` + `src/crop/cropfixed_runtime.f90` sub-arc as a single dispatch (the two files share heavily). Each has one `use variables` block.
>
> Read both files' use-lists. Most symbols likely already on `state%crop%fixed` or `state%cfg%crop%fixed` (the crop_config_global rotation has fixed-crop subschemas). Class A/B rebinds dominant.
>
> **Sub-arc commits**:
>
> 1. **Associate refactor** (both files): sub-record aliases. Byte-identical.
> 2. **Per-cluster retirement** (size by cluster — typically 2-4 commits).
> 3. **Final cleanup** (both files): drop use-list lines, strip stale tags.
>
> **Verification**: same gates.
>
> **Output contract**: same format. Note any symbols still consumed by cropgrowth.f90 (most likely some flags).

- [ ] **Step 2: Review and update ledger**

---

## Task 7: `cropgrass_init.f90` + `cropgrass_runtime.f90` sub-arc

**Goal:** Make both `cropgrass_*.f90` files `use variables`-free.

**Files (allow-list):**
- `src/crop/cropgrass_init.f90`, `src/crop/cropgrass_runtime.f90`
- `src/state/crop_grass_state.f90`, `src/state/crop_common_state.f90`
- Standard legacy-cleanup files

- [ ] **Step 1: Dispatch subagent**

Subagent prompt:

> Execute the `src/crop/cropgrass_init.f90` + `src/crop/cropgrass_runtime.f90` sub-arc. cropgrass_init has 1 use-block; cropgrass_runtime has 1 use-block at module level + has its existing associate pattern at lines 132, 167 (per-rotation config access).
>
> Symbol clusters: ET-related (cropwofost shares these), grass-growth params (already on state%crop%grass via prior arc), mowing/grazing state.
>
> Same sub-arc structure as prior tasks. Notable: cropgrass_runtime.f90:167 already uses `crop_config_global%rotation_grass(state%crop%common%icrop)` — preserve this pattern; the use-variables block is for OTHER symbols.
>
> **Verification + output**: same as prior tasks.

- [ ] **Step 2: Review and update ledger**

---

## Task 8: `cropwofost_init.f90` + `cropwofost_runtime.f90` sub-arc

**Goal:** Make both `cropwofost_*.f90` files `use variables`-free. Largest sub-arc by lines (~2500).

**Files (allow-list):**
- `src/crop/cropwofost_init.f90`, `src/crop/cropwofost_runtime.f90`
- `src/state/crop_wofost_state.f90`, `src/state/crop_common_state.f90`
- Standard legacy-cleanup files

- [ ] **Step 1: Dispatch subagent**

Subagent prompt:

> Execute the `src/crop/cropwofost_init.f90` + `src/crop/cropwofost_runtime.f90` sub-arc. cropwofost_init has 2 use-blocks; cropwofost_runtime has 1 + existing associate at line 161.
>
> Largest pair (~2500 lines combined). Many symbols already on `state%crop%wofost`; expect mostly Class A/B with a few Class C for runtime state that hasn't been migrated yet.
>
> Cluster expectations: assimilation/biomass, leaf-area development, nutrient stress (delegates to wofostnut.f90 — already clean), CO2 (coordinate with helpers' FacCO2 cluster — if helpers already moved co2tables to `state%cfg%X`, just reuse), partitioning, harvest.
>
> Same sub-arc structure: 1 associate refactor commit + 4-6 cluster commits + 1 final cleanup. Possibly more — adjust as needed; report each in your output.
>
> **Verification + output**: same.

- [ ] **Step 2: Review and update ledger**

---

## Task 9: `cropgrowth.f90` dispatcher sub-arc

**Goal:** Make `src/crop/cropgrowth.f90` `use variables`-free. Dispatcher cleanup.

**Files (allow-list):**
- `src/crop/cropgrowth.f90`
- Standard legacy-cleanup files (for any final retirements)

- [ ] **Step 1: Dispatch subagent**

Subagent prompt:

> Execute the `src/crop/cropgrowth.f90` dispatcher sub-arc. 1 use-block at the top of `CropGrowth(task, tsoil, state)` around line ~67. Many of the symbols imported here are now retired or canonical thanks to the prior sub-arcs.
>
> Symbol classes likely after prior sub-arcs:
>
> - All `flCrop*` flags → `state%crop%common%flCrop*` (Class A rebind).
> - `icrop`, `daycrop`, `cropstart`, `cropend`, `swcrp`, `flCropNut`, `flHarvestDay` → `state%crop%common%X` (Class A).
> - `noddrz` → already retired by Task 3.
> - Nutrient flags `nlue`, `anlv`, `anst`, `nmxlv`, etc. → `state%crop%common%nutrient` or `state%crop%wofost%nutrient`.
> - Germination/sow/prep clusters → already migrated by Task 5 (helpers).
> - `dummy_tsoil_cg_ => tsoil` → drop the renaming. Verify `grep -rn "\\btsoil\\b" src/` shows no consumer outside state/heat — if so, drop the declaration in `variables.f90` too.
>
> **Sub-arc commits**:
>
> 1. **Associate refactor**: replace existing `tc_t1900`/`tc_daynr`/`at_tavd` per-field aliases with `time`/`atmo` sub-record. Byte-identical.
>
> 2. **Final cluster retirements**: cut over the remaining bare-name reads (most are flCrop* and crop-common state). Drop the legacy writes at lines 171 and 399-400 if Task 3 hasn't already.
>
> 3. **Final cleanup**: drop the `use variables` block entirely. Strip all stale `[SS-GR-CROPRT B*]`, `[GR-CROPWS B*]`, `[SS-GR-FINAL B*]`, etc. tags throughout the file. Verify `grep -n "use variables\\|Use Variables" src/crop/cropgrowth.f90` returns zero hits.
>
> **Verification + output**: same as prior tasks.

- [ ] **Step 2: Review and update ledger (should be empty after this task)**

---

## Task 10: Final audit

**Goal:** Verify the arc is complete and ledger is empty.

- [ ] **Step 1: Run repository-wide audit**

Dispatch a read-only subagent (Explore type) with this prompt:

> Audit the crop `use variables` sweep. Run these checks and report:
>
> ```bash
> grep -rn "use variables\\|Use Variables" src/crop/ --include="*.f90" | grep -v "^[^:]*:[^:]*:.*!.*$"
> ```
> Expected: zero lines (only comment-stripped output).
>
> ```bash
> grep -c "use variables" src/crop/*.f90
> ```
> Expected: every file shows 0 or only comment matches.
>
> ```bash
> git log --oneline 743a7ed..HEAD
> ```
> Expected: 35-50 commits since the spec commit.
>
> ```bash
> pixi run check-fast
> ```
> Expected: 748 pFUnit + 4/4 regression — byte-identical.
>
> Then read `docs/superpowers/plans/2026-05-25-crop-symbol-ledger.md` and confirm every entry is marked "retired" or moved to canonical state/config home with no leftovers.
>
> Report any deviations.

- [ ] **Step 2: Final commit — update memory + close arc**

```bash
# Update memory with arc completion
# Append to MEMORY.md a one-line entry referencing the new memory file
```

Create a project memory file at
`/home/zawadzkim/.claude/projects/-home-zawadzkim-Code-swap/memory/project_crop_use_variables_sweep_complete_2026-05-25.md`
summarizing:

- Total commits, files affected, symbols retired.
- New state sub-records added (`state%crop%oxygen`, irrigation, etc.).
- New config sub-records added.
- Any dormant modules created in `src/crop/dormant/`.
- Cross-file ledger final state (should be empty).
- Lessons learned for the next subsystem sweep.

Add a single-line index entry to `MEMORY.md`.

- [ ] **Step 3: Push (optional, on user request only)**

Branch is `development`. Do not push without explicit user instruction.

---

## Self-review notes

- Every task has its own subagent prompt that is self-contained (file allow-list, symbol scope, commit plan, verification gate, output contract). No "see Task N" cross-references.
- Cross-file symbol resolution mechanism (the ledger) is documented in the Conventions section and referenced from each Task.
- The `noddrz` last-consumer retirement is explicitly assigned to Task 3 (rootextraction) with the supporting Task 2 (irrigation) only doing read-site cutover.
- The `ParamVG` mutator design choice is an explicit decision gate inside Task 1 — subagent stops and asks before proceeding.
- State-schema commits flagged for `rm -rf builddir` clean-rebuild.
- Stale-tag stripping is in every file's final cleanup commit.
- Final audit (Task 10) verifies arc completion.
