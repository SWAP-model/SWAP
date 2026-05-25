# Crop Subsystem `use variables` Sweep — Design

**Date**: 2026-05-25
**Status**: design approved, awaiting implementation plan

## 1. Goal

Make every active file in `src/crop/` `use variables`-free. Replace the
existing per-field associate aliases (`cw_qrot`, `tc_t1900`, `at_ptra`,
…) with sub-record aliases (`soil`, `time`, `atmo`, …) matching the
convention settled in the soil-subsystem sweep (commits
`bd1a67f..2756abe`, 2026-05-24). Every retired legacy global ends with
either:

- a state field (created in the same commit if missing), or
- a config field (with TOML reader wiring added if missing), or
- explicit retirement documented in the commit (orphan, retired-zero,
  or dormant move).

No deferrals. Cross-file symbols are resolved when first encountered
(via a running ledger maintained by the orchestrator), not put on a
TODO list. Stale retirement-round tags from prior phases (`[SS-GR-*]`,
`[GR-CROPWS B*]`, `[SS-GR-CROPRT *]`, `[GR-SOL …]`, `! DEFERRED:`,
etc.) are stripped in the final cleanup commit for each file since
the symbols they reference will by then be either retired or
canonical.

## 2. Scope

Twelve files in `src/crop/` currently import `use variables` (one or
more `use variables, only: …` blocks each):

| Group | Files |
|---|---|
| Leaves | `tillage.f90`, `irrigation.f90` |
| Compute leaves | `rootextraction.f90`, `oxygenstress.f90` |
| Helpers | `cropgrowth_helpers.f90` (5 use-blocks) |
| Crop-type runtimes | `cropfixed_init.f90`, `cropfixed_runtime.f90`, `cropgrass_init.f90`, `cropgrass_runtime.f90`, `cropwofost_init.f90` (2 use-blocks), `cropwofost_runtime.f90` |
| Dispatcher | `cropgrowth.f90` |

**Already clean** (no `use variables`): `wofostnut.f90`, the
`wofost_soil_*.f90` set (8 files), `management_soil.f90`, and
`crop_config_global.f90`. These are not touched.

**Out of scope**:

- Refactoring `crop_config_global.f90` itself (it is consumed; not
  reshaped).
- Refactors to `src/state/crop_*.f90` beyond field additions.
- Anything in `wofost_soil_*.f90` (already clean).
- The legacy reader path (`readswap.f90`) and the legacy
  `variables.f90` module beyond declaration removal.

## 3. Approach

### 3.1 Sub-arc structure per file

Each file's sub-arc is a self-contained sequence of commits with
`pixi run check-fast` byte-identical between every commit.

1. **Pre-flight (no commit)** — read the file's `use variables`
   blocks; classify every imported symbol as one of:

   - **A — state-rebind**: the state field exists; switch reads
     to `state%X%Y` (e.g., `icrop` → `state%crop%common%icrop`).
   - **B — config-direct**: lives in `state%cfg%X`; read directly
     via a `cfg_X` associate alias.
   - **C — state-migrate**: needs a new state field; add to the
     most natural sub-record in this commit.
   - **D — config-migrate**: needs a new config field + TOML reader
     wiring + validator.
   - **E — retired-zero / orphan**: zero live readers or always-zero;
     drop with an inline `! [GR-CROP …] retired — orphan` breadcrumb.
   - **F — dormant**: a legacy capability with no live dispatch site;
     move the body to `src/crop/dormant/` with a reactivation
     checklist and stub-error fatalerr at the (now-unreachable) call
     site, mirroring the soil sweep's `sptabulated.f90` retirement.

2. **Associate refactor commit** — rewrite the file's associate
   block to sub-record style (`crop => state%crop`,
   `soil => state%soilwater`, `time => state%timecontrol`,
   `atmo => state%atmosphere`, `heat => state%heat`,
   `mesh => state%mesh`, `drai => state%drainage`,
   `surf => state%surfacewater`, plus `cfg_X => state%cfg%X` as
   needed). Bulk-rewrite per-field-prefix reads
   (`cw_qrot` → `soil%qrot`, `tc_t1900` → `time%t1900`,
   `at_ptra` → `atmo%ptra`, etc.). Bulk-rewrite bare-name reads of
   already-homed symbols (Class A) to the same sub-record form.
   **No state-schema or legacy-global changes in this commit.** This
   produces a byte-identical refactor that becomes the foundation
   for the cluster commits that follow.

3. **Per-symbol-cluster commits** — group the remaining symbols into
   logical clusters (e.g., germination, nutrients, biomass, CO2,
   irrigation-schedule). Each cluster is one commit that:

   - Adds new state/config fields if needed (Class C/D), including
     TOML reader updates for D.
   - Cuts over consumers within this file to the new home.
   - Drops the legacy globals from the use-list of this file.
   - If this file was the last consumer of a legacy global, also
     drops the declaration from `variables.f90`,
     `state/legacy_state.f90`, the zero-fill from
     `core/initialize.f90`, and any legacy mirror write in
     `io/toml/config_to_variables.f90`. Otherwise leaves them alive
     with a breadcrumb (`! still consumed by <other-file>, retire
     in that sub-arc`).
   - State-schema or config-schema commits trigger
     `rm -rf builddir` before `check-fast` per the lesson at
     `feedback_state_schema_clean_rebuild.md`.

4. **Final cleanup commit** — drop the file's `use variables, only:`
   line entirely (or keep only a one-line deprecation comment).
   Strip stale retirement-round tags from comments in the file body.
   Confirm `grep -n "use variables" <file>` returns zero hits.

### 3.2 Subagent prompt template

Each sub-arc is dispatched as one subagent run. The prompt is
self-contained and includes:

- **File path** and the canonical state-and-config homes the
  subagent is allowed to write to.
- **Symbol inventory** — every symbol from the file's `use variables`
  blocks, with proposed classification (A–F) and the target home or
  retirement breadcrumb.
- **Commit plan** — N+2 commits: associate refactor, N cluster
  commits, final cleanup.
- **Cross-file boundary contract**:
  - The subagent MAY drop a legacy global's declaration only when
    this file is its LAST consumer (verified by `grep -rn` against
    `src/`).
  - For not-last-consumer cases, the subagent retires only the read
    sites in this file and leaves the declaration alive with a
    breadcrumb naming the remaining consumers.
  - The subagent MUST NOT edit files outside the prompt's
    allow-list. Files allowed for edits: this file, the file's
    cross-file ledger entries (`variables.f90`,
    `state/legacy_state.f90`, `core/initialize.f90`,
    `io/toml/config_to_variables.f90`), state-schema files for
    Class C migrations, and config-schema/reader files for Class D
    migrations.
- **Verification gates** — every commit must pass
  `pixi run check-fast` (748 pFUnit + 4 regression cases) with
  byte-identical results; state/config schema changes require
  `rm -rf builddir` before `check-fast`.
- **Output contract** — the subagent's final message must include:
  - List of commits made (hashes + one-line summaries).
  - Symbol disposition table (one row per imported symbol; final
    classification + home + retirement status).
  - Cross-file leftovers: legacy globals NOT retired in this sub-arc
    because another crop file still reads them, named consumer per
    leftover.
  - Any deviations from the prompt's commit plan.

### 3.3 Cross-file ledger

The orchestrator (me, in the implementation phase) maintains a
running symbol ledger: each entry is one legacy global with its
current consumer set inside `src/crop/`. When a subagent reports a
"cross-file leftover", the ledger updates. When the last consumer's
sub-arc completes, the next sub-arc that touches that symbol's
declaration site retires the declaration.

This avoids the "all consumers must agree before any can land"
deadlock: each sub-arc retires what it can and bumps the ledger;
the cleanup naturally happens at the last hop.

## 4. File ordering & per-file sub-arc preview

Bottom-up. Each row is one sequential subagent dispatch.

| # | File(s) | Lines | Sub-arc commits (est.) | Notable items |
|---|---|---|---|---|
| 1 | `tillage.f90` | 577 | 2-3 | `swsolu`, `SwDiscrvert`, `ParamVG` mutator decision (see §5.2). Likely creates `state%soilwater%vg_params_layer(:)`. |
| 2 | `irrigation.f90` | 486 | 3-4 | `irrigevent`, `dt_SSDI_event`, SSDI scheduling state, `gird`, `wclos`/`wcmes`/`wchis`, `isua`, `swcirrthres`, `noddrz` read sites. New: `state%cfg%crop%irrigation` sub-record. |
| 3 | `rootextraction.f90` | 1012 | 4-5 | Three subroutine-scoped use-blocks. `criterhr`, `flhydrlift`, `kroot`, `kstem`, `oxygenintercept`, `oxygenslope`, `rootcoefa`, `rooteff`, `rootradius`, `rxylem`, `stephr`, `swfrost`, `twilt`, `wiltpoint`. **Last-consumer retirement of `noddrz`** (after files 2 and 5 cut over their reads). |
| 4 | `oxygenstress.f90` | 1812 | 4-5 | `c_mroot`, `f_senes`, `q10_root`, `q10_microbial`, `shape_factor_rootr`, `specific_resp_humus`, `max_resp_factor`, `c_top`, `rid`, `w_root_ss`, `tsoil` staging buffer. **New sub-record: `state%crop%oxygen`** (config snapshots + persistent O2 SAVE state currently in `module O2_pars`). |
| 5 | `cropgrowth_helpers.f90` | 762 | 5-6 | Five separate use-blocks, one per subroutine (`cropoutput`, `ArableLandGerm`, `FacCO2`, `update_rootdistribution`, `sumttd`). Germination params, CO2 tables (`co2year`, `mayrs`, `co2ppm`, `co2amaxtb`, `co2efftb`, `co2tratb`), `gwrt`/`wrtmin`, `tsumdepth`/`tsumtemp`/`tsumtime`. New homes: `state%cfg%crop%germination`, `state%cfg%crop%co2`, `state%cfg%crop%wofost%rootgrowth`. |
| 6 | `cropfixed_init.f90` + `cropfixed_runtime.f90` | 180 + 264 | 4-5 | Fixed-crop tabulated params. Mostly Class A/B rebinds via `state%cfg%crop%fixed` and `state%crop%fixed`. |
| 7 | `cropgrass_init.f90` + `cropgrass_runtime.f90` | 498 + 1356 | 5-6 | ET-related, grass-growth params. Symbols mostly already on `state%crop%grass`. |
| 8 | `cropwofost_init.f90` + `cropwofost_runtime.f90` | combined ~2500 | 6-8 | Largest pair. Wofost params + runtime (CO2, biomass, nutrients). Many Class A rebinds via `state%crop%wofost`. |
| 9 | `cropgrowth.f90` | 903 | 3-4 | Dispatcher. By this point most symbols are already retired. Associate refactor + cleanup pass + tag strip. |

Total estimate: **~35-45 commits across ~9 sub-arcs**.

## 5. Special items

### 5.1 `noddrz` last-consumer retirement

`noddrz` (compartment number at root-zone bottom) was partially
retired during the soil sweep (`6bdba67`). Remaining consumers in
`src/crop/`:

- `irrigation.f90` — loop bound + 4 conditionals at lines 171-179.
- `rootextraction.f90` — ~20 reads across 3 subroutines.
- `cropgrowth_helpers.f90` — line 486 (loop bound in
  `update_rootdistribution`).
- `cropgrowth.f90` — writers at lines 171 and 399-400, both already
  dual-write `state%crop%common%noddrz`.

After sub-arcs 2 (irrigation), 3 (rootextraction), 5 (helpers) cut
over reads, the rootextraction sub-arc (largest consumer; sub-arc 3)
will retire the legacy global declaration in its final cleanup
commit — since by then files 1, 2, and the soil sweep have already
stopped reading, and files 4, 6-9 don't read it. Cropgrowth's
dispatcher writers are dropped in sub-arc 9 (the dispatcher itself
ships only the state write).

### 5.2 `ParamVG` mutator design choice

`tillage.f90:242-264` mutates the per-layer Van Genuchten parameter
matrix `ParamVG(N, lay)` on each tillage event (rescaling theta_s,
K_sat, alpha, n_par from the new bulk density). The soil sweep
retired all OTHER `paramvg` consumers but left the legacy global
alive specifically for this mutator.

Two options for the tillage sub-arc (sub-arc 1):

**(a) Mutable copy on state** (proposed). Add
`state%soilwater%vg_params_layer(:)` (an array of
`vanGenuchten_params_t` records, one per soil-physical layer).
Populated by `soilwater_init` from `state%cfg%soil%hydraulics`.
Tillage mutates this; after each event, the existing per-node
`vg_params(:)` is rebuilt from it. This is the principled
long-term move and finally retires the legacy `paramvg(21, maho)`.

**(b) Keep paramvg alive permanently as tillage scratch space**.
Conservative fallback if (a) explodes scope.

The tillage subagent presents both at the start of the sub-arc;
the user picks. Default is (a).

### 5.3 `state%crop%oxygen` sub-record (new)

The `oxygenstress.f90` sub-arc adds `state%crop%oxygen` as a new
sub-record on `crop_state_t`. It hosts:

- The persistent O2 SAVE state currently in `module O2_pars`
  (`o2_d_soil_term1`, `o2_d_soil_term2`, `o2_gfp100`,
  `o2_capac_term`, `o2_nmin1`, `o2_mplus1`, `o2_ini_stress`).
- Config snapshots from `state%cfg%crop%X` for the
  oxygenstress-specific params (`c_mroot`, `f_senes`, `q10_root`,
  `q10_microbial`, `shape_factor_rootr`, `specific_resp_humus`,
  `max_resp_factor`, `c_top`).

Allocation gated on the swoxygen feature flag if applicable.

### 5.4 `tsoil` config-staging buffer

Three crop files use the renaming pattern
`dummy_tsoil_X_ => tsoil` (verified via grep):
`cropgrowth.f90`, `cropgrowth_helpers.f90` (twice — `ArableLandGerm`
and `sumttd`), and `cropgrass_runtime.f90`. The renaming imports a
stale config-staging buffer that has no live reader after the heat
subsystem's tsoil migration (per memory
`project_atmosphere_arc_complete_2026-05-23.md`). Each sub-arc drops
the renaming locally in its associate-refactor commit.
The `tsoil` declaration itself remains in `variables.f90` only if
another non-crop reader exists — verified at the dispatcher sub-arc
(#9) with a final `grep -rn "\\btsoil\\b" src/` against non-comment
hits; if zero remain outside `variables.f90`/`legacy_state.f90`/
`initialize.f90`/`config_to_variables.f90`, the declaration is also
dropped.

## 6. Risks

1. **Dual-write detection traps**: many crop globals are written to
   BOTH the legacy global and the state field today (e.g.,
   `cropgrowth.f90:171 state%crop%common%noddrz = noddrz`).
   Subagents drop the legacy write *and* verify no remaining
   bare-name reader before retirement. Mitigation: each subagent
   pre-flight greps `\\b<symbol>\\b` across all of `src/` to confirm.

2. **Cross-file ledger drift**: if a "last consumer" is misjudged,
   the symbol retires too early and a later file fails to compile.
   Mitigation: subagent verifies `grep -rn` for non-comment,
   non-declaration hits returns zero before dropping declarations.

3. **TOML schema gates**: some imported symbols guard dormant
   compute branches. If the branch is dead at runtime, prefer
   retired-zero + fatalerr stub (Class E) over migration. Subagent
   flags case-by-case.

4. **`ParamVG` mutator scope** (see §5.2): option (a) is a
   non-trivial state-schema change.

5. **Subagent context bleed**: subagents see their prompt only. The
   orchestrator must pass each subagent a snapshot of the
   cross-file ledger relevant to that file.

## 7. Verification gates

Every commit:

- `pixi run check-fast` — 748 pFUnit + 4 regression cases
  (`hupselbrook`, `surfacewater`, `grassgrowth`, `salinitystress`).
  Byte-identical.
- State-schema or config-schema commits: mandatory
  `rm -rf builddir` clean rebuild beforehand, per
  `feedback_state_schema_clean_rebuild.md`.

Every sub-arc's final commit:

- `grep -n "use variables\\|Use Variables" <file>` returns zero
  non-comment hits.

The whole arc's final-final commit:

- `grep -rn "use variables\\|Use Variables" src/crop/ --include="*.f90"`
  returns zero non-comment hits.
- variables.f90 line count audit: confirm the retired symbol set's
  declarations are gone.

## 8. Out of scope (re-stated)

- `crop_config_global.f90` refactor.
- `wofost_soil_*.f90` (already clean).
- `src/state/crop_*.f90` beyond field additions.
- `readswap.f90` and the legacy reader path.
- `variables.f90` beyond declaration removal of symbols this arc
  retires.
- Symbol retirements outside `src/crop/`'s ownership (e.g., heat,
  drainage symbols imported through `tsoil` or similar).

## 9. Success criteria

1. Every file in `src/crop/` matching `grep -l "use variables"` is
   `use variables`-free at arc end.
2. Every commit in the arc passes `pixi run check-fast` byte-identical.
3. Every retired legacy global has a clear home (state, config) or
   documented retirement (orphan, dormant) in its retirement commit.
4. The associate-pattern convention across `src/crop/` matches the
   soil sweep: sub-record aliases (`soil`, `time`, `atmo`, …) and
   no remaining `cw_X`/`tc_X`/`at_X` per-field prefix aliases.
5. Stale retirement-round tags from prior phases stripped in each
   file's final cleanup commit.
6. Cross-file ledger reaches empty (no symbol left dangling because
   "another file still reads it").
