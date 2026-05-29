# ADR 0047 — Post-cfg polish arc

**Status:** Accepted (2026-05-29)
**Closes:** W8 + W9 + W10 + W12 + W5 (interpretation A) from the 2026-05-28 model-setup-layer analysis. Plus obsolete-tooling and dormant-file cleanup.

## Context

After the major refactoring chain (typed CSV → orchestrator dissolved → `state%cfg` retired — ADRs 0044/0045/0046), the analysis's remaining items were a mix of small architectural cleanups, dedup work, and end-of-strangler-pattern housekeeping. This arc bundles them into one branch — none of them individually justified an arc, but coherent together as "post-cfg polish."

W6 (defensive guards in soilwater) and W7 (CAPI path resolution) deferred per user direction.

## Decisions

### W9 — TOML reader dedup
Extracted `src/io/toml/toml_array_helpers.f90` with canonical `read_table_2d` / `read_real_array_1d` / `read_int_array_1d` / `read_real_pair_array`. Removed 5 verbatim copies of `read_table_2d` and 3+ variants of `read_array_1d` across the TOML readers. Two specialised variants (`read_real_array` scalar-fallback in solute, `read_array_1d_int` with count param in cropgrass) intentionally kept local — they have genuinely different interfaces. Net **−335 lines**.

### W8 — pathwork unification
`config%general%pathwork` becomes the canonical anchor for relative paths. Defaults to `directory_of(main_toml_path)` (was hardcoded `'./'`). `pathatm/pathcrop/pathdrain` default to `pathwork` (was `'./'` each). TOML subfile loaders (`read_crop_toml`, `read_drainage_toml`) drop their `base_path` parameter — they read `config%general%pathwork` directly (`read_general_toml` runs first, so pathwork is populated by then). Two competing conventions replaced by one. Test fixtures updated: removed explicit `[general.paths]` blocks that hardcoded `work = "./"` (resolved against project root, not case dir, post-unification).

### W5 (interpretation A) — extract crop subfile loader
`.crp.toml` file loading + dispatch moved from `read_crop_toml.f90` into a dedicated `src/io/toml/load_crop_rotation_files.f90`. `read_crop_toml.f90` becomes pure TOML parsing (no `toml_load`, no `path_helpers` dep, no `base_path` param — shrinks 128 → 83 lines). The new loader iterates `config%crop%rotation_file(:)`, resolves paths via `config%general%pathwork` (post W8 convention), and dispatches to `read_cropfixed/grass/wofost_toml`. Called from `load_swap_config` after `apply_section_readers`. Cleanly separates I/O from pure parsing.

### W10 — nutrients WSN dual-write retirement
The 8 dual-written `Wofost_Soil_Declarations` globals (`Amend`/`MatNum`/`VolaFrac`/`TimeAmend`/`NuAmend`/`iamend`/`namend`/`isme`) retired in favour of `state%nutrients` fields. `Wofost_SoilAmendents()` signature expanded to take `state`; `management_soil.f90`'s amendment-trigger block uses `state%nutrients%timeamend/isme`. This was **the only remaining "state init writes to globals" pattern** in the codebase. Irrigation defense-in-depth duplicate validation (separate item) explicitly NOT touched per user direction.

### W12 — `legacy_state` deletion
Recount after the cfg arc revealed every field on `legacy_state_t` had zero live readers. The entire `src/state/legacy_state.f90` file (379 lines), the `state%legacy` composition field on `swap_state_t`, and all build registrations DELETED. The Sweep 3 transitional bag (introduced 2026-05-22 to ease the variables.f90 retirement) is closed.

### Cleanup — obsolete scripts + dormant jongvanlier
4 obsolete python scripts deleted:
* `audit_variables.py` — audited the now-deleted `variables.f90`
* `output_parity_check.py` — checked the now-deleted `.crp`/`.snw` writers
* `retire_global.py` — helper for the now-complete GR-DVS arc
* `sweep1_config_side.py` — walker for the now-deleted `config_to_variables.f90`

`show_consumers.py` kept (generic symbol finder; docstring updated).

`src/crop/dormant/jongvanlier.f90` deleted — uncompiled museum file with a stale `state%cfg` reference. Stub in `rootextraction.f90` updated to remove the file path reference (points to git history instead).

## Consequences

+ **The strangler pattern is fully retired end-to-end.** `variables.f90` → `legacy_state_t` → both gone. `swap_state_t` contains only typed subsystem state records. No transitional bag remains.
+ **The last "state init writes to globals" pattern is gone** (W10). Nutrients no longer dual-writes to WSN globals; `Wofost_SoilAmendents` takes `state`.
+ **Single path convention** (pathwork as canonical anchor with TOML-dir default). User can still override at each level via `general.paths.{work,atmosphere,crop,drain}`.
+ **I/O cleanly separated from pure parsing** in the crop subsystem. `read_crop_toml` is a pure TOML parser; `load_crop_rotation_files` owns file I/O.
+ **Build tree is leaner.** Net **~−700 lines** across the arc (−335 from W9 dedup, −379 from legacy_state, plus smaller deltas from script deletion).

- Test fixtures had to drop explicit `work = "./"` overrides — minor change to 8 `swap.toml` files in `tests/swap-cases` submodule (separate commit there).

## Verification

- check-fast PASS: 4/4 byte-identical.
- check-full PASS: 5/5 non-xfail byte-identical (`winter` + `soilhysteresis` xfails pre-existing).
- 806 pFUnit tests pass (unchanged count — most changes were refactor-only, no new test surface).
- 6 commits on branch `polish-arc` (plus 1 separate commit on the `tests/swap-cases` submodule for the test-fixture `[general.paths]` block removal).

## What remains from the original analysis

- **W6**: defensive `if (.not. allocated())` guards in soilwater — user will handle.
- **W7**: CAPI hardcoded `base_dir = './'` — deferred (BMI/CAPI is its own arc).
- **W11**: `crop_irrigation_state_init` size — re-assessed during csv-families, no longer worth splitting.
- Irrigation defense-in-depth duplicate validation — deferred per user direction.

All major architectural items from the original analysis are now resolved. Anything remaining is either explicitly deferred or already addressed indirectly.
