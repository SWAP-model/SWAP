---
title: "Phase 4c-a — Cross-File TOML + Fixed and Grass Crops"
author: Mateusz Zawadzki
date: 2026-04-25
status: draft
---

# Phase 4c-a: Cross-File TOML + Fixed and Grass Crops

The first half of Phase 4c. Brings the new TOML pipeline to the point where four of the five non-macropore regression cases (excluding hupselbrook's WOFOST entry and the salinitystress case, both deferred to Phase 4c-b) reach full field-by-field parity against the legacy `variables` globals — including the `.dra` and `.crp` file contents that Phase 4b deferred.

## Background

Phase 4b matured `.swp` parity for hupselbrook and surfaced three legacy quirks that were workaround-papered:

- `tstart`/`tend` parity tolerance widened to 1.5 days because legacy ttutil uses 1-based JDN-since-1900 while the new TOML helper uses 0-based.
- `swmacro` parity assertion dropped because `readswap.f90:20` declares a local `integer SwMacro` that shadows the `use variables` import, so the module global is never written.
- `nrlevs` parity assertion dropped because `readswap.f90` hard-codes `nrlevs = 1` whenever `dramet /= 3`.

Phase 4c-a fixes all three at the source. The user accepted the policy that "small targeted main-program changes are preferred over carrying complex workarounds" — these three fit that brief.

The Phase 4b parity test asserted only `.swp` fields. Each case's drainage and crop content lives in side files (`.dra`, `.crp`) that the legacy reader follows from inside `readswap`. Phase 4b deferred cross-file TOML loading; Phase 4c-a delivers it. The new loader follows explicit `file = "..."` references the same way the legacy reader follows path conventions.

The crop work splits across two sub-phases because the WOFOST type 2 model is much larger (~150+ params with ~30 tables) than the type 1 (fixed, ~80 params) and type 3 (grass, ~110 params) types. Phase 4c-a delivers types 1 and 3 only; Phase 4c-b adds WOFOST.

## Scope

### In

- **Cross-file TOML loading**, opt-in per section: `[drainage].file = "..."` and `[[crop.rotation]].file = "..."`. Loader follows references and merges contents. Inline still works (Phase 4b shape is preserved as a fallback).
- **Two new config types** in `src/config/`: `cropfixed_config_t` (~80 fields) and `cropgrass_config_t` (~110 fields). Each declares its own complete field set (no shared base type — see Decision section).
- **Two new section readers** in `src/io/toml/`: `read_cropfixed_toml.f90` and `read_cropgrass_toml.f90`.
- **Path-resolution helper** `path_helpers.f90` in `src/io/toml/` exposing `resolve_relative_path(base, rel)` and `directory_of(path)`.
- **Extensions to** `crop_config_t` (parallel arrays `rotation_fixed(:)` and `rotation_grass(:)`), `read_drainage_toml.f90` (file-following), `read_crop_toml.f90` (per-entry file-following + dispatch by type), `load_swap_config.f90` (passes `base_path` to file-following readers).
- **Three small targeted source fixes**:
  - `parse_date_to_days1900` aligned to legacy 1-based JDN.
  - `integer SwMacro` local removed from `readswap.f90:20`.
  - `nrlevs` mirror added to `drainage_config_t%finalize`.
- **Per-case TOML files** under `tests/swap-cases/toml/<case>/`:
  - `1.hupselbrook/`: existing `swap.toml` extended with cross-file references; new `swap.dra.toml`; new `maizes.crp.toml` and `grassd.crp.toml`; placeholder `potatod.crp.toml` for 4c-b.
  - `2.grassgrowth/`, `4.oxygenstress/`, `6.surfacewater/`: existing `swap.toml` extended; new `swap.dra.toml` and matching `*.crp.toml`.
- **Per-case parity tests**: new `test_grassgrowth_parity.pf`, `test_oxygenstress_parity.pf`, `test_surfacewater_parity.pf`. Existing `test_hupselbrook_parity.pf` extended with type 1 + type 3 sub-section assertions. The Phase 4b workaround tolerances are tightened (date parity to `1.0d-9`); the dropped `swmacro` and `nrlevs` assertions are re-enabled.
- **Schema doc refresh**: `docs/configuration-schema.md` reconciled to Phase 4b/4c-a reality.

### Out

- **WOFOST (type 2)** — `cropwofost_config_t` and parity for hupselbrook's potatod entry + the entire 5.salinitystress case. Phase 4c-b.
- **Macropore case** (3.macroporeflow) — out for all crop parity work in 4c-x.
- **Wiring the new TOML path into `swap.f90` execution** — still legacy-readswap-only at runtime. Phase 4 proper.
- **Bottom_boundary, heat, solute config types** — Phase 4d or later.
- **Mid-rotation crop parity** — only the LAST rotation entry's loaded `variables` globals can be asserted (legacy reader sequentially overwrites globals during `readswap`). Documented as a 4c-a limitation; mitigation deferred.
- **Retiring `readswap.f90`** — Phase 4 proper.

## Architecture

### Single-entry recursive loader

`load_swap_config(path, config, errors)` reads `path`, dispatches to section readers, and the file-following readers (drainage, crop) recursively load and merge their referenced files.

```
load_swap_config(path)  e.g. "tests/swap-cases/toml/2.grassgrowth/swap.toml"
   ├─ read swap.toml root
   ├─ read_general_toml      → config%general
   ├─ read_simulation_toml   → config%simulation
   ├─ read_meteorology_toml  → config%meteo
   ├─ read_drainage_toml     → if [drainage].file present, load <swap.toml's dir>/<file>
   │                            then read [drainage] from THAT root → config%drain
   │                          (inline path still works as fallback)
   ├─ read_soil_toml         → config%soil
   ├─ read_crop_toml         → reads [crop] + [[crop.rotation]] entries
   │                            for each entry with file= and type=:
   │                              load <swap.toml's dir>/<file>
   │                              dispatch by type:
   │                                type 1 → read_cropfixed_toml → config%crop%rotation_fixed(i)
   │                                type 3 → read_cropgrass_toml → config%crop%rotation_grass(i)
   │                                type 2 → (4c-b) read_cropwofost_toml
   ↓
config (fully populated)
   ↓
config%validate(errors)
config%finalize(errors)        nrlevs mirror happens here
errors%abort_if_fatal()
```

### Path resolution

Every `file = "..."` is resolved RELATIVE to the directory containing the file that referenced it. `swap.toml`'s references resolve against `swap.toml`'s directory; sub-file references would resolve against their containing file (no nested references in 4c-a — `.dra.toml` and `.crp.toml` are leaves).

Implementation:
- `directory_of(path)` returns the directory part (everything up to and including the last `/`).
- `resolve_relative_path(base, rel)`:
  - If `rel` is absolute (starts with `/`), return it.
  - Otherwise concatenate `base` + `rel`.
  - Phase 4c-a only handles same-directory and below-directory paths. `../` may appear without crashing (just becomes part of the string), but case authors should avoid it.

## Module design

### `src/config/cropfixed_config.f90`

`cropfixed_config_t` — ~80 fields grouped by category:

- **Phenology** (`idev`, `lcc` if `idev=1`)
- **Light & growth** (`kdif`, `kdir`, `eff`)
- **Crop factor & height tables** (`cftb`, `chtb`)
- **Root growth** (`rdi`, `rri`, `rdc`, `rdctb`)
- **Water stress (Feddes)** (`hlim1`...`hlim4`, `adcrh`, `adcrl`, `rsc`)
- **Oxygen stress** (`swhydrlift`, `oxstress_type`, plus model-specific)
- **Salinity stress** (`ecmax`, `ecslop`)
- **Interception** (`cofab`)
- **Sowing/planting** (`crpsta`)

Type-bound `validate` checks ranges/enums. Type-bound `finalize` is a no-op for now; added if parity reveals a need.

### `src/config/cropgrass_config.f90`

`cropgrass_config_t` — ~110 fields. Extends the fixed-crop categories with grass-specific:
- **Mowing schedule** (`swharv`, `nmow`, `dates_mowing`, `lai_after_mow`)
- **Grazing** (`swgraz`, `nstart_graz`, `nstop_graz`, livestock parameters)
- **Fertilizer table** (per-cycle N and timing)

Same validate/finalize pattern.

Both types are independent declarations. Field-name overlap (e.g., `hlim1`) is accepted for clarity; deduplication (shared base type) is YAGNI for Phase 4c-a and reconsidered when Phase 4c-b lands the third type.

### `crop_config_t` extension

```fortran
type :: crop_config_t
   integer :: swcrop = 0
   real(real64),       allocatable :: rotation_start(:)        ! existing
   real(real64),       allocatable :: rotation_end(:)
   character(len=256), allocatable :: rotation_file(:)
   integer,            allocatable :: rotation_type(:)
   type(cropfixed_config_t), allocatable :: rotation_fixed(:)   ! NEW
   type(cropgrass_config_t), allocatable :: rotation_grass(:)   ! NEW
contains
   procedure :: validate => crop_config_validate                ! extended
   procedure :: finalize => crop_config_finalize
end type
```

`rotation_fixed` and `rotation_grass` allocated to length N (rotation length) when the rotation table is parsed. For each `i`, only the entry corresponding to `rotation_type(i)` is meaningful; others stay default-initialized.

`crop_config_validate` extension calls `self%rotation_fixed(i)%validate(errors)` for each `i` where `rotation_type(i) == 1`, and `self%rotation_grass(i)%validate(errors)` where `rotation_type(i) == 3`.

### Cross-file loader pattern

`read_drainage_toml(doc_root, config, errors, base_path)` — `base_path` is optional. If `[drainage].file` is present and `base_path` provided, the reader loads the referenced file and reads `[drainage]` from THAT root via the inner helper `read_drainage_inner(table_or_root, config, errors)`. Inline path: `read_drainage_inner` reads from the original `[drainage]` table.

Same shape for `read_crop_toml`. The rotation walker examines each entry's `file` key and dispatches by type to the matching crop reader. Crop readers (`read_cropfixed_toml` / `read_cropgrass_toml`) take a TOML root pointer (the loaded `.crp.toml` document) and a destination `cropfixed_config_t` or `cropgrass_config_t`.

### `load_swap_config` change

```fortran
subroutine load_swap_config(path, config, errors)
   character(len=:), allocatable :: base_dir
   ! ... load doc ...
   base_dir = directory_of(trim(path))
   call read_drainage_toml(doc_ptr, config%drain, errors, base_path=base_dir)
   call read_crop_toml    (doc_ptr, config%crop,  errors, base_path=base_dir)
   ! other section readers don't take base_path — inline-only
end subroutine
```

## Targeted source fixes

### Fix 1: `parse_date_to_days1900` align to legacy 1-based

`src/io/toml/toml_field_helpers.f90`:

```diff
- jd1900 = 2415021
+ jd1900 = 2415020
```

After fix: `2002-01-01` → 37257 (matches legacy `tstart`). Existing `test_parse_date_to_days1900` constant updates from `37255` to `37257`. Existing round-trip test (`test_hupselbrook_roundtrip`) recomputes; the `tstart`/`tend` tolerance in `test_hupselbrook_parity.pf` tightens from `1.5d0` back to `1.0d-9`.

### Fix 2: `swmacro` shadow removal

`src/io/readswap.f90:20`:

```diff
-      integer SwMacro
```

After fix: `use variables` resolves `swmacro` to the module global; `rdsinr('swmacro', ...)` writes the global; downstream code (and parity tests) see the user's actual value. Verified by `check-full` 6/6 green (the macropore case still uses `flmacropore` for runtime branching, unaffected).

### Fix 3: `nrlevs` mirror in `drainage_config_t%finalize`

`src/config/drainage_config.f90`:

```diff
   subroutine drainage_config_finalize(self, errors)
      class(drainage_config_t), intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
+     ! Mirror legacy convention from readswap.f90: single-level drainage
+     ! methods (dramet 1 or 2) clobber nrlevs to 1 regardless of input.
+     ! Matching this is required for parity. The validator already restricts
+     ! nrlevs to [0, 5]; this finalize clobber lands AFTER validate.
+     if (self%dramet /= 3) self%nrlevs = 1
   end subroutine
```

Validator behavior is unchanged — it still permits `nrlevs in [0, 5]`. The finalize clobber for `dramet /= 3` lands after validate, which matches the legacy reader's "read-then-overwrite" sequence. The dedicated unit test asserts both: (a) when `dramet=2` and any input nrlevs, post-finalize nrlevs is 1; (b) when `dramet=3`, finalize leaves nrlevs alone.

New unit test in `tests/unit/config/test_drainage_config.pf` for the rule. Hupselbrook's `swap.toml` updated to `nrlevs = 0` (truthful to its `dramet=2` config); finalize normalizes to 1; parity assertion succeeds.

## Per-case TOML authoring

Inside the `tests/swap-cases` submodule (branch `main`):

| Case | swap.toml | swap.dra.toml | *.crp.toml | Notes |
|---|---|---|---|---|
| 1.hupselbrook | extend with file refs | new | maizes.crp.toml + grassd.crp.toml + potatod.crp.toml stub | Type 2 entry's stub is a placeholder; Phase 4c-b fills it. |
| 2.grassgrowth | extend | new | grassd.crp.toml | Fully covered (type 3 only). |
| 4.oxygenstress | extend | new | grassd.crp.toml (per-case copy; values may differ) | Fully covered. |
| 6.surfacewater | extend | new | grass.crp.toml (type 1) | Fully covered. |
| 5.salinitystress | NOT touched | — | — | Whole case deferred to 4c-b. |

Hand-authored from the matching legacy files. Estimated effort: ~1-2 hours per `.crp.toml` × ~5 unique crop files; ~30 minutes per `.dra.toml` × 4 cases.

## Test plan

### Unit coverage (extends existing patterns)

- `cropfixed_config_t`: ~5-8 validator tests (range, enum, required-presence on a representative sample of fields; not all 80).
- `cropgrass_config_t`: ~6-10 validator tests covering shared + grass-specific fields.
- `read_cropfixed_toml`: happy + missing-section + malformed (3 tests).
- `read_cropgrass_toml`: happy + missing-section + malformed (3 tests).
- `path_helpers`: 4 tests covering `directory_of` and `resolve_relative_path` (absolute, relative, edge cases).
- `drainage_config_t%finalize` nrlevs mirror: 2 tests (when applies, when doesn't).

### Cross-file integration

- One smoke test per in-scope case loads + validates + finalizes from `swap.toml` end-to-end (cross-file references followed). Pattern from Phase 4a's smoke test, extended to confirm `rotation_fixed(i)` or `rotation_grass(i)` is populated for relevant entries.

### Parity tests (the headline)

One pFUnit test file per case, modeled on `test_hupselbrook_parity.pf`:
- `test_grassgrowth_parity.pf` (full type 3 parity, ~10 fields beyond `.swp`)
- `test_oxygenstress_parity.pf` (full type 3 parity)
- `test_surfacewater_parity.pf` (full type 1 parity)
- `test_hupselbrook_parity.pf` extended: type 1 (maizes) + type 3 (grassd) sub-section assertions. Type 2 (potatod) stays unasserted with comment.

Each parity test:
- chdirs into case dir, stages template, calls `readswap()`, chdirs back.
- Calls `load_swap_config` against the case's TOML.
- Asserts `variables%foo == config%section%field` for ~30-50 fields, of which ~10-20 are crop-specific (the new contribution).
- Cropdata parity asserts the LAST rotation entry's loaded data only (legacy globals state).

Test count delta: ~6 tests/case × 3 cases + 2 new tests on hupselbrook + ~22 new unit tests ≈ +42 tests. End count around 151.

## Schema doc refresh

`docs/configuration-schema.md` reconciled in one focused commit:
- `[crop.rotation]` documented as `[[array-of-tables]]` (matches Phase 4b reader).
- `[output]` keys nested under `[simulation.output]` (matches Phase 4b reader).
- `[meteorology].alt`/`altw` documented at the section root, not under `[meteorology.evapotranspiration]`.
- New "cross-file references" subsection covering `[drainage].file` and `[[crop.rotation]].file` plus `cropfixed`/`cropgrass` schema sketches.
- Cross-link to `docs/toml-format-guide.md`.

## Exit criteria

1. `pixi run -e test check-full` 6/6 green.
2. `pixi run -e test test-pfunit` clean — ~151 tests passing.
3. `parse_date_to_days1900` aligned to legacy 1-based; round-trip test passes; Phase 4b parity tolerance narrowed back to `1.0d-9`.
4. `swmacro` shadow removed; parity assertion on `swmacro` re-enabled and passing in hupselbrook test.
5. `drainage_config_t%finalize` mirrors `nrlevs=1` when `dramet /= 3`; dedicated unit test passes; hupselbrook parity asserts `nrlevs` and passes.
6. Cross-file loader follows `[drainage].file` and `[[crop.rotation]].file` references; load + validate + finalize succeed for the 4 in-scope cases.
7. Per-case parity tests pass for cases 2, 4, 6 (full); case 1 (partial — type 1 + type 3 only).
8. `docs/configuration-schema.md` refreshed; FORD builds clean.
9. Coverage tracked + recorded; no per-domain drop vs Phase 4b baseline.
10. `main` fast-forwarded to `development` locally; tag `rescue/phase-4c-a-crop-fixed-grass`. No push.

## Risks and mitigations

| Risk | Mitigation |
|---|---|
| Legacy `readcropfixed`/`readcropgrass` write to globals not in `variables` (e.g., common blocks, shared structures) — parity assertions can't reach them | Same as Phase 4b: per-field grep of the legacy reader. Drop assertion + log as Phase 4d cleanup if a global isn't reachable. |
| Cross-file path resolution surprises (`../`, absolute paths, symlinks) | `resolve_relative_path` handles only same-dir and below. Author all `.crp.toml` in `swap.toml`'s directory. |
| `swmacro` removal breaks something downstream that relied on the local defaulting to 0 | `check-full` is the gate. Macropore case (sets `swmacro=1`) is the crucial regression. |
| Per-entry crop parity asserts only on LAST rotation entry | Documented as 4c-a limitation. Mid-rotation parity needs a "load-and-snapshot" wrapper, deferred. |
| `.crp.toml` authoring is tedious | Audit doc serves as field checklist. Hand-author one .crp.toml first as template. |
| Hupselbrook's potatod (type 2) stub is incomplete and confuses readers | Stub is a single-line `# Phase 4c-b will populate`. Document in case README. |

## Follow-on (Phase 4c-b and beyond)

After Phase 4c-a tags:
- **Phase 4c-b**: `cropwofost_config_t` (~150+ fields, ~30 tables); `read_cropwofost_toml`; `potatod.crp.toml` populated; salinitystress case authored end-to-end; full parity for cases 1 + 5; tag `rescue/phase-4c-b-wofost`.
- **Phase 4d**: `bottom_boundary_config_t`, `heat_config_t`, `solute_config_t` for the legacy fields that don't yet have homes. Per-case extensions for 1, 4, 5.
- **Phase 4 proper**: wire new TOML path into `swap.f90` execution; retire `readswap.f90` and ttutil.

## Definition of done for this phase

- Tag `rescue/phase-4c-a-crop-fixed-grass` applied locally.
- Every exit criterion above verified.
- `docs/` complete and rendered by FORD.
- All in-scope `.swp.template` / `.dra` / `.crp` files have `.toml` equivalents under `tests/swap-cases/toml/<case>/`.
- Three legacy quirks (date convention, swmacro shadow, nrlevs hardcode) eliminated.
