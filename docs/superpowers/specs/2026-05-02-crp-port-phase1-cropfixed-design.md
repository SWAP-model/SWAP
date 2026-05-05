# `.crp` port — Phase 1 (cropfixed via case 6 surfacewater) — design

**Status:** Draft, awaiting review
**Date:** 2026-05-02
**Predecessors:**
- `docs/superpowers/specs/2026-05-01-swap-dra-port-design.md` (same strangler shape)
- `docs/adr/0015-strangler-narrow-scope-stub-errors.md` (stub-error pattern)
- `docs/adr/0016-per-rotation-crop-config-cache.md` (cache pattern, written alongside this spec)

## Goal

Port case 6 (surfacewater)'s `grass.crp` (type 1, fixed-crop mode) to TOML. The
new executable reads case 6's crop data only from `grass.crp.toml` plus the
typed pipeline; no `.crp` ASCII file is opened on the runtime path for type-1
rotations. Case 6 regression remains 5/5 green throughout. The legacy
`readcropfixed` reader stays alive in `readswap.f90` as a parity-test fixture.

This is **Phase 1 of a 4-phase `.crp` port**:
- **Phase 1** (this spec): cropfixed via case 6.
- **Phase 2**: cropwofost (detailed) via case 5 (salinitystress).
- **Phase 3**: cropgrass via case 4 (oxygenstress) and case 2 (grassgrowth).
- **Phase 4**: case 1 (hupselbrook) — integration test exercising all three
  types in one rotation; legacy fallback in `cropgrowth.f90` is removed here.

Each phase has its own spec, plan, and implementation cycle.

## Non-goals (Phase 1)

- Porting cropwofost or cropgrass — Phases 2/3.
- Porting cases 1, 2, 3, 4, 5 — they continue to run via the legacy reader for
  their respective crop modes.
- Removing legacy `readcropfixed` from `readswap.f90` — kept as parity fixture.
- Fully implementing the runtime semantics for switched-off branches in
  case 6's grass.crp (e.g. SWOXYGEN=2 Bartholomeus, SWDROUGHT=2 De Jong van
  Lier). The schema accepts those values 1:1; the validator rejects them with
  ADR-0015 stub-errors. Case 6 authors all switches at the supported value.

## Status quo

### Schema (partial)

`src/config/cropfixed_config.f90` (~93 lines, Phase 4c-a) covers:
- Phenology: `idev`, `lcc`
- Light: `kdif`, `kdir`
- Crop tables (declared, partially parsed): `cftb`, `chtb`
- Root: `rdi`, `rri`, `rdc`, `rdctb`
- Feddes: `hlim1..hlim4`, `adcrh`, `adcrl`, `rsc`
- Salinity: `ecmax`, `ecslop`
- Interception: `cofab`
- Per-crop irrigation: `schedule` (irrigation_schedule_t)

That's ~17 of legacy `readcropfixed`'s ~50 fields.

### Parser (partial)

`src/io/toml/read_cropfixed_toml.f90` parses the existing schema. Does not
yet handle the missing scalars or the `gctb`/`rdtb` tables.

### Loader (skeleton only)

`src/io/toml/read_crop_toml.f90` parses the `[[crop.rotation]]` array into
parallel arrays `rotation_start`, `rotation_end`, `rotation_file`,
`rotation_type`. The referenced `.crp.toml` files are NOT loaded at this
stage. Per-rotation crop content is invisible to the typed pipeline.

### Adapter (rotation metadata only)

`src/io/toml/config_to_variables.f90` writes `croptype`, `cropstart`,
`cropend`, `cropfil` arrays to module globals. Per-rotation crop fields are
not written — `cropgrowth.f90`'s legacy sub-readers fill them at simulation
time.

### Runtime (legacy reader per rotation)

`src/crop/cropgrowth.f90:368` (`cropfixed(task=1)`) calls `readcropfixed(icrop,
cropfil(icrop), lcc, swhydrlift)` per rotation entry. The reader:
1. Opens `<pathcrop>/<cropfil(icrop)>.crp` from disk.
2. Issues ~50 `rd*` calls reading scalars and tables into ~70 module globals.
3. Builds runtime tables (`cumdens`, root density distribution, etc.) from
   the just-read inputs.

For case 6 (3 rotations, all type=1, all referencing `grass.crp`), the
reader runs 3× during simulation init.

### Test fixtures

- `tests/swap-cases/toml/6.surfacewater/grass.crp` — 202 lines (legacy ASCII).
- `tests/swap-cases/toml/6.surfacewater/grass.crp.toml` — 33 lines, 11 keys
  (skeletal).

## Approach

Mirror the swap.dra port's three-bucket split (validators / finalizers /
runtime-init module), with one new piece: the **per-rotation crop config
cache** (ADR 0016). Loader fans out to `.crp.toml` files at config-load time
and stores the parsed content in `crop_config_t.rotation_cropfixed(:)` (and
parallel arrays for types 2/3 in later phases).

### Unit 1 — Schema 1:1 with legacy `grass.crp`

Extend `cropfixed_config_t` so every field that legacy `readcropfixed` reads
has a corresponding TOML key, regardless of whether the parent switch is on
or off in case 6. The user-stated rationale: pay the schema cost once now so
future phases extending semantic coverage to currently-disabled branches
(SWOXYGEN=2, SWDROUGHT=2, etc.) don't have to extend the schema/parser/adapter
pipeline incrementally per branch.

**Added scalar fields** (Part = section in legacy .crp file):
- Part 0a/b/c: `swprep`, `swsow`, `swgerm`
- Part 0d: `dvsend`, `swharv`
- Part 1 (idev=2 path): `tsumea`, `tsumam`, `tbase`
- Part 3: `swgc`
- Part 4: `swcf`
- Part 10: `swrd`, `swdmi2rd`, `swrdc`
- Part 11: `swoxygen`, `swwrtnonox`
- Part 12: `swdrought`
- Part 13: `swsalinity`
- Part xx: `swcompensate`
- Part 14: `swinter`
- Part 15 (irrigation scheduling): `schedule` (already present)

Plus the deeper-detail fields under SWOXYGEN=2 / SWDROUGHT=2 paths
(`criterhr`, `stephr`, `kroot`, `rxylem`, `kstem`, `swrootradius`,
`root_radiusO2`, `q10_root`, `q10_microbial`, `dry_mat_cont_roots`,
`air_filled_root_por`, `spec_weight_root_tissue`, `var_a`, `f_senes`,
`oxygenslope`, `oxygenintercept`, `swoxygentype`, `wiltpoint`, `taccur`,
`alphacrit`, `dcritrtz`, `aeratecrit`, etc.). Schema accepts them with
defaults of 0.0; validator skips range-checks when the parent switch is at
the disabled value.

**Added flat-array tables** (1-D `real(real64)` arrays of `(dvs, value)`
pairs, mirroring legacy `gctb` storage):
- `gctb` — LAI/SCF table (when swgc=1 or 2)
- `cftb` — crop factor table (when swcf=1, already declared, extend parser)
- `chtb` — crop height table (when swcf=2, already declared, extend parser)
- `cfeictb` — wet-crop factor table (only if swcf=3 → stub-errored)
- `rdtb` — root depth vs DVS (when swrd=1)
- `rdctb` — root density distribution (always required)

**Stub-errored switch values** (validator rejects with ADR-0015
`ERR_VALIDATION_CROSS_FIELD`):
- `swoxygen=2` (Bartholomeus) — case 6 has `swoxygen=0`
- `swdrought=2` (De Jong van Lier) — case 6 has `swdrought=1`
- `swcompensate ∈ {1, 2}` (Jarvis, Walsum) — case 6 has `swcompensate=0`
- `swcf=3` (wet-crop, would activate cfeictb) — case 6 has `swcf=1`
- `swharv=1` (DVS-based timing) — case 6 has `swharv=0`
- `schedule=1` (per-crop irrigation scheduling) — case 6 has `schedule=0`

The active subset for case 6 is: SWPREP=0, SWSOW=0, SWGERM=0, DVSEND=2.0,
SWHARV=0, IDEV=1 (LCC=366), KDIF=0.75, KDIR=0.75, SWGC=1, SWCF=1, SWRD=1,
SWOXYGEN=0, SWWRTNONOX=0, SWDROUGHT=1, SWSALINITY=0, SWCOMPENSATE=0,
SWINTER=1, SCHEDULE=0. All supported.

### Unit 2 — Parser extension

`src/io/toml/read_cropfixed_toml.f90` extended to read the ~30 new scalars
and the 6 tables. New table-reading helper if existing helpers don't cover
1-D `(dvs, value)` paired arrays cleanly.

New file-based wrapper:

```fortran
subroutine read_cropfixed_file_toml(path, config, errors, base_path)
   ! Resolves <base_path>/<path>, opens via toml_load, calls existing
   ! read_cropfixed_toml(doc, config, errors).
end subroutine
```

### Unit 3 — Loader extension

`src/io/toml/read_crop_toml.f90`'s rotation-loop dispatcher allocates and
populates the new cache.

```fortran
type :: crop_config_t
   integer                                    :: swcrop = 0
   real(real64), allocatable                  :: rotation_start(:), rotation_end(:)
   integer,      allocatable                  :: rotation_type(:)
   character(len=:), allocatable              :: rotation_file(:)
   ! Phase 1 addition (parallel arrays per type for Phases 2/3):
   type(cropfixed_config_t), allocatable      :: rotation_cropfixed(:)
end type
```

Loader walks `[[crop.rotation]]` once. Per entry:
- Always populate the metadata arrays.
- Dispatch on `rotation_type(i)`:
  - `1` → call `read_cropfixed_file_toml(rotation_file(i), rotation_cropfixed(i), errors, base_path=case_dir)`
  - `2` → no-op in Phase 1; Phase 2 fills `rotation_cropwofost(i)`.
  - `3` → no-op in Phase 1; Phase 3 fills `rotation_cropgrass(i)`.

Slots not populated by the active phase remain at their default
unallocated/zero-initialized state. The `populated` sentinel (Section
"Sentinel" below) gives the runtime an unambiguous check.

**Path resolution:** `rotation_file(i)` is a basename (e.g. `"grass.crp.toml"`).
Resolved relative to the case working directory (the dir containing
`swap.toml`), matching the CSV companion file convention.

**Error handling:** missing or unparseable `.crp.toml` appends a fatal to
`errors` with context `crop.rotation[i].file = "<path>"`. The loader continues
to surface multiple errors per ADR 0008's collected-error pattern.

### Unit 4 — Sentinel for "this slot is config-driven"

Add to `cropfixed_config_t`:

```fortran
logical :: populated = .false.   ! true once read_cropfixed_toml populates
```

Set to `.true.` at the end of `read_cropfixed_toml` (the in-place reader, not
the file-based wrapper) so it covers both the load-from-file path and any
direct unit-test fixture that uses the in-place reader.

`cropgrowth.f90` checks `crop_config_global%rotation_cropfixed(icrop)%populated`
to dispatch between the new and legacy paths.

### Unit 5 — Runtime init module: `cropfixed_init`

New file: `src/crop/cropfixed_init.f90`. Public sub:

```fortran
subroutine cropfixed_init_from_config(cfg, icrop)
   class(cropfixed_config_t), intent(in) :: cfg
   integer,                   intent(in) :: icrop
end subroutine
```

Two halves:

1. **Config → globals copy.** One-for-one mirror of legacy `readcropfixed`'s
   `rd*` calls. Same defensive control flow (e.g. `if (idev == 1) ... else ...`),
   same module globals targeted (`idev`, `kdif`, `gctb`, `cftb`, `hlim1..hlim4`,
   `rdtb`, `rdctb`, `swoxygen`, `swdrought`, etc.).

2. **Runtime init math.** Verbatim port of legacy `readcropfixed`'s tail:
   - Build `cumdens(:)` from the `rdctb` table (~50 lines of array math).
   - Initialize `dvs`, `tsum`, `daycrop`, `nofd` per-rotation runtime state.
   - Any unit normalizations the legacy reader applied inline.

Defense-in-depth runtime guards (matching the validator stub-errors):

```fortran
if (swoxygen == 2 .or. swdrought == 2 .or. swcompensate /= 0 .or. &
    swcf == 3 .or. swharv == 1 .or. cfg%schedule%schedule == 1) then
   call fatalerr_collected('cropfixed_init', &
      'Unsupported runtime branch reached on the TOML path. The validator ' // &
      'should have caught this earlier.')
end if
```

### Unit 6 — Module-level config reference

New module: `src/crop/crop_config_global.f90`. Provides one public variable:

```fortran
type(crop_config_t), pointer :: crop_config_global => null()
```

`config_to_variables`, at the end of its crop block, sets
`crop_config_global => config%crop`. `cropgrowth.f90` reads through this
pointer to access the per-rotation cache.

**Transitional element.** See "Teardown" below.

### Unit 7 — Wiring change in `cropgrowth.f90:368`

Before:
```fortran
call readcropfixed (icrop, cropfil(icrop), lcc, swhydrlift)
```

After:
```fortran
if (associated(crop_config_global) .and. &
    allocated(crop_config_global%rotation_cropfixed) .and. &
    crop_config_global%rotation_cropfixed(icrop)%populated) then
   call cropfixed_init_from_config(crop_config_global%rotation_cropfixed(icrop), icrop)
else
   call readcropfixed (icrop, cropfil(icrop), lcc, swhydrlift)   ! transitional
end if
```

The `else` branch is **transitional** (see Teardown).

## Teardown plan (strangler-fig hygiene)

Per the project's strangler-fig discipline, every transitional element this
phase introduces must have a documented teardown trigger and replacement.

| Element | Introduced for | Teardown trigger | What gets removed | Replacement |
|---|---|---|---|---|
| `else call readcropfixed(...)` fallback in `cropgrowth.f90:368` | Phase 1 ships before Phases 2/3; type-2/3 rotations still need the legacy reader during the multi-phase transition | Phase 4 (hupselbrook) lands all three types | The `if/else` dispatch around all three legacy reader calls (`readcropfixed`, `readwofost`, `readgrass`) | Unconditional `*_init_from_config` calls; the `if (associated(...))` guard becomes unnecessary because every rotation has a populated cache slot |
| Same fallback pattern around `readwofost` (Phase 2) and `readgrass` (Phase 3) | Same | Phase 4 | Same | Same |
| Legacy `readcropfixed`/`readwofost`/`readgrass` subroutines in `readswap.f90` (~2200 lines) | Parity-test fixtures (per ADR 0015) | A future cleanup task that rewrites parity tests to compare against fixture values rather than legacy-reader output | The three subroutines | Fixture-based parity tests using stored expected values |
| `crop_config_global` module-level pointer in `src/crop/crop_config_global_mod` | Bridge from typed config (held in `core/swap.f90` local scope) to legacy runtime subs that don't yet take config as an argument | ADR 0016's "config + state passing" direction lands (post-Phase 4, separate spec) | The `crop_config_global` module entirely | `cropfixed_init_from_config` and downstream computation subs take `config` and `state` as explicit arguments; threaded from `core/swap.f90` |
| `populated :: logical` sentinel on `cropfixed_config_t` (and the same field on `cropwofost_config_t`/`cropgrass_config_t` in Phases 2/3) | Runtime needs an unambiguous "this slot is config-driven, not a default-initialized placeholder" check while the legacy fallback is still reachable | Same as `crop_config_global` (config-passing direction) | The `populated` field | Dispatch happens at the call site based on `rotation_type(icrop)`; no runtime sentinel needed because there is no fallback |
| `rotation_file(:)` legacy `.crp` files left on disk in `tests/swap-cases/<N>.<case>/` (the legacy case dir) | Legacy executable still reads them | Legacy executable retirement (out of current scope; tracked under "supported subset" docs) | The legacy case dir contents | None — the legacy executable retirement is an independent decision |

The `tests/swap-cases/toml/6.surfacewater/grass.crp` file (in the TOML case
dir) is **not transitional** — Phase 1 deletes it as Task 11. The legacy
copy at `tests/swap-cases/6.surfacewater/grass.crp` stays for the legacy
executable and the surfacewater parity test fixture.

## Test plan

### New unit tests

- `tests/unit/config/test_cropfixed_config.pf` (extend) — range checks for
  the ~30 new scalar fields; table-shape checks for `gctb`/`cftb`/`chtb`/
  `rdtb`/`rdctb`; stub-error rejection of `swoxygen=2`, `swdrought=2`,
  `swcompensate ∈ {1,2}`, `swcf=3`, `swharv=1`, `schedule=1`.
- `tests/unit/io/toml/test_read_cropfixed_toml.pf` (extend) — round-trip a
  fixture covering all new fields and tables.
- `tests/unit/io/toml/test_load_swap_config.pf` (extend) — assert that loading
  case 6's `swap.toml` populates `rotation_cropfixed(1..3)` with
  `populated=.true.`, all three pointing at content from the same parsed
  `grass.crp.toml`.
- `tests/unit/crop/test_cropfixed_init.pf` (new dir + file) — assert that
  `cropfixed_init_from_config` writes the right values to `variables`
  module globals (idev, kdif, hlim*, gctb, cumdens, etc.). Hand-computed
  `cumdens` assertion for case 6's rdctb (similar to the sttab assertion
  in the swap.dra port).

### Parity test

Extend or add `tests/unit/io/toml/test_surfacewater_crop_parity.pf` (or
fold into existing `test_surfacewater_parity.pf`) asserting that
`cropfixed_init_from_config(rotation_cropfixed(i))` produces the same
module-global state as `read_legacy_cropfixed` on `grass.crp` for each
of the 3 rotations.

### Regression

Case 6 (surfacewater) must remain green. Smoke test (Task 10) renames
`grass.crp` to `grass.crp.disabled` and re-runs regression to confirm the
runtime path is fully decoupled.

### Acceptance gate

- `pixi run test-pfunit` → all green.
- `pixi run regression` → 5/5 cases green.
- `tests/swap-cases/toml/6.surfacewater/grass.crp` → does not exist.
- `git grep "call readcropfixed" src/` → returns the legacy fallback in
  `cropgrowth.f90:368` only (the `else` branch). Other reachability gone.
- `git grep "call readcropfixed" tests/` → still returns parity-test
  references via `legacy_crop_helper`.
- `crop_config_global` exists, populated by `config_to_variables`, read by
  `cropgrowth.f90`.

## Implementation tasks (working draft)

To be elaborated in `writing-plans`:

1. Audit `readcropfixed` (readswap.f90:2037-2519) line-by-line into
   READ/VALIDATE/NORMALIZE/RUNTIME/GUARDED buckets. Docs only.
2. Extend `cropfixed_config_t` schema 1:1 with grass.crp. Add stub-error
   validators. Add `populated` sentinel.
3. Extend `read_cropfixed_toml.f90` parser for all new fields + tables. Set
   `populated = .true.` at end.
4. Add `read_cropfixed_file_toml` wrapper for file-based loading.
5. Extend `read_crop_toml.f90` to dispatch on `rotation_type` and populate
   `rotation_cropfixed(:)` per type-1 entry.
6. Author case 6's `grass.crp.toml` 1:1 with legacy values. Submodule pair
   commit.
7. Create `src/crop/cropfixed_init.f90` with `cropfixed_init_from_config`
   (config-to-globals copy + runtime init math + defense-in-depth guards).
8. Create `src/crop/crop_config_global_mod` module-level pointer; wire
   `config_to_variables` to set it.
9. Wire `cropgrowth.f90:368` to dispatch on the `populated` sentinel; legacy
   fallback in the `else` branch.
10. Smoke test: rename `grass.crp` → `grass.crp.disabled`, run regression,
    confirm 5/5, restore.
11. Delete `grass.crp` from the case-6 TOML dir. Submodule pair commit.
12. Write/extend ADR 0016 (cache pattern + config-passing direction +
    teardown plan).
13. Update `docs/csv-companion-files.md`'s path-resolution section (note the
    `.crp.toml` resolution rule). Update `docs/configuration-schema.md` if
    its `[crop]` / `[[crop.rotation]]` table is missing the cropfixed-fields
    documentation.

## Risk register

- **Hidden invariants in `readcropfixed`'s tail-init math.** The `cumdens`
  build and root-density-distribution interpolation are subtle. Risk
  mitigation: line-by-line audit (Task 1) before extracting to
  `cropfixed_init`, plus a hand-computed `cumdens` assertion in unit tests.
- **Case 6's grass.crp is the simplest cropfixed case** (all-disabled
  switches). The runtime branches that exercise the disabled paths are not
  tested here. Phases 2/3 (different cases, different active branches) will
  surface invariants this phase doesn't. Mitigation: when a future case
  authors a switch we currently stub-error, the work to lift the stub is
  focused and well-bounded by the schema being already in place.
- **`crop_config_global` ordering risk.** If something tries to read
  `crop_config_global` before `config_to_variables` runs (or after the
  config goes out of scope), the pointer is null/dangling. Mitigation:
  document the lifecycle in the new module's header comment; the only
  reader (`cropgrowth.f90`) checks `associated(...)` defensively. The
  `core/swap.f90` flow guarantees the config outlives the simulation
  (it's a local in the same enclosing scope as `SurfaceWater(1)`).
