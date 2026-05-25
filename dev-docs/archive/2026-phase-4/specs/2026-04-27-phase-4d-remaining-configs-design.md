---
title: "Phase 4d — Bottom Boundary, Heat, Irrigation, Solute Configs + Grass Mowing/Grazing"
author: Mateusz Zawadzki
date: 2026-04-27
status: draft
---

# Phase 4d: Bottom Boundary, Heat, Irrigation, Solute Configs + Grass Mowing/Grazing

The final phase of the configuration-modernization arc. After 4d closes, every section of every regression case's legacy `.swp` / `.dra` / `.crp` files maps to a TOML field with a parity test. The TOML pipeline becomes the canonical "what the case wants" representation; the legacy `readswap.f90` runtime path stays in place for now (its retirement is Phase 4 proper).

## Background

Phase 4c-b closed with five non-macropore cases at full crop parity, plus the new `legacy_crop_helper` module that lets parity tests call legacy crop sub-readers directly. But four legacy `.swp`/`.dra` sections remain unmodelled in `swap_config_t`:

- **Bottom boundary** — SWAP's `SWBOTB` switch with seven branches (SWBOTB=1..8). The schema currently captures none. Case 5 (salinitystress) uses SWBOTB=3 (Cauchy with regional aquifer head); the audit (Task 9 of 4c-b) flagged this as the largest 4d gap.
- **Heat transport** — `SWHEA` switch + frost params + soil-temperature initialization. Case 4 (oxygenstress) and case 5 use it.
- **Irrigation** — Two halves: fixed-irrigation at `.swp` level (`SWIRFIX`, `IRGFIL`, fixed-event table) and scheduling at `.crp` level (`SCHEDULE`, `TCS`, `DCS`, threshold/depth tables). Phase 4c-b deferred the `.crp` half explicitly.
- **Solute** — `SWSOLU` switch + decomposition, dispersion, root uptake, salinity diffusion. Case 5 uses solute transport with `SWSALINITY=1` for Maas-Hoffman.

Plus one Phase 4c-a deferral:
- **Grass per-event tables** — `cropgrass_config_t` only handles `swharv=0` paths today. Per-event mowing/grazing tables (`swharv=1` DM-threshold, `swharv=2` fixed-date table, plus grazing schedule) exist in case 2 (grassgrowth) and case 4 (oxygenstress) `.crp` files but are currently authored around with `swharv=0`. 4d completes the type-3 grass coverage.

After 4d, no regression case has any "READING NOTE: deferred to Phase 4d" markers in its parity test.

## Scope

### In

- **Five new config types** in `src/config/` (alphabetical):
  - `bottom_boundary_config_t` (~40 fields across 7 SWBOTB branches).
  - `heat_config_t` (~15 fields: switches, conduct/capacity options, frost params, initial soil-temperature table).
  - `irrigation_config_t` (~30 fields top-level + nested `irrigation_schedule_t` sub-type for per-crop scheduling — Option A per user decision).
  - `solute_config_t` (~25 fields: solute transport switches, dispersion, decomposition, root uptake, salinity).
- **One config-type extension** to existing 4c-a code:
  - `cropgrass_config_t` — extend with mowing-event and grazing-event tables (2D allocatables of dates+amounts, similar shape to AFGEN tables).
- **Five new section readers** in `src/io/toml/`:
  - `read_bottom_boundary_toml.f90`, `read_heat_toml.f90`, `read_irrigation_toml.f90`, `read_solute_toml.f90`, plus extension to `read_cropgrass_toml.f90` for the new tables.
- **Wiring** in `swap_config_t` (`config%bottom_boundary`, `config%heat`, `config%irrigation`, `config%solute`) and `load_swap_config.f90` (call new readers).
- **Per-crop nested irrigation schedule**: `cropwofost_config_t` and `cropfixed_config_t` get a `schedule :: irrigation_schedule_t` field. `cropgrass_config_t` only if a regression case uses irrigation scheduling on grass — verify via audit step.
- **Per-case TOML extensions** under `tests/swap-cases/toml/` for cases that exercise the new configs:
  - All 5 non-macropore cases get `[bottom_boundary]`, `[heat]`, `[irrigation]`, `[solute]` sections. Most will be in the SWBOTB=2 or SWBOTB=8 simple branches; case 5 gets SWBOTB=3 with the regional-aquifer params.
  - Cases 2 (grassgrowth) and 4 (oxygenstress) get extended `grassd.crp.toml` with mowing/grazing tables.
  - Cases that use scheduled irrigation (TBD via audit step) get `[irrigation_schedule]` blocks inside their per-crop TOMLs.
- **Parity-test extensions** to the four 4c-a parity tests + the 4c-b salinitystress test, asserting the new fields against `variables%` globals.
- **Schema doc extension**: `docs/configuration-schema.md` gains four new sections + a grass-mowing extension.
- **Coverage rebaseline + tag** `rescue/phase-4d-remaining-configs`.

### Out

- **Macropore case** (3.macroporeflow) — out for all parity work in 4-x.
- **Wiring new TOML path into `swap.f90` execution** — still legacy-readswap-only at runtime. Phase 4 proper.
- **Nitrogen-use (LINTUL4) section in WOFOST** — none of the regression cases activate it.
- **Python bindings, DLL multi-instance** — out of scope for the rescue.
- **`variables` module retirement** — Phase 4 item 1.

## Critical decisions

### D1 — Irrigation config shape (user-confirmed Option A)

One `irrigation_config_t` at the top level of `swap_config_t` for the `.swp` fixed-irrigation half:
```fortran
type :: irrigation_config_t
   integer                       :: swirfix
   character(len=:), allocatable :: irgfil       ! when swirfix=1
   real(real64),     allocatable :: fixed_events(:,:)  ! Date / depth / conc / type
contains
   procedure :: validate => irrigation_config_validate
   procedure :: finalize => irrigation_config_finalize
end type irrigation_config_t
```

Plus a sub-type `irrigation_schedule_t` nested inside per-crop config types for the `.crp` scheduling half:
```fortran
type :: irrigation_schedule_t
   integer      :: schedule          ! 0/1
   integer      :: startirr_day, startirr_month
   integer      :: endirr_day,   endirr_month
   real(real64) :: cirrs
   integer      :: isuas             ! 0=sprinkling, 1=surface
   integer      :: tcs               ! 1..8 timing-criterion switch
   integer      :: dcs               ! 1=back-to-FC, 2=fixed depth
   real(real64), allocatable :: trel_table(:,:)   ! TCS=1
   real(real64), allocatable :: raw_table(:,:)    ! TCS=2
   real(real64), allocatable :: taw_table(:,:)    ! TCS=3
   real(real64), allocatable :: dwa_table(:,:)    ! TCS=4
   real(real64) :: irgthreshold     ! TCS=6
   real(real64), allocatable :: hcri_table(:,:)   ! TCS=7
   real(real64), allocatable :: tcri_table(:,:)   ! TCS=8
   real(real64) :: dcrit             ! TCS=7 or 8
   integer      :: swcirrthres
   real(real64) :: cirrthres, perirrsurp
   integer      :: tcsfix
   integer      :: irgdayfix
   real(real64) :: phfieldcapacity
   real(real64), allocatable :: di_table(:,:)     ! DCS=1
   real(real64) :: raithreshold
   real(real64), allocatable :: fid_table(:,:)    ! DCS=2
   integer      :: dcslim
   real(real64) :: irgdepmin, irgdepmax
contains
   procedure :: validate => irrigation_schedule_validate
   procedure :: finalize => irrigation_schedule_finalize
end type irrigation_schedule_t
```

Each per-crop config type (`cropfixed_config_t`, `cropwofost_config_t`, `cropgrass_config_t`) gains a `schedule :: irrigation_schedule_t` field. The validator gates schedule fields on `schedule%schedule == 1`.

**Why this matters**: mirrors the legacy split exactly (case 1 hupselbrook has fixed irrigation in `.swp` for one rotation entry and scheduling in `.crp` for another). One TOML file per side, no duplication.

### D2 — Bottom boundary `SWBOTB` branching

`SWBOTB` has 8 legal values:
1. Pressure head as function of time
2. Flux as function of time
3. Cauchy: regional bottom flux ∝ (gw_level - aquifer_level), with aquifer params
4. Free drainage (no inflow)
5. Pressure head specified at fixed depth below surface (no time series)
6. Free outflow at depth = unsaturated zone bottom
7. Constant zero flux
8. Soil moisture content at the bottom is fixed at field capacity

Most regression cases use one of {2, 4, 6, 7, 8}. Case 5 uses 3.

The schema captures all 8 branches in one `bottom_boundary_config_t` with conditional fields:
- Always: `swbotb` enum.
- For SWBOTB=1: `swc_table(:,:)` (date / pressure-head pairs).
- For SWBOTB=2: `qbot_table(:,:)`.
- For SWBOTB=3: `shape`, `hdrain`, `rimlay`, `aqave`, `aqamp`, `aqomeg`, plus an optional `cofqha` table.
- For SWBOTB=5: `hbot`, `rhobot`.
- Other branches: scalar params per legacy spec.

Validator gates table presence on `swbotb` value.

### D3 — Heat transport — `SWHEA` and conduct/capacity options

```fortran
type :: heat_config_t
   integer :: swhea          ! 0/1
   integer :: swcalt         ! 0/1 — calculation type
   integer :: swtopbhea      ! top-boundary type
   integer :: swbotbhea      ! bottom-boundary type
   ! conduction/capacity options:
   real(real64), allocatable :: psand(:)      ! per soil layer
   real(real64), allocatable :: pclay(:)
   real(real64), allocatable :: porg(:)
   ! initial soil-temperature table:
   real(real64), allocatable :: tsoil(:,:)    ! depth / temperature
   ! frost params:
   real(real64) :: tfroststa, tfrostend
contains
   procedure :: validate, finalize
end type heat_config_t
```

When `swhea=0`, none of the dependent fields are required.

### D4 — Solute — `SWSOLU` and the salinity-Maas-Hoffman cross-tie

```fortran
type :: solute_config_t
   integer :: swsolu         ! 0/1
   integer :: swbotbc        ! solute boundary at bottom
   real(real64) :: cdrain    ! drain-water concentration
   real(real64) :: cseep     ! seepage concentration
   real(real64) :: tscf      ! transpiration stream concentration factor
   real(real64) :: ldis      ! dispersion length
   ! root uptake:
   real(real64) :: rtheta    ! water content threshold
   real(real64) :: bexp      ! exponent
   ! salinity (already in cropwofost%salinity, but the .swp side has globals):
   real(real64) :: ecmax, ecslop      ! when crop-side swsalinity unset; cross-checked against per-crop
   integer      :: swsoltyp  ! solute type
   integer      :: swdc      ! decomposition switch
   real(real64), allocatable :: pertabsolu(:,:)  ! decomposition table per layer
contains
   procedure :: validate, finalize
end type solute_config_t
```

The crop-side salinity stress (`config%crop%rotation_wofost(N)%salinity%saltmax` etc.) and the solute-side parameters (`config%solute%cdrain` etc.) are SEPARATE; both load from independent legacy globals. Validators don't cross-check — that would couple unrelated branches.

### D5 — Grass mowing/grazing event tables

`cropgrass_config_t` extension:
```fortran
type, extends(cropgrass_config_t_base) :: cropgrass_config_t
   ! ... existing fields ...
   ! Phase 4d additions:
   integer      :: nmow
   real(real64), allocatable :: mowing_dates(:)        ! day-of-year per event
   real(real64), allocatable :: mowing_heights(:)      ! optional, leave allocated for swdmmow=1
   integer      :: swdmmow                             ! 0=use heights, 1=use DM threshold
   real(real64) :: dmharvest, daylastharvest, dmlastharvest
   integer      :: maxdaymow
   integer      :: nstart_graz, nstop_graz
   integer      :: maxdaygrz
   real(real64) :: dmgrazing
   integer      :: swdmgrz
   real(real64) :: lsdb_default
   ! ... rest unchanged ...
end type
```

Validator: when `swharv=1`, expect `nmow >= 1` and `mowing_dates` allocated; when `swharv=2`, expect a fixed-date table; when `swharv=0`, all of these are optional.

The 4c-a `cropgrass_config_t` already has `nmow`, `mowing_dates`, `swharv` — Phase 4d *extends* the validator branches and adds `swdmmow` / DM-threshold handling.

### D6 — Per-case TOML rollout strategy

The four new sections (`[bottom_boundary]`, `[heat]`, `[irrigation]`, `[solute]`) are added to ALL non-macropore cases simultaneously. For cases that don't exercise a particular section, the section still has a switch (`swhea=0`, `swsolu=0`, `swirfix=0`, `swbotb=4`) but no payload.

This keeps every case's TOML self-contained — no "case 5 has solute, case 1 doesn't" dichotomy. The TOML schema is uniform; the values vary.

### D7 — Audit step before authoring

Before any code lands, an audit task (Task 1) reads each case's `.swp` and relevant `.crp` files to record exactly which `swbotb`, `swhea`, `swsolu`, `swirfix` values each case sets, plus which crops use scheduled irrigation. The audit doc gets committed and drives Tasks 2-N. Mirrors Phase 4c-b Task 9's pattern.

## Critical files

### Create

- `src/config/bottom_boundary_config.f90`
- `src/config/heat_config.f90`
- `src/config/irrigation_config.f90` (top-level + nested schedule sub-type)
- `src/config/solute_config.f90`
- `src/io/toml/read_bottom_boundary_toml.f90`
- `src/io/toml/read_heat_toml.f90`
- `src/io/toml/read_irrigation_toml.f90`
- `src/io/toml/read_solute_toml.f90`
- `tests/unit/config/test_bottom_boundary_config.pf`
- `tests/unit/config/test_heat_config.pf`
- `tests/unit/config/test_irrigation_config.pf`
- `tests/unit/config/test_solute_config.pf`
- `tests/unit/io/toml/test_read_bottom_boundary_toml.pf`
- `tests/unit/io/toml/test_read_heat_toml.pf`
- `tests/unit/io/toml/test_read_irrigation_toml.pf`
- `tests/unit/io/toml/test_read_solute_toml.pf`
- `docs/phase-4d-case-audit.md` (audit doc)

### Modify

- `src/config/swap_config.f90` — add `bottom_boundary`, `heat`, `irrigation`, `solute` fields.
- `src/config/cropgrass_config.f90` — extend with per-event tables.
- `src/config/cropwofost_config.f90` — add `schedule :: irrigation_schedule_t`.
- `src/config/cropfixed_config.f90` — add `schedule :: irrigation_schedule_t`.
- `src/io/toml/read_cropgrass_toml.f90` — extend for new tables.
- `src/io/toml/read_cropwofost_toml.f90` — extend for `[irrigation_schedule]` block.
- `src/io/toml/read_cropfixed_toml.f90` — extend for `[irrigation_schedule]` block.
- `src/io/toml/load_swap_config.f90` — call new top-level readers.
- All five non-macropore parity tests (`test_hupselbrook_parity.pf`, `test_grassgrowth_parity.pf`, `test_oxygenstress_parity.pf`, `test_surfacewater_parity.pf`, `test_salinitystress_parity.pf`) — extend with bottom_boundary/heat/irrigation/solute assertions.
- `tests/unit/meson.build` — register new sources + tests.
- `tests/unit/testSuites.inc` — register new pFUnit suites.
- `meson.build` — register new `src/` files in `sources`.
- `docs/configuration-schema.md` — four new sections + grass extension.

### Submodule edits (`tests/swap-cases/toml/`)

For each of cases 1, 2, 4, 5, 6: extend `swap.toml` with `[bottom_boundary]`, `[heat]`, `[irrigation]`, `[solute]` sections. For cases 2 and 4: extend `grassd.crp.toml` with mowing/grazing tables. For cases that use scheduled irrigation (verify via Task 1 audit): add `[irrigation_schedule]` blocks to per-crop TOMLs.

## Tasks (high-level — full breakdown in plan file)

### Part A — Audit + foundation

1. **Case audit**. Read all five non-macropore cases' `.swp`, `.dra`, and `.crp` files. Record `swbotb`, `swhea`, `swsolu`, `swirfix`, `schedule` values per case. Document in `docs/phase-4d-case-audit.md`. Identify which crops use scheduled irrigation. No code change.

### Part B — Bottom boundary

2. `bottom_boundary_config_t` skeleton.
3. Validators + per-sub-type pFUnit tests.
4. `read_bottom_boundary_toml`.
5. Wire into `swap_config_t` + `load_swap_config`.

### Part C — Heat

6. `heat_config_t` skeleton + validators.
7. `read_heat_toml` + tests.
8. Wire into `swap_config_t`.

### Part D — Irrigation

9. `irrigation_config_t` (top-level) + `irrigation_schedule_t` (nested sub-type) — both in `irrigation_config.f90`.
10. Validators + tests.
11. `read_irrigation_toml` (top-level) + extension to per-crop readers (`read_cropwofost_toml`, `read_cropfixed_toml`, `read_cropgrass_toml`) for the `[irrigation_schedule]` block.
12. Wire into `swap_config_t`.

### Part E — Solute

13. `solute_config_t` + validators + reader + tests.
14. Wire into `swap_config_t`.

### Part F — Grass mowing/grazing extension

15. Extend `cropgrass_config_t` with per-event tables. Validator branches for `swharv=0/1/2`.
16. Extend `read_cropgrass_toml` for the new tables.

### Part G — Per-case TOML rollout

17. For each non-macropore case: extend `swap.toml` with the four new sections. Submodule commits + outer-repo bumps.
18. Cases 2 and 4: extend `grassd.crp.toml` with mowing/grazing tables.
19. Cases with scheduled irrigation (per Task 1 audit): extend per-crop TOMLs with `[irrigation_schedule]`.

### Part H — Parity-test extensions

20. Extend each parity test with assertions for the new sections. Use `read_legacy_*` wrappers from `legacy_crop_helper` where the legacy reader populates `variables%` (which it does for solute, heat, bottom_boundary, irrigation — these are all `readswap`-side, not crop-sub-reader-side).

### Part I — Closeout

21. Schema doc extension.
22. Coverage rebaseline.
23. check-full final log + tag `rescue/phase-4d-remaining-configs`.

## Risks & mitigations

| Risk | Mitigation |
|---|---|
| Legacy `readswap` populates `variables%` for some new fields conditionally on switch state — assertions need to gate on the switch | Mirror Phase 4c-b's READING NOTE pattern; gate parity assertions on `swbotb`, `swhea`, etc. |
| `bottom_boundary_config_t` SWBOTB=3 has the most fields; the others are stubs — validator may bloat | Group SWBOTB=3 fields into a sub-section `[bottom_boundary.cauchy]`; other SWBOTB branches stay flat |
| Case 5 has SWBOTB=3 with regional aquifer head/amplitude/period params — small risk of legacy-reader-vs-TOML mismatch on units | Audit step (Task 1) records exact legacy values; parity test catches mismatches |
| Per-event mowing tables in case 2 may overlap with case 2's existing `swharv=0` workaround in `cropgrass_config_t` | Task 15 explicitly removes the workaround when extending; commit is "extend cropgrass for per-event tables (Phase 4d Task 15)" |
| Irrigation scheduling block inside each per-crop type bloats `cropwofost_config_t` further | Sub-type stays small (`irrigation_schedule_t` is independent); validator branches gate cleanly |
| Per-case TOML rollout (Tasks 17-19) hits the submodule-commit permission issue from 4c-b | Same workaround: subagent writes files only, user batch-commits at end |
| 4c-b parity tests need updating after Tasks 5/8/12/14/16 land — easily missed | Task 20 explicitly walks all five parity tests in one focused pass |
| Coverage stays flat at 53.1% (gcovr issue) | Pre-existing, not a 4d regression — document as before |

## Verification

After each Part: `pixi run -e test test-pfunit` green.
After each Part touching `src/` or main code: `pixi run -e test check-fast` green (4/4 in <90s).
After all per-case TOML rollout (Parts G + H): `pixi run -e test check-full` 6/6 green.

Final closeout (Part I):
```
pixi run -e test check-full
pixi run -e test test-pfunit
pixi run -e coverage coverage-report
git tag rescue/phase-4d-remaining-configs
```

No push.

## File count summary

- Create: ~20 files (5 config_t modules, 5 readers, 10 test files, 1 audit doc).
- Modify: ~12 files (swap_config + 3 crop configs + 3 crop readers + load_swap + 5 parity tests + meson + docs).
- LoC estimate: ~3000-3500 lines added (mostly mechanical: field declarations, validator stanzas, readers, parity assertions).
- Submodule edits: 5 swap.toml extensions + 2 grassd.crp.toml extensions + N per-crop irrigation_schedule additions (TBD via Task 1 audit).

## Phase 4 proper preview (informational; out of 4d scope)

After 4d closes, the only remaining work to retire `readswap.f90` is the runtime path: rewire `swap.f90`'s execution to read from `swap_config_t` instead of `variables` globals. This is Phase 4 item 1 — a substantial refactor that builds on top of all 4a/b/c/d schema work. The rescue arc treats it as a separate phase because it changes physics-path code.
