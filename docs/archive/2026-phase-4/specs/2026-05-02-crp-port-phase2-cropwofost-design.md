# `.crp` port — Phase 2 (cropwofost via case 5 salinitystress) — design

**Status:** Draft, awaiting review
**Date:** 2026-05-02
**Predecessors:**
- `docs/superpowers/specs/2026-05-01-swap-dra-port-design.md` (strangler-fig shape)
- `docs/adr/0015-strangler-narrow-scope-stub-errors.md` (stub-error pattern)
- `docs/adr/0016-per-rotation-crop-config-cache.md` (cache pattern, per-rotation arrays)
- `docs/superpowers/specs/2026-05-02-crp-port-phase1-cropfixed-design.md` (Phase 1; mirror structure here)

## Goal

Port case 5 (salinitystress)'s `potatod.crp` (type 2, WOFOST detailed-crop mode)
to TOML. The new executable reads case 5's crop data only from
`potatod.crp.toml` plus the typed pipeline; no `.crp` ASCII file is opened
on the runtime path for type-2 rotations. Case 5 regression remains 5/5 green
throughout. The legacy `readwofost` reader stays alive in `readswap.f90` as a
parity-test fixture.

This is **Phase 2 of a 4-phase `.crp` port**:
- **Phase 1** (separate spec/plan): cropfixed via case 6 surfacewater.
- **Phase 2** (this spec): cropwofost via case 5 salinitystress.
- **Phase 3**: cropgrass via case 4 (oxygenstress) and case 2 (grassgrowth).
- **Phase 4**: case 1 (hupselbrook) — integration test exercising all three
  types in one rotation; legacy fallback in `cropgrowth.f90` is removed here.

## Non-goals (Phase 2)

- Porting cropfixed or cropgrass — Phases 1/3.
- Porting cases 1, 2, 3, 4 — they continue running via the legacy reader for
  their respective crop modes.
- Removing legacy `readwofost` from `readswap.f90` — kept as parity fixture.
- Fully implementing runtime semantics for every branch in `readwofost`
  (soybean variant `swsoybean=1`, bulb crop `swbulb=1`, drought De Jong van
  Lier `swdrought=2`, oxygen Bartholomeus `swoxygen=2`, CO2 correction
  `swco2=1`, irrigation scheduling `schedule=1`, N-P-K nutrient block
  `flCropNut=.true.`, `swinter=2/3`). Per ADR 0015 these are stub-errored at
  validate time. Case 5 authors all switches at the supported value.
- Adding `crop_config_global_mod` — that module is introduced by Phase 1 and
  reused here without modification.

## Status quo

### Schema (comprehensive — already at Phase 4c-b level)

`src/config/cropwofost_config.f90` (866 lines) is the most complete of the
three crop-mode schemas. It covers all 21 sub-types of `cropwofost_config_t`:
`wofost_preparation_t`, `wofost_sowing_t`, `wofost_germination_t`,
`wofost_harvest_t`, `wofost_cropfactor_t`, `wofost_phenology_t`,
`wofost_initial_t`, `wofost_greenarea_t`, `wofost_assimilation_t`,
`wofost_conversion_t`, `wofost_respiration_t`, `wofost_partitioning_t`,
`wofost_death_t`, `wofost_root_t`, `wofost_oxygen_stress_t`,
`wofost_drought_stress_t`, `wofost_salinity_t`, `wofost_compensate_t`,
`wofost_interception_t`, `wofost_co2_t`, `wofost_management_t`, plus
`irrigation_schedule_t`. Validators and finalizers are stubs (return early).

**Schema gap vs `readwofost`:** The schema is already 1:1 with the `readwofost`
field list for the supported branches. The key finding from the audit (Task 1
of the implementation plan) will be to confirm:
- No missing fields for the case-5-active branches.
- The soybean-variant fields (`mg`, `dvsi`, `dvrmax1`, `dvrmax2`, `tmaxdvr`,
  `tmindvr`, `toptdvr`, `flrfphotoveg`, `flphenodayl`, `popt`, `pcrt`) are
  absent from the schema — they should be added for schema 1:1 completeness,
  then stub-errored because `swsoybean=0` in all current test cases.
- The bulb-crop fields (`swbulb`, `fbltb`, `pld`, `plwti`, `remoc`) are absent
  from the schema — same treatment.
- The vernalization table `vernrtb` (used when `idsl=2`) is stored under
  `phenology` as `verntb` in the schema — verify name reconciliation vs the
  legacy `vernrtb` variable.
- The N-P-K nutrient block fields (`LRNR`, `LSNR`, `NLAI`, etc.) are read
  from the same `.crp` file in the legacy reader's task=1 block via a second
  `rdinit` call (lines ~994-1021). These are not in the schema. They must be
  added to `cropwofost_config_t` for schema 1:1, with a stub-error when
  `flCropNut=.true.` (not yet supported in the TOML pipeline; case 5 has
  `flCropNut=.false.`).
- The `swrdc` field (switch for root density development, default 0 in
  `readwofost`) is absent from the schema — add as an integer field.

### Schema-gap summary (quantitative)

| Gap | Fields | Case 5 value | Action |
|---|---|---|---|
| Soybean variant | `swsoybean`, `mg`, `dvsi`, `dvrmax1`, `dvrmax2`, `tmaxdvr`, `tmindvr`, `toptdvr`, `flrfphotoveg`, `flphenodayl`, `popt`, `pcrt` | `swsoybean=0` (absent key) | Add to schema; `swsoybean=1` → stub-error |
| Bulb crops | `swbulb`, `fbltb(:,:)`, `pld`, `plwti`, `remoc` | `swbulb=0` (absent key) | Add to schema; `swbulb=1` → stub-error |
| Root density switch | `swrdc` | `0` (hard-coded default in legacy) | Add to schema; schema accepts `0` or `1`; `swrdc=1` stub-error deferred |
| N-P-K nutrient | ~18 nutrient scalars + `nmxlv(:)` table | `flCropNut=.false.` | Add nutrient sub-type to schema; `flCropNut=.true.` → stub-error |
| Irrigation schedule | `schedule` top-level switch | `schedule=0` | Already in schema via `irrigation_schedule_t`; `schedule=1` → stub-error |
| `swco2=1` CO2 tables | `co2amaxtb`, `co2efftb`, `co2tratb` + `atmofil` | `swco2=0` | Already in schema; `swco2=1` → stub-error |
| `swdrought=2` | `wiltpoint`, `kstem`, `rxylem`, `rootradius`, `kroot`, `rootcoefa`, `rooteff`, `stephr`, `criterhr`, `taccur`, `swhydrlift` | `swdrought=1` | Already in schema (under `drought_stress`); `swdrought=2` → stub-error |
| `swoxygen=2` | Bartholomeus fields; `swoxygentype`, `swtopsub`, `nrstaring` | `swoxygen=1` | Already in schema (under `oxygen_stress`); `swoxygen=2` → stub-error |
| `swinter=2/3` | Gash / storage-cap tables | `swinter=1` | Schema has `swinter ∈ {0,1,2}`; `swinter=2` → stub-error |
| `swharv=1` | DVS-based harvest timing | `swharv=0` | Already in schema; `swharv=1` → stub-error |
| `swcompensate=1/2` | `alphacrit`, `dcritrtz` | `swcompensate=0` | Already in schema; `swcompensate ∈ {1,2}` → stub-error |

### Parser (comprehensive — already at Phase 4c-b level)

`src/io/toml/read_cropwofost_toml.f90` (302 lines) reads all 21 sub-sections
via `get_table` + scalar helpers + `read_table_2d`. The coverage matches the
schema: every field present in `cropwofost_config_t` is parsed. Parser gaps
mirror the schema gaps above (missing soybean, bulb-crop, `swrdc`, nutrient
fields — all to be added in Phase 2 Task 2).

The `irrigation_schedule` section is already read via
`read_irrigation_schedule_from_section`. The parser sets no `populated`
sentinel; that is added in Phase 2 Task 2.

### Loader (already wired for type=2)

`src/io/toml/read_crop_toml.f90:116-117` already dispatches `rotation_type=2`
to `read_cropwofost_toml` and sets `rotation_loaded(i) = .true.`. The
`rotation_wofost(:)` array on `crop_config_t` is allocated at line 57 and
populated at line 116. **No changes required to the loader for Phase 2.**

This is a key difference from Phase 1 (which had to build the entire loader
infrastructure): the loader wiring for type=2 is already in place.

### Adapter (rotation metadata only)

`config_to_variables` writes `croptype`, `cropstart`, `cropend`, `cropfil`
arrays. Per-rotation WOFOST fields are not yet written — `cropgrowth.f90`'s
legacy `readwofost` fills them at simulation time. The `crop_config_global`
pointer will be set at the end of `config_to_variables`'s crop block —
Phase 1 introduces this wiring, Phase 2 reuses it unchanged.

### Runtime (legacy reader per rotation)

`src/crop/cropgrowth.f90:984` (`wofost(task=1)`) calls
```fortran
call readwofost(icrop, cropfil(icrop), swhydrlift, swsoybean, mg, dvsi,
                dvrmax1, dvrmax2, flrfphotoveg, tmaxdvr, tmindvr, toptdvr,
                popt, pcrt, flphenodayl, FraDeceasedLvToSoil)
```
For case 5 (4 rotations, all type=2, all referencing `potatod.crp`), the
reader runs 4× during simulation init.

The `readwofost` signature carries 16 arguments — noticeably more complex than
Phase 1's `readcropfixed(icrop, cropfil, lcc, swhydrlift)` (4 arguments). The
additional arguments are soybean-variant parameters passed back by reference.
In `cropwofost_init_from_config`, the soybean-variant path is stub-errored, so
these output arguments can be set to safe defaults.

### Test fixtures

- `tests/swap-cases/toml/5.salinitystress/potatod.crp` — 418 lines (legacy ASCII).
- `tests/swap-cases/toml/5.salinitystress/potatod.crp.toml` — 216 lines (partial).
  Missing: `[irrigation_schedule]` section (case 5 has `schedule=0`, so the
  parser's default is fine; but the section should be present for completeness).

An existing parity test (`test_salinitystress_parity_potatod_wofost`) already
covers the schema round-trip comparison between `readwofost` output and
`rotation_wofost(1)` config fields. Phase 2 extends this test with assertions
covering the new `populated` sentinel and the `cropwofost_init_from_config`
output.

## Qualitative differences from Phase 1 (cropfixed)

1. **Schema and parser are already comprehensive.** Phase 1 had to extend
   `cropfixed_config_t` by ~30 scalar fields and 6 tables. Phase 2's schema
   is already at 866 lines covering all major sections. The Task 1 audit may
   confirm that only soybean, bulb-crop, `swrdc`, and nutrient fields are
   missing — a much narrower gap than Phase 1.

2. **The loader dispatch is already wired.** Phase 1 had to build `read_crop_toml.f90`'s
   type-1 branch from scratch. Phase 2 inherits the working type-2 branch at
   `read_crop_toml.f90:116-117`.

3. **`crop_config_global_mod` is reused from Phase 1.** Phase 1 creates the
   module. Phase 2 adds a `use` statement to `cropgrowth.f90`'s `wofost`
   subroutine and adds the dispatch branch. No new module needed.

4. **`readwofost` takes 16 arguments** (vs `readcropfixed`'s 4). The dispatch
   wiring in `wofost(task=1)` must preserve the existing `call readwofost`
   in the `else` branch with all 16 arguments. The soybean-variant output
   arguments (`mg`, `dvsi`, `dvrmax1`, `dvrmax2`, `tmaxdvr`, `tmindvr`,
   `toptdvr`, `flrfphotoveg`, `flphenodayl`, `popt`, `pcrt`) must be set to
   safe-default values by `cropwofost_init_from_config` even though the soybean
   branch is stub-errored — they are used later in the `wofost(task=2+)` loop.

5. **Runtime init is qualitatively more complex.** `readwofost`'s tail includes:
   - The `cumdens` computation (same as Phase 1, but conditioned on `swdrought=1`).
   - DVS/tsum/crop state initialization (same structure).
   - An N-P-K nutrient sub-reader block (gated on `flCropNut`; stub-errored).
   - A CO2 file reader block (gated on `swco2=1`; stub-errored).
   - A `.END`-file restart block (gated on `swinco=3`; case 5 uses `swinco=3`
     so this block IS exercised — see Risk register).
   - An irrigation scheduling init call (gated on `schedule=1`; stub-errored).

6. **The `.END`-file restart is exercised by case 5.** Case 5 has `swinco=3`
   (read initial conditions from a previous-run state file). Inside
   `readwofost`, the restart block (lines 3123-3227 of `readswap.f90`) reads
   crop state from an `.END` file only when `t1900 - tstart < 1e-3` and
   `swinco=3` and `|t1900 - cropstart(icrop)| >= 1e-3`. This condition means
   the block fires for the first time-step of the simulation if the crop season
   started before the simulation start. In practice for case 5's 2012 season
   (crop starts 2012-04-17, simulation starts 2012-01-01), the crop start is
   after the simulation start, so the condition evaluates to false and the block
   is skipped. Nevertheless, `cropwofost_init_from_config` must correctly
   replicate the conditional skip, with a defense-in-depth stub-error if the
   block would fire.

## Approach

Mirror the Phase 1 three-bucket split (validators / runtime-init module /
dispatch wiring), adapted to cropwofost's larger schema and more complex
`readwofost` argument list.

### Unit 1 — Schema gap-fill (soybean, bulb-crop, swrdc, nutrient, populated)

Extend `cropwofost_config_t` with:
- A new sub-type `wofost_soybean_t` on the top-level (or flat fields under
  phenology); adds `swsoybean`, `mg`, `dvsi`, `dvrmax1`, `dvrmax2`, `tmaxdvr`,
  `tmindvr`, `toptdvr`, `flrfphotoveg :: logical`, `flphenodayl :: logical`,
  `popt`, `pcrt`. The validator stub-errors `swsoybean=1`.
- A new sub-type `wofost_bulb_t` on the top-level; adds `swbulb`, `fbltb(:,:)`,
  `pld`, `plwti`, `remoc`. The validator stub-errors `swbulb=1`.
- Field `swrdc :: integer = 0` on `wofost_root_t`. Accepted value is 0 only
  for now; `swrdc=1` → stub-error.
- A new sub-type `wofost_nutrient_t` on the top-level; adds `flcropnut ::
  logical = .false.` as the switch, plus the ~18 nutrient scalars and
  `nmxlv(30)` table. The validator stub-errors `flcropnut = .true.`.
- `populated :: logical = .false.` on `cropwofost_config_t`. Set to `.true.`
  at the end of `read_cropwofost_toml`.

**Stub-errored switch values** (validator rejects with ADR-0015
`ERR_VALIDATION_CROSS_FIELD`):
- `swsoybean=1` — case 5 has `swsoybean=0` (absent key, defaults to 0)
- `swbulb=1` — case 5 has `swbulb=0` (absent key)
- `swrdc=1` — case 5 has `swrdc=0` (absent key)
- `flcropnut=.true.` — case 5 has `flCropNut=.false.`
- `schedule=1` (irrigation scheduling) — case 5 has `schedule=0`
- `swco2=1` — case 5 has `swco2=0`
- `swdrought=2` — case 5 has `swdrought=1`
- `swoxygen=2` — case 5 has `swoxygen=1`
- `swinter=2` — case 5 has `swinter=1`
- `swcompensate ∈ {1,2}` — case 5 has `swcompensate=0`
- `swharv=1` — case 5 has `swharv=0`
- `swsalinity=2` (osmotic head) — case 5 has `swsalinity=1` (Maas-Hoffman,
  which IS supported)

The active subset for case 5 is: SWPREP=0, SWSOW=0, SWGERM=0, DVSEND=3.0,
SWHARV=0, SWCF=2 (chtb + albedo + rsc + rsw), IDSL=0, swsoybean=0,
TSUMEA=150, TSUMAM=1550, TDWI=75, LAIEM=0.0589, RGRLAI=0.012, SPA=0, SSA=0,
SPAN=37, TBASE=2, KDIF=1.0, KDIR=0.75, EFF=0.45, CVL=0.72, CVO=0.85,
CVR=0.72, CVS=0.69, Q10=2.0, RML=0.03, RMO=0.0045, RMR=0.01, RMS=0.015,
PERDL=0.03, SWRD=2 (rdi+rri+rdc+swdmi2rd), SWDMI2RD=1, RDCTB (2-row),
SWOXYGEN=1 (Feddes, hlim1/hlim2u/hlim2l), SWWRTNONOX=1, AERATECRIT=0.5,
SWDROUGHT=1 (hlim3h/hlim3l/hlim4/adcrh/adcrl), SWSALINITY=1 (saltmax+saltslope),
SWCOMPENSATE=0, SWINTER=1 (cofab=0.25), SWCO2=0, FraHarLosOrm_lv/st/so,
FraDeceasedLvToSoil=0.3, SWPOTRELMF=2, RELMF=0.8, SCHEDULE=0. All supported.

### Unit 2 — Parser gap-fill and `populated` sentinel

Extend `read_cropwofost_toml.f90` to parse the new sub-sections (`soybean`,
`bulb`, `root.swrdc`, `nutrient`) and set `config%populated = .true.` at
the very end of `read_cropwofost_toml`. The `irrigation_schedule` section
is already parsed; no changes needed there.

Complete `potatod.crp.toml` by adding the `[irrigation_schedule]` section
(case 5 `schedule=0`; the section was absent).

### Unit 3 — Loader: no changes required

`read_crop_toml.f90` already populates `rotation_wofost(i)` and sets
`rotation_loaded(i) = .true.` for type=2 entries. The `populated` sentinel
on the sub-config is set by the parser, not the loader. **No changes needed.**

### Unit 4 — Runtime init module: `cropwofost_init`

New file: `src/crop/cropwofost_init.f90`. Public sub:

```fortran
subroutine cropwofost_init_from_config(cfg, icrop)
   use crop_config_global, only: crop_config_global
   class(cropwofost_config_t), intent(in) :: cfg
   integer,                    intent(in) :: icrop
end subroutine
```

Two halves (same structure as Phase 1's `cropfixed_init_from_config`):

1. **Config → globals copy.** One-for-one mirror of `readwofost`'s `rd*` call
   sequence for the supported subset. Targets the same `variables`-module
   globals: `swcf`, `chtb`, `cftb`, `albedo`, `rsc`, `rsw`, `idsl`, `tsumea`,
   `tsumam`, `dtsmtb`, `tdwi`, `laiem`, `rgrlai`, `slatb`, `spa`, `ssa`,
   `span`, `tbase`, `kdif`, `kdir`, `eff`, `amaxtb`, `tmpftb`, `tmnftb`,
   `cvl`, `cvo`, `cvr`, `cvs`, `q10`, `rml`, `rmo`, `rmr`, `rms`, `rfsetb`,
   `frtb`, `fltb`, `fstb`, `fotb`, `perdl`, `rdrrtb`, `rdrstb`, `swoxygen`,
   `swwrtnonox`, `aeratecrit`, `hlim1`, `hlim2u`, `hlim2l`, `swdrought`,
   `hlim3h`, `hlim3l`, `hlim4`, `adcrh`, `adcrl`, `swsalinity`, `saltmax`,
   `saltslope`, `salthead`, `swcompensate`, `swstressor`, `swinter`, `cofab`,
   `swrd`, `rdi`, `rri`, `rdc`, `swdmi2rd`, `rdctb`, `dvsend`, `swharv`,
   `swrdc`, `schedule`, `relmf`, `swpotrelmf`, `FraDeceasedLvToSoil`,
   `fraharlosorm_lv`, `fraharlosorm_st`, `fraharlosorm_so`, `swco2`, `flco2`.
   Soybean output variables (`mg`, `dvsi`, etc.) set to safe defaults.

2. **Runtime init math.** Verbatim port of `readwofost`'s tail (~lines 3085-3225):
   - `cumdens` computation from `rdctb` (same as Phase 1, active when `swdrought=1`).
   - DVS/tsum/crop state initialization.
   - `.END`-file restart conditional (replicated, but stub-errored if the
     block would fire for an unexpected configuration).
   - N-P-K block skipped entirely (stub-error if `flCropNut=.true.`).
   - CO2 file block skipped (stub-error if `swco2=1`).
   - Irrigation scheduling init call skipped (stub-error if `schedule=1`).

Defense-in-depth runtime guards at entry:

```fortran
if (cfg%soybean%swsoybean == 1 .or. cfg%bulb%swbulb == 1 .or. &
    cfg%nutrient%flcropnut .or. cfg%co2%swco2 == 1 .or. &
    cfg%schedule%schedule == 1 .or. cfg%drought_stress%swdrought == 2 .or. &
    cfg%oxygen_stress%swoxygen == 2 .or. cfg%interception%swinter == 2 .or. &
    cfg%compensate%swcompensate /= 0 .or. cfg%harvest%swharv == 1 .or. &
    cfg%root%swrdc == 1) then
   call fatalerr_collected('cropwofost_init', &
      'Unsupported runtime branch reached on the TOML path. ' // &
      'The validator should have caught this earlier.')
end if
```

### Unit 5 — Module-level config reference (Phase 1 artefact, reused)

`src/crop/crop_config_global.f90` — created by Phase 1. Phase 2 adds a `use`
statement to `src/crop/cropgrowth.f90`'s `wofost` subroutine and registers
`src/crop/cropwofost_init.f90` in `meson.build`. No modifications to the
global module itself.

### Unit 6 — Wiring change in `cropgrowth.f90:984`

Before:
```fortran
call readwofost(icrop, cropfil(icrop), swhydrlift, swsoybean, mg, dvsi, &
                dvrmax1, dvrmax2, flrfphotoveg, tmaxdvr, tmindvr, toptdvr, &
                popt, pcrt, flphenodayl, FraDeceasedLvToSoil)
```

After:
```fortran
if (associated(crop_config_global) .and. &
    allocated(crop_config_global%rotation_wofost) .and. &
    crop_config_global%rotation_wofost(icrop)%populated) then
   call cropwofost_init_from_config( &
           crop_config_global%rotation_wofost(icrop), icrop)
else
   call readwofost(icrop, cropfil(icrop), swhydrlift, swsoybean, mg, dvsi, &  ! transitional
                   dvrmax1, dvrmax2, flrfphotoveg, tmaxdvr, tmindvr, toptdvr, &
                   popt, pcrt, flphenodayl, FraDeceasedLvToSoil)
end if
```

The `else` branch is **transitional** (see Teardown). The full 16-argument
call must be preserved verbatim in the `else` branch.

## Teardown plan (strangler-fig hygiene)

Per the project's strangler-fig discipline, every transitional element
introduced by this phase must have a documented teardown trigger.

| Element | Introduced for | Teardown trigger | What gets removed | Replacement |
|---|---|---|---|---|
| `else call readwofost(...)` fallback in `cropgrowth.f90:984` | Phase 2 ships before Phases 3/4; type-1/3 rotations still need the legacy reader | Phase 4 (hupselbrook) lands all three types | The `if/else` dispatch around all three legacy reader calls | Unconditional `*_init_from_config` calls |
| `readwofost` subroutine in `readswap.f90` (~916 lines) | Parity-test fixture (per ADR 0015) | Future cleanup rewrites parity tests to fixture values | The subroutine itself | Fixture-based parity tests using stored expected values |
| `populated :: logical` on `cropwofost_config_t` | Runtime dispatch sentinel while legacy fallback is still reachable | Same as `crop_config_global` (config-passing direction, post-Phase 4) | The `populated` field | Dispatch at call site on `rotation_type(icrop)`; no sentinel needed |
| `crop_config_global` module pointer (Phase 1 artefact) | Bridge from typed config to legacy globals | ADR 0016 Part C (config-passing direction) lands | The `crop_config_global_mod` module | Explicit config+state threading |
| `tests/swap-cases/toml/5.salinitystress/potatod.crp` (legacy ASCII in TOML case dir) | Legacy ASCII left in TOML dir during transition | Task 9 of Phase 2 plan | The file | None — the TOML executable reads only `potatod.crp.toml` |
| Soybean / bulb / nutrient / CO2 stub-errors in `cropwofost_config_validate` | Narrow scope for Phase 2 | A future case authors one of those branches | The stub-error for that branch | Full implementation of the branch (schema already 1:1) |

The `tests/swap-cases/5.salinitystress/potatod.crp` in the legacy case dir
is **not transitional** — it stays for the legacy executable and the parity
test fixture.

## Test plan

### New unit tests

- `tests/unit/config/test_cropwofost_config.pf` (extend) — add stub-error
  rejection tests for `swsoybean=1`, `swbulb=1`, `flcropnut=.true.`,
  `swco2=1`, `schedule=1`, `swdrought=2`, `swoxygen=2`, `swinter=2`,
  `swcompensate=1`, `swcompensate=2`, `swharv=1`, `swsalinity=2`. Add a
  "case 5 supported values pass" test with all active switches at their
  case-5 values.
- `tests/unit/io/toml/test_read_cropwofost_toml.pf` (extend) — round-trip
  the full `potatod.crp.toml` fixture; assert `populated=.true.` at end.
- `tests/unit/io/toml/test_load_swap_config.pf` (extend) — assert that
  loading case 5's `swap.toml` populates all 4 `rotation_wofost(1..4)`
  slots with `populated=.true.` and `rotation_loaded(1..4) = .true.`.
- `tests/unit/crop/test_cropwofost_init.pf` (new file) — assert that
  `cropwofost_init_from_config` writes correct values to `variables` module
  globals for case 5's potato data: `swrd`, `rdctb`, `cumdens` (hand-computed),
  `swsalinity`, `saltmax`, `saltslope`, `swoxygen`, `hlim1`, `hlim2u`,
  `hlim2l`, `swdrought`, `hlim3h`, `hlim3l`, `hlim4`, `adcrh`, `adcrl`,
  `swinter`, `cofab`, `swco2`, `relmf`.

### Parity test (extend existing)

`tests/unit/io/toml/test_salinitystress_parity.pf` already contains
`test_salinitystress_parity_potatod_wofost`, which compares `readwofost`
output globals to `rotation_wofost(1)` schema fields. Phase 2 adds a
companion assertion block (`test_salinitystress_wofost_init_parity`) that:
1. Calls `cropwofost_init_from_config(rotation_wofost(1), 1)`.
2. Asserts the same `variables`-module globals as those set by `readwofost`.
3. Adds assertions for `cumdens(1..202)` hand-computed from case 5's
   2-row `rdctb = [[0.0, 1.0], [1.0, 0.0]]` (linear ramp → triangular
   integration → known analytic result: `cumdens(i) = (i-1)^2 / 100^2` for
   the normalized form).

### Regression

Case 5 (salinitystress) must remain green. The smoke test (Task 10) renames
`potatod.crp` to `potatod.crp.disabled` in the TOML case dir and re-runs
regression to confirm 5/5, restoring afterward.

### Acceptance gate

- `pixi run test-pfunit` → all green.
- `pixi run regression` → 5/5 cases green.
- `tests/swap-cases/toml/5.salinitystress/potatod.crp` → does not exist.
- `git grep "call readwofost" src/` → returns the legacy fallback in
  `cropgrowth.f90:984` only (the `else` branch). No other reachability.
- `git grep "call readwofost" tests/` → returns parity-test references
  via `legacy_crop_helper`.
- `crop_config_global` populated by `config_to_variables` (from Phase 1),
  read by `cropgrowth.f90`'s `wofost` subroutine.
- `cropwofost_init_from_config` in `src/crop/cropwofost_init.f90` exists,
  registered in `meson.build`.

## Implementation tasks (working draft — see plan for full detail)

1. Audit `readwofost` (readswap.f90:2520-3435) line-by-line into
   READ/VALIDATE/NORMALIZE/RUNTIME/GUARDED buckets. Docs only.
2. Extend `cropwofost_config_t` with missing fields (soybean, bulb, swrdc,
   nutrient); add `populated` sentinel; add stub-error validators.
3. Extend `read_cropwofost_toml.f90` parser for new fields; set
   `populated = .true.` at end.
4. Complete `potatod.crp.toml` with `[irrigation_schedule]` section.
   Submodule pair commit.
5. Create `src/crop/cropwofost_init.f90` with `cropwofost_init_from_config`
   (config-to-globals copy + `cumdens` math + defense-in-depth guards).
6. Wire `cropgrowth.f90:984` to dispatch on the `populated` sentinel;
   legacy fallback in the `else` branch.
7. Smoke test: rename `potatod.crp` → `potatod.crp.disabled`, run regression
   5/5, restore.
8. Delete `potatod.crp` from the case-5 TOML dir. Submodule pair commit.

## Risk register

- **`.END`-file restart block gated on `swinco=3`.** Case 5 has `swinco=3`
  but the block only fires when `t1900 - tstart < 1e-3` AND `swinco=3` AND
  `|t1900 - cropstart(icrop)| >= 1e-3`. For case 5 the 2012 crop starts
  after the simulation start (April vs January) so the condition is `false`
  and the block is skipped. However `cropwofost_init_from_config` must
  replicate the conditional precisely, including for rotations 2/3/4 (2013,
  2014, 2015) which also start mid-year. Mitigation: the smoke test
  (`pixi run regression` with `potatod.crp.disabled`) validates all 4
  rotations produce identical output to the legacy path.
- **16-argument `readwofost` in the `else` fallback.** The `else` branch
  must preserve all 16 arguments verbatim — `swsoybean`, `mg`, `dvsi`,
  `dvrmax1`, `dvrmax2`, `flrfphotoveg`, `tmaxdvr`, `tmindvr`, `toptdvr`,
  `popt`, `pcrt`, `flphenodayl` are local variables in the `wofost` subroutine.
  The dispatch wiring must not disturb their declaration or their values on
  the legacy path. Mitigation: read `cropgrowth.f90:882-986` carefully before
  inserting the dispatch block.
- **`cumdens` computation depends on `swdrought`.** The computation only runs
  when `swdrought=1`. Case 5 has `swdrought=1`, so this is covered. The
  hand-computed assertion for the parity test uses case 5's rdctb = [[0.0,
  1.0], [1.0, 0.0]] (linear ramp): the integrated density over [0,1] is 0.5,
  so the normalized cumdens values at 101 points follow a quadratic shape.
  Mitigation: derive the expected values analytically and assert them in the
  unit test before writing the init sub.
- **`wofost_soybean_t` fields are Fortran `logical` and `real(8)` mixed.**
  The `flrfphotoveg` and `flphenodayl` legacy variables are `logical` in the
  `wofost` subroutine scope. The TOML schema should model them as `logical`
  fields. Mitigation: the schema-1:1 audit (Task 1) clarifies the exact types
  before writing the schema extension.
- **N-P-K nutrient block shares the same `.crp` file handle.** Legacy
  `readwofost` closes the `.crp` file before re-opening it for N-P-K via
  a second `rdinit` call. `cropwofost_init_from_config` never opens any
  file, so there is no file-handle concern. The nutrient block is simply
  skipped with a defense-in-depth guard. Risk is low; documented for clarity.
- **`potatod.crp.toml` is shared across 4 rotations.** All 4 rotation entries
  in case 5's `swap.toml` reference the same file. The loader reads it 4× —
  once per entry. Each slot in `rotation_wofost(1..4)` gets an independent
  populated copy. This is correct and consistent with Phase 1's pattern.
  Mitigation: the `test_load_swap_config.pf` extension asserts all 4 slots
  are populated.
