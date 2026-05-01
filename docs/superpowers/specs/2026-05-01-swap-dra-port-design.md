# swap.dra → TOML port — design

**Status:** Draft, awaiting review
**Date:** 2026-05-01
**Predecessors:** `docs/superpowers/specs/2026-05-01-swap-ini-port-design.md` (same pattern)

## Goal

Eliminate the legacy `swap.dra` ASCII file from the TOML pipeline. After this
work the new executable reads only `swap.toml`, `swap.dra.toml`, and CSV
companions. The legacy executable continues to read `swap.dra` from its own
case directory (`tests/swap-cases/<N>.<case>/`) and is unaffected.

The terminal state:
- `tests/swap-cases/toml/6.surfacewater/swap.dra` is deleted.
- `surfacewater.f90:56` no longer calls `rddre`.
- `rddre` is removed from `readswap.f90`, along with its upstream `drfil` /
  `pathdrain` reads.
- `tests/swap-cases/toml/6.surfacewater` runs end-to-end against the regression
  reference output via the new executable.

## Non-goals

- Supporting `.dra` branches that no TOML case exercises (`swsrf=3`,
  `swsec=1`, `swqhr=2`, `nrman2>0`). These are stubbed with explicit fatal
  errors so a future case authoring those switches gets a clear message
  rather than silent miscalculation. Full coverage lands when a test case
  demands it.
- Changing the legacy executable's behaviour. It continues to read
  `swap.dra` from `tests/swap-cases/<N>.<case>/`.
- Reworking the per-level drainage table for cases that don't have it
  (1.hupselbrook through 5.salinitystress already author what they need).

## Status quo

### Schema and parsers (already in place)

- `drainage_config_t` (Phase 4c-a) — basic drainage scalars, per-level
  arrays (`swdtyp`, `zbotdr`, `drares`, `infres`, `L`, `gwlinf`, `rdrain`,
  `rinfi`, `rentry`, `rexit`, `widthr`, `taludr`), `cofani`,
  `surface_runoff` sub-section.
- `surface_water_config_t` (Phase 4f-prep) — `swsrf`, `swsec`, `wlact`,
  `osswlm`, `nmper`, per-period arrays (`impend`, `swman`, `wscap`,
  `wldip`, `intwl`), `swqhr`, `sofcu`, weir arrays (`hbweir`, `alphaw`,
  `betaw`). Validator covers per-period sizes and enums; finalize already
  applies the `alphaw *= 8.64 * 100^(1-betaw) / sofcu` normalization.
- `read_drainage_toml.f90` parses `[drainage.levels]` array-of-tables —
  but currently only reads `swdtyp`, `swallo`, `zbotdr`, `drares`,
  `infres`, `L`, `owltab_file`. The surface-water-extended fields
  (`gwlinf`, `rdrain`, `rinfi`, `rentry`, `rexit`, `widthr`, `taludr`)
  are declared on the type but not yet parsed. **Gap to close in this
  port.**

### Adapter wiring (already in place)

`config_to_variables.f90` writes the same module globals that `rddre`
overwrites: `swsrf`, `osswlm`, `nmper`, `swqhr`, `impend`, `swman`,
`wscap`, `wldip`, `intwl`, `hbweir`, `alphaw`, `betaw`, plus
`drainage.surface_runoff` cluster. The adapter does not currently
populate `wls1`, `wlstar`, `sttab`, `swstini`, `swst`, `wlsbak`, or
`numadj` — those are computed inside `rddre`.

### TOML content (already in place)

Case 6's `swap.toml` already has:
- `[surface_water]` (swsrf=2, swsec=2, wlact, osswlm, nmper=28, swqhr=1, sofcu).
- `[surface_water.management]` (28-row `impend`/`swman`/`wscap`/`wldip`/`intwl`).
- `[surface_water.weir]` (28-row `hbweir`/`alphaw`/`betaw`).

Case 6's `swap.dra.toml` has the basic part 0 only (`swdivd`, `swdislay`,
`altcu`, `nrlevs=0`). The `[[drainage.levels]]` table is missing — the
file flags it as "Phase 4d skipped" in a header comment.

### The blocker: `rddre`

`readswap.f90:4273-4915` is a 642-line subroutine called from
`drainage/surfacewater.f90:56`. It opens `swap.dra`, parses 62
`rd*` calls, and runs initialization math afterwards. Even on the
TOML path it executes at simulation init and clobbers the
adapter-populated globals from the file on disk. That's why
`swap.dra` cannot simply be deleted — `rddre` would fail to open it.

The 642 lines decompose into three buckets:

1. **File reads + range checks** — 62 `rd*` calls. Goes away with the port.
2. **Cross-section validation** — weir crest vs deepest channel bottom,
   target-level vs supply consistency, period-bookkeeping (nmper sums),
   q-h table consistency (4d), automatic-weir cross-checks (4e1/4e2).
3. **Coordinate normalization** — `hbweir[i] -= altcu`, `alphaw` formula
   (already in finalize), `wls1 = wlact - altcu`, `wldip[i] = abs(...)`,
   `hdepth[i] = -abs(...)` (4e1).
4. **Runtime numerical init** — `sttab` storage table (`widthr`/`taludr`/
   `zbotdr` math), `swstini = swstlev(wls1)`, `wlsbak(1:4) = 0`,
   `numadj = 0`, initial `wls1`/`wlp1` interpolation from `wlstab`/
   `wlptab` (only swsec=1 / swsrf=3 path).

The case-6 subset exercises only buckets 2/3/4 corresponding to swsrf=2 +
swsec=2 + swqhr=1 + all swman=1. Buckets for swsrf=3 / swsec=1 / swqhr=2 /
nrman2>0 don't run for case 6.

## Approach

Three units, each with one purpose, replacing `rddre`:

### Unit 1 — extend cross-section validation

Push (2) into config validators.

- `surface_water_config_validate` (already exists): tighten existing
  per-period checks; add `wldip[i] >= 0` if not already covered.
- `swap_config_validate`: add cross-section rules (need both drainage
  and surface_water in scope):
  - `hbweir[i] > zbotdr(1 + nrpri)` (weir crest above deepest channel
    bottom of secondary system)
  - `hbweir[i] - wldip[i] > zbotdr(1 + nrpri) + 1e-4` when
    `swman[i] == 1 .and. wscap[i] > 1e-7` (target level above bottom
    when supply is attempted)
  - `nrlevs >= 1 + nrpri` when `swsrf == 2`

For now `nrpri = 0` (swsrf=2 has no primary system); when swsrf=3 lands,
`nrpri = 1` and the same rules generalize.

The unimplemented branches (swsrf=3, swsec=1, swqhr=2, nrman2>0) get a
single guard at the top of the cross-section block:

```fortran
if (self%surface_water%swsrf == 3 .or. self%surface_water%swsec == 1 .or. &
    self%surface_water%swqhr == 2 .or. count(self%surface_water%swman == 2) > 0) then
   call errors%append(ERR_VALIDATION_NOT_SUPPORTED, &
      'surface_water: swsrf=3, swsec=1, swqhr=2, and swman=2 are not yet ' // &
      'supported in the TOML pipeline. Use the legacy executable for ' // &
      'cases requiring these branches.', 'swap_config')
end if
```

### Unit 2 — extend finalizers for coordinate normalization

Push (3) into config finalizers. `surface_water_config_finalize` already
has the alphaw precedent.

- `hbweir[i] -= altcu` (per-period). `altcu` is in `drainage_config_t`,
  not `surface_water_config_t` — so this becomes a finalize step in
  `swap_config_finalize` (cross-section), wrapping
  `surface_water%finalize` with an extra pass that consumes
  `drain%altcu`.
- `wldip[i] = abs(wldip[i])` — pure within-config; goes in
  `surface_water_config_finalize`.

After Unit 2, the typed config carries values in the same coordinate
system that `rddre` produced. The adapter is unchanged.

### Unit 3 — extract runtime init into a new module

Push (4) into a new module. Bucket (1) is deleted.

New file: `src/drainage/surfacewater_init.f90`
New public subroutine: `surfacewater_init(wls1, wlp1)`

Operates on already-populated module globals (filled by
`config_to_variables`). Does:

- `numadj = 0`
- `wlsbak(1:4) = 0.0`
- `wls1 = wlact - altcu` (using the wlact written by adapter; wlact
  itself stays in the typed config — see "Decision log" below)
- `wlstar = wls1`
- `sttab` build (lines 4862-4897 verbatim — pure widthr/taludr/zbotdr math)
- `swstini = swstlev(wls1)`; `swst = swstini`
- For unimplemented branches: a single fatal-error gate matching the
  validator (defense in depth):

```fortran
if (swsrf == 3 .or. swsec == 1 .or. swqhr == 2) then
   call fatalerr_collected('surfacewater_init', &
      'swsrf=3, swsec=1, swqhr=2 not yet supported on the TOML path')
end if
```

`drainage/surfacewater.f90:56` is updated:

```fortran
! before:
call rddre (wls,wlp)

! after:
call surfacewater_init (wls, wlp)
```

### Unit 4 — TOML authoring + parser extension

- Extend `read_drainage_toml.f90`'s `[[drainage.levels]]` loop to also
  parse `gwlinf`, `rdrain`, `rinfi`, `rentry`, `rexit`, `widthr`,
  `taludr` (currently declared on the type but not parsed).
- Author case 6's `[[drainage.levels]]` table from the legacy `swap.dra`
  Part 1 LEV/SWDTYP/L/ZBOTDRE/GWLINF/RDRAIN/RINFI/RENTRY/REXIT/WIDTHR/
  TALUDR table (2 rows).
- Author `cofani = [10.0, 10.0, 10.0]` at `[drainage]` top level
  (matches the 3 soil-physical layers in case 6).
- Bump `nrlevs` from 0 to 2 in case 6's `swap.dra.toml`.
- Strip the "Phase 4d skipped" comment block.

### Unit 5 — delete `rddre` and upstream

Once all callers route through `surfacewater_init`:

- Delete subroutine `rddre` (`readswap.f90:4273-4915`).
- Delete `subroutine checkdate` if only `rddre` calls it. (Verify.)
- Delete `if (swdra .ne. 0) call rdscha ('drfil',drfil)` at
  `readswap.f90:1003` — `drfil` is now dead.
- Delete `drfil` argument propagation at `readswap.f90:1766,1783` and
  the `drfil` field from any common/module that only `rddre` reads.
- Delete `pathdrain` if no other caller. (Verify.)

After this `readswap.f90` shrinks by ~700 lines and no code anywhere
opens `swap.dra`.

### Unit 6 — delete `swap.dra` from case 6

`tests/swap-cases/toml/6.surfacewater/swap.dra` is removed. Submodule
gets a paired commit (inner repo + outer-repo bump). The legacy ASCII
file remains in `tests/swap-cases/6.surfacewater/swap.dra` for the
legacy executable.

`tests/swap-cases/run_case.sh` and the regression harness need no
changes — neither references `swap.dra` directly when running the
TOML executable.

### Unit 7 — docs

- `docs/csv-companion-files.md`: update the "Path resolution and
  staging" section to remove the `swap.dra` (`swdra=2`) parenthetical.
- `docs/configuration-schema.md`: extend `[[drainage.levels]]` row
  examples to cover the seven new fields if not already covered.

## Decision log

**Decided: narrow port scope to the case-6 subset.** The unimplemented
branches (swsrf=3, swsec=1, swqhr=2, nrman2>0) are gated with explicit
fatal errors at validation and runtime. No silent misbehaviour. Future
work expands coverage when a test case authors those switches.

**Decided: hbweir altcu-subtraction lives in `swap_config_finalize`,
not `surface_water_config_finalize`.** It's a cross-section
normalization (needs `drain%altcu`). Mirrors the existing `swap_config`
cross-section validation pattern.

**Decided: `wls1 = wlact - altcu` lives in `surfacewater_init`, not
finalize.** `wls1` is runtime state, not config. The typed config
holds `wlact` in its authored form (relative to reference level).
The runtime sub does the offset translation.

**Decided: schema-level fields for unimplemented branches not added
in this port.** No `[surface_water.qh_table]` for swqhr=2, no
`[surface_water.automatic_weir]` for swman=2, no
`[surface_water.primary_levels]` for swsrf=3. Adding them speculatively
when no test exercises them creates dead schema. Add when a case lands.

## Test plan

### New unit tests

`tests/unit/config/test_surface_water_config.pf` (extend or create):
- `test_surface_water_swsrf3_rejected` — error code matches stub.
- `test_surface_water_swsec1_rejected` — error code matches stub.
- `test_surface_water_swqhr2_rejected` — error code matches stub.
- `test_surface_water_swman2_rejected` — error code matches stub.
- `test_surface_water_finalize_wldip_abs` — negative wldip becomes
  positive after finalize.

`tests/unit/config/test_swap_config.pf` (extend):
- `test_swap_config_hbweir_below_channel_bottom_fails` —
  hbweir + altcu < zbotdr(1+nrpri) trips cross-section error.
- `test_swap_config_target_below_channel_with_supply_fails` —
  hbweir-wldip < zbotdr+1e-4 with swman=1, wscap>0 trips error.
- `test_swap_config_finalize_hbweir_subtracts_altcu` — after finalize,
  hbweir reflects `authored - altcu`.

`tests/unit/io/toml/test_read_drainage_toml.pf` (extend):
- `test_drainage_levels_parses_gwlinf_rdrain_etc` — round-trip the 7
  new per-level fields.

`tests/unit/drainage/test_surfacewater_init.pf` (new):
- `test_sttab_open_channel_two_levels` — direct math test of the
  storage-vs-level table for a known (widthr, taludr, zbotdr) input.
- `test_swstini_matches_swstlev` — initial storage equals `swstlev(wls1)`.

### Regression coverage

`tests/swap-cases/toml/6.surfacewater` runs against
`tests/regression/test_output_regression.py` and produces output
identical to the legacy reference. No change in tolerance, no
change in output schema.

### Acceptance

- `pixi run unit-tests` green.
- `pixi run regression-tests` 5/5 cases green (no perf regression
  beyond noise).
- `git grep rddre src/` returns nothing.
- `git grep swap.dra src/` returns only the legacy `tests/swap-cases/`
  reference.
- `tests/swap-cases/toml/6.surfacewater/swap.dra` does not exist.

## Risk register

- **`sttab` math is subtle**: lines 4862-4897 mix integer-divide trapezium
  geometry with rectangle geometry above soil surface. Direct verbatim
  port (no rewrite) protects parity. The new unit test for sttab uses
  the exact case-6 inputs and reference outputs sampled from a current
  legacy run.
- **Hidden `rddre` invariants**: the audit above is from a single read.
  Implementation Task 1 should re-audit line-by-line and classify each
  line into a bucket before cutting. Document the audit in a one-page
  table inside the implementation plan.
- **`drfil`/`pathdrain` cleanup may cascade**: these globals may be
  read by other paths. Implementation Task 5 verifies with `git grep`
  before deletion; if other readers exist, scope grows or cleanup is
  deferred to a follow-up.
- **Legacy executable build**: removing `rddre` from `readswap.f90`
  may break the legacy executable build target if the TOML pipeline
  shares object files with it. Verify build matrix early.
