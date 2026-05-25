# Port `swap.ini` to TOML + CSV Companions — Design

**Date:** 2026-05-01
**Phase:** 4f cleanup, follow-up to CSV meteo finalization
**Status:** approved

## Goal

Replace the last legacy ASCII initial-conditions file (`swap.ini`) consumed by
the modern TOML pipeline with a structured `[soil.initial]` TOML sub-section
for scalars and per-profile CSV companion files for the z-indexed arrays.
After this lands, the salinitystress case no longer needs a `swap.ini` file
in its TOML directory; the only remaining legacy ASCII companions in the
TOML pathway are the `.crp` crop files and the surfacewater extended-drainage
`swap.dra`.

## Context

The CSV meteo finalization tagged at `rescue/phase-csv-meteo-complete` left
two known stop-gap legacy ASCII files in the TOML case directories:

- `5.salinitystress/swap.ini` — initial conditions for `swinco=3` mode
- `6.surfacewater/swap.dra` — extended drainage for `swdra=2` mode

Of those, `swap.ini` is the smaller and simpler port: it's read in one
contiguous block in the adapter (`src/io/toml/config_to_variables.f90:619–650`)
via TTutil `rdinit`/`rdsdor`/`rdador`/`rdfdor`. The block reads three scalars
(`ssnow`, `slw`, `pond`), one mandatory profile (`z_h`/`h`), and two
optional profiles (`z_Tsoil`/`Tsoil` when heat is active; `z_Cml`/`Cml` when
solute is active).

Inspection of the actual `salinitystress/swap.ini` reveals four additional
scalar fields the current adapter does not consume: `SWIRRIGATE`, `ldwet`,
`dt`, and the 7-element `atmin7` array. These are warm-restart artifacts from
`Result.end` (the file SWAP writes at simulation end). The user has confirmed
these should also move to the typed schema even though the current adapter
ignores them — keeping the TOML schema as the canonical home for all
initial-state inputs.

The `swap.dra` port is deferred to a separate spec (it touches ~30 fields
across 4 sections of the legacy file and pulls in surface-water management
schema).

## Architecture

**Replace the entire `swinco==3 .and. allocated(inifil)` block** in
`config_to_variables.f90` with two read paths:

1. **Scalars and small fixed arrays** → typed slots in a new
   `[soil.initial]` TOML sub-section. The adapter copies them directly to
   the legacy globals (`ssnow`, `slw`, `pond`, `pondini`, `ldwet`, `dt`,
   `atmin7(:)`). `swirrigate` is metadata only — the legacy reader uses it
   locally to control parse flow in `readswap.f90`, but the TOML pipeline
   does not have an equivalent global to update; the slot exists in the
   schema for completeness so the value round-trips, but the adapter
   ignores it.
2. **Z-indexed profile tables** → 1–3 separate CSV companion files read via
   the existing `read_csv_table` reader (real-keyed col 1).

Per-profile CSVs (rather than one wide CSV) follow the legacy semantics:
each profile is independently sized and may have its own z-grid. Each CSV
has a minimal `z,<value>` two-column header, matching the convention
established by other CSV companions (`<case>.gwl.csv`, `<case>.haquif.csv`).

The legacy `[soil].inifil` slot is **removed**. The new `[soil.initial]`
sub-section is the canonical home; there is no fallback or compatibility
shim.

The 4 cases that use `swinco=2` (or `swinco=0`) are unaffected: the new
sub-section's defaults make the entire block a no-op.

## Components

| Area | File | Change |
|---|---|---|
| Schema — type | `src/config/soil_config.f90` | Add nested `type :: soil_initial_t` with 7 scalar fields, 1 fixed-7 real array, 3 file-path strings. Add `type(soil_initial_t) :: initial` member to `soil_config_t`. Remove `inifil` slot. |
| Schema — validator | `src/config/soil_config.f90` | Add `soil_initial_t%validate(parent_swinco, swhea, swcalt, swsolu, errors)` that fires `ERR_VALIDATION_REQUIRED` per gating combination (see Validator below). Wire from `swap_config_t%validate` (which has visibility into all sections). |
| Schema test | `tests/unit/config/test_soil_config.pf` | Add tests for the new fields + the 3 validator rules (each with one positive and one negative case). |
| TOML reader | `src/io/toml/read_soil_toml.f90` | Read the `[soil.initial]` sub-table: 7 scalar fields, 1 array (TOML inline array of length 7), 3 string slots. Drop the `inifil` read. |
| TOML reader test | `tests/unit/io/toml/test_read_soil_toml.pf` | Add a test that loads a `[soil.initial]` fixture and asserts all 11 fields populate. |
| Adapter | `src/io/toml/config_to_variables.f90` | Replace lines 619–650 with: scalar copy + 3 conditional `read_csv_table` calls. Use existing `block` scoping pattern. |
| Cross-section validation | `src/config/swap_config.f90` | Update `swap_config_t%validate` to call `soil_initial%validate(swinco, heat.swhea, heat.swcalt, solute.swsolu, errors)` with the cross-section context. |
| Case data — schema | `tests/swap-cases/toml/5.salinitystress/swap.toml` | Remove `[soil].inifil = "swap.ini"`. Add `[soil.initial]` block populated from `swap.ini` scalars (`swirrigate=1`, `ssnow=0`, `slw=0`, `pond=0`, `ldwet=10.0`, `dt=1.0e-7`, `atmin7=[14.1, 13.2, 13.8, 16.0, 18.4, 16.3, 16.6]`) and three `*_file` slots. |
| Case data — h profile | `tests/swap-cases/toml/5.salinitystress/salinitystress.ini.h.csv` | New file. Header `z,h`. ~180 rows extracted from `swap.ini` `z_h,h` table. |
| Case data — Tsoil profile | `tests/swap-cases/toml/5.salinitystress/salinitystress.ini.tsoil.csv` | New file. Header `z,tsoil`. ~180 rows. |
| Case data — Cml profile | `tests/swap-cases/toml/5.salinitystress/salinitystress.ini.cml.csv` | New file. Header `z,cml`. ~180 rows. |
| Case data — legacy delete | `tests/swap-cases/toml/5.salinitystress/swap.ini` | **Delete from submodule.** No longer required at runtime. |
| Docs | `docs/csv-companion-files.md` | Document the 3 new slots in the schema table; mention `swap.ini` removal in migration history. |
| Docs | `docs/configuration-schema.md` | Document the `[soil.initial]` schema if a `[soil]` section exists in this doc. |

## TOML schema example

```toml
[soil]
swinco = 3                     # (no longer takes inifil; reads [soil.initial])

[soil.initial]
swirrigate = 1
ssnow      = 0.0
slw        = 0.0
pond       = 0.0
ldwet      = 10.0
dt         = 1.0e-7
atmin7     = [14.1, 13.2, 13.8, 16.0, 18.4, 16.3, 16.6]
h_file     = "salinitystress.ini.h.csv"
tsoil_file = "salinitystress.ini.tsoil.csv"
cml_file   = "salinitystress.ini.cml.csv"
```

## CSV file shapes

```
salinitystress.ini.h.csv     →   z,h         (header) + N rows
salinitystress.ini.tsoil.csv →   z,tsoil     (header) + M rows
salinitystress.ini.cml.csv   →   z,cml       (header) + K rows
```

`N`, `M`, `K` are independent. Each CSV is read separately via
`read_csv_table` and produces an independent `(rows, 2)` table that the
adapter unpacks into the corresponding paired arrays.

## Adapter logic (replaces lines 619–650)

```
if (swinco == 3) then
   ssnow      = config%soil%initial%ssnow
   slw        = config%soil%initial%slw
   pond       = config%soil%initial%pond
   pondini    = pond
   ldwet      = config%soil%initial%ldwet
   dt         = config%soil%initial%dt
   atmin7(:)  = config%soil%initial%atmin7(:)
   ! config%soil%initial%swirrigate is read but unused — no TOML-side consumer
   if (swsnow /= 1) ssnow = 0.0d0    ! preserve legacy zeroing rule

   ! Mandatory: initial pressure-head profile
   call read_csv_table(h_file, ['z', 'h'], tbl, errs)
   call errs%abort_if_fatal()
   nhead = size(tbl, 1)
   zi(1:nhead) = tbl(:, 1)
   h(1:nhead)  = tbl(:, 2)

   if (swhea == 1 .and. swcalt == 2) then
      call read_csv_table(tsoil_file, ['z', 'tsoil'], tbl_t, errs_t)
      call errs_t%abort_if_fatal()
      n_t = size(tbl_t, 1)
      zh(1:n_t)    = tbl_t(:, 1)
      tsoil(1:n_t) = tbl_t(:, 2)
   end if

   if (swsolu == 1) then
      call read_csv_table(cml_file, ['z', 'cml'], tbl_c, errs_c)
      call errs_c%abort_if_fatal()
      nconc = size(tbl_c, 1)
      zc(1:nconc) = tbl_c(:, 1)
      cml(1:nconc) = tbl_c(:, 2)
   end if
end if
```

(The pseudocode above shows shape; the actual implementation will use the
project's existing `block ... end block` scoping and error-collection
patterns from the daily/rain/detail meteo CSV pre-loads.)

## Validator

`soil_initial_validate(self, parent_swinco, swhea, swcalt, swsolu, errors)`:

| Condition | Error code | Message |
|---|---|---|
| `swinco=3` and `h_file` empty | `ERR_VALIDATION_REQUIRED` | `soil.initial.h_file required when soil.swinco=3` |
| `swinco=3 ∧ swhea=1 ∧ swcalt=2` and `tsoil_file` empty | `ERR_VALIDATION_REQUIRED` | `soil.initial.tsoil_file required when soil.swinco=3 and heat.swhea=1 and heat.swcalt=2` |
| `swinco=3 ∧ swsolu=1` and `cml_file` empty | `ERR_VALIDATION_REQUIRED` | `soil.initial.cml_file required when soil.swinco=3 and solute.swsolu=1` |

The cross-section context (`swhea`, `swcalt`, `swsolu`) is passed by
`swap_config_t%validate`. The `soil_initial_t` itself does not depend on
other config types.

Range validation for the scalar fields:

| Field | Range |
|---|---|
| `swirrigate` | `[0, 1]` |
| `ssnow` | `[0.0, 1000.0]` |
| `slw` | `[0.0, 1000.0]` |
| `pond` | `[0.0, 100.0]` |
| `ldwet` | `[0.0, 366.0]` |
| `dt` | `[1.0e-12, 1.0]` |
| `atmin7(i)` | `[-50.0, 50.0]` |

Profile CSV content (z monotonic, h within depth range, etc.) is not
validated by the schema validator — that is the reader's job and is
handled at adapter time when `read_csv_table` returns malformed data.

## Data flow

```
swap.toml ──load_swap_config──▶ swap_config_t (with [soil.initial])
                                        │
                                        ▼  swap_config%validate
                              soil.initial%validate(swinco, swhea, swcalt, swsolu)
                              required-file checks fire here
                                        │
                                        ▼  config_to_variables (adapter)
                              if swinco==3:
                                  copy 7 scalars + atmin7(:) → globals
                                  read_csv_table(h_file, ['z','h'])     → zi, h, nhead
                                  read_csv_table(tsoil_file, ...)       → zh, tsoil
                                  read_csv_table(cml_file, ...)         → zc, cml, nconc
```

## Error handling

Validator-time:

| Condition | Error code |
|---|---|
| Required `*_file` slot empty | `ERR_VALIDATION_REQUIRED` |
| Scalar out of range | `ERR_VALIDATION_OUT_OF_RANGE` (existing) |

Reader-time (already covered by `read_csv_table`):

| Condition | Error code |
|---|---|
| File not found | `ERR_IO_OPEN_FAILED` |
| Header mismatch (e.g. `z,Cml` vs expected `z,cml`) | `ERR_PARSE_HEADER_MISMATCH` |
| Numeric parse failure | `ERR_PARSE_TYPE_MISMATCH` |

All errors flow through `error_collection_t`; the run aborts at the first
checkpoint after any error.

## Testing

**Schema validator unit tests** (`tests/unit/config/test_soil_config.pf`):

- One positive case per validator rule: all required slots present,
  validator returns no errors.
- One negative case per validator rule: required slot empty, validator
  emits the right error code with the right field name in the message.
- Range tests for the 7 scalars (one out-of-range case each).

**TOML reader test** (`tests/unit/io/toml/test_read_soil_toml.pf`):

- Load a fixture with the full `[soil.initial]` block populated and
  assert all 11 fields round-trip correctly.

**Regression** (`pixi run -e test check-full`):

- `salinitystress`: must remain green at 1e-2 cm tolerance after the
  migration. The change is purely about *where* the data lives, not the
  numeric content.
- Other 4 cases: must stay green; the new section is dormant for them.

**Manual smoke test:**

```bash
cd tests/swap-cases
./run_case.sh -c salinitystress -e ../../builddir/swap
ls toml/5.salinitystress/   # confirm no swap.ini lingering
```

## Out of scope

- **Surfacewater `swap.dra` port** (extended drainage) — separate spec.
- **Crop `.crp` port** — separate spec, larger scope.
- **Other 4 regression cases** — they use `swinco=0` or `swinco=2`,
  unaffected by this change.
- **Per-array z-grid validation** — z monotonic / depth bounds are
  reader-level concerns, not schema-level.
- **Compatibility shim for legacy `inifil`** — none. The slot is removed
  outright.
- **Removing TTutil `rdsdor`/`rdador`/`rdfdor` from elsewhere** — the
  `.swp` pathway still calls these in other contexts.

## Risks

| Risk | Mitigation |
|---|---|
| Numeric drift in salinitystress regression | The CSV row values must be byte-equal in floating-point representation to what the legacy `rdfdor` parsed from `swap.ini`. Use Fortran's default real-formatter (matches legacy `*` format-spec). Run regression after the migration; any drift is a fixture-extraction bug, not an architectural one. |
| `[soil.initial]` array length for `atmin7` mismatch | TOML inline array. Reader checks `size == 7` and emits a clear error if not. |
| Cross-section validation coupling | `soil_initial%validate` takes the dependent switches as args, not as references to other config types. Keeps the validator decoupled. The composition happens in `swap_config_t%validate`. |
| `dt` default of `0.0` could cause solver divide-by-zero on `swinco=2` cases | The adapter only consumes `dt` when `swinco=3`. For other values of `swinco`, the existing solver-init logic sets `dt`; the new field is ignored. Confirm by inspection of `config_to_variables.f90` and `swap.f90` solver-init. |
| Removing `inifil` breaks other (non-regression) test fixtures | `git grep -n inifil` to confirm scope. Update or delete any stale references. |
