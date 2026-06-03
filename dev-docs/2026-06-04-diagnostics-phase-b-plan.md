# Diagnostics Phase B — TOML block + stray-diagnostic cleanup — Plan

**Goal:** Put a `[logging]` TOML block into the precedence chain (default < TOML < env <
C-API), route the genuine stray diagnostics through the logger, and drop dead debug —
byte-identical, no `swap_state_t` schema change.

**Base:** branch `diagnostics-bcde` off `development` (75ed6ea). Phase A is already in.

**Conventions:** `pixi run -e test check-fast` green before each commit (byte-identical
4/4). New `.pf` registered in `pf_files` + `pfunit_extra_sources` + `testSuites.inc`.
pFUnit gotcha: no inline comments after `@assert`; one assertion per line. Targeted
`git add`. Commit trailer `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`.

Scope notes (decided during planning):
- The `timecontrol_mod.f90` screen writes are gated by `swscre==2` (a deliberate
  "print day numbers to screen" feature) and the per-day one is per-timestep — routing
  them to INFO would violate the per-run-INFO rule. **Leave them.**
- `src/soilwater/dormant/sptabulated.f90` magic-unit writes are dormant (unreached).
  Low value; **leave** (Phase B does not touch dormant code).

---

## Task B1 — `[logging]` TOML overrides in the precedence chain

**Files:** `src/error/diagnostics.f90` (extractor + `init_logging_from_toml`), the 4
entry points, test `tests/unit/error/test_diagnostics.pf`.

Add to `diagnostics_mod`:
- An extractor that pulls a `[logging]` table from a parsed toml doc into a
  `diag_overrides_t` (absent table → empty overrides, no error):

```fortran
   !> Extract [logging] overrides from a parsed TOML document.
   !! Recognised keys: level (string), file (string), to_stdout, to_stderr,
   !! timestamps (booleans). Absent [logging] table → empty overrides.
   subroutine extract_logging_overrides(doc, ov)
      use tomlf, only: toml_table, get_value
      type(toml_table), pointer, intent(in)  :: doc
      type(diag_overrides_t),    intent(out) :: ov
      type(toml_table), pointer     :: sec
      character(len=:), allocatable :: sval
      logical :: bval
      integer :: stat
      call get_value(doc, 'logging', sec, requested=.false.)
      if (.not. associated(sec)) return
      call get_value(sec, 'level', sval, stat=stat)
      if (stat == 0 .and. allocated(sval)) then
         if (len_trim(sval) > 0) then
            ov%has_level = .true.; ov%level = log_level_from_name(sval)
         end if
      end if
      call get_value(sec, 'file', sval, stat=stat)
      if (stat == 0 .and. allocated(sval)) then
         if (len_trim(sval) > 0) then; ov%has_file = .true.; ov%log_file = sval; end if
      end if
      call get_value(sec, 'to_stdout', bval, stat=stat)
      if (stat == 0) then; ov%has_stdout = .true.; ov%to_stdout = bval; end if
      call get_value(sec, 'to_stderr', bval, stat=stat)
      if (stat == 0) then; ov%has_stderr = .true.; ov%to_stderr = bval; end if
      call get_value(sec, 'timestamps', bval, stat=stat)
      if (stat == 0) then; ov%has_timestamps = .true.; ov%timestamps = bval; end if
   end subroutine extract_logging_overrides
```
  (Verify the exact `tomlf` `get_value` optional-key signatures against an existing
  reader, e.g. `src/io/toml/read_output_csv_toml.f90` / `read_tz_z1_z2`; adjust if the
  allocatable-string or logical overloads differ. If a key is absent, `stat /= 0`.)

- Two thin parse helpers + a TOML-aware init, all public:

```fortran
   !> Parse [logging] overrides from a TOML file path (missing file → empty).
   subroutine read_logging_overrides_from_file(path, ov)
      use tomlf, only: toml_table, toml_error, toml_load
      character(len=*),       intent(in)  :: path
      type(diag_overrides_t), intent(out) :: ov
      type(toml_table), allocatable, target :: doc
      type(toml_table), pointer             :: doc_ptr
      type(toml_error), allocatable         :: terr
      call toml_load(doc, path, error=terr)
      if (allocated(terr)) return            ! malformed/missing → no overrides
      doc_ptr => doc
      call extract_logging_overrides(doc_ptr, ov)
   end subroutine read_logging_overrides_from_file

   !> Parse [logging] overrides from a TOML text buffer (in-memory C-API path).
   subroutine read_logging_overrides_from_text(text, ov)
      use tomlf, only: toml_table, toml_error, toml_loads
      character(len=*),       intent(in)  :: text
      type(diag_overrides_t), intent(out) :: ov
      type(toml_table), allocatable, target :: doc
      type(toml_table), pointer             :: doc_ptr
      type(toml_error), allocatable         :: terr
      call toml_loads(doc, text, error=terr)
      if (allocated(terr)) return
      doc_ptr => doc
      call extract_logging_overrides(doc_ptr, ov)
   end subroutine read_logging_overrides_from_text

   !> Resolve base < toml < env and initialise the logger (one resolve, early).
   subroutine init_logging(base, toml_ov)
      type(diagnostics_config_t), intent(in) :: base
      type(diag_overrides_t),     intent(in) :: toml_ov
      type(diagnostics_config_t) :: dcfg
      type(diag_overrides_t)     :: none_ov, env_ov
      call read_env_overrides(env_ov)
      dcfg = resolve_diagnostics_config(base, toml_ov, env_ov, none_ov)
      call apply_logging_config(dcfg)
   end subroutine init_logging
```
  Refactor the existing `init_logging_from_env(base)` to call `init_logging(base,
  none_ov)`; factor the if/else `log_init` block into a private
  `apply_logging_config(dcfg)` used by both. Add the new public names.

**Entry-point wiring** (each reads its TOML overrides and calls `init_logging`):
- `swap_main.f90`: `call read_logging_overrides_from_file('swap.toml', toml_ov); call
  init_logging(default_cli_config(), toml_ov)`.
- `swap_bmi_mod.f90` `bmi_initialize`: `read_logging_overrides_from_file(f_config_file,
  toml_ov)` (after `c_to_f_string`), then `init_logging(default_embedded_config(),
  toml_ov)`.
- `swap_ensemble_mod.f90` `ensemble_init`: `read_logging_overrides_from_file(config_file,
  toml_ov)` then `init_logging(default_embedded_config(), toml_ov)`.
- `swap_capi_mod.f90` `swap_initialize_from_toml_string`: it already has the toml text
  (`f_text`); `read_logging_overrides_from_text(f_text, toml_ov)` then
  `init_logging(default_embedded_config(), toml_ov)` (replace the `capi_init_logging`/
  `init_logging_from_env` call). Make sure the text is available before the call (read
  it from the `buf` first if needed).

**Test** (TDD): add a `@test` that builds a small TOML text with a `[logging]` table and
asserts `extract_logging_overrides` (via `read_logging_overrides_from_text`) sets
`has_level`/level (e.g. `level="debug"`→LOGLEVEL_DEBUG), `has_stdout`/to_stdout=.true.,
and that a text WITHOUT `[logging]` yields all `has_*=.false.`. Build, `check-fast`,
commit.

Precedence test (optional, cheap): a `[logging] level="warn"` text + `SWAP_LOG_LEVEL`
unset → resolved level WARN; documents TOML-over-default.

---

## Task B2 — route stray solver/guard diagnostics through the logger

**Files:** `src/crop/oxygenstress.f90`, `src/atmosphere/interception.f90`,
`src/crop/cropgrowth_helpers.f90`.

These are byte-identical (console/dead output only; never the compared output files).
Run `check-fast` after — must stay 4/4.

1. **`src/crop/oxygenstress.f90`** — four bare `write (*,*) '<solver failure>'`:
   - line ~1379 `'QROMBD Too many steps.'`
   - line ~1418 `'QROMBDtab Too many steps.'`
   - line ~1458 `'NR_POLINTD: DEN = 0'`
   - line ~1566 `'ZBREND exceeding maximum iterations.'`
   Replace each `write (*,*) '<msg>'` with `call log_warn('oxygenstress', '<msg>')`.
   Add `use swap_log, only: log_warn` to the module (or the relevant subroutines) if not
   present. Keep surrounding control flow (the routines still return/stop as before —
   only the message emission changes; if a line is immediately followed by `return`/loop
   exit, leave that).

2. **`src/atmosphere/interception.f90`** lines ~257-259 (in the quarantined `msw1eic`):
   ```fortran
   write(ib,9199) k
   write(*,9199) k
   stop
   ```
   Replace all three with a single fatal through the error system:
   ```fortran
   write(messag,9199) k
   call fatalerr_collected('msw1eic', trim(messag))
   ```
   where `messag` is a `character(len=200)` local (declare if needed) and `9199` is the
   existing format. Add `use error_mod, only: fatalerr_collected`. This removes the last
   `log_unit_handle()` consumer — also delete the now-unused `ib_i4 = log_unit_handle()`
   line (~202) and its `use swap_log, only: log_unit_handle` import IF nothing else in
   the file uses `ib`/`log_unit_handle` (grep first; if `ib` is used elsewhere, leave the
   handle but still convert the stop). Note: `msw1eic` is dormant (swinter=3 path), so
   this is unreached in regression — byte-identical trivially, but keep the format/message
   text identical.

3. **`src/crop/cropgrowth_helpers.f90`** lines ~417-425 — the `! TEMPORARY OUTPUT DELETE`
   block writing to units 777/888. **Delete** the marked block entirely (the two
   `write(777,...)`/`write(888,...)` and their `TEMPORARY OUTPUT DELETE` comment lines and
   any `if`/guard that exists solely to gate them). Confirm units 777/888 are never
   `open`ed anywhere (`grep -n "777\|888" src` — they are not), so deletion changes
   nothing at runtime. Verify the surrounding loop/logic still compiles and the
   subroutine's real work is untouched.

Commit B2 as one or two commits (oxygenstress+interception together; cropgrowth-helpers
delete separate is fine). `check-fast` green each commit.

---

## Task B3 — drop unused log imports / context tidy (small)

**Files:** wherever Phase A/B left an unused `use swap_log, only: log_*` import, plus the
two known dead imports from the review (`src/config/simulation_config.f90`,
`src/soilwater/boundtop.f90` — confirm they still `use` a log routine they never call;
if so, remove the unused name from the `only:` list). Normalise the obviously-divergent
context string `'CropGrowth_Wofost'` → `'wofost'` ONLY in files you already touched (do
not sweep). Build + `check-fast`. Commit.

Keep this task minimal — if a removal causes an "unused" cascade or any ambiguity, skip
that item and note it rather than chasing it.

---

## Phase B done-when
- A `[logging]` TOML block (level/file/to_stdout/to_stderr/timestamps) is honored with
  precedence default < TOML < env < C-API, in all four entry points, unit-tested.
- oxygenstress/interception solver+guard messages go through the logger; the
  TEMPORARY-DELETE 777/888 debug is gone; the `log_unit_handle` consumer is gone (if it
  was the only one).
- `check-fast` 4/4 byte-identical; new test counted in `OK (N)`.
