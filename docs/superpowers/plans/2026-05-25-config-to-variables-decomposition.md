# config_to_variables Decomposition Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Break up the 1631-line `src/io/toml/config_to_variables.f90` adapter by moving each subsystem's seeding logic into a typed `state%X%init(config%Y, …)` method on the corresponding state record; the adapter becomes a thin orchestrator (`seed_state_from_config`), and `swap_init_from_loaded_config` is folded into `swap_init` so the call chain is `swap_main → swap_init → seed_state_from_config`.

**Architecture:** Strangler-pattern terminal cleanup. The TOML readers (`src/io/toml/read_*_toml.f90`) and typed `*_config_t` records stay as-is. We move the **config→state seeding** out of the adapter and onto each subsystem's typed state. Subsystems that already have a type-bound `init` get extended; those with free `<subsystem>_init` subs get promoted to type-bound. The orchestrator dispatches each subsystem's init in dependency order (mesh + timecontrol first, then everything else).

**Tech Stack:** Fortran 2003+ (gfortran, `-std=legacy`), Meson/Ninja build (`pixi run check-fast`), pFUnit unit tests, byte-identical regression-output fixtures.

**Required reading before starting:**
- Spec: `docs/superpowers/specs/2026-05-25-config-to-variables-decomposition-design.md` — the section "Current state — inventory" lists every line range to be moved.
- `src/io/toml/config_to_variables.f90` — the source.
- `src/state/swap_state.f90` — defines `swap_state_t` and the slots `state%timecontrol`, `state%atmosphere`, … (so you know the type names).

**Cross-cutting conventions for every task:**
- Branch is already `development`; commit directly. No new branches.
- Per the `feedback_state_schema_clean_rebuild` memory: **any change that adds/removes fields on a `state_*_t` requires `rm -rf builddir`** before `pixi run check-fast`. Always rebuild from scratch when in doubt.
- Commits use the tag `[GR-SEED 2026-05-25]` so the arc is grep-able later.
- Verification gate at end of every task: `pixi run check-fast` must report **748/748 pFUnit + 4/4 regression byte-identical** (`✓ hupselbrook`, `✓ surfacewater`, `✓ salinitystress`, `✓ grassgrowth`). Annual-stats fixture mismatch = bug, do not commit.
- If a test in `tests/unit/io/toml/test_config_to_variables.pf` references a symbol you retired, update the test in the same commit (precedent set by Phase 6 Step 2, `0a44e09`).

---

## File Structure

**Existing (modified):**
- `src/state/timecontrol_state.f90` — gain new `init(self, config_simulation, config_general)` type-bound proc.
- `src/state/atmosphere_state.f90` — extend existing `init(self, config)`.
- `src/state/solute_state.f90` — gain type-bound `init(self, config_solute, numlay)`; free `solute_init` retired.
- `src/state/nutrients_state.f90` — extend existing `init(self, nlay)` → `init(self, nlay, config_nut, pathwork)`.
- `src/state/crop_irrigation_state.f90` — extend stub `init(self)` → `init(self, config_irrigation, tstart, tend, mesh)`.
- `src/state/surfacewater_state.f90` — extend existing `init(self, config_sw, config_drain, numnod)`.
- `src/state/drainage_state.f90` — gain type-bound `init(self, config_drain, numnod)`; free `drainage_init` retired.
- `src/state/soilwater_state.f90` — gain type-bound `init(self, config_soil, config_bb, numnod, numlay, state_solute, state_mesh)`; free `soilwater_init` retired. Bottom-boundary case dispatch folded in.
- `src/state/tillage_state.f90` — gain type-bound `init(self, config_tillage, tend, numlay)`; free `tillage_init` retired.
- `src/state/crop_state.f90` — extend existing `init(self, crop_cfg)` → `init(self, crop_cfg, pathwork)`.
- `src/io/toml/config_to_variables.f90` → **rename** to `src/io/toml/seed_state_from_config.f90` at Task 11; module renamed to `seed_state_from_config_mod`; subroutine renamed to `seed_state_from_config`.
- `src/core/swap_mod.f90` — `swap_init_from_loaded_config` folded into `swap_init` at Task 12.
- `tests/unit/io/toml/test_config_to_variables.pf` → renamed to `test_seed_state_from_config.pf` at Task 11.
- `meson.build` / `tests/unit/meson.build` — source list entries updated.

**Deleted:**
- The free `<subsystem>_init` subroutines for solute, drainage, soilwater, tillage (their implementations move into the type-bound methods).
- `swap_init_from_loaded_config` (folded into `swap_init`).

**Created:**
- None — every new home is an existing state module.

---

## Task 1: state%timecontrol%init

**Files:**
- Modify: `src/state/timecontrol_state.f90` (add type-bound init + body)
- Modify: `src/io/toml/config_to_variables.f90:67-124` (replace section with single call)
- Test (existing — to update later): `tests/unit/io/toml/test_config_to_variables.pf` (no change yet — `state%timecontrol` assertions already use `state%X` paths)

**Spec sections covered:** "General + simulation" + "Simulation.numerical" + "Output switches".

- [ ] **Step 1: Add type-bound init proc declaration**

Edit `src/state/timecontrol_state.f90`. Inside the `type :: timecontrol_state_t` block (currently ending around line 130), before `end type timecontrol_state_t`, find or add a `contains` line. The end of the type currently looks like:

```fortran
   end type timecontrol_state_t
```

Add a `contains` clause if not present, and one procedure:

```fortran
   contains
      procedure :: init => timecontrol_state_init
   end type timecontrol_state_t
```

- [ ] **Step 2: Implement `timecontrol_state_init` in the same module**

Add to the `contains` section of `timecontrol_state_mod` (the body of the module after `contains` keyword at module level):

```fortran
   !> Seed runtime time-control state from typed config.
   !!
   !! Computes iyear/imonth from tstart via dtdpar() (mirrors readswap.f90:126-128).
   !! Allocates outdat/outdatint to MAOUT and zero-fills (consumed by
   !! timecontrol_advance to gate output dumps). When config_simulation%swmonth==1,
   !! populates outdatint with end-of-month dates via populate_outdatint_monthly
   !! and forces period/swres/swodat = 0 (legacy behaviour).
   subroutine timecontrol_state_init(self, config_simulation, config_general)
      use simulation_config_mod, only: simulation_config_t
      use general_config_mod,    only: general_config_t
      use swap_array_dimensions, only: maout
      class(timecontrol_state_t),  intent(inout) :: self
      type(simulation_config_t),   intent(in)    :: config_simulation
      type(general_config_t),      intent(in)    :: config_general
      integer :: datea_init(6)
      real    :: fsec_init

      ! General — screen-/result-file switches
      self%swscre = config_general%swscre

      ! Simulation — clock window + output cadence + numerical solver
      self%tstart    = config_simulation%tstart
      self%tend      = config_simulation%tend
      self%nprintday = config_simulation%nprintday
      self%period    = config_simulation%period
      self%swres     = config_simulation%swres
      self%swodat    = config_simulation%swodat

      ! Numerical (sub-record)
      self%dt    = config_simulation%numerical%dt
      self%dtmin = config_simulation%numerical%dtmin
      self%dtmax = config_simulation%numerical%dtmax
      self%MaxIt = config_simulation%numerical%MaxIt
      self%msteps = config_simulation%numerical%msteps
      ! Not in schema — defaults match legacy.
      self%MaxIterTime   = 0
      self%flMaxIterTime = .false.

      ! Derive iyear/imonth from tstart (mirrors readswap.f90:126-128).
      call dtdpar(self%tstart + 0.1d0, datea_init, fsec_init)
      self%iyear  = datea_init(1)
      self%imonth = datea_init(2)

      ! Output-date schedules — allocate to legacy cap, zero-fill.
      if (.not. allocated(self%outdat))    allocate(self%outdat(maout))
      if (.not. allocated(self%outdatint)) allocate(self%outdatint(maout))
      self%outdat    = 0.0_real64
      self%outdatint = 0.0_real64

      ! Monthly output: populate end-of-month dates + clobber daily-period switches.
      if (config_simulation%swmonth == 1) then
         call populate_outdatint_monthly(self%tend, self%iyear, self%imonth, &
                                         self%outdatint)
         self%period = 0
         self%swres  = 0
         self%swodat = 0
      end if

      ! Output switches: legacy forced to 0 by ADR 0009.
      self%swheader = 0
   end subroutine timecontrol_state_init
```

`dtdpar` is a free subroutine from `src/core/dtutil.f90` (no `use` needed).

- [ ] **Step 3: Move `populate_outdatint_monthly` from adapter into this module**

Cut the body of `populate_outdatint_monthly` (currently at `src/io/toml/config_to_variables.f90:1122-1155`). Paste it into `src/state/timecontrol_state.f90` immediately after `timecontrol_state_init`, with the same signature:

```fortran
   !> Populate `outdatint(:)` with end-of-month dates between tstart and tend.
   !! Mirrors readswap.f90:181-204 (swmonth==1 branch).
   subroutine populate_outdatint_monthly(tend, iyear, imonth, outdatint)
      real(real64), intent(in)    :: tend
      integer,      intent(in)    :: iyear
      integer,      intent(in)    :: imonth
      real(real64), intent(inout) :: outdatint(:)
      integer :: datea_om(6), i_om
      real    :: fsec_om
      real(real64) :: outdate_om

      datea_om = 0
      datea_om(1) = iyear
      datea_om(2) = imonth
      if (datea_om(2) < 12) then
         datea_om(2) = datea_om(2) + 1
      else
         datea_om(1) = datea_om(1) + 1
         datea_om(2) = 1
      end if
      datea_om(3) = 1
      fsec_om = 0.0
      call dtardp(datea_om, fsec_om, outdate_om)
      i_om = 0
      do while ((outdate_om - 1.0d0) < (tend + 0.1d0))
         i_om = i_om + 1
         outdatint(i_om) = outdate_om - 1.0d0
         if (datea_om(2) < 12) then
            datea_om(2) = datea_om(2) + 1
         else
            datea_om(1) = datea_om(1) + 1
            datea_om(2) = 1
         end if
         call dtardp(datea_om, fsec_om, outdate_om)
      end do
   end subroutine populate_outdatint_monthly
```

Make this `private` to the module (no `public ::` line for it).

- [ ] **Step 4: Replace adapter section with call to new init**

Edit `src/io/toml/config_to_variables.f90`. Locate the General+simulation block (lines ~67-124) plus the Simulation.numerical block (~126-136). Delete every line that writes `state%timecontrol%X = config%X`, the `dtdpar` derivation block, the `outdat`/`outdatint` allocation block, the `populate_outdatint_monthly` call, and the swheader assignment (line ~1112). Replace **all of it** with one line near the top of the body:

```fortran
      call state%timecontrol%init(config%simulation, config%general)
```

Also delete the swheader assignment near the bottom (search for `state%timecontrol%swheader = 0` around line 1112 and delete that line — it's now inside the init).

- [ ] **Step 5: Clean rebuild and verify**

Run:

```bash
rm -rf builddir && pixi run check-fast
```

Expected (final lines):

```
 OK
 (748 tests)
...
✓ hupselbrook: regression ok (annual stats match fixture)
✓ surfacewater: regression ok (annual stats match fixture)
✓ salinitystress: regression ok (annual stats match fixture)
✓ grassgrowth: regression ok (annual stats match fixture)
```

If any test fails, check `popular_outdatint_monthly` was moved (not duplicated), the swheader write was deleted exactly once, and the section in the adapter was fully replaced (no leftover `self%timecontrol%X = …` lines).

- [ ] **Step 6: Commit**

```bash
git add src/state/timecontrol_state.f90 src/io/toml/config_to_variables.f90
git commit -m "$(cat <<'EOF'
refactor(seed): Task 1 — state%timecontrol%init absorbs general+simulation seeding

[GR-SEED 2026-05-25] Move 12 fields of seeding (swscre, tstart, tend,
nprintday, period, swres, swodat, dt, dtmin, dtmax, MaxIt, msteps,
MaxIterTime, flMaxIterTime, swheader, plus iyear/imonth derivation,
plus outdat/outdatint allocation, plus monthly-output populate) out of
config_to_variables.f90 into state%timecontrol%init(config%simulation,
config%general).

populate_outdatint_monthly relocated as a module-private helper of
timecontrol_state_mod.

check-fast 5/5 byte-identical.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: state%atmosphere%init extension

**Files:**
- Modify: `src/state/atmosphere_state.f90` (extend `init` body)
- Modify: `src/io/toml/config_to_variables.f90:137-263` (replace meteo+snow block)

**Spec sections covered:** "Meteorology" (incl. CSV pre-loads + snow scalars + evaporation scalar).

- [ ] **Step 1: Locate the existing atmosphere_state_init**

In `src/state/atmosphere_state.f90`, find `subroutine atmosphere_state_init(self, config)` (around line 282). It currently snapshots scalar config fields onto `self`. We will extend its body — signature stays.

- [ ] **Step 2: Add meteo CSV pre-load logic to atmosphere_state_init**

Inside `atmosphere_state_init`, after the existing scalar copies, insert:

```fortran
      ! [GR-SEED 2026-05-25 Task 2] Meteo CSV caches — pre-load once at init.
      block
         use csv_reader_mod, only: read_csv_table
         use error_mod,      only: error_collection_t
         character(len=300) :: metfile_lc, csvpath
         character(len=9)   :: hdr(9)
         real(8), allocatable :: tbl(:,:)
         type(error_collection_t) :: errs
         integer :: r

         metfile_lc = ''
         if (allocated(config%meteo%metfile))  metfile_lc = config%meteo%metfile
         call lowerc(metfile_lc)

         if (index(trim(metfile_lc), '.csv') > 0) then
            hdr(1) = 'date     '
            hdr(2) = 'rad      '
            hdr(3) = 'tmin     '
            hdr(4) = 'tmax     '
            hdr(5) = 'hum      '
            hdr(6) = 'wind     '
            hdr(7) = 'rain     '
            hdr(8) = 'etref    '
            hdr(9) = 'wet      '
            csvpath = trim(config%general%pathatm) // trim(metfile_lc)
            call read_csv_table(trim(csvpath), hdr, tbl, errs)
            call errs%abort_if_fatal()
            self%nmetcsv = size(tbl, 1)
            if (allocated(self%metcsv_dat)) deallocate(self%metcsv_dat)
            allocate(self%metcsv_dat(self%nmetcsv, 9))
            do r = 1, self%nmetcsv
               self%metcsv_dat(r, :) = tbl(r, :)
            end do
         end if
      end block

      ! Detail meteo CSV (swmetdetail=1 + detail_file allocated).
      if (config%meteo%swmetdetail == 1 .and. allocated(config%meteo%detail_file)) then
         if (len_trim(config%meteo%detail_file) > 0) then
            block
               use csv_reader_mod, only: read_csv_table
               use error_mod,      only: error_collection_t
               character(len=300) :: csvpath
               character(len=8)   :: hdr(7)
               real(8), allocatable :: tbl(:,:)
               type(error_collection_t) :: errs
               integer :: r
               hdr(1) = 'datetime'
               hdr(2) = 'record  '
               hdr(3) = 'rad     '
               hdr(4) = 'temp    '
               hdr(5) = 'hum     '
               hdr(6) = 'wind    '
               hdr(7) = 'rain    '
               csvpath = trim(config%general%pathatm) // trim(config%meteo%detail_file)
               call read_csv_table(trim(csvpath), hdr, tbl, errs)
               call errs%abort_if_fatal()
               self%nmetcsv_det = size(tbl, 1)
               if (allocated(self%metcsv_det)) deallocate(self%metcsv_det)
               allocate(self%metcsv_det(self%nmetcsv_det, 7))
               do r = 1, self%nmetcsv_det
                  self%metcsv_det(r, :) = tbl(r, :)
               end do
            end block
         end if
      end if

      ! Rain events CSV (swrain=3 + rain_events_file allocated).
      if (config%meteo%swrain == 3 .and. allocated(config%meteo%rain_events_file)) then
         if (len_trim(config%meteo%rain_events_file) > 0) then
            block
               use csv_reader_mod, only: read_csv_table
               use error_mod,      only: error_collection_t
               character(len=300) :: csvpath
               character(len=8)   :: hdr(2)
               real(8), allocatable :: tbl(:,:)
               type(error_collection_t) :: errs
               integer :: r
               hdr(1) = 'datetime'
               hdr(2) = 'amount  '
               csvpath = trim(config%general%pathatm) // trim(config%meteo%rain_events_file)
               call read_csv_table(trim(csvpath), hdr, tbl, errs)
               call errs%abort_if_fatal()
               self%nraincsv = size(tbl, 1)
               if (allocated(self%raincsv_dat)) deallocate(self%raincsv_dat)
               allocate(self%raincsv_dat(self%nraincsv, 2))
               do r = 1, self%nraincsv
                  self%raincsv_dat(r, :) = tbl(r, :)
               end do
            end block
         end if
      end if

      ! Snow scalars (already partly here — keep one source of truth).
      self%TePrRain = config%meteo%snow%teprrain
      self%TePrSnow = config%meteo%snow%teprsnow
```

If `TePrRain`/`TePrSnow` lines already exist higher in the body, keep one copy only.

- [ ] **Step 3: Move evaporation scalar (state%crop%cfbs) into crop_state init prep**

The line `state%crop%cfbs = config%meteo%evaporation%cfbs` is currently in the meteo block. cfbs lives on `state%crop`, not on atmosphere. Move that single line into Task 10 (crop init); for now leave the line in the adapter outside the meteo block — DO NOT delete yet. Add a `! [SS-T10]` comment above it to flag the relocation.

- [ ] **Step 4: Replace adapter meteo section with single call**

In `src/io/toml/config_to_variables.f90`, locate the meteo block (around lines 137–263 — section header `! Meteorology (audit: 12 + evaporation + snow)`). The line `call state%atmosphere%init(config)` is **already called downstream from swap_mod**, NOT from this adapter. So our removal here just deletes the seeding work; we do NOT add a new call here.

Verify the call site exists at `src/core/swap_mod.f90:280` (search for `call state%atmosphere%init(config)`). It does. Good.

Now delete from `config_to_variables.f90`: lines 137–263 inclusive — the entire meteorology block. Replace with a single comment line:

```fortran
      ! [GR-SEED 2026-05-25 Task 2] Meteorology/snow/evap seeding moved to
      ! state%atmosphere%init(config) (called from swap_mod) and state%crop init.
```

- [ ] **Step 5: Clean rebuild and verify**

```bash
rm -rf builddir && pixi run check-fast
```

Expected: 748 pFUnit + 4/4 regression byte-identical.

If meteo data appears missing during regression (e.g., `arad(1)` is zero), check that `state%atmosphere%init(config)` is still called from `swap_mod.f90` (line ~280) and that the CSV pre-load runs to completion (add `print *, self%nmetcsv` in init to debug if needed).

- [ ] **Step 6: Commit**

```bash
git add src/state/atmosphere_state.f90 src/io/toml/config_to_variables.f90
git commit -m "$(cat <<'EOF'
refactor(seed): Task 2 — state%atmosphere%init absorbs meteo CSV pre-loads

[GR-SEED 2026-05-25] Meteo CSV caches (metcsv_dat/nmetcsv,
metcsv_det/nmetcsv_det, raincsv_dat/nraincsv) and snow scalars
(TePrRain/TePrSnow) now live inside state%atmosphere%init(config). The
existing call from swap_mod.f90 picks up the new work. Adapter's
meteorology block (lines 137-263) collapses to a single comment.

state%crop%cfbs assignment retained in adapter pending Task 10.

check-fast 5/5 byte-identical.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: state%solute%init (promote + extend)

**Files:**
- Modify: `src/state/solute_state.f90` (add type-bound init, replace free `solute_init`)
- Modify: `src/io/toml/config_to_variables.f90:883-957` (replace solute block)
- Modify: `src/core/swap_mod.f90:382` (`call solute_init(state)` → `call state%solute%init(...)`)

**Spec sections covered:** "Solute".

- [ ] **Step 1: Inspect current solute_init**

```bash
grep -n "subroutine solute_init\|end subroutine solute_init" src/solute/solute.f90
```

Expected: two line numbers bracketing the body of the current free subroutine. Read those lines.

- [ ] **Step 2: Add type-bound init declaration to solute_state**

Edit `src/state/solute_state.f90`. Inside `type :: solute_state_t`, before `end type`, add (or extend) a `contains` block:

```fortran
   contains
      procedure :: init => solute_state_init
   end type solute_state_t
```

- [ ] **Step 3: Implement `solute_state_init`**

Add after the type's `contains` at the module level (with `contains` keyword at module level):

```fortran
   !> Seed solute state from typed config + dimension args.
   !!
   !! Absorbs the per-layer-array seeding + cseeptab flatten previously in
   !! config_to_variables.f90:883-957 and the runtime zero-fill of state
   !! fields from the legacy free `solute_init(state)`.
   subroutine solute_state_init(self, config_solute, numlay)
      use solute_config_mod, only: solute_config_t
      class(solute_state_t),    intent(inout) :: self
      type(solute_config_t),    intent(in)    :: config_solute
      integer,                  intent(in)    :: numlay
      integer :: i

      ! Scalars
      self%swbotbc = config_solute%swbotbc
      self%cdrain  = config_solute%cdrain
      self%tscf    = config_solute%tscf
      self%bexp    = config_solute%bexp
      self%cref    = config_solute%cref
      self%cpre    = config_solute%cpre
      self%ddif    = config_solute%ddif
      self%frexp   = config_solute%frexp
      self%gampar  = config_solute%gampar
      self%daquif  = config_solute%daquif
      self%kfsat   = config_solute%kfsat
      self%decsat  = config_solute%decsat
      self%poros   = config_solute%poros
      self%rtheta  = config_solute%rtheta
      self%swbr    = config_solute%swbr

      ! Per-layer arrays — broadcast scalar to ldis(1) when no array provided.
      if (allocated(config_solute%ldis_array)) then
         do i = 1, min(size(config_solute%ldis_array), size(self%ldis))
            self%ldis(i) = config_solute%ldis_array(i)
         end do
      else
         self%ldis(1) = config_solute%ldis
      end if

      if (allocated(config_solute%kf)) then
         do i = 1, min(size(config_solute%kf), size(self%kf))
            self%kf(i) = config_solute%kf(i)
         end do
      end if
      if (allocated(config_solute%decpot)) then
         do i = 1, min(size(config_solute%decpot), size(self%decpot))
            self%decpot(i) = config_solute%decpot(i)
         end do
      end if
      if (allocated(config_solute%fdepth)) then
         do i = 1, min(size(config_solute%fdepth), size(self%fdepth))
            self%fdepth(i) = config_solute%fdepth(i)
         end do
      end if

      ! 2D cseeptab → interleaved afgen layout: (time, conc) pairs.
      if (allocated(config_solute%cseeptab)) then
         do i = 1, min(size(config_solute%cseeptab, 1), size(self%cseeptab) / 2)
            self%cseeptab(2*i - 1) = config_solute%cseeptab(i, 1)
            self%cseeptab(2*i)     = config_solute%cseeptab(i, 2)
         end do
      end if

      ! Runtime zero-fill (formerly free solute_init body — copy verbatim).
      ! [INSERT here the zero-fill lines from current src/solute/solute.f90
      !  subroutine solute_init body, with `state%solute%X` → `self%X`.]
   end subroutine solute_state_init
```

Open `src/solute/solute.f90`, find `subroutine solute_init(state)`, copy every `state%solute%X = …` line out of its body into the `! Runtime zero-fill` placeholder above, replacing each `state%solute%X` with `self%X`. The free `solute_init` is then deletable.

- [ ] **Step 4: Delete the free solute_init from solute.f90**

In `src/solute/solute.f90`, delete the entire `subroutine solute_init(state) … end subroutine solute_init` block. Also remove `solute_init` from the module's `public ::` list (if exported).

- [ ] **Step 5: Update the caller in swap_mod**

`src/core/swap_mod.f90:382` currently reads:

```fortran
   if (flSolute) call solute_init(state)
```

Replace with:

```fortran
   if (state%cfg%solute%swsolu == 1) call state%solute%init(config%solute, state%mesh%numlay)
```

(There is no `flSolute` left after Phase 6; the predicate is the config switch.) If `flSolute` is referenced and not yet retired, leave the predicate as `flSolute` for now and surface that mismatch as a separate cleanup — but check first:

```bash
grep -n "flSolute" src/core/swap_mod.f90
```

If grep shows any `flSolute` left, replace each with `state%cfg%solute%swsolu == 1` in the same commit.

Also drop the `use solute_mod, only: …, solute_init` import — only `solute` (the timestep proc) stays.

- [ ] **Step 6: Replace adapter solute section with no-op comment**

Delete `src/io/toml/config_to_variables.f90` lines 883–957 (the solute block — section header `! Solute (audit: 8 fields + 14 Phase 0 promoted fields)`). Replace with:

```fortran
      ! [GR-SEED 2026-05-25 Task 3] Solute seeding moved to state%solute%init
      ! (called from swap_mod).
```

- [ ] **Step 7: Clean rebuild and verify**

```bash
rm -rf builddir && pixi run check-fast
```

Expected: 748 pFUnit + 4/4 regression byte-identical.

If `salinitystress` regression drifts, check that the per-layer broadcast logic was preserved (the `else self%ldis(1) = config_solute%ldis` branch is easy to miss).

- [ ] **Step 8: Commit**

```bash
git add src/state/solute_state.f90 src/solute/solute.f90 src/core/swap_mod.f90 src/io/toml/config_to_variables.f90
git commit -m "$(cat <<'EOF'
refactor(seed): Task 3 — state%solute%init promoted to type-bound

[GR-SEED 2026-05-25] Free solute_init(state) retired; replaced with
type-bound state%solute%init(config%solute, numlay). Body absorbs the
solute scalars + per-layer array broadcasts + 2D cseeptab flatten
previously in config_to_variables.f90:883-957, plus the runtime
zero-fill formerly in src/solute/solute.f90.

Caller in swap_mod.f90:382 updated. Adapter's solute block collapses
to a single comment.

check-fast 5/5 byte-identical.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: state%nutrients%init (extend)

**Files:**
- Modify: `src/state/nutrients_state.f90` (extend `init` body, change signature)
- Modify: `src/io/toml/config_to_variables.f90:1357-1499` (relocate helpers)
- Modify: `src/core/swap_mod.f90:134` (update init call args)

**Spec sections covered:** Helpers `apply_nutrients` + `apply_nutrients_events`.

- [ ] **Step 1: Change init signature in nutrients_state**

Edit `src/state/nutrients_state.f90`. Locate `subroutine nutrients_state_init(self, nlay)`. Change the signature and body:

```fortran
   subroutine nutrients_state_init(self, nlay, config_nut, pathwork_in)
      use nutrients_config_mod, only: nutrients_config_t
      class(nutrients_state_t),    intent(inout) :: self
      integer,                     intent(in)    :: nlay
      type(nutrients_config_t),    intent(in)    :: config_nut
      character(len=*),            intent(in)    :: pathwork_in
      ! ... existing nlay-based zero-fill code goes here unchanged ...

      ! Then: nutrient config seeding (moved from apply_nutrients).
      call apply_nutrients_into_state(config_nut, pathwork_in)
   end subroutine nutrients_state_init
```

Keep the existing per-layer zero-fill at the top of the routine (whatever it currently does with `nlay`).

- [ ] **Step 2: Move `apply_nutrients` + `apply_nutrients_events` bodies into nutrients_state_mod**

Cut the two subroutines `apply_nutrients` and `apply_nutrients_events` from `src/io/toml/config_to_variables.f90` (lines ~1357-1499). Paste into `src/state/nutrients_state.f90` after `nutrients_state_init`, renaming `apply_nutrients` to `apply_nutrients_into_state`. Both are `private` to the module.

Adjust uses inside them: drop `use config_to_variables_mod, only: apply_nutrients_events`. The renamed routine should call its sibling directly (now in the same module).

The bodies write to `Wofost_Soil_Declarations` (a separate module). That's intentional — those WSN globals are out of scope for this arc. Keep the `use Wofost_Soil_Declarations, only: …` line at the top of `apply_nutrients_into_state`.

- [ ] **Step 3: Update swap_mod caller**

`src/core/swap_mod.f90:134`:

```fortran
   call state%nutrients%init(state%mesh%numlay)
```

Replace with:

```fortran
   call state%nutrients%init(state%mesh%numlay, config%nutrients, config%general%pathwork)
```

- [ ] **Step 4: Drop `call apply_nutrients` from adapter**

Find in `src/io/toml/config_to_variables.f90` the line `call apply_nutrients(config%nutrients, config%general%pathwork)` (around line 445). Delete it. The init above now does the same work, in the same place where `call state%nutrients%init(...)` is invoked from swap_mod.

Also remove `public :: apply_nutrients` and `public :: apply_nutrients_events` from the module's public list (the `module config_to_variables_mod` declaration block at the top).

- [ ] **Step 5: Update tests if any import the retired publics**

```bash
grep -rn "apply_nutrients\|apply_nutrients_events" tests/
```

If any test uses these directly, update its `use config_to_variables_mod, only: …` line to point at `nutrients_state_mod` and the renamed `apply_nutrients_into_state` (or just call the type-bound init in tests).

- [ ] **Step 6: Clean rebuild and verify**

```bash
rm -rf builddir && pixi run check-fast
```

Expected: 748 pFUnit + 4/4 regression byte-identical.

- [ ] **Step 7: Commit**

```bash
git add src/state/nutrients_state.f90 src/io/toml/config_to_variables.f90 src/core/swap_mod.f90 tests/
git commit -m "$(cat <<'EOF'
refactor(seed): Task 4 — state%nutrients%init absorbs nutrients seeding

[GR-SEED 2026-05-25] state%nutrients%init signature extended:
init(self, nlay) → init(self, nlay, config_nut, pathwork_in). Body
absorbs apply_nutrients + apply_nutrients_events helpers from the
adapter; both relocated into nutrients_state_mod as private helpers
(apply_nutrients_into_state).

Wofost_Soil_Declarations writes preserved (out-of-scope cleanup).
Caller in swap_mod.f90:134 updated.

check-fast 5/5 byte-identical.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: state%crop%irrigation%init (extend)

**Files:**
- Modify: `src/state/crop_irrigation_state.f90` (extend stub `init`)
- Modify: `src/io/toml/config_to_variables.f90:818-881, 1291-1629` (move helpers; drop blocks)
- Modify: `src/core/swap_mod.f90` if/where `state%crop%irrigation%init` is invoked

**Spec sections covered:** Fixed irrigation events block + helpers `apply_irrigation_ssdi`, `apply_ssdi_mode0`, `apply_ssdi_mode1`.

- [ ] **Step 1: Identify current init call site for crop%irrigation%init**

```bash
grep -n "state%crop%irrigation%init\|crop_irrigation_state_init" src/core/swap_mod.f90 src/state/crop_state.f90
```

If the call already happens from `state%crop%init` (the umbrella), Task 5's caller update lives there. If not, a new call must be added to swap_mod. Note which.

- [ ] **Step 2: Extend init signature in crop_irrigation_state**

Edit `src/state/crop_irrigation_state.f90`. Change `subroutine crop_irrigation_state_init(self)` to:

```fortran
   subroutine crop_irrigation_state_init(self, config_irrigation, tstart, tend, mesh)
      use irrigation_config_mod, only: irrigation_config_t
      use mesh_state_mod,        only: mesh_state_t
      class(crop_irrigation_state_t), intent(inout) :: self
      type(irrigation_config_t),      intent(in)    :: config_irrigation
      real(real64),                   intent(in)    :: tstart, tend
      type(mesh_state_t),             intent(in)    :: mesh

      ! Keep existing default-init stub behaviour (defaults set on type decl).

      ! Switch + counter
      self%swssdi = config_irrigation%swssdi
      self%nirri  = 1

      ! Fixed-irrigation events seeding (formerly adapter:818-881)
      call seed_fixed_irrigation(self, config_irrigation, tend)

      ! SSDI seeding (formerly apply_irrigation_ssdi)
      if (config_irrigation%swssdi == 1) then
         call apply_ssdi_seed(self, config_irrigation%ssdi, tstart, tend, mesh)
      end if
   end subroutine crop_irrigation_state_init
```

- [ ] **Step 3: Move fixed-irrigation block into a private helper**

In `src/state/crop_irrigation_state.f90` add (as private):

```fortran
   !> Seed fixed-irrigation event arrays. Mirrors the inline + CSV branches
   !! of the legacy adapter (config_to_variables.f90:818-881).
   subroutine seed_fixed_irrigation(self, config_irrigation, tend)
      use irrigation_config_mod, only: irrigation_config_t
      ! [INSERT here the body from current config_to_variables.f90 lines 818-881,
      !  with destinations switched from state%crop%irrigation%X → self%X,
      !  and the surrounding `if (allocated(config%irrigation%fixed_events)) then`
      !  preserved verbatim.]
   end subroutine seed_fixed_irrigation
```

Read `src/io/toml/config_to_variables.f90:818-881` to copy the inline + CSV branches verbatim into the helper.

- [ ] **Step 4: Move `apply_irrigation_ssdi` + `apply_ssdi_mode0` + `apply_ssdi_mode1` into crop_irrigation_state_mod**

Cut the three subroutines from `config_to_variables.f90` (lines ~1291-1354 for `apply_irrigation_ssdi`, 1502-1592 for `apply_ssdi_mode0`, 1599-1629 for `apply_ssdi_mode1`). Paste into `crop_irrigation_state_mod` as private helpers, renaming:

- `apply_irrigation_ssdi` → `apply_ssdi_seed(self, ssdi, tstart, tend, mesh)`
- The mode0/mode1 helpers keep their names but become private to the module.

Adjust the dispatch inside `apply_ssdi_seed`:

```fortran
   subroutine apply_ssdi_seed(self, ssdi, tstart, tend, mesh)
      use irrigation_config_mod, only: irrigation_ssdi_t
      use mesh_state_mod,        only: mesh_state_t
      type(crop_irrigation_state_t), intent(inout) :: self
      type(irrigation_ssdi_t),       intent(in)    :: ssdi
      real(real64),                  intent(in)    :: tstart, tend
      type(mesh_state_t),            intent(in)    :: mesh
      ! ncomp resolution code from current adapter — write into self instead of state.

      if (ssdi%schedule == 0) then
         call apply_ssdi_mode0(ssdi, ncomp, tstart, tend, self, pathwork_in)
      else
         call apply_ssdi_mode1(ssdi, ncomp, self)
      end if
   end subroutine apply_ssdi_seed
```

The signature change: `pathwork_in` must be threaded through if mode0 needs it for CSV reads. Pass it as a 6th arg to `apply_ssdi_seed` if so, and onward from `crop_irrigation_state_init`. Check `apply_ssdi_mode0` body — if it does CSV reads with `pathwork_in`, thread the arg.

- [ ] **Step 5: Update the crop%irrigation%init call**

In whatever file currently calls `state%crop%irrigation%init(...)` (likely `state%crop%init` in `src/state/crop_state.f90` — check), update the call to pass the new args:

```fortran
   call self%irrigation%init(config%irrigation, state%timecontrol%tstart, &
                             state%timecontrol%tend, state%mesh)
```

If the call must be moved from one place to another (e.g., from swap_mod to crop_state.init), reflect that in the commit message.

- [ ] **Step 6: Drop apply_irrigation_ssdi call from adapter**

In `src/io/toml/config_to_variables.f90`, find `if (config%irrigation%swssdi == 1) call apply_irrigation_ssdi(...)` (around line 437). Delete it — the init now handles this.

Also delete the fixed-irrigation block (lines 818-881) — replace with comment:

```fortran
      ! [GR-SEED 2026-05-25 Task 5] Fixed-irrigation seeding moved to
      ! state%crop%irrigation%init (called from state%crop%init).
```

Remove `apply_irrigation_ssdi` from the module's `public ::` list. The test `tests/unit/io/toml/test_apply_irrigation_ssdi.pf` imports `apply_irrigation_ssdi` — update its `use config_to_variables_mod, only: apply_irrigation_ssdi` to import from the new home (or re-target the test to call `state%crop%irrigation%init` directly).

- [ ] **Step 7: Clean rebuild and verify**

```bash
rm -rf builddir && pixi run check-fast
```

Expected: 748 pFUnit + 4/4 regression byte-identical, including the SSDI tests in `test_apply_irrigation_ssdi.pf`.

If SSDI mode0/mode1 tests fail, check the renamed-helper imports in the test file and the `pathwork_in` threading.

- [ ] **Step 8: Commit**

```bash
git add src/state/crop_irrigation_state.f90 src/io/toml/config_to_variables.f90 src/core/swap_mod.f90 src/state/crop_state.f90 tests/unit/io/toml/test_apply_irrigation_ssdi.pf
git commit -m "$(cat <<'EOF'
refactor(seed): Task 5 — state%crop%irrigation%init absorbs irrigation seeding

[GR-SEED 2026-05-25] state%crop%irrigation%init extended to
init(self, config_irrigation, tstart, tend, mesh). Absorbs:
- fixed-irrigation event seeding (inline + CSV branches);
- SSDI seeding (apply_irrigation_ssdi + mode0 + mode1) as private
  helpers (apply_ssdi_seed, apply_ssdi_mode0, apply_ssdi_mode1).

Adapter's irrigation block (818-881) and SSDI helpers (1291-1629)
deleted. test_apply_irrigation_ssdi.pf rewired to import from
crop_irrigation_state_mod.

check-fast 5/5 byte-identical.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: state%surfacewater%init (extend)

**Files:**
- Modify: `src/state/surfacewater_state.f90` (extend init body)
- Modify: `src/io/toml/config_to_variables.f90:959-1046` (drop block)

**Spec sections covered:** "Surface water (management periods)".

- [ ] **Step 1: Extend surfacewater_state_init body**

The existing signature is `subroutine surfacewater_state_init(self, config_sw, config_drain, numnod)`. Keep it. Inside the body, after existing seeding, append:

```fortran
      ! [GR-SEED 2026-05-25 Task 6] Management-period arrays from config_sw.
      block
         use swap_array_dimensions, only: mamp
         integer :: i

         self%osswlm = config_sw%osswlm
         self%nmper  = config_sw%nmper
         self%swqhr  = config_sw%swqhr

         if (allocated(config_sw%impend)) then
            if (.not. allocated(self%impend))  allocate(self%impend(mamp))
            do i = 1, min(size(config_sw%impend), size(self%impend))
               self%impend(i) = config_sw%impend(i)
            end do
         end if
         if (allocated(config_sw%swman)) then
            if (.not. allocated(self%swman))   allocate(self%swman(mamp))
            do i = 1, min(size(config_sw%swman), size(self%swman))
               self%swman(i) = config_sw%swman(i)
            end do
         end if
         if (allocated(config_sw%wscap)) then
            if (.not. allocated(self%wscap))   allocate(self%wscap(mamp))
            do i = 1, min(size(config_sw%wscap), size(self%wscap))
               self%wscap(i) = config_sw%wscap(i)
            end do
         end if
         if (allocated(config_sw%wldip)) then
            if (.not. allocated(self%wldip))   allocate(self%wldip(mamp))
            do i = 1, min(size(config_sw%wldip), size(self%wldip))
               self%wldip(i) = config_sw%wldip(i)
            end do
         end if
         if (allocated(config_sw%intwl)) then
            if (.not. allocated(self%intwl))   allocate(self%intwl(mamp))
            do i = 1, min(size(config_sw%intwl), size(self%intwl))
               self%intwl(i) = config_sw%intwl(i)
            end do
         end if
         if (allocated(config_sw%hbweir)) then
            if (.not. allocated(self%hbweir))  allocate(self%hbweir(mamp))
            do i = 1, min(size(config_sw%hbweir), size(self%hbweir))
               self%hbweir(i) = config_sw%hbweir(i)
            end do
         end if
         ! alphaw, betaw if present in schema
         if (allocated(config_sw%alphaw)) then
            if (.not. allocated(self%alphaw))  allocate(self%alphaw(mamp))
            do i = 1, min(size(config_sw%alphaw), size(self%alphaw))
               self%alphaw(i) = config_sw%alphaw(i)
            end do
         end if
         if (allocated(config_sw%betaw)) then
            if (.not. allocated(self%betaw))   allocate(self%betaw(mamp))
            do i = 1, min(size(config_sw%betaw), size(self%betaw))
               self%betaw(i) = config_sw%betaw(i)
            end do
         end if
      end block
```

- [ ] **Step 2: Drop adapter surface-water block**

In `src/io/toml/config_to_variables.f90`, delete lines ~959-1046 (the surface_water block). Replace with:

```fortran
      ! [GR-SEED 2026-05-25 Task 6] Surface-water management-period seeding
      ! moved to state%surfacewater%init (called from swap_mod).
```

- [ ] **Step 3: Clean rebuild and verify**

```bash
rm -rf builddir && pixi run check-fast
```

Expected: 748 pFUnit + 4/4 regression byte-identical, including the `surfacewater` regression case (which is the most relevant here).

- [ ] **Step 4: Commit**

```bash
git add src/state/surfacewater_state.f90 src/io/toml/config_to_variables.f90
git commit -m "$(cat <<'EOF'
refactor(seed): Task 6 — surfacewater management-periods seeding into init

[GR-SEED 2026-05-25] state%surfacewater%init absorbs the 8 management-
period arrays (impend, swman, wscap, wldip, intwl, hbweir, alphaw,
betaw) plus osswlm/nmper/swqhr scalars. Adapter's surface-water block
(959-1046) collapses to one comment.

check-fast 5/5 byte-identical.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: state%drainage%init (promote + extend)

**Files:**
- Modify: `src/state/drainage_state.f90` (add type-bound init)
- Modify: `src/drainage/drainage.f90` (retire free `drainage_init`)
- Modify: `src/io/toml/config_to_variables.f90:265-418` (drop drainage block)
- Modify: `src/core/swap_mod.f90` (update `call drainage_init(state, config)`)

**Spec sections covered:** "Drainage" (incl. DRAMET=2 dispatch, swdtyp/swallo arrays, owltab CSV, surface_runoff sub-section).

- [ ] **Step 1: Inspect current drainage_init**

```bash
grep -n "subroutine drainage_init\|end subroutine drainage_init" src/drainage/drainage.f90
```

Read the lines bracketed by these results. Note the variables it sets (drainl, wetper, etc.).

- [ ] **Step 2: Add type-bound init declaration to drainage_state**

Edit `src/state/drainage_state.f90`. Inside `type :: drainage_state_t`, add:

```fortran
   contains
      procedure :: init => drainage_state_init
   end type drainage_state_t
```

- [ ] **Step 3: Implement drainage_state_init**

In the module's `contains` section, add `drainage_state_init`. The body has three pieces:

```fortran
   subroutine drainage_state_init(self, config_drain, numnod)
      use drainage_config_mod,   only: drainage_config_t
      use swap_array_dimensions, only: madr, maowl
      class(drainage_state_t),  intent(inout) :: self
      type(drainage_config_t),  intent(in)    :: config_drain
      integer,                  intent(in)    :: numnod
      integer :: i

      ! Piece A: scalars + DRAMET=2 ipos chain (from adapter:265-310)
      self%dramet   = config_drain%dramet
      self%swdislay = config_drain%swdislay
      self%basegw   = config_drain%basegw
      self%entres   = config_drain%entres
      if (config_drain%dramet == 2) then
         self%ipos  = config_drain%ipos
         self%khtop = config_drain%khtop
         if (config_drain%ipos >= 3) then
            self%khbot = config_drain%khbot
            self%zintf = config_drain%zintf
         end if
         if (config_drain%ipos >= 4) then
            self%kvtop = config_drain%kvtop
            self%kvbot = config_drain%kvbot
         end if
         if (config_drain%ipos == 5) then
            self%geofac = config_drain%geofac
         end if
      end if

      ! Piece B: per-level swdtyp + swallo arrays
      if (allocated(config_drain%swdtyp)) then
         if (.not. allocated(self%swdtyp)) then
            allocate(self%swdtyp(size(config_drain%swdtyp)))
            self%swdtyp = 0
         end if
         do i = 1, size(config_drain%swdtyp)
            self%swdtyp(i) = config_drain%swdtyp(i)
         end do
      end if
      if (allocated(config_drain%swallo)) then
         if (.not. allocated(self%swallo)) then
            allocate(self%swallo(size(config_drain%swallo)))
            self%swallo = 0
         end if
         do i = 1, size(config_drain%swallo)
            self%swallo(i) = config_drain%swallo(i)
         end do
      end if
      self%swliminf = config_drain%swliminf

      ! Piece C: owltab CSV pre-load (from adapter:330-355)
      if (allocated(config_drain%owltab_file)) then
         block
            use iso_fortran_env, only: real64
            use csv_reader_mod, only: read_csv_table
            use error_mod, only: error_collection_t
            real(real64), allocatable :: csv_table(:,:)
            type(error_collection_t)  :: csv_errs
            integer :: lev, nrows, k
            character(len=6) :: hdr(2)
            hdr(1) = 'date  '
            hdr(2) = 'level '
            if (.not. allocated(self%owltab)) then
               allocate(self%owltab(config_drain%nrlevs, 2*MAOWL))
               self%owltab = 0.0_real64
            end if
            do lev = 1, size(config_drain%owltab_file)
               if (len_trim(config_drain%owltab_file(lev)) == 0) cycle
               call read_csv_table(trim(config_drain%owltab_file(lev)), hdr, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               if (csv_errs%count() == 0) then
                  nrows = size(csv_table, 1)
                  self%nowltab(lev) = nrows
                  do k = 1, nrows
                     self%owltab(lev, 2*k-1) = csv_table(k, 1)
                     self%owltab(lev, 2*k)   = csv_table(k, 2)
                  end do
               end if
            end do
         end block
      end if

      ! Piece D: surface_runoff sub-section (from adapter:370-418)
      ! ... copy the swtopdislay/ftopdislay/NumLevRapDra block from the adapter
      ! and replace state%drainage%X → self%X.

      ! Piece E: legacy drainage_init body (drainl(1), wetper(1), etc. — copy
      ! from src/drainage/drainage.f90's current drainage_init body).
   end subroutine drainage_state_init
```

Open `src/drainage/drainage.f90`, find `drainage_init`, and copy the body lines (which set up runtime state arrays like `drainl`, `wetper`, allocations) into Piece E above, replacing `state%drainage%X` with `self%X`. Also open `src/io/toml/config_to_variables.f90:370-418` (the `surface_runoff` block) and copy that into Piece D verbatim, replacing `state%drainage%X` with `self%X` and adding `numnod` as the `MADR` cap source where needed.

- [ ] **Step 4: Delete free drainage_init**

In `src/drainage/drainage.f90`, delete the `subroutine drainage_init(state, config) … end subroutine drainage_init` block. Remove `drainage_init` from the module's `public ::` list.

- [ ] **Step 5: Update swap_mod caller**

`src/core/swap_mod.f90:358` (currently `call drainage_init(state, config)`):

```fortran
   call state%drainage%init(config%drain, state%mesh%numnod)
```

Drop `use drainage_mod, only: drainage_init` import.

- [ ] **Step 6: Drop adapter drainage block**

Delete `src/io/toml/config_to_variables.f90` lines 265-418 (entire drainage block — section headers `! Drainage (audit: 20 fields ...)`). Replace with:

```fortran
      ! [GR-SEED 2026-05-25 Task 7] Drainage seeding moved to state%drainage%init
      ! (called from swap_mod).
      ! state%surfacewater%swdra is seeded from config%drain%swdra in
      ! state%surfacewater%init.
```

Note: `state%surfacewater%swdra = config%drain%swdra` (line 251 originally) should already be folded into Task 6's surfacewater init. If not, fold it in now and commit that as part of this task.

- [ ] **Step 7: Clean rebuild and verify**

```bash
rm -rf builddir && pixi run check-fast
```

Expected: 748 pFUnit + 4/4 regression byte-identical.

If `surfacewater` regression drifts, the swdra seeding might be misordered (state%surfacewater%init must run before any surfacewater step that reads swdra).

- [ ] **Step 8: Commit**

```bash
git add src/state/drainage_state.f90 src/drainage/drainage.f90 src/io/toml/config_to_variables.f90 src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
refactor(seed): Task 7 — state%drainage%init promoted to type-bound

[GR-SEED 2026-05-25] Free drainage_init(state, config) retired; new
state%drainage%init(config_drain, numnod) type-bound proc absorbs:
- DRAMET=2 ipos chain + scalar seeding (basegw/entres/swdislay);
- per-level swdtyp/swallo arrays;
- owltab CSV pre-load + allocation;
- surface_runoff sub-section (swtopdislay/ftopdislay/NumLevRapDra);
- legacy drainage_init runtime zero-fill body.

Adapter's drainage block (265-418) collapses to one comment.

check-fast 5/5 byte-identical.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: state%soilwater%init (promote + extend; absorbs bottom_boundary)

**Files:**
- Modify: `src/state/soilwater_state.f90` (add type-bound init, absorb free `soilwater_init`)
- Modify: `src/io/toml/config_to_variables.f90:420-795` (drop soil + bottom_boundary blocks)
- Modify: `src/core/swap_mod.f90` (update caller)

**Spec sections covered:** "Soil" + "Bottom boundary" (largest task).

- [ ] **Step 1: Inspect current soilwater_init**

```bash
grep -n "subroutine soilwater_init\|end subroutine soilwater_init" src/state/soilwater_state.f90
```

Read lines and note current signature `(self, numnod, numlay)`.

- [ ] **Step 2: Add type-bound init wrapper**

Edit `src/state/soilwater_state.f90`. Add a `contains` block inside `type :: soilwater_state_t` with:

```fortran
   contains
      procedure :: init => soilwater_state_init_v2
```

(Use a temporary suffix `_v2` if there's a name collision with the existing free sub — we'll rename to `init` at the end of this task once the free sub is fully gone.)

Implement the new init with extended signature:

```fortran
   subroutine soilwater_state_init_v2(self, config_soil, config_bb, numnod, numlay)
      use soil_config_mod,           only: soil_config_t
      use bottom_boundary_config_mod, only: bottom_boundary_config_t
      class(soilwater_state_t),       intent(inout) :: self
      type(soil_config_t),            intent(in)    :: config_soil
      type(bottom_boundary_config_t), intent(in)    :: config_bb
      integer,                        intent(in)    :: numnod
      integer,                        intent(in)    :: numlay

      ! Piece A: current soilwater_init body (the runtime zero-fill / alloc
      ! based on numnod/numlay) — copy verbatim from the existing free sub.

      ! Piece B: soil seeding from adapter:420-595 (swsophy/swinco/bdens/swfrost,
      ! initial-conditions h_file CSV in swinco=3 branch, etc.).
      ! Replace state%soilwater%X with self%X.

      ! Piece C: bottom_boundary case dispatch from adapter:597-795.
      ! Folded in here because all destinations are self%X (gwltab, qbotab,
      ! haqtab, hbotab).
      call seed_bottom_boundary(self, config_bb)
   end subroutine soilwater_state_init_v2
```

Copy Piece A from the current free `soilwater_init` body. Copy Piece B from `config_to_variables.f90:420-595`. Piece C is delegated to a new private helper below.

- [ ] **Step 3: Add private `seed_bottom_boundary` helper**

In the same module:

```fortran
   subroutine seed_bottom_boundary(self, config_bb)
      use bottom_boundary_config_mod, only: bottom_boundary_config_t
      type(soilwater_state_t),        intent(inout) :: self
      type(bottom_boundary_config_t), intent(in)    :: config_bb

      select case (config_bb%swbotb)
      case (1)
         ! gwl_file CSV → self%gwltab (copy from adapter:609-640)
      case (2)
         ! qbot2_file CSV → self%qbotab (copy from adapter:642-672)
      case (3)
         ! haquif_file CSV → self%haqtab + optional qbot4 (copy from adapter:679-740)
      case (4)
         ! qhbot_file CSV → self%qbotab (copy from adapter:742-775)
      case (5)
         ! hbot5_file CSV → self%hbotab (copy from adapter:777-795)
      case default
         ! cases 6,7,8 — no seeding
      end select
   end subroutine seed_bottom_boundary
```

Open `src/io/toml/config_to_variables.f90:597-795` and copy each case branch into the matching `case ()` in the helper, replacing `state%soilwater%X` with `self%X`.

- [ ] **Step 4: Delete free `soilwater_init`**

After Piece A is folded into `soilwater_state_init_v2`, delete the entire `subroutine soilwater_init(...) … end subroutine soilwater_init` block. Rename `soilwater_state_init_v2` back to `soilwater_state_init` and the procedure binding to `procedure :: init => soilwater_state_init`. Remove `soilwater_init` from any `public ::` list.

- [ ] **Step 5: Update swap_mod caller**

`src/core/swap_mod.f90:133` (currently `call soilwater_init(state%soilwater, state%mesh%numnod, state%mesh%numlay)`):

```fortran
   call state%soilwater%init(config%soil, config%bottom_boundary, &
                              state%mesh%numnod, state%mesh%numlay)
```

Drop `use soilwater_state_mod, only: soilwater_init` import.

- [ ] **Step 6: Drop adapter soil + bottom_boundary blocks**

Delete `src/io/toml/config_to_variables.f90` lines 420-795 entirely. Replace with:

```fortran
      ! [GR-SEED 2026-05-25 Task 8] Soil + bottom_boundary seeding moved to
      ! state%soilwater%init (called from swap_mod).
```

- [ ] **Step 7: Clean rebuild and verify**

```bash
rm -rf builddir && pixi run check-fast
```

Expected: 748 pFUnit + 4/4 regression byte-identical. This is the highest-risk task — if any regression case fails, the most likely culprits are:

1. Init ordering — `state%soilwater%init` must run AFTER `state%mesh` is built (so `numnod`/`numlay` are valid) but BEFORE any compute call that reads `self%gwltab` etc.
2. The swinco=3 warm-restart h_file CSV reader was previously also writing into `config%soil%initial%z_init` (in-place mutation). If you preserved that, mark `config_soil` as `intent(inout)`.

- [ ] **Step 8: Commit**

```bash
git add src/state/soilwater_state.f90 src/io/toml/config_to_variables.f90 src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
refactor(seed): Task 8 — state%soilwater%init absorbs soil + bottom_boundary

[GR-SEED 2026-05-25] Free soilwater_init retired; new type-bound
state%soilwater%init(config_soil, config_bb, numnod, numlay) absorbs:
- the legacy runtime zero-fill / allocation logic from soilwater_init;
- soil seeding (swsophy/swinco/bdens/swfrost + swinco=3 warm restart
  h_file CSV) from adapter:420-595;
- bottom_boundary case dispatch (swbotb=1..5 CSV pre-loads) from
  adapter:597-795 as private seed_bottom_boundary helper.

Adapter's soil + bottom_boundary blocks (420-795) collapse to one
comment.

check-fast 5/5 byte-identical.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: state%tillage%init (promote + extend)

**Files:**
- Modify: `src/state/tillage_state.f90` (add type-bound init)
- Modify: `src/crop/tillage.f90` or `src/state/tillage_state.f90` (retire free `tillage_init`)
- Modify: `src/io/toml/config_to_variables.f90:1202-1288, 425` (drop helper + adapter call)
- Modify: `src/core/swap_mod.f90:340` (update caller)

**Spec sections covered:** Helper `apply_soil_tillage`.

- [ ] **Step 1: Inspect current tillage_init**

```bash
grep -n "subroutine tillage_init\|end subroutine tillage_init" src/crop/tillage.f90 src/state/tillage_state.f90
```

Note the file and the current signature `(self, numlay)`.

- [ ] **Step 2: Add type-bound init declaration**

Edit `src/state/tillage_state.f90`. Add `contains / procedure :: init => tillage_state_init` to the type. Implement:

```fortran
   subroutine tillage_state_init(self, config_tillage, tend, numlay)
      use soil_config_mod, only: soil_tillage_t  ! verify exact type name
      class(tillage_state_t),  intent(inout) :: self
      type(soil_tillage_t),    intent(in)    :: config_tillage
      real(real64),            intent(in)    :: tend
      integer,                 intent(in)    :: numlay

      ! Piece A: current tillage_init body (zero-fill / nlay-based allocation).
      ! Copy verbatim, replacing self%X.

      ! Piece B: apply_soil_tillage body — config_tillage → self.
      ! Copy from config_to_variables.f90:1202-1288 with arg renames.
   end subroutine tillage_state_init
```

Look up the exact name and module of `soil_tillage_t` with:

```bash
grep -n "type :: soil_tillage_t\|type :: tillage_config_t" src/config/soil_config.f90
```

Use whichever is correct (likely `soil_tillage_t` based on Phase 6's `apply_soil_tillage(config%soil%tillage, …)`).

- [ ] **Step 3: Move `apply_soil_tillage` body into tillage_state_init**

Open `src/io/toml/config_to_variables.f90:1202-1288`. Copy the body into Piece B above, replacing `state%tillage%X` with `self%X` and adjusting the signature variables (`tillage`, `tend` → match init args).

The helper `parse_iso_date_to_days1900` is used inside this body. **Relocate it** to `src/io/toml/toml_field_helpers_mod.f90` (a shared utilities home). The helper is small (~14 lines). Add it as a public sub in that module so the init can `use toml_field_helpers_mod, only: parse_iso_date_to_days1900`.

- [ ] **Step 4: Delete free `tillage_init`**

Delete the `subroutine tillage_init(...) … end subroutine tillage_init` block from its current location (likely `src/crop/tillage.f90` or `src/state/tillage_state.f90`). Remove from `public ::` list.

- [ ] **Step 5: Update swap_mod caller**

`src/core/swap_mod.f90:340` (currently `call tillage_init(state%tillage, state%mesh%numlay)`):

```fortran
   call state%tillage%init(config%soil%tillage, state%timecontrol%tend, state%mesh%numlay)
```

- [ ] **Step 6: Drop apply_soil_tillage call from adapter**

In `src/io/toml/config_to_variables.f90`, find:

```fortran
      if (config%soil%swtill == 1) call apply_soil_tillage(config%soil%tillage, state%timecontrol%tend, state)
```

Delete it. Also delete the `subroutine apply_soil_tillage … end subroutine apply_soil_tillage` block (lines 1202-1288). Remove `apply_soil_tillage` from the module's `public ::` list.

Replace with:

```fortran
      ! [GR-SEED 2026-05-25 Task 9] Tillage seeding moved to state%tillage%init
      ! (called from swap_mod). The (swtill == 1) gate now lives in the init.
```

The init should now check `if (config_tillage%swtill /= 1) return` near the top so the gate semantics are preserved.

- [ ] **Step 7: Update tests**

```bash
grep -rn "apply_soil_tillage" tests/
```

Update any test (e.g., `tests/unit/io/toml/test_apply_soil_tillage.pf`) to call `state%tillage%init(...)` directly instead of `apply_soil_tillage`. The test file likely needs to construct a `state%timecontrol` to access `tend` — match the precedent in `test_apply_irrigation_ssdi.pf`.

- [ ] **Step 8: Clean rebuild and verify**

```bash
rm -rf builddir && pixi run check-fast
```

Expected: 748 pFUnit + 4/4 regression byte-identical.

- [ ] **Step 9: Commit**

```bash
git add src/state/tillage_state.f90 src/crop/tillage.f90 src/io/toml/toml_field_helpers.f90 src/core/swap_mod.f90 src/io/toml/config_to_variables.f90 tests/
git commit -m "$(cat <<'EOF'
refactor(seed): Task 9 — state%tillage%init promoted; apply_soil_tillage folded in

[GR-SEED 2026-05-25] Free tillage_init retired; new type-bound
state%tillage%init(config_tillage, tend, numlay) absorbs apply_soil_tillage
from the adapter (config_to_variables.f90:1202-1288). The swtill==1 gate
moves into the init body. parse_iso_date_to_days1900 helper relocated to
toml_field_helpers_mod for shared use.

check-fast 5/5 byte-identical.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: state%crop%init (extend)

**Files:**
- Modify: `src/state/crop_state.f90` (extend init)
- Modify: `src/io/toml/config_to_variables.f90:1048-1087, 250-251` (drop crop + evap.cfbs)
- Modify: `src/core/swap_mod.f90` (update caller if needed)

**Spec sections covered:** "Crop" + the deferred `state%crop%cfbs` from Task 2.

- [ ] **Step 1: Inspect current crop_state_init signature**

```bash
grep -n "subroutine crop_state_init" src/state/crop_state.f90
```

Note signature `(self, crop_cfg)`.

- [ ] **Step 2: Extend crop_state_init**

Edit `src/state/crop_state.f90`. Change signature to:

```fortran
   subroutine crop_state_init(self, crop_cfg, meteo_cfg, pathwork_in)
      use crop_config_mod,         only: crop_config_t
      use meteorology_config_mod,  only: meteorology_config_t
      class(crop_state_t),         intent(inout) :: self
      type(crop_config_t), target, intent(in)    :: crop_cfg
      type(meteorology_config_t),  intent(in)    :: meteo_cfg
      character(len=*),            intent(in)    :: pathwork_in
      integer :: i

      ! existing init body — preserved

      ! [Task 10] flCropReadFile / flCropOpenFile + croptype + cfbs.
      if (crop_cfg%swcrop == 1) then
         self%common%flCropReadFile = .true.
         self%common%flCropOpenFile = .true.
      end if

      if (allocated(crop_cfg%rotation_type)) then
         if (.not. allocated(self%common%croptype)) then
            allocate(self%common%croptype(size(crop_cfg%rotation_type)))
         end if
         do i = 1, size(crop_cfg%rotation_type)
            self%common%croptype(i) = crop_cfg%rotation_type(i)
         end do
      end if

      ! cfbs — deferred from Task 2 (evaporation belongs to crop transpiration setup).
      self%cfbs = meteo_cfg%evaporation%cfbs

      ! crop_config_global legacy pointer — set here as the canonical seeding step.
      crop_config_global => crop_cfg
   end subroutine crop_state_init
```

`crop_config_global` is a module-level pointer in `crop_config_mod`. Read its declaration first to confirm:

```bash
grep -n "crop_config_global" src/config/crop_config.f90
```

If not visible from `crop_state_mod`, add `use crop_config_mod, only: crop_config_global` to the head of `crop_state_init`.

- [ ] **Step 3: Update caller of crop init**

Find where `state%crop%init(...)` is called:

```bash
grep -n "state%crop%init\|crop_state_init" src/core/swap_mod.f90 src/io/toml/config_to_variables.f90
```

If the call lives in `swap_init_from_loaded_config`, update its args to:

```fortran
   call state%crop%init(config%crop, config%meteo, config%general%pathwork)
```

If the call is currently inside `config_to_variables` (with a different shape), relocate the call to `swap_mod` and align args.

- [ ] **Step 4: Drop adapter crop + evap.cfbs lines**

In `src/io/toml/config_to_variables.f90`:
- Delete the meteo evaporation `state%crop%cfbs = config%meteo%evaporation%cfbs` line (left from Task 2).
- Delete lines 1048-1087 (the crop block — section header `! Crop ... Phase 0`).

Replace with:

```fortran
      ! [GR-SEED 2026-05-25 Task 10] Crop seeding (incl. evaporation%cfbs and
      ! the crop_config_global legacy pointer) moved to state%crop%init.
```

- [ ] **Step 5: Clean rebuild and verify**

```bash
rm -rf builddir && pixi run check-fast
```

Expected: 748 pFUnit + 4/4 regression byte-identical.

- [ ] **Step 6: Commit**

```bash
git add src/state/crop_state.f90 src/io/toml/config_to_variables.f90 src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
refactor(seed): Task 10 — state%crop%init absorbs crop seeding + cfbs + global pointer

[GR-SEED 2026-05-25] state%crop%init signature extended:
init(self, crop_cfg) → init(self, crop_cfg, meteo_cfg, pathwork_in).
Absorbs:
- flCropReadFile/flCropOpenFile flags (swcrop==1 gate);
- rotation_type → state%crop%common%croptype copy;
- state%crop%cfbs from config%meteo%evaporation%cfbs (deferred from Task 2);
- crop_config_global legacy pointer assignment.

Adapter's crop block (1048-1087) and the deferred cfbs line collapse
to one comment.

check-fast 5/5 byte-identical.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 11: Rename to seed_state_from_config

**Files:**
- Rename: `src/io/toml/config_to_variables.f90` → `src/io/toml/seed_state_from_config.f90`
- Modify: the renamed file (module + subroutine names)
- Rename: `tests/unit/io/toml/test_config_to_variables.pf` → `tests/unit/io/toml/test_seed_state_from_config.pf`
- Modify: `src/core/swap_mod.f90` (update `use` + call)
- Modify: `meson.build`, `tests/unit/meson.build`

**Spec sections covered:** Task 11.

- [ ] **Step 1: Verify the file is now thin**

```bash
wc -l src/io/toml/config_to_variables.f90
```

Expected: < 100 lines (most of the body is now comments + the orchestrator).

If still large, identify which Task's section wasn't fully cut over and complete that work first.

- [ ] **Step 2: Audit the orchestrator body**

```bash
sed -n '/subroutine config_to_variables/,/end subroutine config_to_variables/p' src/io/toml/config_to_variables.f90
```

It should now look like:

```fortran
   subroutine config_to_variables(config, state)
      type(swap_config_t), target, intent(inout) :: config
      type(swap_state_t),          intent(inout) :: state

      state%cfg => config

      ! All subsystem inits are called from swap_mod (where state is already
      ! fully wired). This routine is intentionally a near-no-op stub now —
      ! retained only so the swap_mod caller continues to compile through
      ! the rename in Task 11.
   end subroutine config_to_variables
```

If there are leftover writes, finish moving them first (revisit the earlier task that should have absorbed them).

- [ ] **Step 3: Rename the file with git mv**

```bash
git mv src/io/toml/config_to_variables.f90 src/io/toml/seed_state_from_config.f90
```

- [ ] **Step 4: Edit module + subroutine names inside the renamed file**

```fortran
! At the top:
module seed_state_from_config_mod
   use iso_fortran_env, only: real64
   use swap_config_mod, only: swap_config_t
   implicit none
   private

   public :: seed_state_from_config

contains

   subroutine seed_state_from_config(config, state)
      type(swap_config_t), target, intent(inout) :: config
      type(swap_state_t),          intent(inout) :: state

      state%cfg => config
      ! Subsystem inits dispatched from swap_mod (see Task 12 for end-state).
   end subroutine seed_state_from_config

end module seed_state_from_config_mod
```

Also remove any leftover `public ::` entries (the `apply_*` helpers are all gone now).

- [ ] **Step 5: Update the test file similarly**

```bash
git mv tests/unit/io/toml/test_config_to_variables.pf tests/unit/io/toml/test_seed_state_from_config.pf
```

Open the renamed test, update:
- `use config_to_variables_mod, only: config_to_variables` → `use seed_state_from_config_mod, only: seed_state_from_config`
- Any `call config_to_variables(config, state)` → `call seed_state_from_config(config, state)`

- [ ] **Step 6: Update swap_mod**

In `src/core/swap_mod.f90`, change:
- `use config_to_variables_mod, only: config_to_variables` → `use seed_state_from_config_mod, only: seed_state_from_config`
- `call config_to_variables(config, state)` → `call seed_state_from_config(config, state)`

- [ ] **Step 7: Update meson source lists**

In `meson.build`, find `'src/io/toml/config_to_variables.f90'` and rename to `'src/io/toml/seed_state_from_config.f90'`.

In `tests/unit/meson.build`, do the same plus update `test_config_to_variables.F90` references to `test_seed_state_from_config.F90` (the pFUnit-generated file follows the test source name).

```bash
grep -rn "config_to_variables" meson.build tests/unit/meson.build
```

Replace every remaining hit.

- [ ] **Step 8: Clean rebuild and verify**

```bash
rm -rf builddir && pixi run check-fast
```

Expected: 748 pFUnit + 4/4 regression byte-identical. If any file errors due to stale `use config_to_variables_mod` imports:

```bash
grep -rn "config_to_variables" src/ tests/
```

Update every remaining hit.

- [ ] **Step 9: Commit**

```bash
git add -A
git commit -m "$(cat <<'EOF'
refactor(seed): Task 11 — rename config_to_variables → seed_state_from_config

[GR-SEED 2026-05-25] After Tasks 1-10 emptied the adapter into per-
subsystem state%X%init methods, the file is renamed to reflect its
actual current role:

  src/io/toml/config_to_variables.f90    → src/io/toml/seed_state_from_config.f90
  module config_to_variables_mod         → module seed_state_from_config_mod
  subroutine config_to_variables         → subroutine seed_state_from_config
  tests/unit/io/toml/test_*.pf           → test_seed_state_from_config.pf

meson.build entries updated. The single caller in swap_mod.f90 updated.

check-fast 5/5 byte-identical.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 12: Flatten swap_init call chain

**Files:**
- Modify: `src/core/swap_mod.f90` (fold `swap_init_from_loaded_config` into `swap_init`)
- Modify: `src/state/swap_state.f90` if `swap_init_from_loaded_config` is publicly exported

**Spec sections covered:** "Call-chain flattening (Task 12)".

- [ ] **Step 1: Audit public exports**

```bash
grep -n "swap_init_from_loaded_config" src/ tests/ -r
```

If anything outside `swap_mod.f90` itself calls or imports `swap_init_from_loaded_config`, list those — they'll need updates (test fixtures, BMI shim, etc.).

- [ ] **Step 2: Inline the wrapper body into swap_init**

In `src/core/swap_mod.f90`, locate `subroutine swap_init(config_file, state, config)` (around line 19) and `subroutine swap_init_from_loaded_config(state, config)` (around line 53). Copy the **entire body** of `swap_init_from_loaded_config` into `swap_init` immediately after the `errors%abort_if_fatal()` call.

The result inside `swap_init` should be:

```fortran
   subroutine swap_init(config_file, state, config)
      use load_swap_config_mod, only: load_swap_config
      use seed_state_from_config_mod, only: seed_state_from_config
      ! ... other uses migrated from swap_init_from_loaded_config ...
      character(len=*),            intent(in)  :: config_file
      type(swap_state_t),          intent(out) :: state
      type(swap_config_t), target, intent(out) :: config

      ! Phase 1: load + validate + finalize TOML
      block
         use error_mod, only: error_collection_t
         type(error_collection_t) :: errors
         call load_swap_config(config_file, config, errors)
         call config%validate(errors)
         call config%finalize(errors)
         call errors%abort_if_fatal()
      end block

      ! Phase 2: seed state from config (pure config→state, no compute)
      call seed_state_from_config(config, state)
      ! [Task 12] Note: most subsystem state%X%init calls also live here for
      ! now; ideally Task 11's seed_state_from_config orchestrates them.

      ! Migrated from swap_init_from_loaded_config: itertime, timecontrol_init,
      ! CalcGrid, all state%X%init calls not yet inside seed_state_from_config,
      ! SoilWater(1), SwapOutput(1), SoilWaterOutput(1), log_info.
      ! [INSERT ALL LINES FROM CURRENT swap_init_from_loaded_config BODY HERE]
   end subroutine swap_init
```

The exact list of "Phase 3" calls (timecontrol_init, CalcGrid, SoilWater(1), etc.) comes from the current body of `swap_init_from_loaded_config`. Copy them verbatim.

- [ ] **Step 3: Delete `swap_init_from_loaded_config`**

After the body has been copied into `swap_init`, delete the entire `subroutine swap_init_from_loaded_config(state, config) … end subroutine swap_init_from_loaded_config` block from `swap_mod.f90`.

Remove `swap_init_from_loaded_config` from the module's `public ::` list (line 15).

- [ ] **Step 4: Update any external callers**

If Step 1 found external callers (e.g., a BMI shim), update them to call `swap_init(...)` instead. If the BMI shim needs an entry point that takes an already-loaded config (rather than a path), keep a separate `swap_init_from_loaded_config` *as a thin two-call wrapper around `seed_state_from_config + Phase 3 inline body*`, but only if there's a real external caller.

- [ ] **Step 5: Clean rebuild and verify**

```bash
rm -rf builddir && pixi run check-fast
```

Expected: 748 pFUnit + 4/4 regression byte-identical.

- [ ] **Step 6: Commit**

```bash
git add src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
refactor(seed): Task 12 — flatten swap_init; retire swap_init_from_loaded_config

[GR-SEED 2026-05-25] swap_init_from_loaded_config wrapper (~330 lines)
folded into swap_init. The end-state call chain is two-deep:

  swap_main → swap_init(config_file, state, config)
                ├─ load + validate + finalize config
                ├─ seed_state_from_config(config, state)   ← Layer A
                └─ [Layer B: timecontrol_init, CalcGrid, SoilWater(1),
                              SwapOutput(1), SoilWaterOutput(1), log_info]

swap_init_from_loaded_config removed from module public list. Phases
are now explicit in swap_init's body.

check-fast 5/5 byte-identical (748 pFUnit + 4/4 regression).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Post-arc audit

After Task 12 commits, run the following to confirm the arc is complete:

```bash
wc -l src/io/toml/seed_state_from_config.f90       # expect: < 50
wc -l src/core/swap_mod.f90                         # expect: shrunk by ~330L
grep -n "swap_init_from_loaded_config" src/ tests/ -r  # expect: 0 hits
grep -rn "config_to_variables" src/ tests/             # expect: 0 hits except comments
grep -rn "^\s*use variables" src/ tests/               # expect: 0 hits except dormant + comments
rm -rf builddir && pixi run check-fast                 # expect: 748 + 4/4 byte-identical
```

If all four succeed, save an arc-completion memory:

```bash
cat > ~/.claude/projects/-home-zawadzkim-Code-swap/memory/project_seed_state_arc_2026-05-25.md <<'EOF'
---
name: project-seed-state-arc-2026-05-25
description: config_to_variables decomposition COMPLETE — adapter retired into per-subsystem state%X%init methods; swap_init call chain flattened
metadata:
  type: project
---

config_to_variables.f90 decomposition arc — **COMPLETE** 2026-05-25.

12 tasks shipped. Adapter (was 1631 lines) became a ~30-line orchestrator
(`seed_state_from_config`). Every subsystem's seeding lives on its own
typed state%X%init(config%Y, …) method. swap_init_from_loaded_config
retired; call chain is now `swap_main → swap_init → seed_state_from_config`.

Strangler pattern truly complete.
EOF
```

Update `MEMORY.md` with one-line entry pointing to the file.

---

## Self-review checklist (filled in by plan author)

**Spec coverage:**
- ✓ All 12 spec tasks have a corresponding plan task.
- ✓ All 12 main-routine sections in the spec inventory are absorbed (general+simulation+output_switches → Task 1; meteorology → Task 2; drainage → Task 7; soil + bottom_boundary → Task 8; irrigation → Task 5; solute → Task 3; surface_water → Task 6; heat → already done in earlier phases, noted in spec; crop → Task 10).
- ✓ All 5 helper subroutines have relocation targets (populate_outdatint_monthly → Task 1; apply_soil_tillage → Task 9; apply_irrigation_ssdi + mode0/1 → Task 5; apply_nutrients + apply_nutrients_events → Task 4).
- ✓ Open question on bottom_boundary resolved: folded into state%soilwater%init via `seed_bottom_boundary` private helper (Task 8).
- ✓ Open question on crop_config_global resolved: set inside state%crop%init (Task 10).
- ⚠ Open question on `parse_iso_date_to_days1900`: relocated to `toml_field_helpers_mod` in Task 9.

**Type consistency:**
- ✓ Signatures match between tasks (e.g., `state%mesh` arg in Task 5 matches `state%mesh%numnod` arg in Tasks 4/7).
- ✓ `pathwork_in` threading is consistent across SSDI helpers and crop init.

**Placeholder scan:**
- ⚠ Tasks 3, 7, 8, 9 contain `[INSERT here the body from current …]` placeholders. These are deliberate — the implementer has to read the source file to fetch the exact body. The placeholder is bracketed with the exact line range and a clear instruction.
