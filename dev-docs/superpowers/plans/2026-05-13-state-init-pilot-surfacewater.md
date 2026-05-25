# State Init Pilot — Surfacewater Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the `surfacewater_init` procedure (currently buried inside `SurfaceWater(task=1)`) with a type-bound `state%surfacewater%init(config_sw, config_drain, numnod)` procedure on the state module, called explicitly from `swap_init`'s S4 pass.

**Architecture:** Single subsystem pilot for the state-init consolidation arc. The init body absorbs the existing `surfacewater_init.f90` (148 LoC) plus the post-init block in `SurfaceWater(task=1)` (hwlman/vtair zero + ZDraBas init). The body reads typed config slices instead of legacy globals, inlines the `wls1_init = wlact - altcu` computation (retiring one transient-buffer indirection), and is unit-testable in isolation. The old `surfacewater_init.f90` module is deleted; `SurfaceWater(task=1)` case becomes a stub.

**Tech Stack:** Fortran 2008+ free-form, meson + ninja build (`pixi run build-linux`), pFUnit unit tests (`pixi run test-pfunit`), Python regression harness (`pixi run check-fast` for everyday gate, `check-full` for end-of-arc).

**Reference spec:** `docs/superpowers/specs/2026-05-13-state-init-pilot-surfacewater-design.md`. Read before starting.

**Per-task verification gate (mandatory):** Per project workflow convention (memory: `feedback_per_task_regression_gate.md`), every task's "Verify" step ends with **all three**:

1. `pixi run build-linux` — no new warnings
2. `pixi run test-pfunit` — all tests pass
3. `pixi run check-fast` — regression byte-identical (4 cases: hupselbrook, surfacewater, salinitystress, grassgrowth)

For tasks 1-4 the runtime behavior is unchanged (new procedure exists but is uncalled), so `check-fast` must be byte-identical to baseline. Task 5 is where the call switches; same gate.

---

## File Structure

### Files modified

- `src/state/surfacewater_state.f90` — add `init` type-bound procedure with full body (Tasks 1-4).
- `src/utils/surfacewaterutils.f90` — extract `swstlev_from_table(sttab, wlev)` helper from existing `swstlev(state, wlev)` (Task 3).
- `tests/unit/state/test_surfacewater_state.pf` — add new tests covering `init` behavior (Tasks 1-4).
- `src/drainage/surfacewater.f90` — stub `SurfaceWater(task=1)` case(1) to a `return` (Task 5).
- `src/core/swap_mod.f90` — add `if (flSurfaceWater) call state%surfacewater%init(...)` to the S4 pass (Task 5).
- `meson.build` — remove `'src/drainage/surfacewater_init.f90'` from the source list (Task 5).

### Files deleted

- `src/drainage/surfacewater_init.f90` (Task 5).

### Files NOT touched

- `src/io/toml/config_to_variables.f90` — explicitly out of scope. The orphan `wls1_init` write at line 1133 stays for the separate follow-up arc.
- `src/state/variables.f90` (or wherever `wls1_init` global lives) — same reason.
- All other `_init` procedures (`soilwater_init`, `atmosphere_init`, etc.) — separate arc.

---

## Task 1: Add `init` type-bound procedure with allocation logic + first test

**Files:**
- Modify: `src/state/surfacewater_state.f90`
- Modify: `tests/unit/state/test_surfacewater_state.pf`

### Step 1.1: Add the type-bound procedure declaration

- [ ] **Add `init` to the `surfacewater_state_t` `contains` block.**

In `src/state/surfacewater_state.f90`, locate the `contains` block inside `type :: surfacewater_state_t` (after the cohort flattening this currently has three `procedure ::` lines for `reset_intermediate`, `reset_cumulative_drainage`, `reset_cumulative_reservoir`). Add an `init` procedure line at the top:

```fortran
   contains
      procedure :: init                        => surfacewater_state_init
      procedure :: reset_intermediate          => surfacewater_reset_intermediate
      procedure :: reset_cumulative_drainage   => surfacewater_reset_cumulative_drainage
      procedure :: reset_cumulative_reservoir  => surfacewater_reset_cumulative_reservoir
   end type surfacewater_state_t
```

### Step 1.2: Add the `use` imports needed by the init body

- [ ] **At the top of the module (after `use, intrinsic :: iso_fortran_env`), add:**

```fortran
   use surface_water_config_mod, only: surface_water_config_t
   use drainage_config_mod,      only: drainage_config_t
   use surfacewater_utils,       only: swstlev_from_table
   use error_mod,                only: fatalerr_collected
```

(The `swstlev_from_table` helper is extracted in Task 3 — for Task 1 it's not yet used so the import will be unused-warning until Task 3. To avoid the warning during the iterative tasks, defer that single import to Task 3 and add only the other three here.)

Concretely for Task 1, add:

```fortran
   use surface_water_config_mod, only: surface_water_config_t
   use drainage_config_mod,      only: drainage_config_t
   use error_mod,                only: fatalerr_collected
```

### Step 1.3: Add the `surfacewater_state_init` subroutine (allocation only)

- [ ] **Add the procedure body to the module's `contains` section** (alongside `surfacewater_reset_intermediate` etc., before `end module`):

```fortran
   !> One-time runtime initialization for surfacewater state.
   !! Replaces the surviving math from legacy `rddre` (now retired from
   !! readswap) plus the post-init block inside `SurfaceWater(task=1)`.
   !!
   !! Scope: swsrf=2, swsec=2, swqhr=1, swman=1, drainage.altcu=0 only.
   !! Other branches are guarded with fatalerr_collected (defense in
   !! depth — surface_water_config_validate rejects them upstream too).
   subroutine surfacewater_state_init(self, config_sw, config_drain, numnod)
      class(surfacewater_state_t),  intent(inout) :: self
      type(surface_water_config_t), intent(in)    :: config_sw
      type(drainage_config_t),      intent(in)    :: config_drain
      integer,                      intent(in)    :: numnod

      ! Defensive guards mirroring surface_water_config_validate.
      ! swman is a fixed-size array (dimensioned mamp); slice 1:nmper.
      if (config_sw%swsrf == 3 .or. config_sw%swsec == 1 .or. config_sw%swqhr == 2) then
         call fatalerr_collected('surfacewater_state_init', &
            'swsrf=3, swsec=1, or swqhr=2 not supported on the TOML path')
         return
      end if
      if (allocated(config_sw%swman)) then
         if (any(config_sw%swman(1:config_sw%nmper) == 2)) then
            call fatalerr_collected('surfacewater_state_init', &
               'swman=2 (automatic weir) not supported on the TOML path')
            return
         end if
      end if

      ! Allocate per-level arrays.
      allocate(self%cqdrain    (config_drain%nrlevs));            self%cqdrain    = 0.0_real64
      allocate(self%cqdrainin  (config_drain%nrlevs));            self%cqdrainin  = 0.0_real64
      allocate(self%cqdrainout (config_drain%nrlevs));            self%cqdrainout = 0.0_real64
      allocate(self%inqdra     (config_drain%nrlevs, numnod));    self%inqdra     = 0.0_real64
      allocate(self%inqdra_in  (config_drain%nrlevs, numnod));    self%inqdra_in  = 0.0_real64
      allocate(self%inqdra_out (config_drain%nrlevs, numnod));    self%inqdra_out = 0.0_real64

   end subroutine surfacewater_state_init
```

### Step 1.4: Add the first pFUnit test

- [ ] **Add this test to `tests/unit/state/test_surfacewater_state.pf` (append at the end of the file, after the existing tests):**

```fortran
@test
subroutine test_surfacewater_init_allocates_cohort_arrays()
   use funit
   use iso_fortran_env, only: real64
   use surfacewater_state_mod, only: surfacewater_state_t
   use surface_water_config_mod, only: surface_water_config_t
   use drainage_config_mod, only: drainage_config_t
   type(surfacewater_state_t)   :: sw
   type(surface_water_config_t) :: cfg_sw
   type(drainage_config_t)      :: cfg_drain
   integer, parameter :: nrlevs = 3
   integer, parameter :: numnod = 50

   ! Minimal valid config for swsrf=2, swsec=2, swqhr=1, no automatic weir.
   cfg_sw%swsrf = 2
   cfg_sw%swsec = 2
   cfg_sw%swqhr = 1
   cfg_sw%nmper = 1
   allocate(cfg_sw%swman(1));  cfg_sw%swman = 1
   cfg_sw%wlact = -50.0_real64

   cfg_drain%nrlevs = nrlevs
   cfg_drain%altcu  = 0.0_real64
   allocate(cfg_drain%zbotdr(nrlevs));   cfg_drain%zbotdr = [-100.0_real64, -120.0_real64, -150.0_real64]
   allocate(cfg_drain%swdtyp(nrlevs));   cfg_drain%swdtyp = 0     ! open channels
   allocate(cfg_drain%widthr(nrlevs));   cfg_drain%widthr = 1.0_real64
   allocate(cfg_drain%taludr(nrlevs));   cfg_drain%taludr = 2.0_real64
   allocate(cfg_drain%l(nrlevs));        cfg_drain%l      = 1000.0_real64

   call sw%init(cfg_sw, cfg_drain, numnod)

   @assertTrue(allocated(sw%cqdrain))
   @assertTrue(allocated(sw%cqdrainin))
   @assertTrue(allocated(sw%cqdrainout))
   @assertTrue(allocated(sw%inqdra))
   @assertTrue(allocated(sw%inqdra_in))
   @assertTrue(allocated(sw%inqdra_out))
   @assertEqual(nrlevs, size(sw%cqdrain))
   @assertEqual([nrlevs, numnod], shape(sw%inqdra))
   @assertEqual([nrlevs, numnod], shape(sw%inqdra_in))
   @assertEqual([nrlevs, numnod], shape(sw%inqdra_out))
end subroutine test_surfacewater_init_allocates_cohort_arrays
```

### Step 1.5: Verify

- [ ] **Run verification gate:**

```bash
cd /home/zawadzkim/Code/swap
pixi run build-linux
pixi run test-pfunit
pixi run check-fast
```

Expected:
- Build: clean.
- pFUnit: all tests pass (including the new `test_surfacewater_init_allocates_cohort_arrays`). Test count increases by 1 from baseline 733 → 734.
- check-fast: 4/4 regression cases byte-identical (the new procedure exists but is uncalled, so runtime behavior is unchanged).

### Step 1.6: Commit

- [ ] **Commit:**

```bash
cd /home/zawadzkim/Code/swap
git add src/state/surfacewater_state.f90 tests/unit/state/test_surfacewater_state.pf

git commit -m "$(cat <<'EOF'
refactor(state): add surfacewater_state_t%init stub with allocation

First task of the state-init pilot arc (spec
docs/superpowers/specs/2026-05-13-state-init-pilot-surfacewater-design.md).
Adds a type-bound `init(config_sw, config_drain, numnod)` procedure on
surfacewater_state_t. Body currently does allocation only (cohort
per-level arrays). New pFUnit test covers the allocation behavior.
Procedure is not yet called from anywhere; coexists with the existing
surfacewater_init in src/drainage/surfacewater_init.f90.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Add zero defaults + config seed (L1 + L2)

**Files:**
- Modify: `src/state/surfacewater_state.f90`
- Modify: `tests/unit/state/test_surfacewater_state.pf`

### Step 2.1: Extend the init body with zero defaults and config seed

- [ ] **In `src/state/surfacewater_state.f90`, inside `surfacewater_state_init`, insert after the defensive-guard block (before the allocation block) the following lines:**

```fortran
      ! ---- L1: zero defaults ----
      self%numadj = 0
      self%wlsbak = 0.0_real64

      ! ---- L2: config-derived seeds ----
      ! wls1 = wlact - altcu (legacy rddre line; altcu=0 enforced by drainage_config_validate
      ! so this equals wlact). Inlined here — retires the wls1_init transient buffer
      ! indirection at the call site.
      self%wls    = config_sw%wlact - config_drain%altcu
      self%wlstar = self%wls
```

The order is: defensive guards → L1 zero defaults → L2 config seed → L0 allocation. (Final body order will be cleaned up in Task 4; for now keep edits additive.)

### Step 2.2: Note on `wlp` global write

The legacy `surfacewater_init` body also wrote `wlp = 0.0d0` (a legacy global). The spec keeps that legacy write retained (other physics still reads `wlp` via `use variables`). Add it to the init body:

- [ ] **Add at the top of the init body, after the `intent(in)` declarations and after the defensive-guard block, before L1:**

```fortran
      ! Legacy global write retained for now: bocodre reads wlp via `use variables`
      ! for the primary surface water level (not surfacewater-state owned). A separate
      ! arc will migrate readers; this write stays until then.
      block
         use variables, only: wlp
         wlp = 0.0_real64  ! swsrf=2 has no primary system
      end block
```

(Using an internal `block` keeps the legacy-global dependency scoped to the one line that needs it.)

### Step 2.3: Add tests for the zero defaults and config seed

- [ ] **Append to `tests/unit/state/test_surfacewater_state.pf`:**

```fortran
@test
subroutine test_surfacewater_init_sets_wls_from_config()
   use funit
   use iso_fortran_env, only: real64
   use surfacewater_state_mod, only: surfacewater_state_t
   use surface_water_config_mod, only: surface_water_config_t
   use drainage_config_mod, only: drainage_config_t
   type(surfacewater_state_t)   :: sw
   type(surface_water_config_t) :: cfg_sw
   type(drainage_config_t)      :: cfg_drain
   integer, parameter :: nrlevs = 1
   integer, parameter :: numnod = 10

   cfg_sw%swsrf = 2;  cfg_sw%swsec = 2;  cfg_sw%swqhr = 1
   cfg_sw%nmper = 1
   allocate(cfg_sw%swman(1));  cfg_sw%swman = 1
   cfg_sw%wlact = -75.0_real64

   cfg_drain%nrlevs = nrlevs
   cfg_drain%altcu  = 0.0_real64
   allocate(cfg_drain%zbotdr(nrlevs));  cfg_drain%zbotdr = -100.0_real64
   allocate(cfg_drain%swdtyp(nrlevs));  cfg_drain%swdtyp = 0
   allocate(cfg_drain%widthr(nrlevs));  cfg_drain%widthr = 1.0_real64
   allocate(cfg_drain%taludr(nrlevs));  cfg_drain%taludr = 2.0_real64
   allocate(cfg_drain%l(nrlevs));       cfg_drain%l      = 1000.0_real64

   call sw%init(cfg_sw, cfg_drain, numnod)

   @assertEqual(-75.0_real64, sw%wls,    1.0e-12_real64)
   @assertEqual(-75.0_real64, sw%wlstar, 1.0e-12_real64)
   @assertEqual(0, sw%numadj)
   @assertEqual(0.0_real64, sum(sw%wlsbak), 1.0e-12_real64)
end subroutine test_surfacewater_init_sets_wls_from_config
```

### Step 2.4: Verify

- [ ] **Run verification gate:**

```bash
cd /home/zawadzkim/Code/swap
pixi run build-linux
pixi run test-pfunit
pixi run check-fast
```

Expected:
- Build: clean.
- pFUnit: all pass. Test count 734 → 735.
- check-fast: byte-identical (new code still uncalled).

### Step 2.5: Commit

- [ ] **Commit:**

```bash
cd /home/zawadzkim/Code/swap
git add src/state/surfacewater_state.f90 tests/unit/state/test_surfacewater_state.pf

git commit -m "$(cat <<'EOF'
refactor(state): add zero defaults and config seeds to surfacewater init

Extends surfacewater_state_t%init with L1 (zero defaults: numadj,
wlsbak) and L2 (config-derived seeds: wls, wlstar, plus legacy wlp
write retained pending separate reader-migration arc). The wls1
formula (wlact - altcu) is inlined directly from typed config —
retiring the wls1_init transient buffer at the call site once the
hoist lands in Task 5.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Extract `swstlev_from_table` helper, add sttab build + swst computation

The existing `swstlev(state, wlev)` in `src/utils/surfacewaterutils.f90` takes a full `swap_state_t`, which doesn't match a state-bound init (which has access only to `self`, the `surfacewater_state_t`). Extract a pure-math helper `swstlev_from_table(sttab, wlev)` and make `swstlev(state, wlev)` a thin wrapper. Then use the helper in init.

**Files:**
- Modify: `src/utils/surfacewaterutils.f90`
- Modify: `src/state/surfacewater_state.f90`
- Modify: `tests/unit/state/test_surfacewater_state.pf`

### Step 3.1: Extract `swstlev_from_table` helper in surfacewater_utils

- [ ] **In `src/utils/surfacewaterutils.f90`, update the `public ::` line to export the new helper:**

```fortran
   public :: wlevst, swstlev, swstlev_from_table, qhtab, runoff
```

- [ ] **Add the new helper at the bottom of the `contains` section (before `end module`):**

```fortran
   !> Pure-table-lookup variant of `swstlev` — operates directly on the storage
   !! table without needing a full `swap_state_t`. Used by `surfacewater_state_t%init`
   !! where only the partially-populated state is available.
   function swstlev_from_table(sttab, wlev) result(swstlev_r)
      implicit none
      real(real64), intent(in) :: sttab(22, 2)
      real(real64), intent(in) :: wlev
      real(real64) :: swstlev_r

      integer :: i
      real(real64) :: dwl
      character(len=200) :: messag

      if (wlev < sttab(22,1)) then
         messag = 'Surface water storage below bottom of table'
         call fatalerr_collected('swstlev_from_table', messag)
      end if
      if (wlev > sttab(1,1)) then
         messag = 'Surface water storage above top of table'
         call fatalerr_collected('swstlev_from_table', messag)
      end if

      i = 0
      do
         i = i + 1
         if (wlev >= sttab(i+1,1) .and. wlev <= sttab(i,1)) exit
      end do

      dwl = (wlev - sttab(i+1,1)) / (sttab(i,1) - sttab(i+1,1))
      swstlev_r = sttab(i+1,2) + dwl * (sttab(i,2) - sttab(i+1,2))
   end function swstlev_from_table
```

- [ ] **Replace the body of the existing `swstlev(state, wlev)` to delegate to the new helper:**

Find the function (around line 102-138) and replace its body so the function becomes a thin wrapper:

```fortran
   function swstlev(state, wlev) result(swstlev_r)
      implicit none
      type(swap_state_t), intent(in) :: state
      real(real64),       intent(in) :: wlev
      real(real64) :: swstlev_r

      swstlev_r = swstlev_from_table(state%surfacewater%sttab, wlev)
   end function swstlev
```

### Step 3.2: Add the import in surfacewater_state.f90

- [ ] **In `src/state/surfacewater_state.f90`, add to the use-list:**

```fortran
   use surfacewater_utils, only: swstlev_from_table
```

### Step 3.3: Add sttab build + swst computation to the init body

- [ ] **In `surfacewater_state_init`, after the allocation block, append:**

```fortran
      ! ---- L3: readswap-style shape math (sttab + swst init) ----
      !
      ! sttab(:,1): water-level rows.
      !   Row 1 = +100cm above soil surface (top).
      !   Row 2 = 0cm (soil surface).
      !   Rows 3..22: divide [0, zbotdr(1+nrpri)] into 20 compartments.
      ! For swsrf=2 (no primary system) nrpri=0, so zbotdr index is 1.
      block
         integer :: i, ilev, nrpri
         real(real64) :: wdepth, wvolum, wbreadth

         nrpri = 0

         self%sttab(1, 1) = 100.0_real64
         self%sttab(2, 1) =   0.0_real64
         do i = 3, 22
            self%sttab(i, 1) = config_drain%zbotdr(1 + nrpri) * real(i - 2, real64) / 20.0_real64
         end do

         ! sttab(:,2): storage volume per unit area (cm), summed across
         ! open-channel levels (swdtyp=0). Verbatim port from legacy rddre
         ! (readswap.f90:4878-4897). l(:) is in centimetres (D6 conversion
         ! done at TOML read time).
         do i = 1, 22
            self%sttab(i, 2) = 0.0_real64
            do ilev = 1 + nrpri, config_drain%nrlevs
               if (config_drain%swdtyp(ilev) == 0 .and. self%sttab(i, 1) > config_drain%zbotdr(ilev)) then
                  if (self%sttab(i, 1) <= 0.0_real64) then
                     ! Trapezium below soil surface
                     wdepth = self%sttab(i, 1) - config_drain%zbotdr(ilev)
                     wvolum = wdepth * (config_drain%widthr(ilev) + wdepth / config_drain%taludr(ilev))
                  else
                     ! Trapezium up to surface, plus rectangle above
                     wdepth   = -config_drain%zbotdr(ilev)
                     wvolum   = wdepth * (config_drain%widthr(ilev) + wdepth / config_drain%taludr(ilev))
                     wbreadth = config_drain%widthr(ilev) + 2.0_real64 * wdepth / config_drain%taludr(ilev)
                     wdepth   = self%sttab(i, 1)
                     wvolum   = wvolum + wbreadth * wdepth
                  end if
                  self%sttab(i, 2) = self%sttab(i, 2) + wvolum / config_drain%l(ilev)
               end if
            end do
         end do
      end block

      ! Initial storage state derived from sttab + wls.
      self%swstini = swstlev_from_table(self%sttab, self%wls)
      self%swst    = self%swstini
```

### Step 3.4: Add tests for sttab + swst

- [ ] **Append to `tests/unit/state/test_surfacewater_state.pf`:**

```fortran
@test
subroutine test_surfacewater_init_builds_sttab()
   use funit
   use iso_fortran_env, only: real64
   use surfacewater_state_mod, only: surfacewater_state_t
   use surface_water_config_mod, only: surface_water_config_t
   use drainage_config_mod, only: drainage_config_t
   type(surfacewater_state_t)   :: sw
   type(surface_water_config_t) :: cfg_sw
   type(drainage_config_t)      :: cfg_drain
   integer, parameter :: nrlevs = 1
   integer, parameter :: numnod = 10
   real(real64), parameter :: zbot = -100.0_real64
   integer :: i

   cfg_sw%swsrf = 2;  cfg_sw%swsec = 2;  cfg_sw%swqhr = 1
   cfg_sw%nmper = 1
   allocate(cfg_sw%swman(1));  cfg_sw%swman = 1
   cfg_sw%wlact = -50.0_real64

   cfg_drain%nrlevs = nrlevs
   cfg_drain%altcu  = 0.0_real64
   allocate(cfg_drain%zbotdr(nrlevs));  cfg_drain%zbotdr = zbot
   allocate(cfg_drain%swdtyp(nrlevs));  cfg_drain%swdtyp = 0
   allocate(cfg_drain%widthr(nrlevs));  cfg_drain%widthr = 1.0_real64
   allocate(cfg_drain%taludr(nrlevs));  cfg_drain%taludr = 2.0_real64
   allocate(cfg_drain%l(nrlevs));       cfg_drain%l      = 1000.0_real64

   call sw%init(cfg_sw, cfg_drain, numnod)

   ! Row 1 = +100, row 2 = 0, rows 3..22 divide [0, zbot] = [0, -100] into 20 compartments.
   @assertEqual( 100.0_real64, sw%sttab(1, 1), 1.0e-12_real64)
   @assertEqual(   0.0_real64, sw%sttab(2, 1), 1.0e-12_real64)
   @assertEqual(  -5.0_real64, sw%sttab(3, 1), 1.0e-12_real64)   ! zbot * (3-2)/20 = -100*1/20 = -5
   @assertEqual(-100.0_real64, sw%sttab(22, 1), 1.0e-12_real64)  ! zbot * (22-2)/20 = -100

   ! Storage column must be monotonically decreasing as depth increases
   ! (more depth = more volume from the table-row perspective).
   do i = 1, 21
      @assertTrue(sw%sttab(i, 2) >= sw%sttab(i+1, 2))
   end do
end subroutine test_surfacewater_init_builds_sttab

@test
subroutine test_surfacewater_init_computes_swst_from_sttab()
   use funit
   use iso_fortran_env, only: real64
   use surfacewater_state_mod, only: surfacewater_state_t
   use surface_water_config_mod, only: surface_water_config_t
   use drainage_config_mod, only: drainage_config_t
   use surfacewater_utils, only: swstlev_from_table
   type(surfacewater_state_t)   :: sw
   type(surface_water_config_t) :: cfg_sw
   type(drainage_config_t)      :: cfg_drain
   integer, parameter :: nrlevs = 1
   integer, parameter :: numnod = 10

   cfg_sw%swsrf = 2;  cfg_sw%swsec = 2;  cfg_sw%swqhr = 1
   cfg_sw%nmper = 1
   allocate(cfg_sw%swman(1));  cfg_sw%swman = 1
   cfg_sw%wlact = -50.0_real64

   cfg_drain%nrlevs = nrlevs
   cfg_drain%altcu  = 0.0_real64
   allocate(cfg_drain%zbotdr(nrlevs));  cfg_drain%zbotdr = -100.0_real64
   allocate(cfg_drain%swdtyp(nrlevs));  cfg_drain%swdtyp = 0
   allocate(cfg_drain%widthr(nrlevs));  cfg_drain%widthr = 1.0_real64
   allocate(cfg_drain%taludr(nrlevs));  cfg_drain%taludr = 2.0_real64
   allocate(cfg_drain%l(nrlevs));       cfg_drain%l      = 1000.0_real64

   call sw%init(cfg_sw, cfg_drain, numnod)

   ! swst and swstini must equal the table-lookup at wls=-50.
   @assertEqual(swstlev_from_table(sw%sttab, sw%wls), sw%swst,    1.0e-12_real64)
   @assertEqual(sw%swst, sw%swstini, 1.0e-12_real64)
end subroutine test_surfacewater_init_computes_swst_from_sttab
```

### Step 3.5: Verify

- [ ] **Run verification gate:**

```bash
cd /home/zawadzkim/Code/swap
pixi run build-linux
pixi run test-pfunit
pixi run check-fast
```

Expected:
- Build: clean.
- pFUnit: all pass. Test count 735 → 737.
- check-fast: byte-identical. The `swstlev_from_table` extraction is a pure refactor (the existing `swstlev` body delegates with identical semantics), and the init procedure is still uncalled.

### Step 3.6: Commit

- [ ] **Commit:**

```bash
cd /home/zawadzkim/Code/swap
git add src/state/surfacewater_state.f90 \
        src/utils/surfacewaterutils.f90 \
        tests/unit/state/test_surfacewater_state.pf

git commit -m "$(cat <<'EOF'
refactor(state): add sttab build + swst init to surfacewater init

Extracts a pure-math swstlev_from_table(sttab, wlev) helper from the
existing swstlev(state, wlev), so the state-bound init (which has only
self : surfacewater_state_t, not a full swap_state_t) can compute
initial storage. swstlev becomes a thin wrapper for backwards
compatibility. surfacewater_state_t%init now builds the open-channel
storage table sttab from drainage geometry and computes swst/swstini
from sttab + wls. Two new pFUnit tests cover sttab shape and swst
correctness.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Add post-init block (hwlman/vtair zero + ZDraBas)

This is the block currently inside `SurfaceWater(task=1)` at `src/drainage/surfacewater.f90:68-93`. After the hoist in Task 5, the call site will be gone, so this work must be absorbed into init now.

**Files:**
- Modify: `src/state/surfacewater_state.f90`
- Modify: `tests/unit/state/test_surfacewater_state.pf`

### Step 4.1: Add hwlman / vtair zero defaults + ZDraBas block

- [ ] **In `surfacewater_state_init`, after the L3 block from Task 3, append:**

```fortran
      ! ---- Post-init shape math (was in SurfaceWater(task=1) post-call block) ----
      ! hwlman/vtair: only output reads them; default to zero here.
      self%hwlman = 0.0_real64
      self%vtair  = 0.0_real64

      ! ZDraBas: macropore drainage basis. Only initialized once (flInitDraBas guards
      ! subsequent runs but in the pilot init always runs once at startup so we
      ! unconditionally set it here and clear the flag).
      !
      ! NOTE: NumLevRapDra and wlstab/maowl/tc_t1900 are legacy globals; macropore
      ! is retired (ADR 0040), so NumLevRapDra is permanently 0 and the swdtyp(0)
      ! branch is unreachable. For swsec=2 (the only branch we support on the TOML
      ! path; the validator enforces this), ZDraBas := wlstar — which we already
      ! seeded in the L2 step.
      if (config_sw%swsec == 2) then
         self%ZDraBas = self%wlstar
      end if
      ! swsec=1 branch (legacy wlstab table lookup) and swdtyp(NumLevRapDra)=1 branch
      ! (drain tube — needs zbotdr(NumLevRapDra) where NumLevRapDra=0 under
      ! macropore-retire) are both UNREACHABLE on the TOML path. Validator rejects
      ! swsec=1; macropore retirement zeros NumLevRapDra. Defensive guard added:
      if (config_sw%swsec /= 2) then
         call fatalerr_collected('surfacewater_state_init', &
            'swsec /= 2 not supported on the TOML path (validator rejects)')
         return
      end if

      self%flInitDraBas = .false.
```

### Step 4.2: Add test for ZDraBas

- [ ] **Append to `tests/unit/state/test_surfacewater_state.pf`:**

```fortran
@test
subroutine test_surfacewater_init_sets_zdrabas_for_swsec_2()
   use funit
   use iso_fortran_env, only: real64
   use surfacewater_state_mod, only: surfacewater_state_t
   use surface_water_config_mod, only: surface_water_config_t
   use drainage_config_mod, only: drainage_config_t
   type(surfacewater_state_t)   :: sw
   type(surface_water_config_t) :: cfg_sw
   type(drainage_config_t)      :: cfg_drain
   integer, parameter :: nrlevs = 1
   integer, parameter :: numnod = 10

   cfg_sw%swsrf = 2;  cfg_sw%swsec = 2;  cfg_sw%swqhr = 1
   cfg_sw%nmper = 1
   allocate(cfg_sw%swman(1));  cfg_sw%swman = 1
   cfg_sw%wlact = -25.0_real64

   cfg_drain%nrlevs = nrlevs
   cfg_drain%altcu  = 0.0_real64
   allocate(cfg_drain%zbotdr(nrlevs));  cfg_drain%zbotdr = -100.0_real64
   allocate(cfg_drain%swdtyp(nrlevs));  cfg_drain%swdtyp = 0
   allocate(cfg_drain%widthr(nrlevs));  cfg_drain%widthr = 1.0_real64
   allocate(cfg_drain%taludr(nrlevs));  cfg_drain%taludr = 2.0_real64
   allocate(cfg_drain%l(nrlevs));       cfg_drain%l      = 1000.0_real64

   call sw%init(cfg_sw, cfg_drain, numnod)

   ! swsec=2 → ZDraBas = wlstar = wls = -25.
   @assertEqual(-25.0_real64, sw%ZDraBas, 1.0e-12_real64)
   @assertEqual(0.0_real64, sw%hwlman, 1.0e-12_real64)
   @assertEqual(0.0_real64, sw%vtair,  1.0e-12_real64)
   @assertFalse(sw%flInitDraBas)
end subroutine test_surfacewater_init_sets_zdrabas_for_swsec_2
```

### Step 4.3: Verify

- [ ] **Run verification gate:**

```bash
cd /home/zawadzkim/Code/swap
pixi run build-linux
pixi run test-pfunit
pixi run check-fast
```

Expected:
- Build: clean.
- pFUnit: all pass. Test count 737 → 738.
- check-fast: byte-identical.

### Step 4.4: Commit

- [ ] **Commit:**

```bash
cd /home/zawadzkim/Code/swap
git add src/state/surfacewater_state.f90 tests/unit/state/test_surfacewater_state.pf

git commit -m "$(cat <<'EOF'
refactor(state): add post-init block (hwlman/vtair/ZDraBas) to surfacewater init

Absorbs the work currently in SurfaceWater(task=1)'s post-call block:
hwlman and vtair zero defaults, and ZDraBas initialization. For the
TOML pipeline scope (swsrf=2, swsec=2; validator rejects others;
macropore retired per ADR 0040), ZDraBas reduces to wlstar. Defensive
guard added against swsec/=2.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Hoist caller, stub SurfaceWater(task=1) case(1), delete surfacewater_init.f90, final check-full

This is the only task that changes runtime behavior. After it, `check-fast` must still be byte-identical (state is initialized to the same values via the new path).

**Files:**
- Modify: `src/core/swap_mod.f90`
- Modify: `src/drainage/surfacewater.f90`
- Modify: `meson.build`
- Delete: `src/drainage/surfacewater_init.f90`

### Step 5.1: Add the call in swap_init

- [ ] **In `src/core/swap_mod.f90`, locate the line:**

```fortran
   if (flSurfaceWater) call SurfaceWater(1, state, request_smaller_dt)
```

(currently around line 161 — find with `grep -n "SurfaceWater(1" src/core/swap_mod.f90`)

- [ ] **Replace it with the hoisted init call followed by the (soon-to-be-stub) SurfaceWater(1) call:**

```fortran
   ! State init S4: hoisted out of SurfaceWater(task=1) per the 2026-05-13
   ! state-init pilot (spec docs/superpowers/specs/2026-05-13-state-init-pilot-surfacewater-design.md).
   if (flSurfaceWater) call state%surfacewater%init(config%surface_water, config%drain, numnod)
   if (flSurfaceWater) call SurfaceWater(1, state, request_smaller_dt)   ! case(1) is now a stub
```

`numnod` is already in the use-list (line ~24 of swap_mod.f90). `config%surface_water` and `config%drain` are accessible since `config` is the type-`swap_config_t` argument to `swap_init`.

### Step 5.2: Stub SurfaceWater(task=1) case(1)

- [ ] **In `src/drainage/surfacewater.f90`, locate the `case (1)` block (around line 61-95).**

- [ ] **Replace the entire `case (1)` block body with:**

```fortran
      case (1)

! === initialization — HOISTED ==========================================
! State init for surfacewater was hoisted to swap_init's S4 pass on 2026-05-13
! (spec docs/superpowers/specs/2026-05-13-state-init-pilot-surfacewater-design.md).
! state%surfacewater%init(...) runs directly from swap_init before this Phase-1 task.
! case(1) is kept as a no-op stub to preserve the task=1 dispatcher signature;
! removal of the dispatcher case is a separate cleanup follow-up.

      return
```

(Remove the `call surfacewater_init(state)` line, the hwlman/vtair zero lines, the entire `if (state%surfacewater%flInitDraBas)` block, and the trailing `return`. The block becomes just the comment + `return`.)

- [ ] **Remove the now-unused `use surfacewater_init_mod, only: surfacewater_init` line** from `src/drainage/surfacewater.f90` (around line 38).

### Step 5.3: Delete `src/drainage/surfacewater_init.f90`

- [ ] **Delete the file:**

```bash
cd /home/zawadzkim/Code/swap
rm src/drainage/surfacewater_init.f90
```

### Step 5.4: Update meson.build

- [ ] **Remove the line `'src/drainage/surfacewater_init.f90',` from `meson.build`** (around line 156). Use:

```bash
cd /home/zawadzkim/Code/swap
grep -n "surfacewater_init.f90" meson.build
# Open meson.build and delete the matching line, e.g.:
sed -i "/'src\/drainage\/surfacewater_init.f90',/d" meson.build
grep -n "surfacewater_init" meson.build || echo "meson.build clean"
```

### Step 5.5: Audit grep

- [ ] **Confirm no remaining references to the deleted module:**

```bash
cd /home/zawadzkim/Code/swap
grep -rn 'surfacewater_init_mod\|surfacewater_init(' src/ tests/ || echo "audit clean"
```

Expected: no output (the module name and the procedure call have been retired from production code).

### Step 5.6: Verify

- [ ] **Run verification gate (this is the critical step — runtime behavior should be identical):**

```bash
cd /home/zawadzkim/Code/swap
pixi run build-linux
pixi run test-pfunit
pixi run check-fast
```

Expected:
- Build: clean. No "undefined reference" or "missing module" errors.
- pFUnit: all 738 tests pass.
- check-fast: 4/4 regression cases **byte-identical** to baseline. Any drift here is a real bug — STOP and investigate before committing. The most likely drift sources are (a) field-order or call-order divergence between old surfacewater_init and new init body, (b) a missed legacy-global write (e.g., `wlp`), (c) sttab indexing or interpolation mismatch.

### Step 5.7: Run the full check (end-of-arc gate)

- [ ] **Run check-full to cover all 6 regression cases (including the slower 2 omitted by check-fast):**

```bash
cd /home/zawadzkim/Code/swap
pixi run check-full
```

Expected: all 6 cases byte-identical.

### Step 5.8: Commit

- [ ] **Commit:**

```bash
cd /home/zawadzkim/Code/swap
git add src/core/swap_mod.f90 \
        src/drainage/surfacewater.f90 \
        meson.build
git rm src/drainage/surfacewater_init.f90

git commit -m "$(cat <<'EOF'
refactor(state): hoist surfacewater init to swap_init S4 pass

Completes the state-init pilot for surfacewater (spec
docs/superpowers/specs/2026-05-13-state-init-pilot-surfacewater-design.md).

Changes:
- swap_init: adds explicit `if (flSurfaceWater) call state%surfacewater%init(
  config%surface_water, config%drain, numnod)` to the S4 pass, immediately
  before the existing Phase-1 SurfaceWater(1, ...) call.
- SurfaceWater(task=1) case(1): collapses to a no-op stub. The old
  `call surfacewater_init(state)` plus the post-call hwlman/vtair/ZDraBas
  block are now subsumed by state%surfacewater%init.
- src/drainage/surfacewater_init.f90: deleted (148 LoC retired). meson.build
  source list updated.
- src/drainage/surfacewater.f90: drops the `use surfacewater_init_mod` import.

`check-full` byte-identical across all 6 regression cases. Establishes
the top-level-orchestration + state-bound-init pattern for the remaining
6 subsystems' future migrations.

The wls1_init transient buffer at config_to_variables.f90:1133 is now
orphaned (no reader); its removal awaits the separate adapter-retirement
arc.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Plan Self-Review

(Run after writing the complete plan; fix issues inline before handoff.)

**Spec coverage check:**

- Spec §"Decision: target architecture" / state-bound `init` pattern → Tasks 1-4 build the procedure incrementally.
- Spec §"Pilot scope — In scope":
  - Define `surfacewater_state_t%init` → Tasks 1-4 ✓
  - Body absorbs surfacewater_init.f90 + post-init block → Tasks 1-4 collectively ✓
  - Reads typed config slices, not legacy globals → Task 1 onwards ✓ (with `wlp` legacy write exception called out in Task 2)
  - Inline `wlact - altcu` (retire transient buffer) → Task 2 ✓
  - Hoist call site → Task 5 ✓
  - Retire `wlp = 0` legacy write decision → Task 2 retains it inside a local `block` with explanatory comment ✓
  - Delete `src/drainage/surfacewater_init.f90` → Task 5 ✓
  - Stub SurfaceWater(task=1) case(1) → Task 5 ✓
- Spec §"Architecture details — Module dependencies" → Tasks 1-3 add the imports incrementally ✓
- Spec §"Validation defenses" (the two `fatalerr_collected` guards) → Task 1's body ✓
- Spec §"Tests" — five listed tests:
  - `test_surfacewater_init_allocates_cohort_arrays` → Task 1 ✓
  - `test_surfacewater_init_sets_wls_from_config` → Task 2 ✓
  - `test_surfacewater_init_builds_sttab` → Task 3 ✓
  - `test_surfacewater_init_computes_swst_from_sttab` → Task 3 ✓
  - `test_surfacewater_init_sets_zdrabas_for_swsec_2` → Task 4 ✓
- Spec §"Verification gates" — build + pFUnit + check-fast per task; check-full at end → Tasks 1-5 each have all three; Task 5 also runs check-full ✓

No spec section uncovered.

**Placeholder scan:** no "TBD", "implement later", or "similar to". All code blocks are concrete.

**Type consistency:**
- `surfacewater_state_init` (implementation) ↔ `init` (type-bound name) ↔ `state%surfacewater%init(...)` (caller) — consistent across all tasks. ✓
- `swstlev_from_table(sttab, wlev)` — signature consistent in surfacewater_utils (Task 3.1), use in init body (Task 3.3), use in test (Task 3.4). ✓
- `config_sw`, `config_drain`, `numnod` argument names — consistent in init signature, call site, and all tests. ✓
- `config_sw%swsrf`, `config_sw%swsec`, `config_sw%swqhr`, `config_sw%nmper`, `config_sw%swman`, `config_sw%wlact` — used in init body and tests; field names match `surface_water_config_t`. ✓
- `config_drain%nrlevs`, `config_drain%altcu`, `config_drain%zbotdr`, `config_drain%swdtyp`, `config_drain%widthr`, `config_drain%taludr`, `config_drain%l` — used in init body and tests; field names match `drainage_config_t`. ✓

**Potentially questionable assumption** I'm flagging here: the `block` syntax with `use variables, only: wlp` inside the init body (Task 2 Step 2.2) is rare in this codebase. If the implementer subagent finds a compilation issue with that pattern, the alternative is to import `wlp` at module scope and add a comment explaining the temporary necessity. Either works; the block form keeps the scope tight but the module-scope form is more conventional Fortran.

**One latent concern about Task 4's defensive guard:**

The added `if (config_sw%swsec /= 2)` guard rejects any case where swsec/=2. This is consistent with the spec (the existing surfacewater_init already had the same guard implicit in its else-branch structure — the only branch that DOES anything is swsec=2 under macropore-retired conditions). The validator enforces this upstream. But if a downstream test case ever exercises swsec=1 deliberately (e.g., a future feature unlock), this guard becomes an active rejection. Acceptable for the pilot — the validator's current behavior is to reject these branches, so the defensive fatalerr is correct as long as the validator's policy holds.
