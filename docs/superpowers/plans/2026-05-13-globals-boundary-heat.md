# GR-BH Implementation Plan — Boundary + Heat Globals Retirement + Mesh Extraction

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Retire `use variables` from `boundbottom.f90`, `boundtop.f90`, `frozencond.f90`, `temperature.f90` (4 target files) by extracting a `state%mesh` subrecord, extending `state%soilwater` with 12 fields, extending `state%drainage` with 10 fields, sweeping every other reader across the codebase (~30 files), and deleting the migrated globals from `variables.f90`. Relocate `heat_init` to a type-bound `state%heat%init(config%heat, numnod)` per the surfacewater_state pilot pattern.

**Architecture:** Three sub-phases inside one arc. Phase A is additive schema + dual-write (no reader changes). Phase B cuts the 4 target readers over to state/config and retires `heat_init` to type-bound form. Phase C sweeps every remaining mesh / layer-flat / drainage-geo reader and deletes the bare globals from `variables.f90` (compile-time surfaces any missed reader). Each task ends with the full regression gate; each phase ends with `check-full` 5/5 byte-for-byte.

**Tech Stack:** Fortran 2018, gfortran, Meson, pFUnit 4.15, pixi, Python regression tests.

**Spec:** `docs/superpowers/specs/2026-05-13-globals-boundary-heat-design.md`

**Builds on:** GR-UTILS (Arc 1, completed 2026-05-12 — commit `26123ad`).

**Verification gate per task (VG)** — per memory `feedback_state_schema_clean_rebuild.md` + `feedback_per_task_regression_gate.md`:

```bash
rm -rf builddir && pixi run build-linux \
  && pixi run test-pfunit \
  && pixi run -e test python tests/regression/test_output_regression.py \
       hupselbrook surfacewater salinitystress grassgrowth
```

Expected: clean build, pFUnit all-pass, regression **4/4 byte-for-byte**. Any deviation = task NOT complete; fix or escalate. Do not commit on red.

**Phase close gate:** `pixi run check-full` 5/5 byte-for-byte.

**Subagent dispatch convention:** every implementer prompt MUST include this VG verbatim. The clean rebuild (`rm -rf builddir`) is mandatory because Meson's incremental build does not propagate `.mod` deps across the `swap_modern` / `swap_legacy` static-library boundary when state schema changes.

---

## File Structure

**New module:**
- `src/state/mesh_state.f90` — `mesh_state_t` (numnod + 6 allocatable mesh arrays + type-bound `init`)

**Extended modules:**
- `src/state/swap_state.f90` — adds `type(mesh_state_t) :: mesh` field
- `src/state/soilwater_state.f90` — +8 layer-flat fields, +4 runtime scalars
- `src/state/drainage_state.f90` — +10 drainage geo/switch fields
- `src/state/heat_state.f90` — adds `procedure :: init`

**Migrated readers (full `use variables` retirement):**
- `src/heat/frozencond.f90` (2 `use variables` sites)
- `src/heat/temperature.f90` (2 sites — bare + `only: NumNod` helper)
- `src/boundary/boundbottom.f90` (1 bare site)
- `src/boundary/boundtop.f90` (2 sites — bare in `boundtop` + `only:` in `PONDRUNOFF`; retains one narrow `only: nird` deferred to Arc 8)

**Codebase sweep (mesh / layer-flat / drainage-geo references swap to state):**
- `src/core/`: `initialize.f90`, `timecontrol_mod.f90`, `swap_bmi_mod.f90`, `swap_mod.f90`, `swap_capi_mod.f90`
- `src/atmosphere/`: `et.f90`, `interception.f90`, `meteoday.f90`, `meteodt.f90`
- `src/io/`: `readmeteo.f90`, `swap_csv_output.f90`, `swapoutput.f90`
- `src/soil/`: `soilgrid.f90`, `soilhydraulics.f90`, `waterbalance.f90`
- `src/drainage/`: `drainage.f90`, `surfacewater.f90`, `divdra.f90`
- `src/crop/`: `cropfixed_init.f90`, `cropgrass_init.f90`, `cropwofost_init.f90`, `wofost_soil_parameters.f90`, `cropgrowth.f90`, `irrigation.f90`, `oxygenstress.f90`, `rootextraction.f90`, `tillage.f90`

**Populator + adapter shrinks:**
- `src/io/toml/config_to_variables.f90` — adds Phase A dual-writes; loses legacy global writes in Phase C
- `src/core/initialize.f90` — loses zero-fills for deleted globals
- `src/core/variables.f90` — declarations of mesh + layer-flats + drainage-geo deleted

**Build:**
- `meson.build` — add `src/state/mesh_state.f90` to legacy sources list, **before** `src/state/swap_state.f90`
- `tests/unit/meson.build` — add to `pfunit_extra_sources`, same ordering

---

## Phase A — Schema + Dual-Write (Tasks 1–10)

Additive only. After Phase A: legacy globals untouched, state mirrors populated, no reader uses state yet, byte-for-byte trivial.

### Task 1: Pre-flight baseline

**Files:** No code changes; capture starting state.

- [ ] **Step 1: Clean rebuild.** Run: `rm -rf builddir && pixi run build-linux` — must succeed.
- [ ] **Step 2: pFUnit baseline.** Run: `pixi run test-pfunit` — record total pass count.
- [ ] **Step 3: check-full baseline.** Run: `pixi run check-full` — 5/5 pass.
- [ ] **Step 4: Empty marker commit.**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-bh): pre-flight baseline — boundary+heat globals retirement begins

check-full 5/5. Baseline locked before mesh extraction + soilwater/drainage
state extensions + 4-file reader migration + codebase sweep.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Create `mesh_state_mod`

**Files:**
- Create: `src/state/mesh_state.f90`
- Modify: `src/state/swap_state.f90` (add `type(mesh_state_t) :: mesh` field; add `use mesh_state_mod`)
- Modify: `meson.build` (add new file to legacy sources, **before** `src/state/swap_state.f90`)
- Modify: `tests/unit/meson.build` (add to `pfunit_extra_sources`, same ordering)

- [ ] **Step 1: Create the new module.**

```fortran
!> @file mesh_state.f90
!! SS-GR-BH: typed mesh / vertical-discretization state record.
!! Replaces the bare globals `numnod`, `layer(:)`, `dz(:)`, `z(:)`,
!! `disnod(:)`, `ztopcp(:)`, `zbotcp(:)` from `variables.f90`.
!! Populated once by `mesh_init` from `config_to_variables.f90`;
!! treated as read-only thereafter.
module mesh_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: mesh_state_t

   type :: mesh_state_t
      integer :: numnod = 0
      integer,      allocatable :: layer(:)
      real(real64), allocatable :: dz(:)
      real(real64), allocatable :: z(:)
      real(real64), allocatable :: disnod(:)
      real(real64), allocatable :: ztopcp(:)
      real(real64), allocatable :: zbotcp(:)
   contains
      procedure :: init => mesh_init
   end type mesh_state_t

contains

   subroutine mesh_init(self, numnod_in, dz_in, z_in, disnod_in, &
                        ztopcp_in, zbotcp_in, layer_in)
      class(mesh_state_t), intent(inout) :: self
      integer,             intent(in)    :: numnod_in
      real(real64),        intent(in)    :: dz_in(:), z_in(:), disnod_in(:), &
                                            ztopcp_in(:), zbotcp_in(:)
      integer,             intent(in)    :: layer_in(:)

      self%numnod = numnod_in

      if (allocated(self%dz))     deallocate(self%dz)
      if (allocated(self%z))      deallocate(self%z)
      if (allocated(self%disnod)) deallocate(self%disnod)
      if (allocated(self%ztopcp)) deallocate(self%ztopcp)
      if (allocated(self%zbotcp)) deallocate(self%zbotcp)
      if (allocated(self%layer))  deallocate(self%layer)

      allocate(self%dz(numnod_in))
      allocate(self%z(numnod_in))
      allocate(self%disnod(numnod_in+1))
      allocate(self%ztopcp(numnod_in))
      allocate(self%zbotcp(numnod_in))
      allocate(self%layer(numnod_in))

      self%dz(:)     = dz_in(1:numnod_in)
      self%z(:)      = z_in(1:numnod_in)
      self%disnod(:) = disnod_in(1:numnod_in+1)
      self%ztopcp(:) = ztopcp_in(1:numnod_in)
      self%zbotcp(:) = zbotcp_in(1:numnod_in)
      self%layer(:)  = layer_in(1:numnod_in)
   end subroutine mesh_init

end module mesh_state_mod
```

**Note:** `disnod` is sized `numnod+1` per legacy convention (`disnod(numnod+1)` is referenced in `boundbottom.f90:180`). All other mesh arrays are `numnod`-sized.

- [ ] **Step 2: Wire into `swap_state_mod`.** Edit `src/state/swap_state.f90`: add `use mesh_state_mod, only: mesh_state_t` near other state-record use lines (line ~17), and add `type(mesh_state_t) :: mesh` near the other subrecord fields (after `type(timecontrol_state_t) :: timecontrol`).

- [ ] **Step 3: Add to `meson.build`.** Find the legacy sources list and add `'src/state/mesh_state.f90',` **immediately before** `'src/state/swap_state.f90'`. The dependency ordering matters — `swap_state_mod` `use`s `mesh_state_mod`.

- [ ] **Step 4: Add to `tests/unit/meson.build`.** Find `pfunit_extra_sources` and add `'src/state/mesh_state.f90',` before `'src/state/swap_state.f90'`.

- [ ] **Step 5: Run VG.** Expected: clean build, pFUnit pass, regression 4/4.

- [ ] **Step 6: Commit.**

```bash
git add src/state/mesh_state.f90 src/state/swap_state.f90 meson.build tests/unit/meson.build
git commit -m "$(cat <<'EOF'
schema(gr-bh): add mesh_state_mod + state%mesh field

New typed mesh subrecord — numnod + layer(:) + dz(:)/z(:)/disnod(:)/
ztopcp(:)/zbotcp(:). disnod sized numnod+1 per legacy convention.
Allocated/copied via state%mesh%init(...). No reader uses it yet
(allocated but unused) — schema is additive.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: `state%mesh` dual-write populator

**Files:** `src/io/toml/config_to_variables.f90`

Locate the section where mesh globals (`numnod`, `dz`, `z`, `disnod`, `ztopcp`, `zbotcp`, `layer`) are populated. Add a `state%mesh%init(...)` call **immediately after** those globals are fully written.

- [ ] **Step 1: Find the mesh population block.** Run: `grep -n "numnod\s*=\|disnod(numnod" src/io/toml/config_to_variables.f90 | head -10`. Identify the line where `numnod`, `dz`, `disnod`, `z`, `ztopcp`, `zbotcp`, `layer` are all populated.

- [ ] **Step 2: Add the dual-write call** immediately after the mesh globals are written:

```fortran
   ! [SS-GR-BH A3] dual-write: populate state%mesh alongside legacy mesh globals
   call state%mesh%init(numnod, dz, z, disnod, ztopcp, zbotcp, layer)
```

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/io/toml/config_to_variables.f90
git commit -m "$(cat <<'EOF'
schema(gr-bh): dual-write state%mesh in config_to_variables

state%mesh%init(...) called after legacy mesh globals are populated.
state and legacy globals carry bitwise-identical mesh data; no readers
have migrated yet.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: Extend `soilwater_state_t` — layer flats (8 fields)

**Files:** `src/state/soilwater_state.f90`

Add 8 allocatable layer-flat fields. These are `maho`-sized layer-indexed arrays (one entry per soil layer, not per node).

- [ ] **Step 1: Add fields.** Edit `src/state/soilwater_state.f90`. Locate the existing layer-indexed field block (e.g. near where `vg_params(:)` lives) and add:

```fortran
      ! [SS-GR-BH A4] layer-flat fields — maho-sized, one per soil layer
      real(real64), allocatable :: ksatexm(:)    !! layer Ksat (examined extension)
      real(real64), allocatable :: ksatfit(:)    !! layer fitted Ksat
      real(real64), allocatable :: cofani(:)     !! layer anisotropy coefficient
      logical                   :: flksatexm = .false.  !! global flag: Ksatexm present in input
      real(real64), allocatable :: orgmat(:)     !! layer gravimetric organic matter
      real(real64), allocatable :: psand(:)      !! layer sand fraction
      real(real64), allocatable :: psilt(:)      !! layer silt fraction
      real(real64), allocatable :: pclay(:)      !! layer clay fraction
```

- [ ] **Step 2: Run VG** (full regression — schema-only change, all 4/4 must still pass).

- [ ] **Step 3: Commit.**

```bash
git add src/state/soilwater_state.f90
git commit -m "$(cat <<'EOF'
schema(gr-bh): extend soilwater_state_t with 8 layer flats

ksatexm/ksatfit/cofani/orgmat/psand/psilt/pclay (maho-sized allocatables)
+ flksatexm logical scalar. Unpopulated until Task 6 (dual-write).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: Extend `soilwater_state_t` — runtime scalars (4 fields)

**Files:** `src/state/soilwater_state.f90`

Add 4 runtime scalars that today live as bare globals shared across boundtop/PONDRUNOFF/boundbottom.

- [ ] **Step 1: Add fields.** Edit `src/state/soilwater_state.f90`. Locate the scalar field block and add:

```fortran
      ! [SS-GR-BH A5] runtime scalars formerly bare globals (boundtop/PONDRUNOFF/boundbottom)
      real(real64) :: q0     = 0.0_real64    !! surface flux (precip + runon - reva) [cm/d]
      real(real64) :: k1max  = 0.0_real64    !! max conductivity at z=0 [cm/d]
      real(real64) :: H0max  = 0.0_real64    !! max ponding pre-runoff [cm]
      integer      :: swbotb_runtime = 0      !! runtime-overridable bottom-boundary switch
```

`swbotb_runtime` is seeded in Task 7 from `config%bottom_boundary%swbotb`. Legacy code mutates the bare global `swbotb` to `-2` at runtime when bottom node goes oven-dry; the runtime copy lives on state. Config value stays immutable.

- [ ] **Step 2: Run VG.**

- [ ] **Step 3: Commit.**

```bash
git add src/state/soilwater_state.f90
git commit -m "$(cat <<'EOF'
schema(gr-bh): extend soilwater_state_t with 4 runtime scalars

q0/k1max/H0max (cross-subroutine state previously held in bare globals)
+ swbotb_runtime (legacy mutates global swbotb=-2 at runtime; this is
now the runtime-mutable copy).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: Dual-write soilwater layer flats

**Files:** `src/io/toml/config_to_variables.f90`

Mirror the 8 layer-flat fields into `state%soilwater` after legacy globals are populated.

- [ ] **Step 1: Find the layer-flat population block.** Run: `grep -n "ksatexm\s*=\s*config\|ksatfit\s*=\s*config\|cofani(\s*[0-9i]" src/io/toml/config_to_variables.f90 | head -20`.

- [ ] **Step 2: After each block where legacy globals are populated, add dual-write.** Use `move_alloc` for arrays (cleaner than allocate+copy) or explicit allocate+copy. Add these blocks after the corresponding legacy assignments:

```fortran
   ! [SS-GR-BH A6] dual-write: state%soilwater layer flats
   if (allocated(ksatexm)) then
      if (allocated(state%soilwater%ksatexm)) deallocate(state%soilwater%ksatexm)
      allocate(state%soilwater%ksatexm(size(ksatexm)))
      state%soilwater%ksatexm = ksatexm
   end if
   if (allocated(ksatfit)) then
      if (allocated(state%soilwater%ksatfit)) deallocate(state%soilwater%ksatfit)
      allocate(state%soilwater%ksatfit(size(ksatfit)))
      state%soilwater%ksatfit = ksatfit
   end if
   if (allocated(cofani)) then
      if (allocated(state%soilwater%cofani)) deallocate(state%soilwater%cofani)
      allocate(state%soilwater%cofani(size(cofani)))
      state%soilwater%cofani = cofani
   end if
   state%soilwater%flksatexm = flksatexm
   if (allocated(orgmat)) then
      if (allocated(state%soilwater%orgmat)) deallocate(state%soilwater%orgmat)
      allocate(state%soilwater%orgmat(size(orgmat)))
      state%soilwater%orgmat = orgmat
   end if
   if (allocated(psand)) then
      if (allocated(state%soilwater%psand)) deallocate(state%soilwater%psand)
      allocate(state%soilwater%psand(size(psand)))
      state%soilwater%psand = psand
   end if
   if (allocated(psilt)) then
      if (allocated(state%soilwater%psilt)) deallocate(state%soilwater%psilt)
      allocate(state%soilwater%psilt(size(psilt)))
      state%soilwater%psilt = psilt
   end if
   if (allocated(pclay)) then
      if (allocated(state%soilwater%pclay)) deallocate(state%soilwater%pclay)
      allocate(state%soilwater%pclay(size(pclay)))
      state%soilwater%pclay = pclay
   end if
```

**Placement note:** Some legacy fields are non-allocatable fixed-size arrays (e.g., `orgmat(maho)`). For those, drop the `if (allocated(...))` guards and use `allocate(state%soilwater%X(maho))`, then copy element-wise. Inspect each declaration in `src/core/variables.f90` (lines 800–1115) before writing the dual-write block.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/io/toml/config_to_variables.f90
git commit -m "$(cat <<'EOF'
schema(gr-bh): dual-write soilwater layer flats

state%soilwater%(ksatexm|ksatfit|cofani|flksatexm|orgmat|psand|psilt|pclay)
populated alongside legacy globals. State carries bitwise-identical layer
data; no readers have migrated yet.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 7: Dual-write soilwater runtime scalars

**Files:** `src/io/toml/config_to_variables.f90`

Seed the 4 runtime scalars in `state%soilwater`.

- [ ] **Step 1: Locate** the section where `swbotb` is assigned from `config%bottom_boundary%swbotb` (~line 755).

- [ ] **Step 2: Add seeding.** Immediately after `swbotb = config%bottom_boundary%swbotb`:

```fortran
   ! [SS-GR-BH A7] dual-write: seed state%soilwater runtime scalars
   state%soilwater%swbotb_runtime = config%bottom_boundary%swbotb
   state%soilwater%q0    = 0.0_real64
   state%soilwater%k1max = 0.0_real64
   state%soilwater%H0max = 0.0_real64
```

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/io/toml/config_to_variables.f90
git commit -m "$(cat <<'EOF'
schema(gr-bh): seed soilwater runtime scalars in adapter

swbotb_runtime seeded from config%bottom_boundary%swbotb; q0/k1max/H0max
zero-initialized. State runtime scalars carry the same values as legacy
globals at simulation start.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 8: Extend `drainage_state_t` (10 fields)

**Files:** `src/state/drainage_state.f90`

Add 10 fields for drainage geometry + switches.

- [ ] **Step 1: Add fields** to the type body. Place near existing drainage fields:

```fortran
      ! [SS-GR-BH A8] drainage geometry + switches formerly in variables.f90
      integer :: nrlevs      = 0
      integer :: swdivd      = 0
      integer :: swnrsrf     = 0
      integer :: swtopnrsrf  = 0
      integer :: swdivdinf   = 0
      real(real64) :: FacDpthInf = 0.0_real64
      real(real64), allocatable :: L(:)       !! drainage spacing per level [cm]
      real(real64), allocatable :: zbotdr(:)  !! drainage depth per level [cm]
      real(real64), allocatable :: owltab(:)  !! open-water level table
```

**Note:** `swdra` already lives on `state%surfacewater%swdra` (added in GR-UTILS Task 4). Check via `grep -n "swdra" src/state/surfacewater_state.f90` — if present there, the drainage_state copy is redundant. Decision: **swdra stays on state%surfacewater** (per Arc 1 precedent), so do NOT add it to drainage_state. Adjust the field count in commit message.

- [ ] **Step 2: Run VG.**

- [ ] **Step 3: Commit.**

```bash
git add src/state/drainage_state.f90
git commit -m "$(cat <<'EOF'
schema(gr-bh): extend drainage_state_t with 9 fields

nrlevs/swdivd/swnrsrf/swtopnrsrf/swdivdinf/FacDpthInf scalars +
L(:)/zbotdr(:)/owltab(:) allocatables. (swdra stays on state%surfacewater
per GR-UTILS precedent.) Unpopulated until Task 9 (dual-write).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 9: Dual-write drainage geometry

**Files:** `src/io/toml/config_to_variables.f90`

Mirror drainage geometry into `state%drainage` after legacy globals are populated.

- [ ] **Step 1: Find** drainage population block. Run: `grep -n "nrlevs\s*=\s*config\|swdivd\s*=\s*config\|zbotdr.*=\s*config" src/io/toml/config_to_variables.f90 | head -10`.

- [ ] **Step 2: Add dual-write block** after legacy assigns:

```fortran
   ! [SS-GR-BH A9] dual-write: state%drainage geometry + switches
   state%drainage%nrlevs     = nrlevs
   state%drainage%swdivd     = swdivd
   state%drainage%swnrsrf    = swnrsrf
   state%drainage%swtopnrsrf = swtopnrsrf
   state%drainage%swdivdinf  = swdivdinf
   state%drainage%FacDpthInf = FacDpthInf
   if (allocated(state%drainage%L))      deallocate(state%drainage%L)
   if (allocated(state%drainage%zbotdr)) deallocate(state%drainage%zbotdr)
   if (allocated(state%drainage%owltab)) deallocate(state%drainage%owltab)
   allocate(state%drainage%L(size(L)))
   allocate(state%drainage%zbotdr(size(zbotdr)))
   allocate(state%drainage%owltab(size(owltab)))
   state%drainage%L      = L
   state%drainage%zbotdr = zbotdr
   state%drainage%owltab = owltab
```

**Per-field caveat:** check declarations in `src/core/variables.f90` for each — some may be fixed-size (`zbotdr(madr)` etc.) rather than allocatable. Adjust `allocate(...)` calls accordingly: allocate to the dimension parameter from `swap_array_dimensions` when source is fixed-size.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/io/toml/config_to_variables.f90
git commit -m "$(cat <<'EOF'
schema(gr-bh): dual-write drainage geo in adapter

state%drainage%(nrlevs|swdivd|swnrsrf|swtopnrsrf|swdivdinf|FacDpthInf|
L|zbotdr|owltab) populated alongside legacy globals. No readers migrated
yet.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 10: Phase A close — check-full 5/5

**Files:** None (verification only).

- [ ] **Step 1: Clean rebuild.** `rm -rf builddir && pixi run build-linux`.
- [ ] **Step 2: pFUnit.** `pixi run test-pfunit` — all-pass.
- [ ] **Step 3: check-full.** `pixi run check-full` — **5/5 byte-for-byte**.
- [ ] **Step 4: Sanity greps.**
  - `grep -n "state%mesh%init" src/io/toml/config_to_variables.f90` — ≥1 hit.
  - `grep -n "state%soilwater%ksatexm\s*=\|state%soilwater%orgmat\s*=" src/io/toml/config_to_variables.f90` — ≥2 hits.
  - `grep -n "state%drainage%nrlevs\s*=" src/io/toml/config_to_variables.f90` — ≥1 hit.
- [ ] **Step 5: Phase A close marker.**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-bh): Phase A complete — schema + dual-write landed

state%mesh, state%soilwater (+12 fields), state%drainage (+9 fields)
populated alongside legacy globals. check-full 5/5 byte-for-byte.
No readers have migrated yet — every legacy global still load-bearing.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Phase B — Reader Cutover + heat_init Relocation (Tasks 11–23)

The 4 target files drop `use variables` and read state/config exclusively. heat_init moves to type-bound. After Phase B: `grep "^[[:space:]]*use variables" src/heat/ src/boundary/` returns at most one expected hit (`only: nird` in boundtop, deferred to Arc 8).

### Task 11: Relocate `heat_init` to type-bound `state%heat%init`

**Files:**
- Modify: `src/state/heat_state.f90` — add `procedure :: init`, contains block, and the relocated subroutine
- Modify: `src/heat/temperature.f90` — delete `heat_init` (lines ~535–560) and remove from `public ::` list
- Modify: `src/core/swap_mod.f90` — change call site from `call heat_init(state)` to `call state%heat%init(config%heat, state%mesh%numnod)`

- [ ] **Step 1: Add type-bound procedure to `heat_state_mod`.**

Edit `src/state/heat_state.f90`:

```fortran
   ! ... existing type :: heat_state_t with fields ...
      character(len=32), allocatable :: output_columns(:)
      integer                        :: output_n_cols = 0

   contains
      procedure :: init => heat_state_init
   end type heat_state_t

contains

   subroutine heat_state_init(self, heat_cfg, numnod_in)
      use, intrinsic :: iso_fortran_env, only: real64
      use heat_config_mod, only: heat_config_t
      class(heat_state_t),    intent(inout) :: self
      type(heat_config_t),    intent(in)    :: heat_cfg
      integer,                intent(in)    :: numnod_in

      ! heat_cfg currently unused; signature reserved for future seed migration.
      if (.not. allocated(self%tsoil))   allocate(self%tsoil(numnod_in))
      if (.not. allocated(self%heacap))  allocate(self%heacap(numnod_in))
      if (.not. allocated(self%heacon))  allocate(self%heacon(numnod_in))
      if (.not. allocated(self%rfcp))    allocate(self%rfcp(numnod_in))
      if (.not. allocated(self%fquartz)) allocate(self%fquartz(numnod_in))
      if (.not. allocated(self%fclay))   allocate(self%fclay(numnod_in))
      if (.not. allocated(self%forg))    allocate(self%forg(numnod_in))

      self%tsoil   = 0.0_real64
      self%heacap  = 0.0_real64
      self%heacon  = 0.0_real64
      self%rfcp    = 1.0_real64    ! NB: 1.0 not 0.0 — matches legacy initial value
      self%fquartz = 0.0_real64
      self%fclay   = 0.0_real64
      self%forg    = 0.0_real64
   end subroutine heat_state_init

end module heat_state_mod
```

Add `use heat_config_mod, only: heat_config_t` to the imports near the top of the module.

- [ ] **Step 2: Delete `heat_init` from `temperature.f90`** (lines ~535–560). Also remove `heat_init` from the `public ::` line at the top of `temperature_mod` (around line 41). Replace `public :: temperature, devries, heat_init` with `public :: temperature, devries`.

- [ ] **Step 3: Update call site in `swap_mod.f90`.** Find the existing call (likely `call heat_init(state)`):

```fortran
   ! Before:
   call heat_init(state)
   ! After:
   call state%heat%init(config%heat, state%mesh%numnod)
```

Also remove the `use temperature_mod, only: heat_init` if it exists (or narrow it to `temperature, devries`).

- [ ] **Step 4: Run VG.**

- [ ] **Step 5: Commit.**

```bash
git add src/state/heat_state.f90 src/heat/temperature.f90 src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
refactor(gr-bh): relocate heat_init to type-bound state%heat%init

Mirrors surfacewater_state pilot pattern. Signature becomes
state%heat%init(config%heat, numnod) — config%heat reserved for future
seed migration. Call site in swap_mod.f90 updated.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 12: Migrate `FrozenCond`

**Files:** `src/heat/frozencond.f90` (lines ~70–155, the `FrozenCond` subroutine)

`FrozenCond` reads: `numnod`, `swfrost`, `tfroststa`, `tfrostend`, `z`, `disnod` from `use variables`.

- [ ] **Step 1: Replace the `use variables` line at line ~71** with:

```fortran
    use swap_state_mod, only: swap_state_t
    use swap_config_mod, only: swap_config_t
```

(Remove `use variables`. Keep the existing `use swap_state_mod` line — merge if duplicated.)

- [ ] **Step 2: Add `config` to the signature.**

```fortran
   subroutine FrozenCond(state, config)
     use swap_state_mod, only: swap_state_t
     use swap_config_mod, only: swap_config_t
     implicit none
     type(swap_state_t),  intent(inout) :: state
     type(swap_config_t), intent(in)    :: config
```

- [ ] **Step 3: Replace references inside the body:**

| Was | Now |
|---|---|
| `numnod` | `state%mesh%numnod` |
| `swfrost` | `config%soil%frost%swfrost` |
| `tfroststa` | `config%heat%tfroststa` |
| `tfrostend` | `config%heat%tfrostend` |
| `z(...)` | `state%mesh%z(...)` |
| `disnod(...)` | `state%mesh%disnod(...)` |

Use `replace_all` carefully — `z` is a single-character name and may collide. Use a context-aware grep/sed: `grep -n "\bz(" src/heat/frozencond.f90` first, then substitute.

- [ ] **Step 4: Update caller** in `swap_mod.f90`. Find `call FrozenCond(state)` and change to `call FrozenCond(state, config)`.

- [ ] **Step 5: Run VG.**

- [ ] **Step 6: Commit.**

```bash
git add src/heat/frozencond.f90 src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
refactor(gr-bh): FrozenCond — drop use variables (Phase B reader cutover)

Mesh refs → state%mesh%X. swfrost → config%soil%frost; tfroststa/end →
config%heat. Signature gains config arg.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 13: Migrate `FrozenBounds`

**Files:** `src/heat/frozencond.f90` (lines ~193–337, the `FrozenBounds` subroutine)

`FrozenBounds` reads: `swdra`, `macp` (array dim, parameter), `nrlevs`, `zbotdr`, `layer`, `ksatexm`, `ksatfit`, `cofani`, `dz`, `swdivd`, `Swdivdinf`, `Swnrsrf`, `SwTopnrsrf`, `FacDpthInf`, `owltab`, `L`.

- [ ] **Step 1: Replace the `use variables` at line ~194** with state/config imports.

```fortran
    use swap_state_mod, only: swap_state_t
    use swap_config_mod, only: swap_config_t
    use distribute_drainage, only: DIVDRA
    use swap_array_dimensions, only: macp
```

Add `config` to the signature (`type(swap_config_t), intent(in) :: config`).

- [ ] **Step 2: Replace references:**

| Was | Now |
|---|---|
| `swdra` | `state%surfacewater%swdra` |
| `nrlevs` | `state%drainage%nrlevs` |
| `zbotdr(level)` | `state%drainage%zbotdr(level)` |
| `layer(node)` | `state%mesh%layer(node)` |
| `ksatexm(...)` | `state%soilwater%ksatexm(...)` |
| `ksatfit(...)` | `state%soilwater%ksatfit(...)` |
| `cofani(...)` | `state%soilwater%cofani(...)` |
| `dz(...)` | `state%mesh%dz(...)` |
| `swdivd` | `state%drainage%swdivd` |
| `Swdivdinf` | `state%drainage%swdivdinf` |
| `Swnrsrf` | `state%drainage%swnrsrf` |
| `SwTopnrsrf` | `state%drainage%swtopnrsrf` |
| `FacDpthInf` | `state%drainage%FacDpthInf` |
| `owltab` | `state%drainage%owltab` |
| `L` | `state%drainage%L` |

**DIVDRA call:** existing signature passes `dz`, `L`, `owltab`, `Zbotdr`, `FacDpthInf` positionally. Update the call to read from state instead:

```fortran
call divdra (state%mesh%numnod, state%drainage%nrlevs, state%mesh%dz, &
             ksatcp, ksatcp, sw_fluseksatexm, &
             layercp, cofanicp, ztop, state%drainage%L, qdrain, qdra, &
             state%drainage%swdivdinf, state%drainage%swnrsrf, &
             state%drainage%swtopnrsrf, state%drainage%zbotdr, &
             tc_dt, state%drainage%FacDpthInf, state%drainage%owltab, tc_t1900)
```

(DIVDRA's own `use variables` migration is deferred to Task 30 / Phase C C5.)

- [ ] **Step 3: Update caller** in `swap_mod.f90`. `call FrozenBounds(state)` → `call FrozenBounds(state, config)`.

- [ ] **Step 4: Run VG.**

- [ ] **Step 5: Commit.**

```bash
git add src/heat/frozencond.f90 src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
refactor(gr-bh): FrozenBounds — drop use variables

Mesh/drainage/layer-flat refs → state%(mesh|drainage|soilwater)%X.
DIVDRA call now reads state%drainage geometry; DIVDRA-internal migration
deferred to Phase C C5 (drainage cluster).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 14: Verify `frozencond.f90` clean

**Files:** None (verification only).

- [ ] **Step 1: Grep.** Run: `grep -n "^[[:space:]]*use variables" src/heat/frozencond.f90`. **Expected:** empty.
- [ ] **Step 2: Run VG.**
- [ ] **Step 3: No commit needed** (verification gate only).

If grep returns ANY hit, return to Task 12/13 and finish migration. Do not advance.

---

### Task 15: Migrate `temperature`

**Files:** `src/heat/temperature.f90` (the `temperature` subroutine, lines ~82–279)

Symbols read via `use variables` at line ~83: `macp`, `swcalt`, `swinco`, `nheat`, `tsoil(init)`, `zh`, `numnod`, `z`, `layer`, `orgmat`, `psand`, `psilt`, `pclay`, `tmean`, `tampli`, `ddamp`, `timref`, `swtopbhea`, `swbotbhea`, `temtoptab`, `tembtab`, `mabbc`, `Tav`, `atav`, `dz`, `disnod`.

- [ ] **Step 1: Replace `use variables` at line 83** with explicit imports:

```fortran
      use swap_state_mod,        only: swap_state_t
      use swap_config_mod,       only: swap_config_t
      use array_utils,           only: afgen
      use numericalsolvers_mod,  only: tridag
      use swap_array_dimensions, only: macp, mabbc
```

- [ ] **Step 2: Add `config` to signature.**

```fortran
   subroutine temperature(task, state, config)
     ! ... imports above ...
     implicit none
     integer,             intent(in)    :: task
     type(swap_state_t),  intent(inout) :: state
     type(swap_config_t), intent(in)    :: config
```

- [ ] **Step 3: Replace references** inside the body:

| Was | Now |
|---|---|
| `numnod` | `state%mesh%numnod` |
| `z(i)` | `state%mesh%z(i)` |
| `dz(i)` | `state%mesh%dz(i)` |
| `disnod(i)` | `state%mesh%disnod(i)` |
| `layer(i)` | `state%mesh%layer(i)` |
| `orgmat(lay)` | `state%soilwater%orgmat(lay)` |
| `psand(lay)`, `psilt(lay)`, `pclay(lay)` | `state%soilwater%psand/psilt/pclay(lay)` |
| `swcalt`, `swinco` | `config%heat%swcalt`, `config%simulation%swinco` |
| `nheat`, `tsoil(init)`, `zh(init)` | `config%heat%nheat`, `config%heat%tsoil_init`, `config%heat%zh` (check actual names in `heat_config_mod`) |
| `tmean`, `tampli`, `ddamp`, `timref` | `config%heat%tmean`, `config%heat%tampli`, `config%heat%ddamp`, `config%heat%timref` |
| `swtopbhea`, `swbotbhea` | `config%heat%swtopbhea`, `config%heat%swbotbhea` |
| `temtoptab`, `tembtab` | `config%heat%temtoptab`, `config%heat%tembtab` |
| `Tav` | `state%atmosphere%Tav` |
| `atav(state%timecontrol%wrecord)` | `state%atmosphere%atav(state%timecontrol%wrecord)` (verify `atav` is on atmosphere_state) |

**Caveat — config field names:** Before substitution, verify each name in `src/config/heat_config.f90`. Names may differ slightly (e.g., `tsoil_init` may be `tsoil(:)` on the heat_config). Use `grep -n "tmean\|tampli\|ddamp\|nheat\|temtoptab" src/config/heat_config.f90` to confirm.

**Caveat — `Tav`/`atav` location:** Run `grep -n "Tav\b\|atav\b" src/state/atmosphere_state.f90` to confirm. If absent, they're still on variables.f90 — leave as `use variables, only: Tav, atav` and note as deferred to Arc 4 (atmosphere). Update plan close-criteria for Task 17 accordingly.

- [ ] **Step 4: Update `devries` call.** Existing call: `call devries(theave,heacap_loc,heacnd,ht_fquartz,ht_fclay,ht_forg,sw_thetas)`. Task 16 migrates `Devries` itself; for now leave the call as-is — it works because `sw_thetas` is already the typed alias.

- [ ] **Step 5: Update caller** in `swap_mod.f90`. `call temperature(task, state)` → `call temperature(task, state, config)`.

- [ ] **Step 6: Run VG.**

- [ ] **Step 7: Commit.**

```bash
git add src/heat/temperature.f90 src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
refactor(gr-bh): temperature() — drop use variables

Mesh refs → state%mesh%X. Heat config switches and tables → config%heat;
swinco → config%simulation. Soil composition → state%soilwater%
(orgmat/psand/psilt/pclay). Signature gains config arg.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 16: Migrate `Devries` helper

**Files:** `src/heat/temperature.f90` (the `Devries` subroutine, lines ~351–529)

Line 352 has `use variables, only: NumNod`. Replace with explicit `numnod` arg.

- [ ] **Step 1: Drop the `use variables` import** at line 352.

- [ ] **Step 2: Add `numnod` as explicit arg.**

```fortran
   subroutine Devries (numnod_in, theta, HeaCap, HeaCon, fquartz_in, fclay_in, forg_in, thetas_in)
     use swap_array_dimensions, only: macp
     implicit none
     integer, intent(in) :: numnod_in
     real(8) theta(macp)
     ! ... rest unchanged ...
```

Replace internal `NumNod` references with `numnod_in`.

- [ ] **Step 3: Update caller** in `temperature` (Task 15 left it unchanged):

```fortran
   call devries(state%mesh%numnod, theave, heacap_loc, heacnd, &
                ht_fquartz, ht_fclay, ht_forg, sw_thetas)
```

- [ ] **Step 4: Run VG.**

- [ ] **Step 5: Commit.**

```bash
git add src/heat/temperature.f90
git commit -m "$(cat <<'EOF'
refactor(gr-bh): Devries — drop use variables, explicit numnod arg

Pure helper now has fully explicit signature. NumNod passed by caller
from state%mesh%numnod.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 17: Verify `temperature.f90` clean

**Files:** None (verification only).

- [ ] **Step 1: Grep.** Run: `grep -n "^[[:space:]]*use variables" src/heat/temperature.f90`. **Expected:** empty (unless Task 15 left a documented `use variables, only: Tav, atav` deferral — verify the deferral is captured in Phase B close criteria).
- [ ] **Step 2: Run VG.**
- [ ] **Step 3: No commit.**

---

### Task 18: Migrate `BoundBottom`

**Files:** `src/boundary/boundbottom.f90`

Symbols read via `use variables` at line 12: `swbotb` (mutated!), `gwltab`, `mabbc`, `sw2`, `sinave`, `sinamp`, `sinmax`, `qbotab`, `sw3`, `aqave`, `aqamp`, `aqper`, `aqtmax`, `haqtab`, `SwBotb3ResVert`, `hdrain`, `shape`, `rimlay`, `sw4`, `swqhbot`, `cofqha`, `cofqhb`, `swcofqhc`, `cofqhc`, `hbotab`, `logf`.

- [ ] **Step 1: Replace `use variables` at line 12** with explicit imports:

```fortran
    use swap_state_mod,        only: swap_state_t
    use swap_config_mod,       only: swap_config_t
    use swap_log,              only: log_debug, to_str
    use swap_array_dimensions, only: mabbc
```

- [ ] **Step 2: Add `config` to signature** at line 55:

```fortran
    subroutine BoundBottom(state, config)
        use array_utils, only: afgen
        use soilhydraulics_utils, only: watcon, hconduc
        implicit none
        type(swap_state_t),  intent(inout) :: state
        type(swap_config_t), intent(in)    :: config
```

- [ ] **Step 3: Replace references:**

| Was | Now |
|---|---|
| `swbotb` (read AND write) | `state%soilwater%swbotb_runtime` |
| `gwltab`, `qbotab`, `haqtab`, `hbotab` | `config%bottom_boundary%(gwltab/qbotab/haqtab/hbotab)` |
| `sw2`, `sw3`, `sw4` | `config%bottom_boundary%(sw2/sw3/sw4)` |
| `sinave`, `sinamp`, `sinmax` | `config%bottom_boundary%(sinave/sinamp/sinmax)` |
| `aqave`, `aqamp`, `aqper`, `aqtmax` | `config%bottom_boundary%(aqave/aqamp/aqper/aqtmax)` |
| `SwBotb3ResVert` | `config%bottom_boundary%SwBotb3ResVert` |
| `hdrain`, `shape`, `rimlay` | `config%bottom_boundary%(hdrain/shape/rimlay)` |
| `swqhbot`, `cofqha`, `cofqhb`, `swcofqhc`, `cofqhc` | `config%bottom_boundary%(swqhbot/cofqha/cofqhb/swcofqhc/cofqhc)` |
| `numnod`, `dz(...)`, `ztopcp(...)`, `zbotcp(...)` | `state%mesh%(numnod/dz/ztopcp/zbotcp)` |
| `logf` | leave as-is — `logf` is the log unit number, used by `call warn(...)`. Add `use variables, only: logf` deferred to Arc 9 if no state placement is obvious. Note this in commit message. |

**Caveat — config field names:** verify in `src/config/bottom_boundary_config.f90`. Names may have different casing.

**Caveat — `swbotb` mutation:** the legacy code does `swbotb = -2` at line 102. The new code does `state%soilwater%swbotb_runtime = -2`. This is the key state-vs-config separation; the **immutable** `config%bottom_boundary%swbotb` is not touched after init.

- [ ] **Step 4: Update caller** in `swap_mod.f90`. `call BoundBottom(state)` → `call BoundBottom(state, config)`.

- [ ] **Step 5: Run VG.**

- [ ] **Step 6: Commit.**

```bash
git add src/boundary/boundbottom.f90 src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
refactor(gr-bh): BoundBottom — drop use variables

All bottom-boundary config switches/tables → config%bottom_boundary%X.
Mesh refs → state%mesh%X. swbotb mutation → state%soilwater%swbotb_runtime
(immutable config%bottom_boundary%swbotb preserved). logf stays via
narrow use-variables-only — deferred to Arc 9.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 19: Verify `boundbottom.f90` clean

**Files:** None (verification only).

- [ ] **Step 1: Grep.** `grep -n "^[[:space:]]*use variables" src/boundary/boundbottom.f90`. **Expected:** at most one hit, the `use variables, only: logf` deferred import.
- [ ] **Step 2: Run VG.**
- [ ] **Step 3: No commit.**

---

### Task 20: Migrate `boundtop`

**Files:** `src/boundary/boundtop.f90` (the `boundtop` subroutine, lines ~19–190)

Symbols read via bare `use variables` at line 20: `q0`, `k1max`, `flrunon`, `runonarr`, `dz`, `disnod`, `swkmean`, `swredu`, `nird`, `ArMpSs`, `swdra`, `ksatexm`, `ksatfit`, `layer`, `pondmx`.

- [ ] **Step 1: Replace `use variables` at line 20** with explicit imports:

```fortran
      use swap_state_mod,  only: swap_state_t
      use swap_config_mod, only: swap_config_t
      use swap_log,        only: log_debug, to_str
      use surfacewater_utils, only: runoff
      use variables,       only: nird   ! [SS-GR-BH B10] DEFERRED to Arc 8 (irrigation cluster)
```

- [ ] **Step 2: Add `config` to signature.**

```fortran
   subroutine boundtop(state, config)
   use soilhydraulics_utils, only: watcon, hconduc, hcomean
   implicit none
   type(swap_state_t),  intent(inout) :: state
   type(swap_config_t), intent(in)    :: config
```

- [ ] **Step 3: Replace references:**

| Was | Now |
|---|---|
| `q0`, `k1max` | `state%soilwater%(q0/k1max)` |
| `flrunon` | `config%simulation%flrunon` (verify in simulation_config — may be in general_config) |
| `runonarr(...)` | `config%simulation%runonarr(...)` |
| `dz(...)`, `disnod(...)`, `layer(...)` | `state%mesh%(dz/disnod/layer)` |
| `swkmean`, `swredu` | `config%soil%swkmean`, `config%soil%swredu` (verify) |
| `nird` | `nird` (unchanged — deferred via narrow `use variables, only: nird`) |
| `ArMpSs` | **DELETE the line `ArMpSs = 0.d0`** (legacy zero-write; macropore retired). Then replace the read at line 145: `q0 = (state%atmosphere%nraidt+nird+state%atmosphere%melt)*(1.0d0-ArMpSs) + ...` → `q0 = (state%atmosphere%nraidt+nird+state%atmosphere%melt) + state%soilwater%runon - state%soilwater%reva`. (The `* (1.0d0-0.0d0)` factor is identity; bit-equivalent.) |
| `swdra` | `state%surfacewater%swdra` (already on state from GR-UTILS) |
| `ksatexm(layer(1))`, `ksatfit(layer(1))` | `state%soilwater%ksatexm(state%mesh%layer(1))`, `state%soilwater%ksatfit(state%mesh%layer(1))` |
| `pondmx` | `state%surfacewater%pondmx` (already on state from GR-UTILS) |

**Caveat — `state%soilwater%q0` write:** the new code writes `state%soilwater%q0 = (state%atmosphere%nraidt + ...)` instead of the bare global `q0 =`. PONDRUNOFF (Task 21) reads this back via state. Cross-subroutine state crosses cleanly.

- [ ] **Step 4: Update caller** in `swap_mod.f90`. `call boundtop(state)` → `call boundtop(state, config)`.

- [ ] **Step 5: Run VG.**

- [ ] **Step 6: Commit.**

```bash
git add src/boundary/boundtop.f90 src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
refactor(gr-bh): boundtop — drop use variables (modulo nird)

q0/k1max → state%soilwater. ArMpSs assignment deleted (macropore retired,
factor *(1-0)=1 is identity — bit-equivalent). Mesh refs → state%mesh.
Layer flats → state%soilwater. Config switches (flrunon/runonarr/swkmean/
swredu) → config. nird retained via narrow use-only — deferred to Arc 8.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 21: Migrate `PONDRUNOFF`

**Files:** `src/boundary/boundtop.f90` (the `PONDRUNOFF` subroutine, lines ~194–290)

`use variables, only:` at line 205: `swdra`, `disnod`, `H0max`, `k1max`, `pondmx`, `q0`, `rsro`, `rsroexp`, `swpondmx`, `pondmxtab`.

- [ ] **Step 1: Replace the import** at line 205 with state/config-only:

```fortran
      use swap_state_mod,  only: swap_state_t
      use swap_config_mod, only: swap_config_t
      use array_utils,     only: afgen
      use surfacewater_utils, only: runoff
      use swap_array_dimensions, only: mairg
```

- [ ] **Step 2: Add `config` to signature.**

```fortran
   subroutine pondrunoff(state, config)
     ! ... imports ...
     implicit none
     type(swap_state_t),  intent(inout) :: state
     type(swap_config_t), intent(in)    :: config
```

- [ ] **Step 3: Replace references:**

| Was | Now |
|---|---|
| `swdra` | `state%surfacewater%swdra` |
| `disnod(1)` | `state%mesh%disnod(1)` |
| `H0max` | `state%soilwater%H0max` |
| `k1max` | `state%soilwater%k1max` |
| `pondmx` | `state%surfacewater%pondmx` |
| `q0` | `state%soilwater%q0` |
| `rsro` | `state%surfacewater%rsro` |
| `rsroexp` | `state%surfacewater%rsroexp` |
| `swpondmx` | `config%surface_water%swpondmx` (verify) |
| `pondmxtab` | `config%surface_water%pondmxtab` (verify) |

**Note — `pondmx` write:** the legacy code at line 226 writes `pondmx = afgen(...)`. This is a **runtime mutation** of a config-loaded array. Decision: the runtime-mutable `pondmx` lives on `state%surfacewater%pondmx` (already there per GR-UTILS Task 4). The config copy stays immutable. Update the assignment to `state%surfacewater%pondmx = afgen(config%surface_water%pondmxtab, 2*mairg, t1900+dt)`.

- [ ] **Step 4: Update caller** in `boundtop` (inside the same file). Find `call PONDRUNOFF(state)` and change to `call PONDRUNOFF(state, config)`.

- [ ] **Step 5: Also update** any other callers of PONDRUNOFF — grep `grep -rn "call.*PONDRUNOFF\|call.*pondrunoff" src/` to find them.

- [ ] **Step 6: Run VG.**

- [ ] **Step 7: Commit.**

```bash
git add src/boundary/boundtop.f90 src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
refactor(gr-bh): PONDRUNOFF — drop use variables

All refs → state (q0/k1max/H0max/pondmx/rsro/rsroexp/swdra on
state%(soilwater|surfacewater); disnod on state%mesh) or config
(swpondmx/pondmxtab). Cross-subroutine state with boundtop now flows
purely through state.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 22: Verify `boundtop.f90` clean

**Files:** None (verification only).

- [ ] **Step 1: Grep.** `grep -n "^[[:space:]]*use variables" src/boundary/boundtop.f90`. **Expected:** exactly one hit, the deferred `use variables, only: nird`.
- [ ] **Step 2: Run VG.**
- [ ] **Step 3: No commit.**

---

### Task 23: Phase B close — check-full 5/5

**Files:** None (verification only).

- [ ] **Step 1: Clean rebuild + check-full.** `rm -rf builddir && pixi run build-linux && pixi run check-full` — **5/5 byte-for-byte**.
- [ ] **Step 2: Sanity greps.**
  - `grep -rn "^[[:space:]]*use variables" src/heat/ src/boundary/ | grep -v "only: *nird\|only: *logf"` — empty.
  - `grep -n "heat_init" src/heat/temperature.f90` — empty (relocated).
  - `grep -n "procedure :: init" src/state/heat_state.f90` — ≥1 hit.
- [ ] **Step 3: Phase B close marker.**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-bh): Phase B complete — 4-file reader cutover + heat_init relocated

frozencond, temperature, boundbottom, boundtop migrated off use variables
(boundtop retains only: nird; boundbottom retains only: logf — both
deferred). heat_init now state%heat%init type-bound. check-full 5/5
byte-for-byte.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Phase C — Codebase Mesh Sweep + Global Retirement (Tasks 24–39)

Every other reader of mesh / migrated layer-flat / migrated drainage-geo globals migrates to state. Then bare globals are deleted from `variables.f90` (compile surfaces any missed reader).

**Scope discipline:** Phase C does NOT touch any other `use variables` symbols. Files keep their `use variables` clause for unrelated globals. Surgical replacement only — substituting names, not removing imports of unrelated names.

**Contingency for non-state-bearing files:** if a file reads `numnod` but doesn't take `state`, add `type(swap_state_t), intent(in) :: state` to its signature. All current callers hold state.

### Symbol replacement reference (used by Tasks 24–30)

| Was (bare global) | Now |
|---|---|
| `numnod` | `state%mesh%numnod` |
| `dz`, `dz(i)` | `state%mesh%dz`, `state%mesh%dz(i)` |
| `z`, `z(i)` | `state%mesh%z`, `state%mesh%z(i)` |
| `disnod`, `disnod(i)` | `state%mesh%disnod`, `state%mesh%disnod(i)` |
| `ztopcp`, `ztopcp(i)` | `state%mesh%ztopcp`, `state%mesh%ztopcp(i)` |
| `zbotcp`, `zbotcp(i)` | `state%mesh%zbotcp`, `state%mesh%zbotcp(i)` |
| `layer`, `layer(i)` | `state%mesh%layer`, `state%mesh%layer(i)` |
| `ksatexm`, `ksatfit`, `cofani`, `orgmat`, `psand`, `psilt`, `pclay`, `flksatexm` | `state%soilwater%X` |
| `nrlevs`, `swdivd`, `swnrsrf`, `swtopnrsrf`, `swdivdinf`, `FacDpthInf`, `L`, `zbotdr`, `owltab` | `state%drainage%X` |

### Task 24: Cluster `src/core/`

**Files:** `src/core/initialize.f90`, `src/core/timecontrol_mod.f90`, `src/core/swap_bmi_mod.f90`, `src/core/swap_mod.f90`, `src/core/swap_capi_mod.f90` (verify each takes state).

- [ ] **Step 1: Per-file inventory.** For each file, run: `grep -nw "numnod\|dz\|z\|disnod\|ztopcp\|zbotcp\|layer" <file> | grep -v "state%\|config%\|integer :: \|real(8) :: \|! "`. Note all bare references.
- [ ] **Step 2: Apply substitutions** per the reference table above. Leave `use variables` imports for OTHER symbols intact — only stop importing names that have moved to state. Edit the `use variables, only:` lines to remove migrated names; if the file uses bare `use variables`, leave the import statement alone but rewrite individual reads.
- [ ] **Step 3: Initialize.f90 caveat:** this file zero-fills many arrays at start of simulation. Mesh-array zero-fills (if any) become writes to `state%mesh%X`. But the `state%mesh%init` call in `config_to_variables.f90` happens AFTER initialize.f90 — verify initialize.f90 zero-fills run BEFORE config_to_variables population, otherwise the zero-fill is dead. If dead, leave it as-is for now; Task 35 (C12) handles cleanup.
- [ ] **Step 4: Run VG.**
- [ ] **Step 5: Commit.**

```bash
git add src/core/initialize.f90 src/core/timecontrol_mod.f90 src/core/swap_bmi_mod.f90 src/core/swap_mod.f90 src/core/swap_capi_mod.f90
git commit -m "$(cat <<'EOF'
sweep(gr-bh): src/core/ — mesh refs to state%mesh%X

initialize/timecontrol_mod/swap_bmi_mod/swap_mod/swap_capi_mod migrated
off bare numnod/dz/z/disnod/ztopcp/zbotcp/layer reads. Files retain
use variables for unrelated globals (deferred to other arcs).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 25: Cluster `src/atmosphere/`

**Files:** `src/atmosphere/et.f90`, `src/atmosphere/interception.f90`, `src/atmosphere/meteoday.f90`, `src/atmosphere/meteodt.f90`.

- [ ] **Step 1: Inventory.** `for f in src/atmosphere/{et,interception,meteoday,meteodt}.f90; do echo "=== $f ==="; grep -nw "numnod\|dz\|z\|disnod\|ztopcp\|zbotcp\|layer" "$f" | grep -v "state%\|config%"; done`.
- [ ] **Step 2: Apply substitutions** per reference table.
- [ ] **Step 3: State plumbing.** If any of these subroutines don't take state today, add `type(swap_state_t), intent(in) :: state` to the signature and update callers. Verify by grep on call sites: `grep -rn "call ETref\|call interception\|call meteoday\|call meteodt" src/`.
- [ ] **Step 4: Run VG.**
- [ ] **Step 5: Commit.**

```bash
git add src/atmosphere/et.f90 src/atmosphere/interception.f90 src/atmosphere/meteoday.f90 src/atmosphere/meteodt.f90
git commit -m "$(cat <<'EOF'
sweep(gr-bh): src/atmosphere/ — mesh refs to state%mesh%X

et/interception/meteoday/meteodt migrated off mesh globals. Files retain
use variables for unrelated atmospheric globals (Arc 4 territory).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 26: Cluster `src/io/`

**Files:** `src/io/readmeteo.f90`, `src/io/swap_csv_output.f90`, `src/io/swapoutput.f90`.

- [ ] **Step 1: Inventory.** Same grep pattern as Task 25 over `src/io/{readmeteo,swap_csv_output,swapoutput}.f90`. **Expect:** swapoutput has many sites (~12 `use variables` blocks; many will mention `numnod`).
- [ ] **Step 2: Apply substitutions.** Note that swapoutput.f90 has several `use variables, only:` blocks — each block's `only:` list may need editing to drop migrated names. Also handle `nrlevs` → `state%drainage%nrlevs` where present.
- [ ] **Step 3: Run VG.**
- [ ] **Step 4: Commit.**

```bash
git add src/io/readmeteo.f90 src/io/swap_csv_output.f90 src/io/swapoutput.f90
git commit -m "$(cat <<'EOF'
sweep(gr-bh): src/io/ — mesh + nrlevs refs to state%X

readmeteo/swap_csv_output/swapoutput migrated off mesh globals and
nrlevs (drainage geo). Files retain use variables for unrelated globals
(meteo arrays etc; Arc 4 / Arc 7 territory).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 27: Cluster `src/soil/`

**Files:** `src/soil/soilgrid.f90`, `src/soil/soilhydraulics.f90`, `src/soil/waterbalance.f90`.

- [ ] **Step 1: Inventory.** Grep both mesh and layer-flat names. Soilhydraulics in particular reads `orgmat`/`psand`/`psilt`/`pclay`/`ksatexm`/`ksatfit`/`cofani`/`flksatexm` heavily.
- [ ] **Step 2: Apply substitutions** — mesh and layer-flats. Be especially careful: in `soilhydraulics.f90`, some of these arrays are also being WRITTEN (it's the populator for `state%soilwater%vg_params`). If a value is computed in this file and later read elsewhere, ensure both legacy globals and state mirrors get the new value (write to BOTH until Task 36/C10 retires the global). Actually — since Phase A dual-write copies legacy globals into state mirrors at the END of `config_to_variables`, any runtime computation in `soilhydraulics.f90` that updates legacy globals must ALSO update state mirrors. Otherwise state will drift from legacy.

**Decision:** for any runtime mutation of a migrated global, dual-write to state immediately. Mark such sites with a `! [GR-BH C4 dual-write]` comment for future cleanup.

- [ ] **Step 3: Run VG.**
- [ ] **Step 4: Commit.**

```bash
git add src/soil/soilgrid.f90 src/soil/soilhydraulics.f90 src/soil/waterbalance.f90
git commit -m "$(cat <<'EOF'
sweep(gr-bh): src/soil/ — mesh + layer flats to state%X

soilgrid/soilhydraulics/waterbalance migrated off mesh globals and layer
flats (ksatexm/ksatfit/cofani/orgmat/psand/psilt/pclay/flksatexm).
Runtime mutators in soilhydraulics dual-write to state mirrors to keep
state and legacy in sync until Task 36 retirement.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 28: Cluster `src/drainage/`

**Files:** `src/drainage/drainage.f90`, `src/drainage/surfacewater.f90`, `src/drainage/divdra.f90`.

- [ ] **Step 1: Inventory.** Grep mesh + drainage geo + relevant layer flats (`cofani`, `ksatexm`, `ksatfit` appear here).
- [ ] **Step 2: Apply substitutions.** DIVDRA already receives drainage geo via Task 13's caller change — verify the signature matches (`numnod`, `nrlevs`, `dz`, `ksatcp`, `ksatcp`, `fluseksatexm`, `layercp`, `cofanicp`, `ztop`, `L`, `qdrain`, `qdra`, `swdivdinf`, `swnrsrf`, `swtopnrsrf`, `zbotdr`, `dt`, `FacDpthInf`, `owltab`, `t1900`). DIVDRA's BODY may still read these from `use variables` — migrate now.
- [ ] **Step 3: Run VG.**
- [ ] **Step 4: Commit.**

```bash
git add src/drainage/drainage.f90 src/drainage/surfacewater.f90 src/drainage/divdra.f90
git commit -m "$(cat <<'EOF'
sweep(gr-bh): src/drainage/ — mesh + drainage geo + layer flats

drainage/surfacewater/divdra migrated off mesh, drainage geo (nrlevs/
swdivd/swnrsrf/swtopnrsrf/swdivdinf/FacDpthInf/L/zbotdr/owltab), and
layer flats. divdra signature was already updated by FrozenBounds (Task
13); body now reads from incoming args / state instead of bare globals.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 29: Cluster `src/crop/` part 1 — init files

**Files:** `src/crop/cropfixed_init.f90`, `src/crop/cropgrass_init.f90`, `src/crop/cropwofost_init.f90`, `src/crop/wofost_soil_parameters.f90`.

- [ ] **Step 1: Inventory.** Grep mesh + layer-flat names. These files initialize crop-soil-water relations.
- [ ] **Step 2: Apply substitutions.**
- [ ] **Step 3: Run VG.**
- [ ] **Step 4: Commit.**

```bash
git add src/crop/cropfixed_init.f90 src/crop/cropgrass_init.f90 src/crop/cropwofost_init.f90 src/crop/wofost_soil_parameters.f90
git commit -m "$(cat <<'EOF'
sweep(gr-bh): src/crop/ init files — mesh + layer flats to state

cropfixed_init/cropgrass_init/cropwofost_init/wofost_soil_parameters
migrated off mesh and layer-flat globals. Files retain use variables
for crop globals (Arc 8 territory).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 30: Cluster `src/crop/` part 2 — runtime

**Files:** `src/crop/cropgrowth.f90`, `src/crop/irrigation.f90`, `src/crop/oxygenstress.f90`, `src/crop/rootextraction.f90`, `src/crop/tillage.f90`.

- [ ] **Step 1: Inventory.**
- [ ] **Step 2: Apply substitutions.** Be especially careful with `cropgrowth.f90` — it has the longest `use variables, only:` lists.
- [ ] **Step 3: Run VG.**
- [ ] **Step 4: Commit.**

```bash
git add src/crop/cropgrowth.f90 src/crop/irrigation.f90 src/crop/oxygenstress.f90 src/crop/rootextraction.f90 src/crop/tillage.f90
git commit -m "$(cat <<'EOF'
sweep(gr-bh): src/crop/ runtime — mesh + layer flats to state

cropgrowth/irrigation/oxygenstress/rootextraction/tillage migrated off
mesh and layer-flat globals. Files retain use variables for crop globals
(Arc 8 territory).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 31: Audit pass — find stale reads

**Files:** None (audit + targeted fixes).

- [ ] **Step 1: Audit grep.** Run:

```bash
grep -rnw "numnod\|dz\|disnod\|ztopcp\|zbotcp" src/ --include="*.f90" \
  | grep -v "state%\|config%\|variables.f90\|mesh_state\|! \|use variables\|integer :: \|real(8) :: \|real(real64) :: "
```

**Expected:** empty. Any hit is a missed bare reference — investigate and patch.

- [ ] **Step 2: Layer-flat audit.** Run:

```bash
grep -rnw "ksatexm\|ksatfit\|cofani\|flksatexm\|orgmat\|psand\|psilt\|pclay" src/ --include="*.f90" \
  | grep -v "state%\|config%\|variables.f90\|! \|use variables\|integer :: \|real(8) :: \|real(real64) :: \|logical :: "
```

Expected: empty.

- [ ] **Step 3: Drainage-geo audit.** Run:

```bash
grep -rnw "nrlevs\|swdivd\|swnrsrf\|swtopnrsrf\|swdivdinf\|FacDpthInf\|owltab" src/ --include="*.f90" \
  | grep -v "state%\|config%\|variables.f90\|! \|use variables\|integer :: \|real(8) :: \|real(real64) :: "
```

Expected: empty.

- [ ] **Step 4: Runtime-state audit** (q0, k1max, H0max, swbotb). Run:

```bash
grep -rnw "q0\|k1max\|H0max\|swbotb\b\|ArMpSs" src/ --include="*.f90" \
  | grep -v "state%\|config%\|variables.f90\|! \|use variables\|real(8) :: \|integer :: "
```

Expected: empty (any swbotb hit should be `state%soilwater%swbotb_runtime`).

- [ ] **Step 5: If any hits, patch in one file per commit.** Each commit message: `audit(gr-bh): <file> — patch stale <symbol> reference`.

- [ ] **Step 6: Final pass:** re-run all 4 grep commands above until each returns empty.

- [ ] **Step 7: Run VG** after each patch.

- [ ] **Step 8: Audit close marker** (no code change if no patches needed):

```bash
git commit --allow-empty -m "$(cat <<'EOF'
audit(gr-bh): codebase sweep clean — no stale mesh/layer-flat/drainage-geo reads

All 4 audit greps empty. Ready for Phase C global deletion (Tasks 32-37).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

**This task is a MANDATORY gate.** Do not proceed to Task 32 until all 4 audit greps return empty.

---

### Task 32: Drop mesh writes from `config_to_variables.f90`

**Files:** `src/io/toml/config_to_variables.f90`

Delete legacy global writes for mesh fields. Keep the `state%mesh%init(...)` call.

- [ ] **Step 1: Find legacy mesh writes.** Run: `grep -n "^[[:space:]]*\(numnod\|dz(\|z(\|disnod(\|ztopcp(\|zbotcp(\|layer(\)" src/io/toml/config_to_variables.f90`.

- [ ] **Step 2: Delete each write line** (the legacy `dz(...) = config%X`, `numnod = config%X` lines, etc.).

- [ ] **Step 3: Keep** the `call state%mesh%init(numnod, dz, z, disnod, ztopcp, zbotcp, layer)` call. **BUT:** since the legacy locals `dz`, `z`, etc. are about to be deleted from variables.f90 (Task 35), this call needs the source data from `config` directly. Change it to read from config:

```fortran
   call state%mesh%init(config%soil%mesh%numnod, &
                        config%soil%mesh%dz, config%soil%mesh%z, &
                        config%soil%mesh%disnod, config%soil%mesh%ztopcp, &
                        config%soil%mesh%zbotcp, config%soil%mesh%layer)
```

**Caveat:** verify the actual config path for mesh — it may be `config%soil` or `config%simulation` or a dedicated `config%mesh`. Use `grep -n "numnod\s*=\s*config\|mesh%numnod" src/io/toml/config_to_variables.f90` to find the source.

- [ ] **Step 4: Run VG.**

- [ ] **Step 5: Commit.**

```bash
git add src/io/toml/config_to_variables.f90
git commit -m "$(cat <<'EOF'
retire(gr-bh): drop legacy mesh writes from adapter

numnod/dz/z/disnod/ztopcp/zbotcp/layer no longer written to bare globals.
state%mesh%init now sourced directly from config. Sets up Task 35
deletion of mesh declarations from variables.f90.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 33: Drop layer-flats + runtime-state writes from `config_to_variables.f90`

**Files:** `src/io/toml/config_to_variables.f90`

Delete legacy global writes for layer flats and runtime scalars. Keep state-mirror writes.

- [ ] **Step 1: Find.** `grep -n "^[[:space:]]*\(ksatexm\b\|ksatfit\b\|cofani\b\|flksatexm\b\|orgmat\b\|psand\b\|psilt\b\|pclay\b\|q0\b\|k1max\b\|H0max\b\|ArMpSs\b\|swbotb\b\)" src/io/toml/config_to_variables.f90`.

- [ ] **Step 2: Delete legacy writes**, keep the dual-write blocks added in Tasks 6 and 7. Edit those blocks now to read from `config` directly rather than from locally-declared globals that are about to be deleted.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/io/toml/config_to_variables.f90
git commit -m "$(cat <<'EOF'
retire(gr-bh): drop legacy layer-flat + runtime writes from adapter

ksatexm/ksatfit/cofani/flksatexm/orgmat/psand/psilt/pclay (layer flats),
q0/k1max/H0max/ArMpSs/swbotb (runtime scalars) no longer written to
bare globals. state mirrors populated directly from config. Sets up
Task 36 deletion.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 34: Drop drainage-geo writes from `config_to_variables.f90`

**Files:** `src/io/toml/config_to_variables.f90`

Delete legacy global writes for drainage geo.

- [ ] **Step 1: Find.** `grep -n "^[[:space:]]*\(nrlevs\b\|swdivd\b\|swnrsrf\b\|swtopnrsrf\b\|swdivdinf\b\|FacDpthInf\b\|L(\|zbotdr(\|owltab(\)" src/io/toml/config_to_variables.f90`.

- [ ] **Step 2: Delete legacy writes**, keep state-mirror writes (Task 9). Source state writes from `config%drain%X` directly.

- [ ] **Step 3: Run VG.**

- [ ] **Step 4: Commit.**

```bash
git add src/io/toml/config_to_variables.f90
git commit -m "$(cat <<'EOF'
retire(gr-bh): drop legacy drainage-geo writes from adapter

nrlevs/swdivd/swnrsrf/swtopnrsrf/swdivdinf/FacDpthInf/L/zbotdr/owltab
no longer written to bare globals. state%drainage mirrors populated
directly from config%drain. Sets up Task 37 deletion.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 35: Delete mesh globals from `variables.f90`

**Files:** `src/core/variables.f90`

Delete bare global declarations. The build will surface any reader missed by the audit.

- [ ] **Step 1: Find declarations.** `grep -n "^[[:space:]]*integer[[:space:]]*numnod\b\|^[[:space:]]*real(8)[[:space:]]*dz(\|^[[:space:]]*real(8)[[:space:]]*z(\|^[[:space:]]*real(8)[[:space:]]*disnod(\|^[[:space:]]*real(8)[[:space:]]*ztopcp(\|^[[:space:]]*real(8)[[:space:]]*zbotcp(\|^[[:space:]]*integer[[:space:]]*layer(" src/core/variables.f90`.

- [ ] **Step 2: Delete each declaration line.** Also delete adjacent comment lines describing each.

- [ ] **Step 3: Drop zero-fills in `initialize.f90`.** Run: `grep -n "numnod\|dz(\|z(\|disnod\|ztopcp\|zbotcp\|layer(" src/core/initialize.f90` → delete or update each to write to `state%mesh%X` (decide per case; usually just delete since `state%mesh%init` covers it).

- [ ] **Step 4: Run VG.** If compile fails, the error line surfaces a missed reader — patch in the offending file, re-run.

- [ ] **Step 5: Commit.**

```bash
git add src/core/variables.f90 src/core/initialize.f90
git commit -m "$(cat <<'EOF'
retire(gr-bh): delete mesh globals from variables.f90

numnod/dz(:)/z(:)/disnod(:)/ztopcp(:)/zbotcp(:)/layer(:) declarations
deleted. Zero-fills in initialize.f90 dropped. state%mesh is the sole
home for vertical-discretization data.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 36: Delete layer-flat + runtime globals from `variables.f90`

**Files:** `src/core/variables.f90`

- [ ] **Step 1: Find declarations.** `grep -n "^[[:space:]]*\(real(8)\|integer\|logical\)[[:space:]]*\(ksatexm\b\|ksatfit\b\|cofani\b\|flksatexm\b\|orgmat\b\|psand\b\|psilt\b\|pclay\b\|q0\b\|k1max\b\|H0max\b\|ArMpSs\b\|swbotb\b\)" src/core/variables.f90`.

- [ ] **Step 2: Delete declarations.**

- [ ] **Step 3: Drop zero-fills in `initialize.f90`.**

- [ ] **Step 4: Run VG.**

- [ ] **Step 5: Commit.**

```bash
git add src/core/variables.f90 src/core/initialize.f90
git commit -m "$(cat <<'EOF'
retire(gr-bh): delete layer-flat + runtime globals from variables.f90

ksatexm/ksatfit/cofani/flksatexm/orgmat/psand/psilt/pclay (8 layer flats)
+ q0/k1max/H0max/ArMpSs/swbotb (5 runtime scalars; ArMpSs deletion
completes macropore retirement) deleted. state%soilwater + state mirrors
are the sole home.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 37: Delete drainage-geo globals from `variables.f90`

**Files:** `src/core/variables.f90`

- [ ] **Step 1: Find declarations.** `grep -n "^[[:space:]]*\(real(8)\|integer\)[[:space:]]*\(nrlevs\b\|swdivd\b\|swnrsrf\b\|swtopnrsrf\b\|swdivdinf\b\|FacDpthInf\b\|L(\|zbotdr(\|owltab(\)" src/core/variables.f90`.

- [ ] **Step 2: Delete declarations.**

- [ ] **Step 3: Drop zero-fills in `initialize.f90`.**

- [ ] **Step 4: Run VG.**

- [ ] **Step 5: Commit.**

```bash
git add src/core/variables.f90 src/core/initialize.f90
git commit -m "$(cat <<'EOF'
retire(gr-bh): delete drainage-geo globals from variables.f90

nrlevs/swdivd/swnrsrf/swtopnrsrf/swdivdinf/FacDpthInf/L(:)/zbotdr(:)/
owltab(:) declarations deleted. state%drainage is the sole home for
drainage geometry / switches.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

### Task 38: Sanity-grep retirement

**Files:** None (verification only).

- [ ] **Step 1: Mesh-globals declarations gone.**

```bash
grep -nw "numnod\|dz\|z\|disnod\|ztopcp\|zbotcp\|layer" src/core/variables.f90
```

Expected: empty (or only comment hits).

- [ ] **Step 2: Layer-flat globals gone.**

```bash
grep -nw "ksatexm\|ksatfit\|cofani\|flksatexm\|orgmat\|psand\|psilt\|pclay" src/core/variables.f90
```

Expected: empty.

- [ ] **Step 3: Runtime + drainage-geo globals gone.**

```bash
grep -nw "q0\|k1max\|H0max\|ArMpSs\|swbotb\|nrlevs\|swdivd\|swnrsrf\|swtopnrsrf\|swdivdinf\|FacDpthInf\|owltab" src/core/variables.f90
```

Expected: empty.

- [ ] **Step 4: State-mesh usage spread.** `grep -rl "state%mesh%" src/ | wc -l` ≥ 20.

- [ ] **Step 5: No commit.**

---

### Task 39: Final verification + arc-complete marker

**Files:** None.

- [ ] **Step 1: Clean rebuild + check-full.** `rm -rf builddir && pixi run build-linux && pixi run check-full` — **5/5 byte-for-byte**.

- [ ] **Step 2: BMI + cffi suites.** `pixi run -e test test-bmi && pixi run -e test test-cffi-demo` — both pass.

- [ ] **Step 3: pFUnit.** `pixi run test-pfunit` — all-pass.

- [ ] **Step 4: Final sanity greps.**
  - `grep -n "^[[:space:]]*use variables" src/heat/ src/boundary/ -r | grep -v "only: *nird\|only: *logf"` empty.
  - `grep -rl "state%mesh%" src/ | wc -l` ≥ 25.
  - `grep -rl "state%soilwater%ksatexm\|state%soilwater%orgmat" src/ | wc -l` ≥ 5.
  - `grep -rl "state%drainage%nrlevs\|state%drainage%zbotdr" src/ | wc -l` ≥ 5.

- [ ] **Step 5: Arc-complete marker.**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(gr-bh): GR-BH complete — boundary+heat globals retired, mesh extracted

check-full 5/5 byte-for-byte. BMI + cffi-demo passing. pFUnit all-pass.

Net retirement:
- 4 readers (boundbottom/boundtop/frozencond/temperature) off use variables
  (boundtop retains only: nird → Arc 8; boundbottom retains only: logf
  → Arc 9)
- heat_init relocated to state%heat%init (type-bound, surfacewater_state
  pilot pattern)
- 7 mesh globals retired → state%mesh
- 8 soilwater layer flats + 4 runtime scalars retired → state%soilwater
- 9 drainage geo/switch fields retired → state%drainage
- ArMpSs deleted (macropore retirement completion)
- ~30 codebase readers swept to state references

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Plan Self-Review

**Spec coverage:**
- ✅ Schema additions (mesh_state, soilwater extensions, drainage extensions, heat type-bound init) → Tasks 2, 4, 5, 8, 11
- ✅ Dual-write populator → Tasks 3, 6, 7, 9
- ✅ Phase A close → Task 10
- ✅ Heat reader migration → Tasks 12, 13, 15, 16
- ✅ Boundary reader migration → Tasks 18, 20, 21
- ✅ Verification tasks → Tasks 14, 17, 19, 22
- ✅ Phase B close → Task 23
- ✅ Codebase mesh sweep clusters → Tasks 24–30
- ✅ Audit pass → Task 31
- ✅ Populator legacy-write retirement → Tasks 32, 33, 34
- ✅ Global declaration deletion → Tasks 35, 36, 37
- ✅ Phase C close + arc-complete → Tasks 38, 39

**Placeholder scan:** No "TBD" / "implement later" / "similar to Task N". Code blocks present in every code step. Audit checklist concrete with exact greps. Caveats explicit per task where field-naming verification is required by the implementer.

**Type consistency:** `mesh_state_t` field names consistent across spec + Task 2 + Task 3 dual-write + Task 32 retirement source. `swbotb_runtime` naming consistent. `state%mesh%init` signature consistent.

**Known soft spots** (flagged inline, not bugs):
- Several config field names (`config%heat%tmean`, `config%simulation%flrunon`, etc.) are noted as "verify via grep before substitution" — implementer subagent confirms in each task. Names may differ slightly from spec table.
- `Tav`/`atav` placement (state%atmosphere vs variables) is a verify-before-migrate caveat in Task 15. If absent on atmosphere_state, retained via narrow `use variables, only:` deferral.
- `swdra` placement: Task 8 notes it lives on state%surfacewater per GR-UTILS precedent; drainage_state does NOT duplicate.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-05-13-globals-boundary-heat.md`. Two execution options:

1. **Subagent-Driven (recommended)** — Fresh subagent per task, two-stage review (spec compliance + code quality) between each. Fast iteration, isolated context per task. Standard pattern for this arc series.
2. **Inline Execution** — Execute tasks in this session via executing-plans, batch execution with checkpoints. Higher main-context pressure.

Which approach?
