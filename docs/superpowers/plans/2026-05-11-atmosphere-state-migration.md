# Atmosphere State-Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Carve ~33 atmosphere-owned globals into `state%atmosphere` (NEW state record); introduce two cohort sub-records (`atmosphere_intermediate_t`, `atmosphere_cumulative_t`) under `flzerointr`/`flzerocumu` gates; preserve byte-identical check-full at every commit.

**Architecture:** Flat top-level scalars (22 fields) + 2 cohort sub-records (8 intermediate + 10 cumulative = 18 cohort fields). Each cohort has a type-bound `reset()` procedure replacing 3 scattered reset blocks. `atmosphere_init(state%atmosphere)` wires at `swap.f90:188` after `soilwater_init`. All 7 atmosphere files in scope — including `meteoday.f90` and `meteodt.f90` (scope correction from discovery exclusion plan); the future atmosphere REFACTOR arc handles orchestration restructuring. ASSOCIATE prefix `at_` used in large compute bodies.

**Tech Stack:** Fortran 2008+ derived types with type-bound procedures, ASSOCIATE, pixi+meson, pFUnit, check-full regression.

**Spec:** `docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md`
**Discovery:** `docs/superpowers/specs/2026-05-11-state-migration-atmosphere-discovery.md`
**Branch:** all commits go on `refactor/atmosphere-state` (already created from `development`).

---

## Lessons-learned applied

- **Dual-write before cutover.** Every global write gets a paired `state%atmosphere%X = …` mirror before readers are cut over. Drop globals only in Phase 2 after all readers are migrated.
- **Verify-state-not-reports.** Run `check-full` output comparison, not just pFUnit, before every commit. pFUnit alone misses global-default regressions (see `feedback_verify_before_committing.md`).
- **ASSOCIATE prefix `at_`.** Shadow-safe aliasing in heavy compute bodies. Follows the `bw_` (boundary), `ht_` (heat), `dr_` (drainage), `sl_` (solute), `cw_` (crop-water) convention.
- **Grep field shapes during pre-flight.** Before deleting a global, run both `use variables.*\b<var>\b` AND raw-symbol greps. Both must return zero hits.
- **Compile-driven Phase 2.7 is expected.** After dropping dual-writes, expect 3–6 hidden readers. Plan for iterative fixup commits.
- **Cohort sub-record paths are one level deeper.** `state%atmosphere%intr%igrai`, `state%atmosphere%cumu%cgrai`. ASSOCIATE blocks at call sites keep local names short.
- **Cumulatives have a single canonical reset owner.** The double-reset of `igrai`/`inrai` (waterbalance + meteoday) collapses to `state%atmosphere%intr%reset()` at one site. Do not partially remove duplicate resets; coordinate A-1.8/A-2.2 to remove both together.

---

## Cohort field allocation

### Flat top-level (22 scalars)

**Instantaneous (11):** `peva`, `ptra`, `empreva`, `melt`, `subl`, `slw`, `ssnow`, `snowinco`, `graidt`, `nraidt`, `aintcdt`

**Per-day (9):** `grai`, `nraida`, `atmdem`, `pevaday`, `ptraday`, `gsnow`, `snrai`, `fprecnosnow`, `sicact`

**Per-event / tillage-reset (3) [also flat top-level]:** `ldwet`, `spev`, `saev`

### `atmosphere_intermediate_t` (8 fields, `flzerointr`-reset)

`igrai`, `inrai`, `ipeva`, `iptra`, `ievap`, `igsnow`, `isubl`, `isnrai`

### `atmosphere_cumulative_t` (10 fields, `flzerocumu`-reset)

`cgrai`, `cnrai`, `caintc`, `cpeva`, `cptra`, `cevap`, `cgsnow`, `csubl`, `csnrai`, `cmelt`

**Total: 40 fields in one record (22 flat + 8 intr + 10 cumu).**

---

## File structure

**New:**
- `src/state/atmosphere_state.f90` — `atmosphere_state_t`, `atmosphere_intermediate_t`, `atmosphere_cumulative_t`, `atmosphere_init` (A-1.1)
- `tests/unit/state/test_atmosphere_state.pf` — pFUnit lifecycle tests (A-1.1)
- `docs/adr/0037-state-migration-atmosphere.md` — ADR (A-2.7)

**Modified:**
- `src/state/swap_state.f90` — add `atmosphere_state_t` component (A-1.2)
- `src/core/swap.f90` — wire `atmosphere_init` at line 188; snow.f90 callers at 216/298 (A-1.2, A-1.3)
- `src/atmosphere/snow.f90` — signature promotion + dual-write (A-1.3)
- `src/atmosphere/et.f90` — dual-write empreva/ldwet/spev/saev/peva via `reduceva`; new `state` arg (A-1.4)
- `src/atmosphere/interception.f90` — dual-write nraida/sicact; new `state` arg to DivIntercep (A-1.5)
- `src/atmosphere/precipitation.f90` — dual-write grai/gsnow/snrai/fprecnosnow + ssnow mutation; new `state` arg (A-1.6)
- `src/atmosphere/snow.f90` — dual-write ssnow/melt/subl/slw/snowinco + cohort fields (covered by A-1.3)
- `src/atmosphere/meteoday.f90` — dual-write peva/ptra/atmdem/pevaday/ptraday/grai; cohort reset dual-patch (A-1.7)
- `src/atmosphere/meteodt.f90` — dual-write peva/ptra (A-1.8)
- `src/soil/soilhydraulics.f90` — cohort reset migration; melt/ldwet/spev/saev readers cutover (A-1.9, A-2.1)
- `src/soil/waterbalance.f90` — igrai/inrai/ipeva/iptra accumulator paths + duplicate-reset removal (A-2.2)
- `src/crop/rootextraction.f90` — ptra/atmdem readers cutover (A-2.3)
- `src/crop/cropgrowth.f90` — ptra readers cutover (A-2.3)
- `src/crop/tillage.f90` — nraida readers cutover (A-2.4)
- `src/boundary/boundtop.f90` — peva/empreva/melt readers cutover (A-2.4)
- `src/heat/temperature.f90` — ssnow readers cutover (A-2.4)
- `src/solute/solute.f90` + `src/solute/agetracer.f90` — nird-adjacent reads (verify only) (A-2.4)
- `src/io/swapoutput.f90` — full cumulative + intermediate + ssnow/snowinco/sicact output cutover (A-2.5)
- `src/io/swap_csv_output.f90` — EPOT/TPOT/SSNOW + dstor cutover (A-2.5)
- `src/core/variables.f90` — retire ~33 globals with `[SS-ATM]` markers (A-2.6)
- `src/core/initialize.f90` — drop redundant zero-inits for retired globals (A-2.6)
- `docs/superpowers/specs/state-migration-playbook.md` — append atmosphere lessons (A-2.7)

---

# Phase 0 — Config gaps

### Task A-0.1: Audit `swsublim` typed-config coverage + `spev`/`saev` adapter seeding

**Design decision:** D12
**Files:** read-only audit; close gaps if found.

`swsublim` is used at `snow.f90:104` to gate the sublimation path but its config-coverage status was not confirmed in the discovery. `spev` and `saev` get initial values (reset on tillage) from soil config. Both are Phase 0 candidates (discovery Section 6).

- [ ] **Grep `swsublim` in config:**
  ```bash
  grep -rn "swsublim" src/config/ src/io/toml/ --include="*.f90"
  ```
  If absent from `meteorology_config.snow` or the TOML adapter, add default and seeding (0 = no sublimation). If present, confirm adapter seeds the variable.
- [ ] **Grep `spev`/`saev` initial seeding:**
  ```bash
  grep -rn "\bspev\b\|\bsaev\b" src/config/ src/io/toml/ --include="*.f90"
  ```
  Both should be seeded from `soil_config.initial` fields (confirmed in discovery Section 2a). If absent, add config fields with `0.0` defaults. If present, document audit result.
- [ ] **If gaps found:** close each with a minimal config addition + adapter seed line. Each gap = one commit.
- [ ] **If no gaps:** one documentation-only commit records the audit.

```bash
git commit -m "$(cat <<'EOF'
docs(config): SS-ATM Phase 0 A-0.1 — audit swsublim + spev/saev config coverage

Audit of meteorology_config.snow (swsublim) and soil_config.initial
(spev/saev) adapter seeding. [FILL: found N gaps / found 0 gaps].
Any gaps closed; no-gap result documented.

Spec: docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md (ADR 0037)
EOF
)"
```

---

# Phase 1 — State type + init + dual-write

### Task A-1.1: Create `src/state/atmosphere_state.f90`; pFUnit lifecycle tests

**Design decisions:** D2, D3, D4
**Files:**
- Create: `src/state/atmosphere_state.f90`
- Create: `tests/unit/state/test_atmosphere_state.pf`

Define `atmosphere_intermediate_t` (8 fields + `reset()` type-bound procedure), `atmosphere_cumulative_t` (10 fields + `reset()` type-bound procedure), and `atmosphere_state_t` (22 flat scalars + `type(atmosphere_intermediate_t) :: intr` + `type(atmosphere_cumulative_t) :: cumu`). Add `atmosphere_init(atm)` subroutine that zeros all 22 flat scalars (cohort fields are already zero-defaulted at declaration). No `numnod`/`nlay` parameters — all fields are scalars.

Field layout to implement:

```fortran
type :: atmosphere_intermediate_t
   real(real64) :: igrai  = 0.0_real64
   real(real64) :: inrai  = 0.0_real64
   real(real64) :: ipeva  = 0.0_real64
   real(real64) :: iptra  = 0.0_real64
   real(real64) :: ievap  = 0.0_real64
   real(real64) :: igsnow = 0.0_real64
   real(real64) :: isubl  = 0.0_real64
   real(real64) :: isnrai = 0.0_real64
contains
   procedure :: reset => atmosphere_intermediate_reset
end type

type :: atmosphere_cumulative_t
   real(real64) :: cgrai  = 0.0_real64
   real(real64) :: cnrai  = 0.0_real64
   real(real64) :: caintc = 0.0_real64
   real(real64) :: cpeva  = 0.0_real64
   real(real64) :: cptra  = 0.0_real64
   real(real64) :: cevap  = 0.0_real64
   real(real64) :: cgsnow = 0.0_real64
   real(real64) :: csubl  = 0.0_real64
   real(real64) :: csnrai = 0.0_real64
   real(real64) :: cmelt  = 0.0_real64
contains
   procedure :: reset => atmosphere_cumulative_reset
end type
```

- [ ] **Write `atmosphere_state.f90`** with both cohort types, their `reset()` implementations, `atmosphere_state_t`, and `atmosphere_init`. All fields zero-defaulted. Module: `atmosphere_state_mod`. Public exports: `atmosphere_state_t`, `atmosphere_intermediate_t`, `atmosphere_cumulative_t`, `atmosphere_init`.
- [ ] **Write `test_atmosphere_state.pf`:**
  - `@test` all 22 flat scalars default to `0.0_real64` on a fresh `atmosphere_state_t()`.
  - `@test` `intr%reset()` zeroes all 8 fields (seed non-zero values first, then call, then assert zero).
  - `@test` `cumu%reset()` zeroes all 10 fields (same pattern).
  - `@test` two independent `atmosphere_state_t` instances share no state.
- [ ] **Add to meson build** — add `atmosphere_state.f90` to the state library sources; add the test file to pFUnit test sources.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
feat(state): SS-ATM Phase 1 A-1.1 — atmosphere_state_t with two cohorts

New atmosphere_state_t: 22 flat scalars (11 instantaneous + 9 per-day +
3 per-event) + atmosphere_intermediate_t (8 fields, flzerointr-reset) +
atmosphere_cumulative_t (10 fields, flzerocumu-reset). Each cohort has
a type-bound reset() procedure. atmosphere_init zeros all scalars.
First arc to debut the ADR 0033 cohort pattern at type-creation.
pFUnit 4 new tests green. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md (ADR 0037)
EOF
)"
```

---

### Task A-1.2: Aggregate `state%atmosphere` in `swap_state.f90`; wire `atmosphere_init` at `swap.f90:188`

**Design decisions:** D3, D4
**Files:**
- Modify: `src/state/swap_state.f90`
- Modify: `src/core/swap.f90`

- [ ] **Add to `swap_state.f90`:**
  - `use atmosphere_state_mod, only: atmosphere_state_t`
  - `type(atmosphere_state_t) :: atmosphere` component in `swap_state_t`.
- [ ] **Add `atmosphere_init` call in `swap.f90`** at line 188 (immediately after the existing `call soilwater_init(state%soilwater, numnod, numlay)`). Add `use atmosphere_state_mod, only: atmosphere_init` to `swap.f90`'s use block if not already imported via `swap_state_mod`.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full` (no compute paths changed yet).

```bash
git commit -m "$(cat <<'EOF'
feat(state): SS-ATM Phase 1 A-1.2 — aggregate state%atmosphere; wire atmosphere_init

swap_state_t gains type(atmosphere_state_t) :: atmosphere. atmosphere_init
called at swap.f90:188 after soilwater_init. No fields yet assigned — init
establishes zero baseline. pFUnit green. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md (ADR 0037)
EOF
)"
```

---

### Task A-1.3: `snow.f90` signature promotion + dual-write

**Design decisions:** D9, D1
**Files:**
- Modify: `src/atmosphere/snow.f90`
- Modify: `src/core/swap.f90` (callers at lines 216, 298 — verify no change needed since both already pass full `state`)

Discovery hazard #7: `snow.f90:106–109` reads `peva`, writes `subl = peva`, then zeros `peva` and `empreva`. The current `optional intent(in)` prevents writes. Promotion to `intent(inout)` enables the state write-backs.

Dual-write targets in `snow.f90`:
- `ssnow` (lines 74, 76, 122, 141, 156, 165)
- `snowinco` (lines 74, 99)
- `slw` (lines 141, 144, 153, 156, 166)
- `melt` (lines 123, 138, 157, 163)
- `subl` (lines 103, 106, 124, 164)
- `peva` and `empreva` zero-writes (lines 108–109) → `state%atmosphere%peva = 0.0_real64`; `state%atmosphere%empreva = 0.0_real64`
- `gsnow` (lines 121, 123, 141, 171) — read only in snow; written by precipitation.f90 (covered in A-1.6)
- `snrai` (lines 131, 132, 144, 173, 177) — read only in snow; written by precipitation.f90 (covered in A-1.6)
- Cohort fields written by task=1 (init): `snowinco = ssnow` (line 74) → no intermediate/cumulative touches at init
- Cohort fields written by task=2 (per-day): `igsnow` (lines 87, 171), `isubl` (lines 88, 172), `isnrai` (lines 89, 173) → `state%atmosphere%intr%igsnow = igsnow` etc.; `cgsnow` (lines 95, 174), `csubl` (lines 96, 175), `cmelt` (lines 98, 176), `csnrai` (lines 97, 177) → `state%atmosphere%cumu%cgsnow = cgsnow` etc.
- **Inline reset block at snow.f90:87–98 (flzerointr/flzerocumu gates):** add dual-patch here: after `igsnow = 0.0` add `state%atmosphere%intr%igsnow = 0.0_real64`; same for the other 3 intr fields and 4 cumu fields. The inline resets remain until Phase 2 fully removes them (A-2.1).

- [ ] **Promote `state` arg** in `snow` from `optional intent(in)` to `intent(inout)`. Update the dummy arg declaration. Verify `swap.f90:216` and `swap.f90:298` already pass full `state` (they do per discovery Section 1 entry-point table).
- [ ] **Add ASSOCIATE block** at top of `snow` task=2 body: `associate(at => state%atmosphere)`.
- [ ] **Install dual-writes** for all write sites listed above. Use `at%ssnow = ssnow` form inside the ASSOCIATE.
- [ ] **Dual-patch the inline reset block** at lines 87–98 (both `flzerointr` and `flzerocumu` blocks).
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-ATM Phase 1 A-1.3 — snow.f90 signature promotion + dual-write

state arg promoted from optional intent(in) to intent(inout). Dual-writes
installed for ssnow, snowinco, slw, melt, subl, peva/empreva zero-paths,
igsnow/isubl/isnrai (intermediate cohort), cgsnow/csubl/csnrai/cmelt
(cumulative cohort). Inline reset block dual-patched. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md (ADR 0037)
EOF
)"
```

---

### Task A-1.4: `et.f90` dual-write — `reduceva` state arg + 5 owned fields

**Design decisions:** D7, D8
**Files:**
- Modify: `src/atmosphere/et.f90`

`reduceva` is the only subroutine in `et.f90` that touches global state. It reads/writes `empreva`, `ldwet`, `spev`, `saev`, and reads `peva`. Discovery Section 2a write sites:
- `empreva`: lines 485, 515, 543, 569, 573, 638 (via `black_reduction`/`boesten_stroosnijder_reduction`)
- `ldwet`: lines 497, 498, 505, 511, 639
- `spev`: lines 556, 581, 583, 640
- `saev`: lines 558, 562, 565, 576, 641
- `peva`: read at line 637 (the `if (pond > POND_THRESHOLD)` path — pond is deferred)
- `snow.f90:108` zeroes `empreva` via legacy global (already handled in A-1.3)

- [ ] **Add `state` arg to `reduceva`:** `subroutine reduceva(task, nrai, state)` with `type(swap_state_t), intent(inout) :: state`. Update callers: `meteoday.f90:808`, `meteodt.f90:358`, `meteodt.f90:448` (additive — pass state).
- [ ] **Add `use swap_state_mod` to `et.f90`** if not already present.
- [ ] **Add ASSOCIATE** at top of `reduceva`: `associate(at => state%atmosphere)`.
- [ ] **Dual-write `empreva`, `ldwet`, `spev`, `saev`** at every write site inside `reduceva`, `black_reduction`, and `boesten_stroosnijder_reduction`. Form: `empreva = X; at%empreva = X` (or `at%empreva = empreva` immediately after).
- [ ] **Add comment marker at `pond` read** (line 637): `! [SS-ATM] reads legacy pond — soil-water-core arc migrates`.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-ATM Phase 1 A-1.4 — et.f90 reduceva dual-write (5 atmosphere fields)

reduceva gains state intent(inout) arg. Dual-writes for empreva, ldwet,
spev, saev (all write sites in reduceva/black_reduction/boesten). peva
read unchanged (legacy global). pond read marked [SS-ATM] deferral.
Callers meteoday:808, meteodt:358+448 updated. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md (ADR 0037)
EOF
)"
```

---

### Task A-1.5: `interception.f90` dual-write — `nraida`, `sicact`

**Design decisions:** D1
**Files:**
- Modify: `src/atmosphere/interception.f90`

`DivIntercep` writes `nraida` (lines 403, 407, 410) and `nird` (lines 404, 408, 411 — co-write, NOT migrated per D10). `ruttervw` writes `sicact` via `msw1eic` return (line 189). `eintc` and `aintc` are local/returned — not in variables.f90, not migrated.

- [ ] **Add `state` arg to `DivIntercep`:** `subroutine DivIntercep(aintc, state)` with `type(swap_state_t), intent(inout) :: state`. Update callers: `meteoday.f90:530`, `meteoday.f90:653`.
- [ ] **Dual-write `nraida`** at lines 403, 407, 410: `nraida = X; state%atmosphere%nraida = X`.
- [ ] **Add `state` arg to the `ruttervw`/`msw1eic` call chain** (or directly to the `ruttervw` call at `meteoday.f90:649`) so `sicact` can be dual-written. `sicact` is written via `msw1eic` return at `interception.f90:189`. Dual-write at that return point: `sicact = msw1eic(...); state%atmosphere%sicact = sicact`.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-ATM Phase 1 A-1.5 — interception.f90 dual-write (nraida, sicact)

DivIntercep gains state intent(inout) arg. Dual-writes for nraida (3 sites)
and sicact (ruttervw return path, 1 site). nird co-write retained as legacy
(irrigation-owned, deferred per D10). Callers meteoday:530+653 updated.
check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md (ADR 0037)
EOF
)"
```

---

### Task A-1.6: `precipitation.f90` dual-write — `grai`, `gsnow`, `snrai`, `fprecnosnow` + `ssnow` mutation

**Design decisions:** D1, D11
**Files:**
- Modify: `src/atmosphere/precipitation.f90`

Write sites (discovery Section 2a):
- `grai`: lines 73, 117, 119
- `gsnow`: lines 79, 83, 88, 108, 126
- `snrai`: lines 93, 95, 109, 128
- `fprecnosnow`: lines 101, 103, 110, 129
- `ssnow = 0.0d0` mutation at line 127 (swmetdetail=1 path; keep verbatim per D11)

- [ ] **Add `state` arg to `PartitionPrecipitation`:** `subroutine PartitionPrecipitation(..., state)` with `type(swap_state_t), intent(inout) :: state`. Update caller: `meteoday.f90:369`.
- [ ] **Dual-write all 4+1 fields** at every write site. `ssnow` mutation at line 127: add `state%atmosphere%ssnow = 0.0_real64` immediately after `ssnow = 0.0d0` (D11 — keep legacy mutation verbatim).
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-ATM Phase 1 A-1.6 — precipitation.f90 dual-write

PartitionPrecipitation gains state intent(inout) arg. Dual-writes for
grai (3 sites), gsnow (5), snrai (4), fprecnosnow (4). ssnow=0 mutation
at line 127 (swmetdetail=1) dual-patched verbatim per D11. Caller
meteoday:369 updated. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md (ADR 0037)
EOF
)"
```

---

### Task A-1.7: `meteoday.f90` dual-write — `peva`, `ptra`, `atmdem`, day fields, cohort resets

**Design decisions:** D1
**Files:**
- Modify: `src/atmosphere/meteoday.f90`

This is the largest single dual-write task. `meteoday.f90` (884 LoC) is now a full home routine, not an excluded file. Write sites from discovery (authoritative inventory):

**`peva` writes:** lines 706, 708, 713, 715, 722, 725, 727, 735, 737

**`ptra` writes:** lines 743, 745, 747, 751, 752, 757, 864

**`atmdem` writes:** lines 816, 857, 860

**`pevaday`/`ptraday`** writes (ETSine path): in `ProcessMeteoDay` / `ETSine` helper — grep for `pevaday =` and `ptraday =` before starting.

**`grai` write:** line 295 (load from `arai`, mm→cm convert)

**`aintcdt` writes:** grep for `aintcdt =` in meteoday.

**Cohort reset blocks (dual-patch):**
- `igrai = 0.0d0` at line 414 → add `state%atmosphere%intr%igrai = 0.0_real64`
- `inrai = 0.0d0` at line 415 → add `state%atmosphere%intr%inrai = 0.0_real64`
- `cgrai = 0.0d0` at line 420 → add `state%atmosphere%cumu%cgrai = 0.0_real64`
- `cnrai = 0.0d0` at line 421 → add `state%atmosphere%cumu%cnrai = 0.0_real64`
- `caintc = 0.0d0` at line 422 → add `state%atmosphere%cumu%caintc = 0.0_real64`

NOTE: the full cohort `reset()` calls do NOT replace the inline resets yet — that is Phase 2 (A-2.2). The dual-patch here keeps the legacy zeros AND adds state writes. This preserves byte-identical behaviour while establishing the state path.

- [ ] **Add `state` arg** to `ProcessMeteoDay` (and `ReadMeteoDay` if it writes any owned globals — check). Signature: `subroutine ProcessMeteoDay(state)` with `type(swap_state_t), intent(inout) :: state`. Update caller at `swap.f90:285`.
- [ ] **Add `state` threading to `ResetMetFlx`** (the reset subroutine) if it is a separate subroutine — check. Pass state down to it.
- [ ] **Add ASSOCIATE block** at top of `ProcessMeteoDay`: `associate(at => state%atmosphere)`.
- [ ] **Dual-write all `peva` write sites** (9 sites listed above): `peva = X; at%peva = X`.
- [ ] **Dual-write all `ptra` write sites** (7 sites): same pattern.
- [ ] **Dual-write all `atmdem` write sites** (3 sites): same pattern.
- [ ] **Dual-write `pevaday`/`ptraday`** (ETSine path — grep to confirm line numbers).
- [ ] **Dual-write `grai`** at line 295.
- [ ] **Dual-write `aintcdt`** (grep to confirm line numbers).
- [ ] **Dual-patch the cohort reset block** at lines 412–423 (5 fields as listed above).
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-ATM Phase 1 A-1.7 — meteoday.f90 dual-write

ProcessMeteoDay gains state intent(inout). Dual-writes for peva (9 sites),
ptra (7 sites), atmdem (3 sites), pevaday/ptraday (ETSine path), grai
(line 295), aintcdt. Cohort reset block (lines 412-423) dual-patched:
igrai/inrai/cgrai/cnrai/caintc now written to state intr/cumu alongside
legacy zeros. Inline resets remain (removed in A-2.2). check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md (ADR 0037)
EOF
)"
```

---

### Task A-1.8: `meteodt.f90` dual-write — `peva`, `ptra`

**Design decisions:** D1
**Files:**
- Modify: `src/atmosphere/meteodt.f90`

`meteodt.f90` (453 LoC) writes `peva` at lines 349 and 444, and `ptra` at lines 348 and 445. Also calls `reduceva` (lines 358, 448) which already takes `state` after A-1.4. The `ETSine` path writes `pevaday`/`ptraday` — check whether these are in `meteodt.f90` or `meteoday.f90` (discovery indicates meteodt:444–445 write peva/ptra under ETSine; confirm exact pevaday/ptraday owners).

- [ ] **Add `state` arg to `MeteoDT`:** `subroutine MeteoDT(state)` with `type(swap_state_t), intent(inout) :: state`. Update caller at `swap.f90:291`.
- [ ] **Thread state into `ETSine`** helper if it writes pevaday/ptraday (check meteodt.f90 content).
- [ ] **Add ASSOCIATE** at top of `MeteoDT`: `associate(at => state%atmosphere)`.
- [ ] **Dual-write `peva`** at lines 349, 444: `peva = X; at%peva = X`.
- [ ] **Dual-write `ptra`** at lines 348, 445: `ptra = X; at%ptra = X`.
- [ ] **Dual-write `pevaday`/`ptraday`** if written here (confirm via grep).
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-ATM Phase 1 A-1.8 — meteodt.f90 dual-write (peva, ptra)

MeteoDT gains state intent(inout). Dual-writes for peva (lines 349+444)
and ptra (lines 348+445). pevaday/ptraday dual-written in ETSine path
if owned here. reduceva callers at 358+448 already pass state (A-1.4).
swap.f90:291 caller updated. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md (ADR 0037)
EOF
)"
```

---

### Task A-1.9: `soilhydraulics.f90` cohort reset dual-patch + `ldwet`/`spev`/`saev`/`nird` reset sites

**Design decisions:** D5
**Files:**
- Modify: `src/soil/soilhydraulics.f90`

`soilhydraulics.f90` is both a co-writer (resets cpeva/cptra/cevap at 1140–1142; ipeva/iptra/ievap at 1112–1114; ldwet/spev/saev at 873–875; nird=0 at 871) and a future reader (melt reads deferred to A-2.1). This task dual-patches the reset sites; the full cohort `reset()` calls and legacy-reset removal happen in Phase 2 (A-2.1).

Dual-patch pattern for each reset field:
- `cpeva = 0.0d0` at line 1140 → add `state%atmosphere%cumu%cpeva = 0.0_real64`
- `cptra = 0.0d0` at line 1141 → add `state%atmosphere%cumu%cptra = 0.0_real64`
- `cevap = 0.0d0` at line 1142 → add `state%atmosphere%cumu%cevap = 0.0_real64`
- `ipeva = 0.0d0` at line 1113 → add `state%atmosphere%intr%ipeva = 0.0_real64`
- `iptra = 0.0d0` at line 1112 → add `state%atmosphere%intr%iptra = 0.0_real64`
- `ievap = 0.0d0` at line 1114 → add `state%atmosphere%intr%ievap = 0.0_real64`
- `ldwet` reset at line 873: `ldwet = X; state%atmosphere%ldwet = X`
- `spev` reset at line 874: `spev = 0.0d0; state%atmosphere%spev = 0.0_real64`
- `saev` reset at line 875: `saev = 0.0d0; state%atmosphere%saev = 0.0_real64`

- [ ] **Confirm `state` is already in scope** in the soilhydraulics routines containing these lines (post-ADR-0036 — should be true for SoilWater(1) body and the tillage-reset block).
- [ ] **Add ASSOCIATE if needed**: `associate(at => state%atmosphere)`.
- [ ] **Dual-patch all 9 reset sites** as listed above. Do NOT remove the legacy resets yet.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-ATM Phase 1 A-1.9 — soilhydraulics.f90 cohort reset dual-patch

Dual-patch for 9 reset sites: ipeva/iptra/ievap (intermediate cohort,
lines 1112-1114), cpeva/cptra/cevap (cumulative cohort, lines 1140-1142),
ldwet/spev/saev (tillage reset, lines 873-875). Legacy resets remain;
state%atmosphere now mirrors them. Full cohort reset() and removal of
legacy inline resets happen in A-2.1. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md (ADR 0037)
EOF
)"
```

---

# Phase 2 — Reader cutover + retire

### Task A-2.1: `soilhydraulics.f90` reader cutover — `melt` reads + cohort reset migration

**Design decisions:** D5, D1
**Files:**
- Modify: `src/soil/soilhydraulics.f90`

`melt` is read at lines 109, 639, 649, 651 in the top-flux equation (net inflow includes melt). Migrate these 4 reads from the legacy global to `state%atmosphere%melt`. Also: replace the 6 legacy-inline-reset lines (dual-patched in A-1.9) with clean `call state%atmosphere%intr%reset()` and `call state%atmosphere%cumu%reset()` invocations under their respective flags. Remove the legacy zero assignments.

- [ ] **Migrate 4 `melt` read sites** — `melt` → `state%atmosphere%melt`. Drop `melt` from `use variables, only:` (or from the `use variables` clause) in the relevant routines.
- [ ] **Replace cpeva/cptra/cevap reset block (lines 1140–1142)** with `call state%atmosphere%cumu%reset()` under the existing `flzerocumu` gate. The `cumu%reset()` zeroes all 10 cumulative fields including the 3 previously set here.
- [ ] **Replace ipeva/iptra/ievap reset block (lines 1112–1114)** with `call state%atmosphere%intr%reset()` under the existing `flzerointr` gate. The `intr%reset()` zeroes all 8 intermediate fields.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-ATM Phase 2 A-2.1 — soilhydraulics: melt reads + cohort reset()

4 melt read sites migrated to state%atmosphere%melt. Legacy inline reset
blocks for ipeva/iptra/ievap and cpeva/cptra/cevap replaced with type-bound
reset() calls (intr%reset() and cumu%reset()) — 6 lines of scattered zeros
collapse to 2 cohort resets. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md (ADR 0037)
EOF
)"
```

---

### Task A-2.2: `waterbalance.f90` reader cutover + duplicate-reset removal (`igrai`/`inrai`)

**Design decisions:** D6, D1
**Files:**
- Modify: `src/soil/waterbalance.f90`

`waterbalance.f90` accumulates atmosphere fields (discovery Section 3.3: ~14 sites). It also has the double-reset of `igrai`/`inrai` at lines 388–389. After A-1.7 dual-patched meteoday's reset block, and A-2.1 landed the canonical `intr%reset()` in soilhydraulics, the duplicate zero-writes at waterbalance:388–389 can be safely removed.

Reads to migrate:
- `ipeva` accumulator path (waterbalance.f90:470): `ipeva = ipeva + X` → `state%atmosphere%intr%ipeva = state%atmosphere%intr%ipeva + X`
- `iptra` (line 469), `igrai` (line 476), `inrai` (line 478), `ievap` (line ~471)
- `cpeva` (line 499), `cptra` (line 500), `cevap` (line 501), `cgrai` (line 512), `cnrai` (line 513), `caintc` (line 510)
- `aintcdt` read at lines 467, 510 (as source for caintc accumulation)

Also: `ptra` is read at lines 396 and 442 for mass-balance output — migrate to `state%atmosphere%ptra`.

- [ ] **Confirm `state` in scope** in `integral` and `checkmassbal` — grep for `intent.*state` in waterbalance.f90.
- [ ] **Migrate all accumulator READ sides** to `state%atmosphere%intr%X` and `state%atmosphere%cumu%X` paths. The WRITE sides of these accumulations already write into state if dual-writes are active (they do — but the write side in waterbalance reads the legacy global and writes to the legacy global; after cutover, read from state and write to state).
- [ ] **Remove duplicate `igrai = 0.0d0` at line 388** and `inrai = 0.0d0` at line 389`. These are now owned by `state%atmosphere%intr%reset()` at soilhydraulics (A-2.1). Also remove from meteoday:414–415 (additive removal since meteoday's `ResetMetFlx` block was dual-patched in A-1.7).
- [ ] **Migrate `ptra` reads** at lines 396, 442.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-ATM Phase 2 A-2.2 — waterbalance: accumulator paths + igrai/inrai dedup

~14 accumulator read/write sites migrated to state%atmosphere%intr/cumu.
igrai/inrai double-reset eliminated: waterbalance:388-389 removed;
meteoday:414-415 (ResetMetFlx) removed. Single owner: intr%reset() in
soilhydraulics (A-2.1). ptra reads at lines 396+442 migrated. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md (ADR 0037)
EOF
)"
```

---

### Task A-2.3: `rootextraction.f90` + `cropgrowth.f90` reader cutover — `ptra`, `atmdem`

**Design decisions:** D1
**Files:**
- Modify: `src/crop/rootextraction.f90`
- Modify: `src/crop/cropgrowth.f90`

`rootextraction.f90` has ~28 read sites for `ptra` (25 sites) and `atmdem` (3 sites). `cropgrowth.f90` has ~9 `ptra` read sites. Both already take `state` per ADR 0036 — pure reader cutover, no signature work.

- [ ] **List all `ptra` read sites in rootextraction.f90** — `grep -n "\bptra\b" src/crop/rootextraction.f90 | grep -v "state%\|!"`. ~25 expected.
- [ ] **List all `atmdem` read sites** — same grep. ~3 expected.
- [ ] **Migrate all reads** to `state%atmosphere%ptra` and `state%atmosphere%atmdem`. Drop from `use variables, only:` clauses.
- [ ] **List all `ptra` read sites in cropgrowth.f90** (lines 669, 672, 715, 719, 722, 1712, 1715, 3004, 3007 per discovery). Migrate to `state%atmosphere%ptra`.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-ATM Phase 2 A-2.3 — rootextraction/cropgrowth: ptra/atmdem cutover

rootextraction.f90: ~28 ptra/atmdem read sites migrated to state%atmosphere.
cropgrowth.f90: ~9 ptra read sites migrated. Both files already state-plumbed
(ADR 0036) — pure reader cutover, no signature changes. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md (ADR 0037)
EOF
)"
```

---

### Task A-2.4: `boundtop.f90`, `temperature.f90`, `tillage.f90` reader cutover

**Design decisions:** D1
**Files:**
- Modify: `src/boundary/boundtop.f90`
- Modify: `src/heat/temperature.f90`
- Modify: `src/crop/tillage.f90`

Small files with few sites:
- `boundtop.f90`: `peva` at line 127, `empreva` at line 129, `melt` at line 137 and 180 (4 sites total). Already state-plumbed (ADR 0035).
- `temperature.f90`: `ssnow` at lines 165, 170 (2 sites). Already state-plumbed (ADR 0034).
- `tillage.f90`: `nraida` at lines 177, 178, 181, 370, 371 (5 sites). State plumbing: `tillage.f90` takes `state` per ADR 0020.

- [ ] **Migrate boundtop.f90 reads** — `peva` → `state%atmosphere%peva`; `empreva` → `state%atmosphere%empreva`; `melt` (2 sites) → `state%atmosphere%melt`. Drop from use clauses.
- [ ] **Migrate temperature.f90 reads** — `ssnow` (2 sites) → `state%atmosphere%ssnow`. Drop from use clauses.
- [ ] **Migrate tillage.f90 reads** — `nraida` (5 sites) → `state%atmosphere%nraida`. Drop from use clauses.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-ATM Phase 2 A-2.4 — boundtop/temperature/tillage reader cutover

boundtop.f90: peva+empreva+melt (4 sites). temperature.f90: ssnow (2 sites).
tillage.f90: nraida (5 sites). All read from state%atmosphere. All files
already state-plumbed — pure reader cutover. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md (ADR 0037)
EOF
)"
```

---

### Task A-2.5: `swapoutput.f90` + `swap_csv_output.f90` output cutover (~36 sites)

**Design decisions:** D1
**Files:**
- Modify: `src/io/swapoutput.f90`
- Modify: `src/io/swap_csv_output.f90`

This is the heaviest output task in this arc. `swapoutput.f90` has ~30 sites; `swap_csv_output.f90` has ~6.

`swapoutput.f90` read clusters:
- `.bal` cumulative composite block (lines 212–344): `cgrai`, `cnrai`, `caintc`, `cpeva`, `cptra`, `cevap`, `cgsnow`, `csubl`, `csnrai`, `cmelt`
- `.inc` intermediate block (lines 381–558): `igrai`, `inrai`, `ipeva`, `iptra`, `ievap`, `igsnow`, `isubl`, `isnrai`
- Snow output block: `ssnow`, `snowinco`, `gsnow`, `snrai`, `melt`, `subl`
- `sicact` at lines 1068, 1233
- Mini-sim writeback check: `swapoutput.f90:3745–3820` — confirm whether any atmosphere field (peva, ptra, ssnow, etc.) is snapshotted/restored. If so, retarget the restore write to state (per playbook lesson boundary #6).

`swap_csv_output.f90` clusters:
- EPOT column: `ipeva` (line ~8 or ~237)
- TPOT column: `iptra` (line ~12 or ~244)
- SSNOW column: `ssnow` (line ~256)
- dstor computation: may read `ssnow` or cumulatives (check)

- [ ] **Grep all atmosphere fields in swapoutput.f90:**
  ```bash
  grep -n "cgrai\|cnrai\|caintc\|cpeva\|cptra\|cevap\|cgsnow\|csubl\|csnrai\|cmelt\|igrai\|inrai\|ipeva\|iptra\|ievap\|igsnow\|isubl\|isnrai\|ssnow\|snowinco\|gsnow\|snrai\|melt\|subl\|sicact\|peva\|ptra\|empreva\|atmdem" src/io/swapoutput.f90 | grep -v "state%\|!" | head -60
  ```
- [ ] **Check mini-sim writeback block** at lines 3745–3820 for any atmosphere fields. Retarget any restore writes to `state%atmosphere%X` if found.
- [ ] **Migrate all output reads** to the appropriate `state%atmosphere%X`, `state%atmosphere%intr%X`, `state%atmosphere%cumu%X` paths. Use ASSOCIATE: `associate(at => state%atmosphere)` to keep lines compact.
- [ ] **Migrate swap_csv_output.f90 reads** (EPOT, TPOT, SSNOW, dstor).
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-ATM Phase 2 A-2.5 — swapoutput/swap_csv_output reader cutover

swapoutput.f90: ~30 sites migrated (.bal cumulative block, .inc intermediate
block, snow output, sicact). swap_csv_output.f90: ~6 sites (EPOT/TPOT/SSNOW/
dstor). Mini-sim writeback at 3745-3820 checked [and retargeted if needed].
ASSOCIATE at_ block keeps lines compact. check-full 5/5.

Spec: docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md (ADR 0037)
EOF
)"
```

---

### Task A-2.6: Drop dual-writes + retire ~33 globals from `variables.f90`; compile-driven fixups

**Design decisions:** D1, D13
**Files:**
- Modify: `src/atmosphere/et.f90`
- Modify: `src/atmosphere/interception.f90`
- Modify: `src/atmosphere/precipitation.f90`
- Modify: `src/atmosphere/snow.f90`
- Modify: `src/atmosphere/meteoday.f90`
- Modify: `src/atmosphere/meteodt.f90`
- Modify: `src/soil/soilhydraulics.f90`
- Modify: `src/core/variables.f90`
- Modify: `src/core/initialize.f90`

This task has the highest compile-driven discovery risk. After dropping dual-writes, expect 3–6 hidden readers to surface. Each compile-driven fixup is its own commit before finalizing the retirement commit.

Pre-flight — for each of the ~33 globals, run the two-grep test before deletion:
```bash
# For each field, e.g. peva, ptra, ssnow, etc.:
grep -rEn "use variables.*\bpeva\b" src/ --include="*.f90" | grep -v "variables.f90\|initialize.f90"
grep -rEn "\bpeva\b" src/ --include="*.f90" | grep -v "variables.f90\|initialize.f90\|state%atmosphere\|config%" | grep -v "!"
# Also check ASSOCIATE alias forms:
grep -rEn "=> state%atmosphere" src/ --include="*.f90"
```
Both greps must return zero hits before deletion. Investigate any remaining hit.

- [ ] **Drop dual-write legacy-global write side** in all 6 home-tree files (et.f90 reduceva path; interception.f90 DivIntercep; precipitation.f90 PartitionPrecipitation; snow.f90 all task branches; meteoday.f90 ProcessMeteoDay + ResetMetFlx; meteodt.f90 MeteoDT). Also remove the now-redundant legacy inline resets that survived as dual-patched pairs in meteoday and soilhydraulics (the `if (flzerointr) X=0 / if (flzerocumu) X=0` blocks adjacent to the cohort resets).
- [ ] **Compile** — `pixi run -e build` — fix compiler errors one by one. Each fix is its own commit labeled `refactor(state): SS-ATM compile-fix N — <symbol> in <file>`.
- [ ] **Comment out ~33 globals in `variables.f90`** with provenance markers:
  ```fortran
  ! real(8) peva  ! [SS-ATM] retired 2026-05-11 — moved to state%atmosphere%peva (ADR 0037)
  ```
  Full list: peva, pevaday, ptra, ptraday, empreva, atmdem, ldwet, spev, saev, grai, graidt, nraida, nraidt, aintcdt, fprecnosnow, sicact, ssnow, snowinco, slw, melt, subl, gsnow, snrai, igrai, inrai, ipeva, iptra, ievap, igsnow, isubl, isnrai, cgrai, cnrai, caintc, cpeva, cptra, cevap, cgsnow, csubl, csnrai, cmelt.
- [ ] **Drop zero-init lines from `initialize.f90`** for all retired globals (bulk zero-init block that becomes redundant after state defaults).
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full` (integration gate — must be byte-identical).

```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-ATM Phase 2 A-2.6 — retire ~33 atmosphere globals; drop dual-writes

Drop dual-writes from all 7 home-tree files. state%atmosphere is now
authoritative for all atmosphere-owned fields. Legacy globals commented
out in variables.f90 with [SS-ATM] provenance markers. Redundant zero-inits
dropped from initialize.f90. Compile-driven fixups resolved (N hidden readers
surfaced). pFUnit green. check-full 5/5 byte-identical.

Spec: docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md (ADR 0037)
EOF
)"
```

---

### Task A-2.7: ADR 0037 + playbook update + merge prep

**Design decisions:** all
**Files:**
- Create: `docs/adr/0037-state-migration-atmosphere.md`
- Modify: `docs/superpowers/specs/state-migration-playbook.md`
- Modify: `docs/adr/index.md` (add ADR 0037 entry)

- [ ] **Author ADR 0037** — status accepted; migration #7; two-cohort layout; meteoday/meteodt scope correction rationale; `atmosphere_init` (scalars only); snow.f90 signature promotion; igrai/inrai double-reset resolution (D6); ~33 globals retired; pond and nird/gird deferrals; compile-driven fixup count. Cross-reference discovery, design, plan, ADRs 0030–0036.
- [ ] **Update playbook** — append atmosphere-arc lessons. Key new patterns: (1) cohort pattern debuted at type-creation (not retrofitted); (2) scope correction mid-discovery (excluded-file → home-routine) and its effect on arc shape; (3) double-reset elimination via cohort ownership transfer; (4) multi-file cohort reset consolidation (3 files → 1 cohort reset() call site).
- [ ] **Update `docs/adr/index.md`** with ADR 0037 entry.
- [ ] **Final PASS** — `pixi run -e test test-pfunit && pixi run check-full`.

```bash
git commit -m "$(cat <<'EOF'
docs(adr): ADR 0037 — atmosphere state-type migration; playbook update

ADR 0037: atmosphere as migration #7; atmosphere_state_t with two cohorts
(intermediate 8 fields, cumulative 10 fields) + 22 flat scalars; first arc
to debut ADR 0033 cohort pattern at type-creation. meteoday/meteodt
included as full home routines (scope correction). igrai/inrai double-reset
eliminated. ~33 globals retired. pond/nird/gird deferred.

Playbook: new atmosphere lessons appended (cohort-at-creation, scope
correction, double-reset elimination, multi-file reset consolidation).

Spec: docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md (ADR 0037)
EOF
)"
```

---

## Self-Review Notes

- **A-1.7 (meteoday.f90) is the highest-risk Phase 1 task.** The file is 884 LoC with 9 peva write sites + 7 ptra write sites + 3 atmdem write sites spread across multiple subroutines (`ReadMeteoDay`, `ProcessMeteoDay`, `ResetMetFlx`, `ETSine`, etc.). Use the discovery's authoritative line numbers as anchors; do NOT rely solely on grep order. Missing one dual-write site produces a subtle regression where the state lags one day.

- **A-1.8 (meteodt.f90) — confirm pevaday/ptraday ownership.** Discovery indicates meteodt writes peva/ptra at lines 348/349 and 444/445 under the ETSine path. But pevaday/ptraday may be written in meteoday's `ETSine` helper OR in meteodt — the discovery is ambiguous. Before A-1.8, grep both files: `grep -n "\bpevaday\b\|\bptraday\b" src/atmosphere/meteoday.f90 src/atmosphere/meteodt.f90`. Whichever file writes them gets the dual-write in that task; update the task plan accordingly.

- **A-2.2 (igrai/inrai double-reset removal) must coordinate with A-1.7.** A-1.7 dual-patched meteoday:414–415 (the ResetMetFlx block); A-2.2 removes the duplicate from both waterbalance:388–389 AND meteoday:414–415. Confirm A-1.7 has landed before starting A-2.2's removal step. Do not remove until the cohort `intr%reset()` at soilhydraulics (A-2.1) is active.

- **A-2.5 (swapoutput.f90) — ASSOCIATE alias hazard.** Discovery playbook lesson (ADR 0033 Phase B finding): swapoutput.f90 may already use `ASSOCIATE(sw => state%surfacewater)` or `ASSOCIATE(sl => state%solute)` blocks. If it adds `at => state%atmosphere`, confirm there is no name collision with other subsystem aliases. Also check whether any existing ASSOCIATE block aliases a field that overlaps with atmosphere field names (e.g., an existing `at =>` for something else).

- **Compile-driven Phase 2.7 estimate (3–6 fixups).** Most likely surprise locations: `swapoutput.f90` use clauses (the `.bal` cumulative block imports may reference atmosphere fields via bare-name after old use-clauses, not just via the ASSOCIATE path); `macroporeoutput.f90` (discovery noted commented-out `ssnow` reference — may still have a live import); `config_to_variables.f90` for `ssnow`/`ldwet` initial seeds (Cat 4 init-seed pattern from playbook). Plan for up to 6 fixup commits before declaring A-2.6 done.

- **Per-event fields `ldwet`/`spev`/`saev` have a THIRD reset site** beyond the two in soilhydraulics (tillage reset at lines 873–875) and the tillage.f90 caller. The config-initial values are seeded from `soil_config.initial` → `config_to_variables.f90:565`. After migration, `config_to_variables.f90` must write `state%atmosphere%ldwet` (and `spev`/`saev`) instead of the legacy globals. This is a Cat 4 init-seed that will surface as a compile error during A-2.6 if not caught earlier. Check during A-0.1 pre-flight and add the dual-write to `config_to_variables.f90` during A-1.9 or A-2.4.

- **`snowinco` write at snow.f90:74 and read at snow.f90:99** — `snowinco` is written as a snapshot at task=1 init (`snowinco = ssnow` or `ssnow = snowinco` depending on `swinco`) and read as a baseline at the cumulative-reset-period start. After migration: `state%atmosphere%snowinco = state%atmosphere%ssnow` and the reverse. Ensure the ordering of dual-writes preserves the read-before-write or write-before-read semantics at these two lines. Check snow.f90:74 and 99 carefully during A-1.3.
