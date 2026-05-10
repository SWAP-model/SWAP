# Boundary-Conditions State-Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Carve 12 instantaneous-scalar boundary fields into `state%soilwater` (boundary subset); close 7–8 Phase 0 config gaps in `bottom_boundary_config_t`; preserve byte-identical check-full at every commit.

**Architecture:** Flat `soilwater_state_t` (initial form), aggregated under `swap_state_t`. `boundtop` and `BoundBottom` already take `state`; the arc is field-carve + co-writer dual-writes + reader migration, not signature surgery.

**Tech Stack:** Fortran 2008, ASSOCIATE, pixi+meson, pFUnit, check-full regression.

**Spec:** `docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md`
**Discovery:** `docs/superpowers/specs/2026-05-10-state-migration-boundary-discovery.md`
**Branch:** all commits go on `refactor/boundary-state` (new from `development`).

---

## Lessons-learned applied

- **Dual-write before cutover.** Every global write gets a paired `state%soilwater%X = …` mirror before readers are cut over. Drop globals only in Phase 2 after all readers are migrated.
- **Verify-state-not-reports.** Run `check-full` output comparison, not just pFUnit, before every commit. `pFUnit` alone misses global-default regressions.
- **ASSOCIATE prefix `bw_`.** Shadow-safe aliasing in compute bodies that also import `use variables`. Follows the `ht_` (heat), `dr_` (drainage), `sl_` (solute) convention.
- **Grep field shapes during pre-flight.** Before deleting a global, grep both `use variables.*\b<var>\b` AND raw-symbol forms. Both must be empty.
- **Defer hazards explicitly.** `pond`, `gwl`, `kmean` entries are deferred to soil-water-core per D5–D7. Do not migrate them here even if they appear near a boundary write site.

---

## File Structure

**Created:**
- `src/state/soilwater_state.f90` — `soilwater_state_t` definition (Phase 1 B-1.1)
- `tests/unit/state/test_soilwater_state.pf` — pFUnit tests (Phase 1 B-1.1)
- `tests/unit/config/test_boundary_phase0_fields.pf` — pFUnit config tests (Phase 0)
- `docs/adr/0035-state-migration-boundary.md` — finalized at end (Phase 2 B-2.7)

**Modified:**
- `src/config/bottom_boundary_config.f90` — add 8 Phase 0 fields (B-0.1–B-0.3)
- `src/io/toml/read_bottom_boundary_toml.f90` (or equivalent reader) — add Phase 0 reads (B-0.1–B-0.3)
- `src/io/toml/config_to_variables.f90` — adapter populates 8 legacy globals (B-0.4)
- `src/state/swap_state.f90` — add `type(soilwater_state_t) :: soilwater` (B-1.1)
- `src/core/swap.f90` — `soilwater_init(state)` call at line 183 (B-1.2)
- `src/boundary/boundtop.f90` — dual-write 7 top-boundary fields + ASSOCIATE bw_ (B-1.3)
- `src/boundary/boundbottom.f90` — dual-write 5 bottom-boundary fields + ASSOCIATE bw_ (B-1.3)
- `src/soil/soilhydraulics.f90` — dual-write qbot (12 sites) + qtop (1 site) (B-1.4); reader cutover (B-2.1)
- `src/heat/frozencond.f90` — dual-write qbot (3 sites in FrozenBounds) (B-1.4); cutover (B-2.5)
- `src/soil/waterbalance.f90` — reader cutover reva/runots/qbot (B-2.2)
- `src/solute/solute.f90` — reader cutover qtop/qbot/runots (B-2.3)
- `src/solute/agetracer.f90` — reader cutover qtop/qbot (B-2.3)
- `src/drainage/surfacewater.f90` — reader cutover runots (B-2.4)
- `src/drainage/drainage.f90` — reader cutover runots (B-2.4)
- `src/utils/surfacewaterutils.f90` — reader cutover (B-2.4)
- `src/io/swapoutput.f90` — output cutover + mini-sim writeback retarget (B-2.6)
- `src/io/swap_csv_output.f90` — output cutover (B-2.6)
- `src/core/variables.f90` — comment out 12 retired globals (B-2.7)
- `tests/unit/testSuites.inc`, `tests/unit/meson.build` — register new suites
- `meson.build` or equivalent source list — add `src/state/soilwater_state.f90`

---

# Phase 0 — Config gap closure

### Task B-0.1: Add `sinmax` / `sinamp` / `sinave` to `bottom_boundary_config_t`

**Design decisions:** D13
**Files:**
- Modify: `src/config/bottom_boundary_config.f90`
- Modify: `src/io/toml/read_bottom_boundary_toml.f90` (or equivalent reader — verify filename)
- Create: `tests/unit/config/test_boundary_phase0_fields.pf` (start the file here; extend in B-0.2 and B-0.3)
- Modify: `tests/unit/testSuites.inc`, `tests/unit/meson.build`

These three scalars support the `swbotb=2 .and. sw2=1` sine-wave bottom-flux branch (`boundbottom.f90:104`). Discovery Section 6 / `variables.f90` lines 947–949.

- [ ] **Inspect existing `bottom_boundary_config_t`** — read `src/config/bottom_boundary_config.f90` lines 1–80; note existing field patterns and validator style.
- [ ] **Write failing pFUnit tests** — in `test_boundary_phase0_fields.pf`, add `@test` stubs for `sinmax`, `sinamp`, `sinave`: default value (0.0), in-range, out-of-range. Register suite in `testSuites.inc`. Confirm FAIL with `pixi run -e test test-pfunit`.
- [ ] **Add fields + validators** to `bottom_boundary_config_t`: `real(real64) :: sinmax = 0.0_real64`, `sinamp`, `sinave`. Gate validator: active when `swbotb == 2 .and. sw2 == 1`.
- [ ] **Add TOML reads** — `get_optional_real_with_default` for each of the three fields.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
feat(config): SS-BND Phase 0 B-0.1 — sinmax/sinamp/sinave in bottom_boundary_config_t

Promotes three swbotb=2 sine-wave scalars (sinmax, sinamp, sinave)
from legacy variables.f90 into bottom_boundary_config_t with typed
defaults, validators, and TOML reader. Closes the swbotb=2 config
gap identified in the boundary-state discovery (Section 6).

check-full 5/5 byte-identical. Spec: docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md
EOF
)"
```

---

### Task B-0.2: Add `cofqha` / `cofqhb` / `cofqhc` / `swcofqhc` to `bottom_boundary_config_t`

**Design decisions:** D13
**Files:**
- Modify: `src/config/bottom_boundary_config.f90`
- Modify: TOML reader (same file as B-0.1)
- Modify: `tests/unit/config/test_boundary_phase0_fields.pf` (extend)

These four parameters support the `swbotb=4 .and. swqhbot=1` exponential q(h) branch (`boundbottom.f90:152,153`). Note: `drainage` config has a `cofqha_table` for surface runoff — a different field; do not confuse. Discovery Section 6 / `variables.f90` lines 761–763.

- [ ] **Write failing tests** — extend `test_boundary_phase0_fields.pf` with stubs for `cofqha`, `cofqhb`, `cofqhc` (real scalars) and `swcofqhc` (integer switch). Confirm FAIL.
- [ ] **Add fields** to `bottom_boundary_config_t`: `real(real64) :: cofqha = 0.0_real64`, `cofqhb`, `cofqhc`, `integer :: swcofqhc = 0`. Gate validators: active when `swbotb == 4 .and. swqhbot == 1`.
- [ ] **Add TOML reads** — optional real reads for cofqha/b/c; integer read for swcofqhc.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
feat(config): SS-BND Phase 0 B-0.2 — cofqha/b/c + swcofqhc in bottom_boundary_config_t

Promotes four swbotb=4 exponential q(h) parameters from legacy globals
into bottom_boundary_config_t. Distinct from drainage config's cofqha_table
(surface-runoff path). Spec: docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md
EOF
)"
```

---

### Task B-0.3: Add `hplate` to `bottom_boundary_config_t`

**Design decisions:** D13
**Files:**
- Modify: `src/config/bottom_boundary_config.f90`
- Modify: TOML reader
- Modify: `tests/unit/config/test_boundary_phase0_fields.pf` (extend)

`hplate` is used at `soilhydraulics.f90:271,557` for the `swbotb=8` lysimeter bottom boundary. Semantically a bottom-boundary parameter; belongs in `bottom_boundary_config_t`. Discovery Section 6 / `variables.f90` line 813.

- [ ] **Write failing test** — add stub for `hplate` (real scalar, default 0.0). Confirm FAIL.
- [ ] **Add field + validator** — `real(real64) :: hplate = 0.0_real64`; gate: `swbotb == 8`. Range check (negative ok — it is a pressure head).
- [ ] **Add TOML read** — `get_optional_real_with_default`.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
feat(config): SS-BND Phase 0 B-0.3 — hplate in bottom_boundary_config_t

Promotes swbotb=8 lysimeter parameter hplate from legacy variables.f90
into bottom_boundary_config_t. Used by soilhydraulics.f90:271,557.
Spec: docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md
EOF
)"
```

---

### Task B-0.4: Adapter — `config_to_variables.f90` populates all 8 new globals

**Design decisions:** D13
**Files:**
- Modify: `src/io/toml/config_to_variables.f90` — the `apply_bottom_boundary` (or equivalent) adapter section

Verify that the adapter function that handles `bottom_boundary_config_t` now populates all 8 newly-added legacy globals (`sinmax`, `sinamp`, `sinave`, `cofqha`, `cofqhb`, `cofqhc`, `swcofqhc`, `hplate`). Cross-check each global with its `variables.f90` line (discovery Section 6 table) and confirm it appears on the adapter's write side.

- [ ] **Grep adapter for existing bottom-boundary section** — `grep -n "bottom_boundary\|swbotb\|sinmax\|cofqha\|hplate" src/io/toml/config_to_variables.f90 | head -30`.
- [ ] **Add adapter assignments** for all 8 new fields — one assignment per field, mirroring the existing pattern (e.g., `sinmax = config%bottom_boundary%sinmax`).
- [ ] **Cross-check default values** against legacy `.swp` default or `variables.f90` initializer — confirm no silent-default deviation.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
feat(config): SS-BND Phase 0 B-0.4 — adapter wires 8 new bottom-boundary globals

config_to_variables.f90 now populates sinmax/sinamp/sinave/cofqha/cofqhb/
cofqhc/swcofqhc/hplate from bottom_boundary_config_t. Closes the full
Phase 0 adapter gap for the boundary-state arc. Spec: docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md
EOF
)"
```

---

### Task B-0.5: Regression-case audit — swbotb=2/4/8 coverage check

**Design decisions:** D13
**Files:** read-only audit (no code changes expected)

Grep TOML fixtures for which `swbotb` values are exercised. If swbotb=2 (sine), swbotb=4 (exp q(h)), or swbotb=8 (lysimeter) are absent from all five check-full cases, document as a known-issue. Do not add new fixtures in this arc.

- [ ] **Grep TOML fixtures** — `grep -rn "swbotb" tests/` (or wherever the TOML regression cases live). Note which values appear.
- [ ] **Grep legacy `.swp` fixtures if present** — confirm no parallel test set uses uncovered modes.
- [ ] **Document findings** — if any of swbotb=2/4/8 are uncovered, add a comment to `docs/adr/0035-state-migration-boundary.md` (stub may not yet exist) or to this plan's self-review section. No new fixture required this arc.
- [ ] **No commit needed** unless a trivially fixable gap is found. If commit is needed:
```bash
git commit -m "$(cat <<'EOF'
docs: SS-BND Phase 0 B-0.5 — regression-coverage audit result

Documents which swbotb modes are exercised in check-full regression.
Gaps (if any) noted as known-issues in ADR 0035 stub.
Spec: docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md
EOF
)"
```

---

# Phase 1 — State type, init, dual-write

### Task B-1.1: Create `soilwater_state_t` + add to `swap_state_t`

**Design decisions:** D2, D3, D8, D9
**Files:**
- Create: `src/state/soilwater_state.f90`
- Create: `tests/unit/state/test_soilwater_state.pf`
- Modify: `src/state/swap_state.f90`
- Modify: `tests/unit/testSuites.inc`, `tests/unit/meson.build`, meson source list

12-field flat type (design doc D2). No cohort sub-records. All scalars with zero/false defaults.

- [ ] **Write failing pFUnit tests** — default scalar values, logical defaults, independent-instance isolation. Confirm FAIL.
- [ ] **Create `soilwater_state_mod`** in `src/state/soilwater_state.f90` — 12 fields per D2 layout. File-header comment references ADR 0035 and the design spec.
- [ ] **Add `type(soilwater_state_t) :: soilwater` to `swap_state_t`** — after `heat` field; add `use soilwater_state_mod` import. Update meson source list.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
feat(state): SS-BND Phase 1 B-1.1 — soilwater_state_t flat type + swap_state aggregation

12-field flat soilwater_state_t (boundary subset): qtop, reva, hsurf,
runots, QMpLatSs, ftoph, FlRunoff, qbot, qbot_nonfrozen, hbot, gwlinp,
deepgw. All instantaneous scalars; no cohort sub-records (matches ADR 0034
heat precedent). Aggregated as state%soilwater in swap_state_t.

Spec: docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md
EOF
)"
```

---

### Task B-1.2: Add `soilwater_init(state)` and wire at `swap.f90:183`

**Design decisions:** D8
**Files:**
- Modify: `src/state/soilwater_state.f90` (add `soilwater_init` subroutine, or new `src/state/soilwater_init.f90`)
- Modify: `src/core/swap.f90`

All 12 fields are scalars — allocation is trivial (no `allocate` calls needed; default initializers handle it). `soilwater_init` may simply zero all fields and document the call-site contract. Placement at `swap.f90:183` (after `CalcGrid()`, before `DoTillage(1)`) ensures the state record is ready before any boundary code runs.

- [ ] **Inspect `swap.f90` lines 180–195** — confirm `CalcGrid()` at 182 and `DoTillage(1)` at 184; identify exact insertion point.
- [ ] **Implement `soilwater_init(state)`** — sets all 12 scalar fields to their zero/false defaults. Public subroutine in `soilwater_state_mod` (or a dedicated `soilwater_init_mod`).
- [ ] **Wire in `swap.f90`** — add `use soilwater_state_mod, only: soilwater_init` (or via `swap_state_mod`); insert `call soilwater_init(state)` between `CalcGrid()` and `DoTillage(1)`.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
feat(state): SS-BND Phase 1 B-1.2 — soilwater_init wired at swap.f90:183

soilwater_init(state) zeros all 12 soilwater_state_t scalar fields.
Called after CalcGrid() and before DoTillage(1) — earliest safe point
per D8 (forward-compatible with future per-node array additions).

Spec: docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md
EOF
)"
```

---

### Task B-1.3: Dual-write in boundary home tree (`boundtop`, `boundbottom`, `PONDRUNOFF`)

**Design decisions:** D1, D4, D10, D11
**Files:**
- Modify: `src/boundary/boundtop.f90`
- Modify: `src/boundary/boundbottom.f90`

For each of the 12 owned fields, find all write sites in the home tree and add a paired `state%soilwater%X = X` write immediately after each legacy global write. Use ASSOCIATE with `bw_` prefix at the top of each routine's compute body. Do NOT remove legacy global writes yet (that is Phase 2 Task B-2.7).

Fields to dual-write per file:
- `boundtop.f90`: `qtop`, `reva`, `hsurf`, `ftoph`, `FlRunoff`, `QMpLatSs`; `PONDRUNOFF` section: `runots` (discovery Section 2a write-site inventory).
- `boundbottom.f90`: `qbot` (11 sites), `qbot_nonfrozen` (1 site), `hbot` (1 site), `gwlinp` (1 site), `deepgw` (2 sites).

- [ ] **Grep write sites in home tree** — `grep -n "qtop\|reva\|hsurf\|ftoph\|FlRunoff\|QMpLatSs\|runots\|qbot\|hbot\|gwlinp\|deepgw" src/boundary/boundtop.f90 src/boundary/boundbottom.f90 | grep -v "!" | head -60`. Compare with discovery Section 2a table.
- [ ] **Add ASSOCIATE block** at top of each relevant subroutine body — `associate(bw_qbot => state%soilwater%qbot, bw_qtop => state%soilwater%qtop, …)`.
- [ ] **Add dual-writes** after each global write — `bw_X = X` (or `state%soilwater%X = X` directly).
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full` (must be byte-identical: dual-writes are additive).
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-BND Phase 1 B-1.3 — dual-write 12 fields in boundary home tree

boundtop.f90 and boundbottom.f90 (PONDRUNOFF included) write both
legacy global and state%soilwater%X for all 12 owned boundary fields.
ASSOCIATE bw_ prefix. Legacy globals unchanged — check-full byte-identical.

Spec: docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md
EOF
)"
```

---

### Task B-1.4: Dual-write at co-writer sites — `soilhydraulics.f90` and `frozencond.f90`

**Design decisions:** D1, D4
**Files:**
- Modify: `src/soil/soilhydraulics.f90`
- Modify: `src/heat/frozencond.f90`

`soilhydraulics.f90` co-writes `qbot` at 12 sites (lines 122, 252, 254, 257, 266, 274, 534, 537, 541, 552, 560, 721) and `qtop` at 1 site (line 637) inside `headcalc`. `frozencond.f90:FrozenBounds` co-writes `qbot` at 3 sites (lines 216, 243, 287). Both files already take `state` in scope. Do NOT touch `pond`, `gwlinp`, `hbot` writes in soilhydraulics here — only `qbot` and `qtop`.

- [ ] **Confirm state in scope** — `grep -n "intent.*state\|swap_state_t" src/soil/soilhydraulics.f90 src/heat/frozencond.f90 | head`. Both must show state as a dummy arg in relevant subroutines.
- [ ] **Add dual-writes in soilhydraulics** — at each of the 12 `qbot` write sites and 1 `qtop` write site, add `state%soilwater%qbot = qbot` (or via ASSOCIATE) immediately after the legacy write.
- [ ] **Add dual-writes in frozencond** — at `FrozenBounds` lines 216, 243, 287, add `state%soilwater%qbot = qbot`.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-BND Phase 1 B-1.4 — dual-write qbot/qtop at co-writer sites

soilhydraulics.f90 (12 qbot sites + 1 qtop site in headcalc) and
frozencond.f90:FrozenBounds (3 qbot sites) write both legacy global
and state%soilwater%X. Both files already take state. Legacy globals
unchanged — check-full byte-identical.

Spec: docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md
EOF
)"
```

---

# Phase 2 — Reader cutover and global retirement

### Task B-2.1: `soilhydraulics.f90` reader cutover

**Design decisions:** D1, D4
**Files:**
- Modify: `src/soil/soilhydraulics.f90`

Migrate all reads of `qtop`, `qbot`, `reva`, `hsurf`, `ftoph`, `FlRunoff`, `hbot`, `gwlinp`, `deepgw` in `soilhydraulics.f90` to `state%soilwater%X`. Discovery Section 3 reader index: ~25 read sites. State is already in scope throughout (`headcalc(state)`, `soilwater(task, state)`). Drop migrated symbols from `use variables, only:` clauses after all sites are updated.

- [ ] **List all read sites** — `grep -n "qtop\|reva\|hsurf\|ftoph\|FlRunoff\|hbot\|gwlinp\|deepgw" src/soil/soilhydraulics.f90 | grep -v "!" | grep -v "state%soilwater"` — compare with discovery table.
- [ ] **Migrate reads** — replace each bare-name read with `state%soilwater%X`. Use ASSOCIATE or direct access consistently.
- [ ] **Drop from `use variables, only:`** — verify no remaining bare-name reads.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-BND Phase 2 B-2.1 — soilhydraulics reads state%soilwater for 9 fields

All ~25 read sites for qtop/qbot/reva/hsurf/ftoph/FlRunoff/hbot/gwlinp/
deepgw in soilhydraulics.f90 (headcalc + SoilWater) now read state%soilwater.
Dropped from use variables. State authoritative for reads; dual-write still active.

Spec: docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md
EOF
)"
```

---

### Task B-2.2: `waterbalance.f90` reader cutover

**Design decisions:** D1, D6
**Files:**
- Modify: `src/soil/waterbalance.f90`

Migrate reads of `reva`, `runots`, `qbot` in `waterbalance.f90` (discovery Section 3: `integral`, `calcgwl`, `checkmassbal` — ~10 sites). Note: `gwl` is deferred (D6); do not touch `calcgwl`'s `gwl` write. Verify `state` is in scope in `integral` and `checkmassbal`; if not, add it and update callers.

- [ ] **Confirm state scope** — `grep -n "intent.*state\|subroutine integral\|subroutine checkmassbal" src/soil/waterbalance.f90 | head`. If missing, trace call chain from `swap.f90` to identify plumbing needed.
- [ ] **Migrate reva/runots/qbot reads** — replace bare names with `state%soilwater%X` in `integral` and `checkmassbal`. Do not touch `gwl` writes in `calcgwl`.
- [ ] **Update `use variables, only:`** — drop migrated symbols.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-BND Phase 2 B-2.2 — waterbalance reads state%soilwater (reva/runots/qbot)

integral() and checkmassbal() in waterbalance.f90 read reva, runots,
qbot from state%soilwater. calcgwl gwl writes untouched (gwl deferred
to soil-water-core per D6).

Spec: docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md
EOF
)"
```

---

### Task B-2.3: `solute.f90` + `agetracer.f90` reader cutover

**Design decisions:** D1
**Files:**
- Modify: `src/solute/solute.f90`
- Modify: `src/solute/agetracer.f90`

Migrate reads of `qtop`, `qbot`, `runots` (and `pond` if present — leave `pond` as legacy per D5). Discovery Section 3: ~7 read sites in `solute.f90`, 1–2 in `agetracer.f90`. Both already take `state` from the solute migration arc.

- [ ] **List read sites** — `grep -n "qtop\|qbot\|runots\|pond" src/solute/solute.f90 src/solute/agetracer.f90 | grep -v "state%\|!" | head -20`. Confirm `pond` reads remain as legacy.
- [ ] **Migrate qtop/qbot/runots reads** — `state%soilwater%X`. Leave `pond` reads on bare global.
- [ ] **Update `use variables, only:`** — drop migrated symbols only.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-BND Phase 2 B-2.3 — solute/agetracer reads state%soilwater

solute.f90 and agetracer.f90 read qtop/qbot/runots from state%soilwater.
pond reads left on legacy global (deferred per D5).

Spec: docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md
EOF
)"
```

---

### Task B-2.4: `surfacewater.f90`, `drainage.f90`, `surfacewaterutils.f90` reader cutover

**Design decisions:** D1, D5
**Files:**
- Modify: `src/drainage/surfacewater.f90`
- Modify: `src/drainage/drainage.f90`
- Modify: `src/utils/surfacewaterutils.f90`

`surfacewater.f90` has ~10 `runots` read sites (reservoir balance). `drainage.f90` has ~5 sites (runots indirect). `surfacewaterutils.f90` uses `pond` via `qhtab` — leave `pond` as legacy per D5. Verify `state` is in scope in all three.

- [ ] **List read sites** — `grep -n "runots\|qtop\|qbot\|pond" src/drainage/surfacewater.f90 src/drainage/drainage.f90 src/utils/surfacewaterutils.f90 | grep -v "state%\|!" | head -30`. Note which fields to migrate vs leave.
- [ ] **Migrate runots reads** in surfacewater.f90 and drainage.f90 — `state%soilwater%runots`.
- [ ] **Leave `pond` reads on legacy global** in all three files (D5).
- [ ] **Update `use variables, only:`** — drop `runots` from clauses where fully migrated.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-BND Phase 2 B-2.4 — surfacewater/drainage reads state%soilwater (runots)

surfacewater.f90 and drainage.f90 read runots from state%soilwater.
pond reads in surfacewater.f90 and surfacewaterutils.f90 left on legacy
global (deferred to soil-water-core per D5).

Spec: docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md
EOF
)"
```

---

### Task B-2.5: `frozencond.f90:FrozenBounds` cutover — `qbot` + `qbot_nonfrozen`

**Design decisions:** D9
**Files:**
- Modify: `src/heat/frozencond.f90`

`FrozenBounds` reads `qbot_nonfrozen` at line 216 (sole read) and reads/writes `qbot` at lines 216, 243, 287. State is already in scope. Cut both reads and the co-writes from legacy globals to `state%soilwater%qbot` and `state%soilwater%qbot_nonfrozen`.

- [ ] **Confirm write sites** — `grep -n "qbot\|qbot_nonfrozen" src/heat/frozencond.f90 | grep -v "!"`. Cross-check with discovery Section 2a and co-writers table.
- [ ] **Cut over reads** — `state%soilwater%qbot`, `state%soilwater%qbot_nonfrozen`.
- [ ] **Cut over co-writes** — the FrozenBounds frost-override writes to `state%soilwater%qbot` (dual-write from B-1.4 becomes state-only write).
- [ ] **Drop from `use variables, only:`**.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-BND Phase 2 B-2.5 — frozencond:FrozenBounds uses state%soilwater

FrozenBounds reads and writes qbot/qbot_nonfrozen via state%soilwater.
Dual-write dropped; state-only writes. qbot_nonfrozen retained in state
for symmetry (D9).

Spec: docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md
EOF
)"
```

---

### Task B-2.6: `swapoutput.f90` + `swap_csv_output.f90` cutover + mini-sim writeback retarget

**Design decisions:** D11
**Files:**
- Modify: `src/io/swapoutput.f90`
- Modify: `src/io/swap_csv_output.f90`

Output cutover: migrate all read sites for the 12 owned fields in output routines to `state%soilwater%X`. Discovery Section 3 Cat 1 reads: ~15 sites in `swapoutput.f90`, ~6 in `swap_csv_output.f90`.

Mini-sim writeback retarget (discovery Hazard #1): `swapoutput.f90:3745–3820` snapshots `qbot`/`gwl`/`pond` before a perturbation mini-sim, then writes back. After this task, the `qbot` writeback targets `state%soilwater%qbot`; `gwl` and `pond` writebacks stay on legacy globals (deferred per D5–D6).

- [ ] **List output read sites** — `grep -n "qtop\|qbot\|reva\|hsurf\|ftoph\|FlRunoff\|runots\|QMpLatSs\|hbot\|gwlinp\|deepgw\|qbot_nonfrozen" src/io/swapoutput.f90 src/io/swap_csv_output.f90 | grep -v "state%\|!" | head -40`.
- [ ] **Migrate output reads** — `state%soilwater%X` for all 12 fields.
- [ ] **Retarget mini-sim writeback** — in `swapoutput.f90:3818` (qbot writeback): change from bare `qbot = qbot_snapshot` to `state%soilwater%qbot = qbot_snapshot`. Leave `gwl`/`pond` writebacks on legacy.
- [ ] **Drop from `use variables, only:`** for migrated symbols.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-BND Phase 2 B-2.6 — output reads state%soilwater; mini-sim qbot retargeted

swapoutput.f90 and swap_csv_output.f90 read all 12 boundary fields from
state%soilwater. Mini-sim writeback at swapoutput.f90:3818 retargets qbot
to state%soilwater%qbot; gwl/pond writebacks remain on legacy globals
(deferred per D5-D6, same pattern as drainage mini-sim resolution).

Spec: docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md
EOF
)"
```

---

### Task B-2.7: Drop dual-writes + retire 12 globals from `variables.f90`

**Design decisions:** D1, D12
**Files:**
- Modify: `src/boundary/boundtop.f90`
- Modify: `src/boundary/boundbottom.f90`
- Modify: `src/soil/soilhydraulics.f90`
- Modify: `src/heat/frozencond.f90`
- Modify: `src/core/variables.f90`
- Create stub: `docs/adr/0035-state-migration-boundary.md` (initial draft)

Pre-flight: for each of the 12 fields, run both grep forms (see below) — both must be empty before any deletion.

- [ ] **Pre-flight grep for each field** — for `X` in `qtop qbot qbot_nonfrozen hbot gwlinp deepgw reva hsurf ftoph runots FlRunoff QMpLatSs`:
  ```bash
  grep -rEn "use variables.*\b${X}\b" src/ --include="*.f90" | grep -v "variables.f90\|initialize.f90"
  grep -rEn "\b${X}\b" src/ --include="*.f90" | grep -v "variables.f90\|initialize.f90\|state%soilwater\|config%" | grep -v "!"
  ```
  Both must return zero hits. Investigate any remaining hits before proceeding.
- [ ] **Drop dual-write legacy-global writes** in `boundtop.f90`, `boundbottom.f90`, `soilhydraulics.f90`, `frozencond.f90` — remove the bare-name global write side; keep `state%soilwater%X =` writes.
- [ ] **Comment out 12 globals in `variables.f90`** — with provenance markers:
  ```fortran
  ! real(8) qtop    ! [SS-BND] retired 2026-05-10 — moved to state%soilwater%qtop (ADR 0035)
  ```
- [ ] **Drop corresponding zero-init lines in `initialize.f90`** if present.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full` (integration gate — must be byte-identical).
- [ ] **Author ADR 0035 stub** — create `docs/adr/0035-state-migration-boundary.md` with status accepted, migration #5, key decisions (D1–D13), consequences, known residuals (pond/gwl/kmean deferred; swbotb=-2 mutation; any swbotb=2/4/8 coverage gaps from B-0.5). Reference discovery, design, plan.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-BND Phase 2 B-2.7 — retire 12 boundary globals; ADR 0035

Drop dual-writes; state%soilwater is now authoritative for all 12
boundary fields. Legacy globals commented out in variables.f90 with
[SS-BND] provenance markers. swbotb=-2 runtime mutation retained as
legacy config mutation (D12). pond/gwl/kmean deferred to soil-water-core.

pFUnit green. check-full 5/5 byte-identical.
ADR: docs/adr/0035-state-migration-boundary.md

Spec: docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md
EOF
)"
```

---

## Self-Review Notes

- **Spec coverage:** D1 (scope) → all tasks. D2 (flat shape) → B-1.1. D3 (aggregation) → B-1.1. D4 (argument-threading, bw_ ASSOCIATE) → B-1.3/B-1.4. D5 (pond deferred) → explicitly called out in B-2.3/B-2.4/B-2.6/B-2.7. D6 (gwl deferred) → B-2.2/B-2.7. D7 (kmean deferred) → pre-flight check in B-2.7. D8 (soilwater_init placement) → B-1.2. D9 (qbot_nonfrozen retained) → B-2.5. D10 (reva boundary-owned) → included in home-tree dual-write B-1.3. D11 (mini-sim preserved/retargeted) → B-2.6. D12 (swbotb=-2 kept) → documented in B-2.7 commit and ADR 0035 stub. D13 (Phase 0 config gaps) → B-0.1–B-0.5.
- **Task that may need splitting:** B-2.1 (`soilhydraulics.f90` ~25 read sites) is the highest-density single-file task. If soilhydraulics has grown since discovery or if there are co-write entanglements, the implementer should split into `headcalc` reads vs `SoilWater(1)` init reads as two sub-commits.
- **B-0.5 (regression audit):** if swbotb=2/4/8 are all uncovered, the audit produces only a documentation artifact. No test fixture is required this arc — confirm this is acceptable before proceeding to Phase 1.
- **`waterbalance.f90` state scope (B-2.2):** `calcgwl()` is a no-arg subroutine; `integral()` may or may not already take `state`. Verify before B-2.2; if plumbing is needed, it adds ~30 min of work not scoped in the heat-plan template.
