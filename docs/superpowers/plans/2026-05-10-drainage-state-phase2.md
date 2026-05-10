# Drainage State Migration — Phase 2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Complete the drainage state migration. After Phase 2: zero drainage-owned globals remain in `variables.f90`, every reader of those variables across the codebase reads from `state%drainage`, the `divdra` call sites pass typed state directly, and the pre-existing `geofac` adapter ordering bug is fixed.

**Architecture:** Migrate the remaining qdrain/qdra cross-subsystem readers (frozencond, surfacewater, waterbalance) to `state%drainage`. Update `divdra` callers to pass `state%drainage%qdra` / `state%drainage%qdrain` (the assumed-shape signature from Phase 1 Task 4 enables this). Drop drainage compute's dual-write of `qdra` and `qdrain`. Remove 6 owned globals from `variables.f90`. Fix `geofac` ordering bug at `config_to_variables.f90:431` (one-line gate).

**Tech Stack:** Fortran 2008, gfortran, meson + ninja + pixi, pFUnit, check-full byte-identical regression.

**Spec:** `docs/superpowers/specs/2026-05-10-state-migration-drainage-design.md` (sections D6, D7, D9 + Phase 2 phasing)

**Branch:** all commits go on `refactor/surfacewater-state` (the umbrella migration branch).

---

## Lessons-learned applied from Phase 1

- **Cross-subsystem reader inventory enumerated upfront.** Phase 1's discovery missed compute readers of `qdrd` (WLEVBAL/WBALLEV in surfacewater.f90 read it globally), which Task 7's check-full failure surfaced. Phase 2 plan starts with a complete pre-flight inventory.
- **Dual-write drop is the integration gate.** Each task migrates one consumer; dual-write drop (Task 4) is the proof that all consumers are migrated.
- **`divdra` already accepts assumed-shape arrays** (Phase 1 Task 4) — Phase 2 just switches call sites from passing legacy globals to passing typed state.

---

## Pre-flight inventory (all remaining readers/writers of qdra and qdrain)

`qdra` is already fully migrated read-side (Phase 1 Task 3 renamed all consumers to `state%drainage%qdra`). The only remaining work for `qdra` is dropping the legacy global write in compute and updating `divdra` callers.

`qdrain` is still largely on globals:

| File | qdrain reference | Type | Phase 2 task |
|---|---|---|---|
| `src/heat/frozencond.f90:251, 258, 265, 269, 276, 289, 294, 303` | reads + writes (FrozenBounds) | compute | Task 1 |
| `src/drainage/surfacewater.f90:143` | passed to divdra | call site | Task 3 |
| `src/drainage/surfacewater.f90:177, 219, 227` | reads (sumqdr, qdra fill, qdrtot accumulation) | compute | Task 2 |
| `src/soil/waterbalance.f90:513, 514, 516, 517, 519` | reads (cqdrain accumulation in integral) | compute | Task 2 |
| `src/drainage/drainage.f90:174` | `use variables, only: …, qdrain, …` in bocodrb | static read of drain geometry — KEEP (this is a config-time read, not the runtime per-level array) | (none — legitimate) |

For `qdra`:

| File | qdra reference | Type | Phase 2 task |
|---|---|---|---|
| `src/drainage/surfacewater.f90:143, 219` | passed to divdra (143); written from qdrain (219) | call site / dual-write | Tasks 3, 4 |
| `src/heat/frozencond.f90:276, 294` | passed to divdra (276); read in qdrain accumulation (294) | call site / read | Tasks 1, 3 |
| `src/drainage/drainage.f90:511` | passed to divdra | call site | Task 3 |

After Tasks 1-3: every `qdrain` and `qdra` consumer reads from `state%drainage`. After Task 4: drainage compute no longer writes legacy globals. After Task 5: the 6 globals are deleted from `variables.f90`.

---

## File Structure

**Modified:**
- `src/heat/frozencond.f90` — `FrozenBounds` reads + writes `qdrain` and `qdra` from/to `state%drainage` (Task 1)
- `src/drainage/surfacewater.f90` — `WLEVBAL` `qdrain` reads → `state%drainage%qdrain`; pass typed state to divdra (Tasks 2, 3)
- `src/soil/waterbalance.f90` — `integral`'s `qdrain` reads → `state%drainage%qdrain` (Task 2)
- `src/drainage/drainage.f90` — `bocodre` writes `state%drainage%qdrain` only (drops dual-write); pass typed state to divdra (Tasks 3, 4)
- `src/drainage/divdra.f90` — no change (assumed-shape from Phase 1 Task 4 already accepts both globals and typed state)
- `src/core/variables.f90` — delete 6 owned globals (Task 5)
- `src/core/initialize.f90` — drop initialization lines (Task 5)
- `src/io/toml/config_to_variables.f90:431` — gate `geofac` write on `ipos /= 5` (Task 5)
- `docs/adr/0031-state-migration-drainage.md` — Phase 2 outcomes finalized (Task 6)
- `docs/adr/index.md` — add ADR 0031 entry (Task 6)
- `docs/superpowers/specs/2026-05-09-state-migration-drainage-discovery.md` — Section 3.5 Phase 1 lessons-learned addendum (Task 6)

---

### Task 1: Migrate `frozencond.f90:FrozenBounds` qdrain and qdra reads/writes to `state%drainage`

`FrozenBounds` is a post-processor that runs after drainage compute and before SoilWater(2). It modifies `qdrain` (zeroes drains below frozen zone, redistributes flux to deepest unfrozen level) and then calls `DIVDRA` to recompute `qdra`. It also sums `qdrain` into `state%surfacewater%qdrtot`.

After this task: `FrozenBounds` reads/writes `state%drainage%qdrain` and `state%drainage%qdra`. Legacy globals are no longer touched by frozencond.

**Files:**
- Modify: `src/heat/frozencond.f90`

- [ ] **Step 1: Inspect FrozenBounds**

```bash
grep -n "subroutine FrozenBounds\|subroutine FrozenCond\|qdrain\|qdra\b" src/heat/frozencond.f90 | head -20
sed -n '240,310p' src/heat/frozencond.f90
```

Identify the qdrain write/read sites and the DIVDRA call.

- [ ] **Step 2: Confirm state argument is in scope**

Phase 2 of surfacewater migration plumbed state through frozencond. Verify:

```bash
grep -n "state\b" src/heat/frozencond.f90 | head
```

`FrozenBounds` should already have `type(swap_state_t), intent(inout) :: state` as an argument. If not, plumb it from the caller (likely `swap.f90`).

- [ ] **Step 3: Migrate qdrain references**

Inside `FrozenBounds`, replace:

```fortran
qdrain(level) = 0.0d0          ! line 251
qdratot = qdratot + qdrain(level)   ! line 258
qdrain(leveldeepest) = qbot         ! line 265
qdrain(level) = qdrain(level) * (1.0d0 + qbot/qdratot)   ! line 269
qdrain(level) = 0.0d0          ! line 289
qdrain(level) = qdrain(level) + qdra(level,node)         ! line 294
state%surfacewater%qdrtot = state%surfacewater%qdrtot + qdrain(level)   ! line 303
```

with `state%drainage%qdrain(...)` substitutions. Use ASSOCIATE if it cleans things up:

```fortran
associate(qdrain => state%drainage%qdrain, qdra => state%drainage%qdra)
   ! body unchanged
end associate
```

- [ ] **Step 4: Update the DIVDRA call (line ~276)**

Currently passes legacy globals `qdrain, qdra`. Change to:

```fortran
call divdra(..., state%drainage%qdrain, state%drainage%qdra, ...)
```

(The assumed-shape signature from Phase 1 Task 4 accepts these directly.)

- [ ] **Step 5: Drop qdrain and qdra from `use variables, only:` clause**

```bash
grep -n "use variables" src/heat/frozencond.f90 | head
```

Remove `qdrain` and `qdra` from the only-list if present.

- [ ] **Step 6: Build and verify**

Run: `pixi run check-full`
Expected: 5/5 byte-identical. Drainage compute still dual-writes globals (Task 4 hasn't dropped that yet), so frozencond's transition from reading globals to reading state is invisible to other consumers.

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: zero failures.

- [ ] **Step 7: Commit**

```bash
git add src/heat/frozencond.f90
git commit -m "$(cat <<'EOF'
refactor(state): SS-DRST Phase 2 Task 1 — frozencond reads qdrain/qdra from state

FrozenBounds now reads and writes state%drainage%qdrain and
state%drainage%qdra instead of legacy globals. The DIVDRA call
inside FrozenBounds passes typed state directly.

Spec: docs/superpowers/specs/2026-05-10-state-migration-drainage-design.md
EOF
)"
```

---

### Task 2: Migrate surfacewater.f90 + waterbalance.f90 qdrain reads to `state%drainage`

Two files, one task — both read `qdrain` for compute purposes and need the same rename.

**Files:**
- Modify: `src/drainage/surfacewater.f90` — qdrain reads in WLEVBAL (sumqdr, qdra fill, qdrtot accumulation)
- Modify: `src/soil/waterbalance.f90` — qdrain reads in integral (cqdrain accumulation)

- [ ] **Step 1: Inspect surfacewater.f90 qdrain reads**

```bash
grep -n "\bqdrain\b" src/drainage/surfacewater.f90 | head -10
```

The non-call-site reads are at lines 177, 219, 227 (per pre-flight inventory). Line 143 is a divdra call site — Task 3 handles that.

- [ ] **Step 2: Migrate surfacewater.f90 qdrain reads**

Each read of `qdrain(level)` becomes `state%drainage%qdrain(level)`. Use ASSOCIATE if helpful:

```fortran
associate(qdrain => state%drainage%qdrain)
   ! body unchanged
end associate
```

Drop `qdrain` from any `use variables, only:` clause in surfacewater.f90 if no longer referenced (the divdra call at line 143 passes legacy `qdrain` for now — Task 3 changes that — so keep the import until Task 3 lands; OR migrate line 143 in this task too if cleaner).

**Recommendation:** migrate line 143 here (pass `state%drainage%qdrain` and `state%drainage%qdra` to divdra) since you're already in the file. Then drop both from the `use variables, only:` list.

- [ ] **Step 3: Migrate waterbalance.f90 qdrain reads in integral**

```bash
grep -n "\bqdrain\b" src/soil/waterbalance.f90 | head
```

Lines 513, 514, 516, 517, 519 — all reads of `qdrain(level)` for cqdrain* accumulation. Replace each with `state%drainage%qdrain(level)`.

The accumulation pattern:
```fortran
if (qdrain(level).lt.0.0d0) then
   state%surfacewater%cqdrainin(level) = state%surfacewater%cqdrainin(level) - qdrain(level)*dt
else if (qdrain(level).gt.0.0d0) then
   state%surfacewater%cqdrainout(level) = state%surfacewater%cqdrainout(level) + qdrain(level)*dt
state%surfacewater%cqdrain(level) = state%surfacewater%cqdrain(level) + qdrain(level)*dt
```

becomes:
```fortran
if (state%drainage%qdrain(level).lt.0.0d0) then
   state%surfacewater%cqdrainin(level) = state%surfacewater%cqdrainin(level) - state%drainage%qdrain(level)*dt
else if (state%drainage%qdrain(level).gt.0.0d0) then
   state%surfacewater%cqdrainout(level) = state%surfacewater%cqdrainout(level) + state%drainage%qdrain(level)*dt
state%surfacewater%cqdrain(level) = state%surfacewater%cqdrain(level) + state%drainage%qdrain(level)*dt
```

Or with ASSOCIATE for readability:
```fortran
associate(qdrain => state%drainage%qdrain)
   ! body unchanged
end associate
```

Drop `qdrain` from `use variables, only:` in waterbalance.f90 if no longer referenced.

- [ ] **Step 4: Build and verify**

```bash
pixi run check-full
pixi run -e test test-pfunit 2>&1 | tail -10
```

Expected: 5/5 byte-identical, zero pFUnit failures.

- [ ] **Step 5: Commit**

```bash
git add src/drainage/surfacewater.f90 src/soil/waterbalance.f90
git commit -m "$(cat <<'EOF'
refactor(state): SS-DRST Phase 2 Task 2 — surfacewater + waterbalance read qdrain from state

WLEVBAL (surfacewater.f90) and integral (waterbalance.f90) now read
qdrain from state%drainage instead of the legacy global. The divdra
call at surfacewater.f90:143 passes typed state directly (the
assumed-shape signature from Phase 1 Task 4 accepts it).

Spec: docs/superpowers/specs/2026-05-10-state-migration-drainage-design.md
EOF
)"
```

---

### Task 3: Migrate remaining `divdra` call sites to pass typed state

After Tasks 1, 2, the only remaining legacy-global call sites for `divdra` are inside `drainage.f90` itself. This task switches drainage's own divdra call to pass typed state, completing the divdra migration.

**Files:**
- Modify: `src/drainage/drainage.f90` (the `divdra` call inside `drainage()`)

- [ ] **Step 1: Find divdra call sites still passing legacy globals**

```bash
grep -rEn "\bcall divdra\b" src/ --include="*.f90"
```

Verify Tasks 1, 2 already updated frozencond and surfacewater call sites. The remaining one(s) should be in `drainage.f90`.

- [ ] **Step 2: Migrate drainage.f90's divdra call**

Find the call (likely around line 511 per inventory):
```fortran
call divdra(..., qdrain, qdra, ...)
```

Change to:
```fortran
call divdra(..., state%drainage%qdrain, state%drainage%qdra, ...)
```

`state` is already an argument of `Drainage` (added during surfacewater Phase 2 Task 3).

- [ ] **Step 3: Drop qdrain/qdra from drainage.f90's `use variables, only:` clauses**

For each subroutine in drainage.f90 that references `qdrain` or `qdra` only via the migrated call site, drop them from the only-list. Verify no other references remain in that subroutine.

- [ ] **Step 4: Build and verify**

```bash
pixi run check-full
pixi run -e test test-pfunit 2>&1 | tail -10
```

Expected: 5/5 byte-identical.

- [ ] **Step 5: Commit**

```bash
git add src/drainage/drainage.f90
git commit -m "$(cat <<'EOF'
refactor(state): SS-DRST Phase 2 Task 3 — divdra call sites pass typed state

The divdra call inside drainage() now passes state%drainage%qdrain
and state%drainage%qdra directly (assumed-shape signature from
Phase 1 Task 4). After this commit, no divdra call site passes
legacy globals — Task 4 can drop the dual-write.

Spec: docs/superpowers/specs/2026-05-10-state-migration-drainage-design.md
EOF
)"
```

---

### Task 4: Drop drainage compute's qdrain and qdra dual-write — INTEGRATION GATE

After Tasks 1-3 every consumer reads from `state%drainage`. Drainage compute can stop writing the legacy globals. This is the integration gate for Phase 2.

**Files:**
- Modify: `src/drainage/drainage.f90` — remove legacy global writes for qdrain and qdra
- Modify: `src/drainage/surfacewater.f90` — drop the `qdra(level,numnod) = qdrain(level)` legacy write at line 219 (Task 2 may have already migrated this)
- Modify: tests if applicable

- [ ] **Step 1: Find dual-write sites in drainage.f90**

```bash
grep -nE "^\s+(qdrain|qdra)\s*[=\(]" src/drainage/drainage.f90 | head -30
grep -nE "^\s+state%drainage%(qdrain|qdra)\s*[=\(]" src/drainage/drainage.f90 | head -30
```

Identify pairs: `<global> = <RHS>` paired with `state%drainage%<global> = <RHS>`. Delete the legacy global write; keep the state write.

- [ ] **Step 2: Find any remaining writes elsewhere**

```bash
grep -rEn "^\s+qdrain\s*[=\(]|^\s+qdra\s*[=\(]" src/ --include="*.f90" | grep -v "state%"
```

Should return empty or only state-side writes. If a non-state global write remains in a file Tasks 1-3 didn't touch, investigate (probably a missed test fixture).

- [ ] **Step 3: Trim use variables**

In each file modified by Tasks 1-3 + this task, drop `qdrain` and `qdra` from `use variables, only:` clauses where no longer referenced. The only legitimate remaining import is `bocodrb` reading `qdrain` static drain geometry at drainage.f90:174 — verify by inspection that it's actually `qdrain` the array (not the runtime per-level flux). If it's a static read (e.g., for Hooghoudt geometry), keep it — Task 5 handles the variables.f90 declaration.

- [ ] **Step 4: Build and verify — THE INTEGRATION GATE**

```bash
pixi run check-full
```

Expected: `Results: 5 passed, 0 failed`. **This is the proof Phase 2's reader migration is complete** — no consumer needs the legacy globals.

```bash
pixi run -e test test-pfunit 2>&1 | tail -10
```

Expected: zero failures.

If a case fails: a missed reader. Find with:
```bash
grep -rEn "use variables.*\b(qdrain|qdra)\b" src/ --include="*.f90" | grep -v "src/core/variables.f90"
```

The hit is the unmigrated reader. Migrate or restore that file's dual-write.

- [ ] **Step 5: Commit**

```bash
git add src/drainage/drainage.f90 src/drainage/surfacewater.f90
# Plus any other files modified
git commit -m "$(cat <<'EOF'
refactor(state): SS-DRST Phase 2 Task 4 — drop qdra/qdrain dual-write — state authoritative

After Tasks 1-3 migrated all readers (frozencond, surfacewater
WLEVBAL, waterbalance integral, all divdra call sites), drainage
compute no longer writes qdra and qdrain to legacy globals. State
is authoritative for all 6 drainage-owned fields.

Integration gate: check-full byte-identical confirms no consumer
needs the globals. Task 5 deletes the declarations from
variables.f90.

Spec: docs/superpowers/specs/2026-05-10-state-migration-drainage-design.md
EOF
)"
```

---

### Task 5: Delete 6 owned globals from variables.f90 + fix geofac adapter

Final cleanup. Remove the 6 declarations + matching initialize.f90 lines + fix the `geofac` adapter ordering bug at `config_to_variables.f90:431` (one-line gate).

**Files:**
- Modify: `src/core/variables.f90` — delete (or comment out) `qdrd`, `qdrain(Madr)`, `qdra(Madr,macp)`, `drainl(Madr)`, `wetper(Madr)`, `ztopdislay(Madr)` declarations
- Modify: `src/core/initialize.f90` — drop the corresponding zero-init lines
- Modify: `src/io/toml/config_to_variables.f90:431` — gate `geofac` write on `config%drain%ipos /= 5`

- [ ] **Step 1: Confirm zero remaining global readers/writers**

```bash
for v in qdrd drainl wetper ztopdislay qdrain qdra; do
  count=$(grep -rEn "\b$v\b" src/ --include="*.f90" 2>/dev/null \
    | grep -v "src/core/variables.f90" \
    | grep -v "src/core/initialize.f90" \
    | grep -v "src/state/" \
    | grep -v "state%" \
    | grep -v "config%drain" \
    | wc -l)
  if [ "$count" -gt 0 ]; then
    echo "$count  $v — STILL HAS REFS"
  fi
done
```

Expected: zero counts. The exception is `bocodrb` in drainage.f90:174 which reads `qdrain` for static drain geometry (Hooghoudt formula). Inspect that reference: if it's reading `qdrain(level)` at a *static* moment (geometry-only, not the runtime flux), it's safe — but the `qdrain` global is shared between static config-derived values and runtime flux. Investigate whether the static read can come from `state%drainage%qdrain` (after `drainage_init` seeds it from globals) or from `config%drain%...` directly.

If `bocodrb`'s read is truly using runtime flux values (not static geometry), the implementer should fail Task 5 here and re-open a smaller task to migrate `bocodrb`. If it's static, keep `qdrain` import in that one place — but the global declaration in variables.f90 must remain. **Investigate carefully.**

- [ ] **Step 2: Delete or comment out the 6 declarations**

Match the surfacewater pattern from Phase 1: comment out (not delete) with provenance markers:

```fortran
! real(8)   qdrd               ! Moved to drainage_state_t%qdrd (ADR 0031)
! real(8)   drainl(Madr)       ! Moved to drainage_state_t%drainl
! real(8)   qdrain(Madr)       ! Moved to drainage_state_t%qdrain
! real(8)   qdra(Madr,macp)    ! Moved to drainage_state_t%qdra
! real(8)   wetper(Madr)       ! Moved to drainage_state_t%wetper
! real(8)   ztopdislay(Madr)   ! Moved to drainage_state_t%ztopdislay
```

Find each declaration (per Phase 1 Task 8 verification grep output: lines 784, 880, 881, 967, 976, 1294 area). Comment them out with the provenance comment.

If `bocodrb`'s `qdrain` static-geometry read genuinely needs the global, KEEP `qdrain` declared (don't comment it out) and document the exception in the commit message.

- [ ] **Step 3: Drop initialization lines from initialize.f90**

```bash
grep -nE "^\s+(qdrd|drainl|wetper|ztopdislay|qdrain|qdra)\s*=" src/core/initialize.f90
```

Each `<varname> = 0.0` (or similar) initializer for migrated variables: delete or comment out. Keep `qdrain` if Step 1's investigation determined it must stay.

- [ ] **Step 4: Fix geofac adapter ordering**

In `src/io/toml/config_to_variables.f90` around line 431, find:

```fortran
geofac       = config%drain%surface_runoff%geofac
```

Wrap with the ipos guard:

```fortran
if (config%drain%ipos /= 5) then
   geofac       = config%drain%surface_runoff%geofac
end if
```

Add a comment explaining the conflict (line 306 already wrote `geofac` for `ipos==5`; the unconditional line 431 was overwriting it).

- [ ] **Step 5: Build and verify**

```bash
pixi run check-full
pixi run -e test test-pfunit 2>&1 | tail -10
```

Expected: 5/5 byte-identical, zero pFUnit failures.

If a case fails on the geofac fix: the test fixture might exercise `ipos==5` AND have the surface_runoff.geofac set differently — in which case the fix is correct (the ipos==5 value should win) but the fixture's expected output reflects the old bug. Investigate by inspecting which case has `ipos == 5`.

- [ ] **Step 6: Commit**

```bash
git add src/core/variables.f90 src/core/initialize.f90 src/io/toml/config_to_variables.f90
git commit -m "$(cat <<'EOF'
refactor(state): SS-DRST Phase 2 Task 5 — delete 6 drainage globals + geofac fix

variables.f90 sheds the drainage section: qdrd, qdrain, qdra,
drainl, wetper, ztopdislay declarations commented out with
provenance markers (state%drainage%* is the live home for all
drainage-owned data). initialize.f90 drops the matching zero-init
lines.

Inline fix: config_to_variables.f90:431 gates the surface-runoff
geofac write on ipos /= 5 — the ipos==5 path's geofac (set at line
306) was being unconditionally overwritten by surface_runoff%geofac.
Two distinct TOML fields map to one legacy global; the gate
preserves both intended behaviors.

Spec: docs/superpowers/specs/2026-05-10-state-migration-drainage-design.md
EOF
)"
```

---

### Task 6: ADR 0031 finalize + index + verification

**Files:**
- Modify: `docs/adr/0031-state-migration-drainage.md` — finalize Phase 2 outcomes
- Modify: `docs/adr/index.md` — add ADR 0031 entry
- Modify: `docs/superpowers/specs/2026-05-09-state-migration-drainage-discovery.md` — append Phase 1 lessons-learned (Section 3.5 categorization should include compute readers)

- [ ] **Step 1: Update ADR 0031 Consequences**

In `docs/adr/0031-state-migration-drainage.md`, replace `Phase 2 consequences (to be filled after completion): TBD` with:

```markdown
**Phase 2 outcomes (completed YYYY-MM-DD):**

- 6 drainage-owned globals removed from variables.f90 and
  initialize.f90 (qdrd, qdrain, qdra, drainl, wetper, ztopdislay).
  [If `qdrain` was kept for bocodrb's static-geometry read, document
  the exception.]
- Cross-subsystem reader migrations:
  - frozencond.FrozenBounds reads/writes state%drainage%{qdrain,qdra}
  - surfacewater.WLEVBAL reads state%drainage%qdrain
  - waterbalance.integral reads state%drainage%qdrain
  - All divdra call sites pass typed state directly
- geofac adapter ordering bug fixed: ipos=5 path now gates the
  surface_runoff geofac overwrite.
- Tests: 622 pFUnit (no count change in Phase 2 — only consumer
  migrations and global deletions). check-full 5/5 byte-identical
  at every commit.

**Architectural learnings (Phase 2):**

- divdra's assumed-shape signature (Phase 1 Task 4) was a clean
  unblock — Phase 2 just switched call sites without further
  refactoring. The investment paid off.
- The geofac schema duplication (two TOML fields mapping to one
  legacy global) was fixed inline in the migration arc; the deeper
  schema reconciliation (do we need both fields?) is documented
  as future work.
- bocodrb's static-geometry read of qdrain [if applicable] is the
  only remaining tie to the global — see commit XYZ for details.
```

- [ ] **Step 2: Update ADR index**

Append to `docs/adr/index.md` (after ADR 0030's entry):

```markdown
- [ADR 0031 — Drainage state-type migration](0031-state-migration-drainage.html) — Second subsystem migrated off variables.f90 globals into typed `drainage_state_t` aggregated under `swap_state_t`. qdra and qdrain (previously misclassified into surfacewater_state_t per ADR 0030) moved to drainage_state_t. divdra modernized to assumed-shape arrays. 6 owned globals (qdrd, qdrain, qdra, drainl, wetper, ztopdislay) removed from variables.f90. geofac adapter ordering bug fixed inline (ipos=5 path no longer overwritten by surface_runoff field).
```

- [ ] **Step 3: Add Phase 1 lessons-learned addendum to discovery doc**

Append to `docs/superpowers/specs/2026-05-09-state-migration-drainage-discovery.md` Section 3.5:

```markdown
### Phase 1 lessons-learned addendum (added 2026-05-10)

The original Section 3.5 categorization missed compute readers of
`qdrd`. Specifically, `WLEVBAL` and `WBALLEV` in `surfacewater.f90`
read `qdrd` from globals — they are *compute* code, not output, but
they consume drainage-owned state. Phase 1 Task 7's check-full
failure (~0.36-0.91 cm GWL drift) caught this miss; the readers
were migrated as part of that task's recovery.

**Future subsystem-migration discoveries should categorize Section
3.5 readers by intent (compute vs output vs config), not just by
file location.** A "non-home-tree compute" reader is just as
load-bearing as an "output" reader.

The grep pattern that catches compute readers:
\`\`\`bash
grep -rEn "use variables.*\b<owned-var>\b" src/ --include="*.f90" \\
  | grep -v "src/core/variables.f90 ; src/state/ ; src/io/" \\
  | grep -v "src/<home-tree>"
\`\`\`

Anything matching is a compute reader that needs migration before
the dual-write can be dropped.
```

- [ ] **Step 4: Final verification**

```bash
pixi run -e test test-pfunit 2>&1 | tail -10
pixi run check-full
```

Expected: zero pFUnit failures, 5/5 check-full.

```bash
# Confirm zero drainage-owned globals readable from variables.f90:
grep -rEn "use variables.*\b(qdrd|drainl|wetper|ztopdislay|qdra)\b" src/ --include="*.f90" | grep -v "src/core/variables.f90"
```

Should return empty (or just `qdrain` if Step 1's investigation determined that import must stay).

- [ ] **Step 5: Commit**

```bash
git add docs/adr/0031-state-migration-drainage.md docs/adr/index.md \
        docs/superpowers/specs/2026-05-09-state-migration-drainage-discovery.md
git commit -m "$(cat <<'EOF'
docs(adr+spec): SS-DRST Phase 2 Task 6 — ADR 0031 final, discovery Section 3.5 addendum

ADR 0031 Phase 2 outcomes finalized (6 globals removed, geofac fix,
divdra call-site migrations, all consumer reads from state%drainage).
Index entry added.

Discovery doc gains a Phase 1 lessons-learned addendum to Section
3.5: future subsystem migrations should categorize readers by
intent (compute vs output vs config), not just by file location —
a missed compute reader of qdrd in WLEVBAL/WBALLEV caused a
check-full regression in Phase 1 Task 7.

Drainage migration complete. Branch refactor/surfacewater-state
remains the umbrella for the next subsystem migration.

Spec: docs/superpowers/specs/2026-05-10-state-migration-drainage-design.md
EOF
)"
```

---

## Self-Review Notes

- **Spec coverage:** D6 (l(Madr) m→cm) — done in Phase 1 Task 9 of *surfacewater* arc, not this. D7 (qdrain rule to drainage / ztopdislay) — Task 7 of surfacewater Phase 2 already did the qdrain rule; ztopdislay turned out to have no active write site (Phase 1 Task 5 finding). D9 (geofac fix) — Task 5. Phase 2 phasing — all tasks.
- **Lessons-learned applied:** Pre-flight inventory enumerated all qdrain/qdra readers and writers upfront so scope is locked. The Phase 1 lesson (compute readers, not just output) is captured in the discovery doc addendum.
- **Risks:** Task 4 is the integration gate. If a missed reader exists, check-full catches it. Task 5's `bocodrb` static-geometry exception may force `qdrain` to stay declared — document the exception clearly.
- **Branch policy:** All commits on `refactor/surfacewater-state`. No merge to development until all subsystems land.
