# Surface-Water State Migration — Phase 2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete the surface-water state migration. After Phase 2: zero owned globals remain in `variables.f90` (the 29 surface-water-owned variables plus `fldecdt` are deleted), every reader of those variables across the codebase reads from `state%surfacewater`, the cross-subsystem ownership hazards from spec D6/D7 are resolved, and `swapoutput.f90`'s dead `SurfaceWater(2)` call is removed.

**Architecture:** Threaded `state` argument propagates further into the codebase — every Phase 1 dual-write site has a corresponding cross-subsystem reader that this phase migrates. Each subsystem we touch (drainage, macropore, heat, solute, crop/management_soil, soilhydraulics, waterbalance) gains `state` as an `intent(in)` argument on the routines that read surface-water-owned variables; their bodies switch from `use variables` to `state%surfacewater`. We do **not** migrate those subsystems' own globals — only their reads of surface-water-owned variables. Their own state migrations are future arcs.

**Tech Stack:** Fortran 2008, gfortran, meson + ninja + pixi, pFUnit, check-full byte-identical regression.

**Spec:** `docs/superpowers/specs/2026-05-09-state-migration-surfacewater-design.md` (sections D5, D6, D7, D9 + Phase 2 phasing)

**Branch:** all commits go on `refactor/surfacewater-state` (Phase 1 already on this branch). After Phase 2 verification is green, the branch is ready for review and merge to `development`.

---

## Lessons Learned from Phase 1

These shape how Phase 2 is structured:

1. **Cross-subsystem readers are bigger than the discovery doc cataloged.** The Phase 1 discovery doc enumerated 29 *owned* globals but did not enumerate which *other* subsystems read them. Phase 1 found these only at execution time, surfacing as "kept dual-write because XYZ reads it." Phase 2 starts with that inventory upfront (section "Cross-subsystem reader inventory" below) so scope is locked in. Future subsystem-migration discovery docs should include this section.

2. **Dual-write transitional pattern was the safety net.** Every Phase 1 task ended with check-full byte-identical because compute wrote both state and globals during the transition. Phase 2 inverts: we drop dual-writes one variable family at a time as their readers migrate. check-full at every task is still the gate.

3. **`use variables` shadows ASSOCIATE bare-name aliases.** Phase 1 had to use `sw_*` prefix on ASSOCIATE binds in `WLEVBAL`/`WBALLEV` because `use variables` already brought `wls` etc. into scope. Phase 2 uses the same convention. When migrating cross-subsystem readers, drop the `use variables, only:` symbol and bind via ASSOCIATE without prefix (no shadow if the symbol is no longer imported).

4. **Hidden consumers exist.** Phase 1 discovered `OutputModflow` reading surface-water state via a SAVE-local clone (`state_om`). Watch for similar shims in Phase 2; explicitly grep before declaring scope complete.

5. **Discovery template upgrade.** Phase 2 Task 12 updates the surfacewater discovery doc with a "Section 3.5: external readers of owned globals" — the inventory below — so the template propagates the lesson to future subsystem-migration discoveries.

---

## Cross-subsystem reader inventory

The 16 surface-water-owned globals that still have dual-writes after Phase 1, with the files outside the surface-water home tree that read them:

| Owned global | External reader files | Phase 2 task |
|---|---|---|
| `wls` | `src/drainage/drainage.f90` | Task 3 |
| `swst` | `src/drainage/drainage.f90` | Task 3 |
| `vtair` | `src/soil/soilhydraulics.f90` | Task 6 |
| `cqdra` | `src/drainage/drainage.f90`, `src/soil/waterbalance.f90` | Tasks 3, 7 |
| `ZDraBas` | `src/drainage/drainage.f90`, `src/macropore/macrorate.f90`, `src/macropore/macropore.f90` | Tasks 3, 4 |
| `iqdra` | `src/drainage/drainage.f90`, `src/soil/waterbalance.f90` | Tasks 3, 7 |
| `qdrtot` | `src/heat/frozencond.f90`, `src/soil/waterbalance.f90`, `src/drainage/drainage.f90`, `src/solute/solute.f90` | Tasks 3, 5, 7 |
| `flInitDraBas` | `src/drainage/drainage.f90`, `src/macropore/macropore.f90` | Tasks 3, 4 |
| `imper` | `src/drainage/drainage.f90`, `src/config/swap_config.f90` (validation only — read-only at config time, fine) | Task 3 |
| `sttab` | `src/config/drainage_config.f90` (config-time read, fine) | Task 2 |
| `cqdrain` | `src/soil/waterbalance.f90`, `src/drainage/drainage.f90` | Tasks 3, 7 |
| `cqdrainin` | `src/drainage/drainage.f90`, `src/soil/waterbalance.f90` | Tasks 3, 7 |
| `cqdrainout` | `src/drainage/drainage.f90`, `src/soil/waterbalance.f90` | Tasks 3, 7 |
| `qdra` | `src/heat/frozencond.f90`, `src/drainage/drainage.f90`, `src/soil/waterbalance.f90`, `src/soil/soilhydraulics.f90`, `src/solute/solute.f90` | Tasks 3, 5, 6, 7 |
| `inqdra` | `src/drainage/drainage.f90`, `src/soil/soilgrid.f90`, `src/soil/waterbalance.f90`, `src/crop/management_soil.f90` | Tasks 3, 5, 7 |
| `inqdra_in` | `src/drainage/drainage.f90`, `src/soil/waterbalance.f90` | Tasks 3, 7 |
| `inqdra_out` | `src/drainage/drainage.f90`, `src/soil/waterbalance.f90` | Tasks 3, 7 |

Plus `fldecdt` (separate from the 29; replaced by `intent(out) :: request_smaller_dt`):
| `fldecdt` | `src/core/swap.f90` (3 reads), `src/core/timecontrol.f90` (4 reads), `src/soil/soilhydraulics.f90` (1 read) | Task 1 |

Files NOT in the inventory but already migrated by Phase 1 (no Phase 2 work):
`src/io/swap_csv_output.f90`, `src/io/swapoutput.f90`, `src/io/macroporeoutput.f90`, `src/utils/surfacewaterutils.f90` (Task 2 finishes this).

`src/config/swap_config.f90` and `src/config/drainage_config.f90` references are config-time validation reads of `imper` and `sttab` — these are not runtime reads, do not affect migration ordering, and are addressed inside their respective Phase 2 tasks.

---

## File Structure

**Modified (compute):**
- `src/core/swap.f90` — read `request_smaller_dt`, write to local `fldecdt` (Task 1)
- `src/core/timecontrol.f90` — drop `fldecdt` reads; take a parameter or read from a local-scope-passed flag (Task 1)
- `src/utils/surfacewaterutils.f90` — `wlevst`, `swstlev` take `state` (Task 2)
- `src/drainage/drainage.f90` — `bocodre` and any other surface-water consumers take `state`; relocate `qdrain` zero rule from surfacewater.f90 (Task 3)
- `src/macropore/macrorate.f90`, `src/macropore/macropore.f90` — take `state` (Task 4)
- `src/heat/frozencond.f90`, `src/solute/solute.f90`, `src/soil/soilgrid.f90`, `src/crop/management_soil.f90` — take `state` (Task 5)
- `src/soil/soilhydraulics.f90` — drop `vtair` global read; switch to `state%surfacewater%vtair` (Task 6) — note: `state` already plumbed through soilhydraulics by Phase 1 Task 5, just adding the field reference
- `src/soil/waterbalance.f90` — drop dual-write of surface-water-owned variables in `integral`; reads still supported via state (Task 7)

**Modified (other):**
- `src/drainage/surfacewater.f90` — drop the 16 remaining post-ASSOCIATE dual-writes from `WLEVBAL`/`WBALLEV` after Tasks 3–6 land their consumers (Task 7); remove `qdrain` zero rule (Task 3 moves it to drainage); finalize ASSOCIATE blocks
- `src/drainage/surfacewater_init.f90` — drop remaining dual-writes (sttab, etc.) once consumers migrate (Task 2)
- `src/io/swapoutput.f90` — remove the dead `SurfaceWater(2)` sensitivity loop (Task 10); `OutputModflow` `state_om` re-init scheme adjusted (Task 8)
- `src/config/drainage_config.f90` — `l_unit` field or m→cm at TOML load (Task 9)
- `src/io/toml/read_drainage_toml.f90` — m→cm at load (Task 9)
- `src/core/variables.f90` — delete the 29 owned globals + `fldecdt` (Task 11)
- `src/core/initialize.f90` — drop initialization lines for deleted globals (Task 11)

**Documentation:**
- `docs/superpowers/specs/2026-05-08-state-migration-surfacewater-discovery.md` — add Section 3.5 with the inventory (Task 12)
- `docs/adr/0030-state-migration-surfacewater-pilot.md` — finalize Consequences section (Task 12)
- `docs/adr/index.md` — add ADR 0030 entry (Task 12)

---

### Task 1: fldecdt full migration

`fldecdt` is the only Phase-1 transitional global *not* mirrored in `surfacewater_state_t`. It's a one-bit signal ("decrease timestep"). Phase 1 set `request_smaller_dt = .true.` (intent-out from `SurfaceWater(2)`) AND wrote `fldecdt = .true.` (legacy). Phase 1's call site in `swap.f90` reads `fldecdt` to gate further compute — it does NOT yet read `request_smaller_dt`.

This task: at the call sites, propagate `request_smaller_dt` into a swap_main-local `fldecdt`-equivalent flag, migrate `timecontrol.f90`'s `fldecdt` readers to that flag, then drop the global write from `WLEVBAL`/`WBALLEV` and delete `fldecdt` from `variables.f90`.

**Order matters:** all readers must migrate BEFORE the global is removed. We do this in one task because the surface area is small (3 files, ~10 sites total).

**Files:**
- Modify: `src/core/swap.f90`
- Modify: `src/core/timecontrol.f90`
- Modify: `src/soil/soilhydraulics.f90` (1 fldecdt read)
- Modify: `src/drainage/surfacewater.f90` (drop `fldecdt = .true.` writes from WLEVBAL)
- Modify: `src/core/variables.f90` (delete `logical fldecdt` declaration — moved to Task 11 if dependencies block it)

- [ ] **Step 1: Inventory all fldecdt sites**

```bash
grep -nE "\bfldecdt\b" src/ -r --include="*.f90" 2>/dev/null
```

Confirm the inventory:
- `src/drainage/surfacewater.f90` — write sites in WLEVBAL (typically 2)
- `src/core/swap.f90` — 3 read sites (lines ~286, 290, 294, 297 — guard expressions on subsystem calls)
- `src/core/timecontrol.f90` — 4 sites (lines 62, 70, 71, 584, 593 — reads + a couple of resets)
- `src/soil/soilhydraulics.f90` — 1 read

Document the actual list as found; the line numbers above may have drifted.

- [ ] **Step 2: Add a request_smaller_dt-driven local in swap_main**

In `src/core/swap.f90`, inside `swap_main`, the existing call sequence is:

```fortran
if (.not.fldecdt .and. flSurfaceWater) call SurfaceWater(2, state, request_smaller_dt)
...
if (.not.fldecdt) call SoilWater(2, state)
...
if (.not.fldecdt .and. flSurfaceWater) call SurfaceWater(3, state, request_smaller_dt)
```

After the `SurfaceWater(2)` call, add:
```fortran
if (request_smaller_dt) fldecdt = .true.
```

This bridges the new typed signal to the legacy global. Now `request_smaller_dt` actually drives the gates. Verify the sequencing — the `if (.not.fldecdt)` checks must come AFTER this propagation.

- [ ] **Step 3: Migrate timecontrol.f90's fldecdt readers**

Read `src/core/timecontrol.f90` lines 60–80 and 580–600. The four reads control timestep reduction. Two paths:
- **Path A (recommended):** keep `fldecdt` as a temporary local-or-module variable inside `timecontrol`, fed by `swap_main`'s setting from Step 2. After this task, `fldecdt` exists only inside `swap.f90`/`timecontrol.f90` — not in `variables.f90`.
- **Path B:** thread `request_smaller_dt` as a parameter through `timecontrol` calls. More invasive.

Choose Path A: add `logical, save :: fldecdt = .false.` at module scope in `timecontrol.f90` (or `swap.f90`'s common block), drop the `use variables, only: fldecdt` from both files. The variable is no longer in `variables.f90`'s globals.

If `timecontrol.f90` already takes `state` (Phase 1 may have plumbed it), an even cleaner path: pass `request_smaller_dt` as an argument. Check the current signature.

- [ ] **Step 4: Migrate the soilhydraulics.f90 fldecdt read**

Find the single read in `src/soil/soilhydraulics.f90`. If it's a sub-routine that already takes `state` (Phase 1 plumbed `state` through `soilwater`/`headcalc`/`PONDRUNOFF`), pass the flag through as an additional argument or via the local module-level `fldecdt`. Match Path A from Step 3.

- [ ] **Step 5: Drop fldecdt = .true. write from WLEVBAL**

In `src/drainage/surfacewater.f90`, delete the legacy `fldecdt = .true.` lines paired with `request_smaller_dt = .true.`. Keep the `request_smaller_dt = .true.` writes — they are now load-bearing.

- [ ] **Step 6: Delete fldecdt from variables.f90**

In `src/core/variables.f90`, remove the `logical fldecdt` declaration (typically line ~93). Also drop any initialization in `initialize.f90`.

- [ ] **Step 7: Build and run check-full**

Run: `pixi run check-full`
Expected: `Results: 5 passed, 0 failed`. Byte-identical because the timestep-reduction signal is preserved end-to-end via the new path.

If any case fails, the new path doesn't reach a reader that the old path did. Find the missing migration via `grep -rn "fldecdt"` (now returns no globals — but might find a stale comment or dropped read).

- [ ] **Step 8: Run pFUnit**

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: zero failures, count = 616.

- [ ] **Step 9: Commit**

```bash
git add src/core/swap.f90 src/core/timecontrol.f90 src/soil/soilhydraulics.f90 \
        src/drainage/surfacewater.f90 src/core/variables.f90 src/core/initialize.f90
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWST Phase 2 Task 1 — fldecdt migrated end-to-end

request_smaller_dt (intent(out) on SurfaceWater) is now propagated
into a local fldecdt-equivalent flag in swap_main, consumed by
timecontrol.f90 and soilhydraulics.f90. Legacy global fldecdt
deleted from variables.f90.

Spec: docs/superpowers/specs/2026-05-09-state-migration-surfacewater-design.md
EOF
)"
```

---

### Task 2: surfacewaterutils.f90 (`wlevst`, `swstlev`) take state

These are the only surface-water-home-file functions still on globals — they read `sttab` (and possibly `wls`/`swst`) via module-level `use variables`. Migrate to take `state` and read from `state%surfacewater%sttab`.

`drainage_config.f90` reads `sttab` at config-validation time only (it's not a runtime read; it's during validate or finalize). That doesn't need state — sttab is computed by `surfacewater_init` AFTER config validation runs, so the config-time read is reading the un-initialized sttab (which is fine — it's just reading the array's allocation status or default values). Verify by inspection; if it's a real value-dependent read at config time, escalate.

**Files:**
- Modify: `src/utils/surfacewaterutils.f90` — `wlevst`, `swstlev` take `state`
- Modify: callers of `wlevst`/`swstlev` — `surfacewater_init.f90`, `surfacewater.f90`
- Modify: `src/config/drainage_config.f90` — verify the `sttab` read; either confirm it's safe (no migration) or migrate (escalate to Task 9 territory)
- Modify: `src/drainage/surfacewater_init.f90` and `src/drainage/surfacewater.f90` — drop the `sttab` global dual-write once `wlevst`/`swstlev` no longer read globally

- [ ] **Step 1: Locate wlevst and swstlev callers**

```bash
grep -rn "\b\(wlevst\|swstlev\)\b" src/ --include="*.f90"
```

- [ ] **Step 2: Add state argument to both functions**

```fortran
function wlevst(state, swstn) result(wlev)
   use swap_state_mod, only: swap_state_t
   type(swap_state_t), intent(in) :: state
   real(real64),       intent(in) :: swstn
   real(real64) :: wlev
   ! ...body uses state%surfacewater%sttab via ASSOCIATE...
end function

function swstlev(state, wlev) result(swstn)
   ! similar
end function
```

Drop `sttab` from `use variables, only:` clause if it was there.

- [ ] **Step 3: Update callers**

Pass `state` to every `wlevst` / `swstlev` call site. They're called from `surfacewater.f90` (in WLEVBAL/WBALLEV — the ASSOCIATE block context already has state), `surfacewater_init.f90`, and elsewhere. Update each.

- [ ] **Step 4: Drop sttab dual-write**

In `surfacewater_init.f90`, remove the legacy global `sttab(i, j) = ...` writes (keep only `state%surfacewater%sttab(i, j) = ...`). Trim `use variables, only:` to remove `sttab` if no longer referenced.

- [ ] **Step 5: Verify drainage_config.f90's sttab read**

Read `src/config/drainage_config.f90` for `sttab` references. If it's a config-validate call that reads sttab as state (post-init), that's a bug in the codebase. If it's a config-time validation that doesn't depend on init having run, it's fine. Document the finding.

- [ ] **Step 6: Build and run**

```bash
pixi run check-full
pixi run -e test test-pfunit 2>&1 | tail -10
```

Expected: 5/5 byte-identical, zero pFUnit failures.

- [ ] **Step 7: Commit**

```bash
git add src/utils/surfacewaterutils.f90 src/drainage/surfacewater.f90 \
        src/drainage/surfacewater_init.f90 src/config/drainage_config.f90
git commit -m "refactor(state): SS-SWST Phase 2 Task 2 — wlevst/swstlev take state, drop sttab global"
```

---

### Task 3: drainage.f90 (`bocodre` + others) take state; relocate qdrain rule

`drainage.f90` is the largest cross-subsystem reader. It reads `wls`, `swst`, `cqdra`, `ZDraBas`, `iqdra`, `qdrtot`, `flInitDraBas`, `imper`, `cqdrain`, `cqdrainin`, `cqdrainout`, `qdra`, `inqdra`, `inqdra_in`, `inqdra_out` — 15 of the surface-water-owned variables.

This task threads `state` into `drainage.f90`'s public entry point(s) (likely `Drainage()` or `bocodre()`), updates internal references to use `state%surfacewater`, drops those reads from `use variables, only:`, and **moves the `qdrain(level) = 0.0` zeroing rule (when `gwl > 998`) from surfacewater.f90 into bocodre** per spec D7.

**Files:**
- Modify: `src/drainage/drainage.f90`
- Modify: `src/drainage/surfacewater.f90` — remove the qdrain zero rule (lines ~134-138, see the Phase 1 cleanup comment); drop the now-unused dual-writes that drainage.f90 was the consumer for
- Modify: `src/core/swap.f90` — `Drainage()` call site needs to pass `state`

- [ ] **Step 1: Locate drainage.f90 surfacewater reads**

```bash
grep -nE "\b(wls|swst|cqdra|ZDraBas|iqdra|qdrtot|flInitDraBas|imper|cqdrain|cqdrainin|cqdrainout|qdra|inqdra|inqdra_in|inqdra_out)\b" src/drainage/drainage.f90
```

Tally the read sites and identify the subroutines that need `state`.

- [ ] **Step 2: Add state to drainage.f90 subroutines**

Top-level entry point `Drainage()` (or whatever it's called — discover via the call from `swap.f90`) gains `state` as `intent(inout)`. Internal subroutines (`bocodre`, `divdra`, etc.) that read surface-water state get `state` as `intent(inout)` too. Body uses ASSOCIATE without prefix:

```fortran
associate(wls          => state%surfacewater%wls,          &
          swst         => state%surfacewater%swst,         &
          cqdra        => state%surfacewater%cqdra,        &
          ZDraBas      => state%surfacewater%ZDraBas,      &
          iqdra        => state%surfacewater%iqdra,        &
          ! ... etc
          )
   ! body unchanged
end associate
```

The bare-name aliasing works here because we **drop those names from `use variables, only:`** in the same edit. No shadow.

- [ ] **Step 3: Move qdrain zero rule into bocodre**

Inside `bocodre`, at the top (or wherever the loop body begins), add:

```fortran
! Surface water above wet limit: drainage halts (relocated from
! surfacewater.f90 per ADR 0030 / spec D7).
if (gwl > 998.0d0) then
   do level = 1, nrlevs
      qdrain(level) = 0.0d0
   end do
   return  ! or break, depending on bocodre's structure
end if
```

`gwl` is the relevant condition; `nrlevs` and `qdrain(level)` are local-or-imported names in `drainage.f90`'s scope. Verify that this rule, placed at the top of `bocodre`, has the same observable effect as the previous in-place rule in `SurfaceWater(2)`.

In `surfacewater.f90`, **delete** the original rule (lines ~133-139 in the current file — find them via the Phase 1 cleanup comment marker).

- [ ] **Step 4: Drop surface-water dual-writes that drainage.f90 was the only consumer for**

For each of the 15 variables in the inventory: if drainage.f90 was the ONLY external reader (per the inventory table at the top of this plan), drop the dual-write of that variable from `WLEVBAL`/`WBALLEV` post-ASSOCIATE blocks now. If it has other readers (e.g., `cqdra` is also read by `waterbalance.f90`), keep the dual-write — Task 7 drops it.

Drainage-only readers to drop dual-writes for in this task: `wls`, `swst`, `imper` (also read by `swap_config.f90` validation, but that's config-time only).

Multi-reader vars (keep dual-write):
- `cqdra`, `iqdra`, `qdrtot`, `cqdrain`, `cqdrainin`, `cqdrainout` — also read by `waterbalance.f90` (Task 7)
- `qdra`, `inqdra` — also read by other subsystems (Tasks 5, 6, 7)
- `ZDraBas`, `flInitDraBas` — also read by `macropore.f90`/`macrorate.f90` (Task 4)

- [ ] **Step 5: Update Drainage() caller in swap.f90**

```fortran
if (fldrain) call Drainage(state)
```

- [ ] **Step 6: Build, check-full, pFUnit**

```bash
pixi run check-full
pixi run -e test test-pfunit 2>&1 | tail -10
```

Expected: 5/5 byte-identical (the qdrain rule relocation must produce identical output — same condition, same effect, just executed from `bocodre` instead of in-place in `SurfaceWater(2)`), zero pFUnit failures.

- [ ] **Step 7: Commit**

```bash
git add src/drainage/drainage.f90 src/drainage/surfacewater.f90 src/core/swap.f90
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWST Phase 2 Task 3 — drainage.f90 takes state, qdrain rule relocated

bocodre and the rest of drainage.f90 now read surface-water-owned
fields from state%surfacewater. The 'qdrain(:) = 0 when gwl > 998'
rule (spec D7) moves from surfacewater.f90:134-139 to the top of
bocodre. wls/swst/imper dual-writes dropped (drainage was the sole
consumer); other multi-reader vars (cqdra, qdrtot, etc.) keep
dual-write until Task 7.

Spec: docs/superpowers/specs/2026-05-09-state-migration-surfacewater-design.md
EOF
)"
```

---

### Task 4: macropore subsystem (macrorate.f90, macropore.f90) takes state

These read `ZDraBas` and `flInitDraBas`. Smaller scope than drainage.

**Files:**
- Modify: `src/macropore/macrorate.f90`, `src/macropore/macropore.f90`
- Modify: callers (likely in `swap.f90`)
- Modify: `src/drainage/surfacewater.f90` — drop `ZDraBas`, `flInitDraBas` dual-writes (no remaining external readers after this task)

- [ ] **Step 1: Locate the reads**

```bash
grep -nE "\b(ZDraBas|flInitDraBas)\b" src/macropore/
```

- [ ] **Step 2: Add state to the macropore subroutines that read these**

Public entry points (find via `grep -rn "call <macropore-routine-name>" src/`) get `state` as `intent(inout)`. Internal references switch to `state%surfacewater%ZDraBas`, `state%surfacewater%flInitDraBas` via ASSOCIATE; drop these symbols from `use variables, only:`.

- [ ] **Step 3: Update callers**

Trace from `swap.f90` down. Add `state` to each call.

- [ ] **Step 4: Drop ZDraBas, flInitDraBas dual-writes from surfacewater.f90**

In `WLEVBAL` post-ASSOCIATE block, remove `ZDraBas = state%surfacewater%ZDraBas` and `flInitDraBas = state%surfacewater%flInitDraBas`.

- [ ] **Step 5: Build, check-full, pFUnit**

Expected: 5/5 byte-identical, zero pFUnit failures.

- [ ] **Step 6: Commit**

```bash
git add src/macropore/macrorate.f90 src/macropore/macropore.f90 src/drainage/surfacewater.f90 src/core/swap.f90
git commit -m "refactor(state): SS-SWST Phase 2 Task 4 — macropore takes state for ZDraBas/flInitDraBas"
```

---

### Task 5: heat / solute / soilgrid / management_soil take state

Smaller cross-subsystem readers, lumped into one task because each is a single-call-site change.

**Files:**
- `src/heat/frozencond.f90` — reads `qdrtot`, `qdra`
- `src/solute/solute.f90` — reads `qdrtot`, `qdra`
- `src/soil/soilgrid.f90` — reads `inqdra`
- `src/crop/management_soil.f90` — reads `inqdra`

- [ ] **Step 1: Locate read sites in each file**

```bash
grep -nE "\b(qdrtot|qdra|inqdra)\b" src/heat/frozencond.f90 src/solute/solute.f90 \
                                       src/soil/soilgrid.f90 src/crop/management_soil.f90
```

- [ ] **Step 2: For each file, identify the public subroutine entry point**

Find via `grep -rn "call <routine>" src/` for each.

- [ ] **Step 3: Add state to the entry point and read from state%surfacewater**

Pattern is the same as Task 4: `state` gets added as `intent(inout)` (or `intent(in)` if the routine doesn't write surfacewater state — likely intent(in) for these readers); ASSOCIATE binds the field; drop the symbol from `use variables, only:`.

- [ ] **Step 4: Update callers**

Trace from `swap.f90` (or wherever the routine is called from) and propagate `state`.

- [ ] **Step 5: Build, check-full, pFUnit**

Expected: 5/5 byte-identical.

- [ ] **Step 6: Commit**

```bash
git add src/heat/frozencond.f90 src/solute/solute.f90 src/soil/soilgrid.f90 \
        src/crop/management_soil.f90 src/core/swap.f90
git commit -m "refactor(state): SS-SWST Phase 2 Task 5 — heat/solute/soilgrid/crop read surfacewater from state"
```

---

### Task 6: soilhydraulics.f90 vtair migration

`soilhydraulics.f90` already takes `state` (Phase 1 plumbed it through `soilwater`/`headcalc`/`PONDRUNOFF`). One field still goes through `use variables`: `vtair`. Switch it to `state%surfacewater%vtair` and drop from the `use` clause. Same with `qdra` if it appears there.

**Files:**
- Modify: `src/soil/soilhydraulics.f90`
- Modify: `src/drainage/surfacewater.f90` — drop `vtair` dual-write (no other external readers)

- [ ] **Step 1: Find vtair (and qdra, if any) reads in soilhydraulics.f90**

```bash
grep -nE "\b(vtair|qdra)\b" src/soil/soilhydraulics.f90
```

- [ ] **Step 2: Switch reads to state**

ASSOCIATE block (or direct `state%surfacewater%vtair` if the read is single-site).

- [ ] **Step 3: Drop vtair, qdra (if any) from `use variables, only:` clause**

- [ ] **Step 4: Drop vtair dual-write from surfacewater.f90**

The `vtair = state%surfacewater%vtair` line (if it exists) in WLEVBAL post-ASSOCIATE.

- [ ] **Step 5: Build, check-full, pFUnit**

- [ ] **Step 6: Commit**

```bash
git add src/soil/soilhydraulics.f90 src/drainage/surfacewater.f90
git commit -m "refactor(state): SS-SWST Phase 2 Task 6 — soilhydraulics reads vtair from state"
```

---

### Task 7: waterbalance.f90 integral — drop dual-write

`integral` was modified in Phase 1's Task 5-fix to dual-write `iqdra`/`cqdra`/`cqdrain`/`inqdra`/`inqdra_in`/`inqdra_out`/`cqdrainin`/`cqdrainout`/`qdrtot`/`qdra` to state. After Tasks 3+5+6 land their consumers, the legacy globals are no longer read. Drop the dual-write from `integral`.

This is a one-file task, but it's the integration test for everything before — if Tasks 3/5/6 missed a consumer, this step's check-full will fail.

**Files:**
- Modify: `src/soil/waterbalance.f90`

- [ ] **Step 1: Find all dual-write blocks in integral**

```bash
grep -nE "(state%surfacewater%[a-zA-Z_]+ = |cqdra |iqdra |cqdrain |inqdra |qdrtot |qdra )" src/soil/waterbalance.f90 | head -40
```

Identify each `<global> = ...` line that has a paired `state%surfacewater%<global> = ...` directly above or below.

- [ ] **Step 2: Drop the legacy global writes**

For each pair, remove the `<global> = ...` line, keep the `state%surfacewater%<global> = ...` line.

- [ ] **Step 3: Trim use variables, only:**

Remove the surface-water-owned symbols from `use variables, only:` in `integral` (and any helpers).

- [ ] **Step 4: Build, check-full, pFUnit**

Expected: 5/5 byte-identical. **If a case fails, it means a consumer in Tasks 3/5/6 was missed** — find it via grep:

```bash
grep -rn "use variables.*\b(iqdra|cqdra|cqdrain|inqdra|inqdra_in|inqdra_out|cqdrainin|cqdrainout|qdrtot|qdra)\b" src/ --include="*.f90" | grep -v "src/core/variables.f90\|src/state/\|src/drainage/surfacewater\|src/io/\|src/utils/surfacewaterutils\|src/soil/waterbalance"
```

If it returns anything, those are unmigrated readers. Either migrate them in this task or roll back to keep the dual-write for that variable.

- [ ] **Step 5: Commit**

```bash
git add src/soil/waterbalance.f90
git commit -m "refactor(state): SS-SWST Phase 2 Task 7 — waterbalance.integral drops dual-write"
```

---

### Task 8: OutputModflow re-evaluation

`OutputModflow` in `swapoutput.f90` has a SAVE-local `state_om` (Phase 1 transitional shim). It calls `surfacewater_init(state_om)` once, then runs its own `SurfaceWater(2/3, state_om)` mini-loop for diagnostic sensitivity output. With Phase 2 dropping global dual-writes, `state_om` and `state` (main) start identical (both initialized by `surfacewater_init`) but evolve independently. The question: is that the correct semantic, or should `state_om` be reset from `state` (main) at the start of each mini-loop iteration?

Investigation needed: read the `OutputModflow` body (in `src/io/swapoutput.f90`, find via `grep -n "subroutine OutputModflow"`) and confirm whether the diagnostic loop is intended to:
- **(A) Run a fresh mini-simulation** from the same initial state — the current `state_om` SAVE-local approach is correct.
- **(B) Run a one-step what-if** from main state — `state_om` should be re-copied from `state` at every diagnostic iteration.

Decide by reading the comments and the loop body.

**Files:**
- Modify: `src/io/swapoutput.f90` — `OutputModflow`

- [ ] **Step 1: Read OutputModflow body**

```bash
grep -n "subroutine OutputModflow" src/io/swapoutput.f90
sed -n '<line>,<line+150>p' src/io/swapoutput.f90
```

Form a hypothesis about the intended semantic.

- [ ] **Step 2: Implement chosen approach**

If (A): no change needed — current behavior is correct, just confirm with a comment.

If (B): at the top of each mini-loop iteration, `state_om%surfacewater = state%surfacewater` (deep-copy via Fortran intrinsic assignment for derived types). `state` must be threaded into `OutputModflow` as an additional argument. Update its caller in `swap.f90`.

- [ ] **Step 3: Build, check-full, pFUnit**

Note: check-full may not exercise OutputModflow if `swmodflow` is off in the regression cases. If so, pFUnit must cover this — write a minimal test if no existing test exercises the path.

- [ ] **Step 4: Commit**

```bash
git add src/io/swapoutput.f90 src/core/swap.f90
git commit -m "refactor(state): SS-SWST Phase 2 Task 8 — OutputModflow state_om re-evaluation"
```

---

### Task 9: l(Madr) m→cm conversion to drainage config load

Spec D6 hazard. Currently `surfacewater_init` mutates `l(Madr)` in-place from m → cm. This is a drainage-owned variable — surface-water shouldn't be mutating it.

**Files:**
- Modify: `src/io/toml/read_drainage_toml.f90` — convert `l` from m to cm at TOML-read time
- Modify: `src/config/drainage_config.f90` — declare that `l` is in cm (update doc / unit comment)
- Modify: `src/drainage/surfacewater_init.f90` — remove the m→cm conversion

- [ ] **Step 1: Locate the conversion**

```bash
grep -n "l(.*) = " src/drainage/surfacewater_init.f90 | grep -i "100\|cm\|m "
```

Should find a line like `l(i) = l(i) * 100.0_real64` or similar.

- [ ] **Step 2: Find where l is loaded from TOML**

```bash
grep -n "config%drainage.*l\b\|drainage%l\b" src/io/toml/read_drainage_toml.f90
```

- [ ] **Step 3: Move the conversion to TOML-load**

In `read_drainage_toml.f90`, immediately after reading `l` from TOML, convert m → cm. Annotate with a comment that the legacy `surfacewater_init` conversion is removed.

- [ ] **Step 4: Update drainage_config.f90 documentation**

Update the inline comment for `l` in the config struct to specify "cm (converted from TOML's m at load time)".

- [ ] **Step 5: Remove the conversion from surfacewater_init.f90**

Delete the in-place mutation lines.

- [ ] **Step 6: Build, check-full, pFUnit**

Expected: 5/5 byte-identical. Conversion happens once at load instead of once at init — same value, different timing. Tests must cover the TOML reader's conversion (write a fixture-based test if missing).

- [ ] **Step 7: Commit**

```bash
git add src/io/toml/read_drainage_toml.f90 src/config/drainage_config.f90 src/drainage/surfacewater_init.f90 tests/unit/io/toml/...
git commit -m "refactor(config): SS-SWST Phase 2 Task 9 — l(Madr) m→cm at TOML load (spec D6)"
```

---

### Task 10: Remove dead SurfaceWater(2) call from swapoutput.f90

Spec D9 hazard 3. Inside `swapoutput.f90`, around lines 3639/3651 (per discovery doc), there's a `SurfaceWater(2, ...)` call inside a sensitivity loop. Per the Phase 1 spec, this is dead code in the modern flow — confirm by inspection and delete.

**Files:**
- Modify: `src/io/swapoutput.f90`

- [ ] **Step 1: Locate the call**

```bash
grep -nE "call (SurfaceWater|surfacewater)" src/io/swapoutput.f90
```

The discovery doc said lines 3639, 3651. Find the actual lines.

- [ ] **Step 2: Confirm it's dead**

Read the surrounding code. Verify the loop is dead in the modern flow (not reachable by any code path that's not stub-errored or removed). The discovery doc said "leftover from interactive-debug code" — verify this claim against the current code.

If it IS reachable and live, this is more than a deletion — it's a real consumer that needs migration. Escalate.

- [ ] **Step 3: Delete the loop or call**

If purely dead: delete. If part of a live path that's just gated off: depends — may be a no-op deletion, may need preserving the gate but changing the body.

- [ ] **Step 4: Build, check-full, pFUnit**

Expected: 5/5 byte-identical (deleting dead code can't change observable behavior).

- [ ] **Step 5: Commit**

```bash
git add src/io/swapoutput.f90
git commit -m "refactor(state): SS-SWST Phase 2 Task 10 — remove dead SurfaceWater(2) call from swapoutput (spec D9)"
```

---

### Task 11: Remove the 29 owned globals from variables.f90

The cleanup. After Tasks 1–10, no compute reads or writes any of the 29 surface-water-owned variables from `variables.f90`. Confirm via grep, then delete the declarations. `fldecdt` was already deleted in Task 1.

**Files:**
- Modify: `src/core/variables.f90`
- Modify: `src/core/initialize.f90` — drop initialization lines for the deleted globals

- [ ] **Step 1: Confirm zero readers/writers**

```bash
for v in wls wlstar swst swstini hwlman vtair wlsold cqdrd cwsupp cwout cqdra ZDraBas iqdra qdrtot overfl flInitDraBas imper numadj wlsbak sttab cqdrain cqdrainin cqdrainout qdra inqdra inqdra_in inqdra_out; do
  count=$(grep -rEn "\b$v\b" src/ --include="*.f90" 2>/dev/null \
    | grep -v "src/core/variables.f90" \
    | grep -v "src/core/initialize.f90" \
    | grep -v "src/state/" \
    | wc -l)
  echo "$count  $v"
done | sort -rn
```

Expected: every count is 0. If any are nonzero, that variable still has a reader — go back to whichever Task should have migrated it.

- [ ] **Step 2: Delete declarations from variables.f90**

Find the declaration lines (per discovery doc: lines 1255–1273 area, plus `vtair` at line 960). Delete each.

- [ ] **Step 3: Drop initialization lines from initialize.f90**

Find via `grep -n "<varname>" src/core/initialize.f90` for each. Delete or comment out.

- [ ] **Step 4: Build, check-full, pFUnit**

Expected: 5/5 byte-identical, zero pFUnit failures. **This is the moment of truth** — if the build fails, a consumer was missed. Restore the declaration for that variable, find the consumer, migrate it.

- [ ] **Step 5: Variables.f90 LoC stats**

```bash
wc -l src/core/variables.f90
```

Compare to the pre-Phase-1 count (1314). Phase 2 should reduce by ~25–30 lines (29 deletions plus a couple of declaration-block trims).

- [ ] **Step 6: Commit**

```bash
git add src/core/variables.f90 src/core/initialize.f90
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWST Phase 2 Task 11 — remove 29 surface-water-owned globals

variables.f90 sheds the surface-water section. Every reader and
writer of the 29 owned variables now goes through state%surfacewater.
This is the proof that Phase 1+2 of the surface-water state migration
is complete: typed state is the sole authoritative location.

Spec: docs/superpowers/specs/2026-05-09-state-migration-surfacewater-design.md
EOF
)"
```

---

### Task 12: ADR 0030 finalize + discovery template upgrade + final verification

Last task. Lock in the documentation, update the playbook artifact for future subsystems, run the final integration check.

**Files:**
- Modify: `docs/adr/0030-state-migration-surfacewater-pilot.md` — finalize Consequences section
- Modify: `docs/adr/index.md` — add ADR 0030 entry
- Modify: `docs/superpowers/specs/2026-05-08-state-migration-surfacewater-discovery.md` — add Section 3.5 (external readers of owned globals)

- [ ] **Step 1: Update ADR 0030 Consequences**

Replace the "Phase 2 consequences (to be filled after completion): TBD" with concrete outcomes:
- Variables.f90 LoC delta (e.g., 1314 → 1283).
- Number of cross-subsystem readers migrated.
- Names of consequence files (the 13 from the inventory at the top of this plan).
- Any architectural learnings (e.g., the qdrain rule relocation pattern, OutputModflow handling).

- [ ] **Step 2: Add ADR 0030 to docs/adr/index.md**

Match the house style (per ADR 0028's index entry):

```markdown
- [ADR 0030 — Surface-water state-type migration (pilot)](0030-state-migration-surfacewater-pilot.html) — First subsystem migrated off variables.f90 globals: surfacewater_state_t holds the 27 runtime fields, swap_state_t aggregates per-subsystem state. Argument-threaded from swap_main; ASSOCIATE in compute. fldecdt replaced by intent(out) request_smaller_dt. qdrain rule relocated to drainage's bocodre. l(Madr) m→cm moved to drainage TOML load. Establishes the migration playbook for subsequent subsystems.
```

- [ ] **Step 3: Add Section 3.5 to surfacewater discovery doc**

Insert after Section 3 (Borrowed globals) the inventory from this plan as a new section:

```markdown
## 3.5 External readers of owned globals (cross-subsystem reader inventory)

Variables this subsystem WRITES that other subsystems READ. Discovered during Phase 1 of the migration; not present in the original discovery (lessons-learned for future discoveries).

| Owned variable | External reader files |
|---|---|
| (table from this plan) |
```

This makes the discovery template more complete; future subsystem migrations should populate this section during discovery, not during execution.

- [ ] **Step 4: Final verification**

```bash
pixi run -e test test-pfunit 2>&1 | tail -10
pixi run check-full
```

Expected:
- pFUnit: zero failures.
- check-full: 5/5 byte-identical.

- [ ] **Step 5: Confirm zero owned globals remain in variables.f90**

```bash
grep -nE "^\s+(real|integer|logical|character).*\b(wls|wlstar|swst|swstini|hwlman|vtair|wlsold|cqdrd|cwsupp|cwout|cqdra|ZDraBas|iqdra|qdrtot|overfl|flInitDraBas|imper|numadj|wlsbak|sttab|cqdrain|cqdrainin|cqdrainout|qdra|inqdra|inqdra_in|inqdra_out|fldecdt)\b" src/core/variables.f90
```

Expected: empty output. (Some symbol names may collide with other subsystems' identically-named variables — investigate any matches case-by-case.)

- [ ] **Step 6: Commit**

```bash
git add docs/adr/0030-state-migration-surfacewater-pilot.md docs/adr/index.md \
        docs/superpowers/specs/2026-05-08-state-migration-surfacewater-discovery.md
git commit -m "$(cat <<'EOF'
docs(adr+spec): SS-SWST Phase 2 Task 12 — ADR 0030 final, discovery template upgrade

ADR 0030 Consequences section finalized with concrete Phase 2
outcomes (variables.f90 LoC delta, cross-subsystem readers
migrated, architectural learnings).

Surfacewater discovery doc gains Section 3.5 (external readers of
owned globals) — populated retroactively from Phase 1 findings.
This becomes part of the migration playbook for subsequent
subsystems: future discoveries enumerate this section upfront, not
at execution time.

Branch refactor/surfacewater-state is now ready for review and
merge to development.

Spec: docs/superpowers/specs/2026-05-09-state-migration-surfacewater-design.md
EOF
)"
```

---

## Self-Review Notes

**Spec coverage:**
- D5 (fldecdt as intent(out)) → Task 1 fully migrates and removes the global.
- D6 (l(Madr) m→cm) → Task 9.
- D7 (qdrain rule to drainage) → Task 3.
- D9 (output coupling: dead SurfaceWater(2) call removal) → Task 10. Output reads were Phase 1's job.
- Phase 2 phasing: "Remove fldecdt, the 29 owned globals" → Tasks 1 + 11.
- Cross-subsystem reader migrations (NOT in spec, lessons-learned addition) → Tasks 2-7.
- OutputModflow re-evaluation (lessons-learned addition) → Task 8.
- ADR finalize → Task 12.

**Lessons-learned applied:**
- Cross-subsystem reader inventory enumerated upfront (avoid Phase-1-style mid-execution surprises).
- Dual-write drop ordering: per-task, only after the consumer is migrated.
- Discovery template gains Section 3.5.

**Risks:**
- Task 7 (drop waterbalance dual-write) is the integration test for Tasks 3 + 5 + 6. If those tasks miss a consumer, Task 7's check-full reveals it. This is the architectural integration gate.
- Task 11 (remove globals) is the second integration gate. If anything was missed, the build fails.
- OutputModflow (Task 8) may not be exercised by check-full. If pFUnit doesn't cover it either, this risk goes uncaught — flag as DONE_WITH_CONCERNS if no test exists.

**Branch policy:** Phase 2 commits go on `refactor/surfacewater-state`. After Task 12, the branch is ready for review and merge to `development`.
