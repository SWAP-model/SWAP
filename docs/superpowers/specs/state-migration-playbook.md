# State Migration Playbook

This document captures the patterns established by the first three subsystem state migrations (surfacewater ADR 0030, drainage ADR 0031, solute ADR 0032) plus the cumulative-reset cohort refactor (ADR 0033). Future subsystem migrations bake these patterns in from the start.

**Status:** living document. Append patterns and gotchas as new migrations surface them.

---

## Phase shape

Every subsystem migration follows the same four-phase shape. Each phase produces a versioned artifact:

1. **Discovery** (read-only fact-finding) → `docs/superpowers/specs/YYYY-MM-DD-state-migration-<subsystem>-discovery.md`
2. **Design** (decisions: state-type fields, public API, hazards resolved) → `docs/superpowers/specs/YYYY-MM-DD-state-migration-<subsystem>-design.md`
3. **Plan** (per-task decomposition for subagent execution) → `docs/superpowers/plans/YYYY-MM-DD-<subsystem>-state-phase{N}.md`
4. **Execution** → commits + ADR

Larger subsystems (>15 owned globals) typically split into Phase 1 (state type + threading + dual-write + drop) and Phase 2 (cross-subsystem reader migration + global removal). Smaller ones can ship in a single overarching plan.

---

## Discovery template — Section 2 categorization (mandatory upfront)

Section 2 of each discovery doc enumerates the subsystem's owned globals (variables in `variables.f90` that the subsystem mutates). **Classify each owned global by reset cadence** — this directly informs the cohort sub-record design in the state type.

| Reset cadence | Description | Cohort target | Examples |
|---|---|---|---|
| **Instantaneous** | Reset every step unconditionally OR computed fresh each step (no `flzero*` gating) | Stays as flat field on `<subsystem>_state_t` | `dt`, `csurf`, `cml(:)` (recomputed daily), `isqbot` |
| **Intermediate** | Reset on `flzerointr` (intermediate-period rollover, typically daily/sub-output) | `<subsystem>_intermediate_t` | `iqdra`, `imsqprec`, `inqdra(:,:)` |
| **Cumulative** | Reset on `flzerocumu` (cumulative-period rollover, typically yearly/output-period) | `<subsystem>_cumulative_t` (or partitioned by activity gate — see below) | `cqdra`, `sqdra`, `samini`, `cqdrain(:)` |

### Cohort partitioning by activity gate (added 2026-05-10 from SS-CRR Phase A correction)

When cumulative fields in one subsystem's state are accumulated under different activity flags, **partition the cohort by flag** — not by name prefix or data type. Surfacewater's cumulative cohort is the worked example:

| Cohort | Fields | Activity gate | Reset owner |
|---|---|---|---|
| `surfacewater_drainage_cumulative_t` | `cqdra`, `cqdrain(:)`, `cqdrainin(:)`, `cqdrainout(:)` | `fldrain` (active under `swdra=1` OR `swdra=2`) | `Drainage()` (and `SurfaceWater(2)` when `swdra=2`) |
| `surfacewater_reservoir_cumulative_t` | `cqdrd`, `cwsupp`, `cwout` | `flSurfaceWater` (active only under `swdra=2`) | `SurfaceWater(2)` only |

The owner of each cohort is the subsystem whose activity flag gates that cohort's accumulation, and only the owner calls `reset()`. The original mistake (Task A5) was a single combined cohort with an inline comment papering over the fact that drainage zeros a strict subset; the type system can express the contract directly via partitioning.

**Discovery checklist:** during Section 2 categorization of cumulative fields, also note the activity flag(s) that gate each field's accumulation paths. If two fields in the same subsystem have different gates, they belong in different cohorts. See ADR 0033 for the full pattern.

The classification grid in Section 2 should look like:

```markdown
| Variable | Type | Cadence | Cohort | Description |
|----------|------|---------|--------|-------------|
| wls      | real | step    | flat   | surface water level |
| iqdra    | real | inter   | intermediate | per-step lateral drainage total |
| cqdra    | real | cumu    | cumulative | cumulative lateral drainage |
```

If a global is reset by *multiple* call sites with non-identical subsets (the surfacewater/drainage finding from ADR 0033 Phase A), flag it in Section 8 (hazards) — the cohort can't enforce subset-reset symmetry; it's a per-call-site policy.

---

## State-type design pattern

Each subsystem state record:

```fortran
type :: <subsystem>_state_t
   ! Instantaneous fields (flat)
   real(real64) :: <field1> = 0.0_real64
   real(real64), allocatable :: <field2>(:)
   ! …

   ! Cohort sub-records — one per reset cadence
   type(<subsystem>_intermediate_t) :: intermediate
   type(<subsystem>_cumulative_t)   :: cumulative
end type
```

Each cohort declares the type-bound `reset()` procedure as its **only** method:

```fortran
type :: <subsystem>_intermediate_t
   real(real64) :: <field> = 0.0_real64
   ! …
contains
   procedure :: reset => <subsystem>_intermediate_reset
end type

! In the module's contains block:
subroutine <subsystem>_intermediate_reset(self)
   class(<subsystem>_intermediate_t), intent(inout) :: self
   self%<field> = 0.0_real64
   ! Use `if (allocated(self%<arr>)) self%<arr> = 0.0_real64` for allocatables.
end subroutine
```

**Direct field writes for accumulation. No setters, no `add_*` methods, no auto-resets.** Only `reset()` is type-bound.

Special non-zero rebases (e.g. `samini = sampro` in solute) stay inline at the flag-gated call site, immediately after the corresponding `reset()` call, with an explanatory comment.

---

## Reset block pattern

Inline reset blocks collapse to two flag-gated `reset()` calls plus any inline rebase:

```fortran
if (flzerointr) call state%<subsystem>%intermediate%reset()
if (flzerocumu) then
   call state%<subsystem>%cumulative%reset()
   ! Optional: inline non-zero rebase, e.g. samini = sampro for solute.
end if
```

### Asymmetric reset hazard (ADR 0033 Phase A finding)

If multiple call sites reset cohort fields and the subsets differ, calling the full cohort `reset()` from every site can corrupt mid-accumulation values at the wrong site. The pattern doesn't enforce full-reset symmetry.

**Resolution:** call the full cohort `reset()` from the canonical owner (typically the latest call in the timestep call chain), and keep element-by-element subset zeroing at intermediate sites with a comment pointing to the ADR section explaining the asymmetry.

---

## Cross-subsystem reader inventory (Section 3.5 — categorization framework)

For each owned global, catalog readers OUTSIDE the home tree by intent:

1. **Output readers** — output routines (e.g., `outdrf`, `outbal`, `set_values`). Most obvious.
2. **Compute readers** — physics routines in OTHER subsystems consuming this subsystem's state (e.g., `frozencond.FrozenBounds` reads drainage's `qdrain`).
3. **Working-buffer reads** — code paths using the legacy global as scratch space.
4. **Init-routine seed reads** — init routines reading legacy globals to populate state at startup.
5. **Call-site argument reads** — reads via routine arguments (e.g., `divdra` taking `qdra` as arg).

The grep templates that catch each category:

```bash
# Categories 1+2 (use variables clauses):
grep -rEn "use variables.*\b<owned-var>\b" src/ --include="*.f90" \
  | grep -v "src/core/variables.f90 ; src/state/ ; src/<home-tree>"

# Categories 3+4+5 (raw symbol references after migration):
grep -rEn "\b<owned-var>\b" src/ --include="*.f90" \
  | grep -v "src/core/variables.f90" \
  | grep -v "src/core/initialize.f90" \
  | grep -v "src/state/" \
  | grep -v "state%" \
  | grep -v "config%"

# ALSO: aliased ASSOCIATE forms — drainage Phase B B3 found `sl => state%solute` blocks
# that the canonical grep missed. Look for any `<alias> => state%<subsystem>%` pattern.
grep -rEn "=> state%<subsystem>%" src/ --include="*.f90"
```

**Verification sequencing:** before deleting a legacy global declaration, run BOTH greps + the alias-form grep. Any non-zero result is an unmigrated reader that must be handled in the same task as the deletion.

---

## Argument threading

Threading `state` from `swap_main` is the entry-point pattern. Subsystem entry points take `state` as `intent(inout)`. ASSOCIATE blocks inside compute bodies bind frequently-used fields:

```fortran
associate(sw => state%surfacewater)
   ! Body uses sw%X — reads and writes work as for the type
end associate
```

When the subsystem also has `use variables, only: …` for borrowed globals (or even bare `use variables`), bare-name aliases inside ASSOCIATE shadow the imports. Use a 2-character prefix to avoid ambiguity:

| Subsystem | Prefix |
|---|---|
| surfacewater | `sw_` |
| drainage | `dr_` |
| solute | `sl_` |

Where there's no ambiguity (i.e., `use variables` is narrowed to non-overlapping symbols), bare-name aliases are cleaner.

---

## Init-routine pattern

Each subsystem has a `<subsystem>_init(state)` subroutine in its home tree that:

1. Allocates per-node / per-level arrays in `state%<subsystem>` from the relevant dimensions (`numnod`, `nrlevs`, etc.).
2. Seeds `state%<subsystem>` from the config-time-populated legacy globals (the two-stage seeding pattern).

Call from `swap_main` once at run init, after `config_to_variables` and before the first per-step compute call. Order matters when a subsystem's state seeds another's input.

---

## Lifecycle invariants

Each migration arc must maintain these invariants:

- **check-full byte-identical at every commit.** Not just at the end of a phase. The dual-write transitional pattern is the safety net; verification at every commit is mandatory.
- **pFUnit zero failures.**
- **No physics changes.** Migration moves data + reset logic; bodies of compute routines should be semantically unchanged.
- **Single commit for the cohort-nesting migration** (Task A3 / B3 pattern). Partial state will not compile.

---

## Common gotchas (lessons learned)

1. **`use Variables` wildcard imports.** Hard to clean up cleanly because narrowing requires enumerating dozens of names. Defer to subsystem reactivation arc when feasible (see solute Task 9 / agetracer.f90).
2. **Compute readers in init routines.** Init routines may read legacy globals to seed state — these reads stay until the subsystem fully retires the global. See drainage_init's `bocodrb`-static-geometry exception (ADR 0031).
3. **MatricFlux-style optional state args.** When a routine is called from many call chains (some pre-migration, some post-migration), accept `state` as `optional intent(in)` so init-time callers don't need to be touched.
4. **Hidden ASSOCIATE shadow names.** Searches for `state%<subsystem>%X` miss aliased forms (`sl => state%solute`, `sw => state%surfacewater`, etc.). Always grep for both.
5. **Mass-balance rebase semantics.** `samini = sampro` in solute is physics, not cohort policy. Stays inline at the call site immediately after `cumulative%reset()`. Cohort `reset()` zeroes `samini` first; the rebase anchors it to the new period's baseline.
6. **Output writeback hazards.** Output routines that mutate state (e.g., year-end `swstini = swst` in `outswb`) are anti-patterns. Move to dedicated callbacks (`<subsystem>_year_reset(state%<subsystem>)`) invoked from `swap_main`.

---

---

## Lessons from heat migration (ADR 0034, added 2026-05-10)

These five lessons generalize beyond heat and apply to future subsystem migrations.

1. **Init-order trap — use `allocated()` guards for state fields consumed before `<subsystem>_init` runs.** In `swap.f90`'s startup sequence, `SoilWater(1, state)` runs before `heat_init`. `SoilWater` resets `state%heat%rfcp(:) = 1.0` every Richards sub-step. Because `heat_init` allocates `rfcp` later, the first Richards call sees an unallocated array and segfaults without the guard. Fix: `if (allocated(state%heat%rfcp)) state%heat%rfcp(:) = 1.0_real64`. **General rule:** during discovery, note the startup call order in `swap.f90`/`swap_main`. Any state field that an earlier-running subsystem reads or writes needs an `allocated()` guard at those call sites.

2. **Non-module routines + state plumbing — use non-optional dummy arrays + use-rename instead of optional+interface.** Passing typed `state` into external (non-module) routines as an `optional` dummy arg requires an explicit interface block at the call site. For routines with large implicit interfaces this is impractical. Alternative: pass a **non-optional** array slice (`tsoil(:)`) as an explicit dummy arg — gfortran resolves this without an explicit interface — and exclude the homonymous global from scope via `use Variables, dummy_X_ => tsoil`. The rename trick cleanly eliminates exactly one global from the namespace without touching unrelated imports. Used for `ArableLandGerm`, `sumttd`, and `grass` in the heat arc.

3. **Output-side parameter sweeps are not writeback hazards — fix is a local scratch array.** A `call <physics>(synthetic_input, owned_field, …)` inside an output routine where `owned_field` is an out-arg looks like a stale-state writeback hazard (playbook gotcha #6). It may actually be a **parameter sweep**: synthetic inputs fed through the physics function to produce a range of output columns (e.g., `heacap` at pF0, pF1, …, pF4.2). In that case the owned field should never have been modified in the first place — the fix is a local scratch array, not relocating the compute to the physics side. Inspect the call's purpose before choosing a fix: if the inputs are synthetic (not the actual current-step state), it is a parameter sweep.

4. **Audit compute auxiliaries for hidden global reads — devries-style.** When migrating a compute routine, also audit the auxiliary subroutines it calls. They may import globals via `use variables, only:` even if the top-level compute was clean. In the heat arc, `devries` (called from `outheapar`'s parameter sweep) read composition arrays from the global namespace. Task 8 extended devries' signature to take those arrays explicitly — eliminating the hidden dependency. **General rule:** after migrating a compute entry point, grep for `use variables` in every helper it calls, including helpers in utility modules (`soilhydraulicsutils`, etc.).

5. **Flat state-type layout when no fields are flag-gated cumulatives.** The cohort sub-record pattern (ADR 0033) is only needed when owned fields accumulate under `flzero*` gates. If all owned fields are instantaneous (recomputed or overwritten every step), a plain flat `*_state_t` is correct — no `intermediate` or `cumulative` sub-records. Heat is the first migration to use this layout. **General rule:** during discovery Section 2 categorization, if every owned global is marked "instantaneous," the cohort pattern does not apply and can be omitted from the design spec.

---

---

## Lessons from boundary migration (ADR 0035, added 2026-05-10)

Boundary was migration #5 and the first coupling-surface arc of the soil-water decomposition. Eight lessons generalize to future arcs.

1. **Coupling-surface decomposition — decompose large subsystems by coupling surface, not by phase.** When a subsystem owns ~100 globals (soil-water), one arc is too large. Decompose by the surface each arc touches: boundary fields first, then crop-uptake, then atmosphere, then the Richards-interior core. Each arc carves only the fields crossing its surface; the final core arc handles internal state. This minimizes blast radius per commit and produces clean ADR narratives. The anti-pattern is splitting by implementation phase (Phase 0 / Phase 1 / Phase 2) within one giant arc — that keeps blast radius high regardless.

2. **Defer fields with un-plumbed co-writers.** If a field has co-writers outside the home tree that have no `state` argument yet (here: `pond` co-written by `tillage.f90`; `gwl` co-written by `calcgwl()` in `waterbalance.f90`), defer the field rather than expanding the arc to plumb those routines. Each arc stays focused; the co-writer plumbing happens naturally in the arc that owns those routines. Trying to migrate `pond` here would have forced tillage state-plumbing, an init-order guard, and an ASSOCIATE block into a boundary-focused commit sequence.

3. **Pre-existing state-arg windfall — check before scoping.** When previous arcs (here: SS-HEAT Task 9) already plumbed `state` into the home-tree entry points, the current arc reduces to field-carve + reader cutover with no signature surgery. Before scoping discovery, grep all home-tree entry points for `intent(in)` or `intent(inout)` state args. A windfall can cut estimated scope by 30–50%.

4. **Pure-reader files can drop the legacy `use variables, only:` import in Phase 2.** For files that only READ a field and never write it, the global can be removed from the `use variables, only:` clause as soon as reads are migrated to `state%X`. There is no Phase 2.7 dual-write to drop. Mixed reader-writer files (co-write sites) must retain the import until Phase 2.7 cleanup. Identifying pure-reader vs. mixed files during discovery simplifies the Phase 2 task breakdown.

5. **Compile-driven Phase 2.7 is now an expected step, not a surprise.** After dropping dual-writes, the compiler surfaces remaining global references that the grep-based read-only discovery missed. In this arc: 4 hidden readers surfaced (macropore `MACROINTEGRAL`, `waterbalance.f90:calcgwl` fallback branch, vestigial `hbot` in `config_to_variables.f90`, `soilhydraulics.f90` diagnostics block). Expect 4–8 readers per arc; plan for 4–8 compiler-driven fix-up commits. The pattern is: drop dual-write → compile → fix reader → compile again → iterate.

6. **Mini-sim writeback retargeting is a recurring pattern — flag it during discovery.** Output routines that snapshot state, run a parameter-perturbation mini-sim, then restore (drainage's mini-sim was first; boundary's `swapoutput.f90:3745–3820` is the second) need their restore writes retargeted to `state%X` when the field migrates. The pattern is now established; during discovery, grep `swapoutput.f90` for any snapshot+restore block touching the subsystem's owned fields and flag it as a Cat 3 hazard. The fix is always mechanical (retarget the writeback line); the risk is forgetting to check.

7. **Validator error-code distinction — `ERR_VALIDATION_OUT_OF_RANGE` (300) vs `ERR_VALIDATION_ENUM` (301).** When adding pFUnit tests for typed-config validators, distinguish range-check errors from enum-check errors. `swcofqhc` is an integer switch (enum-style), so its out-of-value test expects `ERR_VALIDATION_ENUM` (301), not `ERR_VALIDATION_OUT_OF_RANGE` (300). Caught during B-0.2. General rule: check whether the validator uses `validate_range` or `validate_enum` and match the expected error code in the test.

8. **Init-order trap (recap from heat, confirmed again) — `soilwater_init` before `DoTillage(1)`.** `soilwater_init` was placed at `swap.f90:183`, immediately after `CalcGrid()` and before `DoTillage(1)`, matching the pattern established in the heat arc. This is necessary for forward compatibility: later arcs (crop-uptake) will add per-node arrays to `soilwater_state_t` that must be allocated before tillage runs. The general rule (from heat lesson #1) generalizes: every typed-state init runs early in the `CalcGrid` → `DoTillage(1)` → `SoilWater(1)` sequence; subsequent arcs slot their inits into the same position naturally.

---

---

## Lessons from crop-uptake migration (ADR 0036, added 2026-05-11)

Crop-uptake was migration #6 and the second coupling-surface arc of the soil-water decomposition. Six lessons generalize to future arcs.

1. **`soilwater_init` signature-bump precedent — bump once, all subsequent arcs benefit.** When `state%X_init` starts as `(state)` for scalar-only arcs (boundary) and a later arc needs per-node arrays, bump the signature to `(state, numnod, nlay)` in that arc. Do not try to retrofit per-arc with different dimension arguments — one clean signature change covers all subsequent arcs. The boundary design anticipated this explicitly ("later arcs will add per-node arrays"); crop-uptake landed it exactly as planned. Atmosphere and soil-water-core arcs inherit the stable signature without further surgery.

2. **Move init to the state-aware dispatcher, not the leaf subroutines.** When a global is initialized inside subroutines that lack state plumbing (here: `MatricFlux(1)` called from CropFixed/Wofost/Grass task=1 init paths), and threading state into those leaf subroutines is expensive, move the init to the dispatcher that already has state (`CropGrowth(1)` in this arc). The dispatcher becomes the canonical init site; the leaf paths lose their co-write. This generalizes boundary's lesson #2 (defer un-plumbed co-writers): when the dispatcher is state-plumbed but the leaves are not, the dispatcher is the natural migration point. Note: if the init also depends on values computed after state init (as `mfluxtable` depended on `SoilHydraulics(1)` outputs), allocate in `soilwater_init` but build in the dispatcher.

3. **Optional-state-arg for 5+ hidden readers in large compute subroutines.** When compile-driven discovery surfaces many readers spread across large subroutines (here: `flWrtNonox` at 5 sites across `CropFixed(3)`, `CropWofost(3)`, `CropGrass(3)/(3)/(3)` — three large compute blocks where full state threading would require touching hundreds of lines), add `optional intent(in) :: state` to the subroutine and thread state from the dispatcher. The bare-name fallback (read legacy global when state absent) handles call chains without state. Performance is irrelevant (flag reads). This is symmetric with the heat arc's `MatricFlux` optional-state pattern (playbook gotcha #3) and confirms the pattern generalizes beyond single-callee cases.

4. **Stub-errored config paths are zero-risk migration.** Fields gated by config switches that stub-error at TOML parse time (here: `swdrought=2` gates all 14 JvL fields; `swcalt=1` in the heat arc gated the analytical-method path) cannot be exercised at runtime through TOML. Their migration is byte-identical-free by construction — the JvL compute block can change freely without regression impact. The work is code-review hygiene, not regression risk. Document this explicitly in the ADR so future reviewers do not demand integration test coverage for un-exercisable paths.

5. **Naming-prefix discriminator for cohort ownership (`q*` / `i*` / `c*`).** When a subsystem has multiple temporal variants of the same field (here: `qrot` instantaneous, `iqrot` period-sum, `cqrot` cumulative, `qpotrot_day` per-day), the prefix discriminates arc ownership. Crop-uptake owns the `q*` instantaneous fields; soil-water-core owns the `i*` and `c*` accumulators and `*_day` intermediates. Apply this taxonomy during discovery Section 2 categorization to cleanly separate arc scope without ambiguity. When in doubt: who resets it (`flzerointr` / `flzerocumu` gates) determines the owner.

6. **Compile-iteration count as arc-maturity signal.** The number of compile-driven Phase 2.7 iterations tracks how much the pattern has matured: boundary needed 8 iterations (first coupling-surface arc; many surprises); crop-uptake needed 2 iterations (second arc; prior lessons applied). Each successive arc benefits from the playbook accumulating. When scoping a new arc, estimate compile-iteration count by counting how many state-arg windfalls accumulated from prior arcs — more windfalls means fewer hidden readers means fewer iterations.

---

---

## Lessons from atmosphere migration (ADR 0037, added 2026-05-11)

Atmosphere was migration #7 and the third coupling-surface arc of the soil-water decomposition. Six lessons generalize to future arcs.

1. **Pre-flight dual-write coverage check — verify non-zero writes before reader cutover.** Before cutting over any reader in Phase 2, confirm each field has at least one non-zero dual-write site (beyond zero-init and `reset()` calls):

   ```bash
   grep -rn "state%<sub>%<field>\s*=" src/ | grep -v "0\.0\|0_real64"
   ```

   Must return at least one result. If only zero-init or `reset()` writes are found, dual-write coverage is missing for non-zero accumulation paths — fix in Phase 1 first. **Caught by A-2.1's `nraidt`/`aintcdt` regression**: accumulator dual-writes in `waterbalance.f90` had not been added before reader cutover, causing the regression. A-2.2 and A-2.5 applied this check successfully.

2. **Multi-owner cumulative reset → cohort `reset()` consolidation.** When several files reset overlapping cohort fields under the same `if (flzero*)` block (atmosphere had `meteoday.f90` + `snow.f90` + `soilhydraulics.f90` all resetting cumulative cohort fields at 3 separate inline blocks), consolidate to a single `call state%X%cumu%reset()` at the canonical reset site (the first occurrence in the timestep, typically the per-day init). Drop all redundant scattered zero-writes. Confirms ADR 0033's cohort pattern value beyond the surfacewater pilot.

3. **Cohort-pattern debut in new records — introduce at type-creation, not as retrofit.** When a state record is introduced from scratch and the subsystem has both cumulative and intermediate fields sharing reset gates, design the cohorts into the type from the start (`atmosphere_state_t` is the worked example). surfacewater retrofitted cohorts onto an existing flat type (ADR 0033 Phases A+B); atmosphere introduced them at type-creation, which is cleaner. Future new state records with `flzero*`-gated fields should follow the atmosphere precedent.

4. **Co-writer accumulators of cohort fields surface dual-write gaps during Phase 2.** When Phase 2 reader cutover causes regressions, the cause is often that a routine OUTSIDE the home tree ACCUMULATES cohort fields (here: `waterbalance.f90:UpdateWaterBalance` accumulates `igrai`, `inrai`, `ipeva`, `iptra`, `cpeva`, `cptra`, `cevap`, `caintc`, `cgrai`, `cnrai` — 10 atmosphere intermediate+cumulative fields). The home-tree dual-writes cover the compute-write sites; co-writer accumulator sites need separate dual-write additions. During Phase 1 design, explicitly inventory accumulate-write sites in `waterbalance.f90` and other external co-writers — do not rely on grep of the home tree alone.

5. **Architectural state-arg additions during Phase 2.6 are a natural cost, not scope creep.** Compile-driven discovery surfaces routines that need `state` for the first time when hidden global references are exposed. In this arc: `DoTillage`, `CNmethod`, `checkmassbal`, and `Consolidate_Bdens` all gained `state` args during Phase 2.6. These routines had no state plumbing because no prior migration needed them; atmosphere's wider coupling surface is what surfaces them. Budget 4–8 such additions per arc proportional to the number of external reader files. They are structural improvements, not rework.

6. **Field-count progression as scope indicator for arc decomposition.** Heat: 11 fields. Boundary: 12. Crop-uptake: 23. Atmosphere: 40. Larger arcs benefit disproportionately from coupling-surface decomposition — a single soil-water mega-arc would have been unmanageable. Atmosphere's cohort layout and dual 7-file home tree would not have been tractable under the boundary/heat flat precedent. When a subsystem's field-count exceeds ~25 AND the fields span multiple reset cadences, coupling-surface decomposition is the right default, even if it produces 4+ sequential arcs.

---

---

---

## Lessons from soil-water-core migration (ADR 0038, added 2026-05-11)

Soil-water core was migration #8 and the FINAL coupling-surface arc of the soil-water decomposition. Six lessons generalize beyond this arc.

1. **Compile-driven retirement (Strategy B) — the master method for dropping legacy dual-writes.** When the time comes to drop the legacy half of Phase 1 dual-writes and retire the globals, do NOT attempt to manually identify which file to touch first. The which-file-first ordering problem becomes intractable for arcs with 40+ globals spread across 20+ files:

   - Comment out all owned globals in `variables.f90` first (with provenance markers).
   - The compiler enumerates every remaining legacy reference as a precise compile error.
   - Each error is an exact instruction: drop the legacy write line, OR migrate a missed read to state, OR drop a use-clause import.
   - Iterate: fix the batch of errors, recompile, fix the next batch. Typically 20–30 passes for a large arc.
   - One commit resolves everything.

   This eliminates all guesswork and avoids the "partially-dropped dual-write" state that caused NR-divergence hangs in two prior attempts. Cost: one longer iterative-compile commit. Benefit: structurally unambiguous correctness — the compiler is the exhaustive inventory.

2. **Hidden NR-coupled readers — the canary pattern.** When dropping legacy writes causes an NR-divergence hang (hupselbrook hangs at 99.9% CPU with no output, no error), the bug is: a utility routine called from INSIDE the Newton-Raphson iteration reads a dropped legacy global and gets stale zero-init values. In soil-water, `soilhydraulicsutils.f90` (`hconduc`/`dhconduc`/`watcon`/`moiscap`) and `WC_K_models_04_11.f90` both read `cofgen` and `fluseksatexm` from inside `headcalc`'s NR loop. Once the legacy writes stopped, those reads returned zeros → Richards diverged → hang.

   **General rule:** before dropping legacy writes for any array that feeds the Richards solver, grep for `use variables, only: <array>` in all utility modules called from `headcalc`, `fluxes`, or any tight physics loop. Treat hupselbrook as the canary — if it hangs, the root cause is always a stale NR-coupled read, not a logic bug.

3. **Module-level pointer pattern for utility modules with many callers.** When a utility module (`soilhydraulicsutils.f90`, `WC_K_models_04_11.f90`) provides functions called from 80+ sites and reads a soil-physical-property array (`cofgen`, `fluseksatexm`), threading `state` into all callers is prohibitively expensive. Cleaner: add module-level pointers in the utility module; have the home file (`soilhydraulics.f90`) call a `bind_state_targets` helper once at init time that points those pointers at the state arrays. The utility functions read through the pointers without any signature change; all 81 call sites stay untouched.

   This generalizes: any utility module with a flat `use variables, only: <field>` and 30+ call sites is a pointer-pattern candidate. Threading state is correct but expensive; pointer binding is pragmatic for deeply-shared utilities.

4. **Transient buffer pattern for pre-init config-adapter state writes.** `config_to_variables` runs before any `*_init` call in the startup sequence, so state allocatables are not yet allocated. Use module-level transient buffers in the adapter: capture the values during config-to-variables, then drain them into state after `soilwater_init` completes. This is the same pattern as atmosphere A-2.6's `ssnow`/`ldwet` handling, now confirmed general:

   **General rule:** any config-adapter or init-time writer that runs before `<subsystem>_init` needs either (a) an `allocated()` guard at the write site (heat arc lesson #1) or (b) a transient buffer if the field is a per-node allocatable that cannot be written without allocation. Allocatable scalars can use the `allocated()` guard; per-node arrays need the buffer.

5. **Twin-cohort two-gate-reset — ADR 0033 mild extension.** A single cohort sub-record can carry fields with TWO different reset gates without splitting into two separate types. `soilwater_intermediate_t` demonstrates: `reset()` under `flzerointr` zeroes all 22 fields; `reset_per_day()` under `flDayStart` zeroes only the 8 per-day fields nested inside the same cohort. This is cleaner than adding a separate `soilwater_per_day_t` — it keeps `soilwater_state_t`'s shape parallel to `atmosphere_state_t` (two cohort components, not three) while expressing both reset cadences in the type system.

   **General rule:** when discovery identifies per-day fields that naturally belong to the intermediate reset family (they accumulate within a day and are rolled up into the intermediate period), fold them into `<subsystem>_intermediate_t` with a second `reset_per_day()` method. Only split into a separate type if the per-day fields have a fundamentally different ownership story.

6. **4-arc coupling-surface decomposition validated end-to-end.** The original soil-water mega-discovery (~100 globals, 897 lines) was rejected as a single arc. Four coupling-surface arcs shipped successively: boundary 12 + crop-uptake 23 + atmosphere 40 + core 74 = 149 fields across 4 commits, each byte-identical from Phase 1 through retirement.

   Key measures of success: (a) each arc scoped below 80 owned globals, keeping blast radius manageable; (b) each arc's Phase 1 windfall grew because prior arcs plumbed state into more entry points; (c) no rework was required between arcs — the coupling-surface interfaces were stable; (d) hupselbrook stayed green at 5/5 throughout all four arcs.

   **Generalizes as:** when a subsystem estimates >50 owned globals, look for natural coupling surfaces to decompose — the boundary between callers (top/bottom, crop-sink, atmosphere, core-Richards) is the right cut, not an arbitrary file-count split. The final core arc benefits most from prior arcs' state-arg windfalls.

---

---

## Lessons from tillage migration (ADR 0039, added 2026-05-12)

Tillage was migration #9 — the smallest arc to date (6 tasks, 2 compile passes, T-4 reader
cutover a no-op). Three lessons generalize to future arcs.

1. **Config-constant vs runtime-state distinction — split the inventory before scoping.**
   When retiring globals from a subsystem with mixed legacy storage, classify each owned
   global by mutability, not just naming convention:
   - Fields that are **immutable after init** (config-table copies, switch parameters,
     event-table arrays) are config-constants. If a typed config home already exists, leave
     them as legacy globals until a config-consolidation arc retires them. Migrating them to
     a runtime state record would duplicate typed config.
   - Fields that are **per-event computed, mutated, or accumulated** during the simulation
     are genuine runtime state and belong in `<subsystem>_state_t`.
   Apply this split during discovery Section 2 categorization — it sets the scope for the
   entire arc. Tillage debuted this pattern: 18 of 31 `till_*` globals were config-constants
   (kept legacy, Groups A+B); only 13 runtime-state fields migrated (Groups C+D+E).

2. **Cross-subsystem side-effect globals stay legacy until a dedicated ownership arc.**
   When a subsystem writes globals that are read by multiple OTHER subsystems for coupling
   (here: tillage writes `Bdens` and `ParamVG`; `soilhydraulics.f90`, `solute.f90`, and
   `oxygenstress.f90` read them), those fields cannot be absorbed into the writer's state
   record without also migrating all downstream readers. Cross-subsystem ownership
   clarification is its own arc. Leave such fields as legacy globals with a comment noting
   the deferred ownership work. Trying to absorb them into the writer's record forces
   unrelated downstream plumbing and expands blast radius beyond the arc's focus.

3. **Self-contained subsystem hint — flag zero-external-readers during discovery.**
   If discovery shows zero external compute readers of the subsystem's owned globals (outside
   the home file and co-writer), the reader cutover task (T-4 in the standard plan) may be
   a no-op. Flag this during discovery so the plan can drop or stub that task and reduce
   overall arc scope. Tillage was the first arc to surface this: the "8 external sites"
   initially found during discovery turned out to be Group A+B config-constants in
   `config_to_variables.f90` (out of scope), not runtime-state readers. Strategy B (T-5)
   confirmed: 2 compile passes, zero external reader fixes needed.

---

## Reference ADRs

- ADR 0030 — Surface-water state-type migration (pilot)
- ADR 0031 — Drainage state-type migration
- ADR 0032 — Solute state-type migration
- ADR 0033 — Cumulative reset cohorts
- ADR 0034 — Heat subsystem state-type migration
- ADR 0035 — Boundary subsystem state-type migration (first coupling-surface arc)
- ADR 0036 — Crop water uptake state-type migration (second coupling-surface arc)
- ADR 0037 — Atmosphere subsystem state-type migration (third coupling-surface arc)
- ADR 0038 — Soil-water core state-type migration (FINAL coupling-surface arc)
- ADR 0039 — Tillage state-type migration (smallest arc; config-constant vs runtime-state distinction)
