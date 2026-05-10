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
| **Cumulative** | Reset on `flzerocumu` (cumulative-period rollover, typically yearly/output-period) | `<subsystem>_cumulative_t` | `cqdra`, `sqdra`, `samini`, `cqdrain(:)` |

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

## Reference ADRs

- ADR 0030 — Surface-water state-type migration (pilot)
- ADR 0031 — Drainage state-type migration
- ADR 0032 — Solute state-type migration
- ADR 0033 — Cumulative reset cohorts
