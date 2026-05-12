# Tillage State-Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Carve 13 runtime-state `till_*` globals into `state%tillage` (NEW state record); fix `set_iTill` tautological-condition bug; preserve byte-identical check-full at every commit.

**Architecture:** Flat `tillage_state_t` with 7 per-layer arrays (Group C) + 3 per-step scalars (Group D) + 3 init-once geometry/cursor fields (Group E). New `tillage_init` wired at `swap.f90` between `atmosphere_init` (line 214) and `DoTillage(1)` (line 227). Pre-init pattern: `config_to_variables` writes 18 config-constant `till_*` legacy globals (unchanged); `DoTillage(1)` seeds `state%tillage` Group E fields after `tillage_init` allocates Group C arrays. Strategy B compile-driven retirement for the 13 migrated globals.

**Tech Stack:** Fortran 2008+, ASSOCIATE, pixi+meson, pFUnit, check-full regression.

**Spec:** `docs/superpowers/specs/2026-05-12-state-migration-tillage-design.md`
**Discovery:** `docs/superpowers/specs/2026-05-12-state-migration-tillage-discovery.md`
**Branch:** all commits go on `refactor/tillage-state`.

---

## Lessons-learned applied

- **Strategy B compile-driven retirement (ADR 0038 lesson #1).** Comment out all 13 migrated globals in `variables.f90`; compiler enumerates every remaining reference. No guesswork, no partial-drop state. One iterative-compile commit.
- **Pre-flight dual-write coverage check (ADR 0037 atmosphere lesson #1).** Before cutting over any reader, confirm `grep "state%tillage%X\s*=" src/ | grep -v "0\.0\|0_real64"` returns at least one hit per field.
- **H-4 fix safety verification (D9).** Run `grep -rn "flTillage\s*=\s*.true." tests/` before T-3 commit. If all regression cases have `flTillage=false`, the H-4 fix is byte-identical by construction.
- **Grep field shapes during pre-flight.** Before retiring a global in T-5, grep both `use variables.*\btill_X\b` AND raw-symbol forms. Both must return zero non-state hits.
- **ASSOCIATE prefix `tl_`.** Shadow-safe aliasing inside `tillage.f90` compute bodies. Follows `ht_` (heat), `bw_` (boundary), `cw_` (crop-water), `at_` (atmosphere) convention.

---

## File structure

**Created:**
- `src/state/tillage_state.f90` — `tillage_state_t` + `tillage_init` (T-1, T-2)
- `tests/unit/state/test_tillage_state.pf` — pFUnit tests (T-1)

**Modified:**
- `src/state/swap_state.f90` — add `type(tillage_state_t) :: tillage` (T-1)
- `src/core/swap.f90` — `tillage_init` call between `atmosphere_init` and `DoTillage(1)` (T-2)
- `src/crop/tillage.f90` — dual-write + H-4 fix + reader cutover + drop dual-writes (T-3, T-4)
- `src/io/toml/config_to_variables.f90` — Group C/D/E reads replaced by `state%tillage%X` (T-4; only if it reads Group C/D/E fields — discovery confirms it does not; verify)
- `src/core/variables.f90` — comment out 13 retired globals (T-5)
- `tests/unit/testSuites.inc`, `tests/unit/meson.build` — register new test suite (T-1)
- meson source list — add `src/state/tillage_state.f90` (T-1)
- `docs/adr/0039-state-migration-tillage.md` — final ADR (T-6)

---

## Task T-1: Create `tillage_state_t` + add to `swap_state_t`

**Design decisions:** D2, D3, D4
**Files:**
- Create: `src/state/tillage_state.f90`
- Create: `tests/unit/state/test_tillage_state.pf`
- Modify: `src/state/swap_state.f90`
- Modify: `tests/unit/testSuites.inc`, `tests/unit/meson.build`, meson source list

13-field flat type. No `tillage_init` wiring yet (that is T-2). The 7 Group C fields are `allocatable` with no default initializer (allocatables cannot have scalar initializers). The 6 non-allocatable fields carry zero/false defaults.

- [ ] **Write failing pFUnit tests** — in `test_tillage_state.pf`, add `@test` stubs:
  - `test_tillage_state_scalar_defaults` — verify `sumDWC`, `sumAvail1`, `sumAvail2`, `MaxNumSoilHo`, `MaxNumSoilCP`, `iTill` all zero on fresh instance.
  - `test_tillage_state_allocatables_not_allocated` — verify all 7 per-layer arrays are `not allocated` before `tillage_init`.
  - `test_tillage_state_independent_instances` — two instances mutate independently.
  - Register suite. Confirm FAIL with `pixi run -e test test-pfunit`.
- [ ] **Create `tillage_state_mod`** in `src/state/tillage_state.f90`:
  ```fortran
  module tillage_state_mod
     use, intrinsic :: iso_fortran_env, only: real64
     implicit none
     private
     public :: tillage_state_t

     type :: tillage_state_t
        ! Group C — per-layer working state (allocated by tillage_init to NumLay)
        real(real64), allocatable :: Rho_tillage(:)  !! post-event target bulk density per layer
        real(real64), allocatable :: Rho_cons(:)     !! consolidated density target per layer
        real(real64), allocatable :: Rho_last(:)     !! density at start of current step per layer
        real(real64), allocatable :: K_R_cons(:)     !! consolidation rate per layer
        real(real64), allocatable :: Rho_match(:)    !! matching-point density per layer
        real(real64), allocatable :: N_match(:)      !! matching-point n-value per layer
        real(real64), allocatable :: Slope_match(:)  !! slope at matching point per layer
        ! Group D — per-step working scalars (written by Adapt_WC_H each call)
        real(real64) :: sumDWC    = 0.0_real64  !! sum of water content changes
        real(real64) :: sumAvail1 = 0.0_real64  !! available pore space (wetting)
        real(real64) :: sumAvail2 = 0.0_real64  !! available water (draining)
        ! Group E — init-once geometry/cursor (set by det_MNSH / set_iTill at DoTillage(1))
        integer :: MaxNumSoilHo = 0  !! horizon count up to Max_Z_tillage
        integer :: MaxNumSoilCP = 0  !! node count up to Max_Z_tillage
        integer :: iTill        = 0  !! current event-table cursor index
     end type tillage_state_t

  end module tillage_state_mod
  ```
- [ ] **Add `state%tillage` to `swap_state_t`** — add `use tillage_state_mod, only: tillage_state_t` and `type(tillage_state_t) :: tillage` field. Update meson source list.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
feat(state): SS-TIL T-1 — tillage_state_t flat type + swap_state aggregation

13-field flat tillage_state_t: 7 per-layer allocatables (Group C),
3 per-step scalars (Group D), 3 init-once geometry/cursor integers
(Group E). No cohorts (no flzero* gating anywhere in tillage).
Aggregated as state%tillage in swap_state_t.

ADR: docs/adr/0039-state-migration-tillage.md (stub)
Spec: docs/superpowers/specs/2026-05-12-state-migration-tillage-design.md
EOF
)"
```

---

## Task T-2: Add `tillage_init` + wire in `swap.f90`

**Design decisions:** D7, D8
**Files:**
- Modify: `src/state/tillage_state.f90` — add `tillage_init` public subroutine
- Modify: `src/core/swap.f90` — call `tillage_init` between `atmosphere_init` (line 214) and `DoTillage(1)` (line 227)

`tillage_init(state, numlay)` allocates the 7 Group C per-layer arrays to `numlay` and zeroes the Group D scalars. Group E fields (`MaxNumSoilHo`, `MaxNumSoilCP`, `iTill`) are left at their zero defaults; `DoTillage(1)` sets them via `det_MNSH` and `set_iTill`.

- [ ] **Inspect `swap.f90` lines 210–230** — confirm `atmosphere_init` at line 214 and `DoTillage(1)` at line 227; identify exact insertion point (after `atmosphere_init`, inside the `if (flTillage)` guard or outside).
- [ ] **Implement `tillage_init(state, numlay)`** — public subroutine in `tillage_state_mod`. Allocates `state%Rho_tillage(numlay)`, etc. for all 7 Group C arrays. Sets Group D scalars to 0.0. Does NOT set Group E fields (deferred to DoTillage(1)).
  - Note: place the call outside the `if (flTillage)` guard — allocation is unconditional; DoTillage(1) only runs if flTillage. This matches `soilwater_init` precedent.
- [ ] **Add pFUnit allocation lifecycle test** — `test_tillage_init_allocates_numlay`: call `tillage_init(t, 5)`, verify all 7 arrays have `size == 5` and are allocated.
- [ ] **Wire in `swap.f90`** — add `use tillage_state_mod, only: tillage_init` (or via `swap_state_mod` if exported there); insert `call tillage_init(state%tillage, numlay)` after `atmosphere_init` call.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
feat(state): SS-TIL T-2 — tillage_init allocates 7 per-layer arrays; wired in swap.f90

tillage_init(state%tillage, numlay) allocates Group C per-layer arrays
to NumLay. Inserted at swap.f90 after atmosphere_init (line 214) and
before DoTillage(1) (line 227). Group E geometry/cursor fields set later
by det_MNSH/set_iTill inside DoTillage(1).

check-full 5/5 byte-identical.
Spec: docs/superpowers/specs/2026-05-12-state-migration-tillage-design.md
EOF
)"
```

---

## Task T-3: Dual-write all 13 fields in `tillage.f90` + fix H-4

**Design decisions:** D8, D9, D10
**Files:**
- Modify: `src/crop/tillage.f90`

Dual-write phase: add `state%tillage%X = X` mirrors for all 13 runtime fields at every write site inside `tillage.f90`. Legacy global writes stay alive. Use ASSOCIATE with `tl_` prefix in compute bodies. Also fix the H-4 `set_iTill` tautological-condition bug.

- [ ] **Pre-flight: H-4 safety check** — run `grep -rn "swtill\|flTillage" tests/` and `grep -rn "flTillage" /home/zawadzkim/Code/swap/tests/swap-cases/ 2>/dev/null || echo "no swap-cases"`. If all regression cases have `flTillage=false` (likely), the H-4 fix is byte-identical by construction. Document finding in commit message.
- [ ] **List all write sites for 13 fields** — `grep -n "Rho_tillage\|Rho_cons\|Rho_last\|K_R_cons\|Rho_match\|N_match\|Slope_match\|sumDWC\|sumAvail1\|sumAvail2\|MaxNumSoilHo\|MaxNumSoilCP\|iTill" src/crop/tillage.f90 | grep -v "state%\|!"`. Cross-reference discovery Section 1 subroutine table (DoTillage(1,2), Adapt_WC_H, Consolidate_Bdens, Change_Tillage_Info, set_iTill, det_MNSH).
- [ ] **Add ASSOCIATE `tl_` block** at the top of each relevant subroutine's compute body:
  ```fortran
  associate( &
     tl_Rho_tillage => state%tillage%Rho_tillage, &
     tl_Rho_cons    => state%tillage%Rho_cons,    &
     ! ... all 13 fields
  )
  ```
- [ ] **Add dual-writes** — after each legacy global write (e.g., `Rho_tillage(i) = ...`), add `tl_Rho_tillage(i) = Rho_tillage(i)`. For Group E: `state%tillage%MaxNumSoilHo = MaxNumSoilHo`, etc. after `det_MNSH` and `set_iTill` compute their values.
- [ ] **Fix H-4 in `set_iTill`** — change the tautological condition:
  ```fortran
  ! BEFORE (bug — same index on both sides; branch always false):
  if (t1900 >= Date_tillage(i-1) .and. t1900 < Date_tillage(i-1)) then
  ! AFTER (fix — second index advances to i):
  if (t1900 >= Date_tillage(i-1) .and. t1900 < Date_tillage(i)) then
  ```
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`. If H-4 fix causes regression (only possible if any check-full case uses multi-event tillage), isolate the fix to a separate commit after investigating.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-TIL T-3 — dual-write 13 tillage fields; fix set_iTill H-4 bug

All 13 till_* runtime-state fields (Groups C+D+E) dual-written in
tillage.f90 (DoTillage(1,2), Adapt_WC_H, Consolidate_Bdens,
Change_Tillage_Info, set_iTill, det_MNSH). ASSOCIATE tl_ prefix.
Legacy globals unchanged — dual-write is additive.

H-4 fix: set_iTill tautological condition corrected (Date_tillage(i-1)
< Date_tillage(i-1) → < Date_tillage(i)). Loop now advances iTill past
event 1. Byte-identical on check-full (flTillage=false in all 5 cases).

Spec: docs/superpowers/specs/2026-05-12-state-migration-tillage-design.md
EOF
)"
```

---

## Task T-4: Reader cutover — external read sites

**Design decisions:** D1, D8
**Files:**
- Modify: `src/crop/tillage.f90` — internal reads of Group C/D/E fields switch to `state%tillage%X`
- Verify: `src/io/toml/config_to_variables.f90` — reads only Group A+B fields; no Group C/D/E reads (discovery confirmed); no change expected

Discovery found all external read/write sites of Group C/D/E are inside `tillage.f90` itself (the subsystem is self-contained). `config_to_variables.f90` only writes Groups A+B config-constants. `swap.f90` passes `state` as an arg; no direct `till_*` reads.

- [ ] **Pre-flight dual-write coverage check** — for each of the 13 fields, run:
  ```bash
  grep -rn "state%tillage%X\s*=" src/ | grep -v "0\.0\|0_real64"
  ```
  Must return at least one non-zero write per field. If any field only has zero-init writes, go back to T-3 and add the missing dual-write.
- [ ] **List internal reads in `tillage.f90`** — `grep -n "Rho_tillage\|Rho_cons\|Rho_last\|K_R_cons\|Rho_match\|N_match\|Slope_match\|sumDWC\|sumAvail1\|sumAvail2\|MaxNumSoilHo\|MaxNumSoilCP\b\|iTill\b" src/crop/tillage.f90 | grep -v "state%\|!.*retired"`. Categorize each as write vs. read.
- [ ] **Replace reads with `state%tillage%X`** — update ASSOCIATE blocks if needed so reads pull from `tl_X => state%tillage%X`. All compute reads in DoTillage(2), Adapt_WC_H, Consolidate_Bdens, Change_Tillage_Info that use Group C/D/E values must now read from `state%tillage`.
- [ ] **Check external files** — `grep -rn "till_Rho_tillage\|till_sumDWC\|till_MaxNumSoilHo\|till_iTill" src/ --include="*.f90" | grep -v "variables.f90\|tillage.f90\|state%"`. Expected: zero hits outside home file and variables.f90. If any appear, migrate them.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-TIL T-4 — reader cutover; tillage.f90 reads state%tillage for 13 fields

All internal reads of Group C/D/E till_* fields in tillage.f90 now
read from state%tillage%X via ASSOCIATE tl_ aliases. Dual-writes still
active (legacy globals still written). State authoritative for reads.

No external readers found outside tillage.f90 (subsystem self-contained
per discovery Section 3).

Spec: docs/superpowers/specs/2026-05-12-state-migration-tillage-design.md
EOF
)"
```

---

## Task T-5: Strategy B compile-driven retirement

**Design decisions:** D10
**Files:**
- Modify: `src/core/variables.f90` — comment out 13 `till_*` runtime-state globals
- Modify: `src/crop/tillage.f90` — drop dual-write legacy-global sides; state-only writes remain

Apply Strategy B from ADR 0038 (soil-water core playbook lesson #1): comment out all 13 globals first, let the compiler enumerate every remaining reference, fix systematically.

- [ ] **Pre-flight — enumerate legacy write sites** — `grep -n "Rho_tillage\|Rho_cons\|Rho_last\|K_R_cons\|Rho_match\|N_match\|Slope_match\|till_sumDWC\|till_sumAvail1\|till_sumAvail2\|till_MaxNumSoilHo\|till_MaxNumSoilCP\|till_iTill" src/crop/tillage.f90 | grep -v "state%"`. These are the dual-write legacy sides to drop.
- [ ] **Comment out 13 globals in `variables.f90`** — with provenance markers:
  ```fortran
  ! real(8), allocatable :: till_Rho_tillage(:)  ! [SS-TIL] retired 2026-05-12 — moved to state%tillage%Rho_tillage (ADR 0039)
  ```
  Do all 13 in one edit pass (7 allocatables in Group C, 3 scalars in Group D, 3 integers in Group E).
- [ ] **First compile** — `pixi run build 2>&1 | grep "Error\|error:" | head -40`. Expect errors in `tillage.f90` for the 13 legacy write sides. May also surface 3–5 hidden readers in other files not enumerated by discovery.
- [ ] **Iterate: fix each batch of compile errors** — for each error:
  - In `tillage.f90`: remove the legacy-global half of each dual-write (keep `tl_X = ...` or `state%tillage%X = ...`). Drop `till_*` from `use variables, only:` clauses where fully retired.
  - In other files (if any surface): migrate remaining reads to `state%tillage%X`. If a file lacks `state` in scope, thread it or use the optional-state-arg pattern (playbook gotcha #3).
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full` (integration gate — must be byte-identical).
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
refactor(state): SS-TIL T-5 — retire 13 till_* runtime-state globals (Strategy B)

Comment out 13 till_* runtime-state fields in variables.f90 with
[SS-TIL] provenance markers. Compile-driven discovery resolved N
hidden readers (see commit body). state%tillage is now authoritative
for all Group C/D/E fields.

Groups A+B (18 config-constant till_* globals) remain in variables.f90
(future config-consolidation arc). Bdens/ParamVG untouched (D6).

pFUnit green. check-full 5/5 byte-identical.
ADR: docs/adr/0039-state-migration-tillage.md
Spec: docs/superpowers/specs/2026-05-12-state-migration-tillage-design.md
EOF
)"
```

---

## Task T-6: ADR 0039 + playbook update + merge prep

**Design decisions:** all
**Files:**
- Create: `docs/adr/0039-state-migration-tillage.md`
- Modify: `docs/superpowers/specs/state-migration-playbook.md` — append 1-2 new lessons

Author ADR 0039 and append any new playbook lessons surfaced by this arc.

- [ ] **Author `docs/adr/0039-state-migration-tillage.md`** — follow ADR 0035 structure (status, date, migration #, context, decisions D1–D10, deferred items, consequences, known residuals, references). Key points to capture:
  - Migration #9; flat `tillage_state_t`; smallest arc to date.
  - Groups A+B deferral rationale (already in typed `soil_tillage_t`).
  - Bdens/ParamVG stay-legacy rationale (cross-subsystem ownership).
  - H-4 `set_iTill` fix.
  - Strategy B compile-driven retirement; N hidden readers found.
  - Phase 0 zero — ADR 0021 windfall.
  - check-full 5/5 byte-identical; pFUnit N passing, 0 failures.
- [ ] **Append playbook lessons** — consider adding to `state-migration-playbook.md` under "Lessons from tillage migration (ADR 0039)":
  - **Cross-subsystem side-effect globals:** when a subsystem writes globals consumed by multiple other subsystems for side-effect coupling (here: Bdens/ParamVG written by tillage, read by soilhydraulics/solute/oxygenstress), defer cross-subsystem ownership clarification; stay legacy. Moving such globals into the writing subsystem's state would silently sever the reading subsystems.
  - **Config-constant vs runtime-state distinction:** when retiring globals, classify them first. Constants already in a typed config record (Groups A+B in `soil_tillage_t`) should not be duplicated in `state%X` — they belong to a config-consolidation arc. Only retire the genuine runtime state (Groups C+D+E) in the state-migration arc.
- [ ] **Verify PASS** — `pixi run -e test test-pfunit && pixi run check-full`.
- [ ] **Commit:**
```bash
git commit -m "$(cat <<'EOF'
docs(adr): SS-TIL T-6 — ADR 0039 tillage state-migration + playbook update

ADR 0039: migration #9 complete. Flat tillage_state_t (13 runtime fields);
Groups A+B (18 config-constants) deferred to config-consolidation arc;
Bdens/ParamVG stay legacy (cross-subsystem coupling). H-4 set_iTill
tautological-condition fix. Strategy B compile-driven retirement.

Playbook: two new lessons — cross-subsystem side-effect globals; config-
constant vs runtime-state distinction.

pFUnit N passing. check-full 5/5 byte-identical. Smallest arc to date.
Spec: docs/superpowers/specs/2026-05-12-state-migration-tillage-design.md
EOF
)"
```

---

## Self-Review Notes

- **H-4 fix byte-identical assumption:** depends on all 5 check-full cases having `flTillage=false`. If any case has multi-event tillage tables and `flTillage=true`, the corrected `set_iTill` would advance `iTill` further and produce different output. The implementer MUST verify with `grep -rn "swtill\s*=\s*1\|flTillage" tests/` before committing T-3. If any case uses tillage, run check-full before AND after the H-4 fix separately.

- **T-3 ASSOCIATE scope:** the 7 Group C arrays are allocatables. The ASSOCIATE construct works correctly with allocatable targets — ensure `state%tillage` is argument-associated through `state` (which it is, since `DoTillage` takes `state intent(inout)`). No `allocated()` guard needed at dual-write sites because `tillage_init` runs unconditionally before `DoTillage(1)`.

- **T-5 hidden reader estimate:** discovery found zero external compute readers of Group C/D/E fields. The 8 external "read/write sites" in discovery Section 3 are: `swap.f90` (3 call-site lines — not direct `till_*` reads), `config_to_variables.f90` (Groups A+B only, not in scope), `read_soil_tillage_toml.f90` (typed config record, not `till_*` globals). Net hidden readers outside `tillage.f90` proper may be as low as 0. Expect compile-driven discovery to find 0–3 surprises, mostly cleanup of `use variables, only: till_*` import lines.

- **`till_iTill` mutable cursor (H-4):** this field (Group E) is mutated at runtime (`iTill = iTill + 1` at tillage.f90:149). The dual-write at T-3 must include this increment site: after `iTill = iTill + 1`, add `state%tillage%iTill = iTill`. The reader cutover at T-4 must also ensure that reads of `iTill` in the event-dispatch logic pull from `state%tillage%iTill` not the legacy global.

- **`tillage_init` placement — outside `if (flTillage)` guard:** wiring the call unconditionally (not inside the flTillage guard) matches `soilwater_init` precedent and avoids a conditional-allocation footgun. The 7 per-layer arrays will be allocated even when `flTillage=false`; this is negligible memory overhead and simplifies downstream `allocated()` checks.

- **Task split rationale:** T-3 and T-4 are kept separate to isolate dual-write (additive, byte-identical by construction) from reader cutover (removes legacy-global reads, creates first dependency on `state%tillage` for correctness). This follows the established boundary/atmosphere pattern. T-5 is isolated because Strategy B requires the compiler to enumerate — mixing it with T-4 would obscure which compiler error came from which operation.
