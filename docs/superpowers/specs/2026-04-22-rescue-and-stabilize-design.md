---
title: "Rescue & Stabilize: SWAP Refactor Recovery Plan"
author: Mateusz Zawadzki
date: 2026-04-22
status: draft
---

# Rescue & Stabilize

A workflow spec for bringing the SWAP modernization branch back to a clean, tested, documented state from which the compartment-based state refactor, Python bindings, multicore, and GPU work can proceed in separate follow-on specs.

## Background

The SWAP modernization started by lifting the legacy Fortran 77–era code (`swap_org/source_swp_4.2.0/`) into modern Fortran with explicit state management, TOML configuration, and pFUnit-based testing. The design direction — `config_t` / `initial_t` / `state_t` separation, ASSOCIATE-based local aliasing, aggregator pattern in `src/core/swap_state.f90` — is correct. The problem is drift: new work started before prior phases were closed, several modules have duplicate or orphaned WIP files, I/O has dual legacy/TOML paths coexisting, and the current tree does not compile cleanly.

The last commit at which the model builds and regression tests pass is `e256bc0` ("refactor: cleaning up tests and docs"). The drift between that commit and the current work branch is not worth salvaging commit-by-commit: it is small in volume and did not fix any known issue. The clean-slate path — reset to `e256bc0` and re-derive the forward work under discipline — is both simpler and safer than cherry-picking.

## Scope

This spec covers the work needed to turn the green baseline at `e256bc0` into a state where every pre-condition in the (soon-to-be-rewritten) compartment-state design is genuinely satisfied:

- Clean repository layout, one obvious build path, no root-level clutter.
- Canonical documentation covering architecture, state management, configuration schema, build/test, and contribution, rendered by FORD.
- A tested codebase: every `*_state_t` has lifecycle tests, every TOML reader has fixture tests, every pure physics routine has a unit test or characterization test.
- Every module has `config_t` / `initial_t` / `state_t` separation.
- The main program uses `init(state, config, initial)` everywhere.
- Every process function has `config` as `intent(in)`.
- Full regression suite green at current tolerances.

## Non-goals

- No new physics, features, or performance work. The macropore regression time (344s) and the pre-existing MOWDM deviation in `oxygenstress` are **accepted as-is** for the duration of this spec.
- No cherry-picking from the drifted work; it is archived read-only.
- No compiler migration. `gfortran` stays the only supported compiler until this spec is done.
- No compartment-based state refactor (follow-on spec).
- No multicore, Python bindings, or GPU work (follow-on specs).
- No Intel/ifx re-enablement (follow-on spec).
- No changes to physics behavior: all regression cases must match current fixtures (within current tolerances) at every checkpoint.

## Test protocol

Two tiers, used deliberately throughout this spec and onward as a permanent project convention:

| Tier | Cases | Budget | Use |
|---|---|---|---|
| **Fast** | hupselbrook, surfacewater, salinitystress, grassgrowth | ~67s | Every local iteration, every commit before phase tag |
| **Full** | Fast + oxygenstress (~178s) + macropore (~344s) | ~590s | Before every phase tag, and before any merge to `main` |

Both tiers are exposed as pixi tasks:

- `pixi run check-fast` — compiles, runs pFUnit, runs fast regression set.
- `pixi run check-full` — same plus the two slow cases.

## Git workflow for this spec

- **Local-only.** Nothing is pushed to `origin/main` for the duration of this spec.
- **Two branches kept alive:** `main` and `development`. `development` is where rescue-phase work lands directly (no per-change branches for Phase 0–3). `main` only moves when a phase is complete and full-test is green.
- **Starting in Phase 4**, per-change feature branches branch off `development`, merge back when green.
- **Discarded branches.** `swaplib` and `swaplib-simple` are deleted after confirming nothing unique lives on them. The current drifted work branch is renamed to `archive/wip-drifted` and kept locally for reference only.
- **Legacy SWAP preservation.** `swap_org/` (the reference SWAP 4.2.0 tree currently untracked on the machine) is captured as an orphan branch `legacy/swap-4.2.0` so the reference implementation is inside this repository rather than a loose directory.
- **Tagging.** Each phase exit: `rescue/phase-N-<name>`, plus `rescue/complete` at the end.

## Phases

### Phase 0 — Baseline reset

**Goal:** Unambiguous git topology anchored on the green commit, with the legacy SWAP preserved in-repo.

Steps:

1. Inventory all local and remote branches; confirm that `swaplib` and `swaplib-simple` have no unique unmerged work (diff each against `main` and `e256bc0`). If anything unique turns up, branch it to `archive/<name>` first.
2. Create `archive/wip-drifted` from the current drifted HEAD (whatever non-`main` branch the user was last working on).
3. Create `archive/main-pre-rescue` from current `main` tip (`7587ca3`), purely as a safety net.
4. Reset local `main` to `e256bc0`. **Do not push.** `origin/main` is left untouched for the duration of this spec.
5. Create `development` from `main`. All Phase 0–3 commits go directly to `development`.
6. Create `legacy/swap-4.2.0` as an orphan branch capturing a curated snapshot of the reference implementation.

   **License status (verified 2026-04-22):** SWAP 4.2.0 is distributed under GPL v2; TTUTIL 4.27 under LGPL v2.1. Redistribution is explicitly permitted, so the orphan branch may be pushed to any remote. The GPL v2 `LICENSE` file from `swap_org/source_swp_4.2.0/` is preserved at the root of the orphan branch as required by GPL §1.

   **Curation list.** The orphan branch contains only what is useful for reference and reproducibility:

   | Keep | Omit |
   |---|---|
   | `source_swp_4.2.0/` (59 Fortran files + `LICENSE`) | `bin/` (object files, build artifacts) |
   | `source_ttutil_4.27/` (source form; skip the `.zip` if both exist) | `linux/swap420` (Intel-compiled binary — useless elsewhere, bloat) |
   | `cases/` (test-case inputs: `.swp`, `.dra`, `.crp`, `.met`) | Any generated outputs inside `cases/*` |
   | `doc/` (manuals) | `.zip` archives (pick source form instead) |
   | `compiler_settings/` | |
   | `license/`, `readme_4.2.0.txt`, `xdata/` | |

   The curation is for repository size and portability, not license compliance.

   **Modernization repo license note.** Because the modernized tree in `src/` is a derivative work of the GPL v2 original, the modernization repo itself is bound by GPL v2. Verify that a `LICENSE` file exists at repository root declaring GPL v2 (or GPL v2+); add one if missing. Further implications are recorded in `docs/adr/0005-gpl-v2-inheritance.md` during Phase 2.
7. Delete `swaplib` and `swaplib-simple`.
8. Verify: `pixi run check-fast` green; `pixi run check-full` green at recorded tolerances (MOWDM oxygenstress deviation accepted; macropore ~344s accepted).
9. Commit this spec to `development` (`docs/superpowers/specs/2026-04-22-rescue-and-stabilize-design.md`).
10. Tag `rescue/phase-0-baseline`.

**Exit criterion:** `main` = `development` = `e256bc0` (plus the spec commit on `development`); `legacy/swap-4.2.0` orphan branch exists; `swaplib` and `swaplib-simple` gone; full-test green at documented tolerances.

### Phase 1 — Repository & infrastructure hygiene

**Goal:** One obvious way to do each build/test operation. Root directory is clean. Dependencies self-explanatory.

At baseline `e256bc0`, several of the drift-era problems do not yet exist (no `fpm.toml`, no `fpm_install/`, no `.bak` files in `src/crop/`, no `subroutines_to_inspect.md`, no `switch_to_gfort.md`). This phase is correspondingly smaller than it would have been starting from the drifted HEAD.

Targets:

| Issue | Action |
|---|---|
| Two meson build directories (`builddir` vs `builddir_gfortran`) | Unify into a single `builddir` with `enable_pfunit=true` as default. Update pixi tasks so there is exactly one `_configure` task and exactly one `build` task. |
| `tests/pFUnit/` full git clone in tree | Vendor pFUnit as a meson subproject (matching how `test-drive`, `toml-f`, `ttutil` are vendored) or pin it through pixi. Remove the local `.git` checkout. |
| `tests/swap-cases` submodule at detached HEAD | Pin to a tracked tag or branch; document update process in `docs/build-and-test.md`. |
| Hardcoded reference binary path (`swap_org/linux/swap420`) | Make path configurable through a pixi variable or task argument; document the override. |
| Missing `check-fast` / `check-full` pixi tasks | Add both. `check-fast` becomes the everyday command; `check-full` gates phase tags. |
| Optional pre-commit hook | Install a hook that runs `fprettify` on staged `.f90` files and `pixi run check-fast`. Bypassable with `--no-verify` but discouraged. |
| `tests/unit/` mixes `.pf` and standalone `.f90` files | Consolidate to pFUnit-only: one suite file per source module, all under `tests/unit/<domain>/test_<module>_suite.pf`. Regenerate `testSuites.inc` and add a meson check that `testSuites.inc` matches the set of `*_suite.pf` files on disk. |

**Stdlib / fpm note.** At baseline, neither `fortran-stdlib` nor `fpm` is referenced by `meson.build`. Both were drift-era additions. If a future phase decides we do need `stdlib`, bring it in via meson subproject or pixi dependency — not via a parallel fpm build.

**Exit criterion:** One build directory; `pixi run check-fast` completes under 90 seconds; `pixi run check-full` under 15 minutes; all tests under `tests/unit/**/*_suite.pf`; `testSuites.inc` regenerated and consistency-checked; root directory contains no WIP notes. Tag `rescue/phase-1-infra`.

### Phase 2 — Canonical documentation

**Goal:** A new contributor or future author can understand the project from `docs/` alone. No drift-era confusion.

FORD drives documentation generation. The project file is `docs.md` at repo root, and static pages live under `docs/`. FORD-generated API reference goes to `docs/api/` (generated output — not hand-edited). Phase 2 does not change the FORD driver; it rewrites the hand-authored pages `docs/` reads.

**Step 1 — Purge.** Delete the drift-era and duplicate pages:

- `docs/refactoring.md`
- `docs/refactoring-diary.md`
- `docs/next-refactoring-task.md` (its forward-looking content is rewritten as the follow-on compartment-states spec, not preserved here)
- `docs/newschema.md`
- `docs/state_management_pattern.md` (merged into the new `state-management.md`)
- `docs/state_management_architecture.md` (merged into the new `architecture.md` + `state-management.md`)
- `docs/style_guide.md` (duplicate of `code-style-guide.md`; keep one)

**Step 2 — Write canonical pages.** Structure `docs/` (hand-authored pages only; `docs/api/` and `docs/public/` are FORD-managed) as:

```
docs/
├── index.md                    Landing page (FORD root; links to the rest)
├── architecture.md             Three-phase control flow; aggregator state pattern; I/O layer; dependency graph
├── state-management.md         Canonical pattern: config_t / initial_t / state_t; ASSOCIATE; lifecycle
├── configuration-schema.md     TOML schema reference, per section
├── build-and-test.md           How to build; fast-vs-full test protocol; pFUnit conventions; regression harness
├── dependency-management.md    toml-f, ttutil, test-drive, pFUnit; how each is vendored; how to bump
├── code-style.md               Fortran 2008 conventions, naming, intent, module layout (single doc)
├── contributing.md             Branch model for this phase, commit convention, phase tags
├── adr/
│   ├── 0001-gfortran-first.md
│   ├── 0002-toml-over-namelist.md
│   ├── 0003-aggregator-state-over-compartments-for-now.md
│   ├── 0004-pfunit-for-unit-tests.md
│   └── 0005-gpl-v2-inheritance.md
└── superpowers/
    └── specs/                  This spec and future spec files
```

Each architecture decision record (ADR) is short: context, decision, consequences.

**Step 3 — Update `docs.md`** (the FORD project file) so its page ordering and `page_dir` settings match the new layout.

**Step 4 — Per-subdirectory READMEs.** Every `src/<domain>/` receives a `README.md` with exactly three sections:

1. *Responsibility* — one paragraph: what this subdir owns.
2. *Public interface* — list of public types and procedures intended for consumers.
3. *Dependencies* — which other `src/` subdirs it uses.

These READMEs are the single source of truth for "what does this directory do" and are the first thing a reader encounters. They stay short; detail lives in `docs/architecture.md` and source code comments.

**Step 5 — Agent definition rewrite.** `.github/agents/swap-fortran.agent.md` is rewritten to reflect this spec: pre-compartment gates, fast-vs-full test protocol, no legacy physics changes, pointer to this spec file. `swap-fortran-stage2.agent.md` and `fortran-repetition-hunter.agent.md` are deleted unless actively used.

**Exit criterion:** All top-level docs and per-subdir READMEs committed; `pixi run -e docs docs-build` succeeds and produces a complete FORD site; `docs/index.md` has working links to every top-level page and every ADR. Tag `rescue/phase-2-docs`.

### Phase 3 — Test coverage audit & expansion

**Goal:** Enough test coverage that Phase 4's TDD work is actually possible, with the crop consolidation and I/O consolidation safely protected.

**Step 1 — Coverage audit.** Add `gcov` / `lcov` to the `test` pixi feature. Run once; commit the baseline numbers (not the HTML) to `docs/coverage-baseline.md`.

**Step 2 — Categorize and fill.** For each `src/<domain>/`:

- **State lifecycle tests.** Every `*_state_t` needs `init` / `populate` / `finalize` tests.
- **I/O tests.** Every TOML reader needs a happy-path fixture test, at least one malformed-input test, and a round-trip test where applicable.
- **Physics unit tests.** Every *pure* routine gets a reference-value test against hand-computed expected output. Impure routines get *characterization* tests that lock in current behavior — not necessarily correct behavior — so Phase 4 refactors cannot silently change results.
- **Integration tests.** The four fast regression cases remain the top-level behavioral gate; they are the spine, not an afterthought.

**Step 3 — Target.** Coverage is **tracked, not gated**, in this phase. The target is "enough tests that Phase 4 refactors cannot break anything silently"; concrete percentages are a secondary indicator. A reasonable rough aim is ≥50% line coverage project-wide by phase exit, with `src/io/` and `src/core/` meaningfully higher because Phase 4 will churn them most, but the judgment call on coverage-vs-time rests with the author.

**Step 4 — Test organization pass.** By the end of this phase `tests/unit/` should mirror `src/` one-to-one: every subdirectory under `src/` has a matching subdirectory under `tests/unit/` with one `*_suite.pf` file per source module.

**Crop note.** Do not fix the crop duplication mess in Phase 3. Instead, write characterization tests that lock in whatever `cropfixed.f90`, `cropgrass.f90`, and `cropwofost.f90` currently do. The consolidation happens in Phase 4 under test protection.

**Exit criterion:** Coverage baseline recorded; `tests/unit/` mirrors `src/`; every `*_state_t` and every TOML reader has tests; crop characterization tests committed; full-test still green. Tag `rescue/phase-3-coverage`.

### Phase 4 — TDD fix-and-clean of main code

**Goal:** Close every pre-compartment gate under test protection. The development branch gets per-change feature branches for the first time.

Workflow changes starting here:

- Each ordered item below becomes a feature branch off `development` (e.g., `phase4/crop-consolidation`).
- Each feature branch merges back to `development` only when `pixi run check-full` is green.
- `development` merges fast-forward to `main` at the end of each numbered item if and only if `check-full` is green on `development`.

Ordered work items:

1. **Purge legacy `use variables` residue.** Enumerate every file still importing the legacy globals module; replace each with explicit state-type parameters. Start with leaf modules (`src/utils/`, `src/atmosphere/`) and work upward.
2. **Boundary conditions modernization.** Port `boundbottom.f90` and `boundtop.f90` from F77 style to the state-type / `intent(in)` pattern. `boundary_state_t` already exists; wire it through. Characterization tests from Phase 3 protect the change.
3. **I/O consolidation.** Commit to TOML as the single input path. Delete the legacy fixed-format readers (`readswap.f90`, `readdra.f90`, and the crop readers that are superseded by `readcrop_toml.f90`). Update `src/core/swap.f90` to import only the TOML readers. Update any regression case inputs still in fixed-format.
4. **Crop consolidation.** A single dispatcher in `cropgrowth.f90` delegating to one of three clean model files (`cropfixed.f90`, `cropgrass.f90`, `cropwofost.f90`), each contributing its own `*_state_t` fragment wired into `crop_state_t`. One reader (`readcrop_toml.f90`). Delete everything else in `src/crop/`. Characterization tests from Phase 3 protect the physics.
5. **Utils cleanup.** Convert F77-style collections (`sharedexchange.f90`, `sharedsimulation.f90`, `arrayutils.f90`, `ioutils.f90`) to modern modules with `private` defaults and `intent` attributes. Delete anything no longer referenced.
6. **Error handling.** Replace `fatalerr` calls and scattered status codes with a consistent pattern — either a populated `swap_error_t` flow or `error stop` plus `error_unit` writes (author picks one; documented in `docs/code-style.md`). Integrate throughout; `src/error/error.f90` stops being a stub.
7. **Pre-compartment gate closure sweep.** Verify each gate explicitly:
   - Every module has clean `config_t` / `initial_t` / `state_t` separation.
   - The TOML parser populates both `config` and `initial`.
   - The main program uses `init(state, config, initial)` universally.
   - Every process function has `config` as `intent(in)`.
   - Full regression green.

**Why not further splits.** Module interdependence in this codebase means splitting any of items 1–7 into smaller parallel pieces creates more merge friction than it saves. Keep them sequential.

**Exit criterion:** All seven items merged into `development` and fast-forwarded to `main`. Full-test green. Coverage not worse than Phase 3 baseline. Every gate above verified manually and recorded in `docs/architecture.md`. Tag `rescue/phase-4-clean` and `rescue/complete`.

## Risks and mitigations

| Risk | Mitigation |
|---|---|
| Something in `archive/wip-drifted` turns out to be needed | Cherry-pick into a Phase 4 feature branch with tests first; do not revive as "historical code" untested. |
| Phase 3 characterization tests lock in *incorrect* physics behavior | When Phase 4 reveals incorrect behavior, fix the test in a separate commit with a rationale comment, then fix the code. Never silently change both at once. |
| Full-test runtime grows as tests are added | Fast-test remains the <90s contract. Full-test budget is recorded at Phase 1 exit and revisited only with explicit author sign-off. |
| Crop consolidation (item 4.4) surfaces bugs that were masked by the duplication | Accept the scope growth; do not defer bug fixes past Phase 4. If a fix requires a physics decision the author cannot make quickly, revert the consolidation for that model and leave a `TODO` issue. |
| FORD page structure conflicts with proposed `docs/` layout | Phase 2 Step 3 explicitly adjusts `docs.md` to the new layout. If FORD has a hard constraint incompatible with the tree shown, adjust the tree rather than FORD config. |
| Author loses momentum during Phase 2 doc rewrite | Keep Phase 2 deliverables concrete (eight files + per-subdir READMEs + four ADRs) and time-boxed. If a page takes more than half a day, it is probably trying to say too much; split or cut. |

## Out of scope (follow-on specs)

1. **State-sync performance fix** — the 344s macropore regression. Separate spec to land between Phase 4 and the compartment-state refactor.
2. **Compartment-based state refactor** — replaces `docs/next-refactoring-task.md`. Introduces `ponding_state_t`, `canopy_state_t`, `snow_state_t`, `soil_state_t`, `saturated_zone_state_t`, `macropore_state_t`, `crop_state_t`, `drainage_state_t`, `irrigation_state_t`, `surface_water_state_t`, `boundary_reader_state_t`, `meteo_forcing_t`, `meteo_reader_state_t`. Flux functions become explicit `compute_*(from, to, flux)` triples.
3. **Python bindings + multicore** — likely `iso_c_binding` layer plus a thin Python wrapper (pyswap direction), then OpenMP across independent state instances.
4. **Vectorization + GPU** — targeting the spatial loops over compartments.
5. **Intel / ifx re-enablement** — restore the compile flags archived during the gfortran-only phase.

## Definition of done for this spec

- `main` tagged `rescue/complete` at the end of Phase 4.
- Every pre-compartment gate from the "Scope" section above verified.
- `pixi run check-full` green at documented tolerances (MOWDM oxygenstress deviation still accepted; macropore runtime still within its tracked envelope).
- `docs/` complete and rendered by FORD.
- `archive/wip-drifted`, `archive/main-pre-rescue`, and `legacy/swap-4.2.0` preserved in the repository.
- Zero `*.bak`, zero orphan `_new` files, zero root-level WIP notes, one meson build directory.
