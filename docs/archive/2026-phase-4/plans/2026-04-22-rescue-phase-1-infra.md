# Rescue Phase 1 — Repository & Infrastructure Hygiene Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the pFUnit target actually build and run green, lock the compiler to gfortran, consolidate build/test entry points into one obvious path, and record the two policy decisions (compiler, build-dir topology) as ADRs — so Phase 3 (test coverage) has a stable foundation to work on.

**Architecture:** Tight scope, ruthless cut of nice-to-haves. Remove all Intel/ifx branches from `meson.build`, reconcile `tests/unit/meson.build` with source files that actually exist at the `e256bc0` baseline (deleting drift-era forward-declarations rather than creating stubs), unify the two separate meson build directories into one with pFUnit enabled by default, and introduce two orthogonal pixi tasks (`check-fast`, `check-full`) that become the permanent iteration/verification boundary.

**Tech Stack:** gfortran, meson, ninja, pixi, pFUnit 4.15, python 3.11 (for regression harness), git.

**Spec:** [docs/superpowers/specs/2026-04-22-rescue-and-stabilize-design.md](../specs/2026-04-22-rescue-and-stabilize-design.md) — Phase 1

**Starting state at execution:**
- On `development` branch, HEAD = `3a65acf` (Phase 0 tail: swap_org removal + reference binaries preservation). `main` fast-forwarded to the same commit. `origin/main` at `7587ca3`, untouched.
- Tag `rescue/phase-0-baseline` at `ad4411d`.
- Single `builddir/` at root (created by Phase 0 verification build; gitignored).
- pFUnit target broken: `pixi run test-pfunit` fails meson configure because `tests/unit/meson.build` references source files that don't exist at baseline.
- Full regression passes 5/6 at documented tolerances (~343s macropore, MOWDM deviation on oxygenstress).

**End state after this plan:**
- `meson.build` locked to gfortran (`-ffree-line-length-none`, `-std=legacy`, no Intel code paths); Intel/ifx branches removed.
- `tests/unit/meson.build` references only sources that exist at baseline; pFUnit configure succeeds; at least one pFUnit test runs and passes.
- A single `builddir/` used for both the production `swap` binary and pFUnit tests (`enable_pfunit` defaults to `true`); pixi's `_configure-pfunit` and `builddir_gfortran` eliminated.
- Two new pixi tasks: `check-fast` (build + pFUnit + 4 fast regression cases, ≤90s target) and `check-full` (same plus oxygenstress and macropore, ~10 min).
- `.github/agents/swap-fortran.agent.md` rewritten to reflect this spec; `swap-fortran-stage2.agent.md` and `fortran-repetition-hunter.agent.md` deleted.
- Two ADRs: `docs/adr/0001-gfortran-first.md` (compiler policy) and `docs/adr/0002-single-builddir.md` (build topology).
- Tag `rescue/phase-1-infra` on `main` = `development`.

**Deferred out of Phase 1 (tracked for Phase 1.5 or Phase 2):**
- Full pFUnit re-vendoring as a meson subproject (only a short `tests/pFUnit/README.md` note in this plan).
- Pre-commit hook installation (spec marks optional; author can add manually anytime).
- Auto-regeneration check for `testSuites.inc` (manual regeneration in this plan; auto-check deferred to Phase 3 when coverage tooling lands).
- Root `LICENSE` GPL v2 correction — Phase 2 ADR 0005.
- `tests/swap-cases` submodule dirtying-on-regression behavior documentation — Phase 2 `docs/build-and-test.md`.

**Conventions:**
- One concern per commit. Conventional-commit-style messages (`chore:`, `fix:`, `docs:`, `refactor:`, `feat:`).
- Every task ends by verifying green state with either `pixi run build-linux`, `pixi run test-pfunit` (until Task 6 replaces it), or `pixi run check-fast` (after Task 6 introduces it).
- `git status --short` before each commit — nothing unexpected staged.
- No pushes to any remote during Phase 1.

---

## Task 1: Lock compiler to gfortran

**Goal:** Remove all Intel/ifx code paths from `meson.build` so there is exactly one way to build: gfortran via the pixi-managed toolchain. This decision is committed in Task 8's ADR; this task just executes it in the build config.

**Files:**
- Modify: `/home/zawadzkim/Code/swap/meson.build`

### Step 1: Inspect the current compiler-branching structure

Run:
```bash
grep -nE "is_intel|is_gcc|fc.get_id" /home/zawadzkim/Code/swap/meson.build
```

Expected: lines referencing `is_intel`, `is_gcc`, `fc.get_id()` at known locations (roughly lines 13–20 and 77–97 in the current file). If the line numbers diverge significantly from what this plan expects (±5 lines), read the file first and adjust the edits accordingly.

### Step 2: Replace the compiler-detection block with a gfortran-only guard

Find this block in `meson.build` (around lines 13–20):

```
fc = meson.get_compiler('fortran')
fs = import('fs')

src_inc = include_directories(
    'src', 'src/core', 'src/io', 'src/soil', 'src/atmosphere',
    'src/crop', 'src/drainage', 'src/boundary', 'src/macropore',
    'src/solute', 'src/heat', 'src/utils'
)

# Platform/compiler detection
is_windows = host_machine.system() == 'windows'
is_intel = fc.get_id() == 'intel'
is_gcc = fc.get_id() == 'gcc'
```

Replace with:

```
fc = meson.get_compiler('fortran')
fs = import('fs')

src_inc = include_directories(
    'src', 'src/core', 'src/io', 'src/soil', 'src/atmosphere',
    'src/crop', 'src/drainage', 'src/boundary', 'src/macropore',
    'src/solute', 'src/heat', 'src/utils'
)

# Platform detection. The rescue spec locks the supported compiler
# to GCC/gfortran (see docs/adr/0001-gfortran-first.md). Reject any
# other compiler up front rather than silently producing a different
# binary.
is_windows = host_machine.system() == 'windows'

if fc.get_id() != 'gcc'
    error('Unsupported Fortran compiler: @0@. Only gfortran is supported during the rescue (see docs/adr/0001-gfortran-first.md).'.format(fc.get_id()))
endif
```

### Step 3: Replace the compiler-flags block with gfortran-only flags

Find the block (around lines 22–46 previously — the `if is_intel ... elif is_gcc ... endif` section):

```
if is_intel
    add_project_arguments([
        '-O2', '-fpscomp', 'general', '-warn', 'declarations',
        '-warn', 'unused', '-warn', 'uncalled', '-warn', 'interfaces',
        '-init=zero', '-fpe0', '-fp-model', 'source'
    ], language: 'fortran')
    
    if not is_windows
        add_project_link_arguments(['-static-intel'], language: 'fortran')
    endif

elif is_gcc
    gcc_flags = [
        '-O2', '-ffree-line-length-none', '-Wno-line-truncation',
        '-Wno-compare-reals', '-Wno-missing-include-dirs',
        '-Wno-unused-variable', '-Wno-unused-dummy-argument', '-std=legacy'
    ]
    
    if is_windows
        gcc_flags += ['-mconsole', '-static-libgfortran', '-static-libgcc']
        add_project_link_arguments(['-static', '-mconsole'], language: 'fortran')
    endif
    
    add_project_arguments(gcc_flags, language: 'fortran')
endif
```

Replace with:

```
gfortran_flags = [
    '-O2',
    '-ffree-line-length-none',
    '-Wno-line-truncation',
    '-Wno-compare-reals',
    '-Wno-missing-include-dirs',
    '-Wno-unused-variable',
    '-Wno-unused-dummy-argument',
    '-std=legacy',
]

if is_windows
    gfortran_flags += ['-mconsole', '-static-libgfortran', '-static-libgcc']
    add_project_link_arguments(['-static', '-mconsole'], language: 'fortran')
endif

add_project_arguments(gfortran_flags, language: 'fortran')
```

### Step 4: Simplify the ttutil dependency block

Find the block in `meson.build` starting at the `# External dependency: ttutil` section (roughly line 137–108):

```
# ============================================================================
# External dependency: ttutil
# ============================================================================
if is_gcc
    # Use ttutil as a Meson subproject (built from source with gfortran)
    ttutil_prj = subproject('ttutil', default_options: ['default_library=static'])
    ttutil_dep = ttutil_prj.get_variable('ttutil_dep')
else
    # Use prebuilt library for Intel / Windows
    ttutil_version = '4.2.7'
    ttutil_lib_dir = join_paths(meson.current_build_dir(), 'lib')

    if is_windows and not is_gcc
        ttutil_filename = 'libttutil@0@-windows.a'.format(ttutil_version)
    elif is_windows and is_gcc
        ttutil_filename = 'libttutil@0@-mingw.a'.format(ttutil_version)
    else
        ttutil_filename = 'libttutil@0@-linux.a'.format(ttutil_version)
    endif

    ttutil_lib_path = join_paths(ttutil_lib_dir, 'libttutil.a')
    ttutil_url = 'https://github.com/SWAP-model/ttutil/releases/download/v@0@/@1@'.format(
        ttutil_version, ttutil_filename)

    if not fs.exists(ttutil_lib_path)
        run_command('mkdir', '-p', ttutil_lib_dir, check: true)
        run_command('curl', '-L', ttutil_url, '-o', ttutil_lib_path, check: true)
    endif

    ttutil_dep = declare_dependency(link_args: [ttutil_lib_path])
endif
```

Replace with:

```
# ============================================================================
# External dependency: ttutil
# ============================================================================
# Vendored as a meson subproject and built from source with gfortran.
# The Intel/ifx download-prebuilt-binary path has been removed; see
# docs/adr/0001-gfortran-first.md.
ttutil_prj = subproject('ttutil', default_options: ['default_library=static'])
ttutil_dep = ttutil_prj.get_variable('ttutil_dep')
```

### Step 5: Simplify the final-executable link_args

Find the block near the bottom of `meson.build`:

```
# ============================================================================
# Main executable
# ============================================================================
if (is_windows and is_gcc) or not is_windows
    link_args = ['-static']
else
    link_args = []
endif
```

Since we are now gcc-only, simplify to:

```
# ============================================================================
# Main executable
# ============================================================================
# Always statically link with gfortran — see docs/adr/0001-gfortran-first.md.
link_args = ['-static']
```

### Step 6: Verify that no other references to `is_intel` remain

Run:
```bash
grep -nE "is_intel|fpscomp|init=zero|fp-model" /home/zawadzkim/Code/swap/meson.build
```

Expected: zero matches. If any remain, delete them — they are Intel-specific leftovers.

### Step 7: Rebuild and verify the production binary

Run:
```bash
rm -rf builddir && pixi run build-linux 2>&1 | tail -20
ls -la builddir/swap
file builddir/swap
```

Expected:
- Meson configure succeeds; notes the compiler as `GCC gfortran` (not `Intel ifx`).
- Ninja compile succeeds; produces `builddir/swap`.
- `file` identifies the binary as an ELF 64-bit statically linked executable.

If configure fails complaining about the compiler, run:
```bash
which gfortran
gfortran --version
```
Confirm gfortran is on `$PATH` inside pixi. If not, investigate the pixi environment — this is a blocker and must be resolved before proceeding.

### Step 8: Run fast regression (4 cases) to confirm behavior is preserved

Run:
```bash
time pixi run regression 2>&1 | tee /tmp/phase1-task1-regression.log | tail -30
```

Expected: same pattern as Phase 0 baseline — 5/6 pass, `oxygenstress` fails only on MOWDM, total ~10 min. If the failure pattern changes in any other way, STOP — the compiler change has broken physics and needs investigation.

### Step 9: Commit

Run:
```bash
git add meson.build
git status --short
git commit -m "chore(build): lock fortran compiler to gfortran only

Remove the Intel/ifx code paths from meson.build, including the
prebuilt-binary download fallback for ttutil. Any non-GCC compiler
now fails configure with a clear error pointing at the rescue ADR.

Consolidates three things:
- Compiler-detection block replaced with an explicit guard.
- Flags block collapsed to a single gfortran_flags list with a
  windows branch for static linking.
- ttutil dependency always resolves via the meson subproject; the
  Intel/Windows prebuilt-download path is removed.
- Final executable always gets -static link_args.

Regression suite still passes 5/6 at documented tolerances
(oxygenstress MOWDM and macropore ~343s pre-existing). No physics
change. Full rationale in docs/adr/0001-gfortran-first.md (added
in task 8 of this plan)."
```

---

## Task 2: Reconcile tests/unit/meson.build with baseline sources

**Goal:** Make `pixi run test-pfunit` configure succeed by editing `tests/unit/meson.build` to reference only source files that actually exist at `e256bc0`. Delete orphan test .f90 and .pf files whose source dependencies do not exist.

**Files:**
- Modify: `/home/zawadzkim/Code/swap/tests/unit/meson.build`
- Delete: `/home/zawadzkim/Code/swap/tests/unit/io/test_readcrop_toml.f90`
- Delete: `/home/zawadzkim/Code/swap/tests/unit/io/test_readcrop_sectionreaders.f90`
- Delete: `/home/zawadzkim/Code/swap/tests/unit/io/test_readdra_toml.f90`
- Possibly modify: `/home/zawadzkim/Code/swap/tests/unit/testSuites.inc` (in Task 4)

### Step 1: Inventory existing source files at baseline

Run these to produce a ground-truth map the subagent consults during editing:

```bash
echo "=== src/core ==="; ls /home/zawadzkim/Code/swap/src/core/ | sort
echo "=== src/io ==="; ls /home/zawadzkim/Code/swap/src/io/ | sort
echo "=== src/atmosphere ==="; ls /home/zawadzkim/Code/swap/src/atmosphere/ | sort
echo "=== src/heat ==="; ls /home/zawadzkim/Code/swap/src/heat/ | sort
echo "=== src/crop ==="; ls /home/zawadzkim/Code/swap/src/crop/ | sort
echo "=== src/boundary ==="; ls /home/zawadzkim/Code/swap/src/boundary/ | sort
echo "=== src/drainage ==="; ls /home/zawadzkim/Code/swap/src/drainage/ | sort
echo "=== src/solute ==="; ls /home/zawadzkim/Code/swap/src/solute/ | sort
echo "=== src/soil ==="; ls /home/zawadzkim/Code/swap/src/soil/ | sort
echo "=== src/macropore ==="; ls /home/zawadzkim/Code/swap/src/macropore/ | sort
echo "=== src/utils ==="; ls /home/zawadzkim/Code/swap/src/utils/ | sort
echo "=== tests/unit/io ==="; ls /home/zawadzkim/Code/swap/tests/unit/io/
echo "=== tests/unit/state_lifecycles ==="; ls /home/zawadzkim/Code/swap/tests/unit/state_lifecycles/ 2>/dev/null
echo "=== tests/unit/utils ==="; ls /home/zawadzkim/Code/swap/tests/unit/utils/ 2>/dev/null
```

Record the output. The expected baseline inventory is (from prior reconnaissance):

- `src/io/`: `macroporeoutput.f90`, `readdrainagetoml.f90`, `readmeteo.f90`, `readswap.f90`, `readswaptoml.f90`, `swap_csv_output.f90`, `swapoutput.f90`.
- State modules present: `atmosphere_state.f90`, `boundary_state.f90`, `drainage_state.f90`, `surfacewater_state.f90`, `heat_state.f90`, `solute_state.f90`, `soil_state.f90`, `macropore_state.f90`, `cropgrowth_state.f90`.
- `src/core/`: includes `swap_state_mod.f90` (NOT `swap_state.f90`), `variables.f90`, `swap.f90`, `swap_main.f90`, plus others.
- **Not present at baseline**: `readcrop_toml.f90`, `readdra_toml.f90` (old name), `readswap_toml.f90` (old name), `time_state.f90`, `snow_state.f90`, `frozencond_state.f90`, `irrigation_state.f90`, `tillage_state.f90`, `oxygenstress_state.f90`, `management_soil_state.f90`, `rootextraction_state.f90`, `wofost_state.f90`, `simulation_config.f90`, `output_config.f90`, `swap_config.f90`, `swap_io.f90`, `swap_state.f90`.

If the actual listing diverges from this — for example if some of the "not present" modules turn out to exist — adjust the edits below accordingly. The principle is: reference only what exists.

### Step 2: Replace `tests/unit/meson.build` with a baseline-accurate version

Overwrite the file entirely with the following content. This rewrite: (a) keeps only sources that exist at baseline, (b) drops the three orphan standalone test executables, (c) shrinks the pFUnit source list to what the surviving `.pf` files actually need, (d) adjusts reader source names from the drift-era (`readdra_toml`, `readswap_toml`) to the baseline names (`readdrainagetoml`, `readswaptoml`).

```meson
# ============================================================================
# Unit test build configuration
# ============================================================================
# Only sources that exist at the rescue baseline (commit e256bc0) are
# referenced here. Drift-era forward declarations for state modules and
# TOML readers that were never actually committed have been removed.
#
# When Phase 4 re-adds the missing state modules and readers, this file
# grows back; the Phase 4 plan explicitly adds them as part of each
# module's TDD sub-task.
# ============================================================================

test_base_sources = [
    '../../src/core/arrays.f90',
    '../../src/core/constants.f90',
    '../../src/core/swap_log.f90',
    '../../src/utils/arrayutils.f90',
    '../../src/atmosphere/atmosphere_state.f90',
    '../../src/boundary/boundary_state.f90',
    '../../src/drainage/drainage_state.f90',
    '../../src/drainage/surfacewater_state.f90',
    '../../src/heat/heat_state.f90',
    '../../src/solute/solute_state.f90',
    '../../src/soil/soil_state.f90',
    '../../src/macropore/macropore_state.f90',
    '../../src/crop/cropgrowth_state.f90',
]

# ============================================================================
# pFUnit test suite (gfortran only)
# ============================================================================
if get_option('enable_pfunit')
    python3 = import('python').find_installation('python3')

    pfunit_root_result = run_command('sh', '-c', 'echo ${PFUNIT_ROOT:-}', check: false)
    pfunit_root = pfunit_root_result.stdout().strip()

    if pfunit_root == ''
        pfunit_root = join_paths(meson.current_source_dir(),
            '../../tests/pFUnit/build/install_gfortran/PFUNIT-4.15')
        if not fs.exists(pfunit_root)
            error('pFUnit not found. Set PFUNIT_ROOT to a PFUNIT-4.15 install, or build one under tests/pFUnit/build/install_gfortran/PFUNIT-4.15. See tests/pFUnit/README.md.')
        endif
    endif

    message('Using pFUnit from: ' + pfunit_root)

    pfunit_dep = dependency('PFUNIT',
        method: 'cmake',
        required: true,
        cmake_args: ['-DCMAKE_PREFIX_PATH=' + pfunit_root]
    )

    funitproc = join_paths(pfunit_root, 'bin', 'funitproc')
    if not fs.exists(funitproc)
        error('funitproc not found at: ' + funitproc)
    endif

    pfunit_pp = generator(python3,
        output: '@BASENAME@.F90',
        arguments: [funitproc, '@INPUT@', '@OUTPUT@']
    )

    driver_src = join_paths(pfunit_root, 'include', 'driver.F90')

    # pFUnit suites that exist at baseline and whose referenced source
    # modules also exist. Suites relying on TOML readers or state modules
    # that are not present are commented out with a Phase 4 reference.
    pfunit_base_sources = test_base_sources + [
        '../../src/utils/tomlutils.f90',
        '../../src/io/readdrainagetoml.f90',
        '../../src/io/readswaptoml.f90',
    ]

    pf_files = files(
        'io/test_readswaptoml_full.pf',
        'io/test_sectionreaders_suite.pf',
        'state_lifecycles/test_state_lifecycles_suite.pf',
    )
    pp_sources = pfunit_pp.process(pf_files)

    unit_tests = executable('unit-swap-tests',
        sources: [pp_sources, driver_src] + pfunit_base_sources,
        dependencies: [pfunit_dep, ttutil_dep, tomlf_dep],
        include_directories: src_inc,
        fortran_args: ['-D_TEST_SUITES="testSuites.inc"'],
        install: false
    )

    test('unit-swap-tests', unit_tests,
        suite: 'unit-pfunit',
        workdir: meson.project_source_root()
    )
endif
```

If `src/utils/tomlutils.f90` does not exist at baseline — verify with `ls src/utils/tomlutils.f90` — remove that line. If the existing `.pf` files reference a module that isn't in the reduced source list, the pFUnit configure will fail later and the subagent adjusts iteratively.

### Step 3: Delete orphan standalone test .f90 files

Run:
```bash
cd /home/zawadzkim/Code/swap
rm tests/unit/io/test_readcrop_toml.f90
rm tests/unit/io/test_readcrop_sectionreaders.f90
rm tests/unit/io/test_readdra_toml.f90
ls tests/unit/io/
```

Expected: after removal, `tests/unit/io/` contains only `test_readswaptoml_full.pf` and `test_sectionreaders_suite.pf` (the two `.pf` files that survive).

### Step 4: Verify meson configure now succeeds with pFUnit enabled

Run:
```bash
cd /home/zawadzkim/Code/swap
rm -rf builddir
pixi run _configure-pfunit 2>&1 | tail -40
```

Expected: meson configure succeeds. The message `Using pFUnit from: ...` appears. No `ERROR: File ... does not exist` lines.

If configure fails because one of the `.pf` files USE-imports a module not in `pfunit_base_sources`, inspect the failing `.pf` file, add the missing module to `pfunit_base_sources` if its source exists at baseline, or DELETE the `.pf` file entry from `pf_files` if it depends on missing source. Re-run configure until green. Each adjustment is a single edit — keep iterating.

### Step 5: Compile the pFUnit test executable

Run:
```bash
cd /home/zawadzkim/Code/swap
pixi run _compile-pfunit 2>&1 | tail -40
```

Expected: compiles cleanly. Some warnings are acceptable; no errors. Produces an executable at `builddir_gfortran/tests/unit/unit-swap-tests` (or similar).

If a `.pf` file's compiled output references a symbol that doesn't link, the unresolved symbol name tells you which module is missing from `pfunit_base_sources`. Add it or drop the `.pf` file from the list.

### Step 6: Run the pFUnit suite

Run:
```bash
cd /home/zawadzkim/Code/swap
pixi run test-pfunit 2>&1 | tail -60
```

Expected: at least one suite runs. Pass/fail of individual tests is acceptable — what matters is the process gets past configure and compile to actual execution. Record the observed output for the ADR 0002 commit message in Task 10.

If individual tests FAIL (not configure-fail), that's a Phase 3 task. Note the failures and move on — Phase 1's exit is "tests build and run", not "tests all pass".

### Step 7: Commit

Run:
```bash
cd /home/zawadzkim/Code/swap
git add tests/unit/meson.build
git rm tests/unit/io/test_readcrop_toml.f90 \
       tests/unit/io/test_readcrop_sectionreaders.f90 \
       tests/unit/io/test_readdra_toml.f90
git status --short
git commit -m "fix(tests): reconcile tests/unit/meson.build with baseline sources

The unit-test meson config was aspirational — it listed ~15 source
files that were never actually committed at e256bc0. This caused
pixi run test-pfunit to fail at meson configure, blocking Phase 3
coverage work.

Rewrite the file to reference only sources present at the rescue
baseline:
- test_base_sources reduced to nine state modules + arrays +
  constants + swap_log + arrayutils.
- pFUnit source list shrunk similarly; renames readdra_toml ->
  readdrainagetoml and readswap_toml -> readswaptoml to match the
  actual file names at baseline.
- Three standalone readcrop_* / readdra_toml test .f90 files
  deleted; they tested sources that do not exist.

Phase 4 will re-add the missing state modules and readers as part
of each module's TDD sub-task, and the corresponding test entries
come back then."
```

---

## Task 3: Regenerate testSuites.inc to match surviving pFUnit suites

**Goal:** Ensure `tests/unit/testSuites.inc` lists exactly the test suites that survived Task 2. pFUnit's test driver reads this file to register suites.

**Files:**
- Modify: `/home/zawadzkim/Code/swap/tests/unit/testSuites.inc`

### Step 1: Inspect current testSuites.inc

Run:
```bash
cat /home/zawadzkim/Code/swap/tests/unit/testSuites.inc
```

Expected (current): three `ADD_TEST_SUITE(...)` lines for `test_readswaptoml_full_suite`, `test_sectionreaders_suite_suite`, `test_state_lifecycles_suite_suite`.

### Step 2: Determine which suites the surviving `.pf` files define

Each `.pf` file defines a module whose name ends in `_suite`. The `ADD_TEST_SUITE` name is the module name. Run:

```bash
grep -l "module\s" /home/zawadzkim/Code/swap/tests/unit/**/*.pf 2>/dev/null
grep -E "^\s*module\s+" /home/zawadzkim/Code/swap/tests/unit/**/*.pf 2>/dev/null
```

Expected: modules named `test_readswaptoml_full_suite`, `test_sectionreaders_suite_suite`, `test_state_lifecycles_suite_suite` (possibly also `test_utils_suite_suite` if `tests/unit/utils/` contains a `.pf` file).

If any `.pf` file found here is NOT listed in the `pf_files` block of `tests/unit/meson.build` as edited in Task 2, either add it there or delete it from disk — it must not be orphan.

### Step 3: Write testSuites.inc matching the surviving suites

If the four-or-fewer surviving suites are the three already in the file — no change needed. Verify and move on.

If `tests/unit/utils/test_utils_suite.pf` exists and is referenced in `tests/unit/meson.build`, add its `ADD_TEST_SUITE` line:

```
ADD_TEST_SUITE(test_readswaptoml_full_suite)
ADD_TEST_SUITE(test_sectionreaders_suite_suite)
ADD_TEST_SUITE(test_state_lifecycles_suite_suite)
ADD_TEST_SUITE(test_utils_suite_suite)
```

Otherwise keep the current three-line content.

### Step 4: Rebuild and rerun pFUnit to confirm suites register

Run:
```bash
cd /home/zawadzkim/Code/swap
rm -rf builddir_gfortran
pixi run test-pfunit 2>&1 | tail -30
```

Expected: the configure + compile + run flow completes. Every registered suite runs at least one test (or reports 0 tests cleanly).

### Step 5: Commit (only if testSuites.inc changed)

If Step 3 made changes:
```bash
git add tests/unit/testSuites.inc
git commit -m "fix(tests): regenerate testSuites.inc for surviving pFUnit suites

Aligns the test-driver suite registration with the pf_files list
in tests/unit/meson.build after Task 2 reconciled the latter with
baseline sources."
```

If no change was needed, skip the commit. Do NOT create an empty commit.

---

## Task 4: Unify meson build directories

**Goal:** Replace the two separate build directories (`builddir` for production, `builddir_gfortran` for tests) with a single `builddir` that always builds the production binary and optionally the pFUnit tests. Change `enable_pfunit` default to `true`.

**Files:**
- Modify: `/home/zawadzkim/Code/swap/meson_options.txt`
- Modify: `/home/zawadzkim/Code/swap/pixi.toml` (setup, build, and test tasks)

### Step 1: Flip the pFUnit default to true

Run:
```bash
cat /home/zawadzkim/Code/swap/meson_options.txt
```

Expected current content:
```
option('enable_pfunit', type: 'boolean', value: false, description: 'Build pFUnit-based unit tests')
option('enable_unit_tests', type: 'boolean', value: false, description: 'Build standalone unit test executables')
```

Replace with:
```
option('enable_pfunit', type: 'boolean', value: true, description: 'Build pFUnit-based unit tests (default: on).')
option('enable_unit_tests', type: 'boolean', value: false, description: 'Build standalone unit test executables (default: off; kept for occasional non-pFUnit tests).')
```

### Step 2: Unify the pixi setup / build / test task tree

Read the current pixi.toml tasks block to locate the surrounding context:

```bash
sed -n '70,110p' /home/zawadzkim/Code/swap/pixi.toml
```

The current task graph has (roughly):
- `_fetch-deps` → `_configure-swap` (production, no pFUnit, `builddir`) → `build-linux`.
- `_fetch-deps` → `_configure-pfunit` (tests, `builddir_gfortran`) → `_compile-pfunit` → `test-pfunit`.

Replace this two-track graph with a single-track one. Edit `pixi.toml` so the tasks section becomes (preserve surrounding block comments and RUN TASKS / LINT TASKS / DOC TASKS sections unchanged):

```toml
#============================================================================
# BUILD TASKS — single builddir; pFUnit on by default.
#============================================================================
_fetch-deps    = { cmd = "fpm install --prefix fpm_install", outputs = ["fpm_install"] }
_configure     = { cmd = "meson setup builddir", depends-on = ["_fetch-deps"] }
build-linux    = { cmd = "meson compile -C builddir", depends-on = ["_configure"] }

# Tests compile as part of the main build now. test-pfunit runs them.
test-pfunit    = { cmd = "meson test -C builddir --suite unit-pfunit --print-errorlogs", depends-on = ["build-linux"] }
```

The `_configure-swap`, `_configure-pfunit`, and `_compile-pfunit` task entries are removed entirely. `builddir_gfortran/` will never be created again.

Any other task that referenced `builddir_gfortran` must be updated to `builddir`. Search and fix:

```bash
grep -n "builddir_gfortran" /home/zawadzkim/Code/swap/pixi.toml
```

Expected after Step 2 edits: zero matches. If any remain, edit them to `builddir` one by one.

### Step 3: Update the `PFUNIT_ROOT` env var if it references `builddir_gfortran`

Run:
```bash
grep -n "PFUNIT_ROOT\|pFUnit" /home/zawadzkim/Code/swap/pixi.toml
```

If `PFUNIT_ROOT` points at a path that only made sense in the old two-track layout, leave it as-is IF the path is actually the pFUnit *install* dir (not the meson build dir). `tests/pFUnit/build/install_gfortran/PFUNIT-4.15` is the install dir — that's fine, leave it. Only change it if it points into `builddir_gfortran/`.

### Step 4: Wipe old build dirs and verify

Run:
```bash
cd /home/zawadzkim/Code/swap
rm -rf builddir builddir_gfortran
pixi run build-linux 2>&1 | tail -20
```

Expected: configure completes with pFUnit enabled (visible in meson's output), compile succeeds, `builddir/swap` exists.

```bash
ls builddir/swap
ls builddir | head
```

Expected: `swap` binary plus test executables inside `builddir/`.

### Step 5: Run pFUnit via the unified task

```bash
pixi run test-pfunit 2>&1 | tail -30
```

Expected: at least one test suite runs. If Step 5 of Task 2 captured failing tests, they still fail here — that's fine; Phase 1 doesn't try to fix Phase 3's work.

### Step 6: Run fast regression

```bash
time pixi run regression 2>&1 | tail -20
```

Expected: 5/6 pass, MOWDM deviation on oxygenstress unchanged, ~10 min total. Physics still intact.

### Step 7: Commit

```bash
git add meson_options.txt pixi.toml
git status --short
git commit -m "refactor(build): unify on a single builddir

Phase 0 inherited two meson build directories — builddir for
production and builddir_gfortran for pFUnit tests — with separate
configure steps and a convention that tests get built only when
an operator remembered to run a test-specific task. This was a
foot-gun.

One builddir now builds everything. meson_options.enable_pfunit
defaults to true; pixi collapses _configure-swap and
_configure-pfunit into a single _configure task. test-pfunit
just runs the relevant suite in the one build tree.

builddir_gfortran is no longer created anywhere. Rationale and
migration notes in docs/adr/0002-single-builddir.md (added in
task 8)."
```

---

## Task 5: Add check-fast and check-full pixi tasks

**Goal:** Two canonical verification commands that Phase 2–4 work reference. `check-fast` is the everyday iteration gate; `check-full` is the phase-end gate.

**Files:**
- Modify: `/home/zawadzkim/Code/swap/pixi.toml` (tasks + add fast/full regression split)
- Possibly: `/home/zawadzkim/Code/swap/tests/regression/test_output_regression.py` — inspect to determine how to filter cases.

### Step 1: Understand how the regression harness selects cases

Run:
```bash
head -80 /home/zawadzkim/Code/swap/tests/regression/test_output_regression.py
```

Look for how cases are collected (likely from `tests/swap-cases/*/`) and whether there's a CLI flag or env-var to filter.

If the harness accepts a `--cases` / `-k` pytest filter, use that. If it iterates over every subdirectory regardless, Step 2 adds filtering.

### Step 2: Add a case-filtering option to the regression harness

If the harness already supports filtering, skip this step.

Otherwise, edit the harness to accept a `--cases` argument or equivalent. The minimal addition:

```python
# near the existing CLI argparse block, or at the top if none exists
import argparse
import sys

FAST_CASES = {"hupselbrook", "surfacewater", "salinitystress", "grassgrowth"}
ALL_CASES = FAST_CASES | {"oxygenstress", "macropore"}

def _parse_argv(argv):
    p = argparse.ArgumentParser()
    p.add_argument("--cases", choices=["fast", "full"], default="full",
                   help="fast = hupselbrook+surfacewater+salinitystress+grassgrowth; full = all six.")
    return p.parse_args(argv)
```

Use the parsed selection to filter whatever list of cases the harness iterates over. Do NOT change output-comparison logic.

If the harness is more complex than a simple loop, or if modifying Python is out of scope for the subagent, fall back to a wrapper: add a small pixi-level shell that exports `SWAP_REGRESSION_CASES=fast` and have the harness read that env var at startup. Either approach is acceptable — pick the one with fewer code changes.

### Step 3: Add the pixi task definitions

Edit `pixi.toml` to add (in the RUN TASKS or a new VERIFICATION TASKS section):

```toml
#============================================================================
# VERIFICATION TASKS — the canonical gates.
# check-fast: everyday iteration. ~90s budget (build + pFUnit + 4 fast cases).
# check-full: phase-end gate. ~10 min budget (check-fast + oxygenstress + macropore).
#============================================================================
check-fast = { cmd = "python tests/regression/test_output_regression.py --cases fast", depends-on = ["build-linux", "test-pfunit"] }
check-full = { cmd = "python tests/regression/test_output_regression.py --cases full", depends-on = ["build-linux", "test-pfunit"] }
```

If Step 2 took the env-var fallback, use `SWAP_REGRESSION_CASES=fast python ...` instead.

### Step 4: Run check-fast and verify budget

```bash
cd /home/zawadzkim/Code/swap
time pixi run check-fast 2>&1 | tail -40
```

Expected:
- Build + pFUnit + 4 fast regression cases all run.
- Wall time well under 90 seconds. If over, re-examine whether build is truly incremental or whether pFUnit ran unnecessary tests.
- Regression summary shows 4/4 passing (fast cases don't include the MOWDM deviation case).

### Step 5: Run check-full and confirm baseline equivalence

```bash
time pixi run check-full 2>&1 | tail -40
```

Expected:
- All 6 cases run; 5 pass; oxygenstress fails only on MOWDM.
- Wall time ~10 minutes.

### Step 6: Commit

```bash
git add pixi.toml tests/regression/test_output_regression.py
git status --short
git commit -m "feat(tests): add check-fast and check-full pixi gates

These two tasks become the canonical verification boundary for the
rest of the rescue and ongoing work:

- check-fast: build + pFUnit + four fast regression cases
  (hupselbrook, surfacewater, salinitystress, grassgrowth).
  Target budget < 90s; runs before every phase commit.
- check-full: same plus oxygenstress and macropore. Target ~10min;
  runs before every phase tag and before any merge to main.

Regression harness accepts --cases {fast,full} (or
SWAP_REGRESSION_CASES env var) to filter the case set."
```

---

## Task 6: Minimal pFUnit documentation

**Goal:** Leave a single-screen README at `tests/pFUnit/` explaining the current state, how to rebuild if the install is lost, and that full vendoring (as a meson subproject) is a deferred Phase 1.5/2 task. Do NOT remove the existing install.

**Files:**
- Create: `/home/zawadzkim/Code/swap/tests/pFUnit/README.md`

### Step 1: Inspect current state

```bash
ls /home/zawadzkim/Code/swap/tests/pFUnit/ | head
du -sh /home/zawadzkim/Code/swap/tests/pFUnit/
test -e /home/zawadzkim/Code/swap/tests/pFUnit/.git && echo "has .git"
test -e /home/zawadzkim/Code/swap/tests/pFUnit/build/install_gfortran/PFUNIT-4.15 && echo "has PFUNIT-4.15 install"
```

Expected: full git checkout + a built install at the expected path. Record the total size.

### Step 2: Write the README

Create `tests/pFUnit/README.md` with this content:

```markdown
# pFUnit 4.15 (vendored — minimal state for the rescue)

This directory is a local checkout + build of the Goddard pFUnit 4.15 release. It is consumed by `tests/unit/meson.build` via the `PFUNIT_ROOT` environment variable (set by pixi) pointing at `tests/pFUnit/build/install_gfortran/PFUNIT-4.15`.

## Why it lives here

During the rescue (Phase 1), the install is kept in-tree so unit tests build without re-fetching pFUnit on every fresh clone. Full vendoring as a meson subproject — matching how `test-drive`, `toml-f`, and `ttutil` are handled — is a deferred item (Phase 1.5 or Phase 2).

## If the install is ever lost

Either:

1. Delete `build/` and rebuild from the local clone:

       cd tests/pFUnit
       mkdir build && cd build
       cmake -DCMAKE_INSTALL_PREFIX=install_gfortran ..
       make -j
       make install

2. Or replace this directory entirely with a fresh clone of https://github.com/Goddard-Fortran-Ecosystem/pFUnit at the `v4.15` tag, then rebuild as above.

The meson config at `tests/unit/meson.build` expects the install at the path above and errors with a pointer back here if it is missing.

## Why this is not a submodule

Historical: it was dropped in as a raw clone before the submodule convention was established elsewhere in the repo. Converting it to a proper submodule or meson subproject is tracked as deferred work.
```

### Step 3: Commit

```bash
git add tests/pFUnit/README.md
git status --short
git commit -m "docs(tests): document pFUnit local install + deferred vendoring

The rescue Phase 1 leaves tests/pFUnit/ as a local checkout + built
install. Real vendoring (meson subproject like test-drive / toml-f /
ttutil) is deferred to Phase 1.5 or 2. README captures the current
state and rebuild steps so the install can be reproduced without
ambient knowledge."
```

---

## Task 7: Clean up `.github/agents/`

**Goal:** Keep exactly one agent definition, updated to reflect this spec. Delete the drift-era stage-2 and repetition-hunter definitions.

**Files:**
- Modify: `/home/zawadzkim/Code/swap/.github/agents/swap-fortran.agent.md`
- Delete: `/home/zawadzkim/Code/swap/.github/agents/swap-fortran-stage2.agent.md`
- Delete: `/home/zawadzkim/Code/swap/.github/agents/fortran-repetition-hunter.agent.md` (if present)

### Step 1: Inspect what's there

```bash
ls /home/zawadzkim/Code/swap/.github/agents/
```

Expected: at least `swap-fortran.agent.md` plus `swap-fortran-stage2.agent.md` and possibly `fortran-repetition-hunter.agent.md`.

### Step 2: Rewrite `swap-fortran.agent.md`

Replace its content with:

```markdown
---
name: swap-fortran
description: Fortran agent for the SWAP modernization repo — operates under the rescue spec.
tools: [read, edit, search, bash, todo]
---

# swap-fortran agent

Use this agent to work on the SWAP modernization. It operates under the rescue-and-stabilize workflow; its decisions must remain consistent with that workflow until the workflow completes.

## Required context before acting

Read these before making any change:

1. `docs/superpowers/specs/2026-04-22-rescue-and-stabilize-design.md` — the rescue spec.
2. `docs/superpowers/specs/2026-04-22-baseline-record.md` — what was observed at the green baseline `e256bc0`.
3. The current Phase's plan under `docs/superpowers/plans/`.
4. `docs/adr/` — every ADR (architecture decision record). They record the non-obvious choices.

## Non-negotiables

- **Compiler**: gfortran only. Do not re-introduce Intel/ifx code paths. See `docs/adr/0001-gfortran-first.md`.
- **Build**: a single `builddir/` with `enable_pfunit=true` by default. See `docs/adr/0002-single-builddir.md`.
- **Verification**: every commit must pass `pixi run check-fast`. Phase tags also require `pixi run check-full` to be green at documented tolerances.
- **Accepted tolerances**: `oxygenstress` MOWDM deviation and `macropore` ~343 s runtime are pre-existing conditions recorded in the baseline record; do not treat them as regressions, but do not "fix" them in this rescue either.
- **Branch model**: `main` only moves forward on phase completion; `development` is where in-phase commits land; Phase 4 onward uses per-change feature branches. Nothing is pushed to `origin/main` during the rescue.

## When in doubt

Stop and ask. A partial, honest report beats a confident-but-wrong commit.
```

### Step 3: Delete the other agent files

```bash
cd /home/zawadzkim/Code/swap
test -e .github/agents/swap-fortran-stage2.agent.md && rm .github/agents/swap-fortran-stage2.agent.md
test -e .github/agents/fortran-repetition-hunter.agent.md && rm .github/agents/fortran-repetition-hunter.agent.md
ls .github/agents/
```

Expected: only `swap-fortran.agent.md` remains.

### Step 4: Commit

```bash
git add .github/agents/
git status --short
git commit -m "chore(agents): consolidate .github/agents to a single rescue-aware definition

Delete swap-fortran-stage2.agent.md and fortran-repetition-hunter.
agent.md — drift-era definitions that conflict with the rescue
workflow. Keep swap-fortran.agent.md, rewritten to reference the
rescue spec, baseline record, ADR directory, fast/full test
protocol, and branch model."
```

---

## Task 8: Write ADR 0001 (gfortran-first) and ADR 0002 (single-builddir)

**Goal:** Document the two policy decisions made in Tasks 1 and 4 as permanent architecture decision records.

**Files:**
- Create: `/home/zawadzkim/Code/swap/docs/adr/0001-gfortran-first.md`
- Create: `/home/zawadzkim/Code/swap/docs/adr/0002-single-builddir.md`

### Step 1: Ensure the adr directory exists

```bash
mkdir -p /home/zawadzkim/Code/swap/docs/adr
```

### Step 2: Write ADR 0001

Create `docs/adr/0001-gfortran-first.md` with this exact content:

```markdown
# ADR 0001 — gfortran-first (rescue phase)

Status: accepted (2026-04-22, during rescue Phase 1)

## Context

At the rescue baseline commit `e256bc0`, `meson.build` branched between Intel Fortran (`ifx` / `ifort`) and GCC (`gfortran`), with a separate prebuilt-binary download path for ttutil under Intel and a source-subproject path under GCC. The baseline verification (Phase 0) inadvertently used `ifx` because the author's shell had sourced the Intel oneAPI environment; the drift-era survey had assumed gfortran-only. Two supported compilers means two sets of flags, two build behaviours, two failure modes to investigate when a regression surfaces.

The rescue's explicit ordering is: *working car first, then fine-tuned bolid*. Reducing compiler variance removes one axis of uncertainty while physics and test-coverage debt are being paid down.

## Decision

During the rescue (Phases 1–4), the only officially supported Fortran compiler is **GCC / gfortran** as provided by pixi. `meson.build` rejects any other compiler at configure time with an explicit error message.

## Consequences

- **Positive**: deterministic builds. One flag set. ttutil always resolves via the meson subproject. No Intel-specific download-retry / corporate-proxy path.
- **Positive**: simpler meson.build — the Intel branch (~30 lines including the ttutil download logic) is gone.
- **Negative**: the Intel binary has historically been faster at some physics hot paths. A performance regression follow-on spec (between Phase 4 and the compartment refactor) will re-evaluate.
- **Negative**: users who habitually source Intel oneAPI in their shell will see a clear configure error. The error message points at this ADR.

## Re-enablement

Intel/ifx support may return after Phase 4 exit, once the codebase is stabilised. The criteria for re-enablement:
1. All pre-compartment gates from the rescue spec satisfied.
2. A dedicated follow-on spec (Phase 4.x or separate) scoped to restore Intel with regression green on both compilers.
3. The regression harness extended to exercise both compilers in CI.

Reverting this ADR is a deliberate, test-gated act — not an environment-variable flip.
```

### Step 3: Write ADR 0002

Create `docs/adr/0002-single-builddir.md` with this content:

```markdown
# ADR 0002 — Single `builddir/` with pFUnit on by default

Status: accepted (2026-04-22, during rescue Phase 1)

## Context

Phase 0 inherited two separate meson build directories:

- `builddir/` — production: `enable_pfunit=false`. Produces the `swap` binary.
- `builddir_gfortran/` — tests: `enable_pfunit=true`. Produces the pFUnit test binary.

Separate directories mean every compilation unit was built twice, every `meson setup` ran twice, and an operator who built with `pixi run build-linux` silently didn't build tests. Two pixi configure tasks (`_configure-swap`, `_configure-pfunit`) encoded this split in the task graph.

## Decision

A single `builddir/` used for everything. `meson_options.enable_pfunit` defaults to `true`. `pixi.toml` has one `_configure` task, one `build-linux` (alias for `meson compile -C builddir`), and one `test-pfunit` (`meson test -C builddir --suite unit-pfunit`).

## Consequences

- **Positive**: half the compile time for any change that isn't test-only.
- **Positive**: one obvious way to do each operation. New contributors don't need to know about the split.
- **Positive**: pFUnit tests build by default, so a configure regression in test land is caught immediately instead of when someone remembers to run the test-specific task.
- **Negative**: pFUnit becomes a required dependency at configure time — if the install at `tests/pFUnit/build/install_gfortran/PFUNIT-4.15` is missing, configure fails unless `PFUNIT_ROOT` is set or `enable_pfunit` is explicitly overridden to `false`. This is acceptable for the rescue since `tests/pFUnit/` is vendored in-tree (see `tests/pFUnit/README.md`).

## Override for exotic cases

Pass `-Denable_pfunit=false` to `meson setup` to build the production binary only. This is intended for CI jobs that specifically want to exercise the non-test path, not for everyday use.
```

### Step 4: Commit

```bash
git add docs/adr/0001-gfortran-first.md docs/adr/0002-single-builddir.md
git status --short
git commit -m "docs(adr): ADR 0001 (gfortran-first) and ADR 0002 (single builddir)

Lock in the two policy decisions made during Phase 1:

0001: gfortran is the only supported compiler until Phase 4 exit.
0002: one builddir, pFUnit on by default, unified pixi task graph.

Each ADR records context, decision, consequences, and an explicit
re-enablement / override path so future work knows how to undo
the decision cleanly if that ever becomes the right move."
```

---

## Task 9: Run `pixi run check-full`, confirm green, and tag `rescue/phase-1-infra`

**Goal:** Close Phase 1 with a clean verification and a tagged milestone. Fast-forward `main` to `development`.

**Files:** none modified.

### Step 1: Clean build from scratch and run full verification

```bash
cd /home/zawadzkim/Code/swap
rm -rf builddir
time pixi run check-full 2>&1 | tee /tmp/phase1-check-full.log | tail -40
```

Expected: full verification pipeline completes. Build green. pFUnit runs (some individual test pass/fail is acceptable). Regression: 5/6 pass with the known MOWDM deviation only.

If `pixi run check-full` fails in a way that wasn't present at Phase 0 baseline, STOP — investigate the regression before proceeding to the tag. The compiler change in Task 1 is the most likely source of any new failure; revert or adjust.

### Step 2: Record the Phase 1 check-full log

```bash
cp /tmp/phase1-check-full.log /home/zawadzkim/Code/swap/tests/regression/baselines/phase-1-infra.log
git add tests/regression/baselines/phase-1-infra.log
git commit -m "docs(baselines): record phase 1 check-full output"
```

### Step 3: Fast-forward main

```bash
cd /home/zawadzkim/Code/swap
git log --oneline -1 main
git log --oneline -1 development
git branch -f main development
git log --oneline -1 main
git log --oneline -1 origin/main
```

Expected: local `main` moves to the tip of `development`. `origin/main` unchanged at `7587ca3`.

### Step 4: Tag `rescue/phase-1-infra`

```bash
git tag -a rescue/phase-1-infra main -m "$(cat <<'EOF'
Phase 1 complete: repo & infrastructure hygiene.

- meson.build locked to gfortran (ADR 0001); Intel/ifx paths removed.
- Single builddir/; enable_pfunit on by default (ADR 0002).
- tests/unit/meson.build reconciled with baseline sources; orphan
  test .f90 files deleted; testSuites.inc regenerated.
- check-fast and check-full pixi tasks added as canonical gates.
- tests/pFUnit/ local install documented; full vendoring deferred.
- .github/agents/ consolidated to swap-fortran.agent.md only.

check-full green at documented tolerances
(oxygenstress MOWDM and macropore ~343s accepted, unchanged from
Phase 0 baseline).

Next: Phase 2 (canonical documentation). ADR 0005 will address the
root LICENSE LGPL v2.1 vs GPL v2 finding carried from Phase 0.
EOF
)"
git tag -l 'rescue/*'
git rev-parse rescue/phase-1-infra^{commit}
git rev-parse main
```

Expected:
- Tag exists.
- `rescue/phase-1-infra^{commit}` and `main` resolve to the same SHA.

### Step 5: Final state summary

```bash
echo "=== branches ==="
git branch
git branch --list 'archive/*'
git branch --list 'legacy/*'
echo "=== tags ==="
git tag -l 'rescue/*'
echo "=== HEAD / main / origin ==="
git log --oneline -1 main
git log --oneline -1 development
git log --oneline -1 origin/main
```

Expected:
- Local branches: `main`, `development` (current), `legacy/swap-4.2.0`, four `archive/*`.
- Tags: `rescue/phase-0-baseline`, `rescue/phase-1-infra`.
- `main = development`. `origin/main = 7587ca3` (untouched).

---

## Plan self-review

**Spec coverage.** Phase 1 of the spec lists eight issue / action rows. Mapping:

| Spec row | Plan task |
|---|---|
| Two meson build directories → unify | Task 4 |
| `fpm.toml` ambiguous role | N/A at baseline (drift-era only; Phase 0 already confirmed absent) |
| pFUnit built locally → vendor or pin | Task 6 (minimal — leave install, document; full vendoring deferred) |
| `tests/swap-cases` submodule at detached HEAD | Deferred (already declares `branch=main`; Phase 2 handles the dirtying-on-regression doc) |
| Hardcoded reference binary path | Already resolved in Phase 0 tail commit (`3a65acf`) |
| Missing `check-fast` / `check-full` | Task 5 |
| Optional pre-commit hook | Deferred (spec explicitly marks optional) |
| Tests mixing `.pf` and standalone `.f90` | Task 2 step 3 (deletes the orphan `.f90` files; pFUnit-only survivors) |
| `testSuites.inc` drift | Task 3 (manual regen; automated check deferred) |

Plus the three Phase 0 findings that the spec said Phase 1 would address:
- pFUnit broken at baseline → Task 2 (primary).
- Compiler policy → Task 1 + Task 8 ADR 0001.
- `.github/agents/` drift-era files → Task 7.

Deferred items are each recorded at the top of the plan under "Deferred out of Phase 1".

**Placeholder scan.** No "TBD", "TODO", "implement later", "similar to Task N", or description-without-code steps. Every command, every code block, every edit has its full content.

One place that depends on runtime behavior: Task 2 Step 4 and Step 5 include "iterate until configure/compile succeeds" because the exact residual set of missing modules can't be predicted without running meson against the baseline. This is labeled explicitly — it's not a placeholder, it's an honest instruction to iterate within a bounded loop.

**Type consistency.** Task names used across the plan:
- `builddir` as the single canonical build directory (used in Tasks 4, 5, 9). Never `builddir_gfortran` in any task after Task 4.
- pixi task names: `_fetch-deps`, `_configure`, `build-linux`, `test-pfunit`, `check-fast`, `check-full`, `regression`. Consistent from Task 4 onward.
- `enable_pfunit` option name. Consistent everywhere.
- `tests/regression/baselines/phase-N-*.log` naming. Consistent with Phase 0's `phase-0-baseline.log`.
- `rescue/phase-N-name` tag naming. Consistent.

**Scope discipline.** Phase 1 produces a working, testable baseline (Task 9 Step 1 proves this). Every deferred item has an explicit phase target (1.5, 2, 3, or 4). No task is open-ended.
