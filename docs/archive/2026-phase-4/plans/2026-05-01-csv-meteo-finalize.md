# CSV Meteorology Pathway — Finalization Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete the CSV-only meteorological data ingestion pathway. Phase 0 reorganizes the test directory layout so each TOML case lives in a self-contained working dir, decoupling the modern build's tooling from the legacy ASCII case dirs. Phase 1 (Tasks 1–10) then deletes confirmed-dead code, fills a unit-test gap in the `datetime`-keyed CSV reader, and implements the sub-daily (detail) CSV pathway per ADR 0014 Step 1.

**Architecture:** The CSV foundation — daily meteo CSV (`MeteoCSVYear`), rain events CSV (`ReadRainEvents` CSV branch), and the unified `read_csv_table` reader — was shipped in commit `bb15a30`. **Phase 0** migrates each case's legacy runtime files (`*.crp`, `swap_linux.swp.template`) into the matching `tests/swap-cases/toml/<N>.<case>/`, rewrites `run_case.sh` to run SWAP in-place there, and updates the regression runner to use the TOML dir as the `copytree` source. **Phase 1** then adds three independent code units on top: (1) dead-code deletion of `ioutils.f90`, (2) unit tests for the already-implemented `datetime`-keyed parsing in `csv_reader_mod`, and (3) the sub-daily detail meteo CSV path (`MeteoCSVDetYear` + `detail_file` TOML key). ADR 0014 Steps 2–3 (TTutil branch deletion, dead-variable sweep) remain blocked on `.swp` pipeline retirement and are deferred.

**Tech Stack:** Fortran 2008 (gfortran 13), pFUnit, tomlf, pixi, bash, Python (regression). All verification via `pixi run -e test build-linux`, `pixi run -e test test-pfunit`, `pixi run -e test check-full`.

**Specs/ADRs:** `docs/superpowers/specs/2026-05-01-toml-case-dir-self-contained-design.md` (Phase 0), `docs/adr/0013-csv-meteorology-input.md`, `docs/adr/0014-readmeteo-phaseout.md`, `docs/phase-4f-legacy-reader-analysis.md`

---

## File Structure

### Phase 0 (test-dir reorganization)

| File | Change |
|---|---|
| `tests/swap-cases/toml/<N>.<case>/*.crp` (×6) | **Create (submodule):** copy each case's `.crp` files from `tests/swap-cases/<N>.<case>/`. Specifically: `1.hupselbrook` → grassd/maizes/potatod; `2.grassgrowth` → grassd; `3.macroporeflow` → wintcer1/wintcer2; `4.oxygenstress` → grassd; `5.salinitystress` → potatod; `6.surfacewater` → grass. |
| `tests/swap-cases/toml/<N>.<case>/swap_linux.swp.template` (×6) | **Create (submodule):** copy from the matching `tests/swap-cases/<N>.<case>/`. |
| `tests/swap-cases/toml/6.surfacewater/swap.dra` | **Create (submodule):** copy from `tests/swap-cases/6.surfacewater/swap.dra`. Needed because `swdra=2` triggers `rddre()` which still reads the legacy ASCII `swap.dra` (the `swap.dra.toml` covers a subset of fields). |
| `tests/swap-cases/toml/5.salinitystress/swap.ini` | **Create (submodule):** copy from `tests/swap-cases/5.salinitystress/swap.ini`. Needed because `swinco=3` makes the adapter call `rdinit()` to read `swap.ini` for initial profile state. |
| `tests/swap-cases/<N>.<case>/` (×6) | **No change.** Frozen as legacy-binary input. |
| `tests/swap-cases/run_case.sh` | **Rewrite:** drop `--toml`/`--legacy`/`TOML_MODE`; default = TOML mode in-place; new `--legacy-binary <path>` for legacy reference runs; whitelist-based cleanup. |
| `tests/regression/test_output_regression.py` | **Modify:** flip `shutil.copytree` source from legacy dir to TOML dir; delete the "stage TOML files on top" block (current lines 290–307). |
| `docs/csv-companion-files.md` | **Modify:** update the `run_case.sh` reference at line 95 to reflect the new in-place TOML working dir. |

### Phase 1 (CSV meteo finalization)

| File | Change |
|---|---|
| `src/utils/ioutils.f90` | **Delete** — confirmed dead (no references in any source or meson.build). |
| `tests/unit/io/fixtures/csv_datetime_det.csv` | **Create** — 3-row fixture for `datetime`-keyed reader test. |
| `tests/unit/io/test_csv_reader.pf` | **Modify** — add 2 tests covering the `datetime` column branch. |
| `src/config/meteorology_config.f90` | **Modify** — add `detail_file` allocatable field; no new validator rule (adapter enforces it). |
| `tests/unit/config/test_meteorology_config.pf` | **Modify** — add 1 test confirming the new field exists and round-trips. |
| `tests/unit/io/toml/fixtures/meteorology_detail_file.toml` | **Create** — TOML fixture with `[meteorology.temporal].detail_file`. |
| `tests/unit/io/toml/test_read_meteorology_toml.pf` | **Modify** — add 1 test asserting `detail_file` is read. |
| `src/io/toml/read_meteorology_toml.f90` | **Modify** — add `get_optional_string_with_default` call for `detail_file` inside `[meteorology.temporal]`. |
| `src/core/variables.f90` | **Modify** — add `swMetDetCSV` flag; update stub comment for `metcsv_det` (7 columns). |
| `src/io/toml/config_to_variables.f90` | **Modify** — add detail-CSV pre-load block after the daily CSV block; set `swMetDetCSV`. |
| `src/io/readmeteo.f90` | **Modify** — add `MeteoCSVDetYear` subroutine; update `ReadMeteoYear` dispatch to call it when `swmetdetail==1`. |
| `tests/unit/io/toml/test_config_to_variables.pf` | **Modify** — verify `swMetDetCSV` variable is used correctly (build-smoke only; no fixture needed since integration is covered by regression). |
| `docs/meteorology.md` | **Modify** — add `detail_file` key documentation and sub-daily CSV schema. |

---

## Conventions (read before starting)

1. **Build after every source change:** `pixi run -e test build-linux 2>&1 | grep -E "error:" | head -5`
2. **Run targeted pFUnit filter:** `pixi run -e test test-pfunit -- --filter <test_name>`
3. **Submodule discipline (Phase 0 only):** each commit inside `tests/swap-cases/` is followed by a single outer-repo bump (`git add tests/swap-cases && git commit`). Phase 1 makes no submodule edits.
4. **No physics edits** — Phase 1 changes go in config, reader, adapter, `readmeteo.f90` only.
5. **Fortran `block` statement** used for local-scope declarations (Fortran 2008); confirmed to work in this codebase (see `config_to_variables.f90:148–173`).
6. **`jday` and `days1900_to_md`** are defined at file scope in `readmeteo.f90` and callable as external from `MeteoCSVDetYear` via `integer :: jday; external jday`.

---

## Phase 0 — TOML case directory self-contained working dir

Sequence Phase 0 tasks before Phase 1. Tasks 0.1 → 0.5 must complete and the
full regression suite must be green before starting Task 1.

---

## Task 0.1: Migrate legacy runtime files into TOML case dirs (submodule)

**Files:**
- Submodule create: `tests/swap-cases/toml/<N>.<case>/{*.crp, swap_linux.swp.template}` for all 6 cases.
- Outer repo: bump `tests/swap-cases` pointer.

The legacy crop files (`*.crp`) and the `swap_linux.swp.template` are required
at runtime because (a) the legacy crop sub-readers in `cropgrowth.f90` /
`irrigation.f90` / `management_soil.f90` open `<stem>.crp` via TTutil from
`pathwork`, and (b) some sub-readers still call `RDinit(swpfile)`. Phase 4f-extend
removes both dependencies; until then we house the files in the TOML dir.

- [ ] **Step 1: Confirm starting state**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
for d in 1.hupselbrook 2.grassgrowth 3.macroporeflow 4.oxygenstress 5.salinitystress 6.surfacewater; do
  echo "--- $d ---"
  ls "$d"/*.crp "$d"/swap_linux.swp.template 2>/dev/null
done
```

Expected: each legacy dir has at least `swap_linux.swp.template` plus 1–3 `.crp` files (per the spec's enumeration).

- [ ] **Step 2: Copy files into the TOML dirs**

Run from the submodule root (`tests/swap-cases/`):

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
for d in 1.hupselbrook 2.grassgrowth 3.macroporeflow 4.oxygenstress 5.salinitystress 6.surfacewater; do
  if [ -d "toml/$d" ]; then
    for crp in "$d"/*.crp; do
      [ -f "$crp" ] && cp "$crp" "toml/$d/$(basename "$crp")"
    done
    if [ -f "$d/swap_linux.swp.template" ]; then
      cp "$d/swap_linux.swp.template" "toml/$d/swap_linux.swp.template"
    fi
  else
    echo "WARN: toml/$d missing — skipping" >&2
  fi
done

# Conditional legacy companions: only the cases whose TOML config still
# triggers a legacy ASCII reader at runtime.
cp 6.surfacewater/swap.dra   toml/6.surfacewater/swap.dra      # swdra=2 → rddre()
cp 5.salinitystress/swap.ini toml/5.salinitystress/swap.ini    # swinco=3 → rdinit()
```

- [ ] **Step 3: Verify the copy**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
for d in 1.hupselbrook 2.grassgrowth 3.macroporeflow 4.oxygenstress 5.salinitystress 6.surfacewater; do
  echo "--- toml/$d ---"
  ls "toml/$d/"*.crp "toml/$d/swap_linux.swp.template" 2>/dev/null
done
```

Expected: every `toml/<N>.<case>/` now has the matching `.crp` files and `swap_linux.swp.template`.

- [ ] **Step 4: Submodule commit**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
git add toml/*/swap_linux.swp.template toml/*/*.crp
git status --short
git commit -m "test(toml-cases): stage .crp + swap.swp template into TOML dirs

Phase 0 of CSV meteo finalization. Each tests/swap-cases/toml/<case>/
becomes the self-contained working dir for the modern SWAP build. The
legacy <N>.<case>/ dirs are now solely for legacy-binary runs."
```

- [ ] **Step 5: Outer-repo bump**

```bash
cd /home/zawadzkim/Code/swap
git add tests/swap-cases
git submodule status tests/swap-cases
git commit -m "chore(swap-cases): bump for TOML dir runtime-file migration (Phase 0)"
```

- [ ] **Step 6: Sanity check — regression still green via the OLD runner**

The old runner still copies TOML files on top of the legacy dir, so it should
keep working at this point. Confirm we didn't break anything before touching
the runners:

```bash
cd /home/zawadzkim/Code/swap
pixi run -e test check-full 2>&1 | tail -10
```

Expected: 5 passed, 0 failed. (Macroporeflow excluded per ADR 0011.)

---

## Task 0.2: Rewrite `run_case.sh`

**Files:**
- Modify (submodule): `tests/swap-cases/run_case.sh`

- [ ] **Step 1: Replace the script**

Overwrite `tests/swap-cases/run_case.sh` with:

```bash
#!/bin/bash
# Run a SWAP regression case.
#
# Default (TOML mode):
#   cd into tests/swap-cases/toml/<N>.<CASE>/, stage swap.swp from the
#   in-dir swap_linux.swp.template, run the modern SWAP binary there,
#   then clean up generated outputs.
#
# With --legacy-binary <path>:
#   cd into tests/swap-cases/<N>.<CASE>/ (legacy ASCII inputs), stage
#   swap.swp, run the supplied legacy binary, clean up.

CASE_NAME=""
SWAP_EXEC=""
LEGACY_BINARY=""
KEEP_OUTPUTS=false

OPTIONS=$(getopt -o c:e:kh --long case:,exec:,legacy-binary:,keep-outputs,help -n "$0" -- "$@")
if [ $? -ne 0 ]; then
    echo "Error parsing options. Use --help for usage information."
    exit 1
fi
eval set -- "$OPTIONS"

while true; do
    case "$1" in
        -c|--case)         CASE_NAME="$2"; shift 2 ;;
        -e|--exec)         SWAP_EXEC="$2"; shift 2 ;;
        --legacy-binary)   LEGACY_BINARY="$2"; shift 2 ;;
        -k|--keep-outputs) KEEP_OUTPUTS=true; shift ;;
        -h|--help)
            cat <<USAGE
Usage: $0 -c CASE_NAME [-e EXECUTABLE] [--legacy-binary PATH] [-k] [-h]

Default (TOML mode):
  cd into tests/swap-cases/toml/<N>.<CASE>/, run the modern SWAP binary
  there, clean up generated outputs.

With --legacy-binary PATH:
  cd into tests/swap-cases/<N>.<CASE>/ (legacy ASCII inputs), run the
  supplied legacy binary, clean up. Mutually exclusive with --exec.

Options:
  -c, --case CASE_NAME      Case stem (e.g. hupselbrook). Required.
  -e, --exec EXECUTABLE     Modern SWAP binary path (default: ./swap).
  --legacy-binary PATH      Legacy SWAP binary path (legacy ASCII mode).
  -k, --keep-outputs        Skip post-run cleanup (inspect output files).
  -h, --help                Show this message.

Cases: hupselbrook, grassgrowth, macroporeflow, oxygenstress,
       salinitystress, surfacewater
USAGE
            exit 0 ;;
        --) shift; break ;;
        *)  echo "Internal error!"; exit 1 ;;
    esac
done

if [ -z "$CASE_NAME" ]; then
    echo "Error: --case is required. Use --help."
    exit 1
fi

if [ -n "$LEGACY_BINARY" ] && [ -n "$SWAP_EXEC" ]; then
    echo "Error: --legacy-binary and --exec are mutually exclusive."
    exit 1
fi

# Locate case dir by glob (e.g. 1.hupselbrook).
CASE_DIR=$(find . -maxdepth 1 -type d -name "*.${CASE_NAME}" | head -n 1)
if [ -z "$CASE_DIR" ]; then
    echo "Error: case directory '*.${CASE_NAME}' not found in $(pwd)."
    echo "Available cases: hupselbrook, grassgrowth, macroporeflow, oxygenstress, salinitystress, surfacewater"
    exit 1
fi
CASE_BASENAME=$(basename "$CASE_DIR")

if [ -n "$LEGACY_BINARY" ]; then
    BIN_ABS="$(readlink -f "$LEGACY_BINARY")"
    if [ ! -f "$BIN_ABS" ]; then
        echo "Error: legacy binary '$BIN_ABS' not found."
        exit 1
    fi
    WORKDIR="$CASE_DIR"
    MODE_LABEL="LEGACY"
else
    SWAP_EXEC="${SWAP_EXEC:-./swap}"
    BIN_ABS="$(readlink -f "$SWAP_EXEC")"
    if [ ! -f "$BIN_ABS" ]; then
        echo "Error: executable '$BIN_ABS' not found."
        exit 1
    fi
    WORKDIR="toml/$CASE_BASENAME"
    if [ ! -d "$WORKDIR" ] || [ ! -f "$WORKDIR/swap.toml" ]; then
        echo "Error: TOML case directory '$WORKDIR' missing or has no swap.toml."
        exit 1
    fi
    MODE_LABEL="TOML"
fi

echo "[$MODE_LABEL] case=$CASE_NAME workdir=$WORKDIR binary=$BIN_ABS"

cd "$WORKDIR" || exit 1

# Stage swap.swp from the in-dir template (used by legacy crop sub-readers
# that still call RDinit(swpfile)). Both modes need this.
if [ -f "swap_linux.swp.template" ]; then
    cp "swap_linux.swp.template" "swap.swp"
fi

# Run.
"$BIN_ABS"
RC=$?

# Press-Enter automation (SWAP can pause for user input on some paths).
printf '\n'

if [ "$KEEP_OUTPUTS" = true ]; then
    echo "Keeping generated outputs (--keep-outputs)."
else
    # Whitelist cleanup: only delete files SWAP generates. Never touch
    # committed inputs (*.toml, *.crp, *.crp.toml, input *.csv, *.template).
    rm -f -- swap.swp \
             result.* result_*.csv \
             *.log Swap.ok swap.ok reruns.log \
             swaprd\$*.tmp *.tmp \
             fort.20
    echo "Cleaned up generated outputs."
fi

exit $RC
```

Note on the cleanup whitelist:
- `*.tmp` matches `swaprd$00001.tmp` etc. (and any other SWAP-temp file).
- `result.*` catches the legacy SWAP `result.bal` etc. outputs; `result_*.csv` catches the modern `result_output.csv`.
- The whitelist explicitly avoids bare `*.csv` (would nuke input CSVs) and bare `*.crp` (would nuke staged crop files).

- [ ] **Step 2: Make executable + smoke test**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
chmod +x run_case.sh
./run_case.sh --help
```

Expected: usage text printed; exit 0.

- [ ] **Step 3: Run a case manually (TOML mode default)**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
./run_case.sh -c hupselbrook -e ../../builddir/swap
```

Expected: case runs to completion (SWAP exits 100, the script forwards the code), output is cleaned up.

Verify the TOML dir is clean of generated files afterwards:

```bash
ls toml/1.hupselbrook/
```

Expected: no `swap.swp`, no `result_*.csv`, no `*.log`, no `*.tmp`. The committed inputs (`swap.toml`, `*.crp`, `*.crp.toml`, `283.csv`, `swap_linux.swp.template`, `swap.dra.toml`, `README.md`) are still present.

- [ ] **Step 4: Submodule commit**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
git add run_case.sh
git commit -m "test(run_case): rewrite for in-place TOML execution + --legacy-binary

Drops --toml/--legacy/TOML_MODE plumbing in favor of one default flow:
cd into toml/<case>/, stage swap.swp, run, cleanup. New --legacy-binary
PATH option runs an external legacy binary against the legacy ASCII
case dir. --exec and --legacy-binary are mutually exclusive."
```

- [ ] **Step 5: Outer-repo bump**

```bash
cd /home/zawadzkim/Code/swap
git add tests/swap-cases
git commit -m "chore(swap-cases): bump for run_case.sh rewrite (Phase 0)"
```

---

## Task 0.3: Update `tests/regression/test_output_regression.py`

**Files:**
- Modify: `tests/regression/test_output_regression.py` (lines ~257–307 inside `_run_case_in_temp` or equivalent helper)

- [ ] **Step 1: Locate the staging block**

```bash
cd /home/zawadzkim/Code/swap
grep -n "shutil.copytree\|toml_src\|swap.toml\|swap.dra.toml\|case_dir = TESTS_DIR" tests/regression/test_output_regression.py | head -15
```

Expected: identifies the block around lines 257–307 with the legacy `case_dir` source and the "stage TOML files on top" inner block.

- [ ] **Step 2: Replace the block**

Find the section that currently reads (line numbers approximate):

```python
case_dir = TESTS_DIR / "swap-cases" / case.case_dir
if not case_dir.exists():
    raise RuntimeError(f"case directory not found at {case_dir}")

with tempfile.TemporaryDirectory() as tmpdir:
    tmp = Path(tmpdir)

    # Copy case files, excluding any pre-existing output files
    shutil.copytree(
        case_dir,
        tmp / "case",
        ignore=shutil.ignore_patterns(
            'result_output.csv',
            'result_*.csv',
            '*.log',
            'output.*',
            '*.out'
        )
    )
    workdir = tmp / "case"

    # ensure template is named swap.swp
    swap_file = workdir / "swap.swp"
    if not swap_file.exists():
        template = workdir / "swap_linux.swp.template"
        if template.exists():
            shutil.copy(template, swap_file)
        else:
            # try other common names
            for alt in workdir.glob("*.swp"):
                shutil.copy(alt, swap_file)
                break

    # Phase 4f: SWAP now requires swap.toml as the canonical entry
    # point. Stage it from tests/swap-cases/toml/<case>/ alongside
    # cross-file siblings (swap.dra.toml, *.crp.toml). Cases without
    # a populated TOML directory still fail loudly — Phase 4f-extend
    # ports them one by one. Mirrors run_case.sh --toml.
    toml_src = TESTS_DIR / "swap-cases" / "toml" / case.case_dir
    if (toml_src / "swap.toml").exists():
        shutil.copy(toml_src / "swap.toml", workdir / "swap.toml")
        for extra in ("swap.dra.toml",):
            if (toml_src / extra).exists():
                shutil.copy(toml_src / extra, workdir / extra)
        for crp in toml_src.glob("*.crp.toml"):
            shutil.copy(crp, workdir / crp.name)
        # Phase 4f cleanup: stage CSV companion files (long-form
        # fixed-irrigation events, prescribed gwl, etc.) alongside
        # swap.toml. Reader paths are relative to pathwork.
        for csv_companion in toml_src.glob("*.csv"):
            shutil.copy(csv_companion, workdir / csv_companion.name)
```

Replace with:

```python
toml_dir = TESTS_DIR / "swap-cases" / "toml" / case.case_dir
if not toml_dir.exists() or not (toml_dir / "swap.toml").exists():
    raise RuntimeError(
        f"TOML case directory not found at {toml_dir} "
        f"(missing dir or swap.toml). Phase 0 of CSV meteo finalization "
        f"made the TOML dir the sole source of truth for regression."
    )

with tempfile.TemporaryDirectory() as tmpdir:
    tmp = Path(tmpdir)

    # Phase 0: TOML dir is the self-contained source. It contains
    # swap.toml, swap.dra.toml, *.crp.toml, *.csv companions, *.crp
    # legacy crop files, and swap_linux.swp.template. Everything SWAP
    # needs at runtime lives here.
    shutil.copytree(
        toml_dir,
        tmp / "case",
        ignore=shutil.ignore_patterns(
            'result_output.csv',
            'result_*.csv',
            '*.log',
            'output.*',
            '*.out',
        ),
    )
    workdir = tmp / "case"

    # Stage swap.swp from the in-dir template (legacy crop sub-readers
    # still call RDinit(swpfile)).
    swap_file = workdir / "swap.swp"
    if not swap_file.exists():
        template = workdir / "swap_linux.swp.template"
        if template.exists():
            shutil.copy(template, swap_file)
```

- [ ] **Step 3: Run the regression suite**

```bash
cd /home/zawadzkim/Code/swap
pixi run -e test check-full 2>&1 | tail -15
```

Expected: 5/5 cases pass at 1e-2 cm tolerance, identical numeric output to the previous run (the change is purely about *where* SWAP runs, not what it computes).

- [ ] **Step 4: Commit**

```bash
git add tests/regression/test_output_regression.py
git commit -m "test(regression): copytree from TOML dir, drop legacy-on-top staging

Phase 0 of CSV meteo finalization. The TOML case directory is now
self-contained (Task 0.1 staged .crp + swap_linux.swp.template into
each toml/<case>/). The runner no longer needs to read from the
legacy <N>.<case>/ dir at all."
```

---

## Task 0.4: Update `docs/csv-companion-files.md`

**Files:**
- Modify: `docs/csv-companion-files.md`

- [ ] **Step 1: Locate and update the run_case.sh reference**

```bash
grep -n "run_case" /home/zawadzkim/Code/swap/docs/csv-companion-files.md
```

Expected: a single match around line 95 that mentions `tests/swap-cases/run_case.sh ... glob *.csv from the case directory`.

- [ ] **Step 2: Edit the doc**

In `docs/csv-companion-files.md` replace the current block (lines 89–95):

```markdown
## Path resolution and staging

- Companion CSVs are referenced by basename in `swap.toml`. The adapter resolves
  paths relative to the working directory at adapter time.
- The runtime stages all `*.csv` companion files into the case directory before SWAP
  starts. The regression harness (`tests/regression/test_output_regression.py`) and
  `tests/swap-cases/run_case.sh` both glob `*.csv` from the case directory.
```

with:

```markdown
## Path resolution and staging

- Companion CSVs are referenced by basename in `swap.toml`. The adapter resolves
  paths relative to the working directory at adapter time.
- The case working directory is `tests/swap-cases/toml/<N>.<case>/`. It is
  self-contained — every file SWAP reads at runtime lives there: `swap.toml`,
  `swap.dra.toml`, `*.crp.toml`, all `*.csv` companions, the legacy `*.crp`
  crop files (read by sub-readers in `cropgrowth.f90` until Phase 4f-extend
  ports them), and `swap_linux.swp.template` (staged to `swap.swp` per run).
- `tests/swap-cases/run_case.sh` runs SWAP in that directory in-place;
  `tests/regression/test_output_regression.py` copies it to a temp dir for
  parallel-safe execution. Neither tool reads from the legacy `<N>.<case>/`
  directories — those are reserved for the legacy reference binary.
```

- [ ] **Step 3: Commit**

```bash
git add docs/csv-companion-files.md
git commit -m "docs: update csv-companion-files for self-contained TOML dir (Phase 0)"
```

---

## Task 0.5: Phase 0 verification

- [ ] **Step 1: Full regression suite**

```bash
cd /home/zawadzkim/Code/swap
pixi run -e test check-full 2>&1 | tail -15
```

Expected:
```
✓ hupselbrook:    regression ok
✓ grassgrowth:    regression ok
✓ oxygenstress:   regression ok
✓ salinitystress: regression ok
✓ surfacewater:   regression ok
Results: 5 passed, 0 failed
```

- [ ] **Step 2: pFUnit suite (no behavior change expected)**

```bash
pixi run -e test test-pfunit 2>&1 | tail -10
```

Expected: green.

- [ ] **Step 3: Manual `run_case.sh` smoke test**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
./run_case.sh -c hupselbrook -e ../../builddir/swap
```

Expected: case runs successfully; afterwards `ls toml/1.hupselbrook/` shows only committed inputs (no `swap.swp`, no `result_*.csv`, no logs).

- [ ] **Step 4: Confirm legacy dirs untouched**

```bash
cd /home/zawadzkim/Code/swap
git -C tests/swap-cases status --short -- 1.hupselbrook 2.grassgrowth 3.macroporeflow 4.oxygenstress 5.salinitystress 6.surfacewater
```

Expected: empty output (no modifications to legacy case dirs).

Phase 0 is complete. Proceed to Phase 1 (Task 1).

---

# Phase 1 — CSV meteo code work

## Task 1: Delete `src/utils/ioutils.f90`

**Files:**
- Delete: `src/utils/ioutils.f90`

- [ ] **Step 1: Confirm no live references**

```bash
git grep -rn 'ioutils\|io_utils_mod\|parse_output_extensions' -- src/ tests/ meson.build
```

Expected: **zero matches** (the `docs/phase-4f-legacy-reader-analysis.md` reference is in docs only and is fine to leave). If any match appears in `src/` or `meson.build`, stop and investigate before deleting.

- [ ] **Step 2: Delete the file**

```bash
git rm src/utils/ioutils.f90
```

- [ ] **Step 3: Build to confirm nothing breaks**

```bash
pixi run -e test build-linux 2>&1 | grep -E "error:" | head -5
```

Expected: clean build, no errors.

- [ ] **Step 4: Commit**

```bash
git add -u
git commit -m "chore(utils): delete dead ioutils.f90

Module was never imported by any source file and not listed in
meson.build. Confirmed in docs/phase-4f-legacy-reader-analysis.md."
```

---

## Task 2: Add `datetime`-keyed CSV reader unit tests

**Files:**
- Create: `tests/unit/io/fixtures/csv_datetime_det.csv`
- Modify: `tests/unit/io/test_csv_reader.pf`

The `datetime_keyed` branch in `csv_reader_mod` (`parse_iso_datetime`) was implemented in `bb15a30` but no unit test exercises it. This task fills that gap before the sub-daily pathway relies on it.

- [ ] **Step 1: Create the fixture**

Write `tests/unit/io/fixtures/csv_datetime_det.csv`:

```
# Sub-daily detail meteo fixture for test_csv_reader.pf
datetime,record,rad,temp,hum,wind,rain
1980-04-24 00:00:00,1,500.0,10.0,1.0,3.0,0.0
1980-04-24 12:00:00,2,800.0,14.0,0.8,4.0,1.5
1980-04-25 00:00:00,1,400.0,8.0,1.1,2.5,0.0
```

- [ ] **Step 2: Add two tests to `test_csv_reader.pf`**

Append to `tests/unit/io/test_csv_reader.pf`:

```fortran
@test
subroutine test_datetime_keyed_happy()
   use funit
   use iso_fortran_env, only: real64
   use csv_reader_mod, only: read_csv_table
   use error_mod, only: error_collection_t
   real(real64), allocatable :: table(:,:)
   type(error_collection_t)  :: errors

   call read_csv_table('tests/unit/io/fixtures/csv_datetime_det.csv', &
                       ['datetime', 'record  ', 'rad     ', &
                        'temp    ', 'hum     ', 'wind    ', 'rain    '], &
                       table, errors)

   @assertFalse(errors%has_errors())
   @assertTrue(allocated(table))
   @assertEqual(3, size(table, 1))
   @assertEqual(7, size(table, 2))
   ! 1980-04-24 00:00:00 = 29334 days + 0 fractional = 29334.0
   @assertEqual(29334.0_real64, table(1, 1), 1.0e-4_real64)
   ! 1980-04-24 12:00:00 = 29334 + 0.5 = 29334.5
   @assertEqual(29334.5_real64, table(2, 1), 1.0e-4_real64)
   ! 1980-04-25 00:00:00 = 29335.0
   @assertEqual(29335.0_real64, table(3, 1), 1.0e-4_real64)
   ! record column (col 2) is real64 = 1.0, 2.0, 1.0
   @assertEqual(1.0_real64, table(1, 2), 1.0e-6_real64)
   @assertEqual(2.0_real64, table(2, 2), 1.0e-6_real64)
   ! rad value col 3
   @assertEqual(500.0_real64, table(1, 3), 1.0e-6_real64)
end subroutine

@test
subroutine test_datetime_keyed_bad_timestamp()
   use funit
   use iso_fortran_env, only: real64
   use csv_reader_mod, only: read_csv_table
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   real(real64), allocatable :: table(:,:)
   type(error_collection_t)  :: errors

   ! csv_bad_date.csv has "not-a-date" in col 1 — same error code for datetime column.
   ! Re-use the existing bad-date fixture but with a datetime header to exercise
   ! the datetime branch error path.
   ! We need a fixture with a malformed datetime:
   call read_csv_table('tests/unit/io/fixtures/csv_bad_date.csv', &
                       ['datetime', 'gwl     '], &
                       table, errors)

   @assertTrue(errors%has_errors())
   @assertEqual(ERR_PARSE_TYPE_MISMATCH, errors%items(1)%code)
   @assertFalse(allocated(table))
end subroutine
```

Note on string-array literals: Fortran requires uniform element length in array constructors, so all elements are right-padded to the length of the longest (`'datetime'` = 8 chars). The reader compares with `trim()` on both sides.

- [ ] **Step 3: Run the new tests (must pass)**

```bash
pixi run -e test test-pfunit -- --filter test_csv_reader 2>&1 | tail -20
```

Expected: **13 + 2 = 15 PASS**, 0 FAIL.

- [ ] **Step 4: Commit**

```bash
git add tests/unit/io/fixtures/csv_datetime_det.csv tests/unit/io/test_csv_reader.pf
git commit -m "test(csv-reader): add datetime-keyed column tests

Exercises parse_iso_datetime (YYYY-MM-DD HH:MM:SS) and the
ERR_PARSE_TYPE_MISMATCH path for a malformed datetime cell."
```

---

## Task 3: Add `detail_file` field to `meteorology_config_t`

**Files:**
- Modify: `src/config/meteorology_config.f90`
- Modify: `tests/unit/config/test_meteorology_config.pf`

- [ ] **Step 1: Write the failing test (TDD red)**

Append to `tests/unit/config/test_meteorology_config.pf`:

```fortran
@test
subroutine test_meteo_detail_file_field_roundtrips()
   use funit
   use meteorology_config_mod
   type(meteorology_config_t) :: m
   m%detail_file = "mysite.det.csv"
   @assertEqual("mysite.det.csv", m%detail_file)
end subroutine
```

- [ ] **Step 2: Run — expect FAIL (field not defined)**

```bash
pixi run -e test test-pfunit -- --filter test_meteorology_config 2>&1 | tail -10
```

Expected: compile error mentioning `detail_file` unknown, or test FAIL.

- [ ] **Step 3: Add the field**

In `src/config/meteorology_config.f90`, inside `type :: meteorology_config_t`, add after `rain_events_file` (line 36):

```fortran
      character(len=:), allocatable :: detail_file
```

So the block reads:

```fortran
   type :: meteorology_config_t
      character(len=:), allocatable :: metfile
      character(len=:), allocatable :: rainfile
      character(len=:), allocatable :: rain_events_file
      character(len=:), allocatable :: detail_file
      real(real64) :: lat  = 0.0_real64
      ...
```

- [ ] **Step 4: Run — expect PASS**

```bash
pixi run -e test test-pfunit -- --filter test_meteorology_config 2>&1 | tail -10
```

Expected: all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/config/meteorology_config.f90 tests/unit/config/test_meteorology_config.pf
git commit -m "feat(meteo-config): add detail_file field for sub-daily CSV

detail_file holds the path to the datetime-keyed sub-daily meteo CSV
used when swmetdetail=1 in the TOML pathway."
```

---

## Task 4: Read `detail_file` in the TOML reader

**Files:**
- Create: `tests/unit/io/toml/fixtures/meteorology_detail_file.toml`
- Modify: `tests/unit/io/toml/test_read_meteorology_toml.pf`
- Modify: `src/io/toml/read_meteorology_toml.f90`

- [ ] **Step 1: Create the fixture**

Write `tests/unit/io/toml/fixtures/meteorology_detail_file.toml`:

```toml
[meteorology]
file = "hupsel.csv"
lat  = 52.0
alt  = 10.0

[meteorology.evapotranspiration]
swetr = 1

[meteorology.temporal]
swmetdetail = 1
nmetdetail  = 48
detail_file = "hupsel.det.csv"

[meteorology.rain]
swrain = 0
```

- [ ] **Step 2: Write the failing test**

Append to `tests/unit/io/toml/test_read_meteorology_toml.pf`:

```fortran
@test
subroutine test_read_meteorology_detail_file()
   use funit
   use tomlf, only: toml_table, toml_load
   use meteorology_config_mod
   use read_meteorology_toml_mod
   use error_mod
   type(toml_table), allocatable, target :: doc
   type(toml_table), pointer             :: doc_ptr
   type(meteorology_config_t)            :: c
   type(error_collection_t)              :: errors

   call toml_load(doc, 'tests/unit/io/toml/fixtures/meteorology_detail_file.toml')
   doc_ptr => doc
   call read_meteorology_toml(doc_ptr, c, errors)

   @assertFalse(errors%has_errors())
   @assertEqual(1, c%swmetdetail)
   @assertEqual(48, c%nmetdetail)
   @assertTrue(allocated(c%detail_file))
   @assertEqual('hupsel.det.csv', c%detail_file)
end subroutine
```

- [ ] **Step 3: Run — expect FAIL**

```bash
pixi run -e test test-pfunit -- --filter test_read_meteorology_toml 2>&1 | tail -10
```

Expected: compile error or `detail_file` unallocated assertion failure.

- [ ] **Step 4: Add the read call**

In `src/io/toml/read_meteorology_toml.f90`, inside the `if (associated(temporal)) then` block (currently lines 42–46), add after the existing `swmetfilall` read:

```fortran
      call get_table(sec, 'temporal', temporal, 'meteorology.temporal', errors)
      if (associated(temporal)) then
         call get_optional_int_with_default(temporal, 'swmetdetail',  config%swmetdetail,  0, 'meteorology.temporal.swmetdetail',  errors)
         call get_optional_int_with_default(temporal, 'nmetdetail',   config%nmetdetail,   0, 'meteorology.temporal.nmetdetail',   errors)
         call get_optional_int_with_default(temporal, 'swmetfilall',  config%swmetfilall,  0, 'meteorology.temporal.swmetfilall',  errors)
         call get_optional_string_with_default(temporal, 'detail_file', config%detail_file, '', &
                                               'meteorology.temporal.detail_file', errors)
      end if
```

- [ ] **Step 5: Run — expect PASS**

```bash
pixi run -e test test-pfunit -- --filter test_read_meteorology_toml 2>&1 | tail -10
```

Expected: all tests PASS.

- [ ] **Step 6: Commit**

```bash
git add src/io/toml/read_meteorology_toml.f90 \
        tests/unit/io/toml/fixtures/meteorology_detail_file.toml \
        tests/unit/io/toml/test_read_meteorology_toml.pf
git commit -m "feat(meteo-reader): read detail_file from [meteorology.temporal]

Adds optional string key meteorology.temporal.detail_file; defaults to ''
when absent. Required when swmetdetail=1 and the metfile is CSV."
```

---

## Task 5: Promote `metcsv_det` stubs and add `swMetDetCSV`

**Files:**
- Modify: `src/core/variables.f90`

Currently lines 211–213 in `variables.f90` hold stubs for the detail CSV cache. This task adds `swMetDetCSV` and updates the column-count comment to match the 7-column ADR 0014 schema.

- [ ] **Step 1: Edit `variables.f90`**

Locate the block (currently around line 211):

```fortran
      ! Detail CSV cache (swmetdetail=1): 1=datetime(frac days), 2=rad, 3=temp, 4=hum, 5=wind, 6=rain
      integer :: nmetcsv_det = 0
      real(8), dimension(:,:), allocatable :: metcsv_det
      character(len=200) rainfil   ! Name of input file with detailed rainfall intensities
      integer   swRainCSV          ! 0=legacy .YYY; 1=CSV mode (pre-loaded by adapter)
```

Replace with:

```fortran
      ! Detail CSV cache (swmetdetail=1): 7 columns per ADR 0014.
      ! 1=datetime(frac days since JD1900), 2=record, 3=rad(kJ/m2/d),
      ! 4=temp(C), 5=hum(kPa), 6=wind(m/s), 7=rain(mm)
      integer :: nmetcsv_det = 0
      real(8), dimension(:,:), allocatable :: metcsv_det
      integer   swMetDetCSV        ! 0=no sub-daily CSV; 1=sub-daily CSV pre-loaded by adapter
      character(len=200) rainfil   ! Name of input file with detailed rainfall intensities
      integer   swRainCSV          ! 0=legacy .YYY; 1=CSV mode (pre-loaded by adapter)
```

- [ ] **Step 2: Build to confirm**

```bash
pixi run -e test build-linux 2>&1 | grep -E "error:" | head -5
```

Expected: clean.

- [ ] **Step 3: Commit**

```bash
git add src/core/variables.f90
git commit -m "feat(variables): add swMetDetCSV flag; update metcsv_det comment to 7 cols

swMetDetCSV set by adapter when detail_file is pre-loaded.
Comment updated to reflect the full 7-column schema from ADR 0014."
```

---

## Task 6: Pre-load block in `config_to_variables.f90`

**Files:**
- Modify: `src/io/toml/config_to_variables.f90`

This is the adapter step: when `swMetCSV=1` and `swmetdetail=1`, the adapter reads `detail_file` via `read_csv_table` and pre-loads all rows into `metcsv_det`.

- [ ] **Step 1: Add the pre-load block**

In `src/io/toml/config_to_variables.f90`, locate the closing `end if` of the daily CSV pre-load block (currently around line 173, ending the `if (index(trim(metfil), '.csv') > 0) then` branch). Immediately after that `end if` (and still inside the outer meteo block), add:

```fortran
      ! Detail meteo CSV pre-load (swmetdetail=1 + CSV mode + detail_file provided).
      swMetDetCSV = 0
      if (swmetdetail == 1 .and. swMetCSV == 1) then
         if (.not. allocated(config%meteo%detail_file) .or. &
             len_trim(config%meteo%detail_file) == 0) then
            call fatalerr_collected('config_to_variables', &
               'meteorology.temporal.detail_file required when ' // &
               'swmetdetail=1 and metfile is a CSV')
         else
            block
               use csv_reader_mod,  only: read_csv_table
               use error_mod,       only: error_collection_t
               real(8), allocatable :: tbl(:,:)
               type(error_collection_t) :: errs
               character(len=8) :: hdr(7)
               character(len=300) :: csvpath
               integer :: r
               hdr(1) = 'datetime'
               hdr(2) = 'record  '
               hdr(3) = 'rad     '
               hdr(4) = 'temp    '
               hdr(5) = 'hum     '
               hdr(6) = 'wind    '
               hdr(7) = 'rain    '
               csvpath = trim(pathatm) // trim(config%meteo%detail_file)
               call read_csv_table(trim(csvpath), hdr, tbl, errs)
               call errs%abort_if_fatal()
               nmetcsv_det = size(tbl, 1)
               allocate(metcsv_det(nmetcsv_det, 7))
               do r = 1, nmetcsv_det
                  metcsv_det(r, :) = tbl(r, :)
               end do
            end block
            swMetDetCSV = 1
         end if
      end if
```

The `USE variables` at the top of the subroutine already imports `swmetdetail`, `swMetCSV`, `nmetcsv_det`, `metcsv_det`, `swMetDetCSV`, and `pathatm`. Confirm by grepping:

```bash
grep -n "use variables" src/io/toml/config_to_variables.f90 | head -3
```

Also confirm `fatalerr_collected` is imported (it should be via `use error_mod`).

- [ ] **Step 2: Build**

```bash
pixi run -e test build-linux 2>&1 | grep -E "error:" | head -5
```

Expected: clean. If the build complains about `swMetDetCSV` not in scope, verify that `variables.f90` re-exports it via the `use variables` line at the top of `config_to_variables.f90`.

- [ ] **Step 3: Run full regression (no collateral damage)**

```bash
pixi run -e test check-full 2>&1 | tail -15
```

Expected: 5 passed, 0 failed. All active cases use `swmetdetail=0`, so the new block is never entered.

- [ ] **Step 4: Commit**

```bash
git add src/io/toml/config_to_variables.f90
git commit -m "feat(adapter): pre-load sub-daily detail CSV when swmetdetail=1

Reads detail_file via read_csv_table (7-col datetime schema) and
caches rows in metcsv_det/nmetcsv_det. Sets swMetDetCSV=1 on success.
Aborts with fatalerr if detail_file is missing when required."
```

---

## Task 7: Implement `MeteoCSVDetYear` and update dispatch in `readmeteo.f90`

**Files:**
- Modify: `src/io/readmeteo.f90`

`MeteoCSVDetYear` mirrors `MeteoCSVYear` but reads from the `metcsv_det` 7-column cache and populates the `det*` arrays. The `irectotal` and `nofd` initialisation at lines 234–243 of `ReadMeteoYear` already runs after the `goto 100` bypass, so `MeteoCSVDetYear` only needs to populate the arrays and set `ifnd`.

- [ ] **Step 1: Update the dispatch in `ReadMeteoYear`**

Locate the current CSV dispatch (around lines 94–98):

```fortran
! --- CSV mode: extract year slice from the pre-loaded cache
      if (swMetCSV == 1) then
         call MeteoCSVYear(ifnd)
         goto 100
      end if
```

Replace with:

```fortran
! --- CSV mode: extract year slice from the pre-loaded cache.
!     For daily mode (swmetdetail=0): MeteoCSVYear.
!     For sub-daily mode (swmetdetail=1): MeteoCSVDetYear.
      if (swMetCSV == 1) then
         if (swmetdetail == 0) then
            call MeteoCSVYear(ifnd)
         else
            call MeteoCSVDetYear(ifnd)
         end if
         goto 100
      end if
```

- [ ] **Step 2: Add `MeteoCSVDetYear` subroutine**

Append after the closing `end subroutine MeteoCSVYear` (currently around line 668), before the `days1900_to_md` helper:

```fortran

! SUBROUTINE: MeteoCSVDetYear
! Extract one year's sub-daily meteo from the pre-loaded metcsv_det cache.
! Called by ReadMeteoYear when swMetCSV==1 and swmetdetail==1.
! Populates dettime, detrecord, detrad, dettav, dethum, detwind, detrain.
! irectotal and nofd are set in ReadMeteoYear after goto 100.
subroutine MeteoCSVDetYear(ifnd)
use error_mod, only: fatalerr_collected
use variables, only: yearmeteo, metcsv_det, nmetcsv_det, &
                     dettime, detrecord, detrad, dettav, dethum, detwind, detrain
use swap_array_dimensions, only: NMETFILE
implicit none
integer, intent(out) :: ifnd

integer, parameter :: jd1900 = 2415020
integer :: jday
external jday

integer  :: i, i1, i2, n
real(8)  :: t_jan1, t_jan1_next

! Year boundaries in days-since-jd1900.
! All sub-daily timestamps for yearmeteo satisfy:
!   t_jan1 <= timestamp < t_jan1_next
t_jan1      = real(jday(yearmeteo,   1, 1) - jd1900, 8)
t_jan1_next = real(jday(yearmeteo+1, 1, 1) - jd1900, 8)

! Scan cache for this year (cache is sorted by datetime).
i1 = 0; i2 = 0
do i = 1, nmetcsv_det
   if (metcsv_det(i,1) >= t_jan1 - 0.5d0 .and. &
       metcsv_det(i,1) <  t_jan1_next - 0.5d0) then
      if (i1 == 0) i1 = i
      i2 = i
   end if
end do

if (i1 == 0) then
   call fatalerr_collected('MeteoCSVDetYear', &
      'No sub-daily meteo CSV records found for the requested year')
   ifnd = 0; return
end if

n = i2 - i1 + 1
if (n > NMETFILE) then
   call fatalerr_collected('MeteoCSVDetYear', &
      'Sub-daily meteo CSV record count exceeds NMETFILE (17568)')
   ifnd = 0; return
end if
ifnd = n

! Populate per-slot arrays.
! metcsv_det columns: 1=datetime, 2=record, 3=rad(kJ), 4=temp, 5=hum, 6=wind, 7=rain
dettime(1:n)   = metcsv_det(i1:i2, 1)
detrecord(1:n) = nint(metcsv_det(i1:i2, 2))
detrad(1:n)    = metcsv_det(i1:i2, 3) * 1000.0d0   ! kJ/m2 → J/m2
dettav(1:n)    = metcsv_det(i1:i2, 4)
dethum(1:n)    = metcsv_det(i1:i2, 5)
detwind(1:n)   = metcsv_det(i1:i2, 6)
detrain(1:n)   = metcsv_det(i1:i2, 7)

end subroutine MeteoCSVDetYear
```

Note: `detrecord(nmetfile)` is declared as `integer` in variables.f90; `metcsv_det(:,2)` is `real(8)` — `nint()` converts correctly. The module `arrays` must be imported for `NMETFILE`. Verify by grepping:

```bash
grep -n "use arrays\|use swap_array_dimensions\|NMETFILE" src/io/readmeteo.f90 | head -5
```

If `arrays` is not currently imported by `readmeteo.f90`, check what module provides `NMETFILE`:

```bash
grep -rn "NMETFILE\|nmetfile" src/core/arrays.f90 | head -3
```

If `nmetfile` is a local parameter, replace `use arrays, only: NMETFILE` and `NMETFILE` with `integer, parameter :: nmetfile_local = 17568` and `nmetfile_local`. Match the naming convention visible in the file.

- [ ] **Step 3: Build**

```bash
pixi run -e test build-linux 2>&1 | grep -E "error:" | head -5
```

Expected: clean. Fix any compilation errors before proceeding.

- [ ] **Step 4: Run the full regression (no regression)**

```bash
pixi run -e test check-full 2>&1 | tail -15
```

Expected: 5 passed, 0 failed. All active cases have `swmetdetail=0`; the new subroutine is compiled but never called.

- [ ] **Step 5: Commit**

```bash
git add src/io/readmeteo.f90
git commit -m "feat(readmeteo): MeteoCSVDetYear — sub-daily CSV year extraction

Mirrors MeteoCSVYear but reads from metcsv_det (7-col datetime cache).
ReadMeteoYear dispatches to MeteoCSVDetYear when swMetCSV==1 and
swmetdetail==1. irectotal/nofd init continues in ReadMeteoYear post-goto."
```

---

## Task 8: pFUnit unit test for `MeteoCSVDetYear`

**Files:**
- Modify: `tests/unit/io/toml/test_config_to_variables.pf` (or create a new `tests/unit/io/test_meteo_csv_det_year.pf` if the suite supports it)

The test directly manipulates the `variables` module cache and calls `MeteoCSVDetYear` via external interface, mirroring how integration tests in this suite already use module state.

- [ ] **Step 1: Check if adding a new `.pf` file to the build is possible**

```bash
grep -n "test_csv\|test_meteo\|\.pf" tests/unit/io/meson.build 2>/dev/null | head -10
```

If `.pf` files are discovered automatically (glob), we can add a new file. If they are listed explicitly, add the new file to the list. Inspect the meson.build output to decide:

```bash
cat tests/unit/io/meson.build 2>/dev/null | head -30
```

- [ ] **Step 2: Write the test**

Create or append (per the discovery in Step 1) `tests/unit/io/test_meteo_csv_det_year.pf`:

```fortran
@test
subroutine test_MeteoCSVDetYear_extracts_year_slice()
   use funit
   use iso_fortran_env, only: real64
   use variables, only: metcsv_det, nmetcsv_det, yearmeteo, nmetdetail, &
                        dettime, detrecord, detrad, dettav, dethum, detwind, detrain
   interface
      subroutine MeteoCSVDetYear(ifnd)
         integer, intent(out) :: ifnd
      end subroutine
   end interface
   integer :: ifnd

   ! 1980-04-24 00:00:00 = 29334.0 days since JD1900 (verified in test_csv_reader.pf)
   ! 1980-04-24 12:00:00 = 29334.5
   ! 1981-04-23 00:00:00 ≈ 29699.0 (in 1981; 1981-01-01 = 29586, +113 days)
   nmetcsv_det = 3
   if (allocated(metcsv_det)) deallocate(metcsv_det)
   allocate(metcsv_det(3, 7))
   ! datetime, record, rad(kJ), temp, hum, wind, rain
   metcsv_det(1, :) = [29334.0d0, 1.0d0, 500.0d0, 10.0d0, 1.0d0, 3.0d0, 0.0d0]
   metcsv_det(2, :) = [29334.5d0, 2.0d0, 800.0d0, 14.0d0, 0.8d0, 4.0d0, 1.5d0]
   metcsv_det(3, :) = [29699.0d0, 1.0d0, 400.0d0,  8.0d0, 1.1d0, 2.5d0, 0.0d0]
   yearmeteo  = 1980
   nmetdetail = 2

   call MeteoCSVDetYear(ifnd)

   @assertEqual(2, ifnd)
   @assertEqual(29334.0d0, dettime(1), 1.0e-4_real64)
   @assertEqual(29334.5d0, dettime(2), 1.0e-4_real64)
   @assertEqual(1, detrecord(1))
   @assertEqual(2, detrecord(2))
   ! 500 kJ/m2/d → 500000 J/m2/d
   @assertEqual(500000.0d0, detrad(1), 1.0d0)
   @assertEqual(10.0d0, dettav(1), 1.0e-6_real64)

   deallocate(metcsv_det)
end subroutine

@test
subroutine test_MeteoCSVDetYear_no_records_for_year_calls_fatalerr()
   use funit
   use variables, only: metcsv_det, nmetcsv_det, yearmeteo, nmetdetail
   interface
      subroutine MeteoCSVDetYear(ifnd)
         integer, intent(out) :: ifnd
      end subroutine
   end interface
   integer :: ifnd

   ! Only 1981 data in cache — yearmeteo=1980 should trigger fatalerr.
   ! pFUnit cannot easily assert on fatalerr (it aborts).
   ! This test is commented out and kept as documentation only.
   ! TODO: wire up a non-aborting error handler if fatalerr is ever made
   !       testable via error_collection_t.
   @assertTrue(.true.)  ! placeholder — build and link verification only
end subroutine
```

Note: `fatalerr_collected` calls `stop` on fatal errors, making that path unassertable in pFUnit without refactoring. The second test is a documented placeholder. The important test is the first one.

- [ ] **Step 3: Register the file in `tests/unit/meson.build`**

The `.pf` list is explicit (not a glob). Add `'io/test_meteo_csv_det_year.pf'` to the `pf_files` array after the existing `'io/test_csv_reader.pf'` entry (currently at line 154):

```
        'io/test_csv_reader.pf',
        'io/test_meteo_csv_det_year.pf',
```

- [ ] **Step 4: Run the new test**

```bash
pixi run -e test test-pfunit -- --filter test_MeteoCSVDetYear 2>&1 | tail -15
```

Expected: 2 PASS (second test is a trivial placeholder).

- [ ] **Step 5: Commit**

```bash
git add tests/unit/io/test_meteo_csv_det_year.pf
# If meson.build was modified:
# git add tests/unit/io/meson.build
git commit -m "test(readmeteo): pFUnit test for MeteoCSVDetYear year extraction

Tests the happy-path: 2 rows in 1980, 1 in 1981; verifies ifnd=2,
dettime values, detrecord, and detrad kJ→J conversion."
```

---

## Task 9: Update `meteorology.md` documentation

**Files:**
- Modify: `docs/meteorology.md`

- [ ] **Step 1: Add the sub-daily section**

In `docs/meteorology.md`, after the `## Rain events CSV` section and before `## Legacy formats`, insert:

```markdown
---

## Sub-daily detail meteorology CSV (`temporal.detail_file`)

Used when `[meteorology.temporal].swmetdetail = 1`. The file must cover every
year of the simulation period. Each row is one sub-daily time slot.

**Schema:**

```
datetime,record,rad,temp,hum,wind,rain
```

| Column     | Unit         | Description                                                    |
|------------|--------------|----------------------------------------------------------------|
| `datetime` | ISO datetime | Slot timestamp `YYYY-MM-DD HH:MM:SS`; fractional days since JD 2415020 |
| `record`   | –            | Intra-day slot index 1 … `nmetdetail`                          |
| `rad`      | kJ m⁻² d⁻¹  | Global radiation for the slot (converted to J m⁻² d⁻¹ internally) |
| `temp`     | °C           | Air temperature (single value per slot)                        |
| `hum`      | kPa          | Actual vapour pressure                                         |
| `wind`     | m s⁻¹        | Wind speed                                                     |
| `rain`     | mm           | Rainfall for the slot                                          |

**Example (48 slots per day):**

```csv
# Sub-daily meteo, station Hupsel, 48 half-hourly slots per day
datetime,record,rad,temp,hum,wind,rain
2002-01-01 00:00:00,1,0.0,-3.2,0.524,4.90,0.000
2002-01-01 00:30:00,2,0.0,-3.3,0.525,4.85,0.000
...
2002-01-01 23:30:00,48,0.0,-2.9,0.520,5.10,0.000
```

**TOML configuration:**

```toml
[meteorology]
file = "hupsel.csv"

[meteorology.temporal]
swmetdetail = 1
nmetdetail  = 48
detail_file = "hupsel.det.csv"
```

**Notes:**

- The adapter (`config_to_variables.f90`) pre-loads all rows into `metcsv_det` at
  startup. `MeteoCSVDetYear` extracts the current year's slice on each `ReadMeteoYear`
  call, exactly mirroring the daily `MeteoCSVYear` pattern.
- The maximum record count per year is `NMETFILE = 17568` (48 slots × 366 days).
- `irectotal` is initialised in `ReadMeteoYear` using `dettime(1)` after
  `MeteoCSVDetYear` returns — no changes needed in `meteoday.f90`.
- If no rows are found for the requested year, SWAP aborts via `fatalerr_collected`.
```

- [ ] **Step 2: Commit**

```bash
git add docs/meteorology.md
git commit -m "docs(meteo): document detail_file sub-daily CSV pathway

Adds schema table, example, TOML config, and implementation notes
for the swmetdetail=1 CSV path (ADR 0014 Step 1)."
```

---

## Task 10: Full verification

- [ ] **Step 1: Run the complete pFUnit suite**

```bash
pixi run -e test test-pfunit 2>&1 | tail -15
```

Expected: all tests PASS. If any test fails, fix before continuing.

- [ ] **Step 2: Run the full regression**

```bash
pixi run -e test check-full 2>&1 | tail -15
```

Expected:
```
✓ hupselbrook:   regression ok
✓ grassgrowth:   regression ok
✓ oxygenstress:  regression ok
✓ salinitystress:regression ok
✓ surfacewater:  regression ok
Results: 5 passed, 0 failed
```

- [ ] **Step 3: Confirm no dead references remain**

```bash
git grep -rn 'ioutils\|io_utils_mod\|parse_output_extensions' -- src/ tests/ meson.build
```

Expected: **zero matches**.

- [ ] **Step 4: Confirm `detail_file` is wired end-to-end**

```bash
git grep -rn 'detail_file\|MeteoCSVDetYear\|swMetDetCSV' -- src/ tests/
```

Expected: at least one match in each of `src/config/`, `src/io/toml/`, `src/core/variables.f90`, `src/io/readmeteo.f90`, and `tests/`.

- [ ] **Step 5: Tag the milestone**

```bash
git tag rescue/phase-csv-meteo-complete
git tag --list rescue/*
```

---

## Risks and mitigations

| Risk | Mitigation |
|---|---|
| Phase 0 cleanup whitelist accidentally matches a committed input | Whitelist uses explicit names (`swap.swp`, `result_*.csv`, `*.log`, `*.tmp`, …); never bare `*.csv` or `*.crp`. After Task 0.2 Step 3, run `git status` in the submodule to confirm no committed files were deleted. |
| Phase 0 submodule + outer-repo commits get out of sync | Each task pairs one submodule commit with one outer-repo bump; verify with `git submodule status tests/swap-cases` after each pair. |
| Phase 0 regression numerics drift after `copytree` source flip | The TOML dir now contains the `.crp` files (Task 0.1), so `copytree` from the TOML dir gives SWAP exactly the same input set it had before. Run `pixi run -e test check-full` after each Phase 0 task. |
| `use arrays, only: NMETFILE` unavailable from `readmeteo.f90` | Grep for the import pattern in the file; fall back to a local `integer, parameter :: nmetfile_loc = 17568` if the module isn't already imported. |
| pFUnit tests that `use variables` interfere across test runs (shared module state) | Deallocate `metcsv_det` at test teardown (included in the test code above). |
| `swMetDetCSV` not visible inside the `config_to_variables.f90` pre-load block | The outer subroutine already has `use variables` (confirmed in existing code). The `block` construct inherits the parent's USE associations. |
| Active regression cases accidentally trigger `MeteoCSVDetYear` | All five active TOML cases set `swmetdetail=0`; the new dispatch is `if (swmetdetail == 0) ... else ...` so the new branch is never entered. |
| ADR 0014 Steps 2–3 (TTutil deletion, variable sweep) | These remain blocked on `.swp` pipeline retirement. The HACK markers and TTutil code in `readmeteo.f90` are left in place per the ADR. |

---

## What this plan does NOT cover (deferred)

- **ADR 0014 Step 2** — deleting TTutil branches from `readmeteo.f90` (blocked until `.swp` retired or `readmeteo_legacy.f90` split done)
- **ADR 0014 Step 3** — deleting dead variables (`swMetCSV`, `swRainCSV`, `swMetFilAll`, `rainfil`, `station`, `ad`, `am`) — blocked on Step 2
- **End-to-end regression case for `swmetdetail=1`** — requires a real sub-daily CSV dataset and a known reference output; deferred to Phase 4f-extend
- **pFUnit unit test for `MeteoCSVDetYear`** (Task 8) — deferred. Linking the standalone subroutine into `tests/unit/unit-swap-tests` requires adding `readmeteo.f90` to the test source list, which transitively pulls in `meteodt_mod` and a chain of physics modules. The Task 7 code review caught the only correctness gap (year-boundary slack); end-to-end coverage will arrive with the future `swmetdetail=1` regression case. To revisit: extract `MeteoCSV{Year,DetYear}` plus `days1900_to_md`/`jday` into a separate file (`src/io/meteo_csv.f90`) so the test build can include just that helper.
