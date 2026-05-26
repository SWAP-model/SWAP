# Regression Coverage Expansion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add clustered regression test cases that activate `SW*` model options dormant across the existing 6 cases, each proving the modern gfortran build reproduces the SWAP 4.2.0 reference (`swap420`).

**Architecture:** Approach B (subsystem clusters). Each new case is a legacy-ASCII + TOML pair; reference fixture from `swap420`, gfortran baseline from the modern build, cross-checked within `TOL=1e-2`; registered in `test_output_regression.py`'s `CASES`. Phase 0 triage gates which clusters are viable.

**Tech Stack:** Fortran (SWAP), TOML configs, Python regression harness (`tests/regression/test_output_regression.py`), pixi tasks (`swap-ref`, `check-full`), gfortran build in `builddir/`.

---

## Task 1: Phase 0 — Feasibility triage

**Files:**
- Create: `dev-docs/superpowers/plans/2026-05-27-option-triage.md` (the triage table)

- [ ] **Step 1: Verify swap420 runs and inspect its CSV columns**

```bash
cd /home/zawadzkim/Code/swap
ldd tests/reference/swap420 | grep -i "not found" || echo "deps OK"
pixi run swap-ref -c hupselbrook -k 2>&1 | tail -20
head -3 tests/swap-cases/1.hupselbrook/result_output.csv   # column header line
```
Expected: swap420 exits, `result_output.csv` produced. Record its header columns.

- [ ] **Step 2: For each uncovered option, classify support in the modern build**

For every option in the spec's Category A/B lists, grep the modern source (excluding
`/toml/` for the core count) and the TOML reader:
```bash
for sw in swsnow swfrost swhyst swsophy swdc swdrought swcalt swsalinity swinter \
          swgc swqhbot swqhr swco2 swnrsrf swdislay swgerm swharv swsow swprep \
          swrdc swtsum swdivide swcrop swrootradius swkmean swdmgrz swdmmow \
          swlossgrz swlossmow swkimpl swsublim swgraz swman swtill swrain \
          swmetdetail swetsine swinco swbotbc swbotbhea swsrf swsec swallo \
          swcompensate swcirrthres; do
  core=$(grep -rliE "$sw" src --include=*.f90 | grep -v '/toml/' | wc -l)
  toml=$(grep -rliE "$sw" src/io/toml --include=*.f90 | wc -l)
  printf "%-14s core=%-3s toml=%-3s\n" "$sw" "$core" "$toml"
done
```
Mark `core=0` as **unsupported** (exclude). Spot-check `core=1` (e.g. `swdc`) by
reading the one hit to confirm it is wired into a compute path, not read-and-ignored.

- [ ] **Step 3: Confirm which output columns each cluster's physics moves**

```bash
grep -niE "snow|ssnow|tfrost|hysteresis|cropfra|gcrop" src/io/toml/write_swap_config.f90 \
  scripts/output_parity_check.py 2>/dev/null
# Identify CSV column names emitted for snow/frost/hysteresis/cover-fraction etc.
```
Confirm whether `swap420`'s CSV header (Step 1) carries those columns. Columns present
only in the modern build → mark the case for golden-master fallback on those columns.

- [ ] **Step 4: Write the triage table and finalize the cluster set**

Write `dev-docs/superpowers/plans/2026-05-27-option-triage.md` with one row per option:
`option | support (yes/no/marginal) | base case | moved columns | swap420 emits? | cluster`.
Drop unsupported options. Split/merge clusters whose options proved incompatible.

- [ ] **Step 5: Commit**

```bash
git add dev-docs/superpowers/plans/2026-05-27-option-triage.md
git commit -m "docs(plan): Phase 0 option-support triage for regression expansion"
```

---

## Task 2: Phase 1 — soil-hysteresis case (recipe proof, end-to-end)

Proves the full recipe on one cluster before fan-out. Base = hupselbrook (swap420
supports hysteresis; no modern-only columns).

**Files:**
- Create: `tests/swap-cases/7.soilhysteresis/` (clone of `1.hupselbrook/`, ASCII)
- Create: `tests/swap-cases/toml/7.soilhysteresis/` (clone of `toml/1.hupselbrook/`)
- Create: `tests/regression/soilhysteresis_expected.json` (swap420 reference)
- Create: `tests/regression/soilhysteresis_expected_gfortran.json` (modern baseline)
- Modify: `tests/regression/test_output_regression.py` (add `CaseConfig`)

- [ ] **Step 1: Clone the base case in both formats**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
cp -r 1.hupselbrook 7.soilhysteresis
cp -r toml/1.hupselbrook toml/7.soilhysteresis
# remove any stale generated outputs from the clone
rm -f 7.soilhysteresis/*.tmp 7.soilhysteresis/*_swap.log 7.soilhysteresis/swap.swp \
      7.soilhysteresis/fort.20 7.soilhysteresis/reruns.log 7.soilhysteresis/swap_debug.log
rm -f toml/7.soilhysteresis/result_output.csv toml/7.soilhysteresis/Swap.ok
```

- [ ] **Step 2: Activate SWHYST in the ASCII template**

In `7.soilhysteresis/swap_linux.swp.template`, set `SWHYST = 1` and add the companion
hysteresis parameters. Read the template's own commented block for the exact parameter
names (`TAU`, and per-layer `ALFAW`); copy the canonical example values documented there
or in the upstream `legacy/swap-4.2.0` sample. Set `SWKIMPL = 1` and pick one alternate
`SWKMEAN` (e.g. `2`). Mirror the same in the staged `swap.swp` is automatic (run_case
copies the template).

- [ ] **Step 3: Activate SWHYST in the TOML twin**

In `toml/7.soilhysteresis/swap.toml`, set the equivalent keys (`swhyst`, `tau`,
per-layer `alfaw`, `swkimpl`, `swkmean`). Use `read_soil_toml.f90` to confirm the exact
TOML key names the modern reader expects.

- [ ] **Step 4: Generate the swap420 reference fixture**

```bash
cd /home/zawadzkim/Code/swap
pixi run swap-ref -c soilhysteresis -k
python3 - <<'PY'  # aggregate using the harness's own aggregate()
import sys; sys.path.insert(0,'tests/regression')
from test_output_regression import aggregate
import json
a,t,m = aggregate(__import__('pathlib').Path('tests/swap-cases/7.soilhysteresis/result_output.csv'),
                  ["RAIN","INTERC","RUNOFF","EPOT","EACT","DRAINAGE","QBOTTOM","TPOT","TACT","DSTOR"],
                  ["GWL"], [])
json.dump({"years":a,"total":t,"mean":m}, open('tests/regression/soilhysteresis_expected.json','w'),
          indent=2, sort_keys=True)
print("wrote reference")
PY
```
Expected: a `soilhysteresis_expected.json` with non-zero annual blocks.

- [ ] **Step 5: Add the CaseConfig and generate the gfortran baseline**

In `tests/regression/test_output_regression.py`, add to `CASES`:
```python
    "soilhysteresis": CaseConfig(
        name="soilhysteresis",
        case_dir="7.soilhysteresis",
        fixture="soilhysteresis_expected_gfortran.json",
        flux_vars=["RAIN", "INTERC", "RUNOFF", "EPOT", "EACT",
                   "DRAINAGE", "QBOTTOM", "TPOT", "TACT", "DSTOR"],
        state_vars=["GWL"],
    ),
```
Then:
```bash
pixi run build-linux   # rebuild if any src touched (not expected here)
python3 tests/regression/test_output_regression.py --regenerate-fixtures soilhysteresis
```
Expected: `soilhysteresis_expected_gfortran.json` written.

- [ ] **Step 6: Cross-check gate — reference vs gfortran agree within TOL**

```bash
python3 - <<'PY'
import json, math
a=json.load(open('tests/regression/soilhysteresis_expected.json'))
b=json.load(open('tests/regression/soilhysteresis_expected_gfortran.json'))
bad=[]
for blk in ("years","total","mean"):
    ea,eb=a.get(blk,{}),b.get(blk,{})
    def walk(x,y,path=""):
        if isinstance(x,dict):
            for k,v in x.items(): walk(v,y.get(k,{}),f"{path}.{k}")
        else:
            if not (isinstance(y,(int,float)) and math.isclose(x,y,abs_tol=1e-2)):
                bad.append((path,x,y))
    walk(ea,eb,blk)
print("DIVERGENCES:",bad if bad else "none (within 1e-2)")
PY
```
Expected: `none`. If divergences appear, confirm hysteresis actually engaged (compare
to base hupselbrook fixture — values must differ from base), then decide: misconfig vs
real divergence → document in `INVESTIGATION_NOTES.md`.

- [ ] **Step 7: Run the regression for the new case + full check**

```bash
python3 tests/regression/test_output_regression.py soilhysteresis
pixi run check-full
```
Expected: `✓ soilhysteresis: regression ok`; check-full green.

- [ ] **Step 8: Commit**

```bash
git add tests/swap-cases/7.soilhysteresis tests/swap-cases/toml/7.soilhysteresis \
        tests/regression/soilhysteresis_expected.json \
        tests/regression/soilhysteresis_expected_gfortran.json \
        tests/regression/test_output_regression.py
git commit -m "test(regression): add soil-hysteresis case (SWHYST/SWKIMPL/SWKMEAN)"
```

> NOTE: `tests/swap-cases/` is a git submodule. Commit inside it first, then bump the
> superproject pointer. Confirm during execution and adjust the `git add` paths.

---

## Tasks 3–N: Cluster fan-out

Each surviving cluster from the Task 1 triage table becomes one task following the
**exact recipe of Task 2** (clone both formats → activate switches + companion blocks →
swap420 reference → CaseConfig + gfortran baseline → cross-check gate → regression +
check-full → commit). Per-cluster specifics below; companion-block values and asserted
columns are taken from the triage table and the existing cases at execution time.

### Task 3: winter (`SWSNOW`, `SWFROST`, `SWSUBLIM`)
- Base: a case with cold meteo (hupselbrook 283 or grass 260 — choose the one with
  sub-zero temps). Asserted columns: add `SNOW`/`SSNOW` (+ frost indicator) **iff**
  swap420 emits them (triage Step 3); else golden-master those columns.

### Task 4: soil-tabulated (`SWSOPHY`)
- Base: hupselbrook. Replaces analytic Mualem–van Genuchten with tabulated soil-physical
  functions → own case (incompatible with hysteresis). Needs the soil-table input block.

### Task 5: drought-vanlier (`SWDROUGHT=2`, `SWRDC`, `SWROOTRADIUS=1`, `SWDIVIDE=0`)
- Base: oxygenstress. De Jong van Lier needs its parameter block + root radius input.
  Assert `TREDDRY`/root-extraction columns.

### Task 6: crop-calendar (`SWSOW`, `SWPREP`, `SWHARV`, `SWGERM=1`, `SWCO2`, `SWGC=2`)
- Base: a fixed-crop case (maize/potato from hupselbrook). `SWGC=2` swaps `LAI`→`GCTB`.
  Assert crop dry-matter / DVS columns.

### Task 7: grass-management (`SWDMGRZ=1`/`SWDMMOW=1`, `SWLOSSGRZ`, `SWLOSSMOW`, `SWTSUM` alt)
- Base: grassgrowth. Assert `PGRASSDM`/`GRASSDM`/`PMOWDM`/`MOWDM` (cumulative).

### Task 8: salinity-osmotic (`SWSALINITY=2`, `SWDC` if supported, `SWBOTBC=2`)
- Base: salinitystress. Assert `TREDSOL`/`CONC[...]` columns.

### Task 9: extended-drainage (`SWNRSRF`, `SWDISLAY`, `SWSEC=1`, `SWSRF` alt, `SWALLO=2`, `SWQHBOT=1`, `SWQHR=2`, `SWBOTBHEA=2`)
- Base: surfacewater. Assert `GWL`/`POND`/drainage-flux columns. Split if options conflict.

### Task 10: meteo-detail (`SWRAIN=1`, `SWMETDETAIL=1`, `SWETSINE`, `SWINCO=1`)
- Base: hupselbrook. `SWRAIN=1` needs mean-intensity column in meteo; `SWMETDETAIL=1`
  needs sub-daily meteo. Assert `RAIN`/`RUNOFF`/`INTERC`.

Each: one commit on `development`, both input dirs + both fixtures + `CASES` entry; run
the new case + `check-fast` before committing (per-task regression gate).

---

## Self-review notes

- **Spec coverage:** Task 1 ⇒ Phase 0 triage; Task 2 ⇒ Phase 1 recipe proof; Tasks 3–10
  ⇒ the 9 provisional clusters (winter, soil-tabulated, drought-vanlier, crop-calendar,
  grass-management, salinity-osmotic, extended-drainage, meteo-detail + Task 2's
  soil-hysteresis). Cross-check gate, fixture pair, column-selection, and golden-master
  fallback all appear in Task 2 and are referenced by the fan-out. ✓
- **Submodule caveat:** `tests/swap-cases` is a submodule (per repo `.git` file). Commit
  order (submodule then superproject pointer) is flagged in Task 2 and applies to all
  case-creating tasks.
- **Placeholder note:** Companion parameter-block *values* are intentionally sourced at
  execution time from the existing cases / upstream samples rather than transcribed
  here, because they are large per-layer tables; the plan specifies exactly where each
  comes from. This is a config-generation arc, not algorithm code.
