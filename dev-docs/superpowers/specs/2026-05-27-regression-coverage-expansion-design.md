# Regression Coverage Expansion — Activating Dormant Model Options

**Date:** 2026-05-27
**Status:** Approved (brainstorming complete)
**Branch:** `development`

## Problem

The SWAP regression suite (`tests/regression/test_output_regression.py`) runs only
5 of the 6 historical cases (macropore is retired per ADR 0040). Across those cases,
a large fraction of the model's `SW*` option switches are never activated, or are
activated only in one of their several valid modes. Whole subsystems (snow, frost,
hysteresis, tabulated soil physics, De Jong van Lier drought, osmotic salinity, …)
run in **no** regression case, so a modernization change that silently breaks them
would not be caught.

A switch-coverage audit of the six legacy cases' input files
(`swap_linux.swp.template`, `*.crp`, `swap.dra`) against the full mode domains
documented in the template comment blocks produced three buckets:

- **Category A — off in all cases** (~24 switches with real alternatives never run):
  `SWSNOW`, `SWFROST`, `SWHYST`, `SWSOPHY`, `SWDC`, `SWSP`, `SWNRSRF`, `SWDISLAY`,
  `SWCO2`, `SWCIRRTHRES`, `SWRDC`, `SWHARV`, `SWSOW`, `SWPREP`, `SWLOSSGRZ`,
  `SWLOSSMOW`, `SWKIMPL`, `SWPONDMX`, `SWBOTB3RESVERT`, plus output-only toggles.
- **Category B — multi-mode, only some modes exercised** (~27 switches): e.g.
  `SWDROUGHT=2` (van Lier), `SWCALT=1`, `SWINCO=1`, `SWGC=2`, `SWGERM=1`, `SWINTER=2`,
  `SWKMEAN=2..6`, `SWQHBOT=1`, `SWQHR=2`, `SWRAIN=1`, `SWREDU=0`, `SWROOTRADIUS=1`,
  `SWSALINITY=2`, `SWSEC=1`, `SWSRF=1/3`, `SWSTRESSOR=2/4/5`, `SWTOPBHEA=2`,
  `SWTSUM=0/2`, `SWDMGRZ=1`, `SWDMMOW=1`, `SWDIVIDE=0`, `SWCROP=0`, `SWALLO*=2`,
  `SWBOTBC=2`, `SWBOTBHEA=2`, `SWCOMPENSATE=2`.
- **Category C — read but never set in any case** (default-only): `SWHYDRLIFT`,
  `SWCAPRISE`, `SWGRAZ`, `SWMAN`, `SWTILL`, `SWSOYBEAN`, `SWSSDI`, `SWSUBLIM`, …

## Goal

Expand the regression suite to exercise these dormant options by adding new test
cases that activate them, **grouped into coherent subsystem clusters** (Approach B).
Each new case proves the modern (gfortran) build reproduces the upstream SWAP 4.2.0
reference (`tests/reference/swap420`) for the newly-activated physics.

Non-goal: re-introducing physics the rescue deleted (macropore). Options the modern
build does not implement are documented as known gaps, not built.

## Constraints discovered

1. **Every case is a pair.** The reference fixture (`<name>_expected.json`) is
   produced by `swap420`, which reads **legacy ASCII** input
   (`swap_linux.swp.template` + `*.crp` + `swap.dra`). The regression itself runs the
   **modern** build (`builddir/swap`) on the **TOML twin**
   (`tests/swap-cases/toml/<N>.<name>/`), producing `<name>_expected_gfortran.json`.
   Both must exist and be configuration-equivalent.

2. **Activating a switch pulls in its companion parameter block.** "Minimal edit" is
   misleading: `SWHYST=1` needs hysteresis parameters per soil layer, `SWSNOW=1` needs
   snow parameters + cold meteo, `SWGC=2` needs a `GCTB` table instead of `LAI`, etc.

3. **Not every option is reproducible.** The modern build is a rescue branch with
   deleted physics. Confirmed: `SWSP` (solute adsorption) has **no** modern
   implementation (`grep src` → core=0, toml=0); macropore is deleted (ADR 0040);
   `SWDC` is marginal (one reference site each) and must be verified as wired before
   use. A go/no-go triage is mandatory.

## Workflow per case (the recipe)

1. Clone the nearest existing case in **both** formats (legacy dir + `toml/` dir).
2. Activate the cluster's switches and add each switch's companion parameter block in
   both the ASCII inputs and the TOML twin, kept equivalent.
3. Reference: `pixi run swap-ref -c <name>` (runs `swap420` on the legacy dir) →
   aggregate its `result_output.csv` → write `<name>_expected.json`.
4. gfortran baseline: run `builddir/swap` on the TOML twin (or
   `test_output_regression.py --regenerate-fixtures <name>`) → `<name>_expected_gfortran.json`.
5. **Cross-check gate:** `_expected.json` and `_expected_gfortran.json` must agree
   within `TOL = 1e-2` on the asserted columns. Agreement *is* the fidelity proof.
   A genuine divergence is surfaced and triaged like the existing macropore/MOWDM
   notes in `INVESTIGATION_NOTES.md` — not silently tolerated.
6. Register a `CaseConfig` in `test_output_regression.py` whose
   `flux_vars`/`state_vars`/`cumul_vars` include the **columns the new physics moves**
   (e.g. snow → `SNOW`, `SSNOW`), so the fixture asserts the feature rather than an
   unchanged water balance.

### Fallback (Approach C, documented per use)

If `swap420` cannot run a scenario, or its `result_output.csv` lacks a column the
modern build added (snow/crop parity columns), that column/case falls back to
**golden-master**: the modern build's own output is the reference, catching future
drift only. Each fallback is documented inline in the `CaseConfig` and the spec's
triage table.

## Phase 0 — Triage (gates everything)

Classify every uncovered option as **supported** / **unsupported** / **marginal**
in the modern build, map each to its nearest base case, and record the output
column(s) it moves and whether `swap420` can emit them. Output: a triage table that
finalizes the cluster set. Also verify `swap420` actually runs on this machine and
what columns its CSV carries.

## Phase 1 — Recipe proof (one cluster)

Build the **soil-hysteresis** cluster end-to-end (single base = hupselbrook,
`swap420` definitely supports it, no modern-only columns): clone both formats,
activate `SWHYST` (+ `SWKMEAN` alt, `SWKIMPL`), generate both fixtures, confirm the
cross-check gate, register in `CASES`, run the full regression green. This validates
the per-case cost and the gate before fan-out.

## Phase 2+ — Cluster fan-out

Provisional clusters (triage may split/merge/drop):

| Cluster | Base case | Switches activated |
|---|---|---|
| winter | hupselbrook/grass + cold meteo | `SWSNOW`, `SWFROST`, `SWSUBLIM` |
| soil-hysteresis | hupselbrook | `SWHYST`, `SWKMEAN`(alt), `SWKIMPL` |
| soil-tabulated | hupselbrook | `SWSOPHY` (own case — replaces analytic functions) |
| drought-vanlier | oxygenstress | `SWDROUGHT=2`, `SWRDC`, `SWROOTRADIUS=1`, `SWDIVIDE=0` |
| crop-calendar | fixed-crop | `SWSOW`, `SWPREP`, `SWHARV`, `SWGERM=1`, `SWCO2`, `SWGC=2` |
| grass-management | grassgrowth | `SWDMGRZ=1`/`SWDMMOW=1`, `SWLOSSGRZ`, `SWLOSSMOW`, `SWTSUM` alt, `SWGRAZ/SWMAN/SWTILL` |
| salinity-osmotic | salinitystress | `SWSALINITY=2`, `SWDC`(if supported), `SWBOTBC=2` |
| extended-drainage | surfacewater | `SWNRSRF`, `SWDISLAY`, `SWSEC=1`, `SWSRF` alt, `SWALLO=2`, `SWQHBOT=1`, `SWQHR=2`, `SWBOTBHEA=2` |
| meteo-detail | hupselbrook | `SWRAIN=1`, `SWMETDETAIL=1`, `SWETSINE`, `SWINCO=1` |

Each cluster: one case-pair, one commit on `development`, both input dirs + both
fixtures + the `CASES` entry together. Each cluster runs its own new case +
`check-fast` before its commit (per the per-task regression-gate convention), not a
deferred end-of-arc gate.

## Testing / acceptance

- Each new case is green in `test_output_regression.py` against its
  `_expected_gfortran.json`.
- Each new case's `_expected.json` (swap420) and `_expected_gfortran.json` agree
  within `TOL`, or the divergence is documented in `INVESTIGATION_NOTES.md`.
- `pixi run check-full` stays green across the arc.
- The asserted columns for each case demonstrably change when the activated switch is
  toggled off (the fixture tests the feature, not just survival).

## Out of scope

- Macropore and other deleted physics.
- `SWSP` (no modern implementation).
- Refactoring the regression harness beyond adding `CaseConfig` entries and any small
  helper needed for new-column aggregation.
