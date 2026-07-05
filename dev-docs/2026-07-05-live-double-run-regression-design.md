# Live double-run regression — design & plan

**Status:** Accepted (2026-07-05) — decisions locked; implementation pending.
**Relates to:** ADR 0054 (test-suite consolidation & pytest standardization);
supersedes the *stored-fixture* mechanism it describes. Will be recorded by a
new superseding ADR (0055) when Phase 2 lands.

**Locked decisions (2026-07-05):**
- Relocate the oracle binary to `swap-testcases` (do it, not defer — Phase 3 is
  in scope).
- Drop all three fixture helper scripts (`regen_reference.py`,
  `regen_expected.py`, `_switch_validate.py`) — obsolete once comparison is live.
- Make the reference **pluggable** (a selectable descriptor, default
  `swap420gf`), so a future SWAP release can serve as the reference too.
- Keep the ifx `*_expected.json` for provenance (CLAUDE.md); drop the vestigial
  `*_expected_gfortran.json`. No committed oracle cache (ephemeral only, if ever).
**Motivating discussion:** the "where do the `*_reference_gf.json` fixtures
belong" problem — a fixture is a function of *(oracle binary × case inputs)*,
the inputs live in `swap-testcases`, and there is no non-awkward home for the
stored outputs. MODFLOW dissolves the analogous problem by never storing
expected outputs: it runs a reference build and the current build and compares
live.

## Summary

Stop storing `*_reference_gf.json`. For each regression case, run **both** the
modern build (on the TOML inputs) **and** the `swap420gf` 4.2.0 oracle (on the
legacy ASCII inputs) at test time, aggregate both with the existing
`aggregate()`, and compare live with the existing `compare()` / `TOL`. The
"expected" side is regenerated every run instead of loaded from disk.

This is a small change — we already run the oracle-and-aggregate path in
`regen_reference.py`; we just call it at test time instead of writing JSON — and
it removes an entire category of "where does this live / is it stale / which
repo owns it" problems, because **there are no fixtures**.

## Why this works for SWAP (and the one caveat from MODFLOW)

MODFLOW's live comparison is safe because it compares with a **tolerance**
(`htol`), so compiler/FP drift vs. the last release doesn't fail. SWAP's rule is
"byte-identical" — but our regression already aggregates daily output to annual
stats and compares with `TOL = 1e-2` (2 decimals), i.e. it is *already* a
tolerance comparison on aggregates, not a raw byte diff. And our oracle
(`swap420gf`) is a single **deterministic** static binary, so re-running it
reproduces the same aggregates every time. So the live path reuses the exact
same `aggregate()` + `compare()` + `TOL`; only the *source* of the expected side
changes (fresh oracle run instead of stored JSON).

## Design

### Comparison flow (per case, one pytest test)

```
run modern build on toml/        -> (m_annual, m_totals, m_means)   # RuntimeError => pending_restore xfail
run reference on ref.input_variant -> (r_annual, r_totals, r_means) # the live "expected"
compare({years:r_annual, total:r_totals, mean:r_means},
        m_annual, m_totals, m_means)                                # AssertionError => mismatch / known_divergence xfail
```

Run the modern build **first** so a `pending_restore` case (modern fatal-errors
→ `RuntimeError`) short-circuits to xfail before a reference run is spent.

### Pluggable reference (do not hardcode the oracle)

The "expected" side is produced by a selectable **reference descriptor**, not a
hardcoded `swap420gf` path:

```python
class Reference(NamedTuple):
    name: str            # "swap420gf"
    binary: Path         # resolved from swap-testcases/oracle/ (env-overridable)
    input_variant: str   # "legacy" | "toml"  — which case subdir the ref reads
    rc_ok: tuple         # accepted exit codes, e.g. (0, 100)
```

- A small `REFERENCES` registry (like `CASES`), default `swap420gf`
  (`legacy`, rc {0,100}), selected via `SWAP_REGRESSION_REF` env / pytest option.
- `run_reference_and_aggregate(case, ref)` stages `case/<ref.input_variant>/`,
  runs `ref.binary`, accepts `ref.rc_ok`, aggregates with the case's var sets.
- Adding a reference is **config, not code**: drop a binary into
  `swap-testcases/oracle/` + one registry entry.
- Payoff: the *same* machinery gives two regression modes — **modern-vs-4.2.0
  oracle** (default, `legacy` inputs) and eventually **current-vs-last-released
  SWAP** (a future reference reading `toml` inputs, MODFLOW's actual pattern).
  Every case dir already ships both `legacy/` and `toml/`, so both variants
  exist per case.
- Build the seam now; populate only the `swap420gf` entry.

### What changes

- `regression_harness.py` gains the `Reference` type + `REFERENCES` registry and
  `run_reference_and_aggregate(case, ref)` — generalized from
  `regen_reference.py::_run_reference_and_aggregate` (stage `case/<input_variant>/`
  into a temp dir, run `ref.binary`, accept `ref.rc_ok`, aggregate
  `result_output.csv` with the case's var sets). The active reference is resolved
  once (`SWAP_REGRESSION_REF`, default `swap420gf`).
- `test_output_regression.py::test_regression` runs both sides live and compares
  (no `load_fixture`). The `known_divergence` → xfail(`AssertionError`) and
  `pending_restore` → xfail(`RuntimeError`) markers are unchanged and still
  correct: `winter` still diverges live → xfail; `swdrought2` modern still
  fatal-errors → xfail.
- Delete the 20 `*_reference_gf.json` files and the vestigial
  `*_expected_gfortran.json` diagnostics (keep the ifx `*_expected.json`).
- Delete all three helper scripts — `regen_reference.py` (its oracle-run logic
  now lives in the harness), `regen_expected.py` (modern snapshots — redundant
  once comparison is live), and `_switch_validate.py` (a hand-rolled
  reference-vs-modern comparison the harness now does generically).

### What stays

- `CASES` registry, `aggregate()`, `compare()`, `TOL`, the modern
  `run_and_aggregate()` — all unchanged. The registry (var selection,
  `known_divergence`, `pending_restore`) is the modern engine's *judgment* and
  stays in SWAP.
- The pytest/xdist wiring, `fast` marker, pixi tasks — unchanged.
- `swap420gf` + `build_swap420gf.sh` (see oracle-location decision below).

### Oracle / reference location

References live in `swap-testcases/oracle/`, next to the inputs they consume —
"reference dataset + reference engine(s)" in one repo. `swap420gf` +
`build_swap420gf.sh` move there (with the 4.2.0 source provenance for full
self-containment). SWAP resolves `ref.binary` from the pinned `swap-testcases`
checkout, exactly as it already resolves inputs (`SWAP_TESTCASES_PATH` /
`TESTCASES_REF`); `SWAP_REGRESSION_REF` selects which reference. Future
references (e.g. a released modern SWAP) drop into the same `oracle/` dir.

To keep the risky steps separable, Phase 1 lands the live comparison with the
binary still resolvable from `tests/reference/` (via an env default), and
Phase 3 performs the actual relocation + pin bump — but the relocation IS in
scope, not deferred.

## Edge cases & risks

- **Oracle determinism** — `swap420gf` is a fixed binary; we already generated
  fixtures from it deterministically, so live runs reproduce. Low risk. (One
  guard: the aggregation already drops non-numeric fields and rounds to 2 dp, so
  any incidental run metadata in the CSV is ignored.)
- **Oracle fails to run for a registered case** — propagates as `RuntimeError`
  and *fails* the test (correct: it's an infra problem), except on
  `pending_restore` cases, which fail at the modern step first and never reach
  the oracle. No silent masking.
- **`swsalinity=2`** stays a non-case (swap420gf SIGSEGVs on it — no oracle
  possible); unchanged.
- **Compute** — each case now runs two binaries. `check-fast` (4 cases)
  ~2.6s → ~5s; `check-full` (18 cases, xdist-parallel) ~47s → ~80–90s. Within
  the documented budgets. Oracle runs parallelize across cases exactly like the
  modern runs already do.
- **Oracle as hard test-time dependency** — it must be present on every run
  (not just when regenerating). It already is (tracked in `tests/reference`, or
  checked out from `swap-testcases` under option B). MODFLOW pays the same price.

## Phased plan

**Phase 1 — live comparison + pluggable reference (binary still in `tests/reference`).**
1. Add `Reference` + `REFERENCES` (default `swap420gf` → `legacy`, rc {0,100})
   and `run_reference_and_aggregate(case, ref)` to `regression_harness.py`
   (generalize `regen_reference.py`'s oracle runner). Resolve the active
   reference from `SWAP_REGRESSION_REF`; resolve `ref.binary` via an env default
   that currently points at `tests/reference/swap420gf`.
2. Rewrite `test_output_regression.py::test_regression` to run modern +
   reference live and compare; drop `load_fixture`.
3. Session fixture fails fast if the active reference binary is missing (as we
   already do for the modern build).
4. Gate: `pixi run -e test check-full` — expect the same 18 pass / 2 xfail
   (`winter`, `swdrought2`). Confirm no case regresses.

**Phase 2 — retire stored fixtures & all helper scripts.**
5. `git rm` the 20 `*_reference_gf.json` and the `*_expected_gfortran.json`
   diagnostics (keep the ifx `*_expected.json`).
6. Delete `regen_reference.py`, `regen_expected.py`, `_switch_validate.py`.
7. Update docs: CLAUDE.md "Fixtures" bullet (no stored reference; the reference
   is the live comparand), `tests/reference/README.md`, and write **ADR 0055**
   superseding the stored-fixture decision.
8. Gate: `check-full` green; CI green on `development`.

**Phase 3 — relocate the reference into `swap-testcases`.**
9. Move `swap420gf` + `build_swap420gf.sh` (+ 4.2.0 source provenance) into
   `swap-testcases/oracle/`; `ref.binary` resolves from the pinned checkout
   (`SWAP_TESTCASES_PATH`), env-overridable.
10. `swap-testcases` release + `TESTCASES_REF` bump; CI env already provides the
    checkout. `tests/reference/` is then empty of binaries (README may stay as a
    pointer, or the dir is removed).

## Resolved decisions

All prior open questions are settled (see the Locked decisions header):
relocate the reference to `swap-testcases` (Phase 3, in scope); drop all three
helper scripts; make the reference pluggable; keep ifx `*_expected.json`, drop
`*_expected_gfortran.json`; accept ~80–90s `check-full` with no committed cache.

## Rollback

Every phase is a separate reviewed commit on `development`. Phase 1 is additive
(the stored fixtures still exist and could be re-consulted); the point of no
easy return is Phase 2's `git rm`, recoverable from history if needed.
