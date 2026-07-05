# ADR 0055 — Live double-run regression (no stored fixtures)

**Status:** Accepted — landed 2026-07-05
**Arc:** Test-suite review (follows ADR 0054)
**Supersedes:** the *stored-fixture* mechanism described in ADR 0054 (the
`*_reference_gf.json` golden files) and the "fixtures stay in SWAP" decision
recorded in the repo-ecosystem plan.
**Design doc:** `dev-docs/2026-07-05-live-double-run-regression-design.md`.

## Context

The byte-identical regression compared the modern build's aggregated output
against stored `*_reference_gf.json` fixtures produced by the `swap420gf` 4.2.0
oracle. Once the case *inputs* moved to the `swap-testcases` sibling repo
(ADR 0054), the fixtures had no coherent home: a reference fixture is a function
of *(oracle binary × case inputs)*, the inputs live in `swap-testcases`, and the
oracle-derived outputs are neither an assertion about the modern engine (so the
"keep them with the bisectable modern code" rationale doesn't apply — they never
change when modern code changes) nor a natural resident of a frozen 4.2.0
source tag (which can't host a fixture set that grows with every new case).

MODFLOW dissolves the analogous problem by never storing expected outputs: it
runs a reference build and the current build and compares live (tolerance-based,
`htol`).

## Decision

Replace stored fixtures with a **live double-run**. For each case the suite runs
two engines in isolated temp dirs and compares their annual aggregates within
the existing `TOL`:

- the **modern** build on the TOML inputs, and
- a **reference** engine on its input variant.

The reference is **pluggable** — a `Reference` descriptor
`{name, binary, input_variant, rc_ok}` in a `REFERENCES` registry, selected by
`SWAP_REGRESSION_REF` (default `swap420gf`, reading `legacy` inputs). Adding a
reference is config, not code; in particular a future *released SWAP* can serve
as the reference (reading `toml` inputs) — MODFLOW's current-vs-last-release
pattern — with one registry entry. The comparison reuses the unchanged
`aggregate()` / `compare()` / `TOL`; only the *source* of the expected side
changes (a fresh reference run instead of a JSON load). This works under SWAP's
byte-identity rule because the oracle is deterministic and the comparison is
already a 2-decimal aggregate tolerance, not a raw byte diff.

Consequently:

- All `*_reference_gf.json` (and the vestigial `*_expected_gfortran.json`
  diagnostics) are deleted. The historical ifx `*_expected.json` snapshots are
  kept as provenance only (never compared against).
- The three fixture helper scripts (`regen_reference.py`, `regen_expected.py`,
  `_switch_validate.py`) are deleted — their logic is either absorbed into the
  harness or made redundant by live comparison.
- The reference binary is **not committed anywhere**. It is published as a
  release asset by a new archival repo, **`SWAP-model/swap-4.2.0`** (the pristine
  4.2.0 source + a multi-compiler Meson build; the native-gfortran build
  reproduces the previous `swap420gf` byte-for-byte), and downloaded on demand
  into a gitignored cache — mirroring how the 4.2.0 build already pulls
  `libttutil` from the `SWAP-model/ttutil` releases. The `Reference` descriptor
  carries `{repo, tag, asset}`; the release tag is the pin. `tests/reference/`
  is deleted.

## Consequences

- **The fixture-location problem is gone** — nothing to store, locate, version,
  or assign to a repo. The three repos separate cleanly: `swap-4.2.0` owns the
  *reference engine* (source + released binaries), `swap-testcases` owns the
  *data* (case inputs), and SWAP owns the *judgment* (the `CASES` registry:
  which vars to compare, `known_divergence`, `pending_restore`).
- **Cost:** each case runs two binaries. `check-fast` ~2.6s→~5s; `check-full`
  ~47s→~60s (xdist absorbs most of the doubling). The oracle is a ~7 MB
  download, cached after first use (CI can cache the dir).
- **The reference is a hard test-time dependency** — downloaded from the pinned
  `swap-4.2.0` release every fresh run (concurrency-safe under xdist via a
  lock), as it is for MODFLOW. `SWAP_REFERENCE_BIN` overrides with a local path
  for offline/dev work; a download failure fails the session fast with a clear
  message.
- **No "regenerate fixtures" workflow** — a physics change either matches the
  reference or fails; you fix code, never a golden file. Divergences are `xfail`
  markers on the registry.
- Verified: full suite 18 passed / 2 xfailed (`winter`, `swdrought2`) —
  identical to the stored-fixture result.
