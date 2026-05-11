# Known Issue: SS-CRP Phase 0 — swdrought=2 / swoxygen=2 / swcompensate / swfrost Coverage Gap

**Date:** 2026-05-10
**Status:** known-issue
**Related arc:** SS-CRP (crop water uptake state migration)
**Filed at:** C-0.1 (read-only audit; no code changes)

---

## Context

The crop water uptake arc (`rootextraction.f90`) carves 22 fields into
`state%soilwater`. The JvL path (`swdrought=2`, lines 313–760) and oxygen
stress path (`swoxygen=2`, `oxygenstress.f90`) are large; zero check-full
regression coverage of these paths is a known risk.

Four switches gate the code paths most relevant to the 22 migrated fields:

| Switch | Values in scope | Physical meaning |
|---|---|---|
| `swdrought` | 1=Feddes (1978), **2=De Jong van Lier** | Drought-stress model. `swdrought=2` activates the entire JvL compute path: `JongvanLier`, `JongvanLierLoop`, `MatricFlux` Task 2, and the `mfluxtable` / `hroot` / `hleaf` / `mflux` / `mroot` / `rootrho` / `rootphi` / `rmax` / `alpJvLier` / `Tactual` / `Hxylem` fields — 14 of the 22 owned fields are `swdrought=2`-only. |
| `swoxygen` | 0=none, 1=Feddes, **2=Bartholomeus** | Oxygen-stress model. `swoxygen=2` activates `OxygenStress(node, rwu_factor, state)` with the full Bartholomeus physical sub-model — a separate 1800-line computation in `oxygenstress.f90`. `swoxygen=1` uses Feddes reduction factors already embedded in `rootextraction.f90`. |
| `swcompensate` | 0=none, **1=Jarvis** (1989), 2=Walsum (stub) | Compensation for reduced uptake. `swcompensate=1` activates a re-distribution loop (rootextraction lines ~230–280) that can increase `qrot` in non-stressed nodes, potentially redistributing all primary arrays (`qrot`, `qredwet`, `qreddry`, `qredsol`, `qredfrs`). |
| `swfrost` | 0=none, **1=active** | Frost reduction. Gates `alpfrs = 0.0d0` when `state%heat%tsoil(node) < 0.0` (rootextraction line 193). Writes into `qredfrs` / `qredfrssum`. Depends on `state%heat%tsoil` already plumbed from SS-HEAT Task 6. |

---

## check-full Coverage (Audited C-0.1)

Switch values from `tests/swap-cases/toml/<case>/*.crp.toml`.
`swfrost` is in `swap.toml` (soil section); none of the 6 cases set it
(Fortran default = 0, confirmed in `soil_config.f90:56`).

| Case | swdrought | swoxygen | swcompensate | swfrost |
|---|---|---|---|---|
| 1.hupselbrook | 1 (Feddes) | 1 (Feddes) | **1 (Jarvis)** on grassd; 0 on potatod/maizes | 0 (default) |
| 2.grassgrowth | 1 (Feddes) | 1 (Feddes) | 0 | 0 (default) |
| 3.macroporeflow | 0 (default) | 0 (default) | 0 (default) | 0 (default) |
| 4.oxygenstress | 1 (Feddes) | **2 (Bartholomeus)** | **1 (Jarvis)** | 0 (default) |
| 5.salinitystress | 1 (Feddes) | 1 (Feddes) | 0 | 0 (default) |
| 6.surfacewater | 1 (Feddes) | 0 (none) | 0 | 0 (default) |

**Summary:**

- `swdrought=2` (JvL): **zero check-full coverage** across all 6 cases.
- `swoxygen=2` (Bartholomeus): covered by **4.oxygenstress** only.
- `swcompensate=1` (Jarvis): covered by **1.hupselbrook (grassd)** and **4.oxygenstress**.
- `swfrost=1`: **zero check-full coverage** across all 6 cases.

---

## Stub-Error Status

The config validators enforce the following at parse time (as of this audit):

| Switch + value | cropfixed | cropgrass | cropwofost |
|---|---|---|---|
| `swdrought=2` | stub-errored | stub-errored | stub-errored |
| `swoxygen=2` | stub-errored | supported (swoxygentype=2 stub) | stub-errored |
| `swcompensate=2` (Walsum) | (no check found) | stub-errored | (stub noted in code) |

`swdrought=2` is **globally stub-errored** across all three crop types.
This means no TOML simulation can currently reach the JvL compute path
at runtime — the config validator will reject the input before execution.

---

## Risk Assessment

**swdrought=2 (JvL path): LOW runtime risk, MODERATE migration risk.**
Because `swdrought=2` is stub-errored in all three crop configs, no
regression case can exercise this path. The 14 JvL-gated fields
(`mfluxtable`, `hroot`, `hleaf`, `mflux`, `mroot`, `rootrho`, `rootphi`,
`rmax`, `Hxylem`, `alpJvLier`, `Tactual`, and the 3 init co-writes in
`cropgrowth.f90`) will be dual-written and cutover without any integration
check. A silent semantic error in the JvL dual-write (C-1.4) or in the
`MatricFlux(1)` relocation (C-1.2) would not be caught by check-full.
Manual inspection of write sites is the only safeguard during this arc.

**swoxygen=2 (Bartholomeus): LOW risk.** Covered by 4.oxygenstress.
`OxygenStress(node, rwu_factor, state)` is already state-plumbed from
SS-HEAT; no fields owned by this arc are written inside `oxygenstress.f90`.

**swcompensate=1 (Jarvis): LOW risk.** Covered by both 1.hupselbrook
(grassd) and 4.oxygenstress. The Jarvis compensation re-distribution loop
exercises `qrot`, `qredwet`, `qreddry`, `qredsol`, `qredfrs` at lines
~230–295; these are primary-path fields covered by every case with crop.

**swfrost=1: LOW runtime risk, NEGLIGIBLE migration risk.**
`swfrost=1` gates only a single scalar zero-assignment (`alpfrs = 0.0d0`,
line 193) and propagates into `qredfrs` / `qredfrssum`. These are primary
path fields written regardless. The frost path itself only suppresses uptake;
no JvL-exclusive fields are touched. Runtime risk is low.

**Overall regression posture for this arc:**
The byte-identical check-full gate protects all primary-path fields through
5 of 6 cases. The JvL sub-path (14 fields, `swdrought=2`) has no integration
coverage but is blocked from production use by stub-errors. Migration
correctness relies on code-review and manual inspection of C-1.2 / C-1.3 /
C-1.4 write sites.

---

## Recommended Fixture Authoring

Out of scope for this arc. Deferred per plan C-0.1 (analogous to the
2026-05-08 nutrient-regression-case deferral and the 2026-05-10 boundary
swbotb=2/4/8 deferral). Prerequisite: stub-errors for `swdrought=2` must be
lifted in all three crop configs before a regression fixture can be authored.

Suggested future cases:

- `tests/swap-cases/toml/7.jvldrought/` — `swdrought=2`, cropgrass or cropwofost
- `tests/swap-cases/toml/8.frostuptake/` — `swfrost=1`, any crop type

Each requires reference output from a known-good binary (non-TOML legacy
run) before inclusion in the check-full harness.

---

## Cross-References

| Artifact | Path |
|---|---|
| Discovery | `docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-discovery.md` |
| Design | `docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md` |
| Plan | `docs/superpowers/plans/2026-05-10-crop-uptake-state-migration.md` |
| Boundary gap (analogous) | `docs/superpowers/specs/2026-05-10-boundary-phase0-coverage-gap.md` |
| Successor ADR | `docs/adr/0036-state-migration-crop-uptake.md` (drafted at C-2.6) |
