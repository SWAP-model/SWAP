# Regression fixtures — investigation notes

**Status:** open, bumped to a follow-on spec
**Opened:** 2026-04-22 (Rescue Phase 1, Task 1)

---

## 2026-05-27 — Reference basis switched to gfortran-4.2.0; hysteresis regression found

**Reference basis change.** The regression now compares the modern build against
`tests/reference/swap420gf` — the *unmodified* SWAP 4.2.0 source recompiled with the
modern build's gfortran flags (`-O2 -ffree-line-length-none -std=legacy -finit-local-zero`,
**no source edits**). Previously the harness compared the modern build against its own
golden-master snapshot (`*_expected_gfortran.json`); it now compares against
`*_reference_gf.json` produced by `regen_reference.py`. This makes the suite a *physics*
fidelity check rather than a self-consistency check, because the compiler is held constant.

Verification that compiler is not a confound: on all five pre-existing cases, the Intel
`swap420` and gfortran `swap420gf` builds of 4.2.0 agree to the fixtures' 2-decimal
precision. So any modern-vs-`swap420gf` divergence is a genuine code difference.

**New finding — hysteresis (`SWHYST=1`) regression.** The new `soilhysteresis` case
(clone of hupselbrook, hysteresis on) reveals that the modern hysteresis path diverges
from 4.2.0: up to **1.86 cm GWL** (2003), **1.46 cm** on total DRAINAGE, ~1.3 cm DSTOR.
This is **physics, not compiler**: 4.2.0-ifx ≡ 4.2.0-gfortran (0.000 drift) on this case,
while 4.2.0-gfortran vs modern-gfortran diverges on 18/54 aggregated values.

The `hysteresis` subroutine itself (`src/soil/soilhydraulics.f90`) is a faithful
line-by-line transcription of 4.2.0's `hysteresis.f90`. The case is registered with
`known_divergence=` so the harness reports it as an expected divergence (xfail).

### ROOT CAUSE (traced 2026-05-27) — adaptive-`dt` desync, NOT hysteresis

Full systematic trace (against `swap420gf`, base hupselbrook with daily output):

1. **Not hysteresis.** Driver, `prhead`, `moiscap`, `indeks` init, `hm1` save/reset, and
   the timestep call order are all faithful transcriptions.
2. **Not the compiler.** `swap420` (ifx) ≡ `swap420gf` (gfortran) to 6 decimals every day.
3. **Localized to GWL/DRAINAGE/DSTOR.** In the base (no-hysteresis) case these diverge
   ~1e-5/day from day 1 (max ~0.26 cm daily); `RAIN/EPOT/EACT/INTERC/RUNOFF/QBOTTOM/TPOT`
   match exactly. So the seed is in the saturated-zone / lateral-drainage path, and it
   exists in **every** case — averaged away in the 5 aggregated regression fixtures,
   accumulated to ~1.9 cm only under hysteresis threshold amplification.
4. **Fixed `dt` ⇒ bit-identical.** With `DTMIN=DTMAX` fixed, modern and `swap420gf` are
   bit-identical (modern just has one extra no-op step at the front). The divergence
   appears **only with adaptive `dt`**. So all spatial/flux/hysteresis code is faithful;
   the seed is the adaptive timestep controller.
5. **First-step `numbit` desync.** On the very first timestep (tcum=0, dt=2e-4, identical
   hydrostatic init), `swap420gf` solves in `numbit=4` (gwl moves −75.0→−75.055) while
   modern solves in `numbit=1` (gwl unchanged). The `dt` controller (faithful to 4.2.0)
   doubles `dt` when `numbit≤3`, so modern doubles `dt` prematurely and the two builds
   walk different timestep sequences thereafter → the pervasive ~1e-5 drift.
6. **Init is identical and correct.** The initial `h`-profile is bit-identical between the
   builds AND at exact hydrostatic equilibrium (`h = gwl − z`, deviation 0.0). So the
   `numbit` difference is **not** an init bug; it is a sub-tolerance difference in the
   first solve's residual flipping the `Fmax < CritDevBalCp (1e-6)` convergence threshold.

**Conclusion:** the modern build is faithful to 4.2.0 in all physics and (with fixed `dt`)
bit-identical. The divergence is numerical-threshold sensitivity in the *adaptive-`dt`
first-step convergence*, amplified only under hysteresis. Fixing it to match 4.2.0 would
mean perturbing the core solver to replicate a sub-1e-6 4.2.0 behaviour — high risk to the
5 passing cases for a <2 cm effect. Decision pending; documented as xfail meanwhile.

**New finding — frost-path drift (`8.winter`, SWSNOW+SWFROST+SWSUBLIM).** The winter case
shows the **snow path reproduces 4.2.0 exactly** (SNOW peak 0.574 cm in both builds), but
the **frost path (SWFROST=1) drifts minutely**: 2003 DRAINAGE/RUNOFF differ by ~0.03–0.04
cm (total DRAINAGE 74.10 vs 74.14). Much smaller than hysteresis, but real (base case has
0.000 drift). Suspected in the frost soil-water-flow reduction factor. Registered xfail.

Note on output columns: the modern build sources its CSV column list from TOML
`[output.csv] inlist`, **not** from the staged `swap.swp` `INLIST_CSV` (which only
`swap420gf` reads). New-physics columns (e.g. `snow`) must be added to BOTH for an
apples-to-apples comparison.

---

## What changed

The rescue committed to **gfortran-only** during Phases 1–4 (see `docs/adr/0001-gfortran-first.md` — written in Task 8 of the Phase 1 plan). Phase 0's baseline regression numbers were recorded under **ifx** (Intel), because `pixi.toml` silently hardcoded `FC=ifx` in the production configure tasks and `meson.build` had an Intel-specific flag path including `-init=zero`.

When the compiler was swapped to gfortran (plus `-finit-local-zero` to match Intel's zero-initialization semantics), five of six cases matched the existing ifx-based fixtures cleanly. Two cases diverged:

- **`macropore`**: `DRAINAGE` drifts by ~21 at year 1998 (out of ~100), ~3 at year 1999, with corresponding small `GWL` deltas (0.3–2.0). This is large enough to be a real behavioral difference, not floating-point noise.
- **`oxygenstress`**: the pre-existing `MOWDM` deviation (max 85 in year 1995) was already visible under ifx; same magnitude and same year under gfortran — so this deviation is NOT compiler-driven. It is a pre-existing physics or fixture issue that has been tolerated from well before the rescue.

## What we did

Rather than picking a side, we kept both reference sets:

- `*_expected.json` — historical ifx-produced values. **Unchanged.** These document what the upstream-Intel-compiled reference produced.
- `*_expected_gfortran.json` — regenerated from the gfortran+finit-local-zero build on 2026-04-22. These are what the regression harness actually compares against under the current compiler policy.

The harness `CASES` dict points at the `_gfortran` files only. A future compiler-policy change would add a selector — not in scope for the rescue.

## What needs investigation (not part of the rescue)

The open question is whether the **macropore DRAINAGE divergence** is a physics problem or a compiler-flag artefact. Specifically:

1. Does gfortran `-finit-local-zero` cover every initialization path that ifx's `-init=zero` covered? (`SAVE`d module variables, allocated arrays, derived-type components, and COMMON blocks are all separately-controllable dimensions.)
2. Do additional FP-model flags (`-ffp-contract=off`, `-fno-unsafe-math-optimizations`, `-fno-fast-math`) move gfortran's output closer to ifx's? If so, the divergence is numerical compiler choice. If not, something in the macropore code is genuinely compiler-sensitive (stale pointer, undefined-order-of-evaluation arithmetic, etc.).
3. For oxygenstress MOWDM: is the 1995-spike a crop-growth/mowing-schedule physics bug that the fixture has always masked, or a numerical artefact of a scheduled event? Since gfortran and ifx both hit the same number, this is almost certainly physics, not compiler.

None of these block the rescue. They are tracked here as open items to revisit during Phase 4 module cleanup or as a dedicated physics-audit follow-on spec.

## Reproducing the ifx reference (if needed for an investigation)

The legacy SWAP 4.2.0 Intel-compiled Linux binary is preserved at `tests/reference/swap420`. Running it via `pixi run swap-ref` generates output comparable to the historical `*_expected.json` fixtures (modulo the MOWDM deviation the fixtures themselves encode).

---

## 2026-05-31 — cropfixed (type-1) switch re-enablement: two divergences surfaced

Re-enabling the switches that legacy `readcropfixed` supports but the modern
cropfixed TOML path stub-errored (ADR 0015 shortcuts). Most were over-broad
guards over intact compute: **swrd=2, swharv=1, swcompensate=1/2, swcf=3,
swinter=2** all reproduce `swap420gf` byte-for-byte on a maizes (type-1)
variant and were shipped. Two switches surfaced genuine modern-vs-4.2.0
divergences and were **kept gated / flagged**:

### swsalinity=1 (Maas-Hoffman) — kept gated

The reduction kernel (`rootextraction.f90`, `if swsalinity==1 ...`) is identical
to legacy and the config/state plumbing already exists, but enabling it for a
simple crop diverges from `swap420gf`: on hupselbrook maize (SWSOLU=1,
saltmax=3.0, saltslope=0.1) **TACT differs ~3.4 cm/yr**, with matching shifts in
DRAINAGE/GWL.

Diagnostic: with **saltslope=0** (the branch executes but produces zero
reduction) modern ≡ legacy byte-for-byte. So enabling the switch does **not**
perturb init — the divergence appears only when the reduction is non-zero.
Salinity is the only stress whose magnitude reads the solute concentration
`sol%cml`, so it closes a feedback loop *salinity → reduced uptake → cml →
salinity*. Drought/oxygen depend on pressure head, not solute — which is why the
baseline (those stresses + solute, salinity off) matches but this does not. The
likely cause is a timestep-ordering difference in when `cml` is refreshed
relative to the salinity evaluation, introduced in the solute refactor and
invisible until `swsalinity=1` (the only in-`rootextraction` consumer of `cml`).
Note the WOFOST `salinitystress` regression case passes — so either its
concentrations stay below threshold or its orchestration path refreshes `cml`
differently from the simple-crop path; not yet isolated.

### swcf=1/3 + swetr=0 (crop-factor ET under Penman-Monteith) — pre-existing 0.01 cm GWL

swcf=3 (wet-crop factor) was shipped: byte-identical on every water flux and all
ET terms vs `swap420gf`. The lone difference is the **annual-average GWL, off by
one output unit (0.01 cm)**. This is **pre-existing and not introduced by the
wet-crop work**: the already-allowed **swcf=1** path shows the identical 0.01 cm
GWL delta on the same maize case. It is a floating-point rounding property of the
crop-factor ET path combined with Penman-Monteith (swetr=0) — a combination no
regression case exercises (cases pair swcf=2 with swetr=0, or swcf=1 with
swetr=1). All water-balance fluxes are byte-identical; only the GWL daily
interpolation rounds differently. Tracked as an open FP-sensitivity item, not a
blocker.

Validation tool: `tests/regression/_switch_validate.py` runs a maizes (type-1)
variant through both `swap420gf` (legacy ASCII) and the modern build (TOML) per
switch and compares the harness's aggregated flux/state vars.
