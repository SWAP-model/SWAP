# Regression fixtures — investigation notes

**Status:** open, bumped to a follow-on spec
**Opened:** 2026-04-22 (Rescue Phase 1, Task 1)

---

## 2026-05-27 — Reference basis switched to gfortran-4.2.0; hysteresis regression found

**Reference basis change.** The regression compares the modern build against
`swap420gf` — the *unmodified* SWAP 4.2.0 source recompiled with the
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

The legacy SWAP 4.2.0 Intel-compiled Linux binary (`swap420`) is published as a release asset of `SWAP-model/swap-4.2.0`. Running it on a case's legacy inputs generates output comparable to the historical `*_expected.json` fixtures (modulo the MOWDM deviation the fixtures themselves encode).

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

### 2026-05-31 addendum — grass/wofost restoration pass

The same gated branches exist in cropgrass and cropwofost; restored the ones with
intact shared compute, all byte-identical vs swap420gf on potatod/grassd:
- **cropwofost**: swharv=1, swcompensate=1/2 (init was missing the alphacrit/
  dcritrtz copies), swinter=2 (split the schema gashtb into the atmosphere arrays).
- **cropgrass**: swcompensate=2 (Walsum; Jarvis=1 already worked), swinter=2
  (added a gashtb schema field + reader + atmosphere split).

Still gated, consistent with the cropfixed findings:
- **swsalinity** (grass): left gated — same salinity→cml feedback divergence as
  cropfixed (wofost swsalinity=1 already works, so it is simple/grass-path
  specific; not isolated). swsalinity=2 stays gated everywhere (needs swdrought=2).
- **swcf=3** (grass stub-errored; wofost enum-excluded to [1,2] though legacy
  reads swcf=1..3): NOT restored in this pass. It is the same wet-crop-factor
  (cfeic) branch as cropfixed swcf=3 and so inherits the pre-existing crop-factor
  + Penman-Monteith 0.01 cm GWL artifact; it also needs per-module cfeictb
  plumbing. RESTORED for both (per-module cfeictb plumbing; the LAI-indexed
  runtime cfeic lookup was already intact). Both validated byte-identical incl.
  GWL on grassd/potatod — so the cropfixed maize 0.01 cm GWL artifact is
  case-specific, NOT inherent to swcf=3.

Net after this pass: the only remaining gated branch with intact-but-divergent
compute is swsalinity=1 (cropfixed + grass) — the cml-feedback issue above.
Everything else still gated has genuinely deleted/dormant compute (Tier C/D):
swdrought=2, swinter=3, swoxygen=2 inputs/repro, swsalinity=2.

---

## 2026-05-31 — Tier C/D (deleted/dormant compute) restoration scoping

Investigation of the four remaining gated branches whose compute was deleted or
made dormant (vs. the Tier-A/B switches, which only needed re-wiring). All
recovery paths verified. **Common blocker: no TOML regression case exercises any
of these, so each restoration needs a bespoke byte-identical case** (and
swoxygen=2 needs numerical heat, swcalt=2).

### Cluster 1 — swoxygen=2 (Bartholomeus) — LEAST work, kernel proven
- `OxygenStress` kernel (src/crop/oxygenstress.f90) is intact AND compiled AND
  already regression-tested: the `4.oxygenstress` case runs **grass swoxygen=2 /
  swoxygentype=1 byte-identical** (cropgrass only gates swoxygentype==2).
- Kernel branches on croptype (oxygenstress.f90:198,241): croptype 2/3 compute
  `max_resp_factor` internally via GET_MAX_RESP_FACTOR (no input tables); croptype 1
  needs `w_root_ss` from a `wrtb` table.
- → **wofost swoxygen=2/type1 = guard removal** (validator cropwofost_config:781,
  init cropwofost_init:125). **cropfixed swoxygen=2/type1 = add wrtb/mrftb static
  tables** (hard-coded 0 at cropfixed_runtime:135,139,216,220) + guard removal.
- `swoxygentype=2` (reproduction fns): dormant `src/crop/dormant/oxygenrepro.f90`
  (not in build), dispatch stub at rootextraction.f90:158. Needs build inclusion +
  OxygenSlope/OxygenIntercept state + oxygen_dat wiring. No reference case exists.

### Cluster 2 — swinter=3 (adapted-Rutter storage interception)
- `msw1eic` + `ruttervw` deleted in **f653aed** ("drop dead msw1eic Gash kernel");
  recover via `git show f653aed^:src/atmosphere/interception.f90`.
- State (`sicact`/`siccapact`/`fimin`) still on atmosphere_state. Runtime stubs
  hard-code `siccapact=0` (cropfixed_runtime:130,211 + wofost/grass analogues);
  DivIntercep guarded with `swinter.ne.3` (meteo_orchestrator:326); wet-fraction
  branches reference an `eintc` that's never set. CAVEAT: legacy msw1eic used
  real(4) (MetaSWAP interop) — must port to real64 carefully for byte-identity.
- Shared in meteo_orchestrator (crop-unaware, runs for active crop); the per-crop
  runtimes only set siccapact = siccaplai*lai.

### Cluster 3 — swdrought=2 (De Jong van Lier) + swsalinity=2 — MOST work
- `jongvanlier.f90` (JongvanLier + JongvanLierLoop, nested Newton-Raphson on
  hleaf/Hxylem) deleted in **5c82f0a**; recover via
  `git show 5c82f0a^:src/crop/dormant/jongvanlier.f90` (already modernized to
  state%, but imports ~9 bare globals: kroot/rxylem/kstem/rootradius/rootcoefa/
  rooteff/wiltpoint/stephr/criterhr — must finish migrating to state).
- `MatricFlux`/`matric_flux`/`matricflux_build_table` are STILL LIVE in
  rootextraction.f90 (called from cropgrowth.f90 when swdrought=2); the swsalinity=2
  osmotic-head correction (`hosm = salthead*cml`) lives inside matric_flux and is
  complete. So **swsalinity=2 comes for free with swdrought=2** (only reachable
  there). State fields (mflux/mroot/hroot/mfluxtable/Tactual/alpJvLier/twilt) exist.
- Config params declared as stubs (cropfixed_config:111-127) but not TOML-read.
- Highest byte-identical risk: Newton-Raphson iteration count is tolerance-sensitive
  (CriterHr/StepHr).

**Recommended order: 1 → 2 → 3** (ascending effort/risk). Each needs its own
heat/feature-enabled regression case authored or converted from a legacy ASCII
crop (legacy/swap-4.2.0:xdata/crops/* use SWDROUGHT=2/SWINTER=3).

### 2026-06-01 — wofost swoxygen=2 (Bartholomeus type-1): gated, 0.01 cm TACT divergence

Cluster-1 attempt. The plumbing was completed and VERIFIED (swoxygentype field +
reader + the swoxygen==2 init block mirroring the proven cropgrass path; every
kernel input wired: q10, rmr→c_mroot, rfsetb→f_senes, q10_microbial,
specific_resp_humus, srl, swrootradius, root_radiusO2). Yet potatod flipped to
swoxygen=2 diverges from swap420gf by **0.01 cm on annual TACT** (two years),
marginally over the 1e-2 harness tolerance; GWL and all other vars match.

Decisive diagnostics:
- A full all-columns diff of the UNMODIFIED hupselbrook run (legacy vs modern) is
  byte-identical — so potatod biomass is NOT pre-diverged; the 0.01 appears ONLY
  when swoxygen=2 is enabled.
- The grass (type-3) oxygen path IS byte-identical (4.oxygenstress passes). So the
  kernel is exact for type 3 but the **type-2 (wofost) path was never validated**.

→ This is the same class as swsalinity=1: a feedback-coupled stress (oxygen →
reduced uptake → biomass → max_resp_factor → oxygen) amplifying a sub-threshold
FP difference in the type-2 kernel branch. REVERTED; wofost swoxygen=2 stays
gated pending root-cause. Note the non-feedback restorations (swcompensate,
swinter=2, swcf=3) were all byte-identical — the divergences cluster on the
feedback-coupled stresses (salinity, oxygen, and — high risk — the JvL drought
Newton-Raphson still to come).

### 2026-06-01 — permanent crop-switch regression cases + swsalinity=2 SIGSEGV

Added standalone regression cases (hupselbrook + one setting, legacy+TOML, under
tests/regression/cases/) with swap420gf fixtures, registered via the new local=
/pending_restore= CaseConfig fields. PASS: swrd2, swharv1, swcompensate1/2,
swinter2, swcf3 (wofost+grass). xfail: swcf3_maize (known 0.01cm GWL), and
restoration targets swsalinity1, swoxygen2, swinter3, swdrought2 (pending_restore).

**swsalinity=2 (osmotic head): no case — swap420gf itself SIGSEGVs.** Authoring a
maizes swdrought=2 + swsalinity=2 case, the LEGACY reference crashes in
`matricflux_`/`jongvanlier_` (SIGSEGV) on the hupselbrook config. So 4.2.0 itself
cannot run osmotic-head salinity here and no oracle fixture can be produced —
swsalinity=2 is omitted from the suite. (swdrought=2 alone runs fine.)

### 2026-06-01 — swinter=3 (adapted-Rutter) RESTORED (Tier-C cluster 2)

Recovered the deleted `msw1eic` + `ruttervw` kernel (verbatim from f653aed^, real(4)
MetaSWAP ODE kept for byte-identity) into src/atmosphere/interception.f90, restored
the DAILY swinter=3 branch in meteo_orchestrator (Section 5: gctp, ruttervw,
DivIntercep), and re-wired the inputs: cropfixed_config reads fimin/siccaplai,
cropfixed_init plumbs fimin→atmosphere + siccaplai→crop state, cropfixed_runtime
computes siccapact = siccaplai*lai (was hard-zeroed). swinter=3 validator/init
guards removed.

**Faithful — bit-identical to 4.2.0 under fixed dt** (the /tmp fixed-dt probe: all
output columns, worst diff 0.0). Under adaptive dt it diverges ~2.3 cm INTERC/GWL on
the maize year, because the Rutter kernel integrates its canopy ODE over `dt`
DIRECTLY (legacy ruttervw uses the global adaptive `dt`, not the daily `dttp=1`), so
it is acutely sensitive to the SAME adaptive-dt first-step desync that drives the
hysteresis/winter xfails. Registered as known_divergence, not pending_restore.

Scope: only the DAILY meteo path is restored (hupselbrook is SWMETDETAIL=0). The
sub-daily swinter=3 branch (siccaptb afgen) in the SWMETDETAIL=1 orchestrator is
still stubbed — a follow-up when a sub-daily case needs it.

### 2026-06-01 — swdrought=2 (De Jong van Lier): migration done, runtime HANGS (WIP)

Attempted Tier-C restoration. The kernel migration is COMPLETE and COMPILES:
recovered jongvanlier.f90 (JongvanLier + JongvanLierLoop) from 5c82f0a^, migrated
its `use variables` reads (kroot/kstem/rxylem/rootradius/rootcoefa/rooteff/stephr/
criterhr/wiltpoint) to new state%crop%common scalars + the retired
`state%cfg%simulation%numerical%taccur` to state%crop%common%taccur, and MERGED the
two subroutines into rootextraction_mod (siblings of matric_flux — avoids the
jongvanlier↔rootextraction module cycle). Wired config reader + validation + init
plumbing; dispatch `call JongvanLier(state)`; matricflux_build_table already live.

BLOCKER: the modern run HANGS (>200s; legacy swap420gf runs fast). Signature is
dt-collapse — the restored Newton-Raphson produces wrong transpiration → soil
solver fails → adaptive dt shrinks toward zero → effectively unbounded timesteps.
Two contributing facts: (a) the counter>1000 caps use fatalerr_collected (collects
+ continues) not legacy's fatalerr (stops) — added `exit` but still hangs, so it's
the dt-collapse not those loops; (b) a value bug in the migration makes JvL diverge
from legacy. Needs side-by-side instrumented debugging (compare JvL hleaf/Tactual/
mflux per call legacy-vs-modern) — a dedicated session.

REVERTED to keep the suite green; swdrought2 stays a pending_restore target. The
full (compiling) migration is preserved at dev-docs/wip/swdrought2-jvl-restoration.patch
(git apply to resume). NB: swap420gf itself SIGSEGVs on swsalinity=2 (the JvL
osmotic-head sibling), so that path is doubly fragile.

### 2026-06-01 (cont.) — swdrought=2: root-caused the hang (wiltpoint bug), now a PERF wall

Instrumented the hang. **Two findings:**

1. **REAL BUG FOUND (the hypothesis was right).** `matricflux_build_table` /
   `matric_flux` / the runtime `twilt` used `state%crop%common%hlim4` as a stand-in
   for `wiltpoint` — a migration shortcut whose comment said *"wiltpoint legacy
   global is always 0.0 on TOML"*. Once wiltpoint IS plumbed (-20000 vs hlim4's
   -8000), the matric-flux table was cut off at the wrong pressure head, so the JvL
   Newton-Raphson got inconsistent matric flux in the -8000..-20000 range and never
   converged → dt held normal but the solver spun. Fixed: use
   `state%crop%common%wiltpoint` at all three swdrought=2 sites (NOT the Feddes
   swdrought=1 hlim4 uses). Debug confirmed: wiltpoint=-20000, kstem/kroot/taccur all
   correct; JvL then converges to sane values (alp=1.0, hleaf=-54, qrosum=ptra).

2. **NEW BLOCKER — ~1000x perf regression.** With convergence fixed, the sim
   PROGRESSES correctly (dt normal ~0.02-0.04, sane JvL output) but is
   catastrophically slow: ~10s/sim-day → would take hours. **Legacy swap420gf runs
   the same case in 0.70 s.** So modern's JvL path is ~1000x slower — state-record
   indirection (state%soilwater%mfluxtable(lay,count), state%mesh%layer(node), …) in
   the matric_flux hot loop (called ~thousands of times/JvL call) vs legacy's flat
   module globals, likely amplified by an inner-loop convergence-iteration
   difference. Needs profiling + caching the hot state fields into locals — a
   distinct optimization task.

Status: the swdrought=2 restoration is now CORRECT but perf-blocked. Reverted to
keep the suite green (swdrought2 stays pending_restore). Updated patch (incl. the
wiltpoint fix) at dev-docs/wip/swdrought2-jvl-restoration.patch.

### 2026-06-01 (cont.) — swdrought=2 perf: it's the Richards solver, NOT the JvL

Profiled the ~1000x slowdown (perf unavailable — used cpu_time + stub isolation):

- **The JvL kernel is CHEAP.** Direct cpu_time: 137 JvL calls = 182 us total
  (~1 us/call); only ~15-21 matric_flux calls/JvL-call; matricflux_build_table runs
  exactly ONCE. So the JvL is NOT the cost.
- **The cost is the per-timestep Richards solve.** Each outer step is ~0.146 s
  wall vs ~us without swdrought=2, but <1 us of that is JvL. An early-return JvL
  stub made the sim finish in 0.45 s — BUT that stub also zeroed soil%qrot, so it
  removed the solver's *response* to the JvL qrot, not the JvL cost.
- **Mechanism.** RootExtraction (→JvL) is called ONCE per outer step (same as
  legacy, swap.f90:209), computing qrot held fixed through the Richards solve. The
  modern `headcalc` does not converge well with the JvL qrot → trips the
  `time%fldtreduce` dt-reduction loop (swap_mod.f90) → many soilwater_step retries
  per outer step → slow. dt itself reads normal (0.02-0.04) at the outer level; the
  spin is in the inner reduce/retry. Legacy runs the identical case in 0.70 s, so
  legacy's headcalc converges with the same qrot.
- **Conclusion.** This is the SAME adaptive-dt / Richards convergence-threshold
  sensitivity that drives the hysteresis/winter known_divergences — here amplified
  by the stiff JvL coupling into a *performance* spiral. NOT a leftover sync /
  dual-write (checked JongvanLier/JongvanLierLoop/matric_flux — clean Newton-Raphson,
  no array copies). Fixing it means touching core headcalc convergence (high risk to
  the 5 passing cases), so it is out of scope for a mechanical restoration.

Net: swdrought=2 is restored + correct (JvL converges to sane values after the
wiltpoint fix) but perf-blocked by the Richards-convergence sensitivity. Patch
(incl. wiltpoint fix) preserved at dev-docs/wip/swdrought2-jvl-restoration.patch.

### 2026-06-01 (cont.) — swdrought=2 perf ROOT CAUSE: solute sub-stepping, not the JvL/Richards

Full profiling chain (cpu_time + marker isolation; perf unavailable). Each layer
ruled out the previous suspect:
- JvL kernel: CHEAP (137 calls = 182 us cpu, ~1us/call). NOT it.
- headcalc (Richards): converges in numbit=2, called ~once/step. NOT it.
- soilwater_step / dt-reduction loop: fldecdt never fires; <2000 calls. NOT it.
- **solute_step: THE CULPRIT.** Markers around the post-Richards components showed
  solute_step ENTERED 3350x but EXITED 3349x — the sim hangs INSIDE solute_step.

Mechanism: solute transport sub-steps with a von Neumann dispersion stability limit
(solute.f90 prepare_solute_dispersion): `dtsolu = min(dt, dz^2/(2*dispr))`,
`dispr = diffus + ldis*|q|/theta`, then `do while (dt-tcumsol>1e-8)` with
`dtsolu = max(dtsolu, dtmin)`. Under swdrought=2 the JvL extracts water toward the
(correct, post-fix) wiltpoint=-20000, so a node's theta/soil flux q drives `dispr`
large → `dtsolu` collapses to dtmin → the loop runs up to dt/dtmin (~40000)
sub-steps PER timestep → the ~1000x.

Legacy comparison: the legacy solute loop is BYTE-FOR-BYTE the same formula+clamp
(solute.f90:110-121) yet runs the whole case in 0.70s — so legacy's `dispr` does
NOT collapse. The difference is modern's soil water `q`/`theta` at the critical
node vs 4.2.0 — i.e. the SAME modern-vs-4.2.0 soil-water numerical (adaptive-dt /
Richards convergence-threshold) difference behind the hysteresis/winter
known_divergences, here amplified through the solute stability criterion into a
performance collapse rather than a value drift.

CONCLUSION: NOT a leftover sync/dual-write (JvL/matric_flux/headcalc/solute loops
are clean and identical to legacy). The fix is the core soil-water parity work
(make modern q/theta match 4.2.0) — high risk to the 5 passing cases, out of scope.
swdrought=2 stays restored-but-perf-blocked; patch at
dev-docs/wip/swdrought2-jvl-restoration.patch.

---

## 2026-06-11 — ROOT CAUSE FOUND & FIXED: adaptive-dt desync = drainage first-step no-op

The "adaptive-dt threshold desync" behind the hysteresis / winter / swinter3 /
swcf3_maize known-divergences (and the swdrought2 perf collapse) is **NOT** an
irreducible sub-1e-6 convergence-threshold sensitivity. It is a concrete
**ordering bug in `src/drainage/drainage.f90`**, now fixed.

### Mechanism

`drainage()` is called once per timestep from `swap_run_step`. It contained:

```
if (surf%flInitDraBas) then      ! one-time macropore ZDraBas setup
   ... set ZDraBas ; flInitDraBas = .false.
else                              ! normal path: bocodrb/divdra -> qdra
   ... compute drainage ...
end if
```

`flInitDraBas` defaults `.true.` and is only cleared inside the **heavy
surfacewater init**, which is gated on `swsec==2` (i.e. swdra=2 surface-water
cases). For **basic-drainage cases (swdra=1 — the entire hupselbrook family)**
the heavy init never runs, so `flInitDraBas` is still `.true.` at the first
timestep. The first `drainage()` call therefore took the **init-only** branch
and computed **zero `qdra`**. With no drainage sink, the first Richards solve
sees a residual of 0 (the profile is at hydrostatic equilibrium) → converges in
`numbit=1` with the groundwater level unchanged.

The adaptive-dt controller doubles `dt` whenever `numbit ≤ 3`. So the wasted
no-op first step (numbit=1) **doubled `dt` prematurely**, and every subsequent
step ran on a `dt` sequence offset from 4.2.0 → the pervasive ~1e-5/day drift,
amplified to ~1.9 cm GWL under hysteresis, ~0.04 cm under frost, and into a
~1000× performance spiral under swdrought=2 (via the solute stability clamp).

Legacy SWAP 4.2.0 did the `ZDraBas` setup in a **separate init-phase call**
(`Drainage(task=1)`) *before* the time loop; its first per-step call
(`task=2`) computed real drainage. The strangler refactor collapsed both into
the per-step routine but made them mutually exclusive.

### Proof (instrumented trace, swcf3 case, first timestep)

| build | step-1 `qdra` | step-1 `numbit` | step-1 GWL |
|-------|--------------|-----------------|------------|
| before fix | 0.0 | 1 | −75.000000 (no move) |
| after fix  | 0.0715 | 4 | **−75.055108** |

`swap420gf` solves the first step in `numbit=4`, GWL −75.0 → **−75.055** (this
exact value was already documented in the 2026-05-27 trace, point 5). The fixed
modern build reproduces −75.055108 — i.e. it now matches 4.2.0's first step.

### Fix

Make the `flInitDraBas` block a standalone one-time init that **falls through**
to the (now unconditional) drainage computation on the same call, so the first
stepping call does init **and** computes `qdra`, matching legacy's
init-then-step order. At the first call `time%t1900 == tstart` (timecontrol
advances after drainage), identical to legacy's init-time value, so `ZDraBas` is
bit-identical; for basic drainage `ZDraBas` is dead anyway (macropore-only,
retired by ADR 0040). Surgical: only the first `drainage()` call of swdra=1
cases changes behaviour; swdra=2 cases (flInitDraBas already cleared) and
no-drainage cases are untouched.

### Verified

- `swinter3`: xfail → **xpass** (was the adaptive-dt amplifier case).
- `swcf3_maize`: xfail → **xpass** (the "0.01 cm GWL FP artifact" was this desync,
  not floating-point rounding).
- All 6 byte-identical local switch cases still pass (fix preserves byte-identity).
- `soilhysteresis` + `winter` share the identical root cause and are expected to
  xpass, but could not be re-run here (private `tests/swap-cases` submodule
  unavailable). Re-run `check-full` with swap-cases present and remove their
  `known_divergence` flags once confirmed.

### 2026-06-11 follow-up — swdrought2 perf wall is very likely now unblocked

The 2026-06-01 swdrought2 investigation concluded its ~1000x perf collapse was
the SAME adaptive-dt / Richards-convergence sensitivity behind the other
known-divergences, amplified through the solute stability clamp, and that the
fix was "the core soil-water parity work (make modern q/theta match 4.2.0)."

The drainage first-step fix above IS that parity fix — it makes the modern dt
sequence (and therefore q/theta) match 4.2.0 from step 1. So the solute
sub-stepping should no longer collapse `dtsolu` to dtmin under swdrought=2.
**Recommend re-attempting the swdrought2 restoration in a dedicated session**
now that the root cause is fixed. NB: dev-docs/wip/swdrought2-jvl-restoration.patch
no longer applies cleanly (rootextraction.f90 moved in the ADR-0048 feature-first
reorg, `src/crop/rootextraction.f90`); recover jongvanlier.f90 from 5c82f0a^ and
re-port against the current tree rather than `git apply`. The wiltpoint fix
(use state%crop%common%wiltpoint, not hlim4, at the three swdrought=2 matric-flux
sites) is the key correctness fix from that patch to carry forward.

---

## 2026-06-11 — regression suite made self-contained; hysteresis/snow confirmed; frost isolated

**swap-cases toml/ tree is lost.** The `tests/swap-cases` submodule pinned a
commit the upstream remote no longer contains; current `main` has only legacy
ASCII cases (1-6), no `toml/`. The modern TOML base cases were generated locally
and never pushed. Reconstructed base **hupselbrook** TOML from the local switch-
case clones (each patches only specific `.crp` files): swap.toml/dra/csv/template
+ maizes (from swcf3, unpatched) + potatod/grassd (from swrd2, unpatched). Modern
on the reconstruction matches the committed `hupselbrook_reference_gf.json`
byte-for-byte → reconstruction validated.

**Suite is now self-contained.** Promoted hupselbrook + the reconstructed cases
to LOCAL cases under `tests/regression/cases/` so the whole suite no longer
depends on the (private, history-rewritten) submodule. `gen_switch_cases.py` now
sources the local hupselbrook. Retired the 4 orphaned base cases
(grassgrowth/oxygenstress/salinitystress/surfacewater) whose modern TOML is lost
and which need from-scratch conversion (legacy .met meteo, swbotb=3, swdra=2,
swinco=3); legacy ASCII + fixtures preserved.

**soilhysteresis (SWHYST=1) — now byte-identical.** Confirmed the drainage
first-step fix resolves the hysteresis adaptive-dt divergence. The residual that
remained in the first reconstruction was a missing `tau = 0.2` in the modern TOML
(legacy `.swp` has `TAU = 0.2`; modern silently defaults `tau` to 0.0 when
swhyst≠0 — a degenerate value that flips the wetting/drying scanning curve every
step). With `tau = 0.2` the case passes. The hysteresis kernel is faithful.
NB robustness gap: modern should validate/require `tau>0` when `swhyst≠0`.

**snow (SWSNOW=1) — byte-identical.** Snow accumulation/melt path reproduces
4.2.0 exactly (new local `snow` case, SNOWINCO=0).

**winter / frost (SWFROST=1) — genuine divergence, NOT the drainage bug.**
Isolation: snow-only passes byte-identical; frost-only fails (DRAINAGE/RUNOFF/GWL
diverge ~0.5-0.8 cm in the cold grass year, and the magnitude GROWS with frost
intensity — SNOWINCO=0 frost-only shows ~0.8 cm vs the ~0.04 cm the old snow+
frost case showed, because snow insulates the soil and limits frost). The frost
reduction code (FrozenCond/FrozenBounds, src/heat/frozencond.f90) reads faithful.
Since base hupselbrook (swhea=1/swcalt=2 numerical heat) is byte-identical, the
most likely cause is a heat<->frost feedback (frost reduces flow → different
theta → different thermal conductivity → different tsoil near the freezing
threshold → different rfcp) amplifying a sub-threshold difference — same CLASS as
the salinity/oxygen feedback divergences, not the (fixed) drainage first-step bug.
Registered `winter` as a known_divergence for a dedicated session. Needs legacy-
vs-modern tsoil/rfcp per-node instrumentation to confirm.

### 2026-06-11 — NEW BUG surfaced: swbotb=7 (free drainage) blows up in some years

While adding bottom-boundary coverage, a hupselbrook + SWBOTB=7 (free drainage)
case revealed a modern-only bug: 2002/2003 match swap420gf to ~0.18 cm, but 2004
(the grass year, deep/dry GWL) EXPLODES — QBOTTOM 1099 vs 13.9, RUNOFF 999.62 vs
0.0 (999 is a SWAP below-profile sentinel), DRAINAGE 109 vs 26. The legacy oracle
is well-behaved; only modern diverges. Signature: with free drainage the column
drains until GWL drops below the profile and the modern code appears to feed the
999 sentinel into runoff/qbottom arithmetic instead of clamping as 4.2.0 does.
NOT investigated/fixed (distinct from the drainage first-step bug). The case was
dropped (not registered). Reproduce: hupselbrook with `swbotb = 7` in swap.toml
and `SWBOTB = 7` in the .swp. Likely in the swbotb=7 branch of headcalc_residual
(src/soilwater/soilhydraulics.f90:702) or calcgwl's below-profile handling.
A dedicated-session item; bottom-boundary cases (swbotb 1/2/3/5/7/8) remain
uncovered by the regression suite.

---

## 2026-06-11 — 4th lost case (salinitystress) reconstructed; tsoil_file gap fixed

**ENGINE BUG FIXED — swinco=3 + numerical heat was broken.** The validator
required `soil.initial.tsoil_file` for swinco=3 + swhea=1 + swcalt=2, but NO code
ever loaded it into `heat%tsoil` (temperature_seed's tsoil-init branch is guarded
`swinco .ne. 3`). So swinco=3 numerical-heat runs left tsoil unseeded. Fix: let
the existing `[heat].tsoil_init` table seed tsoil for swinco=3 too (afgen at each
node; a full per-node table reproduces the legacy swap.ini warm-restart exactly,
since the table depths equal the node centres), and relax the validator to accept
`tsoil_file` OR `tsoil_init`. Also added `*.irg` to regen_reference's
LEGACY_INPUTS (was missing → swap420gf couldn't run swirgfil=1 cases).

**salinitystress reconstructed, runs, ~3% residual (NOT byte-identical).**
Full swinco=3 warm-restart: h_init.csv + cml_init.csv seeded per-node (numnod=195),
atmosphere ldwet/atmin7 in [soil.initial], tsoil via the 195-row tsoil_init table.
swbotb=3 sinus aquifer head, dramet=3 2-level drainage, wofost potato (= hupselbrook
potatod + dvsend=3/swgerm=0/salinity-on/relmf=0.8/one extra frtb row), 585 fixed
irrigation events (irrig.csv), Maas-Hoffman salinity.

Diagnosis of the residual (CWSO/CPWSO ~3%, CONC ~0.1, from year 1):
- **NOT salinity.** saltslope=0 gives the IDENTICAL divergence (CPWSO is potential,
  salinity-independent). So this is not the salinity→cml feedback class.
- **NOT warm-restart incompleteness.** numnod=195; h/cml seeded per-node directly,
  tsoil afgen-exact.
- **NOT the general solute code.** Adding CONC/CWSO assertions to hupselbrook shows
  CWSO/CPWSO byte-identical (0.0); CONC diverges only transiently (final matches).
- **Specific to this case's feature mix** (swinco=3 warm-restart + swbotb=3 Cauchy
  + 585 irrigation-solute events). Needs targeted instrumentation (compare modern
  vs swap420gf cml per node over the first days). Registered known_divergence.

**Net: all 4 lost base cases reconstructed & in the regression — grassgrowth,
oxygenstress, surfacewater byte-identical; salinitystress runs with a documented
~3% solute residual.**

---

## 2026-06-11 — salinitystress deep dive (legacy-vs-modern instrumented diff)

Built a logging copy of the legacy SWAP 4.2.0 (worktree from the `v4.2.0` tag +
TTUTIL source from github SWAP-model/ttutil, compiled with the swap420gf flags).
Verified byte-identical to the swap420gf oracle on salinitystress (only the
timestamp comment differs), then added matching per-node / per-step logging to
both the legacy and modern solute paths to locate the ~3% divergence.

Findings (in order):
1. **Initial profiles byte-identical.** Dumped z/h/theta/cml per node right after
   solute init: max|diff| = 0.0 across all 195 nodes. The swinco=3 warm restart
   (h_init/cml_init per-node, tsoil via the table) is exact.
2. **The water diverges too, slowly.** GWL drifts from ~0 to ~3 cm over 4 years;
   DRAINAGE/QBOTTOM ~0.07/0.03 cm. RAIN/IRRIG identical. So the solute ~3% is
   DOWNSTREAM of a slow water (GWL) drift, not a solute-transport bug per se
   (confirmed: hupselbrook's CWSO/CPWSO are byte-identical when asserted).
3. **The drift is an adaptive-dt sequence desync.** Per-step dump of (t, dt,
   numbit, gwl): steps 0-2 byte-identical; at step ~3 the modern dt DOUBLES
   (1e-6 -> 2e-6) while legacy holds 1e-6 one step longer — even though **numbit
   is identical (2) on both**. So it is the dt-controller's event timing
   (tEvent/dtEvent/flprevious evolution), not a convergence/numbit difference.
   Day-1 step counts: legacy 62, modern 58.
4. **Ruled out:** salinity (saltslope=0 gives the identical divergence), rain
   intensity (swrain=0 still diverges), init, and a numbit difference. The
   timecontrol dt/event code itself diffs as a faithful transcription.
5. **swinco=3-specific.** The other reconstructed cases (grassgrowth/oxygenstress/
   surfacewater) share swrain=2 but are swinco=2 and are byte-identical.

**Real bug fixed along the way:** the swinco=3 initial timestep. Legacy clamps
the restart dt to dtmin (its per-step `dt = max(dt, dtmin)`), so a sub-dtmin
swap.ini dt (1e-7 < dtmin 1e-6) becomes dtmin on step 1. The modern set
`time%dt = config%soil%initial%dt` unclamped, running the first step 10x too
small. Fixed in swap_mod.f90 (clamp to [dtmin, dtmax]) — steps 0-2 now match
4.2.0 exactly. It does not by itself close the aggregate gap (the step-~3 event-
timing desync dominates), but it removes a genuine discrepancy.

**Status:** the residual is a sub-threshold dt-controller event-timing sensitivity
in the swinco=3 path — the same hard adaptive-dt class as the winter/frost
divergence. Not a pinpointable transcription bug (init + numbit + controller code
all match). Closing it needs instrumenting tEvent/dtEvent/flprevious *inside* both
timecontrol modules step-by-step to find the first variable that diverges.

### 2026-06-11 (cont.) — salinitystress: 2nd bug fixed (dtprevious), residual is FP

Continued the legacy-vs-modern dt-controller trace by logging tEvent/dtEvent/
flprevious/dtprevious/numbit/gwl inside BOTH timecontrol modules:

**2nd control-flow bug FIXED — swinco=3 dtprevious.** Entry-state trace showed
step-0 `dtprevious`: legacy 1e-6 vs modern 2e-4. Legacy clamps the sub-dtmin
restart dt to dtmin and carries that as dtprevious; modern's timecontrol_init
took the `dt<dtmin -> sqrt(dtmin*dtmax)` branch leaving dtprevious=2e-4. On step 0
the event-limit branch sets `dt = dtprevious`, so modern's 2e-4 spuriously
triggered an event-limit -> different flprevious/flTnext -> the dt sequence
desynced from step 1. Fix: in swap_mod.f90 set BOTH time%dt and time%dtprevious
to the clamped restart dt. After this, the dt sequence AND numbit are
byte-for-byte identical to 4.2.0 for ~50 days (GWL firstdiff moved from row 1 to
row ~51).

**Residual = floating-point ordering.** At ~day 50 the GWL diverges by <1e-7 with
**identical dt, numbit, AND flprevious** on both builds. Same control flow + same
inputs + different output = a pure FP operation-ordering difference in the solve
(the modernization changed FP order via state-record indirection / associate /
loop restructuring). On this sensitive config (195 nodes, 1cm top compartments,
4-year run, swbotb=3 Cauchy, numerical heat) it accumulates to ~2.9 cm GWL and
~3% solute by year 3-4. The simpler byte-identical cases don't accumulate past
the 2-decimal tolerance.

**Conclusion:** the two control-flow bugs are fixed (genuine faithfulness gains);
the remaining ~3% is the irreducible FP-ordering limit, not a transcription bug.
Closing it would require matching the exact FP operation order of 4.2.0 in the
refactored solve — out of scope and high-risk. salinitystress stays a documented
known_divergence. Both byte-identical 4.2.0 control-flow bugs were found only by
the instrumented legacy build (worktree from v4.2.0 + TTUTIL source).

### 2026-06-11 (cont.) — salinitystress RESOLVED: byte-identical, NOT FP-irreducible

The "irreducible FP-ordering" conclusion above was WRONG. Re-instrumented the
legacy-vs-modern diff (rebuilt the v4.2.0+TTUTIL gfortran oracle, verified
byte-identical to swap420gf, added per-column daily diffs + tz-profile + reduceva
traces). The ~3% residual was THREE input-transcription bugs in the reconstructed
TOML/state mapping, all now fixed. The case is byte-identical to swap420gf across
all 8 asserted columns × 1461 days (and the full 195-node H/WC/CONC tz-profile).

1. **solute.ldis scalar broadcast bug (src/state/solute_state.f90).** The scalar
   `ldis = 5.0` shorthand only seeded `self%ldis(1)`; every deeper soil layer ran
   with ldis=0 (zero mechanical dispersion below layer 1). The legacy .swp LDIS is
   a per-layer column (`LDIS = 5.0 5.0`), and read_solute_toml.f90's own comment
   already says the scalar form "broadcasts to all layers" — the adapter just
   wasn't doing it. Fixed: `self%ldis(:) = config_solute%ldis`. This was masked in
   every other case because only salinitystress asserts CONC directly. Surfaced as
   a day-1 solute-profile divergence growing to −75% at mid-profile.

2. **swredu wrong model (cases/salinitystress/toml/swap.toml).** Authored
   `swredu = 1` (Black) + cofredbl/rsigni; the legacy .swp is `SWREDU = 2`
   (Boesten–Stroosnijder, `COFREDBO = 0.54`). Wrong soil-evaporation reduction
   model → EACT diverged from day 49 (first dry-down after the wet start).
   Confirmed via matched reduceva(task=1) traces: ldwet/dt path identical, only
   the model branch differed. Fixed to swredu=2 + cofredbo=0.54.

3. **FRTB dropped row (cases/salinitystress/toml/potatod.crp.toml).** The legacy
   potato FRTB has 5 rows incl. `1.27 0.0`; the TOML kept only 4
   (`0.00 0.2 / 1.00 0.2 / 1.36 0.0 / 2.00 0.0`), so the root fraction declined
   over [1.00,1.36] instead of [1.00,1.27], leaving FR≈0.05 at DVS 1.27 (legacy 0).
   More dry matter to roots, less to storage organs → CPWSO/CWSO (POTENTIAL, hence
   salinity/water-independent) ran ~2–4% low from tuber initiation. Restored the
   `1.27 0.0` row → CPWSO/CWSO byte-identical.

The two earlier swinco=3 control-flow fixes (initial dt + dtprevious clamp) remain
correct and load-bearing for the dt/numbit sync; they were necessary but not
sufficient. Lesson: a "few-percent residual that grows from year 1" on a
reconstructed case is far more likely a transcription bug than FP ordering —
diff EVERY asserted column against the oracle (not just the aggregate) and split
soil/crop/solute before reaching for "irreducible." known_divergence removed;
case registered as a normal pass.

## 2026-06-11 (cont.) — swsalinity1 + swoxygen2 RESTORED (byte-identical)

Reviewing the three "extremely fast" (0.02-0.03s) xfail cases: all three were
honest input-validation gates (the binary fatal-errors at swap_init with an
explicit "not yet supported in the TOML pipeline" message), NOT silent passes.
Two of the three turned out to be cheap, real recoveries:

**swsalinity1 (cropfixed maize, swsalinity=1 Maas-Hoffman) — RESTORED.** The
2026-05-31 "salinity→cml feedback divergence, irreducible" diagnosis was WRONG.
It rested on "salinitystress passes" as evidence the feedback was fine — but
salinitystress only "passed" because it was flagged known_divergence. Once
salinitystress was fixed to genuinely byte-identical (incl. TREDSOL, the salinity
reduction), it PROVED the salinity→cml path is exact. The cropfixed divergence was
the SAME solute.ldis broadcast bug: salinity is the only stress reading sol%cml,
so zero dispersion below layer 1 → wrong cml → wrong Maas-Hoffman reduction →
TACT drift (~3.4 cm/yr). The ldis fix corrects cml → byte-identical. Lifted the
two gates (cropfixed_config_validate + cropfixed_init runtime guard); state
plumbing (saltmax/saltslope/salthead) was already wired. swsalinity=2 stays gated.

**swoxygen2 (wofost potato, swoxygen=2 Bartholomeus) — RESTORED.** Never wired,
only stub-errored. The Bartholomeus kernel (oxygenstress.f90) was already live for
the grass path (oxygenstress case). cropwofost_init now copies the Bartholomeus
fields to state%crop%oxygen / %common (mirroring cropgrass_init) and sets
swoxygentype=1 — legacy readwofost's default; WOFOST crps have no SWOXYGENTYPE
field, and swoxygentype=2 (reproduction functions) is the grass-only alternative,
still dormant in rootextraction.f90. Lifted the config + init gates. Byte-identical
to swap420gf (TACT carries the oxygen-stress signal). Updated the two now-stale
"rejected" config unit tests to "accepted" (assertTrue→assertFalse on the gate
message).

**Grass swsalinity left gated** (no grass-salinity regression case to prove it),
but its comment is corrected: it is very likely fixed by the same ldis change,
since it uses the identical rootextraction.f90 kernel. Lift once a case exists.

**swdrought2 (De Jong van Lier) — still the ONE genuine port.** Compute deleted
(5c82f0a). The preserved restoration patch (dev-docs/wip/swdrought2-jvl-
restoration.patch) no longer applies — it predates the crop feature-folder reorg.
Its documented perf-block (Richards/adaptive-dt convergence spiral, 2026-06-01)
predates the drainage first-step no-op fix (2026-06-11) that cleared the SAME
adaptive-dt class for hysteresis/winter/swinter3, so the wall is very likely now
gone. Prime candidate for a dedicated session: re-port the 579-line kernel against
the current layout, then re-check perf. Left as pending_restore.

Suite: 18 passed / 0 failed / 2 xfail (winter frost-feedback, swdrought2). The
winter case genuinely runs (~2.3s) and diverges in the frost path — a separate
real heat↔frost numerical item, not a gate.

## 2026-06-11 (cont.) — swdrought=2 (De Jong van Lier) restoration: spiral bug found, perf residual

Attempted the full JvL restoration (the deleted compute, not just a gate). The
preserved WIP patch (dev-docs/wip/swdrought2-jvl-restoration.patch) no longer
applies, but re-applying it against the current layout was mostly mechanical:

- Kernel (JongvanLier + JongvanLierLoop) merged into rootextraction_mod (siblings
  of matric_flux, to avoid a module cycle), `use variables` reads aliased to
  state%crop%common via associate. All soilwater JvL fields (hleaf/Hxylem/mflux/
  mroot/hroot/rootrho/rootphi/rmax/qrosum/alpJvLier) already exist; only the
  crop_common scalar params (kroot/kstem/rxylem/rootradius/rootcoefa/rooteff/
  stephr/criterhr/wiltpoint/taccur) needed adding. Config/init/reader/dispatch
  wired; gen_switch_cases.py emits the JvL .crp/.toml params. It COMPILES and runs
  the physics.

**KEY FINDING — a real bug behind the "perf wall".** The 2026-06-01 conclusion
("Richards/adaptive-dt convergence spiral, needs core-headcalc changes, out of
scope") was — like the salinitystress/swsalinity1 "irreducible" calls — partly a
missed bug. `cropfixed_runtime` (and grass/wofost runtime) compute the
wilting-point water content as `twilt(i) = watcon(crop%common%hlim4, ...)`, but
legacy is `twilt(i) = watcon(i, wiltpoint)` — the JvL LEAF wilting head
(-20000), NOT the Feddes hlim4 (-8000). `twilt` is read ONLY by the JvL qmax
limit, so the wrong (too-wet) twilt starved qmax, which fed a bad qrot into the
Richards solve and spiralled the dt-reduction retry loop INFINITELY (the run hung,
0 outer steps). Fixing twilt to use `wiltpoint` removed the spiral entirely:
cumretry went 0, the run progresses normally pre-crop (124 days in ~30s).

**RESIDUAL — still perf-blocked, but differently.** After the twilt fix the run
no longer hangs, but it is still too slow to finish (a 9-month window doesn't
complete in 7 min): during the crop transpiration season dt stays tiny
(~2.5e-3) with cumretry=0 — i.e. the adaptive-dt CONTROLLER itself picks small dt
(high numbit), not the forced-reduction loop. Legacy runs the full 3yr in 0.69s,
so its headcalc converges fast with the same qrot. Two candidates remain: (a) the
modern JvL qrot still differs subtly from legacy (faithfulness — fixable), or (b)
genuine core-headcalc convergence sensitivity (the prior session's hypothesis).
Distinguishing needs an instrumented-legacy per-step qrot/dt/numbit comparison
(the method that cracked salinitystress) — a dedicated session.

**DISPOSITION.** Reverted the restoration to keep the tree clean (the project
convention is that dormant compute is uncompiled) and the suite green; swdrought2
stays pending_restore (gated, fast-fails). The twilt bug + the "spiral was a bug,
not core-headcalc" finding are the resumption lead. NEXT SESSION: re-apply the
restoration (kernel into rootextraction_mod + the config/init/state/reader wiring
+ gen params, all worked out above), apply the twilt fix (watcon(wiltpoint,...) in
cropfixed/grass/wofost runtime), then instrument modern-vs-legacy qrot per step to
close the dt residual.
