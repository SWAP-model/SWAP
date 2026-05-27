# Solute pure-kernels pilot — design

- **Date:** 2026-05-27
- **Module under work:** `src/solute/solute.f90` (`solute_mod`)
- **Status:** approved design, pending implementation plan
- **Type:** readability/testability refactor + latent-bug fix (pilot)

## Motivation

`solute_step` is a single ~230-line subroutine that interleaves orchestration
(timestep loop, cohort resets, boundary handling) with physics formulas
(Freundlich isotherm inversion, decomposition factors, dispersion). The physics
is buried in the loop body and has **zero unit tests** — `tests/unit/solute/`
exists but contains only a `.gitkeep`. All current solute tests cover
config/state/TOML parsing, never the transport math.

The codebase already has a proven shape for testable physics:
`src/atmosphere/interception.f90` and `src/atmosphere/et.f90` expose `pure`
functions with explicit scalar/array args and no `state` dependency
("Pure function — explicit args, no state dependency"), with table lookups
hoisted to the caller. This pilot extends that pattern into `solute` and
evaluates whether it should be transposed to other subsystems.

## Finding that reshaped the work (latent bug)

While scoping the refactor we confirmed a latent correctness bug:

`solute_step` reads five derived coefficient arrays — `ddiffwcs`, `bdenskf`,
`bdenskfcref`, `bdenskfsatporos`, `decpotfdepth` — that it **never assigns**.
They are declared as locals in both `solute_seed` and `solute_step`, computed
(correctly) only in `solute_seed`, then discarded when `solute_seed` returns.
`solute_step` therefore reads uninitialised local arrays.

Origin: the task-dispatch split (commit `52ac7ff`, "solute — drop task
dispatch, `solute_seed`/`step(state)`"). In the pre-split `solute(task, state)`
these arrays were computed in the `task=1` branch and read in `task=2/3`; the
split turned them into per-call locals in two *separate* procedures, severing
the data flow.

**Why it is currently masked.** The only solute-active regression case is
`salinitystress` (conservative salt). There:

- `ddif = 0` → correct `ddiffwcs = 0`
- `decpot = 0` → correct `decpotfdepth = 0`
- `kf = 0` → correct `bdenskf = bdenskfcref = 0`
- `bdenskfsatporos` is read only inside the `swbr == 1` aquifer-breakthrough
  branch, which this case does not enter.

gfortran zero-inits (or stack-reuses) the automatics, and zero equals the
physically-correct value here, so `salinitystress` passes byte-identical
against the gfortran-4.2.0 reference. Verified: `regression salinitystress`
→ "regression ok (annual stats match fixture)".

The bug only bites a case with **sorption (`kf>0`), decay (`decpot>0`),
molecular diffusion (`ddif>0`), or aquifer breakthrough (`swbr=1`)** — none
present in the current suite. These five coefficients are all time-invariant
soil/solute properties, so the correct design is to compute them once and
persist them in typed state, which both fixes the bug and makes the seed→step
data flow explicit.

## Goal

Make `solute_step` read as a sequence of named operations, and make its physics
independently unit-testable, by extracting:

1. derived coefficient computation (anchor; also fixes the latent bug),
2. the Freundlich isotherm inversion,
3. the decomposition factors,

into `pure` functions in a dedicated kernels module — using `solute` as a pilot
to decide whether the pattern is worth transposing to `soilwater`/`crop`/`heat`.

## Architecture

### New module: `src/solute/solute_kernels.f90`

`module solute_kernels_mod` — `pure` functions only, explicit args, no
`swap_state_t` dependency. Mirrors `interception_mod`/`et.f90`. `solute_mod`
`use`s it; pFUnit suites `use solute_kernels_mod` directly.

### Extraction 1 — derived coefficients (anchor + bug fix)

State schema change in `src/state/solute_state.f90`:

- Add five `real(real64), allocatable :: name(:)` fields to `solute_state_t`:
  `bdenskf`, `bdenskfcref`, `bdenskfsatporos`, `ddiffwcs`, `decpotfdepth`,
  each documented in the existing per-node comment style.
- Allocate them (sized `numnod`) in `solute_state_init`, alongside the existing
  `cml`/`cmsy` allocation (init allocates; seed fills — the established
  two-phase pattern).

Kernels (Shape B — per-coefficient `elemental` functions, chosen over a single
builder for test granularity, call-site readability, and a strictly state/mesh-
free contract; closest to the `interception.f90` precedent). In
`solute_kernels_mod`:

- `elemental function bdenskf_coeff(bdens, kf) result(v)` — `bdens*kf`
- `elemental function bdenskfsatporos_coeff(bdens, kfsat, poros) result(v)`
  — `bdens*kfsat + poros` (scalars `kfsat`/`poros` broadcast over `bdens`)
- `elemental function ddiffwcs_coeff(ddif, thetsl) result(v)` — `ddif/thetsl**2`
- `elemental function decpotfdepth_coeff(decpot, fdepth) result(v)`
  — `decpot*fdepth`

The trivial derived product is **not** wrapped: `solute_seed` computes
`sol%bdenskfcref = sol%bdenskf * cref` as a plain line after `bdenskf`.

Layer→node gathering of layer-indexed soil/solute properties happens at the
call site (in `solute_seed`) — the elemental kernels only ever see per-node (or
scalar) values, so they stay genuinely state- and mesh-free. `solute_seed`
assigns each result into `state%solute` (e.g.
`sol%bdenskf = bdenskf_coeff(bdens_node, kf_node)`), reading as a narrative
sequence of named formulas.

`solute_step` reads `sol%bdenskf(i)` etc. The dead local declarations are
**deleted** from both routines.

Byte-identical on the current suite (correct values equal the masked zeros in
`salinitystress`); the regression gate remains a valid safety net.

### Extraction 2 — Freundlich isotherm inversion

`pure function solute_cml_from_cmsy(cmsy, theta, bdenskf, frexp, cref, cml_guess) result(cml)`
covering the three branches at lines 227–244:

- `cmsy < vsmall` → `cml = 0`
- `|frexp - 1| < 0.001` (linear) → `cml = cmsy / (theta + bdenskf)`
- otherwise → fixed-point iteration to relative tolerance `rer`, seeded from
  `cml_guess` (clamped to `vsmall`).

`vsmall`/`rer` remain module parameters; pass or re-declare in the kernel as
appropriate. `solute_step` calls this per node, replacing the inline loop.

### Extraction 3 — decomposition factors

Small `pure` functions for the temperature factor (inputs `tsoil`, `gampar`,
`fl_temperature` logical), the moisture factor (inputs `theta`, `rtheta`,
`bexp`), and the transformation rate `ctrans` (composing `decact`, `theta`,
`cml`, `bdenskfcref`, `cref`, `frexp`). `solute_step` composes them.

## Data flow (after)

```
solute_state_init   : allocate cml, cmsy, + 5 coefficient arrays (numnod)
solute_seed         : fill cml profile; gather layer→node; call coeff builder
                      → store bdenskf/bdenskfcref/bdenskfsatporos/ddiffwcs/
                        decpotfdepth into state%solute
solute_step (loop)  : read sol%<coeff>(i); per node call
                      solute_cml_from_cmsy(...) and decomposition kernels
```

## Testing

Populate `tests/unit/solute/` with the first solute physics unit tests:

- **Coefficients:** per-function scalar golden values (`bdenskf_coeff(2,3)==6`
  etc.); explicit zero-input cases (`kf=0`, `decpot=0`, `ddif=0`) that reproduce
  today's masked-to-zero behavior; a non-trivial `bdenskfsatporos_coeff`
  (`poros>0`) case that the current code never exercises; `ddiffwcs_coeff`
  exercised with `thetsl>0`.
- **Freundlich:** linear branch (`frexp≈1`), sub-`vsmall` zeroing, a
  convergence case with a hand-computed fixed point, and `frexp≠1` with known
  inputs.
- **Decomposition:** `fl_temperature=.false.` → `ftemp=0`; the `tsoil<35`
  exponential and the `tsoil≥35` cap; `ftheta` min-clamp at 1.

Register every new suite in `tests/unit/testSuites.inc` and verify the
pFUnit `OK (N tests)` count rises (an unregistered `.pf` compiles but never
runs).

## Verification

- `check-fast` after **each** extraction — byte-identical gate (per the
  per-task regression-gate convention; do not defer to an end-of-arc gate).
- `check-full` before each commit.
- `rm -rf builddir` clean rebuild after the `solute_state.f90` schema change
  (incremental Meson builds do not propagate `.mod` deps across the
  swap_modern↔swap_legacy boundary).
- One commit per extraction.

## Pilot evaluation criteria

After completion, decide on transposition to `soilwater`/`crop`/`heat` from:

1. **Readability:** does the `solute_step` per-compartment loop now read as a
   sequence of named operations rather than inline formulas?
2. **Defect leverage:** the pattern already surfaced one latent bug while
   scoping; did the kernels prevent/expose anything else?
3. **Cost/benefit:** diff size and ceremony (new module, state fields, tests)
   versus the clarity/testability gain.
4. **Generality:** does the `pure` + explicit-args shape survive the messier
   state coupling in `soilwater`/`crop`, or does `solute`'s relatively clean
   per-node structure flatter the pattern?

## Out of scope (YAGNI)

- Splitting `solute_step` into stage subroutines (surface flux / loop /
  boundary / balance). Possible follow-up; not this pilot.
- `solute_seed`'s profile-initialisation logic beyond invoking the builder.
- Activating new regression cases for sorption/decay/diffusion/breakthrough
  (separate coverage arc — but worth doing soon, since it would have caught
  this bug).
- Performance work; renaming the coefficient fields (keep existing names to
  minimise churn).
