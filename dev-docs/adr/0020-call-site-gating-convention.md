---
title: "ADR 0020 — Call-site gating convention for optional subsystems"
date: 2026-05-05
status: accepted
---

# ADR 0020: Call-site gating convention for optional subsystems

## Context

The SWAP runtime composes ~10 optional subsystems on top of the core
water-balance solver: crop nutrients, macropore physics, surface water,
hysteresis, salinity, snow, frost, oxygen stress, transpiration
reduction, etc. Each is opt-in: a TOML switch (`swcropnut`, `swmacro`,
`swsurfwat`, …) decides whether the subsystem participates in a given
run.

The dominant pattern in `src/core/swap.f90` is **call-site gating**:

```fortran
if (flCropNut)   call SoilManagement(*)
if (flMacroPore) call MacroPore(*)
if (flSurfaceWater) call surfacewater(*)
```

Each gate flag (`flCropNut`, `flMacroPore`, …) lives in
`src/core/variables.f90`, defaults to `.false.`, and is set by the
TOML→variables adapter (`config_to_variables.f90`) from the
corresponding `config%<group>%sw<x>` integer. The subsystem body
itself assumes it was only invoked when the gate is true; it does not
re-check.

Two outliers persisted through Phase 4f:

- `Read_Tillage` (`src/crop/tillage.f90`) carried an internal
  `if (swtill /= 1) return` self-check.
- `SSDI_irrigation(1)` (`src/crop/irrigation.f90`) carried an internal
  `if (swssdi == 0) return` self-check.

Both were unconditionally invoked from `swap.f90` /
`timecontrol.f90`, with the subsystem body deciding whether to do
real work. This split convention made the code harder to audit:

1. SS-9/SS-10 audits had to trace the self-check to know whether a
   subsystem was actually wired.
2. SS-10.5 added validator stub-errors rejecting `swtill=1` and
   `swssdi=1`, treating both features as "deferred." That was wrong —
   the legacy TTutil-based parameter readers for both subsystems are
   still operational; they simply require the right call-site wiring.
3. SS-C (legacy reader physical deletion) needs the call-site gating
   in place before it can safely strip `case(0)` short-circuits from
   the subsystem bodies.

## Decision

**Codify call-site gating as the single convention** for optional
subsystems in the SWAP runtime. Concretely:

1. Each optional subsystem has a `logical :: flX = .false.` flag in
   `src/core/variables.f90`.
2. The TOML→variables adapter sets `flX = (config%group%swX == 1)`.
3. Every call site in the runtime is wrapped with
   `if (flX) call X(...)`.
4. The subsystem body assumes the flag is true on entry. It does not
   self-check; it does not maintain disable-defaults inline.
5. Validators check `swX` is a member of its enum (e.g. `[0, 1]`).
   They do not stub-error legitimate values.

Apply this convention to `DoTillage` and `SSDI_irrigation` to bring
them in line:

- Add `flTillage` and `flSSDI` to `variables.f90`.
- Adapter populates them from `config%soil%swtill` and
  `config%irrigation%swssdi`.
- Wrap the 3 `call DoTillage(*)` sites in `swap.f90` and the 3
  `call SSDI_irrigation(*)` sites (2 in `swap.f90`, 1 in
  `timecontrol.f90`) with the flag tests.
- Drop the self-checks from `Read_Tillage` and `SSDI_irrigation(1)`.
- Remove the SS-10.5 stub-error blocks for `swtill=1` and
  `swssdi=1` from the validators.

## Consequences

**Audit:** A grep for `if (flX)` in `swap.f90` and `timecontrol.f90`
gives a complete inventory of the optional subsystems wired into the
runtime. Subsystem bodies are no longer split between "real work" and
"disable defaults"; bodies are pure.

**Validators:** `swtill=1` and `swssdi=1` become accepted enum
values. The TOML pipeline now activates the legacy TTutil-based
parameter readers for both subsystems via the call-site gate. No
behaviour change for any existing regression case (all use
`swtill=0`, `swssdi=0`).

**Future work (ADR 0021 candidate):** This ADR does **not** decide
whether `Read_Tillage` and `SSDI_irrigation(1)` keep their TTutil-
based reading of staged `swap.swp`. They currently do, because no
schema slots exist for the tillage block (Date_tillage, Z_tillage,
I_tillage, Type_tillage, iType_Tillage, TAB_Rho_*, TAB_K_R_cons,
TAB_N_match, …) or the SSDI block (ssdi_file companion table,
ssdi_date, ssdi_amount, …). Porting those parameter blocks to TOML
schema is a separate future arc tracked as ADR 0021.

**Test infrastructure:** The two SS-10.5 stub-error tests
(`test_soil_swtill_one_stub_errors`,
`test_irrigation_swssdi_one_stub_errors`) are inverted to assert the
acceptance of `=1`. The `flTillage`/`flSSDI` plumbing itself is
exercised by the unit-test gate (validator passes) and the regression
gate (all 5 cases continue to pass with flags `.false.`).

**Code-shape:** Subsystem headers in `tillage.f90` and
`irrigation.f90` reference this ADR to document the convention shift.

## Related

- ADR 0007 — config parse-validate-finalize-adapter pipeline (where
  `config_to_variables.f90` lives).
- ADR 0015 — Strangler narrow-scope stub-errors (which the SS-10.5
  swtill/swssdi blocks were a misuse of; this ADR retires them).
- ADR 0019 — Legacy readers retired from runtime; retained as parity
  fixtures (the broader Phase 4f-extend context).
- ADR 0021 (candidate) — TOML port of the tillage and SSDI
  parameter blocks.
