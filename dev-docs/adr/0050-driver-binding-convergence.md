# ADR 0050 — Driver / binding convergence (one library, registry-backed)

**Status:** Accepted — sub-arc 1 landed 2026-07-03; sub-arc 2 planned
**Arc:** T1-G′ (Tier-1, from the 2026-07-03 modernization review)
**Relates to:** the review's §7 and the driver-convergence discussion. Does not
supersede an existing ADR (the two-library split was a meson artifact documented
only in build comments, never an ADR).

## Context

SWAP exposes three ways to drive the engine, diverging only at the *backing
store* (all bottom out in `swap_init`/`swap_run_step`/`swap_close` on an explicit
`(state, config)`):

- CLI `swap` (static) — local `state,config`.
- `libswap_bmi.so` — CSDMS BMI + pragmatic CAPI over the singleton `capi_state`.
- `libswap_xmi.so` — XMI over the ensemble `columns(:)`.

Two shared libraries exist **only** because `swap_bmi_mod` and `swap_xmi_mod` both
export 14 identical C symbols (`initialize`/`update`/`finalize`/…) backed by
different stores — a hard link-time collision. Costs: `c_to_f_string`/
`f_to_c_string` duplicated ×3; per-variable `select case` switches hand-maintained
in each facade (≈6 edit sites per new variable, no cross-check); three variable
vocabularies.

MODFLOW 6 is the reference answer: **one** `libmf6` exposing BMI + XMI over one
memory-managed store, a thin `mf6` CLI built from the same core, and variables
resolved generically through a registry (the memory manager) rather than
per-variable code. Python drives it through one binding (`xmipy`/`modflowapi`).

## Decision

Converge SWAP onto the same shape, in two byte-behavior-preserving sub-arcs:

1. **Variable registry (sub-arc 1 — landed).** Introduce `var_registry_t`
   (`src/bindings/swap_var_registry_mod.f90`): a namespaced table mapping a
   variable name to `{units, rank, readable/settable, scale, pointer-into-state}`,
   with generic accessors. Facades resolve variables by lookup instead of
   switches. Namespace tags (`NS_BMI`, `NS_CAPI`, …) keep each facade's
   enumerate/count/name-list operations scoped to its own vocabulary, so the
   change is byte-for-byte compatible with the switches it replaces. The BMI
   facade and the CAPI `swap_view_array` are the first two consumers; the registry
   is built after `swap_init` in both singleton init paths.

2. **One `libswap.so` on the ensemble backing store (sub-arc 2 — planned).** Port
   the CAPI rich accessors + in-memory init + a per-timestep step onto the
   ensemble; make the BMI facade a 1-column ensemble; retire the singleton so the
   14 colliding symbols have a single ensemble-backed definition; merge the two
   `shared_library` targets into `libswap.so`; move XMI `get_value_ptr` and the
   C-string helpers onto the shared registry; update every consumer to load the
   one library. **The CLI stays static** (built from the same core — the drift
   guarantee comes from same-source/same-build, not from dynamic linking).

The registry is the *index-over-existing-state* form (pointers into `swap_state_t`
bound post-init), not MF6's stronger *manager-owns-memory* form; that keeps every
allocation site untouched and is sufficient for the BMI/XMI surface. Deepening
toward the MF6 form is a later option if exchanges need to address arbitrary
internal variables by string.

## Consequences

- Adding an exposed variable becomes one `add_*` line in `build_variable_registry`,
  visible to every metadata call and (after sub-arc 2) every facade — no more
  6-edit-sites-with-no-cross-check.
- Sub-arc 1 is LoC-roughly-neutral for a single facade (the win is the
  single-source-of-truth and the substrate); the LoC reduction and the two-library
  collapse land in sub-arc 2.
- Verified byte-behavior-preserving: BMI `hello_swap` + CAPI `run_ensemble` (driven
  against the self-contained regression case, since the meson `bmi`/`cffi-demo`
  suites need the un-checked-out `tests/swap-cases` submodule), `check-fast` 4/4
  byte-identical, XMI smoke OK, pFUnit 833 → 838.
- The `f_to_c_string` copies were deliberately NOT merged in sub-arc 1: BMI's is
  bounded by `max_len`, XMI's is unbounded — merging is a behavior change to be
  reconciled during the sub-arc-2 facade unification, not a mechanical move.
</content>
