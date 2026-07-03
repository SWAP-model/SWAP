# ADR 0050 — Driver / binding convergence (one library, registry-backed)

**Status:** Accepted — sub-arcs 1 and 2 landed 2026-07-03 (see Addendum)
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

## Addendum (2026-07-03) — sub-arc 2 landed, with one design revision

**Revision.** The original sub-arc-2 sketch proposed "BMI `update()` = one day"
in the merged library. That was a behavior change dressed as a default and was
dropped. The landed design is **mode-dependent lifecycle semantics**:
`initialize()` detects the `ensemble.txt` sidecar next to the config — present →
coupled ensemble mode (former XMI behavior exactly: `update()` = one day,
`get_time_step()` = 1.0, item counts 1/2, `get_var_nbytes` = 8·ncol); absent →
single-column standalone mode (former BMI behavior exactly: `update()` = one
Richards substep, `get_time_step()` = dt, registry-driven metadata). Both
existing consumer populations observe bit-identical behavior.

**Storage unification** happened by pointer indirection, not a port:
`capi_state`/`capi_config` are now pointers bound to the ensemble's
`columns(1)`/`configs(1)` (the ensemble arrays gained `target`), so every CAPI
accessor compiled unchanged while `swap_ensemble_mod` became the sole owner of
state storage. Pre-init `swap_set_headless` (a real consumer call order) is
buffered and applied at bind time. `swap_c_strings_mod` (2a) holds the shared
string helpers; the generic `f_to_c_string` preserves both bounded/unbounded
contracts by arity.

**Symbol unification.** `swap_bmi_mod` is the sole definer of the 14
formerly-colliding C names; `swap_xmi_mod` retains only the XMI-specific verbs
(`prepare_time_step`/`solve`/`get_value_ptr`/`get_version`/rank/shape/
`get_last_bmi_error`, with last-error shared via the ensemble). Two arity
conflicts were resolved toward the XMI (unbounded) forms after verifying the
callers from source: xmipy's `get_value_ptr` calls `get_var_type(name, buf)`
(2-arg), and no BMI consumer calls `get_var_type`/`get_component_name` at all;
`initialize` keeps its 2-arg form with `n` unread (xmipy passes 1 arg).

**Build/consumers.** One `shared_library('libswap', name_prefix: '')` →
`libswap.so`; `libswap_bmi.so`/`libswap_xmi.so` and the `swap_modern_xmi`
static lib deleted; all consumer references updated (meson test registrations,
`tests/coupling/_defaults.py`, `imod_coupler.toml` + its `run_coupled.py`
rewrite string, `prototype/*.py`, docstrings). The CLI remains static-linked
from the same core.

**Verified** on the merged library: BMI hello_swap + CAPI run_ensemble
(bit-identical values, substep semantics intact), XMI smoke + SWAP↔MODFLOW6
coupled smoke (bit-identical exchange values / water table), `check-fast` 4/4
byte-identical, pFUnit 838.

**Deferred to later arcs:** XMI `get_value_ptr` onto the registry (an `NS_XMI`
namespace over the exchange arrays); handle-based multi-instance (per-instance
ensemble + error state); threaded ensemble (gated on the Tier-3 crop
de-globaling, T1-E).
</content>
