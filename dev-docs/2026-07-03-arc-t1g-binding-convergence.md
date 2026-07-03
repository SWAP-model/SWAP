---
title: "Arc T1-G′ — driver / binding convergence (spec + plan)"
date: 2026-07-03
status: in-progress
tags: [tier1, bindings, bmi, xmi, capi, registry, coupling]
parent: 2026-07-03-modernization-review.md
---

# Arc T1-G′ — driver / binding convergence

From the [modernization review](2026-07-03-modernization-review.md) §7 and the
driver-convergence discussion. Goal: move SWAP toward the MODFLOW 6 shape — **one
shared library exposing BMI + XMI over one backing store, driven by a single
variable registry**, with the thin CLI kept static-linked.

## Where we are today (verified)

Three drive paths diverge at the *backing store*, not the engine (all bottom out
in `swap_init`/`swap_run_step`/`swap_close` on an explicit `(state, config)`):

| Path | Library | Backs onto | Variable vocabulary |
|---|---|---|---|
| CLI | `swap` (static) | local `state,config` (`swap_main.f90`) | — |
| BMI + CAPI | `libswap_bmi.so` | singleton `capi_state` | CSDMS (`soil_water_content`) + raw (`theta`) |
| XMI | `libswap_xmi.so` | ensemble `columns(:)` | exchange (`gwl`,`qbot_volume`) |

Two shared libraries exist **only** because `swap_bmi_mod` and `swap_xmi_mod` both
export the C symbols `initialize`/`update`/`finalize`/… backed by different stores
(`meson.build:65,311,333`). Costs: `c_to_f_string`/`f_to_c_string` duplicated ×3;
per-variable `select case` switches hand-maintained in each facade with no
cross-check; three vocabularies.

The ensemble is day-oriented and coupling-specific: `ensemble_step_day()` only
(no per-timestep step), file+`ensemble.txt` init only (no in-memory TOML), and
none of the CAPI rich accessors (results view, water balance, meteo buffer,
in-memory config). So making the ensemble the sole backing store is a real port.

## Decomposition — two shippable sub-arcs

The two goals are separable, and the registry de-risks the merge, so:

### Sub-arc 1 — variable registry (BMI + CAPI array views)  *(this session)*
Behavior-preserving substrate; **no library or backing-store change**; all
existing binding tests stay green. Two facades on one registry from day one, so
the single-source-of-truth is demonstrated (not just set up).

- New `swap_var_registry_mod`: a `var_registry_t` holding, per `real64` variable,
  `{name, units, is_array, readable, settable, scale, pointer-into-state}`, plus
  generic accessors (find / read-into-`dest` / c-pointer+size / units / nbytes /
  packed name-list).
- `build_variable_registry(reg, state)` registers, binding pointers into the
  (post-init) state arrays:
  - the **10 BMI variables** (8 readable + 2 settable BCs) under their CSDMS
    names — encoding the exact current semantics: array vs scalar, the
    `recharge = -qbot` sign flip (`scale = -1`), the `groundwater_level_imposed →
    h(numnod)` BC target, the `dest(m+1:n)=0` fill;
  - the **4 CAPI array views** (`theta`, `h`, `tsoil`, `inqrot`) under their raw
    names (aliases onto the same fields where they overlap).
- `swap_bmi_mod` variable + metadata accessors (`get_value_double`,
  `set_value_double`, `get_var_units`, `get_var_nbytes`, `get_output_var_names`,
  `get_input_var_names`, item counts) and `swap_capi_mod`'s `swap_view_array`
  become registry-driven; the hand-written switches are deleted. Registry is built
  after `swap_init` in both singleton init paths (`bmi_initialize`, CAPI
  `swap_initialize_from_toml_string`).
- **Deliberately out of sub-arc 1** (sub-arc 2, when the libraries merge and it
  pays off): the C-string-helper dedup — the three `f_to_c_string` copies have
  *different* contracts (BMI bounded by `max_len`, XMI unbounded), so merging them
  is a behavior change, not a move; `swap_get_scalar` (heterogeneous int/logical
  types that don't fit a `real64`-pointer registry without type tags); and the
  XMI `get_value_ptr` facade.

**Rationale for stopping here first:** byte-behavior-preserving and test-covered
(`get_value_double` via bmi, `swap_view_array` via cffi-demo), so safe to land
immediately; and it forces the per-variable contract into one place — the
precondition for pointing a merged library at it in sub-arc 2 without re-deriving
semantics.

### Sub-arc 2 — one `libswap.so` on the ensemble backing store  *(in progress)*
The structural move, decomposed into verifiable steps. **The full test path is
runnable in this environment** — `libmf6.so` + `mf6` + flopy are present, so the
SWAP↔MODFLOW6 coupling smoke test runs — so every step below is verifiable across
BMI + CAPI + XMI + coupling + check-fast + pFUnit (a major de-risk vs. the earlier
assumption the coupling path was blind).

- **2a — shared C-string helpers (DONE 2026-07-03).** `swap_c_strings_mod`
  (modern_core): one `c_to_f_string` + a generic `f_to_c_string` resolving BMI's
  bounded (3-arg) and XMI's unbounded (2-arg) forms by argument count, so every
  call site is unchanged. The three per-facade copies deleted. Verified across all
  six paths byte-behavior-preserving.
- **2b — ensemble gains uncoupled single-column mode + in-memory init + a
  per-substep step.** Additive extension of `swap_ensemble_mod` (an `ncol=1`,
  `flcoupled_gwl=.false.` path that skips the gwl-injection/re-equilibration that
  today's `ensemble_step_day` assumes a MODFLOW driver supplies; a `swap_run_step`
  substep entry). Existing coupled path untouched.
- **2c — port CAPI accessors + registry onto ensemble column 1.** The results
  view, water balance, `swap_view_array`, in-memory config attach, meteo buffer,
  and the registry build operate on `columns(1)`.
- **2d — retire the singleton; unify the colliding lifecycle symbols.**
  *(Revised 2026-07-03 — the earlier "BMI `update()` = one day" proposal was a
  behavior change dressed as a default and is dropped.)* The unified lifecycle is
  **mode-dependent**: `initialize()` detects the `ensemble.txt` sidecar — present
  → coupled ensemble mode, lifecycle behaves exactly as XMI today (`update()` =
  one day via `ensemble_step_day`, `get_time_step()` = 1.0, item counts 1/2);
  absent → single-column mode, lifecycle behaves exactly as BMI today
  (`update()` = one Richards substep, `get_time_step()` = dt, registry-driven
  counts). Both existing consumer populations see bit-for-bit unchanged behavior.
  Storage unifies via pointer indirection, not a port: `capi_state`/`capi_config`
  become pointers bound to `columns(1)`/`configs(1)` (ensemble arrays gain
  `target`), so every CAPI accessor compiles unchanged while the ensemble becomes
  the sole owner. Signature reconciliation for the 3 colliding metadata calls
  whose BMI/XMI arities differ (`get_var_type`, `get_component_name`): adopt the
  2-arg unbounded XMI form — xmipy is their only caller; no BMI consumer calls
  them. `swap_set_headless` before init (orchestrator's call order) is buffered
  and applied at init, since the pointer is unassociated pre-init.
- **2e — merge the meson targets → `libswap.so`; update consumers.** Collapse
  `modern_sources`/`modern_xmi_sources`; one `shared_library('swap')`; update the
  hardcoded `libswap_bmi.so`/`libswap_xmi.so` paths in `prototype/`,
  `tests/coupling/_defaults.py`, `imod_coupler.toml`, and the meson test
  registrations. CLI stays static (built from the same core).

## Sub-arc 2 outcome (2026-07-03) — DONE (2b–2e landed together)

One `libswap.so` now exposes BMI + CAPI + XMI; `libswap_bmi.so`/`libswap_xmi.so`
are gone. Implementation as revised:

- **Ensemble is the sole storage owner.** `swap_ensemble_mod` gained an
  uncoupled single-column mode (`ensemble_allocate_single` /
  `ensemble_init_single`), `target` columns/configs, pointer getters, the
  `ensemble.txt` sidecar reader, and a facade-shared `ensemble_last_error`.
- **The singleton became pointers.** `capi_state`/`capi_config` are pointers
  bound to `columns(1)`/`configs(1)` (`capi_bind_first_column`), so every CAPI
  accessor compiled unchanged. `swap_set_headless` before init is buffered and
  applied at bind (the orchestrator's call order).
- **Mode-aware lifecycle in one definer.** `swap_bmi_mod` holds the 14
  formerly-colliding C names: `initialize()` detects the sidecar (present →
  coupled ensemble init; absent → single-column); `update`/`update_until`/
  `get_time_step`/item-counts/`get_var_nbytes` branch on mode, reproducing the
  former BMI and XMI behaviors exactly. `get_var_type`/`get_component_name`
  adopt the 2-/1-arg unbounded XMI arity (verified from xmipy source: its
  `get_value_ptr` calls `get_var_type` with 2 args; no BMI consumer calls
  either). `initialize` keeps the unused-`n` 2-arg form (xmipy passes 1 arg;
  `n` is never read). `swap_xmi_mod` shrank to the XMI-specific verbs only.
- **Code-verified consumer contracts** (per "don't trust the docs"): the
  coupled driver's full call set on the SWAP library is initialize /
  prepare_time_step / prepare_solve / solve / finalize_solve /
  finalize_time_step / finalize / get_version / get_value_ptr(+rank/type/
  shape); `report_timing_totals` is an xmipy python-side timer (never a lib
  symbol); the coupled driver does NOT call `get_current_time`/`get_time_step`
  on SWAP (an earlier agent report claimed otherwise — wrong).

**Verification (all six paths, merged library):** BMI hello_swap and CAPI
run_ensemble bit-identical values (theta[0]=0.2712560347028374, dt=0.0002 —
substep semantics preserved); XMI smoke qbot_volume bit-identical; SWAP↔MODFLOW6
coupled smoke water table bit-identical; check-fast 4/4 byte-identical; pFUnit
838. Remaining for a later arc (recorded in ADR 0050): re-home XMI
`get_value_ptr` onto the registry (NS_XMI), handle-based multi-instance.

## Plan for sub-arc 1 (each step ends with the binding tests + check-fast green)

1. `swap_c_strings_mod` + switch all three facades to it (dedup). Build; run
   `test-bmi`, `test-cffi-demo`, `test-xmi`, `test-coupling`.
2. `swap_var_registry_mod` (type + accessors) + `build_bmi_registry`. Add
   `capi_registry` singleton in `swap_capi_mod`; build it in both init paths.
3. Convert `swap_bmi_mod` accessors to registry-driven; delete the switches.
   Clean rebuild (facade/state-adjacent) + full binding-test sweep + `check-fast`.
4. Commit.

## Gate commands
`pixi run -e test test-bmi test-cffi-demo test-xmi test-coupling check-fast`
(bmi/cffi exercise the singleton BMI/CAPI surface; xmi/coupling exercise the XMI
library — both must stay green since sub-arc 1 leaves both libraries in place).

## Sub-arc 1 outcome (2026-07-03) — DONE

Landed `src/bindings/swap_var_registry_mod.f90` (namespaced `var_registry_t` +
`build_variable_registry`); `swap_bmi_mod` variable/metadata accessors and
`swap_capi_mod` `swap_view_array` are registry-driven (the hand-written switches
deleted — `swap_bmi_mod` −132/+66 lines). Registry built after `swap_init` in both
singleton init paths. Added 5 pFUnit tests (`tests/unit/bindings/test_var_registry.pf`).

**Verification** (the meson `bmi`/`cffi-demo` suites can't run here — they point at
the un-checked-out `tests/swap-cases` submodule — so the two scripts were driven
directly against the self-contained `hupselbrook` regression case):
- BMI `hello_swap.py` (registry `get_value_double` + unknown-name reject): OK.
- CAPI `run_ensemble.py` (registry `swap_view_array`, in-memory init): OK.
- `check-fast`: 4/4 byte-identical (exe unaffected).
- XMI smoke (`libswap_xmi.so`, unchanged facade, rebuilt with the registry in
  `modern_core`): OK.
- pFUnit: 833 → **838**.

Byte-behavior preserved throughout. Two consumers (BMI + CAPI) now share one
registry — the substrate for sub-arc 2 (one `libswap.so`, XMI onto the registry,
backing-store unification).
</content>
