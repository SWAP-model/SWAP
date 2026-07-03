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

### Sub-arc 2 — one `libswap.so` on the ensemble backing store  *(next spec)*
The structural move. Port the CAPI accessors + in-memory init + a per-timestep
step onto the ensemble; point the BMI facade at a 1-column ensemble; collapse
`modern_sources`/`modern_xmi_sources` so `swap_bmi_mod` + `swap_xmi_mod` no longer
export colliding names (one lifecycle, ensemble-backed); merge the two
`shared_library` targets into `libswap.so`; update every consumer
(`tests/{bmi,cffi-demo,coupling}`, `prototype/`) to load the one library. Keep the
CLI static. Semantic reconciliation to resolve: BMI `update()` = one timestep vs
XMI `solve()` = one day.

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
