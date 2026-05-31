# src/bindings

## Responsibility

C-ABI facades that expose SWAP to external callers (Python, imod_coupler).

- **`swap_capi_mod.f90`** — singleton C-API backing the standalone BMI library.
- **`swap_bmi_mod.f90`** — BMI lifecycle/time/info over the singleton; drives
  `libswap_bmi.so`.
- **`swap_xmi_mod.f90`** — ensemble-backed coupling C-ABI; drives
  `libswap_xmi.so` (the library imod_coupler / xmipy loads).
- **`bmi_constants_mod.f90`** — BMI string-length and status constants.

## Build note

`swap_capi/bmi` and `swap_xmi` export **overlapping** `bind(C)` names
(`initialize`/`update`/`finalize`/…) and therefore **must not** be linked into
the same shared library. `libswap_bmi.so` is built from `swap_capi`+`swap_bmi`;
`libswap_xmi.so` from `swap_xmi`. See `meson.build` (`modern_sources` vs
`modern_xmi_sources`).

## Dependencies

Wraps `driver/`. Depends on `core/`, `state/`, `config/`.
