---
title: "ADR 0021 — Tillage parameters ported to [soil.tillage] TOML block"
date: 2026-05-06
status: accepted
---

# ADR 0021: Tillage parameters ported to `[soil.tillage]` TOML block

## Context

ADR 0019's "Update 2026-05-06" closed the umbrella physical deletion
of the legacy fixed-format reader, but two TTutil-based reads
survived as foot-noted exceptions: `Read_Tillage` and
`read_ssdi_input`. Both still opened the staged `swap.swp` via
`RDinit` to load parameters that had no TOML schema yet.

This ADR ports the tillage subsystem; SSDI follows in ADR 0022.

## Decision

Tillage parameters move to a new `[soil.tillage]` sub-table. The
`swtill` switch stays at `[soil].swtill` (least disruption to
existing TOML cases). Per-event and per-type tables are inline TOML
arrays-of-tables (`[[soil.tillage.events]]`, `[[soil.tillage.types]]`)
— small, tightly coupled, and benefit from same-file readability.

`Read_Tillage` is deleted. A new helper `apply_soil_tillage` in
`config_to_variables.f90` populates the legacy globals
(`till_Date_tillage`, `till_Z_tillage`, `till_Type_Tillage`,
`till_iType_Tillage`, `till_TAB_Rho_cons`, `till_TAB_Rho_tillage`,
`till_TAB_K_R_cons`, optional `till_TAB_Rho_match` /
`till_TAB_N_match`, plus `till_Ntill`, `till_Ntypes`,
`till_Max_Z_tillage`, `till_iTT1`, `till_iTT2`,
`till_i_n_model`, `till_iRedist`) at config-load time.

## Schema

```toml
[soil]
swtill = 1   # unchanged location

[soil.tillage]
i_n_model = 2   # 1..3
iRedist   = 2   # 0..2

[[soil.tillage.events]]
date      = 2003-04-15
z         = 30.0
intensity = 1.0
type_id   = 1

[[soil.tillage.types]]
id          = 1
rho_cons    = 1500.0
rho_tillage = 1100.0
k_R         = 0.05
# rho_match / N_match required only when i_n_model = 3
```

When `swtill = 0`, the `[soil.tillage]` block can be omitted.

## Consequences

- `swpfile` global pointer no longer read by `tillage.f90`. SSDI
  is the only remaining `swpfile` reader; ADR 0022 will retire it.
- `Read_Tillage` deleted; `apply_soil_tillage` is the sole
  tillage-init path.
- Pattern reusable for SSDI (next ADR).
- `flCropNut` un-stub-erring + nutrient reactivation remains a
  separate arc.

## Tests

pFUnit suites under `tests/unit/io/toml/`:
- `test_read_soil_tillage_toml`
- `test_soil_tillage_validate`
- `test_apply_soil_tillage`

No `swtill = 1` regression case is added in this ADR; deferred
until a known-good legacy comparator is available.

## Acceptance

- `grep -n "subroutine Read_Tillage\b" src/` → no matches.
- `grep -n "swpfile" src/crop/tillage.f90` → no matches.
- `pixi run -e test test-pfunit` → `Ok: 1, Fail: 0` (555 tests,
  1 disabled).
- `pixi run -e test check-full` → `5 passed, 0 failed`.
