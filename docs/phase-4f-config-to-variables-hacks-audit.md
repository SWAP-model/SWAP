# Phase 4f-extend HACK Slots Audit — `config_to_variables.f90`

**Date:** 2026-05-03
**File audited:** `src/io/toml/config_to_variables.f90`
**Branch:** development (post-Phase 4f .crp port, all 6 regression cases bit-identical)

**Update 2026-05-04:** Bucket A slots (lines 294, 303, 330, 1230) resolved in
commit `8f4f5a3` (feat: bucket-A schema additions — rsigni, cfevappond, outfil, drfil).

**Update 2026-05-03:** Bucket B slots 5/6/7 resolved:
- `rdmax` (line 1144) — commit `a85faa8`
- `swliminf` (line 484) — commit `66559f2`
- `swredu` (line 280) — commit `ef7447d`

---

## Summary

12 real `HACK Phase 4f-extend` slots were catalogued at the original audit. 4
have since been resolved (Bucket A → Bucket E). 8 remain open.

| Bucket | Label | Count (original) | Count (remaining) |
|--------|-------|-----------------|-------------------|
| A | Trivial schema add | 4 | 0 (all resolved) |
| B | Schema add + cross-section validator | 2 | 0 (slots 5/6/7 resolved 2026-05-03) |
| C | Schema redesign / non-trivial | 1 | 1 |
| D | Will-be-irrelevant after config-passing refactor (ADR 0016 Part C) | 3 | 3 |
| E | Already-resolved (stale comment) | 2 | 9 |

Two entries (line 1236 / CSV block, line 1156 / rdmax) are
partially resolved: a schema field exists but hardcoded companion values
remain; they are graded B and B respectively on that basis.

---

## Full Slot Table

| # | Line | Field / Global | Bucket | What's hacked | What's needed | Affected case(s) | Effort |
|---|------|----------------|--------|---------------|---------------|------------------|--------|
| 1 | 280 | `swredu` | E | ~~Inferred from `cofredbo != 0.35`~~ **Resolved in `ef7447d`.** Added `meteorology.evaporation.swredu` (integer, default 1) to `meteorology_evaporation_t`; parser + validator (enum [1,2] + cross-field swredu=2 requires cofredbo≠0.35) + adapter wired. Case 5 authors `swredu = 2`. | — | — | — |
| 2 | 294 | `rsigni` | E | ~~Hardcoded `0.5` (cm)~~ **Resolved in `8f4f5a3`.** Added `meteorology.evaporation.rsigni` (real, default 0.5) to `meteo_evaporation_config_t`; parser + validator + adapter wired. | — | — | — |
| 3 | 303 | `cfevappond` | E | ~~Hardcoded `1.25`~~ **Resolved in `8f4f5a3`.** Added `meteorology.evaporation.cfevappond` (real, default 1.25) to `meteo_evaporation_config_t`; parser + validator + adapter wired. | — | — | — |
| 4 | 330 | `drfil` | E | ~~Hardcoded `'swap'`~~ **Resolved in `8f4f5a3`.** Added `drainage.drfil` (string, default `'swap'`) to `drainage_config_t`; parser + adapter wired. | — | — | — |
| 5 | 484 | `swliminf` | E | ~~Hardcoded `1` when `dramet=3`~~ **Resolved in `66559f2`.** Added `drainage.swliminf` (integer, default 0) to `drainage_config_t`; parser + validator (enum [0,1] + cross-field swliminf=1 requires dramet=3) + adapter wired. Cases 2, 3, 4, 5 (all dramet=3) author `swliminf = 1`. | — | — | — |
| 6 | 716 | `iHWCKmodel(:)` | C | Fixed to `1` (uni-modal van Genuchten) for all soil layers. Legacy reader allows per-layer override (values 1–11), but no schema slot exists. | Requires a per-layer field (array) in `soil_hydraulics_config_t`, schema validation, and potentially a lookup table / enum for the 11 model codes — this is a multi-model hydraulic dispatch, not a simple scalar. | No current regression case authors a value other than 1. Risk is latent. | L |
| 7 | 743 | `paramvg(10,:)` / `ksatexm` threshold branch | B | Always sets `paramvg(10,i) = -999.0` (sentinel = no threshold-Ksat). The legacy path computes `relsatthr`/`ksatthr` when `ksatexm(i) > ksatfit(i)`. All current regression cases author `ksatexm == ksatfit`, so `flksatexm` would be `.false.` in the legacy path — the sentinel is therefore correct for all current cases, but semantically the schema already reads `ksatexm` without wiring the threshold logic. | Port the `flksatexm` branch (readswap.f90:802–815): when `ksatexm(i) > ksatfit(i)`, compute `relsatthr(i)`, `ksatthr(i)` and set `paramvg(10,i) = ksatexm(i)`. Update `cofgen(11,:)` / `cofgen(12,:)`. | No current case triggers this (all have `ksatexm == ksatfit`). Future cases with structured-Ksat reduction need this. | M |
| 8 | 1144 | `rdmax` (`RDS`) | E | ~~Hardcoded `200.0` cm~~ **Resolved in `a85faa8`.** Added `crop.rdmax` (real64, default 200.0 cm) to `crop_config_t`; parser + validator (range [1,5000] cm) + adapter wired. Cases 3 (320.0), 5 (100.0), 6 (60.0) author non-default values; cases 1/2/4 use the default. | — | — | — |
| 9 | 1191 | `cropfil(:)` suffix stripping | D | Strips `.crp.toml` / `.toml` suffixes from `rotation_file(i)` before writing to legacy `cropfil(i)`, because the legacy per-rotation readers (still called via `cropgrowth.f90 ArableLandGerm`) expect a bare stem. | This hack goes away when ADR 0016 Part C lands: once all three crop-mode readers (`readcropfixed`, `readwofost`, `readgrass`) are fully retired and replaced by `read_crop*_toml` init functions, `cropfil` can pass through unchanged (or be eliminated). Until then this is load-bearing. | All TOML crop cases (1, 2, 4, 5, 6). Must remain until legacy reader retirement is complete. | — (defer) |
| 10 | 1230 | `outfil` | E | ~~Hardcoded `'result'`~~ **Resolved in `8f4f5a3`.** Added `general.outfil` (string, default `'result'`) to `general_config_t`; parser + adapter wired. | — | — | — |
| 11 | 1236 | `swcsv` + `InList_csv` / `swcsv_tz` / `InList_csv_tz` | B | `swcsv` is hardcoded to 1 (CSV always on). `InList_csv` is partially resolved — `general.inlist_csv` schema field exists and is used when authored; otherwise falls back to a hardcoded default. `swcsv_tz` is hardcoded to 0 and `InList_csv_tz` to `'wc,h,conc'`. | The partial resolution of `InList_csv` (already done) leaves 3 remaining gaps: (a) `swcsv` (the enable switch, always 1); (b) `swcsv_tz` (depth-profile CSV switch, always 0); (c) `InList_csv_tz` (depth-profile column list, hardcoded). These belong in a `[output.csv]` section. Cross-section concern: the output section does not yet exist in the schema. | All 6 TOML cases are affected — they silently inherit `swcsv=1`. Cases 2, 4, 5, 6 override `inlist_csv`; none yet need `swcsv_tz=1`. | M |
| 12 | 1260 | `swpfile` + `logf` (legacy I/O state) | D | Sets `swpfile = 'swap.swp'` and opens `swap_swap.log` so that unported per-crop / drainage legacy readers (`RDinit(unit, logf, swpfile)`) can open the right input file. | This entire block goes away reader-by-reader as Phase 4f retires the legacy ASCII readers. Fixing it in isolation (e.g. making `swpfile` a typed-config field) would be pointless — the fix is to retire the readers that need it. | All TOML crop cases currently call `ArableLandGerm` / crop-mode readers that depend on this. | — (defer) |

---

## Notes on Specific Slots

### Slot 7 (line 743, ksatexm threshold)
All 6 active regression cases author `ksatexm == ksatfit` element-wise, making
`flksatexm` logically false in the legacy path. The `-999` sentinel is therefore
bit-identical to what the legacy path would produce. This is not Bucket E
(already-resolved) because the underlying branch is genuinely unported and would
produce wrong `cofgen(11,:)` / `cofgen(12,:)` for any case with
`ksatexm > ksatfit`. Bucket B is assigned because the port is straightforward
once the schema already reads `ksatexm`.

### Slot 8 (line 1156, rdmax)
Hardcoded `200.0` passes all 6 cases because crop configs keep `rdc`/`rdtb`
within the per-case soil profile depth. However, the actual per-case `RDS`
values are: 200 (cases 1/2/4), 100 (case 5), 60 (case 6), 320 (case 3 —
excluded). Cases 5 and 6 are latently wrong; any future case with
`rdtb_max > min(rds, 200)` would silently compute incorrect rooting depths.

### Slot 11 (line 1236, CSV output)
Partially resolved: `inlist_csv` has a schema slot and is correctly dispatched.
Only `swcsv` (hardcoded 1), `swcsv_tz` (hardcoded 0), and `InList_csv_tz`
(hardcoded `'wc,h,conc'`) remain open. Graded B (not A) because the output
section does not yet exist in the schema — adding `[output.csv]` requires
defining the section type and registering it in the top-level `swap_config_t`.

### Slots 9 and 12 (lines 1191 and 1260)
Both are intrinsically tied to the legacy-reader-retirement roadmap. Attempting
to fix either in isolation would require retaining the complexity but in a new
location. They should be closed as a side-effect of retiring
`readcropfixed` / `readwofost` / `readgrass` / `rddre` from `cropgrowth.f90`.

---

## Recommended Order

### Quick wins (Bucket A) — COMPLETE
All four Bucket A slots resolved in `8f4f5a3`:
- `rsigni` (line 294) → `meteorology.evaporation.rsigni`
- `cfevappond` (line 303) → `meteorology.evaporation.cfevappond`
- `outfil` (line 1230) → `general.outfil`
- `drfil` (line 330) → `drainage.drfil`

### Medium priority (Bucket B) — COMPLETE
All three Bucket B slots resolved in commits `a85faa8` / `66559f2` / `ef7447d`:
- `rdmax` (line 1144) → `crop.rdmax`
- `swliminf` (line 484) → `drainage.swliminf`
- `swredu` (line 280) → `meteorology.evaporation.swredu`

### Remaining medium priority (Bucket B, ~M)
8. **CSV output block (line 1236)** — Add `[output.csv]` section type to
   `general_config_t`. Three sub-fields: `swcsv` (int 0/1), `swcsv_tz`
   (int 0/1), `InList_csv_tz` (string). All 6 case TOMLs need the section.
   Dependency: schema section type must exist before the parser and adapter
   can be wired.

### Defer (Bucket C and D)

- **`iHWCKmodel` (line 716, Bucket C)** — No regression case triggers
  non-uni-modal hydraulics. Defer until a case requiring models 2–11 is
  added to the test suite. Needs a spec for the 11-model dispatch.

- **`cropfil` suffix stripping (line 1191, Bucket D)** — Load-bearing
  until the legacy crop readers are fully retired (ADR 0016 Part C). Do
  not address in isolation.

- **`swpfile` + `logf` (line 1260, Bucket D)** — Goes away reader-by-reader
  as legacy ASCII readers are retired. No value in fixing in isolation.
