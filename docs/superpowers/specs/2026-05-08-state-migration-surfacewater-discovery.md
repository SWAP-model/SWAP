# Subsystem Migration Discovery: Surface Water

**Date:** 2026-05-08
**Status:** discovery (read-only inventory)
**Pilot for:** state-type migration playbook (first subsystem migration arc)

## 1. Source files

| File | LoC | Role |
|------|-----|------|
| `src/drainage/surfacewater.f90` | 681 | Dispatcher + leaf compute (SurfaceWater task 1/2/3; WLEVBAL; WBALLEV) |
| `src/drainage/surfacewater_init.f90` | 133 | Init — builds sttab, converts L m→cm, zeroes wlsbak/numadj, sets wls/wlstar |
| `src/utils/surfacewaterutils.f90` | 222 | Utility — wlevst, swstlev, qhtab, runoff (pure-reads from module globals) |
| `src/config/surface_water_config.f90` | 196 | Typed config struct + validate + finalize |
| `src/io/toml/read_surface_water_toml.f90` | 201 | TOML reader for `[surface_water]` section |
| `tests/unit/drainage/test_surfacewater_init.pf` | 93 | pFUnit — init state zeroing, sttab geometry, storage row |
| `tests/unit/config/test_surface_water_config.pf` | 355 | pFUnit — config validate, finalize (alphaw normalisation, wldip abs) |
| `tests/unit/io/toml/test_read_surface_water_toml.pf` | — | pFUnit — TOML reader coverage (missing section, full populate, length mismatch) |
| `tests/unit/io/toml/test_surfacewater_parity.pf` | — | pFUnit — literal-parity against case-6 legacy globals |

**Total home-file LoC (four subsystem-home files + TOML reader): 1433**

---

## 2. Owned globals (this subsystem writes)

Variables declared in `src/core/variables.f90` that the subsystem assigns to. Listed with `variables.f90` line, type, inline comment, and write sites.

| Variable | Type | `variables.f90` line | Description | Write sites (file:line) |
|----------|------|----------------------|-------------|--------------------------|
| `wls` | `real(8)` | 1271 | Surface water level (cm) | `surfacewater.f90:448,457,471,480,506,558` (WLEVBAL); `surfacewater.f90:651` (WBALLEV) |
| `wlstar` | `real(8)` | 1263 | Target surface water level (cm) | `surfacewater.f90:91` (init via out-param); `surfacewater.f90:361,406,407,409` (WLEVBAL); `surfacewater_init.f90:91` |
| `swst` | `real(8)` | 1270 | Surface water storage per unit area (cm) | `surfacewater.f90:447,454,470,479,507,559` (WLEVBAL); `surfacewater.f90:655` (WBALLEV); `surfacewater_init.f90:130` |
| `swstini` | `real(8)` | 1270 | Initial surface water storage (cm) | `surfacewater_init.f90:129` |
| `sttab(22,2)` | `real(8)` | 1269 | Level-storage table (col 1 = level cm, col 2 = storage cm) | `surfacewater_init.f90:96,99,108,121,122,123` |
| `wlsbak(4)` | `real(8)` | 1270 | Circular buffer of last 4 surface water levels (oscillation detection) | `surfacewater_init.f90:84,85`; `surfacewater.f90:574,575,576,577` (WLEVBAL) |
| `numadj` | `integer` | 1255 | Counter of target-level adjustments | `surfacewater_init.f90:82`; `surfacewater.f90:415` (WLEVBAL) |
| `hwlman` | `real(8)` | 1272 | Pressure head used for target level computation (cm) | `surfacewater.f90:59` (task=1); `surfacewater.f90:394` (WLEVBAL) |
| `vtair` | `real(8)` | 960 | Total air volume in soil column (cm) | `surfacewater.f90:60` (task=1); `surfacewater.f90:380,382,383` (WLEVBAL) |
| `overfl` | `logical` | 1273 | Flag: automatic weir overflow occurred | `surfacewater.f90:336,508,510` (WLEVBAL) |
| `imper` | `integer` | 1256 | Current surface water management period index | `surfacewater.f90:342,343` (WLEVBAL) |
| `cqdrd` | `real(8)` | 1271 | Cumulative drainage into surface water reservoir (cm) | `surfacewater.f90:103` (task=2 reset); `surfacewater.f90:598,674` (WLEVBAL, WBALLEV) |
| `cwsupp` | `real(8)` | 1271 | Cumulative external supply to reservoir (cm) | `surfacewater.f90:104` (task=2 reset); `surfacewater.f90:599,675` (WLEVBAL, WBALLEV) |
| `cwout` | `real(8)` | 1271 | Cumulative outflow from reservoir (cm) | `surfacewater.f90:105` (task=2 reset); `surfacewater.f90:600,676` (WLEVBAL, WBALLEV) |
| `wlsold` | `real(8)` | 1272 | Surface water level at previous time step (cm) | `surfacewater.f90:648` (WBALLEV) |
| `cqdra` | `real(8)` | 768 | Cumulative lateral drainage, all levels (cm) | `surfacewater.f90:106` (task=2 reset) |
| `cqdrain(Madr)` | `real(8)` | 769 | Cumulative drainage per level (cm) | `surfacewater.f90:108` (task=2 reset) |
| `cqdrainin(Madr)` | `real(8)` | 770 | Cumulative infiltration per level (cm) | `surfacewater.f90:109` (task=2 reset) |
| `cqdrainout(Madr)` | `real(8)` | 771 | Cumulative drainage outflow per level (cm) | `surfacewater.f90:110` (task=2 reset) |
| `inqdra(Madr,macp)` | `real(8)` | 826 | Intermediate lateral drainage per level/compartment (cm) | `surfacewater.f90:93` (task=2 reset) |
| `inqdra_in(Madr,macp)` | `real(8)` | 827 | Intermediate infiltration per level/compartment (cm) | `surfacewater.f90:94` (task=2 reset) |
| `inqdra_out(Madr,macp)` | `real(8)` | 827 | Intermediate drainage outflow per level/compartment (cm) | `surfacewater.f90:95` (task=2 reset) |
| `iqdra` | `real(8)` | 837 | Intermediate lateral drainage total (cm) | `surfacewater.f90:98` (task=2 reset) |
| `qdra(Madr,macp)` | `real(8)` | 877 | Lateral drainage flux per level/compartment (cm/d) | `surfacewater.f90:129,178,181,208` (task=2, partitioning) |
| `qdrtot` | `real(8)` | 881 | Total lateral drainage flux (cm/d) | `surfacewater.f90:212,214` (task=2) |
| `ZDraBas` | `real(8)` | 1212 | Drainage basis level for rapid macropore drainage (cm) | `surfacewater.f90:70,72,74` (task=1, macropore branch) |
| `flInitDraBas` | `logical` | 1214 | Flag: init drainage basis for rapid drainage (macropore) | `surfacewater.f90:77` (task=1, set to `.false.`) |
| `l(Madr)` | `real(8)` | 859 | Drain spacing (converted from m to cm in init) | `surfacewater_init.f90:79` (m→cm conversion) |
| `fldecdt` | `logical` | 93 | Flag: decrease time step | `surfacewater.f90:569,584` (WLEVBAL, ponding/oscillation check) |

**Note on `qdrain(Madr)`:** The subsystem sets `qdrain(level) = 0.0d0` at line 117–120 of `surfacewater.f90` (task=2, when `gwl > 998`). However, `qdrain` is also written by `drainage.f90` (bocodre) — it is a shared variable, owned primarily by drainage but conditionally overridden here. This is a cross-subsystem hazard (see Section 8).

---

## 3. Borrowed globals (this subsystem reads, owned elsewhere)

Variables read by the subsystem but not written by it (or written only by another module).

| Variable | Best-guess owner subsystem | Read sites (file:line) |
|----------|----------------------------|--------------------------|
| `gwl` | water-balance / soil-water | `surfacewater.f90:115,375` (WLEVBAL task=2, WLEVBAL) |
| `dt` | time-control | `surfacewater.f90:224,406,433,467,493,526,545,569,598,599,600` (WLEVBAL); `surfacewater.f90:651,660,674,675,676` (WBALLEV); `surfacewaterutils.f90:13` |
| `t1900` | time-control | `surfacewater.f90:72,138,353` (task=1, task=2); `surfacewater.f90:583,589` (WLEVBAL); `surfacewater.f90:651` (WBALLEV) |
| `T` | time-control | `surfacewater.f90:369` (WLEVBAL) |
| `tcum` | time-control | `surfacewater.f90:372` (WLEVBAL) |
| `flzerointr` | time-control | `surfacewater.f90:90` (task=2) |
| `flzerocumu` | time-control | `surfacewater.f90:102` (task=2) |
| `fldtmin` | time-control | `surfacewater.f90:583` (WLEVBAL) |
| `logf` | output/control | `surfacewater.f90:591` (WLEVBAL) |
| `swscre` | output/control | `surfacewater.f90:591` (WLEVBAL) |
| `NUMNOD` | soil-water init | `surfacewater.f90:381` (WLEVBAL) |
| `THETAS(macp)` | soil-water | `surfacewater.f90:382` (WLEVBAL) |
| `THETA(macp)` | soil-water | `surfacewater.f90:382` (WLEVBAL) |
| `DZ(macp)` | soil-water | `surfacewater.f90:137,382` (task=2, WLEVBAL) |
| `H(macp)` | soil-water | `surfacewater.f90:390,394` (WLEVBAL) |
| `pond` | boundary / water-balance | `surfacewater.f90:568` (WLEVBAL); `surfacewaterutils.f90:13` |
| `pondmx` | boundary / water-balance | `surfacewater.f90:568` (WLEVBAL); `surfacewaterutils.f90:13` |
| `rsro` | boundary / water-balance | `surfacewater.f90:569` (WLEVBAL); `surfacewaterutils.f90:13` |
| `rsroexp` | boundary / water-balance | `surfacewaterutils.f90:13` |
| `runots` | boundary (owned by `boundtop.f90`) | `surfacewater.f90:433,460,468,479,493,526,545,560` (WLEVBAL); `surfacewater.f90:660` (WBALLEV) |
| `QRapDra` | macropore | `surfacewater.f90:433,468,493,526,545` (WLEVBAL); `surfacewater.f90:660` (WBALLEV) |
| `qdrd` | drainage (`drainage.f90` / bocodre) | `surfacewater.f90:433,468,493,526,545,598` (WLEVBAL); `surfacewater.f90:660,674` (WBALLEV) |
| `nrlevs` | soil-water / drainage init | `surfacewater.f90:107,112,127,207,213` (task=2); `surfacewater_init.f90:78,109` |
| `numnod` | soil-water init | `surfacewater.f90:91,127` (task=2) |
| `swdivd` | drainage config | `surfacewater.f90:133` (task=2) |
| `swdislay` | drainage config | `surfacewater.f90:142,150` (task=2) |
| `swtopdislay(Madr)` | drainage config | `surfacewater.f90:144,152,175` (task=2) |
| `zTopDisLay(Madr)` | drainage config | `surfacewater.f90:145,156,161,178` (task=2) |
| `fTopDisLay(Madr)` | drainage config | `surfacewater.f90:145,146` (task=2) |
| `FacDpthInf` | drainage config | `surfacewater.f90:138` (task=2, DIVDRA call) |
| `ksatfit(maho)` | soil-water | `surfacewater.f90:137` (task=2, DIVDRA call) |
| `ksatexm(maho)` | soil-water | `surfacewater.f90:137` (task=2, DIVDRA call) |
| `fluseksatexm(macp)` | soil-water | `surfacewater.f90:137` (task=2, DIVDRA call) |
| `layer(macp)` | soil-water | `surfacewater.f90:137` (task=2, DIVDRA call) |
| `cofani(maho)` | soil-water | `surfacewater.f90:137` (task=2, DIVDRA call) |
| `owltab(Madr,2*maowl)` | surface-water config (legacy reader) | `surfacewater.f90:72,139` (task=1 macropore; task=2 DIVDRA call) |
| `zbotdr(Madr)` | drainage config | `surfacewater.f90:70,71,423,448` (task=1, WLEVBAL); `surfacewater_init.f90:99,110,117` |
| `widthr(Madr)` | drainage config | `surfacewater_init.f90:114,119` |
| `taludr(Madr)` | drainage config | `surfacewater_init.f90:114,119` |
| `swdtyp(Madr)` | drainage config | `surfacewater.f90:69` (task=1); `surfacewater_init.f90:110` |
| `wlstab(2*mawls)` | surface-water config (legacy reader) | `surfacewater.f90:72` (task=1); `surfacewater.f90:651` (WBALLEV) |
| `wlptab(2*mawlp)` | surface-water config (legacy reader) | `surfacewater.f90:224` (task=3, swsrf=3 branch) |
| `wls1_init` | surface-water config (TOML adapter) | `surfacewater_init.f90:89` |
| `NRPRI` | surface-water config (legacy reader) | `surfacewater.f90:423,448` (WLEVBAL) |
| `nmper` | surface-water config | `surfacewater.f90:346` (WLEVBAL); `surfacewater_init.f90:65` |
| `swman(mamp)` | surface-water config | `surfacewater.f90:356,488,513` (WLEVBAL); `surfacewater_init.f90:65,66` |
| `impend(mamp)` | surface-water config | `surfacewater.f90:353` (WLEVBAL) |
| `intwl(mamp)` | surface-water config | `surfacewater.f90:369` (WLEVBAL) |
| `hbweir(mamp)` | surface-water config | `surfacewater.f90:361,498,519,534,541` (WLEVBAL) |
| `alphaw(mamp)` | surface-water config | `surfacewater.f90:499,541` (WLEVBAL) |
| `betaw(mamp)` | surface-water config | `surfacewater.f90:499,541,542` (WLEVBAL) |
| `wscap(mamp)` | surface-water config | `surfacewater.f90:425` (WLEVBAL) |
| `wldip(mamp)` | surface-water config | `surfacewater.f90:422` (WLEVBAL) |
| `wlsman(mamp,mamte)` | surface-water config (legacy reader only, not in `surface_water_config_t`) | `surfacewater.f90:396` (WLEVBAL) |
| `gwlcrit(mamp,mamte)` | surface-water config (legacy reader only) | `surfacewater.f90:375` (WLEVBAL) |
| `nphase(mamp)` | surface-water config (legacy reader only) | `surfacewater.f90:374` (WLEVBAL) |
| `nodhd(mamp)` | surface-water config (legacy reader only) | `surfacewater.f90:390,394` (WLEVBAL) |
| `dropr(mamp*mamte)` | surface-water config (legacy reader only) | `surfacewater.f90:405,406` (WLEVBAL) |
| `VCRIT(mamp,mamte)` | surface-water config (legacy reader only) | `surfacewater.f90:385` (WLEVBAL) |
| `HCRIT(mamp,mamte)` | surface-water config (legacy reader only) | `surfacewater.f90:390` (WLEVBAL) |
| `SWQHR` | surface-water config | `surfacewater.f90:497,518,521` (WLEVBAL); `surfacewater_init.f90:60` |
| `QQHTAB(mamp,mamte)` | surface-water config (legacy reader only) | `surfacewater.f90:522` (WLEVBAL); `surfacewaterutils.f90:13` |
| `hqhtab(mamp,mamte)` | surface-water config (legacy reader only) | `surfacewaterutils.f90:13` |
| `osswlm` | surface-water config | `surfacewater.f90:582` (WLEVBAL) |
| `swsrf` | surface-water config | `surfacewater.f90:223` (task=3); `surfacewater_init.f90:60` |
| `swsec` | surface-water config | `surfacewater.f90:226,229` (task=3); `surfacewater_init.f90:60` |
| `swdra` | drainage config | `surfacewaterutils.f90:13` |
| `sttab(22,2)` | owned by this subsystem (init), read by utils | `surfacewaterutils.f90:13` (module-level import) |
| `imper` | owned by this subsystem (WLEVBAL), read by utils | `surfacewaterutils.f90:13` (module-level import for qhtab) |
| `NumLevRapDra` | macropore | `surfacewater.f90:64` (task=1) |
| `flInitDraBas` | owned by this subsystem, read here too | `surfacewater.f90:63` (task=1 branch check) |
| `WLSOLD` | owned by this subsystem (WBALLEV writes it) | read back in WBALLEV same subroutine |

---

## 4. Entry points (subroutines called from outside this subsystem)

| Subroutine (this subsystem) | Called from (file:line) | Purpose / lifecycle stage |
|-----------------------------|--------------------------|---------------------------|
| `SurfaceWater(1)` | `src/core/swap.f90:174` | Init — read/build drainage tables, zero buffers (init) |
| `SurfaceWater(2)` | `src/core/swap.f90:283` | Compute lateral drainage fluxes and distribute over compartments (per-step) |
| `SurfaceWater(2)` | `src/io/swapoutput.f90:3639` | Same, called from storage/recharge output subroutine (per-step, inside sensitivity loop) |
| `SurfaceWater(3)` | `src/core/swap.f90:290` | Compute surface water balance (WLEVBAL or WBALLEV) (per-step) |
| `SurfaceWater(3)` | `src/io/swapoutput.f90:3651` | Same, from storage/recharge output (per-step, inside sensitivity loop) |
| `runoff()` | `src/boundary/boundtop.f90:22` (module-level use); `boundtop.f90:253,266,281,299` | Calculate surface runoff from ponded water (per-step, in PONDRUNOFF) |

**Note:** `SurfaceWaterOutput` (task 1/2/3) is defined in `src/io/swapoutput.f90` — it is not part of the four subsystem-home files. It calls `outdrf` and `outswb`, which are also in `swapoutput.f90`. Those are output-layer files, not owned by the subsystem itself.

---

## 5. Internal call graph

```
SurfaceWater(task)          [surfacewater.f90:24]
├── (task=1) surfacewater_init(wls, wlp)
│   │   [surfacewater_init.f90:44]
│   ├── → external: fatalerr_collected (guard checks)
│   └── swstlev(wls1)       [surfacewaterutils.f90:93]
│       └── (reads sttab; terminates on bounds error)
│
├── (task=1, macropore branch) afgen(wlstab, …, t1900)
│   └── → external: array_utils/afgen
│
├── (task=2) bocodre(dh)
│   └── → external: drainage_mod/bocodre [drainage.f90:508]
│
├── (task=2, swdivd=1) DIVDRA(…)
│   └── → external: distribute_drainage/DIVDRA
│
├── (task=3, swsrf=3) afgen(wlptab, …, t1900)
│   └── → external: array_utils/afgen
│
├── (task=3, swsec=2) WLEVBAL()
│   │   [surfacewater.f90:241]
│   ├── swstlev(wlstar)     [surfacewaterutils.f90:93]
│   ├── wlevst(swst)        [surfacewaterutils.f90:40]
│   ├── qhtab(wlstar)       [surfacewaterutils.f90:142]
│   │   └── (reads hqhtab, qqhtab, imper)
│   ├── dtdpst(…)           → external: date utility
│   └── warn(…)             → external: output utility
│
└── (task=3, swsec=1) WBALLEV()
        [surfacewater.f90:608]
    ├── AFGEN(WLSTAB, …)    → external: array_utils/afgen
    └── swstlev(wls)        [surfacewaterutils.f90:93]
```

`surfacewater_utils` functions (`wlevst`, `swstlev`, `qhtab`, `runoff`) are leaf routines — they call only `fatalerr_collected` and perform arithmetic on module-level globals. `runoff()` is additionally called from `boundtop.f90` (external to subsystem).

---

## 6. Output coupling

| Output writer (file:subroutine) | Reads which fields (subsystem-owned globals) | Destination file |
|----------------------------------|-----------------------------------------------|-----------------|
| `src/io/swapoutput.f90:outdrf` | `cqdrain(1..nrlevs)`, `cqdrd`, `crunoff` (water-balance owned), `cQMpOutDrRap` (macropore owned) | `*.drf` |
| `src/io/swapoutput.f90:outswb` | `wlstar`, `wls`, `swst`, `swstini`, `cqdrd`, `cwsupp`, `cwout`, `overfl`, `hwlman`, `vtair`, `numadj`, `imper`, `hbweir(imper)` | `*.swb`, `*.man` |
| `src/io/swap_csv_output.f90:csv_out` | `iqdra`, `inqdra(1:nrlevs,:)`, `inqdra_in(1:nrlevs,:)`, `inqdra_out(1:nrlevs,:)` | user-specified CSV (INLIST_CSV) |

**Note on `outswb` writing `swstini`:** `swapoutput.f90:3330` contains `if (daynr.eq.nrOfDays) swstini = swst` inside `outswb`. This is the output layer reaching back into the subsystem's state to update an owned variable as a year-reset operation. This is a bidirectional coupling that needs design-phase resolution (see Section 8, item 4).

---

## 7. Config inputs

Config fields fed from `[surface_water]` TOML → `surface_water_config_t` → adapter in `config_to_variables.f90` (lines 1029–1080).

| TOML key | Typed-config field path | Global variable(s) set (adapter line) | Compute-time read |
|----------|--------------------------|---------------------------------------|-------------------|
| `surface_water.swsrf` | `config%surface_water%swsrf` | `swsrf` (1029) | `surfacewater.f90:223`; `surfacewater_init.f90:60` |
| `surface_water.swsec` | `config%surface_water%swsec` | `swsec` (1030) | `surfacewater.f90:226,229`; `surfacewater_init.f90:60` |
| `surface_water.wlact` | `config%surface_water%wlact` | `wls1_init` via `wlact - altcu` (1036) | `surfacewater_init.f90:89` |
| `surface_water.osswlm` | `config%surface_water%osswlm` | `osswlm` (1037) | `surfacewater.f90:582` (WLEVBAL oscillation threshold) |
| `surface_water.nmper` | `config%surface_water%nmper` | `nmper` (1038) | `surfacewater.f90:346`; `surfacewater_init.f90:65` |
| `surface_water.swqhr` | `config%surface_water%swqhr` | `swqhr` (1039) | `surfacewater.f90:497,518,521`; `surfacewater_init.f90:60` |
| `surface_water.sofcu` | `config%surface_water%sofcu` | (consumed in `finalize`, normalises `alphaw`; not written to a global directly) | — |
| `surface_water.management.impend` | `config%surface_water%impend(:)` | `impend(1:nmper)` (1041–1044) | `surfacewater.f90:353` |
| `surface_water.management.swman` | `config%surface_water%swman(:)` | `swman(1:nmper)` (1046–1049) | `surfacewater.f90:356,488,513`; `surfacewater_init.f90:65` |
| `surface_water.management.wscap` | `config%surface_water%wscap(:)` | `wscap(1:nmper)` (1051–1054) | `surfacewater.f90:425` |
| `surface_water.management.wldip` | `config%surface_water%wldip(:)` | `wldip(1:nmper)` (1056–1059); finalize applies `abs()` before copy | `surfacewater.f90:422` |
| `surface_water.management.intwl` | `config%surface_water%intwl(:)` | `intwl(1:nmper)` (1061–1064) | `surfacewater.f90:369` |
| `surface_water.weir.hbweir` | `config%surface_water%hbweir(:)` | `hbweir(1:nmper)` (1068–1071) | `surfacewater.f90:361,498,519,534,541` |
| `surface_water.weir.alphaw` | `config%surface_water%alphaw(:)` | `alphaw(1:nmper)` (1073–1076); finalize applies `8.64 * 100^(1-betaw) / sofcu` before copy | `surfacewater.f90:499,541` |
| `surface_water.weir.betaw` | `config%surface_water%betaw(:)` | `betaw(1:nmper)` (1078–1081) | `surfacewater.f90:499,541,542` |

**Surface-water config fields NOT yet in `surface_water_config_t` (legacy reader only):**
`wlsman`, `gwlcrit`, `nphase`, `nodhd`, `dropr`, `vcrit`, `hcrit`, `hqhtab`, `qqhtab`, `owltab`, `wlstab`, `wlptab`, `NRPRI`, `nrsec`. These are read by WLEVBAL and are required for `swman=2` (automatic weir), `swsec=1` (input water level), and `swqhr=2` (QH-table discharge) — all currently stub-errored on the TOML path.

---

## 8. Open questions / cross-subsystem hazards

1. **`qdrain(Madr)` dual ownership** (`surfacewater.f90:117`; `drainage.f90:178,187,216,231,247,280,282,283,296,297`). `bocodre` in `drainage.f90` is the primary writer of `qdrain`. The surfacewater subsystem conditionally overrides it with zeros (`qdrain(level) = 0.0d0` when `gwl > 998`). A clean `surfacewater_state_t` cannot own `qdrain` exclusively; it needs to be in a shared "drainage-flux state" type or the override must be restructured. Needs design-phase resolution.

2. **`runoff()` in `surfacewater_utils` is called by `boundtop.f90`** (`boundtop.f90:253,266,281,299`). The function reads `swdra`, `pond`, `pondmx`, `rsro`, `rsroexp`, `wls`, `swst`, `dt` via module-level `use variables`. After migration, if `wls` and `swst` move into `surfacewater_state_t`, the function signature must change (threaded argument or accessor). This coupling crosses the subsystem boundary into the boundary/top condition subsystem.

3. **`imper` is both owned (WLEVBAL writes it) and re-read by `surfacewater_utils:qhtab`** via the module-level `use variables`. The current module-level import in `surfacewaterutils.f90:13` (`use variables, only: sttab, imper, …`) means the utility functions are statefully coupled to the global `imper`. Post-migration, `qhtab` must receive `imper` as an argument or `surfacewater_utils` must import from `surfacewater_state_t`.

4. **`outswb` (output layer) writes `swstini`** (`swapoutput.f90:3330`): `if (daynr.eq.nrOfDays) swstini = swst`. This is the output layer mutating a subsystem-owned state variable as a year-reset. After migration, the output layer must not write directly into `surfacewater_state_t`; the year-reset should be triggered via a callback or an explicit reset entry point in the subsystem.

5. **`SurfaceWater(2)` is called from `swapoutput.f90:3639`** inside a storage/recharge perturbation loop (for computing storage-discharge sensitivity). This is a compute call from within an output file — an architectural smell that will need consideration when deciding what signature `SurfaceWater(task)` acquires post-migration.

6. **`vtair` is declared in the soil-water section of `variables.f90` (line 960)** but written exclusively by WLEVBAL (`surfacewater.f90:380,382`). The inline comment says "Total depth of air in the soil column (L)" — a soil-water concept, but the surface water subsystem computes it for the purpose of target-level determination. Ownership needs explicit assignment in the design phase.

7. **`fldecdt` is both a time-control variable and written by WLEVBAL** (`surfacewater.f90:569,584`). The timestep-reduction flag is set when surface water oscillations or ponding exceed thresholds. Post-migration, the subsystem state type must be able to signal the main loop to reduce timestep, either via a flag field in the state type or a returned logical.

8. **`wlsman`, `gwlcrit`, `nphase`, `nodhd`, `dropr`, `vcrit`, `hcrit`, `hqhtab`, `qqhtab` are NOT in `surface_water_config_t`** and are only populated by the legacy `readswap`/`rddre` reader. These are required for `swman=2`, `swsec=1`, and `swqhr=2` branches — all stub-errored on the TOML path. Before those branches are enabled, these fields must be added to `surface_water_config_t`. This is a pre-requisite for full migration.

9. **`wls` is read by `boundtop.f90`** indirectly through `surfacewater_utils:runoff()` (which uses `wls` from module-level `use variables`). Water balance and surface runoff both depend on the current surface water level. This is a mutual coupling: `wls` → `runoff()` → `runots` → WLEVBAL reads `runots`. The dependency chain `wls → runots → wls` must be checked for ordering hazards in the migration.

10. **`l(Madr)` unit conversion in `surfacewater_init`** (`surfacewater_init.f90:79`): The init module converts `l(i)` from metres to centimetres in-place on the global array. This is a destructive side-effect on a drainage-owned global, guarded by the module comment (only for `dramet=0`/`swdra=2` paths). After migration, this conversion must not mutate a shared global; it should happen either at config-load time or within the init by reading from a typed config field.

---

## 3.5 External readers of owned globals

**Added retroactively at end of Phase 2** (lessons-learned from execution: the original discovery cataloged owned writes but missed external reads; future subsystem migrations populate this section during discovery, not during execution).

For each surface-water-owned global, the files OUTSIDE this subsystem's home tree that read it. This determines the scope of cross-subsystem migration work in Phase 2.

| Owned global | External reader files |
|---|---|
| `wls` | `src/drainage/drainage.f90` (via `bocodre`); also indirectly via `surfacewater_utils:runoff()` from `boundtop.f90` |
| `swst` | `src/drainage/drainage.f90`; also indirectly via `runoff()` |
| `imper` | `src/drainage/drainage.f90`; `src/utils/surfacewaterutils.f90:qhtab` |
| `cqdra` | `src/drainage/drainage.f90`, `src/soil/waterbalance.f90:integral` |
| `ZDraBas` | `src/drainage/drainage.f90`, `src/macropore/macrorate.f90` (`RAPIDDRAIN`), `src/macropore/macropore.f90` (`MACROINIT`) |
| `iqdra` | `src/drainage/drainage.f90`, `src/soil/waterbalance.f90:integral`, `src/io/swap_csv_output.f90:set_values`, `src/io/swapoutput.f90:outwba`/`outinc`/`OutputModflow`/`csv_write` |
| `qdrtot` | `src/heat/frozencond.f90:FrozenBounds`, `src/soil/waterbalance.f90:integral`+`fluxes`, `src/drainage/drainage.f90`, `src/solute/solute.f90` |
| `flInitDraBas` | `src/drainage/drainage.f90`, `src/macropore/macropore.f90:MACROINIT` |
| `sttab` | `src/utils/surfacewaterutils.f90:wlevst`/`swstlev`; `src/config/drainage_config.f90` (comment-only, safe) |
| `cqdrain(Madr)` | `src/soil/waterbalance.f90:integral`, `src/drainage/drainage.f90`, `src/io/swapoutput.f90:outdrf`/`outbal`/`outblc` |
| `cqdrainin(Madr)` | `src/drainage/drainage.f90`, `src/soil/waterbalance.f90:integral`, `src/io/swapoutput.f90:outblc` |
| `cqdrainout(Madr)` | `src/drainage/drainage.f90`, `src/soil/waterbalance.f90:integral`, `src/io/swapoutput.f90:outblc` |
| `qdra(Madr,macp)` | `src/heat/frozencond.f90:FrozenBounds`, `src/drainage/drainage.f90` (passed to `divdra` as explicit-shape), `src/soil/waterbalance.f90:integral`+`fluxes`, `src/soil/soilhydraulics.f90:headcalc`, `src/solute/solute.f90` |
| `inqdra(Madr,macp)` | `src/drainage/drainage.f90`, `src/soil/soilgrid.f90:ConvertDiscrVert`, `src/soil/waterbalance.f90:integral`, `src/crop/management_soil.f90:SoilManagement` |
| `inqdra_in(Madr,macp)` | `src/drainage/drainage.f90`, `src/soil/waterbalance.f90:integral` |
| `inqdra_out(Madr,macp)` | `src/drainage/drainage.f90`, `src/soil/waterbalance.f90:integral` |
| `cqdrd`, `cwsupp`, `cwout` | `src/io/swapoutput.f90:outdrf` only (output-side, migrated in Phase 1) |
| `wlstar`, `wlsold`, `wlsbak`, `numadj`, `vtair`, `hwlman`, `swstini`, `qdrtot`, `overfl` | (no live external readers in Phase 1; owned-only or already migrated) |

Plus the special case (not in surfacewater_state_t):
| `fldecdt` | `src/core/swap.f90`, `src/core/timecontrol.f90`, `src/soil/soilhydraulics.f90:headcalc`, `src/io/swapoutput.f90:OutputModflow` |

**Implication for migration scope:** Phase 2 had to migrate 13 cross-subsystem files in addition to the home tree. Future subsystem-migration discoveries should populate this section during the read-only investigation, BEFORE writing the design doc, so the design phase can plan owner-rule relocations and reader migrations explicitly.

---

## 9. Test surface

| Test file | What it covers |
|-----------|----------------|
| `tests/unit/drainage/test_surfacewater_init.pf` | `surfacewater_init`: state zeroing (numadj, wlsbak), sttab depth rows, sttab storage (trapezoid geometry), wls1/wlp1 out-params |
| `tests/unit/config/test_surface_water_config.pf` | `surface_water_config_t` validate: swsrf enum, swsec=1 rejected, swsec=2 minimal valid, nmper range, array size mismatch, swman invalid enum, swqhr=2 rejected, swman=2 rejected, swsrf=3 rejected; finalize: alphaw normalisation, wldip abs |
| `tests/unit/io/toml/test_read_surface_water_toml.pf` | `read_surface_water_toml`: missing section silent, switch-only configs, full populate (swsec=2 + swqhr=1), per-period array length mismatch, impend date decoding |
| `tests/unit/io/toml/test_surfacewater_parity.pf` | Literal-parity of `config%surface_water` fields against case-6 legacy globals (swsrf, swsec, wlact, osswlm, nmper, swqhr, sofcu, impend[1..28], swman, wscap, wldip, intwl, hbweir, alphaw, betaw) |
| `tests/regression/surfacewater_expected.json` + `surfacewater_expected_gfortran.json` | End-to-end regression: GWL and POND at year-end for 1997/1998/1999 from case-6 (Wildenborch) |
| `tests/swap-cases/toml/6.surfacewater/` | Full integration case: swdra=2, swsrf=2, swsec=2, swqhr=1, 28 management periods, 2 drainage levels |

**Coverage gap:** There are no pFUnit tests for the compute routines WLEVBAL and WBALLEV in isolation, or for `surfacewater_utils` (wlevst, swstlev, qhtab, runoff). Coverage of the compute path relies entirely on the end-to-end regression case.

---

## 10. Summary statistics

- **Total LoC across four home files + TOML reader:** 1433 (681 + 133 + 222 + 196 + 201)
- **Owned globals (Section 2):** 29 variables
- **Borrowed globals (Section 3):** 54 variables
- **Entry points (Section 4):** 6 call sites (2 unique public subroutines: `SurfaceWater(task)` and `runoff()`)
- **Output coupling sites (Section 6):** 3 writers (`outdrf`, `outswb`, `csv_out`)
