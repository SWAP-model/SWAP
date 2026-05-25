# readcropfixed Line-by-Line Audit

**File:** `src/io/readswap.f90`, lines 2037–2517  
**Subroutine:** `readcropfixed(icrop, crpfil, lcc, swhydrlift)`  
**Date of audit:** 2026-05-01  
**Purpose:** Classification of every non-blank, non-comment line in the ~480-line `readcropfixed` subroutine to inform the `.crp` → TOML port (Tasks 2–7). The subroutine reads the legacy `*.crp` ASCII file, validates its contents, applies hardcoded defaults, and builds runtime data structures (notably the cumulative root density function). After the port, it is replaced by:

- Config validators (VALIDATE-bucket lines)
- Config finalizers / defaults (NORMALIZE-bucket lines)
- A new `cropfixed_init` module (RUNTIME-bucket lines)
- Stub/error checks for unimplemented branches (GUARDED-bucket lines)

Lines that fall into multiple buckets appear in multiple rows — this is intentional and useful.

---

## Classification Table

| Lines         | Bucket    | Notes                                                                                                                |
| ------------- | --------- | -------------------------------------------------------------------------------------------------------------------- |
| 2037          | READ      | Subroutine signature `readcropfixed(icrop,crpfil,lcc,swhydrlift)`                                                   |
| 2044-2053     | READ      | `use variables, only: ...` — imports all module-level crop state variables                                          |
| 2054-2055     | READ      | `use irrigation_mod, only: irrigation`; `use oxygenstress_mod, only: oxygen_dat`                                    |
| 2056-2057     | READ      | `use array_utils, only: afgen`; `use swap_array_dimensions, only: magrs, mayrs`                                    |
| 2058          | READ      | `implicit none`                                                                                                      |
| 2060-2061     | READ      | Dummy argument declarations: `icrop,lcc,swhydrlift`, `crpfil`                                                       |
| 2064-2071     | READ      | Local variable declarations (`crp,getun2,ifnd,i,ini`, `swIrrigate`, `sum,depth,rootdis`, temp arrays, `message,filnam`, `rdinqr`) |
| 2075          | READ      | `filnam = trim(pathcrop)//trim(crpfil)//'.crp'` — build file path                                                   |
| 2076          | READ      | `crp = getun2(10,90,2)` — get free file unit number                                                                 |
| 2077          | READ      | `call rdinit(crp,logf,filnam)` — open .crp file and initialise parser                                               |
| 2080          | READ      | `call rdsinr('idev',1,2,idev)` — read phenology switch                                                              |
| 2081-2082     | READ      | `if (idev.eq.1)`: `call rdsinr('lcc',1,366,lcc)` — fixed crop cycle length                                         |
| 2083-2087     | READ      | `elseif (idev.eq.2)`: read `tsumea`, `tsumam`, `tbase` — temperature-sum phenology                                  |
| 2090-2091     | READ      | `call rdsdor('kdif',...)`, `call rdsdor('kdir',...)` — light extinction coefficients                                |
| 2094          | READ      | `call rdsinr('swgc',1,2,swgc)` — LAI vs soil cover fraction switch                                                  |
| 2095-2096     | READ      | `if (swgc.eq.1)`: `call rdador('gctb',0,12,gctb,...)` — LAI table                                                  |
| 2097-2104     | READ      | `elseif (swgc.eq.2)`: `call rdador('gctb',0,2,gctb,...)` — soil cover fraction table                               |
| 2099-2104     | VALIDATE  | Loop over gctb even elements: `if (gctb(i).gt.1.0)` fatalerr — SoilCover must be ≤ 1 (swgc=2)                      |
| 2108          | READ      | `call rdsinr('swcf',1,3,swcf)` — crop factor / crop height switch                                                   |
| 2111-2115     | VALIDATE  | `if (swetr.eq.1 .and. swcf.eq.2)` fatalerr — ETref requires swcf=1 or 3, not swcf=2                                |
| 2117-2125     | READ      | `if (swcf.eq.1 .or. swcf.eq.3)`: read dvs array and cf array into `cftb`                                           |
| 2122-2125     | RUNTIME   | Loop: pack `dvsinput`/`cfinput` into interleaved `cftb` array                                                       |
| 2126-2135     | GUARDED   | `if (swcf.eq.3)`: read additional dvs+cfw arrays into `cfeictb` (wet-crop factor). Port scope is swcf=1; swcf=3 is guarded. |
| 2136          | NORMALIZE | `chtb = -99.99d0` — sentinel fill for unused crop-height table when swcf=1 or 3                                     |
| 2137-2146     | READ      | `else` (swcf=2): read dvs+ch arrays into `chtb`                                                                     |
| 2142-2145     | RUNTIME   | Loop: pack `dvsinput`/`chinput` into interleaved `chtb`                                                             |
| 2146          | NORMALIZE | `cftb = -99.99d0` — sentinel fill for unused crop-factor table when swcf=2                                          |
| 2150-2154     | NORMALIZE | `if (swcf.eq.1 .or. swcf.eq.3)`: `albedo = 0.23d0`, `rsc = 70.0d0`, `rsw = 0.0d0` — hardcoded ETref standard defaults |
| 2155-2160     | READ      | `else` (swcf=2): `call rdsdor('albedo',...)`, `call rdsdor('rsc',...)`, `call rdsdor('rsw',...)` — crop-specific values |
| 2165          | NORMALIZE | `swrdc = 0` — default: root density distribution does not change in time                                             |
| 2168          | READ      | `call rdador('rdctb',0,100,rdctb,22,ifnd)` — root density distribution table (always read, 11 pairs)                |
| 2171          | NORMALIZE | `swrd = 1` — default: root extension depends on development stage                                                    |
| 2172-2174     | READ      | `if (rdinqr('swrd'))`: `call rdsinr('swrd',1,3,swrd)` — optional override of root extension switch                  |
| 2177-2179     | READ      | `if (swrd.eq.1)`: `call rdador('rdtb',0,1000,rdtb,...)` — rooting depth vs DVS table                               |
| 2182-2192     | GUARDED   | `elseif (swrd.eq.2)`: read `rdi`, `rri`, `rdc`, optional `swdmi2rd` — max-daily-increase root extension. Not in port scope. |
| 2195-2200     | GUARDED   | `elseif (swrd.eq.3)` fatalerr — root biomass extension not supported by simple cropgrowth; always fatal here        |
| 2204          | NORMALIZE | `dvsend = 2.0d0` — default end-of-season DVS                                                                        |
| 2205-2207     | READ      | `if (rdinqr('dvsend'))`: `call rdsdor('dvsend',0,3,dvsend)` — optional override                                     |
| 2208          | NORMALIZE | `swharv = 0` — default: no harvest switch                                                                            |
| 2209-2211     | READ      | `if (rdinqr('swharv'))`: `call rdsinr('swharv',0,1,swharv)` — optional harvest timing switch                        |
| 2209-2211     | GUARDED   | swharv=1 activates harvest scheduling (not in port scope; swharv=0 is the default path)                              |
| 2214          | NORMALIZE | `swoxygen = 1` — default oxygen stress model (Feddes)                                                               |
| 2215-2216     | READ      | `if (rdinqr('swoxygen'))`: `call rdsinr('swoxygen',0,2,swoxygen)` — optional override                               |
| 2218-2223     | VALIDATE  | `if (swoxygen.eq.2 .and. (swhea.eq.0 .or. swcalt.eq.1))` fatalerr — physical O2 stress needs numerical heat flow    |
| 2225-2230     | VALIDATE  | `if (swoxygen.eq.2 .and. bdens(1).lt.100)` fatalerr — physical O2 stress needs realistic bulk density               |
| 2233-2238     | READ      | `if (swoxygen.eq.1)`: read `hlim1`, `hlim2u`, `hlim2l` — Feddes anaerobiosis pressure heads                         |
| 2240-2278     | GUARDED   | `if (swoxygen.eq.2)` block — Bartholomeus physical O2 stress. Entire block is out of port scope.                    |
| 2243          | GUARDED   | `swoxygentype = 1` — default physical O2 sub-model (swoxygen=2 only)                                                |
| 2244-2246     | GUARDED   | Optional read of `swoxygentype` (swoxygen=2 only)                                                                    |
| 2248-2270     | GUARDED   | `if (swoxygentype.eq.1)`: read q10_root, q10_microbial, specific_resp_humus, c_mroot, srl, f_senes, wrtb, mrftb, swrootradius, and swrootradius-dependent radius params |
| 2271-2277     | GUARDED   | `else` (swoxygentype=2): read SwTopSub, NrStaring; call `oxygen_dat(...)` — reproduction function approach           |
| 2282          | NORMALIZE | `swWrtNonox = 0` — default: no non-oxidative root senescence switch                                                  |
| 2283-2285     | READ      | `if (rdinqr('swWrtNonox'))`: `call rdsinr('swWrtNonox',0,1,swWrtNonox)` — optional override                         |
| 2287          | NORMALIZE | `aeratecrit = 0.0001d0` — default aeration criterion                                                                 |
| 2288-2290     | READ      | `if (swWrtNonox.eq.1)`: `call rdsdor('aeratecrit',...)` — read aeration criterion threshold                          |
| 2293          | NORMALIZE | `swdrought = 1` — default: drought stress according to Feddes                                                        |
| 2294-2296     | READ      | `if (rdinqr('swdrought'))`: `call rdsinr('swdrought',1,2,swdrought)` — optional override                             |
| 2298-2305     | READ      | `if (swdrought.eq.1)`: read `hlim3h`, `hlim3l`, `hlim4`, `adcrh`, `adcrl` — Feddes drought parameters               |
| 2307-2320     | GUARDED   | `else` (swdrought=2): read De Jong van Lier drought params: `wiltpoint`, `kstem`, `rxylem`, `rootradius`, `kroot`, `rootcoefa`, `swhydrlift`, `rooteff`, `stephr`, `criterhr`, `taccur` |
| 2324-2347     | GUARDED   | `if (flsolute)` block — salt stress section. swsalinity∈{1,2} both guarded for port.                                |
| 2326-2328     | GUARDED   | `if (rdinqr('swsalinity'))`: read `swsalinity` (flsolute=true only)                                                  |
| 2330-2334     | GUARDED   | `if (swsalinity.eq.1)`: read `saltmax`, `saltslope` (Maas-Hoffman function)                                          |
| 2335-2345     | GUARDED   | `elseif (swsalinity.eq.2)`: read `salthead`; validate swdrought must be 2 with swsalinity=2                          |
| 2339-2345     | VALIDATE  | `if (swdrought.eq.1)` inside swsalinity=2 block: fatalerr — osmotic head salt stress requires De Jong van Lier       |
| 2350-2351     | NORMALIZE | `swcompensate = 0; swjarvis = 0` — defaults for compensation switches                                                |
| 2352-2353     | READ      | `if (rdinqr('swcompensate'))`: `call rdsinr('swcompensate',0,2,swcompensate)`                                        |
| 2354-2367     | READ      | `else if (rdinqr('swjarvis'))` legacy path: warn, read `swjarvis`, derive `swcompensate=1`                           |
| 2356-2357     | VALIDATE  | `call warn(...)` — SWJARVIS is deprecated; migrate to SWCOMPENSATE                                                   |
| 2361-2362     | VALIDATE  | `call warn(...)` — SWJARVIS applied to all stresses (when swjarvis≠4)                                                |
| 2364          | NORMALIZE | `swcompensate = 1` — coerce legacy swjarvis > 0 to new swcompensate=1                                                |
| 2370-2375     | VALIDATE  | `if (swcompensate.gt.0 .and. swdrought.eq.2)` fatalerr — compensation not allowed with De Jong van Lier drought      |
| 2378          | NORMALIZE | `swstressor = 1` — default: compensate all stressors                                                                 |
| 2379-2383     | READ      | `if (swcompensate.gt.0 .and. rdinqr('swstressor'))`: read `swstressor`                                               |
| 2386-2390     | READ      | `if (swcompensate.eq.1)`: `alphacrit = 1.0d0` default; read `alphacrit`                                              |
| 2386-2390     | NORMALIZE | `alphacrit = 1.0d0` — default critical stress index for Jarvis compensation                                          |
| 2393-2398     | GUARDED   | `elseif (swcompensate.eq.2)`: `dcritrtz = 0.0d0` default; read `dcritrtz` — Walsum compensation. Not in port scope. |
| 2402          | READ      | `call rdsinr('swinter',0,3,swinter)` — interception model switch                                                     |
| 2403-2404     | READ      | `if (swinter.eq.1)`: read `cofab` — simple von Hoyningen-Huene coefficient                                           |
| 2405-2423     | GUARDED   | `else if (swinter.eq.2)`: read t, pfree, pstem, scanopy, avprec, avevap; pack into 5 lookup tables. Not in port scope. |
| 2424-2427     | GUARDED   | `else if (swinter.eq.3)`: read `fimin`, `siccaplai`. Not in port scope.                                              |
| 2430          | READ      | `call rdsinr('schedule',0,1,schedule)` — irrigation scheduling switch                                                |
| 2432-2437     | GUARDED   | `if (schedule.eq.1 .and. swdrought.eq.2)`: re-read hlim3h, hlim3l, hlim4 for scheduling override. Scheduling is guarded (port scope is schedule=0). |
| 2440          | READ      | `close(crp)` — close the .crp file                                                                                   |
| 2443-2446     | GUARDED   | `if (schedule.eq.1)`: `call irrigation(1)` + optional `IrrigationOutput(1)` — scheduling initialisation. Guarded.   |
| 2449-2478     | RUNTIME   | `if (swdrought.eq.1)` block: build `cumdens` array — normalized cumulative root density function over 101 depth points |
| 2454-2458     | RUNTIME   | Inner loop i=0,100: populate `rootdis` using `afgen(rdctb,22,depth)` at 0.01-spaced depths                           |
| 2461-2464     | RUNTIME   | Copy depth values (odd indices) from `rootdis` into `cumdens`                                                         |
| 2465-2472     | RUNTIME   | Trapezoidal integration loop i=4,202,2: accumulate `sum` and store in `cumdens` even indices                          |
| 2475-2477     | RUNTIME   | Normalization loop i=2,202,2: `cumdens(i) = cumdens(i) / sum`                                                         |
| 2481-2514     | GUARDED   | `if (t1900 - tstart .lt. 1d-3 .and. swinco.eq.3 .and. ...)` — warm-restart from `.END` file (swinco=3). Entire block guarded. |
| 2485-2486     | GUARDED   | `ini = getun2(...)`; `call rdinit(ini,logf,inifil)` — open `.END` file (swinco=3 only)                               |
| 2489-2497     | GUARDED   | Read `swIrrigate` and conditionally `dayfix` from `.END` file (swinco=3 only)                                        |
| 2490-2494     | GUARDED   | Derive `flIrrigate` boolean from `swIrrigate` integer (swinco=3 only)                                                |
| 2500-2503     | GUARDED   | Read `rd`, `rdpot`, optionally `sicact` (swinco=3 only)                                                              |
| 2506-2509     | GUARDED   | Read `nofd`, `daycrop`, `dvs`, `tsum` from `.END` file (swinco=3 only)                                               |
| 2512          | GUARDED   | `close(ini)` — close `.END` file (swinco=3 only)                                                                     |
| 2516          | READ      | `return`                                                                                                              |
| 2517          | READ      | `end` — end of subroutine                                                                                             |

---

## Summary

### Lines per bucket

| Bucket    | Approx. line count | Description                                                                                   |
| --------- | ------------------ | --------------------------------------------------------------------------------------------- |
| READ      | ~120 lines         | `rd*` calls, file open/close, declarations, use statements                                    |
| VALIDATE  | ~40 lines          | `fatalerr` / `warn` calls and surrounding `if` logic                                          |
| NORMALIZE | ~20 lines          | Hardcoded defaults (`albedo=0.23`, `rsc=70`, `rdc*=-99.99`), default switches, sentinel fills |
| RUNTIME   | ~35 lines          | cumdens build (post-close, no file I/O); table packing loops for cftb/chtb                    |
| GUARDED   | ~100 lines         | swcf=3 (wet crop factor), swrd=2/3, swoxygen=2, swsalinity∈{1,2}, swcompensate=2, swinter∈{2,3}, schedule=1, swinco=3 warm restart |

Total: ~480 lines of substantive content (2037–2517).

### Call count summary

- **READ:** ~45 `rd*` calls (rdsinr, rdsdor, rdador, rdfdor, rdinqr, rdinit) plus 3 file operations (getun2, rdinit×2, close×2)
- **VALIDATE:** 7 `fatalerr` calls + 3 `warn` calls
- **NORMALIZE:** 10 distinct default assignments (swrdc, swrd, dvsend, swharv, swoxygen, swWrtNonox, aeratecrit, swdrought, swcompensate/swjarvis, swstressor, alphacrit; plus chtb/cftb sentinel fills, albedo/rsc/rsw ETref defaults)
- **RUNTIME:** 1 initialization action: cumdens build (3 loops); 2 table-packing loops (cftb, chtb)
- **GUARDED:** 8 conditional branches: swcf=3, swrd=2, swrd=3, swoxygen=2, swsalinity∈{1,2}, swcompensate=2, swinter∈{2,3}, schedule=1, swinco=3

---

## Implementation Guidance

### 1. Many switches have hardcoded defaults via `rdinqr` guards

Unlike `rddre`, almost every optional parameter in `readcropfixed` uses the pattern:

```fortran
swXxx = <default>
if (rdinqr('swXxx')) then
  call rdsinr ('swXxx', ...)
endif
```

This means the TOML path must preserve these defaults in the config finalizer or schema. The defaults are: `swrdc=0`, `swrd=1`, `dvsend=2.0`, `swharv=0`, `swoxygen=1`, `swWrtNonox=0`, `aeratecrit=1e-4`, `swdrought=1`, `swcompensate=0`, `swstressor=1`, `alphacrit=1.0`. All must be documented in the TOML schema and applied in the finalizer.

### 2. albedo/rsc/rsw are hardcoded for swcf=1 or swcf=3

Lines 2152-2154 silently override any user-supplied values with `albedo=0.23`, `rsc=70.0`, `rsw=0.0` whenever `swcf=1` or `swcf=3`. This is a NORMALIZE action, not a default — it happens *after* any read, so user values are discarded. The TOML config validator should either refuse these fields when swcf≠2, or the finalizer must enforce the override and warn if user-supplied values differ.

### 3. cftb/chtb sentinel fill is bidirectional

When `swcf=1 or 3`, `chtb` is filled with `-99.99` (line 2136). When `swcf=2`, `cftb` is filled with `-99.99` (line 2146). These are runtime invariants that downstream code uses as "not-applicable" markers. The `cropfixed_init` module (Task 6) must reproduce this fill; it is not a config-time concern.

### 4. cumdens is built after file close (lines 2449-2478)

The normalized cumulative root density function is computed *after* `close(crp)` at line 2440. It depends only on `rdctb` and `swdrought`. This is a pure post-read computation and belongs in `cropfixed_init`. It is only computed when `swdrought=1` — for `swdrought=2` the downstream code uses a different mechanism (De Jong van Lier does not use cumdens).

### 5. swoxygen=2 block is large but entirely guarded (~38 lines, 2240-2278)

The Bartholomeus physical oxygen stress model has two sub-variants (swoxygentype=1: physical; swoxygentype=2: reproduction functions via `oxygen_dat`). Both are out of port scope. The stub should refuse swoxygen=2 at config validation time (the validator already exists for the `fatalerr` cases at 2218-2230, but those only guard swcf/swhea/bdens combinations; swoxygen=2 itself needs an explicit "not supported" error in the new validator).

### 6. swinco=3 warm-restart block reads a second file (the `.END` file)

Lines 2481-2514 conditionally open `inifil` (the `.END` restart file from a previous run) and read runtime state: `swIrrigate`, `dayfix`, `rd`, `rdpot`, `sicact`, `nofd`, `daycrop`, `dvs`, `tsum`. This is **not** crop parameter config — it is runtime state persistence. The TOML port does not need to replicate this; the `.END`/`swinco=3` path should be refused by the config validator (`swinco=3` is already documented as guarded in the swap-config module).

### 7. The swjarvis legacy path (lines 2354-2367) is a backward-compat shim

It reads the deprecated `SWJARVIS` key, warns, and coerces to `swcompensate=1`. In TOML the shim should not be replicated — the TOML schema should only expose `swcompensate`. If a user has an old `.crp` with `swjarvis`, they must migrate manually; the validator should reject unknown keys.

### 8. swsalinity is gated on the runtime flag `flsolute`

The entire salt stress section (lines 2324-2347) is inside `if (flsolute)`. The `flsolute` flag comes from the `.swp` file (not the `.crp` file), meaning salt stress configuration is conditional on a cross-file runtime dependency. The TOML crop validator cannot check this alone — either the cross-file validator must gate swsalinity checks on the companion `.swp` config, or the finalizer must apply the gate.

### 9. schedule=1 re-reads hlim3h/hlim3l/hlim4 (lines 2432-2437)

When both `schedule=1` and `swdrought=2`, the three Feddes limiting pressure heads are re-read from a dedicated irrigation block in the `.crp` file. These are the **same** variable names as the Feddes drought block (lines 2301-2303), but appear in a different position in the file, possibly with different values. The TOML schema should treat these as separate named fields (e.g. `irrigation_hlim3h` vs `hlim3h`) or use a sub-table, to avoid ambiguity. Since schedule=1 is guarded in the port, this can be deferred.

### 10. swrd=3 is always a fatal error (line 2195-2199)

`readcropfixed` immediately calls `fatalerr` for swrd=3. This means swrd=3 was never validly supported by the simple cropgrowth module. The TOML validator should reject swrd=3 as an invalid value (range 1–2) rather than as an "unsupported" port path.

### 11. interception table packing (swinter=2, lines 2412-2423) stores 5 parallel tables

When `swinter=2`, 6 arrays (t, pfree, pstem, scanopy, avprec, avevap) are read and then repacked into 5 interleaved lookup tables (`pfreetb`, `pstemtb`, `scanopytb`, `avprectb`, `avevaptb`) using the `tinter` array as the common x-axis. The `cropfixed_init` module, if it ever supports swinter=2, must reproduce this packing exactly.
