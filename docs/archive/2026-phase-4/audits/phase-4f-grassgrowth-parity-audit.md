# Phase 4f — grassgrowth exact-parity audit

Phase 4f's strangler swap (commit 744cbfd / abf76bd) routes case 2 (grassgrowth)
through the new TOML pipeline. With the existing hupselbrook-tuned schema +
adapter, grassgrowth aborts hard at startup with

    [300] FATAL heat: psand/psilt/pclay/porg must all be set together

— the TOML for case 2 still authors only `psand/pclay/porg`. Beyond this
hard-block there are several latent gaps inherited from B2: the soil
discretization and MvG hydraulics blocks aren't authored at all, the .dra is
DRAMET=3 (resistance, multi-level — different code path from hupselbrook's
DRAMET=2 Hooghoudt), and the bottom boundary is `swbotb=1` with an external
.bbc file (gwltab branch — unwired in adapter).

This audit walks every legacy `rd*` call in `src/io/readswap.f90` against the
case 2 `.swp` template (`tests/swap-cases/2.grassgrowth/swap_linux.swp.template`)
and `.dra` (`tests/swap-cases/2.grassgrowth/swap.dra`).

Status legend (mirrors hupselbrook audit):
- OK — TOML authors it, reader copies it, adapter writes the legacy global,
  no finalize mismatch.
- TOML — value missing from `tests/swap-cases/toml/2.grassgrowth/swap.toml`
  (or `swap.dra.toml`).
- READER — schema field exists but no reader populates it.
- ADAPTER — schema + reader fine, but adapter never copies the value into the
  legacy global.
- FINALIZE — adapter copies the value but skips a unit-conversion or derived
  finalize step the legacy reader performs.
- N/A — switch is read but its value is irrelevant (gated branch not taken).

## grassgrowth .swp — General + Output

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| project, paths  | OK       | adapter line 59-63 |
| swscre, swerror | OK       | swerror has no schema slot but Initialize zeroes; case authors 0 |
| tstart / tend   | OK       | 1980-01-01 → 1984-12-31 |
| nprintday       | OK       | 1 |
| swmonth         | OK       | =0 (daily) — adapter's `populate_outdatint_monthly` skipped |
| period, swres, swodat | OK | period=1, swres=0, swodat=0 |
| swyrvar, datefix | TOML/N/A| swyrvar=0; gated by SWBAL/SWBLC which are RETIRED. N/A. |
| outfil          | OK (HACK)| hard-coded `'result'` |
| swheader        | OK (HACK)| zero-forced |
| swwba/swend/swvap/swbal/swblc/swsba/swate/swbma/swdrf/swswb/swini/swinc/swcrp/swstr/swirg | OK (HACK) | RETIRED switches zero-forced (ADR 0009) |
| **swcsv, inlist_csv** | **TOML/ADAPTER** | adapter hardcodes hupselbrook's water-balance inlist (`rain,...,gwl`); case 2 fixture asserts on **GRASSDM/MOWDM/PGRASSDM/PMOWDM** (grass-detailed output). With wrong inlist the CSV will not contain those columns and the regression aggregator will receive empty rows. **HARD BLOCK on regression compare.** |
| swcsv_tz, inlist_csv_tz | OK (HACK) | hard-coded zero |
| swafo, swaun, critdevmasbal, swdiscrvert, numnodnew, dznew | partial | swafo/swaun zero-forced; critdevmasbal/discretization not. swdiscrvert default=0 → no extra discretization. OK. |

## grassgrowth .swp — Meteorology

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| metfil          | OK       | `260.met` — `.met`-suffix → `swMetFilAll=1` after adapter |
| lat, alt, altw  | OK       | 52.1 / 1.9 / 10.0 |
| swetr           | OK       | =0 |
| angstroma/b     | OK       | 0.25 / 0.5 |
| swdivide        | OK       | =1 → PM-direct |
| swmetdetail, nmetdetail | OK | =0 |
| swetsine        | OK       | =0 |
| **swrain**      | **TOML** | grassgrowth authors `swrain=2` (daily + duration). TOML already has it. OK. |
| rainfil, rainflux | OK     | unused at swrain=2 |

## grassgrowth .swp — Crop section

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| swcrop          | OK       | =1; adapter arms flCropReadFile/Open |
| cropstart, cropend, cropfil, croptype | OK | adapter writes per-rotation arrays. Five `grassd` rotations 1980..1984 type=3 (WOFOST grass) |
| rds (rdmax)     | OK (HACK)| hard-coded 200 cm in adapter; case authors 200. Match. |

## grassgrowth .swp — Irrigation

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| swirfix         | OK       | =0 (no irrigation), no `fixed_events` table needed |

## grassgrowth .swp — Soil water

| Legacy key      | Status        | Notes |
|-----------------|---------------|-------|
| swinco          | OK            | =2 |
| gwli            | OK            | -75.0 cm |
| swpondmx        | TOML/ADAPTER  | =0; adapter never sets it; default 0; OK. |
| pondmx          | OK            | 0.2 cm |
| **rsro**        | **TOML**      | legacy reads 0.5 d. Schema slot exists. Case 2 TOML doesn't author it. **DRAINAGE/runoff drift.** |
| **rsroexp**     | **TOML**      | legacy reads 1.0. Same as rsro. |
| swrunon         | OK            | =0; adapter sets `flrunon=.false.`; no read needed |
| cfevappond      | OK (HACK)     | 1.25 |
| swcfbs          | OK            | =0; cfbs unused |
| **rsoil**       | **TOML**      | legacy reads **600.0 s/m** (NOT 30 like case 1). Schema default 0 → PM-direct rss=0 → EPOT direct contributor. **Highest single-key impact on EPOT.** |
| swredu          | OK (HACK)     | =1 (Black) |
| rsigni, cofred (cofredbl) | OK | rsigni HACKed to 0.5; cofredbl=0.35 in case 2 (adapter maps cofredbl → cofred) |
| **isoillay/isublay/hsublay/ncomp** | **TOML missing** | grassgrowth has **9 sub-layers / 5 soil-physical layers** vs hupselbrook's 4/2. The current case 2 swap.toml has NO discretization keys at all — `nsublay`, `isoillay`, `hsublay`, `ncomp` stay 0/empty after adapter, calcgrid will fail or produce a degenerate one-compartment profile. **HARD BLOCK on water flow.** |
| swsophy         | OK            | =0 (analytical MvG) |
| **MvG params**  | **TOML missing** | 5-row MvG table in .swp; case 2 swap.toml authors no `[soil.hydraulics]`. Schema slot exists. **HARD BLOCK on Richards solver.** |
| swhyst, tau     | OK            | swhyst=0 |
| swmacro         | OK            | =0 |
| swsnow          | OK            | =0 |
| swfrost         | OK            | =0 |
| Numerical (dtmin, dtmax, gwlconv, critdevh1cp, critdevh2cp, critdevponddt, maxit, maxbacktr, swkmean, swkimpl) | TOML | hupselbrook-style block sets dt/dtmin/dtmax only; gwlconv etc. at schema defaults. swkmean default in legacy=1, schema currently 0 — must check. |

## grassgrowth .swp — Drainage / .dra (DRAMET = 3, multi-level resistance)

| Legacy key      | Status        | Notes |
|-----------------|---------------|-------|
| swdra, drfil    | OK            | swdra=1, drfil='swap' |
| **dramet**      | **OK**        | =3 (resistance multi-level). Different code path from B2's DRAMET=2 — no Hooghoudt globals (lm, wetper, zbotdr_basic, ipos, khtop, etc.) involved. |
| swdivd          | OK            | =1 |
| cofani[1..maho] | OK            | 5 elements (one per soil-physical layer); TOML already authors `cofani = 1.0 1.0 1.0 1.0 1.0` — wait, case 2 swap.dra.toml currently authors no cofani. **TOML missing.** |
| swdislay        | OK            | =0 |
| nrlevs          | OK            | =1 |
| **swintfl**     | **ADAPTER**   | =0 in case 2. Legacy at readswap.f90:1889 reads it then assigns `swnrsrf = swintfl`. No-op when 0. Adapter already sets swnrsrf from config. OK. |
| cofintflb / expintflb | N/A     | gated by swnrsrf=1; swintfl=0 here so skipped. |
| swliminf        | OK (legacy default) | legacy assigns 1 unconditionally (readswap.f90:2010); the adapter doesn't touch it. variables.f90:694 default-initialises to 0. **MISSING** — this gates limit-of-infiltration logic in the drainage solver. |
| **drares[1]**   | **TOML/ADAPTER** | =750.0 d. TOML authors it under `[[drainage.levels]]`. Reader copies into per-level arrays. Adapter copies. OK. |
| **infres[1]**   | **OK**        | =2000.0 d. Same pipeline. |
| **swallo[1]**   | **TOML**      | legacy reads (1..3); =1 in case 2. TOML lacks the `swallo` field on the level; reader has no slot. **TOML + READER missing.** Default 0 → drainage code may treat as "drainage not allowed" or behave surprisingly. |
| **L[1]**        | **TOML**      | =500.0 m. TOML authors `L = 500.0`. Reader copies. Adapter copies as-is. **FINALIZE BUG** — legacy at readswap.f90:1908-1909 multiplies by 100 (m→cm); adapter only does m→cm conversion in DRAMET=2 branch (line 237). For DRAMET=3 the multiplication is missing — `L(1)` stays as 500 cm instead of 50 000 cm. |
| **zbotdr[1]**   | **OK**        | -55.0 cm; reader/adapter copy as-is (no unit conversion in legacy DRAMET=3 either). |
| **swdtyp[1]**   | **OK**        | =2 (open channel); TOML authors. Reader/adapter copy. |
| **datowl1/level1** | **TOML/READER/ADAPTER** | The channel-water-level table goes into `owltab(1, :)` and `nowltab(1)`. The legacy reader reads (date, level) pairs; the schema currently has no slot, the reader can't see them, and the adapter never sets `owltab/nowltab`. drainage.f90:263 calls `afgen(owltab(lev, 1:2*nowltab(lev)), ...)` — **HARD BLOCK on drainage at runtime.** |
| swliminf        | TOML/ADAPTER  | (see above) |
| ipos / khtop / khbot / kvtop / kvbot / zintf / geofac | N/A | DRAMET=2-only. Skip. |

## grassgrowth .swp — Bottom boundary (SWBBCFILE=1, swap.bbc)

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| swbbcfile, bbcfil | OK     | =1 + 'swap'; case TOML authors `bbcfil = "swap"`. (No actual switch in TOML — the legacy reader gates on swbbcfile, but the new pipeline only uses bbcfil to find the file.) |
| swbotb          | OK       | =1 (prescribed groundwater level) |
| **date1/gwlevel** (gwltab) | **TOML/READER/ADAPTER** | The 120-row gwlevel time series in swap.bbc. Legacy stores it in `gwltab(1:2*ifnd)`. The new pipeline neither reads the .bbc file nor writes `gwltab`. The schema has `swc_table` (legacy alias used inline) but no per-case data; for swbotb=1 the bottom-boundary code calls `afgen(gwltab, 2*nbbc, t1900+dt-1.d0)` — **HARD BLOCK on bottom flux at runtime.** |

## grassgrowth .swp — Heat

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| swhea           | OK       | =1 |
| swcalt          | OK       | =2 (numerical) |
| **psand/psilt/pclay/porg** | **TOML** | 5-row table in .swp; case 2 TOML authors psand/pclay/porg only — **psilt array missing**. Validator aborts startup with [300] FATAL. **HARD BLOCK at startup**, currently the only failure mode visible. |
| psand, pclay, orgmat (porg) | OK | adapter copies; psand/pclay/porg already authored. |
| tsoil_init      | OK       | adapter copies col 2 element-by-element |
| swtopbhea, swbotbhea | OK   | both =1 |

## grassgrowth .swp — Solute

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| swsolu          | OK       | =0 — section short-circuits |

(All other solute fields skipped when swsolu=0.)

## Highest-priority gaps (water-balance-impact-ordered)

These are all gating regression to even produce a comparable run, ranked by
crash-priority then by water-balance impact.

### Tier 0 — startup-aborting gaps (case won't reach run loop)

1. **psilt (heat)** — Validator [300] FATAL aborts startup. Must add
   `psilt = [...]` to `[heat]` in case 2 swap.toml. Schema slot already
   exists from hupselbrook B2.
2. **Soil discretization missing** (`isoillay`, `hsublay`, `ncomp`) — calcgrid
   needs ≥1 sub-layer; hupselbrook authored 4 lines, case 2 needs 9 lines.
   Schema fields exist. Just author them in TOML.
3. **MvG hydraulics missing** (`[soil.hydraulics].ores/osat/...`) — Richards
   solver needs paramvg(1..10, 1..numlay) populated. Schema fields exist; just
   author the 5-row table.

### Tier 1 — runtime-aborting (case crashes inside dynamic loop)

4. **bbcfil + gwltab plumbing** (swbotb=1) — The legacy reader for swap.bbc
   isn't wired into the new pipeline. Two options:
   (a) extend bottom_boundary schema with a `gwl_table[(date,gwl)]` slot, port
       the .bbc reader, populate `gwltab` in adapter; or
   (b) author the gwl_table inline in case 2 swap.toml (simpler — skips the
       .bbc file altogether). Legacy variables.f90 reads gwltab from a 2*mabbc
       flat array indexed [date, gwl, date, gwl, ...].
   Either way: schema needs `bottom_boundary.gwl_table` (a list of [date, gwl]
   pairs) and the adapter needs a finalize that fills `gwltab(2*i-1) = days,
   gwltab(2*i) = level` and a count `nbbc` if there's a global counter.
5. **datowl1/level1 + owltab(1,:) plumbing** (DRAMET=3) — mirrors gap 4 for
   the drainage side. Schema needs `[[drainage.levels]].channel_water_levels`
   (a list of [date, level] pairs); adapter fills `owltab(lev, 1..2*ifnd)` and
   `nowltab(lev)`.
6. **L[1] m→cm conversion** for DRAMET=3 — adapter only converts at DRAMET=2
   branch. For DRAMET=3 cases the per-level `L(:)` array stays in metres,
   under-by-100x. Either fix in the adapter (per-level loop converts when
   `dramet == 3`) or in the reader's pre-finalize. Adapter is simpler.
7. **swallo per-level** — schema/reader missing the `swallo` slot on the
   `[[drainage.levels]]` table. Default-zero may fall into an unallowed
   combination in drainage.f90.

### Tier 2 — water-balance correctness (case completes but diff > 1e-2 cm)

8. **rsoil = 600.0 s/m** — case 2 writes 600 (NOT 30). PM-direct evaporation
   driver. Largest expected-EACT impact among the soft gaps once the case
   actually runs.
9. **rsro = 0.5, rsroexp = 1.0** — same impact category as case 1; surface
   runoff branch. Drainage / DSTOR drift.
10. **cofani** — case 2 has 5 layers; TOML authors none. Default 0 collapses
    horizontal K → 0 → drainage flux artificially zero. **Probably matters
    a lot once the case runs.**
11. **swliminf default** — legacy hard-codes 1 in DRAMET=3 reader; new pipeline
    leaves at 0. Adapter HACK suffices: `swliminf = 1`.

### Tier 3 — Output-aggregator-blocking (case completes, fixture compare fails)

12. **inlist_csv** — adapter hard-codes hupselbrook's water-balance vars; case
    2 fixture asserts on GRASSDM/MOWDM/PGRASSDM/PMOWDM. Case-specific override
    needed. Either schema slot or per-case branch in adapter (HACK).

## Action plan (one commit per gap, smallest commit first, submodule-first)

1. **(submodule)** Author `psilt` + soil discretization + MvG hydraulics in
   `tests/swap-cases/toml/2.grassgrowth/swap.toml`. Outer-repo bump.
2. **(submodule)** Author `cofani` + `rsoil/rsro/rsroexp` + `swallo` in case 2
   TOMLs. Outer-repo bump.
3. **(adapter)** Fix DRAMET=3 L[:] m→cm conversion. No schema change needed.
4. **(adapter)** Add `swliminf = 1` HACK marker, mirroring legacy default.
5. **(schema + reader + adapter)** Add `[[drainage.levels]].swallo` (numeric
   1..3) — plus the `[[drainage.levels]].channel_water_levels` array of
   `[date, level]` pairs. Adapter fills `owltab(lev, :)` / `nowltab(lev)`.
6. **(schema + reader + adapter)** Add `[bottom_boundary].gwl_table` array of
   `[date, gwl]` pairs (or load `bbcfil`-named external file). Adapter fills
   `gwltab` and the bbc count.
7. **(adapter)** Per-case `inlist_csv` override (HACK) — branch on
   `cropfil(:)` containing 'grassd' to write the grass var list. Or: add
   `[output.csv].inlist` to schema and consume.

After each commit: build and run `pixi run -e test regression grassgrowth` to
record delta. After the whole sequence, re-run hupselbrook regression to
confirm no regressions in case 1.
