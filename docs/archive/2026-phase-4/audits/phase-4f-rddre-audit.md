# rddre Line-by-Line Audit

**File:** `src/io/readswap.f90`, lines 4273–4912  
**Subroutine:** `rddre(wls1, wlp1)`  
**Date of audit:** 2026-05-01  
**Purpose:** Classification of every non-blank, non-comment line in the 642-line `rddre` subroutine to inform the swap.dra → TOML port (Tasks 4–8). The subroutine reads the legacy `swap.dra` ASCII file, validates its contents, applies coordinate/unit conversions, and builds runtime data structures. After the port, it is replaced by:

- Config validators (VALIDATE-bucket lines)
- Config finalizers (NORMALIZE-bucket lines)
- A new `surfacewater_init` module (RUNTIME-bucket lines)
- Stub/error checks for unimplemented branches (GUARDED-bucket lines)

Lines that fall into multiple buckets appear in multiple rows — this is intentional and useful.

---

## Classification Table

| Lines         | Bucket    | Notes                                                                                     |
| ------------- | --------- | ----------------------------------------------------------------------------------------- |
| 4273          | READ      | Subroutine signature `rddre(wls1,wlp1)`                                                   |
| 4295-4303     | READ      | `use` declarations (variables, surfacewater_utils, array_utils, swap_array_dimensions)    |
| 4304          | READ      | `implicit none`                                                                            |
| 4307          | READ      | `real(8) wls1, wlp1` — dummy args; become outputs of `surfacewater_init`                 |
| 4310-4324     | READ      | Local variable declarations (dra, level, itab, imper*, flweir, etc.)                     |
| 4329          | READ      | Build `filnam` = pathdrain // drfil // '.dra'                                             |
| 4330          | READ      | `dra = getun2(10,90,2)` — get free unit number                                            |
| 4331          | READ      | `call rdinit(dra,logf,filnam)` — open swap.dra and init parser                           |
| 4334          | READ      | `call rdsinr('swdivd',0,1,swdivd)` — read swdivd                                        |
| 4335-4339     | VALIDATE  | `if (swdivd.eq.0)` warn about numerical instability (call warn, not fatalerr)            |
| 4341-4353     | READ      | `if (swdivd.eq.1)`: conditionally read swdivdinf, FacDpthInf, cofani (scalar or array)  |
| 4356          | READ      | `call rdsdor('altcu',...)` — read altcu (altitude of control unit)                       |
| 4361          | READ      | `call rdsinr('nrsrf',1,Madr,nrlevs)` — read number of drainage levels                   |
| 4362-4367     | VALIDATE  | `if (nrlevs.gt.5)` warn about limited output (call warn)                                 |
| 4370-4379     | READ      | `if (swdivd.eq.1)`: read swdislay, conditionally read swtopdislay/ztopdislay or ftopdislay |
| 4383-4394     | READ      | Read per-level arrays: lev, swdtyp, l, zbotdre, gwlinf, rdrain, rinfi, rentry, rexit, widthr, taludr |
| 4397-4434     | NORMALIZE | Loop i=1,nrlevs: `l(i) = l(i)*100.0` (m→cm), `zbotdr(i) -= altcu` (absolute→relative)  |
| 4401-4407     | NORMALIZE | `if (swdivd=1 .and. swdislay=1)`: `ztopdislay(i) -= altcu`                               |
| 4403-4406     | VALIDATE  | `if (ztopdislay(i).lt.zbotdr(i))` fatalerr — ztopdislay cannot be below zbotdre          |
| 4414-4417     | VALIDATE  | `if (level(i).ne.i)` fatalerr — drainage level index not consistent                     |
| 4418-4423     | VALIDATE  | `if (i.gt.1 .and. zbotdr(i).lt.zbotdr(i-1))` fatalerr — levels not ordered (deepest first) |
| 4424-4427     | VALIDATE  | `if (gwlinf(i).gt.zbotdr(i))` fatalerr — gwlinf must be below zbotdr                    |
| 4428-4432     | VALIDATE  | `if (swdtyp(i).eq.0 .and. widthr(i).lt.small)` fatalerr — widthr > 0 for open channels  |
| 4437          | READ      | `call rdsinr('swnrsrf',0,2,swnrsrf)` — read swnrsrf                                     |
| 4438-4444     | READ      | `if (swnrsrf.eq.1)`: read rsurfdeep, rsurfshallow; `if (swnrsrf.eq.2)`: read cofintfl, expintfl |
| 4448-4450     | READ      | `if (swdivd.eq.1 .and. swnrsrf.gt.0)`: read SwTopnrsrf                                  |
| 4454          | READ      | `call rdsinr('swsrf',1,3,swsrf)` — read swsrf (surface water regime switch)             |
| 4455-4460     | RUNTIME   | Set nrpri/nrsec based on swsrf (`if swsrf.eq.3: nrpri=1; if swsrf.eq.2: nrpri=0; nrsec=nrlevs-nrpri`) — module-level state, not a `rd*` call |
| 4461-4465     | GUARDED   | `if (swsrf.eq.1)`: no surface water system — close and return. Our port targets swsrf=2; swsrf=1 closes early. |
| 4467-4470     | VALIDATE  | `if (swdtyp(1+nrpri).ne.0)` fatalerr — deepest secondary level must be open channel     |
| 4475-4500     | GUARDED   | Entire `if (swsrf.eq.3)` block — primary system wlptab read + wlp1 init. Port scope is swsrf=2 only. |
| 4477-4479     | GUARDED   | Init wlptab to zero (swsrf=3 only)                                                       |
| 4481-4484     | GUARDED   | `call rdatim('date1',...)` + `checkdate(...)` — read primary system date series (swsrf=3) |
| 4485-4490     | GUARDED   | `call rdfdor('wlp',...)` + loop: convert wlp to wlptab (altcu-relative), store dates     |
| 4492          | GUARDED   | `wlp1 = afgen(wlptab,2*mawlp,t1900-1.d0)` — set initial primary wl via interpolation   |
| 4504          | READ      | `call rdsinr('swsec',1,2,swsec)` — read swsec (secondary system management switch)      |
| 4508-4529     | GUARDED   | Entire `if (swsec.eq.1)` block — secondary wl is prescribed (wlstab read + wls1 init via afgen). Port scope is swsec=2. |
| 4513-4515     | GUARDED   | Init wlstab to zero (swsec=1 only)                                                       |
| 4517-4520     | GUARDED   | `call rdatim('date2',...)` + `checkdate(...)` — read secondary system date series (swsec=1) |
| 4522-4527     | GUARDED   | `call rdfdor('wls',...)` + loop: convert wls to wlstab (altcu-relative), store dates     |
| 4529          | GUARDED   | `wls1 = afgen(wlstab,2*mawls,t1900-1.d0)` — set initial secondary wl via interpolation (swsec=1) |
| 4534-4537     | READ      | `elseif (swsec.eq.2)`: read wlact, convert to altcu-relative: `wls1 = wls1 - altcu`     |
| 4536-4537     | NORMALIZE | `wls1 = wls1 - altcu` — make initial secondary wl relative to control unit              |
| 4538          | READ      | `call rdsdor('osswlm',0.0,10.0,osswlm)` — read allowed overshoot of target wl           |
| 4542          | READ      | `call rdsinr('nmper',1,mamp,nmper)` — read number of management periods                 |
| 4543          | RUNTIME   | `wlstar = wls1` — initialize target level to initial water level                        |
| 4545-4550     | READ      | Read per-period arrays: imper_4b, impend, swman, wscap, wldip, intwl                    |
| 4553-4554     | RUNTIME   | `nrman1=0; nrman2=0` — counter initializations derived from swman array, not `rd*` calls |
| 4555-4573     | VALIDATE  | Loop over nmper: validate intwl.ge.1 when swman=2; compute wldip=abs(wldip); count nrman1/nrman2; validate swman range; validate nrman1+nrman2=nmper |
| 4560          | NORMALIZE | `wldip(imper) = abs(wldip(imper))` — normalize wldip to positive value                  |
| 4556-4559     | VALIDATE  | `if (swman(imper).eq.2 .and. intwl(imper).lt.1)` fatalerr — intwl must be >= 1         |
| 4561-4568     | READ      | Count nrman1 / nrman2 based on swman; validate swman not out-of-range                    |
| 4565-4568     | VALIDATE  | `if (swman(imper).ne.1 .and. .ne.2)` fatalerr — swman out of range                     |
| 4570-4573     | VALIDATE  | `if ((nrman1+nrman2).ne.nmper)` fatalerr — period count mismatch                        |
| 4576          | READ      | `call rdsinr('swqhr',1,2,swqhr)` — read type of discharge relationship                  |
| 4582          | READ      | `call rdsdor('sofcu',0.1,100000,sofcu)` — read surface area of control unit             |
| 4584          | READ      | `call rdfinr('imper_4c',1,nmper,imper_4c,mamp,nmper)` — read period indices             |
| 4585-4589     | READ      | Compute zb=min(zbotdr(1),zbotdr(2)); read hbweir array and alphaw, betaw arrays         |
| 4593-4594     | READ      | Init imperb=0, imperi=0 (loop control for 4c duplicate-check)                            |
| 4596          | NORMALIZE | `hbweir(imper) = hbweir(imper) - altcu` — weir crest level: absolute → relative        |
| 4598-4599     | NORMALIZE | `alphaw(imper) = alphaw(imper) * (8.64*100^(1-betaw)/sofcu)` — unit/formula normalization |
| 4600-4604     | VALIDATE  | `if (hbweir(imper).lt.zbotdr(1+nrpri))` fatalerr — weir crest below channel bottom      |
| 4607-4613     | VALIDATE  | `if (swman=1 .and. wscap>1e-7 .and. (hbweir-wldip).lt.(zbotdr+1e-4))` fatalerr — supply impossible below zbotdr |
| 4615-4621     | VALIDATE  | `if (imper.ne.imperb)` / else fatalerr — imper_4c must be unique                        |
| 4625-4628     | VALIDATE  | `if (imperi.ne.nmper)` fatalerr — part 4c record count must equal nmper                 |
| 4630-4716     | GUARDED   | `elseif (swqhr.eq.2)` block (part 4d) — q-h table discharge relation. Not in port scope. |
| 4637-4642     | GUARDED   | Init arrays nqh, flweir, exists, flzero (swqhr=2 only)                                  |
| 4645-4649     | GUARDED   | Read imper_4d, imptab, htab, qtab (swqhr=2 only)                                        |
| 4652-4655     | GUARDED   | Convert: `hqhtab(imper,itab) = hhtab(i) - altcu`; store qqhtab (swqhr=2 only)           |
| 4658-4700     | GUARDED   | Validate q-h table: nqh count, level above zbotdr, unique imper, hqhtab first=altcu+100, descending h/q (swqhr=2 only) |
| 4703-4706     | GUARDED   | `if (imperi.ne.nmper)` fatalerr — period count for 4d (swqhr=2 only)                   |
| 4709-4714     | GUARDED   | `do imper=1,nmper`: if qqhtab not going to zero fatalerr (swqhr=2 only)                 |
| 4720-4859     | GUARDED   | `if (nrman2.gt.0)` block (parts 4e1+4e2) — automatic weir management. Not in port scope. |
| 4727-4729     | GUARDED   | Read imper_4e1, dropr, hdepth arrays (nrman2>0 only)                                    |
| 4731-4755     | GUARDED   | Loop: hdepth=abs, find compartment node (nodhd), validate swman=2, unique imper (nrman2>0) |
| 4756-4759     | GUARDED   | `if (imperi.ne.nrman2)` fatalerr — 4e1 period count mismatch (nrman2>0)                |
| 4762-4768     | GUARDED   | Init nphase, imperb, imperi for 4e2 (nrman2>0)                                          |
| 4771-4778     | GUARDED   | Read imper_4e2, impphase, wlsman, gwlcrit, hcrit, vcrit (nrman2>0)                     |
| 4781-4786     | GUARDED   | Convert: `wlsman -= altcu`; store gwlcrit, hcrit, vcrit (nrman2>0)                     |
| 4789-4807     | GUARDED   | Validate 4e2: wlsman above zbotdr, swman=2 consistency, unique imper counting (nrman2>0) |
| 4809-4814     | GUARDED   | `if (imperi.ne.nrman2)` fatalerr — 4e2 period count mismatch (nrman2>0)                |
| 4817-4859     | GUARDED   | Consistency checks: gwlcrit(1)=0, vcrit(1)=0, hcrit(1)=0; hbweir within 1cm of wlsman; wlsman/gwlcrit/hcrit/vcrit monotonicity across phases (nrman2>0); closing `endif` + structural comment at 4857-4859 |
| 4860          | RUNTIME   | `numadj = 0` — init counter for target level adjustments                                |
| 4866-4867     | RUNTIME   | `sttab(1,1)=100.0; sttab(2,1)=0.0` — set top two levels of storage table               |
| 4868-4872     | RUNTIME   | Loop i=3,22: `sttab(i,1) = zbotdr(1+nrpri)*(i-2)/20.0` — divide depth to 20 compartments |
| 4878-4897     | RUNTIME   | Double loop i=1,22 / ilev=1+nrpri,nrlevs: compute sttab(i,2) as surface storage per unit area (open channels only, trapezoidal for below surface, rectangular contribution for ponding above) |
| 4900          | RUNTIME   | `swstini = swstlev(wls1)` — compute initial storage from initial water level             |
| 4901          | RUNTIME   | `swst = swstini` — set current storage = initial storage                                 |
| 4904-4906     | RUNTIME   | Loop i=1,4: `wlsbak(i) = 0.0` — zero-init 4-timestep water level memory array          |
| 4909          | READ      | `CLOSE(DRA)` — close swap.dra file                                                       |
| 4911          | READ      | `RETURN`                                                                                  |
| 4912          | READ      | `END` — end of subroutine                                                                 |

---

## Summary

### Lines per bucket

| Bucket   | Approx. line-range count | Description                                                     |
| -------- | ------------------------ | --------------------------------------------------------------- |
| READ     | ~212 lines               | `rd*` calls, file open/close, local declarations, use statements, array init for reading |
| VALIDATE | ~85 lines                | `fatalerr` / `warn` calls and surrounding `if` logic            |
| NORMALIZE| ~20 lines                | Unit and coordinate conversions (altcu subtraction, cm conversion, alphaw formula, abs(wldip)) |
| RUNTIME  | ~33 lines                | sttab build, swstini/swst init, wlsbak zero-init, wlstar init, nrpri/nrsec assignment, nrman1/nrman2 init |
| GUARDED  | ~184 lines               | swsrf=3 primary wlp table, swsec=1 prescribed secondary wl, swqhr=2 q-h table, nrman2>0 automatic weir (incl. closing structural lines at 4857-4859) |

Total: 640 lines of substantive content (4273–4912, minus ~2 blank header lines at the top).

### Call count summary

- **READ:** ~40 `rd*` calls (rdsinr, rdsdor, rdfdor, rdfinr, rdatim, rdainr, rdftim, rdinqr, rdinar)
- **VALIDATE:** ~28 `fatalerr` calls + 2 `warn` calls
- **NORMALIZE:** 5 distinct conversion expressions (l×100, zbotdr−altcu, ztopdislay−altcu, hbweir−altcu, alphaw×unit-factor, wls1−altcu, wldip=abs, wlsman−altcu)
- **RUNTIME:** 5 initialization actions (wlstar, numadj, sttab levels, sttab storage, swstini/swst, wlsbak)
- **GUARDED:** 4 conditional branches: swsrf=3, swsec=1, swqhr=2, nrman2>0

---

## Key Findings for Implementation Tasks

### 1. altcu is the master coordinate reference

Every depth/level read from swap.dra is stored in absolute altitude (relative to some datum) and must have `altcu` subtracted to convert to soil-surface-relative coordinates. This affects: `zbotdr`, `ztopdislay`, `hbweir`, `wlp`/`wlptab`, `wls`/`wlstab`, `wlact` (→`wls1`), `wlsman`, `hqhtab`. The TOML path must either apply these conversions in the finalizer or store the altcu-relative value directly in the TOML.

### 2. alphaw conversion depends on sofcu and betaw

The alphaw normalization at lines 4598-4599 (`alphaw *= 8.64 * 100^(1-betaw) / sofcu`) is a compound formula that must be reproduced exactly in the finalizer. It is NOT a simple unit scale — it folds in sofcu (surface area of control unit). This means sofcu must be finalized before alphaw, or both must be finalized together in the same pass.

### 3. wls1 is initialized differently for swsec=1 vs swsec=2

- `swsec=1`: `wls1 = afgen(wlstab, 2*mawls, t1900-1.0)` — interpolation from a prescribed table (GUARDED).
- `swsec=2`: `wls1 = wlact - altcu` — read from wlact field (port scope). The `surfacewater_init` module only needs to handle swsec=2.

### 4. wlp1 is only computed for swsrf=3 (GUARDED)

`wlp1 = afgen(wlptab, 2*mawlp, t1900-1.0)`. The port targets swsrf=2; wlp1 is always zero/unused. The `surfacewater_init` signature should still accept wlp1 as an output for API compatibility, but can set it to zero or skip it.

### 5. swsrf=1 causes an early close+return at line 4463-4464

The validator for swsrf must catch swsrf=1 as a "no surface water" path. The port path (swsrf=2) never hits the early return. The stub should emit a clean error for swsrf=1 if called on the new path. The early-return skip also means none of swsec/swqhr/nmper are read when swsrf=1 — validators must reflect this gating.

### 6. sttab geometry — key invariant

The storage table `sttab(22,2)` is built using `l(ilev)` which was already converted to cm (×100) at line 4399. If `l` is stored in meters in the new TOML config (raw value), the `surfacewater_init` module must apply the ×100 factor before using it in sttab computation, or the finalizer must store l in cm.

### 7. nrman2>0 block is the most complex guarded section

Lines 4720-4859 (~140 lines) cover multi-phase automatic weir management. This is deeply intertwined with the numnod/dz soil column state (line 4738-4742 determines `nodhd` by walking the soil nodes). This cannot be stubbed with a simple "not supported" message in the config validator alone — it requires runtime node-depth data unavailable at config-read time. The stub at validator level should refuse swman=2 entirely.

### 8. zb = min(zbotdr(1),zbotdr(2)) dependency

At line 4585, `zb = min(zbotdr(1), zbotdr(2))` is computed before reading hbweir. This means the hbweir read-range bound depends on already-converted zbotdr values. In the TOML path, hbweir is read directly as altcu-relative and validated post-finalize, so this dependency is resolved naturally.

### 9. imper uniqueness checks use a running imperb variable

Parts 4c, 4d, 4e1, 4e2 each use an `imperb` running tracker to verify that management period indices are unique and sequential. This pattern must be replicated in the validator for swqhr=1 (part 4c) when implementing Task 5.

### 10. Double-loop storage computation (lines 4878-4897) handles ponding

The sttab(i,2) computation branches on `sttab(i,1) <= 0` (below surface) vs above. When the water level is above the soil surface (ponding, sttab level = 100 cm), the channel is treated as a rectangle of width `wbreadth` above the trapezoid. This must be reproduced exactly in `surfacewater_init`.
