---
title: "crop_config_global Retirement — Tier-3 sub-project 1"
date: 2026-05-27
status: approved
branch: development
parent: "Tier-3 crop-cluster decomposition (multi-instance readiness)"
---

# crop_config_global Retirement (Tier-3 sub-project 1)

## Context

Tier-3 (crop-cluster) work is decomposed into sequenced sub-projects (each its
own spec → plan → arc), because the dispatch retirement and the multi-instance
`SAVE`-elimination are coupled and total ~135+ shared-state variables across
~5500 lines:

1. **`crop_config_global` retirement** ← THIS SPEC (independent, mechanical)
2. Crop dispatch (`CropGrowth`/`CropFixed`/`Wofost`/`Grass`/`ArableLandGerm`)
   + cross-phase `SAVE` → `state%crop%{fixed,grass,wofost}`
3. `SoilManagement` dispatch + WSN nutrient module-`SAVE`
   (`wofost_soil_declarations` ~100 vars, `cropwofost_init_mod` 34 vars) →
   `state%nutrients`
4. Output file-handle `SAVE`s (nba/om1/om2) + read-only gauss constants

This sub-project removes the first named multi-instance blocker and establishes
the `state%cfg%crop` access pattern the dispatch arc (sub-project 2) will reuse.

## Problem

`crop_config_global` is a module-global pointer (`src/crop/crop_config_global.f90:23`,
`type(crop_config_t), pointer :: crop_config_global => null()`, in module
`crop_config_global_mod`, flagged transitional under ADR 0016). It is associated
once in `crop_state_init` (`src/state/crop_state.f90:85`:
`crop_config_global => crop_cfg`) and read at ~27–32 runtime sites across the
crop runtime files. As a module global reassigned per `crop_state_init`, it
breaks concurrent multi-instance execution: initialise instance A then B and the
pointer aims at B's config, so A's crop routines read B's rotations.

## Key facts (verified)

- **Identity:** `state%crop%init(config%crop, …)` (`swap_mod.f90:111`) makes
  `crop_cfg` == `config%crop`; with `state%cfg => config` (`swap_mod.f90:102`),
  `crop_config_global` and `state%cfg%crop` are the **same object**. Both are
  wired in `swap_init_body` before any crop runtime call → identical at every
  read site. Substituting `crop_config_global%X` → `state%cfg%crop%X` is
  semantically identical.
- **Read-only at runtime:** the only writes to `rotation_loaded` (the lazy-load
  cache flag) are in `src/io/toml/read_crop_toml.f90` at config-load time,
  through the real `config%`, never through `crop_config_global`. The runtime
  reads it into a local `use_cache`. No config mutation through the pointer.
- **Idiom present:** 4 of 5 reader files already open
  `associate (crop_cfg => state%cfg%crop)`:
  `cropgrowth.f90:89`, `cropfixed_runtime.f90:55`, `cropgrass_runtime.f90:101`,
  `cropwofost_runtime.f90:171`. Only `cropgrowth_helpers.f90` lacks it.

## Design

**Access form — reuse the existing `crop_cfg` associate alias.** In the 4 files
that already alias `crop_cfg => state%cfg%crop`, rewrite each
`crop_config_global%X` → `crop_cfg%X`, and rewrite any
`associate (… => crop_config_global%…)` target → `crop_cfg%…`. This matches the
established sub-record-associate idiom and keeps reads concise.

For `cropgrowth_helpers.f90` (3 reads in the `ArableLandGerm` prep/sow/germ
helper procedures, no existing alias): add `crop_cfg => state%cfg%crop` to each
helper's `associate` and use `crop_cfg%…`. (`state` is in scope at all three
sites.)

**Then delete the global:** once readers hit zero, remove
`crop_config_global => crop_cfg` (`crop_state.f90:85`), the
`crop_config_global_mod` module file (`src/crop/crop_config_global.f90`), and
every now-dead `use crop_config_global_mod` / `use …, only: crop_config_global`
import (including the one in `crop_state.f90`).

**Read sites by file (anchors; grep authoritatively at implementation time):**

| File | reads | existing `crop_cfg` alias |
|---|---|---|
| `src/crop/cropgrowth.f90` | ~16 (187–210, 356–357) | yes (`:89`) |
| `src/crop/cropgrass_runtime.f90` | 5 (118–131) | yes (`:101`) |
| `src/crop/cropfixed_runtime.f90` | 4 (73–80) | yes (`:55`) |
| `src/crop/cropwofost_runtime.f90` | 4 (189–196) | yes (`:171`) |
| `src/crop/cropgrowth_helpers.f90` | 3 (123, 153, 188) | no — add it |

## Scope guard

- No changes to `read_crop_toml.f90` config-load writes.
- No dispatch changes (those are sub-project 2).
- No other subsystem touched.
- The associate scope must cover each read site — verify the existing
  `crop_cfg` alias is in scope where reads occur; if a read sits outside the
  associate block, extend the block or spell out `state%cfg%crop%…` there.

## Sequencing & verification

One commit per reader file (5), each gated by byte-identical `pixi run -e test
check-fast` (4/4 clean; `soilhysteresis`/`winter` are known pre-existing xfails
not in the fast set). A final 6th commit deletes the module file + the
`=>` set-site + dead imports once `grep -rn crop_config_global src/` shows zero
code references. `check-full` at arc end (the 2 known xfails expected).

Because the pointer is read-only and identical to `state%cfg%crop`, every
commit is behavior-preserving by construction.

## Success criteria

- `grep -rn "crop_config_global" src/ --include=*.f90` → zero matches (no module,
  no pointer, no imports, no reads; comments may be cleaned but are harmless).
- `src/crop/crop_config_global.f90` deleted.
- `check-full`: 5 passed / 0 failed / 2 known xfails — zero new regressions.
