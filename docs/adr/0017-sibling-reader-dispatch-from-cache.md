---
title: "ADR 0017 — Sibling-reader dispatch from per-rotation cache"
date: 2026-05-02
status: accepted
---

# ADR 0017: Sibling-reader dispatch from per-rotation cache

## Context

The `.crp` port (Phases 1-4 of Phase 4f, planned in 2026-05-02) replaces the
legacy `readcropfixed`/`readwofost`/`readgrass` readers with a typed-pipeline
+ runtime-init pattern (ADR 0016: per-rotation cache). The original Phase 1
plan assumed a single legacy reader per crop type, called once at the
per-rotation init in `cropgrowth.f90`. The wiring task replaced that one
call with a dispatch on `crop_config_global%rotation_loaded(icrop)`:
cache-hit calls a typed-config init module; cache-miss falls back to the
legacy reader.

Phase 1 Task 8's smoke test (renaming `grass.crp` and re-running regression)
revealed an unspoken assumption: the per-crop-type runtime path actually
calls **multiple** legacy readers from the same `.crp` file, each scoped to
a sub-aspect of the crop:

- `readcropfixed` — main scalars + tables (read by `cropfixed(task=1)`).
- `readarablelandgerm` — preparation/sowing/germination switches (read by
  `ArableLandGerm(task=1)`, which runs in `CropGrowth(task=1)` *before*
  `cropfixed(task=1)`).
- `irrigation` — per-crop scheduling (read by `irrigation(task=1)`,
  conditionally called from inside legacy `readcropfixed` when `schedule=1`).

Each of these openers does the same boilerplate: build the path
(`pathcrop // cropfil(icrop) // '.crp'`), open the file, parse a subset of
keys, close. We call them **sibling readers**: they share the file but
parse disjoint key sets.

The Task 8 smoke test failed because removing the `.crp` from disk broke
`readarablelandgerm`'s open even though `cropfixed`'s open had been
correctly bypassed. The pattern repeats for every sibling reader.

## Decision

When a sibling reader is discovered to open the same legacy file we are
porting, apply the **same cache-driven dispatch** at its callsite that we
applied to the primary reader. Do NOT re-open the `.crp.toml` from disk a
second time; consume the already-parsed `cropfixed_config_t` (or
analogous typed config) from `crop_config_global%rotation_*(icrop)`.

Specifically:

1. **Extend the typed config schema** to cover the sibling reader's fields
   (if not already covered by the schema-1:1 decision in ADR 0015 +
   Phase 1's spec). Per the schema-1:1 philosophy, the relevant fields are
   usually already present as members of `cropfixed_config_t` /
   `cropwofost_config_t` / `cropgrass_config_t`, just unused by the primary
   init module.
2. **Extend the runtime-init module** (or add a sibling sub) to copy those
   fields to the legacy module globals the sibling reader writes to.
3. **Dispatch at the sibling callsite**: replace
   `call <sibling_reader>(...)` with the same `block`-scoped guard used in
   the primary dispatch:
   ```fortran
   block
      use crop_config_global_mod, only: crop_config_global
      use <init_mod>, only: <init_sub>
      logical :: use_cache
      ! same associated/allocated/bounds/loaded chain
      if (use_cache) then
         call <init_sub>(...)
      else
         call <legacy_sibling_reader>(...)   ! transitional fallback
      end if
   end block
   ```
4. **Document the sibling in the runtime audit**. The audit table for the
   primary reader already exists; add a one-line cross-reference for the
   sibling: "ALSO: `readarablelandgerm` (readswap.f90:3248) reads
   `swPrep/swSow/swGerm` from the same .crp; bypassed by cache dispatch
   in Phase 1 Task 7-fix."

## Consequences

**Positive:**

- Each sibling-reader fix is mechanical: same pattern, different reader,
  different fields. No new architecture per sibling.
- The cache loaded by `read_crop_toml.f90` is the single source of truth.
  No data lives in two places.
- The teardown story is the same: when the legacy fallback is removed
  (end of Phase 4 of the relevant `.crp` phase), all the sibling-reader
  `else` branches go away in lock-step with the primary.

**Negative:**

- The sibling discovery is reactive: the smoke test (rename file,
  re-run regression) is what surfaces them. This means the audit task at
  the start of each `.crp` phase needs to grep the codebase for ALL
  `pathcrop // ... // '.crp'` opens, not just inside the primary reader.
  Update the audit checklist accordingly.
- Each sibling fix is a small commit (TDD: failing smoke → extend init →
  dispatch → smoke green). For Phase 1, this is the
  `readarablelandgerm`/`ArableLandGerm(1)` patch that this ADR documents.
  For Phases 2 and 3, expect equivalent siblings (e.g.
  `irrigation(task=1)` if the case authors `schedule=1`).

**Neutral:**

- The init module gains sibling-specific responsibilities. As long as the
  module's docstring lists which legacy readers it replaces, this stays
  comprehensible.

## When to apply this ADR

Whenever the Task-8-style smoke test (or any later regression) reveals a
`File does not exist` error from a `pathcrop // ... // '.crp'` open after
the primary reader's dispatch has been wired. The signature is: the
runtime crashes during simulation init with `FOPENG: File does not exist`
+ a path ending in `.crp`. The fix is mechanical (per the Decision above),
not architectural.

## Pre-flight audit checklist

For each `.crp` port phase, the audit task (Task 1 in Phase 1's plan) must
include:

```bash
grep -n "trim(pathcrop)\|'.crp'" src/io/readswap.f90 src/crop/*.f90 \
   | grep -v "outfil\|outputfile" | sort -u
```

Every open found in this grep must either be:
- The primary reader (replaced by the runtime-init module).
- A sibling reader (dispatched at its callsite per this ADR).
- Already known to be reached only via a stub-errored switch (so it
  cannot run when the cache is loaded — verify the gate, not just the
  intent).

## Pattern reference

First documented application: Phase 1 of the .crp port, Task 7-fix. The
primary dispatch in `cropgrowth.f90:397` (cropfixed) was extended with a
sibling dispatch in `cropgrowth.f90:92` (`call ArableLandGerm(1)`) covering
the `readarablelandgerm` open at `readswap.f90:3273`. The cropfixed runtime
init module gained the four `flCrop*` flag assignments derived from
`cfg%swprep/swsow/swgerm`.

Future expected applications:
- Phase 1 Task 7-fix2 (if `irrigation(task=1)` discovers schedule=1 in any
  case): dispatch around `irrigation(task=1)` callsite from
  `cropgrowth.f90` and `readcropfixed`.
- Phase 2 (cropwofost): predict similar sibling readers around `wofost`
  initialization. Audit grep needed.
- Phase 3 (cropgrass): same.
