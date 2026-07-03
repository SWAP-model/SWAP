---
title: "Arc T1-B — dead-code & dormant-clarity sweep (spec + plan)"
date: 2026-07-03
status: in-progress
tags: [tier1, dead-code, dormant, hygiene, byte-identical]
parent: 2026-07-03-modernization-review.md
---

# Arc T1-B — dead-code & dormant-clarity sweep

First Tier-1 arc from the [modernization review](2026-07-03-modernization-review.md).
The review's §4 estimated "~9 K lines of dead code to remove." **Discovery
revised that down sharply, and that revision is itself the arc's main value:**
most of the "dead" surface is either intentional or actually live.

## Findings from discovery (verified, not reported)

1. **The `dormant/` modules (~9 K lines) are a deliberate museum, not dead
   code.** All three `dormant/` dirs carry a `README.md` documenting each
   module and its reactivation checklist; meson has tombstone comments pointing
   at them; they are build-excluded (zero binary cost). The `crop/dormant/README.md`
   even has a "Why not delete?" section. **Disposition: keep.** Deleting them
   would destroy curated future-work seeds (the `legacy/swap-4.2.0` branch
   preserves only the *raw* original, not the modern-integrated restoration
   notes).

2. **The museum docs have drifted.** `crop/dormant/README.md` still lists
   `jongvanlier.f90`, which ADR 0047 (W12) deleted. Restoration steps across all
   three READMEs still say "drop the `use variables` blocks" / "restore the
   retired `variables.f90` declarations" — but `variables.f90` was deleted; the
   modern targets are `state%…` records. **Disposition: refresh to match reality.**

3. **10 of the 12 "dead state fields" in review §4 are LIVE.** A full-tree grep
   (reads *and* writes, including `associate` aliases, output writers, bindings)
   found: `evp` read at `soilhydraulics.f90:101`; `runonarr`/`flrunon` at
   `boundtop.f90:53-54`; `alfaw_layer` at `soilhydraulics.f90:1349`; `ksatexm`
   read in ~8 sites (drainage, surfacewater, boundtop, frozencond,
   soilhydraulics(utils), tillage); `cofani` in drainage/surfacewater/frozencond;
   `flksatexm` at `tillage.f90:50`; `bpegwl` read in the `waterbalance.f90`
   perched-gwl loop; `psilt`/`pclay` in tillage + temperature. **Disposition:
   keep, and record here so they are not re-attempted.**

4. **Only two fields are genuinely dead:** `pegwl` and `npegwl`
   (`soilwater_state.f90:191,193`) — written in `waterbalance.f90` but read
   nowhere (not in compute, output, or bindings). Their sibling `bpegwl` (the
   node index) is the only live output of that block.

5. **`stepnr()`** (`arrayutils.f90:56-90`, exported at `:10`) — a complete public
   function with **no caller anywhere** in `src/` or `tests/`. Genuinely dead.

6. **The macropore "retired-zero placeholders" are mostly gone.** Post-strangler
   they are tombstone comments plus a couple of load-bearing always-zero locals
   that preserve byte-identity (`ArMpSs`), and `FrArMtrx` — which is *live*
   (read in the output writer, `csv_output.f90:843`). **Disposition: out of
   scope; nothing safely removable.**

## Scope

**In:**
- Delete `stepnr()` (definition + public export).
- Delete `pegwl` / `npegwl` fields and their write-only assignments.
- Refresh the three `dormant/README.md` files to match the actual files and the
  post-`variables.f90` reality.
- Record the "these 10 fields are LIVE" correction (this doc + review §4 update).
- Write an ADR formalizing the dormant-module preservation policy.

**Out:** deleting any `dormant/` module; macropore scaffolding; `level()` (kept —
live caller at `waterbalance.f90:87` computes `gwl`); the TOML/registry/binding
work (later arcs).

## Byte-identity analysis for the removals

- `stepnr` — no callers → removal cannot change output.
- `npegwl` writes (`waterbalance.f90:132,148,152`): RHS are `node` / `1` / `-1`
  (no side effects).
- `pegwl` writes (`:133,141,143,146`): RHS are `level(...)`, `min(...)`,
  `soil%pond`, `0`. `level()` takes `state` as `intent(in)` → **pure**, no
  mutation; it is still called at `:87` for the live `gwl`, so the function
  stays. Removing the `:133` call (whose result only fed the dead `pegwl`) is
  byte-identical *and* drops a ~59 K-calls/run pure evaluation as a minor bonus.
- The surrounding control flow (`flsat`, `node`, `nodheq1`, `bpegwl`) is left
  untouched, so the block still computes the live `bpegwl` bit-for-bit.

## Plan (each code step ends `check-fast` green; state edit ⇒ clean rebuild)

1. **stepnr** — drop from the `public::` list and delete the function body.
2. **pegwl/npegwl** — delete the two declarations in `soilwater_state.f90`; delete
   the seven assignments in `waterbalance.f90` (lines 132,133,141,143,146,148,152).
   Clean rebuild (`src/state/` change) + `check-fast` → byte-identical is the proof.
3. **Museum READMEs** — refresh all three (remove `jongvanlier`; retarget
   `variables.f90` restoration steps to `state%…`).
4. **ADR** — new ADR "dormant-module preservation pattern"; add to ADR index.
5. **Review-doc update** — correct §4 (dead-field list → the 2 real ones; note the
   museum is intentional).

Prerequisite (separate concern, its own commit): **T1-A CI gate** — a push/PR
workflow running `check-fast` + `test-pfunit`.
</content>
