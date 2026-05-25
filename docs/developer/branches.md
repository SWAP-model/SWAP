---
title: Branch convention
date: 2026-05-05
---

# Branch convention

| Branch | Role | Push target |
|---|---|---|
| `main` | **Modernization release line.** Only merge from `development` when a version is release-worthy. | Tagged releases only. |
| `development` | **Active day-to-day work.** Topic branches feed into here. | Pushed continuously. |
| `v4.2.0` tag | **Legacy SWAP 4.2.0 reference snapshot** (commit `7587ca3`). pyswap and other downstream consumers pin to this tag for the legacy artifact. | Frozen. |

## Workflow

```
topic-branch        ← optional short-lived feature work
    │
    ▼ git merge --no-ff (or fast-forward)
development         ← daily work; push freely to origin/development
    │
    ▼ git merge --no-ff   (only when a version is release-worthy)
main                ← release line; tagged when a version ships
```

- New work happens on `development` (or topic branches that merge into
  it).
- Push `development` to `origin/development` whenever you want a remote
  backup or want CI to run.
- **Do not merge `development` → `main` until a version is
  release-worthy.** Each merge to `main` should correspond to a
  tagged release.

## Why this convention

The modernization (Phase 4 / Phase 4f-extend, completed 2026-05-05;
see [`PHASE-4-MODERNIZATION-SUMMARY.md`](PHASE-4-MODERNIZATION-SUMMARY.md))
landed ~460 commits that retire the legacy fixed-format reader chain.
The first push of those commits accidentally landed on `origin/main`
on 2026-05-05 13:40 UTC+2; rather than reverting, we accepted `main`
as the new modernization line and pinned pyswap (and other consumers)
to the `v4.2.0` tag for legacy artifacts.

Going forward, `main` is reserved for releases. `development` carries
the in-progress work and is the default push target.

## Tags

The `rescue/*` and `docs/*` tags annotate the modernization journey.
The capstone summary
([`PHASE-4-MODERNIZATION-SUMMARY.md`](PHASE-4-MODERNIZATION-SUMMARY.md))
lists each tag with its scope. The `v4.2.0` tag points at the
pre-modernization legacy SHA (`7587ca3`) and is the canonical anchor
for downstream consumers that need the legacy reference.

## Future plan: fork into a separate repository

Once the modernization is independently consumable (Python bindings,
release artifacts, documentation site), the `main` branch will be
forked into a new repository — a "modern SWAP" repo separate from
this one (which becomes a legacy/archive holding pen). At that point:

- New repo: hosts modernization as its `main`.
- This repo: continues hosting the legacy 4.2.0 reference (via the
  `v4.2.0` tag) for archival / pyswap-style consumers.

Until then, both histories live here.

## Local convention

```
local/main          ← mirrors origin/main (modernization, release line)
local/development   ← mirrors origin/development (daily work)
```

`git fetch && git pull --ff-only` on either is safe. Merge from
`development` to `main` is a deliberate, infrequent action gated on
release-worthiness.
