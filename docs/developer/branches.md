---
title: Branch convention
date: 2026-05-05
updated: 2026-06-26
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
see `dev-docs/phase-4-modernization-summary.md` in the repository)
landed ~460 commits that retire the legacy fixed-format reader chain.
The first push of those commits accidentally landed on `origin/main`
on 2026-05-05 13:40 UTC+2; rather than reverting, we accepted `main`
as the new modernization line and pinned pyswap (and other consumers)
to the `v4.2.0` tag for legacy artifacts.

Going forward, `main` is reserved for releases. `development` carries
the in-progress work and is the default push target.

A short-lived `origin/modern` branch once existed to hold the
modernization line separately from a legacy `main`, on the assumption
that pushing to `main` would trigger executable rebuilds. That premise
no longer holds: CI (`.github/workflows/ci.yaml`) fires **only on
`v*` version-tag pushes**, never on branch pushes. `origin/modern` had
also drifted into an exact duplicate of `origin/main`, so it was
deleted on 2026-06-26. The active branch set is simply `main` +
`development`, with the `v4.2.0` tag as the legacy anchor.

## Tags

The `rescue/*` and `docs/*` tags annotate the modernization journey.
The capstone summary
(`dev-docs/phase-4-modernization-summary.md` in the repository)
lists each tag with its scope. The `v4.2.0` tag points at the
pre-modernization legacy SHA (`7587ca3`) and is the canonical anchor
for downstream consumers that need the legacy reference.

## On forking into a separate repository (deferred)

An earlier plan proposed eventually forking the modernization into a
new "modern SWAP" repository and leaving this one as a legacy/archive
holding pen for the 4.2.0 reference. As of 2026-06-26 that split is
**not planned**: this repository *is* the modern SWAP repo. The whole
tree is the modernization line; the only legacy artifact retained here
is the frozen `v4.2.0` tag, which downstream consumers (pyswap, etc.)
pin to directly. No separate legacy repo is maintained.

If a fork is ever revisited, it would be a deliberate future decision
recorded in a new ADR — not an assumed end-state of this convention.

## Local convention

```
local/main          ← mirrors origin/main (modernization, release line)
local/development   ← mirrors origin/development (daily work)
```

`git fetch && git pull --ff-only` on either is safe. Merge from
`development` to `main` is a deliberate, infrequent action gated on
release-worthiness.
