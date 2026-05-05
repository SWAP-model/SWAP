---
title: Branch convention
date: 2026-05-05
---

# Branch convention

This repository hosts two parallel histories on the remote:

| Branch | What it holds | Consumers |
|---|---|---|
| `main` | **Legacy SWAP 4.2.0 reference implementation.** Frozen against accidental modernization changes. | `pyswap` (hardcoded `git clone main` for build artifacts); GitHub Actions legacy build workflow. |
| `modern` | **Modernized runtime: TOML-only input pipeline, typed config schema, CSV outputs.** Active development line. | Internal modernization work; future consumers once they migrate. |

## Why two branches

The modernization (Phase 4 / Phase 4f-extend, completed 2026-05-05;
see [`PHASE-4-MODERNIZATION-SUMMARY.md`](PHASE-4-MODERNIZATION-SUMMARY.md))
landed ~460 commits of changes that retire the legacy fixed-format
reader chain. External consumers — `pyswap` and the GitHub Actions
legacy artifact build — depend on `main` pointing at the SWAP 4.2.0
codebase. Pushing the modernized history to `main` would break those
consumers without warning.

`modern` exists so the modernization work has a remote home (backed
up, browsable, taggable) without disturbing the consumers of `main`.

## Day-to-day workflow

```
development        ← short-lived feature work; merge upward
   │
   ▼ git merge --no-ff
modern             ← active modernization branch (push here, never to main)
   │
   ▼ (eventual cutover, see below)
main               ← currently legacy v4.2.0; future home of the modern branch once consumers migrate
```

- New work happens on `development` (or topic branches).
- When stable, merge into `modern` with `--no-ff` (preserves the
  feature branch in history).
- Push `modern` to `origin/modern`.
- Never push directly to `origin/main`.

## Tags

The `rescue/*` and `docs/*` tags annotate the modernization journey
on the `modern` branch. The capstone summary
([`PHASE-4-MODERNIZATION-SUMMARY.md`](PHASE-4-MODERNIZATION-SUMMARY.md))
lists each tag with its scope.

## Future plan: fork into a separate repository

Once the modernization is independently consumable (Python bindings,
release artifacts, documentation site), the `modern` branch will be
forked into a new repository — separating "modern SWAP" cleanly from
the SWAP 4.2.0 legacy reference. At that point:

- New repo: hosts `modern` as its `main`.
- This repo: continues hosting the legacy 4.2.0 reference for
  archival / pyswap-style consumers.

Until then, the two-branch arrangement keeps both histories live in
one place.

## Local convention

After the initial setup:

```
local/main       ← mirrors origin/main (legacy)
local/modern     ← tracks origin/modern (modernization)
local/development ← topic branches feed into modern
```

`git fetch && git pull --ff-only` on either branch is safe — neither
diverges from its remote.
