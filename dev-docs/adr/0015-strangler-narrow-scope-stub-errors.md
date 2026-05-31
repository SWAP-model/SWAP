---
title: "ADR 0015 — Strangler narrow-scope stub-errors for unsupported branches"
date: 2026-05-02
status: accepted
---

# ADR 0015: Strangler narrow-scope stub-errors for unsupported branches

## Context

The Phase 4f port effort replaces legacy ASCII readers (`swap.ini`, `swap.dra`,
`*.crp`, …) with the typed-config + adapter pipeline (ADR 0007). Each legacy
reader spans a wide schema with many branches the user can author via legacy
switches (`swsrf=3`, `swsec=1`, `swqhr=2`, `swirfix`, multiple `swcrop` modes,
…). The TOML test cases exercise only a small subset of those branches:

| Reader     | Branches in case set                | Branches in legacy reader        |
| ---------- | ----------------------------------- | -------------------------------- |
| `swap.ini` | `swinco=3`                          | `swinco ∈ {1,2,3}` + heat / solute combinations |
| `swap.dra` | `swdra=2 + swsrf=2 + swqhr=1`       | `swdra ∈ {0,1,2}` × `swsrf ∈ {1,2,3}` × `swqhr ∈ {1,2}` × `swman ∈ {1,2}` × … |
| `*.crp`    | TBD (pending ADR-bound brainstorm)  | three crop modes (simple / detailed / WOFOST / grass) × many sub-modes |

Two strategies for handling the gap between authored test coverage and the
legacy schema's full breadth:

1. **Full coverage upfront**: extend the typed schema, the validator, the
   finalizer, the adapter, and the runtime path for every legacy branch
   before deleting any legacy code. Maximises future-proofing; very expensive
   per port; pulls in plumbing work (e.g. the altcu coordinate-system
   threading through `zbotdr`, `hbweir`, and `wls1`) that may never be
   exercised by a test case.

2. **Narrow scope with stub-errors**: port only what the test cases exercise.
   For unsupported branches, add a stub-error in the appropriate validator
   that rejects the configuration with a clear message ("not yet supported in
   the TOML pipeline; use the legacy executable"). Defer the plumbing until
   a case authors the branch.

The swap.ini and swap.dra ports both used strategy (2). It worked: each port
shipped in 1-2 sessions, regression stayed green, and the deferred branches
were trivially identifiable from the validator stub-error messages.

## Decision

For Phase 4f legacy-reader ports, the default is **strategy (2): narrow
scope with stub-errors**.

A port may extend to additional branches when:
- A TOML test case authors the branch, or
- A clear, near-term consumer (downstream Python bindings, an active
  external user) needs the branch, or
- The cost of adding the branch is comparable to the cost of writing the
  stub-error.

Otherwise, the port:

- Extends the typed config, validator, finalizer, adapter, and runtime
  module(s) for the branches the test cases exercise.
- Adds explicit stub-error validators for unsupported branches, rejecting
  them at validate time with `ERR_VALIDATION_CROSS_FIELD` and a message of
  the form: `"<section>.<field>=<value> not yet supported in the TOML
  pipeline; use the legacy executable."`
- Adds defense-in-depth runtime guards (`fatalerr_collected`) inside any
  new runtime module that mirrors the validator stub-errors.
- Keeps the legacy reader (e.g. `rddre`, `rddrb`, the `*.crp` sub-readers)
  alive in `readswap.f90` / `cropgrowth.f90` for the parity test fixtures.
  The reader is never called from the runtime path; only the parity tests
  exercise it.
- Documents the deferred branches in the port's design spec
  and / or the migration history section of `docs/csv-companion-files.md`.

## Consequences

**Positive:**

- Each port ships predictably in 1-2 sessions instead of 4-6.
- The validator surface becomes the documentation of what's supported. A
  new contributor can `grep "not yet supported" src/config/` to find every
  deferred branch and its location.
- Speculative plumbing for branches no one uses doesn't exist, so it
  doesn't bit-rot or carry hidden bugs.
- The legacy executable remains the canonical fallback for users with cases
  outside the supported subset. No functionality is lost — only the entry
  point is constrained.

**Negative:**

- Users authoring a TOML case with a deferred branch see a stub-error and
  must either (a) refactor their case to fit the supported subset, (b) use
  the legacy executable, or (c) wait for someone to extend the port. This
  is a real friction cost that grows with deferred surface area.
- The TOML executable is "almost feature-complete" rather than fully
  feature-complete. Documentation must be precise about the supported
  subset to avoid misleading users.
- A future port that wants to lift a stub-error must do focused
  archaeology — figuring out what plumbing the legacy reader applied
  inside its branch (coordinate normalizations, cross-section validations,
  runtime initializations) and where it should live in the typed
  pipeline. This is the work the original narrow-scope decision deferred.

**Neutral:**

- The migration history of each port records which branches were deferred,
  giving future implementers a clear backlog if they want to broaden
  coverage. The swap.dra port deferred `swsrf=3 / swsec=1 / swqhr=2 /
  swman=2 / altcu /= 0 / dramet /= 0`. The swap.ini port deferred nothing
  (it covered its full schema).

## Pattern reference

Implementations of this ADR:

- `src/config/surface_water_config.f90` — stub-errors for `swsrf=3`,
  `swsec=1`, `swqhr=2`, `swman=2` (Phase 4f swap.dra port).
- `src/config/drainage_config.f90` — stub-errors for `altcu /= 0` and
  `swdra=2 + dramet /= 0` (Phase 4f swap.dra port).
- `src/drainage/surfacewater_init.f90` — defense-in-depth runtime guards
  matching the validator stub-errors above.

When porting `*.crp`, expect the brainstorm to surface a similar table:
which crop modes / sub-modes does the test set exercise? What are the
remaining branches, and at what level (mode / sub-mode / per-table) do
they get stub-errored? The answer shapes the port's scope.
