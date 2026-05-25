---
title: Dependency management
author: SWAP modernization team
---

# Dependency management

## Two layers

The project's dependencies split into two layers, managed by different tools:

1. **pixi** manages the runtime toolchain — Python interpreter, Meson build
   driver, regression test harness, Fortran formatter, documentation
   generator. One-command install (`pixi install`) on every developer
   workstation; cross-platform by declaration (`platforms = ["linux-64",
   "win-64"]` in `pixi.toml`). The toolchain is pinned via `pixi.lock`, so
   every contributor solves to the same versions.
2. **meson subprojects** handle Fortran library dependencies that ship as
   source. The `.wrap` files under `subprojects/` point at the pinned
   upstream sources; Meson fetches and builds them as part of the main
   configure step.

The two layers are orthogonal: pixi does not know about subprojects, and
Meson does not know about pixi. Keeping them separate means the Fortran
build is compiler-environment-agnostic and the toolchain is reproducible
without system-level admin.

## pixi dependencies

All entries below are read directly from `pixi.toml`. Roles are the
observable use inside the repository, not marketing descriptions.

### Base `[dependencies]`

| Package | Version pin | Role |
|---|---|---|
| `python` | `3.11.*` | Interpreter used by Meson and by ad-hoc scripts in the repo. |
| `meson` | `==1.7` | Build driver; invoked via `meson setup builddir` / `meson compile -C builddir`. |

Note: `gfortran` is **not** in `pixi.toml`. The build task sets
`FC=gfortran` (see `_configure` in `[tasks]`) and expects `gfortran` to be
available on `PATH`. This is deliberate on Linux — pixi-provided compiler
toolchains have historically added friction on mixed-distro developer
setups — but it means a bare developer machine needs a system `gfortran`
before `pixi run build-linux` will succeed. See ADR 0001 on the
gfortran-first stance.

### `[pypi-dependencies]` (base)

| Package | Version pin | Role |
|---|---|---|
| `fprettify` | `>=0.3.7, <0.4` | Fortran source formatter; invoked via `pixi run lint`. |

### `test` feature

Activated by `pixi run -e test …`.

| Package | Version pin | Role |
|---|---|---|
| `python` | `3.11.*` | Regression harness runtime. |
| `pytest` | `>=8.0` | Test-runner library; used by the regression harness. |
| `pytest-xdist` | `>=3.5` | Parallel test execution across regression cases. |
| `pytest-timeout` | `>=2.2` | Per-test timeout enforcement; prevents a hung case from blocking the gate. |
| `pandas` | `<3.0` | CSV aggregation inside `tests/regression/test_output_regression.py`. |
| `pyswap` (pypi) | `>=0.3.9` | Regression harness dependency (case-file parsing / SWAP output helpers). |

### `docs` feature

Activated by `pixi run -e docs …`.

| Package | Version pin | Role |
|---|---|---|
| `ford` | `>=7.0` | Generates Fortran API reference from source comments. |
| `python` | `3.11.*` | FORD runtime. |

## Meson subprojects

Each entry below corresponds to one file under `subprojects/`. The
populated subdirectories (`subprojects/ttutil/`, `subprojects/toml-f/`,
`subprojects/test-drive/`) are present after a Meson configure step —
they are the checked-out copies Meson fetched, and are tracked as
untracked working-tree noise rather than committed code.

### `ttutil`

`subprojects/ttutil.wrap`:

```
[wrap-git]
directory = ttutil
url = https://github.com/SWAP-model/ttutil.git
revision = v4.2.7
```

Legacy utility library carried forward from the original SWAP
distribution and now hosted at `SWAP-model/ttutil`. Provides helpers the
core Fortran code depends on at configure time. Required by the main
`swap` binary. Built as a static library from source under
`subprojects/ttutil/` whenever `meson setup builddir` runs. Pin is a
named tag (`v4.2.7`), so `meson subprojects update` is deterministic.

### `toml-f`

`subprojects/toml-f.wrap`:

```
[wrap-git]
directory = toml-f
url = https://github.com/toml-f/toml-f
revision = head
```

Pure-Fortran TOML parser from `https://github.com/toml-f/toml-f`.
Required by the TOML configuration loaders in the `src/io/` tree (the
`.toml` variants of the readers). Chosen over YAML-based alternatives
because it is pure Fortran — no C, no Python runtime — for config
parsing inside the Fortran binary.

**Flag:** `revision = head` is **not a pin.** It resolves to whatever
`HEAD` is on the upstream default branch at the moment Meson fetches.
This is a reproducibility hazard and should be replaced with a tagged
revision as part of Phase 3 / dependency hygiene. Left as-is here to
document the current state, not to endorse it.

### `test-drive`

`subprojects/test-drive.wrap`:

```
[wrap-redirect]
filename = toml-f/subprojects/test-drive.wrap
```

The top-level wrap is a **redirect** to toml-f's own wrap at
`subprojects/toml-f/subprojects/test-drive.wrap`, which pins
`revision = v0.4.0` from `https://github.com/fortran-lang/test-drive.git`.

Lightweight Fortran unit-test helper. Present transitively because
toml-f's own test suite uses it. pFUnit is SWAP's primary test harness;
test-drive is vendored as a side effect of toml-f and is not used by any
live SWAP tests today.

## pFUnit

pFUnit sits outside `subprojects/` and is treated separately. It lives at
`tests/pFUnit` as a **gitlink** (mode `160000`, SHA
`581cd32be937df1922ce50fb06ddbc5522e6ce66` at the time of writing), not
as a Meson wrap and not as a configured Git submodule. The outer repo
records the SHA but there is no `.gitmodules` entry to re-clone it
automatically.

See `docs/build-and-test.md` for the gitlink peculiarity, the
rebuild-into-`install_gfortran/PFUNIT-4.15` step, and the
`PFUNIT_ROOT` activation variable that `pixi.toml` exports. Bumping
pFUnit is a manual rebuild-and-re-pin (see below), not a wrap-file bump,
precisely because it is not a wrap.

## Bumping versions

### Bumping a pixi dep

1. Edit `pixi.toml` (base `[dependencies]` or the relevant
   `[feature.<name>.dependencies]` block).
2. Run `pixi install` (base) or `pixi install -e test` / `-e docs` to
   resolve and update `pixi.lock`.
3. Run `pixi run -e test check-full` and confirm green.
4. Commit **both** `pixi.toml` and `pixi.lock` in the same commit.

### Bumping a subproject

1. Edit `subprojects/<name>.wrap` — change the `revision = …` line to the
   new tag or commit SHA. For `toml-f` specifically, the first bump
   should be `head` → a named tag (see Flag above).
2. Remove the checked-out copy so Meson re-fetches:
   `rm -rf subprojects/<name>/`.
3. Reconfigure: `meson setup builddir --reconfigure`.
4. Run `pixi run -e test check-full` to catch any behaviour change.
5. Commit the `.wrap` file.

### Bumping pFUnit

1. Rebuild per `docs/build-and-test.md` — either from the existing local
   clone at a new tag, or from a fresh clone alongside the repo.
2. Stage the new gitlink SHA: `git add tests/pFUnit`.
3. Run `pixi run -e test check-full` (which runs `test-pfunit`).
4. Commit.

Converting pFUnit into a proper Git submodule or a Meson subproject is a
deferred chore — not required to bump the pin.

## Dependency discipline during the rescue

One variable at a time. New dependencies are **not** added during the
rescue (Phases 1–4) unless they are load-bearing for the phase's exit
criterion. Proposals for new deps wait until after rescue/complete.

The rescue is already managing a large number of moving pieces:
compiler migration, builddir consolidation, fixture regeneration, test
reconciliation, and pFUnit re-population. Adding new deps mid-rescue
compounds the matrix of things that can break and makes bisection
harder. Remove-or-replace changes to existing deps are likewise deferred
unless the current version actively blocks the phase.
