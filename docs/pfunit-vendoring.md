# pFUnit 4.15 vendoring — current state and rebuild steps

This is a **rescue-phase** note describing the baseline state of the pFUnit test-harness dependency as of commit `rescue/phase-1-infra`. Phase 2 (`docs/build-and-test.md`) absorbs this content into a broader build-and-test guide.

## Where pFUnit lives

`tests/pFUnit/` is a local checkout + build of the Goddard pFUnit 4.15 release. The CMake-built install tree at `tests/pFUnit/build/install_gfortran/PFUNIT-4.15` is what meson consumes — its path is surfaced via the `PFUNIT_ROOT` environment variable set in pixi's `[activation]` block.

## How the outer repo tracks it

**This is weird and worth knowing about.** `tests/pFUnit/` is a **gitlink** (git tree entry mode `160000`) pinned at commit `581cd32be937df1922ce50fb06ddbc5522e6ce66` in the outer repo, but it is **not** declared in `.gitmodules`. The result:

- `git clone` of the outer repo does NOT automatically populate `tests/pFUnit/`.
- `git submodule update` has no entry to act on.
- `git add tests/pFUnit/<anything>` fails with "Pathspec is in submodule".
- The working-tree contents are whatever was checked out when pFUnit was first dropped in.

This is a rescue-era inheritance, not a design decision. Converting it to one of: (a) a properly-declared submodule, (b) a meson subproject via `.wrap`, or (c) plain in-tree contents tracked by the outer repo — is **deferred** to Phase 1.5 or Phase 2. For now, the install exists on author machines and nothing active depends on populating it afresh.

## If the install is ever lost

Two options:

1. From the existing local clone:

       cd tests/pFUnit
       mkdir build && cd build
       cmake -DCMAKE_INSTALL_PREFIX=install_gfortran ..
       make -j
       make install

2. From scratch — clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit at tag `v4.15`, place it at `tests/pFUnit/`, then rebuild as above.

The meson config at `tests/unit/meson.build` expects the install at `tests/pFUnit/build/install_gfortran/PFUNIT-4.15` and errors with a pointer back to this doc if the path is missing.

## Why this hasn't been fixed during rescue Phase 1

Phase 1's MVP scope explicitly defers full pFUnit vendoring — the install works, Phase 4 testing isn't on the critical path yet, and the effort to convert to a proper meson subproject (pFUnit is CMake-based, needs a CMake-dependency wrap) is disproportionate to the remaining rescue work. Flagged for Phase 1.5 or Phase 2.
