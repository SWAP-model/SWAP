# Reference binaries — SWAP 4.2.0

Pre-compiled binary of the unmodified upstream SWAP 4.2.0 implementation, used by
the regression harness to compare the modernization's output against the oracle.

## Files

| File | Target | Built | Notes |
|---|---|---|---|
| `swap420gf` | Linux x86-64, **statically linked** | 2026-05-27 | **gfortran** compile of the pristine 4.2.0 source, using the modern build's flags. Self-contained. **This is the regression reference.** |

The former Intel-compiled `swap420` (Linux) and `swap.exe` (Windows) binaries were
removed 2026-07-05: the project is gfortran-only (ADR 0001), nothing in the build or
tests consumed them, and their only references were four already-dead pixi run tasks.
They remain recoverable from git history and rebuildable from the `legacy/swap-4.2.0`
branch if ever needed.

### Why a gfortran build of 4.2.0?

The modern build is gfortran, and the original 4.2.0 reference was Intel-compiled. To
make the regression a *physics* comparison rather than an *Intel-vs-gfortran* one,
`swap420gf` rebuilds the **unmodified** 4.2.0 source with the modern build's exact flags
(`-O2 -ffree-line-length-none -std=legacy -finit-local-zero`). On every case the Intel
and gfortran 4.2.0 builds agreed to the fixture's 2-decimal precision (verified
2026-05-27), confirming the compiler is not a source of drift — so any
modern-vs-`swap420gf` divergence is a genuine code difference. **No source edits were
needed** to compile 4.2.0 under gfortran.

Build recipe: `build_swap420gf.sh` (run after `git worktree add /tmp/swap-legacy legacy/swap-4.2.0`).

## License

SWAP 4.2.0 is distributed under the **GNU General Public License Version 2**. A small number of files (TTUTIL427.LIB) are under LGPL v2.1. See the full text in the `legacy/swap-4.2.0` orphan branch (`source_swp_4.2.0/LICENSE`).

## Source

The full source tree that produced these binaries is preserved as the `legacy/swap-4.2.0` orphan branch in this repository. To inspect it:

```bash
git worktree add /tmp/swap-legacy legacy/swap-4.2.0
```

## Usage

`swap420gf` reads **legacy ASCII** inputs (`.swp/.crp/.dra`). The regression harness and
`regen_reference.py` invoke it directly against the legacy case tree (case *inputs* live
in the `swap-testcases` sibling repo; see `tests/regression/`), so there is normally no
need to run it by hand.

To (re)generate the regression reference fixtures from `swap420gf` for every registered
case:

```bash
python3 tests/regression/regen_reference.py            # all cases
python3 tests/regression/regen_reference.py hupselbrook # one case
```

This writes `<case>_reference_gf.json` next to the other fixtures; the regression harness
compares the modern build's output against these.

## Portability notes

`swap420` is dynamically linked and built against an older glibc. On modern Linux distributions it should still run but requires the usual shared libraries (`libm`, `libpthread`, `libc`). If it fails to start, inspect `ldd tests/reference/swap420` to see which libraries are unresolved. Static relinking is out of scope here; rebuild from `legacy/swap-4.2.0` if a different target is needed.
