# Reference binaries — SWAP 4.2.0

Pre-compiled binaries of the unmodified upstream SWAP 4.2.0 implementation. Used by regression tasks that compare the modernization's output against the reference.

## Files

| File | Target | Built | Notes |
|---|---|---|---|
| `swap420` | Linux x86-64, dynamically linked (glibc 2.6.32+) | 2021-10-06 | Intel Fortran compile from upstream. ELF 64-bit LSB, not stripped. |
| `swap.exe` | Windows x86-64 | 2021-10-06 | Intel Fortran compile from upstream. |
| `swap420gf` | Linux x86-64, **statically linked** | 2026-05-27 | **gfortran** compile of the *same* pristine 4.2.0 source, using the modern build's flags. Self-contained. **This is the regression reference.** |

### Why a gfortran build of 4.2.0?

`swap420` is Intel-compiled; the modern build is gfortran. To make the regression a
*physics* comparison rather than an *Intel-vs-gfortran* comparison, `swap420gf` rebuilds
the **unmodified** 4.2.0 source with the modern build's exact flags
(`-O2 -ffree-line-length-none -std=legacy -finit-local-zero`). On every existing case the
Intel and gfortran 4.2.0 builds agree to the fixture's 2-decimal precision (verified
2026-05-27), confirming the compiler is not a source of drift — so any modern-vs-`swap420gf`
divergence is a genuine code difference. **No source edits were needed** to compile 4.2.0
under gfortran.

Build recipe: `build_swap420gf.sh` (run after `git worktree add /tmp/swap-legacy legacy/swap-4.2.0`).

## License

SWAP 4.2.0 is distributed under the **GNU General Public License Version 2**. A small number of files (TTUTIL427.LIB) are under LGPL v2.1. See the full text in the `legacy/swap-4.2.0` orphan branch (`source_swp_4.2.0/LICENSE`).

## Source

The full source tree that produced these binaries is preserved as the `legacy/swap-4.2.0` orphan branch in this repository. To inspect it:

```bash
git worktree add /tmp/swap-legacy legacy/swap-4.2.0
```

## Usage

A 4.2.0 binary reads **legacy ASCII** inputs (`.swp/.crp/.dra`), so it must be run with
`run_case.sh --legacy-binary` (which `cd`s into the legacy `<N>.<case>/` dir), **not**
`--exec` (which selects TOML mode — wrong for a 4.2.0 binary):

```bash
cd tests/swap-cases
./run_case.sh -c hupselbrook --legacy-binary ../../tests/reference/swap420gf -k
```

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
