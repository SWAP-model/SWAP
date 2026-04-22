# Reference binaries — SWAP 4.2.0

Pre-compiled binaries of the unmodified upstream SWAP 4.2.0 implementation. Used by regression tasks that compare the modernization's output against the reference.

## Files

| File | Target | Built | Notes |
|---|---|---|---|
| `swap420` | Linux x86-64, dynamically linked (glibc 2.6.32+) | 2021-10-06 | Intel Fortran compile from upstream. ELF 64-bit LSB, not stripped. |
| `swap.exe` | Windows x86-64 | 2021-10-06 | Intel Fortran compile from upstream. |

## License

SWAP 4.2.0 is distributed under the **GNU General Public License Version 2**. A small number of files (TTUTIL427.LIB) are under LGPL v2.1. See the full text in the `legacy/swap-4.2.0` orphan branch (`source_swp_4.2.0/LICENSE`).

## Source

The full source tree that produced these binaries is preserved as the `legacy/swap-4.2.0` orphan branch in this repository. To inspect it:

```bash
git worktree add /tmp/swap-legacy legacy/swap-4.2.0
```

## Usage

Regression against these reference binaries is wired into pixi:

```bash
pixi run swap-ref       # Linux: run a case with tests/reference/swap420
pixi run swap-win-ref   # Windows: run a case with tests/reference/swap.exe
```

Both tasks invoke `tests/swap-cases/run_case.sh` (or `.cmd`) with the reference binary as `--exec`. Case selection and other flags follow the normal `run_case.sh` conventions.

## Portability notes

`swap420` is dynamically linked and built against an older glibc. On modern Linux distributions it should still run but requires the usual shared libraries (`libm`, `libpthread`, `libc`). If it fails to start, inspect `ldd tests/reference/swap420` to see which libraries are unresolved. Static relinking is out of scope here; rebuild from `legacy/swap-4.2.0` if a different target is needed.
