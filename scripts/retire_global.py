#!/usr/bin/env python3
"""Per-symbol globals migration helper for GR-DVS arc.

Usage:
    python scripts/retire_global.py <symbol> <state_path> [--has-t-suffix]

Steps:
1. sed across crop-runtime subroutine ranges (state in scope).
2. Fix use-clauses (replace 'state%...%<symbol>' import with just the comma).
3. Delete self-assignment lines (X = X).
4. Patch variables.f90 declaration.
5. Patch initialize.f90 zero-init.
6. Patch swap_mod.f90 dual-write (X = X assignment + use clause).

Does NOT run the build or commit — caller should verify and commit.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# Subroutine ranges where module-level globals appear and `state` is in scope.
# Outside these ranges, the same symbol name may be a local parameter — leave alone.
RANGES = {
    'src/crop/cropwofost_runtime.f90': (19, 1282),
    'src/crop/cropgrass_runtime.f90': (18, 1430),
    'src/crop/cropfixed_runtime.f90': (15, 274),
    'src/crop/cropgrowth.f90': (20, 610),
    'src/crop/cropgrowth_helpers.f90': (28, 380),
    'src/crop/rootextraction.f90': (31, 1008),
    'src/crop/oxygenstress.f90': (96, 660),
    'src/crop/irrigation.f90': (33, 488),
    'src/atmosphere/meteoday.f90': (1, 9999),
    'src/atmosphere/interception.f90': (1, 9999),
    'src/soil/waterbalance.f90': (1, 9999),
    'src/solute/agetracer.f90': (1, 9999),
    'src/solute/solute.f90': (1, 9999),
    'src/drainage/drainage.f90': (1, 9999),
    'src/drainage/surfacewater.f90': (1, 9999),
    'src/atmosphere/snow.f90': (1, 9999),
    'src/heat/temperature.f90': (1, 9999),
    'src/heat/frozencond.f90': (1, 9999),
}


def _read_lines(path: Path):
    """Read file preserving line endings. Returns (lines, eol) where eol is '\r\n' or '\n'."""
    data = path.read_bytes().decode('utf-8')
    if '\r\n' in data:
        return data.split('\r\n'), '\r\n'
    return data.split('\n'), '\n'


def _write_lines(path: Path, lines, eol):
    path.write_bytes(eol.join(lines).encode('utf-8'))


def sed_migrate(symbol: str, state_path: str, has_t_suffix: bool):
    excl = '[^a-zA-Z0-9_t]' if has_t_suffix else '[^a-zA-Z0-9_]'
    for rel, (lo, hi) in RANGES.items():
        path = REPO / rel
        if not path.exists():
            continue
        pattern = (
            f'{lo},{hi}{{ /^\\s*!/!s/([^%a-zA-Z0-9_]|^){symbol}'
            f'({excl}|$)/\\1{state_path}\\2/gI }}'
        )
        subprocess.run(['sed', '-i', '-E', pattern, str(path)], check=True)


def fix_use_clauses(symbol: str, state_path: str):
    """Remove the mangled 'state%...%<sym>' from use-variables import lists.

    Walk lines, tracking continuation state via trailing '&'. A use clause
    starts at `use variables` and continues until the first line that does
    NOT end with '&' (ignoring trailing comments).
    """
    sp_re = re.escape(state_path)

    def strip_state_path(line: str) -> str:
        # Strip the import token in any position
        new = line
        new = re.sub(r',\s*' + sp_re + r'\s*(?=,)', ',', new)
        new = re.sub(r',\s*' + sp_re + r'\s*(?=\s*&)', '', new)
        new = re.sub(r',\s*' + sp_re + r'\s*$', '', new)
        new = re.sub(sp_re + r'\s*,\s*', '', new)
        new = re.sub(r'\b' + sp_re + r'\b', '', new)
        new = re.sub(r',\s*,', ',', new)
        return new

    files_to_check = list(RANGES.keys()) + ['src/core/swap_mod.f90']
    for rel in files_to_check:
        path = REPO / rel
        if not path.exists():
            continue
        lines, eol = _read_lines(path)
        in_use = False
        out_lines = []
        for line in lines:
            stripped = line.lstrip()
            if stripped.startswith('use variables') or stripped.startswith('use Variables'):
                in_use = True
            if in_use and state_path in line:
                line = strip_state_path(line)
            code_part = line.split('!', 1)[0]
            code_stripped = code_part.rstrip()
            ends_amp = code_stripped.endswith('&')
            out_lines.append(line)
            if in_use and not ends_amp:
                in_use = False
        # NOTE: trailing-comma cleanup left to manual fix-up after build error.
        # Auto-fixing is brittle (comments, multi-line patterns).
        text = eol.join(out_lines)
        # Final pass: fix trailing commas right before a continuation that has no symbols
        # Pattern: ",\s*&\s*(![^\n]*)?\n\s*&" → "  &" then on a body line that ends with ", &" but the
        # following continuation is empty.
        # Simpler heuristic: drop bare "," that immediately precedes "&" at end of code on a line
        # whose own list has had all entries stripped.
        original = path.read_bytes().decode('utf-8')
        if text != original:
            path.write_bytes(text.encode('utf-8'))


def delete_self_assignments(state_path: str):
    """Remove lines where 'state%X = state%X' (no RHS computation)."""
    # Match: state%path = state%path  (optionally with trailing comment)
    pat = re.compile(
        r'^(\s*)' + re.escape(state_path) +
        r'\s*=\s*' + re.escape(state_path) +
        r'\s*(![^\n]*)?$'
    )
    for rel in list(RANGES.keys()) + [
        'src/core/swap_mod.f90',
        'src/crop/cropgrowth.f90',
    ]:
        path = REPO / rel
        if not path.exists():
            continue
        lines, eol = _read_lines(path)
        out_lines = []
        for line in lines:
            if pat.match(line.rstrip('\r')):
                continue  # drop
            out_lines.append(line)
        new_text = eol.join(out_lines)
        original = path.read_bytes().decode('utf-8')
        if new_text != original:
            path.write_bytes(new_text.encode('utf-8'))


def patch_variables_decl(symbol: str, state_path: str):
    """Replace real(8) <symbol> ... line with a tombstone comment."""
    path = REPO / 'src/core/variables.f90'
    text = path.read_bytes().decode('utf-8')
    new_text = re.sub(
        r'^(\s*)real\(8\)\s+' + re.escape(symbol) + r'\s+!\s*([^\n]*)$',
        rf'\1! [GR-CROP-DVS] {symbol} retired — see {state_path}',
        text, flags=re.MULTILINE
    )
    if new_text != text:
        path.write_bytes(new_text.encode('utf-8'))


def patch_initialize_zero(symbol: str):
    """Drop the '<symbol> = 0.0d0' line in initialize.f90."""
    path = REPO / 'src/core/initialize.f90'
    text = path.read_bytes().decode('utf-8')
    new_text = re.sub(
        r'^\s*' + re.escape(symbol) + r'\s*=\s*0\.0?d0\s*\n',
        '', text, flags=re.MULTILINE
    )
    if new_text != text:
        path.write_bytes(new_text.encode('utf-8'))


def patch_swap_mod_dual_write(symbol: str, state_path: str):
    """Drop 'state%X = symbol' dual-write in swap_mod and drop from its use clause."""
    path = REPO / 'src/core/swap_mod.f90'
    text = path.read_bytes().decode('utf-8')
    # Drop dual-write assignment
    text = re.sub(
        r'^\s*' + re.escape(state_path) + r'\s*=\s*' + re.escape(symbol) + r'\s*\n',
        '', text, flags=re.MULTILINE
    )
    # Drop from use clause (handle both ', sym,' and 'sym,')
    text = re.sub(r',\s*' + re.escape(symbol) + r'\b', '', text)
    text = re.sub(r'\b' + re.escape(symbol) + r'\s*,\s*', '', text)
    path.write_bytes(text.encode('utf-8'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('symbol')
    ap.add_argument('state_path')
    ap.add_argument('--has-t-suffix', action='store_true',
                    help='Add t to exclusion class (e.g. for symbols where suffix-t versions exist like tagp/tagpt)')
    args = ap.parse_args()

    print(f'Migrating {args.symbol} → {args.state_path}')
    sed_migrate(args.symbol, args.state_path, args.has_t_suffix)
    fix_use_clauses(args.symbol, args.state_path)
    delete_self_assignments(args.state_path)
    patch_variables_decl(args.symbol, args.state_path)
    patch_initialize_zero(args.symbol)
    patch_swap_mod_dual_write(args.symbol, args.state_path)
    print('Done. Verify with: pixi run -e test build-linux && pixi run -e test check-fast')


if __name__ == '__main__':
    main()
