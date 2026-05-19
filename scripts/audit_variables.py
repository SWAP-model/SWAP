#!/usr/bin/env python3
"""Audit src/core/variables.f90 active declarations.

For each declared symbol, classify as:
- state-homed: a state type has a field of this name (case-insensitive)
- config-side-adapter: assigned from `config%X%Y` in config_to_variables.f90
- config-side-init: assigned from `cfg%X` in any crop*init.f90 file
- runtime-no-home: written by runtime code but no state field exists yet
- unknown: not assigned anywhere reachable from active code

Outputs a markdown table to stdout.
"""

import re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

VARIABLES = REPO / "src/core/variables.f90"
STATE_DIR = REPO / "src/state"
CONFIG_ADAPTER = REPO / "src/io/toml/config_to_variables.f90"
INIT_FILES = [
    REPO / "src/crop/cropfixed_init.f90",
    REPO / "src/crop/cropwofost_init.f90",
    REPO / "src/crop/cropgrass_init.f90",
]


def list_active_symbols():
    """Return list of (type, symbol) for active (non-tombstoned) decls."""
    out = []
    for line in VARIABLES.read_text().splitlines():
        stripped = line.lstrip()
        if stripped.startswith("!"):
            continue
        # Match: real(8) sym, integer sym, logical sym, character...
        m = re.match(
            r"\s+(real\(\d+\)|real\(real64\)|integer(?:\(\d+\))?|logical|character\([^)]+\)|character\*\d+)\s+(\w+)",
            line,
        )
        if m:
            t, sym = m.group(1), m.group(2)
            out.append((t, sym))
    # Dedupe while preserving first-seen
    seen = set()
    uniq = []
    for t, s in out:
        if s.lower() in seen:
            continue
        seen.add(s.lower())
        uniq.append((t, s))
    return uniq


def state_homes():
    """Map lowercase symbol -> state filename basename."""
    homes = {}
    for f in STATE_DIR.glob("*.f90"):
        text = f.read_text()
        for m in re.finditer(
            r"::\s+(\w+)\s*(?:\([^)]+\))?\s*=", text
        ):
            sym = m.group(1)
            homes.setdefault(sym.lower(), f.name)
    return homes


def adapter_writes():
    """Return set of lowercase symbols written from config%X in the adapter."""
    text = CONFIG_ADAPTER.read_text()
    syms = set()
    for m in re.finditer(
        r"^\s*(\w+)(?:\([^)]*\))?\s*=\s*config%", text, flags=re.MULTILINE
    ):
        syms.add(m.group(1).lower())
    # Also catch `call copy_table(config%X, sym)` patterns
    for m in re.finditer(
        r"call\s+copy_table\s*\(\s*config%[^,]+,\s*(\w+)", text, flags=re.IGNORECASE
    ):
        syms.add(m.group(1).lower())
    return syms


def init_writes():
    """Return set of lowercase symbols written from cfg% in any init file."""
    syms = set()
    for f in INIT_FILES:
        if not f.exists():
            continue
        text = f.read_text()
        for m in re.finditer(
            r"^\s*(\w+)(?:\([^)]*\))?\s*=\s*cfg%", text, flags=re.MULTILINE
        ):
            syms.add(m.group(1).lower())
        for m in re.finditer(
            r"call\s+copy_table\s*\(\s*cfg%[^,]+,\s*(\w+)", text, flags=re.IGNORECASE
        ):
            syms.add(m.group(1).lower())
    return syms


def find_readers(sym: str):
    """grep for any non-comment reference to <sym> outside skip patterns.
    Returns list of (file, line) tuples.
    """
    pat = rf"\b{re.escape(sym)}\b"
    cmd = ["grep", "-rnEi", pat, "--include=*.f90", "src/"]
    refs = []
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO, check=False)
    except Exception:
        return refs
    for line in r.stdout.splitlines():
        # filename:line:content
        parts = line.split(":", 2)
        if len(parts) < 3:
            continue
        fpath, lineno, content = parts[0], parts[1], parts[2]
        if "swap_legacy" in fpath:
            continue
        if "/state/" in fpath:
            continue
        if "/config/" in fpath:
            continue
        if "variables.f90" in fpath:
            continue
        # Skip pure comment lines
        code = content.split("!", 1)[0]
        if not code.strip():
            continue
        # Skip the declaration in init pattern (cfg%X assignments are config-side,
        # not real consumers)
        refs.append((fpath, lineno))
    return refs


def reader_subroutines(sym: str):
    """Approximate: count distinct subroutines that reference the symbol by
    scanning back from each ref to the nearest `subroutine NAME` line."""
    pat = rf"\b{re.escape(sym)}\b"
    cmd = ["grep", "-rnEi", pat, "--include=*.f90", "src/"]
    subs = set()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO, check=False)
    except Exception:
        return subs
    file_lines = {}
    for line in r.stdout.splitlines():
        parts = line.split(":", 2)
        if len(parts) < 3:
            continue
        fpath, lineno, content = parts[0], int(parts[1]), parts[2]
        if "swap_legacy" in fpath or "/state/" in fpath or "/config/" in fpath:
            continue
        if "variables.f90" in fpath or "initialize.f90" in fpath:
            continue
        # comment-only skip
        code = content.split("!", 1)[0]
        if not code.strip():
            continue
        if fpath not in file_lines:
            try:
                file_lines[fpath] = (REPO / fpath).read_text().splitlines()
            except Exception:
                continue
        lines = file_lines[fpath]
        # Walk back to find enclosing subroutine
        for j in range(min(lineno - 1, len(lines) - 1), -1, -1):
            m = re.match(r"\s*(subroutine|function)\s+(\w+)", lines[j], re.IGNORECASE)
            if m:
                subs.add(f"{fpath}:{m.group(2).lower()}")
                break
            # Stop at end of previous subroutine (we'd go too far)
            if re.match(r"\s*end\s+(subroutine|function)", lines[j], re.IGNORECASE) and j < lineno - 1:
                break
    return subs


def find_writers(sym: str):
    """grep for any write to <sym>; return list of (file, kind).
    Kinds: assign, read-stmt, copy-table-dst, output-arg.
    """
    writers = []
    pat = rf"\b{re.escape(sym)}\b\s*(\([^)]*\))?\s*="
    cmd = ["grep", "-rnE", "-i", pat, "--include=*.f90", "src/"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO, check=False)
    except Exception:
        return writers
    for line in r.stdout.splitlines():
        # Skip declaration sites and config field declarations
        if "::" in line.split("!")[0]:
            continue
        # Skip comparisons (== or .eq.)
        before_eq = line.split("=", 1)[0]
        after_eq = line.split("=", 1)[1] if "=" in line else ""
        if after_eq.startswith("="):  # ==
            continue
        # Skip if line starts with ! comment
        code = line.split(":", 2)[-1].lstrip()
        if code.startswith("!"):
            continue
        writers.append(line)
    # Also check copy_table destinations
    pat2 = rf"call\s+copy_table\s*\([^,)]+,\s*{re.escape(sym)}\b"
    cmd2 = ["grep", "-rnEi", pat2, "--include=*.f90", "src/"]
    try:
        r2 = subprocess.run(cmd2, capture_output=True, text=True, cwd=REPO, check=False)
        for line in r2.stdout.splitlines():
            writers.append(line + "  [copy_table-dst]")
    except Exception:
        pass
    return writers


def written_in_runtime(sym: str):
    """True if there's a non-adapter, non-init, non-state writer."""
    for w in find_writers(sym):
        if "swap_legacy" in w:
            continue
        if "/state/" in w:
            continue
        if "/config/" in w:  # config types
            continue
        if "initialize.f90" in w:
            continue
        if "variables.f90" in w:
            continue
        if "config_to_variables" in w:
            continue
        if any(init.name in w for init in INIT_FILES):
            continue
        return True
    return False


def written_by_adapter_or_init(sym: str):
    for w in find_writers(sym):
        if "config_to_variables" in w:
            return True
        if any(init.name in w for init in INIT_FILES):
            return True
    return False


def any_writer(sym: str):
    return len(find_writers(sym)) > 0


def has_runtime_reader(sym: str):
    """True if symbol is read (referenced) outside variables.f90 and config files."""
    pat = rf"\b{re.escape(sym)}\b"
    cmd = ["grep", "-rlEi", pat, "--include=*.f90", "src/"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO, check=False)
    except Exception:
        return False
    for line in r.stdout.splitlines():
        if "swap_legacy" in line:
            continue
        if "/state/" in line:
            continue
        if "/config/" in line:
            continue
        if "variables.f90" in line:
            continue
        if "initialize.f90" in line:
            continue
        if "config_to_variables" in line:
            continue
        if any(init.name in line for init in INIT_FILES):
            continue
        return True
    return False


def classify(sym: str, homes: dict, adapter: set, init: set):
    lc = sym.lower()
    state_home = homes.get(lc)
    has_runtime_writer = written_in_runtime(sym)
    config_writer = written_by_adapter_or_init(sym)
    has_any_writer = any_writer(sym)

    # state-homed takes precedence
    if state_home and has_runtime_writer:
        return "state-homed-runtime"
    if state_home and config_writer:
        return "state-homed-config-only"
    if state_home and not has_any_writer:
        return "state-homed-dead"
    if state_home:
        return "state-homed-unknown"
    # no state home
    if config_writer and not has_runtime_writer:
        return "config-side"
    if has_runtime_writer:
        return "runtime-no-home"
    # No writer found anywhere active. Two sub-cases:
    if has_runtime_reader(sym):
        return "orphan-read-only"  # read but never written — potential bug or config-gap
    return "dead-no-readers"


def reader_files(sym: str):
    """Return set of files where sym appears in non-comment code (excluding
    skip patterns: legacy, state, config, variables.f90, initialize.f90)."""
    pat = rf"\b{re.escape(sym)}\b"
    cmd = ["grep", "-rlEi", pat, "--include=*.f90", "src/"]
    files = set()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO, check=False)
    except Exception:
        return files
    for fpath in r.stdout.splitlines():
        if "swap_legacy" in fpath: continue
        if "/state/" in fpath: continue
        if "/config/" in fpath: continue
        if "variables.f90" in fpath: continue
        if "initialize.f90" in fpath: continue
        if "config_to_variables" in fpath: continue
        files.add(fpath)
    return files


def writer_files(sym: str):
    files = set()
    for w in find_writers(sym):
        parts = w.split(":", 1)
        if not parts: continue
        fpath = parts[0]
        if "swap_legacy" in fpath: continue
        if "/state/" in fpath: continue
        if "/config/" in fpath: continue
        if "variables.f90" in fpath: continue
        if "initialize.f90" in fpath: continue
        if "config_to_variables" in fpath: continue
        files.add(fpath)
    return files


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--locality", action="store_true",
                    help="Include reader-file-count, writer-file-count, sub-count columns")
    args = ap.parse_args()

    actives = list_active_symbols()
    homes = state_homes()
    adapter = adapter_writes()
    init = init_writes()

    rows = []
    for t, sym in actives:
        cat = classify(sym, homes, adapter, init)
        home = homes.get(sym.lower(), "")
        extra = {}
        if args.locality:
            rfiles = reader_files(sym)
            wfiles = writer_files(sym)
            subs = reader_subroutines(sym)
            extra["nread_files"] = len(rfiles)
            extra["nwrite_files"] = len(wfiles)
            extra["nsubs"] = len(subs)
            extra["rfiles"] = sorted(rfiles)
            extra["subs"] = sorted(subs)
        rows.append((sym, t, cat, home, extra))

    # Print summary
    from collections import Counter
    counts = Counter(r[2] for r in rows)
    print("# variables.f90 audit")
    print()
    print(f"Total active declarations: {len(rows)}")
    print()
    for cat, n in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"- **{cat}**: {n}")
    print()

    if args.locality:
        # Single-sub candidate count
        single_sub = [r for r in rows if r[4].get("nsubs", 999) == 1
                       and r[4].get("nread_files", 999) == 1
                       and r[4].get("nwrite_files", 999) <= 1]
        print(f"## Locality candidates")
        print(f"Symbols touched by exactly 1 subroutine in 1 file: {len(single_sub)}")
        print()
        print("| Symbol | Type | Category | Sub | File |")
        print("|---|---|---|---|---|")
        for sym, t, cat, home, extra in sorted(single_sub, key=lambda r: r[0].lower()):
            subs = list(extra.get("subs", []))
            files = list(extra.get("rfiles", []))
            print(f"| `{sym}` | {t} | {cat} | {subs[0] if subs else ''} | {files[0] if files else ''} |")
        print()

    print("## Per-symbol table")
    print()
    if args.locality:
        print("| Symbol | Type | Category | State home | NReadFiles | NWriteFiles | NSubs |")
        print("|---|---|---|---|---|---|---|")
        for sym, t, cat, home, extra in sorted(rows, key=lambda r: (r[2], r[0].lower())):
            print(f"| `{sym}` | {t} | {cat} | {home} | {extra.get('nread_files','')} | {extra.get('nwrite_files','')} | {extra.get('nsubs','')} |")
    else:
        print("| Symbol | Type | Category | State home |")
        print("|---|---|---|---|")
        for sym, t, cat, home, _ in sorted(rows, key=lambda r: (r[2], r[0].lower())):
            print(f"| `{sym}` | {t} | {cat} | {home} |")


if __name__ == "__main__":
    main()
