#!/usr/bin/env python3
"""show_consumers.py — locate every reference to a symbol across src/ and tests/.

Categorizes each reference so a human can manually edit each site with
informed eyes. Does NOT modify any files.

Usage:
    python3 scripts/show_consumers.py <SYMBOL>

Output sections:
  1) declaration  — type declaration lines
  2) zero-init    — zero-initialisation in initialize.f90
  3) adapter      — legacy config_to_variables adapter writes (historical; bucket
                    will be empty now that config_to_variables.f90 is gone)
  4) init         — crop*init.f90 lines (X = cfg%X)
  5) use-clauses  — files that `use variables, only: ..., X, ...` (historical;
                    bucket will be empty now that variables.f90 is gone)
  6) signatures   — subroutine/function lines where sym appears as a dummy arg
  7) locals       — files where sym is declared as a local var or intent arg
  8) reads/writes — every other reference, grouped by file with line context
  9) tests        — references in tests/

Rule of thumb for manual editing:
  - Items in (6) and (7) are LOCAL to that subroutine — DO NOT replace.
  - Items in (8) become `state%<path>%X` (or relevant state path),
    EXCEPT inside subroutines that appear in (6)/(7) for that file.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

SKIP_DIR_PATTERNS = ["swap_legacy", "/state/", "/config/", "subprojects/"]


def is_skipped(fpath: str) -> bool:
    return any(p in fpath for p in SKIP_DIR_PATTERNS)


def _is_use_clause_continuation(fpath: str, lineno: int) -> bool:
    """True if this line is inside a `use variables` clause via `&` chain."""
    path = REPO / fpath
    if not path.exists():
        return False
    try:
        text_lines = path.read_text(errors="ignore").splitlines()
    except Exception:
        return False
    j = lineno - 1
    while j > 0:
        code = text_lines[j - 1].split("!", 1)[0].rstrip()
        if not code:
            j -= 1
            continue
        if re.match(r"\s*use\s+variables\b", text_lines[j - 1], re.IGNORECASE):
            return True
        # Found a non-empty code line that's NOT a use clause and doesn't continue
        if not code.endswith("&"):
            return False
        j -= 1
    return False


def grep_lines(sym: str):
    """Return list of (file, line_num, line_text) for case-insensitive matches."""
    pat = rf"\b{re.escape(sym)}\b"
    cmd = ["grep", "-rnEi", pat, "--include=*.f90", "--include=*.pf", "src/", "tests/"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO, check=False)
    except Exception:
        return []
    out = []
    for line in r.stdout.splitlines():
        parts = line.split(":", 2)
        if len(parts) < 3:
            continue
        fpath, lineno, text = parts[0], parts[1], parts[2]
        if is_skipped(fpath):
            continue
        try:
            lineno_i = int(lineno)
        except ValueError:
            continue
        out.append((fpath, lineno_i, text))
    return out


def classify(refs, sym: str):
    """Group references by category. Returns dict of category -> list of (file, line, text)."""
    buckets = {
        "decl": [],          # declaration in variables.f90
        "init": [],          # X = cfg%X in crop*init.f90
        "adapter": [],       # X = config%... in config_to_variables.f90
        "zero_init": [],     # X = 0.X in initialize.f90
        "use_clause": [],    # use variables, only: ... X ...
        "signature": [],     # subroutine NAME(... X ...)
        "local_decl": [],    # real(8) X / intent(...) :: X (local var)
        "test": [],          # tests/ — for awareness
        "comment": [],       # comment-only mention (informational)
        "read_write": [],    # everything else — these are the migration sites
    }

    sp_re = re.compile(rf"\b{re.escape(sym)}\b", re.IGNORECASE)

    for fpath, lineno, text in refs:
        # Comment-only check
        code, _, _ = text.partition("!")
        if not sp_re.search(code):
            buckets["comment"].append((fpath, lineno, text))
            continue

        # variables.f90 declaration
        if "variables.f90" in fpath:
            if re.match(r"\s*(real|integer|logical|character)", code, re.IGNORECASE):
                buckets["decl"].append((fpath, lineno, text))
            else:
                buckets["read_write"].append((fpath, lineno, text))
            continue

        # initialize.f90 zero-init
        if "initialize.f90" in fpath and re.match(rf"\s*{re.escape(sym)}\s*=\s*[0.dDeE+\-]+", code, re.IGNORECASE):
            buckets["zero_init"].append((fpath, lineno, text))
            continue

        # config_to_variables.f90 adapter writes
        if "config_to_variables" in fpath and re.match(rf"\s*{re.escape(sym)}\s*=\s*config%", code, re.IGNORECASE):
            buckets["adapter"].append((fpath, lineno, text))
            continue

        # crop*_init.f90 writes (X = cfg%X)
        if "_init.f90" in fpath and re.match(rf"\s*{re.escape(sym)}\s*=\s*cfg%", code, re.IGNORECASE):
            buckets["init"].append((fpath, lineno, text))
            continue

        # Use clauses (start)
        if re.search(r"use\s+variables\b", code, re.IGNORECASE):
            buckets["use_clause"].append((fpath, lineno, text))
            continue

        # Use clause continuation: walk back to see if previous code lines
        # chain back to a `use variables` line through unbroken `&` continuations.
        # If so, treat this as a use-clause line too.
        if _is_use_clause_continuation(fpath, lineno):
            buckets["use_clause"].append((fpath, lineno, text))
            continue

        # Subroutine / function signature line
        if re.match(r"\s*(pure\s+|elemental\s+|recursive\s+)*(subroutine|function)\s+\w+\s*\(", code, re.IGNORECASE):
            buckets["signature"].append((fpath, lineno, text))
            continue

        # Local declaration (intent / type-decl with `::`)
        if "::" in code and re.search(rf"::\s*[^!]*\b{re.escape(sym)}\b", code, re.IGNORECASE):
            buckets["local_decl"].append((fpath, lineno, text))
            continue

        # Tests
        if fpath.startswith("tests/"):
            buckets["test"].append((fpath, lineno, text))
            continue

        # Everything else: read/write at runtime
        buckets["read_write"].append((fpath, lineno, text))

    return buckets


def find_continuation_signature_lines(refs, sym: str):
    """Some signatures span continuations: 'subroutine X(a,b, &\n  &  sym, c)'.
    Detect refs whose previous code line ends with `&` AND is inside a sig.

    Returns set of (file, line) tuples to be re-classified as signature.
    """
    by_file = {}
    for f, ln, t in refs:
        by_file.setdefault(f, []).append(ln)
    augmented = set()
    for fpath, lines in by_file.items():
        path = REPO / fpath
        if not path.exists():
            continue
        try:
            text_lines = path.read_text(errors="ignore").splitlines()
        except Exception:
            continue
        for ln in lines:
            # Walk back from ln-1 checking each code line ends with `&`.
            # If we hit a `use ...` statement, stop — we're inside a use clause,
            # not a signature.
            j = ln - 1
            while j > 0:
                code = text_lines[j - 1].split("!", 1)[0].rstrip()
                if not code:
                    j -= 1
                    continue
                # Encountered a `use` keyword starting a statement — abort walk
                if re.match(r"\s*use\s+", text_lines[j - 1], re.IGNORECASE):
                    break
                # Encountered a subroutine/function signature — classify as sig
                if re.match(r"\s*(pure\s+|elemental\s+|recursive\s+)*(subroutine|function)\s+\w+\s*\(",
                            text_lines[j - 1], re.IGNORECASE):
                    augmented.add((fpath, ln))
                    break
                if code.endswith("&"):
                    j -= 1
                else:
                    break
    return augmented


def print_section(title, items, *, with_path=True, indent=2):
    if not items:
        return
    print(f"\n=== {title} ({len(items)}) ===")
    last_file = None
    for fpath, lineno, text in items:
        if with_path and fpath != last_file:
            print(f"  {fpath}")
            last_file = fpath
        prefix = " " * indent if with_path else ""
        text_trim = text[:200]
        print(f"  {prefix}{lineno}: {text_trim}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("symbol")
    args = ap.parse_args()
    sym = args.symbol

    refs = grep_lines(sym)
    if not refs:
        print(f"No references found for {sym}")
        return 0

    buckets = classify(refs, sym)

    # Promote continuation-signature lines from read_write to signature
    cont_sig = find_continuation_signature_lines(buckets["read_write"], sym)
    promoted = []
    keep = []
    for ref in buckets["read_write"]:
        if (ref[0], ref[1]) in cont_sig:
            promoted.append(ref)
        else:
            keep.append(ref)
    buckets["signature"].extend(promoted)
    buckets["read_write"] = keep

    print(f"=== Symbol: {sym} ===")
    print(f"Total references: {len(refs)}")

    print_section("1) DECLARATION (delete after all readers migrated)", buckets["decl"])
    print_section("2) ZERO-INIT (delete with decl)", buckets["zero_init"])
    print_section("3) ADAPTER write (rewrite or delete)", buckets["adapter"])
    print_section("4) INIT write (`X = cfg%X` → `state%path%X = cfg%X`)", buckets["init"])
    print_section("5) USE-CLAUSE (drop sym from each `use variables` list)", buckets["use_clause"])
    print_section("6) SIGNATURE (LOCAL ARG — do NOT replace)", buckets["signature"])
    print_section("7) LOCAL DECL (LOCAL VAR — do NOT replace)", buckets["local_decl"])
    print_section("8) READ/WRITE (migrate each — check enclosing sub for local-shadow!)", buckets["read_write"])
    print_section("9) TESTS (update assertions/imports)", buckets["test"])
    if buckets["comment"]:
        print(f"\n=== 10) COMMENT-ONLY ({len(buckets['comment'])}) — informational, may be stale ===")

    # Summary
    print(f"\n=== SUMMARY ===")
    print(f"  read/write sites to migrate: {len(buckets['read_write'])}")
    print(f"  use-clause imports to drop:  {len(buckets['use_clause'])}")
    print(f"  init writes to redirect:     {len(buckets['init'])}")
    if buckets["signature"] or buckets["local_decl"]:
        print(f"  ⚠  LOCAL SHADOWS:")
        print(f"     signatures: {len(buckets['signature'])}")
        print(f"     local decls: {len(buckets['local_decl'])}")
        print(f"     → DO NOT replace sym inside these subroutines.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
