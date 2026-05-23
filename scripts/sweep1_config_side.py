#!/usr/bin/env python3
"""Sweep 1: retire truly-global config-side symbols via state%cfg.

For every line in config_to_variables.f90 of the form `<sym> = <config_path>`:
  1. Skip if config_path is not a pure `config%...` rvalue (transforms/derived)
  2. Skip if symbol has special handling (numlay, rsoil, etc.)
  3. For every consumer file that does `use variables, only: <sym>`:
     - drop sym from the use-clause
     - replace bare reads of sym with state%cfg%<path>
  4. Delete the adapter line in config_to_variables.f90
  5. Tombstone the decl in variables.f90
  6. Delete the zero-init in initialize.f90

The script preserves CRLF/LF line endings.

Usage:
    python3 scripts/sweep1_config_side.py [--dry-run]
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
ADAPTER = REPO / "src/io/toml/config_to_variables.f90"
VARIABLES = REPO / "src/core/variables.f90"
INITIALIZE = REPO / "src/core/initialize.f90"

# Symbols to skip — non-trivial copies or already handled
SKIP_SYMS = {
    "numlay",        # derived: numlay = config%soil%isoillay(size(...))
    "rsoil",         # non-config: assigned from a literal in some paths
    "swsrf", "swsec", # pilot — already done
    "swdra", "swsophy", "pondmx", "rsro", "rsroexp",   # already done
    "tepRrain", "tepRsnow",  # already done
}


def read_lf(path: Path):
    """Return (text_lf, eol) where eol is '\\r\\n' or '\\n'."""
    data = path.read_bytes().decode("utf-8")
    if "\r\n" in data:
        return data.replace("\r\n", "\n"), "\r\n"
    return data, "\n"


def write_lf(path: Path, text_lf: str, eol: str):
    out = text_lf.replace("\n", eol) if eol == "\r\n" else text_lf
    path.write_bytes(out.encode("utf-8"))


def find_pure_copy_assignments():
    """Scan adapter for `X = config%path` lines. Return list of (sym, config_path, line_num).
    Pure copies only — RHS is exactly `config%...` with no expressions/calls.
    """
    text, _ = read_lf(ADAPTER)
    results = []
    pat = re.compile(
        r"^(\s*)(\w+)\s*=\s*(config%[\w%]+)\s*(?:!.*)?$",
        flags=re.MULTILINE,
    )
    for m in pat.finditer(text):
        sym = m.group(2)
        path = m.group(3)
        if sym.lower() in {s.lower() for s in SKIP_SYMS}:
            continue
        # Ensure path is config%... terminating with a field name (no parens/calls)
        if "(" in path or ")" in path:
            continue
        # Skip if path equals exactly `config%X%sym` (sym at end) — we still proceed,
        # but make sure path is well-formed
        results.append((sym, path))
    return results


def find_consumer_files(sym: str):
    """Files that `use variables` and reference sym."""
    pat = rf"\b{re.escape(sym)}\b"
    cmd = ["grep", "-rlEi", pat, "--include=*.f90", "src/"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO, check=False)
    except Exception:
        return []
    files = []
    for fpath in r.stdout.splitlines():
        if "swap_legacy" in fpath: continue
        if "/state/" in fpath: continue
        if "/config/" in fpath: continue
        if "variables.f90" in fpath: continue
        if "initialize.f90" in fpath: continue
        if "config_to_variables" in fpath: continue
        files.append(fpath)
    return files


def find_test_files(sym: str):
    pat = rf"\b{re.escape(sym)}\b"
    cmd = ["grep", "-rlEi", pat, "--include=*.pf", "--include=*.f90", "tests/"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO, check=False)
    except Exception:
        return []
    return [f for f in r.stdout.splitlines() if "swap_legacy" not in f]


def strip_from_use_clause(text_lf: str, sym: str) -> str:
    """Within `use variables, only: ...` continuations, drop sym.

    Walks lines tracking whether we are inside a `use variables` clause.
    Within the clause, removes occurrences of sym in all positions:
      - leading on a continuation line: `         sym, &`
      - leading on a continuation line ending the list: `         sym`
      - middle: `, sym,` / `, sym &`
      - after `only:`: `only: sym, ...` / `only: sym ...`
    """
    lines = text_lf.split("\n")
    out = []
    in_use = False
    sp = re.escape(sym)
    for line in lines:
        code, sep, comment = line.partition("!")
        stripped = code.lstrip()
        if re.match(r"use\s+variables\b", stripped, re.IGNORECASE):
            in_use = True
        if in_use:
            new = code
            # 1) leading-on-continuation: "<spaces>sym," or "<spaces>sym &" or "<spaces>sym\n"
            new = re.sub(rf"^(\s+){sp}\s*,\s*", r"\1", new, flags=re.IGNORECASE)
            new = re.sub(rf"^(\s+){sp}\s*(?=\s*(?:&|$))", r"\1", new, flags=re.IGNORECASE)
            # 2) middle
            new = re.sub(rf",\s*{sp}\s*,", ",", new, flags=re.IGNORECASE)
            new = re.sub(rf",\s*{sp}\s*(?=\s*(?:&|$))", "", new, flags=re.IGNORECASE)
            # 3) after only:
            new = re.sub(rf"(only\s*:\s*){sp}\s*,\s*", r"\1", new, flags=re.IGNORECASE)
            new = re.sub(rf"(only\s*:\s*){sp}\s*(?=\s*(?:&|$))", r"\1", new, flags=re.IGNORECASE)
            # 4) right after a continuation marker
            new = re.sub(rf"(&\s+){sp}\s*,\s*", r"\1", new, flags=re.IGNORECASE)
            new = re.sub(rf"(&\s+){sp}\s*(?=\s*(?:&|$))", r"\1", new, flags=re.IGNORECASE)
            code = new
            line = code + sep + comment
        # Track end of use-clause: only count code-bearing lines
        code_rstrip = code.rstrip()
        if in_use and code_rstrip != "" and not code_rstrip.endswith("&"):
            in_use = False
        out.append(line)
    return "\n".join(out)


def _clause_contains_sym(lines, start_idx, word_pat):
    """Starting at start_idx (a `use variables` line), collect the multi-line
    clause and return (uses_sym, last_idx)."""
    n = len(lines)
    use_text = [lines[start_idx]]
    code = lines[start_idx].split("!", 1)[0].rstrip()
    k = start_idx
    while code.endswith("&") and k + 1 < n:
        k += 1
        use_text.append(lines[k])
        next_code = lines[k].split("!", 1)[0].rstrip()
        if next_code == "":
            code = "&"  # comment-only; treat as continuation
        else:
            code = next_code
    joined = " ".join(use_text)
    return bool(word_pat.search(joined)), k


def find_importing_subs(text_lf: str, sym: str):
    """Return list of (start_line, end_line) for subroutines that import sym
    via a `use variables, only: ... sym ...` clause — either at module scope
    (applies to all subs in the module) or in the sub's body."""
    lines = text_lf.split("\n")
    n = len(lines)
    out = []
    sp = re.escape(sym)
    word_pat = re.compile(rf"\b{sp}\b", re.IGNORECASE)

    # Build a list of subroutine boundaries (assume flat — no nesting inside a sub)
    sub_starts = []
    for idx, line in enumerate(lines):
        if re.match(r"\s*(pure\s+)?(elemental\s+)?(recursive\s+)?(subroutine|function)\s+\w+\s*\(",
                    line, re.IGNORECASE):
            sub_starts.append(idx)

    # Check for MODULE-level use clauses (between `module X` and first `contains`/sub).
    # If found and they import sym, treat all subs as importing.
    module_imports_sym = False
    first_sub = sub_starts[0] if sub_starts else n
    for idx in range(first_sub):
        stripped = lines[idx].lstrip()
        if re.match(r"use\s+variables\s*(?:$|!)", stripped, re.IGNORECASE):
            # Blanket use Variables (module scope) — imports everything
            module_imports_sym = True
            break
        if re.match(r"use\s+variables\b", stripped, re.IGNORECASE):
            uses, _ = _clause_contains_sym(lines, idx, word_pat)
            if uses:
                module_imports_sym = True
                break

    # Walk through subroutines.
    i = 0
    while i < n:
        # Match start of subroutine or function definition (with `(`)
        m = re.match(r"\s*(pure\s+)?(elemental\s+)?(recursive\s+)?(subroutine|function)\s+(\w+)\s*\(",
                     lines[i], re.IGNORECASE)
        if m:
            sub_start = i
            sub_name = m.group(5)
            # End is right before the next subroutine or end-of-file
            next_starts = [s for s in sub_starts if s > i]
            sub_end = (next_starts[0] - 1) if next_starts else n - 1

            # Module-level import wins: mark every sub as importing
            if module_imports_sym:
                out.append((sub_start, sub_end))
                i = sub_end + 1
                continue

            # Within [sub_start, sub_end], look for use variables clauses that contain sym
            k = sub_start
            in_use = False
            use_text = []
            uses_sym = False
            while k <= sub_end:
                line = lines[k]
                stripped = line.lstrip()
                # Blanket `use variables` (no only:) imports everything
                m_blanket = re.match(r"use\s+variables\s*(?:$|!)", stripped, re.IGNORECASE)
                if m_blanket:
                    uses_sym = True
                    break
                if re.match(r"use\s+variables\b", stripped, re.IGNORECASE):
                    in_use = True
                if in_use:
                    use_text.append(line)
                    code = line.split("!", 1)[0].rstrip()
                    if code == "":
                        k += 1
                        continue
                    if not code.endswith("&"):
                        joined = " ".join(use_text)
                        if word_pat.search(joined):
                            uses_sym = True
                            break
                        use_text = []
                        in_use = False
                k += 1
            if uses_sym:
                out.append((sub_start, sub_end))
            i = sub_end + 1
        else:
            i += 1
    return out


def replace_bare_reads_in_ranges(text_lf: str, sym: str, replacement: str, ranges) -> str:
    """Replace bare `sym` with `replacement` within the given subroutine ranges.
    Skips use-clauses and declarations.
    """
    sp = re.escape(sym)
    pat = re.compile(rf"(?<![%a-zA-Z0-9_]){sp}(?![a-zA-Z0-9_%])", flags=re.IGNORECASE)
    lines = text_lf.split("\n")
    if not ranges:
        return text_lf
    # Apply replacement within each range. Track use-clause spans (multi-line)
    # and skip them entirely. A use-clause ends at the first NON-COMMENT, NON-EMPTY
    # code-line whose code part does NOT end with `&`.
    for (start, end) in ranges:
        in_use_block = False
        for k in range(start, end + 1):
            line = lines[k]
            code, sep, comment = line.partition("!")
            stripped = code.lstrip()
            # Detect use-clause start
            if re.match(r"use\s+", stripped, re.IGNORECASE):
                in_use_block = True
                continue
            if in_use_block:
                code_rstrip = code.rstrip()
                # Comment-only or blank line within the use clause — stay in block
                if code_rstrip == "":
                    continue
                # Code present; check if continuation
                if not code_rstrip.endswith("&"):
                    in_use_block = False
                continue
            # Declaration lines (real(8) :: x, intent(in) :: y) — skip
            if "::" in code:
                continue
            lines[k] = pat.sub(replacement, code) + sep + comment
    return "\n".join(lines)


def is_use_clause_line(line: str) -> bool:
    return bool(re.match(r"\s*use\s+variables\b", line, re.IGNORECASE))


def fix_use_clause_artifacts(text_lf: str) -> str:
    """Heuristic fixups after use-clause edits:
    - "only: ," → "only:"
    - " ,  ," → " , "
    - Lines whose entire code is just `&` (orphan continuations) → delete
    - Empty `use variables, only:` clauses → delete the whole clause
    """
    text_lf = re.sub(r"only\s*:\s*,\s*", "only: ", text_lf, flags=re.IGNORECASE)
    text_lf = re.sub(r",\s*,", ",", text_lf)
    # Drop orphan `&` lines
    new_lines = []
    for line in text_lf.split("\n"):
        code, sep, comment = line.partition("!")
        if code.strip() == "&":
            continue
        new_lines.append(line)
    text_lf = "\n".join(new_lines)

    # Drop empty `use variables, only:` clauses.
    # The clause might span multiple lines (continuations); detect by scanning
    # forward until non-continuation, non-comment, non-empty line. If the entire
    # clause body has no symbols, delete all lines making up the clause.
    lines = text_lf.split("\n")
    out = []
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.lstrip()
        # Only handle `use variables, only:` form — never delete blanket `use variables`.
        m = re.match(r"use\s+variables\b\s*,\s*only\s*:", stripped, re.IGNORECASE)
        if m:
            # Build clause by collecting lines while previous code link is open.
            clause_idxs = [i]
            j = i
            code_at_j = lines[i].split("!", 1)[0].rstrip()
            chain_open = code_at_j.endswith("&")
            while chain_open and j + 1 < len(lines):
                code_next = lines[j + 1].split("!", 1)[0].rstrip()
                if code_next == "":
                    # Comment-only line — part of clause (& chain not broken)
                    j += 1
                    clause_idxs.append(j)
                    continue
                # Real code line. If it starts a NEW statement (use/implicit/
                # type-declaration), the clause was actually broken upstream —
                # don't include this line in the clause.
                stripped_next = code_next.lstrip()
                if re.match(
                    r"(use|implicit|integer|real|logical|character|type|select|"
                    r"if|else|endif|end\s|do|module|subroutine|function|"
                    r"return|call|associate)\b",
                    stripped_next,
                    re.IGNORECASE,
                ):
                    break
                # Real code line: include only if chain still open
                if chain_open:
                    j += 1
                    clause_idxs.append(j)
                    chain_open = code_next.endswith("&")
                else:
                    break
            joined = " ".join(lines[k].split("!", 1)[0] for k in clause_idxs)
            content = re.sub(r"use\s+variables\b\s*,\s*only\s*:", "", joined, flags=re.IGNORECASE)
            content = content.replace("&", "").replace(",", "").strip()
            if not content:
                # Empty clause — delete all clause lines
                i = j + 1
                continue
        out.append(line)
        i += 1
    return "\n".join(out)


def migrate_symbol(sym: str, config_path: str, dry_run: bool = False):
    """Run the full migration for one symbol. Returns (touched_files, ok)."""
    state_path = config_path.replace("config%", "state%cfg%")
    print(f"  → {sym}: {config_path} → {state_path}")
    consumers = find_consumer_files(sym)
    tests = find_test_files(sym)
    touched = []

    for fpath in consumers:
        path = REPO / fpath
        if not path.exists():
            continue
        text_lf, eol = read_lf(path)
        original = text_lf
        # IMPORTANT ORDER: detect importing subs on the original text first,
        # then strip from use clause, then replace reads in those sub ranges.
        ranges = find_importing_subs(text_lf, sym)
        text_lf = strip_from_use_clause(text_lf, sym)
        text_lf = replace_bare_reads_in_ranges(text_lf, sym, state_path, ranges)
        text_lf = fix_use_clause_artifacts(text_lf)
        if text_lf != original:
            if not dry_run:
                write_lf(path, text_lf, eol)
            touched.append(fpath)

    # Test files: just patch use-clause + replace reads (less strict — tests
    # might also use the legacy global as the assertion target).
    for fpath in tests:
        path = REPO / fpath
        if not path.exists():
            continue
        text_lf, eol = read_lf(path)
        original = text_lf
        ranges = find_importing_subs(text_lf, sym)
        text_lf = strip_from_use_clause(text_lf, sym)
        text_lf = replace_bare_reads_in_ranges(text_lf, sym, config_path, ranges)
        text_lf = fix_use_clause_artifacts(text_lf)
        if text_lf != original:
            if not dry_run:
                write_lf(path, text_lf, eol)
            touched.append(fpath)

    # Delete adapter line
    text_lf, eol = read_lf(ADAPTER)
    pat = re.compile(
        rf"^(\s*){re.escape(sym)}\s*=\s*{re.escape(config_path)}\s*(?:!.*)?$\n?",
        flags=re.MULTILINE | re.IGNORECASE,
    )
    new = pat.sub("", text_lf)
    # Also replace bare reads of `sym` in the adapter's own conditional logic
    # with `config_path`. The adapter has config in scope.
    sp = re.escape(sym)
    bare_pat = re.compile(rf"(?<![%a-zA-Z0-9_]){sp}(?![a-zA-Z0-9_%])", flags=re.IGNORECASE)
    new_lines = []
    for line in new.split("\n"):
        code, sep, comment = line.partition("!")
        if "::" not in code:
            code = bare_pat.sub(config_path, code)
        new_lines.append(code + sep + comment)
    new = "\n".join(new_lines)
    if new != text_lf:
        if not dry_run:
            write_lf(ADAPTER, new, eol)
        touched.append("src/io/toml/config_to_variables.f90")

    # Tombstone decl in variables.f90
    text_lf, eol = read_lf(VARIABLES)
    # Match: `      real(8)   sym ! comment` OR `      integer sym(array) ! ...`
    pat = re.compile(
        rf"^(\s+)(real\(8\)|integer|logical|character\([^)]+\))(\s+){re.escape(sym)}(\([^)]*\))?(\s*!.*)?$",
        flags=re.MULTILINE | re.IGNORECASE,
    )
    new = pat.sub(rf"\1! [GR-CROP-DVS] {sym} retired — see {state_path}", text_lf)
    if new != text_lf:
        if not dry_run:
            write_lf(VARIABLES, new, eol)
        touched.append("src/core/variables.f90")
    # Also try: `real(8) sym, sym2, ...` — multi-symbol decl
    # Just handle removing one entry from a list; harder edge case, skip for now.

    # Delete zero-init in initialize.f90
    text_lf, eol = read_lf(INITIALIZE)
    pat = re.compile(
        rf"^\s+{re.escape(sym)}\s*=\s*[^\n]*\n",
        flags=re.MULTILINE | re.IGNORECASE,
    )
    new = pat.sub("", text_lf)
    if new != text_lf:
        if not dry_run:
            write_lf(INITIALIZE, new, eol)
        touched.append("src/core/initialize.f90")

    return touched, True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--limit", type=int, default=0, help="Only process first N symbols")
    ap.add_argument("--only", action="append", default=[], help="Only process named symbol")
    args = ap.parse_args()

    pairs = find_pure_copy_assignments()
    if args.only:
        only = {s.lower() for s in args.only}
        pairs = [p for p in pairs if p[0].lower() in only]
    if args.limit:
        pairs = pairs[: args.limit]

    print(f"Found {len(pairs)} pure-copy adapter lines to migrate")
    print()

    for sym, path in pairs:
        touched, ok = migrate_symbol(sym, path, dry_run=args.dry_run)
        print(f"     {len(touched)} files touched")
    print()
    print("Done. Build with: pixi run -e test build-linux")


if __name__ == "__main__":
    main()
