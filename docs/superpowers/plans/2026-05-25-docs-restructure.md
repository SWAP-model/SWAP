# Docs site restructure — Phase A+B Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Spec:** `docs/superpowers/specs/2026-05-25-docs-restructure-design.md` (this file moves to `dev-docs/superpowers/plans/` during Task 2; orchestrator tracks).

**Goal:** Reorganise the FORD documentation site so the public surface holds only Tutorial / Reference / Developer / API content, with an Atlas-styled navbar and inline typeahead search. All dev-only material (ADRs, specs, plans, archive, status reports) moves to `dev-docs/`.

**Architecture:** Two interleaved tracks. (1) **Files:** `git mv` dev material to `dev-docs/`, reorganise user/contributor pages into `docs/reference/` and `docs/developer/`, add section landers, rewrite `docs/index.md`. (2) **Templates:** new `docs/templates/base.html` overriding FORD's Bootstrap navbar with the Atlas palette plus a lunr.js-powered typeahead search bound to FORD's existing `search_database.json`. Intermediate verification gate after Task 8 confirms the structural cleanup before template work begins.

**Tech Stack:** FORD 6.x (Jinja2 templates), Python 3.11, lunr.js (CDN, already used by FORD), pixi tasks (`docs-rebuild`, `docs-serve`), Bootstrap 5 (kept for FORD inner pages), Atlas palette (Spectral/Inter/JetBrains Mono fonts).

**Verification gates:**
- After **Task 8**: `pixi run -e docs docs-rebuild` exits 0 with zero `Warning:` lines about missing titles, skipped directories, or unknown metadata keys.
- After **Task 11**: smoke test of served site + search index integrity check + internal-link audit all pass.

**Commit convention:** development-only commits; co-author trailer required (matches repo convention in `CLAUDE.md`).

---

## File map

**Created:**
- `dev-docs/README.md` — index of relocated dev material
- `docs/tutorial/index.md` — stub
- `docs/reference/index.md` — section lander
- `docs/developer/index.md` — section lander
- `docs/templates/base.html` — custom Atlas navbar + typeahead

**Modified:**
- `docs.md` — remove `./docs/superpowers` from `exclude_dir`
- `docs/index.md` — full rewrite (drop Phase-4 framing)
- `docs/reference/meteorology.md` — add YAML frontmatter with `title:`
- `docs/reference/csv-companion-files.md` — add YAML frontmatter with `title:`

**Deleted:**
- `docs/templates/_index.html` — stale stock template backup
- `swap_debug.log` — untracked debug artefact at repo root

**Moved (via `git mv`):**
- `docs/adr/` → `dev-docs/adr/` (42 files as a tree)
- `docs/superpowers/` → `dev-docs/superpowers/` (specs + plans + notes; the spec and this plan move with it)
- `docs/archive/` → `dev-docs/archive/` (2026-phase-4)
- `docs/PHASE-4-MODERNIZATION-SUMMARY.md` → `dev-docs/phase-4-modernization-summary.md`
- `docs/coverage-baseline.md` → `dev-docs/coverage-baseline.md`
- `docs/configuration-schema.md` → `docs/reference/configuration-schema.md`
- `docs/toml-format-guide.md` → `docs/reference/toml-format-guide.md`
- `docs/meteorology.md` → `docs/reference/meteorology.md`
- `docs/csv-companion-files.md` → `docs/reference/csv-companion-files.md`
- `docs/architecture.md` → `docs/developer/architecture.md`
- `docs/contributing.md` → `docs/developer/contributing.md`
- `docs/code-style.md` → `docs/developer/code-style.md`
- `docs/build-and-test.md` → `docs/developer/build-and-test.md`
- `docs/state-management.md` → `docs/developer/state-management.md`
- `docs/error-handling.md` → `docs/developer/error-handling.md`
- `docs/logging.md` → `docs/developer/logging.md`
- `docs/validation.md` → `docs/developer/validation.md`
- `docs/branches.md` → `docs/developer/branches.md`
- `docs/dependency-management.md` → `docs/developer/dependency-management.md`

---

## Task 1: Create `dev-docs/` and move the ADR tree

**Files:**
- Create dir: `dev-docs/`
- Move: `docs/adr/` → `dev-docs/adr/`

- [ ] **Step 1: Create the dev-docs/ directory and move ADRs**

```bash
mkdir -p dev-docs
git mv docs/adr dev-docs/adr
```

- [ ] **Step 2: Verify the move**

```bash
ls dev-docs/adr/ | wc -l
ls docs/adr 2>&1
```

Expected:
- First command prints `42` (or whatever the current ADR count is — adjust if a new ADR was added).
- Second command: `ls: cannot access 'docs/adr': No such file or directory`.

- [ ] **Step 3: Commit**

```bash
git commit -m "$(cat <<'EOF'
docs(restructure): move ADRs from docs/ to dev-docs/ — phase A

ADRs are dev-only history; they should not appear on the public docs
site. Whole tree moves together so sibling cross-references stay valid.

Spec: docs/superpowers/specs/2026-05-25-docs-restructure-design.md

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Move remaining internal dev material

**Files:**
- Move: `docs/superpowers/` → `dev-docs/superpowers/`
- Move: `docs/archive/` → `dev-docs/archive/`
- Move: `docs/PHASE-4-MODERNIZATION-SUMMARY.md` → `dev-docs/phase-4-modernization-summary.md`
- Move: `docs/coverage-baseline.md` → `dev-docs/coverage-baseline.md`

> **Note:** `docs/superpowers/` includes both this plan file and its spec. After this task, this plan lives at `dev-docs/superpowers/plans/2026-05-25-docs-restructure.md`. The orchestrator (or next subagent invocation) must reference the new path.

- [ ] **Step 1: Move the directories and individual files**

```bash
git mv docs/superpowers dev-docs/superpowers
git mv docs/archive dev-docs/archive
git mv docs/PHASE-4-MODERNIZATION-SUMMARY.md dev-docs/phase-4-modernization-summary.md
git mv docs/coverage-baseline.md dev-docs/coverage-baseline.md
```

- [ ] **Step 2: Verify**

```bash
ls dev-docs/
ls docs/superpowers docs/archive docs/PHASE-4-MODERNIZATION-SUMMARY.md docs/coverage-baseline.md 2>&1
```

Expected:
- First command lists at least: `adr archive phase-4-modernization-summary.md coverage-baseline.md superpowers`.
- Second command: all four targets report "No such file or directory".

- [ ] **Step 3: Commit**

```bash
git commit -m "$(cat <<'EOF'
docs(restructure): move superpowers/archive/status reports to dev-docs/

Specs, plans, notes, the Phase-4 archive, the modernization summary, and
the coverage baseline are dev-only artifacts. Out of the public docs
tree but preserved in repo.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Move user-facing reference pages and fix frontmatters

**Files:**
- Create dir: `docs/reference/`
- Move: `docs/configuration-schema.md` → `docs/reference/configuration-schema.md`
- Move: `docs/toml-format-guide.md` → `docs/reference/toml-format-guide.md`
- Move: `docs/meteorology.md` → `docs/reference/meteorology.md` (frontmatter fix)
- Move: `docs/csv-companion-files.md` → `docs/reference/csv-companion-files.md` (frontmatter fix)

- [ ] **Step 1: Create the directory and move the four pages**

```bash
mkdir -p docs/reference
git mv docs/configuration-schema.md docs/reference/configuration-schema.md
git mv docs/toml-format-guide.md docs/reference/toml-format-guide.md
git mv docs/meteorology.md docs/reference/meteorology.md
git mv docs/csv-companion-files.md docs/reference/csv-companion-files.md
```

- [ ] **Step 2: Add YAML frontmatter to meteorology.md**

The current first line of `docs/reference/meteorology.md` is `# Meteorology`. Prepend YAML frontmatter so FORD accepts it. Use the Edit tool with:

- **old_string:**
  ```
  # Meteorology
  ```
- **new_string:**
  ```
  ---
  title: Meteorology
  ---

  # Meteorology
  ```

- [ ] **Step 3: Add YAML frontmatter to csv-companion-files.md**

The current first line of `docs/reference/csv-companion-files.md` is `# CSV companion files`. Use the Edit tool:

- **old_string:**
  ```
  # CSV companion files
  ```
- **new_string:**
  ```
  ---
  title: CSV companion files
  ---

  # CSV companion files
  ```

- [ ] **Step 4: Verify frontmatter is in place**

```bash
head -4 docs/reference/meteorology.md docs/reference/csv-companion-files.md
```

Expected: both files start with `---\ntitle: ...\n---`.

Also check all four reference pages have a title:

```bash
grep -L "^title:" docs/reference/*.md
```

Expected: empty output (every file has a title).

- [ ] **Step 5: Commit**

```bash
git add docs/reference/
git commit -m "$(cat <<'EOF'
docs(restructure): move user reference pages to docs/reference/ + fix frontmatters

Adds YAML title frontmatter to meteorology.md and csv-companion-files.md
which FORD requires (was emitting parse-error warnings each build).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Move contributor pages to `docs/developer/`

**Files:**
- Create dir: `docs/developer/`
- Move ten contributor-facing pages.

- [ ] **Step 1: Create the directory and move the pages**

```bash
mkdir -p docs/developer
git mv docs/architecture.md          docs/developer/architecture.md
git mv docs/contributing.md          docs/developer/contributing.md
git mv docs/code-style.md            docs/developer/code-style.md
git mv docs/build-and-test.md        docs/developer/build-and-test.md
git mv docs/state-management.md      docs/developer/state-management.md
git mv docs/error-handling.md        docs/developer/error-handling.md
git mv docs/logging.md               docs/developer/logging.md
git mv docs/validation.md            docs/developer/validation.md
git mv docs/branches.md              docs/developer/branches.md
git mv docs/dependency-management.md docs/developer/dependency-management.md
```

- [ ] **Step 2: Verify and check frontmatter titles**

```bash
ls docs/developer/ | wc -l
grep -L "^title:" docs/developer/*.md
```

Expected:
- First command prints `10`.
- Second command may print page filenames missing a `title:` line. If any are listed, fix them by prepending YAML frontmatter the same way as in Task 3 (Step 2/3), deriving the title from the existing top-level `# Heading`.

- [ ] **Step 3: Verify only intended files remain at docs/ root**

```bash
find docs -maxdepth 1 -type f -name "*.md"
```

Expected output (order may vary):
```
docs/index.md
```

If any other `.md` shows up at this depth, it was missed — move it now (probably to `docs/developer/`) before committing.

- [ ] **Step 4: Commit**

```bash
git add docs/developer/
git commit -m "$(cat <<'EOF'
docs(restructure): move contributor pages to docs/developer/

architecture, contributing, code-style, build-and-test, state-management,
error-handling, logging, validation, branches, dependency-management.
Wired into the public site under the new Developer tab.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Write `dev-docs/README.md`

**Files:**
- Create: `dev-docs/README.md`

- [ ] **Step 1: Create the file**

Write `dev-docs/README.md` with this exact content:

```markdown
# Dev documentation (internal)

This directory holds development-only material that used to live under `docs/`
but is no longer published on the FORD-generated public site. It's kept here
so contributors can still find it via repository search and editor file
navigation.

The public site source is at `docs/`. The FORD config (`docs.md`) does not
scan this directory.

## Contents

### Architecture decisions

- `adr/` — 42+ Architecture Decision Records. Each ADR captures a non-obvious
  choice made during the modernization arcs and the reasoning behind it. The
  ADR carries custom metadata keys (`status:`, `supersedes:`,
  `supersedes-portion-of:`) that FORD doesn't recognise — that's why these
  files don't belong on the public site.

### Plans, specs, and notes

- `superpowers/specs/` — design documents brainstormed before implementation.
- `superpowers/plans/` — implementation plans corresponding to each spec.
- `superpowers/notes/` — ad-hoc per-arc notes.

### Phase 4 archive

- `archive/2026-phase-4/` — frozen specs, plans, and per-reader audits from
  the Phase 4 / Phase 4f-extend modernization (completed 2026-05-05; tags
  `rescue/phase-4f-extend-complete` and `rescue/phase-4f-extend-followups`).

### Historical reports

- `phase-4-modernization-summary.md` — capstone summary of the Phase 4 arc.
- `coverage-baseline.md` — Phase 3 line/branch coverage snapshot. Coverage
  is tracked, not gated (see `adr/0006-coverage-tracked-not-gated.md`).
```

- [ ] **Step 2: Commit**

```bash
git add dev-docs/README.md
git commit -m "$(cat <<'EOF'
docs(restructure): add dev-docs/README.md index

One-screen guide to what lives under dev-docs/ — ADRs, specs/plans,
archive, historical reports — so future contributors can find moved
content quickly.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Create section landing pages

**Files:**
- Create: `docs/tutorial/index.md`
- Create: `docs/reference/index.md`
- Create: `docs/developer/index.md`

- [ ] **Step 1: Create `docs/tutorial/index.md`**

```bash
mkdir -p docs/tutorial
```

Write `docs/tutorial/index.md` with:

```markdown
---
title: Tutorial
---

# Tutorial

The end-to-end SWAP tutorial is in progress. It will walk through one full
simulation — preparing input files, authoring `swap.toml`, running the model,
and reading outputs.

In the meantime, head over to the [Reference](../reference/index.html) section
for per-topic pages on configuration, meteorology, and input formats.
```

- [ ] **Step 2: Create `docs/reference/index.md`**

Write `docs/reference/index.md` with:

```markdown
---
title: Reference
---

# Reference

Per-topic reference for authoring a SWAP run.

- **[Configuration schema](configuration-schema.html)** — the full TOML key
  reference, organised by section.
- **[TOML format guide](toml-format-guide.html)** — conventions, examples,
  and gotchas for authoring `swap.toml`.
- **[Meteorology](meteorology.html)** — `[meteorology]` section options,
  reference evapotranspiration methods, daily meteo CSV format.
- **[CSV companion files](csv-companion-files.html)** — when to use a CSV
  alongside `swap.toml`, expected schemas, the `read_csv_table` adapter.
```

- [ ] **Step 3: Create `docs/developer/index.md`**

Write `docs/developer/index.md` with:

```markdown
---
title: Developer
---

# Developer

For contributors to the SWAP modernization.

## Orientation

- **[Architecture](architecture.html)** — three-phase control flow, state
  aggregation, module boundaries.
- **[Contributing](contributing.html)** — branch model, commit conventions,
  verification gates.
- **[Build and test](build-and-test.html)** — building locally, fast/full
  test split, regression harness.
- **[Code style](code-style.html)** — Fortran 2008 conventions matching
  `.fprettify.rc`.

## Subsystems

- **[State management](state-management.html)** — `config_t` / `initial_t`
  / `state_t` lifecycle and the ASSOCIATE pattern.
- **[Error handling](error-handling.html)** — error collection, calling
  convention, abort checkpoint.
- **[Logging](logging.html)** — `swap_log` facility, levels, format.
- **[Validation](validation.html)** — primitive checks, section validators,
  how to add new ones.

## Operations

- **[Branches](branches.html)** — branch model and naming.
- **[Dependency management](dependency-management.html)** — subprojects,
  pixi deps, pFUnit, how to bump.
```

- [ ] **Step 4: Verify and commit**

```bash
ls docs/tutorial/ docs/reference/ docs/developer/
git add docs/tutorial/ docs/reference/index.md docs/developer/index.md
git commit -m "$(cat <<'EOF'
docs(restructure): add section landers for tutorial/reference/developer

Tutorial is a stub pointing to Reference until Phase C populates it.
Reference and Developer are short overviews with bulleted page lists
grouped by purpose.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Update `docs.md` and delete stale files

**Files:**
- Modify: `docs.md` (remove one line from `exclude_dir`)
- Delete: `docs/templates/_index.html`
- Delete: `swap_debug.log`

- [ ] **Step 1: Edit `docs.md`**

Currently `docs.md` line 12-13 reads:
```
exclude_dir: ./subprojects
             ./docs/superpowers
```

Use the Edit tool to remove the second line:

- **old_string:**
  ```
  exclude_dir: ./subprojects
               ./docs/superpowers
  ```
- **new_string:**
  ```
  exclude_dir: ./subprojects
  ```

- [ ] **Step 2: Delete stale files**

```bash
git rm docs/templates/_index.html
rm -f swap_debug.log
```

(`swap_debug.log` is untracked, so it's `rm` not `git rm`.)

- [ ] **Step 3: Rewrite `docs/index.md`**

The current `docs/index.md` is built around Phase 4 framing and links to the old flat structure. Replace its body entirely. Use the Read tool first to confirm current content, then Write the new content:

```markdown
---
title: SWAP documentation
author: SWAP Team
---

# SWAP — Soil Water Atmosphere Plant

SWAP is a one-dimensional vertical simulation model for transport processes
(water, heat, solutes) in the soil–plant–atmosphere continuum at field scale.
This site is the home of the modernization of SWAP 4.2.0.

## For users

- **[Tutorial](tutorial/index.html)** — work through a full SWAP run end to end.
- **[Reference](reference/index.html)** — per-topic pages on configuration,
  meteorology, and input formats.

## For contributors

- **[Developer](developer/index.html)** — architecture, contributing guide,
  build/test, subsystem docs.

## API reference

The Fortran source is browsable at [api/index.html](api/index.html). Generated
by FORD; rebuild with `pixi run -e docs docs-build`.
```

- [ ] **Step 4: Verify and commit**

```bash
git status --short
```

Expected: `docs.md`, `docs/index.md`, and `docs/templates/_index.html` (deleted) are all in the staged or unstaged sections.

```bash
git add docs.md docs/index.md docs/templates/_index.html
git commit -m "$(cat <<'EOF'
docs(restructure): wire FORD config + landing page to new layout

- docs.md: drop ./docs/superpowers from exclude_dir (no longer present)
- docs/index.md: full rewrite, drops Phase-4 framing and ADR list,
  points at Tutorial / Reference / Developer / API
- Delete docs/templates/_index.html (stale stock-template backup)

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Intermediate verification — clean structural build

This is the **Phase A done gate**. No template work proceeds until the build is clean.

- [ ] **Step 1: Run the docs build**

```bash
pixi run -e docs docs-rebuild 2>&1 | tee /tmp/docs-build.log | tail -30
```

Expected: exit code 0, last line `Browse the generated documentation: ...`.

- [ ] **Step 2: Check for forbidden warnings**

```bash
grep -E "no title metadata|Skipping creating page|Ignoring unknown Ford metadata key" /tmp/docs-build.log
```

Expected: **empty output**. If any line matches, the cleanup is incomplete — read the warning, find the offending file, fix it (add frontmatter or move it out of `docs/`), then re-run from Step 1.

- [ ] **Step 3: Check the generated search index for leaks**

```bash
pixi run -e docs python -c "
import json, pathlib
raw = pathlib.Path('docs/api/search/search_database.json').read_text()
pages = json.loads(raw.removeprefix('var tipuesearch = '))['pages']
urls  = [p['url'] for p in pages]
for forbidden in ('/adr/', '/superpowers/', '/archive/'):
    leaked = [u for u in urls if forbidden in u]
    assert not leaked, f'{forbidden!r} leaked into search index: {leaked[:3]}'
assert any('/reference/' in u for u in urls), 'reference pages missing from index'
assert any('/developer/' in u for u in urls), 'developer pages missing from index'
print(f'OK: {len(pages)} pages indexed')
"
```

Expected: `OK: <N> pages indexed`. If any assertion fires, fix the offending move/exclude and re-build.

- [ ] **Step 4: Manual smoke (no commit yet)**

```bash
pixi run -e docs docs-serve &
SERVE_PID=$!
sleep 2
curl -sf http://localhost:8001/index.html > /dev/null && echo "landing OK"
curl -sf http://localhost:8001/page/reference/meteorology.html > /dev/null && echo "reference page OK"
curl -sf http://localhost:8001/page/developer/architecture.html > /dev/null && echo "developer page OK"
curl -s -o /dev/null -w "%{http_code}\n" http://localhost:8001/page/adr/0001-gfortran-first.html | grep -q 404 && echo "ADR correctly 404"
kill $SERVE_PID 2>/dev/null
```

Expected output:
```
landing OK
reference page OK
developer page OK
ADR correctly absent
```

- [ ] **Step 5: Mark the milestone**

No new commit (no source changes since Task 7). Tag the current commit locally so we can find the gate later:

```bash
git tag docs-restructure-phase-a-clean
```

> **Stop here for orchestrator review** before proceeding to template work.

---

## Task 9: Write `docs/templates/base.html` (Atlas navbar, no search JS yet)

**Files:**
- Create: `docs/templates/base.html`

- [ ] **Step 1: Read the existing index.html and stock base.html for reference**

```bash
sed -n '1,50p' .pixi/envs/docs/lib/python3.11/site-packages/ford/templates/base.html
```

This is the file we're replacing — keep the `<head>` block (Bootstrap CSS, MathJax, FORD's local.css, etc.) intact and only restyle the `<body>` chrome.

- [ ] **Step 2: Write the new base.html**

Create `docs/templates/base.html` with the following content:

```html
{# ============================================================
   SWAP · FORD base template (V2 · Atlas)
   Replaces FORD's stock base.html. Provides the Atlas-styled
   navbar + inline typeahead search on every generated page.
   ============================================================ #}
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    {% if summary %}
      <meta name="description" content="{{ summary|striptags }}">
    {% else %}
      <meta name="description" content="Documentation for {{ project }}">
    {% endif %}
    <meta name="author" content="{{ author }}" >
    <link rel="icon" href="{{ project_url }}/favicon.png">

    <title>{% block title %} {{ project }} {% endblock title %}</title>

    <!-- Bootstrap (kept for FORD inner-page grid classes) -->
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet"
          integrity="sha384-QWTKZyjpPEjISv5WaRU9OFeRpok6YctnYmDr5pNlyT2bRjXh0JMhjY6hW+ALEwIH" crossorigin="anonymous">
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/js/bootstrap.bundle.min.js"
            integrity="sha384-YvpcrYf0tY3lHB60NNkmXc5s9fDVZLESaAA55NDzOxhy9GkcIdslK1eN7N6jIeHz" crossorigin="anonymous"></script>
    <!-- Font Awesome -->
    <link href="{{ project_url }}/css/fontawesome.min.css" rel="stylesheet">
    <link href="{{ project_url }}/css/brands.min.css" rel="stylesheet">
    <link href="{{ project_url }}/css/regular.min.css" rel="stylesheet">
    <link href="{{ project_url }}/css/solid.min.css" rel="stylesheet">
    <link href="{{ project_url }}/css/v4-font-face.min.css" rel="stylesheet">
    <link href="{{ project_url }}/css/v4-shims.min.css" rel="stylesheet">
    <!-- MathJax -->
    <script type="text/x-mathjax-config">
      MathJax.Hub.Config({ TeX: { equationNumbers: { autoNumber: "AMS" } } });
    </script>
    {% if mathjax_config %}
      <script src="{{ project_url }}/js/MathJax-config/{{ path.basename(mathjax_config) }}"></script>
    {% endif %}
    <script src="https://cdn.jsdelivr.net/npm/mathjax@2.7.9/MathJax.js?config=TeX-AMS-MML_HTMLorMML" async
            integrity="sha256-DViIOMYdwlM/axqoGDPeUyf0urLoHMN4QACBKyB58Uw=" crossorigin="anonymous"></script>
    <!-- FORD-provided styles -->
    <link href="{{ project_url }}/css/local.css" rel="stylesheet">
    <link href="{{ project_url }}/css/pygments.css" rel="stylesheet">
    {% if css %}
      <link href="{{ project_url }}/css/user.css" rel="stylesheet">
    {% endif %}
    <script src="{{ project_url }}/js/svg-pan-zoom.min.js"></script>

    <!-- Atlas palette + fonts (inherited by all child pages) -->
    <link href="https://fonts.googleapis.com/css2?family=Spectral:ital,wght@0,400;0,500;1,400;1,500&family=Inter:wght@400;500;600&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
    <style>
      :root {
        --paper:      #f6f1e7;
        --paper-alt:  #efe8d8;
        --card:       #fdfaf3;
        --ink:        #1d1916;
        --ink-soft:   #3a322b;
        --muted:      #7a6f60;
        --faint:      #a89e8c;
        --rule:       #d9cfb8;
        --rule-soft:  #e6dcc4;
        --accent:     #b85535;
        --serif: 'Spectral', 'Source Serif Pro', Georgia, 'Times New Roman', serif;
        --sans:  'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', system-ui, sans-serif;
        --mono:  'JetBrains Mono', 'IBM Plex Mono', ui-monospace, Menlo, Consolas, monospace;
      }
      body { background: var(--paper); color: var(--ink); font-family: var(--sans); }

      /* ---------- Atlas navbar ---------- */
      .swap-nav {
        position: sticky; top: 0; z-index: 100;
        display: flex; align-items: center; gap: 24px;
        padding: 14px 32px;
        background: var(--paper);
        border-bottom: 1px solid var(--rule);
        font-family: var(--sans);
      }
      .swap-nav__brand {
        font-family: var(--serif); font-size: 20px; font-weight: 500;
        color: var(--ink); text-decoration: none; letter-spacing: -0.01em;
      }
      .swap-nav__brand:hover { color: var(--accent); }
      .swap-nav__links {
        list-style: none; margin: 0; padding: 0;
        display: flex; align-items: center; gap: 4px;
      }
      .swap-nav__links a, .swap-nav__dropdown > button {
        font-family: var(--sans); font-size: 14px; font-weight: 500;
        color: var(--ink-soft); text-decoration: none;
        padding: 6px 12px; border: none; background: transparent;
        border-radius: 2px; cursor: pointer;
      }
      .swap-nav__links a:hover, .swap-nav__dropdown > button:hover {
        color: var(--accent); background: var(--paper-alt);
      }
      .swap-nav__divider {
        width: 1px; height: 16px; background: var(--rule); margin: 0 4px;
      }
      .swap-nav__dropdown { position: relative; }
      .swap-nav__menu {
        position: absolute; top: 100%; left: 0; margin-top: 4px;
        background: var(--card); border: 1px solid var(--rule);
        border-radius: 2px; min-width: 180px; padding: 6px 0;
        display: none; box-shadow: 0 4px 12px rgba(0,0,0,0.05);
      }
      .swap-nav__menu a {
        display: block; padding: 8px 16px; font-size: 13px;
        color: var(--ink-soft); text-decoration: none;
      }
      .swap-nav__menu a:hover { background: var(--paper-alt); color: var(--accent); }
      .swap-nav__dropdown.is-open .swap-nav__menu { display: block; }

      .swap-nav__search {
        margin-left: auto; position: relative;
        display: flex; align-items: center; gap: 8px;
        padding: 6px 12px;
        background: var(--card); border: 1px solid var(--rule);
        border-radius: 2px; min-width: 280px;
      }
      .swap-nav__search:focus-within { border-color: var(--accent); }
      .swap-nav__search input {
        flex: 1; border: none; outline: none; background: transparent;
        font-family: var(--sans); font-size: 14px; color: var(--ink);
      }
      .swap-nav__search input::placeholder { color: var(--muted); }
      .swap-nav__search kbd {
        font-family: var(--mono); font-size: 11px; color: var(--muted);
        border: 1px solid var(--rule); border-radius: 3px;
        padding: 2px 6px; background: var(--paper);
      }
      #swap-search-results {
        position: absolute; top: 100%; left: 0; right: 0; margin-top: 4px;
        background: var(--card); border: 1px solid var(--rule);
        border-radius: 2px; max-height: 360px; overflow-y: auto;
        box-shadow: 0 4px 12px rgba(0,0,0,0.05);
      }
      #swap-search-results a {
        display: block; padding: 10px 14px;
        color: var(--ink); text-decoration: none;
        border-bottom: 1px solid var(--rule-soft);
        font-size: 13px;
      }
      #swap-search-results a:last-child { border-bottom: none; }
      #swap-search-results a:hover, #swap-search-results a.is-active {
        background: var(--paper-alt); color: var(--accent);
      }
      #swap-search-results .swap-search-title { font-weight: 500; }
      #swap-search-results .swap-search-snippet {
        font-size: 12px; color: var(--muted); margin-top: 2px;
      }
      #swap-search-results .swap-search-empty {
        padding: 12px 14px; color: var(--muted); font-size: 13px;
      }

      /* ---------- Footer ---------- */
      .swap-footer {
        margin-top: 64px; padding: 24px 32px;
        border-top: 1px solid var(--rule);
        background: var(--paper-alt);
        font-family: var(--mono); font-size: 11px; color: var(--muted);
        letter-spacing: 0.06em;
        display: flex; justify-content: space-between;
      }
      .swap-footer a { color: inherit; text-decoration: none; border-bottom: 1px solid var(--rule); }
      .swap-footer a:hover { color: var(--accent); border-bottom-color: var(--accent); }

      @media (max-width: 800px) {
        .swap-nav { flex-wrap: wrap; padding: 12px 16px; gap: 12px; }
        .swap-nav__search { order: 99; width: 100%; min-width: 0; }
        .swap-footer { flex-direction: column; gap: 8px; padding: 16px; }
      }
    </style>
  </head>

  <body>
    <nav class="swap-nav">
      <a class="swap-nav__brand" href="{{ project_url }}/index.html">{{ project }}{% if version %} <small style="font-size:12px;color:var(--muted);">{{ version }}</small>{% endif %}</a>
      <ul class="swap-nav__links">
        <li><a href="{{ project_url }}/page/tutorial/index.html">Tutorial</a></li>
        <li><a href="{{ project_url }}/page/reference/index.html">Reference</a></li>
        <li><a href="{{ project_url }}/page/developer/index.html">Developer</a></li>
        <li class="swap-nav__divider"></li>
        <li class="swap-nav__dropdown" id="swap-api-dropdown">
          <button type="button" aria-haspopup="true" aria-expanded="false">API ▾</button>
          <div class="swap-nav__menu" role="menu">
            {% if project.modules %}<a href="{{ project_url }}/lists/modules.html">Modules</a>{% endif %}
            {% if project.procedures %}<a href="{{ project_url }}/lists/procedures.html">Procedures</a>{% endif %}
            {% if project.types %}<a href="{{ project_url }}/lists/types.html">Derived Types</a>{% endif %}
            {% if project.absinterfaces %}<a href="{{ project_url }}/lists/absint.html">Abstract Interfaces</a>{% endif %}
            {% if project.namelists %}<a href="{{ project_url }}/lists/namelists.html">Namelists</a>{% endif %}
            {% if incl_src %}<a href="{{ project_url }}/lists/files.html">Source Files</a>{% endif %}
          </div>
        </li>
      </ul>
      {% if search %}
      <div class="swap-nav__search">
        <input type="search" id="swap-search-input" placeholder="Search the docs…" autocomplete="off">
        <kbd>/</kbd>
        <div id="swap-search-results" hidden></div>
      </div>
      {% endif %}
    </nav>

    <div class="container">
      {% block body %}
      {% endblock body %}
    </div>

    <footer class="swap-footer">
      <span>
        {{ project }}{% if author %} · developed by {{ author }}{% endif %}
        · &copy; {{ year }}{% if license %} {{ license }}{% endif %}
      </span>
      <span>
        Generated by <a href="https://github.com/Fortran-FOSS-Programmers/ford">FORD</a>
        {% if print_creation_date %}· {{ creation_date }}{% endif %}
      </span>
    </footer>

    <script>
      // API dropdown toggle
      (function () {
        var dd = document.getElementById('swap-api-dropdown');
        if (!dd) return;
        var btn = dd.querySelector('button');
        btn.addEventListener('click', function (e) {
          e.stopPropagation();
          var open = dd.classList.toggle('is-open');
          btn.setAttribute('aria-expanded', open ? 'true' : 'false');
        });
        document.addEventListener('click', function () {
          dd.classList.remove('is-open');
          btn.setAttribute('aria-expanded', 'false');
        });
      })();
    </script>

  </body>
</html>
```

- [ ] **Step 3: Rebuild and confirm the navbar renders**

```bash
pixi run -e docs docs-rebuild 2>&1 | tail -5
```

Expected: exit code 0, last line `Browse the generated documentation: ...`.

```bash
grep -q "swap-nav__brand" docs/api/index.html && echo "navbar rendered on landing"
grep -q "swap-nav__brand" docs/api/page/reference/meteorology.html && echo "navbar rendered on reference page"
grep -q "swap-nav__brand" docs/api/lists/modules.html && echo "navbar rendered on modules page"
```

Expected: all three lines print. If any is missing, the template isn't being picked up — re-check that `html_template_dir: ./docs/templates` is still in `docs.md` and `docs/templates/base.html` exists.

- [ ] **Step 4: Commit**

```bash
git add docs/templates/base.html
git commit -m "$(cat <<'EOF'
docs(templates): add Atlas-styled base.html overriding FORD navbar

Custom navbar with brand + Tutorial/Reference/Developer/API dropdown,
inline search input slot (typeahead JS in next commit), Atlas palette
inherited from index.html. Bootstrap CSS kept in <head> so FORD's
inner-page grid classes still render.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: Add typeahead search to `docs/templates/base.html`

**Files:**
- Modify: `docs/templates/base.html` (append lunr.js script + search JS)

- [ ] **Step 1: Add lunr.js CDN and search JS before `</body>`**

Use the Edit tool on `docs/templates/base.html`. Locate the existing closing block:

- **old_string:**
  ```
    </script>

  </body>
  </html>
  ```
- **new_string:**
  ```
    </script>

    {% if search %}
    <script src="https://unpkg.com/lunr@2.3.9/lunr.min.js"
            integrity="sha384-jHjQy43yPI0lvtq28QhFAehv4OmtPq2VxlfFwkQbvjbflnpKAEhmJgsg5SADQyrf"
            crossorigin="anonymous"></script>
    <script>
      // Inline typeahead search backed by FORD's existing search_database.json.
      (function () {
        var input   = document.getElementById('swap-search-input');
        var results = document.getElementById('swap-search-results');
        if (!input || !results) return;

        var idx = null, pageMap = null, loaded = false, loading = false;
        var dbUrl = '{{ project_url }}/search/search_database.json';

        function loadDb() {
          if (loaded || loading) return;
          loading = true;
          fetch(dbUrl).then(function (r) { return r.text(); }).then(function (txt) {
            var data = JSON.parse(txt.replace(/^var tipuesearch\s*=\s*/, ''));
            pageMap = data.pages;
            idx = lunr(function () {
              this.ref('id');
              this.field('title', { boost: 10 });
              this.field('text');
              this.field('tags', { boost: 5 });
              var self = this;
              pageMap.forEach(function (p, i) {
                self.add({
                  id:    i,
                  title: p.title || '',
                  text:  p.text  || '',
                  tags:  p.tags  || ''
                });
              });
            });
            loaded = true; loading = false;
            // If the user typed something while we were loading, re-query.
            if (input.value) run(input.value);
          }).catch(function (err) {
            loading = false;
            console.error('search db load failed:', err);
          });
        }

        function escapeHtml(s) {
          return String(s).replace(/[&<>"']/g, function (c) {
            return ({ '&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'})[c];
          });
        }

        function snippet(text, query) {
          if (!text) return '';
          var i = text.toLowerCase().indexOf(query.toLowerCase());
          if (i < 0) return text.slice(0, 100) + (text.length > 100 ? '…' : '');
          var start = Math.max(0, i - 40);
          var end   = Math.min(text.length, i + query.length + 60);
          return (start > 0 ? '…' : '') + text.slice(start, end) + (end < text.length ? '…' : '');
        }

        function run(q) {
          if (!loaded) return;
          q = q.trim();
          if (!q) { hide(); return; }
          var hits;
          try { hits = idx.search(q + '*'); }
          catch (e) { hits = []; }
          hits = hits.slice(0, 8);
          if (!hits.length) {
            results.innerHTML = '<div class="swap-search-empty">No matches for "' + escapeHtml(q) + '"</div>';
          } else {
            results.innerHTML = hits.map(function (h, i) {
              var p = pageMap[parseInt(h.ref, 10)];
              return '<a href="' + escapeHtml(p.url) + '"' + (i === 0 ? ' class="is-active"' : '') + '>' +
                       '<div class="swap-search-title">' + escapeHtml(p.title || p.url) + '</div>' +
                       '<div class="swap-search-snippet">' + escapeHtml(snippet(p.text, q)) + '</div>' +
                     '</a>';
            }).join('');
          }
          results.hidden = false;
        }

        function hide() { results.hidden = true; }

        input.addEventListener('focus',  loadDb);
        input.addEventListener('input',  function () { run(input.value); });
        input.addEventListener('keydown', function (e) {
          if (e.key === 'Enter') {
            var first = results.querySelector('a');
            if (first) { e.preventDefault(); window.location.href = first.href; }
          } else if (e.key === 'Escape') {
            hide(); input.blur();
          }
        });
        document.addEventListener('click', function (e) {
          if (!results.contains(e.target) && e.target !== input) hide();
        });
        document.addEventListener('keydown', function (e) {
          if (e.key === '/' && !/^(input|textarea|select)$/i.test(document.activeElement.tagName)) {
            e.preventDefault();
            input.focus();
          }
        });
      })();
    </script>
    {% endif %}

  </body>
  </html>
  ```

- [ ] **Step 2: Rebuild**

```bash
pixi run -e docs docs-rebuild 2>&1 | tail -5
```

Expected: exit code 0.

- [ ] **Step 3: Verify lunr is loaded on at least the landing page and a deep inner page**

```bash
grep -l "lunr@2.3.9" docs/api/index.html docs/api/page/reference/meteorology.html docs/api/lists/modules.html
```

Expected: all three filenames printed.

- [ ] **Step 4: Commit**

```bash
git add docs/templates/base.html
git commit -m "$(cat <<'EOF'
docs(templates): inline typeahead search with lunr.js

Lazy-loads search_database.json on first input focus, builds a lunr
index, renders top-8 hits in a dropdown beneath the navbar input.
Keyboard: Enter goes to first result, '/' focuses the input from
anywhere, Escape closes the dropdown.

Reuses FORD's existing search_database.json — no build-config changes.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 11: Final verification

This is the **Phase A+B done gate**.

- [ ] **Step 1: Clean build**

```bash
pixi run -e docs docs-rebuild 2>&1 | tee /tmp/docs-build-final.log | tail -5
```

Expected: exit code 0.

- [ ] **Step 2: Zero forbidden warnings**

```bash
grep -E "no title metadata|Skipping creating page|Ignoring unknown Ford metadata key|AttributeError" /tmp/docs-build-final.log
```

Expected: empty output.

- [ ] **Step 3: Search index sanity**

```bash
pixi run -e docs python -c "
import json, pathlib
raw = pathlib.Path('docs/api/search/search_database.json').read_text()
pages = json.loads(raw.removeprefix('var tipuesearch = '))['pages']
urls  = [p['url'] for p in pages]
for forbidden in ('/adr/', '/superpowers/', '/archive/', 'phase-4-modernization-summary', 'coverage-baseline'):
    leaked = [u for u in urls if forbidden in u]
    assert not leaked, f'{forbidden!r} leaked: {leaked[:3]}'
assert any('/reference/' in u for u in urls), 'reference missing'
assert any('/developer/' in u for u in urls), 'developer missing'
assert any('/tutorial/' in u for u in urls), 'tutorial missing'
print(f'OK: {len(pages)} pages indexed; no dev material leaked')
"
```

Expected: `OK: <N> pages indexed; no dev material leaked`.

- [ ] **Step 4: Served-site smoke test**

```bash
pixi run -e docs docs-serve &
SERVE_PID=$!
sleep 2

# Landing page returns 200 and contains the Atlas brand
curl -sf http://localhost:8001/index.html | grep -q "swap-nav__brand" && echo "landing nav OK"

# Reference page renders with the new nav
curl -sf http://localhost:8001/page/reference/meteorology.html | grep -q "swap-nav__brand" && echo "reference nav OK"

# Developer page renders with the new nav
curl -sf http://localhost:8001/page/developer/architecture.html | grep -q "swap-nav__brand" && echo "developer nav OK"

# Module page renders with the new nav (inner FORD page inheritance)
MOD=$(ls docs/api/module/ 2>/dev/null | head -1)
[ -n "$MOD" ] && curl -sf "http://localhost:8001/module/$MOD" | grep -q "swap-nav__brand" && echo "module nav OK"

# Search assets are loaded on inner pages
curl -sf http://localhost:8001/page/reference/meteorology.html | grep -q "lunr@2.3.9" && echo "search JS loaded inner page"

# Old dev paths return 404
curl -s -o /dev/null -w "%{http_code}\n" http://localhost:8001/page/adr/0001-gfortran-first.html | grep -q 404 && echo "ADR correctly 404"

kill $SERVE_PID 2>/dev/null
wait $SERVE_PID 2>/dev/null
```

Expected output (order may vary slightly):
```
landing nav OK
reference nav OK
developer nav OK
module nav OK
search JS loaded inner page
ADR correctly 404
```

If any line is missing, investigate before the final commit. Common causes:
- `module nav OK` missing → `docs/api/module/` is empty (FORD didn't index any modules) or path differs; check `ls docs/api/module/` and adjust.
- `search JS loaded` missing → `search:` is not `true` in `docs.md`; check the `{% if search %}` gating.

- [ ] **Step 5: Internal-link audit**

```bash
grep -rn "docs/adr\|docs/superpowers\|docs/archive\|PHASE-4-MODERNIZATION-SUMMARY" docs/ src/ README.md CLAUDE.md 2>/dev/null | grep -v "docs/api/" | grep -v "dev-docs/"
```

Expected: empty output. (Matches under `docs/api/` are FORD-generated; matches under `dev-docs/` are internal cross-references between moved files, both expected.)

- [ ] **Step 6: Manual visual check (orchestrator only — not subagent)**

Open `http://localhost:8001/` in a browser. Confirm:
- Landing page Atlas styling renders correctly.
- Navbar shows: brand, Tutorial / Reference / Developer / API ▾, search input with `/` kbd hint.
- Click "Reference" → reference lander renders, navbar persists.
- Type "meteorology" in the search → dropdown shows the page within ~200ms, click navigates.
- Press `/` on any page → search input is focused.
- Hover "API ▾" → dropdown shows Modules / Procedures / Derived Types / Source Files.

- [ ] **Step 7: Final tag**

```bash
git tag docs-restructure-phase-b-complete
```

No new commit (no source changes since Task 10).

- [ ] **Step 8: Summary report**

Print a final report (orchestrator step, no commit):

```bash
echo "=== Docs restructure — Phase A+B complete ==="
echo "Commits added this arc:"
git log --oneline docs-restructure-phase-a-clean^..docs-restructure-phase-b-complete 2>/dev/null || \
  git log --oneline -20
echo
echo "Build: $(pixi run -e docs docs-rebuild 2>&1 | tail -1)"
echo "Tag (Phase A): docs-restructure-phase-a-clean"
echo "Tag (Phase B): docs-restructure-phase-b-complete"
```

---

## Out of scope (Phase C — separate spec)

- Tutorial narrative authoring (currently a one-paragraph stub).
- Per-page TOC sidebars for current-section navigation.
- Customising FORD's API-side templates (`mod_page.html`, `proc_page.html`, `type_page.html`).
- Versioned docs (multiple SWAP releases on one site).
- ADR publishing pipeline (if we later decide some ADRs deserve a public home).
