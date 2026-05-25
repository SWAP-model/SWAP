# Docs site restructure — Phase A+B (structure + templates)

**Date:** 2026-05-25
**Status:** Approved by user (brainstorming complete); awaiting implementation plan.
**Scope:** Phases A (file reorganisation) and B (template overrides) of a larger
docs cleanup. Phase C (Tutorial + Reference content authoring) is a separate
future spec.

---

## Problem

The FORD-generated documentation site has accumulated cruft over the rescue
arcs:

- Dev-only material lives in `docs/` alongside user-facing pages:
  42 ADRs, frozen Phase-4 archive, `superpowers/` specs/plans/notes,
  status reports (`PHASE-4-MODERNIZATION-SUMMARY.md`, `coverage-baseline.md`).
- Two pages (`csv-companion-files.md`, `meteorology.md`) lack the YAML
  `title:` frontmatter FORD requires, producing parse errors on every build.
- ADRs carry custom metadata keys (`status:`, `supersedes-portion-of:`)
  FORD doesn't recognise, producing 30+ warnings per build.
- `docs/index.md` still frames the project around "Phase 4" (which completed
  2026-05-05) and links to only ADRs 0001–0006 of the 42 that exist.
- A new Atlas-style `index.html` template was added, but FORD's stock
  Bootstrap-dark navbar still wraps every page, creating visual dissonance.
- The site search submits to a dedicated results page — feels dated next
  to the new landing page.
- `docs/templates/_index.html` is a stale backup of the stock FORD template.
- `swap_debug.log` is an untracked debug artefact at repo root.

Underlying cause: the docs site has no taxonomy. Dev artefacts, contributor
guides, and user references all sit at the same depth.

## Goals

1. The public site contains only user-facing and contributor-facing content.
2. Three clear top-level audiences: **Tutorial**, **Reference**, **Developer**,
   plus the FORD-generated **API**.
3. `pixi run -e docs docs-rebuild` exits 0 with **zero** FORD warnings about
   missing titles, skipped directories, or unknown metadata keys.
4. The navbar and search bar match the new Atlas palette and feel modern
   (inline typeahead, `/` shortcut).
5. Dev material (ADRs, specs, plans, archive, status reports) is preserved
   in the repo at a known location, just out of FORD's scan path.

## Non-goals

- **Authoring new tutorial prose.** A stub `docs/tutorial/index.md` is enough
  for Phase A+B. The narrative walkthrough is Phase C.
- **Rewriting existing page bodies.** Frontmatter fixes only.
- **Changing FORD's API-side templates** (`mod_page.html`, `proc_page.html`,
  etc.). Inner pages keep FORD defaults; only the chrome (`base.html`)
  changes.
- **Dropping Bootstrap CSS.** FORD's generated inner pages use Bootstrap
  grid classes — Bootstrap stays in `<head>`, just isn't used by the new
  navbar.
- **Modifying `pixi.toml` docs tasks.** Already verified: `docs-build`,
  `docs-clean`, `docs-serve`, `docs-open` reference only `docs/api/`
  (FORD output), no source paths.

---

## Design

### 1. Directory layout

```
swap/
├── docs/                          ← public site source (FORD reads this)
│   ├── index.md                   (rewritten; no Phase-4 framing)
│   ├── tutorial/
│   │   └── index.md               (stub; Phase C populates)
│   ├── reference/                 (user-facing)
│   │   ├── index.md               (section lander)
│   │   ├── configuration-schema.md
│   │   ├── toml-format-guide.md
│   │   ├── meteorology.md         (+ title frontmatter)
│   │   └── csv-companion-files.md (+ title frontmatter)
│   ├── developer/                 (contributor-facing)
│   │   ├── index.md               (section lander)
│   │   ├── architecture.md
│   │   ├── contributing.md
│   │   ├── code-style.md
│   │   ├── build-and-test.md
│   │   ├── state-management.md
│   │   ├── error-handling.md
│   │   ├── logging.md
│   │   ├── validation.md
│   │   ├── branches.md
│   │   └── dependency-management.md
│   ├── templates/
│   │   ├── base.html              (NEW — Atlas navbar + typeahead search)
│   │   └── index.html             (existing Atlas template, unchanged)
│   └── public/                    (favicon, images — unchanged)
│
├── dev-docs/                      ← NEW; internal only; NOT in FORD's scan
│   ├── README.md                  (NEW; index of what lives here)
│   ├── adr/                       (all 42 ADRs, moved as a tree)
│   ├── superpowers/               (specs, plans, notes — moved as a tree;
│   │                               this design doc moves here as part of
│   │                               implementation)
│   ├── archive/2026-phase-4/      (moved as a tree)
│   ├── phase-4-modernization-summary.md
│   └── coverage-baseline.md
│
└── docs.md                        ← FORD config; one targeted edit
```

**Files deleted (not moved):**

- `docs/templates/_index.html` — stale backup of the stock FORD index.
- `swap_debug.log` — untracked debug artefact at repo root.

**`docs.md` changes:**

- Remove the `./docs/superpowers` entry from `exclude_dir`. After the move,
  no subdirectory of `docs/` needs to be excluded.
- No other edits. FORD walks `docs/tutorial/`, `docs/reference/`,
  `docs/developer/` automatically.

**Cross-reference updates:**

- `docs/index.md` is rewritten end-to-end — drops the ADR list, the Phase-4
  framing, and the dev-page links. Points to Tutorial / Reference / Developer
  / API only.
- ADR-internal cross-references (ADRs that link to sibling ADRs) keep working
  because the whole tree moves together.
- No source code references any of the moved paths — verified via
  `grep -rn "docs/adr\|docs/superpowers\|docs/archive" src/ pixi.toml README.md`
  before the move (implementation step).

### 2. Frontmatter normalisation

Every relocated `.md` page is checked for a YAML frontmatter block with a
`title:` field. Two known offenders today:

- `docs/csv-companion-files.md` → moved to `docs/reference/csv-companion-files.md`
  with `title: CSV companion files` prepended.
- `docs/meteorology.md` → moved to `docs/reference/meteorology.md`
  with `title: Meteorology` prepended.

The implementation plan includes a `grep -L "^title:" <moved-md-files>` check
to catch any other missing-title pages before the build.

### 3. Custom `docs/templates/base.html`

Replaces FORD's stock `base.html`. Inherits the Jinja2 contract (extends
nothing; defines `{% block title %}` and `{% block body %}`) so FORD's
inner-page templates (`mod_page.html`, `proc_page.html`, etc.) keep working.

**`<head>` section — kept from FORD's stock:**

- Bootstrap CSS bundle (used by FORD's inner pages)
- Bootstrap JS bundle
- FontAwesome (icons in FORD inner pages)
- MathJax (equations)
- FORD's `local.css` and `pygments.css`
- Favicon, meta description, viewport

**`<head>` section — added:**

- Atlas palette CSS variables on `:root` (extracted from `index.html` so
  every page inherits them).
- lunr.js `<script>` (CDN; already used by FORD on `search.html`).
- Lazy loader for `search/search_database.json` (deferred to first focus).

**Navbar — replaces the Bootstrap dark navbar entirely:**

```html
<nav class="swap-nav">
  <a class="swap-nav__brand" href="{{ project_url }}/index.html">{{ project }}</a>
  <ul class="swap-nav__links">
    <li><a href="{{ project_url }}/tutorial/index.html">Tutorial</a></li>
    <li><a href="{{ project_url }}/reference/index.html">Reference</a></li>
    <li><a href="{{ project_url }}/developer/index.html">Developer</a></li>
    <li class="swap-nav__divider"></li>
    <li class="swap-nav__dropdown">
      <button type="button" aria-haspopup="true" aria-expanded="false">API ▾</button>
      <div class="swap-nav__menu" role="menu">
        {% if project.modules %}<a href="{{ project_url }}/lists/modules.html">Modules</a>{% endif %}
        {% if project.procedures %}<a href="{{ project_url }}/lists/procedures.html">Procedures</a>{% endif %}
        {% if project.types %}<a href="{{ project_url }}/lists/types.html">Derived Types</a>{% endif %}
        {% if incl_src %}<a href="{{ project_url }}/lists/files.html">Source Files</a>{% endif %}
      </div>
    </li>
  </ul>
  <div class="swap-nav__search">
    <input type="search" id="swap-search-input" placeholder="Search…" autocomplete="off">
    <kbd>/</kbd>
    <div id="swap-search-results" hidden></div>
  </div>
</nav>
```

**Styling — scoped under `.swap-nav`** so it never collides with Bootstrap
classes on inner pages. Uses the Atlas variables (`--paper`, `--ink`,
`--accent`, `--rule`, etc.) — same palette as `index.html`.

**Typeahead search behaviour:**

1. On first focus of `#swap-search-input`, fetch
   `{{ project_url }}/search/search_database.json` once (it's already
   generated by FORD's `Tipue_Search_JSON_Generator`).
2. The JSON declares `var tipuesearch = { pages: [...] }` — same format
   FORD's own `search.html` consumes.
3. Build a lunr index over `pages[*].{title, text, tags}`.
4. On `input` event, query the index and render the top 8 hits as a
   dropdown of `<a>` links beneath the input.
5. **Enter** → navigate to the first hit. **Click** → navigate to that hit.
   **Escape** or outside-click → close dropdown.
6. **`/` keypress** (when no other input is focused) → focus the search
   input. Mirrors GitHub / Algolia / FORD's own keyboard shortcut convention.
7. FORD's `search.html` remains unchanged as a fallback for users who
   reach it via a deep link.

**Footer — keep FORD's footer markup, restyle minimally** (font-family,
colours) to match the Atlas palette. No content changes.

### 4. Section landing pages

Three small Markdown files are created. Each is a one-screen overview with
a bulleted list of the section's pages. No new prose beyond a sentence or
two of context per page.

- `docs/tutorial/index.md` — title + one sentence: "The end-to-end tutorial
  is in progress. For now, see the **Reference** section."
- `docs/reference/index.md` — title + list of the four reference pages with
  one-line summaries (copy summaries from existing page headings).
- `docs/developer/index.md` — title + list of the ten developer pages,
  grouped: "Orientation" (architecture, contributing, build-and-test,
  code-style), "Subsystems" (state-management, error-handling, logging,
  validation), "Operations" (branches, dependency-management).

### 5. `docs/index.md` rewrite

Replace the current Phase-4-centric content with three short sections:

1. **What is SWAP** — one paragraph (steal/trim from the current summary).
2. **For users** — links to `tutorial/index.html` and `reference/index.html`.
3. **For contributors** — links to `developer/index.html`.
4. **API reference** — link to `api/index.html`.

No mention of phases, rescue arcs, or ADRs in the public landing page.

### 6. `dev-docs/README.md`

A one-screen index of the moved content so future contributors can find
it. Headings: ADRs, Specs & plans (superpowers), Archive (Phase 4),
Historical reports. Each section is a bulleted list with one-line
descriptions of what's inside.

---

## Verification

Phase A+B is "done" when **all four** of these pass:

1. **Build is clean.** `pixi run -e docs docs-rebuild` exits 0, with zero
   `Warning:` lines containing any of:
   - `"no title metadata"`
   - `"Skipping creating page for"`
   - `"Ignoring unknown Ford metadata key"`

2. **Search index has the right pages.** A short Python check:
   ```python
   import json, pathlib
   raw = pathlib.Path("docs/api/search/search_database.json").read_text()
   pages = json.loads(raw.removeprefix("var tipuesearch = "))["pages"]
   urls  = [p["url"] for p in pages]
   assert not any("/adr/" in u for u in urls),          "ADR leaked"
   assert not any("/superpowers/" in u for u in urls),  "superpowers leaked"
   assert not any("/archive/" in u for u in urls),      "archive leaked"
   assert any("/reference/" in u for u in urls),        "reference missing"
   assert any("/developer/" in u for u in urls),        "developer missing"
   ```

3. **Manual smoke test of the served site** (`pixi run -e docs docs-serve`,
   browser to `http://localhost:8001/`):
   - Landing page renders with Atlas styling.
   - Navbar shows: brand, Tutorial / Reference / Developer / API ▾, search input.
   - Typing "meteorology" in the search input shows a typeahead dropdown
     within ~200 ms; Enter navigates to the reference page.
   - Pressing `/` on any page focuses the search input.
   - Clicking into a module page from API ▾ still shows the new navbar
     (template inheritance works on inner pages).

4. **Internal-link integrity.**
   `grep -rn "docs/adr\|docs/superpowers\|docs/archive" docs/ src/ README.md CLAUDE.md`
   returns nothing inside `docs/` source files. (`docs/api/` matches are
   FORD-generated and ignored.)

---

## Risks and mitigations

| Risk | Mitigation |
|---|---|
| Moving ADRs/specs breaks internal cross-references between them | All dev material moves as a tree together; relative links between siblings stay valid. Only `docs/index.md` outbound links to the moved tree need editing. |
| Custom `base.html` breaks FORD's inner page rendering | Keep all `<head>` assets and preserve the `{% block body %}` / `{% block title %}` contract. Only the `<nav>` markup changes. Verified by browsing a module page in the smoke test. |
| lunr.js bundle on every page slows initial load | Defer-load lunr (`<script defer>` from CDN) and lazy-fetch the search database on first input focus, not on page load. |
| Atlas CSS variables only defined in `index.html`'s scoped block | Move the `:root` palette declarations into `base.html`'s `<head>`. `index.html`'s scoped `.swap-index` block stays as-is and inherits. |
| `docs.md` config drift | Single targeted edit: remove `./docs/superpowers` from `exclude_dir`. Diff stays one line + context. |
| Rollback needed | All changes are within `docs/`, `dev-docs/`, and `docs.md`. `git checkout HEAD -- docs/ docs.md && git clean -fd dev-docs/` reverts cleanly. |

---

## Out of scope (deferred to Phase C)

- Authoring the actual end-to-end tutorial narrative.
- Reorganising or rewriting reference page bodies (only frontmatter changes
  here).
- Adding a left-sidebar table of contents for the current section.
- Customising FORD's API-side templates (`mod_page.html`, `proc_page.html`).
- Adding versioned docs (multiple SWAP releases on the same site).

These are tracked for a future spec; nothing in Phase A+B precludes them.
