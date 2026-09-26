# AGENTS.md

Academic writing repo (research on the Ellipsoid Method / convex optimization). Sources of truth are Pandoc-Markdown files; PDFs/HTML are generated from them. There is **no application code, no test suite, and no CI** — edits are academic prose where citations and cross-references are load-bearing. Preserve them exactly.

## Documents (don't confuse them)

- `ell-review.md` — **main paper** (article). Canonical draft. Metadata in `ell-review.yaml`.
- `ell-review2.md` — parallel/alternate draft of the same paper. **Diverges in cross-ref label style** (see below).
- `ellipsoid.md` — beamer slide deck "(II)" → `ell.pdf` (built via `make slides`; the committed `ell.pdf` is a stale 2018 deck and does not match this source).
- `cutting_plane.md` — beamer slide deck "(I)" → `cutting_plane.pdf`.
- `MLE.md` / `MLE.tex` — unrelated older paper (intra-die spatial correlation, IEICE class). `MLE.tex` cannot build here: it needs `ieice.cls`, which is not in the repo.
- `try.md`, `survey.md`, `fir.md`, `min_cycle_ratio.md`, `optimal_scaling.md`, `multiplierless.md`, `profit-max.md`, `notes_ai.md`, `ell_review_ai.md`, `ell-review-cn.md` — notes / Q&A / AI summaries / a Chinese translation. Not build targets.
- `Makefile` — **the build entry point**: `make paper` / `make html` / `make slides` / `make clean`.
- `note.md` — raw-command cheat-sheet (the exact pandoc invocations).
- `ellipsoid.files/` — figures. The paper references the `.pdf` versions (`![...](ellipsoid.files/...pdf){width="80%"}`); `.svg` twins exist for HTML. `refs/*.pdf` are source references, not build outputs.

## Build

Requires pandoc + pandoc-crossref + a LaTeX engine (MiKTeX/TeX Live). Prefer the `Makefile`:

```powershell
make paper     # -> ell-review.pdf
make html      # -> ell-review.html  (needs a local katex/ directory)
make slides    # -> cutting_plane.pdf, ell.pdf  (uses xelatex)
make clean
```

The paper command it runs is:

```powershell
pandoc -F pandoc-crossref --lua-filter=secspacing.lua --citeproc -s -t latex -N --reference-links --shift-heading-level-by=-1 --csl=applied-mathematics-letters.csl ell-review.yaml latex.yaml crossref.yaml ell-review.md -o ell-review.pdf
```

Why each non-obvious flag (each fixed a real breakage — see issue #4):

- `--shift-heading-level-by=-1`: body headings start at `##` (the title is the H1). Without it pandoc emits no `\section` and numbering degrades to `0.1`, `0.2`, …. Do not remove it unless the headings are re-leveled.
- `crossref.yaml` uses `cref:false`: `siamltex.cls` redefines `\label`/`\refstepcounter`, which breaks cleveref. `cref:true` renders **every** cross-reference as `??`. Don't flip it back without patching the class.
- `secspacing.lua` fixes the cosmetic side effect of `cref:false`: pandoc-crossref joins the section symbol and the number with a non-breaking space (`§~\ref{…}`), which prints as `§ 4.2`. The filter runs **after** `-F pandoc-crossref` and rewrites the emitted `Str "§~"` to `§\ref{…}`, so references read `§4.2`. Order matters — keep `--lua-filter=secspacing.lua` after `-F pandoc-crossref`.
- `ell-review.md` references the existing `.pdf` figures, so no SVG converter (`rsvg-convert`/`inkscape`) is needed for the paper. The slide decks contain emoji and therefore build with `xelatex` (emoji render as missing glyphs until a font is configured).
- `applied-mathematics-letters.csl` is **self-contained** (its Elsevier parent style is inlined), so `--citeproc` works offline. Don't reintroduce `rel="independent-parent"`.

The `*.yaml` files are **pandoc metadata**, not app config: `latex.yaml`/`beamer.yaml` (class options), `crossref.yaml` (pandoc-crossref), plus per-doc `ell-review.yaml`. HTML builds need a local `katex/` directory (gitignored, not shipped).

## Writing conventions (enforced by tooling)

- Pandoc-Markdown: `$...$` / `$$...$$` math, citations `[@key]` / `@key` / `[@k1; @k2]`, cross-refs `@sec:...`, heading labels `{#sec:cutting_plane}`.
- **Heading-label separator is file-specific**: `ell-review.md` uses underscores (`#sec:cutting_plane`), `ell-review2.md` uses hyphens (`#sec:cutting-plane`). Match the file you edit — a mixed separator silently breaks `@sec:` links.
- Each document declares its own `bibliography:` list in YAML front matter. Cite only keys present in a bib already listed (or update the list). Many `*.bib` files exist.
- The trailing `## References {-}` heading is intentional (pandoc fills it).
- `ell-review-cn.md` content is Chinese; keep it.

## Lint / verification

No tests. Two markdown linters are configured and they disagree: `.markdownlint.json` (line length 500) and `rumdl.toml` (line length 1200, currently untracked). Recent history applies prettier, so keep formatting consistent. To verify an edit, run `make paper` (or the matching `note.md` command) and check the output.

## Gotchas

- PDFs, `.docx`, `.pptx`, and some `*.tex` are committed build artifacts. Don't hand-edit generated files.
- Stale generated outputs `ell-review.tex`, `main-diff.tex`, `main-diff.pdf` are **untracked + gitignored** (removed from the index). Don't re-add them.
- Editor backups (`*~`, `*.un~`) are gitignored; `ell-review.md~` / `.ell-review.md.un~` are untracked but still on disk.
- `envconfig.sh` is legacy Gitpod/conda setup (references a nonexistent `Config.py`, Python 3.6); `.gitpod.yml` only `chmod +x`es it — it is not a build script.
- `temp.*` is gitignored yet `temp4.docx`/`temp4.pptx` are tracked (the glob doesn't match them).

## Git

- Default branch `master`; remote branches `revision`, `imgbot`, `luk036/gitpod-setup`.
- Work often starts dirty: many `.md` files modified and `rumdl.toml` untracked. Run `git status` before assuming a clean tree.
- Commits are small prose edits ("fix typo", "improve the paper"). Match that style; do not commit unless asked.
