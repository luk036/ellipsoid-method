# AGENTS.md

Academic writing repo (research on the Ellipsoid Method / convex optimization). Sources of truth are Pandoc-Markdown files; PDFs/HTML are generated from them. There is **no application code, no test suite, and no CI** — edits are academic prose where citations and cross-references are load-bearing. Preserve them exactly.

## Documents (don't confuse them)

- `ell-review.md` — **main paper** (article). Canonical draft. Metadata in `ell-review.yaml`.
- `ell-review2.md` — parallel/alternate draft of the same paper. **Diverges in cross-ref label style** (see below).
- `ellipsoid.md` — beamer slide deck "(II)" → `ell.pdf`.
- `cutting_plane.md` — beamer slide deck "(I)" → `cutting_plane.pdf`.
- `MLE.md` / `MLE.tex` — unrelated older paper (intra-die spatial correlation, IEICE class). `MLE.tex` cannot build here: it needs `ieice.cls`, which is not in the repo.
- `try.md`, `survey.md`, `fir.md`, `min_cycle_ratio.md`, `optimal_scaling.md`, `multiplierless.md`, `profit-max.md`, `notes_ai.md`, `ell_review_ai.md`, `ell-review-cn.md` — notes / Q&A / AI summaries / a Chinese translation. Not build targets.
- `note.md` — **the build cheat-sheet**. Read it before building; it holds the exact pandoc invocations.
- `ellipsoid.files/` — SVG/PDF/PNG figures referenced by markdown (`![...](ellipsoid.files/...svg){width="80%"}`). `refs/*.pdf` are source references, not build outputs.

## Build

pandoc + pandoc-crossref + MiKTeX are installed on Windows. Do **not** invent commands; use `note.md`. Representative examples:

```powershell
# Article PDF (main paper)
pandoc -F pandoc-crossref --citeproc -s -t latex -N --reference-links --csl=applied-mathematics-letters.csl ell-review.yaml latex.yaml crossref.yaml ell-review.md -o ell-review.pdf

# Beamer slides
pandoc -F pandoc-crossref -s -t beamer --toc --natbib --reference-links --csl=applied-mathematics-letters.csl beamer.yaml cutting_plane.md -o temp.tex
```

The `*.yaml` files are **pandoc metadata**, not app config: `latex.yaml`/`beamer.yaml` (class options), `crossref.yaml`/`crossref_2.yaml` (pandoc-crossref), plus per-doc `ell-review.yaml`. HTML builds use `--katex=katex/`; the `katex/` directory is gitignored and absent, so you must supply it.

## Writing conventions (enforced by tooling)

- Pandoc-Markdown: `$...$` / `$$...$$` math, citations `[@key]` / `@key` / `[@k1; @k2]`, cross-refs `@sec:...`, heading labels `{#sec:cutting_plane}`.
- **Heading-label separator is file-specific**: `ell-review.md` uses underscores (`#sec:cutting_plane`), `ell-review2.md` uses hyphens (`#sec:cutting-plane`). Match the file you edit — a mixed separator silently breaks `@sec:` links.
- Each document declares its own `bibliography:` list in YAML front matter. Cite only keys present in a bib already listed (or update the list). Many `*.bib` files exist.
- The trailing `## References {-}` heading is intentional (pandoc fills it).
- `ell-review-cn.md` content is Chinese; keep it.

## Lint / verification

No tests. Two markdown linters are configured and they disagree: `.markdownlint.json` (line length 500) and `rumdl.toml` (line length 1200, currently untracked). Recent history applies prettier, so keep formatting consistent. To verify an edit, rebuild the affected doc with the matching `note.md` command and check the output.

## Gotchas

- PDFs, `.docx`, `.pptx`, and some `*.tex` are committed build artifacts. Don't hand-edit generated files.
- `ell-review.tex` is a **stale beamer output** — it does not match the current `ell-review.md`/`ell-review.yaml`. `main-diff.tex`/`main-diff.pdf` are 2019 `latexdiff` output (reference `/tmp/...` paths) — also stale.
- Editor backups are committed and **not** gitignored: `ell-review.md~`, `.ell-review.md.un~`. Leave them alone unless asked.
- `envconfig.sh` is legacy Gitpod/conda setup (references a nonexistent `Config.py`, Python 3.6); `.gitpod.yml` only `chmod +x`es it — it is not a build script.
- `temp.*` is gitignored yet `temp4.docx`/`temp4.pptx` are tracked (the glob doesn't match them).

## Git

- Default branch `master`; remote branches `revision`, `imgbot`, `luk036/gitpod-setup`.
- Work often starts dirty: many `.md` files modified and `rumdl.toml` untracked. Run `git status` before assuming a clean tree.
- Commits are small prose edits ("fix typo", "improve the paper"). Match that style; do not commit unless asked.
