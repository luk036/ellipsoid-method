# Build entry point for the ellipsoid-method documents.
#
# Requires: pandoc, pandoc-crossref and a LaTeX engine (MiKTeX/TeX Live) with pdflatex.
# Citations use applied-mathematics-letters.csl, which is self-contained (no network).
#
# Policy notes (see AGENTS.md / issue #4):
#   * crossref.yaml uses cref:false because siamltex.cls overrides \label/\refstepcounter
#     and breaks cleveref (cref:true would render every cross-reference as "??").
#   * --shift-heading-level-by=-1 because the body headings start at "##" (title is the H1);
#     without it pandoc emits no \section and numbering degrades to 0.1, 0.2, ...
#   * ell-review.md references the existing .pdf figures, so no SVG converter is needed.
#   * secspacing.lua removes the non-breaking space pandoc-crossref puts in section
#     references ("§~\ref{...}"), so "§4.2" is printed rather than "§ 4.2".

PANDOC   := pandoc
CROSSREF := pandoc-crossref
CSL      := applied-mathematics-letters.csl

PAPER_FLAGS := -F $(CROSSREF) --lua-filter=secspacing.lua --citeproc -s -t latex -N --reference-links \
               --shift-heading-level-by=-1 --csl=$(CSL)
PAPER_META  := ell-review.yaml latex.yaml crossref.yaml

.PHONY: all paper multiplierless html slides clean

all: paper

# --- Main paper (SIAM article) ---------------------------------------------
paper: ell-review.pdf

ell-review.pdf: ell-review.md $(PAPER_META) $(CSL) secspacing.lua
	$(PANDOC) $(PAPER_FLAGS) $(PAPER_META) ell-review.md -o $@

# --- Multiplierless FIR paper ----------------------------------------------
# Shares latex.yaml and crossref.yaml with the main paper; same PAPER_FLAGS.
# Only the per-document metadata (multiplierless.yaml) differs.
multiplierless: multiplierless.pdf

multiplierless.pdf: multiplierless.md multiplierless.yaml latex.yaml crossref.yaml $(CSL) secspacing.lua
	$(PANDOC) $(PAPER_FLAGS) multiplierless.yaml latex.yaml crossref.yaml multiplierless.md -o $@

# Requires a local katex/ directory (gitignored, not shipped).
html: ell-review.html

ell-review.html: ell-review.md $(PAPER_META) $(CSL)
	$(PANDOC) -F $(CROSSREF) --citeproc -s -t html -N --katex=katex/ \
	          --reference-links --csl=$(CSL) $(PAPER_META) ell-review.md -o $@

# --- Beamer slides ----------------------------------------------------------
# Slides contain emoji, so they need xelatex (pdflatex rejects non-ASCII).
# Remaining "Missing character" warnings for emoji are expected until the
# source cleanup tracked in issue #4 (P2).
slides: cutting_plane.pdf ell.pdf

cutting_plane.pdf: cutting_plane.md beamer.yaml
	$(PANDOC) -F $(CROSSREF) -s -t beamer --toc --natbib --reference-links \
	          --pdf-engine=xelatex --csl=$(CSL) beamer.yaml cutting_plane.md -o $@

ell.pdf: ellipsoid.md beamer.yaml
	$(PANDOC) -s -t beamer --natbib --toc --pdf-engine=xelatex beamer.yaml ellipsoid.md -o $@

# --- Housekeeping -----------------------------------------------------------
# The trailing `; true` forces make to run these through a shell; without it make
# may exec `rm` directly, which fails on Windows (`rm` is not an .exe there).
clean:
	-@rm -f temp.tex temp.pdf temp.md temp.html temp.txt ; true
	-@rm -f *.aux *.bbl *.blg *.log *.out *.nav *.snm *.vrb *.toc ; true
