# 👉 Note

# Preferred: use the Makefile (see Makefile).
#   make paper   -> ell-review.pdf
#   make html    -> ell-review.html   (requires a local katex/ directory)
#   make slides  -> cutting_plane.pdf ell.pdf   (uses xelatex)
#   make clean
#
# Raw commands below, for reference.
#
# Main paper -> ell-review.pdf.
# crossref.yaml uses cref:false: siamltex.cls redefines \label/\refstepcounter,
#   which breaks cleveref (cref:true renders every cross-reference as "??").
# --shift-heading-level-by=-1 : body headings in ell-review.md start at "##".
# --lua-filter=secspacing.lua : drop the non-breaking space pandoc-crossref puts
#   in section references, so they print as "§4.2" instead of "§ 4.2".
pandoc -F pandoc-crossref --lua-filter=secspacing.lua --citeproc -s -t latex -N --reference-links --shift-heading-level-by=-1 --csl=siam-numeric.csl ell-review.yaml latex.yaml crossref.yaml ell-review.md -o ell-review.pdf

pandoc -F pandoc-crossref --citeproc -s -t html -N --katex=katex/ --reference-links --csl=siam-numeric.csl ell-review.yaml latex.yaml crossref.yaml ell-review.md -o ell-review.html

pandoc -F pandoc-crossref -s -t beamer --toc --natbib --reference-links --pdf-engine=xelatex --csl=siam-numeric.csl beamer.yaml cutting_plane.md -o cutting_plane.pdf

pandoc -F pandoc-crossref -s -t html --katex=katex/ --toc --natbib --reference-links --csl=siam-numeric.csl beamer.yaml cutting_plane.md -o cutting_plane.html

pandoc -s --wrap=preserve ell-review.md -o temp.md

pandoc -s -t beamer --natbib --toc --pdf-engine=xelatex beamer.yaml ellipsoid.md -o ell.pdf

pandoc -F pandoc-crossref --citeproc -s -t latex -N latex.yaml crossref.yaml ell-review.md -o temp.tex

pandoc -F pandoc-crossref --citeproc -s -t html -N --katex=katex/ crossref.yaml ell-review.md -o temp.html

\usepackage{subfig}
\AtBeginDocument{%
\renewcommand*\figurename{Figure}
\renewcommand*\tablename{Table}
}
\AtBeginDocument{%
\renewcommand*\listfigurename{\#\#
List of Figures}
\renewcommand*\listtablename{\#\#
List of Tables}
}
\usepackage{float}
\floatstyle{ruled}
\makeatletter
\@ifundefined{c@chapter}{\newfloat{codelisting}{h}{lop}}{\newfloat{codelisting}{h}{lop}{[}chapter{]}}
\makeatother
\floatname{codelisting}{Listing}
\newcommand\*\listoflistings{\listof{codelisting}{List
of Listings}}
\usepackage{cleveref}
\crefname{figure}{Fig.}{Fig.}
\crefname{table}{Table}{Table}
\crefname{equation}{Eq.}{Eq.}
\crefname{listing}{Listing}{Listing}
\crefname{section}{§}{§}
\Crefname{figure}{Fig.}{Fig.}
\Crefname{table}{Table}{Table}
\Crefname{equation}{Eq.}{Eq.}
\Crefname{listing}{Listing}{Listing}
\Crefname{section}{§}{§}
\makeatletter
\crefname{codelisting}{\cref@listing@name}{\cref@listing@name@plural}
\Crefname{codelisting}{\Cref@listing@name}{\Cref@listing@name@plural}
\makeatother
