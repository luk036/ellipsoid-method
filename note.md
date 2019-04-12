pandoc -s --wrap=preserve ell-review.md -o temp.md

pandoc -F pandoc-crossref -F pandoc-citeproc -s -t latex -N latex.yaml crossref.yaml ell-review.md -o temp.tex

pandoc -F pandoc-crossref -F pandoc-citeproc -s -t html -N --katex=katex/ crossref.yaml ell-review.md -o temp.html

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
\newcommand*\listoflistings{\listof{codelisting}{List
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
