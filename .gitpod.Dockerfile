FROM gitpod/workspace-full:latest

USER root
# Install util tools.
RUN apt-get update \
 && apt-get install -y \
  apt-utils \
  sudo \
  git \
  less \
  texlive \
  texlive-science \
  texlive-latex-extra \
  latexmk \
  chktex \
  latexdiff \
  ktikz \
  wget

# Give back control
USER root

# Cleaning
RUN apt-get clean
