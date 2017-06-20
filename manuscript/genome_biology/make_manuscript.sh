#!/usr/bin/env sh
pdflatex -shell-escape bmc_article.tex;
bibtex bmc_article;
pdflatex -shell-escape bmc_article.tex;
pdflatex -shell-escape bmc_article.tex;
mv bmc_article.log bmc_article.aux bmc_article.pdf bmc_article.bbl bmc_article.blg compile_dir;
