#!/bin/bash
## Senan Hogan-Hennessy, 11 February 2025
## Master bash script, for all simulation and data evidence for the paper:
## Causal Mediation in Natural Experiments.

# Note:
# "R CMD BATCH" runs an Rscript, and logs output in a corresponding *.Rout file.

## Enact simulation evidence, and figures.
cd simulations
# Run the first one.
R CMD BATCH --no-save roy-sim.R
# Present the results
R CMD BATCH --no-save plot-sim.R
# Go back to base folder.
cd ..

## Compile the paper.

# First append DOIs to each bib entry, by my custom tool,
# and point to the cleaned version in my working paper.
cd text/sections
~/venv/bin/python3 ~/Dropbox/latex-templates/bib-edit.py 07-bibliography.bib
bibOld="bibliography.bib"
bibNew="bibliography-doi.bib"
sed -i -e "s/$bibOld/$bibNew/g" ../paper.tex
cd ..

# Push TeX files through latexmk to get a pdf, then clean the intermed files.
latexmk -pdf paper.tex
latexmk -c
# Take the work-in-progress version, and declare as the base version.
cp paper.pdf ../mediation-natural-experiment-2025.pdf
# Put the bib file back to how it was.
sed -i -e "s/$bibNew/$bibOld/g" paper.tex
cd ..
