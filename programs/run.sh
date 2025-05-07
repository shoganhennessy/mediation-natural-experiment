#!/bin/bash
## Senan Hogan-Hennessy, 11 February 2025
## Master bash script, for all simulation and data evidence for the paper:
## Causal Mediation in Natural Experiments.

# Note:
# "R CMD BATCH" runs an Rscript, and logs output in a corresponding *.Rout file.



## Enact simulation evidence, and figures.
# May need to run on the COrnell BioHPC/econ cluster node.
# ssh ecco01@XXX.biohpc.cornell.edu
cd simulations
# Run the Heckman correction simulation, and plot its results
R CMD BATCH --no-save heckit-sim.R
R CMD BATCH --no-save heckit-plot.R
# Run the semi-parametric simulation, and plot its results
R CMD BATCH --no-save semiparametric-sim.R
R CMD BATCH --no-save semiparametric-plot.R
# Go back to base folder.
cd ..

## Compile the paper.
cd text
# Adjust DOIs in each bib entry (by my own custom tool)
# point to the cleaned bib file in my working paper.
cd sections
# ~/venv/bin/python3 ~/Dropbox/latex-templates/bib-edit.py 07-bibliography.bib
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

## Compile the presentation, and declare a version of the (30min) presentation.
cd presentation
# Remove the pauses for display version
presOld="dvipsnames"
presNew="dvipsnames,handout"
sed -i -e "s/$presOld/$presNew/g" presentation.tex
latexmk -pdf presentation.tex
latexmk -c
cp presentation.pdf ../presentation-2025.pdf
# Put it back to how it was.
sed -i -e "s/$presNew/$presOld/g" presentation.tex
