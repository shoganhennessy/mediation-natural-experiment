#!/bin/bash
## Senan Hogan-Hennessy, 25 September 2023
## Master bash script, for scripts showing IV without exclusion estimator.
# Note "R CMD BATCH" generates .Rout files showing consol output of each script.

# Folder for simulation evidence
cd simulations
# Run the first one.
R CMD BATCH --no-save roy-sim.R
# Present the results
R CMD BATCH --no-save plot-sim.R
# Go back to base folder.
cd ..

# Compile paper, with outputs from R scripts above as inputs for TeX files
cd ../text
latexmk -pdf paper.tex
latexmk -c
cp paper.pdf ../mediation-natural-experiment.pdf
