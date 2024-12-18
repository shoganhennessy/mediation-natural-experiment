#!/bin/bash
## Senan Hogan-Hennessy, 25 September 2023
## Master bash script, for scripts showing IV without exclusion estimator.
# Note "R CMD BATCH" generates .Rout files showing consol output of each script.

# Folder for simulation evidence
cd simulations
# Run the first one.
R CMD BATCH --no-save mediation-simulation.R
# Go back to base folder.
cd ..

# Build data sets of relevance
cd data-build
# Build a panel of HRS data
R CMD BATCH --no-save hrs-panel.R
# Combine panel of HRS data with genetic educ scores
R CMD BATCH --no-save hrs-build.R
# Collapse the panel of HRS data, keeping genetic educ scores
R CMD BATCH --no-save hrs-collapse.R
# Go back to base folder.
cd ..

# Statistical analysis code
cd data-analyse
# Summarise HRS data.
R CMD BATCH --no-save hrs-summarise.R
# Analsyse HRS data, with the invalid instrument educ score.
R CMD BATCH --no-save hrs-educ-iv.R
# Go back to base folder.
cd ..

# Compile paper, with outputs from R scripts above as inputs for TeX files
cd ../text
latexmk -pdf paper.tex
latexmk -c
cp paper.pdf ../mr-education.pdf
