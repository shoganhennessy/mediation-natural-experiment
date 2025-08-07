# mediation-natural-experiment

Senan Hogan-Hennessy, most recent update 7 August 2025.

## Causal Mediation in Natural Experiments

Natural experiments are a cornerstone of applied economics, providing settings for estimating causal effects with a compelling argument for treatment randomisation, but give little indication of the mechanisms behind causal effects.
Causal Mediation (CM) provides a framework to analyse mechanisms by identifying the average direct and indirect effects (CM effects), yet conventional CM methods require the relevant mediator is as-good-as-randomly assigned.
When people choose the mediator based on costs and benefits (whether to visit a doctor, to attend university, etc.), this assumption fails and conventional CM analyses are at risk of bias.
I propose a control function strategy that uses instrumental variation in mediator take-up costs, delivering unbiased direct and indirect effects when selection is driven by unobserved gains.
The method identifies CM effects via the marginal effect of the mediator, with parametric or semi-parametric estimation that is simple to implement in two stages.
Applying these methods to the Oregon Health Insurance Experiment reveals a substantial portion of the Medicaid lottery's effect on self-reported health and happiness flows through increased healthcare usage --- an effect that a conventional CM analysis would mistake.
This approach gives applied researchers an alternative method to estimate CM effects when an initial treatment is quasi-randomly assigned, but the mediator is not, as is common in natural experiments.

- [mediation-natural-experiment-2025.pdf](https://raw.githubusercontent.com/shoganhennessy/mediation-natural-experiment/main/mediation-natural-experiment-2025.pdf) is the latest version of the working paper.
- [presentation-2025.pdf](https://raw.githubusercontent.com/shoganhennessy/mediation-natural-experiment/main/presentation-2025.pdf) is the latest version of the associated slides.

## Replication

This folder is the replication package for the paper, with statistical simulation and data analysis of the replication package for the Oregon Health Insurance Experiment (hosted at [ICSPR, and not here](https://doi.org/10.3886/ICPSR34314.v3)) using the programming language *R*, and associated packages.

A master bash file calls all code to build analysis data from raw files, produce the analysis files, and compile the final paper in one-go.

The replicator should expect the code to run for 4-6 hours, using at maximum 8GB of RAM.


## Analysis

Folder "programs/" contains multiple analysis, all using the *R* language, to run simulations and analyse  data.

- "simulations/" contains *R* scripts that simulate data, then apply new CM methods to estimate CM effects in spite of extreme unobserved selection-into-mediator (i.e., the Roy model).

- "data-build/" contains *R* scripts that extract ICSPR data from a raw folder; a replicator can run these files as is after downloading the [Oregon Health Insurance Experiment replication package](https://doi.org/10.3886/ICPSR34314.v3) for ICSPR.

- "data-analyse/" contains *R* scripts that statistically analyse, producing visualisations and tables presented in the paper.

## Text

Folder "text/" contains all files regarding the paper, written in LaTeX.

## Presentation

Folder "presentation/" contains all files regarding the presentation, written in Beamer.
