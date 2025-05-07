#!/usr/bin/R
## Senan Hogan-Hennessy, 7 May 2025
## Plot the estimates for indirect + direct effects, with semi-parametric CF.
# Show the date:
print(format(Sys.time(), "%H:%M %Z %A, %d %B %Y"))

## Load libraries
# Functions for data manipulation and visualisation
library(tidyverse)
# Library for better colour choice.
library(ggthemes)
# Library for equations in plots
library(latex2exp)

## Set up the R environment
set.seed(47)
# Define number of digits in tables and graphs
digits.no <- 3
# Define where output files go.
output.folder <- file.path("sim-output")
# Set the options for the plot sizes, in saving ggplot output.
fig.height <- 10
fig.width <- fig.height
# Define colours to use.
colour.list <- c("orange", "#2ca02c") # Green


################################################################################
## Load the simulated data.  (saved as separate file in advance.)

# Load data from DGP-strapping 10,000 times.
sim.data <- read_csv(file.path(output.folder, "dist-semiparametric-data.csv"))
print(sim.data)
#View(head(sim.data)


################################################################################
## Plot the dist of the ADE and AIE estimates, around truth.

# ADE estimates, by type.
direct_dist.plot <- sim.data %>%
    ggplot() +
    # Dist of OLS estimates.
    geom_density(aes(x = (ols_direct_effect - truth_direct_effect),
        y = after_stat(density)),
        colour = "black", fill = colour.list[1], alpha = 0.75) +
    annotate("text", colour = colour.list[1],
        x = 0.5, y = 4,
        fontface = "bold",
        label = ("Unadjusted"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = colour.list[1],
        x = 0.5, y = 3.875,
        xend = -0.25, yend = 3,
        linewidth = 0.75, curvature = -0.25,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    # Dist of CF estimates.
    geom_density(aes(x = (cf_direct_effect - truth_direct_effect),
        y = after_stat(density)),
        colour = "black", fill = colour.list[2], alpha = 0.75) +
    annotate("text", colour = colour.list[2],
        x = 0.5, y = 2.5,
        fontface = "bold",
        label = ("Semi-parametric CF"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = colour.list[2],
        x = 0.5, y = 2.375,
        xend = 0.15, yend = 1.6,
        linewidth = 0.75, curvature = -0.25,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    # Truth value
    geom_vline(xintercept = 0,
        colour = "black", linetype = "dashed", linewidth = 1) +
    # Other presentation options
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = TeX("Estimate$-$ True Value"),,
        breaks = seq(-1.0, 1.0, by = 0.25),
        limits = c(-1.0, 1.0)) +
    scale_y_continuous(expand = c(0, 0),
        name = "", limits = c(0, 5.1)) +
    ggtitle("Density") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(output.folder, "semiparametric-direct-dist.png"),
    plot = direct_dist.plot,
    units = "cm", width = fig.width, height = fig.height)

# AIE estimates, by type.
indirect_dist.plot <- sim.data %>%
    ggplot() +
    # Dist of OLS estimates.
    geom_density(aes(x = ols_indirect_effect - truth_indirect_effect,
        y = after_stat(density)),
        colour = "black", fill = colour.list[1], alpha = 0.75) +
    # Dist of CF estimates.
    geom_density(aes(x = cf_indirect_effect - truth_indirect_effect,
        y = after_stat(density)),
        colour = "black", fill = colour.list[2], alpha = 0.75) +
    # Truth value
    geom_vline(xintercept = 0,
        colour = "black", linetype = "dashed", linewidth = 1) +
    # Other presentation options
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = TeX("Estimate$-$ True Value"),
        breaks = seq(-1.0, 1.0, by = 0.25),
        limits = c(-1.0, 1.0)) +
    scale_y_continuous(expand = c(0, 0),
        name = "", limits = c(0, 5.1)) +
    ggtitle("Density") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(output.folder, "semiparametric-indirect-dist.png"),
    plot = indirect_dist.plot,
    units = "cm", width = fig.width, height = fig.height)
