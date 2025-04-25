#!/usr/bin/R
## Senan Hogan-Hennessy, 20 Jan 2025
## Plot the system for indirect + direct effects, with Heckman correction.
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
colour.list <- c("orange", "blue")


################################################################################
## Load the simulated data.  (saved as separate file in advance.)

# Load data from DGP-strapping 10,000 times.
sim.data <- read_csv(file.path(output.folder, "dist-heckit-data.csv"))
print(sim.data)


################################################################################
## Plot the dist of the ADE and AIE estimates, around truth.

# ADE estimates, by type.
direct_dist.plot <- sim.data %>%
    ggplot() +
    # Dist of OLS estimates.
    geom_density(aes(x = ols_direct_effect, y = after_stat(density)),
        colour = "black", fill = colour.list[1], alpha = 0.75) +
    annotate("text", colour = colour.list[1],
        x = 0.25, y = 3,
        fontface = "bold",
        label = ("Unadjusted"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = colour.list[1],
        x = 0.35, y = 2.75,
        xend = 0.75, yend = 2.5,
        linewidth = 0.75, curvature = 0.25,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    # Dist of CF estimates.
    geom_density(aes(x = cf_direct_effect, y = after_stat(density)),
        colour = "black", fill = colour.list[2], alpha = 0.75) +
    annotate("text", colour = colour.list[2],
        x = 2, y = 2,
        fontface = "bold",
        label = ("Selection model"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = colour.list[2],
        x = 1.875, y = 1.75,
        xend = 1.5, yend = 1.25,
        linewidth = 0.75, curvature = -0.25,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    # Truth value
    geom_vline(xintercept = mean(sim.data$truth_direct_effect),
        colour = "black", linetype = "dashed", linewidth = 1) +
    annotate("text", colour = "black",
        x = 2, y = 3.5,
        fontface = "bold",
        label = ("Truth"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = "black",
        x = 1.875, y = 3.5,
        xend = 1.625, yend = 3.25,
        linewidth = 0.75,
        curvature = -0.25,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    # Other presentation options
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Estimate",
        breaks = seq(0, 2.5, by = 0.25),
        limits = c(0, 2.5)) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        limits = c(0, 4.5)) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(output.folder, "heckit-direct-dist.png"),
    plot = direct_dist.plot,
    units = "cm", width = fig.width, height = fig.height)

# AIE estimates, by type.
indirect_dist.plot <- sim.data %>%
    ggplot() +
    # Dist of OLS estimates.
    geom_density(aes(x = ols_indirect_effect, y = after_stat(density)),
        colour = "black", fill = colour.list[1], alpha = 0.75) +
    # Dist of CF estimates.
    geom_density(aes(x = cf_indirect_effect, y = after_stat(density)),
        colour = "black", fill = colour.list[2], alpha = 0.75) +
    # Truth value
    geom_vline(xintercept = mean(sim.data$truth_indirect_effect),
        colour = "black", linetype = "dashed", linewidth = 1) +
    # Other presentation options
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = "Estimate",
        breaks = seq(0, 2.5, by = 0.25),
        limits = c(0, 2.5)) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        limits = c(0, 4.5)) +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(output.folder, "heckit-indirect-dist.png"),
    plot = indirect_dist.plot,
    units = "cm", width = fig.width, height = fig.height)
