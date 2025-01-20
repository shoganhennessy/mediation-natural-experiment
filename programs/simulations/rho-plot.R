#!/usr/bin/R
## Senan Hogan-Hennessy, 20 Jan 2025
## Plot the system for indirect + direct effects, with Roy selection.
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

# Load the data file for varying correlated errors.
rho.data <- read_csv(file.path(output.folder, "rho-sim-data.csv"))
print(names(rho.data))


################################################################################
## Plot the Direct effect estimates, by OLS + CF

# Plot the bias in direct effect est vs rho, with different sigma values
directeffect_bias.plot <- rho.data %>%
    ggplot(aes(x = rho)) +
    # OLS est + 95 % CI
    geom_point(aes(y = ols_direct_effect), colour = colour.list[1]) +
    geom_ribbon(aes(ymin = ols_direct_effect_low, ymax = ols_direct_effect_up),
        fill = colour.list[1], alpha = 0.2) +
    annotate("text", colour = colour.list[1],
        x = -0.5, y = 1.65,
        fontface = "bold",
        label = ("OLS"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = colour.list[1],
        x = -0.4, y = 1.75,
        xend = -0.25, yend = 1.9,
        linewidth = 0.75,
        curvature = 0.25,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    # CF est + 95 % CI
    geom_point(aes(y = cf_direct_effect), colour = colour.list[2]) +
    geom_ribbon(aes(ymin = cf_direct_effect_low, ymax = cf_direct_effect_up),
        fill = colour.list[2], alpha = 0.2) +
    annotate("text", colour = colour.list[2],
        x = -0.5, y = 1,
        fontface = "bold",
        label = ("Control function"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = colour.list[2],
        x = -0.5, y = 1.125,
        xend = -0.375, yend = 1.3,
        linewidth = 0.75,
        curvature = 0.25,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    # Truth:
    geom_hline(aes(yintercept = mean(truth_direct_effect)),
        colour = "black", linetype = "dashed", linewidth = 0.75) +
    annotate("text", colour = "black",
        x = 0.65, y = 0.75,
        fontface = "bold",
        label = ("Truth"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = "black",
        x = 0.6, y = 0.875,
        xend = 0.4, yend = 1.4,
        linewidth = 0.75,
        curvature = -0.25,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    # Presentation options
    theme_bw() +
    scale_x_continuous(
        name = TeX("$\\rho$"),
        expand = c(0, 0),
        breaks = seq(-1, 1, by = 0.2),
        limits = c(-1.05, 1.05)) +
    scale_y_continuous(name = "",
        breaks = seq(0, 2.25, by = 0.25),
        limits = c(0, 2.25),
        expand = c(0.01, 0.01)) +
    ggtitle("Estimate") +
    theme(plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0.25, 0.25), "mm"))
# Save this plot
ggsave(file.path(output.folder, "rho-directeffect-bias.png"),
    plot = directeffect_bias.plot, dpi = 300,
    units = "cm", width = fig.width, height = fig.height)


################################################################################
## Plot the Indirect effect estimates, by OLS + CF.

# Plot the bias in indirect effect est vs rho, with different sigma values
indirecteffect_bias.plot <- rho.data %>%
    ggplot(aes(x = rho)) +
    # OLS est + 95 % CI
    geom_point(aes(y = ols_indirect_effect), colour = colour.list[1]) +
    geom_ribbon(aes(ymin = ols_indirect_effect_low, ymax = ols_indirect_effect_up),
        fill = colour.list[1], alpha = 0.2) +
    # CF est + 95 % CI
    geom_point(aes(y = cf_indirect_effect), colour = colour.list[2]) +
    geom_ribbon(aes(ymin = cf_indirect_effect_low, ymax = cf_indirect_effect_up),
        fill = colour.list[2], alpha = 0.2) +
    # Truth:
    geom_hline(aes(yintercept = mean(truth_indirect_effect)),
        colour = "black", linetype = "dashed", linewidth = 0.75) +
    # Presentation options
    theme_bw() +
    scale_x_continuous(
        name = TeX("$\\rho$"),
        expand = c(0, 0),
        breaks = seq(-1, 1, by = 0.2),
        limits = c(-1.05, 1.05)) +
    scale_y_continuous(name = "",
        breaks = seq(0, 2.25, by = 0.25),
        limits = c(0, 2.25),
        expand = c(0.01, 0.01)) +
    ggtitle("Estimate") +
    theme(plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0.25, 0.25), "mm"),
        legend.position = "bottom")
# Save this plot
ggsave(file.path(output.folder, "rho-indirecteffect-bias.png"),
    plot = indirecteffect_bias.plot, dpi = 300,
    units = "cm", width = fig.width, height = fig.height)
