#!/usr/bin/R
## Senan Hogan-Hennessy, 20 Jan 2025
## Plot the estimates for indirect + direct effects, with Heckman correction.
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
print(head(sim.data))
#View(sim.data)


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
        label = ("Parametric CF"),
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
ggsave(file.path(output.folder, "heckit-direct-dist.png"),
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
ggsave(file.path(output.folder, "heckit-indirect-dist.png"),
    plot = indirect_dist.plot,
    units = "cm", width = fig.width, height = fig.height)


#! Below is unfinished.
quit("no")
################################################################################
## Plot the CF effect estimates, by OLS + CF,
## different $\sigma = \sigma_1 / \sigma_0$ values.

# Plot the bias in direct effect est vs sigma
sigma_directeffect_bias.plot <- sigma.data %>%
    ggplot(aes(x = sigma)) +
    # OLS est + 95 % CI
    geom_point(aes(y = (ols_direct_effect - truth_direct_effect)),
        colour = colour.list[1]) +
    geom_ribbon(aes(
        ymin = (ols_direct_effect_low - truth_direct_effect),
        ymax = (ols_direct_effect_up - truth_direct_effect)),
        fill = colour.list[1], alpha = 0.2) +
    # CF est + 95 % CI
    geom_point(aes(y = (cf_direct_effect - truth_direct_effect)),
        colour = colour.list[2]) +
    geom_ribbon(aes(
        ymin = (cf_direct_effect_low - truth_direct_effect),
        ymax = (cf_direct_effect_up - truth_direct_effect)),
        fill = colour.list[2], alpha = 0.2) +
    # Truth:
    geom_line(aes(y = (truth_direct_effect - truth_direct_effect)),
        colour = "black", linetype = "dashed", linewidth = 1) +
    # Presentation options
    theme_bw() +
    scale_x_continuous(
        name = TeX("$\\sigma$"),
        expand = c(0, 0),
        breaks = seq(0, 2, by = 0.25),
        limits = c(-0.02, 2.02)) +
    scale_y_continuous(name = "",
        breaks = seq(-1.5, 1.5, by = 0.25),
        limits = c(-1.5, 1.5),
        expand = c(0.01, 0.01)) +
    ggtitle("Estimate") +
    theme(plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 1.5, 0.25), "mm"),
        axis.title.x = element_text(vjust = -0.25))
# Save this plot
ggsave(file.path(output.folder, "sigma-directeffect-bias.png"),
    plot = sigma_directeffect_bias.plot, dpi = 300,
    units = "cm", width = fig.width, height = fig.height)

# Plot the bias in indirect effect est vs sigma
sigma_indirecteffect_bias.plot <- sigma.data %>%
    ggplot(aes(x = sigma)) +
    # OLS est + 95 % CI
    geom_point(aes(y = ols_indirect_effect - truth_indirect_effect),
        colour = colour.list[1]) +
    geom_ribbon(aes(
        ymin = (ols_indirect_effect_low - truth_indirect_effect),
        ymax = (ols_indirect_effect_up - truth_indirect_effect)),
        fill = colour.list[1], alpha = 0.2) +
    # CF est + 95 % CI
    geom_point(aes(y = cf_indirect_effect- truth_indirect_effect),
        colour = colour.list[2]) +
    geom_ribbon(aes(
        ymin = (cf_indirect_effect_low - truth_indirect_effect),
        ymax = (cf_indirect_effect_up - truth_indirect_effect)),
        fill = colour.list[2], alpha = 0.2) +
    # Truth:
    geom_line(aes(y = truth_indirect_effect - truth_indirect_effect),
        colour = "black", linetype = "dashed", linewidth = 1) +
    # Presentation options
    theme_bw() +
    scale_x_continuous(
        name = TeX("$\\sigma$"),
        expand = c(0, 0),
        breaks = seq(0, 2, by = 0.25),
        limits = c(-0.02, 2.02)) +
    scale_y_continuous(name = "",
        breaks = seq(-1.5, 1.5, by = 0.25),
        limits = c(-1.5, 1.5),
        expand = c(0.01, 0.01)) +
    ggtitle("Estimate") +
    theme(plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 1.5, 0.25), "mm"),
        axis.title.x = element_text(vjust = -0.25))
# Save this plot
ggsave(file.path(output.folder, "sigma-indirecteffect-bias.png"),
    plot = sigma_indirecteffect_bias.plot, dpi = 300,
    units = "cm", width = fig.width, height = fig.height)


################################################################################
## Plot the Direct effect estimates, by OLS + CF, different $\rho$ values.

# Plot the bias in direct effect est vs rho
rho_directeffect_bias.plot <- rho.data %>%
    ggplot(aes(x = rho)) +
    # OLS est + 95 % CI
    geom_point(aes(y = ols_direct_effect), colour = colour.list[1]) +
    geom_ribbon(aes(ymin = ols_direct_effect_low, ymax = ols_direct_effect_up),
        fill = colour.list[1], alpha = 0.2) +
    annotate("text", colour = colour.list[1],
        x = 0.15, y = 0.2,
        fontface = "bold",
        label = ("OLS"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = colour.list[1],
        x = 0.3, y = 0.25,
        xend = 0.5, yend = 0.5,
        linewidth = 0.75,
        curvature = 0.25,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    # CF est + 95 % CI
    geom_point(aes(y = cf_direct_effect), colour = colour.list[2]) +
    geom_ribbon(aes(ymin = cf_direct_effect_low, ymax = cf_direct_effect_up),
        fill = colour.list[2], alpha = 0.2) +
    annotate("text", colour = colour.list[2],
        x = -0.5, y = 2.00,
        fontface = "bold",
        label = ("Parametric CF"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = colour.list[2],
        x = -0.375, y = 1.95,
        xend = -0.25, yend = 1.55,
        linewidth = 0.75,
        curvature = -0.125,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    # Truth:
    geom_line(aes(y = (truth_direct_effect)),
        colour = "black", linetype = "dashed", linewidth = 1) +
    annotate("text", colour = "black",
        x = 0.65, y = 1.8,
        fontface = "bold",
        label = ("Truth"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = "black",
        x = 0.6, y = 1.75,
        xend = 0.475, yend = 1.4875,
        linewidth = 0.75, curvature = 0.125,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    # Presentation options
    theme_bw() +
    scale_x_continuous(
        name = TeX("$\\rho$"),
        expand = c(0, 0),
        breaks = seq(-1, 1, by = 0.25),
        limits = c(-1.025, 1.025)) +
    scale_y_continuous(name = "",
        breaks = seq(0, 2.5, by = 0.25),
        limits = c(0, 2.5),
        expand = c(0.01, 0.01)) +
    ggtitle("Estimate") +
    theme(plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 1.5, 0.25), "mm"),
        axis.title.x = element_text(vjust = -0.25))
# Save this plot
ggsave(file.path(output.folder, "rho-directeffect-bias.png"),
    plot = rho_directeffect_bias.plot, dpi = 300,
    units = "cm", width = fig.width, height = fig.height)

# Plot the bias in indirect effect est vs rho
rho_indirecteffect_bias.plot <- rho.data %>%
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
    geom_line(aes(y = (truth_indirect_effect)),
        colour = "black", linetype = "dashed", linewidth = 1) +
    # Presentation options
    theme_bw() +
    scale_x_continuous(
        name = TeX("$\\rho$"),
        expand = c(0, 0),
        breaks = seq(-1, 1, by = 0.25),
        limits = c(-1.025, 1.025)) +
    scale_y_continuous(name = "",
        breaks = seq(0, 2.5, by = 0.25),
        limits = c(0, 2.5),
        expand = c(0.01, 0.01)) +
    ggtitle("Estimate") +
    theme(plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 1.5, 0.25), "mm"),
        axis.title.x = element_text(vjust = -0.25))
# Save this plot
ggsave(file.path(output.folder, "rho-indirecteffect-bias.png"),
    plot = rho_indirecteffect_bias.plot, dpi = 300,
    units = "cm", width = fig.width, height = fig.height)


################################################################################
## Plot the CF effect estimates, by OLS + CF, different $\sigma_1$ values.

# Plot the bias in direct effect est vs sigma_1
sigma_1_directeffect_bias.plot <- sigma_1.data %>%
    ggplot(aes(x = sigma_1)) +
    # OLS est + 95 % CI
    geom_point(aes(y = ols_direct_effect), colour = colour.list[1]) +
    geom_ribbon(aes(ymin = ols_direct_effect_low,
        ymax = ols_direct_effect_up),
        fill = colour.list[1], alpha = 0.2) +
    # CF est + 95 % CI
    geom_point(aes(y = cf_direct_effect), colour = colour.list[2]) +
    geom_ribbon(aes(ymin = cf_direct_effect_low, ymax = cf_direct_effect_up),
        fill = colour.list[2], alpha = 0.2) +
    # Truth:
    geom_line(aes(y = truth_direct_effect),
        colour = "black", linetype = "dashed", linewidth = 1) +
    # Presentation options
    theme_bw() +
    scale_x_continuous(
        name = TeX("$\\sigma_1$"),
        expand = c(0, 0),
        breaks = seq(0, 2, by = 0.25),
        limits = c(-0.02, 2.02)) +
    scale_y_continuous(name = "",
        breaks = seq(0, 2.5, by = 0.25),
        limits = c(0, 2.5),
        expand = c(0.01, 0.01)) +
    ggtitle("Estimate") +
    theme(plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 1.5, 0.25), "mm"),
        axis.title.x = element_text(vjust = -0.25))
# Save this plot
ggsave(file.path(output.folder, "sigma-1-directeffect-bias.png"),
    plot = sigma_1_directeffect_bias.plot, dpi = 300,
    units = "cm", width = fig.width, height = fig.height)


# Plot the bias in indirect effect est vs sigma_1
sigma_1_indirecteffect_bias.plot <- sigma_1.data %>%
    ggplot(aes(x = sigma_1)) +
    # OLS est + 95 % CI
    geom_point(aes(y = ols_indirect_effect), colour = colour.list[1]) +
    geom_ribbon(aes(ymin = ols_indirect_effect_low,
        ymax = ols_indirect_effect_up),
        fill = colour.list[1], alpha = 0.2) +
    # CF est + 95 % CI
    geom_point(aes(y = cf_indirect_effect), colour = colour.list[2]) +
    geom_ribbon(aes(ymin = cf_indirect_effect_low, ymax = cf_indirect_effect_up),
        fill = colour.list[2], alpha = 0.2) +
    # Truth:
    geom_line(aes(y = truth_indirect_effect),
        colour = "black", linetype = "dashed", linewidth = 1) +
    # Presentation options
    theme_bw() +
    scale_x_continuous(
        name = TeX("$\\sigma_1$"),
        expand = c(0, 0),
        breaks = seq(0, 2, by = 0.25),
        limits = c(-0.02, 2.02)) +
    scale_y_continuous(name = "",
        breaks = seq(0, 2.5, by = 0.25),
        limits = c(0, 2.5),
        expand = c(0.01, 0.01)) +
    ggtitle("Estimate") +
    theme(plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 1.5, 0.25), "mm"),
        axis.title.x = element_text(vjust = -0.25))
# Save this plot
ggsave(file.path(output.folder, "sigma-1-indirecteffect-bias.png"),
    plot = sigma_1_indirecteffect_bias.plot, dpi = 300,
    units = "cm", width = fig.width, height = fig.height)
