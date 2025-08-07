#!/usr/bin/R
## Senan Hogan-Hennessy, 7 May 2025
## Plot the estimates for indirect + direct effects, with semiparametric CF.
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
# Define where input files come form + output files go.
input.folder <- file.path("sim-output")
output.folder <- file.path("..", "..", "text", "sections", "figures")
# Set the options for the plot sizes, in saving ggplot output.
fig.height <- 10
fig.width <- fig.height
# Define colours to use.
colour.list <- c("orange", "blue", "#2ca02c") # Orange, blue, Green


################################################################################
## Load the simulated data.  (saved as separate file in advance.)

# Load data from normal DGP-strapping 10,000 times.
normal.data <- read_csv(file.path(input.folder, "normal-cf-data.csv"))
print(normal.data)

# Load data from uniform DGP-strapping 10,000 times.
uniform.data <- read_csv(file.path(input.folder, "uniform-cf-data.csv"))
print(uniform.data)

# Load data from different Corr(U_0, U_1) parameter values
rho.data <- read_csv(file.path(input.folder, "rho-cf-data.csv"))
print(rho.data)

# Load data from different sd(U_1) parameter values, relative to sd(U_0) = 1
sigma_1.data <- read_csv(file.path(input.folder, "sigma1-cf-data.csv"))
print(sigma_1.data)


################################################################################
## Normally dist errors: Plot dist of ADE and AIE estimates, around truth.

# ADE estimates, by type.
direct_dist.plot <- normal.data %>%
    ggplot() +
    # Dist of OLS estimates.
    geom_density(aes(x = (ols_direct_effect - truth_direct_effect),
        y = after_stat(density)), adjust = 2,
        colour = "black", fill = colour.list[1], alpha = 0.75) +
    annotate("text", colour = colour.list[1],
        x = 0.5, y = 7,
        fontface = "bold",
        label = ("Unadjusted"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = colour.list[1],
        x = 0.5, y = 6.875,
        xend = -0.3, yend = 5.25,
        linewidth = 0.75, curvature = -0.25,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    # Dist of CF estimates.
    geom_density(aes(x = (heckit_direct_effect - truth_direct_effect),
        y = after_stat(density)), adjust = 2,
        colour = "black", fill = colour.list[2], alpha = 0.75) +
    annotate("text", colour = colour.list[2],
        x = 0.5, y = 2.5,
        fontface = "bold",
        label = ("Parametric CF"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = colour.list[2],
        x = 0.5, y = 2.375,
        xend = 0.25, yend = 1.25,
        linewidth = 0.75, curvature = -0.25,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    # Truth value
    geom_vline(xintercept = 0,
        colour = "black", linetype = "dashed", linewidth = 1) +
    # Other presentation options
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = TeX("Estimate $-$ True Value"),
        breaks = seq(-1.0, 1.0, by = 0.25),
        limits = 0.81 * c(-1, 1)) +
    scale_y_continuous(expand = c(0, 0), name = "",
        breaks = seq(0, 10, by = 1), limits = c(0, 8.5)) +
    ggtitle("Density") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(output.folder, "normal-direct-dist.png"),
    plot = direct_dist.plot,
    units = "cm", width = fig.width, height = fig.height)

# AIE estimates, by type.
indirect_dist.plot <- normal.data %>%
    ggplot() +
    # Dist of OLS estimates.
    geom_density(aes(x = ols_indirect_effect - truth_indirect_effect,
        y = after_stat(density)), adjust = 2,
        colour = "black", fill = colour.list[1], alpha = 0.75) +
    # Dist of CF estimates.
    geom_density(aes(x = heckit_indirect_effect - truth_indirect_effect,
        y = after_stat(density)), adjust = 2, adjust = 2,
        colour = "black", fill = colour.list[2], alpha = 0.75) +
    # Truth value
    geom_vline(xintercept = 0,
        colour = "black", linetype = "dashed", linewidth = 1) +
    # Other presentation options
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = TeX("Estimate $-$ True Value"),
        breaks = seq(-1.0, 1.0, by = 0.25)         ,
        limits = 0.81 * c(-1, 1)) +
    scale_y_continuous(expand = c(0, 0), name = "",
        breaks = seq(0, 10, by = 1), limits = c(0, 8.5)) +
    ggtitle("Density") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(output.folder, "normal-indirect-dist.png"),
    plot = indirect_dist.plot,
    units = "cm", width = fig.width, height = fig.height)


################################################################################
## Uniform dist errors: Plot dist of ADE and AIE estimates, around truth.

# ADE estimates, by type.
direct_dist.plot <- uniform.data %>%
    ggplot() +
    # Dist of OLS estimates.
    geom_density(aes(x = (ols_direct_effect - truth_direct_effect),
        y = after_stat(density)), adjust = 2,
        colour = "black", fill = colour.list[1], alpha = 0.75) +
    annotate("text", colour = colour.list[1],
        x = 0.5, y = 7,
        fontface = "bold",
        label = ("Unadjusted"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = colour.list[1],
        x = 0.5, y = 6.875,
        xend = -0.3, yend = 5.25,
        linewidth = 0.75, curvature = -0.25,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    # Dist of Heckit estimates.
    geom_density(aes(x = (heckit_direct_effect - truth_direct_effect),
        y = after_stat(density)), adjust = 2,
        colour = "black", fill = colour.list[2], alpha = 0.75) +
    # Dist of CF estimates.
    geom_density(aes(x = (cf_direct_effect - truth_direct_effect),
        y = after_stat(density)), adjust = 2,
        colour = "black", fill = colour.list[3], alpha = 0.75) +
    annotate("text", colour = colour.list[3],
        x = 0.4, y = 4,
        fontface = "bold",
        label = ("Semi-parametric CF"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = colour.list[3],
        x = 0.5, y = 3.875,
        xend = 0.25, yend = 0.75,
        linewidth = 0.75, curvature = -0.25,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    # Truth value
    geom_vline(xintercept = 0,
        colour = "black", linetype = "dashed", linewidth = 1) +
    # Other presentation options
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = TeX("Estimate $-$ True Value"),
        breaks = seq(-1.0, 1.0, by = 0.25)         ,
        limits = 0.81 * c(-1.0, 1.0)) +
    scale_y_continuous(expand = c(0, 0), name = "",
        breaks = seq(0, 100, by = 1),
        limits = c(0, 8.5), oob = scales::rescale_none) +
    ggtitle("Density") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(output.folder, "uniform-direct-dist.png"),
    plot = direct_dist.plot,
    units = "cm", width = fig.width, height = fig.height)

# AIE estimates, by type.
indirect_dist.plot <- uniform.data %>%
    ggplot() +
    # Dist of OLS estimates.
    geom_density(aes(x = ols_indirect_effect - truth_indirect_effect,
        y = after_stat(density)), adjust = 2,
        colour = "black", fill = colour.list[1], alpha = 0.75) +
    # Dist of Heckit CF estimates.
    geom_density(aes(x = heckit_indirect_effect - truth_indirect_effect,
        y = after_stat(density)), adjust = 2,
        colour = "black", fill = colour.list[2], alpha = 0.75) +
    # Dist of semi-parametric CF estimates.
    geom_density(aes(x = cf_indirect_effect - truth_indirect_effect,
        y = after_stat(density)), adjust = 2,
        colour = "black", fill = colour.list[3], alpha = 0.75) +
    # Truth value
    geom_vline(xintercept = 0,
        colour = "black", linetype = "dashed", linewidth = 1) +
    # Other presentation options
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        name = TeX("Estimate $-$ True Value"),
        breaks = seq(-1.0, 1.0, by = 0.25)         ,
        limits = 0.81 * c(-1.0, 1.0)) +
    scale_y_continuous(expand = c(0, 0), name = "",
        breaks = seq(0, 100, by = 1),
        limits = c(0, 8.5), oob = scales::rescale_none) +
    ggtitle("Density") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))
# Save this plot
ggsave(file.path(output.folder, "uniform-indirect-dist.png"),
    plot = indirect_dist.plot,
    units = "cm", width = fig.width, height = fig.height)


################################################################################
## Plot the estimates, by OLS + CF, different rho = Corr(U_0, U_1) values.

# Plot the bias in direct effect est vs rho
rho_directeffect_bias.plot <- rho.data %>%
    ggplot(aes(x = rho)) +
    # OLS est + 95 % CI
    geom_point(aes(y = ols_direct_effect), colour = colour.list[1]) +
    geom_ribbon(aes(ymin = ols_direct_effect_low, ymax = ols_direct_effect_up),
        fill = colour.list[1], alpha = 0.2) +
    geom_line(aes(y = (ols_direct_effect_low)), alpha = 0.5,
        colour = colour.list[1], linetype = "dashed") +
    geom_line(aes(y = (ols_direct_effect_up)), alpha = 0.5,
        colour = colour.list[1], linetype = "dashed") +
    annotate("text", colour = colour.list[1],
        x = 0.25, y = 0.25,
        fontface = "bold",
        label = ("Unadjusted"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = colour.list[1],
        x = 0.55, y = 0.325,
        xend = 0.75, yend = 0.55,
        linewidth = 0.75,
        curvature = 0.25,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    # CF est + 95 % CI
    geom_point(aes(y = cf_direct_effect), colour = colour.list[3]) +
    geom_ribbon(aes(ymin = cf_direct_effect_low, ymax = cf_direct_effect_up),
        fill = colour.list[3], alpha = 0.2) +
    geom_line(aes(y = (cf_direct_effect_low)), alpha = 0.5,
        colour = colour.list[3], linetype = "dashed") +
    geom_line(aes(y = (cf_direct_effect_up)), alpha = 0.5,
        colour = colour.list[3], linetype = "dashed") +
    annotate("text", colour = colour.list[3],
        x = -0.5, y = 2.25,
        fontface = "bold",
        label = ("Semi-parametric CF"),
        size = 4.25, hjust = 0.5, vjust = 0) +
    annotate("curve", colour = colour.list[3],
        x = -0.65, y = 2.2,
        xend = -0.55, yend = 1.65,
        linewidth = 0.75,
        curvature = 0.125,
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
        x = 0.65, y = 1.75,
        xend = 0.7, yend = 1.45,
        linewidth = 0.75, curvature = -0.125,
        arrow = arrow(length = unit(0.25, 'cm'))) +
    # Presentation options
    theme_bw() +
    scale_x_continuous(
        name = TeX(r"(Corr$(U_{i,0}, U_{i,1})$)"),
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
    geom_line(aes(y = (ols_indirect_effect_low)), alpha = 0.5,
        colour = colour.list[1], linetype = "dashed") +
    geom_line(aes(y = (ols_indirect_effect_up)), alpha = 0.5,
        colour = colour.list[1], linetype = "dashed") +
    # CF est + 95 % CI
    geom_point(aes(y = cf_indirect_effect), colour = colour.list[3]) +
    geom_ribbon(aes(ymin = cf_indirect_effect_low, ymax = cf_indirect_effect_up),
        fill = colour.list[3], alpha = 0.2) +
    geom_line(aes(y = (cf_indirect_effect_low)), alpha = 0.5,
        colour = colour.list[3], linetype = "dashed") +
    geom_line(aes(y = (cf_indirect_effect_up)), alpha = 0.5,
        colour = colour.list[3], linetype = "dashed") +
    # Truth:
    geom_line(aes(y = (truth_indirect_effect)),
        colour = "black", linetype = "dashed", linewidth = 1) +
    # Presentation options
    theme_bw() +
    scale_x_continuous(
        name = TeX(r"(Corr$(U_{i,0}, U_{i,1})$)"),
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
    geom_ribbon(aes(ymin = ols_direct_effect_low, ymax = ols_direct_effect_up),
        fill = colour.list[1], alpha = 0.2) +
    geom_line(aes(y = (ols_direct_effect_low)), alpha = 0.5,
        colour = colour.list[1], linetype = "dashed") +
    geom_line(aes(y = (ols_direct_effect_up)), alpha = 0.5,
        colour = colour.list[1], linetype = "dashed") +
    # CF est + 95 % CI
    geom_point(aes(y = cf_direct_effect), colour = colour.list[3]) +
    geom_ribbon(aes(ymin = cf_direct_effect_low, ymax = cf_direct_effect_up),
        fill = colour.list[3], alpha = 0.2) +
    geom_line(aes(y = (cf_direct_effect_low)), alpha = 0.5,
        colour = colour.list[3], linetype = "dashed") +
    geom_line(aes(y = (cf_direct_effect_up)), alpha = 0.5,
        colour = colour.list[3], linetype = "dashed") +
    # Truth:
    geom_line(aes(y = truth_direct_effect),
        colour = "black", linetype = "dashed", linewidth = 1) +
    # Presentation options
    theme_bw() +
    scale_x_continuous(
        name = TeX(r"(Var$(U_{i,1})^{0.5}$)"),
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
ggsave(file.path(output.folder, "sigma1-directeffect-bias.png"),
    plot = sigma_1_directeffect_bias.plot, dpi = 300,
    units = "cm", width = fig.width, height = fig.height)

# Plot the bias in indirect effect est vs sigma_1
sigma_1_indirecteffect_bias.plot <- sigma_1.data %>%
    ggplot(aes(x = sigma_1)) +
    # Ols est + 95% CI
    geom_point(aes(y = ols_indirect_effect), colour = colour.list[1]) +
    geom_ribbon(aes(ymin = ols_indirect_effect_low, ymax = ols_indirect_effect_up),
        fill = colour.list[1], alpha = 0.2) +
    geom_line(aes(y = (ols_indirect_effect_low)), alpha = 0.5,
        colour = colour.list[1], linetype = "dashed") +
    geom_line(aes(y = (ols_indirect_effect_up)), alpha = 0.5,
        colour = colour.list[1], linetype = "dashed") +
    # CF est + 95 % CI
    geom_point(aes(y = cf_indirect_effect), colour = colour.list[3]) +
    geom_ribbon(aes(ymin = cf_indirect_effect_low, ymax = cf_indirect_effect_up),
        fill = colour.list[3], alpha = 0.2) +
    geom_line(aes(y = (cf_indirect_effect_low)), alpha = 0.5,
        colour = colour.list[3], linetype = "dashed") +
    geom_line(aes(y = (cf_indirect_effect_up)), alpha = 0.5,
        colour = colour.list[3], linetype = "dashed") +
    # Truth:
    geom_line(aes(y = truth_indirect_effect),
        colour = "black", linetype = "dashed", linewidth = 1) +
    # Presentation options
    theme_bw() +
    scale_x_continuous(
        name = TeX(r"(Var$(U_{i,1})^{0.5}$)"),
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
ggsave(file.path(output.folder, "sigma1-indirecteffect-bias.png"),
    plot = sigma_1_indirecteffect_bias.plot, dpi = 300,
    units = "cm", width = fig.width, height = fig.height)
