#!/usr/bin/R
## Senan Hogan-Hennessy, 16 Jan 2025
## SimulaPlotte the system for indirect + direct effects, with Roy selection.
# Show the date:
print(format(Sys.time(), "%H:%M %Z %A, %d %B %Y"))

## Load libraries
# Functions for data manipulation and visualisation
library(tidyverse)
# Library for better colour choice.
library(ggthemes)
# Library for equations in plots
library(latex2exp)

# Define colours to use.
colour.list <- c("orange", "blue")


################################################################################
## Load the simulated data.  (saved as separate file in advance.)

# Load the data file.
rho.data <- read_csv(file.path(output.folder, "rho-sim-data.csv"))
print(names(rho.data))


################################################################################
## Plot the Direct effect estimates, by OLS + CF

# Plot the bias in direct effect est vs rho, with different sigma values
directeffect_bias.plot <- rho.data %>%
    ggplot(aes(x = rho)) +
    geom_line(aes(y = truth_direct_effect),
        colour = "black", linetype = "dashed") +
    geom_point(aes(y = ols_direct_effect), colour = colour.list[1]) +
    geom_point(aes(y = cf_direct_effect), colour = colour.list[2]) +
    theme_bw() +
    scale_x_continuous(
        name = TeX("$\\rho$"),
        expand = c(0, 0),
        breaks = seq(-1, 1, by = 0.2),
        limits = c(-1.05, 1.05)) +
    scale_y_continuous(name = "") +#,
        #breaks = seq(-1, 1, by = 0.1),
        #limits = c(0, 2),
        #expand = c(0, 0)) +
    ggtitle("Estimate") +
    theme(plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0.25, 0.25), "mm"),
        legend.position = "bottom")
# Save this plot
ggsave(file.path(output.folder, "directeffect-bias.png"),
    plot = directeffect_bias.plot, dpi = 300,
    units = "cm", width = fig.width, height = fig.height)


################################################################################
## Plot the Indirect effect estimates, by OLS + CF

# Plot the bias in indirect effect est vs rho, with different sigma values
indirecteffect_bias.plot <- rho.data %>%
    ggplot(aes(x = rho)) +
    geom_line(aes(y = truth_indirect_effect),
        colour = "black", linetype = "dashed") +
    geom_point(aes(y = ols_indirect_effect), colour = colour.list[1]) +
    geom_point(aes(y = cf_indirect_effect), colour = colour.list[2]) +
    theme_bw() +
    scale_x_continuous(
        name = TeX("$\\rho$"),
        expand = c(0, 0),
        breaks = seq(-1, 1, by = 0.2),
        limits = c(-1.05, 1.05)) +
    scale_y_continuous(name = "") +#,
        #breaks = seq(-1, 1, by = 0.1),
        #limits = c(0, 2),
        #expand = c(0, 0)) +
    ggtitle("Estimate") +
    theme(plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0.25, 0.25), "mm"),
        legend.position = "bottom")
# Save this plot
ggsave(file.path(output.folder, "indirecteffect-bias.png"),
    plot = indirecteffect_bias.plot, dpi = 300,
    units = "cm", width = fig.width, height = fig.height)
