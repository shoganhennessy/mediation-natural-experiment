#!/usr/bin/R
## Senan Hogan-Hennessy, 11 May 2026.
## Script to Visualise MTE in the Oregon Health Insurance experiment.
print(Sys.time())
set.seed(47)

## Packages:
renv::load()
# functions for data manipulation and visualisation
library(tidyverse)
# Semi-parametric splines
library(mgcv)
# Library for equations in plots
library(latex2exp)


# Define folder paths (1) input data (2) clean data.
data.folder <- file.path("..", "..", "data", "oregon-lottery-icspr")
figures.folder <- file.path("..", "..", "text", "sections", "figures")
tables.folder <- file.path("..", "..", "text", "sections", "tables")
presentation.folder <- file.path("..", "..", "presentation",
    "presentation-files", "figures")
# Size of figures.
fig.width <- 15
fig.height <- (2 / 3) * fig.width
fig.width <- fig.height
# List of 3 default colours.
colour.list <- c(
    "#1f77b4", # Blue
    "#2ca02c", # Green
    "#d62728") # Red


################################################################################
## Load the Oregon Health Insurance Experiment replication data.

# Load the pre-cleaned Oregon Health data.
analysis.data <- data.folder %>%
    file.path("cleaned-oregon-data.csv") %>%
    read_csv()


################################################################################
## Estimate the first-stage for the MTE model

# This estimates the MTE of any healthcare utilisation, D_i, on outcome Y_i.
# In the CM notation:
#   Z_i = lottery_iv
#   D_i = any_healthcare
#   Y_i = Y_health or Y_happy
#   mediator IV = initial_health_location

mte.data <- analysis.data %>%
    select(
        Y_health,
        Y_happy,
        any_healthcare,
        lottery_iv,
        hh_size,
        initial_health_location,
        initial_dia_diagnosis,
        initial_ast_diagnosis,
        initial_hbp_diagnosis,
        initial_emp_diagnosis,
        initial_chf_diagnosis,
        initial_dep_diagnosis,
        survey_weight) %>%
    drop_na()

mte.data$hh_size <- factor(mte.data$hh_size)
mte.data$initial_health_location <- factor(
    mte.data$initial_health_location)

# First stage for the mediator:
#   D_i = any_healthcare
#   excluded mediator IV = initial_health_location
mediator.firststage.reg <- glm(any_healthcare ~
    1 + lottery_iv + hh_size + initial_health_location +
    initial_dia_diagnosis +
    initial_ast_diagnosis +
    initial_hbp_diagnosis +
    initial_emp_diagnosis +
    initial_chf_diagnosis +
    initial_dep_diagnosis,
    family = quasibinomial(link = "probit"),
    weights = survey_weight,
    data = mte.data)
mte.data$pi.est <- fitted(mediator.firststage.reg)
print(summary(mediator.firststage.reg))
summary(mte.data$pi.est)

## Plot the estimated first-stage propensity score.
pi.plot <- mte.data %>%
    mutate(healthcare_label =
        ifelse(any_healthcare == 0,
            r"(Lost lottery, π(0; X))",
            r"(Won lottery, π(1; X))")) %>%
    ggplot() +
    geom_histogram(
        aes(x = pi.est, weight = survey_weight, fill = healthcare_label),
        bins = 10,
        colour = "black",
        alpha = 0.75) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0),
        limits = c(0.25, 1.05),
        #name = TeX(r"(Unobserved Cost to Healthcare Utilisation, \textit{$U_i$})"),
        name = TeX(r"(Estimated Healthcare Utilisation Propensity, \textit{$pi(z; \textbf{X}_i)$})"),
        breaks = seq(0, 1, by = 0.2)) +
    scale_y_continuous(expand = c(0.001, 0, 0.1, 0),
        name = "") +
    ggtitle(TeX(r"(Observations, \textit{$n$})")) +
    theme(legend.position = "none",
        plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm")) +
    facet_wrap(~ healthcare_label)
# Save the resulting figure.
ggsave(file.path(figures.folder, "oregon-pi-est.png"),
    plot = pi.plot,
    units = "cm", width = fig.width, height = fig.height)


################################################################################
## Local IV MTE curve for healthcare utilisation.

# Local Iv package
library(localIV)
localiv.bandwidth <- 0.2

# Estimate MTE by Local IV
health.mte <- mte(
    selection = any_healthcare ~
        1 + lottery_iv + hh_size + initial_health_location +
        initial_dia_diagnosis +
        initial_ast_diagnosis +
        initial_hbp_diagnosis +
        initial_emp_diagnosis +
        initial_chf_diagnosis +
        initial_dep_diagnosis,
    outcome = Y_health ~ 
        1 + lottery_iv + hh_size +
        initial_dia_diagnosis +
        initial_ast_diagnosis +
        initial_hbp_diagnosis +
        initial_emp_diagnosis +
        initial_chf_diagnosis +
        initial_dep_diagnosis,
    method = "localIV",
    data = mte.data,
    bw = localiv.bandwidth)
# fitted model
health_mte.data <- mte_at(u = seq(0.25, 0.95, 0.01), model = health.mte)

# Estimate MTE by Local IV
happy.mte <- mte(
    selection = any_healthcare ~
        1 + lottery_iv + hh_size + initial_health_location +
        initial_dia_diagnosis +
        initial_ast_diagnosis +
        initial_hbp_diagnosis +
        initial_emp_diagnosis +
        initial_chf_diagnosis +
        initial_dep_diagnosis,
    outcome = Y_happy ~ 
        1 + lottery_iv + hh_size +
        initial_dia_diagnosis +
        initial_ast_diagnosis +
        initial_hbp_diagnosis +
        initial_emp_diagnosis +
        initial_chf_diagnosis +
        initial_dep_diagnosis,
    method = "localIV",
    data = mte.data,
    bw = localiv.bandwidth)
# fitted model
happy_mte.data <- mte_at(u = seq(0.25, 0.95, 0.01), model = happy.mte)

# Collect the outcomes in the grids.
mte.health.grid <- health_mte.data %>%
    mutate(outcome_name = "Health overall good?")
mte.happy.grid <- happy_mte.data %>%
    mutate(outcome_name = "Happy overall?")
mte_plot.data <- bind_rows(mte.health.grid, mte.happy.grid)

# Calculate the extrapolation region
pi0 <- mean(mte.data$any_healthcare[mte.data$lottery_iv == 0])
pi1 <- mean(mte.data$any_healthcare[mte.data$lottery_iv == 1])

# Plot the MTE curves.
mte.plot <- mte_plot.data %>%
    mutate(mte.fill = ifelse(pi0 <= u & u <= pi1, value, NA)) %>%
    ggplot(aes(x = u)) +
    geom_line(aes(y = value, colour = outcome_name), linewidth = 1.5) +
    geom_ribbon(aes(ymin = 0, ymax = mte.fill, fill = outcome_name), alpha = 0.35) +
    #geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    scale_x_continuous(expand = c(0.01, 0.01),
        limits = c(0.3, 1),
        name = TeX(r"(Unobserved Cost to Healthcare Utilisation, \textit{$U_i$})"),
        breaks = seq(0, 1, by = 0.2)) +
    scale_y_continuous(expand = c(0.001, 0, 0.1, 0),
        limits = c(0, 0.75),
        breaks = seq(0, 1, by = 0.1),
        name = "") +
    ggtitle(TeX(r"(Local IV Point Estimate)")) +
    theme(legend.position = c(0.75, 0.96),
        legend.title = element_blank(),
        plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))

# Save the resulting figure.
ggsave(file.path(figures.folder, "oregon-mte-est.png"),
    plot = mte.plot,
    units = "cm", width = fig.width, height = fig.height)
