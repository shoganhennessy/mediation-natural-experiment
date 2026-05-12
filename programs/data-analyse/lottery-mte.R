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
    1 + lottery_iv * hh_size * initial_health_location,
    family = quasibinomial(link = "logit"),
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
    scale_x_continuous(expand = c(0.01, 0.01),
        limits = c(0.4, 0.85),
        #name = TeX(r"(Unobserved Cost to Healthcare Utilisation, \textit{$U_i$})"),
        name = TeX(r"(Estimated Healthcare Utilisation Propensity, \textit{$pi(z; \textbf{X}_i)$})"),
        breaks = seq(0, 1, by = 0.1)) +
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

library(localIV)

mod <- mte(
    selection = any_healthcare ~ 1 + lottery_iv + hh_size +
        initial_health_location,
    outcome = Y_happy ~ 1 + lottery_iv + hh_size,
    method = "localIV",
    data = mte.data,
    bw = 0.1)
# fitted propensity score model
summary(mod$ps_model)
mte_vals <- mte_at(u = seq(0.3, 0.95, 0.01), model = mod)

# SHow the Local IV MTE(u) curve
localiv.plot <- ggplot(mte_vals, aes(x = u, y = value)) +
  geom_line(size = 1) +
  xlab("Latent Resistance U") +
  ylab("Estimates of MTE at Average values of X") +
  theme_minimal(base_size = 14)


################################################################################
## Estimate a restricted MTE curve for healthcare utilisation.

# Estimate GAM outcome equation for health well-being
health.gam <- gam(Y_health ~ 1 + lottery_iv * hh_size +
    s(pi.est, bs = "cr"),
    weights = survey_weight,
    method = "REML",
    data = mte.data)
print(summary(health.gam))
plot(health.gam, pages = 1, shade = TRUE)

# Estimate GAM outcome equation for happiness.
happy.gam <- gam(Y_happy ~ 1 + lottery_iv * hh_size +
    s(pi.est, bs = "cr"),
    weights = survey_weight,
    method = "REML",
    data = mte.data)
print(summary(happy.gam))
plot(happy.gam, pages = 1, shade = TRUE)

# Recover MTE curve from numerical derivative of the GAM.
alpha <- 0.05
p.lower <- quantile(mte.data$pi.est, alpha / 2, na.rm = TRUE)
p.upper <- quantile(mte.data$pi.est, 1 - alpha / 2, na.rm = TRUE)
p.grid <- seq(p.lower, p.upper, length.out = 100)

# Small perturbation for central finite differences.
eps <- 0.001 * (p.upper - p.lower)
mte.health.grid <- data.frame(
    pi.est = p.grid,
    mte = NA_real_)
mte.happy.grid <- data.frame(
    pi.est = p.grid,
    mte = NA_real_)

# Operate the grid
for (j in 1:NROW(p.grid)){
    data.plus <- mte.data
    data.minus <- mte.data
    data.plus$pi.est <- p.grid[j] + eps
    data.minus$pi.est <- p.grid[j] - eps
    # Project the pi predictions around the epsilon neighbourhood
    health.plus <- predict(health.gam,
        newdata = data.plus,
        type = "response")
    health.minus <- predict(health.gam,
        newdata = data.minus,
        type = "response")
    happy.plus <- predict(happy.gam,
        newdata = data.plus,
        type = "response")
    happy.minus <- predict(happy.gam,
        newdata = data.minus,
        type = "response")
    # Calculate the implied slope.
    mte.health.grid$mte[j] <- weighted.mean(
        (health.plus - health.minus) / (2 * eps),
        w = mte.data$survey_weight,
        na.rm = TRUE)
    mte.happy.grid$mte[j] <- weighted.mean(
        (happy.plus - happy.minus) / (2 * eps),
        w = mte.data$survey_weight,
        na.rm = TRUE)
}

# Collect the outcomes in the grids.
mte.health.grid <- mte.health.grid %>%
    mutate(outcome_name = "Health overall good?")
mte.happy.grid <- mte.happy.grid %>%
    mutate(outcome_name = "Happy overall?")
mte.gam.plot.data <- bind_rows(
    mte.health.grid,
    mte.happy.grid)

# Calculate the extrapolation region
pi0 <- mean(mte.data$any_healthcare[mte.data$lottery_iv == 0])
pi1 <- mean(mte.data$any_healthcare[mte.data$lottery_iv == 1])

# Plot the GAM-based MTE curves.
mte.plot <- mte.gam.plot.data %>%
    mutate(mte.fill = ifelse(pi0 <= pi.est & pi.est <= pi1, mte, NA)) %>%
    ggplot(aes(x = pi.est)) +
    geom_line(aes(y = mte, colour = outcome_name), linewidth = 1.5) +
    geom_ribbon(aes(ymin = 0, ymax = mte.fill, fill = outcome_name), alpha = 0.35) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    scale_x_continuous(expand = c(0.01, 0.01),
        limits = c(0.4, 0.825),
        name = TeX(r"(Unobserved Cost to Healthcare Utilisation, \textit{$U_i$})"),
        breaks = seq(0, 1, by = 0.1)) +
    scale_y_continuous(expand = c(0.001, 0, 0.1, 0),
        breaks = seq(0, 1, by = 0.05),
        name = "") +
    ggtitle(TeX(r"(MTE estimate, \textit{E$[ Y_i(1) - Y_i(0) \; | \; U_i \, ]$})")) +
    theme(legend.position = c(0.25, 0.75),
        legend.title = element_blank(),
        plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0, 0), "mm"))

# Save the resulting figure.
ggsave(file.path(figures.folder, "oregon-mte-est.png"),
    plot = mte.plot,
    units = "cm", width = fig.width, height = fig.height)
