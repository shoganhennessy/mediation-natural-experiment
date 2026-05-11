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

# Factorise the relevant variables.
analysis.data$hh_size <- factor(analysis.data$hh_size)
analysis.data$initial_health_location <- factor(
    analysis.data$initial_health_location)


################################################################################
## Estimate a restricted MTE curve for healthcare utilisation.

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
        survey_weight) %>%
    drop_na()

mte.data$hh_size <- factor(mte.data$hh_size)
mte.data$initial_health_location <- factor(
    mte.data$initial_health_location)

# First stage for the mediator:
#   D_i = any_healthcare
#   excluded mediator IV = initial_health_location
mediator.firststage.reg <- glm(any_healthcare ~ lottery_iv +
        hh_size + initial_health_location,
    family = quasibinomial(link = "logit"),
    weights = survey_weight,
    data = mte.data)

mte.data$pi.est <- fitted(mediator.firststage.reg)

print(summary(mediator.firststage.reg))
summary(mte.data$pi.est)
hist(mte.data$pi.est)

## Estimate GAM outcome equation for happiness.
happy.gam <- gam(Y_happy ~ lottery_iv + hh_size +
        s(pi.est, bs = "cr"),
    weights = survey_weight,
    method = "REML",
    data = mte.data)

print(summary(happy.gam))
plot(happy.gam, pages = 1, shade = TRUE)


health.gam <- gam(Y_health ~ lottery_iv + hh_size + s(pi.est, bs = "cr"),
    weights = survey_weight,
    method = "REML",
    data = mte.data)
print(summary(health.gam))
plot(health.gam, pages = 1, shade = TRUE)


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

for (j in 1:NROW(p.grid)){

    data.plus <- mte.data
    data.minus <- mte.data

    data.plus$pi.est <- p.grid[j] + eps
    data.minus$pi.est <- p.grid[j] - eps

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

    mte.health.grid$mte[j] <- weighted.mean(
        (health.plus - health.minus) / (2 * eps),
        w = mte.data$survey_weight,
        na.rm = TRUE)

    mte.happy.grid$mte[j] <- weighted.mean(
        (happy.plus - happy.minus) / (2 * eps),
        w = mte.data$survey_weight,
        na.rm = TRUE)
}

mte.health.grid <- mte.health.grid %>%
    mutate(outcome_name = "Health overall good?")

mte.happy.grid <- mte.happy.grid %>%
    mutate(outcome_name = "Happy overall?")

# Plot the GAM-based MTE curves.
mte.gam.plot.data <- bind_rows(
    mte.health.grid,
    mte.happy.grid)

mte.gam.plot <- mte.gam.plot.data %>%
    ggplot(aes(x = pi.est, y = mte, colour = outcome_name)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_line(size = 1.1) +
    theme_bw() +
    scale_colour_manual("", values = colour.list[c(1, 3)]) +
    scale_x_continuous(name = TeX(r"(Estimated healthcare propensity, $\hat{\pi}$)")) +
    scale_y_continuous(name = "Estimated MTE of healthcare utilisation") +
    ggtitle(TeX(r"(GAM-based restricted MTE curve.)")) +
    theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 0, 0, 0), "mm"))

mte.gam.plot

ggsave(file.path(figures.folder, "oregon-gam-mte-shape.png"),
    plot = mte.gam.plot,
    units = "cm", width = fig.width, height = fig.height)
