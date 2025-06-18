#!/usr/bin/R
## Senan Hogan-Hennessy, 14 June 2025.
## Script to to extract relevant data from Oregon Health Insurance rep package.
print(Sys.time())
set.seed(47)

## Packages:
# functions for data manipulation and visualisation
library(tidyverse)

# Define folder paths (1) input data (2) clean data.
data.folder <- file.path("..", "..", "data", "oregon-lottery-icspr")
figures.folder <- file.path("..", "..", "text", "sections", "figures")
presentation.folder <- file.path("..", "..", "presentation",
    "presentation-files", "figures")
# Size of figures.
fig.width <- 15
fig.height <- (2 / 3) * fig.width
presentation.width <- 15
presentation.height <- (2 / 3) * presentation.width
# List of 3 default colours.
colour.list <- c("orange", "blue")
    #"#1f77b4", # Blue
    #"#2ca02c", # Green
    #"#d62728") # Red

################################################################################
## Load the Oregon Health Insurance Experiment replication data.

# Load the pre-cleaned Oregon Health data.
oregon.data <- data.folder %>%
    file.path("cleaned-oregon-data.csv") %>%
    read_csv()

################################################################################
## Validate the lottery encouragement instrument for having health insurance.

# Show lottery IV -> insurance is strong.
iv_firststage.reg <- lm(any_insurance ~ 1 + lottery_iv,
    data = oregon.data)
print(summary(iv_firststage.reg))

# Show lottery IV -> insurance is strong, conditional on household size.
iv_firststage.reg <- lm(any_insurance ~ 1 + lottery_iv * factor(hh_size),
    data = oregon.data)
print(summary(iv_firststage.reg))

# Calculate Pr(Z_iv = 1 | household size), the known instrument prop score.
iv_prop.reg <- lm(lottery_iv ~ 0 + factor(hh_size),
    data = oregon.data)
print(summary(iv_prop.reg))


################################################################################
## Create a figure showing the happinness + subjective health effects.

## Get relevant data.
analysis.data <- oregon.data %>%
    select(health_level_survey, happiness_level_survey, lottery_iv,
        any_insurance, health_needs_met,
        any_doc_visits, any_hospital_visits, hh_size) %>%
    drop_na()

#! Test: bootstrap
#analysis.data <- analysis.data %>%
#    sample_frac(replace = TRUE)

# Health survey outcome -> Overall health is good or better (2 is fair, 1 is poor).
Y_health <- as.integer(analysis.data$health_level_survey >= 3)
# Happiness outcome -> very or pretty overall happiness (1 == not happy)
Y_happy <- as.integer(analysis.data$happiness_level_survey <= 2)
# Lottery instrument for health insurance.
Z_iv <- analysis.data$lottery_iv
# Treatment -> have any health insurance insurance.
Z <- analysis.data$any_insurance
# Mediator -> health_needs_met
D <- analysis.data$health_needs_met
#! Test: DOctor visits
D <- as.integer((analysis.data$any_doc_visits +
    analysis.data$any_hospital_visits) > 0)

## Use the Abadie (2003) kappa weights to get E[ Y(d') | lottery complier ]
est_probZ_iv <- glm(Z_iv ~ 1 + factor(analysis.data$hh_size))
hat_probZ_iv <- est_probZ_iv$fitted
kappa_1 <- Z * ((Z_iv - hat_probZ_iv) / ((1 - hat_probZ_iv) * hat_probZ_iv))
kappa_0 <- (1 - Z) * (((1 - Z_iv) - (1 - hat_probZ_iv)) / ((1 - hat_probZ_iv) * hat_probZ_iv))
kappa <- kappa_0 * (1 - hat_probZ_iv) + kappa_1 * hat_probZ_iv

## Estimate mean outcomes among (un)insured lottery compliers.
# Health insurance effect of the lottery (among entire population).
Z_0_complier      <- mean(Z[Z_iv == 0])
Z_1_complier      <- mean(Z[Z_iv == 1])
Z_effect_complier <- Z_1_complier - Z_0_complier
print(Z_effect_complier)
# Health needs met?
D_0_complier      <- mean(kappa_0 * D) / mean(kappa_0)
D_1_complier      <- mean(kappa_1 * D) / mean(kappa_1)
D_effect_complier <- D_1_complier - D_0_complier
print(D_effect_complier)
# Health overall good?
Y_health_0_complier      <- mean(kappa_0 * Y_health) / mean(kappa_0)
Y_health_1_complier      <- mean(kappa_1 * Y_health) / mean(kappa_1)
Y_health_effect_complier <- Y_health_1_complier - Y_health_0_complier
print(Y_health_effect_complier)
# Happy overall?
Y_happy_0_complier      <- mean(kappa_0 * Y_happy) / mean(kappa_0)
Y_happy_1_complier      <- mean(kappa_1 * Y_happy) / mean(kappa_1)
Y_happy_effect_complier <- Y_happy_1_complier - Y_happy_0_complier
print(Y_happy_effect_complier)


################################################################################
## Bar chart of the health insurance effects

# Name of the outcome variables (in order)
outcome_name.list <- c("Health insured?",
    "Survey: \nAny healthcare visits?",
    "Survey: \nHealth overall good?",
    "Survey: \nHappy overall?")

# Get a dataframe of the relevant effects.
complier.data <- data.frame(
    Z_iv = c("0", "1"),
    outcome_value = c(
        Z_0_complier, Z_1_complier,
        D_0_complier, D_1_complier,
        Y_health_0_complier, Y_health_1_complier,
        Y_happy_0_complier, Y_happy_1_complier),
    outcome_name = c(
        rep(outcome_name.list[1], 2),
        rep(outcome_name.list[2], 2),
        rep(outcome_name.list[3], 2),
        rep(outcome_name.list[4], 2)))

# Show a barchart of just the lottery effects.
lottery.plot <- complier.data %>%
    mutate(outcome_value =
        ifelse(outcome_name != outcome_name.list[1], 0, outcome_value)) %>%
    ggplot() +
    geom_bar(aes(fill = Z_iv, x = outcome_name, y = outcome_value),
        colour = 1, position = "dodge", stat = "identity") +
    theme_bw() +
    scale_x_discrete(name = "") +
    scale_fill_manual("", values = colour.list) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        limits = c(0, 1.025),
        breaks = seq(0, 1, by = 0.1)) +
    ggtitle("Mean Outcome") +
    guides(fill = guide_legend(ncol = 2)) + 
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        legend.position = c(0.9125, 1.05),
        plot.margin = unit(c(0.5, 3, 0, 0), "mm")) +
    # Add a caliper noting the lottery effects
    ggbrace::stat_brace(
        data = data.frame(x = c(0.6, 1.4), y = c(0.75, 0.85)), aes(x, y),
        size = 1, colour = "black") +
    annotate("text", x = 1, y = 0.94,
        label = ("Lottery effect"),
        size = 4, hjust = 0.5, vjust = 0) +
    # Label the effect sizes.
    annotate("text", x = 0.8, y = Z_0_complier + Z_effect_complier / 2,
        label = paste0("+ ", round(100 * Z_effect_complier), "%"),
        size = 4, hjust = 0.5, vjust = 0.5,
        fontface = "bold", colour = colour.list[2])
# Save this plot
ggsave(file.path(presentation.folder, "lottery-effects.png"),
    plot = lottery.plot,
    units = "cm", width = presentation.width, height = presentation.height)

# Full barchart
complier.plot <- complier.data %>%
    ggplot() +
    geom_bar(aes(fill = Z_iv, x = outcome_name, y = outcome_value),
        colour = 1, position = "dodge", stat = "identity") +
    theme_bw() +
    scale_x_discrete(name = "") +
    scale_fill_manual("", values = colour.list) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        limits = c(0, 1.025),
        breaks = seq(0, 1, by = 0.1)) +
    ggtitle("Mean Outcome") +
    guides(fill = guide_legend(ncol = 2)) + 
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        legend.position = c(0.9125, 1.05),
        plot.margin = unit(c(0.5, 3, 0, 0), "mm")) +
    # Add a caliper noting the lottery effects
    ggbrace::stat_brace(
        data = data.frame(x = c(0.6, 1.4), y = c(0.75, 0.85)), aes(x, y),
        size = 1, colour = "black") +
    annotate("text", x = 1, y = 0.94,
        label = ("Lottery effect"),
        size = 4, hjust = 0.5, vjust = 0) +
    # Add a caliper noting the health insurance effects (lottery compliers).
    ggbrace::stat_brace(
        data = data.frame(x = c(1.6, 4.4), y = c(0.725, 0.8125)), aes(x, y),
        size = 1, colour = "black") +
    annotate("text", x = 4.25, y = 0.94,
        label = ("Health insurance effect (lottery compliers)"),
        size = 4, hjust = 1, vjust = 0) +
    # Label the effect sizes.
    annotate("text", x = 0.8, y = Z_0_complier + Z_effect_complier / 2,
        label = paste0("+ ", round(100 * Z_effect_complier), "%"),
        size = 4, hjust = 0.5, vjust = 0.5,
        fontface = "bold", colour = colour.list[2]) +
    annotate("text", x = 1.8, y = D_0_complier + D_effect_complier / 2,
        label = paste0("+ ", round(100 * D_effect_complier), "%"),
        size = 4, hjust = 0.5, vjust = 0.5,
        fontface = "bold", colour = colour.list[2]) +
    annotate("text", x = 2.8, y = Y_happy_0_complier + Y_happy_effect_complier / 2,
        label = paste0("+ ", round(100 * Y_happy_effect_complier), "%"),
        size = 4, hjust = 0.5, vjust = 0.5,
        fontface = "bold", colour = colour.list[2]) +
    annotate("text", x = 3.8, y = Y_health_0_complier + Y_health_effect_complier / 2,
        label = paste0("+ ", round(100 * Y_health_effect_complier), "%"),
        size = 4, hjust = 0.5, vjust = 0.5,
        fontface = "bold", colour = colour.list[2])

# Save this plot
ggsave(file.path(presentation.folder, "insurance-effects.png"),
    plot = complier.plot,
    units = "cm", width = presentation.width, height = presentation.height)
ggsave(file.path(figures.folder, "insurance-effects.png"),
    plot = complier.plot,
    units = "cm", width = fig.width, height = fig.height)

# Code from LARF package for estimating OLS with possibly negative kappa weights.
X <- 1
solve ( t(cbind(Z,X) * kappa) %*% cbind(Z,X)) %*% t(cbind(Z,X) * kappa)  %*% Y_health



lm(any_insurance ~ 1 + lottery_iv * factor(hh_size), data = oregon.data) %>%
    summary() 

# Effect of insurance on healthcare (mediator)
lm(health_needs_met ~ 1 + lottery_iv * factor(hh_size), data = oregon.data) %>%
    summary() 
ivreg::ivreg(health_needs_met ~ 1 + any_insurance | 1 + lottery_iv * factor(hh_size),
    data = oregon.data) %>%
    summary()
# Effect on resulting happiness + health (ATE)
ivreg::ivreg(I(health_level_survey >= 3) ~ 1 + any_insurance | 1 + lottery_iv * factor(hh_size),
    data = oregon.data) %>%
    summary()
ivreg::ivreg(I(happiness_level_survey <= 2) ~ 1 + any_insurance | 1 + lottery_iv * factor(hh_size),
    data = oregon.data) %>%
    summary()
# TODO: Use the Abadaie (2003) estimator, to simplify the later mediation.


# Informal mechanism, direct effect (controlling for healthcare satisfaction)
ivreg::ivreg(I(health_level_survey >= 3) ~ 1 + any_insurance + health_needs_met
    | 1 + lottery_iv * factor(hh_size) + health_needs_met,
    data = oregon.data) %>%
    summary()
ivreg::ivreg(I(happiness_level_survey <= 2) ~ 1 + any_insurance + health_needs_met
    | 1 + lottery_iv * factor(hh_size) + health_needs_met,
    data = oregon.data) %>%
    summary()

# Using the cost of provider as instruments for the healthcare usage.
# -> Direct effect of insurance on health/happiness now zero.
ivreg::ivreg(I(health_level_survey >= 3)
    ~ 1 + any_insurance + health_needs_met
    | 1 + lottery_iv * factor(hh_size) + factor(usual_health_location),
    data = oregon.data) %>%
    summary()
ivreg::ivreg(I(happiness_level_survey <= 2)
    ~ 1 + any_insurance + health_needs_met
    | 1 + lottery_iv * factor(hh_size) + factor(usual_health_location),
    data = oregon.data) %>%
    summary()
