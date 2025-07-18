#!/usr/bin/R
## Senan Hogan-Hennessy, 14 June 2025.
## Script to to extract relevant data from Oregon Health Insurance rep package.
print(Sys.time())
set.seed(47)

## Packages:
# functions for data manipulation and visualisation
library(tidyverse)
# Functions for bootstrapping.
library(boot)
# Library for better colour choice.
library(ggthemes)
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
presentation.width <- 15
presentation.height <- (2 / 3) * presentation.width
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
analysis.data$usual_health_location <- factor(
    analysis.data$usual_health_location)


################################################################################
## Validate the lottery encouragement instrument for having health insurance.

# Show lottery IV -> insurance is strong.
iv_firststage.reg <- lm(any_insurance ~ 1 + lottery_iv,
    weights = survey_weight,
    data = analysis.data)
print(summary(iv_firststage.reg))

# Show lottery IV -> insurance is strong, conditional on household size.
iv_firststage.reg <- lm(any_insurance ~ 1 + lottery_iv * factor(hh_size),
    weights = survey_weight,
    data = analysis.data)
print(summary(iv_firststage.reg))

# Calculate Pr(Z_iv = 1 | household size), the known instrument prop score.
iv_prop.reg <- lm(lottery_iv ~ 0 + factor(hh_size),
    weights = survey_weight,
    data = analysis.data)
print(summary(iv_prop.reg))


################################################################################
## Create a figure showing the happiness + subjective health effects.

# Use the Abadie (2003) kappa weights to get E[ Y(z') | lottery complier ]
abadie.late <- function(data, outcome, Z, Z_iv, control_iv,
        indices = NULL){
    # Bootstrap sample, if indices provided.
    if (!is.null(indices)){
        data <- data[indices,]
    }
    # Get the relevant variables.
    outcome <- data[[outcome]]
    Z          <- data[[Z]]
    Z_iv       <- data[[Z_iv]]
    control_iv <- data[[control_iv]]
    # Estimate the kappa weighting
    est_probZ_iv <- glm(Z_iv ~ 1 + control_iv)
    hat_probZ_iv <- est_probZ_iv$fitted
    kappa_0 <- (1 - Z) * (((1 - Z_iv) - (1 - hat_probZ_iv)) / 
        ((1 - hat_probZ_iv) * hat_probZ_iv))
    kappa_1 <- Z * ((Z_iv - hat_probZ_iv) / (
        (1 - hat_probZ_iv) * hat_probZ_iv))
    kappa <- kappa_0 * (1 - hat_probZ_iv) + kappa_1 * hat_probZ_iv
    # Calculate the complier levels for z' = 0, 1, and resulting effect.
    outcome_0_complier      <- mean(kappa_0 * outcome) / mean(kappa_0)
    outcome_1_complier      <- mean(kappa_1 * outcome) / mean(kappa_1)
    outcome_effect_complier <- outcome_1_complier - outcome_0_complier
    return(c(outcome_0_complier, outcome_1_complier, outcome_effect_complier))
}


## Estimate mean outcomes among (un)insured lottery compliers, with boot SEs
boot.samples <- 10^3
# Health insurance effect of the lottery (among entire population).
Z.late <- boot(statistic = abadie.late, R = boot.samples,
    data = analysis.data,
    outcome = "any_insurance", Z = "lottery_iv",
    Z_iv = "lottery_iv", control_iv = "hh_size")
Z_0_complier      <- mean(Z.late$t[, 1])
Z_1_complier      <- mean(Z.late$t[, 2])
Z_effect_complier <- mean(Z.late$t[, 3])
Z_0_complier.se      <- sd(Z.late$t[, 1])
Z_1_complier.se      <- sd(Z.late$t[, 2])
Z_effect_complier.se <- sd(Z.late$t[, 3])
print(c(Z_effect_complier, Z_effect_complier.se))
# Mediator: healthcare utilisation (among lottery compliers).
D.late <- boot(statistic = abadie.late, R = boot.samples,
    data = analysis.data,
    outcome = "any_healthcare", Z = "any_insurance",
    Z_iv = "lottery_iv", control_iv = "hh_size")
D_0_complier      <- mean(D.late$t[, 1])
D_1_complier      <- mean(D.late$t[, 2])
D_effect_complier <- mean(D.late$t[, 3])
D_0_complier.se      <- sd(D.late$t[, 1])
D_1_complier.se      <- sd(D.late$t[, 2])
D_effect_complier.se <- sd(D.late$t[, 3])
print(c(D_effect_complier, D_effect_complier.se))
# Outcome: Health overall good? (among lottery compliers).
Y_health.late <- boot(statistic = abadie.late, R = boot.samples,
    data = analysis.data,
    outcome = "Y_health", Z = "any_insurance",
    Z_iv = "lottery_iv", control_iv = "hh_size")
Y_health_0_complier      <- mean(Y_health.late$t[, 1])
Y_health_1_complier      <- mean(Y_health.late$t[, 2])
Y_health_effect_complier <- mean(Y_health.late$t[, 3])
Y_health_0_complier.se      <- sd(Y_health.late$t[, 1])
Y_health_1_complier.se      <- sd(Y_health.late$t[, 2])
Y_health_effect_complier.se <- sd(Y_health.late$t[, 3])
print(c(Y_health_effect_complier, Y_health_effect_complier.se))
# Outcome: Happy overall?  (among lottery compliers).
Y_happy.late <- boot(statistic = abadie.late, R = boot.samples,
    data = analysis.data,
    outcome = "Y_happy", Z = "any_insurance",
    Z_iv = "lottery_iv", control_iv = "hh_size")
Y_happy_0_complier      <- mean(Y_happy.late$t[, 1])
Y_happy_1_complier      <- mean(Y_happy.late$t[, 2])
Y_happy_effect_complier <- mean(Y_happy.late$t[, 3])
Y_happy_0_complier.se      <- sd(Y_happy.late$t[, 1])
Y_happy_1_complier.se      <- sd(Y_happy.late$t[, 2])
Y_happy_effect_complier.se <- sd(Y_happy.late$t[, 3])
print(c(Y_happy_effect_complier, Y_happy_effect_complier.se))


################################################################################
## Bar chart of the health insurance effects

# Name of the outcome variables (in order)
outcome_name.list <- c(
    "0               1\nHealth insured?",
    "0               1\nAny use of healthcare?",
    "0               1\nSurvey: \nHealth overall good?",
    "0               1\nSurvey: \nHappy overall?")

# Get a dataframe of the relevant effects.
complier.data <- data.frame(
    Z_iv = c("0", "1"),
    outcome_value = c(
        Z_0_complier, Z_1_complier,
        D_0_complier, D_1_complier,
        Y_health_0_complier, Y_health_1_complier,
        Y_happy_0_complier, Y_happy_1_complier),
    #ci_lower = c(
    #    Z_0_lower, Z_1_lower,
    #    D_0_lower, D_1_lower,
    #    Y_health_0_lower, Y_health_1_lower, 
    #    Y_happy_0_lower, Y_happy_1_lower),
    #ci_upper = c(
    #    Z_0_upper, Z_1_upper,
    #    D_0_upper, D_1_upper,
    #    Y_health_0_upper, Y_health_1_upper, 
    #    Y_happy_0_upper, Y_happy_1_upper),
    outcome_name = c(
        rep(outcome_name.list[1], 2),
        rep(outcome_name.list[2], 2),
        rep(outcome_name.list[3], 2),
        rep(outcome_name.list[4], 2)))

# Full barchart
complier.plot <- complier.data %>%
    ggplot() +
    geom_bar(aes(group = Z_iv,
        fill = outcome_name, x = outcome_name, y = outcome_value),
        colour = 1, position = "dodge", stat = "identity") +
    theme_bw() +
    scale_x_discrete(name = "", limits = outcome_name.list) +
    scale_fill_manual("", values = colour.list[c(2, 1, 3, 3)]) +
    scale_y_continuous(expand = c(0, 0),
        name = "",
        limits = c(0, 1.025),
        breaks = seq(0, 1, by = 0.1)) +
    ggtitle(TeX(r"(Mean Outcome, for each $z' =0,1$.)")) +
    theme(legend.position = "none",
        plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 0, -2.5, 0), "mm")) +
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
        label = paste0("+ ", round(Z_effect_complier, 2),
            "\n(", round(Z_effect_complier.se, 2), ")"),
        size = 4, hjust = 0.5, vjust = 0.5,
        fontface = "bold", colour = colour.list[1]) +
    annotate("text", x = 1.8, y = D_0_complier + D_effect_complier / 2,
        label = paste0("+ ", round(D_effect_complier, 2),
            "\n(", round(D_effect_complier.se, 2), ")"),
        size = 4, hjust = 0.5, vjust = 0.5,
        fontface = "bold", colour = colour.list[2]) +
    annotate("text", x = 2.8, y = Y_health_0_complier + Y_health_effect_complier / 2,
        label = paste0("+ ", round(Y_health_effect_complier, 2),
            "\n(", round(Y_health_effect_complier.se, 2), ")"),
        size = 4, hjust = 0.5, vjust = 0.5,
        fontface = "bold", colour = colour.list[3])  +
    annotate("text", x = 3.8, y = Y_happy_0_complier + Y_happy_effect_complier / 2,
        label = paste0("+ ", round(Y_happy_effect_complier, 2),
            "\n(", round(Y_happy_effect_complier.se, 2), ")"),
        size = 4, hjust = 0.5, vjust = 0.5,
        fontface = "bold", colour = colour.list[3])

# Save this plot
ggsave(file.path(presentation.folder, "insurance-effects.png"),
    plot = complier.plot,
    units = "cm", width = presentation.width, height = presentation.height)
ggsave(file.path(figures.folder, "insurance-effects.png"),
    plot = complier.plot,
    units = "cm", width = fig.width, height = fig.height)
