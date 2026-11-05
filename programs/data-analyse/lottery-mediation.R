#!/usr/bin/R
## Senan Hogan-Hennessy, 5 November 2026.
## Estimate CM effects in the Oregon Health Insurance Experiment.
## Updated to identify AIE(0) from the Z = 0 conditional-outcome curve,
## then estimate ADE(1) = ATE - AIE(0).
print(Sys.time())
set.seed(47)

## Load packages, in the RENV environment
renv::load()
# Functions for data manipulation and visualisation.
library(tidyverse)
# Functions for bootstrapping and applied regressions.
library(boot)
library(fixest)
# Library for better colour choice.
library(ggthemes)
# Library for equations in plots
library(latex2exp)
# Package for LaTeX tables
library(xtable)

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
presentation.height <- (7 / 12) * presentation.width
# Number of digits to round to.
digits.no <- 2
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

# List all the observed control variables.
diagnosis.list <- c(
    "initial_dia_diagnosis",
    "initial_ast_diagnosis",
    "initial_hbp_diagnosis",
    "initial_emp_diagnosis",
    "initial_chf_diagnosis",
    "initial_dep_diagnosis")
controls.list <- paste(c("hh_size", diagnosis.list), collapse = " + ")
nocontrols.list <- "hh_size"


################################################################################
## Define the estimation functions.

# Create counterfactual treatment copies for prediction.
treatment.copies <- function(example.data){
    example_Z0.data <- example.data
    example_Z1.data <- example.data
    example_Z0.data$Z <- 0
    example_Z1.data$Z <- 1
    return(list(Z0 = example_Z0.data, Z1 = example_Z1.data))
}

# Summarise whether endpoint propensity scores lie in the observed Z = 0 range.
# With the categorical IV, the restricted polynomial and rank condition provide
# identification. These quantities remain useful warnings about extrapolation.
support.summary <- function(pi_realised_Z0, pi_0_Z0, pi_1_Z0){
    support.lower <- min(pi_realised_Z0, na.rm = TRUE)
    support.upper <- max(pi_realised_Z0, na.rm = TRUE)
    endpoint.lower <- pmin(pi_0_Z0, pi_1_Z0)
    endpoint.upper <- pmax(pi_0_Z0, pi_1_Z0)
    return(list(
        support.lower = support.lower,
        support.upper = support.upper,
        share.pi0.in.support = mean(pi_0_Z0 >= support.lower &
            pi_0_Z0 <= support.upper, na.rm = TRUE),
        share.pi1.in.support = mean(pi_1_Z0 >= support.lower &
            pi_1_Z0 <= support.upper, na.rm = TRUE),
        share.interval.in.support = mean(endpoint.lower >= support.lower &
            endpoint.upper <= support.upper, na.rm = TRUE)))
}

# Estimate mediation without correcting for selection into the mediator.
mediate.unadjusted <- function(Y, Z, D, X_iv, data,
        X_minus = NULL, indices = NULL){
    # Bootstrap sample, if indices are provided.
    if (!is.null(indices)){
        data <- data[indices, ]
    }
    # Make the names consistent.
    data$Y <- data[[Y]]
    data$Z <- data[[Z]]
    data$D <- data[[D]]
    # Define the estimating equations.
    total.formula <- formula(paste0("Y ~ 1 + Z + ", X_minus))
    firststage.formula <- formula(paste0(
        "D ~ 1 + Z * (", X_iv, ") + ", X_minus))
    secondstage.formula <- formula(paste0(
        "Y ~ 1 + Z * D + ", X_minus))
    # Counterfactual copies.
    copies <- treatment.copies(data)
    input_Z0.data <- copies$Z0
    input_Z1.data <- copies$Z1
    input_D0.data <- data
    input_D1.data <- data
    input_D0.data$D <- 0
    input_D1.data$D <- 1
    # Count estimated coefficients.
    count.coef <- 0
    # 0. Total effect regression.
    totaleffect.reg <- lm(total.formula, data = data)
    totaleffect.est <- mean(
        predict(totaleffect.reg, newdata = input_Z1.data) -
        predict(totaleffect.reg, newdata = input_Z0.data),
        na.rm = TRUE)
    count.coef <- count.coef + length(totaleffect.reg$coefficients)
    # 1. Mediator first stage.
    firststage.reg <- glm(
        firststage.formula,
        family = binomial(link = "probit"),
        data = data)
    pi_0.est <- predict(firststage.reg, newdata = input_Z0.data,
        type = "response")
    pi_1.est <- predict(firststage.reg, newdata = input_Z1.data,
        type = "response")
    pi.bar <- mean(pi_1.est - pi_0.est, na.rm = TRUE)
    count.coef <- count.coef + length(firststage.reg$coefficients)
    # 2. Naive outcome regression.
    secondstage.reg <- lm(secondstage.formula, data = data)
    count.coef <- count.coef + length(secondstage.reg$coefficients)
    mediator.effect.0 <- predict(
        secondstage.reg, newdata = input_D1.data) - predict(
            secondstage.reg, newdata = input_D0.data)
    # 3. Conventional plug-in AIE(0), and ADE(1) as the remainder.
    aie0.est <- mean(
        (pi_1.est - pi_0.est) * mediator.effect.0,
        na.rm = TRUE)
    ade1.est <- totaleffect.est - aie0.est
    # Return the estimates in the established order.
    output.list <- c(
        pi.bar,
        totaleffect.est,
        ade1.est,
        aie0.est,
        count.coef)
    return(output.list)
}

# Estimate AIE(0) from the Z = 0 restricted conditional-outcome curve,
# then estimate ADE(1) = ATE - AIE(0).
mediate.semiparametric <- function(Y, Z, D, X_iv, data,
        X_minus = NULL, indices = NULL,
        support.warning = FALSE, return.models = FALSE){
    ### The estimator follows the new identification argument:
    ## 1. Estimate the mediator propensity score on the full sample.
    ## 2. Estimate m(p;0,X^-) in the Z = 0 subsample.
    ## 3. Predict m at pi(1;X) and pi(0;X) for Z = 0 observations.
    ## 4. Average endpoint differences to estimate AIE(0).
    ## 5. Estimate ADE(1) = ATE - AIE(0).

    # Bootstrap sample, if indices are provided.
    if (!is.null(indices)){
        data <- data[indices, ]
    }
    # Make the names consistent.
    data$Y <- data[[Y]]
    data$Z <- data[[Z]]
    data$D <- data[[D]]
    # Define equations.
    total.formula <- formula(paste0("Y ~ 1 + Z + ", X_minus))
    firststage.formula <- formula(paste0(
        "D ~ 1 + Z * (", X_iv, ") + ", X_minus))
    outcome.formula <- formula(paste0(
        "Y ~ 1 + ", X_minus,
        " + poly(pi.est, degree = 5, raw = TRUE, simple = TRUE)"))
    # Counterfactual treatment copies.
    copies <- treatment.copies(data)
    input_Z0.data <- copies$Z0
    input_Z1.data <- copies$Z1
    # Count estimated coefficients.
    count.coef <- 0
    # 1. Estimate the total effect from quasi-random treatment assignment.
    totaleffect.reg <- lm(total.formula, data = data)
    totaleffect.est <- mean(
        predict(totaleffect.reg, newdata = input_Z1.data) -
            predict(totaleffect.reg, newdata = input_Z0.data),
                na.rm = TRUE)
    count.coef <- count.coef + length(totaleffect.reg$coefficients)
    # 2. Estimate the mediator propensity score.
    mte_firststage.reg <- glm(
        firststage.formula,
        family = binomial(link = "probit"),
        data = data)
    data$pi.est <- predict(mte_firststage.reg, type = "response")
    pi_0.est <- predict(mte_firststage.reg,
        newdata = input_Z0.data, type = "response")
    pi_1.est <- predict(mte_firststage.reg,
        newdata = input_Z1.data, type = "response")
    pi.bar <- mean(pi_1.est - pi_0.est, na.rm = TRUE)
    count.coef <- count.coef + length(mte_firststage.reg$coefficients)
    # 3. Restrict outcome-curve estimation and standardisation to Z = 0.
    outcome_Z0.data <- data[data$Z == 0, , drop = FALSE]
    endpoint_Z0.data <- outcome_Z0.data
    endpoint_Z0.data$pi_0.est <- pi_0.est[data$Z == 0]
    endpoint_Z0.data$pi_1.est <- pi_1.est[data$Z == 0]
    # Support warning.
    support.est <- support.summary(
        pi_realised_Z0 = outcome_Z0.data$pi.est,
        pi_0_Z0 = endpoint_Z0.data$pi_0.est,
        pi_1_Z0 = endpoint_Z0.data$pi_1.est)
    if (support.warning == TRUE &&
            support.est$share.interval.in.support < 1){
        warning(paste0(
            round(100 * (1 - support.est$share.interval.in.support), 2),
            "% of Z = 0 observations have at least one endpoint outside the ",
            "observed Z = 0 propensity-score range."))
    }
    # 4. Estimate m(p;0,X^-) in the Z = 0 subsample.
    mte_outcome_Z0.reg <- lm(outcome.formula, data = outcome_Z0.data)
    count.coef <- count.coef + length(mte_outcome_Z0.reg$coefficients)
    # Full-column-rank diagnostic for the maintained restricted model.
    outcome.design <- model.matrix(mte_outcome_Z0.reg)
    rank.diagnostic <- c(rank = qr(outcome.design)$rank,
        columns = ncol(outcome.design))
    # 5. Evaluate the same Z = 0 curve at both endpoint propensities.
    predict_pi0.data <- endpoint_Z0.data
    predict_pi1.data <- endpoint_Z0.data
    predict_pi0.data$pi.est <- endpoint_Z0.data$pi_0.est
    predict_pi1.data$pi.est <- endpoint_Z0.data$pi_1.est
    m_pi0.est <- predict(mte_outcome_Z0.reg, newdata = predict_pi0.data,
        type = "response")
    m_pi1.est <- predict(mte_outcome_Z0.reg, newdata = predict_pi1.data,
        type = "response")
    # 6. Estimate AIE(0), and ADE(1) as the remainder from the ATE.
    conditional.aie0.est <- m_pi1.est - m_pi0.est
    aie0.est <- mean(conditional.aie0.est, na.rm = TRUE)
    ade1.est <- totaleffect.est - aie0.est
    # Return models and diagnostics for one-shot inspection, but retain the
    # compact numeric output expected by the bootstrap and reporting code.
    if (return.models == TRUE){
        return(list(
            estimates = c(
                pi.bar,
                totaleffect.est,
                ade1.est,
                aie0.est,
                count.coef),
            support = support.est,
            rank = rank.diagnostic,
            firststage.reg = mte_firststage.reg,
            outcome.reg = mte_outcome_Z0.reg,
            conditional.aie = conditional.aie0.est))
    }
    output.list <- c(
        pi.bar,
        totaleffect.est,
        ade1.est,
        aie0.est,
        count.coef)
    return(output.list)
}

# Define a function to bootstrap the complete estimation procedure.
mediate.bootstrap <- function(Y, Z, D, X_iv, data,
        X_minus = NULL, type = "semi-parametric", boot.reps = 10){
    # Define an empty data frame.
    boot.data <- data.frame(matrix(ncol = 4, nrow = 0))
    names(boot.data) <- c(
        "First-stage", "ATE", "ADE", "AIE")
    for (i in seq_len(boot.reps)){
        if ((boot.reps >= 100) && ((100 * i / boot.reps) %% 5 == 0)){
            cat(paste0(i, " out of ", boot.reps, ", ",
                100 * (i / boot.reps), "% done.", "\n"))
        }
        boot.indices <- sample(
            seq_len(nrow(data)), nrow(data), replace = TRUE)
        if (type == "semi-parametric"){
            point.est <- mediate.semiparametric(
                Y, Z, D, X_iv, data,
                X_minus = X_minus, indices = boot.indices)
        }
        else if (type == "unadjusted"){
            point.est <- mediate.unadjusted(
                Y, Z, D, X_iv, data,
                X_minus = X_minus, indices = boot.indices)
        }
        else {
            stop(paste0("The type option only takes values of ",
                'c("semi-parametric", "unadjusted").'))
        }
        boot.data[i, ] <- point.est[1:4]
    }
    return(boot.data)
}

# Define a wrapper around the point estimates and bootstrap.
mediate.selection <- function(Y, Z, D, X_iv, X_minus, data,
        type = "semi-parametric", boot.reps = 10){
    # Calculate the point estimates.
    if (type == "semi-parametric"){
        point.est <- mediate.semiparametric(
            Y, Z, D, X_iv, data, X_minus = X_minus)
    }
    else if (type == "unadjusted"){
        point.est <- mediate.unadjusted(
            Y, Z, D, X_iv, data, X_minus = X_minus)
    }
    else {
        stop(paste0("The type option only takes values of ",
            'c("semi-parametric", "unadjusted").'))
    }
    count.coef <- point.est[5]
    # Calculate SEs by non-parametric bootstrap.
    if (!is.null(boot.reps)){
        if (boot.reps < 500){
            print(paste0("Attempting to bootstrap with fewer than 500 reps.",
                " Are you sure? This is likely not enough for convergence."))
        }
        point.boot <- mediate.bootstrap(
            Y = Y, Z = Z, D = D, X_minus = X_minus, X_iv = X_iv,
            data = data, type = type, boot.reps = boot.reps)
        point.se <- as.matrix(c(
            sd(100 * point.boot$"First-stage", na.rm = TRUE),
            sd(100 * point.boot$"ATE", na.rm = TRUE),
            sd(100 * point.boot$"ADE", na.rm = TRUE),
            sd(100 * point.boot$"AIE", na.rm = TRUE),
            sd(point.boot$"AIE" / point.boot$"ATE", na.rm = TRUE)))
    }
    else {
        point.se <- as.matrix(rep(NA, 5))
    }
    # Report output.
    point.est <- as.matrix(c(
        100 * point.est[1:4],
        point.est[4] / point.est[2]))
    tratio <- as.matrix(point.est / point.se)
    ptratio <- as.matrix(2 * pt(
        abs(tratio),
        df = nrow(data) - count.coef,
        lower.tail = FALSE))
    out <- list(
        coefficients = point.est,
        SE = point.se,
        tratio = tratio,
        ptratio = ptratio,
        type = type,
        variables = paste(Z, D, Y, sep = ", "),
        boot.reps = boot.reps)
    rownames(out$coefficients) <-
        c("First-stage", "ATE", "ADE(1)", "AIE(0)", "Proportion, AIE(0) / ATE")
    rownames(out$SE) <- rownames(out$coefficients)
    rownames(out$tratio) <- rownames(out$coefficients)
    rownames(out$ptratio) <- rownames(out$coefficients)
    class(out) <- "mediate.selection"
    return(out)
}

# Print applied to the function.
print.mediate.selection <- function(x, digits = 4, ...){
    cat("Treatment, Mediator, Outcome: \n")
    cat(x$variables)
    cat("\n")
    est <- cbind(x$coefficients, x$SE)
    colnames(est) <- c("Coefficients", "SE")
    cat(paste0("\n", x$type, " Estimates, With SEs from ",
        x$boot.reps, " bootstrap replications."))
    cat("\n\n")
    print.default(format(est, digits = digits), quote = FALSE)
}

# Apply the summary function, to get a presentable output.
summary.mediate.selection <- function(object, ...){
    TAP <- cbind(
        Estimate = coef(object),
        SE = object$SE,
        ptratio = object$ptratio)
    colnames(TAP) <- c("Estimate", "SE", "P")
    res <- list(variables = object$variables, coefficients = TAP)
    class(res) <- "summary.larf"
    return(res)
}

# Presentable summary, via printing.
print.summary.mediate.selection <- function(x, digits = 4, ...){
    cat("Treatment, Mediator, Outcome: \n")
    cat(x$variables)
    cat("\n")
    print.default(round(x$coefficients, digits = digits), quote = FALSE)
}


################################################################################
## One-shot tests of the updated estimators.

unadjusted.test <- mediate.unadjusted(
    Y = "Y_health", Z = "lottery_iv", D = "any_healthcare",
    X_iv = "initial_health_location",
    X_minus = controls.list,
    data = analysis.data)
print(unadjusted.test)

semiparametric.test <- mediate.semiparametric(
    Y = "Y_health", Z = "lottery_iv", D = "any_healthcare",
    X_iv = "initial_health_location",
    X_minus = controls.list,
    data = analysis.data,
    support.warning = TRUE,
    return.models = TRUE)
print(semiparametric.test$estimates)
print(semiparametric.test$support)
print(semiparametric.test$rank)
print(summary(semiparametric.test$firststage.reg))
print(summary(semiparametric.test$outcome.reg))


################################################################################
## Show the regular location is a strong IV for healthcare visits.

# Note values in usual health location:
initial_health_location.list <- c(
    "1. Private clinic",
    "2. Public clinic",
    "3. Hospital clinic",
    "4. Hospital A&E",
    "5. Urgent care",
    "6. Other clinic",
    "7. No regular")

# Show that the health location influences healthcare take-up
location.data <- analysis.data %>%
    mutate(initial_health_location_name =
        ifelse(initial_health_location == 1, initial_health_location.list[1],
        ifelse(initial_health_location == 2, initial_health_location.list[2],
        ifelse(initial_health_location == 3, initial_health_location.list[3],
        ifelse(initial_health_location == 4, initial_health_location.list[4],
        ifelse(initial_health_location == 5, initial_health_location.list[5],
        ifelse(initial_health_location == 6, initial_health_location.list[6],
        ifelse(initial_health_location == 7, initial_health_location.list[7],
            "WARNING")))))))) %>%
    group_by(initial_health_location_name) %>%
    summarise(
        any_healthcare_mean = mean(any_healthcare, na.rm = TRUE),
        any_healthcare_sd = sd(any_healthcare, na.rm = TRUE),
        count = n(),
        .groups = "drop") %>%
    mutate(any_healthcare_se = any_healthcare_sd / sqrt(count),
        any_healthcare_mean_lower = any_healthcare_mean -
            1.96 * any_healthcare_se,
        any_healthcare_mean_upper = any_healthcare_mean +
            1.96 * any_healthcare_se)

location.plot <- location.data %>%
    ggplot(aes(x = initial_health_location_name)) +
    geom_bar(aes(y = any_healthcare_mean, fill = initial_health_location_name),
        colour = "black", stat = "identity") +
    geom_errorbar(
        aes(ymin = any_healthcare_mean_lower,
            ymax = any_healthcare_mean_upper),
        width = 0.2, position = position_dodge(0.9)) +
    coord_flip() +
    theme_bw() +
    scale_x_discrete(name = "", limits = initial_health_location.list[7:1]) +
    scale_y_continuous(expand = c(0, 0),
        name = TeX(r"(Visited healthcare in following 12 months, $\Pr(\,D_i=1\mid X_i^{IV}\,)$)"),
        limits = c(0, 1.025), oob = scales::rescale_none,
        breaks = seq(0, 1, by = 0.1)) +
    ggtitle("Usual Healthcare Location") +
    theme(legend.position = "none",
        axis.text.y = element_text(hjust = 0),
        plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 2, 0, 0), "mm"))

ggsave(file.path(figures.folder, "location-effects.png"),
    plot = location.plot,
    units = "cm", width = fig.width, height = fig.height)

# First-stage relevance diagnostics.
print(mean(analysis.data$any_healthcare, na.rm = TRUE))
location.reg <- lm(
    formula(paste0("any_healthcare ~ 1 + factor(initial_health_location)",
        " + lottery_iv * ", nocontrols.list)),
    data = analysis.data)
print(summary(location.reg))
iv.list <- paste0("factor(initial_health_location)", 2:7)
print(car::linearHypothesis(location.reg, test = "F", iv.list))
print(100 * table(
    analysis.data$initial_health_location,
    exclude = NULL) / NROW(analysis.data))


################################################################################
## Examine associations between healthcare use and outcomes.

## Library for better colour choice and additional graphing.
library(ggthemes)
# remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)

outcome_name.list <- "0               1\n" %>%
    paste0(c("Health overall good?", "Happy overall?"))
healthcare.data <- analysis.data %>%
    group_by(any_healthcare) %>%
    summarise(
        health_mean = mean(Y_health),
        health_se = sd(Y_health) / sqrt(n() - 1),
        happy_mean = mean(Y_happy),
        happy_se = sd(Y_happy) / sqrt(n() - 1),
        .groups = "drop") %>%
    pivot_longer(cols = c(health_mean, health_se, happy_mean, happy_se),
        names_to = c("outcome_name", ".value"),
        names_pattern = "(.*)_(mean|se)") %>%
    mutate(ci_lower = mean - 1.96 * se,
        ci_upper = mean + 1.96 * se) %>%
    mutate(outcome_name =
        ifelse(outcome_name == "health", outcome_name.list[1],
            ifelse(outcome_name == "happy", outcome_name.list[2], "")))

#Plot it.
healthcare.plot <- healthcare.data %>%
    ggplot(aes(x = outcome_name, group = any_healthcare)) +
    geom_col_pattern(
        aes(pattern = any_healthcare, fill = outcome_name, y = mean),
        position = "dodge",
        pattern_angle = 45,
        pattern_density = 0.01,
        pattern_spacing = 0.1,
        pattern_fill = "black",
        colour = "black") +
    geom_errorbar(aes(x = outcome_name, ymin = ci_lower, ymax = ci_upper),
        size = 0.5, alpha = 0.5,
        stat = "identity", position = position_dodge(0.9), width = 1 / 3) +
    theme_bw() +
    scale_x_discrete(name = "Outcome, attending healthcare or not.",
        limits = outcome_name.list) +
    scale_fill_manual("", values = colour.list[c(2, 2)]) +
    scale_y_continuous(expand = c(0, 0), name = "",
        limits = c(0.275, 0.825), oob = scales::rescale_none,
        breaks = seq(0, 1, by = 0.1)) +
    theme(legend.position = "none",
        plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 0, 0, 0), "mm"))
# Save this plot for the presentation.
ggsave(file.path(figures.folder, "healthcare-ols.png"),
    plot = healthcare.plot,
    units = "cm", width = 0.8 * fig.width, height = 0.65 * fig.height)

healthcare.presentation.plot <- healthcare.data %>%
    filter(outcome_name != "0               1\nHealth overall good?") %>%
    ggplot(aes(x = outcome_name, group = any_healthcare)) +
    geom_col_pattern(
        aes(pattern = any_healthcare, fill = outcome_name, y = mean),
        position = "dodge",
        pattern_angle = 45,
        pattern_density = 0.01,
        pattern_spacing = 0.1,
        pattern_fill = "black",
        colour = "black") +
    geom_errorbar(aes(x = outcome_name, ymin = ci_lower, ymax = ci_upper),
        size = 0.5, alpha = 0.5,
        stat = "identity", position = position_dodge(0.9), width = 1 / 3) +
    theme_bw() +
    scale_x_discrete(name = "Outcome, attending healthcare or not.") +
    scale_fill_manual("", values = colour.list[c(2, 2)]) +
    scale_y_continuous(expand = c(0, 0), name = "",
        limits = c(0.275, 0.825), oob = scales::rescale_none,
        breaks = seq(0, 1, by = 0.1)) +
    theme(legend.position = "none",
        plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 0, 0, 0), "mm"))

ggsave(file.path(figures.folder, "healthcare-ols-presentation.png"),
    plot = healthcare.presentation.plot,
    units = "cm", width = 0.5 * fig.width, height = 0.65 * fig.height)

# OLS correlations.
health.reg <- lm(
    formula(paste0("Y_health ~ 1 + any_healthcare + ", nocontrols.list)),
    data = analysis.data)
print(summary(health.reg))
happy.reg <- lm(
    formula(paste0("Y_happy ~ 1 + any_healthcare + ", nocontrols.list)),
    data = analysis.data)
print(summary(happy.reg))

health_controls.reg <- lm(
    formula(paste0("Y_health ~ 1 + any_healthcare + ", controls.list)),
    data = analysis.data)
print(summary(health_controls.reg))
happy_controls.reg <- lm(
    formula(paste0("Y_happy ~ 1 + any_healthcare + ", controls.list)),
    data = analysis.data)
print(summary(happy_controls.reg))

# IV associations.
health.iv <- feols(
    formula(paste0("Y_health ~ 1 + ", nocontrols.list,
        " | any_healthcare ~ factor(initial_health_location)")),
    data = analysis.data)
print(summary(health.iv))
happy.iv <- feols(
    formula(paste0("Y_happy ~ 1 + ", nocontrols.list,
        " | any_healthcare ~ factor(initial_health_location)")),
    data = analysis.data)
print(summary(happy.iv))

# Implied average mediator effect among mediator compliers.
mediate.est <- mediate.selection(
    Y = "Y_health", Z = "lottery_iv", D = "any_healthcare",
    X_iv = "initial_health_location",
    X_minus = nocontrols.list,
    boot.reps = NULL,
    type = "semi-parametric",
    data = analysis.data)
print(coeftable(mediate.est)["AIE(0)", "Estimate"] /
    coeftable(mediate.est)["First-stage", "Estimate"])
mediate.est <- mediate.selection(
    Y = "Y_happy", Z = "lottery_iv", D = "any_healthcare",
    X_iv = "initial_health_location",
    X_minus = nocontrols.list,
    boot.reps = NULL,
    type = "semi-parametric",
    data = analysis.data)
print(coeftable(mediate.est)["AIE(0)", "Estimate"] /
    coeftable(mediate.est)["First-stage", "Estimate"])


################################################################################
## Estimate the CM effects with the conventional and MTE methods.

# State how many bootstrap replications are needed.
boot.reps <- 10^3

# Extract effects for plotting.
effects.extract <- function(mediate.est, model.name){
    reg.summary <- summary(mediate.est)
    total.est <- coeftable(reg.summary)["ATE", "Estimate"]
    total.ci.upper <- total.est + 1.96 * coeftable(reg.summary)["ATE", "SE"]
    total.ci.lower <- total.est - 1.96 * coeftable(reg.summary)["ATE", "SE"]
    direct.est <- coeftable(reg.summary)["ADE(1)", "Estimate"]
    direct.ci.upper <- direct.est + 1.96 * coeftable(reg.summary)["ADE(1)", "SE"]
    direct.ci.lower <- direct.est - 1.96 * coeftable(reg.summary)["ADE(1)", "SE"]
    indirect.est <- coeftable(reg.summary)["AIE(0)", "Estimate"]
    indirect.ci.upper <- indirect.est + 1.96 * coeftable(reg.summary)["AIE(0)", "SE"]
    indirect.ci.lower <- indirect.est - 1.96 * coeftable(reg.summary)["AIE(0)", "SE"]
    permediated.est <- coeftable(reg.summary)[
        "Proportion, AIE(0) / ATE", "Estimate"]
    permediated.ci.upper <- permediated.est + 1.96 * coeftable(reg.summary)[
        "Proportion, AIE(0) / ATE", "SE"]
    permediated.ci.lower <- permediated.est - 1.96 * coeftable(reg.summary)[
        "Proportion, AIE(0) / ATE", "SE"]
    data.return <- data.frame(
        effect = c("Total", "Direct", "Indirect", "Percent Mediated"),
        pointest = c(total.est, direct.est, indirect.est, permediated.est),
        upperest = c(total.ci.upper, direct.ci.upper,
            indirect.ci.upper, permediated.ci.upper),
        lowerest = c(total.ci.lower, direct.ci.lower,
            indirect.ci.lower, permediated.ci.lower)) %>%
        mutate(model = model.name)
    return(data.return)
}

# Panel A: Self-reported health.
unadjusted.health.est <- mediate.selection(
    Y = "Y_health", Z = "lottery_iv", D = "any_healthcare",
    X_iv = "initial_health_location", X_minus = controls.list,
    type = "unadjusted", data = analysis.data, boot.reps = boot.reps)
print(unadjusted.health.est)
semiparametric.health.est <- mediate.selection(
    Y = "Y_health", Z = "lottery_iv", D = "any_healthcare",
    X_iv = "initial_health_location", X_minus = controls.list,
    type = "semi-parametric", data = analysis.data, boot.reps = boot.reps)
print(semiparametric.health.est)

# Extract for the LaTeX table.
panelA.data <- data.frame(
    unadjusted_point = coef(summary(unadjusted.health.est))[, "Estimate"],
    unadjusted_se = coef(summary(unadjusted.health.est))[, "SE"],
    semiparametric_point = coef(summary(semiparametric.health.est))[, "Estimate"],
    semiparametric_se = coef(summary(semiparametric.health.est))[, "SE"])

Y_health.data <- rbind(
    effects.extract(unadjusted.health.est, "Conventional"),
    effects.extract(semiparametric.health.est, "Semi-parametric MTE"))

panelA.table <- panelA.data %>%
    signif(digits.no) %>%
    format(scientific = FALSE)
panelA.table$unadjusted_se <- paste0("(", panelA.table$unadjusted_se, ")")
panelA.table$semiparametric_se <- paste0(
    "(", panelA.table$semiparametric_se, ")")
panelA.table <- data.frame(t(panelA.table))
panelA.table <- cbind(model = c(
    "Conventional CM", "", "Semi-parametric MTE", ""), panelA.table)

# Save the LaTeX table
panelA.table %>%
    xtable() %>%
    print(
        digits = digits.no,
        sanitize.colnames.function = identity,
        sanitize.text.function = identity,
        NA.string = " ",
        include.colnames = FALSE,
        include.rownames = FALSE,
        only.contents = TRUE,
        hline.after = NULL,
        format.args = list(big.mark = ","),
        file = file.path(tables.folder, "cm-oregon-health.tex"))

# Panel B: Self-reported happiness.
unadjusted.happy.est <- mediate.selection(
    Y = "Y_happy", Z = "lottery_iv", D = "any_healthcare",
    X_iv = "initial_health_location", X_minus = controls.list,
    type = "unadjusted", data = analysis.data, boot.reps = boot.reps)
print(unadjusted.happy.est)
semiparametric.happy.est <- mediate.selection(
    Y = "Y_happy", Z = "lottery_iv", D = "any_healthcare",
    X_iv = "initial_health_location", X_minus = controls.list,
    type = "semi-parametric", data = analysis.data, boot.reps = boot.reps)
print(semiparametric.happy.est)

panelB.data <- data.frame(
    unadjusted_point = coef(summary(unadjusted.happy.est))[, "Estimate"],
    unadjusted_se = coef(summary(unadjusted.happy.est))[, "SE"],
    semiparametric_point = coef(summary(semiparametric.happy.est))[, "Estimate"],
    semiparametric_se = coef(summary(semiparametric.happy.est))[, "SE"])

Y_happy.data <- rbind(
    effects.extract(unadjusted.happy.est, "Conventional"),
    effects.extract(semiparametric.happy.est, "Semi-parametric MTE"))

# Extract for the LaTeX table.
panelB.table <- panelB.data %>%
    signif(digits.no) %>%
    format(scientific = FALSE)
panelB.table$unadjusted_se <- paste0("(", panelB.table$unadjusted_se, ")")
panelB.table$semiparametric_se <- paste0(
    "(", panelB.table$semiparametric_se, ")")
panelB.table <- data.frame(t(panelB.table))
panelB.table <- cbind(model = c(
    "Conventional CM", "", "Semi-parametric MTE", ""), panelB.table)

print(panelB.table)

panelB.table %>%
    xtable() %>%
    print(
        digits = digits.no,
        sanitize.colnames.function = identity,
        sanitize.text.function = identity,
        NA.string = " ",
        include.colnames = FALSE,
        include.rownames = FALSE,
        only.contents = TRUE,
        hline.after = NULL,
        format.args = list(big.mark = ","),
        file = file.path(tables.folder, "cm-oregon-happy.tex"))


################################################################################
## Plot the results.

health_mediation.plot <- Y_health.data %>%
    filter(effect != "Percent Mediated") %>%
    ggplot(aes(
        fill = factor(effect, levels = c("Total", "Direct", "Indirect")),
        x = model)) +
    geom_bar(aes(y = pointest),
        stat = "identity", position = "dodge", colour = "black") +
    geom_errorbar(aes(ymin = lowerest, ymax = upperest),
        size = 2 / 3,
        stat = "identity", position = position_dodge(0.9), width = 1 / 3) +
    theme_bw() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_discrete(
        name = "",
        limits = c("Total", "Direct", "Indirect"),
        type = colour.list[c(3, 1, 2)]) +
    scale_x_discrete(
        name = "",
        limits = c("Conventional", "Semi-parametric MTE")) +
    scale_y_continuous(expand = c(0, 0, 0.01, 0),
        limits = c(-0.5, 9.75), breaks = seq(-10, 10, by = 1),
        oob = scales::rescale_none, name = "") +
    ggtitle("Estimate, percent effect on subjective health") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 3, 0.25, 0), "mm"),
        legend.position = c(0.8, 1.025),
        legend.direction = "horizontal")

ggsave(file.path(figures.folder, "mediation-health.png"),
    plot = health_mediation.plot,
    units = "cm",
    width = presentation.width, height = presentation.height)

happy_mediation.plot <- Y_happy.data %>%
    filter(effect != "Percent Mediated") %>%
    ggplot(aes(
        fill = factor(effect, levels = c("Total", "Direct", "Indirect")),
        x = model)) +
    geom_bar(aes(y = pointest),
        stat = "identity", position = "dodge", colour = "black") +
    geom_errorbar(aes(ymin = lowerest, ymax = upperest),
        size = 2 / 3,
        stat = "identity", position = position_dodge(0.9), width = 1 / 3) +
    theme_bw() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_discrete(
        name = "",
        limits = c("Total", "Direct", "Indirect"),
        type = colour.list[c(3, 1, 2)]) +
    scale_x_discrete(
        name = "",
        limits = c("Conventional", "Semi-parametric MTE")) +
    scale_y_continuous(expand = c(0, 0, 0.01, 0),
        limits = c(-0.5, 9.75), breaks = seq(-10, 10, by = 1),
        oob = scales::rescale_none, name = "") +
    ggtitle("Estimate, percent effect on subjective well-being") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 3, 0.25, 0), "mm"),
        legend.position = c(0.8, 1.025),
        legend.direction = "horizontal")

ggsave(file.path(figures.folder, "mediation-happy.png"),
    plot = happy_mediation.plot,
    units = "cm",
    width = presentation.width, height = presentation.height)


################################################################################
## Placeholder plots showing only conventional estimates.

health_mediation.placeholder <- Y_health.data %>%
    filter(effect != "Percent Mediated") %>%
    mutate(
        pointest = ifelse(model == "Conventional", pointest, NA),
        lowerest = ifelse(model == "Conventional", lowerest, NA),
        upperest = ifelse(model == "Conventional", upperest, NA)) %>%
    ggplot(aes(
        fill = factor(effect, levels = c("Total", "Direct", "Indirect")),
        x = model)) +
    geom_bar(aes(y = pointest),
        stat = "identity", position = "dodge", colour = "black",
        na.rm = TRUE) +
    geom_errorbar(aes(ymin = lowerest, ymax = upperest),
        size = 2 / 3,
        stat = "identity", position = position_dodge(0.9), width = 1 / 3,
        na.rm = TRUE) +
    theme_bw() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_discrete(
        name = "",
        limits = c("Total", "Direct", "Indirect"),
        type = colour.list[c(3, 1, 2)]) +
    scale_x_discrete(
        name = "",
        limits = c("Conventional", "Semi-parametric MTE")) +
    scale_y_continuous(expand = c(0, 0, 0.01, 0),
        limits = c(-0.5, 9.75), breaks = seq(-10, 10, by = 1),
        oob = scales::rescale_none, name = "") +
    ggtitle("Estimate, percent effect on subjective health") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 3, 0.25, 0), "mm"),
        legend.position = c(0.8, 1.025),
        legend.direction = "horizontal")

ggsave(file.path(figures.folder, "mediation-health-placeholder.png"),
    plot = health_mediation.placeholder,
    units = "cm",
    width = presentation.width, height = presentation.height)

happy_mediation.placeholder <- Y_happy.data %>%
    filter(effect != "Percent Mediated") %>%
    mutate(
        pointest = ifelse(model == "Conventional", pointest, NA),
        lowerest = ifelse(model == "Conventional", lowerest, NA),
        upperest = ifelse(model == "Conventional", upperest, NA)) %>%
    ggplot(aes(
        fill = factor(effect, levels = c("Total", "Direct", "Indirect")),
        x = model)) +
    geom_bar(aes(y = pointest),
        stat = "identity", position = "dodge", colour = "black",
        na.rm = TRUE) +
    geom_errorbar(aes(ymin = lowerest, ymax = upperest),
        size = 2 / 3,
        stat = "identity", position = position_dodge(0.9), width = 1 / 3,
        na.rm = TRUE) +
    theme_bw() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_discrete(
        name = "",
        limits = c("Total", "Direct", "Indirect"),
        type = colour.list[c(3, 1, 2)]) +
    scale_x_discrete(
        name = "",
        limits = c("Conventional", "Semi-parametric MTE")) +
    scale_y_continuous(expand = c(0, 0, 0.01, 0),
        limits = c(-0.5, 9.75), breaks = seq(-10, 10, by = 1),
        oob = scales::rescale_none, name = "") +
    ggtitle("Estimate, percent effect on subjective well-being") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 3, 0.25, 0), "mm"),
        legend.position = c(0.8, 1.025),
        legend.direction = "horizontal")

ggsave(file.path(figures.folder, "mediation-happy-placeholder.png"),
    plot = happy_mediation.placeholder,
    units = "cm",
    width = presentation.width, height = presentation.height)
