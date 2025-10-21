#!/usr/bin/R
## Senan Hogan-Hennessy, 11 July 2025.
## Script to to extract relevant data from Oregon Health Insurance rep package.
print(Sys.time())
set.seed(47)

## Packages:
# functions for data manipulation and visualisation
library(tidyverse)
# Functions for bootstrapping.
library(boot)
library(fixest)
# Package forsemi-parametric MTE by splines.
library(mgcv)
library(splines)
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

# List all the observed control variables
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
## Define the functions to use.

# Estimate the values, given a first and second-stages
estimated.values <- function(firststage.reg, secondstage.reg, totaleffect.reg,
    data, complier.adjustment = NULL){
    ### Inputs:
    ## data, a data frame simulated from above.
    input_Z0.data <- data
    input_Z1.data <- data
    input_D0.data <- data
    input_D1.data <- data
    input_Z0.data$Z <- 0
    input_Z1.data$Z <- 1
    input_D0.data$D <- 0
    input_D1.data$D <- 1
    # calculate the first-stage by prediction
    firststage.est <- predict(
        firststage.reg, newdata = input_Z1.data, type = "response") - predict(
            firststage.reg, newdata = input_Z0.data, type = "response")
    # Calculate the total effect estimate by prediction.
    totaleffect.est <- predict(
        totaleffect.reg, newdata = input_Z1.data) - predict(
            totaleffect.reg, newdata = input_Z0.data)
    # calculate the second-stage direct effect
    direct.est <- predict(
        secondstage.reg, newdata = input_Z1.data) - predict(
            secondstage.reg, newdata = input_Z0.data)
    # calculate the second-stage (controlled) indirect effect
    indirect.est <- predict(
        secondstage.reg, newdata = input_D1.data) - predict(
            secondstage.reg, newdata = input_D0.data)
    # Add the Kline Walters (2019) IV-type complier adjustment (provided externally).
    if (!is.null(complier.adjustment)) {
        indirect.est <- indirect.est + complier.adjustment
    }
    # Return the mean estimates.
    output.list <- list(
        "first-stage"     = mean(firststage.est, na.rm = TRUE),
        "total-effect"    = mean(totaleffect.est, na.rm = TRUE),
        "direct-effect"   = mean(direct.est, na.rm = TRUE),
        "indirect-effect" = mean(firststage.est * indirect.est, na.rm = TRUE))
    # Return the output.list
    return(output.list)
}

# Define a function to estimate mediation, in two-stages (no adjustment).
mediate.unadjusted <- function(Y, Z, D, X_iv, data,
    X_minus = NULL, control_iv = NULL, indices = NULL){
    # Bootstrap sample, if indices provided.
    if (!is.null(indices)){
        data <- data[indices, ]
    }
    # Make the names consistent
    data$Y <- data[[Y]]
    data$Z <- data[[Z]]
    data$D <- data[[D]]
    total.formula <- formula(paste0("Y ~ 1 + Z + ", X_minus, " + ", X_iv))
    firststage.formula <- formula(paste0("D ~ 1 + Z * (", X_iv, ") +", X_minus))
    secondstage.formula <- formula(paste0("Y ~ 1 + Z * D + ", X_minus, " + ", X_iv))
    # Start counting coefficients estimated.
    count.coef <- 0
    # 0. Total effect regression.
    totaleffect.reg <- lm(total.formula, data = data)
    count.coef <- count.coef + length(totaleffect.reg$coefficients)
    # 1. Regular first-stage (well identified).
    firststage.reg <- lm(firststage.formula, data = data)
    count.coef <- count.coef + length(firststage.reg$coefficients)
    # 2. Estimate second-stage (naive case has no MTEs).
    unadjusted_secondstage.reg <- lm(secondstage.formula, data = data)
    count.coef <- count.coef + length(unadjusted_secondstage.reg$coefficients)
    # Compile the estimates.
    unadjusted.est <- estimated.values(
        firststage.reg, unadjusted_secondstage.reg,
        totaleffect.reg, data,
        complier.adjustment = NULL)
    # Return the off-setting estimates.
    output.list <- c(
        unadjusted.est$"first-stage",
        unadjusted.est$"total-effect",
        unadjusted.est$"direct-effect",
        unadjusted.est$"indirect-effect",
        count.coef)
    return(output.list)
}

#! Test it out, with the specification I want.
mediate.est <- mediate.unadjusted(
    Y = "Y_health", Z = "lottery_iv", D = "any_insurance",
    X_iv = "initial_health_location",
    X_minus = controls.list,
    analysis.data)# %>% sample_frac(prop = 1, replace = TRUE))
print(mediate.est)

## Define a function to Heckman selection correct mediation est, in two-stages.
# First the Mills ratio selection correction term(s)
mills.ratio <- function(pi.est){
    # Inv Mills ratio, taking as input the estimated mediator propensity.
    return(dnorm(qnorm(pi.est)) / pnorm(qnorm(pi.est)))
    }
# The two-stage mediation model.
mediate.heckit <- function(Y, Z, D, X_iv, data,
    X_minus = NULL, indices = NULL){
    # Bootstrap sample, if indices provided.
    if (!is.null(indices)){
        data <- data[indices, ]
    }
    # Make the names consistent
    data$Y <- data[[Y]]
    data$Z <- data[[Z]]
    data$D <- data[[D]]
    total.formula <- formula(paste0("Y ~ 1 + Z + ", X_minus))
    firststage.formula <- formula(paste0("D ~ 1 + Z * (", X_iv, ") +", X_minus))
    secondstage.formula <- formula(paste0("Y ~ 1 + Z * D + ", X_minus,
        " + lambda_0 + lambda_1"))
    # Start counting coefficients estimated.
    count.coef <- 0
    # 0. Total effect regression.
    totaleffect.reg <- lm(total.formula, data = data)
    count.coef <- count.coef + length(totaleffect.reg$coefficients)
    # 1. Probit first-stage (well identified).
    heckit_firststage.reg <- glm(firststage.formula,
        family = binomial(link = "probit"),
        data = data)
    count.coef <- count.coef + length(heckit_firststage.reg$coefficients)
    # 2. Define the MTEs --- for assumed N(0,1) dist.
    pi.est <- predict(heckit_firststage.reg, type = "response")
    data$lambda_0 <- (1 - data$D) * mills.ratio(pi.est) * (
        - pi.est / (1 - pi.est))
    data$lambda_1 <- data$D * mills.ratio(pi.est)
    # 3. Estimate second-stage, including the MTEs.
    heckit_secondstage.reg <- lm(secondstage.formula, data = data)
    count.coef <- count.coef + length(heckit_secondstage.reg$coefficients)
    # Compensate complier difference in AIE, by Kline Walters (2019) IV-type adjustment.
    input_Z0.data <- data
    input_Z1.data <- data
    input_Z0.data$Z <- 0
    input_Z1.data$Z <- 1
    pi_0.est <- predict(heckit_firststage.reg, newdata = input_Z0.data)
    pi_1.est <- predict(heckit_firststage.reg, newdata = input_Z1.data)
    Gamma.big <-  (pi_1.est * mills.ratio(pi_1.est)
        - pi_0.est * mills.ratio(pi_0.est)) / (pi_1.est - pi_0.est)
    rho_0 <- as.numeric(coef(heckit_secondstage.reg)["lambda_0"])
    rho_1 <- as.numeric(coef(heckit_secondstage.reg)["lambda_1"])
    complier.adjustment <- (rho_1 - rho_0) * Gamma.big
    # Compile the estimates.
    heckit.est <- estimated.values(
        heckit_firststage.reg, heckit_secondstage.reg,
        totaleffect.reg, data,
        complier.adjustment = complier.adjustment)
    # Return the off-setting estimates.
    output.list <- c(
        heckit.est$"first-stage",
        heckit.est$"total-effect",
        heckit.est$"direct-effect",
        heckit.est$"indirect-effect",
        count.coef)
    return(output.list)
}

#! Test it out, with the specification I want.
mediate.est <- mediate.heckit(
    Y = "Y_health", Z = "lottery_iv", D = "any_healthcare",
    X_iv = "initial_health_location",
    X_minus = controls.list,
    data = analysis.data)# %>% sample_frac(prop = 1, replace = TRUE))
print(mediate.est)

# Define a function to two-stage semi-parametric MTE for CM effects.
mediate.semiparametric <- function(Y, Z, D, X_iv, data,
    X_minus = NULL, indices = NULL){
    # Bootstrap sample, if indices provided.
    if (!is.null(indices)){
        data <- data[indices, ]
    }
    # Make the names consistent
    data$Y <- data[[Y]]
    data$Z <- data[[Z]]
    data$D <- data[[D]]
    total.formula <- formula(paste0("Y ~ 1 + Z +", X_minus))
    firststage.formula <- formula(paste0("D ~ 1 + Z * (", X_iv, ") +", X_minus))
    secondstage.formula <- formula(paste0("Y ~ 1 + Z + ", X_minus,
        " + s(pi.est, bs = 'cr')"))
    # Get relevant columns for imputation.
    input_Z0.data <- data
    input_Z1.data <- data
    input_D0.data <- data
    input_D1.data <- data
    input_Z0.data$Z <- 0
    input_Z1.data$Z <- 1
    input_D0.data$D <- 0
    input_D1.data$D <- 1
    # Start counting coefficients estimated.
    count.coef <- 0
    # 1. Total effect regression.
    totaleffect.reg <- lm(total.formula, data = data)
    totaleffect.est <- mean(predict(
        totaleffect.reg, newdata = input_Z1.data) - predict(
            totaleffect.reg, newdata = input_Z0.data))
    count.coef <- count.coef + length(totaleffect.reg$coefficients)
    # 2. Semi-parametric first-stage
    mte_firststage.reg <- glm(firststage.formula,
        family = binomial(link = "probit"),
        data = data)
    data$pi.est <- predict(mte_firststage.reg, type = "response")
    pi_0.est <- predict(mte_firststage.reg, newdata = input_Z0.data, type = "response")
    pi_1.est <- predict(mte_firststage.reg, newdata = input_Z1.data, type = "response")
    pi.bar <- mean(pi_1.est - pi_0.est)
    count.coef <- count.coef + length(mte_firststage.reg$coefficients)
    # Calculate the levels of pi, accounting for few values in the IV.
    distinct_mte.values <- min(
        length(unique(data$pi.est)) - 2, as.integer(nrow(data) / 2000))
    # 3. Semi-parametric series estimation of the second-stage.
    mte_secondstage_D0.reg <- gam(secondstage.formula,
        method = "REML", data = data, subset = (D == 0))
    mte_secondstage_D1.reg <- gam(secondstage.formula,
        method = "REML", data = data, subset = (D == 1))
    count.coef <- count.coef + length(mte_secondstage_D0.reg$coefficients)
    count.coef <- count.coef + length(mte_secondstage_D1.reg$coefficients)
    # 4. Compose the CM effects from this object.
    D_0 <- 1 - mean(data$D)
    D_1 <- 1 - D_0
    # 4.1 ADE point estimate, from the MTE model.
    gammma.est <- coef(mte_secondstage_D0.reg)["Z"]
    delta_plus.est <- coef(mte_secondstage_D1.reg)["Z"]
    ade.est <- as.numeric(D_0 * gammma.est + D_1 * delta_plus.est)
    # 4.2 AIE by using ADE estimate, relative to ATE.
    # (Avoiding semi-parametric extrapolation, see notes on ATE comparison)
    delta.est <- delta_plus.est - gammma.est
    ade_Z0.est <- gammma.est + delta.est * mean(
        data$D[data$Z == 0])
    ade_Z1.est <- gammma.est + delta.est * mean(
        data$D[data$Z == 1])
    aie.est <- (totaleffect.est - mean(
        (1 - data$Z) * ade_Z1.est + (data$Z) * ade_Z0.est))
    # Return the estimates.
    output.list <- c(
        pi.bar,
        totaleffect.est,
        ade.est,
        aie.est,
        count.coef)
    return(output.list)
}

#! Test it out, with the specification I want.
mediate.est <- mediate.semiparametric(
    Y = "Y_health", Z = "lottery_iv", D = "any_healthcare",
    X_iv = "initial_health_location",
    X_minus = controls.list,
    data = analysis.data)# %>% sample_frac(prop = 1, replace = TRUE))
print(mediate.est)

# Define a function to bootstrap.
mediate.bootstrap <- function(Y, Z, D, X_iv, data,
        X_minus = NULL, type = "parametric", boot.reps = 10){
    # Define an empty data.frame.
    boot.data <- data.frame(matrix(ncol = 4, nrow = 0))
    names(boot.data) <- c(
        "First-stage", "ATE", "ADE", "AIE")
    j <- 1
    for (i in 1:boot.reps){
        if ((boot.reps >= 100) & ((100 * i / boot.reps) %% 5 == 0)){
            cat(paste0(i, " out of ", boot.reps, ", ", 100 * (i / boot.reps),
                "% done.", "\n"))
        }
        boot.indices <- sample(1:nrow(data), nrow(data), replace = TRUE)
        if (type == "parametric"){
            point.est <- mediate.heckit(Y, Z, D, X_iv, data,
                X_minus = X_minus, indices = boot.indices)
            }
            else if (type == "semi-parametric"){
                point.est <- mediate.semiparametric(
                    Y, Z, D, X_iv, data,
                    X_minus = X_minus, indices = boot.indices)
            }
            else if (type == "unadjusted"){
                point.est <- mediate.unadjusted(Y, Z, D, X_iv, data,
                    X_minus = X_minus, indices = boot.indices)
            }
            else {
                stop(paste0("The type option only takes values of ",
                    'c("parametric", "semi-parametric", "unadjusted").'))
            }
        boot.data[i, ] <- point.est
    }
    return(boot.data)
}

#! Test it out.
mediate.boot <- mediate.bootstrap(
    Y = "Y_health", Z = "lottery_iv", D = "any_healthcare",
    X_iv = "initial_health_location",
    X_minus = controls.list,
    data = analysis.data, type = "unadjusted", boot.reps = 10)
print(mediate.boot)

## Define a function to wrap around all the others.
mediate.selection <- function(Y, Z, D, X_iv, X_minus, data,
    type = "parametric", boot.reps = 10){
    # Calculate the point estimates.
    if (type == "parametric"){
        point.est <- mediate.heckit(Y, Z, D, X_iv, data, X_minus = X_minus)
    }
    else if (type == "semi-parametric"){
        point.est <- mediate.semiparametric(Y, Z, D, X_iv, data, X_minus = X_minus)
    }
    else if (type == "unadjusted"){
        point.est <- mediate.unadjusted(Y, Z, D, X_iv, data, X_minus = X_minus)
    }
    else {
        stop(paste0("The type option only takes values of ",
            'c("parametric", "semi-parametric", "unadjusted").'))
    }
    count.coef <- point.est[5]
    # Calculate the SEs by a non-parametric bootstrap.
    if (!is.null(boot.reps)){
        if (boot.reps < 500){
            print(paste0("Attempting to bootstrap with fewer than 500 reps.",
                "  Are you sure?  This is likely not enough for convergence."))
        }
        point.boot <- mediate.bootstrap(
            Y = Y, Z = Z, D = D, X_minus = X_minus, X_iv = X_iv,
            data = data, type = type, boot.reps = boot.reps)
        point.se <- as.matrix(c(
            sd(100 * point.boot$"First-stage"),
            sd(100 * point.boot$"ATE"),
            sd(100 * point.boot$"ADE"),
            sd(100 * point.boot$"AIE"),
            sd(point.boot$"AIE" / point.boot$"ATE")))
    }
    else {
        point.se <- as.matrix(c(
            NA,
            NA,
            NA,
            NA,
            NA))
    }
    # Report output
    point.est <- as.matrix(c(100 * point.est[1:4], point.est[4] / point.est[2]))
    tratio <- as.matrix(point.est / point.se)
    ptratio <- as.matrix(2 * pt(abs(tratio),
        df = nrow(data) - count.coef, lower.tail = FALSE))
    # Preapred object to putput.
    out <- list(
        coefficients = point.est,
        SE = point.se,
        tratio = tratio,
        ptratio = ptratio,
        type = type,
        variables = paste(Z, D, Y, sep = ", "),
        boot.reps = boot.reps)
    rownames(out$coefficients) <-
        c("First-stage", "ATE", "ADE", "AIE", "Proportion, AIE / ATE")
    rownames(out$SE)      <-rownames(out$coefficients)
    rownames(out$tratio)  <-rownames(out$coefficients)
    rownames(out$ptratio) <-rownames(out$coefficients)
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
    cat(paste0("\n", x$type, " Estimates, With SEs from ", x$boot.reps, " bootstrap replications."))
    cat("\n\n")
    print.default(format(est, digits = digits), quote = FALSE)
}

# Apply the summary function, to get a presentable output.
summary.mediate.selection <- function(object, ...){
    TAP <- cbind(
        Estimate = coef(object),
        SE = object$SE,
        ptratio = object$ptratio )
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
        count = n()) %>%
    ungroup() %>%
    mutate(any_healthcare_se = any_healthcare_sd / (count^(0.5)),
        any_healthcare_mean_lower = any_healthcare_mean - 1.96 * any_healthcare_se,
        any_healthcare_mean_upper = any_healthcare_mean + 1.96 * any_healthcare_se)

# Draw the horizontal bar chart
location.plot <- location.data %>%
    ggplot(aes(x = initial_health_location_name)) +
    # Bars of the mean
    geom_bar(aes(y = any_healthcare_mean, fill = initial_health_location_name),
        colour = "black", stat = "identity") +
    # Error bars
    geom_errorbar(
        aes(ymin = any_healthcare_mean_lower, ymax = any_healthcare_mean_upper),
        width = 0.2, position = position_dodge(0.9)) +
    # Make horizontal and format.
    coord_flip() +
    theme_bw() +
    scale_x_discrete(name = "", limits = initial_health_location.list[7:1]) +
    scale_y_continuous(expand = c(0, 0),
        name = TeX(r"(Visited healthcare in following 12 months, $\Pr( \,D_i = 1 \, | \, X^{IV}_i \,)$)"),
        limits = c(0, 1.025), oob = scales::rescale_none,
        breaks = seq(0, 1, by = 0.1)) +
    ggtitle("Usual Healthcare Location") +
    theme(legend.position = "none",
        axis.text.y = element_text(hjust = 0),
        plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 2, 0, 0), "mm"))
# Save the plot.
ggsave(file.path(figures.folder, "location-effects.png"),
    plot = location.plot,
    units = "cm", width = fig.width, height = fig.height)

# Show the OLS correlation between D (mediator) and Y (outcome.)
library(margins)
health.reg <- lm(
    formula(paste0("Y_health ~ 1 + any_healthcare + ", nocontrols.list)),
    data = analysis.data)
print(summary(health.reg))
happy.reg <- lm(
    formula(paste0("Y_happy ~ 1 + any_healthcare + ", nocontrols.list)),
    data = analysis.data)
print(summary(happy.reg))
# Show the IV effect between D (mediator) and Y (outcome.)
library(fixest)
health.iv <- feols(
    formula(paste0("Y_health ~ 1 + ", nocontrols.list,
        "| any_healthcare ~ factor(initial_health_location)")),
    data = analysis.data)
print(summary(health.iv))
happy.iv <- feols(
    formula(paste0("Y_happy ~ 1 + ", nocontrols.list,
    "| any_healthcare ~ factor(initial_health_location)")),
    data = analysis.data)
print(summary(happy.iv))

# Show the structural estimate for mediator complier's D -> Y effect.
mediate.est <- mediate.selection(
    Y = "Y_health", Z = "lottery_iv", D = "any_healthcare",
    X_iv = "initial_health_location",
    X_minus = controls.list,
    boot.reps = NULL,
    type = "parametric",
    data = analysis.data)
print(coeftable(mediate.est)["AIE", "Estimate"] /
    coeftable(mediate.est)["First-stage", "Estimate"])

# Get the F statistic, for location -> D (healthcare.)
library(car)
print(mean(analysis.data$any_healthcare, na.rm = TRUE))
location.reg <- lm(
    formula(paste0("any_healthcare ~ 1 + factor(initial_health_location)",
        " + lottery_iv * ", nocontrols.list)),
    data = analysis.data)
print(summary(location.reg))
iv.list <- paste0("factor(initial_health_location)", 2:7)
print(car::linearHypothesis(location.reg, test = "F", iv.list))


################################################################################
## Estimate the CM effects with my methods.

# State how many bootstrap replications are needed.
boot.reps <- 10^3

## Panel A: Self-reported healthiness.
# Naive selection-on-observables
unadjusted.est <- mediate.selection(
    Y = "Y_health", Z = "lottery_iv", D = "any_healthcare",
    X_iv = "initial_health_location", X_minus = controls.list,
    type = "unadjusted",
    data = analysis.data,
    boot.reps = boot.reps)
print(unadjusted.est)
# Parametric MTE
parametric.est <- mediate.selection(
    Y = "Y_health", Z = "lottery_iv", D = "any_healthcare",
    X_iv = "initial_health_location", X_minus = controls.list,
    type = "parametric",
    data = analysis.data,
    boot.reps = boot.reps)
print(parametric.est)
# Semi-parametric MTE
semiparametric.est <- mediate.selection(
    Y = "Y_health", Z = "lottery_iv", D = "any_healthcare",
    X_iv = "initial_health_location", X_minus = controls.list,
    type = "semi-parametric",
    data = analysis.data,
    boot.reps = boot.reps)
print(semiparametric.est)

# Show the inferred indirect effect among mediator compliers.
print(100 * (coef(summary(unadjusted.est))["AIE", "Estimate"] / 100) /
    (coef(summary(unadjusted.est))["First-stage", "Estimate"] / 100))
print(100 * (coef(summary(parametric.est))["AIE", "Estimate"] / 100) /
    (coef(summary(parametric.est))["First-stage", "Estimate"] / 100))
print(100 * (coef(summary(semiparametric.est))["AIE", "Estimate"] / 100) /
    (coef(summary(semiparametric.est))["First-stage", "Estimate"] / 100))

# Extract the relevant estimates.
panelA.data <- data.frame(
    unadjusted_point = coef(summary(unadjusted.est))[, "Estimate"],
    unadjusted_se = coef(summary(unadjusted.est))[, "SE"],
    parametric_point = coef(summary(parametric.est))[, "Estimate"],
    parametric_se = coef(summary(parametric.est))[, "SE"],
    semiparametric_point = coef(summary(semiparametric.est))[, "Estimate"],
    semiparametric_se = coef(summary(semiparametric.est))[, "SE"])

# Save the estimates in data.
effects.extract <- function(mediate.est, model.name){
    # Compile the mediation regression results.
    reg.summary <- summary(mediate.est)
    # Get the total effect estimates.
    total.est      <- coeftable(reg.summary)["ATE", "Estimate"]
    total.ci.upper <- total.est + 1.96 * coeftable(reg.summary)["ATE", "SE"]
    total.ci.lower <- total.est - 1.96 * coeftable(reg.summary)["ATE", "SE"]
    # Get the direct effect estimates.
    direct.est       <- coeftable(reg.summary)["ADE", "Estimate"]
    direct.ci.upper  <- direct.est + 1.96 * coeftable(reg.summary)["ADE", "SE"]
    direct.ci.lower  <- direct.est - 1.96 * coeftable(reg.summary)["ADE", "SE"]
    # Get the indirect effect estimates.
    indirect.est       <- coeftable(reg.summary)["AIE", "Estimate"]
    indirect.ci.upper  <- indirect.est + 1.96 * coeftable(reg.summary)["AIE", "SE"]
    indirect.ci.lower  <- indirect.est - 1.96 * coeftable(reg.summary)["AIE", "SE"]
    # Get the percent mediated estimates.
    permediated.est      <- coeftable(reg.summary)["Proportion, AIE / ATE", "Estimate"]
    permediated.ci.upper <- permediated.est + 1.96 * coeftable(reg.summary)["Proportion, AIE / ATE", "SE"]
    permediated.ci.lower <- permediated.est - 1.96 * coeftable(reg.summary)["Proportion, AIE / ATE", "SE"]
    # Put it all into a dataframe.
    data.return <- data.frame(
        effect = c("Total", "Direct", "Indirect", "Percent Mediated"),
        pointest = c(total.est, direct.est, indirect.est, permediated.est),
        upperest = c(total.ci.upper, direct.ci.upper, indirect.ci.upper, permediated.ci.upper),
        lowerest = c(total.ci.lower, direct.ci.lower, indirect.ci.lower, permediated.ci.lower))
    # Label it with the model name
    data.return <- data.return %>% mutate(model = model.name)
    return(data.return)
}

# Collect estimates for plotting.
Y_health.data <- rbind(
    effects.extract(unadjusted.est,     "Conventional"),
    effects.extract(parametric.est,     "Parametric MTE"),
    effects.extract(semiparametric.est, "Semi-parametric MTE"))

# Clean up the data for a table.
panelA.table <- panelA.data %>%
    signif(digits.no) %>%
    format(scientific = FALSE)

panelA.table$unadjusted_se <- panelA.table$unadjusted_se %>% paste0("(", ., ")")
panelA.table$parametric_se <- panelA.table$parametric_se %>% paste0("(", ., ")")
panelA.table$semiparametric_se <- panelA.table$semiparametric_se %>% paste0("(", ., ")")
panelA.table <- data.frame(t(panelA.table))
# Add on the first column
panelA.table$model <- c(
    "Unadjusted", "", "Parametric MTE", "", "Semi-parametric MTE", "")
panelA.table <- panelA.table[c(6, 1:5)]

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


## Panel B: Self-reported happiness.
# Naive selection-on-observables
unadjusted.est <- mediate.selection(
    Y = "Y_happy", Z = "lottery_iv", D = "any_healthcare",
    X_iv = "initial_health_location", X_minus = controls.list,
    type = "unadjusted",
    data = analysis.data,
    boot.reps = boot.reps)
print(unadjusted.est)
# Parametric MTE
parametric.est <- mediate.selection(
    Y = "Y_happy", Z = "lottery_iv", D = "any_healthcare",
    X_iv = "initial_health_location", X_minus = controls.list,
    type = "parametric",
    data = analysis.data,
    boot.reps = boot.reps)
print(parametric.est)
# Semi-parametric MTE
semiparametric.est <- mediate.selection(
    Y = "Y_happy", Z = "lottery_iv", D = "any_healthcare",
    X_iv = "initial_health_location", X_minus = controls.list,
    type = "semi-parametric",
    data = analysis.data,
    boot.reps = boot.reps)
print(semiparametric.est)

# Extract the relevant figures.
panelB.data <- data.frame(
    unadjusted_point = coef(summary(unadjusted.est))[, "Estimate"],
    unadjusted_se = coef(summary(unadjusted.est))[, "SE"],
    parametric_point = coef(summary(parametric.est))[, "Estimate"],
    parametric_se = coef(summary(parametric.est))[, "SE"],
    semiparametric_point = coef(summary(semiparametric.est))[, "Estimate"],
    semiparametric_se = coef(summary(semiparametric.est))[, "SE"])

# Save for a plot.
Y_happy.data <- rbind(
    effects.extract(unadjusted.est,     "Conventional"),
    effects.extract(parametric.est,     "Parametric MTE"),
    effects.extract(semiparametric.est, "Semi-parametric MTE"))

# Show the inferred controlled indirect effect.
print(100 * (coef(summary(unadjusted.est))["AIE", "Estimate"] / 100) /
    (coef(summary(unadjusted.est))["First-stage", "Estimate"] / 100))
print(100 * (coef(summary(parametric.est))["AIE", "Estimate"] / 100) /
    (coef(summary(parametric.est))["First-stage", "Estimate"] / 100))
print(100 * (coef(summary(semiparametric.est))["AIE", "Estimate"] / 100) /
    (coef(summary(semiparametric.est))["First-stage", "Estimate"] / 100))

# Clean up the data.
panelB.table <- panelB.data %>%
    signif(digits.no) %>%
    format(scientific = FALSE)

panelB.table$unadjusted_se <- panelB.table$unadjusted_se %>% paste0("(", ., ")")
panelB.table$parametric_se <- panelB.table$parametric_se %>% paste0("(", ., ")")
panelB.table$semiparametric_se <- panelB.table$semiparametric_se %>% paste0("(", ., ")")
panelB.table <- data.frame(t(panelB.table))
# Add on the first column
panelB.table$model <- c(
    "Unadjusted", "", "Parametric MTE", "", "Semi-parametric MTE", "")
panelB.table <- panelB.table[c(6, 1:5)]

# Save the LaTeX table
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
## Plot the results

# Plot health results in a bar chart.
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
        limits = c("Conventional", "Parametric MTE", "Semi-parametric MTE")) +
    scale_y_continuous(expand = c(0, 0, 0.01, 0),
        limits = c(-0.1, 8), breaks = seq(-10, 10, by = 1),
        oob = scales::rescale_none,
        name = "") +
    ggtitle("Estimate, percent effect on subjective health") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 3, 0.25, 0), "mm"),
        legend.position = c(0.66, 0.9375),
        legend.direction = "horizontal")
# Save this file
ggsave(file.path(figures.folder, "mediation-health.png"),
    plot = health_mediation.plot,
    units = "cm",
    width = presentation.width, height = presentation.height)

# Plot happiness results in a bar chart.
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
        limits = c("Conventional", "Parametric MTE", "Semi-parametric MTE")) +
    scale_y_continuous(expand = c(0, 0, 0.01, 0),
        limits = c(-0.1, 8), breaks = seq(-10, 10, by = 1),
        oob = scales::rescale_none,
        name = "") +
    ggtitle("Estimate, percent effect on subjective well-being") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 3, 0.25, 0), "mm"),
        legend.position = c(0.66, 0.9375),
        legend.direction = "horizontal")
# Save this file
ggsave(file.path(figures.folder, "mediation-happy.png"),
    plot = happy_mediation.plot,
    units = "cm",
    width = presentation.width, height = presentation.height)


################################################################################
## Same plot, but with empty plots for the MTE parts.

# Plot health results in a bar chart.
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
        limits = c("Conventional", "Parametric MTE", "Semi-parametric MTE")) +
    scale_y_continuous(expand = c(0, 0, 0.01, 0),
        limits = c(-0.1, 8), breaks = seq(-10, 10, by = 1),
        oob = scales::rescale_none,
        name = "") +
    ggtitle("Estimate, percent effect on subjective health") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 3, 0.25, 0), "mm"),
        legend.position = c(0.66, 0.9375),
        legend.direction = "horizontal")

# Save placeholder version
ggsave(file.path(figures.folder, "mediation-health-placeholder.png"),
    plot = health_mediation.placeholder,
    units = "cm",
    width = presentation.width, height = presentation.height)

# Plot happiness results in a bar chart.
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
        limits = c("Conventional", "Parametric MTE", "Semi-parametric MTE")) +
    scale_y_continuous(expand = c(0, 0, 0.01, 0),
        limits = c(-0.1, 8), breaks = seq(-10, 10, by = 1),
        oob = scales::rescale_none,
        name = "") +
    ggtitle("Estimate, percent effect on subjective well-being") +
    theme(plot.title = element_text(size = rel(1), hjust = 0),
        plot.title.position = "plot",
        plot.margin = unit(c(0, 3, 0.25, 0), "mm"),
        legend.position = c(0.66, 0.9375),
        legend.direction = "horizontal")

# Save placeholder version
ggsave(file.path(figures.folder, "mediation-happy-placeholder.png"),
    plot = happy_mediation.placeholder,
    units = "cm",
    width = presentation.width, height = presentation.height)
