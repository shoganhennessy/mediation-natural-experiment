#!/usr/bin/R
## Senan Hogan-Hennessy, 24 October 2026
## Identifying ADE + AIE with a semi-parametric MTE interpolation
## Identifying AIE(0) from the Z = 0 conditional-outcome curve,
## Updated simulation code for the semi-parametric endpoint estimator.

# Show the date:
print(format(Sys.time(), "%H:%M %Z %A, %d %B %Y"))

## Load packages, in the RENV environment
renv::load()
# Functions for data manipulation and visualisation
library(tidyverse)
# Package for distributions to sample from.
library(MASS)
# Package for semi-parametric outcome-curve estimation by splines.
library(mgcv)


## Set up the R environment
set.seed(47)
# Define number of digits in tables and graphs
digits.no <- 3
# Define where output files go.
output.folder <- file.path("sim-output")
# Define the sample size to work with.
sample.N <- 5 * 10^3


################################################################################
## Define a function to simulate data in the triangular system.

# Define a function to simulate all observed + unobserved data
roy.data <- function(rho, sigma_0, sigma_1, sigma_C,
        error.dist = "normal", sample.size = sample.N){
    ### Inputs:
    ## rho \in [-1, +1] measuring correlation between U_0, U_1.
    ## sigma_0 >= 0 measuring standard deviation of U_0.
    ## sigma_1 >= 0 measuring standard deviation of U_1.
    ## sigma_C >= 0 measuring standard deviation of U_C.
    ## sample.size: integer, representing output sample size (i.e., N).

    # First covariate, \vec X_i^-
    X_minus <- 4 + rnorm(sample.size, mean = 0, sd = 1)
    # Second covariate, the mediator instrument.
    X_IV <- 4 * runif(sample.size, -1, 1)

    # Simulate the five-dimensional joint error vector ordered as
    # (U(0,0), U(0,1), U(1,0), U(1,1), U_C).
    # The singular covariance matrix imposes U(0,d) = U(1,d), d = 0,1.
    mu.vector <- rep(0, 5)
    sigma.matrix <- matrix(c(
        sigma_0^2,                 rho * sigma_0 * sigma_1, sigma_0^2,                 rho * sigma_0 * sigma_1, 0,
        rho * sigma_0 * sigma_1,   sigma_1^2,                 rho * sigma_0 * sigma_1, sigma_1^2,                 0,
        sigma_0^2,                 rho * sigma_0 * sigma_1, sigma_0^2,                 rho * sigma_0 * sigma_1, 0,
        rho * sigma_0 * sigma_1,   sigma_1^2,                 rho * sigma_0 * sigma_1, sigma_1^2,                 0,
        0,                         0,                         0,                         0,                         sigma_C^2),
        nrow = 5, ncol = 5, byrow = TRUE)
    U_all <- mvrnorm(n = sample.size, mu = mu.vector, Sigma = sigma.matrix,
        empirical = FALSE, tol = 1e-8)
    if (error.dist == "normal"){
        U_0_0 <- U_all[, 1]
        U_0_1 <- U_all[, 2]
        U_1_0 <- U_all[, 3]
        U_1_1 <- U_all[, 4]
        U_C <- U_all[, 5]
    }
    else if (error.dist == "uniform"){
        # Uniform errors centred on zero, preserving the target SDs.
        U_0_0 <- (2 * pnorm(U_all[, 1], 0, sigma_0) - 1) * sigma_0 * sqrt(3)
        U_0_1 <- (2 * pnorm(U_all[, 2], 0, sigma_1) - 1) * sigma_1 * sqrt(3)
        U_1_0 <- (2 * pnorm(U_all[, 3], 0, sigma_0) - 1) * sigma_0 * sqrt(3)
        U_1_1 <- (2 * pnorm(U_all[, 4], 0, sigma_1) - 1) * sigma_1 * sqrt(3)
        U_C <- (2 * pnorm(U_all[, 5], 0, sigma_C) - 1) * sigma_C * sqrt(3)
    }
    else {
        stop("Choose an error.dist from `normal`, `uniform`.")
    }
    # Define the mean potential outcomes.
    mu_outcome_z_d_X <- function(z, d, x_minus){
        return(x_minus + z + d + z * d)
    }
    mu_cost_z_X <- function(z, x_minus, x_iv){
        # Lower mediator costs under treatment increase mediator take-up.
        return(-3 * z + x_minus - x_iv)
    }
    # Potential outcomes: Y_i(z,d) = mu_d(z;X_i) + U_i(z,d).
    Y_0_0 <- mu_outcome_z_d_X(0, 0, X_minus) + U_0_0
    Y_0_1 <- mu_outcome_z_d_X(0, 1, X_minus) + U_0_1
    Y_1_0 <- mu_outcome_z_d_X(1, 0, X_minus) + U_1_0
    Y_1_1 <- mu_outcome_z_d_X(1, 1, X_minus) + U_1_1
    # Roy selection: D_i(z) = 1{Y_i(z,1) - Y_i(z,0) >= C_i(z)}.
    D_0 <- as.integer(Y_0_1 - Y_0_0 >=
        mu_cost_z_X(0, X_minus, X_IV) + U_C)
    D_1 <- as.integer(Y_1_1 - Y_1_0 >=
        mu_cost_z_X(1, X_minus, X_IV) + U_C)
    # Observed treatment, mediator, and outcome.
    probZ <- 0.5
    Z <- rbinom(sample.size, 1, probZ)
    D <- Z * D_1 + (1 - Z) * D_0
    Y <- Z * D * Y_1_1 +
        Z * (1 - D) * Y_1_0 +
        (1 - Z) * D * Y_0_1 +
        (1 - Z) * (1 - D) * Y_0_0
    combined.data <- data.frame(
        # Observed data
        Z, D, Y, X_minus, X_IV,
        # Potential mediator choices and outcomes
        D_0, D_1,
        Y_0_0, Y_0_1, Y_1_0, Y_1_1,
        # Unobserved errors
        U_0_0, U_0_1, U_1_0, U_1_1, U_C)
    return(combined.data)
}


################################################################################
## Define a function to show the theoretical values for the data.

theoretical.values <- function(sim.data, digits.no = 3, print.truth = FALSE){
    ### Inputs:
    ## sim.data, a data frame simulated from roy.data().
    D <- sim.data$D
    D_0 <- sim.data$D_0
    D_1 <- sim.data$D_1
    Y_0_0 <- sim.data$Y_0_0
    Y_0_1 <- sim.data$Y_0_1
    Y_1_0 <- sim.data$Y_1_0
    Y_1_1 <- sim.data$Y_1_1
    # First stage.
    first_stage <- D_1 - D_0
    average_first_stage <- mean(first_stage)
    # ATE = E[Y(1,D(1)) - Y(0,D(0))].
    Y_1_D1 <- D_1 * Y_1_1 + (1 - D_1) * Y_1_0
    Y_0_D0 <- D_0 * Y_0_1 + (1 - D_0) * Y_0_0
    total_effect <- Y_1_D1 - Y_0_D0
    average_total_effect <- mean(total_effect)
    # AIE(0) = E[Y(0,D(1)) - Y(0,D(0))].
    Y_0_D1 <- D_1 * Y_0_1 + (1 - D_1) * Y_0_0
    indirect_effect_0 <- Y_0_D1 - Y_0_D0
    average_indirect_effect <- mean(indirect_effect_0)
    # ADE(1) = E[Y(1,D(1)) - Y(0,D(1))].
    direct_effect_1 <- Y_1_D1 - Y_0_D1
    average_direct_effect <- mean(direct_effect_1)
    if (print.truth == TRUE){
        print("Here is a summary of the (unobserved) true effects:")
        print("Mediator principal strata, table(D_0, D_1) / N:")
        print(table(D_0, D_1) / NROW(sim.data))
        print("Pr(D = 1):")
        print(mean(D))
        print(paste0("The ATE: ",
            round(average_total_effect, digits.no)))
        print(paste0("The average first-stage: ",
            round(average_first_stage, digits.no)))
        print(paste0("AIE(0): ",
            round(average_indirect_effect, digits.no)))
        print(paste0("ADE(1): ",
            round(average_direct_effect, digits.no)))
        print(paste0("ATE - ADE(1) - AIE(0): ",
            round(average_total_effect - average_direct_effect -
                average_indirect_effect, digits.no)))
    }
    output.list <- list(
        average_first_stage     = average_first_stage,
        average_total_effect    = average_total_effect,
        average_direct_effect   = average_direct_effect,
        average_indirect_effect = average_indirect_effect)
    return(output.list)
}


################################################################################
## Define helper functions for prediction and support diagnostics.

# Create counterfactual treatment copies for prediction.
treatment.copies <- function(example.data){
    example_Z0.data <- example.data
    example_Z1.data <- example.data
    example_Z0.data$Z <- 0
    example_Z1.data$Z <- 1
    return(list(Z0 = example_Z0.data, Z1 = example_Z1.data))
}

# Summarise whether endpoint propensity scores lie in the observed Z = 0 support.
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


################################################################################
## Conventional OLS comparison estimator.

mediate.ols <- function(example.data){
    copies <- treatment.copies(example.data)
    example_Z0.data <- copies$Z0
    example_Z1.data <- copies$Z1
    # Total effect from quasi-random treatment assignment.
    totaleffect.reg <- lm(Y ~ 1 + Z + X_minus, data = example.data)
    ate.est <- mean(predict(totaleffect.reg, newdata = example_Z1.data) -
        predict(totaleffect.reg, newdata = example_Z0.data))
    # Naive mediation regression, retained only as a comparison.
    secondstage.reg <- lm(Y ~ 1 + Z * D + X_minus, data = example.data)
    example_Z0_D0.data <- example_Z0.data
    example_Z0_D1.data <- example_Z0.data
    example_Z0_D0.data$D <- 0
    example_Z0_D1.data$D <- 1
    mediator.effect.0 <- predict(secondstage.reg, newdata = example_Z0_D1.data) -
        predict(secondstage.reg, newdata = example_Z0_D0.data)
    # First-stage eestimation
    firststage.reg <- glm(D ~ 1 + Z + X_IV + X_minus,
        family = binomial(link = "probit"), data = example.data)
    pi_0.est <- predict(firststage.reg, newdata = example_Z0.data,
        type = "response")
    pi_1.est <- predict(firststage.reg, newdata = example_Z1.data,
        type = "response")
    # Conventional plug-in indirect effect, shown as a biased benchmark.
    aie.est <- mean((pi_1.est - pi_0.est) * mediator.effect.0)
    ade.est <- ate.est - aie.est
    return(list(
        "first-stage"     = mean(pi_1.est - pi_0.est),
        "total-effect"    = ate.est,
        "direct-effect"   = ade.est,
        "indirect-effect" = aie.est))
}


################################################################################
## Semi-parametric endpoint estimator of AIE(0) and ADE(1).

mediate.semiparametric <- function(example.data,
        spline.k = 10, spline.bs = "cr", support.warning = FALSE){
    ### The estimator follows the paper's new identification argument:
    ## 1. Estimate the mediator propensity score on the full sample.
    ## 2. Estimate m(p;0,X^-) in the Z = 0 subsample with mgcv.
    ## 3. Predict m at pi(1;X) and pi(0;X) for Z = 0 observations.
    ## 4. Average endpoint differences to estimate AIE(0).
    ## 5. Estimate ADE(1) = ATE - AIE(0).
    copies <- treatment.copies(example.data)
    example_Z0.data <- copies$Z0
    example_Z1.data <- copies$Z1
    # 1. Estimate the total effect from quasi-random treatment assignment.
    totaleffect.reg <- lm(Y ~ 1 + Z + X_minus, data = example.data)
    ate.est <- mean(predict(totaleffect.reg, newdata = example_Z1.data) -
        predict(totaleffect.reg, newdata = example_Z0.data))
    # 2. Estimate the mediator propensity score.
    # The probit specification is correctly specified for the normal DGP.
    cf_firststage.reg <- glm(D ~ 1 + Z + X_IV + X_minus,
        family = binomial(link = "probit"), data = example.data)
    example.data$pi.est <- predict(cf_firststage.reg, type = "response")
    pi_0.est <- predict(cf_firststage.reg, newdata = example_Z0.data,
        type = "response")
    pi_1.est <- predict(cf_firststage.reg, newdata = example_Z1.data,
        type = "response")
    # 3. Restrict outcome-curve estimation and standardisation to Z = 0.
    outcome_Z0.data <- example.data[example.data$Z == 0, , drop = FALSE]
    endpoint_Z0.data <- outcome_Z0.data
    endpoint_Z0.data$pi_0.est <- pi_0.est[example.data$Z == 0]
    endpoint_Z0.data$pi_1.est <- pi_1.est[example.data$Z == 0]
    # Check that the endpoint propensities are covered by the realised Z = 0
    # propensity-score support. The estimator does not trim or clamp endpoints,
    # because doing so would conceal unsupported extrapolation.
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
    # 4. Estimate m(p;0,X^-) directly in the Z = 0 subsample.
    # Use a cubic regression spline with REML smoothing selection.
    cf_outcome_Z0.reg <- gam(
        Y ~ 1 + X_minus + s(pi.est, bs = spline.bs, k = spline.k),
        method = "REML", data = outcome_Z0.data)
    # 5. Evaluate the same Z = 0 outcome curve at both endpoint propensities.
    predict_pi0.data <- endpoint_Z0.data
    predict_pi1.data <- endpoint_Z0.data
    predict_pi0.data$pi.est <- endpoint_Z0.data$pi_0.est
    predict_pi1.data$pi.est <- endpoint_Z0.data$pi_1.est
    m_pi0.est <- predict(cf_outcome_Z0.reg, newdata = predict_pi0.data,
        type = "response")
    m_pi1.est <- predict(cf_outcome_Z0.reg, newdata = predict_pi1.data,
        type = "response")
    conditional.aie0.est <- m_pi1.est - m_pi0.est
    aie0.est <- mean(conditional.aie0.est, na.rm = TRUE)
    ade1.est <- ate.est - aie0.est
    # Return all the objects
    output.list <- list(
        "first-stage"     = mean(
            endpoint_Z0.data$pi_1.est - endpoint_Z0.data$pi_0.est, na.rm = TRUE),
        "total-effect"    = ate.est,
        "direct-effect"   = ade1.est,
        "indirect-effect" = aie0.est,
        "support"         = support.est,
        "firststage.reg"  = cf_firststage.reg,
        "outcome.reg"     = cf_outcome_Z0.reg,
        "conditional-aie" = conditional.aie0.est)
    return(output.list)
}


################################################################################
## Define a function to automate DGP re-simulation or bootstrap estimation.

estimated.loop <- function(boot.reps, example.data,
        bootstrap = TRUE, print.progress = FALSE, error.dist = "normal",
        # DGP parameters
        rho, sigma_0, sigma_1, sigma_C,
        spline.k = 10, spline.bs = "cr") {
    # 1. Naive OLS benchmark.
    ols_total_effect <- c()
    ols_direct_effect <- c()
    ols_indirect_effect <- c()
    # 2. Semi-parametric endpoint estimator.
    cf_total_effect <- c()
    cf_direct_effect <- c()
    cf_indirect_effect <- c()
    cf_support_share <- c()
    # Truth values.
    truth_total_effect <- c()
    truth_direct_effect <- c()
    truth_indirect_effect <- c()
    # Get the truth estimate.
    truth.est <- theoretical.values(example.data)
    # Get the point estimates.
    for (i in seq_len(boot.reps)){
        if (bootstrap == TRUE) {
            boot.indicies <- sample(seq_len(NROW(example.data)),
                NROW(example.data), replace = TRUE)
            boot.data <- example.data[boot.indicies, ]
        }
        else if (bootstrap == FALSE){
            boot.data <- roy.data(rho, sigma_0, sigma_1, sigma_C,
                error.dist = error.dist)
            truth.est <- theoretical.values(boot.data)
        }
        else {
            stop("The `bootstrap` option only takes values of TRUE or FALSE.")
        }
        if (print.progress == TRUE &&
                (i == 1 || i == boot.reps || i %% max(1, floor(boot.reps / 100)) == 0)){
            print(paste0(i, " out of ", boot.reps, ", ",
                round(100 * i / boot.reps), "% done."))
        }
        ols.est <- tryCatch(mediate.ols(boot.data),
            error = function(e) NULL)
        cf.est <- tryCatch(mediate.semiparametric(boot.data,
            spline.k = spline.k, spline.bs = spline.bs),
            error = function(e) NULL)
        truth_total_effect[i] <- truth.est$average_total_effect
        truth_direct_effect[i] <- truth.est$average_direct_effect
        truth_indirect_effect[i] <- truth.est$average_indirect_effect
        if (!is.null(ols.est)){
            ols_total_effect[i] <- ols.est$"total-effect"
            ols_direct_effect[i] <- ols.est$"direct-effect"
            ols_indirect_effect[i] <- ols.est$"indirect-effect"
        }
        else {
            ols_total_effect[i] <- NA
            ols_direct_effect[i] <- NA
            ols_indirect_effect[i] <- NA
        }
        if (!is.null(cf.est)){
            cf_total_effect[i] <- cf.est$"total-effect"
            cf_direct_effect[i] <- cf.est$"direct-effect"
            cf_indirect_effect[i] <- cf.est$"indirect-effect"
            cf_support_share[i] <- cf.est$"support"$share.interval.in.support
        }
        else {
            cf_total_effect[i] <- NA
            cf_direct_effect[i] <- NA
            cf_indirect_effect[i] <- NA
            cf_support_share[i] <- NA
        }
    }
    # Return the output.
    output.list <- list()
    output.list$data <- data.frame(
        truth_total_effect = truth_total_effect,
        truth_direct_effect = truth_direct_effect,
        truth_indirect_effect = truth_indirect_effect,
        ols_total_effect = ols_total_effect,
        ols_direct_effect = ols_direct_effect,
        ols_indirect_effect = ols_indirect_effect,
        cf_total_effect = cf_total_effect,
        cf_direct_effect = cf_direct_effect,
        cf_indirect_effect = cf_indirect_effect,
        cf_support_share = cf_support_share)
    summarise.est <- function(x){
        c(mean = mean(x, na.rm = TRUE),
            se = sd(x, na.rm = TRUE),
            low = as.numeric(quantile(x, 0.025, na.rm = TRUE)),
            up = as.numeric(quantile(x, 0.975, na.rm = TRUE)))
    }
    ols.direct.summary <- summarise.est(ols_direct_effect)
    ols.indirect.summary <- summarise.est(ols_indirect_effect)
    cf.direct.summary <- summarise.est(cf_direct_effect)
    cf.indirect.summary <- summarise.est(cf_indirect_effect)
    output.list$estimates <- data.frame(
        truth_direct_effect = mean(truth_direct_effect, na.rm = TRUE),
        truth_indirect_effect = mean(truth_indirect_effect, na.rm = TRUE),
        ols_direct_effect = ols.direct.summary["mean"],
        ols_direct_effect_se = ols.direct.summary["se"],
        ols_direct_effect_low = ols.direct.summary["low"],
        ols_direct_effect_up = ols.direct.summary["up"],
        ols_indirect_effect = ols.indirect.summary["mean"],
        ols_indirect_effect_se = ols.indirect.summary["se"],
        ols_indirect_effect_low = ols.indirect.summary["low"],
        ols_indirect_effect_up = ols.indirect.summary["up"],
        cf_direct_effect = cf.direct.summary["mean"],
        cf_direct_effect_se = cf.direct.summary["se"],
        cf_direct_effect_low = cf.direct.summary["low"],
        cf_direct_effect_up = cf.direct.summary["up"],
        cf_indirect_effect = cf.indirect.summary["mean"],
        cf_indirect_effect_se = cf.indirect.summary["se"],
        cf_indirect_effect_low = cf.indirect.summary["low"],
        cf_indirect_effect_up = cf.indirect.summary["up"],
        cf_support_share = mean(cf_support_share, na.rm = TRUE))
    return(output.list)
}


################################################################################
## Compare estimation methods in a one-shot simulation.

rho <- 0.5
sigma_0 <- 1
sigma_1 <- 1.5
sigma_C <- 0.5

# Get the simulated data.
simulated.data <- roy.data(rho, sigma_0, sigma_1, sigma_C)
print(theoretical.values(simulated.data, print.truth = TRUE))

# Show how the naive mediation estimator is biased.
print(mediate.ols(simulated.data))

# Show the new endpoint estimator of AIE(0), then ADE(1) = ATE - AIE(0).
cf.one.shot <- mediate.semiparametric(simulated.data,
    spline.k = 10, spline.bs = "cr", support.warning = TRUE)
print(summary(cf.one.shot$outcome.reg))


################################################################################
## Repeated-DGP simulation for normally distributed errors.

# Large simulated population for approximate truth values.
roy.data(rho, sigma_0, sigma_1, sigma_C, sample.size = 10^6,
    error.dist = "normal") %>%
    theoretical.values(print.truth = TRUE) %>%
    print()

normal.data <- roy.data(rho, sigma_0, sigma_1, sigma_C,
    error.dist = "normal")

# Change to 10^4 for the final simulation. A smaller value is useful for testing.
sim.reps <- 10^3
normal.est <- estimated.loop(sim.reps, normal.data,
    bootstrap = FALSE, print.progress = TRUE, error.dist = "normal",
    rho, sigma_0, sigma_1, sigma_C,
    spline.k = 10, spline.bs = "cr")

normal.est$data %>%
    write_csv(file.path(output.folder, "normal-endpoint-cf-data.csv"))
normal.est$estimates %>%
    write_csv(file.path(output.folder, "normal-endpoint-cf-estimates.csv"))


################################################################################
## Optional: bootstrap inference for one simulated sample.

# Set boot.reps to the desired number for final inference.
boot.reps <- 10^3
normal.boot <- estimated.loop(boot.reps, normal.data,
    bootstrap = TRUE, print.progress = TRUE, error.dist = "normal",
    rho, sigma_0, sigma_1, sigma_C,
    spline.k = 10, spline.bs = "cr")

normal.boot$data %>%
    write_csv(file.path(output.folder, "normal-endpoint-cf-bootstrap.csv"))
normal.boot$estimates %>%
    write_csv(file.path(output.folder, "normal-endpoint-cf-bootstrap-estimates.csv"))


################################################################################
## Plot repeated-DGP estimation errors.

plot.data <- normal.est$data %>%
    transmute(
        OLS_AIE0 = ols_indirect_effect - truth_indirect_effect,
        MTE_AIE0 = cf_indirect_effect - truth_indirect_effect,
        OLS_ADE1 = ols_direct_effect - truth_direct_effect,
        MTE_ADE1 = cf_direct_effect - truth_direct_effect) %>%
    pivot_longer(cols = everything(), names_to = "estimand_method",
        values_to = "error") %>%
    separate(estimand_method, into = c("method", "estimand"), sep = "_") %>%
    mutate(
        method = recode(method,
            "OLS" = "Conventional CM",
            "MTE" = "Semi-parametric MTE"),
        estimand = recode(estimand,
            "AIE0" = "AIE(0)",
            "ADE1" = "ADE(1)"))

error.plot <- ggplot(plot.data, aes(x = error, fill = method, colour = method)) +
    geom_density(alpha = 0.25, na.rm = TRUE) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
    facet_wrap(~ estimand, scales = "free") +
    labs(x = "Estimate minus true effect", y = "Density",
        fill = NULL, colour = NULL) +
    theme_minimal()

print(error.plot)
#ggsave(file.path(output.folder, "normal-endpoint-cf-density.png"),
#    error.plot, height = 6, width = 10, dpi = 300)
