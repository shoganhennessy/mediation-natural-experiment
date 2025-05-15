#!/usr/bin/R
## Senan Hogan-Hennessy, 18 April 2025
## Identifying ADE + AIE with a parametric control function
## see Hogan-Hennessy (2025).

# Show the date:
print(format(Sys.time(), "%H:%M %Z %A, %d %B %Y"))

## Load libraries
# Functions for data manipulation and visualisation
library(tidyverse)
# Causal medation package, Imai Keele Yamamoto (2010)
library(mediation)
# Package for classical selection estimators (i.e., MLE)
library(sampleSelection)
# Package for more distributions to sample from.
library(MASS)
# Package for semi-parametric regressor, splines by bs(.).
library(splines)

## Set up the R environment
#set.seed(47)
# Define number of digits in tables and graphs
digits.no <- 3
# Define where output files go.
output.folder <- file.path("sim-output")
# Set the options for the plot sizes, in saving ggplot output.
fig.height <- 10
fig.width <- fig.height

# Define the sample size to work with.
sample.N <- 10^3


################################################################################
## Define a function to simulate data in the triangular system.

# Define a function to simulate all observed + unobserved data 
simulate.data <- function(rho, sigma_0, sigma_1, sigma_C,
        sample.size = sample.N){
    ### Inputs:
    ## X, a matrix of covariates, continuous or binary values.
    ## rho \in [-1, +1] measuring correlation between U_0, U_1.
    ## sigma_0 >= 0 measuring standard deviation of U_0.
    ## sigma_1 >= 0 measuring standard deviation of U_1.
    ## sigma_C >= 0 measuring standard deviation of U_C.
    ## sample.size: integer, representing output sample size (i.e., N).
    # First covariate (\vec X_i^-)
    X_minus <- 4 + rnorm(sample.size, mean = 0, sd = 1)
    # Second covariate (instrument for the control function).
    X_IV <- 2 * runif(sample.size, -1, 1)
    # Simulate the unobserved error terms.
    U_all <- mvrnorm(n = sample.size, mu = c(0, 0, 0),
        Sigma = matrix(c(
            sigma_0^2,               rho * sigma_0 * sigma_1, 0,
            rho * sigma_0 * sigma_1, sigma_1^2,               0,
            0,                       0,                       sigma_C^2),
                ncol = 3),
        empirical = FALSE)
    U_0 <- U_all[, 1]
    U_1 <- U_all[, 2]
    U_C <- U_all[, 3]
    # Define the mean potential outcomes.
    mu_outcome_z_d_X <- function(z, d, x_minus){
        return(x_minus + (z + d + z * d))
    }
    mu_cost_z_X <- function(z, x_minus, x_iv){
        return(- 3 * z + x_minus - x_iv)
    }
    # Mean outcomes: Y_i(z, d) = mu_d(z; X_i) + U_D
    Y_0_0 <- mu_outcome_z_d_X(0, 0, X_minus) + U_0
    Y_0_1 <- mu_outcome_z_d_X(0, 1, X_minus) + U_1
    Y_1_0 <- mu_outcome_z_d_X(1, 0, X_minus) + U_0
    Y_1_1 <- mu_outcome_z_d_X(1, 1, X_minus) + U_1
    # Roy selection: D_i(z) = 1{ Y(z, 1) - Y(z, 0) >= C_i }
    D_0 <- as.integer(Y_0_1 - Y_0_0 >= mu_cost_z_X(0, X_minus, X_IV) + U_C)
    D_1 <- as.integer(Y_1_1 - Y_1_0 >= mu_cost_z_X(1, X_minus, X_IV) + U_C)
    # Observed outcomes: treatment Z
    probZ <- 0.5
    Z <- rbinom(sample.size, 1, probZ)
    # Observed outcomes: mediator D
    D <- (Z * D_1) + ((1 - Z) * D_0)
    # Observed outcomes: outcome Y
    Y <- (Z * D * Y_1_1) +
        (Z * (1 - D) * Y_1_0) +
        ((1 - Z) * D * Y_0_1) +
        ((1 - Z) * (1 - D) * Y_0_0)
    # Put these data to a coherent data frame.
    combined.data <- data.frame(
        # Observed data
        Z, D, Y,  X_minus, X_IV,
        # Unobserved, potential outcomes and compliance.
        D_0, D_1,
        Y_0_0, Y_0_1, Y_1_0, Y_1_1,
        #mu0_Z0_X, mu1_Z0_X, mu0_Z1_X, mu1_Z1_X, muC_Z0_X, muC_Z1_X, 
        U_0, U_1, U_C)
    # Return the simulated data as a data frame.
    return(combined.data)
}


################################################################################
## Define a function to show the theoretical values for the data.
theoretical.values <- function(sim.data, digits.no = 3, print.truth = FALSE){
    ### Inputs:
    ## sim.data, a data frame simulated from above.
    # Extract the potentials from simulated data.
    Z <- sim.data$Z
    D <- sim.data$D
    Y <- sim.data$Y
    X_minus <- sim.data$X_minus
    X_IV <- sim.data$X_IV
    D_0 <- sim.data$D_0
    D_1 <- sim.data$D_1
    Y_0_0 <- sim.data$Y_0_0
    Y_0_1 <- sim.data$Y_0_1
    Y_1_0 <- sim.data$Y_1_0
    Y_1_1 <- sim.data$Y_1_1
    U_0 <- sim.data$U_0
    U_1 <- sim.data$U_1
    U_C <- sim.data$U_C
    # Get the true first-stage effects
    first_stage <- D_1 - D_0
    average_first_stage <- mean(first_stage)
    # Get the theoretical total effect/reduced form/ATE
    total_effect <-
        (Y_1_1 - Y_0_0) * (D_1 == 1 & D_0 == 0) +
        (Y_1_1 - Y_0_1) * (D_1 == 1 & D_0 == 1) +
        (Y_1_0 - Y_0_0) * (D_1 == 0 & D_0 == 0)
    average_total_effect <- mean(total_effect)
    # Get the theoretical indirect effect.
    indirect_effect <-
        (Z * (Y_1_1 - Y_1_0) + (1 - Z) * (Y_0_1 - Y_0_0))
    controlled_indirect_effect <- mean(indirect_effect)
    average_indirect_effect <- mean((D_1 == 1 & D_0 == 0) * indirect_effect)
    # Get the theoretical direct effect.
    direct_effect <- (D * (Y_1_1 - Y_0_1) + (1 - D) * (Y_1_0 - Y_0_0))
    average_direct_effect <- mean(direct_effect)
    # Show the real underlying values.
    if (print.truth == TRUE){
        print("Here is a summary of the (unobserved) true effects:")
        # Show how many ATs, NTs, Compliers in terms of D_i(Z) for Z = 0, 1.
        print("How many mediator compliers in the system?")
        print(table(D_1, D_0) / NROW(sim.data))
        print("How many actually took the mediator, i.e. Pr(D = 1)?")
        print(mean(D))
        # Show the real treatment effects
        print(paste0(c("The average total effect:",    as.numeric(average_total_effect))))
        print(paste0(c("The average first-stage:",     as.numeric(average_first_stage))))
        print(paste0(c("The average direct effect:",   as.numeric(average_direct_effect))))
        print(paste0(c("The average indirect effect:", as.numeric(average_indirect_effect))))
    }
    # Define a named list to return
    output.list <- list(
        average_first_stage     = average_first_stage,
        average_total_effect    = average_total_effect,
        average_direct_effect   = average_direct_effect,
        average_indirect_effect = average_indirect_effect)
    # Return the output.list
    return(output.list)
}


################################################################################
## Define a function to estimate mediation, given the first + second-stages.

# Estimate the values, given a first and second-stages
estimated.values <- function(firststage.reg, secondstage.reg, example.data,
    complier.adjustment = NULL){
    ### Inputs:
    ## example.data, a data frame simulated from above.
    # calculate the first-stage by prediction
    firststage.est <- predict(
        firststage.reg, newdata = mutate(example.data, Z = 1), type = "response") - predict(
            firststage.reg, newdata = mutate(example.data, Z = 0), type = "response")
    # calculate the second-stage direct effect
    direct.est <- predict(
        secondstage.reg, newdata = mutate(example.data, Z = 1)) -
        predict(secondstage.reg, newdata = mutate(example.data, Z = 0))
    # calculate the second-stage (controlled) indirect effect
    indirect.est <- predict(
        secondstage.reg, newdata = mutate(example.data, D = 1)) -
        predict(secondstage.reg, newdata = mutate(example.data, D = 0))
    # Add the Kline Walters (2019) IV-type complier adjustment (provided external)
    if (!is.null(complier.adjustment)) {
        indirect.est <- indirect.est + complier.adjustment
    }
    # Return the mean estimates.
    output.list <- list(
        "first-stage"     = mean(firststage.est, na.rm = TRUE),
        "direct-effect"   = mean(direct.est, na.rm = TRUE),
        "indirect-effect" = mean(firststage.est * indirect.est, na.rm = TRUE))
    # Return the output.list
    return(output.list)
}

# Define a function to Heckman selection correct mediation est, in two-stages.
mediate.heckit <- function(example.data){
    # 1. Probit first-stage (well identified).
    cf_firststage.reg <- glm(D ~ (1 + Z) + X_IV + X_minus,
        family = binomial(link = "probit"),
        data = example.data)
    # 2. Define the control functions --- with assumed N(0,1) dist.
    lambda_1.fun <- function(p){
        # Inv Mills ratio, taking as input the estimated mediator propensity.
        return(dnorm(qnorm(p)) / pnorm(qnorm(p)))
    }
    P <- predict(cf_firststage.reg, type = "response")
    example.data$lambda_0 <- (1 - example.data$D) * lambda_1.fun(P) * (- P / (1 - P))
    example.data$lambda_1 <- example.data$D * lambda_1.fun(P)
    # 3. Estimate second-stage, including the CFs.
    cf_secondstage.reg <- lm(Y ~ (1 + Z * D) + X_minus + lambda_0 + lambda_1,
        data = example.data)
    # Compensate complier difference in AIE, by Kline Walters (2019) IV-type adjustment.
    P_0 <- predict(cf_firststage.reg, newdata = mutate(example.data, Z = 0), type = "response")
    P_1 <- predict(cf_firststage.reg, newdata = mutate(example.data, Z = 1), type = "response")
    Gamma.big <-  (P_1 * lambda_1.fun(P_1) - P_0 * lambda_1.fun(P_0)) / (P_1 - P_0)
    rho_0 <- coef(cf_secondstage.reg)["lambda_0"]
    rho_1 <- coef(cf_secondstage.reg)["lambda_1"]
    add.term <- (rho_1 - rho_0) * Gamma.big
    # Return the first and second-stages, and the complier compensating term.
    output.list <- list(
        "first-stage"         = cf_firststage.reg,
        "second-stage"        = cf_secondstage.reg,
        "complier-adjustment" = add.term,
        "heckit-data"         = example.data)
    return(output.list)
}

# Bootstrap the estimates.
estimated.loop <- function(boot.reps, example.data,
        bootstrap = TRUE, print.progress = FALSE,
        # Default data parameters
        rho = 0.5, sigma_0 = 1, sigma_1 = 2, sigma_C = 2) {
    # Define lists the will be returned:
    # 2. Naive OLS.
    ols_direct_effect <- c()
    ols_indirect_effect <- c()
    # 3. Control function.
    cf_direct_effect <- c()
    cf_indirect_effect <- c()
    # More in the future ....
    # Calculate the truth values, given the input data
    truth_direct_effect <- c()
    truth_indirect_effect <- c()
    truth.est <- theoretical.values(example.data)
    ## Loop across the bootstraps values.
    for (i in seq(1, boot.reps)){
        # If bootstrapping, just resample from provided data.
        if (bootstrap == TRUE) {
            boot.indicies <- sample(
                seq(1, NROW(example.data)), NROW(example.data), replace = TRUE)
            boot.data <- example.data[boot.indicies, ]
        }
        # If a regular re-simulation, get new data.
        else if (bootstrap == FALSE){
            boot.data <- simulate.data(0.5, 1, 3, 2)
            # Update the truth values to the newest simulated data, if so.
            truth.est <- theoretical.values(boot.data)
        }
        else {stop("The `bootstrap' option only takes values of TRUE or FALSE.")}
        # Print, if want the consol output of how far we are.
        if (print.progress == TRUE){
            if ((100 * (i / boot.reps)) %% 1 == 0) {
                print(paste0(i, " out of ", boot.reps, ", ", 100 * (i / boot.reps), "% done."))
            }
        }
        # Now get the mediation effects, by different approaches.
        # 2. OLS estimate of second-stage
        ols_firststage.reg <- lm(D ~ (1 + Z) + X_minus + X_IV, data = boot.data)
        ols_secondstage.reg <- lm(Y ~ 1 + Z * D + X_minus, data = boot.data)
        ols.est <- estimated.values(ols_firststage.reg, ols_secondstage.reg,
            boot.data)
        # 3. Heckman-style selection-into-mediator model estimates.
        heckit.reg <- mediate.heckit(boot.data)
        cf.est <- estimated.values(
            heckit.reg$"first-stage",
            heckit.reg$"second-stage",
            heckit.reg$"heckit-data",
            complier.adjustment = heckit.reg$"complier-adjustment")
        # Save the outputs.
        truth_direct_effect[i]   <- truth.est$average_direct_effect
        truth_indirect_effect[i] <- truth.est$average_indirect_effect
        ols_direct_effect[i]     <- ols.est$`direct-effect`
        ols_indirect_effect[i]   <- ols.est$`indirect-effect`
        cf_direct_effect[i]      <- cf.est$`direct-effect`
        cf_indirect_effect[i]    <- cf.est$`indirect-effect`
    }
    # Return the bootstrap data.
    output.list <- list()
    output.list$data <- data.frame(
        truth_direct_effect   = truth_direct_effect,
        ols_direct_effect     = ols_direct_effect,
        cf_direct_effect      = cf_direct_effect,
        truth_indirect_effect = truth_indirect_effect,
        ols_indirect_effect   = ols_indirect_effect,
        cf_indirect_effect    = cf_indirect_effect)
    # Calculate the needed statistics, to return
    output.list$estimates <- data.frame(
        # Truth
        truth_direct_effect     = as.numeric(mean(truth_direct_effect)),
        truth_indirect_effect   = as.numeric(mean(truth_indirect_effect)),
        # OLS mean, and the 95% confidence intervals
        ols_direct_effect       = as.numeric(mean(ols_direct_effect)),
        ols_direct_effect_se    = as.numeric(sd(ols_direct_effect)),
        ols_direct_effect_up    = as.numeric(quantile(ols_direct_effect,
            probs = 0.975, na.rm = TRUE)),
        ols_direct_effect_low   = as.numeric(quantile(ols_direct_effect,
            probs = 0.025, na.rm = TRUE)),
        ols_indirect_effect     = as.numeric(mean(ols_indirect_effect)),
        ols_indirect_effect_se  = as.numeric(sd(ols_indirect_effect)),
        ols_indirect_effect_up  = as.numeric(quantile(ols_indirect_effect,
            probs = 0.975, na.rm = TRUE)),
        ols_indirect_effect_low = as.numeric(quantile(ols_indirect_effect,
            probs = 0.025, na.rm = TRUE)),
        # Control Fun mean, and the 95% confidence intervals
        cf_direct_effect        = as.numeric(mean(cf_direct_effect)),
        cf_direct_effect_se     = as.numeric(sd(cf_direct_effect)),
        cf_direct_effect_up     = as.numeric(quantile(cf_direct_effect,
            probs = 0.975, na.rm = TRUE)),
        cf_direct_effect_low    = as.numeric(quantile(cf_direct_effect,
            probs = 0.025, na.rm = TRUE)),
        cf_indirect_effect      = as.numeric(mean(cf_indirect_effect)),
        cf_indirect_effect_se   = as.numeric(sd(cf_indirect_effect)),
        cf_indirect_effect_up   = as.numeric(quantile(cf_indirect_effect,
            probs = 0.975, na.rm = TRUE)),
        cf_indirect_effect_low  = as.numeric(quantile(cf_indirect_effect,
            probs = 0.025, na.rm = TRUE)))
    return(output.list)
}


################################################################################
## Compare estimation methods, in one simulation.

## Simulate the data with given rho, sigma_0, sigma_1, sigma_C values.
rho <- 0.5
sigma_0 <- 1
sigma_1 <- 3
sigma_C <- 2
simulated.data <- simulate.data(rho, sigma_0, sigma_1, sigma_C)
# SHow the theoretical direct + indirect values
print(theoretical.values(simulated.data, print.truth = TRUE))

# Show that the regression specification holds exactly, if specified correctly.
true_firststage.reg <- glm(D ~ (1 + Z) + X_IV + X_minus +
        U_C + U_0 + U_1,
    family = binomial(link = "probit"),
    data = simulated.data)
true_secondstage.reg <- lm(Y ~ (-1 + Z * D) + X_minus +
    # including the unobserved errors:
    I((1 - D) * U_0) + I(D * U_1),
    data = simulated.data)
print(theoretical.values(simulated.data))
print(estimated.values(true_firststage.reg, true_secondstage.reg, simulated.data))
# See how the first and second-stages are perfect:
print(summary(true_firststage.reg))
print(summary(true_secondstage.reg))
#! Note: this automatically includes the complier term in the indirect effect
#!       without full access to U_0, U_1 this is not automatic.

# Show how the OLS result gives a bias result (if rho != 0),
# using the same second-stage for both direct + indirect estimates.
ols_firststage.reg <- lm(D ~ (1 + Z) + X_minus + X_IV, data = simulated.data)
ols_secondstage.reg <- lm(Y ~ 1 + Z * D + X_minus, data = simulated.data)
print(theoretical.values(simulated.data))
print(estimated.values(ols_firststage.reg, ols_secondstage.reg, simulated.data))


################################################################################
## Testing out a Heckman selection model.

# Show how a control function gets it correct, in 2 steps.
cf_firststage.reg <- glm(D ~ (1 + Z) + X_IV + X_minus,
    family = binomial(link = "probit"),
    data = simulated.data)
# Second-stage, with Heckman normal errors
lambda_1.fun <- function(p){
    return(dnorm(qnorm(p)) / pnorm(qnorm(p)))
}
P <- predict(cf_firststage.reg, type = "response")
simulated.data$lambda_0 <- (1 - simulated.data$D) * lambda_1.fun(P) * (- P / (1 - P))
simulated.data$lambda_1 <- simulated.data$D * lambda_1.fun(P)
# Second-stage, including the control functions lambda_0, lambda_1
cf_secondstage.reg <- lm(Y ~ (1 + Z * D) + X_minus + lambda_0 + lambda_1,
    data = simulated.data)
print(summary(cf_secondstage.reg))
# Show how it gets these correct (with complier adjustment).
print(theoretical.values(simulated.data))
print(estimated.values(cf_firststage.reg, cf_secondstage.reg, simulated.data))
print(estimated.values(cf_firststage.reg, cf_secondstage.reg, simulated.data,
    complier.adjustment = mediate.heckit(simulated.data)$"complier-adjustment"))


################################################################################
## Testing out the Semi-parametric approach.

# 1. Non-parametric first-stage.
cf_firststage.reg <- gam(D ~ 1 + Z + s(X_IV) + s(X_minus),
    family = binomial, data = simulated.data)
print(summary(cf_firststage.reg))
P <- predict(cf_firststage.reg, type = "response")
P_0 <- predict(cf_firststage.reg, newdata = mutate(simulated.data, Z = 0), type = "response")
P_1 <- predict(cf_firststage.reg, newdata = mutate(simulated.data, Z = 1), type = "response")
# Second-stage, with semi-parametric CF.
simulated.data$intercept <- 1
simulated.data$P <- P

# (1) Estimate \lambda_1 in D = 1 sample.
cf_D1_semi.reg <- lm(Y ~ (0 + intercept + Z) + X_minus +
    bs(P, knots = seq(0, 1, by = 0.1), intercept = TRUE),
    data = filter(simulated.data, D == 1))
print(summary(cf_D1_semi.reg))
simulated.data$hat_lambda_1 <- predict(cf_D1_semi.reg,
    newdata = mutate(simulated.data, intercept = 0, Z = 0, X_minus = 0))
# (2) Use the CFs in the estimation
# Compose the cross-estimates of \tilde \lambda_1(pi)
simulated.data$lambda_0 <- (1 - simulated.data$D) * (-P / (1 - P)) * simulated.data$hat_lambda_1
simulated.data$lambda_1 <- simulated.data$D * simulated.data$hat_lambda_1
# Second-stage, including the control functions lambda_0, lambda_1
cf_secondstage.reg <- lm(I(Y - lambda_1) ~ (1 + Z * D) + X_minus + lambda_0,
    data = simulated.data)
print(summary(cf_secondstage.reg))
# Take the \tilde\rho = rho_0 / rho_1 estimate.
tilde_rho <- coef(cf_secondstage.reg)["lambda_0"]
# Lastly, the complier adjustment.
hat_lambda_1_P_0 <- predict(cf_D1_semi.reg,
    newdata = mutate(simulated.data, intercept = 0, Z = 0, X_minus = 0, P = P_0))
hat_lambda_1_P_1 <- predict(cf_D1_semi.reg,
    newdata = mutate(simulated.data, intercept = 0, Z = 0, X_minus = 0, P = P_1))
complier.adjustment <- (1 - tilde_rho) * (
    P_1 * hat_lambda_1_P_1 - P_0 * hat_lambda_1_P_0) / (P_1 - P_0)

print(estimated.values(cf_firststage.reg, cf_secondstage.reg, simulated.data,
    complier.adjustment = complier.adjustment))


################################################################################
## Plot bootstrap results for one DGP, repeatedly drawn

# Base data to input
simulated.data <- simulate.data(rho, sigma_0, sigma_1, sigma_C)

# Get bootstrapped point est for the CF approach
sim.reps <- 10^4
sim.est <- estimated.loop(sim.reps, simulated.data,
    bootstrap = FALSE, print.progress = TRUE)
sim.data <- sim.est$data

## Save the repeated DGPs' point estimates as separate data.
sim.data %>% write_csv(file.path(output.folder, "dist-heckit-data.csv"))
