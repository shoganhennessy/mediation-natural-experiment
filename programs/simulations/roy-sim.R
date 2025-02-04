#!/usr/bin/R
## Senan Hogan-Hennessy, 16 Jan 2025
## Simulate the system for indirect + direct effects, with Roy selection.
# Show the date:
print(format(Sys.time(), "%H:%M %Z %A, %d %B %Y"))

## Load libraries
# Functions for data manipulation and visualisation
library(tidyverse)
# Library for better colour choice.
library(ggthemes)
# Library for equations in plots
library(latex2exp)
# Causal medation package, Imai Keele Yamamoto (2010)
library(mediation)
# Package for classical selection estimators (i.e., MLE)
library(sampleSelection)
# Package for semi-parametric regressor
library(splines)

## Set up the R environment
set.seed(47)
# Define number of digits in tables and graphs
digits.no <- 3
# Define where output files go.
output.folder <- file.path("sim-output")
# Set the options for the plot sizes, in saving ggplot output.
fig.height <- 10
fig.width <- fig.height

# Define the sample size to work with.
sample.N <- 10^4


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
    ## sample.size: integer, representing output sample size i.e., N.
    # First covariate (\vec X_i^-)
    X_minus <- 4 + rnorm(sample.size, mean = 0, sd = 1)
    # Second covariate (instrument for the control function).
    X_IV <- rbinom(sample.size, 1, 1 / 2)
    # Simulate the unobserved error terms.
    U_both <- MASS::mvrnorm(
        n = sample.size,
        mu = c(0, 0, 0),
        Sigma = matrix(c(
            sigma_0^2,                rho * sigma_0 * sigma_1,  0,
            rho * sigma_0 * sigma_1,  sigma_1^2,                0,
            0,                        0,                        sigma_C^2), ncol = 3),
        empirical = FALSE)
    U_0 <- U_both[, 1]
    U_1 <- U_both[, 2]
    U_C <- U_both[, 3]
    # Define the mean potential outcomes.
    mu_outcome_z_d_X <- function(z, d, x_minus){
        return(x_minus + (z + d + z * d))
    }
    mu_cost_z_X <- function(z, x_minus, x_iv){
        return(- 3 * z + x_minus - x_iv)
    }
    # Y_i(Z, D) = mu_D(Z; X_i) + U_D
    Y_0_0 <- mu_outcome_z_d_X(0, 0, X_minus) + U_0
    Y_0_1 <- mu_outcome_z_d_X(0, 1, X_minus) + U_1
    Y_1_0 <- mu_outcome_z_d_X(1, 0, X_minus) + U_0
    Y_1_1 <- mu_outcome_z_d_X(1, 1, X_minus) + U_1
    # D_i(Z)= 1{ Y(Z, 1) - Y(Z, 0) >= C_i }
    D_0 <- as.integer(Y_0_1 - Y_0_0 >= mu_cost_z_X(0, X_minus, X_IV) + U_C)
    D_1 <- as.integer(Y_1_1 - Y_1_0 >= mu_cost_z_X(1, X_minus, X_IV) + U_C)
    # Generate the individual effects (direct + indirect)
    probZ <- 0.5
    Z <- rbinom(sample.size, 1, probZ)
    # Observed outcomes: D, Y
    D <- (Z * D_1) + ((1 - Z) * D_0)
    # Generate the list of observed outcomes
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
    # Get the theoretical total effect/reduced form/ITT
    total_effect <-
        (Y_1_1 - Y_0_0) * (D_1 == 1 & D_0 == 0) +
        (Y_1_1 - Y_0_1) * (D_1 == 1 & D_0 == 1) +
        (Y_1_0 - Y_0_0) * (D_1 == 0 & D_0 == 0)
    average_total_effect <- mean(total_effect)
    # Get the theoretical indirect effect.
    indirect_effect <-
        (Z * (Y_1_1 - Y_1_0) + (1 - Z) * (Y_0_1 - Y_0_0)) * (D_1 == 1 & D_0 == 0)
    average_indirect_effect <- mean(indirect_effect)
    # Get the theoretical direct effect.
    direct_effect <- #(D * (Y_1_1 - Y_0_1) + (1 - D) * (Y_1_0 - Y_0_0))
        (D * (Y_1_1 - Y_0_1) + (1 - D) * (Y_1_0 - Y_0_0)) * (D_1 == 1 & D_0 == 0) +
        (Y_1_1 - Y_0_1) * (D_1 == 1 & D_0 == 1) +
        (Y_1_0 - Y_0_0) * (D_1 == 0 & D_0 == 0)
    average_direct_effect <- mean(direct_effect)
    # Show the real underlying values.
    if (print.truth == TRUE){
        print("Here is a summary of the (unobserved) true effects:")
        # Show how many ATs, NTs, Compliers in terms of D_i(Z) for Z = 0, 1.
        print("How many compliers in the system?")
        print(table(D_1, D_0) / NROW(sim.data))
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
estimated.values <- function(firststage.reg, secondstage.reg, example.data){
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
    # calculate the second-stage indirect effect
    indirect.est <- predict(
        secondstage.reg, newdata = mutate(example.data, D = 1)) -
        predict(secondstage.reg, newdata = mutate(example.data, D = 0))
    # Add on the K_0, K_1 conditional on D_i = 0, 1 respectively for compliers.
    # Estimate the kappa-weight
    hat_probZ <- 0.5
    kappa_1 <- example.data$D * ((example.data$Z - hat_probZ) / (
        (1 - hat_probZ) * hat_probZ))
    kappa_0 <- (1 - example.data$D) * (((1 - example.data$Z) - (1 - hat_probZ)) / (
        (1 - hat_probZ) * hat_probZ))
    kappa.weight <- kappa_1 * hat_probZ + kappa_0 * (1 - hat_probZ)
    # Calculate the term to add on.
    add.term <- (weighted.mean(predict(secondstage.reg, newdata = mutate(
        filter(example.data, D == 1), Z = 0, D = 0, X_minus = 0, K_0 = 0)),
            kappa.weight[example.data$D == 1])
        - weighted.mean(predict(secondstage.reg, newdata = mutate(
            filter(example.data, D == 0), Z = 0, D = 0, X_minus = 0, K_1 = 0)),
                kappa.weight[example.data$D == 0]))
    indirect.est <- indirect.est + add.term
    # Return the mean estimates.
    output.list <- list(
        "first-stage"     = mean(firststage.est),
        "direct-effect"   = mean(direct.est),
        "indirect-effect" = mean(firststage.est * indirect.est))
    # Return the output.list
    return(output.list)
}

# Bootstrap the estimates.
estimated.loop <- function(boot.reps, example.data,
        sample.size = sample.N, bootstrap = TRUE){
    # Define lists the will be returned:
    # 2. Naive OLS.
    ols_direct_effect <- c()
    ols_indirect_effect <- c()
    # 3. Control function.
    cf_direct_effect <- c()
    cf_indirect_effect <- c()
    # More in the future ....
    ## Loop across the bootstraps values.
    for (i in seq(1, boot.reps)){
        # If bootstrapping, just resample from provided data.
        if (bootstrap == TRUE) {
            boot.indicies <- sample(seq(1, sample.size), sample.size, replace = TRUE)
            boot.data <- example.data[boot.indicies, ]
        }
        # If a regular re-simulation, get new data.
        else if (bootstrap == FALSE){
            boot.data <- simulate.data(0.5, 1, 2, 0.25)
            if ((100 * (i / boot.reps)) %% 1 == 0) {
                print(paste0(i, " out of ", boot.reps, ", ", 100 * (i / boot.reps), "% done."))
            }
        }
        else {stop("The `bootstrap' option only takes values of TRUE or FALSE.")}
        # Calculate the truth values, given the simulated data
        truth.est <- theoretical.values(example.data)
        # Now get the mediation effects, by different approaches.
        # 2. OLS estimate of second-stage
        ols_firststage.reg <- lm(D ~ (1 + Z) + X_minus + X_IV, data = boot.data)
        ols_secondstage.reg <- lm(Y ~ 1 + Z * D + X_minus, data = boot.data)
        ols.est <- estimated.values(ols_firststage.reg, ols_secondstage.reg, boot.data)
        # 3. Control Function estimates.
        cf_firststage.reg <- lm(D ~ (1 + Z) * X_IV *
            bs(X_minus, df = 20, intercept = TRUE),
            data = boot.data)
        boot.data$D_hat <- cf_firststage.reg$fitted
        control.fun <- ecdf(boot.data$D_hat)
        boot.data$K <- control.fun(boot.data$D_hat)
        boot.data$K_0 <- (1 - boot.data$D) * boot.data$K
        boot.data$K_1 <- boot.data$D * boot.data$K
        cf_secondstage.reg <- lm(
            Y ~ (1 + Z * D) + X_minus +
            bs(K_0, knots = seq(0, 1, by = 0.025), intercept = FALSE) +
            bs(K_1, knots = seq(0, 1, by = 0.025), intercept = FALSE),
            #bs(K, knots = seq(-1, 1, by = 0.025), intercept = FALSE),
            data = boot.data)
        cf.est <- estimated.values(cf_firststage.reg, cf_secondstage.reg, boot.data)
        # Save the outputs.
        ols_direct_effect[i]   <- ols.est$`direct-effect`
        ols_indirect_effect[i] <- ols.est$`indirect-effect`
        cf_direct_effect[i]    <- cf.est$`direct-effect`
        cf_indirect_effect[i]  <- cf.est$`indirect-effect`
    }
    # Return the bootstrap data.
    output.list <- list()
    output.list$data <- data.frame(
        truth_direct_effect   = as.numeric(truth.est$average_direct_effect),
        ols_direct_effect     = ols_direct_effect,
        cf_direct_effect      = cf_direct_effect,
        truth_indirect_effect = as.numeric(truth.est$average_indirect_effect),
        ols_indirect_effect   = ols_indirect_effect,
        cf_indirect_effect    = cf_indirect_effect
    )
    # Calculate the needed statistics, to return
    output.list$estimates <- data.frame(
        # Truth
        truth_direct_effect     = as.numeric(truth.est$average_direct_effect),
        truth_indirect_effect   = as.numeric(truth.est$average_indirect_effect),
        # OLS mean, and the 95% confidence intervals
        ols_direct_effect       = as.numeric(mean(ols_direct_effect)),
        ols_direct_effect_se    = as.numeric(sd(ols_direct_effect)),
        ols_direct_effect_up    = as.numeric(quantile(ols_direct_effect, probs = 0.975)),
        ols_direct_effect_low   = as.numeric(quantile(ols_direct_effect, probs = 0.025)),
        ols_indirect_effect     = as.numeric(mean(ols_indirect_effect)),
        ols_indirect_effect_se  = as.numeric(sd(ols_indirect_effect)),
        ols_indirect_effect_up  = as.numeric(quantile(ols_indirect_effect, probs = 0.975)),
        ols_indirect_effect_low = as.numeric(quantile(ols_indirect_effect, probs = 0.025)),
        # Control Fun mean, and the 95% confidence intervals
        cf_direct_effect        = as.numeric(mean(cf_direct_effect)),
        cf_direct_effect_se     = as.numeric(sd(cf_direct_effect)),
        cf_direct_effect_up     = as.numeric(quantile(cf_direct_effect, probs = 0.975)),
        cf_direct_effect_low    = as.numeric(quantile(cf_direct_effect, probs = 0.025)),
        cf_indirect_effect      = as.numeric(mean(cf_indirect_effect)),
        cf_indirect_effect_se   = as.numeric(sd(cf_indirect_effect)),
        cf_indirect_effect_up   = as.numeric(quantile(cf_indirect_effect, probs = 0.975)),
        cf_indirect_effect_low  = as.numeric(quantile(cf_indirect_effect, probs = 0.025))
    )
    return(output.list)
}


################################################################################
## Compare estimation methods, in one simulation.

## Simulate the data: rho, sigma_0, sigma_1, sigma_C = 0.5, 1, 2, 1.
simulated.data <- simulate.data(0.5, 1, 2, 0.25)
# SHow the theoretical direct + indirect values
print(theoretical.values(simulated.data, print.truth = TRUE))

# Show that the regression specification holds exactly (after debiasing outcome).
true_firststage.reg <- glm(D ~ (1 + Z) + X_IV +
    bs(X_minus, df = 5, intercept = TRUE),
    family = binomial(link = "probit"), data = simulated.data)
simulated.data$K_0 <- (1 - simulated.data$D) * simulated.data$U_0
simulated.data$K_1 <- simulated.data$D * simulated.data$U_1
true_secondstage.reg <- lm(Y ~ (1 + Z * D) + X_minus +
    # including the unobserved errors: U_0 + (D * U_0) + (D * U_1),
    K_0 + K_1,
    data = simulated.data)
print(theoretical.values(simulated.data))
print(estimated.values(true_firststage.reg, true_secondstage.reg, simulated.data))
# This is because the truth second-stage is perfect:
print(summary(true_secondstage.reg))

# Show how the OLS result gives a bias result (if rho != 0)
ols_firststage.reg <- lm(D ~ (1 + Z) + X_minus + X_IV, data = simulated.data)
ols_secondstage.reg <- lm(Y ~ 1 + Z * D + X_minus, data = simulated.data)
print(theoretical.values(simulated.data))
print(estimated.values(ols_firststage.reg, ols_secondstage.reg, simulated.data))

# Show how (unknown) control function gets it correct, in 2 steps (with splines)
# -> Use the estimated CDF as the control function (as in Imbens Newey 2009).
cf_firststage.reg <- lm(D ~ (1 + Z) * X_IV *
    bs(X_minus, df = 20, intercept = TRUE),
    data = simulated.data)
simulated.data$D_hat <- cf_firststage.reg$fitted
control.fun <- ecdf(simulated.data$D_hat)
simulated.data$K <- control.fun(simulated.data$D_hat)
simulated.data$K_0 <- (1 - simulated.data$D) * simulated.data$K
simulated.data$K_1 <- simulated.data$D * simulated.data$K
cf_secondstage.reg <- lm(Y ~ (1 + Z * D) + X_minus +
    bs(K_0, knots = seq(0, 1, by = 0.025), intercept = FALSE) +
    bs(K_1, knots = seq(0, 1, by = 0.025), intercept = FALSE),
    #bs(K, knots = seq(-1, 1, by = 0.05), intercept = FALSE),
    data = simulated.data)
print(summary(cf_secondstage.reg))
print(theoretical.values(simulated.data))
print(estimated.values(cf_firststage.reg, cf_secondstage.reg, simulated.data))

#! Test, add on the K_0 and K_1 conditional on D_i = 0,1 respectively in truth
firststage.est <- predict(
    true_firststage.reg, newdata = mutate(simulated.data, Z = 1), type = "response") - predict(
        true_firststage.reg, newdata = mutate(simulated.data, Z = 0), type = "response")
# calculate the second-stage indirect effect
indirect.est <- predict(
    true_secondstage.reg, newdata = mutate(simulated.data, D = 1)) -
    predict(true_secondstage.reg, newdata = mutate(simulated.data, D = 0))
add.term <- (mean(simulated.data$K_1[simulated.data$D_0 == 0 & simulated.data$D_1 == 1 & simulated.data$D == 1])
    - mean(simulated.data$K_0[simulated.data$D_0 == 0 & simulated.data$D_1 == 1 & simulated.data$D == 1]))
print(mean(firststage.est * (indirect.est + add.term)))
#! Do the same thing, but kappa weighted to the compliers inside the CF estimate.
firststage.est <- predict(
    cf_firststage.reg, newdata = mutate(simulated.data, Z = 1), type = "response") - predict(
        cf_firststage.reg, newdata = mutate(simulated.data, Z = 0), type = "response")
# calculate the second-stage indirect effect
indirect.est <- predict(
    cf_secondstage.reg, newdata = mutate(simulated.data, D = 1)) -
    predict(cf_secondstage.reg, newdata = mutate(simulated.data, D = 0))
# Estimate the kappa-weight
hat_probZ <- 0.5
kappa_1 <- simulated.data$D * ((simulated.data$Z - hat_probZ) / (
    (1 - hat_probZ) * hat_probZ))
kappa_0 <- (1 - simulated.data$D) * (((1 - simulated.data$Z) - (1 - hat_probZ)) / (
    (1 - hat_probZ) * hat_probZ))
kappa.weight <- kappa_1 * hat_probZ + kappa_0 * (1 - hat_probZ)
# Calculate the term to add on.
add.term <- (weighted.mean(simulated.data$K_1[simulated.data$D == 1],
        kappa.weight[simulated.data$D == 1])
    - weighted.mean(simulated.data$K_0[simulated.data$D == 0],
        kappa.weight[simulated.data$D == 0]))
mean(firststage.est * (indirect.est + add.term))

#! Test: show the ADE given Z_i = 1, and similar for AIE
firststage.est <- predict(
    true_firststage.reg, newdata = mutate(simulated.data, Z = 1), type = "response") - predict(
        true_firststage.reg, newdata = mutate(simulated.data, Z = 0), type = "response")
# calculate the second-stage direct effect
direct.est <- predict(
    true_secondstage.reg, newdata = mutate(simulated.data, Z = 1)) -
    predict(true_secondstage.reg, newdata = mutate(simulated.data, Z = 0))
# calculate the second-stage indirect effect
indirect.est <- predict(
    true_secondstage.reg, newdata = mutate(simulated.data, D = 1)) -
    predict(true_secondstage.reg, newdata = mutate(simulated.data, D = 0))
# Return the mean estimates.
print(theoretical.values(simulated.data))
print(c("first-stage",     mean(firststage.est)))
print(c("direct-effect",   mean(direct.est)))
print(c("indirect-effect", mean(firststage.est * indirect.est)))
# Get the treated versions
Z <- simulated.data$Z
D <- simulated.data$D
Y <- simulated.data$Y
D_0 <- simulated.data$D_0
D_1 <- simulated.data$D_1
Y_0_0 <- simulated.data$Y_0_0
Y_0_1 <- simulated.data$Y_0_1
Y_1_0 <- simulated.data$Y_1_0
Y_1_1 <- simulated.data$Y_1_1
# AIE, estimae vs the group differences term.
print(c("indirect-effect", mean(firststage.est * indirect.est)))
indirect_effect <- (Z * (Y_1_1 - Y_1_0) + (1 - Z) * (Y_0_1 - Y_0_0))
print(mean(indirect_effect * (D_1 != D_0)))
print(mean(indirect_effect * mean(D_1 != D_0)))
print(estimated.values(true_firststage.reg, true_secondstage.reg, simulated.data)$`indirect-effect`)

# ADE, estimate vs group differences term.
print(c("direct-effect",   mean(direct.est)))
direct_effect <- (D * (Y_1_1 - Y_0_1) + (1 - D) * (Y_1_0 - Y_0_0))
print(mean(direct_effect))
print(mean(direct_effect[Z == 1]))

#! Test: note the difference between AIE and LAIE (i.e., group differences term).
# show gains to D, on average
print(mean(simulated.data$Z * (simulated.data$Y_1_1 - simulated.data$Y_1_0) +
    (1 - simulated.data$Z) * (simulated.data$Y_0_1 - simulated.data$Y_0_0)))
# show gains to D, among compliers
print(mean((simulated.data$Z * (simulated.data$Y_1_1 - simulated.data$Y_1_0) +
    (1 - simulated.data$Z) * (simulated.data$Y_0_1 - simulated.data$Y_0_0)) * (
        simulated.data$D_1 == 1 & simulated.data$D_0 == 0)))

#! Test, get compliers correct
#! -> needs a different weighting scheme
#! -> Abadie (2003) kappa weights do not get it correct.
library(mgcv)
cf_firststage.reg <- lm(D ~ (1 + Z) * X_IV * bs(X_minus, df = 20, intercept = TRUE),
    data = simulated.data)
complier.weights <- predict(
    cf_firststage.reg, newdata = mutate(simulated.data, Z = 1), type = "response") - predict(
        cf_firststage.reg, newdata = mutate(simulated.data, Z = 0), type = "response")
complier.weights[complier.weights < 0] <- 0
complier_secondstage.reg <- lm(Y ~ (1 + Z * D) + X_minus +
    bs(K, knots = seq(-1, 1, by = 0.05)),
    weights = complier.weights,
    data = simulated.data)
print(theoretical.values(simulated.data))
print(estimated.values(cf_firststage.reg, complier_secondstage.reg, simulated.data))

#! Test: including the K_0 and K_1 terms in complier gains
simulated.data$K_0 <- (1 - simulated.data$D) * (1 - cf_firststage.reg$residuals)
simulated.data$K_1 <- simulated.data$D * cf_firststage.reg$residuals
complier_secondstage.reg <- lm(Y ~ (1 + Z * D) + X_minus +
    bs(K_0, knots = seq(-0.05, 1.05, by = 0.05)) +
    bs(K_1, knots = seq(-0.05, 1.05, by = 0.05)),
    data = simulated.data)

# Get the returns estimate correct.
print(mean((simulated.data$Y_0_1 - simulated.data$Y_0_0)[
        simulated.data$Z == 0]))
print(mean((simulated.data$Y_0_1 - simulated.data$Y_0_0)[
        simulated.data$Z == 0 & simulated.data$D_1 == 1 & simulated.data$D_0 == 0]))
print(mean((predict(
    cf_firststage.reg, newdata = mutate(simulated.data, Z = 1), type = "response") - predict(
        cf_firststage.reg, newdata = mutate(simulated.data, Z = 0), type = "response")) * (
    predict(complier_secondstage.reg,
        newdata = mutate(simulated.data, Z = 0, D = 1, K_0 = 0, K_1 = K_1)) -
    predict(complier_secondstage.reg,
        newdata = mutate(simulated.data, Z = 0, D = 0, K_0 = K_0, K_1 = 0)))))
print(mean(predict(complier_secondstage.reg,
        newdata = mutate(simulated.data, Z = 0, D = 1, K_0 = 0, K_1 = 0)) -
    predict(complier_secondstage.reg,
        newdata = mutate(simulated.data, Z = 0, D = 0, K_0 = 0, K_1 = 0))))

# calculate the second-stage indirect effect
firststage.est <- predict(
    cf_firststage.reg, newdata = mutate(simulated.data, Z = 1), type = "response") - predict(
        cf_firststage.reg, newdata = mutate(simulated.data, Z = 0), type = "response")
indirect.est <- predict(
    complier_secondstage.reg, newdata = mutate(simulated.data, D = 1, K_0 = 0)) -
        predict(complier_secondstage.reg, newdata = mutate(simulated.data, D = 0, K_1 = 0))
direct.est <- predict(
    complier_secondstage.reg, newdata = mutate(simulated.data, Z = 1)) -
        predict(complier_secondstage.reg, newdata = mutate(simulated.data, Z = 0))
# Show the means
print(theoretical.values(simulated.data))
print(mean(direct.est))
print(mean(firststage.est * indirect.est))
# Compared the un-weighted version (which has group differences bias).
print(estimated.values(cf_firststage.reg, complier_secondstage.reg, simulated.data))
print(estimated.values(cf_firststage.reg, cf_secondstage.reg, simulated.data))


################################################################################
## Plot bootstrap results for one DGP

# Base data to test out.
simulated.data <- simulate.data(0.5, 1, 2, 0.25)

# Get bootstrapped point est for the CF approach
boot.reps <- 10^2
boot.est <- estimated.loop(boot.reps, simulated.data, bootstrap = FALSE)
boot.data <- boot.est$data

## Save the bootstrapped point estimates data.
boot.data %>% write_csv(file.path(output.folder, "boot-sim-data-test.csv"))
exit.

################################################################################
## Compare estimation methods, across different sigma values.

# Define an empty dataframe, to start adding to.
boot.values <- estimated.loop(1, simulated.data)$estimates
boot.values$sigma <- NA
sigma.data <- boot.values[0, ]
print(sigma.data)
# Define values in rho \in [-1, 1] to go across
sigma.values <- seq(0, 2, by = 0.25)
# Define the number of boot reps for each
boot.reps <- 10^3
i <- 0

# Start the sigma loop
for (sigma in sigma.values){
    i <- i + 1
    # Simulate the data: rho, sigma_0, sigma_1, sigma_C
    sigma_sim.data <- simulate.data(0.5, 1, sigma, 0.25)
    # Get the truth + estimates + bootstrapped SEs, and save rho value
    sigma.boot <- estimated.loop(boot.reps, simulated.data,
        bootstrap = TRUE)$estimates
    sigma.boot$sigma <- sigma
    # Add to the dataframe.
    sigma.data[i, ] <- sigma.boot
    # SHow far we are.
    print(paste0(sigma, " in [0, 2], ", 100 * i / length(sigma.values), "% done."))
    gc()
}
# Save the output data.
sigma.data %>% write_csv(file.path(output.folder, "sigma-sim-data.csv"))


################################################################################
## Compare estimation methods, across different rho values.

# Define an empty dataframe, to start adding to.
boot.values <- estimated.loop(1, simulated.data)$estimates
boot.values$rho <- NA
rho.data <- boot.values[0, ]
print(rho.data)
# Define values in rho \in [-1, 1] to go across
rho.values <- seq(-1, 1, by = 0.25)
# Define the number of boot reps for each
boot.reps <- 10^3
i <- 0

# Start the rho loop
for (rho in rho.values){
    i <- i + 1
    # Simulate the data: rho, sigma_0, sigma_1, sigma_C
    rho_sim.data <- simulate.data(rho, 1, 2, 0.25)
    # Get the truth + estimates + bootstrapped SEs, and save rho value
    rho.boot  <- estimated.loop(boot.reps, simulated.data,
        bootstrap = TRUE)$estimates
    rho.boot$rho <- rho
    # Add to the dataframe.
    rho.data[i, ] <- rho.boot
    # Show how far we 
    print(paste0(rho, " in [-1, 1], ", 100 * i / length(rho.values), "% done."))
    gc()
}

## Save the output data.
rho.data %>% write_csv(file.path(output.folder, "rho-sim-data.csv"))
