#!/usr/bin/R
## Senan Hogan-Hennessy, 16 Jan 2025
## Simulate the system for indirect + direct effects, with Roy selection.
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
    ## sample.size: integer, representing output sample size (i.e., N).
    # First covariate (\vec X_i^-)
    X_minus <- 4 + rnorm(sample.size, mean = 0, sd = 1)
    # Second covariate (instrument for the control function).
    X_IV <- runif(sample.size, 0, 1) # rbinom(sample.size, 1, 1 / 2)
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
        (Z * (Y_1_1 - Y_1_0) + (1 - Z) * (Y_0_1 - Y_0_0)) * (D_1 == 1 & D_0 == 0)
    average_indirect_effect <- mean(indirect_effect)
    # Get the theoretical direct effect.
    direct_effect <- (D * (Y_1_1 - Y_0_1) + (1 - D) * (Y_1_0 - Y_0_0))
        #(D * (Y_1_1 - Y_0_1) + (1 - D) * (Y_1_0 - Y_0_0)) * (D_1 == 1 & D_0 == 0) +
        #(Y_1_1 - Y_0_1) * (D_1 == 1 & D_0 == 1) +
        #(Y_1_0 - Y_0_0) * (D_1 == 0 & D_0 == 0)
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
estimated.values <- function(
    firststage.reg, direct.reg, indirect.reg, example.data){
    ### Inputs:
    ## example.data, a data frame simulated from above.
    # calculate the first-stage by prediction
    firststage.est <- predict(
        firststage.reg, newdata = mutate(example.data, Z = 1), type = "response") - predict(
            firststage.reg, newdata = mutate(example.data, Z = 0), type = "response")
    # calculate the second-stage direct effect
    direct.est <- predict(
        direct.reg, newdata = mutate(example.data, Z = 1)) -
        predict(direct.reg, newdata = mutate(example.data, Z = 0))
    # calculate the second-stage indirect effect
    indirect.est <- predict(
        indirect.reg, newdata = mutate(example.data, D = 1)) -
        predict(indirect.reg, newdata = mutate(example.data, D = 0))
    # Add on the K_0, K_1 conditional on D_i = 0, 1 respectively for compliers.
    # # Estimate the kappa-weight
    # hat_probZ <- lm(Z ~ 1 + X_IV * bs(X_minus, df = 10, intercept = TRUE),
    #     data = example.data)$fitted
    # kappa_1 <- example.data$D * ((example.data$Z - hat_probZ) / (
    #     (1 - hat_probZ) * hat_probZ))
    # kappa_0 <- (1 - example.data$D) * (((1 - example.data$Z) - (1 - hat_probZ)) / (
    #     (1 - hat_probZ) * hat_probZ))
    # kappa.weight <- kappa_1 * hat_probZ + kappa_0 * (1 - hat_probZ)
    # # Calculate the term to add on.
    # K_0 <- weighted.mean(predict(secondstage.reg, newdata = mutate(
    #     example.data, Z = 0, D = 0, X_minus = 0)),
    #     kappa.weight * (1 - example.data$D), na.rm = TRUE)
    # K_1 <- weighted.mean(predict(secondstage.reg, newdata = mutate(
    #     example.data, Z = 0, D = 0, X_minus = 0)),
    #     kappa.weight * example.data$D, na.rm = TRUE)
    # indirect.est <- indirect.est + (K_1 - K_0)
    # Return the mean estimates.
    output.list <- list(
        "first-stage"     = mean(firststage.est, na.rm = TRUE),
        "direct-effect"   = mean(direct.est, na.rm = TRUE),
        "indirect-effect" = mean(firststage.est * indirect.est, na.rm = TRUE))
    # Return the output.list
    return(output.list)
}

# Define a function to cross-fit the semi-parametric control function.
cf_crossfit_mediate <- function(example.data){
    # 1. Calculate the split in half, two cross-fit samples
    example.size <- NROW(example.data)
    cross.index <- sample(seq(1, example.size),
        size = 0.5 * example.size, replace = FALSE)
    firstcross.data <- example.data[cross.index,]
    secondcross.data <- example.data[-cross.index,]
    # 2. calculate the CF model in the first cross sample.
    firstcross_firststage.reg <- lm(D ~ (1 + Z) * X_IV *
        bs(X_minus, df = 10, intercept = TRUE),
        data = firstcross.data)
    firstcross.data$K <- firstcross_firststage.reg$residuals
    #TODO: remove the separated direct + indirect effects.
    #firstcross_direct.reg <- lm(Y ~ (1 + Z * D) + X_minus +
    #    bs(K, knots = seq(-1, 1, by = 0.05), intercept = FALSE),
    #    data = firstcross.data)
    firstcross.data$K_0 <- (1 - firstcross.data$D) * firstcross.data$K
    firstcross.data$K_1 <- (firstcross.data$D) * firstcross.data$K
    firstcross_indirect.reg <- lm(Y ~ (1 + Z * D) + X_minus +
        #bs(K_0, knots = seq(0, 1, by = 0.05), intercept = TRUE) +
        #bs(K_1, knots = seq(0, 1, by = 0.05), intercept = TRUE),
        bs(K, knots = seq(-1, 1, by = 0.05), intercept = TRUE),
        data = firstcross.data)
    # 3. calculate the CF model on the second cross sample.
    secondcross_firststage.reg <- lm(D ~ (1 + Z) * X_IV *
        bs(X_minus, df = 10, intercept = TRUE),
        data = secondcross.data)
    secondcross.data$K <- secondcross_firststage.reg$residuals
    #secondcross_direct.reg <- lm(
    #    Y ~ (1 + Z * D) + X_minus +
    #    bs(K, knots = seq(-1, 1, by = 0.05), intercept = FALSE),
    #    data = secondcross.data)
    secondcross.data$K_0 <- (1 - secondcross.data$D) * secondcross.data$K
    secondcross.data$K_1 <- (secondcross.data$D) * secondcross.data$K
    secondcross_indirect.reg <- lm(
        Y ~ (1 + Z * D) + X_minus +
        #bs(K_0, knots = seq(0, 1, by = 0.05), intercept = TRUE) +
        #bs(K_1, knots = seq(0, 1, by = 0.05), intercept = TRUE),
        bs(K, knots = seq(-1, 1, by = 0.05), intercept = TRUE),
        data = secondcross.data)
    # 4. Predict the estimate on the opposite data.
    firstcross.est <- estimated.values(firstcross_firststage.reg,
        firstcross_indirect.reg, firstcross_indirect.reg, secondcross.data)
    secondcross.est <- estimated.values(secondcross_firststage.reg,
        secondcross_indirect.reg, secondcross_indirect.reg, firstcross.data)
    # Resturn the averaged estimates.
    output.list <- list(
        "first-stage"     = mean(c(
            firstcross.est$`first-stage`, secondcross.est$`first-stage`), na.rm = FALSE),
        "direct-effect"   = mean(c(
            firstcross.est$`direct-effect`, secondcross.est$`direct-effect`), na.rm = FALSE),
        "indirect-effect" = mean(c(
            firstcross.est$`indirect-effect`, secondcross.est$`indirect-effect`), na.rm = FALSE))
    # Return the output.list
    return(output.list)
}

# Bootstrap the estimates.
estimated.loop <- function(boot.reps, example.data,
    bootstrap = TRUE, print.progress == FALSE){
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
            boot.indicies <- sample(
                seq(1, NROW(example.data)), NROW(example.data), replace = TRUE)
            boot.data <- example.data[boot.indicies, ]
        }
        # If a regular re-simulation, get new data.
        else if (bootstrap == FALSE){
            boot.data <- simulate.data(0.5, 1, 2, 0.5)
        }
        else {stop("The `bootstrap' option only takes values of TRUE or FALSE.")}
        # Print, if want the consol output of how far we are.
        if (print.progress == TRUE){
            if ((100 * (i / boot.reps)) %% 1 == 0) {
                print(paste0(i, " out of ", boot.reps, ", ", 100 * (i / boot.reps), "% done."))
            }
        }
        # Calculate the truth values, given the simulated data
        truth.est <- theoretical.values(example.data)
        # Now get the mediation effects, by different approaches.
        # 2. OLS estimate of second-stage
        ols_firststage.reg <- lm(D ~ (1 + Z) + X_minus + X_IV, data = boot.data)
        ols_secondstage.reg <- lm(Y ~ 1 + Z * D + X_minus, data = boot.data)
        ols.est <- estimated.values(ols_firststage.reg,
            ols_secondstage.reg, ols_secondstage.reg, boot.data)
        # 3. Control Function estimates.
        cf.est <- cf_crossfit_mediate(boot.data)
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
            probs = 0.025, na.rm = TRUE))
    )
    return(output.list)
}


################################################################################
## Compare estimation methods, in one simulation.

## Simulate the data: rho, sigma_0, sigma_1, sigma_C = 0.5, 1, 2, 0.5
simulated.data <- simulate.data(0.5, 1, 2, 0.5)
# SHow the theoretical direct + indirect values
print(theoretical.values(simulated.data, print.truth = TRUE))

# Show that the regression specification holds exactly (after debiasing outcome).
true_firststage.reg <- glm(D ~ (1 + Z) + X_IV + X_minus +
    U_C, #U_0 + U_1,
    family = binomial(link = "probit"),
    data = simulated.data)
true_secondstage.reg <- lm(Y ~ (1 + Z * D) + X_minus +
    # including the unobserved errors:
    U_0 + (D * U_0) + (D * U_1),
    data = simulated.data)
print(theoretical.values(simulated.data))
print(estimated.values(true_firststage.reg,
    true_secondstage.reg, true_secondstage.reg, simulated.data))
# See how the first and second-stages are perfect:
print(summary(true_firststage.reg))
print(summary(true_secondstage.reg))

# Show how the OLS result gives a bias result (if rho != 0),
# using the same second-stage for both direct + indirect estimates.
ols_firststage.reg <- lm(D ~ (1 + Z) + X_minus + X_IV, data = simulated.data)
ols_secondstage.reg <- lm(Y ~ 1 + Z * D + X_minus, data = simulated.data)
print(theoretical.values(simulated.data))
print(estimated.values(ols_firststage.reg,
    ols_secondstage.reg, ols_secondstage.reg, simulated.data))

# Show how (unknown) control function gets it correct, in 2 steps (with splines)
cf_firststage.reg <- lm(D ~ (1 + Z) * X_IV  *
    bs(X_minus, df = 10, intercept = TRUE),
    data = simulated.data)
# Direct estimate.
simulated.data$K <- cf_firststage.reg$residuals
cf_direct.reg <- lm(Y ~ (1 + Z * D) + X_minus +
    bs(K, knots = seq(-1, 1, by = 0.05), intercept = FALSE),
    data = simulated.data)
# Indirect estimate.
simulated.data$K_0 <- (1 - simulated.data$D) * simulated.data$K
simulated.data$K_1 <- (simulated.data$D) * simulated.data$K
cf_indirect.reg <- lm(Y ~ (1 + Z * D) + X_minus +
    bs(K_0, knots = seq(0, 1, by = 0.05), intercept = TRUE) +
    bs(K_1, knots = seq(0, 1, by = 0.05), intercept = TRUE),
    data = simulated.data)
print(theoretical.values(simulated.data))
print(estimated.values(cf_firststage.reg,
    cf_direct.reg, cf_indirect.reg, simulated.data))
# Cross-fit the CF approach to avoid over-fitting bias in the semi-para step.
print(theoretical.values(simulated.data))
print(cf_crossfit_mediate(simulated.data))

#! Test: Imbens Newey (2009) conditional CDF as the control function.
cf_firststage.reg <- lm(D ~ (1 + Z) * X_IV *
    bs(X_minus, df = 20, intercept = TRUE),
    data = simulated.data)
control.fun <- ecdf(cf_firststage.reg$fitted)
simulated.data$K <- as.numeric(control.fun(cf_firststage.reg$fitted))
cdf_direct.reg <- lm(Y ~ (1 + Z * D) + X_minus +
    bs(K, knots = seq(0, 1, by = 0.05), intercept = TRUE),
    data = simulated.data)
simulated.data$K_0 <- (1 - simulated.data$D) * (1 - simulated.data$K)
simulated.data$K_1 <- simulated.data$D * simulated.data$K
cdf_indirect.reg <- lm(Y ~ (1 + Z * D) + X_minus +
    bs(K_0, knots = seq(0, 1, by = 0.05), intercept = TRUE)+
    bs(K_1, knots = seq(0, 1, by = 0.05), intercept = TRUE),
    data = simulated.data)
print(theoretical.values(simulated.data))
print(estimated.values(cf_firststage.reg,
    cdf_direct.reg, cdf_indirect.reg, simulated.data))

#! Test, add on the K_0 and K_1 conditional on D_i = 0,1 respectively in truth
firststage.est <- predict(
    true_firststage.reg, newdata = mutate(simulated.data, Z = 1), type = "response") - predict(
        true_firststage.reg, newdata = mutate(simulated.data, Z = 0), type = "response")
# calculate the second-stage indirect effect
indirect.est <- predict(
    true_secondstage.reg, newdata = mutate(simulated.data, D = 1)) -
    predict(true_secondstage.reg, newdata = mutate(simulated.data, D = 0))
add.term <- weighted.mean(simulated.data$U_1 - simulated.data$U_0,
    simulated.data$D_0 == 0 & simulated.data$D_1 == 1)
mean(firststage.est * (indirect.est + add.term))
print(theoretical.values(simulated.data)$average_indirect_effect)

#! Do the same thing, but kappa weighted to the compliers inside the CF estimate.
# Estimate the kappa-weight.
hat_probZ <- lm(Z ~ 1 * X_IV * bs(X_minus, df = 10, intercept = TRUE),
    data = simulated.data)$fitted
kappa_1 <- simulated.data$D * ((simulated.data$Z - hat_probZ) / (
    (1 - hat_probZ) * hat_probZ))
kappa_0 <- (1 - simulated.data$D) * (((1 - simulated.data$Z) - (1 - hat_probZ)) / (
    (1 - hat_probZ) * hat_probZ))
kappa.weight <- kappa_1 * hat_probZ + kappa_0 * (1 - hat_probZ)
# Calculate the term to add on.
errors <- predict(cdf_direct.reg, newdata = mutate(simulated.data,
    Z = 0, D = 0, X_minus = 0))
add.term <- weighted.mean(errors, simulated.data$D * kappa.weight) -
    weighted.mean(errors, (1 - simulated.data$D) * kappa.weight)
add.term <- weighted.mean(simulated.data$D * errors -
    (1 - simulated.data$D) * errors, kappa.weight)
mean(firststage.est * (indirect.est + add.term))


################################################################################
## Plot bootstrap results for one DGP, repeatedly drawn

# Base data to test out.
simulated.data <- simulate.data(0.5, 1, 2, 0.5)

# Get bootstrapped point est for the CF approach
sim.reps <- 10^4
sim.est <- estimated.loop(sim.reps, simulated.data,
    bootstrap = FALSE, print.progress = TRUE)
sim.data <- sim.est$data

## Save the repated DGPs' point estimates as separate data.
sim.data %>% write_csv(file.path(output.folder, "boot-sim-data.csv"))


################################################################################
## Compare estimation methods, across different sigma values.

# Define an empty dataframe, to start adding to.
boot.values <- estimated.loop(1, simulated.data)$estimates
boot.values$sigma <- NA
sigma.data <- boot.values[0, ]
# Define values in rho \in [-1, 1] to go across
sigma.values <- seq(0, 2, by = 0.25)
# Define the number of boot reps for each
boot.reps <- 10^3
i <- 0

# Start the sigma loop
for (sigma in sigma.values){
    # Simulate the data: rho, sigma_0, sigma_1, sigma_C
    sigma_sim.data <- simulate.data(0.5, sigma, 2 * sigma, 0.5)
    # Get the truth + estimates + bootstrapped SEs, and save rho value
    sigma.boot <- estimated.loop(boot.reps, sigma_sim.data,
        bootstrap = FALSE)$estimates
    sigma.boot$sigma <- sigma
    # Add to the dataframe.
    i <- i + 1
    sigma.data[i, ] <- sigma.boot
    # SHow far we are.
    print(paste0("simga = ", sigma,
        " in [0, 2], ", 100 * i / length(sigma.values), "% done."))
    gc()
}
# Save the output data.
sigma.data %>% write_csv(file.path(output.folder, "sigma-sim-data.csv"))


################################################################################
## Compare estimation methods, across different sigma_1 values.

# Define an empty dataframe, to start adding to.
boot.values <- estimated.loop(1, simulated.data)$estimates
boot.values$sigma_1 <- NA
sigma_1.data <- boot.values[0, ]
# Define values in rho \in [-1, 1] to go across
sigma_1.values <- seq(0, 2, by = 0.25)
# Define the number of boot reps for each
boot.reps <- 10^3
i <- 0

# Start the sigma_1 loop
for (sigma_1 in sigma_1.values){
    i <- i + 1
    # Simulate the data: rho, sigma_0, sigma_1, sigma_C
    sigma_1_sim.data <- simulate.data(0.5, 1, sigma_1, 0.5)
    # Get the truth + estimates + bootstrapped SEs, and save rho value
    sigma_1.boot <- estimated.loop(boot.reps, sigma_1_sim.data,
        bootstrap = FALSE)$estimates
    sigma_1.boot$sigma_1 <- sigma_1
    # Add to the dataframe.
    sigma_1.data[i, ] <- sigma_1.boot
    # SHow far we are.
    print(paste0("sigma_1 = ", sigma_1,
        " in [0, 2], ", 100 * i / length(sigma_1.values), "% done."))
    gc()
}
# Save the output data.
sigma_1.data %>% write_csv(file.path(output.folder, "sigma-1-sim-data.csv"))


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
    rho_sim.data <- simulate.data(rho, 1, 2, 0.5)
    # Get the truth + estimates + bootstrapped SEs, and save rho value
    rho.boot  <- estimated.loop(boot.reps, rho_sim.data,
        bootstrap = FALSE)$estimates
    rho.boot$rho <- rho
    # Add to the dataframe.
    rho.data[i, ] <- rho.boot
    # Show how far we 
    print(paste0("rho = ", rho,
        " in [-1, 1], ", 100 * i / length(rho.values), "% done."))
    gc()
}

## Save the output data.
rho.data %>% write_csv(file.path(output.folder, "rho-sim-data.csv"))
