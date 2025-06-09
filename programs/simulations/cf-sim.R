#!/usr/bin/R
## Senan Hogan-Hennessy, 6 May 2025
## Identifying ADE + AIE with a sem-parametric control function
## see Hogan-Hennessy (2025).

# Show the date:
print(format(Sys.time(), "%H:%M %Z %A, %d %B %Y"))

## Load libraries
# Functions for data manipulation and visualisation
library(tidyverse)
# Package for more distributions to sample from.
library(MASS)
# Package forsemi-parametric CF, motivated by MTEs.
library(ivmte)
# Semi-parametric binary choice model (for mediator first-stage).
# (Klein Spady 1993, https://github.com/henrykye/semiBRM)
library(semiBRM)


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
sample.N <- 10^3


################################################################################
## Define a function to simulate data in the triangular system.

# Define a function to simulate all observed + unobserved data 
roy.data <- function(rho, sigma_0, sigma_1, sigma_C,
        error.dist = "normal", sample.size = sample.N){
    ### Inputs:
    ## X, a matrix of covariates, continuous or binary values.
    ## rho \in [-1, +1] measuring correlation between U_0, U_1.
    ## sigma_0 >= 0 measuring standard deviation of U_0.
    ## sigma_1 >= 0 measuring standard deviation of U_1.
    ## sigma_C >= 0 measuring standard deviation of U_C.
    ## sample.size: integer, representing output sample size (i.e., N).
    # First covariate, \vec X_i^-
    X_minus <- 4 + rnorm(sample.size, mean = 0, sd = 1)
    # Second covariate (instrument for the control function).
    X_IV <- runif(sample.size, -1, 1)
    # Simulate the unobserved error terms.
    mu.vector <- c(0, 0, 0)
    sigma.matrix <- matrix(c(
        sigma_0^2,               rho * sigma_0 * sigma_1, 0,
        rho * sigma_0 * sigma_1, sigma_1^2,               0,
        0,                       0,                       sigma_C^2),
        ncol = 3)
    # Simulate a joint normal, from which a prob integral transform can get
    # any other (potentially correlated) distribution. 
    U_all <- mvrnorm(n = sample.size, mu = mu.vector, sigma.matrix,
        empirical = FALSE)
    # Get the desired distribution.
    if (error.dist == "normal"){
        U_0 <- U_all[, 1]
        U_1 <- U_all[, 2]
        U_C <- U_all[, 3]
    } 
    else if (error.dist == "uniform"){
        # uniform errors (centred on zero), with SD = sigma_0,1,C.
        U_0 <- (2 * pnorm(U_all[, 1], mean = 0, sd = sigma_0) - 1) * sigma_0 * sqrt(3)
        U_1 <- (2 * pnorm(U_all[, 2], mean = 0, sd = sigma_1) - 1) * sigma_1 * sqrt(3)
        U_C <- (2 * pnorm(U_all[, 3], mean = 0, sd = sigma_C) - 1) * sigma_C * sqrt(3)
    }
    else { stop("Choose an error.dist from `normal`, `uniform`.")
    }
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
        print(paste0(c("Cov(V, U_1):", cov(U_C - U_1 + U_0, U_0))))
        print(paste0(c("Cov(V, U_0):", cov(U_C - U_1 + U_0, U_1))))
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
estimated.values <- function(firststage.reg, secondstage.reg, totaleffect.reg,
    example.data, complier.adjustment = NULL){
    ### Inputs:
    ## example.data, a data frame simulated from above.
    # calculate the first-stage by prediction
    firststage.est <- predict(
        firststage.reg, newdata = mutate(example.data, Z = 1),
            type = "response") - predict(
                firststage.reg, newdata = mutate(example.data, Z = 0),
                    type = "response")
    # Calculate the total effect estimate by prediction.
    totaleffect.est <- predict(
        totaleffect.reg, newdata = mutate(example.data, Z = 1)) - predict(
            totaleffect.reg, newdata = mutate(example.data, Z = 0))
    # calculate the second-stage direct effect
    direct.est <- predict(
        secondstage.reg, newdata = mutate(example.data, Z = 1)) - predict(
            secondstage.reg, newdata = mutate(example.data, Z = 0))
    # calculate the second-stage (controlled) indirect effect
    indirect.est <- predict(
        secondstage.reg, newdata = mutate(example.data, D = 1)) - predict(
            secondstage.reg, newdata = mutate(example.data, D = 0))
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


################################################################################
## Two mediation estimation functions (1) Parametric CF (2) Semi-parametric CF.

# Define a function to Heckman selection correct mediation est, in two-stages.
mediate.heckit <- function(example.data){
    # 0. Total effect regression.
    totaleffect.reg <- lm(Y ~ 1 + Z + X_minus,
        data = example.data)
    # 1. Probit first-stage (well identified).
    heckit_firststage.reg <- glm(D ~ (1 + Z) + X_IV + X_minus,
        family = binomial(link = "probit"),
        data = example.data)
    # 2. Define the CFs --- for assumed N(0,1) dist.
    lambda_1.fun <- function(pi.est){
        # Inv Mills ratio, taking as input the estimated mediator propensity.
        return(dnorm(qnorm(pi.est)) / pnorm(qnorm(pi.est)))
    }
    pi.est <- predict(heckit_firststage.reg, type = "response")
    example.data$lambda_0 <- (1 - example.data$D) * lambda_1.fun(pi.est) * (
        - pi.est / (1 - pi.est))
    example.data$lambda_1 <- example.data$D * lambda_1.fun(pi.est)
    # 3. Estimate second-stage, including the CFs.
    heckit_secondstage.reg <- lm(Y ~ (1 + Z * D) + X_minus + lambda_0 + lambda_1,
        data = example.data)
    # Compensate complier difference in AIE, by Kline Walters (2019) IV-type adjustment.
    pi_0.est <- predict(heckit_firststage.reg, newdata = mutate(example.data, Z = 0), type = "response")
    pi_1.est <- predict(heckit_firststage.reg, newdata = mutate(example.data, Z = 1), type = "response")
    Gamma.big <-  (pi_1.est * lambda_1.fun(pi_1.est) - pi_0.est * lambda_1.fun(pi_0.est)) / (
        pi_1.est - pi_0.est)
    rho_0 <- coef(summary(heckit_secondstage.reg))["lambda_0", "Estimate"]
    rho_1 <- coef(summary(heckit_secondstage.reg))["lambda_1", "Estimate"]
    complier.adjustment <- (rho_1 - rho_0) * Gamma.big
    # Compile the estimates.
    heckit.est <- estimated.values(heckit_firststage.reg, heckit_secondstage.reg,
        totaleffect.reg, example.data,
        complier.adjustment = complier.adjustment)
    # Return the off-setting estimates.
    output.list <- list(
        "first-stage"     = heckit.est$"first-stage",
        "total-effect"    = heckit.est$"total-effect",
        "direct-effect"   = heckit.est$"direct-effect",
        "indirect-effect" = heckit.est$"indirect-effect")
    return(output.list)
}

# Define a function to two-stage semi-parametric CF for CM effects.
mediate.semiparametric <- function(example.data){
    #example.data <- simulated.data
    # 0. Total effect regression.
    totaleffect.reg <- lm(Y ~ 1 + Z + X_minus,
        data = example.data)
    totaleffect.est <- mean(predict(
        totaleffect.reg, newdata = mutate(example.data, Z = 1)) - predict(
            totaleffect.reg, newdata = mutate(example.data, Z = 0)))
    # 2. Use ivmte for semi-parametric second-stage, with a polynomial.
    example.data$pi_0 <- mean(example.data$D[example.data$Z == 0])
    example.data$pi_1 <- mean(example.data$D[example.data$Z == 1])
    cf_secondstage.reg <- ivmte(
        propensity = D ~ (1 + Z) + X_IV + X_minus,
        outcome =  "Y",
        m0 = ~ 1 + Z + X_minus + u + I(u^2) + I(u^3) + I(u^4) + I(u^5),
        m1 = ~ 1 + Z + X_minus + u + I(u^2) + I(u^3) + I(u^4) + I(u^5),
        target = "late", late.from = c(Z = 0), late.to = c(Z = 1),
        # target = "genlate", genlate.lb = pi_0, genlate.ub = pi_0,
        solver = "rmosek",
        point = TRUE,
        data = example.data)
    # 3. Compose the CM effects from this object.
    D_0 <- 1 - mean(example.data$D)
    D_1 <- 1 - D_0
    pi.bar <- (mean(cf_secondstage.reg$propensity$phat[example.data$Z == 1])
        - mean(cf_secondstage.reg$propensity$phat[example.data$Z == 0]))
    # 3.1 ADE point estimate, from the CF model.
    ade.est <- as.numeric(D_0 * cf_secondstage.reg$mtr.coef["[m0]Z"]
        + D_1 * cf_secondstage.reg$mtr.coef["[m1]Z"])
    # 3.2 AIE estimate by the IVMTE extrapolation across semi-parametric CFs.
    # aie.est <- pi.bar * cf_secondstage.reg$point.estimate
    # 3.2.1 AIE by using ADE estimate, relative to ATE.
    # (Avoiding semi-parametric extrapolation, see notes on ATE comparison)
    gammma.est <- cf_secondstage.reg$mtr.coef["[m0]Z"]
    delta.est <- cf_secondstage.reg$mtr.coef["[m1]Z"] - gammma.est
    ade_Z0.est <- gammma.est + delta.est * mean(
        example.data$D[example.data$Z == 0])
    ade_Z1.est <- gammma.est + delta.est * mean(
        example.data$D[example.data$Z == 1])
    aie.est <- totaleffect.est - mean(
        (1 - example.data$Z) * ade_Z1.est + example.data$Z * ade_Z0.est)
    # Return the estimates.
    output.list <- list(
        "first-stage"     = pi.bar,
        "total-effect"    = totaleffect.est,
        "direct-effect"   = ade.est,
        "indirect-effect" = aie.est)
    return(output.list)
}


################################################################################
## Define a function to automate the DGP generation, or bootstrapping models.

# Define a loop.
estimated.loop <- function(boot.reps, example.data,
        bootstrap = TRUE, print.progress = FALSE, error.dist = "normal",
        # DGP parameters
        rho, sigma_0, sigma_1, sigma_C) {
    # Define lists the will be returned:
    # 1. Naive OLS (fine for total effect, but not CM effects).
    ols_total_effect <- c()
    ols_direct_effect <- c()
    ols_indirect_effect <- c()
    # 2. Heckman selection model adjustment
    heckit_direct_effect <- c()
    heckit_indirect_effect <- c()
    # 3. Semi-parametric Control function.
    cf_direct_effect <- c()
    cf_indirect_effect <- c()
    # Calculate the truth values, given the input data
    truth_total_effect <- c()
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
        # If a regular DGP re-simulation, get new data.
        else if (bootstrap == FALSE){
            boot.data <- roy.data(rho, sigma_0, sigma_1, sigma_C,
                error.dist = error.dist)
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
        # First, get the total effect, by a standard regression.
        totaleffect.reg <- lm(Y ~ 1 + Z + X_minus, data = boot.data)
        # Now get the mediation effects, by different approaches.
        # 2. OLS estimate of second-stage
        ols_firststage.reg <- lm(D ~ (1 + Z) + X_minus + X_IV, data = boot.data)
        ols_secondstage.reg <- lm(Y ~ 1 + Z * D + X_minus, data = boot.data)
        ols.est <- estimated.values(ols_firststage.reg, ols_secondstage.reg,
            totaleffect.reg, boot.data)
        # 3. Heckman-style selection-into-mediator model estimates.
        heckit.est <- mediate.heckit(boot.data)
        # 4. CF / MTE semi-parametric selection-into-mediator model estimates.
        cf.est <- mediate.semiparametric(boot.data)
        # Save the outputs.
        truth_total_effect[i]     <- truth.est$average_total_effect
        truth_direct_effect[i]    <- truth.est$average_direct_effect
        truth_indirect_effect[i]  <- truth.est$average_indirect_effect
        ols_total_effect[i]       <- ols.est$"total-effect"
        ols_direct_effect[i]      <- ols.est$"direct-effect"
        ols_indirect_effect[i]    <- ols.est$"indirect-effect"
        heckit_direct_effect[i]   <- heckit.est$"direct-effect"
        heckit_indirect_effect[i] <- heckit.est$"indirect-effect"
        cf_direct_effect[i]       <- cf.est$"direct-effect"
        cf_indirect_effect[i]     <- cf.est$"indirect-effect"
    }
    # Return the bootstrap data.
    output.list <- list()
    output.list$data <- data.frame(
        # Total effect
        truth_total_effect     = truth_total_effect,
        ols_total_effect       = ols_total_effect,
        # Direct effect
        truth_direct_effect    = truth_direct_effect,
        ols_direct_effect      = ols_direct_effect,
        heckit_direct_effect   = heckit_direct_effect,
        cf_direct_effect       = cf_direct_effect,
        # Indirect effect
        truth_indirect_effect  = truth_indirect_effect,
        ols_indirect_effect    = ols_indirect_effect,
        heckit_indirect_effect = heckit_indirect_effect,
        cf_indirect_effect     = cf_indirect_effect)
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
        heckit_direct_effect        = as.numeric(mean(heckit_direct_effect)),
        heckit_direct_effect_se     = as.numeric(sd(heckit_direct_effect)),
        heckit_direct_effect_up     = as.numeric(quantile(heckit_direct_effect,
            probs = 0.975, na.rm = TRUE)),
        heckit_direct_effect_low    = as.numeric(quantile(heckit_direct_effect,
            probs = 0.025, na.rm = TRUE)),
        heckit_indirect_effect      = as.numeric(mean(heckit_indirect_effect)),
        heckit_indirect_effect_se   = as.numeric(sd(heckit_indirect_effect)),
        heckit_indirect_effect_up   = as.numeric(quantile(heckit_indirect_effect,
            probs = 0.975, na.rm = TRUE)),
        heckit_indirect_effect_low  = as.numeric(quantile(heckit_indirect_effect,
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
sigma_1 <- 1.5
sigma_C <- 0.5
simulated.data <- roy.data(rho, sigma_0, sigma_1, sigma_C)
# SHow the theoretical direct + indirect values
print(theoretical.values(simulated.data, print.truth = TRUE))

# Show that the regression specification holds exactly, if specified correctly.
true_firststage.reg <- glm(D ~ (1 + Z) + X_IV + X_minus + U_C + U_0 + U_1,
    family = binomial(link = "probit"),
    data = simulated.data)
true_secondstage.reg <- lm(Y ~ (1 + Z * D) + X_minus +
    # including the unobserved errors:
    D * U_0 + D * U_1,
    data = simulated.data)
true_total.reg <- lm(Y ~ (1 + Z) + X_minus, data = simulated.data)
print(theoretical.values(simulated.data))
print(estimated.values(true_firststage.reg, true_secondstage.reg,
    true_total.reg, simulated.data))
# See how the first and second-stages are perfect:
print(summary(true_firststage.reg))
print(summary(true_secondstage.reg))
#! Note: this automatically includes the complier term in the indirect effect
#!       if U_0, U_1 unobserved, this is not automatic.

# Show how the OLS result gives a bias result (if rho != 0),
# using the same second-stage for both direct + indirect estimates.
ols_firststage.reg <- lm(D ~ (1 + Z) + X_minus + X_IV, data = simulated.data)
ols_secondstage.reg <- lm(Y ~ 1 + Z * D + X_minus, data = simulated.data)
print(theoretical.values(simulated.data))
print(estimated.values(ols_firststage.reg, ols_secondstage.reg,
    true_total.reg, simulated.data))

# Show how the CF approaches get it correct (with complier adjustment).
print(theoretical.values(simulated.data))
print(mediate.heckit(simulated.data))
print(mediate.semiparametric(simulated.data))

# One-shot showing how the uniform errors screw up the heckman selection model.
uniform.data <- roy.data(rho, sigma_0, sigma_1, sigma_C, error.dist = "uniform")
print(theoretical.values(uniform.data))
print(mediate.heckit(uniform.data))
print(mediate.semiparametric(uniform.data))


################################################################################
## Plot DGP-strap results for normally distributed errors, repeatedly drawn
# Note, using the same input parameters as before, to keep comparable.

# First, show us what huge sample gives for the truth values.
roy.data(rho, sigma_0, sigma_1, sigma_C, sample.size = 10^6,
    error.dist = "normal") %>%
    theoretical.values(print.truth = TRUE) %>%
    print()

# Base data to input.
normal.data <- roy.data(rho, sigma_0, sigma_1, sigma_C, error.dist = "normal")

# Get bootstrapped point est for the CF approach
sim.reps <- 10^4
normal.est <- estimated.loop(sim.reps, normal.data,
    bootstrap = FALSE, print.progress = TRUE,  error.dist = "normal",
    rho, sigma_0, sigma_1, sigma_C)

# Save the repeated DGPs' point estimates as separate data.
normal.est$data %>% write_csv(file.path(output.folder, "normal-cf-data.csv"))


################################################################################
## Plot DGP-strap results for uniform distributed errors, repeatedly drawn
# Note, using the same input parameters as before, to keep comparable.

# Base data to input.
uniform.data <- roy.data(rho, sigma_0, sigma_1, sigma_C, error.dist = "uniform")

# Get bootstrapped point est for the CF approach
uniform.est <- estimated.loop(sim.reps, uniform.data,
    bootstrap = FALSE, print.progress = TRUE, error.dist = "uniform",
    rho, sigma_0, sigma_1, sigma_C)

# Save the repeated DGPs' point estimates as separate data.
uniform.est$data %>% write_csv(file.path(output.folder, "uniform-cf-data.csv"))


################################################################################
## Compare estimation methods, across different rho values.

# Define an empty dataframe, to start adding to.
rho_normal.data <- roy.data(rho, sigma_0, sigma_1, sigma_C,
    error.dist = "normal")
rho.data <- estimated.loop(1, rho_normal.data,
    bootstrap = TRUE, print.progress = FALSE, error.dist = "normal",
    rho, sigma_0, sigma_1, sigma_C)$estimates
rho.values$rho <- NA
rho.data <- rho.values$data[0, ]
# Define values in rho \in [-1, 1] to go across
rho.values <- seq(-1, 1, by = 0.25)
# Define the number of boot reps for each
boot.reps <- 1000
i <- 0

# Start the rho loop
for (rho.value in rho.values){
    i <- i + 1
    # Simulate the data: rho, sigma_0, sigma_1, sigma_C
    rho_sim.data <- roy.data(rho.value, sigma_0, sigma_1, sigma_C,
        error.dist = "normal")
    # Get the truth + estimates + bootstrapped SEs, and save rho value
    rho.boot <- estimated.loop(boot.reps, rho_sim.data,
        bootstrap = TRUE, print.progress = FALSE, error.dist = "normal",
        rho.value, sigma_0, sigma_1, sigma_C)$estimates
    # Add to the dataframe.
    rho.boot$rho <- rho.value
    rho.data[i, ] <- rho.boot
    # SHow far we are.
    print(paste0(rho.value, " in [-1, 1], ",
        100 * i / length(rho.values), "% done."))
    gc()
}

## Save the output data.
rho.data %>% write_csv(file.path(output.folder, "rho-cf-data.csv"))


################################################################################
## Compare estimation methods, across different sd(U_1) values.

# Define an empty dataframe, to start adding to.
sigma_1_normal.data <- roy.data(rho, sigma_0, sigma_1, sigma_C,
    error.dist = "normal")
sigma_1.values <- estimated.loop(1, sigma_1_normal.data,
    bootstrap = TRUE, print.progress = FALSE, error.dist = "normal",
    rho, sigma_0, sigma_1, sigma_C)
sigma_1.values$sigma_1 <- NA
sigma_1.data <- sigma_1.values$data[0, ]
# Define values in sigma_1 in [0, 2] to go across
sigma_1.values <- seq(0, 2, by = 0.25)
# Define the number of boot reps for each
boot.reps <- 10^3
i <- 0

# Start the sigma_1 loop
for (sigma_1.value in sigma_1.values){
    i <- i + 1
    # Simulate the data: rho, sigma_0, sigma_1, sigma_C
    sigma_1_sim.data <- roy.data(rho, sigma_0, sigma_1.value, sigma_C,
        error.dist = "normal")
    # Get the truth + estimates + bootstrapped SEs, and save rho value
    sigma_1.boot <- estimated.loop(boot.reps, sigma_1_sim.data,
        bootstrap = TRUE, print.progress = FALSE, error.dist = "normal",
        rho, sigma_0, sigma_1.value, sigma_C)$estimates
    # Add to the dataframe.
    sigma_1.boot$sigma_1 <- sigma_1.value
    sigma_1.data[i, ] <- sigma_1.boot
    # Show how far we are.
    print(paste0(sigma_1.value, " in [0, 2], ",
        100 * i / length(sigma_1.values), "% done."))
    gc()
}

## Save the output data.
sigma_1.data %>% write_csv(file.path(output.folder, "sigma1-cf-data.csv"))
