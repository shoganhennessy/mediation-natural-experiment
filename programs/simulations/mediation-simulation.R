#!/usr/bin/R
## Senan Hogan-Hennessy, 04 October 2023
## Simulate the system for indirect + direct effects
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

## Set up the R environment
set.seed(47)
# Define number of digits in tables and graphs
digits.no <- 3
# Define where output files go.
output.folder <- file.path("sim-output")
# Set the options for the plot sizes, in saving ggplot output.
fig.height <- 10
fig.width <- fig.height


################################################################################
## Define a function to simulate data in the triangular system.
simulate.data <- function(rho, sigma_0, sigma_1, sigma_C, sample.size = 10^4){
    ### Inputs:
    ## X, a matrix of covariates, continuous or binary values.
    ## rho \in [-1, +1] measuring correlation between U_0, U_1.
    ## sigma_0 >= 0 measuring standard deviation of U_0.
    ## sigma_1 >= 0 measuring standard deviation of U_1.
    ## sample.size: integer, representing output sample size i.e., N.
    
    # Generate matrix of covariates, along which treatment intensity varies.
    X <- 5 + rnorm(sample.size)
    # Put in terms of a linear model, with \beta_1 - \beta_0 > 0 educ gains.
    # Y_i(Z_i, 0) = \beta_0 \vec X_i + \gamma_0 Z_i + U_0
    # Y_i(Z_i, 1) = \beta_1 \vec X_i + \gamma_1 Z_i + U_0
    beta_0 <- 1
    beta_1 <- 2
    gamma_0 <- 0.5
    gamma_1 <- 1
    mu0_Z0_X <- beta_0 * X
    mu1_Z0_X <- beta_1 * X
    mu0_Z1_X <- mu0_Z0_X + gamma_0
    mu1_Z1_X <- mu1_Z0_X + gamma_1
    muC_Z0_X <- rep(5, sample.size)
    muC_Z1_X <- 0.75 * muC_Z0_X
    # Unobserved selection, between D = 1 and D = 0
    # Assume: U_1, U_0 ~ BivarNormal(rho, mu_0 = mu_1 = 0, sigma_0, sigma_1)
    # For simplicity, independent costs with var = 0.1
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
    #! Test, add another instrument for education/treatment uptake cost.
    probCost <- 0.5
    instCost <- rbinom(sample.size, 1, probCost)
    U_C <- U_C + (instCost - 1) * sigma_C
    # Compose potential outcomes from the observed + unobserved factors.
    # Y_i(Z, D) = mu_D(Z; X_i) + U_D
    Y_0_0 <- mu0_Z0_X + U_0
    Y_0_1 <- mu1_Z0_X + U_1
    Y_1_0 <- mu0_Z1_X + U_0
    Y_1_1 <- mu1_Z1_X + U_1
    # Model compliance as Roy selection --- gains as latent index/prop score.
    # D(Z) = 1{ Y_i(Z, D = 0) <= Y(Z, D_i = 1) - costs }
    #      = 1{ mu_0(Z; X_i) + U_0 <= mu_1(Z; X_i) + U_1 - mu_C(Z; X_i) }
    D_0 <- D_1 <- replicate(0, sample.size)
    D_0 <- as.integer(Y_0_0 <= Y_0_1 - muC_Z0_X - U_C)
    D_1 <- as.integer(Y_1_0 <= Y_1_1 - muC_Z1_X - U_C)
    # Label compliance
    firststage.designation <- replicate(sample.size, "", "vector")
    firststage.designation[D_1 == 0 & D_0 == 0] <- "Never-taker"
    firststage.designation[D_1 == 1 & D_0 == 1] <- "Always-taker"
    firststage.designation[D_1 == 1 & D_0 == 0] <- "Complier"
    firststage.designation[D_1 == 0 & D_0 == 1] <- "Defier"
    # Show how many ATs, NTs, Compliers in terms of D_i(Z) for Z = 0, 1.
    print("How many compliers in the sample?")
    print(table(firststage.designation) / length(firststage.designation))
    # Generate the list of a random instrument.
    probZ <- 0.5
    Z <- rbinom(sample.size, 1, probZ)
    # Generate the list of observed binary treatment
    D <- (Z * D_1) + ((1 - Z) * D_0)
    # Generate the list of observed outcomes
    Y <- (Z * D * Y_1_1) +
        (Z * (1 - D) * Y_1_0) +
        ((1 - Z) * D * Y_0_1) +
        ((1 - Z) * (1 - D) * Y_0_0)
    # Put these data to a coherent data frame.
    combined.data <- tibble(data.frame(
        # Observed data
        Z, D, Y, X,
        # Unobserved, potential outcomes and compliance.
        D_0, D_1,
        Y_0_0, Y_0_1, Y_1_0, Y_1_1,
        mu0_Z0_X, mu1_Z0_X, mu0_Z1_X, mu1_Z1_X, muC_Z0_X, muC_Z1_X, 
        U_0, U_1, U_C, instCost,
        firststage.designation))
    # Return the simulated data as a data frame.
    return(combined.data)
}


################################################################################
## Define a function to show the theoretical values for the data.
theoretical.values <- function(simulated.data, digits.no = 3){
    ### Inputs:
    ## simulated.data, a data frame simulated from above.
    # Extract the potentials from simulated data.
    Z <- simulated.data$Z
    D <- simulated.data$D
    Y <- simulated.data$Y
    X <- simulated.data$X
    D_0 <- simulated.data$D_0
    D_1 <- simulated.data$D_1
    Y_0_0 <- simulated.data$Y_0_0
    Y_0_1 <- simulated.data$Y_0_1
    Y_1_0 <- simulated.data$Y_1_0
    Y_1_1 <- simulated.data$Y_1_1
    firststage.designation <- simulated.data$firststage.designation
    mu0_Z0_X <- simulated.data$mu0_Z0_X
    mu1_Z0_X <- simulated.data$mu1_Z0_X
    mu0_Z1_X <- simulated.data$mu0_Z1_X
    mu1_Z1_X <- simulated.data$mu1_Z1_X
    muC_Z0_X <- simulated.data$muC_Z0_X
    muC_Z1_X <- simulated.data$muC_Z1_X
    U_0 <- simulated.data$U_0
    U_1 <- simulated.data$U_1
    U_C <- simulated.data$U_C
    # Get the true first-stage effects
    first_stage <- D_1 - D_0
    average_first_stage <- mean(first_stage)
    # Get the theoretical total effect/reduced form/ITT
    total_effect <- (Y_1_1 - Y_0_0) * (firststage.designation == "Complier") +
        (Y_1_1 - Y_0_1) * (firststage.designation == "Always-taker") +
        (Y_1_0 - Y_0_0) * (firststage.designation == "Never-taker")
    average_total_effect <- mean(total_effect)
    # Get the theoretical indirect effect, and complier returns to education 
    indirect_effect <- (Z * (Y_1_1 - Y_1_0) + (1 - Z) * (Y_0_1 - Y_0_0)) * (
            firststage.designation == "Complier") +
        (0) * (firststage.designation == "Always-taker") +
        (0) * (firststage.designation == "Never-taker")
    average_indirect_effect <- mean(indirect_effect)
    late <- mean(indirect_effect[firststage.designation == "Complier"])
    # Get the theoretical direct effect.
    direct_effect <- (Z * (Y_1_1 - Y_0_1) + (1 - Z) * (Y_1_0 - Y_0_0)) * (
            firststage.designation == "Complier") +
        (Y_1_1 - Y_0_1) * (firststage.designation == "Always-taker") +
        (Y_1_0 - Y_0_0) * (firststage.designation == "Never-taker")
    average_direct_effect <- mean(direct_effect)
    # Variance in compliance
    var_D <- var(D)
    var_D1_D0 <- var(D_1 - D_0)
    ## True value of population compliance explained by X_i
    # 1. Notice that we can rewrite D_i(1) as
    # D_1 = 1{mu0_Z1_X + U_0 <= mu1_Z1_X + U_1} = 1{U_0-U_1 <= mu1_Z1_X - mu0_Z1_X}
    # Since (U_0 - U_1) is normally distributed, then the conditional mean given X
    prob_D1_X <- pnorm((mu1_Z1_X - mu0_Z1_X - muC_Z1_X) / sd(U_0 - U_1 - U_C))
    # 2. By similar logic for D_i(0)
    # D_0 = 1{mu0_Z0_X + U_0 <= mu1_Z0_X + U_1} = 1{U_0-U_1 <= mu1_Z0_X - mu0_Z0_X}
    prob_D0_X <- pnorm((mu1_Z0_X - mu0_Z0_X - muC_Z0_X) / sd(U_0 - U_1 - U_C))
    cate_D_x <- prob_D1_X - prob_D0_X
    var_cate_D_x <- var(cate_D_x)
    ## Compare to the true value of population mediation. explained by X_i
    D_X <- Z * (mu1_Z1_X - mu0_Z1_X - muC_Z1_X) + (1 - Z) * (mu1_Z0_X - mu0_Z0_X)
    prob_D_X <- pnorm(D_X / sd(U_0 - U_1 - U_C))
    var_prob_D_X <- var(prob_D_X)
    # Define a named list to return
    output.list <- list(
        average_first_stage     = average_first_stage,
        average_total_effect    = average_total_effect,
        average_indirect_effect = average_indirect_effect,
        average_direct_effect   = average_direct_effect,
        complier_indirect_effect= late,
        var_D                   = var_D,
        var_prob_D_X            = var_prob_D_X,
        var_D1_D0               = var_D1_D0,
        var_cate_D_x            = var_cate_D_x)
    # Return the output.list
    return(output.list)
}


################################################################################
## Define a function to estimate mediation, assuming sequential ig. only on X_i.
estimated.values <- function(simulated.data, digits.no = 3){
    ### Inputs:
    ## simulated.data, a data frame simulated from above.
    # Define the first stage / mediation model
    first_stage.reg <- lm(D ~ 1 + Z + X, data = simulated.data)
    # Total effect
    total_stage.reg <- lm(Y ~ 1 + Z + X, data = simulated.data)
    # Define the second stage
    second_stage.reg <- lm(Y ~ 1 + Z + D + Z:D + X, data = simulated.data)
    # Estimate the mechanism model.
    mediation.reg <- mediate(first_stage.reg, second_stage.reg,
        treat = "Z", mediator = "D",
        robustSE = FALSE, sims = 2)
    # Get the point estimates of the direct + indirect effects.
    est_first_stage     <- as.numeric(coef(first_stage.reg )["Z"])
    est_total_effect    <- as.numeric(coef(total_stage.reg)["Z"])
    est_indirect_effect <- mediation.reg$d.avg
    est_direct_effect   <- mediation.reg$z.avg
    # Define a named list to return
    output.list <- list(
        est_first_stage     = est_first_stage,
        est_total_effect    = est_total_effect,
        est_indirect_effect = est_indirect_effect,
        est_direct_effect   = est_direct_effect)
    # Return the output.list
    return(output.list)
}


################################################################################
## Simulate data, for an example.

# Simulate the data, with rho, sigma_0, sigma_1, sigma_C = 3/4, 1, 2, 1
example.data <- simulate.data(3/4, 1, 2, 1)

# Show the theoretical values
print("True treatment effects in the simulation:")
print(theoretical.values(example.data)[
    c("average_first_stage",
        "average_total_effect",
        "average_indirect_effect",
        "average_direct_effect")])

# Plot the Potential Outcomes
outcomes.plot <- example.data %>%
    ggplot(aes(x = X)) +
    geom_point(aes(y = Y_1_1, colour = "Y(1, 1)")) +
    geom_point(aes(y = Y_1_0, colour = "Y(1, 0)")) +
    geom_point(aes(y = Y_0_1, colour = "Y(0, 1)")) +
    geom_point(aes(y = Y_0_0, colour = "Y(0, 0)")) +
    geom_point(aes(y = Y,     colour = "Y"), colour = "black", alpha = 0.25) +
    theme_bw() +
    scale_x_continuous(name = TeX("$X_i$"),
        expand = c(0, 0)) +
    scale_y_continuous(name = "",
        expand = c(0, 0),
        limits = c(0, max(example.data$Y) + 1)) +
    scale_colour_tableau() +
    ggtitle(TeX("$\\Y_i(Z, D_i(Z))$")) +
    theme(plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0.25, 0), "mm"),
        legend.title = element_blank(),
        legend.position = "bottom")
# Save this plot
ggsave(file.path(output.folder, "outcomes-plot.png"),
    plot = outcomes.plot, dpi = 300,
    units = "cm", width = fig.width, height = fig.height)


################################################################################
## Test a control function estimator on the simulated data.
# Package for classical selection estimators (i.e., MLE)
library(sampleSelection)
# Generalised Random Forests, https://grf-labs.github.io/grf/
library(grf)
# DML estimates of the direct and indirect effects.
library(causalweight)

# Simulate the data, with rho, sigma_0, sigma_1, sigma_C = 3/4, 1, 2, 1
example.data <- simulate.data(3/4, 1, 2, 1.5)
Z <- example.data$Z
D <- example.data$D
Y <- example.data$Y
X <- example.data$X
D_0 <- example.data$D_0
D_1 <- example.data$D_1
Y_0_0 <- example.data$Y_0_0
Y_0_1 <- example.data$Y_0_1
Y_1_0 <- example.data$Y_1_0
Y_1_1 <- example.data$Y_1_1
firststage.designation <- example.data$firststage.designation
mu0_Z0_X <- example.data$mu0_Z0_X
mu1_Z0_X <- example.data$mu1_Z0_X
mu0_Z1_X <- example.data$mu0_Z1_X
mu1_Z1_X <- example.data$mu1_Z1_X
muC_Z0_X <- example.data$muC_Z0_X
muC_Z1_X <- example.data$muC_Z1_X
U_0 <- example.data$U_0
U_1 <- example.data$U_1
U_C <- example.data$U_C
instCost <- example.data$instCost

# Show the theoretical values
print("True treatment effects in the simulation:")
print(theoretical.values(example.data)[
    c("average_first_stage",
        "average_total_effect",
        "average_indirect_effect",
        "average_direct_effect")])

# Perfectly get the reduced form equation, using unobserved control functions.
summary(lm(Y ~ (1 + Z * D) * X))
summary(lm(Y ~ (1 + Z * D) * X + U_0 + I(D * (U_1 - U_0))))
summary(lm(Y ~ (1 + Z * D) * X + D : U_1 + I(1 - D) : U_0))
# SImilarly in mediation model, though U_0, U_1 unobserved.
mediation.reg <- mediate(
    lm(D ~ (1 + Z) * poly(X, 3)),
    lm(Y ~ (1 + Z * D) * poly(X, 3) + U_0 + I(D * (U_1 - U_0))),
    treat = "Z", mediator = "D",
    robustSE = FALSE, sims = 300)
print(summary(mediation.reg))
theoretical.values(example.data)[
    c("average_first_stage",
        "average_total_effect",
        "average_indirect_effect",
        "average_direct_effect")]

## Calculate the direct + indirect effect point estimates by expectations.
# Get the total effect by prediction
totaleffect.reg <- lm(Y ~ (1 + Z) * X, data = example.data)
mean(predict(totaleffect.reg, newdata = mutate(example.data, Z = 1)) -
    predict(totaleffect.reg, newdata = mutate(example.data, Z = 0)))
# Get the first-stage by prediction
firststage.reg <- probit(D ~ (1 + Z) * X, data = example.data)
mean(predict(firststage.reg, newdata = mutate(example.data, Z = 1), type = "response") -
    predict(firststage.reg, newdata = mutate(example.data, Z = 0), type = "response"))
# Get the second-stage by prediction
secondstage.reg <- lm(Y ~ (1 + Z * D) * X, data = example.data)
# Indirect effect, E[ Y_i(Z, 1) - Y_i(Z, 0) \times (D_i(1) - D_i(0)) ]
mean((predict(firststage.reg, newdata = mutate(example.data, Z = 1), type = "response") -
    predict(firststage.reg, newdata = mutate(example.data, Z = 0), type = "response")) * (
        predict(secondstage.reg, newdata = mutate(example.data, D = 1)) -
        predict(secondstage.reg, newdata = mutate(example.data, D = 0))))
# Direct Effect.
mean(predict(secondstage.reg, newdata = mutate(example.data, Z = 1)) -
    predict(secondstage.reg, newdata = mutate(example.data, Z = 0)))


#! Unfinished: the cost instrument switching the analysis to cost instrument compliers,
#!             Changing the av effects to LATES w.r.t. random cost variation.
#!             and not average effects.
#TODO: implement the scond instrument analysis in full, identifying those complier effects.

# Note: above estimates can be corrected by adding the control function.
# Unobserved part to add on, + U_0 + I(D * (U_1 - U_0)), see H-H (2024).
# Or estimated with a Heckman (switching) selection approach.
firststage.eq    <- (D ~ (1 + Z * instCost) * poly(X, 3))
outcomeY_Z_0.eq  <- (Y ~ (1 + Z) * poly(X, 3))
outcomeY_Z_1.eq  <- (Y ~ (1 + Z) * poly(X, 3))
selection.reg <- selection(
    firststage.eq,
    list(outcomeY_Z_0.eq, outcomeY_Z_1.eq),
    data = example.data,
    method = "ml")
# Show real estimands.
print(theoretical.values(example.data)[
    c("average_first_stage",
        "average_total_effect",
        "average_indirect_effect",
        "average_direct_effect")])
# Direct Effect, E[ Y_i(1, D) - Y_i(0, D)].
print(mean(predict(selection.reg, newdata = mutate(example.data, Z = 1)) -
    predict(selection.reg, newdata = mutate(example.data, Z = 0))))
# Indirect effect, E[ Y_i(Z, 1) - Y_i(Z, 0) \times (D_i(1) - D_i(0)) ]
firststage.reg <- probit(firststage.eq, data = example.data)
print(mean((
    predict(firststage.reg, newdata = mutate(example.data, Z = 1), type = "response") -
    predict(firststage.reg, newdata = mutate(example.data, Z = 0), type = "response")) * (
        predict(selection.reg, newdata = mutate(example.data, D = 1))[, 2] -
        predict(selection.reg, newdata = mutate(example.data, D = 0))[, 1])))

# Get the CF by hand, assuming Bivar Normal
# https://github.com/dssg/learning/blob/master/Causality/demo-of-iv-and-control-function-methods.r
errors <- predict(firststage.reg)
InvMills_D1 <- ifelse(D == 1,  dnorm(-errors) / (1-pnorm(-errors)), 0)
InvMills_D0 <- ifelse(D == 0, -dnorm(-errors) / (  pnorm(-errors)), 0)
example.data$controlfun  <- ifelse(D == 1, InvMills_D1, InvMills_D0)
example.data$errors  <- predict(firststage.reg, type = "response")
secondstage.reg <- lm(Y ~ (1 + Z * D + controlfun) * (X),
    data = example.data)
print(theoretical.values(example.data)[
    c("average_first_stage",
        "average_total_effect",
        "average_indirect_effect",
        "average_direct_effect")])
# Indirect effect, E[ Y_i(Z, 1) - Y_i(Z, 0) \times (D_i(1) - D_i(0)) ]
mean((predict(firststage.reg, newdata = mutate(example.data, Z = 1), type = "response") -
    predict(firststage.reg, newdata = mutate(example.data, Z = 0), type = "response")) * (
        predict(secondstage.reg, newdata = mutate(example.data, D = 1)) -
        predict(secondstage.reg, newdata = mutate(example.data, D = 0))))
# Direct Effect.
mean(predict(secondstage.reg, newdata = mutate(example.data, Z = 1)) -
    predict(secondstage.reg, newdata = mutate(example.data, Z = 0)))

#! Test: Try it with non-linear control fun, with cost instrument.
cf_firststage.reg <- probit(D ~ (1 + Z) * (1 + poly(X, 3) + poly(U_C, 3)),
    data = example.data)
example.data$controlfun <- (
    example.data$D - predict(cf_firststage.reg, type = "response"))
secondstage.reg <- lm(Y ~ ((1 + Z * D)) * poly(X, 4) + poly(controlfun, 4),
    data = example.data)
summary(secondstage.reg)
print(theoretical.values(example.data)[
    c("average_first_stage",
        "average_total_effect",
        "average_indirect_effect",
        "average_direct_effect")])
# Indirect effect, E[ Y_i(Z, 1) - Y_i(Z, 0) \times (D_i(1) - D_i(0)) ]
mean((predict(firststage.reg, newdata = mutate(example.data, Z = 1), type = "response") -
    predict(firststage.reg, newdata = mutate(example.data, Z = 0), type = "response")) * (
        predict(secondstage.reg, newdata = mutate(example.data, D = 1)) -
        predict(secondstage.reg, newdata = mutate(example.data, D = 0))))
# Direct Effect.
mean(predict(secondstage.reg, newdata = mutate(example.data, Z = 1)) -
    predict(secondstage.reg, newdata = mutate(example.data, Z = 0)))


## Estimate the mediation with a semi-parametric, CF approach.
# CF needs an external instrument to avoid distributional assumptions.
controlfun.forest <- probability_forest(
    X = cbind(example.data$Z, example.data$X, example.data$instCost),
    Y = as.factor(example.data$D))
D.est <- predict(firststage.forest)$predictions[, 2]
controlfun.est <- example.data$D - D.est

#! Open question: does the Abadie Imbens debiased matching estimator work
#!                well for the mediation + control function estimation?
#!                what about a DML approach, where CF is only used in 2nd-stage?
#!                Could one use an external instrument to est the CF, and then
#!                test the seq ig assumption (with CF) using Huber's method?
#!                In a similar sense, the new DR selection model could be used 
#!                in a CF setting.

# Define the first stage / mediation model
first_stage.reg <- lm(D ~ (1 + Z) * poly(X, 4))
# Define the second stage
second_stage.reg <- lm(Y ~ (1 + Z * D) * poly(X, 4) + poly(controlfun.est, 4))
# Estimate the mechanism model.
mediation.reg <- mediate(first_stage.reg, second_stage.reg,
    treat = "Z", mediator = "D",
    robustSE = FALSE, sims = 100)
print(theoretical.values(example.data)[
    c("average_first_stage",
        "average_total_effect",
        "average_indirect_effect",
        "average_direct_effect")])
print(summary(mediation.reg))

#! Comment: the 2nd stage OLS gets it exactly correct.
#!          However, mediation package does not bootstrap for Control function
#!          calculation in the first-stage -> needs to do so as part of my code.
#!          Note: this does not work well only useing non-linear forests,
#!          presumably because of poor overlap.


################################################################################
## Sequence of new IV estimator to try for this 

# 1. Generate a binary instrument for D, say tuition cost instrument, \tilde D.
# 2. Simulate the system with Roy selection Z -> D -> Y, and a direct effect
# 3. Calculate the Abadie (2003) weights for \tilde D -> D
# 4. Run the CM analysis with these weights, and see if the results correctly
#    estimate the direct + indirect effects local to \tilde D -> D compliers.

# Hand-code Weighted least Squares solution (Abadie 2003 Section 4.2.1)
hat_probZ <- glm(Z ~ 1 + poly(X, 3), family = binomial(probit))$fitted
k_1 <- D * ((Z - hat_probZ) / ((1 - hat_probZ) * hat_probZ))
k_0 <- (1 - D) * (((1 - Z) - (1 - hat_probZ)) / ((1 - hat_probZ) * hat_probZ))
k <- k_1 * hat_probZ + k_0 * (1 - hat_probZ)
# Weighted least squares with (possibly) negative weights k
solve(t(cbind(D, 1, X) * k) %*% cbind(D, 1, X)) %*% t(cbind(D, 1, X) * k) %*% as.matrix(Y)


# library(mgcv)
# firstStep <- gam(e401k ~ s(inc) + s(age) + s(agesq) + marr + s(fsize), 
# data=c401k, family=binomial(link = "probit"))
# zProb <- firstStep$fitted
# est2<- larf(nettfa ~ inc + age + agesq + marr + fsize, treatment = p401k, 
# instrument = e401k, data = c401k, zProb = zProb)





################################################################################
## Simulate these data across different parameter values.

#TODO: Make this simulation go across fewer parameters, 
#TODO: but average the point estimate across 100 different simulated data files for each parameter.

# Define a list of rho, sigma_0, sigma_1 values to loop across.
sim.length    <- 10^2
rho.max       <- 0.95
rho.list      <- seq(from = -rho.max, to = rho.max, length.out = 10)
sigma_0.value <- 2
sigma_1.list  <- c(0, 0.25, 0.5, 1, 2, 4)
sigma_C.value <- sigma_0.value
# Generate empty lists to go across.
rho_value.list <- c()
sigma_0_value.list <- c()
sigma_1_value.list <- c()
average_total_effect.list <- c()
average_indirect_effect.list <- c()
average_direct_effect.list <- c()
est_first_stage.list <- c()
est_total_effect.list <- c()
est_indirect_effect.list <- c()
est_direct_effect.list <- c()
var_D.list <- c()
var_prob_D_X.list <- c()
var_D1_D0.list <- c()
var_cate_D_x.list <- c()
# Perform the simulation across these values.
for (rho.value in rho.list){
    for (sigma_1.value in sigma_1.list){
        # Simulate the data, with given rho, sigma_0, sigma_1
        print(paste(rho.value, sigma_0.value, sigma_1.value))
        sim.data <- simulate.data(
            rho.value, sigma_0.value, sigma_1.value, sigma_C.value,
                sample.size = 10^4)
        # Generate the real values, and estimated values.
        theoretical.list <- theoretical.values(sim.data)
        estimated.list   <- estimated.values(sim.data)
        # Store the values.
        rho_value.list               <- c(rho.value,                                rho_value.list)
        sigma_0_value.list           <- c(sigma_0.value,                            sigma_0_value.list)
        sigma_1_value.list           <- c(sigma_1.value,                            sigma_1_value.list)
        sigma_C_value.list           <- c(sigma_C.value,                            sigma_1_value.list)
        average_total_effect.list    <- c(theoretical.list$average_total_effect,    average_total_effect.list)
        est_total_effect.list        <- c(estimated.list$est_total_effect,          est_total_effect.list)
        average_indirect_effect.list <- c(theoretical.list$average_indirect_effect, average_indirect_effect.list)
        est_indirect_effect.list     <- c(estimated.list$est_indirect_effect,       est_indirect_effect.list)
        average_direct_effect.list   <- c(theoretical.list$average_direct_effect,   average_direct_effect.list)
        est_direct_effect.list       <- c(estimated.list$est_direct_effect,         est_direct_effect.list)
        var_D.list                   <- c(theoretical.list$var_D,                   var_D.list)
        var_prob_D_X.list            <- c(theoretical.list$var_prob_D_X,            var_prob_D_X.list)
        var_D1_D0.list               <- c(theoretical.list$var_D1_D0,               var_D1_D0.list)
        var_cate_D_x.list            <- c(theoretical.list$var_cate_D_x,            var_cate_D_x.list)
    }
}


# Put everything together into data/
parameterised.data <- data.frame(
    rho                     = rho_value.list,
    sigma_0                 = sigma_0_value.list,
    sigma_1                 = sigma_1_value.list,
    average_total_effect    = average_total_effect.list,
    est_total_effect        = est_total_effect.list,
    average_indirect_effect = average_indirect_effect.list,
    est_indirect_effect     = est_indirect_effect.list,
    average_direct_effect   = average_direct_effect.list,
    est_direct_effect       = est_direct_effect.list,
    var_D                   = var_D.list,
    var_prob_D_X            = var_prob_D_X.list,
    var_D1_D0               = var_D1_D0.list,
    var_cate_D_x            = var_cate_D_x.list)


################################################################################
## Plot the bias in the estimates, using linear CM methods (Imai+ 2010).

# Plot the bias in reduced form effect est vs rho, with different sigma values.
reducedform_bias.plot <- parameterised.data %>%
    ggplot(aes(x = rho,
        y = ((est_total_effect - average_total_effect
            ) / average_total_effect),
        colour = factor(sigma_1 / sigma_0))) +
    geom_jitter(alpha = 2 / 3) +
    geom_smooth(aes(group = sigma_1 / sigma_0),
        se = FALSE, method = "loess", formula = y ~ x) +
    theme_bw() +
    scale_x_continuous(
        name = TeX("$\\rho$"),
        expand = c(0, 0),
        breaks = seq(-1, 1, by = 0.2),
        limits = c(-1, 1)) +
    scale_y_continuous(name = "",
        #breaks = seq(0, 0.1, by = 0.02),
        #limits = c(0, 0.1),
        expand = c(0, 0)) +
    ggtitle(TeX(paste0(
        "Percent Bias, given $sigma_{0} = ",
            sigma_0.value, "$"))) +
    guides(colour = guide_legend(
        TeX("$\\sigma_1 / \\sigma_0$"), nrow = 1)) +
    theme(plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0.25, 0), "mm"),
        legend.position = "bottom")
# Save this plot
ggsave(file.path(output.folder, "reducedform-bias.png"),
    plot = reducedform_bias.plot, dpi = 300,
    units = "cm", width = fig.width, height = fig.height)

# Plot the bias in direct effect est vs rho, with different sigma values
directeffect_bias.plot <- parameterised.data %>%
    ggplot(aes(x = rho,
        y = ((est_direct_effect - average_direct_effect
            ) / average_direct_effect),
        colour = factor(sigma_1 / sigma_0))) +
    geom_jitter(alpha = 2 / 3) +
    geom_smooth(aes(group = sigma_1 / sigma_0),
        se = FALSE, method = "loess", formula = y ~ x) +
    theme_bw() +
    scale_x_continuous(
        name = TeX("$\\rho$"),
        expand = c(0, 0),
        breaks = seq(-1, 1, by = 0.2),
        limits = c(-1, 1)) +
    scale_y_continuous(name = "",
        #breaks = seq(-1, 1, by = 0.1),
        #limits = c(-1, 1),
        expand = c(0, 0)) +
    ggtitle(TeX(paste0(
        "Percent Bias, given $sigma_0 = ",
            sigma_0.value, "$"))) +
    guides(colour = guide_legend(
        TeX("$\\sigma_1 / \\sigma_0$"), nrow = 1)) +
    theme(plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0.25, 0.25), "mm"),
        legend.position = "bottom")
# Save this plot
ggsave(file.path(output.folder, "directeffect-bias.png"),
    plot = directeffect_bias.plot, dpi = 300,
    units = "cm", width = fig.width, height = fig.height)


# Plot the bias in indirect effect est vs rho, with different sigma values
indirecteffect_bias.plot <- parameterised.data %>%
    ggplot(aes(x = rho,
        y = ((est_indirect_effect - average_indirect_effect
            ) / average_indirect_effect),
        colour = factor(sigma_1 / sigma_0))) +
    geom_jitter(alpha = 0.5) +
    geom_smooth(aes(group = sigma_1 / sigma_0),
        se = FALSE, method = "loess", formula = y ~ x) +
    theme_bw() +
    scale_x_continuous(
        expand = c(0, 0),
        name = TeX("$\\rho$"),
        breaks = seq(-1, 1, by = 0.2),
        limits = c(-1, 1)) +
    scale_y_continuous(name = "",
        #breaks = seq(-1, 1, by = 0.1),
        #limits = c(-1, 1),
        expand = c(0, 0)) +
    ggtitle(TeX(paste0(
        "Percent Bias, given $sigma_0 = ",
            sigma_0.value, "$"))) +
    guides(colour = guide_legend(
        TeX("$\\sigma_1 / \\sigma_0$"), nrow = 1)) +
    theme(plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0.25, 0.25), "mm"),
        legend.position = "bottom")
# Save this plot
ggsave(file.path(output.folder, "indirecteffect-bias.png"),
    plot = indirecteffect_bias.plot, dpi = 300,
    units = "cm", width = fig.width, height = fig.height)


################################################################################
## Plot the real $R^2_U$ values

# $R^2_U$, where $R^2_X$ is based on $Var(D_i(1) - D_i(0))$
R2_U_var_D1_D0.plot <- parameterised.data %>%
    ggplot(aes(x = rho, y = var_cate_D_x / var_D1_D0,
        colour = factor(sigma_1 / sigma_0))) +
    geom_jitter(alpha = 0.5) +
    geom_smooth(aes(group = sigma_1 / sigma_0),
        se = FALSE, method = "loess", formula = y ~ x) +
    theme_bw() +
    scale_x_continuous(
        expand = c(0, 0),
        name = TeX("$\\rho$"),
        breaks = seq(-1, 1, by = 0.2),
        limits = c(-1, 1)) +
    scale_y_continuous(name = "",
        #breaks = seq(0, 0.5, by = 0.025),
        #limits = c(0, 0.2),
        expand = c(0, 0)) +
    ggtitle(TeX(paste0(
        "$R^2_U$, given $sigma_0 = ",
            sigma_0.value, "$"))) +
    guides(colour = guide_legend(
        TeX("$\\sigma_1 / \\sigma_0$"), nrow = 1)) +
    theme(plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0.25, 0), "mm"),
        legend.position = "bottom")
# Save this plot
ggsave(file.path(output.folder, "R2Uvar-D1D0.png"),
    plot = R2_U_var_D1_D0.plot, dpi = 300,
    units = "cm", width = fig.width, height = fig.height)

# $R^2_U$, where $R^2_X$ is based on $Var(D_i)$
R2_U_var_D.plot <- parameterised.data %>%
    ggplot(aes(x = rho, y = var_prob_D_X / var_D,
        colour = factor(sigma_1 / sigma_0))) +
    geom_jitter(alpha = 0.5) +
    geom_smooth(aes(group = sigma_1 / sigma_0),
        se = FALSE, method = "loess", formula = y ~ x) +
    theme_bw() +
    scale_x_continuous(name = TeX("$\\rho$"),
        breaks = seq(-1, 1, by = 0.2),
        limits = c(-1, 1),
        expand = c(0, 0)) +
    scale_y_continuous(name = "",
        #breaks = seq(0, 0.5, by = 0.025),
        #limits = c(0, 0.2),
        expand = c(0, 0)) +
    ggtitle(TeX(paste0(
        "$R^2_U$, given $sigma_0 = ",
            sigma_0.value, "$"))) +
    guides(colour = guide_legend(
        TeX("$\\sigma_1 / \\sigma_0$"), nrow = 1)) +
    theme(plot.title = element_text(hjust = 0, size = rel(1)),
        plot.title.position = "plot",
        plot.margin = unit(c(0.5, 3, 0.25, 0), "mm"),
        legend.position = "bottom")
# Save this plot
ggsave(file.path(output.folder, "R2Uvar-D.png"),
    plot = R2_U_var_D.plot, dpi = 300,
    units = "cm", width = fig.width, height = fig.height)
