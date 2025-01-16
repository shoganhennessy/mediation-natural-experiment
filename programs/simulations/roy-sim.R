#!/usr/bin/R
## Senan Hogan-Hennessy, 16 Jan October 2025
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
# Package for classical selection estimators (i.e., MLE)
library(sampleSelection)

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
sample.size <- 10^4

################################################################################
## First Simulation: triangular system, all assumptions good.

## Simulate the covariates (observed + unobserved).
# First covariate (\vec X_i^-)
X_minus <- rnorm(sample.size, mean = 4, sd = 1)
# Second covariate (instrument for the control fun).
X_IV <- rbinom(sample.size, 1, 1 / 2)
# Simulate the unobserved error terms
rho <- 3 / 4
sigma_0 <- 1
sigma_1 <- 2 * sigma_0
sigma_C <- 1
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
# Define the potential outcomes.
mu_outcome <- function(d, z, x_minus){
    return(x_minus + (z + d + z * d))
}
mu_cost <- function(z, x_minus, x_iv){
    return(- 3 * z + x_minus * x_iv)
}
# Y_i(Z, D) = mu_D(Z; X_i) + U_D
Y_0_0 <- mu_outcome(0, 0, X_minus) + U_0
Y_0_1 <- mu_outcome(0, 1, X_minus) + U_1
Y_1_0 <- mu_outcome(1, 0, X_minus) + U_0
Y_1_1 <- mu_outcome(1, 1, X_minus) + U_1
# D_i(Z)= 1{ Y(Z, 1) - Y(Z, 0) >= C_i }
D_0 <- as.integer(Y_0_1 - Y_0_0 >= mu_cost(0, X_minus, X_IV) + U_C)
D_1 <- as.integer(Y_1_1 - Y_1_0 >= mu_cost(1, X_minus, X_IV) + U_C)
# Generate the individual effects (direct + indirect)
probZ <- 0.5
Z <- rbinom(sample.size, 1, probZ)
direct_effects <-
    (Z * (Y_1_1 - Y_0_1) + (1 - Z) * (Y_1_0 - Y_0_0)) * (D_1 == 1 & D_0 == 0) +
    (Y_1_1 - Y_0_1) * (D_1 == 1 & D_0 == 1) +
    (Y_1_0 - Y_0_0) * (D_1 == 0 & D_0 == 0)
indirect_effects <-
    (Z * (Y_1_1 - Y_1_0) + (1 - Z) * (Y_0_1 - Y_0_0)) * (D_1 != D_0)
# Observed outcomes: D, Y
D <- (Z * D_1) + ((1 - Z) * D_0)
# Generate the list of observed outcomes
Y <- (Z * D * Y_1_1) +
    (Z * (1 - D) * Y_1_0) +
    ((1 - Z) * D * Y_0_1) +
    ((1 - Z) * (1 - D) * Y_0_0)
# Save the dataframe
example.data <- data.frame(Z, D, Y, X_minus, X_IV)

## Show how the system operates:
# How many mediator D_i compliers (w.r.t. Z)?  No defiers.
print(table(D_1, D_0) / sample.size)
# Show that the regression specification holds exactly (with correlated error).
# True result: Y = 0 + 1 D + 1 Z + 1 Z D + X_minus + (1 - D) U_0 - D U_1
# Error term is a function of D, so ensure tractability by subtracting from Y
print(summary(lm(
    I(Y - (1 - D) * U_0 - D * U_1) ~ 1 + D + Z + Z:D + X_minus)))
# SHow how the OLS result gives a bias result (if rho != 0)
print(summary(lm(Y ~ 1 + D + Z + Z:D + X_minus)))

# SHow the exact same thing, for CM effects
print(paste0(c("ADE: ", mean(direct_effects))))
print(paste0(c("AIE: ", mean(indirect_effects))))
#print(paste0(c("First-stage: ", mean(D_1 - D_0))))

# True model (including unobserved error terms):
firststage.reg <- lm(D ~ (1 + Z) + X_minus * X_IV)
secondstage.reg <-  lm(I(Y - (1 - D) * U_0 - D * U_1) ~ (1 + Z * D) + X_minus)
print(paste0(c("ADE est: ",
    mean(coef(secondstage.reg)["Z"] + D * coef(secondstage.reg)["Z:D"]))))
print(paste0(c("AIE est: ", mean(coef(firststage.reg)["Z"] * (
    coef(secondstage.reg)["D"] + Z * coef(secondstage.reg)["Z:D"])))))

# Standard OLS approach fails
firststage.reg <- lm(D ~ (1 + Z) + X_minus * X_IV, data = example.data)
secondstage.reg <-  lm(Y ~ (1 + Z * D) + X_minus, data = example.data)
print(paste0(c("ADE est: ",
    mean(coef(secondstage.reg)["Z"] + D * coef(secondstage.reg)["Z:D"]))))
print(paste0(c("AIE est: ", mean(coef(firststage.reg)["Z"] * (
    coef(secondstage.reg)["D"] + Z * coef(secondstage.reg)["Z:D"])))))


# CF approach, with a probit
controlfun.reg <- probit(D ~ (1 + Z) * X_minus * X_IV, data = example.data)
example.data$K <- example.data$D - predict(controlfun.reg, type = "response")
# controlfun.forest <- grf::probability_forest(X = cbind(Z, X_minus, X_IV), Y = as.factor(D))
# example.data$K <- D - predict(controlfun.forest)$predictions[, 2]
firststage.reg <- lm(D ~ (1 + Z) * poly(X_minus, 5) * X_IV, data = example.data)
secondstage.reg <-  lm(Y ~ (1 + Z * D) + X_minus + poly(K, 5), data = example.data)
mean(predict(secondstage.reg, newdata = mutate(example.data, Z = 1)) -
    predict(secondstage.reg, newdata = mutate(example.data, Z = 0)))
print(paste0(c("ADE truth: ", mean(direct_effects))))
mean((predict(firststage.reg, newdata = mutate(example.data, Z = 1), type = "response") -
    predict(firststage.reg, newdata = mutate(example.data, Z = 0), type = "response")) * (
        predict(secondstage.reg, newdata = mutate(example.data, D = 1)) -
        predict(secondstage.reg, newdata = mutate(example.data, D = 0))))
print(paste0(c("AIE truth: ", mean(indirect_effects))))







controlfun.forest <- grf::probability_forest(X = cbind(Z, X_minus, X_IV), Y = as.factor(D))
example.data$K <- D - predict(controlfun.forest)$predictions[, 2]
print(summary(lm(Y ~ (1 + Z * D) + X_minus + poly(K, 3))))

K_0 <- (1 - D) * predict(controlfun.forest)$predictions[, 1]
K_1 <- D * predict(controlfun.forest)$predictions[, 2]
print(summary(lm(Y ~ (1 + Z * D) + X_minus + poly(K_0, 3) + poly(K_1, 3))))


# K <- D * K + (1 - D) * (1 - K)
K <- D - K
print(summary(lm(Y ~ (1 + Z * D) + poly(X_minus, 5) + poly(K, 5))))



firststage.reg <- lm(D ~ (1 + Z) + X_minus * X_IV)
secondstage.reg <- lm(Y ~ (1 + Z * D) + X_minus + poly(K, 3))
print(paste0(c("ADE est: ",
    mean(coef(secondstage.reg)["Z"] + D * coef(secondstage.reg)["Z:D"]),
    "truth: ", mean(direct_effects))))
print(paste0(c("AIE est: ", mean(coef(firststage.reg)["Z"] * (
    coef(secondstage.reg)["D"] + Z * coef(secondstage.reg)["Z:D"])),
    "truth: ", mean(indirect_effects))))

print(summary(lm(D ~ (1 + Z) + X_minus * X_IV)))
print(mean(D_1 - D_0))
