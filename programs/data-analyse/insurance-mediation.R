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


################################################################################ Define the functions to use.

# Estimate the values, given a first and second-stages
estimated.values <- function(firststage.reg, secondstage.reg, totaleffect.reg,
    input.data, complier.adjustment = NULL){
    ### Inputs:
    ## input.data, a data frame simulated from above.
    input_Z0.data <- input.data
    input_Z1.data <- input.data
    input_D0.data <- input.data
    input_D1.data <- input.data
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

# Define a function to Heckman selection correct mediation est, in two-stages.
mediate.naive <- function(Y, Z, D, X_iv, X_minus, input.data){
    # Make the names consistent
    input.data$Y <- input.data[[Y]]
    input.data$Z <- input.data[[Z]]
    input.data$D <- input.data[[D]]
    input.data$X_iv <- input.data[[X_iv]]
    input.data$X_minus <- input.data[[X_minus]]
    # 0. Total effect regression.
    totaleffect.reg <- lm(Y ~ 1 + Z + X_minus,
        data = input.data)
    # 1. Probit first-stage (well identified).
    firststage.reg <- glm(D ~ (1 + Z) + X_iv + X_minus,
        family = binomial(link = "binomial"),
        data = input.data)
    # 2. Estimate second-stage, including the CFs.
    naive_secondstage.reg <- lm(Y ~ (1 + Z * D) + X_minus,
        data = input.data)
    # Compile the estimates.
    naive.est <- estimated.values(
        firststage.reg, naive_secondstage.reg,
        totaleffect.reg, input.data,
        complier.adjustment = NULL)
    # Return the off-setting estimates.
    output.list <- c(
        naive.est$"first-stage",
        naive.est$"total-effect",
        naive.est$"direct-effect",
        naive.est$"indirect-effect")
    return(output.list)
}

# Define a function to Heckman selection correct mediation est, in two-stages.
mediate.heckit <- function(Y, Z, D, X_iv, X_minus, input.data){
    # Make the names consistent
    input.data$Y <- input.data[[Y]]
    input.data$Z <- input.data[[Z]]
    input.data$D <- input.data[[D]]
    input.data$X_iv <- input.data[[X_iv]]
    input.data$X_minus <- input.data[[X_minus]]
    # Get relevant columns for imputation.
    input_Z0.data <- input.data
    input_Z1.data <- input.data
    input_D0.data <- input.data
    input_D1.data <- input.data
    input_Z0.data$Z <- 0
    input_Z1.data$Z <- 1
    input_D0.data$D <- 0
    input_D1.data$D <- 1
    # 0. Total effect regression.
    totaleffect.reg <- lm(Y ~ 1 + Z + X_minus,
        data = input.data)
    # 1. Probit first-stage (well identified).
    heckit_firststage.reg <- glm(D ~ (1 + Z) + X_iv + X_minus,
        family = binomial(link = "probit"),
        data = input.data)
    # 2. Define the CFs --- for assumed N(0,1) dist.
    lambda_1.fun <- function(pi.est){
        # Inv Mills ratio, taking as input the estimated mediator propensity.
        return(dnorm(qnorm(pi.est)) / pnorm(qnorm(pi.est)))
    }
    pi.est <- predict(heckit_firststage.reg, type = "response")
    input.data$lambda_0 <- (1 - input.data$D) * lambda_1.fun(pi.est) * (
        - pi.est / (1 - pi.est))
    input.data$lambda_1 <- input.data$D * lambda_1.fun(pi.est)
    # 3. Estimate second-stage, including the CFs.
    heckit_secondstage.reg <- lm(Y ~ (1 + Z * D) + X_minus + lambda_0 + lambda_1,
        data = input.data)
    # Compensate complier difference in AIE, by Kline Walters (2019) IV-type adjustment.
    pi_0.est <- predict(heckit_firststage.reg,
        newdata = input_Z0.data, type = "response")
    pi_1.est <- predict(heckit_firststage.reg,
        newdata = input_Z1.data, type = "response")
    Gamma.big <-  (pi_1.est * lambda_1.fun(pi_1.est)
        - pi_0.est * lambda_1.fun(pi_0.est)) / (pi_1.est - pi_0.est)
    rho_0 <- coef(summary(heckit_secondstage.reg))["lambda_0", "Estimate"]
    rho_1 <- coef(summary(heckit_secondstage.reg))["lambda_1", "Estimate"]
    complier.adjustment <- (rho_1 - rho_0) * Gamma.big
    # Compile the estimates.
    heckit.est <- estimated.values(
        heckit_firststage.reg, heckit_secondstage.reg,
        totaleffect.reg, input.data,
        complier.adjustment = complier.adjustment)
    # Return the off-setting estimates.
    output.list <- c(
        heckit.est$"first-stage",
        heckit.est$"total-effect",
        heckit.est$"direct-effect",
        heckit.est$"indirect-effect")
    return(output.list)
}

# Define a function to two-stage semi-parametric CF for CM effects.
mediate.semiparametric <- function(Y, Z, D, X_iv, X_minus, input.data){
    # Make the names consistent
    input.data$Y <- input.data[[Y]]
    input.data$Z <- input.data[[Z]]
    input.data$D <- input.data[[D]]
    input.data$X_iv <- input.data[[X_iv]]
    input.data$X_minus <- input.data[[X_minus]]
    # Get relevant columns for imputation.
    input_Z0.data <- input.data
    input_Z1.data <- input.data
    input_D0.data <- input.data
    input_D1.data <- input.data
    input_Z0.data$Z <- 0
    input_Z1.data$Z <- 1
    input_D0.data$D <- 0
    input_D1.data$D <- 1
    # 1. Total effect regression.
    totaleffect.reg <- lm(Y ~ 1 + Z + X_minus,
        data = input.data)
    totaleffect.est <- mean(predict(
        totaleffect.reg, newdata = input_Z1.data) - predict(
            totaleffect.reg, newdata = input_Z0.data))
    # 2. Semi-parametric first-stage
    cf_firststage.reg <- semiBRM(D ~ 1 + Z + X_iv + X_minus,
        data = input.data, control = list(iterlim = 100))
    input.data$pi.est <- predict(cf_firststage.reg, type = "response")$prob
    pi_0.est <- predict(cf_firststage.reg,
        newdata = input_Z0.data, type = "response")$prob
    pi_1.est <- predict(cf_firststage.reg,
        newdata = input_Z1.data, type = "response")$prob
    pi.bar <- mean(pi_1.est - pi_0.est)
    # 3. Semi-parametric series estimation of the second-stage.
    cf_secondstage_D0.reg <- gam(Y ~ 1 + Z + X_minus + s(pi.est, bs = "cr"),
        method = "REML",
        data = input.data, subset = (D == 0))
    cf_secondstage_D1.reg <- gam(Y ~ 1 + Z + X_minus + s(pi.est, bs = "cr"),
        method = "REML",
        data = input.data, subset = (D == 1))
    # 4. Compose the CM effects from this object.
    D_0 <- 1 - mean(input.data$D)
    D_1 <- 1 - D_0
    # 4.1 ADE point estimate, from the CF model.
    gammma.est <- coef(cf_secondstage_D0.reg)["Z"]
    delta_plus.est <- coef(cf_secondstage_D1.reg)["Z"]
    ade.est <- as.numeric(D_0 * gammma.est + D_1 * delta_plus.est)
    # 4.2 AIE by using ADE estimate, relative to ATE.
    # (Avoiding semi-parametric extrapolation, see notes on ATE comparison)
    delta.est <- delta_plus.est - gammma.est
    ade_Z0.est <- gammma.est + delta.est * mean(
        input.data$D[input.data$Z == 0])
    ade_Z1.est <- gammma.est + delta.est * mean(
        input.data$D[input.data$Z == 1])
    aie.est <- (totaleffect.est - mean(
        (1 - input.data$Z) * ade_Z1.est + (input.data$Z) * ade_Z0.est))
    # Return the estimates.
    output.list <- c(
        pi.bar,
        totaleffect.est,
        ade.est,
        aie.est)
    return(output.list)
}

##TODO Define a function to use these functions,
##TODO and provide a summary object (with booted SEs).

mediate.selection <- function(Y, Z, D, X_iv, X_minus, input.data,
    type = "parametric", boot.reps = NULL){
    # If no bootstrapping, just give the point estimates.
    if (is.null(boot.reps)){
        if (type == "parametric"){
            point.est <- mediate.heckit(Y, Z, D, X_iv, X_minus, input.data)
        }
        else if (type == "semi-parametric"){
            point.est <- mediate.semiparametric(
                Y, Z, D, X_iv, X_minus, input.data)
        }
        else if (type == "naive"){
            point.est <- mediate.naive(Y, Z, D, X_iv, X_minus, input.data)
        }
        else {
            stop(paste0("The type option only takes values of ",
                'c("parametric", "semi-parametric", "naive").'))}
        return(point.est)
    }
    # IF bootstrapping for SEs, then calculated via the bootstrap.

}


## Test on the Oregon data.
analysis.data$intercept <- 0
mediate.est <- mediate.selection(
    Y = "Y_happy", Z = "any_insurance", D = "any_healthcare",
    X_iv = "usual_health_location", X_minus = "intercept"
    #TODO: Z_iv = "lottery_iv", 
    #TODO: control_iv = "hh_size"
)
    
    Y, Z, D, X_iv, X_minus, input.data,


# Show that the health location influences healthcare take-up
lm(any_healthcare ~ 0 + factor(usual_health_location), data = analysis.data) %>%
    summary()
# Note values in usual health location:
# 1 private clinic
# 2 public clinic
# 3 hospital-based clinic
# 4 hospital ER
# 5 urgent care clinic
# 6 other place
# 7 don't have usual place
#! TEST: is it related to the IV?  Answer: not significantly.
lm(lottery_iv ~ 0 + factor(usual_health_location), data = analysis.data) %>%
    summary()






quit("no")
# Code from LARF package for estimating OLS with (possibly negative) k weights.
X <- 1
solve ( t(cbind(Z,X) * kappa) %*% cbind(Z,X)) %*% t(cbind(Z,X) * kappa)  %*% Y_health
