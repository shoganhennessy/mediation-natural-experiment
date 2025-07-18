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
# Package forsemi-parametric CF by splines.
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
presentation.height <- (2 / 3) * presentation.width
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
analysis.data$usual_health_location <- factor(
    analysis.data$usual_health_location)


################################################################################
## Define the functions to use.

################################################################################
## Weighted-LS *or* Probit with (possibly) negative weight
################################################################################
wls_negweight <- function(formula, w_input, data,
    family = c("gaussian", "probit"),
    control = list(maxit = 100, reltol = 1e-8)) {
        family <- match.arg(family)
        ## model frame / matrices
        mf <- model.frame(formula, data, na.action = na.omit)
        y  <- model.response(mf)
        X  <- model.matrix(formula, mf)
        w  <- data[[w_input]]
        if (family == "gaussian") {
        ############# OLS / WLS ####################################################
        Xw   <- X * w                # apply weight
        XtWX <- t(Xw) %*% X
        XtWY <- t(Xw) %*% y
        beta <- as.vector(solve(XtWX, XtWY))

        fitted <- as.vector(X %*% beta)
        resid  <- y - fitted
        df     <- nrow(X) - ncol(X)
        sigma2 <- sum(resid^2) / df
        vcovB  <- sigma2 * solve(XtWX)

  } else {
    ############# Probit MLE with weight #####################################
    if (!all(y %in% 0:1))
      stop("For probit you need a 0/1 outcome.")

    nll <- function(beta, X, y, w, eps = 1e-12) {           # –log-lik
      eta <- X %*% beta
      p   <- pnorm(eta)
      # keep inside (0,1) to avoid log(0)
      p   <- pmin(pmax(p, eps), 1 - eps)
      -sum(w * (y * log(p) + (1 - y) * log(1 - p)))
    }

    ## start values: unweighted probit via glm (allows only +ve weight)
    start  <- stats::glm.fit(X, y, family = binomial(link = probit))$coefficients
    opt    <- optim(par   = start,
                    fn    = nll,
                    X = X, y = y, w = w,
                    method = "BFGS",
                    control = list(maxit = control$maxit,
                                   reltol = control$reltol),
                    hessian = TRUE)

    if (opt$convergence != 0)
      warning("optim() did not fully converge")

    beta   <- as.vector(opt$par)
    Hinv   <- tryCatch(solve(opt$hessian),
                       error = function(e) { matrix(NA, ncol(X), ncol(X)) })
    vcovB  <- Hinv
    fitted <- as.vector(pnorm(X %*% beta))
    resid  <- y - fitted
    df     <- nrow(X) - ncol(X)
    sigma2 <- NA                                # not defined for probit
  }

  #### assemble return object ##################################################
  out <- list(coefficients   = beta,
              residuals      = resid,
              fitted.values  = fitted,
              vcov           = vcovB,
              sigma          = sqrt(sigma2),
              df.residual    = df,
              formula        = formula,
              call           = match.call(),
              terms          = attr(mf, "terms"),
              family         = family)

  class(out) <- "wls_negweight"
  out
}

###############################################################################
## predict method
###############################################################################
predict.wls_negweight <- function(object, newdata = NULL,
                                   type = c("response", "link"), ...) {

  type <- match.arg(type)

  ## 1.  Compute the linear predictor  η  (“link” scale)
  if (is.null(newdata)) {
    # from the fitted object itself
    if (object$family == "gaussian") {
      eta <- as.vector(object$fitted.values)      # for OLS link==response
    } else {
      eta <- as.vector(qnorm(object$fitted.values)) # inverse link for probit
    }
  } else {
    # for new data
    Xnew <- model.matrix(delete.response(object$terms), newdata)
    eta  <- as.vector(Xnew %*% object$coefficients)
  }

  ## 2.  Return either link or response scale
  if (object$family == "gaussian") {
    return(eta)                    # same for "response" and "link"
  } else {                         # probit
    if (type == "link") return(eta)
    return(pnorm(eta))             # response = probability
  }
}

###############################################################################
## summary method (t- or z-tests as appropriate)
###############################################################################
summary.wls_negweight <- function(object, ...) {
  se  <- sqrt(diag(object$vcov))
  tst <- object$coefficients / se
  if (object$family == "gaussian") {
    pval <- 2 * pt(-abs(tst), df = object$df.residual)
    stat <- "t value"
  } else {                         # probit uses z statistics
    pval <- 2 * pnorm(-abs(tst))
    stat <- "z value"
  }

  tab <- cbind(Estimate = object$coefficients,
               `Std. Error` = se,
               setNames(tst, stat),
               `Pr(>|t|)` = pval)

  cat("\nCall:\n")
  print(object$call)
  cat("\nCoefficients:\n")
  print(tab, digits = 4)
  if (object$family == "gaussian")
    cat("\nResidual SE:", round(object$sigma, 4),
        "on", object$df.residual, "DF\n")
  invisible(tab)
}

###############################################################################
## print method
###############################################################################
print.wls_negweight <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients, ...)
  invisible(x)
}

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

# Define a function to Heckman selection correct mediation est, in two-stages.
mediate.unadjusted <- function(Y, Z, D, X_iv, X_minus, data,
    Z_iv = NULL, control_iv = NULL, indices = NULL){
    # Bootstrap sample, if indices provided.
    if (!is.null(indices)){
        data <- data[indices, ]
    }
    # Make the names consistent
    data$Y <- data[[Y]]
    data$Z <- data[[Z]]
    data$D <- data[[D]]
    data$X_iv <- data[[X_iv]]
    data$X_minus <- data[[X_minus]]
    # Start counting coefficients estimated.
    count.coef <- 0
    # If given a treatment IV, use that.
    if (!is.null(Z_iv)){
        data$Z_iv <- data[[Z_iv]]
        if (!is.null(control_iv)){
            data$control_iv <- data[[control_iv]]}
        else { data$control_iv <- 0}
        # Calculate the kappa.weight weight
        est_probZ_iv <- glm(Z_iv ~ 1 + control_iv, data = data)
        hat_probZ_iv <- est_probZ_iv$fitted
        kappa.weight_0 <- (1 - data$Z) * (((1 - data$Z_iv
            ) - (1 - hat_probZ_iv)) / ( 
                (1 - hat_probZ_iv) * hat_probZ_iv))
        kappa.weight_1 <- data$Z * ((data$Z_iv - hat_probZ_iv) / (
            (1 - hat_probZ_iv) * hat_probZ_iv))
        data$kappa.weight <- kappa.weight_0 * (
            1 - hat_probZ_iv) + kappa.weight_1 * hat_probZ_iv
        count.coef <- count.coef + length(est_probZ_iv$coefficients)
    }
    # If no treatment IV, set the IV re-weighting to unity.
    else { data$kappa.weight <- 1 }
    # 0. Total effect regression.
    totaleffect.reg <- wls_negweight(Y ~ 1 + Z + X_minus,
        "kappa.weight", family = "gaussian",
        data = data)
    count.coef <- count.coef + length(totaleffect.reg$coefficients)
    # 1. Regular first-stage (well identified).
    firststage.reg <- wls_negweight(D ~ 1 + Z * X_iv + X_minus,
        "kappa.weight", family = "probit",
        data = data)
    count.coef <- count.coef + length(firststage.reg$coefficients)
    # 2. Estimate second-stage (naive case has no CFs).
    unadjusted_secondstage.reg <- wls_negweight(Y ~ (1 + Z * D) + X_minus,
        "kappa.weight", family = "gaussian",
        data = data)
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
    Y = "Y_health", Z = "any_insurance", D = "any_healthcare",
    X_iv = "usual_health_location", X_minus = "hh_size",
    Z_iv = "lottery_iv", control_iv = "hh_size",
    analysis.data)
print(mediate.est)


# Define a function to Heckman selection correct mediation est, in two-stages.
lambda_1.fun <- function(pi.est){
        # Inv Mills ratio, taking as input the estimated mediator propensity.
        return(dnorm(qnorm(pi.est)) / pnorm(qnorm(pi.est)))
    }
# THe actual function.
mediate.heckit <- function(Y, Z, D, X_iv, X_minus, data,
    Z_iv = NULL, control_iv = NULL, indices = NULL){
    # Bootstrap sample, if indices provided.
    if (!is.null(indices)){
        data <- data[indices, ]
    }
    # Make the names consistent
    data$Y <- data[[Y]]
    data$Z <- data[[Z]]
    data$D <- data[[D]]
    data$X_iv <- data[[X_iv]]
    data$X_minus <- data[[X_minus]]
    # Start counting coefficients estimated.
    count.coef <- 0
    # If given a treatment IV, use that.
    if (!is.null(Z_iv)){
        data$Z_iv <- data[[Z_iv]]
        if (!is.null(control_iv)){
            data$control_iv <- data[[control_iv]] }
        else { data$control_iv <- 0}
        # Calculate the kappa.weight weight
        est_probZ_iv <- glm(Z_iv ~ 1 + control_iv, data = data)
        hat_probZ_iv <- est_probZ_iv$fitted
        kappa.weight_0 <- (1 - data$Z) * (((1 - data$Z_iv
            ) - (1 - hat_probZ_iv)) / ( 
                (1 - hat_probZ_iv) * hat_probZ_iv))
        kappa.weight_1 <- data$Z * ((data$Z_iv - hat_probZ_iv) / (
            (1 - hat_probZ_iv) * hat_probZ_iv))
        data$kappa.weight <- kappa.weight_0 * (
            1 - hat_probZ_iv) + kappa.weight_1 * hat_probZ_iv
        count.coef <- count.coef + length(est_probZ_iv$coefficients)
    }
    # If no treatment IV, set the IV re-weighting to unity.
    else { data$kappa.weight <- 1 }
    # 0. Total effect regression.
    totaleffect.reg <- wls_negweight(Y ~ 1 + Z + X_minus,
        "kappa.weight", family = "gaussian",
        data = data)
    count.coef <- count.coef + length(totaleffect.reg$coefficients)
    # 1. Probit first-stage (well identified).
    heckit_firststage.reg <- wls_negweight(D ~ 1 + Z + X_iv + X_minus,
        "kappa.weight", family = "probit",
        data = data)
    count.coef <- count.coef + length(heckit_firststage.reg$coefficients)
    # 2. Define the CFs --- for assumed N(0,1) dist.
    pi.est <- predict(heckit_firststage.reg)
    data$lambda_0 <- (1 - data$D) * lambda_1.fun(pi.est) * (
        - pi.est / (1 - pi.est))
    data$lambda_1 <- data$D * lambda_1.fun(pi.est)
    # 3. Estimate second-stage, including the CFs.
    heckit_secondstage.reg <- wls_negweight(Y ~ 1 + Z + D + Z:D + X_minus +
        lambda_0 + lambda_1,
        "kappa.weight", family = "gaussian",
        data = data)
    count.coef <- count.coef + length(heckit_secondstage.reg$coefficients)
    # Compensate complier difference in AIE, by Kline Walters (2019) IV-type adjustment.
    input_Z0.data <- data
    input_Z1.data <- data
    input_Z0.data$Z <- 0
    input_Z1.data$Z <- 1
    pi_0.est <- predict(heckit_firststage.reg,
        newdata = input_Z0.data)
    pi_1.est <- predict(heckit_firststage.reg,
        newdata = input_Z1.data)
    Gamma.big <-  (pi_1.est * lambda_1.fun(pi_1.est)
        - pi_0.est * lambda_1.fun(pi_0.est)) / (pi_1.est - pi_0.est)
    rho_0 <- coef(heckit_secondstage.reg)[6] #as.numeric(coef(heckit_secondstage.reg)["lambda_0"])
    rho_1 <- coef(heckit_secondstage.reg)[7] #as.numeric(coef(heckit_secondstage.reg)["lambda_1"])
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
    Y = "Y_health", Z = "any_insurance", D = "any_healthcare",
    X_iv = "usual_health_location", X_minus = "hh_size",
    Z_iv = "lottery_iv", control_iv = "hh_size",
    analysis.data)
print(mediate.est)


# Define a function to two-stage semi-parametric CF for CM effects.
mediate.semiparametric <- function(Y, Z, D, X_iv, X_minus, data,
    Z_iv = NULL, control_iv = NULL, indices = NULL){
    # Bootstrap sample, if indices provided.
    if (!is.null(indices)){
        data <- data[indices, ]
    }
    # Make the names consistent
    data$Y <- data[[Y]]
    data$Z <- data[[Z]]
    data$D <- data[[D]]
    data$X_iv <- data[[X_iv]]
    data$X_minus <- data[[X_minus]]
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
    # If given a treatment IV, use that.
    if (!is.null(Z_iv)){
        data$Z_iv <- data[[Z_iv]]
        if (!is.null(control_iv)){
            data$control_iv <- data[[control_iv]]}
        else { data$control_iv <- 0}
        # Calculate the kappa.weight weights
        est_probZ_iv <- glm(Z_iv ~ 1 + control_iv, data = data)
        hat_probZ_iv <- est_probZ_iv$fitted
        kappa.weight_0 <- (1 - data$Z) * (((1 - data$Z_iv
            ) - (1 - hat_probZ_iv)) / ( 
                (1 - hat_probZ_iv) * hat_probZ_iv))
        kappa.weight_1 <- data$Z * ((data$Z_iv - hat_probZ_iv) / (
            (1 - hat_probZ_iv) * hat_probZ_iv))
        data$kappa.weight <- kappa.weight_0 * (
            1 - hat_probZ_iv) + kappa.weight_1 * hat_probZ_iv
        count.coef <- count.coef + length(est_probZ_iv$coefficients)
    }
    # If no treatment IV, set the IV re-weighting to unity.
    else { data$kappa.weight <- 1 }
    # 1. Total effect regression.
    totaleffect.reg <- wls_negweight(Y ~ 1 + Z + X_minus,
        "kappa.weight", family = "gaussian",
        data = data)
    totaleffect.est <- mean(predict(
        totaleffect.reg, newdata = input_Z1.data) - predict(
            totaleffect.reg, newdata = input_Z0.data))
    count.coef <- count.coef + length(totaleffect.reg$coefficients)
    # 2. Semi-parametric first-stage
    cf_firststage.reg <- wls_negweight(D ~ 1 + Z * X_iv + X_minus,
        "kappa.weight", family = "probit",
        data = data)
    data$pi.est <- predict(cf_firststage.reg)
    pi_0.est <- predict(cf_firststage.reg, newdata = input_Z0.data)
    pi_1.est <- predict(cf_firststage.reg, newdata = input_Z1.data)
    pi.bar <- mean(pi_1.est - pi_0.est)
    count.coef <- count.coef + length(cf_firststage.reg$coefficients)
    # Calculate the levels of pi, accounting for few values in the IV.
    distinct_cf.values <- min(
        length(unique(data$pi.est)) - 2, as.integer(nrow(data) / 2000))
    # 3. Semi-parametric series estimation of the second-stage.
    cf_secondstage_D0.reg <- wls_negweight(Y ~ 1 + Z + X_minus  +
        bs(pi.est, df = distinct_cf.values),
        "kappa.weight", family = "gaussian",
        data = data[data$D == 0, ])
    cf_secondstage_D1.reg <- wls_negweight(Y ~ 1 + Z + X_minus +
        bs(pi.est, df = distinct_cf.values),
        "kappa.weight", family = "gaussian",
        data = data[data$D == 1, ])
    count.coef <- count.coef + length(cf_secondstage_D0.reg$coefficients)
    count.coef <- count.coef + length(cf_secondstage_D1.reg$coefficients)
    # 4. Compose the CM effects from this object.
    D_0 <- 1 - mean(data$D)
    D_1 <- 1 - D_0
    # 4.1 ADE point estimate, from the CF model.
    gammma.est <- coef(cf_secondstage_D0.reg)[2]#["Z"]
    delta_plus.est <- coef(cf_secondstage_D1.reg)[2]#["Z"]
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
    Y = "Y_health", Z = "any_insurance", D = "any_healthcare",
    X_iv = "usual_health_location", X_minus = "hh_size",
    Z_iv = "lottery_iv", control_iv = "hh_size",
    analysis.data)
print(mediate.est)

# Define a function to bootstrap.
mediate.bootstrap <- function(Y, Z, D, X_iv, X_minus, data,
        Z_iv = NULL, control_iv = NULL,
        type = "parametric", boot.reps = 10){
    # Define an empty data.frame.
    boot.data <- data.frame(matrix(ncol = 4, nrow = 0))
    names(boot.data) <- c(
        "First-stage", "ATE, Total", "ADE", "AIE")
    j <- 1
    for (i in 1:boot.reps){
        print(i)
        boot.indices <- sample(1:nrow(data), nrow(data), replace = TRUE)
        tryCatch({
            error
            if (type == "parametric"){
                point.est <- mediate.heckit(Y, Z, D, X_iv, X_minus, data,
                    Z_iv = Z_iv, control_iv = control_iv, indices = boot.indices)
                }
                else if (type == "semi-parametric"){
                    point.est <- mediate.semiparametric(
                        Y, Z, D, X_iv, X_minus, data,
                        Z_iv = Z_iv, control_iv = control_iv, indices = boot.indices)
                }
                else if (type == "unadjusted"){
                    point.est <- mediate.unadjusted(Y, Z, D, X_iv, X_minus, data,
                        Z_iv = Z_iv, control_iv = control_iv, indices = boot.indices)
                }
                else {
                    stop(paste0("The type option only takes values of ",
                        'c("parametric", "semi-parametric", "unadjusted").'))
                }
            boot.data[j, ] <- point.est 
            j <- j + 1
        }, error = function(msg){
            print(paste0("Missed ", j))
        })
    }
    return(boot.data)
}

#! Test it out.
mediate.boot <- mediate.bootstrap(
    Y = "Y_health", Z = "any_insurance", D = "any_healthcare",
    X_iv = "usual_health_location", X_minus = "hh_size",
    Z_iv = "lottery_iv", control_iv = "hh_size",
    data = analysis.data, type = "parametric", boot.reps = 10)
print(mediate.boot)


## Define a function to wrap around all the others.
mediate.selection <- function(Y, Z, D, X_iv, X_minus, data,
    Z_iv = NULL, control_iv = NULL,
    type = "parametric", boot.reps = 10){
    # Calculate the point estimates.
    if (type == "parametric"){
        point.est <- mediate.heckit(Y, Z, D, X_iv, X_minus, data,
            Z_iv = Z_iv, control_iv = control_iv)
    }
    else if (type == "semi-parametric"){
        point.est <- mediate.semiparametric(
            Y, Z, D, X_iv, X_minus, data,
                Z_iv = Z_iv, control_iv = control_iv)
    }
    else if (type == "unadjusted"){
        point.est <- mediate.unadjusted(Y, Z, D, X_iv, X_minus, data,
            Z_iv = Z_iv, control_iv = control_iv)
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
        point.boot <- mediate.bootstrap(Y, Z, D, X_iv, X_minus, data,
            Z_iv = Z_iv, control_iv = control_iv,
            type = type, boot.reps = boot.reps)
    }
    # Report output
    point.est <- as.matrix(c(point.est[1:4], point.est[4] / point.est[2]))
    point.se <- as.matrix(c(
        sd(point.boot$"First-stage"),
        sd(point.boot$"ATE, Total"),
        sd(point.boot$"ADE"),
        sd(point.boot$"AIE"),
        sd(point.boot$"AIE" / point.boot$"ATE, Total")))
    tratio <- as.matrix(point.est / point.se)
    ptratio <- as.matrix(2 * pt(abs(tratio),
        df = nrow(data) - count.coef, lower.tail = FALSE))
    # Preapred object to putput.
    out <- list(
        coefficients = point.est,
        SE = point.se,
        tratio = tratio,
        ptratio = ptratio,
        variables = paste(Z, D, Y, sep = ", "),
        boot.reps = boot.reps)
    rownames(out$coefficients) <-
        c("First-stage", "ATE, Total", "ADE", "AIE", "Proportion, AIE / ATE")
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
    cat(paste0("\nEstimates, With SEs from ", x$boot.reps, " bootstrap replications."))
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

# Show that the health location influences healthcare take-up
location.data <- analysis.data %>%
    group_by(usual_health_location) %>%
    summarise(
        any_healthcare_mean = mean(any_healthcare, na.rm = TRUE),
        any_healthcare_sd = sd(any_healthcare, na.rm = TRUE),
        count = n()) %>%
    ungroup() %>%
    mutate(any_healthcare_se = any_healthcare_sd / (count^(0.5)))
# Note values in usual health location:
# 1 private clinic
# 2 public clinic
# 3 hospital-based clinic
# 4 hospital ER
# 5 urgent care clinic
# 6 other place
# 7 don't have usual place


#TODO: code this as a simple a figure, justifying the instrument.
#TODO: panel (a) first-stage, (b) second-stage effect on health + happiness.


################################################################################
## Estimate the CM effects with my methods.

# State how many bootstrap replications are needed.
boot.reps <- 10^3

## Panel A: Self-reported healthiness.
# Naive selection-on-observables
unadjusted.est <- mediate.selection(
    Y = "Y_health", Z = "any_insurance", D = "any_healthcare",
    X_iv = "usual_health_location", X_minus = "hh_size",
    Z_iv = "lottery_iv", control_iv = "hh_size",
    type = "unadjusted",
    data = analysis.data,
    boot.reps = boot.reps)
print(unadjusted.est)
# Parametric CF
parametric.est <- mediate.selection(
    Y = "Y_health", Z = "any_insurance", D = "any_healthcare",
    X_iv = "usual_health_location", X_minus = "hh_size",
    Z_iv = "lottery_iv", control_iv = "hh_size",
    type = "parametric",
    data = analysis.data,
    boot.reps = boot.reps)
print(parametric.est)
# Semi-parametric CF
semiparametric.est <- mediate.selection(
    Y = "Y_health", Z = "any_insurance", D = "any_healthcare",
    X_iv = "usual_health_location", X_minus = "hh_size",
    Z_iv = "lottery_iv", control_iv = "hh_size",
    type = "semi-parametric",
    data = analysis.data,
    boot.reps = boot.reps)
print(semiparametric.est)

# Extract the relevant figures.
panelA.data <- data.frame(
    unadjusted_point = coef(summary(unadjusted.est))[, "Estimate"],
    unadjusted_se = coef(summary(unadjusted.est))[, "SE"],
    parametric_point = coef(summary(parametric.est))[, "Estimate"],
    parametric_se = coef(summary(parametric.est))[, "SE"],
    semiparametric_point = coef(summary(semiparametric.est))[, "Estimate"],
    semiparametric_se = coef(summary(semiparametric.est))[, "SE"])

# Clean up the data.
panelA.data <- panelA.data %>%
    signif(digits.no) %>%
    format(scientific = FALSE) %>%
    as.character()
panelA.data$unadjusted_se <- panelA.data$unadjusted_se %>% paste0("(", ., ")")
panelA.data$parametric_se <- panelA.data$parametric_se %>% paste0("(", ., ")")
panelA.data$semiparametric_se <- panelA.data$semiparametric_se %>% paste0("(", ., ")")
panelA.data <- data.frame(t(panelA.data))
# Add on the first column
panelA.data$model <- c(
    "Unadjusted", "", "Parametric CF", "", "Semi-parametric CF", "")
panelA.data <- panelA.data[c(6, 1:5)]

# Save the LaTeX table
panelA.data %>%
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
    Y = "Y_happy", Z = "any_insurance", D = "any_healthcare",
    X_iv = "usual_health_location", X_minus = "hh_size",
    Z_iv = "lottery_iv", control_iv = "hh_size",
    type = "unadjusted",
    data = analysis.data,
    boot.reps = boot.reps)
# Parametric CF
parametric.est <- mediate.selection(
    Y = "Y_happy", Z = "any_insurance", D = "any_healthcare",
    X_iv = "usual_health_location", X_minus = "hh_size",
    Z_iv = "lottery_iv", control_iv = "hh_size",
    type = "parametric",
    data = analysis.data,
    boot.reps = boot.reps)
# Semi-parametric CF
semiparametric.est <- mediate.selection(
    Y = "Y_happy", Z = "any_insurance", D = "any_healthcare",
    X_iv = "usual_health_location", X_minus = "hh_size",
    Z_iv = "lottery_iv", control_iv = "hh_size",
    type = "semi-parametric",
    data = analysis.data,
    boot.reps = boot.reps)

# Extract the relevant figures.
panelB.data <- data.frame(
    unadjusted_point = coef(summary(unadjusted.est))[, "Estimate"],
    unadjusted_se = coef(summary(unadjusted.est))[, "SE"],
    parametric_point = coef(summary(parametric.est))[, "Estimate"],
    parametric_se = coef(summary(parametric.est))[, "SE"],
    semiparametric_point = coef(summary(semiparametric.est))[, "Estimate"],
    semiparametric_se = coef(summary(semiparametric.est))[, "SE"])

# Clean up the data.
panelB.data <- panelB.data %>%
    signif(digits.no) %>%
    format(scientific = FALSE) %>%
    as.character()
panelB.data$unadjusted_se <- panelB.data$unadjusted_se %>% paste0("(", ., ")")
panelB.data$parametric_se <- panelB.data$parametric_se %>% paste0("(", ., ")")
panelB.data$semiparametric_se <- panelB.data$semiparametric_se %>% paste0("(", ., ")")
panelB.data <- data.frame(t(panelB.data))
# Add on the first column
panelB.data$model <- c(
    "Unadjusted", "", "Parametric CF", "", "Semi-parametric CF", "")
panelB.data <- panelB.data[c(6, 1:5)]

# Save the LaTeX table
panelB.data %>%
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
