## -----------------------------------------------------------------------------
## Description: This script performs estimation of cause-specific absolute risks
##              when using a case-cohort with 
##              - exhaustive sampling of the cases with an event of interest
##              (y1) or with a competing event (y2)
##              - a fixed baseline hazard for y1 and y2
##              - design weights 
##              - calibrated weights  
##
##              Calibration of the design weights is against 
##              - (i)+(ii)+(iii) when focusing on y1 only
##              - (i)+(ii)+(iii) when focusing on y1 and y2
##              - (i)+(ii) when focusing on y1 only
##
##              We compare the robust variance estimate and the variance 
##              estimate that properly accounts for the case-cohort sampling, as
##              proposed in Etievant and Gail (2024)
## -----------------------------------------------------------------------------
## Required Package: dplyr, ggplot2, grid, gridExtra, gtable, nnet, parallel, 
## survival, xtable
## -----------------------------------------------------------------------------
## Required Functions: help.functions.R
## -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Read functions and packages --------------------------------------------------

library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(nnet)
library(parallel)
library(survival)
library(xtable)

source("help.functions.R")

# ------------------------------------------------------------------------------
# Different scenarios (with different combinations of parameters) to be 
# investigated over simulations ------------------------------------------------

parameters  <- NULL 
prob.y1     <- c(0.02, 0.05, 0.1) 
prob.y2     <- c(0.05, 0.1, 0.2)
n           <- c(5000)
noncases.per.case <- c(1, 2)
part        <- 10
for (l in 1:length(prob.y1)) {
  for (k in 1:length(prob.y2)) {
    for (m in 1:length(n)) {
      for (i in 1:length(noncases.per.case)) {
        for (j in 1:part) {
          parameters <- rbind(parameters, c(prob.y1[l], prob.y2[k], n[m], 
                                            noncases.per.case[i], j))
        }}}}}

colnames(parameters) <- c("Prob.y1", "Prob.y2", "n", "noncases.per.case", 
                          "part")

# ------------------------------------------------------------------------------
# Fixed parameters -------------------------------------------------------------

alpha11  <- 0.05 
alpha21  <- -0.35
alpha12  <- -0.25 
alpha22  <- 0.5
x11.cat  <- 1:3
x12.cat  <- 1:3
Prob.x21 <- matrix(c(0.7, 0.05, 0.25, 0.45, 0.2, 0.35, 0.40, 0.3, 0.3), nrow = 3,
                  byrow = TRUE)
Prob.x22 <- matrix(c(0.6, 0.15, 0.25, 0.5, 0.2, 0.3, 0.3, 0.4, 0.3), nrow = 3,
                   byrow = TRUE)
beta11  <- - 0.2 
beta21  <- 0.25
beta31  <- - 0.4
beta12  <- - 0.1 
beta22  <- 0
beta32  <- 0.6
names.beta.y1  <- paste0("beta", c(1,2,3), "1")
names.beta.y2  <- paste0("beta", c(1,3), "2")
time    <- 10 

shape.y1  <- 1
shape.y2  <- 1

t1      <- 0
t2      <- 8
x1      <- c(-1, 1, -0.6) # covariate profiles for the pure and absolute risks
x2      <- c(0.5, 0, 1)   # covariate profiles for the pure and absolute risks

densityX1X3.y1 <- function (x) {
  return(1 / sqrt(2 * pi) * exp(-1/2 * (x - (beta11 + beta31 * alpha11))^2))
}
densityX1X3.y2 <- function (x) {
  return(1 / sqrt(2 * pi) * exp(-1/2 * (x - (beta12 + beta32 * alpha11))^2))
}
densityX1X3.y1 <- Vectorize(densityX1X3.y1)
E.exp.y1       <- exp(1/2 * beta31^2 + 1/2 * (beta11 + beta31 * alpha11)^2) * 
  (integrate(f = densityX1X3.y1, lower = -10, upper = -2)$value * 
     sum(exp((beta21 + beta31 * alpha21) * 0:2) * Prob.x21[1,]) + 
     integrate(densityX1X3.y1, -2, 1)$value * 
     sum(exp((beta21 + beta31 * alpha21) * 0:2) * Prob.x21[2,]) + 
     integrate(densityX1X3.y1, 1, 10)$value * 
     sum(exp((beta21 + beta31 * alpha21) * 0:2) * Prob.x21[3,])) # E(exp(beta11 * X1 + beta21 * X21 + beta31 * X3))
densityX1X3.y2 <- Vectorize(densityX1X3.y2)
E.exp.y2       <- exp(1/2 * beta32^2 + 1/2 * (beta12 + beta32 * alpha12)^2) * 
  (integrate(f = densityX1X3.y2, lower = -10, upper = -1)$value * 
     sum(exp((beta22 + beta32 * alpha22) * 0:2) * Prob.x22[1,]) + 
     integrate(densityX1X3.y2, -1, 2)$value * 
     sum(exp((beta22 + beta32 * alpha22) * 0:2) * Prob.x22[2,]) + 
     integrate(densityX1X3.y2, 2, 10)$value * 
     sum(exp((beta22 + beta32 * alpha22) * 0:2) * Prob.x22[3,])) # E(exp(beta12 * X1 + beta32 * X3))

Nreplic <- 5000 / 10 # we break down the replications in 10

set.seed(12345)
seed <- round(abs(rnorm(10) * 10^5))

Onerun <- function (p) {
  
  Prob.y1           <- parameters[p, 1] # not accounting for censoring and competing events
  Prob.y2           <- parameters[p, 2] # not accounting for censoring and competing events
  n                 <- parameters[p, 3] 
  K                 <- parameters[p, 4] # number of non-cases we wish to sample for each case
  part              <- parameters[p, 5]
  
  set.seed(seed[part])
  
  names(n) <- NULL
  
  # Time-fixed cause-specific baseline hazards ---------------------------------
  bh.y1             <- Prob.y1 / E.exp.y1 / time # (time-fixed) baseline hazard for event 1
  names(bh.y1)      <- NULL

  bh.y2             <- Prob.y2 / E.exp.y2 / time # (time-fixed) baseline hazard for event 2
  names(bh.y2)      <- NULL

  # True values of the pure risk and absolute risk in (t1, t2] and for covariate
  # profiles x -----------------------------------------------------------------
  
  AR.y1.int <- function (t) {
    bh.y1 * exp(x1 %*% c(beta11, beta21, beta31)) * 
      exp(- (bh.y1 * exp(x1 %*% c(beta11, beta21, beta31)) * (t^(shape.y1) - t1^(shape.y1)) + 
               bh.y2 * exp(x2 %*% c(beta12, beta22, beta32)) * (t^(shape.y2) - t1^(shape.y2)))) *
      shape.y1 * t^(shape.y1 - 1)
  }
  AR.y1.int <- Vectorize(AR.y1.int)
  AR.y1 <- integrate(f = AR.y1.int, lower = t1, upper = t2)$value
  
  PR.y1             <- c(1 - exp(- bh.y1 * exp(x1 %*% c(beta11, beta21, beta31)) * t2^(shape.y1)))
  
  res               <- NULL
  
  for (nrep in 1:Nreplic) {
    
    # --------------------------------------------------------------------------
    # Generation of the cohort -------------------------------------------------
    
    # Generation of the covariates ---------------------------------------------
    
    X11      <- rnorm(n)
    X11.cat  <- 1 * (X11 < -2) + 2 *((X11 >= -2)&(X11 < 1)) + 3 * (X11 >= 1)
    X21      <- rep(NA, n)
    for (i in 1:length(x11.cat)) {
      X21[which(X11.cat == x11.cat[i])] <- sample(0:2, 
                                                  sum(X11.cat == x11.cat[i]), 
                                                  replace = TRUE, 
                                                  prob = Prob.x21[i,])
    }
    mu31     <- alpha11 * X11 + alpha21 * X21
    X31      <- mapply(mu31, FUN = X.generation)
    
    X12      <- rnorm(n)
    X12.cat  <- 1 * (X12 < -1) + 2 *((X12 >= -1)&(X12 < 2)) + 3 * (X12 >= 2)
    X22      <- rep(NA, n)
    for (i in 1:length(x12.cat)) {
      X22[which(X12.cat == x12.cat[i])] <- sample(0:2, 
                                                  sum(X12.cat == x12.cat[i]), 
                                                  replace = TRUE, 
                                                  prob = Prob.x22[i,])
    }
    mu32    <- alpha12 * X12 + alpha22 * X22
    X32     <- mapply(mu32, FUN = X.generation)
    
    # Generation of the proxy variables for the covariates ---------------------
    
    X11.proxy  <- X11 + rnorm(n, sd = 0.75) 
    X31.proxy  <- X31 + rnorm(n, sd = 0.75)
    
    X12.proxy  <- X12 + rnorm(n, sd = 0.75) 
    X32.proxy  <- X32 + rnorm(n, sd = 0.75)
 
    
    # Generation of the entry times, event times, and times to independent 
    # censoring ----------------------------------------------------------------
    
    entry.times     <- runif(n, min = 0, max = 5) 

    U               <- runif(n)
    S <- function (j, t) {
      exp(- t^shape.y1 * bh.y1 * exp(beta11 * X11[j] + beta21 * X21[j] + 
                                       beta31 * X31[j]) - 
            t^shape.y2 * bh.y2 * exp(beta12 * X12[j] + beta22 * X21[j] + 
                                       beta32 * X32[j])) - U[j]
    }
    y.times   <- sapply(1:n, function(j) uniroot(S, c(0.0001, 1000), 
                                                 tol = 10^(-5), 
                                                 extendInt ="yes", j = j)$root)
    
    U               <- runif(n)
    death.times     <- (-log(U) / (0.02 / 10)) # 2% pure censoring at 10 years
    
    
    # Observed times and case indicators ---------------------------------------
    
    times           <- pmin(death.times, y.times, 10 - entry.times) 
    status          <- rep(0, n)
    
    prob = ((shape.y1 * y.times^(shape.y1 - 1) * bh.y1 * exp(beta11 * X11 + 
                                                               beta21 * X21 + 
                                                               beta31 * X31)) / 
              (shape.y1 * y.times^(shape.y1 - 1) * bh.y1 * exp(beta11 * X11 + 
                                                                 beta21 * X21 + 
                                                                 beta31 * X31) + 
                 shape.y2 * y.times^(shape.y2 - 1) * bh.y2 * exp(beta12 * X12 + 
                                                                   beta22 * X22 + 
                                                                   beta32 * X32)))
    bin             <- rbinom(n = sum(times == y.times), size = 1, 
                              prob = prob[times == y.times])
    status[times == y.times] <- 2 - (bin == 1)
    status1         <- 1 * (status == 1)
    status2         <- 1 * (status == 2)

    # Number of cases for the two event types ----------------------------------
    n.y1      <- sum(status1 == 1)
    n.y2      <- sum(status2 == 1)
    n.y0      <- n - n.y1 - n.y2
    n.strata  <- c(n.y0, n.y1, n.y2) 
    
    m         <- round(K * n.y1 / (1 - n.y1 / n))

    cohort    <- cbind(id = 1:n, X11, X21, X31, X12, X22, X32, X11.proxy, 
                       X31.proxy, X12.proxy, X32.proxy, times, status, status1, 
                       status2, n, m, n.strata = n.strata[status + 1], 
                       n.y1 = n.y1, n.y2 = n.y2)
    cohort    <- as.data.frame(cohort)
    
    
    # --------------------------------------------------------------------------
    # Sampling of the case cohort ----------------------------------------------
    
    # Fixed number of individuals that will be sampled to obtain the case-cohort
    
    indiv.sampled         <- sample(cohort$id, size = m, replace = FALSE)
    phase2                <- rep(0, n)
    phase2[indiv.sampled] <- 1
    subcohort             <- phase2
    phase2[cohort$status1 == 1] <- 1 # we add all the remaining cases for event type 1
    phase2[cohort$status2 == 1] <- 1 # and for event type 2
    
    cohort$weights <- 1
    cohort$weights[(cohort$status == 0)] <- n / m # 1 / sampling probability
    
    cohort <- cbind(cohort, subcohort, phase2)
    
    casecohort <- cohort
    casecohort <- casecohort[-which((phase2 == 0)), ]
    
    
    # --------------------------------------------------------------------------
    # Estimation using the whole cohort ----------------------------------------
    
    mod.cohort.y1  <- coxph(Surv(times, status1) ~ X11 + X21 + X31, data = cohort, 
                            robust = TRUE) # treating event 2 as independent censoring
    
    mod.cohort.y2  <- coxph(Surv(times, status2) ~ X12 + X32, data = cohort, 
                            robust = TRUE) # treating event 1 as independent censoring
    
    est.y1.cohort <- influences(mod = mod.cohort.y1, x = c(-1, 1, -0.6), t1 = 0, 
                                t2 = 8)
    est.y2.cohort <- influences(mod = mod.cohort.y2, x = c(0.5, 1), t1 = 0, 
                                t2 = 8, competing = TRUE, 
                                event.times.y1 = est.y1.cohort$event.times)
    est.risk.cohort <- risk.influences(est.y1.cohort, est.y2.cohort)
    
    beta.y1.hat.cohort <- mod.cohort.y1$coefficients
    beta.y2.hat.cohort <- mod.cohort.y2$coefficients
    
    AR.y1.hat.cohort <- est.risk.cohort$AR.y1
    PR.y1.hat.cohort <- est.risk.cohort$PR.y1
    
    robvar.beta.y1.hat.cohort   <- robustvariance(est.y1.cohort$infl.beta) # equal to diag(mod.cohort.y1$var)
    robvar.beta.y2.hat.cohort   <- robustvariance(est.y2.cohort$infl.beta) # equal to diag(mod.cohort.y2$var)
    robvar.AR.y1.hat.cohort     <- robustvariance(est.risk.cohort$infl.AR.y1)
    robvar.PR.y1.hat.cohort     <- robustvariance(est.risk.cohort$infl.PR.y1)
    
    est.y2.cohort   <- influences(mod = mod.cohort.y2, x = c(0.5, 1), t1 = 0, 
                                t2 = 8)
    n.events.y2     <- length(est.y2.cohort$Lambda0.t12.hat)
    Lambda0.t12.y2.hat.cohort <- est.y2.cohort$Lambda0.t12.hat[n.events.y2]
    robvar.Lambda0.t12.y2.hat.cohort   <- robustvariance(as.matrix(est.y2.cohort$infl.Lambda0.t12[,n.events.y2]))
    
    res.cohort <- cbind("Cohort", t(beta.y1.hat.cohort), t(beta.y2.hat.cohort),
                        Lambda0.t12.y2.hat.cohort,
                        AR.y1.hat.cohort, PR.y1.hat.cohort, 
                        t(robvar.beta.y1.hat.cohort), 
                        t(robvar.beta.y2.hat.cohort),
                        robvar.Lambda0.t12.y2.hat.cohort,
                        robvar.AR.y1.hat.cohort, robvar.PR.y1.hat.cohort)
    
    
    # --------------------------------------------------------------------------
    # Estimation using the case-cohort with design weights ---------------------
    
    mod.casecohort.y1  <- coxph(Surv(times, status1) ~ X11 + X21 + X31, 
                                data = casecohort, robust = TRUE, 
                                weights = weights)
    
    mod.casecohort.y2  <- coxph(Surv(times, status2) ~ X12 + X32, 
                                data = casecohort,
                                robust = TRUE, weights = weights)
    
    est.y1.casecohort <- influences(mod = mod.casecohort.y1, x = c(-1, 1, -0.6), 
                                    t1 = 0, t2 = 8)
    est.y2.casecohort <- influences(mod = mod.casecohort.y2, x = c(0.5, 1), 
                                    t1 = 0, t2 = 8, competing = TRUE, 
                                    event.times.y1 = est.y1.casecohort$event.times)
    est.risk.casecohort <- risk.influences(est.y1.casecohort, est.y2.casecohort)
    
    beta.y1.hat.casecohort <- mod.casecohort.y1$coefficients
    beta.y2.hat.casecohort <- mod.casecohort.y2$coefficients
     
    AR.y1.hat.casecohort <- est.risk.casecohort$AR.y1
    PR.y1.hat.casecohort <- est.risk.casecohort$PR.y1
    
    robvar.beta.y1.hat.casecohort   <- robustvariance(est.y1.casecohort$infl.beta) # equal to diag(mod.casecohort.y1$var)
    robvar.beta.y2.hat.casecohort   <- robustvariance(est.y2.casecohort$infl.beta) # equal to diag(mod.casecohort.y2$var)
    robvar.AR.y1.hat.casecohort     <- robustvariance(est.risk.casecohort$infl.AR.y1)
    robvar.PR.y1.hat.casecohort     <- robustvariance(est.risk.casecohort$infl.PR.y1)
    
    var.beta.y1.hat.casecohort      <- variance(casecohort, est.y1.casecohort$infl.beta) 
    var.beta.y2.hat.casecohort      <- variance(casecohort, est.y2.casecohort$infl.beta) 
    var.AR.y1.hat.casecohort        <- variance(casecohort, est.risk.casecohort$infl.AR.y1) 
    var.PR.y1.hat.casecohort        <- variance(casecohort, est.risk.casecohort$infl.PR.y1) 
    
    est.y2.casecohort   <- influences(mod = mod.casecohort.y2, x = c(0.5, 1), t1 = 0, 
                                  t2 = 8)
    n.events.y2 <- length(est.y2.casecohort$Lambda0.t12.hat)
    Lambda0.t12.y2.hat.casecohort <- est.y2.casecohort$Lambda0.t12.hat[n.events.y2]
    robvar.Lambda0.t12.y2.hat.casecohort   <- robustvariance(est.y2.casecohort$infl.Lambda0.t12[,n.events.y2]) 
    var.Lambda0.t12.y2.hat.casecohort      <- variance(casecohort = casecohort, 
                                                                  infl = as.matrix(est.y2.casecohort$infl.Lambda0.t12[,n.events.y2])) 

    res.casecohort.rob <- cbind("CaseCohort Robust", t(beta.y1.hat.casecohort), 
                            t(beta.y2.hat.casecohort), 
                            Lambda0.t12.y2.hat.casecohort,
                            AR.y1.hat.casecohort, 
                            PR.y1.hat.casecohort, 
                            t(robvar.beta.y1.hat.casecohort),
                            t(robvar.beta.y2.hat.casecohort),
                            robvar.Lambda0.t12.y2.hat.casecohort,
                            robvar.AR.y1.hat.casecohort, 
                            robvar.PR.y1.hat.casecohort)
    
    res.casecohort <- cbind("CaseCohort", t(beta.y1.hat.casecohort), 
                            t(beta.y2.hat.casecohort), 
                            Lambda0.t12.y2.hat.casecohort,
                            AR.y1.hat.casecohort, 
                            PR.y1.hat.casecohort, t(var.beta.y1.hat.casecohort),
                            t(var.beta.y2.hat.casecohort),
                            var.Lambda0.t12.y2.hat.casecohort,
                            var.AR.y1.hat.casecohort, var.PR.y1.hat.casecohort)
    
    
    # --------------------------------------------------------------------------
    # Estimation using the case-cohort with calibrated weights -----------------
    
    # Imputing the phase-two covariates on the entire cohort -------------------
    
    mod.pred.X31    <- lm(X31 ~ X31.proxy + X11.proxy + X21, data = casecohort, 
                          weights = weights)
    cohort$X31.pred <- predict(mod.pred.X31, 
                               newdata = as.data.frame(cohort %>% 
                                                         select(X31.proxy, 
                                                                X11.proxy, 
                                                                X21)), 
                               type = "response")
    mod.pred.X11    <- lm(X11 ~ X11.proxy, data = casecohort, weights = weights)
    cohort$X11.pred <- predict(mod.pred.X11, 
                               newdata = as.data.frame(cohort %>% 
                                                         select(X11.proxy)), 
                               type = "response")
    
    mod.pred.X32    <- lm(X32 ~ X32.proxy + X12.proxy + X22, data = casecohort, 
                          weights = weights)
    cohort$X32.pred <- predict(mod.pred.X32, 
                               newdata = as.data.frame(cohort %>% 
                                                         select(X32.proxy, 
                                                                X12.proxy, 
                                                                X22)), 
                               type = "response")
    mod.pred.X12    <- lm(X12 ~ X12.proxy, data = casecohort, weights = weights)
    cohort$X12.pred <- predict(mod.pred.X12, 
                               newdata = as.data.frame(cohort %>% 
                                                         select(X12.proxy)), 
                               type = "response")
    
    # Running the coxph model on the imputed cohort data -----------------------
    
    mod.imputedcohort.y1 <- coxph(Surv(times, status1) ~ X11.pred + X21 + 
                                    X31.pred, data = cohort, robust = TRUE)
    
    mod.imputedcohort.y2 <- coxph(Surv(times, status2) ~ X12.pred + X32.pred, 
                                  data = cohort, robust = TRUE)
    
    # Building the auxiliary variables proposed by Breslow et al. (2019) -------
    
    A.y1.Breslow <- residuals(mod.imputedcohort.y1, type = "dfbeta", 
                              weighted = T)
    
    A.y2.Breslow <- residuals(mod.imputedcohort.y2, type = "dfbeta", 
                              weighted = T)
    
    # Building the auxiliary variables proposed by Shin et al. (2020) ----------
    
    time.on.study <- pmax(pmin(t2, cohort$times) - t1, 0)
    
    A             <- as.matrix(cbind(1, A.y1.Breslow))
    total         <- colSums(A)
    indiv.phase2  <- rownames(casecohort)
    calibration1  <- calibration(A.phase2 = A[indiv.phase2,], 
                                 design.weights = casecohort$weights, 
                                 total = total, eta0 = rep(0, ncol(A)), 
                                 niter.max = 10^3, epsilon.stop = 10^(-10))
    casecohort$weights.calib3 <- c(calibration1$calibrated.weights)
    
    mod.calib1.y1 <- coxph(Surv(times, status1) ~ X11 + X21 + X31, 
                           data = casecohort, weight = weights.calib3, 
                           id = id, robust = TRUE)
    
    A.Shin.y1 <- c(time.on.study * exp(mod.calib1.y1$coefficients %*% 
                                         t(cbind(cohort$X11.pred, cohort$X21, 
                                                 cohort$X31.pred))))
    
    # Calibrating the design weights against the auxiliary variables proposed by
    # Breslow et al. (2009) and Shin et al. (2020), focusing on event type 1 ---
    
    A             <- as.matrix(cbind(1, A.y1.Breslow, A.Shin.y1))
    total         <- colSums(A)
    calibration2  <- NA
    casecohort$weights.calib <- NA
    calibration2  <- calibration(A.phase2 = A[indiv.phase2,], 
                                 design.weights = casecohort$weights, 
                                 total = total, eta0 = rep(0, ncol(A)), 
                                 niter.max = 10^3, epsilon.stop = 10^(-10))
    casecohort$weights.calib <- c(calibration2$calibrated.weights)
    
    if(is.na(sum(casecohort$weights.calib))){
      A             <- as.matrix(cbind(1, A.y1.Breslow))
      total         <- colSums(A)
      calibration2  <- calibration(A.phase2 = A[indiv.phase2,], 
                                   design.weights = casecohort$weights, 
                                   total = total, eta0 = rep(0, ncol(A)), 
                                   niter.max = 10^3, epsilon.stop = 10^(-10))
      casecohort$weights.calib <- c(calibration2$calibrated.weights)
    }
    
    mod.calib.y1            <- coxph(Surv(times, status1) ~ X11 + X21 + X31, 
                                     data = casecohort, 
                                     weight = weights.calib, id = id, 
                                     robust = TRUE)
    
    mod.calib.y2            <- coxph(Surv(times, status2) ~ X12 + X32, 
                                     data = casecohort, 
                                     weight = weights.calib, id = id, 
                                     robust = TRUE)
    
    est.y1.calib <- influences.calib(mod = mod.calib.y1, A = A,
                                     x = c(-1, 1, -0.6), t1 = 0, t2 = 8)
    
    est.y2.calib <- influences.calib(mod = mod.calib.y2, A = A, x = c(0.5, 1), 
                                     t1 = 0, t2 = 8, competing = TRUE, 
                                     event.times.y1 = est.y1.casecohort$event.times)
    
    est.risk.calib <- risk.influences.calib(est.y1.calib, est.y2.calib)
    
    beta.y1.hat.calib <- mod.calib.y1$coefficients
    beta.y2.hat.calib <- mod.calib.y2$coefficients
    
    AR.y1.hat.calib <- est.risk.calib$AR.y1
    PR.y1.hat.calib <- est.risk.calib$PR.y1
    
    robvar.beta.y1.hat.calib   <- robustvariance(est.y1.calib$infl1.beta + est.y1.calib$infl2.beta) 
    robvar.beta.y2.hat.calib   <- robustvariance(est.y2.calib$infl1.beta + est.y2.calib$infl2.beta) 
    robvar.AR.y1.hat.calib     <- robustvariance(est.risk.calib$infl1.AR.y1 + est.risk.calib$infl2.AR.y1)
    robvar.PR.y1.hat.calib     <- robustvariance(est.risk.calib$infl1.PR.y1 + est.risk.calib$infl2.PR.y1)
    
    var.beta.y1.hat.calib      <- variance.calib(casecohort, 
                                                 est.y1.calib$infl1.beta,
                                                 est.y1.calib$infl2.beta,
                                                 cohort) 
    var.beta.y2.hat.calib      <- variance.calib(casecohort, 
                                                 est.y2.calib$infl1.beta,
                                                 est.y2.calib$infl2.beta,
                                                 cohort) 
    var.AR.y1.hat.calib        <- variance.calib(casecohort, 
                                                 est.risk.calib$infl1.AR.y1,
                                                 est.risk.calib$infl2.AR.y1,
                                                 cohort) 
    var.PR.y1.hat.calib        <- variance.calib(casecohort, 
                                                 est.risk.calib$infl1.PR.y1,
                                                 est.risk.calib$infl2.PR.y1,
                                                 cohort) 
    
    est.y2.calib   <- influences.calib(mod = mod.calib.y2, A = A, x = c(0.5, 1), 
                                       t1 = 0, t2 = 8)
    n.events.y2 <- length(est.y2.calib$Lambda0.t12.hat)
    Lambda0.t12.y2.hat.calib <- est.y2.calib$Lambda0.t12.hat[n.events.y2]
    robvar.Lambda0.t12.y2.hat.calib   <- robustvariance(est.y2.calib$infl1.Lambda0.t12[,n.events.y2] + est.y2.calib$infl2.Lambda0.t12[,n.events.y2]) 
    var.Lambda0.t12.y2.hat.calib      <- variance.calib(casecohort, 
                                                        as.matrix(est.y2.calib$infl1.Lambda0.t12[,n.events.y2]),
                                                        as.matrix(est.y2.calib$infl2.Lambda0.t12[,n.events.y2]),
                                                        cohort) 
    
    res.calib.rob <- cbind("CaseCohort Robust Calibration", t(beta.y1.hat.calib), 
                                t(beta.y2.hat.calib), 
                                Lambda0.t12.y2.hat.calib,
                                AR.y1.hat.calib, 
                                PR.y1.hat.calib, 
                                t(robvar.beta.y1.hat.calib),
                                t(robvar.beta.y2.hat.calib),
                                robvar.Lambda0.t12.y2.hat.calib,
                                robvar.AR.y1.hat.calib, 
                                robvar.PR.y1.hat.calib)
    
    res.calib <- cbind("CaseCohort Calibration", t(beta.y1.hat.calib), 
                            t(beta.y2.hat.calib), 
                            Lambda0.t12.y2.hat.calib,
                            AR.y1.hat.calib, 
                            PR.y1.hat.calib, t(var.beta.y1.hat.calib),
                            t(var.beta.y2.hat.calib),
                            var.Lambda0.t12.y2.hat.calib,
                            var.AR.y1.hat.calib, var.PR.y1.hat.calib)
    
    # --------------------------------------------------------------------------
    # Estimation using the case-cohort with calibrated weights, with calibration
    # against the two event types ----------------------------------------------
    
    # Building the auxiliary variables proposed by Shin et al. (2020) ----------
    
    A             <- as.matrix(cbind(1, A.y1.Breslow , A.y2.Breslow))
    total         <- colSums(A)
    indiv.phase2  <- rownames(casecohort)
    calibration2  <- NA
    casecohort$weights.calib1 <- NA
    calibration1  <- calibration(A.phase2 = A[indiv.phase2,], 
                                 design.weights = casecohort$weights, 
                                 total = total, eta0 = rep(0, ncol(A)), 
                                 niter.max = 10^3, epsilon.stop = 10^(-10))
    casecohort$weights.calib1 <- c(calibration1$calibrated.weights)
    
    if(is.na(sum(casecohort$weights.calib1))){
      A             <- as.matrix(cbind(1, A.y1.Breslow))
      total         <- colSums(A)
      calibration1  <- calibration(A.phase2 = A[indiv.phase2,], 
                                   design.weights = casecohort$weights, 
                                   total = total, eta0 = rep(0, ncol(A)), 
                                   niter.max = 10^3, epsilon.stop = 10^(-10))
      casecohort$weights.calib1 <- c(calibration1$calibrated.weights)
    }
    
    mod.calib1.y1 <- coxph(Surv(times, status1) ~ X11 + X21 + X31, 
                           data = casecohort, weight = weights.calib1, 
                           id = id, robust = TRUE)
    
    A.Shin.y1 <- c(time.on.study * exp(mod.calib1.y1$coefficients %*% 
                                         t(cbind(cohort$X11.pred, cohort$X21, 
                                                 cohort$X31.pred))))
    
    mod.calib1.y2 <- coxph(Surv(times, status2) ~ X12 + X32, 
                           data = casecohort, weight = weights.calib1, 
                           id = id, robust = TRUE)
    
    A.Shin.y2 <- c(time.on.study * exp(mod.calib1.y2$coefficients %*% 
                                         t(cbind(cohort$X12.pred, 
                                                 cohort$X32.pred))))
    
    # Calibrating the design weights against the auxiliary variables proposed by
    # Breslow et al. (2009) and Shin et al. (2020), focusing on event type 1 and
    # event type 2 -------------------------------------------------------------
    
    A             <- as.matrix(cbind(1, A.y1.Breslow, A.y2.Breslow, A.Shin.y1, 
                                     A.Shin.y2))
    total         <- colSums(A)
    calibration2  <- NA
    casecohort$weights.calib2 <- NA
    calibration2  <- calibration(A.phase2 = A[indiv.phase2,], 
                                 design.weights = casecohort$weights, 
                                 total = total, eta0 = rep(0, ncol(A)), 
                                 niter.max = 10^3, epsilon.stop = 10^(-10))
    casecohort$weights.calib2 <- c(calibration2$calibrated.weights)
    
    if(is.na(sum(casecohort$weights.calib2))){
      A             <- as.matrix(cbind(1, A.y1.Breslow, A.y2.Breslow, 
                                       A.Shin.y1))
      total         <- colSums(A)
      calibration2  <- calibration(A.phase2 = A[indiv.phase2,], 
                                   design.weights = casecohort$weights, 
                                   total = total, eta0 = rep(0, ncol(A)), 
                                   niter.max = 10^3, epsilon.stop = 10^(-10))
      casecohort$weights.calib2 <- c(calibration2$calibrated.weights)
    }
    if(is.na(sum(casecohort$weights.calib2))){
      A             <- as.matrix(cbind(1, A.y1.Breslow, A.y2.Breslow))
      total         <- colSums(A)
      calibration2  <- calibration(A.phase2 = A[indiv.phase2,], 
                                   design.weights = casecohort$weights, 
                                   total = total, eta0 = rep(0, ncol(A)), 
                                   niter.max = 10^3, epsilon.stop = 10^(-10))
      casecohort$weights.calib2 <- c(calibration2$calibrated.weights)
    }
    if(is.na(sum(casecohort$weights.calib2))){
      A             <- as.matrix(cbind(1, A.y1.Breslow))
      total         <- colSums(A)
      calibration2  <- calibration(A.phase2 = A[indiv.phase2,], 
                                   design.weights = casecohort$weights, 
                                   total = total, eta0 = rep(0, ncol(A)), 
                                   niter.max = 10^3, epsilon.stop = 10^(-10))
      casecohort$weights.calib2 <- c(calibration2$calibrated.weights)
    }
    
    mod.calib2.y1            <- coxph(Surv(times, status1) ~ X11 + X21 + X31, 
                                     data = casecohort, 
                                     weight = weights.calib2, id = id, 
                                     robust = TRUE)
    
    mod.calib2.y2            <- coxph(Surv(times, status2) ~ X12 + X32, 
                                     data = casecohort, 
                                     weight = weights.calib2, id = id, 
                                     robust = TRUE)
    
    est.y1.calib2 <- influences.calib(mod = mod.calib2.y1, A = A,
                                     x = c(-1, 1, -0.6), t1 = 0, t2 = 8)
    
    est.y2.calib2 <- influences.calib(mod = mod.calib2.y2, A = A, x = c(0.5, 1), 
                                     t1 = 0, t2 = 8, competing = TRUE, 
                                     event.times.y1 = est.y1.casecohort$event.times)
    
    est.risk.calib2 <- risk.influences.calib(est.y1.calib2, est.y2.calib2)
    
    beta.y1.hat.calib2 <- mod.calib2.y1$coefficients
    beta.y2.hat.calib2 <- mod.calib2.y2$coefficients
    
    AR.y1.hat.calib2 <- est.risk.calib2$AR.y1
    PR.y1.hat.calib2 <- est.risk.calib2$PR.y1
    
    robvar.beta.y1.hat.calib2   <- robustvariance(est.y1.calib2$infl1.beta + est.y1.calib2$infl2.beta) 
    robvar.beta.y2.hat.calib2   <- robustvariance(est.y2.calib2$infl1.beta + est.y2.calib2$infl2.beta) 
    robvar.AR.y1.hat.calib2     <- robustvariance(est.risk.calib2$infl1.AR.y1 + est.risk.calib2$infl2.AR.y1)
    robvar.PR.y1.hat.calib2     <- robustvariance(est.risk.calib2$infl1.PR.y1 + est.risk.calib2$infl2.PR.y1)
    
    var.beta.y1.hat.calib2      <- variance.calib(casecohort, 
                                                 est.y1.calib2$infl1.beta,
                                                 est.y1.calib2$infl2.beta,
                                                 cohort) 
    var.beta.y2.hat.calib2      <- variance.calib(casecohort, 
                                                 est.y2.calib2$infl1.beta,
                                                 est.y2.calib2$infl2.beta,
                                                 cohort) 
    var.AR.y1.hat.calib2        <- variance.calib(casecohort, 
                                                 est.risk.calib2$infl1.AR.y1,
                                                 est.risk.calib2$infl2.AR.y1,
                                                 cohort) 
    var.PR.y1.hat.calib2        <- variance.calib(casecohort, 
                                                 est.risk.calib2$infl1.PR.y1,
                                                 est.risk.calib2$infl2.PR.y1,
                                                 cohort) 
    
    est.y2.calib2   <- influences.calib(mod = mod.calib2.y2, A = A, x = c(0.5, 1), 
                                       t1 = 0, t2 = 8)
    n.events.y2 <- length(est.y2.calib2$Lambda0.t12.hat)
    Lambda0.t12.y2.hat.calib2 <- est.y2.calib2$Lambda0.t12.hat[n.events.y2]
    robvar.Lambda0.t12.y2.hat.calib2   <- robustvariance(est.y2.calib2$infl1.Lambda0.t12[,n.events.y2] + est.y2.calib2$infl2.Lambda0.t12[,n.events.y2]) 
    var.Lambda0.t12.y2.hat.calib2      <- variance.calib(casecohort, 
                                                        as.matrix(est.y2.calib2$infl1.Lambda0.t12[,n.events.y2]),
                                                        as.matrix(est.y2.calib2$infl2.Lambda0.t12[,n.events.y2]),
                                                        cohort) 
    
    res.calib2.rob <- cbind("CaseCohort Robust Calibration 2", t(beta.y1.hat.calib2), 
                           t(beta.y2.hat.calib2), 
                           Lambda0.t12.y2.hat.calib2,
                           AR.y1.hat.calib2, 
                           PR.y1.hat.calib2, 
                           t(robvar.beta.y1.hat.calib2),
                           t(robvar.beta.y2.hat.calib2),
                           robvar.Lambda0.t12.y2.hat.calib2,
                           robvar.AR.y1.hat.calib2, 
                           robvar.PR.y1.hat.calib2)
    
    res.calib2 <- cbind("CaseCohort Calibration 2", t(beta.y1.hat.calib2), 
                       t(beta.y2.hat.calib2), 
                       Lambda0.t12.y2.hat.calib2,
                       AR.y1.hat.calib2, 
                       PR.y1.hat.calib2, t(var.beta.y1.hat.calib2),
                       t(var.beta.y2.hat.calib2),
                       var.Lambda0.t12.y2.hat.calib2,
                       var.AR.y1.hat.calib2, var.PR.y1.hat.calib2)
    
    # --------------------------------------------------------------------------
    # Estimation using the case-cohort with calibrated weights, with calibration
    # against the event type of interest and Breslow (i) only ------------------
    
    # Building the auxiliary variables proposed by Shin et al. (2020) ----------
    
    A             <- as.matrix(cbind(1, A.y1.Breslow))
    total         <- colSums(A)
 
    mod.calib3.y1            <- coxph(Surv(times, status1) ~ X11 + X21 + X31, 
                                      data = casecohort, 
                                      weight = weights.calib3, id = id, 
                                      robust = TRUE)
    
    mod.calib3.y2            <- coxph(Surv(times, status2) ~ X12 + X32, 
                                      data = casecohort, 
                                      weight = weights.calib3, id = id, 
                                      robust = TRUE)
    
    est.y1.calib3 <- influences.calib(mod = mod.calib3.y1, A = A,
                                      x = c(-1, 1, -0.6), t1 = 0, t2 = 8)
    
    est.y2.calib3 <- influences.calib(mod = mod.calib3.y2, A = A, x = c(0.5, 1), 
                                      t1 = 0, t2 = 8, competing = TRUE, 
                                      event.times.y1 = est.y1.casecohort$event.times)
    
    est.risk.calib3 <- risk.influences.calib(est.y1.calib3, est.y2.calib3)
    
    beta.y1.hat.calib3 <- mod.calib3.y1$coefficients
    beta.y2.hat.calib3 <- mod.calib3.y2$coefficients
    
    AR.y1.hat.calib3 <- est.risk.calib3$AR.y1
    PR.y1.hat.calib3 <- est.risk.calib3$PR.y1
    
    robvar.beta.y1.hat.calib3   <- robustvariance(est.y1.calib3$infl1.beta + est.y1.calib3$infl2.beta) 
    robvar.beta.y2.hat.calib3   <- robustvariance(est.y2.calib3$infl1.beta + est.y2.calib3$infl2.beta) 
    robvar.AR.y1.hat.calib3     <- robustvariance(est.risk.calib3$infl1.AR.y1 + est.risk.calib3$infl2.AR.y1)
    robvar.PR.y1.hat.calib3     <- robustvariance(est.risk.calib3$infl1.PR.y1 + est.risk.calib3$infl2.PR.y1)
    
    var.beta.y1.hat.calib3      <- variance.calib(casecohort, 
                                                  est.y1.calib3$infl1.beta,
                                                  est.y1.calib3$infl2.beta,
                                                  cohort) 
    var.beta.y2.hat.calib3      <- variance.calib(casecohort, 
                                                  est.y2.calib3$infl1.beta,
                                                  est.y2.calib3$infl2.beta,
                                                  cohort) 
    var.AR.y1.hat.calib3        <- variance.calib(casecohort, 
                                                  est.risk.calib3$infl1.AR.y1,
                                                  est.risk.calib3$infl2.AR.y1,
                                                  cohort) 
    var.PR.y1.hat.calib3        <- variance.calib(casecohort, 
                                                  est.risk.calib3$infl1.PR.y1,
                                                  est.risk.calib3$infl2.PR.y1,
                                                  cohort) 
    
    est.y2.calib3   <- influences.calib(mod = mod.calib3.y2, A = A, x = c(0.5, 1), 
                                        t1 = 0, t2 = 8)
    n.events.y2 <- length(est.y2.calib3$Lambda0.t12.hat)
    Lambda0.t12.y2.hat.calib3 <- est.y2.calib3$Lambda0.t12.hat[n.events.y2]
    robvar.Lambda0.t12.y2.hat.calib3   <- robustvariance(est.y2.calib3$infl1.Lambda0.t12[,n.events.y2] + est.y2.calib3$infl2.Lambda0.t12[,n.events.y2]) 
    var.Lambda0.t12.y2.hat.calib3      <- variance.calib(casecohort, 
                                                         as.matrix(est.y2.calib3$infl1.Lambda0.t12[,n.events.y2]),
                                                         as.matrix(est.y2.calib3$infl2.Lambda0.t12[,n.events.y2]),
                                                         cohort) 
    
    res.calib3.rob <- cbind("CaseCohort Robust Calibration 3", t(beta.y1.hat.calib3), 
                            t(beta.y2.hat.calib3), 
                            Lambda0.t12.y2.hat.calib3,
                            AR.y1.hat.calib3, 
                            PR.y1.hat.calib3, 
                            t(robvar.beta.y1.hat.calib3),
                            t(robvar.beta.y2.hat.calib3),
                            robvar.Lambda0.t12.y2.hat.calib3,
                            robvar.AR.y1.hat.calib3, 
                            robvar.PR.y1.hat.calib3)
    
    res.calib3 <- cbind("CaseCohort Calibration 3", t(beta.y1.hat.calib3), 
                        t(beta.y2.hat.calib3), 
                        Lambda0.t12.y2.hat.calib3,
                        AR.y1.hat.calib3, 
                        PR.y1.hat.calib3, t(var.beta.y1.hat.calib3),
                        t(var.beta.y2.hat.calib3),
                        var.Lambda0.t12.y2.hat.calib3,
                        var.AR.y1.hat.calib3, var.PR.y1.hat.calib3)
    
    
    recap <- rbind(res.casecohort, res.casecohort.rob, res.calib, 
                   res.calib.rob, res.calib2, res.calib2.rob, res.calib3, 
                   res.calib3.rob, res.cohort)
    colnames(recap) <- c("Method", paste0(names.beta.y1, ".hat"),
                         paste0(names.beta.y2, ".hat"), "Lambda0.t12.y2.hat",
                         "AR.y1.hat", 
                         "PR.y1.hat", paste0("var.", names.beta.y1, ".hat"),
                         paste0("var.", names.beta.y2, ".hat"), 
                         "var.Lambda0.t12.y2.hat",
                         "var.AR.y1.hat", "var.PR.y1.hat")
    
    res <- rbind(res, cbind(recap, n, K, n.y0 = n.strata[1], 
                            n.y1 = n.strata[2], n.y2 = n.strata[3], 
                            m = m, ncasecohort = nrow(casecohort), 
                            Prob.y1, Prob.y2, beta11, beta21, beta31, beta12, 
                            beta32, AR.y1, PR.y1))
    row.names(res) <- NULL
  }
  myfile  <- paste0("SimulationResults_exhaustive-n", n, "-Prob.y1", Prob.y1, 
                    "-Prob.y2", Prob.y2, "-K", K, "Part", part, ".Rdata")
  save(res, file = myfile)
}

# ------------------------------------------------------------------------------
# Running the simulations in parallel for all the scenarios --------------------

# result <- mclapply(1:nrow(parameters), Onerun, mc.cores = 64)

# ------------------------------------------------------------------------------
# Combining the simulation results ---------------------------------------------

RES <- NULL
for (p in 1:nrow(parameters)) {
  
  Prob.y1           <- parameters[p, 1] # not accounting for censoring and competing events
  Prob.y2           <- parameters[p, 2] # not accounting for censoring and competing events
  n                 <- parameters[p, 3] 
  K                 <- parameters[p, 4] # number of non-cases we wish to sample for each case
  part              <- parameters[p, 5]
  
  load(paste0("SimulationResults_exhaustive-n", n, "-Prob.y1", Prob.y1, 
              "-Prob.y2", Prob.y2, "-K", K, "Part", part, ".Rdata"))
  RES <- rbind(RES, res)
  
}
RECAP           <- as.data.frame(RES)
ColNames        <- colnames(RECAP[,c(2:33)])
RECAP[ColNames] <- sapply(RECAP[ColNames], as.numeric)
RECAP$Method    <- as.factor(RECAP$Method)
RECAP$n         <- as.factor(RECAP$n)
RECAP$K         <- as.factor(RECAP$K)
RECAP$Prob.y1   <- as.factor(RECAP$Prob.y1)
RECAP$Prob.y2   <- as.factor(RECAP$Prob.y2)
myfile  <- paste0("SimulationResults_exhaustive.Rdata")
save(RECAP, file = myfile)

# ------------------------------------------------------------------------------
# Details of the results  ------------------------------------------------------

Nreplic         <- 5000
coverage.left   <- round(0.95 - 1.96 * sqrt(0.95 * 0.05 / Nreplic), digits = 3)
coverage.right  <- round(0.95 + 1.96 * sqrt(0.95 * 0.05 / Nreplic), digits = 3)

parameters <- parameters[,-5]
parameters <- parameters[which(duplicated.matrix(parameters)==FALSE),]

list.methods <- c("Cohort", "CaseCohort", "CaseCohort Robust", 
                  "CaseCohort Calibration", "CaseCohort Robust Calibration", 
                  "CaseCohort Calibration 2", "CaseCohort Robust Calibration 2",
                  "CaseCohort Calibration 3", "CaseCohort Robust Calibration 3")
list.methods.col <- c("Cohort", "CC", "CCRobust", "CCCalib", "CCCalibRobust",
                      "CCCalib2", "CCCalib2Robust", "CCCalib3", 
                      "CCCalib3Robust")
n.method <- length(list.methods)

details.AR.y1     <- NULL
details.PR.y1     <- NULL
details.logAR.y1     <- NULL
details.logPR.y1     <- NULL
for (i in 1:nrow(parameters)) {
  
  RECAP1 <- RECAP[((i-1) * (n.method * Nreplic) + 1):(i * (n.method * Nreplic)), ]
  
  beta11  <- RECAP1[1,]$beta11
  beta21  <- RECAP1[1,]$beta21
  beta31  <- RECAP1[1,]$beta31
  beta12  <- RECAP1[1,]$beta12
  beta32  <- RECAP1[1,]$beta32
 
  AR.y1   <- RECAP1[1,]$AR.y1
  PR.y1   <- RECAP1[1,]$PR.y1
  n.y0 <- mean(RECAP1$n.y0)
  n.y1 <- mean(RECAP1$n.y1)
  n.y2 <- mean(RECAP1$n.y2)
  m <- mean(RECAP1$m)
  ncasecohort     <- mean(RECAP1$ncasecohort)

  logAR.y1              <- log(AR.y1)
  RECAP1$logAR.y1.hat   <- log(RECAP1$AR.y1.hat)
  logPR.y1              <- log(PR.y1)
  RECAP1$logPR.y1.hat   <- log(RECAP1$PR.y1.hat)
  
  nameColResults <- function (name) {
    a <- NULL
    for (j in 1:length(list.methods)) {
      a <- c(a, paste0(name, ".", list.methods.col[j]))
    }
    return(a)
  }
  
  # Bias and mean estimate
  mean.AR.y1  <- NULL
  mean.PR.y1  <- NULL
  mean.logAR.y1  <- NULL
  mean.logPR.y1  <- NULL
  mean.logLambda0.t12.y2  <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    mean.AR.y1 <- cbind(mean.AR.y1, mean((RECAP1$AR.y1.hat[which(RECAP1$Method == method)])))
    mean.PR.y1 <- cbind(mean.PR.y1, mean((RECAP1$PR.y1.hat[which(RECAP1$Method == method)])))

    mean.logAR.y1 <- cbind(mean.logAR.y1, mean((RECAP1$logAR.y1.hat[which(RECAP1$Method == method)])))
    mean.logPR.y1 <- cbind(mean.logPR.y1, mean((RECAP1$logPR.y1.hat[which(RECAP1$Method == method)])))
  }

  mean.AR.y1 <- as.data.frame(mean.AR.y1)
  colnames(mean.AR.y1) <- nameColResults("mean.AR.y1")
  mean.PR.y1 <- as.data.frame(mean.PR.y1)
  colnames(mean.PR.y1) <- nameColResults("mean.PR.y1")

  mean.logAR.y1 <- as.data.frame(mean.logAR.y1)
  colnames(mean.logAR.y1) <- nameColResults("mean.logAR.y1")
  mean.logPR.y1 <- as.data.frame(mean.logPR.y1)
  colnames(mean.logPR.y1) <- nameColResults("mean.logPR.y1")
  
  bias.AR.y1 <- mean.AR.y1 - AR.y1
  colnames(bias.AR.y1) <- nameColResults("bias.AR.y1")
  bias.PR.y1 <- mean.PR.y1 - PR.y1
  colnames(bias.PR.y1) <- nameColResults("bias.PR.y1")
  
   bias.logAR.y1 <- mean.logAR.y1 - logAR.y1
  colnames(bias.logAR.y1) <- nameColResults("bias.logAR.y1")
  bias.logPR.y1 <- mean.logPR.y1 - logPR.y1
  colnames(bias.logPR.y1) <- nameColResults("bias.logPR.y1")
  
  # Relative bias
  relative.bias.AR.y1 <- abs(bias.AR.y1 / AR.y1)
  colnames(relative.bias.AR.y1) <- nameColResults("relative.bias.AR.y1")
  relative.bias.PR.y1 <- abs(bias.PR.y1 / PR.y1)
  colnames(relative.bias.PR.y1) <- nameColResults("relative.bias.PR.y1")
 
  relative.bias.logAR.y1 <- abs(bias.logAR.y1 / logAR.y1)
  colnames(relative.bias.logAR.y1) <- nameColResults("relative.bias.logAR.y1")
  relative.bias.logPR.y1 <- abs(bias.logPR.y1 / logPR.y1)
  colnames(relative.bias.logPR.y1) <- nameColResults("relative.bias.logPR.y1")
  
  # Empirical variance
  empir.var.AR.y1  <- NULL
  empir.var.PR.y1  <- NULL
  empir.var.logAR.y1  <- NULL
  empir.var.logPR.y1  <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    empir.var.AR.y1 <- cbind(empir.var.AR.y1, var(RECAP1$AR.y1.hat[which(RECAP1$Method == method)]))
    empir.var.PR.y1 <- cbind(empir.var.PR.y1, var(RECAP1$PR.y1.hat[which(RECAP1$Method == method)]))
    
    empir.var.logAR.y1 <- cbind(empir.var.logAR.y1, var(RECAP1$logAR.y1.hat[which(RECAP1$Method == method)]))
    empir.var.logPR.y1 <- cbind(empir.var.logPR.y1, var(RECAP1$logPR.y1.hat[which(RECAP1$Method == method)]))
  }
  empir.var.AR.y1 <- as.data.frame(empir.var.AR.y1)
  colnames(empir.var.AR.y1) <- nameColResults("empir.var.AR.y1")
  empir.var.PR.y1 <- as.data.frame(empir.var.PR.y1)
  colnames(empir.var.PR.y1) <- nameColResults("empir.var.PR.y1")
  
  empir.var.logAR.y1 <- as.data.frame(empir.var.logAR.y1)
  colnames(empir.var.logAR.y1) <- nameColResults("empir.var.logAR.y1")
  empir.var.logPR.y1 <- as.data.frame(empir.var.logPR.y1)
  colnames(empir.var.logPR.y1) <- nameColResults("empir.var.logPR.y1")

  ## Variance estimation for the log-transformed estimates
  var.logAR.y1.hat <- NULL
  var.logPR.y1.hat <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    var.logAR.y1.hat <- cbind(var.logAR.y1.hat, ((RECAP1$var.AR.y1.hat[which(RECAP1$Method == method)]) / (RECAP1$AR.y1.hat[which(RECAP1$Method == method)])^2))
    var.logPR.y1.hat <- cbind(var.logPR.y1.hat, ((RECAP1$var.PR.y1.hat[which(RECAP1$Method == method)]) / (RECAP1$PR.y1.hat[which(RECAP1$Method == method)])^2))
  }
  var.logAR.y1.hat <- as.data.frame(var.logAR.y1.hat)
  colnames(var.logAR.y1.hat) <- nameColResults("var.logAR.y1.hat")
  var.logPR.y1.hat <- as.data.frame(var.logPR.y1.hat)
  colnames(var.logPR.y1.hat) <- nameColResults("var.logPR.y1.hat")

  # Mean variance
  mean.var.AR.y1  <- NULL
  mean.var.PR.y1  <- NULL
  mean.var.logAR.y1  <- NULL
  mean.var.logPR.y1  <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    mean.var.AR.y1 <- cbind(mean.var.AR.y1, mean(RECAP1$var.AR.y1.hat[which(RECAP1$Method == method)]))
    mean.var.PR.y1 <- cbind(mean.var.PR.y1, mean(RECAP1$var.PR.y1.hat[which(RECAP1$Method == method)]))
    
    mean.var.logAR.y1 <- cbind(mean.var.logAR.y1, mean(var.logAR.y1.hat[,j]))
    mean.var.logPR.y1 <- cbind(mean.var.logPR.y1, mean(var.logPR.y1.hat[,j]))
  }
  mean.var.AR.y1 <- as.data.frame(mean.var.AR.y1)
  colnames(mean.var.AR.y1) <- nameColResults("mean.var.AR.y1")
  mean.var.PR.y1 <- as.data.frame(mean.var.PR.y1)
  colnames(mean.var.PR.y1) <- nameColResults("mean.var.PR.y1")
  
  mean.var.logAR.y1 <- as.data.frame(mean.var.logAR.y1)
  colnames(mean.var.logAR.y1) <- nameColResults("mean.var.logAR.y1")
  mean.var.logPR.y1 <- as.data.frame(mean.var.logPR.y1)
  colnames(mean.var.logPR.y1) <- nameColResults("mean.var.logPR.y1")
  
  # Relative difference between empirical variance and mean variance estimates
  relative.diffvar.AR.y1 <- abs(empir.var.AR.y1 - mean.var.AR.y1) / mean.var.AR.y1
  colnames(relative.diffvar.AR.y1) <- nameColResults("relative.diffvar.AR.y1")
  relative.diffvar.PR.y1 <- abs(empir.var.PR.y1 - mean.var.PR.y1) / mean.var.PR.y1
  colnames(relative.diffvar.PR.y1) <- nameColResults("relative.diffvar.PR.y1")
  
   relative.diffvar.logAR.y1 <- abs(empir.var.logAR.y1 - mean.var.logAR.y1) / mean.var.logAR.y1
  colnames(relative.diffvar.logAR.y1) <- nameColResults("relative.diffvar.logAR.y1")
  relative.diffvar.logPR.y1 <- abs(empir.var.logPR.y1 - mean.var.logPR.y1) / mean.var.logPR.y1
  colnames(relative.diffvar.logPR.y1) <- nameColResults("relative.diffvar.logPR.y1")
  
  # Construction of the CIs
  RECAP1$CI.left.AR.y1.hat <-  RECAP1$beta11
  RECAP1$CI.right.AR.y1.hat <- RECAP1$beta11
  RECAP1$CI.left.PR.y1.hat <-  RECAP1$beta11
  RECAP1$CI.right.PR.y1.hat <- RECAP1$beta11
  RECAP1$CI.left.logAR.y1.hat <-  RECAP1$beta11
  RECAP1$CI.right.logAR.y1.hat <- RECAP1$beta11
  RECAP1$CI.left.logPR.y1.hat <-  RECAP1$beta11
  RECAP1$CI.right.logPR.y1.hat <- RECAP1$beta11
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    RECAP1$CI.left.AR.y1.hat[which(RECAP1$Method == method)] <-    RECAP1$AR.y1.hat[which(RECAP1$Method == method)] - qnorm(0.975) * (sqrt(RECAP1$var.AR.y1.hat[which(RECAP1$Method == method)]))
    RECAP1$CI.right.AR.y1.hat[which(RECAP1$Method == method)] <-    RECAP1$AR.y1.hat[which(RECAP1$Method == method)] + qnorm(0.975) * (sqrt(RECAP1$var.AR.y1.hat[which(RECAP1$Method == method)]))
    RECAP1$CI.left.PR.y1.hat[which(RECAP1$Method == method)] <-    RECAP1$PR.y1.hat[which(RECAP1$Method == method)] - qnorm(0.975) * (sqrt(RECAP1$var.PR.y1.hat[which(RECAP1$Method == method)]))
    RECAP1$CI.right.PR.y1.hat[which(RECAP1$Method == method)] <-    RECAP1$PR.y1.hat[which(RECAP1$Method == method)] + qnorm(0.975) * (sqrt(RECAP1$var.PR.y1.hat[which(RECAP1$Method == method)]))
    
    RECAP1$CI.left.logAR.y1.hat[which(RECAP1$Method == method)] <-    RECAP1$logAR.y1.hat[which(RECAP1$Method == method)] - qnorm(0.975) * (sqrt(var.logAR.y1.hat[,j]))
    RECAP1$CI.right.logAR.y1.hat[which(RECAP1$Method == method)] <-    RECAP1$logAR.y1.hat[which(RECAP1$Method == method)] + qnorm(0.975) * (sqrt(var.logAR.y1.hat[,j]))
    RECAP1$CI.left.logPR.y1.hat[which(RECAP1$Method == method)] <-    RECAP1$logPR.y1.hat[which(RECAP1$Method == method)] - qnorm(0.975) * (sqrt(var.logPR.y1.hat[,j]))
    RECAP1$CI.right.logPR.y1.hat[which(RECAP1$Method == method)] <-    RECAP1$logPR.y1.hat[which(RECAP1$Method == method)] + qnorm(0.975) * (sqrt(var.logPR.y1.hat[,j]))
     }
  
  # Coverage of the CIs
  coverage.AR.y1  <- NULL
  coverage.PR.y1  <- NULL
  coverage.logAR.y1  <- NULL
  coverage.logPR.y1  <- NULL
  for (j in 1:length(list.methods)) {
    method <- list.methods[j]
    coverage.AR.y1 <- cbind(coverage.AR.y1, mean((RECAP1$CI.left.AR.y1.hat[which(RECAP1$Method == method)] < AR.y1)&(RECAP1$CI.right.AR.y1.hat[which(RECAP1$Method == method)] > AR.y1)))
    coverage.PR.y1 <- cbind(coverage.PR.y1, mean((RECAP1$CI.left.PR.y1.hat[which(RECAP1$Method == method)] < PR.y1)&(RECAP1$CI.right.PR.y1.hat[which(RECAP1$Method == method)] > PR.y1)))
    
    coverage.logAR.y1 <- cbind(coverage.logAR.y1, mean((RECAP1$CI.left.logAR.y1.hat[which(RECAP1$Method == method)] < logAR.y1)&(RECAP1$CI.right.logAR.y1.hat[which(RECAP1$Method == method)] > logAR.y1)))
    coverage.logPR.y1 <- cbind(coverage.logPR.y1, mean((RECAP1$CI.left.logPR.y1.hat[which(RECAP1$Method == method)] < logPR.y1)&(RECAP1$CI.right.logPR.y1.hat[which(RECAP1$Method == method)] > logPR.y1)))
   }
  coverage.AR.y1 <- as.data.frame(coverage.AR.y1)
  colnames(coverage.AR.y1) <- nameColResults("coverage.AR.y1")
  coverage.PR.y1 <- as.data.frame(coverage.PR.y1)
  colnames(coverage.PR.y1) <- nameColResults("coverage.PR.y1")
  
  coverage.logAR.y1 <- as.data.frame(coverage.logAR.y1)
  colnames(coverage.logAR.y1) <- nameColResults("coverage.logAR.y1")
  coverage.logPR.y1 <- as.data.frame(coverage.logPR.y1)
  colnames(coverage.logPR.y1) <- nameColResults("coverage.logPR.y1")  
  
  # nominal coverage
  nominal.coverage.AR.y1  <- NULL
  nominal.coverage.PR.y1  <- NULL
  nominal.coverage.logAR.y1  <- NULL
  nominal.coverage.logPR.y1  <- NULL
  for (j in 1:length(list.methods)) {
     nominal.coverage.AR.y1 <- cbind(nominal.coverage.AR.y1, 
                                     (coverage.AR.y1[j] > coverage.left) & (coverage.AR.y1[j] < coverage.right))
    nominal.coverage.PR.y1 <- cbind(nominal.coverage.PR.y1, 
                                    (coverage.PR.y1[j] > coverage.left) & (coverage.PR.y1[j] < coverage.right))

     nominal.coverage.logAR.y1 <- cbind(nominal.coverage.logAR.y1, 
                                    (coverage.logAR.y1[j] > coverage.left) & (coverage.logAR.y1[j] < coverage.right))
    nominal.coverage.logPR.y1 <- cbind(nominal.coverage.logPR.y1, 
                                    (coverage.logPR.y1[j] > coverage.left) & (coverage.logPR.y1[j] < coverage.right))  
    }
  nominal.coverage.AR.y1 <- as.data.frame(nominal.coverage.AR.y1)
  colnames(nominal.coverage.AR.y1) <- nameColResults("nominal.coverage.AR.y1")
  nominal.coverage.PR.y1 <- as.data.frame(nominal.coverage.PR.y1)
  colnames(nominal.coverage.PR.y1) <- nameColResults("nominal.coverage.PR.y1")

  nominal.coverage.logAR.y1 <- as.data.frame(nominal.coverage.logAR.y1)
  colnames(nominal.coverage.logAR.y1) <- nameColResults("nominal.coverage.logAR.y1")
  nominal.coverage.logPR.y1 <- as.data.frame(nominal.coverage.logPR.y1)
  colnames(nominal.coverage.logPR.y1) <- nameColResults("nominal.coverage.logPR.y1")  
  
  # efficiency gain
  efficiency.AR.y1  <- as.numeric(empir.var.AR.y1[1]) / empir.var.AR.y1
  efficiency.PR.y1  <- as.numeric(empir.var.PR.y1[1]) / empir.var.PR.y1
  efficiency.logAR.y1  <- as.numeric(empir.var.logAR.y1[1]) / empir.var.logAR.y1
  efficiency.logPR.y1  <- as.numeric(empir.var.logPR.y1[1]) / empir.var.logPR.y1

  efficiency.AR.y1 <- as.data.frame(efficiency.AR.y1)
  colnames(efficiency.AR.y1) <- nameColResults("efficiency.AR.y1")
  efficiency.PR.y1 <- as.data.frame(efficiency.PR.y1)
  colnames(efficiency.PR.y1) <- nameColResults("efficiency.PR.y1")
  
  efficiency.logAR.y1 <- as.data.frame(efficiency.logAR.y1)
  colnames(efficiency.logAR.y1) <- nameColResults("efficiency.logAR.y1")
  efficiency.logPR.y1 <- as.data.frame(efficiency.logPR.y1)
  colnames(efficiency.logPR.y1) <- nameColResults("efficiency.logPR.y1")  

  details.AR.y1   <- rbind(details.AR.y1, c(bias.AR.y1, relative.bias.AR.y1,
                                            empir.var.AR.y1, 
                                            mean.var.AR.y1, 
                                            relative.diffvar.AR.y1,
                                            coverage.AR.y1,
                                            nominal.coverage.AR.y1,
                                            efficiency.AR.y1,
                                            beta11 = beta11, 
                                            beta21 = beta21, beta31 = beta31, 
                                            beta12 = beta12, beta32 = beta32,
                                              AR.y1 = AR.y1, PR.y1 = PR.y1,
                                            logAR.y1 = logAR.y1, 
                                            logPR.y1 = logPR.y1,
                                            n = as.character(RECAP1[1,]$n), 
                                            K = as.character(RECAP1[1,]$K), 
                                            Prob.y1 = as.character(RECAP1[1,]$Prob.y1), 
                                            Prob.y2 = as.character(RECAP1[1,]$Prob.y2), 
                                            ncasecohort = ncasecohort,
                                            mean.n.y0 = n.y0,
                                            mean.n.y1 = n.y1,
                                            mean.n.y2 = n.y2,
                                            mean.m = m))
  
  details.PR.y1   <- rbind(details.PR.y1, c(bias.PR.y1, relative.bias.PR.y1,
                                            empir.var.PR.y1,
                                            mean.var.PR.y1, 
                                            relative.diffvar.PR.y1,
                                            coverage.PR.y1,
                                            nominal.coverage.PR.y1,
                                            efficiency.PR.y1,
                                            beta11 = beta11, 
                                            beta21 = beta21, beta31 = beta31, 
                                            beta12 = beta12, beta32 = beta32,
                                            AR.y1 = AR.y1, PR.y1 = PR.y1,
                                            logAR.y1 = logAR.y1, 
                                            logPR.y1 = logPR.y1,
                                            n = as.character(RECAP1[1,]$n), 
                                            K = as.character(RECAP1[1,]$K), 
                                            Prob.y1 = as.character(RECAP1[1,]$Prob.y1), 
                                            Prob.y2 = as.character(RECAP1[1,]$Prob.y2), 
                                            ncasecohort = ncasecohort,
                                            mean.n.y0 = n.y0,
                                            mean.n.y1 = n.y1,
                                            mean.n.y2 = n.y2,
                                            mean.m = m))
  
  details.logAR.y1   <- rbind(details.logAR.y1, c(bias.logAR.y1, 
                                                  relative.bias.logAR.y1,
                                                  empir.var.logAR.y1, 
                                            mean.var.logAR.y1, 
                                            relative.diffvar.logAR.y1,
                                            coverage.logAR.y1,
                                            nominal.coverage.logAR.y1,
                                            efficiency.logAR.y1,
                                            beta11 = beta11, 
                                            beta21 = beta21, beta31 = beta31, 
                                            beta12 = beta12, beta32 = beta32,
                                              AR.y1 = AR.y1, PR.y1 = PR.y1,
                                            logAR.y1 = logAR.y1, 
                                            logPR.y1 = logPR.y1,
                                            n = as.character(RECAP1[1,]$n), 
                                            K = as.character(RECAP1[1,]$K), 
                                            Prob.y1 = as.character(RECAP1[1,]$Prob.y1), 
                                            Prob.y2 = as.character(RECAP1[1,]$Prob.y2), 
                                            ncasecohort = ncasecohort,
                                            mean.n.y0 = n.y0,
                                            mean.n.y1 = n.y1,
                                            mean.n.y2 = n.y2,
                                            mean.m = m))
  
  details.logPR.y1   <- rbind(details.logPR.y1, c(bias.logPR.y1, 
                                                  relative.bias.logPR.y1,
                                                  empir.var.logPR.y1, 
                                                  mean.var.logPR.y1, 
                                                  relative.diffvar.logPR.y1,
                                                  coverage.logPR.y1,
                                                  nominal.coverage.logPR.y1,
                                                  efficiency.logPR.y1,
                                                  beta11 = beta11, 
                                                  beta21 = beta21, beta31 = beta31, 
                                                  beta12 = beta12, beta32 = beta32,
                                                  AR.y1 = AR.y1, PR.y1 = PR.y1,
                                                  logAR.y1 = logAR.y1, 
                                                  logPR.y1 = logPR.y1,
                                                  n = as.character(RECAP1[1,]$n), 
                                                  K = as.character(RECAP1[1,]$K), 
                                                  Prob.y1 = as.character(RECAP1[1,]$Prob.y1), 
                                                  Prob.y2 = as.character(RECAP1[1,]$Prob.y2), 
                                                  ncasecohort = ncasecohort,
                                                  mean.n.y0 = n.y0,
                                                  mean.n.y1 = n.y1,
                                                  mean.n.y2 = n.y2,
                                                  mean.m = m))
}
details.AR.y1 <- as.data.frame(details.AR.y1)
details.PR.y1 <- as.data.frame(details.PR.y1)
details.logAR.y1 <- as.data.frame(details.logAR.y1)
details.logPR.y1 <- as.data.frame(details.logPR.y1)

simul.name = "_exhausive"

save(details.AR.y1, file = paste0("details.AR.y1", simul.name, ".Rdata"))
col.names.AR.y1 <- colnames(details.AR.y1[,c(1:90)])
details.AR.y1[col.names.AR.y1] <- sapply(details.AR.y1[col.names.AR.y1], 
                                           as.numeric)
write.csv(details.AR.y1, file = paste0("details.AR.y1", simul.name, ".csv"))

save(details.PR.y1, file = paste0("details.PR.y1", simul.name, ".Rdata"))
col.names.PR.y1 <- colnames(details.PR.y1[,c(1:90)])
details.PR.y1[col.names.PR.y1] <- sapply(details.PR.y1[col.names.PR.y1], 
                                         as.numeric)
write.csv(details.PR.y1, file = paste0("details.PR.y1", simul.name, ".csv"))

save(details.logAR.y1, file = paste0("details.logAR.y1", simul.name, ".Rdata"))
col.names.logAR.y1 <- colnames(details.logAR.y1[,c(1:90)])
details.logAR.y1[col.names.logAR.y1] <- sapply(details.logAR.y1[col.names.logAR.y1], 
                                         as.numeric)
write.csv(details.logAR.y1, file = paste0("details.logAR.y1", simul.name, ".csv"))

save(details.logPR.y1, file = paste0("details.logPR.y1", simul.name, ".Rdata"))
col.names.logPR.y1 <- colnames(details.logPR.y1[,c(1:90)])
details.logPR.y1[col.names.logPR.y1] <- sapply(details.logPR.y1[col.names.logPR.y1], 
                                         as.numeric)
write.csv(details.logPR.y1, file = paste0("details.logPR.y1", simul.name, ".csv"))


## Print Latex tables for manuscript -------------------------------------------

### No calibration
a <- cbind(details.logAR.y1[details.logAR.y1$K == 1, c(84, 85, 80, 2, 20, 29:30, 47:48, 65)]) # case cohort with K = 1, no calibration
b <- cbind(details.logAR.y1[details.logAR.y1$K == 2, c(84, 85, 80, 2, 20, 29:30, 47:48, 65)]) # case cohort with K = 2, no calibration

a_scientific        <- a
a_scientific[,4:7]  <- data.frame(lapply(a[,4:7], function(x) format(x, scientific = TRUE, digits = 2)))
a_scientific[,10]    <- format(a[,10], digits = 2)
a                   <- data.frame(lapply(a_scientific , function(x) format(x, digits = 3)))
print(xtable(a), include.rownames=FALSE) # K = 1, no calibration

b_scientific        <- b 
b_scientific[,4:7]  <- data.frame(lapply(b[,4:7], function(x) format(x, scientific = TRUE, digits = 2)))
b_scientific[,10]    <- format(b[,10], digits = 2)
b                   <- data.frame(lapply(b_scientific , function(x) format(x, digits = 3)))
print(xtable(b), include.rownames=FALSE) # K = 2, no calibration

### Calibration
a <- cbind(details.logAR.y1[details.logAR.y1$K == 1, c(84, 85, 80, 4, 22, 31:32, 49:50, 67, 69, 71)]) # case cohort with K = 1 and calibration
b <- cbind(details.logAR.y1[details.logAR.y1$K == 2, c(84, 85, 80, 4, 22, 31:32, 49:50, 67, 69, 71)]) # case cohort with K = 2 and calibration

a_scientific        <- a
a_scientific[,4:7]  <- data.frame(lapply(a[,4:7], function(x) format(x, scientific = TRUE, digits = 2)))
a_scientific[,10:12]  <- data.frame(lapply(a[,10:12], function(x) format(x, digits = 2)))
a                   <- data.frame(lapply(a_scientific , function(x) format(x, digits = 3)))
print(xtable(a), include.rownames=FALSE) # K = 1 and calibration

b_scientific        <- b 
b_scientific[,4:7]  <- data.frame(lapply(b[,4:7], function(x) format(x, scientific = TRUE, digits = 2)))
b_scientific[,10:12]<- data.frame(lapply(b[,10:12], function(x) format(x, digits = 2)))
b                   <- data.frame(lapply(b_scientific , function(x) format(x, digits = 3)))
print(xtable(b), include.rownames=FALSE) # K = 2 and calibration

### Cohort
d <- cbind(details.logAR.y1[details.logAR.y1$K == 1, c(84, 85, 80, 1, 19, 28, 46)])

d_scientific        <- d
d_scientific[,4:6]  <- data.frame(lapply(d[,4:6], function(x) format(x, scientific = TRUE, digits = 2)))
d                   <- data.frame(lapply(d_scientific , function(x) format(x, digits = 3)))
print(xtable(d), include.rownames=FALSE) # K = 1

### Subcohort and case-cohort size
a <- cbind(details.logAR.y1[, c(84, 85, 83, 90, 86)]) 
a[,4:5] <- data.frame(lapply(a[,4:5] , function(x) format(x, digits = 1)))
print(xtable(a), include.rownames=FALSE) 

