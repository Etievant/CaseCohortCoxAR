## -----------------------------------------------------------------------------
## Description: This script performs estimation of cause-specific absolute risks
##              in PLCO when using a case-cohort with 
##              - stratified sampling based on case status
##              - sampling of only a fraction (m.y1) of the cases with an event 
##              of interest (y1) and only a fraction (m.y2) of the cases with a 
##              competing event (y2). m.y0 and m.y2 are based on on n.y0 and 
##              n.y2, respectively
##              - design weights 
##
##              We compare the robust variance estimate and the variance 
##              estimate that properly accounts for the case-cohort sampling, as
##              proposed in Etievant and Gail (2024)
##
##              Results are displayed in Web Appendix D. For the analyses 
##              displayed in the Main Document (with m.y0 and m.y2 both based on
##               n.y1), see Rscript in PLCO analysis.R 
## -----------------------------------------------------------------------------
## Required Package: dplyr, parallel, survival, xtable
## -----------------------------------------------------------------------------
## Required Functions: help.functions.R
## -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Read functions and packages --------------------------------------------------

library(dplyr)
library(survival)
library(parallel)
library(xtable)

source("help.functions.R")

# ------------------------------------------------------------------------------
# Load data set ----------------------------------------------------------------

load("entirecohortPLCO_competingrisks.Rdata")
n     <- nrow(cohort)

# ------------------------------------------------------------------------------
# Specify fixed parameters -----------------------------------------------------

# Time interval for the absolute risk, e.g., [0, 15] 
t1  <- 0
t2  <- 10

# Covariate profile for the absolute and pure risks 
x.y1 <- c(2.5, 1, 8)
x.y2 <- c(1, 0, 28, 1, 0, 1, 0, 0, 1, 1, 60)

names.x.y1  <- c("x66","x7","x8")
names.x.y2  <- c("x1","x2","x5","x9","x10","x11","x12","x13","x14","x16","x18")

# x1 current cigarette smoker
# x2 former cigarette smoker
# x5 BMI at 50 yo
# x66 log psa level from most recent test prior to diagnosis
# x7 had prostatectomy
# x8 gleason score from prostatectomy or biopsy
# x9 indicator of being non-Hispanic White
# x10 indicator of being non-Hispanic Black
# x11 had hypertension
# x12 had stroke
# x13 had heart attack
# x14 had diabetes
# x16 had pulmonary morbidity (chronic bronchitis or emphysema)
# x18 age at diagnosis

# ------------------------------------------------------------------------------
# Estimation using the entire cohort -------------------------------------------

mod.cohort.y1  <- coxph(Surv(yearsonstudy, status1) ~ x66 + x7 + 
                          x8, data = cohort, robust = TRUE, id = id) # treating event 2 as independent censoring

mod.cohort.y2  <- coxph(Surv(yearsonstudy, status2) ~ x1 + x2 + x5 + x9 + 
                          x10 + x11 + x12 + x13 + x14 + x16 + x18, 
                        data = cohort, robust = TRUE, id = id) # treating event 1 as independent censoring

est.y1.cohort <- influences(mod = mod.cohort.y1, x = x.y1, t1 = t1, 
                            t2 = t2)
est.y2.cohort <- influences(mod = mod.cohort.y2, x = x.y2, t1 = t1, 
                            t2 = t2, competing = TRUE, 
                            event.times.y1 = est.y1.cohort$event.times)

est.risk.cohort <- risk.influences(est.y1.cohort, est.y2.cohort)

beta.y1.hat.cohort <- mod.cohort.y1$coefficients
beta.y2.hat.cohort <- mod.cohort.y2$coefficients

AR.y1.hat.cohort <- est.risk.cohort$AR.y1
PR.y1.hat.cohort <- est.risk.cohort$PR.y1

logAR.y1.hat.cohort <- log(AR.y1.hat.cohort)
logPR.y1.hat.cohort <- log(PR.y1.hat.cohort)

robvar.beta.y1.hat.cohort   <- robustvariance(est.y1.cohort$infl.beta) # equal to diag(mod.cohort.y1$var)
robvar.beta.y2.hat.cohort   <- robustvariance(est.y2.cohort$infl.beta) # equal to diag(mod.cohort.y2$var)
robvar.AR.y1.hat.cohort     <- robustvariance(est.risk.cohort$infl.AR.y1)
robvar.PR.y1.hat.cohort     <- robustvariance(est.risk.cohort$infl.PR.y1)

robvar.logAR.y1.hat.cohort <- robvar.AR.y1.hat.cohort / (AR.y1.hat.cohort)^2
robvar.logPR.y1.hat.cohort <- robvar.PR.y1.hat.cohort / (PR.y1.hat.cohort)^2

# CI.left.AR.y1.hat.cohort.rob    <- AR.y1.hat.cohort - qnorm(0.975) * (sqrt(robvar.AR.y1.hat.cohort))
# CI.right.AR.y1.hat.cohort.rob   <- AR.y1.hat.cohort + qnorm(0.975) * (sqrt(robvar.AR.y1.hat.cohort))
# CI.left.PR.y1.hat.cohort.rob    <- PR.y1.hat.cohort - qnorm(0.975) * (sqrt(robvar.PR.y1.hat.cohort))
# CI.right.PR.y1.hat.cohort.rob   <- PR.y1.hat.cohort + qnorm(0.975) * (sqrt(robvar.PR.y1.hat.cohort))
# 
# CI.left.logAR.y1.hat.cohort.rob    <- logAR.y1.hat.cohort - qnorm(0.975) * (sqrt(robvar.logAR.y1.hat.cohort))
# CI.right.logAR.y1.hat.cohort.rob   <- logAR.y1.hat.cohort + qnorm(0.975) * (sqrt(robvar.logAR.y1.hat.cohort))
# CI.left.logPR.y1.hat.cohort.rob    <- logPR.y1.hat.cohort - qnorm(0.975) * (sqrt(robvar.logPR.y1.hat.cohort))
# CI.right.logPR.y1.hat.cohort.rob   <- logPR.y1.hat.cohort + qnorm(0.975) * (sqrt(robvar.logPR.y1.hat.cohort))

# --------------------------------------------------------------------------
# Sampling of the case cohort, when m.y1 = 157, m.y0 = 200, and m.y2 = 196 -----

set.seed(12346)

# Number of individuals to be sampled from the 3 strata based on case-status

n.y1            <- sum(cohort$status == 1)
n.y2            <- sum(cohort$status == 2)
n.y0            <- sum(cohort$status == 0)
n.strata        <- c(n.y0, n.y1, n.y2) 

m.y1            <- ceiling(n.y1 / 2)
m.y0            <- ceiling(n.y0 / 10)
m.y2            <- ceiling(n.y2 / 10)
m.strata        <- c(m.y0, m.y1, m.y2)

cohort          <- cbind(cohort, n.strata = n.strata[cohort$status + 1], 
                         m.strata = m.strata[cohort$status + 1])

indiv.sampled <- sample(cohort$id[cohort$status == 0], size = m.y0, 
                        replace = FALSE)
indiv.sampled <- c(indiv.sampled, sample(cohort$id[cohort$status == 2], 
                                         size = m.y2, replace = FALSE))
indiv.sampled <- c(indiv.sampled, sample(cohort$id[cohort$status == 1], 
                                         size = m.y1, replace = FALSE))

phase2                <- rep(0, n)
names(phase2)         <- cohort$id
phase2[indiv.sampled] <- 1

cohort$weights <- n.strata[cohort$status + 1] / m.strata[cohort$status + 1]

cohort <- cbind(cohort, phase2)

casecohort <- cohort
casecohort <- casecohort[-which((phase2 == 0)), ]


# ------------------------------------------------------------------------------
# At risk indicator matrix and counting process matrix for the phase-one data at
# all event times, even for cases not in phase-two -----------------------------

mod.cohort.y1.detail   <- coxph.detail(mod.cohort.y1, riskmat = TRUE) 
mod.cohort.y2.detail   <- coxph.detail(mod.cohort.y2, riskmat = TRUE)   

riskmat.phase1.y1      <- mod.cohort.y1.detail$riskmat
riskmat.phase1.y2      <- mod.cohort.y2.detail$riskmat

rownames(riskmat.phase1.y1)  <- cohort[, "id"]
rownames(riskmat.phase1.y2)  <- cohort[, "id"]

observed.times.phase1.y1 <- apply(riskmat.phase1.y1, 1,
                                  function(v) {which.max(cumsum(v))})
observed.times.phase1.y2 <- apply(riskmat.phase1.y2, 1,
                                  function(v) {which.max(cumsum(v))})

dNt.phase1.y1            <- matrix(0, nrow(riskmat.phase1.y1), 
                                   ncol(riskmat.phase1.y1))
dNt.phase1.y1[cbind(1:nrow(riskmat.phase1.y1), observed.times.phase1.y1)] <- 1
dNt.phase1.y1            <- sweep(dNt.phase1.y1, 1, cohort$status1, "*")
colnames(dNt.phase1.y1)  <- colnames(riskmat.phase1.y1)
rownames(dNt.phase1.y1)  <- rownames(riskmat.phase1.y1)

dNt.phase1.y2            <- matrix(0, nrow(riskmat.phase1.y2), 
                                   ncol(riskmat.phase1.y2))
dNt.phase1.y2[cbind(1:nrow(riskmat.phase1.y2), observed.times.phase1.y2)] <- 1
dNt.phase1.y2            <- sweep(dNt.phase1.y2, 1, cohort$status2, "*")
colnames(dNt.phase1.y2)  <- colnames(riskmat.phase1.y2)
rownames(dNt.phase1.y2)  <- rownames(riskmat.phase1.y2)


# ------------------------------------------------------------------------------
# Estimation using the case-cohort with design weights, when the numerator of 
# the Breslow estimator is not weighted ----------------------------------------

mod.casecohort.y1  <- coxph(Surv(yearsonstudy, status1) ~ x66 +  x7 + 
                              x8, data = casecohort, robust = TRUE, 
                            id = id, weights = weights) # treating event 2 as independent censoring

mod.casecohort.y2  <- coxph(Surv(yearsonstudy, status2) ~ x1 + x2 + x5 + x9 + 
                              x10 + x11 + x12 + x13 + x14 + x16 + x18, 
                            data = casecohort, robust = TRUE, id = id,
                            weights = weights) # treating event 1 as independent censoring

est.y1.casecohort <- influences.generalized(mod = mod.casecohort.y1, x = x.y1, 
                                            t1 = t1, t2 = t2, 
                                            Breslow.all = TRUE,
                                            riskmat.phase1 = riskmat.phase1.y1, 
                                            dNt.phase1 = dNt.phase1.y1)

est.y2.casecohort <- influences.generalized(mod = mod.casecohort.y2, x = x.y2, 
                                            t1 = t1, t2 = t2, competing = TRUE, 
                                            event.times.y1 = est.y1.casecohort$event.times,
                                            riskmat.phase1 = riskmat.phase1.y2, 
                                            dNt.phase1 = dNt.phase1.y2, 
                                            Breslow.all = TRUE)

est.risk.casecohort <- risk.influences.generalized(est.y1.casecohort, 
                                                   est.y2.casecohort)

beta.y1.hat.casecohort <- mod.casecohort.y1$coefficients
beta.y2.hat.casecohort <- mod.casecohort.y2$coefficients

AR.y1.hat.casecohort <- est.risk.casecohort$AR.y1
PR.y1.hat.casecohort <- est.risk.casecohort$PR.y1

logAR.y1.hat.casecohort <- log(AR.y1.hat.casecohort)
logPR.y1.hat.casecohort <- log(PR.y1.hat.casecohort)

robvar.beta.y1.hat.casecohort   <- robustvariance(est.y1.casecohort$infl1.beta + est.y1.casecohort$infl2.beta) # equal to diag(mod.casecohort.y1$var)
robvar.beta.y2.hat.casecohort   <- robustvariance(est.y2.casecohort$infl1.beta + est.y2.casecohort$infl2.beta) # equal to diag(mod.casecohort.y2$var)
robvar.AR.y1.hat.casecohort     <- robustvariance(est.risk.casecohort$infl1.AR.y1 + est.risk.casecohort$infl2.AR.y1)
robvar.PR.y1.hat.casecohort     <- robustvariance(est.risk.casecohort$infl1.PR.y1 + est.risk.casecohort$infl2.PR.y1)

var.beta.y1.hat.casecohort      <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.y1.casecohort$infl1.beta, 
                                                        infl2 = est.y1.casecohort$infl2.beta, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 
var.beta.y2.hat.casecohort      <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.y2.casecohort$infl1.beta, 
                                                        infl2 = est.y2.casecohort$infl2.beta, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 
var.AR.y1.hat.casecohort        <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.risk.casecohort$infl1.AR.y1, 
                                                        infl2 = est.risk.casecohort$infl2.AR.y1, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 
var.PR.y1.hat.casecohort        <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.risk.casecohort$infl1.PR.y1, 
                                                        infl2 = est.risk.casecohort$infl2.PR.y1, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 

var.logAR.y1.hat.casecohort <- var.AR.y1.hat.casecohort / (AR.y1.hat.casecohort)^2
var.logPR.y1.hat.casecohort <- var.PR.y1.hat.casecohort / (PR.y1.hat.casecohort)^2

robvar.logAR.y1.hat.casecohort <- robvar.AR.y1.hat.casecohort / (AR.y1.hat.casecohort)^2
robvar.logPR.y1.hat.casecohort <- robvar.PR.y1.hat.casecohort / (PR.y1.hat.casecohort)^2

# print the results in a table

AR.y1.hat.cohort <- round(AR.y1.hat.cohort, digits = 3)
PR.y1.hat.cohort <- round(PR.y1.hat.cohort, digits = 3)

AR.y1.hat.casecohort <- round(AR.y1.hat.casecohort, digits = 3)
PR.y1.hat.casecohort <- round(PR.y1.hat.casecohort, digits = 3)

robvar.AR.y1.hat.cohort <- round(robvar.AR.y1.hat.cohort, digits = 6)
robvar.PR.y1.hat.cohort <- round(robvar.PR.y1.hat.cohort, digits = 6)

robvar.AR.y1.hat.casecohort <- round(robvar.AR.y1.hat.casecohort, digits = 5)
robvar.PR.y1.hat.casecohort <- round(robvar.PR.y1.hat.casecohort, digits = 5)

var.AR.y1.hat.casecohort <- round(var.AR.y1.hat.casecohort, digits = 5)
var.PR.y1.hat.casecohort <- round(var.PR.y1.hat.casecohort, digits = 5)

AR.est <- c(AR.y1.hat.casecohort, var.AR.y1.hat.casecohort,
            robvar.AR.y1.hat.casecohort)

PR.est <- c(PR.y1.hat.casecohort, var.PR.y1.hat.casecohort,
            robvar.PR.y1.hat.casecohort)

# ------------------------------------------------------------------------------
# Estimation using the case-cohort with design weights and when the numerator of
# the Breslow estimator is weighted --------------------------------------------

est.y1.casecohort2 <- influences(mod = mod.casecohort.y1, x = x.y1, t1 = t1, 
                                 t2 = t2)

est.y2.casecohort2 <- influences(mod = mod.casecohort.y2, x = x.y2, t1 = t1, 
                                 t2 = t2, competing = TRUE, 
                                 event.times.y1 = est.y1.casecohort2$event.times)

est.risk.casecohort2 <- risk.influences(est.y1.casecohort2, est.y2.casecohort2)

AR.y1.hat.casecohort2 <- est.risk.casecohort2$AR.y1
PR.y1.hat.casecohort2 <- est.risk.casecohort2$PR.y1

robvar.beta.y1.hat.casecohort2   <- robustvariance(est.y1.casecohort2$infl.beta) # equal to diag(mod.casecohort.y1$var)
robvar.beta.y2.hat.casecohort2   <- robustvariance(est.y2.casecohort2$infl.beta) # equal to diag(mod.casecohort.y2$var)
robvar.AR.y1.hat.casecohort2     <- robustvariance(est.risk.casecohort2$infl.AR.y1)
robvar.PR.y1.hat.casecohort2     <- robustvariance(est.risk.casecohort2$infl.PR.y1)

var.beta.y1.hat.casecohort2      <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.y1.casecohort2$infl.beta, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 
var.beta.y2.hat.casecohort2      <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.y2.casecohort2$infl.beta, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 
var.AR.y1.hat.casecohort2        <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.risk.casecohort2$infl.AR.y1, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 
var.PR.y1.hat.casecohort2        <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.risk.casecohort2$infl.PR.y1, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 

# ------------------------------------------------------------------------------
# Print the results in a TeX table ---------------------------------------------

AR.y1.hat.casecohort2 <- round(AR.y1.hat.casecohort2, digits = 3)
PR.y1.hat.casecohort2 <- round(PR.y1.hat.casecohort2, digits = 3)

robvar.AR.y1.hat.casecohort2 <- round(robvar.AR.y1.hat.casecohort2, digits = 5)
robvar.PR.y1.hat.casecohort2 <- round(robvar.PR.y1.hat.casecohort2, digits = 5)

var.AR.y1.hat.casecohort2 <- round(var.AR.y1.hat.casecohort2, digits = 5)
var.PR.y1.hat.casecohort2 <- round(var.PR.y1.hat.casecohort2, digits = 5)

AR.est2 <- c(AR.y1.hat.casecohort2, 
            var.AR.y1.hat.casecohort2,robvar.AR.y1.hat.casecohort2)

PR.est2 <- c(PR.y1.hat.casecohort2, 
            var.PR.y1.hat.casecohort2, robvar.PR.y1.hat.casecohort2)

print(xtable(cbind(rbind(AR.est, PR.est), rbind(AR.est2, PR.est2)), 
             digits = -2), include.rownames = FALSE)


# ------------------------------------------------------------------------------
# Using another covariate profile for the absolute and pure risks --------------
x.y1 <- c(1.6, 1, 7)
x.y2 <- c(0, 1, 25, 0, 1, 0, 1, 0, 1, 0, 55)

names.x.y1  <- c("x66","x7","x8")
names.x.y2  <- c("x1","x2","x5","x9","x10","x11","x12","x13","x14","x16","x18")

# --------------------------------------------------------------------------
# Estimation using the case-cohort with design weights, when the numerator of 
# the Breslow estimator is not weighted ----------------------------------------

mod.casecohort.y1  <- coxph(Surv(yearsonstudy, status1) ~ x66 +  x7 + 
                              x8, data = casecohort, robust = TRUE, 
                            id = id, weights = weights) # treating event 2 as independent censoring

mod.casecohort.y2  <- coxph(Surv(yearsonstudy, status2) ~ x1 + x2 + x5 + x9 + 
                              x10 + x11 + x12 + x13 + x14 + x16 + x18, 
                            data = casecohort, robust = TRUE, id = id,
                            weights = weights) # treating event 1 as independent censoring

est.y1.casecohort <- influences.generalized(mod = mod.casecohort.y1, x = x.y1, 
                                            t1 = t1, t2 = t2, 
                                            Breslow.all = TRUE,
                                            riskmat.phase1 = riskmat.phase1.y1, 
                                            dNt.phase1 = dNt.phase1.y1)

est.y2.casecohort <- influences.generalized(mod = mod.casecohort.y2, x = x.y2, 
                                            t1 = t1, t2 = t2, competing = TRUE, 
                                            event.times.y1 = est.y1.casecohort$event.times,
                                            riskmat.phase1 = riskmat.phase1.y2, 
                                            dNt.phase1 = dNt.phase1.y2, 
                                            Breslow.all = TRUE)

est.risk.casecohort <- risk.influences.generalized(est.y1.casecohort, 
                                                   est.y2.casecohort)

beta.y1.hat.casecohort <- mod.casecohort.y1$coefficients
beta.y2.hat.casecohort <- mod.casecohort.y2$coefficients

AR.y1.hat.casecohort <- est.risk.casecohort$AR.y1
PR.y1.hat.casecohort <- est.risk.casecohort$PR.y1

logAR.y1.hat.casecohort <- log(AR.y1.hat.casecohort)
logPR.y1.hat.casecohort <- log(PR.y1.hat.casecohort)

robvar.beta.y1.hat.casecohort   <- robustvariance(est.y1.casecohort$infl1.beta + est.y1.casecohort$infl2.beta) # equal to diag(mod.casecohort.y1$var)
robvar.beta.y2.hat.casecohort   <- robustvariance(est.y2.casecohort$infl1.beta + est.y2.casecohort$infl2.beta) # equal to diag(mod.casecohort.y2$var)
robvar.AR.y1.hat.casecohort     <- robustvariance(est.risk.casecohort$infl1.AR.y1 + est.risk.casecohort$infl2.AR.y1)
robvar.PR.y1.hat.casecohort     <- robustvariance(est.risk.casecohort$infl1.PR.y1 + est.risk.casecohort$infl2.PR.y1)

var.beta.y1.hat.casecohort      <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.y1.casecohort$infl1.beta, 
                                                        infl2 = est.y1.casecohort$infl2.beta, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 
var.beta.y2.hat.casecohort      <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.y2.casecohort$infl1.beta, 
                                                        infl2 = est.y2.casecohort$infl2.beta, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 
var.AR.y1.hat.casecohort        <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.risk.casecohort$infl1.AR.y1, 
                                                        infl2 = est.risk.casecohort$infl2.AR.y1, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 
var.PR.y1.hat.casecohort        <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.risk.casecohort$infl1.PR.y1, 
                                                        infl2 = est.risk.casecohort$infl2.PR.y1, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 

var.logAR.y1.hat.casecohort <- var.AR.y1.hat.casecohort / (AR.y1.hat.casecohort)^2
var.logPR.y1.hat.casecohort <- var.PR.y1.hat.casecohort / (PR.y1.hat.casecohort)^2

robvar.logAR.y1.hat.casecohort <- robvar.AR.y1.hat.casecohort / (AR.y1.hat.casecohort)^2
robvar.logPR.y1.hat.casecohort <- robvar.PR.y1.hat.casecohort / (PR.y1.hat.casecohort)^2

# print the results in a table

AR.y1.hat.casecohort <- round(AR.y1.hat.casecohort, digits = 3)
PR.y1.hat.casecohort <- round(PR.y1.hat.casecohort, digits = 3)

robvar.AR.y1.hat.casecohort <- round(robvar.AR.y1.hat.casecohort, digits = 6)
robvar.PR.y1.hat.casecohort <- round(robvar.PR.y1.hat.casecohort, digits = 6)

var.AR.y1.hat.casecohort <- round(var.AR.y1.hat.casecohort, digits = 6)
var.PR.y1.hat.casecohort <- round(var.PR.y1.hat.casecohort, digits = 6)

AR.est <- c(AR.y1.hat.casecohort, var.AR.y1.hat.casecohort,
            robvar.AR.y1.hat.casecohort)

PR.est <- c(PR.y1.hat.casecohort, var.PR.y1.hat.casecohort,
            robvar.PR.y1.hat.casecohort)

# ------------------------------------------------------------------------------
# Estimation using the case-cohort with design weights and when the numerator of
# the Breslow estimator is weighted --------------------------------------------

est.y1.casecohort2 <- influences(mod = mod.casecohort.y1, x = x.y1, t1 = t1, 
                                 t2 = t2)

est.y2.casecohort2 <- influences(mod = mod.casecohort.y2, x = x.y2, t1 = t1, 
                                 t2 = t2, competing = TRUE, 
                                 event.times.y1 = est.y1.casecohort2$event.times)

est.risk.casecohort2 <- risk.influences(est.y1.casecohort2, est.y2.casecohort2)

AR.y1.hat.casecohort2 <- est.risk.casecohort2$AR.y1
PR.y1.hat.casecohort2 <- est.risk.casecohort2$PR.y1

robvar.beta.y1.hat.casecohort2   <- robustvariance(est.y1.casecohort2$infl.beta) # equal to diag(mod.casecohort.y1$var)
robvar.beta.y2.hat.casecohort2   <- robustvariance(est.y2.casecohort2$infl.beta) # equal to diag(mod.casecohort.y2$var)
robvar.AR.y1.hat.casecohort2     <- robustvariance(est.risk.casecohort2$infl.AR.y1)
robvar.PR.y1.hat.casecohort2     <- robustvariance(est.risk.casecohort2$infl.PR.y1)

var.beta.y1.hat.casecohort2      <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.y1.casecohort2$infl.beta, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 
var.beta.y2.hat.casecohort2      <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.y2.casecohort2$infl.beta, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 
var.AR.y1.hat.casecohort2        <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.risk.casecohort2$infl.AR.y1, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 
var.PR.y1.hat.casecohort2        <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.risk.casecohort2$infl.PR.y1, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 

# ------------------------------------------------------------------------------
# Print the results in a TeX table ---------------------------------------------

AR.y1.hat.casecohort2 <- round(AR.y1.hat.casecohort2, digits = 3)
PR.y1.hat.casecohort2 <- round(PR.y1.hat.casecohort2, digits = 3)

robvar.AR.y1.hat.casecohort2 <- round(robvar.AR.y1.hat.casecohort2, digits = 6)
robvar.PR.y1.hat.casecohort2 <- round(robvar.PR.y1.hat.casecohort2, digits = 6)

var.AR.y1.hat.casecohort2 <- round(var.AR.y1.hat.casecohort2, digits = 6)
var.PR.y1.hat.casecohort2 <- round(var.PR.y1.hat.casecohort2, digits = 6)

AR.est2 <- c(AR.y1.hat.casecohort2, 
             var.AR.y1.hat.casecohort2,robvar.AR.y1.hat.casecohort2)

PR.est2 <- c(PR.y1.hat.casecohort2, 
             var.PR.y1.hat.casecohort2, robvar.PR.y1.hat.casecohort2)

print(xtable(cbind(rbind(AR.est, PR.est), rbind(AR.est2, PR.est2)), 
             digits = -2), include.rownames = FALSE)



# ------------------------------------------------------------------------------
# Using another covariate profile for the absolute and pure risks --------------
x.y1 <- c(0, 0, 5)
x.y2 <- c(0, 1, 32, 0, 1, 0, 1, 0, 1, 0, 60)

names.x.y1  <- c("x66","x7","x8")
names.x.y2  <- c("x1","x2","x5","x9","x10","x11","x12","x13","x14","x16","x18")

# --------------------------------------------------------------------------
# Estimation using the case-cohort with design weights, when the numerator of 
# the Breslow estimator is not weighted ----------------------------------------

mod.casecohort.y1  <- coxph(Surv(yearsonstudy, status1) ~ x66 +  x7 + 
                              x8, data = casecohort, robust = TRUE, 
                            id = id, weights = weights) # treating event 2 as independent censoring

mod.casecohort.y2  <- coxph(Surv(yearsonstudy, status2) ~ x1 + x2 + x5 + x9 + 
                              x10 + x11 + x12 + x13 + x14 + x16 + x18, 
                            data = casecohort, robust = TRUE, id = id,
                            weights = weights) # treating event 1 as independent censoring

est.y1.casecohort <- influences.generalized(mod = mod.casecohort.y1, x = x.y1, 
                                            t1 = t1, t2 = t2, 
                                            Breslow.all = TRUE,
                                            riskmat.phase1 = riskmat.phase1.y1, 
                                            dNt.phase1 = dNt.phase1.y1)

est.y2.casecohort <- influences.generalized(mod = mod.casecohort.y2, x = x.y2, 
                                            t1 = t1, t2 = t2, competing = TRUE, 
                                            event.times.y1 = est.y1.casecohort$event.times,
                                            riskmat.phase1 = riskmat.phase1.y2, 
                                            dNt.phase1 = dNt.phase1.y2, 
                                            Breslow.all = TRUE)

est.risk.casecohort <- risk.influences.generalized(est.y1.casecohort, 
                                                   est.y2.casecohort)

beta.y1.hat.casecohort <- mod.casecohort.y1$coefficients
beta.y2.hat.casecohort <- mod.casecohort.y2$coefficients

AR.y1.hat.casecohort <- est.risk.casecohort$AR.y1
PR.y1.hat.casecohort <- est.risk.casecohort$PR.y1

logAR.y1.hat.casecohort <- log(AR.y1.hat.casecohort)
logPR.y1.hat.casecohort <- log(PR.y1.hat.casecohort)

robvar.beta.y1.hat.casecohort   <- robustvariance(est.y1.casecohort$infl1.beta + est.y1.casecohort$infl2.beta) # equal to diag(mod.casecohort.y1$var)
robvar.beta.y2.hat.casecohort   <- robustvariance(est.y2.casecohort$infl1.beta + est.y2.casecohort$infl2.beta) # equal to diag(mod.casecohort.y2$var)
robvar.AR.y1.hat.casecohort     <- robustvariance(est.risk.casecohort$infl1.AR.y1 + est.risk.casecohort$infl2.AR.y1)
robvar.PR.y1.hat.casecohort     <- robustvariance(est.risk.casecohort$infl1.PR.y1 + est.risk.casecohort$infl2.PR.y1)

var.beta.y1.hat.casecohort      <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.y1.casecohort$infl1.beta, 
                                                        infl2 = est.y1.casecohort$infl2.beta, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 
var.beta.y2.hat.casecohort      <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.y2.casecohort$infl1.beta, 
                                                        infl2 = est.y2.casecohort$infl2.beta, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 
var.AR.y1.hat.casecohort        <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.risk.casecohort$infl1.AR.y1, 
                                                        infl2 = est.risk.casecohort$infl2.AR.y1, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 
var.PR.y1.hat.casecohort        <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.risk.casecohort$infl1.PR.y1, 
                                                        infl2 = est.risk.casecohort$infl2.PR.y1, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 

var.logAR.y1.hat.casecohort <- var.AR.y1.hat.casecohort / (AR.y1.hat.casecohort)^2
var.logPR.y1.hat.casecohort <- var.PR.y1.hat.casecohort / (PR.y1.hat.casecohort)^2

robvar.logAR.y1.hat.casecohort <- robvar.AR.y1.hat.casecohort / (AR.y1.hat.casecohort)^2
robvar.logPR.y1.hat.casecohort <- robvar.PR.y1.hat.casecohort / (PR.y1.hat.casecohort)^2

# print the results in a table

AR.y1.hat.casecohort <- round(AR.y1.hat.casecohort, digits = 4)
PR.y1.hat.casecohort <- round(PR.y1.hat.casecohort, digits = 4)

robvar.AR.y1.hat.casecohort <- round(robvar.AR.y1.hat.casecohort, digits = 7)
robvar.PR.y1.hat.casecohort <- round(robvar.PR.y1.hat.casecohort, digits = 7)

var.AR.y1.hat.casecohort <- round(var.AR.y1.hat.casecohort, digits = 7)
var.PR.y1.hat.casecohort <- round(var.PR.y1.hat.casecohort, digits = 7)

AR.est <- c(AR.y1.hat.casecohort, var.AR.y1.hat.casecohort,
            robvar.AR.y1.hat.casecohort)

PR.est <- c(PR.y1.hat.casecohort, var.PR.y1.hat.casecohort,
            robvar.PR.y1.hat.casecohort)

# ------------------------------------------------------------------------------
# Estimation using the case-cohort with design weights and when the numerator of
# the Breslow estimator is weighted --------------------------------------------

est.y1.casecohort2 <- influences(mod = mod.casecohort.y1, x = x.y1, t1 = t1, 
                                 t2 = t2)

est.y2.casecohort2 <- influences(mod = mod.casecohort.y2, x = x.y2, t1 = t1, 
                                 t2 = t2, competing = TRUE, 
                                 event.times.y1 = est.y1.casecohort2$event.times)

est.risk.casecohort2 <- risk.influences(est.y1.casecohort2, est.y2.casecohort2)

AR.y1.hat.casecohort2 <- est.risk.casecohort2$AR.y1
PR.y1.hat.casecohort2 <- est.risk.casecohort2$PR.y1

robvar.beta.y1.hat.casecohort2   <- robustvariance(est.y1.casecohort2$infl.beta) # equal to diag(mod.casecohort.y1$var)
robvar.beta.y2.hat.casecohort2   <- robustvariance(est.y2.casecohort2$infl.beta) # equal to diag(mod.casecohort.y2$var)
robvar.AR.y1.hat.casecohort2     <- robustvariance(est.risk.casecohort2$infl.AR.y1)
robvar.PR.y1.hat.casecohort2     <- robustvariance(est.risk.casecohort2$infl.PR.y1)

var.beta.y1.hat.casecohort2      <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.y1.casecohort2$infl.beta, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 
var.beta.y2.hat.casecohort2      <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.y2.casecohort2$infl.beta, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 
var.AR.y1.hat.casecohort2        <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.risk.casecohort2$infl.AR.y1, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 
var.PR.y1.hat.casecohort2        <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.risk.casecohort2$infl.PR.y1, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 

# ------------------------------------------------------------------------------
# Print the results in a TeX table ---------------------------------------------

AR.y1.hat.casecohort2 <- round(AR.y1.hat.casecohort2, digits = 4)
PR.y1.hat.casecohort2 <- round(PR.y1.hat.casecohort2, digits = 4)

robvar.AR.y1.hat.casecohort2 <- round(robvar.AR.y1.hat.casecohort2, digits = 7)
robvar.PR.y1.hat.casecohort2 <- round(robvar.PR.y1.hat.casecohort2, digits = 7)

var.AR.y1.hat.casecohort2 <- round(var.AR.y1.hat.casecohort2, digits = 7)
var.PR.y1.hat.casecohort2 <- round(var.PR.y1.hat.casecohort2, digits = 7)

AR.est2 <- c(AR.y1.hat.casecohort2, 
             var.AR.y1.hat.casecohort2,robvar.AR.y1.hat.casecohort2)

PR.est2 <- c(PR.y1.hat.casecohort2, 
             var.PR.y1.hat.casecohort2, robvar.PR.y1.hat.casecohort2)

print(xtable(cbind(rbind(AR.est, PR.est), rbind(AR.est2, PR.est2)), 
             digits = -2), include.rownames = FALSE)



# --------------------------------------------------------------------------
# Sampling of the case cohort, when m.y1 = 313, m.y0 = 400, and m.y2 = 391 -----

set.seed(12346)

# Number of individuals to be sampled from the 3 strata based on case-status

n.y1            <- sum(cohort$status == 1)
n.y2            <- sum(cohort$status == 2)
n.y0            <- sum(cohort$status == 0)
n.strata        <- c(n.y0, n.y1, n.y2) 

m.y1            <- ceiling(n.y1)
m.y0            <- ceiling(n.y0 / 5)
m.y2            <- ceiling(n.y2 / 5)
m.strata        <- c(m.y0, m.y1, m.y2)

cohort          <- cohort[,-c(30:33)]
cohort          <- cbind(cohort, n.strata = n.strata[cohort$status + 1], 
                         m.strata = m.strata[cohort$status + 1])

indiv.sampled <- sample(cohort$id[cohort$status == 0], size = m.y0, 
                        replace = FALSE)
indiv.sampled <- c(indiv.sampled, sample(cohort$id[cohort$status == 2], 
                                         size = m.y2, replace = FALSE))
indiv.sampled <- c(indiv.sampled, sample(cohort$id[cohort$status == 1], 
                                         size = m.y1, replace = FALSE))

phase2                <- rep(0, n)
names(phase2)         <- cohort$id
phase2[indiv.sampled] <- 1

cohort$weights <- n.strata[cohort$status + 1] / m.strata[cohort$status + 1]

cohort <- cbind(cohort, phase2)

casecohort <- cohort
casecohort <- casecohort[-which((phase2 == 0)), ]

# ------------------------------------------------------------------------------
# Covariate profile for the absolute and pure risks ----------------------------
x.y1 <- c(2.5, 1, 8)
x.y2 <- c(1, 0, 28, 1, 0, 1, 0, 0, 1, 1, 60)

names.x.y1  <- c("x66","x7","x8")
# ------------------------------------------------------------------------------
# Estimation using the case-cohort with design weights, when the numerator of 
# the Breslow estimator is not weighted ----------------------------------------

mod.casecohort.y1  <- coxph(Surv(yearsonstudy, status1) ~ x66 +  x7 + 
                              x8, data = casecohort, robust = TRUE, 
                            id = id, weights = weights) # treating event 2 as independent censoring

mod.casecohort.y2  <- coxph(Surv(yearsonstudy, status2) ~ x1 + x2 + x5 + x9 + 
                              x10 + x11 + x12 + x13 + x14 + x16 + x18, 
                            data = casecohort, robust = TRUE, id = id,
                            weights = weights) # treating event 1 as independent censoring

est.y1.casecohort <- influences.generalized(mod = mod.casecohort.y1, x = x.y1, 
                                            t1 = t1, t2 = t2, 
                                            Breslow.all = TRUE,
                                            riskmat.phase1 = riskmat.phase1.y1, 
                                            dNt.phase1 = dNt.phase1.y1)

est.y2.casecohort <- influences.generalized(mod = mod.casecohort.y2, x = x.y2, 
                                            t1 = t1, t2 = t2, competing = TRUE, 
                                            event.times.y1 = est.y1.casecohort$event.times,
                                            riskmat.phase1 = riskmat.phase1.y2, 
                                            dNt.phase1 = dNt.phase1.y2, 
                                            Breslow.all = TRUE)

est.risk.casecohort <- risk.influences.generalized(est.y1.casecohort, 
                                                   est.y2.casecohort)

beta.y1.hat.casecohort <- mod.casecohort.y1$coefficients
beta.y2.hat.casecohort <- mod.casecohort.y2$coefficients

AR.y1.hat.casecohort <- est.risk.casecohort$AR.y1
PR.y1.hat.casecohort <- est.risk.casecohort$PR.y1

logAR.y1.hat.casecohort <- log(AR.y1.hat.casecohort)
logPR.y1.hat.casecohort <- log(PR.y1.hat.casecohort)

robvar.beta.y1.hat.casecohort   <- robustvariance(est.y1.casecohort$infl1.beta + est.y1.casecohort$infl2.beta) # equal to diag(mod.casecohort.y1$var)
robvar.beta.y2.hat.casecohort   <- robustvariance(est.y2.casecohort$infl1.beta + est.y2.casecohort$infl2.beta) # equal to diag(mod.casecohort.y2$var)
robvar.AR.y1.hat.casecohort     <- robustvariance(est.risk.casecohort$infl1.AR.y1 + est.risk.casecohort$infl2.AR.y1)
robvar.PR.y1.hat.casecohort     <- robustvariance(est.risk.casecohort$infl1.PR.y1 + est.risk.casecohort$infl2.PR.y1)

var.beta.y1.hat.casecohort      <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.y1.casecohort$infl1.beta, 
                                                        infl2 = est.y1.casecohort$infl2.beta, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 
var.beta.y2.hat.casecohort      <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.y2.casecohort$infl1.beta, 
                                                        infl2 = est.y2.casecohort$infl2.beta, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 
var.AR.y1.hat.casecohort        <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.risk.casecohort$infl1.AR.y1, 
                                                        infl2 = est.risk.casecohort$infl2.AR.y1, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 
var.PR.y1.hat.casecohort        <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.risk.casecohort$infl1.PR.y1, 
                                                        infl2 = est.risk.casecohort$infl2.PR.y1, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 

var.logAR.y1.hat.casecohort <- var.AR.y1.hat.casecohort / (AR.y1.hat.casecohort)^2
var.logPR.y1.hat.casecohort <- var.PR.y1.hat.casecohort / (PR.y1.hat.casecohort)^2

robvar.logAR.y1.hat.casecohort <- robvar.AR.y1.hat.casecohort / (AR.y1.hat.casecohort)^2
robvar.logPR.y1.hat.casecohort <- robvar.PR.y1.hat.casecohort / (PR.y1.hat.casecohort)^2

# print the results in a table

AR.y1.hat.cohort <- round(AR.y1.hat.cohort, digits = 3)
PR.y1.hat.cohort <- round(PR.y1.hat.cohort, digits = 3)

AR.y1.hat.casecohort <- round(AR.y1.hat.casecohort, digits = 3)
PR.y1.hat.casecohort <- round(PR.y1.hat.casecohort, digits = 3)

robvar.AR.y1.hat.casecohort <- round(robvar.AR.y1.hat.casecohort, digits = 6)
robvar.PR.y1.hat.casecohort <- round(robvar.PR.y1.hat.casecohort, digits = 5)

var.AR.y1.hat.casecohort <- round(var.AR.y1.hat.casecohort, digits = 6)
var.PR.y1.hat.casecohort <- round(var.PR.y1.hat.casecohort, digits = 5)

AR.est <- c(AR.y1.hat.casecohort, var.AR.y1.hat.casecohort,
            robvar.AR.y1.hat.casecohort)

PR.est <- c(PR.y1.hat.casecohort, var.PR.y1.hat.casecohort,
            robvar.PR.y1.hat.casecohort)

# ------------------------------------------------------------------------------
# Estimation using the case-cohort with design weights and when the numerator of
# the Breslow estimator is weighted --------------------------------------------

est.y1.casecohort2 <- influences(mod = mod.casecohort.y1, x = x.y1, t1 = t1, 
                                 t2 = t2)

est.y2.casecohort2 <- influences(mod = mod.casecohort.y2, x = x.y2, t1 = t1, 
                                 t2 = t2, competing = TRUE, 
                                 event.times.y1 = est.y1.casecohort2$event.times)

est.risk.casecohort2 <- risk.influences(est.y1.casecohort2, est.y2.casecohort2)

AR.y1.hat.casecohort2 <- est.risk.casecohort2$AR.y1
PR.y1.hat.casecohort2 <- est.risk.casecohort2$PR.y1

robvar.beta.y1.hat.casecohort2   <- robustvariance(est.y1.casecohort2$infl.beta) # equal to diag(mod.casecohort.y1$var)
robvar.beta.y2.hat.casecohort2   <- robustvariance(est.y2.casecohort2$infl.beta) # equal to diag(mod.casecohort.y2$var)
robvar.AR.y1.hat.casecohort2     <- robustvariance(est.risk.casecohort2$infl.AR.y1)
robvar.PR.y1.hat.casecohort2     <- robustvariance(est.risk.casecohort2$infl.PR.y1)

var.beta.y1.hat.casecohort2      <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.y1.casecohort2$infl.beta, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 
var.beta.y2.hat.casecohort2      <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.y2.casecohort2$infl.beta, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 
var.AR.y1.hat.casecohort2        <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.risk.casecohort2$infl.AR.y1, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 
var.PR.y1.hat.casecohort2        <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.risk.casecohort2$infl.PR.y1, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 

# ------------------------------------------------------------------------------
# Print the results in a TeX table ---------------------------------------------

AR.y1.hat.casecohort2 <- round(AR.y1.hat.casecohort2, digits = 4)
PR.y1.hat.casecohort2 <- round(PR.y1.hat.casecohort2, digits = 4)

robvar.AR.y1.hat.casecohort2 <- round(robvar.AR.y1.hat.casecohort2, digits = 6)
robvar.PR.y1.hat.casecohort2 <- round(robvar.PR.y1.hat.casecohort2, digits = 5)

var.AR.y1.hat.casecohort2 <- round(var.AR.y1.hat.casecohort2, digits = 6)
var.PR.y1.hat.casecohort2 <- round(var.PR.y1.hat.casecohort2, digits = 5)

AR.est2 <- c(AR.y1.hat.casecohort2, 
             var.AR.y1.hat.casecohort2,robvar.AR.y1.hat.casecohort2)

PR.est2 <- c(PR.y1.hat.casecohort2, 
             var.PR.y1.hat.casecohort2, robvar.PR.y1.hat.casecohort2)

print(xtable(cbind(rbind(AR.est, PR.est), rbind(AR.est2, PR.est2)), 
             digits = -2), include.rownames = FALSE)


# ------------------------------------------------------------------------------
# Using another covariate profile for the absolute and pure risks --------------
x.y1 <- c(1.6, 1, 7)
x.y2 <- c(0, 1, 25, 0, 1, 0, 1, 0, 1, 0, 55)

names.x.y1  <- c("x66","x7","x8")
names.x.y2  <- c("x1","x2","x5","x9","x10","x11","x12","x13","x14","x16","x18")

# --------------------------------------------------------------------------
# Estimation using the case-cohort with design weights, when the numerator of 
# the Breslow estimator is not weighted ----------------------------------------

mod.casecohort.y1  <- coxph(Surv(yearsonstudy, status1) ~ x66 +  x7 + 
                              x8, data = casecohort, robust = TRUE, 
                            id = id, weights = weights) # treating event 2 as independent censoring

mod.casecohort.y2  <- coxph(Surv(yearsonstudy, status2) ~ x1 + x2 + x5 + x9 + 
                              x10 + x11 + x12 + x13 + x14 + x16 + x18, 
                            data = casecohort, robust = TRUE, id = id,
                            weights = weights) # treating event 1 as independent censoring

est.y1.casecohort <- influences.generalized(mod = mod.casecohort.y1, x = x.y1, 
                                            t1 = t1, t2 = t2, 
                                            Breslow.all = TRUE,
                                            riskmat.phase1 = riskmat.phase1.y1, 
                                            dNt.phase1 = dNt.phase1.y1)

est.y2.casecohort <- influences.generalized(mod = mod.casecohort.y2, x = x.y2, 
                                            t1 = t1, t2 = t2, competing = TRUE, 
                                            event.times.y1 = est.y1.casecohort$event.times,
                                            riskmat.phase1 = riskmat.phase1.y2, 
                                            dNt.phase1 = dNt.phase1.y2, 
                                            Breslow.all = TRUE)

est.risk.casecohort <- risk.influences.generalized(est.y1.casecohort, 
                                                   est.y2.casecohort)

beta.y1.hat.casecohort <- mod.casecohort.y1$coefficients
beta.y2.hat.casecohort <- mod.casecohort.y2$coefficients

AR.y1.hat.casecohort <- est.risk.casecohort$AR.y1
PR.y1.hat.casecohort <- est.risk.casecohort$PR.y1

logAR.y1.hat.casecohort <- log(AR.y1.hat.casecohort)
logPR.y1.hat.casecohort <- log(PR.y1.hat.casecohort)

robvar.beta.y1.hat.casecohort   <- robustvariance(est.y1.casecohort$infl1.beta + est.y1.casecohort$infl2.beta) # equal to diag(mod.casecohort.y1$var)
robvar.beta.y2.hat.casecohort   <- robustvariance(est.y2.casecohort$infl1.beta + est.y2.casecohort$infl2.beta) # equal to diag(mod.casecohort.y2$var)
robvar.AR.y1.hat.casecohort     <- robustvariance(est.risk.casecohort$infl1.AR.y1 + est.risk.casecohort$infl2.AR.y1)
robvar.PR.y1.hat.casecohort     <- robustvariance(est.risk.casecohort$infl1.PR.y1 + est.risk.casecohort$infl2.PR.y1)

var.beta.y1.hat.casecohort      <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.y1.casecohort$infl1.beta, 
                                                        infl2 = est.y1.casecohort$infl2.beta, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 
var.beta.y2.hat.casecohort      <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.y2.casecohort$infl1.beta, 
                                                        infl2 = est.y2.casecohort$infl2.beta, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 
var.AR.y1.hat.casecohort        <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.risk.casecohort$infl1.AR.y1, 
                                                        infl2 = est.risk.casecohort$infl2.AR.y1, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 
var.PR.y1.hat.casecohort        <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.risk.casecohort$infl1.PR.y1, 
                                                        infl2 = est.risk.casecohort$infl2.PR.y1, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 

var.logAR.y1.hat.casecohort <- var.AR.y1.hat.casecohort / (AR.y1.hat.casecohort)^2
var.logPR.y1.hat.casecohort <- var.PR.y1.hat.casecohort / (PR.y1.hat.casecohort)^2

robvar.logAR.y1.hat.casecohort <- robvar.AR.y1.hat.casecohort / (AR.y1.hat.casecohort)^2
robvar.logPR.y1.hat.casecohort <- robvar.PR.y1.hat.casecohort / (PR.y1.hat.casecohort)^2

# print the results in a table

AR.y1.hat.casecohort <- round(AR.y1.hat.casecohort, digits = 3)
PR.y1.hat.casecohort <- round(PR.y1.hat.casecohort, digits = 3)

robvar.AR.y1.hat.casecohort <- round(robvar.AR.y1.hat.casecohort, digits = 7)
robvar.PR.y1.hat.casecohort <- round(robvar.PR.y1.hat.casecohort, digits = 7)

var.AR.y1.hat.casecohort <- round(var.AR.y1.hat.casecohort, digits = 7)
var.PR.y1.hat.casecohort <- round(var.PR.y1.hat.casecohort, digits = 7)

AR.est <- c(AR.y1.hat.casecohort, var.AR.y1.hat.casecohort,
            robvar.AR.y1.hat.casecohort)

PR.est <- c(PR.y1.hat.casecohort, var.PR.y1.hat.casecohort,
            robvar.PR.y1.hat.casecohort)

# ------------------------------------------------------------------------------
# Estimation using the case-cohort with design weights and when the numerator of
# the Breslow estimator is weighted --------------------------------------------

est.y1.casecohort2 <- influences(mod = mod.casecohort.y1, x = x.y1, t1 = t1, 
                                 t2 = t2)

est.y2.casecohort2 <- influences(mod = mod.casecohort.y2, x = x.y2, t1 = t1, 
                                 t2 = t2, competing = TRUE, 
                                 event.times.y1 = est.y1.casecohort2$event.times)

est.risk.casecohort2 <- risk.influences(est.y1.casecohort2, est.y2.casecohort2)

AR.y1.hat.casecohort2 <- est.risk.casecohort2$AR.y1
PR.y1.hat.casecohort2 <- est.risk.casecohort2$PR.y1

robvar.beta.y1.hat.casecohort2   <- robustvariance(est.y1.casecohort2$infl.beta) # equal to diag(mod.casecohort.y1$var)
robvar.beta.y2.hat.casecohort2   <- robustvariance(est.y2.casecohort2$infl.beta) # equal to diag(mod.casecohort.y2$var)
robvar.AR.y1.hat.casecohort2     <- robustvariance(est.risk.casecohort2$infl.AR.y1)
robvar.PR.y1.hat.casecohort2     <- robustvariance(est.risk.casecohort2$infl.PR.y1)

var.beta.y1.hat.casecohort2      <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.y1.casecohort2$infl.beta, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 
var.beta.y2.hat.casecohort2      <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.y2.casecohort2$infl.beta, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 
var.AR.y1.hat.casecohort2        <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.risk.casecohort2$infl.AR.y1, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 
var.PR.y1.hat.casecohort2        <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.risk.casecohort2$infl.PR.y1, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 

# ------------------------------------------------------------------------------
# Print the results in a TeX table ---------------------------------------------

AR.y1.hat.casecohort2 <- round(AR.y1.hat.casecohort2, digits = 3)
PR.y1.hat.casecohort2 <- round(PR.y1.hat.casecohort2, digits = 3)

robvar.AR.y1.hat.casecohort2 <- round(robvar.AR.y1.hat.casecohort2, digits = 7)
robvar.PR.y1.hat.casecohort2 <- round(robvar.PR.y1.hat.casecohort2, digits = 7)

var.AR.y1.hat.casecohort2 <- round(var.AR.y1.hat.casecohort2, digits = 7)
var.PR.y1.hat.casecohort2 <- round(var.PR.y1.hat.casecohort2, digits = 7)

AR.est2 <- c(AR.y1.hat.casecohort2, 
             var.AR.y1.hat.casecohort2,robvar.AR.y1.hat.casecohort2)

PR.est2 <- c(PR.y1.hat.casecohort2, 
             var.PR.y1.hat.casecohort2, robvar.PR.y1.hat.casecohort2)

print(xtable(cbind(rbind(AR.est, PR.est), rbind(AR.est2, PR.est2)), 
             digits = -2), include.rownames = FALSE)



# ------------------------------------------------------------------------------
# Using another covariate profile for the absolute and pure risks --------------
x.y1 <- c(0, 0, 5)
x.y2 <- c(0, 1, 32, 0, 1, 0, 1, 0, 1, 0, 60)

names.x.y1  <- c("x66","x7","x8")
names.x.y2  <- c("x1","x2","x5","x9","x10","x11","x12","x13","x14","x16","x18")

# --------------------------------------------------------------------------
# Estimation using the case-cohort with design weights, when the numerator of 
# the Breslow estimator is not weighted ----------------------------------------

mod.casecohort.y1  <- coxph(Surv(yearsonstudy, status1) ~ x66 +  x7 + 
                              x8, data = casecohort, robust = TRUE, 
                            id = id, weights = weights) # treating event 2 as independent censoring

mod.casecohort.y2  <- coxph(Surv(yearsonstudy, status2) ~ x1 + x2 + x5 + x9 + 
                              x10 + x11 + x12 + x13 + x14 + x16 + x18, 
                            data = casecohort, robust = TRUE, id = id,
                            weights = weights) # treating event 1 as independent censoring

est.y1.casecohort <- influences.generalized(mod = mod.casecohort.y1, x = x.y1, 
                                            t1 = t1, t2 = t2, 
                                            Breslow.all = TRUE,
                                            riskmat.phase1 = riskmat.phase1.y1, 
                                            dNt.phase1 = dNt.phase1.y1)

est.y2.casecohort <- influences.generalized(mod = mod.casecohort.y2, x = x.y2, 
                                            t1 = t1, t2 = t2, competing = TRUE, 
                                            event.times.y1 = est.y1.casecohort$event.times,
                                            riskmat.phase1 = riskmat.phase1.y2, 
                                            dNt.phase1 = dNt.phase1.y2, 
                                            Breslow.all = TRUE)

est.risk.casecohort <- risk.influences.generalized(est.y1.casecohort, 
                                                   est.y2.casecohort)

beta.y1.hat.casecohort <- mod.casecohort.y1$coefficients
beta.y2.hat.casecohort <- mod.casecohort.y2$coefficients

AR.y1.hat.casecohort <- est.risk.casecohort$AR.y1
PR.y1.hat.casecohort <- est.risk.casecohort$PR.y1

logAR.y1.hat.casecohort <- log(AR.y1.hat.casecohort)
logPR.y1.hat.casecohort <- log(PR.y1.hat.casecohort)

robvar.beta.y1.hat.casecohort   <- robustvariance(est.y1.casecohort$infl1.beta + est.y1.casecohort$infl2.beta) # equal to diag(mod.casecohort.y1$var)
robvar.beta.y2.hat.casecohort   <- robustvariance(est.y2.casecohort$infl1.beta + est.y2.casecohort$infl2.beta) # equal to diag(mod.casecohort.y2$var)
robvar.AR.y1.hat.casecohort     <- robustvariance(est.risk.casecohort$infl1.AR.y1 + est.risk.casecohort$infl2.AR.y1)
robvar.PR.y1.hat.casecohort     <- robustvariance(est.risk.casecohort$infl1.PR.y1 + est.risk.casecohort$infl2.PR.y1)

var.beta.y1.hat.casecohort      <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.y1.casecohort$infl1.beta, 
                                                        infl2 = est.y1.casecohort$infl2.beta, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 
var.beta.y2.hat.casecohort      <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.y2.casecohort$infl1.beta, 
                                                        infl2 = est.y2.casecohort$infl2.beta, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 
var.AR.y1.hat.casecohort        <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.risk.casecohort$infl1.AR.y1, 
                                                        infl2 = est.risk.casecohort$infl2.AR.y1, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 
var.PR.y1.hat.casecohort        <- variance.generalized(casecohort = casecohort, 
                                                        infl1 = est.risk.casecohort$infl1.PR.y1, 
                                                        infl2 = est.risk.casecohort$infl2.PR.y1, 
                                                        cohort = cohort, 
                                                        Breslow.all = TRUE) 

var.logAR.y1.hat.casecohort <- var.AR.y1.hat.casecohort / (AR.y1.hat.casecohort)^2
var.logPR.y1.hat.casecohort <- var.PR.y1.hat.casecohort / (PR.y1.hat.casecohort)^2

robvar.logAR.y1.hat.casecohort <- robvar.AR.y1.hat.casecohort / (AR.y1.hat.casecohort)^2
robvar.logPR.y1.hat.casecohort <- robvar.PR.y1.hat.casecohort / (PR.y1.hat.casecohort)^2

# print the results in a table

AR.y1.hat.casecohort <- round(AR.y1.hat.casecohort, digits = 4)
PR.y1.hat.casecohort <- round(PR.y1.hat.casecohort, digits = 4)

robvar.AR.y1.hat.casecohort <- round(robvar.AR.y1.hat.casecohort, digits = 8)
robvar.PR.y1.hat.casecohort <- round(robvar.PR.y1.hat.casecohort, digits = 8)

var.AR.y1.hat.casecohort <- round(var.AR.y1.hat.casecohort, digits = 8)
var.PR.y1.hat.casecohort <- round(var.PR.y1.hat.casecohort, digits = 8)

AR.est <- c(AR.y1.hat.casecohort, var.AR.y1.hat.casecohort,
            robvar.AR.y1.hat.casecohort)

PR.est <- c(PR.y1.hat.casecohort, var.PR.y1.hat.casecohort,
            robvar.PR.y1.hat.casecohort)

# ------------------------------------------------------------------------------
# Estimation using the case-cohort with design weights and when the numerator of
# the Breslow estimator is weighted --------------------------------------------

est.y1.casecohort2 <- influences(mod = mod.casecohort.y1, x = x.y1, t1 = t1, 
                                 t2 = t2)

est.y2.casecohort2 <- influences(mod = mod.casecohort.y2, x = x.y2, t1 = t1, 
                                 t2 = t2, competing = TRUE, 
                                 event.times.y1 = est.y1.casecohort2$event.times)

est.risk.casecohort2 <- risk.influences(est.y1.casecohort2, est.y2.casecohort2)

AR.y1.hat.casecohort2 <- est.risk.casecohort2$AR.y1
PR.y1.hat.casecohort2 <- est.risk.casecohort2$PR.y1

robvar.beta.y1.hat.casecohort2   <- robustvariance(est.y1.casecohort2$infl.beta) # equal to diag(mod.casecohort.y1$var)
robvar.beta.y2.hat.casecohort2   <- robustvariance(est.y2.casecohort2$infl.beta) # equal to diag(mod.casecohort.y2$var)
robvar.AR.y1.hat.casecohort2     <- robustvariance(est.risk.casecohort2$infl.AR.y1)
robvar.PR.y1.hat.casecohort2     <- robustvariance(est.risk.casecohort2$infl.PR.y1)

var.beta.y1.hat.casecohort2      <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.y1.casecohort2$infl.beta, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 
var.beta.y2.hat.casecohort2      <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.y2.casecohort2$infl.beta, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 
var.AR.y1.hat.casecohort2        <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.risk.casecohort2$infl.AR.y1, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 
var.PR.y1.hat.casecohort2        <- variance.generalized(casecohort = casecohort, 
                                                         infl = est.risk.casecohort2$infl.PR.y1, 
                                                         cohort = cohort, 
                                                         Breslow.all = FALSE) 

# ------------------------------------------------------------------------------
# Print the results in a TeX table ---------------------------------------------

AR.y1.hat.casecohort2 <- round(AR.y1.hat.casecohort2, digits = 4)
PR.y1.hat.casecohort2 <- round(PR.y1.hat.casecohort2, digits = 4)

robvar.AR.y1.hat.casecohort2 <- round(robvar.AR.y1.hat.casecohort2, digits = 8)
robvar.PR.y1.hat.casecohort2 <- round(robvar.PR.y1.hat.casecohort2, digits = 8)

var.AR.y1.hat.casecohort2 <- round(var.AR.y1.hat.casecohort2, digits = 8)
var.PR.y1.hat.casecohort2 <- round(var.PR.y1.hat.casecohort2, digits = 8)

AR.est2 <- c(AR.y1.hat.casecohort2, 
             var.AR.y1.hat.casecohort2,robvar.AR.y1.hat.casecohort2)

PR.est2 <- c(PR.y1.hat.casecohort2, 
             var.PR.y1.hat.casecohort2, robvar.PR.y1.hat.casecohort2)

print(xtable(cbind(rbind(AR.est, PR.est), rbind(AR.est2, PR.est2)), 
             digits = -2), include.rownames = FALSE)


