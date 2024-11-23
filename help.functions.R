## -----------------------------------------------------------------------------
## Function: influences()
## -----------------------------------------------------------------------------
## Description: This function estimates the influences on the cause-specific 
##              log-relative hazard, cause-specific baseline hazard at each 
##              unique event time, as well as cause-specific cumulative baseline
##              hazard until each event time of the event of primary interest in
##              interval (t1,t2]. This function should be used with design
##              weights under the case-cohort design with exhaustive sampling of 
##              cases, or under the case-cohort design with stratified sampling 
##              based on case status when using only the event times of the
##              cases in the case-cohort (Breslow.weight). If the event type is 
##              not that of primary interest, interval (t1,t2] and the event 
##              times of the primary event of interest need to be provided. Also 
##              returns parameters estimates
## -----------------------------------------------------------------------------
## Arguments:
##
##  mod             a cox model object, result of function coxph 
##
##  x               vector of length p, specifying the covariate profile to be
##                  considered for the pure risk. Default is (0,...,0)
##
##  competing       is the event type that of primary interest or a competing 
##                  event. If competing = TRUE, the arguments below need to be
##                  provided. Default is FALSE
##
##  t1              left bound of the time interval considered for the 
##                  cumulative baseline hazard and pure risk. If 
##                  competing = FALSE, the default is the first event time. If 
##                  competing = TRUE, this argument needs to be provided
##                  
##  t2              right bound of the time interval considered for the 
##                  cumulative baseline hazard and pure risk. If 
##                  competing = FALSE, the default is the last event time. If 
##                  competing = TRUE, this argument needs to be provided
##
##  event.times.y1  vector with the event times for the primary event of  
##                  interest (y1). If competing = TRUE, this argument needs to
##                  be provided
## -----------------------------------------------------------------------------

influences  <- function (mod, competing = NULL, event.times.y1 = NULL, 
                         t1 = NULL, t2 = NULL, x = NULL) {
  
  if (is.null(competing)) {
    competing <- FALSE
  }
  
  # ----------------------------------------------------------------------------
  # Quantities needed for the influences ---------------------------------------
  
  mod.detail        <- coxph.detail(mod, riskmat = T)
  riskmat           <- mod.detail$riskmat 
  number.times      <- ncol(riskmat) 
  n                 <- nrow(riskmat) 
  observed.times    <- apply(riskmat, 1, function(v) {which.max(cumsum(v))}) 
  X                 <- model.matrix(mod) 
  
  weights           <- mod$weights # weights used for the fit, w_i
  if (is.null(weights)) { weights <- 1 } # when using the whole cohort, w_i = 1
  
  beta.hat          <- mod$coefficients
  
  p                 <- length(beta.hat)
  if ((length(x) != p) | (is.null(x))) {x <- rep(0, p)} # if no covariate profile provided, use 0 as reference level
  exp.x             <- c(exp(x %*% beta.hat))
  
  exp.X.weighted    <- weights * exp(X %*% beta.hat)
  Y.exp.X.weighted  <- riskmat * matrix(exp.X.weighted, nrow = n, 
                                        ncol = number.times, byrow = FALSE) 
  S0t               <- t(riskmat) %*% (exp.X.weighted) 
  S1t               <- t(riskmat) %*% (X * matrix(exp.X.weighted, nrow = n,
                                                  ncol = p, byrow = FALSE)) 
  dNt               <- matrix(0, n, number.times) 
  dNt[cbind(1:n, observed.times)] <- 1 
  dNt1              <- dNt
  dNt               <- dNt * matrix(mod$y[, ncol(mod$y)], nrow(riskmat), 
                                    number.times, byrow = FALSE) 
  dNt.weighted      <- dNt * matrix(weights, nrow = n, ncol = number.times, byrow = FALSE) 
  
  # ----------------------------------------------------------------------------
  # Estimation of the cause-specific baseline hazard at each unique event time -
  
  lambda0.t.hat         <- t(colSums(dNt.weighted) / S0t) 
  
  # ----------------------------------------------------------------------------
  # Estimation of the cause-specific "survival function" at each unique event 
  # time -----------------------------------------------------------------------
  
  if(competing == FALSE){
    if (is.null(t1)) { t1 <- floor(min(mod.detail$time)) }
    if (is.null(t2)) { t2 <- floor(max(mod.detail$time)) }
    t1t2.times        <- which((t1 < mod.detail$time) & (mod.detail$time <= t2)) 
    
    lambda0.t12.hat   <- lambda0.t.hat[t1t2.times] # baseline hazard for each of the event times in (t1, t2]
    Lambda0.t12.hat   <- cumsum(lambda0.t12.hat) # cumulative baseline hazard in time interval (t1, t] for each of the event times t in (t1, t2]
    S.t12.hat         <- exp(- Lambda0.t12.hat * exp.x) # vector with the S(t) / S(t1) for each of the event times t in (t1, t2]
  
  }else{
    if (is.null(t1)) { stop("t1 should be the same as for the event of primary 
                            interest") }
    if (is.null(t2)) { stop("t2 should be the same as for the event of primary 
                            interest") }
    if (is.null(event.times.y1)) { stop("the event times for the event of 
                                        primary interest should be provided") }
    
    Lambda0.t12.hat   <- sapply(event.times.y1, 
                                function(t) {sum(lambda0.t.hat[which((t1 < mod.detail$time) & 
                                                                       (mod.detail$time <= t))])}) # cumulative baseline hazard in time interval (t1, t] for each of the event times t of the primary event of interest in (t1, t2]
    S.t12.hat         <- exp(- Lambda0.t12.hat * exp.x)
  }
  
  # --------------------------------------------------------------------------
  # Computation of the influences on the cause-specific log-relative hazard --
  
  infl.beta <- residuals(mod, type = "dfbeta", weighted = T)
  
  # --------------------------------------------------------------------------
  # Computation of the influences on the cause-specific the baseline hazard at
  # each unique event time ---------------------------------------------------
  
  infl.lambda0.t    <- (dNt.weighted - (Y.exp.X.weighted + infl.beta %*% 
                                          t(S1t)) * matrix(lambda0.t.hat, 
                                                           nrow = n, 
                                                           ncol = number.times,
                                                           byrow = TRUE)) /
    matrix(S0t, nrow = n, ncol = number.times, byrow = TRUE)
  
  if(competing == FALSE){
    infl.lambda0.t12    <- infl.lambda0.t[,t1t2.times]
    infl.Lambda0.t12    <- t(apply(infl.lambda0.t12, 1, cumsum))
    
    return(list(beta.hat = beta.hat, exp.x = exp.x, 
                lambda0.t.hat = lambda0.t.hat, competing = competing,
                event.times.y1 = mod.detail$time[t1t2.times], t1 = t1, 
                t2 = t2, x = x, lambda0.t12.hat = lambda0.t12.hat, 
                Lambda0.t12.hat = Lambda0.t12.hat, S.t12.hat = S.t12.hat,
                infl.beta = infl.beta, infl.lambda0.t = infl.lambda0.t, 
                infl.lambda0.t12 = infl.lambda0.t12,
                infl.Lambda0.t12 = infl.Lambda0.t12))
  }else{
    infl.Lambda0.t12   <- sapply(event.times.y1, 
                                 function(t) {rowSums(as.matrix(infl.lambda0.t[,which((t1 < mod.detail$time) & 
                                                                                        (mod.detail$time <= t))]))})
    return(list(beta.hat = beta.hat, exp.x = exp.x, 
                lambda0.t.hat = lambda0.t.hat, competing = competing, 
                event.times.y1 = event.times.y1, t1 = t1, t2 = t2, x = x, 
                Lambda0.t12.hat = Lambda0.t12.hat, S.t12.hat = S.t12.hat, 
                infl.beta = infl.beta, infl.lambda0.t = infl.lambda0.t,
                infl.Lambda0.t12 = infl.Lambda0.t12))
  }
}



## -----------------------------------------------------------------------------
## Function: influences.generalized()
## -----------------------------------------------------------------------------
## Description: This function estimates the influences on the cause-specific 
##              log-relative hazard, cause-specific baseline hazard at each 
##              unique event time, as well as cause-specific cumulative baseline
##              hazard until each event time of the event of primary interest in
##              interval (t1,t2]. This function should be used with design
##              weights under the case-cohort design with stratified sampling 
##              based on case status when using the event event times of all the
##              cases in the cohort (Breslow.all). If the event type is not that 
##              of primary interest, interval (t1,t2] and the event times of the 
##              primary event of interest need to be provided. Because the event 
##              times of all the cases in the cohort are used in the Breslow 
##              estimator, even if only a fraction of the cases are included
##              in the case-cohort, the phase-one influences and the phase-two 
##              influences are returned. Also returns parameters estimates
## -----------------------------------------------------------------------------
## Arguments:
##
##  mod             a cox model object, result of function coxph 
##
##  x               vector of length p, specifying the covariate profile to be 
##                  considered for the pure risk. Default is (0,...,0)
##
##  competing       is the event type that of primary interest or a competing 
##                  event. If competing = TRUE, arguments t1, t2, and
##                  event.times.y1 below need to be provided. Default is FALSE
##
##  t1              left bound of the time interval considered for the 
##                  cumulative baseline hazard and pure risk. If 
##                  competing = FALSE, the default is the first event time. If 
##                  competing = TRUE, this argument needs to be provided
##                  
##  t2              right bound of the time interval considered for the 
##                  cumulative baseline hazard and pure risk. If 
##                  competing = FALSE, the default is the last event time. If 
##                  competing = TRUE, this argument needs to be provided
##
##  event.times.y1  vector with the event times for the primary event of  
##                  interest (y1). If competing = TRUE, this argument needs to
##                  be provided
##
##  Breslow.all     was only a fraction of the cases included in the case-cohort 
##                  and were the event times of all the cases in the cohort used
##                  in the numerator of the Breslow estimator? If 
##                  Breslow.all = TRUE, the arguments below need to be provided. 
##                  Default is TRUE. If Breslow.all = FALSE, 
##                  influences.generalized() is equivalent to using influences()
##
##  riskmat.phase1  at risk matrix for the cohort (phase-one data) at all the  
##                  cases event times, even that not in phase-two. If 
##                  Breslow.all = TRUE, this argument needs to be provided
##
##  status.phase1   vector indicating the case status in the case-cohort (phase- 
##                  one data). If Breslow.all = TRUE and dNt.phase1 = NULL, this 
##                  argument needs to be provided
##
##  dNt.phase1      counting process matrix for failures in the cohort (phase-
##                  one data). If Breslow.all = TRUE and status.phase1 = NULL, 
##                  this argument needs to be provided
## -----------------------------------------------------------------------------

influences.generalized <- function (mod, competing = NULL, 
                                    event.times.y1 = NULL, Breslow.all = NULL, 
                                    riskmat.phase1 = NULL, dNt.phase1 = NULL, 
                                    status.phase1 = NULL, t1 = NULL, t2 = NULL, 
                                    x = NULL) {
  
  if (is.null(competing)) {
    competing <- FALSE
  }
  
  if (is.null(Breslow.all)) {
    Breslow.all <- TRUE
  }
  
  # If only a fraction of the cases is used ------------------------------------
  
  if (Breslow.all == TRUE) {
    if (is.null(riskmat.phase1)) {
      stop("The at risk matrix for the phase-one data at all of the cases event 
           times need to be provided")
    } else {
      if (is.null(dNt.phase1)&is.null(status.phase1)) {
        stop("The status for the phase-one data, or the counting process matrix 
             for failures in the phase-one data, need to be provided")
      } else {
        if (is.null(dNt.phase1)) { 
          
          # If not provided, compute the counting process matrix for failures in 
          # the phase-one data -------------------------------------------------
          
          observed.times.phase1 <- apply(riskmat.phase1, 1,
                                         function(v) {which.max(cumsum(v))})
          dNt.phase1            <- matrix(0, nrow(riskmat.phase1), 
                                          ncol(riskmat.phase1))
          dNt.phase1[cbind(1:nrow(riskmat.phase1), observed.times.phase1)] <- 1
          dNt.phase1            <- sweep(dNt.phase1, 1, status.phase1, "*") 
        }
        
        # ----------------------------------------------------------------------
        # Quantities needed for the influences ---------------------------------
        
        mod.detail            <- coxph.detail(mod, riskmat = TRUE)
        X                     <- model.matrix(mod) 
        riskmat               <- mod.detail$riskmat
        n.phase2              <- nrow(X)
        n                     <- nrow(riskmat.phase1)
        number.times.phase1   <- ncol(riskmat.phase1)
        indiv.phase2          <- row.names(X)
        riskmat.phase2        <- riskmat.phase1[indiv.phase2, ] # use the cases' 
        # actual failure time, known even for cases not in phase-two
        cases.times           <- as.numeric(colnames(riskmat.phase1)) # cases failure times, even 
        # the one who are not in phase-two
        
        weights           <- mod$weights # weights used for the fit, w_i
        if (is.null(weights)) { weights <- 1 } # when using the whole cohort, w_i = 1
        
        beta.hat          <- mod$coefficients
        
        p                 <- length(beta.hat)
        if ((length(x) != p) | (is.null(x))) {x <- rep(0, p)} # if no covariate profile provided, use 0 as reference level
        exp.x             <- c(exp(x %*% beta.hat))
        
        exp.X.weighted    <- weights * exp(X %*% beta.hat)
  
        Y.exp.X.weighted.casestimes   <- sweep(riskmat.phase2, 1, 
                                               exp.X.weighted, "*")
        
        S0t.casestimes  <- t(riskmat.phase2) %*% (exp.X.weighted)  # for any
        # event time t (even the one of cases not in phase-two)
        S1t.casestimes  <- t(riskmat.phase2) %*% (sweep(X, 1, 
                                                        exp.X.weighted, "*"))
  
        # ----------------------------------------------------------------------
        # Estimation of the cause-specific baseline hazard at each unique event 
        # time -----------------------------------------------------------------
        
        lambda0.t.hat         <- t(colSums(dNt.phase1) / S0t.casestimes) 
        
        # ----------------------------------------------------------------------
        # Estimation of the cause-specific "survival function" at each unique 
        # event time -----------------------------------------------------------
        
        if(competing == FALSE){
          if (is.null(t1)) { t1 <- floor(min(mod.detail$time)) }
          if (is.null(t2)) { t2 <- floor(max(mod.detail$time)) }
          t1t2.times        <- which((t1 < cases.times) & (cases.times <= t2)) 
          
          lambda0.t12.hat   <- lambda0.t.hat[t1t2.times] # baseline hazard for each of the event times in (t1, t2]
          Lambda0.t12.hat   <- cumsum(lambda0.t12.hat) # cumulative baseline hazard in time interval (t1, t] for each of the event times in (t1, t2]
          S.t12.hat         <- exp(- Lambda0.t12.hat * exp.x) # vector with the S(t) / S(t1) for t in (t1, t2]
          
        }else{
          if (is.null(t1)) { stop("t1 should be the same as for the event of primary 
                            interest") }
          if (is.null(t2)) { stop("t2 should be the same as for the event of primary 
                            interest") }
          if (is.null(event.times.y1)) { stop("the event times for the event of 
                                        primary interest should be provided") }
          
          Lambda0.t12.hat   <- sapply(event.times.y1, 
                                      function(t) {sum(lambda0.t.hat[which((t1 < cases.times) & 
                                                                                             (cases.times <= t))])}) # cumulative baseline hazard in time interval (t1, t] for each of the event times t of the primary event of interest in (t1, t2]
          S.t12.hat         <- exp(- Lambda0.t12.hat * exp.x)
        }
        
        infl1.beta      <- matrix(0, nrow = nrow(riskmat.phase1), ncol = p)
        rownames(infl1.beta) <- rownames(riskmat.phase1)
        infl2.beta      <- infl1.beta
        infl2.beta[indiv.phase2,] <- residuals(mod, type = "dfbeta", weighted = T)
        
        infl1.lambda0.t     <- dNt.phase1 / matrix(S0t.casestimes, nrow = n, 
                                                  ncol = number.times.phase1, 
                                                  byrow = TRUE)
        infl2.lambda0.t     <- matrix(0, nrow = nrow(infl1.lambda0.t), 
                                      ncol = ncol(infl1.lambda0.t))
        rownames(infl2.lambda0.t) <- rownames(riskmat.phase1)
        infl2.lambda0.t[indiv.phase2,]    <- (- (Y.exp.X.weighted.casestimes + 
                                                   infl2.beta[indiv.phase2,] %*% 
                                                   t(S1t.casestimes)) * 
          matrix(lambda0.t.hat, nrow = n.phase2, ncol = number.times.phase1,
                     byrow = TRUE)) /
          matrix(S0t.casestimes, nrow = n.phase2, ncol = number.times.phase1, 
                 byrow = TRUE)
        
        if(competing == FALSE){
          
          infl1.lambda0.t12     <- infl1.lambda0.t[,t1t2.times]
          infl1.Lambda0.t12     <- t(apply(infl1.lambda0.t12, 1, cumsum))
          
          infl2.lambda0.t12     <- infl2.lambda0.t[,t1t2.times]
          infl2.Lambda0.t12     <- t(apply(infl2.lambda0.t12, 1, cumsum))
          
          return(list(beta.hat = beta.hat, exp.x = exp.x, 
                      lambda0.t.hat = lambda0.t.hat, competing = competing,
                      event.times.y1 = cases.times[t1t2.times], t1 = t1, 
                      t2 = t2, x = x, lambda0.t12.hat = lambda0.t12.hat, 
                      Lambda0.t12.hat = Lambda0.t12.hat, S.t12.hat = S.t12.hat,
                      infl1.beta = infl1.beta, infl2.beta = infl2.beta, 
                      infl1.lambda0.t = infl1.lambda0.t,
                      infl2.lambda0.t = infl2.lambda0.t,
                      infl1.lambda0.t12 = infl1.lambda0.t12,
                      infl2.lambda0.t12 = infl2.lambda0.t12,
                      infl1.Lambda0.t12 = infl1.Lambda0.t12,
                      infl2.Lambda0.t12 = infl2.Lambda0.t12,
                      Breslow.all = Breslow.all))
        }else{
          infl1.Lambda0.t12   <- sapply(event.times.y1, 
                                        function(t) {rowSums(as.matrix(infl1.lambda0.t[,which((t1 < cases.times) & 
                                                                                                                (cases.times <= t))]))})
          infl2.Lambda0.t12   <- sapply(event.times.y1, 
                                        function(t) {rowSums(as.matrix(infl2.lambda0.t[,which((t1 < cases.times) & 
                                                                                                                (cases.times <= t))]))})
          return(list(beta.hat = beta.hat, exp.x = exp.x, 
                      lambda0.t.hat = lambda0.t.hat, competing = competing, 
                      event.times.y1 = event.times.y1, t1 = t1, t2 = t2, x = x, 
                      Lambda0.t12.hat = Lambda0.t12.hat, S.t12.hat = S.t12.hat, 
                      infl1.beta = infl1.beta, infl2.beta = infl2.beta,
                      infl1.lambda0.t = infl1.lambda0.t,
                      infl2.lambda0.t = infl2.lambda0.t,
                      infl1.Lambda0.t12 = infl1.Lambda0.t12,
                      infl2.Lambda0.t12 = infl2.Lambda0.t12,
                      Breslow.all = Breslow.all))
        }
      }
    }
  }else{
    
    # --------------------------------------------------------------------------
    # Quantities needed for the influences -------------------------------------
    
    mod.detail        <- coxph.detail(mod, riskmat = T)
    riskmat           <- mod.detail$riskmat 
    number.times      <- ncol(riskmat) 
    n                 <- nrow(riskmat) 
    observed.times    <- apply(riskmat, 1, function(v) {which.max(cumsum(v))}) 
    X                 <- model.matrix(mod) 

    weights           <- mod$weights # weights used for the fit, w_i
    if (is.null(weights)) { weights <- 1 } # when using the whole cohort, w_i = 1
    
    beta.hat          <- mod$coefficients
    
    p                 <- length(beta.hat)
    if ((length(x) != p) | (is.null(x))) {x <- rep(0, p)} # if no covariate profile provided, use 0 as reference level
    exp.x             <- c(exp(x %*% beta.hat))
    
    exp.X.weighted    <- weights * exp(X %*% beta.hat)
    Y.exp.X.weighted  <- riskmat * matrix(exp.X.weighted, nrow = n, 
                                          ncol = number.times, byrow = FALSE) 
    S0t               <- t(riskmat) %*% (exp.X.weighted) 
    S1t               <- t(riskmat) %*% (X * matrix(exp.X.weighted, nrow = n,
                                                    ncol = p, byrow = FALSE)) 
    dNt               <- matrix(0, n, number.times) 
    dNt[cbind(1:n, observed.times)] <- 1 
    dNt1 = dNt
    dNt               <- dNt * matrix(mod$y[, ncol(mod$y)], nrow(riskmat), 
                                      number.times, byrow = FALSE) 
    dNt.weighted      <- dNt * matrix(weights, nrow = n, ncol = number.times, 
                                      byrow = FALSE) 
    
    # --------------------------------------------------------------------------
    # Estimation of the cause-specific baseline hazard at each unique event time

    lambda0.t.hat         <- t(colSums(dNt.weighted) / S0t) 
    
    # --------------------------------------------------------------------------
    # Estimation of the cause-specific "survival function" at each unique event 
    # time ---------------------------------------------------------------------
    
    if(competing == FALSE){
      if (is.null(t1)) { t1 <- floor(min(mod.detail$time)) }
      if (is.null(t2)) { t2 <- floor(max(mod.detail$time)) }
      t1t2.times        <- which((t1 < mod.detail$time) & (mod.detail$time <= t2)) 
      
      lambda0.t12.hat   <- lambda0.t.hat[t1t2.times] # baseline hazard for each of the event times in (t1, t2]
      Lambda0.t12.hat   <- cumsum(lambda0.t12.hat) # cumulative baseline hazard in time interval (t1, t] for each of the event times t in (t1, t2]
      S.t12.hat         <- exp(- Lambda0.t12.hat * exp.x) # vector with the S(t) / S(t1) for each of the event times t in (t1, t2]
    }else{
      if (is.null(t1)) { stop("t1 should be the same as for the event of primary 
                            interest") }
      if (is.null(t2)) { stop("t2 should be the same as for the event of primary 
                            interest") }
      if (is.null(event.times.y1)) { stop("the event times for the event of 
                                        primary interest should be provided") }
      
      Lambda0.t12.hat   <- sapply(event.times.y1, 
                                  function(t) {sum(lambda0.t.hat[which((t1 < mod.detail$time) & 
                                                                         (mod.detail$time <= t))])}) # cumulative baseline hazard in time interval (t1, t] for each of the event times t of the primary event of interest in (t1, t2]
      S.t12.hat         <- exp(- Lambda0.t12.hat * exp.x)
    }
    
    # --------------------------------------------------------------------------
    # Computation of the influences on the cause-specific log-relative hazard --
    
    infl.beta <- residuals(mod, type = "dfbeta", weighted = T)
    
    # --------------------------------------------------------------------------
    # Computation of the influences on the cause-specific the baseline hazard at
    # each unique event time ---------------------------------------------------
    
    infl.lambda0.t    <- (dNt.weighted - (Y.exp.X.weighted + infl.beta %*% t(S1t)) * 
                            matrix(lambda0.t.hat, nrow = n, ncol = number.times,
                                   byrow = TRUE)) /
      matrix(S0t, nrow = n, ncol = number.times, byrow = TRUE)
    
    if(competing == FALSE){
      infl.lambda0.t12    <- infl.lambda0.t[,t1t2.times]
      infl.Lambda0.t12    <- t(apply(infl.lambda0.t12, 1, cumsum))
      
      return(list(beta.hat = beta.hat, exp.x = exp.x, 
                  lambda0.t.hat = lambda0.t.hat, competing = competing,
                  event.times.y1 = mod.detail$time[t1t2.times], t1 = t1, 
                  t2 = t2, x = x, lambda0.t12.hat = lambda0.t12.hat, 
                  Lambda0.t12.hat = Lambda0.t12.hat, S.t12.hat = S.t12.hat,
                  infl.beta = infl.beta, infl.lambda0.t = infl.lambda0.t, 
                  infl.lambda0.t12 = infl.lambda0.t12,
                  infl.Lambda0.t12 = infl.Lambda0.t12,
                  Breslow.all = Breslow.all))
    }else{
      infl.Lambda0.t12   <- sapply(event.times.y1, 
                                   function(t) {rowSums(as.matrix(infl.lambda0.t[,which((t1 < mod.detail$time) & 
                                                                                                          (mod.detail$time <= t))]))})
      return(list(beta.hat = beta.hat, exp.x = exp.x, 
                  lambda0.t.hat = lambda0.t.hat, competing = competing, 
                  event.times.y1 = event.times.y1, t1 = t1, t2 = t2, x = x, 
                  Lambda0.t12.hat = Lambda0.t12.hat, S.t12.hat = S.t12.hat, 
                  infl.beta = infl.beta, infl.lambda0.t = infl.lambda0.t,
                  infl.Lambda0.t12 = infl.Lambda0.t12,
                  Breslow.all = Breslow.all))
    }
  }
}



## -----------------------------------------------------------------------------
## Function: influences.calib()
## -----------------------------------------------------------------------------
## Description: This function estimates the influences on the cause-specific 
##              log-relative hazard, cause-specific baseline hazard at each 
##              unique event time, as well as cause-specific cumulative baseline
##              hazard until each event time of the event of primary interest in
##              interval (t1,t2]. This function should be used with calibrated
##              weights under the case-cohort design with exhaustive sampling of 
##              cases. If the event type is not that of primary interest,
##              interval (t1,t2] and the event times of the primary event
##              event of interest need to be provided. Because calibrated 
##              weights are used, the phase-one influences and the phase-two 
##              influences. are returned. Also returns parameters estimates
## -----------------------------------------------------------------------------
## Arguments:
##
##  mod             a cox model object, result of function coxph 
##
##  x               vector of length p, specifying the covariate profile 
##                  considered for the pure risk. Default is (0,...,0)
##
##  A               matrix with the values of the auxiliary variables used for 
##                  the calibration of the design weights in the whole cohort
##
##  competing       is the event type that of primary interest or a competing 
##                  event. If competing = TRUE, the arguments below need to be
##                  provided. Default is FALSE
##
##  t1              left bound of the time interval considered for the 
##                  cumulative baseline hazard and pure risk. If 
##                  competing = FALSE, the default is the first event time. If 
##                  competing = TRUE, this argument needs to be provided
##                  
##  t2              right bound of the time interval considered for the 
##                  cumulative baseline hazard and pure risk. If 
##                  competing = FALSE, the default is the last event time. If 
##                  competing = TRUE, this argument needs to be provided
##
##  event.times.y1  vector with the event times for the primary event of  
##                  interest (y1). If competing = TRUE, this argument needs to
##                  be provided
## -----------------------------------------------------------------------------

influences.calib  <- function (mod, A, competing = NULL, event.times.y1 = NULL,
                               t1 = NULL, t2 = NULL, x = NULL) {
  
  if (is.null(competing)) {
    competing <- FALSE
  }
  
  # ----------------------------------------------------------------------------
  # Quantities needed for the influences ---------------------------------------
  
  mod.detail        <- coxph.detail(mod, riskmat = T)
  riskmat           <- mod.detail$riskmat 
  number.times      <- ncol(riskmat) 
  n                 <- nrow(riskmat) 
  X                 <- model.matrix(mod) 
  p                 <- ncol(X) 
  weights           <- mod$weights # weights used for the fit, w_i
  if (is.null(weights)) { weights <- 1 } # when using the whole cohort, w_i = 1
  
  beta.hat          <- mod$coefficients
  
  p                 <- length(beta.hat)
  if ((length(x) != p) | (is.null(x))) {x <- rep(0, p)} # if no covariate profile provided, use 0 as reference level
  exp.x             <- c(exp(x %*% beta.hat))
  
  exp.X.weighted    <- weights * exp(X %*% beta.hat)
  Y.exp.X.weighted  <- riskmat * matrix(exp.X.weighted, nrow = n, 
                                        ncol = number.times, byrow = FALSE) 
  S0t               <- t(riskmat) %*% (exp.X.weighted) 
  S1t               <- t(riskmat) %*% (X * matrix(exp.X.weighted, nrow = n,
                                                  ncol = p, byrow = FALSE)) 
  observed.times    <- apply(riskmat, 1, function(v) {which.max(cumsum(v))}) 
  dNt               <- matrix(0, n, number.times) 
  dNt[cbind(1:n, observed.times)] <- 1 
  dNt               <- dNt * matrix(mod$y[, ncol(mod$y)], nrow(riskmat), 
                                    number.times, byrow = FALSE) 
  dNt.weighted      <- dNt * matrix(weights, nrow = n, ncol = number.times, 
                                    byrow = FALSE) 
  
  infomat.indiv     <- mod.detail$imat  
  infomat           <- apply(X = infomat.indiv, MARGIN = c(1,2), FUN = sum)
  
  # ----------------------------------------------------------------------------
  # Estimation of the cause-specific baseline hazard at each unique event time -
  
  lambda0.t.hat         <- t(colSums(dNt) / S0t) 
  
  # ----------------------------------------------------------------------------
  # Estimation of the cause-specific "survival function" at each unique event 
  # time -----------------------------------------------------------------------
  
  if(competing == FALSE){
    if (is.null(t1)) { t1 <- floor(min(mod.detail$time)) }
    if (is.null(t2)) { t2 <- floor(max(mod.detail$time)) }
    t1t2.times        <- which((t1 < mod.detail$time) & (mod.detail$time <= t2)) 
    
    lambda0.t12.hat   <- lambda0.t.hat[t1t2.times] # baseline hazard for each of the event times in (t1, t2]
    Lambda0.t12.hat   <- cumsum(lambda0.t12.hat) # cumulative baseline hazard in time interval (t1, t] for each of the event times t in (t1, t2]
    S.t12.hat         <- exp(- Lambda0.t12.hat * exp.x) # vector with the S(t) / S(t1) for each of the event times t in (t1, t2]
  }else{
    if (is.null(t1)) { stop("t1 should be the same as for the event of primary 
                            interest") }
    if (is.null(t2)) { stop("t2 should be the same as for the event of primary 
                            interest") }
    if (is.null(event.times.y1)) { stop("the event times for the event of 
                                        primary interest should be provided") }
    
    Lambda0.t12.hat   <- sapply(event.times.y1, 
                                function(t) {sum(lambda0.t.hat[which((t1 < mod.detail$time) & 
                                                                       (mod.detail$time <= t))])}) # cumulative baseline hazard in time interval (t1, t] for each of the event times t of the primary event of interest in (t1, t2]
    S.t12.hat         <- exp(- Lambda0.t12.hat * exp.x)
  }
  
  # ----------------------------------------------------------------------------
  # Computation of the influences for the Lagrangian multipliers, eta ----------
  
  indiv.phase2      <- row.names(X)
  A.phase2          <- A[indiv.phase2, ]
  q                 <- ncol(A.phase2)
  
  A.weighted         <- sweep(A.phase2, 1, weights, "*")
  AA.weighted        <- array(NA, dim = c(q, q, n))
  for(i in 1:nrow(A.phase2)) {
    AA.weighted[,, i] <- tcrossprod(A.phase2[i,], A.weighted[i,])
  }
  sum.AA.weighted    <- apply(X = AA.weighted, MARGIN = c(1,2), FUN = sum)
  sum.AA.weighted.inv <- solve(sum.AA.weighted)
  
  infl1.eta         <- A %*% sum.AA.weighted.inv
  infl2.eta         <- 0 * A
  infl2.eta[indiv.phase2,] <- - A.weighted %*% sum.AA.weighted.inv
  infl.eta          <- infl1.eta + infl2.eta 
  
  # ----------------------------------------------------------------------------
  # Computation of the influences for log-relative hazard, beta ----------------
  
  XA              <- array(NA, dim = c(p, q, n))
  for(i in 1:n) {
    XA[,, i]      <- tcrossprod(X[i,], A.phase2[i,])
  }
  drond.G1t.eta   <- array(NA, dim = c(p, q, number.times))
  drond.G0t.eta   <- array(NA, dim = c(q, number.times))
  drond.S1t.eta   <- array(NA, dim = c(p, q, number.times))
  drond.S0t.eta   <- array(NA, dim = c(q, number.times))
  for(t in 1:number.times) { 
    drond.G1t.eta[,, t] <- apply(X = sweep(XA, 3, dNt.weighted[, t], "*"), 
                                 MARGIN = c(1,2), FUN = sum)
    drond.G0t.eta[, t]  <- colSums(A.phase2 * matrix(dNt.weighted[, t], 
                                                     nrow = n, ncol = q, 
                                                     byrow = FALSE))
    drond.S1t.eta[,, t] <- apply(X = sweep(XA, 3, Y.exp.X.weighted[, t],
                                           "*"), MARGIN = c(1,2), FUN = sum)
    drond.S0t.eta[, t]  <- colSums(A.phase2 * matrix(Y.exp.X.weighted[, t], 
                                                     nrow = n, ncol = q, 
                                                     byrow = FALSE))
  }
  drond.Ut.eta    <- array(NA, dim = c(p, q, number.times))
  X.expect        <- S1t / matrix(S0t, nrow = number.times, ncol = p, 
                                  byrow = FALSE)
  for(t in 1:number.times) { 
    drond.Ut.eta[,, t] <- drond.G1t.eta[,, t] - tcrossprod(X.expect[t, ],
                                                           drond.G0t.eta[, t]) -
      (colSums(dNt.weighted) / S0t)[t] * drond.S1t.eta[,, t] +
      (colSums(dNt.weighted) / S0t ^ 2)[t] * 
      tcrossprod(S1t[t, ], drond.S0t.eta[, t])
  }
  drond.U.eta     <- apply(X = drond.Ut.eta, MARGIN = c(1,2), FUN = sum)
  infl1.beta      <- infl1.eta %*% t(drond.U.eta) %*% solve(infomat)
  score           <- array(NA, dim = c(n, p, number.times))
  for(i in 1:n) {
    score[i,,]    <- - sweep(t(X.expect), 1, X[i,], "-")
  }
  for(t in 1:number.times) {
    dMt           <- dNt.weighted[, t] - riskmat[, t] *
      exp.X.weighted * (colSums(dNt.weighted)/S0t)[t]
    score[,, t]   <- score[,, t] * matrix(dMt, nrow = n, ncol = p, 
                                          byrow = FALSE)
  }
  score.beta      <- apply(X = score, MARGIN = c(1,2), FUN = sum) 
  infl2.beta      <- infl2.eta %*% t(drond.U.eta) %*% solve(infomat)
  infl2.beta[indiv.phase2,] <- infl2.beta[indiv.phase2,] +
    score.beta %*% solve(infomat) 
  
  # --------------------------------------------------------------------------
  # Computation of the influences on the cause-specific the baseline hazard at
  # each unique event time ---------------------------------------------------
  
  infl1.lambda0.t <- (- (infl1.beta %*% t(S1t) + infl1.eta %*% drond.S0t.eta) *
                        matrix(lambda0.t.hat, nrow = nrow(A), 
                               ncol = number.times, byrow = TRUE)) /
    matrix(S0t, nrow = nrow(A), ncol = number.times, 
           byrow = TRUE) 
  infl2.lambda0.t <- 0 * infl1.lambda0.t
  infl2.lambda0.t[indiv.phase2,] <- (dNt - (Y.exp.X.weighted + 
                                              infl2.beta[indiv.phase2,] %*% 
                                              t(S1t) + infl2.eta[indiv.phase2,] %*%
                                              drond.S0t.eta) * 
                                       matrix(lambda0.t.hat, nrow = n, 
                                              ncol = number.times, byrow = TRUE)) /
    matrix(S0t, nrow = n, ncol = number.times,
           byrow = TRUE) 
  
  if(competing == FALSE){
    infl1.lambda0.t12     <- infl1.lambda0.t[,t1t2.times]
    infl1.Lambda0.t12     <- t(apply(infl1.lambda0.t12, 1, cumsum))
    
    infl2.lambda0.t12     <- infl2.lambda0.t[,t1t2.times]
    infl2.Lambda0.t12     <- t(apply(infl2.lambda0.t12, 1, cumsum))
    
    return(list(beta.hat = beta.hat, exp.x = exp.x, 
                lambda0.t.hat = lambda0.t.hat, competing = competing,
                event.times.y1 = mod.detail$time[t1t2.times], t1 = t1, 
                t2 = t2, x = x, lambda0.t12.hat = lambda0.t12.hat, 
                Lambda0.t12.hat = Lambda0.t12.hat, S.t12.hat = S.t12.hat,
                infl1.beta = infl1.beta, infl2.beta = infl2.beta, 
                infl1.lambda0.t = infl1.lambda0.t,
                infl2.lambda0.t = infl2.lambda0.t,
                infl1.lambda0.t12 = infl1.lambda0.t12,
                infl2.lambda0.t12 = infl2.lambda0.t12,
                infl1.Lambda0.t12 = infl1.Lambda0.t12,
                infl2.Lambda0.t12 = infl2.Lambda0.t12))
  }else{
    infl1.Lambda0.t12   <- sapply(event.times.y1, function(t) {rowSums(as.matrix(infl1.lambda0.t[,which((t1 < mod.detail$time) & 
                                                                                                          (mod.detail$time <= t))]))})
    infl2.Lambda0.t12   <- sapply(event.times.y1, function(t) {rowSums(as.matrix(infl2.lambda0.t[,which((t1 < mod.detail$time) & 
                                                                                                          (mod.detail$time <= t))]))})
    return(list(beta.hat = beta.hat, exp.x = exp.x, 
                lambda0.t.hat = lambda0.t.hat, competing = competing, 
                event.times.y1 = event.times.y1, t1 = t1, t2 = t2, x = x, 
                Lambda0.t12.hat = Lambda0.t12.hat, S.t12.hat = S.t12.hat, 
                infl1.beta = infl1.beta, infl2.beta = infl2.beta,
                infl1.lambda0.t = infl1.lambda0.t,
                infl2.lambda0.t = infl2.lambda0.t,
                infl1.Lambda0.t12 = infl1.Lambda0.t12,
                infl2.Lambda0.t12 = infl2.Lambda0.t12))
  }
}



## -----------------------------------------------------------------------------
## Function: influences.generalized.calib()
## -----------------------------------------------------------------------------
## Description: This function estimates the influences on the cause-specific 
##              log-relative hazard, cause-specific baseline hazard at each 
##              unique event time, as well as cause-specific cumulative baseline
##              hazard until each event time of the event of primary interest in
##              interval (t1,t2]. This function should be used with calibrated
##              weights under the case-cohort design with stratified sampling 
##              based on case status (Breslow.all or Breslow.weight). If the 
##              event type is not that of primary interest, interval (t1,t2] and 
##              the event times of the primary event of interest need to be 
##              provided. Because calibrated weights are used, the phase-one 
##              influences and the phase-two influences are returned. Also 
##              returns parameters estimates
## -----------------------------------------------------------------------------
## Arguments:
##
##  mod             a cox model object, result of function coxph 
##
##  x               vector of length p, specifying the covariate profile 
##                  considered for the pure risk. Default is (0,...,0)
##
##  A               matrix with the values of the auxiliary variables used for 
##                  the calibration of the weights in the whole cohort
##
##  competing       is the event type that of primary interest or a competing 
##                  event. If competing = TRUE, arguments t1, t2, and
##                  event.times.y1 below need to be provided. Default is FALSE
##                  
##  t1              left bound of the time interval considered for the 
##                  cumulative baseline hazard and pure risk. If 
##                  competing = FALSE, the default is the first event time. If 
##                  competing = TRUE, this argument needs to be provided
##                  
##  t2              right bound of the time interval considered for the 
##                  cumulative baseline hazard and pure risk. If 
##                  competing = FALSE, the default is the last event time. If 
##                  competing = TRUE, this argument needs to be provided
##
##  event.times.y1  vector with the event times for the primary event of  
##                  interest (y1). If competing = TRUE, this argument needs to
##                  be provided
##
##  fraction.y      was only a fraction of the cases included in the case-cohort? 
##                  If fraction.y = TRUE, the arguments below need to be 
##                  provided. If fraction.y = FALSE, 
##                  influences.generalized.calib() is equivalent to using 
##                  influences.calib()
##
##  Breslow.all     was only a fraction of the cases included in the case-cohort 
##                  and were the event times of all the cases in the cohort used
##                  in the numerator of the Breslow estimator? If 
##                  Breslow.all = TRUE, the arguments below need to be provided. 
##                  Default is TRUE
##
##  riskmat.phase1  at risk matrix for the cohort (phase-one data) at all the  
##                  cases event times, even that not in phase-two. If 
##                  Breslow.all = TRUE, this argument needs to be provided
##
##  status.phase1   vector indicating the case status in the case-cohort (phase- 
##                  one data). If Breslow.all = TRUE and dNt.phase1 = NULL, this 
##                  argument needs to be provided
##
##  dNt.phase1      counting process matrix for failures in the cohort (phase-
##                  one data). If Breslow.all = TRUE and status.phase1 = NULL, 
##                  this argument needs to be provided
## -----------------------------------------------------------------------------

influences.generalized.calib <- function (mod, A, competing = NULL, 
                                          event.times.y1 = NULL, 
                                          Breslow.all = NULL, 
                                          riskmat.phase1 = NULL, 
                                          dNt.phase1 = NULL, 
                                          status.phase1 = NULL, t1 = NULL, 
                                          t2 = NULL, x = NULL, 
                                          fraction.y = NULL) {
  
  if (is.null(competing)) {
    competing <- FALSE
  }
  
  if (is.null(Breslow.all)) {
    Breslow.all <- TRUE
  }
  
  if (is.null(fraction.y)) {
    fraction.y <- TRUE
  }
  
  if (fraction.y == FALSE & Breslow.all == TRUE) {
    stop("Breslow.all can only be used in the setting where only a fraction of 
         the cases included in the case-cohort")
  }
  
  # If only a fraction of the cases is used ------------------------------------
  
  if (fraction.y == TRUE) {
    if (is.null(riskmat.phase1)) {
      stop("The at risk matrix for the phase-one data at all of the cases event 
           times need to be provided")
    } else {
      if (is.null(dNt.phase1)&is.null(status.phase1)) {
        stop("The status for the phase-one data, or the counting process matrix 
             for failures in the phase-one data, need to be provided")
      } else {
        if (is.null(dNt.phase1)) { 
          
          # If not provided, compute the counting process matrix for failures in 
          # the phase-one data -------------------------------------------------
          
          observed.times.phase1 <- apply(riskmat.phase1, 1,
                                         function(v) {which.max(cumsum(v))})
          dNt.phase1            <- matrix(0, nrow(riskmat.phase1), 
                                          ncol(riskmat.phase1))
          dNt.phase1[cbind(1:nrow(riskmat.phase1), observed.times.phase1)] <- 1
          dNt.phase1            <- sweep(dNt.phase1, 1, status.phase1, "*") 
        }
        
        if(Breslow.all == TRUE){ # we use all of the cases' failure time, 
          # even for cases not in phase-two, in the numerator of the Breslow 
          # estimator
          
          # --------------------------------------------------------------------
          # Quantities needed for the influences -------------------------------
          
          mod.detail            <- coxph.detail(mod, riskmat = TRUE)
          X                     <- model.matrix(mod) 
          riskmat               <- mod.detail$riskmat
          number.times.phase2   <- ncol(riskmat)
          n.phase2              <- nrow(X)
          n                     <- nrow(riskmat.phase1)
          number.times.phase1   <- ncol(riskmat.phase1)
          indiv.phase2          <- row.names(X)
          riskmat.phase2        <- riskmat.phase1[indiv.phase2, ] # use the cases' 
          # actual failure time, known even for cases not in phase-two
          cases.times           <- as.numeric(colnames(riskmat.phase1)) # cases failure times, even 
          # the one who are not in phase-two
          
          weights           <- mod$weights # weights used for the fit, w_i
          if (is.null(weights)) { weights <- 1 } # when using the whole cohort, w_i = 1
          
          beta.hat          <- mod$coefficients
          
          p                 <- length(beta.hat)
          if ((length(x) != p) | (is.null(x))) {x <- rep(0, p)} # if no covariate profile provided, use 0 as reference level
          exp.x             <- c(exp(x %*% beta.hat))
          
          exp.X.weighted    <- weights * exp(X %*% beta.hat)
          Y.exp.X.weighted.casestimes   <- sweep(riskmat.phase2, 1, 
                                                 exp.X.weighted, "*")
          
          observed.times    <- apply(riskmat, 1, function(v) {which.max(cumsum(v))}) 
          dNt               <- matrix(0, n.phase2, number.times.phase2) 
          dNt[cbind(1:nrow(riskmat), observed.times)] <- 1 
          dNt               <- dNt * matrix(mod$y[, ncol(mod$y)], nrow(riskmat), 
                                            number.times.phase2, byrow = FALSE) 
          dNt.weighted      <- dNt * matrix(weights, nrow = n.phase2, 
                                            ncol = number.times.phase2,  byrow = FALSE) 
          infomat.indiv     <- mod.detail$imat  
          infomat           <- apply(X = infomat.indiv, MARGIN = c(1,2), FUN = sum)
          
          Y.exp.X.weighted  <- sweep(riskmat, 1, exp.X.weighted, "*") 
          
          S0t   <- t(riskmat) %*% (exp.X.weighted) 
          S1t   <- t(riskmat) %*% (sweep(X, 1, exp.X.weighted, "*")) 
          S1t   <- t(riskmat) %*% (X * matrix(exp.X.weighted, nrow = n.phase2, ncol = p, byrow = FALSE)) 
          
          X.expect <- S1t / matrix(S0t, nrow = number.times.phase2, ncol = p, byrow = FALSE)
          
          S0t.casestimes  <- t(riskmat.phase2) %*% (exp.X.weighted)  # for any
          # event time t (even the one of cases not in phase-two)
          S1t.casestimes  <- t(riskmat.phase2) %*% (sweep(X, 1, exp.X.weighted, "*"))
          
          # --------------------------------------------------------------------
          # Estimation of the cause-specific baseline hazard at each unique event 
          # time ---------------------------------------------------------------
          
          lambda0.t.hat         <- t(colSums(dNt.phase1) / S0t.casestimes) 
          
          # --------------------------------------------------------------------
          # Estimation of the cause-specific "survival function" at each unique 
          # event time ---------------------------------------------------------
          
          if(competing == FALSE){
            if (is.null(t1)) { t1 <- floor(min(mod.detail$time)) }
            if (is.null(t2)) { t2 <- floor(max(mod.detail$time)) }
            t1t2.times        <- which((t1 < cases.times) & (cases.times <= t2)) 
            
            lambda0.t12.hat   <- lambda0.t.hat[t1t2.times] # baseline hazard for each of the event times in (t1, t2]
            Lambda0.t12.hat   <- cumsum(lambda0.t12.hat) # cumulative baseline hazard in time interval (t1, t] for each of the event times in (t1, t2]
            S.t12.hat         <- exp(- Lambda0.t12.hat * exp.x) # vector with the S(t) / S(t1) for t in (t1, t2]
            
          }else{
            if (is.null(t1)) { stop("t1 should be the same as for the event of primary 
                            interest") }
            if (is.null(t2)) { stop("t2 should be the same as for the event of primary 
                            interest") }
            if (is.null(event.times.y1)) { stop("the event times for the event of 
                                        primary interest should be provided") }
            
            Lambda0.t12.hat   <- sapply(event.times.y1, 
                                        function(t) {sum(lambda0.t.hat[which((t1 <  cases.times) & 
                                                                               (cases.times <= t))])}) # cumulative baseline hazard in time interval (t1, t] for each of the event times of the primary event of interest in (t1, t2]
            S.t12.hat         <- exp(- Lambda0.t12.hat * exp.x)
          }
          
          # ----------------------------------------------------------------------------
          # Computation of the influences for the Lagrangian multipliers, eta ----------
          
          A.phase2          <- A[indiv.phase2, ]
          q                 <- ncol(A.phase2)
          
          A.weighted         <- sweep(A.phase2, 1, weights, "*")
          AA.weighted        <- array(NA, dim = c(q, q, n.phase2))
          for(i in 1:nrow(A.phase2)) {
            AA.weighted[,, i] <- tcrossprod(A.phase2[i,], A.weighted[i,])
          }
          sum.AA.weighted    <- apply(X = AA.weighted, MARGIN = c(1,2), FUN = sum)
          sum.AA.weighted.inv <- solve(sum.AA.weighted)
          
          infl1.eta         <- A %*% sum.AA.weighted.inv
          infl2.eta         <- 0 * A
          infl2.eta[indiv.phase2,] <- - A.weighted %*% sum.AA.weighted.inv
          
          # ----------------------------------------------------------------------------
          # Computation of the influences for log-relative hazard, beta ----------------
          
          XA              <- array(NA, dim = c(p, q, n.phase2))
          for(i in 1:n.phase2) {
            XA[,, i]      <- tcrossprod(X[i,], A.phase2[i,])
          }
          drond.G1t.eta   <- array(NA, dim = c(p, q, number.times.phase2))
          drond.G0t.eta   <- array(NA, dim = c(q, number.times.phase2))
          drond.S1t.eta   <- array(NA, dim = c(p, q, number.times.phase2))
          drond.S0t.eta   <- array(NA, dim = c(q, number.times.phase2))
          for(t in 1:number.times.phase2) { 
            drond.G1t.eta[,, t] <- apply(X = sweep(XA, 3, dNt.weighted[, t], "*"), 
                                         MARGIN = c(1,2), FUN = sum)
            drond.G0t.eta[, t]  <- colSums(A.phase2 * matrix(dNt.weighted[, t], 
                                                             nrow = n.phase2, ncol = q, 
                                                             byrow = FALSE))
            drond.S1t.eta[,, t] <- apply(X = sweep(XA, 3, Y.exp.X.weighted[, t],
                                                   "*"), MARGIN = c(1,2), FUN = sum)
            drond.S0t.eta[, t]  <- colSums(A.phase2 * matrix(Y.exp.X.weighted[, t], 
                                                             nrow = n.phase2, 
                                                             ncol = q, 
                                                             byrow = FALSE))
          }
          
          drond.S1t.casestimes.eta   <- array(NA, dim = c(p, q, number.times.phase1))
          drond.S0t.casestimes.eta   <- array(NA, dim = c(q, number.times.phase1))
          for(t in 1:number.times.phase1) { 
            drond.S1t.casestimes.eta[,, t] <- apply(X = sweep(XA, 3, Y.exp.X.weighted.casestimes [, t],
                                                              "*"), MARGIN = c(1,2), FUN = sum)
            drond.S0t.casestimes.eta[, t]  <- colSums(A.phase2 * matrix(Y.exp.X.weighted.casestimes [, t], 
                                                                        nrow = n.phase2, 
                                                                        ncol = q, 
                                                                        byrow = FALSE))
          }
          
          drond.Ut.eta    <- array(NA, dim = c(p, q, number.times.phase2))
          X.expect        <- S1t / matrix(S0t, nrow = number.times.phase2, 
                                          ncol = p, byrow = FALSE)
          for(t in 1:number.times.phase2) { 
            drond.Ut.eta[,, t] <- drond.G1t.eta[,, t] - tcrossprod(X.expect[t, ],
                                                                   drond.G0t.eta[, t]) -
              (colSums(dNt.weighted) / S0t)[t] * drond.S1t.eta[,, t] +
              (colSums(dNt.weighted) / S0t ^ 2)[t] * 
              tcrossprod(S1t[t, ], drond.S0t.eta[, t])
          }
          drond.U.eta     <- apply(X = drond.Ut.eta, MARGIN = c(1,2), FUN = sum)
          infl1.beta      <- infl1.eta %*% t(drond.U.eta) %*% solve(infomat)
          score           <- array(NA, dim = c(n.phase2, p, number.times.phase2))
          for(i in 1:n.phase2) {
            score[i,,]    <- - sweep(t(X.expect), 1, X[i,], "-")
          }
          for(t in 1:number.times.phase2) {
            dMt           <- dNt.weighted[, t] - riskmat[, t] *
              exp.X.weighted * (colSums(dNt.weighted)/S0t)[t]
            score[,, t]   <- score[,, t] * matrix(dMt, nrow = n.phase2, ncol = p, 
                                                  byrow = FALSE)
          }
          score.beta      <- apply(X = score, MARGIN = c(1,2), FUN = sum) 
          infl2.beta      <- infl2.eta %*% t(drond.U.eta) %*% solve(infomat)
          infl2.beta[indiv.phase2,] <- infl2.beta[indiv.phase2,] +
            score.beta %*% solve(infomat) 
          
          # --------------------------------------------------------------------------
          # Computation of the influences on the cause-specific the baseline hazard at
          # each unique event time ---------------------------------------------------
          
          infl1.lambda0.t <- dNt.phase1 / matrix(S0t.casestimes, nrow = n, 
                                                 ncol = number.times.phase1, 
                                                 byrow = TRUE) + 
            (- (infl1.beta %*% t(S1t.casestimes) + infl1.eta %*% drond.S0t.casestimes.eta) *
               matrix(lambda0.t.hat, nrow = nrow(A), 
                      ncol = number.times.phase1, byrow = TRUE)) /
            matrix(S0t.casestimes, nrow = nrow(A), ncol = number.times.phase1, 
                   byrow = TRUE) 
          
          infl2.lambda0.t <- 0 * infl1.lambda0.t
          infl2.lambda0.t[indiv.phase2,] <- (- (Y.exp.X.weighted.casestimes + 
                                                  infl2.beta[indiv.phase2,] %*% 
                                                  t(S1t.casestimes) + infl2.eta[indiv.phase2,] %*%
                                                  drond.S0t.casestimes.eta) * 
                                               matrix(lambda0.t.hat, nrow = n.phase2, 
                                                      ncol = number.times.phase1, byrow = TRUE)) /
            matrix(S0t.casestimes, nrow = n.phase2, ncol = number.times.phase1,
                   byrow = TRUE) 
          
          if(competing == FALSE){
            infl1.lambda0.t12     <- infl1.lambda0.t[,t1t2.times]
            infl1.Lambda0.t12     <- t(apply(infl1.lambda0.t12, 1, cumsum))
            
            infl2.lambda0.t12     <- infl2.lambda0.t[,t1t2.times]
            infl2.Lambda0.t12     <- t(apply(infl2.lambda0.t12, 1, cumsum))
            
            return(list(beta.hat = beta.hat, exp.x = exp.x, 
                        lambda0.t.hat = lambda0.t.hat, competing = competing,
                        event.times.y1 = cases.times[t1t2.times], t1 = t1, 
                        t2 = t2, x = x, lambda0.t12.hat = lambda0.t12.hat, 
                        Lambda0.t12.hat = Lambda0.t12.hat, S.t12.hat = S.t12.hat,
                        infl1.beta = infl1.beta, infl2.beta = infl2.beta, 
                        infl1.lambda0.t = infl1.lambda0.t,
                        infl2.lambda0.t = infl2.lambda0.t,
                        infl1.lambda0.t12 = infl1.lambda0.t12,
                        infl2.lambda0.t12 = infl2.lambda0.t12,
                        infl1.Lambda0.t12 = infl1.Lambda0.t12,
                        infl2.Lambda0.t12 = infl2.Lambda0.t12))
          }else{
            infl1.Lambda0.t12   <- sapply(event.times.y1, function(t) {rowSums(as.matrix(infl1.lambda0.t[,which((t1 < cases.times) & 
                                                                                                                  (cases.times <= t))]))})
            infl2.Lambda0.t12   <- sapply(event.times.y1, function(t) {rowSums(as.matrix(infl2.lambda0.t[,which((t1 < cases.times) & 
                                                                                                                  (cases.times <= t))]))})
            return(list(beta.hat = beta.hat, exp.x = exp.x, 
                        lambda0.t.hat = lambda0.t.hat, competing = competing, 
                        event.times.y1 = event.times.y1, t1 = t1, t2 = t2, x = x, 
                        Lambda0.t12.hat = Lambda0.t12.hat, S.t12.hat = S.t12.hat, 
                        infl1.beta = infl1.beta, infl2.beta = infl2.beta,
                        infl1.lambda0.t = infl1.lambda0.t,
                        infl2.lambda0.t = infl2.lambda0.t,
                        infl1.Lambda0.t12 = infl1.Lambda0.t12,
                        infl2.Lambda0.t12 = infl2.Lambda0.t12))
          }
          
        }else{ # we only use the event times of cases in phase-two and weight 
          # the numerator of the Breslow estimator
          
          # ----------------------------------------------------------------------------
          # Quantities needed for the influences ---------------------------------------
          
          mod.detail        <- coxph.detail(mod, riskmat = T)
          riskmat           <- mod.detail$riskmat 
          number.times      <- ncol(riskmat) 
          n                 <- nrow(riskmat) 
          X                 <- model.matrix(mod) 
          p                 <- ncol(X) 
          weights           <- mod$weights # weights used for the fit, w_i
          if (is.null(weights)) { weights <- 1 } # when using the whole cohort, w_i = 1
          
          beta.hat          <- mod$coefficients
          
          p                 <- length(beta.hat)
          if ((length(x) != p) | (is.null(x))) {x <- rep(0, p)} # if no covariate profile provided, use 0 as reference level
          exp.x             <- c(exp(x %*% beta.hat))
          
          exp.X.weighted    <- weights * exp(X %*% beta.hat)
          Y.exp.X.weighted  <- riskmat * matrix(exp.X.weighted, nrow = n, 
                                                ncol = number.times, byrow = FALSE) 
          S0t               <- t(riskmat) %*% (exp.X.weighted) 
          S1t               <- t(riskmat) %*% (X * matrix(exp.X.weighted, nrow = n,
                                                          ncol = p, byrow = FALSE)) 
          observed.times    <- apply(riskmat, 1, function(v) {which.max(cumsum(v))}) 
          dNt               <- matrix(0, n, number.times) 
          dNt[cbind(1:n, observed.times)] <- 1 
          dNt               <- dNt * matrix(mod$y[, ncol(mod$y)], nrow(riskmat), 
                                            number.times, byrow = FALSE) 
          dNt.weighted      <- dNt * matrix(weights, nrow = n, ncol = number.times, 
                                            byrow = FALSE) 
          
          infomat.indiv     <- mod.detail$imat  
          infomat           <- apply(X = infomat.indiv, MARGIN = c(1,2), FUN = sum)
          
          # ----------------------------------------------------------------------------
          # Estimation of the cause-specific baseline hazard at each unique event time -
          
          lambda0.t.hat         <- t(colSums(dNt.weighted) / S0t) 
          
          # ----------------------------------------------------------------------------
          # Estimation of the cause-specific "survival function" at each unique event 
          # time -----------------------------------------------------------------------
          
          if(competing == FALSE){
            if (is.null(t1)) { t1 <- floor(min(mod.detail$time)) }
            if (is.null(t2)) { t2 <- floor(max(mod.detail$time)) }
            t1t2.times        <- which((t1 < mod.detail$time) & (mod.detail$time <= t2)) 
            
            lambda0.t12.hat   <- lambda0.t.hat[t1t2.times] # baseline hazard for each of the event times in (t1, t2]
            Lambda0.t12.hat   <- cumsum(lambda0.t12.hat) # cumulative baseline hazard in time interval (t1, t] for each of the event times t in (t1, t2]
            S.t12.hat         <- exp(- Lambda0.t12.hat * exp.x) # vector with the S(t) / S(t1) for each of the event times t in (t1, t2]
          }else{
            if (is.null(t1)) { stop("t1 should be the same as for the event of primary 
                            interest") }
            if (is.null(t2)) { stop("t2 should be the same as for the event of primary 
                            interest") }
            if (is.null(event.times.y1)) { stop("the event times for the event of 
                                        primary interest should be provided") }
            
            Lambda0.t12.hat   <- sapply(event.times.y1, 
                                        function(t) {sum(lambda0.t.hat[which((t1 < mod.detail$time) & 
                                                                               (mod.detail$time <= t))])}) # cumulative baseline hazard in time interval (t1, t] for each of the event times of the primary event of interest in (t1, t2]
            S.t12.hat         <- exp(- Lambda0.t12.hat * exp.x)
          }
          
          # ----------------------------------------------------------------------------
          # Computation of the influences for the Lagrangian multipliers, eta ----------
          
          indiv.phase2      <- row.names(X)
          A.phase2          <- A[indiv.phase2, ]
          q                 <- ncol(A.phase2)
          
          A.weighted         <- sweep(A.phase2, 1, weights, "*")
          AA.weighted        <- array(NA, dim = c(q, q, n))
          for(i in 1:nrow(A.phase2)) {
            AA.weighted[,, i] <- tcrossprod(A.phase2[i,], A.weighted[i,])
          }
          sum.AA.weighted    <- apply(X = AA.weighted, MARGIN = c(1,2), FUN = sum)
          sum.AA.weighted.inv <- solve(sum.AA.weighted)
          
          infl1.eta         <- A %*% sum.AA.weighted.inv
          infl2.eta         <- 0 * A
          infl2.eta[indiv.phase2,] <- - A.weighted %*% sum.AA.weighted.inv
          infl.eta          <- infl1.eta + infl2.eta 
          
          # ----------------------------------------------------------------------------
          # Computation of the influences for log-relative hazard, beta ----------------
          
          XA              <- array(NA, dim = c(p, q, n))
          for(i in 1:n) {
            XA[,, i]      <- tcrossprod(X[i,], A.phase2[i,])
          }
          drond.G1t.eta   <- array(NA, dim = c(p, q, number.times))
          drond.G0t.eta   <- array(NA, dim = c(q, number.times))
          drond.S1t.eta   <- array(NA, dim = c(p, q, number.times))
          drond.S0t.eta   <- array(NA, dim = c(q, number.times))
          for(t in 1:number.times) { 
            drond.G1t.eta[,, t] <- apply(X = sweep(XA, 3, dNt.weighted[, t], "*"), 
                                         MARGIN = c(1,2), FUN = sum)
            drond.G0t.eta[, t]  <- colSums(A.phase2 * matrix(dNt.weighted[, t], 
                                                             nrow = n, ncol = q, 
                                                             byrow = FALSE))
            drond.S1t.eta[,, t] <- apply(X = sweep(XA, 3, Y.exp.X.weighted[, t],
                                                   "*"), MARGIN = c(1,2), FUN = sum)
            drond.S0t.eta[, t]  <- colSums(A.phase2 * matrix(Y.exp.X.weighted[, t], 
                                                             nrow = n, ncol = q, 
                                                             byrow = FALSE))
          }
          drond.Ut.eta    <- array(NA, dim = c(p, q, number.times))
          X.expect        <- S1t / matrix(S0t, nrow = number.times, ncol = p, 
                                          byrow = FALSE)
          for(t in 1:number.times) { 
            drond.Ut.eta[,, t] <- drond.G1t.eta[,, t] - tcrossprod(X.expect[t, ],
                                                                   drond.G0t.eta[, t]) -
              (colSums(dNt.weighted) / S0t)[t] * drond.S1t.eta[,, t] +
              (colSums(dNt.weighted) / S0t ^ 2)[t] * 
              tcrossprod(S1t[t, ], drond.S0t.eta[, t])
          }
          drond.U.eta     <- apply(X = drond.Ut.eta, MARGIN = c(1,2), FUN = sum)
          infl1.beta      <- infl1.eta %*% t(drond.U.eta) %*% solve(infomat)
          score           <- array(NA, dim = c(n, p, number.times))
          for(i in 1:n) {
            score[i,,]    <- - sweep(t(X.expect), 1, X[i,], "-")
          }
          for(t in 1:number.times) {
            dMt           <- dNt.weighted[, t] - riskmat[, t] *
              exp.X.weighted * (colSums(dNt.weighted)/S0t)[t]
            score[,, t]   <- score[,, t] * matrix(dMt, nrow = n, ncol = p, 
                                                  byrow = FALSE)
          }
          score.beta      <- apply(X = score, MARGIN = c(1,2), FUN = sum) 
          infl2.beta      <- infl2.eta %*% t(drond.U.eta) %*% solve(infomat)
          infl2.beta[indiv.phase2,] <- infl2.beta[indiv.phase2,] +
            score.beta %*% solve(infomat) 
          
          # --------------------------------------------------------------------------
          # Computation of the influences on the cause-specific the baseline hazard at
          # each unique event time ---------------------------------------------------
          
          infl1.lambda0.t <- (- (infl1.beta %*% t(S1t) + infl1.eta %*% drond.S0t.eta) *
                                matrix(lambda0.t.hat, nrow = nrow(A), 
                                       ncol = number.times, byrow = TRUE)) /
            matrix(S0t, nrow = nrow(A), ncol = number.times, 
                   byrow = TRUE) 
          infl2.lambda0.t <- 0 * infl1.lambda0.t
          infl2.lambda0.t[indiv.phase2,] <- (dNt.weighted - (Y.exp.X.weighted + 
                                                               infl2.beta[indiv.phase2,] %*% 
                                                               t(S1t) + infl2.eta[indiv.phase2,] %*%
                                                               drond.S0t.eta) * 
                                               matrix(lambda0.t.hat, nrow = n, 
                                                      ncol = number.times, byrow = TRUE)) /
            matrix(S0t, nrow = n, ncol = number.times,
                   byrow = TRUE) 
          
          if(competing == FALSE){
            infl1.lambda0.t12     <- infl1.lambda0.t[,t1t2.times]
            infl1.Lambda0.t12     <- t(apply(infl1.lambda0.t12, 1, cumsum))
            
            infl2.lambda0.t12     <- infl2.lambda0.t[,t1t2.times]
            infl2.Lambda0.t12     <- t(apply(infl2.lambda0.t12, 1, cumsum))
            
            return(list(beta.hat = beta.hat, exp.x = exp.x, 
                        lambda0.t.hat = lambda0.t.hat, competing = competing,
                        event.times.y1 = mod.detail$time[t1t2.times], t1 = t1, 
                        t2 = t2, x = x, lambda0.t12.hat = lambda0.t12.hat, 
                        Lambda0.t12.hat = Lambda0.t12.hat, S.t12.hat = S.t12.hat,
                        infl1.beta = infl1.beta, infl2.beta = infl2.beta, 
                        infl1.lambda0.t = infl1.lambda0.t,
                        infl2.lambda0.t = infl2.lambda0.t,
                        infl1.lambda0.t12 = infl1.lambda0.t12,
                        infl2.lambda0.t12 = infl2.lambda0.t12,
                        infl1.Lambda0.t12 = infl1.Lambda0.t12,
                        infl2.Lambda0.t12 = infl2.Lambda0.t12))
          }else{
            infl1.Lambda0.t12   <- sapply(event.times.y1, function(t) {rowSums(as.matrix(infl1.lambda0.t[,which((t1 < mod.detail$time) & 
                                                                                                                  (mod.detail$time <= t))]))})
            infl2.Lambda0.t12   <- sapply(event.times.y1, function(t) {rowSums(as.matrix(infl2.lambda0.t[,which((t1 < mod.detail$time) & 
                                                                                                                  (mod.detail$time <= t))]))})
            
            return(list(beta.hat = beta.hat, exp.x = exp.x, 
                        lambda0.t.hat = lambda0.t.hat, competing = competing, 
                        event.times.y1 = event.times.y1, t1 = t1, t2 = t2, x = x, 
                        Lambda0.t12.hat = Lambda0.t12.hat, S.t12.hat = S.t12.hat, 
                        infl1.beta = infl1.beta, infl2.beta = infl2.beta,
                        infl1.lambda0.t = infl1.lambda0.t,
                        infl2.lambda0.t = infl2.lambda0.t,
                        infl1.Lambda0.t12 = infl1.Lambda0.t12,
                        infl2.Lambda0.t12 = infl2.Lambda0.t12))
          }
        }
      }
    }
    
  }else{
    
    # ----------------------------------------------------------------------------
    # Quantities needed for the influences ---------------------------------------
    
    mod.detail        <- coxph.detail(mod, riskmat = T)
    riskmat           <- mod.detail$riskmat 
    number.times      <- ncol(riskmat) 
    n                 <- nrow(riskmat) 
    X                 <- model.matrix(mod) 
    p                 <- ncol(X) 
    weights           <- mod$weights # weights used for the fit, w_i
    if (is.null(weights)) { weights <- 1 } # when using the whole cohort, w_i = 1
    
    beta.hat          <- mod$coefficients
    
    p                 <- length(beta.hat)
    if ((length(x) != p) | (is.null(x))) {x <- rep(0, p)} # if no covariate profile provided, use 0 as reference level
    exp.x             <- c(exp(x %*% beta.hat))
    
    exp.X.weighted    <- weights * exp(X %*% beta.hat)
    Y.exp.X.weighted  <- riskmat * matrix(exp.X.weighted, nrow = n, 
                                          ncol = number.times, byrow = FALSE) 
    S0t               <- t(riskmat) %*% (exp.X.weighted) 
    S1t               <- t(riskmat) %*% (X * matrix(exp.X.weighted, nrow = n,
                                                    ncol = p, byrow = FALSE)) 
    observed.times    <- apply(riskmat, 1, function(v) {which.max(cumsum(v))}) 
    dNt               <- matrix(0, n, number.times) 
    dNt[cbind(1:n, observed.times)] <- 1 
    dNt               <- dNt * matrix(mod$y[, ncol(mod$y)], nrow(riskmat), 
                                      number.times, byrow = FALSE) 
    dNt.weighted      <- dNt * matrix(weights, nrow = n, ncol = number.times, 
                                      byrow = FALSE) 
    
    infomat.indiv     <- mod.detail$imat  
    infomat           <- apply(X = infomat.indiv, MARGIN = c(1,2), FUN = sum)
    
    # ----------------------------------------------------------------------------
    # Estimation of the cause-specific baseline hazard at each unique event time -
    
    lambda0.t.hat         <- t(colSums(dNt) / S0t) 
    
    # ----------------------------------------------------------------------------
    # Estimation of the cause-specific "survival function" at each unique event 
    # time -----------------------------------------------------------------------
    
    if(competing == FALSE){
      if (is.null(t1)) { t1 <- floor(min(mod.detail$time)) }
      if (is.null(t2)) { t2 <- floor(max(mod.detail$time)) }
      t1t2.times        <- which((t1 < mod.detail$time) & (mod.detail$time <= t2)) 
      
      lambda0.t12.hat   <- lambda0.t.hat[t1t2.times] # baseline hazard for each of the event times in (t1, t2]
      Lambda0.t12.hat   <- cumsum(lambda0.t12.hat) # cumulative baseline hazard in time interval (t1, t] for each of the event times t in (t1, t2]
      S.t12.hat         <- exp(- Lambda0.t12.hat * exp.x) # vector with the S(t) / S(t1) for each of the event times t in (t1, t2]
    }else{
      if (is.null(t1)) { stop("t1 should be the same as for the event of primary 
                            interest") }
      if (is.null(t2)) { stop("t2 should be the same as for the event of primary 
                            interest") }
      if (is.null(event.times.y1)) { stop("the event times for the event of 
                                        primary interest should be provided") }
      
      Lambda0.t12.hat   <- sapply(event.times.y1, 
                                  function(t) {sum(lambda0.t.hat[which((t1 < mod.detail$time) & 
                                                                         (mod.detail$time <= t))])}) # cumulative baseline hazard in time interval (t1, t] for each of the event times of the primary event of interest in (t1, t2]
      S.t12.hat         <- exp(- Lambda0.t12.hat * exp.x)
    }
    
    # ----------------------------------------------------------------------------
    # Computation of the influences for the Lagrangian multipliers, eta ----------
    
    indiv.phase2      <- row.names(X)
    A.phase2          <- A[indiv.phase2, ]
    q                 <- ncol(A.phase2)
    
    A.weighted         <- sweep(A.phase2, 1, weights, "*")
    AA.weighted        <- array(NA, dim = c(q, q, n))
    for(i in 1:nrow(A.phase2)) {
      AA.weighted[,, i] <- tcrossprod(A.phase2[i,], A.weighted[i,])
    }
    sum.AA.weighted    <- apply(X = AA.weighted, MARGIN = c(1,2), FUN = sum)
    sum.AA.weighted.inv <- solve(sum.AA.weighted)
    
    infl1.eta         <- A %*% sum.AA.weighted.inv
    infl2.eta         <- 0 * A
    infl2.eta[indiv.phase2,] <- - A.weighted %*% sum.AA.weighted.inv
    infl.eta          <- infl1.eta + infl2.eta 
    
    # ----------------------------------------------------------------------------
    # Computation of the influences for log-relative hazard, beta ----------------
    
    XA              <- array(NA, dim = c(p, q, n))
    for(i in 1:n) {
      XA[,, i]      <- tcrossprod(X[i,], A.phase2[i,])
    }
    drond.G1t.eta   <- array(NA, dim = c(p, q, number.times))
    drond.G0t.eta   <- array(NA, dim = c(q, number.times))
    drond.S1t.eta   <- array(NA, dim = c(p, q, number.times))
    drond.S0t.eta   <- array(NA, dim = c(q, number.times))
    for(t in 1:number.times) { 
      drond.G1t.eta[,, t] <- apply(X = sweep(XA, 3, dNt.weighted[, t], "*"), 
                                   MARGIN = c(1,2), FUN = sum)
      drond.G0t.eta[, t]  <- colSums(A.phase2 * matrix(dNt.weighted[, t], 
                                                       nrow = n, ncol = q, 
                                                       byrow = FALSE))
      drond.S1t.eta[,, t] <- apply(X = sweep(XA, 3, Y.exp.X.weighted[, t],
                                             "*"), MARGIN = c(1,2), FUN = sum)
      drond.S0t.eta[, t]  <- colSums(A.phase2 * matrix(Y.exp.X.weighted[, t], 
                                                       nrow = n, ncol = q, 
                                                       byrow = FALSE))
    }
    drond.Ut.eta    <- array(NA, dim = c(p, q, number.times))
    X.expect        <- S1t / matrix(S0t, nrow = number.times, ncol = p, 
                                    byrow = FALSE)
    for(t in 1:number.times) { 
      drond.Ut.eta[,, t] <- drond.G1t.eta[,, t] - tcrossprod(X.expect[t, ],
                                                             drond.G0t.eta[, t]) -
        (colSums(dNt.weighted) / S0t)[t] * drond.S1t.eta[,, t] +
        (colSums(dNt.weighted) / S0t ^ 2)[t] * 
        tcrossprod(S1t[t, ], drond.S0t.eta[, t])
    }
    drond.U.eta     <- apply(X = drond.Ut.eta, MARGIN = c(1,2), FUN = sum)
    infl1.beta      <- infl1.eta %*% t(drond.U.eta) %*% solve(infomat)
    score           <- array(NA, dim = c(n, p, number.times))
    for(i in 1:n) {
      score[i,,]    <- - sweep(t(X.expect), 1, X[i,], "-")
    }
    for(t in 1:number.times) {
      dMt           <- dNt.weighted[, t] - riskmat[, t] *
        exp.X.weighted * (colSums(dNt.weighted)/S0t)[t]
      score[,, t]   <- score[,, t] * matrix(dMt, nrow = n, ncol = p, 
                                            byrow = FALSE)
    }
    score.beta      <- apply(X = score, MARGIN = c(1,2), FUN = sum) 
    infl2.beta      <- infl2.eta %*% t(drond.U.eta) %*% solve(infomat)
    infl2.beta[indiv.phase2,] <- infl2.beta[indiv.phase2,] +
      score.beta %*% solve(infomat) 
    
    # --------------------------------------------------------------------------
    # Computation of the influences on the cause-specific the baseline hazard at
    # each unique event time ---------------------------------------------------
    
    infl1.lambda0.t <- (- (infl1.beta %*% t(S1t) + infl1.eta %*% drond.S0t.eta) *
                          matrix(lambda0.t.hat, nrow = nrow(A), 
                                 ncol = number.times, byrow = TRUE)) /
      matrix(S0t, nrow = nrow(A), ncol = number.times, 
             byrow = TRUE) 
    infl2.lambda0.t <- 0 * infl1.lambda0.t
    infl2.lambda0.t[indiv.phase2,] <- (dNt - (Y.exp.X.weighted + 
                                                infl2.beta[indiv.phase2,] %*% 
                                                t(S1t) + infl2.eta[indiv.phase2,] %*%
                                                drond.S0t.eta) * 
                                         matrix(lambda0.t.hat, nrow = n, 
                                                ncol = number.times, byrow = TRUE)) /
      matrix(S0t, nrow = n, ncol = number.times,
             byrow = TRUE) 
    
    if(competing == FALSE){
      infl1.lambda0.t12     <- infl1.lambda0.t[,t1t2.times]
      infl1.Lambda0.t12     <- t(apply(infl1.lambda0.t12, 1, cumsum))
      
      infl2.lambda0.t12     <- infl2.lambda0.t[,t1t2.times]
      infl2.Lambda0.t12     <- t(apply(infl2.lambda0.t12, 1, cumsum))
      
      return(list(beta.hat = beta.hat, exp.x = exp.x, 
                  lambda0.t.hat = lambda0.t.hat, competing = competing,
                  event.times.y1 = mod.detail$time[t1t2.times], t1 = t1, 
                  t2 = t2, x = x, lambda0.t12.hat = lambda0.t12.hat, 
                  Lambda0.t12.hat = Lambda0.t12.hat, S.t12.hat = S.t12.hat,
                  infl1.beta = infl1.beta, infl2.beta = infl2.beta, 
                  infl1.lambda0.t = infl1.lambda0.t,
                  infl2.lambda0.t = infl2.lambda0.t,
                  infl1.lambda0.t12 = infl1.lambda0.t12,
                  infl2.lambda0.t12 = infl2.lambda0.t12,
                  infl1.Lambda0.t12 = infl1.Lambda0.t12,
                  infl2.Lambda0.t12 = infl2.Lambda0.t12))
    }else{
      infl1.Lambda0.t12   <- sapply(event.times.y1, function(t) {rowSums(as.matrix(infl1.lambda0.t[,which((t1 < mod.detail$time) & 
                                                                                                            (mod.detail$time <= t))]))})
      infl2.Lambda0.t12   <- sapply(event.times.y1, function(t) {rowSums(as.matrix(infl2.lambda0.t[,which((t1 < mod.detail$time) & 
                                                                                                            (mod.detail$time <= t))]))})
      return(list(beta.hat = beta.hat, exp.x = exp.x, 
                  lambda0.t.hat = lambda0.t.hat, competing = competing, 
                  event.times.y1 = event.times.y1, t1 = t1, t2 = t2, x = x, 
                  Lambda0.t12.hat = Lambda0.t12.hat, S.t12.hat = S.t12.hat, 
                  infl1.beta = infl1.beta, infl2.beta = infl2.beta,
                  infl1.lambda0.t = infl1.lambda0.t,
                  infl2.lambda0.t = infl2.lambda0.t,
                  infl1.Lambda0.t12 = infl1.Lambda0.t12,
                  infl2.Lambda0.t12 = infl2.Lambda0.t12))
    }
  }
}



## -----------------------------------------------------------------------------
## Function: risk.estimation()
## -----------------------------------------------------------------------------
## Description: This function estimates the absolute risk and pure risk for the 
##              primary event of interest. The covariate profile and time or age
##              interval for which the absolute risk and pure risk are to be
##              estimated, are the ones that were provided when estimating the
##              log-relative hazard, baseline hazard and cumulative baseline
##              hazard 
## -----------------------------------------------------------------------------
## Arguments:
##
##  est.y1          return object from function influences() or
##                  influences.generalized(), for the primary event of interest
##
##  est.y2          return object from function influences() or
##                  influences.generalized(), for the competing event
## ----------------------------------------------------------------------------

risk.estimation <- function(est.y1, est.y2) {
  
  AR.y1   <- est.y1$exp.x * sum(est.y1$lambda0.t12.hat * est.y1$S.t12.hat * 
                                  est.y2$S.t12.hat)
  
  PR.y1   <- 1- exp(- est.y1$exp.x * sum(est.y1$lambda0.t12.hat))
  
  return(list(AR.y1 = AR.y1, PR.y1 = PR.y1))
}



## -----------------------------------------------------------------------------
## Function: risk.influences()
## -----------------------------------------------------------------------------
## Description: This function estimates the influences on the absolute risk and 
##              pure risk for the primary event of interest. It should be used
##              with design weights under the case-cohort design with exhaustive 
##              sampling of cases, or under the case-cohort design with 
##              stratified sampling based on case status when using only the 
##              event times of the cases in the case-cohort (Breslow.weight). 
##              Also returns parameters estimates. The covariate profile and 
##              time or age interval for which the absolute risk and pure risk 
##              are to be estimated, are the ones that were provided when 
##              estimating the log-relative hazard, baseline hazard and 
##              cumulative baseline hazard 
## -----------------------------------------------------------------------------
## Arguments:
##
##  est.y1          return object from function influences, for the primary 
##                  event of interest
##
##  est.y2          return object from function influences, for the competing
##                  event
## ----------------------------------------------------------------------------

risk.influences <- function(est.y1, est.y2) {
  
  AR.y1       <- est.y1$exp.x * sum(est.y1$lambda0.t12.hat * est.y1$S.t12.hat * 
                                      est.y2$S.t12.hat)
  
  PR.y1       <- 1- exp(- est.y1$exp.x * sum(est.y1$lambda0.t12.hat))
  
  infl.AR.y1  <- as.matrix(rowSums(est.y1$exp.x * 
                                     sweep(((est.y1$infl.beta %*% 
                                               est.y1$x %*% 
                                               est.y1$lambda0.t12.hat) + 
                                              est.y1$infl.lambda0.t12 - 
                                              sweep(((est.y1$exp.x * 
                                                        (est.y1$infl.beta %*% 
                                                           est.y1$x %*% 
                                                           t(est.y1$Lambda0.t12.hat) + 
                                                           est.y1$infl.Lambda0.t12)) + 
                                                       (est.y2$exp.x * 
                                                          (est.y2$infl.beta %*% 
                                                             est.y2$x %*% 
                                                             t(est.y2$Lambda0.t12.hat) + 
                                                             est.y2$infl.Lambda0.t12))), 2,
                                                    est.y1$lambda0.t12.hat, "*")), 2, 
                                           (est.y1$S.t12.hat * est.y2$S.t12.hat), "*" )))
  
  infl.PR.y1    <- ((1 - PR.y1) * est.y1$exp.x * sum(est.y1$lambda0.t12.hat)) * 
    (est.y1$infl.beta %*% est.y1$x) + 
    ((1 - PR.y1) * est.y1$exp.x) * rowSums(est.y1$infl.lambda0.t12)
  
  return(list(AR.y1 = AR.y1, PR.y1 = PR.y1, infl.AR.y1 = infl.AR.y1,
              infl.PR.y1 = infl.PR.y1))
}



## -----------------------------------------------------------------------------
## Function: risk.influences.generalized()
## -----------------------------------------------------------------------------
## Description: This function estimates the influences on the absolute risk and 
##              pure risk for the primary event of interest. It should be used 
##              with design weights under the case-cohort design with stratified 
##              sampling based on case status when using the event times of all 
##              the cases in the cohort (Breslow.all). Also returns parameters 
##              estimates. The covariate profile and time or age interval for 
##              which the absolute risk and pure risk are to be estimated, are 
##              the ones that were provided when estimating the log-relative 
##              hazard, baseline hazard and cumulative baseline hazard 
## -----------------------------------------------------------------------------
## Arguments:
##
##  est.y1          return object from function influences, for the primary 
##                  event of interest
##
##  est.y2          return object from function influences, for the competing
##                  event
## ----------------------------------------------------------------------------

risk.influences.generalized  <- function(est.y1, est.y2) {
  
  AR.y1       <- est.y1$exp.x * sum(est.y1$lambda0.t12.hat * est.y1$S.t12.hat * 
                                      est.y2$S.t12.hat)
  
  PR.y1       <- 1- exp(- est.y1$exp.x * sum(est.y1$lambda0.t12.hat))
  
  
  infl1.AR.y1  <- as.matrix(rowSums(est.y1$exp.x * 
                                      sweep(((est.y1$infl1.beta %*% 
                                                est.y1$x %*% 
                                                est.y1$lambda0.t12.hat) + 
                                               est.y1$infl1.lambda0.t12 - 
                                               sweep(((est.y1$exp.x * 
                                                         (est.y1$infl1.beta %*% 
                                                            est.y1$x %*% 
                                                            t(est.y1$Lambda0.t12.hat) + 
                                                            est.y1$infl1.Lambda0.t12)) + 
                                                        (est.y2$exp.x * 
                                                           (est.y2$infl1.beta %*% 
                                                              est.y2$x %*% 
                                                              t(est.y2$Lambda0.t12.hat) + 
                                                              est.y2$infl1.Lambda0.t12))), 2,
                                                     est.y1$lambda0.t12.hat, "*")), 2, 
                                            (est.y1$S.t12.hat * est.y2$S.t12.hat), "*" )))
  
  infl1.PR.y1    <- ((1 - PR.y1) * est.y1$exp.x * sum(est.y1$lambda0.t12.hat)) * 
    (est.y1$infl1.beta %*% est.y1$x) + 
    ((1 - PR.y1) * est.y1$exp.x) * rowSums(est.y1$infl1.lambda0.t12)
  
  infl2.AR.y1  <- as.matrix(rowSums(est.y1$exp.x * 
                                      sweep(((est.y1$infl2.beta %*% 
                                                est.y1$x %*% 
                                                est.y1$lambda0.t12.hat) + 
                                               est.y1$infl2.lambda0.t12 - 
                                               sweep(((est.y1$exp.x * 
                                                         (est.y1$infl2.beta %*% 
                                                            est.y1$x %*% 
                                                            t(est.y1$Lambda0.t12.hat) + 
                                                            est.y1$infl2.Lambda0.t12)) + 
                                                        (est.y2$exp.x * 
                                                           (est.y2$infl2.beta %*% 
                                                              est.y2$x %*% 
                                                              t(est.y2$Lambda0.t12.hat) + 
                                                              est.y2$infl2.Lambda0.t12))), 2,
                                                     est.y1$lambda0.t12.hat, "*")), 2, 
                                            (est.y1$S.t12.hat * est.y2$S.t12.hat), "*" )))
  
  infl2.PR.y1    <- ((1 - PR.y1) * est.y1$exp.x * sum(est.y1$lambda0.t12.hat)) * 
    (est.y1$infl2.beta %*% est.y1$x) + 
    ((1 - PR.y1) * est.y1$exp.x) * rowSums(est.y1$infl2.lambda0.t12)
  
  return(list(AR.y1 = AR.y1, PR.y1 = PR.y1, infl1.AR.y1 = infl1.AR.y1,
              infl1.PR.y1 = infl1.PR.y1, infl2.AR.y1 = infl2.AR.y1,
              infl2.PR.y1 = infl2.PR.y1))
}



## -----------------------------------------------------------------------------
## Function: risk.influences.calib()
## -----------------------------------------------------------------------------
## Description: This function estimates the influences on the absolute risk and 
##              pure risk for the primary event of interest. It should be used 
##              with calibrated weights. Also returns parameters estimates. The 
##              covariate profile and time or age interval for which the 
##              absolute risk and pure risk are to be estimated, are the ones 
##              that were provided when estimating the log-relative hazard, 
##              baseline hazard and cumulative baseline hazard 
## -----------------------------------------------------------------------------
## Arguments:
##
##  est.y1          return object from function influences, for the primary 
##                  event of interest
##
##  est.y2          return object from function influences, for the competing
##                  event
## ----------------------------------------------------------------------------

risk.influences.calib  <- function(est.y1, est.y2) {
  
  AR.y1       <- est.y1$exp.x * sum(est.y1$lambda0.t12.hat * est.y1$S.t12.hat * 
                                      est.y2$S.t12.hat)
  
  PR.y1       <- 1- exp(- est.y1$exp.x * sum(est.y1$lambda0.t12.hat))
  
  
  infl1.AR.y1  <- as.matrix(rowSums(est.y1$exp.x * 
                                      sweep(((est.y1$infl1.beta %*% 
                                                est.y1$x %*% 
                                                est.y1$lambda0.t12.hat) + 
                                               est.y1$infl1.lambda0.t12 - 
                                               sweep(((est.y1$exp.x * 
                                                         (est.y1$infl1.beta %*% 
                                                            est.y1$x %*% 
                                                            t(est.y1$Lambda0.t12.hat) + 
                                                            est.y1$infl1.Lambda0.t12)) + 
                                                        (est.y2$exp.x * 
                                                           (est.y2$infl1.beta %*% 
                                                              est.y2$x %*% 
                                                              t(est.y2$Lambda0.t12.hat) + 
                                                              est.y2$infl1.Lambda0.t12))), 2,
                                                     est.y1$lambda0.t12.hat, "*")), 2, 
                                            (est.y1$S.t12.hat * est.y2$S.t12.hat), "*" )))
  
  infl1.PR.y1    <- ((1 - PR.y1) * est.y1$exp.x * sum(est.y1$lambda0.t12.hat)) * 
    (est.y1$infl1.beta %*% est.y1$x) + 
    ((1 - PR.y1) * est.y1$exp.x) * rowSums(est.y1$infl1.lambda0.t12)
  
  infl2.AR.y1  <- as.matrix(rowSums(est.y1$exp.x * 
                                      sweep(((est.y1$infl2.beta %*% 
                                                est.y1$x %*% 
                                                est.y1$lambda0.t12.hat) + 
                                               est.y1$infl2.lambda0.t12 - 
                                               sweep(((est.y1$exp.x * 
                                                         (est.y1$infl2.beta %*% 
                                                            est.y1$x %*% 
                                                            t(est.y1$Lambda0.t12.hat) + 
                                                            est.y1$infl2.Lambda0.t12)) + 
                                                        (est.y2$exp.x * 
                                                           (est.y2$infl2.beta %*% 
                                                              est.y2$x %*% 
                                                              t(est.y2$Lambda0.t12.hat) + 
                                                              est.y2$infl2.Lambda0.t12))), 2,
                                                     est.y1$lambda0.t12.hat, "*")), 2, 
                                            (est.y1$S.t12.hat * est.y2$S.t12.hat), "*" )))
  
  infl2.PR.y1    <- ((1 - PR.y1) * est.y1$exp.x * sum(est.y1$lambda0.t12.hat)) * 
    (est.y1$infl2.beta %*% est.y1$x) + 
    ((1 - PR.y1) * est.y1$exp.x) * rowSums(est.y1$infl2.lambda0.t12)
  
  return(list(AR.y1 = AR.y1, PR.y1 = PR.y1, infl1.AR.y1 = infl1.AR.y1,
              infl1.PR.y1 = infl1.PR.y1, infl2.AR.y1 = infl2.AR.y1,
              infl2.PR.y1 = infl2.PR.y1))
}



## -----------------------------------------------------------------------------
## Function: robustvariance()
## -----------------------------------------------------------------------------
## Description: This function returns the robust variance estimate, i.e., 
##              computes the sum of the squared influences. The overall 
##              influences should be provided
## -----------------------------------------------------------------------------
## Argument:
##
##  infl        overall influences on parameters such as log-relative hazard, 
##              covariate specific pure risk or covariate specific absolute risk 
## -----------------------------------------------------------------------------

robustvariance <- function (infl) {
  
  return(robust.var = diag(crossprod(infl)))
  
}



## -----------------------------------------------------------------------------
## Function: variance()
## -----------------------------------------------------------------------------
## Description: This function returns the influence-based variance estimate that 
##              follows the complete variance decomposition. It should be used
##              with design weights under the case-cohort design with exhaustive 
##              sampling of cases. The overall influences should be provided
## -----------------------------------------------------------------------------
## Arguments:
##
## casecohort       data frame with id (identification) and weights (design 
##                  weights) for each individual in the case cohort data. It 
##                  should also contain status (case status), m (subcohort size)
##                  and n (cohort size)
##
## infl             overall influences on a parameter such as log-relative 
##                  hazard, cumulative baseline hazard or covariate specific 
##                  pure-risk 
## -----------------------------------------------------------------------------

variance <- function(casecohort, infl) {
  
  if (is.null(casecohort$weights)) {
    print("the subjects weights must be provided")
    break
  }
  
  noncases          <- which(casecohort$status == 0)
  n.noncases        <- length(noncases)
  m                 <- casecohort$m[1]
  n                 <- casecohort$n[1]
  m.n               <- m / n
  n.m               <- n / m
  
  Product.sigma.w  <- matrix(1 - m * (n - 1) / 
                               (n * (m - 1)), 
                             nrow = n.noncases, ncol = n.noncases)
  diag(Product.sigma.w) <- 1 - m.n
  
  phase1.var  <- n / (n - 1) * crossprod((infl / sqrt(casecohort$weights)))
  phase2.var  <- with(casecohort, t(infl[status == 0,]) %*% 
                        Product.sigma.w %*% infl[status == 0,])
  var         <- diag(phase1.var + phase2.var)
  
  return(variance = var)
  
}



## -----------------------------------------------------------------------------
## Function: variance.calib()
## -----------------------------------------------------------------------------
## Description: This function returns the influence-based variance estimate that 
##              follows the complete variance decomposition. It should be used 
##              with calibrated weights under the case-cohort design with 
##              exhaustive sampling of cases
## -----------------------------------------------------------------------------
## Arguments:
##
## casecohort       data frame with id (identification) and weights (design 
##                  weights) for each individual in the case cohort data. It 
##                  should also contain status (case status), m (subcohort size)
##                  and n (cohort size)
##
## cohort           data frame with status (case status) and subcohort 
##                  (subcohort sampling indicators) for each individual
##
## infl1            phase-one influences on a parameter such as log-relative 
##                  hazard, cumulative baseline hazard or covariate specific 
##                  pure-risk 
##
## infl2            phase-two influences
## -----------------------------------------------------------------------------

variance.calib <- function(casecohort, infl1, infl2, cohort, 
                           variance.phase2 = NULL) {
  
  if (is.null(variance.phase2)) {
    variance.phase2 <- FALSE
  }
  
  noncases          <- which(casecohort$status == 0)
  n.noncases        <- length(noncases)
  m                 <- casecohort$m[1]
  n                 <- casecohort$n[1]
  m.n               <- m / n
  n.m               <- n / m
  
  product.sigma  <- matrix(1 - m * (n - 1) / (n * (m - 1)), 
                             nrow = n.noncases, ncol = n.noncases)
  diag(product.sigma) <- 1 - m.n
  
  product.sigma[which(is.na(product.sigma))] <- 0
  
  omega                 <- rep(1, n)
  omega[casecohort$id]  <- casecohort$weights
  
  phase1.var  <- n / (n - 1) * colSums((infl1) ^ 2 + 2 * 
                                           (infl1) * infl2 + 
                                           (infl2 / sqrt(omega)) ^ 2)
  
  phase2.var    <- with(cohort, 
                        t(infl2[(status == 0) & (subcohort == TRUE),]) %*% 
                          product.sigma %*% 
                          infl2[(status == 0) & (subcohort == TRUE),])
  var           <- phase1.var + diag(phase2.var)
  
  if (variance.phase2 == FALSE) {
    return(variance = var)
  } else {
    return(list(variance = var, variance.phase2 = diag(phase2.var)))
  }
}



## -----------------------------------------------------------------------------
## Function: variance.generalized()
## -----------------------------------------------------------------------------
## Description: This function returns the influence-based variance estimate that 
##              follows the complete variance decomposition. It should be used 
##              with design weights or calibrated weights under the case-cohort 
##              design with stratified sampling based on case status
## -----------------------------------------------------------------------------
## Arguments:
##
## n                number of individuals in the whole cohort
##
## casecohort       data frame with id (identification) and weights (design 
##                  weights) for each individual in the case cohort data. It 
##                  should also contain status (case status, used as strata), 
##                  strata.m (number of sampled individuals in the stratum) and 
##                  strata.n (stratum size)
##
## cohort           data frame with status (case status, used as strata), id
##                  (identification), and phase2 (indicator of being in the 
##                  case-cohort, i.e., phase-two sample) for each individual
##
## variance.phase2  should the phase-two component of the variance be returned?
##
## Breslow.all      was only a fraction of the cases included in the case-cohort 
##                  and were the event times of all the cases in the cohort used
##                  in the numerator of the Breslow estimator? If 
##                  Breslow.all = TRUE, arguments infl1 and infl2 below need to 
##                  be provided. If Breslow.all = FALSE and calibrated = FALSE, 
##                  the last argument below needs to be provided. Default is 
##                  TRUE
##
## calibrated       were calibrated weights used for estimation? If 
##                  calibrated = TRUE, arguments infl1 and infl2 below need to 
##                  be provided. If Breslow.all = FALSE and calibrated = FALSE, 
##                  the last argument below needs to be provided
##
## infl1            phase-one influences on a parameter such as log-relative 
##                  hazard, cumulative baseline hazard or covariate specific 
##                  pure-risk. Should be provided if Breslow.all = TRUE or if
##                  calibrated = TRUE
##
## infl2            phase-two influences. Should be provided if 
##                  Breslow.all = TRUE or if calibrated = TRUE
##
## infl             overall influences on a parameter such as log-relative 
##                  hazard, cumulative baseline hazard or covariate specific 
##                  pure-risk. Should be provided if Breslow.all = FALSE and 
##                  calibrated = FALSE
## -----------------------------------------------------------------------------

variance.generalized <- function(casecohort, infl1 = NULL, infl2 = NULL,
                                 infl = NULL, cohort, Breslow.all = NULL,
                                 calibrated = NULL,
                                 variance.phase2 = NULL) {
  
  if (is.null(Breslow.all)) {
    Breslow.all <- FALSE
  }
  
  if (is.null(calibrated)) {
    calibrated <- FALSE
  }
  
  if (is.null(variance.phase2)) {
    variance.phase2 <- FALSE
  }
  
  n                     <- nrow(cohort)
  ncasecohort           <- nrow(casecohort)
  
  casecohort.W          <- casecohort$status
  casecohort.m          <- casecohort$m.strata
  casecohort.n          <- casecohort$n.strata
  casecohort.m.n        <- casecohort.m / casecohort.n
  casecohort.n.m        <- casecohort.n / casecohort.m
  
  product.sigma.w <- function (indiv) { 
    return((casecohort.W == casecohort.W[indiv]) * (1 - casecohort.m.n * 
             (casecohort.n - 1) / (casecohort.m - 1)))
  }
  Product.sigma.w                      <- do.call(rbind, lapply(1:ncasecohort, 
                                                                product.sigma.w))
  diag(Product.sigma.w)                <- 1 - casecohort.m.n
  Product.sigma.w[which(is.na(Product.sigma.w))] <- 0
  
  if (Breslow.all == TRUE | calibrated == TRUE) {
    
    omega                 <- rep(1, n)
    names(omega)          <- cohort$id
    omega[casecohort$id]  <- casecohort$weights
    
    phase1.var  <- n / (n - 1) * colSums((infl1) ^ 2 + 2 * 
                                             (infl1) * infl2 + 
                                             (infl2 / sqrt(omega)) ^ 2)
    
    phase2.var    <- with(cohort, t(infl2[(phase2 == TRUE),]) %*% 
                            Product.sigma.w %*% infl2[(phase2 == TRUE),])
    
    var           <- phase1.var + diag(phase2.var)
    
  } else {
    
    phase1.var  <- n / (n - 1) * crossprod((infl / sqrt(casecohort$weights)))
    
    phase2.var    <- with(casecohort, t(infl) %*% Product.sigma.w %*% infl)

    var           <- diag(phase1.var + phase2.var)
  }

  if (variance.phase2 == FALSE) {
    return(variance = var)
  } else {
    return(list(variance = var, variance.phase2 = diag(phase2.var)))
  }
}



## -----------------------------------------------------------------------------
## Function: conf.interval()
## -----------------------------------------------------------------------------
## Description: This function computes the 95% confidence interval of a 
##              parameter estimate and its variance estimate, assuming normality
## -----------------------------------------------------------------------------
## Arguments:
##
##  param.estimate  parameter estimate
##                    
##  var.estimate    variance estimate for the parameter
## -----------------------------------------------------------------------------

conf.interval <- function (param.estimate, var.estimate) {
  
  return(c(param.estimate - sqrt(var.estimate) * qnorm(0.975), 
           param.estimate + sqrt(var.estimate) * qnorm(0.975)))
  
}



## -----------------------------------------------------------------------------
## Function: calibration()
## -----------------------------------------------------------------------------
## Description: This function calibrates the design weights using the raking 
##              procedure. It matches the weighted totals of the auxiliary 
##              variables in the stratified case-cohort (with calibrated 
##              weights), to the un-weighted auxiliary variables totals in the 
##              whole cohort. The Newton Raphson method is used to solve the 
##              optimization problem
## -----------------------------------------------------------------------------
## Arguments:
##
##  A               matrix with the values of the auxiliary variables to be used 
##                  for the calibration of the weights in the case-cohort 
##                  (phase-two data)
##
##  design.weights  design weights to be calibrated
##
##  total           un-weighted auxiliary variables totals in the whole cohort
##
##  eta0            vector of initial/possible values for eta, to be used as 
##                  seed in the iterative procedure. Default is (0,...,0)
##
##  niter.max       maximum number of iterations for the Newton Raphson method. 
##                  Default is 10^4 iterations
##
##  epsilon.stop    threshold for the difference between the estimated weighted 
##                  total and the total in the whole cohort. If this difference 
##                  is less than the value of epsilon.stop, no more iterations 
##                  will be performed. Default is 10^(-10)
## -----------------------------------------------------------------------------

calibration <- function (A.phase2, design.weights, total, eta0 = NULL, 
                         niter.max = NULL, epsilon.stop = NULL) {
  
  A.phase2            <- as.matrix(A.phase2)
  
  n                   <- nrow(A.phase2)
  q                   <- ncol(A.phase2)
  
  if ((is.null(eta0)) || (length(eta0) != q)) {
    eta0          <- rep(0, q)
    print("eta0 has not been filled out. Default is 0. ")
  }
  if ((is.null(epsilon.stop)) || (epsilon.stop <= 0)) {
    epsilon.stop  <- 10^(-10)
    print("epsilon.stop has not been filled out. Default is 10^(-10). ")
  }
  if ((is.null(niter.max)) || (niter.max <= 0)) {
    niter.max     <- 10^4
    print("niter.max has not been filled out. Default is 10^4 iterations.")
  }
  
  for (niter in 1:niter.max) {
    
    calibrated.weights  <- design.weights * exp(A.phase2 %*% eta0)
    A.calweighted       <- sweep(A.phase2, 1, calibrated.weights, "*")
    f                   <- colSums(A.calweighted) - total
    
    epsilon             <- max(abs(f))
    if(epsilon < epsilon.stop){
      break
    }
    
    AA.calweighted      <- array(NA, dim = c(q, q, n))
    for (i in 1:n) {
      AA.calweighted[,, i] <- tcrossprod(A.phase2[i,], A.calweighted[i,])
    }
    drond.f.eta0        <- apply(X = AA.calweighted, MARGIN = c(1,2), FUN = sum)
    
    if(is.numeric(try(solve(drond.f.eta0), silent = TRUE))){
      eta0  <- eta0 - c(f %*% solve(drond.f.eta0))
    }else{
      eta0  <- NA
    }
    
    if(is.na(sum(eta0)) == TRUE){
      break
    }
  }
  
  if(is.na(sum(eta0)) == TRUE){
    calibrated.weights    <- NA
    A.calweighted         <- NA
    estimated.total       <- NA
  }else{
    calibrated.weights    <- design.weights * exp(A.phase2 %*% eta0)
    A.calweighted         <- sweep(A.phase2, 1, calibrated.weights, "*")
    estimated.total       <- colSums(A.calweighted)
  }
  
  return(list(eta.hat = eta0, calibrated.weights = calibrated.weights, 
              estimated.total = estimated.total))
}



## -----------------------------------------------------------------------------
## Function: X.generation()
## -----------------------------------------------------------------------------
## Description:   This function generates a gaussian covariate value, with unit
##                variance and given mean
## -----------------------------------------------------------------------------
## Arguments:
##
##  mu                    mean of the covariate
## -----------------------------------------------------------------------------

X.generation <- function (mu) { return(rnorm(n = 1, mean = mu, sd = 1)) }

