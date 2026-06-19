### function "cal_c()" calculates Concordance using concordant vs discordant pairs empirical estimator
## marker = \hat\eta_i
## Stime = survival times
## status = event status (0 = censored, 1 = event)
cal_c <- function(marker, Stime, status){
  utimes <- sort(unique(Stime[status==1]))
  num <- denom <- 0
  for(ut in seq_along(utimes)){
    ## current time
    ti    <- utimes[ut]
    ## subjects who experienced an event at current time
    inx_i <- which(Stime == ti & status == 1)
    ## subjects with observed times beyond event current time
    inx_j <- which(Stime > ti)
    ## number of "cases" and "controls" at time current time
    n_case_t    <- length(inx_i)
    n_control_t <- length(inx_j)
    for(i in seq_along(inx_i)){
      num <- num + sum((marker[inx_j] > marker[inx_i[i]])) + 0.5 * sum(marker[inx_j] == marker[inx_i[i]])
    }
    denom <- denom + n_case_t * n_control_t
  }
  1 - num / denom
}

### function "cal_stime()" returns the survival times estimates
cal_stime <- function(fit, data, tgrid = seq(0, 10, len = 1000), family = "cox.ph"){
  if (family == "cox.ph") {
    n <- nrow(data)
    nS_pred <- length(tgrid)
    df_pred <- data[rep(1:n, each = nS_pred), ]
    df_pred$Y <- rep(tgrid, n)
    S_i <- matrix(predict(fit, newdata = df_pred, type = "response"), nrow = n, ncol = nS_pred, byrow = TRUE)
    # ## calculate estimated baseline hazard
    # t0 <- rev(fit$family$data$tr) # observed event times
    # H0_hat <- rev(fit$family$data$h) # Breslow estimator of the cumulative hazard function
    # H0_fit <- scam(log(H0_hat + 1e-8) ~ s(t0, bs = "mpi") - 1) # smooth while imposing non-decreasing shape constraints
    # #H0_fit <- scam(H0_hat ~ s(t0, bs = "mpi") - 1)
    # ## evaluate smoothed H0 on fine grid
    # H0_prd <- exp(predict(H0_fit, newdata = data.frame(t0 = tgrid)))
    # #H0_prd <- pmax(0, predict(H0_fit, newdata = data.frame(t0 = tgrid)))
    # ## estimate eta from estimated functional surface of real data
    # eta_i <- matrix(predict(fit, data), ncol = 1)
    # ## calculate survival times
    # S_i <- exp(-(exp(eta_i) %*% H0_prd))
  } else if (family == "lognormal"){
    lp <- predict(fit, data, type = "response")
    scale <- as.numeric(gsub(".*\\(([^)]+)\\).*", "\\1", fit$family$family))
    S_i <- outer(lp, tgrid, function(lp_i, tgrid_j) pnorm((log(tgrid_j) - lp_i) / scale, lower.tail = FALSE))
  }
  return(S_i)
}

### function "choose_ub_unif()" finds the upper bound of Uniform(0, ub) censoring
# that achieves a target censoring rate.
# P(censored) = P(C < T) where C ~ Uniform(0, ub).
# For a given set of true event times T_true, solves for ub such that
#   mean(pmin(T_true, ub) / ub) = censor_rate
# using the identity P(C < T) = E[pmin(T, ub) / ub].
choose_ub_unif <- function(T_true, censor_rate, tol = 1e-8) {
  stopifnot(censor_rate > 0, censor_rate < 1)
  g <- function(ub) mean(pmin(T_true, ub) / ub) - censor_rate
  lo <- max(tol, min(T_true[T_true > 0], na.rm = TRUE) * 1e-6)
  hi <- max(T_true, na.rm = TRUE) * 2 + 1
  if (g(lo) < 0) lo <- tol
  while (g(hi) > 0) hi <- hi * 2
  uniroot(g, lower = lo, upper = hi)$root
}

### function "cal_IPCW_IBS()" calculates IPCW integrated Brier score
# Evaluates at observed event times in the test set.
# Integration: trapezoidal rule for first interval, right Riemann sum thereafter,
# normalized by max observed event time. Censoring distribution estimated by KM on training data.
# Arguments:
#   S_mat       : n_test x length(tgrid) survival probability matrix
#   time_test   : observed times in test set
#   event_test  : event indicators in test set (1 = event, 0 = censored)
#   time_train  : observed times in training set (for KM censoring estimate)
#   event_train : event indicators in training set
#   tgrid       : time grid on which S_mat is evaluated (should start > 0)
cal_IPCW_IBS <- function(S_mat, time_test, event_test, time_train, event_train, tgrid) {
  n_test  <- length(time_test)
  km_cens <- survfit(Surv(time_train, 1 - event_train) ~ 1)
  get_G <- function(t_vals) {
    sf     <- stepfun(km_cens$time, c(1, km_cens$surv))
    G_vals <- sf(t_vals)
    G_vals[G_vals == 0] <- min(km_cens$surv[km_cens$surv > 0], na.rm = TRUE)
    G_vals
  }
  utimes <- sort(unique(time_test[event_test == 1]))
  if (length(utimes) == 0) return(NA_real_)
  Brier_t <- vapply(utimes, function(ti) {
    col_idx <- max(which(tgrid <= ti))
    S_ti    <- S_mat[, col_idx]
    idx1    <- (time_test <= ti) & (event_test == 1)
    idx2    <- time_test > ti
    (sum(S_ti[idx1]^2 / get_G(time_test[idx1])) +
       sum((1 - S_ti[idx2])^2 / get_G(ti))) / n_test
  }, numeric(1))
  iBrier <- 0.5 * utimes[1] * Brier_t[1]
  if (length(utimes) > 1)
    iBrier <- iBrier + sum(Brier_t[-1] * diff(utimes))
  iBrier / max(utimes)
}

### function "cal_Brier()" calculates Brier's score under right censored survival setting
cal_Brier <- function(S, Stime, status, tgrid = seq(0, 10, len=1000)){
  # get unique ordered survival times
  ut_train <- unique(time_train[event_train == 1])
  ut_train <- ut_train[order(ut_train)]
  ut_test <- unique(time_test[event_test == 1])
  ut_test <- ut_test[order(ut_test)]
  
  # derive the KM estimate of the censoring time
  ut_train_censor <- unique(time_train[event_train == 0])
  ut_train_censor <- ut_train_censor[order(ut_train_censor)]
  ut_test_censor <- unique(time_test[event_test == 0])
  ut_test_censor <- ut_test_censor[order(ut_test_censor)]
  # get KM estimates of censoring time for train data
  KM_train_censor <- unique(survfit(Surv(time_train,1-event_train) ~ 1)$surv)
  # get KM estimates of censoring time for test data where test and training overlap
  # impute linearly on the log scale for survival times which are not in the training dataset
  KM_test_censor  <- rep(NA, length(ut_test_censor))
  for(i in seq_along(ut_test_censor)){
    # if test data survival time less than the minimum observed event time in the training dataset
    # impute survival time on the linearly log scale
    if(ut_test_censor[i] < min(ut_train_censor)){
      KM_test_censor[i] <- exp(log(1) - (log(1) - log(KM_train_censor[1])) / (ut_train_censor[1]-0) * (ut_test_censor[1]-0))
    }
    # if test data survival time within the observed range of event times in the training dataset,
    # either use the KM estimate (if the training time is in the test times), or impute linearly on the log scale
    if(ut_test_censor[i] >= min(ut_train_censor) & ut_test_censor[i] < max(ut_train_censor)){
      inx_l <- max(which(ut_train_censor <= ut_test_censor[i]))
      inx_r <- min(which(ut_train_censor > ut_test_censor[i]))
      st_l <- KM_train_censor[inx_l]
      st_r <- KM_train_censor[inx_r]
      t_l <- ut_train_censor[inx_l]
      t_r <- ut_train_censor[inx_r]
      KM_test_censor[i] <- exp(log(st_l) - (log(st_l) - log(st_r))/(t_r -t_l)  * (ut_test_censor[i] - t_l))
    }
    # if the test data survival time is beyond the observe range of event times in the training dataset,
    # use the last observed survival probability
    if(ut_test_censor[i] >= max(ut_train_censor)){
      KM_test_censor[i] <- min(KM_train_censor)
    }
  }
  
  iBrier <- 0
  utimes <- sort(unique(Stime[status == 1]))
  if(max(utimes) == 10){
    utimes <- utimes[-length(utimes)]
  }
  Brier.t <- rep(NA, length(utimes))
  for(ut in seq_along(utimes)){
    ## current time
    ti <- utimes[ut]
    ## subjects with observed times beyond event current time
    inx_j <- which(Stime > ti)
    ## subjects who experienced an event before or at current time
    inx_i <- which(Stime <= ti & status==1)
    ## extract survival probability of all subjects at this time
    surv.prob <- S[,max(which(tgrid <= ti))]
    tmp <- 0
    for(i in 1:nrow(S)){
      G <- 1
      if(i %in% inx_j){
        if(max(which(ut_test_censor < ti)) > 0){
          G <- KM_test_censor[max(which(ut_test_censor < ti))]
        }
        tmp <- tmp + (1-surv.prob[i])^2/G
      }else if(i %in% inx_i){
        if(max(which(ut_test_censor < Stime[i])) > 0){
          G <- KM_test_censor[max(which(ut_test_censor < Stime[i]))]
        }
        tmp <- tmp + (0-surv.prob[i])^2/G
      }
    }
    Brier.t[ut] <- tmp/nrow(S)
    if(ut == 1){
      iBrier <- iBrier + 0.5*utimes[ut]*Brier.t[ut]
    }else{
      iBrier <- iBrier + Brier.t[ut]*(utimes[ut] - utimes[ut-1])
    }
  }
  iBrier <- iBrier/max(utimes)
  return(iBrier)
}
