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
