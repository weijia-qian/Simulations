### Calculate cross-validated c-index and Brier's score
# cross validated c-index and Brier's score for each model
## introduce variables in the cross validation
set.seed(135) 
n_folds <- 10 # number of folds
### variables that store the c-index and Brier's score in each calculation
AUC_norm <- AUC_cox <- AUC_sieve <- rep(0, n_folds)
Brier_norm <- Brier_cox <- Brier_sieve <- rep(0, n_folds)
### use "createFolds" function from "caret" package to divide real data into n_folds folds
fold <- createFolds(sim_data$data$Y, k = n_folds, list = FALSE)

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
    inx_i <- which(Stime == ti & status==1)
    ## subjects with observed times beyond event current time
    inx_j <- which(Stime > ti)
    ## number of "cases" and "controls" at time current time
    n_case_t    <- length(inx_i)
    n_control_t <- length(inx_j)
    for(i in seq_along(inx_i)){
      num   <- num + sum( (marker[inx_j] > marker[inx_i[i]] ) ) + 0.5*sum(marker[inx_j] == marker[inx_i[i]])
    }
    denom <- denom + n_case_t*n_control_t
  }
  1-num/denom
}

### function "cal_stime()" returns the survival times estimates on the test data in cross-validated prediction
cal_stime <- function(fit, tgrid_sim = seq(0, 10, len=1000), family = 'cox.ph'){
  if (family == 'cox.ph') {
    ## calculate estimated baseline hazard
    t0 <- rev(fit$family$data$tr)
    H0_hat <- rev(fit$family$data$h) ## estimate for the baseline cumulative hazard function \int_0^t \lambda_0(s)ds
    H0_fit <- gam(H0_hat ~ s(t0, pc = 0)-1) ## smooth the cumulative hazard
    ## evaluate smoothed H0 on fine grid
    H0_prd <- pmax(0, predict(H0_fit, newdata = data.frame(t0 = tgrid_sim))) ## truncate cumulative hazard at 0
    ## estimate eta from estimated functional surface of real data
    eta_i <- matrix(predict(fit, test_data), ncol = 1)
    ## calcualte survival times
    S_i <- exp(-(exp(eta_i) %*% H0_prd))
  } else if (family == 'lognormal'){
    lp <- fit$linear.predictors
    scale <- fit$scale
    S_i <- outer(lp, tgrid_sim, function(lp_i, tgrid_sim_j) pnorm((log(tgrid_sim_j) - lp_i) / scale, lower.tail = FALSE))
  } else if (family == 'loglogistic'){
    lp <- fit$linear.predictors
    scale <- fit$scale
    S_i <- outer(lp, tgrid_sim, function(lp_i, tgrid_sim_j) 1 - 1 / (1 + (exp(lp_i) / tgrid_sim_j)^(1 / scale)))
  }
  return(S_i)
}

### function "cal_Brier()" calculates Brier's score under right censored survival setting
cal_Brier <- function(S, Stime, status, tgrid = seq(0, 10, len=1000)){
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

## run cross validation
for(k in 1:n_folds){
  ## split train and test data
  train_data <- sim_data$data[-which(fold == k),]
  test_data <- sim_data$data[which(fold == k),]
  
  ## survival time and status
  time_test   <- test_data$Y
  event_test <- test_data$delta
  time_train   <- train_data$Y
  event_train <- train_data$delta
  
  ## get unique ordered survival times
  ut_train <- unique(time_train[event_train == 1])
  ut_train <- ut_train[order(ut_train)]
  ut_test <- unique(time_test[event_test == 1])
  ut_test <- ut_test[order(ut_test)]
  
  ## get KM estimates of survival function for train data
  KM_train <- unique(survfit(Surv(time_train,event_train)~1)$surv)
  ## get KM estimates of survival function for test data where test and training overlap
  ## impute linearly on the log scale for survival times which are not in the training dataset
  KM_test  <- rep(NA, length(ut_test))
  for(i in seq_along(ut_test)){
    ## if test data survival time less than the minimum observed event time in the training dataset
    ## impute survival time on the linearly log scale
    if(ut_test[i] < min(ut_train)){
      KM_test[i] <- exp(log(1) - (log(1) - log(KM_train[1])) / (ut_train[1]-0) * (ut_test[1]-0))
    }
    ## if test data survival time within the observed range of event times in the training dataset,
    ## either use the KM estimate (if the training time is in the test times), or impute linearly on the log scale
    if(ut_test[i] >= min(ut_train) & ut_test[i] < max(ut_train)){
      inx_l <- max(which(ut_train <= ut_test[i]))
      inx_r <- min(which(ut_train > ut_test[i]))
      st_l <- KM_train[inx_l]
      st_r <- KM_train[inx_r]
      t_l <- ut_train[inx_l]
      t_r <- ut_train[inx_r]
      KM_test[i] <- exp(log(st_l) - (log(st_l) - log(st_r))/(t_r -t_l)  * (ut_test[i] - t_l))
    }
    ## if the test data survival time is beyond the observe range of event times in the training dataset,
    ## use the last observed survival probability
    if(ut_test[i] >= max(ut_train)){
      KM_test[i] <- min(KM_train)
    }
  }
  
  ## derive the KM estimate of the censoring time
  ut_train_censor <- unique(time_train[event_train == 0])
  ut_train_censor <- ut_train_censor[order(ut_train_censor)]
  ut_test_censor <- unique(time_test[event_test == 0])
  ut_test_censor <- ut_test_censor[order(ut_test_censor)]
  ## get KM estimates of censoring time for train data
  KM_train_censor <- unique(survfit(Surv(time_train,1-event_train)~1)$surv)
  ## get KM estimates of censoring time for test data where test and training overlap
  ## impute linearly on the log scale for survival times which are not in the training dataset
  KM_test_censor  <- rep(NA, length(ut_test_censor))
  for(i in seq_along(ut_test_censor)){
    ## if test data survival time less than the minimum observed event time in the training dataset
    ## impute survival time on the linearly log scale
    if(ut_test_censor[i] < min(ut_train_censor)){
      KM_test_censor[i] <- exp(log(1) - (log(1) - log(KM_train_censor[1])) / (ut_train_censor[1]-0) * (ut_test_censor[1]-0))
    }
    ## if test data survival time within the observed range of event times in the training dataset,
    ## either use the KM estimate (if the training time is in the test times), or impute linearly on the log scale
    if(ut_test_censor[i] >= min(ut_train_censor) & ut_test_censor[i] < max(ut_train_censor)){
      inx_l <- max(which(ut_train_censor <= ut_test_censor[i]))
      inx_r <- min(which(ut_train_censor > ut_test_censor[i]))
      st_l <- KM_train_censor[inx_l]
      st_r <- KM_train_censor[inx_r]
      t_l <- ut_train_censor[inx_l]
      t_r <- ut_train_censor[inx_r]
      KM_test_censor[i] <- exp(log(st_l) - (log(st_l) - log(st_r))/(t_r -t_l)  * (ut_test_censor[i] - t_l))
    }
    ## if the test data survival time is beyond the observe range of event times in the training dataset,
    ## use the last observed survival probability
    if(ut_test_censor[i] >= max(ut_train_censor)){
      KM_test_censor[i] <- min(KM_train_censor)
    }
  }
  
  
  ## obtain estimated log hazard (excluding baseline hazard)
  eta_norm <- -fit.norm$linear.predictors[which(fold == k)]
  eta_cox <- rowSums(predict(fit.cox, test_data, type="terms"))
  
  ## derive cross-validated c-index
  AUC_norm[k] <- cal_c(marker = eta_norm, Stime = time_test, status = event_test)
  AUC_cox[k] <- cal_c(marker = eta_cox, Stime = time_test, status = event_test)

  ## calculate the cross-validated brier score
  S_norm <- cal_stime(fit = fit.norm, family = 'lognormal')
  S_cox <- cal_stime(fit = fit.cox, family = 'cox.ph')
  
  tmax <- round(max(sim_data$data$t))
  tvec <- seq(0, tmax, length.out = 1000)
  
  Brier_norm[k] <- cal_Brier(S_norm, Stime = time_test, status = event_test, tgrid = tvec)
  Brier_cox[k] <- cal_Brier(S_cox, Stime = time_test, status = event_test, tgrid = tvec)

}
