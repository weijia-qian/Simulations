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

### function "cal_stime()" returns the survival times estimates
cal_stime <- function(fit, data, tgrid = seq(0, 10, len = 1000), family = "cox.ph"){
  if (family == "cox.ph") {
    ## calculate estimated baseline hazard
    t0 <- rev(fit$family$data$tr) # observed event times
    H0_hat <- rev(fit$family$data$h) # Breslow estimator of the cumulative hazard function
    H0_fit <- scam(H0_hat ~ s(t0, bs = "mpi") - 1) # smooth while imposing non-decreasing shape constraints
    #H0_fit <- gam(H0_hat ~ s(t0, pc = 0)-1)
    ## evaluate smoothed H0 on fine grid
    H0_prd <- pmax(0, predict(H0_fit, newdata = data.frame(t0 = tgrid)))
    ## estimate eta from estimated functional surface of real data
    eta_i <- matrix(predict(fit, data), ncol = 1)
    ## calculate survival times
    S_i <- exp(-(exp(eta_i) %*% H0_prd))
  } else if (family == "lognormal"){
    lp <- predict(fit, data, type = "response")
    scale <- fit$scale
    S_i <- outer(lp, tgrid, function(lp_i, tgrid_j) pnorm((log(tgrid_j) - lp_i) / scale, lower.tail = FALSE))
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
