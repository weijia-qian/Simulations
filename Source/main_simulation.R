####################################################################
# Weijia Qian
# November 2024
#
# This file simulates survival data under different data generation mechanisms
# and fits functional AFT and Cox models
####################################################################

#suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(fda))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(refund))
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(library(tidyverse))

wd = getwd()

if(substring(wd, 2, 6) == "Users"){
  doLocal = TRUE
}else{
  doLocal = FALSE
}

###############################################################
## define or source functions used in code below
###############################################################
source(here("Source", "simulate_AFT.R"))
source(here("Source", "simulate_Cox.R"))
source(here("Source", "est_sieve.R"))
source(here("Source", "calc_auc_brier.R"))
source(here("Source", "utils_summary.R"))
source(here("Source", "fit_penalized_aft.R"))
source(here("Software", "sourceFTTM.R"))
load(here("Source", "dat_func.Rdata")) # load real data

###############################################################
## set simulation design elements
###############################################################
family = c("lognormal", "loglogistic", "cox.ph")
n = c(100, 200, 500)
#nS = c(50, 100, 500)
nS = c(100)
#beta_type = c('simple', 'complex')
beta_type = c('complex')
b = c(0.1)
seed_start = 1000
N_iter = 500

params = expand.grid(seed_start = seed_start,
                     family = family,
                     n = n,
                     nS = nS,
                     beta_type = beta_type,
                     b = b)

## define number of simulations and parameter scenarios
if(doLocal) {
  scenario = 1
  N_iter = 30
}else{
  # defined from batch script params
  scenario <- as.numeric(commandArgs(trailingOnly=TRUE))
}

n = params$n[scenario]
family = params$family[scenario]
nS = params$nS[scenario]
beta_type = params$beta_type[scenario]
b = params$b[scenario]
SEED.START = params$seed_start[scenario]

###############################################################
## run simulations
###############################################################
results = vector("list", length = N_iter)

# simulate a test dataset
if(family %in% c("lognormal", "loglogistic")){
  sim_data_test <- simulate_AFT(family = family, n = n, nS = nS, beta_type = beta_type, b = b, seed = SEED.START)
}else{
  sim_data_test <- simulate_Cox(n = n, nS = nS, beta_type = beta_type, seed = SEED.START)
}

for(iter in 1:N_iter){
  print(iter)
  # set seed
  seed.iter = (SEED.START - 1) * N_iter + iter

  # simulate data
  if(family %in% c("lognormal", "loglogistic")){
    sim_data <- simulate_AFT(family = family, n = n, nS = nS, beta_type = beta_type, b = b, seed = seed.iter)
  }else{
    sim_data <- simulate_Cox(n = n, nS = nS, beta_type = beta_type, seed = seed.iter)
  }

  ###############################################################
  ## fit functional AFT and Cox model
  ###############################################################
  # linear functional log-normal AFT model
  time.norm <- as.numeric(system.time(fit.norm <- gam(logY ~ 1 + s(S, by = X_L, bs = "ps", k = 30), 
                                                      family = cnorm(), data = sim_data$data))[3])
  # linear functional Cox model
  time.cox <- as.numeric(system.time(fit.cox <- gam(Y ~ s(S, by = X_L, bs = "ps", k = 30), data = sim_data$data, 
                                                        weights = delta, family = cox.ph))[3])
  
  # my function
  tic()
  lambda_grid <- exp(seq(log(1000), log(10000), length.out = 500))
  model <- "lognormal"
  best_lambda <- optimize_lambda(Y, delta, X, data = sim_data$data, family = model, lambda_grid)
  fit.faft <- optimize_AFT(Y, delta, X, data = sim_data$data, family = model, lambda = best_lambda, se = TRUE)
  time_stamp <- toc(quiet = TRUE)
  time_faft <- time_stamp$toc - time_stamp$tic
  
  
  ###############################################################
  ## Out-of-sample Harrell’s C-index and Brier score
  ###############################################################
  
  # survival time and status
  time_train <- sim_data$data$Y
  event_train <- sim_data$data$delta
  time_test <- sim_data_test$data$Y
  event_test <- sim_data_test$data$delta

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

  ## obtain linear predictors
  eta_norm <- predict(fit.norm, sim_data_test$data, type = "response")
  eta_cox <- rowSums(predict(fit.cox, sim_data_test$data, type = "terms"))
  eta_faft <- fit.faft$lp
  
  ## calculate c-index
  AUC_norm <- cal_c(marker = -eta_norm, Stime = time_test, status = event_test)
  AUC_cox <- cal_c(marker = eta_cox, Stime = time_test, status = event_test)
  AUC_faft <- cal_c(marker = -eta_faft, Stime = time_test, status = event_test)
  
  ## calculate the brier score
  tmax_test <- round(quantile(sim_data_test$data$t, 0.99))
  tgrid_test <- seq(0, tmax_test, length.out = 1000)
  
  S_norm <- cal_stime(fit = fit.norm, data = sim_data_test$data, tgrid = tgrid_test, family = 'lognormal')
  S_cox <- cal_stime(fit = fit.cox, data = sim_data_test$data, tgrid = tgrid_test, family = 'cox.ph')
  lp_faft <- as.numeric(fit.faft$lp)
  scale_faft <- fit.faft$b_hat
  S_faft <- outer(lp_faft, tgrid_test, function(lp_faft_i, tgrid_test_j) pnorm((log(tgrid_test_j) - lp_faft_i) / scale_faft, 
                                                                                 lower.tail = FALSE))
  
  Brier_norm <- cal_Brier(S_norm, Stime = time_test, status = event_test, tgrid = tgrid_test)
  Brier_cox <- cal_Brier(S_cox, Stime = time_test, status = event_test, tgrid = tgrid_test)
  Brier_faft <- cal_Brier(S_faft, Stime = time_test, status = event_test, tgrid = tgrid_test)
  
  ###############################################################
  ## pointwise squared errors, pointwise CIs and CMA CIs for estimated beta
  ###############################################################
  # calculate pointwise squared errors
  sgrid <- seq(0, 1, len = nS)
  df_pred <- data.frame(S = sgrid, X_L = 1)
  
  coef.true <- sim_data$coefficients$beta1
  coef.est.norm <- predict(fit.norm, newdata = df_pred, type = "terms", se.fit = TRUE)
  coef.est.cox <- predict(fit.cox, newdata = df_pred, type = "terms", se.fit = TRUE)
  
  if (family %in% c("lognormal", "loglogistic")) {
    se.coef.norm <- (coef.true - coef.est.norm[[1]])^2
    se.coef.cox <- (-coef.true - coef.est.cox[[1]])^2
    se.coef.faft <- (coef.true - fit.faft$beta1_hat)^2
  } else {
    se.coef.norm <- (-coef.true - coef.est.norm[[1]])^2
    se.coef.cox <- (coef.true - coef.est.cox[[1]])^2
    se.coef.faft <- (-coef.true - fit.faft$beta1_hat)^2
  }
  
  # calculate CMA CIs
  cma.coef.norm <- get_CMA(fit.norm)
  cma.coef.cox <- get_CMA(fit.cox)
  
  # fit AFT model using sieve algorithm
  #sieve.results <- est_sieve(data = sim_data$data)

  # summary
  df_coef <- data.frame(true_coef = coef.true,                           # true coefficient function
                        est_coef_norm = as.numeric(coef.est.norm[[1]]),  # estimated coefficient function
                        se_coef_norm = as.numeric(se.coef.norm),         # pointwise squared error
                        lb_coef_norm = as.numeric(coef.est.norm[[1]] - qnorm(0.975) * coef.est.norm[[2]]), # pointwise CI lower bound
                        ub_coef_norm = as.numeric(coef.est.norm[[1]] + qnorm(0.975) * coef.est.norm[[2]]), # pointwise CI upper bound
                        cma_lb_coef_norm = as.numeric(cma.coef.norm[[1]]), # CMA CI lower bound
                        cma_ub_coef_norm = as.numeric(cma.coef.norm[[2]]), # CMA CI upper bound
                        est_coef_cox = as.numeric(coef.est.cox[[1]]), 
                        se_coef_cox = as.numeric(se.coef.cox), 
                        #est_coef_sieve = as.numeric(sieve.results[[1]]),
                        #se_coef_sieve = as.numeric(sieve.results[[2]]),
                        lb_coef_cox = as.numeric(coef.est.cox[[1]] - qnorm(0.975) * coef.est.cox[[2]]), 
                        ub_coef_cox = as.numeric(coef.est.cox[[1]] + qnorm(0.975) * coef.est.cox[[2]]),
                        cma_lb_coef_cox = as.numeric(cma.coef.cox[[1]]),
                        cma_ub_coef_cox = as.numeric(cma.coef.cox[[2]]),
                        est_coef_faft = fit.faft$beta1_hat, 
                        se_coef_faft = se.coef.faft, 
                        lb_coef_faft = fit.faft$beta1_ci_lower, 
                        ub_coef_faft = fit.faft$beta1_ci_upper) %>%
    mutate(cover_coef_norm = (true_coef > lb_coef_norm) & (true_coef < ub_coef_norm),
           cover_coef_cox = (true_coef > lb_coef_cox) & (true_coef < ub_coef_cox),
           cover_coef_faft = (true_coef > lb_coef_faft) & (true_coef < ub_coef_faft),
           cover_cma_coef_norm = (true_coef > cma_lb_coef_norm) & (true_coef < cma_ub_coef_norm),
           cover_cma_coef_cox = (true_coef > cma_lb_coef_cox) & (true_coef < cma_ub_coef_cox))
  
  ###############################################################
  ## pointwise squared errors for survival function
  ###############################################################
  # set the time grid to evaluate survival function
  if (family == "cox.ph"){
    tmax <- 300 # consistent with 'tmax' in simulate_Cox()
    tgrid <- seq(0, tmax, len = 1000) 
  } else {
    tmax <- round(quantile(sim_data$data$t, 0.99))
    tgrid <- seq(0, tmax, len = 1000)
    lp = sim_data$data$lp
    scale = sim_data$coefficients$b[1]
  }
  
  # true survival function
  if (family == "lognormal") {
    S_true <- outer(lp, tgrid, function(lp_i, tgrid_j) pnorm((log(tgrid_j) - lp_i) / scale, lower.tail = FALSE))
  } else if (family == "loglogistic") {
    S_true <- outer(lp, tgrid, function(lp_i, tgrid_j) 1 - 1 / (1 + (exp(lp_i) / tgrid_j)^(1 / scale)))
  } else {
    S_true = sim_data$data$Si
  }
  
  # estimated survival function
  S_norm <- cal_stime(fit = fit.norm, data = sim_data$data, tgrid = tgrid, family = 'lognormal')
  S_cox <- cal_stime(fit = fit.cox, data = sim_data$data, tgrid = tgrid, family = 'cox.ph')
  S_faft <- outer(lp_faft, tgrid, function(lp_faft_i, tgrid_j) pnorm((log(tgrid_j) - lp_faft_i) / scale_faft, 
                                                                                 lower.tail = FALSE))
  
  # calculate pointwise squared error and MISE
  df_surv <- data.frame(time = tgrid,
                        #true_surv = 1 - p.true,
                        #est_surv_norm = 1 - p.est.norm,
                        #est_surv_cox = 1 - p.est.cox,
                        se_surv_norm = colMeans((S_true - S_norm)^2),
                        se_surv_cox = colMeans((S_true - S_cox)^2),
                        se_surv_faft = colMeans((S_true - S_faft)^2))
  
  df_info <- data.frame(scenario = scenario,
                        iter = iter,
                        seed = seed.iter,
                        family = family,
                        n = n,
                        nS = nS,
                        beta_type = beta_type,
                        b = b,
                        censor_rate = 1 - mean(sim_data$data$delta),
                        AUC_norm,
                        AUC_cox,
                        AUC_faft,
                        Brier_norm,
                        Brier_cox,
                        Brier_faft,
                        time_norm = time.norm,
                        #time_sieve = sieve.results[[3]],
                        time_cox = time.cox,
                        time_faft = time_faft)
  
  res <- list(info = df_info, coef = df_coef, surv = df_surv)

  results[[iter]] = res

} # end for loop

# record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())
dir.create(file.path(here::here("Output"), Date), showWarnings = FALSE)

filename = paste0(here::here("Output", Date), "/", scenario, ".RDA")
save(results,
     file = filename)

###############################################################
## end sim
###############################################################


