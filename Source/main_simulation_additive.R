####################################################################
# Weijia Qian
# November 2024
#
# This file simulates survival data under different data generation mechanisms
# and fits functional AFT and Cox models
####################################################################

#suppressPackageStartupMessages(library(caret))
#suppressPackageStartupMessages(library(fda))
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
source(here("Source", "simulate_AFT2.R"))
source(here("Source", "calc_auc_brier.R"))
#source(here("Source", "utils_summary.R"))
#source(here("Source", "fit_penalized_aft.R"))
#source(here("Software", "sourceFTTM.R"))
load(here("Source", "dat_func.Rdata")) # load real data

###############################################################
## set simulation design elements
###############################################################
family = c("lognormal", "loglogistic")
n = c(100, 200, 500)
nS = c(50, 100, 500)
beta_type = c('linear')
b = c(0.5)
seed_start = 1000
N_iter = 100

params = expand.grid(seed_start = seed_start,
                     family = family,
                     n = n,
                     nS = nS,
                     beta_type = beta_type,
                     b = b)

## define number of simulations and parameter scenarios
if(doLocal) {
  scenario = 1
  N_iter = 2
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
  sim_data_test <- simulate_AFT2(family = family, n = n, nS = nS, beta_type = beta_type, b = b, u = 150, seed = SEED.START)
}else{
  sim_data_test <- simulate_afcm(n = n, nS = nS, beta_type = beta_type, seed = SEED.START)
}

for(iter in 1:N_iter){
  print(iter)
  # set seed
  seed.iter = (SEED.START - 1) * N_iter + iter

  # simulate data
  if(family %in% c("lognormal", "loglogistic")){
    sim_data <- simulate_AFT2(family = family, n = n, nS = nS, beta_type = beta_type, b = b, u = 150, seed = seed.iter)
  }else{
    sim_data <- simulate_afcm(n = n, nS = nS, beta_type = beta_type, seed = seed.iter)
  }

  ###############################################################
  ## fit functional AFT and Cox model
  ###############################################################
  # additive functional log-normal AFT model
  tic()
  fit.aaft <- gam(logY ~ ti(X, S, by = L, bs = c("cr", "cr"), k = c(20, 20), mc = c(TRUE, FALSE)), family = cnorm(), data = sim_data$data)
  time_stamp <- toc(quiet = TRUE)
  time.aaft <- time_stamp$toc - time_stamp$tic
  
  # additive functional Cox model
  tic()
  fit.afcm <- gam(Y ~ ti(X, S, by = L, bs = c("cr", "cr"), k = c(20, 20), mc = c(TRUE, FALSE)), weights = delta, family = cox.ph(), data = sim_data$data)
  time_stamp <- toc(quiet = TRUE)
  time.afcm <- time_stamp$toc - time_stamp$tic
  
  # linear functional log-normal AFT model
  tic()
  fit.laft <- gam(logY ~ 1 + s(S, by = X_L, bs = "ps", k = 30), family = cnorm(), data = sim_data$data)
  time_stamp <- toc(quiet = TRUE)
  time.laft <- time_stamp$toc - time_stamp$tic
  
  # linear functional Cox model
  tic()
  fit.lfcm <- gam(Y ~ s(S, by = X_L, bs = "ps", k = 30), weights = delta, family = cox.ph(), data = sim_data$data)
  time_stamp <- toc(quiet = TRUE)
  time.lfcm <- time_stamp$toc - time_stamp$tic
  
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
  eta_aaft <- predict(fit.aaft, sim_data_test$data, type = "response")
  eta_afcm <- rowSums(predict(fit.afcm, sim_data_test$data, type = "terms"))
  eta_laft <- predict(fit.laft, sim_data_test$data, type = "response")
  eta_lfcm <- rowSums(predict(fit.lfcm, sim_data_test$data, type = "terms"))
  
  ## calculate c-index
  AUC_aaft <- cal_c(marker = -eta_aaft, Stime = time_test, status = event_test)
  AUC_afcm <- cal_c(marker = eta_afcm, Stime = time_test, status = event_test)
  AUC_laft <- cal_c(marker = -eta_laft, Stime = time_test, status = event_test)
  AUC_lfcm <- cal_c(marker = eta_lfcm, Stime = time_test, status = event_test)
  
  ## calculate the brier score
  tmax_test <- round(quantile(sim_data_test$data$t, 0.99))
  tgrid_test <- seq(0, tmax_test, length.out = 1000)
  
  S_aaft <- cal_stime(fit = fit.aaft, data = sim_data_test$data, tgrid = tgrid_test, family = 'lognormal')
  S_afcm <- cal_stime(fit = fit.afcm, data = sim_data_test$data, tgrid = tgrid_test, family = 'cox.ph')
  S_laft <- cal_stime(fit = fit.laft, data = sim_data_test$data, tgrid = tgrid_test, family = 'lognormal')
  S_lfcm <- cal_stime(fit = fit.lfcm, data = sim_data_test$data, tgrid = tgrid_test, family = 'cox.ph')
  
  Brier_aaft <- cal_Brier(S_aaft, Stime = time_test, status = event_test, tgrid = tgrid_test)
  Brier_afcm <- cal_Brier(S_afcm, Stime = time_test, status = event_test, tgrid = tgrid_test)
  Brier_laft <- cal_Brier(S_laft, Stime = time_test, status = event_test, tgrid = tgrid_test)
  Brier_lfcm <- cal_Brier(S_lfcm, Stime = time_test, status = event_test, tgrid = tgrid_test)
  
  ###############################################################
  ## pointwise squared errors, pointwise CIs and CMA CIs for estimated beta
  ###############################################################
  # calculate pointwise squared errors
  xgrid <- as.matrix(seq(-40, 10, len = nS), ncol = 1)
  sgrid <- as.matrix(seq(0, 1, len = 500), ncol = 1)
  coef.true <- t(-0.5 * xgrid %*% t(sgrid))
  
  df_pred <- expand.grid(X = xgrid, S = sgrid)
  df_pred$L = 1
  pred_smooth_aaft <- predict(fit.aaft, newdata = df_pred, type = "terms", terms = "ti(X,S):L")
  coef.est.aaft <- matrix(pred_smooth_aaft[,"ti(X,S):L"], ncol = nS, nrow = 500, byrow = TRUE)
  pred_smooth_afcm <- predict(fit.afcm, newdata = df_pred, type = "terms", terms = "ti(X,S):L")
  coef.est.afcm <- matrix(pred_smooth_afcm[,"ti(X,S):L"], ncol = nS, nrow = 500, byrow = TRUE)
  
  # coef.true <- sim_data$coefficients$beta1
  # coef.est.aaft <- predict(fit.aaft, newdata = df_pred, type = "terms", se.fit = TRUE)
  # coef.est.afcm <- predict(fit.afcm, newdata = df_pred, type = "terms", se.fit = TRUE)
  
  if (family %in% c("lognormal", "loglogistic")) {
    se.coef.aaft <- (coef.true - coef.est.aaft)^2
    se.coef.afcm <- (-coef.true - coef.est.afcm)^2
  } else {
    se.coef.aaft <- (-coef.true - coef.est.aaft[[1]])^2
    se.coef.afcm <- (coef.true - coef.est.afcm[[1]])^2
  }
  

  # summary
  df_coef <- data.frame(true_coef = I(coef.true),                           # true coefficient function
                        est_coef_aaft = I(coef.est.aaft),  # estimated coefficient function
                        se_coef_aaft = rowMeans(se.coef.aaft),         # pointwise squared error
                        #lb_coef_aaft = as.numeric(coef.est.aaft[[1]] - qnorm(0.975) * coef.est.aaft[[2]]), # pointwise CI lower bound
                        #ub_coef_aaft = as.numeric(coef.est.aaft[[1]] + qnorm(0.975) * coef.est.aaft[[2]]), # pointwise CI upper bound
                        #cma_lb_coef_aaft = as.numeric(cma.coef.aaft[[1]]), # CMA CI lower bound
                        #cma_ub_coef_aaft = as.numeric(cma.coef.aaft[[2]]), # CMA CI upper bound
                        est_coef_afcm = I(coef.est.afcm),
                        se_coef_afcm = rowMeans(se.coef.afcm))
  
  ###############################################################
  ## pointwise squared errors for survival function
  ###############################################################
  # set the time grid to evaluate survival function
  if (family == "cox.ph"){
    tmax <- 300 # consistent with 'tmax' in simulate_afcm()
    tgrid <- seq(0, tmax, len = 1000) 
  } else {
    tmax <- round(quantile(sim_data$data$t, 0.99))
    tgrid <- seq(0, tmax, len = 1000)
    lp = sim_data$data$lp
    scale = sim_data$b
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
  S_aaft <- cal_stime(fit = fit.aaft, data = sim_data$data, tgrid = tgrid, family = 'lognormal')
  S_afcm <- cal_stime(fit = fit.afcm, data = sim_data$data, tgrid = tgrid, family = 'cox.ph')
  S_laft <- cal_stime(fit = fit.laft, data = sim_data$data, tgrid = tgrid, family = 'lognormal')
  S_lfcm <- cal_stime(fit = fit.lfcm, data = sim_data$data, tgrid = tgrid, family = 'cox.ph')
  # S_faft <- outer(lp_faft, tgrid, function(lp_faft_i, tgrid_j) pnorm((log(tgrid_j) - lp_faft_i) / scale_faft, 
  #                                                                                lower.tail = FALSE))
  
  # calculate pointwise squared error and MISE
  df_surv <- data.frame(time = tgrid,
                        #true_surv = 1 - p.true,
                        #est_surv_aaft = 1 - p.est.aaft,
                        #est_surv_afcm = 1 - p.est.afcm,
                        se_surv_aaft = colMeans((S_true - S_aaft)^2),
                        se_surv_afcm = colMeans((S_true - S_afcm)^2),
                        se_surv_laft = colMeans((S_true - S_laft)^2),
                        se_surv_lfcm = colMeans((S_true - S_lfcm)^2))
  
  df_info <- data.frame(scenario = scenario,
                        iter = iter,
                        seed = seed.iter,
                        family = family,
                        n = n,
                        nS = nS,
                        beta_type = beta_type,
                        b = b,
                        censor_rate = 1 - mean(sim_data$data$delta),
                        AUC_aaft,
                        AUC_afcm,
                        AUC_laft,
                        AUC_lfcm,
                        Brier_aaft,
                        Brier_afcm,
                        Brier_laft,
                        Brier_lfcm,
                        time_aaft = time.aaft,
                        time_afcm = time.afcm,
                        time_laft = time.laft,
                        time_lfcm = time.lfcm)
  
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


