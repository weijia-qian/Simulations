####################################################################
# Weijia Qian
# November 2024
#
# This file simulates survival data under different data generation mechanisms
# and fits functional AFT and Cox models
####################################################################

suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(fda))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(refund))
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(survival))
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
source(here::here("Source", "simulate_AFT.R"))
source(here::here("Source", "simulate_Cox.R"))
source(here::here("Source", "est_sieve.R"))
source(here::here("Source", "calc_auc_brier.R"))
source(here::here("Source", "utils_summary.R"))
load(here::here("Source", "dat_func.Rdata")) # load real data

###############################################################
## set simulation design elements
###############################################################
n = c(100, 200, 500)
family = c("lognormal", "loglogistic", "cox.ph")
nS = c(50, 100, 500)
k = c(6, 8)
seed_start = 1000
N_iter = 500

params = expand.grid(seed_start = seed_start,
                     n = n,
                     family = family,
                     nS = nS,
                     k = k)
params = params[!(params$family == "cox.ph" & params$k == 8),]

## define number of simulations and parameter scenarios
if(doLocal) {
  scenario = 7
  N_iter = 2
}else{
  # defined from batch script params
  scenario <- as.numeric(commandArgs(trailingOnly=TRUE))
}

n = params$n[scenario]
family = params$family[scenario]
nS = params$nS[scenario]
k = params$k[scenario]
SEED.START = params$seed_start[scenario]

###############################################################
## run simulations
###############################################################
results = vector("list", length = N_iter)
for(iter in 1:N_iter){
  # set seed
  seed.iter = (SEED.START - 1)*N_iter + iter

  # simulate data
  if(family %in% c("lognormal", "loglogistic")){
    sim_data <- simulate_AFT(family = family, n = n, nS = nS, k = k, seed = seed.iter)
  }else{
    sim_data <- simulate_Cox(n = n, nS = nS, seed = seed.iter)
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
  
  ###############################################################
  ## Harrell’s C-index and Brier score
  ###############################################################
  #set.seed(seed.iter) 
  #n_folds <- 10 # number of folds
  ### variables that store the c-index and Brier's score in each calculation
  #AUC_norm <- AUC_cox <- rep(0, n_folds)
  #Brier_norm <- Brier_cox <- rep(0, n_folds)
  ### use "createFolds" function from "caret" package to divide real data into n_folds folds
  #fold <- createFolds(sim_data$data$Y, k = n_folds, list = FALSE)
  
  ## survival time and status
  time_test <- sim_data$data$Y
  event_test <- sim_data$data$delta

  ## get unique ordered survival times
  ut_test <- unique(time_test[event_test == 1])
  ut_test <- ut_test[order(ut_test)]
  
  ## derive the KM estimate of the censoring time
  ut_test_censor <- unique(time_test[event_test == 0])
  ut_test_censor <- ut_test_censor[order(ut_test_censor)]
  ## get KM estimates of censoring time for test data
  KM_test_censor <- unique(survfit(Surv(time_test,1-event_test)~1)$surv)

  ## obtain linear predictors
  eta_norm <- -fit.norm$linear.predictors
  #eta_cox <- rowSums(predict(fit.cox, test_data, type="terms"))
  eta_cox <- fit.cox$linear.predictors
  
  ## calculate c-index
  AUC_norm <- cal_c(marker = eta_norm, Stime = time_test, status = event_test)
  AUC_cox <- cal_c(marker = eta_cox, Stime = time_test, status = event_test)
  
  ## calculate the brier score
  tmax <- round(quantile(sim_data$data$t, 0.95))
  tvec <- seq(0, tmax, length.out = 1000)
  
  S_norm <- cal_stime(fit = fit.norm, tgrid = tvec, family = 'lognormal')
  S_cox <- cal_stime(fit = fit.cox, tgrid = tvec, family = 'cox.ph')
  
  Brier_norm <- cal_Brier(S_norm, Stime = time_test, status = event_test, tgrid = tvec)
  Brier_cox <- cal_Brier(S_cox, Stime = time_test, status = event_test, tgrid = tvec)
  
  ###############################################################
  ## pointwise squared errors, pointwise CIs and CMA CIs for estimated beta
  ###############################################################
  # calculate pointwise squared errors
  svec <- seq(0, 1, len = nS)
  df_pred <- data.frame(S = svec, X_L = 1)
  coef.true <- sim_data$coefficients$beta1
  coef.est.norm <- predict(fit.norm, newdata = df_pred, type = "terms", se.fit = TRUE)
  coef.est.cox <- predict(fit.cox, newdata = df_pred, type = "terms", se.fit = TRUE)
  if (family %in% c("lognormal", "loglogistic")) {
    se.coef.norm <- (coef.true - coef.est.norm[[1]])^2
    se.coef.cox <- (-coef.true - coef.est.cox[[1]])^2
  } else {
    se.coef.norm <- (-coef.true - coef.est.norm[[1]])^2
    se.coef.cox <- (coef.true - coef.est.cox[[1]])^2
  }
  #mise.coef.norm <- mean(se.coef.norm)
  #mise.coef.cox <- mean(se.coef.cox)
  
  # calculate CMA CIs
  cma.coef.norm <- get_CMA(fit.norm)
  cma.coef.cox <- get_CMA(fit.cox)
  
  # fit AFT model using sieve algorithm
  sieve.results <- est_sieve(data = sim_data$data)

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
                        lb_coef_cox = as.numeric(coef.est.cox[[1]] - qnorm(0.975) * coef.est.cox[[2]]), 
                        ub_coef_cox = as.numeric(coef.est.cox[[1]] + qnorm(0.975) * coef.est.cox[[2]]),
                        cma_lb_coef_cox = as.numeric(cma.coef.cox[[1]]),
                        cma_ub_coef_cox = as.numeric(cma.coef.cox[[2]]),
                        est_coef_sieve = as.numeric(sieve.results[[1]]),
                        se_coef_sieve = as.numeric(sieve.results[[2]])) %>%
    mutate(cover_coef_norm = (true_coef > lb_coef_norm) & (true_coef < ub_coef_norm),
           cover_coef_cox = (true_coef > lb_coef_cox) & (true_coef < ub_coef_cox),
           cover_cma_coef_norm = (true_coef > cma_lb_coef_norm) & (true_coef < cma_ub_coef_norm),
           cover_cma_coef_cox = (true_coef > cma_lb_coef_cox) & (true_coef < cma_ub_coef_cox))
  
  ###############################################################
  ## pointwise squared errors for survival function
  ###############################################################
  # true survival function
  if (family == "lognormal") {
    lp = sim_data$data$lp
    scale = sim_data$coefficients$b[1]
    S_true <- outer(lp, tvec, function(lp_i, tvec_j) pnorm((log(tvec_j) - lp_i) / scale, lower.tail = FALSE))
  } else if (family == "loglogistic") {
    lp = sim_data$data$lp
    scale = sim_data$coefficients$b[1]
    S_true <- outer(lp, tvec, function(lp_i, tvec_j) 1 - 1 / (1 + (exp(lp_i) / tvec_j)^(1 / scale)))
  } else {
    S_true = sim_data$data$Si
  }
  
  # calculate pointwise squared error and MISE
  df_surv <- data.frame(time = tvec,
                        #true_surv = 1 - p.true,
                        #est_surv_norm = 1 - p.est.norm,
                        #est_surv_cox = 1 - p.est.cox,
                        se_surv_norm = colMeans((S_true - S_norm)^2),
                        se_surv_cox = colMeans((S_true - S_cox)^2))
  #mise.surv.norm <- mean(se.surv.norm)
  #mise.surv.cox <- mean(se.surv.cox)
  
  df_info <- data.frame(scenario = scenario,
                        iter = iter,
                        seed = seed.iter,
                        n = n,
                        family = family,
                        nS = nS,
                        k = k,
                        AUC_norm,
                        AUC_cox,
                        Brier_norm,
                        Brier_cox,
                        time_norm = time.norm,
                        time_cox = time.cox,
                        time_sieve = sieve.results[[3]])
  
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


