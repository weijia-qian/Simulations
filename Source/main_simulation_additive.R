####################################################################
# Weijia Qian
# November 2024
#
# This file simulates survival data under different data generation mechanisms
# and fits functional AFT and Cox models
####################################################################

suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(refund))
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
source(here("Source", "simulate_Cox2.R"))
source(here("Source", "calc_auc_brier.R"))
load(here("Source", "dat_func.Rdata")) # load real data

###############################################################
## set simulation design elements
###############################################################
family = c("lognormal", "loglogistic", "cox.ph")
n = c(200)
nS = c(50, 100, 500)
beta_type = c('nonlinear1')
b = c(0.5)
seed_start = 1234
N_iter = 100

params = expand.grid(family = family,
                     n = n,
                     nS = nS,
                     beta_type = beta_type,
                     b = b)

## define number of simulations and parameter scenarios
if(doLocal) {
  scenario = 2
  N_iter = 1
}else{
  # defined from batch script params
  scenario <- as.numeric(commandArgs(trailingOnly=TRUE))
}

n = params$n[scenario]
family = params$family[scenario]
nS = params$nS[scenario]
beta_type = params$beta_type[scenario]
b = params$b[scenario]

###############################################################
## run simulations
###############################################################
results = vector("list", length = N_iter)

# simulate a test dataset
if (family %in% c("lognormal", "loglogistic")) {
  sim_data_test <- simulate_AFT2(family = family, n = n, nS = nS, beta_type = beta_type, b = b, u = 2000, seed = seed_start)
} else {
  sim_data_test <- simulate_Cox2(n = n, nS = nS, beta_type = beta_type, tmax = 2000, u = 2000, seed = seed_start)
}

for (iter in 1:N_iter) {
  print(iter)
  # set seed
  seed.iter = (seed_start - 1) * N_iter + iter

  # simulate data
  if (family %in% c("lognormal", "loglogistic")) {
    sim_data <- simulate_AFT2(family = family, n = n, nS = nS, beta_type = beta_type, b = b, u = 2000, seed = seed.iter)
  } else {
    sim_data <- simulate_Cox2(n = n, nS = nS, beta_type = beta_type, tmax = 2000, u = 2000, seed = seed.iter)
  }

  ###############################################################
  ## fit functional AFT and Cox model
  ###############################################################
  # additive functional log-normal AFT model
  tic()
  fit.aaft <- gam(logY ~ ti(X, S, by = L, bs = c("cr", "cr"), k = c(20, 20), mc = c(TRUE, FALSE)), family = cnorm(), data = sim_data$data)
  time_stamp <- toc(quiet = TRUE)
  time_aaft <- time_stamp$toc - time_stamp$tic
  
  # additive functional Cox model
  tic()
  fit.afcm <- gam(Y ~ ti(X, S, by = L, bs = c("cr", "cr"), k = c(20, 20), mc = c(TRUE, FALSE)), weights = delta, family = cox.ph(), data = sim_data$data)
  time_stamp <- toc(quiet = TRUE)
  time_afcm <- time_stamp$toc - time_stamp$tic
  
  # linear functional log-normal AFT model
  tic()
  fit.laft <- gam(logY ~ 1 + s(S, by = X_L, bs = "ps", k = 20), family = cnorm(), data = sim_data$data)
  time_stamp <- toc(quiet = TRUE)
  time_laft <- time_stamp$toc - time_stamp$tic
  
  # linear functional Cox model
  tic()
  fit.lfcm <- gam(Y ~ s(S, by = X_L, bs = "ps", k = 20), weights = delta, family = cox.ph(), data = sim_data$data)
  time_stamp <- toc(quiet = TRUE)
  time_lfcm <- time_stamp$toc - time_stamp$tic
  
  ###############################################################
  ## Out-of-sample Harrell’s C-index and Brier score
  ###############################################################
  # survival time and status
  time_train <- sim_data$data$Y
  event_train <- sim_data$data$delta
  time_test <- sim_data_test$data$Y
  event_test <- sim_data_test$data$delta

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
  
  # set the time grid to evaluate survival function
  tmax <- 2000 # consistent with 'tmax' in simulate_Cox2()
  tgrid <- seq(0, tmax, len = 1000) 
  
  # estimated survival function
  S_aaft_test <- cal_stime(fit = fit.aaft, data = sim_data_test$data, tgrid = tgrid, family = 'lognormal')
  S_afcm_test <- cal_stime(fit = fit.afcm, data = sim_data_test$data, tgrid = tgrid, family = 'cox.ph')
  S_laft_test <- cal_stime(fit = fit.laft, data = sim_data_test$data, tgrid = tgrid, family = 'lognormal')
  S_lfcm_test <- cal_stime(fit = fit.lfcm, data = sim_data_test$data, tgrid = tgrid, family = 'cox.ph')
  
  ## calculate brier score
  Brier_aaft <- cal_Brier(S_aaft_test, Stime = time_test, status = event_test, tgrid = tgrid)
  Brier_afcm <- cal_Brier(S_afcm_test, Stime = time_test, status = event_test, tgrid = tgrid)
  Brier_laft <- cal_Brier(S_laft_test, Stime = time_test, status = event_test, tgrid = tgrid)
  Brier_lfcm <- cal_Brier(S_lfcm_test, Stime = time_test, status = event_test, tgrid = tgrid)
  
  ###############################################################
  ## pointwise squared errors of coefficient surface
  ###############################################################
  # set the grid to evaluate coefficient surface
  kX = kS = 500
  xgrid <- seq(-40, 10, len = kX)
  sgrid <- seq(0, 1, len = kS)
  df_coef <- expand.grid(X = xgrid, S = sgrid)
  df_coef$L = 1
  
  # true coefficient surface
  beta <- sim_data$beta
  df_coef$true_coef <- beta(df_coef$X, df_coef$S)
  
  # estimated coefficient surface
  df_coef$est_coef_aaft <- as.numeric(predict(fit.aaft, newdata = df_coef, type = "terms", terms = "ti(X,S):L"))
  df_coef$est_coef_afcm <- as.numeric(predict(fit.afcm, newdata = df_coef, type = "terms", terms = "ti(X,S):L"))
  
  # calculate pointwise squared errors
  df_coef$se_coef_aaft <- (df_coef$true_coef - df_coef$est_coef_aaft)^2
  df_coef$se_coef_afcm <- (df_coef$true_coef - df_coef$est_coef_afcm)^2

  ###############################################################
  ## pointwise squared errors ofr survival function
  ###############################################################
  if (family %in% c("lognormal", "loglogistic")) {
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

  # calculate pointwise squared errors
  df_surv <- data.frame(time = tgrid,
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
                        time_aaft,
                        time_afcm,
                        time_laft,
                        time_lfcm)
  
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


