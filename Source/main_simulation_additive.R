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
suppressPackageStartupMessages(library(survAUC))
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
source(here("Source", "fit_penalized_aft.R"))
load(here("Source", "dat_func.Rdata")) # load real data

###############################################################
## set simulation design elements
###############################################################
family = c("lognormal", "cox.ph")
n = c(100, 200, 500)
nS = c(100)
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
  scenario = 1
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
  
  res <- tryCatch({
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
  
  # linear functional log-normal AFT model (mgcv)
  tic()
  fit.laft <- gam(logY ~ 1 + s(S, by = X_L, bs = "ps", k = 20), family = cnorm(), data = sim_data$data)
  time_stamp <- toc(quiet = TRUE)
  time_laft <- time_stamp$toc - time_stamp$tic
  
  # lognormal lfAFT
  lambda_grid <- exp(seq(log(1000), log(10000), length.out = 500))
  model <- "lognormal"
  tic()
  best_lambda <- optimize_lambda(lambda_grid = lambda_grid, data = sim_data$data, y = "Y", delta = "delta",
                                 x = "X", x_as_regex = FALSE, family = model)$lambda
  fit.laft2 <- optimize_AFT(data = sim_data$data, y = "Y", delta = "delta", x = "X", x_as_regex = FALSE, 
                           family = model, lambda = best_lambda, se = TRUE, bootstrap = FALSE)
  time_stamp <- toc(quiet = TRUE)
  time_laft2 <- time_stamp$toc - time_stamp$tic
  
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
  eta_laft2 <- predict_AFT(fit.laft2, sim_data_test$data)
  eta_lfcm <- rowSums(predict(fit.lfcm, sim_data_test$data, type = "terms"))
  
  ## Surv objects for UnoC (Uno's IPCW C-index, consistent with real_data_analysis.Rmd)
  Surv_train <- Surv(time_train, event_train)
  Surv_test  <- Surv(time_test,  event_test)

  ## calculate C-index (Uno's IPCW; sign-flip so higher lpnew = higher risk)
  # Note: eta_laft2 is overwritten below as as.numeric(); compute C-index first.
  Cindex_aaft  <- UnoC(Surv_train, Surv_test, lpnew = -eta_aaft)
  Cindex_afcm  <- UnoC(Surv_train, Surv_test, lpnew =  eta_afcm)
  Cindex_laft  <- UnoC(Surv_train, Surv_test, lpnew = -eta_laft)
  Cindex_laft2 <- UnoC(Surv_train, Surv_test, lpnew = -eta_laft2)
  Cindex_lfcm  <- UnoC(Surv_train, Surv_test, lpnew =  eta_lfcm)

  # set the time grid to evaluate survival function
  tmax <- 2000 # consistent with 'tmax' in simulate_Cox2()
  nS_pred <- 500
  tgrid <- seq(0.1, tmax, len = nS_pred)  # start at 0.1 to avoid log(0) in AFT survival

  # estimated survival function
  S_aaft_test <- cal_stime(fit = fit.aaft, data = sim_data_test$data, tgrid = tgrid, family = 'lognormal')
  df_pred <- sim_data_test$data[rep(1:n, each = nS_pred), ]
  df_pred$Y <- rep(tgrid, n)
  S_afcm_test <- matrix(predict(fit.afcm, newdata = df_pred, type = "response"), nrow = n, ncol = nS_pred, byrow = TRUE)
  S_laft_test <- cal_stime(fit = fit.laft, data = sim_data_test$data, tgrid = tgrid, family = 'lognormal')
  S_lfcm_test <- matrix(predict(fit.lfcm, newdata = df_pred, type = "response"), nrow = n, ncol = nS_pred, byrow = TRUE)

  eta_laft2 <- as.numeric(eta_laft2)
  scale_laft2 <- fit.laft2$b_hat
  S_laft2_test <- outer(eta_laft2, tgrid,
                        function(eta_laft2_i, tgrid_j) pnorm((log(tgrid_j) - eta_laft2_i) / scale_laft2, lower.tail = FALSE))

  ## calculate IPCW integrated Brier score (consistent with real_data_analysis.Rmd)
  Brier_aaft  <- cal_IPCW_IBS(S_aaft_test,  time_test, event_test, time_train, event_train, tgrid)
  Brier_afcm  <- cal_IPCW_IBS(S_afcm_test,  time_test, event_test, time_train, event_train, tgrid)
  Brier_laft  <- cal_IPCW_IBS(S_laft_test,  time_test, event_test, time_train, event_train, tgrid)
  Brier_laft2 <- cal_IPCW_IBS(S_laft2_test, time_test, event_test, time_train, event_train, tgrid)
  Brier_lfcm  <- cal_IPCW_IBS(S_lfcm_test,  time_test, event_test, time_train, event_train, tgrid)
  
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
  df_pred <- sim_data$data[rep(1:n, each = nS_pred), ]
  df_pred$Y <- rep(tgrid, n)
  S_afcm <- matrix(predict(fit.afcm, newdata = df_pred, type = "response"), nrow = n, ncol = nS_pred, byrow = TRUE)
  S_laft <- cal_stime(fit = fit.laft, data = sim_data$data, tgrid = tgrid, family = 'lognormal')
  S_lfcm <- matrix(predict(fit.lfcm, newdata = df_pred, type = "response"), nrow = n, ncol = nS_pred, byrow = TRUE)
  
  lp_laft2 <- as.numeric(fit.laft2$lp)
  S_laft2 <- outer(lp_laft2, tgrid, 
                  function(lp_laft2_i, tgrid_j) pnorm((log(tgrid_j) - lp_laft2_i) / scale_laft2, lower.tail = FALSE))

  # calculate pointwise squared errors
  df_surv <- data.frame(time = tgrid,
                        se_surv_aaft = colMeans((S_true - S_aaft)^2),
                        se_surv_afcm = colMeans((S_true - S_afcm)^2),
                        se_surv_laft = colMeans((S_true - S_laft)^2),
                        se_surv_laft2 = colMeans((S_true - S_laft2)^2),
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
                        # C-index (Uno's IPCW, out-of-sample)
                        Cindex_aaft,
                        Cindex_afcm,
                        Cindex_laft,
                        Cindex_laft2,
                        Cindex_lfcm,
                        # IPCW integrated Brier score (out-of-sample)
                        Brier_aaft,
                        Brier_afcm,
                        Brier_laft,
                        Brier_laft2,
                        Brier_lfcm,
                        time_aaft,
                        time_afcm,
                        time_laft,
                        time_laft2,
                        time_lfcm)
  
  list(
    info = df_info,
    coef = df_coef,
    surv = df_surv
  )
  
  }, error = function(e) {
    warning(sprintf(
      "Iteration %d skipped due to error:\n  %s",
      iter, e$message
    ))
    NULL
  })
  
  ## only save non-NULL results
  if (!is.null(res)) {
    results[[iter]] <- res
  }
  
} # end for loop

# drop NULL entries
results <- Filter(Negate(is.null), results)

# record date for analysis; create directory for results
Date = paste0(gsub("-", "", Sys.Date()), "_additive")
dir.create(file.path(here::here("Output"), Date), showWarnings = FALSE)

filename = paste0(here::here("Output", Date), "/", scenario, ".RDA")
save(results,
     file = filename)

###############################################################
## end sim
###############################################################


