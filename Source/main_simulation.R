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
source(here("Source", "simulate_AFT.R"))
source(here("Source", "simulate_Cox.R"))
source(here("Source", "calc_auc_brier.R"))
source(here("Source", "utils_summary.R"))
source(here("Source", "fit_penalized_aft.R"))
source(here("Source", "FTTM.R"))
load(here("Source", "dat_func.Rdata")) # load real data

###############################################################
## set simulation design elements
###############################################################
family = c("lognormal", "loglogistic", "cox.ph")
n = c(100, 200, 500)
nS = c(50, 100, 500)
beta_type = c('simple')
b = c(0.5)
N_iter = 500

params = expand.grid(family = family,
                     n = n,
                     nS = nS,
                     beta_type = beta_type,
                     b = b)

## define number of simulations and parameter scenarios
if(doLocal) {
  scenario = 3
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

###############################################################
## run simulations
###############################################################
results = vector("list", length = N_iter)

# Pre-draw per-iteration seeds (reproducible given the R session seed at this point)
seeds <- sample.int(1e8, N_iter)

# simulate a fixed test dataset — seed derived from scenario for reproducibility
# across runs without depending on seed_start
if(family %in% c("lognormal", "loglogistic")){
  sim_data_test <- simulate_AFT(family = family, n = n, nS = nS, beta_type = beta_type, b = b,
                                seed = scenario * 1000L + 99L)
}else{
  sim_data_test <- simulate_Cox(n = n, nS = nS, beta_type = beta_type,
                                seed = scenario * 1000L + 99L)
}

for(iter in 1:N_iter){
  print(iter)
  # set seed for this iteration
  seed.iter <- seeds[[iter]]

  # simulate data
  if(family %in% c("lognormal", "loglogistic")){
    sim_data <- simulate_AFT(family = family, n = n, nS = nS, beta_type = beta_type, b = b, seed = seed.iter)
  }else{
    sim_data <- simulate_Cox(n = n, nS = nS, beta_type = beta_type, seed = seed.iter)
  }

  res <- tryCatch({
  ###############################################################
  ## fit functional AFT and Cox model
  ###############################################################
  # linear functional lognormal AFT model
  tic()
  fit.norm <- gam(logY ~ 1 + s(S, by = X_L, bs = "ps", k = 20), family = cnorm(), data = sim_data$data)
  time_stamp <- toc(quiet = TRUE)
  time_norm <- time_stamp$toc - time_stamp$tic
  
  # linear functional Cox model
  tic()
  fit.cox <- gam(Y ~ s(S, by = X_L, bs = "ps", k = 20), weights = delta, family = cox.ph(), data = sim_data$data)
  time_stamp <- toc(quiet = TRUE)
  time_cox <- time_stamp$toc - time_stamp$tic
  
  # lognormal lfAFT
  tic()
  lambda_grid <- exp(seq(log(1000), log(10000), length.out = 500))
  model <- "lognormal"
  best_lambda <- optimize_lambda(lambda_grid = lambda_grid, data = sim_data$data, y = "Y", delta = "delta",
                                 x = "X", x_as_regex = FALSE, family = model)$lambda
  fit.faft <- optimize_AFT(data = sim_data$data, y = "Y", delta = "delta", x = "X", x_as_regex = FALSE, 
                           family = model, lambda = best_lambda, se = TRUE, bootstrap = FALSE)
  time_stamp <- toc(quiet = TRUE)
  time_faft <- time_stamp$toc - time_stamp$tic
  
  # loglogistic lfAFT
  tic()
  model <- "loglogistic"
  best_lambda <- optimize_lambda(lambda_grid = lambda_grid, data = sim_data$data, y = "Y", delta = "delta",
                                 x = "X", x_as_regex = FALSE, family = model)$lambda
  fit.faft2 <- optimize_AFT(data = sim_data$data, y = "Y", delta = "delta", x = "X", x_as_regex = FALSE, 
                           family = model, lambda = best_lambda, se = TRUE, bootstrap = FALSE)
  time_stamp <- toc(quiet = TRUE)
  time_faft2 <- time_stamp$toc - time_stamp$tic
  
  # Functional Time Transformation Model (FTTM)
  # tic()
  # S <- seq(0, 1, len = nS)
  # T_obs = sim_data$data$Y # observed survival time
  # delta = sim_data$data$delta # censoring indicator
  # n = length(T_obs)
  # Xmat = matrix(1, nrow = n, ncol = 1) # intercept
  # Xs = as.matrix(sim_data$data$X) # functional covariate matrix
  # tau <- ceiling(max(T_obs))
  # p1 <- ncol(Xmat)
  # fit.fttm <- fit_FTTM(data = sim_data$data)
  # time_stamp <- toc(quiet = TRUE)
  # time_fttm <- time_stamp$toc - time_stamp$tic

  ###############################################################
  ## Out-of-sample Harrell’s C-index and Brier score
  ###############################################################
  # survival time and status
  time_train <- sim_data$data$Y
  event_train <- sim_data$data$delta
  time_test <- sim_data_test$data$Y
  event_test <- sim_data_test$data$delta

  ## obtain linear predictors
  eta_norm <- predict(fit.norm, sim_data_test$data, type = "response")
  eta_cox <- rowSums(predict(fit.cox, sim_data_test$data, type = "terms"))
  eta_faft <- predict_AFT(fit.faft, sim_data_test$data)
  eta_faft2 <- predict_AFT(fit.faft2, sim_data_test$data)
  # eta_fttm <- predict_FTTM(fit.fttm, sim_data_test$data)
  
  ## Surv objects for UnoC (Uno's IPCW C-index, consistent with real_data_analysis.Rmd)
  Surv_train <- Surv(time_train, event_train)
  Surv_test  <- Surv(time_test,  event_test)

  ## calculate Uno's IPCW C-index (sign-flip so higher lpnew = higher risk)
  Cindex_norm  <- UnoC(Surv_train, Surv_test, lpnew = -eta_norm)
  Cindex_cox   <- UnoC(Surv_train, Surv_test, lpnew =  eta_cox)
  Cindex_faft  <- UnoC(Surv_train, Surv_test, lpnew = -eta_faft)
  Cindex_faft2 <- UnoC(Surv_train, Surv_test, lpnew = -eta_faft2)
  # Cindex_fttm <- UnoC(Surv_train, Surv_test, lpnew = -eta_fttm)

  ## calculate Harrell's C-index (no IPCW; sign-flip so higher marker = higher risk)
  Harrell_norm  <- cal_c(marker = -eta_norm,  Stime = time_test, status = event_test)
  Harrell_cox   <- cal_c(marker =  eta_cox,   Stime = time_test, status = event_test)
  Harrell_faft  <- cal_c(marker = -eta_faft,  Stime = time_test, status = event_test)
  Harrell_faft2 <- cal_c(marker = -eta_faft2, Stime = time_test, status = event_test)
  # Harrell_fttm <- cal_c(marker = -eta_fttm, Stime = time_test, status = event_test)

  ## calculate the IPCW Brier score (consistent with real_data_analysis.Rmd)
  tmax <- 120
  nS_pred <- 500
  tgrid <- seq(0.1, tmax, len = nS_pred)  # start at 0.1 to avoid log(0) in AFT survival
  
  S_norm <- cal_stime(fit = fit.norm, data = sim_data_test$data, tgrid = tgrid, family = 'lognormal')
  # S_cox <- cal_stime(fit = fit.cox, data = sim_data_test$data, tgrid = tgrid, family = 'cox.ph')
  df_pred <- sim_data_test$data[rep(1:n, each = nS_pred), ]
  df_pred$Y <- rep(tgrid, n)
  S_cox <- matrix(predict(fit.cox, newdata = df_pred, type = "response"), nrow = n, ncol = nS_pred, byrow = TRUE)
  
  eta_faft <- as.numeric(eta_faft)
  scale_faft <- fit.faft$b_hat
  S_faft <- outer(eta_faft, tgrid, 
                  function(eta_faft_i, tgrid_j) pnorm((log(tgrid_j) - eta_faft_i) / scale_faft, lower.tail = FALSE))
  
  eta_faft2 <- as.numeric(eta_faft2)
  scale_faft2 <- fit.faft2$b_hat
  S_faft2 <- outer(eta_faft2, tgrid, 
                   function(eta_faft2_i, tgrid_j) 1 / (1 + exp((log(tgrid_j) - eta_faft2_i) / scale_faft2)))
  
  Brier_norm  <- cal_IPCW_IBS(S_norm,  time_test, event_test, time_train, event_train, tgrid)
  Brier_cox   <- cal_IPCW_IBS(S_cox,   time_test, event_test, time_train, event_train, tgrid)
  Brier_faft  <- cal_IPCW_IBS(S_faft,  time_test, event_test, time_train, event_train, tgrid)
  Brier_faft2 <- cal_IPCW_IBS(S_faft2, time_test, event_test, time_train, event_train, tgrid)
  
  ###############################################################
  ## pointwise squared errors, pointwise CIs and CMA CIs for estimated beta
  ###############################################################
  # calculate pointwise squared errors
  sgrid <- seq(0, 1, len = nS)
  df_pred <- data.frame(S = sgrid, X_L = 1)
  
  coef.true <- sim_data$coefficients$beta1
  coef.est.norm <- predict(fit.norm, newdata = df_pred, type = "terms", se.fit = TRUE)
  coef.est.cox <- predict(fit.cox, newdata = df_pred, type = "terms", se.fit = TRUE)
  
  se.coef.norm <- (coef.true - coef.est.norm[[1]])^2
  se.coef.cox <- (coef.true - coef.est.cox[[1]])^2
  se.coef.faft <- (coef.true - fit.faft$beta1_hat)^2
  se.coef.faft2 <- (coef.true - fit.faft2$beta1_hat)^2
  # se.coef.fttm <- (coef.true - fit.fttm$beta1_hat)^2
  
  # calculate CMA CIs
  # cma.coef.norm <- get_CMA(fit.norm)
  # cma.coef.cox <- get_CMA(fit.cox)
  
  # fit AFT model using sieve algorithm
  #sieve.results <- est_sieve(data = sim_data$data)

  # summary
  df_coef <- data.frame(true_coef = coef.true,                           # true coefficient function
                        est_coef_norm = as.numeric(coef.est.norm[[1]]),  # estimated coefficient function
                        se_coef_norm = as.numeric(se.coef.norm),         # pointwise squared error
                        lb_coef_norm = as.numeric(coef.est.norm[[1]] - qnorm(0.975) * coef.est.norm[[2]]), # pointwise CI lower bound
                        ub_coef_norm = as.numeric(coef.est.norm[[1]] + qnorm(0.975) * coef.est.norm[[2]]), # pointwise CI upper bound
                        # cma_lb_coef_norm = as.numeric(cma.coef.norm[[1]]), # CMA CI lower bound
                        # cma_ub_coef_norm = as.numeric(cma.coef.norm[[2]]), # CMA CI upper bound
                        est_coef_cox = as.numeric(coef.est.cox[[1]]), 
                        se_coef_cox = as.numeric(se.coef.cox), 
                        #est_coef_sieve = as.numeric(sieve.results[[1]]),
                        #se_coef_sieve = as.numeric(sieve.results[[2]]),
                        lb_coef_cox = as.numeric(coef.est.cox[[1]] - qnorm(0.975) * coef.est.cox[[2]]), 
                        ub_coef_cox = as.numeric(coef.est.cox[[1]] + qnorm(0.975) * coef.est.cox[[2]]),
                        # cma_lb_coef_cox = as.numeric(cma.coef.cox[[1]]),
                        # cma_ub_coef_cox = as.numeric(cma.coef.cox[[2]]),
                        est_coef_faft = fit.faft$beta1_hat, 
                        se_coef_faft = se.coef.faft, 
                        lb_coef_faft = fit.faft$beta1_ci_lower, 
                        ub_coef_faft = fit.faft$beta1_ci_upper,
                        est_coef_faft2 = fit.faft2$beta1_hat, 
                        se_coef_faft2 = se.coef.faft2, 
                        lb_coef_faft2 = fit.faft2$beta1_ci_lower, 
                        ub_coef_faft2 = fit.faft2$beta1_ci_upper
                        ) %>%
    mutate(cover_coef_norm = (true_coef > lb_coef_norm) & (true_coef < ub_coef_norm),
           cover_coef_cox = (true_coef > lb_coef_cox) & (true_coef < ub_coef_cox),
           cover_coef_faft = (true_coef > lb_coef_faft) & (true_coef < ub_coef_faft),
           cover_coef_faft2 = (true_coef > lb_coef_faft2) & (true_coef < ub_coef_faft2)
           # cover_cma_coef_norm = (true_coef > cma_lb_coef_norm) & (true_coef < cma_ub_coef_norm),
           # cover_cma_coef_cox = (true_coef > cma_lb_coef_cox) & (true_coef < cma_ub_coef_cox)
           )
  
  ###############################################################
  ## pointwise squared errors for survival function
  ###############################################################
  # set the time grid to evaluate survival function
  if (family %in% c("lognormal", "loglogistic")) {
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
  # S_cox <- cal_stime(fit = fit.cox, data = sim_data$data, tgrid = tgrid, family = 'cox.ph')
  df_pred <- sim_data$data[rep(1:n, each = nS_pred), ]
  df_pred$Y <- rep(tgrid, n)
  S_cox <- matrix(predict(fit.cox, newdata = df_pred, type = "response"), nrow = n, ncol = nS_pred, byrow = TRUE)
  lp_faft <- as.numeric(fit.faft$lp)
  S_faft <- outer(lp_faft, tgrid, 
                  function(lp_faft_i, tgrid_j) pnorm((log(tgrid_j) - lp_faft_i) / scale_faft, lower.tail = FALSE))
  lp_faft2 <- as.numeric(fit.faft2$lp)
  S_faft2 <- outer(lp_faft2, tgrid, 
                   function(lp_faft2_i, tgrid_j) 1 - 1 / (1 + (exp(lp_faft2_i) / tgrid_j)^(1 / scale_faft2)))
  
  # calculate pointwise squared error and MISE
  df_surv <- data.frame(time = tgrid,
                        #true_surv = 1 - p.true,
                        #est_surv_norm = 1 - p.est.norm,
                        #est_surv_cox = 1 - p.est.cox,
                        se_surv_norm = colMeans((S_true - S_norm)^2),
                        se_surv_cox = colMeans((S_true - S_cox)^2),
                        se_surv_faft = colMeans((S_true - S_faft)^2),
                        se_surv_faft2 = colMeans((S_true - S_faft2)^2)
                        )
  
  df_info <- data.frame(scenario = scenario,
                        iter = iter,
                        seed = seed.iter,
                        family = family,
                        n = n,
                        nS = nS,
                        beta_type = beta_type,
                        b = b,
                        censor_rate = 1 - mean(sim_data$data$delta),
                        # Uno's IPCW C-index (out-of-sample)
                        Cindex_norm,
                        Cindex_cox,
                        Cindex_faft,
                        Cindex_faft2,
                        # Harrell's C-index (out-of-sample, no IPCW correction)
                        Harrell_norm,
                        Harrell_cox,
                        Harrell_faft,
                        Harrell_faft2,
                        # IPCW integrated Brier score (out-of-sample)
                        Brier_norm,
                        Brier_cox,
                        Brier_faft,
                        Brier_faft2,
                        time_norm,
                        time_cox,
                        time_faft,
                        time_faft2
                        )
  
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
Date = gsub("-", "", Sys.Date())
dir.create(file.path(here::here("Output"), Date), showWarnings = FALSE)

filename = paste0(here::here("Output", Date), "/", scenario, ".RDA")
save(results,
     file = filename)

###############################################################
## end sim
###############################################################


