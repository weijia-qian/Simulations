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
source(here("Source", "calc_auc_brier.R"))
source(here("Source", "utils_summary.R"))
source(here("Source", "fit_penalized_aft.R"))
source(here("Source", "FTTM.R"))
load(here("Source", "dat_func.Rdata")) # load real data

###############################################################
## set simulation design elements
###############################################################
family = c("lognormal", "loglogistic", "cox.ph")
# n = c(100, 200, 500)
n = c(100, 200, 500)
nS = c(100)
# nS = c(50, 100, 500)
beta_type = c('simple')
b = c(0.5)
seed_start = 1000
N_iter = 100

params = expand.grid(family = family,
                     n = n,
                     nS = nS,
                     beta_type = beta_type,
                     b = b)

params <- params %>%
  filter(family == "lognormal" | n == 100)

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

for(iter in 1:N_iter){
  print(iter)
  # set seed
  seed.iter = (seed_start - 1) * N_iter + iter

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
  lambda_grid <- exp(seq(log(1000), log(10000), length.out = 500))
  model <- "lognormal"
  tic()
  best_lambda <- optimize_lambda(lambda_grid = lambda_grid, data = sim_data$data, y = "Y", delta = "delta",
                                 x = "X", x_as_regex = FALSE, family = model)$lambda
  fit.faft <- optimize_AFT(data = sim_data$data, y = "Y", delta = "delta", x = "X", x_as_regex = FALSE, 
                           family = model, lambda = best_lambda, se = TRUE, bootstrap = FALSE)
  time_stamp <- toc(quiet = TRUE)
  time_faft <- time_stamp$toc - time_stamp$tic
  
  # loglogistic lfAFT
  model <- "loglogistic"
  tic()
  best_lambda <- optimize_lambda(lambda_grid = lambda_grid, data = sim_data$data, y = "Y", delta = "delta",
                                 x = "X", x_as_regex = FALSE, family = model)$lambda
  fit.faft2 <- optimize_AFT(data = sim_data$data, y = "Y", delta = "delta", x = "X", x_as_regex = FALSE, 
                           family = model, lambda = best_lambda, se = TRUE, bootstrap = FALSE)
  time_stamp <- toc(quiet = TRUE)
  time_faft2 <- time_stamp$toc - time_stamp$tic
  
  # Functional Time Transformation Model (FTTM)
  tic()
  S <- seq(0, 1, len = nS)
  T_obs = sim_data$data$Y # observed survival time
  delta = sim_data$data$delta # censoring indicator
  n = length(T_obs)
  Xmat = matrix(1, nrow = n, ncol = 1) # intercept
  Xs = as.matrix(sim_data$data$X) # functional covariate matrix
  tau <- ceiling(max(T_obs))
  p1 <- ncol(Xmat)
  fit.fttm <- fit_FTTM(data = sim_data$data)
  time_stamp <- toc(quiet = TRUE)
  time_fttm <- time_stamp$toc - time_stamp$tic


  df_info <- data.frame(scenario = scenario,
                        iter = iter,
                        seed = seed.iter,
                        family = family,
                        n = n,
                        nS = nS,
                        beta_type = beta_type,
                        b = b,
                        censor_rate = 1 - mean(sim_data$data$delta),
                        time_norm,
                        time_cox,
                        time_faft,
                        time_faft2,
                        time_fttm
                        )
  
  list(
    info = df_info
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
Date = paste0(gsub("-", "", Sys.Date()), "_fttm")
dir.create(file.path(here::here("Output"), Date), showWarnings = FALSE)

filename = paste0(here::here("Output", Date), "/", scenario, ".RDA")
save(results,
     file = filename)

###############################################################
## end sim
###############################################################


