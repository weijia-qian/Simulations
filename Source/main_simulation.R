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
suppressPackageStartupMessages(library(refund))
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(splines))
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
source(here::here("source", "simulate_AFT.R"))
source(here::here("source", "simulate_Cox.R"))
source(here::here("source", "utils_summary.R"))

###############################################################
## set simulation design elements
###############################################################

n = c(100, 200, 500)
family = c("lognormal", "loglogistic", "cox.ph")
nS = c(50, 100, 500)
seed_start = 1000
N_iter = 50

params = expand.grid(seed_start = seed_start,
                     n = n,
                     family = family,
                     nS = nS)

## record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())
dir.create(file.path(here::here("Output"), Date), showWarnings = FALSE)

## define number of simulations and parameter scenarioå
if(doLocal) {
  scenario = 1
  N_iter = 50
}else{
  # defined from batch script params
  scenario <- as.numeric(commandArgs(trailingOnly=TRUE))
}

###############################################################
## start simulation code
###############################################################

###############################################################
## set simulation design elements
###############################################################
# load real data
load(here::here("source", "dat_func.Rdata"))

n = params$n[scenario]
family = params$family[scenario]
nS = params$nS[scenario]
SEED.START = params$seed_start[scenario]

results = vector("list", length = N_iter)
for(iter in 1:N_iter){
  # set seed
  seed.iter = (SEED.START - 1)*N_iter + iter
  #set.seed(seed.iter)

  # simulate data
  if(family %in% c("lognormal", "loglogistic")){
    sim_data <- simulate_AFT(family = family, n = n, nS = nS, seed = seed.iter)
  }else{
    sim_data <- simulate_Cox(n = n, nS = nS, seed = seed.iter)
  }

  ###############################################################
  ## fit functional AFT and Cox model
  ###############################################################
  
  ### linear functional log-normal AFT model
  time.norm <- as.numeric(system.time(fit.norm <- gam(logY ~ 1 + s(S, by = X_L, bs = "ps", k = 20), 
                                                      family = cnorm(), data = sim_data$data))[3])
  ### linear functional Cox model
  time.cox <- as.numeric(system.time(fit.cox <- gam(Y ~ s(S, by = X_L, bs = "ps", k = 20), data = sim_data$data, 
                                                        weights = delta, family = cox.ph))[3])
  
  ### calculate pointwise squared error and MISE of estimated beta(s)
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
                        cma_ub_coef_cox = as.numeric(cma.coef.cox[[2]])) %>%
    mutate(cover_coef_norm = (true_coef > lb_coef_norm) & (true_coef < ub_coef_norm),
           cover_coef_cox = (true_coef > lb_coef_cox) & (true_coef < ub_coef_cox),
           cover_cma_coef_norm = (true_coef > cma_lb_coef_norm) & (true_coef < cma_ub_coef_norm),
           cover_cma_coef_cox = (true_coef > cma_lb_coef_cox) & (true_coef < cma_ub_coef_cox))
  
  ### calculate pointwise squared error and MISE of estimated survival function
  tmax <- round(max(sim_data$data$t))
  tvec <- seq(0, tmax, length.out = 1000)
  
  # function to calculate P(T < t|X) for AFT models
  get_cdf_AFT <- function(time, lp, scale, family = "loglogistic") {
    if (family == "loglogistic") {
      p <- 1 / (1 + (exp(lp) / time)^(1 / scale))
    } else if (family == "lognormal") {
      z <- (log(time) - lp) / scale
      p <- pnorm(z)
    } 
    return(p)
  }
  
  # true survival function
  if (family %in% c("lognormal", "loglogistic")) {
    p.true <- matrix(nrow = n, ncol = length(tvec))
    for (i in 1:n) {
      for (j in 1:length(tvec)) {
        p.true[i,j] <- get_cdf_AFT(time = tvec[j], lp = sim_data$data$lp[i], scale = sim_data$coefficients$b[1], family = family)
      }
    }
  } else {
    p.true = 1 - sim_data$data$Si
  }

  # estimated survival function from AFT model
  p.est.norm <- matrix(nrow = n, ncol = length(tvec))
  for (i in 1:nrow(p.est.norm)) {
    for (j in 1:ncol(p.est.norm)) {
      p.est.norm[i,j] <- get_cdf_AFT(time = tvec[j], lp = fit.norm$linear.predictors[i], scale = fit.norm$scale, family = family)
    }
  }

  # estimated survival function from Cox model
  t0 <- rev(fit.cox$family$data$tr)
  H0_hat <- rev(fit.cox$family$data$h) 
  H0_fit <- scam(H0_hat ~ s(t0, bs = "mpi") - 1) 
  H0_prd <- pmax(0, predict(H0_fit, newdata = data.frame(t0 = tvec)))
  eta_i <- matrix(fit.cox$linear.predictors, ncol = 1)
  Si <- exp(-(exp(eta_i) %*% H0_prd))
  p.est.cox <- 1 - Si
  
  # calculate pointwise squared error and MISE
  df_surv <- data.frame(time = tvec,
                        #true_surv = 1 - p.true,
                        #est_surv_norm = 1 - p.est.norm,
                        #est_surv_cox = 1 - p.est.cox,
                        se_surv_norm = colMeans((p.true - p.est.norm)^2),
                        se_surv_cox = colMeans((p.true - p.est.cox)^2))
  #mise.surv.norm <- mean(se.surv.norm)
  #mise.surv.cox <- mean(se.surv.cox)
  
  df_info <- data.frame(scenario = scenario,
                        iter = iter,
                        seed = seed.iter,
                        n = n,
                        family = family,
                        nS = nS,
                        time_norm = time.norm,
                        time_cox = time.cox)
  
  res <- list(info = df_info, coef = df_coef, surv = df_surv)

  results[[iter]] = res

} # end for loop


filename = paste0(here::here("Output", Date), "/", scenario, ".RDA")
save(results,
     file = filename)

###############################################################
## end sim
###############################################################


