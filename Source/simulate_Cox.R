### Simulate survival data from Cox model

simulate_Cox = function(data = dat_func,
                        beta_type = "simple",
                        n = 500,
                        npc = 5,
                        tmax = 300,
                        range_s = 1,
                        nS = 401,
                        u = 140,
                        seed = 916){
  
  ### 1. Simulate new curves Xi(s)
  # FPCA on real data
  df_wide <- pivot_wider(subset(data, select = -seconds), names_from = frame, values_from = percent_change, 
                         names_prefix = "percent_change_")
  s <- unique(data$seconds)   # original grid
  sgrid <- seq(0, range_s, length.out = nS)  # simulation grid
  if (nS <= length(s)) {
    # downsampling the original functional domain
    matrix_wide <- as.matrix(df_wide[, -(1:8)])[, seq(1, length(s), length.out = nS)] 
  } else {
    # interpolate each row to have nS evenly spaced points
    matrix_wide <- apply(as.matrix(df_wide[, -(1:8)]), 1, function(row) {
      approx(seq(1, length(s)), row, xout = seq(1, length(s), length.out = nS))$y
    })
    matrix_wide <- t(matrix_wide)
  }
  fpca.results <- fpca.face(matrix_wide, argvals = sgrid, pve = 0.99)
  
  # re-scaled to appear on the scale of the original functions 
  Phi <- sqrt(nS) * fpca.results$efunctions
  eigenvalues <- fpca.results$evalues / nS
  
  # simulate new scores from MVN
  set.seed(seed)
  sim_scores <- mvrnorm(n = n, mu = rep(0, npc), Sigma = diag(eigenvalues[1:npc]))
  mu <- matrix(rep(fpca.results$mu, n), nrow = n, byrow = TRUE)
  
  # generate new curves based on FPCs and new scores
  sim_curves <- mu + sim_scores %*% t(Phi[, 1:npc])
  sim_curves[,1] <- 0 # asgridign the initial value to zero
  
  ### 2. Fit a linear functional Cox model on real data
  df_wide$X = I(matrix_wide)
  # quadrature weights for Riemann integration
  lvec <- matrix(1 / nS, nS, 1) 
  # matrix containing quadrature weights for all participants
  L <- kronecker(matrix(1, nrow(df_wide), 1), t(lvec))
  # matrix containing functional domain values
  S <- kronecker(matrix(1, nrow(df_wide), 1), t(sgrid))
  df_wide$S <- I(S)
  # pointwise product of w_i(s) and the quadrature weights, "L"
  df_wide$X_L <- I(df_wide$X * L)
  fit <- gam(time_delay ~ s(S, by = X_L, bs = "bs", k = 30), data = df_wide, 
             weights = use_num, family = cox.ph)

  ### 3. Estimate the cumulative baseline hazard
  # derive raw estimates 
  t0 <- rev(fit$family$data$tr) 
  H0_hat <- rev(fit$family$data$h) 
  # smooth while imposing non-decreasing shape constraints 
  H0_fit <- scam(H0_hat ~ s(t0, bs = "mpi") - 1) 
  # set the time grid to evaluate cumulative baseline hazard 
  tgrid <- seq(0, tmax, len = 1000) 
  # derive final estimates on the grid 
  H0_prd <- pmax(0, predict(H0_fit, newdata = data.frame(t0 = tgrid)))

  ### 4. Estimate the linear predictor 
  df_sim <- data.frame(X = I(sim_curves), 
                       L = I(matrix(1 / nS, ncol = nS, nrow = n)), 
                       S = I(matrix(sgrid, ncol = nS, nrow = n, byrow = TRUE)))
  df_sim$X_L = I(df_sim$X * df_sim$L)
  if (beta_type == "simple"){
    beta <- function(s) 0.3 - (s - 0.2)^2
    eta_i <- sim_curves %*% beta(sgrid) * (range_s / nS)
  } else {
    #beta <- function(s) -0.3 + cos(10 * s)
    eta_i <- predict(fit, newdata = df_sim, type = "terms")
  }
  
  ### 5. Estimate the survival function 
  Si <- exp(-(exp(eta_i) %*% H0_prd))

  ### 6. Simulate survival times 
  set.seed(seed)
  U <- runif(n) 
  Ti <- rep(NA, n) 
  for(i in 1:n){ 
    if(all(Si[i,] > U[i])){Ti[i] <- max(tgrid) + 1} 
    else{Ti[i] <- tgrid[min(which(Si[i,] < U[i]))]} 
  }
  summary(Ti)

  ### 7. Simulate censoring times from uniform (0, u)
  set.seed(seed)
  Ci <- runif(n, 60, u)
  Yi <- pmin(Ci, Ti) # observed time to event 
  di <- as.numeric(Ti <= Ci) # binary event indicator

  # set up the AFT response
  logY <- cbind(log(Yi), log(Yi))
  logY[di == 0, 2] <- Inf # right censoring

  # save simulated data
  sim_data_wide = data.frame(ID = seq(1:n), 
                       Y = Yi,
                       t = Ti,
                       C = Ci,
                       delta = di,
                       X = I(sim_curves),
                       X_L = I(df_sim$X_L),
                       S = I(df_sim$S),
                       lp = eta_i,
                       logY = I(logY),
                       Si = I(Si))
  
  # save true coefficient functions
  if (beta_type == "simple"){
    df_coef = data.frame(time = sgrid,
                         beta1 = beta(sgrid)
                         )
  } else {
    df_coef = data.frame(time = sgrid,
                         beta1 = as.numeric(predict(fit, newdata = data.frame(S = sgrid, X_L = 1), type = "terms"))
                         )
  }
  
  return(list(data = sim_data_wide, coefficients = df_coef, family = "cox.ph"))

}
