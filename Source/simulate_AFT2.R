simulate_AFT2 = function(data = dat_func, 
                        family = "lognormal",
                        n = 500,
                        npc = 5,
                        smax = 1,
                        nS = 401,
                        beta_type = "linear",
                        beta_0 = 0.5,
                        b = 0.1,
                        u = 1600,
                        seed = 916){
  
  # FPCA on real data
  df_wide <- pivot_wider(subset(data, select = -seconds), names_from = frame, values_from = percent_change, 
                         names_prefix = "percent_change_")
  s <- unique(data$seconds)   # original grid
  sgrid <- seq(0, smax, length.out = nS)  # simulation grid
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
  sim_curves[, 1] <- 0 # assign the initial value to zero
  #sim_curves <- sweep(sim_curves, 2, colMeans(sim_curves)) ## mean center simulated functional covariates
  
  # observed points on the functional domain
  sgrid <- seq(0, smax, length.out = nS)
  
  ## define pre-specified functional coefficient surface F(.,.)
  if(beta_type == "linear"){
    beta <- function(X, s) -0.5*X*s
  }else if(beta_type == "const"){
    beta <- function(X, s) -X*0.1
  }else if(beta_type == "nonlinear1"){
    beta <- function(X, s) (0.1*X^2)*s
  }else if(beta_type == "nonlinear2"){
    beta <- function(X, s) sin(-X*s)
  }

  # numeric integral over the grid given pre-specified F
  num_int <- as.matrix(apply(sim_curves, 1, function(x) sum(beta(X = x, s = sgrid))* (smax / nS)), ncol = 1)
  
  # simulate z_i from logistic(0,1)
  if (family == "loglogistic"){
    z <- rlogis(n)
  } else if (family == "lognormal"){
    #z <- rnorm(n)
    # use Box-Muller transform
    u1 <- runif(n)
    u2 <- runif(n)
    z <- sqrt(-2 * log(u1)) * cos(2 * pi * u2)
  } else {
    stop('Invalid input!')
  }
  
  # obtain true survival times (T)
  lp <- rep(beta_0, n) + num_int
  t <- exp(lp + b * z)
  
  # simulate censoring times from uniform (0, u)
  C <- runif(n, 0, u)
  
  # obtain observed survival times (Y)
  Y <- pmin(t, C)
  
  # obtain censoring indicator
  delta <- ifelse(Y == t, 1, 0)

  # set up the AFT response
  logY <- cbind(log(Y), log(Y))
  logY[delta == 0, 2] <- Inf # right censoring
  
  # set up data structure for mgcv fitting
  lvec <- matrix(1 / nS, nS, 1) # quadrature weights for Riemann integration
  L <- kronecker(matrix(1, n, 1), t(lvec)) # matrix containing quadrature weights for all participants
  S <- kronecker(matrix(1, n, 1), t(sgrid)) # matrix containing functional domain values
  
  # save simulated data
  sim_data_wide = data.frame(ID = seq(1:n), 
                             Y = Y,
                             t = t,
                             C = C,
                             delta = delta,
                             X = I(sim_curves),
                             L = I(L),
                             X_L = I(sim_curves * L),
                             S = I(S),
                             lp = lp,
                             logY = I(logY))
  
  # save true coefficient functions
  # df_coef = data.frame(time = sgrid,
  #                      beta0 = rep(beta_0, nS),
  #                      beta1 = beta(sgrid),
  #                      b = b)
  
  return(list(data = sim_data_wide, beta_type = beta_type, family = family, beta_0 = beta_0, beta = beta, b = b))
  
}