##### Function to fit the penalized AFT model #####
optimize_AFT <- function(Y,      # survival time
                         delta,  # censoring indicator
                         X,      # matrix of functional covariate
                         data,   # name of dataset
                         family, # "lognormal" or "loglogistic"
                         k = 20, # number of spline basis to construct beta1(s)
                         lambda,  # smoothing parameter
                         se = FALSE # confidence intervals for parameters
                         ) {
  
  # Extract elements
  Y <- data$Y
  delta <- data$delta
  X <- data$X
  
  # Generate spline basis matrix
  nS <- ncol(X) # dimension of the functional predictor
  s <- seq(0, 1, length.out = nS)  # domain of s
  B <- bs(s, df = k)  # basis matrix (nS x k)
  C <- cbind(1, X %*% B / nS) # design matrix (n x (k+1))
  
  # Initial guesses for beta and b
  init_params <- c(rep(0, k+1), 1)
  
  # Construct the penalty matrix
  Pen <- penalty_matrix(kp = k+1, nS = nS, a = 0.001)
  
  # Optimization
  fit <- optim(
    par = init_params,
    fn = penalized_loglik,
    gr = penalized_score,
    method = "BFGS",
    Y = Y,
    delta = delta,
    X = C,
    family = family,
    lambda = lambda,
    Pen = Pen,
    control = list(maxit = 2000))
  
  # Extract estimates
  beta0_hat <- fit$par[1]
  beta1_hat <- as.numeric(fit$par[2:(k+1)] %*% t(B))
  b_hat <- fit$par[k+2]
  mu_hat <- C %*% fit$par[1:(k+1)]
  #RSS <- sum((log(Y) - mu_hat)^2)
  
  # Compute GCV
  beta_hat <- fit$par[1:(k+1)]
  GCV <- compute_gcv(Y, delta, C, mu_hat, b_hat, lambda, Pen, family)
  
  if (se == TRUE) {
    # Covariance matrix of beta
    hessian <- optimHess(fit$par, fn = penalized_loglik, gr = penalized_score,               
                         Y = Y, delta = delta, X = C, family = family, lambda = lambda, Pen = Pen)
    cov_beta <- solve(hessian)
    
    # Compute standard error for beta0, beta1, b
    se_beta0 <- sqrt(diag(cov_beta)[1])
    se_beta1 <- sqrt(rowSums((B %*% cov_beta[2:(k+1), 2:(k+1)]) * B))
    se_b <- sqrt(diag(cov_beta)[k+2])
    
    # Confidence intervals
    beta0_ci_lower <- beta0_hat - qnorm(0.975) * se_beta0
    beta0_ci_upper <- beta0_hat + qnorm(0.975) * se_beta0
    beta1_ci_lower <- beta1_hat - qnorm(0.975) * se_beta1
    beta1_ci_upper <- beta1_hat + qnorm(0.975) * se_beta1
    b_ci_lower <- b_hat - qnorm(0.975) * se_b
    b_ci_upper <- b_hat + qnorm(0.975) * se_b
    
    return(list(beta0_hat = beta0_hat, beta1_hat = beta1_hat, b_hat = b_hat, lp = mu_hat, GCV = GCV,
                family = family, lambda = lambda, beta0_ci_lower = beta0_ci_lower, beta0_ci_upper = beta0_ci_upper,
                beta1_ci_lower = beta1_ci_lower, beta1_ci_upper = beta1_ci_upper, b_ci_lower = b_ci_lower, b_ci_upper = b_ci_upper))
  } else {
    
    return(list(beta0_hat = beta0_hat, beta1_hat = beta1_hat, b_hat = b_hat, lp = mu_hat, GCV = GCV,
                family = family, lambda = lambda))
  }
}

##### Function to construct the penalty matrix #####
penalty_matrix <- function(kp, nS, a){
  D <- nS
  s <- seq(0, 1, length.out = nS)
  spline_basis <- bs(s, df = kp, intercept = TRUE)
  diff2 <- matrix(rep(c(1, -2, 1, rep(0, D-2)), D - 2)[1:((D-2) * D)], D-2, D, byrow = TRUE)
  P2 <- t(spline_basis) %*% t(diff2) %*% diff2 %*% spline_basis # not full rank
  Pen <- a * diag(kp) + (1-a) * P2
  Pen[1, ] <- 0 # no penalty on the intercept
  Pen[, 1] <- 0 # no penalty on the intercept
  return(Pen)
}

##### Function to compute the negative log-likelihood #####
penalized_loglik <- function(params, Y, delta, X, family, lambda, Pen) {
  # Extract parameters
  beta_coef <- params[-length(params)]
  b <- params[length(params)]
  
  # Ensure b > 0
  if (b <= 0) {
    #print("b <= 0 detected")
    return(Inf)
  }
  
  mu <- X %*% beta_coef
  z <- (log(Y) - mu) / b
  
  if (family == "lognormal"){
    log_f <- delta * (-log(Y) - log(b) - 0.5 * z^2 - 0.5 * log(2 * pi))
    log_S <- (1 - delta) * log(1 - pnorm(z))
  } else if (family == "loglogistic"){
    p_z <- 1 / (1 + exp(-z)) # logistic cdf
    log_f <- delta * log((1 / (b * Y)) * p_z * (1 - p_z))
    log_S <- (1 - delta) * log(1 - p_z)
    # # guard against zeros in log()
    # f <- pmax(f, .Machine$double.eps)
    # S <- pmax(S, .Machine$double.eps)
  }
  
  # Penalized log-likelihood
  penalty <- lambda * crossprod(beta_coef, Pen) %*% beta_coef
  loglik <- sum(log_f + log_S) - penalty
  
  # Check for invalid log-likelihood
  if (is.nan(loglik) | is.infinite(loglik)) {
    #print("Invalid log-likelihood detected")
    return(Inf)
  }
  
  return(-loglik)  # negate for minimization
}

##### Function to compute the gradient (score equations) #####
penalized_score <- function(params, Y, delta, X, family, lambda, Pen) {
  # Extract parameters
  beta_coef <- params[-length(params)]
  b <- params[length(params)]

  mu <- X %*% beta_coef
  z <- (log(Y) - mu) / b
  
  if (family == "lognormal"){
    # Log-normal components
    f_z <- dnorm(z) # PDF of z
    S_z <- pnorm(-z)  # survival function of z
    
    # Score for beta
    score_beta <- t(X) %*% (delta * z / b - (1 - delta) * (f_z / (S_z * b))) - 2 * lambda * Pen %*% beta_coef
    
    # Score for b
    score_b <- sum(delta * (-1 / b + z^2 / b) + (1 - delta) * f_z * z / (b * S_z))
    
  }
  else if (family == "loglogistic"){
    # Logistic CDF
    p_z <- 1 / (1 + exp(-z))
    
    # Score for beta
    score_beta <-  crossprod(X, delta * (2 * p_z - 1) + (1 - delta) * p_z) / b - 2 * lambda * Pen %*% beta_coef
    
    # Score for b
    score_b <- sum(- delta + delta * (2*p_z - 1) * z + (1 - delta) * p_z * z) / b
  }

  return(-c(score_beta, score_b)) # negative score for minimization
}


##### Function to find optimal lambda using GCV #####
optimize_lambda <- function(Y, delta, X, data, family, lambda_grid) {
  gcv_values <- sapply(lambda_grid, function(lambda) optimize_AFT(Y, delta, X, data, family, lambda = lambda)[5])
  optimal_lambda <- lambda_grid[which.min(gcv_values)]
  # search lambda on a finer grid
  # new_lambda_grid <- seq(optimal_lambda - 100, optimal_lambda + 100, by = 1)
  # new_gcv_values <- sapply(new_lambda_grid, function(lambda) optimize_AFT(Y, delta, X, data, family, lambda = lambda)[5])
  # optimal_lambda <- new_lambda_grid[which.min(new_gcv_values)]
  return(optimal_lambda)
}

##### Function to compute GCV value #####
compute_gcv <- function(Y, delta, C, mu_hat, b_hat, lambda, Pen, family = c("lognormal", "loglogistic")) {
  family <- match.arg(family)
  n <- length(Y)
  
  # Log-likelihood
  loglik <- switch(family,
                   "lognormal" = {
                     z <- (log(Y) - mu_hat) / b_hat
                     log_f <- dnorm(z, log = TRUE) - log(Y * b_hat)
                     log_S <- pnorm(z, lower.tail = FALSE, log.p = TRUE)
                     sum(delta * log_f + (1 - delta) * log_S)
                   },
                   "loglogistic" = {
                     z <- (log(Y) - mu_hat) / b_hat
                     log_f <- log(1 / b_hat) - log(Y) - 2 * log(1 + exp(-z))
                     log_S <- -log(1 + exp(z))
                     sum(delta * log_f + (1 - delta) * log_S)
                   }
  )
  
  # Weights for DF calculation (approximate second derivatives)
  w <- switch(family,
              "lognormal" = {
                z <- (log(Y) - mu_hat) / b_hat
                phi_z <- dnorm(z)
                S_z <- pnorm(z, lower.tail = FALSE)
                w <- numeric(n)
                w[delta == 1] <- 1 / b_hat^2
                w[delta == 0] <- (phi_z[delta == 0]^2) / (S_z[delta == 0]^2 * b_hat^2)
                w
              },
              "loglogistic" = {
                z <- (log(Y) - mu_hat) / b_hat
                # Approximate weight as derivative of log-hazard:
                p <- exp(z) / (1 + exp(z))  # derivative of log(1 + exp(z)) is logistic
                w <- numeric(n)
                w[delta == 1] <- 1 / b_hat^2  # constant approx
                w[delta == 0] <- (p[delta == 0]^2) / b_hat^2  # rough approx
                w
              }
  )
  
  # Degrees of freedom
  W <- diag(w)
  XtWX <- t(C) %*% W %*% C
  H <- solve(XtWX + lambda * Pen, XtWX)
  df <- sum(diag(H))
  
  # GCV
  gcv <- -loglik / (1 - df / n)^2
  return(gcv)
}

##### Function to return the linear predictor for new obs #####
predict_AFT <- function(fit, newdata) {
  
  # extract new X
  X_new <- newdata$X
  if (is.null(X_new)) stop("`newdata` must have an element named `X`.")
  if (!is.matrix(X_new)) X_new <- as.matrix(X_new)
  
  # extract AFT parameters
  beta0 <- fit$beta0_hat
  beta1 <- fit$beta1_hat
  
  # check dimensions
  nS <- ncol(X_new)
  if (length(beta1) != nS) {
    warning(sprintf(
      "Length of beta1_hat (%d) does not match # columns of newdata$X (%d).",
      length(beta1), nS
    ))
  }
  
  # compute linear predictor
  mu_new <- as.vector(beta0 + (X_new %*% beta1) / nS)
  
  # preserve row names if present
  if (!is.null(rownames(X_new))) {
    names(mu_new) <- rownames(X_new)
  }
  
  return(mu_new)
}