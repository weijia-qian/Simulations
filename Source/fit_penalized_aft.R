##### Function to fit the penalized AFT model #####
optimize_AFT <- function(Y,      # survival time
                         delta,  # censoring indicator
                         X,      # matrix of functional covariate
                         Z = NULL, # scalar covaraites
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
  
  # Optional scalar covariate matrix
  if (!is.null(Z)) {
    if (is.character(Z)) {
      Z_mat <- as.matrix(data[, Z, drop = FALSE])
    } else {
      Z_mat <- as.matrix(Z)
    }
  } else {
    Z_mat <- NULL
  }
  
  # Generate spline basis matrix for functional covariate
  nS <- ncol(X)
  s <- seq(0, 1, length.out = nS) 
  B <- bs(s, df = k)  # nS x k
  X_mat <- X %*% B / nS # n x k
  
  # Design matrix C: intercept + scalar covariates + functional covariate
  C <- cbind(1, Z_mat, X_mat)
  p <- ncol(C)
  
  # Initial guesses for beta and b
  init_params <- c(rep(0, p), 1)
  
  # Construct the penalty matrix
  nZ <- if (is.null(Z_mat)) 0 else ncol(Z_mat)
  Pen <- matrix(0, nrow = p, ncol = p)
  Pen_block <- penalty_matrix(kp = k, nS = nS, a = 0.001)
  Pen[(2 + nZ):p, (2 + nZ):p] <- Pen_block
  
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
  coef_all <- fit$par[1:p]
  beta0_hat <- coef_all[1]
  betaZ_hat <- if (nZ > 0) coef_all[2:(1 + nZ)] else NULL
  beta1_hat <- as.numeric(coef_all[(2 + nZ):p] %*% t(B))  # evaluate beta(s)
  b_hat <- fit$par[p + 1]
  mu_hat <- C %*% coef_all

  # beta0_hat <- fit$par[1]
  # beta1_hat <- as.numeric(fit$par[2:(k+1)] %*% t(B))
  # b_hat <- fit$par[k+2]
  # mu_hat <- C %*% fit$par[1:(k+1)]
  
  # Compute GCV
  GCV <- compute_gcv(Y, delta, C, mu_hat, b_hat, lambda, Pen, family)
  
  if (se == TRUE) {
    # Covariance matrix of beta
    hessian <- optimHess(fit$par, fn = penalized_loglik, gr = penalized_score,               
                         Y = Y, delta = delta, X = C, family = family, lambda = lambda, Pen = Pen)
    cov_beta <- solve(hessian)
    
    # Compute standard error for beta0, beta1, b
    se_beta0 <- sqrt(diag(cov_beta)[1])
    se_betaZ <- if (nZ > 0) sqrt(diag(cov_beta)[2:(1 + nZ)]) else NULL
    se_beta1 <- sqrt(rowSums((B %*% cov_beta[(2 + nZ):p, (2 + nZ):p]) * B))
    se_b <- sqrt(diag(cov_beta)[p + 1])
    
    # se_beta1 <- sqrt(rowSums((B %*% cov_beta[2:(k+1), 2:(k+1)]) * B))
    # se_b <- sqrt(diag(cov_beta)[k+2])
    
    # Confidence intervals
    beta0_ci_lower <- beta0_hat - qnorm(0.975) * se_beta0
    beta0_ci_upper <- beta0_hat + qnorm(0.975) * se_beta0
    betaZ_ci_lower <- betaZ_hat - qnorm(0.975) * se_betaZ
    betaZ_ci_upper <- betaZ_hat + qnorm(0.975) * se_betaZ
    beta1_ci_lower <- beta1_hat - qnorm(0.975) * se_beta1
    beta1_ci_upper <- beta1_hat + qnorm(0.975) * se_beta1
    b_ci_lower <- b_hat - qnorm(0.975) * se_b
    b_ci_upper <- b_hat + qnorm(0.975) * se_b
    
    return(list(beta0_hat = beta0_hat, betaZ_hat = betaZ_hat, beta1_hat = beta1_hat, b_hat = b_hat, lp = mu_hat, GCV = GCV,
                family = family, lambda = lambda, Z_names = Z,
                beta0_ci_lower = beta0_ci_lower, beta0_ci_upper = beta0_ci_upper,
                betaZ_ci_lower = betaZ_ci_lower, betaZ_ci_upper = betaZ_ci_upper,
                beta1_ci_lower = beta1_ci_lower, beta1_ci_upper = beta1_ci_upper, b_ci_lower = b_ci_lower, b_ci_upper = b_ci_upper))
  } else {
    
    return(list(beta0_hat = beta0_hat, betaZ_hat = betaZ_hat, beta1_hat = beta1_hat, b_hat = b_hat, lp = mu_hat, GCV = GCV,
                family = family, lambda = lambda, Z_names = Z))
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
    log_S <- (1 - delta) * pnorm(z, lower.tail = FALSE, log.p = TRUE)
    # log_S <- (1 - delta) * log(1 - pnorm(z))
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
    
  } else if (family == "loglogistic"){
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
optimize_lambda <- function(Y, delta, X, Z = NULL, data, family, lambda_grid) {
  gcv_values <- sapply(lambda_grid, function(lambda) optimize_AFT(Y, delta, X, Z, data, family, lambda = lambda)$GCV)
  optimal_lambda <- lambda_grid[which.min(gcv_values)]
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
                     #log_f <- -log(Y) - log(b_hat) - 0.5 * z^2 - 0.5 * log(2 * pi)
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
  
  # Clip weights for stability
  w <- pmin(pmax(w, 1e-6), 1e6)
  
  # Degrees of freedom
  W <- diag(w)
  XtWX <- t(C) %*% W %*% C
  df <- tryCatch({
    H <- solve(XtWX + lambda * Pen, XtWX)
    sum(diag(H))
  }, error = function(e) {
    warning("Matrix inversion failed during GCV calculation.")
    return(NA)
  })
  # Guard against invalid df
  if (is.na(df) || df <= 0 || df / n > 0.95) {
    return(Inf)
  }
  
  # GCV
  gcv <- -loglik / (1 - df / n)^2
  return(gcv)
}

##### Function to return the linear predictor for new obs #####
predict_AFT <- function(fit, newdata) {
  
  # extract new covariates
  X_new <- newdata$X
  if (is.null(X_new)) stop("`newdata` must have an element named `X`.")
  if (!is.matrix(X_new)) X_new <- as.matrix(X_new)
  
  # extract AFT parameters
  beta0 <- fit$beta0_hat
  beta1 <- fit$beta1_hat
  betaZ <- fit$betaZ_hat
  
  # check dimensions
  nS <- ncol(X_new)
  if (length(beta1) != nS) {
    warning(sprintf(
      "Length of beta1_hat (%d) does not match # columns of newdata$X (%d).",
      length(beta1), nS
    ))
  }
  
  # check scalar covariates
  Z_names <- fit$Z_names
  if (!is.null(betaZ)) {
    if (is.null(Z_names)) stop("Scalar covariates estimated but Z_names not found in `fit`.")
    if (!all(Z_names %in% names(newdata))) {
      missing <- Z_names[!Z_names %in% names(newdata)]
      stop("The following scalar covariates are missing in newdata: ", paste(missing, collapse = ", "))
    }
    Z_mat <- as.matrix(newdata[, Z_names, drop = FALSE])
    # compute linear predictor
    mu_new <- as.vector(beta0 + Z_mat %*% betaZ + (X_new %*% beta1) / nS)
  } else {
    mu_new <- as.vector(beta0 + (X_new %*% beta1) / nS)
  }

  
  
  # preserve row names if present
  if (!is.null(rownames(X_new))) {
    names(mu_new) <- rownames(X_new)
  }
  
  return(mu_new)
}
