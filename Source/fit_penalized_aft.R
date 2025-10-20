#' Penalized AFT with easy column selection + optional bootstrap CIs
#'
#' @param data   data.frame containing all variables
#' @param y      character(1). Column name of survival time (Y)
#' @param delta  character(1). Column name of censoring indicator (1=observed, 0=censored)
#' @param x      character. Either:
#'               - vector of exact column names for the functional covariate matrix, OR
#'               - a single regex/prefix used to grep matching columns (set x_as_regex=TRUE)
#' @param z      NULL or character vector of scalar covariate names in `data`
#' @param family "lognormal" or "loglogistic"
#' @param k      integer. # of spline basis for beta1(s)
#' @param lambda numeric. smoothing parameter
#' @param s_grid NULL or numeric vector of length = #functional grid points (optional).
#'               If NULL, assumes an equally spaced grid in [0,1].
#' @param basis     character; one of c("bs","ns","ps","tp","cr","cp")
#' @param basis_args list; forwarded to basis constructor (e.g., degree=3 for 'bs',
#'                   m=c(2,2) for 'ps', etc.)
#' @param x_as_regex logical. If TRUE and length(x)==1, select X columns via grep(x, names(data))
#' @param se         logical. Wald SE / CI via Hessian
#' @param bootstrap  logical. Row-resampling bootstrap percentile CIs
#' @param B          integer. # bootstrap replicates
#' @param boot_seed  NULL or integer seed for reproducibility
#'
#' @return list with estimates, optional Wald and bootstrap CIs, and meta info
# Requires: penalized_loglik(), penalized_score(), penalty_matrix(), compute_gcv()
optimize_AFT <- function(data,
                         y,                 # column name of survival time
                         delta,             # column name of censoring indicator (1=obs, 0=cens)
                         x,                 # vector of X colnames OR single regex/prefix (set x_as_regex=TRUE)
                         z = NULL,          # optional scalar covariate names
                         family = c("lognormal", "loglogistic"),
                         k = 20,            # # of spline basis for beta(s) (df for bs/ns)
                         lambda,            # smoothing parameter
                         s_grid = NULL,     # optional functional grid (length = ncol(X))
                         basis = c("bs", "ns"),  # splines::bs/ns
                         basis_args = list(), # forwarded to splines::bs/ns
                         x_as_regex = FALSE,
                         se = FALSE,        # Wald SE / CIs via Hessian
                         bootstrap = FALSE, # bootstrap percentile CIs
                         B = 500,           # bootstrap reps
                         boot_seed = NULL
) {  
  
  family <- match.arg(family)
  basis  <- match.arg(basis)
  
  stopifnot(is.data.frame(data))
  stopifnot(is.character(y) && length(y) == 1L)
  stopifnot(is.character(delta) && length(delta) == 1L)
  stopifnot(is.character(x) && length(x) >= 1L)
  
  # ---- pull Y, Delta ----
  if (!all(c(y, delta) %in% names(data)))
    stop("`y` and/or `delta` not found in `data`.")
  Y     <- data[[y]]
  Delta <- data[[delta]]
  
  # ---- select X (functional, wide) ----
  if (length(x) == 1L && isTRUE(x_as_regex)) {
    x_cols <- grep(x, names(data), value = TRUE)
    if (!length(x_cols)) stop("No X columns matched: ", x)
  } else {
    if (!all(x %in% names(data)))
      stop("Some X columns not found: ",
           paste(setdiff(x, names(data)), collapse = ", "))
    x_cols <- x
  }
  X  <- as.matrix(data[, x_cols, drop = FALSE])
  n  <- nrow(X)
  nS <- ncol(X)
  
  # ---- Z (scalar) ----
  if (!is.null(z)) {
    if (!all(z %in% names(data)))
      stop("Some Z columns not found: ",
           paste(setdiff(z, names(data)), collapse = ", "))
    Z_mat   <- as.matrix(data[, z, drop = FALSE])
    Z_names <- colnames(Z_mat)
  } else {
    Z_mat <- NULL
    Z_names <- NULL
  }
  nZ <- if (is.null(Z_mat)) 0L else ncol(Z_mat)
  
  # ---- s-grid ----
  if (is.null(s_grid)) {
    s <- seq(0, 1, length.out = nS)
  } else {
    if (length(s_grid) != nS) stop("`s_grid` must have length = ncol(X).")
    s <- as.numeric(s_grid)
  }
  
  # ---- build bs/ns basis matrix on s (df = k) ----
  #    basis_args is forwarded (e.g., degree=3 for bs)
  if (basis == "bs") {
    Bspl <- do.call(splines::bs, c(list(x = s, df = k), basis_args))
  } else {
    Bspl <- do.call(splines::ns, c(list(x = s, df = k), basis_args))
  }
  kb <- ncol(Bspl)  # should equal k, but keep it general
  
  # ---- projected functional design (Riemann average) ----
  X_mat <- X %*% Bspl / nS         # n x kb
  
  # ---- design matrix: intercept + Z + X_mat ----
  C <- cbind(1, Z_mat, X_mat)
  p <- ncol(C)
  
  # ---- penalty matrix on functional block ----
  Pen <- matrix(0, nrow = p, ncol = p)
  Pen_block <- penalty_matrix(kp = kb, nS = nS, a = 0.001)
  Pen[(2 + nZ):p, (2 + nZ):p] <- Pen_block
  
  # ---- optimization ----
  init_params <- c(rep(0, p), 1)   # last is scale parameter
  fit <- optim(par     = init_params,
               fn      = penalized_loglik,
               gr      = penalized_score,
               method  = "BFGS",
               Y       = Y,
               delta   = Delta,
               X       = C,
               family  = family,
               lambda  = lambda,
               Pen     = Pen,
               control = list(maxit = 2000))
  
  # ---- estimates ----
  coef_all  <- fit$par[1:p]
  beta0_hat <- coef_all[1]
  betaZ_hat <- if (nZ > 0) coef_all[2:(1 + nZ)] else NULL
  beta1_hat <- as.numeric(coef_all[(2 + nZ):p] %*% t(Bspl))  # length nS
  b_hat     <- fit$par[p + 1]
  mu_hat    <- as.vector(C %*% coef_all)
  
  # ---- GCV ----
  GCV <- compute_gcv(Y, Delta, C, mu_hat, b_hat, lambda, Pen, family)
  
  out <- list(
    beta0_hat = beta0_hat,
    betaZ_hat = betaZ_hat,
    beta1_hat = beta1_hat,
    b_hat     = b_hat,
    lp        = mu_hat,
    GCV       = GCV,
    family    = family,
    lambda    = lambda,
    Z_names   = Z_names,
    x_names   = x_cols,
    s_grid    = s,
    basis     = basis,
    kbasis    = kb
  )
  
  # ---- Wald SEs / CIs (optional) ----
  if (isTRUE(se)) {
    hessian  <- optimHess(fit$par, fn = penalized_loglik, gr = penalized_score,
                          Y = Y, delta = Delta, X = C, family = family,
                          lambda = lambda, Pen = Pen)
    cov_beta <- tryCatch(solve(hessian), error = function(e) NULL)
    
    if (is.null(cov_beta)) {
      warning("Hessian inversion failed; Wald SEs unavailable.")
    } else {
      se_beta0 <- sqrt(diag(cov_beta)[1])
      se_betaZ <- if (nZ > 0) sqrt(diag(cov_beta)[2:(1 + nZ)]) else NULL
      se_beta1 <- sqrt(rowSums((Bspl %*% cov_beta[(2 + nZ):p, (2 + nZ):p]) * Bspl))
      se_b     <- sqrt(diag(cov_beta)[p + 1])
      
      z <- qnorm(0.975)
      out$beta0_ci_lower <- beta0_hat - z * se_beta0
      out$beta0_ci_upper <- beta0_hat + z * se_beta0
      
      if (nZ > 0) {
        out$betaZ_ci_lower <- as.numeric(betaZ_hat - z * se_betaZ)
        out$betaZ_ci_upper <- as.numeric(betaZ_hat + z * se_betaZ)
        names(out$betaZ_ci_lower) <- Z_names
        names(out$betaZ_ci_upper) <- Z_names
      }
      
      out$beta1_ci_lower <- beta1_hat - z * se_beta1
      out$beta1_ci_upper <- beta1_hat + z * se_beta1
      out$b_ci_lower     <- b_hat - z * se_b
      out$b_ci_upper     <- b_hat + z * se_b
      
      out$se_beta0 <- se_beta0
      out$se_betaZ <- se_betaZ
      out$se_beta1 <- se_beta1
      out$se_b     <- se_b
    }
  }
  
  # ---- Bootstrap CIs (optional) ----
  if (isTRUE(bootstrap)) {
    if (!is.null(boot_seed)) set.seed(boot_seed)
    
    beta0_boot <- numeric(B)
    b_boot     <- numeric(B)
    beta1_boot <- matrix(NA_real_, nrow = B, ncol = nS)
    betaZ_boot <- if (nZ > 0) matrix(NA_real_, nrow = B, ncol = nZ) else NULL
    
    for (bb in seq_len(B)) {
      idx <- sample.int(n, size = n, replace = TRUE)
      
      Yb     <- Y[idx]
      Deltab <- Delta[idx]
      Xb     <- X[idx, , drop = FALSE]
      Zb     <- if (nZ > 0) Z_mat[idx, , drop = FALSE] else NULL
      
      # reuse same Bspl and penalty
      X_mat_b <- Xb %*% Bspl / nS
      Cb      <- cbind(1, Zb, X_mat_b)
      
      fit_b <- tryCatch(
        optim(par     = init_params,
              fn      = penalized_loglik,
              gr      = penalized_score,
              method  = "BFGS",
              Y       = Yb,
              delta   = Deltab,
              X       = Cb,
              family  = family,
              lambda  = lambda,
              Pen     = Pen,
              control = list(maxit = 2000)),
        error = function(e) NULL
      )
      
      if (is.null(fit_b) || fit_b$convergence != 0 || any(!is.finite(fit_b$par))) next
      
      coef_b <- fit_b$par[1:p]
      beta0_boot[bb] <- coef_b[1]
      if (nZ > 0) betaZ_boot[bb, ] <- coef_b[2:(1 + nZ)]
      beta1_boot[bb, ] <- as.numeric(coef_b[(2 + nZ):p] %*% t(Bspl))
      b_boot[bb]      <- fit_b$par[p + 1]
    }
    
    qfun <- function(x) stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    
    out$bootstrap_info <- list(
      B    = B,
      used = sum(is.finite(beta0_boot) & is.finite(b_boot))
    )
    
    beta1_ci <- apply(beta1_boot, 2, qfun)
    out$beta1_boot_ci_lower <- as.numeric(beta1_ci[1, ])
    out$beta1_boot_ci_upper <- as.numeric(beta1_ci[2, ])
    
    b_ci <- qfun(b_boot)
    out$b_boot_ci_lower <- b_ci[1]; out$b_boot_ci_upper <- b_ci[2]
    
    beta0_ci <- qfun(beta0_boot)
    out$beta0_boot_ci_lower <- beta0_ci[1]; out$beta0_boot_ci_upper <- beta0_ci[2]
    
    if (nZ > 0) {
      Z_ci <- apply(betaZ_boot, 2, qfun)
      out$betaZ_boot_ci_lower <- as.numeric(Z_ci[1, ])
      out$betaZ_boot_ci_upper <- as.numeric(Z_ci[2, ])
      names(out$betaZ_boot_ci_lower) <- Z_names
      names(out$betaZ_boot_ci_upper) <- Z_names
    }
  }
  
  out
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
  if (b <= 0 || !is.finite(b)) return(Inf)
  
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
  }
  
  # Penalized log-likelihood
  penalty <- lambda * crossprod(beta_coef, Pen) %*% beta_coef
  loglik <- sum(log_f + log_S) - penalty
  
  # Return negative log-likelihood for minimization
  if (!is.finite(loglik)) return(Inf)
  return(-loglik)
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
    score_beta <- t(X) %*% (delta * z / b + (1 - delta) * (f_z / (S_z * b))) - 2 * lambda * Pen %*% beta_coef
    
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
optimize_lambda <- function(lambda_grid,
                            data,
                            y,
                            delta,
                            x,
                            z = NULL,
                            family = c("lognormal", "loglogistic"),
                            k = 20,
                            s_grid = NULL,
                            basis = c("bs", "ns"),
                            basis_args = list(),
                            x_as_regex = FALSE) {
  
  family <- match.arg(family)
  basis  <- match.arg(basis)
  
  if (length(lambda_grid) < 1L) stop("`lambda_grid` must have length >= 1.")
  if (!is.numeric(lambda_grid)) stop("`lambda_grid` must be numeric.")
  
  # Evaluate GCV for each lambda (ignore se/bootstrap during tuning)
  gcv_values <- sapply(lambda_grid, function(lam) {
    fit_i <- try(
      optimize_AFT(
        data = data, y = y, delta = delta, x = x, z = z,
        family = family, k = k, lambda = lam,
        s_grid = s_grid, x_as_regex = x_as_regex,
        basis = basis, basis_args = basis_args,
        se = FALSE, bootstrap = FALSE,   # speed: force off during tuning
      ),
      silent = TRUE
    )
    if (inherits(fit_i, "try-error") || is.null(fit_i$GCV) || !is.finite(fit_i$GCV)) {
      return(NA_real_)
    }
    fit_i$GCV
  })
  
  # Pick the best finite GCV
  finite_idx <- which(is.finite(gcv_values))
  if (length(finite_idx) == 0L) {
    stop("All GCV evaluations failed or returned non-finite values.")
  }
  best_idx <- finite_idx[which.min(gcv_values[finite_idx])]
  optimal_lambda <- lambda_grid[best_idx]
  
  # Named vector for convenience
  names(gcv_values) <- format(lambda_grid, digits = 6)
  
  list(
    lambda = optimal_lambda,
    gcv = gcv_values,
    grid = lambda_grid
  )
}

# optimize_lambda <- function(Y, delta, X, Z = NULL, data, family, lambda_grid) {
#   gcv_values <- sapply(lambda_grid, function(lambda) optimize_AFT(Y, delta, X, Z, data, family, lambda = lambda)$GCV)
#   optimal_lambda <- lambda_grid[which.min(gcv_values)]
#   return(optimal_lambda)
# }

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
