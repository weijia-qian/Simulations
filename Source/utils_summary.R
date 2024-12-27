get_CMA <- function(fit){
  # call predict.gam
  coef_est <- predict(fit, newdata = df_pred, type = "terms", se.fit = TRUE) 
  # extract point estimates and pointwise standard errors 
  beta_hat <- coef_est$fit[, 1] 
  se_beta_hat <- coef_est$se.fit[,1]
  # get the design matrix associated with sgrid
  lpmat <- predict(fit, newdata = df_pred, type = "lpmatrix")
  # get the column indices associated with the functional term
  inx_beta <- which(grepl("s\\(S\\):X_L\\.[0-9]+", dimnames(lpmat)[[2]]))
  # get the design matrix A associated with beta(s)
  Bmat <- lpmat[,inx_beta]
  # get Var spline coefficients (beta_sp) associated with beta(u,s)
  beta_sp <- coef(fit)[inx_beta] 
  Vbeta_sp <- vcov(fit)[inx_beta, inx_beta]
  # number of bootstrap samples (B)
  nboot <- 1e4 
  # set up container for bootstrap 
  beta_mat_boot <- matrix(NA, nboot, length(sgrid))
  #Do the bootstrap 
  for(i in 1:nboot){ 
    beta_sp_i <- MASS::mvrnorm(n = 1, mu = beta_sp, Sigma = Vbeta_sp) 
    beta_mat_boot[i,] <- Bmat %*% beta_sp_i 
  }
  #Find the max statistic 
  dvec <- apply(beta_mat_boot, 1, function(x) max(abs(x - beta_hat) / se_beta_hat))
  #Get 95% global confidence band 
  Z_global <- quantile(dvec, 0.95)
  beta_hat_LB_global <- beta_hat - Z_global * se_beta_hat
  beta_hat_UB_global <- beta_hat + Z_global * se_beta_hat
  
  return(list(beta_hat_LB_global, beta_hat_UB_global))
}
