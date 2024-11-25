### sieve maximum likelihood approach by Liu et al.

est_sieve <- function(data = sim_data$data, nbasis_g = 5, nbasis_beta = 5, lower_bound = -12, upper_bound = 5){
  optim_bound <- 1000 # an option in optimization procedure not related to simulation setup
  factr <- 1e12
  a <- lower_bound
  b <- upper_bound
  rangeval_g <- c(a,b) 
  basisobj_g <- create.bspline.basis(rangeval_g, nbasis_g)
  
  sim_itr <- 1
  itr <- 1
  param_hat <- matrix(0, nrow = sim_itr, ncol = nbasis_beta + nbasis_g) # storing estimates for accessing the performance
  param_hat_se <- matrix(0, nrow =sim_itr, ncol = nbasis_beta + nbasis_g) # storing estimates for accessing the performance
  
  beta_norm_c <- rep(0, sim_itr)
  beta_norm_2 <- rep(0, sim_itr)
  g_norm_c <- rep(0, sim_itr)
  g_norm_A <- rep(0, sim_itr)
  beta_cp <- matrix(NA, nrow = 10, ncol = sim_itr)
  
  Y <- as.numeric(data$logY[, 1])
  delta <- data$delta
  z_centered <- as.matrix(scale(data$X, center = T, scale = F))
  z <- vector("list", n)
  
  for(i in 1:n){
    z[[i]] <- splinefun(z_centered[i, ])
  }
  
  # Create B-spline basis for beta(s) estimation
  rangeval_beta <- c(0, 1)
  basisobj_beta <- create.bspline.basis(rangeval_beta, nbasis_beta) # nbasis_beta is set before the loop
  
  # Calculate int phi_k(s)z_i(s)ds
  gen_z_phi <- function(z, basisobj, k){
    function(x){
      z(x) * eval.basis(x, basisobj)[, k]
    }
  }
  
  # placeholder of the integral of phi(s) * Z(s)
  # where phi(s) is the basis function of beta(s)
  # please refer to  implementation.pdf for the definition
  int_phi_z <- matrix(0, nbasis_beta, n)
  
  for (i in 1:nbasis_beta) {
    for(j in 1:n) {
      int_phi_z[i,j] <- integrate(gen_z_phi(z[[j]], basisobj_beta, i), 0, 1)$value
    }
  }
  
  ############# INITIALIZING ESTIMATION #################
  time_start <- Sys.time()
  beta_hat <- rep(0, nbasis_beta)
  g_hat <- rep(0, nbasis_g)
  
  # placeholder of the gradient of parameters
  S_theta <- rep(0, nbasis_beta + nbasis_g)
  
  prev_loss <- 10000
  
  while(TRUE){
    
    # placeholder of the gradient of parameters
    S_theta <- rep(0, nbasis_beta + nbasis_g)
    
    # objective function (refer to implementation.pdf)
    fr_a <- function(x){
      
      # alpha_hat <- x[1:3]
      beta_hat <- x[0 + (1:nbasis_beta)]
      
      # calculate r = yi - mu(Ui,theta)
      r <- Y - c(beta_hat %*% int_phi_z)
      r <- pmin(pmax(r, lower_bound + 0.1), upper_bound - 0.1)
      
      g_hat <- x[nbasis_beta + (1:nbasis_g)]   
      
      #  print(c(min(r),max(r), alpha_hat, g_hat))
      sum_g_psi <- c(eval.basis(r, basisobj_g) %*% g_hat)
      
      # generate integrant
      gen_integ_k <- function(basisobj_g, g_hat){
        function(s){
          exp(eval.basis(s, basisobj_g) %*% g_hat)
        }
      }
      
      int_exp_phi_I <- rep(0, n)
      integral <- gen_integ_k(basisobj_g, g_hat)
      
      for(i in 1:n)
        int_exp_phi_I[i] <- integrate(integral,a,r[i])$value   
      
      # you may uncomment below for showing intermediate result
      # print(c(min(r),max(r), alpha_hat, beta_hat, g_hat, -mean(sum_g_psi * delta - int_exp_phi_I)))
      # print(c(alpha_hat, -mean(sum_g_psi * delta - int_exp_phi_I)))
      
      -mean(sum_g_psi * delta - int_exp_phi_I)
    }
    
    # first derivative (refer to implementation.pdf)
    grr_a <- function(x){
      
      beta_hat <- x[0 + (1:nbasis_beta)]
      
      # calculate r = yi - mu(Ui,theta)
      r <- Y - (c(beta_hat %*% int_phi_z))
      r <- pmin(pmax(r, lower_bound + 0.1), upper_bound - 0.1)
      g_hat <- x[nbasis_beta + (1:nbasis_g)]   
      
      #  sum_g_psi_pp <- c(eval.basis(r,basisobj_g,Lfdobj=2) %*% g_hat)
      sum_g_psi_p <- c(eval.basis(r, basisobj_g, Lfdobj = 1) %*% g_hat)
      tmp <- eval.basis(r, basisobj_g) %*% g_hat
      sum_g_psi <- c(tmp)
      
      ##################################################################
      # calculate the first derivate dl/d(alpha1), dl/d(beta), dl/d(g) #
      ##################################################################
      tmp1 <- -sum_g_psi_p * delta
      tmp2 <- exp(tmp)
      tmp3 <- tmp1 + as.vector(tmp2)
      
      S_theta[0 + (1:nbasis_beta)] <- colMeans(t(int_phi_z) * tmp3)
      
      # dl/dg
      # generate integrant
      gen_integ_k <- function(basisobj_g, g_hat,k){
        function(s){
          tmp <- eval.basis(s, basisobj_g)
          exp(tmp %*% g_hat) * tmp[, k]
        }
      }
      
      int_exp_phi_I <- matrix(0, n, nbasis_g)
      
      for(i in 1:nbasis_g){
        integral <- gen_integ_k(basisobj_g, g_hat, i)
        for(j in 1:n)
          int_exp_phi_I[j,i] <- integrate(integral, a, r[j])$value   
      }
      
      S_theta[nbasis_beta + (1:nbasis_g)] <- colMeans(eval.basis(r, basisobj_g) * delta - int_exp_phi_I)
      
      -S_theta # returned value
    }
    
    param_hat_tmp <- optim(c(beta_hat, g_hat), fr_a, grr_a, 
                           method = 'L-BFGS-B', upper = optim_bound, lower = -optim_bound, 
                           control=list(factr = factr, maxit = 1000))$par
    
    # alpha_hat <- param_hat_tmp[1:2]
    beta_hat <- param_hat_tmp[1:nbasis_beta]
    g_hat    <- param_hat_tmp[nbasis_beta + (1:nbasis_g)]
    
    # placeholder of the gradient of parameters
    S_theta <- rep(0, nbasis_g)
    
    # objective function (refer to implementation.pdf)
    fr_b <- function(x){
      
      # calculate r = yi - mu(Ui,theta)
      r <- Y - c(beta_hat %*% int_phi_z)
      r <- pmin(pmax(r,lower_bound + 0.1),upper_bound - 0.1)
      
      g_hat <- x[1:nbasis_g]   
      
      sum_g_psi <- c(eval.basis(r, basisobj_g) %*% g_hat)
      
      # generate integrant
      gen_integ_k <- function(basisobj_g, g_hat){
        function(s){
          exp(eval.basis(s, basisobj_g) %*% g_hat)
        }
      }
      
      int_exp_phi_I <- rep(0, n)
      integral <- gen_integ_k(basisobj_g, g_hat)
      
      for(i in 1:n){
        int_exp_phi_I[i] <- integrate(integral, a, r[i])$value   
      }
      
      # you may uncomment below for showing intermediate result
      # print(c(min(r),max(r), alpha_hat, beta_hat, g_hat, -mean(sum_g_psi * delta - int_exp_phi_I)))
      # print(c(alpha_hat, -mean(sum_g_psi * delta - int_exp_phi_I)))
      
      -mean(sum_g_psi * delta - int_exp_phi_I)
    }
    
    # first derivative (refer to implementation.pdf)
    grr_b <- function(x){
      
      # calculate r = yi - mu(Ui,theta)
      r <- Y - (c(beta_hat %*% int_phi_z))
      r <- pmin(pmax(r, lower_bound + 0.1), upper_bound - 0.1)
      g_hat <- x[1:nbasis_g]   
      
      sum_g_psi_p <- c(eval.basis(r, basisobj_g, Lfdobj = 1) %*% g_hat)
      tmp <- eval.basis(r, basisobj_g) %*% g_hat
      sum_g_psi <- c(tmp)
      
      ##################################################################
      # calculate the first derivate dl/d(alpha1), dl/d(beta), dl/d(g) #
      ##################################################################
      tmp1 <- -sum_g_psi_p * delta
      tmp2 <- exp(tmp)
      tmp3 <- tmp1 + as.vector(tmp2)
      
      # # dl/d(alpha1)
      # S_theta[1] <- mean(x1 * tmp3)
      # 
      # # dl/d(alpha2)
      # S_theta[2] <- mean(x2 * tmp3)
      
      # dl/dg
      # generate integrant
      gen_integ_k <- function(basisobj_g, g_hat,k){
        function(s){
          tmp <- eval.basis(s, basisobj_g)
          exp(tmp %*% g_hat) * tmp[, k]
        }
      }
      
      int_exp_phi_I <- matrix(0, n,nbasis_g)
      
      for(i in 1:nbasis_g){
        integral <- gen_integ_k(basisobj_g, g_hat,i)
        for(j in 1:n)
          int_exp_phi_I[j,i] <- integrate(integral, a, r[j])$value   
      }
      
      S_theta[1:nbasis_g] <- colMeans(eval.basis(r, basisobj_g) * delta - int_exp_phi_I)
      
      -S_theta # returned value
    }
    
    # perform optimization on beta
    param_hat_tmp <- optim(c(g_hat), fr_b, grr_b, 
                           method = 'L-BFGS-B', upper = optim_bound, lower = -optim_bound, 
                           control = list(factr = factr, maxit = 1000))$par
    
    g_hat     <- param_hat_tmp[1:nbasis_g]
    # placeholder of the gradient of parameters
    S_theta <- rep(0, nbasis_beta)
    
    # objective function (refer to implementation.pdf)
    fr_g <- function(x){
      beta_hat <- x[(1:nbasis_beta)]
      
      # calculate r = yi - mu(Ui,theta)
      r <- Y - (c(beta_hat %*% int_phi_z))
      r <- pmin(pmax(r, lower_bound + 0.1), upper_bound - 0.1)
      
      sum_g_psi <- c(eval.basis(r, basisobj_g) %*% g_hat)
      
      # generate integrant
      gen_integ_k <- function(basisobj_g, g_hat){
        function(s){
          exp(eval.basis(s, basisobj_g) %*% g_hat)
        }
      }
      
      int_exp_phi_I <- rep(0,n)
      integral <- gen_integ_k(basisobj_g, g_hat)
      
      for(i in 1:n)
        int_exp_phi_I[i] <- integrate(integral, a, r[i])$value   
      
      # you may uncomment below for showing intermediate result
      # print(c(min(r),max(r), alpha_hat, beta_hat, g_hat, -mean(sum_g_psi * delta - int_exp_phi_I)))
      # print(c(alpha_hat, -mean(sum_g_psi * delta - int_exp_phi_I)))
      
      -mean(sum_g_psi * delta - int_exp_phi_I)
    }
    
    # first deriviative (refer to implementation.pdf)
    grr_g <- function(x){
      
      beta_hat <- x[(1:nbasis_beta)]
      
      # calculate r = yi - mu(Ui,theta)
      r <- Y - (c(beta_hat %*% int_phi_z))
      r <- pmin(pmax(r, lower_bound + 0.1), upper_bound - 0.1)
      sum_g_psi_p <- c(eval.basis(r, basisobj_g, Lfdobj = 1) %*% g_hat)
      tmp <- eval.basis(r,basisobj_g) %*% g_hat
      sum_g_psi <- c(tmp)
      
      ##################################################################
      # calculate the first derivate dl/d(alpha1), dl/d(beta), dl/d(g) #
      ##################################################################
      tmp1 <- -sum_g_psi_p * delta
      tmp2 <- exp(tmp)
      tmp3 <- tmp1 + as.vector(tmp2)
      
      # # dl/d(alpha0)
      # S_theta[1] <- mean(x1 * tmp3)
      # 
      # # dl/d(alpha1)
      # S_theta[2] <- mean(x2 * tmp3)
      
      S_theta[1:nbasis_beta] <- colMeans(t(int_phi_z) * tmp3)
      
      -S_theta # returned value
      
    }
    
    # perform optimization on g
    check <- optim(c(beta_hat), fr_g, grr_g, 
                   method = 'L-BFGS-B', upper = optim_bound, lower = -optim_bound, 
                   control = list(factr = factr, maxit = 1000))
    
    param_hat_tmp <- check$par
    cur_loss <- check$value
    
    beta_hat <- param_hat_tmp[1:nbasis_beta]
    
    if(abs(cur_loss - prev_loss) < 1e-4) break
    prev_loss <- cur_loss
    
  }
  
  param_hat[itr, ] <- c(beta_hat, g_hat)
  time_end <- Sys.time()
  time.sieve <- as.numeric(time_end-time_start)
  
  #write.csv(param_hat, file = paste0("../sieve/sim_data/my_param",n,"_",error_distribution,"_", tau, "_", nbasis_beta, ".csv"))
  
  #####################################
  # calculate Hessian for SE          #
  #####################################
  H_theta <- matrix(0, nbasis_beta + nbasis_g, nbasis_beta + nbasis_g)
  beta_hat <- param_hat[itr, 1:nbasis_beta]
  
  r <- Y - (c(beta_hat %*% int_phi_z))
  r <- pmin(pmax(r, lower_bound + 0.1), upper_bound - 0.1)
  g_hat <- param_hat[itr, nbasis_beta + (1:nbasis_g)] 
  
  sum_g_psi_pp <- c(eval.basis(r, basisobj_g, Lfdobj=2) %*% g_hat)
  sum_g_psi_p <- c(eval.basis(r, basisobj_g, Lfdobj=1) %*% g_hat)
  
  tmp1 <- (sum_g_psi_pp * delta - exp(eval.basis(r, basisobj_g) %*% g_hat) * sum_g_psi_p)
  
  # # d2l/d(alpha1)d(alpha1)
  # H_theta[1,1] <- mean(x1*x1*tmp1)
  # 
  # # d2l/d(alpha3)d(alpha3)
  # H_theta[2,2] <- mean(x2*x2*tmp1)
  # 
  # # d2l/d(alpha1)d(alpha2)
  # H_theta[1,2] <- mean(x1*x2*tmp1)
  # H_theta[2,1] <- H_theta[1,2]
  # 
  # # d2l/d(alpha1)d(beta)
  # H_theta[1, alpha_no+(1:nbasis_beta)] <- colMeans(t(int_phi_z) *x1*as.vector(tmp1))
  # H_theta[alpha_no + (1:nbasis_beta), 1] <- H_theta[1,alpha_no+(1:nbasis_beta)]
  # 
  # # d2l/d(alpha2)d(beta)
  # H_theta[2,alpha_no+(1:nbasis_beta)] <- colMeans(t(int_phi_z)*x2*as.vector(tmp1))
  # H_theta[alpha_no+(1:nbasis_beta),2] <- H_theta[2,alpha_no+(1:nbasis_beta)]
  
  # d2l/d(beta)d(beta)
  diag(H_theta)[1:nbasis_beta] <- colMeans(t(int_phi_z)*t(int_phi_z)*as.vector(tmp1))
  
  for(i in 1:(nbasis_beta - 1)){
    for(j in (i+1):nbasis_beta){
      H_theta[i,j] <- mean(int_phi_z[i, ] * int_phi_z[j, ] * tmp1)
      H_theta[j,i] <- H_theta[i,j]
    }
  }
  
  tmp2 <- -eval.basis(r, basisobj_g, Lfdobj = 1) * delta + as.vector(exp(eval.basis(r, basisobj_g) %*% g_hat)) * eval.basis(r, basisobj_g)
  
  # # d2l/d(alpha1)d(g)
  # H_theta[1,alpha_no+nbasis_beta+(1:nbasis_g)] <-  colMeans(x1*tmp2)
  # H_theta[alpha_no+nbasis_beta+(1:nbasis_g),1] <- H_theta[1,alpha_no+nbasis_beta+(1:nbasis_g)]
  # 
  # # d2l/d(alpha2)d(g)
  # H_theta[2,alpha_no+nbasis_beta+(1:nbasis_g)] <-  colMeans(x2*tmp2)
  # H_theta[alpha_no+nbasis_beta+(1:nbasis_g),2] <- H_theta[2,alpha_no+nbasis_beta+(1:nbasis_g)]
  
  # d2l/d(beta)d(g)
  for(i in 1:nbasis_beta){
    H_theta[i, nbasis_beta+(1:nbasis_g)] <- colMeans(c(int_phi_z[i, ]) * tmp2)
    H_theta[nbasis_beta + (1:nbasis_g), i] <- H_theta[i, nbasis_beta + (1:nbasis_g)]
  }
  
  # d2l/dgdg
  # generate integrant
  gen_integ_k_l <- function(basisobj_g, g_hat, k, l){
    function(s){
      exp(eval.basis(s, basisobj_g) %*% g_hat) * eval.basis(s, basisobj_g)[, k] * eval.basis(s, basisobj_g)[, l]
    }
  }
  
  int_exp_phi_I <- matrix(0, nbasis_g, nbasis_g)
  
  for(i in 1:nbasis_g){
    for(j in i:nbasis_g){
      integral <- gen_integ_k_l(basisobj_g, g_hat, i, j)
      for(k in 1:n)
        int_exp_phi_I[i,j] <- int_exp_phi_I[i,j] + integrate(integral, a, r[k])$value   
      int_exp_phi_I[j,i] <- int_exp_phi_I[i,j]
    }
  }
  
  H_theta[nbasis_beta + (1:nbasis_g), nbasis_beta + (1:nbasis_g)] <- -int_exp_phi_I/n
  tryCatch(                
    
    # Specifying expression
    expr = {                      
      
      param_hat_se[itr, ] <- sqrt(diag(solve(-H_theta))/n)
      
      #write.csv(param_hat_se, file = paste0("../FAFT/sim_data/my_param",n,"_",error_distribution,"_", tau, "_", nbasis_beta, "_se.csv"))
      
      #########################
      # Simultaneous CI
      #########################
      ss <- seq(0, 1, .01)
      gg <- eval.basis(ss, basisobj_beta)
      var_cov <- solve(-H_theta)[1:nbasis_beta, 1:nbasis_beta]
      
      mean_vec <- rep(0, length(ss))
      cov_mat <- gg %*% var_cov %*% t(gg)
      cov_mat <- cov_mat/(sqrt(diag(cov_mat)) %*% t(sqrt(diag(cov_mat))))
      
      Q <- mvrnorm(n = 10000,
                   mu = mean_vec,
                   Sigma = cov_mat)
      
      Q_max <- rep(0,10000)
      for(i in 1:nrow(Q)){
        Q_max[i] <- max(abs(Q[i, ]))
      }
      
      Q_n_a <- as.numeric(quantile(Q_max, c(0.95)))
      
      xi <- seq(0, .5, 0.05)
      for(k in 1:length(xi)){
        ss <- seq(xi[k], 1-xi[k], .01)
        gg <- eval.basis(ss, basisobj_beta)
        var_cov <- solve(-H_theta)[1:nbasis_beta, 1:nbasis_beta]
        mean_vec <- rep(0, length(ss))
        cov_mat <- gg %*% var_cov %*% t(gg)
        cov_mat <- cov_mat/(sqrt(diag(cov_mat))%*% t(sqrt(diag(cov_mat))))
        
        beta_hat_s <- eval.basis(ss, basisobj_beta) %*% param_hat[itr, 1:nbasis_beta]
        beta_cp[k, itr] <- (max(abs((beta_hat_s - beta(ss))/sqrt(diag(gg %*% var_cov %*% t(gg))/n))) <= Q_n_a)
      }    
    },
    
    # Specifying error message
    error = function(e){
      print("There was an error message.")
    },

    warning = function(w){
      print("There was a warning message.")
    },

    finally = {
      print("finally Executed")
    }
  )
  
  #################g C-norm###########################
  # r_hat <- Y - c(beta_hat %*% int_phi_z)
  # r_hat <- pmin(pmax(r_hat, lower_bound + 0.1), upper_bound - 0.1)
  # g_hat <- eval.basis(r_hat, basisobj_g) %*% param_hat[itr, nbasis_beta + (1:nbasis_g)]
  # 
  # r <- Y - (alpha1*x1 + alpha2*x2 + z_times_beta)
  # r <- pmin(pmax(r,lower_bound+0.1),upper_bound-0.1)
  # 
  # if(error_distribution=="weibull_1"){
  #   
  #   g_true <- log(dweibull(exp(r),shape=weibull_shape)*exp(r)/pweibull(exp(r),shape=weibull_shape, lower.tail=FALSE))
  #   
  # }else if(error_distribution=="gauss_mix_0.5"){
  #   
  #   g_true <- log(dmyMix(r)/pmyMix(r, lower.tail=FALSE))
  #   
  # }else if(error_distribution=="extrm"){
  #   
  #   g_true <- log(devd(r, scale = extrm_scale)/(1-pevd(r, scale = extrm_scale)))
  # }
  # 
  # g_norm_A[itr] <- sqrt(mean(delta*(g_hat - g_true)^2))
  
  #################beta C-norm###########################
  svec <- seq(0, 1, len = nS)
  coef.est.sieve <- eval.basis(svec, basisobj_beta) %*% param_hat[itr, 1:nbasis_beta]
  se.coef.sieve <- (coef.true - coef.est.sieve)^2
  
  return(list(coef.est.sieve, se.coef.sieve, time.sieve))
  
  # ds <- 0.0001
  # ss <- seq(0,1,ds)
  # beta_hat <- eval.basis(ss, basisobj_beta) %*% param_hat[itr, 1:nbasis_beta]
  # beta_true <- beta(ss)
  # d_beta <- beta_hat - beta_true
  # beta_norm_c[itr] <- sum(d_beta*ds)^2
  # 
  # for(i in 2:nbasis_X){
  #   beta_norm_c[itr] <- beta_norm_c[itr] + sum(cos((i-1)*pi*ss) * d_beta * ds)^2 * 2 / i
  # }
  # beta_norm_c[itr] <- sqrt(beta_norm_c[itr]/3*(aa^2))
  # beta_norm_2[itr] <- sqrt(sum(d_beta^2 * ds))
  ##################################################
  
  
  # if(itr > 1){
  #   print('Beta CP')
  #   # print(paste(mean(beta_cp0[1:itr]),mean(beta_cp1[1:itr]),mean(beta_cp2[1:itr]),mean(beta_cp3[1:itr]),mean(beta_cp4[1:itr])))
  #   print(rowMeans(beta_cp[, 1:itr], na.rm = TRUE))
  #   
  #   print('current norm_c and norm_2')
  #   print(paste(beta_norm_c[itr],beta_norm_2[itr],g_norm_c[itr],g_norm_A[itr]))
  #   print(paste(mean(beta_norm_c[1:itr]),mean(beta_norm_2[1:itr]),mean(g_norm_c[1:itr],na.rm=TRUE),mean(g_norm_A[1:itr],na.rm=TRUE)))
  #   # print(paste(mean(beta_norm_c[1:itr]^2),mean(beta_norm_2[1:itr]^2),mean(g_norm_c[1:itr]^2),mean(g_norm_2[1:itr]^2)))
  #   
  #   write.csv(c(colMeans(abs(param_hat[1:itr,1:2] - c(alpha1,alpha2))),
  #               colMeans(param_hat_se[1:itr,1:2], na.rm = TRUE),
  #               mean(abs(param_hat[1:itr,1] - alpha1) > 1.96*param_hat_se[1:itr,1], na.rm = TRUE),
  #               mean(abs(param_hat[1:itr,2] - alpha2) > 1.96*param_hat_se[1:itr,2], na.rm = TRUE),
  #               mean(beta_norm_c[1:itr]),mean(beta_norm_2[1:itr]),mean(g_norm_c[1:itr],na.rm=TRUE),mean(g_norm_A[1:itr],na.rm=TRUE),
  #               rowMeans(beta_cp[,1:itr],na.rm=TRUE)), 
  #             file = paste0("../FAFT/sim_data/curve",n,"_",error_distribution,"_", tau, "_", nbasis_beta, "_alpha.csv"))
  # }
  
  }
