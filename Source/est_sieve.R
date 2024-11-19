### sieve maximum likelihood approach by Liu et al.
# install.packages('fda')
# library(fda)

# This setup follows:
#
# A SIEVE M-THEOREM FOR BUNDLED PARAMETERS IN
# SEMIPARAMETRIC MODELS, WITH APPLICATION TO THE
# EFFICIENT ESTIMATION IN A LINEAR MODEL FOR
# CENSORED DATA
# Y = 2 + X1 + X2 + integral of beta(s)*z(s)
simulation <- function(Y, delta, n, nbasis_g, nbasis_beta, tau, error_distribution, lower_bound, upper_bound){
  
  extrm_scale <- 2
  weibull_shape <- 0.6
  aa <- 3
  alpha_no <- 2
  
  nbasis_X <- 50
  set.seed(2024)
  optim_bound <- 1000 # an option in optimization procedure not related to simulation setup
  factr <- 1e12
  
  
  alpha1 <- 1
  alpha2 <- 1
  
  
  a <- lower_bound
  b <- upper_bound
  rangeval_g <- c(a,b) 
  basisobj_g <- create.bspline.basis(rangeval_g, nbasis_g)
  
  
  param_hat <- matrix(0,nrow=sim_itr,ncol=2+nbasis_beta+nbasis_g) # storing estimates for accessing the performance
  param_hat_se <- matrix(0,nrow=sim_itr,ncol=2+nbasis_beta+nbasis_g) # storing estimates for accessing the performance
  
  beta_norm_c <- rep(0,sim_itr)
  beta_norm_2 <- rep(0,sim_itr)
  g_norm_c <- rep(0,sim_itr)
  g_norm_A <- rep(0,sim_itr)
  
  beta_cp <- matrix(NA,nrow=10,ncol=sim_itr)
  

    ############# START GENERATING DATA (delta, X1, X2, Z, Censored T) #################
    # x1 <- rbinom(n,1,0.5) 
    x1 <- rbinom(n,1,0.5) - 0.5
    # x2 <- pmin(pmax(rnorm(n,0,sqrt(0.5)),-2),2)
    x2 <- pmin(pmax(rnorm(n,0,0.5),-2),2)
    # generating functional covariates according to
    # JOS Optimal Estimation for the Functional Cox Model (P.12) 
    # by Simeng Qu, Jane-Ling Wang and Xiao Wang
    gen_z <- function(v)
    {
      u <- runif(50,-3,3)
      
      function(s){
        
        k <- 1
        xi <- (-1)^(k+1)*k^(-v/2)
        z <- xi*u[k]*rep(1,length(s))
        
        for(k in 2:50){
          
          xi <- (-1)^(k+1)*k^(-v/2)
          z <- z + xi*u[k]*sqrt(2)*cos((k-1)*pi*s) * ((1/(1+exp(-1000*s))-0.5)*2)* ((1/(1+exp(-1000*(-(s-1))))-0.5)*2)
        }
        z
      }
    }
    
    # product of Z(s) and beta(s)
    gen_z_beta <- function(z,beta){
      
      function(s){
        
        z(s)*beta(s)
      }
    }
    
    z <- vector("list",n)
    z_times_beta <- rep(0,n) # vector of integral of Z(s) * beta(s) 
    
    for(i in 1:n){
      z[[i]] <- gen_z(v)
      z_times_beta[i] <- integrate(gen_z_beta(z[[i]],beta),0,1)$value
    }
    
    
    # error term (normal distribution)
    # e <- rnorm(n,0,1)
    if(error_distribution == "weibull_1"){
      
      # e <- log(rexp(n, rate=1)) # (i) exp(e) ~ exponential(1)
      e <- log(rweibull(n, shape=weibull_shape))
      
    }else if(error_distribution == "gauss_mix_0.5"){
      
      e <- rmyMix(n) # (i) exp(e) ~ exponential(1)
      
    }else if(error_distribution == "extrm"){
      
      e <- revd(n, scale = extrm_scale) # (i) exp(e) ~ exponential(1)
    }
    
    # event time
    fail_T <- alpha1*x1 + alpha2*x2 + z_times_beta + e
    
    print(paste0("itr: ", itr))
    censoring_T <- log(runif(n,0,tau))  # censoring time (uniform from 0 to tau)
    Y <- pmin(fail_T,censoring_T)  # observation time
    delta <- fail_T < censoring_T  # delta
    
    
    Y <- sim_data$data$logY
    delta <- sim_data$data$delta
    
    # Create B-spline basis for beta(s) estimation
    rangeval_beta <- c(0,1)
    basisobj_beta <- create.bspline.basis(rangeval_beta, nbasis_beta) # nbasis_beta is set before the loop
    
    # Calculate int phi_k(s)z_i(s)ds
    gen_z_phi <- function(z, basisobj, k){
      function(s){
        z(s)*eval.basis(s,basisobj)[,k]
      }
    }
    
    # placeholder of the integral of phi(s) * Z(s)
    # where phi(s) is the basis function of beta(s)
    # please refer to  implementation.pdf for the definition
    int_phi_z <- matrix(0, nbasis_beta, n)
    
    for(i in 1:nbasis_beta)
      for(j in 1:n)
        int_phi_z[i,j] <- integrate(gen_z_phi(z[[j]],basisobj_beta,i),0,1)$value        
    
    
    
    # Define parameters
    s <- seq(0, 1, length.out = 500)  # equally spaced time points on [0,1]
    df <- 5  # degrees of freedom (number of basis functions)
    
    # Generate B-spline basis matrix
    basis <- bs(s, df = df, intercept = TRUE)  # Returns a matrix of size (500 x df)
    
    my_int_phi_z <- matrix(0, n, nbasis_beta)
    for (i in 1:n) {
      # Element-wise product of f_i(s) with each basis function phi_k(s)
      my_int_phi_z[i, ] <- colSums(basis * sim_data$data$X[i, ] * (1 / 500))  # 1/500 is delta_s for trapezoidal rule
    }
    int_phi_z <- t(my_int_phi_z)
    
    ############# INITIALIZING ESTIMATION #################
#    alpha_hat <- rep(0,2)
    beta_hat <- rep(0,nbasis_beta)
    g_hat     <- rep(0,nbasis_g)
    
    # placeholder of the gradient of parameters
    S_theta <- rep(0, nbasis_beta+nbasis_g)
    
    # param_hat[itr,] <- c(alpha_hat, beta_hat, g_hat)
    prev_loss <- 10000
    
    while(TRUE){
      
      # placeholder of the gradient of parameters
      S_theta <- rep(0, nbasis_beta + nbasis_g)
      
      
      # objective function (refer to implementation.pdf)
      fr_a <- function(x){
        
        # alpha_hat <- x[1:3]
        beta_hat <- x[0+(1:nbasis_beta)]
        
        # calculate r = yi - mu(Ui,theta)
        r <- Y - (c(beta_hat %*% int_phi_z))
        r <- pmin(pmax(r,lower_bound+0.1),upper_bound-0.1)
        # print(paste0('r: ',min(r),' ',max(r)))
        # print(c(alpha_hat,beta_hat,g_hat))
        
        g_hat <- x[0+nbasis_beta+(1:nbasis_g)]   
        
        #  print(c(min(r),max(r), alpha_hat, g_hat))
        sum_g_psi <- c(eval.basis(r,basisobj_g) %*% g_hat)
        
        # generate integrant
        gen_integ_k <- function(basisobj_g,g_hat){
          function(s){
            exp(eval.basis(s,basisobj_g) %*% g_hat)
          }
        }
        
        int_exp_phi_I <- rep(0,n)
        integral <- gen_integ_k(basisobj_g,g_hat)
        
        for(i in 1:n)
          int_exp_phi_I[i] <- integrate(integral,a,r[i])$value   
        
        # you may uncomment below for showing intermediate result
        # print(c(min(r),max(r), alpha_hat, beta_hat, g_hat, -mean(sum_g_psi * delta - int_exp_phi_I)))
        # print(c(alpha_hat, -mean(sum_g_psi * delta - int_exp_phi_I)))
        
        -mean(sum_g_psi * delta - int_exp_phi_I)
      }
      
      # first derivative (refer to implementation.pdf)
      grr_a <- function(x){
        
        # alpha_hat <- x[1:3]
        beta_hat <- x[0+(1:nbasis_beta)]
        
        # calculate r = yi - mu(Ui,theta)
        r <- Y - (c(beta_hat %*% int_phi_z))
        r <- pmin(pmax(r,lower_bound+0.1),upper_bound-0.1)
        g_hat <- x[0+nbasis_beta+(1:nbasis_g)]   
        
        
        #  sum_g_psi_pp <- c(eval.basis(r,basisobj_g,Lfdobj=2) %*% g_hat)
        sum_g_psi_p <- c(eval.basis(r,basisobj_g,Lfdobj=1) %*% g_hat)
        tmp <- eval.basis(r,basisobj_g) %*% g_hat
        sum_g_psi <- c(tmp)
        
        ##################################################################
        # calculate the first derivate dl/d(alpha1), dl/d(beta), dl/d(g) #
        ##################################################################
        tmp1 <- -sum_g_psi_p * delta
        tmp2 <- exp(tmp)
        tmp3 <- tmp1 + as.vector(tmp2)
        
        S_theta[0+(1:nbasis_beta)] <- colMeans(t(int_phi_z) * tmp3)
        
        # dl/dg
        # generate integrant
        gen_integ_k <- function(basisobj_g,g_hat,k){
          function(s){
            tmp <- eval.basis(s,basisobj_g)
            exp(tmp %*% g_hat) * tmp[,k]
          }
        }
        
        int_exp_phi_I <- matrix(0,n,nbasis_g)
        
        for(i in 1:nbasis_g){
          integral <- gen_integ_k(basisobj_g,g_hat,i)
          for(j in 1:n)
            int_exp_phi_I[j,i] <- integrate(integral,a,r[j])$value   
        }
        
        S_theta[0+nbasis_beta+(1:nbasis_g)] <- colMeans(eval.basis(r,basisobj_g) * delta - int_exp_phi_I)
        
        -S_theta # returned value
        
      }
      
      param_hat_tmp <- optim(c(beta_hat,g_hat), fr_a, grr_a, 
                             method='L-BFGS-B', upper=optim_bound, lower=-optim_bound, 
                             control=list(factr =factr, maxit=1000))$par
      
      # alpha_hat <- param_hat_tmp[1:2]
      beta_hat <- param_hat_tmp[0+(1:nbasis_beta)]
      g_hat     <- param_hat_tmp[0+nbasis_beta+(1:nbasis_g)]
      
      # placeholder of the gradient of parameters
      S_theta <- rep(0,2+0+nbasis_g)
      
      # objective function (refer to implementation.pdf)
      fr_b <- function(x){
        
#        alpha_hat <- x[1:2]
        # beta_hat <- x[0+(1:nbasis_beta)]
        
        # calculate r = yi - mu(Ui,theta)
        r <- Y - (c(beta_hat %*% int_phi_z))
        r <- pmin(pmax(r,lower_bound+0.1),upper_bound-0.1)
        # print(paste0('r: ',min(r),' ',max(r)))
        # print(c(alpha_hat,beta_hat,g_hat))
        
        g_hat <- x[0+(1:nbasis_g)]   
        
        #  print(c(min(r),max(r), alpha_hat, g_hat))
        sum_g_psi <- c(eval.basis(r,basisobj_g) %*% g_hat)
        
        # generate integrant
        gen_integ_k <- function(basisobj_g,g_hat){
          function(s){
            exp(eval.basis(s,basisobj_g) %*% g_hat)
          }
        }
        
        int_exp_phi_I <- rep(0,n)
        integral <- gen_integ_k(basisobj_g,g_hat)
        
        for(i in 1:n)
          int_exp_phi_I[i] <- integrate(integral,a,r[i])$value   
        
        # you may uncomment below for showing intermediate result
        # print(c(min(r),max(r), alpha_hat, beta_hat, g_hat, -mean(sum_g_psi * delta - int_exp_phi_I)))
        # print(c(alpha_hat, -mean(sum_g_psi * delta - int_exp_phi_I)))
        
        -mean(sum_g_psi * delta - int_exp_phi_I)
      }
      
      # first derivative (refer to implementation.pdf)
      grr_b <- function(x){
        
#        alpha_hat <- x[1:2]
        # beta_hat <- x[3+(1:nbasis_beta)]
        
        # calculate r = yi - mu(Ui,theta)
        r <- Y - (c(beta_hat %*% int_phi_z))
        r <- pmin(pmax(r,lower_bound+0.1),upper_bound-0.1)
        g_hat <- x[0+(1:nbasis_g)]   
        
        #  sum_g_psi_pp <- c(eval.basis(r,basisobj_g,Lfdobj=2) %*% g_hat)
        sum_g_psi_p <- c(eval.basis(r,basisobj_g,Lfdobj=1) %*% g_hat)
        tmp <- eval.basis(r,basisobj_g) %*% g_hat
        sum_g_psi <- c(tmp)
        
        ##################################################################
        # calculate the first derivate dl/d(alpha1), dl/d(beta), dl/d(g) #
        ##################################################################
        tmp1 <- -sum_g_psi_p * delta
        tmp2 <- exp(tmp)
        tmp3 <- tmp1 + as.vector(tmp2)
        
        # dl/d(alpha1)
        S_theta[1] <- mean(x1 * tmp3)
        
        # dl/d(alpha2)
        S_theta[2] <- mean(x2 * tmp3)
        
        # dl/dg
        # generate integrant
        gen_integ_k <- function(basisobj_g,g_hat,k){
          function(s){
            tmp <- eval.basis(s,basisobj_g)
            exp(tmp %*% g_hat) * tmp[,k]
          }
        }
        
        int_exp_phi_I <- matrix(0,n,nbasis_g)
        
        for(i in 1:nbasis_g){
          integral <- gen_integ_k(basisobj_g,g_hat,i)
          for(j in 1:n)
            int_exp_phi_I[j,i] <- integrate(integral,a,r[j])$value   
        }
        
        S_theta[0+(1:nbasis_g)] <- colMeans(eval.basis(r,basisobj_g) * delta - int_exp_phi_I)
        
        -S_theta # returned value
      }
      
      # perform optimization on beta
      # param_hat_tmp <- optim(c(alpha_hat,g_hat), fr_b, grr_b, 
      #                        method='L-BFGS-B', upper=optim_bound, lower=-optim_bound)$par
      param_hat_tmp <- optim(c(alpha_hat,g_hat), fr_b, grr_b, 
                             method='L-BFGS-B', upper=optim_bound, lower=-optim_bound, 
                             control=list(factr =factr, maxit=1000))$par
      
      
      
#      alpha_hat <- param_hat_tmp[1:2]
      g_hat     <- param_hat_tmp[0+(1:nbasis_g)]
      # placeholder of the gradient of parameters
      S_theta <- rep(0,nbasis_beta+0)
      
      # objective function (refer to implementation.pdf)
      fr_g <- function(x){
        
#        alpha_hat <- x[1:2]
        beta_hat <- x[(1:nbasis_beta)]
        
        # calculate r = yi - mu(Ui,theta)
        r <- Y - (c(beta_hat %*% int_phi_z))
        r <- pmin(pmax(r,lower_bound+0.1),upper_bound-0.1)
        # print(paste0('r: ',min(r),' ',max(r)))
        # print(c(alpha_hat,beta_hat,g_hat))
        
        #  print(c(min(r),max(r), alpha_hat, g_hat))
        sum_g_psi <- c(eval.basis(r,basisobj_g) %*% g_hat)
        
        # generate integrant
        gen_integ_k <- function(basisobj_g,g_hat){
          function(s){
            exp(eval.basis(s,basisobj_g) %*% g_hat)
          }
        }
        
        int_exp_phi_I <- rep(0,n)
        integral <- gen_integ_k(basisobj_g,g_hat)
        
        for(i in 1:n)
          int_exp_phi_I[i] <- integrate(integral,a,r[i])$value   
        
        # you may uncomment below for showing intermediate result
        # print(c(min(r),max(r), alpha_hat, beta_hat, g_hat, -mean(sum_g_psi * delta - int_exp_phi_I)))
        # print(c(alpha_hat, -mean(sum_g_psi * delta - int_exp_phi_I)))
        print(paste("r: ",c(min(r),max(r))))
        
        -mean(sum_g_psi * delta - int_exp_phi_I)
      }
      
      # first deriviative (refer to implementation.pdf)
      grr_g <- function(x){
        
#        alpha_hat <- x[1:2]
        beta_hat <- x[(1:nbasis_beta)]
        
        # calculate r = yi - mu(Ui,theta)
        r <- Y - (c(beta_hat %*% int_phi_z))
        r <- pmin(pmax(r,lower_bound+0.1),upper_bound-0.1)
        sum_g_psi_p <- c(eval.basis(r,basisobj_g,Lfdobj=1) %*% g_hat)
        tmp <- eval.basis(r,basisobj_g) %*% g_hat
        sum_g_psi <- c(tmp)
        
        ##################################################################
        # calculate the first derivate dl/d(alpha1), dl/d(beta), dl/d(g) #
        ##################################################################
        tmp1 <- -sum_g_psi_p * delta
        tmp2 <- exp(tmp)
        tmp3 <- tmp1 + as.vector(tmp2)
        
        # dl/d(alpha0)
 #       S_theta[1] <- mean(x1 * tmp3)
        
        # dl/d(alpha1)
 #       S_theta[2] <- mean(x2 * tmp3)
        
        S_theta[(1:nbasis_beta)] <- colMeans(t(int_phi_z) * tmp3)
        
        -S_theta # returned value
        
      }
      
      # perform optimization on g
      # check <- optim(c(alpha_hat,beta_hat), fr_g, grr_g, 
      #                method='L-BFGS-B', upper=optim_bound, lower=-optim_bound)
      check <- optim(c(alpha_hat,beta_hat), fr_g, grr_g, 
                     method='L-BFGS-B', upper=optim_bound, lower=-optim_bound, 
                     control=list(factr =factr, maxit=1000))
      
      param_hat_tmp <- check$par
      cur_loss <- check$value
      print(c(prev_loss,cur_loss))
      
#      alpha_hat <- param_hat_tmp[1:2]
      beta_hat <- param_hat_tmp[(1:nbasis_beta)]
      
      if(abs(cur_loss-prev_loss) < 1e-4) break
      prev_loss <- cur_loss
      
      # print(c(alpha_hat, beta_hat, g_hat))
    }
    
    param_hat[itr,] <- c(alpha_hat, beta_hat, g_hat)
    
    # show current estimation
    print(param_hat[itr,])
    
    if(itr > 1){
      # show average estimation so far
      print(c(colMeans(param_hat[1:itr,]),itr))
    }
    
    write.csv(param_hat, file = paste0("../FAFT/sim_data/my_param",n,"_",error_distribution,"_", tau, "_", nbasis_beta, ".csv"))
    
    #####################################
    # calculate Hessian for SE          #
    #####################################
    H_theta <- matrix(0, nbasis_beta+nbasis_g, nbasis_beta+nbasis_g)
#    alpha_hat <- param_hat[itr,1:alpha_no]
    beta_hat <- param_hat[itr,(1:nbasis_beta)]
    
    r <- Y - (c(beta_hat %*% int_phi_z))
    r <- pmin(pmax(r,lower_bound+0.1),upper_bound-0.1)
    g_hat <- param_hat[itr,nbasis_beta+(1:nbasis_g)] 
    
    sum_g_psi_pp <- c(eval.basis(r,basisobj_g,Lfdobj=2) %*% g_hat)
    sum_g_psi_p <- c(eval.basis(r,basisobj_g,Lfdobj=1) %*% g_hat)
    
    tmp1 <- (sum_g_psi_pp * delta  - exp(eval.basis(r,basisobj_g) %*% g_hat) * sum_g_psi_p)
    
    # d2l/d(alpha1)d(alpha1)
#    H_theta[1,1] <- mean(x1*x1*tmp1)
    
    # d2l/d(alpha3)d(alpha3)
#    H_theta[2,2] <- mean(x2*x2*tmp1)
    
    # d2l/d(alpha1)d(alpha2)
#    H_theta[1,2] <- mean(x1*x2*tmp1)
#    H_theta[2,1] <- H_theta[1,2]
    
    # d2l/d(alpha1)d(beta)
    alpha_no=0
    H_theta[1,alpha_no+(1:nbasis_beta)] <- colMeans(t(int_phi_z)*x1*as.vector(tmp1))
    H_theta[alpha_no+(1:nbasis_beta),1] <- H_theta[1,alpha_no+(1:nbasis_beta)]
    
    # d2l/d(alpha2)d(beta)
    H_theta[2,alpha_no+(1:nbasis_beta)] <- colMeans(t(int_phi_z)*x2*as.vector(tmp1))
    H_theta[alpha_no+(1:nbasis_beta),2] <- H_theta[2,alpha_no+(1:nbasis_beta)]
    
    # d2l/d(beta)d(beta)
    diag(H_theta)[alpha_no+(1:nbasis_beta)] <- colMeans(t(int_phi_z)*t(int_phi_z)*as.vector(tmp1))
    
    for(i in 1:(nbasis_beta-1))
      for(j in (i+1):nbasis_beta){
        H_theta[alpha_no+i,alpha_no+j] <- mean(int_phi_z[i,] * int_phi_z[j,]*tmp1 )
        H_theta[alpha_no+j,alpha_no+i] <- H_theta[alpha_no+i,alpha_no+j]
      }
    
    tmp2 <- -eval.basis(r,basisobj_g,Lfdobj=1) * delta + as.vector(exp(eval.basis(r,basisobj_g) %*% g_hat)) * eval.basis(r,basisobj_g)
    
    # d2l/d(alpha1)d(g)
    H_theta[1,alpha_no+nbasis_beta+(1:nbasis_g)] <-  colMeans(x1*tmp2)
    H_theta[alpha_no+nbasis_beta+(1:nbasis_g),1] <- H_theta[1,alpha_no+nbasis_beta+(1:nbasis_g)]
    
    # d2l/d(alpha2)d(g)
    H_theta[2,alpha_no+nbasis_beta+(1:nbasis_g)] <-  colMeans(x2*tmp2)
    H_theta[alpha_no+nbasis_beta+(1:nbasis_g),2] <- H_theta[2,alpha_no+nbasis_beta+(1:nbasis_g)]
    
    # d2l/d(beta)d(g)
    for(i in 1:nbasis_beta){
      H_theta[alpha_no+i,alpha_no+nbasis_beta+(1:nbasis_g)] <-  colMeans(c(int_phi_z[i,])*tmp2)
      H_theta[alpha_no+nbasis_beta+(1:nbasis_g),alpha_no+i] <- H_theta[alpha_no+i,alpha_no+nbasis_beta+(1:nbasis_g)]
    }
    
    # d2l/dgdg
    # generate integrant
    gen_integ_k_l <- function(basisobj_g,g_hat,k,l){
      
      function(s){
        exp(eval.basis(s,basisobj_g) %*% g_hat) * eval.basis(s,basisobj_g)[,k] * eval.basis(s,basisobj_g)[,l]
      }
    }
    
    int_exp_phi_I <- matrix(0,nbasis_g,nbasis_g)
    
    for(i in 1:nbasis_g){
      for(j in i:nbasis_g){
        integral <- gen_integ_k_l(basisobj_g,g_hat,i,j)
        for(k in 1:n)
          int_exp_phi_I[i,j] <- int_exp_phi_I[i,j] + integrate(integral,a,r[k])$value   
        
        int_exp_phi_I[j,i] <- int_exp_phi_I[i,j]
      }
      
    }
    
    H_theta[alpha_no+nbasis_beta+(1:nbasis_g),alpha_no+nbasis_beta+(1:nbasis_g)] <- -int_exp_phi_I/n
    tryCatch(                
      
      # Specifying expression
      expr = {                      
        
        param_hat_se[itr,] <- sqrt(diag(solve(-H_theta))/n)
        
        print(param_hat_se[itr,1:2])
        
        write.csv(param_hat_se, file = paste0("../FAFT/sim_data/my_param",n,"_",error_distribution,"_", tau, "_", nbasis_beta, "_se.csv"))
        
        #########################
        # Simultaneous CI
        #########################
        # 
        ss <- seq(0,1,.01)
        # ss <- seq(0,1,.01)
        gg <- eval.basis(ss,basisobj_beta)
        var_cov <- solve(-H_theta)[alpha_no+(1:nbasis_beta),alpha_no+(1:nbasis_beta)]
        
        mean_vec <- rep(0,length(ss))
        cov_mat <- gg %*% var_cov %*% t(gg)
        cov_mat <- cov_mat/(sqrt(diag(cov_mat))%*% t(sqrt(diag(cov_mat))))
        
        Q <- mvrnorm(n = 10000,
                     mu = mean_vec,
                     Sigma = cov_mat)
        
        Q_max <- rep(0,10000)
        for(i in 1:nrow(Q)){
          Q_max[i] <- max(abs(Q[i,]))
        }
        
        Q_n_a <- as.numeric(quantile(Q_max, c(0.95))) #2.79
        
        xi <- seq(0,.5,0.05)
        for(k in 1:length(xi)){
          ss <- seq(xi[k],1-xi[k],.01)
          gg <- eval.basis(ss,basisobj_beta)
          var_cov <- solve(-H_theta)[alpha_no+(1:nbasis_beta),alpha_no+(1:nbasis_beta)]
          mean_vec <- rep(0,length(ss))
          cov_mat <- gg %*% var_cov %*% t(gg)
          cov_mat <- cov_mat/(sqrt(diag(cov_mat))%*% t(sqrt(diag(cov_mat))))
          
          beta_hat_s <- eval.basis(ss,basisobj_beta) %*% param_hat[itr,alpha_no+(1:nbasis_beta)]
          beta_cp[k,itr] <- (max(abs((beta_hat_s - beta(ss))/sqrt(diag(gg %*% var_cov %*% t(gg))/n))) <= Q_n_a)
        }    
        
        # print("Everything was fine.")
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
    
    if(itr > 1){
      # jpeg(paste0("beta",n,"_",bound,".jpg"))
      # show the pointwise average of estimated beta(s)
      ss <- seq(0,1,0.01)
      plot(ss,eval.basis(ss,basisobj_beta) %*% colMeans(param_hat[1:itr,2+(1:nbasis_beta)]),
           xlab='s',
           ylab='beta(s)',
           type="l",
           col='red')
      lines(ss,beta(ss),col='blue',lty=2)
      legend(0, -1.4, legend=c("Estimate", "True"),
             col=c("red", "blue"), lty=1:2)
      write.csv(cbind(ss,beta(ss),
                      eval.basis(ss,basisobj_beta) %*% colMeans(param_hat[1:itr,2+(1:nbasis_beta)])), 
                file = paste0("../FAFT/sim_data/curve",n,"_",error_distribution,"_", tau, "_", nbasis_beta, "_beta.csv"))
      
      # dev.off()
      
      # jpeg(paste0("g",n,"_",bound,".jpg"))
      # show the pointwise average of estimated g(s)
      ss <- seq(lower_bound,upper_bound,0.01)
      if(error_distribution=="weibull_1"){
        plot(ss,log(dweibull(exp(ss),shape=weibull_shape)*exp(ss)/pweibull(exp(ss),shape=weibull_shape, lower.tail=FALSE)), type='l', col='blue')
        write.csv(cbind(ss,log(dweibull(exp(ss),shape=weibull_shape)*exp(ss)/pweibull(exp(ss),shape=weibull_shape, lower.tail=FALSE)),
                        eval.basis(ss,basisobj_g) %*% colMeans(param_hat[1:itr,2+nbasis_beta+(1:nbasis_g)])),
                  file = paste0("../FAFT/sim_data/curve",n,"_",error_distribution,"_", tau, "_", nbasis_beta, "_g.csv"))
        
      }else if(error_distribution=="gauss_mix_0.5"){
        
        plot(ss,log(dmyMix(ss)/pmyMix(ss, lower.tail=FALSE)), type='l', col='blue')
        write.csv(cbind(ss,log(dmyMix(ss)/pmyMix(ss, lower.tail=FALSE)),
                        eval.basis(ss,basisobj_g) %*% colMeans(param_hat[1:itr,2+nbasis_beta+(1:nbasis_g)])),
                  file = paste0("../FAFT/sim_data/curve",n,"_",error_distribution,"_", tau, "_", nbasis_beta, "_g.csv"))
        
      }else if(error_distribution=="extrm"){
        
        plot(ss,log(devd(ss, scale = extrm_scale)/(1-pevd(ss, scale = extrm_scale))), type='l', col='blue')
        print(cbind(ss,log(devd(ss, scale = extrm_scale)/(1-pevd(ss, scale = extrm_scale))),
                    eval.basis(ss,basisobj_g) %*% colMeans(param_hat[1:itr,2+nbasis_beta+(1:nbasis_g)])))
        write.csv(cbind(ss,log(devd(ss, scale = extrm_scale)/(1-pevd(ss, scale = extrm_scale))),
                        eval.basis(ss,basisobj_g) %*% colMeans(param_hat[1:itr,2+nbasis_beta+(1:nbasis_g)])),
                  file = paste0("../FAFT/sim_data/curve",n,"_",error_distribution,"_", tau, "_", nbasis_beta, "_g.csv"))
        
      }
      
      lines(ss,eval.basis(ss,basisobj_g) %*% colMeans(param_hat[1:itr,2+nbasis_beta+(1:nbasis_g)]),
            xlab='s',
            ylab='g(s)',
            type='l',
            lty=2,col='red')
      legend(10, -150, legend=c("Estimate", "True"),
             col=c("red", "blue"), lty=2:1)
      # dev.off()
      
      
      
      
    }
    
    #################g C-norm###########################
    r_hat <- Y - (alpha_hat[1]*x1 + alpha_hat[2]*x2 + c(beta_hat %*% int_phi_z))
    r_hat <- pmin(pmax(r_hat,lower_bound+0.1),upper_bound-0.1)
    g_hat <- eval.basis(r_hat,basisobj_g) %*% param_hat[itr,2+nbasis_beta+(1:nbasis_g)]
    
    r <- Y - (alpha1*x1 + alpha2*x2 + z_times_beta)
    r <- pmin(pmax(r,lower_bound+0.1),upper_bound-0.1)
    
    if(error_distribution=="weibull_1"){
      
      g_true <- log(dweibull(exp(r),shape=weibull_shape)*exp(r)/pweibull(exp(r),shape=weibull_shape, lower.tail=FALSE))
      
    }else if(error_distribution=="gauss_mix_0.5"){
      
      g_true <- log(dmyMix(r)/pmyMix(r, lower.tail=FALSE))
      
    }else if(error_distribution=="extrm"){
      
      g_true <- log(devd(r, scale = extrm_scale)/(1-pevd(r, scale = extrm_scale)))
    }
    
    g_norm_A[itr] <- sqrt(mean(delta*(g_hat - g_true)^2))
    print("normA")
    print(g_norm_A[itr])
    
    ##################################################
    
    #################beta C-norm###########################
    ds <- 0.0001
    ss <- seq(0,1,ds)
    beta_hat <- eval.basis(ss,basisobj_beta) %*% param_hat[itr,2+(1:nbasis_beta)]
    beta_true <- beta(ss)
    d_beta <- beta_hat - beta_true
    beta_norm_c[itr] <- sum(d_beta*ds)^2
    
    for(i in 2:nbasis_X){
      beta_norm_c[itr] <- beta_norm_c[itr] + sum(cos((i-1)*pi*ss) * d_beta * ds)^2 * 2 / i 
    }
    beta_norm_c[itr] <- sqrt(beta_norm_c[itr]/3*(aa^2))
    beta_norm_2[itr] <- sqrt(sum(d_beta^2 * ds))
    # beta_norm_2[itr] <- sqrt(sum(abs(d_beta) * ds))
    ##################################################
    
    
    if(itr > 1){
      print('Beta CP')
      # print(paste(mean(beta_cp0[1:itr]),mean(beta_cp1[1:itr]),mean(beta_cp2[1:itr]),mean(beta_cp3[1:itr]),mean(beta_cp4[1:itr])))
      print(rowMeans(beta_cp[,1:itr],na.rm=TRUE))
      
      print('BIAS1 BIAS2')
      # print(abs(param_hat[itr,1:2] - c(alpha1,alpha2)))
      print(colMeans(abs(param_hat[1:itr,1:2] - c(alpha1,alpha2))))
      
      print('ESE1 and ESE2')
      print(colMeans(param_hat_se[1:itr,1:2], na.rm = TRUE))
      print(paste(mean(param_hat_se[1:itr,1], na.rm = TRUE),
                  mean(param_hat_se[1:itr,2], na.rm = TRUE)))
      
      print('current CP1 and CP2')
      print(paste(mean(abs(param_hat[1:itr,1] - alpha1) > 1.96*param_hat_se[1:itr,1], na.rm = TRUE),
                  mean(abs(param_hat[1:itr,2] - alpha2) > 1.96*param_hat_se[1:itr,2], na.rm = TRUE)))
      
      print('current norm_c and norm_2')
      print(paste(beta_norm_c[itr],beta_norm_2[itr],g_norm_c[itr],g_norm_A[itr]))
      print(paste(mean(beta_norm_c[1:itr]),mean(beta_norm_2[1:itr]),mean(g_norm_c[1:itr],na.rm=TRUE),mean(g_norm_A[1:itr],na.rm=TRUE)))
      # print(paste(mean(beta_norm_c[1:itr]^2),mean(beta_norm_2[1:itr]^2),mean(g_norm_c[1:itr]^2),mean(g_norm_2[1:itr]^2)))
      
      write.csv(c(colMeans(abs(param_hat[1:itr,1:2] - c(alpha1,alpha2))),
                  colMeans(param_hat_se[1:itr,1:2], na.rm = TRUE),
                  mean(abs(param_hat[1:itr,1] - alpha1) > 1.96*param_hat_se[1:itr,1], na.rm = TRUE),
                  mean(abs(param_hat[1:itr,2] - alpha2) > 1.96*param_hat_se[1:itr,2], na.rm = TRUE),
                  mean(beta_norm_c[1:itr]),mean(beta_norm_2[1:itr]),mean(g_norm_c[1:itr],na.rm=TRUE),mean(g_norm_A[1:itr],na.rm=TRUE),
                  rowMeans(beta_cp[,1:itr],na.rm=TRUE)), 
                file = paste0("../FAFT/sim_data/curve",n,"_",error_distribution,"_", tau, "_", nbasis_beta, "_alpha.csv"))
    }
  }
  
}
