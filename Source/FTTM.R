fit_FTTM <- function(data,
                     N0grid = c(4, 7, 10, 13), # N0: order of polynomial basis to estimate H(t)
                     N1grid = c(3, 5, 7, 9),   # N1: order of Bernstein basis to estimate beta(s)
                     rgrid = seq(0, 4, length = 10) # r: logarithmic error distribution class
                     ){
  
  # T_obs = data$Y # observed survival time
  # delta = data$delta # censoring indicator
  # n = length(T_obs)
  # Xmat = matrix(1, nrow = n, ncol = 1) # intercept
  # Xs = as.matrix(data$X) # functional covariate matrix
  
  # tau <- ceiling(max(T_obs))
  # p1 <- ncol(Xmat)
  
  # S <- seq(0, 1, length = ncol(Xs))

  #r=0 #Cox specification
  gridN <- expand.grid(N0 = N0grid, N1 = N1grid, r = rgrid)
  gridN <- as.matrix(gridN)
  AICres <- apply(gridN, 1, ChooseAICval)  
  AICvec <- c()
  for (k in 1:nrow(gridN)){
    AICvec[k]<-AICres[[k]]$AICval
  }
  indmin <- which.min(AICvec)
  Nrvecsel <- gridN[indmin,] # optimal N with minimum AIC
  N0 <- as.numeric(Nrvecsel[1])
  N1 <- as.numeric(Nrvecsel[2])
  r <- as.numeric(Nrvecsel[3])
  
  finalest <- estFTTM.final(Nrvecsel)
  
  psiest <- -finalest$psiest
  SEpar <- finalest$SEpar
  VCOVpar <- finalest$VCOVpar
  
  betaest <- psiest[1:p1] # estimated beta0
  #true beta#
  #simdata$coefficients$beta0[1]
  ##############
  
  etaest <- psiest[(p1+1):(p1+N0+1)]
  gammaest <- c()
  gammaest[1] <- etaest[1]
  for (k in 2: length(etaest)){
    gammaest[k]<-gammaest[k-1]+exp(etaest[k])
  }
  
  thetaest <- psiest[(p1+N0+1+1):(p1+N0+1+N1+1)]
  betas_tt_est <- beta_s_estfun(S, thetaest, N1) # estimated beta1(s)
  
  ###covariance matrix of parameters##
  VCOVparfunc <- VCOVpar[(p1+N0+1+1):(p1+N0+1+N1+1),(p1+N0+1+1):(p1+N0+1+N1+1)]
  basismat <- NULL
  for(k in 0:N1){
    basismat <- cbind(basismat, dbeta((S), shape1 = (k+1), shape2 = N1-k+1))
  }
  basismat <- (1/(N1+1)) * basismat
  vcov_betas <- basismat%*%VCOVparfunc%*%t(basismat)
  SE_betas <- sqrt(diag(vcov_betas))
  
  #asymptotic 95% point-wise confidence interval
  lcb_betas <- betas_tt_est-1.96*SE_betas
  ucb_betas <- betas_tt_est+1.96*SE_betas
  
  #linear predictor
  lp = Xmat %*% betaest + Xs %*% betas_tt_est / ncol(Xs)
  
  return(list(beta0_hat = betaest, beta1_hat = betas_tt_est, beta1_ci_lower = lcb_betas, beta1_ci_upper = ucb_betas, 
              lp = lp, gamma_hat = gammaest, N0 = N0, N1 = N1, r = r, tau = tau))
  
}

predict_FTTM <- function(fit, newdata){
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



# # estimated H(t) vs true H(t)
# plot(tgrid, H0, type = "l", col = "red", lwd = 1.5,
#      xlab = "t", ylab = "Value",
#      main = expression(paste("Comparison of ", log(t), " and ", hat(H)(t))))
# lines(tgrid, log(tgrid), col = "blue", lwd = 1.5)
# legend("topright", legend = c(expression(log(t)), expression(hat(H)(t))),
#        col = c("blue", "red"), lty = 1, lwd = 1.5)
# 
# 
# plot(S,betaS,col=1,type="l",ylim=c(-2,3), main = "Simple true beta(s)") 
# lines(S,betas_tt_est,col=4)
# lines(S,lcb_betas,col=4,lty=3)
# lines(S,ucb_betas,col=4,lty=3) #CI looks good
# abline(h=0)

G<-function(t,r=0)
{if(r==0)
{
  val<-t  
}
  if(r>0)
  {val<-log(1+r*t)/r}
  val
}

F_epsilon_bar<-function(t,r=0)
{exp(-G(exp(t),r=r))
}



f_epsilon<-function(t,r=0)
{
  if(r==0)
  {
    val<-(exp(t))*(exp(-exp(t)))  
  }
  
  if(r>0)
  {val<- (exp(-G(exp(t),r)))*(exp(t))*(1/(1+r*exp(t)))}
  val
}

#need to define AIC,BIC, to choose N0,N1

H_baseline<-function(t,gamma,N0,Tau)
{
  basismat <- NULL
  for(k in 0:N0){
    basismat<- cbind(basismat, dbeta((t/Tau), shape1=(k+1), shape2=N0-k+1))
  }
  #basis should add up to 1
  basismat<-(1/(N0+1))*basismat
  val<-as.numeric(basismat%*%gamma)
  val
}

H_baseline_der<-function(t,gamma,N0,Tau)
{
  dgamma<-diff(gamma)
  nbas<-N0-1
  basismat <- NULL
  for(k in 0:nbas){
    basismat<- cbind(basismat, dbeta((t/Tau), shape1=(k+1), shape2=nbas-k+1))
  }
  #basis should add up to 1
  basismat<-(1/(nbas+1))*basismat
  val<-as.numeric(basismat%*%dgamma)
  val*N0
}


beta_s_estfun<-function(s,theta,N1)
{
  basismat <- NULL
  for(k in 0:N1){
    basismat<- cbind(basismat, dbeta((s), shape1=(k+1), shape2=N1-k+1))
  }
  #basis should add up to 1
  basismat<-(1/(N1+1))*basismat
  val<-as.numeric(basismat%*%theta)
  val
}



S_N_fun<-function(t,beta,gamma,theta,Xvec,Xf,S,tau,r=0,N0,N1)
{val1<-H_baseline(t,gamma,N0,tau)
val2<-as.numeric(Xvec%*%beta)
basismat <- NULL
for(k in 0:N1){
  basismat<- cbind(basismat, dbeta((S), shape1=(k+1), shape2=N1-k+1))
}
#basis should add up to 1
basismat<-(1/(N1+1))*basismat
val3w<-(Xf%*%basismat)*(1/length(S))
val3<-as.numeric(val3w%*%theta)
totval<-val1+val2+val3
F_epsilon_bar(totval,r=r)
}


f_N_fun<-function(t,beta,gamma,theta,Xvec,Xf,S,tau,r=0,N0,N1)
{val1<-H_baseline(t,gamma,N0,tau)
val2<-as.numeric(Xvec%*%beta)
basismat <- NULL
for(k in 0:N1){
  basismat<- cbind(basismat, dbeta((S), shape1=(k+1), shape2=N1-k+1))
}
#basis should add up to 1
basismat<-(1/(N1+1))*basismat
val3w<-(Xf%*%basismat)*(1/length(S))
val3<-as.numeric(val3w%*%theta)
totval<-val1+val2+val3
f_epsilon(totval,r=r)*H_baseline_der(t,gamma,N0,tau)*(1/tau)
}

loglikfun2<-function(psi,T_obs,delta,Xmat,Xs,S,tau,r,N0,N1)
{ 
  p1<-ncol(Xmat)
  beta<-psi[1:p1]
  #gamma<-psi[(p1+1):(p1+N0+1)] 
  ##define gamma from eta
  eta<-psi[(p1+1):(p1+N0+1)]
  gamma<-c()
  gamma[1]<-eta[1]
  for (k in 2: length(eta))
  {gamma[k]<-gamma[k-1]+exp(eta[k])
  }
  
  theta<-psi[(p1+N0+1+1):(p1+N0+1+N1+1)]
  #ll<-c()
  temp1<-delta*log(f_N_fun(T_obs,beta,gamma,theta,Xvec=Xmat,Xf=Xs,S,tau,r=r,N0,N1))
  temp2<-(1-delta)*log(S_N_fun(T_obs,beta,gamma,theta,Xvec=Xmat,Xf=Xs,S,tau,r=r,N0,N1))
  ll<-temp1+temp2
  loglikval<-sum(ll) 
  objval<--loglikval
  if(is.finite(objval))
  {return(objval)}
  if(!is.finite(objval))
  {return(9999999)}
  
}

#have to start at good value

choose.init<-function(T_obs,delta,Xmat,Xs,S,N1)
{basismat <- NULL
for(k in 0:N1){
  basismat<- cbind(basismat, dbeta((S), shape1=(k+1), shape2=N1-k+1))
}
#basis should add up to 1
basismat<-(1/(N1+1))*basismat
val3w<-(Xs%*%basismat)*(1/length(S)) #these are w
library(mgcv)
fit<-gam(T_obs~Xmat+val3w,family=cox.ph(),weights=delta)
p1<-ncol(Xmat)
betainit<-as.numeric(coef(fit)[1:p1])
thetainit<-as.numeric(coef(fit)[(p1+1):(p1+1+N1)])
result<-list(betainit=betainit,thetainit=thetainit)
return(result)
}


ICfunc<-function(psiest,N0,N1,r)
{neglval<-2*loglikfun2(psiest,T_obs,delta,Xmat,Xs,S,tau,r=r,N0,N1)
AIC<-neglval+2*(p1+N0+N1+1)
return(AIC)
}

ChooseAICval<-function(Nrvec)
{
  N0<-as.numeric(Nrvec[1])
  N1<-as.numeric(Nrvec[2])
  r<-as.numeric(Nrvec[3])
  library(mgcv)
  initval<-choose.init(T_obs=T_obs,delta=delta,Xmat=Xmat,Xs=Xs,S=S,N1=N1)#N0 not needed here
  print(initval)
  betainit<-initval$betainit
  thetainit<-initval$thetainit
  gammainit<- sort((rnorm((N0+1),0,1)))#c(-0.5,0.1,0.3,0.5) #has to be increasing #H(t)=beta0+logt
  #need to give better values by understanding H#
  etainit<-c()
  etainit[1]<-gammainit[1]
  for(k in 2: length(gammainit))
  {etainit[k]<-log(gammainit[k]-gammainit[k-1]) #these cannot be all same.. chose baseline hazard differently
  }
  psi<-c(betainit,etainit,thetainit)
  out <- optim(psi,loglikfun2,T_obs=T_obs,delta=delta,Xmat=Xmat,Xs=Xs,S=S,tau=tau,r=r,N0=N0,N1=N1, method = "BFGS",control=list(trace=1,maxit=500)) #inc maxit
  psiest<-out$par
  AICval<-ICfunc(psiest,N0,N1,r)
  #BICval  
  result<-list(AICval=AICval,psiest=psiest)
  #print(Nrvec)
  return(result)
}


estFTTM.final<-function(Nrvec)
{ N0<-as.numeric(Nrvec[1])
N1<-as.numeric(Nrvec[2])
r<-as.numeric(Nrvec[3])
library(mgcv)
initval<-choose.init(T_obs,delta,Xmat,Xs,S,N1=N1)#N0 not needed here
print(initval)
betainit<-initval$betainit
thetainit<-initval$thetainit
gammainit<- sort((rnorm((N0+1),0,1)))#c(-0.5,0.1,0.3,0.5) #has to be increasing #H(t)=beta0+logt
etainit<-c()
etainit[1]<-gammainit[1]
for(k in 2: length(gammainit))
{etainit[k]<-log(gammainit[k]-gammainit[k-1]) 
}
psi<-c(betainit,etainit,thetainit)
out <- optim(psi,loglikfun2,T_obs=T_obs,delta=delta,Xmat=Xmat,Xs=Xs,S=S,tau=tau,r=r,N0=N0,N1=N1, method = "BFGS",control=list(trace=1,maxit=500),hessian = TRUE) #inc maxit
psiest<-out$par
H<-out$hessian   
library(Matrix)
nH<-nearPD(H)$mat
VCOVpar<-solve(nH)
SEpar<-sqrt(diag(VCOVpar)) 
AICval<-ICfunc(psiest,N0,N1,r)
result<-list(AICval=AICval,psiest=psiest,SEpar=SEpar,VCOVpar=VCOVpar)
#print(Nrvec)
return(result)
}