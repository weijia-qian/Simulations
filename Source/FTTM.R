# simdata <- simulate_AFT(family = "lognormal", n = 500, nS = 500, beta_type = "simple, beta_0 = 0.5,
#                         b = 0.5)
simdata <- readRDS("simdata.rds")
source("sourceFTTM.R")

T_obs=simdata$data$Y # observed survival time
delta=simdata$data$delta # censoring indicator
n=length(T_obs)
Xmat=matrix(1, nrow = n, ncol = 1) # intercept
Xs=simdata$data$X # functional covariate matrix
betaS=simdata$coefficients$beta1*-1 # true beta(s)

tau<-ceiling(max(T_obs))
p1<-ncol(Xmat)

S <- seq(0, 1, length = ncol(Xs))

N0grid<-c(4, 7, 10, 13)
N1grid<-c(3, 5, 7, 9) #can make big given resources
#r=0 #Cox specification
rgrid<-seq(0, 4, length = 5)
gridN<-expand.grid(N0=N0grid,N1=N1grid,r=rgrid)
gridN<-as.matrix(gridN)
AICres<-apply(gridN,1,ChooseAICval)  
AICvec<-c()
for (k in 1:nrow(gridN))
{AICvec[k]<-AICres[[k]]$AICval
}
indmin<-which.min(AICvec)
Nrvecsel<-gridN[indmin,] #optimal N with minimum AIC
N0<-as.numeric(Nrvecsel[1])
N1<-as.numeric(Nrvecsel[2])
r<-as.numeric(Nrvecsel[3])

finalest<-estFTTM.final(Nrvecsel)

psiest<-finalest$psiest
SEpar<-finalest$SEpar
VCOVpar<-finalest$VCOVpar

betaest<-psiest[1:p1] 
### estimated beta ###
betaest
#true beta#
simdata$coefficients$beta0[1]*-1
##############

etaest<-psiest[(p1+1):(p1+N0+1)]
gammaest<-c()
gammaest[1]<-etaest[1]
for (k in 2: length(etaest))
{gammaest[k]<-gammaest[k-1]+exp(etaest[k])
}
tgrid <- seq(0, 120, length = 200)
H_est <- H_baseline(t = tgrid, gamma = gammaest, N0, Tau = tau)

# estimated H(t) vs true H(t)
plot(tgrid, H_est, type = "l", col = "red", lwd = 1.5,
     xlab = "t", ylab = "Value",
     main = expression(paste("Comparison of ", log(t), " and ", hat(H)(t))))
lines(tgrid, log(tgrid), col = "blue", lwd = 1.5)
legend("topright", legend = c(expression(log(t)), expression(hat(H)(t))),
       col = c("blue", "red"), lty = 1, lwd = 1.5)

thetaest<-psiest[(p1+N0+1+1):(p1+N0+1+N1+1)]
betas_tt_est<-beta_s_estfun(S,thetaest,N1)
###covariance matrix of parameters##
VCOVparfunc<-VCOVpar[(p1+N0+1+1):(p1+N0+1+N1+1),(p1+N0+1+1):(p1+N0+1+N1+1)]
basismat <- NULL
for(k in 0:N1){
  basismat<- cbind(basismat, dbeta((S), shape1=(k+1), shape2=N1-k+1))
}
basismat<-(1/(N1+1))*basismat
vcov_betas<-basismat%*%VCOVparfunc%*%t(basismat)
SE_betas<-sqrt(diag(vcov_betas))

#asymptotic 95% point-wise confidence interval
lcb_betas<-betas_tt_est-1.96*SE_betas
ucb_betas<-betas_tt_est+1.96*SE_betas

plot(S,betaS,col=1,type="l",ylim=c(-2,3), main = "Simple true beta(s)") 
lines(S,betas_tt_est,col=4)
lines(S,lcb_betas,col=4,lty=3)
lines(S,ucb_betas,col=4,lty=3) #CI looks good
abline(h=0)

