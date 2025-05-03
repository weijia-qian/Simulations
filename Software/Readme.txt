1. FTTM.html show illustration of the proposed FTTM
   for simulation scenarios A1,A2 in the paper.
 
2. FTTM.Rmd file contains the R markdown code for the same.

3. Source the sourcFTTM.R file which loads the following functions used for estimation.

i) G function, where the survival function of error epsilon is given by exp(-G(e^t)).
#@t= time
#@r= choice of r in the logarithmic error distribution class
Use: G(t,r=r)
Value: numeric value of G(t,r=r)

ii) F_epsilon_bar function, the survival function of error epsilon.
#@t= time
#@r= choice of r in the logarithmic error distribution class
Use: F_epsilon_bar(t,r=r)
Value: numeric value of F_epsilon_bar(t,r)


iii) f_epsilon function, the density of error epsilon.
#@t= time
#@r= choice of r in the logarithmic error distribution class
Use: f_epsilon(t,r=r)
Value: numeric value of f_epsilon(t,r)

iv) H_baseline: The trasnformation function H(.).
#@t= time
#@gamma= basis coefficients
#@N0= order of Bernstein polynomial basis functions for H
#@Tau= Defined as in paper (maximum observed survival time)
Use: H_baseline(t,gamma,N0,Tau)
Value: numeric value of H function at t/Tau.

v) H_baseline_der: Derivative of the trasnformation function H(.).
#@t= time
#@gamma= basis coefficients
#@N0= order of Bernstein polynomial basis functions
#@Tau= Defined as in paper (maximum observed survival time)
Use: H_baseline_der(t,gamma,N0,Tau)
Value: numeric value of derivative H function at t/Tau

vi) beta_s_estfun: functional coefficient beta(s).
#@s= functional argument
#@theta= basis coefficients
#@N1= order of Bernstein polynomial basis functions
Use: beta_s_estfun(s,theta,N1)
Value: numeric value of functional coefficient beta(s).

vii) choose.init: choosing initial values of the model parameters.
#@T_obs= observed event times
#@delta= event indicators
#@Xmat= scalar covariate matrix
#@Xs= functional covariate matrix (n*m) observed on a grid
#@S= grid on which functional covariate is observed
#@N1= order of Bernstein polynomial basis functions
Use: choose.init(T_obs,delta,Xmat,Xs,S,N1)
Value:
$betainit: initial values for scalar parameters
$thetainit:initial values for basis coeffcients corresponding to slope beta(s)


viii) loglikfun2: observed log-likelihood (negative).
#@psi= basis coefficients
#@T_obs= observed event times
#@delta= event indicators
#@Xmat= scalar covariate matrix
#@Xs= functional covariate matrix (n*m) observed on a grid
#@S= grid on which functional covariate is observed
#@Tau= Defined as in paper (maximum observed survival time)
#@r= choice of r in the logarithmic error distribution class
#@N0= order of Bernstein polynomial basis functions for H
#@N1= order of Bernstein polynomial basis functions
Use: loglikfun2(psi,T_obs,delta,Xmat,Xs,S,tau,r,N0,N1)
Value: value of the negative log-likelihood

ix) ICfunc: function for calculating AIC.
#@psiest= basis coefficients
#@N0= order of Bernstein polynomial basis functions for H
#@N1= order of Bernstein polynomial basis functions
#@r= choice of r in the logarithmic error distribution class
Use: ICfunc(psiest,N0,N1,r)
Value: AIC value


x) ChooseAICval: Funtion for estimation, returns AIC value and estimated basis coefficients.
#@Nrvec= A vector with following three elements (in this order)
N0= order of Bernstein polynomial basis functions for H
N1= order of Bernstein polynomial basis functions
r= choice of r in the logarithmic error distribution class
Use: ChooseAICval(Nrvec)
Value: 
$AICval= Value of AIC at the MLE
$psiest= MLE of the basis coefficients

xi) estFTTM.final: Funtion for estimation and inference, returns AIC value, estimated basis coefficients and their estimated variance covariance matrix.
#@Nrvec= A vector with following three elements (in this order)
N0= order of Bernstein polynomial basis functions for H
N1= order of Bernstein polynomial basis functions
r= choice of r in the logarithmic error distribution class
Use: estFTTM.final(Nrvec)
Value: 
$AICval= Value of AIC at the MLE
$psiest= MLE of the basis coefficients
$SEpar= Standard error of the estimated basis coefficients
$VCOVpar= Variance covariance matrix of the estimated basis coefficients

















