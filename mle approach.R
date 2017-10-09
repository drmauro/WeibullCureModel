######## Weibull cure rate regression model 
time = data$time
delta = data$delta
n=length(time)
loglike_weibull = function(varp)
{
  b1=varp[1]
  b2=varp[2]
  b3=varp[3]
  b4=varp[4]
  b5=varp[5]
  b6=varp[6]
  
  ## Weibull parameters linked with one covariates y
  a   = exp(varp[1])*exp(y*varp[2])
  lam = exp(varp[3])*exp(y*varp[4])
  
  ## cure fraction linked with one covariates y
  pz = (exp(varp[5])*exp(y*varp[6]))/(1+exp(varp[5])*exp(y*varp[6]))time 
  
  fpop = (1-pz)*dweibull(time, shape=a, scale = lam, log = FALSE)
  Spop = pz+ (1-pz)*pweibull(time,a,lam,lower.tail = F , log = FALSE ) 
  loglike = sum(delta*log(fpop) +(1-delta)*log(Spop))
  loglike
}

#Test different starting points.
ini = c(-1,-1,-1,-1,-1,-1)
loglike_weibull(ini) 
res_weibull = optim(ini, loglike_weibull, method="BFGS", control= list(fnscale=-1),hessian=F )
res_weibull