
# Calculating Abar for given gs(w), which is in the form of gs(w) = p1*w^p2 + p3*w^p4 where p1 to p4 are given
f2 <- function(par,
               ca, k, MAP,
               Vcmax=50, cp=30, Km=703, Rd=1, VPD=0.02, p=43200, LAI=1, nZ=0.5,
               a=1.6, l=1.8e-5, MDP=MAP/365, gamma=1/((MDP/k)/1000)*nZ, h=l*a*LAI/nZ*p){
  
  p1 <- abs(par[1]) # we require p1 be nonnegativ
  p2 <- par[2]
  p3 <- abs(par[3]) # we require p3 be nonnegativ
  p4 <- par[4]
  
  g <- function(w){p1*w^p2+p3*w^p4}
  A <- function(w){LAI*1/2*(Vcmax+(Km+ca)*g(w)-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*g(w)+((ca+Km)*g(w)+Rd)^2-2*Rd*Vcmax)^(1/2))}
  Ev <- function(w){h*VPD*g(w)}
  rEv <- function(w){1/Ev(w)}
  integralEv <- Vectorize(function(w){integrate(rEv, 0, w)$value})
  fnoc <- function(w){1/Ev(w)*exp(-gamma*w+k*integralEv(w))}
  integralfnoc <- integrate(fnoc, 0, 1)$value
  c <- 1/(integralfnoc+1/k)
  f <- function(w){c*fnoc(w)}
  fA <- function(w){f(w)*A(w)}
  averA <- integrate(fA, 0, 1)$value
  return(averA)
}

f2_wrapper <- function(par){
  ca <- ca
  k <- k
  MAP <- MAP
  res <- f2(par, ca, k, MAP)
  return(res)
}
