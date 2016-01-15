
#Derived variables
dvsf <- function(ca, k, MAP,
                q,
                Vcmax=50, cp=30, Km=703, Rd=1, VPD=0.02, p=43200, LAI=1, nZ=0.5,
                a=1.6, l=1.8e-5, MDP=MAP/365, gamma=1/((MDP/k)/1000)*nZ, h=l*a*LAI/nZ*p){
  
  f1 <- function(w){
    w <- w
    O <- function(g){
      integrand <- function(g){(4*h*(cp+Km)*LAI*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*g)/(((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))*(ca*k*LAI*Rd+k*Km*LAI*Rd-ca*k*LAI*Vcmax+2*cp*k*LAI*Vcmax+k*Km*LAI*Vcmax+gamma*h*LAI*Rd^2*VPD-2*gamma*h*LAI*Rd*Vcmax*VPD+gamma*h*LAI*Vcmax^2*VPD+ca^2*k*LAI*g+2*ca*k*Km*LAI*g+k*Km^2*LAI*g+ca*gamma*h*LAI*Rd*VPD*g+gamma*h*Km*LAI*Rd*VPD*g-ca*gamma*h*LAI*Vcmax*VPD*g+2*cp*gamma*h*LAI*Vcmax*VPD*g+gamma*h*Km*LAI*Vcmax*VPD*g-ca*k*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))-k*Km*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+2*gamma*h*q*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+gamma*h*LAI*Rd*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))-gamma*h*LAI*Vcmax*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))))}
      integrate(integrand, 0, g)$value-w
    }
    u <- uniroot(O, c(0,1), tol=1e-10, extendInt = "yes")
    return(u$root)
  }
  
  f2 <- Vectorize(f1)
  
  A <- function(w){LAI*1/2*(Vcmax+(Km+ca)*f2(w)-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*f2(w)+((ca+Km)*f2(w)+Rd)^2-2*Rd*Vcmax)^(1/2))}
  Ev <- function(w){h*VPD*f2(w)}
  rEv <- function(w){1/Ev(w)}
  integralEv <- Vectorize(function(w){integrate(rEv, 0, w)$value})
  fnoc <- function(w){1/Ev(w)*exp(-gamma*w+k*integralEv(w))}
  integralfnoc <- integrate(fnoc, 0, 1)$value
  c <- 1/(integralfnoc+1/k)
  f <- function(w){c*fnoc(w)}
  
  #derived variables
  frw <- function(w){f(w)*w*100}
  averrw <- integrate(frw, 0, 1)$value
  fEv <- function(w){f(w)*Ev(w)}
  averEv <- integrate(fEv, 0, 1)$value
  fg1 <- function(w){sqrt(VPD*100)*(f2(w)*ca/A(w)-1)*f(w)}
  averg1 <- integrate(fg1, 0, 1)$value
  fcica <- function(w){(ca-A(w)/f2(w))*f(w)/ca}
  avercica <- integrate(fcica, 0, 1)$value
  p0 <- c/k*100
  lambda <- function(w){(2*h*VPD)/(-((2*Vcmax*(-ca+2*cp+Km)+2*(ca+Km)*(f2(w)*(ca+Km)+Rd))/(2*sqrt(2*f2(w)*Vcmax*(-ca+2*cp+Km)+(f2(w)*(ca+Km)+Rd)^2-2*Rd*Vcmax+Vcmax^2)))+ca+Km)}
  averlambda <- integrate(lambda, 0, 1)$value
  
  res <- c(averrw, averEv, averg1, avercica, p0, averlambda, averEv*nZ*1000/MDP)
  return(res)
}
