
source("R/Single power.r")
source("R/Double power.r")

# gw: the optimal gs(w), following Eq. 14
gw <- function(averA,
               ca, k, MAP,
               Vcmax=50, cp=30, Km=703, Rd=1, VPD=0.02, p=43200, LAI=1, nZ=0.5,
               a=1.6, l=1.8e-5, MDP=MAP/365, gamma=1/((MDP/k)/1000)*nZ, h=l*a*LAI/nZ*p){
  
  f1 <- function(w){
    w <- w
    O <- function(g){
      integrand <- function(g){(4*h*(cp+Km)*LAI*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*g)/(((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))*(ca*k*LAI*Rd+k*Km*LAI*Rd-ca*k*LAI*Vcmax+2*cp*k*LAI*Vcmax+k*Km*LAI*Vcmax+gamma*h*LAI*Rd^2*VPD-2*gamma*h*LAI*Rd*Vcmax*VPD+gamma*h*LAI*Vcmax^2*VPD+ca^2*k*LAI*g+2*ca*k*Km*LAI*g+k*Km^2*LAI*g+ca*gamma*h*LAI*Rd*VPD*g+gamma*h*Km*LAI*Rd*VPD*g-ca*gamma*h*LAI*Vcmax*VPD*g+2*cp*gamma*h*LAI*Vcmax*VPD*g+gamma*h*Km*LAI*Vcmax*VPD*g-ca*k*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))-k*Km*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+2*gamma*h*averA*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+gamma*h*LAI*Rd*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))-gamma*h*LAI*Vcmax*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))))}
      integrate(integrand, 0, g)$value-w
    }
    u <- uniroot(O, c(0,1), tol=1e-10, extendInt = "yes")
    return(u$root)
  }
  
  f2 <- Vectorize(f1)
  
  g <- vector("numeric", length=length(w))
  g <- f2(w=w)

  return(g)
}

ca <- 400
k <- 0.025
MAP <- 1095

averA <- dvs[dvs$ca==ca & dvs$k==k & dvs$MAP==MAP, "A"]
w <- seq(0, 1, length=1e3)
g <- gw(averA, ca, k, MAP)

#single power
int1 <- c(0.5, 0.5)
# gs(w) as a single power function with optimized parameters
par1 <- optimx(int1, f1_wrapper, itnmax=5000, method="BFGS", control=list(maximize=T))
gw1 <- function(w){abs(par1$p1)*w^par1$p2}

#double power
int2 <- c(0.1, 4, 0.1, 0.5)
# gs(w) as a double power function with optimized parameters
par2 <- optimx(int2, f2_wrapper, itnmax=5000, method="BFGS", control=list(maximize=T))
gw2 <- function(w){abs(par2$p1)*w^par2$p2+abs(par2$p3)*w^par2$p4}

# Plots
windows(8, 6)
par(mgp = c(2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(3.5, 4, 1, 2))

# Optimal
plot(w*100, g, type = "l",
     xlab = expression(italic(w)~"(%)"), 
     ylab = expression(italic(g[s])~(mol~m^-2~s^-1)),
     xlim = c(0, 100), ylim = c(0, 0.2),
     cex.lab = 1.3,
     col = "red"
)
#single power
data1 <- gw1(w)
lines(w*100, data1, lty = 1, col = "darkgreen")
#double power
data2 <- gw2(w)
lines(w*100, data2, lty = 2, col = "blue")

Legend = c("Optimal function", "Single-term power function", "Double-term power function")
legend("topleft", Legend, lty = c(1, 1, 2), col = c("red", "darkgreen", "blue"))

box()

dev.copy2pdf(file = "output/figures/Figure 3.pdf")
