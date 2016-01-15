
# a simulation model: generating data of daily soil water content,
# stomatal conductance, and photosynthesis rate under stochastic rainfall
f <- function(averA,
              ca, k, MAP,
              Vcmax=50, cp=30, Km=703, Rd=1, VPD=0.02, p=43200, LAI=1, nZ=0.5,
              a=1.6, l=1.8e-5, MDP=MAP/365, gamma=1/((MDP/k)/1000)*nZ, h=l*a*LAI/nZ*p){
  
  n <- 1e5 #the number of rainfall events
  RW0 <- 0.01 #initial soil water content
  
  set.seed(1)
  
  #fAD: generating AD (a time series of rainfall)
  fAD <- function(n){
    # LD[i] - the length of the interval between (i-1)th and ith rainfall events
    LD <- rexp(n, k)
    # DR[i] - the day when the ith rainfall event occurs
    DR <- floor(cumsum(LD))
    # LS - the length of simulation run
    LS <- tail(DR, 1)
    #AR[i] <- the amount of ith rainfall event
    AR <- rexp(n, gamma)
    #AD[i] - the amount of rainfall on the ith day
    AD <- rep(0, length=LS)
    for (i in 1:n){
      AD[DR[i]] <- AD[DR[i]]+AR[i]
    }
    return(AD)
  }
  
  AD <- fAD(n)
  LS <- length(AD)
  
  gsw1 <- function(w){
    w <- w
    O <- function(g){
      integrand <- function(g){(4*h*(cp+Km)*LAI*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*g)/(((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))*(ca*k*LAI*Rd+k*Km*LAI*Rd-ca*k*LAI*Vcmax+2*cp*k*LAI*Vcmax+k*Km*LAI*Vcmax+gamma*h*LAI*Rd^2*VPD-2*gamma*h*LAI*Rd*Vcmax*VPD+gamma*h*LAI*Vcmax^2*VPD+ca^2*k*LAI*g+2*ca*k*Km*LAI*g+k*Km^2*LAI*g+ca*gamma*h*LAI*Rd*VPD*g+gamma*h*Km*LAI*Rd*VPD*g-ca*gamma*h*LAI*Vcmax*VPD*g+2*cp*gamma*h*LAI*Vcmax*VPD*g+gamma*h*Km*LAI*Vcmax*VPD*g-ca*k*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))-k*Km*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+2*gamma*h*averA*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+gamma*h*LAI*Rd*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))-gamma*h*LAI*Vcmax*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))))}
      integrate(integrand, 0, g)$value-w
    }
    u <- uniroot(O, c(0,1), tol=1e-10, extendInt = "yes")
    return(u$root)
  }
  
  gsw <- Vectorize(gsw1)
  
  #g[i] - average stomatal conductance on the ith day
  g <- vector("numeric")
  g[1] <- 0
  # A[i] - photosynthesis rate on the ith day; μmol per second
  A <- vector("numeric")
  A[1] <- 0
  # E[i] - the amount of soil water transpired on ith day  
  E <- vector("numeric")
  E[1] <- 0
  # I[i] - the amount (on relative scale) of rain water infiltrating the soil on ith day
  I <- vector("numeric")
  I[1] <- 0
  # dRWdt[i] - the change of relative soil water content on ith day
  dRWdt <- vector("numeric")
  dRWdt[1] <- 0
  # RW[i] - relative soil water content on ith day
  RW <- vector("numeric")
  RW[1] <- RW0
  
  # Simulation, stops after 1000 days
  for(i in 2:LS){
    g[i] <- gsw(RW[i-1])
    A[i] <- (LAI*1/2*(Vcmax+(Km+ca)*g[i]-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*g[i]+((ca+Km)*g[i]+Rd)^2-2*Rd*Vcmax)^(1/2)))
    E[i] <- h*VPD*g[i]
    #sum of infiltration and current RW must be smaller than 1
    I[i] <- min(AD[i], 1-RW[i-1])
    dRWdt[i] <- I[i]-E[i]
    x <- RW[i-1]+dRWdt[i]
    RW[i] <- ifelse(x<0, 0, x)
    if (i>=1000) break
  }
  res <- data.frame(w=RW, g=g, A=A)
  return(res)
}

ca <- 400
k <- 0.025
MAP <- 1095

averA <- dvs[dvs$ca==ca & dvs$k==k & dvs$MAP==MAP, "A"]
res <- f(averA, ca, k, MAP)
write.csv(res, "output/data/Figure 1.csv", row.names = FALSE)

#plots
windows(8, 6)
Cols <- c("red", "darkgreen", "blue")
par(mar=c(3, 6.3, 1, 2.8), lwd = 2)
plot(res$w[2:1000]*100, type="s",
     xlab = "", ylab = "",
     xlim = c(0, 1000), ylim = c(0, 100),
     col = Cols[3],
     axes = FALSE
)

axis(2, ylim = c(0, 100), pos = 0, lwd = 2)
mtext(expression(italic(w~"(%)")),side = 2,line = 0.9, cex = 1.3)
text(x = 500, y= 19, labels = expression(italic(w)), col = Cols[3], cex = 1.3)

par(new=TRUE)

plot(res$A[2:1000], type="s",
     xlab = "", ylab = "",
     xlim = c(0, 1000), ylim = c(0, 15),
     col = Cols[1],
     axes = FALSE
)

axis(2, ylim = c(0, 10), pos = -130, lwd = 2)
mtext(expression(italic(A)~(mu*mol~m^-2~s^-1)),side = 2,line = 4.5, cex = 1.3)
text(x = 500, y= 11.7, labels = expression(italic(A)), col = Cols[1], cex = 1.3)

par(new=TRUE)

plot(res$g[2:1000], type="s",
     xlab = "", ylab = "",
     xlim = c(0, 1000), ylim = c(0, 0.15),
     col = Cols[2],
     axes = FALSE
)

axis(4, ylim = c(0, 0.03), pos = 1000, lwd = 2)
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)),side = 4,line = 1.5, cex = 1.3)
text(x = 500, y= 0.089, labels = expression(italic(g[s])), col = Cols[2], cex = 1.3)

axis(1, xlim = c(0, 1000), pos = 0, lwd = 2)
mtext("days",side = 1,line = 1, cex = 1.3)

dev.copy2pdf(file = "output/figures/Figure 1.pdf")
