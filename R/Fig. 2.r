
# a simulation model: generating data of daily soil water content, stomatal
# conductance, and marginal water cost during a prolonged drought period
f <- function(averA,
              ca, k, MAP,
              w0=1,
              Vcmax=50, cp=30, Km=703, Rd=1, VPD=0.02, p=43200, LAI=1, nZ=0.5,
              a=1.6, l=1.8e-5, MDP=MAP/365, gamma=1/((MDP/k)/1000)*nZ, h=l*a*LAI/nZ*p){
  
  gsw1 <- function(w){
    w <- w
    O <- function(g){
      integrand <- function(g){(4*h*(cp+Km)*LAI*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*g)/(((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))*(ca*k*LAI*Rd+k*Km*LAI*Rd-ca*k*LAI*Vcmax+2*cp*k*LAI*Vcmax+k*Km*LAI*Vcmax+gamma*h*LAI*Rd^2*VPD-2*gamma*h*LAI*Rd*Vcmax*VPD+gamma*h*LAI*Vcmax^2*VPD+ca^2*k*LAI*g+2*ca*k*Km*LAI*g+k*Km^2*LAI*g+ca*gamma*h*LAI*Rd*VPD*g+gamma*h*Km*LAI*Rd*VPD*g-ca*gamma*h*LAI*Vcmax*VPD*g+2*cp*gamma*h*LAI*Vcmax*VPD*g+gamma*h*Km*LAI*Vcmax*VPD*g-ca*k*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))-k*Km*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+2*gamma*h*averA*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+gamma*h*LAI*Rd*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))-gamma*h*LAI*Vcmax*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))))}
      integrate(integrand, 0, g)$value-w
    }
    u <- uniroot(O, c(0,1), tol=1e-10, extendInt = "yes")
    return(u$root)
  }
  
  gsw2 <- Vectorize(gsw1)
  
  g <- vector("numeric")
  lambda <- vector("numeric")
  E <- vector("numeric")
  w <- vector("numeric")
  w[1] <- w0
  day <- vector("numeric")
  day[1] <- 0
  
  for(i in 2:1e5){
    g[i] <- gsw2(w[i-1])
    lambda[i] <- (2*1.6*VPD)/(-((2*Vcmax*(-ca+2*cp+Km)+2*(ca+Km)*(g[i]*(ca+Km)+Rd))/(2*sqrt(2*g[i]*Vcmax*(-ca+2*cp+Km)+(g[i]*(ca+Km)+Rd)^2-2*Rd*Vcmax+Vcmax^2)))+ca+Km)
    E[i] <- h*VPD*g[i]
    w[i] <- w[i-1]-E[i]
    day[i] <- day[i-1]+1
    if (w[i]<=0) break
  }
  
  res <- data.frame(day=day[-1], g=g[-1], w=w[-1], lambda=lambda[-1])
  
  return(res)
}

ca <- 400
k <- 0.025
MAP <- 1095

averA <- dvs[dvs$ca==ca & dvs$k==k & dvs$MAP==MAP, "A"]
res <- f(averA, ca, k, MAP)
write.csv(res, "output/data/Figure 2.csv", row.names = FALSE)

#plots
data <- head(res, -1)
colnames(data) <- c("day", "g", "w", "lambda")

windows(8, 6)
par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(3.3, 7, 1, 3.7), mfrow=c(1,1))

#g
plot(data$day, data$g, type = "l",
     xlim = c(1, tail(data$day, 1)), ylim = c(0, 0.2),
     xaxt = "n",
     xlab = NA, ylab = NA,
     cex.lab = 1.3,
     col = "darkgreen")

mtext(expression(italic(g[s])~(mol~m^-2~s^-1)),side = 2,line = 1.9, cex = 1.3)
text(145, 0.08, expression(italic(g[s])), adj = c(0, 1), col = "darkgreen", cex = 1.3)

#lambda
par(new=TRUE)

plot(data$day, data$lambda*10^6, type="l",
     xlim = c(1, tail(data$day, 1)), ylim = c(0, 3000),
     xaxt = "n", yaxt = "n",
     xlab = NA, ylab = NA,
     cex.lab = 1.3,
     col = "orange"
)

axis(4, ylim = c(0, 3000), pos = tail(data$day, 1), lwd = 2)
mtext(expression(italic(lambda)~(mol~(H[2]*O)/mol~(CO[2]))),side = 4,line = 2.4, cex = 1.3)
text(145, 470, expression(italic(lambda)), adj = c(0, 1), col = "orange", cex = 1.3)

#w
par(new=TRUE)

plot(data$day, data$w*100, type = "l",
     xlim = c(1, tail(data$day, 1)), ylim = c(0, 100),
     xaxt = "n", yaxt = "n",
     xlab = NA, ylab = NA,
     cex.lab = 1.3,
     col = "blue")

axis(2, ylim = c(0, 100), pos = -37, lwd = 2)
mtext(expression(italic(w~"(%)")),side = 2,line = 5.5, cex = 1.3)
text(145, 30, expression(italic(w)), adj = c(0, 1), col = "blue", cex = 1.3)

axis(1, xlim = c(0, tail(data$day, 1)), pos = 0, lwd = 2)
mtext("days",side = 1,line = 2, cex = 1.3)
box()

dev.copy2pdf(file = "output/figures/Figure 2.pdf")
