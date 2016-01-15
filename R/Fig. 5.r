
# calculating gsmax (i.e., gs(1)) and linearity
f <- function(averA,
              ca, k, MAP,
              Vcmax=50, cp=30, Km=703, Rd=1, VPD=0.02, p=43200, LAI=1, nZ=0.5,
              a=1.6, l=1.8e-5, MDP=MAP/365, gamma=1/((MDP/k)/1000)*nZ, h=l*a*LAI/nZ*p){
  
  integrand <- function(g){(4*h*(cp+Km)*LAI*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*g)/(((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))*(ca*k*LAI*Rd+k*Km*LAI*Rd-ca*k*LAI*Vcmax+2*cp*k*LAI*Vcmax+k*Km*LAI*Vcmax+gamma*h*LAI*Rd^2*VPD-2*gamma*h*LAI*Rd*Vcmax*VPD+gamma*h*LAI*Vcmax^2*VPD+ca^2*k*LAI*g+2*ca*k*Km*LAI*g+k*Km^2*LAI*g+ca*gamma*h*LAI*Rd*VPD*g+gamma*h*Km*LAI*Rd*VPD*g-ca*gamma*h*LAI*Vcmax*VPD*g+2*cp*gamma*h*LAI*Vcmax*VPD*g+gamma*h*Km*LAI*Vcmax*VPD*g-ca*k*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))-k*Km*LAI*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+2*gamma*h*averA*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))+gamma*h*LAI*Rd*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))-gamma*h*LAI*Vcmax*VPD*sqrt((Rd-Vcmax+(ca-Km)*g)^2+4*g*(Km*Rd+cp*Vcmax+ca*Km*g))))}
  slope <- function(g)1/integrand(g)
  f1 <- function(w){
    w <- w
    O <- function(g){integrate(integrand, 0, g)$value-w}
    u <- uniroot(O, c(0,1), tol=1e-10, extendInt = "yes")
    return(u$root)
  }
  
  f2 <- Vectorize(f1)
  
  gsmax <- f2(1)
  slopemin <- optim(0.5, slope, method="Brent", lower=0, upper=gsmax)$value
  linearity <- slopemin/gsmax
  
  return(c(gsmax, linearity))
}

#environmental conditions
ca <- c(400)
k <- c(0.025, 0.1)
MAP <- seq(0.5, 10, by=0.5)*365
env <- as.vector(expand.grid(ca, k, MAP))

res1 <- matrix(nrow=nrow(env), ncol=2)

for(i in 1:nrow(env)){
  res1[i, ] <- f(ca=env[i, 1], k=env[i, 2], MAP=env[i, 3], averA=dvs[dvs$ca==env[i, 1] & dvs$k==env[i, 2] & dvs$MAP==env[i, 3], "A"])
}

res <- cbind(env, res1)
colnames(res) <- c("ca", "k", "MAP", "gsmax", "linearity") 
write.csv(res, "output/data/Figure 5.csv", row.names = FALSE)

Cols <- c("red", "blue")

# plots
windows(8, 12)
par(mgp = c(2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(2, 4, 1.5, 2), mfrow=c(2,1))

plotBy(gsmax ~ MAP | k, data=res,
       type='l', legend = FALSE, legendwhere="topleft",
       xlab = NA, 
       ylab = expression(italic(g[smax])~(mol~m^-2~s^-1)),
       xlim = c(0, 4000),
       ylim = c(0, 0.6),
       xaxt = "n",
       cex.lab = 1.3,
       col = Cols)

axis(1, xlim = c(0, 4000), pos = 0, lwd = 2, labels = FALSE)
text(40, 0.6*0.99, "a", adj = c(0, 1), col = "black", cex = 1.3)

par(mar=c(3.8, 4, 0, 2))

plotBy(linearity ~ MAP | k, data=res,
       type='l', legend = FALSE, legendwhere="topleft",
       xlab = NA, 
       ylab = expression(Linearity~of~italic(g[s])(italic(w))),
       xlim = c(0, 4000),
       ylim = c(0, 1),
       xaxt = "n",
       cex.lab = 1.3,
       col = Cols)

axis(1, xlim = c(0, 4000), pos = 0, lwd = 2)
mtext(expression(MAP~(mm~year^-1)),side = 1,line = 2.5, cex = 1.3)
text(40, 1*0.99, "b", adj = c(0, 1), col = "black", cex = 1.3)
box()

legend("bottomright", expression(italic(k==0.025), italic(k==0.1)),
       col = c("red","blue"), lty=c(1, 1), lwd=c(2, 2))

dev.copy2pdf(file = "output/figures/Figure 5.pdf")
