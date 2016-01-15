
#a
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

w <- seq(0,1,length=1e4)

ca <- 400
Low <- 365
High <- 3650
Infrequent <- 0.025
Frequent <- 0.1

if(!file.exists("output/data/Figure 4a.csv")){
  # L and I
  gwLI <- gw(ca=ca, k=Infrequent, MAP=Low, averA=dvs[dvs$ca==ca & dvs$k==Infrequent & dvs$MAP==Low, "A"])
  
  # L and F
  gwLF <- gw(ca=ca, k=Frequent, MAP=Low, averA=dvs[dvs$ca==ca & dvs$k==Frequent & dvs$MAP==Low, "A"])
  
  # H and I
  gwHI <- gw(ca=ca, k=Infrequent, MAP=High, averA=dvs[dvs$ca==ca & dvs$k==Infrequent & dvs$MAP==High, "A"])
  
  # H and F
  gwHF <- gw(ca=ca, k=Frequent, MAP=High, averA=dvs[dvs$ca==ca & dvs$k==Frequent & dvs$MAP==High, "A"])
  
  gwres <- data.frame(w=w, LI=gwLI, LF=gwLF, HI=gwHI, HF=gwHF)
  write.csv(gwres, "output/data/Figure 4a.csv", row.names=FALSE)
} else {
  gwres <- read.csv("output/data/Figure 4a.csv")
}

#b
# Eq. 4
pdff <- function(averA,
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
  
  #PDF
  Ev <- function(w){h*VPD*f2(w)}
  rEv <- function(w){1/Ev(w)}
  integralEv <- Vectorize(function(w){integrate(rEv, 0, w)$value})
  fnoc <- function(w){1/Ev(w)*exp(-gamma*w+k*integralEv(w))}
  integralfnoc <- integrate(fnoc, 0, 1)$value
  c <- 1/(integralfnoc+1/k)
  f <- function(w){c*fnoc(w)}
  pdf <- f(w2)
  return(pdf)
}

w2 <- seq(1/500, 1, length=500)

ca <- 400
Low <- 365
High <- 3650
Infrequent <- 0.025
Frequent <- 0.1


if(!file.exists("output/data/Figure 4b.csv")){
  # L and I
  pdfLI <- pdff(ca=ca, k=Infrequent, MAP=Low, averA=dvs[dvs$ca==ca & dvs$k==Infrequent & dvs$MAP==Low, "A"])
  
  # L and F
  pdfLF <- pdff(ca=ca, k=Frequent, MAP=Low, averA=dvs[dvs$ca==ca & dvs$k==Frequent & dvs$MAP==Low, "A"])
  
  # H and I
  pdfHI <- pdff(ca=ca, k=Infrequent, MAP=High, averA=dvs[dvs$ca==ca & dvs$k==Infrequent & dvs$MAP==High, "A"])
  
  # H and F
  pdfHF <- pdff(ca=ca, k=Frequent, MAP=High, averA=dvs[dvs$ca==ca & dvs$k==Frequent & dvs$MAP==High, "A"])
  
  pdfres <- data.frame(w=w2, LI=pdfLI, LF=pdfLF, HI=pdfHI, HF=pdfHF)
  write.csv(pdfres, "output/data/Figure 4b.csv", row.names=FALSE)
} else {
  pdfres <- read.csv("output/data/Figure 4b.csv")
}

#plots
#a
Cols <- c("red","red", "blue", "blue")
windows(8, 12)
par(mgp=c(2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(2, 4, 1.5, 2), mfrow=c(2,1))

with(gwres, {
  plot(w*100, LI, type="l", 
       xaxt = "n",
       xlab=NA, 
       ylab=expression(italic(g[s])~(mol~m^-2~s^-1)),
       xlim=c(0, 100),
       ylim=c(0, 0.6),
       cex.lab=1.3,
       col=Cols[1])
  lines(w*100, HI, lty=1,
        col=Cols[3])
  lines(w*100, LF, lty=2, 
        col=Cols[2])
  lines(w*100, HF, lty=2,
        col=Cols[4])
})

axis(1, xlim=c(0, 100), pos=0, lwd=2, labels=FALSE)
text(96.5, 0.6*0.99, "a", adj=c(0, 1), col="black", cex=1.3)

Legend <- c("Low MAP & Low frequency", "Low MAP & High frequency", "High MAP & Low frequency", "High MAP & High frequency")
legend("topleft", Legend, col=Cols, lty=c(1, 2, 1, 2), lwd=c(2, 2, 2, 2))

box()

#b
par(mar=c(3.5, 4, 0, 2))
with(pdfres, {
  plot(w*100, LI, col=Cols[1], type="l", 
       xaxt="n",
       xlab=NA, 
       ylab="Probability Density",
       xlim=c(0,100),
       ylim=c(0, 2),
       cex.lab=1.3)
  lines(w*100, HI, col=Cols[3], lty=1)
  lines(w*100, LF, col=Cols[2], lty=2)
  lines(w*100, HF, col=Cols[4], lty=2)
})

axis(1, xlim=c(0, 4000), pos=0, lwd=2)
mtext(expression(italic(w)~"(%)"),side=1,line=2.3, cex=1.3)
text(96.5, 2*0.99, "b", adj=c(0, 1), col = "black", cex=1.3)

box()

dev.copy2pdf(file = "output/figures/Figure 4.pdf")
