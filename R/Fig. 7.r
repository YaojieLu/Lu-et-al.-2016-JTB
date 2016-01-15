
data1 <- read.csv("data/NPP.csv")
colnames(data1) <- c("MAP", "TNPP")
MAP <- data1$MAP
TNPP <- data1$TNPP

# 3000
nls1 <- nls(TNPP ~ a*MAP^b/exp(c*MAP), start=list(a = 0.5346608789, b = 1.0787801564, c= 0.0003562646), data=data1)
p <- predict_nls(nls1, from=min(MAP), to=max(MAP), interval="confidence")

data <- dvs
data$ca <- as.factor(data$ca)
data$k <- as.factor(data$k)
Cols <- c("red","darkgreen","blue")

# Figure 5
windows(8, 6)
par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(4, 4, 1, 4), mfrow=c(1,1))
# average A - average daily rainfall
plotBy(data$A ~ data$MAP | k, data=subset(data, ca == 400),
       type='l', legend = FALSE, legendwhere="topleft",
       xlim = c(0, 3000),
       ylim = c(0, 15),
       xlab = NA,
       ylab = expression(bar(italic(A))~(mu*mol~m^-2~s^-1)),
       xaxt = "n",
       yaxt = "n",
       cex.lab = 1.3,
       col = Cols)

axis(1, xlim = c(0, 3000), pos = 0, lwd = 2)
axis(2, xlim = c(0, 15), pos = 0, lwd = 2, at = c(0, 5, 10, 15))
mtext(expression(MAP~(mm~year^-1)),side = 1,line = 2.5, cex = 1.3)

par(new=TRUE)

p1 <- data.frame(x = p$x, pred = p$pred, lwr = p$lwr, upr = p$upr)
p2 <- subset(p1, x<=2500)

plot(MAP, TNPP, panel.first = {
  addpoly(p2$x, p2$lwr, p2$upr)
  lines(p2$x, p2$pred)
},
type = "n",
xlim = c(0, 3000),
ylim = c(0, 1420),
xlab = "", ylab = "",
axes = FALSE)

rect(0, 0, 182.5, 1420, col = "white", border = NA)
rect(2480.152, 0, 3000, 1420, col = "white", border = NA)

legend("bottomright", expression(italic(k==0.025), italic(k==0.05), italic(k==0.1), Del~Grosso~italic(et~al.)~2008),
       col = c("red","darkgreen","blue", "black"), lty=c(1, 1, 1, 1), lwd=c(2, 2, 2, 3))

axis(4, ylim = c(0, 1200), pos = 3000, lwd = 2, at = c(0, 500, 1000))
mtext(expression(TNPP~(gC~m^-2~yr^-1)),side = 4,line = 2.8, cex = 1.3)

box()

dev.copy2pdf(file = "output/figures/Figure 7.pdf")