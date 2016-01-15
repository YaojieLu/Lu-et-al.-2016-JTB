
# Lin et al. 2015
g1data <- read.csv("data/g1.csv")
colnames(g1data) <- c("MAP", "g1")

# Prentice et al. 2011 (simplified)
data1 <- read.csv("data/delta13.csv")
colnames(data1) <- c("MAP", "delta13")
MAP <- data1$MAP
cica <- (-8.3-data1$delta13-4.4)/22.6
nls1 <- nls(cica ~ a+b*MAP, start=list(a = 0.411991, b = 0.000619469), data=data1)
p <- predict_nls(nls1, from=min(MAP), to=max(MAP), interval="confidence")

data <- dvs
data$ca <- as.factor(data$ca)
data$k <- as.factor(data$k)
Cols <- c("red","darkgreen","blue")

#plots
windows(8, 18)
par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(1.7, 4, 2.6, 2), mfrow=c(3,1))

# ci/ca - average daily rainfall
plot(MAP, cica, panel.first = {
  addpoly(p$x, p$lwr, p$upr)
  lines(p$x, p$pred)
},
type = "n",
xlim = c(0, 2000),
ylim = c(0, 1),
xlab = NA, ylab = expression(italic(c[i]/c[a])),
xaxt = "n",
cex.lab = 1.3)

plotBy(cica ~ MAP | k, data=subset(data, ca == 400 & MAP <= 1825),
       type='l', legend = FALSE, legendwhere="topleft",
       xlim = c(0, 2000),
       ylim = c(0, 1),
       xlab = NA, ylab = expression(italic(c[i]/c[a])),
       xaxt = "n",
       cex.lab = 1.3,
       col = Cols, add = TRUE)

axis(1, xlim = c(0, 100), pos = 0, lwd = 2, labels = FALSE)
text(1900, 1*0.99, "a", adj = c(0, 1), col = "black", cex = 1.3)

legend("bottomright", expression(italic(k==0.025), italic(k==0.05), italic(k==0.1), Prentice~italic(et~al.)~2011),
       col = c("red","darkgreen","blue", "black"), lty=c(1, 1, 1, 1), lwd = c(2, 2, 2, 2))

# lambda
par(mar=c(3, 4, 1.3, 2))
plotBy(lambda ~ MAP | k, data=subset(data, ca == 400 & MAP <= 1825),
       type='l', legend = FALSE, legendwhere="topleft",
       xlim = c(0, 2000),
       ylim = c(0, 0.006),
       xlab = NA, ylab = expression(italic(lambda)~(mol~(H[2]*O)/mol~(CO[2]))),
       xaxt = "n",
       yaxt = "n",
       cex.lab = 1.3,
       col = Cols)

axis(1, xlim = c(0, 100), pos = 0, lwd = 2, labels = FALSE)
axis(2, ylim = c(0, 0.006), pos = 0, lwd = 2, at = c(0, 0.002, 0.004, 0.006))
text(1900, 0.006*0.99, "b", adj = c(0, 1), col = "black", cex = 1.3)

legend("bottomright", expression(italic(k==0.025), italic(k==0.05), italic(k==0.1)),
       col = c("red","darkgreen","blue"), lty=c(1, 1, 1), lwd = c(2, 2, 2))

# g1 - average daily rainfall
par(mar=c(4.3, 4, 0, 2))
plot(g1data$MAP, g1data$g1,
     xlim = c(0, 2000), ylim = c(0, 15),
     xaxt = "n", yaxt = "n", pch = c(19),
     xlab = NA, ylab = NA)

par(new=TRUE)

plotBy(g1 ~ MAP | k, data=subset(data, ca == 400 & MAP <= 1825),
       type='l', legend = FALSE, legendwhere="topleft",
       xlim = c(0, 2000),
       ylim = c(0, 15),
       xlab = NA, ylab = expression(bar(italic(g[1]))),
       xaxt = "n",
       yaxt = "n",
       cex.lab = 1.3,
       col = Cols)

axis(1, xlim = c(0, 2000), pos = 0, lwd = 2)
axis(2, ylim = c(0, 15), pos = 0, lwd = 2, at = c(0, 5, 10, 15))
mtext(expression(MAP~(mm~year^-1)),side = 1,line = 3.3, cex = 1.3)
text(1900, 15*0.99, "c", adj = c(0, 1), col = "black", cex = 1.3)

legend("topleft", expression(italic(k==0.025), italic(k==0.05), italic(k==0.1), Lin~italic(et~al.)~2015),
       col = c("red","darkgreen","blue", "black"), lty=c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), lwd = c(2, 2, 2, NA))

box()

dev.copy2pdf(file = "output/figures/Figure 9.pdf")
