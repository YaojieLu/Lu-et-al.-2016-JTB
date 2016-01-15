
#data from Zhang et al. 2001
MAP <- c(
  1260, 1819, 645, 1639, 2978, 995, 1469, 1494, 1272, 1451, 1602, 596, 596, 701, 721, 788, 788, 836, 959, 1128, 1212, 1017, 2049, 2080, 1023, 1153, 1113, 1768, 2015, 1309, 704, 2648, 2425, 2851, 1983, 1800, 1750, 1854, 1611, 1193, 1207, 1333, 813, 813, 639, 457, 457, 1153, 1113, 2098, 2081, 1460, 1100, 579, 597
  )
E <- c(
  1115, 1363, 622, 1190, 1483, 898, 1204, 1200, 1127, 1199, 1295, 552, 515, 683, 696, 701, 706, 714, 863, 879, 939, 956, 1328, 1291, 938, 859, 823, 1059, 1154, 1033, 631, 1311, 1440, 1481, 1265, 1145, 1019, 1016, 1063, 1168, 1122, 1180, 727, 726, 568, 437, 439, 860, 823, 1394, 1295, 990, 910, 467, 442
  )
MAPE <- data.frame(MAP = MAP, E = E)
# the functional form follows Eq. 6 Zhang et al. 2001
nls1 <- nls(E ~ MAP*((1+w*E0/MAP)/(1+w*E0/MAP+MAP/E0)), start=list(w = 3.822, E0 = 1128.347), data=MAPE)
p <- predict_nls(nls1, from=min(MAP), to=max(MAP), interval="confidence")

data <- dvs
data$ca <- as.factor(data$ca)
data$k <- as.factor(data$k)
Cols <- c("red","darkgreen","blue")

# Figure 6
windows(8, 12)
par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(2, 4, 1.8, 2), mfrow=c(2,1))
# average E - average daily rainfall

plot(MAP, E/365, panel.first = {
  addpoly(p$x, p$lwr/365, p$upr/365)
  lines(p$x, p$pred/365)
},
type = "n",
xlim = c(0, 3500),
ylim = c(0, 8),
xlab = NA, ylab = expression(bar(italic(E))~(mm~day^-1)),
xaxt = "n",
cex.lab = 1.3)

plotBy(data$E*500 ~ data$MAP | k, data=subset(data, ca == 400),
       type='l', legend = FALSE, legendwhere="topleft",
       xlim = c(0, 3500),
       ylim = c(0, 8),
       xlab = NA, ylab = expression(bar(italic(E))~(mm~day^-1)),
       xaxt = "n",
       cex.lab = 1.3,
       col = Cols, add = T)

rect(0, 0, 457, 8, col = "white", border = NA)
rect(2978, 0, 3500, 8, col = "white", border = NA)

axis(1, xlim = c(0, 3500), pos = 0, lwd = 2, labels = FALSE)
text(40, 8*0.99, "a", adj = c(0, 1), col = "black", cex = 1.3)
box()

# E/P - average daily rainfall
par(mar=c(3.8, 4, 0, 2))

plot(MAP, E/MAP, panel.first = {
  addpoly(p$x, p$lwr/p$x, p$upr/p$x)
  lines(p$x, p$pred/p$x)
},
type = "n",
xlim = c(0, 3500),
ylim = c(0, 1),
xlab = NA, ylab = expression(bar(italic(E))/MAP),
xaxt = "n",
cex.lab = 1.3)

plotBy(EP ~ MAP | k, data=subset(data, ca == 400),
       type='l', legend = FALSE, legendwhere="topleft",
       xlim = c(0, 3500),
       ylim = c(0, 1),
       xlab = NA, ylab = expression(bar(italic(E))/MAP),
       xaxt = "n",
       cex.lab = 1.3,
       col = Cols, add = T)

rect(0, 0, 457, 1, col = "white", border = NA)
rect(2978, 0, 3500, 1, col = "white", border = NA)

legend("bottomleft", expression(italic(k==0.025), italic(k==0.05), italic(k==0.1), Zhang~italic(et~al.)~2001),
       col = c("red","darkgreen","blue", "black"), lty=c(1, 1, 1, 1), lwd=c(2, 2, 2, 3))

axis(1, xlim = c(0, 3500), pos = 0, lwd = 2)
mtext(expression(MAP~(mm~year^-1)),side = 1,line = 2.5, cex = 1.3)
text(40, 1*0.99, "b", adj = c(0, 1), col = "black", cex = 1.3)
box()

dev.copy2pdf(file = "output/figures/Figure 6.pdf")
