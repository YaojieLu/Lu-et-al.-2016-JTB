
data <- dvs
data$ca <- as.factor(data$ca)
data$k <- as.factor(data$k)
Cols <- c("red","darkgreen","blue")

#plots
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 1.5, 2), mfrow=c(1,1))
# the effect (ratio) of elevated ca on A & E
ca400medium <- subset(data, ca==400 & k==0.05)
ca800medium <- subset(data, ca==800 & k==0.05)
a <- length(ca400medium$MAP)
caA <- data.frame(MAP=numeric(length=a), ratio=numeric(length=a), stringsAsFactors=FALSE)
caA$MAP <- ca400medium$MAP
caA$ratio <- ca800medium$A/ca400medium$A
caE <- data.frame(MAP=numeric(length=a), ratio=numeric(length=a), stringsAsFactors=FALSE)
caE$MAP <- ca400medium$MAP
caE$ratio <- ca800medium$E/ca400medium$E

plot(caA$MAP, caA$ratio, type="l",
     xlab="", ylab=expression(Ratio),
     xlim=c(0, 4000), ylim=c(0.5, 2.5),
     xaxt="n",     
     col=c("red"), cex.lab=1.3
)

lines(caA$MAP, caE$ratio, lty=1, col=c("blue"))

abline(h=1, col=c("black"), lwd=2, lty=2)

axis(1, xlim=c(0, 4000), pos=0.5, lwd=2)
mtext(expression(MAP~(mm~year^-1)),side=1,line=2.5, cex=1.3)

legend("topright", expression(italic(bar(A)[c[a]==800]/bar(A)[c[a]==400]), italic(bar(E)[c[a]==800]/bar(E)[c[a]==400])),
       col=c("red", "blue"), lty=c(1, 1), lwd=c(2, 2))

box()

dev.copy2pdf(file="output/figures/Figure 8.pdf")
