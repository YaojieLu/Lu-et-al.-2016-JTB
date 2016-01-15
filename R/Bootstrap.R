
source("Bootstrap predict_nls.R")
source("Bootstrap functions.R")

f <- function(t, k=0.5)1 - exp(-k*t)
x <- seq(0,12,length=101)
y <- f(x) + rnorm(101,sd=0.1)
dfr <- data.frame(x=x,y=y)
nls1 <- nls(y ~ b0 + A*(1-exp(-k*x)), start=list(b0=0,A=1,k=1), data=dfr)

p <- predict_nls(nls1, from=min(x), to=max(x), interval="confidence")

plot(x,y)
with(p, {
  lines(x, pred, lty=5)
  lines(x, lwr, lty=4, col="red")
  lines(x, upr, lty=4, col="red")
})


plot(x,y, panel.first={
  addpoly(p$x, p$lwr, p$upr)
  lines(p$x, p$pred)
})
