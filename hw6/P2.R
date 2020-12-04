library(pracma)
library(fpCompare)

set.seed(3)
par(mfrow = c(5, 2))
for(i in 1:10) {
  # i不可再用?
  par(pch=19)
  A=randn(10)-2*diag(10)
  # log_t=seq(-4.6,3,by=0.1)
  # t=exp(log_t)
  t=seq(0,20,by=0.1)
  f=function(x) {
    return(norm(exp(A*x)))
  }
  log_f_t=log(sapply(t,f))
  plot(t, log_f_t,cex=0.5)
  abline(b=max(Re(eig(A))),a=0,lty="dashed",col=2)
}

