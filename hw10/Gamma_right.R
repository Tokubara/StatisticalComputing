set.seed(543) 
lambda <- c(.6, .25, .15) #rate is 1/(2lambda) 
m<-2000
lam <- sample(lambda, size = 2000, replace = TRUE) 
y <- rgamma(m, shape = .5, rate = 1/(2*lam)) 
N <- 10000 #max. number of iterations 
L <- c(.5, .4, .1) #initial est. for lambdas 
tol <- .Machine$double.eps^0.5 
L.old <- L + 1 
for (j in 1:N) { 
  f1 <- dgamma(y, shape=1/2, rate=1/(2*L[1])) 
  f2 <- dgamma(y, shape=1/2, rate=1/(2*L[2])) 
  f3 <- dgamma(y, shape=1/2, rate=1/(2*L[3])) 
  py <- f1 / (f1 + f2 + f3) #posterior prob y from 1 
  qy <- f2 / (f1 + f2 + f3) #posterior prob y from 2 
  ry <- f3 / (f1 + f2 + f3) #posterior prob y from 3 
  mu1 <- sum(y * py) / sum(py) #update means 
  mu2 <- sum(y * qy) / sum(qy) 
  mu3 <- sum(y * ry) / sum(ry) 
  L <- c(mu1, mu2, mu3) #update lambdas 
  L <- L / sum(L) 
  if (sum(abs(L - L.old)/L.old) < tol) break 
  L.old <- L 
}
