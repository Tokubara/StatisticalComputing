# 题目是, 已知Y是混合分布, 混合概率都是1/3, 每一个分布分别是Ga(1/2,1/(2*λj)). 现在需要估出λj.

# 造数据
set.seed(543)
m <- 2000
lambda <- c(.6, .25, .15) #rate is 1/(2lambda)
lam <- sample(lambda, size = 2000, replace = TRUE)
y <- rgamma(m, shape = .5, rate = 1 / (2 * lam))

N <- 10000 #max. number of iterations 
L <- c(.5, .4, .1) #initial est. for lambdas 
tol <- .Machine$double.eps ^ 0.5

test_debug.f<-function() {
  iter_counter = 0
  while (iter_counter < N) {
    L.old = L
    iter_counter = iter_counter + 1
    # 先计算出L
    likelihood.m = sapply(L, function(lambda) { dgamma(y, 0.5, 1 / (2 * lambda)) })
    likelihood.m=likelihood.m/rowSums(likelihood.m)
    L = apply(likelihood.m, MARGIN = 2, function(col) dot(col, y) / sum(col))
    L=L/sum(L)
    if (sum(abs(L.old - L) / L.old) < tol) {
      break
    }
  } 
  return(list(lambda=L, iter_num=iter_counter))
}

res=test_debug.f()
# 用的是相对变化, 那么需要存下L.old
# .t=cbind(c(1,2),c(3,4))
# .t/rowSums(.t)


