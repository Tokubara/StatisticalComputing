library(pracma)
sigma.svd <- numeric(30)
sigma.square <- numeric(30)
for(m in 1:30) {
  A <- matrix(0, m, m)
  diag(A) <- 0.1
  A[upper.tri(A)] <- 1
  i = m
  sigma.svd[i] = svd(A)$d[m]
  sigma.square[i] = sqrt(min(eigen(t(A)%*%A)$values))
}


sigma.svd - sigma.square
col.square = rep(1, 30)
col.square[is.nan(sigma.square)] = 2
plot(1:30, log(sigma.svd), pch = 19, col = 3,ylim=c(-50,1))
lines(1:30, type="p",ifelse(is.nan(sigma.square), 0, log(sigma.square)), pch = 16, col = col.square)





