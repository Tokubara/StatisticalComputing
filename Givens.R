library(fpCompare)

Givens<-function(A) {
  m = nrow(A)
  n = ncol(A)
  if (n > m) {
    # A满足行数>=列数
    stop("n>m, stop")
  }
  Q = diag(m)
  for (i in 1:n) {
    xi=A[i,i]
    for(j in (i+1):n) {
      if(j>n) break
      xj=A[j,i]
      .t = sqrt(xi ^ 2 + xj ^ 2)
      .c=xi/.t
      .s=xj/.t
      G=cbind(c(.c,-.s),c(.s,.c))
      # 用这个Q去作用
      A[c(i, j), (i):n] = G%*% A[c(i, j), (i):n] # fix:一开始,没有对A[,i]作用, 这会让A[,i]得不到改变的
      Q[, c(i, j)] = Q[, c(i, j)]%*%t(G)
    }
  }
  return(list(Q = Q, R = A))
}

A = randi(20,10,8)
res=Givens(A)
Q=res$Q
R=res$R
(A)%==%(Q%*%R)

