library(fpCompare)

A=cbind(c(1,2,1),c(1,0,1))
old_A=A

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
      A[c(i, j), (i + 1):n] = G%*% A[c(i, j), (i + 1):n]
      Q[, c(i, j)] = Q[, c(i, j)]%*%t(G)
    }
    A[(i+1):m,i]=0
  }
  # return(list(Q = Q, R = A))
