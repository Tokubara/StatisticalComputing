library(fpCompare)

householder<-function(A) {
  m = nrow(A)
  n = ncol(A)
  if (n > m) {
    # A满足行数>=列数
    stop("n>m, stop")
  }
  Q=diag(m)
  for(i in 1:n) {
    qi=A[i:m,i]
    .t = norm(qi, "2")
    if((.t)%==%(abs(qi[1]))) {
      next
    }
    qi[1] = qi[1] + sign(qi[1]) * .t
    qi=qi/norm(qi,"2") # 单位化, 但好像也未必需要
    A[i, i] = -sign(qi[1]) * .t
    A[(i+1):m,i]=0
    for(j in (i+1):n) {
      # 对其它列作用
      if(j>n) break # fix:if(j>n) break
      aj = A[i:m, j]
      .t=dot(qi,aj)
      aj=aj-2*.t*qi # fix:aj=aj-.t*qi
      A[i:m, j]=aj
    }
    H=diag(m)
    H[i:m,i:m]=diag(m-i+1)-2*qi%*%t(qi)
    Q=Q%*%t(H)
  }
  return(list(Q=Q,R=A))
}

householder.test<-function(A=NULL) {
  if(is.null(A)) {
    A=randi(20,10,8)
  }
  res=householder(A)
  Q=res$Q
  R=res$R
  print(all((Q%*%R)%==%(A)))
}

# householder.test()
