library(pracma)
grad <-function(x) { # 已测试
  return(c(  2*(x[1]-1)+400*x[1]*(x[1]^2-x[2]),200*(x[2]-x[1]^2)  ))
}
Hesse <- function(x) { # 已测试
  return(matrix(c(2+400*(3*x[1]^2-x[2]),-400*x[1],-400*x[1],200),nrow=2))
}
.rk<-function(x) {
  -grad(x)
}
#它们都是接受向量为参数

maxit = 1000000
tol=1e-6

xk=c(2,3)
rk = .rk(xk)
dk=rk
n=2
iter=0
while(sum(dk^2)>tol) {
  if(iter>maxit) { 
    cat("reach maxit, exit\n")
    break
  }
  for (k in range(n)) {
    A = Hesse(xk)
    denom = as.vector((dk %*% A %*% dk))
    ak = dot(dk, rk) / denom
    xk = xk + ak * dk
    if(k<n) {
      rk = .rk(xk)
      bk = as.vector((dk %*% A %*% rk)) / denom
      dk = rk - bk * dk
    }
  }
  iter=iter+1
  rk=.rk(xk)
  dk=.rk(xk)
  # cat(xk, sum(dk ^ 2),"\n")
}

