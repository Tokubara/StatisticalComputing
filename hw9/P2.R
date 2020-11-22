grad <-function(x) {
  return(c(  2*(x[1]-1)+400*x[1]*(x[1]^2-x[2]),200*(x[2]-x[1]^2)  ))
}
Hesse <- function(x) {
  return(matrix(c(2+400*(3*x[1]^2-x[2]),-400*x[1],-400*x[1],200),nrow=2))
}
.rk<-function(x) {
  -grad(x)
}
#它们都是接受向量为参数

xk=c(2,3)
rk = .rk(xk)
dk=rk
n=2
for(k in range(10*n)) {
  A=Hesse(xk)
  denom = (dk %*% A %*% dk)
  ak = dot(dk, rk) / denom
  xk=x0+ak*dk
  rk=.rk(xk)
  bk = (dk %*% A %*% rk) / denom
  dk=rk-bk*dk
}

