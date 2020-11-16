f<-function(x) { # x是向量
  return((10*x[1]^2+x[2]^2)/2)
}
grad_f <- function(x) {
  return(c(10*x[1],x[2]))
}

x_=c(2,3)
tol1=1e-4
tol2=1e-6
iter=0
MAX_ITER=500

while(T) {
  if(iter>MAX_ITER) {
    stop("not found solution")
  }
  grad_=grad_f(x_)
  step=0.01  # 这里需要改成switch
  diff = step * grad_
  x_=x_-diff
  if (sqrt(sum(diff ^ 2)) < tol1 * (sqrt(sum(x_ ^ 2)) + tol2)) { # 感觉判停准则不是很理想, 因为step就限制了不可能很大
    print("x=")
    print(x)
    print("f")
    print(f(x))
    break
  }
  iter=iter+1
}