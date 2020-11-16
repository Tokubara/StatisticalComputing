f<-function(x) { # x是向量
  return((10*x[1]^2+x[2]^2)/2)
}
grad_f <- function(x) {
  return(c(10*x[1],x[2]))
}

printf <- function(...) cat(sprintf(...))

x_=c(2,3)
tol1=1e-4
tol2=1e-6
iter=1
MAX_ITER=5000

x1_plot=numeric(MAX_ITER)
x2_plot = numeric(MAX_ITER)

exact_step<-function(x_,grad_) {
  return((10*grad_[1]*x_[1]+x_[2]*grad_[2])/(10*grad_[1]^2+grad_[2]^2))
}

backtracking_step <- function(x, a=0.3,b=0.6) {
  step=1
  grad_ = grad_f(x)
  f.x = f(x) # x在这个过程中是不更新的, grad也不更新
  max_iter=300
  iter=0
  const = a * sum(grad_ ^ 2)
  while(T) {
    iter=iter+1
    if(iter>max_iter) {
      stop("not found step")
    }
    x_new = x - grad_ * step
    f.x_new = f(x_new)
    if(f.x-f.x_new>=const*step) {
      return(step)
    } else { # 虽然else多余
      step=step*b
    }
  }
}
# method=""
while(T) {
  if(iter>MAX_ITER) {
    stop("not found solution")
  }
  grad_=grad_f(x_)
  # step=0.01 # 定步长
  # step = exact_step(x_, grad_)
  step = backtracking_step(x_)
  diff = step * grad_
  x_=x_-diff
  x1_plot[iter]=x_[1]
  x2_plot[iter]=x_[2]
  if (sqrt(sum(diff ^ 2)) < tol1 * (sqrt(sum(x_ ^ 2)) + tol2)) { # 感觉判停准则不是很理想, 因为step就限制了不可能很大
    print("x=")
    print(x_)
    print("f")
    print(f(x_))
    print("iter=")
    print(iter)
    break
  }
  iter=iter+1
}
plot(x1_plot,x2_plot)
