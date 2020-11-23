library(pracma)
f <- function(x) {
  (1 - x[1]) ^ 2 + 100 * (x[2] - x[1] * x[1]) ^ 2
}
grad <-function(x) { # 已测试
  return(c(  2*(x[1]-1)+400*x[1]*(x[1]^2-x[2]),200*(x[2]-x[1]^2)  ))
}
Hesse <- function(x) { # 已测试
  return(matrix(c(2+400*(3*x[1]^2-x[2]),-400*x[1],-400*x[1],200),nrow=2))
}

conjugate_gradient <- function(xk, f, grad, Hesse, maxit = 100000, tol = 1e-6) {
  .rk <- function(x) {
    -grad(x)
  }
  #它们都是接受向量为参数
  rk = .rk(xk)
  dk = rk
  n = length(xk)
  iter = 0
  while (sum(dk ^ 2) > tol) {
    if (iter > maxit) {
      cat("reach maxit, exit\n")
      break
    }
    for (k in range(n)) {
      A = Hesse(xk)
      denom = as.vector((dk %*% A %*% dk))
      ak = dot(dk, rk) / denom
      xk = xk + ak * dk
      if (k < n) {
        rk = .rk(xk)
        bk = as.vector((dk %*% A %*% rk)) / denom
        dk = rk - bk * dk
      }
    }
    iter = iter + 1
    rk = .rk(xk)
    dk = .rk(xk)
    # cat(xk, sum(dk ^ 2),"\n")
  }
  return(list(x = xk, f.x = f(xk), iter = iter))
}
conjugate_gradient(c(2,3),f,grad,Hesse)

Newton <- function(xk,f,grad, Hesse, maxit = 50000) {
  iter = 0
  dx = 1
  while (sum(dx ^ 2) > tol) {
    if (iter > maxit) {
      cat("reach maxit, exit\n")
      break
    }
    H = Hesse(xk)
    g = grad(xk)
    dx = solve(H) %*% g
    xk = xk - dx # 如果这里太慢了, 可以考虑直接代入公式.
    iter = iter + 1
  }
  return(list(x = xk, f.x = f(xk), iter = iter))
}

xk = c(2, 3)
Newton(c(2,3),f,grad,Hesse)

gradient_descent<-function(xk,f,grad) {
  iter = 0
  backtracking_step <- function(x, a = 0.3, b = 0.6) {
    step = 1
    grad_ = grad(x)
    f.x = f(x) # x在这个过程中是不更新的, grad也不更新
    max_iter = 300
    iter = 0
    const = a * sum(grad_ ^ 2)
    while (T) {
      iter = iter + 1
      if (iter > max_iter) {
        # cat("not found step\n")
        return(-1)
      }
      x_new = x - grad_ * step
      f.x_new = f(x_new)
      if (f.x - f.x_new >= const * step) {
        return(step)
      } else {
        # 虽然else多余
        step = step * b
      }
    }
  }
  while (T) {
    if (iter > maxit) {
      cat("reach maxit, exit\n")
      break
    }
    grad_ = grad(xk)
    step = backtracking_step(xk)
    if (step < 0) {
      # 修改了判停准则
      break
    }
    diff = step * grad_
    xk = xk - diff
    iter = iter + 1
  }
  return(list(x=xk,f.x=f(xk),iter=iter))
}

gradient_descent(c(2,3),f,grad)



