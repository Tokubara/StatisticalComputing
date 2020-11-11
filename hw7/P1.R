library(fpCompare)

d <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44, 3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75, 0.27, 43.21)

l_f <- function(theta) {
  sum(-log(pi * (1 + (d - theta) ^ 2)))
}

neg_l_f <- function(theta) {
  sum(log(pi * (1 + (d - theta) ^ 2)))
}

# 作图
theta = seq(-10, 50, by = 0.05)
l_c = sapply(theta, l_f)
plot(x = theta, y = l_c, pch = 19)

dl_f <- function(x) {
  # ppt上设的是0.01, 太小了确实有抵消的问题
  (l_f(x + dx) - l_f(x - dx)) / (2 * dx)
}

ddl_f <- function(x) {
  (dl_f(x + dx) - dl_f(x - dx)) / (2 * dx)
}

# 二分法
bisection <- function(l_start, r_start, tol = 1E-12) {
  l = l_start
  r = r_start
  ITER_UPPER = 5000
  iter_n = 0
  l_val = dl_f(l) # >0
  r_val = dl_f(r) # <0
  while (TRUE) {
    mid = (l + r) / 2
    iter_n = iter_n + 1
    if (iter_n > ITER_UPPER) {
      stop('not found solution')
    }
    mid_val = dl_f(mid)
    if (abs(mid_val) < tol) {
      break
    } else if (mid_val > 0) {
      l = mid
    } else {
      r = mid
    }
  }
  return(mid)
}

bisection(-1, 1, 1E-16)

newton <- function(x0 = 1, N = 50, tol = 1E-12, dx = 0.01, fisher_info = FALSE) {
  iter <- integer(1)
  while (TRUE) {
    ddx <- ddl_f(x0)
    if (ddx %==% 0) {
      stop("导数为0")
    }
    dx <- dl_f(x0)
    x1 <- x0 - dx / ifelse(fisher_info, -0.5 * length(d), ddx) # bug:fisher-scoring-method 没有加上-
    iter = iter + 1
    if ((abs(x1 - x0) < tol)) {
      return(c(x1, iter))
    }
    x0 = x1 # bug:忘了更新x0
    if (iter > N) {
      stop("cannot find solution")
    }
  }
}

newton_start <- c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38)
sapply(newton_start, function(x) try(newton(x, N = 50)))
sapply(newton_start, function(x) try(newton(x, fisher_info = T, N = 50000)))

newton(mean(d))

newton(-0.2, N = 500000, fisher_info = T)

abs(1 + 0.1 * ddl_f(-0.1922947))
