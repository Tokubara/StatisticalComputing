vandermonde <- function(v, n) {
  # 向量为v, 长度为n, 得到的是大小为length(v)*n的矩阵
  return(t(sapply(v, function(x) x ^ (0:(n - 1)))))
}

A = vandermonde(seq(0, 1, len = 25), 15)
cond(A) # 36302854020
cond(t(A) %*% A) # 9.577515e+17

P2 <- function(m) {
  # 参数m表示mode
  res = QR(A, m)
  Q = res$Q
  R = res$R
  n = ncol(Q)
  # browser()
  return(c(cond(A - Q %*% R), cond(diag(n) - t(Q) %*% Q)))
}

m = c("M", "C", "G", "H")
P2(m[1]) # rank < n, stop
P2(m[2]) # Inf 1.060788e+16
P2(m[3]) # 2.674211e+13 Inf
P2(m[4]) # 76.18462 39.93274