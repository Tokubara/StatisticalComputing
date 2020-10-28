pacman::p_boot()
pacman::p_load(fpCompare, pracma)

.f <- function(p, i) {
  # 为什么需要这个函数, 因为as.integer(names(which.max(abs())))
  if (all(p == 0)) {
    return(U[, i])
  } else {
    return(U[-p, i])
  }
}

MGS <- function(A) {
  tol = 1e-5 # 误差限
  m = nrow(A)
  n = ncol(A)
  if (n > m) {
    # A满足行数>=列数
    stop("n>m, stop")
  }
  Q = A
  R = matrix(0, n, n)
  for (i in 1:n) {
    qi = Q[, i]
    norm2 = norm(qi, "2")
    if (abs(norm2) < tol) {
      stop("rank < n, stop") # 如果列不满秩, 报错
    }
    R[i, i] = norm2
    qi = qi / norm2
    Q[, i] = qi
    for (j in (i + 1):n) {
      # 没办法, 这样写很不美观
      if (j > n) break
      # 去改变所有
      R[i, j] = dot(Q[, j], qi)
      Q[, j] = Q[, j] - qi * R[i, j]
    }
  }
  return(list(Q = Q, R = R))
}

CGS <- function(A) {
  tol = 1e-5 # 误差限
  m = nrow(A)
  n = ncol(A)
  if (n > m) {
    # A满足行数>=列数
    stop("n>m, stop")
  }
  Q = A
  R = matrix(0, n, n)
  for (i in 1:n) {
    qi = Q[, i]
    for (j in 1:(i - 1)) {
      # bug:1:(i-1)会出问题
      if (j > i - 1 || j < 1) break
      # 去改变所有
      R[j, i] = dot(Q[, j], Q[, i])
      qi = qi - Q[, j] * R[j, i]
    }
    norm2 = norm(qi, "2")
    if (abs(norm2) < tol) {
      stop("rank < n, stop") # 如果列不满秩, 报错
    }
    R[i, i] = norm2
    qi = qi / norm2
    Q[, i] = qi
  }
  return(list(Q = Q, R = R))
}

A = matrix(c(1,4,4,2,4,6,2,2,4),nrow=3)
n = nrow(A)
N_C = 1:n # 为了表明这是个常量, 而且是向量
pivoting=TRUE

p=integer(n) # 虽然只需要n-1
# 为debug好看, 先保证A没有修改
L=diag(n)
U=A
rownames(U)=1:3
if(n!=ncol(A)) {
  stop("A is not n*n")
}
for(i in 1:(n-1)) {
  # i=2 # 调试
  # 因为一共要进行n-1次
  if(pivoting) {
    p[i] = as.integer(names(which.max(abs(.f(p,i))))) # 找出第i列绝对值最大者 # fix:p[i]=which.max(abs(A[,i])) A没有改变, 应该是U # ifelse(all(p == 0), U[, i], U[-p, i])不行, 因为它会去掉names
  } else {
    p[i]=i
  }
  qi = U[p[i], (i + 1):n] # perf:qi = A[p[i], i:n]:因为对角元不参与运算 # fix:这里应该是U而不是A
  aii=qi[1] # 就是这个对角元
  for (j in N_C[-p]) {
    qij=U[j,i]/aii
    U[j, (i + 1):n] = U[j, (i + 1):n] - qij * qi
    L[j,i]=qij
  }
}





# test函数似乎可以继续用
test <- function(mode = c("c", "m")) {
  # c表示CGS, m表示MGS. # 0表示测试QR=A以及Q'Q=I, 1表示与R的QR作对比
  mode = match.arg(mode)
  if (mode == "c") {
    f = CGS
  } else {
    f = MGS
  }
  # 测试不满秩
  # A = randi(2, 3, 3)
  # A[, 3] = 2 * A[, 1] - 3 * A[, 2] # 使它不满秩
  # f(A) # 报错:不满秩
  repeat {
    # 一直执行就是没出错
    A = randi(20, 10, 9)
    res = f(A)
    Q = res$Q
    R = res$R
    n = ncol(A)
    r_res = qr(A)
    # browser()
    stopifnot((A %==% (Q %*% R)), ((t(Q) %*% Q) %==% diag(n)), (abs(qr.Q(r_res)) %==% abs(Q)), (abs(qr.R(r_res)) %==% abs(R)))
    print("Right")
  }
}
