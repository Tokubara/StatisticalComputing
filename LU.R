library(pracma)
library(fpCompare)
LU <- function(A, pivoting = TRUE) {
  .f <- function() {
    # 这个函数的作用是返回, 已知p,i的情况下, 去掉p的那些行, 看第i列哪一行的元素最大, 返回这一行
    # 为什么需要多写这么个函数, 因为R关于-在矩阵index中不一致的表现, p=c(0,0,0),那么a[-p]得到的为空
    if (i == 1) {
      # 用不了ifelse, 因为很遗憾, ifelse会去掉names属性
      .t = A[, i]
    } else {
      .t = A[-p, i]
    }
    return(as.integer(names(which.max(abs(.t)))))
  }
  # 下面进行循环前的初始化
  n = nrow(A)
  if (n != ncol(A)) {
    stop("A is not n*n") # 必须对方阵进行
  }
  N_C = 1:n # 大小表明这是个常向量
  p = integer(n)
  # 这里A没有修改
  # L = matrix(0,n,n)
  # A = A
  rownames(A) = 1:n
  for (i in 1:(n - 1)) {
    # 因为一共要进行n-1次
    p[i] = ifelse(pivoting, .f(), i)
    aii = A[p[i], i] # 就是这个对角元
    if (aii == 0) {
      # 如果是部分选主元, 这种情况很难发生, 这主要是写给不选主元的
      stop("can't be decomposed")
    }
    qi = A[p[i], (i + 1):n] # perf:qi = A[p[i], i:n]:因为对角元不参与运算 # fix:这里应该是U而不是A
    for (j in N_C[-p]) {
      A[j, i] = A[j, i] / aii
      A[j, (i + 1):n] = A[j, (i + 1):n] - A[j, i] * qi
    }
  }
  p[n] = N_C[-p] # 补全p, 以便得到L,U,P
  P = matrix(0, n, n)
  P[cbind(1:n, p)] = 1
  U=matrix(0,n,n)
  L=diag(n)
  for(i in 1:n) {
    U[i, i:n] = A[p[i], i:n]
    L[i, 1:(i - 1)] = A[p[i], 1:(i - 1)]
  }
  return(list(L = L, U = U, P = P))
}
comment(LU) <-"如果pivoting为TRUE, 部分选主元, 否则全选主元, 此时P返回为I"

test_LU<-function(A=NULL) {
  # 检验
  if(is.null(A)) {
    A=randi(20,10)
  }
  res=LU(A)
  L=res$L
  U=res$U
  P=res$P
  print(all((P%*%A)%==%(L%*%U)))
}
