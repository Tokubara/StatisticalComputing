library(pracma)
library(fpCompare)
LU <- function(A, pivoting = TRUE) {
  .f <- function() {
    # 这个函数的作用是返回, 已知p,i的情况下, 去掉p的那些行,  看第i列哪一行的元素最大, 返回这一行
    if (i==1) { # 用不了ifelse, 因为很遗憾, ifelse会去掉names属性
      .t=U[, i]
    } else {
      .t=U[-p, i]
    }
    return(as.integer(names(which.max(abs(.t)))))
  }
  n = nrow(A)
  if (n != ncol(A)) {
    stop("A is not n*n") # 必须对方阵进行
  }
  N_C = 1:n # 大小表明这是个常向量
  p = integer(n)
  # 这里A没有修改
  L = matrix(0,n,n)
  U = A
  rownames(U) = 1:n
  for (i in 1:(n - 1)) {
    # 因为一共要进行n-1次
    if (pivoting) {
      p[i] = .f() # 找出第i列绝对值最大者 # fix:p[i]=which.max(abs(A[,i])) A没有改变, 应该是U # ifelse(all(p == 0), U[, i], U[-p, i])不行, 因为它会去掉names
    } else {
      p[i] = i
    }
    aii = U[p[i],i] # 就是这个对角元
    if (aii == 0) {
      stop("can't be decomposed")
    }
    qi = U[p[i], (i + 1):n] # perf:qi = A[p[i], i:n]:因为对角元不参与运算 # fix:这里应该是U而不是A
    for (j in N_C[-p]) {
      qij = U[j, i] / aii
      U[j, (i + 1):n] = U[j, (i + 1):n] - qij * qi
      L[j, i] = qij
    }
  }
  p[n] = N_C[-p] # 补全p, 以便得到L,U,P
  # 得到U
  U = U[p,]
  U[lower.tri(U)] <- 0
  # 得到L
  L = L[p,]
  diag(L) <- 1
  L[upper.tri(L)] <- 0
  # 得到P
  P = matrix(0, n, n)
  P[cbind(1:n, p)] = 1
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

