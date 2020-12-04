A = rbind(c(0, 2), c(3, 1))
res = LU(A) # 选主元的情况
b = res$P %*% c(4, 4)
.t = forwardsolve(res$L, b, upper.tri = FALSE)
x.p = forwardsolve(res$U, .t, upper.tri = TRUE)

res = LU(A, FALSE) # 不选主元的情况

