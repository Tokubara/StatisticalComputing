m = 3
A = randi(10, m, m)
B = matrix(0, 2 * m, 2 * m)
B[1:m, (m + 1):(2 * m)] = t(A)
B[(m + 1):(2 * m), 1:m] = A
svd_A = svd(A) # 13.606468  5.698691  2.321411
eig_B = eigen(B) # 13.606468   5.698691   2.321411  -2.321411  -5.698691 -13.606468
# 正数的验证
.x = c(svd_A$v[, 1], svd_A$u[, 1])
B %*% .x - svd_A$d[1] * .x
# 负数的验证
.x = c(svd_A$v[, 1], - svd_A$u[, 1])
B %*% .x + svd_A$d[1] * .x
X = rbind(cbind(svd_A$v, svd_A$v), cbind(svd_A$u, -svd_A$u))
L = diag(c(svd_A$d, -svd_A$d))
B-X%*%L%*%solve(X)