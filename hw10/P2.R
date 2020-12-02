library(pracma)
# 设定真实参数
len_rv = 2000
mu.preset = c(-20, 0, 30)
sigma2.preset = c(6, 1, 12)
sigma.preset = sqrt(sigma2.preset)
p.preset = c(0.2, 0.3, 0.5)

# 生成数据点
set.seed(17)
rv.latent = sample(x = 1:3, size = len_rv, prob = p.preset, replace = TRUE)
rv.norm = rnorm(len_rv, mu.preset[rv.latent], sigma.preset[rv.latent])

# 设定初始参数
mu.t = c(-10, 0, 10)
sigma2.t = c(6, 3, 9)
sigma.t = sqrt(sigma2.t)
p.t = c(1, 1, 1) / 3

# 进行迭代
for (k_iter in 1:500) {
  lij.matrix = cbind(dnorm(rv.norm, mu.t[1], sigma.t[1]) * p.t[1], dnorm(rv.norm, mu.t[2], sigma.t[2]) * p.t[2], dnorm(rv.norm, mu.t[3], sigma.t[3]) * p.t[3]) # m表示是个矩阵
  lij.matrix = t(apply(lij.matrix, MARGIN = 1, function(x) x / sum(x)))
  lij.colsum = colSums(lij.matrix)
  mu.t = apply(lij.matrix, MARGIN = 2, function(x) dot(x, rv.norm) / sum(x))
  for (j_percol in 1:3) {
    sigma2.t[j_percol] = dot(lij.matrix[, j_percol], (rv.norm - mu.t[j_percol]) ^ 2) / lij.colsum[j_percol]
  }
  sigma.t = sqrt(sigma2.t)
  p.t = lij.colsum / len_rv
}
