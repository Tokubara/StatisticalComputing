library(pracma)

gmm_norm_3 <- function(rv.norm, mu.t = c(0, 0, 0), sigma2.t = c(1, 1, 1), p.t = c(1, 1, 1) / 3, max_iter = 5000, tol = 1e-6) {
  # 之后计算需要用到的量
  len_rv = length(rv.norm)
  sigma.t = sqrt(sigma2.t)
  iter_counter = 0
  while (iter_counter < max_iter) {
    # 计算条件概率
    mu.t.old=mu.t
    lij.matrix = sapply(1:3, function(i) { dnorm(rv.norm, mu.t[i], sigma.t[i]) *p.t[i] }) # m表示是个矩阵
    lij.matrix = lij.matrix/rowSums(lij.matrix)
    stopifnot("underflow" = !any(is.na(lij.matrix)), "underflow" = !any(is.nan(lij.matrix)))
    lij.colsum = colSums(lij.matrix)
    # mu.t.old = mu.t # 保存用于判停
    # 更新μ
    mu.t = apply(lij.matrix, MARGIN = 2, function(col) dot(col, rv.norm) / sum(col))
    # 更新σ2
    for (j_percol in 1:3) {
      sigma2.t[j_percol] = dot(lij.matrix[, j_percol], (rv.norm - mu.t[j_percol]) ^ 2) / lij.colsum[j_percol]
    }
    sigma.t = sqrt(sigma2.t)
    # 更新p
    p.t = lij.colsum / len_rv
    if(abs(sum(p.t)-1)>=tol) {
      print(p.t)
      stop("sum of p!=1")
    }
    # 判停
    iter_counter = iter_counter + 1
    print(mu.t)
    # print(iter_counter)
    # print(sum((mu.t.old - mu.t) ^ 2))
    # if (sqrt(sum((mu.t.old - mu.t) ^ 2)) / sqrt(sum(mu.t.old^2))<tol) {
    #   break
    # }
  }
  return(list(mu = mu.t, sigma2 = sigma2.t, sigma = sigma.t, p = p.t, iter_num = iter_counter))
}

# 设定真实参数
len_rv = 2000
mu.preset = c(3, 1, 7)
sigma2.preset = c(2, 1, 3)
sigma.preset = sqrt(sigma2.preset)
p.preset = c(0.2, 0.3, 0.5)

# 生成数据点
set.seed(17)
rv.latent = sample(x = 1:3, size = len_rv, prob = p.preset, replace = TRUE)
rv.norm = rnorm(len_rv, mu.preset[rv.latent], sigma.preset[rv.latent])

res = gmm_norm_3(rv.norm,c(2,2,2))
