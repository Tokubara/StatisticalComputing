library(pracma)

cut_n_parts <- function(x, part_num) {
  # 表示将x从小到大等分为3部分, 返回每一部分的均值
  stopifnot(is.null(dim(x)))
  # 排序
  x = sort(x)
  len_x = length(x)
  per_part_len <- floor(len_x / part_num)
  # 计算每一段的index
  index <- seq(from = 1, len = part_num, by = per_part_len)
  index[part_num + 1] = len_x + 1
  # stopifnot(all(start_index<=end_index))
  # 得到均值
  per_part <- lapply(1:part_num, function(i) {(x[index[i]:(index[i + 1] - 1)]) })
  return(per_part)
}
# cut_n_parts(1:10, 3)
# cut_n_parts(1:2000, 4)


gmm_norm <- function(rv.norm, n, mu.t = NULL, sigma.t = NULL, p.t = NULL, max_iter = 5000, tol = 1e-6) {
  cut_rv.norm = cut_n_parts(rv.norm, n)
  if (is.null(mu.t)) {
    mu.t = sapply(cut_rv.norm, mean)
  }
  if (is.null(sigma.t)) {
    sigma.t = sapply(cut_rv.norm, sd)
  }
  p.t = rep(1 / n, n)
  # 之后计算需要用到的量
  len_rv = length(rv.norm)
  sigma2.t = sigma.t ^ 2
  iter_counter = 0
  while (iter_counter < max_iter) {
    # 计算条件概率
    mu.t.old = mu.t
    lij.matrix = sapply(1:n, function(i) { dnorm(rv.norm, mu.t[i], sigma.t[i]) * p.t[i] }) # m表示是个矩阵
    lij.matrix = lij.matrix / rowSums(lij.matrix)
    stopifnot("underflow" = !any(is.na(lij.matrix)), "underflow" = !any(is.nan(lij.matrix)))
    lij.colsum = colSums(lij.matrix)
    # mu.t.old = mu.t # 保存用于判停
    # 更新μ
    mu.t = apply(lij.matrix, MARGIN = 2, function(col) dot(col, rv.norm) / sum(col))
    # 更新σ2
    for (j_percol in 1:n) {
      sigma2.t[j_percol] = dot(lij.matrix[, j_percol], (rv.norm - mu.t[j_percol]) ^ 2) / lij.colsum[j_percol]
    }
    sigma.t = sqrt(sigma2.t)
    # 更新p
    p.t = lij.colsum / len_rv
    stopifnot("sum of p!=1" = abs(sum(p.t) - 1) <= tol)
    # 判停
    iter_counter = iter_counter + 1
    if (sqrt(sum((mu.t.old - mu.t) ^ 2)) / sqrt(sum(mu.t.old ^ 2)) < tol) {
      break
    }
  }
  return(list(mu = mu.t, sigma2 = sigma2.t, sigma = sigma.t, p = p.t, iter_num = iter_counter))
}


# 设定真实参数
test.gmm_norm <- function(len_rv, mu.preset, sigma2.preset, p.preset, max_iter=5000) {
  # 生成数据点
  stopifnot("sum of p not equals 1" = near(sum(p.preset), 1), length(mu.preset) == length(sigma2.preset), length(mu.preset) == length(p.preset))
  set.seed(17)
  n = length(mu.preset)
  rv.latent = sample(x = 1:n, size = len_rv, prob = p.preset, replace = TRUE)
  sigma.preset = sqrt(sigma2.preset)
  rv.norm = rnorm(len_rv, mu.preset[rv.latent], sigma.preset)
  res = gmm_norm(rv.norm, n = n,max_iter = max_iter)
  return(res)
}

res = test.gmm_norm(10000, c(-1, 2, 4), c(1, 1, 1), c(0.2, 0.3, 0.5), max_iter=10000)

