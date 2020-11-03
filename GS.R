MGS <- function(A) {
  m = nrow(A)
  n = ncol(A)
  if (n > m) {
    # A满足行数>=列数, 这里有问题
    stop("n>m, stop")
  }
  Q = A
  R = matrix(0, n, n)
  for (i in 1:n) {
    qi = Q[, i]
    norm2 = norm(qi, "2")
    if (abs(norm2) %==% 0) {
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
    if (abs(norm2) %==% 0) {
      stop("rank < n, stop") # 如果列不满秩, 报错
    }
    R[i, i] = norm2
    qi = qi / norm2
    Q[, i] = qi
  }
  return(list(Q = Q, R = R))
}