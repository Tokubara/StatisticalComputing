
QR <- function(A, mode = c("Househoulder","MGS", "CGS", "Givens")) {
  mode=match.arg(mode)
  f=switch(mode,
    MGS=MGS,
    CGS=CGS,
    Househoulder = householder,
    Givens=Givens
  )
  return(f(A))
}

QR.test <- function(A = NULL,mode="C") {
  if (is.null(A)) {
    A = randi(20, 10, 8)
  }
  # browser()
  res = QR(A,mode=mode)
  Q = res$Q
  R = res$R
  print(all((Q %*% R) %==% (A)))
}
