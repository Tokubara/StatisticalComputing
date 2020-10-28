library(fpCompare)
A = randi(10, 3, 3)
old_A=A

m = nrow(A)
n = ncol(A)

if (n > m) {
  # A满足行数>=列数
  stop("n>m, stop")
}
Q=diag(m)
for(i in 1:n) {
  qi=A[i:m,i]
  .t = norm(qi, "2")
  if((.t)%==%(abs(qi[1]))) {
    next
  }
  qi[1] = qi[1] + sign(qi[1]) * .t
  qi=qi/norm(qi,"2") # 单位化, 但好像也未必需要
  A[i, i] = sign(qi[1]) * .t
  A[(i+1):m,i]=0
  for(j in (i+1):n) {
    # 对其它列作用
    aj = A[i:m, j]
    if(j>=n) break
    .t=dot(qi,aj)
    aj=aj-2*.t*qi
    A[i:m, j]=aj
  }
  H=diag(m)
  H[i:m,i:m]=diag(m-i+1)-2*qi%*%t(qi)
  Q=Q%*%t(H)
}
R=A
old_A-Q%*%R