A=hilb(10)
res.lu=LU(A)
b=randi(10,10,1)
.t=forwardsolve(res.lu$L,res.lu$P%*%b,upper.tri=FALSE)
x.lu=forwardsolve(res.lu$U,.t,upper.tri = TRUE)

res.qr=QR(A)
Q=res.qr$Q
R=res.qr$R
.b=t(Q)%*%b
x.qr = forwardsolve(R, .b, upper.tri = TRUE)

cond(b-A%*%x.lu)
cond(b - A %*% x.qr)

