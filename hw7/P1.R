d <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44, 3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75, 0.27, 43.21)

l_f<-function(x,theta) {
  sum(-log(pi*(1+(x-theta)^2)))
}

theta = seq(-1, 1, by = 0.01)
l_c = sapply(theta, l_f,d)
plot(x=theta,y=l_c,pch=19)

# 二分法
l=-1
r=1
mid=(l+r)/2

d_l_f<-function(l_f,x,dx=0.01) { # ppt上设的是0.01, 太小了确实有抵消的问题
  (l_f(d,x+dx)-l_f(d,x-dx))/(2*dx)
}

f_mid = d_l_f(, mid)

while()