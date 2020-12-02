library(dplyr)
library(pracma)
# 设定真实参数
n=2000
mu=c(-20,0,30)
sigma2=c(6,1,12)
sigma=sqrt(sigma2)
p=c(0.2,0.3,0.5)

# 生成数据点
set.seed(17)
ry=sample(x=1:3,size=n,prob=p,replace = TRUE)
rx=rnorm(n,mu[ry],sigma[ry])

# 设定初始参数
mu_=c(-10,0,10)
sigma2_=c(6,3,9)
sigma_=sqrt(sigma2_)
p_=c(1,1,1)/3

# 进行迭代
for(k in 1:500) {
  am = cbind(dnorm(rx, mu_[1], sigma_[1]) * p_[1], dnorm(rx, mu_[2], sigma_[2]) * p_[2], dnorm(rx, mu_[3], sigma_[3]) * p_[3]) # m表示是个矩阵
  am=t(apply(am,MARGIN=1,function(x)x/sum(x)))
  am_sum = colSums(am)
  mu_=apply(am,MARGIN=2,function(x)dot(x,rx)/sum(x))
  for(j in 1:3) {
    sigma2_[j] = dot(am[, j], (rx - mu_[j]) ^ 2) / am_sum[j]
  }
  sigma_=sqrt(sigma2_)
  p_=am_sum/n
}
