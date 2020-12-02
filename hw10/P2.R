library(dplyr)
library(pracma)
# 设定真实参数
n=2000
mu=c(-1.2,0,1.6)
sigma=c(0.1,0.1,0.1)
p=c(0.2,0.4,0.4)

# 生成数据点
set.seed(17)
ry=sample(1:3,n,replace = TRUE)
rx=rnorm(n,mu[ry],sigma[ry])

# 设定初始参数
mu_=c(-1,0,1)
sigma_=c(0.3,0.3,0.3)
p_=c(1,1,1)/3

# 进行迭代
for(k in 1:5000) {
am = cbind(dnorm(x, mu_[1], sigma_[1]) * p_[1], dnorm(x, mu_[2], sigma_[2]) * p_[2], dnorm(x, mu_[3], sigma_[3]) * p_[3]) # m表示是个矩阵
am=t(apply(am,MARGIN=1,function(x)x/sum(x)))
mu_=apply(am,MARGIN=2,function(x)dot(x,rx)/sum(x))
for(j in 1:3) {
  sigma_[j] = dot(am[, j], (rx - mu_[j]) ^ 2) / sum(am[,j])
}
p_=colSums(am)/n
}
