# 此函数在其它函数之前执行, 用于载入包和运算符
library(pracma)
library(fpCompare)
`:` <- function(a, b) {
  if (b < a) {
    return(numeric(0))
  } else {
    return(seq(a, b))
  }
}
source("GS.R")
source("Givens.R")
source("Householder.R")
source("QR.R")
source("LU.R")