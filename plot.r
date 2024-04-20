#nolint start
library(MASS)

n <- 100
p <- 300
set.seed(1234)

corr_x <- matrix(rep(0,p*p),p,p)
for (i in 1:p) {
  for (j in 1:p) {
    corr_x[i,j] <- 0.5^(abs(i-j))
  }
}
x <- mvrnorm(n, numeric(p), corr_x)
x[,1] <- 1
beta <- rep(0,p)
beta[2:4] <- 2
beta[1] <- 1
y <- x %*% beta + rnorm(n, 0, sqrt(3))

index <- c(1:p)
BIC_arr <- c()  #BIC array of all possible combination
min_BIC <- Inf  #global minimal BIC
M <- c()        #set of index that satisfies the inequality conditon
BIC_M <- 0      #BIC of M

EBIC_arr <- c()
min_EBIC <- Inf
E_M <- c()
EBIC_M <- 0

for(i in 1:p){
  combi <- combn(index,i) 
  print(i)
  for(j in 1:ncol(combi)){
    temp <- combi[,j]
    fit <- lm(y ~ x[,temp] - 1)
    BIC_arr <- c(BIC_arr, BIC(fit))
    EBIC_arr <- c(EBIC_arr, (BIC_arr[length(BIC_arr)]) + 2 * lchoose(p, length(temp)))
    if(min_BIC > min(BIC_arr)){
      min_BIC <- min(BIC_arr)
      M_hat <- temp
    }
    if(min_EBIC > min(EBIC_arr)){
      min_EBIC <- min(EBIC_arr)
      EBIC_M_hat <- temp
    }
  }
}
min_M <- M_hat
EBIC_min_M <- EBIC_M_hat

L <- length(M_hat) - 1
for(i in 1:L){
  combi <- combn(M_hat,i)
  for(j in 1:ncol(combi)){
    temp <- combi[,j]
    fit <- lm(y ~ x[,temp] - 1)
    BIC_temp <- BIC(fit)
    #print(temp)
    #print(BIC_temp)
    #print(min_BIC - BIC_temp)
    if(abs(min_BIC - BIC_temp) < 2 * ((length(M_hat) - length(temp)) * log(p))){
      M <- temp
      BIC_M <- BIC_temp
    if(length(M) < length(min_M)){
      min_M <- M
      min_BIC_M <- BIC_M
    }
    }
  }
}

beta_hat <- function(col){
    return(((t(x[,col]) %*% x[,col])^(-1)) %*% t(x[,col]) %*% y)
  }
  sigma_hat <- function(col){
    return(((t(y - x[,col] %*% beta) %*% (y - x[,col]))^2) / n)
  }
  vector_norm <- function(col, hat_beta){
    return(as.numeric(t(y - x[,col] %*% hat_beta) %*% (y - x[,col] %*% hat_beta)))
  }
  likelihood <- function(sigma, col, hat_sigma, vector){
    return(log(1 / (2 * pi * hat_sigma)^(n/2) * exp(-1 / (2 * hat_sigma) * vector)))
  }
  func_BIC <- function(likelihood, col){
    return(-2 * likelihood + ncol(col) * log(n))
  }

L<-length(BIC_arr)
true_BIC <- BIC(lm(y ~ x[,c(1,2,3,4)] - 1))
plot((1:L)[true_condition==1],BIC_arr[true_condition==1],pch=3,col=2,ylim=c(min(BIC_arr),max(BIC_arr)),xlim=c(0,L))
points((1:L)[true_condition!=1],BIC_arr[true_condition!=1],pch=4,col=4)
abline(h=BIC_arr[which(true_model==1)])
points(1,true_BIC,pch = 7, col = 19)
#nolint end