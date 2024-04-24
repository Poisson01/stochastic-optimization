#nolint start
library(MASS)
n <- 100 #50, 100, 150, 200, ...
p <- 1000
g <- 1
m <- 1000
times <- 1
TNR_arr <- matrix(0,times,3)
TPR_arr <- matrix(0,times,3)
model_selected_arr <- rep(0,times)
for(aa in 1:times){
  TNR_denominator <- p
  TPR_denominator <- 0
  set.seed(1000 + aa)
  corr_x <- matrix(rep(0,p*p),p,p)
  for (i in 1:p) {
    for (j in 1:p) {
      corr_x[i,j] <- 0.5^(abs(i-j))
    }
  }
  x <- mvrnorm(n, numeric(p), corr_x)
  x[,1] <- 1
  beta <- rep(0,p)
  beta[2:7] <- 2
  beta[1] <- 1
  y <- x %*% beta + rnorm(n, 0, sqrt(3))
  for(i in beta){
    if(i){
      TPR_denominator <- TPR_denominator + 1
      TNR_denominator <- TNR_denominator - 1
    }
  }
  TRUE_values <- which(beta != 0)
  FALSE_values <- which(beta == 0)
  #EBIC
  all_index <- c(1:p)
  init <- sample(2:p, 1)
  s <- c(1); s_min <- s
  fit <- lm(y ~ x[, s] - 1)
  s_EBIC_s <- BIC(fit) + 2 * lchoose(p, length(s) + 1)
  EBIC_s_min <- s_EBIC_s
  for (iter in 1:m) {
    for(i in 2:p){
      plus_s <- rep(FALSE,p); minus_s <- rep(FALSE,p)
      plus_s[s] <- TRUE; minus_s[s] <- TRUE
      if(plus_s[i]==TRUE){
        EBIC_plus <- s_EBIC_s
        minus_s[i] <- FALSE
        minus_fit <- lm(y ~ x[,minus_s] - 1)
        EBIC_minus <- BIC(minus_fit) + 2 * lchoose(p, sum(minus_s) + 1)
      }else{
        EBIC_minus <- s_EBIC_s
        plus_s[i] <- TRUE
        plus_fit <- lm(y ~ x[,plus_s] - 1)
        EBIC_plus <- BIC(plus_fit) + 2 * lchoose(p, sum(plus_s) + 1)
      }
      if(sum(plus_s) <= n^(2/3)){
        if(EBIC_s_min > EBIC_plus){
          EBIC_s_min <- EBIC_plus
          s_min <- all_index[plus_s]
        }
      }
      if(EBIC_s_min > EBIC_minus){
        EBIC_s_min <- EBIC_minus
        s_min <- all_index[minus_s]
      }
      prob <- 1/(1 + exp((-0.5 * EBIC_minus) + (0.5 * EBIC_plus)))
      z <- rbinom(1,1,prob)
      if((sum(plus_s) > n^(2/3)) & (z == 1)){
        
      }else{
        if (z == 1) {
          s <- all_index[plus_s]
          s_EBIC_s <- EBIC_plus
        }else{
          s <- all_index[minus_s]
          s_EBIC_s <- EBIC_minus
        }
      }
      #print(s)
    }
    #if(iter %% 100 == 0){
      #cat("EBIC",iter,"\n")
    #}
  }
  EBIC_TPR_numerator <- sum(s_min %in% TRUE_values)
  EBIC_TNR_numerator <- TNR_denominator - sum(s_min %in% FALSE_values)
  TPR_arr[aa,1] <- EBIC_TPR_numerator / TPR_denominator
  TNR_arr[aa,1] <- EBIC_TNR_numerator / TNR_denominator
  
  #BIC 1st
  s <- c(1); M_hat <- s
  fit <- lm(y ~ x[, s] - 1)
  s_BIC_s <- BIC(fit)
  M_hat_BIC <- s_BIC_s
  for (iter in 1:m) {
    for(i in 2:p){
      plus_s <- rep(FALSE,p); minus_s <- rep(FALSE,p)
      plus_s[s] <- TRUE; minus_s[s] <- TRUE
      if(plus_s[i]==TRUE){
        BIC_plus <- s_BIC_s
        minus_s[i] <- FALSE
        minus_fit <- lm(y ~ x[,minus_s] - 1)
        BIC_minus <- BIC(minus_fit)
      }else{
        BIC_minus <- s_BIC_s
        plus_s[i] <- TRUE
        plus_fit <- lm(y ~ x[,plus_s] - 1)
        BIC_plus <- BIC(plus_fit)
      }
      if(sum(plus_s) <= n^(2/3)){
        if(M_hat_BIC > BIC_plus){
          M_hat_BIC <- BIC_plus
          M_hat <- all_index[plus_s]
        }
      }
      if(M_hat_BIC > BIC_minus){
        M_hat_BIC <- BIC_minus
        M_hat <- all_index[minus_s]
      }
      prob <- 1/(1 + exp((-0.5 * BIC_minus) + (0.5 * BIC_plus)))
      z <- rbinom(1,1,prob)
      if((sum(plus_s) > n^(2/3)) & (z == 1)){
        
      }else{
        if (z == 1) {
          s <- all_index[plus_s]
          s_BIC_s <- BIC_plus
        }else{
          s <- all_index[minus_s]
          s_BIC_s <- BIC_minus
        }
      }
    }
    #if(iter %% 100 == 0){
      #cat("BIC 1st",iter,"\n")
    #}
  }
  BIC1st_TPR_numerator <- sum(M_hat %in% TRUE_values)
  BIC1st_TNR_numerator <- TNR_denominator - sum(M_hat %in% FALSE_values)
  TPR_arr[aa,2] <- BIC1st_TPR_numerator / TPR_denominator
  TNR_arr[aa,2] <- BIC1st_TNR_numerator / TNR_denominator

  # BIC 2nd in paper
  M_hat <- sort(M_hat); min_M <- M_hat; M <- c(); BIC_M <- 0; temp <- 2
  M_hat_beta <- solve(t(x[,M_hat]) %*% x[,M_hat]) %*% t(x[,M_hat]) %*% y
  M_hat_sigma <- (t(y - x[,M_hat] %*% M_hat_beta) %*% (y - x[,M_hat] %*% M_hat_beta)) / n
  visited_BIC <- c()
  model_selected <- FALSE
  M_hat_BIC <- -2 * sum(dnorm(y, x[,M_hat] %*% M_hat_beta, sqrt(M_hat_sigma), log = TRUE)) #+ (length(M_hat)) * log(n)
  #visited_s <- matrix(0,m,p)
  s <- c(1)
  s_beta_hat <- (solve(t(x[,s]) %*% x[,s])) %*% t(x[,s]) %*% y
  s_sigma_hat <- (t(y - x[,s] %*% s_beta_hat) %*% (y - x[,s] %*% s_beta_hat)) / n
  s_BIC <-  -2 * sum(dnorm(y, x[,s] %*% s_beta_hat, sqrt(s_sigma_hat), log = TRUE)) #+ (length(s)) * log(n)
  s_fit <- lm(y ~ x[,s] - 1)
  s_EBIC <- BIC(s_fit) + 2 * lchoose(p, length(s) + 1)
  for(iter in 1:m){
    for(i in M_hat[-1]){
      plus_s <- rep(FALSE, p); minus_s <- rep(FALSE, p)
      plus_s[s] <- TRUE; minus_s[s] <- TRUE
      if(plus_s[i] == TRUE){
        BIC_plus <- s_BIC
        EBIC_plus <- s_EBIC
        minus_s[i] <- FALSE
        minus_beta_hat <- (solve(t(x[,minus_s]) %*% x[,minus_s])) %*% t(x[,minus_s]) %*% y
        minus_sigma_hat <- (t(y - x[,minus_s] %*% minus_beta_hat) %*% (y - x[,minus_s] %*% minus_beta_hat)) / n
        BIC_minus <-  -2 * sum(dnorm(y, x[,minus_s] %*% minus_beta_hat, sqrt(minus_sigma_hat), log = TRUE)) #+ (sum(minus_s)) * log(n)
        minus_fit <- lm(y ~ x[,minus_s] - 1)
        EBIC_minus <- BIC(minus_fit) + 2 * lchoose(p, sum(minus_s) + 1)
        visited_BIC <- append(visited_BIC,BIC_minus)
      }else {
        BIC_minus <- s_BIC
        EBIC_minus <- s_EBIC
        plus_s[i] <- TRUE
        plus_beta_hat <- (solve(t(x[,plus_s]) %*% x[,plus_s])) %*% t(x[,plus_s]) %*% y
        plus_sigma_hat <- (t(y - x[,plus_s] %*% plus_beta_hat) %*% (y - x[,plus_s] %*% plus_beta_hat)) / n
        BIC_plus <-  -2 * sum(dnorm(y, x[,plus_s] %*% plus_beta_hat, sqrt(plus_sigma_hat), log = TRUE)) #+ (sum(plus_s)) * log(n)
        plus_fit <- lm(y ~ x[,plus_s] - 1)
        EBIC_plus <- BIC(plus_fit) + 2 * lchoose(p, sum(plus_s) + 1)
        visited_BIC <- append(visited_BIC,BIC_plus)
      }
      if(sum(plus_s) <= n^(2/3)){
        p_s <- (length(M_hat) - sum(plus_s))
        plus_cut_off <- n * 2 #2 * p_s * (log(p)) #+ 2 * p_s * log(log(p)) + p_s
        if(abs(M_hat_BIC - BIC_plus) <= plus_cut_off){
          M <- all_index[plus_s]
          BIC_M <- BIC_plus
          if(length(M) < length(min_M)){
            min_M <- M
            min_BIC_M <- BIC_M
          }
        }
      }
      m_s <- (length(M_hat) - sum(minus_s))
      minus_cut_off <- n * 2 #2 * m_s * (log(p)) #+ 2 * m_s * log(log(p)) + m_s
      if(abs(M_hat_BIC - BIC_minus) <= minus_cut_off){
        M <- all_index[minus_s]
        BIC_M <- BIC_minus
        if(length(M) < length(min_M)){
          min_M <- M
          min_BIC_M <- BIC_M
        }
      }
      prob <- 1/(1 + exp((-0.5 / temp * EBIC_minus) + (0.5 / temp * EBIC_plus)))
      z <- rbinom(1, 1, prob)
      if((sum(plus_s) > n ^ (2 / 3)) & (z == 1)){

      }else {
        if(z == 1){
          s <- all_index[plus_s]
          s_EBIC <- EBIC_plus
          s_BIC <- BIC_plus
        }else {
            s <- all_index[minus_s]
            s_EBIC <- EBIC_minus
            s_BIC <- BIC_minus
        }
      }
      #print(s)
      if (identical(as.numeric(s),c(1,2,3,4,5,6,7))) {
        model_selected <- TRUE
      }
    }
    #if(iter %% 100 == 0){
    #  cat("BIC 2nd",iter,"\n")
    #}
  }
  model_selected_arr[aa] <- model_selected
  print(aa)
  #cat("EBIC",s_min,"\n",EBIC_s_min,"\n")
  #cat("BIC 1st",M_hat,"\n",M_hat_BIC,"\n")
  #cat("BIC 2nd",min_M,"\n",min_BIC_M,"\n")
BIC2nd_TPR_numerator <- sum(min_M %in% TRUE_values)
BIC2nd_TNR_numerator <- TNR_denominator - sum(min_M %in% FALSE_values)
TPR_arr[aa,3] <- BIC2nd_TPR_numerator / TPR_denominator
TNR_arr[aa,3] <- BIC2nd_TNR_numerator / TNR_denominator
}

mean(TNR_arr[,1])
mean(TNR_arr[,3])
mean(TPR_arr[,1])
mean(TPR_arr[,3])
model_selected_arr
boxplot(TNR_arr[,1],TNR_arr[,3])
boxplot(TPR_arr[,1],TPR_arr[,3])

for(i in 1:times){
  cat("EBIC", i, "th TPR:", TPR_arr[i,1],"\n")
  cat("BIC first stochatic", i, "th TPR:", TPR_arr[i,2],"\n")
  cat("BIC second stochatic", i, "th TPR:", TPR_arr[i,3],"\n")
  cat("EBIC", i, "th TNR:", TNR_arr[i,1],"\n")
  cat("BIC first stochatic", i, "th TNR:", TNR_arr[i,2],"\n")
  cat("BIC second stochatic", i, "th TNR:", TNR_arr[i,3],"\n")
}
#nolint end