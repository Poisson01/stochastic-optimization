#nolint start
library(MASS)
#define data
n <- 100
p <- 1000
g <- 1
times <- 100
win_greedy <- 0
win_stochastic <- 0
win_both <- 0
EBIC_arr <- matrix(0, times, 2)
greedy_index_arr <- matrix(FALSE, times, p)
stochastic_index_arr <- matrix(FALSE, times, p)
time_greedy_arr <- rep(0, times)
time_stochastic_arr <- rep(0, times)

for(temp in 1:times){
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
  beta[1] <-1
  y <- x %*% beta + rnorm(n, 0, sqrt(3))
  assign(paste0("x",temp),x)
  assign(paste0("y",temp),y)
  
#Greedy algorithm
  greedy_start <- Sys.time()
  init <- sample(2:p,1)
  s <- c(1,init)
  fit <- lm(y ~ x[,s] - 1)
  EBIC_g <- BIC(fit) + 2 * g * lchoose(p, length(s))
  loop <- TRUE
  while (loop == TRUE) {
    min_EBIC <- EBIC_g
    loop <- FALSE
    for(i in s){
      st <- sort(union(1,setdiff(s, i)))
      fit <- lm(y ~ x[,st] - 1)
      EBIC <- BIC(fit) + 2 * g * lchoose(p, length(st))
      if (min_EBIC > EBIC) {
        min_EBIC <- EBIC
        st_plus <- st
      }
    }
    for(i in setdiff(c(2:p), s)){
      st <- sort(c(s, i))
      fit <- lm(y ~ x[,st] - 1)
      EBIC <- BIC(fit) + 2 * g * lchoose(p, length(st))
      if (min_EBIC > EBIC){
        min_EBIC <- EBIC
        st_plus <- st
      }
    }
    if(EBIC_g > min_EBIC){
      EBIC_g <- min_EBIC
      s <- st_plus
      loop <- TRUE
    } 
    if(length(s) > n * 0.5){
      loop <- FALSE
    }
  }
  greedy_end <- Sys.time()
  time_greedy_arr[temp] <- (greedy_end - greedy_start)
  greedy_s <- s
  greedy_EBIC <- EBIC_g
  EBIC_arr[temp, 1] <- greedy_EBIC
  bool_vector <- rep(FALSE, p)
  bool_vector[greedy_s] <- TRUE
  greedy_index_arr[temp,] <- bool_vector
  
#stochastic algorithm
  stochastic_start <- Sys.time()
  M <- 300
  all_index <- c(1:p)
  init <- sample(2:p, 1)
  s <- c(1,init); s_min <- s
  fit <- lm(y ~ x[, s] - 1)
  EBIC_ss <- BIC(fit) + 2 * g * lchoose(p, length(s))
  EBIC_s <- EBIC_ss
  for (iter in 1:M) {
    for(i in 2:p){
      plus_s <- rep(FALSE,p); minus_s <- rep(FALSE,p)
      plus_s[s] <- TRUE; minus_s[s] <- TRUE
      if(plus_s[i]==TRUE){
        EBIC_plus <- EBIC_ss
        minus_s[i] <- FALSE
        minus_fit <- lm(y ~ x[,minus_s] - 1)
        EBIC_minus <- BIC(minus_fit) + 2 * g * lchoose(p, sum(minus_s))
      }else{
        EBIC_minus <- EBIC_ss
        plus_s[i] <- TRUE
        plus_fit <- lm(y ~ x[,plus_s] - 1)
        EBIC_plus <- BIC(plus_fit) + 2 * g * lchoose(p, sum(plus_s))
      }
      if(EBIC_s > EBIC_plus ){
        EBIC_s <- EBIC_plus 
        s_min <- all_index[plus_s]
      }
      if(EBIC_s > EBIC_minus ){
        EBIC_s <- EBIC_minus 
        s_min <- all_index[minus_s]
      }
      
      prob <- 1/(1 + exp((-0.5 * EBIC_minus) + (0.5 * EBIC_plus)))
      z <- rbinom(1,1,prob)
      if((sum(plus_s) > 0.5 * n) & (z == 1)){
      
      }else{
      if (z == 1) {
        s <- all_index[plus_s]
        EBIC_ss <- EBIC_plus
      }else{
        s <- all_index[minus_s]
        EBIC_ss <- EBIC_minus
      }
      }
    }
  }
  stochastic_end <- Sys.time()
  time_stochastic_arr[temp] <- (stochastic_end - stochastic_start)
  stochastic_s <- s_min
  stochastic_EBIC <- EBIC_s
  EBIC_arr[temp,2] <- stochastic_EBIC
  bool_vector <- rep(FALSE, p)
  bool_vector[stochastic_s] <- TRUE
  stochastic_index_arr[temp,] <- bool_vector
  print(temp)
}

 for(i in 1:times){
  if(EBIC_arr[i,1] < EBIC_arr[i,2]){
    win_greedy <- win_greedy + 1
  }else if (EBIC_arr[i,1] > EBIC_arr[i,2]) {
     win_stochastic <- win_stochastic + 1
  } else {
    win_both <- win_both + 1
  }
}
cat("greedy performed better", win_greedy, "times.")
cat("stochastic performed better", win_stochastic, "times.")
cat("same performed", win_both, "times.")
cat("greedy average running time: ", mean(time_greedy_arr))
cat("stochastic average running time: ", mean(time_stochastic_arr))
for(i in 1:times){
  if(EBIC_arr[i,1] < EBIC_arr[i,2]){
    cat("greedy performed better",i,"th iteration\n")
    cat("greedy: ",EBIC_arr[i,1],"/ stochastic: ",EBIC_arr[i,2],"\n")
    cat("greedy: ",all_index[greedy_index_arr[i,]],"/ stochastic: ",all_index[stochastic_index_arr[i,]],"\n")
    cat("\n")
  }
}
for(i in 1:times){
  cat("greedy", i, "th select model: ", all_index[greedy_index_arr[i,]],"\n")
  cat("stochastic", i, "th select model: ", all_index[stochastic_index_arr[i,]],"\n")
  cat(i,"th EBIC: ",EBIC_arr[i,],"\n")
  if(EBIC_arr[i,1] < EBIC_arr[i,2]){
    cat("greedy win\n")
  }else if(EBIC_arr[i,1] > EBIC_arr[i,2]){
    cat("stochastic win\n")
  }
  cat("\n")
}
EBIC_200_arr <- EBIC_arr
greedy_200_index_arr <- greedy_index_arr
stochastic_200_index_arr <- stochastic_index_arr

EBIC_arr <- EBIC_100_arr
greedy_index_arr <- greedy_100_index_arr
stochastic_index_arr <- stochastic_100_index_arr

win_greedy <- 0 
win_stochastic <- 0
win_both <- 0

#nolint end