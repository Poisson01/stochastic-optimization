#nolint start
library(MASS)
#define data
n <- 150
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
  beta[1:7] <- 1.5
  y <- 1 + x %*% beta + rnorm(n, 0, sqrt(3))
  
#Greedy algorithm
  greedy_start <- Sys.time()
  init <- sample(2:p,1)
  s <- c(1,init)
  fit <- lm(y ~ x[,s] - 1)
  EBIC_s <- BIC(fit) + 2 * g * lchoose(p, length(s))
  roop <- TRUE
  while (roop == TRUE) {
    min_EBIC <- EBIC_s
    roop <- FALSE
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
    if(EBIC_s > min_EBIC){
      EBIC_s <- min_EBIC
      s <- st_plus
      roop <- TRUE
    } 
    if(length(s) > n * 0.6){
      roop <- FALSE
    }
  }
  greedy_end <- Sys.time()
  time_greedy_arr[temp] <- (greedy_end - greedy_start)
  greedy_s <- s
  greedy_EBIC <- EBIC_s
  EBIC_arr[temp, 1] <- greedy_EBIC
  bool_vector <- rep(FALSE, p)
  bool_vector[greedy_s] <- TRUE
  greedy_index_arr[temp,] <- bool_vector
  
#stochastic algorithm
  stochastic_start <- Sys.time()
  M <- 100
  all_index <- c(1:p)
  init <- sample(2:p, 1)
  s <- c(1,init); s_min <- s
  fit <- lm(y ~ x[, s] - 1)
  s_EBIC_s <- BIC(fit) + 2 * g * lchoose(p, length(s))
  EBIC_s_min <- s_EBIC_s
  for (iter in 1:M) {
    for(i in 2:p){
      plus_s <- rep(FALSE,p); minus_s <- rep(FALSE,p)
      plus_s[s] <- TRUE; minus_s[s] <- TRUE
      if(plus_s[i]==TRUE){
        EBIC_plus <- s_EBIC_s
        minus_s[i] <- FALSE
        minus_fit <- lm(y ~ x[,minus_s] - 1)
        EBIC_minus <- BIC(minus_fit) + 2 * g * lchoose(p, sum(minus_s))
      }else{
        EBIC_minus <- s_EBIC_s
        plus_s[i] <- TRUE
        plus_fit <- lm(y ~ x[,plus_s] - 1)
        EBIC_plus <- BIC(plus_fit) + 2 * g * lchoose(p, sum(plus_s))
      }
      #minus_s[i] <- FALSE
      #minus_fit <- lm(y ~ x[,minus_s] - 1)
      #EBIC_minus <- BIC(minus_fit) + 2 * lchoose(p, sum(minus_s))
      prob <- 1/(1 + exp((-0.5 * EBIC_minus) + (0.5 * EBIC_plus)))
      z <- rbinom(1,1,prob)
      if((sum(plus_s) > 0.6 * n) & (z == 1)){
      
      }else{
      if (z == 1) {
        s <- all_index[plus_s]
        s_EBIC_s <- EBIC_plus
      }else{
        s <- all_index[minus_s]
        s_EBIC_s <- EBIC_minus
      }
      #fit <- lm(y ~ x[,s] - 1)
      #s_EBIC_s <- BIC(fit) + 2 * lchoose(p, length(s))
      if(EBIC_s_min > s_EBIC_s){
        EBIC_s_min <- s_EBIC_s
        s_min <- s
        }
      }
    }
    #cat("---------",iter,"---------\n")
    #print("s")
    #print(s)
    #print("local_EBIC")
    #print(s_EBIC_s)
    #print("s_min")
    #print(s_min)
    #print("global_EBIC")
    #print(EBIC_s_min)
  }
  stochastic_end <- Sys.time()
  time_stochastic_arr[temp] <- (stochastic_end - stochastic_start)
  stochastic_s <- s_min
  stochastic_EBIC <- EBIC_s_min
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
    cat("greedy performed better ",i,"th iteration\n")
    cat("greedy: ",EBIC_arr[i,1],"stochastic: ",EBIC_arr[i,2],"\n")
    cat("greedy: ",all_index[greedy_index_arr[i,]],"stochastic: ",all_index[stochastic_index_arr[i,]],"\n")
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
EBIC_150_arr <- EBIC_arr
greedy_150_index_arr <- greedy_index_arr
stochastic_150_index_arr <- stochastic_index_arr

#nolint end