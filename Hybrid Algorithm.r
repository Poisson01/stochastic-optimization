#nolint start
library(MASS)
#define data
n <- 100
p <- 1000
g <- 1
times <- 100
win_greedy <- 0
win_stochastic <- 0
win_hybrid <- 0
all_same <- 0
EBIC_arr <- matrix(0, times, 3)
greedy_index_arr <- matrix(FALSE, times, p)
stochastic_index_arr <- matrix(FALSE, times, p)
hybrid_index_arr <- matrix(FALSE, times, p)
time_greedy_arr <- rep(0, times)
time_stochastic_arr <- rep(0, times)
time_hybrid_arr <- rep(0, times)

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
  beta[1] <- 1
  y <- x %*% beta + rnorm(n, 0, sqrt(3))
  
#Greedy algorithm
  greedy_start <- Sys.time()
  init <- sample(2:p,1)
  s <- c(1)
  fit <- lm(y ~ x[,s] - 1)
  EBIC_s <- BIC(fit) + 2 * g * lchoose(p, length(s))
  loop <- TRUE
  while (loop == TRUE) {
    min_EBIC <- EBIC_s
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
    if(EBIC_s > min_EBIC){
      EBIC_s <- min_EBIC
      s <- st_plus
      loop <- TRUE
    } 
    if(length(s) > n * 0.6){
      loop <- FALSE
    }
  }
  print("greedy end")
  greedy_end <- Sys.time()
  time_greedy_arr[temp] <- difftime(Sys.time(),greedy_start,units = "secs")
  greedy_s <- s
  greedy_EBIC <- EBIC_s
  EBIC_arr[temp, 1] <- greedy_EBIC
  bool_vector <- rep(FALSE, p)
  bool_vector[greedy_s] <- TRUE
  greedy_index_arr[temp,] <- bool_vector
  
#stochastic algorithm
  stochastic_start <- Sys.time()
  M <- 1000
  all_index <- c(1:p)
  init <- sample(2:p, 1)
  s <- c(1); s_min <- s
  fit <- lm(y ~ x[, s] - 1)
  s_EBIC_s <- BIC(fit) + 2 * lchoose(p, length(s))
  EBIC_s_min <- s_EBIC_s
  for (iter in 1:M) {
    for(i in 2:p){
      plus_s <- rep(FALSE,p); minus_s <- rep(FALSE,p)
      plus_s[s] <- TRUE; minus_s[s] <- TRUE
      if(plus_s[i]==TRUE){
        EBIC_plus <- s_EBIC_s
        minus_s[i] <- FALSE
        minus_fit <- lm(y ~ x[,minus_s] - 1)
        EBIC_minus <- BIC(minus_fit) + 2 * lchoose(p, sum(minus_s))
      }else{
        EBIC_minus <- s_EBIC_s
        plus_s[i] <- TRUE
        plus_fit <- lm(y ~ x[,plus_s] - 1)
        EBIC_plus <- BIC(plus_fit) + 2 * lchoose(p, sum(plus_s))
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
    }
  }
  print("stochastic end")
  stochastic_end <- Sys.time()
  time_stochastic_arr[temp] <- difftime(Sys.time(),stochastic_start,units = "secs")
  stochastic_s <- s_min
  stochastic_EBIC <- EBIC_s_min
  EBIC_arr[temp,2] <- stochastic_EBIC  
  bool_vector <- rep(FALSE, p)
  bool_vector[stochastic_s] <- TRUE
  stochastic_index_arr[temp,] <- bool_vector

  #hybrid algorithm
  hybrid_start <- Sys.time()
  M <- 100
  init <- sample(2:p, 1)
  s_h <- c(1)
  fit <- lm(y ~ x[, s_h] - 1)
  EBIC_h <- BIC(fit) + 2 * g * lchoose(p, length(s_h))
  global_loop <- TRUE; loop <- TRUE
  while(global_loop){
    global_loop <- FALSE
    while(loop){
      EBIC_g <- EBIC_h
      loop <- FALSE
      for(i in s_h){
        st <- sort(union(1, setdiff(s_h, 1)))
        fit <- lm(y ~ x[, st] - 1)
        EBIC <- BIC(fit) + 2 * g * lchoose(p, length(st))
        if(EBIC_g > EBIC){
          EBIC_g <- EBIC
          s_g <- st
        }
      }
      for(i in setdiff(c(2:p), s_h)){
        st <- sort(c(s_h, i))
        fit <- lm(y ~ x[, st] - 1)
        EBIC <- BIC(fit) + 2 * g * lchoose(p, length(st))
        if(EBIC_g > EBIC){
          EBIC_g <- EBIC
          s_g <- st
        }
      }
      if(EBIC_h > EBIC_g){
        EBIC_h <- EBIC_g
        s_h <- s_g
        loop <- TRUE
        global_loop <- TRUE
      }
      if(length(s_h) > n^(2/3)){
        loop <- FALSE
      }
    }
    s <- s_h;
    EBIC_s <- EBIC_h
    condition <- FALSE
    for(iter in 1:M){
      for(i in 2:p){
        plus_s <- rep(FALSE, p); minus_s <- rep(FALSE, p)
        plus_s[s] <- TRUE; minus_s[s] <- TRUE
        if(plus_s[i] == TRUE){
          EBIC_plus <- EBIC_s
          minus_s[i] <-FALSE
          minus_fit <- lm(y ~ x[,minus_s] - 1)
          EBIC_minus <- BIC(minus_fit) + 2 * g * lchoose(p, sum(minus_s))
        }else {
           EBIC_minus <- EBIC_s
           plus_s[i] <- TRUE
           plus_fit <- lm(y ~ x[,plus_s] - 1)
           EBIC_plus <- BIC(plus_fit) + 2 * g * lchoose(p, sum(plus_s))
        }
        if(sum(plus_s) <= n^(2/3)){
          if(EBIC_h > EBIC_plus){
            s_h <- all_index[plus_s]
            EBIC_h <- EBIC_plus
            global_loop <- TRUE
            condition <- TRUE
            break
          }
        }
        if(EBIC_h > EBIC_minus){
          s_h <- all_index[minus_s]
          EBIC_h <- EBIC_minus
          global_loop <- TRUE
          condition <- TRUE
          break
        }
        prob <- 1 / (1 + exp((-0.5 * EBIC_minus) + (0.5 * EBIC_plus)))
        z <- rbinom(1, 1, prob)
        if((sum(plus_s) > n^(2/3)) & (z == 1)){

        }else {
           if(z == 1){
            s <- all_index[plus_s]
            EBIC_s <- EBIC_plus
           }else {
            s <- all_index[minus_s]
            EBIC_s <- EBIC_minus
           }
        }
      }
      if(iter == M){
        global_loop <- FALSE
      }
      if(condition){
        break
      }
    }
  }
  print("hybrid end")
  hybrid_end <- Sys.time()
  time_hybrid_arr[temp] <- difftime(Sys.time(),hybrid_start, units = "secs")
  hybrid_s <- s_h
  hybrid_EBIC <- EBIC_h
  EBIC_arr[temp, 3] <- hybrid_EBIC
  bool_vector <- rep(FALSE, p)
  bool_vector[hybrid_s] <- TRUE
  hybrid_index_arr[temp, ] <- bool_vector
  print(temp)
}

for(i in 1:times){
  if((EBIC_arr[i, 1] == EBIC_arr[i, 2]) & (EBIC_arr[i, 1] == EBIC_arr[i, 3]) & (EBIC_arr[i, 2] == EBIC_arr[i, 3])){
    all_same <- all_same + 1
  }
}
for(i in 1:times){
  min_value <- min(EBIC_arr[i,])
  if(EBIC_arr[i, 1] == min_value){
    win_greedy <- win_greedy + 1
  }if(EBIC_arr[i, 2] == min_value){
    win_stochastic <- win_stochastic + 1
  }if(EBIC_arr[i, 3] == min_value){
    win_hybrid <- win_hybrid + 1
  }
}

cat("greedy the best performed ratio: ", win_greedy / times)
cat("stochastic the best performed ratio: ", win_stochastic / times)
cat("hybrid the best performed ratio: ", win_hybrid / times)
cat("same performed ratio: ", all_same)
cat("greedy average running time: ", mean(time_greedy_arr))
cat("stochastic average running time: ", mean(time_stochastic_arr))
cat("hybrid average running time: ", mean(time_hybrid_arr))
for(i in 1:times){
  cat("greedy", i, "th select model: ", all_index[greedy_index_arr[i,]],"\n")
  cat("stochastic", i, "th select model: ", all_index[stochastic_index_arr[i,]],"\n")
  cat("hybrid", i, "th select model: ", all_index[hybrid_index_arr[i,]],"\n")
  cat(i,"th EBIC: ",EBIC_arr[i,],"\n")
  cat("\n")
}
#nolint end