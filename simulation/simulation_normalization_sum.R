library(tidyverse)
library(doParallel)

cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

Trials <- 3000
Sample <- 3000
sigma_dr <- 0
th <- 1
dud <- seq(0, 90, by = 15)
results_all <- list()
sigma_signal <- 0.05

results_all <- foreach (d = 1:length(dud), .combine = "rbind") %do% {# loop over conditions (values of dud-alternative: 0 (equivalent to 2AFC) to 90 (equal to the second largest value))
  # initialize
  Stimulus <- c(100, 90, dud[d]) / 20000
  theta <- Stimulus / (sigma_dr + sum(Stimulus))
  Nalt <- sum(Stimulus != 0)
  # 15 sec.
  foreach (i = 1:Trials, .combine = "rbind") %do% { # loop over trials
    #set accumulator
    accum <- data.frame(a1 = c(0, rep(c(theta * Stimulus)[1], Sample)), 
                        a2 = c(0, rep(c(theta * Stimulus)[2], Sample)), 
                        a3 = c(0, rep(c(theta * Stimulus)[3], Sample))) 
    accum[-1, 1] <- accum[-1, 1] + rnorm(Sample, 0, sigma_signal)
    accum[-1, 2] <- accum[-1, 2] + rnorm(Sample, 0, sigma_signal)
    accum[-1, 3] <- accum[-1, 3] + rnorm(Sample, 0, sigma_signal)
    accum <- cumsum(accum)
    if (dud[d] == 0){# when 2AFC
      accum <- accum[, 1:2]
    }
    
    # calculate the differences between accumulators
    delta <- foreach (j = 1:Nalt, .combine = "cbind") %dopar% {
      accum[, j] - apply(as.data.frame(accum[, -j]), 1, max)
    }
    rownames(delta) <- 1:nrow(delta)
    
    rt_cand <- vector()
    for (j in 1:Nalt){
      rt_cand[j] <- which(delta[, j] > th)[1]
    }
    
    chooseDud <- 0
    corr <- 0
    if (sum(!is.na(rt_cand)) == 0){
      rt <- NA
      choice <- NA
      conf <- NA
    } else {
      rt <- min(rt_cand, na.rm = TRUE)
      choice <- which(rt_cand == rt)
      if (length(choice) != 1){
        choice <- sample(choice, 1)
      }
      conf <- accum[rt, choice]
      if (choice == 1){
        corr <- 1
      }
      if (choice == 3){
        chooseDud <- 1
      }
    }
    
    c(dud[d], rt, choice, conf, corr, chooseDud)
  }
}

rownames(results_all) <- 1:nrow(results_all)
colnames(results_all) <- c("condition", "rt", "choice", "conf", "corr", "chooseDud")
results_all <- as.data.frame(results_all)
######
accuracy <- vector()
conf_all <- vector()
rt_mean <- vector()
for (i in 1:length(dud)) {
  result <- na.omit(filter(results_all, condition == dud[i]))
  #result$conf <- if_else(result$choice == 1, result$conf, -result$conf)
  accuracy[i] <- sum(result$corr) / nrow(result)
  conf_all[i] <- mean(result$conf)
  rt_mean[i] <- mean(result$rt)
}

df <- data.frame(condition = dud, acc = accuracy, conf = conf_all, rt = rt_mean)
ggplot(df, aes(x = condition, y = acc)) + geom_line()
ggplot(df, aes(x = condition, y = conf)) + geom_line()
ggplot(df, aes(x = condition, y = rt)) + geom_line()
