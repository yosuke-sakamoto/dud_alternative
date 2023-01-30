# import library
library(doParallel)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tictoc)

cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

dr    <- 0.0001
sigma <- 0.06
n_dt  <- 250
w     <- -0.0022

Trial <- 10000
th <- 1.2
Sample <- 3000
condition <- c(0, 27, 45, 63, 76, 86)
simdat <- c()

tic()
for (j in condition) {
    
    d <- foreach(i = 1:Trial, .combine = "rbind", .packages = c("magrittr", "dplyr")) %dopar% {
        Stimulus <- c(100, 90, j)
        Stimulus <- Stimulus / (1 + w * sum(Stimulus))
        Stimulus <- Stimulus * dr
        nAlt <- sum(Stimulus != 0)
        
        accum <- data.frame(a1 = c(0, rep(Stimulus[1], Sample)), 
                            a2 = c(0, rep(Stimulus[2], Sample)), 
                            a3 = c(0, rep(Stimulus[3], Sample)))
        accum[-1, 1] <- accum[-1, 1] + rnorm(Sample, 0, sigma)
        accum[-1, 2] <- accum[-1, 2] + rnorm(Sample, 0, sigma)
        accum[-1, 3] <- accum[-1, 3] + rnorm(Sample, 0, sigma)
        accum <- cumsum(accum)
        accum <- accum[, 1:nAlt]
        accum <- mutate(accum, m1 = apply(accum, 1, max), 
                        m2 = apply(accum, 1, function(x){return(max(x[-which.max(x)]))}), 
                        delta = m1 - m2) # 時刻を変えた時に刺激が入れ替わる
        simrt <- which(accum$delta > th)[1]
        choice <- max.col(accum[, 1:nAlt])[simrt]
        
        pre <- ifelse(simrt < 300, simrt, 300)
        post <- 300
        
        c(simrt, choice, accum$a1[simrt], accum$a2[simrt], accum$m1[simrt], accum$m2[simrt], 
          accum$delta[simrt], accum$delta[simrt + pre], accum$delta[simrt + post])
    }
    
    d <- as.data.frame(d)
    colnames(d) <- c("rt", "choice", "a1", "a2", "m1", "m2", "delta", "pre_delta", "post_delta")
    d$rt   <- d$rt + n_dt
    d$rt[d$rt > 3000] <- NA
    d <- na.omit(d)
    d$rt <- d$rt/1000 #####
    d$condition <- j
    simdat <- rbind(simdat, d)
    
}
toc()

# model visualization
simdat %>%
    complete(choice, condition) %>%
    ggplot() + geom_histogram(aes(x = rt, fill = factor(choice))) + facet_wrap(. ~ factor(condition) + factor(choice), nrow = 6)

simdat %>%
    group_by(choice, condition) %>%
    mutate(n1 = n()) %>%
    ungroup(choice, condition) %>%
    group_by(condition) %>%
    mutate(n2 = n()) %>%
    select(n1, n2, choice, condition) %>%
    distinct() %>%
    mutate(p = n1/n2) %>%
    ggplot() + geom_line(aes(x = condition, y = p, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86))

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_rt = mean(rt)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_rt, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86))

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_m1 = mean(m1)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_m1, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86))

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_delta = mean(delta)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_delta, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86))

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_pre_delta = mean(pre_delta)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_pre_delta, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86))

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_post_delta = mean(post_delta)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_post_delta, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86))