# import
library(doParallel)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tictoc)
library(cowplot)

cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

dr    <- 0.0001
sigma <- 0.06
n_dt  <- 250
w     <- -0.0022

Trial <- 10000
th <- 1.2
threshold <- rnorm(Trial, 1.2, 0.1) # jitter for threshold
Sample <- 3000
condition <- c(0, 27, 45, 63, 76, 86)
simdat <- c()

tic()
for (j in condition) {
    
    d <- foreach(i = 1:Trial, .combine = "rbind", .packages = c("magrittr", "dplyr")) %dopar% {
        # th <- threshold[i] # jitter for threshold
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
        accum <- accum[, 1:nAlt]
        accum <- cumsum(accum)[-1, ]
        accum <- mutate(accum, delta = apply(accum, 1, max) - apply(accum, 1, function(x){return(max(x[-which.max(x)]))}))
        
        simrt  <- which(accum$delta > th)[1]
        choice <- max.col(accum[, 1:nAlt])[simrt]
        third  <- max.col(-accum[, 1:nAlt])[simrt]
        second <- ifelse(nAlt == 3, setdiff(c(1, 2, 3), c(choice, third)), third)
        
        pre <- ifelse(simrt <= 300, simrt - 1, 300)
        post <- 300
        
        as.numeric(c(simrt, choice,
                     accum[simrt, choice],        accum[simrt, second],        accum[simrt, third],
                     accum[simrt - pre, choice],  accum[simrt - pre, second],  accum[simrt - pre, third],                                      
                     accum[simrt + post, choice], accum[simrt + post, second], accum[simrt + post, third]))                                                                   
        
    }
    
    d <- as.data.frame(d)
    colnames(d) <- c("rt", "choice", "m1", "m2", "m3", "m1_pre", "m2_pre", "m3_pre", "m1_post", "m2_post", "m3_post")
    d$rt   <- d$rt + n_dt
    #d$rt[d$rt > 3000] <- NA
    #d <- na.omit(d)
    d$condition <- j
    simdat <- rbind(simdat, d)
    
}
toc()

# model visualization
simdat %>% mutate(delta = m1 - m2, delta_pre = m1_pre - m2_pre, delta_post = m1_post - m2_post) -> simdat
simdat %>%
    complete(choice, condition) %>%
    ggplot() + geom_histogram(aes(x = rt, fill = factor(choice))) + facet_wrap(. ~ factor(condition) + factor(choice), nrow = 6) -> p1
p1

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
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) -> p2
p2

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_rt = mean(rt)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_rt, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) -> p3
p3

### at choice
simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_m1 = mean(m1)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_m1, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylim(0, 8) -> p4
p4

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_m2 = mean(m2)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_m2, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylim(0, 8) -> p5
p5

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_m3 = mean(m3)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_m3, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylim(0, 8) -> p6
p6

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_delta = mean(delta)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_delta, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) -> p7
p7

### before choice
simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_m1_pre = mean(m1_pre)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_m1_pre, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylim(0, 8) -> p8
p8

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_m2_pre = mean(m2_pre)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_m2_pre, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylim(0, 8) -> p9
p9

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_m3_pre = mean(m3_pre)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_m3_pre, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylim(0, 8) -> p10
p10

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_delta_pre = mean(delta_pre)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_delta_pre, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) -> p11
p11

### post choice
simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_m1_post = mean(m1_post)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_m1_post, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) -> p12
p12

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_m2_post = mean(m2_post)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_m2_post, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) -> p13
p13

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_m3_post = mean(m3_post)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_m3_post, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) -> p14
p14

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_delta_post = mean(delta_post)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_delta_post, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) -> p15
p15