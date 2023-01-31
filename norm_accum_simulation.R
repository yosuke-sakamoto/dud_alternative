# import ２AFCは別ルートで統合
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

Trial <- 3000
th <- 1.2
Sample <- 3000
condition <- c(27, 45, 63, 76, 86)
simdat <- c()

tic()
threshold <- rnorm(Trial, 1.2, 0.1) # jitter for threshold
for (j in condition) {
    
    d <- foreach(i = 1:Trial, .combine = "rbind", .packages = c("magrittr", "dplyr")) %dopar% {
        #th <- threshold[i] # jitter for threshold
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
        accum <- cumsum(accum)[-1, ]
       # accum <- accum[, 1:nAlt]
        accum %>% mutate(accum, first =  max.col(accum), third = max.col(-accum)) %>%
            rowwise() %>% mutate(second = setdiff(c(1, 2, 3), c(first, third))) %>%
            rowwise() %>% mutate(delta = pmax(a1, a2, a3) - sort(c(a1, a2, a3), decreasing = TRUE)[2]) -> accum
        
        simrt <- which(accum$delta > th)[1]
        choice <- max.col(accum[, 1:nAlt])[simrt]
        
        pre <- ifelse(simrt < 300, simrt - 1, 300)
        post <- 300
        
        as.numeric(c(simrt, choice,
          accum[simrt, accum$first[simrt]], accum[simrt, accum$second[simrt]], accum[simrt, accum$third[simrt]],
          accum[simrt - pre, accum$first[simrt]], accum[simrt - pre, accum$second[simrt]], accum[simrt - pre, accum$third[simrt]],                                      
          accum[simrt + post, accum$first[simrt]], accum[simrt + post, accum$second[simrt]], accum[simrt + post, accum$third[simrt]]))                                                                   
          
    }
    
    d <- as.data.frame(d)
    colnames(d) <- c("rt", "choice", "a1", "a2", "a3", "a1_pre", "a2_pre", "a3_pre", "a1_post", "a2_post", "a3_post")
    d$rt   <- d$rt + n_dt
    d$rt[d$rt > 3000] <- NA
    d <- na.omit(d)
    d$rt <- d$rt/1000
    d$condition <- j
    simdat <- rbind(simdat, d)
    
}
toc()

# model visualization
simdat %>% mutate(delta = a1 - a2, delta_pre = a1_pre - a2_pre, delta_post = a1_post - a2_post) -> simdat
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

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_a1 = mean(a1)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_a1, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylim(0, 8) -> p4
p4

simdat %>%
    group_by(choice, condition) %>%
    ggplot(aes(x = condition, y = a1, color = factor(choice))) + geom_point() +
    stat_summary(fun.y = "mean", geom = "line") +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylim(0, 8) -> p4
p4

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_a2 = mean(a2)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_a2, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylim(0, 8) -> p5
p5

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_a3 = mean(a3)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_a3, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylim(0, 8) -> p6
p6

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_delta = mean(delta)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_delta, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) -> p6
p6

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_delta_pre = mean(delta_pre)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_delta_pre, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) -> p7
p7

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_delta_post = mean(delta_post)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_delta_post, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) -> p8
p8 #平均だけではわからないのでgrid searchでモデリングする

plot(simdat$a1, simdat$a2_pre) # 3つの関係を図示
plot(simdat$a1, simdat$a2_post)