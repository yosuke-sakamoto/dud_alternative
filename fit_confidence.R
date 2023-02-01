### fitting confidence data
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

df <- read.csv("exp1_data_behavior.csv", header = TRUE)
df <- na.omit(df)
df$chosenItem <- ifelse(df$chosenItem == "target", 1, 
                        ifelse(df$chosenItem == "distractor", 2, 3))
df01 <- filter(df, subj == "sub01" & chosenItem != 3)

fit_conf <- function(simdat, empdat, r) {
    
    guess <- c(0, 2.5/(0.15 + r), 4/(0.15 + r))
    
    params <- list("r" = r,
                   "m1_correct" =   filter(simdat, choice == 1)$m1_post,
                   "m2_correct" =   filter(simdat, choice == 1)$m2_post,
                   "m1_incorrect" = filter(simdat, choice == 2)$m1_post,
                   "m2_incorrect" = filter(simdat, choice == 2)$m2_post,
                   "conf_correct" = filter(empdat, chosenItem == 1)$conf,
                   "conf_incorrect" = filter(empdat, chosenItem == 2)$conf)
    
    fit <- optim(par = guess, fn = fit_conf_ks, gr = NULL, method = "L-BFGS-B", parameters = params,
                 lower = c(-5, -5, 0),
                 upper = c(20, 20, 20),
                 control = list("maxit" = 100000, "parscale" = c(10, 10, 10)))
    
    est <- data.frame(c1 = fit$par[1], c2 = fit$par[2], c3 = fit$par[3], ks = fit$value)
    return(est)
}

fit_conf_ks <- function(x, parameters) {
    c1 <- x[1] 
    c2 <- x[2]
    c3 <- x[3]
    r  <- parameters$r 
    
    cv_incorrect <- parameters$m1_incorrect - r * parameters$m2_incorrect
    cv_correct <- parameters$m1_correct - r * parameters$m2_correct
    simconf_incorrect <- c(rep(1, sum(cv_incorrect < c1)), rep(2, sum(cv_incorrect >= c1 & cv_incorrect < c2)), 
                           rep(3, sum(cv_incorrect >= c2 & cv_incorrect < c3)), rep(4, sum(cv_incorrect >= c3)))
    simconf_correct <-  c(rep(1, sum(cv_correct < c1)), rep(2, sum(cv_correct >= c1 & cv_correct < c2)), 
                          rep(3, sum(cv_correct >= c2 & cv_correct < c3)), rep(4, sum(cv_correct >= c3)))
    simconf <- c(simconf_correct, -1 * simconf_incorrect)
    empconf <- c(parameters$conf_correct, -1 * parameters$conf_incorrect)
    
    ks <- ks.test(simconf, empconf)$statistic # Kolmogorov-Smirnov test
    
    if (is.nan(ks)) {
        ks <- Inf
    }
    
    # define constraints
    if (x[1] < x[2] & x[2] < x[3]) {
        return(ks)
    } else {
        return(1)
    }
}

fit_conf_ll <- function(x, parameters) {
    c1 <- x[1] 
    c2 <- x[2]
    c3 <- x[3]
    r  <- parameters$r 
    
    cv_incorrect <- parameters$m1_incorrect - r * parameters$m2_incorrect
    cv_correct <- parameters$m1_correct - r * parameters$m2_correct
    
    simconf_incorrect <- c(sum(cv_incorrect < c1), sum(cv_incorrect >= c1 & cv_incorrect < c2), 
                           sum(cv_incorrect >= c2 & cv_incorrect < c3), sum(cv_incorrect >= c3))
    simconf_correct <- c(sum(cv_correct < c1), sum(cv_correct >= c1 & cv_correct < c2), 
                         sum(cv_correct >= c2 & cv_correct < c3), sum(cv_correct >= c3))
    sim_incorrect_rate <- simconf_incorrect + 0.5 / sum(simconf_incorrect)
    sim_correct_rate <- simconf_correct + 0.5 / sum(simconf_correct)
    
    
    empconf_incorrect <- parameters$conf_incorrect
    empconf_correct <- parameters$conf_correct
    emp_incorrect <- c(sum(empconf_incorrect == 1), sum(empconf_incorrect == 2), sum(empconf_incorrect == 3), sum(empconf_incorrect == 4))
    emp_correct <- c(sum(empconf_correct == 1), sum(empconf_correct == 2), sum(empconf_correct == 3), sum(empconf_correct == 4))
    
    ll <- sum(emp_incorrect * log(sim_incorrect_rate) + emp_correct * log(sim_correct_rate))
    
    if (is.nan(ll)) {
        ll <- -9999
    }
    
    ll <- -ll
    
    # define constraints
    if (x[1] < x[2] & x[2] < x[3]) {
        return(ll)
    } else {
        return(9999)
    }
}

# grid search for weight (weightによってcvのスケールが大きく変わるため)
weight <- seq(0.9, 1, 0.01)
fit <- foreach(r = weight, .combine = "rbind", .packages = c("magrittr", "dplyr")) %dopar% {
    f <- fit_conf(simdat = simdat, empdat = df01, r = r)
    f$r <- r
    print(f)
}
fit

# model prediction
f <- fit[which.min(fit$ks), ] # winning model
f <- fit[11, ]

c1 <- as.numeric(f[1])
c2 <- as.numeric(f[2])
c3 <- as.numeric(f[3])
r  <- as.numeric(f[5])

simdat %>%
    mutate(cv = m1_post - r * m2_post, 
           simconf = ifelse(cv < c1, 1,
                            ifelse(cv >= c1 & cv < c2, 2,
                                   ifelse(cv >= c2 & cv < c3, 3, 4)))) %>%
    group_by(choice, condition) %>%
    summarise(mean_simconf = mean(simconf)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_simconf, color = factor(choice))) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86))