# import library
library(tidyverse)
library(doParallel)
library(dplyr)
library(tictoc)

cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

# import data
df_behav <- read.csv("exp1_data_behavior.csv", header = TRUE)
df_behav <- na.omit(df_behav)


chisq <- function(par, dat){
  
  dr    <- par[1]
  sigma <- par[2]
  n_dt  <- par[3]
  w     <- par[4]
  
  rt_dat <- dat$rt * 1000
  rt_q <- quantile(dat$rt, probs = c(0.1, 0.3, 0.5, 0.7, 0.9)) * 1000
  Trial <- nrow(dat)
  th <- 1
  Sample <- 3000
  
  d <- foreach (i = 1:Trial, .combine = "rbind", .packages = "tidyverse") %dopar% {
    Stimulus <- c(100, 90, dat$val3[i])
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
                    delta = m1 - m2)
    simrt <- which(accum$delta > th)[1]
    choice <- max.col(accum[, 1:nAlt])[simrt]
    
    c(simrt, choice)
  }
  
  d <- na.omit(as.data.frame(d))
  colnames(d) <- c("rt", "choice")
  d$rt <- d$rt + n_dt
  d_corr <- filter(d, choice == 1)
  d_incorr <- filter(d, choice != 1)
  
  simrt_corr <- c(sum(d_corr$rt < rt_q[1]), 
                  sum(d_corr$rt > rt_q[2] & d_corr$rt < rt_q[3]), 
                  sum(d_corr$rt > rt_q[3] & d_corr$rt < rt_q[4]),
                  sum(d_corr$rt > rt_q[4] & d_corr$rt < rt_q[5]),
                  sum(d_corr$rt > rt_q[5]))
  simrt_incorr <- c(sum(d_incorr$rt < rt_q[1]), 
                    sum(d_incorr$rt > rt_q[2] & d_incorr$rt < rt_q[3]), 
                    sum(d_incorr$rt > rt_q[3] & d_incorr$rt < rt_q[4]),
                    sum(d_incorr$rt > rt_q[4] & d_incorr$rt < rt_q[5]),
                    sum(d_incorr$rt > rt_q[5]))
  tab_simrt <- matrix(c(simrt_corr, simrt_incorr), nrow = 2, byrow = TRUE)
  
  rt_corr <- c(sum(rt_dat < rt_q[1]), 
               sum(rt_dat > rt_q[2] & rt_dat < rt_q[3]), 
               sum(rt_dat > rt_q[3] & rt_dat < rt_q[4]),
               sum(rt_dat > rt_q[4] & rt_dat < rt_q[5]),
               sum(rt_dat > rt_q[5]))
  rt_incorr <- c(sum(rt_dat < rt_q[1]), 
                 sum(rt_dat > rt_q[2] & rt_dat < rt_q[3]), 
                 sum(rt_dat > rt_q[3] & rt_dat < rt_q[4]),
                 sum(rt_dat > rt_q[4] & rt_dat < rt_q[5]),
                 sum(rt_dat > rt_q[5]))
  tab_rt <- matrix(c(rt_corr, rt_incorr), nrow = 2, byrow = TRUE)
  
  chisq <- sum((tab_simrt - tab_rt) ^ 2 / tab_rt)
  
  return (chisq)
}

#maximize objective function
fit_dnm <- function(dat_behav){
  init_par <- c(0.00001, 0.05, 300, 0.01)
  
  fit <- suppressWarnings(optim(par = init_par, fn = chisq, gr = NULL, 
                                method = "L-BFGS-B", dat = dat_behav,
                                lower = c(10e-9, 10e-5, 10, -5), upper = c(10e-2, 1, 1000, 5),
                                control = list("maxit" = 100000, "parscale" = c(10, 1, 10, 1))))
  
  est <- data.frame(dr = fit$par[1], sigma = fit$par[2], n_dt = fit$par[3], w = fit$par[4], chisq = fit$value)
  return(est)
}

#individual fittings
fit_indiv <- data.frame()
for (i in unique(df_behav$subj)){
  df_behav_indiv <- filter(df_behav, subj == i)
  fit <- fit_dnm(df_behav_indiv)
  fit_indiv <- bind_rows(fit_indiv, fit)
}