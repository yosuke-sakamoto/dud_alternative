# import library
library(tidyverse)
library(doParallel)
library(dplyr)
library(tictoc)

cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

# import data
df <- read.csv("exp1_data_behavior.csv", header = TRUE)
df <- na.omit(df)
df$chosenItem <- ifelse(df$chosenItem == "target", 1, 
                        ifelse(df$chosenItem == "distractor", 2, 3))

chisq <- function(par, dat1){
  
  dr    <- par[1]
  sigma <- par[2]
  n_dt  <- par[3]
  w     <- par[4]
  
  dat1$rt <- dat1$rt * 100
  rt_q <- quantile(dat1$rt, probs = c(0.1, 0.3, 0.5, 0.7, 0.9))
  th <- 1
  Sample <- 300
  chisq_stats <- c()
  
  for (c in unique(dat1$condition)){
    
    dat2 <- filter(dat1, condition == c)
    dudVal <- unique(dat2$val3)
    
    d <- foreach (i = 1:nrow(dat2), .combine = "rbind", .packages = "tidyverse") %dopar% {
      Stimulus <- c(100, 90, dudVal)
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
    
    d <- as.data.frame(d)
    colnames(d) <- c("rt", "choice")
    d$rt   <- d$rt + n_dt
    d$rt[d$rt > 300] <- NA
    d <- na.omit(d)
    
    d_tar  <- filter(d, choice == 1)
    d_dist <- filter(d, choice == 2)
    d_dud  <- filter(d, choice == 3)
    
    simrt_tar  <- c(sum(d_tar$rt <= rt_q[1]), 
                    sum(d_tar$rt > rt_q[1] & d_tar$rt <= rt_q[2]), 
                    sum(d_tar$rt > rt_q[2] & d_tar$rt <= rt_q[3]), 
                    sum(d_tar$rt > rt_q[3] & d_tar$rt <= rt_q[4]),
                    sum(d_tar$rt > rt_q[4] & d_tar$rt <= rt_q[5]),
                    sum(d_tar$rt > rt_q[5]))
    simrt_dist <- c(sum(d_dist$rt <= rt_q[1]), 
                    sum(d_dist$rt > rt_q[1] & d_dist$rt <= rt_q[2]), 
                    sum(d_dist$rt > rt_q[2] & d_dist$rt <= rt_q[3]), 
                    sum(d_dist$rt > rt_q[3] & d_dist$rt <= rt_q[4]),
                    sum(d_dist$rt > rt_q[4] & d_dist$rt <= rt_q[5]),
                    sum(d_dist$rt > rt_q[5]))
    simrt_dud  <- c(sum(d_dud$rt <= rt_q[1]), 
                    sum(d_dud$rt > rt_q[1] & d_dud$rt <= rt_q[2]), 
                    sum(d_dud$rt > rt_q[2] & d_dud$rt <= rt_q[3]), 
                    sum(d_dud$rt > rt_q[3] & d_dud$rt <= rt_q[4]),
                    sum(d_dud$rt > rt_q[4] & d_dud$rt <= rt_q[5]),
                    sum(d_dud$rt > rt_q[5]))
    tab_simrt <- matrix(c(simrt_tar, simrt_dist, simrt_dud), nrow = 3, byrow = TRUE)
    
    dat_tar  <- filter(dat2, chosenItem == 1)
    dat_dist <- filter(dat2, chosenItem == 2)
    dat_dud  <- filter(dat2, chosenItem == 3)
    
    rt_tar <- c(sum(dat_tar$rt <= rt_q[1]), 
                sum(dat_tar$rt > rt_q[1] & dat_tar$rt <= rt_q[2]), 
                sum(dat_tar$rt > rt_q[2] & dat_tar$rt <= rt_q[3]), 
                sum(dat_tar$rt > rt_q[3] & dat_tar$rt <= rt_q[4]),
                sum(dat_tar$rt > rt_q[4] & dat_tar$rt <= rt_q[5]),
                sum(dat_tar$rt > rt_q[5]))
    rt_dist <- c(sum(dat_dist$rt <= rt_q[1]), 
                 sum(dat_dist$rt > rt_q[1] & dat_dist$rt <= rt_q[2]), 
                 sum(dat_dist$rt > rt_q[2] & dat_dist$rt <= rt_q[3]), 
                 sum(dat_dist$rt > rt_q[3] & dat_dist$rt <= rt_q[4]),
                 sum(dat_dist$rt > rt_q[4] & dat_dist$rt <= rt_q[5]),
                 sum(dat_dist$rt > rt_q[5]))
    rt_dud  <- c(sum(dat_dud$rt <= rt_q[1]), 
                 sum(dat_dud$rt > rt_q[1] & dat_dud$rt <= rt_q[2]), 
                 sum(dat_dud$rt > rt_q[2] & dat_dud$rt <= rt_q[3]), 
                 sum(dat_dud$rt > rt_q[3] & dat_dud$rt <= rt_q[4]),
                 sum(dat_dud$rt > rt_q[4] & dat_dud$rt <= rt_q[5]),
                 sum(dat_dud$rt > rt_q[5]))
    tab_rt <- matrix(c(rt_tar, rt_dist, rt_dud), nrow = 3, byrow = TRUE)
    
    chisq_cond <- (tab_simrt - tab_rt) ^ 2 / tab_rt
    chisq_cond <- chisq_cond[is.finite(chisq_cond)]
    
    chisq_stats <- c(chisq_stats, sum(chisq_cond))
  }
  
  return (sum(chisq_stats))
}

#maximize objective function
fit_dnm <- function(dat_behav){
  init_par <- c(0.001, 0.05, 20, 0.01)
  
  fit <- suppressWarnings(optim(par = init_par, fn = chisq, gr = NULL, 
                                method = "L-BFGS-B", dat1 = dat_behav,
                                lower = c(10e-6, 10e-4, 5, -0.003), upper = c(10e-1, 1, 100, 5),
                                control = list("maxit" = 100000, "parscale" = c(10, 10, 10, 1))))
  
  est <- data.frame(dr = fit$par[1], sigma = fit$par[2], n_dt = fit$par[3], w = fit$par[4], chisq = fit$value)
  return(est)
}

#individual fittings
fit_indiv <- data.frame()
for (i in unique(df$subj)){
  df_indiv <- filter(df, subj == i)
  fit <- fit_dnm(df_indiv)
  fit_indiv <- bind_rows(fit_indiv, fit)
}