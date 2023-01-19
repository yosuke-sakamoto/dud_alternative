#+ message = F
library(tidyverse)
library(data.table)


#'# data loading
dat <- fread("exp1_data_behavior.csv", sep = ",", header = T)
#dat <- subset(dat, dat$participant != "sub01")
dat <- na.omit(dat)


#'# functions
fs_dnm <- function(val3, corr) {
    guess <- c(0.5, 2.22 * 10e-3)
    params <- list("val3" = val3, "corr" = corr) # input vector
    
    fit <- suppressWarnings(optim(par = guess, fn = deviance, gr = NULL, method = "L-BFGS-B", parameters = params,
                                  lower = c(0.001, -5), upper = c(20, 5),
                                  control = list("maxit" = 100000,
                                                 "parscale" = c(0.004, 0.0002))))
    
    est <- data.frame(sH = fit$par[1], w = fit$par[2], dev = fit$value)
    return(est)
}

deviance <- function(x, parameters) {
    val3 <- parameters$val3
    corr <- parameters$corr
    
    sH = x[1]
    w = x[2]
    
    Vcu <- as.data.frame(cbind(100, 90, val3))/100
    colnames(Vcu) <- c("V1", "V2", "V3")
    Vcu <- mutate(Vcu, meanVcu = ifelse(V3 > 0, (V1 + V2 + V3)/3, (V1 + V2)/2))
    Vcu <- mutate(Vcu, normalizer = sH + w * meanVcu)
    M <- Vcu[, 1:3]/Vcu$normalizer
    M <- mutate(M, P = exp(V1) / (exp(V1) + exp(V2) + exp(V3)))
    M <- mutate(M, P = ifelse(P < .01, .01, P)) #adjustments to avoid punishing models too much for very unlikely predictions
    M <- mutate(M, P = ifelse(P > .99, .99, P)) #adjustments to avoid punishing models too much for very unlikely predictions
    M$corr <- corr
    
    dev = -2 * sum(M$corr * log(M$P) + (1 - M$corr) * log(1 - M$P))
    return(dev)
}


#'# individual fittings
fits <- c()

for (sub in unique(dat$subj)) {
    d <- subset(dat, dat$subj == sub)
    val3 <- d$val3
    corr <- d$corr
    f <- tryCatch(fs_dnm(val3, corr), error = function(e){f = cbind(NA, NA, NA)}) 
    colnames(f) <- c("sH", "w", "dev")
    fits <- try(rbind(fits, cbind(f, sub)))
    fits <- na.omit(fits)
    fits$sH <- as.numeric(fits$sH)
    fits$w <- as.numeric(fits$w)
    fits$dev <- as.numeric(fits$dev)
}


#'# visualization 
prediction <- c()
m1 <- mean(c(1, 0.9))
m2 <- mean(c(1, 0.9, 0.27))
m3 <- mean(c(1, 0.9, 0.45))
m4 <- mean(c(1, 0.9, 0.63))
m5 <- mean(c(1, 0.9, 0.76))
m6 <- mean(c(1, 0.9, 0.86))

for (subj in unique(fits$sub)) {
    f <- subset(fits, fits$sub == subj)
    
    pred <- c(
        exp(1 / (f[1] + f[2] * m1)) / ( exp(1 / (f[1] + f[2] * m1)) + exp(0.9 / (f[1] + f[2] * m1)) + exp(0 / (f[1] + f[2] * m1)) ),
        exp(1 / (f[1] + f[2] * m2)) / ( exp(1 / (f[1] + f[2] * m2)) + exp(0.9 / (f[1] + f[2] * m2)) + exp(0.27 / (f[1] + f[2] * m2)) ),
        exp(1 / (f[1] + f[2] * m3)) / ( exp(1 / (f[1] + f[2] * m3)) + exp(0.9 / (f[1] + f[2] * m3)) + exp(0.45 / (f[1] + f[2] * m3)) ),
        exp(1 / (f[1] + f[2] * m4)) / ( exp(1 / (f[1] + f[2] * m4)) + exp(0.9 / (f[1] + f[2] * m4)) + exp(0.63 / (f[1] + f[2] * m4)) ),
        exp(1 / (f[1] + f[2] * m5)) / ( exp(1 / (f[1] + f[2] * m5)) + exp(0.9 / (f[1] + f[2] * m5)) + exp(0.76 / (f[1] + f[2] * m5)) ),
        exp(1 / (f[1] + f[2] * m6)) / ( exp(1 / (f[1] + f[2] * m6)) + exp(0.9 / (f[1] + f[2] * m6)) + exp(0.86 / (f[1] + f[2] * m6)) ))
    
    pred <- as.numeric(pred)
    pred <- as.data.frame(pred)
    prediction <- rbind(prediction, pred)
}

dat %>%
    group_by(subj, val3) %>%
    summarise(Accuracy = mean(corr)) -> acc

acc <- subset(acc, acc$subj %in% unique(fits$sub))
acc <- cbind(acc, prediction)
acc$pred <- as.numeric(acc$pred)

mean_dnm <- ggplot(acc) + geom_point(aes(x = val3, y = Accuracy, color = subj)) +
    geom_line(mapping = aes(x = val3, y = pred, color = subj)) + ylim(0.49, 0.9) + xlab("Dud stimulus") + ggtitle("mean_dnm") + guides(color = guide_legend(title = NULL))
mean_dnm

fits_mean_dnm <- fits
fits_mean_dnm <- mutate(fits_mean_dnm, aic = dev + 4)