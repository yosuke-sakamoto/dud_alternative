#+ message = F
library(tidyverse)
library(data.table)


#'# data loading
dat <- fread("exp1_data_behavior.csv", sep = ",", header = T)
dat <- mutate(dat, Corr = ifelse(ChosenITM == CorrectITM, 1, 0))
#dat <- subset(dat, dat$participant != "sub01")
dat <- na.omit(dat)


#'# functions
fs_logit <- function(Val3, Corr) {
    guess <- 0.5
    params <- list("Val3" = Val3, "Corr" = Corr) # input vector
    
    fit <- suppressWarnings(optim(par = guess, fn = deviance, gr = NULL, method = "L-BFGS-B", parameters = params,
                                  upper = 20,
                                  control = list("maxit" = 100000)))
    
    est <- data.frame(sH = fit$par[1], w = 0, dev = fit$value)
    return(est)
}

deviance <- function(x, parameters) {
    Val3 <- parameters$Val3
    Corr <- parameters$Corr
    
    sH = x[1]
    w = 0
    
    Vcu <- as.data.frame(cbind(100, 90, Val3))/100
    colnames(Vcu) <- c("V1", "V2", "V3")
    Vcu <- mutate(Vcu, sumVcu = V1 + V2 + V3)
    Vcu <- mutate(Vcu, normalizer = sH + w * sumVcu)
    M <- Vcu[, 1:3]/Vcu$normalizer
    M <- mutate(M, P = exp(V1) / (exp(V1) + exp(V2) + exp(V3)))
    M <- mutate(M, P = ifelse(P < .01, .01, P)) #adjustments to avoid punishing models too much for very unlikely predictions
    M <- mutate(M, P = ifelse(P > .99, .99, P)) #adjustments to avoid punishing models too much for very unlikely predictions
    M$Corr <- Corr
    
    dev = -2 * sum(M$Corr * log(M$P) + (1 - M$Corr) * log(1 - M$P))
    return(dev)
}


#'# individual fittings
fits <- c()

for (sub in unique(dat$participant)) {
    d <- subset(dat, dat$participant == sub)
    Val3 <- d$Val3
    Corr <- d$Corr
    f <- tryCatch(fs_logit(Val3, Corr), error = function(e){f = cbind(NA, NA, NA)}) 
    colnames(f) <- c("sH", "w", "dev")
    fits <- try(rbind(fits, cbind(f, sub)))
    fits <- na.omit(fits)
    fits$sH <- as.numeric(fits$sH)
    fits$w <- as.numeric(fits$w)
    fits$dev <- as.numeric(fits$dev)
}


#'# visualization 
prediction <- c()

for (subj in unique(fits$sub)) {
    f <- subset(fits, fits$sub == subj)
    
    pred <- c(
        exp(1 / (f[1] + f[2] * 1.90)) / ( exp(1 / (f[1] + f[2] * 1.90)) + exp(0.9 / (f[1] + f[2] * 1.90)) + exp(0 / (f[1] + f[2] * 1.90)) ),
        exp(1 / (f[1] + f[2] * 2.17)) / ( exp(1 / (f[1] + f[2] * 2.17)) + exp(0.9 / (f[1] + f[2] * 2.17)) + exp(0.27 / (f[1] + f[2] * 2.17)) ),
        exp(1 / (f[1] + f[2] * 2.35)) / ( exp(1 / (f[1] + f[2] * 2.35)) + exp(0.9 / (f[1] + f[2] * 2.35)) + exp(0.45 / (f[1] + f[2] * 2.35)) ),
        exp(1 / (f[1] + f[2] * 2.53)) / ( exp(1 / (f[1] + f[2] * 2.53)) + exp(0.9 / (f[1] + f[2] * 2.53)) + exp(0.63 / (f[1] + f[2] * 2.53)) ),
        exp(1 / (f[1] + f[2] * 2.66)) / ( exp(1 / (f[1] + f[2] * 2.66)) + exp(0.9 / (f[1] + f[2] * 2.66)) + exp(0.76 / (f[1] + f[2] * 2.66)) ),
        exp(1 / (f[1] + f[2] * 2.76)) / ( exp(1 / (f[1] + f[2] * 2.76)) + exp(0.9 / (f[1] + f[2] * 2.76)) + exp(0.86 / (f[1] + f[2] * 2.76)) ))
    
    pred <- as.numeric(pred)
    pred <- as.data.frame(cbind(pred, subj))
    prediction <- rbind(prediction, pred)
}

dat %>%
    group_by(participant, Val3) %>%
    summarise(Accuracy = mean(Corr)) -> acc

acc <- subset(acc, acc$participant %in% unique(fits$sub))
acc <- cbind(acc, prediction)
acc$pred <- as.numeric(acc$pred)

logit <- ggplot(acc) + geom_point(aes(x = Val3, y = Accuracy, color = participant)) +
    geom_line(mapping = aes(x = Val3, y = pred, color = participant)) + ylim(0.49, 0.9) + ggtitle("logit")
logit

fits_logit <- fits
fits_logit <- mutate(fits_logit, aic = dev + 2)