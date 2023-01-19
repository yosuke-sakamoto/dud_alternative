#+ message = F
library(tidyverse)
library(data.table)
library(evd)


#'# data loading
dat <- fread("exp1_data_behavior.csv", sep = ",", header = T)
#dat <- subset(dat, dat$participant != "sub01")
dat <- na.omit(dat)


#'# functions
fs_logit <- function(val3, corr) {
    guess <- 0.5
    params <- list("val3" = val3, "corr" = corr) # input vector
    
    fit <- suppressWarnings(optim(par = guess, fn = deviance, gr = NULL, method = "L-BFGS-B", parameters = params,
                                  upper = 20,
                                  control = list("maxit" = 100000)))
    
    est <- data.frame(sH = fit$par[1], w = 0, dev = fit$value)
    return(est)
}

deviance <- function(x, parameters) {
    val3 <- parameters$val3
    corr <- parameters$corr
    
    sH = x[1]
    w = 0
    
    Vcu <- as.data.frame(cbind(100, 90, val3))/100
    colnames(Vcu) <- c("V1", "V2", "V3")
    Vcu <- mutate(Vcu, sumVcu = V1 + V2 + V3)
    Vcu <- mutate(Vcu, normalizer = sH + w * sumVcu)
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
    f <- tryCatch(fs_logit(val3, corr), error = function(e){f = cbind(NA, NA, NA)}) 
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
    
    n1 <- exp(1 / (f[1] + f[2] * 1.90)) + exp(0.9 / (f[1] + f[2] * 1.90)) + exp(0 / (f[1] + f[2] * 1.90))
    n2 <- exp(1 / (f[1] + f[2] * 2.17)) + exp(0.9 / (f[1] + f[2] * 2.17)) + exp(0.27 / (f[1] + f[2] * 2.17))
    n3 <- exp(1 / (f[1] + f[2] * 2.35)) + exp(0.9 / (f[1] + f[2] * 2.35)) + exp(0.45 / (f[1] + f[2] * 2.35))
    n4 <- exp(1 / (f[1] + f[2] * 2.53)) + exp(0.9 / (f[1] + f[2] * 2.53)) + exp(0.63 / (f[1] + f[2] * 2.53))
    n5 <- exp(1 / (f[1] + f[2] * 2.66)) + exp(0.9 / (f[1] + f[2] * 2.66)) + exp(0.76 / (f[1] + f[2] * 2.66))
    n6 <- exp(1 / (f[1] + f[2] * 2.76)) + exp(0.9 / (f[1] + f[2] * 2.76)) + exp(0.86 / (f[1] + f[2] * 2.76))           
    
    # softmax transformation                                                                      
    pred <- c(
        exp(1 / (f[1] + f[2] * 1.90)) / n1,
        exp(1 / (f[1] + f[2] * 2.17)) / n2,
        exp(1 / (f[1] + f[2] * 2.35)) / n3,
        exp(1 / (f[1] + f[2] * 2.53)) / n4,
        exp(1 / (f[1] + f[2] * 2.66)) / n5,
        exp(1 / (f[1] + f[2] * 2.76)) / n6)
    
    # mean evidence for gumbel sdt
    evs <- rbind(
        cbind( 1 / (f[1] + f[2] * 1.90), 0.9 / (f[1] + f[2] * 1.90), 0 / (f[1] + f[2] * 1.90) ),
        cbind( 1 / (f[1] + f[2] * 2.17), 0.9 / (f[1] + f[2] * 2.17), 0.27 / (f[1] + f[2] * 2.17) ),
        cbind( 1 / (f[1] + f[2] * 2.35), 0.9 / (f[1] + f[2] * 2.35), 0.45 / (f[1] + f[2] * 2.35) ),
        cbind( 1 / (f[1] + f[2] * 2.53), 0.9 / (f[1] + f[2] * 2.53), 0.63 / (f[1] + f[2] * 2.53) ),
        cbind( 1 / (f[1] + f[2] * 2.66), 0.9 / (f[1] + f[2] * 2.66), 0.76 / (f[1] + f[2] * 2.66) ),
        cbind( 1 / (f[1] + f[2] * 2.76), 0.9 / (f[1] + f[2] * 2.76), 0.86 / (f[1] + f[2] * 2.76) ))
    
    pred <- as.numeric(pred)
    pred <- as.data.frame(pred)
    pred$id <- subj
    pred$condition <- c(0, 27, 45, 63, 76, 86)
    prediction <- rbind(prediction, cbind(pred, evs))
    
}


#'# accuracy
dat %>%
    group_by(subj, val3) %>%
    summarise(Accuracy = mean(corr)) -> acc

acc <- subset(acc, acc$subj %in% unique(fits$sub))
acc$pred <- prediction$pred

logit <- ggplot(acc) + geom_point(aes(x = val3, y = Accuracy, color = subj)) +
    geom_line(mapping = aes(x = val3, y = pred, color = subj)) + ylim(0.49, 0.9) + xlab("Dud stimulus") + ggtitle("logit")  + guides(color = guide_legend(title = NULL))
logit

fits_logit <- fits
fits_logit <- mutate(fits_logit, aic = dev + 4)


#'# confidence
colnames(prediction) <- c("pred", "id", "condition", "ev1", "ev2", "ev3")
prediction %>% arrange(condition) -> prediction
conf_dat <- c()

for (i in 1:length(unique(prediction$id))) {
    ev1 <- cbind(rgumbel(10000, loc = prediction$ev1[i], scale = 1), "ev1")
    ev2 <- cbind(rgumbel(10000, loc = prediction$ev2[i], scale = 1), "ev2")
    ev3 <- cbind(rep(NA, 10000), "ev3")
    
    choice <- as.data.frame(cbind(as.numeric(ev1[, 1]), as.numeric(ev2[, 1]), as.numeric(ev3[, 1])))
    choice$max <- max.col(choice[, 1:2])
    choice %>% rowwise() %>% mutate(., dce = max(V1, V2)) -> choice
    choice %>% rowwise() %>% mutate(., be = dce - sum(V1 + V2) + dce) -> choice
    choice %>% rowwise() %>% mutate(., be2 = dce - sum(V1 + V2) + dce) -> choice
    choice %>% group_by(max) %>% summarise(mean_dce = mean(dce), mean_be = mean(be), mean_be2 = mean(be2), n = n()) -> df
    df$subj <- prediction$id[i]
    df$condition <- prediction$condition[i]
    conf_dat <- rbind(conf_dat, df)
}

for (i in 11:nrow(prediction)) {
    ev1 <- cbind(rgumbel(10000, loc = prediction$ev1[i], scale = 1), "ev1")
    ev2 <- cbind(rgumbel(10000, loc = prediction$ev2[i], scale = 1), "ev2")
    ev3 <- cbind(rgumbel(10000, loc = prediction$ev3[i], scale = 1), "ev3")
    
    choice <- as.data.frame(cbind(as.numeric(ev1[, 1]), as.numeric(ev2[, 1]), as.numeric(ev3[, 1])))
    choice$max <- max.col(choice)
    choice %>% rowwise() %>% mutate(., dce = max(V1, V2, V3)) -> choice
    choice %>% rowwise() %>% mutate(., be = dce - (sum(V1 + V2 + V3) - dce)/2) -> choice
    choice %>% rowwise() %>% mutate(., be2 = dce - sum(V1 + V2 + V3) + dce + min(V1, V2, V3)) -> choice
    choice %>% group_by(max) %>% summarise(mean_dce = mean(dce), mean_be = mean(be), mean_be2 = mean(be2), n = n()) -> df
    df$subj <- prediction$id[i]
    df$condition <- prediction$condition[i]
    conf_dat <- rbind(conf_dat, df)
}


conf_dat <- mutate(conf_dat, max = ifelse(max == 1, "target", ifelse(max == 2, "distractor", "dud")))
conf_dat$max <- factor(conf_dat$max, levels = c("target", "distractor", "dud"))

ggplot(conf_dat) + geom_line(aes(x = condition, y = mean_be, color = subj)) + 
    stat_summary(fun = "mean", geom = "line", aes(x = condition, y = mean_be), size = 1) +
    facet_grid(. ~ max) + ggtitle("logit") -> logit_be
logit_be

ggplot(conf_dat) + geom_line(aes(x = condition, y = mean_be2, color = subj)) + 
    stat_summary(fun = "mean", geom = "line", aes(x = condition, y = mean_be2), size = 1) +
    facet_grid(. ~ max) + ggtitle("logit") -> logit_be2
logit_be2

ggplot(conf_dat) + geom_line(aes(x = condition, y = mean_dce, color = subj)) + 
    stat_summary(fun = "mean", geom = "line", aes(x = condition, y = mean_dce), size = 1) +
    facet_grid(. ~ max) + ggtitle("logit") -> logit_dce
logit_dce