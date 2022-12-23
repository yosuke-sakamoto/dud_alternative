#+ message = F
library(tidyverse)

#' 多項ロジットモデル (P1 + P2 + P3 = 1)  

simdat <- c()
sigma <- c(0, 0.2) # semisaturation parameter in the divisive normalization model
omega <- c(-0.05, 0.05)
Val3 <- c(0, 27, 45, 63, 76, 86)

for (sH in sigma) {
    for (w in omega) {
        
        p1 <- c(
            exp(1 / (sH + w * 1.90)) / ( exp(1 / (sH + w * 1.90)) + exp(0.9 / (sH + w * 1.90)) + exp(0 / (sH + w * 1.90)) ),
            exp(1 / (sH + w * 2.17)) / ( exp(1 / (sH + w * 2.17)) + exp(0.9 / (sH + w * 2.17)) + exp(0.27 / (sH + w * 2.17)) ),
            exp(1 / (sH + w * 2.35)) / ( exp(1 / (sH + w * 2.35)) + exp(0.9 / (sH + w * 2.35)) + exp(0.45 / (sH + w * 2.35)) ),
            exp(1 / (sH + w * 2.53)) / ( exp(1 / (sH + w * 2.53)) + exp(0.9 / (sH + w * 2.53)) + exp(0.63 / (sH + w * 2.53)) ),
            exp(1 / (sH + w * 2.66)) / ( exp(1 / (sH + w * 2.66)) + exp(0.9 / (sH + w * 2.66)) + exp(0.76 / (sH + w * 2.66)) ),
            exp(1 / (sH + w * 2.76)) / ( exp(1 / (sH + w * 2.76)) + exp(0.9 / (sH + w * 2.76)) + exp(0.86 / (sH + w * 2.76)) ))
        
        p2 <- c(
            exp(0.9 / (sH + w * 1.90)) / ( exp(1 / (sH + w * 1.90)) + exp(0.9 / (sH + w * 1.90)) + exp(0 / (sH + w * 1.90)) ),
            exp(0.9 / (sH + w * 2.17)) / ( exp(1 / (sH + w * 2.17)) + exp(0.9 / (sH + w * 2.17)) + exp(0.27 / (sH + w * 2.17)) ),
            exp(0.9 / (sH + w * 2.35)) / ( exp(1 / (sH + w * 2.35)) + exp(0.9 / (sH + w * 2.35)) + exp(0.45 / (sH + w * 2.35)) ),
            exp(0.9 / (sH + w * 2.53)) / ( exp(1 / (sH + w * 2.53)) + exp(0.9 / (sH + w * 2.53)) + exp(0.63 / (sH + w * 2.53)) ),
            exp(0.9 / (sH + w * 2.66)) / ( exp(1 / (sH + w * 2.66)) + exp(0.9 / (sH + w * 2.66)) + exp(0.76 / (sH + w * 2.66)) ),
            exp(0.9 / (sH + w * 2.76)) / ( exp(1 / (sH + w * 2.76)) + exp(0.9 / (sH + w * 2.76)) + exp(0.86 / (sH + w * 2.76)) ))
        
        p3 <- c(
            exp(0    / (sH + w * 1.90)) / ( exp(1 / (sH + w * 1.90)) + exp(0.9 / (sH + w * 1.90)) + exp(0 / (sH + w * 1.90)) ),
            exp(0.27 / (sH + w * 2.17)) / ( exp(1 / (sH + w * 2.17)) + exp(0.9 / (sH + w * 2.17)) + exp(0.27 / (sH + w * 2.17)) ),
            exp(0.45 / (sH + w * 2.35)) / ( exp(1 / (sH + w * 2.35)) + exp(0.9 / (sH + w * 2.35)) + exp(0.45 / (sH + w * 2.35)) ),
            exp(0.63 / (sH + w * 2.53)) / ( exp(1 / (sH + w * 2.53)) + exp(0.9 / (sH + w * 2.53)) + exp(0.63 / (sH + w * 2.53)) ),
            exp(0.76 / (sH + w * 2.66)) / ( exp(1 / (sH + w * 2.66)) + exp(0.9 / (sH + w * 2.66)) + exp(0.76 / (sH + w * 2.66)) ),
            exp(0.86 / (sH + w * 2.76)) / ( exp(1 / (sH + w * 2.76)) + exp(0.9 / (sH + w * 2.76)) + exp(0.86 / (sH + w * 2.76)) ))
        
        sim <- as.data.frame(cbind(p1, p2, p3, Val3, sH, w))
        sim <- mutate(sim, sum = p1 + p2 + p3)
        sim <- mutate(sim, relacc = p1 / (p1 + p2))
        sim <- mutate(sim, odds_p1p2 = p1/p2)
        simdat <- rbind(simdat, sim)
    }
}

ggplot(simdat) + geom_line(aes(x = Val3, y = p1, linetype = factor(sH), color = factor(w))) + ylab("P1 / (P1 + P2 + P3)") + ylim(0, 1)
ggplot(simdat) + geom_line(aes(x = Val3, y = relacc, linetype = factor(sH), color = factor(w))) + ylab("P1 / (P1 + P2)") + ylim(0, 1)
ggplot(simdat) + geom_line(aes(x = Val3, y = odds_p1p2, linetype = factor(sH), color = factor(w))) + ylab("P1 / P2")