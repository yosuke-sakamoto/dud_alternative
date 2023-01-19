library(evd)
library(tidyverse)

fits_logit
fits_dmn
fits_dnm2


#'# 
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
    prediction <- rbind(prediction, cbind(pred, evs))
}

colnames(prediction) <- c("c_pred", "ev1", "ev2", "ev3")


#'# gumbel sdt
ev1 <- cbind(rgumbel(10000, loc = prediction$ev1[1], scale = 1), "ev1")
ev2 <- cbind(rgumbel(10000, loc = prediction$ev2[1], scale = 1), "ev2")
ev3 <- cbind(rgumbel(10000, loc = prediction$ev3[1], scale = 1), "ev3")
g_sdt <- as.data.frame(rbind(ev1, ev2, ev3))
g_sdt$V1 <- as.numeric(g_sdt$V1)
g_sdt$V2 <- factor(g_sdt$V2, levels = c("ev1", "ev2", "ev3"))
ggplot(g_sdt) + geom_density(aes(x = V1, color = V2))

choice <- as.data.frame(cbind(as.numeric(ev1[, 1]), as.numeric(ev2[, 1]), as.numeric(ev3[, 1])))
choice$max <- max.col(choice)
table(choice$max)