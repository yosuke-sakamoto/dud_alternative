#+ message = F
library(tidyverse)
library(data.table)
library(lme4)
library(lmerTest)
library(car)
library(multcomp)
library(sjPlot)

#'# data loading
d <- fread("exp1_data_behavior.csv", sep = ",", header = T)
d <- mutate(d, Corr = ifelse(ChosenITM == CorrectITM, 1, 0))
d$ChoiceRT <- as.numeric(d$ChoiceRT)
d <- na.omit(d) # 67 trials omitted
# d <- subset(d, d$participant != "sub01") # accuracy below mean - 2sd at Val3 = 63

d1 <- rbind(subset(d, d$UVal < 90 & d$ChosenITM == "up"),
            subset(d, d$LVal < 90 & d$ChosenITM == "left"),
            subset(d, d$RVal < 90 & d$ChosenITM == "right")) # dud chosen

d2 <- rbind(subset(d, d$UVal < 90 & d$ChosenITM != "up"),
            subset(d, d$LVal < 90 & d$ChosenITM != "left"),
            subset(d, d$RVal < 90 & d$ChosenITM != "right")) # dud not chosen

table(d1$Condition)
table(d2$Condition)

#'# type2 roc (aggregated)
roc2 <- c()
auc2 <- c()

for (i in unique(d2$Condition)) {
    
    dat <- subset(d2, d2$Condition == i)
    
    s1_c4 <- sum(dat$Corr == 1 & dat$Conf == 4)
    s1_c3 <- sum(dat$Corr == 1 & dat$Conf == 3)
    s1_c2 <- sum(dat$Corr == 1 & dat$Conf == 2)
    s1_c1 <- sum(dat$Corr == 1 & dat$Conf == 1)
    s2_c4 <- sum(dat$Corr == 0 & dat$Conf == 4)
    s2_c3 <- sum(dat$Corr == 0 & dat$Conf == 3)
    s2_c2 <- sum(dat$Corr == 0 & dat$Conf == 2)
    s2_c1 <- sum(dat$Corr == 0 & dat$Conf == 1)
    
    nr_s1 <- c(s1_c4, s1_c3, s1_c2, s1_c1)
    nr_s2 <- c(s2_c4, s2_c3, s2_c2, s2_c1)
    
    nhit <- cumsum(nr_s1)
    nfa <- cumsum(nr_s2)
    
    phit <- nhit/sum(nr_s1)
    pfa <- nfa/sum(nr_s2)
    
    roc <- data.frame(phit, pfa)[-4,]
    roc$Condition <- i
    roc2 <- rbind(roc2, roc)
    
    auc <- phit[1] * pfa[1]/2
    for (n in 1:3) {
        auc <- auc + (phit[n] + phit[n + 1])*(pfa[n + 1] - pfa[n])/2
    }
    auc <- data.frame(auc, i)
    auc2 <- rbind(auc2, auc)
}


p1 <- ggplot(roc2, aes(x = pfa, y = phit, color = factor(Condition))) + geom_point() + geom_line() +
    xlim(0, 1) + ylim(0, 1) + xlab("Type2 false alarm rate") + 
    ylab("Type2 hit rate") + ggtitle("Type-2 ROC")
p1

p2 <- ggplot(auc2, aes(x = factor(i), y = auc)) + geom_point() +
    ylim(0.4, 0.9) + xlab("Condition") + ylab("Type2 AUC") + ggtitle("Type-2 AUC")
p2


#'# type2 roc (individual)
iroc2 <- c()
iauc2 <- c()

for (i in unique(d2$participant)) {
    
    idat <- subset(d2, d2$participant == i)
    
    for (j in unique(idat$Condition)) {
        
        jdat <- subset(idat, idat$Condition == j)
        
        s1_c4 <- sum(jdat$Corr == 1 & jdat$Conf == 4)
        s1_c3 <- sum(jdat$Corr == 1 & jdat$Conf == 3)
        s1_c2 <- sum(jdat$Corr == 1 & jdat$Conf == 2)
        s1_c1 <- sum(jdat$Corr == 1 & jdat$Conf == 1)
        s2_c4 <- sum(jdat$Corr == 0 & jdat$Conf == 4)
        s2_c3 <- sum(jdat$Corr == 0 & jdat$Conf == 3)
        s2_c2 <- sum(jdat$Corr == 0 & jdat$Conf == 2)
        s2_c1 <- sum(jdat$Corr == 0 & jdat$Conf == 1)
        
        nr_s1 <- c(s1_c4, s1_c3, s1_c2, s1_c1)
        nr_s2 <- c(s2_c4, s2_c3, s2_c2, s2_c1)
        
        nhit <- cumsum(nr_s1)
        nfa <- cumsum(nr_s2)
        
        phit <- nhit/sum(nr_s1)
        pfa <- nfa/sum(nr_s2)
        
        roc <- data.frame(phit, pfa)[-4,]
        roc$id <- i
        roc$Condition <- j
        roc <- mutate(roc, order = row_number())
        iroc2 <- rbind(iroc2, roc)
        
        auc <- phit[1] * pfa[1]/2
        for (n in 1:3) {
            auc <- auc + (phit[n] + phit[n + 1])*(pfa[n + 1] - pfa[n])/2
        }
        auc <- data.frame(auc, i, j)
        iauc2 <- rbind(iauc2, auc)
    }
}
        
p3 <- ggplot(iroc2, aes(x = pfa, y = phit, color = factor(Condition))) + geom_point() + geom_line() +
    facet_wrap( ~ id) + xlim(0, 1) + ylim(0, 1) + xlab("Type2 false alarm rate") + 
    ylab("Type2 hit rate") + ggtitle("Type-2 ROC")
p3

p4 <- ggplot(iauc2, aes(x = factor(j), y = auc)) + geom_point() +
    facet_wrap(. ~ i) + ylim(0.4, 0.9) + xlab("Condition") + ylab("Type2 AUC") + ggtitle("Type-2 AUC")
p4

p5 <- iroc2 %>%
    group_by(Condition, order) %>%
    summarise(mhit = mean(phit), mfa = mean(pfa), sehit = sd(phit)/sqrt(10), sefa = sd(pfa)/sqrt(10)) %>%
    ggplot(., aes(x = mfa, y = mhit, color = factor(Condition))) + geom_point() + geom_line() +
    geom_errorbarh(aes(xmax = mfa + sefa, xmin = mfa - sefa)) + xlim(0, 1) + ylim(0, 1) +
    geom_errorbar(aes(ymax = mhit + sehit, ymin = mhit - sehit)) + xlab("Type2 false alarm rate") + 
    ylab("Type2 hit rate") + ggtitle("Type-2 ROC")
p5

p6 <- ggplot(iauc2, aes(x = factor(j), y = auc, color = factor(j))) + geom_point() +
    stat_summary(fun.y = "mean", geom = "crossbar", position = position_dodge(width = .9)) +
    ylim(0.4, 0.9) + xlab("Condition") + ylab("Type2 AUC") + ggtitle("Type-2 AUC") + guides(color = F)
p6