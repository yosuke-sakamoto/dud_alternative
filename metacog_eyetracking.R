#+ message = F
library(tidyverse)
library(data.table)
library(ggh4x)
library(lme4)
library(car)
library(ggeffects)
library(doParallel)

cores <- getOption("mc.cores", detectCores()) # for parallel computation
cl <- makeCluster(cores)
registerDoParallel(cl)


#'# data loading
files <- dir("exp1_data_eyetracking_individual_saccades", ".*csv$", full.names = TRUE) 
dat <- c()

for (i in 1:length(files)) {
    d <- fread(files[i], header = TRUE)
    d$id <- i  # numeric variable for subject id
    dat <- rbind(dat, d)
}

dat <- subset(dat, dat$event == "fixation" & dat$Conf != 0)
dat <- mutate(dat, target = apply(dat, 1, function(x){names(d)[14 + which.max(x[15:17])]}))
dat <- mutate(dat, dud = apply(dat, 1, function(x){names(d)[14 + which.min(x[15:17])]}))
dat$target <- fct_recode(dat$target, up = "UVal", left = "LVal", right = "RVal")
dat$dud <- fct_recode(dat$dud, up = "UVal", left = "LVal", right = "RVal")
dat <- mutate(dat, fix = ifelse(item == target, "target", 
                                ifelse(item == dud, "dud",
                                       ifelse(item == "other", "other", "distractor"))))
dat <- mutate(dat, choice = ifelse(ChosenITM == target, "correct", 
                                   ifelse(ChosenITM == dud, "dud", "distractor")))


#'# assign nFix within a trial
df <- foreach(i = unique(dat$id), .packages = "tidyverse") %dopar% {
    df1 <- c()
    df2 <- subset(dat, dat$id == i)
    for (j in unique(df2$trial)) {
        df1 <- bind_rows(df1, df2 %>% filter(trial == j) %>% mutate(, nFix = row_number()))
    }
    print(df1)
}
dat2 <- c()
for (i in unique(dat$id)) {
    dat2 <- rbind(dat2, df[[i]])
}


#'# proportion of target fixation frequency within a trial
# dat2 <- subset(dat2, dat2$nFix != 1) # first fixation is almost always on the screen center
# dat2 <- subset(dat2, dat2$nFix < 5) # subsetting trials
# dat2 <- subset(dat2, dat2$fromTrialBegin < 1)

df <- foreach(i = unique(dat2$id), .packages = "tidyverse") %dopar% {
    df1 <- c()
    df2 <- subset(dat2, dat2$id == i)
    for (j in unique(df2$trial)) {
        df2 %>% filter(, trial == j) -> d
        df1 <- rbind(df1, cbind(nrow(d),
                                nrow(subset(d, d$fix == "target")),
                                nrow(subset(d, d$fix == "target"))/nrow(d),
                                nrow(subset(d, d$fix == "target"))/nrow(subset(d, d$fix == "target" | d$fix == "distractor")),
                                d$Condition, d$Conf, d$id, d$choice)[1, ])
    }
    print(df1)
}

dat3 <- c()
for (i in 1:length(unique(dat$id))) {
    dat3 <- rbind(dat3, df[[i]])
}

dat3 <- as.data.frame(dat3)
colnames(dat3) <- c("n_fix", "n_target_fix", "p_target_fix", "p_target_fix2", "Condition", "Conf", "id", "choice")
dat3$id <- as.numeric(dat3$id)


#'# type2 roc (aggregated)
dat3$n_fix <- as.numeric(dat3$n_fix)
dat3$n_target_fix <- as.numeric(dat3$n_target_fix)
dat3$p_target_fix <- as.numeric(dat3$p_target_fix)
dat3$p_target_fix2 <- as.numeric(dat3$p_target_fix2)
dat3$Conf <- as.numeric(dat3$Conf)

q <- quantile(dat3$p_target_fix)
roc2 <- c()
auc2 <- c()

for (i in unique(dat3$Condition)) {
    
    df <- subset(dat3, dat3$Condition == i)
    
    s1_c4 <- sum(df$choice == "correct" & df$p_target_fix > q[4])
    s1_c3 <- sum(df$choice == "correct" & df$p_target_fix > q[3] & df$p_target_fix <= q[4])
    s1_c2 <- sum(df$choice == "correct" & df$p_target_fix > q[2] & df$p_target_fix <= q[3])
    s1_c1 <- sum(df$choice == "correct" & df$p_target_fix <= q[2])
    
    s2_c4 <- sum(df$choice != "correct" & df$p_target_fix > q[4])
    s2_c3 <- sum(df$choice != "correct" & df$p_target_fix > q[3] & df$p_target_fix <= q[4])
    s2_c2 <- sum(df$choice != "correct" & df$p_target_fix > q[2] & df$p_target_fix <= q[3])
    s2_c1 <- sum(df$choice != "correct" & df$p_target_fix <= q[2])
    
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
    ylab("Type2 hit rate") + ggtitle("Type-2 ROC (p_target_fix)")
p1

p2 <- ggplot(auc2, aes(x = factor(i), y = auc)) + geom_point() +
    ylim(0.4, 0.9) + xlab("Condition") + ylab("Type2 AUC") + ggtitle("Type-2 AUC (p_target_fix)")
p2


#'# type2 roc (individual)
iroc2 <- c()
iauc2 <- c()

for (i in unique(dat3$id)) {
    
    idat <- subset(dat3, dat3$id == i)
    q <- quantile(idat$p_target_fix)
    
    for (j in unique(idat$Condition)) {
        
        jdat <- subset(idat, idat$Condition == j)
        
        s1_c4 <- sum(jdat$choice == "correct" & jdat$p_target_fix > q[4])
        s1_c3 <- sum(jdat$choice == "correct" & jdat$p_target_fix > q[3] & jdat$p_target_fix <= q[4])
        s1_c2 <- sum(jdat$choice == "correct" & jdat$p_target_fix > q[2] & jdat$p_target_fix <= q[3])
        s1_c1 <- sum(jdat$choice == "correct" & jdat$p_target_fix <= q[2])
        
        s2_c4 <- sum(jdat$choice != "correct" & jdat$p_target_fix > q[4])
        s2_c3 <- sum(jdat$choice != "correct" & jdat$p_target_fix > q[3] & jdat$p_target_fix <= q[4])
        s2_c2 <- sum(jdat$choice != "correct" & jdat$p_target_fix > q[2] & jdat$p_target_fix <= q[3])
        s2_c1 <- sum(jdat$choice != "correct" & jdat$p_target_fix <= q[2])
        
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
    ylab("Type2 hit rate") + ggtitle("Type-2 ROC (p_target_fix)")
p3

p4 <- ggplot(iauc2, aes(x = factor(j), y = auc)) + geom_point() +
    facet_wrap(. ~ i) + ylim(0.4, 0.9) + xlab("Condition") + ylab("Type2 AUC") + ggtitle("Type-2 AUC (p_target_fix)")
p4

p5 <- iroc2 %>%
    group_by(Condition, order) %>%
    summarise(mhit = mean(phit), mfa = mean(pfa), sehit = sd(phit)/sqrt(10), sefa = sd(pfa)/sqrt(10)) %>%
    ggplot(., aes(x = mfa, y = mhit, color = factor(Condition))) + geom_point() + geom_line() +
    geom_errorbarh(aes(xmax = mfa + sefa, xmin = mfa - sefa)) + xlim(0, 1) + ylim(0, 1) +
    geom_errorbar(aes(ymax = mhit + sehit, ymin = mhit - sehit)) + xlab("Type2 false alarm rate") + 
    ylab("Type2 hit rate") + ggtitle("Type-2 ROC (p_target_fix)")
p5

p6 <- ggplot(iauc2, aes(x = factor(j), y = auc, color = factor(j))) + geom_point() +
    stat_summary(fun.y = "mean", geom = "crossbar", position = position_dodge(width = .9)) +
    ylim(0.4, 0.9) + xlab("Condition") + ylab("Type2 AUC") + ggtitle("Type-2 AUC (p_target_fix)") + guides(color = F)
p6



#'# correlation with confidence
dat4 <- subset(dat3, dat3$choice != "dud")

p7 <- ggplot(dat4, aes(x = n_fix, y = Conf, color = factor(choice))) + geom_jitter(alpha = 0.2, height = 0.2) +
    stat_smooth(method = "lm", size = 1.5) + 
    facet_wrap( ~ Condition) + ylab("Confidence")
p7

p8 <- ggplot(dat4, aes(x = n_target_fix, y = Conf, color = factor(choice))) + geom_jitter(alpha = 0.2, height = 0.2) +
    stat_smooth(method = "lm", size = 1.5) + 
    facet_wrap( ~ Condition) + ylab("Confidence")
p8

p9 <- ggplot(dat4, aes(x = p_target_fix, y = Conf, color = factor(choice))) + geom_jitter(alpha = 0.2, height = 0.2) +
    stat_smooth(method = "lm", size = 1.5) + 
    facet_wrap( ~ Condition) + ylab("Confidence")
p9