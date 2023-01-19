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
dat <- fread("exp1_data_eyetracking.csv", header = TRUE)
dat <- subset(dat, dat$event == "fixation" & dat$conf != 0) # conf should be NA
dat$condition <- as.factor(dat$condition)
dat$fixItem <- factor(dat$fixItem, levels = c("target", "distractor", "dud", "other"))
dat$chosenItem <- factor(dat$chosenItem, levels = c("target", "distractor", "dud"))
bdat <- fread("exp1_data_behavior.csv", header = TRUE)


#'# subject-wise fixation plot
plot1 <- foreach(i = unique(dat$subj), .packages = c("tidyverse", "ggh4x")) %dopar% {
    subset(dat, dat$subj == i) %>%
        ggplot() + geom_point(aes(x = x, y = y, size = dur, color = condition), alpha = 0.3) + 
        facet_nested(. ~ targetPos + dudPos) + ggtitle(i) -> p
    print(p)
}
plot1


#'# Position-based fixation frequency
dat %>%
    group_by(condition, targetPos, dudPos, fixItem, id) %>%
    summarise(n = n()) -> freq

plot2 <- foreach(i = unique(freq$Condition), .packages = c("tidyverse", "ggh4x")) %dopar% {
    p <- ggplot(subset(freq, freq$Condition == i)) + geom_violin(aes(x = item, y = n, color = item)) + geom_point(aes(x = item, y = n)) + 
        facet_nested(. ~ target + dud) + 
        scale_x_discrete(limits = c("up", "left", "right", "noFix")) +
        scale_color_discrete(limits = c("up", "left", "right", "noFix")) + ggtitle(i)
    print(p)
}
plot2


#'# total fixation frequency
dat %>%
    group_by(condition, subj) %>%
    summarise(n = n(), fpt = n()/288) -> ff1
ggplot(ff1) + geom_violin(aes(x = "", y = fpt)) + geom_point(aes(x = "", y = fpt, color = subj)) + 
    xlab("") + ylab("Mean fixations per trial") + facet_wrap(. ~ condition)

#'# condition-wise fixation frequency
dat %>%
    group_by(fixItem, condition, subj) %>%
    summarise(n = n(), fpt = n()/288) -> ff2
ggplot(ff2) + geom_violin(aes(x = fixItem, y = fpt)) + geom_point(aes(x = fixItem, y = fpt)) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylab("Mean fixations per trial") + facet_wrap(. ~ condition)

#'# condition-wise fixation proportion (試行ワイズでカウントしていない)
ff2$totalFixWithinCondition <- c(ff1$n, ff1$n, ff1$n[11:60], ff1$n)
ff2$pFix <- ff2$n/ff2$totalFixWithinCondition 
ggplot(ff2) + geom_violin(aes(x = fixItem, y = pFix)) + geom_point(aes(x = fixItem, y = pFix)) + 
        theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylab("Total fixation proportion") + facet_wrap(. ~ condition)    

 
#'# total fixation frequency (AOI only) 参加者間の比較・個人差
dat %>%
    filter(., fixItem != "other") %>%
    group_by(condition, subj) %>%
    summarise(n = n(), fpt = n()/288) -> ff3
ggplot(ff3) + geom_violin(aes(x = "", y = fpt)) + geom_point(aes(x = "", y = fpt, color = subj)) + 
    xlab("") + ylab("Mean fixations per trial (AOI only)") + facet_wrap(. ~ condition)

#'# condition-wise fixation frequency (AOI only)   
dat %>%
    filter(., fixItem != "other") %>%
    group_by(fixItem, condition, subj) %>%
    summarise(n = n(), fpt = n()/288) -> ff4
ggplot(ff4) + geom_violin(aes(x = fixItem, y = fpt)) + geom_point(aes(x = fixItem, y = fpt)) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylab("Mean fixation per trial (AOI only)") + facet_wrap(. ~ condition)

#'# condition-wise fixation proportion (AOI only)   
ff4$totalFixWithinCondition <- c(ff3$n, ff1$n, ff3$n[11:60])
ff4$pFix <- ff4$n/ff4$totalFixWithinCondition 
ggplot(ff4) + geom_violin(aes(x = fixItem, y = pFix)) + geom_point(aes(x = fixItem, y = pFix)) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylab("Total fixation proportion (AOI only)") + facet_wrap(. ~ condition)    


#'# total fixation frequency (target and distractor only)   
dat %>%
    filter(., fixItem != "other" & fixItem != "dud") %>%
    group_by(condition, subj) %>%
    summarise(n = n(), fpt = n()/288) -> ff5
ggplot(ff5) + geom_violin(aes(x = "", y = fpt)) + geom_point(aes(x = "", y = fpt)) + 
    xlab("") + ylab("fixation frequency (Mean fixations per trial (target and distractor only)") + facet_wrap(. ~ condition)

#'# condition-wise fixation frequency (target and distractor only)   
dat %>%
    filter(., fixItem != "other" & fixItem != "dud") %>%
    group_by(fixItem, condition, subj) %>%
    summarise(n = n(), fpt = n()/288) -> ff6
ggplot(ff6) + geom_violin(aes(x = fixItem, y = fpt)) + geom_point(aes(x = fixItem, y = fpt)) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylab("Mean fixations per trial (target and distractor only)") + facet_wrap(. ~ condition)

#'# condition-wise fixation proportion (target and distractor only)   
ff6$totalFixWithinCondition <- c(ff5$n, ff5$n)
ff6$pFix <- ff6$n/ff6$totalFixWithinCondition 
ggplot(ff6) + geom_violin(aes(x = fixItem, y = pFix)) + geom_point(aes(x = fixItem, y = pFix)) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylab("Relative fixation proportion (target and distractor only)") + facet_wrap(. ~ condition)   


#'# Stimulus-based fixation frequency (choice considered)
dat %>%
    group_by(condition, chosenItem, subj) %>%
    summarise(n = n()) -> ff7
bdat %>%
    group_by(condition, chosenItem, subj) %>%
    summarise(n = n()) -> bf7
ff7 %>%
    ungroup() %>% complete(condition, chosenItem, subj) %>%
    ggplot() + geom_violin(aes(x = chosenItem, y = n)) + geom_point(aes(x = chosenItem, y = n)) +
    ylab("fixation frequency") + facet_wrap(. ~ condition)

dat %>%
    group_by(condition, fixItem, chosenItem, subj) %>%
    summarise(n = n()) -> ff8    
ff8 %>%
    ungroup() %>% complete(condition, fixItem, chosenItem, subj) %>%
    ggplot() + geom_violin(aes(x = fixItem, y = n)) + geom_point(aes(x = fixItem, y = n)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylab("fixation frequency") + facet_nested(. ~ chosenItem + condition)



#'# Stimulus-based fixation frequency (confidence considered)
dat %>%
    group_by(Condition, fix, choice, Conf, id) %>% 
    summarise(n = n()) %>% ungroup() %>% complete(Condition, fix, Conf) -> fixation3

fixation3$fix <- as.factor(fixation3$fix)
fixation3$fix <- factor(fixation3$fix, levels = c("target", "distractor", "dud", "noFix"))

# correct choice
subset(fixation3, fixation3$choice == "correct") %>%
    ungroup() %>% complete(Condition, fix, choice) %>%
    ggplot() + geom_point(position = position_dodge(width = .8), aes(x = fix, y = n, color = Condition)) +
    facet_grid(Conf ~ .) +
    stat_summary(fun.y = "mean", geom = "crossbar", position = position_dodge(width = .8), 
                 mapping = aes(x = fix, y = n, color = Condition)) + ggtitle("fixation frequency in correct trials")


#'# Stimulus-based fixation frequency (choice considered)
dat %>%
    group_by(Condition, fix, choice, id) %>%
    summarise(n = n()) -> fixation2

fixation2$fix <- as.factor(fixation2$fix)
fixation2$fix <- factor(fixation2$fix, levels = c("target", "distractor", "dud", "noFix"))

# correct choice
subset(fixation2, fixation2$choice == "correct") %>%
    ungroup() %>% complete(Condition, fix, choice) %>%
    ggplot() + geom_point(position = position_dodge(width = .8), aes(x = fix, y = n, color = Condition)) +
    stat_summary(fun.y = "mean", geom = "crossbar", position = position_dodge(width = .8), 
                 mapping = aes(x = fix, y = n, color = Condition)) + ggtitle("fixation frequency in correct trials")

# distractor choice
subset(fixation2, fixation2$choice == "distractor") %>%
    ungroup() %>% complete(Condition, fix, choice) %>%
    ggplot() + geom_point(position = position_dodge(width = .8), aes(x = fix, y = n, color = Condition)) +
    stat_summary(fun.y = "mean", geom = "crossbar", position = position_dodge(width = .8), 
                 mapping = aes(x = fix, y = n, color = Condition)) + ggtitle("fixation frequency in distractor-chosen trials")

# dud choice
subset(fixation2, fixation2$choice == "dud") %>%
    ungroup() %>% complete(Condition, fix, choice) %>%
    ggplot() + geom_point(position = position_dodge(width = .8), aes(x = fix, y = n, color = Condition)) +
    stat_summary(fun.y = "mean", geom = "crossbar", position = position_dodge(width = .8), 
                 mapping = aes(x = fix, y = n, color = Condition)) + ggtitle("fixation frequency in dud-chosen trials")



#'# Stimulus-based fixation proportion
fixation_prop1 <- c()

for (i in unique(fixation$id)) {
    for (j in unique(fixation$Condition)) {
        d <- subset(fixation,  fixation$id == i & fixation$Condition == j)
        s <- sum(d$n)
        fixation_prop1 <- rbind(fixation_prop1, mutate(d, prop = n/s))
    }
}

ggplot(fixation_prop1, aes(x = fix, y = prop, color = Condition)) + geom_violin() + 
    geom_point(position = position_dodge(width = .9)) +
    stat_summary(fun.y = "mean", geom = "crossbar", position = position_dodge(width = .9)) 


# noFix excluded)
fixation2 <- subset(fixation, fixation$fix != "noFix")
fixation_prop2 <- c()

for (i in unique(fixation2$id)) {
    for (j in unique(fixation2$Condition)) {
        d <- subset(fixation2,  fixation2$id == i & fixation2$Condition == j)
        s <- sum(d$n)
        fixation_prop2 <- rbind(fixation_prop2, mutate(d, prop = n/s))
    }
}

ggplot(fixation_prop2, aes(x = fix, y = prop, color = Condition)) + geom_violin() + 
    geom_point(position = position_dodge(width = .9)) +
    stat_summary(fun.y = "mean", geom = "crossbar", position = position_dodge(width = .9)) 

# dud excluded
fixation3 <- subset(fixation, fixation$fix == "target" | fixation$fix == "distractor")
fixation_prop3 <- c()

for (i in unique(fixation3$id)) {
    for (j in unique(fixation3$Condition)) {
        d <- subset(fixation3,  fixation3$id == i & fixation3$Condition == j)
        s <- sum(d$n)
        fixation_prop3 <- rbind(fixation_prop3, mutate(d, prop = n/s))
    }
}

ggplot(fixation_prop3, aes(x = fix, y = prop, color = Condition)) + geom_violin() + 
    geom_point(position = position_dodge(width = .9)) +
    stat_summary(fun.y = "mean", geom = "crossbar", position = position_dodge(width = .9)) 


#'# fixation dynamics
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


dat2 %>%
    group_by(id, fix, nFix) %>%
    summarise(n = n()) -> fixs

fix_prop <- c()
for (i in unique(fixs$id)) {
    for (j in 1:20) {
        fixs %>% filter(, id == i & nFix == j) -> d
        s <- sum(d$n)
        fix_prop <- rbind(fix_prop, mutate(d, prop = n/s))
    }
}

ggplot(fix_prop, aes(x = nFix, y = prop, color = fix)) + geom_point(position = position_dodge(width = .9)) +
    stat_summary(fun.y = "mean", geom = "line", position = position_dodge(width = .9))


fix_prop2 <- c()
for (i in unique(fixs$id)) {
    for (j in 1:20) {
        fixs %>% filter(, id == i & nFix == j & fix != "noFix") -> d
        s <- sum(d$n)
        fix_prop2 <- rbind(fix_prop2, mutate(d, prop = n/s))
    }
}

ggplot(fix_prop2, aes(x = nFix, y = prop, color = fix)) + geom_point(position = position_dodge(width = .9)) +
    stat_summary(fun.y = "mean", geom = "line", position = position_dodge(width = .9))