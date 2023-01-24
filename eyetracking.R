#+ message = F
# 参加者間の比較・個人差 確信度と判断の一次元化
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
    group_by(condition, targetPos, dudPos, fixItem, subj) %>%
    summarise(n = n()) %>%
    ungroup() -> freq # ungroup() is necessary for plotting dud fix data

plot2 <- foreach(i = unique(freq$condition), .packages = c("tidyverse", "ggh4x")) %dopar% {
    p <- ggplot(subset(freq, freq$condition == i)) + 
        geom_violin(aes(x = fixItem, y = n, color = fixItem)) + 
        geom_point(aes(x = fixItem, y = n, color = fixItem)) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        facet_nested(. ~ targetPos + dudPos) + ggtitle(i)
    print(p)
}
plot2


#'# All trials
dat %>%
    group_by(subj, condition) %>%
    mutate(n_trials = n_distinct(trial), sum_fixations = n()) %>%
    group_by(fixItem, condition, subj) %>%
    mutate(n_fixations = n(), fpt = sum_fixations/n_trials, cfpt = n()/n_trials) %>%
    select(sum_fixations, n_fixations, n_trials, fpt, cfpt) %>%
    distinct() -> fd1

# total fixation frequency
ggplot(fd1) + geom_violin(aes(x = "", y = fpt)) + geom_point(aes(x = "", y = fpt, color = subj)) + 
    xlab("") + ylab("Mean fixations per trial") + facet_wrap(. ~ condition)

# condition-wise fixation frequency
ggplot(fd1) + geom_violin(aes(x = fixItem, y = cfpt)) + geom_point(aes(x = fixItem, y = cfpt)) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylab("Mean fixations per trial") + facet_wrap(. ~ condition)


#'# AOI only
# total fixation frequency (AOI only) other以外へのfixationがなかった試行は除かれる
dat %>%
    filter(., fixItem != "other") %>%
    group_by(subj, condition) %>%
    mutate(n_trials = n_distinct(trial), sum_fixations = n()) %>%
    group_by(fixItem, condition, subj) %>%
    mutate(n_fixations = n(), fpt = sum_fixations/n_trials, cfpt = n()/n_trials) %>%
    select(sum_fixations, n_fixations, n_trials, fpt, cfpt) %>%
    distinct() -> fd2

# total fixation frequency (AOI only)
ggplot(fd2) + geom_violin(aes(x = "", y = fpt)) + geom_point(aes(x = "", y = fpt, color = subj)) + 
    xlab("") + ylab("Mean fixations per trial (AOI only)") + facet_wrap(. ~ condition)

# condition-wise fixation frequency (AOI only)   
ggplot(fd2) + geom_violin(aes(x = fixItem, y = cfpt)) + geom_point(aes(x = fixItem, y = cfpt)) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylab("Mean fixation per trial (AOI only)") + facet_wrap(. ~ condition)


#'# Target and distractor only   
dat %>%
    filter(., fixItem != "other" & fixItem != "dud") %>%
    group_by(subj, condition) %>%
    mutate(n_trials = n_distinct(trial), sum_fixations = n()) %>%
    group_by(fixItem, condition, subj) %>%
    mutate(n_fixations = n(), fpt = sum_fixations/n_trials, cfpt = n()/n_trials) %>%
    select(sum_fixations, n_fixations, n_trials, fpt, cfpt) %>%
    distinct() -> fd3

# total fixation frequency 
ggplot(fd3) + geom_violin(aes(x = "", y = fpt)) + geom_point(aes(x = "", y = fpt)) + 
    xlab("") + ylab("Mean fixations per trial (target and distractor only)") + facet_wrap(. ~ condition)

# condition-wise fixation frequency (target and distractor only)   
ggplot(fd3) + geom_violin(aes(x = fixItem, y = cfpt)) + geom_point(aes(x = fixItem, y = cfpt)) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylab("Mean fixations per trial (target and distractor only)") + facet_wrap(. ~ condition)


#'# Stimulus-based fixation frequency (choice considered)
dat %>%
    group_by(subj, condition, chosenItem) %>%
    mutate(n_trials = n_distinct(trial), sum_fixations = n()) %>%
    group_by(fixItem, condition, subj, chosenItem) %>%
    mutate(n_fixations = n(), fpt = sum_fixations/n_trials, cfpt = n()/n_trials) %>%
    select(sum_fixations, n_fixations, n_trials, fpt, cfpt) %>%
    distinct() -> fd4

ggplot(fd4) + geom_violin(aes(x = chosenItem, y = fpt)) + geom_point(aes(x = chosenItem, y = fpt)) +
    ylab("Mean fixations per trial") + facet_wrap(. ~ condition)

ggplot(fd4) + geom_violin(aes(x = chosenItem, y = cfpt, color = fixItem)) + 
    geom_point(aes(x = chosenItem, y = cfpt, color = fixItem), position = position_dodge(width = 0.85)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylim(0, 2.5) + ylab("Mean fixations per trial") + facet_wrap(. ~ condition)


#'# Stimulus-based fixation frequency (confidence considered)
dat %>%
    distinct(subj, condition, conf, trial, nFix, nFix_target) %>%
    group_by(subj, condition, conf) %>%
    summarise(fpt = mean(nFix), tfpt = mean(nFix_target)) %>%
    ungroup(subj, condition, conf) -> fd5

# total fixations
ggplot(fd5) + geom_violin(aes(x = factor(conf), y = fpt)) + geom_point(aes(x = factor(conf), y = fpt)) +
    ylab("Mean fixations per trial") + facet_wrap(. ~ condition)

# total fixations
ggplot(fd5) + geom_violin(aes(x = factor(conf), y = tfpt)) + geom_point(aes(x = factor(conf), y = tfpt)) +
    ylab("Mean target fixations per trial") + facet_wrap(. ~ condition)


#'# Stimulus-based fixation frequency (choice and confidence considered)
dat %>%
    distinct(subj, condition, chosenItem, conf, trial, nFix, nFix_target) %>%
    group_by(subj, condition, chosenItem, conf) %>%
    summarise(fpt = mean(nFix), tfpt = mean(nFix_target)) %>%
    ungroup(subj, condition, chosenItem, conf) %>%
    complete(subj, condition, chosenItem, conf) -> fd6
fd6$fpt[is.na(fd6$fpt)] <- 0
fd6$tfpt[is.na(fd6$tfpt)] <- 0

# total fixations
ggplot(fd6) + geom_violin(aes(x = factor(conf), y = fpt, color = chosenItem)) + 
    geom_point(aes(x = factor(conf), y = fpt, color = chosenItem), position = position_dodge(width = 0.85)) +
    ylab("Mean fixations per trial") + facet_wrap(. ~ condition)

# total target fixations
ggplot(fd6) + geom_violin(aes(x = factor(conf), y = tfpt, color = chosenItem)) + 
    geom_point(aes(x = factor(conf), y = tfpt, color = chosenItem), position = position_dodge(width = 0.85)) +
    ylab("Mean target fixations per trial") + ylim(0, 2.5) + facet_wrap(. ~ condition)


#'# fixation dynamics
prop <- foreach(i = 1:7, .packages = c("tidyverse", "ggh4x")) %dopar% {
    dat %>% filter(event == "fixation" & countFix <= i) %>%
        group_by(subj, condition) %>%
        mutate(totalFix = n()) %>%
        group_by(subj, fixItem, condition) %>%
        mutate(fix = n(), pFix = n()/totalFix) %>%
        select(subj, fixItem, fix, totalFix, pFix, condition) -> df
    df$i <- i
    print(distinct(df))
}

p_dat <- c()
for (i in 1:7) {
    p_dat <- rbind(p_dat, prop[[i]])
}

p_dat %>%
    ungroup(subj, fixItem, i) %>%
    complete(subj, fixItem, i) -> p_dat
p_dat$fix[is.na(p_dat$fix)] <- 0
p_dat$pFix[is.na(p_dat$pFix)] <- 0

ggplot(p_dat, aes(x = i, y = pFix, color = fixItem)) + geom_point() +
    stat_summary(fun.y = "mean", geom = "line", position = position_dodge(width = .9)) +
    scale_x_continuous(breaks = seq(2, 7, 1), limits = c(1.5, 7.5))
    + ylab("Cumulative fixation proportion") + facet_wrap(. ~ condition)


#'# fixation dynamics (exclude other fixations)
prop <- foreach(i = 1:7, .packages = c("tidyverse", "ggh4x")) %dopar% {
    dat %>% filter(event == "fixation" & fixItem != "other" & countFix <= i) %>%
        group_by(subj, condition) %>%
        mutate(totalFix = n()) %>%
        group_by(subj, fixItem, condition) %>%
        mutate(fix = n(), pFix = n()/totalFix) %>%
        select(subj, fixItem, fix, totalFix, pFix, condition) -> df
    df$i <- i
    print(distinct(df))
}

p_dat <- c()
for (i in 1:7) {
    p_dat <- rbind(p_dat, prop[[i]])
}

p_dat %>%
    ungroup(subj, fixItem, i) %>%
    complete(subj, fixItem, i) -> p_dat
p_dat$fix[is.na(p_dat$fix)] <- 0
p_dat$pFix[is.na(p_dat$pFix)] <- 0
p_dat <- subset(p_dat, p_dat$fixItem != "other")

ggplot(p_dat, aes(x = i, y = pFix, color = fixItem)) + geom_point() +
    stat_summary(fun.y = "mean", geom = "line") +
    scale_x_continuous(breaks = seq(2, 7, 1), limits = c(1.5, 7.5)) + ylim(0, 0.7) + 
    xlab("countFix") + ylab("Cumulative fixation proportion") + facet_wrap(. ~ condition)


#'# fixation dynamics (exclude other and dud fixations)
prop <- foreach(i = 1:7, .packages = c("tidyverse", "ggh4x")) %dopar% {
    dat %>% filter(event == "fixation" & fixItem != "other" & fixItem != "dud" & countFix <= i) %>%
        group_by(subj, condition) %>%
        mutate(totalFix = n()) %>%
        group_by(subj, fixItem, condition) %>%
        mutate(fix = n(), pFix = n()/totalFix) %>%
        select(subj, fixItem, fix, totalFix, pFix, condition) -> df
    df$i <- i
    print(distinct(df))
}

p_dat <- c()
for (i in 1:7) {
    p_dat <- rbind(p_dat, prop[[i]])
}

p_dat %>%
    ungroup(subj, fixItem, i) %>%
    complete(subj, fixItem, i) -> p_dat
p_dat$fix[is.na(p_dat$fix)] <- 0
p_dat$pFix[is.na(p_dat$pFix)] <- 0
p_dat <- subset(p_dat, p_dat$fixItem != "other" & p_dat$fixItem != "dud")

ggplot(p_dat, aes(x = i, y = pFix, color = fixItem)) + geom_point() +
    stat_summary(fun.y = "mean", geom = "line") +
    scale_x_continuous(breaks = seq(2, 7, 1), limits = c(1.5, 7.5)) + ylim(0, 0.7) + 
    xlab("countFix") + ylab("Cumulative fixation proportion") + facet_wrap(. ~ condition)

ggplot(subset(p_dat, p_dat$i != 1), aes(x = as.numeric(as.character(condition)), y = pFix, color = fixItem)) + geom_point() +
    stat_summary(fun.y = "mean", geom = "line") +
    xlab("Condition") + ylab("Cumulative fixation proportion") + facet_wrap(. ~ i)