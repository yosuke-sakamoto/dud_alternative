#+ message = F
# 参加者間の比較・個人差 確信度と判断の一次元化
# sub03 almost always gives conf of 4

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
dat %>% group_by(subj) %>% mutate(conf_normalized = scale(conf)) -> dat # subject-wise normalization of confidence 
dat %>% mutate(q_dur_distractor = ifelse(dur_distractor > quantile(dur_distractor)[4], 0.75,
                                  ifelse(dur_distractor <= quantile(dur_distractor)[4] & dur_distractor > quantile(dur_distractor)[3], 0.5,
                                  ifelse(dur_distractor <= quantile(dur_distractor)[3] & dur_distractor > quantile(dur_distractor)[2], 0.25, 0)))) -> dat
dat %>% mutate(tdDurationRatio = dur_target/dur_distractor) -> dat
                                                              

#'# subject-wise fixation plot
plot1 <- foreach(i = unique(dat$subj), .packages = c("tidyverse", "ggh4x")) %dopar% {
    subset(dat, dat$subj == i) %>%
        ggplot() + geom_point(aes(x = x, y = y, size = dur, color = condition), alpha = 0.3) + 
        facet_nested(. ~ targetPos + dudPos) + ggtitle(i) %>% print()
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
        facet_nested(. ~ targetPos + dudPos) + ggtitle(i) %>% print()
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
p_dat <- foreach(i = 1:7, .combine = rbind, .packages = "tidyverse") %dopar% {
    dat %>% filter(event == "fixation" & countFix <= i) %>%
        group_by(subj, condition) %>%
        mutate(totalFix = n()) %>%
        group_by(subj, fixItem, condition) %>%
        mutate(fix = n(), pFix = n()/totalFix) %>%
        select(subj, fixItem, fix, totalFix, pFix, condition) -> df
    df$i <- i
    print(distinct(df))
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


#'# nFix_target, nFix_distractorの両者で反応正誤を説明
hist(dat$nFix_target)
hist(dat$nFix_distractor)
cor(dat$nFix_target, dat$nFix_distractor)

# condition aggregated
ggplot(subset(dat, dat$nFix_target <= 3 & dat$nFix_distractor <= 3),
       aes(x = nFix_target, , y = corr, color = factor(nFix_distractor))) + 
    geom_count(position = position_dodge(width = 0.3)) +
    stat_summary(fun.y = "mean", geom = "line") +
    scale_x_continuous(breaks = seq(0, 3, 1), limits = c(-0.5, 3.5))

f1 <- glm(corr ~ nFix_target * nFix_distractor, family = binomial, 
          data = subset(dat, dat$nFix_target <= 3 & dat$nFix_distractor <= 3))
summary(f1)
Anova(f1)
plot(ggpredict(f1, terms = c("nFix_target", "nFix_distractor")))

f2 <- glm(corr ~ nFix_target * factor(nFix_distractor), family = binomial, 
          data = subset(dat, dat$nFix_target <= 3 & dat$nFix_distractor <= 3))
summary(f2)
Anova(f2)
plot(ggpredict(f2, terms = c("nFix_target", "nFix_distractor")))


# condition separated
ggplot(subset(dat, dat$nFix_target <= 3 & dat$nFix_distractor <= 3),
       aes(x = nFix_target, , y = corr, color = factor(nFix_distractor))) + 
    geom_count(position = position_dodge(width = 1.2)) +
    stat_summary(fun.y = "mean", geom = "line") +
    scale_x_continuous(breaks = seq(0, 3, 1), limits = c(-0.5, 3.5)) + facet_wrap(. ~ condition)

f3 <- glm(corr ~ nFix_target * nFix_distractor * factor(condition), family = binomial, 
          data = subset(dat, dat$nFix_target <= 3 & dat$nFix_distractor <= 3))
summary(f3)
Anova(f3)
plot(ggpredict(f3, terms = c("nFix_target", "nFix_distractor", "condition")))

f4 <- glm(corr ~ nFix_target * factor(nFix_distractor) * factor(condition), family = binomial, 
          data = subset(dat, dat$nFix_target <= 3 & dat$nFix_distractor <= 3))
summary(f4)
Anova(f4)
plot(ggpredict(f4, terms = c("nFix_target", "nFix_distractor", "condition")))


#'# nFix_target, nFix_distractorの両者で確信度を説明
hist(dat$nFix_target)
hist(dat$nFix_distractor)
hist(dat$conf)

# condition aggregated
ggplot(subset(dat, dat$nFix_target <= 3 & dat$nFix_distractor <= 3),
       aes(x = nFix_target, , y = conf, color = factor(nFix_distractor))) + 
    geom_count(position = position_dodge(width = 0.3)) +
    stat_summary(fun.y = "mean", geom = "line") +
    scale_x_continuous(breaks = seq(0, 3, 1), limits = c(-0.5, 3.5)) + facet_wrap(. ~ chosenItem)
f1 <- glm(corr ~ nFix_target * nFix_distractor, family = binomial, 
          data = subset(dat, dat$nFix_target <= 3 & dat$nFix_distractor <= 3))
summary(f1)
Anova(f1)
plot(ggpredict(f1, terms = c("nFix_target", "nFix_distractor")))

f2 <- glm(corr ~ nFix_target * factor(nFix_distractor), family = binomial, 
          data = subset(dat, dat$nFix_target <= 3 & dat$nFix_distractor <= 3))
summary(f2)
Anova(f2)
plot(ggpredict(f2, terms = c("nFix_target", "nFix_distractor")))


# condition separated
ggplot(subset(dat, dat$nFix_target <= 3 & dat$nFix_distractor <= 3),
       aes(x = nFix_target, , y = corr, color = factor(nFix_distractor))) + 
    geom_count(position = position_dodge(width = 1.2)) +
    stat_summary(fun.y = "mean", geom = "line") +
    scale_x_continuous(breaks = seq(0, 3, 1), limits = c(-0.5, 3.5)) + facet_wrap(. ~ condition)

f3 <- glm(corr ~ nFix_target * nFix_distractor * factor(condition), family = binomial, 
          data = subset(dat, dat$nFix_target <= 3 & dat$nFix_distractor <= 3))
summary(f3)
Anova(f3)
plot(ggpredict(f3, terms = c("nFix_target", "nFix_distractor", "condition")))

f4 <- glm(corr ~ nFix_target * factor(nFix_distractor) * factor(condition), family = binomial, 
          data = subset(dat, dat$nFix_target <= 3 & dat$nFix_distractor <= 3))
summary(f4)
Anova(f4)
plot(ggpredict(f4, terms = c("nFix_target", "nFix_distractor", "condition")))


#'# nFix_target, nFix_distractorの両者で標準化された確信度を説明
hist(dat$nFix_target)
hist(dat$nFix_distractor)
hist(dat$conf_normalized)

# condition aggregated
ggplot(subset(dat, dat$nFix_target <= 3 & dat$nFix_distractor <= 3 & dat$subj != "sub03"), # sub03 almost always gives conf of 4
       aes(x = nFix_target, , y = conf_normalized, color = factor(nFix_distractor))) + 
    geom_count(position = position_dodge(width = 0.3)) +
    stat_summary(fun.y = "mean", geom = "line") +
    scale_x_continuous(breaks = seq(0, 3, 1), limits = c(-0.5, 3.5)) + 
    ylim(-3, 3) + facet_wrap(. ~ chosenItem)

f5 <- lm(conf_normalized ~ nFix_target * nFix_distractor * chosenItem, 
          data = subset(dat, dat$nFix_target <= 3 & dat$nFix_distractor <= 3 & 
                            dat$chosenItem != "dud" & dat$subj != "sub03"))
summary(f5)
Anova(f5)
plot(ggpredict(f5, terms = c("nFix_target", "nFix_distractor", "chosenItem")))

f6 <- lm(conf_normalized ~ nFix_target * factor(nFix_distractor) * chosenItem, 
         data = subset(dat, dat$nFix_target <= 3 & dat$nFix_distractor <= 3 & 
                           dat$chosenItem != "dud" & dat$subj != "sub03"))
summary(f6)
Anova(f6)
plot(ggpredict(f6, terms = c("nFix_target", "nFix_distractor", "chosenItem")))


# condition separated
ggplot(subset(dat, dat$nFix_target <= 3 & dat$nFix_distractor <= 3 & 
                  dat$chosenItem != "dud" & dat$subj != "sub03"), # sub03 almost always gives conf of 4
       aes(x = nFix_target, , y = conf_normalized, color = factor(nFix_distractor))) + 
    geom_count(position = position_dodge(width = 0.3)) +
    stat_summary(fun.y = "mean", geom = "line") +
    scale_x_continuous(breaks = seq(0, 3, 1), limits = c(-0.5, 3.5)) + 
    ylim(-3, 3) + facet_nested(. ~ chosenItem + condition)

f7 <- lm(conf_normalized ~ nFix_target * nFix_distractor * chosenItem * factor(condition), 
         data = subset(dat, dat$nFix_target <= 3 & dat$nFix_distractor <= 3 & 
                           dat$chosenItem != "dud" & dat$subj != "sub03"))
summary(f7)
Anova(f7)
plot(ggpredict(f7, terms = c("nFix_target", "nFix_distractor", "chosenItem", "condition")))

f8 <- lm(conf_normalized ~ nFix_target * factor(nFix_distractor) * chosenItem * factor(condition), 
         data = subset(dat, dat$nFix_target <= 3 & dat$nFix_distractor <= 3 & 
                           dat$chosenItem != "dud" & dat$subj != "sub03"))
summary(f8)
Anova(f8)
plot(ggpredict(f8, terms = c("nFix_target", "nFix_distractor", "chosenItem", "condition")))


#'# dur_target, dur_distractorの両者で反応正誤を説明
hist(dat$dur_target)
hist(dat$dur_distractor)
plot(dat$dur_target, dat$dur_distractor)
cor(dat$dur_target, dat$dur_distractor)

# condition aggregated
ggplot(dat, aes(x = dur_target, y = corr, color = factor(q_dur_distractor))) + 
    geom_count(alpha = 0.5) + stat_smooth() +
    scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1))
f9 <- glm(corr ~ dur_target * dur_distractor, family = binomial, data = dat)
summary(f9)
Anova(f9)
plot(ggpredict(f9, terms = c("dur_target", "dur_distractor")))

f10 <- glm(corr ~ dur_target * factor(q_dur_distractor), data = dat)
summary(f10)
Anova(f10)
plot(ggpredict(f10, terms = c("dur_target", "q_dur_distractor")))


# condition separated
ggplot(dat, aes(x = dur_target, y = corr, color = factor(q_dur_distractor))) + 
    geom_count(alpha = 0.5) + stat_smooth() +
    scale_x_continuous(breaks = seq(0, 0.6, 0.2), limits = c(0, 0.6)) + facet_wrap(. ~ condition)

f11 <- glm(corr ~ dur_target * dur_distractor * factor(condition), family = binomial, data = dat)
summary(f11)
Anova(f11)
plot(ggpredict(f11, terms = c("dur_target", "dur_distractor", "condition")))

f12 <- glm(corr ~ dur_target * factor(q_dur_distractor) * factor(condition), data = dat)
summary(f12)
Anova(f12)
plot(ggpredict(f12, terms = c("dur_target", "q_dur_distractor", "condition")))

#'# dur_target, dur_distractorの両者で標準化された確信度を説明
hist(dat$dur_target)
hist(dat$dur_distractor)
hist(dat$conf_normalized)

# condition aggregated
ggplot(subset(dat, dat$conf_normalized > -3), aes(x = dur_target, y = conf_normalized, color = factor(q_dur_distractor))) + 
    geom_count(alpha = 0.5) + stat_smooth(size = 1.2) +
    scale_x_continuous(breaks = seq(0, 0.6, 0.2), limits = c(0, 0.6)) + ylim(-3, 3) + facet_wrap(. ~ chosenItem)

f13 <- lm(conf_normalized ~ dur_target * dur_distractor * chosenItem, 
         data = subset(dat, dat$conf_normalized > -3 & 
                           dat$chosenItem != "dud" & dat$subj != "sub03"))
summary(f13)
Anova(f13)
plot(ggpredict(f13, terms = c("dur_target", "dur_distractor", "chosenItem")))

f14 <- lm(conf_normalized ~ dur_target * factor(q_dur_distractor) * chosenItem, 
         data = subset(dat, dat$conf_normalized > -3 & 
                           dat$chosenItem != "dud" & dat$subj != "sub03"))
summary(f14)
Anova(f14)
plot(ggpredict(f14, terms = c("dur_target", "q_dur_distractor", "chosenItem")))


# condition separated
ggplot(subset(dat, dat$conf_normalized > -3 & dat$chosenItem != "dud"), aes(x = dur_target, y = conf_normalized, color = factor(q_dur_distractor))) + 
    geom_count(alpha = 0.5) + stat_smooth(size = 1.2) +
    scale_x_continuous(breaks = seq(0, 0.6, 0.2), limits = c(0, 0.6)) + ylim(-3, 3) + facet_nested(. ~ chosenItem + condition)

f15 <- lm(conf_normalized ~ dur_target * dur_distractor * chosenItem * factor(condition), 
          data = subset(dat, dat$conf_normalized > -3 & 
                            dat$chosenItem != "dud" & dat$subj != "sub03"))
summary(f15)
Anova(f15)
plot(ggpredict(f15, terms = c("dur_target", "dur_distractor", "chosenItem", "condition")))

f16 <- lm(conf_normalized ~ dur_target * factor(q_dur_distractor) * chosenItem * factor(condition), 
          data = subset(dat, dat$conf_normalized > -3 & 
                            dat$chosenItem != "dud" & dat$subj != "sub03"))
summary(f16)
Anova(f16)
plot(ggpredict(f16, terms = c("dur_target", "q_dur_distractor", "chosenItem", "condition")))


#'# first fixation item
dat %>%
    group_by(subj, condition, chosenItem, firstFixItem) %>%
    summarise(n = n()) %>%
    ungroup(subj, condition, chosenItem) %>%
    complete(subj, condition, chosenItem) %>%
    ggplot(., aes(x = as.numeric(as.character(condition)), y = n, color = firstFixItem)) + 
    geom_point() + stat_summary(fun.y = "mean", geom = "line") + ggtitle("Number of first fixation")


#'# first fixation item (choice considered)
dat %>%
    group_by(subj, condition, chosenItem, firstFixItem) %>%
    summarise(n = n()) %>%
    ungroup(subj, condition, chosenItem) %>%
    complete(subj, condition, chosenItem) %>%
    ggplot(., aes(x = as.numeric(as.character(condition)), y = n, color = firstFixItem)) + 
    geom_point() + stat_summary(fun.y = "mean", geom = "line") + ggtitle("Number of first fixation") + facet_wrap(. ~ chosenItem)