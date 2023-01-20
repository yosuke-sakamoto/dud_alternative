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


#'# data loading 参加者間の比較・個人差
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
    mutate(., proportion = n_fixations/sum_fixations) %>%
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
    mutate(., proportion = n_fixations/sum_fixations) %>%
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
    mutate(., proportion = n_fixations/sum_fixations) %>%
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
    mutate(., proportion = n_fixations/sum_fixations) %>%
    distinct() -> fd4

ggplot(fd4) + geom_violin(aes(x = chosenItem, y = fpt)) + geom_point(aes(x = chosenItem, y = fpt)) +
    ylab("Mean fixations per trial") + facet_wrap(. ~ condition)

ggplot(fd4) + geom_violin(aes(x = chosenItem, y = cfpt, color = fixItem)) + 
    geom_point(aes(x = chosenItem, y = cfpt, color = fixItem), position = position_dodge(width = 0.85)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylim(0, 2.5) + ylab("Mean fixations per trial") + facet_wrap(. ~ condition)


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