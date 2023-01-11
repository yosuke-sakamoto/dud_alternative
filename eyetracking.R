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
dat <- mutate(dat, target = apply(dat, 1, function(x){names(d)[13 + which.max(x[14:16])]}))
dat <- mutate(dat, dud = apply(dat, 1, function(x){names(d)[13 + which.min(x[14:16])]}))
dat$Condition <- as.factor(dat$Condition)
dat$target <- fct_recode(dat$target, up = "UVal", left = "LVal", right = "RVal")
dat$dud <- fct_recode(dat$dud, up = "UVal", left = "LVal", right = "RVal")
dat <- mutate(dat, fix = ifelse(item == target, "target", 
                         ifelse(item == dud, "dud",
                         ifelse(item == "noFix", "noFix", "distractor"))))
dat <- mutate(dat, choice = ifelse(ChosenITM == target, "correct", 
                            ifelse(ChosenITM == dud, "dud", "distractor")))


#'# subject-wise fixation plot
plot1 <- foreach(i = 1:length(files), .packages = c("tidyverse", "ggh4x")) %dopar% {
    subset(dat, dat$id == i) %>%
    ggplot() + geom_point(aes(x = x, y = y, size = dur, color = Condition), alpha = 0.3) + 
        facet_nested(. ~ target + dud) + ggtitle(i) -> p
    print(p)
}
plot1


#'# Position-based fixation frequency
dat %>%
    group_by(Condition, target, dud, item, id) %>%
    summarise(n = n()) -> freq

plot2 <- foreach(i = unique(freq$Condition), .packages = c("tidyverse", "ggh4x")) %dopar% {
    p <- ggplot(subset(freq, freq$Condition == i)) + geom_violin(aes(x = item, y = n, color = item)) + geom_point(aes(x = item, y = n)) + 
        facet_nested(. ~ target + dud) + 
        scale_x_discrete(limits = c("up", "left", "right", "noFix")) +
        scale_color_discrete(limits = c("up", "left", "right", "noFix")) + ggtitle(i)
    print(p)
}
plot2


#'# Stimulus-based fixation frequency
dat %>%
    group_by(Condition, fix, id) %>%
    summarise(n = n()) -> fixation

fixation$fix <- as.factor(fixation$fix)
fixation$fix <- factor(fixation$fix, levels = c("target", "distractor", "dud", "noFix"))

ggplot(fixation, aes(x = fix, y = n, color = Condition)) + geom_violin() + 
    geom_point(position = position_dodge(width = .9)) +
    stat_summary(fun.y = "mean", geom = "crossbar", position = position_dodge(width = .9)) 

f1 <- lmer(n ~ fix * Condition + (1 + fix + Condition | id), data = fixation)
summary(f1)
Anova(f1)
plot(ggpredict(f1, terms = c("fix", "Condition"), type = "fe"))


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