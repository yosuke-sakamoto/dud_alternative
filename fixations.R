#+ message = F
library(tidyverse)
library(data.table)
library(ggh4x)
library(lme4)
library(car)
library(ggeffects)


#'# data loading
files <- dir("exp1_data_eyetracking_individual_saccades", ".*csv$", full.names = TRUE) 
dat <- c()

for (i in 1:length(files)) {
    d <- fread(files[i], header = TRUE)
    d$id <- i  # numeric variable for subject id
    dat <- rbind(dat, d)
}

dat <- subset(dat, dat$event == "fixation")
dat <- mutate(dat, target = apply(dat, 1, function(x){names(d)[13 + which.max(x[14:16])]}))
dat <- mutate(dat, dud = apply(dat, 1, function(x){names(d)[13 + which.min(x[14:16])]}))
dat$Condition <- as.factor(dat$Condition)
dat$target <- fct_recode(dat$target, up = "UVal", left = "LVal", right = "RVal")
dat$dud <- fct_recode(dat$dud, up = "UVal", left = "LVal", right = "RVal")
dat <- mutate(dat, fix = ifelse(item == target, "target", 
                         ifelse(item == dud, "dud",
                         ifelse(item == "noFix", "noFix", "distractor"))))


#'# subject-wise fixation plot
for (i in 1:length(files)) {
    d <- subset(dat, dat$id == i)
    p <- ggplot(d) + geom_point(aes(x = x, y = y, size = dur, color = Condition), alpha = 0.3) + 
    facet_nested(. ~ target + dud) + ggtitle(i)
    print(p)
}


#'# subject-wise fixation frequency
dat %>%
    group_by(Condition, target, dud, item, id) %>%
    summarise(n = n()) -> freq

for (i in dat$Condition) {
    p <- ggplot(subset(freq, freq$Condition == i)) + geom_violin(aes(x = item, y = n, color = item)) + geom_point(aes(x = item, y = n)) + 
        facet_nested(. ~ target + dud) + 
        scale_x_discrete(limits = c("up", "left", "right", "noFix")) +
        scale_color_discrete(limits = c("up", "left", "right", "noFix")) + ggtitle(i)
    print(p)
}


#'# subject-wise fixation frequency
dat %>%
    group_by(Condition, fix, id) %>%
    summarise(n = n()) -> fixation

f1 <- lmer(n ~ fix * Condition + (1 + fix + Condition | id), data = fixation)
summary(f1)
plot(ggpredict(f1, terms = c("fix", "Condition"), type = "fe"))