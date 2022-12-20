#+ message = F
library(tidyverse)
library(data.table)
library(lme4)
library(lmerTest)
library(car)
library(multcomp)
library(sjPlot)

#'# data loading
d <- fread("exp1_data_behavior.txt", sep = ",", header = T)
d <- mutate(d, Corr = ifelse(ChosenITM == CorrectITM, 1, 0))
d$ChoiceRT <- as.numeric(d$ChoiceRT)
d <- na.omit(d) # 67 trials omitted
d <- subset(d, d$participant != "sub01") # accuracy below mean - 2sd at Val3 = 63


#'# visualization
d %>% group_by(Val3, participant) %>%
    summarise(Accuracy = mean(Corr)) %>%
    ggplot() + geom_line(aes(x = Val3, y = Accuracy, color = participant)) +
    stat_summary(fun = "mean", geom = "line", aes(x = Val3, y = Accuracy), size = 1)

d %>% group_by(Val3, participant, Corr) %>%
    summarise(RT = mean(ChoiceRT)) %>%
    ggplot() + geom_line(aes(x = Val3, y = RT, color = participant)) + facet_grid(Corr ~ .) +
    stat_summary(fun = "mean", geom = "line", aes(x = Val3, y = RT), size = 1) + facet_grid(Corr ~ .)

d %>% group_by(Val3, participant, Corr) %>%
    summarise(Confidence = mean(Conf)) %>%
    ggplot() + geom_line(aes(x = Val3, y = Confidence, color = participant)) + facet_grid(Corr ~ .) +
    stat_summary(fun = "mean", geom = "line", aes(x = Val3, y = Confidence), size = 1) + facet_grid(Corr ~ .)


#'# statistical test
d$Val3 <- as.factor(d$Val3)
f1 <- glmer(Corr ~ Val3 + (1 + Val3|participant), data = d, family = binomial)
summary(f1)
Anova(f1)
# summary(glht(f1, mcp(Val3 = "Tukey")))
plot_model(f1, type = "pred", terms = "Val3", pred.type = "fe") 

d$Corr <- as.factor(d$Corr)
f2 <- lmer(ChoiceRT ~ Val3 * Corr + (1 + Val3 + Corr|participant), data = d)
summary(f2)
Anova(f2)
difflsmeans(f2) # uncorrected linear contrasts 
# summary(glht(f2, mcp(Val3 = "Tukey")))
plot_model(f2, type = "pred", terms = c("Val3", "Corr"), pred.type = "fe") 
plot_model(f2, type = "diag", terms = c("Val3", "Corr"), pred.type = "fe")

f3 <- lmer(Conf ~ Val3 * Corr + (1 + Val3 + Corr|participant), data = d)
summary(f3)
Anova(f3)
difflsmeans(f3) # uncorrected linear contrasts
# summary(glht(f3, mcp(Val3 = "Tukey")))
plot_model(f3, type = "pred", terms = c("Val3", "Corr"), pred.type = "fe") 
plot_model(f3, type = "diag", terms = c("Val3", "Corr"), pred.type = "fe")