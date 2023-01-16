library(lme4)
library(lmerTest)
library(car)
library(ggeffects)

data <- subset(dat3, dat3$choice != "dud")
data$Condition <- as.factor(data$Condition)


#'# number of fixations
f1 <- lm(Conf ~ n_fix * choice * Condition * factor(id), data = data)
summary(f1)
Anova(f1)
plot(ggpredict(f1, terms = c("n_fix", "choice", "Condition"))) 


#'# number of target fixations
f2 <- lm(Conf ~ n_target_fix * choice * Condition * factor(id), data = data)
summary(f2)
Anova(f2)
plot(ggpredict(f2, terms = c("n_target_fix", "choice", "Condition"))) 


#'# number of target fixations
f3 <- lm(Conf ~ p_target_fix * choice * Condition * factor(id), data = data)
summary(f3)
Anova(f3)
plot(ggpredict(f3, terms = c("p_target_fix", "choice", "Condition"))) 