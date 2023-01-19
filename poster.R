#+ message = F
library(tidyverse)
library(data.table)
library(lme4)
library(lmerTest)
library(car)
library(sjPlot)
library(ggsignif)
library(cowplot)

theme_set(theme_classic(base_size = 8, base_family = "Helvetica")) 


#'# data loading
d <- fread("exp1_data_behavior.csv", sep = ",", header = T)
d <- mutate(d, Corr = ifelse(ChosenITM == CorrectITM, 1, 0))
d <- mutate(d, Correct = ifelse(ChosenITM == CorrectITM, "correct", "incorrect"))
d$Correct <- as.factor(d$Correct)
d$ChoiceRT <- as.numeric(d$ChoiceRT)
d$nVal3 <- as.numeric(as.character(d$Val3))
d <- na.omit(d) # 67 trials omitted
# d <- subset(d, d$participant != "sub01") # accuracy below mean - 2sd at Val3 = 63

d1 <- rbind(subset(d, d$UVal < 90 & d$ChosenITM == "up"),
            subset(d, d$LVal < 90 & d$ChosenITM == "left"),
            subset(d, d$RVal < 90 & d$ChosenITM == "right")) # subset dud chosen

d2 <- rbind(subset(d, d$UVal < 90 & d$ChosenITM != "up"),
            subset(d, d$LVal < 90 & d$ChosenITM != "left"),
            subset(d, d$RVal < 90 & d$ChosenITM != "right")) # subset dud not chosen


#'# visualization for whole data
d %>% group_by(Val3, participant) %>%
    summarise(Accuracy = mean(Corr)) %>%
    ggplot() + geom_line(aes(x = Val3, y = Accuracy, color = participant)) + xlab("Number of dud stimulus") +
    stat_summary(fun = "mean", geom = "line", aes(x = Val3, y = Accuracy), size = 1) + guides(color = F)

d %>% group_by(Val3, participant, Correct) %>%
    summarise(RT = mean(ChoiceRT)) %>%
    ggplot() + geom_line(aes(x = Val3, y = RT, color = participant))  + xlab("Number of dud stimulus") +
    stat_summary(fun = "mean", geom = "line", aes(x = Val3, y = RT), size = 1) + facet_grid(Correct ~ .) + guides(color = F)

d %>% group_by(Val3, participant, Correct) %>%
    summarise(Confidence = mean(Conf)) %>%
    ggplot() + geom_line(aes(x = Val3, y = Confidence, color = participant))  + xlab("Number of dud stimulus") +
    stat_summary(fun = "mean", geom = "line", aes(x = Val3, y = Confidence), size = 1) + facet_grid(Correct ~ .)


d$Val3 <- as.factor(d$Val3)
f1 <- glmer(Corr ~ Val3 + (1 + Val3|participant), data = d, family = binomial)
summary(f1)
Anova(f1)
# summary(glht(f1, mcp(Val3 = "Tukey")))
ap <- p.adjust(as.data.frame(summary(f1)[10])[2:6, 4], method = "holm") # adjusted p values for correct choices
pv <- data.frame(x = c(0, 0, 0, 0, 0),
                  y = c(0.84, 0.85, 0.86, 0.87, 0.88),
                  xend = c(27, 45, 63, 76, 86),
                  yend = c(0.84, 0.85, 0.86, 0.87, 0.88),
                  values = (as.character(round(ap, 3))))
p1 <- plot_model(f1, type = "pred", terms = "Val3", pred.type = "fe") + 
    xlab("Number of dud stimulus") +  ylab("Accuracy") + ggtitle("") + ylim(0.68, 0.9) +
    geom_signif(data = pv, aes(x = x, xend = xend, y = y, yend = yend, annotation = values, textsize = 2), 
                stat = "identity", tip_length = 1) + guides(color = F)
p1


f2 <- glmer(Corr ~ I(nVal3) + I(nVal3^2) + (1 + I(nVal3) + I(nVal3-2)|participant), data = d, family = binomial)
summary(f2)
Anova(f2)
p2 <- plot_model(f2, type = "pred", terms = "nVal3 [all]", pred.type = "fe") + 
    xlab("Number of dud stimulus") +  ylab("Accuracy") + ggtitle("")
p2


f3 <- lmer(ChoiceRT ~ Val3 * Correct + (1 + Val3 + Correct|participant), data = d)
summary(f3)
Anova(f3)

dl3 <- difflsmeans(f3) # uncorrected linear contrasts 
dl3[17:21,] # p values for correct choices
ap1 <- p.adjust(dl3[17:21,][, 7] , method = "holm") # adjusted p values for correct choices
dl3[68:72,] # p values for incorrect choices
ap2 <- p.adjust(dl3[68:72,][, 7] , method = "holm") # adjusted p values for incorrect choices

pv1 <- data.frame(x = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                 y = c(0.75, 0.7, 0.65, 0.6, 0.55, 1.45, 1.5, 1.55, 1.6, 1.65),
                 xend = c(27, 45, 63, 76, 86, 27, 45, 63, 76, 86),
                 yend = c(0.75, 0.7, 0.65, 0.6, 0.55, 1.45, 1.5, 1.55, 1.6, 1.65),
                 group_col = c(rep("correct", 5), rep("incorrect", 5)),
                 values = (as.character(c(round(ap1, 3), round(ap2, 3)))))

p3 <- plot_model(f3, type = "pred", terms = c("Val3", "Correct"), pred.type = "fe") +
    xlab("Number of dud stimulus") +  ylab("Response time") + ggtitle("") + ylim(0.55, 1.7) + 
    geom_signif(data = pv1, aes(x = x, xend = xend, y = y, yend = yend, annotation = values, color = group_col, textsize = 2), 
                stat = "identity", tip_length = 1) + guides(color = F)
p3


f4 <- lmer(Conf ~ Val3 * Correct + (1 + Val3 + Correct|participant), data = d)
summary(f4)
Anova(f4)
difflsmeans(f4) # uncorrected linear contrasts

dl4 <- difflsmeans(f4) # uncorrected linear contrasts 
dl4[17:21,] # p values for correct choices
ap1 <- p.adjust(dl4[17:21,][, 7] , method = "holm") # adjusted p values for correct choices
dl4[68:72,] # p values for incorrect choices
ap2 <- p.adjust(dl4[68:72,][, 7] , method = "holm") # adjusted p values for incorrect choices

pv2 <- data.frame(x = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                 y = c(3.5, 3.6, 3.7, 3.8, 3.9, 1.9, 1.8, 1.7, 1.6, 1.5),
                 xend = c(27, 45, 63, 76, 86, 27, 45, 63, 76, 86),
                 yend = c(3.5, 3.6, 3.7, 3.8, 3.9, 1.9, 1.8, 1.7, 1.6, 1.5),
                 group_col = c(rep("correct", 5), rep("incorrect", 5)),
                 values = (as.character(c(round(ap1, 3), round(ap2, 3)))))

p4 <- plot_model(f4, type = "pred", terms = c("Val3", "Correct"), pred.type = "fe") +
    xlab("Number of dud stimulus") +  ylab("Confidence") + ggtitle("") + ylim(1.5, 4) + 
    geom_signif(data = pv2, aes(x = x, xend = xend, y = y, yend = yend, annotation = values, color = group_col, textsize = 2),  
                stat = "identity") + guides(color = F)
p4


#'# target vs. distractor comparison (test of IIA)
d1 %>% group_by(Val3, participant) %>%
    summarise(dud_choice = n()) %>%
    ggplot() + geom_line(aes(x = Val3, y = dud_choice, color = participant)) +
    stat_summary(fun = "mean", geom = "line", aes(x = Val3, y = dud_choice), size = 1) +
    ylab("Number of dud choice") + ylim(0, 40)

d2 %>% group_by(Val3, participant) %>%
    summarise(Accuracy = mean(Corr)) %>%
    ggplot() + geom_line(aes(x = Val3, y = Accuracy, color = participant)) +
    stat_summary(fun = "mean", geom = "line", aes(x = Val3, y = Accuracy), size = 1) + ylab("Relative accuracy (dud excluded)")

d2 %>% group_by(Val3, participant, Corr) %>%
    summarise(RT = mean(ChoiceRT)) %>%
    ggplot() + geom_line(aes(x = Val3, y = RT, color = participant)) +
    stat_summary(fun = "mean", geom = "line", aes(x = Val3, y = RT), size = 1) + facet_grid(Corr ~ .)

d2 %>% group_by(Val3, participant, Corr) %>%
    summarise(Confidence = mean(Conf)) %>%
    ggplot() + geom_line(aes(x = Val3, y = Confidence, color = participant)) +
    stat_summary(fun = "mean", geom = "line", aes(x = Val3, y = Confidence), size = 1) + facet_grid(Corr ~ .)


d2$Val3 <- as.factor(d2$Val3)
f5 <- glmer(Corr ~ Val3 + (1 + Val3|participant), data = d2, family = binomial)
summary(f5)
Anova(f5)
ap <- p.adjust(as.data.frame(summary(f5)[10])[2:6, 4], method = "holm") # adjusted p values for correct choices
pv <- data.frame(x = c(0, 0, 0, 0, 0),
                 y = c(0.84, 0.85, 0.86, 0.87, 0.88),
                 xend = c(27, 45, 63, 76, 86),
                 yend = c(0.84, 0.85, 0.86, 0.87, 0.88),
                 values = (as.character(round(ap, 3))))
p5 <- plot_model(f5, type = "pred", terms = "Val3", pred.type = "fe") + 
    xlab("Number of dud stimulus") +  ylab("Relative accuracy (dud excluded)") + ggtitle("") + ylim(0.68, 0.9) +
    geom_signif(data = pv, aes(x = x, xend = xend, y = y, yend = yend, annotation = values, textsize = 2), 
                stat = "identity", tip_length = 1) + guides(color = F)
p5


d2$Val3 <- as.numeric(as.character(d2$Val3))/100
f6 <- glmer(Corr ~ Val3 + (1 + Val3|participant), data = d2, family = binomial(link = probit))
summary(f6)
Anova(f6)


fig1 <- plot_grid(p1, p5, p3, p4, labels = c("a", "b", "c", "d"), align = "h", scale = 0.99, vjust = 0.99)
ggsave(file = "fig1.jpg", plot = fig1, dpi = 500, width = 6, height = 6)


#'# normalization model