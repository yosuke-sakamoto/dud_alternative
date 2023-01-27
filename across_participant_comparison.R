#'# Stimulus-based fixation frequency (choice considered)
# どのような参加者がdud effectを示しやすいか
# 参加者間の比較・個人差
# 確信度と判断の一次元化
# sub03 almost always gives conf of 4

dat %>% group_by(fixItem, condition, subj, chosenItem) %>%
    summarize(m_pupil = mean(pupil)) %>%
    distinct() -> fd4

ggplot(fd4) + geom_violin(aes(x = chosenItem, y = m_pupil)) + geom_point(aes(x = chosenItem, y = m_pupil)) +
    ylab("Mean fixations per trial") + facet_wrap(. ~ condition)

ggplot(fd4) + geom_violin(aes(x = chosenItem, y = m_pupil, color = fixItem)) + 
    geom_point(aes(x = chosenItem, y = m_pupil, color = fixItem), position = position_dodge(width = 0.85)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylim(2, 6) + ylab("Mean pupil size") + facet_wrap(. ~ condition)
