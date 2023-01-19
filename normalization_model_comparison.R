#+ message = F
library(tidyverse)
library(cowplot)

theme_set(theme_classic(base_size = 8, base_family = "Helvetica")) 


#'# model fitting
#+ message = F
#' logit: multinomial logit without normalization  
#' dnm:   nomalization model under 0 < w < 10  
#' dnm2:  nomalization model under -5 < w < 5  

source("fs_logit.R")
source("fs_dnm.R")
source("fs_dnm2.R")

fits_logit$model <- "logit"
fits_dnm$model <- "dnm"
fits_dnm2$model <- "dnm2"


#'# fitting curves
#+ message = F
g1 <- plot_grid(logit + theme(legend.position = "none"),
                dnm   + theme(legend.position = "none"),
                dnm2  + theme(legend.position = "none"), 
                get_legend(logit + guides(color = guide_legend(ncol = 2))),
                labels = c("a", "b", "c", NA), align = "vh", scale = 0.99, vjust = 0.87)
g1


#'# model comparison
#+ message = F
fittings <- as.data.frame(rbind(fits_logit, fits_dnm, fits_dnm2))
fittings %>%
    group_by(model) %>%
    summarise(summed_aic = sum(aic)) %>%
    ggplot() + geom_point(aes(x = model, y = summed_aic)) + 
    scale_x_discrete(limits = c("logit", "dnm", "dnm2")) + ylab("Summed AIC") + ylim(18100, 18150) -> g2
g2

#'# parameter distribution
#+ message = F
# dnm2
t.test(fits_dnm2$w, mu = 0)

counts <- rbind(fits_dnm, fits_dnm2)
ggplot(counts) + geom_histogram(aes(x = w, fill = model), alpha = 0.4, position = "identity") + labs(fill = "") +
    xlim(-0.1, 0.1) + ylim(0, 10) + theme(legend.position = c(0.85, 0.8), legend.key.size = unit(0.2, 'cm')) + ylab("Count") -> g3
g3

fig2 <- plot_grid(logit + theme(legend.position = "none"),
          dnm   + theme(legend.position = "none"),
          dnm2  + theme(legend.position = "none"), 
          g2, g3, get_legend(logit + guides(color = guide_legend(ncol = 2, title = NULL))),
          labels = c("a", "b", "c", "d", "e", NA), align = "vh", scale = 0.99, vjust = 0.87)

ggsave(file = "fig2.jpg", plot = fig2, dpi = 500, width = 4.5, height = 3)