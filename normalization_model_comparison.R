#+ message = F
library(tidyverse)
library(cowplot)

#'# model fitting
#+ message = F
#' logit: multinomial logit without normalization  
#' dnm:   nomalization model under 0 < w < 10  
#' dnm2:   nomalization model under -5 < w < 5  

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
                get_legend(logit),
                labels = c("a", "b", "c", NA), align = "vh", scale = 0.99, vjust = 0.87)
g1


#'# model comparison
#+ message = F
fittings <- as.data.frame(rbind(fits_logit, fits_dnm, fits_dnm2))
fittings %>%
    group_by(model) %>%
    summarise(summed_aic = sum(aic)) %>%
    ggplot() + geom_point(aes(x = model, y = summed_aic))


#'# parameter distribution
#+ message = F
# dnm 
ggplot(fits_dnm) + geom_histogram(aes(x = w)) + xlim(-0.1, 0.1)
t.test(fits_dnm$w, mu = 0)

# dnm2
ggplot(fits_dnm2) + geom_histogram(aes(x = w)) + xlim(-0.1, 0.1)
t.test(fits_dnm2$w, mu = 0)