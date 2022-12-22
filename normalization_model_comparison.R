#+ message = F
library(tidyverse)

fits_logit$model <- "logit"
fits_dnm$model <- "dnm"
fits_dnm2$model <- "dnm2"

fittings <- rbind(fits_logit, fits_dnm, fits_dnm2)
fittings %>%
    group_by(model) %>%
    summarise(summed_aic = sum(aic))