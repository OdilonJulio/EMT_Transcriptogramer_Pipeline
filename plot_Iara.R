library(transcriptogramer)
library(tidyverse)

load("odilon_data/SomenteTranscriptogramas.RData")

map_dfr(ls()[grep("^t", ls())], function(x) {
  
  #normalizar
  transcriptogram <- apply(get(x, envir = .GlobalEnv)@transcriptogramS2[, -c(1,2)], 2, function(x) {x / sum(x)})
  
  data.frame(
    Protein = get(x, envir = .GlobalEnv)@transcriptogramS2[,1],
    Position = get(x, envir = .GlobalEnv)@transcriptogramS2[,2],
    mean = apply(transcriptogram, 1, mean, na.rm = T),
    sd = apply(transcriptogram, 1, sd,  na.rm = T),
    treatment = x
  ) -> res
  return(res)
  
}) -> df

ggplot(df, aes(x = Position, y = mean, col = treatment, group = treatment)) +
  geom_line() +
  facet_wrap(~treatment)