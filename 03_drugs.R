library(tidyverse)

read_csv('data/grmetrics.csv', col_types=cols()) %>%
    pull(agent) %>% unique() %>%
    map_chr(~str_c('"', .x, '"')) %>%
    str_flatten(',\n  ') %>%
    str_c('params.drugs = [\n  ', ., '\n]\n') %>%
    cat(file="output/drugs.config")

