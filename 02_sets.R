library(tidyverse)

vAll <- scan("data/all-genes.txt", what=character())

GS <- list.files("input", pattern="*.txt", full.names=TRUE) %>%
    set_names() %>% map(scan, what=character())

map(GS, ~.x[duplicated(.x)])
map(GS, setdiff, vAll)

