library(tidyverse)

vAll <- scan("data/all-genes.txt", what=character())
core282 <- scan("sets/core282.txt", what=character())

GS <- read_tsv("data/pub-sets.tsv", col_types=cols()) %>%
    select(-URL) %>%
    mutate(across(`Gene Set`, str_split, ", "),
           across(PMID, partial(str_c, "PMID"))) %>%
    with(set_names(`Gene Set`, PMID))

nDiff <- map(GS, setdiff, vAll) %>% map_int(length)
stopifnot(all(nDiff == 0))

iwalk(GS, ~cat(.x, sep="\n",
               file=str_c("sets/", .y, ".txt")))
