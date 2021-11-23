library(tidyverse)

GS <- read_tsv("data/pub-sets.tsv", col_types=cols()) %>%
    select(-URL) %>%
    mutate(across(`Gene Set`, str_split, ", "),
           across(PMID, partial(str_c, "PMID:"))) %>%
    with(set_names(`Gene Set`, PMID))

