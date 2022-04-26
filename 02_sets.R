library(tidyverse)

lit <- read_csv("data/pmids_for_drugs_tidy.csv", col_types = cols()) %>%
    mutate(PMID = str_c("PMID:", PMID))

# Count the number of unique drugs mentioned in each PMID
unq_cnt <- group_by(lit, PMID) %>%
    summarize(nDrug = length(unique(Drug)), .groups = "drop")

# Top 50 papers for each drug
res <- inner_join(lit, unq_cnt, by = "PMID") %>%
    select(Drug, PMID, nDrug, Evidence, Year, nCite) %>%
    group_by(Drug) %>%
    arrange(nDrug, desc(Evidence), desc(Year), desc(nCite)) %>%
    slice_head(n = 50) %>%
    ungroup()

write_csv(res, "output/lit-cand.csv")
