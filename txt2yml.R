library(tidyverse)

# Identify all mined signatures
vfn <- c(
    list.files("data/sigs/lit", full.names = TRUE),
    list.files("data/sigs/univ", full.names = TRUE) %>%
        keep(~grepl("PMID", .x))
)
vfn <- set_names(vfn, str_sub(basename(vfn), 1, -5))

# Load the signatures
sigs <- vfn %>%
    map(scan, what = character()) %>%
    map(~list(Signature = .x))

# Load all relevant metadata
metad <- read_csv("data/agents_metadata.csv")
metal <- read_tsv("data/lit-meta.tsv") %>%
    select(-Title, -URL)

# Combine signatures with metadata
res <- metal %>%
    mutate(across(PMID, ~str_c("PMID", .x))) %>%
    split(., .$PMID) %>%
    map(select, -PMID) %>%
    map(as.list) %>%
    list_merge(!!!sigs)
