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
    map(~list(signature = .x))

# Load all relevant metadata
meta <- read_tsv("data/lit-meta.tsv") %>%
    select(-Title, -URL, -Size) %>%
    rename_with(str_to_lower, -PMID) %>%
    rename_with(~str_replace_all(.x, " ", "-")) %>%
    mutate(across(PMID, ~str_c("PMID", .x))) %>%
    mutate(across(type, recode, Universal = "Prior-based"))

# Combine signatures with metadata
res <- meta %>%
    split(., .$PMID) %>%
    map(select, -PMID) %>%
    map(as.list) %>%
    list_merge(!!!sigs)

yaml::write_yaml(res, "literature.yml")
