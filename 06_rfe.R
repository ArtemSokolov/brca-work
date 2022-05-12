library(tidyverse)

synapser::synLogin(Sys.getenv("SYN_USER"), Sys.getenv("SYN_TOKEN"))

# Idenify and download relevant files
v <- synExtra::synGlob("syn30291818", "*", "*_rfe.csv")
fns <- synExtra::synDownloader("data/rfe")(v)

# Loads a file and retrieves a signature
get_sig <- function(fn) {
    read_csv(fn, col_types = cols()) %>%
        filter(is.na(drop_iteration)) %>%
        pull(features)
}

# Load and process all files
res <- tibble(fn = fns) %>%
    mutate(
        drug = str_replace(basename(fn), "_rfe.csv", ""),
        path = map(fn, get_sig)) %>%
        with(set_names(path, drug))

jsonlite::toJSON(res, pretty = TRUE) %>%
    cat(file = "output/BRCA_rfesigs.json")
    