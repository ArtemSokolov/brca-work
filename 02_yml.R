library(tidyverse)

`%nin%` <- function(x, y) {
  union(setdiff(x, y), setdiff(y, x))
}

# Load the original signatures
df0 <- yaml::read_yaml("literature.yml") %>%
  map(modify_at, "signature", list) %>%
  map(~tibble(!!!.x)) %>%
  bind_rows(.id = "PMID")

# Load the updated mappings
df1 <- yaml::read_yaml("literature_mapped.yml") %>%
  map(modify_at, vars(ends_with("mappings")), enframe) %>%
  map(modify_at, vars(ends_with("mappings")), unnest, value) %>%
  map(modify_at, vars(starts_with("signature")), list) %>%
  map(~tibble(!!!.x)) %>%
  bind_rows(.id = "PMID")

# Verify that old fields are unmodified
dff <- full_join(df0, df1, by = "PMID")
names(df0) %>%
  setdiff("PMID") %>%
  map_lgl(~identical(dff[[str_c(., ".x")]], dff[[str_c(., ".y")]]))

# Verify that the mappings span the same gene domain
map2_lgl(df1$signature_gene_mappings, df1$signature, ~all(.x$name %in% .y))
map2_lgl(df1$signature_protein_mappings, df1$signature, ~all(.x$name %in% .y))

# Verfiy that the mappings match the target gene/protein space
rnaseq <- read_csv("data/rnaseq_log2rpkm.csv")
map_lgl(df1$signature_gene_mappings, ~all(.x$value %in% rnaseq$gene_name))

ms <- read_csv("data/mass_spec.csv")
map_lgl(df1$signature_protein_mappings, ~all(.x$value %in% ms$Gene_Symbol))
