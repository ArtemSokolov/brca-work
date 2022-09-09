library(tidyverse)

sigs <- yaml::read_yaml("literature_mapped.yml") %>%
  map(`[`, c("signature_gene_mappings", "signature_protein_mappings")) %>%
  map(modify, enframe) %>%
  map(modify, unnest, value) %>%
  map(modify, pull, value)

rnasig <- map(sigs, `[[`, "signature_gene_mappings")
mssig <- map(sigs, `[[`, "signature_protein_mappings") %>%
  keep(~length(.x) > 5)

dir.create("output/lit-gene")
dir.create("output/lit-protein")

iwalk(rnasig, ~cat(.x, sep="\n", file=str_c("output/lit-gene/", .y, ".txt")))
iwalk(mssig, ~cat(.x, sep="\n", file=str_c("output/lit-protein/", .y, ".txt")))
