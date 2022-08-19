library(tidyverse)

f <- function(v) str_split(v, ",", simplify = TRUE)[, 1]

s2 <- openxlsx::read.xlsx("PMID27048245-TableS2.xlsx") %>%
    select(Gene = Genes.names, log2FC = `log.2.Res./.Par`) %>%
    mutate(across(Gene, f)) %>%
    tibble()

s3 <- openxlsx::read.xlsx("PMID27048245-TableS3.xlsx") %>%
    select(Gene = Genes.names, log2FC = `ratios.Res./.Par`) %>%
    mutate(across(Gene, f)) %>%
    tibble()

g2 <- filter(s2, abs(log2FC) > 2)$Gene
g3 <- filter(s3, abs(log2FC) > 2)$Gene

intersect(g2, g3)
