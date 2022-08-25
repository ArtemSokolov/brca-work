library(tidyverse)

f <- function(sheet) {
    openxlsx::read.xlsx("PMID31462705-TableS1.xlsx", sheet) %>%
        tibble() %>%
        select(Gene = 1, pval = 2, log2FC = 3) %>%
        mutate(padj = p.adjust(pval, method = "bonferroni"))
}

skbr3 <- f("SKBR3 THZ1, p<0.05")
bt474 <- f("BT474 THZ1, p<0.05")

v1 <- filter(skbr3, abs(log2FC) > 2, padj < 0.05)$Gene
v2 <- filter(bt474, abs(log2FC) > 2, padj < 0.05)$Gene

intersect(v1, v2)
