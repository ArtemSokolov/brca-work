library(tidyverse)

## Metadata for additional universal signatures
meta_univ <- tibble(
    Signature = c("core250", "pam50"),
    Size      = c(250, 50),
    Type      = "Universal"
)

## Combine with the metadata for publishes signatures
meta <- read_tsv("data/lit-meta.tsv",
        col_types = cols(PMID = col_character())
    ) %>%
    mutate(
        Signature = str_c("PMID", PMID),
        across(Drug, recode,
            `fulvestrant+ribociclib` = "ribociclib",
            `hormonal therapy` = "hormonal"),
        Type = ifelse(Type != "Universal", "Experimental", "Universal"),
        Of = ifelse(Type == "Universal", Modality, Drug),
        across(Of, recode, `Gene Essentiality` = "Essentiality")
    ) %>%
    select(Signature, Of, Type, Size) %>%
    bind_rows(meta_univ)

## Background models for RNAseq performance
bk <- read_csv("output/BK-models.csv", col_types = cols()) %>%
    filter(Modality == "RNA") %>%
    select(-Modality)

## Exclude some of the experimental drugs to improve interpretability
dexcl <- c(
    "BSJ-01-175",
    "BSJ-03-124",
    "E17",
    "FMF-04-107-2",
    "FMF-04-112-1",
    "THZ-P1-2",
    "THZ-P1-2R",
    "ZZ1-33B"
)

## Load AUC values, match up against metadata and background models
## Compute p values associated with each AUC measure
x <- read_csv("data/aucs/rnasig-aucs.csv", col_types = cols()) %>%
    inner_join(bk, by = "Drug") %>%
    inner_join(meta, by = "Signature") %>%
    mutate(
        AUCe = Size * Slope + Intercept,
        pval  = pmap_dbl(list(AUC, AUCe, SD), pnorm, lower.tail = FALSE)
    ) %>%
    select(Type, Signature, Of, Generic, Class, AUC, pval, `F-statistic`) %>%
    filter(!(Generic %in% dexcl))

## Harmonic mean to aggregate p values for a single per-signature metric
hmean <- function(x) {
    length(x) / sum(1 / x)
}

hmp <- group_by(x, Signature, Of, Type) %>%
    summarize(across(pval, hmean), .groups = "drop") %>%
    mutate(
        Generic = "Harmonic\nmean p-val",
        Class = ""
    )

## Supplementary computations for the figure
y <- bind_rows(x, hmp, .id = "Category") %>%
    mutate(
        nlogp = -log10(pval),
        Label = ifelse(
            (pval < 0.05) | (Category == 2),
            str_sub(as.character(round(pval, 3)), 2),
            ""
        ),
        Signature = ifelse(
            is.na(Of),
            Signature,
            str_c(Signature, "\n(", Of, ")")
        ),
        Highlight = ifelse(Of == Generic, "yes", "no"),
        Class = recode(
            Class,
            `BCL2 family` = "BCL2",
            `DNA damage` = "DNA dmg"
            ),
        Class = factor(Class, c(sort(unique(Class))[-1], ""))
    )

## Short-hand for bold element_text of desired size
etxt <- function(s, ...) {
    ggplot2::element_text(size = s, face = "bold", ...)
}

## Plot all-by-all hits
pal <- c("#F7F7F7", rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu"))[4:7])
vs <- c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005)
gg <- ggplot(y, aes(y = Signature, x = Generic, fill = nlogp)) +
    theme_bw() + geom_tile() +
    geom_tile(
        data = filter(y, Highlight == "yes"),
        color = "black", size = 1
    ) +
    geom_text(aes(label = Label), color = "black", angle = 90) +
    scale_fill_gradientn(
        name = "p-value", colors = pal, limits = c(0, 3.1),
        labels = vs, breaks = -log10(vs)
    ) +
    xlab("Drug") +
    facet_grid(Type~Class, space = "free", scales = "free") +
    theme(axis.text.x  = etxt(12, angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y  = etxt(12), axis.title = etxt(14),
          legend.text = etxt(12), legend.title = etxt(14),
          legend.key.height = unit(3, "line"),
          strip.text.y = element_blank(),
          strip.text.x = etxt(12),
          strip.background = element_blank())
ggsave("plots/05-hits.png", gg, width = 17, height = 9)
