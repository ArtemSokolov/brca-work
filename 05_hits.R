library(tidyverse)

## Harmonic mean to aggregate p values for a single per-signature metric
hmean <- function(x) {
    length(x) / sum(1 / x)
}

split2 <- function(.df, .col) {
    split(.df, .df[[.col]]) |> purrr::map(dplyr::select, -dplyr::all_of(.col))
}

## Background models for RNAseq performance
bk <- read_csv("output/BK-models.csv", col_types = cols()) |>
  select(-Drug) |> rename(Drug = Generic) |> split2("Modality")

## Load metadata
meta <- yaml::read_yaml("literature.yml") %>%
  map(modify_at, "signature", list) %>%
  map(~tibble(!!!.x)) %>%
  bind_rows(.id = "PMID") %>%
  mutate(
    size = map_int(signature, length),
    across(drug, recode, `hormonal therapy` = "hormonal"),
    across(modality, recode,
      `Gene Essentiality` = "Essentiality",
      PCR = "Targeted PCR")
  ) %>%
  select(Signature = PMID, Of = drug, Type = type,
    Modality = modality, Size = size)

## Load AUC values, match up against metadata and background models
## Compute p values associated with each AUC measure
auc <- read_csv("data/aucs/lit-gene.csv", col_types = cols())
x <- inner_join(auc, bk$RNA, by = "Drug") %>%
    inner_join(meta, by = "Signature") %>%
    mutate(
        AUCe = Size * Slope + Intercept,
        pval  = pmap_dbl(list(AUC, AUCe, SD), pnorm, lower.tail = FALSE)
    ) %>%
    select(Type, Signature, Of, Modality,
      Drug, Class, AUC, pval, `F-statistic`) %>%
    filter(Drug %in% Of)

hmp <- group_by(x, Signature, Of, Type) %>%
    summarize(across(pval, hmean), .groups = "drop") %>%
    mutate(
        Drug = "Harmonic\nmean p-val",
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
            str_c(Signature, " (", Of, ")")
        ),
        Highlight = ifelse(Of == Drug, "yes", "no"),
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
gg <- ggplot(y, aes(y = Signature, x = Drug, fill = nlogp)) +
    theme_bw() + geom_tile() +
    geom_tile(
        data = filter(y, Highlight == "yes"),
        color = "black", size = 1
    ) +
    geom_text(aes(label = Label), color = "black", angle = 0) +
    scale_fill_gradientn(
        name = "p-value", colors = pal, limits = c(0, 3.1),
        labels = vs, breaks = -log10(vs)
    ) +
    xlab("Drug") +
    facet_grid(Type~Class, space = "free", scales = "free") +
    theme(axis.text.x  = etxt(12, angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y  = etxt(12), axis.title = etxt(14),
          legend.text = etxt(12, angle = 90, hjust = 1),
          legend.title = etxt(14),
          legend.key.width = unit(3, "line"),
          legend.position = "bottom",
          strip.text.y = etxt(12),
          strip.text.x = etxt(12),
          strip.background = element_blank())

# Save the figure
pfx <- str_c("plots/Figure-Lit-", Sys.Date())
ggsave(str_c(pfx, ".png"), gg, width = 15, height = 12)
ggsave(str_c(pfx, ".pdf"), gg, width = 15, height = 12)
