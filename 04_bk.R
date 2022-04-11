library(tidyverse)

## Fit simple linear regression models
fitLM <- function(nf, auc) {
    b1 <- cov(nf, auc) / var(nf)
    b0 <- mean(auc) - b1 * mean(nf)
    s  <- sd(auc - (b0 + nf * b1))

    tibble(Slope = b1, Intercept = b0, SD = s)
}

## Identify the 95% CI of AUC for a given number of features
CI95 <- function(nf, mdl) {
    tibble(nFeats = nf) %>%
        mutate(AUC   = mdl$Slope * nFeats + mdl$Intercept,
               Upper = AUC + 1.96 * mdl$SD,
               Lower = AUC - 1.96 * mdl$SD)
}

## Reads a background file and isolates #Features and Iteration index
read_bg <- function(fn)
    read_csv(fn, col_types = cols()) %>%
        mutate(tokens = str_split(Signature, "-"),
               nFeats = map_int(tokens, ~as.integer(.x[2])),
               Iter   = map_int(tokens, ~as.integer(.x[3])),
               tokens = NULL)

## Load background files
XRNA <- read_bg("data/aucs/rnabg-aucs.csv")
XMS  <- read_bg("data/aucs/msbg-aucs.csv")
XX   <- bind_rows(RNA = XRNA, MS = XMS, .id = "Modality")

## Other data files
MT <- read_csv("data/agents_metadata.csv", col_types = cols()) %>%
    mutate(agent = str_split(agent, "/", simplify = TRUE)[, 1]) %>%
    select(Drug = agent, Generic = generic_name)

F <- read_csv("data/anova.csv", col_types = cols()) %>%
    select(Generic = generic_name, Class = drug_class, `F-statistic` = F_GR_AOC)

## Fit simple linear models to each drug/platform pair
M <- group_by(XX, Modality, Drug) %>%
    summarize(across(c(nFeats, AUC), list), .groups = "drop") %>%
    mutate(Model = map2(nFeats, AUC, fitLM)) %>%
    select(Modality, Drug, Model) %>%
    unnest(Model) %>%
    left_join(MT, by = "Drug") %>%
    inner_join(F, by = "Generic")
write_csv(M, "output/BK-models.csv")

## Okabe Ito palette and a void plot
pal <- c("black", "#E69F00", "#56B4E9",
    "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7", "#999999")
void <- ggplot() + theme_void()

## Plot a model for a subset of drugs
plotDrugs <- function(X, M, drugs) {
    Z  <- dplyr::inner_join(X, M, by = "Drug") %>%
        dplyr::filter(Generic %in% drugs)
    Md <- M %>%
        filter(Generic %in% drugs) %>%
        nest(Model = c(Slope, Intercept, SD)) %>%
        mutate(Fit = map(Model, partial(CI95, unique(Z$nFeats)))) %>%
        unnest(Fit)

    ggplot(Z, aes(x = nFeats, y = AUC)) + theme_bw() +
        geom_line(data = Md, color = "red") +
        geom_errorbar(data = Md, aes(ymin = Lower, ymax = Upper),
                      color = "red", alpha = 0.4,
                      width = 10, linetype = "dashed") +
        geom_point() + facet_wrap(~Generic, nrow = 1) +
        xlab("Number of Features") +
        theme(strip.background = element_blank())
}

## Plot property p1 against p2 for all models
## Explicitly highlights drugs in h
plotModels <- function(M, p1, p2, h) {
    M1 <- M %>% dplyr::mutate(Label = ifelse(Generic %in% h, Generic, ""))
    ggplot(M1, aes(x = {{p1}}, y = {{p2}}, color = Class)) +
        theme_bw() + geom_point() + scale_color_manual(values = pal) +
        ggrepel::geom_text_repel(aes(label = Label), show.legend = FALSE)
}

## Main figure -- results in the RNA space
vh <- c("alpelisib", "Pin1-3", "ipatasertib", "everolimus")
MRNA <- M %>% filter(Modality == "RNA")
rsq <- with(MRNA, cor(`F-statistic`, Intercept)^2) %>%
    round(3) %>% str_c("~italic(r)^2 == ", ., "")

## Individual panels
pA <- plotDrugs(XRNA, MRNA, vh)
pB <- plotModels(MRNA, Slope, Intercept, vh) + guides(color = "none")
pC <- plotModels(MRNA, `F-statistic`, Intercept, vh) +
    annotate("text", 1.2, 0.5, hjust = 0, vjust = 1, label = rsq, parse = TRUE)

## Assemble the main figure
f1 <- cowplot::plot_grid(pB, void, pC, ncol = 3, labels = c("b", "", "c"),
                         rel_widths = c(0.75, 0.02, 1), label_size = 20)
ff <- cowplot::plot_grid(pA, NULL, f1, ncol = 1, labels = c("a", "", ""),
                         rel_heights = c(1, 0.02, 1.7), label_size = 20)
pfx <- str_c("plots/Figure-BK-", Sys.Date())
str_c(pfx, ".png") %>% ggsave(ff, width = 10, height = 6)
str_c(pfx, ".pdf") %>% ggsave(ff, width = 10, height = 6)

## Supplement -- intercept vs. standard deviation
ggBKSD <- plotModels(MRNA, Intercept, SD, vh) + ylab("Standard Deviation")
ggplot2::ggsave("plots/Suppl-BKSD.png", ggBKSD, width = 6, height = 4)

## Supplement -- RNA vs Mass Spec space
Mvs <- select(M, Modality, Drug, Generic, Intercept, Class) %>%
    pivot_wider(names_from = "Modality", values_from = "Intercept")
rsq_vs <- with(Mvs, cor(RNA, MS)^2) %>% round(3) %>% 
    str_c("~italic(r)^2 == ", ., "")

ggvs <- plotModels(Mvs, RNA, MS, c()) +
    xlab("Intercept (RNAseq)") + ylab("Intercept (Mass Spec)") +
    annotate("text", 0.7, 0.5, hjust = 0, vjust = 1,
        label = rsq_vs, parse = TRUE)
ggplot2::ggsave("plots/Suppl-BK-RNA-v-MS.png", ggvs, width = 6, height = 4)
