library(tidyverse)

## Fit simple linear regression models
fitLM <- function(nf, auc) {
    b1 <- cov(nf,auc) / var(nf)
    b0 <- mean(auc) - b1*mean(nf)
    s  <- sd(auc - (b0 + nf*b1))

    tibble(Slope=b1, Intercept=b0, SD=s)
}

## Identify the 95% CI of AUC for a given number of features
CI95 <- function(nf, mdl) {
    tibble(nFeats=nf) %>%
        mutate(AUC   = mdl$Slope * nFeats + mdl$Intercept,
               Upper = AUC + 1.96*mdl$SD,
               Lower = AUC - 1.96*mdl$SD)
}

## Plot a model for a subset of drugs
plotDrug <- function(X, M, drugs) {
    Z  <- inner_join(X, M, by="Drug") %>% filter(Generic %in% drugs)
    Md <- M %>% filter(Generic %in% drugs) %>%
        nest(Model=c(Slope,Intercept,SD)) %>%
        mutate(Fit = map(Model, partial(CI95, unique(Z$nFeats)))) %>%
        unnest(Fit)

    ggplot(Z, aes(x=nFeats, y=AUC)) + theme_bw() +
        geom_line(data=Md, color='red') +
        geom_errorbar(data=Md, aes(ymin=Lower, ymax=Upper),
                      color='red', alpha=0.4,
                      width=10, linetype='dashed') +
        geom_point() + facet_wrap(~Generic, nrow=1) +
        xlab("Number of Features") +
        theme(strip.background = element_blank())
}

## Reads a background file and isolates #Features and Iteration index
read_bg <- function(fn)
    read_csv(fn, col_types=cols()) %>%
        mutate(tokens = str_split(Signature, "-"),
               nFeats = map_int(tokens, ~as.integer(.x[2])),
               Iter   = map_int(tokens, ~as.integer(.x[3])),
               tokens = NULL)

## Fetch all background files
## synapser::synLogin()
## synExtra::synDownloader("aucs")("syn27221488", "syn27221487")

## Load background files
XRNA <- read_bg("aucs/rnabg-aucs.csv")
XMS  <- read_bg("aucs/msbg-aucs.csv")

## Quick correlation check
bind_rows(RNA=XRNA, MS=XMS, .id="Modality") %>%
    group_by(Modality, Drug) %>%
    summarize(across(AUC, mean), .groups="drop") %>%
    spread(Modality, AUC) %>%
    with(cor(MS, RNA, method="pearson"))

## Other data files
MT <- read_csv("data/agents_metadata_nov30-2021.csv", col_types=cols()) %>%
    mutate(agent = str_split(agent, "/", simplify=TRUE)[,1]) %>%
    select(Drug = agent, Generic = generic_name)

F <- read_csv("data/anova.csv", col_types=cols()) %>%
    select(Generic = generic_name, Class = drug_class, `F-statistic` = F_GR_AOC)

## Fit simple linear models to each drug
M <- XRNA %>% group_by(Drug) %>%
    summarize( across(c(nFeats, AUC), list) ) %>%
    mutate(Model = map2(nFeats, AUC, fitLM)) %>%
    select(Drug, Model) %>% unnest(Model) %>%
    left_join(MT, by="Drug") %>% inner_join(F, by="Generic")

## Plot an overview
vh <- c("alpelisib", "Pin1-3", "ipatasertib", "everolimus")
pal <- c("black", colorblindr::palette_OkabeIto)
M1 <- M %>% mutate( Label=ifelse(Generic %in% vh, Generic, "") )
gg1 <- ggplot(M1, aes(x=Slope, y=Intercept, color=Class)) +
    theme_bw() + geom_point() +
    scale_color_manual(values=pal, guide='none') +
    ggrepel::geom_text_repel(aes(label=Label), show.legend=FALSE)
gg2 <- ggplot(M1, aes(x=Intercept, y=SD, color=Class)) +
    theme_bw() + geom_point() +
    scale_color_manual(values=pal) +
    ylab("Standard Deviation") +
    ggrepel::geom_text_repel(aes(label=Label), show.legend=FALSE)

## Plot a handful of examples
gg3 <- plotDrug(XRNA, M, vh)

rsq <- with(M1, cor(`F-statistic`, Intercept)^2) %>%
    round(3) %>% str_c("~italic(r)^2 == ", ., "")
             
gg4 <- ggplot(M1, aes(x=`F-statistic`, y=Intercept, color=Class)) +
    ylab("Intercept") + xlab("Correlation w/ subtype") +
    theme_bw() + geom_point() + scale_color_manual(values=pal) +
    ggrepel::geom_text_repel(aes(label=Label), show.legend=FALSE) +
    annotate('text', 1.2, 0.5, hjust=0, vjust=1, label=rsq, parse=TRUE)

## Put everything together
void <- ggplot() + theme_void()
f1 <- cowplot::plot_grid( gg1, void, gg4, ncol=3, labels=c("b", "", "c"),
                          rel_widths=c(0.75,0.02,1), label_size=20 )

ff <- cowplot::plot_grid( gg3, NULL, f1, ncol=1, labels=c("a","",""),
                         rel_heights=c(1,0.02,1.7), label_size=20 )
##    cowplot::draw_plot_label( c("a", "b"), c(0, 0.4) )

str_c("plots/Figure-BK-", Sys.Date(), ".png") %>% ggsave(ff, width=10, height=6)
str_c("plots/Figure-BK-", Sys.Date(), ".pdf") %>% ggsave(ff, width=10, height=6)

##ggplot2::ggsave("plots/04-bias-noise.png", gg2, width=6, height=4)
##ggplot2::ggsave("plots/04-bias-cor.png", gg4, width=6, height=4)

write_csv(M, "output/BK-models.csv")
