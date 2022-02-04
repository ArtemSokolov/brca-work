library(tidyverse)

## Fit simple linear regression models
fitLM <- function(nf, auc) {
    b1 <- cov(nf,auc) / var(nf)
    b0 <- mean(auc) - b1*mean(nf)
    s  <- sd(auc - (b0 + nf*b1))

    tibble(Slope=b1, Bias=b0, Noise=s)
}

## Identify the 95% CI of AUC for a given number of features
CI95 <- function(nf, mdl) {
    tibble(nFeats=nf) %>%
        mutate(AUC   = mdl$Slope * nFeats + mdl$Bias,
               Upper = AUC + 1.96*mdl$Noise,
               Lower = AUC - 1.96*mdl$Noise)
}

## Plot a model for a subset of drugs
plotDrug <- function(X, M, drugs) {
    Z  <- inner_join(X, M, by="Drug") %>% filter(Generic %in% drugs)
    Md <- M %>% filter(Generic %in% drugs) %>%
        nest(Model=c(Slope,Bias,Noise)) %>%
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

## Fetch all background files
##library(synExtra)
##synapser::synLogin()
##syn <- synDownloader("aucs/rnabg")
##fns <- synChildren("syn26487075") %>% syn()

## RNAseq background
X <- tibble(fn = list.files("aucs/rnabg", full.names=TRUE)) %>%
    mutate(fnb    = basename(fn),
           tokens = str_split(fnb, "-"),
           nFeats = map_int(tokens, ~as.integer(.x[2])),
           Iter   = map_int(tokens, ~as.integer(.x[3])),
           Data   = map(fn, read_csv, col_types=cols())) %>%
    select(-fn, -fnb, -tokens) %>% unnest(Data)

## MassSpec background
XMS <- tibble(fn = list.files("aucs/msbg", full.names=TRUE)) %>%
    mutate(fnb    = basename(fn),
           tokens = str_split(fnb, "-"),
           nFeats = map_int(tokens, ~as.integer(.x[2])),
           Iter   = map_int(tokens, ~as.integer(.x[3])),
           Data   = map(fn, read_csv, col_types=cols())) %>%
    select(-fn, -fnb, -tokens) %>% unnest(Data)

## Quick correlation check
inner_join(select(X, Drug, AUC), select(XMS, Drug, AUC),
           by="Drug", suffix=c("_rna", "_ms")) %>%
    group_by(Drug) %>%
    summarize_at( c("AUC_rna", "AUC_ms"), mean ) %>%
    with( cor(AUC_rna, AUC_ms, method="pearson") )
    
## Other data files
MT <- read_csv("data/agents_metadata_nov30-2021.csv", col_types=cols()) %>%
    mutate(agent = str_split(agent, "/", simplify=TRUE)[,1]) %>%
    select(Drug = agent, Generic = generic_name)

F <- read_csv("data/anova.csv", col_types=cols()) %>%
    select(Generic = generic_name, Class = drug_class, `F-statistic` = F_GR_AOC)

## Fit simple linear models to each drug
M <- X %>% group_by(Drug) %>%
    summarize( across(c(nFeats, AUC), list) ) %>%
    mutate(Model = map2(nFeats, AUC, fitLM)) %>%
    select(Drug, Model) %>% unnest(Model) %>%
    left_join(MT, by="Drug") %>% inner_join(F, by="Generic")

## Plot an overview
vh <- c("alpelisib", "Pin1-3", "ipatasertib", "everolimus")
pal <- c("black", colorblindr::palette_OkabeIto)
M1 <- M %>% mutate( Label=ifelse(Generic %in% vh, Generic, "") )
gg1 <- ggplot(M1, aes(x=Slope, y=Bias, color=Class)) +
    theme_bw() + geom_point() +
    scale_color_manual(values=pal, guide='none') +
    ylab("Intercept") +
    ggrepel::geom_text_repel(aes(label=Label), show.legend=FALSE)
gg2 <- ggplot(M1, aes(x=Bias, y=Noise, color=Class)) +
    theme_bw() + geom_point() +
    scale_color_manual(values=pal) +
    xlab("Intercept") +
    ggrepel::geom_text_repel(aes(label=Label), show.legend=FALSE)

## Plot a handful of examples
gg3 <- plotDrug(X, M, vh)

rsq <- with(M1, cor(`F-statistic`, Bias)^2) %>%
    round(3) %>% str_c("~italic(r)^2 == ", ., "")
             
gg4 <- ggplot(M1, aes(x=`F-statistic`, y=Bias, color=Class)) +
    ylab("Intercept") + xlab("Correlation w/ subtype") +
    theme_bw() + geom_point() + scale_color_manual(values=pal) +
    annotate('text', 1.2, 0.5, hjust=0, vjust=1, label=rsq, parse=TRUE)


gg <- gridExtra::grid.arrange(gg1, gg2, gg3,
                              layout_matrix=rbind(c(1,2),
                                                  c(3,3)),
                              heights=c(1.7,1), widths=c(0.8,1))

ggsave("plots/04-bkmodel.png", gg, width=10, height=6)

ggplot2::ggsave("plots/04-bias-noise.png", gg2, width=6, height=4)
ggplot2::ggsave("plots/04-bias-cor.png", gg4, width=6, height=4)

write_csv(M, "output/BK-models.csv")
