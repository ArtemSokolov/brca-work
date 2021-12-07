library(tidyverse)
library(synExtra)

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
    Z <- filter(X, Drug %in% drugs)
    Md <- M %>% filter(Drug %in% drugs) %>%
        nest(Model=c(Slope,Bias,Noise)) %>%
        mutate(Fit = map(Model, partial(CI95, unique(Z$nFeats)))) %>%
        unnest(Fit)

    ggplot(Z, aes(x=nFeats, y=AUC)) + theme_bw() +
        geom_line(data=Md, color='red') +
        geom_errorbar(data=Md, aes(ymin=Lower, ymax=Upper),
                      color='red', alpha=0.4,
                      width=10, linetype='dashed') +
        geom_point() + facet_wrap(~Drug, nrow=1) +
        xlab("Number of Features") +
        theme(strip.background = element_blank())
}

## Fetch all background files
##synapser::synLogin()
##syn <- synDownloader("data/syn")
##fns <- synChildren("syn26487075") %>% syn()
fns <- list.files("data/syn", full.names=TRUE)

## Load everything
X <- tibble(fn = fns) %>%
    mutate(fnb    = basename(fn),
           tokens = str_split(fnb, "-"),
           nFeats = map_int(tokens, ~as.integer(.x[2])),
           Iter   = map_int(tokens, ~as.integer(.x[3])),
           Data   = map(fn, read_csv, col_types=cols())) %>%
    select(-fn, -fnb, -tokens) %>% unnest(Data)

MT <- read_csv("data/agents_metadata_nov30-2021.csv", col_types=cols()) %>%
    mutate(agent = str_split(agent, "/", simplify=TRUE)[,1]) %>%
    select(Drug = agent, Class=drug_class)

## Fit simple linear models to each drug
M <- X %>% group_by(Drug) %>%
    summarize( across(c(nFeats, AUC), list) ) %>%
    mutate(Model = map2(nFeats, AUC, fitLM)) %>%
    select(Drug, Model) %>% unnest(Model) %>%
    left_join(MT, by="Drug")

## Plot an overview
vh <- c("Alpelisib", "Pin1-3", "Ipatasertib", "Everolimus")
pal <- c("black", ggthemes::few_pal()(8), "darkgray")
M1 <- M %>% mutate( Label=ifelse(Drug %in% vh, Drug, "") )
gg1 <- ggplot(M1, aes(x=Slope, y=Bias, color=Class)) +
    theme_bw() + geom_point() +
    scale_color_manual(values=pal, guide=FALSE) +
    ggrepel::geom_text_repel(aes(label=Label), show.legend=FALSE)
gg2 <- ggplot(M1, aes(x=Bias, y=Noise, color=Class)) +
    theme_bw() + geom_point() +
    scale_color_manual(values=pal) +
    ggrepel::geom_text_repel(aes(label=Label), show.legend=FALSE)

## Plot a handful of examples
gg3 <- plotDrug(X, M, vh)

gg <- gridExtra::grid.arrange(gg1, gg2, gg3,
                              layout_matrix=rbind(c(1,2),
                                                  c(3,3)),
                              heights=c(1.7,1), widths=c(0.8,1))

ggsave("plots/bkmodel.png", gg, width=10, height=6)
write_csv(M, "output/BK-models.csv")
