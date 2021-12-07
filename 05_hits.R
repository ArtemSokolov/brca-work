library(tidyverse)
library(synExtra)

## Fetch relevant files
## synapser::synLogin()
## synSet <- synDownloader("data/sets")
## synAUC <- synDownloader("data/aucs")
## fnsSet <- synChildren("syn26529258") %>% synSet() %>% c(synSet("syn26486678"))
## fnsAUC <- synChildren("syn26529259") %>% synAUC() %>% c(synAUC("syn26529261"))
fnsSet <- list.files("data/sets", full.names=TRUE)
fnsAUC <- list.files("data/aucs", full.names=TRUE)

## Load everything and match up results against input sets and background models
MBK <- read_csv("output/BK-models.csv", col_types=cols())
S <- tibble(fn = fnsSet) %>%
    mutate(Name   = str_split(basename(fn), "\\.", simplify=TRUE)[,1],
           Set    = map(fn, scan, what=character(), quiet=TRUE),
           nFeats = map_int(Set, length)) %>% select(-fn)
X <- tibble(fn = fnsAUC) %>%
    mutate(Name = str_split(basename(fn), "-", simplify=TRUE)[,1],
           AUC  = map(fn, read_csv, col_types=cols())) %>%
    inner_join(S, by="Name") %>% select(-fn, -Set) %>%
    unnest(AUC) %>% inner_join(MBK, by="Drug")

## Compute the p values associated with each AUC measure
P <- X %>% mutate(AUCe  = nFeats*Slope + Bias,
                  Slope = NULL, Bias = NULL,
                  pval  = pmap_dbl(list(AUC, AUCe, Noise), pnorm, lower.tail=FALSE))
write_csv(P, "output/summary.csv" )
P <- P %>% mutate(nlog10p = -log10(pval))

P0 <- P %>% mutate(Label = ifelse(pval < 0.01, Drug, ""))
vs <- c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005)
pal <- c("black", ggthemes::few_pal()(8))
ggplot(P0, aes(x=AUC, y=nlog10p, color=Name)) + theme_bw() +
    geom_point() +
    scale_y_continuous(breaks=-log10(vs), labels=as.character(vs), name="p value") +
    scale_color_manual(values=pal) +
    ggrepel::geom_text_repel(aes(label=Label), show.legend=FALSE) +
    ggsave("plots/05-summary.png", width=8, height=6)

P1 <- P %>% filter(Name == "core282")
ggplot(P1, aes(x=Drug, y=nlog10p)) + theme_bw() +
    geom_point() + facet_wrap(~Name) +
    scale_y_continuous(breaks=-log10(vs), labels=as.character(vs), name="p value") +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
          strip.background = element_blank()) +
    ggsave("plots/05-core282.png", width=10, height=6)

hmean <- function(x) {
    length(x) / sum( 1/x )
}

pal <- c("#F7F7F7", rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))[4:7])
##pal <- RColorBrewer::brewer.pal( 9, "YlOrBr" ) %>% colorRampPalette
HMP <- P %>% group_by( Name, Class ) %>%
    summarize(hmp = hmean(pval), .groups="drop") %>%
    mutate(nlogp = -log10(hmp),
           Label = ifelse(hmp < 0.05, as.character(round(hmp, 3)), ""))

vs <- c(0.5, 0.1, 0.05, 0.01)
ggplot( HMP, aes(x=Name, y=Class, fill=nlogp) ) +
    theme_bw() + geom_tile() +
    geom_text(aes(label=Label)) +
    scale_fill_gradientn(colors=pal, breaks=-log10(vs),
                         labels=vs, name="HMP") +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
    ggsave( "plots/05-hmp.png", width=9, height=7 )

