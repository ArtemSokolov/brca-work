library(tidyverse)

## Fetch relevant files
## library(synExtra)
## synapser::synLogin()
## synLit <- synDownloader("aucs/lit")
## fnsCore <- synChildren("syn26529408") %>% synDownloader("aucs/core")()
## fnsLit  <- synChildren("syn26529259") %>% synDownloader("aucs/lit")()
fnsCore <- list.files("aucs/core", full.names=TRUE)
fnsLit  <- list.files("aucs/lit", full.names=TRUE)
fnsSet  <- list.files("sets", full.names=TRUE)
stopifnot(length(fnsCore) + length(fnsLit) == length(fnsSet))

## Metadata for literature signatures
LM <- read_tsv("data/lit-meta.tsv", col_types=cols(PMID=col_character())) %>%
    select( Name=PMID, SigOf=Drug ) %>%
    mutate( across(SigOf, recode, `fulvestrant+ribociclib` = "fulv+ribo") )

## Load everything and match up results against input sets and background models
MBK <- read_csv("output/BK-models.csv", col_types=cols())
S <- tibble(fn = fnsSet) %>%
    mutate(Name   = str_split(basename(fn), "\\.", simplify=TRUE)[,1],
           Set    = map(fn, scan, what=character(), quiet=TRUE),
           nFeats = map_int(Set, length)) %>% select(-fn)
X <- list( Core = fnsCore, Literature = fnsLit ) %>%
    enframe("Type", "fn") %>% unnest(fn) %>%
    mutate(Name = str_split(basename(fn), "-", simplify=TRUE)[,1],
           AUC  = map(fn, read_csv, col_types=cols())) %>%
    inner_join(S, by="Name") %>% select(-fn, -Set) %>%
    unnest(AUC) %>% inner_join(MBK, by="Drug") %>%
    left_join(LM, by="Name")

## Compute the p values associated with each AUC measure
P <- X %>% mutate(AUCe  = nFeats*Slope + Bias,
                  pval  = pmap_dbl(list(AUC, AUCe, Noise), pnorm, lower.tail=FALSE)) %>%
    select(Type, Name, SigOf, Generic, Class, nFeats, `F-statistic`, AUCe, Noise, AUC, pval)
write_csv(P, "output/summary.csv" )

P <- P %>% mutate(nlog10p = -log10(pval),
                  Highlight = ifelse(pval < 0.01, "yes", "no")) %>%
    group_by(Highlight) %>% mutate(Index = 1:n()) %>% ungroup() %>%
    mutate(Label=ifelse(Highlight == "yes", as.character(Index), ""))

vs <- c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005)
gg1 <- ggplot(P, aes(x=AUC, y=nlog10p, color=Highlight)) + theme_bw() +
    geom_point(size=0.5) +
    scale_y_continuous(breaks=-log10(vs), labels=as.character(vs), name="p value") +
    facet_wrap(~Type, nrow=1) +
    theme(strip.background=element_blank()) +
    scale_color_manual(values=c("yes"="red", "no"="black"), guide='none') +
    ggrepel::geom_text_repel(aes(label=Label), show.legend=FALSE, size=2, segment.size=0.2)

ggsave("plots/05-summary.png", gg1, width=6, height=3)

hmean <- function(x) {
    length(x) / sum( 1/x )
}

pal <- c("#F7F7F7", rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))[4:7])
##pal <- RColorBrewer::brewer.pal( 9, "YlOrBr" ) %>% colorRampPalette
HMP <- P %>% group_by( Name, Class ) %>%
    summarize(hmp = hmean(pval), .groups="drop") %>%
    mutate(nlogp = -log10(hmp),
           Label = ifelse(hmp < 0.05, as.character(round(hmp, 3)), ""))

vs <- c(0.5, 0.1, 0.05, 0.01, 0.005)
gg2 <- ggplot( HMP, aes(x=Name, y=Class, fill=nlogp) ) +
    theme_bw() + geom_tile() +
    geom_text(aes(label=Label)) +
    scale_fill_gradientn(colors=pal, breaks=-log10(vs),
                         labels=vs, name="HMP") +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave( "plots/05-hmp.png", gg2, width=9, height=6 )

