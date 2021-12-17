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
    mutate(Name = str_c("PMID", PMID)) %>%
    select(Name, SigOf=Drug) %>%
    mutate(across(SigOf, recode, `fulvestrant+ribociclib` = "ribociclib"))

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
    mutate(Type = ifelse(Name == "PMID25501949", "Core", Type)) %>%
    unnest(AUC) %>% inner_join(MBK, by="Drug") %>%
    left_join(LM, by="Name")

## Compute the p values associated with each AUC measure
P <- X %>% mutate(AUCe  = nFeats*Slope + Bias,
                  pval  = pmap_dbl(list(AUC, AUCe, Noise), pnorm, lower.tail=FALSE)) %>%
    select(Type, Name, SigOf, Generic, Class, nFeats, `F-statistic`, AUCe, Noise, AUC, pval)
write_csv(P, "output/summary.csv" )

## Harmonic mean to aggregate p values for a single per-signature metric
hmean <- function(x) {
    length(x) / sum( 1/x )
}
HMP <- P %>% group_by(Name, SigOf, Type) %>%
    summarize( across(pval, hmean), .groups="drop" ) %>%
    mutate(Generic = "Harmonic\nmean p-val")

## Supplementary computations for the figure
PLT <- bind_rows(P, HMP, .id="Category") %>%
    mutate(nlogp = -log10(pval),
           Label = ifelse((pval < 0.05) | (Category == 2),
                          as.character(round(pval,3)), ""),
           Signature = ifelse(is.na(SigOf), Name, str_c(Name, "\n(", SigOf, ")")),
           Highlight = ifelse(SigOf == Generic, "yes", "no"))

## Plot all-by-all hits
pal <- c("#F7F7F7", rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))[4:7])
vs <- c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005)
gg <- ggplot(PLT, aes(x=Signature, y=Generic, fill=nlogp)) +
    theme_bw() + geom_tile() +
    geom_tile(data=filter(PLT, Highlight=="yes"), color="black", size=1) +
    geom_text(aes(label=Label), color="black") +
    scale_fill_gradientn(name="p-value", colors=pal,
                         labels=vs, breaks=-log10(vs)) +
    ylab("Drug") +
    facet_grid(Category~Type, space="free", scales="free") +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
          strip.text = element_blank(),
          strip.background = element_blank())
ggsave("plots/05-all.png", gg, width=9, height=10)
