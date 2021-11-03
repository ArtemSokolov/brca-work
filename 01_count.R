library(tidyverse)

PMID <- jsonlite::read_json('data/pmids_for_drugs.json') %>%
    map(unlist) %>% enframe(name='agent', value='pmids')

dcmap <- c(Chemotherapy="Chemo",
           `BCL2 family`="BCL2")

X <- read_csv('data/agents_metadata.csv', col_types=cols()) %>%
    select( -library_plate, -nominal_target ) %>%
    left_join(PMID, by='agent') %>%
    mutate(n = map_int(pmids, length))

Z <- X %>%
    arrange(desc(n)) %>%
    mutate(agent = factor(agent,agent),
           drug_class = recode(drug_class, !!!dcmap),
           lbl = ifelse(n >= 100, as.character(n), ""),
           lbl = str_c(lbl, "  "),
           n = pmin(n,100))

fplot <- function(ZZ) {
    ggplot(ZZ, aes(x=agent, y=n)) +
        theme_minimal() +
        geom_bar(stat='identity', color='gray', fill='white') +
        geom_text(aes(label=lbl), angle=90, size=3, hjust=1) +
        facet_grid(~drug_class, scales='free', space='free') +
        scale_y_continuous(breaks=seq(0,100,by=25),
                           labels=c("0","25","50","75","100+"),
                           name="# Papers") +
        xlab("Agent") +
        theme(panel.grid.major.x = element_blank(),
              axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
}

fplot(Z) + ggsave("pmid-count.png", width=14, height=4)
fplot(filter(Z, n>0)) + ggsave("pmid-count-nz.png", width=11, height=4)

Yc <- PMID %>% unnest(pmids) %>% count(pmids)
Y <- PMID %>% unnest(pmids) %>% inner_join(Yc, by="pmids") %>%
    mutate_at("pmids", as.integer)

