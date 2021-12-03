library(tidyverse)

M <- read_csv('data/agents-metadata.csv', col_types=cols()) %>%
    select(Drug = agent, Class = drug_class )

PMID <- jsonlite::read_json('data/pmids_for_drugs.json') %>%
    keep(~length(.x) > 0) %>%
    modify_depth(2, modify_if, is.null, ~as.integer(NA)) %>%
    modify_depth(2, as_tibble) %>%
    map(enframe, "PMID") %>% map(unnest, value) %>%
    enframe("Drug") %>% unnest(value) %>%
    select(Drug, PMID, Year=year, nCite = citation_count,
           PubChem = pubchem_support, MeSH = mesh_support,
           Grounding = grounding_support) %>%
    left_join(M, by="Drug")

plotCounts <- function() {
    X <- PMID %>% group_by(Drug) %>%
        summarize(nPMID = length(PMID),
                  Class = unique(Class)) %>%
        arrange(Class, desc(nPMID)) %>%
        mutate(Drug = factor(Drug,Drug),
               lbl = ifelse(nPMID >= 100, as.character(nPMID), ""),
               lbl = str_c(lbl, "  "),
               nPMID = pmin(nPMID,100))

    pal <- c("white", ggthemes::few_pal()(8), "darkgray")
    ggplot(X, aes(x=Drug, y=nPMID)) +
        theme_minimal() +
        geom_bar(stat='identity', color='gray', aes(fill=Class)) +
        geom_text(aes(label=lbl), angle=90, size=3, hjust=1) +
        scale_y_continuous(breaks=seq(0,100,by=25),
                           labels=c("0","25","50","75","100+"),
                           name="# Papers") +
        scale_fill_manual(values=pal) +
        xlab("Agent") +
        theme(panel.grid.major.x = element_blank(),
              axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
        ggsave("plots/pmid-count.png", width=10, height=4)
}
