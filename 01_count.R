library(tidyverse)

#M <- read_csv('data/agents_metadata.csv', col_types=cols()) %>%
#    select(Drug = agent, Class = drug_class )
#
#PMID <- jsonlite::read_json('data/pmids_for_drugs.json') %>%
#    keep(~length(.x) > 0) %>%
#    modify_depth(2, modify_if, is.null, ~as.integer(NA)) %>%
#    modify_depth(2, as_tibble) %>%
#    map(enframe, "PMID") %>% map(unnest, value) %>%
#    enframe("Drug") %>% unnest(value) %>%
#    select(Drug, PMID, Year=year, nCite = citation_count,
#           PubChem = pubchem_support, MeSH = mesh_support,
#           Grounding = grounding_support) %>%
#    mutate(Evidence = PubChem + MeSH + Grounding) %>%
#    left_join(M, by="Drug")
#
#write_csv(PMID, 'data/pmids_for_drugs_tidy.csv')

PMID <- read_csv("data/pmids_for_drugs_tidy.csv")

plotCounts <- function(.df, nThresh=100) {
    X <- .df %>% group_by(Drug) %>%
        summarize(nPMID = length(PMID),
                  Class = unique(Class)) %>%
        arrange(Class, desc(nPMID)) %>%
        mutate(Drug = factor(Drug, Drug),
               lbl = ifelse(nPMID >= nThresh, as.character(nPMID), ""),
               lbl = str_c(lbl, "  "),
               nPMID = pmin(nPMID, nThresh))

    pal <- c("white", ggthemes::few_pal()(8), "darkgray") %>%
        set_names(unique(PMID$Class))
    ggplot(X, aes(x = Drug, y = nPMID)) +
        theme_minimal() +
        geom_bar(stat = "identity", color = "gray", aes(fill = Class)) +
        geom_text(aes(label=lbl), angle=90, size=3, hjust=1) +
        scale_y_continuous(breaks=seq(0,100,by=25),
                           labels=c("0","25","50","75","100+"),
                           name="# Papers") +
        scale_fill_manual(values=pal) +
        xlab("Agent") +
        theme(panel.grid.major.x = element_blank(),
              axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
}

plotCounts(PMID) %>%
    ggplot2::ggsave("plots/01-pmid-count.png", ., width = 10, height = 4)
filter(PMID, Evidence == 3) %>%
    plotCounts(50) %>%
    ggplot2::ggsave("plots/01-pmid-count3.png", ., width = 5, height = 4)

plotCitations <- function() {
    X <- PMID %>% filter(!is.na(nCite)) %>%
        select(PMID, Year, nCite) %>%
        distinct() %>%
        mutate( Label = ifelse(nCite>500, PMID, "") )

    gg1 <- ggplot(X, aes(x=Year, y=nCite)) +
        theme_minimal() + geom_point() +
        ggrepel::geom_text_repel(aes(label=Label), nudge_y=50) +
        ylab("Number of citations")
    ggplot2::ggsave("plots/01-cite-count.png", gg1, width=8, height=6)

    nTop <- c("2021"=0, "2020"=1, "2019"=1, "2018"=3, "2017"=3, "2016"=5,
              "2015"=2, "2014"=3, "2014"=3, "2013"=3, "2012"=6, "2011"=2, "2010"=3)
    X1 <- X %>% filter( Year >= 2010 ) %>%
        group_by(Year) %>%
        mutate(Rank = length(nCite)-rank(nCite)+1) %>%
        ungroup() %>%
        mutate(nCite = pmin(nCite,700),
               Label = ifelse(Rank <= nTop[as.character(Year)], PMID, ""),
               Highlight = ifelse(Label == "", "no", "yes"))
    
    gg2 <- ggplot(X1, aes(x=nCite, y=Year, color=Highlight)) +
        theme_minimal() + geom_point() +
        ggrepel::geom_text_repel(aes(label=Label), size=3, nudge_y=-0.25) +
        scale_y_continuous(breaks=seq(2010,2021),
                           labels=seq(2010,2021)) +
        scale_x_continuous(breaks=seq(0,700,by=100),
                           labels=c(as.character(seq(0,600,by=100)), "700+"),
                           name="Number of citations") +
        scale_color_manual(values=c("yes"="red", "no"="black"), guide='none') +
        theme(panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank())

    ggplot2::ggsave("plots/01-cite-count-2010+.png", gg2, width=8, height=6)    
}

plotCitations()

lit <- PMID %>% mutate(PMID = str_c("PMID:", PMID))

# Count the number of unique drugs mentioned in each PMID
unq_cnt <- group_by(lit, PMID) %>%
    summarize(nDrug = length(unique(Drug)), .groups = "drop")

# Top 50 papers for each drug
res <- inner_join(lit, unq_cnt, by = "PMID") %>%
    select(Drug, PMID, nDrug, Evidence, Year, nCite) %>%
    group_by(Drug) %>%
    arrange(nDrug, desc(Evidence), desc(Year), desc(nCite)) %>%
    slice_head(n = 50) %>%
    ungroup()

write_csv(res, "output/lit-cand.csv")