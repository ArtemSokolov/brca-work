library(tidyverse)

GC <- jsonlite::read_json('data/gene_counts.json') %>%
    map(enframe, 'Gene', 'Count') %>%
    enframe('PMID', 'Counts') %>%
    unnest(Counts) %>%
    unnest(Count)

N <- GC %>% group_by(Gene) %>%
    summarize(Total   = sum(Count),
              nPapers = length(Count))

ggplot(N, aes(x=nPapers, y=Total)) +
    theme_bw() + geom_point() +
    scale_y_log10(name="Total # of mentions") +
    scale_x_log10(name="# papers with at least 1 mention") +
    ggsave("total-counts.png")

##GC %>% group_by(Gene) %>% mutate(nPapers = length(Count)) %>%
## ungroup %>% filter( nPapers == 1 ) %>% arrange(desc(Count))
