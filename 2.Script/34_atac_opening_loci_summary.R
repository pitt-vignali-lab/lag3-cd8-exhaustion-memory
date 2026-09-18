# ----------------------------------------------------------------------------
# Fig 4 | Summarise loci opening and closing between subsets
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/openingLociSum.R
# Input    : 1.Data/atac/expressionMatrixBatchEffFixed.RDS
# Output   : 3.Output/atac/
# Upstream : 23_atac_de_peaks_for_viz.R
# ----------------------------------------------------------------------------

library(tidyverse)
library(ggrepel)

annoDf <- data.table::fread('1.Data/atac/MergedReplicate/consensus_peaks.mRp.clN.annotatePeaks.txt') %>% 
  rename("interval_id" = 1)
boolDf <- data.table::fread('1.Data/atac/MergedReplicate/consensus_peaks.mRp.clN.boolean.txt')

supDf <- merge(annoDf[, c(1:4, 9:11)],
               boolDf[, c(4, 7:14)], 
               by = 'interval_id') %>% 
  select(-interval_id)

names(supDf) <- gsub('.mRp.clN.bool', '_Open', names(supDf))
write.csv(supDf, '3.Output/atac/allOpenningRegion.csv')

condis <- names(boolDf) 
condis <- condis[- grep('LAG3', condis)]  
condis <- condis[grep('bool', condis)]

# [grep('bool', names(boolDf))][-grep('LAG3', names(boolDf))]
  
allcondiDf <- data.frame()

for (i in condis) {
  subBool <- boolDf %>% 
    dplyr::select("interval_id", i) %>% 
    rename('readout' = 2) %>% 
    mutate(condi = gsub('.mRp.clN.bool', '', i)) %>%
    merge(annoDf %>% select("interval_id", 'Gene Name'),
          by = "interval_id") %>% 
    filter(readout) %>% 
    group_by(condi, `Gene Name`) %>% 
    summarise(n = n()) %>% 
    rename('GeneName' = 2)
    
   assign(paste0('p', gsub('.mRp.clN.bool', '', i)),
          ggplot(subBool %>% 
                 filter(n > 10),
               aes(x = reorder(GeneName, -n), y = n)) +
     geom_point())
     
    allcondiDf <- rbind(subBool, allcondiDf)
}

allcondiDf <- allcondiDf %>% 
  group_by(condi) %>% 
  arrange(desc(n)) %>% 
  mutate(ord = row_number())

labDf <- allcondiDf %>% 
  group_by(condi) %>% 
  top_n(50, -ord)

allcondiDf %>% 
  filter(ord <= 25) %>% 
  
ggplot(aes(y = reorder(ord, n), x = n)) +
  # geom_bar(stat = 'identity') +
  coord_cartesian(clip = "off") +
  geom_text(aes(label = GeneName),
            hjust = 0,
            nudge_x = 2, 
            size = 5) +
  xlim(10, 75) +
  geom_point(col = 'blue') +
  facet_grid(~factor(condi, levels = gsub('.mRp.clN.bool', '', condis)[c(3:4, 1:2, 5, 7, 6)])) +
  # coord_flip() +
  labs(x = 'Number of accessible chromatin regions',
       y = '') +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 15),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line.x = element_line()) 


ggsave('3.Output/atac/top25OpenningRegion.png', width = 20, height = 5, units = 'in')

