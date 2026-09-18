# ----------------------------------------------------------------------------
# Fig 4 | Genomic distribution of the opening loci
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/openingLociSumPie.R
# Input    : 1.Data/atac/MergedReplicate/consensus_peaks.mRp.clN.annotatePeaks.txt
# Output   : 3.Output/atac/
# Upstream : 34_atac_opening_loci_summary.R
# ----------------------------------------------------------------------------

library(tidyverse)
library(ggrepel)
library(export)

df <- read.table('1.Data/atac/MergedReplicate/mrep_peak_annotation-plot-1.tsv',
                 header = T) %>% 
  pivot_longer(names_to = 'region', cols = 2:6) %>% 
  group_by(X.Category.) %>% 
  # summarise(sum = sum(value))
  mutate(region = factor(region, levels = c('Intergenic', 'TTS', 'exon', 'intron', 'promoter.TSS')),
         perc = value/sum(value))


df %>% 
  ggplot() +
  geom_bar(aes(x = perc, y = X.Category., fill = region),
           width = 1) +
  coord_polar(theta = "x") +
  facet_grid(~X.Category.)

df %>% 
  filter(X.Category. != "TUMOR_LAG3_SP") %>% 
  mutate(X.Category. = factor(X.Category., 
                              levels = c('NDLN_DN', 'NDLN_TDT_SP',
                                         'DLN_DN', 'DLN_TDT_SP',
                                         'TUMOR_DN', "TUMOR_TDT_SP", 'TUMOR_DP'))) %>% 
  arrange(region) %>% 
  mutate(prop = value / sum(value) *100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop ) %>% 
  ggplot(aes(x="", y=prop, fill=rev(region))) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  facet_wrap(~X.Category., nrow = 2) +
  geom_text_repel(aes(y = ypos, label = paste0(round(prop, 2), '%')), color = "gray20", size=3) +
  theme_void()

graph2ppt(file='3.Output/atac/PieRegion.pptx', width=8, height=5)
ggsave('3.Output/atac/PieRegion.png', width = 8, height = 5, units = 'in')
