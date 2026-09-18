# ----------------------------------------------------------------------------
# Fig 3 | Marker gene expression across tdT subsets
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/2.2.D11D19MarkerGeneVisulization.R
# Input    : 1.Data/scrna/D11D19DatWithMonocleAnnot.RDS
# Output   : 3.Output/scrna/D11D19MarkerGenes*.png
# Upstream : 08_fig3_partition_analysis.R
# ----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(cowplot)
library(export)
library(COSG)

dat <- readRDS('1.Data/scrna/D11D19DatWithMonocleAnnot.RDS')

Idents(dat) <- 'anno'

cosgDf <- cosg(dat)

write.csv(cosgDf$names, '3.Output/scrna/D11D19Dataset/cosgGeneDf.csv')

mkDf <- rbind(
  data.frame(id = 'C1: Precursor T cells',
             marker = c('Tcf12')),
  data.frame(id = 'C2: Stem-like memory',
             marker = c('Sell', 'Lef1', 'Ccr7')),
  data.frame(id = 'C3: Central memory',
             marker = c('Bcl2', 'Cxcr3', 'Sell')),
  data.frame(id = 'C4: Long-lived memory',
             marker = c('Lef1', 'Klf2', 'Ccr7', 'Sell')),
  data.frame(id = 'C5: Biosynthetic',
             marker = c('Myc','Fbl', 'Ppan')),
  data.frame(id = 'C6: IFN-induced',
             marker = c('Ifit1','Ifit3','Irf7')),
  data.frame(id = 'C7: Effector memmory',
             marker = c('Cxcr3','Ccl5','Gzmk')),
  data.frame(id = 'C8: Memory precursor',
             marker = c('Slamf6', 'Cd9', 'Itgb1')), # SLAMF6 marks stem-like/precursor CD8 states; CD9 and ITGB1 capture the adhesion/tissue-interaction/TRM-primed
  data.frame(id = 'C9: Co-stimulating',
             marker = c('Tnfrsf4','Tnfrsf18','Tnfrsf25')),
  data.frame(id = 'C10: SLAMF6 progenitor-memory',
             marker = c('Tcf7', 'Slamf6', 'Id3')), # stem/memory programming (Tcf7/Id3) with quiescent longevity/maintenance (Myb).
  data.frame(id = 'C11: Activated',
             marker = c('Mcm2','Mcm3','Mcm5','Mcm7','E2f1')),
  data.frame(id = 'C12: Late activated',
             marker = c('Bcl3','Socs3','Elovl6')),
  data.frame(id = 'C13: Effector',
             marker = c('Gzmb','Gzma','Gzmk','Ccl5','Fasl')),
  data.frame(id = 'C14: Terminally Exhausted',
             marker = c('Tox','Pdcd1','Lag3','Havcr2','Tnfrsf9')),
  data.frame(id = 'C15: Proliferative Exhausted',
             marker = c('Tox','Mki67','Top2a','Cdk1'))
)


plotDf <- DotPlot(dat, features = unique(mkDf$marker))$data %>% 
  mutate(id = gsub('C4: Stress response', 'C4: Long-lived memory', id)) %>% 
  mutate(id = gsub('C10: Long-lived memory', 'C10: SLAMF6 progenitor-memory', id)) %>% 
  merge(mkDf, by.x = 'features.plot', by.y = 'marker') %>% 
  mutate(id.x = factor(gsub(':.*', '', id.x),
                       levels = c(paste0('C', 1:15))),
         id.y = factor(id.y, levels = c("C1: Precursor T cells",
                                 "C2: Stem-like memory",
                                 "C3: Central memory",
                                 "C4: Long-lived memory",
                                 "C5: Biosynthetic",
                                 "C6: IFN-induced",
                                 "C7: Effector memmory",
                                 "C8: Memory precursor",
                                 "C9: Co-stimulating",
                                 "C10: SLAMF6 progenitor-memory",
                                 "C11: Activated",
                                 "C12: Late activated",
                                 "C13: Effector",
                                 "C14: Terminally Exhausted",
                                 "C15: Proliferative Exhausted")))

mkp <- plotDf %>% 
  ggplot(aes(x = id.x, y = features.plot)) +
  geom_point(aes(size = pct.exp, 
                 fill = avg.exp.scaled),
             shape = 21, stroke = 0.1, col = 'black') +
  scale_size_continuous(range = c(0.5, 2)) +
  scale_fill_gradient2(high = 'tomato', low = 'blue', mid = 'white', midpoint = 0) +
  facet_grid(id.y~., scales = "free", space = "free") +
  theme_bw() +
  labs(x = "clusters",
       y = 'Marker genes',
       fill = 'Average Expression',
       size = 'Percent Expression') +
  theme(strip.text.y = element_text(angle = 0),
        panel.spacing = unit(0.1, "lines"), 
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        text = element_text(size = 8, face = 2),
        legend.position = 'bottom')

mkp
ggsave('3.Output/scrna/D11D19MarkerGenes2.png', width = 8, height = 12, units = 'in')         
graph2ppt(mkp, '3.Output/scrna/D11D19MarkerGenes.pptx', width = 3.5, height = 6)

mkp +
  labs(x = '', 
y = '') +
  coord_fixed() +
  theme(axis.text = element_blank(),
        strip.text.y = element_blank(),
        plot.background = element_rect(fill = 'transparent', color = NA),
      legend.position = 'none')
ggsave('3.Output/scrna/D11D19MarkerGenes_clen.png', width = 3, height = 5.2, units = 'in',
       bg = "transparent")  
