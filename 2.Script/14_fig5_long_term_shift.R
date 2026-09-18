# ----------------------------------------------------------------------------
# Fig 5 | Long-term shift of tdT subsets in the periphery
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/4.1.longTermShift.R
# Input    : 1.Data/scrna/
# Output   : 3.Output/scrna/D11D19D60CD44.png
# Upstream : 07_scrna_clustering_final.R
# ----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(cowplot)
library(vegan)
library(presto)
library(COSG)
library(ggalluvial)
library(ggforce)
library(export)

dat <- readRDS('1.Data/scrna/clstedDatRes08FilteredWithMonocleAnnot.RDS')

datSub <- subset(dat, time %in% c('D11', 'D19', 'D60') & genotype == 'CD44+tdT+') 

datSub <- subset(datSub, idents = 5, invert = T)

datSub$newclusters <- recode(datSub$seurat_clusters,
                             `1` = 2, `2` = 4, `3` = 8, `4` = 3, 
                             `5` = 13, `6` = 5, `7` = 1, `8` = 9, 
                             `9` = 14, `10` = 10, `11` = 11, 
                             `12` = 7, `13` = 15, `14` = 6, 
                             `15` = 12)
unique( paste0('C', datSub$newclusters, ': ', 
               gsub('C.*: ', '', datSub$anno)))

datSub$anno <- paste0('C', datSub$newclusters, ': ', 
                      gsub('C.*: ', '', datSub$anno))

DimPlot(datSub, group.by = 'anno', label = T) 

Idents(datSub) <- 'newclusters'

datSub$anno <- recode(datSub$anno,
                      "C2: naïve" = 'C2: Stem-like memory', 
                      "C8: Stem-like memory" = 'C8: Memory precursor', 
                      "C5: Antigent-presenting" = 'C5: Biosynthetic',
                      "C4: Stress response" = "C4: Long-lived memory",
                      "C10: Long-lived memory" = "C10: SLAMF6 progenitor-memory"
                    ) 

datSub$anno <- factor(datSub$anno, 
                      levels = c("C1: Precursor T cells",
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
                                 "C15: Proliferative Exhausted"))

DimPlot(datSub, group.by = 'anno', label = T) +
  coord_fixed()

table(datSub$anno, datSub$time) %>% 
  data.frame() %>% 
  group_by(Var2) %>% 
  mutate(ratio = 100 * Freq/sum(Freq)) %>% 
  ggplot(aes(x = Var2, y = ratio, fill = Var1)) +
  geom_bar(stat = 'identity', color = 'white', width = .5, alpha= .7) + 
  geom_flow(aes(alluvium = Var1), alpha= .7,
            curve_type = "sigmoid", col = 'white', 
            # “linear”, “cubic”, “quintic”, “sine”, “arctangent”, “sigmoid”
            width = .5) + 
  # scale_fill_manual(values = c("#A6A6A6", "#F4A6B0", "darkgreen", "#D73027", "#762A83")) +
  # facet_grid(~ part) +
  coord_fixed(ratio = 1/25) +
  theme_bw() +
  theme(text = element_text(size = 8, face = 'bold', color = 'black'),
        legend.position = 'bottom')

table(datSub$anno, datSub$time) %>% 
  data.frame() %>% 
  group_by(Var2) %>% 
  mutate(ratio = 100 * Freq/sum(Freq)) %>% 
ggplot(aes(x = as.numeric(gsub('D', '', Var2)), y = ratio,  alluvium = Var1)) +
  geom_alluvium(aes(fill = Var1), col = 'white', 
                   width = 0.6 ,
                   alpha = 2/3, decreasing = NA,
                   curve_type = "sigmoid") +
  labs(y = 'Percentage (%)',
       x = 'Time',
      fill = 'Clusters') +
  scale_x_continuous(expand=c(0,0), breaks = c(11, 19, 60), labels = c('D11', 'D19', 'D60')) +
  scale_y_continuous(expand=c(0,0)) +
    theme( # remove the vertical grid lines
           panel.grid.major.x = element_blank() ,
           legend.position = 'bottom', 
           # explicitly set the horizontal lines (or they will disappear too)
           panel.grid.major.y = element_blank() ,
           text = element_text(size = 6, face = 2)
    )

graph2ppt(file = '3.Output/scrna/D11D19D60CD44.pptx', width = 12, height = 10)

table(datSub$anno, datSub$time) %>% 
  data.frame() %>% 
  group_by(Var2) %>% 
  mutate(ratio = 100 * Freq/sum(Freq)) %>% 
ggplot(aes(x = as.numeric(gsub('D', '', Var2)), y = ratio,  alluvium = Var1)) +
  geom_alluvium(aes(fill = Var1), col = 'white', 
                   width = 0.6 ,
                   alpha = 0.8, decreasing = NA,
                   curve_type = "sigmoid") +
  labs(x = '', y = '') +
  scale_x_continuous(expand=c(0,0), breaks = c(11, 19, 60), labels = c('D11', 'D19', 'D60')) +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text = element_blank(),
        strip.text.y = element_blank(),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.grid = element_blank(), 
      legend.position = 'none')

ggsave(file = '3.Output/scrna/D11D19D60CD44.png', width = 3, height = 3.5)
