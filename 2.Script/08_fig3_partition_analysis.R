# ----------------------------------------------------------------------------
# Fig 3 | Partition and cell distribution of tdT subsets at d11 and d19
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/2.1.D11D19partitionAnalysis.R
# Input    : 1.Data/scrna/D11D19DatWithMonocleAnnot.RDS
# Output   : 3.Output/scrna/D11D19Dataset/cellDistributionInPartition.png
#            3.Output/scrna/D11D19Dataset/umap.png
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

datSub <- subset(dat, time %in% c('D11', 'D19')) 

monoClst <- readRDS('3.Output/scrna/3.Pseudotime/monoclePartitions.RDS')

datSub <- AddMetaData(datSub, monoClst[colnames(datSub)], col.name = 'monoclePart')

DimPlot(datSub, group.by = 'monoclePart', split.by = 'time', ncol = 2) +
  coord_fixed()

Idents(datSub) <- 'monoclePart'

deG <- wilcoxauc(datSub)

mkG <- cosg(datSub)

# P 1
# 
# Key markers: Sell, Il7r, Lef1, Satb1, Klf2, Ccr7, Foxp1, Tcf7
# 
# Hallmarks of naïve/central memory CD8 T cells — lymphoid homing (SELL, CCR7), survival (BCL2), transcriptional regulation for quiescence and differentiation (TCF7, LEF1, KLF2, FOXP1).
# Suggested name:
#   
#   Naïve/central memory CD8 program
# 
# P 2
# 
# Key markers: Slamf6, Itgb1, Itgae, Cxcr3, Ptger4, Ikzf2, Cd40lg, Tnfsf8, Ltb
# 
# Indicates activated/migratory CD8 T cells — costimulation (SLAMF6, CD40LG), adhesion/integrins (ITGB1, ITGAE), trafficking chemokine receptors (CXCR3), early activation TFs (IKZF2).
# Suggested name:
#   
#   Activated/migratory CD8 program
# 
# P 3
# 
# Key markers: Gzmb, Gzmk, Klrc1, Klrk1, Fasl, Prf1, Pdcd1, Tigit, Havcr2, Ctla4, Ccl5
# 
# Strong cytotoxic effector + exhaustion profile — cytolytic machinery (PRF1, GZMB, GZMK), NK-like receptors (KLRC1, KLRK1), inhibitory checkpoints (PDCD1, TIGIT, HAVCR2, CTLA4).
# Suggested name:
#   
#   Cytotoxic effector/exhausted CD8 program

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

mkG <- cosg(datSub)

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

DimPlot(datSub, group.by = 'anno', label = T) 
saveRDS(datSub, '1.Data/scrna/D11D19DatWithMonocleAnnot.RDS')

###################################### plotting #####################################################################
datSub <- readRDS('1.Data/scrna/D11D19DatWithMonocleAnnot.RDS')
datSub$timeGeno <- paste(datSub$time, datSub$genotype)
DimPlot(datSub, group.by = 'newclusters', label = T) +
  coord_fixed()

table(datSub$timeGeno, datSub$monoclePart) %>% 
  data.frame() %>% 
  separate(Var1, into = c('time', 'genotype'), sep = ' ') %>% 
  dplyr::rename(part = Var2) %>% 
  group_by(part, time) %>% 
  mutate(ratio = 100 * Freq/sum(Freq),
         part = recode(part,
                       `1` = "Naïve/central memory",
                       `2` = "Activated/migratory",
                       `3` = "Cytotoxic effector/exhausted"),
         genotype = factor(genotype,
                           levels = c("CD44+tdT-",
                                      "CD44+tdT+", 
                                      "DN", 
                                      "LAG3-tdT+",
                                      "LAG3+tdT+"))) %>% 
  ggplot(aes(x = time, y = ratio, fill = genotype)) +
  geom_bar(stat = 'identity', position = "stack", width = .5, col = 'gray50') +
  geom_flow(aes(alluvium = genotype), alpha= .5,
            curve_type = "sigmoid", col = 'gray50', 
            # “linear”, “cubic”, “quintic”, “sine”, “arctangent”, “sigmoid”
            width = .5) + 
  scale_fill_manual(values = c("#A6A6A6", "#F4A6B0", "darkgreen", "#D73027", "#762A83")) +
  facet_grid(~ part) +
  coord_fixed(ratio = 1/25) +
  theme_bw() +
  theme(text = element_text(size = 8, face = 'bold', color = 'black'),
legend.position = 'none')

ggsave("3.Output/scrna/D11D19Dataset/cellDistributionInPartition.png", width = 8, height = 5)
graph2ppt(file="3.Output/scrna/D11D19Dataset/cellDistributionInPartition.pptx", 
          width = 6, height = 2.5)

plotDat <- data.frame(monoclePart = datSub$monoclePart, 
                      time = datSub$time, 
                      anno = datSub$anno,
                      Embeddings(datSub, 'umap')) 

plotDat %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  # 1) filled outer region per group
  stat_density_2d_filled(
    aes(fill = monoclePart,
        col = monoclePart),
    contour_var = "ndensity",
    breaks = c(0.05, Inf),   # one region: ndensity >= 0.05 (tweak this)
    n = 400,                 # higher grid to avoid blocky edges
    adjust = 0.75,            # smooth a bit; tweak per your data
    alpha = 0.25,
    show.legend = FALSE
  ) +
  geom_point(aes(col = anno), size = 0.5, alpha = 0.5) + 
  geom_label(data = plotDat %>% 
               group_by(anno) %>% 
               summarise(UMAP_1 = median(UMAP_1),
                         UMAP_2 = median(UMAP_2)),
             aes(UMAP_1, UMAP_2, label = gsub('\\:.*', '', anno), 
                 fill = anno), 
             alpha = 0.3,
             show.legend = FALSE) +
  scale_fill_manual(values = c('#cabfe7', '#d488f6', '#ad07f9', 
                               scales::hue_pal()(15))) +
  scale_color_manual(values = c('#cabfe7', '#d488f6', '#ad07f9', 
                                scales::hue_pal()(15))) +
  xlim(-9, 10) +
  ylim(-14, 11.5)  +
  coord_fixed()

p <- ggplot()

labs15 <- c(
  "C1: Precursor T cells","C2: Stem-like memory","C3: Central memory","C4: Long-lived memory",
  "C5: Biosynthetic","C6: IFN-induced","C7: Effector memmory","C8: Memory precursor",
  "C9: Co-stimulating","C10: SLAMF6 progenitor-memory","C11: Activated","C12: Late activated",
  "C13: Effector","C14: Terminally Exhausted","C15: Proliferative Exhausted"
)

for (i in labs15) {
  p <- p + geom_point(data = plotDat |> filter(anno == i), 
                      aes(UMAP_1, UMAP_2, fill = anno), col = 'black',
                    shape = 21, size = 1.5, stroke = 1) +
           geom_point(data = plotDat |> filter(anno == i), 
                      aes(UMAP_1, UMAP_2, col = anno),
                    shape = 16, size = 1)
}

# palette for 15 + extras (named so mapping is stable)
pal <- c(scales::hue_pal()(15), '#cabfe7', '#d488f6', '#ad07f9')
names(pal) <- c(labs15, "1", "2", "3")

p +
  # 1) filled outer region per group
  stat_density_2d_filled(
    data = plotDat,
    aes(UMAP_1, UMAP_2, 
      fill = monoclePart,
        col = monoclePart),
    contour_var = "ndensity",
    breaks = c(0.05, Inf),   # one region: ndensity >= 0.05 (tweak this)
    n = 400,                 # higher grid to avoid blocky edges
    adjust = 0.75,            # smooth a bit; tweak per your data
    alpha = 0.25,
    show.legend = FALSE
  ) +
  # geom_label(data = plotDat %>% 
  #              group_by(anno) %>% 
  #              summarise(UMAP_1 = median(UMAP_1),
  #                        UMAP_2 = median(UMAP_2)),
  #            aes(UMAP_1, UMAP_2, label = gsub('\\:.*', '', anno), 
  #                fill = anno), 
  #            alpha = 0.3,
  #            show.legend = FALSE, fontface = 2, size = 4) +
  scale_fill_manual(values = pal,
                    breaks = labs15,
                    drop = F) +
  scale_color_manual(values = pal,
                    breaks = labs15,
                    drop = F) +
  guides(
  fill = guide_legend(override.aes=list(shape = 21, 
                                              size = 5, 
                                              exclude = c("1","2","3")))
) + 
  labs(fill = '', color = '') +
  xlim(-9, 10) +
  ylim(-14, 11.5)  +
  coord_fixed() +
  theme_void() +
  theme(text = element_text(size = 15, face = 'bold', color = 'black'),
legend.position = 'none') 

ggsave('3.Output/scrna/D11D19Dataset/umap_nolabel.png', width = 8, height = 8)

########################
library("ggdensity")

labDf <- plotDat %>% 
  group_by(anno) %>% 
  summarise(UMAP_1 = mean(UMAP_1),
            UMAP_2 = mean(UMAP_2))

# one clean outline per remaining cluster
ggplot(plotDat, aes(UMAP_1, UMAP_2))  +
  geom_hdr() +
  stat_density_2d(data = plotDat_filt, 
    aes(UMAP_1, UMAP_2, group = anno, color = anno),
    contour_var = "ndensity",
    lty = 5,
    lwd = 0.3, 
    breaks = c(0.1),   # outer envelope; tweak 0.01–0.05
    adjust = 1.2,
    show.legend = FALSE
  ) +
  # geom_text(data = labDf, 
  #   aes(label = gsub(':.*', '', anno), col = anno),
  #   fontface = 2, size = 3) +
  facet_grid(~time) +
  xlim(-9, 10) +
  ylim(-15, 11.5) +
  coord_fixed() +
  theme_void() +
  theme(text = element_text(size = 12, face = 2),
        legend.position = 'none')

graph2ppt(file="3.Output/scrna/D11D19Dataset/cellDistributionDensity.pptx", 
          width = 6, height = 3)
