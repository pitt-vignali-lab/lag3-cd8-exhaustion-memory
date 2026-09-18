# ----------------------------------------------------------------------------
# Analysis | Render the tumour-set velocity summary plots
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/veloTumorSet/plotDf.R
# Input    : 2.Script/veloTumorSet/
# Output   : 3.Output/scrna/
# Upstream : 18_velocity_tumor_prepare.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(COSG)
library(ggpubr)

obj <- readRDS('integrated_tumorSet.rds')

obj <- FindClusters(obj, resolution = 0.35, cluster.name = "tumor_clusters")

table(obj$tumor_clusters, obj$anno)


DimPlot(obj, group.by = "tumor_clusters", alpha = 0.5, label = T) +
  coord_fixed() +
  theme_void()


mkG <- cosg(JoinLayers(obj))

# Cluster 0: TCF7/IL7R/KLF2 signature → early effector–memory precursors
# Cluster 1: MKI67/Top2a/Ccnb genes → cycling/proliferating
# Cluster 2: metabolic/translation-intermediate → “Intermediate Activation CD8” is accurate
# Cluster 3: LAG3/OX40/4-1BB/TIGIT/HAVCR2 → exhausted transition → “Exhausted” is valid
# Cluster 4: GZMB/PRF1/IFN-high → “Cytotoxic Effector” is spot-on

tumor_anno <- c(
  '0' = 'Memory Precursors',
  '1' = 'Proliferating',
  '2' = 'Intermediate Activation',
  '3' = 'Exhausted',
  '4' = 'Cytotoxic Effector'
)
t <- data.frame(tumor_anno = tumor_anno[as.character(obj$tumor_clusters)])
rownames(t) <- names(obj$tumor_clusters)

obj <- AddMetaData(obj, t)

DimPlot(obj, group.by = "tumor_anno", alpha = 0.5, label = T) +
  coord_fixed() +
  theme_void()

plotDf <- Embeddings(obj, 'umap') %>%
  data.frame(genotype = obj$genotype, anno = obj$tumor_anno)

plotDf %>%
  ggplot(aes(x = umap_1, y = umap_2, col = genotype)) +
  geom_point(size = 0.5) +
  scale_color_manual(
    values = c('LAG3+tdT+' = '#6a5acd', 'LAG3-tdT+' = '#e74c3c')
  ) +
  coord_fixed() +
  theme_void()

ggsave('umap_genotype.png', width = 4, height = 4)

ggplot(plotDf, aes(x = umap_1, y = umap_2, col = anno)) +
  geom_point(alpha = 1, size = 0.5) +
  scale_color_brewer(palette = 'Set2') +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom')

ggsave('umap_anno.png', width = 6, height = 4)

table(anno = obj$tumor_anno, genotype = obj$genotype, rep = obj$rep) %>%
  data.frame() %>%
  group_by(rep, genotype) %>% # Calculate proportion within each rep AND genotype
  mutate(Freq = Freq / sum(Freq)) %>%
  ggplot(aes(
    x = factor(
      gsub(' ', ' \n', anno),
      levels = c(
        'Memory \nPrecursors',
        'Proliferating',
        'Intermediate \nActivation',
        'Exhausted',
        'Cytotoxic \nEffector'
      )
    ),
    y = Freq,
    col = genotype
  )) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
    size = 0.6
  ) +
  stat_compare_means(aes(label = after_stat(p.signif))) +
  scale_color_manual(
    values = c('LAG3+tdT+' = '#6a5acd', 'LAG3-tdT+' = '#e74c3c')
  ) +
  labs(x = NULL, y = "Proportion within genotype") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 5, face = 2, color = 'black'),
    legend.position = 'none'
  )

export::graph2ppt(
  file = 'genotype_percentage_by_rep.pptx',
  width = 2,
  height = 1.6
)
