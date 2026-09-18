# ----------------------------------------------------------------------------
# Setup | Batch correction by reciprocal PCA and Harmony
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/1.3.batch_effects_removal.R
# Input    : 1.Data/scrna/datInt.RDS
# Output   : 1.Data/scrna/datIntegRpca.RDS
#            1.Data/scrna/datIntegHarmony.RDS
# Upstream : 02_scrna_qc_batch_detect.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)

dat <- readRDS("1.Data/scrna/datInt.RDS")

# Harmony ----
library(harmony)

datInt <- dat %>%
  RunHarmony(group.by.vars = "batch2", plot_convergence = TRUE)

# UMAP and clustering with harmonized PCs
datInt <- RunUMAP(datInt, reduction='harmony', dims = 1:30)
DimPlot(datInt, split.by = "batch2") +
  coord_fixed(ratio = 1)

saveRDS(datInt, "1.Data/scrna/datIntegHarmony.RDS")

# figures ----

plt <- DimPlot(datInt, label = T, label.size = 8) + 
  coord_fixed(ratio = 1) +
  NoLegend()

# x-range
round(diff(layer_scales(plt)$x$range$range)/2, 1)
# y-range
round(diff(layer_scales(plt)$y$range$range)/2, 1)

plotDf <- data.frame(datInt$batch2, datInt$genotype, 
                     Embeddings(datInt, reduction = "umap"))

names(plotDf) <- c('batch', 'genotype', 'x', 'y')

ggplot(plotDf, aes(x, y)) +
  geom_point(col = "gray20") +
  coord_fixed(ratio = 1) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("3.Output/scrna/1.DatInit/harmonyOverall.png", width = 9.3, height = 11.3,
       units = "in", dpi = 300, bg = "white")

ggplot(plotDf, aes(x, y, col = batch)) +
  geom_point() +
  facet_grid(~ batch)+
  coord_fixed(ratio = 1) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("3.Output/scrna/1.DatInit/harmonyBatch.png", width = 14, height = 8.5,
       units = "in", dpi = 300, bg = "white")


ggplot(plotDf, aes(x, y, col = genotype)) +
  geom_point()+
  coord_fixed(ratio = 1) +
  facet_grid(batch ~ genotype)+
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("3.Output/scrna/1.DatInit/harmonyBatchGenotype.png", width = 27.9, height = 11.3,
       units = "in", dpi = 300, bg = "white")

DimPlot(datRpca, group.by = "batch2") +
  coord_fixed() +
DimPlot(datRpca, group.by = "genotype")+
  coord_fixed()
ggsave("3.Output/scrna/1.DatInit/batch_effects_removal_Rpca.png",
       width = 9, height = 4, units = "in", dpi = 300, bg = "white")

DimPlot(datHarm, group.by = "rep") +
  coord_fixed() +
  DimPlot(datHarm, group.by = "genotype")+
  coord_fixed()
ggsave("3.Output/scrna/1.DatInit/batch_effects_removal_Harmony.png",
       width = 9, height = 4, units = "in", dpi = 300, bg = "white")


# Tested but not used ----
## Reciprocal PCA (RPCA) with CD44-CD62L+tdT- group 
dat$batch3 <- paste0(dat$batch2, dat$genotype)

datList <- SplitObject(dat, split.by = "batch3")

datList <- lapply(datList, function(x){
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = datList)

datList <- lapply(X = datList, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = datList, 
                                  reference = c(6, 12), 
                                  reduction = "rpca",
                                  dims = 1:50)
datInt <- IntegrateData(anchorset = anchors, dims = 1:50)

datInt <- ScaleData(datInt, verbose = FALSE)

datInt <- RunPCA(datInt, verbose = FALSE)

datInt <- RunUMAP(datInt, dims = 1:50, verbose = FALSE)

saveRDS(datInt, "1.Data/scrna/datIntegRpcaWithNaive.RDS")

# Reciprocal PCA (RPCA) 
datList <- SplitObject(dat, split.by = "batch2")

datList <- lapply(datList, function(x){
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = datList)

datList <- lapply(X = datList, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = datList, 
                                  # reference = c(1, 2), 
                                  reduction = "rpca",
                                  dims = 1:50)

datInt <- IntegrateData(anchorset = anchors, dims = 1:50)

datInt <- ScaleData(datInt, verbose = FALSE)

datInt <- RunPCA(datInt, verbose = FALSE)

datInt <- RunUMAP(datInt, dims = 1:50, verbose = FALSE)

saveRDS(datInt, "1.Data/scrna/datIntegRpca.RDS")

