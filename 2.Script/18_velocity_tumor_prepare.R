# ----------------------------------------------------------------------------
# Analysis | Export embeddings for the tumour-set velocity run
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/veloTumorSet/1.init.R
# Input    : 1.Data/scrna/
# Output   : 2.Script/veloTumorSet/
# Upstream : 07_scrna_clustering_final.R
# ----------------------------------------------------------------------------

library(Seurat)
library(scales)
library(tidyverse)
library(SeuratData)
library(SeuratWrappers)

dat <- readRDS("1.Data/scrna/D11D19DatWithMonocleAnnot.RDS")

datSub <- subset(
  dat,
  tissue == 'Tumor' & genotype %in% c('LAG3-tdT+', 'LAG3+tdT+')
)

DimPlot(datSub, group.by = "genotype") + coord_fixed()
# datSub <- subset(dat, time %in% c('D19') & genotype %in% c('LAG3-tdT+', 'CD44+tdT+'))
obj <- datSub
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)

obj <- IntegrateLayers(
  object = obj,
  method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)

obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.2, cluster.name = "clusters")

obj <- RunUMAP(
  obj,
  reduction = "integrated.mnn",
  dims = 1:15,
  reduction.name = "umap"
)

DimPlot(obj, group.by = "genotype", alpha = 0.5) +
  scale_color_manual(
    values = c("LAG3-tdT+" = "#D73027", "LAG3+tdT+" = "#800080")
  ) +
  coord_fixed() +
  theme_void()

ggsave("figures/umap_genotype.png", width = 4, height = 4)

saveRDS(obj, file = "integrated_tumorSet.rds")

datSub <- subset(obj, time == 'D11')

write.csv(Cells(datSub), file = "cellID.csv", row.names = FALSE)

# umap
Embeddings(datSub, reduction = "umap") %>%
  write.csv(file = "umap.csv")

# cluster
datSub$clust <- gsub(':.*', '', datSub$genotype)

datSub$clust %>%
  write.csv(file = "clusters.csv")

# color collection
# data.frame(id = gsub(':.*', '', levels(datSub$anno)),
#            col = hue_pal()(length(levels(datSub$anno)))) %>%
#   arrange(id) %>%
#   write.csv(file = "2.Script/veloD11D19/col.csv")

data.frame(id = c('LAG3+tdT+', 'LAG3-tdT+'), col = c("#800080", "#D73027")) %>%
  arrange(id) %>%
  write.csv(file = "col.csv")

# figure size
plt <- DimPlot(datSub) +
  coord_fixed()

# x-range
x = round(diff(layer_scales(plt)$x$range$range), 1)

# y-range
y = round(diff(layer_scales(plt)$y$range$range), 1)

write.csv(data.frame(x, y), file = "coord.csv")

# for (sub in c('overall')) {
#
#   if (sub == 'overall') {
#     datSub <- dat} else {
#       datSub <- subset(dat, sample == sub)}
#
# DimPlot(dat) +
#   coord_fixed()

# cell barcodes

# old ----
# newClst <- read.csv("1.Data/scrna/TCR_defined_clusters.csv", row.names = 1)
# names(newClst) <- "newClst"

# dat <- AddMetaData(dat, newClst)

# dat <- subset(dat, sample == "PML CD8")
# Idents(dat) <- "newClst"

# DimPlot(dat)

# write.csv(Cells(dat), file = "cellID.csv", row.names = FALSE)

# embd <- Embeddings(dat, reduction = "umap")
# write.csv(embd, file = "umap.csv")

# clst <- Idents(dat) # dat@meta.data["newClst"]
# write.csv(clst, file = "clusters.csv")

# # Idents(dat) <- factor(Idents(dat), levels = 0:10)

# hue_pal()(length(table(Idents(dat))))

# plt <- DimPlot(dat) +
#   coord_fixed()

# plt <- DimPlot(subset(dat, idents = c(0, 7))) +
#   coord_fixed()

# # x-range
# round(diff(layer_scales(plt)$x$range$range), 1)

# # y-range
# round(diff(layer_scales(plt)$y$range$range), 1)

# head(Cells(dat))

# plt <- DimPlot(subset(dat, idents = c(0, 7))) +
#   coord_fixed()
