# ----------------------------------------------------------------------------
# Fig 3 | Export embeddings and clusters for RNA velocity
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/veloD11D19/1.init.R
# Input    : 1.Data/scrna/D11D19DatWithMonocleAnnot.RDS
# Output   : 2.Script/veloD11D19/coord.csv
#            2.Script/veloD11D19/clusters.csv
# Upstream : 08_fig3_partition_analysis.R
# Note     : Spliced / unspliced counts from Velocyto
# ----------------------------------------------------------------------------

library(Seurat)
library(scales)
library(tidyverse)

dat <- readRDS("1.Data/scrna/D11D19DatWithMonocleAnnot.RDS")

datSub <- subset(dat, time %in% c('D19') & genotype %in% c('LAG3-tdT+', 'CD44+tdT+')) 

DimPlot(datSub) + coord_fixed()

write.csv(Cells(datSub), 
          file = "2.Script/veloD11D19/cellID.csv", 
          row.names = FALSE)

# umap 
Embeddings(datSub, reduction = "umap") %>% 
  write.csv(file = "2.Script/veloD11D19/umap.csv")

# cluster
datSub$clust <- gsub(':.*', '', datSub$genotype)

datSub$clust %>% 
  write.csv(file = "2.Script/veloD11D19/clusters.csv")

# color collection 
# data.frame(id = gsub(':.*', '', levels(datSub$anno)), 
#            col = hue_pal()(length(levels(datSub$anno)))) %>% 
#   arrange(id) %>% 
#   write.csv(file = "2.Script/veloD11D19/col.csv")

data.frame(id = c('CD44+tdT+', 'LAG3-tdT+'), 
           col = c("#F4A6B0", "#D73027")) %>% 
  arrange(id) %>% 
  write.csv(file = "2.Script/veloD11D19/col.csv")

# figure size
plt <- DimPlot(datSub) +
  coord_fixed()

# x-range
x = round(diff(layer_scales(plt)$x$range$range), 1)

# y-range
y = round(diff(layer_scales(plt)$y$range$range), 1)

write.csv(data.frame(x, y), 
          file =  "2.Script/veloD11D19/coord.csv")


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
