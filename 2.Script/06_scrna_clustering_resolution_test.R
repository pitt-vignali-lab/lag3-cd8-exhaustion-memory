# ----------------------------------------------------------------------------
# Setup | Scan clustering resolutions
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/1.6.clustringResTesting.R
# Input    : 1.Data/scrna/datIntCellCycleScore2.RDS
# Output   : 3.Output/scrna/
# Upstream : 04_scrna_cell_cycle_scoring.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(presto)
library(clustree)
library(gridExtra)

if(!dir.exists("3.Output/scrna/2.Clustering")){
  dir.create("3.Output/scrna/2.Clustering")
}

dat <- readRDS("1.Data/scrna/datIntegHarmony.RDS")

DimPlot(dat) +
  coord_fixed()

# cell cycle scoring ----
dat <- CellCycleScoring(dat, 
                        s.features = stringr::str_to_sentence(cc.genes.updated.2019$s.genes), 
                        g2m.features = stringr::str_to_sentence(cc.genes.updated.2019$g2m.genes), 
                        set.ident = TRUE)
DimPlot(dat, group.by = "Phase") + coord_fixed()

DimPlot(dat, split.by = "Phase") + coord_fixed()

# clustering ---- 
dat %>% RunHarmony("batch2", theta = 2, plot_convergence = TRUE) # , epsilon.cluster=-Inf, epsilon.harmony=-Inf)

dat <- FindNeighbors(dat, reduction = "harmony", dims = 1:35, force.recalc = T)
datTest <- FindClusters(dat, resolution = 0.8, group.singletons = T,
                        method = 1, graph.name = 'RNA_snn') 

datTest <- RunUMAP(datTest, graph = 'RNA_snn', metric = 'euclidean')

DimPlot(datTest, label = T) +
  theme_bw() +
  coord_fixed() +NoLegend()

DimPlot(datTest, split.by = "seurat_clusters") +
  theme_bw() +
  coord_fixed() + NoLegend()

DimPlot(datTest, group.by = "Phase") + coord_fixed()
DimPlot(datTest, group.by = "Phase", split.by = "Phase") + coord_fixed()
# png("3.Output/scrna/2.Clustering/cell_distribution_in_each_cluster.png", 
#     width = 6, height = 12, units = "in", res = 300)
clustree(datTest)
# dev.off()

plotDf <- data.frame(Embeddings(datTest, reduction = "umap"),
                     datTest$genotype,
                     datTest$tissue,
                     datTest$time,
                     datTest@meta.data[c(55:65)])

ggplot(plotDf, aes(UMAP_1, UMAP_2, col = RNA_snn_res.0.8)) +
  geom_point() +
  facet_grid(datTest.tissue ~ datTest.time + datTest.genotype) +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(paste0("3.Output/scrna/2.Clustering/res0.8_cluster_dist.png"),
       width = 25, units = "in", height = 8, dpi = 300, 
       bg = "white")

for (i in 7:16) {
  
  id <- names(plotDf)[i]
  wilcoxauc(datTest, group_by = id) %>%
  write.csv(paste0("3.Output/scrna/2.Clustering/", id, "_marker_gene.csv"))
  
  png(paste0("3.Output/scrna/2.Clustering/", id,"_cluster_cell_freq.png"), 
      width = 20, height = 2, units = "in", res = 300)
  dfTab <- plotDf %>% rename(cluster = i) %>% select(cluster) %>% table %>% as.data.frame() %>% t
  pTab <- tableGrob(dfTab)
  grid.arrange(pTab)
  dev.off()
  
  
  p1 <- ggplot(plotDf %>% rename(cluster = i), 
               aes(UMAP_1, UMAP_2, col = cluster)) +
    geom_point(size = 0.5)+
    geom_text(data = plotDf %>% rename(cluster = i) %>%
                group_by(cluster) %>% 
                summarise(x = mean(UMAP_1),
                          y = mean(UMAP_2)),
              aes(x, y, label = cluster), col = "black") +
    coord_fixed() +
    theme_bw() +
    theme(legend.position = "none")
  
  p2 <- ggplot(plotDf %>% rename(cluster = i), 
               aes(UMAP_1, UMAP_2, col = cluster)) +
    geom_point(size = 0.5) +
    geom_text(data = plotDf %>% rename(cluster = i) %>%
                group_by(cluster) %>% 
                summarise(x = mean(UMAP_1),
                          y = mean(UMAP_2)),
              aes(x, y, label = cluster), col = "black") +
    facet_wrap(datTest.genotype ~ datTest.time, nrow = 3) +
    coord_fixed() +
    theme_bw() +
    theme(legend.position = "none")
  
  p3 <- ggplot(plotDf %>% rename(cluster = i), 
               aes(UMAP_1, UMAP_2, col = cluster)) +
    geom_point(size = 0.5) +
    facet_graph(~ cluster)
  p1 + p2
  ggsave(paste0("3.Output/scrna/2.Clustering/", id,"_cluster_dist.png"),
         width = 15, units = "in", height = 8, dpi = 300, 
         bg = "white")
  
  cat(names(plotDf)[i], "is done \n")
}



