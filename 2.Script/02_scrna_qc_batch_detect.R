# ----------------------------------------------------------------------------
# Setup | Quality control and batch effect assessment
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/1.2.QC_batch_effects_detect.R
# Input    : 1.Data/scrna/datInt.RDS
# Output   : 3.Output/scrna/
# Upstream : 01_scrna_data_init.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)

if(!dir.exists("3.Output/scrna/1.DatInit")){
  dir.create("3.Output/scrna/1.DatInit")
}

# merge samples ----
files <- list.files("1.Data/scrna/initializedDat/")

for (file in files) {
  assign(gsub(".RDS", "", file), readRDS(paste0("1.Data/scrna/initializedDat/", file)))
}

rm(list = ls()[- grep("GL_",ls())])

projs <- ls()

dat <- merge(eval(parse(text = projs[1])), 
                  y=sapply(projs[-1], function(x) eval(parse(text = x))), 
                  add.cell.ids = projs,
                  project = "ALLCELLS")

rm(list = ls()[ls() != "dat"])

# QC ----

dat$log10GenesPerUMI <- log10(dat$nFeature_RNA) / log10(dat$nCount_RNA)
dat$mitoRatio <- PercentageFeatureSet(object = dat, pattern = "^mt-")/100
dat$ribonRatio <- PercentageFeatureSet(dat, pattern = "^Rp[sl]")/100
dat$hbRatio <- PercentageFeatureSet(dat, pattern = "^Hb[^(p)]")/100

# Visualize the number of cell counts per sample
dat@meta.data %>% 
  ggplot(aes(x=rep, fill=rep)) + 
  geom_bar()+
  geom_label(stat='count', aes(label = ..count..), 
             colour = "black", fill = "White", size = 3, alpha = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        legend.position = "none") +
  ggtitle("NCells by replicates")
ggsave("3.Output/scrna/1.DatInit/cell_number_distr_by_replicates.png", 
       width = 8, height = 5, bg = "white")


dat@meta.data %>% 
  ggplot(aes(x=genotype, fill=genotype)) + 
  geom_bar() +
  geom_label(stat='count', aes(label = ..count..), 
             colour = "black", fill = "White", size = 3, alpha = 0.8) +
  theme_classic() +
  ggtitle("NCells") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        legend.position = "none") 
ggsave("3.Output/scrna/1.DatInit/cell_number_distr_by_genotype.png", 
       width = 8, height = 5, bg = "white")

# Visualize the number UMIs/transcripts per cell
dat@meta.data %>% 
  ggplot(aes(x = nCount_RNA, color=rep, fill= rep)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 300) + 
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        legend.position = c(0.9, 0.8)) +
  ggtitle("nCount_RNA")
ggsave("3.Output/scrna/1.DatInit/cell_UMIs_distr_by_replicates.png", 
       width = 8, height = 5, bg = "white")


# Visualize the distribution of genes detected per cell
dat@meta.data %>% 
  ggplot(aes(color=rep, x=nFeature_RNA, fill= rep)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 250) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nFeature_RNA")
ggsave("3.Output/scrna/1.DatInit/genes_per_cell_distr.png", 
       width = 8, height = 5, bg = "white")

dat@meta.data %>% 
  ggplot(aes(y=log10(nFeature_RNA), x= rep, fill= rep)) + 
  geom_boxplot(alpha = 0.2) + 
  # scale_x_log10() + 
  theme_classic() +
  ylab("lgo10(nGene)") +
  geom_vline(xintercept = 300) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nFeature_RNA")
ggsave("3.Output/scrna/1.DatInit/genes_per_cell_box.png", 
       width = 8, height = 5, bg = "white")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs 
dat@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 800) +
  geom_hline(yintercept = 500) +
  facet_wrap(~rep) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCount_RNA vs. nFeature_RNA")
ggsave("3.Output/scrna/1.DatInit/UMIs_vs_genes detected.png", 
       width = 8, height = 5, bg = "white")


# Visualize the distribution of mitochondrial gene expression detected per cell
dat@meta.data %>% 
  ggplot(aes(color=rep, x=mitoRatio, fill=rep)) + 
  geom_density(alpha = 0.15) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.5) +
  theme(plot.title = element_text(hjust=0.2, face="bold")) +
  ggtitle("Mitochondrial counts ratio")
ggsave("3.Output/scrna/1.DatInit/mitochondrialCountsRatio.png", width = 8, height = 5, bg = "white")

dat@meta.data %>% 
  ggplot(aes(color=rep, x=ribonRatio, fill=rep)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.15) +
  theme(plot.title = element_text(hjust=0.2, face="bold")) +
  ggtitle("Ribosome counts ratio")
ggsave("3.Output/scrna/1.DatInit/ribosomeCountsRatio.png", 
       width = 8, height = 5, bg = "white")


dat@meta.data %>% 
  ggplot(aes(color=rep, x=hbRatio, fill=rep)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 2e-3) +
  theme(plot.title = element_text(hjust=0.2, face="bold")) +
  ggtitle("Hemoglobin counts ratio")
ggsave("3.Output/scrna/1.DatInit/hemoglobinCountsRatio.png", 
       width = 8, height = 5, bg = "white")

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
dat@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = rep, fill=rep)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = c(0.8)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Complexity")
ggsave("3.Output/scrna/1.DatInit/complexity.png", 
       width = 8, height = 5, bg = "white")

# Filter out low quality reads using selected thresholds - these will change with experiment
datFiltered <- subset(x = dat, 
                      subset= nCount_RNA >= 300 & 
                        nFeature_RNA >= 250 &
                        log10GenesPerUMI > 0.80 &
                        mitoRatio < 0.2)

# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = datFiltered, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keepGenes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filteredCts <- counts[keepGenes, ]

# Reassign to filtered Seurat object
datFiltered <- datFiltered[, colnames(datFiltered) %in% colnames(filteredCts)]

rm(counts, nonzero, keepGenes, filteredCts)

datFiltered@meta.data %>% 
  ggplot(aes(x=genotype, fill = rep)) + 
  geom_bar() +
  geom_label(stat='count', aes(label = ..count..), 
             colour = "black", fill = "White", size = 3, alpha = 0.8) +
  theme_classic() +
  ggtitle("NCells") +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold")) 

ggsave("3.Output/scrna/1.DatInit/cell_number_distr_split_after_filtering.png", 
       width = 8, height = 8, bg = "white")

# batch effects checking

datFiltered <- NormalizeData(datFiltered)

# memory.limit(9999999999)
datFiltered <- ScaleData(datFiltered, features = rownames(datFiltered))

datFiltered <- FindVariableFeatures(datFiltered)

datFiltered <- RunPCA(datFiltered, npcs = 50)

ElbowPlot(datFiltered, ndims = 50)
ggsave("3.Output/scrna/1.DatInit/elbowPlotBeforeBatEffRm.png", width = 15, height = 4, bg = "white")

# Determine percent of variation associated with each PC
pct <- datFiltered[["pca"]]@stdev / sum(datFiltered[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# Minimum of the two calculation
pcs <- min(co1, co2)
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))
# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
ggsave("3.Output/scrna/1.DatInit/elbowPCSelectionBeforeBatEffRm.png", width = 15, height = 4, bg = "white")

# datFiltered <- RunTSNE(datFiltered, dims = 1:co2)

datFiltered <- RunUMAP(datFiltered, dims = 1:co2)

datFiltered <- RunTSNE(datFiltered, dims = 1:co2)

cowplot::plot_grid(ncol = 3,
                   DimPlot(datFiltered, reduction = "pca", group.by = "rep")+
                     NoAxes()+ggtitle("PCA") + coord_fixed(ratio = 1),
                   DimPlot(datFiltered, reduction = "tsne", group.by = "rep")+
                     NoAxes()+ggtitle("TSNE") + coord_fixed(ratio = 1),
                   DimPlot(datFiltered, reduction = "umap", group.by = "rep")+
                     NoAxes()+ggtitle("UMAP") + coord_fixed(ratio = 1)
)

ggsave("3.Output/scrna/1.DatInit/batchEffChecking_by_rep.png", width = 16, height = 5)

cowplot::plot_grid(ncol = 3,
                   DimPlot(datFiltered, reduction = "pca", group.by = "batch")+
                     NoAxes()+ggtitle("PCA") + coord_fixed(ratio = 1),
                   DimPlot(datFiltered, reduction = "tsne", group.by = "batch")+
                     NoAxes()+ggtitle("TSNE") + coord_fixed(ratio = 1),
                   DimPlot(datFiltered, reduction = "umap", group.by = "batch")+
                     NoAxes()+ggtitle("UMAP") + coord_fixed(ratio = 1)
)

ggsave("3.Output/scrna/1.DatInit/batchEffChecking_by_batch.png", width = 16, height = 5)

#  information supplement and fixation 
datFiltered$genotype <- gsub(" ", "", datFiltered$genotype)

datFiltered$genotype[datFiltered$genotype == "CD44-tdT-CD62L+"] <- "CD44-CD62L+tdT-"

datFiltered$batch2 <- ifelse(datFiltered$batch %in% c("Exp. 083121" ,"Exp. 091521"),
                     "batch_1", "batch_2")
# figures 
cowplot::plot_grid(ncol = 3,
                   DimPlot(datFiltered, reduction = "pca", group.by = "batch2")+
                     NoAxes()+ggtitle("PCA") + coord_fixed(ratio = 1),
                   DimPlot(datFiltered, reduction = "tsne", group.by = "batch2")+
                     NoAxes()+ggtitle("TSNE") + coord_fixed(ratio = 1),
                   DimPlot(datFiltered, reduction = "umap", group.by = "batch2")+
                     NoAxes()+ggtitle("UMAP") + coord_fixed(ratio = 1)
)

ggsave("3.Output/scrna/1.DatInit/batchEffChecking_by_batch2.png", width = 16, height = 5)

saveRDS(datFiltered, "1.Data/scrna/datInt.RDS")

unlink("1.Data/scrna/initializedDat/", recursive = TRUE)
