# ----------------------------------------------------------------------------
# Setup | Read CellRanger matrices and build the Seurat object
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/1.1.datInit.R
# Input    : 1.Data/scrna/ (CellRanger output)
#            1.Data/scrna/hashSettings_042722.csv
# Output   : 1.Data/scrna/datInt.RDS
# Note     : Samples resolved from TotalSeq hashtags
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)

dir.create("1.Data/scrna/initializedDat/")

reps <- list.files("1.Data/scrna/")[1:6]

for (rep in reps) {
  
  annoDf <- read.csv("1.Data/scrna/hashSettings_042722.csv")
  
  annoDf <- annoDf[annoDf$rep == rep,]
  
  hashtags <- annoDf$hashtag
  
  # gene expression matrix
  datRaw <- Read10X(paste0("1.Data/scrna/", rep, "/Gex/filtered_feature_bc_matrix/"))
  colnames(datRaw) <- gsub("-1", "", colnames(datRaw))
  
  # hashtag 
  citeRaw <- Read10X(paste0("1.Data/scrna/", rep, "/CITEseq/umi_count/"), 
                     gene.column = 1) 
  citeRaw <- citeRaw[grep("m_", rownames(citeRaw))[hashtags],]
  
  # Select cell barcodes detected by both RNA and HTO In the example datasets we have already
  # filtered the cells for you, but perform this step for clarity.
  cellsUse <- intersect(colnames(datRaw), colnames(citeRaw))
  
  # Subset RNA and HTO counts by joint cell barcodes
  datRaw <- datRaw[, cellsUse]
  citeRaw <- as.matrix(citeRaw[, cellsUse])
  
  # Confirm that the HTO have the correct names
  rownames(citeRaw) <- paste0("HTO_", hashtags)
  
  # Setup Seurat object
  datTag <- CreateSeuratObject(counts = datRaw)
  
  # Normalize RNA data with log normalization
  datTag <- NormalizeData(datTag)
  
  # Find and scale variable features
  datTag <- FindVariableFeatures(datTag, selection.method = "mean.var.plot")
  datTag <- ScaleData(datTag, features = VariableFeatures(datTag))
  
  # Add HTO data as a new assay independent from RNA
  datTag[["HTO"]] <- CreateAssayObject(counts = citeRaw)
  
  # Normalize HTO data with centered log-ratio (CLR) transformation
  datTag <- NormalizeData(datTag, assay = "HTO", normalization.method = "CLR")
  
  # If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid clustering function for large applications You can also play with additional parameters (see documentation for HTODemux()) to adjust the threshold for classification Here we are using the default settings
  
  datTag <- HTODemux(datTag, assay = "HTO", positive.quantile = 1-10E-15)
  
  table(datTag$HTO_classification.global)
  
  # Extract the singlets
  dat <- subset(datTag, HTO_classification.global == "Singlet")
  
  # add annotation ----
  dat <- AddMetaData(object = dat,
                     metadata = dat@meta.data %>% select("HTO_classification") %>% 
                       rownames_to_column(var = "cellID") %>% 
                       merge(annoDf %>% 
                               mutate(hashtag = paste0("HTO-", hashtag)), 
                             by.x = "HTO_classification",
                             by.y = "hashtag") %>% 
                       column_to_rownames(var = "cellID"))
  
  # Add tcr data ----
  tcr <- read.csv(paste0("1.Data/scrna/", rep, "/TCR/filtered_contig_annotations.csv"))
  tcr$barcode <- gsub("-1", "", tcr$barcode)
  tcr <- tcr[!duplicated(tcr$barcode), ]
  rownames(tcr) <- tcr$barcode
  
  dat <- AddMetaData(object = dat, metadata = tcr)
  
  saveRDS(dat, paste0("1.Data/scrna/initializedDat/", rep, ".RDS"))
  
  gc()
}


