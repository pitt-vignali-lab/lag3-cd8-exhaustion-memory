# ----------------------------------------------------------------------------
# Setup | Cell cycle scoring and regression
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/1.4.cell_cycle_scoring.R
# Input    : 1.Data/scrna/datIntegRpca.RDS
# Output   : 1.Data/scrna/datIntCellCycleScore1.RDS
#            1.Data/scrna/datIntCellCycleScore2.RDS
# Upstream : 03_scrna_batch_correction.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)

dat <- readRDS("1.Data/scrna/datInt.RDS")
dat <- CellCycleScoring(dat, 
                        s.features = stringr::str_to_sentence(cc.genes.updated.2019$s.genes), 
                        g2m.features = stringr::str_to_sentence(cc.genes.updated.2019$g2m.genes), 
                        set.ident = TRUE)
dat <- RunPCA(dat, 
              features = c(stringr::str_to_sentence(cc.genes.updated.2019$s.genes),
                           stringr::str_to_sentence(cc.genes.updated.2019$g2m.genes)))

DimPlot(dat, reduction = "pca") + coord_fixed()
ggsave("3.Output/scrna/1.DatInit/cellcyclescroingBeforeOverall.png", 
       width = 5, height = 5, units = "in", dpi = 300)

DimPlot(dat, reduction = "pca", split.by = "Phase") + coord_fixed()
ggsave("3.Output/scrna/1.DatInit/cellcyclescroingBeforeSplit.png", 
       width = 12, height = 4, units = "in", dpi = 300)


dat <- ScaleData(dat, vars.to.regress = c("S.Score", "G2M.Score"), 
                 features = rownames(dat))
saveRDS(dat, "1.Data/scrna/datIntCellCycleScore1.RDS")

dat$CC.Difference <- dat$S.Score - dat$G2M.Score
dat <- ScaleData(dat, vars.to.regress = "CC.Difference", features = rownames(dat))
saveRDS(dat, "1.Data/scrna/datIntCellCycleScore2.RDS")

dat <- readRDS("1.Data/scrna/datIntCellCycleScore2.RDS")
dat <- RunPCA(dat, 
              features = c(stringr::str_to_sentence(cc.genes.updated.2019$s.genes),
                           stringr::str_to_sentence(cc.genes.updated.2019$g2m.genes)))