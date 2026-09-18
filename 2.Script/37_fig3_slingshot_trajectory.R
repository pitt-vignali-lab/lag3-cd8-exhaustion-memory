# ----------------------------------------------------------------------------
# Fig 3 | Slingshot developmental trajectories across tdT subsets
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/slingshot.RE.R
# Input    : 1.Data/scrna/D11D19DatWithMonocleAnnot.RDS
# Output   : 3.Output/scrna/
# Upstream : 08_fig3_partition_analysis.R
# Note     : Trajectories predicted across ranked cell groups in tumour and periphery
# ----------------------------------------------------------------------------

library(slingshot)
library(SingleCellExperiment)
library(uwot)
library(mclust, quietly = TRUE)
library(RColorBrewer)
library(grDevices)

# means <- rbind(
#   # non-DE genes
#   matrix(rep(rep(c(0.1,0.5,1,2,3), each = 300),100),
#          ncol = 300, byrow = TRUE),
#   # early deactivation
#   matrix(rep(exp(atan( ((300:1)-200)/50 )),50), ncol = 300, byrow = TRUE),
#   # late deactivation
#   matrix(rep(exp(atan( ((300:1)-100)/50 )),50), ncol = 300, byrow = TRUE),
#   # early activation
#   matrix(rep(exp(atan( ((1:300)-100)/50 )),50), ncol = 300, byrow = TRUE),
#   # late activation
#   matrix(rep(exp(atan( ((1:300)-200)/50 )),50), ncol = 300, byrow = TRUE),
#   # transient
#   matrix(rep(exp(atan( c((1:100)/33, rep(3,100), (100:1)/33) )),50), 
#          ncol = 300, byrow = TRUE)
# )
# 
# counts <- apply(means,2,function(cell_means){
#   total <- rnbinom(1, mu = 7500, size = 4)
#   rmultinom(1, total, cell_means)
# })
# 
# rownames(counts) <- paste0('G',1:750)
# colnames(counts) <- paste0('c',1:300)

expBatRem <- readRDS( '1.Data/atac/expressionMatrixBatchEffFixed.RDS')
colnames(expBatRem) <- gsub('\\_REP[0-9].mLb.clN.sorted.bam', '', colnames(expBatRem))

sce <- SingleCellExperiment(assays = List(counts = expBatRem))


FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}

assays(sce)$norm <- FQnorm(assays(sce)$counts)

pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')

plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)

cl1 <- Mclust(rd1)$classification
colData(sce)$GMM <- cl1


plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

cl2 <- kmeans(rd1, centers = 4)$cluster
colData(sce)$kmeans <- cl2

plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

sce <- slingshot(sce)

summary(sce$slingPseudotime_1)


colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

library(tradeSeq)

# fit negative binomial GAM
sce <- fitGAM(sce)

# test for dynamic expression
ATres <- associationTest(sce)

topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce$GMM[pst.ord]

heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

lin1 <- getLineages(rd, cl, start.clus = '1')
lin1