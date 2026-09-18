# ----------------------------------------------------------------------------
# Setup | The same for the intratumoural samples only
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/1.2.deInitializationTumor.R
# Input    : 1.Data/atac/consensus_peaks.mRp.clN.featureCounts.txt
# Output   : 1.Data/atac/1.deseqDdsFileTumor.RDS
# Upstream : 21_atac_deseq_init.R
# ----------------------------------------------------------------------------

library(sva)
library(DESeq2)
library(tidyverse)
library(export)

output <- read.table(
  '1.Data/atac/consensus_peaks.mLb.clN.featureCounts.txt',
  header = T
)
hskp <- read.csv('3.Output/atac/1.housekeepingPeaks.csv')
rownames(output) <- output$Geneid

output <- output %>%
  filter(!Geneid %in% hskp$x)

# get expression matrix
expdat <- output[-c(1:6)]
# saveRDS(expdat, '1.Data/atac/expressionMatrix.RDS')

coldat <- read.csv('1.Data/atac/expSet.csv')
# saveRDS(coldat, '1.Data/atac/coldat.RDS')
coldat <- coldat[-grep('LAG3_SP', coldat$group), ]
coldat <- coldat[grep('TUMOR', coldat$group), ]

expdat <- expdat[, coldat$fileName]
expBatRem <- ComBat_seq(counts = as.matrix(expdat), batch = coldat$batch)
# saveRDS(expBatRem, '1.Data/atac/expressionMatrixBatchEffFixed.RDS')

dds <- DESeqDataSetFromMatrix(
  countData = expBatRem,
  # expdat[rownames(coldat)]  expBatRem
  colData = coldat,
  design = ~ batch + group
)


# dn #00B050 sp #C00000 dp #7030A0
dds <- DESeq(dds)
saveRDS(dds, '1.Data/atac/1.deseqDdsFileTumor.RDS')

vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup = c('group'), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

smpls <- c('LAG3^-~tdT^-', 'LAG3^-~tdT^+', 'LAG3^+~tdT^+')

pcaData %>%
  mutate(
    group = gsub('TDT_|TUMOR_', '', group),
    group = factor(group, levels = c('DN', 'SP', 'DP'))
  ) %>%
  ggplot(aes(PC1, PC2, col = group)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(
    values = c('#00B050', '#C00000', '#7030A0'),
    labels = c(
      bquote(bold(LAG3^-{} ~ tdT^-{})),
      bquote(bold(LAG3^-{} ~ tdT^+{})),
      bquote(bold(LAG3^+{} ~ tdT^+{}))
    )
  ) +
  labs(col = '') +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = 'bottom', text = element_text(size = 15, face = 2))

ggsave(
  '3.Output/atac/1.pcaPlotTumor.png',
  width = 8,
  height = 4,
  units = 'in',
  bg = 'transparent'
)
