# ----------------------------------------------------------------------------
# Fig 4 | Call differentially accessible peaks for track visualisation
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/2.0.DEpeaksCallingForPeakVisz.R
# Input    : 1.Data/atac/expressionMatrixBatchEffFixed.RDS
# Output   : 3.Output/atac/
# Upstream : 21_atac_deseq_init.R
# Note     : DESeq2 on the consensus peak matrix
# ----------------------------------------------------------------------------

library(tidyverse)
library(rtracklayer)
library(GenomicRanges)

# Colors for your 4 CD8 conditions
wherryDat <- list.files(
  "1.Data/atac/RefGSE86797/merged_replicate/bigwig",
  pattern = "bigWig$",
  full.names = TRUE
)

boolDf <- data.table::fread(
  "1.Data/atac/Merged_replicate/consensus_peaks.mRp.clN.boolean.txt"
)
anno <- data.table::fread(
  "1.Data/atac/Merged_replicate/consensus_peaks.mRp.clN.annotatePeaks.txt"
) %>%
  dplyr::rename('peak_id' = 1)

myDat <- list.files(
  "1.Data/atac/Merged_replicate/bigwig",
  pattern = "bigWig$",
  full.names = TRUE
)
conds <- c("DLN_TDT_SP", "NDLN_TDT_SP", "TUMOR_DP", "TUMOR_TDT_SP")

geneRange <- anno %>% filter(`Gene Name` == "Tox")
boolRange <- boolDf %>% filter(interval_id %in% geneRange$peak_id)


# Method 2: Extract for all peaks in boolRange
all_regions <- GRanges(
  seqnames = boolRange$chr,
  ranges = IRanges(start = boolRange$start, end = boolRange$end)
)


# Import signal for all regions
bw_signal_all <- import(
  "1.Data/atac/Merged_replicate/bigwig/DLN_TDT_SP.mRp.clN.bigWig",
  which = all_regions,
  as = "GRanges"
)

# Convert GRanges to data frame for ggplot
bw_df <- as.data.frame(bw_signal_all) %>%
  mutate(
    position = start,
    signal = score
  )
