# ----------------------------------------------------------------------------
# Fig 4 | Joint PCA of tdT populations with public reference ATAC datasets
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/2.0.tdt_ref_joint_PCA.R
# Input    : 1.Data/atac/RefGSE86797/merged_library/macs2/consensus_peaks.mLb.clN.featureCounts.txt
# Output   : 3.Output/atac/
# Upstream : 21_atac_deseq_init.R
# Note     : Reference peaks merged with GenomicRanges. Public data from GSE86797
# ----------------------------------------------------------------------------

library(GenomicRanges)
library(DESeq2)
library(export)
library(tidyverse)
library(sva)

# DE analysis in our dataset ####
tdt <- read.table(
  '1.Data/atac/consensus_peaks.mLb.clN.featureCounts.txt',
  header = T
) %>%
  column_to_rownames(var = "Geneid")

ref <- read.table(
  '1.Data/atac/RefGSE86797/merged_library/macs2/consensus_peaks.mLb.clN.featureCounts.txt',
  header = T
) %>%
  column_to_rownames(var = "Geneid")

# Extract just coordinate + count columns
key_cols <- c("Chr", "Start", "End")
tdt_counts <- tdt %>% dplyr::select(all_of(key_cols), contains(".bam"))
ref_counts <- ref %>% dplyr::select(all_of(key_cols), contains(".bam"))


gr_tdt <- GRanges(
  seqnames = tdt_counts$Chr,
  ranges = IRanges(tdt_counts$Start, tdt_counts$End)
)
gr_ref <- GRanges(
  seqnames = ref_counts$Chr,
  ranges = IRanges(ref_counts$Start, ref_counts$End)
)

hits <- findOverlaps(gr_tdt, gr_ref, ignore.strand = TRUE)
# pick unique best matches per tdt peak
hits <- hits[!duplicated(queryHits(hits)) & !duplicated(subjectHits(hits))]

merged <- cbind(
  tdt_counts[queryHits(hits), ],
  ref_counts[subjectHits(hits), -(1:3)] # drop redundant Chr/Start/End
)

# Build count matrix
count_mat <- merged %>%
  dplyr::select(ends_with(".bam")) %>%
  as.matrix()

rownames(count_mat) <- paste(merged$Chr, merged$Start, merged$End, sep = "_")

coldat <- read.csv('1.Data/atac/expSet.csv')

refcol <- colnames(count_mat)[!colnames(count_mat) %in% coldat$fileName]

coldat <- rbind(
  coldat,
  data.frame(
    group = gsub('_.*', '', refcol),
    fastq_1 = '',
    replicate = gsub('.*_REP|.mLb.*', '', refcol),
    batch = 3,
    fileName = refcol
  )
)

coldat <- coldat[-grep('LAG3_SP|DN', coldat$group), ]

expdat <- count_mat[, coldat$fileName]

expBatRem <- ComBat_seq(counts = as.matrix(expdat), batch = coldat$batch)

dds <- DESeqDataSetFromMatrix(
  countData = expBatRem,
  # expdat[rownames(coldat)]  expBatRem
  colData = coldat,
  design = ~ batch + group
)

dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup = c('group'), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# smpls <- c('LAG3^-~tdT^-', 'LAG3^-~tdT^+', 'LAG3^+~tdT^+')

pcaData %>%
  mutate(
    group = factor(
      group,
      levels = c(
        'NDLN_TDT_SP',
        'DLN_TDT_SP',
        'TUMOR_DP',
        'TUMOR_TDT_SP',
        'Naive',
        'Effector',
        'Memory',
        'Exhausted'
      )
    )
  ) %>%
  ggplot(aes(PC1, PC2, fill = group)) +
  geom_point(size = 3, shape = 21, alpha = 0.6) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_fill_manual(
    values = c(
      'pink',
      'tomato',
      '#7030A0',
      '#C00000',
      "gray",
      "deepskyblue2",
      "steelblue3",
      "dodgerblue4"
    ),
    labels = c(
      bquote(bold(LAG3^-{} ~ tdT^+{} ~ NDLN)),
      bquote(bold(LAG3^-{} ~ tdT^+{} ~ DLN)),
      bquote(bold(LAG3^+{} ~ tdT^+{} ~ Tumor)),
      bquote(bold(LAG3^-{} ~ tdT^+{} ~ Tumor)),
      'Ref_Naive',
      'Ref_Effector',
      'Ref_Memory',
      'Ref_Exhausted'
    )
  ) +
  labs(fill = '') +
  # coord_fixed()  +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = 'bottom',
    text = element_text(size = 10, face = 2)
  )

graph2ppt(file = "3.Output/atac/joint_PCA_plot.ppt", width = 4, height = 4)

group_dist <- pcaData %>%
  dplyr::select(name, group, PC1, PC2) %>%
  expand_grid(
    pcaData %>%
      dplyr::select(name2 = name, group2 = group, PC1_2 = PC1, PC2_2 = PC2)
  ) %>%
  dplyr::filter(name != name2) %>%
  mutate(distance = sqrt((PC1 - PC1_2)^2 + (PC2 - PC2_2)^2)) %>%
  group_by(group, group2) %>%
  summarize(mean_distance = mean(distance), .groups = "drop") %>%
  pivot_wider(names_from = group2, values_from = mean_distance)

# 4. Test which TDT samples are closest to reference states
ref_groups <- c("Naive", "Effector", "Memory", "Exhausted")
tdt_groups <- c("NDLN_TDT_SP", "DLN_TDT_SP", "TUMOR_TDT_SP", "TUMOR_DP")

similarity_scores <- expand_grid(
  tdt = tdt_groups,
  ref = ref_groups
) %>%
  rowwise() %>%
  mutate(
    avg_distance = {
      tdt_samples <- pcaData %>% dplyr::filter(group == tdt)
      ref_samples <- pcaData %>% dplyr::filter(group == ref)

      dists <- expand_grid(
        tdt_samples %>% dplyr::select(PC1_t = PC1, PC2_t = PC2),
        ref_samples %>% dplyr::select(PC1_r = PC1, PC2_r = PC2)
      ) %>%
        mutate(dist = sqrt((PC1_t - PC1_r)^2 + (PC2_t - PC2_r)^2))

      mean(dists$dist)
    }
  ) %>%
  ungroup()

# Find closest reference for each TDT sample
closest_ref <- similarity_scores %>%
  group_by(tdt) %>%
  slice_min(avg_distance, n = 1)

group_centroids <- pcaData %>%
  group_by(group) %>%
  summarize(
    PC1_center = mean(PC1),
    PC2_center = mean(PC2),
    .groups = "drop"
  )

tdt_centroids <- group_centroids %>%
  dplyr::filter(
    group %in% c("NDLN_TDT_SP", "DLN_TDT_SP", "TUMOR_TDT_SP", "TUMOR_DP")
  )

ref_centroids <- group_centroids %>%
  dplyr::filter(group %in% c("Naive", "Effector", "Memory", "Exhausted")) %>%
  arrange(PC1_center)

tdt_arrow_data <- closest_ref %>%
  left_join(tdt_centroids, by = c("tdt" = "group")) %>%
  rename(x = PC1_center, y = PC2_center) %>%
  left_join(ref_centroids, by = c("ref" = "group")) %>%
  rename(xend = PC1_center, yend = PC2_center) %>%
  mutate(
    # Calculate midpoint for label placement
    mid_x = (x + xend) / 2,
    mid_y = (y + yend) / 2,
    # Round distance for display
    distance_label = round(avg_distance, 2)
  )

pcaData %>%
  # dplyr::filter(
  #   group %in% c("NDLN_TDT_SP", "DLN_TDT_SP", "TUMOR_TDT_SP", "TUMOR_DP")
  # ) %>%
  mutate(
    group = factor(
      group,
      levels = c(
        'NDLN_TDT_SP',
        'DLN_TDT_SP',
        'TUMOR_DP',
        'TUMOR_TDT_SP',
        'Naive',
        'Effector',
        'Memory',
        'Exhausted'
      )
    )
  ) %>%
  ggplot(aes(PC1, PC2)) +
  # Add arrows to reference states (red) ONLY
  geom_segment(
    data = tdt_arrow_data,
    linetype = "dashed",
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
    color = "red",
    linewidth = 0.8,
    alpha = 0.6
  ) +
  # Add distance labels on arrows
  geom_label(
    data = tdt_arrow_data,
    aes(x = mid_x, y = mid_y, label = distance_label),
    size = 2,
    color = "red",
    fill = "white",
    alpha = 0.8,
    fontface = "bold",
    label.padding = unit(0.15, "lines")
  ) +
  geom_point(aes(fill = group), size = 3, shape = 21, alpha = 1) +
  geom_label(
    data = group_centroids,
    aes(x = PC1_center, y = PC2_center, label = group),
    size = 2,
    alpha = 0.6,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      'pink',
      'tomato',
      '#7030A0',
      '#C00000',
      "gray",
      "deepskyblue2",
      "steelblue3",
      "dodgerblue4"
    ),
    labels = c(
      bquote(bold(LAG3^-{} ~ tdT^+{} ~ NDLN)),
      bquote(bold(LAG3^-{} ~ tdT^+{} ~ DLN)),
      bquote(bold(LAG3^+{} ~ tdT^+{} ~ Tumor)),
      bquote(bold(LAG3^-{} ~ tdT^+{} ~ Tumor)),
      'Ref_Naive',
      'Ref_Effector',
      'Ref_Memory',
      'Ref_Exhausted'
    )
  ) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    # legend.position = 'none',
    text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9)
  )

graph2ppt(file = "3.Output/atac/joint_PCA_plot.ppt", width = 3.65, height = 3.65)
