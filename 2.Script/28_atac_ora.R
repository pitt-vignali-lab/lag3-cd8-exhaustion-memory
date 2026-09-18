# ----------------------------------------------------------------------------
# Fig 4 | Over-representation analysis of differentially accessible regions
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/2.4.ora.R
# Input    : 1.Data/atac/peakAnnotationRgreat.RDS
# Output   : 1.Data/atac/fucntionalEnrichmentResult.RDS
# Upstream : 25_atac_highly_variable_peaks.R
# ----------------------------------------------------------------------------

library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
library(export)
library(GenomicRanges)
library(GenomicFeatures)
library(clusterProfiler)
library(org.Mm.eg.db)
library(export)

# orgs that are not all opening or all closed are collected 
peakBool <- data.table::fread('1.Data/atac/MergedReplicate/consensus_peaks.mRp.clN.boolean.txt')%>% 
  as.data.frame()

hskp <- read.csv('3.Output/atac/1.housekeepingPeaks.csv') 

annoDf <- data.table::fread('1.Data/atac/MergedReplicate/consensus_peaks.mRp.clN.annotatePeaks.txt') %>% 
  dplyr::rename('interval_id' = 1)

peakBool <- peakBool %>% 
  filter(!interval_id %in% hskp$x) %>% 
  merge(annoDf, by = 'interval_id')

condi <- unique(gsub('.mRp.*', '', names(peakBool)))[c(7:12, 14)]

goList <- lapply(condi, function(x){
  selecCol <- c(1:3, grep(paste0(x, '|Gene Name'), 
                          colnames(peakBool)))
  
  geneSet <- na.omit(peakBool[,selecCol]) %>% 
    pull(`Gene Name`) %>% 
    unique()
  
  go <- enrichGO(geneSet, 
                 OrgDb = org.Mm.eg.db, 
                 ont="BP", 
                 pvalueCutoff = 0.05, 
                 keyType = "SYMBOL")
  
  res <- go@result %>% 
    filter(p.adjust < 0.05) %>% 
    pull(ID)
  
  return(res)
})


names(goList) <- condi

saveRDS(goList, '3.Output/atac/2.ORA_jaccard_similarity.RDS')
# Example enrichment results
enrichment_results <- readRDS('3.Output/atac/2.ORA_jaccard_similarity.RDS')

# Convert to a binary matrix
all_pathways <- unique(unlist(enrichment_results))
binary_matrix <- sapply(enrichment_results, function(x) as.integer(all_pathways %in% x))
rownames(binary_matrix) <- all_pathways

# Calculate Jaccard Similarity
jaccard_similarity <- function(x, y) {
  intersection <- sum(x & y)
  union <- sum(x | y)
  return(intersection / union)
}

similarity_matrix <- outer(1:ncol(binary_matrix), 1:ncol(binary_matrix), Vectorize(function(i, j) {
  jaccard_similarity(binary_matrix[, i], binary_matrix[, j])
}))

# Set row and column names
colnames(similarity_matrix) <- rownames(similarity_matrix) <- names(enrichment_results)
# library(ComplexHeatmap)
# library(circlize)
# 
# col_fun = colorRamp2(c(0.5, 0.75, 1), c("white", "yellow", "red"))
# 
# Heatmap(
#   similarity_matrix, 
#   col = col_fun,
#   rect_gp = gpar(col = "white", lwd = 2),
#   row_dend_reorder = TRUE
# )

# Convert to a data frame for visualization
similarity_matrix[lower.tri(similarity_matrix, diag = T)] <- NA
similarity_df <- reshape2::melt(similarity_matrix) %>% 
  # mutate(Var1 = factor(Var1, levels = c("NDLN_DN", "NDLN_TDT_SP", 
  #                                       "DLN_DN", "DLN_TDT_SP",  
  #                                       "TUMOR_DN", "TUMOR_DP", "TUMOR_TDT_SP")),
  #        Var2 = factor(Var2, levels = c("NDLN_DN", "NDLN_TDT_SP", 
  #                                       "DLN_DN", "DLN_TDT_SP",  
  #                                       "TUMOR_DN", "TUMOR_DP", "TUMOR_TDT_SP"))) %>% 
  na.omit()

# Plot the similarity matrix
jsMp <- ggplot(similarity_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), angle = 45, size = 8/2.8) +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "Jaccard Similarity Matrix", x = "Condition", y = "Condition", fill = "Similarity") +
  theme_void() +
  scale_x_discrete(position = "top") +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0, size = 8, face = 2),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 8, face = 2)) +
  coord_fixed()

graph2ppt(jsMp, '3.Output/atac/2.pseaSimilarity.pptx', height = 4, width = 4)

# Save the similarity matrix
write.csv(similarity_matrix, file = "similarity_matrix.csv")

selecCol <- c(1:3, grep(paste0(condi[1], '|Gene Name'), colnames(peakBool)))
peakBoolSub1 <- na.omit(peakBool[,selecCol]) %>% 
  pull(`Gene Name`) %>% 
  unique()


go <- enrichGO(peakBoolSub1, OrgDb = org.Mm.eg.db, ont="BP", pvalueCutoff = 0.05, keyType = "SYMBOL")
goSub <- gofilter(go, level = 2)


selecCol <- c(1:3, grep('TUMOR_TDT_SP', colnames(peakBool)))
peakBoolSub2 <- na.omit(peakBool[,selecCol])

# Convert to GRanges objects
peaks_gr <- GRanges(seqnames = peakBoolSub1$chr, ranges = IRanges(start = peakBoolSub1$start, end = peakBoolSub1$end))
reference_gr <- GRanges(seqnames = peakBoolSub2$chr, ranges = IRanges(start = peakBoolSub2$start, end = peakBoolSub2$end))

# Find overlaps between peaks and reference annotations
overlaps <- findOverlaps(peaks_gr, reference_gr)

# Calculate observed overlap count
observed_overlap <- length(overlaps)


# Define a function to shuffle peak locations
shuffle_peaks <- function(peaks_gr, genome_size) {
  # Generate random starts
  random_starts <- sample(1:genome_size, length(peaks_gr), replace = TRUE)
  # Generate shuffled peaks
  shuffled_peaks <- GRanges(seqnames = seqnames(peaks_gr), ranges = IRanges(start = random_starts, width = width(peaks_gr)))
  return(shuffled_peaks)
}

# Get the total genome size (sum of all chromosomes lengths)
# Adjust this based on the actual genome you are working with
genome_size <- length(peaks_gr)

# Perform permutation test to assess significance
n_permutations <- 1000

permuted_overlaps <- replicate(n_permutations, {
  permuted_peaks <- shuffle_peaks(peaks_gr, genome_size)
  length(findOverlaps(permuted_peaks, reference_gr))
})

# Calculate p-value
p_value <- sum(permuted_overlaps >= observed_overlap) / n_permutations

# Print results
cat("Observed overlap:", observed_overlap, "\n")
cat("P-value:", p_value, "\n")

# Visualization (optional)
hist(permuted_overlaps, main = "Permutation Test Results", xlab = "Overlap Count")
abline(v = observed_overlap, col = "red", lwd = 2)
