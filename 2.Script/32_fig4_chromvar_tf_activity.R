# ----------------------------------------------------------------------------
# Fig 4 | Transcription factor motif activity
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/chromVAR_TF.R
# Input    : 1.Data/atac/consensus_peaks.mRp.clN.featureCounts.txt
# Output   : 3.Output/atac/
# Upstream : 21_atac_deseq_init.R
# Note     : chromVAR on BSgenome.Mmusculus.UCSC.mm10 with JASPAR2020 vertebrate motifs
#            GC bias corrected
# ----------------------------------------------------------------------------

# --- 1. Load Libraries ---
library(chromVAR)
library(SummarizedExperiment)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10) # Mouse Genome
library(JASPAR2020)
library(TFBSTools)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(motifmatchr)
library(tidyverse)
library(ComplexHeatmap)
library(export)
# --- Load the mouse reference genome (mm10 assembly) ---

genome <- BSgenome.Mmusculus.UCSC.mm10


# --- Load mouse TF motifs from the JASPAR database ---
opts <- list(species = "Mus musculus", collection = "CORE")
motifs <- getMatrixSet(JASPAR2020, opts)

# Example placeholder PFM (replace with your real matrix!)
tox_pfm <- matrix(
  c(
    5,  1,  2, 12,  3,  1,
    3, 12,  1,  2,  4,  2,
    2,  1, 12,  1,  6,  3,
    10, 6,  5,  5,  7, 11
  ),
  nrow = 4, byrow = TRUE,
  dimnames = list(c("A","C","G","T"), NULL)
)

# Background frequencies must be named!
bg_vec <- c(A=0.25, C=0.25, G=0.25, T=0.25)

tox_pf <- PFMatrix(
  ID = "TOX_custom",
  name = "TOX",
  matrixClass = "Unknown",
  strand = "+",
  bg = bg_vec,
  profileMatrix = tox_pfm
)

# Add into your motif list
motifs[["TOX"]] <- tox_pf

# --- 1. Load the featureCounts output file ---
# Use `comment.char = "#"` to skip the header line that starts with #
fc_filepath <- "1.Data/atac/consensus_peaks.mRp.clN.featureCounts.txt"
fc_data <- read.table(fc_filepath, header = TRUE, sep = "\t", comment.char = "#")

peaks_gr <- GRanges(seqnames = fc_data$Chr,
                    ranges = IRanges(start = fc_data$Start, end = fc_data$End))

counts_matrix <- as.matrix(fc_data[, grepl('DLN_TDT|NDLN_TDT|TUMOR_DP|TUMOR_TDT', colnames(fc_data))])

# fc_data <- readRDS('1.Data/atac/expressionMatrixBatchEffFixed.RDS')
# an_data <- read.table('1.Data/atac/MergedReplicate/consensus_peaks.mRp.clN.annotatePeaks.txt', sep = '\t', header = TRUE) %>% 
#   dplyr::rename('peak_id' = 1) %>% 
#   column_to_rownames('peak_id')

# an_data <- an_data[rownames(fc_data), ]

# # --- 2. Extract Peak Coordinates (as a GRanges object) ---
# peaks_gr <- GRanges(seqnames = an_data$Chr,
#                     ranges = IRanges(start = an_data$Start, end = an_data$End))

# --- 3. Extract the Counts Matrix ---
# The counts start from the 7th column to the end of the file.
# counts_matrix <- as.matrix(fc_data)

# --- 4. Define Your Experimental Design (colData) ---
# IMPORTANT: You MUST customize this section to match your experiment.
# Create a data frame where each row corresponds to a sample column in your counts matrix.
sample_info <- data.frame(
  row.names = colnames(counts_matrix),
  # I am creating example groups based on your filenames. Please adjust!
  group = sub("(_REP\\d+)\\.mLb\\.clN\\.sorted\\.bam$", "", colnames(counts_matrix))
)

# --- 5. Construct the SummarizedExperiment Object ---
se <- SummarizedExperiment(
    assays = list(counts = counts_matrix),
    rowRanges = peaks_gr,
    colData = sample_info
)

# Check the object to make sure it looks correct
print(se)

# --- 1. Add GC content bias correction ---
se <- addGCBias(se, genome = genome)

# --- 2. Find motif occurrences in the peaks ---
motif_ix <- matchMotifs(motifs, se, genome = genome)

# --- 3. Compute Deviations ---
# This is the main step that calculates TF activity scores
dev <- computeDeviations(object = se, annotations = motif_ix)

# Extract the deviation scores
deviation_scores <- deviationScores(dev)

# Define the groups to compare (change these to your groups of interest)
# group1_samples <- rownames(colData(dev))[which(colData(dev)$group == "TUMOR_DP")]
# group2_samples <- rownames(colData(dev))[which(colData(dev)$group == "NDLN_DN")]

# Run a t-test for each TF to compare the two groups
# results_ttest <- t(apply(deviation_scores, 1, function(row) {
#   t_res <- t.test(row[group1_samples], row[group2_samples])
#   return(c(mean_diff = t_res$estimate[1] - t_res$estimate[2], p_value = t_res$p.value))
# }))
# results_ttest <- as.data.frame(results_ttest)
# results_ttest$p_adjust <- p.adjust(results_ttest$p_value, method = "BH")
# results_ttest$tf_name <- TFBSTools::name(motifs)[rownames(results_ttest)]

# # View top results
# head(results_ttest %>% arrange(p_adjust))

# # Visualize with a Volcano Plot
# ggplot(results_ttest %>% 
#   rename(mean_diff = 1), aes(x = mean_diff, y = -log10(p_adjust))) +
#   geom_point(aes(color = p_adjust < 0.05 & abs(mean_diff) > 1), alpha = 0.7) +
#   geom_text(aes(label=ifelse(p_adjust<0.01 & abs(mean_diff)>1.5, tf_name,'')),hjust=0,vjust=0, size=3) +
#   theme_minimal() +
#   labs(title = "TF Activity: TUMOR_DP vs NDLN_DN",
#        x = "Difference in Mean Deviation (TUMOR_DP - NDLN_DN)")

# Find the top 100 most variable TFs
variability <- computeVariability(dev, bootstrap_samples = 10000)
# top_tfs <- head(arrange(variability, desc(variability)), 100)$name

# Get the deviation scores for just these top TFs
# top_scores <- deviation_scores[top_tfs, ]
# rownames(top_scores) <- TFBSTools::name(motifs)[top_tfs] # Use TF names instead of IDs

scoreMat <- deviation_scores %>% 
  merge(variability %>% select(name), by = 0) %>% 
  select(-Row.names) %>% 
  column_to_rownames('name')

rownames(scoreMat)

cd8_relevant_tfs <- c(
  # Exhaustion
  "NFATC2", "NR4A2", "Nr2f6", "ARNT::HIF1A", "Smad2::Smad3", "Smad4", 'TOX', 
  
  # Effector
  "Stat4", "Stat5a", "Stat5b", "Stat5a::Stat5b",  
  
  # Memory
  "Foxo1", "RUNX1"
)

# Generate the heatmap
# pheatmap(scoreMat[cd8_relevant_tfs, grepl('DLN_TDT|NDLN_TDT|TUMOR_DP|TUMOR_TDT', colnames(scoreMat))],
#          main = "Top 100 Variable TF Activities",
#          show_colnames = TRUE,
#          fontsize_row = 7)


mat <- scoreMat[cd8_relevant_tfs, grepl('DLN_TDT|NDLN_TDT|TUMOR_DP|TUMOR_TDT', colnames(scoreMat))]

colnames(mat) <- gsub('.mLb.clN.sorted.bam', '', colnames(mat))

library(preprocessCore)
mat_scaled <- normalize.quantiles(as.matrix(mat))
rownames(mat_scaled) <- rownames(mat)
colnames(mat_scaled) <- colnames(mat)

# mat_scaled <- t(scale(t(mat)))
# mat_scaled <- t(apply(mat, 1, function(x) (x - min(x)) / (max(x) - min(x))))

Heatmap(mat_scaled,
    column_split = gsub('_REP[0-9]', '', colnames(mat)),
    row_split = rep(c('Exhaustion', 'Effector', 'Memory'), c(7, 4, 2)), 
    width = ncol(mat)*unit(10, "mm"), 
    height = nrow(mat)*unit(10, "mm"),
    rect_gp = gpar(col = "gray20", lwd = 1))

graph2ppt(file = '3.Output/atac/TFheatmap2.pptx', width = 8, height = 8)
