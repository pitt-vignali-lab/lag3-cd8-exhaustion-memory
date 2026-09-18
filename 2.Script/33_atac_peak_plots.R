# ----------------------------------------------------------------------------
# Fig 4 | ggplot rendering of selected peak regions
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/ggplotPeaks.R
# Input    : 1.Data/atac/Merged_replicate/bigwig/
# Output   : 3.Output/atac/
# Upstream : 30_fig4_peak_visualization.R
# ----------------------------------------------------------------------------

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(rtracklayer)
library(tidyverse)


if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# -----------------------------
# 2. Define Genomic Region
# -----------------------------
# Get gene coordinates
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
gene_gr <- genes(txdb, filter = list(gene_id = gene_entrez))

if (length(gene_gr) == 0) {
  stop("Gene ID not found.")
}

annoDf <- read.table(
  file = "1.Data/atac/Merged_replicate/consensus_peaks.mRp.clN.annotatePeaks.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  dplyr::rename(PeakID = 1)

# Define query region with padding
query_gr <- GRanges(
  # seqnames = seqnames(gene_gr),
  #   ranges = IRanges(start = start(gene_gr) - left_pad,
  #                    end   = end(gene_gr) + right_pad)

  annoDf %>%
    filter(Gene.Name == 'Pdcd1') %>%
    select(Chr, Start, End)
)

# -----------------------------
# 3. Prepare File List
# -----------------------------
# A. In-house files
inhouse_dir <- "1.Data/atac/Merged_replicate/bigwig"
conds <- c("DLN_TDT_SP", "NDLN_TDT_SP", "TUMOR_DP", "TUMOR_TDT_SP")
inhouse_paths <- file.path(inhouse_dir, paste0(conds, ".mRp.clN.bigWig"))
names(inhouse_paths) <- conds

# B. Wherry Ref files
wherry_paths <- list.files(
  "1.Data/atac/RefGSE86797/merged_replicate/bigwig",
  pattern = "bigWig$",
  full.names = TRUE
)
names(wherry_paths) <- gsub(".bigWig", "", basename(wherry_paths))

# Combine all files to process
all_bw_files <- c(inhouse_paths, wherry_paths)

# -----------------------------
# 4. Extract & Normalize Data
# -----------------------------
# Helper function to process one file
process_bw <- function(filepath, condition_name, region) {
  if (!file.exists(filepath)) {
    warning(paste("File not found:", filepath))
    return(NULL)
  }

  # Import only the specific region (Fast)
  gr <- import(filepath, selection = BigWigSelection(region))

  # Normalize Score (0 to 1 based on max in this region)
  max_score <- max(gr$score, na.rm = TRUE)
  if (!is.finite(max_score) || max_score <= 0) {
    max_score <- 1
  }
  gr$score <- gr$score / max_score

  # Convert to Data Frame
  df <- as.data.frame(gr) %>%
    dplyr::select(seqnames, start, end, score) %>%
    dplyr::rename(chr = seqnames) %>%
    dplyr::mutate(
      gene = gene_symbol,
      condition = condition_name,
      source = ifelse(condition_name %in% conds, "inhouse", "wherry")
    )

  return(df)
}

# Run extraction across all files
peak_values_df <- map2_dfr(
  all_bw_files,
  names(all_bw_files),
  ~ process_bw(.x, .y, query_gr)
)

# 1. Prepare the annotation data
# Make sure column names match your loaded data (e.g. "Gene Name" vs "Gene.Name")
# 1. Prepare annotation
# Uses "Gene.Name" and "Start"/"End" as shown in your colnames

print(colnames(annoDf))

# 1. Identify which condition is on the top row
# (It is usually the first level of the factor. You can also type the string manually like "DLN_TDT_SP")
top_condition <- levels(factor(peak_values_df$condition))[1]

# 2. Add this 'condition' column to your annotation data
# This tells ggplot: "Only put these labels in the panel belonging to top_condition"
posAnno <- annoDf %>%
  filter(Gene.Name == 'Pdcd1') %>%
  mutate(
    pos = (Start + End) / 2,
    dtt = paste0(round(Distance.to.TSS / 1000, 2), 'KB \n▾'),
    posy = row_number() %% 2 == 0
  )

posAnno_top <- posAnno %>%
  mutate(condition = top_condition)

# 3. Plot
ggplot(
  peak_values_df %>% filter(condition == 'DLN_TDT_SP'),
  aes(x = start, xend = end, y = score, color = condition)
) +
  geom_segment(aes(y = 0, yend = score), linewidth = 1) +

  # Use the new 'posAnno_top' dataframe
  geom_text(
    data = posAnno_top,
    aes(x = pos, y = ifelse(posy, 1.1, 1.3), label = dtt),
    inherit.aes = FALSE, # Keep this FALSE to avoid errors
    angle = 0,
    hjust = 0.5,
    size = 3
  ) +

  facet_grid(condition ~ .) +
  ylim(0, 1.5) +
  labs(
    title = "Tox Peak Values",
    x = "Genomic Coordinates",
    y = "Peak Score"
  ) +
  theme_minimal()

ggplot(
  peak_values_df,
  aes(x = start, xend = end, y = score, color = condition)
) +
  geom_segment(aes(y = 0, yend = score), linewidth = 1)
