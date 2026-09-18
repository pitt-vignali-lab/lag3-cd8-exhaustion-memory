# ----------------------------------------------------------------------------
# Fig 4 | Accessibility at the Tox locus
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/peakvisTox.R
# Input    : 1.Data/atac/Merged_replicate/bigwig/
# Output   : 3.Output/atac/
# Upstream : 30_fig4_peak_visualization.R
# ----------------------------------------------------------------------------

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(Gviz)
library(GenomicFeatures)
library(rtracklayer)
library(GenomicRanges)
library(karyoploteR)
library(tidyverse)
library(data.table)


# -----------------------------
# Config
# -----------------------------

# Colors for your 4 CD8 conditions
colPan <- c("pink", "pink", "red", "purple")
names(colPan) <- c("NDLN_TDT_SP", "DLN_TDT_SP", "TUMOR_TDT_SP", "TUMOR_DP")

# Wherry ref bigWigs
wherryDat <- list.files(
  "1.Data/atac/RefGSE86797/",
  pattern = "bigWig$",
  full.names = TRUE
)
# keep their visual order if you like; example from your code:
ord <- c(4, 3, 1, 2)
ord <- ord[ord <= length(wherryDat)]
if (length(ord)) {
  wherryDat <- wherryDat[ord]
}

wherryCol <- c("gray", "dodgerblue1", "dodgerblue2", "dodgerblue4")
# will remap to basenames later

# Boolean peak file (for highlights)
bool_path <- "1.Data/atac/MergedReplicate/consensus_peaks.mRp.clN.boolean.txt"
stopifnot(file.exists(bool_path))
boolDf <- data.table::fread(bool_path)

# Output dir
out_dir <- "3.Output"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Gene(s) to plot (symbol -> Entrez ID)
geneIdDf <- data.frame(
  symbols = c("Tox"),
  id = c("252838"),
  stringsAsFactors = FALSE
)

# Region padding around gene (if annotations missing)
left_pad <- 40000
right_pad <- 50000

# Your 4 in-house bigWigs (must exist under this folder)
inhouse_dir <- "1.Data/atac/MergedReplicate/bigwig"
conds <- c("DLN_TDT_SP", "NDLN_TDT_SP", "TUMOR_DP", "TUMOR_TDT_SP")

# -----------------------------
# Helpers
# -----------------------------

safe_import_bw <- function(path, chr, from, to) {
  import(path, as = "GRanges", which = GRanges(chr, IRanges(from, to)))
}

normalize_gr_list <- function(gr_list) {
  allScore <- unlist(lapply(gr_list, function(gr) gr$score))
  maxScore <- suppressWarnings(max(allScore, na.rm = TRUE))
  if (!is.finite(maxScore) || maxScore <= 0) {
    maxScore <- 1
  }
  lapply(gr_list, function(gr) {
    gr$score <- gr$score / maxScore
    gr
  })
}

# Make a tidy highlight data.frame for a given condition column name in boolDfSub
make_highlight_df <- function(boolDfSub, cond_name) {
  if (!cond_name %in% names(boolDfSub)) {
    return(NULL)
  }
  t <- boolDfSub %>%
    dplyr::select(chr, start, end, all_of(cond_name)) %>%
    dplyr::rename(value = 4) %>%
    dplyr::filter(!is.na(value), value > 0)
  if (nrow(t) == 0) {
    return(NULL)
  }
  t
}

# Convert any chromosome formats in boolDf to "chrN" to match Gviz
normalize_chr_col <- function(df) {
  if (!"chr" %in% names(df)) {
    return(df)
  }
  chr <- df$chr
  # try to coerce to character "chrN"
  chr_char <- as.character(chr)
  # if already "chr", keep; else prefix "chr"
  chr_char <- ifelse(grepl("^chr", chr_char), chr_char, paste0("chr", chr_char))
  df$chr <- chr_char
  df
}

# -----------------------------
# Tracks shared across genes
# -----------------------------

gtrack <- GenomeAxisTrack(
  lwd = 0.4,
  col = "black",
  fontcolor = "black",
  fontsize = 10
)
displayPars(gtrack) <- list(size = 2)

GRList <- transcriptsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = "gene")

# -----------------------------
# Main loop over genes
# -----------------------------

for (gene in geneIdDf$symbols) {
  message(">>> Processing gene: ", gene)

  # Extract gene model and locus
  ent_id <- geneIdDf %>%
    dplyr::filter(symbols == gene) %>%
    dplyr::pull(id) %>%
    as.character()
  stopifnot(length(ent_id) == 1)

  if (!ent_id %in% names(GRList)) {
    stop("Entrez ID ", ent_id, " not found in TxDb.")
  }

  infoDf <- data.frame(GRList[[ent_id]])
  chr <- as.character(unique(infoDf$seqnames))[1]
  start <- min(c(infoDf$start, infoDf$end))
  end <- max(c(infoDf$start, infoDf$end))

  # Reference gene track
  refG <- GeneRegionTrack(
    TxDb.Mmusculus.UCSC.mm10.knownGene,
    chromosome = chr,
    start = start,
    end = end,
    geneSymbol = TRUE,
    transcriptAnnotation = "gene",
    name = "Gene",
    col = "black",
    fill = "black",
    fontcolor.group = "black",
    fontsize.group = 17
  )

  # Rewrite displayed ranges to carry the SYMBOL instead of numeric id
  refGDf <- ranges(refG) %>%
    data.frame() %>%
    dplyr::filter(gene == as.numeric(ent_id))
  z <- GRanges(
    refGDf$seqnames,
    IRanges(refGDf$start, refGDf$end),
    strand = refGDf$strand,
    feature = refGDf$feature,
    id = gene,
    exon = refGDf$exon,
    gene = gene,
    symbol = gene,
    transcript = gene
  )
  ranges(refG) <- z
  displayPars(refG) <- list(size = 0.8)

  # Determine plotting window using your annotatePeaks file (if present)
  anno_path <- "1.Data/atac/consensus_peaks.mLb.clN.annotatePeaks.txt"
  if (file.exists(anno_path)) {
    peakAnno <- data.table::fread(anno_path) %>%
      dplyr::rename(Geneid = 1) %>%
      dplyr::filter(`Gene Name` == gene)
  } else {
    peakAnno <- data.frame()
  }

  if (nrow(peakAnno) > 0) {
    gRange <- range(c(peakAnno$Start, peakAnno$End))
    startNew <- min(start - left_pad, gRange[1])
    endNew <- max(end + right_pad, gRange[2])
  } else {
    startNew <- start - left_pad
    endNew <- end + right_pad
  }

  # -------------------------
  # In-house bigWigs -> DataTracks (normalized jointly)
  # -------------------------
  # import region for each condition
  bw_in_list <- list()
  for (i in conds) {
    bw_path <- file.path(inhouse_dir, paste0(i, ".mRp.clN.bigWig"))
    stopifnot(file.exists(bw_path))
    bw_in_list[[i]] <- safe_import_bw(bw_path, chr, startNew, endNew)
  }
  # normalize scores together
  bw_in_list <- normalize_gr_list(bw_in_list)

  # Build DataTracks (histogram)
  peak_tracks <- list()
  for (i in conds) {
    gr <- bw_in_list[[i]]
    nm <- gsub("_TDT", "", i)
    peak_tracks[[nm]] <- DataTrack(
      gr,
      name = gsub("_", "\n", nm),
      type = "histogram",
      fill.histogram = colPan[i],
      col.histogram = colPan[i],
      ylim = c(0, 1),
      aggregation = "mean"
    )
    displayPars(peak_tracks[[nm]]) <- list(
      fontcolor.title = "black",
      fontsize = 10,
      background.title = "white",
      col.axis = "black",
      lwd = 0.4,
      size = 0.8
    )
  }

  # -------------------------
  # Wherry bigWigs -> DataTracks (normalized jointly)
  # -------------------------
  bw_wherry_list <- list()
  for (bw_path in wherryDat) {
    if (!file.exists(bw_path)) {
      warning("Wherry bigWig not found: ", bw_path)
      next
    }
    bw_wherry_list[[basename(bw_path)]] <- safe_import_bw(
      bw_path,
      chr,
      startNew,
      endNew
    )
  }
  bw_wherry_list <- normalize_gr_list(bw_wherry_list)

  wherry_tracks <- list()
  for (i in seq_along(wherryDat)) {
    bw_path <- wherryDat[i]
    bn <- basename(bw_path)
    if (!bn %in% names(bw_wherry_list)) {
      next
    }

    gr <- bw_wherry_list[[bn]]
    wherry_tracks[[bn]] <- DataTrack(
      gr,
      name = gsub(".bigWig", "", bn),
      type = "histogram",
      fill.histogram = wherryCol[i],
      col.histogram = wherryCol[i],
      ylim = c(0, 1),
      aggregation = "mean"
    )
    displayPars(wherry_tracks[[bn]]) <- list(
      fontcolor.title = "black",
      fontsize = 10,
      background.title = "white",
      col.axis = "black",
      lwd = 0.4,
      size = 0.8
    )
  }

  # -------------------------
  # Boolean peaks -> AnnotationTracks for highlights
  # -------------------------
  boolDfSub <- boolDf %>%
    dplyr::filter(chr == !!chr, end >= startNew, start <= endNew)

  boolDfSub <- normalize_chr_col(boolDfSub)

  highlight_tracks <- list()
  for (cond in conds) {
    hl_df <- make_highlight_df(boolDfSub, cond)
    if (!is.null(hl_df) && nrow(hl_df) > 0) {
      nm <- gsub("_TDT", "", cond)
      highlight_tracks[[nm]] <- AnnotationTrack(
        GRanges(hl_df$chr, IRanges(hl_df$start, hl_df$end)),
        name = paste0(nm, "\npeak"),
        fill = colPan[cond],
        col = colPan[cond],
        stacking = "dense"
      )
      displayPars(highlight_tracks[[nm]]) <- list(
        fontcolor.title = "black",
        fontsize = 10,
        background.title = "white",
        size = 0.5
      )
    }
  }

  # -------------------------
  # Plot
  # -------------------------
  pdf_path <- file.path(out_dir, paste0(gene, "_tracks.pdf"))
  pdf(pdf_path, width = 12, height = 10)

  track_list <- c(
    list(gtrack, refG),
    peak_tracks,
    highlight_tracks,
    wherry_tracks
  )

  plotTracks(
    track_list,
    chromosome = chr,
    from = startNew,
    to = endNew,
    transcriptAnnotation = "symbol"
  )

  dev.off()
  message("Saved: ", pdf_path)

  # -------------------------
  # Extract peak values to data frame
  # -------------------------
  peak_values_list <- list()

  # In-house peaks
  for (cond in conds) {
    gr <- bw_in_list[[cond]]
    if (length(gr) > 0) {
      peak_values_list[[cond]] <- data.frame(
        gene = gene,
        condition = cond,
        chr = as.character(seqnames(gr)),
        start = start(gr),
        end = end(gr),
        score = gr$score,
        source = "inhouse",
        stringsAsFactors = FALSE
      )
    }
  }

  # Wherry peaks
  for (bw_name in names(bw_wherry_list)) {
    gr <- bw_wherry_list[[bw_name]]
    if (length(gr) > 0) {
      peak_values_list[[bw_name]] <- data.frame(
        gene = gene,
        condition = gsub(".bigWig|.mRp.clN", "", bw_name),
        chr = as.character(seqnames(gr)),
        start = start(gr),
        end = end(gr),
        score = gr$score,
        source = "wherry",
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine into single data frame for this gene
  peak_values_df <- bind_rows(peak_values_list)

  # Save to file
  peak_csv <- file.path(out_dir, paste0(gene, "_peak_values.csv"))
  write.csv(peak_values_df, peak_csv, row.names = FALSE)
  message("Saved peak values: ", peak_csv)
}

df <- read.csv('3.Output/atac/Tox_peak_values.csv') %>%
  mutate(condition = factor(condition, levels = unique(condition)))

ggplot(df, aes(x = start, xend = end, y = score, color = condition)) +
  geom_segment(aes(y = 0, yend = score), size = 1) +
  facet_grid(condition ~ ., scales = "free_y") +
  ylim(0, 1) +
  labs(
    title = "Toxicity Peak Values",
    x = "Genomic Coordinates",
    y = "Peak Score"
  ) +
  theme_minimal()

# -------------------------
# Compare similarity among conditions
# -------------------------

# Create a matrix of peak scores across genomic bins
# Use a binning approach to align peaks across conditions

# Define bins across the region
bin_size <- 100 # adjust as needed
bins <- seq(startNew, endNew, by = bin_size)

# Create matrix with conditions as columns, bins as rows
score_matrix <- matrix(
  0,
  nrow = length(bins) - 1,
  ncol = length(peak_values_list)
)
colnames(score_matrix) <- names(peak_values_list)

# Fill matrix with scores
for (cond_name in names(peak_values_list)) {
  df_cond <- peak_values_list[[cond_name]]

  for (i in seq_len(nrow(df_cond))) {
    # Find which bins this peak overlaps
    peak_start <- df_cond$start[i]
    peak_end <- df_cond$end[i]
    peak_score <- df_cond$score[i]

    bin_idx <- which(bins[-length(bins)] <= peak_end & bins[-1] >= peak_start)
    score_matrix[bin_idx, cond_name] <- pmax(
      score_matrix[bin_idx, cond_name],
      peak_score
    )
  }
}

# Calculate correlation matrix
cor_matrix <- cor(score_matrix, method = "pearson")

# Save correlation matrix
cor_csv <- file.path(out_dir, paste0(gene, "_correlation_matrix.csv"))
write.csv(cor_matrix, cor_csv)
message("Saved correlation matrix: ", cor_csv)

# Create heatmap
library(pheatmap)

pdf(
  file.path(out_dir, paste0(gene, "_similarity_heatmap.pdf")),
  width = 8,
  height = 8
)
pheatmap(
  cor_matrix,
  main = paste0(gene, " - Condition Similarity"),
  display_numbers = TRUE,
  number_format = "%.2f",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean"
)
dev.off()
message("Saved similarity heatmap")

# Perform hierarchical clustering
hc <- hclust(dist(t(score_matrix)))

pdf(file.path(out_dir, paste0(gene, "_dendrogram.pdf")), width = 10, height = 6)
plot(
  hc,
  main = paste0(gene, " - Condition Clustering"),
  xlab = "Condition",
  sub = ""
)
dev.off()
message("Saved dendrogram")
