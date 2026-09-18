# ----------------------------------------------------------------------------
# Fig 4 | Genome browser style coverage tracks
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/library(TxDb.Mmusculus.UCSC.mm10.R
# Input    : 1.Data/atac/Merged_replicate/bigwig/*.bigWig
# Output   : 3.Output/atac/
# Upstream : 23_atac_de_peaks_for_viz.R
# Note     : Gviz and karyoploteR on TxDb.Mmusculus.UCSC.mm10.knownGene. The original filename was a save accident
# ----------------------------------------------------------------------------

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(Gviz)
library(GenomicFeatures)
library(rtracklayer)
library(karyoploteR)
library(tidyverse)
library(cowplot)
library(export)

# 1. Setup Colors and Files ----
colPan <- c('green', 'tomato', 'green', 'tomato', 'green', 'red', 'purple')
names(colPan) <- c(
  'NDLN_DN',
  "NDLN_TDT_SP",
  "DLN_DN",
  "DLN_TDT_SP",
  "TUMOR_DN",
  "TUMOR_TDT_SP",
  "TUMOR_DP"
)

wherryDat <- list.files(
  '1.Data/atac/RefGSE86797/merged_replicate/bigwig',
  full.names = T
)[c(4, 3, 1, 2)]
wherryCol <- c('gray', 'royalblue1', 'royalblue4', 'blue4')
names(wherryCol) <- wherryDat

boolDf <- data.table::fread(
  '1.Data/atac/Merged_replicate/consensus_peaks.mRp.clN.boolean.txt'
)

gtrack <- GenomeAxisTrack(
  lwd = 0.4,
  col = 'black',
  fontcolor = "black",
  fontsize = 10
)
displayPars(gtrack) <- list(size = 2)

# 2. Target Genes ----
geneIdDf <- data.frame(
  symbols = c('Tox', 'Pdcd1'),
  id = c('252838', '18566')
)

# GRList generation
GRList <- transcriptsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = "gene")

# 3. Main Loop ----
for (gene in geneIdDf$symbols) {
  # Reset optional tracks for each loop iteration
  markerTrack <- NULL

  # Get gene info
  infoDf <- data.frame(GRList[
    geneIdDf %>% filter(symbols == gene) %>% pull(id)
  ][[1]])

  chr = unique(infoDf$seqnames)
  start = min(c(infoDf$start, infoDf$end))
  end = max(c(infoDf$start, infoDf$end))

  # Setup Reference Gene Track
  refG <- GeneRegionTrack(
    TxDb.Mmusculus.UCSC.mm10.knownGene,
    chromosome = chr,
    start = start,
    end = end,
    geneSymbol = TRUE,
    transcriptAnnotation = "gene",
    name = "Gene",
    col = 'black',
    fill = 'black',
    fontcolor.group = 'black',
    fontsize.group = 17,
    fontfamaily.group = 'Arial'
  )

  ezid <- as.numeric(geneIdDf[geneIdDf$symbols == gene, ]$id)
  refGDf <- ranges(refG) %>%
    data.frame() %>%
    filter(gene == ezid)

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

  # 4. Window Calculation & Marker Logic ----
  if (gene == 'Pdcd1') {
    # range for Pdcd1
    startNew <- start - 5000
    endNew <- end + 24000

    # --- Triangle Marker Logic ---
    gStrand <- as.character(unique(infoDf$strand)[1])

    # Calculate position of -23.8k based on strand direction
    if (gStrand == "+") {
      # If positive strand, upstream is roughly start - 23800
      markPos <- start - 23800
    } else {
      # If negative strand, upstream is roughly end + 23800
      markPos <- end + 23800
    }

    # Create the triangle track
    markerTrack <- DataTrack(
      data = 1, # Dummy y-value
      start = markPos,
      width = 500, # Width of the marker
      chromosome = chr,
      genome = "mm10",
      name = "-23.8k",
      type = "p", # points
      pch = 25, # 25 = Inverted Triangle (▼)
      cex = 1.5, # Size of the triangle
      col = "red", # Border color
      fill = "red", # Fill color
      ylim = c(0, 1),
      showAxis = FALSE,
      background.title = "transparent",
      col.title = "red",
      rotation.title = 0,
      size = 0.5 # Height of the track
    )
  } else if (gene == 'Izumo1r') {
    startNew <- start - 10000
    endNew <- end + 25000
  } else if (gene == 'Adam10') {
    startNew <- start - 25000
    endNew <- end + 25000
  } else if (gene == 'Tox') {
    peakAnno <- data.table::fread(
      '1.Data/atac/Merged_replicate/consensus_peaks.mRp.clN.annotatePeaks.txt'
    ) %>%
      rename(Geneid = 1) %>%
      filter(`Gene Name` == gene)

    gRange <- range(c(peakAnno$Start, peakAnno$End))
    startNew <- start - 150000
    endNew <- end + 50000
  } else if (gene == 'Stat3') {
    peakAnno <- data.table::fread(
      '1.Data/atac/Merged_replicate/consensus_peaks.mRp.clN.annotatePeaks.txt'
    ) %>%
      rename(Geneid = 1) %>%
      filter(`Gene Name` == gene)

    gRange <- range(c(peakAnno$Start, peakAnno$End))
    startNew <- start - 40000
    endNew <- gRange[2]
  } else {
    peakAnno <- data.table::fread(
      '1.Data/atac/Merged_replicate/consensus_peaks.mRp.clN.annotatePeaks.txt'
    ) %>%
      rename(Geneid = 1) %>%
      filter(`Gene Name` == gene)

    gRange <- range(c(peakAnno$Start, peakAnno$End))
    startNew <- ifelse(gRange[1] > start, start - 30000, gRange[1])
    endNew <- ifelse(gRange[2] < end, end + 30000, gRange[2])
  }

  # 5. Peak Import (Main Data) ----
  allScore <- c()
  for (i in c(
    "DLN_DN",
    "DLN_TDT_SP",
    "NDLN_DN",
    "NDLN_TDT_SP",
    "TUMOR_DN",
    "TUMOR_DP",
    "TUMOR_TDT_SP"
  )) {
    t <- import(
      paste0("1.Data/atac/Merged_replicate/bigwig/", i, ".mRp.clN.bigWig"),
      as = "GRanges",
      which = GRanges(c(chr), IRanges(startNew, endNew))
    )
    allScore <- c(allScore, t$score)
  }

  for (i in c(
    "DLN_DN",
    "DLN_TDT_SP",
    "NDLN_DN",
    "NDLN_TDT_SP",
    "TUMOR_DN",
    "TUMOR_DP",
    "TUMOR_TDT_SP"
  )) {
    t <- import(
      paste0("1.Data/atac/Merged_replicate/bigwig/", i, ".mRp.clN.bigWig"),
      as = "GRanges",
      which = GRanges(c(chr), IRanges(startNew, endNew))
    )

    t$score <- t$score / max(allScore)
    assign(
      paste('peak', gsub('_TDT', '', i), sep = '_'),
      DataTrack(
        t,
        name = gsub('_', '\n', gsub('_TDT', '', i)),
        fill.histogram = colPan[i],
        col.histogram = colPan[i],
        ylim = c(0, 1),
        aggregation = "mean",
        fontsize = 12,
        size = 3
      )
    )
  }

  # 6. Highlight Track Creation ----
  targetChr <- chr
  boolDfSub <- boolDf[, c(1:3, 23:30)] %>%
    filter(chr == targetChr, start > startNew, end < endNew)

  names(boolDfSub) <- gsub('.mRp.clN.qval|_LAG3', '', names(boolDfSub))

  for (i in c(
    "DLN_DN",
    "DLN_TDT_SP",
    "NDLN_DN",
    "NDLN_TDT_SP",
    "TUMOR_DN",
    "TUMOR_DP",
    "TUMOR_TDT_SP"
  )) {
    t <- boolDfSub %>%
      mutate(chr = as.numeric(gsub('chr', '', targetChr))) %>%
      select(chr, start, end, all_of(i)) %>%
      rename('value' = 4) %>%
      filter(!is.na(value))

    if (nrow(t) == 0) {
      assign(
        paste('ht', gsub('_TDT', '', i), sep = '_'),
        eval(parse(text = paste0('peak_', gsub('_TDT', '', i))))
      )
    } else {
      assign(
        paste('ht', gsub('_TDT', '', i), sep = '_'),
        HighlightTrack(
          trackList = eval(parse(text = paste0('peak_', gsub('_TDT', '', i)))),
          start = t$start,
          width = t$end - t$start,
          col = NA
        )
      )
    }
  }

  # 7. Wherry Data Import ----
  allScore2 <- c()
  for (i in wherryDat) {
    condi <- gsub('1.Data/atac/RefGSE86797/merged_replicate|.mRp.clN.bigWig', '', i)
    t <- import(
      paste0(i),
      as = "GRanges",
      which = GRanges(c(chr), IRanges(startNew, endNew))
    )
    allScore2 <- c(allScore2, t$score)
  }

  for (i in wherryDat) {
    condi <- gsub(
      '1.Data/atac/RefGSE86797/merged_replicate/bigwig/|.mRp.clN.bigWig',
      '',
      i
    )
    cat(condi, '\n')
    t <- import(
      paste0(i),
      as = "GRanges",
      which = GRanges(c(chr), IRanges(startNew, endNew))
    )
    t$score <- t$score / max(allScore2)

    assign(
      paste('peak', condi, sep = '_'),
      DataTrack(
        t,
        name = paste('Ref', condi, sep = ':\n'),
        fill.histogram = wherryCol[i],
        col.histogram = wherryCol[i],
        ylim = c(0, 1),
        aggregation = "mean",
        fontsize = 12,
        size = 2.5
      )
    )
  }

  # 8. Define Plot Lists ----

  # Define base list for No Highlight (Include markerTrack, remove later if NULL)
  trackList_NoHt <- list(
    gtrack,
    markerTrack, # Only appears if gene is Pdcd1
    # peak_NDLN_DN,   # User had this commented out
    peak_NDLN_SP,
    # peak_DLN_DN,    # User had this commented out
    peak_DLN_SP,
    # peak_TUMOR_DN,  # User had this commented out
    peak_TUMOR_DP,
    peak_TUMOR_SP,
    peak_Naive,
    peak_Memory,
    peak_Effector,
    peak_Exhausted,
    refG
  )
  # Filter out NULLs (important for markerTrack)
  trackList_NoHt <- trackList_NoHt[!sapply(trackList_NoHt, is.null)]

  # Define base list for Highlight
  trackList_Ht <- list(
    gtrack,
    markerTrack, # Only appears if gene is Pdcd1
    # ht_NDLN_DN,
    ht_NDLN_SP,
    # ht_DLN_DN,
    ht_DLN_SP,
    # ht_TUMOR_DN,
    ht_TUMOR_DP,
    ht_TUMOR_SP,
    peak_Naive,
    peak_Memory,
    peak_Effector,
    peak_Exhausted,
    refG
  )
  # Filter out NULLs
  trackList_Ht <- trackList_Ht[!sapply(trackList_Ht, is.null)]

  # 9. Plotting ----

  # Figure without highlight
  png(
    paste0('3.Output/atac/PeakViz/noht_', gene, 'fomal_PeakVizLong.png'),
    width = 4,
    height = 3,
    units = 'in',
    res = 300
  )
  plotTracks(
    trackList_NoHt,
    type = 'hist',
    col.title = 'black',
    col.axis = 'black',
    col.line = 'black',
    fontcolor = 'black',
    windowSize = 2500,
    min.height = 25,
    coverageHeight = 1,
    minCoverageHeight = 0.5,
    background.title = NA
  )
  dev.off()

  # Figure with highlight
  png(
    paste0('3.Output/atac/PeakViz/ht_', gene, 'PeakVizLong2.png'),
    width = 4,
    height = 3,
    units = 'in',
    res = 300
  )
  plotTracks(
    trackList_Ht,
    type = 'hist',
    col.title = 'black',
    col.axis = 'black',
    col.line = 'black',
    fontcolor = 'black',
    windowSize = 2500,
    min.height = 25,
    coverageHeight = 1,
    minCoverageHeight = 0.5,
    background.title = NA
  )
  dev.off()

  gc()
}
