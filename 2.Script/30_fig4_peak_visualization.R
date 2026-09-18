# ----------------------------------------------------------------------------
# Fig 4 | Peak-level visualisation across the subsets
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/3.peakViz.R
# Input    : 1.Data/atac/Merged_replicate/bigwig/
#            1.Data/atac/MergedReplicate/consensus_peaks.mRp.clN.annotatePeaks.txt
# Output   : 3.Output/atac/
# Upstream : 23_atac_de_peaks_for_viz.R
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

# target genes details ----
geneIdDf <- data.frame(
  symbols = c('Tox', 'Pdcd1'),
  id = c('252838', '18566')
)

# memory set
# symbols = c('Tcf7', 'Il7r', 'Bcl2', 'Eomes', 'Id3', 'Bach2', 'Ccr7', 'Sell', 'Cd44', 'Itgae', 'Cd69'),
# id = c('21414', '16197', '12043', '13813', '15903', '12014', '12775', '20343', '12505', '16408', '12515')

GRList <- transcriptsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = "gene")


# gene = 'Tox'

for (gene in geneIdDf$symbols) {
  infoDf <- data.frame(GRList[
    geneIdDf %>% filter(symbols == gene) %>% pull(id)
  ][[1]])

  chr = unique(infoDf$seqnames)
  start = min(c(infoDf$start, infoDf$end))
  end = max(c(infoDf$start, infoDf$end))

  # reference gene track ----
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
  # z <- ranges(refG)
  # z$symbol <- gene
  # z$transcript <- gene
  # z$gene <- gene
  # ranges(refG) <- z
  # displayPars(refG) <- list(size=0.8)

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

  if (gene == 'Pdcd1') {
    # range for Pdcd1
    startNew <- start - 5000
    endNew <- end + 24000
  } else if (gene == 'Izumo1r') {
    # range for Pdcd1
    startNew <- start - 10000
    endNew <- end + 25000
  } else if (gene == 'Adam10') {
    # range for Pdcd1
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
    # range for the rest
    peakAnno <- data.table::fread(
      '1.Data/atac/Merged_replicate/consensus_peaks.mRp.clN.annotatePeaks.txt'
    ) %>%
      rename(Geneid = 1) %>%
      filter(`Gene Name` == gene)

    gRange <- range(c(peakAnno$Start, peakAnno$End))

    startNew <- ifelse(gRange[1] > start, start - 30000, gRange[1])
    endNew <- ifelse(gRange[2] < end, end + 30000, gRange[2])
  }

  # peaks collect ----

  # "p", "l", "b", "a", "s", "g", "r", "S", "smooth", "polygon", "horizon",
  # "histogram", "mountain", "h", "boxplot", "gradient", "heatmap", "confint"

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
      # AlignmentsTrack DataTrack
      DataTrack(
        t,
        name = gsub('_', '\n', gsub('_TDT', '', i)),
        fill.histogram = colPan[i],
        col.histogram = colPan[i],
        ylim = c(0, 1), # pd1: 1.6 # Tox 1
        aggregation = "mean",
        fontsize = 12,
        size = 3
      )
    )
  }

  # add highlight

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

  # add wherry data
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
      # AlignmentsTrack DataTrack
      DataTrack(
        t,
        name = paste('Ref', condi, sep = ':\n'),
        fill.histogram = wherryCol[i],
        col.histogram = wherryCol[i],
        ylim = c(0, 1), # pd1: 1.6 # Tox 1
        aggregation = "mean",
        fontsize = 12,
        size = 2.5
      )
    )
  }

  # figure without hightlight ----

  png(
    paste0('3.Output/atac/PeakViz/noht_', gene, 'fomal_PeakVizLong.png'),
    width = 4,
    height = 3,
    units = 'in',
    res = 300
  )
  plotTracks(
    list(
      gtrack,
      # peak_NDLN_DN,
      peak_NDLN_SP,
      # peak_DLN_DN,
      peak_DLN_SP,
      # peak_TUMOR_DN,
      peak_TUMOR_DP,
      peak_TUMOR_SP,
      peak_Naive,
      peak_Memory,
      peak_Effector,
      peak_Exhausted,
      refG
    ),
    type = 'hist',
    col.title = 'black',
    col.axis = 'black',
    col.line = 'black',
    fontcolor = 'black',
    windowSize = 2500,
    min.height = 25,
    coverageHeight = 0.08,
    minCoverageHeight = 0,
    coverageHeight = 1,
    minCoverageHeight = 0.5,
    background.title = NA
  )

  dev.off()

  # figure with hightlight ----
  png(
    paste0('3.Output/atac/PeakViz/ht_', gene, 'PeakVizLong.png'),
    width = 4,
    height = 3,
    units = 'in',
    res = 300
  )
  plotTracks(
    list(
      gtrack,
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
    ),
    type = 'hist',
    col.title = 'black',
    col.axis = 'black',
    col.line = 'black',
    fontcolor = 'black',
    windowSize = 2500,
    min.height = 25,
    coverageHeight = 0.08,
    minCoverageHeight = 0,
    coverageHeight = 1,
    minCoverageHeight = 0.5,
    background.title = NA
  )

  dev.off()
  gc()
}
