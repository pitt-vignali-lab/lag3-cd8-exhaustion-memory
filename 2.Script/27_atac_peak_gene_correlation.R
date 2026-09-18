# ----------------------------------------------------------------------------
# Analysis | Correlate peak accessibility with gene expression
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/2.3.sigDeGCoorelationTest.R
# Input    : 1.Data/atac/
#            1.Data/scrna/
# Output   : 3.Output/atac/
# Upstream : 25_atac_highly_variable_peaks.R
# ----------------------------------------------------------------------------

library(rtracklayer)
library(tidyverse)
library(Hmisc)
library(corrplot)
library(psych)
library(ggcorrplot2)
library(export)

# significantly expressed genes ----
deG <- unique(unlist(sapply(
  list.files('3.Output/atac/DiffExp/', full.names = T)[c(1:3, 5)],
  \(x) read.csv(x) %>% filter(padj < 0.05) %>% pull(Row.names)
)))

deG <- unique(unlist(sapply(
  list.files('3.Output/atac/DiffExp/', full.names = T)[c(1, 3, 5, 8)],
  \(x) read.csv(x) %>% filter(padj < 0.05) %>% pull(Row.names)
)))


colDat <- readRDS('1.Data/atac/coldat.RDS') %>%
  filter(!grepl('LAG3_SP', group)) %>%
  mutate(group = gsub('_TDT', '', group))

# original Exp ----
# expDat <- read.table('1.Data/atac/consensus_peaks.mLb.clN.featureCounts.txt',
#                      header = T) %>%
#   column_to_rownames('Geneid') %>%
#   dplyr::select(-c(1:5))

# modified Exp ----
expDat <- readRDS('1.Data/atac/expressionMatrixBatchEffFixed(highVariable).RDS')


transDat <- do.call(
  cbind,
  lapply(unique(colDat$group), function(x) {
    rowMeans(expDat[, colDat[colDat$group == x, ]$fileName])
  })
)
colnames(transDat) <- unique(colDat$group)

ct <- corr.test(
  transDat[deG, c(3, 4, 1, 2, 5, 6, 7)] %>%
    as.data.frame() %>%
    rename_with(~ gsub("_", "\n", .x, fixed = TRUE)),
  use = "pairwise",
  method = "spearman",
  adjust = "BH",
  alpha = .05,
  ci = TRUE,
  minlength = 5,
  normal = TRUE
)

ct <- corr.test(
  transDat[deG, c(4, 2, 6, 7)] %>%
    as.data.frame() %>%
    rename_with(~ gsub("_", "\n", .x, fixed = TRUE)),
  use = "pairwise",
  method = "spearman",
  adjust = "BH",
  alpha = .05,
  ci = TRUE,
  minlength = 5,
  normal = TRUE
)

corr <- ct$r
p.mat <- ct$p

# Visualize the correlation matrix
# --------------------------------
# method = "circle" (default)

corrplot.mixed(corr, lower.col = "black", number.cex = 0.7)


t <- ggcorrplot.mixed(
  corr,
  upper = "ellipse",
  lower = "number",
  p.mat = p.mat,
  # tl.col="black",
  insig = "label_sig",
  sig.lvl = c(0.05, 0.01, 0.001)
)

t[["theme"]][["text"]][["face"]] <- 2
t[["theme"]][["text"]][["color"]] <- 'black'

t +
  scale_fill_gradientn(
    colours = rev(t[["plot_env"]][["col"]]),
    limits = c(-1, 1),
    guide = guide_colorbar(
      direction = "vertical",
      title = "",
      nbin = 1000,
      ticks.colour = "black",
      frame.colour = "black",
      barwidth = 1.5,
      barheight = 15
    )
  ) +
  scale_colour_gradientn(
    colours = rev(t[["plot_env"]][["col"]]),
    limits = c(-1, 1),
    guide = guide_colorbar(
      direction = "vertical",
      title = "",
      nbin = 1000,
      ticks.colour = "black",
      frame.colour = "black",
      barwidth = 1.5,
      barheight = 15
    )
  ) +
  coord_fixed() +
  theme(legend.position = "right", text = element_text(size = 15, face = 2))

ggsave(
  '3.Output/atac/sigDeGCoorelation(highVariable)_new.png',
  width = 5,
  height = 4.5,
  units = 'in'
)
graph2ppt(
  file = "3.Output/atac/sigDeGCoorelation(highVariable)_new.pptx",
  width = 5,
  height = 4.5
)
