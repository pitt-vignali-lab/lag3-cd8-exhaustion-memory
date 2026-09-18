# ----------------------------------------------------------------------------
# Fig 5 | Chord diagram with a Pdcd1 expression ring
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/3.2.D11D19TCRSharingChordAdv1_withPD1Exp.R
# Input    : 1.Data/scrna/TCR_seqNtrust_alphabetapaired.csv
#            1.Data/scrna/D11D19DatWithMonocleAnnot.RDS
# Output   : 3.Output/scrna/D11D19Dataset/
# Upstream : 12_fig5_tcr_sharing_chord.R
# ----------------------------------------------------------------------------

library(scales)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(Seurat)

df <- read.csv('3.Output/scrna/D11D19TcrProf.csv')

sectors <- c("CD44+tdT+", "LAG3-tdT+", "LAG3+tdT+")

# Process both time points
process_timepoint <- function(time_point) {
  dfsub <- df %>%
    filter(
      condi %in%
        paste0(time_point, c('_LAG3-tdT+', '_LAG3+tdT+', '_CD44+tdT+')),
      tcrDef == 'High expansion'
    ) %>%
    group_by(cdr3) %>%
    mutate(
      cmn = length(unique(genotype)),
      cmnGt = paste(sort(unique(genotype)), collapse = '_'),
      genotype = factor(genotype, levels = sectors)
    ) %>%
    group_by(genotype) %>%
    arrange(-cmn, cmnGt, cdr3) %>%
    mutate(ord = 1:n()) %>%
    group_by(cdr3, genotype) %>%
    summarise(
      start = min(ord),
      end = max(ord),
      cmn = first(cmn),
      cmnGt = first(cmnGt),
      .groups = 'drop'
    ) %>%
    select(cdr3, genotype, cmn, cmnGt, start, end) %>%
    unique() %>%
    mutate(
      genotype = factor(genotype, levels = sectors),
      x1 = pmax(0, start - 1L),
      x2 = end
    )

  df2 <- df %>%
    filter(
      time == time_point,
      tcrDef == 'High expansion',
      genotype %in% sectors
    ) %>%
    select(cdr3, Row.names, genotype) %>%
    left_join(dfsub, by = c('cdr3', 'genotype')) %>%
    mutate(
      cat = paste0(genotype, ifelse(genotype != cmnGt, '_Share', '_Unique'))
    ) %>%
    select(Row.names, cat, cdr3) %>%
    column_to_rownames('Row.names')

  list(dfsub = dfsub, df2 = df2)
}

# Load Seurat data
dat <- readRDS('1.Data/scrna/D11D19DatWithMonocleAnnot.RDS')

# Process both timepoints
d11_data <- process_timepoint('D11')
d19_data <- process_timepoint('D19')

# Get expression data for both timepoints
get_expression_data <- function(time_point, df2_data, dfsub_data) {
  datSub <- subset(dat, cells = rownames(df2_data))
  datSub <- AddMetaData(datSub, df2_data)
  Idents(datSub) <- 'cat'

  dp1 <- DotPlot(datSub, features = c('Pdcd1'))$data

  dfplot <- df %>%
    filter(
      time == time_point,
      tcrDef == 'High expansion',
      genotype %in% sectors
    ) %>%
    select(cdr3, Row.names, genotype) %>%
    left_join(dfsub_data, by = c('cdr3', 'genotype')) %>%
    mutate(
      cat = paste0(genotype, ifelse(genotype != cmnGt, '_Share', '_Unique'))
    ) %>%
    merge(dp1 %>% select(id, avg.exp.scaled), by.x = 'cat', by.y = 'id')

  dfplot
}

dfplot_d11 <- get_expression_data('D11', d11_data$df2, d11_data$dfsub)
dfplot_d19 <- get_expression_data('D19', d19_data$df2, d19_data$dfsub)

# Calculate unified expression range across both timepoints
exp_range_unified <- range(
  c(dfplot_d11$avg.exp.scaled, dfplot_d19$avg.exp.scaled),
  na.rm = TRUE
)

# Colors
genotype_cols <- c(
  "CD44+tdT+" = "#ff6347",
  "LAG3-tdT+" = "#ff0000",
  "LAG3+tdT+" = "#800080"
)

cmnGt_cols <- setNames(
  alpha(
    c(
      '#ff6347',
      '#ff0000',
      '#800080',
      '#332D56',
      '#596FB7',
      '#71C0BB',
      '#E3EEB2'
    ),
    0.7
  ),
  c(
    "CD44+tdT+",
    "LAG3-tdT+",
    "LAG3+tdT+",
    "CD44+tdT+_LAG3-tdT+",
    "CD44+tdT+_LAG3+tdT+",
    "LAG3-tdT+_LAG3+tdT+",
    "CD44+tdT+_LAG3-tdT+_LAG3+tdT+"
  )
)

# Function to create circos plot
create_circos_plot <- function(dfsub, dfplot, exp_range) {
  xlim_df <- dfsub %>%
    group_by(genotype) %>%
    summarise(n = max(x2), .groups = "drop") %>%
    arrange(genotype)

  circos.clear()
  circos.par(start.degree = 90, gap.degree = 1, cell.padding = c(0, 0, 0, 0))

  circos.initialize(factors = xlim_df$genotype, xlim = cbind(0, xlim_df$n))
  set_track_gap(mm_h(0.2))

  # Track 1: Genotype labels
  circos.track(
    ylim = c(0, 1),
    bg.col = genotype_cols[as.character(xlim_df$genotype)],
    bg.border = NA,
    track.height = 0.05,
    panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, 3, CELL_META$sector.index)
    }
  )

  # Track 2: Expression heatmap (using unified range)
  circos.track(
    ylim = c(0, 1),
    bg.border = NA,
    track.height = 0.12,
    panel.fun = function(x, y) {
      sec <- CELL_META$sector.index
      dat <- dfplot %>% filter(genotype == sec)
      if (nrow(dat) == 0) {
        return(invisible(NULL))
      }

      for (k in seq_len(nrow(dat))) {
        xl <- max(CELL_META$xlim[1], dat$x1[k])
        xr <- min(CELL_META$xlim[2], dat$x2[k])

        if (xl < xr) {
          exp_val <- dat$avg.exp.scaled[k]
          if (!is.na(exp_val)) {
            color_idx <- round(rescale(
              exp_val,
              to = c(1, 100),
              from = exp_range
            ))
            color_idx <- max(1, min(100, color_idx))

            circos.rect(
              xl,
              0,
              xr,
              1,
              col = exp_colors[color_idx],
              border = 'white',
              lwd = 0.3 # Thinner border
            )
          }
        }
      }
    }
  )

  # Track 3: cmnGt intervals
  circos.track(
    ylim = c(0, 1),
    bg.border = NA,
    track.height = 0.05,
    panel.fun = function(x, y) {
      sec <- CELL_META$sector.index
      dat <- dfplot %>% filter(genotype == sec)
      if (nrow(dat) == 0) {
        return(invisible(NULL))
      }

      for (k in seq_len(nrow(dat))) {
        xl <- max(CELL_META$xlim[1], dat$x1[k])
        xr <- min(CELL_META$xlim[2], dat$x2[k])
        if (xl < xr) {
          circos.rect(
            xl,
            0,
            xr,
            1,
            col = cmnGt_cols[[dat$cmnGt[k]]],
            border = 'white',
            lwd = 0.3 # Thinner border
          )
        }
      }
    }
  )

  # Links for shared clonotypes
  cdr_shared <- dfsub %>% filter(cmn > 1) %>% pull(cdr3) %>% unique()
  for (cdr in cdr_shared) {
    sub <- dfsub %>% filter(cdr3 == cdr)
    if (dplyr::n_distinct(sub$genotype) < 2) {
      next
    }
    cmb <- combn(as.character(sub$genotype), 2)
    col_link <- adjustcolor(cmnGt_cols[[unique(sub$cmnGt)[1]]], alpha.f = 0.8)

    for (j in seq_len(ncol(cmb))) {
      g1 <- cmb[1, j]
      g2 <- cmb[2, j]
      x1 <- c(
        sub %>% filter(genotype == g1) %>% pull(x1),
        sub %>% filter(genotype == g1) %>% pull(x2)
      )
      x2 <- c(
        sub %>% filter(genotype == g2) %>% pull(x1),
        sub %>% filter(genotype == g2) %>% pull(x2)
      )
      circos.link(g1, x1, g2, x2, col = col_link, border = "white", lwd = 0.2) # Thinner border
    }
  }
}

exp_colors <- colorRampPalette(c("#E0FFE0", '#00a300', "#003300"))(100)

# Create unified legend
lgd <- Legend(
  col_fun = colorRamp2(
    breaks = seq(exp_range_unified[1], exp_range_unified[2], length.out = 100),
    colors = exp_colors
  ),
  title = "Pdcd1\nExpression",
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 8),
  direction = "horizontal",
  legend_width = unit(2.5, "cm")
)

# Plot D11
create_circos_plot(d11_data$dfsub, dfplot_d11, exp_range_unified)
draw(lgd, x = unit(0.2, "npc"), y = unit(0.1, "npc"))
export::graph2ppt(
  file = "3.Output/scrna/D11D19Dataset/d11cmnTCR.ppt",
  width = 5,
  height = 5,
  bg = 'transparent'
)
export::graph2png(
  file = "3.Output/scrna/D11D19Dataset/d11cmnTCR.png",
  width = 2.5,
  height = 2.5,
  bg = 'transparent'
)

# Plot D19
create_circos_plot(d19_data$dfsub, dfplot_d19, exp_range_unified)
draw(lgd, x = unit(0.2, "npc"), y = unit(0.1, "npc"))
export::graph2ppt(
  file = "3.Output/scrna/D11D19Dataset/d19cmnTCR.ppt",
  width = 5,
  height = 5,
  bg = 'transparent'
)
export::graph2png(
  file = "3.Output/scrna/D11D19Dataset/d19cmnTCR.png",
  width = 5,
  height = 5,
  bg = 'transparent'
)
