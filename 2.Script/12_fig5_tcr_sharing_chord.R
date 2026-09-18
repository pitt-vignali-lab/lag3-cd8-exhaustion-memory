# ----------------------------------------------------------------------------
# Fig 5 | Chord diagram of clonotype sharing between tdT subsets
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/3.2.D11D19TCRSharingChordAdv1.R
# Input    : 1.Data/scrna/TCR_seqNtrust_alphabetapaired.csv
# Output   : 3.Output/scrna/D11D19Dataset/d11cmnTCR.png
#            3.Output/scrna/D11D19Dataset/d19cmnTCR.png
# Upstream : 11_fig5_tcr_analysis.R
# Note     : Drawn with circlize. This panel is also the repository logo
# ----------------------------------------------------------------------------

library(circlize)
library(tidyverse)
library(cowplot)

# Arrange the layout
# layout(matrix(1, 1, 1)) 

df <- read.csv('3.Output/scrna/D11D19TcrProf.csv') 

sectors <- c("CD44+tdT+", "LAG3-tdT+", "LAG3+tdT+")

dfsub <- df %>% 
  filter(condi %in% paste0('D11', c('_LAG3-tdT+', '_LAG3+tdT+', '_CD44+tdT+')),
         tcrDef == 'High expansion') %>% 
  group_by(cdr3) %>% 
  mutate(cmn = length(unique(genotype)),
         cmnGt = paste(sort(unique(genotype)), collapse = '_'),
         genotype = factor(genotype, levels = c('CD44+tdT+', 'LAG3-tdT+', 'LAG3+tdT+'))) %>% 
  group_by(genotype) %>% 
  arrange(-cmn, cmnGt, cdr3) %>% 
  mutate(ord = 1: n()) %>% 
  group_by(cdr3, genotype) %>% 
  summarise(start = min(ord),
            end = max(ord),
            cmn = first(cmn), 
            cmnGt = first(cmnGt),
            .groups = 'drop') %>% 
  select(cdr3, genotype, cmn, cmnGt, start, end) %>% 
  unique() %>%
  mutate(
    genotype = factor(genotype, levels = sectors),
    x1 = pmax(0, start - 1L),   # left boundary for intervals
    x2 = end
  )

# xlim per sector
xlim_df <- dfsub %>%
  group_by(genotype) %>%
  summarise(n = max(x2), .groups = "drop") %>%
  arrange(genotype)

# --- colors ---
genotype_cols <- c("CD44+tdT+" = "#ff6347",
                   "LAG3-tdT+" = "#ff0000",
                   "LAG3+tdT+" = "#800080")

# distinct cmnGt groups -> assign colors
# cmnGt_levels <- sort(unique(dfsub$cmnGt))
cmnGt_cols <- setNames(
  alpha(c('#ff6347', '#ff0000', '#800080', '#332D56', '#596FB7', '#71C0BB', '#E3EEB2'), 0.7),
  c("CD44+tdT+", "LAG3-tdT+", "LAG3+tdT+", "CD44+tdT+_LAG3-tdT+", "CD44+tdT+_LAG3+tdT+", "LAG3-tdT+_LAG3+tdT+", "CD44+tdT+_LAG3-tdT+_LAG3+tdT+")
)

# --- initialize ---
circos.clear()
circos.par(start.degree = 90, 
           gap.degree = 1,
           # track.margin = c(0, 0),
           cell.padding = c(0, 0, 0, 0)
           ) 

circos.initialize(factors = xlim_df$genotype,
                  xlim    = cbind(0, xlim_df$n))
set_track_gap(mm_h(0.2))
# --- Track 1: genotype strip (background by genotype) ---
circos.track(
  ylim = c(0, 1),
  bg.col    = genotype_cols[as.character(xlim_df$genotype)],
  bg.border = NA,
  track.height = 0.05,
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter, 
      3, # CELL_META$ycenter, 
      CELL_META$sector.index
    )
    # optional axis per sector:
    # circos.axis(h = 0, major.at = pretty(CELL_META$xlim), labels.cex = 0.5)
  }
)

# --- Track 2: cmnGt intervals within each genotype sector ---
# Paint each (x1, x2) interval by its cmnGt color
circos.track(
  ylim = c(0, 1),
  bg.border = NA,
  track.height = 0.035,
  panel.fun = function(x, y) {
    sec <- CELL_META$sector.index
    dat <- dfsub %>% filter(genotype == sec)
    if (nrow(dat) == 0) return(invisible(NULL))
    for (k in seq_len(nrow(dat))) {
      # clamp just in case:
      xl <- max(CELL_META$xlim[1], dat$x1[k])
      xr <- min(CELL_META$xlim[2], dat$x2[k])
      if (xl < xr) {
        circos.rect(xl, 0, xr, 1,
                    col = cmnGt_cols[[ dat$cmnGt[k] ]],
                    border = 'white')
      }
    }
  }
)

# --- (Optional) Draw links between sectors, colored by cmnGt ---
# Only for shared clonotypes if you want:
cdr_shared <- dfsub %>% filter(cmn > 1) %>% pull(cdr3) %>% unique()
for (cdr in cdr_shared) {
  sub <- dfsub %>% filter(cdr3 == cdr)
  if (dplyr::n_distinct(sub$genotype) < 2) next
  cmb <- combn(as.character(sub$genotype), 2)
  # Use the group color once per cdr3 (assumes one cmnGt per cdr3; tweak if not)
  col_link <- adjustcolor(cmnGt_cols[[ unique(sub$cmnGt)[1] ]], alpha.f = 0.8)

  for (j in seq_len(ncol(cmb))) {
    g1 <- cmb[1, j]; g2 <- cmb[2, j]
    x1 <- c(sub %>% filter(genotype == g1) %>% pull(x1),
            sub %>% filter(genotype == g1) %>% pull(x2))
    x2 <- c(sub %>% filter(genotype == g2) %>% pull(x1),
            sub %>% filter(genotype == g2) %>% pull(x2))
    circos.link(
      g1, x1, g2, x2,
      col = col_link,
      border = "white"
    )
  }
}

export::graph2ppt(
  file   = "3.Output/scrna/D11D19Dataset/d11cmnTCR.pptx",
  width  = 4,     # inches
  height = 4,
  append = FALSE  # TRUE to add more slides later
)

export::graph2png(
  file   = "3.Output/scrna/D11D19Dataset/d19cmnTCR.png",
  width  = 25,     # inches
  height = 25,
  bg = 'transparent'
)

# (Optional) legend for cmnGt
# if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
#   lgd <- ComplexHeatmap::Legend(title = "cmnGt", labels = names(cmnGt_cols),
#                                 legend_gp = grid::gpar(fill = unname(cmnGt_cols)))
#   grid::grid.draw(lgd)
# }
