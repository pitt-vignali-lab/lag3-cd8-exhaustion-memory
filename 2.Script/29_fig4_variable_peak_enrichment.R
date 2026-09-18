# ----------------------------------------------------------------------------
# Fig 4 | Functional enrichment of the highly variable peaks
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/3.1.HighvariablePeakFuncEnrich.R
# Input    : 1.Data/atac/expressionMatrixBatchEffFixed(highVariable).RDS
# Output   : 3.Output/atac/
# Upstream : 25_atac_highly_variable_peaks.R
# ----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(DESeq2)
library(GSVA)
library(msigdbr)
library(pheatmap)

# 1) Inputs from nf-core/atacseq
ann <- read_tsv("1.Data/atac/consensus_peaks.mLb.clN.annotatePeaks.txt")      # HOMER annotations

# ---- choose ONE representative peak per gene (closest-to-TSS + tie-breakers) ----
# Adjust 'gene_col' to what you have: e.g., "Gene Name" or "Entrez ID" or "Nearest PromoterID"
gene_col <- if ("Gene Name" %in% names(ann)) "Gene Name" else "Entrez ID"

ann2 <- ann %>%
  transmute(
    peak_id = .data[[names(ann)[1]]],
    gene    = .data[[gene_col]],
    dist    = suppressWarnings(as.numeric(`Distance to TSS`)),
    annot   = Annotation,
    score   = suppressWarnings(as.numeric(`Peak Score`))
  ) %>%
  filter(!is.na(gene), gene != "") %>%
  mutate(
    abs_dist   = abs(dist),
    annot_rank = case_when(
      grepl("promoter|TSS", annot, ignore.case=TRUE) ~ 1L,
      grepl("exon", annot, ignore.case=TRUE)         ~ 2L,
      grepl("intron", annot, ignore.case=TRUE)       ~ 3L,
      grepl("downstream", annot, ignore.case=TRUE)   ~ 4L,
      TRUE                                           ~ 5L
    )
  )

rep_peaks <- ann2 %>%
  group_by(gene) %>%
  arrange(abs_dist, annot_rank, desc(score), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(gene, peak_id)

# 2) Build a gene × sample count matrix from featureCounts (replicate-aware)
fc <- readRDS('1.Data/atac/expressionMatrixBatchEffFixed(highVariable).RDS')  %>%      # counts per peak per sample
  data.frame() %>% 
  rownames_to_column('Geneid')

sample_cols <- setdiff(names(fc), c("Geneid","Chr","Start","End","Strand","Length"))

gene_cts <- fc %>%
  semi_join(rep_peaks, by = c("Geneid" = "peak_id")) %>%
  left_join(rep_peaks, by = c("Geneid" = "peak_id")) %>%
  select(gene, all_of(sample_cols)) %>%
  group_by(gene) %>% summarise(across(all_of(sample_cols), sum), .groups="drop") %>%
  tibble::column_to_rownames("gene") %>% as.matrix()

dds <- DESeqDataSetFromMatrix(round(gene_cts), DataFrame(row.names=colnames(gene_cts)), ~ 1)
dds <- estimateSizeFactors(dds)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
mat <- assay(vsd)                      # genes × samples

raw_names  <- colnames(mat)
base_names <- sub("\\.mLb\\.clN\\.sorted\\.bam$", "", raw_names)

mm <- str_match(base_names, "^(DLN|NDLN|TUMOR)_(DN|DP|TDT(?:_SP)?)_REP(\\d+)$")
stopifnot(all(!is.na(mm[,1])))

meta <- tibble(
  sample_raw  = raw_names,
  sample_base = base_names,
  tissue      = mm[,2],
  subset_raw  = mm[,3],
  rep         = as.integer(mm[,4])
) %>%
  mutate(
    # alias TDT_SP -> TDT
    subset = recode(subset_raw, "TDT_SP" = "TDT", .default = subset_raw),
    group  = paste(tissue, subset, sep = "_")
  )

# Compute group means (per gene)
idx_by_group <- split(seq_len(ncol(mat)), meta$group)
mat_groupmean <- sapply(idx_by_group, function(ix) {
  rowMeans(mat[, ix, drop = FALSE], na.rm = TRUE)
})

# order groups for readability
desired_order <- c("NDLN_DN","NDLN_TDT","DLN_DN","DLN_TDT","TUMOR_DN","TUMOR_TDT","TUMOR_DP")
present_order <- intersect(desired_order, colnames(mat_groupmean))
mat_groupmean <- mat_groupmean[, present_order, drop = FALSE]

# 4A) Per-sample pathway activity (GSVA/ssGSEA) — NO CONTROL REQUIRED
# Use gene symbols in msigdbr if your 'gene' are symbols; otherwise switch to entrez_gene.
msig <- msigdbr(species = "Mus musculus", category = "C7")  # "Homo sapiens" for human
paths <- split(msig$gene_symbol, msig$gs_name)

# 4B) Optional “global” enrichment (interpret as ‘generally open promoters’)
# Rank genes by mean normalized accessibility and run pre-ranked GSEA.
library(fgsea)

# 1) Z-score each gene across groups (row-wise scale)
gm_z <- t(scale(t(as.matrix(mat_groupmean)), center = TRUE, scale = TRUE))

# 2) Run fgsea per group with stable settings
register(SerialParam())  # or use nproc = 1 inside fgsea()
set.seed(1)

fg <- lapply(colnames(gm_z), function(grp) {
  stats <- gm_z[, grp]
  # break any residual ties without changing order
  stats <- stats + rnorm(length(stats), sd = 1e-8)
  stats <- sort(stats, decreasing = TRUE)

  fgsea(pathways = paths,
        stats    = stats,
        minSize  = 15,
        maxSize  = 500,
        scoreType = "std",      # now valid because we have pos & neg
        nproc     = 1)          # avoids serialize/connection errors
})

names(fg) <- colnames(gm_z)

saveRDS(fg, '1.Data/atac/fucntionalEnrichmentResult.RDS')

fg$NDLN_TDT %>% 
  filter(padj < 0.05)

fg$DLN_TDT %>% 
  filter(padj < 0.05) %>% 
  filter(grepl('CD8', pathway))

fg$TUMOR_TDT %>% 
  filter(padj < 0.05) %>% 
  filter(grepl('CD8', pathway))

fg$TUMOR_DP %>% 
  filter(padj < 0.05) %>% 
  filter(grepl('CD8', pathway))

fg$DLN_TDT %>% 
  filter(padj < 0.05) %>% 
  filter(grepl('CD8', pathway)) %>%
  arrange(padj) %>%
  mutate(theme = case_when(
    str_detect(pathway, "NAIVE") ~ "Naive",
    str_detect(pathway, "MEMORY") ~ "Memory",
    str_detect(pathway, "EFFECTOR|ACT") ~ "Effector/Activation",
    str_detect(pathway, "EXHAUST") ~ "Exhaustion",
    str_detect(pathway, "IFN|IFNA|IFNAR|STAT") ~ "Interferon",
    TRUE ~ "Other"
  )) %>%
  group_by(theme) %>% slice_min(padj, n = 3, with_ties = FALSE) %>% ungroup() %>%
  ggplot(aes(x = reorder(pathway, NES), y = NES, fill = theme)) +
  geom_col() + coord_flip() +
  labs(title = "Top enriched immunologic signatures", x = NULL, y = "NES")
library(ggvenn)

ggvenn(list(`DLN_TDT` = fg$DLN_TDT %>% 
  filter(padj < 0.05) %>% 
  filter(grepl('CD8', pathway)) %>% pull(pathway),

`TUMOR_TDT` = fg$TUMOR_TDT %>% 
  filter(padj < 0.05) %>% 
  filter(grepl('CD8', pathway))%>% pull(pathway),

`TUMOR_DP` = fg$TUMOR_DP %>% 
  filter(padj < 0.05) %>% 
  filter(grepl('CD8', pathway))%>% pull(pathway)))+
  theme(text = element_text(face = 2, size = 8))

ggsave('3.Output/atac/fucntionVenn.png', width = 6, height = 6)
export::graph2ppt(file = '3.Output/atac/fucntionVenn.ppt', width = 5, height = 5)

rbind(fg$DLN_TDT %>% 
  filter(padj < 0.05) %>% 
  filter(grepl('CD8', pathway)) %>% 
  top_n(10, -padj) %>% 
  mutate(cat = 'DLN TDT'),

fg$TUMOR_TDT %>% 
  filter(padj < 0.05) %>% 
  filter(grepl('CD8', pathway))%>% 
  top_n(10, -padj)%>% 
  mutate(cat = 'TUMOR TDT'),

fg$TUMOR_DP %>% 
  filter(padj < 0.05) %>% 
  filter(grepl('CD8', pathway))%>% 
  top_n(10, -padj)%>% 
  mutate(cat = 'TUMOR DP')) %>% 
  ggplot(aes(x = NES, y = reorder(pathway, -NES))) +
  geom_bar(stat = 'identity') +
   facet_grid(cat ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(text = element_text(face = 2, size = 8))

export::graph2ppt(file = '3.Output/atac/fucntionBar.ppt', width = 6, height = 5)
