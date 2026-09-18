# ----------------------------------------------------------------------------
# Fig 4 | Identify highly variable accessible regions
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/2.1.allHvgDetection(otherVsDN).R
# Input    : 1.Data/atac/expressionMatrixBatchEffFixed.RDS
# Output   : 1.Data/atac/expressionMatrixBatchEffFixed(highVariable).RDS
#            1.Data/atac/1.deseqDdsFile_highVariableP.RDS
# Upstream : 21_atac_deseq_init.R
# ----------------------------------------------------------------------------

if(! dir.exists("3.Output/atac/DiffExp")){
  dir.create('3.Output/atac/DiffExp')
}

library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)

# peakPosRgreat <- readRDS('1.Data/atac/peakAnnotationRgreat.RDS')
peakPosChpskr <- data.table::fread('1.Data/atac/consensus_peaks.mLb.clN.annotatePeaks.txt') %>% 
  rename(Geneid = 1)


dds <- readRDS('1.Data/atac/1.deseqDdsFile_highVariableP.RDS')

# DNLN
resNDLNSp <- results(dds, contrast=c("group", 
                                     'NDLN_TDT_SP',
                                     'NDLN_DN')) %>% 
  as.data.frame() %>% 
  merge(peakPosChpskr, by.x = 0, by.y = 'Geneid')

# DLN
resDLNSp <- results(dds, contrast=c("group", 
                                    'DLN_TDT_SP',
                                    'DLN_DN'))%>% 
  as.data.frame() %>% 
  merge(data.frame(peakPosChpskr), by.x = 0, by.y = 'Geneid')



# Tumor
resTumorSp <- results(dds, contrast=c("group", 
                                      'TUMOR_TDT_SP',
                                      'TUMOR_DN'))%>% 
  as.data.frame() %>% 
  merge(data.frame(peakPosChpskr), by.x = 0, by.y = 'Geneid')

resTumorDp <- results(dds, contrast=c("group", 
                                      'TUMOR_DP',
                                      'TUMOR_DN'))%>% 
  as.data.frame() %>% 
  merge(data.frame(peakPosChpskr), by.x = 0, by.y = 'Geneid')


write.csv(resNDLNSp, '3.Output/atac/DiffExp/DE_NDLN_SP(highVariableP).csv', row.names = F)
write.csv(resDLNSp, '3.Output/atac/DiffExp/DE_DLN_SP(highVariableP).csv', row.names = F)
write.csv(resTumorSp, '3.Output/atac/DiffExp/DE_Tumor_SP(highVariableP).csv', row.names = F)
write.csv(resTumorDp, '3.Output/atac/DiffExp/DE_Tumor_DP(highVariableP).csv', row.names = F)

# heatmap generation 
coldat <- readRDS('1.Data/atac/coldat.RDS') %>% 
  filter(!grepl('LAG3_SP', group)) %>% 
  mutate(FactorRep = paste(group, replicate, sep = '_')) %>% 
  column_to_rownames('fileName')

vsd <- vst(dds, blind=FALSE)

peakSelect <- c(resNDLNSp %>% 
                  group_by(log2FoldChange > 0) %>% 
                  top_n(10, -padj) %>% 
                  pull(Row.names),
                resDLNSp %>% 
                  group_by(log2FoldChange > 0) %>% 
                  top_n(10, -padj) %>% 
                  pull(Row.names),
                resTumorSp %>% 
                  group_by(log2FoldChange > 0) %>% 
                  top_n(10, -padj) %>% 
                  pull(Row.names),
                resTumorDp %>% 
                  group_by(log2FoldChange > 0) %>% 
                  top_n(10, -padj) %>% 
                  pull(Row.names)) %>% 
  unique()


matSub <- assay(vsd) %>% 
  as.data.frame() %>% 
  rownames_to_column('Geneid') %>%
  filter(Geneid %in% peakSelect) %>% 
  merge(peakPosChpskr %>% 
          select(Geneid, `Gene Name`, `Distance to TSS`) %>% 
          mutate(gDetail = paste0(`Gene Name`, '(', 
                                  ifelse(`Distance to TSS` > 0, '+', ''),
                                  `Distance to TSS`, ')')),
        by = 'Geneid') %>% 
  select(-c(Geneid, `Gene Name`, `Distance to TSS`)) %>% 
  column_to_rownames('gDetail')

matSub <- matSub[, c(5:8, 1:4, 9:12, 17:19, 13:16)]
# colnames(matSub) <- coldat[colnames(matSub),]$FactorRep

htmp <- BBmisc::normalize(matSub, range = c(-2, 2), margin = 1)

colnames(htmp) <- colnames(matSub)


png('3.Output/atac/unbiasedHeatmap(highVariableP).png', 
    width = 17, height = 25, units = 'in', res = 300)
Heatmap(as.matrix(htmp),
        cluster_columns = F,
        name = 'Normalized \nExpression',
        column_split = factor(gsub('_TDT', '', coldat[colnames(matSub),]$group),
                                 levels = c('NDLN_DN', 'NDLN_SP',
                                            'DLN_DN', 'DLN_SP',
                                            'TUMOR_DN', 'TUMOR_SP', 'TUMOR_DP')),
        rect_gp = gpar(col = "gray", lwd = 0.05),
        row_dend_width = unit(5, "cm"), 
        column_title_gp = grid::gpar(fontsize = 18,
                                     fontface = 2),
        
        row_names_gp =  grid::gpar(fontsize = 18,
                                   fontface = 2),
        show_column_names = F, 
        # column_names_gp = grid::gpar(fontsize = 6),
        column_gap = unit(1, "mm"),
        # width = ncol(htmp) * unit(5, 'mm'),
        # height = nrow(htmp) * unit(4, 'mm'),
        column_names_rot = 45)

dev.off()

# selected gene heatmap
peakSC <- peakPosChpskr %>% 
  filter(SYMBOL %in% c('Tcf7', 'Nfkb1', 'IL7r', 'Xcl1', 'Id3', 'Slamf6', 'Bcl6', 'Foxo1', 'Lef1')) %>% 
  pull(peakID)

peakEX <- peakPosChpskr %>% 
  filter(SYMBOL %in% c('Tox', 'Pdcd1', 'Ets1', 'Runx1', 'Runx2', 'Blimp1', 'Irf4', 'Stat5a', 'Entpd1', 'Il10', 'Cd224', 'Gzmb', 'Id', 'Batf')) %>% 
  pull(peakID)

matSub <- assay(vsd) %>% 
  as.data.frame() %>% 
  # rownames_to_column('peakID') %>%
  filter(rownames(.) %in% c(peakSC, peakEX))

htmp <- BBmisc::normalize(matSub, range = c(-2, 2), margin = 2)

colnames(htmp) <- colnames(matSub)

# png('3.Output/atac/unbiasedHeatmap.png', 
#     width = 12, height = 12, units = 'in', res = 300)
Heatmap(as.matrix(htmp),
        cluster_columns = T,
        name = 'Normalized \nExpression',
        column_split = coldat[colnames(matSub),]$Factor,
        row_split = c(rep('Stem cell like', length(peakSC)),
                      rep('Exhaustion', length(peakEX))),
        rect_gp = gpar(col = "gray", lwd = 0.05),
        row_dend_width = unit(3, "cm"), 
        column_title_gp = grid::gpar(fontsize = 10,
                                     fontface = 2),
        show_row_names = F,
        show_column_names = F, 
        # column_names_gp = grid::gpar(fontsize = 6),
        column_gap = unit(1, "mm"),
        # width = ncol(htmp) * unit(5, 'mm'),
        # height = nrow(htmp) * unit(4, 'mm'),
        column_names_rot = 45)

# dev.off()

# common peaks

cmnUP <- Reduce(intersect, list(resNDLNSp %>% 
                                  filter(padj < 0.05, log2FoldChange > 0) %>% 
                                  pull(Row.names),
                                resDLNSp %>% 
                                  filter(padj < 0.05, log2FoldChange > 0) %>% 
                                  pull(Row.names),
                                resTumorSp %>% 
                                  filter(padj < 0.05, log2FoldChange > 0) %>% 
                                  pull(Row.names),
                                resTumorDp %>% 
                                  filter(padj < 0.05, log2FoldChange > 0) %>% 
                                  pull(Row.names)))

cmnUP <- peakPosChpskr %>% filter(peakID %in% cmnUP)

cmnDN <- Reduce(intersect, list(resNDLNSp %>% 
                                  filter(padj < 0.05, log2FoldChange < 0) %>% 
                                  pull(Row.names),
                                resDLNSp %>% 
                                  filter(padj < 0.05, log2FoldChange < 0) %>% 
                                  pull(Row.names),
                                resTumorSp %>% 
                                  filter(padj < 0.05, log2FoldChange < 0) %>% 
                                  pull(Row.names),
                                resTumorDp %>% 
                                  filter(padj < 0.05, log2FoldChange < 0) %>% 
                                  pull(Row.names)))

cmnDN <- peakPosChpskr %>% filter(peakID %in% cmnDN)


t <- data.frame(resDLNSp) %>% 
  filter(padj < 0.05, log2FoldChange < 0) %>% 
  merge(peakPosRgreat, by.x = 0, by.y = 'peakID') %>% 
  merge(data.frame(peakPosChpskr), by.x = 'Row.names', by.y = 'peakID')

intersect(unique(t$SYMBOL), unique(as.character(t$annotated_genes)))

sort(table(t$SYMBOL), decreasing = T) %>% 
  head(10)

sort(table(as.character(t$annotated_genes)), decreasing = T) %>% 
  head(10)



# remove lag3+Tdt-
coldatSub <- coldat %>% 
  separate(Factor, sep = '_',into = c('Tis', 'Typ'), remove = F)


for (tiss in unique(coldatSub$Tis)) {
  
  expdatSub <- expdat[, coldatSub %>% filter(Tis == tiss) %>% rownames()]
  
  coldatT <- coldatSub %>% 
    filter(Tis == tiss) %>% 
    mutate(Factor = factor(Factor))
  
  expBatRem <- ComBat_seq(counts = as.matrix(expdat[rownames(coldatT)]), 
                          batch = coldatT$Batch, 
                          group = coldatT$Factor)
  
  dds <- DESeqDataSetFromMatrix(countData = expBatRem, 
                                # expdat[rownames(coldat)]  expBatRem
                                colData = coldatT,
                                design = ~ Batch + Factor)
  dds <- DESeq(dds)
  
  for (gntp in coldatT %>% 
       filter(Typ != 'DN')  %>% 
       pull(Typ) %>% unique()) {
    
    res <- results(dds, contrast=c("Factor", 
                                   paste(tiss, gntp, sep = '_'), 
                                   paste(tiss, 'DN', sep = '_')))
    
    saveRDS(res, paste0('3.Output/atac/2.1.TumorHVG/DEseqProj_', tiss, gntp, '.RDS'))
    
    data.frame(res) %>% 
      merge(peakPosRgreat, by.x = 0, by.y = 'peakID') %>% 
      data.table::fwrite(paste0('3.Output/atac/2.1.TumorHVG/DE_rGreat_', tiss, gntp, '.csv'))
    
    data.frame(res) %>% 
      merge(data.frame(peakPosChpskr), by.x = 0, by.y = 'peakID') %>% 
      data.table::fwrite(paste0('3.Output/atac/2.1.TumorHVG/DE_chpskr_', tiss, gntp, '.csv'))
    
  }
  
}

############# OLD ################
# tumor particular samples ----
coldatTumor <- coldat[grep('Tumor', rownames(coldat)),]
coldatTumor <- coldatTumor[-grep('tdTNLag3P', coldatTumor$Factor),]

hvgDf <- data.frame()
for(smpl in unique(coldatTumor$Factor)){
  coldatT <- coldatTumor
  
  coldatT$Factor <- as.character(coldatT$Factor)
  coldatT[coldatT$Factor != smpl,]$Factor <- 'Rest'
  
  coldatT$Factor <- as.factor(coldatT$Factor)
  
  expBatRem <- ComBat_seq(counts = as.matrix(expdat[rownames(coldatT)]), 
                          batch = coldatT$Batch, 
                          group = coldatT$Factor)
  
  dds <- DESeqDataSetFromMatrix(countData = expBatRem, 
                                # expdat[rownames(coldat)]  expBatRem
                                colData = coldatT,
                                design = ~ Batch + Factor)
  dds <- DESeq(dds)
  
  res <- results(dds, contrast=c("Factor", smpl, "Rest"))
  
  vsd <- vst(dds, blind=FALSE)
  
  # PCA plot ----
  # pcaData <- plotPCA(vsd, intgroup=c('Factor'), returnData=TRUE)
  # percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  # ggplot(pcaData, aes(PC1, PC2, col = group)) +
  #   geom_point() +
  #   ggrepel::geom_text_repel(aes(label = Factor)) +
  #   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  #   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #   coord_fixed()  +
  #   theme_bw() +
  #   theme(legend.position = 'none')
  # 
  # ggsave(paste0('3.Output/atac/2.TumorHVG/PCA_', smpl, '.png'), 
  #        width = 6, height = 5, units = 'in', bg = 'white')
  
  # genreate overall outputs ----
  resDf <- data.frame(res) %>% 
    # filter(log2FoldChange > 0,
    #        padj < 0.05) %>% 
    mutate(sample = smpl) %>% 
    rownames_to_column('Peak')
  
  hvgDf <- rbind(hvgDf, resDf)
}

hvgDf %>% 
  group_by(sample) %>% 
  filter(log2FoldChange > 0,
         padj < 0.05) %>% 
  summarise(n = n()) %>% 
  mutate(sample = gsub('Tumor_', '', sample)) %>% 
  ggplot(aes(x = sample, y = n)) +
  geom_bar(stat = 'identity') +
  scale_y_continuous(trans='log10') +
  coord_fixed(ratio = 1) +
  theme_bw()

ggsave('3.Output/atac/2.TumorHVG/hvgDist.png', 
       width = 5, height = 5, units = 'in', bg = 'white')

gCont <- data.frame(table(hvgDf$Peak))
table(gCont$Freq)

data.table::fwrite(hvgDf %>% 
                     merge(peakPos[c("peakID", 
                                     "annotated_genes",
                                     "dist_to_TSS")], 
                           by.x = 'Peak', by.y = "peakID"), 
                   '3.Output/atac/2.1.TumorHVG/tumorHvgAll.csv')

data.table::fwrite(hvgDf %>% 
                     group_by(sample) %>% 
                     filter(log2FoldChange > 0,
                            padj < 0.05) %>%
                     merge(peakPos[c("peakID", 
                                     "annotated_genes",
                                     "dist_to_TSS")], 
                           by.x = 'Peak', by.y = "peakID"), 
                   '3.Output/atac/2.1.TumorHVG/tumorHvgSig(0.05).csv')

data.table::fwrite(hvgDf %>% 
                     group_by(sample) %>% 
                     filter(log2FoldChange > 0,
                            padj < 0.1) %>%
                     merge(peakPos[c("peakID", 
                                     "annotated_genes",
                                     "dist_to_TSS")], 
                           by.x = 'Peak', by.y = "peakID"), 
                   '3.Output/atac/2.1.TumorHVG/tumorHvgSig(0.1).csv')

