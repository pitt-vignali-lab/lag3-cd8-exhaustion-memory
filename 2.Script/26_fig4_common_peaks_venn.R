# ----------------------------------------------------------------------------
# Fig 4 | Venn comparison of accessible regions between subsets
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/2.2.commonPeaks(Venn).R
# Input    : 1.Data/atac/expressionMatrixBatchEffFixed.RDS
# Output   : 3.Output/atac/
# Upstream : 25_atac_highly_variable_peaks.R
# ----------------------------------------------------------------------------

library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
library(export)

# orgs that are not all opening or all closed are collected 
peakBool <- data.table::fread('1.Data/atac/MergedReplicate/consensus_peaks.mRp.clN.boolean.txt')
hskp <- read.csv('3.Output/atac/1.housekeepingPeaks.csv')
annoDf <- data.table::fread('1.Data/atac/MergedReplicate/consensus_peaks.mRp.clN.annotatePeaks.txt') %>% 
  rename('interval_id' = 1)

peakBool <- peakBool %>% 
  filter(!interval_id %in% hskp$x)

bool <- names(peakBool)[grep('bool', names(peakBool))]

fc <- names(peakBool)[grep('fc', names(peakBool))]

pval <- names(peakBool)[grep('pval', names(peakBool))]

qval <- names(peakBool)[grep('qval', names(peakBool))]

subDf <- data.frame(peakBool$interval_id,
                    peakBool %>% select(bool[1],
                                        fc[1],
                                        pval[1],
                                        qval[1]))
peakSum <- peakBool %>% 
  select('interval_id', bool, -'TUMOR_LAG3_SP.mRp.clN.bool') %>% 
  column_to_rownames('interval_id') %>% 
  mutate(DLN_TDT = DLN_DN.mRp.clN.bool == DLN_TDT_SP.mRp.clN.bool & DLN_TDT_SP.mRp.clN.bool, # DLN opening region
         NDLN_TDT = NDLN_DN.mRp.clN.bool == NDLN_TDT_SP.mRp.clN.bool & NDLN_TDT_SP.mRp.clN.bool,
         Tumor_TDT = TUMOR_DN.mRp.clN.bool == TUMOR_TDT_SP.mRp.clN.bool & TUMOR_TDT_SP.mRp.clN.bool,
         Tumor_DP = TUMOR_DN.mRp.clN.bool == TUMOR_DP.mRp.clN.bool & TUMOR_DP.mRp.clN.bool)

library(ggvenn)

ggvenn(peakSum %>% 
         select(DLN_TDT, NDLN_TDT, Tumor_TDT, Tumor_DP) %>% 
         filter(rowSums(.) != 0))

graph2ppt(file='3.Output/atac/positivePopulationOCRVenn(highVariableP).pptx', width=8, height=5)

ggsave('3.Output/atac/positivePopulationOCRVenn(highVariableP).png',
       width = 8, height = 5)



commonPeak <- peakSum %>% 
  select(DLN_TDT, NDLN_TDT, Tumor_TDT, Tumor_DP) %>% 
  filter(rowSums(.) == 4) %>% 
  rownames()

commonPeakAnno <- annoDf %>% filter(interval_id %in% commonPeak)

sort(table(commonPeakAnno$`Gene Name`), decreasing = T)


phyDf <- data.frame()

N = nrow(annoDf) # sum(table(annoDf$`Gene Name`))
n = nrow(commonPeakAnno)

for (i in unique(commonPeakAnno$`Gene Name`)) {
  m = table(commonPeakAnno$`Gene Name`)[i]
  M = table(annoDf$`Gene Name`)[i]
  
  phyDf <- rbind(phyDf,
                 data.frame(gene = i,
                            pval = phyper(m,M,(N-M),n,lower.tail=F)))
}

t <- data.frame(gene = unique(commonPeakAnno$`Gene Name`)) %>% 
  mutate(N = nrow(annoDf),
         n = nrow(commonPeakAnno))

t$m <- table(commonPeakAnno$`Gene Name`)[t$gene]
t$M <- table(annoDf$`Gene Name`)[t$gene]

t$pval <- apply(t, 1, function(x) phyper(as.numeric(x['m']), as.numeric(x['M']),(as.numeric(x['N'])-as.numeric(x['M'])),as.numeric(x['n']),lower.tail=F))
t$fdr <- p.adjust(t$pval, method = 'fdr')

apply(t[1, ], 1, \(x) as.numeric(x['m']) + as.numeric(x['n']))



t$pval <- sapply(t$gene, 
       function(x) phyper(table(commonPeakAnno$`Gene Name`)[x],
                          table(annoDf$`Gene Name`)[x],
                          (nrow(annoDf)-table(annoDf$`Gene Name`)[x]),
                          nrow(commonPeakAnno),lower.tail=F))

intersect(t[t$fdr < 0.05, ]$gene, comG)




table(commonPeakAnno$`Gene Name`)['Pdcd1']

length(unique(commonPeakAnno$`Gene Name`))

intersect(comG, unique(commonPeakAnno$`Gene Name`))
