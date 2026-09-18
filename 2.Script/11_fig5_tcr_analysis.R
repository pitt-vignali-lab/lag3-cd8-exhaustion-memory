# ----------------------------------------------------------------------------
# Fig 5 | TCR repertoire of the tdT subsets at d11 and d19
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/3.1.D11D19TCR.R
# Input    : 1.Data/scrna/TCR_seqNtrust_alphabetapaired.csv
#            1.Data/scrna/D11D19DatWithMonocleAnnot.RDS
# Output   : 3.Output/scrna/D11D19Dataset/
# Upstream : 05_scrna_tcr_prepare.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(export)

dat <- readRDS('1.Data/scrna/D11D19DatWithMonocleAnnot.RDS')

# tcr data read in 
tcrDf <- read_csv('1.Data/scrna/TCR_seqNtrust_betaOnly.csv') 

bcdCnt <- tcrDf %>% 
  group_by(barcode) %>% 
  summarise(cont = n())

all(bcdCnt$cont == 1)

# subset d11 d19 data only and add cdr3 data in 
datSub <- AddMetaData(dat, 
                      tcrDf %>% select(barcode, cdr3) %>% 
                        column_to_rownames('barcode'))


# summarise the TCr expension data 
cdrCnt <- dat$cdr3 %>% 
  data.frame() %>% 
  dplyr::rename('cdr3' = 1) %>% 
  filter(!is.na(cdr3)) %>% 
  group_by(cdr3) %>% 
  summarise(cont = n()) %>% 
  merge(tcrDf %>% select(barcode, cdr3),
        by = 'cdr3')


# try find a threshold to define highly/lowly expensing TCR  ----
tcrExp <- sort(table(datSub$cdr3), decreasing = T)

resdf <- data.frame()
for (i in 1:length(unique(tcrExp)) ) {
  ratio = sum(tcrExp[tcrExp >= unique(tcrExp)[i]])/sum(tcrExp)
  
  cont = length(unique(names(tcrExp[tcrExp >= unique(tcrExp)[i]])))
  
  resdf <- rbind(resdf, data.frame(cellRatio = ratio, tcrType = cont, tcrCounts = unique(tcrExp)[i]))
  
}

resdf$d <- resdf$cellRatio/resdf$tcrType # cell ratio contributed by tcr expension number. cells with tcr copy number larger or equals to 719 take contribute 3% cell numbers

resdf$cumd <- cumsum(resdf$d)

x = log10(resdf$tcrType)
resdf$y1 = -0.022*x^4+0.2164*x^3-0.696*x^2+0.7612*x-0.0749
resdf$y2 = -0.0044*x^5 + 0.0541*x^4 - 0.232*x^3 + 0.3806*x^2 - 0.0749*x + 0.0383

ggplot(resdf, aes(x = log10(tcrType), y = cumd)) +
  # geom_line() +
  geom_text(aes(label = tcrCounts), size = 2) +
  geom_hline(yintercept= 0.275, # round(max(resdf$cumd), 3),
             col = 'pink', lty = 2) +
  # geom_vline(xintercept=log10(resdf[resdf$tcrCounts == 10,]$tcrType),
  #            col = 'pink', lty = 2) +
  # geom_smooth(method = 'gam') +
  labs(y = 'Cell ratio accumulation',
       x = expression(bold(log[10]^{TCR~diversity}))) +
  theme_bw()

ggplot(resdf, aes(x = log10(tcrType), y = y1)) +
  # geom_line() +
  geom_text(aes(label = tcrCounts)) +
  geom_hline(yintercept= 0.275, # round(max(resdf$cumd), 3),
             col = 'pink', lty = 2) +
  geom_vline(xintercept=log10(resdf[resdf$tcrCounts == 10,]$tcrType),
             col = 'pink', lty = 2) +
  geom_smooth(method = 'gam') +
  labs(y = 'Cell ratio accumulation',
       x = expression(bold(log[10]^{TCR~diversity}))) +
  theme_bw()


# df for ploting ----
plotDf <- data.frame(Embeddings(datSub, 'umap'),
                     time = datSub$time,
                     genotype = datSub$genotype) %>% 
  mutate(condi = paste(time, genotype, sep = '_')) %>% 
  merge(cdrCnt, by.x = 0, by.y = 'barcode') %>% 
  mutate(tcrDef = ifelse(cont > 9, 'High expansion', 'Low expansion'))

head(plotDf)
write.csv(plotDf, '3.Output/scrna/D11D19TcrProf.csv', row.names = F)
write.csv(plotDf$Row.names, '3.Output/scrna/D11D19allTcrExpCellId.csv', row.names = F)
write.csv(plotDf[plotDf$cont > 9,]$Row.names, '3.Output/scrna/D11D19hiTcrExpCellId.csv', row.names = F)

# the percentage here is overall high tcr with their counts and overall low tcr with their counts
plotDf %>% 
  group_by(condi, tcrDef) %>%
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count),
         condi = factor(condi, levels = c("D11_CD44+tdT-", "D11_CD44+tdT+", "D11_DN", "D11_LAG3+tdT+","D11_LAG3-tdT+",
                                          "D19_CD44+tdT-", "D19_CD44+tdT+", "D19_DN", "D19_LAG3+tdT+","D19_LAG3-tdT+"))) %>% 
  ggplot(aes(x = condi, y = perc * 100, fill = gsub('expansion', '',tcrDef))) +
  geom_bar(stat="identity", col = 'black', linewidth = 0.1) +
  scale_fill_manual(values = c('orange', 'darkblue')) +
  labs(x = 'Conditions', 
       y = 'Percentage',
       fill = 'TCR Expansion') +
  theme_classic() +
  theme(legend.position = 'bottom', 
        axis.text.x = element_text(angle = 30, hjust = 1),
        text = element_text(face = 2))


plotDf %>% 
  group_by(time, genotype, tcrDef) %>%
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count),
         genotype = factor(genotype, levels = c("CD44+tdT-", "CD44+tdT+", "DN", "LAG3+tdT+","LAG3-tdT+")),
        time = factor(time, levels = c('D11', 'D19'))) %>% 
  ggplot(aes(x = genotype, y = perc * 100, fill = gsub('expansion', '',tcrDef))) +
  geom_bar(stat="identity", col = 'black', linewidth = 0.1) +
  scale_fill_manual(values = c('orange', 'darkblue')) +
  labs(x = 'Conditions', 
       y = 'Percentage',
       fill = 'TCR Expansion') +
  facet_grid(~time) +
  theme_classic() +
  theme(legend.position = 'bottom', 
        axis.text.x = element_text(angle = 30, hjust = 1),
        text = element_text(face = 2, size = 10))

export::graph2ppt(
  file   = "3.Output/scrna/D11D19Dataset/TCRAbandunce.pptx",
  width  = 3.5,     # inches
  height = 3,
  append = FALSE  # TRUE to add more slides later
)
