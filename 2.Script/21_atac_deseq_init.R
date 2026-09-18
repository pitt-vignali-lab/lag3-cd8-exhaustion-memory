# ----------------------------------------------------------------------------
# Setup | Build the ATAC-seq DESeq2 object from consensus peak counts
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/1.1.deInitialization.R
# Input    : 1.Data/atac/consensus_peaks.mRp.clN.featureCounts.txt
#            1.Data/atac/coldat.RDS
# Output   : 1.Data/atac/expressionMatrix.RDS
#            1.Data/atac/expressionMatrixBatchEffFixed.RDS
# Note     : Counts batch-corrected with ComBat-seq (sva), which preserves integer counts
# ----------------------------------------------------------------------------

library(sva)
library(DESeq2)
library(tidyverse)
library(export)

output <- read.table('1.Data/atac/consensus_peaks.mLb.clN.featureCounts.txt',
                     header = T)
rownames(output) <- output$Geneid

# get expression matrix
expdat <- output[-c(1:6)]
saveRDS(expdat, '1.Data/atac/expressionMatrix.RDS')


coldat <- read.csv('1.Data/atac/expSet.csv')
saveRDS(coldat, '1.Data/atac/coldat.RDS')
coldat <- coldat[-grep('LAG3_SP', coldat$group),]


expdat <- expdat[, coldat$fileName]
expBatRem <- ComBat_seq(counts=as.matrix(expdat), 
                        batch = coldat$batch)
saveRDS(expBatRem, '1.Data/atac/expressionMatrixBatchEffFixed.RDS')

# housekeeping peak detection 
# Define a small constant
epsilon <- 1e-6

# Define a function to calculate mean and regularized CV
calculate_regularized_cv <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  cv <- sd_x / (mean_x + epsilon)
  return(c(mean_x, cv))
}

# Apply the function to each row of the expression matrix
# Create a data frame with peak IDs and their corresponding CV values
cv_df <- t(apply(expBatRem, 1, calculate_regularized_cv)) %>% 
  as.data.frame() %>% 
  dplyr::rename('mean' = 1, CV = 2) %>% 
  rownames_to_column('PeakID') %>% 
  mutate(logCV = log(CV))

# Create a Q-Q plot to visually inspect normality
qqnorm(cv_df$logCV)
qqline(cv_df$logCV, col = "red")

# Perform the Kolmogorov-Smirnov test for normality
ks_test <- ks.test(cv_df$logCV, 
                   "pnorm", 
                   mean = mean(cv_df$logCV, na.rm = TRUE), 
                   sd = sd(cv_df$logCV, na.rm = TRUE))
print(ks_test) # thre result here indicate the distribution is not normality 

# Asymptotic one-sample Kolmogorov-Smirnov test
# 
# data:  cv_values
# D = 0.035327, p-value < 2.2e-16
# alternative hypothesis: two-sided

# Calculate the first quartile (Q1) and third quartile (Q3) of CV values
Q1 <- quantile(cv_df$logCV, 0.25, na.rm = TRUE)
Q3 <- quantile(cv_df$logCV, 0.75, na.rm = TRUE)
IQR <- Q3 - Q1

# Define the threshold as Q1 - 1.5 * IQR
cv_threshold <- Q1 - 0.5 * IQR

cv_df %>% 
  filter(logCV > cv_threshold) %>% 
  nrow()

ggplot(cv_df, aes(x = logCV)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = cv_threshold, color = "red", linetype = "dashed") +
  labs(title = "Distribution of log Coefficient of Variation (CV)",
       subtitle = 'Threshold(houskeeping) < Q1 - 0.5 * IQR',
       x = "CV",
       y = "Frequency") +
  theme_minimal()
ggsave('3.Output/atac/1.houskeepingPeakSelect.png', width = 6, height = 5, units = 'in')

highVariableP <- cv_df %>% 
  filter(logCV > cv_threshold) %>% 
  pull(PeakID)

saveRDS(expBatRem[highVariableP,], 
        '1.Data/atac/expressionMatrixBatchEffFixed(highVariable).RDS')

write.csv(cv_df %>% 
            filter(logCV < cv_threshold) %>% 
            pull(PeakID),
          '3.Output/atac/1.housekeepingPeaks.csv')

dds <- DESeqDataSetFromMatrix(countData = expBatRem[highVariableP,], 
                              # expdat[rownames(coldat)]  expBatRem
                              colData = coldat,
                              design = ~ batch + group)



dds <- DESeq(dds)
saveRDS(dds, '1.Data/atac/1.deseqDdsFile_highVariableP.RDS')

vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c('group'), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC2, PC1, col = group)) +
  geom_point() +
  # coord_flip() +
  ggrepel::geom_text_repel(aes(label = group),
                           fontface = 2, size = 2) +
  
  scale_color_manual(values = c('#00B050', '#C00000', '#00B050', '#C00000','#00B050', '#7030A0', '#C00000')) +
  ylab(paste0("PC1: ",percentVar[1],"% variance")) +
  xlab(paste0("PC2: ",percentVar[2],"% variance")) +
  
  # coord_fixed()  +
  theme_bw() +
  theme(legend.position = 'none', 
        text = element_text(size = 12, face = 2))

ggsave('3.Output/atac/1.pcaPlot.png', 
       width = 15, height = 8, units = 'in', bg = 'white')

graph2ppt(file="3.Output/atac/1.pcaPlot.pptx", width=2.8, height=4.2) 
