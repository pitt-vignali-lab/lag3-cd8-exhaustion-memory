# ----------------------------------------------------------------------------
# Setup | Combine CellRanger vdj and TRUST4 TCR calls
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/1.5.overallTcrPrepare.R
# Input    : 1.Data/scrna/TCR/filtered_contig_annotations.csv
#            TRUST4 reports
# Output   : 1.Data/scrna/TCR_seqNtrust_alphabetapaired.csv
# Upstream : 01_scrna_data_init.R
# Note     : TRUST4 was added to extend repertoire coverage beyond CellRanger vdj
# ----------------------------------------------------------------------------

library(tidyverse)

# TCR seq ----
tcrSeq <- data.frame()

for (dir in list.files("1.Data/scrna/")[grep("GL_", list.files("1.Data/scrna/"))]) {
  tcrDf <- data.table::fread(paste0("1.Data/scrna/", dir, 
                                    "/TCR/filtered_contig_annotations.csv")) %>% 
    mutate(barcode = gsub('-1', '', barcode)) %>% 
    mutate(barcode = paste(dir, barcode, sep = "_")) %>% 
    select(c("barcode", "is_cell", "contig_id", "high_confidence", "length", 
             "chain", "v_gene", "d_gene", "j_gene", "c_gene", "full_length",
             "productive", "cdr3", "cdr3_nt", "reads", "umis", "raw_clonotype_id", 
             "raw_consensus_id"))
  
  tcrSeq <- rbind(tcrSeq, tcrDf)
  
  cat(dir, " is done! \n")
}

tcrSeq <- tcrSeq %>% 
  filter(productive == TRUE)



# TRUST ----
trustSeq <- data.frame()
for (dir in list.files("1.Data/scrna/")[grep("GL_", list.files("1.Data/scrna/"))]) {
  files <- list.files(paste0("1.Data/scrna/", dir, "/TCR/"))
  
  trustDf <- data.table::fread(paste0("1.Data/scrna/", dir, "/TCR/", 
                                      files[grep("TRUST", files)])) %>% 
    mutate(cid = paste(dir, gsub('\\_.*', '', cid), sep = "_"))
  
  trustSeq <- rbind(trustSeq, trustDf)
  
  cat(dir, " is done! \n")
}

trustSeq$chain <- apply(trustSeq[,c('V', 'D', 'J', 'C')],
                        1,
                        function(x) sort(unique(substr(x, 1, 3)),
                                         decreasing = T)[1])

trustSeq <- trustSeq %>% 
  filter(CDR3aa != "out_of_frame") %>% 
  rename(barcode = 'cid', reads = 'count') %>% 
  separate(V, into = c("v_gene", "v"), sep = "\\*") %>% 
  separate(D, into = c("d_gene", "d"), sep = "\\*") %>% 
  separate(J, into = c("j_gene", "v"), sep = "\\*") %>% 
  separate(C, into = c("c_gene", "v"), sep = "\\*") %>% 
  rename(cdr3 = CDR3aa) %>% 
  rename(cdr3_nt = CDR3nt) %>% 
  filter(chain %in% c('TRA', "TRB")) %>% 
  select(barcode, chain, cdr3, v_gene, d_gene, c_gene, j_gene, cdr3_nt, reads) %>% 
  unique()

# conbination 
tcrAll <- rbind(tcrSeq %>% 
                  select(barcode, chain, cdr3, v_gene, d_gene, 
                         c_gene, j_gene, cdr3_nt, reads) %>% 
                  mutate(source = 'TcrSeq'), 
                trustSeq %>% mutate(source = 'TcrTrust')) %>% 
  mutate(cdr3 = gsub('_', '', cdr3))

# beta chain
tcrBeta <- tcrAll %>% filter(chain == 'TRB')%>% 
  group_by(barcode) %>% 
  top_n(1, reads)  

tcrReddt <- tcrBeta %>% group_by(barcode) %>% summarise(n = n()) %>% filter(n > 1)
tcrBeta <- tcrBeta %>% filter(! barcode %in% tcrReddt$barcode)

tcrBeta %>% group_by(barcode) %>% summarise(n = n()) %>% pull(n) %>% table()

table(tcrBeta$source)

# TCR paired 
tcrPair <- merge(tcrAll %>% filter(chain == 'TRA') %>% 
                   group_by(barcode) %>% 
                   top_n(1, reads) %>% 
                   rename_with(~paste('alpha', .x, sep = "_")),
                 tcrAll %>% filter(chain == 'TRB')%>% 
                   group_by(barcode) %>% 
                   top_n(1, reads) %>%
                   rename_with(~paste('beta', .x, sep = "_")),
                 by.x = 'alpha_barcode', by.y = 'beta_barcode') %>% 
  unique() %>% 
  mutate(pair = paste(alpha_cdr3, beta_cdr3, sep = '_'))

tcrPairReddt <- tcrPair %>% group_by(alpha_barcode) %>% 
  summarise(n = n()) %>% filter(n > 1)

tcrPair <- tcrPair %>% filter(! alpha_barcode %in% tcrPairReddt$alpha_barcode)


# Venn
x <- list(
  paired = tcrPair %>% select(alpha_barcode, pair) %>% 
    unique() %>% pull(alpha_barcode),
  beta = tcrBeta %>% select(barcode, cdr3) %>% 
    unique() %>% pull(barcode)
)

ggvenn::ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)

write.csv(tcrPair, "1.Data/scrna/TCR_seqNtrust_alphabetapaired.csv")

write.csv(tcrBeta, '1.Data/scrna/TCR_seqNtrust_betaOnly.csv')

