# ----------------------------------------------------------------------------
# Setup | Fix the final clustering resolution and cluster identities
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/1.7.clusteringResSet.R
# Input    : 1.Data/scrna/datIntCellCycleScore2.RDS
# Output   : 1.Data/scrna/clstedDatRes08.RDS
#            3.Output/scrna/Res0.8_merker_genes.csv
# Upstream : 06_scrna_clustering_resolution_test.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(presto)
library(clustree)
library(gridExtra)
library(ggdendro)
library(vegan)
library(colorspace)
library(cowplot)
library(harmony)

if(!dir.exists("3.Output/scrna/2.Clustering")){
  dir.create("3.Output/scrna/2.Clustering")
}

dat <- readRDS("1.Data/scrna/datIntegHarmony.RDS")

DimPlot(dat) +
  coord_fixed()

# clustering ---- 
dat %>% RunHarmony("batch2", theta = 2, plot_convergence = TRUE) # , epsilon.cluster=-Inf, epsilon.harmony=-Inf)

dat <- FindNeighbors(dat, reduction = "harmony", dims = 1:35, force.recalc = T)
dat <- FindClusters(dat, resolution = 0.8, group.singletons = T,
                        method = 1, graph.name = 'RNA_snn') 

dat <- RunUMAP(dat, graph = 'RNA_snn', metric = 'euclidean')

## saveRDS(dat, "1.Data/scrna/clstedDatRes08.RDS")

# marker genes ----
mkGs <- wilcoxauc(dat)

write.csv(mkGs, paste0("3.Output/scrna/2.Clustering/", Sys.Date(), "Res0.8_merker_genes.csv"))

# barplot with dendrogram (cluster dist in condition) ----

dat$anno <- paste0(dat$time, ", ", dat$genotype, ": ", dat$tissue)

t = data.frame(prop.table(table(Idents(dat), dat$anno),
                          margin = 2)) %>%  reshape::cast(Var2 ~ Var1) %>% 
  as.data.frame() %>% 
  column_to_rownames(var= "Var2")

hc=hclust(dist(t),method="ward.D2")
hc=reorder(hc,wts=-as.matrix(t)%*%seq(ncol(t))^2) # vegan::reorder.hclust
tree=ggdendro::dendro_data(as.dendrogram(hc),type="rectangle")

p1=ggplot(ggdendro::segment(tree))+
  geom_segment(aes(x=y,y=x,xend=yend,yend=xend),lineend="round",size=1)+
  scale_x_continuous(expand=expansion(add=c(0,.01)))+ # don't crop half of line between top-level nodes
  scale_y_continuous(limits=.5+c(0,nrow(t)),expand=c(0,0))+
  theme(
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    axis.ticks.length=unit(0,"pt"), # remove extra space occupied by ticks
    axis.title=element_blank(),
    panel.background=element_rect(fill="white"),
    panel.grid=element_blank(),
    plot.margin=margin(5,0,5,0)
  ) + scale_x_reverse()

t=t[hc$labels[hc$order],]
t2=data.frame(V1=rownames(t)[row(t)],V2=colnames(t)[col(t)],V3=unname(do.call(c,t)))
lab=round(100*t2$V3)
lab[lab==0]=""

p2=ggplot(t2,aes(x=factor(V1,level=rownames(t)),y=V3, fill=factor(V2, levels = 0:22)))+
  geom_bar(stat="identity",width=1,position=position_fill(reverse=T), 
           col = "gray80", alpha = 0.7)+
  # geom_text(aes(label=lab),position=position_stack(vjust=.5,reverse=T),size=3.5)+
  coord_flip()+
  scale_y_continuous(labels = scales::percent, 
                     expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0))+
  # scale_fill_manual(values=colorspace::hex(HSV(head(seq(0,360,length.out=ncol(t)+1),-1),.5,1)))+
  theme(
    axis.text=element_text(color="black",size=12, face = 2),
    legend.position="none",
    axis.title=element_blank(),
    plot.margin=margin(0,20,0,0),
    panel.background = element_blank()
  ) 

legend <- get_legend(
  # create some space to the left of the legend
  p2 + guides(fill = guide_legend(title = "Cluster",
                                  nrow = 2)) +
    theme(legend.position = "bottom")
)

png(paste0("3.Output/scrna/2.Clustering/", Sys.Date(), "_cell_distr_bar_dendro.png"),
    width = 15, height = 6, units = 'in', res = 300)
plot_grid(plot_grid(p1, p2, rel_widths = c(1, 8), align = "h"), 
          legend, ncol = 1, rel_heights = c(1, .1))
dev.off()

# barplot (genotype dist in cluster) ----

data.frame(prop.table(table(dat$anno, Idents(dat)),
                           margin = 2)) %>%  reshape::cast(Var2 ~ Var1) %>% 
  as.data.frame() %>% 
  column_to_rownames(var= "Var2") %>% 
  write.csv(paste0("3.Output/scrna/2.Clustering/", Sys.Date(), 
                   "_cellClusterDistributionPercentageTans.csv"))

ggplot(t2,aes(x = factor(V2, levels = 0:22),y = V3, fill=V1, col = V1))+
  geom_bar(stat="identity",width=0.85,position=position_fill(reverse=T), 
           alpha = 0.7) +
  scale_y_continuous(labels = scales::percent, 
                     expand=c(0,0)) +
  labs(x = "Clusters", y = "", fill = "Conditions", col = "Conditions") +
  theme_classic() +
  coord_fixed(ratio = 22) +
  theme(axis.line.y = element_line(size = 0.5),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())
  
ggsave(paste0("3.Output/scrna/2.Clustering/", Sys.Date(), "cell_distr_bar_Trans.png"),
       width = 10, height = 7, units = 'in', dpi = 300)
