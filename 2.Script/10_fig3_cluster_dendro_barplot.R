# ----------------------------------------------------------------------------
# Fig 3 | Cluster dendrogram and composition bar plot
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/2.3.dendroBarPlot.R
# Input    : 1.Data/scrna/D11D19DatWithMonocleAnnot.RDS
# Output   : 3.Output/scrna/D11D19Dataset/clusterDistribution_dendro.png
# Upstream : 08_fig3_partition_analysis.R
# ----------------------------------------------------------------------------

# if(!dir.exists('3.Output/scrna/D11D19Dataset')){
#   dir.create('3.Output/scrna/D11D19Dataset')
# }

# library(tidyverse)
# library(Seurat)
# library(cowplot)
# library(vegan)

# datSub <- readRDS('1.Data/scrna/D11D19DatWithMonocleAnnot.RDS')

# plotDf <- data.frame(Embeddings(datSub, reduction = 'umap'),
#                      anno = datSub$anno,
#                      time = datSub$time,
#                      genotype = datSub$genotype)
# # labDf <- plotDf %>% 
# #   group_by(anno) %>% 
# #   summarise(UMAP_1 = median(UMAP_1),
# #             UMAP_2 = median(UMAP_2))

# # ggplot(plotDf, aes(UMAP_1, UMAP_2, col = anno)) +
# #   geom_point(shape = '`', key_glyph = draw_key_point) +
# #   geom_label(data = labDf, 
# #              aes(UMAP_1, UMAP_2, label = gsub('\\:.*', '', anno), fill = anno),
# #              show.legend = F, fontface = 2, col = 'black', alpha = 0.3, size = 5) +
# #   coord_fixed() +
# #   guides(col = guide_legend(ncol = 2, 
# #                             override.aes=list(shape = 16, 
# #                                               size = 5,
# #                                               alpha = 0.6))) +
# #   ggforce::theme_no_axes()+
# #   theme(legend.position = 'bottom',
# #         legend.key = element_blank(),
# #         text = element_text(size = 10, face = 2),
# #         legend.title = element_blank(),
# #         panel.background = element_rect(fill='transparent'),
# #         plot.background = element_rect(fill='transparent', color=NA),
# #         panel.grid.major = element_blank(),
# #         panel.grid.minor = element_blank(),
# #         legend.background = element_rect(fill='transparent'))


# # ggsave('3.Output/scrna/D11D19Dataset/UMAPOverall_blockfont.png', 
# #        width = 4.5, height = 8, units = 'in', 
# #        bg = 'transparent')

# # ggplot(plotDf, aes(UMAP_1, UMAP_2, col = anno)) +
# #   geom_point(shape = '`', key_glyph = draw_key_point) +
# #   # geom_label(data = labDf, 
# #   #            aes(UMAP_1, UMAP_2, label = gsub('\\:.*', '', anno), fill = anno),
# #   #            show.legend = F, fontface = 2, col = 'white', alpha = 0.4) +
# #   facet_grid(time ~ genotype) +
# #   coord_fixed() +
# #   guides(col = guide_legend(ncol = 2, 
# #                             override.aes=list(shape = 16, 
# #                                               size = 5,
# #                                               alpha = 0.6))) +
# #   ggforce::theme_no_axes()+
# #   theme(legend.position = 'none',
# #         text = element_text(size = 10, face = 2),
# #         legend.title = element_blank(),
# #         strip.background = element_blank(),
# #         panel.background = element_rect(fill='transparent'),
# #         plot.background = element_rect(fill='transparent', color=NA),
# #         panel.grid.major = element_blank(),
# #         panel.grid.minor = element_blank(),
# #         legend.background = element_rect(fill='transparent'))

# # ggsave('3.Output/scrna/D11D19Dataset/UMAPSplit.png', 
# #        width = 8, height = 5, units = 'in')

# # cluster shifting 

# # table(datSub$anno, datSub$genotype, datSub$time) %>% 
# #   as.data.frame() %>% # filter(Var2 == 'DN', Var3 == 'D11') %>% pull(Freq) %>% sum()
# #   ggplot(aes(y = Freq, x = Var3, fill = Var1, group = Var1)) +
# #   geom_col(width = 0.5, colour = "gray20",
# #            position = position_fill(reverse=F)) +
# #   geom_area(aes(x = c("D11" = 1.25, "D19" = 1.75)[Var3]), 
# #             position = position_fill(reverse=F), 
# #             colour = "gray20", alpha = 0.3,
# #             outline.type = "full") +
# #   facet_grid(~Var2)+
# #   scale_y_continuous(labels = scales::percent, limits = ) +
# #   coord_fixed(ratio = 6)  +
# #   # scale_color_brewer(palette = 'Set2') +
# #   # scale_fill_brewer(palette = 'Set2')+
# #   labs(x = 'Time point',
# #        y = 'Cluster percentage', 
# #        fill = 'Cluster') +
# #   theme_minimal() +
# #   theme(text = element_text(size = 15, face = 2),
# #         legend.position = 'none', 
# #         panel.background = element_rect(fill='transparent'),
# #         panel.border = element_blank(),
# #         plot.background = element_rect(fill='transparent', color=NA),
# #         panel.grid.major = element_blank(),
# #         panel.grid.minor = element_blank(),
# #         legend.background = element_rect(fill='transparent'))

# # ggsave('3.Output/scrna/D11D19Dataset/BarTumorSetDistri.png',
# #        width = 10, height = 8, units = 'in')

# # cluster distribution 
# datSub$timeGeno <- paste(datSub$time, datSub$genotype)

# t = data.frame(prop.table(table(datSub$anno, datSub$timeGeno),
#                           margin = 2)) %>%  reshape::cast(Var2 ~ Var1) %>% 
#   as.data.frame() %>% 
#   column_to_rownames(var= "Var2")

# # t2 = data.frame(prop.table(table(dat$anno, dat$timeGeno),
# #                           margin = 1)) %>%  reshape::cast(Var2 ~ Var1) %>% 
# #   as.data.frame() %>% 
# #   column_to_rownames(var= "Var2")

# hc=hclust(dist(scale(t)),method="ward.D2")
# hc=reorder(hc, wts=-as.matrix(t)%*%seq(ncol(t))^2) # vegan::reorder.hclust
# tree=ggdendro::dendro_data(as.dendrogram(hc),type="rectangle")

# p1=ggplot(ggdendro::segment(tree))+
#   geom_segment(aes(x=y,y=x,xend=yend,yend=xend),lineend="round",linewidth=1)+
#   scale_x_continuous(expand=expansion(add=c(0,.01)))+ # don't crop half of line between top-level nodes
#   scale_y_continuous(limits=.5+c(0,nrow(t)),expand=c(0,0))+
#   theme(
#     axis.text=element_blank(),
#     axis.ticks=element_blank(),
#     axis.ticks.length=unit(0,"pt"), # remove extra space occupied by ticks
#     axis.title=element_blank(),
#     # panel.background=element_rect(fill="white"),
#     panel.grid=element_blank(),
#     plot.margin=margin(5,0,5,0),
#     panel.background = element_rect(fill='transparent'),
#     plot.background = element_rect(fill='transparent', color=NA),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.background = element_rect(fill='transparent')
#   ) + scale_x_reverse()

# t=t[hc$labels[hc$order],]
# t2=data.frame(V1=rownames(t)[row(t)],V2=colnames(t)[col(t)],V3=unname(do.call(c,t)))
# lab=round(100*t2$V3)
# lab[lab==0]=""

# p2=ggplot(t2,aes(x=factor(V1,level=rownames(t)),
#                  y=V3, fill=factor(V2, 
#                                    levels =  c("C1: Precursor T cells",
#                                  "C2: Stem-like memory",
#                                  "C3: Central memory",
#                                  "C4: Long-lived memory",
#                                  "C5: Biosynthetic",
#                                  "C6: IFN-induced",
#                                  "C7: Effector memmory",
#                                  "C8: Memory precursor",
#                                  "C9: Co-stimulating",
#                                  "C10: SLAMF6 progenitor-memory",
#                                  "C11: Activated",
#                                  "C12: Late activated",
#                                  "C13: Effector",
#                                  "C14: Terminally Exhausted",
#                                  "C15: Proliferative Exhausted"))))+
#   geom_bar(stat="identity",width=1,position=position_fill(reverse=T), 
#            col = "gray50", alpha = 0.7)+
#   # geom_text(aes(label=lab),position=position_stack(vjust=.5,reverse=T),size=3.5)+
#   coord_flip()+
#   scale_y_continuous(labels = scales::percent, 
#                      expand=c(0,0)) +
#   scale_x_discrete(expand=c(0,0))+
#   # scale_fill_manual(values=colorspace::hex(HSV(head(seq(0,360,length.out=ncol(t)+1),-1),.5,1)))+
#   theme(
#     axis.text=element_text(color="black",size=12, face = 2),
#     legend.position="none",
#     legend.text = element_text(face = 2),
#     axis.title=element_blank(),
#     plot.margin=margin(0,20,0,0),
#     # panel.background = element_blank(),
#     panel.background = element_rect(fill='transparent'),
#     plot.background = element_rect(fill='transparent', color=NA),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.background = element_rect(fill='transparent')
#   ) 

# legend <- get_legend(
#   # create some space to the left of the legend
#   p2 + guides(fill = guide_legend(title = "Cluster",
#                                   ncol = 8)) +
#     theme(legend.position = "bottom")
# )

# plot_grid(plot_grid(p1, p2, rel_widths = c(1, 8), align = "h"), 
#           legend, ncol = 1, rel_heights = c(1, .1))

# ggsave('3.Output/scrna/D11D19Dataset/clusterDistribution.png', 
#        width = 15, height = 5, units = 'in')

##### time point distribution ##################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(ggdendro)
library(reshape)
library(cowplot)

datSub <- readRDS("1.Data/scrna/D11D19DatWithMonocleAnnot.RDS")
datSub$timeGeno <- paste(datSub$time, datSub$genotype)
plot_cluster_bar_with_dendro <- function(datSub, target_time) {
  datSub_filtered = subset(datSub, time %in% target_time)#  & genotype != c('CD44+tdT-'))
  t = data.frame(prop.table(table(datSub_filtered$anno, datSub_filtered$timeGeno),
                            margin = 2)) %>%  reshape::cast(Var2 ~ Var1) %>% 
    as.data.frame() %>% 
    column_to_rownames(var= "Var2")
  
  # t2 = data.frame(prop.table(table(dat$anno, dat$timeGeno),
  #                           margin = 1)) %>%  reshape::cast(Var2 ~ Var1) %>% 
  #   as.data.frame() %>% 
  #   column_to_rownames(var= "Var2")
  
  hc=hclust(dist(scale(t+ 1e-5)),method="ward.D2") 
  #  "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
  hc=reorder(hc, wts=-as.matrix(t)%*%seq(ncol(t))^2) # vegan::reorder.hclust
  tree=ggdendro::dendro_data(as.dendrogram(hc),type="rectangle")
  
  p1=ggplot(ggdendro::segment(tree))+
    geom_segment(aes(x=y,y=x,xend=yend,yend=xend), lineend="round", linewidth = 0.5, color = 'black')+
    scale_x_continuous(expand=expansion(add=c(0,.01)))+ # don't crop half of line between top-level nodes
    scale_y_continuous(limits=.5+c(0,nrow(t)),expand=c(0,0))+
    theme(
      axis.text.y=element_blank(),
      axis.text.x =element_text(size=8, face = 2, color = 'white'), 
      axis.ticks=element_blank(),
      text = element_text(color = 'white'), 
      axis.ticks.length=unit(0,"pt"), # remove extra space occupied by ticks
      axis.title=element_blank(),
      # panel.background=element_rect(fill="white"),
      panel.grid=element_blank(),
      plot.margin=margin(5,0,5,0),
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill='transparent')
    ) + scale_x_reverse()
  
  t=t[hc$labels[hc$order],]
  t2=data.frame(V1=rownames(t)[row(t)],V2=colnames(t)[col(t)],V3=unname(do.call(c,t)))
  lab=round(100*t2$V3)
  lab[lab==0]=""
  
  p2=ggplot(t2,aes(x=factor(V1,level=rownames(t)),
                   y=V3, fill=factor(V2, 
                                     levels =  c("C1: Precursor T cells",
                                 "C2: Stem-like memory",
                                 "C3: Central memory",
                                 "C4: Long-lived memory",
                                 "C5: Biosynthetic",
                                 "C6: IFN-induced",
                                 "C7: Effector memmory",
                                 "C8: Memory precursor",
                                 "C9: Co-stimulating",
                                 "C10: SLAMF6 progenitor-memory",
                                 "C11: Activated",
                                 "C12: Late activated",
                                 "C13: Effector",
                                 "C14: Terminally Exhausted",
                                 "C15: Proliferative Exhausted"))))+
    geom_bar(stat="identity",width=1,position=position_fill(reverse=T), 
             col = "black", alpha = 0.7, linewidth = 0.1)+
    # geom_text(aes(label=lab),position=position_stack(vjust=.5,reverse=T),size=3.5)+
    coord_flip()+
    scale_y_continuous(labels = scales::percent, 
                       expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0))+
    # scale_fill_manual(values=colorspace::hex(HSV(head(seq(0,360,length.out=ncol(t)+1),-1),.5,1)))+
    theme(
      axis.text=element_text(size=4.5, face = 2, color = 'white'),
      text = element_text(color = 'white'),
      legend.position="none",
      legend.text = element_text(face = 2, color = 'white'),
      axis.title=element_blank(),
      plot.margin=margin(0,20,0,0),
      # panel.background = element_blank(),
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill='transparent')
    ) 
  
  legend <- get_legend(
    # create some space to the left of the legend
    p2 + guides(fill = guide_legend(title = "Cluster",
                                    ncol = 8)) +
      theme(legend.position = "bottom")
  )
  
  
  return(list(dendrogram = p1, barplot = p2, legend = legend))
}

res_d11 <- plot_cluster_bar_with_dendro(datSub, "D11")
res_d19 <- plot_cluster_bar_with_dendro(datSub, "D19")

# 组合图
library(cowplot)
plot_grid(
  plot_grid(res_d11$dendrogram, res_d11$barplot, nrow = 1, rel_widths = c(0.1, 0.9)),
  panel_spacing = unit(0.1, "lines"),
  plot_grid(res_d19$dendrogram, res_d19$barplot, nrow = 1, rel_widths = c(0.1, 0.9)),
  panel_spacing = unit(0.1, "lines"),
  # legend,
  ncol = 1,
  rel_heights = c(1, 0.05, 1, 0.01)  # adjust to control spacing
)

ggsave('3.Output/scrna/D11D19Dataset/clusterDistribution_dendro.png', 
       width = 4.3, height = 2.3, units = 'in')

export::graph2ppt(file="3.Output/scrna/D11D19Dataset/cellDistribution_dendro.pptx", 
          width = 15, height = 2) 

