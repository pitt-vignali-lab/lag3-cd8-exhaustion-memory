# ----------------------------------------------------------------------------
# Analysis | Lag3 and Adam10 expression across subsets
#
# Paper    : Tumor-reactive LAG3+CD8+ T cells diverge into terminally
#            exhausted cells and long-lived memory T cells
# Original : 2.Script/LAg3Adam10Dot.R
# Input    : 1.Data/scrna/D11D19DatWithMonocleAnnot.RDS
# Output   : 3.Output/scrna/Adam10Dot.png
# Upstream : 08_fig3_partition_analysis.R
# ----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(cowplot)
library(scatterpie)
library(export)
library(presto)

dat <- readRDS('1.Data/scrna/clstedDatRes08FilteredWithMonocleAnnot.RDS')

datSub <- subset(dat, time %in% c('D11', 'D19'))

datSub$timeGeno <- paste(datSub$time, datSub$genotype, sep = '_')

Idents(datSub) <- 'timeGeno'

deG <- wilcoxauc(datSub)

deG %>%
  filter(feature %in% c('Lag3', 'Adam10')) %>%
  separate(group, into = c('Time', 'Genotype'), sep = '_') %>%
  mutate(
    Genotype = factor(
      Genotype,
      levels = c("CD44+tdT-", "CD44+tdT+", "DN", "LAG3+tdT+", "LAG3-tdT+")
    ),
    feature = factor(feature, levels = c('Lag3', 'Adam10')),
    tissue = ifelse(Genotype %in% c("CD44+tdT-", "CD44+tdT+"), 'LN', 'Tumnor')
  ) %>%
  ggplot(aes(x = Genotype, y = Time, fill = logFC, size = pct_in)) +
  geom_point(shape = 21) +
  geom_text(
    aes(
      label = case_when(
        padj > 0.05 ~ '',
        0.05 >= padj & padj > 0.01 ~ '*',
        0.01 >= padj & padj > 0.001 ~ '**',
        0.001 >= padj & padj > 0.0001 ~ '***',
        0.0001 >= padj ~ '****'
      )
    ),
    size = 2.5,
    vjust = 0
  ) +
  scale_fill_gradient2(
    high = 'red',
    low = 'blue',
    mid = 'white',
    midpoint = 0
  ) +
  facet_grid(feature ~ tissue, scales = 'free_x', ) +
  theme_bw() +
  labs(size = 'Percent expressed', fill = 'Average expression') +
  # coord_fixed() +
  theme(
    axis.text = element_text(face = 2, color = 'black'),
    text = element_text(size = 6, face = 2, color = 'black'),
    strip.background = element_blank(),
    legend.position = 'none',
    legend.box = "vertical"
  )

graph2ppt(file = paste0('3.Output/scrna/LAg3Adam10Dot.ppt'), width = 2, height = 2)

####################################################3

target = 'Lag3'

plotDf <- data.frame()


for (geno in unique(datSub$genotype)) {
  for (tim in unique(datSub$time)) {
    plotDf <- rbind(
      plotDf,
      FeaturePlot(
        subset(dat, genotype == geno & time == tim),
        feature = target
      )$data %>%
        mutate(time = tim, genotype = geno) %>%
        rename('target' = 4)
    )
  }
}

plotDf %>%
  mutate(
    genotype = factor(
      genotype,
      levels = c("CD44+tdT-", "CD44+tdT+", "DN", "LAG3+tdT+", "LAG3-tdT+")
    )
  ) %>%
  ggplot(aes(x = genotype, y = target, fill = genotype)) +
  geom_violin() +
  geom_jitter(show.legend = F, size = 0.05, alpha = 0.1) +
  scale_fill_manual(
    values = c('gray30', '#d0292f', '#54bc6e', '#8548a5', '#d0292f')
  ) +
  facet_grid(time ~ .) +
  theme_bw() +
  labs(y = target) +
  coord_fixed() +
  theme(
    legend.position = 'none',
    text = element_text(size = 15, face = 2),
    # axis.text.x = element_text(angle = 30, hjust = 1),
    axis.text = element_blank(),
    strip.background = element_blank(),
    axis.title = element_blank(),
    strip.text = element_blank()
  )


ggsave(
  file = paste0('3.Output/scrna/', target, 'Violin.png'),
  width = 5,
  height = 8
)

graph2ppt(
  file = paste0('3.Output/scrna/', target, 'Violin.ppt'),
  width = 5,
  height = 8
)


target = 'Adam10'
DotPlot(datSub, features = target, group.by = 'timeGeno')$data %>%
  separate(id, into = c('Time', 'Genotype'), sep = '_') %>%
  mutate(
    Genotype = factor(
      Genotype,
      levels = c("CD44+tdT-", "CD44+tdT+", "DN", "LAG3+tdT+", "LAG3-tdT+")
    )
  ) %>%
  ggplot(aes(x = Genotype, y = Time, fill = avg.exp.scaled, size = pct.exp)) +
  geom_point(shape = 21) +
  scale_fill_gradient2(
    high = 'red',
    low = 'blue',
    mid = 'white',
    midpoint = 0
  ) +
  theme_bw() +
  labs(size = 'Percent expressed', fill = 'Average expression') +
  coord_fixed() +
  theme(
    axis.text = element_text(size = 12, face = 2, color = 'black'),
    text = element_text(size = 15, face = 2, color = 'black'),
    legend.position = 'bottom',
    legend.box = "vertical"
  )

ggsave(file = paste0('3.Output/scrna/', target, 'Dot.png'), width = 6, height = 5)


target = 'Lag3'
deG %>%
  filter(feature %in% c(target)) %>%
  separate(group, into = c('Time', 'Genotype'), sep = '_') %>%
  mutate(
    Genotype = factor(
      Genotype,
      levels = c("CD44+tdT-", "CD44+tdT+", "DN", "LAG3+tdT+", "LAG3-tdT+")
    )
  ) %>%
  ggplot(aes(x = Genotype, y = Time, fill = logFC, size = pct_in)) +
  geom_point(shape = 21) +
  geom_text(
    aes(
      label = case_when(
        padj > 0.05 ~ '',
        0.05 >= padj & padj > 0.01 ~ '*',
        0.01 >= padj & padj > 0.001 ~ '**',
        0.001 >= padj & padj > 0.0001 ~ '***',
        0.0001 >= padj ~ '****'
      )
    ),
    size = 10,
    vjust = 0
  ) +
  scale_fill_gradient2(
    high = 'red',
    low = 'blue',
    mid = 'white',
    midpoint = 0
  ) +
  theme_bw() +
  labs(size = 'Percent expressed', fill = 'Average expression') +
  coord_fixed() +
  theme(
    axis.text = element_text(size = 12, face = 2, color = 'black'),
    text = element_text(size = 15, face = 2, color = 'black'),
    legend.position = 'bottom',
    legend.box = "vertical"
  )

ggsave(file = paste0('3.Output/scrna/', target, 'Dot.png'), width = 6, height = 5)


deG %>%
  filter(feature %in% c('Adam10')) %>%
  separate(group, into = c('Time', 'Genotype'), sep = '_') %>%
  mutate(
    Genotype = factor(
      Genotype,
      levels = c("CD44+tdT-", "CD44+tdT+", "DN", "LAG3+tdT+", "LAG3-tdT+")
    )
  ) %>%
  # filter(Time == 'D19') %>%
  ggplot(aes(x = Genotype, y = logFC, )) +
  geom_bar(stat = 'identity') +
  facet_grid(Time ~ .) +
  geom_text(
    aes(y = logFC * 1.5, label = ifelse(padj < 0.05, '*', '')),
    size = 15
  ) +
  # scale_fill_gradient2(high = 'red', low = 'blue', mid = 'white', midpoint = 0) +
  theme_bw() +
  labs(size = 'Percent expressed', fill = 'Average expression') +
  coord_fixed(5) +
  theme(
    axis.text = element_text(size = 12, face = 2, color = 'black'),
    text = element_text(size = 15, face = 2, color = 'black'),
    legend.position = 'bottom',
    legend.box = "vertical"
  )

plotDf$read <- round(plotDf$target, 1)


unique(plotDf$read)

p <- ggplot()

for (i in sort(unique(plotDf$read))) {
  p <- p +
    geom_point(
      data = plotDf %>% filter(read == i),
      aes(x = UMAP_1, y = UMAP_2, col = target),
      size = 0.5
    )
}


p +
  coord_fixed() +
  facet_grid(time ~ genotype) +
  scale_color_gradient(high = 'red', low = 'gray') +
  theme_bw()
ggsave('3.Output/scrna/Lag3Umap.png', width = 10, height = 6, units = 'in')

plotDf %>%
  group_by(time, genotype) %>%
  count(target != 0, target == 0) %>%
  mutate(prop = prop.table(n)) %>%
  # filter(time =='D11', genotype == 'CD44+tdT+') %>%
  ggplot(aes(x = '', y = prop, fill = `target != 0`)) +
  geom_bar(stat = "identity", width = 1, col = 'black') +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c('gray90', 'red')) +
  facet_grid(time ~ genotype) +
  theme_void() +
  theme(
    legend.position = 'none',
    text = element_blank(),
    panel.spacing = unit(1, "lines")
  )
graph2ppt(file = '3.Output/scrna/Lag3Pie.ppt', width = 8, height = 5)

plot_grid(p2, rel_widths = 0.01)

pUmap <- ggplotGrob(p1)
pPie <- ggplotGrob(plot_grid(p2, rel_widths = 0.1))

gt <- gtable(
  widths = unit(rep(1, 15), "null"),
  heights = unit(rep(1, 10), "null")
)

gtable_show_layout(gt)

gt <- gtable_add_grob(gt, pUmap, t = 1, b = 10, l = 1, r = 15)
gt <- gtable_add_grob(gt, pPie, t = 5, l = 2.5, b = 5, r = 2.5)
grid.draw(gt)


p +
  +coord_fixed() +
  facet_grid(time ~ genotype) +
  scale_color_gradient(high = 'red', low = 'gray') +
  theme_bw()

sp <- subset(datSub, genotype %in% c('LAG3-tdT+', 'CD44+tdT+'))

spLag3 <- FeaturePlot(
  sp,
  features = c('Lag3', 'Adam10'),
  split.by = 'time'
)$data %>%
  mutate(cat = 'sp')

dp <- subset(datSub, genotype %in% c('LAG3+tdT+'))

dpLag3 <- FeaturePlot(dp, features = 'Lag3')$data %>%
  mutate(cat = 'dp')

dn <- subset(datSub, genotype %in% c('DN', 'CD44+tdT-'))

dnLag3 <- FeaturePlot(dn, features = 'Lag3')$data %>%
  mutate(cat = 'dn')

allLag3 <- rbind(spLag3, dpLag3, dnLag3) %>%
  mutate(
    ident = factor(
      ident,
      levels = c('CD44+tdT-', 'CD44+tdT+', 'DN', 'LAG3-tdT+', 'LAG3+tdT+')
    )
  )

allLag3$read <- round(allLag3$Lag3, 1)
unique(allLag3$read)

p <- ggplot()

for (i in sort(unique(allLag3$read))) {
  p <- p +
    geom_point(
      data = allLag3 %>% filter(read == i),
      aes(x = UMAP_1, y = UMAP_2, col = Lag3),
      size = 0.5
    )
}


p +
  coord_fixed() +
  facet_grid(~ident) +
  scale_color_gradient(high = 'red', low = 'gray') +
  theme_bw()

allLag3 %>%
  group_by(ident)

VlnPlot(datSub, features = 'Lag3', group.by = 'genotype')
