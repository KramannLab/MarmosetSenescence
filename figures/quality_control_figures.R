# Copyright (c) [2021] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- General QC plots

library(Seurat)
library(ggplot2)
library(viridis)
library(cowplot)
library(dplyr)
indir = '~/Dropbox/UKA/marmoset/'
outdir = '~/Dropbox/Marmoset/figures/qc_plots/'
setwd(indir)

sc = readRDS(file = 'marmoset.cca.integration.filter.reclust.annotated.rds')


# Plot cell count per sample
data = as.data.frame(table(sc$sample))
colnames(data) = c('sample', 'count')

p1 = ggplot(data, aes(x = sample, y = log10(count))) +
  geom_bar(aes(fill = sample), 
           alpha = 1, stat = 'identity') +
  scale_fill_viridis(discrete = TRUE) +
  xlab('') + 
  ylab('Number of valid cells (log10)') +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 8),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = 'none')

# Plot UMI count per sample
data = data.frame(cell = names(Idents(sc)), 
                  UMI = as.numeric(sc$nCount_RNA), 
                  gene = as.numeric(sc$nFeature_RNA),
                  sample = sc$sample,
                  mt = sc$percent.mt)

p2 = ggplot(data, aes(x = sample, y = log10(UMI))) +
  geom_boxplot(aes(color = sample),
               outlier.size = 0.5) +
  geom_violin(aes(color = sample, fill = sample), 
              alpha = 0.5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  xlab('') + 
  ylab('Number of UMIs (log10)') +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.ticks = element_blank(),
        legend.position = 'none')

# Plot gene count per sample
p3 = ggplot(data, aes(x = sample, y = log10(gene))) +
  geom_boxplot(aes(color = sample),
               outlier.size = 0.5) +
  geom_violin(aes(color = sample, fill = sample), 
              alpha = 0.5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  xlab('') + 
  ylab('Number of genes (log10)') +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.ticks = element_blank(),
        legend.position = 'none')

# Plot percentage of mitochondrial read mappeing per sample
p4 = ggplot(data, aes(x = sample, y = mt)) +
  geom_boxplot(aes(color = sample),
               outlier.size = 0.5) +
  geom_violin(aes(color = sample, fill = sample), 
              alpha = 0.5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  xlab('') + 
  ylab('Percent MT mapping') +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.ticks = element_blank(),
        legend.position = 'none')


pdf(file = paste0(outdir, 'per_sample_QC.pdf'),
    height = 8)
plot_grid(p1, p2, p3, p4, 
          ncol = 1, 
          align = 'v', 
          rel_heights = c(1, 1, 1.3))
dev.off()



data = data.frame(cell = names(Idents(sc)), 
                  UMI = as.numeric(sc$nCount_RNA), 
                  gene = as.numeric(sc$nFeature_RNA),
                  sample = sc$sample,
                  mt = sc$percent.mt)




