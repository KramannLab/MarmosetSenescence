# Copyright (c) [2021] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- General QC plots

library(Seurat)
library(ggplot2)
library(viridis)
library(cowplot)
'%ni%' = Negate('%in%')
indir = '~/Dropbox/UKA/marmoset/'
outdir = '~/Dropbox/Marmoset/figures/qc_plots/'
setwd(indir)
source('MarmosetSenescence/sc_source/sc_source.R')

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



#---- Mitochondrial read mapping per age group

pdf(file = paste0(outdir, 'percent_mt_age_group.pdf'),
    height = 4,
    width = 5)
VlnPlot(sc, feature = 'percent.mt',
        pt.size = 0,
        group.by = 'age_group',
        cols = viridis(length(unique(sc$age_group)))) +
  theme(text = element_text(size = 8),
      axis.ticks = element_blank(),
      axis.text.x = element_text(size = 8, angle = 45),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10)) +
  labs(y = 'Percent MT mapping') +
  NoLegend() +
  ggtitle('')
dev.off()



#---- Cell cycle scoring

# Adjust gene sets for marmoset
s.genes = cc.genes$s.genes
s.genes = s.genes[s.genes %in% rownames(sc)]
# Add PRIM1, MLF1IP, RPA2, CCNE2, UBR7, RAD51, BRIP1
s.genes = c(s.genes, c('ENSCJAG00000010662',
                       'CENPU',
                       'ENSCJAG00000009473',
                       'ENSCJAG00000000011',
                       'ENSCJAG00000018283',
                       'ENSCJAG00000020821'))

g2m.genes = cc.genes$g2m.genes
g2m.genes = g2m.genes[g2m.genes %in% rownames(sc)]
# Add CKS2, TACC3,  RANGAP1
g2m.genes = c(g2m.genes, c('ENSCJAG00000007792',
                           'ENSCJAG00000001805',
                           'ENSCJAG00000004844'))

sc = CellCycleScoring(sc, s.features = s.genes,
                      g2m.features = g2m.genes)


pdf(file = paste0(outdir, 'cell_cycle_score.pdf'), width = 5, height = 5)
# S phase score per cell type
VlnPlot(sc, features = 'S.Score', 
        pt.size = 0, 
        group.by = 'integrated_annotations_abbrev', 
        cols = cell.type.abbrev.colors) +
  theme(text = element_text(size = 8),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank()) +
  labs(y = 'S cell cycle score') +
  NoLegend() +
  ggtitle('')
# S phase score per age group
VlnPlot(sc, features = 'S.Score', 
        pt.size = 0, 
        group.by = 'age_group', 
        cols = viridis(length(unique(sc$age_group)))) +
  theme(text = element_text(size = 8),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank()) +
  labs(y = 'S cell cycle score') +
  NoLegend() +
  ggtitle('')
# G2M score per cell type
VlnPlot(sc, features = 'G2M.Score', 
        pt.size = 0, 
        group.by = 'integrated_annotations_abbrev', 
        cols = cell.type.abbrev.colors) +
  theme(text = element_text(size = 8),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank()) +
  labs(y = 'G2M cell cycle score') +
  NoLegend() +
  ggtitle('')
# G2M score per age group
VlnPlot(sc, features = 'G2M.Score', 
        pt.size = 0, 
        group.by = 'age_group', 
        cols = viridis(length(unique(sc$age_group)))) +
  theme(text = element_text(size = 8),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank()) +
  labs(y = 'G2M cell cycle score') +
  NoLegend() +
  ggtitle('')
dev.off()



#---- Bar charts with distribution of cell types in age groups

data = data.frame(cell = colnames(sc), 
                  AgeGroup = sc$age_group,
                  cellType = sc$integrated_annotations_abbrev,
                  cellTypeTop = sc$integrated_annotations_abbrev_top,
                  Sample = sc$sample)


pdf(file = paste0(outdir, 'cell_type_distribution_barcharts.pdf'), width = 15)
# Cell type, stacked age
ggplot(data, aes(x = cellType, fill = as.factor(AgeGroup))) + 
  geom_bar(position = position_fill(reverse = TRUE)) + 
  labs(y = 'Proportion', 
       x = element_blank(), 
       fill = 'Age group') + 
  theme_classic() +
  scale_fill_viridis(discrete = TRUE, option = 'viridis') +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8, 
                                   angle = 90, 
                                   hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))

# Cell type top, stacked age
ggplot(data, aes(x = cellTypeTop, fill = as.factor(AgeGroup))) + 
  geom_bar(position = position_fill(reverse = TRUE)) + 
  labs(y = 'Proportion', 
       x = element_blank(), 
       fill = 'Age group') + 
  theme_classic() +
  scale_fill_viridis(discrete = TRUE, option = 'viridis') +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8, 
                                   angle = 90, 
                                   hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))

# Age group, stacked cell type
ggplot(data, aes(x = AgeGroup, fill = as.factor(cellType))) + 
  geom_bar(position = position_fill(reverse = TRUE)) + 
  labs(y = 'Proportion', 
       x = element_blank(), 
       fill = 'Cell type') + 
  theme_classic() +
  scale_fill_manual(values = cell.type.abbrev.colors) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8, 
                                   angle = 90, 
                                   hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))

# Age group, stacked cell type top
ggplot(data, aes(x = AgeGroup, fill = as.factor(cellTypeTop))) + 
  geom_bar(position = position_fill(reverse = TRUE)) + 
  labs(y = 'Proportion', 
       x = element_blank(), 
       fill = 'Cell type') + 
  theme_classic() +
  scale_fill_manual(values = cell.type.abbrev.top.colors) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8, 
                                   angle = 90, 
                                   hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))

# Per sample, stacked cell type
ggplot(data, aes(x = Sample, fill = as.factor(cellType))) + 
  geom_bar(position = position_fill(reverse = TRUE)) + 
  labs(y = 'Proportion', 
       x = element_blank(), 
       fill = 'Cell type') + 
  theme_classic() +
  scale_fill_manual(values = cell.type.abbrev.colors) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8, 
                                   angle = 90, 
                                   hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))

# Per sample, stacked cell type top
ggplot(data, aes(x = Sample, fill = as.factor(cellTypeTop))) + 
  geom_bar(position = position_fill(reverse = TRUE)) + 
  labs(y = 'Proportion', 
       x = element_blank(), 
       fill = 'Cell type') + 
  theme_classic() +
  scale_fill_manual(values = cell.type.abbrev.top.colors) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8, 
                                   angle = 90, 
                                   hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10))
dev.off()


