# Copyright (c) [2021] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Summarize cell type counts for later use with scCODA

library(Seurat)
library(SeuratDisk)
library(dplyr)
indir = '~/Dropbox/UKA/marmoset/'
setwd(indir)

sc = readRDS(file = 'marmoset.cca.integration.filter.reclust.annotated.rds')

sc$sample_condition = paste0(sc$sample, '_', sc$age_group)

counts = sc@meta.data %>%
  group_by(integrated_annotations_abbrev_top) %>%
  count(sample_condition) %>%
  tidyr::spread(key = integrated_annotations_abbrev_top, 
                value = n) %>%
  as.data.frame()

write.table(counts, file = 'sample_cell_type_counts.csv',
            quote = FALSE,
            sep = ',',
            row.names = FALSE)



#---- Differential abundance testing witl standard Wilcoxon rank sum test

# Test between young and old kidney samples
cell_prop = sc@meta.data %>%
  filter(age_group %in% c('young', 'old')) %>%
  group_by(sample_condition) %>%
  dplyr::count(integrated_annotations_abbrev_top) %>% 
  group_by(sample_condition) %>% 
  mutate(prop = n/sum(n)) %>%
  as.data.frame()

cell_prop$condition = gsub('\\d+_', '', cell_prop$sample_condition)

# Wilcoxon rank sum test per cell type
cell.types = unique(cell_prop$integrated_annotations_abbrev_top)

p.values = c()
for (cell.type in cell.types){
  df = cell_prop %>%
    filter(integrated_annotations_abbrev_top == cell.type)
  
  p.values[[cell.type]] = wilcox.test(df$prop ~ df$condition)$p.value
}
p.values = unlist(p.values)
p.adjust(p.values, method = 'BH')

