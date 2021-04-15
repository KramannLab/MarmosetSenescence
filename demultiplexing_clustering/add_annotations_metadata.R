# Copyright (c) [2021] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Add cell type annotations

library(Seurat)
indir = '~/Dropbox/UKA/marmoset/'
setwd(indir)

sc = readRDS(file = 'marmoset.cca.integration.filter.reclust.rds')
Idents(sc) = 'integrated_snn_res.0.6'


# Cell type annotations
sc = RenameIdents(sc, `0` = 'Proximal tubular cells 1',
				`1`= 'Proximal tubular cells 2',
				`2` = 'Proximal tubular cells 3', 
				`3` = 'Proximal tubular cells 4',
				`4` = 'Proximal tubular cells 5', 
				`5` = 'Loop of Henle',
				`6` = 'Descending thin limb', 
				`7` = 'Proximal tubular cells 6',
				`8` = 'Distal convoluted tubule 1', 
				`9` = 'Fibroblasts 1',
				`10` = 'Distal convoluted tubule 2', 
				`11` = 'Connecting tubule 1',
				`12` = 'Endothelial cells 1', 
				`13` = 'Endothelial cells 2',
				`14` = 'Intercalated A cells', 
				`15` = 'Connecting tubule 2 - Collecting duct',
				`16` = 'Intercalated B cells', 
				`17` = 'Immune cells', 
				`18` = 'Podocytes', 
				`21` = 'Unknown 1',
				`22` = 'Fibroblasts 2', 
				`23` = 'Unknown 2')

sc$integrated_annotations = Idents(sc)


# Abbreviated cell type annotations
sc$integrated_annotations_abbrev = plyr::mapvalues(Idents(sc), 
				from = c('Proximal tubular cells 1', 
						'Proximal tubular cells 2',
						'Proximal tubular cells 3', 
						'Proximal tubular cells 4',
						'Proximal tubular cells 5', 
						'Loop of Henle',
						'Descending thin limb', 
						'Proximal tubular cells 6',
						'Distal convoluted tubule 1', 
						'Fibroblasts 1',
						'Distal convoluted tubule 2', 
						'Connecting tubule 1',
						'Endothelial cells 1', 
						'Endothelial cells 2',
						'Intercalated A cells', 
						'Connecting tubule 2 - Collecting duct',
						'Intercalated B cells', 
						'Immune cells',
						'Podocytes', 
						'Unknown 1',
						'Fibroblasts 2', 
						'Unknown 2'),
				to = c('PT1', 
						'PT2',
						'PT3', 
						'PT4',
						'PT5', 
						'LOH',
						'DTL', 
						'PT6',
						'DCT1', 
						'FIB1',
						'DCT2', 
						'CNT1',
						'EC1', 
						'EC2',
						'IC-A', 
						'CNT2-CD',
						'IC-B', 
						'IM',
						'PD', 
						'Unknown 1',
						'FIB2', 
						'Unknown 2'))


# Abbreviated high level cell type annotations
sc$integrated_annotations_abbrev_top = plyr::mapvalues(sc$integrated_annotations_abbrev, 
				from = c('PT1', 'PT2',
						'PT3', 'PT4',
						'PT5', 'LOH',
						'DTL', 'PT6',
						'DCT1', 'FIB1',
						'DCT2', 'CNT1',
						'EC1', 'EC2',
						'IC-A', 'CNT2-CD',
						'IC-B', 'IM',
						'PD', 'Unknown 1',
						'FIB2', 'Unknown 2'),
				to = c('PT', 'PT',
						'PT', 'PT',
						'PT', 'LOH',
						'DTL', 'PT',
						'DCT', 'FIB',
						'DCT', 'CNT',
						'EC', 'EC',
						'IC', 'CNT',
						'IC', 'IM',
						'PD', 'Unknown',
						'FIB', 'Unknown'))



#---- Add meta data

meta.data = read.table(file = '~/Dropbox/Marmoset/Metadata/marmoset.metadata.txt', 
						sep = '\t', 
						header = TRUE)
rownames(meta.data) = meta.data[,1]
meta.data[,1] = NULL
samples = as.character(sc$sample)
meta.table = meta.data[samples,]
rownames(meta.table) = rownames(sc@meta.data)
sc = AddMetaData(sc, metadata = meta.table)


# Add age groups
sc$age_group = 'young'
sc$age_group[sc$age > 3] = 'young_adult'
sc$age_group[sc$age > 7] = 'adult'
sc$age_group[sc$age > 11] = 'old'
sc$age_group_top = 'young'
sc$age_group_top[sc$age > 7] = 'old'


sc$cell_type_age = paste(sc$integrated_annotations, sc$age_group)
sc$cell_type_age_top = paste(sc$integrated_annotations, sc$age_group_top)



#---- Order annotations

sc$integrated_annotations_abbrev = factor(sc$integrated_annotations_abbrev, 
									levels = names(cell.type.abbrev.colors))

sc$integrated_annotations_abbrev_top = factor(sc$integrated_annotations_abbrev_top, 
									levels = names(cell.type.abbrev.top.colors))

sc$age_group = factor(sc$age_group,
				levels = c('young', 'young_adult', 'adult', 'old'))

sc$age_group_top = factor(sc$age_group_top,
				levels = c('young', 'old'))

sc$kidney_function = factor(sc$kidney_function,
					levels = c('good', 'bad'))


# Save data
saveRDS(sc, file = '../marmoset.cca.integration.filter.reclust.annotated.rds')

