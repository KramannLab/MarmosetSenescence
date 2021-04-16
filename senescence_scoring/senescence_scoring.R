# Copyright (c) [2021] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Senescence scoring with Seurat's AddModuleScore

library(Seurat)
library(ggplot2)
library(viridis)
indir = '~/Dropbox/UKA/marmoset/'
setwd(indir)
source('MarmosetSenescence/sc_source/sc_source.R')

sc = readRDS(file = 'marmoset.cca.integration.filter.reclust.annotated.rds')


# Markers from Kiss et al, 2020 PMID: 32236824 (Modified to be compatible with marmoset)
# Note these gene sets have been compiled for rodents and may not be applicable for primates

senescence.gset = list(senescence_core_genes = c('ENSCJAG00000008048', 'ENSCJAG00000036548', 
								'ENSCJAG00000039587', 'TNFRSF10B', 'CDKN1A', 'DAO'), 
					senescence_effector_genes = c('ENSCJAG00000000983', 'ENSCJAG00000018303', 
								'BRF1', 'MAP2K3', 'MAP2K6', 'SMURF2', 'TGFB1I1', 'ANGPTL2'), 
					SASP_genes = c('CCL2', 'CCL24', 'CCL5', 'CTNNB1', 'CXCL1', 'CXCL10', 
								'CXCL12', 'ENSCJAG00000061765', 'CXCL16', 'HGF', 'HMGB1', 
								'ICAM1', 'IGFBP2', 'IGFBP3', 'IGFBP4', 'IGFBP5', 'IGFBP6', 
								'IGFBP7', 'IL15', 'IL18', 'IL1A', 'IL1B', 'IL2', 'IL6', 
								'ENSCJAG00000059106', 'MMP12', 'MMP13', 'MMP14', 'PGF', 
								'PLAT', 'TIMP2', 'SERPINE1', 'ENSCJAG00000039948',
								'SERPINE2', 'HGF', 'NRG1', 'EREG', 'AREG'))

senescence.gset$senescence_overall = c(senescence.gset$senescence_core_genes,
										senescence.gset$senescence_effector_genes,
										senescence.gset$SASP_genes)


# Compute senescence scores
ctrl_genes = 35
DefaultAssay(sc) = 'RNA'
gsets = names(senescence.gset)

for (gset in gsets){
	features = list(gset = senescence.gset[[gset]])
	sc = AddModuleScore(object = sc, features = features, name = gset, ctrl = ctrl_genes)
}


# Visualise senescence scores
for (gset in gsets){
	pdf(file = paste0(outdir, 'score_', gset, '.pdf'), width = 8)

	# Feature plot
	print(FeaturePlot(sc, feature = paste0(gset, '1'), 
		label = FALSE) +
				scale_colour_gradient2(low = 'blue', 
					mid = 'lightgrey', 
					high = 'red',
					midpoint = 0, 
					limits = c(-0.5,0.5), 
					oob = scales::squish) +
				ggtitle(gset))

	# Violin plots 
	print(VlnPlot(sc, feature = paste0(gset, '1'), 
		pt.size = 0, 
		group.by = 'integrated_annotations_abbrev', 
		cols = cell.type.abbrev.colors) + 
		ggtitle(gset))

	print(VlnPlot(sc, feature = paste0(gset, '1'), 
		pt.size = 0, 
		group.by = 'age_group', 
		cols = viridis(length(unique(sc$age_group)))) +
		ggtitle(gset))

	print(VlnPlot(sc, feature = paste0(gset, '1'), 
		pt.size = 0, 
		group.by = 'integrated_annotations_abbrev_top', 
		cols = cell.type.abbrev.top.colors) +
		ggtitle(gset))

	print(VlnPlot(sc, feature = paste0(gset, '1'), 
		pt.size = 0, 
		group.by = 'integrated_annotations_abbrev_top',
		split.by = 'age_group_top', 
		cols = viridis(length(unique(sc$age_group_top))),
		split.plot = TRUE) +
		ggtitle(gset))

	print(VlnPlot(sc, feature = paste0(gset, '1'), 
		pt.size = 0, 
		group.by = 'age_group_top',
		cols = viridis(length(unique(sc$age_group_top)))) +
		ggtitle(gset))

	print(VlnPlot(sc, feature = paste0(gset, '1'), 
		pt.size = 0, 
		group.by = 'age_group',
		split.by = 'kidney_function', 
		cols = viridis(length(unique(sc$kidney_function)))) + 
		ggtitle(gset))

	dev.off()
}



#---- Visualise senescence scores in doughnut plots

cutoff = 0.1

# All cell types combined
pdf(file = paste0(outdir, 'senescent_cells_doughnut_all.pdf'))
my_doughnut_plot(object = sc, score = 'senescence_core_genes1', 
			cutOff = cutoff, 
			group = 'age_group')
my_doughnut_plot(object = sc, score = 'senescence_effector_genes1', 
			cutOff = cutoff, 
			group = 'age_group')
my_doughnut_plot(object = sc, score = 'SASP_genes1', 
			cutOff = cutoff, 
			group = 'age_group')
my_doughnut_plot(object = sc, score = 'senescence_overall1', 
			cutOff = cutoff, 
			group = 'age_group')
my_doughnut_plot(object = sc, score = 'senescence_core_genes1', 
			cutOff = cutoff, 
			group = 'age_group_top')
my_doughnut_plot(object = sc, score = 'senescence_effector_genes1', 
			cutOff = cutoff, 
			group = 'age_group_top')
my_doughnut_plot(object = sc, score = 'SASP_genes1', 
			cutOff = cutoff, 
			group = 'age_group_top')
my_doughnut_plot(object = sc, score = 'senescence_overall1', 
			cutOff = cutoff, 
			group = 'age_group_top')
dev.off()


# Per cell type
pdf(file = paste0(outdir, 'senescent_cells_doughnut_per_cell_type.pdf'), height = 15, width = 15)
p = my_doughnut_plot(object = sc, score = 'senescence_core_genes1', 
					cutOff = cutoff, 
					group = 'age_group', 
					cellTypes = names(cell.type.abbrev.top.colors))
gridExtra::grid.arrange(grobs = p, ncol = 2)

p = my_doughnut_plot(object = sc, score = 'senescence_effector_genes1',
					cutOff = cutoff, 
					group = 'age_group', 
					cellTypes = names(cell.type.abbrev.top.colors))
gridExtra::grid.arrange(grobs = p, ncol = 2)

p = my_doughnut_plot(object = sc, score = 'SASP_genes1', 
					cutOff = cutoff, 
					group = 'age_group', 
					cellTypes = names(cell.type.abbrev.top.colors))
gridExtra::grid.arrange(grobs = p, ncol = 2)

p = my_doughnut_plot(object = sc, score = 'senescence_overall1', 
					cutOff = cutoff, 
					group = 'age_group', 
					cellTypes = names(cell.type.abbrev.top.colors))
gridExtra::grid.arrange(grobs = p, ncol = 2)

p = my_doughnut_plot(object = sc, score = 'senescence_core_genes1', 
					cutOff = cutoff, 
					group = 'age_group_top', 
					cellTypes = names(cell.type.abbrev.top.colors))
gridExtra::grid.arrange(grobs = p, ncol = 2)

p = my_doughnut_plot(object = sc, score = 'senescence_effector_genes1', 
					cutOff = cutoff, 
					group = 'age_group_top', 
					cellTypes = names(cell.type.abbrev.top.colors))
gridExtra::grid.arrange(grobs = p, ncol = 2)

p = my_doughnut_plot(object = sc, score = 'SASP_genes1', 
					cutOff = cutoff, 
					group = 'age_group_top', 
					cellTypes = names(cell.type.abbrev.top.colors))
gridExtra::grid.arrange(grobs = p, ncol = 2)

p = my_doughnut_plot(object = sc, score = 'senescence_overall1', 
					cutOff = cutoff, 
					group = 'age_group_top', 
					cellTypes = names(cell.type.abbrev.top.colors))
gridExtra::grid.arrange(grobs = p, ncol = 2)
dev.off()



#---- Plot canonical senescence markers (Kiss et al, 2020 PMID: 32236824) 

senescence.markers = c('GLB1', 'ENSCJAG00000051705', 'E2F2', 'IL10', 'IL1B', 
						'ITGAM', 'ITGAX', 'LMNB1', 'PARP14', 'TNF', 'CDKN1B', 'CDKN1A')

pdf(file = paste0(outdir, 'canonical_senescence_markers.pdf'), width = 8)
for (marker in senescence.markers){
	print(VlnPlot(sc, feature = marker, 
		pt.size = 0, 
		group.by = 'integrated_annotations_abbrev', 
		cols = cell.type.abbrev.colors))

	print(VlnPlot(sc, feature = marker, 
		pt.size = 0, 
		group.by = 'integrated_annotations_abbrev_top', 
		cols = cell.type.abbrev.top.colors))

	print(VlnPlot(sc, feature = marker, 
		pt.size = 0, 
		group.by = 'age_group', 
		cols = viridis(length(unique(sc$age_group)))))
}
dev.off()


