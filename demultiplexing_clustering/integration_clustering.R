# Copyright (c) [2021] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Integration of demultiplexed snRNA samples with Seurat's CCA

library(Seurat)
library(rlist)
'%ni%' = Negate('%in%')
options(future.globals.maxSize = 30720*1024^2)
indir = '~/Dropbox/UKA/marmoset/data/'
outdir = '~/Dropbox/UKA/marmoset/'
setwd(indir)
samples = c('FP17', 'FP49', 'FP50',
			'FP51', 'FP52', 'FP53', 
			'FP54', 'FP55')

hto.samples = read.table(file = '../hto.samples.txt', 
						header = TRUE, 
						sep = '\t')
hto.samples$hto = paste0('A0', hto.samples$hto)
rownames(hto.samples) = hto.samples$hto


obj.list = c()
for (sample in samples){
	sc = readRDS(file = paste0(sample, '.rds'))
	sc$hto = sub('-.*', '', sc$hash.ID)
	sc$sample = hto.samples[sc$hto,sample]
	sc = RenameCells(object = sc, add.cell.id = sample)

	sc = NormalizeData(sc, verbose = FALSE)
	sc = FindVariableFeatures(sc, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
	obj.list = list.append(obj.list, sc)
}


res = 0.6
sc.anchors = FindIntegrationAnchors(object.list = obj.list, dims = 1:20, verbose = FALSE)
sc = IntegrateData(anchorset = sc.anchors, dims = 1:20, verbose = FALSE)
DefaultAssay(sc) = 'integrated'
sc = ScaleData(sc, verbose = FALSE)
sc = RunPCA(sc, npcs = 30, verbose = FALSE)
sc = RunUMAP(sc, reduction = 'pca', dims = 1:20, verbose = FALSE)
sc = FindNeighbors(sc, reduction = 'pca', dims = 1:20, verbose = FALSE)
sc = FindClusters(sc, resolution = res, verbose = FALSE)



#---- Re-integrate and re-cluster after removing poor quality cells


# Filter out poor quality cells
Idents(sc) = 'integrated_snn_res.0.6'
bad.clusters = c('3', '14', '21', 
				'25', '26', '27', 
				'29', '30', '31', 
				'34', '35', '36')
sc = subset(sc, subset = integrated_snn_res.0.6 %ni% bad.clusters)


# Re-integrate and -cluster cells
obj.list = SplitObject(sc, split.by = 'orig.ident')

obj.list = lapply(X = obj.list, FUN = function(x) {
    x = NormalizeData(x)
    x = FindVariableFeatures(x, selection.method = 'vst', nfeatures = 2000)
})


sc.anchors = FindIntegrationAnchors(object.list = obj.list, dims = 1:20)
sc = IntegrateData(anchorset = sc.anchors, dims = 1:20)
DefaultAssay(sc) = 'integrated'
sc = ScaleData(sc)
sc = RunPCA(sc, npcs = 30)
sc = RunUMAP(sc, reduction = 'pca', dims = 1:20)
sc = FindNeighbors(sc, reduction = 'pca', dims = 1:20)
sc = FindClusters(sc, resolution = res)



#---- Remove small doublet clusters


bad.clusters = c('19', '20', '24', 
				'25', '26', '27')
sc = subset(sc, subset = integrated_snn_res.0.6 %ni% bad.clusters)

# Save data 
saveRDS(sc, file = 'marmoset.cca.integration.filter.reclust.rds')

