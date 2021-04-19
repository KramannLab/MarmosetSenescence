# Copyright (c) [2021] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Subclustering of fibroblasts and endothelial cells

library(Seurat)
options(future.globals.maxSize = 30720*1024^2)
indir = '~/Dropbox/UKA/marmoset/'
setwd(indir)

# Re-integrate and -cluster cells
sc = readRDS(file = 'marmoset.cca.integration.filter.reclust.annotated.rds')

cell.types = c('FIB', 'EC')
res = 0.6

for (cell.type in cell.types){
  sc.subset = subset(sc, integrated_annotations_abbrev_top %in% cell.type)
  
  # Normalize
  obj.list = SplitObject(sc.subset, split.by = 'orig.ident')
  obj.list = lapply(X = obj.list, FUN = function(x) {
    x = NormalizeData(x, verbose = FALSE)
    x = FindVariableFeatures(x, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
  })
  
  # Re-integrate
  sc.anchors = FindIntegrationAnchors(object.list = obj.list, dims = 1:20, verbose = FALSE)
  sc.subset = IntegrateData(anchorset = sc.anchors, dims = 1:20, verbose = FALSE)
  DefaultAssay(sc.subset) = 'integrated'
  sc.subset = ScaleData(sc.subset, features = rownames(sc.subset), verbose = FALSE)
  sc.subset = RunPCA(sc.subset, npcs = 30, verbose = FALSE)
  sc.subset = RunUMAP(sc.subset, reduction = 'pca', dims = 1:20, verbose = FALSE)
  sc.subset = FindNeighbors(sc.subset, reduction = 'pca', dims = 1:20, verbose = FALSE)
  sc.subset = FindClusters(sc.subset, resolution = res, verbose = FALSE)
  
  # Output
  pdf(file = paste0('integrated_', cell.type, '.pdf'))
  print(DimPlot(sc.subset))
  print(DimPlot(sc.subset, group.by = 'integrated_annotations'))
  print(DimPlot(sc.subset, group.by = 'age_group'))
  dev.off()
  
  # Save data
  saveRDS(sc.subset, file = paste0('integrated.', cell.type, '.rds'))
}

