# Copyright (c) [2021] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Convert Seurat object to AnnData file (scanpy, .h5ad) for scCODA analysis

# Source
#https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html

library(SeuratDisk)
indir = '~/Dropbox/UKA/marmoset/'
setwd(indir)

# Save Seurat object as h5Seurat file
sc = readRDS(file = 'marmoset.cca.integration.filter.reclust.annotated.rds')
filename = 'marmoset.integrated.h5Seurat'
SaveH5Seurat(sc, filename = filename)

# Convert h5Seurat file to AnnData file 
Convert(filename, dest = 'h5ad')

