library(Seurat)
library(dplyr)
sc.dir = '~/Dropbox/Marmoset/data/'
outdir = '~/Dropbox/Marmoset/demultiplexing/'
setwd(sc.dir)
sample = 'FP49'
hto.sample = 'FP42'

# Read count matrix for RNA data and hashtag oligo (HTO) data
sc.data = Read10X(data.dir = paste0('snrna/', sample))
hto.data = Read10X(data.dir = paste0('hashtag/', hto.sample, '/umi_count/'), gene.column = 1)

# Subset for the hashtags actually used
n = 5
hto.data = hto.data[rownames(hto.data)[1:n],]
colnames(hto.data) = paste0(colnames(hto.data), '-1')

# Select cell barcodes that are both detected in RNA and HTO
joint.bcs = intersect(colnames(sc.data), colnames(hto.data))

# Subset RNA and HTO counts by the joint cell barcodes
sc.data = sc.data[,joint.bcs]
hto.data = as.matrix(hto.data[,joint.bcs])

# Create Seurat object
sc.hashtag = CreateSeuratObject(counts = sc.data, project = sample, min.cells = 3, min.features = 200)

# Add HTOs as an independent assay
sc.hashtag[['HTO']] = CreateAssayObject(counts = hto.data)

# Filter out poor quality cells before demultiplexing
DefaultAssay(sc.hashtag) = 'RNA'
sc.hashtag[['percent.mt']] = PercentageFeatureSet(sc.hashtag, pattern = '^MT-')
qc_plot = VlnPlot(sc.hashtag, ncol = 2, pt.size = 0.1,
	features = c('nFeature_RNA', 'nCount_RNA', 'nCount_HTO', 'percent.mt'))
sc.hashtag = subset(sc.hashtag, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)

# Normalize, find variable features and scale data
sc.hashtag = NormalizeData(sc.hashtag, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
sc.hashtag = FindVariableFeatures(sc.hashtag, selection.method = 'mean.var.plot', verbose = FALSE)
sc.hashtag = ScaleData(sc.hashtag, verbose = FALSE, features = VariableFeatures(sc.hashtag))

# Normalize HTO data with centered log-ration (CLR) transformation
sc.hashtag = NormalizeData(sc.hashtag, assay = 'HTO', normalization.method = 'CLR', verbose = FALSE)

# Demultiplex cells based on enrichment of HTOs, k-mean clustering
# Default is k-medoid but k-means seems to work better for our data
sc.hashtag = HTODemux(sc.hashtag, assay = 'HTO', kfunc = 'kmeans', positive.quantile = 0.99)

# Save global HTO classification (negative/doublet/positive)
htos = sc.hashtag$HTO_classification.global %>% 
	data.frame() %>% 
	tibble::rownames_to_column(.)
colnames(htos) = c('cell_barcode', 'HTO.global')
htos %>% group_by(HTO.global) %>%
	summarise(count = n()) %>%
	data.frame() %>%
	write.table(., file = paste0(outdir, sample, '.HTO.classification.global.txt'), 
		sep = '\t', quote = FALSE, row.names = FALSE)

# Visualize enrichment of HTOs
Idents(sc.hashtag) = 'HTO_maxID'
p1 = RidgePlot(sc.hashtag, assay = 'HTO', features = rownames(hto.data), ncol = 1)

# UMI violin plot for singlets, doublets and negatives
Idents(sc.hashtag) = 'HTO_classification.global'
p2 = VlnPlot(sc.hashtag, features = 'nCount_RNA', pt.size = 0.1)

# Remove negatives
sc.hashtag = subset(sc.hashtag, idents = 'Negative', invert = TRUE)

# Visualize doublets and singlets in TSNE embedding
hto.dist.mtx = as.matrix(dist(t(GetAssayData(object = sc.hashtag, assay = 'HTO'))))
sc.hashtag = RunTSNE(sc.hashtag, distance.matrix = hto.dist.mtx, perplexity = 100)
p3 = DimPlot(sc.hashtag)
p4 = HTOHeatmap(sc.hashtag, assay = 'HTO', ncells = 5000, raster = FALSE)

# Extract singlets
sc.hashtag = subset(sc.hashtag, idents = 'Singlet')

# Check for batch effects by running Seurat pipeline
sc.hashtag = FindVariableFeatures(sc.hashtag, selection.method = 'mean.var.plot', verbose = FALSE)
sc.hashtag = ScaleData(sc.hashtag, features = VariableFeatures(sc.hashtag), verbose = FALSE)
sc.hashtag = RunPCA(sc.hashtag, features = VariableFeatures(sc.hashtag), verbose = FALSE)
sc.hashtag = FindNeighbors(sc.hashtag, reduction = 'pca', dims = 1:20, verbose = FALSE)
sc.hashtag = FindClusters(sc.hashtag, resolution = 0.6, verbose = FALSE)
sc.hashtag = RunUMAP(sc.hashtag, reduction = 'pca', dims = 1:20, verbose = FALSE)
p5 = DimPlot(sc.hashtag, group.by = 'HTO_classification')

# QC plot per HTO
qc_plot_hto = VlnPlot(sc.hashtag, group.by = 'HTO_classification', ncol = 2, pt.size = 0,
	features = c('nFeature_RNA', 'nCount_RNA', 'nCount_HTO', 'percent.mt'))

# Save plots 
pdf(file = paste0(outdir, sample, '_demultiplex.pdf'), height = 10, width = 8)
qc_plot
p1
p2
p3
p4
p5
qc_plot_hto
dev.off()

# Save HTO classifications in table
htos = sc.hashtag$HTO_classification %>% 
	data.frame() %>% 
	tibble::rownames_to_column(.)
colnames(htos) = c('cell_barcode', 'HTO')
htos %>% group_by(HTO) %>%
	summarise(count = n()) %>%
	data.frame() %>%
	write.table(., file = paste0(outdir, sample, '.HTO.classification.txt'), 
		sep = '\t', quote = FALSE, row.names = FALSE)

# Save singlets
saveRDS(sc.hashtag, file = paste0(outdir, sample, '.rds'))

