# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


# Cell type coloring scheme
# Colors inspired by https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
cell.type.colors = c('Podocytes' = '#000000',
                    'Proximal tubular cells 1' = '#70036a',
                    'Proximal tubular cells 2' = '#911eb4',
                    'Proximal tubular cells 3' = '#a860a4',
                    'Proximal tubular cells 4' = '#f032e6',
                    'Proximal tubular cells 5' = '#ff8ff9',
                    'Proximal tubular cells 6'= '#e6beff',
                    'Descending thin limb' ='#ffe119',
                    'Loop of Henle' = '#ffd8b1',
                    'Distal convoluted tubule 1' = '#800000',
                    'Distal convoluted tubule 2' = '#9A6324',
                    'Connecting tubule 1' = '#cc0000',
                    'Connecting tubule 2 - Collecting duct' = '#ff6666',
                    'Intercalated A cells' = '#42d4f4',
                    'Intercalated B cells' = '#aaffc3',
                    'Endothelial cells 1' = '#000075',
                    'Endothelial cells 2' = 'blue',
                    'Fibroblasts 1' = '#ff4d04',
                    'Fibroblasts 2' = '#ffae19',
                    'Immune cells' = '#808000',
                    'Unknown 1' = '#999999',
                    'Unknown 2' = '#505050')


cell.type.abbrev.colors = c('PD' = '#000000',
                            'PT1' = '#70036a',
                            'PT2' = '#911eb4',
                            'PT3' = '#a860a4',
                            'PT4' = '#f032e6',
                            'PT5' = '#ff8ff9',
                            'PT6'= '#e6beff',
                            'DTL' ='#ffe119',                            
                            'LOH' = '#ffd8b1',
                            'DCT1' = '#800000',
                            'DCT2' = '#9A6324',
                            'CNT1' = '#cc0000',
                            'CNT2-CD' = '#ff6666',
                            'IC-A' = '#42d4f4',
                            'IC-B' = '#aaffc3',
                            'EC1' = '#000075',
                            'EC2' = 'blue',                                                        
                            'FIB1' = '#ff4d04',
                            'FIB2' = '#ffae19',
                            'IM' = '#808000',
                            'Unknown 1' = '#999999',
                            'Unknown 2' = '#505050')


cell.type.abbrev.top.colors = c('PD' = '#000000',
                            'PT' = '#70036a',
                            'DTL' ='#ffe119',                            
                            'LOH' = '#ffd8b1',
                            'DCT' = '#800000',
                            'CNT' = '#cc0000',
                            'IC' = '#42d4f4',
                            'EC' = '#000075',                                                    
                            'FIB' = '#ff4d04',
                            'IM' = '#808000',
                            'Unknown' = '#999999')




#---- Enhanced DoHeatmap function (modified)
DoHeatmap4 = function(SeuratObject, GSC, assay = 'RNA', ident, show_hr = TRUE, title = 'Expression') {
  suppressPackageStartupMessages(library(ComplexHeatmap))
  
  gg_color_hue = function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  gg_color_hue2 = function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 35, c = 100)[1:n]
  }
  
  mat = SeuratObject@assays[[assay]]@scale.data
  
  GSC = GSC[GSC[,1] %in% rownames(mat),]
  
  genes = GSC[,1]  
  genes.cols = GSC[,2]
  
  cl = as.character(ident)
  
  # Reorder
  ord = order(as.numeric(cl), decreasing = FALSE)
  mat = mat[,ord]
  cl = cl[ord]
  
  cl.cols = setNames(gg_color_hue(length(unique(cl))),unique(as.character(sort(as.numeric(cl)))))
  
  common_genes = intersect(genes,rownames(mat))
  diff_genes = setdiff(genes, rownames(mat))
  
  mat2 = rbind(mat[common_genes,],
                matrix(NA, nrow = length(diff_genes), 
                      ncol = ncol(mat),
                      dimnames = list(diff_genes, colnames(mat))))
  mat2 = mat2[genes,]
  
  hc = HeatmapAnnotation(df = data.frame('cluster' = cl),
                        col = list('cluster' = cl.cols), 
                        show_annotation_name = FALSE,
                        show_legend = FALSE)
  
  hr = rowAnnotation(df = data.frame(markers = genes.cols), 
                      show_annotation_name = FALSE,
                      show_legend = TRUE)
  
  f1 = circlize::colorRamp2(c(-5, 0, +5), c('purple', 'black', 'yellow'))

  hp = Heatmap(mat2, cluster_rows = FALSE, cluster_columns = FALSE,
              col = f1, top_annotation = hc, bottom_annotation = hc, 
              name = 'Expression',
              split = factor(genes.cols, levels = unique(genes.cols)),
              row_title_rot = 0, row_gap = unit(1.5, 'mm'),
              column_split = factor(cl, levels = unique(cl)), 
              column_title_rot = 00, column_gap = unit(1.2, 'mm'),
              show_column_names = FALSE, row_names_side = 'left', 
              column_title_gp = gpar(fontsize = 7),
              row_title_gp = gpar(fontsize = 7), 
              row_names_gp = gpar(fontsize = 2.6))

  if(show_hr) {
    hh = hp + hr 
  } else {
    hh = hp
  }
  return(hh)
}




#---- Genesorter
run_genesorter = function(sc, assay = 'RNA', slot = 'data', write.file = FALSE, out.dir = '.', file.name.prefix = NULL){
  
  # Get specificity score (specScore) and conditional probability of expression (condGeneProb)
  library(genesorteR, quietly = TRUE)
  sg = sortGenes(GetAssayData(sc, assay = assay, slot = slot), Idents(sc))
  return.list = c('specScore' = sg$specScore, 'condGeneProb' = sg$condGeneProb)

  # Write files
  if(write.file){
    if(!dir.exists(file.path(out.dir))) stop('out.dir does not exist')

    # specScore
    specScore = as.data.frame(sg$specScore)
    specScore$gene = rownames(specScore)
    specScore = specScore[,c(ncol(specScore), 1:ncol(specScore)-1)]
    write.table(specScore, 
      file = paste0(out.dir, '/', file.name.prefix, 'specScore.txt'), 
      sep = '\t', quote = FALSE, row.names = FALSE)

    # condGeneProb
    condGeneProb = as.data.frame(sg$condGeneProb)
    condGeneProb$gene = rownames(condGeneProb)
    condGeneProb = condGeneProb[,c(ncol(condGeneProb), 1:ncol(condGeneProb)-1)]
    write.table(condGeneProb, 
      file =  paste0(out.dir, '/', file.name.prefix, 'condGeneProb.txt'), 
      sep = '\t', quote = FALSE, row.names = FALSE)
  }

  return(return.list)
}




#---- PanglaoDB bar chart
panglao_db = function(matrix, species = 'hs', n = 50, out.dir = '.', file.name.prefix = NULL, return.genes = FALSE){
  library(stringr, quietly = TRUE)
  library(ggplot2, quietly = TRUE)
  '%ni%' = Negate('%in%')

  # Check input validity
  if(toupper(species) %ni% c('HS', 'MM')) stop('Only valid species input are "Hs" or "Mm"')
  if(!is.matrix(matrix) & !is(matrix, 'sparseMatrix')) stop('Input has be to a (sparse) matrix')
  if(!dir.exists(file.path(out.dir))) stop('out.dir does not exist')


  # Subset Panglao DB to species of choice
  f = read.table('~/Dropbox/UKA/data/PanglaoDB_markers_27_Mar_2020.tsv',
      sep = '\t', header = TRUE, quote = '')
  species = species
  species.subset = str_detect(f$species, regex(species, ignore_case = TRUE))
  f = f[species.subset,]


  # Get top n genes based on scores in input matrix
  genes = list()
  for (i in 1:ncol(matrix)){
    genes[[i]] = rownames(head(matrix[order(matrix[,i], decreasing = TRUE),], n))
  }


  # Get associated cell types for top genes in each cluster
  cell.types = list()
  for (i in 1:length(genes)){
    cell.types[[i]] = f[f$official.gene.symbol %in% str_to_upper(genes[[i]]),]$cell.type
  }
  names(cell.types) = colnames(matrix)


  # Bar charts of PanglaoDB top 5 cell types per cluster
  pdf(file = paste0(out.dir, '/', file.name.prefix, 'PanglaoDB_annotations.pdf'))
  for (cluster in names(cell.types)){
    # Add annotation counts and weighted counts
    cell.prop = data.frame(cell.type = cell.types[[cluster]])
    cell.prop = cell.prop %>%
      count(cell.type, name = 'count', sort = TRUE) %>%
      mutate(pct = count/nrow(cell.prop))
    cell.prop$cluster = cluster

    # Get top 5 occuring cell types
    top.ct = cell.prop %>% 
      pull(cell.type) %>%
      head(5)

    # Subset plot top 5, order cell.type in plot
    cell.prop = cell.prop[cell.prop$cell.type %in% top.ct,]
    cell.prop$cell.type = factor(cell.prop$cell.type, levels = top.ct)

    print(ggplot(cell.prop, aes(x = cluster, y = pct, fill = cell.type)) + 
      geom_bar(position = 'dodge', stat = 'identity') +
      labs(y = 'Annotation fraction', x = element_blank(), fill = 'PanglaoDB annotation') +
      theme_classic())
  }
  dev.off()
  rm(f)

  if(return.genes){
    return(genes)
  }
}




#---- Correlation heatmap with meta data
correlation_heatmap = function(object, conditionVector, assay = 'RNA', cellTypeColors = cell.type.colors){
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(viridis))

  df = Reduce('rbind', 
              AverageExpression(object,
                features = VariableFeatures(object, assay = assay), 
              assays = assay))
  df.plot = cor(df)

  # Get meta data
  condition = str_extract(colnames(df), paste(conditionVector, collapse = '|'))
  cell.types = gsub(paste0('[[:space:]]?', conditionVector, collapse = '|'), '', colnames(df))

  # Heatmap annotations
  mat.colors = colorRampPalette(rev(brewer.pal(n = 7, name = 'RdYlBu')))(100)
  condition.colors = viridis(length(conditionVector))
  names(condition.colors) = conditionVector
  ha = HeatmapAnnotation(Condition = condition,
                        Celltypes = cell.types,
                        col = list(Condition = condition.colors,
                                  Celltypes = cellTypeColors),
                        show_annotation_name = FALSE)

  # Heatmap
  p = Heatmap(matrix = df.plot,
             col = mat.colors,
             name = 'Pearson',
             show_row_names = FALSE,
             show_column_names = FALSE,
             clustering_distance_rows = 'pearson',
             clustering_distance_columns = 'pearson',
             clustering_method_columns = 'average',
             clustering_method_rows = 'average',
             cluster_rows = cluster_within_group(df.plot, cell.types),
             cluster_columns = cluster_within_group(df.plot, cell.types),
             top_annotation = ha)
  return(p)
  rm(df)
  rm(df.plot)
}




#---- Process extracellular-matrix (ECM) scores provided in Naba et al., 2016 (PMID: 26163349)
processNABA = function(filepath = '../data/NABAgsets.xls') {
  con = file(filepath, 'r')
  naba_gsets = list()
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
  split_line = unlist(strsplit(line, split = '\t'))
  naba_gsets[[split_line[1]]] = split_line[3:length(split_line)]
  }
  close(con)
  return(naba_gsets)
}




#---- Doughnut plot of senescence scores (uses abbreviated top annotations!)
doughnut_plot = function(object, score, cutOff, group, cellTypes = NA){
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(tidyr))
  suppressPackageStartupMessages(library(viridis))

  # Get meta data
  df = object@meta.data

  # Format plot legend title
  format.score = gsub('1', '', score)
  format.score = gsub('_', ' ', stringr::str_to_title(format.score))


  # Generate plots per cell typee

  # Prepare data frame per cell type
  if(all(!is.na(cellTypes))){

    prep_data = function(DF, CT){
      df.subset = DF[DF$integrated_annotations_abbrev_top == CT,]

      # Summarize counts of score > cutoff and score <= cutoff for each group
      data = df.subset %>% 
        group_by(group = UQ(rlang::sym(group))) %>%
        summarize(summary = c(sum(UQ(rlang::sym(score)) > cutOff)/n(),
          sum(UQ(rlang::sym(score)) <= cutOff)/n()),
        score = c(paste0('Score > ', cutOff), paste0('Score <= ', cutOff)),
          label = scales::percent(summary, accuracy = 0.1L),
        ymax = cumsum(summary),
          ymin = c(0, head(ymax, n = -1)),
        labelPosition = (ymax + ymin)/2) %>%
        as.data.frame()
      return(data)
    }

    p = lapply(cellTypes, FUN = function(x){
      ggplot(prep_data(df, x), aes(ymax = ymax, ymin = ymin, xmax = 4, 
        xmin = 3, fill = factor(score))) +
          facet_grid(facets =. ~ group) +
          geom_rect() +
          geom_label(x = 3.5, aes(y = labelPosition, label = label), 
                  show.legend = FALSE, size = 2) +
          scale_fill_brewer(palette = 4) +
          coord_polar(theta = 'y') +
          xlim(c(2, 4)) +
          theme_void()  +
          labs(fill = format.score) +
          ggtitle(x)
      })
  }

  # Generate overall plot

  if(all(is.na(cellTypes))){

    # Summarize counts of score > cutoff and score <= cutoff for each group
    data = df %>% 
      group_by(group = UQ(rlang::sym(group))) %>%
      summarize(summary = c(sum(UQ(rlang::sym(score)) > cutOff)/n(),
        sum(UQ(rlang::sym(score)) <= cutOff)/n()),
       score = c(paste0('Score > ', cutOff), paste0('Score <= ', cutOff)),
        label = scales::percent(summary, accuracy = 0.1L),
       ymax = cumsum(summary),
        ymin = c(0, head(ymax, n = -1)),
       labelPosition = (ymax + ymin)/2) %>%
      as.data.frame()

   p = ggplot(data, 
              aes(ymax = ymax, ymin = ymin, xmax = 4, 
                  xmin = 3, fill = factor(score))) +
        facet_grid(facets =. ~ group) +
        geom_rect() +
        geom_label(x = 3.5, aes(y = labelPosition, label = label), 
                  show.legend = FALSE, size = 2) +
        scale_fill_brewer(palette = 4) +
        coord_polar(theta = 'y') +
        xlim(c(2, 4)) +
        theme_void()  +
        labs(fill = format.score)
  }

return(p)
} 



