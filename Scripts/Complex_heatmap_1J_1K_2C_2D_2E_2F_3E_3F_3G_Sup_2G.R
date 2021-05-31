setwd('H:/lab data/Glioblastoma/Figures/Codes for github')

library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(digest)
library(cluster)

#Loading the annotation files
annotation_core = read.csv('Annotation_core.csv',row.names = 1)
annotation_rim = read.csv('Annotation_rim.csv',row.names = 1)
annotation_inv = read.csv('Annotation_inv.csv',row.names = 1)
annotation_neg = read.csv('Annotation_neg.csv',row.names = 1)
annotation_pos = read.csv('Annotation_pos.csv',row.names = 1)

ann_colors = list(
  Region = c(tumor_core = "#f4f4f4", 
             tumor_rim = "#f9ab8c",
             invasive_region='#ffd66b',
             neg='#72cac7',
             pos='#e26240'))


#Loading the matrix files for complex heatmap
mat = read.csv(file.choose(), row.names = 1) 

core = as.matrix(mat[,1:10])
invasive = as.matrix(mat[,11:20])
negative = as.matrix(mat[,21:30])
positive = as.matrix(mat[,31:40])
rim = as.matrix(mat[,41:50])
Breaks <- seq(min(mat), max(mat))

h1 = ComplexHeatmap::pheatmap(as.matrix(core), 
                              annotation_col = annotation_core,
                              annotation_colors = ann_colors,
                              annotation_names_row = F,
                              border_color = NA, breaks = Breaks)

h2 = ComplexHeatmap::pheatmap(as.matrix(rim),
                              annotation_col = annotation_rim,
                              annotation_colors = ann_colors,
                              annotation_names_row = F,
                              border_color = NA, breaks = Breaks)

h3 = ComplexHeatmap::pheatmap(as.matrix(invasive), 
                              annotation_col = annotation_inv,
                              annotation_colors = ann_colors,
                              annotation_names_row = F,
                              border_color = NA, breaks = Breaks)

h4 = ComplexHeatmap::pheatmap(as.matrix(negative), 
                              annotation_col = annotation_neg,
                              annotation_colors = ann_colors,
                              annotation_names_row = F,
                              border_color = NA, breaks = Breaks)

h5 = ComplexHeatmap::pheatmap(as.matrix(positive), 
                              annotation_col = annotation_pos,
                              annotation_colors = ann_colors,
                              annotation_names_col = F,
                              border_color = NA, breaks = Breaks)


#This will merge the 5 heatmaps in a single page
h1 + h2 + h3 + h4 + h5


