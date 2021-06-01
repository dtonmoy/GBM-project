setwd("H:/lab data/Glioblastoma/Figures/Codes for github")

library(pheatmap)
library(readxl)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)

mat = read.xlsx("6A.xlsx", rowNames = T)

ComplexHeatmap::pheatmap(mat,  show_colnames = T,cluster_rows = T,
         fontsize_row = 10, fontsize_col = 10,
         clustering_distance_cols = 'spearman', 
         clustering_distance_rows = 'spearman',
         cellwidth = 10, cellheight = 10, border_color = 'black',
         color = colorRampPalette(c("#f2ea0f","#ee2829"))(100))
