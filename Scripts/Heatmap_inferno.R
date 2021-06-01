setwd('H:/lab data/Glioblastoma/Figures/Codes for github')


library(readxl)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(openxlsx)

### Figure 1B
mat = read.xlsx("1B.xlsx",rowNames = T)
mat = t(mat)
pheatmap(mat, show_colnames = T, 
         fontsize_row = 10, fontsize_col = 10,
         angle_col = 90,
         cluster_cols = F,
         cluster_rows = F,
         color = inferno(10),
         cellwidth = 10, cellheight = 10, border_color = 'grey')


### Figure 3A
mat = read.csv("3A.csv",row.names = 1)
mat[is.na(mat)] = 0
mat = t(mat)
pheatmap(mat,  show_colnames = T, 
         fontsize_row = 10, fontsize_col = 10,
         angle_col = 90,
         cluster_cols = F,
         cluster_rows = F,
         clustering_method = "complete",
         color = inferno(100),
         cellwidth = 10, cellheight = 10, border_color = 'grey')


### Supplementary figure 2K
mat = read.xlsx("Sup 2K.xlsx", rowNames = T)
mat = t(mat)
pheatmap(mat, show_colnames = T, 
         fontsize_row = 10, fontsize_col = 10,
         angle_col = 90,
         cluster_cols = F,
         cluster_rows = F,
         color = inferno(10),
         cellwidth = 10, cellheight = 10, border_color = 'grey')
