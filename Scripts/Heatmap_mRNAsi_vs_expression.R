setwd('H:/lab data/Glioblastoma/Figures/Codes for github')
library(readxl)
library(pheatmap)
library(RColorBrewer)

#load matrix for heatmap
mat = read.csv(file.choose(),row.names = 1)

#load annotation column for heatmap
annotation_col = read.csv(file.choose(), row.names = 1)

#set annotation colors according to the annotation files
ann_colors = list(
  Marker = c(Yamanaka_factor = "#c6dc6a", 
             Cancer_sc = "#48c3cf",
             Normal_sc='#f47541')
)

#generate heatmap
pheatmap(mat,  show_colnames = T, 
         annotation_col = annotation_col,
         annotation_colors = ann_colors, 
         fontsize_row = 10, fontsize_col = 10,
         angle_col = 90,
         cluster_cols = F,
         cluster_rows = F,
         cellwidth = 10, cellheight = 10, border_color = 'black',
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))