setwd('H:/lab data/Glioblastoma/Figures/Codes for github')
library(readxl)
library(ggplot2)
library(ggthemes)

scat = read_excel('6D.xlsx')
attach(scat)

ggplot(scat, aes(x=Expression, y=mRNAsi)) +
  geom_point(alpha = 0.5, size = 3) +
  xlab("Log2 (TPM+1)") +
  geom_smooth(method=lm , color="Blue", fill="#c5c3c2", se=TRUE) +
  theme_bw() + facet_grid(~Category, scales = "free")
  


# ABCA13 = cor(mRNAsi[1:10],Expression[1:10], method = 'pearson')
# FHL2 = cor(mRNAsi[11:20],Expression[11:20], method = 'pearson')
# IL6 = cor(mRNAsi[21:30],Expression[21:30], method = 'pearson')
