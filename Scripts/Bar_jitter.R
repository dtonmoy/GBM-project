setwd('H:/lab data/Glioblastoma/Figures/Codes for github')

library(ggthemes)
library(ggpubr)
library(ggplot2)
library(tidyverse)
library(readxl)

jitter = read.csv('6B.csv')
attach(jitter)
names(jitter)

level_order = c('Core', 'Rim', 'Inv',	'5ALA-', '5ALA+', 'Classical',
                'Mesenchymal', 'Neural', 'Proneural')

ggplot(jitter, aes(x = factor(Condition, level = level_order),
                  y = StemnessScore, fill = Condition))+
  stat_summary(fun.y='mean', geom='bar', colour="black") +
  scale_fill_manual(values = c("#71cac7","#e46141","#d3942b","#f6f8f7",
                               "#fed66a","#93ab3c",
                               "#31b34a","#ea0d8d","#f9aa8c"))+
  geom_jitter()+
  theme_bw() + facet_grid(~Source, scales = "free", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90))+
  xlab("") + theme(legend.position = "none")
