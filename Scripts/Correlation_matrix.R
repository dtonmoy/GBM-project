setwd('H:/lab data/Glioblastoma/Figures/Codes for github')
library(psych)

##boxplot 1A
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggthemes)

log_tpm = read.csv('Sup 1A.csv')

box = log_tpm %>% select(names(log_tpm)) %>%
  pivot_longer(., cols = c(Core_Sample_27_1,
                           Core_Sample_28_1,
                           Core_Sample_29_1,
                           Core_Sample_30_1,
                           Core_Sample_31_1,
                           Core_Sample_33_2,
                           Core_Sample_34_3,
                           Core_Sample_37_2,
                           Core_Sample_38_1,
                           Core_Sample_40_2,
                           Inv_Sample_27_5,
                           Inv_Sample_28_5,
                           Inv_Sample_29_5,
                           Inv_Sample_30_5,
                           Inv_Sample_31_5,
                           Inv_Sample_33_5,
                           Inv_Sample_34_6,
                           Inv_Sample_37_5,
                           Inv_Sample_38_5,
                           Inv_Sample_40_5,
                           Neg_Sample_27_neg,
                           Neg_Sample_28_neg,
                           Neg_Sample_29_neg,
                           Neg_Sample_30_neg,
                           Neg_Sample_31_neg,
                           Neg_Sample_33_neg,
                           Neg_Sample_34_neg,
                           Neg_Sample_37_neg,
                           Neg_Sample_38_neg,
                           Neg_Sample_40_neg,
                           Pos_Sample_27_pos,
                           Pos_Sample_28_pos,
                           Pos_Sample_29_pos,
                           Pos_Sample_30_pos,
                           Pos_Sample_31_pos,
                           Pos_Sample_33_pos,
                           Pos_Sample_34_pos,
                           Pos_Sample_37_pos,
                           Pos_Sample_38_pos,
                           Pos_Sample_40_pos,
                           Rim_Sample_27_2,
                           Rim_Sample_28_2,
                           Rim_Sample_29_3,
                           Rim_Sample_30_2,
                           Rim_Sample_31_2,
                           Rim_Sample_33_3,
                           Rim_Sample_34_5,
                           Rim_Sample_37_3,
                           Rim_Sample_38_2,
                           Rim_Sample_40_3), 
               names_to = "Var", values_to = "Val")

x = box %>% mutate(Category= ifelse(grepl("Core", Var), "a_Core",
                        ifelse(grepl("Rim", Var), "b_Rim",
                               ifelse(grepl("Inv", Var), "c_Inv",
                                      ifelse(grepl("Neg", Var), "Neg",
                                             ifelse(grepl("Pos", Var), "Pos","Other"))))))

ggplot(x, aes(x = Var, y = Val, fill = Category)) +
geom_boxplot(size = .1) + theme_light()+ facet_grid(~Category, scales = "free")+
  theme(axis.text.x = element_text(angle = 90))



##correlation plot 1B
cor_data = read.csv('Sup 1B.csv')
attach(cor_data)
pairs(cor_data[2:6], pch = 19, lower.panel = NULL)

pairs.panels(cor_data[,2:6], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE, # show correlation ellipses
             col="#69b3a2"
)