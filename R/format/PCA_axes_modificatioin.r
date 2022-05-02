#####################################
setwd('~/Desktop/HMH/')
library(dplyr)
library(ggplot2)
library(DESeq2)
library(edgeR)
library(reshape)

#################################################################
########## renormalized on 2022.01.25

########## normalized peaks all : 34271 peaks ##############
norm.peak.all <- read.csv('rds/cutandrun/re_normalized.11.15/cutandrun.2021.11.17.cleared.peak.log2fc.csv', 
                          header = T, row.names = 1)
norm.peak.all %>% head()
norm.peak.all %>% dim() ## 34271


norm.peak <- norm.peak.all[,grep('read_', colnames(norm.peak.all))]
norm.peak %>% dim()
norm.peak <- norm.peak[,c(3:ncol(norm.peak))]
norm.peak %>% colnames()
norm.peak[1:3,]

pca.tmp <- prcomp(norm.peak,scale. = T)
pca.tmp %>% str()
pca.tmp$x
pca.tmp$sdev
## And we can use sdev to compute variance explained by each Principal Component.
var_explained <- pca.tmp$sdev^2/sum(pca.tmp$sdev^2)
q3 <- 0.3231 ## 3rd Qu
pca.tmp$rotation %>% as.data.frame() %>% 
  ggplot(aes((PC1+q3),PC2, shape= substr(rownames(pca.tmp$rotation), 1,8))) + 
  geom_point(size=3) +labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
                     y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))

## sdev gives the standard deviation of principal component. 
## HHX wants to have 0 in both x and y axis
## when it doesn't have it, re-scale it to have 0 in it.

## this is the case with no 0

