#####################################
### Iowa RNA_seq data processing ### 
setwd('~/Desktop/HMH/rds/')
library(dplyr)
library(ggplot2)
library(DESeq2)
library(edgeR)
library(reshape)

raw.counts.filtered <- read.csv('RNA_seq/Iowa_Tcm Tem irradiation/features.raw.counts.filtered.2over3.2022.03.31.csv',
                                row.names = 1)
raw.counts.filtered %>% dim()
raw.counts.filtered[1:3,]
############################################################################
############################################################################
### DESeq

count.mtx <- raw.counts.filtered
count.mtx %>% head()

## input data : info sheet
info <- data.frame(matrix(nrow = ncol(count.mtx), ncol = 3))
colnames(info) <- c('sample', 'cell_type','condition')
info$sample <- colnames(count.mtx)
info$cell_type <- substr(info$sample,1,3)
info$condition <- substr(info$sample,5,7)
info

###############################################
######## total dds : DESeq ####################
dds <- DESeqDataSetFromMatrix(count.mtx, info, ~condition)
dds <- DESeq(dds)

### vst for PCA for quick view### 
vsd <- vst(dds,blind=TRUE)
plotPCA(vsd, intgroup="sample")


# to remove the dependence of the variance on the mean
# plotPCA(vsd, intgroup="mutation")
pcaplot <- plotPCA(vsd, intgroup="sample", return=T) 
pcaplot$group <- substr(pcaplot$name,1,7)
#using the DESEQ2 plotPCA fxn we can
### drawing plot
theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), 
               axis.line = element_line(colour = "black"))
ggplot(pcaplot, aes(PC1,PC2, label=rownames(pcaplot))) + geom_point(size=3, alpha=0.2) + 
  theme + geom_text() +facet_grid(.~group)



ggplot(pcaplot, aes(PC1,PC2, shape=condition, label=rownames(pcaplot))) + geom_point(size=3, alpha=0.2) + 
  theme + xlab('PC1:58% variance') + ylab('PC2:22% variance') + geom_text(size=3,
                                                                          vjust = 0, nudge_y = 0.5, nudge_x = 0.5) +
  facet_wrap(.~group, ncol=2)
