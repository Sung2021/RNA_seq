#####################################
### Iowa RNA_seq data processing ### 
setwd('~/Desktop/HMH/rds/')
library(dplyr)
library(ggplot2)
library(DESeq2)
library(edgeR)
library(reshape)
### raw data were generated by featurecounts
### read raw counts and save in rds format
list.files('RNA_seq_features/Iowa', full.names = T)
feature_list <- list.files('RNA_seq_features/Iowa', full.names = T, 
                           pattern = 'matrix.txt')
feature_list
length(feature_list)
tmp.int <- read.table(feature_list[1],
                      header = T, row.names = 1)
colnames(tmp.int) <- paste0(strsplit(feature_list[1], split = '/')[[1]][3] %>% substr(1,2),
                            '.count')
tmp.df <- tmp.int
for(i in 2:length(feature_list)){
  sample <- strsplit(feature_list[i], split = '/')[[1]][3] %>% substr(1,2)
  tmp <- read.table(feature_list[i],
                    header = T, col.names = c('gene',paste0(sample,'count')))
  tmp.df <- cbind(tmp.df, tmp %>% select(contains('count')))
}

tmp.df %>% head()
## reorder the columns to 1:12
tmp.df <- tmp.df[,c(1,5:12,2:4)]
## rename the columns
colnames(tmp.df) <- paste0(c(rep('Tcm',3),rep('Tem',3),rep('Tcm',3),rep('Tem',3)),'_',
                        c(rep('ctr',6), rep('irr',6)),
                        rep(seq(1:3),4))
tmp.df %>% head()

tmp.df[c('Tcf7','Ctla4'),]
## save the raw count to csv
write.csv(tmp.df, 'RNA_seq/Iowa_Tcm Tem irradiation/features.raw.counts.2022.03.31.csv')


###############################################################################
##### initial filtering ######
### filtering 

raw.counts <- tmp.df
raw.counts %>% head(30)
raw.counts[c('Tcf7','Icos','Ctla4'),]

raw.counts %>% dim()

## no expression gene removed
table(rowSums(raw.counts) == 0)
raw.counts <- raw.counts[rowSums(raw.counts) != 0,]
raw.counts %>% dim()

summary(rowSums(raw.counts))
boxplot(rowSums(raw.counts)) ## outliers : Ftl1, Tmsb4, Rn45s
hist((rowSums(raw.counts)), breaks = 1000)

raw.counts[rowSums(raw.counts) >= 3e6,]
which(rowSums(raw.counts) >= 3e6)
raw.counts <- raw.counts[-which(rowSums(raw.counts) >= 3e6)[[1]],]
raw.counts %>% dim()

write.csv(raw.counts, 'RNA_seq/Iowa_Tcm Tem irradiation/features.raw.counts.filtered.2022.03.31.csv')

###############################################################################
### filter 2 over 3 in any conditions
tmp <- raw.counts !=0
tmp %>% colnames()

## genes 2 over 3 exp in each condition 
f1 <- tmp[tmp[,c(1:3)] %>% rowSums() >=2,] %>% rownames()
f2 <- tmp[tmp[,c(4:6)] %>% rowSums() >=2,] %>% rownames()
f3 <- tmp[tmp[,c(7:9)] %>% rowSums() >=2,] %>% rownames()
f4 <- tmp[tmp[,c(10:12)] %>% rowSums() >=2,] %>% rownames()

union(union(union(f1,f2),f3),f4) %>% length()

raw.counts.filtered <- raw.counts[union(union(union(f1,f2),f3),f4),]
raw.counts.filtered %>% head()
raw.counts.filtered %>% dim()

raw.counts.filtered %>% 
  write.csv('RNA_seq/Iowa_Tcm Tem irradiation/features.raw.counts.filtered.2over3.2022.03.31.csv')

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
pcaplot <- plotPCA(vsd, intgroup="condition", return=T) 
pcaplot$group <- substr(pcaplot$name,1,7)
#using the DESEQ2 plotPCA fxn we can
### drawing plot
theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), 
               axis.line = element_line(colour = "black"))
ggplot(pcaplot, aes(PC1,PC2, shape=condition, label=rownames(pcaplot))) + geom_point(size=3, alpha=0.2) + 
  theme + xlab('PC1:58% variance') + ylab('PC2:22% variance') + geom_text(size=3,
                                                                          vjust = 0, nudge_y = 0.5, nudge_x = 1)

ggplot(pcaplot, aes(PC1,PC2, shape=condition, label=rownames(pcaplot))) + geom_point(size=3, alpha=0.2) + 
  theme + xlab('PC1:58% variance') + ylab('PC2:22% variance') + geom_text(size=3,
                                                                          vjust = 0, nudge_y = 0.5, nudge_x = 0.5) +
  facet_wrap(.~group, ncol=2)