## raw count 2 over 3 was previously calculated
raw.counts.filtered <- read.csv('RNA_seq/Iowa_Tcm Tem irradiation/features.raw.counts.filtered.2over3.2022.03.31.csv',
                                row.names = 1)
raw.counts.filtered %>% dim()

#################################################
## fpkm was calculated previously
fpkm <- read.csv('RNA_seq/Iowa_Tcm Tem irradiation/irr.fpkm.2022.04.14.csv', 
                 row.names = 1)
                 
c(1:3,4:5,6:8,9:10) ## columns in one group
df.fpkm <- data.frame(matrix(nrow = nrow(fpkm)))
rownames(df.fpkm) <- rownames(fpkm)
df.fpkm %>% head()
df.fpkm[,1] <- (fpkm[,c(1:3)] %>% rowMeans()) >=0.5
df.fpkm[,2] <- (fpkm[,c(4:5)] %>% rowMeans()) >=0.5
df.fpkm[,3] <- (fpkm[,c(6:8)] %>% rowMeans()) >=0.5
df.fpkm[,4] <- (fpkm[,c(9:10)] %>% rowMeans()) >=0.5
df.fpkm <-  df.fpkm[df.fpkm %>% rowSums() >= 1,]
df.fpkm %>% dim()
df.fpkm %>% head()
rownames(df.fpkm)
fpkm.filtered <- fpkm[rownames(df.fpkm),]
fpkm.filtered %>% dim()
fpkm.filtered %>% write.csv('RNA_seq/Iowa_Tcm Tem irradiation/irr.fpkm.over0.5.filtered.8721genes.2022.04.20.csv')


#################################################


### DEseq comparisons
## fpkm filtered genes

raw.counts.filtered[rownames(fpkm.filtered),] %>% is.na() %>% table()
raw.counts.filtered[rownames(fpkm.filtered),] %>% 
  write.csv('RNA_seq/Iowa_Tcm Tem irradiation/irr.rowcount.filtered.8721.2022.04.20.csv')
raw.counts.filtered <- raw.counts.filtered[,-c(grep('Tem_ctr1|Tem_irr1', colnames(raw.counts.filtered)))]
raw.counts.filtered <- raw.counts.filtered[rownames(fpkm.filtered),] ## 8721 genes
raw.counts.filtered %>% dim()
raw.counts.filtered[1:3,]
raw.counts.filtered %>% 
  write.csv('RNA_seq/Iowa_Tcm Tem irradiation/irr.rowcount.filtered.8721.2022.04.20.csv')

count.8721 <- read.csv('RNA_seq/Iowa_Tcm Tem irradiation/irr.rowcount.filtered.8721.2022.04.20.csv',
                       row.names = 1)
count.mtx <- count.8721
count.mtx
count.mtx %>% dim()
colnames(count.mtx)

###### wt tcm/wt tem #################

input1 <- 'Tcm_ctr'
input2 <- 'Tem_ctr'

grep(input1, colnames(count.mtx))
grep(input2, colnames(count.mtx))
count.input <- count.mtx[,c(grep(input1, colnames(count.mtx)),
                            grep(input2, colnames(count.mtx)))]  ## input for count 

info <- data.frame(matrix(nrow = ncol(count.input), ncol = 2))
colnames(info) <- c('sample', 'condition')
info$sample <- colnames(count.input)
info$condition <- substr(info$sample, 1,7)
info$condition <- factor(info$condition, levels = c(input1,input2))
info 
str(info)
dds <- DESeqDataSetFromMatrix(count.input, info, ~ condition)
dim(dds)

dds <- DESeq(dds)
res <- results(dds)
dim(res)
res <- data.frame(res)
is.na(res$padj) %>% table()
res.filtered <- res[!(is.na(res$padj)),]
dim(res.filtered)
res.filtered <- res.filtered %>% 
  filter(log2FoldChange >= log2(2)|log2FoldChange <= -log2(2)) %>% 
  filter(padj < 0.05)
res.filtered %>% dim()
genes.763 <- res.filtered %>% rownames()
# table((count.input[res.filtered %>% rownames(),] %>% rowSums()) > 1)
table((fpkm.filtered[res.filtered %>% rownames(),colnames(count.input)] %>% rowSums()) > 1)

fpkm.filtered[res.filtered %>% rownames(),colnames(count.input)] %>% write.csv('~/Desktop/tmp.csv')
fpkm.filtered[res.filtered %>% rownames(),colnames(count.input)] %>% colMeans()
tmp <- fpkm.filtered[res.filtered %>% rownames(),colnames(count.input)]
tmp[,grep('Tcm', colnames(tmp))] %>% rowMeans() %>% as.data.frame() %>% ggplot(aes(.)) + 
  geom_density() +scale_x_log10() + geom_vline(xintercept = 1, color='red')
tmp[,grep('Tcm', colnames(tmp))] %>% rowMeans() %>% summary()

tmp['Lad1',]
df.tmp <- data.frame(matrix(nrow = nrow(tmp)))
rownames(df.tmp) <- rownames(tmp)
threshold <- 0.5
df.tmp[,1] <- (tmp[,grep('Tcm', colnames(tmp))] %>% rowMeans()) >= threshold
df.tmp[,2] <- (tmp[,grep('Tem', colnames(tmp))] %>% rowMeans()) >= threshold
(df.tmp %>% rowSums() > 0) %>% table()
df.tmp[df.tmp %>% rowSums() > 0,] %>% dim()
df.tmp[df.tmp %>% rowSums() > 0,] %>% rownames()
res.filtered[df.tmp[df.tmp %>% rowSums() > 0,] %>% rownames(),] %>% dim()
table(res.filtered[df.tmp[df.tmp %>% rowSums() > 0,] %>% rownames(),]$log2FoldChange > 0 )

res.threshold.0.5 <- res.filtered[df.tmp[df.tmp %>% rowSums() > 0,] %>% rownames(),] 
res.threshold.0.5 %>% 
  write.csv('RNA_seq/Iowa_Tcm Tem irradiation/wt_tcm_tem.deseq.threshold.0.5.746gene.22.04.20.csv')

res.filtered %>% write.csv('RNA_seq/Iowa_Tcm Tem irradiation/wt_tcm_tem.deseq.genes.filtered.22.04.18.csv')
res.filtered <- read.csv('RNA_seq/Iowa_Tcm Tem irradiation/wt_tcm_tem.deseq.genes.filtered.22.04.18.csv',
                         row.names = 1)

res.filtered %>% ggplot(aes(log2FoldChange, -log10(padj))) + geom_point()

###
genes.1186 <- read.columns('RNA_seq/Iowa_Tcm Tem irradiation/genes.1186.csv')
res[(rownames(res) %in% genes),] %>% 
  write.csv('RNA_seq/Iowa_Tcm Tem irradiation/wt_tcm_tem.deseq.genes.1187.csv')




