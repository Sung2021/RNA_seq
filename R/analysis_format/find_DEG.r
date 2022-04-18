######   #################
input1 <- 'Tcm_irr'
input2 <- 'Tem_irr'

count.input <- count.mtx[,c(grep(input1, colnames(count.mtx)),
                            grep(input2, colnames(count.mtx)))]  ## input for count 
info <- data.frame(matrix(nrow = ncol(count.input), ncol = 2))
colnames(info) <- c('sample', 'condition')
info$sample <- colnames(count.input)
info$condition <- substr(info$sample, 1,7)
info$condition <- factor(info$condition, levels = c(input1,input2))
info 
dds <- DESeqDataSetFromMatrix(count.input, info, ~ condition)
dim(dds)
dds <- DESeq(dds)
res <- results(dds)
dim(res)
res <- data.frame(res)
res %>% write.csv('RNA_seq/Iowa_Tcm Tem irradiation/dds.res.7246.genes.tcm_irr.tem_irr.csv')

## filter DEG based on several criteria
is.na(res$padj) %>% table()
res.filtered <- res[!(is.na(res$padj)),]
dim(res.filtered)
res.filtered <- res.filtered %>% 
  filter(log2FoldChange >= log2(2)|log2FoldChange <= -log2(2)) %>% 
  filter(padj < 0.05)
res.filtered %>% dim()
res.filtered %>% rownames()

## volcano plot
res.filtered %>% ggplot(aes(log2FoldChange, -log10(padj))) + geom_point()+
  geom_vline(xintercept = c(-log2(2), log2(2)), color='red') +
  geom_hline(yintercept = 0.0, color='blue')

