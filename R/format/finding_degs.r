

############################################################


input1 <- 'Tem_ctr'
input2 <- 'Tem_irr'

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
res <- res[!(is.na(res$padj)),]
dim(res)
res.filtered <- res %>% filter(log2FoldChange >= log2(2)|log2FoldChange <= -log2(2)) %>% filter(padj < 0.05)
res.filtered %>% dim()
res.filtered %>% rownames()
table((count.input[res.filtered %>% rownames(),] %>% rowSums()) > 1)

res.filtered %>% write.csv('RNA_seq/Iowa_Tcm Tem irradiation/tem_ctr_tem_irr.deseq.genes.filtered.csv')


res.filtered %>% ggplot(aes(log2FoldChange, -log10(padj))) + geom_point()

