###### to add two colums for Tcm_enr, Tem_enr 
B <- set1[set1$log2FoldChange > 0,] %>% rownames()
A <- set1[set1$log2FoldChange < 0,] %>% rownames()
input.data <- fpkm.filtered
tmp.column <- data.frame(matrix(nrow = nrow(input.data),
                                ncol = 2))
rownames(tmp.column) <- rownames(input.data)
colnames(tmp.column) <- c('A','B')

###### in the pheatmap, Tcm_enr : red (+1), Tem_enr : blue (-1), not in either (0)
tmp.column$A <- 0
tmp.column$B <- 0
### using the ATAC_peak information
tmp.column[rownames(tmp.column) %in% A, 'A'] <- 1 ## Tcm_enr
tmp.column[rownames(tmp.column) %in% B, 'B'] <- -1 ## Tem_enr
colnames(tmp.column) <- c('Tcm_enr','Tem_enr')
pheatmap::pheatmap(tmp.column, cluster_cols = F, cluster_rows = F, 
                   show_rownames = F,
                   color=colorRampPalette(c("blue", "white", "red"))(50))
