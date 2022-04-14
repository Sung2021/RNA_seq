############################################
########### k mean clustering ##############
dim(fpkm.above1.filtered) 
saveRDS(fpkm.above1.filtered , '2021_Qiang/2021.08.16.Qiang_RNA_seq.merged.above1.fpkm.rds')
fpkm.above1.filtered[1:3,]
fpkm.tc_wt[1:3,]
fpkm.te_wt[1:3,]
fpkm.tc_te_wt <- cbind(fpkm.tc_wt, fpkm.te_wt)
dim(fpkm.tc_te_wt)
fpkm.tc_te_wt.deg <- fpkm.tc_te_wt[rownames(set1.df),]
fpkm.tc_te_wt.deg

set.seed(20)
kClust <- kmeans(scaledata, centers=4, nstart = 1000, iter.max = 20)
kClusters <- kClust$cluster

set.seed(123)
fit2 <- kmeans(fpkm.tc_te_wt.deg, centers = 2, nstart = 10)
table(fit2$cluster)
fit2.cluster <- fit2$cluster
fit3 <- kmeans(fpkm.tc_te_wt.deg, centers = 3, nstart = 10)
table(fit3$cluster)
fit3.cluster <- fit3$cluster
fit4 <- kmeans(fpkm.tc_te_wt.deg, centers = 4, nstart = 10)
table(fit4$cluster)
fit4.cluster <- fit4$cluster
fit5 <- kmeans(fpkm.tc_te_wt.deg, centers = 5, nstart = 10)
table(fit5$cluster)
fit5.cluster <- fit5$cluster

fit.all <- cbind(fit2.cluster,
                 fit3.cluster,
                 fit4.cluster,
                 fit5.cluster)
write.csv(fit.all, '2021_Qiang/2021.08.16.Tcm_Tem_wt.fit.all.csv')
fit.all[c('Tcf7','Ctla4','Sell'),]
fit.all[fit.all[,1] == '1',]
table(fit.all[,2])
rownames(fit.all[fit.all[,2] == '1',])
fit.all[fit.all[,4] == '4',]

genes.input <- rownames(fit.all[fit.all[,2] == '1',])
fpkm.input <- scale(fpkm.tc_te_wt.deg[genes.input,])
df.anno <- data.frame(fit.all[genes.input,2])
df.anno[,1] <- factor(df.anno[,1])
colnames(df.anno) <- 'cluster'
pheatmap::pheatmap(fpkm.input, 
                   cluster_rows=F, show_rownames=T,
                   cluster_cols=FALSE, 
                   annotation_row = df.anno, 
                   fontsize_row = 7,
                   border_color = 'NA')

fit3 <- kmeans(fpkm.tc_te_wt.deg, 3)
table(fit3$cluster)
fit4 <- kmeans(fpkm.tc_te_wt.deg, 4)
table(fit4$cluster)
fit5 <- kmeans(fpkm.tc_te_wt.deg, 5)
table(fit5$cluster)

fit2$cluster

###############################################
fit2$cluster[fit2$cluster == '2']

