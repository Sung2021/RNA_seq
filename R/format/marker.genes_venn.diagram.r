## set2 : Tcm.ctr_irr
## set3 : Tem.ctr_irr
table(set2$log2FoldChange > 0 )
table(set3$log2FoldChange > 0 )

Tcm.irr.up <- set2 %>% filter(log2FoldChange >0) %>% rownames()
Tcm.irr.dn <- set2 %>% filter(log2FoldChange <0) %>% rownames()
Tem.irr.up <- set3 %>% filter(log2FoldChange >0) %>% rownames()
Tem.irr.dn <- set3 %>% filter(log2FoldChange <0) %>% rownames()

mks <- list(Tcm.irr.up=Tcm.irr.up,
            Tcm.irr.dn=Tcm.irr.dn,
            Tem.irr.up=Tem.irr.up,
            Tem.irr.dn=Tem.irr.dn)

mks <- list(Tcm.irr.up=Tcm.irr.up,
            Tem.irr.up=Tem.irr.up)

mks <- list(Tcm.irr.dn=Tcm.irr.dn,
            Tem.irr.dn=Tem.irr.dn)


gplots::venn(mks)
venn.output <- gplots::venn(mks)
length(attributes(venn.output)$intersections)

## number of genes in venn diagram
genes <- c()
for(i in 1:length(attributes(venn.output)$intersections)){
  genes <- c(genes,c(attributes(venn.output)$intersections[[i]]))
}
genes %>% length() 
fpkm %>% filter(rownames(.) %in% genes) %>% dim() 
fpkm %>% filter(rownames(.) %in% genes) %>% 
  write.csv('RNA_seq/Iowa_Tcm Tem irradiation/fpkm.275.genes.csv')
fpkm %>% filter(rownames(.) %in% genes) %>% 
  write.csv('RNA_seq/Iowa_Tcm Tem irradiation/fpkm.292.genes.csv')


for(i in 1:length(attributes(venn.output)$intersections)){
  genes <- attributes(venn.output)$intersections[[i]]
  assign(paste0('up_',i), genes)
}
ls(pattern = 'up_')

for(i in 1:length(attributes(venn.output)$intersections)){
  genes <- attributes(venn.output)$intersections[[i]]
  assign(paste0('dn_',i), genes)
}
ls(pattern = 'dn_')

up_1
up_2
up_3
dn_1
dn_2
dn_3

fpkm.filtered %>% dim()
fpkm.filtered$group <- 'NA'
fpkm.filtered[rownames(fpkm.filtered) %in% up_1,]$group <- 'up_1'
fpkm.filtered[rownames(fpkm.filtered) %in% up_2,]$group <- 'up_2'
fpkm.filtered[rownames(fpkm.filtered) %in% up_3,]$group <- 'up_3'
fpkm.filtered[rownames(fpkm.filtered) %in% dn_1,]$group <- 'dn_1'
fpkm.filtered[rownames(fpkm.filtered) %in% dn_2,]$group <- 'dn_2'
fpkm.filtered[rownames(fpkm.filtered) %in% dn_3,]$group <- 'dn_3'
