

## Iowa RNA-seq irradiation
## each condition gene list
gene_list <- read.csv('RNA_seq/Iowa_Tcm Tem irradiation/res_all/gene_list_each_condition.2022.05.04.csv')
gene_list %>% head()

## the column numbers to use
tmp.col <- c(2,5,4)
tmp.col <- c(2,6,4)
tmp.col <- c(1,5,3)
tmp.col <- c(1,6,3)
tmp.col <- c(2,6,5,4)
tmp.col <- c(1,6,5,3)
colnames(gene_list)[tmp.col] ## check the names of column to choose
## to prepare the list as an input for venn diagram
mks <- list()
for (i in seq_along(tmp.col)){
  tmp <- gene_list[tmp.col[i]]
  tmp <- tmp[tmp !=""]
  mks[[i]] <- tmp
}
names(mks) <- colnames(gene_list)[tmp.col]
gplots::venn(mks)
venn.output <- gplots::venn(mks)
attributes(venn.output)$intersections

## check the length of each list
attributes(venn.output)$intersections %>% str()

## save the list file to csv or table: ver1
filename = 'dn_regulated_irr.comparison.to.Tem_Tcm.'
save.date = '2022.05.04'
lapply(attributes(venn.output)$intersections, 
       function(x) write.table(data.frame(x), 
                                paste0('RNA_seq/Iowa_Tcm Tem irradiation/res_all/',
                                       filename,save.date,'.csv'), 
                                append= T, sep=',' ))
## save the list file to csv or table: ver2                         
filename = 'UP_group_comparison.to.Tem_Tcm.'
save.date = '2022.05.04'
attributes(venn.output)$intersections %>% plyr::ldply(.,cbind) %>% 
  write.csv(paste0('RNA_seq/Iowa_Tcm Tem irradiation/res_all/',
                                                           filename,save.date,'.csv'))
                                                           

