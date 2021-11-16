### filter 2 over 3 in any conditions
tmp <- raw.counts !=0
tmp
f1 <- tmp[tmp[,grep('d5.wt', colnames(tmp))] %>% rowSums() >=2,] %>% rownames()
f2 <- tmp[tmp[,grep('d5.mt', colnames(tmp))] %>% rowSums() >=2,] %>% rownames()
f3 <- tmp[tmp[,grep('d9.wt', colnames(tmp))] %>% rowSums() >=2,] %>% rownames()
f4 <- tmp[tmp[,grep('d9.mt', colnames(tmp))] %>% rowSums() >=2,] %>% rownames()
union(f1,f2) %>% length()
union(union(f1,f2),f3) %>% length()
union(union(union(f1,f2),f3),f4) %>% length()

raw.counts.filtered <- raw.counts[union(union(union(f1,f2),f3),f4),]
raw.counts.filtered %>% head()
write.csv(raw.counts.filtered, 'rds/RNA_seq/Tfh_AYAA/Tfh_AYAA.raw.count.2over3.2021.11.15.csv')
