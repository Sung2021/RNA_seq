res.threshold.0.5 %>% ggplot(aes(log2FoldChange, -log10(padj))) +
  geom_point(size=1, alpha=0.5) + theme_classic() +
  geom_vline(xintercept = c(-1,1), color='red') +
  geom_hline(yintercept = -log10(0.05), color='blue')

genes <- c('Stat4', 'Cd28','Fasl','Sell','Il2ra','Dusp2','Cd40',
           'Il12a', 'S100a10','Kdm1a', 'Jun','Id3', 'Il12rb2',
           'Cd8a', 'Klrg1', 'Cd4','Ltbr', 'Pou2f2', 'Ifitm2','Ifitm3',
           'Irf8', 'Cx3cr1','Ccr5', 'Izumo1r', 'Ldlr','Cxcr5', 'Prdm1', 'Icosl',
           'Ifng','Il9r','Kdm6b', 'Ccl4','Id2', 'Fos', 'Irf4', 'Trim13', 'Tmtc4', 'Gzmb', 
           'Myc', 'Hdac7', 'Nr4a1', 'Bcl6','Btla', 'Cd200','Cd74', 'Scd2',
           'Scd1')

res.threshold.0.5$gene <- ''
res.threshold.0.5[rownames(res.threshold.0.5) %in% genes,]$gene <- genes

res.threshold.0.5 %>% ggplot(aes(log2FoldChange, -log10(padj), 
                                 label=gene)) +
  geom_point(size=1, alpha=0.5) + theme_classic() +
  geom_vline(xintercept = c(-1,1), color='red') +
  geom_hline(yintercept = -log10(0.05), color='blue') +
  geom_text(hjust = -0.05, nudge_x = 0.05)


res.threshold.0.5 %>% ggplot(aes(log2FoldChange, -log10(padj), 
                                 label=gene)) +
  geom_point(size=1, alpha=0.2) + theme_classic() +
  geom_vline(xintercept = c(-1,1), color='red') +
  geom_hline(yintercept = -log10(0.05), color='blue') +
  geom_text(hjust = -0.05, nudge_x = 0.05) +
  scale_y_log10() + ylim(c(1,5)) + xlim(c(-5,5))