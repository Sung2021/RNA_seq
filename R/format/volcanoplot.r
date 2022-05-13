## basic plot
## input data 
## DEG result in csv format
## input genes are ones with fpkm > 1 filtering
set1 <- read.csv('csv/Tle_RNA_seq.Tcm.Tem.res.10009.genes.csv', row.names = 1)
set2 <- read.csv('csv/Tle_RNA_seq.Tem.ko.wt.res.10009.genes.csv', row.names = 1)
set3 <- read.csv('csv/Tle_RNA_seq.Tcm.ko.wt.res.10009.genes.csv', row.names = 1)

## function for prepare the data
## prepare the data 
res.filtering <- function(input, fdr=0.05, fc=1.5){
  df <- input %>% 
    select(c('log2FoldChange','pvalue','padj'))
  df$color <- 'grey'
  df[df %>% filter(padj < fdr) %>% filter(log2FoldChange > log2(fc)) %>% 
       filter(padj < 0.05) %>% rownames(), ]$color <- 'red'
  df[df %>% filter(padj < fdr) %>% filter(log2FoldChange < -log2(fc)) %>% 
       filter(padj < 0.05) %>% rownames(), ]$color <- 'blue'
  df$color <- factor(df$color)
  df$size <- 0.1
  df[df$color %in% c('red','blue'),]$size <- 1
  return(df)
}
## function for drawing the plot
## drawing the plot
vlc.plot <- function(input, fdr=0.05, fc=1.5,y=c(-0.5,60), x=c(-6,6)){
  input %>% ggplot(aes(log2FoldChange, -log10(padj), 
                           color=color)) +
    geom_point(size=input$size) + theme_classic() +
    geom_vline(xintercept = c(-log2(fc)), color='blue') +
    geom_vline(xintercept = c(log2(fc)), color='red') + 
    geom_hline(yintercept = -log10(fdr), color='black') +
    scale_color_manual(values=c('blue','grey','red'))+
    ylim(y) + xlim(x)
}

## three plots 
p1 <- vlc.plot(input = res.filtering(set1))
p2 <- vlc.plot(input = res.filtering(set2))
p3 <- vlc.plot(input = res.filtering(set3))
## p1: Tem vs Tcm
## p2 : Tcm wt vs ko
## p3 : Tem wt vs ko

cowplot::plot_grid(p1,p3,p2, ncol = 3)

## number of gene in colors
res.filtering(set1)$color %>% table()
res.filtering(set2)$color %>% table()
res.filtering(set3)$color %>% table()


## volcano plot with gene names label

## input mtx: res
res <- read.csv('res_all/Tcm.ctr_irr.626.genes.2022.04.22.csv',
                row.names = 1)
## add gene name information
res$gene <- '' ## genes unwanted to plot
genes <- c('Stat4', 'Cd28','Fasl','Sell','Il2ra','Dusp2','Cd40',
           'Il12a', 'S100a10','Kdm1a', 'Jun','Id3', 'Il12rb2',
           'Cd8a', 'Klrg1', 'Cd4','Ltbr', 'Pou2f2', 'Ifitm2','Ifitm3',
           'Irf8', 'Cx3cr1','Ccr5', 'Izumo1r', 'Ldlr','Cxcr5', 'Prdm1', 'Icosl',
           'Ifng','Il9r','Kdm6b', 'Ccl4','Id2', 'Fos', 'Irf4', 'Trim13', 'Tmtc4', 'Gzmb', 
           'Myc', 'Hdac7', 'Nr4a1', 'Bcl6','Btla', 'Cd200','Cd74', 'Scd2',
           'Scd1')
res[rownames(res) %in% genes,]$gene <- genes ## genes wanted to plot
## draw the plot
res %>% ggplot(aes(log2FoldChange, -log10(padj), 
                                 label=gene)) +
  geom_point(size=1, alpha=0.5) + theme_classic() +
  geom_vline(xintercept = c(-1,1), color='red') +
  geom_hline(yintercept = -log10(0.05), color='blue') +
  geom_text(hjust = -0.05, nudge_x = 0.05) +
  ylim(c(1,5)) + xlim(c(-5,5))

## additional information
## genes to plot on the figure
genes <- c('Klrg1','Cx3cr1','Gzmb','Gzma','Fasl',
           'Cd226','Prdm1','Itga1','Ifng','Tbx21',
           'Cd8a','Nkg7','Id2','Il2ra','Themis', 
           'Ly6c2','Stat4','Cd28','Klf3','Malat1')

genes <- c('Scd1','Cd74','Id3','Sell', 'Icosl',
           'Cxcr5', 'Irf8', 'Id2')
