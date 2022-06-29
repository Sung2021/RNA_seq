## draft code
## 고쳐야함
# create a matrix with 1 million counts for the
# 2nd and 3rd column, the 1st and 4th have
# half and double the counts, respectively.
m <- matrix(1e6 * rep(c(.125, .25, .25, .5), each=4),
            ncol=4, dimnames=list(1:4,1:4))
mode(m) <- "integer"
se <- SummarizedExperiment(list(counts=m), colData=DataFrame(sample=1:4))
dds <- DESeqDataSet(se, ~ 1)

# create 4 GRanges with lengths: 1, 1, 2, 2.5 Kb
gr1 <- GRanges("chr1",IRanges(1,1000)) # 1kb
gr2 <- GRanges("chr1",IRanges(c(1,1001),c( 500,1500))) # 1kb
gr3 <- GRanges("chr1",IRanges(c(1,1001),c(1000,2000))) # 2kb
gr4 <- GRanges("chr1",IRanges(c(1,1001),c(200,1300))) # 500bp
rowRanges(dds) <- GRangesList(gr1,gr2,gr3,gr4)

# the raw counts
counts(dds)

# the FPM values
fpm(dds)

# the FPKM values
fpkm(dds) ## genomic info가 있어야함


DataFrame(sample=1:ncol(count.mtx))


se <- SummarizedExperiment(as.matrix(count.mtx[genes,]), colData=DataFrame(sample=1:ncol(count.mtx)))
dds <- DESeqDataSet(se, ~ 1)
rownames(dds)


############################################################
allGenes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
allGenes.df <- allGenes %>% data.frame()
allGenes.df[1:3,]
allGenes_ranges <- data.frame(allGenes@ranges)
allGenes_ranges$ENTREZID <- allGenes_ranges$names
allGenes_ranges[1:3,]

library(clusterProfiler)
converted.genes <- clusterProfiler::bitr(allGenes_ranges$ENTREZID, fromType ="ENTREZID"  , 
                                         toType = c("SYMBOL","ENSEMBL"), 
                                         OrgDb = "org.Mm.eg.db")
converted.genes$SYMBOL %>% duplicated() %>% table()
converted.genes[-(converted.genes$SYMBOL %>% duplicated()),] %>% dim()
converted.genes <- converted.genes[-(converted.genes$SYMBOL %>% duplicated()),]
converted.genes[1:3,]

selected_genes <- allGenes_ranges[allGenes_ranges$names %in% converted.genes$ENTREZID,]
selected_genes %>% dim()
selected_genes[1:3,]
## join columns by ENTREZID
gene_info <- left_join(converted.genes, selected_genes, by='ENTREZID')
gene_info %>% dim()
gene_info[1:3,]
intersect(gene_info$names, allGenes.df$gene_id) %>% length()
gene_info %>% dim()
allGenes.df %>% dim()
merge(gene_info[,c('SYMBOL','names')], 
      allGenes.df[,c('seqnames','start','end','width','gene_id')], by.x='names',by.y='gene_id') %>% dim()
df.merge <- merge(gene_info[,c('SYMBOL','names')], 
                  allGenes.df[,c('seqnames','start','end','width','gene_id')], by.x='names',by.y='gene_id')
table(is.na(df.merge))
df.merge[1:3,]
table(duplicated(df.merge$SYMBOL))
table(duplicated(df.merge$names))

df.merge <- df.merge[!(duplicated(df.merge$SYMBOL)),]
df.merge %>% write.csv('RNA_seq/gene_genomic_info.csv')

rownames(df.merge) <- df.merge$SYMBOL
df.merge[1:3,]

table(rownames(count.mtx) %in% df.merge$SYMBOL)
genes <- rownames(count.mtx)[rownames(count.mtx) %in% df.merge$SYMBOL]
count.mtx <- count.mtx[genes,]
count.mtx.gr <- df.merge[genes,][,c(3:6,1:2)] %>% makeGRangesFromDataFrame()

se <- SummarizedExperiment(as.matrix(count.mtx[genes,]), 
                           colData=DataFrame(sample=1:ncol(count.mtx)))
dds <- DESeqDataSet(se, ~ 1)
rowRanges(dds) <- count.mtx.gr

# the FPKM values
fpkm <- fpkm(dds) ## genomic info가 있어야함
fpkm %>% write.csv('RNA_seq/Icos/Icos.72hr.RNA_seq.featurecounts.40bp.trimming.fpkm.22.06.29.csv')
fpkm %>% dim()
