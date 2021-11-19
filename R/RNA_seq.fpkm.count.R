###### gene range information from TxDb database #####
###### gene range information needed for fpkm calculation
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicFeatures)

allGenes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
allGenes_ranges <- data.frame(allGenes@ranges)

## add ENTREZID column from name info
allGenes_ranges$ENTREZID <- allGenes_ranges$names

## genes in count matrix
## convert genes into Entrezid
geneset <- rownames(count.mtx) ## gene list
library(clusterProfiler)
genes_to_convert <- clusterProfiler::bitr(geneset, fromType = "SYMBOL", 
                                          toType = c("ENSEMBL","ENTREZID"), 
                                          OrgDb = "org.Mm.eg.db")
## remove duplicated genes (Symbol)
genes_to_convert <- genes_to_convert[!duplicated(genes_to_convert$SYMBOL),]
rownames(genes_to_convert) <- genes_to_convert$SYMBOL

## subset out the selected genes range information
count.mtx <- readRDS('rds/featurecounts/2021.10.15.RNA_seq.raw.count.rds')

selected_genes <- allGenes_ranges[allGenes_ranges$names %in% genes_to_convert$ENTREZID,]
selected_genes %>% dim()
## two inputs
## genes to convert : symbol, ENTREZID
## selected genes : ENTREZID, location, width 
## join columns by ENTREZID
gene_info <- left_join(genes_to_convert, selected_genes, by='ENTREZID')
gene_info %>% dim()

######### FPKM calculation ##########

fpkm <- RNAAgeCalc::count2FPKM(count.mtx,genelength = gene_info$width,idtype = "SYMBOL")
fpkm %>% head()
