source("https://bioconductor.org/biocLite.R")
biocLite("Rsubread")
library(Rsubread)
#library(DESeq)

setwd("/Volumes/JH SPLASHDR/bam_files2")

counts1<- featureCounts(files=c("accepted_hitsCHS1.bam","accepted_hitsCHS2.bam","accepted_hitsFHS1.bam","accepted_hitsFHS2.bam","accepted_hitsCC1.bam","accepted_hitsCC2.bam","accepted_hitsFC1.bam","accepted_hitsFC2.bam"), annot.ext = "Manduca_OGS.gff", isGTFAnnotationFile = TRUE, GTF.attrType = "Parent", useMetaFeatures = TRUE, allowMultiOverlap = TRUE)

write.table(x=data.frame(counts1$annotation[,c("GeneID","Length")],counts1$counts,stringsAsFactors=FALSE),file="counts1.txt",quote=FALSE,sep="\t",row.names=FALSE)

countsTable1 <- read.delim("counts1.txt",header=TRUE)
rownames(countsTable1) <- countsTable1$GeneID
countsTable1 <- countsTable1[,-1]
countsTable1 <- countsTable1[,-1]
head( countsTable1)

conds<- factor(c("CHS", "CHS", "FHS", "FHS", "CNS", "CNS", "FNS", "FNS"))
countDesign <- data.frame(row.names = colnames(countsTable1), condition= conds)


cds<- newCountDataSet(countsTable,conds)
cds = estimateSizeFactors( cds )
sizeFactors( cds )
head( counts( cds, normalized=TRUE ) )
cds = estimateDispersions( cds )
str( fitInfo(cds))
plotDispEsts( cds )
res = nbinomTest( cds, "CNS", "CHS" )
head(res)
plotMA(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
resSig = res[ res$padj < 0.1, ]

#The most sig dif genes 
head(resSig [order(resSig$pval),])

#The most sig dif genes down regulated
head( resSig[ order( resSig$foldChange, -resSig$baseMean ), ] )

#The most sig dif genes up regulated
head( resSig[ order( -resSig$foldChange, -resSig$baseMean ), ] )

#Full design
cdsFull<- newCountDataSet(countsTable1, countDesign)
cdsFull<- estimateSizeFactors(cdsFull)
cdsFull<- estimateDispersions(cdsFull) 
plotDispEsts(cdsFull)

fit1 = fitNbinomGLMs( cdsFull, count ~  condition )
str(fit1)
head(fit1)

#Using DESeq2
library("DESeq2")
install.packages("gplots")
library("gplots")

ddsFullCountTable<- DESeqDataSetFromMatrix(countData = countsTable, colData = countDesign, design = ~condition)

dds <- DESeq(ddsFullCountTable)

res <- results(dds)
head (res)
mcols(res)$description




#------------  Volcano plots   --------
# standard via DESeq2
plotMA(res, main="RNA-seq analysis of Manduca sexta", ylim=c(-2,2))
plotMA(res, alpha =0.01,  main="RNA-seq analysis of Manduca sexta (strict)", ylim=c(-2,2))

#------------  Quick PCA  --------
rld <- rlog(dds)
print(plotPCA(rld, intgroup="reared"))

library(Heatplus)

# clustering takes a long time so we are limiting to first 100 genes
shortCountTable=countsTable[1:50,]
mm=data.matrix(shortCountTable, )#
heatmap_2(mm, legend=2, col=topo.colors(100))


#Repeating with mRNA instead of exon
#

library(Rsubread)
library(DESeq2)

setwd("/Volumes/JH SPLASHDR/bam_files2")


counts<- featureCounts(files=c("accepted_hitsCHS1.bam","accepted_hitsCHS2.bam","accepted_hitsFHS1.bam","accepted_hitsFHS2.bam","accepted_hitsCC1.bam","accepted_hitsCC2.bam","accepted_hitsFC1.bam","accepted_hitsFC2.bam"), annot.ext = "Manduca_OGS.gff", isGTFAnnotationFile = TRUE, GTF.featureType="gene",GTF.attrType = "Name", useMetaFeatures = TRUE, allowMultiOverlap = TRUE)

write.table(x=data.frame(counts$annotation[,c("GeneID","Length")],counts$counts,stringsAsFactors=FALSE),file="counts.txt",quote=FALSE,sep="\t",row.names=FALSE)

countsTable<- read.delim("counts.txt",header=TRUE)
rownames(countsTable) <- countsTable$GeneID
countsTable <- countsTable[,-1]
countsTable <- countsTable[,-1]
head(countsTable)


countDesign <- data.frame(row.names = colnames(countsTable), reared= c("Constant", "Constant", "Fluctuating", "Fluctuating", "Constant", "Constant", "Fluctuating", "Fluctuating"), shock=c(rep("Heat_shock", 4), rep("No_shock", 4)))

ddsFullCountTable<- DESeqDataSetFromMatrix(countData = countsTable, colData = countDesign, design = ~reared)

dds <- DESeq(ddsFullCountTable)

res <- results(dds)
head (res)
mcols(res)$description
resOrdered <- res[order(res$padj),]
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)

#multifactor analysis
ddsMF <- dds
design(ddsMF) <- formula(~ reared + shock)
ddsMF <- DESeq(ddsMF)
resMF<- results(ddsMF)
head(resMF)
design(ddsMF) <- formula(~ reared + shock)
ddsMF <- DESeq(ddsMF)
resMF <- results(ddsMF)
head(resMF)
resMFType <- results(ddsMF, contrast=c("reared","Constant","Fluctuating"))
head(resMFType)

# **** Note: interactions can be included in the model
design(ddsMF) <- formula(~reared + shock +reared:shock)
#run analysis with interaction model
ddsX <- DESeq(ddsMF)
resultsNames(ddsX)
resMF <- results(ddsMF)
head(resMF)

biocLite("pheatmap")
library("pheatmap")
select <- order(rowMeans(counts(ddsX,normalized=TRUE)),decreasing=TRUE)[1:20]
nt <- normTransform(ddsX)
# defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(ddsX)[,c("reared","shock")])

pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df)

plotPCA(rld, intgroup=c("reared", "shock"))

data <- plotPCA(rld, intgroup=c("reared", "shock"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=shock, shape=reared)) +  geom_point(size=3) + theme_linedraw()+ xlab(paste0("PC1: ",percentVar[1],"% variance")) +  ylab(paste0("PC2: ",percentVar[2],"% variance"))


# ddsXContrast1 <- results(ddsX, contrast =list("rearedFluctuating.shockNo_shock","rearedFluctuating.shockHeat_shock"))
ddsXContrast1Ordered <- res[order(ddsXContrast1$padj),]
head (ddsXContrast1Ordered)
mcols(ddsXContrast1)$description