source("https://bioconductor.org/biocLite.R")
biocLite("Rsubread")
library(Rsubread)
library(DESeq2)
#library(DESeq)

setwd("/Volumes/JH SPLASHDR/bam_files2")

#Using DESeq2
library("DESeq2")
install.packages("gplots")
library("gplots")

#Repeating with mRNA instead of exon

counts<- featureCounts(files=c("accepted_hitsCC1.bam","accepted_hitsCC2.bam","accepted_hitsFC1.bam","accepted_hitsFC2.bam","accepted_hitsCHS1.bam","accepted_hitsCHS2.bam","accepted_hitsFHS1.bam","accepted_hitsFHS2.bam"), annot.ext = "Manduca_OGS.gff", isGTFAnnotationFile = TRUE, GTF.featureType="gene",GTF.attrType = "Name", useMetaFeatures = TRUE, allowMultiOverlap = TRUE)

write.table(x=data.frame(counts$annotation[,c("GeneID","Length")],counts$counts,stringsAsFactors=FALSE),file="counts.txt",quote=FALSE,sep="\t",row.names=FALSE)

countsTable<- read.delim("counts.txt",header=TRUE)
rownames(countsTable) <- countsTable$GeneID
countsTable <- countsTable[,-1]
countsTable <- countsTable[,-1]
head(countsTable)


countDesign <- data.frame(row.names = colnames(countsTable), reared= c(rep("Cons",2),rep("Fluc",2), rep("Cons",2),rep("Fluc",2)),shock=c(rep("NS", 4), rep("HS", 4)), cond=c(rep("CNS",2), rep("FNS",2), rep("CHS",2), rep("FHS",2)))

ddsFullCountTable<- DESeqDataSetFromMatrix(countData = countsTable, colData = countDesign, design = ~cond)

dds <- DESeq(ddsFullCountTable)

#Getting contrasts of interest
res_CNSvCHS <- results(dds, contrast=c("cond", "CNS", "CHS"))
head (res_CNSvCHS)
mcols(res)$description
resultsNames(dds)
resOrderedCNSvCHS <- res_CNSvCHS[order(res_CNSvCHS$padj),]

res_CNSvFNS <- results(dds, contrast=c("cond", "CNS", "FNS"))
head (res_CNSvFNS)
resOrderedCNSvFNS <- res_CNSvFNS[order(res_CNSvFNS$padj),]

res_FNSvFHS <- results(dds, contrast=c("cond", "FNS", "FHS"))
head (res_FNSvFHS)
resOrderedFNSvFHS <- res_FNSvFHS[order(res_FNSvFHS$padj),]

res_CHSvFHS <- results(dds, contrast=c("cond", "CHS", "FHS"))
head (res_CHSvFHS)
resOrderedCHSvFHS <- res_CHSvFHS[order(res_CHSvFHS$padj),]

plotMA(res_CNSvCHS,alpha=0.01, main="RNA-seq analysis of Manduca sexta CNSvCHS", ylim=c(-10,10))

summary(res)
sum(res_CNSvCHS$padj < 0.01, na.rm=TRUE)
sum(res_CNSvFNS$padj <0.01, na.rm=TRUE)
sum(res_FNSvFHS$padj <0.01, na.rm=TRUE)
sum(res_CHSvFHS$padj <0.01, na.rm=TRUE)

#multifactor analysis
design(ddsFullCountTable) <- formula(~reared + shock)
ddsMF <- DESeq(ddsFullCountTable)
resMF<- results(ddsMF)
head(resMF)
resMF_rear <- results(ddsMF, contrast=c("reared","Cons","Fluc"))
head(resMF_rear)
resMF_shock<- results(ddsMF, contrast=c("shock","NS","HS"))
head(resMF_shock)


plotMA(resMF_shock)
plotMA(resMF_rear, ylim=c(-3,3))
# **** Note: interactions can be included in the model
#design(ddsMF) <- formula(~ reared*shock)
#run analysis with interaction model
#ddsX <- DESeq(ddsMF)
#resultsNames(ddsX)
#resMF1 <- results(ddsX)
#head(resMF1)

biocLite("pheatmap")
library("pheatmap")
select <- order(rowMeans(counts(ddsMF,normalized=TRUE)),decreasing=TRUE)[1:20]
nt <- normTransform(ddsMF)
# defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(ddsMF)[,c("reared","shock")])

pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df)

rld <- rlog(ddsMF)
plotPCA(rld, intgroup=c("reared", "shock"))

data <- plotPCA(rld, intgroup=c("reared", "shock"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=shock, shape=reared)) +  geom_point(size=3) + theme_linedraw()+ xlab(paste0("PC1: ",percentVar[1],"% variance")) +  ylab(paste0("PC2: ",percentVar[2],"% variance"))



#Adding gene names to the Msex names
annotat<- read.csv("Msex_OGS2_Comp_Annotat_Vogel.csv")
head(annotat)

clean_annotat<- annotat %>% 
  select(Name, Seq_Description, Hit_desc.,GOs,GO_Accession,Enzyme_Codes,InterProScan)

clean_annotat$Seq_uname<- substr(clean_annotat$Name,1,11)


names(clean_annotat)
dim(clean_annotat)
names(res_CNSvCHS)
res_CNSvCHS$Seq_uname<- rownames(res_CNSvCHS)

jCNSCHS<- right_join(as.data.frame(res_CNSvCHS),clean_annotat,"Seq_uname")
names(jCNSCHS)
head(jCNSCHS)

CNS_CHS_diff<- jCNSCHS %>% 
  filter(padj < 0.01) %>% 
  select(Seq_uname,Seq_Description,padj,log2FoldChange,baseMean) %>% 
  arrange((log2FoldChange))

head(CNS_CHS_diff)
