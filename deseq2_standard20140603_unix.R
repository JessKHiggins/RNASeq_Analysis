######### Standard Pipeline DEseq 2 ###############
# Modified from Love MI, Huber W and Anders S (2014).
# “Moderated estimation of fold change and dispersion for RNA-Seq data with DESeq2.” bioRxiv. http://dx.doi.org/10.1101/002832, http://dx.doi.org/10.1101/002832.
# Tested on R v 3.1
# cdjones 2014
#
# to download and install DEseq2
# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
####################################################

#------------ Load Library ------------
# load deseq lib
library("DESeq2")
#--------------------------------------

######### code snippets ###############
# custom R code for functions

#######################################

#----------Identify data files ----------
# note must be raw counts for the NB to work properly

# path to data matrix
# *** note: that for multiple comparisons it is important that the "base" be the first series of data ***

# set working dir
setwd=("/netscr/cdjones/RNAseq/")  ##  CHANGE cdjones to your ONYEN

# assumes these are in the working directory
datamatrixfilename="results/flydev_matrix.txt"
columndatafilename="results/sampleinfo.csv"

# make a count table
CountTable = read.table( datamatrixfilename, header=FALSE, row.names=1 )
head(CountTable)

# load sample information, which sets up expt'l design
sampleInfo <- read.csv( columndatafilename )
head(sampleInfo)
sampleInfo <- DataFrame( sampleInfo ) #works

# make sure the rownames of the sample info have the sample data
rownames( sampleInfo ) <- sampleInfo$sample
head(sampleInfo)

## set up column names if TransCnt_v02.py and maketables_0.3.rb were used.
colnames(CountTable)<-c(rownames(sampleInfo))
head(CountTable)


#----------- build the data object to do the analysis ------------
# remind ourselves of our exp'l design
head(sampleInfo)

# Make dds object
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = CountTable,  # the RNAseq counts
  colData = sampleInfo,   # the sample and treatment
  design = ~stage)  # the statistical model used, similar to ANOVA style models in R for  Kij = NB({mu}ij; {alpha}i)
                                # **** note:  that the order matters... the last item is by default the contrast althought this can be altered.
                                # later changes alt form: design(ddsFullCountTable) <- formula(~stage + genotype)
ddsFullCountTable # prints the summary of the dds


#----------- Run the model and analysis ------------
# This single function executes several functions in DEseq2 that were separate calls in earlier version
# 1) estimating size factors == "normalize counts"
# 2) estimating dispersions == look across whole data set to estimate a function that estimates {alpha}i  
# 3) gene-wise dispersion estimates == now assign {alpha}geneX specifically to geneX 
# 4) mean-dispersion relationship == describe the relationship between {mu}ij and {alpha}i, recalibrate estimates
# 5) final dispersion estimates  == final {alpha}i, assign to genes. We now have fit the data to the NB framework
# 6) fitting model and testing  == use Wald test in an ANOVA-like framework
# **** note:  The ANOVA-like framework means that you can have multiple factors and levels and test simultaneously

#  load data into the DESeq pipeline and do analysis
dds <- DESeq(ddsFullCountTable)

# ouput in human readable format
res <- results(dds)
head (res)
mcols(res)$description
# *** note:  you can change multiple testing p.adjust  "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none")


# -- some simple filtering and ordering --

# cutoff by fold change
resFoldCutoff <-res[res$log2FoldChange>1,]
head(resFoldCutoff)

# cutoff by multiple testing adjusted p-value
resFoldCutoff <-res[which(res$padj<0.05),]
head(resFoldCutoff)

#sort by p-value
resOrdered <- res[order(res$padj),]
head(resOrdered)

# ouput the data to a file XLS can import
write.table(resOrdered,"export_result_table.txt", col.names=NA, sep="\t")



# -- additional data analysis and comparisons --

# *****note: to try these use the alldata.txt and p_sampleinfo.csv files and rerun

# some data sets may have multiple levels for a factor (eg, gneotype: wt, het, mut)
# using contrasts to fit a coefficient for each
# make sure out model is right
design(ddsFullCountTable) <- formula(~stage + genotype)
#run analysis
dds <- DESeq(ddsFullCountTable)

# pull out contrast of interest
resContrast1 <- results(dds, contrast=c("genotype","het","mut"))
resContrast1Ordered <- res[order(resContrast1$padj),]
head (resContrast1Ordered)
mcols(resContrast1)$description

resContrast2 <- results(dds, contrast=c("genotype","wt","het"))
resContrast2Ordered <- res[order(resContrast2$padj),]
head (resContrast2Ordered)
mcols(resContrast2)$description


# **** Note: interactions can be included in the model
design(ddsFullCountTable) <- formula(~stage + genotype +stage:genotype)
#run analysis with interaction model
ddsX <- DESeq(ddsFullCountTable)
resultsNames(ddsX)

# note that stageA.genotype.wt is the interaction of stage A and genotype wt

# pull out contrast of interest
ddsXContrast1 <- results(ddsX, contrast =list("stageA.genotype.wt","stageL1.genotypemut"))
ddsXContrast1Ordered <- res[order(ddsXContrast1$padj),]
head (ddsXContrast1Ordered)
mcols(ddsXContrast1)$description

################################################
#      -- Graphics and viewing the data  --    #

# some standard representations of RNASeq data
# these will use the dds and res objects from 
# sections above (but not the more exotic 
# comparisons in the later half).

# For some of these the alldata.txt and p_sampleinfo.csv may be a better example.

# if gplots is not installed:   install.packages( pkgs= "gplots" )

#------------ Load Additional Libraries --------
# load 
library("RColorBrewer")
library("gplots")
#--------------------------------------

#------------  Volcano plots   --------
# standard via DESeq2
plotMA(res, main="RNA-seq analysis of Fly Development (loose)", ylim=c(-2,2))

plotMA(res, alpha =0.01,  main="RNA-seq analysis of Fly Development (strict)", ylim=c(-2,2))

#------------  Heat Map   --------
hmcol <- colorRampPalette(brewer.pal(12, "GnBu"))(100)  # GnBu

#  look at raw counts
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]  ## subset
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))

# look at how fit to NB cleans up the signal from the noise in the data
vsd <- varianceStabilizingTransformation(dds)  
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))

#-------- with clustering using a different heatmap function--------\
# if you need Heatplus:
#   source("http://bioconductor.org/biocLite.R")
#   biocLite("Heatplus")
library(Heatplus)

# clustering takes a long time so we are limiting to first 100 genes
shortCountTable=CountTable[1:50,]
mm=data.matrix(shortCountTable, )#
heatmap_2(mm, legend=2, col=topo.colors(100))
# *** note: it is easy enough to load in any table you wish such as filtered results


#------------  Quick PCA  --------
rld <- rlog(dds)
print(plotPCA(rld, intgroup=c("genotype", "stage")))
