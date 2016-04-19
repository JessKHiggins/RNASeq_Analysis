#Analysis of Manudca RNASeq data from CLC genomics workbench, goal is to get some sort of gene name and do gene enrichment analyses.

library(dplyr)
library(ggplot2)

edge_results<- read.csv("EDGE_test.csv", header=TRUE)
names(edge_results)

