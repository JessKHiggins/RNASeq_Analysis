#Analysis of Manudca RNASeq data from CLC genomics workbench, goal is to get some sort of gene name and do gene enrichment analyses.

library(ggplot2)
library(dplyr)

edge_results<- read.csv("EDGE_test.csv", header=TRUE)
names(edge_results)

annotat<- read.csv("Msex_OGS2_Comp_Annotat_Vogel.csv")

#Getting top 50 upregulated for each comparison of interest
CNS_CHS_diff<- edge_results %>% 
  filter(EDGE_CNSvCHS_pvalue< 0.05) %>% 
  select(Seq_name,Description,EDGE_CNSvCHS_pvalue,EDGE_CNSvCHS_Foldchange,EDGE_CNSvCHS_Weighted_diff) %>% 
  arrange(desc(EDGE_CNSvCHS_Foldchange))

head(CNS_CHS_diff)

write.csv(top_n(CNS_CHS_diff,50, EDGE_CNSvCHS_Foldchange),"Top_50_CNSvCHS_up.csv")

CHS_FHS_diff<- edge_results %>% 
  filter(EDGE_CHSvFHS_pvalue< 0.05) %>% 
  select(Seq_name,Description,EDGE_CHSvFHS_pvalue,EDGE_CHSvFHS_Foldchange,EDGE_CHSvFHS_Weighted_diff) %>% 
  arrange(desc(EDGE_CHSvFHS_Foldchange))

head(CHS_FHS_diff)

write.csv(top_n(CHS_FHS_diff,50,EDGE_CHSvFHS_Foldchange),"Top_50_CHSvFHS_up.csv")

FNS_FHS_diff<- edge_results %>% 
  filter(EDGE_FNSvFHS_pvalue< 0.05) %>% 
  select(Seq_name,Description,EDGE_FNSvFHS_pvalue,EDGE_FNSvFHS_Foldchange,EDGE_FNSvFHS_Weighted_diff) %>% 
  arrange(desc(EDGE_FNSvFHS_Foldchange))

head(FNS_FHS_diff)

write.csv(top_n(FNS_FHS_diff,50,EDGE_FNSvFHS_Foldchange),"Top_50_FNSvFHS_up.csv")

CNS_FNS_diff<- edge_results %>% 
  filter(EDGE_CNSvFNS_pvalue< 0.05) %>% 
  select(Seq_name,Description,EDGE_CNSvFNS_pvalue,EDGE_CNSvFNS_Foldchange,EDGE_CNSvFNS_Weighted_diff) %>% 
  arrange(desc(EDGE_CNSvFNS_Foldchange))

head(CNS_FNS_diff)

write.csv(top_n(CNS_FNS_diff,50,EDGE_CNSvFNS_Foldchange),"Top_50_CNSvFNS_up.csv")

CNS_FHS_diff<- edge_results %>% 
  filter(EDGE_CNSvFHS_pvalue< 0.05) %>% 
  select(Seq_name,Description,EDGE_CNSvFHS_pvalue,EDGE_CNSvFHS_Foldchange,EDGE_CNSvFHS_Weighted_diff) %>% 
  arrange(desc(EDGE_CNSvFHS_Foldchange))

head(CNS_FHS_diff)

write.csv(top_n(CNS_FHS_diff,50,EDGE_CNSvFHS_Foldchange),"Top_50_CNSvFHS_up.csv")
  




#Getting Top 50 downregulated for each condition
downCNS_CHS<-arrange(CNS_CHS_diff, EDGE_CNSvCHS_Foldchange) 
write.csv(head(downCNS_CHS, n=50),"Top_50_CNSvCHS_down.csv")

downCNS_FHS<-arrange(CNS_FHS_diff, EDGE_CNSvFHS_Foldchange) 
write.csv(head(downCNS_FHS,n=50),"Top_50_CNSvFHS_down.csv")

downCNS_FNS<-arrange(CNS_FNS_diff, EDGE_CNSvFNS_Foldchange)
write.csv(head(downCNS_FNS,n=50),"Top_50_CNSvFNS_down.csv")

downFNS_FHS<-arrange(FNS_FHS_diff, EDGE_FNSvFHS_Foldchange) 
write.csv(head(downFNS_FHS,n=50),"Top_50_FNSvFHS_down.csv")

downCHS_FHS<-arrange(CHS_FHS_diff, EDGE_CHSvFHS_Foldchange) 
write.csv(head(downCHS_FHS,n=50),"Top_50_CHSvFHS_down.csv")

