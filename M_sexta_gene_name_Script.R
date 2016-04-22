#Analysis of Manudca RNASeq data from CLC genomics workbench, goal is to get some sort of gene name and do gene enrichment analyses.

setwd("~/Google Drive/Kingsolver Lab- Post Doc/M_sexta_RNASeq_Analysis/data")

library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

total_data<- read.csv("Msex_CLC_diff.csv", header=TRUE)
names(total_data)

#edge_results<- read.csv("EDGE_test.csv", header=TRUE)
#names(edge_results)
total_data$Seq_uname<- substr(total_data$Seq_name,1,11)

annotat<- read.csv("Msex_OGS2_Comp_Annotat_Vogel.csv")
head(annotat)

clean_annotat<- annotat %>% 
  select(Name, Seq_Description, Hit_desc.,GOs,GO_Accession,Enzyme_Codes,InterProScan)

clean_annotat$Seq_uname<- substr(clean_annotat$Name,1,11)


names(clean_annotat)
dim(clean_annotat)
names(edge_results)
dim(edge_results)

joined<- right_join(total_data,clean_annotat,"Seq_uname")
names(joined)

#Getting top 50 upregulated for each comparison of interest
CNS_CHS_diff<- joined %>% 
  filter(EDGE_CNSvCHS_pvalue< 0.05) %>% 
  select(Seq_uname,Seq_Description,EDGE_CNSvCHS_pvalue,EDGE_CNSvCHS_Foldchange,EDGE_CNSvCHS_Weighted_diff, CNS_Means,CHS_Means,FNS_Means,FHS_Means) %>% 
  arrange(desc(EDGE_CNSvCHS_Foldchange))

head(CNS_CHS_diff)

write.csv(top_n(CNS_CHS_diff,50, EDGE_CNSvCHS_Foldchange),"Top_50_CNSvCHS_up.csv")

CHS_FHS_diff<- joined %>% 
  filter(EDGE_CHSvFHS_pvalue< 0.05) %>% 
  select(Seq_uname,Seq_Description,EDGE_CHSvFHS_pvalue,EDGE_CHSvFHS_Foldchange,EDGE_CHSvFHS_Weighted_diff,CNS_Means,CHS_Means,FNS_Means,FHS_Means) %>% 
  arrange(desc(EDGE_CHSvFHS_Foldchange))

head(CHS_FHS_diff)

write.csv(top_n(CHS_FHS_diff,50,EDGE_CHSvFHS_Foldchange),"Top_50_CHSvFHS_up.csv")

FNS_FHS_diff<- joined %>% 
  filter(EDGE_FNSvFHS_pvalue< 0.05) %>% 
  select(Seq_uname,Seq_Description,EDGE_FNSvFHS_pvalue,EDGE_FNSvFHS_Foldchange,EDGE_FNSvFHS_Weighted_diff,CNS_Means,CHS_Means,FNS_Means,FHS_Means) %>% 
  arrange(desc(EDGE_FNSvFHS_Foldchange))

head(FNS_FHS_diff)

write.csv(top_n(FNS_FHS_diff,50,EDGE_FNSvFHS_Foldchange),"Top_50_FNSvFHS_up.csv")

CNS_FNS_diff<- joined %>% 
  filter(EDGE_CNSvFNS_pvalue< 0.05) %>% 
  select(Seq_uname,Seq_Description,EDGE_CNSvFNS_pvalue,EDGE_CNSvFNS_Foldchange,EDGE_CNSvFNS_Weighted_diff,CNS_Means,CHS_Means,FNS_Means,FHS_Means) %>% 
  arrange(desc(EDGE_CNSvFNS_Foldchange))

head(CNS_FNS_diff)

write.csv(top_n(CNS_FNS_diff,50,EDGE_CNSvFNS_Foldchange),"Top_50_CNSvFNS_up.csv")

CNS_FHS_diff<- joined %>% 
  filter(EDGE_CNSvFHS_pvalue< 0.05) %>% 
  select(Seq_uname,Seq_Description,EDGE_CNSvFHS_pvalue,EDGE_CNSvFHS_Foldchange,EDGE_CNSvFHS_Weighted_diff,CNS_Means,CHS_Means,FNS_Means,FHS_Means) %>% 
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


#Figuring out how many are turned up vs. down

CNS_CHS_diff %>% 
  filter(EDGE_CNSvCHS_Foldchange>0) %>% 
  count()

FNS_FHS_diff %>% 
  filter(EDGE_FNSvFHS_Foldchange<0) %>% 
  count()

CNS_FNS_diff %>% 
  filter(EDGE_CNSvFNS_Foldchange<0) %>% 
  count()

CHS_FHS_diff %>% 
  filter(EDGE_CHSvFHS_Foldchange<0) %>% 
  count()

#Read in CLC_diff and look at overall levels 
#Means are total gene reads for each condition

#Means that are significant in all conditions

Means_join<- full_join()
names(Means_join)

top1<-top_n(CNS_CHS_diff,50, EDGE_CNSvCHS_Foldchange)
top2<-top_n(CHS_FHS_diff,50,EDGE_CHSvFHS_Foldchange)
top3<- top_n(CNS_FNS_diff,50,EDGE_CNSvFNS_Foldchange)
top4<- top_n(FNS_FHS_diff,50,EDGE_FNSvFHS_Foldchange)

down1<- downCNS_CHS[1:50,]
down2<- downCHS_FHS[1:50,]
down3<- downCNS_FNS[1:50,]
down4<- downFNS_FHS[1:50,]

join1<- full_join(top1,top2)
join2<- full_join(top3,top4)
means_joinup<-full_join(join1,join2)

View(means_join)
names(means_join)

means_join_up<- means_joinup %>% 
  select(Seq_uname,Seq_Description, CNS_Means, CHS_Means, FNS_Means,FHS_Means)

names(means_join_up)

join3<- full_join(down1,down2)
join4<- full_join(down3,down4)
means_joindown<-full_join(join3,join4) 

means_join_down<- means_joindown %>% 
  select(Seq_uname,Seq_Description, CNS_Means, CHS_Means, FNS_Means,FHS_Means)

all_means<- full_join(means_join_up,means_join_down)
head(all_means)
View(all_means)


#Filtering out tryp genes
tryp_means<- all_means %>% 
  filter(grepl("tryp", Seq_Description) )

melt_tryp<-melt(tryp_means, id=c("Seq_uname", "Seq_Description")) %>% 
  rename(Condition=variable, Means=value)


#Filtering out immue genes bye using "cin"

immune_means<- all_means %>% 
  filter(grepl("cin", Seq_Description))

melt_immune<-melt(immune_means, id=c("Seq_uname", "Seq_Description")) %>% 
  rename(Condition=variable, Means=value)

#Filtering out heat shock genes

hs_means<- all_means %>% 
  filter(grepl("heat", Seq_Description))

melt_hs<-melt(hs_means, id=c("Seq_uname", "Seq_Description")) %>% 
  rename(Condition=variable, Means=value)


#Rough plots where each color is a sequence description and each point is a unique gene
ggplot(melt_tryp, aes(x=Condition, y=Means, color=Seq_Description))+
  geom_point()+geom_line(aes(group=Seq_uname))

ggplot(melt_immune, aes(x=Condition, y=Means, color=Seq_Description))+
  geom_point()+geom_line(aes(group=Seq_uname))

ggplot(melt_hs, aes(x=Condition, y=Means, color=Seq_Description))+
  geom_point()+geom_line(aes(group=Seq_uname))

