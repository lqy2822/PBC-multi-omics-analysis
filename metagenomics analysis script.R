##############metadata read in############
#metadata
setwd("D:/PBC multiomics data/metagenomics/data")
metadata=read.csv("PBC metagenomics sample metadata.csv")
row.names(metadata)=metadata$SampleID
metadata$group=as.factor(metadata$group)
metadata$Stage=as.factor(metadata$Stage)
metadata$Antibiotics[262:411]=c(rep("0",150))
metadata$Antibiotics=as.numeric(metadata$Antibiotics)
metadata$Strorage_time=as.numeric(metadata$Strorage_time)
metadata$subjectID=as.factor(metadata$subjectID)
metadata$gender=as.factor(metadata$gender)
metadata$gender[metadata$gender=="F "]<-"F"
metadata$gender[metadata$gender=="M "]<-"M"
metadata$gender=metadata$gender[,drop=TRUE]
metadata$Treatment=as.factor(metadata$Treatment)
metadata$Response_6m=as.factor(metadata$Response_6m)
metadata$Response_1y=as.factor(metadata$Response_1y)
metadata$gp210=as.factor(metadata$gp210)
metadata$Sp100=as.factor(metadata$Sp100)
metadata$Paris.I=as.factor(metadata$Paris.I)
metadata$Barcelona=as.factor(metadata$Barcelona)
metadata$Rotterdam=as.factor(metadata$Rotterdam)
metadata$Toronto=as.factor(metadata$Toronto)
str(metadata$Stage)
colnames(metadata)
metadata[,c(20:31)]=sapply(metadata[,20:31],as.numeric)
str(metadata[PBCbase,"Stage"])

#derive the sampleID in each subgroup
library(dplyr)
HC=row.names(filter(metadata,group=="HC"))
HC_no=c("MXH063","MXH1022","MXH1076","MXH1136","MXH444","MXH448","MXH459",
                "MXH496","MXH497","MXH532","MXH565","MXH779","MXH827","MXH840","MXH966","MXH976",
                "MXH455","MXH605","MXH572")
index=HC %in% HC_no #为逻辑值向量
index=which(index==TRUE)#为位置索引
HC_new=HC[-index]
#remove classical PBCbase patients with antibiotics usage
PBCbase=row.names(filter(metadata,group=="PBCbase"& Antibiotics=="0"))
scPBCbase=row.names(filter(metadata,group=="scPBCbase"))#several scPBC patients took antibiotics
PBC_6m=row.names(filter(metadata,group=="PBC_6m"&Antibiotics=="0"))
PBCpair_6m=substring(PBC_6m,1,6)%>% intersect(.,PBCbase)
PBC_6m_pair=paste(PBCpair_6m,"2",sep="-") %>% intersect(.,PBC_6m)

#baseline samples of UDCA 6month responders and non-responder 
PBCbase_resp_6m=row.names(filter(metadata,group=="PBCbase"& Antibiotics=="0"& Response_6m=="1"))
PBCbase_noresp_6m=row.names(filter(metadata,group=="PBCbase"& Antibiotics=="0"& Response_6m=="0"))
PBCpair_6m_resp=intersect(PBCbase_resp_6m,PBCpair_6m)
PBCpair_6m_noresp=intersect(PBCbase_noresp_6m,PBCpair_6m)


PBCbase_cluster1=filter(cluster.res,cluster=="clt1") %>% row.names()
PBCbase_cluster2=filter(cluster.res,cluster=="clt2") %>% row.names()

#cluster1中response和non-response样本
PBCbase_cluster1_resp=intersect(PBCbase_cluster1,PBCbase_resp_6m)
PBCbase_cluster1_noresp=intersect(PBCbase_cluster1,PBCbase_noresp_6m)
PBCbase_cluster2_resp=intersect(PBCbase_cluster2,PBCbase_resp_6m)
PBCbase_cluster2_noresp=intersect(PBCbase_cluster2,PBCbase_noresp_6m)
  
#longitudinal sample ID:response
#6month sample 
PBC_6m_resp=substring(PBC_6m,1,6) %in% PBCpair_6m_resp %>% which(.==TRUE) %>% PBC_6m[.]
PBC_6m_noresp=substring(PBC_6m,1,6) %in% PBCpair_6m_noresp %>% which(.==TRUE) %>% PBC_6m[.]

#longitudinal sample ID:cluster
PBCpair_6m_cluster1=intersect(PBCpair_6m,PBCbase_cluster1)
PBCpair_6m_cluster2=intersect(PBCpair_6m,PBCbase_cluster2)
PBC_6m_cluster1=substring(PBC_6m,1,6) %in% PBCpair_6m_cluster1 %>% which(.==TRUE) %>% PBC_6m[.]
PBC_6m_cluster2=substring(PBC_6m,1,6) %in% PBCpair_6m_cluster2%>% which(.==TRUE) %>% PBC_6m[.]

#所有分析的样本
all_sample=c(PBCbase,PBC_6m,scPBCbase,HC)

#compare arthrimatic parameters among PBC and HC
wilcox.test(metadata[PBCbase,"age"],metadata[HC,"age"])
summary(metadata[HC,"age"])
chisq.test(metadata[c(PBCbase,HC),"gender"],metadata[c(PBCbase,HC),"group"])$observed

###########phylum data############
setwd("D:/multi-omics data/metagenomics data/re-assembled data/data/phylum")
phylum_rela=read.csv("phylum relative abundance renamed.csv",header = T,row.names=1)
colnames(phylum_rela)=all_name
phylum_rela=phylum_rela[,row.names(metadata)]
phylum_rela[,1:411]=sapply(phylum_rela[,1:411],as.numeric)
colSums(phylum_rela)
#keep phylum with a minimum relative abundance of 0.001% and in at least 1% samples 
phylum_rela_sample=phylum_rela[,all_sample]
phylum_rela_filter=phylum_rela_sample
phylum_rela_filter[phylum_rela_filter<0.00001]=0
phylum_rela_filter[phylum_rela_filter>0.00001]=1
phylum_rela_filter<- phylum_rela_sample[which(rowSums(phylum_rela_filter) >=4),]
#absolute abundance
phylum_abs=read.csv("phylum absolute abundance renamed.csv",header = T,row.names = 1)
colnames(phylum_abs)=all_name
library(compositions)
library(zCompositions)
phylum_abs_clr=cmultRepl(phylum_abs, label=0, method="CZM")  %>% clr(.) #zCompositions package
phylum_abs_clr=t(phylum_abs_clr) %>% as.data.frame()

##########family data###########
setwd("D:/multi-omics data/metagenomics data/re-assembled data/data/family")
family_rela=read.csv("family relative abundance renamed.csv",header = T,row.names=1)
colnames(family_rela)=all_name
family_rela=family_rela[,row.names(metadata)]
family_rela[,1:411]=sapply(family_rela[,1:411],as.numeric)
colSums(family_rela)
#keep family with a minimum relative abundance of 0.001% and in more than 1% samples 
family_rela_sample=family_rela[,all_sample]
family_rela_filter=family_rela_sample
family_rela_filter[family_rela_filter<0.00001]=0
family_rela_filter[family_rela_filter>0.00001]=1
family_rela_filter<- family_rela_sample[which(rowSums(family_rela_filter) >=4),]
#family absolute data
setwd("D:/multi-omics data/metagenomics data/re-assembled data/data/family")
family_abs=read.csv("family abs abundance renamed.csv",header = T,row.names = 1)
colnames(family_abs)=all_name
family_abs_clr=cmultRepl(family_abs, label=0, method="CZM")  %>% clr(.) #zCompositions package
family_abs_clr=t(family_abs_clr) %>% as.data.frame()
###################class data and plot##################
setwd("D:/multi-omics data/metagenomics data/re-assembled data/data/class")
class_abs=read.csv("class absolute abundance renamed.csv",header = T,row.names=1)
class_rela=read.csv("class relative abundance renamed.csv",header = T,row.names = 1)
colnames(class_abs)=all_name
colnames(class_rela)=all_name
#clr transformation
library(compositions)
library(zCompositions)
#all taxa
class_abs_clr=cmultRepl(class_abs, label=0, method="CZM")  %>% clr(.)
class_abs_clr=t(class_abs_clr) %>% as.data.frame()
class_abs_clr=class_abs_clr[row.names(metadata),]
class_abs_clr_scale=scale(class_abs_clr[all_sample,]) %>% as.data.frame()#只对分析的样本进行scale
#keep class with a relative abundance 0.001% in at least 1% samples
class_rela_sample=class_rela[,all_sample]
library(plyr)
class_rela_sample=class_rela_sample %>% mutate_each(funs(. / sum(.)))
colSums(class_rela_sample)
class_rela_filter=class_rela_sample
class_rela_filter[class_rela_filter<0.00001]=0
class_rela_filter[class_rela_filter>0.00001]=1
class_rela_filter<- class_rela_sample[which(rowSums(class_rela_filter) >=4),]
row.names(class_rela_filter)
#计算每个class 在PBC的平均相对丰度
class_mean=apply(class_rela_filter[,PBCbase], 1, mean) %>% as.data.frame(.)
row.names(class_mean)
colnames(class_mean)="class_mean"
class_mean=arrange(class_mean,desc(class_mean))
row.names(class_mean)
#丰度最高的10个class
class_top10=row.names(class_mean)[c(2,1,4,7:13)]
#只绘制PBCbase的top10 known class 变化图
class_rela_filter_plot=class_rela_filter[class_top10,PBCbase] %>% t(.) %>% as.data.frame(.)
class_rela_filter_plot=arrange(class_rela_filter_plot,desc(class_rela_filter_plot$`k__Bacteria;p__Firmicutes;c__Clostridia`))
class_rela_filter_plot=t(class_rela_filter_plot)%>% as.data.frame(.)
class_rela_filter_plot$Taxonomy <- factor(rownames(class_rela_filter_plot), levels = rev(rownames(class_rela_filter_plot)))
class_rela_filter_plot <- melt(class_rela_filter_plot, id = 'Taxonomy')
colnames(class_rela_filter_plot)

p <- ggplot(class_rela_filter_plot, aes(variable, 100 * value, fill = Taxonomy)) +
  geom_col(position = 'stack', width = 0.6) +
  scale_fill_manual(values =  rev(c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', 'gray'))) +
  labs(x = '', y = 'Relative Abundance(%)') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13), legend.title = element_blank(), legend.text = element_text(size = 11))
p
#########order data################
setwd("D:/multi-omics data/metagenomics data/re-assembled data/data/order")
order_rela=read.csv("order relative abundance renamed.csv",header = T,row.names=1)
colnames(order_rela)=all_name
order_rela=order_rela[,row.names(metadata)]
order_rela[,1:411]=sapply(order_rela[,1:411],as.numeric)
colSums(order_rela)
#keep order with a minimum relative abundance of 0.001% and in more than 5% samples 
order_rela_filter=order_rela
order_rela_filter[order_rela_filter<0.00001]=0
order_rela_filter[order_rela_filter>=0.00001]=1
order_rela_filter<- order_rela[which(rowSums(order_rela_filter) >= 4),] 
#absolute abundance
order_abs=read.csv("order absolute abundance renamed.csv",header = T,row.names = 1)
colnames(order_abs)=all_name
order_abs=order_abs[,row.names(metadata)]

order_abs_clr=cmultRepl(order_abs[row.names(order_rela_filter),], label=0, method="CZM")  %>% clr(.) 
#在进行Bayes推断时因为存在 Row(s) containing all zeros/unobserved values were found (check it out using zPatterns)，
#所以此处只进行了满足filter标准taxa的clr转化
order_abs_clr=t(order_abs_clr) %>% as.data.frame()


###########genus data############
#abslute abundance data
setwd("D:/multi-omics data/metagenomics data/re-assembled data/data/genus")
genus_abs=read.csv("genus_abs abundance renamed.csv",header = T,row.names=1)
colnames(genus_abs)=all_name
#clr transformation
library(compositions)
library(zCompositions)
#all taxa
genus_abs_clr=cmultRepl(genus_abs, label=0, method="CZM")  %>% clr(.)
genus_abs_clr=t(genus_abs_clr) %>% as.data.frame()
genus_abs_clr=genus_abs_clr[row.names(metadata),]
genus_abs_clr_scale=scale(genus_abs_clr[all_sample,]) %>% as.data.frame()#只对分析的样本进行scale
#read in relative abundance data
genus_rela=read.csv("genus_relative_abundance_renamed.csv",header = T,row.names = 1)
colnames(genus_rela)=all_name
colSums(genus_rela)
genus_rela=genus_rela[,row.names(metadata)]
#keep genus with a relative abundance 0.001% in at least 1% samples
genus_rela_sample=genus_rela[,all_sample]
genus_rela_sample=genus_rela_sample %>% mutate_each(funs(. / sum(.)))
colSums(genus_rela_sample)
genus_rela_filter=genus_rela_sample
genus_rela_filter[genus_rela_filter<0.00001]=0
genus_rela_filter[genus_rela_filter>0.00001]=1
genus_rela_filter<- genus_rela_sample[which(rowSums(genus_rela_filter) >=4),]

##########species data############
#read absolute abundance file
setwd("D:/multi-omics data/metagenomics data/re-assembled data/data/species")
species_abs=read.csv("species abs renamed.csv",header = T,row.names=1)
colnames(species_abs)=all_name
species_abs_clr=cmultRepl(species_abs, label=0, method="CZM")  %>% clr(.)
species_abs_clr=t(species_abs_clr) %>% as.data.frame()

#relative abundance data
species_rela=read.csv("species rela renamed.csv",header = T,row.names = 1)
colnames(species_rela)=all_name
species_rela=species_rela[,row.names(metadata)]
colSums(species_rela)
#keep species with relative abundance 0.00001 in at least 1% samples
species_rela_sample=species_rela[,all_sample]
species_rela_sample=species_rela_sample %>% mutate_each(funs(. / sum(.)))
colSums(species_rela_sample)
species_rela_filter=species_rela_sample
species_rela_filter[species_rela_filter<0.00001]=0
species_rela_filter[species_rela_filter>0.00001]=1
species_rela_filter<- species_rela_sample[which(rowSums(species_rela_filter) >=4),]

#core species present in more than 95% samples
species_rela_filter_core=species_rela_sample
species_rela_filter_core[species_rela_filter_core!=0]=1
species_rela_filter_core<- species_rela_sample[which(rowSums(species_rela_filter_core) >=344),]
rowSums(species_rela_filter_core)
##############KEGG pathway##############
KEGG_path3_rela=read.csv("KEGG level3 pathway relative abundance.csv",header = T,row.names = 1)
colnames(KEGG_path3_rela)=KEGG_path3_rela[1,]
KEGG_path3_rela=KEGG_path3_rela[-1,]
colnames(KEGG_path3_rela)=all_name
KEGG_path3_rela=KEGG_path3_rela[,row.names(metadata)]
KEGG_path3_rela[,1:411]=sapply(KEGG_path3_rela[,1:411], as.numeric)
colSums(KEGG_path3_rela)
KEGG_path3_rela_filter=KEGG_path3_rela[,all_sample]
KEGG_path3_rela_filter[KEGG_path3_rela_filter<0.001]=0
KEGG_path3_rela_filter[KEGG_path3_rela_filter>0.001]=1
KEGG_path3_rela_filter<- KEGG_path3_rela[,all_sample][which(rowSums(KEGG_path3_rela_filter) >=4),] 


############KO gene#############
KO_gene_abs=read.csv("KO gene absolute abundance.csv",header = T,row.names=1)
colnames(KO_gene_abs)=all_name
KO_gene_abs=KO_gene_abs[,row.names(metadata)]
KO_gene_rela=read.csv("KO gene.csv",header = T,row.names = 1)
colnames(KO_gene_rela)=all_name
KO_gene_rela=KO_gene_rela[,row.names(metadata)]
colSums(KO_gene_rela)
#clr transformation
KO_gene_abs_clr=cmultRepl(KO_gene_abs, label=0, method="CZM")  %>% clr(.) #zCompositions package
KO_gene_abs_clr=t(KO_gene_abs_clr) %>% as.data.frame()
#keep KO_gene with a prevalence in at least more than 1% sample
KO_gene_rela_filter=KO_gene_rela[,all_sample]
gene_presence=vector()
for (i in 1:11276) {
  gene_presence[i]=which(KO_gene_abs[,all_sample][i,]!=0) %>% length(.)
}
gene_presence=as.data.frame(gene_presence)
row.names(gene_presence)=row.names(KO_gene_abs)
index=gene_presence$gene_presence>=4 %>% as.vector() 
index=which(index==TRUE)
KO_gene_rela_filter=KO_gene_rela[index,all_sample]
setwd("D:/PBC multiomics data/analysis V6/metagenomics/data")
write.csv(gene_presence,"gene_presence.csv")
#gene annotation
setwd("D:/PBC multiomics data/unsupervised cluster analysis_renamed/metagenomics/result/KO gene")
library(dplyr)
gene_anno=read.csv("gene_anno.csv",header = T,row.names = 1)#
row.names(KO_gene_filter) %in% gene_anno$KO
row.names(gene_anno)=gene_anno$KO
#只注释筛选后的基因
gene_anno_filter=gene_anno[row.names(KO_gene_rela_filter),]


###############Cazy data############
#Cazy 1 data
#read in data
setwd("D:/multi-omics data/metagenomics data/re-assembled data/data/pathway")
Cazy1_abs=read.csv("CAzy1 absolute abundance.csv",header = T,row.names=1)
library(zCompositions)
library(compositions)
colnames(Cazy1_abs)=all_name
CAzy1_abs_clr= clr(Cazy1_abs) #zCompositions package,no zero in data
CAzy1_abs_clr=t(CAzy1_abs_clr) %>% as.data.frame()

Cazy1_rela=read.csv("CAzy1 relative abundance.csv",header = T,row.names = 1)
colnames(Cazy1_rela)=all_name

##Cazy2 data
Cazy2_abs=read.csv("CAzy2 absolute abundance.csv",header = T,row.names=1)
colnames(Cazy2_abs)=all_name
Cazy2_abs_clr= cmultRepl(Cazy2_abs, label=0, method="CZM")  %>% clr(.) #zCompositions package,no zero in data
Cazy2_abs_clr=t(Cazy2_abs_clr) %>% as.data.frame()

Cazy2_rela=read.csv("CAzy2 relative abundance.csv",header = T,row.names = 1)
colnames(Cazy2_rela)=all_name
#cazy2 gene filter
Cazy2_rela_filter=Cazy2_rela[,all_sample]
Cazy2_rela_filter[Cazy2_rela_filter<0.001]=0
Cazy2_rela_filter[Cazy2_rela_filter>0.001]=1
Cazy2_rela_filter<- Cazy2_rela[which(rowSums(Cazy2_rela_filter) >=4),] 
#cazy2 gene anntation
cazy2_annotation=read.table("CAZy_Annotation.txt",header = T,sep = "\t")
colnames(cazy2_annotation)
cazy2_anno=cazy2_annotation[!duplicated(cazy2_annotation$CAZy_family),] %>% as.data.frame(.)
row.names(cazy2_anno)=cazy2_anno$CAZy_family



############PAM cluster analysis###############
library("factoextra")
library(cluster)
#determine the best number of clustering
#silhouette method
clust.num=fviz_nbclust(genus_abs_clr_scale[PBCbase,row.names(genus_rela_filter)], pam, method = "silhouette")
#wss method
clust.num=fviz_nbclust(genus_abs_clr_scale[PBCbase,row.names(genus_rela_filter)], pam, method = "wss")
#gap_stat
clust.num=fviz_nbclust(genus_abs_clr_scale[PBCbase,row.names(genus_rela_filter)], pam, method = "gap_stat")
plot(clust.num)
#fpc package
library(fpc)
pamk(genus_abs_clr_scale[PBCbase,row.names(genus_rela_filter)],criterion = "ch") %>% .$nc # the best number of cluster
#visualize the Calinski-Harabasz index calculated from pamk function
CH_index=pamk(genus_abs_clr_scale[PBCbase,row.names(genus_rela_filter)],criterion = "ch") %>% .$crit %>% as.data.frame(.)
colnames(CH_index)="CH_index"
CH_index$cluster=as.factor(c(1:10))
CH_index
boxplot(CH_index~cluster,data = CH_index)
#final results: the number of cluster is 2.
pam.res <- pam(genus_abs_clr_scale[PBCbase,row.names(genus_rela_filter)],2)
# Visualize
library("factoextra")
cluster.plot=fviz_cluster(pam.res, data =genus_abs_clr_scale[PBCbase,row.names(genus_rela_filter)],
                          ellipse = FALSE,
                          palette = "jco",
                          repel = FALSE,
                          labelsize = 0,
                          show.clust.cent = "TRUE",
                          ggtheme = theme_minimal())
cluster.plot
cluster.res=as.data.frame(pam.res$clustering)
colnames(cluster.res)="cluster"
cluster.res$cluster=as.factor(cluster.res$cluster)
cluster.res$cluster= c('1'="clt1",'2'="clt2")[ as.factor(cluster.res$cluster)]
cluster.res$cluster=as.factor(cluster.res$cluster)
summary(cluster.res$cluster)

#整合metadata cluster信息
metadata_PBCbase=cbind(metadata[PBCbase,],cluster.res)
colnames(metadata_PBCbase)[39]="cluster"
metadata_PBCbase=cbind(metadata_PBCbase,shannon[row.names(metadata_PBCbase),"shannon"])
colnames(metadata_PBCbase)[40]="shannon"
HC_cluster=c(rep("HC",150)) %>% as.data.frame(.)
row.names(HC_cluster)=HC
colnames(HC_cluster)="cluster"
metadata_PBCHC=rbind(cluster.res,HC_cluster)
colnames(metadata_PBCHC)[39]="cluster"
metadata_PBCHC=cbind(metadata_PBCHC,shannon[c(PBCbase,HC),])
colnames(metadata_PBCHC)[40]="shannon"
metadata_PBCHC$SampleID=row.names(metadata_PBCHC)
metadata_PBCHC=cbind(metadata_PBCHC,metadata[c(PBCbase,HC),"group"])
metadata_PBCHC$group=factor(metadata_PBCHC$group,levels = c("PBCbase","HC"))
colnames(metadata_PBCHC)[4]="group"
metadata_PBCHC$group=metadata_PBCHC$group[,drop=TRUE]
metadata_PBCHC=cbind(metadata[row.names(metadata_PBCHC),],metadata_PBCHC[,"cluster"])

#cluster.res_PBCHC添加NLR信息，重新读入
setwd("D:/PBC multiomics data/analysis V6/metagenomics/data")
cluster.res_PBCHC.NLR=read.csv("cluster.res_PBCHC.csv",header = T,row.names=1)
row.names(metadata_PBCHC)
row.names(cluster.res_PBCHC.NLR)
metadata_PBCHC$NLR=cluster.res_PBCHC.NLR$NLR
metadata_PBCbase$NLR=cluster.res_PBCHC.NLR[1:132,]$NLR


#calculate cluster assocaition with phenotype
#categorial factor
chisq.test(metadata[PBCbase,]$Response_6m,cluster.res[PBCbase,1])$observed 
chisq.test(metadata[PBCbase,]$gp210,cluster.res[,1])$observed
chisq.test(metadata[PBCbase,]$gp210,metadata[PBCbase,]$Response_6m)$observed
#选择有Stage信息的样品
library(dplyr)
stage_sample=filter(metadata,Stage=="early"|Stage=="advanced") %>% row.names(.) %>% intersect(.,PBCbase)
metadata_stage=metadata[stage_sample,]$Stage
metadata_stage=as.data.frame(metadata_stage)
colnames(metadata_stage)="stage"
metadata_stage$stage=as.factor(metadata_stage$stage)
metadata_stage$stage=metadata_stage$stage[,drop=TRUE]
str(metadata_stage$stage)
row.names(metadata_stage)=stage_sample
chisq.test(metadata_stage[stage_sample,1],cluster.res[stage_sample,1])#P=0.5024


ggplot(metadata_PBCbase,aes(x=cluster,y=TB,color=cluster))+
  geom_boxplot()+
  geom_jitter(width = 0.2,size=0.9)+
  scale_color_manual(values = c("#088288","#F19433"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "null")
#TB绘图移出极端值

metadata_PBCbase$TB

library(plyr)

#response
clt_resp=chisq.test(metadata[PBCbase,]$Response_6m,cluster.res[PBCbase,1])$observed %>% as.data.frame()
colnames(clt_resp)=c("Response","cluster","Freq")
clt_resp$percent=c(0.2,0.8,0.41,0.59)
clt_resp$Freq
clt_resp$cluster
library(ggplot2)
clt_resp$percent_rev=c("0.36","0.62","0.64","0.38")
ggplot(clt_resp,aes(x=cluster,y=percent,fill=Response))+
  geom_bar(stat = "identity",position = position_stack(reverse = TRUE))+
  scale_fill_manual(values = c("#85ABD1","#BFBFBF"))+
  geom_text(aes(label=percent),position = position_stack(reverse = TRUE),color="white",vjust=2)+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank())



#gp210
clt_gp210=chisq.test(metadata[PBCbase,]$gp210,cluster.res[PBCbase,1])$observed %>% as.data.frame()
colnames(clt_gp210)=c("gp210","cluster","Freq")
clt_gp210$percent=c(0.71,0.29,0.57,0.43)
clt_gp210$percent_rev=c("67%","38%","33%","62%")
clt_gp210=clt_gp210[c(1,2,4,3),]
ggplot(clt_gp210,aes(x=cluster,y=percent,fill=gp210_new))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#85ABD1","#BFBFBF"))+
  geom_text(aes(label=percent),position = position_stack(),color="white",vjust=2)+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank())


#continual variable

wilcox.test(shannon~cluster,data=metadata_PBCbase) #p=6.1e-13
wilcox.test(ALB~cluster,data=metadata_PBCbase) #p=0.24
wilcox.test(age~cluster,data=metadata_PBCbase) #p=0.013
wilcox.test(IgM~cluster,data=metadata_PBCbase) #p=0.07
wilcox.test(ALT~cluster,data=metadata_PBCbase) #p=0.28
wilcox.test(AST~cluster,data=metadata_PBCbase) #p=0.98
wilcox.test(TB~cluster,data=metadata_PBCbase) #p=0.17
wilcox.test(ALP~cluster,data=metadata_PBCbase) #p=0.23
wilcox.test(age~Response_6m,data=metadata_PBCbase) #p=0.013
summary(metadata_PBCbase$cluster)
lm(NLR~cluster+age,data = metadata_PBCbase) %>% summary(.)
#Association of cluster and response after deconfounding age
glm(Response_6m~cluster+age+ALP+TB,data=metadata_PBCbase,family = "binomial") %>% summary()
glm(gp210~cluster+age,data=metadata_PBCbase,family = "binomial") %>% summary()
library(dplyr)
cor.test(metadata_PBCbase$age,metadata_PBCbase$shannon,method = "spearman")

#############PAM cluster analysis of treated samples#####################
#1.治疗后样本的cluster分析
library("factoextra")
library(cluster)
#determine the best number of clustering
#silhouette method
clust.num_6m=fviz_nbclust(genus_abs_clr_scale[PBC_6m_pair,row.names(genus_rela_filter)], pam, method = "silhouette")
#wss method
clust.num_6m=fviz_nbclust(genus_abs_clr_scale[PBC_6m_pair,row.names(genus_rela_filter)], pam, method = "wss")
#gap_stat
clust.num_6m=fviz_nbclust(genus_abs_clr_scale[PBC_6m_pair,row.names(genus_rela_filter)], pam, method = "gap_stat")

#final results: the number of cluster is 2.
pam.res_6m <- pam(genus_abs_clr_scale[PBC_6m_pair,row.names(genus_rela_filter)],2)
# Visualize
library("factoextra")
cluster.plot_6m=fviz_cluster(pam.res_6m, data =genus_abs_clr_scale[PBC_6m_pair,row.names(genus_rela_filter)],
                          ellipse = FALSE,
                          palette = "jco",
                          repel = FALSE,
                          labelsize = 0,
                          show.clust.cent = "TRUE",
                          ggtheme = theme_minimal())
cluster.plot_6m
cluster.res_6m=as.data.frame(pam.res_6m$clustering)
colnames(cluster.res_6m)="cluster"
cluster.res_6m$cluster=as.factor(cluster.res_6m$cluster)
cluster.res_6m$cluster= c('1'="clt1",'2'="clt2")[ as.factor(cluster.res_6m$cluster)]
cluster.res_6m$cluster=as.factor(cluster.res_6m$cluster)
summary(cluster.res_6m$cluster)

#分析两个cluster的多样性差异
wilcox.test(shannon[PBC_6m_clt1_test,],shannon[PBC_6m_clt2_test,])
boxplot(shannon[PBC_6m_clt1_test,],shannon[PBC_6m_clt2_test,])

PBC_6m_clt1_test=filter(cluster.res_6m,cluster=="clt1") %>% row.names()
PBC_6m_clt2_test=filter(cluster.res_6m,cluster=="clt2") %>% row.names()

#2.analyze the consistency of baseline and treated sample cluster membership
PBC_6m_clt1_test
PBC_6m_cluster1
baseline_cluster=as.data.frame(c(rep("clt1",36),rep("clt2",23)))
row.names(baseline_cluster)=c(PBC_6m_cluster1,PBC_6m_cluster2)
colnames(baseline_cluster)="baseline"

treat_cluster=as.data.frame(c(rep("clt1_treat",41),rep("clt2_treat",18)))
row.names(treat_cluster)=c(PBC_6m_clt1_test,PBC_6m_clt2_test)
colnames(treat_cluster)="treat"

cluster_combine=cbind(baseline_cluster,treat_cluster[row.names(baseline_cluster),])
colnames(cluster_combine)[2]="treat"

chisq.test(cluster_combine[,"baseline"],cluster_combine[,"treat"])$observed #P=0.009

#结论：治疗后样本也存在cluster，且和基线的cluster一致；P=0.009

#3.分析治疗后cluster与临床应答的相关性
substring(row.names(treat_cluster),1,6)
chisq.test(metadata[substring(row.names(treat_cluster),1,6),]$Response_6m,treat_cluster[,1]) #P=0.43

#结论：治疗后cluster与临床应答无相关性

#################alpha diversity HC_new############
#(上文代码中shannon源自此分析)
library(vegan)
shannon=diversity(t(species_rela_filter)) %>% as.data.frame()
colnames(shannon)="shannon"
setwd("D:/PBC multiomics data/analysis V6/multiomics/data")
write.csv(shannon,"shannon.csv")
#with HC_new

wilcox.test(shannon[PBCbase_cluster1,],shannon[HC_new,])
wilcox.test(shannon[HC_new,],shannon[PBCbase,])


#plot according to group HC_new
metadata_PBCHC.new$group=factor(metadata_PBCHC.new$group,levels = c("HC","PBCbase"))
metadata_PBCHC.new=cbind(metadata_PBCHC.new,shannon[row.names(metadata_PBCHC.new),1])
colnames(metadata_PBCHC.new)[5]="shannon"
metadata_PBCHC.new=cbind(metadata_PBCHC.new,metadata_PBCHC[row.names(metadata_PBCHC.new),"cluster"])
colnames(metadata_PBCHC.new)[6]="cluster"
metadata_PBCHC.new$cluster=factor(metadata_PBCHC.new$cluster,levels = c("HC","clt1","clt2"))
#i.boxplot with jitter
ggplot(data = metadata_PBCHC.new,aes(x=group,y=shannon,color=group))+geom_boxplot()+
  geom_jitter(width = 0.2,size=0.9)+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  scale_color_manual(values = c("#3F5BA7","#D72326"))

#ii boxplot without jitter
ggplot(data = metadata_PBCHC.new,aes(x=group,y=shannon,fill=group))+geom_boxplot(color="grey")+
  
  theme(axis.line = element_line(colour = "grey"),panel.background = element_blank(),legend.position = "none")+
  scale_fill_manual(values = c("#3F5BA7","#D72326"))

#plot according to cluster 
#i.without HC
ggplot(data = metadata_PBCbase,aes(x=cluster,y=shannon,fill=cluster))+geom_boxplot()+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  scale_fill_manual(values = c("#088288","#F19433"))
#ii with HC_new
ggplot(data = metadata_PBCHC.new,aes(x=cluster,y=shannon,fill=cluster))+geom_boxplot(color="grey")+
  
  theme(axis.line = element_line(colour = "grey"),panel.background = element_blank(),legend.position = "none")+
  scale_fill_manual(values = c("#3F5BA7","#088288","#F19433"))


############beta diversity analysis##################
#Aitchison distance
library(robCompositions)
library(zCompositions)
library(compositions)
library(vegan)
#Distance calculation
Ait_species_abs=vegdist(species_abs_clr[,row.names(species_rela_filter)],method = "euclidean")%>% as.matrix(.)
Ait_gene_abs=vegdist(KO_gene_abs_clr[,row.names(KO_gene_rela_filter)],method = "euclidean")%>% as.matrix(.)
Ait_KEGG3_abs=vegdist(KEGG_path3_abs_clr[,row.names(KEGG_path3_rela_filter)],method = "euclidean") %>% as.matrix(.)
Ait_genus_abs=vegdist(genus_abs_clr[,row.names(genus_rela_filter)],method = "euclidean") %>% as.matrix(.)

#1.PBC+HC
#Plot the PCoA based on Aitchison distance among PBCbase and HC
pcoa_HC_PBC_Ait<- cmdscale(Ait_genus_abs[c(HC_new,PBCbase),c(HC_new,PBCbase)], k = nrow(metadata[c(HC_new,PBCbase),]) - 1, eig = TRUE, add = TRUE)

#Extract PCoA coordinates
pc12_HC_PBC_Ait<- pcoa_HC_PBC_Ait$points[,1:2]
pc12_HC_PBC_Ait=as.data.frame(pc12_HC_PBC_Ait)
pc_importance_HC_PBC_Ait<-round(pcoa_HC_PBC_Ait$eig/sum(pcoa_HC_PBC_Ait$eig)*100,digits = 2)
pc12_HC_PBC_Ait$SampleID<-row.names(pc12_HC_PBC_Ait)
data_HC_PBC_Ait<-merge(pc12_HC_PBC_Ait,metadata_PBCHC[c(HC_new,PBCbase),],by="SampleID")

#根据group
pcoa_plot_HC_PBC_Ait<-ggplot(data_HC_PBC_Ait,aes(x=V1,y=V2,colour=group)) + 
  geom_point(size=1.5) + #绘制点图并设定大小
  theme_bw() +  #使用黑白主题
  scale_color_manual(values=c("#3F5BA7","#D72326"))+theme_bw()+ #使用黑白主题
  guides(color=guide_legend(title = NULL)) + #去除图例标题
  stat_ellipse(level = 0.8)+  #添加置信椭圆（可以删除掉）
  theme(axis.title.x = element_text(size=15,family="sans"), # 修改x轴标题文本的属性
        axis.title.y = element_text(size=15,family="sans",angle=90), # 修改y轴标题文本的属性
        axis.text.y=element_text(size=12,family="sans"), # 修改x轴刻度标签文本的属性
        axis.text.x=element_text(size=12,family="sans"), # 修改y轴刻度标签文本的属性
        panel.grid = element_blank()#隐藏网格线
  )

pcoa_plot_HC_PBC_Ait<-pcoa_plot_HC_PBC_Ait+xlab(paste("PCo1"," ",pc_importance_HC_PBC_Ait[1],"%",sep = ""))+ylab(paste("PCo2"," ",pc_importance_HC_PBC_Ait[2],"%",sep = ""))
pcoa_plot_HC_PBC_Ait
#2.PBC+HC gene Aitchison
pcoa_HC_PBC_Ait<- cmdscale(Ait_gene_abs[c(HC_new,PBCbase),c(HC_new,PBCbase)], k = nrow(metadata[c(HC_new,PBCbase),]) - 1, eig = TRUE, add = TRUE)

#Extract PCoA coordinates
pc12_HC_PBC_Ait<- pcoa_HC_PBC_Ait$points[,1:2]
pc12_HC_PBC_Ait=as.data.frame(pc12_HC_PBC_Ait)
pc_importance_HC_PBC_Ait<-round(pcoa_HC_PBC_Ait$eig/sum(pcoa_HC_PBC_Ait$eig)*100,digits = 2)
pc12_HC_PBC_Ait$SampleID<-row.names(pc12_HC_PBC_Ait)
data_HC_PBC_Ait<-merge(pc12_HC_PBC_Ait,metadata_PBCHC[c(HC_new,PBCbase),],by="SampleID")

#根据group
library(ggplot2)
pcoa_plot_HC_PBC_Ait<-ggplot(data_HC_PBC_Ait,aes(x=V1,y=V2,colour=cluster)) + 
  geom_point(size=1.5) + #绘制点图并设定大小
  theme_bw() +  #使用黑白主题
  scale_color_manual(values=c("#088288","#F19433","#3F5BA7"))+theme_bw()+ #使用黑白主题
  guides(color=guide_legend(title = NULL)) + #去除图例标题
  stat_ellipse(level = 0.8)+  #添加置信椭圆（可以删除掉）
  theme(axis.title.x = element_text(size=15,family="sans"), # 修改x轴标题文本的属性
        axis.title.y = element_text(size=15,family="sans",angle=90), # 修改y轴标题文本的属性
        axis.text.y=element_text(size=12,family="sans"), # 修改x轴刻度标签文本的属性
        axis.text.x=element_text(size=12,family="sans"), # 修改y轴刻度标签文本的属性
        panel.grid = element_blank()#隐藏网格线
  )

pcoa_plot_HC_PBC_Ait<-pcoa_plot_HC_PBC_Ait+xlab(paste("PCo1"," ",pc_importance_HC_PBC_Ait[1],"%",sep = ""))+ylab(paste("PCo2"," ",pc_importance_HC_PBC_Ait[2],"%",sep = ""))
pcoa_plot_HC_PBC_Ait

#3.PBCbase according to cluster **figure
#Plot the PCoA based on Aitchison distance among PBCbase and HC
pcoa_PBC_Ait<- cmdscale(Ait_genus_abs[PBCbase,PBCbase], k = nrow(metadata[PBCbase,]) - 1, eig = TRUE, add = TRUE)

#Extract PCoA coordinates
pc12_PBC_Ait<- pcoa_PBC_Ait$points[,1:2]
pc12_PBC_Ait=as.data.frame(pc12_PBC_Ait)
pc_importance_PBC_Ait<-round(pcoa_PBC_Ait$eig/sum(pcoa_PBC_Ait$eig)*100,digits = 2)
pc12_PBC_Ait$SampleID<-row.names(pc12_PBC_Ait)
data_PBC_Ait<-merge(pc12_PBC_Ait,metadata_PBCHC[PBCbase,],by="SampleID")

#根据cluster
pcoa_plot_PBC_Ait<-ggplot(data_PBC_Ait,aes(x=V1,y=V2,colour=cluster)) + 
  geom_point(size=1.5) + #绘制点图并设定大小
  theme_bw() +  #使用黑白主题
  scale_color_manual(values=c("#088288","#F19433"))+theme_bw()+ #使用黑白主题
  guides(color=guide_legend(title = NULL)) + #去除图例标题
  
  theme(axis.title.x = element_text(size=15,family="sans"), # 修改x轴标题文本的属性
        axis.title.y = element_text(size=15,family="sans",angle=90), # 修改y轴标题文本的属性
        axis.text.y=element_text(size=12,family="sans"), # 修改x轴刻度标签文本的属性
        axis.text.x=element_text(size=12,family="sans"), # 修改y轴刻度标签文本的属性
        panel.grid = element_blank()#隐藏网格线
  )

pcoa_plot_PBC_Ait<-pcoa_plot_PBC_Ait+xlab(paste("PCo1"," ",pc_importance_PBC_Ait[1],"%",sep = ""))+ylab(paste("PCo2"," ",pc_importance_PBC_Ait[2],"%",sep = ""))
pcoa_plot_PBC_Ait


#permanova analysis
library(vegan)

adonis(Ait_species_abs[c(PBCpair_6m,PBC_6m_pair),c(PBCpair_6m,PBC_6m_pair)]~group,data = metadata[c(PBCpair_6m,PBC_6m_pair),])

################compare distance among clt1/clt2 and HC#############
library(dplyr)
#species level
Ait_species_abs_clt1HC=Ait_species_abs[PBCbase_cluster1,HC] %>% as.numeric() %>% as.data.frame()
Ait_species_abs_clt2HC=Ait_species_abs[PBCbase_cluster2,HC] %>% as.numeric() %>% as.data.frame()
colnames(Ait_species_abs_clt1HC)="distance"
colnames(Ait_species_abs_clt2HC)="distance"
wilcox.test(Ait_species_abs_clt1HC$distance,Ait_species_abs_clt2HC$distance)
median(Ait_species_abs_clt1HC$distance)
median(Ait_species_abs_clt2HC$distance)
Ait_species_abs_PBC_HC=rbind(Ait_species_abs_clt1HC,Ait_species_abs_clt2HC)
Ait_species_abs_PBC_HC$group=c(rep("clt1HC",10800),rep("clt2HC",9000))
library(ggplot2)
ggplot(data = Ait_species_abs_PBC_HC,aes(x=group,y=distance,color=group))+
  geom_boxplot()+geom_jitter(width = 0.2,size=0.3)+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  scale_color_manual(values = c("#FFB733","#B75454"))+
  annotate("text",x=2, y=400,label= "P< 2.2e-16")
  
#genus level HC_new *figure
Ait_genus_abs_clt1HC=Ait_genus_abs[PBCbase_cluster1,HC_new] %>% as.numeric() %>% as.data.frame()
Ait_genus_abs_clt2HC=Ait_genus_abs[PBCbase_cluster2,HC_new] %>% as.numeric() %>% as.data.frame()
colnames(Ait_genus_abs_clt1HC)="distance"
colnames(Ait_genus_abs_clt2HC)="distance"
wilcox.test(Ait_genus_abs_clt1HC$distance,Ait_genus_abs_clt2HC$distance)
median(Ait_genus_abs_clt1HC$distance)
median(Ait_genus_abs_clt2HC$distance)
Ait_genus_abs_PBC_HC=rbind(Ait_genus_abs_clt1HC,Ait_genus_abs_clt2HC)
Ait_genus_abs_PBC_HC$group=c(rep("clt1HC",9301),rep("clt2HC",7991))
library(ggplot2)
ggplot(data = Ait_genus_abs_PBC_HC,aes(x=group,y=distance,color=group))+
  geom_boxplot()+geom_jitter(width = 0.2,size=0.3)+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  scale_color_manual(values = c("#088288","#F19433"))


#KEGG path3 level
Ait_KEGG_path3_abs_clt1HC=Ait_KEGG3_abs[PBCbase_cluster1,HC] %>% as.numeric() %>% as.data.frame()
Ait_KEGG_path3_abs_clt2HC=Ait_KEGG3_abs[PBCbase_cluster2,HC] %>% as.numeric() %>% as.data.frame()
colnames(Ait_KEGG_path3_abs_clt1HC)="distance"
colnames(Ait_KEGG_path3_abs_clt2HC)="distance"
wilcox.test(Ait_KEGG_path3_abs_clt1HC$distance,Ait_KEGG_path3_abs_clt2HC$distance)
median(Ait_KEGG_path3_abs_clt1HC$distance)
median(Ait_KEGG_path3_abs_clt2HC$distance)
Ait_KEGG_path3_abs_PBC_HC=rbind(Ait_KEGG_path3_abs_clt1HC,Ait_KEGG_path3_abs_clt2HC)
Ait_KEGG_path3_abs_PBC_HC$group=c(rep("clt1HC",10800),rep("clt2HC",9000))
library(ggplot2)
ggplot(data = Ait_KEGG_path3_abs_PBC_HC,aes(x=group,y=distance,color=group))+
  geom_boxplot()+geom_jitter(width = 0.2,size=0.3)+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  scale_color_manual(values = c("#FFB733","#B75454"))+
  annotate("text",x=2,y=45,label="P< 2.2e-16")

#output clutser membership data
setwd("D:/PBC multiomics data/analysis V6/metagenomics/data")
write.csv(metadata_PBCHC,"cluster.res_PBCHC.csv")

############PCoA colored based on genus abundance################
#1.only PBC sample,based on Bacteroides and Prevotella genus abundance

#Plot the PCoA based on Aitchison distance among PBCbase
pcoa_PBC_Ait<- cmdscale(Ait_genus_abs[PBCbase,PBCbase], k = nrow(metadata[PBCbase,]) - 1, eig = TRUE, add = TRUE)

#Extract PCoA coordinates
pc12_PBC_Ait<- pcoa_PBC_Ait$points[,1:2]
pc12_PBC_Ait=as.data.frame(pc12_PBC_Ait)
pc_importance_PBC_Ait<-round(pcoa_PBC_Ait$eig/sum(pcoa_PBC_Ait$eig)*100,digits = 2)
pc12_PBC_Ait$SampleID<-row.names(pc12_PBC_Ait)
data_PBC_Ait<-merge(pc12_PBC_Ait,metadata_PBCbase[PBCbase,],by="SampleID")

data_PBC_Ait$Bacteroides=genus_abs_clr[PBCbase,'k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides']
data_PBC_Ait$Veillonella=genus_abs_clr[PBCbase,"k__Bacteria;p__Firmicutes;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella"]
row.names(data_PBC_Ait)=data_PBC_Ait$SampleID


# *figure 
pcoa_plot_PBC_Ait<-ggplot(data_PBC_Ait,aes(x=V1,y=V2,colour=Response_6m)) + 
  geom_point(size=4) + #绘制点图并设定大小
  theme_bw() +  #使用黑白主题
  theme_bw()+ #使用黑白主题
  scale_color_manual(values=c("black","#A6DCE4"))+
  guides(color=guide_legend(title = NULL)) + #去除图例标题
  theme(axis.title.x = element_text(size=15,family="sans"), # 修改x轴标题文本的属性
        axis.title.y = element_text(size=15,family="sans",angle=90), # 修改y轴标题文本的属性
        axis.text.y=element_text(size=12,family="sans"), # 修改x轴刻度标签文本的属性
        axis.text.x=element_text(size=12,family="sans"), # 修改y轴刻度标签文本的属性
        panel.grid = element_blank()#隐藏网格线
  )

pcoa_plot_PBC_Ait<-pcoa_plot_PBC_Ait+xlab(paste("PCo1"," ",pc_importance_PBC_Ait[1],"%",sep = ""))+ylab(paste("PCo2"," ",pc_importance_PBC_Ait[2],"%",sep = ""))
pcoa_plot_PBC_Ait





################PCo1-taxa association *figure#################

wilcox.test(data_PBC_Ait[PBCbase_cluster1,2],data_PBC_Ait[PBCbase_cluster2,2])
#1.genus
colnames(data_PBC_Ait)
row.names(data_PBC_Ait)=data_PBC_Ait$SampleID
PBC_PCo1_genus=cbind(data_PBC_Ait[PBCbase,c(1,2,3)],genus_abs_clr[PBCbase,row.names(genus_rela_filter)])

corR_PCo1_genus=sapply(2, function(x) sapply(4:886,function(y){cor.test(PBC_PCo1_genus[,x],PBC_PCo1_genus[,y],method = "spearman")$estimate}))
corP_PCo1_genus=sapply(2, function(x) sapply(4:886,function(y){cor.test(PBC_PCo1_genus[,x],PBC_PCo1_genus[,y],method = "spearman")$p.value}))
corq_PCo1_genus=p.adjust(corP_PCo1_genus,"fdr") %>% matrix(.,ncol = 1)
row.names(corR_PCo1_genus)=colnames(PBC_PCo1_genus)[4:886]
row.names(corP_PCo1_genus)=colnames(PBC_PCo1_genus)[4:886]
row.names(corq_PCo1_genus)=colnames(PBC_PCo1_genus)[4:886]

genus_mean_PBCbase=apply(genus_rela_filter[,PBCbase], 1, mean) %>% as.data.frame(.)
row.names(genus_mean_PBCbase)=row.names(genus_rela_filter)
colSums(genus_rela_filter)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/PCo1_correlation")
write.csv(corP_PCo1_genus,"corP_PCo1_genus.csv")
write.csv(corq_PCo1_genus,"corq_PCo1_genus.csv")
write.csv(corR_PCo1_genus,"corR_PCo1_genus.csv")
write.csv(genus_mean_PBCbase,"genus_mean_PBCbase.csv")

#2.order
PBC_PCo1_order=cbind(data_PBC_Ait[PBCbase,c(1,2,3)],order_abs_clr[PBCbase,row.names(order_rela_filter)])

corR_PCo1_order=sapply(2, function(x) sapply(4:206,function(y){cor.test(PBC_PCo1_order[,x],PBC_PCo1_order[,y],method = "spearman")$estimate}))
corP_PCo1_order=sapply(2, function(x) sapply(4:206,function(y){cor.test(PBC_PCo1_order[,x],PBC_PCo1_order[,y],method = "spearman")$p.value}))
corq_PCo1_order=p.adjust(corP_PCo1_order,"fdr") %>% matrix(.,ncol = 1)
row.names(corR_PCo1_order)=colnames(PBC_PCo1_order)[4:206]
row.names(corP_PCo1_order)=colnames(PBC_PCo1_order)[4:206]
row.names(corq_PCo1_order)=colnames(PBC_PCo1_order)[4:206]

order_mean_PBCbase=apply(order_rela_filter[,PBCbase], 1, mean) %>% as.data.frame(.)
row.names(order_mean_PBCbase)=row.names(order_rela_filter)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/order/PCo1_correlation")
write.csv(corP_PCo1_order,"corP_PCo1_order.csv")
write.csv(corq_PCo1_order,"corq_PCo1_order.csv")
write.csv(corR_PCo1_order,"corR_PCo1_order.csv")
write.csv(order_mean_PBCbase,"order_mean_PBCbase.csv")
data_PBC_Ait$Bacteroidales_order=order_abs_clr[PBCbase,"k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales"]
data_PBC_Ait$Veillonellales_order=order_abs_clr[PBCbase,"k__Bacteria;p__Firmicutes;c__Negativicutes;o__Veillonellales"]
data_PBC_Ait$Clostridiales_order=order_abs_clr[PBCbase,"k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales"]

#3.species
PBC_PCo1_species=cbind(data_PBC_Ait[PBCbase,c(1,2,3)],species_abs_clr[PBCbase,row.names(species_rela_filter)])

corR_PCo1_species=sapply(2, function(x) sapply(4:3857,function(y){cor.test(PBC_PCo1_species[,x],PBC_PCo1_species[,y],method = "spearman")$estimate}))
corP_PCo1_species=sapply(2, function(x) sapply(4:3857,function(y){cor.test(PBC_PCo1_species[,x],PBC_PCo1_species[,y],method = "spearman")$p.value}))
corq_PCo1_species=p.adjust(corP_PCo1_species,"fdr") %>% matrix(.,ncol = 1)
row.names(corR_PCo1_species)=colnames(PBC_PCo1_species)[4:3857]
row.names(corP_PCo1_species)=colnames(PBC_PCo1_species)[4:3857]
row.names(corq_PCo1_species)=colnames(PBC_PCo1_species)[4:3857]

species_mean_PBCbase=apply(species_rela_filter[,PBCbase],1,mean) %>% as.data.frame(.)
row.names(species_mean_PBCbase)=row.names(species_rela_filter)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/PCo1_correlation")
write.csv(corP_PCo1_species,"corP_PCo1_species.csv")
write.csv(corq_PCo1_species,"corq_PCo1_species.csv")
write.csv(corR_PCo1_species,"corR_PCo1_species.csv")
write.csv(species_mean_PBCbase,"species_mean_PBCbase.csv")

#4.family
PBC_PCo1_family=cbind(data_PBC_Ait[PBCbase,c(1,2,3)],family_abs_clr[PBCbase,row.names(family_rela_filter)])

corR_PCo1_family=sapply(2, function(x) sapply(4:336,function(y){cor.test(PBC_PCo1_family[,x],PBC_PCo1_family[,y],method = "spearman")$estimate}))
corP_PCo1_family=sapply(2, function(x) sapply(4:336,function(y){cor.test(PBC_PCo1_family[,x],PBC_PCo1_family[,y],method = "spearman")$p.value}))
corq_PCo1_family=p.adjust(corP_PCo1_family,"fdr") %>% matrix(.,ncol = 1)
row.names(corR_PCo1_family)=colnames(PBC_PCo1_family)[4:336]
row.names(corP_PCo1_family)=colnames(PBC_PCo1_family)[4:336]
row.names(corq_PCo1_family)=colnames(PBC_PCo1_family)[4:336]

family_mean_PBCbase=apply(family_rela_filter[,PBCbase],1,mean) %>% as.data.frame(.)
row.names(family_mean_PBCbase)=row.names(family_rela_filter)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/family/PCo1_correlation")
write.csv(corP_PCo1_family,"corP_PCo1_family.csv")
write.csv(corq_PCo1_family,"corq_PCo1_family.csv")
write.csv(corR_PCo1_family,"corR_PCo1_family.csv")
write.csv(family_mean_PBCbase,"family_mean_PBCbase.csv")

###figure** 根据PCo1相关性可视化driving taxa or ratio
data_PBC_Ait$Oscillibacter=genus_abs_clr[PBCbase,"k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Oscillospiraceae;g__Oscillibacter"]
data_PBC_Ait$BO_ratio=data_PBC_Ait$Bacteroides-data_PBC_Ait$Oscillibacter
data_PBC_Ait$VO_ratio=data_PBC_Ait$Veillonella-data_PBC_Ait$Oscillibacter
data_PBC_Ait$Oscillospiraceae=family_abs_clr[PBCbase,"k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Oscillospiraceae"]
data_PBC_Ait$Bacteroidaceae=family_abs_clr[PBCbase,"k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae"]
data_PBC_Ait$BO_ratio_family=data_PBC_Ait$Bacteroidaceae-data_PBC_Ait$Oscillospiraceae
data_PBC_Ait$BO_ratio_family=data_PBC_Ait$Bacteroidaceae-data_PBC_Ait$Oscillospiraceae
data_PBC_Ait$Ruminococcaceae=family_abs_clr[PBCbase,"k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae"]
data_PBC_Ait$Ruminococcaceae.spp=genus_abs_clr[PBCbase,"k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__unidentified"]
data_PBC_Ait$Butyricicoccus=genus_abs_clr[PBCbase,"k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Butyricicoccus"]
data_PBC_Ait$Clostridium_leptum=species_abs_clr[PBCbase,"k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__unidentified;s__[Clostridium] leptum"]
data_PBC_Ait$Oscillibacter.spp=species_abs_clr[PBCbase,"k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Oscillospiraceae;g__Oscillibacter;s__Oscillibacter sp."]
data_PBC_Ait$I.butyriciproducens=species_abs_clr[PBCbase,"k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__unidentified;g__Intestinimonas;s__Intestinimonas butyriciproducens"]
data_PBC_Ait$V.dispa=species_abs_clr[PBCbase,"k__Bacteria;p__Firmicutes;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella;s__Veillonella dispar"]
data_PBC_Ait$B.fragilis=species_abs_clr[PBCbase,"k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides fragilis"]
data_PBC_Ait$Ruminococcus_563=species_abs_clr[row.names(data_PBC_Ait),"k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminococcus;s__Ruminococcus sp. CAG:563"]

#配色V1
ggplot(data_PBC_Ait,aes(x=V1,y=V2,colour=B.fragilis)) + 
  geom_point(size=4) + 
  scale_colour_gradientn(colours = c("#F7E120","#BBD631","#80C556","#1CA485","#34608E","#461D5F","#46195B"))+
  theme_bw()  +
  xlab(paste("PCo1"," ",pc_importance_PBC_Ait[1],"%",sep = ""))+ylab(paste("PCo2"," ",pc_importance_PBC_Ait[2],"%",sep = ""))

#配色V2
ggplot(data_PBC_Ait,aes(x=V1,y=V2,colour=BO_ratio_family)) + 
  geom_point(size=4) + 
  scale_colour_gradient2(low = "#4A86BD",mid="white",high ="#E62621")+
  theme_bw() +
  xlab(paste("PCo1"," ",pc_importance_PBC_Ait[1],"%",sep = ""))+ylab(paste("PCo2"," ",pc_importance_PBC_Ait[2],"%",sep = ""))

#clt1 clt2 PCo1 loading boxplot
library(ggplot2)
ggplot(data = data_PBC_Ait,aes(x=cluster,y=V1,fill=cluster))+
  geom_boxplot(color="grey")+
  theme(axis.line = element_line(colour = "grey"),panel.background = element_blank(),legend.position = "none")+
  scale_fill_manual(values = c("#088288","#F19433"))

ggplot(data = data_PBC_Ait,aes(x=cluster,y=V2,fill=cluster))+
  geom_boxplot(color="grey")+
  theme(axis.line = element_line(colour = "grey"),panel.background = element_blank(),legend.position = "none")+
  scale_fill_manual(values = c("#088288","#F19433"))
wilcox.test(data_PBC_Ait[PBCbase_cluster1,"V2"],data_PBC_Ait[PBCbase_cluster2,"V2"])

############沿着V1升序绘制细菌每个样本的柱状图#############
data_PBC_Ait_reorder <- data_PBC_Ait[order(data_PBC_Ait$V1),]
data_PBC_Ait_reorder$SampleID=factor(data_PBC_Ait_reorder$SampleID,levels= unique(data_PBC_Ait_reorder$SampleID))

ggplot(data= data_PBC_Ait_reorder,aes(x=SampleID,y=Ruminococcus_563))+geom_bar(mapping = NULL, data = NULL, stat = "identity",
                                                                                 width=0.9, position="dodge",fill="light blue")

setwd("D:/PBC multiomics data/fecal microbiota transplantation")
write.csv(data_PBC_Ait_reorder,"data_PBC_Ait_reorder.csv")
#绘制Bacteroides与Prevotella单独的变化图
ggplot(data = data_PBC_Ait_reorder,aes(x=cluster,y=Bacteroides,fill=cluster))+geom_boxplot()+
  scale_fill_manual(values=c("#FFB733","#B75454"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")

ggplot(data = data_PBC_Ait_reorder,aes(x=cluster,y=Prevotella,fill=cluster))+geom_boxplot()+
  scale_fill_manual(values=c("#FFB733","#B75454"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")


ggplot(data = data_PBC_Ait_reorder,aes(x=cluster,y=BF_ratio))+geom_boxplot()

#绘制与PCo1相关性>0.8的33个abundant species (mean relative abundance>0.0001)的corR barplot **figure 
#为了体现Clostridial是驱动cluster的taxa
species_mean_PBCbase_0.0001
corR_PCo1_species=as.data.frame(corR_PCo1_species)
colnames(corR_PCo1_species)="corR"
library(dplyr)
PCo1_corR0.8_species=filter(corR_PCo1_species,corR< -0.8) %>% row.names(.)
PCo1_corR0.8_species=intersect(PCo1_corR0.8_species,species_mean_PBCbase_0.0001)

newname=vector()
name=strsplit(PCo1_corR0.8_species,"s__")
for (i in 1:33) {
  newname[i]=name[[i]][2]
}

corR_PCo1_species_0.8=corR_PCo1_species[PCo1_corR0.8_species,1] %>% as.data.frame(.)
row.names(corR_PCo1_species_0.8)=PCo1_corR0.8_species
corR_PCo1_species_0.8$species=newname
colnames(corR_PCo1_species_0.8)[1]="corR"
corR_PCo1_species_0.8=arrange(corR_PCo1_species_0.8,desc(corR_PCo1_species_0.8$corR))
corR_PCo1_species_0.8$species=factor(corR_PCo1_species_0.8$species,levels= unique(corR_PCo1_species_0.8$species))
library(ggplot2)
ggplot(corR_PCo1_species_0.8,aes(x=species,y=-corR))+
  geom_bar(stat = "identity",fill="#D36D2D",width = 0.6)+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  coord_flip()

########Metagenome centroid change after treatment *figure############
#6month HC_new
library(vegan)
metadata_centroid=c(PBCpair_6m_cluster1,PBCpair_6m_cluster2,PBC_6m_cluster1,PBC_6m_cluster2,HC_new) %>% as.data.frame()
metadata_centroid$group=c(rep("clt1.base",36),rep("clt2.base",23),rep("clt1.6m",36),rep("clt2.6m",23),rep("HC",131))
metadata_centroid$group=as.factor(metadata_centroid$group)
Ait_genus_abs_centroid=vegdist(genus_abs_clr[c(PBCpair_6m_cluster1,PBCpair_6m_cluster2,PBC_6m_cluster1,PBC_6m_cluster2,HC_new),row.names(genus_rela_filter)],method = "euclidean") 
Ait_species_abs_centroid=vegdist(species_abs_clr[c(PBCpair_6m_cluster1,PBCpair_6m_cluster2,PBC_6m_cluster1,PBC_6m_cluster2,HC_new),row.names(species_rela_filter)],method = "euclidean")
Ait_KEGG_path3_abs_centroid=vegdist(KEGG_path3_abs_clr[c(PBCpair_6m_cluster1,PBCpair_6m_cluster2,PBC_6m_cluster1,PBC_6m_cluster2,HC_new),row.names(KEGG_path3_rela_filter)],method = "euclidean")
Ait_KO_gene_abs_centroid=vegdist(KO_gene_abs_clr[c(PBCpair_6m_cluster1,PBCpair_6m_cluster2,PBC_6m_cluster1,PBC_6m_cluster2,HC_new),row.names(KO_gene_rela_filter)],method = "euclidean")
Ait_GMM_abs_centroid=vegdist(GMM_abs_clr[c(PBCpair_6m_cluster1,PBCpair_6m_cluster2,PBC_6m_cluster1,PBC_6m_cluster2,HC_new),row.names(GMM_rela_filter)],method = "euclidean")
#1.根据cluster显示
mod=betadisper(Ait_species_abs_centroid,metadata_centroid$group)
df.pcoa = data.frame(x=mod$vectors[,1],y=mod$vectors[,2],group=mod$group)
centroids = data.frame(x=mod$centroids[,1],y=mod$centroids[,2],group = factor(rownames(mod$centroids)))
gg <- merge(df.pcoa,centroids, by ='group', suffixes = c('','.centroid')) 
gg$group  = factor(gg$group,levels=c('clt1.base','clt1.6m','clt2.base','clt2.6m','HC'))
library(ggplot2)
ppcoaUDCA_ctrl <- ggplot(gg) +
  geom_point(aes(x=x,y=y,color=group),size=1,alpha=0.4) + scale_color_manual(values = c('#088288','#ADD5DB','#F19433','#FFD966','#3F5BA7')) +
  geom_segment(aes(x=x.centroid, y =y.centroid, xend=x,yend=y,color=group),alpha=0.2) +
   geom_point(data=centroids,aes(x=x,y=y,color=group),size=3)+
  xlab(paste('PCoA1 (',100*round(mod$eig[1] / sum(mod$eig),2),'%)',sep='')) +
  ylab(paste('PCoA2 (',100*round(mod$eig[2] / sum(mod$eig),2),'%)',sep='')) + 
  theme_bw()+
  theme(legend.position = 'none', panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  ppcoaUDCA_ctrl
  

############比较cluster1/2治疗前后与HC的距离(species level) *figure############  
#1.1 cluster1治疗前后与HC的距离
  Ait_species_abs_centroid_matrix=as.matrix(Ait_species_abs_centroid) #将距离数据转化为矩阵
  Ait_species_abs_clt1.base.HC=Ait_species_abs_centroid_matrix[PBCpair_6m_cluster1,HC_new] %>% as.numeric(.) %>% as.data.frame(.)
  Ait_species_abs_clt1.treat.HC=Ait_species_abs_centroid_matrix[PBC_6m_cluster1,HC_new] %>% as.numeric(.) %>% as.data.frame(.)
  colnames(Ait_species_abs_clt1.base.HC)="distance"
  colnames(Ait_species_abs_clt1.treat.HC)="distance"
  wilcox.test(Ait_species_abs_clt1.base.HC$distance,Ait_species_abs_clt1.treat.HC$distance) #P<2.2e-16
  mean(Ait_species_abs_clt1.base.HC$distance) #mean=215.4802
  mean(Ait_species_abs_clt1.treat.HC$distance) #mean=221.3536
  Ait_species_abs_clt1_HC=rbind(Ait_species_abs_clt1.base.HC,Ait_species_abs_clt1.treat.HC)
  Ait_species_abs_clt1_HC$group=c(rep("clt1.base.HC",4716),rep("clt1.treat.HC",4716))
  library(ggplot2)
  ggplot(data = Ait_species_abs_clt1_HC,aes(x=group,y=distance,color=group))+
    geom_boxplot()+geom_jitter(width = 0.2,size=0.3)+
    theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
    scale_color_manual(values = c("#088288","#ADD5DB"))
#1.2 cluster1同一个体治疗前后的距离
  Ait_species_abs_clt1.base.treatment=Ait_species_abs_centroid_matrix[PBCpair_6m_cluster1,PBC_6m_cluster1]  
  clt1.base.treat=diag(Ait_species_abs_clt1.base.treatment) %>% as.data.frame(.)
  colnames(clt1.base.treat)="distance"
  
#1.3 cluster1治疗前后的距离，不论是否为同一个体
  clt1.base.treat_all=Ait_species_abs_clt1.base.treatment[PBCpair_6m_cluster1,PBC_6m_cluster1] %>% as.numeric() %>% as.data.frame()
  colnames(clt1.base.treat_all)="distance"  
  
#2.1 cluster2治疗前后与HC的距离
  Ait_species_abs_clt2.base.HC=Ait_species_abs_centroid_matrix[PBCpair_6m_cluster2,HC_new] %>% as.numeric() %>% as.data.frame()
  Ait_species_abs_clt2.treat.HC=Ait_species_abs_centroid_matrix[PBC_6m_cluster2,HC_new] %>% as.numeric() %>% as.data.frame()
  colnames(Ait_species_abs_clt2.base.HC)="distance" #mean=241.4884
  colnames(Ait_species_abs_clt2.treat.HC)="distance" #mean=243.0914
  wilcox.test(Ait_species_abs_clt2.base.HC$distance,Ait_species_abs_clt2.treat.HC$distance) #P=0.00232
  mean(Ait_species_abs_clt2.base.HC$distance)
  mean(Ait_species_abs_clt2.treat.HC$distance)
  Ait_species_abs_clt2_HC=rbind(Ait_species_abs_clt2.base.HC,Ait_species_abs_clt2.treat.HC)
  Ait_species_abs_clt2_HC$group=c(rep("clt2.base.HC",3013),rep("clt2.treat.HC",3013))
  library(ggplot2)
  ggplot(data = Ait_species_abs_clt2_HC,aes(x=group,y=distance,color=group))+
    geom_boxplot()+geom_jitter(width = 0.2,size=0.3)+
    theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
    scale_color_manual(values = c("#F19433","#FFD966"))
#2.2 cluster2同一个体治疗前后的距离
Ait_species_abs_clt2.base.treatment=Ait_species_abs_centroid_matrix[PBCpair_6m_cluster2,PBC_6m_cluster2]  
clt2.base.treat=diag(Ait_species_abs_clt2.base.treatment) %>%as.data.frame(.)
colnames(clt2.base.treat)="distance"  

#2.3 cluster2治疗前后的距离，不论是否为同一个体
  clt2.base.treat_all=Ait_species_abs_clt2.base.treatment[PBCpair_6m_cluster2,PBC_6m_cluster2] %>% as.numeric() %>% as.data.frame()
  colnames(clt2.base.treat_all)="distance"  
  
#绘制不同cluster同一样本治疗前后距离变化图 
  clt.base.treat=rbind(clt1.base.treat,clt2.base.treat)
  clt.base.treat$group=c(rep("clt1",36),rep("clt2",23))
  ggplot(data = clt.base.treat,aes(x=group,y=distance,color=group))+
    geom_boxplot()+geom_jitter(width = 0.2,size=0.8)+
    theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
    scale_color_manual(values = c("#088288","#F19433"))
  wilcox.test(clt1.base.treat$distance,clt2.base.treat$distance)
#绘制不同cluster所有治疗前后距离变化图 
  clt.base.treat_all=rbind(clt1.base.treat_all,clt2.base.treat_all)
  clt.base.treat_all$group=c(rep("clt1",1296),rep("clt2",529))
  ggplot(data = clt.base.treat_all,aes(x=group,y=distance,color=group))+
    geom_boxplot()+geom_jitter(width = 0.2,size=0.8)+
    theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
    scale_color_manual(values = c("#088288","#F19433"))

 
   
  
  
 #3.将各组治疗前后与HC的距离合并且绘图
  Ait_species_abs_clt_HC=rbind(Ait_species_abs_clt1_HC,Ait_species_abs_clt2_HC)
  ggplot(data = Ait_species_abs_clt_HC,aes(x=group,y=distance,color=group))+
    geom_boxplot()+geom_jitter(width = 0.2,size=0.8)+
    theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
    scale_color_manual(values = c("#088288","#ADD5DB","#F19433","#FFD966"))

  

  ############比较cluster1/2治疗前后与HC的距离(KO gene level)############ 
  #1.1 cluster1治疗前后与HC的距离
  Ait_KO_gene_abs_centroid_matrix=as.matrix(Ait_KO_gene_abs_centroid) #将距离数据转化为矩阵
  Ait_KO_gene_abs_clt1.base.HC=Ait_KO_gene_abs_centroid_matrix[PBCpair_6m_cluster1,HC_new] %>% as.numeric() %>% as.data.frame()
  Ait_KO_gene_abs_clt1.treat.HC=Ait_KO_gene_abs_centroid_matrix[PBC_6m_cluster1,HC_new] %>% as.numeric() %>% as.data.frame()
  colnames(Ait_KO_gene_abs_clt1.base.HC)="distance"
  colnames(Ait_KO_gene_abs_clt1.treat.HC)="distance"
  wilcox.test(Ait_KO_gene_abs_clt1.base.HC$distance,Ait_KO_gene_abs_clt1.treat.HC$distance) #P<2.2e-16
  mean(Ait_KO_gene_abs_clt1.base.HC$distance) #mean=252.2408
  mean(Ait_KO_gene_abs_clt1.treat.HC$distance) #mean=247.6397
  Ait_KO_gene_abs_clt1_HC=rbind(Ait_KO_gene_abs_clt1.base.HC,Ait_KO_gene_abs_clt1.treat.HC)
  Ait_KO_gene_abs_clt1_HC$group=c(rep("clt1.base.HC",4716),rep("clt1.treat.HC",4716))
  library(ggplot2)
  ggplot(data = Ait_KO_gene_abs_clt1_HC,aes(x=group,y=distance,color=group))+
    geom_boxplot()+geom_jitter(width = 0.2,size=0.3)+
    theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
    scale_color_manual(values = c("#088288","#ADD5DB"))
  #1.2 cluster1治疗前后相同样本的距离
  Ait_KO_gene_abs_clt1.base.treatment=Ait_KO_gene_abs_centroid_matrix[PBCpair_6m_cluster1,PBC_6m_cluster1]  
  clt1.base.treat_KO_gene=diag(Ait_KO_gene_abs_clt1.base.treatment) %>% as.data.frame(.)
  colnames(clt1.base.treat_KO_gene)="distance"
  #2.1 cluster2治疗前后与HC的距离
  Ait_KO_gene_abs_clt2.base.HC=Ait_KO_gene_abs_centroid_matrix[PBCpair_6m_cluster2,HC_new] %>% as.numeric() %>% as.data.frame()
  Ait_KO_gene_abs_clt2.treat.HC=Ait_KO_gene_abs_centroid_matrix[PBC_6m_cluster2,HC_new] %>% as.numeric() %>% as.data.frame()
  colnames(Ait_KO_gene_abs_clt2.base.HC)="distance" #mean=261.0791
  colnames(Ait_KO_gene_abs_clt2.treat.HC)="distance" #mean=255.0801
  wilcox.test(Ait_KO_gene_abs_clt2.base.HC$distance,Ait_KO_gene_abs_clt2.treat.HC$distance) #P=0.00232
  mean(Ait_KO_gene_abs_clt2.base.HC$distance)
  mean(Ait_KO_gene_abs_clt2.treat.HC$distance)
  Ait_KO_gene_abs_clt2_HC=rbind(Ait_KO_gene_abs_clt2.base.HC,Ait_KO_gene_abs_clt2.treat.HC)
  Ait_KO_gene_abs_clt2_HC$group=c(rep("clt2.base.HC",3013),rep("clt2.treat.HC",3013))
  library(ggplot2)
  ggplot(data = Ait_KO_gene_abs_clt2_HC,aes(x=group,y=distance,color=group))+
    geom_boxplot()+geom_jitter(width = 0.2,size=0.3)+
    theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
    scale_color_manual(values = c("#F19433","#FFD966"))
  #2.2 cluster2治疗前后相同样本的距离
  Ait_KO_gene_abs_clt2.base.treatment=Ait_KO_gene_abs_centroid_matrix[PBCpair_6m_cluster2,PBC_6m_cluster2]  
  clt2.base.treat_KO_gene=diag(Ait_KO_gene_abs_clt2.base.treatment) %>%as.data.frame(.)
  colnames(clt2.base.treat_KO_gene)="distance"  
  #绘制各组治疗前后变化图 
  clt.base.treat_KO_gene=rbind(clt1.base.treat_KO_gene,clt2.base.treat_KO_gene)
  clt.base.treat_KO_gene$group=c(rep("clt1",36),rep("clt2",23))
  ggplot(data = clt.base.treat_KO_gene,aes(x=group,y=distance,color=group))+
    geom_boxplot()+geom_jitter(width = 0.2,size=0.8)+
    theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
    scale_color_manual(values = c("#088288","#F19433"))
  wilcox.test(clt1.base.treat_KO_gene$distance,clt2.base.treat_KO_gene$distance)
  
  
  
  
  #1.根据group显示,即显示PBC总体治疗前后的变化
  library(vegan)
  metadata_centroid$group_all=c(rep("PBCbase",59),rep("PBC.6m",59),rep("HC",131))
  metadata_centroid$group_all=as.factor(metadata_centroid$group_all)
  mod=betadisper(Ait_species_abs_centroid,metadata_centroid$group_all)
  df.pcoa = data.frame(x=mod$vectors[,1],y=mod$vectors[,2],group=mod$group)
  centroids = data.frame(x=mod$centroids[,1],y=mod$centroids[,2],group = factor(rownames(mod$centroids)))
  gg <- merge(df.pcoa,centroids, by ='group', suffixes = c('','.centroid')) 
  gg$group  = factor(gg$group ,levels=c('PBCbase','PBC.6m','HC'))
  library(ggplot2)
  ppcoaUDCA_ctrl <- ggplot(gg) +
    geom_point(aes(x=x,y=y,color=group),size=1,alpha=0.4) + scale_color_manual(values = c('#D72326','#F29E9C','#3F5BA7')) +
    geom_segment(aes(x=x.centroid, y =y.centroid, xend=x,yend=y,color=group),alpha=0.2) +
    geom_point(data=centroids,aes(x=x,y=y,color=group),size=3)+
    xlab(paste('PCoA1 (',100*round(mod$eig[1] / sum(mod$eig),2),'%)',sep='')) +
    ylab(paste('PCoA2 (',100*round(mod$eig[2] / sum(mod$eig),2),'%)',sep='')) + 
    theme_bw()+
    theme( legend.position = "none",panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  ppcoaUDCA_ctrl

#############治疗前后PCoA图，无HC####################  
#Aitchison distance
library(robCompositions)
library(zCompositions)
library(compositions)
library(vegan)
#Distance calculation
Ait_species_abs=vegdist(species_abs_clr[,row.names(species_rela_filter)],method = "euclidean")%>% as.matrix(.)
  
#1.clt1
#Plot the PCoA based on Aitchison distance 
pcoa_clt1.base.treat_Ait<- cmdscale(Ait_species_abs[c(PBCpair_6m_cluster1,PBC_6m_cluster1),c(PBCpair_6m_cluster1,PBC_6m_cluster1)], k = nrow(metadata[c(PBCpair_6m_cluster1,PBC_6m_cluster1),]) - 1, eig = TRUE, add = TRUE)
  
#Extract PCoA coordinates
  pc12_clt1.base.treat_Ait<- pcoa_clt1.base.treat_Ait$points[,1:2]
  pc12_clt1.base.treat_Ait=as.data.frame(pc12_clt1.base.treat_Ait)
  pc_importance_clt1.base.treat_Ait<-round(pcoa_clt1.base.treat_Ait$eig/sum(pcoa_clt1.base.treat_Ait$eig)*100,digits = 2)
  pc12_clt1.base.treat_Ait$SampleID<-row.names(pc12_clt1.base.treat_Ait)
  data_clt1.base.treat_Ait<-merge(pc12_clt1.base.treat_Ait,metadata[c(PBCpair_6m_cluster1,PBC_6m_cluster1),],by="SampleID")
  data_clt1.base.treat_Ait$group= data_clt1.base.treat_Ait$group[,drop=TRUE]
  #根据group
  pcoa_plot_clt1.base.treat_Ait<-ggplot(data_clt1.base.treat_Ait,aes(x=V1,y=V2,colour=group)) + 
    geom_point(size=1) + #绘制点图并设定大小
    theme_bw() +  #使用黑白主题
    scale_color_manual(values=c("#ADD5DB","#088288"))+theme_bw()+ #使用黑白主题
    guides(color=guide_legend(title = NULL)) + #去除图例标题
    geom_line(aes(group=subjectID),color="#3E3A39")+
    theme(axis.title.x = element_text(size=15,family="sans"), # 修改x轴标题文本的属性
          axis.title.y = element_text(size=15,family="sans",angle=90), # 修改y轴标题文本的属性
          axis.text.y=element_text(size=12,family="sans"), # 修改x轴刻度标签文本的属性
          axis.text.x=element_text(size=12,family="sans"), # 修改y轴刻度标签文本的属性
          panel.grid = element_blank()#隐藏网格线
    )
  pcoa_plot_clt1.base.treat_Ait<-pcoa_plot_clt1.base.treat_Ait+xlab(paste("PCo1"," ",pc_importance_clt1.base.treat_Ait[1],"%",sep = ""))+ylab(paste("PCo2"," ",pc_importance_clt1.base.treat_Ait[2],"%",sep = ""))
  pcoa_plot_clt1.base.treat_Ait
  
  #1.clt2
  #Plot the PCoA based on Aitchison distance 
  pcoa_clt2.base.treat_Ait<- cmdscale(Ait_species_abs[c(PBCpair_6m_cluster2,PBC_6m_cluster2),c(PBCpair_6m_cluster2,PBC_6m_cluster2)], k = nrow(metadata[c(PBCpair_6m_cluster2,PBC_6m_cluster2),]) - 1, eig = TRUE, add = TRUE)
  
  #Extract PCoA coordinates
  pc12_clt2.base.treat_Ait<- pcoa_clt2.base.treat_Ait$points[,1:2]
  pc12_clt2.base.treat_Ait=as.data.frame(pc12_clt2.base.treat_Ait)
  pc_importance_clt2.base.treat_Ait<-round(pcoa_clt2.base.treat_Ait$eig/sum(pcoa_clt2.base.treat_Ait$eig)*100,digits = 2)
  pc12_clt2.base.treat_Ait$SampleID<-row.names(pc12_clt2.base.treat_Ait)
  data_clt2.base.treat_Ait<-merge(pc12_clt2.base.treat_Ait,metadata[c(PBCpair_6m_cluster2,PBC_6m_cluster2),],by="SampleID")
  data_clt2.base.treat_Ait$group= data_clt2.base.treat_Ait$group[,drop=TRUE]
  #根据group
  pcoa_plot_clt2.base.treat_Ait<-ggplot(data_clt2.base.treat_Ait,aes(x=V1,y=V2,colour=group)) + 
    geom_point(size=1) + #绘制点图并设定大小
    theme_bw() +  #使用黑白主题
    scale_color_manual(values=c("#FFD966","#F19433"))+theme_bw()+ #使用黑白主题
    guides(color=guide_legend(title = NULL)) + #去除图例标题
    geom_line(aes(group=subjectID),color="#3E3A39")+
    theme(axis.title.x = element_text(size=15,family="sans"), # 修改x轴标题文本的属性
          axis.title.y = element_text(size=15,family="sans",angle=90), # 修改y轴标题文本的属性
          axis.text.y=element_text(size=12,family="sans"), # 修改x轴刻度标签文本的属性
          axis.text.x=element_text(size=12,family="sans"), # 修改y轴刻度标签文本的属性
          panel.grid = element_blank()#隐藏网格线
    )
  pcoa_plot_clt2.base.treat_Ait<-pcoa_plot_clt2.base.treat_Ait+xlab(paste("PCo1"," ",pc_importance_clt2.base.treat_Ait[1],"%",sep = ""))+ylab(paste("PCo2"," ",pc_importance_clt2.base.treat_Ait[2],"%",sep = ""))
  pcoa_plot_clt2.base.treat_Ait
  
#############MaAsin 结果读入#################
#1.PBCHC
#1.1 species HC_new
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin")
maaslin_species_PBCHC=read.csv("maaslin_species PBCHC_new.csv",header = T,row.names = 1)
row.names(maaslin_species_PBCHC)=maaslin_species_PBCHC$species
maaslin_species_PBCHC_dif=filter(maaslin_species_PBCHC,qval<0.05) %>% row.names(.) #HC_new 
#1.2 genus HC_new
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin")
maaslin_genus_PBCHC=read.csv("maaslin_genus PBCHC_new.csv",header = T,row.names = 1)
row.names(maaslin_genus_PBCHC)=maaslin_genus_PBCHC$genus
maaslin_genus_PBCHC_dif=filter(maaslin_genus_PBCHC,qval<0.05) %>% row.names(.)
#1.3 KEGG_pathway HC_new
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin")
maaslin_KEGG_path3_PBCHC=read.csv("maaslin_KEGG_path3 PBCHC_new.csv",header = T,row.names = 1)
row.names(maaslin_KEGG_path3_PBCHC)=maaslin_KEGG_path3_PBCHC$Pathway
maaslin_KEGG_path3_PBCHC_dif=filter(maaslin_KEGG_path3_PBCHC,qval<0.2) %>% row.names(.)
#1.4 KO_gene HC_new
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin")
maaslin_KO_gene_PBCHC=read.csv("maaslin_KO_gene PBCHC_new.csv",header = T,row.names = 1)
maaslin_KO_gene_PBCHC_dif=filter(maaslin_KO_gene_PBCHC,qval<0.05) %>% row.names(.)
#1.5 GMM HC_new
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin")
maaslin_GMM_PBCHC=read.csv("maaslin_GMM PBCHC_new.csv",header = T,row.names = 1)
row.names(maaslin_GMM_PBCHC)=maaslin_GMM_PBCHC$feature
maaslin_GMM_PBCHC_dif=filter(maaslin_GMM_PBCHC,qval<0.2) %>% row.names(.)


#在MaAsin and paired analysis script.R中进行paired wilcox 分析
#输出结果
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/paired/PBC6m")
write.csv(res_KO_gene_pair.wilcox,"res_KO_gene_pair.wilcox.csv")
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/paired/PBC6m")
write.csv(res_KEGG_path3_pair.wilcox,"res_KEGG_path3_pair.wilcox.csv")
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/paired/PBC6m")
write.csv(res_species_pair.wilcox,"res_species_paired.wilcox.csv")

#2.clt1HC 
#2.1 species HC_new
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin")
maaslin_species_clt1HC=read.csv("maaslin_species clt1HC_new.csv",header = T,row.names = 1)
row.names(maaslin_species_clt1HC)=maaslin_species_clt1HC$species
maaslin_species_clt1HC_dif=filter(maaslin_species_clt1HC,qval<0.05) %>% row.names(.)
#2.2 genus HC_new
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin")
maaslin_genus_clt1HC=read.csv("maaslin_genus clt1HC_new.csv",header = T,row.names = 1)
row.names(maaslin_genus_clt1HC)=maaslin_genus_clt1HC$genus
maaslin_genus_clt1HC_dif=filter(maaslin_genus_clt1HC,qval<0.05) %>% row.names(.)
#2.3 KEGG_path3 HC_new
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin")
maaslin_KEGG_path3_clt1HC=read.csv("maaslin_KEGG_path3 clt1HC_new.csv",header = T,row.names = 1)
row.names(maaslin_KEGG_path3_clt1HC)=maaslin_KEGG_path3_clt1HC$Pathway
maaslin_KEGG_path3_clt1HC_dif=filter(maaslin_KEGG_path3_clt1HC,qval<0.05) %>% row.names(.)
#2.4 KO_gene HC_new n=362
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin")
maaslin_KO_gene_clt1HC=read.csv("maaslin_KO_gene clt1HC_new.csv",header = T,row.names = 1)
maaslin_KO_gene_clt1HC_dif=filter(maaslin_KO_gene_clt1HC,qval<0.05) %>% row.names(.)
#2.5 GMM HC_new
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin")
maaslin_GMM_clt1HC=read.csv("maaslin_GMM clt1HC_new.csv",header = T,row.names = 1)
row.names(maaslin_GMM_clt1HC)=maaslin_GMM_clt1HC$feature
maaslin_GMM_clt1HC_dif=filter(maaslin_GMM_clt1HC,qval<0.05) %>% row.names(.)

#3.clt1clt2
#3.1 species
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin")
maaslin_species_clt1clt2=read.csv("maaslin_species clt1clt2.csv",header = T,row.names = 1)
row.names(maaslin_species_clt1clt2)=maaslin_species_clt1clt2$species
maaslin_species_clt1clt2_dif=filter(maaslin_species_clt1clt2,qval<0.05) %>% row.names(.)
#2.2 genus
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin")
maaslin_genus_clt1clt2=read.csv("maaslin_genus clt1clt2.csv",header = T,row.names = 1)
row.names(maaslin_genus_clt1clt2)=maaslin_genus_clt1clt2$genus
maaslin_genus_clt1clt2_dif=filter(maaslin_genus_clt1clt2,qval<0.05) %>% row.names(.)
#2.3 KEGG_path3
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin")
maaslin_KEGG_path3_clt1clt2=read.csv("maaslin_KEGG_path3 clt1clt2.csv",header = T,row.names = 1)
row.names(maaslin_KEGG_path3_clt1clt2)=maaslin_KEGG_path3_clt1clt2$pathway
maaslin_KEGG_path3_clt1clt2_dif=filter(maaslin_KEGG_path3_clt1clt2,qval<0.05) %>% row.names(.)
#2.4 KO_gene
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin")
maaslin_KO_gene_clt1clt2=read.csv("maaslin_KO_gene clt1clt2.csv",header = T,row.names = 1)
row.names(maaslin_KO_gene_clt1clt2)=maaslin_KO_gene_clt1clt2$KO
maaslin_KO_gene_clt1clt2_dif=filter(maaslin_KO_gene_clt1clt2,qval<0.05) %>% row.names(.)
#2.5 GMM
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin")
maaslin_GMM_clt1clt2=read.csv("maaslin_GMM clt1clt2.csv",header = T,row.names = 1)
row.names(maaslin_GMM_clt1clt2)=maaslin_GMM_clt1clt2$feature
maaslin_GMM_clt1clt2_dif=filter(maaslin_GMM_clt1clt2,qval<0.05) %>% row.names(.)

#4.clt2HC
#4.1 species HC_new
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin")
maaslin_species_clt2HC=read.csv("maaslin_species clt2HC_new.csv",header = T,row.names = 1)
row.names(maaslin_species_clt2HC)=maaslin_species_clt2HC$species
maaslin_species_clt2HC_dif=filter(maaslin_species_clt2HC,qval<0.05) %>% row.names(.)
#4.2 genus HC_new
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin")
maaslin_genus_clt2HC=read.csv("maaslin_genus clt2HC_new.csv",header = T,row.names = 1)
row.names(maaslin_genus_clt2HC)=maaslin_genus_clt2HC$genus
maaslin_genus_clt2HC_dif=filter(maaslin_genus_clt2HC,qval<0.05) %>% row.names(.)

#4.3 KEGG_path3 HC_new
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin")
maaslin_KEGG_path3_clt2HC=read.csv("maaslin_KEGG_path3 clt2HC_new.csv",header = T,row.names = 1)
row.names(maaslin_KEGG_path3_clt2HC)=maaslin_KEGG_path3_clt2HC$Pathway
maaslin_KEGG_path3_clt2HC_dif=filter(maaslin_KEGG_path3_clt2HC,qval<0.05) %>% row.names(.)

#4.4 KO_gene HC_new n=1904
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin")
maaslin_KO_gene_clt2HC=read.csv("maaslin_KO_gene clt2HC_new.csv",header = T,row.names = 1)
maaslin_KO_gene_clt2HC_dif=filter(maaslin_KO_gene_clt2HC,qval<0.05) %>% row.names(.)
#4.5 GMM HC_new
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin")
maaslin_GMM_clt2HC=read.csv("maaslin_GMM clt2HC_new.csv",header = T,row.names = 1)
row.names(maaslin_GMM_clt2HC)=maaslin_GMM_clt2HC$feature
maaslin_GMM_clt2HC_dif=filter(maaslin_GMM_clt2HC,qval<0.05) %>% row.names(.)

################pair test feature selection for cluster##################
#确定clt1及clt2各自用于pair检验的feature:与HC相比有差异及各自对比有差异的feature
#1.clt1
#1.1 species HC_new
maaslin_species_clt1_pair=c(maaslin_species_clt1HC_dif,maaslin_species_clt1clt2_dif)
index=duplicated(maaslin_species_clt1_pair) 
index=which(index==FALSE)
maaslin_species_clt1_pair=maaslin_species_clt1_pair[index]

#1.2 genus
maaslin_genus_clt1_pair=c(maaslin_genus_clt1HC_dif,maaslin_genus_clt1clt2_dif)
index=duplicated(maaslin_genus_clt1_pair) 
index=which(index==FALSE)
maaslin_genus_clt1_pair=maaslin_genus_clt1_pair[index]

#1.3 KEGG_path3
maaslin_KEGG_path3_clt1_pair=c(maaslin_KEGG_path3_clt1HC_dif,maaslin_KEGG_path3_clt1clt2_dif)
index=duplicated(maaslin_KEGG_path3_clt1_pair) 
index=which(index==FALSE)
maaslin_KEGG_path3_clt1_pair=maaslin_KEGG_path3_clt1_pair[index]

#1.4 kO_gene
maaslin_KO_gene_clt1_pair=c(maaslin_KO_gene_clt1HC_dif,maaslin_KO_gene_clt1clt2_dif)
index=duplicated(maaslin_KO_gene_clt1_pair) 
index=which(index==FALSE)
maaslin_KO_gene_clt1_pair=maaslin_KO_gene_clt1_pair[index]

#1.5 GMM
maaslin_GMM_clt1_pair=c(maaslin_GMM_clt1HC_dif,maaslin_GMM_clt1clt2_dif)
index=duplicated(maaslin_GMM_clt1_pair) 
index=which(index==FALSE)
maaslin_GMM_clt1_pair=maaslin_GMM_clt1_pair[index]

#2.clt2
#2.1 species
maaslin_species_clt2_pair=c(maaslin_species_clt2HC_dif,maaslin_species_clt1clt2_dif)
index=duplicated(maaslin_species_clt2_pair) 
index=which(index==FALSE)
maaslin_species_clt2_pair=maaslin_species_clt2_pair[index]

#2.2 genus
maaslin_genus_clt2_pair=c(maaslin_genus_clt2HC_dif,maaslin_genus_clt1clt2_dif)
index=duplicated(maaslin_genus_clt2_pair) 
index=which(index==FALSE)
maaslin_genus_clt2_pair=maaslin_genus_clt2_pair[index]

#2.3 KEGG_path3
maaslin_KEGG_path3_clt2_pair=c(maaslin_KEGG_path3_clt2HC_dif,maaslin_KEGG_path3_clt1clt2_dif)
index=duplicated(maaslin_KEGG_path3_clt2_pair) 
index=which(index==FALSE)
maaslin_KEGG_path3_clt2_pair=maaslin_KEGG_path3_clt2_pair[index]

#1.4 kO_gene
maaslin_KO_gene_clt2_pair=c(maaslin_KO_gene_clt2HC_dif,maaslin_KO_gene_clt1clt2_dif)
index=duplicated(maaslin_KO_gene_clt2_pair) 
index=which(index==FALSE)
maaslin_KO_gene_clt2_pair=maaslin_KO_gene_clt2_pair[index]

#1.5 GMM
maaslin_GMM_clt2_pair=c(maaslin_GMM_clt2HC_dif,maaslin_GMM_clt1clt2_dif)
index=duplicated(maaslin_GMM_clt2_pair) 
index=which(index==FALSE)
maaslin_GMM_clt2_pair=maaslin_GMM_clt2_pair[index]

############可视化功能变化##############
#1.KEGG_path3
maaslin_KEGG_path3_PBCHC_plot=maaslin_KEGG_path3_PBCHC[maaslin_KEGG_path3_PBCHC_dif,]
library(plyr)
library(forcats)
maaslin_KEGG_path3_PBCHC_plot[which(maaslin_KEGG_path3_PBCHC_plot$coef<0),"change"]="HC"
maaslin_KEGG_path3_PBCHC_plot[which(maaslin_KEGG_path3_PBCHC_plot$coef>0),"change"]="PBCbase"
maaslin_KEGG_path3_PBCHC_plot$Pathway=as.character(maaslin_KEGG_path3_PBCHC_plot$Pathway)
maaslin_KEGG_path3_PBCHC_plot$change=as.factor(maaslin_KEGG_path3_PBCHC_plot$change)
maaslin_KEGG_path3_PBCHC_plot <- maaslin_KEGG_path3_PBCHC_plot %>%  mutate(Pathway = fct_reorder(Pathway,coef))

library(ggplot2)
ggplot(data = maaslin_KEGG_path3_PBCHC_plot,aes(x=Pathway,y=coef,fill=change))+
  geom_bar(stat = "identity",position = "identity",width = 0.2)+
  scale_fill_manual(values = c("#3F5BA7","#D72326"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  theme(legend.position ="none")+
  theme(text = element_text(size = 7))+
  labs(fill="")+
  coord_flip()

#2.Tryptophan metabolism by PBCHC group HC_new
trypto_gene=KO_gene_abs_clr[c(PBCbase,HC_new),c("K01501","K07130","K00252","K00463","K07515","K00626",
                                            "K00274","K00838","K01556","K00543","K01426","K11182","K00825",
                                            "K01593","K01721","K07511","K03781","K01692","K01825","K00164","K00022",
                                            "K00466","K00816","K04103","K00128","K01432","K03782","K00392","K01667","K01782",
                                            "K15067","K00169","K00170","K00171","K00172"
)]
trypto_gene=cbind(trypto_gene,metadata_PBCHC[c(PBCbase,HC_new),"group"])  
colnames(trypto_gene)[36]="group"

trypto_gene_mean=as.data.frame(matrix(0,ncol=1,nrow=2))
str(trypto_gene)
for (i in 1:35) {
  res.gene=aggregate(trypto_gene[,i], by=list(type=trypto_gene$group),mean)
  trypto_gene_mean[,i+1]=res.gene$x
  row.names(trypto_gene_mean)=res.gene$type
}
trypto_gene_mean=trypto_gene_mean[,-1]
colnames(trypto_gene_mean)=c("E3.5.5.1","kynB","gcdH","IDO","HADHA","atoB","MAO","ARO8","kynU",
                             "ASMT","amiE","AOC1","AADAT","DDC","nthA","ECHS1","katE","paaF","fadB","OGDH","HADH",
                             "E1.13.12.3","CCBL","ipdC","ALDH","AFMID","katG","ACMSD",
                             "tnaA","fadJ","amnD","porA","porB","porD","porG")
library(pheatmap)
pheatmap(t(trypto_gene_mean),cluster_rows = T,cluster_cols = F,
         color = colorRampPalette(c("#12528C", "white", "#B7232E"))(50),
         main = "Tryptophan metabolism",fontsize = 7)
#2.1 microbial Tryptophan metabolism by PBCHC group HC_new  #还未完成
trypto_gene_micro=KO_gene_abs_clr[c(PBCbase,HC_new),c("K01501","K07130","K00252","K00463","K07515","K00626",
                                                "K00274","K00838","K01556","K00543","K01426","K11182","K00825",
                                                "K01593","K01721","K07511","K03781","K01692","K01825","K00164","K00022",
                                                "K00466","K00816","K04103","K00128","K01432","K03782","K00392","K01667","K01782",
                                                "K15067","K00169","K00170","K00171","K00172"
)]
trypto_gene_micro=cbind(trypto_gene_micro,metadata_PBCHC[c(PBCbase,HC_new),"group"])  
colnames(trypto_gene)[36]="group"

trypto_gene_mean=as.data.frame(matrix(0,ncol=1,nrow=2))
str(trypto_gene)
for (i in 1:35) {
  res.gene=aggregate(trypto_gene[,i], by=list(type=trypto_gene$group),mean)
  trypto_gene_mean[,i+1]=res.gene$x
  row.names(trypto_gene_mean)=res.gene$type
}
trypto_gene_mean=trypto_gene_mean[,-1]
colnames(trypto_gene_mean)=c("E3.5.5.1","kynB","gcdH","IDO","HADHA","atoB","MAO","ARO8","kynU",
                             "ASMT","amiE","AOC1","AADAT","DDC","nthA","ECHS1","katE","paaF","fadB","OGDH","HADH",
                             "E1.13.12.3","CCBL","ipdC","ALDH","AFMID","katG","ACMSD",
                             "tnaA","fadJ","amnD","porA","porB","porD","porG")
library(pheatmap)
pheatmap(t(trypto_gene_mean),cluster_rows = T,cluster_cols = F,
         color = colorRampPalette(c("#12528C", "white", "#B7232E"))(50),
         main = "Tryptophan metabolism",fontsize = 7)

###########零碎分析###############
library(dplyr)
#比较shannon指数
wilcox.test(shannon[PBCbase_resp_6m,1],shannon[PBCbase_noresp_6m,1])
#responder and noresponder
#PBCbase and HC_new
lm(shannon~Response_6m+age+gender+BMI,data = metadata_PBCbase) %>% summary(.)#校正后不显著
lm(shannon~group+age+gender+BMI,data = metadata_PBCHC[c(PBCbase,HC_new),]) %>% summary(.)#校正后P=0.00723

boxplot(shannon[PBCbase_resp_6m,1],shannon[PBCbase_noresp_6m,1],names = c("responder","no-responder"),ylab="shannon"
        ,col = c("#c9daf8","#e69138"),main = " wilcox P=0.11")
#simpson指数
simpson=diversity(t(species_rela_filter),"simpson") %>% as.data.frame()
wilcox.test(simpson[PBCbase_resp_6m,1],simpson[PBCbase_noresp_6m,1])#P=0.20
boxplot(simpson[PBCbase_resp_6m,1],simpson[PBCbase_noresp_6m,1],names = c("responder","no-responder"),ylab="simpson"
        ,col = c("#c9daf8","#e69138"),main = " wilcox P=0.20")

boxplot(Cazy2_abs_clr[HC,"GH128"],Cazy2_abs_clr[PBCbase_cluster2,"GH128"])

row.names(Cazy1_abs)
#比较cluster储存时间的差异
wilcox.test(metadata[PBCbase_cluster1,"Storage_time"],metadata[PBCbase_cluster2,"Storage_time"])


###########可视化###########
#cazy1
CAzy1_abs_clr_plot=cbind(CAzy1_abs_clr[row.names(metadata_PBCHC),],metadata_PBCHC[,"cluster"])
colnames(CAzy1_abs_clr_plot)[7]="cluster"
CAzy1_abs_clr_plot$cluster=factor(CAzy1_abs_clr_plot$cluster,levels = c("clt1","clt2","HC"))
library(reshape2)
library(reshape)
library(ggplot2)
str(melt(CAzy1_abs_clr_plot))
ggplot(data=melt(CAzy1_abs_clr_plot[,c(4,6,7)]),aes(x=variable,y=value,color=cluster))+geom_boxplot()+
  scale_color_manual(values = c("#FFB733","#B75454","#83AAF0"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank())+
  theme(axis.text.x = element_text(angle = -10,hjust = 0.4,vjust =1))+ #调整x轴文本的方向
  ylab("CLR abundance")+xlab("")+
  annotate("text",x=2,y=1,label="clt1 vs HC/clt2 FDR<0.05 \n clt2vsHC无差异")
#cazy1在PBCbase cluster1治疗前后变化
wilcox.test(CAzy1_abs_clr[PBCpair_6m_cluster1,4],CAzy1_abs_clr[PBC_6m_cluster1,4],paired = T)
wilcox.test(CAzy1_abs_clr[PBCpair_6m_cluster1,6],CAzy1_abs_clr[PBC_6m_cluster1,6],paired = T)
boxplot(CAzy1_abs_clr[PBCpair_6m_cluster1,4],CAzy1_abs_clr[PBC_6m_cluster1,4])
data=data.frame(baseline=CAzy1_abs_clr[PBCpair_6m_cluster1,6],treat=CAzy1_abs_clr[PBC_6m_cluster1,6]) #shannon.base|6m为向量
colnames(CAzy1_abs_clr)
library(ggpubr)
ggpaired(data, cond1 = "baseline", cond2 = "treat",
         color = "condition", palette = "npg",ylab = "Glycoside Hydrolases(CLR)",line.color = "grey",
         title = "cluster1 P=0.017")+
  theme(legend.position = "none") 

KO_gene_abs["K13368",PBCbase_cluster2]
KO_gene_abs["K13368",PBCbase_cluster1]
KO_gene_abs["K13368",HC]

#可视化展示通路的变化
#1.clt1.HC
KEGG_path3_maalin_plot_clt1HC=read.csv("KEGG_path3_clt1HC_for plot.csv",header = T)
row.names(KEGG_path3_maalin_plot_clt1HC)=KEGG_path3_maalin_plot_clt1HC$KEGG_pathway
colnames(KEGG_path3_maalin_plot_clt1HC)
library(plyr)
library(forcats)
KEGG_path3_maalin_plot_clt1HC <- KEGG_path3_maalin_plot_clt1HC %>%  mutate(KEGG_pathway = fct_reorder(KEGG_pathway,coef))

#top20 pathway in clt1 vs HC
ggplot(data = KEGG_path3_maalin_plot_clt1HC,aes(x=KEGG_pathway,y=coef,fill=change))+
  geom_bar(stat = "identity",position = "identity")+
  scale_fill_manual(values = c("#83AAF0","#F2644B"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  theme(legend.position =c(0.3,0.8))+
  labs(fill="")+
  coord_flip()
#2. clt2.HC
#1.clt1.HC
KEGG_path3_maalin_plot_clt2HC=read.csv("KEGG_path3_clt2HC_for plot.csv",header = T)
row.names(KEGG_path3_maalin_plot_clt2HC)=KEGG_path3_maalin_plot_clt2HC$KEGG_pathway
colnames(KEGG_path3_maalin_plot_clt2HC)
library(plyr)
library(forcats)
KEGG_path3_maalin_plot_clt2HC <- KEGG_path3_maalin_plot_clt2HC %>%  mutate(KEGG_pathway = fct_reorder(KEGG_pathway,coef))

#top20 pathway in clt2 vs HC
ggplot(data = KEGG_path3_maalin_plot_clt2HC,aes(x=KEGG_pathway,y=coef,fill=change))+
  geom_bar(stat = "identity",position = "identity")+
  scale_fill_manual(values = c("#F2644B","#83AAF0"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  theme(legend.position =c(0.8,0.3))+
  labs(fill="")+
  coord_flip()

##############代谢通路及相关基因的变化##############
KEGG_path3_abs_clr_plot=cbind(KEGG_path3_abs_clr[row.names(metadata_PBCHC),],metadata_PBCHC[,"cluster"])
colnames(KEGG_path3_abs_clr_plot)[416]="cluster"
KEGG_path3_abs_clr_plot$cluster=factor(KEGG_path3_abs_clr_plot$cluster,levels = c("clt1","clt2","HC"))
#1.steroid biosynthesis pathway
ggplot(data=KEGG_path3_abs_clr_plot,aes(y=KEGG_path3_abs_clr_plot$ko00100__Steroid_biosynthesis,x=cluster,color=cluster))+geom_boxplot()+
  scale_color_manual(values = c("#FFB733","#B75454","#83AAF0"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  ylab("Steroid_biosynthesis (CLR)")+xlab("")+
  annotate("text",x=2,y=2.5,label="MaAslin clt1 vs HC FDR=0.02\n clt2vsHC FDR=2.19E-07\nclt1vs.clt2 FDR=2.23E-11")
#1.1 steroid biosynthesis通路中的基因变化
steroid_gene=KO_gene_abs_clr[c(PBCbase,HC),c("K00559","K08242","K07750","K12298","K00801","K07748")]
steroid_gene=cbind(steroid_gene,metadata_PBCHC[c(PBCbase,HC),"cluster"])
colnames(steroid_gene)[7]="cluster"
library(dplyr)
steroid_gene_mean=as.data.frame(matrix(0,ncol=1,nrow=3))
for (i in 1:6) {
  res.gene=aggregate(steroid_gene[,i], by=list(type=steroid_gene$cluster),mean)
  steroid_gene_mean[,i+1]=res.gene$x
  row.names(steroid_gene_mean)=res.gene$type
}
steroid_gene_mean=steroid_gene_mean[,-1]
colnames(steroid_gene_mean)=c("SMT1","SMT2","SC4MOL","CEL","FDFT1","NSDHL")

pheatmap(t(steroid_gene_mean),cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Steroid biosynthesis")

#2.carotenoid pathway
ggplot(data=KEGG_path3_abs_clr_plot,aes(y=KEGG_path3_abs_clr_plot$ko00906__Carotenoid_biosynthesis,x=cluster,color=cluster))+geom_boxplot()+
  scale_color_manual(values = c("#FFB733","#B75454","#83AAF0"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  ylab("Carotenoid_biosynthesis (CLR)")+xlab("")+
  annotate("text",x=2,y=3.5,label="MaAslin clt1 vs HC ns\n clt2vsHC FDR=5.83E-07\nclt1vs.clt2 FDR=8.04E-06")
#2.1 carotenoid gene 
which(KO_gene_abs["K05296",]!=0) %>% length(.)
carotenoid_gene=KO_gene_abs_clr[c(PBCbase,HC),c("K10212","K02292","K14596","K02293","K09841","K09844",
                                                "K09835","K10210","K09846","K02291","K10211","K15745","K10027")]

carotenoid_gene=cbind(carotenoid_gene,metadata_PBCHC[c(PBCbase,HC),"cluster"])
colnames(carotenoid_gene)[14]="cluster"
library(dplyr)
carotenoid_gene_mean=as.data.frame(matrix(0,ncol=1,nrow=3))
for (i in 1:13) {
  res.gene=aggregate(carotenoid_gene[,i], by=list(type=carotenoid_gene$cluster),mean)
  carotenoid_gene_mean[,i+1]=res.gene$x
  row.names(carotenoid_gene_mean)=res.gene$type
}
carotenoid_gene_mean=carotenoid_gene_mean[,-1]
colnames(carotenoid_gene_mean)=c("K10212","crtO","crtX","crtP","ABA2","crtC",
                                 "crtISO","crtP","crtF","crtB","crtQ","AL1","crtI")

pheatmap(t(carotenoid_gene_mean),cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "carotenoid biosynthesis")
#3. LPS biosynthesis
ggplot(data=KEGG_path3_abs_clr_plot,aes(y=KEGG_path3_abs_clr_plot$ko00540__Lipopolysaccharide_biosynthesis,x=cluster,color=cluster))+geom_boxplot()+
  scale_color_manual(values = c("#FFB733","#B75454","#83AAF0"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  ylab("LPS biosynthesis (CLR)")+xlab("")
#3.1 LPS gene (only some representative,太多了)

LPS_gene=KO_gene_abs_clr[c(PBCbase,HC),c("K06041","K03274","K03269","K02535","K19353",
                                                "K00912","K01627","K02527","K00748","K02841","K02560","K02517")]

LPS_gene=cbind(LPS_gene,metadata_PBCHC[c(PBCbase,HC),"cluster"])
colnames(LPS_gene)[13]="cluster"
library(dplyr)
LPS_gene_mean=as.data.frame(matrix(0,ncol=1,nrow=3))
for (i in 1:12) {
  res.gene=aggregate(LPS_gene[,i], by=list(type=LPS_gene$cluster),mean)
  LPS_gene_mean[,i+1]=res.gene$x
  row.names(LPS_gene_mean)=res.gene$type
}
LPS_gene_mean=LPS_gene_mean[,-1]
colnames(LPS_gene_mean)=c("kdsD","gmhD","lpxH","lpxC",
                                 "eptC","lpxK","kdsA","kdtA","lpxB","waaC","lpxM","lpxL")

pheatmap(t(LPS_gene_mean),cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "LPS biosynthesis")
#4 indole propionic acid  
IPA_gene=KO_gene_abs_clr[c(PBCbase,HC),c("K00169","K00170","K00171","K00172","K13607"
                                         )]
IPA_gene=cbind(IPA_gene,metadata_PBCHC[c(PBCbase,HC),"cluster"])
colnames(IPA_gene)[6]="cluster"
library(dplyr)
IPA_gene_mean=as.data.frame(matrix(0,ncol=1,nrow=3))
str(IPA_gene)
for (i in 1:5) {
  res.gene=aggregate(IPA_gene[,i], by=list(type=IPA_gene$cluster),mean)
  IPA_gene_mean[,i+1]=res.gene$x
  row.names(IPA_gene_mean)=res.gene$type
}
IPA_gene_mean=IPA_gene_mean[,-1]
colnames(IPA_gene_mean)=c("porA","porB","porD","porG","fldA")
pheatmap(t(IPA_gene_mean),cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Indole biosynthesis")


#5.糖苷酶基因变化
glycan_gene=KO_gene_abs_clr[c(PBCbase,HC),c("K01818","K01784","K02564","K01209","K12373","K01206","K15923",
                                            "K05989","K01179","K01178","K01190","K07407","K05349","K01187",
                                            "K21065","K01186","K01811","K01235")]
glycan_gene=cbind(glycan_gene,metadata_PBCHC[c(PBCbase,HC),"cluster"])  
colnames(glycan_gene)[19]="cluster"
str(trypto_gene)
glycan_gene_mean=as.data.frame(matrix(0,ncol=1,nrow=3))
str(glycan_gene)
for (i in 1:18) {
  res.gene=aggregate(glycan_gene[,i], by=list(type=glycan_gene$cluster),mean)
  glycan_gene_mean[,i+1]=res.gene$x
  row.names(glycan_gene_mean)=res.gene$type
}
glycan_gene_mean=glycan_gene_mean[,-1]
colnames(glycan_gene_mean)=c("fucI","galE","nagB","abfA","HEXA_B","FUCA","AXY8","ramA","E3.2.1.4",
                             "E3.2.1.3","lacZ","rafA","bglX","malZ","E3.2.1.197","NEU1","xylS","aguA")

pheatmap(t(glycan_gene_mean),cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "polysaccride metabolism")

#6. steroid hormone biosynthesis
steroid.hormone_gene=KO_gene_abs_clr[c(PBCbase,HC),c("K01130","K00699","K07414","K07424",
                                                     "K13368","K00038","K00071","K00545",
                                                    "K15680","K01131","K12343","K00251","K13370","K05296")]

steroid.hormone_gene=cbind(steroid.hormone_gene,metadata_PBCHC[c(PBCbase,HC),"cluster"])  
colnames(steroid.hormone_gene)[15]="cluster"
str(steroid.hormone_gene$cluster)
steroid.hormone_gene_mean=as.data.frame(matrix(0,ncol=1,nrow=3))
str(steroid.hormone_gene)
for (i in 1:14) {
  res.gene=aggregate(steroid.hormone_gene[,i], by=list(type=steroid.hormone_gene$cluster),mean)
  steroid.hormone_gene_mean[,i+1]=res.gene$x
  row.names(steroid.hormone_gene_mean)=res.gene$type
}
steroid.hormone_gene_mean=steroid.hormone_gene_mean[,-1]
colnames(steroid.hormone_gene_mean)=c("aslA","UGT","CYP2D","CYP3A","HSD17B2","HSD3a","HSD11B2","COMT","HSD11B1","STS",
                                      "SRD5A1","AKR1D1","HSD17B8","HSD3B"
                          )

library(pheatmap)
pheatmap(t(steroid.hormone_gene_mean),cluster_rows = T,cluster_cols = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "steroid hormone biosynthesis")

#6.1 steroid hormone biosynthesis 根据group显示 HC_new
steroid.hormone_gene=steroid.hormone_gene[c(PBCbase,HC_new),]
steroid.hormone_gene=cbind(steroid.hormone_gene,metadata_PBCHC[c(PBCbase,HC_new),"group"])  
colnames(steroid.hormone_gene)[16]="group"
str(steroid.hormone_gene$group)
steroid.hormone_gene_mean_group=as.data.frame(matrix(0,ncol=1,nrow=2))
str(steroid.hormone_gene)
for (i in 1:14) {
  res.gene=aggregate(steroid.hormone_gene[,i], by=list(type=steroid.hormone_gene$group),mean)
  steroid.hormone_gene_mean_group[,i+1]=res.gene$x
  row.names(steroid.hormone_gene_mean_group)=res.gene$type
}
steroid.hormone_gene_mean_group=steroid.hormone_gene_mean_group[,-1]
colnames(steroid.hormone_gene_mean_group)=c("aslA","UGT","CYP2D","CYP3A","HSD17B2","HSD3a","HSD11B2","COMT","HSD11B1","STS",
                                      "SRD5A1","AKR1D1","HSD17B8","HSD3B"
)

library(pheatmap)
pheatmap(t(steroid.hormone_gene_mean_group),cluster_rows = T,cluster_cols = T,scale = "none",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "steroid hormone biosynthesis")


#6.2 steroid hormone biosynthesis_selected
steroid.hormone.select_gene=KO_gene_abs_clr[c(PBCbase,HC),c("K01131","K01130","K12343")]

steroid.hormone.select_gene=cbind(steroid.hormone.select_gene,metadata_PBCHC[c(PBCbase,HC),"cluster"])  
colnames(steroid.hormone.select_gene)[4]="cluster"
str(steroid.hormone.select_gene$cluster)
steroid.hormone.select_gene_mean=as.data.frame(matrix(0,ncol=1,nrow=3))

for (i in 1:3) {
  res.gene=aggregate(steroid.hormone.select_gene[,i], by=list(type=steroid.hormone.select_gene$cluster),mean)
  steroid.hormone.select_gene_mean[,i+1]=res.gene$x
  row.names(steroid.hormone.select_gene_mean)=res.gene$type
}
steroid.hormone.select_gene_mean=steroid.hormone.select_gene_mean[,-1]
colnames(steroid.hormone.select_gene_mean)=c("STS","aslA","SRD5A1")

library(pheatmap)
pheatmap(t(steroid.hormone.select_gene_mean),cluster_rows = T,cluster_cols = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "steroid hormone biosynthesis")


#############物种堆叠柱状图#############
colSums(phylum_rela)
row.names(phylum_filter)
library(microbiomeseq)

phylum_filter_stack.plot=phylum_rela
phylum_filter_stack.plot$mean=apply(phylum_filter_stack.plot,1,mean)
phylum_filter_stack.plot$name=row.names(phylum_rela)
phylum_filter_stack.plot=arrange(phylum_filter_stack.plot,desc(phylum_filter_stack.plot$mean))
row.names(phylum_filter_stack.plot)=phylum_filter_stack.plot$name
phylum_filter_stack.plot.top10=phylum_filter_stack.plot[c(1,2,4,6,8,9,11,12,13,14),-c(412,413)]
phylum_filter_stack.plot.other=phylum_filter_stack.plot[-c(1,2,4,6,8,9,11,12,13,14),-c(412,413)]
phylum_other=apply(phylum_filter_stack.plot.other, 2, sum)
phylum_filter_stack.plot.top10[11,]=phylum_other
colSums(phylum_filter_stack.plot.top10)
phylum_filter_stack.plot.top10_clt1.base=phylum_filter_stack.plot.top10[,PBCbase_cluster1]
#1.cluster1
phylum_filter_stack.plot.top10_clt1.base$Taxonomy <- factor(rownames(phylum_filter_stack.plot.top10_clt1.base), levels = rev(rownames(phylum_filter_stack.plot.top10_clt1.base)))
phylum_filter_stack.plot.top10_clt1.base <- melt(phylum_filter_stack.plot.top10_clt1.base, id = 'Taxonomy')
str(phylum_filter_stack.plot.top10_clt1.base)
p <- ggplot(phylum_filter_stack.plot.top10_clt1.base, aes(x=reorder(variable,value), 100 * value, fill = Taxonomy)) +
  geom_col(position = 'stack', width = 0.6) +

  scale_fill_manual(values =  rev(c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', 'gray'))) +
  labs(x = '', y = 'Relative Abundance(%)') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13), legend.title = element_blank(), legend.text = element_text(size = 11))
#2.cluster2
phylum_filter_stack.plot.top10_clt2.base=phylum_filter_stack.plot.top10[,PBCbase_cluster2]

phylum_filter_stack.plot.top10_clt2.base$Taxonomy <- factor(rownames(phylum_filter_stack.plot.top10_clt2.base), levels = rev(rownames(phylum_filter_stack.plot.top10_clt2.base)))
phylum_filter_stack.plot.top10_clt2.base <- melt(phylum_filter_stack.plot.top10_clt2.base, id = 'Taxonomy')
str(phylum_filter_stack.plot.top10_clt2.base)
ggplot(phylum_filter_stack.plot.top10_clt2.base, aes(x=reorder(variable,value), 100 * value, fill = Taxonomy)) +
  geom_col(position = 'stack', width = 0.6) +
  
  scale_fill_manual(values =  rev(c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', 'gray'))) +
  labs(x = '', y = 'Relative Abundance(%)') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13), legend.title = element_blank(), legend.text = element_text(size = 11))
#3.HC
phylum_filter_stack.plot.top10_HC=phylum_filter_stack.plot.top10[,HC]

phylum_filter_stack.plot.top10_HC$Taxonomy <- factor(rownames(phylum_filter_stack.plot.top10_HC), levels = rev(rownames(phylum_filter_stack.plot.top10_HC)))
phylum_filter_stack.plot.top10_HC <- melt(phylum_filter_stack.plot.top10_HC, id = 'Taxonomy')
str(phylum_filter_stack.plot.top10_HC)
ggplot(phylum_filter_stack.plot.top10_HC, aes(x=reorder(variable,value), 100 * value, fill = Taxonomy)) +
  geom_col(position = 'stack', width = 0.6) +
  scale_fill_manual(values =  rev(c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', 'gray'))) +
  labs(x = '', y = 'Relative Abundance(%)') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13), legend.title = element_blank(), legend.text = element_text(size = 11))+
  title("HC")

############计算基因的presence样本数############
gene_presence=vector()
for (i in 1:11276) {
  gene_presence[i]=which(KO_gene_abs[i,]!=0) %>% length(.)
}
gene_presence=as.data.frame(gene_presence)
row.names(gene_presence)=row.names(KO_gene_abs)
setwd("D:/PBC multiomics data/unsupervised cluster analysis_renamed/metagenomics/result/KO gene")
write.csv(gene_presence,"all_KO_gene presence.csv")
#############计算species的PBCbase presence样本数###############
species_filter_presence=vector()
for (i in 1:3854) {
  species_filter_presence[i]=which(species_abs[row.names(species_rela_filter),PBCbase][i,]!=0) %>% length(.)
}
species_filter_presence=as.data.frame(species_filter_presence)
row.names(species_filter_presence)=row.names(species_rela_filter)
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species")
write.csv(species_filter_presence,"species_filter_presence_PBCbase.csv")
##############特定基因/taxa/GMM的boxplot *figure################
#特定基因 HC_new
KO_gene_abs_clr_plot=cbind(KO_gene_abs_clr[c(HC_new,PBCbase),row.names(KO_gene_rela_filter)],metadata_PBCHC.new[c(HC_new,PBCbase),"group"])
colnames(KO_gene_abs_clr_plot)[9537]="group"
KO_gene_abs_clr_plot=cbind(KO_gene_abs_clr_plot,metadata_PBCHC.new[c(HC_new,PBCbase),"cluster"])
colnames(KO_gene_abs_clr_plot)[9538]="cluster"
KO_gene_abs_clr_plot$cluster=factor(KO_gene_abs_clr_plot$cluster,levels = c("HC","clt1","clt2"))
str(KO_gene_abs_clr_plot$cluster)
library(ggplot2)
#1.按照group展示 **figure
ggplot(data = KO_gene_abs_clr_plot,aes(x=group,y=KO_gene_abs_clr_plot$K01667,fill=group))+geom_boxplot()+
  labs(y="tnaA (CLR)")+
  scale_fill_manual(values = c("#3F5BA7","#D72326"))+
  theme(legend.position = "none")
#2.按照cluster展示
ggplot(data = KO_gene_abs_clr_plot,aes(x=cluster,y=KO_gene_abs_clr_plot$K11181,fill=cluster))+geom_boxplot()+
  labs(y="dsrB (CLR)")+
  scale_fill_manual(values = c("#3F5BA7","#088288","#F19433"))+
  theme(legend.position = "none")

#特定species
species_abs_clr_plot=cbind(species_abs_clr[c(PBCbase,HC_new),row.names(species_rela_filter)],metadata_PBCHC.new[c(PBCbase,HC_new),"group"])
colnames(species_abs_clr_plot)[3855]="group"
species_abs_clr_plot=cbind(species_abs_clr_plot,metadata_PBCHC.new[c(PBCbase,HC_new),"cluster"])
colnames(species_abs_clr_plot)[3856]="cluster"
species_abs_clr_plot$cluster=factor(species_abs_clr_plot$cluster,levels = c("HC","clt1","clt2"))
str(species_abs_clr_plot$cluster)
ggplot(data = species_abs_clr_plot,aes(x=cluster,y=species_abs_clr_plot$"k__Bacteria;p__Proteobacteria;c__Deltaproteobacteria;o__Desulfovibrionales;f__Desulfovibrionaceae;g__Bilophila;s__Bilophila wadsworthia",fill=cluster))+geom_boxplot(color="grey")+
  labs(y="Bilophila wadsworthia (CLR)")+
  scale_fill_manual(values = c("#3F5BA7","#088288","#F19433"))+
  theme(legend.position = "none")

#特定GMM
GMM_abs_clr_plot=cbind(GMM_abs_clr[c(HC_new,PBCbase),row.names(GMM_rela_filter)],metadata_PBCHC.new[c(HC_new,PBCbase),"group"])
colnames(GMM_abs_clr_plot)[130]="group"
GMM_abs_clr_plot=cbind(GMM_abs_clr_plot,metadata_PBCHC.new[c(HC_new,PBCbase),"cluster"])
colnames(GMM_abs_clr_plot)[131]="cluster"
str(GMM_abs_clr_plot$cluster)
#1.display according to group
ggplot(data = GMM_abs_clr_plot,aes(x=group,y=GMM_abs_clr_plot$MF0099,fill=group))+geom_boxplot(color="grey")+
  labs(y="methanol conversion (CLR)")+
  scale_fill_manual(values = c("#3F5BA7","#D72326"))+
  theme(axis.line = element_line(colour = "grey"),panel.background = element_blank(),legend.position = "none")
  
#特定genus
genus_abs_clr_plot=cbind(genus_abs_clr[c(HC_new,PBCbase),row.names(genus_rela_filter)],metadata_PBCHC.new[c(HC_new,PBCbase),"group"])
colnames(genus_abs_clr_plot)[884]="group"
genus_abs_clr_plot=cbind(genus_abs_clr_plot,metadata_PBCHC.new[c(HC_new,PBCbase),"cluster"])
colnames(genus_abs_clr_plot)[885]="cluster"
str(genus_abs_clr_plot$cluster)

#1.display according to group
ggplot(data = genus_abs_clr_plot,aes(x=cluster,y=genus_abs_clr_plot$`k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides`,fill=cluster))+geom_boxplot(color="grey")+
  labs(y="Bacteroides (CLR)")+
  scale_fill_manual(values =  c("#3F5BA7","#088288","#F19433"))+
  theme(axis.line = element_line(colour = "grey"),panel.background = element_blank(),legend.position = "none")


#特定order
order_abs_clr_plot=cbind(order_abs_clr[c(HC_new,PBCbase),row.names(order_rela_filter)],metadata_PBCHC.new[c(HC_new,PBCbase),"group"])
colnames(order_abs_clr_plot)[204]="group"
order_abs_clr_plot=cbind(order_abs_clr_plot,metadata_PBCHC.new[c(HC_new,PBCbase),"cluster"])
colnames(order_abs_clr_plot)[205]="cluster"
str(order_abs_clr_plot$cluster)
#1.display according to group
ggplot(data = order_abs_clr_plot,aes(x=cluster,y=order_abs_clr_plot$`k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales`,fill=cluster))+geom_boxplot(color="grey")+
  labs(y="Clostridiales (CLR)")+
  scale_fill_manual(values =  c("#3F5BA7","#088288","#F19433"))+
  theme(axis.line = element_line(colour = "grey"),panel.background = element_blank(),legend.position = "none")

#特定class
class_abs_clr_plot=cbind(class_abs_clr[c(HC_new,PBCbase),row.names(class_rela_filter)],metadata_PBCHC.new[c(HC_new,PBCbase),"group"])
colnames(class_abs_clr_plot)[119]="group"
class_abs_clr_plot=cbind(class_abs_clr_plot,metadata_PBCHC.new[c(HC_new,PBCbase),"cluster"])
colnames(class_abs_clr_plot)[120]="cluster"
str(class_abs_clr_plot$cluster)
#1.display according to group
ggplot(data = class_abs_clr_plot,aes(x=cluster,y=class_abs_clr_plot$`k__Bacteria;p__Firmicutes;c__Clostridia`,fill=cluster))+geom_boxplot(color="grey")+
  labs(y="Clostridia (CLR)")+
  scale_fill_manual(values =  c("#3F5BA7","#088288","#F19433"))+
  theme(axis.line = element_line(colour = "grey"),panel.background = element_blank(),legend.position = "none")


####################绘制microbial feature变化热图###########
#species===================
#计算species平均相对丰度
species_rela_mean=apply(species_rela_filter, 1, mean) %>% as.data.frame(.)
colnames(species_rela_mean)="mean"

##根据mean relative abundance进行排序
species_rela_mean_dif=species_rela_mean[maaslin_species_PBCHC_dif,] %>% as.data.frame(.)
row.names(species_rela_mean_dif)=maaslin_species_PBCHC_dif
colnames(species_rela_mean_dif)="mean"
species_rela_mean_dif$name=maaslin_species_PBCHC_dif
species_rela_mean_dif=arrange(species_rela_mean_dif,desc(mean))
species_rela_mean_dif[1:20,"name"]

#1.P值排名前30位的species PBC vs HC *figure
species_Top30=maaslin_species_PBCHC_dif[1:30]
species_abs_clr_top30=species_abs_clr[c(HC_new,PBCbase),species_Top30]
newname=vector()
name=strsplit(colnames(species_abs_clr_top30),"s__")
for (i in 1:30) {
  newname[i]=name[[i]][2]
}
newname[8]="Veillonella spp"
newname[15]="Anaerobutyricum spp"
colnames(species_abs_clr_top30)=newname

#1.1按照全样本热图展示
anno_sample=metadata[row.names(species_abs_clr_top30),"group"] %>% as.data.frame(.)
row.names(anno_sample)=row.names(species_abs_clr_top30)
colnames(anno_sample)="group"
anno_sample$group=anno_sample$group[,drop=TRUE]
str(anno_sample)
color3 = colorRampPalette(c("white","#2F5597"))(50)
ann_colors = list(
  group = c(PBCbase="#D72326",HC="#3F5BA7")
)
library(pheatmap)
pheatmap(t(species_abs_clr_top30),  show_colnames = FALSE,
         fontsize = 8,cluster_cols = F,color = color3,annotation_col = anno_sample,annotation_colors = ann_colors,
         border=TRUE,gaps_col = 131,scale = "row",breaks =  seq(0,4,length.out = 110))

pdf(3*7)
#输出文件
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin")
write.csv(species_rela_mean,"species_rela_mean.csv")

#1.2 将PBCHC_new 变化最显著的top30按照barplot展示
maaslin_species_PBCHC_top30=maaslin_species_PBCHC[1:30,]
maaslin_species_PBCHC_top30$species=newname
maaslin_species_PBCHC_top30=arrange(maaslin_species_PBCHC_top30,desc(coef))

library(plyr)
library(forcats)
library(dplyr)
#避免barplot按照字符排名
maaslin_species_PBCHC_top30 <- maaslin_species_PBCHC_top30 %>%  mutate(species = fct_reorder(species,coef))
maaslin_species_PBCHC_top30$group=c(rep("PBCbase",14),rep("HC",16))

ggplot(data = maaslin_species_PBCHC_top30,aes(x=species,y=coef,fill=group))+
  geom_bar(stat = "identity",position = "identity",width = 0.6)+
  scale_fill_manual(values = c("#3F5BA7","#D72326"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  theme(legend.position ="null")+
  theme(text = element_text(size = 10))+
  labs(fill="")+
  coord_flip()

#2.clt1/clt2 top40 decresed and top20 increased  *figure
species_Top60_clt1clt2=maaslin_species_clt1clt2 %>% filter(.,coef<0) %>% row.names(.)%>% .[1:40]
species_Top60_clt1clt2=maaslin_species_clt1clt2 %>% filter(.,coef>0) %>% row.names(.)%>% .[1:20] %>% c(.,species_Top60_clt1clt2)
newname=vector()
name=strsplit(species_Top60_clt1clt2,"s__")
for (i in 1:60) {
  newname[i]=name[[i]][2]
}
newname[8]="Bacteroides spp"
species_abs_clr_top60=species_abs_clr[c(PBCbase_cluster1,PBCbase_cluster2),species_Top60_clt1clt2]
colnames(species_abs_clr_top60)=newname

anno_sample=metadata_PBCbase[row.names(species_abs_clr_top60),"cluster"] %>% as.data.frame(.)
row.names(anno_sample)=row.names(species_abs_clr_top60)
colnames(anno_sample)="cluster"
anno_sample$cluster=anno_sample$cluster[,drop=TRUE]
str(anno_sample)
color3 = colorRampPalette(c("white","#2F5597"))(50)
ann_colors = list(
  cluster = c(clt1="#088288",clt2="#F19433")
)
library(pheatmap)
pheatmap(t(species_abs_clr_top60),  show_colnames = FALSE,
         fontsize = 11,cluster_cols = F,color = color3,annotation_col = anno_sample,annotation_colors = ann_colors,
         border=TRUE,gaps_col =71 ,scale = "row",breaks =  seq(0,4,length.out = 110))

#3.clt1 vs clt2升高最显著的15个以及下降最显著的10个species
species_Top25_clt1clt2=maaslin_species_clt1clt2 %>% filter(.,coef<0) %>% row.names(.)%>% .[1:15]
species_Top25_clt1clt2=maaslin_species_clt1clt2 %>% filter(.,coef>0) %>% row.names(.)%>% .[1:10] %>% c(.,species_Top25_clt1clt2)
newname=vector()
name=strsplit(species_Top25_clt1clt2,"s__")
for (i in 1:25) {
  newname[i]=name[[i]][2]
}
newname[8]="Bacteroides spp"


species_abs_clr_top25=species_abs_clr[c(HC_new,PBCbase_cluster1,PBCbase_cluster2),species_Top25_clt1clt2]
colnames(species_abs_clr_top25)=newname

anno_sample=metadata_PBCHC[row.names(species_abs_clr_top25),"cluster"] %>% as.data.frame(.)
row.names(anno_sample)=row.names(species_abs_clr_top25)
colnames(anno_sample)="cluster"
anno_sample$cluster=anno_sample$cluster[,drop=TRUE]
str(anno_sample)
color3 = colorRampPalette(c("white","#2F5597"))(50)

library(pheatmap)

#i.分两步提取order信息
newname=vector()
name=strsplit(species_Top25_clt1clt2,"o__")
for (i in 1:25) {
  newname[i]=name[[i]][2]
}
newname_v2=vector()
name=strsplit(newname,"f__")
for (i in 1:25) {
  newname_v2[i]=name[[i]][1]
}
newname_v2=as.factor(newname_v2)
summary(newname_v2)

annotation_row = data.frame(
  order = newname_v2
)
row.names(annotation_row)=colnames(species_abs_clr_top25)

ann_colors = list(
  order=c("Bacteroidales;"="#A4CFE0","Clostridiales;"="#1F77B6","Sphingomonadales;"="#B2E088","Veillonellales;"="#329F28","unidentified;"="#D9D9D9"),
  cluster = c(clt1="#088288",clt2="#F19433",HC="#3F5BA7")
)

pheatmap(t(species_abs_clr_top25),  show_colnames = FALSE,
         fontsize = 11,cluster_cols = F,color = color3,annotation_col  = anno_sample,annotation_row =annotation_row ,annotation_colors = ann_colors,
         border=TRUE,gaps_col =c(131,202) ,scale = "row",breaks =  seq(0,4,length.out = 110))

#4.clt1 vs clt2升高最显著的100个以及下降最显著的50个species,只标注部分species名称
species_Top150_clt1clt2=maaslin_species_clt1clt2 %>% filter(.,coef<0) %>% row.names(.)%>% .[1:100]
species_Top150_clt1clt2=maaslin_species_clt1clt2 %>% filter(.,coef>0) %>% row.names(.)%>% .[1:50] %>% c(.,species_Top150_clt1clt2)

species_abs_clr_top150=species_abs_clr[c(HC_new,PBCbase_cluster1,PBCbase_cluster2),species_Top150_clt1clt2]


anno_sample=metadata_PBCHC[row.names(species_abs_clr_top150),"cluster"] %>% as.data.frame(.)
row.names(anno_sample)=row.names(species_abs_clr_top150)
colnames(anno_sample)="cluster"
anno_sample$cluster=anno_sample$cluster[,drop=TRUE]
str(anno_sample)
color3 = colorRampPalette(c("white","#2F5597"))(50)

library(pheatmap)

#i.分两步提取order信息
newname=vector()
name=strsplit(species_Top150_clt1clt2,"o__")
for (i in 1:150) {
  newname[i]=name[[i]][2]
}
newname_v2=vector()
name=strsplit(newname,"f__")
for (i in 1:150) {
  newname_v2[i]=name[[i]][1]
}
newname_v2=as.factor(newname_v2)
summary(newname_v2)

#将species<5个的order及unidentified合并为other
newname_v2=replace(newname_v2,which(newname_v2=="Spirochaetales;"),"unidentified;")
newname_v2=replace(newname_v2,which(newname_v2=="Candidatus Borkfalkiales;"),"unidentified;")
newname_v2=replace(newname_v2,which(newname_v2=="Candidatus Brocadiales;"),"unidentified;")
newname_v2=replace(newname_v2,which(newname_v2=="Eggerthellales;"),"unidentified;")
newname_v2=replace(newname_v2,which(newname_v2=="Enterobacterales;"),"unidentified;")
newname_v2=replace(newname_v2,which(newname_v2=="Erysipelotrichales;"),"unidentified;")
newname_v2=replace(newname_v2,which(newname_v2=="Flavobacteriales;"),"unidentified;")
newname_v2=replace(newname_v2,which(newname_v2=="Fusobacteriales;"),"unidentified;")
newname_v2=replace(newname_v2,which(newname_v2=="Micrococcales;"),"unidentified;")
newname_v2=replace(newname_v2,which(newname_v2=="Pasteurellales;"),"unidentified;")
newname_v2=replace(newname_v2,which(newname_v2=="Rhodospirillales;"),"unidentified;")
newname_v2=replace(newname_v2,which(newname_v2=="Selenomonadales;"),"unidentified;")
newname_v2=replace(newname_v2,which(newname_v2=="Sphingomonadales;"),"unidentified;")
newname_v2=replace(newname_v2,which(newname_v2=="Thermoanaerobacterales;"),"unidentified;")

newname_v2=newname_v2[drop=TRUE] 


annotation_row = data.frame(
  order = newname_v2
)
row.names(annotation_row)=colnames(species_abs_clr_top150)

ann_colors = list(
  order=c("Bacteroidales;"="#A4CFE0","Clostridiales;"="#1F77B6","Veillonellales;"="#329F28","unidentified;"="#D9D9D9"),
  cluster = c(clt1="#088288",clt2="#F19433",HC="#3F5BA7")
)

pheatmap(t(species_abs_clr_top150),  show_colnames = FALSE,show_rownames = FALSE,
         fontsize = 11,cluster_cols = F,color = color3,annotation_col  = anno_sample,annotation_row =annotation_row,annotation_colors = ann_colors,
         border=TRUE,gaps_col =c(131,202) ,scale = "row",breaks =  seq(0,4,length.out = 110),cellwidth = 1,cellheight = 1)




#GMM信息及差异绘图 **figure=============================
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin")
maaslin_GMM_PBCHC=read.csv("maaslin_GMM PBCHC_new.csv",header = T,row.names = 1)
row.names(maaslin_GMM_PBCHC)=maaslin_GMM_PBCHC$feature
colnames(maaslin_GMM_PBCHC)
#1.PBCHC差异的GMM
maaslin_GMM_PBCHC_dif=maaslin_GMM_PBCHC %>% filter(.,qval<0.2) %>% row.names(.)
maaslin_GMM_PBCHC_dif_name=maaslin_GMM_PBCHC[maaslin_GMM_PBCHC_dif,"Name"]
GMM_abs_clr_dif=t(GMM_abs_clr[c(PBCbase,HC_new),maaslin_GMM_PBCHC_dif]) %>% as.data.frame(.)
row.names(GMM_abs_clr_dif)=maaslin_GMM_PBCHC_dif_name
library(pheatmap)
color3 = colorRampPalette(c("white","#661311"))(150)

pheatmap(GMM_abs_clr_dif, show_colnames =F,
         fontsize = 8,cluster_cols = F,cluster_rows = T,color = color3,annotation_col = anno_sample,annotation_colors = ann_colors,
         border=TRUE,scale="row",gaps_col = 132,breaks =  seq(-0.5,2,length.out = 150))

#2.clt1/clt2  有差异的GMM FDR<0.05
maaslin_GMM_clt1clt2_dif_name=maaslin_GMM_clt1clt2[maaslin_GMM_clt1clt2_dif,"Name"]
GMM_abs_clr_clt1clt2_dif=t(GMM_abs_clr[c(HC_new,PBCbase_cluster1,PBCbase_cluster2),maaslin_GMM_clt1clt2_dif]) %>% as.data.frame(.)
row.names(GMM_abs_clr_clt1clt2_dif)=maaslin_GMM_clt1clt2_dif_name
library(pheatmap)
color3 = colorRampPalette(c("white","#2F5597"))(150)

colnames(maaslin_GMM_clt1clt2)

pheatmap(GMM_abs_clr_clt1clt2_dif, show_colnames =F,
         fontsize = 8,cluster_cols = F,cluster_rows = T,color = color3,annotation_col = anno_sample,annotation_colors = ann_colors,
         border=TRUE,scale="row",gaps_col = c(131,202),breaks =  seq(-0.5,2,length.out = 150))

#2.1 clt1/clt2 差异最显著的top20 GMM  全样本热图 *figure
maaslin_GMM_clt1clt2_top20=row.names(maaslin_GMM_clt1clt2)[1:20]
maaslin_GMM_clt1clt2_top20_name=maaslin_GMM_clt1clt2[maaslin_GMM_clt1clt2_top20,"Name"]
GMM_abs_clr_clt1clt2_top20=t(GMM_abs_clr[c(HC_new,PBCbase_cluster1,PBCbase_cluster2),maaslin_GMM_clt1clt2_top20]) %>% as.data.frame(.)
row.names(GMM_abs_clr_clt1clt2_top20)=maaslin_GMM_clt1clt2_top20_name
library(pheatmap)
color3 = colorRampPalette(c("white","#2F5597"))(150)

colnames(maaslin_GMM_clt1clt2)

pheatmap(GMM_abs_clr_clt1clt2_top20, show_colnames =F,
         fontsize = 8,cluster_cols = F,cluster_rows = T,color = color3,annotation_col = anno_sample,annotation_colors = ann_colors,
         border=TRUE,scale="row",gaps_col = c(131,202),breaks =  seq(-0.5,2,length.out = 150),annotation_row = annotation_row)

#2.2 clt1/clt2 差异最显著的top30 GMM 中位值热图 *figure
GMM_clt1clt2_top30=row.names(maaslin_GMM_clt1clt2)[1:30]
GMM_abs_clr_clt1clt2_top30=GMM_abs_clr[c(PBCbase_cluster1,PBCbase_cluster2,HC_new),GMM_clt1clt2_top30]
GMM_clt1clt2_top30_plot=cbind(GMM_abs_clr_clt1clt2_top30,metadata_PBCHC[c(PBCbase_cluster1,PBCbase_cluster2,HC_new),"cluster"])
colnames(GMM_clt1clt2_top30_plot)[31]="cluster"
library(dplyr)
GMM_clt1clt2_top30_median=as.data.frame(matrix(0,ncol=1,nrow=3))
for (i in 1:30) {
  res.gmm=aggregate(GMM_clt1clt2_top30_plot[,i], by=list(type=GMM_clt1clt2_top30_plot$cluster),median)
  GMM_clt1clt2_top30_median[,i+1]=res.gmm$x
  row.names(GMM_clt1clt2_top30_median)=res.gmm$type
}
GMM_clt1clt2_top30_median=GMM_clt1clt2_top30_median[,-1]
colnames(GMM_clt1clt2_top30_median)=maaslin_GMM_clt1clt2$Name[1:30]
library(pheatmap)
pheatmap(t(GMM_clt1clt2_top30_median[c(3,1,2),]),cluster_rows = T,cluster_cols = F,border_color = "white",
         color =  colorRampPalette(c("#259042","white","#82388A" ))(50),
         main = "GMM",cellwidth = 18,cellheight = 10)


#2.3 上述top30 gmm在两两比较时的effect size 热图*figure#
#2.3.1.maasLin coefficient
GMM_clt1clt2_top30_clt1hc.coef=-maaslin_GMM_clt1HC[GMM_clt1clt2_top30,"coef"]
GMM_clt1clt2_top30_clt2hc.coef=-maaslin_GMM_clt2HC[GMM_clt1clt2_top30,"coef"]
GMM_clt1clt2_top30_clt1clt2.coef=maaslin_GMM_clt1clt2[GMM_clt1clt2_top30,"coef"]
GMM_clt1clt2_top30_coef=data.frame(
  clt1.hc=GMM_clt1clt2_top30_clt1hc.coef,
  clt2.hc=GMM_clt1clt2_top30_clt2hc.coef,
  clt1clt2=GMM_clt1clt2_top30_clt1clt2.coef
)
row.names(GMM_clt1clt2_top30_coef)=maaslin_GMM_clt1clt2$Name[1:30]
#2.3.2 上述top30 gmm在两两比较时的 fdr
GMM_clt1clt2_top30_clt1hc.qval=maaslin_GMM_clt1HC[GMM_clt1clt2_top30,"qval"]
GMM_clt1clt2_top30_clt2hc.qval=maaslin_GMM_clt2HC[GMM_clt1clt2_top30,"qval"]
GMM_clt1clt2_top30_clt1clt2.qval=maaslin_GMM_clt1clt2[GMM_clt1clt2_top30,"qval"]
GMM_clt1clt2_top30_qval=data.frame(
  clt1.hc=GMM_clt1clt2_top30_clt1hc.qval,
  clt2.hc=GMM_clt1clt2_top30_clt2hc.qval,
  clt1clt2=GMM_clt1clt2_top30_clt1clt2.qval
)
row.names(GMM_clt1clt2_top30_qval)=maaslin_GMM_clt1clt2$Name[1:30]

#结合中位相对丰度的热图调整顺序
index=c(9,30,21,8,15,5,6,22,3,24,25,12,11,27,13,23,
        10,28,7,17,14,20,1,2,19,29,26,4,16,18)

GMM_clt1clt2_top30_coef_plot=GMM_clt1clt2_top30_coef[index,]
GMM_clt1clt2_top30_qval_plot=GMM_clt1clt2_top30_qval[index,]

GMM_clt1clt2_top30_qval_plot[GMM_clt1clt2_top30_qval_plot<0.0001]="****"
GMM_clt1clt2_top30_qval_plot[GMM_clt1clt2_top30_qval_plot>=0.0001&GMM_clt1clt2_top30_qval_plot<0.001]="***"
GMM_clt1clt2_top30_qval_plot[GMM_clt1clt2_top30_qval_plot>=0.001&GMM_clt1clt2_top30_qval_plot<0.01]="**"
GMM_clt1clt2_top30_qval_plot[GMM_clt1clt2_top30_qval_plot>=0.01&GMM_clt1clt2_top30_qval_plot<0.05]="*"
GMM_clt1clt2_top30_qval_plot[GMM_clt1clt2_top30_qval_plot>=0.05]=""
GMM_clt1clt2_top30_qval_plot[28,2]="***" #因为“***”没有显示，手动调整

pheatmap(GMM_clt1clt2_top30_coef_plot,cluster_rows = F,cluster_cols = F,border_color = "white",
         color =  colorRampPalette(c("#2F5597","white","#890202" ))(50),display_numbers =GMM_clt1clt2_top30_qval_plot ,
         main = "GMM",breaks =  seq(-1.5,5.7,length.out = 110),scale="row",fontsize_number  = 11,cellwidth = 18,cellheight = 10) #按照row进行scale

#############小提琴图展示差异特征数目及系数分布 *figure###############
library(ggplot2) 
library(dplyr)
library(hrbrthemes)
library(viridis)
#1.species
species_violin <- data.frame(
  group=c(rep("clt1clt2",2182),rep("clt2HC",1993),rep("clt1HC",654),rep("PBCHC",294) ),
  MaAslin.coef=c( maaslin_species_clt1clt2[maaslin_species_clt1clt2_dif,"coef"], 
          -maaslin_species_clt2HC[maaslin_species_clt2HC_dif,"coef"], 
           -maaslin_species_clt1HC[maaslin_species_clt1HC_dif,"coef"], 
           maaslin_species_PBCHC[maaslin_species_PBCHC_dif,"coef"] )
)
species_violin$group=factor(species_violin$group,levels = c("PBCHC","clt1HC","clt2HC","clt1clt2"))


#绘制小提琴图像 
  ggplot( data=species_violin,aes(x=group, y=MaAslin.coef, fill=group)) +
    geom_violin(scale = "count",width=1.8) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE,option = "H") +
    theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("species") +
  xlab("")
  
#2.KO gene
  KO_gene_violin <- data.frame(
    group=c(rep("clt1clt2",3062),rep("clt2HC",1904),rep("clt1HC",362),rep("PBCHC",269) ),
    MaAslin.coef=c( maaslin_KO_gene_clt1clt2[maaslin_KO_gene_clt1clt2_dif,"coef"], 
                    -maaslin_KO_gene_clt2HC[maaslin_KO_gene_clt2HC_dif,"coef"], 
                    -maaslin_KO_gene_clt1HC[maaslin_KO_gene_clt1HC_dif,"coef"], 
                    maaslin_KO_gene_PBCHC[maaslin_KO_gene_PBCHC_dif,"coef"] )
  )
  KO_gene_violin$group=factor(KO_gene_violin$group,levels = c("PBCHC","clt1HC","clt2HC","clt1clt2"))
  
  
  #绘制小提琴图像 
  ggplot( data=KO_gene_violin,aes(x=group, y=MaAslin.coef, fill=group)) +
    geom_violin(scale = "count",width=1.8) +
    geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    scale_fill_viridis(discrete = TRUE,option = "G") +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("KO_gene") +
    xlab("")
  
  
  
#################cluster network construction V1########
###1.cluster1
#1.1 genera 筛选标准: 1.clt1clt2差异 2.clt1 mean relab>0.1% 3.非g__Unidentified
#计算cluster1中genus的平均相对丰度
genus_rela_mean_clt1=apply(genus_rela_filter[,PBCbase_cluster1], 1, mean) %>% as.data.frame(.)
colnames(genus_rela_mean_clt1)="mean_relab"
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network")
write.csv(genus_rela_mean_clt1,"genus_rela_mean_clt1.csv")
#平均相对丰度 >0.1% 的genus
genus_rela_clt1_0.001=filter(genus_rela_mean_clt1,mean_relab>0.001) %>% row.names(.)
#1.3 移除g__unidentified
genus_rela_clt1_0.001=genus_rela_clt1_0.001[-grep("g__unidentified",genus_rela_clt1_0.001)]
genus_rela_clt1_0.001=genus_rela_clt1_0.001[-39] #此genus为unclassified

#1.4 clt1 clt2有差异的genus
genus_rela_clt1_0.001=intersect(genus_rela_clt1_0.001,maaslin_genus_clt1clt2_dif) #18个genus
#1.5进行spearman分析
genus_abs_clr_network_clt1=genus_abs_clr[PBCbase_cluster1,genus_rela_clt1_0.001]
corR_genus_genus_clt1=sapply(1:18, function(x) sapply(1:18,function(y){cor.test(genus_abs_clr_network_clt1[,x],genus_abs_clr_network_clt1[,y],method = "spearman")$estimate}))
corP_genus_genus_clt1=sapply(1:18, function(x) sapply(1:18,function(y){cor.test(genus_abs_clr_network_clt1[,x],genus_abs_clr_network_clt1[,y],method = "spearman")$p.value}))
corq_genus_genus_clt1=p.adjust(corP_genus_genus_clt1,"fdr") %>% matrix(.,ncol = 18)

#修改network节点，只保留genus信息
network_genus_clt1_simpli=vector()
name=strsplit(genus_rela_clt1_0.001,"g__")
for (i in 1:18) {
  network_genus_clt1_simpli[i]=name[[i]][2]
}


row.names(corP_genus_genus_clt1)=network_genus_clt1_simpli  
row.names(corq_genus_genus_clt1)=network_genus_clt1_simpli  
row.names(corR_genus_genus_clt1)=network_genus_clt1_simpli  
colnames(corP_genus_genus_clt1)=network_genus_clt1_simpli  
colnames(corq_genus_genus_clt1)=network_genus_clt1_simpli  
colnames(corR_genus_genus_clt1)=network_genus_clt1_simpli  

#选择FDR<0.05的相关性
library(reshape)
library(reshape2)
corR_genus_genus_clt1_melt=melt(corR_genus_genus_clt1)
corq_genus_genus_clt1_melt=melt(corq_genus_genus_clt1)
index=which(corq_genus_genus_clt1_melt$value<0.05) #显著的相关性
corR_genus_genus_clt1_melt_sig=corR_genus_genus_clt1_melt[index,]
colnames(corR_genus_genus_clt1_melt_sig)[3]="corR"
corq_genus_genus_clt1_melt_sig=corq_genus_genus_clt1_melt[index,]
colnames(corq_genus_genus_clt1_melt_sig)[3]="corq"
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network/clt1")
genus_network_clt1=cbind(corR_genus_genus_clt1_melt_sig,corq_genus_genus_clt1_melt_sig)

write.csv(genus_network_clt1,"genus_network_clt1.csv")


#此18个network genus的相对丰度
network_genus_clt1_relab=genus_rela_mean_clt1[genus_rela_clt1_0.001,] %>% as.data.frame(.)
row.names(network_genus_clt1_relab)=network_genus_clt1_simpli
colnames(network_genus_clt1_relab)="relab"  
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network/clt1")
write.csv(network_genus_clt1_relab,"network_genus_clt1_relab.csv")

##2.cluster2
#计算cluster2中genus的平均相对丰度
genus_rela_mean_clt2=apply(genus_rela_filter[,PBCbase_cluster2], 1, mean) %>% as.data.frame(.)
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network")
write.csv(genus_rela_mean_clt2,"genus_rela_mean_clt2.csv")

#平均相对丰度 >0.1% 的genus
colnames(genus_rela_mean_clt2)="mean_relab"
genus_rela_clt2_0.001=filter(genus_rela_mean_clt2,mean_relab>0.001) %>% row.names(.)
#1.3 移除g__unidentified
genus_rela_clt2_0.001=genus_rela_clt2_0.001[-grep("g__unidentified",genus_rela_clt2_0.001)]
genus_rela_clt2_0.001=genus_rela_clt2_0.001[-38] #此genus为unclassified
#1.4 在clt1clt2有差异的genus
genus_rela_clt2_0.001=intersect(genus_rela_clt2_0.001,maaslin_genus_clt1clt2_dif) #13个genus
#1.5进行spearman分析
genus_abs_clr_network_clt2=genus_abs_clr[PBCbase_cluster2,genus_rela_clt2_0.001]
corR_genus_genus_clt2=sapply(1:13, function(x) sapply(1:13,function(y){cor.test(genus_abs_clr_network_clt2[,x],genus_abs_clr_network_clt2[,y],method = "spearman")$estimate}))
corP_genus_genus_clt2=sapply(1:13, function(x) sapply(1:13,function(y){cor.test(genus_abs_clr_network_clt2[,x],genus_abs_clr_network_clt2[,y],method = "spearman")$p.value}))
corq_genus_genus_clt2=p.adjust(corP_genus_genus_clt2,"fdr") %>% matrix(.,ncol = 13)

#修改network节点，只保留genus信息
network_genus_clt2_simpli=vector()
name=strsplit(genus_rela_clt2_0.001,"g__")
for (i in 1:13) {
  network_genus_clt2_simpli[i]=name[[i]][2]
}


row.names(corP_genus_genus_clt2)=network_genus_clt2_simpli  
row.names(corq_genus_genus_clt2)=network_genus_clt2_simpli  
row.names(corR_genus_genus_clt2)=network_genus_clt2_simpli  
colnames(corP_genus_genus_clt2)=network_genus_clt2_simpli  
colnames(corq_genus_genus_clt2)=network_genus_clt2_simpli  
colnames(corR_genus_genus_clt2)=network_genus_clt2_simpli  

#选择FDR<0.05的相关性
library(reshape)
library(reshape2)
corR_genus_genus_clt2_melt=melt(corR_genus_genus_clt2)
corq_genus_genus_clt2_melt=melt(corq_genus_genus_clt2)
index=which(corq_genus_genus_clt2_melt$value<0.05)
corR_genus_genus_clt2_melt_sig=corR_genus_genus_clt2_melt[index,]
colnames(corR_genus_genus_clt2_melt_sig)[3]="corR"
corq_genus_genus_clt2_melt_sig=corq_genus_genus_clt2_melt[index,]
colnames(corq_genus_genus_clt2_melt_sig)[3]="corq"
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network/clt2")
genus_network_clt2=cbind(corR_genus_genus_clt2_melt_sig,corq_genus_genus_clt2_melt_sig)
write.csv(genus_network_clt2,"genus_network_clt2.csv")

#此13个network genus的相对丰度
network_genus_clt2_relab=genus_rela_mean_clt2[genus_rela_clt2_0.001,] %>% as.data.frame(.)
row.names(network_genus_clt2_relab)=network_genus_clt2_simpli
colnames(network_genus_clt2_relab)="relab"  
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network/clt2")
write.csv(network_genus_clt2_relab,"network_genus_clt2_relab.csv")



##############cluster network construction V2 *figure##################
#1. genera 筛选: 1.clt1clt2差异 2.all sample mean relab>0.1% 3.移除order信息不明确的genus
#计算all sample 中genus的平均相对丰度
library(dplyr)
genus_rela_mean_allsample=apply(genus_rela_filter, 1, mean) %>% as.data.frame(.)
colnames(genus_rela_mean_allsample)="mean_relab"
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network v2")
write.csv(genus_rela_mean_allsample,"genus_rela_mean_allsample.csv")
#1.1 平均相对丰度 >0.1% 的genus
genus_rela_all_0.001=filter(genus_rela_mean_allsample,mean_relab>0.001) %>% row.names(.)
#1.2 clt1 clt2有差异的genus
genus_rela_all_0.001_clt.dif=intersect(genus_rela_all_0.001,maaslin_genus_clt1clt2_dif) #21个genus
#1.3 移除order信息不明确的genus
genus_rela_all_0.001_clt.dif=genus_rela_all_0.001_clt.dif[-c(19,21)]

#生成进行网络分析的genus simplified name，只保留genus信息
network_genus_simpli=vector()
name=strsplit(genus_rela_all_0.001_clt.dif,"g__")
for (i in 1:19) {
  network_genus_simpli[i]=name[[i]][2]
}
network_genus_simpli[10]="Oscillospiraceae.spp"
network_genus_simpli[15]="Ruminococcaceae.spp"
network_genus_simpli[16]="Clostridiales.spp"

#2. taxa spearman association
#2.1 clt1 
genus_abs_clr_network_clt1=genus_abs_clr[PBCbase_cluster1,genus_rela_all_0.001_clt.dif]
corR_genus_genus_clt1=sapply(1:19, function(x) sapply(1:19,function(y){cor.test(genus_abs_clr_network_clt1[,x],genus_abs_clr_network_clt1[,y],method = "spearman")$estimate}))
corP_genus_genus_clt1=sapply(1:19, function(x) sapply(1:19,function(y){cor.test(genus_abs_clr_network_clt1[,x],genus_abs_clr_network_clt1[,y],method = "spearman")$p.value}))
corq_genus_genus_clt1=p.adjust(corP_genus_genus_clt1,"fdr") %>% matrix(.,ncol = 19)

row.names(corP_genus_genus_clt1)=network_genus_simpli  
row.names(corq_genus_genus_clt1)=network_genus_simpli  
row.names(corR_genus_genus_clt1)=network_genus_simpli  
colnames(corP_genus_genus_clt1)=network_genus_simpli  
colnames(corq_genus_genus_clt1)=network_genus_simpli  
colnames(corR_genus_genus_clt1)=network_genus_simpli  

#选择FDR<0.05的相关性
library(reshape)
library(reshape2)
corR_genus_genus_clt1_melt=melt(corR_genus_genus_clt1)
corq_genus_genus_clt1_melt=melt(corq_genus_genus_clt1)
index=which(corq_genus_genus_clt1_melt$value<0.05) #显著的相关性
corR_genus_genus_clt1_melt_sig=corR_genus_genus_clt1_melt[index,]
colnames(corR_genus_genus_clt1_melt_sig)[3]="corR"
corq_genus_genus_clt1_melt_sig=corq_genus_genus_clt1_melt[index,]
colnames(corq_genus_genus_clt1_melt_sig)[3]="corq"
genus_network_clt1=cbind(corR_genus_genus_clt1_melt_sig,corq_genus_genus_clt1_melt_sig)
genus_network_clt1=genus_network_clt1[,-c(4,5)]
genus_network_clt1=genus_network_clt1[!duplicated(genus_network_clt1$corR),] %>% as.data.frame(.)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network v2/clt1")
write.csv(genus_network_clt1,"genus_network_clt1.csv")

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network v2/clt1")
write.csv(corP_genus_genus_clt1,"corP_genus_genus_clt1.csv")
write.csv(corq_genus_genus_clt1,"corq_genus_genus_clt1.csv")
write.csv(corR_genus_genus_clt1,"corR_genus_genus_clt1.csv")


#2.2 clt2
genus_abs_clr_network_clt2=genus_abs_clr[PBCbase_cluster2,genus_rela_all_0.001_clt.dif]
corR_genus_genus_clt2=sapply(1:19, function(x) sapply(1:19,function(y){cor.test(genus_abs_clr_network_clt2[,x],genus_abs_clr_network_clt2[,y],method = "spearman")$estimate}))
corP_genus_genus_clt2=sapply(1:19, function(x) sapply(1:19,function(y){cor.test(genus_abs_clr_network_clt2[,x],genus_abs_clr_network_clt2[,y],method = "spearman")$p.value}))
corq_genus_genus_clt2=p.adjust(corP_genus_genus_clt2,"fdr") %>% matrix(.,ncol = 19)

row.names(corP_genus_genus_clt2)=network_genus_simpli  
row.names(corq_genus_genus_clt2)=network_genus_simpli  
row.names(corR_genus_genus_clt2)=network_genus_simpli  
colnames(corP_genus_genus_clt2)=network_genus_simpli  
colnames(corq_genus_genus_clt2)=network_genus_simpli  
colnames(corR_genus_genus_clt2)=network_genus_simpli 

#选择FDR<0.05的相关性
library(reshape)
library(reshape2)
corR_genus_genus_clt2_melt=melt(corR_genus_genus_clt2)
corq_genus_genus_clt2_melt=melt(corq_genus_genus_clt2)
index=which(corq_genus_genus_clt2_melt$value<0.05)
corR_genus_genus_clt2_melt_sig=corR_genus_genus_clt2_melt[index,]
colnames(corR_genus_genus_clt2_melt_sig)[3]="corR"
corq_genus_genus_clt2_melt_sig=corq_genus_genus_clt2_melt[index,]
colnames(corq_genus_genus_clt2_melt_sig)[3]="corq"
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network v2/clt2")
genus_network_clt2=cbind(corR_genus_genus_clt2_melt_sig,corq_genus_genus_clt2_melt_sig)
genus_network_clt2=genus_network_clt2[,-c(4,5)]
genus_network_clt2=genus_network_clt2[!duplicated(genus_network_clt2$corR),] %>% as.data.frame(.)
write.csv(genus_network_clt2,"genus_network_clt2.csv")


setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network v2/clt2")
write.csv(corP_genus_genus_clt2,"corP_genus_genus_clt2.csv")
write.csv(corq_genus_genus_clt2,"corq_genus_genus_clt2.csv")
write.csv(corR_genus_genus_clt2,"corR_genus_genus_clt2.csv")

#生成node节点信息及变化信息

genus_network_node.info = data.frame(
  genus = genus_rela_all_0.001_clt.dif,
  simli.name=network_genus_simpli,
  coef=maaslin_genus_clt1clt2[genus_rela_all_0.001_clt.dif,"coef"]
)
row.names(genus_network_node.info)=network_genus_simpli
write.csv(genus_network_node.info,"genus_network_node.info.csv")
#2.3 PBC all
genus_abs_clr_network_PBCbase=genus_abs_clr[PBCbase,genus_rela_all_0.001_clt.dif]
corR_genus_genus_PBCbase=sapply(1:19, function(x) sapply(1:19,function(y){cor.test(genus_abs_clr_network_PBCbase[,x],genus_abs_clr_network_PBCbase[,y],method = "spearman")$estimate}))
corP_genus_genus_PBCbase=sapply(1:19, function(x) sapply(1:19,function(y){cor.test(genus_abs_clr_network_PBCbase[,x],genus_abs_clr_network_PBCbase[,y],method = "spearman")$p.value}))
corq_genus_genus_PBCbase=p.adjust(corP_genus_genus_PBCbase,"fdr") %>% matrix(.,ncol = 19)

row.names(corP_genus_genus_PBCbase)=network_genus_simpli  
row.names(corq_genus_genus_PBCbase)=network_genus_simpli  
row.names(corR_genus_genus_PBCbase)=network_genus_simpli  
colnames(corP_genus_genus_PBCbase)=network_genus_simpli  
colnames(corq_genus_genus_PBCbase)=network_genus_simpli  
colnames(corR_genus_genus_PBCbase)=network_genus_simpli 
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network v2/PBC")
write.csv(corP_genus_genus_PBCbase,"corP_genus_genus_PBCbase.csv")
write.csv(corq_genus_genus_PBCbase,"corq_genus_genus_PBCbase.csv")
write.csv(corR_genus_genus_PBCbase,"corR_genus_genus_PBCbase.csv")

#2.4 HC

genus_abs_clr_network_HC=genus_abs_clr[HC_new,genus_rela_all_0.001_clt.dif]
corR_genus_genus_HC=sapply(1:19, function(x) sapply(1:19,function(y){cor.test(genus_abs_clr_network_HC[,x],genus_abs_clr_network_HC[,y],method = "spearman")$estimate}))
corP_genus_genus_HC=sapply(1:19, function(x) sapply(1:19,function(y){cor.test(genus_abs_clr_network_HC[,x],genus_abs_clr_network_HC[,y],method = "spearman")$p.value}))
corq_genus_genus_HC=p.adjust(corP_genus_genus_HC,"fdr") %>% matrix(.,ncol = 19)

row.names(corP_genus_genus_HC)=network_genus_simpli  
row.names(corq_genus_genus_HC)=network_genus_simpli  
row.names(corR_genus_genus_HC)=network_genus_simpli  
colnames(corP_genus_genus_HC)=network_genus_simpli  
colnames(corq_genus_genus_HC)=network_genus_simpli  
colnames(corR_genus_genus_HC)=network_genus_simpli 

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network v2/HC")
write.csv(corP_genus_genus_HC,"corP_genus_genus_HC.csv")
write.csv(corq_genus_genus_HC,"corq_genus_genus_HC.csv")
write.csv(corR_genus_genus_HC,"corR_genus_genus_HC.csv")

#选择FDR<0.05的相关性
library(reshape)
library(reshape2)
corR_genus_genus_HC_melt=melt(corR_genus_genus_HC)
corq_genus_genus_HC_melt=melt(corq_genus_genus_HC)
index=which(corq_genus_genus_HC_melt$value<0.05)
corR_genus_genus_HC_melt_sig=corR_genus_genus_HC_melt[index,]
colnames(corR_genus_genus_HC_melt_sig)[3]="corR"
corq_genus_genus_HC_melt_sig=corq_genus_genus_HC_melt[index,]
colnames(corq_genus_genus_HC_melt_sig)[3]="corq"
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network v2/HC")
genus_network_HC=cbind(corR_genus_genus_HC_melt_sig,corq_genus_genus_HC_melt_sig)
genus_network_HC=genus_network_HC[,-c(4,5)]
genus_network_HC=genus_network_HC[!duplicated(genus_network_HC$corR),] %>% as.data.frame(.)
write.csv(genus_network_HC,"genus_network_HC.csv")

#2.5 all baseline sample
genus_abs_clr_network_allbase=genus_abs_clr[c(HC_new,PBCbase),genus_rela_all_0.001_clt.dif]
corR_genus_genus_allbase=sapply(1:19, function(x) sapply(1:19,function(y){cor.test(genus_abs_clr_network_allbase[,x],genus_abs_clr_network_allbase[,y],method = "spearman")$estimate}))
corP_genus_genus_allbase=sapply(1:19, function(x) sapply(1:19,function(y){cor.test(genus_abs_clr_network_allbase[,x],genus_abs_clr_network_allbase[,y],method = "spearman")$p.value}))
corq_genus_genus_allbase=p.adjust(corP_genus_genus_allbase,"fdr") %>% matrix(.,ncol = 19)

row.names(corP_genus_genus_allbase)=network_genus_simpli  
row.names(corq_genus_genus_allbase)=network_genus_simpli  
row.names(corR_genus_genus_allbase)=network_genus_simpli  
colnames(corP_genus_genus_allbase)=network_genus_simpli  
colnames(corq_genus_genus_allbase)=network_genus_simpli  
colnames(corR_genus_genus_allbase)=network_genus_simpli 

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network v2/PBC+HC")
write.csv(corP_genus_genus_allbase,"corP_genus_genus_allbase.csv")
write.csv(corq_genus_genus_allbase,"corq_genus_genus_allbase.csv")
write.csv(corR_genus_genus_allbase,"corR_genus_genus_allbase.csv")


###########cluster network construction V3####################
#1.只选用all_sample平均相对丰度 >0.1% 的genus,data来自network V2
genus_rela_all_0.001
#1.1移除order信息不明确的genus
genus_rela_all_0.001_no.unclassi=genus_rela_all_0.001[-c(39,47:49)]

#1.2生成进行网络分析的genus simplified name，只保留genus信息
network_genus_simpli=vector()
name=strsplit(genus_rela_all_0.001_no.unclassi,"g__")
for (i in 1:45) {
  network_genus_simpli[i]=name[[i]][2]
}
network_genus_simpli[10]="Bacteroidales.spp"
network_genus_simpli[23]="Lachnospiraceae.spp"
network_genus_simpli[25]="Oscillospiraceae.spp"
network_genus_simpli[32]="Ruminococcaceae.spp"
network_genus_simpli[33]="Clostridiales.spp"
network_genus_simpli[45]="Enterobacteriaceae.spp"

#2. taxa spearman association
#2.1 clt1 
genus_abs_clr_network_clt1=genus_abs_clr[PBCbase_cluster1,genus_rela_all_0.001_no.unclassi]
corR_genus_genus_clt1=sapply(1:45, function(x) sapply(1:45,function(y){cor.test(genus_abs_clr_network_clt1[,x],genus_abs_clr_network_clt1[,y],method = "spearman")$estimate}))
corP_genus_genus_clt1=sapply(1:45, function(x) sapply(1:45,function(y){cor.test(genus_abs_clr_network_clt1[,x],genus_abs_clr_network_clt1[,y],method = "spearman")$p.value}))
corq_genus_genus_clt1=p.adjust(corP_genus_genus_clt1,"fdr") %>% matrix(.,ncol = 45)

row.names(corP_genus_genus_clt1)=network_genus_simpli  
row.names(corq_genus_genus_clt1)=network_genus_simpli  
row.names(corR_genus_genus_clt1)=network_genus_simpli  
colnames(corP_genus_genus_clt1)=network_genus_simpli  
colnames(corq_genus_genus_clt1)=network_genus_simpli  
colnames(corR_genus_genus_clt1)=network_genus_simpli  

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network v3/clt1")
write.csv(corP_genus_genus_clt1,"corP_genus_genus_clt1.csv")
write.csv(corq_genus_genus_clt1,"corq_genus_genus_clt1.csv")
write.csv(corR_genus_genus_clt1,"corR_genus_genus_clt1.csv")

#2.2 clt2
genus_abs_clr_network_clt2=genus_abs_clr[PBCbase_cluster2,genus_rela_all_0.001_no.unclassi]
corR_genus_genus_clt2=sapply(1:45, function(x) sapply(1:45,function(y){cor.test(genus_abs_clr_network_clt2[,x],genus_abs_clr_network_clt2[,y],method = "spearman")$estimate}))
corP_genus_genus_clt2=sapply(1:45, function(x) sapply(1:45,function(y){cor.test(genus_abs_clr_network_clt2[,x],genus_abs_clr_network_clt2[,y],method = "spearman")$p.value}))
corq_genus_genus_clt2=p.adjust(corP_genus_genus_clt2,"fdr") %>% matrix(.,ncol = 45)

row.names(corP_genus_genus_clt2)=network_genus_simpli  
row.names(corq_genus_genus_clt2)=network_genus_simpli  
row.names(corR_genus_genus_clt2)=network_genus_simpli  
colnames(corP_genus_genus_clt2)=network_genus_simpli  
colnames(corq_genus_genus_clt2)=network_genus_simpli  
colnames(corR_genus_genus_clt2)=network_genus_simpli 
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network v3/clt2")
write.csv(corP_genus_genus_clt2,"corP_genus_genus_clt2.csv")
write.csv(corq_genus_genus_clt2,"corq_genus_genus_clt2.csv")
write.csv(corR_genus_genus_clt2,"corR_genus_genus_clt2.csv")

#2.3 PBC all
genus_abs_clr_network_PBCbase=genus_abs_clr[PBCbase,genus_rela_all_0.001_no.unclassi]
corR_genus_genus_PBCbase=sapply(1:45, function(x) sapply(1:45,function(y){cor.test(genus_abs_clr_network_PBCbase[,x],genus_abs_clr_network_PBCbase[,y],method = "spearman")$estimate}))
corP_genus_genus_PBCbase=sapply(1:45, function(x) sapply(1:45,function(y){cor.test(genus_abs_clr_network_PBCbase[,x],genus_abs_clr_network_PBCbase[,y],method = "spearman")$p.value}))
corq_genus_genus_PBCbase=p.adjust(corP_genus_genus_PBCbase,"fdr") %>% matrix(.,ncol = 45)

row.names(corP_genus_genus_PBCbase)=network_genus_simpli  
row.names(corq_genus_genus_PBCbase)=network_genus_simpli  
row.names(corR_genus_genus_PBCbase)=network_genus_simpli  
colnames(corP_genus_genus_PBCbase)=network_genus_simpli  
colnames(corq_genus_genus_PBCbase)=network_genus_simpli  
colnames(corR_genus_genus_PBCbase)=network_genus_simpli 
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network v3/PBC")
write.csv(corP_genus_genus_PBCbase,"corP_genus_genus_PBCbase.csv")
write.csv(corq_genus_genus_PBCbase,"corq_genus_genus_PBCbase.csv")
write.csv(corR_genus_genus_PBCbase,"corR_genus_genus_PBCbase.csv")

#2.4 HC

genus_abs_clr_network_HC=genus_abs_clr[HC_new,genus_rela_all_0.001_no.unclassi]
corR_genus_genus_HC=sapply(1:45, function(x) sapply(1:45,function(y){cor.test(genus_abs_clr_network_HC[,x],genus_abs_clr_network_HC[,y],method = "spearman")$estimate}))
corP_genus_genus_HC=sapply(1:45, function(x) sapply(1:45,function(y){cor.test(genus_abs_clr_network_HC[,x],genus_abs_clr_network_HC[,y],method = "spearman")$p.value}))
corq_genus_genus_HC=p.adjust(corP_genus_genus_HC,"fdr") %>% matrix(.,ncol = 45)

row.names(corP_genus_genus_HC)=network_genus_simpli  
row.names(corq_genus_genus_HC)=network_genus_simpli  
row.names(corR_genus_genus_HC)=network_genus_simpli  
colnames(corP_genus_genus_HC)=network_genus_simpli  
colnames(corq_genus_genus_HC)=network_genus_simpli  
colnames(corR_genus_genus_HC)=network_genus_simpli 
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network v3/HC")
write.csv(corP_genus_genus_HC,"corP_genus_genus_HC.csv")
write.csv(corq_genus_genus_HC,"corq_genus_genus_HC.csv")
write.csv(corR_genus_genus_HC,"corR_genus_genus_HC.csv")


#2.5 all baseline sample
genus_abs_clr_network_allbase=genus_abs_clr[c(HC_new,PBCbase),genus_rela_all_0.001_no.unclassi]
corR_genus_genus_allbase=sapply(1:45, function(x) sapply(1:45,function(y){cor.test(genus_abs_clr_network_allbase[,x],genus_abs_clr_network_allbase[,y],method = "spearman")$estimate}))
corP_genus_genus_allbase=sapply(1:45, function(x) sapply(1:45,function(y){cor.test(genus_abs_clr_network_allbase[,x],genus_abs_clr_network_allbase[,y],method = "spearman")$p.value}))
corq_genus_genus_allbase=p.adjust(corP_genus_genus_allbase,"fdr") %>% matrix(.,ncol = 45)

row.names(corP_genus_genus_allbase)=network_genus_simpli  
row.names(corq_genus_genus_allbase)=network_genus_simpli  
row.names(corR_genus_genus_allbase)=network_genus_simpli  
colnames(corP_genus_genus_allbase)=network_genus_simpli  
colnames(corq_genus_genus_allbase)=network_genus_simpli  
colnames(corR_genus_genus_allbase)=network_genus_simpli 

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network v3/PBC+HC")
write.csv(corP_genus_genus_allbase,"corP_genus_genus_allbase.csv")
write.csv(corq_genus_genus_allbase,"corq_genus_genus_allbase.csv")
write.csv(corR_genus_genus_allbase,"corR_genus_genus_allbase.csv")

##############cluster network construction V4##############
###1.cluster1
#1.1 genera 筛选标准: 1.clt1clt2差异 2.clt1 mean relab>0.1% 
#计算cluster1中genus的平均相对丰度
genus_rela_mean_clt1=apply(genus_rela_filter[,PBCbase_cluster1], 1, mean) %>% as.data.frame(.)
colnames(genus_rela_mean_clt1)="mean_relab"
#1.1平均相对丰度 >0.1% 的genus
genus_rela_clt1_0.001=filter(genus_rela_mean_clt1,mean_relab>0.001) %>% row.names(.)
#1.2 clt1 clt2有差异的genus
genus_rela_clt1_0.001=intersect(genus_rela_clt1_0.001,maaslin_genus_clt1clt2_dif) #22个genus
#1.3移除order信息不明确的genus
genus_rela_clt1_0.001=genus_rela_clt1_0.001[-c(20,23,25)]
#1.4进行spearman分析
genus_abs_clr_network_clt1=genus_abs_clr[PBCbase_cluster1,genus_rela_clt1_0.001]
corR_genus_genus_clt1=sapply(1:22, function(x) sapply(1:22,function(y){cor.test(genus_abs_clr_network_clt1[,x],genus_abs_clr_network_clt1[,y],method = "spearman")$estimate}))
corP_genus_genus_clt1=sapply(1:22, function(x) sapply(1:22,function(y){cor.test(genus_abs_clr_network_clt1[,x],genus_abs_clr_network_clt1[,y],method = "spearman")$p.value}))
corq_genus_genus_clt1=p.adjust(corP_genus_genus_clt1,"fdr") %>% matrix(.,ncol = 22)

cor.test(genus_abs_clr[PBCbase,"k__Bacteria;p__Firmicutes;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella"],
genus_abs_clr[PBCbase,"k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Coprococcus"],method = "spearman")

#修改network节点，只保留genus信息
network_genus_clt1_simpli=vector()
name=strsplit(genus_rela_clt1_0.001,"g__")
for (i in 1:22) {
  network_genus_clt1_simpli[i]=name[[i]][2]
}

network_genus_clt1_simpli[9]="Clostridiaceae spp"
network_genus_clt1_simpli[12]="Oscillospiraceae spp"
network_genus_clt1_simpli[18]="Ruminococcaceae spp"
network_genus_clt1_simpli[19]="o__Clostridiales spp"



row.names(corP_genus_genus_clt1)=network_genus_clt1_simpli  
row.names(corq_genus_genus_clt1)=network_genus_clt1_simpli  
row.names(corR_genus_genus_clt1)=network_genus_clt1_simpli  
colnames(corP_genus_genus_clt1)=network_genus_clt1_simpli  
colnames(corq_genus_genus_clt1)=network_genus_clt1_simpli  
colnames(corR_genus_genus_clt1)=network_genus_clt1_simpli  


#选择FDR<0.05的相关性
library(reshape)
library(reshape2)
corR_genus_genus_clt1_melt=melt(corR_genus_genus_clt1)
corq_genus_genus_clt1_melt=melt(corq_genus_genus_clt1)
index=which(corq_genus_genus_clt1_melt$value<0.05) #显著的相关性
corR_genus_genus_clt1_melt_sig=corR_genus_genus_clt1_melt[index,]
colnames(corR_genus_genus_clt1_melt_sig)[3]="corR"
corq_genus_genus_clt1_melt_sig=corq_genus_genus_clt1_melt[index,]
colnames(corq_genus_genus_clt1_melt_sig)[3]="corq"

genus_network_clt1=cbind(corR_genus_genus_clt1_melt_sig,corq_genus_genus_clt1_melt_sig)

#移除自身及互相重复的相关性
colnames(genus_network_clt1)
genus_network_clt1=genus_network_clt1[!duplicated(genus_network_clt1$corR),] %>% as.data.frame(.)
genus_network_clt1=genus_network_clt1[,-c(4,5)]

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network v4/clt1")
write.csv(genus_network_clt1,"genus_network_clt1.csv")
write.csv(corP_genus_genus_clt1,"corP_genus_genus_clt1.csv")
write.csv(corq_genus_genus_clt1,"corq_genus_genus_clt1.csv")
write.csv(corR_genus_genus_clt1,"corR_genus_genus_clt1.csv")

##2.cluster2
#计算cluster2中genus的平均相对丰度
genus_rela_mean_clt2=apply(genus_rela_filter[,PBCbase_cluster2], 1, mean) %>% as.data.frame(.)
#2.1平均相对丰度 >0.1% 的genus
colnames(genus_rela_mean_clt2)="mean_relab"
genus_rela_clt2_0.001=filter(genus_rela_mean_clt2,mean_relab>0.001) %>% row.names(.)
#2.2 在clt1clt2有差异的genus
genus_rela_clt2_0.001=intersect(genus_rela_clt2_0.001,maaslin_genus_clt1clt2_dif) 
#2.3 移除order不明确的物种
genus_rela_clt2_0.001=genus_rela_clt2_0.001[-c(14,16)]  #14个genus

#2.4进行spearman分析
genus_abs_clr_network_clt2=genus_abs_clr[PBCbase_cluster2,genus_rela_clt2_0.001]
corR_genus_genus_clt2=sapply(1:14, function(x) sapply(1:14,function(y){cor.test(genus_abs_clr_network_clt2[,x],genus_abs_clr_network_clt2[,y],method = "spearman")$estimate}))
corP_genus_genus_clt2=sapply(1:14, function(x) sapply(1:14,function(y){cor.test(genus_abs_clr_network_clt2[,x],genus_abs_clr_network_clt2[,y],method = "spearman")$p.value}))
corq_genus_genus_clt2=p.adjust(corP_genus_genus_clt2,"fdr") %>% matrix(.,ncol = 14)

#修改network节点，只保留genus信息
network_genus_clt2_simpli=vector()
name=strsplit(genus_rela_clt2_0.001,"g__")
for (i in 1:14) {
  network_genus_clt2_simpli[i]=name[[i]][2]
}

network_genus_clt2_simpli[11]="Clostridiales spp"
  
row.names(corP_genus_genus_clt2)=network_genus_clt2_simpli  
row.names(corq_genus_genus_clt2)=network_genus_clt2_simpli  
row.names(corR_genus_genus_clt2)=network_genus_clt2_simpli  
colnames(corP_genus_genus_clt2)=network_genus_clt2_simpli  
colnames(corq_genus_genus_clt2)=network_genus_clt2_simpli  
colnames(corR_genus_genus_clt2)=network_genus_clt2_simpli  

#选择FDR<0.05的相关性
library(reshape)
library(reshape2)
corR_genus_genus_clt2_melt=melt(corR_genus_genus_clt2)
corq_genus_genus_clt2_melt=melt(corq_genus_genus_clt2)
index=which(corq_genus_genus_clt2_melt$value<0.05)
corR_genus_genus_clt2_melt_sig=corR_genus_genus_clt2_melt[index,]
colnames(corR_genus_genus_clt2_melt_sig)[3]="corR"
corq_genus_genus_clt2_melt_sig=corq_genus_genus_clt2_melt[index,]
colnames(corq_genus_genus_clt2_melt_sig)[3]="corq"
genus_network_clt2=cbind(corR_genus_genus_clt2_melt_sig,corq_genus_genus_clt2_melt_sig)


#移除自身及互相重复的相关性
colnames(genus_network_clt2)
genus_network_clt2=genus_network_clt2[!duplicated(genus_network_clt2$corR),] %>% as.data.frame(.)
genus_network_clt2=genus_network_clt2[,-c(4,5)]

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/genus network v4/clt2")
write.csv(genus_network_clt2,"genus_network_clt2.csv")
write.csv(corP_genus_genus_clt2,"corP_genus_genus_clt2.csv")
write.csv(corq_genus_genus_clt2,"corq_genus_genus_clt2.csv")
write.csv(corR_genus_genus_clt2,"corR_genus_genus_clt2.csv")

#############PBC or cluster治疗前后变化 连线图 **figure################
#1.alpha多样性
#1.1 PBC总体
data=data.frame(baseline=shannon[PBCpair_6m,1],treat=shannon[PBC_6m_pair,1]) #shannon.base|6m为向量

library(ggpubr)
ggpaired(data, cond1 = "baseline", cond2 = "treat",
         color = "condition", palette = c("#D72326","#F29E9C"),ylab = "shannon",line.color = "grey")+
  theme(legend.position = "none") 

#1.2 cluster1
wilcox.test(shannon[PBCpair_6m_cluster1,1],shannon[PBC_6m_cluster1,1],paired = TRUE)
data=data.frame(baseline=shannon[PBCpair_6m_cluster1,1],treat=shannon[PBC_6m_cluster1,1]) #shannon.base|6m为向量

library(ggpubr)
ggpaired(data, cond1 = "baseline", cond2 = "treat",
         color = "condition", palette = c("#088288","#ADD5DB"),ylab = "shannon",line.color = "grey")+
  theme(legend.position = "none") 

#1.3 cluster2
wilcox.test(shannon[PBCpair_6m_cluster2,1],shannon[PBC_6m_cluster2,1],paired = TRUE)

data=data.frame(baseline=shannon[PBCpair_6m_cluster2,1],treat=shannon[PBC_6m_cluster2,1]) #shannon.base|6m为向量

library(ggpubr)
ggpaired(data, cond1 = "baseline", cond2 = "treat",
         color = "condition", palette = c("#F19433","#FFD966"),ylab = "shannon",line.color = "grey")+
  theme(legend.position = "none") 

#2.microbial gene
#2.2 cluster1
metadata[PBC_6m,clinic]
wilcox.test(metadata[PBCpair_6m_cluster1,"ALP"],metadata[PBC_6m_cluster1,"ALP"],paired = TRUE)

library(ggpubr)
data=data.frame(baseline=KO_gene_abs_clr_plot_treat[PBCpair_6m_cluster1,"K00956"],
                treat=KO_gene_abs_clr_plot_treat[PBC_6m_cluster1,"K00956"]) 
ggpaired(data, cond1 = "baseline", cond2 = "treat",
         color = "condition",palette = c("#088288","#ADD5DB"),ylab = "cysN",line.color = "grey")+
  theme(legend.position = "none") 

#2.3 cluster2

library(ggpubr)
data=data.frame(baseline=KO_gene_abs_clr_plot_treat[PBCpair_6m_cluster2,"K00956"],
                treat=KO_gene_abs_clr_plot_treat[PBC_6m_cluster2,"K00956"]) 
ggpaired(data, cond1 = "baseline", cond2 = "treat",
         color = "condition",palette = c("#F19433","#FFD966"),ylab = "cysN",line.color = "grey")+
      theme(legend.position = "none") 

#2.4 PBCbase
library(ggpubr)
data=data.frame(baseline=KO_gene_abs_clr_plot_treat[PBCpair_6m,"K00955"],
                treat=KO_gene_abs_clr_plot_treat[PBC_6m_pair,"K00955"]) 
ggpaired(data, cond1 = "baseline", cond2 = "treat",
         color = "condition",palette =c("#D72326","#F29E9C"),ylab = "cysN",line.color = "grey")+
  theme(legend.position = "none") 

#############7.cluster治疗前后变化  HC_new *figure################
##I boxplot
#1.KO_gene plot with HC
KO_gene_abs_clr_plot_treat=KO_gene_abs_clr[c(PBCpair_6m_cluster1,PBC_6m_cluster1,PBCpair_6m_cluster2,PBC_6m_cluster2,HC_new),row.names(KO_gene_rela_filter)]
KO_gene_abs_clr_plot_treat$group=c(rep("clt1.base",36),rep("clt1.6m",36),rep("clt2.base",23),rep("clt2.6m",23),rep("HC",131))
row.names(KO_gene_abs_clr_plot_treat)
KO_gene_abs_clr_plot_treat$group=factor(KO_gene_abs_clr_plot_treat$group,levels = c("clt1.base","clt1.6m","clt2.base","clt2.6m","HC"))
library(ggplot2)
ggplot(data = KO_gene_abs_clr_plot_treat,aes(x=group,y=KO_gene_abs_clr_plot_treat$K01023,fill=group))+geom_boxplot()+
  labs(y="assT (CLR)")+
  scale_fill_manual(values = c("#088288","#ADD5DB","#F19433","#FFD966","#3F5BA7"))+
  theme(legend.position = "none")

#1.2 KO_gene plot without HC
KO_gene_abs_clr_plot_treat_noHC=KO_gene_abs_clr[c(PBCpair_6m_cluster1,PBC_6m_cluster1,PBCpair_6m_cluster2,PBC_6m_cluster2),row.names(KO_gene_rela_filter)]
KO_gene_abs_clr_plot_treat_noHC$group=c(rep("clt1.base",36),rep("clt1.6m",36),rep("clt2.base",23),rep("clt2.6m",23))

KO_gene_abs_clr_plot_treat_noHC$group=factor(KO_gene_abs_clr_plot_treat_noHC$group,levels = c("clt1.base","clt1.6m","clt2.base","clt2.6m"))
library(ggplot2)
ggplot(data = KO_gene_abs_clr_plot_treat_noHC,aes(x=group,y=KO_gene_abs_clr_plot_treat_noHC$K01023,fill=group))+geom_boxplot()+
  labs(y="assT (CLR)")+
  scale_fill_manual(values = c("#088288","#ADD5DB","#F19433","#FFD966"))+
  theme(legend.position = "none")

#2.species
species_abs_clr_plot_treat=species_abs_clr[c(PBCpair_6m_cluster1,PBC_6m_cluster1,PBCpair_6m_cluster2,PBC_6m_cluster2,HC_new),row.names(species_rela_filter)]
species_abs_clr_plot_treat$group=c(rep("clt1.base",36),rep("clt1.6m",36),rep("clt2.base",23),rep("clt2.6m",23),rep("HC",131))
row.names(species_abs_clr_plot_treat)
species_abs_clr_plot_treat$group=factor(species_abs_clr_plot_treat$group,levels = c("clt1.base","clt1.6m","clt2.base","clt2.6m","HC"))


#######7.1 筛选clt1治疗前后变化最大的top40 species *heatmap########
#2.1.with HC
library(dplyr)
library(plyr)
res_species_pair.wilcox_clt1$species=row.names(res_species_pair.wilcox_clt1)
res_species_pair.wilcox_clt1=arrange(res_species_pair.wilcox_clt1,-desc(res_species_pair.wilcox_clt1$clt1base_6m.qval))
row.names(res_species_pair.wilcox_clt1)=res_species_pair.wilcox_clt1$species
species_clt1base.6m.dif.top40=row.names(res_species_pair.wilcox_clt1)[1:40]

species_abs_clr_plot_treat_dif=species_abs_clr_plot_treat[,c(species_clt1base.6m.dif.top40,"group")]
species_abs_clr_plot_treat_dif$group=species_abs_clr_plot_treat_dif$group[,drop=TRUE]

species_abs_clr_plot_treat_dif_mean=as.data.frame(matrix(0,ncol=1,nrow=5))

for (i in 1:40) {
  res.species=aggregate(species_abs_clr_plot_treat_dif[,i], by=list(type=species_abs_clr_plot_treat_dif$group),mean)
  species_abs_clr_plot_treat_dif_mean[,i+1]=res.species$x
  row.names(species_abs_clr_plot_treat_dif_mean)=res.species$type
}
species_abs_clr_plot_treat_dif_mean=species_abs_clr_plot_treat_dif_mean[,-1]
colnames(species_abs_clr_plot_treat_dif_mean)=species_clt1base.6m.dif.top40

newname=vector()
name=strsplit(species_clt1base.6m.dif.top40,"s__")
for (i in 1:40) {
  newname[i]=name[[i]][2]
}
colnames(species_abs_clr_plot_treat_dif_mean)=newname


#add order annotation
colnames(species_abs_clr_plot_treat_dif_mean)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/plot/plot_source_data")

#i.分两步提取order信息
newname=vector()
name=strsplit(species_clt1base.6m.dif.top40,"o__")
for (i in 1:40) {
  newname[i]=name[[i]][2]
}
newname_v2=vector()
name=strsplit(newname,"f__")
for (i in 1:40) {
  newname_v2[i]=name[[i]][1]
}

annotation_col = data.frame(
 order = factor(newname_v2)
)

for (i in c(2,8,14,15,22,28)) {
  newname_v2[i]="other"
}
newname_v2[2]="other"
newname_v2=as.factor(newname_v2)
summary(newname_v2)

annotation_col = data.frame(
 order = newname_v2
)
row.names(annotation_col)=colnames(species_abs_clr_plot_treat_dif_mean)
#将order信息中unidentified和other合并为”other“
which(annotation_col$order=="unidentified;") #1,4,113,21
annotation_col[21,1]="other"

ann_colors=list(order=c("Bacteroidales;"="#FCBD7E","Clostridiales;"="#5E8EC2",other="#FDF9CD"))
library(pheatmap)
pheatmap(t(species_abs_clr_plot_treat_dif_mean),cluster_rows = T,cluster_cols = F,scale = "row",
         color = colorRampPalette(c("#2F5597","white","#890202" ))(50),gaps_row = c(2,4),angle_col = "45",border_color = NA,fontsize = 7,
         annotation_row  = annotation_col,annotation_colors = ann_colors,cellwidth = 12,cellheight = 7)

#2.2 without HC

species_abs_clr_plot_treat_dif=species_abs_clr_plot_treat[1:118,c(species_clt1base.6m.dif.top40,"group")]
species_abs_clr_plot_treat_dif$group=species_abs_clr_plot_treat_dif$group[,drop=TRUE]

species_abs_clr_plot_treat_dif_mean=as.data.frame(matrix(0,ncol=1,nrow=4))

for (i in 1:40) {
  res.species=aggregate(species_abs_clr_plot_treat_dif[,i], by=list(type=species_abs_clr_plot_treat_dif$group),mean)
  species_abs_clr_plot_treat_dif_mean[,i+1]=res.species$x
  row.names(species_abs_clr_plot_treat_dif_mean)=res.species$type
}
species_abs_clr_plot_treat_dif_mean=species_abs_clr_plot_treat_dif_mean[,-1]
colnames(species_abs_clr_plot_treat_dif_mean)=species_clt1base.6m.dif.top40

newname=vector()
name=strsplit(species_clt1base.6m.dif.top40,"s__")
for (i in 1:40) {
  newname[i]=name[[i]][2]
}
colnames(species_abs_clr_plot_treat_dif_mean)=newname

annotation_col = data.frame(
  order = newname_v2
)
row.names(annotation_col)=colnames(species_abs_clr_plot_treat_dif_mean)

ann_colors=list(order=c("Bacteroidales;"="#FCBD7E","Clostridiales;"="#5E8EC2",other="#FDF9CD","unidentified;"="#BCE0EC"))
library(pheatmap)
pheatmap(species_abs_clr_plot_treat_dif_mean,cluster_rows = F,cluster_cols = T,scale = "column",
         color = colorRampPalette(c("#2F5597","white","#890202" ))(50),gaps_row = c(2,4),angle_col = "45",border_color = NA,fontsize = 7,
         annotation_col = annotation_col,annotation_colors = ann_colors)


#######7.2筛选clt2治疗前后变化top11 species (FDR<0.2) *heatmap#############
library(dplyr)
library(plyr)
res_species_pair.wilcox_clt2$species=row.names(res_species_pair.wilcox_clt2)
res_species_pair.wilcox_clt2=arrange(res_species_pair.wilcox_clt2,-desc(res_species_pair.wilcox_clt2$clt2base_6m.qval))
row.names(res_species_pair.wilcox_clt2)=res_species_pair.wilcox_clt2$species
species_clt2base.6m.dif.top11=row.names(res_species_pair.wilcox_clt2)[1:11]

species_abs_clr_plot_treat_dif_clt2=species_abs_clr_plot_treat[1:118,c(species_clt2base.6m.dif.top11,"group")]
species_abs_clr_plot_treat_dif_clt2$group=species_abs_clr_plot_treat_dif_clt2$group[,drop=TRUE]

species_abs_clr_plot_treat_dif_clt2_mean=as.data.frame(matrix(0,ncol=1,nrow=4))

for (i in 1:11) {
  res.species=aggregate(species_abs_clr_plot_treat_dif_clt2[,i], by=list(type=species_abs_clr_plot_treat_dif_clt2$group),mean)
  species_abs_clr_plot_treat_dif_clt2_mean[,i+1]=res.species$x
  row.names(species_abs_clr_plot_treat_dif_clt2_mean)=res.species$type
}
species_abs_clr_plot_treat_dif_clt2_mean=species_abs_clr_plot_treat_dif_clt2_mean[,-1]
colnames(species_abs_clr_plot_treat_dif_clt2_mean)=species_clt2base.6m.dif.top11

newname=vector()
name=strsplit(species_clt2base.6m.dif.top11,"s__")
for (i in 1:11) {
  newname[i]=name[[i]][2]
}
colnames(species_abs_clr_plot_treat_dif_clt2_mean)=newname

pheatmap(t(species_abs_clr_plot_treat_dif_clt2_mean),cluster_rows = T,cluster_cols = F,scale = "row",
         color = colorRampPalette(c("#2F5597","white","#890202" ))(50),gaps_row = c(2,4),angle_col = "45",border_color = NA,fontsize = 7,
         cellwidth = 12,cellheight = 7)
dev.off()
dev.new()

#############7.3 筛选PBC总体治疗前后变化最大的top50 species **figure###########
#2.1.with HC
library(dplyr)
library(plyr)
res_species_pair.wilcox$species=row.names(res_species_pair.wilcox)
res_species_pair.wilcox=arrange(res_species_pair.wilcox,-desc(res_species_pair.wilcox$PBCbase_6m.pval))
row.names(res_species_pair.wilcox)=res_species_pair.wilcox$species
species_PBCbase.6m.dif.top50=row.names(res_species_pair.wilcox)[1:50]

species_abs_clr_plot_treat_dif_PBC=species_abs_clr_plot_treat[,c(species_PBCbase.6m.dif.top50,"group")]
species_abs_clr_plot_treat_dif_PBC$group=species_abs_clr_plot_treat_dif_PBC$group[,drop=TRUE]

species_abs_clr_plot_treat_dif_PBC_mean=as.data.frame(matrix(0,ncol=1,nrow=3))

for (i in 1:50) {
  res.species=aggregate(species_abs_clr_plot_treat_dif_PBC[,i], by=list(type=species_abs_clr_plot_treat_dif_PBC$group),mean)
  species_abs_clr_plot_treat_dif_PBC_mean[,i+1]=res.species$x
  row.names(species_abs_clr_plot_treat_dif_PBC_mean)=res.species$type
}
species_abs_clr_plot_treat_dif_PBC_mean=species_abs_clr_plot_treat_dif_PBC_mean[,-1]
colnames(species_abs_clr_plot_treat_dif_PBC_mean)=species_PBCbase.6m.dif.top50

newname=vector()
name=strsplit(species_PBCbase.6m.dif.top50,"s__")
for (i in 1:50) {
  newname[i]=name[[i]][2]
}
colnames(species_abs_clr_plot_treat_dif_PBC_mean)=newname
species_abs_clr_plot_treat_dif_PBC_mean=species_abs_clr_plot_treat_dif_PBC_mean[c(3,2,1),]
#add order annotation

#i.分两步提取order信息
newname_V1=vector()
name=strsplit(species_PBCbase.6m.dif.top50,"o__")
for (i in 1:50) {
  newname_V1[i]=name[[i]][2]
}
newname_v2=vector()
name=strsplit(newname_V1,"f__")
for (i in 1:50) {
  newname_v2[i]=name[[i]][1]
}
#将包含细菌少于5个的class合并为other
which(newname_v2=="unidentified;")
for (i in c(11,21,49,18,41,37,19,22,45,48,27,32,1,5,7,17)) {
  newname_v2[i]="other"
}
newname_v2=as.factor(newname_v2)

summary(newname_v2)

annotation_col = data.frame(
  order = newname_v2
)

row.names(annotation_col)=colnames(species_abs_clr_plot_treat_dif_PBC_mean)

ann_colors=list(order=c("Bacteroidales;"="#FCBD7E","Clostridiales;"="#5E8EC2","Veillonellales;"="#B49800",other="#FDF9CD"))

library(pheatmap)
pheatmap(t(species_abs_clr_plot_treat_dif_PBC_mean),cluster_rows = T,cluster_cols = F,scale = "row",
         color = colorRampPalette(c("#259042", "white", "#82388A" ))(50),gaps_row = c(2,4),angle_col = "45",border_color = "white",fontsize = 7,
         annotation_row  = annotation_col,annotation_colors = ann_colors,cellwidth = 12,cellheight = 7)



###############7.5 genus level plot################
genus_abs_clr_plot_treat=genus_abs_clr[c(PBCpair_6m_cluster1,PBC_6m_cluster1,PBCpair_6m_cluster2,PBC_6m_cluster2,HC_new),row.names(genus_rela_filter)]
genus_abs_clr_plot_treat$group=c(rep("clt1.base",36),rep("clt1.6m",36),rep("clt2.base",23),rep("clt2.6m",23),rep("HC",131))
row.names(genus_abs_clr_plot_treat)
genus_abs_clr_plot_treat$group=factor(genus_abs_clr_plot_treat$group,levels = c("clt1.base","clt1.6m","clt2.base","clt2.6m","HC"))
library(ggplot2)
ggplot(data = genus_abs_clr_plot_treat,aes(x=group,y=genus_abs_clr_plot_treat$"k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides",fill=group))+geom_boxplot()+
  labs(y="Bacteroides (CLR)")+
  scale_fill_manual(values = c("#088288","#ADD5DB","#F19433","#FFD966","#3F5BA7"))+
  theme(legend.position = "none")

##II 治疗前后按照热图展示
#筛选PBC治疗前后变化FDR<0.2且丰度大于0.1%的genus
genus_PBCbase.6m.dif=filter(res_genus_pair.wilcox,PBCbase_6m.qval<0.2 & PBC.abundance>0.001) %>% row.names(.)
genus_abs_clr_plot_treat_dif=genus_abs_clr_plot_treat[1:118,c(genus_PBCbase.6m.dif,"group")]
genus_abs_clr_plot_treat_dif$group=genus_abs_clr_plot_treat_dif$group[,drop=TRUE]

genus_abs_clr_plot_treat_dif_mean=as.data.frame(matrix(0,ncol=1,nrow=4))

for (i in 1:28) {
  res.genus=aggregate(genus_abs_clr_plot_treat_dif[,i], by=list(type=genus_abs_clr_plot_treat_dif$group),mean)
  genus_abs_clr_plot_treat_dif_mean[,i+1]=res.genus$x
  row.names(genus_abs_clr_plot_treat_dif_mean)=res.genus$type
}
genus_abs_clr_plot_treat_dif_mean=genus_abs_clr_plot_treat_dif_mean[,-1]
colnames(genus_abs_clr_plot_treat_dif_mean)=genus_PBCbase.6m.dif

newname=vector()
name=strsplit(genus_PBCbase.6m.dif,"g__")
for (i in 1:28) {
  newname[i]=name[[i]][2]
}
newname[6]="Bacteroidales spp"
newname[14]="Lachnospiraceae spp"
newname[15]="Oscillospiraceae spp"
newname[21]="Ruminococcaceae spp"
newname[22]="Clostridiales spp"
newname[25]="Firmicutes spp"
newname[28]="unidentified spp"
colnames(genus_abs_clr_plot_treat_dif_mean)=newname

library(pheatmap)
pheatmap(genus_abs_clr_plot_treat_dif_mean,cluster_rows = F,cluster_cols = T,scale = "column",
         color = colorRampPalette(c( "white", "firebrick3"))(50),gaps_row = 2,angle_col = "45")



###############计算veillonella genus与clostridial的比值及其临床相关性 *figure######################
Veillonella=genus_rela_filter["k__Bacteria;p__Firmicutes;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella",c(PBCbase,HC_new,PBC_6m_pair)]
Clostridial=order_rela_filter["k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales",c(PBCbase,HC_new,PBC_6m_pair)]
Vei_clos=rbind(Veillonella,Clostridial)
row.names(Vei_clos)
colnames(Vei_clos)
Vei_clos=t(Vei_clos) %>% as.data.frame(.)
Vei_clos$clos_vei_ratio=Vei_clos[,2]/Vei_clos[,1]
wilcox.test(Vei_clos[PBCbase_resp_6m,3],Vei_clos[PBCbase_noresp_6m,3])
wilcox.test(Vei_clos[PBCbase_cluster1,3],Vei_clos[PBCbase_cluster2,3])
wilcox.test(Vei_clos[PBCbase_noresp_6m,3],Vei_clos[HC_new,3])
wilcox.test(Vei_clos[PBCpair_6m,3],Vei_clos[PBC_6m_pair,3],paired = T)

boxplot(Vei_clos[PBCbase_resp_6m,3],Vei_clos[PBCbase_noresp_6m,3])

Vei_clos=cbind(Vei_clos,metadata_PBCHC[row.names(Vei_clos),c("cluster","group","Response_6m")])


#1.按照 response 展示
Vei_clos_response=Vei_clos[c(PBCbase_resp_6m,PBCbase_noresp_6m,HC_new),3] %>% as.data.frame(.)
row.names(Vei_clos_response)=c(PBCbase_resp_6m,PBCbase_noresp_6m,HC_new)
Vei_clos_response$response=c(rep("response",92),rep("no-response",39),rep("HC",131))
colnames(Vei_clos_response)[1]="Clos.Vei.ratio"
Vei_clos_response$response=factor(Vei_clos_response$response,levels = c("HC","response","no-response"))
library(ggplot2)

ggplot(data = Vei_clos_response,aes(x=response,y=Vei_clos_response$Clos.Vei.ratio,fill=response))+geom_boxplot()+
  labs(y="Clost/Veil (Rela)")+
  scale_fill_manual(values = c("#3F5BA7","#A6DCE4","#C7C7C7"))+
  theme(legend.position = "none")

#2.按照cluster展示
Vei_clos_cluster=Vei_clos[c(PBCbase_cluster1,PBCbase_cluster2,HC_new),]
Vei_clos_cluster$cluster=factor(Vei_clos_cluster$cluster,levels = c("HC","clt1","clt2"))

ggplot(data = Vei_clos_cluster,aes(x=cluster,y=Vei_clos_cluster$clos_vei_ratio,fill=cluster))+geom_boxplot()+
  labs(y="Clost/Veil (Rela)")+
  scale_fill_manual(values = c("#3F5BA7","#088288","#F19433"))+
  theme(legend.position = "none")

############## clostridial与Veillonella的相关性计算###################
colnames(Vei_clos)
row.names(Vei_clos)
cor.test(Vei_clos[PBCbase,]$`k__Bacteria;p__Firmicutes;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella`,
         Vei_clos[PBCbase,]$`k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales`,method = "spearman") #P=0.36, rho=-0.079
cor.test(Vei_clos[PBCbase_cluster1,]$`k__Bacteria;p__Firmicutes;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella`,
         Vei_clos[PBCbase_cluster1,]$`k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales`,method = "spearman") #P=0.36, rho=-0.079  #P=0.93, rho=-0.009
cor.test(Vei_clos[PBCbase_cluster2,]$`k__Bacteria;p__Firmicutes;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella`,
         Vei_clos[PBCbase_cluster2,]$`k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales`,method = "spearman") #P=0.36, rho=-0.079  #P=0.89, rho=-0.017
cor.test(Vei_clos[PBCbase_resp_6m,]$`k__Bacteria;p__Firmicutes;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella`,
         Vei_clos[PBCbase_resp_6m,]$`k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales`,method = "spearman") #P=0.36, rho=-0.079  #P=0.89, rho=-0.017
cor.test(Vei_clos[PBCbase_noresp_6m,]$`k__Bacteria;p__Firmicutes;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella`,
         Vei_clos[PBCbase_noresp_6m,]$`k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales`,method = "spearman") #P=0.36, rho=-0.079  #P=0.89, rho=-0.017


###########responder and non-responder 结果可视化*figure##############
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin")
maaslin_species_PBCresp=read.csv("maaslin_species PBCresp.csv",header = T,row.names=1)
row.names(maaslin_species_PBCresp)=maaslin_species_PBCresp$species
library(dplyr)
maaslin_species_PBCresp_dif=filter(maaslin_species_PBCresp,pval<0.01) %>% row.names(.)
maaslin_species_PBCresp_dif[1:20]
maaslin_species_PBCresp_dif_20plot=maaslin_species_PBCresp[1:20,]
species_resp_dif_20=maaslin_species_PBCresp_dif[1:20]
newname=vector()
name=strsplit(species_resp_dif_20,"s__")
for (i in 1:20) {
  newname[i]=name[[i]][2]
}
newname[8]="Veillonella spp"
newname[15]="Anaerobutyricum spp"
row.names(maaslin_species_PBCresp_dif_20plot)=newname
maaslin_species_PBCresp_dif_20plot$species=newname
maaslin_species_PBCresp_dif_20plot <- maaslin_species_PBCresp_dif_20plot %>%  mutate(species = fct_reorder(species,coef))

#展现PBC resp改变最大的20taxa coefficient
library(ggplot2)
ggplot(data = maaslin_species_PBCresp_dif_20plot,aes(x=species,y=coef,fill=value))+
  geom_bar(stat = "identity",position = "identity",width = 0.6)+
  scale_fill_manual(values = c("#85ABD1"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  theme(legend.position =c(0.8,0.7))+
  theme(text = element_text(size = 10))+
  labs(fill="")+
  coord_flip()

#展现PBC resp改变最大且相对丰度>0.0001的taxa
colnames(species_mean_PBCbase)="mean_abundance"
species_mean_PBCbase_0.0001=filter(species_mean_PBCbase,mean_abundance>0.0001) %>% row.names(.)
index=maaslin_species_PBCresp_dif %in% species_mean_PBCbase_0.0001 #为逻辑值向量
index=which(index==TRUE)#为位置索引
maaslin_species_PBCresp_dif_0.0001_top20=maaslin_species_PBCresp_dif[index]
maaslin_species_PBCresp_dif_0.0001_top20=maaslin_species_PBCresp_dif_0.0001_top20[-21]

maaslin_species_PBCresp_dif_0.0001_20plot=maaslin_species_PBCresp[maaslin_species_PBCresp_dif_0.0001_top20,]

newname=vector()
name=strsplit(maaslin_species_PBCresp_dif_0.0001_top20,"s__")
for (i in 1:20) {
  newname[i]=name[[i]][2]
}
newname[4]="Bilophila spp"
newname[8]="Phascolarctobacterium spp"
newname[9]="Unidentified"
newname[15]="Sutterellaceae spp"
row.names(maaslin_species_PBCresp_dif_0.0001_20plot)=newname
maaslin_species_PBCresp_dif_0.0001_20plot$species=newname
library(forcats)
maaslin_species_PBCresp_dif_0.0001_20plot <- maaslin_species_PBCresp_dif_0.0001_20plot %>%  mutate(species = fct_reorder(species,coef))

#展现PBC resp改变最大的20taxa coefficient
library(ggplot2)
ggplot(data = maaslin_species_PBCresp_dif_0.0001_20plot,aes(x=species,y=coef,fill=value))+
  geom_bar(stat = "identity",position = "identity",width = 0.6)+
  scale_fill_manual(values = c("#85ABD1","#BFBFBF"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "null")+
  theme(legend.position ="null")+
  theme(text = element_text(size = 10))+
  labs(fill="")+
  coord_flip()