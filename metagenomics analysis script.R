####main codes for PBC multi-omics study
###for questions, please contact QiaoyanLiu (lqy2822@163.com)
##################################################  load data ####################################################
#metadata
setwd("D:/PBC multiomics data/metagenomics/data")
metadata=read.csv("PBC metagenomics sample metadata.csv")

#derive the sampleID in each subgroup
library(dplyr)
HC=row.names(filter(metadata,group=="HC"))

#remove PBC patients with antibiotics usage
PBCbase=row.names(filter(metadata,group=="PBCbase"& Antibiotics=="0"))

#baseline samples of UDCA 6month responders and non-responder 
PBCbase_resp_6m=row.names(filter(metadata,group=="PBCbase"& Antibiotics=="0"& Response_6m=="1"))
PBCbase_noresp_6m=row.names(filter(metadata,group=="PBCbase"& Antibiotics=="0"& Response_6m=="0"))

PBCbase_cluster1=filter(cluster.res,cluster=="clt1") %>% row.names()
PBCbase_cluster2=filter(cluster.res,cluster=="clt2") %>% row.names()

#response and non-response sample ID in each cluster
PBCbase_cluster1_resp=intersect(PBCbase_cluster1,PBCbase_resp_6m)
PBCbase_cluster1_noresp=intersect(PBCbase_cluster1,PBCbase_noresp_6m)
PBCbase_cluster2_resp=intersect(PBCbase_cluster2,PBCbase_resp_6m)
PBCbase_cluster2_noresp=intersect(PBCbase_cluster2,PBCbase_noresp_6m)


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


###################class data and plot##################
setwd("D:/multi-omics data/metagenomics data/re-assembled data/data/class")
class_abs=read.csv("class absolute abundance renamed.csv",header = T,row.names=1)
class_rela=read.csv("class relative abundance renamed.csv",header = T,row.names = 1)

#clr transformation
library(compositions)
library(zCompositions)
#all taxa
class_abs_clr=cmultRepl(class_abs, label=0, method="CZM")  %>% clr(.)
class_abs_clr=t(class_abs_clr) %>% as.data.frame()
class_abs_clr=class_abs_clr[row.names(metadata),]
class_abs_clr_scale=scale(class_abs_clr) %>% as.data.frame()
#keep class with a relative abundance 0.001% in at least 1% samples
class_rela_sample=class_rela[,all_sample]
library(plyr)
class_rela_sample=class_rela_sample %>% mutate_each(funs(. / sum(.)))
colSums(class_rela_sample)
class_rela_filter=class_rela_sample
class_rela_filter[class_rela_filter<0.00001]=0
class_rela_filter[class_rela_filter>0.00001]=1
class_rela_filter<- class_rela_sample[which(rowSums(class_rela_filter) >=4),]

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

#integrate metadata cluster information
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

#calculate cluster assocaition with phenotype
#categorial factor
chisq.test(metadata[PBCbase,]$Response_6m,cluster.res[PBCbase,1])$observed 
chisq.test(metadata[PBCbase,]$gp210,cluster.res[,1])$observed
chisq.test(metadata[PBCbase,]$gp210,metadata[PBCbase,]$Response_6m)$observed

ggplot(metadata_PBCbase,aes(x=cluster,y=TB,color=cluster))+
  geom_boxplot()+
  geom_jitter(width = 0.2,size=0.9)+
  scale_color_manual(values = c("#088288","#F19433"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "null")


################plot showing cluster and clinical relevance#####################

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




#################alpha diversity############
library(vegan)
shannon=diversity(t(species_rela_filter)) %>% as.data.frame()
colnames(shannon)="shannon"


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

#group plot
pcoa_plot_HC_PBC_Ait<-ggplot(data_HC_PBC_Ait,aes(x=V1,y=V2,colour=group)) + 
  geom_point(size=1.5) +
  theme_bw() +  
  scale_color_manual(values=c("#3F5BA7","#D72326"))+theme_bw()+ 
  guides(color=guide_legend(title = NULL)) + 
  stat_ellipse(level = 0.8)+  
  theme(axis.title.x = element_text(size=15,family="sans"),
        axis.title.y = element_text(size=15,family="sans",angle=90), 
        axis.text.y=element_text(size=12,family="sans"), 
        axis.text.x=element_text(size=12,family="sans"), 
        panel.grid = element_blank()
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



#3.PBCbase according to cluster **figure
#Plot the PCoA based on Aitchison distance among PBCbase and HC
pcoa_PBC_Ait<- cmdscale(Ait_genus_abs[PBCbase,PBCbase], k = nrow(metadata[PBCbase,]) - 1, eig = TRUE, add = TRUE)

#Extract PCoA coordinates
pc12_PBC_Ait<- pcoa_PBC_Ait$points[,1:2]
pc12_PBC_Ait=as.data.frame(pc12_PBC_Ait)
pc_importance_PBC_Ait<-round(pcoa_PBC_Ait$eig/sum(pcoa_PBC_Ait$eig)*100,digits = 2)
pc12_PBC_Ait$SampleID<-row.names(pc12_PBC_Ait)
data_PBC_Ait<-merge(pc12_PBC_Ait,metadata_PBCHC[PBCbase,],by="SampleID")



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
  

################PCo1-taxa association *figure#################

#species
PBC_PCo1_species=cbind(data_PBC_Ait[PBCbase,c(1,2,3)],species_abs_clr[PBCbase,row.names(species_rela_filter)])

corR_PCo1_species=sapply(2, function(x) sapply(4:3857,function(y){cor.test(PBC_PCo1_species[,x],PBC_PCo1_species[,y],method = "spearman")$estimate}))
corP_PCo1_species=sapply(2, function(x) sapply(4:3857,function(y){cor.test(PBC_PCo1_species[,x],PBC_PCo1_species[,y],method = "spearman")$p.value}))
corq_PCo1_species=p.adjust(corP_PCo1_species,"fdr") %>% matrix(.,ncol = 1)
row.names(corR_PCo1_species)=colnames(PBC_PCo1_species)[4:3857]
row.names(corP_PCo1_species)=colnames(PBC_PCo1_species)[4:3857]
row.names(corq_PCo1_species)=colnames(PBC_PCo1_species)[4:3857]

species_mean_PBCbase=apply(species_rela_filter[,PBCbase],1,mean) %>% as.data.frame(.)
row.names(species_mean_PBCbase)=row.names(species_rela_filter)


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

############plot:distribution of taxa abundance along PCo1 #############
data_PBC_Ait_reorder <- data_PBC_Ait[order(data_PBC_Ait$V1),]
data_PBC_Ait_reorder$SampleID=factor(data_PBC_Ait_reorder$SampleID,levels= unique(data_PBC_Ait_reorder$SampleID))

ggplot(data= data_PBC_Ait_reorder,aes(x=SampleID,y=Ruminococcus_563))+geom_bar(mapping = NULL, data = NULL, stat = "identity",
                                                                                 width=0.9, position="dodge",fill="light blue")
ggplot(data = data_PBC_Ait_reorder,aes(x=cluster,y=Bacteroides,fill=cluster))+geom_boxplot()+
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


 
  


  
#############MaAsin result#################
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











