
############## read in data #####################
#metadata
setwd("D:/PBC multiomics data/metabolomics/serum/data")
metadata=read.csv("serum_metadata.csv")
row.names(metadata)=metadata$X
metadata=metadata[,-1]
metadata$group=as.factor(metadata$group)
metadata$Stage=as.factor(metadata$Stage)
metadata$Antibiotics=as.numeric(metadata$Antibiotics)
metadata$subjectID=as.factor(metadata$subjectID)
metadata$gender=as.factor(metadata$gender)
metadata$Response_6m=as.factor(metadata$Response_6m)
metadata[,c(14:25)]=sapply(metadata[,c(14:25)],as.numeric)

#serum metabolome data
serum_MS=read.csv("serum_identified_all.csv")
row.names(serum_MS)=serum_MS[,1]
colnames(serum_MS)=serum_MS[1,]
serum_MS=serum_MS[-1,-1]
metaboname=colnames(serum_MS)
name=row.names(serum_MS)
serum_MS=as.data.frame(lapply(serum_MS,as.numeric))
row.names(serum_MS)=name #"as.numeric" transformation distorts row and colnames
colnames(serum_MS)=metaboname
serum=serum_MS[-1,]
#extract PBC and HC samples 
serum_MS=rbind(serum_MS[1,],serum_MS[row.names(metadata),])
serum=serum[row.names(metadata),]
serum_log2=log2(serum)
#output data in simca-P analysis
write.csv(serum_log2,"serum_PBC.HC_log2.csv")

#load PCA loading points file
setwd("D:/PBC multiomics data/metabolomics/serum/data")
PCA_point=read.table("simca-P PBC_HC.txt",header = T,row.names = 1,sep = "\t")
metadata=cbind(metadata,PCA_point)
colnames(metadata)[33]="PCA.point"
metadata$PCA.point=as.numeric(metadata$PCA.point)

#select metabolites for differential analysis with MS score>0.7
serum_0.7=serum_MS %>% t()%>% as.data.frame()
serum_0.7=filter(serum_0.7,serum_0.7$`MS2 score`>0.7) %>% t()%>% as.data.frame()
serum_0.7=serum_0.7[-1,]
#data transformation
serum_0.7_log2=serum_log2[,colnames(serum_0.7)]
serum_0.7_log2_UV=as.data.frame(scale(serum_0.7_log2))
#read in metabolite class data
serum_class=read.csv("serum_identified_all_v2_class.csv",header = T,row.names=1)
serum_class=serum_class[,1:6]
serum_class$class=as.factor(serum_class$class)
serum_class$SuperClass=as.factor(serum_class$SuperClass)
#simplified class
serum_class_simpli=read.csv("serum_class.csv",header = T,row.names = 1)
serum_class_simpli$class=as.factor(serum_class_simpli$class)
serum_class_simpli$class=serum_class_simpli$class[,drop=TRUE]
row.names(serum_class_simpli)
summary(serum_class_simpli$class)


#library packages
library(ggplot2)
library(dplyr)
library(reshape2)
library(reshape)
library(vegan)
library(tidyr)

#derive the sampleID in each subgroup
library(dplyr)
HC=metadata %>% filter(group=="HC"& Antibiotics=="0") %>% filter(between(PCA.point,-21.6832,21.6832)) %>% row.names()
HC_no=intersect(c("HC_H063","HC_H1022","HC_H1076","HC_H1136","HC_H444","HC_H448","HC_H459",
        "HC_H496","HC_H497","HC_H532","HC_H565","HC_H779","HC_H827","HC_H840","HC_H966","HC_H976",
        "HC_H455","HC_H605","HC_H572"),HC)
index=HC %in% HC_no #为逻辑值向量
index=which(index==TRUE)#为位置索引
HC_new=HC[-index]

#remove classical PBCbase patients with antibiotics usage
#select samples within 95% CI in PCA
PBCbase=metadata %>% filter(group=="PBCbase"& Antibiotics=="0") %>% filter(between(PCA.point,-21.6832,21.6832)) %>% row.names()
scPBCbase=metadata %>% filter(group=="scPBCbase") %>% filter(between(PCA.point,-21.6832,21.6832)) %>% row.names()#several scPBC patients took antibiotics
PBC_6m=metadata %>% filter(group=="PBC_6m") %>% filter(between(PCA.point,-21.6832,21.6832)) %>% row.names()
PBCpair_6m=intersect(sub("T","",x=PBC_6m),PBCbase)#select baseline paired samples with both baseline and 6m
PBC_6m_pair=sub("PBC_","PBC_T",PBCpair_6m)#6m samples with both timepoints
PBCbase_resp=metadata %>% filter(Response_6m=="1") %>% row.names(.) %>% intersect(.,PBCbase)
PBCbase_noresp=metadata %>% filter(Response_6m=="0") %>% row.names(.) %>% intersect(.,PBCbase)
PBC_6m_pair_resp=substring(PBC_6m_pair,6,8) %in% substring(PBCbase_resp,5,7) %>% which(.==TRUE) %>% PBC_6m_pair[.]
PBC_6m_pair_noresp=substring(PBC_6m_pair,6,8) %in% substring(PBCbase_noresp,5,7) %>% which(.==TRUE) %>% PBC_6m_pair[.]
PBCpair_6m_resp=intersect(PBCpair_6m,PBCbase_resp)
PBCpair_6m_noresp=intersect(PBCpair_6m,PBCbase_noresp)

#read in cluster file
setwd("D:/PBC multiomics data/analysis V6/metagenomics/data")
cluster.res_PBCHC=read.csv("cluster.res_PBCHC.csv",header = T,row.names = 1)
#生成不同cluster的sampleID
PBCbase_fecal_cluster1=filter(cluster.res_PBCHC,cluster=="clt1") %>% row.names(.) %>% intersect(.,PBCbase_fecal)
PBCbase_fecal_cluster2=filter(cluster.res_PBCHC,cluster=="clt2") %>% row.names(.) %>% intersect(.,PBCbase_fecal)
PBCbase_cluster1= metadata[PBCbase,"fecal.sample"] %in% PBCbase_fecal_cluster1 %>% which(.==TRUE) %>% row.names(metadata[PBCbase,])[.]
PBCbase_cluster2= metadata[PBCbase,"fecal.sample"] %in% PBCbase_fecal_cluster2 %>% which(.==TRUE) %>% row.names(metadata[PBCbase,])[.]
PBCpair_6m_cluster1=intersect(PBCpair_6m,PBCbase_cluster1)
PBCpair_6m_cluster2=intersect(PBCpair_6m,PBCbase_cluster2)
PBC_6m_pair_cluster1=substring(PBC_6m_pair,6,8) %in% substring(PBCbase_cluster1,5,7) %>% which(.==TRUE) %>% PBC_6m_pair[.]
PBC_6m_pair_cluster2=substring(PBC_6m_pair,6,8) %in% substring(PBCbase_cluster2,5,7) %>% which(.==TRUE) %>% PBC_6m_pair[.]

#generate metadata with cluster information  **HC_new
PBCHC_cluster=c(PBCbase_cluster1,PBCbase_cluster2,HC_new) %>% as.data.frame()
PBCHC_cluster$cluster=c(rep("clt1",65),rep("clt2",46),rep("HC",111))
PBCHC_cluster$cluster=as.factor(PBCHC_cluster$cluster) 
row.names(PBCHC_cluster)=PBCHC_cluster$.
colnames(PBCHC_cluster)[1]="sampleID"
PBCHC_cluster$cluster=as.factor(PBCHC_cluster$cluster)
colnames(metadata)
metadata_PBCHC=metadata[c(PBCbase,HC_new),]
metadata_PBCHC$sampleID=row.names(metadata_PBCHC)
metadata_PBCHC=cbind(metadata_PBCHC[,3:6],PBCHC_cluster[row.names(metadata_PBCHC),])
metadata_PBCHC[,c("sampleID","cluster")]
filter(metadata_PBCHC,cluster=="clt1") %>% row.names(.)
metadata_PBCHC$group=as.factor(metadata_PBCHC$group)
metadata_PBCHC$group=metadata_PBCHC$group[,drop=TRUE]
metadata_PBCHC$gender=as.factor(metadata_PBCHC$gender)
metadata_PBCHC$age=as.numeric(metadata_PBCHC$age)
metadata_PBCHC$BMI=as.numeric(metadata_PBCHC$BMI)
metadata_PBCHC$cluster=as.factor(metadata_PBCHC$cluster)
str(metadata_PBCHC)

#############linear model HC_new################
lm_PBCHC=as.data.frame(matrix(0,ncol=1,nrow=281))
#PBCHC
for (i in 1:ncol(serum_0.7_log2_UV)) {
  lm_PBCHC.fit=lm(serum_0.7_log2_UV[c(PBCbase,HC_new),i]~group+gender+age+BMI,data = metadata_PBCHC)    
  lm_PBCHC$pval_PBCHC[i] <- summary(lm_PBCHC.fit)$coefficients["groupPBCbase","Pr(>|t|)"] 
  lm_PBCHC$tval_PBCHC[i] <- summary(lm_PBCHC.fit)$coefficients["groupPBCbase","t value"]  
  }
lm_PBCHC$qval_PBCHC=p.adjust(lm_PBCHC$pval_PBCHC,"fdr")
lm_PBCHC$PBCbase.abundance=apply(serum_0.7[PBCbase,],2,mean)
lm_PBCHC$HC.abundance=apply(serum_0.7[HC_new,],2,mean)
#responder and non-responder
for (i in 1:ncol(serum_0.7_log2_UV)) {
  lm_PBCHC.fit=lm(serum_0.7_log2_UV[PBCbase,i]~Response_6m+gender+age+BMI,data = metadata[PBCbase,])    
  lm_PBCHC$pval_PBCresp.noresp[i] <- summary(lm_PBCHC.fit)$coefficients["Response_6m1","Pr(>|t|)"] 
  lm_PBCHC$tval_PBCresp.noresp[i] <- summary(lm_PBCHC.fit)$coefficients["Response_6m1","t value"]  
}
lm_PBCHC$qval_PBCresp.noresp=p.adjust(lm_PBCHC$pval_PBCresp.noresp,"fdr")
lm_PBCHC$PBCresp.abundance=apply(serum_0.7[PBCbase_resp,],2,mean)
lm_PBCHC$PBCnoresp.abundance=apply(serum_0.7[PBCbase_noresp,],2,mean)
#scPBCbase and HC
for (i in 1:ncol(serum_0.7_log2_UV)) {
  lm_PBCHC.fit=lm(serum_0.7_log2_UV[c(scPBCbase,HC_new),i]~group+gender+age+BMI,data = metadata[c(scPBCbase,HC_new),])    
  lm_PBCHC$pval_scPBCHC[i] <- summary(lm_PBCHC.fit)$coefficients["groupscPBCbase","Pr(>|t|)"] 
  lm_PBCHC$tval_scPBCHC[i] <- summary(lm_PBCHC.fit)$coefficients["groupscPBCbase","t value"]  
}
lm_PBCHC$qval_scPBCHC=p.adjust(lm_PBCHC$pval_scPBCHC,"fdr")
lm_PBCHC$scPBC.abundance=apply(serum_0.7[scPBCbase,],2,mean)
lm_PBCHC$Hc.abundance=apply(serum_0.7[HC_new,],2,mean)

#clt1HC
for (i in 1:ncol(serum_0.7_log2_UV)) {
   lm_PBCHC.fit=lm(serum_0.7_log2_UV[c(PBCbase_cluster1,HC_new),i]~cluster+gender+age+BMI,data = metadata_PBCHC[c(PBCbase_cluster1,HC_new),])    
   lm_PBCHC$pval_clt1HC[i] <- summary(lm_PBCHC.fit)$coefficients["clusterHC","Pr(>|t|)"]  
   lm_PBCHC$tval_clt1HC[i] <- summary(lm_PBCHC.fit)$coefficients["clusterHC","t value"] 
   }
lm_PBCHC$qval_clt1HC=p.adjust(lm_PBCHC$pval_clt1HC,"fdr")
lm_PBCHC$clt1.abundance=apply(serum_0.7[PBCbase_cluster1,],2,mean)
#clt2HC
for (i in 1:ncol(serum_0.7_log2_UV)) {
     lm_PBCHC.fit=lm(serum_0.7_log2_UV[c(PBCbase_cluster2,HC_new),i]~cluster+gender+age+BMI,data = metadata_PBCHC[c(PBCbase_cluster2,HC_new),])    
     lm_PBCHC$pval_clt2HC[i] <- summary(lm_PBCHC.fit)$coefficients["clusterHC","Pr(>|t|)"]  
     lm_PBCHC$tval_clt2HC[i] <- summary(lm_PBCHC.fit)$coefficients["clusterHC","t value"] 
   }
lm_PBCHC$qval_clt2HC=p.adjust(lm_PBCHC$pval_clt2HC,"fdr")
lm_PBCHC$clt2.abundance=apply(serum_0.7[PBCbase_cluster2,],2,mean)
#clt1clt2
for (i in 1:ncol(serum_0.7_log2_UV)) {
  lm_PBCHC.fit=lm(serum_0.7_log2_UV[PBCbase,i]~cluster+gender+age+BMI,data = metadata_PBCHC[PBCbase,])    
  lm_PBCHC$pval_clt1clt2[i] <- summary(lm_PBCHC.fit)$coefficients["clusterclt2","Pr(>|t|)"]  
  lm_PBCHC$tval_clt1clt2[i] <- summary(lm_PBCHC.fit)$coefficients["clusterclt2","t value"] 
  }
lm_PBCHC$qval_clt1clt2=p.adjust(lm_PBCHC$pval_clt1clt2,"fdr")
row.names(lm_PBCHC)=colnames(serum_0.7_log2_UV)

write.csv(lm_PBCHC,"lm_PBCHC.new.csv")

###############比较不同stage的代谢物####################
PBCbase_early=metadata %>% filter(Stage=="early") %>% row.names(.) %>% intersect(.,PBCbase)
PBCbase_advanced=metadata %>% filter(Stage=="moderately advanced"|Stage=="advanced") %>% row.names(.) %>% intersect(.,PBCbase)

wilcox.test(serum_0.7[PBCbase_early,"Indole-3-propionic acid"],serum_0.7[PBCbase_advanced,"Indole-3-propionic acid"])
colnames(serum_0.7)


#############paired 检验################
#1.确定cluster1/2用于前后比较的feature:与HC及各自对比有差异的feature
#1.1 clt1
library(dplyr)
# HC_new
serum_clt1HC_dif=lm_PBCHC %>% filter(.,qval_clt1HC<0.05) %>% row.names(.)
serum_clt1clt2_dif=lm_PBCHC %>% filter(.,qval_clt1clt2<0.05) %>% row.names(.)

serum_clt1_pair=c(serum_clt1HC_dif,serum_clt1clt2_dif)
index=duplicated(serum_clt1_pair)
index=which(index==FALSE)
serum_clt1_pair=serum_clt1_pair[index]
#1.2 clt2 
#HC_new
serum_clt2HC_dif=lm_PBCHC %>% filter(.,qval_clt2HC<0.05) %>% row.names(.)
serum_clt1clt2_dif=lm_PBCHC %>% filter(.,qval_clt1clt2<0.05) %>% row.names(.)

serum_clt2_pair=c(serum_clt2HC_dif,serum_clt1clt2_dif)
index=duplicated(serum_clt2_pair)
index=which(index==FALSE)
serum_clt2_pair=serum_clt2_pair[index]

#paired wilcox test 检验
#1.1 PBCHC  #select metabolites
res_serum_pair.wilcox=as.data.frame(matrix(0,ncol=1,nrow=173))
for (i in 1:173) {
  res_serum_pair.wilcox$PBCbase_6m.pval[i]=wilcox.test(serum_0.7[PBCpair_6m,serum_PBCHC_dif_lm][,i],
                                                       serum_0.7[PBC_6m_pair,serum_PBCHC_dif_lm][,i],paired = T)$p.value
}
res_serum_pair.wilcox$PBCbase_6m.qval=p.adjust(res_serum_pair.wilcox$PBCbase_6m.pval,"fdr")
res_serum_pair.wilcox$PBCbase.abundance=apply(serum_0.7[PBCpair_6m,serum_PBCHC_dif_lm],2,mean)
res_serum_pair.wilcox$PBC.6m.abundance=apply(serum_0.7[PBC_6m_pair,serum_PBCHC_dif_lm],2,mean)
res_serum_pair.wilcox$change=res_serum_pair.wilcox$PBC.6m.abundance-res_serum_pair.wilcox$PBCbase.abundance
res_serum_pair.wilcox$PBCHC.tval=lm_PBCHC[serum_PBCHC_dif_lm,"tval_PBCHC"]
row.names(res_serum_pair.wilcox)=serum_PBCHC_dif_lm

#1.2 PBCHC  #all metabolites
res_serum_pair.wilcox=as.data.frame(matrix(0,ncol=1,nrow=281))
for (i in 1:281) {
  res_serum_pair.wilcox$PBCbase_6m.pval[i]=wilcox.test(serum_0.7[PBCpair_6m,i],
                                                       serum_0.7[PBC_6m_pair,i],paired = T)$p.value
}
res_serum_pair.wilcox$PBCbase_6m.qval=p.adjust(res_serum_pair.wilcox$PBCbase_6m.pval,"fdr")
res_serum_pair.wilcox$PBCbase.abundance=apply(serum_0.7[PBCpair_6m,],2,mean)
res_serum_pair.wilcox$PBC.6m.abundance=apply(serum_0.7[PBC_6m_pair,],2,mean)
res_serum_pair.wilcox$change=res_serum_pair.wilcox$PBC.6m.abundance-res_serum_pair.wilcox$PBCbase.abundance
res_serum_pair.wilcox$PBCHC.tval=lm_PBCHC[,"tval_PBCHC"]
res_serum_pair.wilcox$PBCHC.pval=lm_PBCHC[,"pval_PBCHC"]
res_serum_pair.wilcox$PBCHC.qval=lm_PBCHC[,"qval_PBCHC"]

row.names(res_serum_pair.wilcox)=colnames(serum_0.7)


#2.1  clt1  #select metabolites
res_serum_pair.wilcox_clt1=as.data.frame(matrix(0,ncol=1,nrow=165))
for (i in 1:165) {
  res_serum_pair.wilcox_clt1$cl1base_6m.pval[i]=wilcox.test(serum_0.7[PBCpair_6m_cluster1,serum_clt1_pair][,i],
                                                           serum_0.7[PBC_6m_pair_cluster1,serum_clt1_pair][,i],paired = T)$p.value
}
res_serum_pair.wilcox_clt1$clt1base_6m.qval=p.adjust(res_serum_pair.wilcox_clt1$clt1base_6m.pval,"fdr")
res_serum_pair.wilcox_clt1$clt1base.abundance=apply(serum_0.7[PBCpair_6m_cluster1,serum_clt1_pair],2,mean)
res_serum_pair.wilcox_clt1$clt1.6m.abundance=apply(serum_0.7[PBC_6m_pair_cluster1,serum_clt1_pair],2,mean)
res_serum_pair.wilcox_clt1$change=res_serum_pair.wilcox_clt1$clt1.6m.abundance-res_serum_pair.wilcox_clt1$clt1base.abundance
res_serum_pair.wilcox_clt1$clt1HC.qval=lm_PBCHC[serum_clt1_pair,"qval_clt1HC"]
res_serum_pair.wilcox_clt1$clt1HC.tval=lm_PBCHC[serum_clt1_pair,"tval_clt1HC"]
res_serum_pair.wilcox_clt1$clt1clt2.qval=lm_PBCHC[serum_clt1_pair,"qval_clt1clt2"]
res_serum_pair.wilcox_clt1$clt1clt2.tval=lm_PBCHC[serum_clt1_pair,"tval_clt1clt2"]
row.names(res_serum_pair.wilcox_clt1)=serum_clt1_pair

#2.2  clt2  #all metabolites
res_serum_pair.wilcox_clt1=as.data.frame(matrix(0,ncol=1,nrow=281))
for (i in 1:281) {
  res_serum_pair.wilcox_clt1$clt1base_6m.pval[i]=wilcox.test(serum_0.7[PBCpair_6m_cluster1,i],
                                                            serum_0.7[PBC_6m_pair_cluster1,i],paired = T)$p.value
}
res_serum_pair.wilcox_clt1$clt1base_6m.qval=p.adjust(res_serum_pair.wilcox_clt1$clt1base_6m.pval,"fdr")
res_serum_pair.wilcox_clt1$clt1base.abundance=apply(serum_0.7[PBCpair_6m_cluster1,],2,mean)
res_serum_pair.wilcox_clt1$clt1.6m.abundance=apply(serum_0.7[PBC_6m_pair_cluster1,],2,mean)
res_serum_pair.wilcox_clt1$change=res_serum_pair.wilcox_clt1$clt1.6m.abundance-res_serum_pair.wilcox_clt1$clt1base.abundance
res_serum_pair.wilcox_clt1$clt1HC.qval=lm_PBCHC[,"qval_clt1HC"]
res_serum_pair.wilcox_clt1$clt1HC.tval=-lm_PBCHC[,"tval_clt1HC"]
res_serum_pair.wilcox_clt1$clt1clt2.qval=lm_PBCHC[,"qval_clt1clt2"]
res_serum_pair.wilcox_clt1$clt1clt2.tval=-lm_PBCHC[,"tval_clt1clt2"]
row.names(res_serum_pair.wilcox_clt1)=colnames(serum_0.7)

#3.1 clt2  #select metabolites
res_serum_pair.wilcox_clt2=as.data.frame(matrix(0,ncol=1,nrow=146))
for (i in 1:146) {
  res_serum_pair.wilcox_clt2$clt2base_6m.pval[i]=wilcox.test(serum_0.7[PBCpair_6m_cluster2,serum_clt2_pair][,i],
                                                            serum_0.7[PBC_6m_pair_cluster2,serum_clt2_pair][,i],paired = T)$p.value
}
res_serum_pair.wilcox_clt2$clt2base_6m.qval=p.adjust(res_serum_pair.wilcox_clt2$clt2base_6m.pval,"fdr")
res_serum_pair.wilcox_clt2$clt2base.abundance=apply(serum_0.7[PBCpair_6m_cluster2,serum_clt2_pair],2,mean)
res_serum_pair.wilcox_clt2$clt2.6m.abundance=apply(serum_0.7[PBC_6m_pair_cluster2,serum_clt2_pair],2,mean)
res_serum_pair.wilcox_clt2$change=res_serum_pair.wilcox_clt2$clt2.6m.abundance-res_serum_pair.wilcox_clt2$clt2base.abundance
res_serum_pair.wilcox_clt2$clt2HC.qval=lm_PBCHC[serum_clt2_pair,"qval_clt2HC"]
res_serum_pair.wilcox_clt2$clt2HC.tval=lm_PBCHC[serum_clt2_pair,"tval_clt2HC"]
res_serum_pair.wilcox_clt2$clt1clt2.qval=lm_PBCHC[serum_clt2_pair,"qval_clt1clt2"]
res_serum_pair.wilcox_clt2$clt1clt2.tval=lm_PBCHC[serum_clt2_pair,"tval_clt1clt2"]
row.names(res_serum_pair.wilcox_clt2)=serum_clt2_pair

#3.2 clt2  #all metabolites
res_serum_pair.wilcox_clt2=as.data.frame(matrix(0,ncol=1,nrow=281))
for (i in 1:281) {
  res_serum_pair.wilcox_clt2$clt2base_6m.pval[i]=wilcox.test(serum_0.7[PBCpair_6m_cluster2,i],
                                                             serum_0.7[PBC_6m_pair_cluster2,i],paired = T)$p.value
}
res_serum_pair.wilcox_clt2$clt2base_6m.qval=p.adjust(res_serum_pair.wilcox_clt2$clt2base_6m.pval,"fdr")
res_serum_pair.wilcox_clt2$clt2base.abundance=apply(serum_0.7[PBCpair_6m_cluster2,],2,mean)
res_serum_pair.wilcox_clt2$clt2.6m.abundance=apply(serum_0.7[PBC_6m_pair_cluster2,],2,mean)
res_serum_pair.wilcox_clt2$change=res_serum_pair.wilcox_clt2$clt2.6m.abundance-res_serum_pair.wilcox_clt2$clt2base.abundance
res_serum_pair.wilcox_clt2$clt2HC.qval=lm_PBCHC[,"qval_clt2HC"]
res_serum_pair.wilcox_clt2$clt2HC.tval=-lm_PBCHC[,"tval_clt2HC"]
res_serum_pair.wilcox_clt2$clt1clt2.qval=lm_PBCHC[,"qval_clt1clt2"]
res_serum_pair.wilcox_clt2$clt1clt2.tval=lm_PBCHC[,"tval_clt1clt2"]
row.names(res_serum_pair.wilcox_clt2)=colnames(serum_0.7)

setwd("D:/PBC multiomics data/analysis V6/metabolomics/serum metabolome/result/paired_new")
write.csv(res_serum_pair.wilcox,"res_serum_pair.wilcox_PBC.csv")
write.csv(res_serum_pair.wilcox_clt1,"res_serum_pair.wilcox_clt1.csv")
write.csv(res_serum_pair.wilcox_clt2,"res_serum_pair.wilcox_clt2.csv")


###########可视化展示:metabolome PCA analysis###################
library("FactoMineR")
library("factoextra")
res.pca <- PCA(serum_0.7_log2_UV[c(PBCbase,HC_new),], graph = FALSE)
eig.val <- get_eigenvalue(res.pca)
eig.val
#visualization
fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = metadata_PBCHC[c(PBCbase,HC_new),]$cluster, # color by groups
             palette = c("#088288","#F19433","#3F5BA7"),
             pointshape=19,
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",
             ggtheme = theme_minimal()
)  
#pdf 3*4
##################metabolome PCoA analysis(before and after) HC_new *figure#####################
library(vegan)
#1.cluster base and 6m
serum_dis=vegdist(serum_0.7_log2_UV[c(PBCpair_6m_cluster1,PBCpair_6m_cluster2,PBC_6m_pair_cluster1,PBC_6m_pair_cluster2,HC_new),],method = "euclidean")
metadata_centroid=c(PBCpair_6m_cluster1,PBCpair_6m_cluster2,PBC_6m_pair_cluster1,PBC_6m_pair_cluster2,HC_new) %>% as.data.frame()
metadata_centroid$group=c(rep("clt1.base",25),rep("clt2.base",15),rep("clt1.6m",25),rep("clt2.6m",15),rep("HC",111))
metadata_centroid$group=as.factor(metadata_centroid$group)
mod=betadisper(serum_dis,metadata_centroid$group)
df.pcoa = data.frame(x=mod$vectors[,1],y=mod$vectors[,2],group=mod$group)
centroids = data.frame(x=mod$centroids[,1],y=mod$centroids[,2],group = factor(rownames(mod$centroids)))
gg <- merge(df.pcoa,centroids, by ='group', suffixes = c('','.centroid')) 
gg$group  = factor(gg$group ,levels=c('clt1.base','clt1.6m','clt2.base','clt2.6m','HC'))
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
#2.all PBC base and 6m
metadata_centroid$group_all=c(rep("PBC.base",40),rep("PBC.6m",40),rep("HC",111))
metadata_centroid$group_all=as.factor(metadata_centroid$group_all)
mod=betadisper(serum_dis,metadata_centroid$group_all)
df.pcoa = data.frame(x=mod$vectors[,1],y=mod$vectors[,2],group=mod$group)
centroids = data.frame(x=mod$centroids[,1],y=mod$centroids[,2],group = factor(rownames(mod$centroids)))
gg <- merge(df.pcoa,centroids, by ='group', suffixes = c('','.centroid')) 
gg$group  = factor(gg$group ,levels=c('PBC.base','PBC.6m','HC'))
library(ggplot2)
ppcoaUDCA_ctrl <- ggplot(gg) +
  geom_point(aes(x=x,y=y,color=group),size=1,alpha=0.4) + scale_color_manual(values =c('#D72326','#F29E9C','#3F5BA7')) +
  geom_segment(aes(x=x.centroid, y =y.centroid, xend=x,yend=y,color=group),alpha=0.2) +
  geom_point(data=centroids,aes(x=x,y=y,color=group),size=3)+
  xlab(paste('PCoA1 (',100*round(mod$eig[1] / sum(mod$eig),2),'%)',sep='')) +
  ylab(paste('PCoA2 (',100*round(mod$eig[2] / sum(mod$eig),2),'%)',sep='')) + 
  theme_bw()+
  theme(legend.position = 'none', panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ppcoaUDCA_ctrl

library(vegan)
serum_dis=as.matrix(serum_dis)
adonis(serum_dis[c(PBCpair_6m,PBC_6m_pair),c(PBCpair_6m,PBC_6m_pair)]~group,data = metadata[c(PBCpair_6m,PBC_6m_pair),])

############比较cluster1/2治疗前后与HC的距离  *figure############  
#1.1 cluster1治疗前后与HC的距离
serum_dis_matrix=as.matrix(serum_dis) #将距离数据转化为矩阵
serum_clt1.base.HC=serum_dis_matrix[PBCpair_6m_cluster1,HC_new] %>% as.numeric(.) %>% as.data.frame(.)
serum_clt1.treat.HC=serum_dis_matrix[PBC_6m_pair_cluster1,HC_new] %>% as.numeric(.) %>% as.data.frame(.)
colnames(serum_clt1.base.HC)="distance"
colnames(serum_clt1.treat.HC)="distance"
wilcox.test(serum_clt1.base.HC$distance,serum_clt1.treat.HC$distance) #7.185e-16
mean(serum_clt1.base.HC$distance) #mean=21.60183
mean(serum_clt1.treat.HC$distance) #mean=22.38577
serum_clt1_HC=rbind(serum_clt1.base.HC,serum_clt1.treat.HC)
serum_clt1_HC$group=c(rep("clt1.base.HC",2775),rep("clt1.treat.HC",2775))
library(ggplot2)
ggplot(data = serum_clt1_HC,aes(x=group,y=distance,color=group))+
  geom_boxplot()+geom_jitter(width = 0.2,size=0.3)+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  scale_color_manual(values = c("#088288","#ADD5DB"))


#2.1 cluster2治疗前后与HC的距离
serum_dis_matrix=as.matrix(serum_dis) #将距离数据转化为矩阵
serum_clt2.base.HC=serum_dis_matrix[PBCpair_6m_cluster2,HC_new] %>% as.numeric(.) %>% as.data.frame(.)
serum_clt2.treat.HC=serum_dis_matrix[PBC_6m_pair_cluster2,HC_new] %>% as.numeric(.) %>% as.data.frame(.)
colnames(serum_clt2.base.HC)="distance"
colnames(serum_clt2.treat.HC)="distance"
wilcox.test(serum_clt2.base.HC$distance,serum_clt2.treat.HC$distance) #0.06845
mean(serum_clt2.base.HC$distance) #mean= 23.63023
mean(serum_clt2.treat.HC$distance) #mean=22.8327
serum_clt2_HC=rbind(serum_clt2.base.HC,serum_clt2.treat.HC)
serum_clt2_HC$group=c(rep("clt2.base.HC",1665),rep("clt2.treat.HC",1665))
library(ggplot2)
ggplot(data = serum_clt2_HC,aes(x=group,y=distance,color=group))+
  geom_boxplot()+geom_jitter(width = 0.2,size=0.3)+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  scale_color_manual(values = c("#F19433","#FFD966"))


#3.将各组治疗前后与HC的距离合并且绘图
serum_clt_HC=rbind(serum_clt1_HC,serum_clt2_HC)
ggplot(data = serum_clt_HC,aes(x=group,y=distance,color=group))+
  geom_boxplot()+geom_jitter(width = 0.2,size=0.8)+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  scale_color_manual(values = c("#088288","#ADD5DB","#F19433","#FFD966"))





############输出文件到metagenomics workspace##############
metadata$group
metadata$serum.sample
#将serum 代谢谱的样品名修改为粪便编号,代谢物加后缀"_serum"
metadata[row.names(serum_0.7),"fecal.sample"]
row.names(serum_0.7)
serum_0.7_for_multiomics=serum_0.7[c(PBCbase,PBC_6m,scPBCbase,HC),]
serum_0.7_log2_UV_for_multiomics=serum_0.7_log2_UV[c(PBCbase,PBC_6m,scPBCbase,HC),]
row.names(serum_0.7_for_multiomics)=metadata[c(PBCbase,PBC_6m,scPBCbase,HC),"fecal.sample"]
row.names(serum_0.7_log2_UV_for_multiomics)=metadata[c(PBCbase,PBC_6m,scPBCbase,HC),"fecal.sample"]
colnames(serum_0.7_for_multiomics)=paste(colnames(serum_0.7_for_multiomics),"serum",sep = "_")
colnames(serum_0.7_log2_UV_for_multiomics)=paste(colnames(serum_0.7_for_multiomics),"serum",sep = "_")
lm_PBCHC_for_multiomics=lm_PBCHC
row.names(lm_PBCHC_for_multiomics)=colnames(serum_0.7_for_multiomics)
#输出文件到
setwd("D:/PBC multiomics data/analysis V6/multiomics/data")
write.csv(serum_0.7_for_multiomics,"serum_0.7_for_multiomics.csv")
write.csv(serum_0.7_log2_UV_for_multiomics,"serum_0.7_log2_UV_for_multiomics.csv")
write.csv(lm_PBCHC_for_multiomics,"lm_PBCHC_serum_for_multiomics.csv")

############计算差异代谢物(lm PBCHC q<0.05)的临床相关性 #HC_new #############
library(dplyr)
serum_PBCHC_dif_lm=filter(lm_PBCHC,qval_PBCHC<0.05) %>% row.names(.)
clinic=c("IgM","IgA","IgG","ALT","AST","ALP","GGT","ALB","TB","CB")
merge_PBCbase=serum_0.7[PBCbase,serum_PBCHC_dif_lm]
merge_PBCbase=cbind(metadata[PBCbase,c("age","gender","BMI",clinic)],merge_PBCbase)

library(PResiduals)
corR_clini.serum=sapply(4:13, function(x) sapply(14:185, function(y) {
  partial_Spearman(merge_PBCbase[,x]|merge_PBCbase[,y]~merge_PBCbase[,"age"]+merge_PBCbase[,"gender"]+merge_PBCbase[,"BMI"],data=merge_PBCbase)$TS$TB$ts
}))
corp_clini.serum=sapply(4:13, function(x) sapply(14:185, function(y) {
  partial_Spearman(merge_PBCbase[,x]|merge_PBCbase[,y]~merge_PBCbase[,"age"]+merge_PBCbase[,"gender"]+merge_PBCbase[,"BMI"],data=merge_PBCbase)$TS$TB$pval
}))
corq_clini.serum=p.adjust(corp_clini.serum,"fdr") %>% matrix(.,ncol = 10)
row.names(corR_clini.serum)=serum_PBCHC_dif_lm
row.names(corp_clini.serum)=serum_PBCHC_dif_lm
row.names(corq_clini.serum)=serum_PBCHC_dif_lm
colnames(corR_clini.serum)=clinic
colnames(corp_clini.serum)=clinic
colnames(corq_clini.serum)=clinic
setwd("D:/PBC multiomics data/analysis V6/metabolomics/serum metabolome/result/HC_new")
write.csv(corp_clini.serum,"corP_clini.serum_new.csv")
write.csv(corq_clini.serum,"corq_clini.serum_new.csv")
write.csv(corR_clini.serum,"corR_clini.serum_new.csv")

#计数显著的相关性数目
index=vector()
for (i in 1:172) {
  medi=corq_clini.serum[i,]<0.05
  index[i]=which(medi==TRUE) %>% length(.)
  
}
sum(index)

############可视化展示-boxplot###########
library(dplyr)
serum_0.7_log2_plot=cbind(serum_0.7_log2,metadata$group)
colnames(serum_0.7_log2_plot)[282]="group"
serum_0.7_log2_plot$group=factor(serum_0.7_log2_plot$group,levels = c("PBCbase","PBC_6m","PBC_1y","HC","scPBCbase","scPBC_6m"))

serum_0.7_log2_plot_PBCHC=serum_0.7_log2_plot[c(HC_new,PBCbase),]
serum_0.7_log2_plot_PBCHC$group=serum_0.7_log2_plot_PBCHC$group[,drop=TRUE]

serum_0.7_log2_plot_PBCHC=cbind(serum_0.7_log2_plot_PBCHC,metadata_PBCHC[c(HC_new,PBCbase),"cluster"])
colnames(serum_0.7_log2_plot_PBCHC)[283]="cluster"
serum_0.7_log2_plot_PBCHC$cluster=factor(serum_0.7_log2_plot_PBCHC$cluster,levels = c("HC","clt1","clt2"))
serum_0.7_log2_plot_PBCHC$group=factor(serum_0.7_log2_plot_PBCHC$group,levels = c("HC","PBCbase"))


library(ggplot2)
#1.所有样本
ggplot(data = serum_0.7_log2_plot,aes(x=group,y=serum_0.7_log2_plot$`Phenyllactic acid`,color=group))+geom_boxplot()+
  geom_jitter(width = 0.2,size=0.9)+
  labs(y="Phenyllactic acid")+
  scale_color_manual(values = c("#F2644B","#87C3DF","#02425A","#89CCA4","#8167AD","#C8B8DA"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank())

#2.PBCbase 和HC_new样本
#2.1 按照cluster分组展示 
ggplot(data =serum_0.7_log2_plot_PBCHC,aes(x=cluster,y=serum_0.7_log2_plot_PBCHC$`Indoxyl sulfate`,color=cluster))+geom_boxplot()+
  geom_jitter(width = 0.2,size=0.9)+
  labs(y="serum indoxyl sulfate",x="")+
  scale_color_manual(values = c("#FFB733","#B75454","#83AAF0"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  annotate("text",x=1.5,y=2,label="lm FDR=0.0318")
 
#2.2 按照group展示 HC_new  **figure
library(ggplot2)
ggplot(data =serum_0.7_log2_plot_PBCHC,aes(x=group,y=serum_0.7_log2_plot_PBCHC$`Tauroursodeoxycholic acid`,fill=group))+geom_boxplot(color="grey")+

  labs(y="TUDCA",x="")+
  scale_fill_manual(values = c("#3F5BA7","#D72326"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")



#2.3 仅包括PBCbase样本，按照cluster展示
serum_0.7_log2_plot_PBCbase=serum_0.7_log2_plot_PBCHC[PBCbase,]
serum_0.7_log2_plot_PBCbase$cluster=serum_0.7_log2_plot_PBCbase$cluster[,drop=TRUE]
ggplot(data =serum_0.7_log2_plot_PBCbase,aes(x=cluster,y=serum_0.7_log2_plot_PBCbase$`Indoxyl sulfate`,fill=cluster))+geom_boxplot(color="grey")+
  
  labs(y="serum Indoxyl sulfate",x="")+
  scale_fill_manual(values = c("#088288","#F19433"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")

#2.4所有样本，按照cluser展示 **figure
ggplot(data =serum_0.7_log2_plot_PBCHC,aes(x=cluster,y=`Indoxyl sulfate`,fill=cluster))+geom_boxplot(color="grey")+
  
  labs(y="serum Indoxyl sulfate",x="")+
  scale_fill_manual(values = c("#3F5BA7","#088288","#F19433"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")


#治疗前后
#boxplot图添加前后的连线
library(ggpubr)
library(ggplot2)
#1.all PBCbase
data=data.frame(baseline=serum_0.7_log2[PBCpair_6m,"Cholesterol sulfate"],treat=serum_0.7_log2[PBC_6m_pair,"Cholesterol sulfate"]) #shannon.base|6m为向量
colnames(data)
ggpaired(data, cond1 = "baseline", cond2 = "treat",
         color = "condition", palette = "npg",ylab = "serum Cholesterol sulfate",line.color = "grey",
         title = "all PBCbase无差异")+
         theme(legend.position = "none")  
#2.cluster1
data=data.frame(baseline=serum_0.7_log2[PBCpair_6m_cluster1,"Dehydroepiandrosterone sulfate"],treat=serum_0.7_log2[PBC_6m_pair_cluster1,"Dehydroepiandrosterone sulfate"]) #shannon.base|6m为向量
colnames(data)
ggpaired(data, cond1 = "baseline", cond2 = "treat",
         color = "condition", palette = "npg",ylab = "serum DHEA-S",line.color = "grey",
         title = "cluster1 P=0.01\n FDR=0.19")+
  theme(legend.position = "none") 
#3.cluster2
data=data.frame(baseline=serum_0.7_log2[PBCpair_6m_cluster2,"Dehydroepiandrosterone sulfate"],treat=serum_0.7_log2[PBC_6m_pair_cluster2,"Dehydroepiandrosterone sulfate"]) #shannon.base|6m为向量
colnames(data)
ggpaired(data, cond1 = "baseline", cond2 = "treat",
         color = "condition", palette = "npg",ylab = "serum DHEA-S",line.color = "grey",
         title = "cluster2 无差异")+
  theme(legend.position = "none") 

#############可视化展示相关性热图 HC_new#############
library(dplyr)
corq_clini.serum_plot=corq_clini.serum
corq_clini.serum_plot=as.data.frame(corq_clini.serum_plot)
#生成corq<0.05的代谢物索引
index=vector()
for (i in 1:172) {
  index[i]=which(corq_clini.serum_plot[i,]<0.05) %>% length(.)
}
index=which(index>0)

corq_clini.serum_plot=corq_clini.serum_plot[index,]
corR_clini.serum_plot=corR_clini.serum[index,]

corq_clini.serum_plot[corq_clini.serum_plot<0.0001]="****"
corq_clini.serum_plot[corq_clini.serum_plot>=0.0001&corq_clini.serum_plot<0.001]="***"
corq_clini.serum_plot[corq_clini.serum_plot>=0.001&corq_clini.serum_plot<0.01]="**"
corq_clini.serum_plot[corq_clini.serum_plot>=0.01&corq_clini.serum_plot<0.05]="*"
corq_clini.serum_plot[corq_clini.serum_plot>=0.05]=""
corq_clini.serum_plot["Tauroursodeoxycholic acid","CB"]="***"
corq_clini.serum_plot["Methyl jasmonate","CB"]="***"
corq_clini.serum_plot["Leucyl-Valine","CB"]="***"
corq_clini.serum_plot["Sambutoxin","ALP"]="***"
#只展示部分代谢物，如下为不展示的代谢物
serum_noclini=c("Succinic acid semialdehyde","Docosahexaenoic acid","Pivampicillin","4-Hydroxy-1H-indole-3-acetonitrile","Methyl jasmonate","gamma-Calacorene","Gentamicin","Sambutoxin","Probucol","Dopexamine","Leucyl-Isoleucine","Oleamide","Neryl rhamnosyl-glucoside")
index=row.names(corq_clini.serum_plot) %in%  serum_noclini
index=which(index==TRUE)
library(pheatmap)
col3 <- colorRampPalette(c("#259042", "white", "#82388A"))
pheatmap(corR_clini.serum_plot[-index,],scale = "none",cluste_rows = T, cluster_cols= T, border=NA, display_numbers =corq_clini.serum_plot[-index,],
         fontsize_number = 12, number_color = "black", cellwidth = 20, cellheight =20,color=col3(100))


#将不显著的代谢物相关系数调整为0
corR_clini.serum_plot_v2=corR_clini.serum_plot
corR_clini.serum_plot_v2[corq_clini.serum[row.names(corR_clini.serum_plot),]>0.05]=0

corR_clini.serum_plot_v2[,1:10]=sapply(corR_clini.serum_plot_v2[,1:10],as.numeric)
corR_clini.serum_plot_v2=as.matrix(corR_clini.serum_plot_v2)
pheatmap(corR_clini.serum_plot_v2[-index,],scale = "none",cluste_rows = T, cluster_cols= T, border=NA, display_numbers =corq_clini.serum_plot[-index,],
         fontsize_number = 12, number_color = "black", cellwidth = 20, cellheight =20,color=col3(100))


##############特定代谢物相关性散点图##################
merge_PBCbase_plot=cbind(serum_0.7_log2[PBCbase,serum_PBCHC_dif_lm],metadata[PBCbase,18:24])
ggplot(data=merge_PBCbase_plot,aes(x=merge_PBCbase_plot[,42],y=merge_PBCbase_plot[,176]))+
  geom_jitter(width = 0.2,size=1,color="#A71E24")+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank())+
  labs(x="succinic acid",y="GGT",text=element_text(size = 8))+
  stat_smooth(method = "lm",se=FALSE,color="grey")


############可视化展示barplot t-value############
#1.PBCHC
lm_PBCHC_plot=lm_PBCHC[serum_PBCHC_dif_lm,]
lm_PBCHC_plot$metabolite=row.names(lm_PBCHC_plot)
lm_PBCHC_plot[which(lm_PBCHC_plot$tval_PBCHC<0),"change"]="HC"
lm_PBCHC_plot[which(lm_PBCHC_plot$tval_PBCHC>0),"change"]="PBCbase"
library(plyr)
library(forcats)
lm_PBCHC_plot <- lm_PBCHC_plot %>%  mutate(metabolite = fct_reorder(metabolite,tval_PBCHC))

#展现PBCHC显著变化的代谢物的tvalue
library(ggplot2)
ggplot(data = lm_PBCHC_plot,aes(x=metabolite,y=tval_PBCHC,fill=change))+
  geom_bar(stat = "identity",position = "identity",width = 0.2)+
  scale_fill_manual(values = c("#83AAF0","#F2644B"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  theme(legend.position =c(0.8,0.7))+
  theme(text = element_text(size = 6))+
  labs(fill="")+
  coord_flip()

#2.PBC resp-noresp nominally相关的代谢物(n=32)
lm_PBCresp_plot=filter(lm_PBCHC,pval_PBCresp.noresp<0.05)
lm_PBCresp_plot$metabolite=row.names(lm_PBCresp_plot)
lm_PBCresp_plot[which(lm_PBCresp_plot$tval_PBCresp.noresp<0),"change"]="UDCA noresponse"
lm_PBCresp_plot[which(lm_PBCresp_plot$tval_PBCresp.noresp>0),"change"]="UDCA response"
library(plyr)
library(forcats)
lm_PBCresp_plot <- lm_PBCresp_plot %>%  mutate(metabolite = fct_reorder(metabolite,tval_PBCresp.noresp))
library(ggplot2)
#展现nominally significant代谢物
ggplot(data = lm_PBCresp_plot,aes(x=metabolite,y=tval_PBCresp.noresp,fill=change))+
  geom_bar(stat = "identity",position = "identity",width = 0.6)+
  scale_fill_manual(values = c("#BFBFBF","#85ABD1"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  theme(legend.position ="null")+
  theme(text = element_text(size = 10))+
  labs(fill="")+
  coord_flip()

#2.scPBCbase HC
#在PBCHC显著改变且在亚临床PBC中nominally改变（P<0.05, FDR<0.15）
lm_scPBCHC_plot=filter(lm_PBCHC,pval_scPBCHC<0.05&qval_PBCHC<0.05)
lm_scPBCHC_plot$metabolite=row.names(lm_scPBCHC_plot)
lm_scPBCHC_plot[which(lm_scPBCHC_plot$tval_scPBCHC<0),"change"]="HC"
lm_scPBCHC_plot[which(lm_scPBCHC_plot$tval_scPBCHC>0),"change"]="scPBCbase"
library(plyr)
library(forcats)
lm_scPBCHC_plot <- lm_scPBCHC_plot %>%  mutate(metabolite = fct_reorder(metabolite,tval_scPBCHC))
library(ggplot2)
#展现nominally significant代谢物
ggplot(data = lm_scPBCHC_plot,aes(x=metabolite,y=tval_scPBCHC,fill=change))+
  geom_bar(stat = "identity",position = "identity",width = 0.2)+
  scale_fill_manual(values = c("#83AAF0","#C8B8DA"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  theme(legend.position =c(0.8,0.3))+
  theme(text = element_text(size = 7))+
  labs(fill="")+
  coord_flip()
############治疗前后可视化 boxplot*figure##################
# 1.boxplot图添加前后的连线
library(ggpubr)
#1.1 all PBCbase
data=data.frame(baseline=serum_0.7_log2[PBCpair_6m,"Succinic acid"],
                treat=serum_0.7_log2[PBC_6m_pair,"Succinic acid"]) 
ggpaired(data, cond1 = "baseline", cond2 = "treat",
         color = "condition", palette = "npg",ylab = "serum Succinic acid",line.color = "grey")+
  theme(legend.position = "none")
#1.2.cluster1
data=data.frame(baseline=serum_0.7_log2[PBCpair_6m_cluster1,"Succinic acid"],
                treat=serum_0.7_log2[PBC_6m_pair_cluster1,"Succinic acid"]) 
ggpaired(data, cond1 = "baseline", cond2 = "treat",
         color = "condition",palette = c("#088288","#ADD5DB"),ylab = "Succinic acid",line.color = "grey")+
  theme(legend.position = "none") 
#1.3.cluster2
data=data.frame(baseline=serum_0.7_log2[PBCpair_6m_cluster2,"Succinic acid"],
                treat=serum_0.7_log2[PBC_6m_pair_cluster2,"Succinic acid"]) 
ggpaired(data, cond1 = "baseline", cond2 = "treat",
         color = "condition",palette = c("#F19433","#FFD966"),ylab = "Succinic acid",line.color = "grey")+
  theme(legend.position = "none") 

#2.all group boxplot  HC_new
serum_0.7_log2_plot_treat=serum_0.7_log2[c(PBCpair_6m_cluster1,PBC_6m_pair_cluster1,PBCpair_6m_cluster2,PBC_6m_pair_cluster2,HC_new),]
serum_0.7_log2_plot_treat$group=c(rep("clt1.base",25),rep("clt1.6m",15),rep("clt2.base",25),rep("clt2.6m",15),rep("HC",111))
serum_0.7_log2_plot_treat$group=factor(serum_0.7_log2_plot_treat$group,levels = c("clt1.base","clt1.6m","clt2.base","clt2.6m","HC"))
library(ggplot2)
ggplot(data = serum_0.7_log2_plot_treat,aes(x=group,y=serum_0.7_log2_plot_treat$`N-Acetyl-L-methionine`,fill=group))+geom_boxplot()+
  labs(y="N-Acetyl-L-methionine")+
  scale_fill_manual(values = c("#088288","#ADD5DB","#F19433","#FFD966","#3F5BA7"))+
  theme(legend.position = "none")

############治疗前后可视化 均值heatmap*figure##################
#PBC治疗前后显著变化(FDR<0.05)的代谢物
#按照PBC整体展示
library(dplyr)
colnames(res_serum_pair.wilcox)
serum_PBCbase.6m_dif=res_serum_pair.wilcox %>% filter(.,PBCbase_6m.qval<0.05) %>% row.names(.)

serum_plot_treat=serum_0.7_log2_plot[c(PBCpair_6m,PBC_6m_pair,HC_new),]
serum_plot_treat_dif=serum_plot_treat[,c(serum_PBCbase.6m_dif,"group")]
serum_plot_treat_dif$group=serum_plot_treat_dif$group[,drop=TRUE]

serum_plot_treat_dif_mean=as.data.frame(matrix(0,ncol=1,nrow=3))

for (i in 1:30) {
  res.serum=aggregate(serum_plot_treat_dif[,i], by=list(type=serum_plot_treat_dif$group),mean)
  serum_plot_treat_dif_mean[,i+1]=res.serum$x
  row.names(serum_plot_treat_dif_mean)=res.serum$type
}
serum_plot_treat_dif_mean=serum_plot_treat_dif_mean[,-1]
colnames(serum_plot_treat_dif_mean)=serum_PBCbase.6m_dif

library(pheatmap)
pheatmap(t(serum_plot_treat_dif_mean),cluster_rows = T,cluster_cols = F,scale = "row",
         color = colorRampPalette(c("#259042", "white", "#82388A" ))(50),gaps_row = c(2,4),angle_col = "45",border_color = "white",fontsize = 7,
         cellwidth = 12,cellheight = 7)


##########pathway analysis及可视化###########
#输出数据到metaboanalyst
serum_0.7_log2UV_for_pathway=serum_0.7_log2_UV[c(PBCbase,HC),serum_PBCHC_dif_lm]
serum_0.7_log2UV_for_pathway=cbind(metadata_PBCHC[c(PBCbase,HC),"group"],serum_0.7_log2UV_for_pathway)
colnames(serum_0.7_log2UV_for_pathway)[1]="label"
serum_0.7_log2UV_for_pathway$label
setwd("D:/PBC multiomics data/analysis V6/metabolomics/serum metabolome/data/pathway analysis")
write.csv(serum_0.7_log2UV_for_pathway,"serum_0.7_log2UV_for_pathway.csv")

#计算不同通路的pathway intensity (with >3 metabolite)
Glycerophospholipid.metabolism=apply(serum[,c("PE(16:0/20:4(5Z,8Z,11Z,14Z))","PC(20:4(8Z,11Z,14Z,17Z)/20:4(8Z,11Z,14Z,17Z))","PC(22:5(7Z,10Z,13Z,16Z,19Z)/20:1(11Z))","PC(22:6(4Z,7Z,10Z,13Z,16Z,19Z)/20:1(11Z))",
                                              "LysoPC(P-18:1(9Z))","LysoPC(P-16:0)","LysoPC(20:3(5Z,8Z,11Z))","Glycerophosphocholine","Choline","Acetylcholine")],1,mean) %>% as.data.frame(.)
Pathway.activity=Glycerophospholipid.metabolism
colnames(Pathway.activity)="Glycerophospholipid metabolism"
Pathway.activity$"Arachidonic acid metabolism"=apply(serum[,c("PC(20:4(8Z,11Z,14Z,17Z)/20:4(8Z,11Z,14Z,17Z))","PC(22:5(7Z,10Z,13Z,16Z,19Z)/20:1(11Z))","PC(22:6(4Z,7Z,10Z,13Z,16Z,19Z)/20:1(11Z))",
                                                              "Arachidonic acid","20-Hydroxyeicosatetraenoic acid","5,6-DHET")],1,mean) 

Pathway.activity$"Primary bile acid biosynthesis"=apply(serum[,c("Glycine","Glycocholic acid(-)","Taurocholic acid","Chenodeoxycholic acid","Taurochenodesoxycholic acid")],1,mean)
Pathway.activity$"Biosynthesis of unsaturated fatty acids"=apply(serum[,c("Docosahexaenoic acid","Alpha-Linolenic acid","Arachidonic acid")],1,mean)
Pathway.activity$"Alanine, aspartate and glutamate metabolism"=apply(serum[,c("L-Aspartic acid","Pyruvic acid","gamma-Aminobutyric acid","Succinic acid semialdehyde","Succinic acid")],1,mean)
Pathway.activity$"Cysteine and methionine metabolism"=apply(serum[,c("L-Methionine","5'-Methylthioadenosine","L-Serine","2-Ketobutyric acid","3-Sulfinoalanine","Pyruvic acid")],1,mean)
Pathway.activity$"Butanoate metabolism"=apply(serum[,c("3-Hydroxybutyric acid","gamma-Aminobutyric acid","Succinic acid semialdehyde","Succinic acid")],1,mean)
Pathway.activity$"Propanoate metabolism"=apply(serum[,c("Succinic acid","2-Ketobutyric acid","2-Hydroxybutyric acid")],1,mean)
Pathway.activity$"Glycine, serine and threonine metabolism"=apply(serum[,c("Choline","Pyruvic acid","Glycine","L-Serine","2-Ketobutyric acid","Hydroxypyruvic acid","L-Threonine(+)")],1,mean)
Pathway.activity$"Valine, leucine and isoleucine biosynthesis"=apply(serum[,c("2-Ketobutyric acid","L-Isoleucine","L-Threonine(+)")],1,mean)
Pathway.activity$"Phenylalanine metabolism"=apply(serum[,c("L-Tyrosine","Hippuric acid","L-Phenylalanine(-)")],1,mean)
Pathway.activity$"Aminoacyl-tRNA biosynthesis"=apply(serum[,c("Glycine","L-Aspartic acid","L-Serine","L-Methionine","L-Isoleucine","(???)-Tryptophan(+)","L-Tyrosine",
                                                              "L-Phenylalanine(-)","L-Threonine(+)")],1,mean)
Pathway.activity$"Sphingolipid metabolism"=apply(serum[,c("L-Serine","Sphinganine","Phytosphingosine","Sphingosine")],1,mean)
Pathway.activity$"Glyoxylate and dicarboxylate metabolism"=apply(serum[,c("L-Serine","Pyruvic acid","Glycine","Hydroxypyruvic acid")],1,mean)
Pathway.activity$"Nicotinate and nicotinamide metabolism"=apply(serum[,c("L-Aspartic acid","1-Methylnicotinamide","N1-Methyl-4-pyridone-3-carboxamide")],1,mean)
row.names(Pathway.activity)=row.names(serum)
#读入metaboanalyst的pathway分析结果
setwd("D:/PBC multiomics data/analysis V6/metabolomics/serum metabolome/data")
pathway_res=read.csv("serum_pathway.csv",header = T,row.names = 1)
pathway_res$FDR=p.adjust(pathway_res$Raw.p,"fdr") #重新计算FDR值
pathway_res$pathway=row.names(pathway_res)
colnames(pathway_res)[4]="-logP"
library(ggplot2)
ggplot(pathway_res, aes(x=pathway, y=Diff_score)) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=Diff_score,color=Diff_score)) +
  geom_point( aes(color=Diff_score, size=pathway_res$`-logP`), alpha=0.9) +
  scale_size_continuous(range=c(2,6))+
  scale_color_gradientn(colours = colorRampPalette(c("#12528C", "white", "#B7232E"),bias=0.3)(100)) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

################根据class进行可视化展示############
#1.读入class分类
setwd("D:/PBC multiomics data/remove sc and 1y/serum/result")
serum_class=read.csv("lm_PBCHC.csv",header = T,row.names = 1)
colnames(serum_class)=serum_class[1,]
row.names(serum_class)=serum_class[,1]
serum_class=serum_class[-1,-1]
serum_class=serum_class[,1:4]
row.names(serum_class)[1:20]
row.names(serum_class)=row.names(lm_PBCHC)
#2.与linear model分析结果合并
serum_class=cbind(serum_class,lm_PBCHC)
#3.绘图
colnames(serum_class)[1:10]
lmt_PBCHC_dif_plot=serum_class[serum_PBCHC_dif_lm,c("class_renamed","tval_PBCHC")]
lmt_PBCHC_dif_plot$class_renamed=reorder(lmt_PBCHC_dif_plot$class_renamed,lmt_PBCHC_dif_plot$tval_PBCHC,FUN = median)
lmt_PBCHC_dif_plot$class_renamed=lmt_PBCHC_dif_plot$class_renamed[,drop=TRUE]
colnames(lmt_PBCHC_dif_plot)=c("class","tvalue")
lmt_PBCHC_dif_plot_serum=lmt_PBCHC_dif_plot
remove(lmt_PBCHC_dif_plot)
library(ggplot2)
ggplot(data=lmt_PBCHC_dif_plot_serum,aes(x=class,y=tvalue,fill=class))+
  geom_boxplot()+
  theme(axis.line = element_line(colour = "grey"),panel.background = element_blank(),legend.position = "none")+
  coord_flip() 

############绘制所有代谢物(n=441)superclass level 饼图#################
serum_class_simpli$SuperClass=as.factor(serum_class_simpli$SuperClass)
superclass_statis=summary(serum_class_simpli$SuperClass)
superclass_statis=as.data.frame(superclass_statis)
colnames(superclass_statis)="superclass"
superclass_statis$name=row.names(superclass_statis)
superclass_statis[1,2]="other"
superclass_statis=arrange(superclass_statis,desc(superclass_statis$superclass))
superclass_statis_new=superclass_statis[1:8,]
setwd("D:/PBC multiomics data/analysis V6/metabolomics/serum metabolome/data")
write.csv(superclass_statis_new,"superclass_statis_new.csv")
superclass_statis_new=read.csv("superclass_statis_new.csv",header = T,row.names = 1)

ggplot(superclass_statis_new, aes(x = '', y =superclass, fill = name)) +
  geom_bar(stat = 'identity', width = 1) +
  coord_polar(theta = 'y') +
  scale_fill_manual(values = c('#FFFFB3','#8DD3C7','#BEBADA', '#FB8072', '#80B1D3', '#FDB462',"#238FD0","#CBC402","#77B829"))
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
  theme( legend.title = element_blank()) +
  labs(x = '', y = '', title = 'Serum metabolome', fill = 'Superclass')
  
#######绘制血清代谢物来源(host, microbiota, co-metabolism等)的柱状图 #####
  setwd("D:/PBC multiomics data/analysis V6/metabolomics/serum metabolome/data")
  serum_origin=read.csv("serum_origin_statistics.csv",header = T,row.names = 1)
  library(ggplot2)
  library(dplyr)
  library(forcats)
  serum_origin <- serum_origin %>%  mutate(name = fct_reorder(name, origin))#按照另一列排列
  
  ggplot(data = serum_origin,aes(x=name,y=origin,fill=name))+
    geom_bar(stat = "identity",position = "identity",width = 0.6)+
    scale_fill_manual(values = c('#F2644B','#87C3DF','#02425A', '#C8B8DA', '#8167AD', '#CADEB8'))+
    theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
    theme(legend.position ="null")+
    theme(text = element_text(size = 10))+
    labs(fill="")+
    coord_flip()
#############小提琴图展示差异特征数目及系数分布###############

  serum_violin <- data.frame(
    group=c(rep("clt1clt2",7),rep("clt2HC",144),rep("clt1HC",158),rep("PBCHC",172) ),
    lm.tval=c( lm_PBCHC[serum_clt1clt2_dif,"tval_clt1clt2" ], 
               -lm_PBCHC[serum_clt2HC_dif,"tval_clt2HC"], 
               -lm_PBCHC[serum_clt1HC_dif,"tval_clt1HC"], 
               lm_PBCHC[serum_PBCHC_dif_lm,"tval_PBCHC"] )
  )
  serum_violin$group=factor(fecal_violin$group,levels = c("clt1clt2","clt2HC","clt1HC","PBCHC"))
  
  ggplot( data=serum_violin,aes(x=group, y=lm.tval, fill=group)) +
    geom_violin(scale = "count",width=1.8) +
    geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    scale_fill_viridis(discrete = TRUE) +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("serum metabolites") +
    xlab("")
