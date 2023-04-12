##############species_serum###################
#读入之前的数据
setwd("D:/PBC multiomics data/remove sc and 1y/multiomics/data")
merge_species_serum_old=read.csv("merge_species_serum.csv",header = T,row.names = 1)
colnames(merge_species_serum_old)=merge_species_serum_old[1,]
merge_species_serum_old=merge_species_serum_old[-1,]

#确定重合的species
library(dplyr)
species_PBCHC_cor.rep=intersect(colnames(merge_species_serum_old)[1:212],species_PBCHC_cor)
species_PBCHC_cor.uni=species_PBCHC_cor %in% species_PBCHC_cor.rep
species_PBCHC_cor.uni=which(species_PBCHC_cor.uni==FALSE)
#结果:species_PBCHC_cor 中184与186是新的，其他均已经算过species_serum correlation

#进行这两个species的计算
library(PResiduals)
merge_species_serum_new=cbind(serum_0.7[c(PBCbase_serum,HC_serum),serum_PBCHC_lm_df],species_abs_clr[c(PBCbase_serum,HC_serum),species_PBCHC_cor[c(184,186)]])
merge_species_serum_new=cbind(merge_species_serum_new,metadata[row.names(merge_species_serum_new),c("age","gender","BMI")])
colnames(merge_species_serum_new)[174]
str(merge_species_serum_new[,176:178]) 

corR_species_serum_new=sapply(174:175, function(x) sapply(1:173, function(y) {
  partial_Spearman(merge_species_serum_new[,x]|merge_species_serum_new[,y]~merge_species_serum_new[,"age"]+merge_species_serum_new[,"gender"]+merge_species_serum_new[,"BMI"],data=merge_species_serum_new)$TS$TB$ts
}))
corp_species_serum_new=sapply(174:175, function(x) sapply(1:173, function(y) {
  partial_Spearman(merge_species_serum_new[,x]|merge_species_serum_new[,y]~merge_species_serum_new[,"age"]+merge_species_serum_new[,"gender"]+merge_species_serum_new[,"BMI"],data=merge_species_serum_new)$TS$TB$pval
}))
corR_species_serum_new=t(corR_species_serum_new)
corp_species_serum_new=t(corp_species_serum_new)
colnames(corR_species_serum_new)=colnames(merge_species_serum_new)[1:173]
colnames(corp_species_serum_new)=colnames(merge_species_serum_new)[1:173]
row.names(corp_species_serum_new)=colnames(merge_species_serum_new)[174:175]
row.names(corR_species_serum_new)=colnames(merge_species_serum_new)[174:175]

##############species_fecal################### 
#
setwd("D:/PBC multiomics data/remove sc and 1y/multiomics/data")
merge_species_fecal_old=read.csv("merge_species_fecal.csv",header = T,row.names = 1)
colnames(merge_species_fecal_old)=merge_species_fecal_old[1,]
merge_species_fecal_old=merge_species_fecal_old[-1,]
colnames(merge_species_fecal_old)[221]
species_PBCHC_cor.rep_fecal=intersect(colnames(merge_species_fecal_old)[1:220],species_PBCHC_cor)
species_PBCHC_cor.uni_fecal=species_PBCHC_cor %in% species_PBCHC_cor.rep_fecal
species_PBCHC_cor.uni_fecal=which(species_PBCHC_cor.uni_fecal==FALSE)
#结果:所有species均不需要重新计算

##############serum_fecal################### 
merge_fecal_serum=cbind(fecal_0.7[c(PBCbase_fecal.serum,HC_fecal.serum),fecal_PBCHC_lm_df],serum_0.7[c(PBCbase_fecal.serum,HC_fecal.serum),serum_PBCHC_lm_df])
merge_fecal_serum=cbind(merge_fecal_serum,metadata[row.names(merge_fecal_serum),c("age","gender","BMI")])
colnames(merge_fecal_serum)[257]
corR_fecal_serum=sapply(1:83, function(x) sapply(84:256, function(y) {
  partial_Spearman(merge_fecal_serum[,x]|merge_fecal_serum[,y]~merge_fecal_serum[,"age"]+merge_fecal_serum[,"gender"]+merge_fecal_serum[,"BMI"],data=merge_fecal_serum)$TS$TB$ts
}))
corp_fecal_serum=sapply(1:83, function(x) sapply(84:256, function(y) {
  partial_Spearman(merge_fecal_serum[,x]|merge_fecal_serum[,y]~merge_fecal_serum[,"age"]+merge_fecal_serum[,"gender"]+merge_fecal_serum[,"BMI"],data=merge_fecal_serum)$TS$TB$pval
}))
colnames(corp_fecal_serum)=colnames(merge_fecal_serum)[1:83]
colnames(corR_fecal_serum)=colnames(merge_fecal_serum)[1:83]
row.names(corp_fecal_serum)=colnames(merge_fecal_serum)[84:256]
row.names(corR_fecal_serum)=colnames(merge_fecal_serum)[84:256]
corq_fecal_serum=p.adjust(corp_fecal_serum,"fdr") %>% matrix(.,ncol = 83)
colnames(corq_fecal_serum)=colnames(merge_fecal_serum)[1:83]
row.names(corq_fecal_serum)=colnames(merge_fecal_serum)[84:256]
library(dplyr)
write.csv(corp_fecal_serum,"corp_fecal_serum.csv")
write.csv(corq_fecal_serum,"corq_fecal_serum.csv")
write.csv(corR_fecal_serum,"corR_fecal_serum.csv")

###########species_clini HC_new *figure######################
#确定差异的species
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin")
maaslin_species_PBCHC=read.csv("maaslin_species PBCHC_new.csv",header = T,row.names = 1)
row.names(maaslin_species_PBCHC)=maaslin_species_PBCHC$species
maaslin_species_PBCHC_dif=filter(maaslin_species_PBCHC,qval<0.05) %>% row.names(.) #HC_new 
#与20% presence species合并
species_PBCHC_cor=intersect(maaslin_species_PBCHC_dif,species_rela_filter_0.2) #n=217
merge_clini_species=cbind(merge_clini_species[,1:10],species_abs_clr[row.names(merge_clini_species),species_PBCHC_cor])
merge_clini_species=cbind(merge_clini_species,metadata[row.names(merge_clini_species),c("age","gender","BMI")])
library(PResiduals)
corR_clini_species=sapply(1:10, function(x) sapply(11:227, function(y) {
  partial_Spearman(merge_clini_species[,x]|merge_clini_species[,y]~merge_clini_species[,"age"]+merge_clini_species[,"gender"]+merge_clini_species[,"BMI"],data=merge_clini_species)$TS$TB$ts
}))
corp_clini_species=sapply(1:10, function(x) sapply(11:227, function(y) {
  partial_Spearman(merge_clini_species[,x]|merge_clini_species[,y]~merge_clini_species[,"age"]+merge_clini_species[,"gender"]+merge_clini_species[,"BMI"],data=merge_clini_species)$TS$TB$pval
}))
corq_clini_species=p.adjust(corp_clini_species,"fdr") %>% matrix(.,ncol=10)
row.names(corp_clini_species)=colnames(merge_clini_species)[11:227]
row.names(corR_clini_species)=colnames(merge_clini_species)[11:227]
colnames(corp_clini_species)=colnames(merge_clini_species)[1:10]
colnames(corR_clini_species)=colnames(merge_clini_species)[1:10]
row.names(corq_clini_species)=colnames(merge_clini_species)[11:227]
colnames(corq_clini_species)=colnames(merge_clini_species)[1:10]
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/partial clinical association_new")
write.csv(corp_clini_species,"corP_clini_species.csv")
write.csv(corq_clini_species,"corq_clini_species.csv")
write.csv(corR_clini_species,"corR_clini_species.csv")
#1.clini-species 绘图 V1
#读入要选择展示的物种
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/partial clinical association_new")
cor_clini_species_plot=read.csv("cor_clini_species_plot.csv",header = T,row.names = 1)
corR_clini_species_plot=corR_clini_species[row.names(cor_clini_species_plot),]
corp_clini_species_plot=corp_clini_species[row.names(cor_clini_species_plot),]

corR_clini_species_plot[corp_clini_species_plot>0.05]=0
corp_clini_species_plot[corp_clini_species_plot<0.05]="*"
corp_clini_species_plot[corp_clini_species_plot>0.05]=""

newname=vector()
name=strsplit(row.names(corR_clini_species_plot),"s__")
for (i in 1:20) {
  newname[i]=name[[i]][2]
}
newname[2]="Veillonella spp"

row.names(corR_clini_species_plot)=newname
row.names(cor_clini_species_plot)=newname
corR_clini_species_plot=as.matrix(corR_clini_species_plot)
corR_clini_species_plot[,1:10]=sapply(corR_clini_species_plot[,1:10],as.numeric)

col3 <- colorRampPalette(c("#12528C", "white", "#B7232E"))

library(pheatmap)
pheatmap(corR_clini_species_plot,scale = "none",cluster_rows = TRUE, cluster_cols= TRUE, border=NA,display_numbers = corp_clini_species_plot,
         fontsize_number = 14, number_color = "white", cellwidth = 15, cellheight =15,color=col3(100))

#根据上述热图对species名称的聚类顺序调整response的热图顺序，使之次序一致(*最新的图中没有添加这个部分)
cor_clini_species_plot=cor_clini_species_plot[c(2,9,1,11,6,15,4,7,3,14,12,18,17,16,10,5,8,13,19,20),]
cor_clini_species_plot=as.data.frame(cor_clini_species_plot)
colnames(cor_clini_species_plot)="response"
row.names(cor_clini_species_plot)=row.names(corR_clini_species_plot)[c(2,9,1,11,6,15,4,7,3,14,12,18,17,16,10,5,8,13,19,20)]
cor_clini_species_plot=as.matrix(cor_clini_species_plot)
pheatmap(cor_clini_species_plot,scale = "none",cluster_rows = FALSE, cluster_cols= FALSE, border=NA,
         fontsize_number = 8, number_color = "black", cellwidth = 15, cellheight =15,color=col3(100))


#clini-species **figrue V2
#与临床指标(IgM,AST,ALT,AKP,GGT,TB,ALB,CB有相关性 FDR<0.1的species)
#上述species选择在D:\PBC multiomics data\analysis V6\metagenomics\result\species\partial clinical association_new\cor_clini_species.xlsx中corq_selection sheet
index=c(6,87,8,15,31,41,1,3,4,22,38,76,100,107,119,151,156,173,182)
corR_clini_species_plot=corR_clini_species[index,]
corp_clini_species_plot=corp_clini_species[index,]
corq_clini_species_plot=corq_clini_species[index,]

corR_clini_species_plot[corp_clini_species_plot>0.05]=0
corq_clini_species_plot[corq_clini_species_plot<0.1]="*"
corq_clini_species_plot[corq_clini_species_plot>0.1&corp_clini_species_plot<0.05]="#"
corq_clini_species_plot[corp_clini_species_plot>0.05]=""
newname=vector()
name=strsplit(row.names(corR_clini_species_plot),"s__")
for (i in 1:19) {
  newname[i]=name[[i]][2]
}
newname[1]="Veillonella spp"
newname[19]="Haemophilus spp"

row.names(corR_clini_species_plot)=newname

corR_clini_species_plot=as.matrix(corR_clini_species_plot)

col3 <- colorRampPalette(c("#12528C", "white", "#B7232E"))
colnames(corR_clini_species_plot)
library(pheatmap)#不展示IgA及IgG的相关性
pheatmap(corR_clini_species_plot[,-c(2,3,7)],scale = "none",cluster_rows = TRUE, cluster_cols= TRUE, border=NA,display_numbers = corq_clini_species_plot[,-c(2,3,7)],
         fontsize_number = 14, number_color = "white", cellwidth = 15, cellheight =15,color=col3(100))

##############pathway_clini#############
#高于20%样本出现的pathway
KEGG_path3_rela_filter_0.2=KEGG_path3_rela_filter
KEGG_path3_rela_filter_0.2[KEGG_path3_rela_filter_0.2<0.001]=0
KEGG_path3_rela_filter_0.2[KEGG_path3_rela_filter_0.2>0.001]=1
KEGG_path3_rela_filter_0.2<- KEGG_path3_rela_filter[which(rowSums(KEGG_path3_rela_filter_0.2) >=73),] %>% row.names(.)
KEGG_path3_PBCHC_cor=intersect(maaslin_KEGG_path3_PBCHC_dif,KEGG_path3_rela_filter_0.2)

merge_clini_KEGG_path3=cbind(metadata_PBCbase[PBCbase,clini],KEGG_path3_abs_clr[PBCbase,KEGG_path3_PBCHC_cor])
merge_clini_KEGG_path3=cbind(merge_clini_KEGG_path3,metadata[PBCbase,c("age","gender","BMI")])

corR_clini_KEGG_path3=sapply(1:10, function(x) sapply(11:39, function(y) {
  partial_Spearman(merge_clini_KEGG_path3[,x]|merge_clini_KEGG_path3[,y]~merge_clini_KEGG_path3[,"age"]+merge_clini_KEGG_path3[,"gender"]+merge_clini_KEGG_path3[,"BMI"],data=merge_clini_KEGG_path3)$TS$TB$ts
}))
corp_clini_KEGG_path3=sapply(1:10, function(x) sapply(11:39, function(y) {
  partial_Spearman(merge_clini_KEGG_path3[,x]|merge_clini_KEGG_path3[,y]~merge_clini_KEGG_path3[,"age"]+merge_clini_KEGG_path3[,"gender"]+merge_clini_KEGG_path3[,"BMI"],data=merge_clini_KEGG_path3)$TS$TB$pval
}))
corq_clini_KEGG_path3=p.adjust(corp_clini_KEGG_path3,"fdr") %>% matrix(.,ncol=10)
row.names(corp_clini_KEGG_path3)=colnames(merge_clini_KEGG_path3)[11:39]
row.names(corR_clini_KEGG_path3)=colnames(merge_clini_KEGG_path3)[11:39]
colnames(corp_clini_KEGG_path3)=colnames(merge_clini_KEGG_path3)[1:10]
colnames(corR_clini_KEGG_path3)=colnames(merge_clini_KEGG_path3)[1:10]
row.names(corq_clini_KEGG_path3)=colnames(merge_clini_KEGG_path3)[11:39]
colnames(corq_clini_KEGG_path3)=colnames(merge_clini_KEGG_path3)[1:10]
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/partial clinical association")
write.csv(corp_clini_KEGG_path3,"corP_clini_KEGG_path3.csv")
write.csv(corq_clini_KEGG_path3,"corq_clini_KEGG_path3.csv")
write.csv(corR_clini_KEGG_path3,"corR_clini_KEGG_path3.csv")

##############pathway serum#################
merge_KEGG_path3_serum=cbind(serum_0.7[c(PBCbase_serum,HC_serum),serum_PBCHC_lm_df],KEGG_path3_abs_clr[c(PBCbase_serum,HC_serum),KEGG_path3_PBCHC_cor])
merge_KEGG_path3_serum=cbind(merge_KEGG_path3_serum,metadata[row.names(merge_KEGG_path3_serum),c("age","gender","BMI")])


corR_KEGG_path3_serum=sapply(174:202, function(x) sapply(1:173, function(y) {
  partial_Spearman(merge_KEGG_path3_serum[,x]|merge_KEGG_path3_serum[,y]~merge_KEGG_path3_serum[,"age"]+merge_KEGG_path3_serum[,"gender"]+merge_KEGG_path3_serum[,"BMI"],data=merge_KEGG_path3_serum)$TS$TB$ts
}))
corp_KEGG_path3_serum=sapply(174:202, function(x) sapply(1:173, function(y) {
  partial_Spearman(merge_KEGG_path3_serum[,x]|merge_KEGG_path3_serum[,y]~merge_KEGG_path3_serum[,"age"]+merge_KEGG_path3_serum[,"gender"]+merge_KEGG_path3_serum[,"BMI"],data=merge_KEGG_path3_serum)$TS$TB$pval
}))
corR_KEGG_path3_serum=t(corR_KEGG_path3_serum)
corp_KEGG_path3_serum=t(corp_KEGG_path3_serum)
colnames(corR_KEGG_path3_serum)=colnames(merge_KEGG_path3_serum)[1:173]
colnames(corp_KEGG_path3_serum)=colnames(merge_KEGG_path3_serum)[1:173]
row.names(corp_KEGG_path3_serum)=colnames(merge_KEGG_path3_serum)[174:202]
row.names(corR_KEGG_path3_serum)=colnames(merge_KEGG_path3_serum)[174:202]


###############pathway fecal#################
merge_KEGG_path3_fecal=cbind(fecal_0.7[c(PBCbase_fecal,HC_fecal),fecal_PBCHC_lm_df],KEGG_path3_abs_clr[c(PBCbase_fecal,HC_fecal),KEGG_path3_PBCHC_cor])
merge_KEGG_path3_fecal=cbind(merge_KEGG_path3_fecal,metadata[row.names(merge_KEGG_path3_fecal),c("age","gender","BMI")])


corR_KEGG_path3_fecal=sapply(84:112, function(x) sapply(1:83, function(y) {
  partial_Spearman(merge_KEGG_path3_fecal[,x]|merge_KEGG_path3_fecal[,y]~merge_KEGG_path3_fecal[,"age"]+merge_KEGG_path3_fecal[,"gender"]+merge_KEGG_path3_fecal[,"BMI"],data=merge_KEGG_path3_fecal)$TS$TB$ts
}))
corp_KEGG_path3_fecal=sapply(84:112, function(x) sapply(1:83, function(y) {
  partial_Spearman(merge_KEGG_path3_fecal[,x]|merge_KEGG_path3_fecal[,y]~merge_KEGG_path3_fecal[,"age"]+merge_KEGG_path3_fecal[,"gender"]+merge_KEGG_path3_fecal[,"BMI"],data=merge_KEGG_path3_fecal)$TS$TB$pval
}))
corR_KEGG_path3_fecal=t(corR_KEGG_path3_fecal)
corp_KEGG_path3_fecal=t(corp_KEGG_path3_fecal)
colnames(corR_KEGG_path3_fecal)=colnames(merge_KEGG_path3_fecal)[1:83]
colnames(corp_KEGG_path3_fecal)=colnames(merge_KEGG_path3_fecal)[1:83]
row.names(corp_KEGG_path3_fecal)=colnames(merge_KEGG_path3_fecal)[84:112]
row.names(corR_KEGG_path3_fecal)=colnames(merge_KEGG_path3_fecal)[84:112]


#############combine with old result##############
setwd("D:/PBC multiomics data/remove sc and 1y/multiomics/result")
#==================================
#1.species_serum数据
#1.1 corp数据
corp_secies_serum_old=read.csv("corp_species_serum_dif.csv",header = T,row.names=1)
colnames(corp_secies_serum_old)=corp_secies_serum_old[1,]
corp_secies_serum_old=corp_secies_serum_old[-1,]
corp_secies_serum_old=corp_secies_serum_old[,species_PBCHC_cor.rep]
corp_secies_serum=rbind(corp_secies_serum_new,t(corp_secies_serum_old))
dim(corp_secies_serum)
colnames(corp_secies_serum)
#1.2.corR数据
corR_secies_serum_old=read.csv("corR_species_serum_dif.csv",header = T,row.names=1)
colnames(corR_secies_serum_old)=corR_secies_serum_old[1,]
corR_secies_serum_old=corR_secies_serum_old[-1,]
corR_secies_serum_old=corR_secies_serum_old[,species_PBCHC_cor.rep]
corR_secies_serum=rbind(corR_secies_serum_new,t(corR_secies_serum_old))
dim(corR_secies_serum)
colnames(corR_secies_serum)
#1.3.corq
corq_secies_serum=p.adjust(corp_secies_serum,"fdr") %>% matrix(.,ncol = 173)
row.names(corq_secies_serum)=row.names(corp_secies_serum)
colnames(corq_secies_serum)=colnames(corp_secies_serum)
#
setwd("D:/PBC multiomics data/analysis V6/multiomics/result/partial correlation")
write.csv(corp_secies_serum,"corp_secies_serum.csv")
write.csv(corq_secies_serum,"corq_secies_serum.csv")
write.csv(corR_secies_serum,"corR_secies_serum.csv")

#2.species_fecal数据
#2.1 corp数据
corp_species_fecal=read.csv("corp_species_fecal_dif.csv",header = T,row.names = 1)
colnames(corp_species_fecal)=corp_species_fecal[1,]
corp_species_fecal=corp_species_fecal[-1,]
corp_species_fecal=corp_species_fecal[,species_PBCHC_cor]
#2.2 corR数据
corR_species_fecal=read.csv("corR_species_fecal_dif.csv",header = T,row.names = 1)
colnames(corR_species_fecal)=corR_species_fecal[1,]
corR_species_fecal=corR_species_fecal[-1,]
corR_species_fecal=corR_species_fecal[,species_PBCHC_cor]
#2.3 corq
corq_species_fecal=p.adjust(as.matrix(corp_species_fecal),"fdr") %>% matrix(.,ncol = 212)
colnames(corq_species_fecal)=colnames(corp_species_fecal)
row.names(corq_species_fecal)=row.names(corp_species_fecal)

setwd("D:/PBC multiomics data/analysis V6/multiomics/result/partial correlation")
write.csv(corp_species_fecal,"corp_species_fecal.csv")
write.csv(corq_species_fecal,"corq_species_fecal.csv")
write.csv(corR_species_fecal,"corR_species_fecal.csv")

###############计算所有代谢物与shannon指数的关联##################
#读入shannon数据
shannon=read.csv("shannon.csv",header = T,row.names = 1)
#HC_new
HC_no=intersect(c("MXH063","MXH1022","MXH1076","MXH1136","MXH444","MXH448","MXH459",
                  "MXH496","MXH497","MXH532","MXH565","MXH779","MXH827","MXH840","MXH966","MXH976",
                  "MXH455","MXH605","MXH572"),HC_serum)
index=HC_serum %in% HC_no #为逻辑值向量
index=which(index==TRUE)#为位置索引
HC_serum_new=HC_serum[-index]
#1.serum_shannon HC_new
merge_shannon_serum=cbind(serum_0.7[c(PBCbase_serum,HC_serum_new),],shannon[c(PBCbase_serum,HC_serum_new),"shannon"])
merge_shannon_serum=cbind(merge_shannon_serum,metadata[c(PBCbase_serum,HC_serum_new),c("age","gender","BMI")])
colnames(merge_shannon_serum)[282]="shannon"
library(PResiduals)
corR_shannon_serum=sapply(282, function(x) sapply(1:281, function(y) {
  partial_Spearman(merge_shannon_serum[,x]|merge_shannon_serum[,y]~merge_shannon_serum[,"age"]+merge_shannon_serum[,"gender"]+merge_shannon_serum[,"BMI"],data=merge_shannon_serum)$TS$TB$ts
}))
corp_shannon_serum=sapply(282, function(x) sapply(1:281, function(y) {
  partial_Spearman(merge_shannon_serum[,x]|merge_shannon_serum[,y]~merge_shannon_serum[,"age"]+merge_shannon_serum[,"gender"]+merge_shannon_serum[,"BMI"],data=merge_shannon_serum)$TS$TB$pval
}))
corq_shannon_serum=p.adjust(corp_shannon_serum,"fdr") %>% matrix(.,ncol = 1)
colnames(corp_shannon_serum)="shannon"
colnames(corq_shannon_serum)="shannon"
colnames(corR_shannon_serum)="shannon"
row.names(corp_shannon_serum)=colnames(merge_shannon_serum)[1:281]
row.names(corq_shannon_serum)=colnames(merge_shannon_serum)[1:281]
row.names(corR_shannon_serum)=colnames(merge_shannon_serum)[1:281]

#2.fecal_shannon
merge_shannon_fecal=cbind(fecal_0.7[PBCbase_fecal,],shannon[PBCbase_fecal,"shannon"])
merge_shannon_fecal=cbind(merge_shannon_fecal,metadata[PBCbase_fecal,c("age","gender","BMI")])
colnames(merge_shannon_fecal)[834]="shannon"

corR_shannon_fecal=sapply(834, function(x) sapply(1:833, function(y) {
  partial_Spearman(merge_shannon_fecal[,x]|merge_shannon_fecal[,y]~merge_shannon_fecal[,"age"]+merge_shannon_fecal[,"gender"]+merge_shannon_fecal[,"BMI"],data=merge_shannon_fecal)$TS$TB$ts
}))
corp_shannon_fecal=sapply(834, function(x) sapply(1:833, function(y) {
  partial_Spearman(merge_shannon_fecal[,x]|merge_shannon_fecal[,y]~merge_shannon_fecal[,"age"]+merge_shannon_fecal[,"gender"]+merge_shannon_fecal[,"BMI"],data=merge_shannon_fecal)$TS$TB$pval
}))
corq_shannon_fecal=p.adjust(corp_shannon_fecal,"fdr") %>% matrix(.,ncol = 1)
colnames(corp_shannon_fecal)="shannon"
colnames(corq_shannon_fecal)="shannon"
colnames(corR_shannon_fecal)="shannon"
row.names(corp_shannon_fecal)=colnames(merge_shannon_fecal)[1:833]
row.names(corq_shannon_fecal)=colnames(merge_shannon_fecal)[1:833]
row.names(corR_shannon_fecal)=colnames(merge_shannon_fecal)[1:833]

#输出结果
setwd("D:/PBC multiomics data/analysis V6/multiomics/result/partial correlation HC_new")
write.csv(corp_shannon_fecal,"corp_shannon_fecal.csv")
write.csv(corq_shannon_fecal,"corq_shannon_fecal.csv")
write.csv(corR_shannon_fecal,"corR_shannon_fecal.csv")
write.csv(corp_shannon_serum,"corp_shannon_serum_new.csv")
write.csv(corq_shannon_serum,"corq_shannon_serum_new.csv")
write.csv(corR_shannon_serum,"corR_shannon_serum_new.csv")

###############计算差异代谢物与shannon指数的关联##################
#linear modeling
merge_shannon_serum_dif=cbind(serum_0.7_log2UV[c(PBCbase_serum,HC_serum),serum_PBCHC_lm_df],shannon[c(PBCbase_serum,HC_serum),"shannon"])
merge_shannon_serum_dif=cbind(merge_shannon_serum_dif,metadata[c(PBCbase_serum,HC_serum),c("age","gender","BMI")])
colnames(merge_shannon_serum_dif)[174]="shannon"
lm(shannon~age+BMI,data = merge_shannon_serum_dif) %>% summary(.)

###############选择显著的feature进行mediation################
#1.species
#1.1 species_clini P<0.01
corp_clini_species
species_sig=vector()
for (i in 1:212) {
  species_sig[i]=which(corp_clini_species[i,]<0.01) %>% length(.)
}
species_sig=which(species_sig>0)
species_sig=row.names(corp_clini_species)[species_sig]
#1.2 species_serum FDR<0.05
corq_secies_serum
species_sig_serum=vector()
for (i in 1:212) {
  species_sig_serum[i]=which(corq_secies_serum[i,]<0.05) %>% length(.)
}
species_sig_serum=which(species_sig_serum>0)
species_sig_serum=row.names(corq_secies_serum)[species_sig_serum]
#1.3 species_fecal FDR<0.05
corq_species_fecal
species_sig_fecal=vector()
for (i in 1:212) {
  species_sig_fecal[i]=which(corq_species_fecal[,i]<0.05) %>% length(.)
}
species_sig_fecal=which(species_sig_fecal>0)
species_sig_fecal=colnames(corq_species_fecal)[species_sig_fecal]

#final species
species_med=intersect(species_sig,species_sig_fecal)
species_med=intersect(species_med,species_sig_serum)

#2.serum
#2.1.serum_clini P<0.05
setwd("D:/PBC multiomics data/analysis V6/metabolomics/serum metabolome/result")
corp_clini_serum=read.csv("corP_clini.serum.csv",header = T,row.names = 1)
colnames(corp_clini_serum)=corp_clini_serum[1,]
corp_clini_serum=corp_clini_serum[-1,]

serum_sig=vector()
for (i in 1:173) {
  serum_sig[i]=which(corp_clini_serum[i,]<0.01) %>% length(.)
}
serum_sig=which(serum_sig>0)
serum_sig=row.names(corp_clini_serum)[serum_sig]

#2.2.serum_species FDR<0.05
serum_sig_species=vector()
for (i in 1:173) {
  serum_sig_species[i]=which(corq_secies_serum[,i]<0.05) %>% length(.)
}
serum_sig_species=which(serum_sig_species>0)
serum_sig_species=colnames(corq_secies_serum)[serum_sig_species]
serum_sig=paste(serum_sig,"serum",sep = "_")

#final serum
serum_med=intersect(serum_sig,serum_sig_species)


#3.fecal
#3.fecal_clini P<0.05
setwd("D:/PBC multiomics data/analysis V6/metabolomics/fecal metabolome/result")
corp_clini_fecal=read.csv("corP_clini.fecal.csv",header = T,row.names=1)
fecal_sig=vector()
for (i in 1:83) {
  fecal_sig[i]=which(corp_clini_fecal[i,]<0.05) %>% length(.)
}
fecal_sig=which(fecal_sig>0)
fecal_sig=row.names(corp_clini_fecal)[fecal_sig]

#4.clinical
#4.1 clinical serum
clini_serum=vector()
for (i in 1:10) {
  clini_serum[i]=which(corp_clini_serum[,i]<0.05) %>% length(.)
}
clini_serum=which(clini_serum>0)
clini_serum=colnames(corp_clini_serum)[clini_serum]
#4.2 clini species
clini_species=vector()
for (i in 1:10) {
  clini_species[i]=which(corp_clini_species[,i]<0.05) %>% length(.)
}
clini_species=which(clini_species>0)
clini_species=colnames(corp_clini_species)[clini_species]

#################mediation analysis####################
library(mediation)
mediation_clini_serum_species=cbind(metadata_PBCbase[PBCbase_serum,clini],serum_0.7_log2UV[PBCbase_serum,serum_med])
mediation_clini_serum_species=cbind(mediation_clini_serum_species,species_abs_clr[PBCbase_serum,species_med]) 

#修改行名
mediation_name=colnames(mediation_clini_serum_species)
mediationID=vector()
for (i in 1:162) {
  mediationID[i]=paste("feature",i,sep = "")
}

colnames(mediation_clini_serum_species)=mediationID
mediation_clini_serum_species=mediation_clini_serum_species[complete.cases(mediation_clini_serum_species),]
mediationID=as.vector(mediationID)
str(mediation_clini_serum_species)
#合并age,gender,BMI数据
mediation_clini_serum_species=cbind(mediation_clini_serum_species,metadata[row.names(mediation_clini_serum_species),c("age","gender","BMI")])
colnames(mediation_clini_serum_species)[163:165]

#====================================
#完整的model
setwd("D:/PBC multiomics data/analysis V6/multiomics/result/mediation analysis")
result=c()
for (i in names(mediation_clini_serum_species)[c(1,4:9)]) {
  for (j in names(mediation_clini_serum_species)[11:57]) {
    for (t in names(mediation_clini_serum_species)[58:162]) {
      newdata=mediation_clini_serum_species[,c(i,j,t,"age","gender","BMI")]
      colnames(newdata)=c("outcome","mediator","treat","age","gender","BMI")
      model1=lm(mediator~treat+age+gender+BMI,data=newdata)
      model2=lm(outcome~mediator+treat+age+gender+BMI,data = newdata)
      model3=mediate(model.m = model1,model.y = model2,treat ="treat",mediator ="mediator",boot=FALSE, sims=1000)
      result=rbind(result,c(i,j,t,summary(model3)$d0,summary(model3)$d0.p,summary(model3)$z0,summary(model3)$z0.p,summary(model3)$n0,summary(model3)$n0.p,summary(model3)$tau.coef,summary(model3)$tau.p))
      
    }
  }
}
write.csv(result,"result.csv")

result1=c()
for (i in names(mediation_clini_serum_species)[4:9]) {
  for (j in names(mediation_clini_serum_species)[40:57]) {
    for (t in names(mediation_clini_serum_species)[58:162]) {
      newdata=mediation_clini_serum_species[,c(i,j,t,"age","gender","BMI")]
      colnames(newdata)=c("outcome","mediator","treat","age","gender","BMI")
      model1=lm(mediator~treat+age+gender+BMI,data=newdata)
      model2=lm(outcome~mediator+treat+age+gender+BMI,data = newdata)
      model3=mediate(model.m = model1,model.y = model2,treat ="treat",mediator ="mediator",boot=FALSE, sims=1000)
      result1=rbind(result1,c(i,j,t,summary(model3)$d0,summary(model3)$d0.p,summary(model3)$z0,summary(model3)$z0.p,summary(model3)$n0,summary(model3)$n0.p,summary(model3)$tau.coef,summary(model3)$tau.p))
      
    }
  }
}
write.csv(result1,"result1.csv")

result2=c()
for (i in names(mediation_clini_serum_species)[7:9]) {
  for (j in names(mediation_clini_serum_species)[57]) {
    for (t in names(mediation_clini_serum_species)[58:162]) {
      newdata=mediation_clini_serum_species[,c(i,j,t,"age","gender","BMI")]
      colnames(newdata)=c("outcome","mediator","treat","age","gender","BMI")
      model1=lm(mediator~treat+age+gender+BMI,data=newdata)
      model2=lm(outcome~mediator+treat+age+gender+BMI,data = newdata)
      model3=mediate(model.m = model1,model.y = model2,treat ="treat",mediator ="mediator",boot=FALSE, sims=1000)
      result2=rbind(result2,c(i,j,t,summary(model3)$d0,summary(model3)$d0.p,summary(model3)$z0,summary(model3)$z0.p,summary(model3)$n0,summary(model3)$n0.p,summary(model3)$tau.coef,summary(model3)$tau.p))
      
    }
  }
}

result3=c()
for (i in names(mediation_clini_serum_species)[5:7]) {
  for (j in names(mediation_clini_serum_species)[11:39]) {
    for (t in names(mediation_clini_serum_species)[58:162]) {
      newdata=mediation_clini_serum_species[,c(i,j,t,"age","gender","BMI")]
      colnames(newdata)=c("outcome","mediator","treat","age","gender","BMI")
      model1=lm(mediator~treat+age+gender+BMI,data=newdata)
      model2=lm(outcome~mediator+treat+age+gender+BMI,data = newdata)
      model3=mediate(model.m = model1,model.y = model2,treat ="treat",mediator ="mediator",boot=FALSE, sims=1000)
      result3=rbind(result3,c(i,j,t,summary(model3)$d0,summary(model3)$d0.p,summary(model3)$z0,summary(model3)$z0.p,summary(model3)$n0,summary(model3)$n0.p,summary(model3)$tau.coef,summary(model3)$tau.p))
      
    }
  }
}
write.csv(result3,"result3.csv")
result4=c()
for (i in names(mediation_clini_serum_species)[8:9]) {
  for (j in names(mediation_clini_serum_species)[11:56]) {
    for (t in names(mediation_clini_serum_species)[58:162]) {
      newdata=mediation_clini_serum_species[,c(i,j,t,"age","gender","BMI")]
      colnames(newdata)=c("outcome","mediator","treat","age","gender","BMI")
      model1=lm(mediator~treat+age+gender+BMI,data=newdata)
      model2=lm(outcome~mediator+treat+age+gender+BMI,data = newdata)
      model3=mediate(model.m = model1,model.y = model2,treat ="treat",mediator ="mediator",boot=FALSE, sims=1000)
      result4=rbind(result4,c(i,j,t,summary(model3)$d0,summary(model3)$d0.p,summary(model3)$z0,summary(model3)$z0.p,summary(model3)$n0,summary(model3)$n0.p,summary(model3)$tau.coef,summary(model3)$tau.p))
      
    }
  }
}

write.csv(result4,"result4.csv")

write.csv(as.data.frame(mediation_name),"mediation name.csv")

#calculate serum associations with CB

###########筛选显著的mediation结果###################
#读入完整的结果
#1.修改行名为三个feature，去除重复计算的结果
result.comple=read.csv("result.complete.csv",header = T,row.names = 1)
result.name=paste(result.comple$species.name,result.comple$metabolite.name,sep = "|")
result.name=paste(result.name,result.comple$clinical.name,sep = "|")
result.comple$name=result.name
result.comple=result.comple[!duplicated(result.comple$name),] %>% as.data.frame(.)
row.names(result.comple)=result.comple$name
#2.选择species_clini P<0.01
library(reshape2)
library(reshape)
library(dplyr)
corp_clini_species_melt=melt(corp_clini_species)
corp_clini_species_melt_sig=corp_clini_species_melt %>% filter(.,value<0.01)
#2.1 species_name
species_clini_sig=corp_clini_species_melt_sig$X1 %>% as.character(.) 
species_clini_sig=species_clini_sig[!duplicated(species_clini_sig)]
#3.species_metabolite FDR<0.05
corq_secies_serum_melt=melt(corq_secies_serum[species_clini_sig,])
corq_secies_serum_melt_sig=corq_secies_serum_melt %>% filter(.,value<0.05)
#3.2 serum meltabolite-name
serum_species_sig=corq_secies_serum_melt_sig$X2 %>% as.character(.) 
serum_species_sig=serum_species_sig[!duplicated(serum_species_sig)]

#4.metabolite-clini
row.names(corp_clini_serum)=rownames(corp_fecal_serum)
corp_clini_serum[,1:10]=sapply(corp_clini_serum[,1:10],as.numeric)
corp_clini_serum=as.matrix(corp_clini_serum)
corp_clini_serum_melt=melt(corp_clini_serum[serum_species_sig,])
corp_clini_serum_melt_sig=filter(corp_clini_serum_melt,value<0.01)
#4.1
serum_clini_sig=corp_clini_serum_melt$X1 %>% as.character(.) 
serum_clini_sig=serum_clini_sig[!duplicated(serum_clini_sig)]

#5.依次倒退
#手动生成tri
write.csv(corp_clini_species_melt_sig,"corp_clini_species_melt_sig.csv")
write.csv(corp_clini_serum_melt_sig,"corp_clini_serum_melt_sig.csv")
#读入手动生成的tri数据(species-metabo-clini)
cor_tri=read.csv("cor tri.csv",header = T,row.names = 1)
cor_tri$species_serum=paste(cor_tri$species,cor_tri$serum,sep = "|")
corq_secies_serum_melt_sig$species_serumv2=paste(corq_secies_serum_melt_sig$X1,corq_secies_serum_melt_sig$X2,sep = "|")
cor_tri_index=cor_tri$species_serum %in% corq_secies_serum_melt_sig$species_serumv2
cor_tri_index=which(cor_tri_index==TRUE)
cor_tri$species_serum_clini=paste(cor_tri$species_serum,cor_tri$clini,sep = "|")
cor_tri_sig=cor_tri[cor_tri_index,]
write.csv(cor_tri_sig,"cor_tri_sig.csv")
#6.互相关联的mediation关系##
result.comple_cortri=result.comple$name %in% cor_tri_sig$species_serum_clini
result.comple_cortri=which(result.comple_cortri==TRUE)
result.comple_cortri_sig=result.comple[result.comple_cortri,]

write.csv(result.comple_cortri_sig,"result.comple_cortri_sig.csv")

#6.species-fecal-clini trio
setwd("D:/PBC multiomics data/analysis V6/multiomics/result/mediation analysis")
cor_tri_fecal=read.csv("cor tri fecal.csv",header = T,row.names = 1)
cor_tri_fecal$species_fecal=paste(cor_tri_fecal$species,cor_tri_fecal$fecal,sep = "|")
#6.1 species_fecal FDR<0.05
library(reshape)
library(reshape2)
library(dplyr)
corq_species_fecal_melt=melt(corq_species_fecal)
corq_species_fecal_melt_sig=filter(corq_species_fecal_melt,value<0.05)
corq_species_fecal_melt_sig$species_fecalv2=paste(corq_species_fecal_melt_sig$Var2,corq_species_fecal_melt_sig$Var1,sep = "|")  
cor_tri_fecal_index=cor_tri_fecal$species_fecal %in% corq_species_fecal_melt_sig$species_fecalv2
cor_tri_fecal_index=which(cor_tri_fecal_index==TRUE)
cor_tri_fecal$species_fecal_clini=paste(cor_tri_fecal$species_fecal,cor_tri_fecal$clini,sep = "|")
cor_tri_fecal_sig=cor_tri_fecal[cor_tri_fecal_index,]
write.csv(cor_tri_fecal_sig,"cor_tri_fecal_sig.csv")

#进行mediation分析的fecal metabolite和species
fecal_med=cor_tri_fecal_sig$fecal[!duplicated(cor_tri_fecal_sig$fecal)]
species_med_fecal=cor_tri_fecal_sig$species[!duplicated(cor_tri_fecal_sig$species)]
clini_med_fecal=cor_tri_fecal_sig$clini[!duplicated(cor_tri_fecal_sig$clini)]

#生成待进行fecal_mediation分析的数据
mediation_clini_fecal_species=cbind(metadata_PBCbase[PBCbase_fecal,clini_med_fecal],fecal_0.7_log2UV[PBCbase_fecal,fecal_med])
mediation_clini_fecal_species=cbind(mediation_clini_fecal_species,species_abs_clr[PBCbase_fecal,species_med_fecal]) 

#修改行名
mediation_name_fecal=colnames(mediation_clini_fecal_species)
mediationID_fecal=vector()
for (i in 1:34) {
  mediationID_fecal[i]=paste("feature",i,sep = "")
}

colnames(mediation_clini_fecal_species)=mediationID_fecal
mediation_clini_fecal_species=mediation_clini_fecal_species[complete.cases(mediation_clini_fecal_species),]
mediationID_fecal=as.vector(mediationID_fecal)
str(mediation_clini_fecal_species)
#合并age,gender,BMI数据
mediation_clini_fecal_species=cbind(mediation_clini_fecal_species,metadata[row.names(mediation_clini_fecal_species),c("age","gender","BMI")])
colnames(mediation_clini_fecal_species)[35:37]

#====================================
#完整的model
setwd("D:/PBC multiomics data/analysis V6/multiomics/result/mediation analysis")
result_fecal=c()
for (i in names(mediation_clini_fecal_species)[c(1:5)]) {
  for (j in names(mediation_clini_fecal_species)[6:17]) {
    for (t in names(mediation_clini_fecal_species)[18:34]) {
      newdata=mediation_clini_fecal_species[,c(i,j,t,"age","gender","BMI")]
      colnames(newdata)=c("outcome","mediator","treat","age","gender","BMI")
      model1=lm(mediator~treat+age+gender+BMI,data=newdata)
      model2=lm(outcome~mediator+treat+age+gender+BMI,data = newdata)
      model3=mediate(model.m = model1,model.y = model2,treat ="treat",mediator ="mediator",boot=FALSE, sims=1000)
      result_fecal=rbind(result_fecal,c(i,j,t,summary(model3)$d0,summary(model3)$d0.p,summary(model3)$z0,summary(model3)$z0.p,summary(model3)$n0,summary(model3)$n0.p,summary(model3)$tau.coef,summary(model3)$tau.p))
      
    }
  }
}
write.csv(result_fecal,"result_fecal.csv")
write.csv(as.data.frame(mediation_name_fecal),"mediation_name_fecal.csv")
#匹配行名后重新读入result_fecal
result_fecal_rename=read.csv("result_fecal.csv",header = T,row.names = 1)
result_fecal_rename$name=paste(result_fecal_rename$species.name,result_fecal_rename$metabolite.name,sep = "|")
result_fecal_rename$name=paste(result_fecal_rename$name,result_fecal_rename$clinical.name,sep = "|")
row.names(result_fecal_rename)=result_fecal_rename$name
#筛选cor_tri_fecal sig
result_fecal_index=result_fecal_rename$name %in% cor_tri_fecal_sig$species_fecal_clini
result_fecal_index=which(result_fecal_index==TRUE)

result_fecal_sig=result_fecal_rename[result_fecal_index,]
write.csv(result_fecal_sig,"result_fecal_sig.csv")

##重新计算血清中CB的关联分析
CB=read.csv("CB.csv",header = T,row.names = 1)
species_CB=CB$feature[1:12]
serum_CB=CB$feature[13:32]
mediation_CB_serum_species=cbind(metadata_PBCbase[PBCbase_serum,"CB"],serum_0.7_log2UV[PBCbase_serum,serum_CB])
mediation_CB_serum_species=cbind(mediation_CB_serum_species,species_abs_clr[PBCbase_serum,species_CB]) 

#修改行名

mediation_CB_serum_species=mediation_CB_serum_species[complete.cases(mediation_CB_serum_species),]
mediation_CB_serum_species=cbind(mediation_CB_serum_species,metadata[row.names(mediation_CB_serum_species),c("age","gender","BMI")])
str(mediation_CB_serum_species)
colnames(mediation_CB_serum_species)[1]="CB"
result_CB=c()
for (i in names(mediation_CB_serum_species)[1]) {
  for (j in names(mediation_CB_serum_species)[2:21]) {
    for (t in names(mediation_CB_serum_species)[22:33]) {
      newdata=mediation_CB_serum_species[,c(i,j,t,"age","gender","BMI")]
      colnames(newdata)=c("outcome","mediator","treat","age","gender","BMI")
      model1=lm(mediator~treat+age+gender+BMI,data=newdata)
      model2=lm(outcome~mediator+treat+age+gender+BMI,data = newdata)
      model3=mediate(model.m = model1,model.y = model2,treat ="treat",mediator ="mediator",boot=FALSE, sims=1000)
      result_CB=rbind(result_CB,c(i,j,t,summary(model3)$d0,summary(model3)$d0.p,summary(model3)$z0,summary(model3)$z0.p,summary(model3)$n0,summary(model3)$n0.p,summary(model3)$tau.coef,summary(model3)$tau.p))
      
    }
  }
}
write.csv(result_CB,"result_CB.csv")

#重命名后读入数据
result_CB_rename=read.csv("result_CB.csv",header = T,row.names = 1)
result_CB_rename$name=paste(result_CB_rename$species,result_CB_rename$metabolite,sep = "|")
result_CB_rename$name=paste(result_CB_rename$name,result_CB_rename$clinical,sep = "|")
result_CB_index=result_CB_rename$name %in% cor_tri_sig$species_serum_clini
result_CB_index=which(result_CB_index==TRUE)
result_CB_sig=result_CB_rename[result_CB_index,]
write.csv(result_CB_sig,"result_CB_sig.csv")

####合并所有trio sig的结果后进行mediation FDR校正###
result_triosig_mediation=read.csv("result triosig mediation.csv",header = T,row.names = 1)
result_triosig_mediation$FDR=p.adjust(result_triosig_mediation$d0.p,"fdr")
write.csv(result_triosig_mediation,"result_triosig_mediation_FDR.csv")


################inverse mediation################
#
species_invmed=result_triosig_mediation$species.name[!duplicated(result_triosig_mediation$species.name)] 
species_invmed=species_invmed[-20]
metabo_invmed=result_triosig_mediation$metabolite.name[!duplicated(result_triosig_mediation$metabolite.name)]
serum_invmed=grepl("serum",metabo_invmed) %>% metabo_invmed[.]
fecal_invmed=metabo_invmed[32:43]
clini_invmed=result_triosig_mediation$clinical.name[!duplicated(result_triosig_mediation$clinical.name)]
clini_invmed=clini_invmed[1:5]

#计算serum metabolite的mediation
Invmediation_serum=mediation_clini_serum_species[,c(clini_invmed,serum_invmed,species_invmed,"age","gender","BMI")]
result_invmed_serum=c()
for (i in names(Invmediation_serum)[1:5]) {
  for (j in names(Invmediation_serum)[6:35]) {
    for (t in names(Invmediation_serum)[36:54]) {
      newdata=Invmediation_serum[,c(i,j,t,"age","gender","BMI")]
      colnames(newdata)=c("treat","mediator","outcome","age","gender","BMI")
      model1=lm(mediator~treat+age+gender+BMI,data=newdata)
      model2=lm(outcome~mediator+treat+age+gender+BMI,data = newdata)
      model3=mediate(model.m = model1,model.y = model2,treat ="treat",mediator ="mediator",boot=FALSE, sims=1000)
      result_invmed_serum=rbind(result_invmed_serum,c(i,j,t,summary(model3)$d0,summary(model3)$d0.p,summary(model3)$z0,summary(model3)$z0.p,summary(model3)$n0,summary(model3)$n0.p,summary(model3)$tau.coef,summary(model3)$tau.p))
      
    }
  }
}

write.csv(result_invmed_serum,"result_invmed_serum.csv")
result_invmed_serum_rename=read.csv("result_invmed_serum.csv",header = T,row.names = 1)

#计算fecal metabolite的mediation
colnames(mediation_clini_fecal_species)=c(mediation_name_fecal,"age","gender","BMI")
Invmediation_fecal_feature=intersect(colnames(mediation_clini_fecal_species),c(clini_invmed,fecal_invmed,species_invmed,"age","gender","BMI"))
Invmediation_fecal=mediation_clini_fecal_species[,Invmediation_fecal_feature]

result_invmed_fecal=c()
for (i in names(Invmediation_fecal)[1:5]) {
  for (j in names(Invmediation_fecal)[6:17]) {
    for (t in names(Invmediation_fecal)[18:34]) {
      newdata=Invmediation_fecal[,c(i,j,t,"age","gender","BMI")]
      colnames(newdata)=c("treat","mediator","outcome","age","gender","BMI")
      model1=lm(mediator~treat+age+gender+BMI,data=newdata)
      model2=lm(outcome~mediator+treat+age+gender+BMI,data = newdata)
      model3=mediate(model.m = model1,model.y = model2,treat ="treat",mediator ="mediator",boot=FALSE, sims=1000)
      result_invmed_fecal=rbind(result_invmed_fecal,c(i,j,t,summary(model3)$d0,summary(model3)$d0.p,summary(model3)$z0,summary(model3)$z0.p,summary(model3)$n0,summary(model3)$n0.p,summary(model3)$tau.coef,summary(model3)$tau.p))
      
    }
  }
}
result_invmed_fecal_rename=as.data.frame(result_invmed_fecal)
colnames(result_invmed_fecal_rename)=colnames(result_invmed_serum_rename)

#合并血清和粪便的inverse mediation结果
result_invmed_serum_fecal=rbind(result_invmed_serum_rename,result_invmed_fecal_rename)
result_invmed_serum_fecal$name=paste(result_invmed_serum_fecal$species.name,result_invmed_serum_fecal$metabolite.name,sep = "|")
result_invmed_serum_fecal$name=paste(result_invmed_serum_fecal$name,result_invmed_serum_fecal$clinical.name,sep = "|")
row.names(result_invmed_serum_fecal)=result_invmed_serum_fecal$name
result_invmed_serum_fecal_repeat=result_invmed_serum_fecal[result_triosig_mediation$name,]

row.names(result_triosig_mediation)=result_triosig_mediation$name
result_pos_inv=cbind(result_triosig_mediation,result_invmed_serum_fecal_repeat)

result_pos_inv$FDR_inv=p.adjust(result_pos_inv[,18],"fdr")

write.csv(result_pos_inv,"result_pos_inv.csv")

############冲击图可视化#########
library(ggplot2)
install.packages("ggalluvial")
library(ggalluvial)

colnames(result_pos_inv)[14:26]=paste(colnames(result_pos_inv[14:26]),"inv",sep="|")
result_pos_inv$species.name_shortend=strsplit(result_pos_inv$species.name,"s__")
for (i in 1:459) {
  result_pos_inv$species.name_shortend[i]=result_pos_inv$species.name_shortend[[i]][2]
}
result_pos_inv[,27]=sapply(result_pos_inv[,27], as.character)

result_pos_inv$species.name_shortend

#选择某些代谢物进行可视化展示()
result_pos_inv_select=result_pos_inv[c(37,76,81,82,97,104,153,174,175,189,198,215,243,340),]

result_pos_inv_select$species.name_shortend[!duplicated(result_pos_inv_select$species.name_shortend)] #6
result_pos_inv_select$metabolite.name[!duplicated(result_pos_inv_select$metabolite.name)] #10
result_pos_inv_select$clinical.name[!duplicated(result_pos_inv_select$clinical.name)]#4

#颜色设置
color_species <- c('#8DD3C7','#8DD3C7', '#8DD3C7', '#8DD3C7', '#8DD3C7','#BEBADA')


color_metabo<- c('#80B1D3', '#FDB462',
                 '#FDB462', '#FDB462', '#EEA9B8', '#FDB462',"#d9ead3",'#e06666','#d9ead3', '#FFED6F')

color_clini <- c('#ead1dc', '#ead1dc', '#d5a6bd', '#a4c2f4')

#具体作图
ggplot(result_pos_inv_select, aes(axis1 = species.name_shortend, axis2 = metabolite.name, axis3 = clinical.name)) +
  geom_flow(aes(fill = metabolite.name)) +  #在这里，将所有连线的颜色按 miRNA 赋值，代表 miRNA 介导的关系
  scale_fill_manual(values = rev(color_metabo)) +
  geom_stratum(fill = c(color_species, color_metabo, color_clini)) +  #类似堆叠柱形图
  geom_text(stat = 'stratum', infer.label = TRUE, size = 2.5) +  #在各个区块中添加文字标签
  scale_x_discrete(limits = c(' species.name_shortend', 'metabolite.name', 'clinical.name')) +
  labs(x = '', y = '') +  #从这儿开始是一些主题参数等的设置
  theme(legend.position = 'none', panel.background = element_blank(),
        line = element_blank(), text = element_blank())


###########GMM analysis####################
#读入GMM数据
setwd("D:/PBC multiomics data/analysis V6/metagenomics/data")
GMM_abs_clr=read.csv("GMM_abs_clr.csv",header = T,row.names = 1)
#读入GMM maaslin 结果
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin")
maaslin_GMM_PBCHC=read.csv("maaslin GMM PBCHC.csv",header = T,row.names = 1)
colnames(maaslin_GMM_PBCHC)
#差异的GMM
maaslin_GMM_PBCHC_dif=maaslin_GMM_PBCHC %>% filter(.,qval<0.2) %>% row.names(.)
maaslin_GMM_PBCHC_dif_name=maaslin_GMM_PBCHC[maaslin_GMM_PBCHC_dif,"name"]
#species和GMM的关联性分析#
merge_species_GMM=cbind(species_abs_clr[c(PBCbase,HC),maaslin_species_PBCHC_dif],GMM_abs_clr[c(PBCbase,HC),maaslin_GMM_PBCHC_dif])
corP_species_GMM=sapply(287:319, function(x) sapply(1:286,function(y){cor.test(merge_species_GMM[,x],merge_species_GMM[,y],method = "spearman")$p.value}))
corR_species_GMM=sapply(287:319, function(x) sapply(1:286,function(y){cor.test(merge_species_GMM[,x],merge_species_GMM[,y],method = "spearman")$estimate}))
corq_species_GMM=p.adjust(corP_species_GMM,"fdr") %>% matrix(.,ncol = 33)
row.names(corP_species_GMM)=maaslin_species_PBCHC_dif
row.names(corq_species_GMM)=maaslin_species_PBCHC_dif
row.names(corR_species_GMM)=maaslin_species_PBCHC_dif
colnames(corR_species_GMM)=maaslin_GMM_PBCHC_dif_name
colnames(corP_species_GMM)=maaslin_GMM_PBCHC_dif_name
colnames(corq_species_GMM)=maaslin_GMM_PBCHC_dif_name
#输出文件
setwd("D:/PBC multiomics data/analysis V6/multiomics/result/nopartial spearman correlation")
write.csv(corP_species_GMM,"corP_species_GMM.csv")
write.csv(corq_species_GMM,"corq_species_GMM.csv")
write.csv(corR_species_GMM,"corR_species_GMM.csv")

################GMM species strong positive correlation visualization###################
corR_species_GMM_plot=corR_species_GMM
corR_species_GMM_plot[corR_species_GMM_plot<0.5]=0
corR_species_GMM_plot[corR_species_GMM_plot>=0.5]=1
#删除没有强正相关的行
species_GMM_cor=rowSums(corR_species_GMM_plot)
index=species_GMM_cor==0 
index=which(index==TRUE)
#删除没有强正相关的列
GMM_species_cor=colSums(corR_species_GMM_plot)
index_GMM=GMM_species_cor==0 
index_GMM=which(index_GMM==TRUE)

#另外不展示都是未知菌种的细菌
unclassfied=c(76,83,194,65,103,104,211,6,13,55,109,127,56,97,114,44,16,148,77,45,48,78,52,129,68,98,10,19,87,49,67,92,
  140,187,188,201,75,124,47,100,184,46,79,28,125,139,160,250,72,182,153,73,96,163,216,158,273,147,165,143,122,112,84,80,66,63,62,58,
  57,53,51,50,17,190,179,167,155,118,14,18,8,15,183,29,170,86,156,43,186,91,89,64)

#保留强正相关的关联
corR_species_GMM_plot=corR_species_GMM_plot[-c(index,unclassfied),-index_GMM]

#行名重命名
newname=vector()
name=strsplit(row.names(corR_species_GMM_plot),"s__")
for (i in 1:28) {
  newname[i]=name[[i]][2]
}
row.names(corR_species_GMM_plot)=newname
colnames(corR_species_GMM_plot)
col3 <- colorRampPalette(c( "white", "#B7232E"))
pheatmap(corR_species_GMM_plot,scale = "none",cluste_rows = F, cluster_cols= F, border="gray",
         fontsize_number = 12, number_color = "black", cellwidth = 20, cellheight =20,color=col3(100),
         angle_col = 45)

##############GMM serum###############

#高于20%样本出现的GMM
GMM_rela_filter_0.2=GMM_rela_filter
GMM_rela_filter_0.2[GMM_rela_filter_0.2<0.001]=0
GMM_rela_filter_0.2[GMM_rela_filter_0.2>0.001]=1
GMM_rela_filter_0.2<- GMM_rela_filter[which(rowSums(GMM_rela_filter_0.2) >=73),] %>% row.names(.)
GMM_PBCHC_cor=intersect(maaslin_GMM_PBCHC_dif,GMM_rela_filter_0.2)

merge_GMM_serum=cbind(serum_0.7[c(PBCbase_serum,HC_serum),serum_PBCHC_lm_df],GMM_abs_clr[c(PBCbase_serum,HC_serum),GMM_PBCHC_cor])
merge_GMM_serum=cbind(merge_GMM_serum,metadata[row.names(merge_GMM_serum),c("age","gender","BMI")])

corR_GMM_serum=sapply(174:199, function(x) sapply(1:173, function(y) {
  partial_Spearman(merge_GMM_serum[,x]|merge_GMM_serum[,y]~merge_GMM_serum[,"age"]+merge_GMM_serum[,"gender"]+merge_GMM_serum[,"BMI"],data=merge_GMM_serum)$TS$TB$ts
}))
corp_GMM_serum=sapply(174:199, function(x) sapply(1:173, function(y) {
  partial_Spearman(merge_GMM_serum[,x]|merge_GMM_serum[,y]~merge_GMM_serum[,"age"]+merge_GMM_serum[,"gender"]+merge_GMM_serum[,"BMI"],data=merge_GMM_serum)$TS$TB$pval
}))
corq_GMM_serum=p.adjust(corp_GMM_serum,"fdr") %>% matrix(.,ncol = 26)
row.names(corp_GMM_serum)=serum_PBCHC_lm_df
row.names(corq_GMM_serum)=serum_PBCHC_lm_df
row.names(corR_GMM_serum)=serum_PBCHC_lm_df
colnames(corp_GMM_serum)=maaslin_GMM_PBCHC[GMM_PBCHC_cor,"Name"]
colnames(corq_GMM_serum)=maaslin_GMM_PBCHC[GMM_PBCHC_cor,"Name"]
colnames(corR_GMM_serum)=maaslin_GMM_PBCHC[GMM_PBCHC_cor,"Name"]
write.csv(corp_GMM_serum,"corp_GMM_serum.csv")
write.csv(corq_GMM_serum,"corq_GMM_serum.csv")
write.csv(corR_GMM_serum,"corR_GMM_serum.csv")

##############GMM fecal###############
merge_GMM_fecal=cbind(fecal_0.7[c(PBCbase_fecal,HC_fecal),fecal_PBCHC_lm_df],GMM_abs_clr[c(PBCbase_fecal,HC_fecal),GMM_PBCHC_cor])
merge_GMM_fecal=cbind(merge_GMM_fecal,metadata[row.names(merge_GMM_fecal),c("age","gender","BMI")])

corR_GMM_fecal=sapply(84:109, function(x) sapply(1:83, function(y) {
  partial_Spearman(merge_GMM_fecal[,x]|merge_GMM_fecal[,y]~merge_GMM_fecal[,"age"]+merge_GMM_fecal[,"gender"]+merge_GMM_fecal[,"BMI"],data=merge_GMM_fecal)$TS$TB$ts
}))
corp_GMM_fecal=sapply(84:109, function(x) sapply(1:83, function(y) {
  partial_Spearman(merge_GMM_fecal[,x]|merge_GMM_fecal[,y]~merge_GMM_fecal[,"age"]+merge_GMM_fecal[,"gender"]+merge_GMM_fecal[,"BMI"],data=merge_GMM_fecal)$TS$TB$pval
}))
corq_GMM_fecal=p.adjust(corp_GMM_fecal,"fdr") %>% matrix(.,ncol = 26)
row.names(corp_GMM_fecal)=fecal_PBCHC_lm_df
row.names(corq_GMM_fecal)=fecal_PBCHC_lm_df
row.names(corR_GMM_fecal)=fecal_PBCHC_lm_df
colnames(corp_GMM_fecal)=maaslin_GMM_PBCHC[GMM_PBCHC_cor,"Name"]
colnames(corq_GMM_fecal)=maaslin_GMM_PBCHC[GMM_PBCHC_cor,"Name"]
colnames(corR_GMM_fecal)=maaslin_GMM_PBCHC[GMM_PBCHC_cor,"Name"]
write.csv(corp_GMM_fecal,"corp_GMM_fecal.csv")
write.csv(corq_GMM_fecal,"corq_GMM_fecal.csv")
write.csv(corR_GMM_fecal,"corR_GMM_fecal.csv")

##########选择显著的相关性构建network##########
library(reshape2)
library(reshape)
#1.species-serum FDR<0.05相关性
corq_species_serum_melt=melt(corq_secies_serum)
corR_species_serum_melt=melt(corR_secies_serum)
index=which(corq_species_serum_melt$value<0.05)
corR_species_serum_melt_sig=corR_species_serum_melt[index,]
cor_species_serum_melt_sig=cbind(corR_species_serum_melt_sig,corq_species_serum_melt[index,3])
colnames(cor_species_serum_melt_sig)[4]="FDR"


#2.species-fecal FDR<0.05相关性
corq_species_fecal_melt=melt(corq_species_fecal)
corR_species_fecal[,1:212]=sapply(corR_species_fecal[,1:212],as.numeric)
corR_species_fecal_matrix=as.matrix(corR_species_fecal)
corR_species_fecal_melt=melt(corR_species_fecal_matrix)
index=which(corq_species_fecal_melt$value<0.05)
corR_species_fecal_melt_sig=corR_species_fecal_melt[index,]
cor_species_fecal_melt_sig=cbind(corR_species_fecal_melt_sig,corq_species_fecal_melt[index,3])
colnames(cor_species_fecal_melt_sig)[4]="FDR"
str(corR_species_fecal_melt_sig)

#3.fecal-serum FDR<0.05
corq_fecal_serum_melt=melt(corq_fecal_serum)
corR_fecal_serum_melt=melt(corR_fecal_serum)
index=which(corq_fecal_serum_melt$value<0.05)
corR_fecal_serum_melt_sig=corR_fecal_serum_melt[index,]
cor_fecal_serum_melt_sig=cbind(corR_fecal_serum_melt_sig,corq_fecal_serum_melt[index,3])
colnames(cor_fecal_serum_melt_sig)[4]="FDR"

#4.合并上述显著相关性，生成cytoscape文件
cor_network=rbind(cor_species_serum_melt_sig,cor_species_fecal_melt_sig,
                  cor_fecal_serum_melt_sig)
setwd("D:/PBC multiomics data/analysis V6/multiomics/result")
write.csv(cor_network,"cor_network_v2.csv")

#===5.feature性质，用于cytoscape节点调整
row.names(maaslin_species_PBCHC)[1:5]
cor_network_feature_species=maaslin_species_PBCHC[species_PBCHC_cor,"PBCHC.coef"]%>%as.data.frame(.)
#重命名species名字
newname_v2=vector()
name=strsplit(species_PBCHC_cor,"s__")
for (i in 1:212) {
  newname_v2[i]=name[[i]][2]
}

cor_network_feature_species$class=c(rep("species",212))
cor_network_feature_species$newname=newname_v2
cor_network_feature_species$name=species_PBCHC_cor

cor_network_feature_serum=serum_res[serum_PBCHC_lm_df,"tval_PBCHC"] %>% as.data.frame(.)
row.names(cor_network_feature_serum)=serum_PBCHC_lm_df
cor_network_feature_serum$class=c(rep("serum",173))
cor_network_feature_serum$newname=serum_PBCHC_lm_df
cor_network_feature_serum$name=serum_PBCHC_lm_df

colnames(cor_network_feature_serum)

cor_network_feature_fecal=fecal_res[fecal_PBCHC_lm_df,"tval_PBCHC"] %>% as.data.frame(.)
row.names(cor_network_feature_fecal)=fecal_PBCHC_lm_df
cor_network_feature_fecal$class=c(rep("fecal",83))
cor_network_feature_fecal$newname=fecal_PBCHC_lm_df
cor_network_feature_fecal$name=fecal_PBCHC_lm_df
colnames(cor_network_feature_fecal)

cor_network_feature=rbind(cor_network_feature_species,cor_network_feature_serum,cor_network_feature_fecal)
cor_network_feature[which(cor_network_feature$.<0),"change"]="HC"
cor_network_feature[which(cor_network_feature$.>0),"change"]="PBCbase"

setwd("D:/PBC multiomics data/analysis V6/multiomics/result")
write.csv(cor_network_feature,"cor_network_feature_v3.csv")

############可视化网络中的feature rank############
#读入cytoscape生成的node rank 文件
setwd("D:/PBC multiomics data/analysis V6/multiomics/result")
node_rank=read.csv("node rank file v2.csv",header = T,row.names = 1)
row.names(node_rank)=node_rank$node_name
node_rank$class=as.factor(node_rank$class)
library(dplyr)
node_rank=arrange(node_rank,desc(MCC))

library(plyr)
library(forcats)
node_rank <- node_rank %>%  mutate(node_name = fct_reorder(node_name,MCC))

#展现top100 feature 按照MCC排序
library(ggplot2)
ggplot(data = node_rank[1:50,],aes(x=node_name,y=MCC,fill=class))+
  geom_bar(stat = "identity",position = "identity",width = 0.2)+
  scale_fill_manual(values = c("#6F6AB0","#EA5B67","#179C75"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  theme(legend.position =c(0.8,0.7))+
  theme(text = element_text(size = 8))+
  labs(fill="")+
  coord_flip()


################进行rmcorr分析治疗前后的相关性#################
#KO_gene和粪便代谢物
PBCpair_6m_fecal=intersect(PBCbase_fecal,PBCpair_6m)

KO_gene_fecal_rmcor=cbind(KO_gene_abs_clr[c(PBCpair_6m,PBC_6m_pair),row.names(KO_gene_rela_filter)],metadata[c(PBCpair_6m,PBC_6m_pair),"subjectID"])
KO_gene_fecal_rmcor=cbind(KO_gene_fecal_rmcor,fecal_0.7_log2UV[c(c(PBCpair_6m,PBC_6m_pair)),])

fecal_Na=which(is.na(KO_gene_fecal_rmcor[,10369]) ==TRUE) #生成缺失粪便代谢物数据的行索引
KO_gene_fecal_rmcor=KO_gene_fecal_rmcor[-c(fecal_Na),]
grep("MXC112",row.names(KO_gene_fecal_rmcor))#移除对应的初治或者治疗后样本:MXC015-2" "MXC050-2" "MXC112-2
KO_gene_fecal_rmcor=KO_gene_fecal_rmcor[-c(8,13,36),]#移除对应的初治样本:MXC015-2" "MXC050-2" "MXC112-2

colnames(KO_gene_fecal_rmcor)[9537]="subjectID"
KO_gene_fecal_rmcor$subjectID=as.factor(KO_gene_fecal_rmcor$subjectID)
KO_gene_fecal_rmcor$subjectID=KO_gene_fecal_rmcor$subjectID[,drop=TRUE]

rmcorr(subjectID,KO_gene_fecal_rmcor[,"K01023"],KO_gene_fecal_rmcor[,"Cholesterol sulfate"],KO_gene_fecal_rmcor)$p

#genus和粪便代谢物
genus_fecal_rmcor=cbind(genus_abs_clr[c(PBCpair_6m,PBC_6m_pair),row.names(genus_rela_filter)],metadata[c(PBCpair_6m,PBC_6m_pair),"subjectID"])
genus_fecal_rmcor=cbind(genus_fecal_rmcor,fecal_0.7_log2UV[c(c(PBCpair_6m,PBC_6m_pair)),])

fecal_Na=which(is.na(genus_fecal_rmcor[,1716]) ==TRUE) #生成缺失粪便代谢物数据的行索引
genus_fecal_rmcor=genus_fecal_rmcor[-c(fecal_Na),]
grep("MXC112",row.names(genus_fecal_rmcor))#移除对应的初治或者治疗后样本:MXC015-2" "MXC050-2" "MXC112-2
genus_fecal_rmcor=genus_fecal_rmcor[-c(8,13,36),]#移除对应的初治样本:MXC015-2" "MXC050-2" "MXC112-2

colnames(genus_fecal_rmcor)[884]="subjectID"
genus_fecal_rmcor$subjectID=as.factor(genus_fecal_rmcor$subjectID)
genus_fecal_rmcor$subjectID=genus_fecal_rmcor$subjectID[,drop=TRUE]

rmcorr(subjectID,genus_fecal_rmcor[,"k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides"],genus_fecal_rmcor[,"Cholesterol sulfate"],genus_fecal_rmcor)

