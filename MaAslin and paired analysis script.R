############1.species#############
#修改名字以与masslin分析的结果match,因为字符的问题;
species_name=row.names(species_rela_filter) %>% as.data.frame()
speciesID=vector()
for (i in 1:3854) {
  speciesID[i]=paste("species",row.names(species_name)[i],sep = "")
}
row.names(species_name)=speciesID

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin")
#1.1 PBCbase and HC_new
species_abs_clr_masslin=species_abs_clr[c(PBCbase,HC_new),row.names(species_rela_filter)]
colnames(species_abs_clr_masslin)=speciesID
write.csv(species_abs_clr_masslin,"species clr PBCbase HC_new masslin.csv")
remove(species_abs_clr_masslin)
metadata_PBCHC.new=metadata_PBCHC[c(PBCbase,HC_new),c("age","gender","BMI","group")]
#metadata
setwd("D:/PBC multiomics data/analysis V6/metagenomics/data")
write.csv(metadata_PBCHC.new,"metadata_PBCHC.new.maaslin.csv")
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/species clr PBCbase HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/data/metadata_PBCHC.new.maaslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/MAaSLIN2_species_PBC_HC.new',
                     fixed_effects = c("age","gender","BMI","group"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin")
write.csv(species_name,"species_name.csv")

#1.2 PBCbase responder and non-responder
species_abs_clr_masslin=species_abs_clr[c(PBCbase_resp_6m,PBCbase_noresp_6m),row.names(species_rela_filter)]
colnames(species_abs_clr_masslin)=speciesID
write.csv(species_abs_clr_masslin,"species clr PBCbase resp noresp masslin.csv")
remove(species_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/species clr PBCbase resp noresp masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove sc and 1y/metagenomics/data/genus/metadata PBCbase resp noresp masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/MAaSLIN2_species_PBCbase_resp_noresp',
                     fixed_effects = c("age","gender","BMI","Response_6m"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)
#1.3 scPBCbase and HC_all
species_abs_clr_masslin=species_abs_clr[c(scPBCbase,HC),row.names(species_rela_filter)]
colnames(species_abs_clr_masslin)=speciesID
write.csv(species_abs_clr_masslin,"species clr scPBCbase HC masslin.csv")
remove(species_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/species clr scPBCbase HC masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove 1y and keep sc/metagenomics/result/species/maaslin/metadata scPBCbase HC.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/MAaSLIN2_species_scPBCbase_HC',
                     fixed_effects = c("age","gender","BMI","group"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)

remove(fit_data)

#1.4 clt1 vs HC_new
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin")
species_abs_clr_masslin=species_abs_clr[c(PBCbase_cluster1,HC_new),row.names(species_rela_filter)]
colnames(species_abs_clr_masslin)=speciesID
write.csv(species_abs_clr_masslin,"species clr clt1 HC_new masslin.csv")
write.csv(metadata_PBCHC[c(PBCbase_cluster1,HC_new),],"metadata clt1 HC_new masslin.csv")
remove(species_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/species clr clt1 HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 HC_new masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/MAaSLIN2_species_clt1_HC.new',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
#1.5 clt2 vs HC_new
species_abs_clr_masslin=species_abs_clr[c(PBCbase_cluster2,HC_new),row.names(species_rela_filter)]
colnames(species_abs_clr_masslin)=speciesID
write.csv(species_abs_clr_masslin,"species clr clt2 HC_new masslin.csv")
write.csv(metadata_PBCHC[c(PBCbase_cluster2,HC_new),],"metadata clt2 HC_new masslin.csv")
remove(species_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/species clr clt2 HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt2 HC_new masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/MAaSLIN2_species_clt2_HC.new',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin")
#1.6 clt1 vs clt2
species_abs_clr_masslin=species_abs_clr[c(PBCbase_cluster1,PBCbase_cluster2),row.names(species_rela_filter)]
colnames(species_abs_clr_masslin)=speciesID
write.csv(species_abs_clr_masslin,"species clr clt1 clt2 masslin.csv")
write.csv(metadata_PBCHC[c(PBCbase_cluster1,PBCbase_cluster2),],"metadata clt1 clt2 masslin.csv")
remove(species_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/species clr clt1 clt2 masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 clt2 masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/MAaSLIN2_species_clt1_clt2',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)

#1.7.1 PBC_pair vs PBC_6m_pair #selective species
res_species_pair.wilcox=as.data.frame(matrix(0,ncol=1,nrow=286))
for (i in 1:286) {
  res_species_pair.wilcox$PBCbase_6m.pval[i]=wilcox.test(species_abs_clr[PBC_6m_pair,maaslin_species_PBCHC_dif][,i],
                                                         species_abs_clr[PBCpair_6m,maaslin_species_PBCHC_dif][,i],paired = T)$p.value
}
res_species_pair.wilcox$PBCbase_6m.qval=p.adjust(res_species_pair.wilcox$PBCbase_6m.pval,"fdr")
res_species_pair.wilcox$PBCbase.abundance=apply(species_rela_filter[maaslin_species_PBCHC_dif,PBCpair_6m],1,mean)
res_species_pair.wilcox$PBC6m.abundance=apply(species_rela_filter[maaslin_species_PBCHC_dif,PBC_6m_pair],1,mean)
res_species_pair.wilcox$change=res_species_pair.wilcox$PBC6m.abundance-res_species_pair.wilcox$PBCbase.abundance
row.names(res_species_pair.wilcox)=maaslin_species_PBCHC_dif

#1.7.2 PBC_pair vs PBC_6m_pair #all species
res_species_pair.wilcox=as.data.frame(matrix(0,ncol=1,nrow=3854))
for (i in 1:3854) {
  res_species_pair.wilcox$PBCbase_6m.pval[i]=wilcox.test(species_abs_clr[PBC_6m_pair,row.names(species_rela_filter)][,i],
                                                         species_abs_clr[PBCpair_6m,row.names(species_rela_filter)][,i],paired = T)$p.value
}
res_species_pair.wilcox$PBCbase_6m.qval=p.adjust(res_species_pair.wilcox$PBCbase_6m.pval,"fdr")
res_species_pair.wilcox$PBCbase.abundance=apply(species_rela_filter[row.names(species_rela_filter),PBCpair_6m],1,mean)
res_species_pair.wilcox$PBC6m.abundance=apply(species_rela_filter[row.names(species_rela_filter),PBC_6m_pair],1,mean)
res_species_pair.wilcox$change=res_species_pair.wilcox$PBC6m.abundance-res_species_pair.wilcox$PBCbase.abundance
res_species_pair.wilcox$PBCHC_coef=maaslin_species_PBCHC[row.names(species_rela_filter),"coef"]
res_species_pair.wilcox$PBCHC_pval=maaslin_species_PBCHC[row.names(species_rela_filter),"pval"]
res_species_pair.wilcox$PBCHC_qval=maaslin_species_PBCHC[row.names(species_rela_filter),"qval"]

row.names(res_species_pair.wilcox)=row.names(species_rela_filter)
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/paired_new")
write.csv(res_species_pair.wilcox,"res_species_pair.wilcox.csv")

#1.8.1 clt1_pair vs clt1_6m_pair #selective species
res_species_pair.wilcox_clt1=as.data.frame(matrix(0,ncol=1,nrow=2372))
for (i in 1:2372) {
  res_species_pair.wilcox_clt1$clt1base_6m.pval[i]=wilcox.test(species_abs_clr[PBCpair_6m_cluster1,maaslin_species_clt1_pair][,i],
                                                         species_abs_clr[PBC_6m_cluster1,maaslin_species_clt1_pair][,i],paired = T)$p.value
}
res_species_pair.wilcox_clt1$clt1base_6m.qval=p.adjust(res_species_pair.wilcox_clt1$clt1base_6m.pval,"fdr")
res_species_pair.wilcox_clt1$clt1base.abundance=apply(species_rela_filter[maaslin_species_clt1_pair,PBCpair_6m_cluster1],1,mean)
res_species_pair.wilcox_clt1$clt1.6m.abundance=apply(species_rela_filter[maaslin_species_clt1_pair,PBC_6m_cluster1],1,mean)
res_species_pair.wilcox_clt1$change=res_species_pair.wilcox_clt1$clt1.6m.abundance-res_species_pair.wilcox_clt1$clt1base.abundance
res_species_pair.wilcox_clt1$clt1HC_coef=maaslin_species_clt1HC[maaslin_species_clt1_pair,"coef"]
res_species_pair.wilcox_clt1$clt1HC_qval=maaslin_species_clt1HC[maaslin_species_clt1_pair,"qval"]
res_species_pair.wilcox_clt1$clt1clt2_coef=maaslin_species_clt1clt2[maaslin_species_clt1_pair,"coef"]
res_species_pair.wilcox_clt1$clt1clt2_qval=maaslin_species_clt1clt2[maaslin_species_clt1_pair,"qval"]
row.names(res_species_pair.wilcox_clt1)=maaslin_species_clt1_pair

#1.8 clt1_pair vs clt1_6m_pair #all species
res_species_pair.wilcox_clt1=as.data.frame(matrix(0,ncol=1,nrow=3854))
for (i in 1:3854) {
  res_species_pair.wilcox_clt1$clt1base_6m.pval[i]=wilcox.test(species_abs_clr[PBCpair_6m_cluster1,row.names(species_rela_filter)][,i],
                                                               species_abs_clr[PBC_6m_cluster1,row.names(species_rela_filter)][,i],paired = T)$p.value
}
res_species_pair.wilcox_clt1$clt1base_6m.qval=p.adjust(res_species_pair.wilcox_clt1$clt1base_6m.pval,"fdr")
res_species_pair.wilcox_clt1$clt1base.abundance=apply(species_rela_filter[row.names(species_rela_filter),PBCpair_6m_cluster1],1,mean)
res_species_pair.wilcox_clt1$clt1.6m.abundance=apply(species_rela_filter[row.names(species_rela_filter),PBC_6m_cluster1],1,mean)
res_species_pair.wilcox_clt1$change=res_species_pair.wilcox_clt1$clt1.6m.abundance-res_species_pair.wilcox_clt1$clt1base.abundance
res_species_pair.wilcox_clt1$clt1HC_coef=-maaslin_species_clt1HC[row.names(species_rela_filter),"coef"]
res_species_pair.wilcox_clt1$clt1HC_qval=maaslin_species_clt1HC[row.names(species_rela_filter),"qval"]
res_species_pair.wilcox_clt1$clt1clt2_coef=-maaslin_species_clt1clt2[row.names(species_rela_filter),"coef"]
res_species_pair.wilcox_clt1$clt1clt2_qval=maaslin_species_clt1clt2[row.names(species_rela_filter),"qval"]
row.names(res_species_pair.wilcox_clt1)=row.names(species_rela_filter)
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/paired_new")
write.csv(res_species_pair.wilcox_clt1,"res_species_pair.wilcox_clt1.csv")



#1.9.1 clt2pair vs clt2_6m_pair #select species
res_species_pair.wilcox_clt2=as.data.frame(matrix(0,ncol=1,nrow=2476))
for (i in 1:2476) {
  res_species_pair.wilcox_clt2$clt2base_6m.pval[i]=wilcox.test(species_abs_clr[PBCpair_6m_cluster2,maaslin_species_clt2_pair][,i],
                                                               species_abs_clr[PBC_6m_cluster2,maaslin_species_clt2_pair][,i],paired = T)$p.value
}
res_species_pair.wilcox_clt2$clt2base_6m.qval=p.adjust(res_species_pair.wilcox_clt2$clt2base_6m.pval,"fdr")
res_species_pair.wilcox_clt2$clt2base.abundance=apply(species_rela_filter[maaslin_species_clt2_pair,PBCpair_6m_cluster2],1,mean)
res_species_pair.wilcox_clt2$clt2.6m.abundance=apply(species_rela_filter[maaslin_species_clt2_pair,PBC_6m_cluster2],1,mean)
res_species_pair.wilcox_clt2$change=res_species_pair.wilcox_clt2$clt2.6m.abundance-res_species_pair.wilcox_clt2$clt2base.abundance
res_species_pair.wilcox_clt2$clt2HC_coef=-maaslin_species_clt2HC[maaslin_species_clt2_pair,"coef"]
res_species_pair.wilcox_clt2$clt2HC_qval=maaslin_species_clt2HC[maaslin_species_clt2_pair,"qval"]
res_species_pair.wilcox_clt2$clt1clt2_coef=maaslin_species_clt1clt2[maaslin_species_clt2_pair,"coef"]
res_species_pair.wilcox_clt2$clt1clt2_qval=maaslin_species_clt1clt2[maaslin_species_clt2_pair,"qval"]
row.names(res_species_pair.wilcox_clt2)=maaslin_species_clt2_pair

#1.9.2 clt2pair vs clt2_6m_pair #all species
res_species_pair.wilcox_clt2=as.data.frame(matrix(0,ncol=1,nrow=3854))
for (i in 1:3854) {
  res_species_pair.wilcox_clt2$clt2base_6m.pval[i]=wilcox.test(species_abs_clr[PBCpair_6m_cluster2,row.names(species_rela_filter)][,i],
                                                               species_abs_clr[PBC_6m_cluster2,row.names(species_rela_filter)][,i],paired = T)$p.value
}
res_species_pair.wilcox_clt2$clt2base_6m.qval=p.adjust(res_species_pair.wilcox_clt2$clt2base_6m.pval,"fdr")
res_species_pair.wilcox_clt2$clt2base.abundance=apply(species_rela_filter[row.names(species_rela_filter),PBCpair_6m_cluster2],1,mean)
res_species_pair.wilcox_clt2$clt2.6m.abundance=apply(species_rela_filter[row.names(species_rela_filter),PBC_6m_cluster2],1,mean)
res_species_pair.wilcox_clt2$change=res_species_pair.wilcox_clt2$clt2.6m.abundance-res_species_pair.wilcox_clt2$clt2base.abundance
res_species_pair.wilcox_clt2$clt2HC_coef=-maaslin_species_clt2HC[row.names(species_rela_filter),"coef"]
res_species_pair.wilcox_clt2$clt2HC_qval=maaslin_species_clt2HC[row.names(species_rela_filter),"qval"]
res_species_pair.wilcox_clt2$clt1clt2_coef=maaslin_species_clt1clt2[row.names(species_rela_filter),"coef"]
res_species_pair.wilcox_clt2$clt1clt2_qval=maaslin_species_clt1clt2[row.names(species_rela_filter),"qval"]
row.names(res_species_pair.wilcox_clt2)=row.names(species_rela_filter)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/paired_new")
write.csv(res_species_pair.wilcox_clt2,"res_species_pair.wilcox_clt2.csv")

#
#1.10 PBC_6m responder and non-responder
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin")
species_abs_clr_masslin=species_abs_clr[c(PBC_6m_noresp,PBC_6m_resp),row.names(species_rela_filter)]
colnames(species_abs_clr_masslin)=speciesID
write.csv(species_abs_clr_masslin,"species clr PBC6m resp noresp masslin.csv")
#metadata
metadata_PBC6m.resp.noresp=metadata[c(PBC_6m_noresp,PBC_6m_resp),c("age","gender","BMI","Response_6m")]
write.csv(metadata_PBC6m.resp.noresp,"metadata_PBC6m.resp.noresp.csv")
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/species clr PBC6m resp noresp masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata_PBC6m.resp.noresp.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/MAaSLIN2_species_PBC6m_resp.noresp.new',
                     fixed_effects = c("age","gender","BMI","Response_6m"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin")
write.csv(species_name,"species_name.csv")

#各组相对丰度
species_rela_mean$PBCbase_resp=apply(species_rela_filter[row.names(species_rela_filter),PBCbase_resp_6m],1,mean)
species_rela_mean$PBCbase_noresp=apply(species_rela_filter[row.names(species_rela_filter),PBCbase_noresp_6m],1,mean)
species_rela_mean$PBC6m_resp=apply(species_rela_filter[row.names(species_rela_filter),PBC_6m_resp],1,mean)
species_rela_mean$PBC6m_noresp=apply(species_rela_filter[row.names(species_rela_filter),PBC_6m_noresp],1,mean)
write.csv(species_rela_mean,"species_rela_mean.csv")
#1.11 PBCpair_6m_resp vs PBC_6m_resp

res_species_pair.wilcox_PBCresp=as.data.frame(matrix(0,ncol=1,nrow=3854))
for (i in 1:3854) {
  res_species_pair.wilcox_PBCresp$PBCbase_6m.resp.pval[i]=wilcox.test(species_abs_clr[PBCpair_6m_resp,row.names(species_rela_filter)][,i],
                                                               species_abs_clr[PBC_6m_resp,row.names(species_rela_filter)][,i],paired = T)$p.value
}
res_species_pair.wilcox_PBCresp$PBCbase_6m.resp.qval=p.adjust(res_species_pair.wilcox_PBCresp$PBCbase_6m.resp.pval,"fdr")
res_species_pair.wilcox_PBCresp$PBCbase.resp.abundance=apply(species_rela_filter[row.names(species_rela_filter),PBCpair_6m_resp],1,mean)
res_species_pair.wilcox_PBCresp$PBC.6m.resp.abundance=apply(species_rela_filter[row.names(species_rela_filter),PBC_6m_resp],1,mean)
res_species_pair.wilcox_PBCresp$change=res_species_pair.wilcox_PBCresp$PBC.6m.resp.abundance-res_species_pair.wilcox_PBCresp$PBCbase.resp.abundance
row.names(res_species_pair.wilcox_PBCresp)=row.names(species_rela_filter)

#1.12 PBCpair_6m_noresp vs PBC_6m_noresp

res_species_pair.wilcox_PBCnoresp=as.data.frame(matrix(0,ncol=1,nrow=3854))
for (i in 1:3854) {
  res_species_pair.wilcox_PBCnoresp$PBCbase_6m.noresp.pval[i]=wilcox.test(species_abs_clr[PBCpair_6m_noresp,row.names(species_rela_filter)][,i],
                                                                      species_abs_clr[PBC_6m_noresp,row.names(species_rela_filter)][,i],paired = T)$p.value
}
res_species_pair.wilcox_PBCnoresp$PBCbase_6m.noresp.qval=p.adjust(res_species_pair.wilcox_PBCnoresp$PBCbase_6m.noresp.pval,"fdr")
res_species_pair.wilcox_PBCnoresp$PBCbase.noresp.abundance=apply(species_rela_filter[row.names(species_rela_filter),PBCpair_6m_noresp],1,mean)
res_species_pair.wilcox_PBCnoresp$PBC.6m.noresp.abundance=apply(species_rela_filter[row.names(species_rela_filter),PBC_6m_noresp],1,mean)
res_species_pair.wilcox_PBCnoresp$change=res_species_pair.wilcox_PBCnoresp$PBC.6m.noresp.abundance-res_species_pair.wilcox_PBCnoresp$PBCbase.noresp.abundance
row.names(res_species_pair.wilcox_PBCnoresp)=row.names(species_rela_filter)

write.csv(res_species_pair.wilcox_PBCnoresp,"res_species_pair.wilcox_PBCnoresp.csv")
write.csv(res_species_pair.wilcox_PBCresp,"res_species_pair.wilcox_PBCresp.csv")

##################2.genus#################
#修改名字以与masslin分析的结果match,因为字符的问题;
genus_name=row.names(genus_rela_filter) %>% as.data.frame()
genusID=vector()
for (i in 1:883) {
  genusID[i]=paste("genus",row.names(genus_name)[i],sep = "")
}
row.names(genus_name)=genusID

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin")
#2.1PBCbase and HC_new
genus_abs_clr_masslin=genus_abs_clr[c(PBCbase,HC_new),row.names(genus_rela_filter)]
colnames(genus_abs_clr_masslin)=genusID
write.csv(genus_abs_clr_masslin,"genus clr PBCbase HC_new masslin.csv")
remove(genus_abs_clr_masslin)

library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin/genus clr PBCbase HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/data/metadata_PBCHC.new.maaslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin/MAaSLIN2_genus_PBC_HC.new',
                     fixed_effects = c("age","gender","BMI","group"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin")
write.csv(genus_name,"genus_name.csv")

#2.2 PBCbase responder and non-responder
genus_abs_clr_masslin=genus_abs_clr[c(PBCbase_resp_6m,PBCbase_noresp_6m),row.names(genus_rela_filter)]
colnames(genus_abs_clr_masslin)=genusID
write.csv(genus_abs_clr_masslin,"genus clr PBCbase resp noresp masslin.csv")
remove(genus_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin/genus clr PBCbase resp noresp masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove sc and 1y/metagenomics/data/genus/metadata PBCbase resp noresp masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin/MAaSLIN2_genus_PBCbase_resp_noresp',
                     fixed_effects = c("age","gender","BMI","Response_6m"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)
#2.3 scPBCbase and HC_all
genus_abs_clr_masslin=genus_abs_clr[c(scPBCbase,HC),row.names(genus_rela_filter)]
colnames(genus_abs_clr_masslin)=genusID
write.csv(genus_abs_clr_masslin,"genus clr scPBCbase HC masslin.csv")
remove(genus_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin/genus clr scPBCbase HC masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove 1y and keep sc/metagenomics/result/species/maaslin/metadata scPBCbase HC.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin/MAaSLIN2_genus_scPBCbase_HC',
                     fixed_effects = c("age","gender","BMI","group"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin")
#2.4 clt1 vs HC_new
genus_abs_clr_masslin=genus_abs_clr[c(PBCbase_cluster1,HC_new),row.names(genus_rela_filter)]
colnames(genus_abs_clr_masslin)=genusID
write.csv(genus_abs_clr_masslin,"genus clr clt1 HC_new masslin.csv")
remove(genus_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin/genus clr clt1 HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 HC_new masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin/MAaSLIN2_genus_clt1_HC.new',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
#2.5 clt2 vs HC_new
genus_abs_clr_masslin=genus_abs_clr[c(PBCbase_cluster2,HC_new),row.names(genus_rela_filter)]
colnames(genus_abs_clr_masslin)=genusID
write.csv(genus_abs_clr_masslin,"genus clr clt2 HC_new masslin.csv")
remove(genus_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin/genus clr clt2 HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt2 HC_new masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin/MAaSLIN2_genus_clt2_HC.new',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin")
#2.6 clt1 vs clt2
genus_abs_clr_masslin=genus_abs_clr[c(PBCbase_cluster1,PBCbase_cluster2),row.names(genus_rela_filter)]
colnames(genus_abs_clr_masslin)=genusID
write.csv(genus_abs_clr_masslin,"genus clr clt1 clt2 masslin.csv")
remove(genus_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin/genus clr clt1 clt2 masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 clt2 masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/genus/maaslin/MAaSLIN2_genus_clt1_clt2',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)

#2.7.1  PBC_pair vs PBC_6m_pair  #select genus
res_genus_pair.wilcox=as.data.frame(matrix(0,ncol=1,nrow=28))
for (i in 1:28) {
  res_genus_pair.wilcox$PBCbase_6m.pval[i]=wilcox.test(genus_abs_clr[PBC_6m_pair,maaslin_genus_PBCHC_dif][,i],
                                                         genus_abs_clr[PBCpair_6m,maaslin_genus_PBCHC_dif][,i],paired = T)$p.value
}
res_genus_pair.wilcox$PBCbase_6m.qval=p.adjust(res_genus_pair.wilcox$PBCbase_6m.pval,"fdr")
res_genus_pair.wilcox$PBCbase.abundance=apply(genus_rela_filter[maaslin_genus_PBCHC_dif,PBCpair_6m],1,mean)
res_genus_pair.wilcox$PBC6m.abundance=apply(genus_rela_filter[maaslin_genus_PBCHC_dif,PBC_6m_pair],1,mean)
res_genus_pair.wilcox$change=res_genus_pair.wilcox$PBC6m.abundance-res_genus_pair.wilcox$PBCbase.abundance
res_genus_pair.wilcox$PBCbase.HC_coef=maaslin_genus_PBCHC[maaslin_genus_PBCHC_dif,"coef"]
res_genus_pair.wilcox$PBCbase.HC_qval=maaslin_genus_PBCHC[maaslin_genus_PBCHC_dif,"qval"]
row.names(res_genus_pair.wilcox)=maaslin_genus_PBCHC_dif

#2.7.2 PBC_pair vs PBC_6m_pair  #all genus
res_genus_pair.wilcox=as.data.frame(matrix(0,ncol=1,nrow=883))
for (i in 1:883) {
  res_genus_pair.wilcox$PBCbase_6m.pval[i]=wilcox.test(genus_abs_clr[PBC_6m_pair,row.names(genus_rela_filter)][,i],
                                                       genus_abs_clr[PBCpair_6m,row.names(genus_rela_filter)][,i],paired = T)$p.value
}
res_genus_pair.wilcox$PBCbase_6m.qval=p.adjust(res_genus_pair.wilcox$PBCbase_6m.pval,"fdr")
res_genus_pair.wilcox$PBCbase.abundance=apply(genus_rela_filter[row.names(genus_rela_filter),PBCpair_6m],1,mean)
res_genus_pair.wilcox$PBC6m.abundance=apply(genus_rela_filter[row.names(genus_rela_filter),PBC_6m_pair],1,mean)
res_genus_pair.wilcox$PBC.abundance=apply(genus_rela_filter[row.names(genus_rela_filter),c(PBC_6m_pair,PBCpair_6m)],1,mean)

res_genus_pair.wilcox$change=res_genus_pair.wilcox$PBC6m.abundance-res_genus_pair.wilcox$PBCbase.abundance
res_genus_pair.wilcox$PBCbase.HC_coef=maaslin_genus_PBCHC[row.names(genus_rela_filter),"coef"]
res_genus_pair.wilcox$PBCbase.HC_pval=maaslin_genus_PBCHC[row.names(genus_rela_filter),"pval"]
res_genus_pair.wilcox$PBCbase.HC_qval=maaslin_genus_PBCHC[row.names(genus_rela_filter),"qval"]

row.names(res_genus_pair.wilcox)=row.names(genus_rela_filter)

#2.8.1 clt1_pair vs clt1_6m_pair  #select genus
res_genus_pair.wilcox_clt1=as.data.frame(matrix(0,ncol=1,nrow=663))
for (i in 1:663) {
  res_genus_pair.wilcox_clt1$clt1base_6m.pval[i]=wilcox.test(genus_abs_clr[PBCpair_6m_cluster1,maaslin_genus_clt1_pair][,i],
                                                               genus_abs_clr[PBC_6m_cluster1,maaslin_genus_clt1_pair][,i],paired = T)$p.value
}
res_genus_pair.wilcox_clt1$clt1base_6m.qval=p.adjust(res_genus_pair.wilcox_clt1$clt1base_6m.pval,"fdr")
res_genus_pair.wilcox_clt1$clt1base.abundance=apply(genus_rela_filter[maaslin_genus_clt1_pair,PBCpair_6m_cluster1],1,mean)
res_genus_pair.wilcox_clt1$clt1.6m.abundance=apply(genus_rela_filter[maaslin_genus_clt1_pair,PBC_6m_cluster1],1,mean)
res_genus_pair.wilcox_clt1$change=res_genus_pair.wilcox_clt1$clt1.6m.abundance-res_genus_pair.wilcox_clt1$clt1base.abundance
res_genus_pair.wilcox_clt1$clt1HC_coef=maaslin_genus_clt1HC[maaslin_genus_clt1_pair,"coef"]
res_genus_pair.wilcox_clt1$clt1HC_qval=maaslin_genus_clt1HC[maaslin_genus_clt1_pair,"qval"]
res_genus_pair.wilcox_clt1$clt1clt2_coef=maaslin_genus_clt1clt2[maaslin_genus_clt1_pair,"coef"]
res_genus_pair.wilcox_clt1$clt1clt2_qval=maaslin_genus_clt1clt2[maaslin_genus_clt1_pair,"qval"]
row.names(res_genus_pair.wilcox_clt1)=maaslin_genus_clt1_pair

#2.8.2 clt1_pair vs clt1_6m_pair  #all genus
res_genus_pair.wilcox_clt1=as.data.frame(matrix(0,ncol=1,nrow=883))
for (i in 1:883) {
  res_genus_pair.wilcox_clt1$clt1base_6m.pval[i]=wilcox.test(genus_abs_clr[PBCpair_6m_cluster1,row.names(genus_rela_filter)][,i],
                                                             genus_abs_clr[PBC_6m_cluster1,row.names(genus_rela_filter)][,i],paired = T)$p.value
}
res_genus_pair.wilcox_clt1$clt1base_6m.qval=p.adjust(res_genus_pair.wilcox_clt1$clt1base_6m.pval,"fdr")
res_genus_pair.wilcox_clt1$clt1base.abundance=apply(genus_rela_filter[row.names(genus_rela_filter),PBCpair_6m_cluster1],1,mean)
res_genus_pair.wilcox_clt1$clt1.6m.abundance=apply(genus_rela_filter[row.names(genus_rela_filter),PBC_6m_cluster1],1,mean)
res_genus_pair.wilcox_clt1$change=res_genus_pair.wilcox_clt1$clt1.6m.abundance-res_genus_pair.wilcox_clt1$clt1base.abundance
res_genus_pair.wilcox_clt1$clt1HC_coef=-maaslin_genus_clt1HC[row.names(genus_rela_filter),"coef"]
res_genus_pair.wilcox_clt1$clt1HC_qval=maaslin_genus_clt1HC[row.names(genus_rela_filter),"qval"]
res_genus_pair.wilcox_clt1$clt1clt2_coef=-maaslin_genus_clt1clt2[row.names(genus_rela_filter),"coef"]
res_genus_pair.wilcox_clt1$clt1clt2_qval=maaslin_genus_clt1clt2[row.names(genus_rela_filter),"qval"]
row.names(res_genus_pair.wilcox_clt1)=row.names(genus_rela_filter)

#2.9.1 clt2pair vs clt2_6m_pair #select genus
res_genus_pair.wilcox_clt2=as.data.frame(matrix(0,ncol=1,nrow=666))
for (i in 1:666) {
  res_genus_pair.wilcox_clt2$clt2base_6m.pval[i]=wilcox.test(genus_abs_clr[PBCpair_6m_cluster2,maaslin_genus_clt2_pair][,i],
                                                               genus_abs_clr[PBC_6m_cluster2,maaslin_genus_clt2_pair][,i],paired = T)$p.value
}
res_genus_pair.wilcox_clt2$clt2base_6m.qval=p.adjust(res_genus_pair.wilcox_clt2$clt2base_6m.pval,"fdr")
res_genus_pair.wilcox_clt2$clt2base.abundance=apply(genus_rela_filter[maaslin_genus_clt2_pair,PBCpair_6m_cluster2],1,mean)
res_genus_pair.wilcox_clt2$clt2.6m.abundance=apply(genus_rela_filter[maaslin_genus_clt2_pair,PBC_6m_cluster2],1,mean)
res_genus_pair.wilcox_clt2$change=res_genus_pair.wilcox_clt2$clt2.6m.abundance-res_genus_pair.wilcox_clt2$clt2base.abundance
res_genus_pair.wilcox_clt2$clt2HC_coef=maaslin_genus_clt2HC[maaslin_genus_clt2_pair,"coef"]
res_genus_pair.wilcox_clt2$clt2HC_qval=maaslin_genus_clt2HC[maaslin_genus_clt2_pair,"qval"]
res_genus_pair.wilcox_clt2$clt1clt2_coef=maaslin_genus_clt1clt2[maaslin_genus_clt2_pair,"coef"]
res_genus_pair.wilcox_clt2$clt1clt2_qval=maaslin_genus_clt1clt2[maaslin_genus_clt2_pair,"qval"]
row.names(res_genus_pair.wilcox_clt2)=maaslin_genus_clt2_pair

#2.9.2 clt2pair vs clt2_6m_pair #all genus
res_genus_pair.wilcox_clt2=as.data.frame(matrix(0,ncol=1,nrow=883))
for (i in 1:883) {
  res_genus_pair.wilcox_clt2$clt2base_6m.pval[i]=wilcox.test(genus_abs_clr[PBCpair_6m_cluster2,row.names(genus_rela_filter)][,i],
                                                             genus_abs_clr[PBC_6m_cluster2,row.names(genus_rela_filter)][,i],paired = T)$p.value
}
res_genus_pair.wilcox_clt2$clt2base_6m.qval=p.adjust(res_genus_pair.wilcox_clt2$clt2base_6m.pval,"fdr")
res_genus_pair.wilcox_clt2$clt2base.abundance=apply(genus_rela_filter[row.names(genus_rela_filter),PBCpair_6m_cluster2],1,mean)
res_genus_pair.wilcox_clt2$clt2.6m.abundance=apply(genus_rela_filter[row.names(genus_rela_filter),PBC_6m_cluster2],1,mean)
res_genus_pair.wilcox_clt2$change=res_genus_pair.wilcox_clt2$clt2.6m.abundance-res_genus_pair.wilcox_clt2$clt2base.abundance
res_genus_pair.wilcox_clt2$clt2HC_coef=-maaslin_genus_clt2HC[row.names(genus_rela_filter),"coef"]
res_genus_pair.wilcox_clt2$clt2HC_qval=maaslin_genus_clt2HC[row.names(genus_rela_filter),"qval"]
res_genus_pair.wilcox_clt2$clt1clt2_coef=maaslin_genus_clt1clt2[row.names(genus_rela_filter),"coef"]
res_genus_pair.wilcox_clt2$clt1clt2_qval=maaslin_genus_clt1clt2[row.names(genus_rela_filter),"qval"]
row.names(res_genus_pair.wilcox_clt2)=row.names(genus_rela_filter)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/genus/paired_new")
write.csv(res_genus_pair.wilcox,"res_genus_pair.wilcox_PBC.csv")
write.csv(res_genus_pair.wilcox_clt1,"res_genus_pair.wilcox_clt1.csv")
write.csv(res_genus_pair.wilcox_clt2,"res_genus_pair.wilcox_clt2.csv")

################3.family###################
#修改名字以与masslin分析的结果match,因为字符的问题;
family_name=row.names(family_rela_filter) %>% as.data.frame()
familyID=vector()
for (i in 1:333) {
  familyID[i]=paste("family",row.names(family_name)[i],sep = "")
}
row.names(family_name)=familyID

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/family/maaslin")
#3.1PBCbase and HC_new
family_abs_clr_masslin=family_abs_clr[c(PBCbase,HC_new),row.names(family_rela_filter)]
colnames(family_abs_clr_masslin)=familyID
write.csv(family_abs_clr_masslin,"family clr PBCbase HC_new masslin.csv")
remove(family_abs_clr_masslin)

library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/family/maaslin/family clr PBCbase HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/data/metadata_PBCHC.new.maaslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/family/maaslin/MAaSLIN2_family_PBC_HC.new',
                     fixed_effects = c("age","gender","BMI","group"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/family/maaslin")
write.csv(family_name,"family_name.csv")

#3.2 PBCbase responder and non-responder
family_abs_clr_masslin=family_abs_clr[c(PBCbase_resp_6m,PBCbase_noresp_6m),row.names(family_rela_filter)]
colnames(family_abs_clr_masslin)=familyID
write.csv(family_abs_clr_masslin,"family clr PBCbase resp noresp masslin.csv")
remove(family_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/family/maaslin/family clr PBCbase resp noresp masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove sc and 1y/metagenomics/data/genus/metadata PBCbase resp noresp masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/family/maaslin/MAaSLIN2_family_PBCbase_resp_noresp',
                     fixed_effects = c("age","gender","BMI","Response_6m"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)
#3.3 scPBCbase and HC_all
family_abs_clr_masslin=family_abs_clr[c(scPBCbase,HC),row.names(family_rela_filter)]
colnames(family_abs_clr_masslin)=familyID
write.csv(family_abs_clr_masslin,"family clr scPBCbase HC masslin.csv")
remove(family_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/family/maaslin/family clr scPBCbase HC masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove 1y and keep sc/metagenomics/result/species/maaslin/metadata scPBCbase HC.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/family/maaslin/MAaSLIN2_family_scPBCbase_HC',
                     fixed_effects = c("age","gender","BMI","group"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)
#3.4 clt1 vs HC_new
family_abs_clr_masslin=family_abs_clr[c(PBCbase_cluster1,HC_new),row.names(family_rela_filter)]
colnames(family_abs_clr_masslin)=familyID
write.csv(family_abs_clr_masslin,"family clr clt1 HC_new masslin.csv")
remove(family_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/family/maaslin/family clr clt1 HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 HC_new masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/family/maaslin/MAaSLIN2_family_clt1_HC.new',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
#3.5 clt2 vs HC_new
family_abs_clr_masslin=family_abs_clr[c(PBCbase_cluster2,HC_new),row.names(family_rela_filter)]
colnames(family_abs_clr_masslin)=familyID
write.csv(family_abs_clr_masslin,"family clr clt2 HC_new masslin.csv")
remove(family_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/family/maaslin/family clr clt2 HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt2 HC_new masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/family/maaslin/MAaSLIN2_family_clt2_HC.new',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/family/maaslin")
#3.6 clt1 vs clt2
family_abs_clr_masslin=family_abs_clr[c(PBCbase_cluster1,PBCbase_cluster2),row.names(family_rela_filter)]
colnames(family_abs_clr_masslin)=familyID
write.csv(family_abs_clr_masslin,"family clr clt1 clt2 masslin.csv")
remove(family_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/family/maaslin/family clr clt1 clt2 masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 clt2 masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/family/maaslin/MAaSLIN2_family_clt1_clt2',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)
################4.phylum##############
#修改名字以与masslin分析的结果match,因为字符的问题;
phylum_name=row.names(phylum_rela_filter) %>% as.data.frame()
phylumID=vector()
for (i in 1:69) {
  phylumID[i]=paste("phylum",row.names(phylum_name)[i],sep = "")
}
row.names(phylum_name)=phylumID

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/phylum/maaslin")
#4.1PBCbase and HC
phylum_abs_clr_masslin=phylum_abs_clr[c(PBCbase,HC),row.names(phylum_rela_filter)]
colnames(phylum_abs_clr_masslin)=phylumID
write.csv(phylum_abs_clr_masslin,"phylum clr PBCbase HC masslin.csv")
remove(phylum_abs_clr_masslin)

library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/phylum/maaslin/phylum clr PBCbase HC masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove sc and 1y/metagenomics/data/KEGG_path3/metadata PBC HC masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/phylum/maaslin/MAaSLIN2_phylum_PBC_HC',
                     fixed_effects = c("age","gender","BMI","group"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/phylum/maaslin")
write.csv(phylum_name,"phylum_name.csv")

#4.2 PBCbase responder and non-responder
phylum_abs_clr_masslin=phylum_abs_clr[c(PBCbase_resp_6m,PBCbase_noresp_6m),row.names(phylum_rela_filter)]
colnames(phylum_abs_clr_masslin)=phylumID
write.csv(phylum_abs_clr_masslin,"phylum clr PBCbase resp noresp masslin.csv")
remove(phylum_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/phylum/maaslin/phylum clr PBCbase resp noresp masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove sc and 1y/metagenomics/data/genus/metadata PBCbase resp noresp masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/phylum/maaslin/MAaSLIN2_phylum_PBCbase_resp_noresp',
                     fixed_effects = c("age","gender","BMI","Response_6m"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)
#4.3 scPBCbase and HC_all
phylum_abs_clr_masslin=phylum_abs_clr[c(scPBCbase,HC),row.names(phylum_rela_filter)]
colnames(phylum_abs_clr_masslin)=phylumID
write.csv(phylum_abs_clr_masslin,"phylum clr scPBCbase HC masslin.csv")
remove(phylum_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/phylum/maaslin/phylum clr scPBCbase HC masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove 1y and keep sc/metagenomics/result/species/maaslin/metadata scPBCbase HC.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/phylum/maaslin/MAaSLIN2_phylum_scPBCbase_HC',
                     fixed_effects = c("age","gender","BMI","group"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)
#4.4 clt1 vs HC
phylum_abs_clr_masslin=phylum_abs_clr[c(PBCbase_cluster1,HC),row.names(phylum_rela_filter)]
colnames(phylum_abs_clr_masslin)=phylumID
write.csv(phylum_abs_clr_masslin,"phylum clr clt1 HC masslin.csv")
remove(phylum_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/phylum/maaslin/phylum clr clt1 HC masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 HC masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/phylum/maaslin/MAaSLIN2_phylum_clt1_HC',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
#4.5 clt2 vs HC
phylum_abs_clr_masslin=phylum_abs_clr[c(PBCbase_cluster2,HC),row.names(phylum_rela_filter)]
colnames(phylum_abs_clr_masslin)=phylumID
write.csv(phylum_abs_clr_masslin,"phylum clr clt2 HC masslin.csv")
remove(phylum_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/phylum/maaslin/phylum clr clt2 HC masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt2 HC masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/phylum/maaslin/MAaSLIN2_phylum_clt2_HC',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/phylum/maaslin")
#4.6 clt1 vs clt2
phylum_abs_clr_masslin=phylum_abs_clr[c(PBCbase_cluster1,PBCbase_cluster2),row.names(phylum_rela_filter)]
colnames(phylum_abs_clr_masslin)=phylumID
write.csv(phylum_abs_clr_masslin,"phylum clr clt1 clt2 masslin.csv")
remove(phylum_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/phylum/maaslin/phylum clr clt1 clt2 masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 clt2 masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/phylum/maaslin/MAaSLIN2_phylum_clt1_clt2',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)

#############5.KEGG_path3###############
#修改名字以与masslin分析的结果match,因为字符的问题;
KEGG_path3_name=row.names(KEGG_path3_rela_filter) %>% as.data.frame()
KEGG_path3ID=vector()
for (i in 1:329) {
  KEGG_path3ID[i]=paste("KEGG_path3",row.names(KEGG_path3_name)[i],sep = "")
}
row.names(KEGG_path3_name)=KEGG_path3ID

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin")

#5.1 PBCbase and HC_new

KEGG_path3_abs_clr_masslin=KEGG_path3_abs_clr[c(PBCbase,HC_new),row.names(KEGG_path3_rela_filter)]
colnames(KEGG_path3_abs_clr_masslin)=KEGG_path3ID
write.csv(KEGG_path3_abs_clr_masslin,"KEGG_path3 clr PBCbase HC_new masslin.csv")
remove(KEGG_path3_abs_clr_masslin)

library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin/KEGG_path3 clr PBCbase HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/data/metadata_PBCHC.new.maaslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin/MAaSLIN2_KEGG_path3_PBC_HC.new',
                     fixed_effects = c("age","gender","BMI","group"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin")
write.csv(KEGG_path3_name,"KEGG_path3_name.csv")

#5.2 PBCbase responder and non-responder
KEGG_path3_abs_clr_masslin=KEGG_path3_abs_clr[c(PBCbase_resp_6m,PBCbase_noresp_6m),row.names(KEGG_path3_rela_filter)]
colnames(KEGG_path3_abs_clr_masslin)=KEGG_path3ID
write.csv(KEGG_path3_abs_clr_masslin,"KEGG_path3 clr PBCbase resp noresp masslin.csv")
remove(KEGG_path3_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin/KEGG_path3 clr PBCbase resp noresp masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove sc and 1y/metagenomics/data/genus/metadata PBCbase resp noresp masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin/MAaSLIN2_KEGG_path3_PBCbase_resp_noresp',
                     fixed_effects = c("age","gender","BMI","Response_6m"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)
#5.3 scPBCbase and HC_all
KEGG_path3_abs_clr_masslin=KEGG_path3_abs_clr[c(scPBCbase,HC),row.names(KEGG_path3_rela_filter)]
colnames(KEGG_path3_abs_clr_masslin)=KEGG_path3ID
write.csv(KEGG_path3_abs_clr_masslin,"KEGG_path3 clr scPBCbase HC masslin.csv")
remove(KEGG_path3_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin/KEGG_path3 clr scPBCbase HC masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove 1y and keep sc/metagenomics/result/species/maaslin/metadata scPBCbase HC.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin/MAaSLIN2_KEGG_path3_scPBCbase_HC',
                     fixed_effects = c("age","gender","BMI","group"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)

#5.4 clt1 vs HC_new
KEGG_path3_abs_clr_masslin=KEGG_path3_abs_clr[c(PBCbase_cluster1,HC_new),row.names(KEGG_path3_rela_filter)]
colnames(KEGG_path3_abs_clr_masslin)=KEGG_path3ID
write.csv(KEGG_path3_abs_clr_masslin,"KEGG_path3 clr clt1 HC_new masslin.csv")
remove(KEGG_path3_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin/KEGG_path3 clr clt1 HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 HC_new masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin/MAaSLIN2_KEGG_path3_clt1_HC.new',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
#5.5 clt2 vs HC_new
KEGG_path3_abs_clr_masslin=KEGG_path3_abs_clr[c(PBCbase_cluster2,HC_new),row.names(KEGG_path3_rela_filter)]
colnames(KEGG_path3_abs_clr_masslin)=KEGG_path3ID
write.csv(KEGG_path3_abs_clr_masslin,"KEGG_path3 clr clt2 HC_new masslin.csv")
remove(KEGG_path3_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin/KEGG_path3 clr clt2 HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt2 HC_new masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin/MAaSLIN2_KEGG_path3_clt2_HC.new',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin")
#5.6 clt1 vs clt2
KEGG_path3_abs_clr_masslin=KEGG_path3_abs_clr[c(PBCbase_cluster1,PBCbase_cluster2),row.names(KEGG_path3_rela_filter)]
colnames(KEGG_path3_abs_clr_masslin)=KEGG_path3ID
write.csv(KEGG_path3_abs_clr_masslin,"KEGG_path3 clr clt1 clt2 masslin.csv")
remove(KEGG_path3_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin/KEGG_path3 clr clt1 clt2 masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 clt2 masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/maaslin/MAaSLIN2_KEGG_path3_clt1_clt2',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)
#5.7.1 PBCpair_6m vs PBC_6m #select pathway
res_KEGG_path3_pair.wilcox=as.data.frame(matrix(0,ncol=1,nrow=52))
for (i in 1:52) {
  res_KEGG_path3_pair.wilcox$PBCbase_6m.pval[i]=wilcox.test(KEGG_path3_abs_clr[PBC_6m_pair,maaslin_KEGG_path3_PBCHC_dif][,i],
                                                         KEGG_path3_abs_clr[PBCpair_6m,maaslin_KEGG_path3_PBCHC_dif][,i],paired = T)$p.value
}
res_KEGG_path3_pair.wilcox$PBCbase_6m.qval=p.adjust(res_KEGG_path3_pair.wilcox$PBCbase_6m.pval,"fdr")
res_KEGG_path3_pair.wilcox$PBCbase.abundance=apply(KEGG_path3_rela_filter[maaslin_KEGG_path3_PBCHC_dif,PBCpair_6m],1,mean)
res_KEGG_path3_pair.wilcox$PBC6m.abundance=apply(KEGG_path3_rela_filter[maaslin_KEGG_path3_PBCHC_dif,PBC_6m_pair],1,mean)
res_KEGG_path3_pair.wilcox$change=res_KEGG_path3_pair.wilcox$PBC6m.abundance-res_KEGG_path3_pair.wilcox$PBCbase.abundance
row.names(res_KEGG_path3_pair.wilcox)=maaslin_KEGG_path3_PBCHC_dif

#5.7.2 PBCpair_6m vs PBC_6m  #all pathway
res_KEGG_path3_pair.wilcox=as.data.frame(matrix(0,ncol=1,nrow=329))
for (i in 1:329) {
  res_KEGG_path3_pair.wilcox$PBCbase_6m.pval[i]=wilcox.test(KEGG_path3_abs_clr[PBC_6m_pair,row.names(KEGG_path3_rela_filter)][,i],
                                                            KEGG_path3_abs_clr[PBCpair_6m,row.names(KEGG_path3_rela_filter)][,i],paired = T)$p.value
}
res_KEGG_path3_pair.wilcox$PBCbase_6m.qval=p.adjust(res_KEGG_path3_pair.wilcox$PBCbase_6m.pval,"fdr")
res_KEGG_path3_pair.wilcox$PBCbase.abundance=apply(KEGG_path3_rela_filter[row.names(KEGG_path3_rela_filter),PBCpair_6m],1,mean)
res_KEGG_path3_pair.wilcox$PBC6m.abundance=apply(KEGG_path3_rela_filter[row.names(KEGG_path3_rela_filter),PBC_6m_pair],1,mean)
res_KEGG_path3_pair.wilcox$change=res_KEGG_path3_pair.wilcox$PBC6m.abundance-res_KEGG_path3_pair.wilcox$PBCbase.abundance
res_KEGG_path3_pair.wilcox$PBCHC_coef=maaslin_KEGG_path3_PBCHC[row.names(KEGG_path3_rela_filter),"coef"]
res_KEGG_path3_pair.wilcox$PBCHC_pval=maaslin_KEGG_path3_PBCHC[row.names(KEGG_path3_rela_filter),"pval"]
res_KEGG_path3_pair.wilcox$PBCHC_qval=maaslin_KEGG_path3_PBCHC[row.names(KEGG_path3_rela_filter),"qval"]

row.names(res_KEGG_path3_pair.wilcox)=row.names(KEGG_path3_rela_filter)

#5.8.1 clt1_pair vs clt1_6m_pair #select pathway
res_KEGG_path3_pair.wilcox_clt1=as.data.frame(matrix(0,ncol=1,nrow=144))
for (i in 1:144) {
  res_KEGG_path3_pair.wilcox_clt1$clt1base_6m.pval[i]=wilcox.test(KEGG_path3_abs_clr[PBCpair_6m_cluster1,maaslin_KEGG_path3_clt1_pair][,i],
                                                             KEGG_path3_abs_clr[PBC_6m_cluster1,maaslin_KEGG_path3_clt1_pair][,i],paired = T)$p.value
}
res_KEGG_path3_pair.wilcox_clt1$clt1base_6m.qval=p.adjust(res_KEGG_path3_pair.wilcox_clt1$clt1base_6m.pval,"fdr")
res_KEGG_path3_pair.wilcox_clt1$clt1base.abundance=apply(KEGG_path3_rela_filter[maaslin_KEGG_path3_clt1_pair,PBCpair_6m_cluster1],1,mean)
res_KEGG_path3_pair.wilcox_clt1$clt1.6m.abundance=apply(KEGG_path3_rela_filter[maaslin_KEGG_path3_clt1_pair,PBC_6m_cluster1],1,mean)
res_KEGG_path3_pair.wilcox_clt1$change=res_KEGG_path3_pair.wilcox_clt1$clt1.6m.abundance-res_KEGG_path3_pair.wilcox_clt1$clt1base.abundance
res_KEGG_path3_pair.wilcox_clt1$clt1HC_coef=maaslin_KEGG_path3_clt1HC[maaslin_KEGG_path3_clt1_pair,"coef"]
res_KEGG_path3_pair.wilcox_clt1$clt1HC_qval=maaslin_KEGG_path3_clt1HC[maaslin_KEGG_path3_clt1_pair,"qval"]
res_KEGG_path3_pair.wilcox_clt1$clt1clt2_coef=maaslin_KEGG_path3_clt1clt2[maaslin_KEGG_path3_clt1_pair,"coef"]
res_KEGG_path3_pair.wilcox_clt1$clt1clt2_qval=maaslin_KEGG_path3_clt1clt2[maaslin_KEGG_path3_clt1_pair,"qval"]
row.names(res_KEGG_path3_pair.wilcox_clt1)=maaslin_KEGG_path3_clt1_pair

#5.8.2 clt1_pair vs clt1_6m_pair  #all pathway
res_KEGG_path3_pair.wilcox_clt1=as.data.frame(matrix(0,ncol=1,nrow=329))
for (i in 1:329) {
  res_KEGG_path3_pair.wilcox_clt1$clt1base_6m.pval[i]=wilcox.test(KEGG_path3_abs_clr[PBCpair_6m_cluster1,row.names(KEGG_path3_rela_filter)][,i],
                                                                  KEGG_path3_abs_clr[PBC_6m_cluster1,row.names(KEGG_path3_rela_filter)][,i],paired = T)$p.value
}
res_KEGG_path3_pair.wilcox_clt1$clt1base_6m.qval=p.adjust(res_KEGG_path3_pair.wilcox_clt1$clt1base_6m.pval,"fdr")
res_KEGG_path3_pair.wilcox_clt1$clt1base.abundance=apply(KEGG_path3_rela_filter[row.names(KEGG_path3_rela_filter),PBCpair_6m_cluster1],1,mean)
res_KEGG_path3_pair.wilcox_clt1$clt1.6m.abundance=apply(KEGG_path3_rela_filter[row.names(KEGG_path3_rela_filter),PBC_6m_cluster1],1,mean)
res_KEGG_path3_pair.wilcox_clt1$change=res_KEGG_path3_pair.wilcox_clt1$clt1.6m.abundance-res_KEGG_path3_pair.wilcox_clt1$clt1base.abundance
res_KEGG_path3_pair.wilcox_clt1$clt1HC_coef=-maaslin_KEGG_path3_clt1HC[row.names(KEGG_path3_rela_filter),"coef"]
res_KEGG_path3_pair.wilcox_clt1$clt1HC_qval=maaslin_KEGG_path3_clt1HC[row.names(KEGG_path3_rela_filter),"qval"]
res_KEGG_path3_pair.wilcox_clt1$clt1clt2_coef=-maaslin_KEGG_path3_clt1clt2[row.names(KEGG_path3_rela_filter),"coef"]
res_KEGG_path3_pair.wilcox_clt1$clt1clt2_qval=maaslin_KEGG_path3_clt1clt2[row.names(KEGG_path3_rela_filter),"qval"]
row.names(res_KEGG_path3_pair.wilcox_clt1)=row.names(KEGG_path3_rela_filter)

#5.9.1 clt2pair vs clt2_6m_pair  #select pathway
res_KEGG_path3_pair.wilcox_clt2=as.data.frame(matrix(0,ncol=1,nrow=135))
for (i in 1:135) {
  res_KEGG_path3_pair.wilcox_clt2$clt2base_6m.pval[i]=wilcox.test(KEGG_path3_abs_clr[PBCpair_6m_cluster2,maaslin_KEGG_path3_clt2_pair][,i],
                                                             KEGG_path3_abs_clr[PBC_6m_cluster2,maaslin_KEGG_path3_clt2_pair][,i],paired = T)$p.value
}
res_KEGG_path3_pair.wilcox_clt2$clt2base_6m.qval=p.adjust(res_KEGG_path3_pair.wilcox_clt2$clt2base_6m.pval,"fdr")
res_KEGG_path3_pair.wilcox_clt2$clt2base.abundance=apply(KEGG_path3_rela_filter[maaslin_KEGG_path3_clt2_pair,PBCpair_6m_cluster2],1,mean)
res_KEGG_path3_pair.wilcox_clt2$clt2.6m.abundance=apply(KEGG_path3_rela_filter[maaslin_KEGG_path3_clt2_pair,PBC_6m_cluster2],1,mean)
res_KEGG_path3_pair.wilcox_clt2$change=res_KEGG_path3_pair.wilcox_clt2$clt2.6m.abundance-res_KEGG_path3_pair.wilcox_clt2$clt2base.abundance
res_KEGG_path3_pair.wilcox_clt2$clt2HC_coef=maaslin_KEGG_path3_clt2HC[maaslin_KEGG_path3_clt2_pair,"coef"]
res_KEGG_path3_pair.wilcox_clt2$clt2HC_qval=maaslin_KEGG_path3_clt2HC[maaslin_KEGG_path3_clt2_pair,"qval"]
res_KEGG_path3_pair.wilcox_clt2$clt1clt2_coef=maaslin_KEGG_path3_clt1clt2[maaslin_KEGG_path3_clt2_pair,"coef"]
res_KEGG_path3_pair.wilcox_clt2$clt1clt2_qval=maaslin_KEGG_path3_clt1clt2[maaslin_KEGG_path3_clt2_pair,"qval"]
row.names(res_KEGG_path3_pair.wilcox_clt2)=maaslin_KEGG_path3_clt2_pair

#5.9.2 clt2pair vs clt2_6m_pair  #all pathway
res_KEGG_path3_pair.wilcox_clt2=as.data.frame(matrix(0,ncol=1,nrow=329))
for (i in 1:329) {
  res_KEGG_path3_pair.wilcox_clt2$clt2base_6m.pval[i]=wilcox.test(KEGG_path3_abs_clr[PBCpair_6m_cluster2,row.names(KEGG_path3_rela_filter)][,i],
                                                                  KEGG_path3_abs_clr[PBC_6m_cluster2,row.names(KEGG_path3_rela_filter)][,i],paired = T)$p.value
}
res_KEGG_path3_pair.wilcox_clt2$clt2base_6m.qval=p.adjust(res_KEGG_path3_pair.wilcox_clt2$clt2base_6m.pval,"fdr")
res_KEGG_path3_pair.wilcox_clt2$clt2base.abundance=apply(KEGG_path3_rela_filter[row.names(KEGG_path3_rela_filter),PBCpair_6m_cluster2],1,mean)
res_KEGG_path3_pair.wilcox_clt2$clt2.6m.abundance=apply(KEGG_path3_rela_filter[row.names(KEGG_path3_rela_filter),PBC_6m_cluster2],1,mean)
res_KEGG_path3_pair.wilcox_clt2$change=res_KEGG_path3_pair.wilcox_clt2$clt2.6m.abundance-res_KEGG_path3_pair.wilcox_clt2$clt2base.abundance
res_KEGG_path3_pair.wilcox_clt2$clt2HC_coef=-maaslin_KEGG_path3_clt2HC[row.names(KEGG_path3_rela_filter),"coef"]
res_KEGG_path3_pair.wilcox_clt2$clt2HC_qval=maaslin_KEGG_path3_clt2HC[row.names(KEGG_path3_rela_filter),"qval"]
res_KEGG_path3_pair.wilcox_clt2$clt1clt2_coef=maaslin_KEGG_path3_clt1clt2[row.names(KEGG_path3_rela_filter),"coef"]
res_KEGG_path3_pair.wilcox_clt2$clt1clt2_qval=maaslin_KEGG_path3_clt1clt2[row.names(KEGG_path3_rela_filter),"qval"]
row.names(res_KEGG_path3_pair.wilcox_clt2)=row.names(KEGG_path3_rela_filter)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KEGG_path3/paired_new")

write.csv(res_KEGG_path3_pair.wilcox,"res_KEGG_path3_pair.wilcox_PBC.csv")
write.csv(res_KEGG_path3_pair.wilcox_clt1,"res_KEGG_path3_pair.wilcox_clt1.csv")
write.csv(res_KEGG_path3_pair.wilcox_clt2,"res_KEGG_path3_pair.wilcox_clt2.csv")

#############6.KO_gene##############
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin")
#6.1 PBCbase and HC_new
KO_gene_abs_clr_masslin=KO_gene_abs_clr[c(PBCbase,HC_new),row.names(KO_gene_rela_filter)]
write.csv(KO_gene_abs_clr_masslin,"KO_gene clr PBCbase HC_new masslin.csv")
remove(KO_gene_abs_clr_masslin)

library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin/KO_gene clr PBCbase HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/data/metadata_PBCHC.new.maaslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin/MAaSLIN2_KO_gene_PBC_HC.new',
                     fixed_effects = c("age","gender","BMI","group"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin")
write.csv(gene_anno_filter,"gene_anno_filter.csv")

#6.2 PBCbase responder and non-responder
KO_gene_abs_clr_masslin=KO_gene_abs_clr[c(PBCbase_resp_6m,PBCbase_noresp_6m),row.names(KO_gene_rela_filter)]
write.csv(KO_gene_abs_clr_masslin,"KO_gene clr PBCbase resp noresp masslin.csv")
remove(KO_gene_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin/KO_gene clr PBCbase resp noresp masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove sc and 1y/metagenomics/data/genus/metadata PBCbase resp noresp masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin/MAaSLIN2_KO_gene_PBCbase_resp_noresp',
                     fixed_effects = c("age","gender","BMI","Response_6m"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)
#6.3 scPBCbase and HC_all
KO_gene_abs_clr_masslin=KO_gene_abs_clr[c(scPBCbase,HC),row.names(KO_gene_rela_filter)]
write.csv(KO_gene_abs_clr_masslin,"KO_gene clr scPBCbase HC masslin.csv")
remove(KO_gene_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin/KO_gene clr scPBCbase HC masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove 1y and keep sc/metagenomics/result/species/maaslin/metadata scPBCbase HC.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin/MAaSLIN2_KO_gene_scPBCbase_HC',
                     fixed_effects = c("age","gender","BMI","group"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)
#6.4 clt1 vs HC_new
KO_gene_abs_clr_masslin=KO_gene_abs_clr[c(PBCbase_cluster1,HC_new),row.names(KO_gene_rela_filter)]
write.csv(KO_gene_abs_clr_masslin,"KO_gene clr clt1 HC_new masslin.csv")
remove(KO_gene_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin/KO_gene clr clt1 HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 HC_new masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin/MAaSLIN2_KO_gene_clt1_HC.new',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
#6.5 clt2 vs HC_new
KO_gene_abs_clr_masslin=KO_gene_abs_clr[c(PBCbase_cluster2,HC_new),row.names(KO_gene_rela_filter)]
write.csv(KO_gene_abs_clr_masslin,"KO_gene clr clt2 HC_new masslin.csv")
remove(KO_gene_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin/KO_gene clr clt2 HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt2 HC_new masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin/MAaSLIN2_KO_gene_clt2_HC.new',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin")
#6.6 clt1 vs clt2
KO_gene_abs_clr_masslin=KO_gene_abs_clr[c(PBCbase_cluster1,PBCbase_cluster2),row.names(KO_gene_rela_filter)]
write.csv(KO_gene_abs_clr_masslin,"KO_gene clr clt1 clt2 masslin.csv")
remove(KO_gene_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin/KO_gene clr clt1 clt2 masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 clt2 masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/maaslin/MAaSLIN2_KO_gene_clt1_clt2',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)
#6.7.1 PBCpair_6m vs PBC_6m  #select gene
res_KO_gene_pair.wilcox=as.data.frame(matrix(0,ncol=1,nrow=1202))
for (i in 1:1202) {
  res_KO_gene_pair.wilcox$PBCbase_6m.pval[i]=wilcox.test(KO_gene_abs_clr[PBC_6m_pair,maaslin_KO_gene_PBCHC_dif][,i],
                                                            KO_gene_abs_clr[PBCpair_6m,maaslin_KO_gene_PBCHC_dif][,i],paired = T)$p.value
}
res_KO_gene_pair.wilcox$PBCbase_6m.qval=p.adjust(res_KO_gene_pair.wilcox$PBCbase_6m.pval,"fdr")
res_KO_gene_pair.wilcox$PBCbase.abs.abundance=apply(KO_gene_abs[maaslin_KO_gene_PBCHC_dif,PBCpair_6m],1,mean)
res_KO_gene_pair.wilcox$PBC6m.abs.abundance=apply(KO_gene_abs[maaslin_KO_gene_PBCHC_dif,PBC_6m_pair],1,mean)
res_KO_gene_pair.wilcox$change=res_KO_gene_pair.wilcox$PBC6m.abs.abundance-res_KO_gene_pair.wilcox$PBCbase.abs.abundance
row.names(res_KO_gene_pair.wilcox)=maaslin_KO_gene_PBCHC_dif

#6.7.2 PBCpair_6m vs PBC_6m #all KO gene
res_KO_gene_pair.wilcox=as.data.frame(matrix(0,ncol=1,nrow=9536))
for (i in 1:9536) {
  res_KO_gene_pair.wilcox$PBCbase_6m.pval[i]=wilcox.test(KO_gene_abs_clr[PBC_6m_pair,row.names(KO_gene_rela_filter)][,i],
                                                         KO_gene_abs_clr[PBCpair_6m,row.names(KO_gene_rela_filter)][,i],paired = T)$p.value
}
res_KO_gene_pair.wilcox$PBCbase_6m.qval=p.adjust(res_KO_gene_pair.wilcox$PBCbase_6m.pval,"fdr")
res_KO_gene_pair.wilcox$PBCbase.abs.abundance=apply(KO_gene_abs[row.names(KO_gene_rela_filter),PBCpair_6m],1,mean)
res_KO_gene_pair.wilcox$PBC6m.abs.abundance=apply(KO_gene_abs[row.names(KO_gene_rela_filter),PBC_6m_pair],1,mean)
res_KO_gene_pair.wilcox$change=res_KO_gene_pair.wilcox$PBC6m.abs.abundance-res_KO_gene_pair.wilcox$PBCbase.abs.abundance
res_KO_gene_pair.wilcox$PBCHC_coef=maaslin_KO_gene_PBCHC[row.names(KO_gene_rela_filter),"coef"]
res_KO_gene_pair.wilcox$PBCHC_pval=maaslin_KO_gene_PBCHC[row.names(KO_gene_rela_filter),"pval"]
res_KO_gene_pair.wilcox$PBCHC_qval=maaslin_KO_gene_PBCHC[row.names(KO_gene_rela_filter),"qval"]
row.names(res_KO_gene_pair.wilcox)=row.names(KO_gene_rela_filter)
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/paired_new")
write.csv(res_KO_gene_pair.wilcox,"res_KO_gene_pair.wilcox_PBC.csv")

#6.8.1 clt1_pair vs clt1_6m_pair #select genes
res_KO_gene_pair.wilcox_clt1=as.data.frame(matrix(0,ncol=1,nrow=3287))
for (i in 1:3287) {
  res_KO_gene_pair.wilcox_clt1$clt1base_6m.pval[i]=wilcox.test(KO_gene_abs_clr[PBCpair_6m_cluster1,maaslin_KO_gene_clt1_pair][,i],
                                                                  KO_gene_abs_clr[PBC_6m_cluster1,maaslin_KO_gene_clt1_pair][,i],paired = T)$p.value
}
res_KO_gene_pair.wilcox_clt1$clt1base_6m.qval=p.adjust(res_KO_gene_pair.wilcox_clt1$clt1base_6m.pval,"fdr")
res_KO_gene_pair.wilcox_clt1$clt1base.abundance=apply(KO_gene_abs[maaslin_KO_gene_clt1_pair,PBCpair_6m_cluster1],1,mean)
res_KO_gene_pair.wilcox_clt1$clt1.6m.abundance=apply(KO_gene_abs[maaslin_KO_gene_clt1_pair,PBC_6m_cluster1],1,mean)
res_KO_gene_pair.wilcox_clt1$change=res_KO_gene_pair.wilcox_clt1$clt1.6m.abundance-res_KO_gene_pair.wilcox_clt1$clt1base.abundance
res_KO_gene_pair.wilcox_clt1$clt1HC_coef=maaslin_KO_gene_clt1HC[maaslin_KO_gene_clt1_pair,"coef"]
res_KO_gene_pair.wilcox_clt1$clt1HC_qval=maaslin_KO_gene_clt1HC[maaslin_KO_gene_clt1_pair,"qval"]
res_KO_gene_pair.wilcox_clt1$clt1clt2_coef=maaslin_KO_gene_clt1clt2[maaslin_KO_gene_clt1_pair,"coef"]
res_KO_gene_pair.wilcox_clt1$clt1clt2_qval=maaslin_KO_gene_clt1clt2[maaslin_KO_gene_clt1_pair,"qval"]
row.names(res_KO_gene_pair.wilcox_clt1)=maaslin_KO_gene_clt1_pair

#6.8.2 clt1_pair vs clt1_6m_pair #all genes
res_KO_gene_pair.wilcox_clt1=as.data.frame(matrix(0,ncol=1,nrow=9536))
for (i in 1:9536) {
  res_KO_gene_pair.wilcox_clt1$clt1base_6m.pval[i]=wilcox.test(KO_gene_abs_clr[PBCpair_6m_cluster1,row.names(KO_gene_rela_filter)][,i],
                                                               KO_gene_abs_clr[PBC_6m_cluster1,row.names(KO_gene_rela_filter)][,i],paired = T)$p.value
}
res_KO_gene_pair.wilcox_clt1$clt1base_6m.qval=p.adjust(res_KO_gene_pair.wilcox_clt1$clt1base_6m.pval,"fdr")
res_KO_gene_pair.wilcox_clt1$clt1base.abundance=apply(KO_gene_abs[row.names(KO_gene_rela_filter),PBCpair_6m_cluster1],1,mean)
res_KO_gene_pair.wilcox_clt1$clt1.6m.abundance=apply(KO_gene_abs[row.names(KO_gene_rela_filter),PBC_6m_cluster1],1,mean)
res_KO_gene_pair.wilcox_clt1$change=res_KO_gene_pair.wilcox_clt1$clt1.6m.abundance-res_KO_gene_pair.wilcox_clt1$clt1base.abundance
res_KO_gene_pair.wilcox_clt1$clt1HC_coef=-maaslin_KO_gene_clt1HC[row.names(KO_gene_rela_filter),"coef"]
res_KO_gene_pair.wilcox_clt1$clt1HC_qval=maaslin_KO_gene_clt1HC[row.names(KO_gene_rela_filter),"qval"]
res_KO_gene_pair.wilcox_clt1$clt1clt2_coef=-maaslin_KO_gene_clt1clt2[row.names(KO_gene_rela_filter),"coef"]
res_KO_gene_pair.wilcox_clt1$clt1clt2_qval=maaslin_KO_gene_clt1clt2[row.names(KO_gene_rela_filter),"qval"]
row.names(res_KO_gene_pair.wilcox_clt1)=row.names(KO_gene_rela_filter)


#6.9.1 clt2pair vs clt2_6m_pair #select genes
res_KO_gene_pair.wilcox_clt2=as.data.frame(matrix(0,ncol=1,nrow=3705))
for (i in 1:3705) {
  res_KO_gene_pair.wilcox_clt2$clt2base_6m.pval[i]=wilcox.test(KO_gene_abs_clr[PBCpair_6m_cluster2,maaslin_KO_gene_clt2_pair][,i],
                                                                  KO_gene_abs_clr[PBC_6m_cluster2,maaslin_KO_gene_clt2_pair][,i],paired = T)$p.value
}
res_KO_gene_pair.wilcox_clt2$clt2base_6m.qval=p.adjust(res_KO_gene_pair.wilcox_clt2$clt2base_6m.pval,"fdr")
res_KO_gene_pair.wilcox_clt2$clt2base.abundance=apply(KO_gene_abs[maaslin_KO_gene_clt2_pair,PBCpair_6m_cluster2],1,mean)
res_KO_gene_pair.wilcox_clt2$clt2.6m.abundance=apply(KO_gene_abs[maaslin_KO_gene_clt2_pair,PBC_6m_cluster2],1,mean)
res_KO_gene_pair.wilcox_clt2$change=res_KO_gene_pair.wilcox_clt2$clt2.6m.abundance-res_KO_gene_pair.wilcox_clt2$clt2base.abundance
res_KO_gene_pair.wilcox_clt2$clt2HC_coef=maaslin_KO_gene_clt2HC[maaslin_KO_gene_clt2_pair,"coef"]
res_KO_gene_pair.wilcox_clt2$clt2HC_qval=maaslin_KO_gene_clt2HC[maaslin_KO_gene_clt2_pair,"qval"]
res_KO_gene_pair.wilcox_clt2$clt1clt2_coef=maaslin_KO_gene_clt1clt2[maaslin_KO_gene_clt2_pair,"coef"]
res_KO_gene_pair.wilcox_clt2$clt1clt2_qval=maaslin_KO_gene_clt1clt2[maaslin_KO_gene_clt2_pair,"qval"]
row.names(res_KO_gene_pair.wilcox_clt2)=maaslin_KO_gene_clt2_pair

#6.9.2 clt2pair vs clt2_6m_pair #all genes
res_KO_gene_pair.wilcox_clt2=as.data.frame(matrix(0,ncol=1,nrow=9536))
for (i in 1:9536) {
  res_KO_gene_pair.wilcox_clt2$clt2base_6m.pval[i]=wilcox.test(KO_gene_abs_clr[PBCpair_6m_cluster2,row.names(KO_gene_rela_filter)][,i],
                                                               KO_gene_abs_clr[PBC_6m_cluster2,row.names(KO_gene_rela_filter)][,i],paired = T)$p.value
}
res_KO_gene_pair.wilcox_clt2$clt2base_6m.qval=p.adjust(res_KO_gene_pair.wilcox_clt2$clt2base_6m.pval,"fdr")
res_KO_gene_pair.wilcox_clt2$clt2base.abundance=apply(KO_gene_abs[row.names(KO_gene_rela_filter),PBCpair_6m_cluster2],1,mean)
res_KO_gene_pair.wilcox_clt2$clt2.6m.abundance=apply(KO_gene_abs[row.names(KO_gene_rela_filter),PBC_6m_cluster2],1,mean)
res_KO_gene_pair.wilcox_clt2$change=res_KO_gene_pair.wilcox_clt2$clt2.6m.abundance-res_KO_gene_pair.wilcox_clt2$clt2base.abundance
res_KO_gene_pair.wilcox_clt2$clt2HC_coef=-maaslin_KO_gene_clt2HC[row.names(KO_gene_rela_filter),"coef"]
res_KO_gene_pair.wilcox_clt2$clt2HC_qval=maaslin_KO_gene_clt2HC[row.names(KO_gene_rela_filter),"qval"]
res_KO_gene_pair.wilcox_clt2$clt1clt2_coef=maaslin_KO_gene_clt1clt2[row.names(KO_gene_rela_filter),"coef"]
res_KO_gene_pair.wilcox_clt2$clt1clt2_qval=maaslin_KO_gene_clt1clt2[row.names(KO_gene_rela_filter),"qval"]
row.names(res_KO_gene_pair.wilcox_clt2)=row.names(KO_gene_rela_filter)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/KO_gene/paired_new")

write.csv(res_KO_gene_pair.wilcox_clt1,"res_KO_gene_pair.wilcox_clt1.csv")
write.csv(res_KO_gene_pair.wilcox_clt2,"res_KO_gene_pair.wilcox_clt2.csv")

################7.core species##################
species_name_core=row.names(species_rela_filter_core) %>% as.data.frame()
speciesID_core=vector()
for (i in 1:1814) {
  speciesID_core[i]=paste("species",row.names(species_name_core)[i],sep = "")
}
row.names(species_name_core)=speciesID_core

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin")
#7.1PBCbase and HC
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/core maaslin")
species_abs_clr_masslin_core=species_abs_clr[c(PBCbase,HC),row.names(species_rela_filter_core)]
colnames(species_abs_clr_masslin_core)=speciesID_core
write.csv(species_abs_clr_masslin_core,"species clr PBCbase HC masslin_core.csv")
remove(species_abs_clr_masslin_core)

library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/core maaslin/species clr PBCbase HC masslin_core.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove sc and 1y/metagenomics/data/KEGG_path3/metadata PBC HC masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/species/core maaslin/MAaSLIN2_species_PBC_HC',
                     fixed_effects = c("age","gender","BMI","group"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/core maaslin")
write.csv(species_name_core,"species_name_core.csv")
#7.2 PBCbase responder and non-responder
species_abs_clr_masslin_core=species_abs_clr[c(PBCbase_resp_6m,PBCbase_noresp_6m),row.names(species_rela_filter_core)]
colnames(species_abs_clr_masslin_core)=speciesID_core
write.csv(species_abs_clr_masslin_core,"species clr PBCbase resp noresp masslin_core.csv")
remove(species_abs_clr_masslin_core)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/core maaslin/species clr PBCbase resp noresp masslin_core.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove sc and 1y/metagenomics/data/genus/metadata PBCbase resp noresp masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/species/core maaslin/MAaSLIN2_species_PBCbase_resp_noresp',
                     fixed_effects = c("age","gender","BMI","Response_6m"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)
#7.3 scPBCbase and HC_all
species_abs_clr_masslin_core=species_abs_clr[c(scPBCbase,HC),row.names(species_rela_filter_core)]
colnames(species_abs_clr_masslin_core)=speciesID_core
write.csv(species_abs_clr_masslin_core,"species clr scPBCbase HC masslin_core.csv")
remove(species_abs_clr_masslin_core)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/core maaslin/species clr scPBCbase HC masslin_core.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove 1y and keep sc/metagenomics/result/species/maaslin/metadata scPBCbase HC.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/species/core maaslin/MAaSLIN2_species_scPBCbase_HC',
                     fixed_effects = c("age","gender","BMI","group"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)

remove(fit_data)

#7.4 clt1 vs HC
species_abs_clr_masslin_core=species_abs_clr[c(PBCbase_cluster1,HC),row.names(species_rela_filter_core)]
colnames(species_abs_clr_masslin_core)=speciesID_core
write.csv(species_abs_clr_masslin_core,"species clr clt1 HC masslin_core.csv")

remove(species_abs_clr_masslin_core)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/core maaslin/species clr clt1 HC masslin_core.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 HC masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/species/core maaslin/MAaSLIN2_species_clt1_HC',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
#7.5 clt2 vs HC
species_abs_clr_masslin_core=species_abs_clr[c(PBCbase_cluster2,HC),row.names(species_rela_filter_core)]
colnames(species_abs_clr_masslin_core)=speciesID_core
write.csv(species_abs_clr_masslin_core,"species clr clt2 HC masslin_core.csv")

remove(species_abs_clr_masslin_core)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/core maaslin/species clr clt2 HC masslin_core.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt2 HC masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/species/core maaslin/MAaSLIN2_species_clt2_HC',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)

#7.6 clt1 vs clt2
species_abs_clr_masslin_core=species_abs_clr[c(PBCbase_cluster1,PBCbase_cluster2),row.names(species_rela_filter_core)]
colnames(species_abs_clr_masslin_core)=speciesID_core
write.csv(species_abs_clr_masslin_core,"species clr clt1 clt2 masslin_core.csv")
remove(species_abs_clr_masslin_core)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/core maaslin/species clr clt1 clt2 masslin_core.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 clt2 masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/species/core maaslin/MAaSLIN2_species_clt1_clt2',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)

############8.metabolic module#################
#1.PBCHC_new
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin")
GMM_abs_clr_maaslin=GMM_abs_clr[c(PBCbase,HC_new),row.names(GMM_rela_filter)]
write.csv(GMM_abs_clr_maaslin,"GMM clr PBCbase HC masslin_new.csv")
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin/GMM clr PBCbase HC masslin_new.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/data/metadata_PBCHC.new.maaslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin/MAaSLIN2_GMM_PBC_HC_filter.new',
                     fixed_effects = c("age","gender","BMI","group"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)

#1.2 PBCbase responder and non-responder
GMM_abs_clr_masslin=GMM_abs_clr[c(PBCbase_resp_6m,PBCbase_noresp_6m),row.names(GMM_rela_filter)]
write.csv(GMM_abs_clr_masslin,"GMM clr PBCbase resp noresp masslin.csv")
remove(GMM_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin/GMM clr PBCbase resp noresp masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove sc and 1y/metagenomics/data/genus/metadata PBCbase resp noresp masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin/MAaSLIN2_GMM_PBCbase_resp_noresp',
                     fixed_effects = c("age","gender","BMI","Response_6m"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)
#1.3 scPBCbase and HC_all
GMM_abs_clr_masslin=GMM_abs_clr[c(scPBCbase,HC),row.names(GMM_rela_filter)]
write.csv(GMM_abs_clr_masslin,"GMM clr scPBCbase HC masslin.csv")
remove(GMM_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin/GMM clr scPBCbase HC masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove 1y and keep sc/metagenomics/result/species/maaslin/metadata scPBCbase HC.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin/MAaSLIN2_GMM_scPBCbase_HC',
                     fixed_effects = c("age","gender","BMI","group"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)

remove(fit_data)

#1.4 clt1 vs HC_new
GMM_abs_clr_masslin=GMM_abs_clr[c(PBCbase_cluster1,HC_new),row.names(GMM_rela_filter)]
write.csv(GMM_abs_clr_masslin,"GMM clr clt1 HC_new masslin.csv")
remove(GMM_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin/GMM clr clt1 HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 HC_new masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin/MAaSLIN2_GMM_clt1_HC.new',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
#1.5 clt2 vs HC_new
GMM_abs_clr_masslin=GMM_abs_clr[c(PBCbase_cluster2,HC_new),row.names(GMM_rela_filter)]

write.csv(GMM_abs_clr_masslin,"GMM clr clt2 HC_new masslin.csv")

remove(GMM_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin/GMM clr clt2 HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt2 HC_new masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin/MAaSLIN2_GMM_clt2_HC.new',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/GMM/maaslin")
#1.6 clt1 vs clt2
GMM_abs_clr_masslin=GMM_abs_clr[c(PBCbase_cluster1,PBCbase_cluster2),row.names(GMM_rela_filter)]

write.csv(GMM_abs_clr_masslin,"GMM clr clt1 clt2 masslin.csv")
remove(GMM_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin/GMM clr clt1 clt2 masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 clt2 masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/maaslin/MAaSLIN2_GMM_clt1_clt2',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)

#8.7.1 PBCpair_6m vs PBC_6m #select GMM
res_GMM_pair.wilcox=as.data.frame(matrix(0,ncol=1,nrow=33))
for (i in 1:33) {
  res_GMM_pair.wilcox$PBCbase_6m.pval[i]=wilcox.test(GMM_abs_clr[PBC_6m_pair,maaslin_GMM_PBCHC_dif][,i],
                                                         GMM_abs_clr[PBCpair_6m,maaslin_GMM_PBCHC_dif][,i],paired = T)$p.value
}
res_GMM_pair.wilcox$PBCbase_6m.qval=p.adjust(res_GMM_pair.wilcox$PBCbase_6m.pval,"fdr")
res_GMM_pair.wilcox$PBCbase.abs.abundance=apply(GMM_abs[maaslin_GMM_PBCHC_dif,PBCpair_6m],1,mean)
res_GMM_pair.wilcox$PBC6m.abs.abundance=apply(GMM_abs[maaslin_GMM_PBCHC_dif,PBC_6m_pair],1,mean)
res_GMM_pair.wilcox$change=res_GMM_pair.wilcox$PBC6m.abs.abundance-res_GMM_pair.wilcox$PBCbase.abs.abundance
res_GMM_pair.wilcox$PBCbase.HC_coef=maaslin_GMM_PBCHC[maaslin_GMM_PBCHC_dif,"coef"]
res_GMM_pair.wilcox$PBCbase.HC_qval=maaslin_GMM_PBCHC[maaslin_GMM_PBCHC_dif,"qval"]

row.names(res_GMM_pair.wilcox)=maaslin_GMM_PBCHC_dif

#8.7.2 PBCpair_6m vs PBC_6m  #all GMM
res_GMM_pair.wilcox=as.data.frame(matrix(0,ncol=1,nrow=129))
for (i in 1:129) {
  res_GMM_pair.wilcox$PBCbase_6m.pval[i]=wilcox.test(GMM_abs_clr[PBC_6m_pair,row.names(GMM_rela_filter)][,i],
                                                     GMM_abs_clr[PBCpair_6m,row.names(GMM_rela_filter)][,i],paired = T)$p.value
}
res_GMM_pair.wilcox$PBCbase_6m.qval=p.adjust(res_GMM_pair.wilcox$PBCbase_6m.pval,"fdr")
res_GMM_pair.wilcox$PBCbase.abs.abundance=apply(GMM_abs[row.names(GMM_rela_filter),PBCpair_6m],1,mean)
res_GMM_pair.wilcox$PBC6m.abs.abundance=apply(GMM_abs[row.names(GMM_rela_filter),PBC_6m_pair],1,mean)
res_GMM_pair.wilcox$change=res_GMM_pair.wilcox$PBC6m.abs.abundance-res_GMM_pair.wilcox$PBCbase.abs.abundance
res_GMM_pair.wilcox$PBCbase.HC_coef=maaslin_GMM_PBCHC[row.names(GMM_rela_filter),"coef"]
res_GMM_pair.wilcox$PBCbase.HC_pval=maaslin_GMM_PBCHC[row.names(GMM_rela_filter),"pval"]
res_GMM_pair.wilcox$PBCbase.HC_qval=maaslin_GMM_PBCHC[row.names(GMM_rela_filter),"qval"]

row.names(res_GMM_pair.wilcox)=row.names(GMM_rela_filter)

#8.8.1 clt1_pair vs clt1_6m_pair #select GMM
res_GMM_pair.wilcox_clt1=as.data.frame(matrix(0,ncol=1,nrow=65))
for (i in 1:65) {
  res_GMM_pair.wilcox_clt1$clt1base_6m.pval[i]=wilcox.test(GMM_abs_clr[PBCpair_6m_cluster1,maaslin_GMM_clt1_pair][,i],
                                                               GMM_abs_clr[PBC_6m_cluster1,maaslin_GMM_clt1_pair][,i],paired = T)$p.value
}
res_GMM_pair.wilcox_clt1$clt1base_6m.qval=p.adjust(res_GMM_pair.wilcox_clt1$clt1base_6m.pval,"fdr")
res_GMM_pair.wilcox_clt1$clt1base.abundance=apply(GMM_abs[maaslin_GMM_clt1_pair,PBCpair_6m_cluster1],1,mean)
res_GMM_pair.wilcox_clt1$clt1.6m.abundance=apply(GMM_abs[maaslin_GMM_clt1_pair,PBC_6m_cluster1],1,mean)
res_GMM_pair.wilcox_clt1$change=res_GMM_pair.wilcox_clt1$clt1.6m.abundance-res_GMM_pair.wilcox_clt1$clt1base.abundance
res_GMM_pair.wilcox_clt1$clt1HC_coef=maaslin_GMM_clt1HC[maaslin_GMM_clt1_pair,"coef"]
res_GMM_pair.wilcox_clt1$clt1HC_qval=maaslin_GMM_clt1HC[maaslin_GMM_clt1_pair,"qval"]
res_GMM_pair.wilcox_clt1$clt1clt2_coef=maaslin_GMM_clt1clt2[maaslin_GMM_clt1_pair,"coef"]
res_GMM_pair.wilcox_clt1$clt1clt2_qval=maaslin_GMM_clt1clt2[maaslin_GMM_clt1_pair,"qval"]
row.names(res_GMM_pair.wilcox_clt1)=maaslin_GMM_clt1_pair

#8.8.2 clt1_pair vs clt1_6m_pair #all GMM
res_GMM_pair.wilcox_clt1=as.data.frame(matrix(0,ncol=1,nrow=129))
for (i in 1:129) {
  res_GMM_pair.wilcox_clt1$clt1base_6m.pval[i]=wilcox.test(GMM_abs_clr[PBCpair_6m_cluster1,row.names(GMM_rela_filter)][,i],
                                                           GMM_abs_clr[PBC_6m_cluster1,row.names(GMM_rela_filter)][,i],paired = T)$p.value
}
res_GMM_pair.wilcox_clt1$clt1base_6m.qval=p.adjust(res_GMM_pair.wilcox_clt1$clt1base_6m.pval,"fdr")
res_GMM_pair.wilcox_clt1$clt1base.abundance=apply(GMM_abs[row.names(GMM_rela_filter),PBCpair_6m_cluster1],1,mean)
res_GMM_pair.wilcox_clt1$clt1.6m.abundance=apply(GMM_abs[row.names(GMM_rela_filter),PBC_6m_cluster1],1,mean)
res_GMM_pair.wilcox_clt1$change=res_GMM_pair.wilcox_clt1$clt1.6m.abundance-res_GMM_pair.wilcox_clt1$clt1base.abundance
res_GMM_pair.wilcox_clt1$clt1HC_coef=-maaslin_GMM_clt1HC[row.names(GMM_rela_filter),"coef"]
res_GMM_pair.wilcox_clt1$clt1HC_qval=maaslin_GMM_clt1HC[row.names(GMM_rela_filter),"qval"]
res_GMM_pair.wilcox_clt1$clt1clt2_coef=-maaslin_GMM_clt1clt2[row.names(GMM_rela_filter),"coef"]
res_GMM_pair.wilcox_clt1$clt1clt2_qval=maaslin_GMM_clt1clt2[row.names(GMM_rela_filter),"qval"]
row.names(res_GMM_pair.wilcox_clt1)=row.names(GMM_rela_filter)


#8.9.1 clt2pair vs clt2_6m_pair  #select GMM
res_GMM_pair.wilcox_clt2=as.data.frame(matrix(0,ncol=1,nrow=60))
for (i in 1:60) {
  res_GMM_pair.wilcox_clt2$clt2base_6m.pval[i]=wilcox.test(GMM_abs_clr[PBCpair_6m_cluster2,maaslin_GMM_clt2_pair][,i],
                                                               GMM_abs_clr[PBC_6m_cluster2,maaslin_GMM_clt2_pair][,i],paired = T)$p.value
}
res_GMM_pair.wilcox_clt2$clt2base_6m.qval=p.adjust(res_GMM_pair.wilcox_clt2$clt2base_6m.pval,"fdr")
res_GMM_pair.wilcox_clt2$clt2base.abundance=apply(GMM_abs[maaslin_GMM_clt2_pair,PBCpair_6m_cluster2],1,mean)
res_GMM_pair.wilcox_clt2$clt2.6m.abundance=apply(GMM_abs[maaslin_GMM_clt2_pair,PBC_6m_cluster2],1,mean)
res_GMM_pair.wilcox_clt2$change=res_GMM_pair.wilcox_clt2$clt2.6m.abundance-res_GMM_pair.wilcox_clt2$clt2base.abundance
res_GMM_pair.wilcox_clt2$clt2HC_coef=maaslin_GMM_clt2HC[maaslin_GMM_clt2_pair,"coef"]
res_GMM_pair.wilcox_clt2$clt2HC_qval=maaslin_GMM_clt2HC[maaslin_GMM_clt2_pair,"qval"]
res_GMM_pair.wilcox_clt2$clt1clt2_coef=maaslin_GMM_clt1clt2[maaslin_GMM_clt2_pair,"coef"]
res_GMM_pair.wilcox_clt2$clt1clt2_qval=maaslin_GMM_clt1clt2[maaslin_GMM_clt2_pair,"qval"]
row.names(res_GMM_pair.wilcox_clt2)=maaslin_GMM_clt2_pair

#8.9.2 clt2pair vs clt2_6m_pair  #all GMM
res_GMM_pair.wilcox_clt2=as.data.frame(matrix(0,ncol=1,nrow=129))
for (i in 1:129) {
  res_GMM_pair.wilcox_clt2$clt2base_6m.pval[i]=wilcox.test(GMM_abs_clr[PBCpair_6m_cluster2,row.names(GMM_rela_filter)][,i],
                                                           GMM_abs_clr[PBC_6m_cluster2,row.names(GMM_rela_filter)][,i],paired = T)$p.value
}
res_GMM_pair.wilcox_clt2$clt2base_6m.qval=p.adjust(res_GMM_pair.wilcox_clt2$clt2base_6m.pval,"fdr")
res_GMM_pair.wilcox_clt2$clt2base.abundance=apply(GMM_abs[row.names(GMM_rela_filter),PBCpair_6m_cluster2],1,mean)
res_GMM_pair.wilcox_clt2$clt2.6m.abundance=apply(GMM_abs[row.names(GMM_rela_filter),PBC_6m_cluster2],1,mean)
res_GMM_pair.wilcox_clt2$change=res_GMM_pair.wilcox_clt2$clt2.6m.abundance-res_GMM_pair.wilcox_clt2$clt2base.abundance
res_GMM_pair.wilcox_clt2$clt2HC_coef=-maaslin_GMM_clt2HC[row.names(GMM_rela_filter),"coef"]
res_GMM_pair.wilcox_clt2$clt2HC_qval=maaslin_GMM_clt2HC[row.names(GMM_rela_filter),"qval"]
res_GMM_pair.wilcox_clt2$clt1clt2_coef=maaslin_GMM_clt1clt2[row.names(GMM_rela_filter),"coef"]
res_GMM_pair.wilcox_clt2$clt1clt2_qval=maaslin_GMM_clt1clt2[row.names(GMM_rela_filter),"qval"]
row.names(res_GMM_pair.wilcox_clt2)=row.names(GMM_rela_filter)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/GMM_v2/paired_new")
write.csv(res_GMM_pair.wilcox,"res_GMM_pair.wilcox_PBC.csv")
write.csv(res_GMM_pair.wilcox_clt1,"res_GMM_pair.wilcox_clt1.csv")
write.csv(res_GMM_pair.wilcox_clt2,"res_GMM_pair.wilcox_clt2.csv")


##############9.class####################
class_name=row.names(class_rela_filter) %>% as.data.frame()
classID=vector()
for (i in 1:118) {
  classID[i]=paste("class",row.names(class_name)[i],sep = "")
}
row.names(class_name)=classID

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/class/maaslin")
#9.1PBCbase and HC_new
class_abs_clr_masslin=class_abs_clr[c(PBCbase,HC_new),row.names(class_rela_filter)]
colnames(class_abs_clr_masslin)=classID
write.csv(class_abs_clr_masslin,"class clr PBCbase HC_new masslin.csv")
remove(class_abs_clr_masslin)

library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/class/maaslin/class clr PBCbase HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/data/metadata_PBCHC.new.maaslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/class/maaslin/MAaSLIN2_class_PBC_HC.new',
                     fixed_effects = c("age","gender","BMI","group"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/class/maaslin")
write.csv(class_name,"class_name.csv")

#9.2 PBCbase responder and non-responder
class_abs_clr_masslin=class_abs_clr[c(PBCbase_resp_6m,PBCbase_noresp_6m),row.names(class_rela_filter)]
colnames(class_abs_clr_masslin)=classID
write.csv(class_abs_clr_masslin,"class clr PBCbase resp noresp masslin.csv")
remove(class_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/class/maaslin/class clr PBCbase resp noresp masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove sc and 1y/metagenomics/data/genus/metadata PBCbase resp noresp masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/class/maaslin/MAaSLIN2_class_PBCbase_resp_noresp',
                     fixed_effects = c("age","gender","BMI","Response_6m"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)
#9.3 scPBCbase and HC_all
class_abs_clr_masslin=class_abs_clr[c(scPBCbase,HC),row.names(class_rela_filter)]
colnames(class_abs_clr_masslin)=classID
write.csv(class_abs_clr_masslin,"class clr scPBCbase HC masslin.csv")
remove(class_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/class/maaslin/class clr scPBCbase HC masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/remove 1y and keep sc/metagenomics/result/species/maaslin/metadata scPBCbase HC.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/class/maaslin/MAaSLIN2_class_scPBCbase_HC',
                     fixed_effects = c("age","gender","BMI","group"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/class/maaslin")
#9.4 clt1 vs HC_new
class_abs_clr_masslin=class_abs_clr[c(PBCbase_cluster1,HC_new),row.names(class_rela_filter)]
colnames(class_abs_clr_masslin)=classID
write.csv(class_abs_clr_masslin,"class clr clt1 HC_new masslin.csv")
remove(class_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/class/maaslin/class clr clt1 HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 HC_new masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/class/maaslin/MAaSLIN2_class_clt1_HC.new',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
#9.5 clt2 vs HC_new
class_abs_clr_masslin=class_abs_clr[c(PBCbase_cluster2,HC_new),row.names(class_rela_filter)]
colnames(class_abs_clr_masslin)=classID
write.csv(class_abs_clr_masslin,"class clr clt2 HC_new masslin.csv")
remove(class_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/class/maaslin/class clr clt2 HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt2 HC_new masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/class/maaslin/MAaSLIN2_class_clt2_HC.new',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/class/maaslin")
#9.6 clt1 vs clt2
class_abs_clr_masslin=class_abs_clr[c(PBCbase_cluster1,PBCbase_cluster2),row.names(class_rela_filter)]
colnames(class_abs_clr_masslin)=classID
write.csv(class_abs_clr_masslin,"class clr clt1 clt2 masslin.csv")
remove(class_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/class/maaslin/class clr clt1 clt2 masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 clt2 masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/class/maaslin/MAaSLIN2_class_clt1_clt2',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)

##################10. order#######################
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/order/maaslin")
order_name=row.names(order_rela_filter) %>% as.data.frame()
orderID=vector()
for (i in 1:203) {
  orderID[i]=paste("order",row.names(order_name)[i],sep = "")
}
row.names(order_name)=orderID
#10.4 clt1 vs HC_new
order_abs_clr_masslin=order_abs_clr[c(PBCbase_cluster1,HC_new),row.names(order_rela_filter)]
colnames(order_abs_clr_masslin)=orderID
write.csv(order_abs_clr_masslin,"order clr clt1 HC_new masslin.csv")
remove(order_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/order/maaslin/order clr clt1 HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 HC_new masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/order/maaslin/MAaSLIN2_order_clt1_HC.new',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
#10.5 clt2 vs HC_new
order_abs_clr_masslin=order_abs_clr[c(PBCbase_cluster2,HC_new),row.names(order_rela_filter)]
colnames(order_abs_clr_masslin)=orderID
write.csv(order_abs_clr_masslin,"order clr clt2 HC_new masslin.csv")
remove(order_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/order/maaslin/order clr clt2 HC_new masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt2 HC_new masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/order/maaslin/MAaSLIN2_order_clt2_HC.new',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)

setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/order/maaslin")
#10.6 clt1 vs clt2
order_abs_clr_masslin=order_abs_clr[c(PBCbase_cluster1,PBCbase_cluster2),row.names(order_rela_filter)]
colnames(order_abs_clr_masslin)=orderID
write.csv(order_abs_clr_masslin,"order clr clt1 clt2 masslin.csv")
remove(order_abs_clr_masslin)
library(Maaslin2)
fit_data <- Maaslin2(input_data = 'D:/PBC multiomics data/analysis V6/metagenomics/result/order/maaslin/order clr clt1 clt2 masslin.txt'
                     ,input_metadata = 'D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin/metadata clt1 clt2 masslin.txt',
                     output='D:/PBC multiomics data/analysis V6/metagenomics/result/order/maaslin/MAaSLIN2_order_clt1_clt2',
                     fixed_effects = c("age","gender","BMI","cluster"),transform = 'NONE',normalization = 'NONE',analysis_method='LM',standardize=FALSE,max_significance = 0.1,plot_heatmap = FALSE,plot_scatter = FALSE)
remove(fit_data)
write.csv(order_name,"order_name.csv")
