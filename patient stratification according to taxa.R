
#workspace: D:/PBC multiomics data/analysis V6/metagenomics/code/metagenomics analysis.RData]
#=========================
#1.Assignment according to BO_ratio_family################
#Dominant taxa groups (Ba or Ru) were assigned from
#family-level relative abundances according to the comparative abundance of
#1.Bacteroidaceae and Oscillospiraceae
PBCbase_Bacteroidaceae=filter(data_PBC_Ait,BO_ratio_family>0) %>% row.names(.)
PBCbase_Oscillospiraceae=filter(data_PBC_Ait,BO_ratio_family<0) %>% row.names(.)
family=c(rep("Bacteroidaceae",55),rep("Oscillospiraceae",77)) %>% as.data.frame(.)
row.names(family)=c(PBCbase_Bacteroidaceae,PBCbase_Oscillospiraceae)
colnames(family)="family"
#2.比较两组的临床指标的差异
metadata_PBCbase=cbind(metadata_PBCbase,family[row.names(metadata_PBCbase),1])
colnames(metadata_PBCbase)[41]="family"

chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"family"]) #P=0.047
chisq.test(metadata_PBCbase[,"gp210"],metadata_PBCbase[,"family"]) #P=0.00947

#=========================
#2.Assignment according to BO_ratio_genus################
#Bacteroides and Oscillibacter
PBCbase_Bacteroides=filter(data_PBC_Ait,BO_ratio>0) %>% row.names(.)
PBCbase_Oscillibacter=filter(data_PBC_Ait,BO_ratio<0) %>% row.names(.)
PBCbase_genus=c(rep("Bacteroides",54),rep("Oscillibacter",78)) %>% as.data.frame(.)
row.names(PBCbase_genus)=c(PBCbase_Bacteroides,PBCbase_Oscillibacter)
colnames(PBCbase_genus)="PBCbase_genus"
#2.比较两组的临床指标的差异
metadata_PBCbase=cbind(metadata_PBCbase,PBCbase_genus[row.names(metadata_PBCbase),1])
colnames(metadata_PBCbase)[42]="PBCbase_genus"

chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"PBCbase_genus"]) #P=0.086
chisq.test(metadata_PBCbase[,"gp210"],metadata_PBCbase[,"PBCbase_genus"]) #P=0.017

quantile(data_PBC_Ait$BO_ratio_family)
filter(data_PBC_Ait,BO_ratio_family>1.2301644) %>% row.names(.)
#=========================
#4.Assignment according to Oscillospiraceae abundance################
quantile(data_PBC_Ait$Oscillospiraceae,probs = c(0.1,0.25,0.5,0.7,0.75,0.9))

PBCbase_Oscillospiraceae_high=filter(data_PBC_Ait,Oscillospiraceae>0.9775213  ) %>% row.names(.)
PBCbase_Oscillospiraceae_low=filter(data_PBC_Ait,Oscillospiraceae<0.9775213  ) %>% row.names(.)
PBCbase_Oscillospiraceae_abundance=c(rep("Oscillospiraceae_high",40),rep("Oscillospiraceae_low",92)) %>% as.data.frame(.)
row.names(PBCbase_Oscillospiraceae_abundance)=c(PBCbase_Oscillospiraceae_high,PBCbase_Oscillospiraceae_low)
colnames(PBCbase_Oscillospiraceae_abundance)="Oscillospiraceae_abundance"
#2.比较两组的临床指标的差异
metadata_PBCbase=cbind(metadata_PBCbase[,1:42],PBCbase_Oscillospiraceae_abundance[row.names(metadata_PBCbase),1])
colnames(metadata_PBCbase)[43]="Oscillospiraceae_abundance"
colnames(metadata_PBCbase)[45]="Oscillibacter_abundance_0.75"
chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"Oscillospiraceae_abundance"]) #P=0.085
#50% percentile: P=0.14 66vs66
#70% percentile: P=0.09 92vs40 
#75% percentile: P=0.057 99vs33

#=========================
#5.Assignment according to Oscillospiracter abundance################
quantile(data_PBC_Ait$Oscillibacter,probs = c(0.1,0.25,0.5,0.7,0.75,0.9))

PBCbase_Oscillibacter_high=filter(data_PBC_Ait,Oscillibacter> 1.0393114 ) %>% row.names(.)
PBCbase_Oscillibacter_low=filter(data_PBC_Ait,Oscillibacter< 1.0393114 ) %>% row.names(.)
PBCbase_Oscillibacter_abundance=c(rep("Oscillibacter_high",40),rep("Oscillibacter_low",92)) %>% as.data.frame(.)
row.names(PBCbase_Oscillibacter_abundance)=c(PBCbase_Oscillibacter_high,PBCbase_Oscillibacter_low)
colnames(PBCbase_Oscillibacter_abundance)="Oscillibacter_abundance"
#2.比较两组的临床指标的差异
metadata_PBCbase=cbind(metadata_PBCbase[,1:45],PBCbase_Oscillibacter_abundance[row.names(metadata_PBCbase),1])
colnames(metadata_PBCbase)[46]="Oscillibacter_abundance_0.7"

chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"Oscillibacter_abundance_0.75"]) 
#50% percentile: P=0.14 66vs66
#70% percentile: P=0.09 92vs40 
#75% percentile: P=0.057 99vs33

#=================
#6.Assignment according to BO_ratio_family value################
#calculate quantile
quantile(data_PBC_Ait$BO_ratio_family)

PBCbase_BO_ratio_family_high=filter(data_PBC_Ait,BO_ratio_family>0) %>% row.names(.) #n=33
PBCbase_BO_ratio_family_low=filter(data_PBC_Ait,BO_ratio_family< 0 ) %>% row.names(.) #n=99
PBCbase_BO_ratio_family_abundance=c(rep("PBCbase_BO_ratio_family_high",55),rep("Oscillospiraceae_low",77)) %>% as.data.frame(.)
row.names(PBCbase_BO_ratio_family_abundance)=c(PBCbase_BO_ratio_family_high,PBCbase_BO_ratio_family_low)
colnames(PBCbase_BO_ratio_family_abundance)="BO_ratio_family_abundance"
#2.比较两组的临床指标的差异
metadata_PBCbase=cbind(metadata_PBCbase[,1:43],PBCbase_BO_ratio_family_abundance[row.names(metadata_PBCbase),1])
colnames(metadata_PBCbase)[44]="BO_ratio_family_abundance"

chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"BO_ratio_family_abundance"]) #P=0.085

#======================
#7.Assignment according to Anaerobutyricum hallii################
data_PBC_Ait$Anaerobutyricum_hallii=t(species_rela_filter)[PBCbase,"k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Anaerobutyricum;s__Anaerobutyricum hallii"]
quantile(data_PBC_Ait$Anaerobutyricum_hallii)
PBCbase_Anaerobutyricum_hallii_high=filter(data_PBC_Ait,Anaerobutyricum_hallii>0.001) %>% row.names(.)
PBCbase_Anaerobutyricum_hallii_low=filter(data_PBC_Ait,Anaerobutyricum_hallii<0.001) %>% row.names(.)
PBCbase_Anaerobutyricum_hallii=c(rep("PBCbase_Anaerobutyricum_hallii_high",36),rep("PBCbase_Anaerobutyricum_hallii_low",96)) %>% as.data.frame(.)
row.names(PBCbase_Anaerobutyricum_hallii)=c(PBCbase_Anaerobutyricum_hallii_high,PBCbase_Anaerobutyricum_hallii_low)
colnames(PBCbase_Anaerobutyricum_hallii)="Anaerobutyricum_hallii_relab"

#2.比较两组的临床指标的差异

metadata_PBCbase=cbind(metadata_PBCbase[,1:46],PBCbase_Anaerobutyricum_hallii[row.names(metadata_PBCbase),1])
colnames(metadata_PBCbase)[47]="Anaerobutyricum_hallii_0.001"

chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"Anaerobutyricum_hallii_0.001"])
chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"Anaerobutyricum_hallii_0.001"])$observed

#======================
#8.Assignment according to [Clostridium]_leptum################
data_PBC_Ait$Clostridium_leptum=t(species_rela_filter)[PBCbase,"k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__unidentified;s__[Clostridium] leptum"]
quantile(data_PBC_Ait$Clostridium_leptum)
PBCbase_Clostridium_leptum_high=filter(data_PBC_Ait,Clostridium_leptum>0.0001) %>% row.names(.)
PBCbase_Clostridium_leptum_low=filter(data_PBC_Ait,Clostridium_leptum<0.0001) %>% row.names(.)
PBCbase_Clostridium_leptum=c(rep("PBCbase_Clostridium_leptum_high",67),rep("PBCbase_Clostridium_leptum_low",65)) %>% as.data.frame(.)
row.names(PBCbase_Clostridium_leptum)=c(PBCbase_Clostridium_leptum_high,PBCbase_Clostridium_leptum_low)
colnames(PBCbase_Clostridium_leptum)="Clostridium_leptum_relab"

#2.比较两组的临床指标的差异

metadata_PBCbase=cbind(metadata_PBCbase[,1:49],PBCbase_Clostridium_leptum[row.names(metadata_PBCbase),1])
colnames(metadata_PBCbase)[50]="Clostridium_leptum_0.0001"

chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"Clostridium_leptum_0.00015"]) #P=0.0235
chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"Clostridium_leptum_0.00015"])$observed

#======================
#9.Assignment according to [Clostridium]_scindens################
data_PBC_Ait$Clostridium_scindens=t(species_rela_filter)[PBCbase,"k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Lachnoclostridium;s__[Clostridium] scindens"]
quantile(data_PBC_Ait$Clostridium_scindens)
PBCbase_Clostridium_scindens_high=filter(data_PBC_Ait,Clostridium_scindens>0.00005) %>% row.names(.)
PBCbase_Clostridium_scindens_low=filter(data_PBC_Ait,Clostridium_scindens<0.00005) %>% row.names(.)
PBCbase_Clostridium_scindens=c(rep("PBCbase_Clostridium_scindens_high",58),rep("PBCbase_Clostridium_scindens_low",74)) %>% as.data.frame(.)
row.names(PBCbase_Clostridium_scindens)=c(PBCbase_Clostridium_scindens_high,PBCbase_Clostridium_scindens_low)
colnames(PBCbase_Clostridium_scindens)="Clostridium_scindens_relab"

#2.比较两组的临床指标的差异

metadata_PBCbase=cbind(metadata_PBCbase[,1:48],PBCbase_Clostridium_scindens[row.names(metadata_PBCbase),1])
colnames(metadata_PBCbase)[49]="Clostridium_scindens_0.00005"

chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"Clostridium_scindens_0.00005"]) #P=0.035
chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"Clostridium_scindens_0.00005"])$observed
