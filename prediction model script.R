#==========================
#预测数据###############
#basic information
HC_no=c("MXH063","MXH1022","MXH1076","MXH1136","MXH444","MXH448","MXH459",
        "MXH496","MXH497","MXH532","MXH565","MXH779","MXH827","MXH840","MXH966","MXH976",
        "MXH455","MXH605","MXH572")
index=HC %in% HC_no #为逻辑值向量
index=which(index==TRUE)#为位置索引
HC_new=HC[-index]
#===================================
##读入species maaslin UDCA response vs. noresponse分析结果
setwd("D:/PBC multiomics data/analysis V6/metagenomics/result/species/maaslin")
maaslin_species_PBCresp=read.csv("maaslin_species PBCresp.csv",header = T,row.names=1)
row.names(maaslin_species_PBCresp)=maaslin_species_PBCresp$species
library(dplyr)
maaslin_species_PBCresp_dif=filter(maaslin_species_PBCresp,pval<0.01) %>% row.names(.)
species_resp_pred=species_abs_clr[c(PBCbase_resp_6m,PBCbase_noresp_6m),maaslin_species_PBCresp_dif]
species_resp_pred=cbind(species_resp_pred,metadata[c(PBCbase_resp_6m,PBCbase_noresp_6m),"Response_6m"])
colnames(species_resp_pred)[164]="response"
#====================================
#血清代谢物的预测数据
PBCbase_serum_resp=intersect(PBCbase_serum,PBCbase_resp_6m)
PBCbase_serum_noresp=intersect(PBCbase_serum,PBCbase_noresp_6m)
serum_PBCresp_lm_df=filter(serum_res,pval_PBCresp.noresp<0.05)%>% row.names(.)
serum_resp_pred=serum_0.7_log2UV[c(PBCbase_serum_resp,PBCbase_serum_noresp),serum_PBCresp_lm_df]
serum_resp_pred=cbind(serum_resp_pred,metadata[c(PBCbase_serum_resp,PBCbase_serum_noresp),"Response_6m"])
colnames(serum_resp_pred)[33]="response"
#====================================
#粪便代谢物的预测数据
PBCbase_fecal_resp=intersect(PBCbase_fecal,PBCbase_resp_6m)
PBCbase_fecal_noresp=intersect(PBCbase_fecal,PBCbase_noresp_6m)
fecal_PBCresp_lm_df=filter(fecal_res,pval_PBCresp.noresp<0.05)%>% row.names(.)
fecal_resp_pred=fecal_0.7_log2UV[c(PBCbase_fecal_resp,PBCbase_fecal_noresp),fecal_PBCresp_lm_df]
fecal_resp_pred=cbind(fecal_resp_pred,metadata[c(PBCbase_fecal_resp,PBCbase_fecal_noresp),"Response_6m"])
colnames(fecal_resp_pred)[48]="response"
#====================================
#联合组学数据进行预测
PBCbase_fecal.serum_resp=intersect(PBCbase_fecal.serum,PBCbase_resp_6m)
PBCbase_fecal.serum_noresp=intersect(PBCbase_fecal.serum,PBCbase_noresp_6m)
combine.sample=c(PBCbase_fecal.serum_resp,PBCbase_fecal.serum_noresp)
combine_resp_pred=cbind(species_resp_pred[combine.sample,1:163],serum_resp_pred[combine.sample,1:32],fecal_resp_pred[combine.sample,1:48])
#=================================
#临床指标的预测数据

clinic_resp_pred=metadata_PBCbase[c(PBCbase_resp_6m,PBCbase_noresp_6m),c(clini,"Response_6m")]

#=====================
#指定genus





############################1.response/no-response prediction############################
###########1.1 random forest analysis###########
###############1.1.1 confusionMatrix method####################
##1.1.1.1 species
library(dplyr)
library(randomForest)
library(pROC)
names(species_resp_pred) <- make.names(names(species_resp_pred))
colnames(species_resp_pred)=c(maaslin_species_PBCresp_dif,"response")
set.seed(1)
#random forest model,为了下面选择importance或者gini排名最靠前的feature
resp_species_rf=randomForest(response~.,data = species_resp_pred,importance=T,proximity=T)
resp_species_rf_votes=as.data.frame(resp_species_rf$votes)
roc_resp_species_rf <- roc(species_resp_pred$response,resp_species_rf_votes$"1")
plot(roc_resp_species_rf, print.thres=TRUE,print.auc=T)
set.seed(1)
#determine the best number of included features
result <- replicate(10, rfcv(species_resp_pred[,1:163], as.factor(species_resp_pred$response)), simplify=FALSE)
error.cv <- sapply(result, "[[", "error.cv")
matplot(result[[1]]$n.var, cbind(rowMeans(error.cv), error.cv), type="l", lwd=c(2, rep(1, ncol(error.cv))), col=1, lty=1, log="x", xlab="Number of variables", ylab="CV Error")
#看每个feature对于模型预测的重要性
resp_species_rf_importance=resp_species_rf$importance %>% as.data.frame(.)
colnames(resp_species_rf_importance)
#修改species名字，因为字符的问题
row.names(resp_species_rf_importance)=maaslin_species_PBCresp_dif

#按照gini指数进行降序排列
resp_species_rf_importance=arrange(resp_species_rf_importance,desc(MeanDecreaseGini))
#排名前41位的species的名称
resp_species_rf_mtry41=row.names(resp_species_rf_importance)[1:41]
#生成新的用于预测的数据，只有这41个species和response信息
species_resp_pred_mtry41=cbind(species_resp_pred[,resp_species_rf_mtry41],species_resp_pred[,164])
colnames(species_resp_pred_mtry41)[42]="response"
#进行模型的预测
names(species_resp_pred_mtry41)=make.names(names(species_resp_pred_mtry41))
resp_species_rf=randomForest(response~.,data = species_resp_pred_mtry41,importance=T,proximity=T,mtry=41)
resp_species_rf_votes=as.data.frame(resp_species_rf$votes)
roc_resp_species_rf <- roc(species_resp_pred_mtry41$response,resp_species_rf_votes$"1")
plot(roc_resp_species_rf, print.thres=TRUE,print.auc=T)
#输出用于建模的species
write.csv(species_resp_pred_mtry41,"species_resp_pred_mtry41.csv")

##1.1.1.2.serum metabolite
colnames(serum_res)
serum_resp.noresp_dif_wilcox=filter(serum_res,wilcox_q_base_resp.noresp<0.2) %>% row.names(.) #11个代谢物
#生成用于预测模型的数据
serum_resp_pred=serum_0.7_log2UV[PBCbase_serum,serum_resp.noresp_dif_wilcox]
serum_resp_pred=cbind(serum_resp_pred,metadata[PBCbase_serum,"Response_6m"])
colnames(serum_resp_pred)[12]="response"
serum_resp_pred=serum_resp_pred[-16,]
names(serum_resp_pred) <- make.names(names(serum_resp_pred))
str(serum_resp_pred)
#进行模型的预测
resp_serum_rf=randomForest(response~.,data = serum_resp_pred,importance=T,proximity=T,mtry=11)
resp_serum_rf_votes=as.data.frame(resp_serum_rf$votes)
roc_resp_serum_rf <- roc(serum_resp_pred$response,resp_serum_rf_votes$"1")
plot(roc_resp_serum_rf, print.thres=TRUE,print.auc=T)

#1.1.1.3.血清代谢物联合species
species_serum_resp_pred=cbind(species_resp_pred_mtry29[row.names(serum_resp_pred),-30],serum_resp_pred)
resp_species_serum_rf=randomForest(response~.,data = species_serum_resp_pred,importance=T,proximity=T,mtry=40)
resp_species_serum_rf_votes=as.data.frame(resp_species_serum_rf$votes)
roc_resp_species_serum_rf <- roc(species_serum_resp_pred$response,resp_species_serum_rf_votes$"1")
plot(roc_resp_species_serum_rf, print.thres=TRUE,print.auc=T)

#1.1.1.4.fecal metabolite(未进行造模，因为fecal metabolite无FDR<0.2的代谢物)
colnames(fecal_res)
fecal_resp.noresp_dif_wilcox=filter(fecal_res,wilcox_q_resp.noresp<0.2) %>% row.names(.) #character(0)

###############1.1.2random forest five-fold cross-validation method####################
library(dplyr)
library(randomForest)
library(pROC)
#species ===============================================================
names(species_resp_pred) <- make.names(names(species_resp_pred))
colnames(species_resp_pred)=c(maaslin_species_PBCresp_dif,"response")
set.seed(1)

randomized_samples = sample(1:nrow(species_resp_pred))
batches = split(randomized_samples,cut(seq_along(randomized_samples),5,labels = F))


train_X = lapply(batches, function(x) species_resp_pred[-x,-164])
test_X = lapply(batches, function(x) species_resp_pred[x,-164])

train_Y = lapply(batches, function(x) species_resp_pred[-x,164])
test_Y = lapply(batches, function(x) species_resp_pred[x,164])

## Build prediction models

#Cross Validation with Manual Fine Tuning
sqtmtry<- round(sqrt(ncol(species_resp_pred) - 1))
rfGrid <- expand.grid(mtry = c(6:26))

all_prediction_models = list()
print(' >> training models')
for(i in 1:5) {
  print(paste0('X-val set ',i,' '))
  
  current_X = train_X[[i]]
  current_Y = train_Y[[i]]
  
  all_prediction_models[[paste0("pred.batch",i)]] = train(x = as.matrix(current_X),
                                                              y = current_Y,
                                                              method = "rf", ntree = 500,tuneGrid = rfGrid,preProcess = c("scale"))
}
# predictions on train set
train.predictions.microb = list()
for (i in 1:5){
  train.predictions.microb[[i]] = 
     predict(all_prediction_models[[paste0("pred.batch",i)]]$finalModel,
            newx = as.matrix(train_X[[i]]))
}
# predictions on test set
# ================================================= 
test.predictions.microb = list()
for (i in 1:5){
  test.predictions.microb[[i]] =
    
    predict(all_prediction_models[[paste0("pred.batch",i)]]$finalModel,
            newdata = as.matrix(test_X[[i]]))
}
## calculate AUCs 
# =================================================
#microbial prediction
library(pROC)
auc.train.microb=list()
for (i in 1:5) {
  auc.train.microb[[i]]=
    auc(roc(train_Y[[i]],ordered(train.predictions.microb[[i]]),direction = "<",levels = c(0,1)))
} 

auc.test.microb=list()
for (i in 1:5) {
  auc.test.microb[[i]]=
    auc(roc(test_Y[[i]],ordered(test.predictions.microb[[i]]),direction = "<",levels = c(0,1)))
}
#Note: ordered(test.predictions.microb[[i]]),针对Error: Predictor must be numeric or ordered.
names(test.predictions.microb[[i]])=NULL
as.numeric(test.predictions.microb[[i]])

#所以随机森林模型的结果不可以采用，预测效果不好

###############1.2 lasso 模型#################

###############1.2.1 species model#################
library(caret)
library(glmnet)
set.seed(1234)
randomized_samples = sample(1:nrow(species_resp_pred))
batches = split(randomized_samples,cut(seq_along(randomized_samples),5,labels = F))


train_X = lapply(batches, function(x) species_resp_pred[-x,-164])
test_X = lapply(batches, function(x) species_resp_pred[x,-164])

train_Y = lapply(batches, function(x) species_resp_pred[-x,164])
test_Y = lapply(batches, function(x) species_resp_pred[x,164])

## Build prediction models
# ===============================================================
all_prediction_models = list()
print(' >> training models')
for(i in 1:5) {
  print(paste0('X-val set ',i,' '))
  
    current_X = train_X[[i]]
    current_Y = train_Y[[i]]
   
    all_prediction_models[[paste0("pred.batch",i)]] = cv.glmnet(x = as.matrix(current_X),
                                                                             y = current_Y,
                                                                             alpha = 1,nfolds = 10,type.measure = "auc",family = "binomial")
  }
## calculate training set and test set predictions
# ========================================= 
# predictions on training set
train.predictions.microb = list()
for (i in 1:5){
  train.predictions.microb[[i]] = 
    # predict(all_prediction_models[[paste0("pred.batch",i,".RData")]],
    #         newx = as.matrix(train_X.microb[[i]]),s="lambda.1se")[,1]
    predict(all_prediction_models[[paste0("pred.batch",i)]],
            newx = as.matrix(train_X[[i]]),s="lambda.1se")[,1]
    
  }
# predictions on test set
# ================================================= 
test.predictions.microb = list()
for (i in 1:5){
  test.predictions.microb[[i]] =
    
    predict(all_prediction_models[[paste0("pred.batch",i)]],
            newx = as.matrix(test_X[[i]]),s="lambda.1se")[,1]
}
## calculate AUCs 
# =================================================
#microbial prediction
library(pROC)
auc.train.microb=list()
 for (i in 1:5) {
  auc.train.microb[[i]]=
   auc(roc(train_Y[[i]],train.predictions.microb[[i]],direction = "<",levels = c(0,1)))
} 

auc.test.microb=list()
 for (i in 1:5) {
   auc.test.microb[[i]]=
    auc(roc(test_Y[[i]],test.predictions.microb[[i]],direction = "<",levels = c(0,1)))
 }
#=============================
#结果输出
setwd("D:/PBC multiomics data/analysis V6/multiomics/result/response prediction/species")
#绘制纳入变量数目,即lambda值与AUC的变化图
for (i in 1:5) {
  png(filename = paste0(i, "_lambda", ".jpg"),width = 300*5,height = 300*5,res = 300)
  plot(all_prediction_models[[paste0("pred.batch",i)]]) 
  
}
dev.off()
plot(all_prediction_models[[paste0("pred.batch",i)]]) 
#绘制每个模型的ROC曲线
plot(roc(test_Y[[i]],test.predictions.microb[[i]],direction = "<",levels = c(0,1)))
for (i in 1:5) {
  png(filename = paste0(i, "_microbe", ".jpg"),width = 300*5,height = 300*5,res = 300)
  plot(roc(test_Y[[i]],test.predictions.microb[[i]],direction = "<",levels = c(0,1)),print.auc=T)
  dev.off()
}
plot(roc(test_Y[[5]],test.predictions.microb[[5]],direction = "<",levels = c(0,1)),print.auc=T)

#
#输出5个预测模型各自纳入的变量及AUC值
all_prediction_microbe_coef= as.data.frame(matrix(0,ncol=1,nrow=164))
for (i in 1:5) {
  all_prediction_microbe_coef[,i+1]=coef(all_prediction_models[[paste0("pred.batch",i)]]$glmnet.fit,exact=F,s=all_prediction_models[[paste0("pred.batch",i)]]$lambda.1se)%>% as.matrix(.)%>% as.data.frame(.)%>%.[,1]
  all_prediction_microbe_coef[165,i+1]=auc.test.microb[[i]]
  }
row.names(all_prediction_microbe_coef)[2:164]=maaslin_species_PBCresp_dif

write.csv(all_prediction_microbe_coef,"all_prediction_microbe_coef.csv")

#==================================================
################1.2.2 serum metabolite model###################

library(caret)
library(glmnet)
set.seed(12306)
randomized_samples_serum = sample(1:nrow(serum_resp_pred))
batches_serum = split(randomized_samples_serum,cut(seq_along(randomized_samples_serum),5,labels = F))


train_X_serum = lapply(batches_serum, function(x) serum_resp_pred[-x,-33])
test_X_serum = lapply(batches_serum, function(x) serum_resp_pred[x,-33])

train_Y_serum = lapply(batches_serum, function(x) serum_resp_pred[-x,33])
test_Y_serum = lapply(batches_serum, function(x) serum_resp_pred[x,33])

## Build prediction models
# ===============================================================
all_prediction_models_serum = list()
print(' >> training models')
for(i in 1:5) {
  print(paste0('   > X-val set ',i,' ...'))
  
  current_X_serum = train_X_serum[[i]]
  current_Y_serum = train_Y_serum[[i]]
  
  #all_prediction_models[["pred.batch",i,".pheno.", j]] = cv.glmnet(as.matrix(complete.X),complete.Y,alpha =1,nfolds = 10,family = "binomial")
  all_prediction_models_serum[[paste0("pred.batch",i)]] = cv.glmnet(x = as.matrix(current_X_serum),
                                                              y = current_Y_serum,
                                                              alpha = 1,nfolds = 5,type.measure = "auc",family = "binomial")
}
## calculate training set and test set predictions
# ========================================= 
# predictions on training set
train.predictions.serum = list()
for (i in 1:5){
  train.predictions.serum[[i]] = 
    # predict(all_prediction_models[[paste0("pred.batch",i,".RData")]],
    #         newx = as.matrix(train_X.microb[[i]]),s="lambda.min")[,1]
    predict(all_prediction_models_serum[[paste0("pred.batch",i)]],
            newx = as.matrix(train_X_serum[[i]]),s="lambda.min")[,1]
  
}
# predictions on test set
# ================================================= 
test.predictions.serum = list()
for (i in 1:5){
  test.predictions.serum[[i]] =
    
    predict(all_prediction_models_serum[[paste0("pred.batch",i)]],
            newx = as.matrix(test_X_serum[[i]]),s="lambda.min")[,1]
}
## calculate AUCs 
# =================================================
#microbial prediction
library(pROC)
auc.train.serum=list()
for (i in 1:5) {
  auc.train.serum[[i]]=
    auc(roc(train_Y_serum[[i]],train.predictions.serum[[i]],direction = "<",levels = c(0,1)))
} 

auc.test.serum=list()
for (i in 1:5) {
  auc.test.serum[[i]]=
    auc(roc(test_Y_serum[[i]],test.predictions.serum[[i]],direction = "<",levels = c(0,1)))
}

auc.test.serum

setwd("D:/PBC multiomics data/analysis V6/multiomics/result/response prediction/serum")
#绘制纳入变量数目,即lambda值与AUC的变化图
for (i in 1:5) {
  png(filename = paste0(i, "_lambda", ".jpg"),width = 300*5,height = 300*5,res = 300)
  plot(all_prediction_models_serum[[paste0("pred.batch",i)]]) 
  dev.off()
}


#绘制每个模型的ROC曲线
plot(roc(test_Y[[5]],test.predictions.microb[[5]],direction = "<",levels = c(0,1)))
for (i in 1:5) {
  png(filename = paste0(i, "_serum.auc", ".jpg"),width = 300*5,height = 300*5,res = 300)
  plot(roc(test_Y_serum[[i]],test.predictions.serum[[i]],direction = "<",levels = c(0,1)),print.auc=T)
  dev.off()
}

plot(roc(test_Y_serum[[3]],test.predictions.serum[[3]],direction = "<",levels = c(0,1)),print.auc=T)

#
#输出5个预测模型各自纳入的变量及AUC值

all_prediction_serum_coef= as.data.frame(matrix(0,ncol=1,nrow=33))
for (i in 1:5) {
  all_prediction_serum_coef[,i+1]=coef(all_prediction_models_serum[[paste0("pred.batch",i)]]$glmnet.fit,exact=F,s=all_prediction_models_serum[[paste0("pred.batch",i)]]$lambda.min)%>% as.matrix(.)%>% as.data.frame(.)%>%.[,1]
  
}
row.names(all_prediction_serum_coef)[2:33]=serum_PBCresp_lm_df
for (i in 1:5) {
  all_prediction_serum_coef[34,i+1]=auc.test.serum[[i]]
}
write.csv(all_prediction_serum_coef,"all_prediction_serum_coef.csv")


#==================================================
################1.2.3 fecal metabolite model###################

library(caret)
library(glmnet)
set.seed(1230)
randomized_samples_fecal = sample(1:nrow(fecal_resp_pred))
batches_fecal = split(randomized_samples_fecal,cut(seq_along(randomized_samples_fecal),5,labels = F))


train_X_fecal = lapply(batches_fecal, function(x) fecal_resp_pred[-x,-48])
test_X_fecal = lapply(batches_fecal, function(x) fecal_resp_pred[x,-48])

train_Y_fecal = lapply(batches_fecal, function(x) fecal_resp_pred[-x,48])
test_Y_fecal = lapply(batches_fecal, function(x) fecal_resp_pred[x,48])

## Build prediction models
# ===============================================================
all_prediction_models_fecal = list()
print(' >> training models')
for(i in 1:5) {
  print(paste0('   > X-val set ',i,' ...'))
  
  current_X_fecal = train_X_fecal[[i]]
  current_Y_fecal = train_Y_fecal[[i]]
  
  #all_prediction_models[["pred.batch",i,".pheno.", j]] = cv.glmnet(as.matrix(complete.X),complete.Y,alpha =1,nfolds = 10,family = "binomial")
  all_prediction_models_fecal[[paste0("pred.batch",i)]] = cv.glmnet(x = as.matrix(current_X_fecal),
                                                                    y = current_Y_fecal,
                                                                    alpha = 1,nfolds = 10,type.measure = "auc",family = "binomial")
}
## calculate training set and test set predictions
# ========================================= 
# predictions on training set
train.predictions.fecal = list()
for (i in 1:5){
  train.predictions.fecal[[i]] = 
    # predict(all_prediction_models[[paste0("pred.batch",i,".RData")]],
    #         newx = as.matrix(train_X.microb[[i]]),s="lambda.1se")[,1]
    predict(all_prediction_models_fecal[[paste0("pred.batch",i)]],
            newx = as.matrix(train_X_fecal[[i]]),s="lambda.1se")[,1]
  
}
# predictions on test set
# ================================================= 
test.predictions.fecal = list()
for (i in 1:5){
  test.predictions.fecal[[i]] =
    
    predict(all_prediction_models_fecal[[paste0("pred.batch",i)]],
            newx = as.matrix(test_X_fecal[[i]]),s="lambda.1se")[,1]
}
## calculate AUCs 
# =================================================
#microbial prediction
library(pROC)
auc.train.fecal=list()
for (i in 1:5) {
  auc.train.fecal[[i]]=
    auc(roc(train_Y_fecal[[i]],train.predictions.fecal[[i]],direction = "<",levels = c(0,1)))
} 

auc.test.fecal=list()
for (i in 1:5) {
  auc.test.fecal[[i]]=
    auc(roc(test_Y_fecal[[i]],test.predictions.fecal[[i]],direction = "<",levels = c(0,1)))
}

auc.test.fecal

setwd("D:/PBC multiomics data/analysis V6/multiomics/result/response prediction/fecal")
#绘制纳入变量数目,即lambda值与AUC的变化图
for (i in 1:5) {
  png(filename = paste0(i, "_lambda", ".jpg"),width = 300*5,height = 300*5,res = 300)
  plot(all_prediction_models_fecal[[paste0("pred.batch",i)]]) 
  dev.off()
}

#绘制每个模型的ROC曲线
plot(roc(test_Y[[5]],test.predictions.microb[[5]],direction = "<",levels = c(0,1)))
for (i in 1:5) {
  png(filename = paste0(i, "_fecal.auc", ".jpg"),width = 300*5,height = 300*5,res = 300)
  plot(roc(test_Y_fecal[[i]],test.predictions.fecal[[i]],direction = "<",levels = c(0,1)),print.auc=T)
  dev.off()
}

plot(roc(test_Y_fecal[[1]],test.predictions.fecal[[1]],direction = "<",levels = c(0,1)),print.auc=T)

#
#输出5个预测模型各自纳入的变量及AUC值

all_prediction_fecal_coef= as.data.frame(matrix(0,ncol=1,nrow=48))
for (i in 1:5) {
  all_prediction_fecal_coef[,i+1]=coef(all_prediction_models_fecal[[paste0("pred.batch",i)]]$glmnet.fit,exact=F,s=all_prediction_models_fecal[[paste0("pred.batch",i)]]$lambda.min)%>% as.matrix(.)%>% as.data.frame(.)%>%.[,1]
  
}
row.names(all_prediction_fecal_coef)[2:48]=fecal_PBCresp_lm_df
for (i in 1:5) {
  all_prediction_fecal_coef[49,i+1]=auc.test.fecal[[i]]
}
write.csv(all_prediction_fecal_coef,"all_prediction_fecal_coef.csv")


#绘制纳入变量数目,即lambda值与AUC的变化图
plot(all_prediction_models_fecal$pred.batch5) 

#==================================================
################1.2.4 combined model###################

library(caret)
library(glmnet)
set.seed(12346)
randomized_samples_combine = sample(1:nrow(combine_resp_pred))
batches_combine = split(randomized_samples_combine,cut(seq_along(randomized_samples_combine),5,labels = F))


train_X_combine = lapply(batches_combine, function(x) combine_resp_pred[-x,-243])
test_X_combine = lapply(batches_combine, function(x) combine_resp_pred[x,-243])

train_Y_combine = lapply(batches_combine, function(x) combine_resp_pred[-x,243])
test_Y_combine = lapply(batches_combine, function(x) combine_resp_pred[x,243])

## Build prediction models
# ===============================================================
all_prediction_models_combine = list()
print(' >> training models')
for(i in 1:5) {
  print(paste0('   > X-val set ',i,' ...'))
  
  current_X_combine = train_X_combine[[i]]
  current_Y_combine = train_Y_combine[[i]]
  
  #all_prediction_models[["pred.batch",i,".pheno.", j]] = cv.glmnet(as.matrix(complete.X),complete.Y,alpha =1,nfolds = 5,family = "binomial")
  all_prediction_models_combine[[paste0("pred.batch",i)]] = cv.glmnet(x = as.matrix(current_X_combine),
                                                                    y = current_Y_combine,
                                                                    alpha = 1,nfolds = 5,type.measure = "auc",family = "binomial")
}
## calculate training set and test set predictions
# ========================================= 
# predictions on training set
train.predictions.combine = list()
for (i in 1:5){
  train.predictions.combine[[i]] = 
    # predict(all_prediction_models[[paste0("pred.batch",i,".RData")]],
    #         newx = as.matrix(train_X.microb[[i]]),s="lambda.min")[,1]
    predict(all_prediction_models_combine[[paste0("pred.batch",i)]],
            newx = as.matrix(train_X_combine[[i]]),s="lambda.min")[,1]
  
}
# predictions on test set
# ================================================= 
test.predictions.combine = list()
for (i in 1:5){
  test.predictions.combine[[i]] =
    
    predict(all_prediction_models_combine[[paste0("pred.batch",i)]],
            newx = as.matrix(test_X_combine[[i]]),s="lambda.min")[,1]
}
## calculate AUCs 
# =================================================
#microbial prediction
library(pROC)
auc.train.combine=list()
for (i in 1:5) {
  auc.train.combine[[i]]=
    auc(roc(train_Y_combine[[i]],train.predictions.combine[[i]],direction = "<",levels = c(0,1)))
} 

auc.test.combine=list()
for (i in 1:5) {
  auc.test.combine[[i]]=
    auc(roc(test_Y_combine[[i]],test.predictions.combine[[i]],direction = "<",levels = c(0,1)))
}

auc.test.combine

setwd("D:/PBC multiomics data/analysis V6/multiomics/result/response prediction/combine")
#绘制纳入变量数目,即lambda值与AUC的变化图
for (i in 1:5) {
  png(filename = paste0(i, "_combine.lambda", ".jpg"),width = 300*5,height = 300*5,res = 300)
  plot(all_prediction_models_combine[[paste0("pred.batch",i)]]) 
  dev.off()
}

#绘制每个模型的ROC曲线
plot(roc(test_Y[[5]],test.predictions.microb[[5]],direction = "<",levels = c(0,1)))
for (i in 1:5) {
  png(filename = paste0(i, "_combine.auc", ".jpg"),width = 300*5,height = 300*5,res = 300)
  plot(roc(test_Y_combine[[i]],test.predictions.combine[[i]],direction = "<",levels = c(0,1)),print.auc=T)
  dev.off()
}

plot(roc(test_Y_combine[[1]],test.predictions.combine[[1]],direction = "<",levels = c(0,1)),print.auc=T)

#============================================
#输出5个预测模型各自纳入的变量及AUC值

all_prediction_combine_coef= as.data.frame(matrix(0,ncol=1,nrow=243))
for (i in 1:5) {
  all_prediction_combine_coef[,i+1]=coef(all_prediction_models_combine[[paste0("pred.batch",i)]]$glmnet.fit,exact=F,s=all_prediction_models_combine[[paste0("pred.batch",i)]]$lambda.min)%>% as.matrix(.)%>% as.data.frame(.)%>%.[,1]
  
}
row.names(all_prediction_combine_coef)[2:243]=colnames(combine_resp_pred)[1:242]

for (i in 1:5) {
  all_prediction_combine_coef[244,i+1]=auc.test.combine[[i]]
}
write.csv(all_prediction_combine_coef,"all_prediction_combine_coef.csv")

#绘制纳入变量数目,即lambda值与AUC的变化图
plot(all_prediction_models_combine$pred.batch3) 
plot(roc(test_Y_combine[[3]],test.predictions.combine[[3]],direction = "<",levels = c(0,1)))


################可视化#################
#1.比较/可视化不同模型的预测效率
fecal_auc=all_prediction_fecal_coef[49,2:6]
serum_auc=all_prediction_serum_coef[34,2:6]
microbe_auc=all_prediction_microbe_coef[165,2:6]
combine_auc=all_prediction_combine_coef[244,2:6]

auc_all=data.frame(microbeAUC=t(microbe_auc[1,]),fecal_AUC=t(fecal_auc[1,]),serum_AUC=t(serum_auc[1,]),combine_AUC=t(combine_auc[1,]))
colnames(auc_all)=c("microbe","fecal","serum","combine")
row.names(auc_all)=c("Fold1","Fold2","Fold3","Fold4","Fold5")
library(reshape)
library(reshape2)
melt(auc_all)
ggplot(data = melt(auc_all),aes(x=variable,y=value,fill=variable))+
  geom_boxplot()+
  scale_fill_manual(values= c("#EF4B4C","#486C87","#765BB9","#C6C6C6"))+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")

#2.可视化在species model fold5中 contribution top5的species
setwd("D:/PBC multiomics data/analysis V6/multiomics/plot/prediction model")
top5_species_contribu=read.csv("top5 species contribution.csv",header = T,row.names = 1)
top5_species_contribu$species=row.names(top5_species_contribu)
top5_species_contribu$species=factor(top5_species_contribu$species,levels= unique(top5_species_contribu$species))
library(ggplot2)
ggplot(top5_species_contribu,aes(x=species,y=contribution))+
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = -10,hjust = 0.4,vjust =1))
  
#3.可视化在combined model中 contribution top10的feature
setwd("D:/PBC multiomics data/analysis V6/multiomics/plot/prediction model")
top10_feature_contribu=read.csv("top10 feature contribution.csv",header = T,row.names = 1)
top10_feature_contribu$feature=row.names(top10_feature_contribu)
top10_feature_contribu$feature=factor(top10_feature_contribu$feature,levels= unique(top10_feature_contribu$feature))
library(ggplot2)
ggplot(top10_feature_contribu,aes(x=feature,y=contribution))+
  geom_bar(stat = "identity",fill="#D36D2D",width = 0.6)+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  coord_flip()
  

###################1.2.5 clinic model####################
clinic_resp_pred_impu=clinic_resp_pred[,1:10]
for (i in 1:10) {
  clinic_resp_pred_impu[,i]=impute(clinic_resp_pred[,i])
}
clinic_resp_pred_impu[,11]=clinic_resp_pred[,11]

str(clinic_resp_pred_impu)
library(caret)
library(glmnet)
set.seed(124)
randomized_samples = sample(1:nrow(clinic_resp_pred_impu))
batches = split(randomized_samples,cut(seq_along(randomized_samples),5,labels = F))


train_X = lapply(batches, function(x) clinic_resp_pred_impu[-x,-11])
test_X = lapply(batches, function(x) clinic_resp_pred_impu[x,-11])

train_Y = lapply(batches, function(x) clinic_resp_pred_impu[-x,11])
test_Y = lapply(batches, function(x) clinic_resp_pred_impu[x,11])

## Build prediction models
# ===============================================================
all_prediction_models = list()
print(' >> training models')
for(i in 1:5) {
  print(paste0('   > X-val set ',i,' ...'))
  
  current_X = train_X[[i]]
  current_Y = train_Y[[i]]
  
  #all_prediction_models[["pred.batch",i,".pheno.", j]] = cv.glmnet(as.matrix(complete.X),complete.Y,alpha = 0.5,nfolds = 10,family = "binomial")
  all_prediction_models[[paste0("pred.batch",i)]] = cv.glmnet(x = as.matrix(current_X),
                                                              y = current_Y,
                                                              alpha = 1,nfolds = 10,type.measure = "auc",family = "binomial")
}
## calculate training set and test set predictions
# ========================================= 
# predictions on training set
train.predictions.clinic = list()
for (i in 1:5){
  train.predictions.clinic[[i]] = 
    # predict(all_prediction_models[[paste0("pred.batch",i,".RData")]],
    #         newx = as.matrix(train_X.microb[[i]]),s="lambda.min")[,1]
    predict(all_prediction_models[[paste0("pred.batch",i)]],
            newx = as.matrix(train_X[[i]]),s="lambda.min")[,1]
  
}
# predictions on test set
# ================================================= 
test.predictions.clinic = list()
for (i in 1:5){
  test.predictions.clinic[[i]] =
    
    predict(all_prediction_models[[paste0("pred.batch",i)]],
            newx = as.matrix(test_X[[i]]),s="lambda.min")[,1]
}
## calculate AUCs 
# =================================================
#microbial prediction
library(pROC)
auc.train.clinic=list()
for (i in 1:5) {
  auc.train.clinic[[i]]=
    auc(roc(train_Y[[i]],train.predictions.clinic[[i]],direction = "<",levels = c(0,1)))
} 

auc.test.clinic=list()
for (i in 1:5) {
  auc.test.clinic[[i]]=
    auc(roc(test_Y[[i]],test.predictions.clinic[[i]],direction = "<",levels = c(0,1)))
}

############1.2.6 clostridial/veillonella###################
Veillonella=genus_rela_filter["k__Bacteria;p__Firmicutes;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella",c(PBCbase,HC_new,PBC_6m_pair)]
Clostridial=order_rela_filter["k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales",c(PBCbase,HC_new,PBC_6m_pair)]
Vei_clos=rbind(Veillonella,Clostridial)
row.names(Vei_clos)
colnames(Vei_clos)
Vei_clos=t(Vei_clos) %>% as.data.frame(.)
Vei_clos$clos_vei_ratio=Vei_clos[,2]/Vei_clos[,1]
#1.clinic+Vei_clos_ratio
clinic.CVratio_resp_pred_impu=cbind(clinic_resp_pred_impu[,1:10],Vei_clos[row.names(clinic_resp_pred_impu),3])
clinic.CVratio_resp_pred_impu=cbind(clinic.CVratio_resp_pred_impu,clinic_resp_pred_impu[,11])
colnames(clinic.CVratio_resp_pred_impu)[11]="Clo_Vei_ratio"
colnames(clinic.CVratio_resp_pred_impu)[12]="response"

library(caret)
library(glmnet)
set.seed(124)
randomized_samples = sample(1:nrow(clinic.CVratio_resp_pred_impu))
batches = split(randomized_samples,cut(seq_along(randomized_samples),5,labels = F))


train_X = lapply(batches, function(x) clinic.CVratio_resp_pred_impu[-x,-12])
test_X = lapply(batches, function(x) clinic.CVratio_resp_pred_impu[x,-12])

train_Y = lapply(batches, function(x) clinic.CVratio_resp_pred_impu[-x,12])
test_Y = lapply(batches, function(x) clinic.CVratio_resp_pred_impu[x,12])

## Build prediction models
all_prediction_models = list()
print(' >> training models')
for(i in 1:5) {
  print(paste0('   > X-val set ',i,' ...'))
  
  current_X = train_X[[i]]
  current_Y = train_Y[[i]]
  
  #all_prediction_models[["pred.batch",i,".pheno.", j]] = cv.glmnet(as.matrix(complete.X),complete.Y,alpha = 0.5,nfolds = 10,family = "binomial")
  all_prediction_models[[paste0("pred.batch",i)]] = cv.glmnet(x = as.matrix(current_X),
                                                              y = current_Y,
                                                              alpha = 1,nfolds = 10,type.measure = "auc",family = "binomial")
}
## calculate training set and test set predictions
# predictions on training set
train.predictions.clinic_ratio = list()
for (i in 1:5){
  train.predictions.clinic_ratio[[i]] = 
    # predict(all_prediction_models[[paste0("pred.batch",i,".RData")]],
    #         newx = as.matrix(train_X.microb[[i]]),s="lambda.min")[,1]
    predict(all_prediction_models[[paste0("pred.batch",i)]],
            newx = as.matrix(train_X[[i]]),s="lambda.min")[,1]
  
}
# predictions on test set
test.predictions.clinic_ratio = list()
for (i in 1:5){
  test.predictions.clinic_ratio[[i]] =
    
    predict(all_prediction_models[[paste0("pred.batch",i)]],
            newx = as.matrix(test_X[[i]]),s="lambda.min")[,1]
}
## calculate AUCs 
#clinic+Clostridla/Veillonella ratio prediction
library(pROC)
auc.train.clinic_ratio=list()
for (i in 1:5) {
  auc.train.clinic_ratio[[i]]=
    auc(roc(train_Y[[i]],train.predictions.clinic_ratio[[i]],direction = "<",levels = c(0,1)))
} 

auc.test.clinic_ratio=list()
for (i in 1:5) {
  auc.test.clinic_ratio[[i]]=
    auc(roc(test_Y[[i]],test.predictions.clinic_ratio[[i]],direction = "<",levels = c(0,1)))
}


#1.3循环计算PCo1_corR>0.8的taxa的预测AUC *figure############
#读入PCo1_corR>0.8的taxa名称
species_corR_over0.8=read.csv("species_corR_over0.8.csv",header = T,row.names = 1)
species_corR_over0.8=species_corR_over0.8$species

##1.根据相对丰度
species_resp_pred_corR=species_rela_filter[species_corR_over0.8,c(PBCbase_resp_6m,PBCbase_noresp_6m)] %>% t(.) %>% as.data.frame(.)
species_resp_pred_corR=cbind(species_resp_pred_corR,metadata[c(PBCbase_resp_6m,PBCbase_noresp_6m),"Response_6m"])
colnames(species_resp_pred_corR)[34]="response"

AUC_species_corR=vector()
for (i in names(species_resp_pred_corR)[1:33]) {
  cl <- makeCluster(2)
  registerDoSNOW(cl)
  set.seed(100)
  newdata=species_resp_pred_corR[,c(i,"response")]
  colnames(newdata)=c("taxa","response")
  opt_cut_b=cutpointr(newdata, taxa, response, boot_runs = 1000,
                      silent = TRUE, allowParallel = TRUE)
  stopCluster(cl)
 
  print(opt_cut_b$AUC)
  
  }
##2.根据CLR
species_resp_pred_corR_clr=species_abs_clr[c(PBCbase_resp_6m,PBCbase_noresp_6m),species_corR_over0.8]
species_resp_pred_corR_clr=cbind(species_resp_pred_corR_clr,metadata[c(PBCbase_resp_6m,PBCbase_noresp_6m),"Response_6m"])
colnames(species_resp_pred_corR_clr)[34]="response"


for (i in names(species_resp_pred_corR_clr)[1:33]) {
  cl <- makeCluster(2)
  registerDoSNOW(cl)
  set.seed(100)
  newdata=species_resp_pred_corR_clr[,c(i,"response")]
  colnames(newdata)=c("taxa","response")
  opt_cut_b=cutpointr(newdata, taxa, response, boot_runs = 1000,
                      silent = TRUE, allowParallel = TRUE)
  stopCluster(cl)
  
  print(opt_cut_b$AUC)
  
}
#可视化33个species的AUC
AUC_species_corR=read.csv("predicting AUC of 33 species.csv",header = T,row.names = 1)
AUC_species_corR$species=row.names(AUC_species_corR)
AUC_species_corR$species=factor(AUC_species_corR$species,levels= unique(AUC_species_corR$species))
library(ggplot2)
ggplot(AUC_species_corR,aes(x=species,y=AUC))+
  geom_bar(stat = "identity",fill="#D36D2D",width = 0.6)+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "none")+
  coord_flip()

#1.3.1重复Clostridium_leptum的区分效能及cutoff *figure################

#1.relative abundance
library("doSNOW")
library("doRNG")
library("tidyverse")
cl <- makeCluster(2)
registerDoSNOW(cl)
set.seed(100)
opt_cut_b=cutpointr(species_resp_pred_corR, "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__unidentified;s__[Clostridium] leptum", response, boot_runs = 1000,
          silent = TRUE, allowParallel = TRUE)
stopCluster(cl)
opt_cut_b %>% select(optimal_cutpoint, data, boot)
summary(opt_cut_b)
plot(opt_cut_b)
library(dplyr)
grep("leptum",colnames(species_resp_pred_corR))
clostrium_leptum=species_resp_pred_corR[,c(25,34)]
colnames(clostrium_leptum)[1]="leptum"
PBCbase_clos.leptum.high=filter(clostrium_leptum,leptum> 1e-04) %>% row.names(.)
PBCbase_clos.leptum.low=filter(clostrium_leptum,leptum< 1e-04) %>% row.names(.)
PBCbase_Clostridium_leptum=c(rep("PBCbase_clos.leptum.high",66),rep("PBCbase_clos.leptum.low",65)) %>% as.data.frame(.)

row.names(PBCbase_Clostridium_leptum)=c(PBCbase_clos.leptum.high,PBCbase_clos.leptum.low)
colnames(PBCbase_Clostridium_leptum)="Clostridium_leptum_rela"

#1.2 比较两组的临床指标的差异
colnames(metadata_PBCbase)
metadata_PBCbase=cbind(metadata_PBCbase[,1:39],PBCbase_Clostridium_leptum[row.names(metadata_PBCbase),1])
colnames(metadata_PBCbase)[40]="Clostridium_leptum_rela"

chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"Clostridium_leptum_rela"]) #P=0.049
chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"Clostridium_leptum_rela"])$observed


#2.CLR abundance
library("doSNOW")
library("doRNG")
library("tidyverse")
cl <- makeCluster(2)
registerDoSNOW(cl)
set.seed(100)
opt_cut_b=cutpointr(species_resp_pred_corR_clr, "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__unidentified;s__[Clostridium] leptum", response, boot_runs = 1000,
                    silent = TRUE, allowParallel = TRUE)
stopCluster(cl)
opt_cut_b %>% select(optimal_cutpoint, data, boot)
summary(opt_cut_b)
plot(opt_cut_b)
?plot_cut_boot
plot_cut_boot(opt_cut_b)

library(dplyr)
grep("leptum",colnames(species_resp_pred_corR))
clostrium_leptum=species_resp_pred_corR_clr[,c(25,34)]
colnames(clostrium_leptum)[1]="leptum"
PBCbase_clos.leptum.high=filter(clostrium_leptum,leptum>  -0.2693  ) %>% row.names(.)
PBCbase_clos.leptum.low=filter(clostrium_leptum,leptum<  -0.2693  ) %>% row.names(.)
PBCbase_Clostridium_leptum=c(rep("PBCbase_clos.leptum.high",79),rep("PBCbase_clos.leptum.low",52)) %>% as.data.frame(.)

row.names(PBCbase_Clostridium_leptum)=c(PBCbase_clos.leptum.high,PBCbase_clos.leptum.low)
colnames(PBCbase_Clostridium_leptum)="Clostridium_leptum_clr"

#1.2 比较两组的临床指标的差异
colnames(metadata_PBCbase)
metadata_PBCbase=cbind(metadata_PBCbase[,1:40],PBCbase_Clostridium_leptum[row.names(metadata_PBCbase),1])
colnames(metadata_PBCbase)[41]="Clostridium_leptum_clr"

chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"Clostridium_leptum_clr"]) #P=0.006
chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"Clostridium_leptum_clr"])$observed

#绘图
#1.donuts plot according to response (only CLR result)
# load library
library(ggplot2)

# Create  data.
#clostridum_high group
data_CL.high <- data.frame(
  category=c("response", "no-response"),
  count=c(63, 16)
)

# Compute percentages
data_CL.high$fraction <- data_CL.high$count / sum(data_CL.high$count)
data_CL.high$fraction=c("80%","20%")
# Compute the cumulative percentages (top of each rectangle)
data_CL.high$ymax <- cumsum(data_CL.high$fraction)

# Compute the bottom of each rectangle
data_CL.high$ymin <- c(0, head(data_CL.high$ymax, n=-1))

# Compute label position
data_CL.high$labelPosition <- (data_CL.high$ymax + data_CL.high$ymin) / 2

# Compute a good label
data_CL.high$label <- paste0(data_CL.high$category,"\n", data_CL.high$fraction)

# Make the plot
ggplot(data_CL.high, aes(ymax=ymax, ymin=ymin, xmax=3, xmin=2, fill=category)) +
  geom_rect() +
  geom_label( x=4, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values = c("#BFBFBF","#85ABD1")) +
  coord_polar(theta="y") +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")

#clostridum_leptum low group
data_CL.low <- data.frame(
  category=c("response", "no-response"),
  count=c(29, 23)
)

# Compute percentages
data_CL.low$fraction <- data_CL.low$count / sum(data_CL.low$count)

# Compute the cumulative percentages (top of each rectangle)
data_CL.low$ymax <- cumsum(data_CL.low$fraction)

# Compute the bottom of each rectangle
data_CL.low$ymin <- c(0, head(data_CL.low$ymax, n=-1))

# Compute label position
data_CL.low$labelPosition <- (data_CL.low$ymax + data_CL.low$ymin) / 2

# Compute a good label
data_CL.low$fraction=c("56%","44%")
data_CL.low$label <- paste0(data_CL.low$category,"\n", data_CL.low$fraction)

# Make the plot
ggplot(data_CL.low, aes(ymax=ymax, ymin=ymin, xmax=3, xmin=2, fill=category)) +
  geom_rect() +
  geom_label( x=4, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values = c("#BFBFBF","#85ABD1")) +
  coord_polar(theta="y") +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")

#2.样本分布 stacked plot
# Create data
sample_CL <- data.frame(
  category=c("C.leptum high", "C.leptum low"),
  count=c(79, 52),
  percentage=c(0.6,0.4),
  group=c("PBCbase","PBCbase")
)

ggplot( sample_CL, aes( x = group, weight = percentage, fill = category))+
  geom_bar( width = 0.2,position = position_stack(reverse = TRUE))+
  scale_fill_manual(values = c("#088288","#F19433"))+theme_bw()


#1.3.2 重复Eubacterium_siraeum的区分效能及cutoff *figure################

#1.relative abundance
library("doSNOW")
library("doRNG")
library("tidyverse")
library(cutpointr)
cl <- makeCluster(2)
registerDoSNOW(cl)
set.seed(100)
opt_cut_b=cutpointr(species_resp_pred_corR, "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__unidentified;s__[Eubacterium] siraeum", response, boot_runs = 1000,
                    silent = TRUE, allowParallel = TRUE)
stopCluster(cl)

summary(opt_cut_b) #cutpoint=1e-04
plot(opt_cut_b)
library(dplyr)
grep("siraeum",colnames(species_resp_pred_corR))
Eubacterium_siraeum=species_resp_pred_corR[,c(29,34)]
colnames(Eubacterium_siraeum)[1]="siraeum"
PBCbase_Eubacterium_siraeum.high=filter(Eubacterium_siraeum,siraeum> 1e-04) %>% row.names(.)
PBCbase_Eubacterium_siraeum.low=filter(Eubacterium_siraeum,siraeum< 1e-04) %>% row.names(.)
PBCbase_Eubacterium_siraeum=c(rep("PBCbase_Eubacterium_siraeum.high",61),rep("PBCbase_Eubacterium_siraeum.low",70)) %>% as.data.frame(.)

row.names(PBCbase_Eubacterium_siraeum)=c(PBCbase_Eubacterium_siraeum.high,PBCbase_Eubacterium_siraeum.low)
colnames(PBCbase_Eubacterium_siraeum)="Eubacterium_siraeum_rela"

#1.2 比较两组的临床指标的差异
colnames(metadata_PBCbase)
metadata_PBCbase=cbind(metadata_PBCbase[,1:41],PBCbase_Eubacterium_siraeum[row.names(metadata_PBCbase),1])
colnames(metadata_PBCbase)[42]="Eubacterium_siraeum_rela"

chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"Eubacterium_siraeum_rela"]) #P=0.07
chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"Eubacterium_siraeum_rela"])$observed


#2.CLR abundance
library("doSNOW")
library("doRNG")
library("tidyverse")
cl <- makeCluster(2)
registerDoSNOW(cl)
set.seed(100)
opt_cut_b_ES_clr=cutpointr(species_resp_pred_corR_clr, "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__unidentified;s__[Eubacterium] siraeum", response, boot_runs = 1000,
                    silent = TRUE, allowParallel = TRUE)
stopCluster(cl)

summary(opt_cut_b_ES_clr) #cutpoint= -0.2507 
plot(opt_cut_b_ES_clr)
library(dplyr)
grep("siraeum",colnames(species_resp_pred_corR_clr))
Eubacterium_siraeum=species_resp_pred_corR_clr[,c(29,34)]
colnames(Eubacterium_siraeum)[1]="siraeum"
PBCbase_Eubacterium_siraeum.high=filter(Eubacterium_siraeum,siraeum>  -0.2507   ) %>% row.names(.)
PBCbase_Eubacterium_siraeum.low=filter(Eubacterium_siraeum,siraeum<  -0.2507   ) %>% row.names(.)
PBCbase_Eubacterium_siraeum=c(rep("PBCbase_Eubacterium_siraeum.high",72),rep("PBCbase_Eubacterium_siraeum.low",59)) %>% as.data.frame(.)

row.names(PBCbase_Eubacterium_siraeum)=c(PBCbase_Eubacterium_siraeum.high,PBCbase_Eubacterium_siraeum.low)
colnames(PBCbase_Eubacterium_siraeum)="Eubacterium_siraeum_clr"

#1.2 比较两组的临床指标的差异
colnames(metadata_PBCbase)
metadata_PBCbase=cbind(metadata_PBCbase[,1:42],PBCbase_Eubacterium_siraeum[row.names(metadata_PBCbase),1])
colnames(metadata_PBCbase)[43]="Eubacterium_siraeum_clr"

chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"Eubacterium_siraeum_clr"]) #P=0.002
chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"Eubacterium_siraeum_clr"])$observed

#绘图
#1.dougnut plot according to UDCA response (only CLR result)
# load library
library(ggplot2)

# Create  data.
#Eubacterium_siraeum_high group
data_ES.high <- data.frame(
  category=c("response", "no-response"),
  count=c(59, 13)
)

# Compute percentages
data_ES.high$fraction <- data_ES.high$count / sum(data_ES.high$count)

# Compute the cumulative percentages (top of each rectangle)
data_ES.high$ymax <- cumsum(data_ES.high$fraction)

# Compute the bottom of each rectangle
data_ES.high$ymin <- c(0, head(data_ES.high$ymax, n=-1))

# Compute label position
data_ES.high$labelPosition <- (data_ES.high$ymax + data_ES.high$ymin) / 2

# Compute a good label
data_ES.high$fraction=c("82%","18%")
data_ES.high$label <- paste0(data_ES.high$category,"\n", data_ES.high$fraction)

# Make the plot
ggplot(data_ES.high, aes(ymax=ymax, ymin=ymin, xmax=3, xmin=2, fill=category)) +
  geom_rect() +
  geom_label( x=4, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values = c("#BFBFBF","#85ABD1")) +
  coord_polar(theta="y") +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")

#Eubacterium_siraeum low group
data_ES.low <- data.frame(
  category=c("response", "no-response"),
  count=c(33, 26)
)

# Compute percentages
data_ES.low$fraction <- data_ES.low$count / sum(data_ES.low$count)

# Compute the cumulative percentages (top of each rectangle)
data_ES.low$ymax <- cumsum(data_ES.low$fraction)

# Compute the bottom of each rectangle
data_ES.low$ymin <- c(0, head(data_ES.low$ymax, n=-1))

# Compute label position
data_ES.low$labelPosition <- (data_ES.low$ymax + data_ES.low$ymin) / 2

# Compute a good label
data_ES.low$fraction=c("56%","44%")
data_ES.low$label <- paste0(data_ES.low$category,"\n", data_ES.low$fraction)

# Make the plot
ggplot(data_ES.low, aes(ymax=ymax, ymin=ymin, xmax=3, xmin=2, fill=category)) +
  geom_rect() +
  geom_label( x=4, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values = c("#BFBFBF","#85ABD1")) +
  coord_polar(theta="y") +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")

#2.样本分布 stacked plot
# Create data
sample_ES <- data.frame(
  category=c("E.siraeum high", "E.siraeum low"),
  count=c(72, 59),
  percentage=c(0.55,0.45),
  group=c("PBCbase","PBCbase")
)

ggplot( sample_ES, aes( x = group, weight = percentage, fill = category))+
  geom_bar( position = position_stack(reverse = TRUE))+
  scale_fill_manual(values = c("#088288","#F19433"))+theme_bw()

#############1.3.3 Ruminococcus spp 的区分效能#####################
library("doSNOW")
library("doRNG")
library(cutpointr)
library("tidyverse")
cl <- makeCluster(2)
registerDoSNOW(cl)
set.seed(100)
opt_cut_b_ES_clr=cutpointr(species_resp_pred_corR_clr, "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminococcus;s__Ruminococcus sp. CAG:563", response, boot_runs = 1000,
                           silent = TRUE, allowParallel = TRUE)
stopCluster(cl)

summary(opt_cut_b_ES_clr) #cutpoint=  1.1573
plot(opt_cut_b_ES_clr)
library(dplyr)
grep("563",colnames(species_resp_pred_corR_clr))
Ruminococcus_563=species_resp_pred_corR_clr[,c(15,34)]
colnames(Ruminococcus_563)[1]="Ruminococcus_563"
PBCbase_Ruminococcus_563.high=filter(Ruminococcus_563,Ruminococcus_563>=  1.1573   ) %>% row.names(.)
PBCbase_Ruminococcus_563.low=filter(Ruminococcus_563,Ruminococcus_563<=  1.1573   ) %>% row.names(.)
PBCbase_Ruminococcus_563=c(rep("PBCbase_Ruminococcus_563.high",54),rep("PBCbase_Ruminococcus_563.low",77)) %>% as.data.frame(.)

row.names(PBCbase_Ruminococcus_563)=c(PBCbase_Ruminococcus_563.high,PBCbase_Ruminococcus_563.low)
colnames(PBCbase_Ruminococcus_563)="Ruminococcus_563_clr"

#1.2 比较两组的临床指标的差异
colnames(metadata_PBCbase)
metadata_PBCbase=cbind(metadata_PBCbase[,1:43],PBCbase_Ruminococcus_563[row.names(metadata_PBCbase),1])
colnames(metadata_PBCbase)[44]="Ruminococcus_563"

chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"Ruminococcus_563"]) #P=0.002
chisq.test(metadata_PBCbase[,"Response_6m"],metadata_PBCbase[,"Ruminococcus_563"])$observed

#绘图

#1.dougnut plot according to UDCA response (only CLR result)
# load library
library(ggplot2)

# Create  data.
#1.Ruminococcus high group
data_Ru.high <- data.frame(
  category=c("response", "no-response"),
  count=c(45, 9)
)
data_Ru.high$fraction <- data_Ru.high$count / sum(data_Ru.high$count)
# Compute the cumulative percentages (top of each rectangle)
data_Ru.high$ymax <- cumsum(data_Ru.high$fraction)
# Compute the bottom of each rectangle
data_Ru.high$ymin <- c(0, head(data_Ru.high$ymax, n=-1))

data_Ru.high$labelPosition <- (data_Ru.high$ymax + data_Ru.high$ymin) / 2
data_Ru.high$label <- paste0(data_Ru.high$category,"\n", data_Ru.high$fraction)
# Make the plot
ggplot(data_Ru.high, aes(ymax=ymax, ymin=ymin, xmax=3, xmin=2, fill=category)) +
  geom_rect() +
  geom_label( x=4, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values = c("#85ABD1","#BFBFBF")) +
  coord_polar(theta="y") +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")


#2.Ruminococcus low group
data_Ru.low <- data.frame(
  category=c("response", "no-response"),
  count=c(47, 30)
)
data_Ru.low$fraction <- data_Ru.low$count / sum(data_Ru.low$count)
# Compute the cumulative percentages (top of each rectangle)
data_Ru.low$ymax <- cumsum(data_Ru.low$fraction)
# Compute the bottom of each rectangle
data_Ru.low$ymin <- c(0, head(data_Ru.low$ymax, n=-1))

data_Ru.low$labelPosition <- (data_Ru.low$ymax + data_Ru.low$ymin) / 2
data_Ru.low$label <- paste0(data_Ru.low$category,"\n", data_Ru.low$fraction)
# Make the plot
ggplot(data_Ru.low, aes(ymax=ymax, ymin=ymin, xmax=3, xmin=2, fill=category)) +
  geom_rect() +
  geom_label( x=4, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values = c("#85ABD1","#BFBFBF")) +
  coord_polar(theta="y") +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")

##########1.3.4 cluster 与C.leptum and E.siraeum的一致性################
PBCbase_strati.1=data.frame(
  C.leptum=c(rep("CL.low",52),rep("CL.high",79))
)
row.names(PBCbase_strati.1)=c(PBCbase_clos.leptum.low,PBCbase_clos.leptum.high)

PBCbase_strati.2=data.frame(
 E.siraeum=c(rep("ES.low",59),rep("ES.high",72))
)
row.names(PBCbase_strati.2)=c(PBCbase_Eubacterium_siraeum.low,PBCbase_Eubacterium_siraeum.high)

PBCbase_strati=cbind(PBCbase_strati.1[,1],PBCbase_strati.2[row.names(PBCbase_strati.1),1])
row.names(PBCbase_strati)=row.names(PBCbase_strati.1)
PBCbase_strati=as.data.frame(PBCbase_strati)
PBCbase_strati=cbind(PBCbase_strati[,1:2],metadata_PBCHC[row.names(PBCbase_strati),"cluster"])
colnames(PBCbase_strati)=c("C.leptum","E.siraeum","cluster")

chisq.test(PBCbase_strati[,"C.leptum"],PBCbase_strati[,"cluster"])$observed
C.leptum_clt=chisq.test(PBCbase_strati[,"C.leptum"],PBCbase_strati[,"cluster"])$observed %>% as.data.frame(.)
colnames(C.leptum_clt)=c("C.leptum","cluster","number")
C.leptum_clt$proportion=c(0.8,0.13,0.2,0.87)
C.leptum_clt=C.leptum_clt[c(2,1,4,3),]
C.leptum_clt$cluster=factor(C.leptum_clt$cluster,levels = c("clt1","clt2"))
C.leptum_clt$C.leptum=factor(C.leptum_clt$C.leptum,levels = c("CL.low","CL.high"))
str(C.leptum_clt)
ggplot(C.leptum_clt,aes(x=C.leptum,y=proportion,fill=cluster))+
  geom_bar(stat = "identity",position = "stack")+
  scale_fill_manual(values = c("#088288","#F19433"))+
  geom_text(aes(label=proportion),position = position_stack(),color="white",vjust=2)+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank())


E.siraeum_clt=chisq.test(PBCbase_strati[,"E.siraeum"],PBCbase_strati[,"cluster"]) $observed %>% as.data.frame(.)
colnames(E.siraeum_clt)=c("E.siraeum","cluster","number")
E.siraeum_clt$proportion=c(0.87,0.12,0.13,0.88)
E.siraeum_clt=E.siraeum_clt[c(2,1,4,3),]
E.siraeum_clt$E.siraeum=factor(E.siraeum_clt$E.siraeum,levels = c("ES.low","ES.high"))
ggplot(E.siraeum_clt,aes(x=E.siraeum,y=proportion,fill=cluster))+
  geom_bar(stat = "identity",position = "stack")+
  scale_fill_manual(values = c("#088288","#F19433"))+
  geom_text(aes(label=proportion),position = position_stack(),color="white",vjust=2)+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank())


############1.4 Response prediction lasso model based on 33 species#################
#============1.4.1 lasso five fold-cross validation#############

library(caret)
library(glmnet)
set.seed(12)
randomized_samples = sample(1:nrow(species_resp_pred_corR_clr))
batches = split(randomized_samples,cut(seq_along(randomized_samples),5,labels = F))


train_X = lapply(batches, function(x) species_resp_pred_corR_clr[-x,-34])
test_X = lapply(batches, function(x) species_resp_pred_corR_clr[x,-34])

train_Y = lapply(batches, function(x) species_resp_pred_corR_clr[-x,34])
test_Y = lapply(batches, function(x) species_resp_pred_corR_clr[x,34])


## Build prediction models
all_prediction_models = list()
print(' >> training models')
for(i in 1:5) {
  print(paste0('X-val set ',i,' '))
  
  current_X = train_X[[i]]
  current_Y = train_Y[[i]]
  
  all_prediction_models[[paste0("pred.batch",i)]] = cv.glmnet(x = as.matrix(current_X),
                                                              y = current_Y,
                                                              alpha = 1,nfolds = 10,type.measure = "auc",family = "binomial")
}
## calculate training set and test set predictions
# predictions on training set
train.predictions.microb = list()
for (i in 1:5){
  train.predictions.microb[[i]] = 
    # predict(all_prediction_models[[paste0("pred.batch",i,".RData")]],
    #         newx = as.matrix(train_X.microb[[i]]),s="lambda.1se")[,1]
    predict(all_prediction_models[[paste0("pred.batch",i)]],
            newx = as.matrix(train_X[[i]]),s="lambda.min")[,1]
  
}
# predictions on test set
test.predictions.microb = list()
for (i in 1:5){
  test.predictions.microb[[i]] =
    
    predict(all_prediction_models[[paste0("pred.batch",i)]],
            newx = as.matrix(test_X[[i]]),s="lambda.min")[,1]
}
## calculate AUCs 
#microbial prediction
library(pROC)
auc.train.microb=list()
for (i in 1:5) {
  auc.train.microb[[i]]=
    auc(roc(train_Y[[i]],train.predictions.microb[[i]],direction = "<",levels = c(0,1)))
} 

auc.test.microb=list()
for (i in 1:5) {
  auc.test.microb[[i]]=
    auc(roc(test_Y[[i]],test.predictions.microb[[i]],direction = "<",levels = c(0,1)))
}


############1.5 Response prediction lasso model based on 5 species#################

species_corR_over0.8_respdif= intersect(species_corR_over0.8,maaslin_species_PBCresp_dif)

species_resp_pred_corR_clr_respdif=species_resp_pred_corR_clr[,c(species_corR_over0.8_respdif,"response")]
#============1.5.1 lasso five fold-cross validation#############

library(caret)
library(glmnet)
set.seed(1234)
randomized_samples = sample(1:nrow(species_resp_pred_corR_clr_respdif))
batches = split(randomized_samples,cut(seq_along(randomized_samples),5,labels = F))


train_X = lapply(batches, function(x) species_resp_pred_corR_clr_respdif[-x,-6])
test_X = lapply(batches, function(x) species_resp_pred_corR_clr_respdif[x,-6])

train_Y = lapply(batches, function(x) species_resp_pred_corR_clr_respdif[-x,6])
test_Y = lapply(batches, function(x) species_resp_pred_corR_clr_respdif[x,6])


## Build prediction models
all_prediction_models = list()
print(' >> training models')
for(i in 1:5) {
  print(paste0('X-val set ',i,' '))
  
  current_X = train_X[[i]]
  current_Y = train_Y[[i]]
  
  all_prediction_models[[paste0("pred.batch",i)]] = cv.glmnet(x = as.matrix(current_X),
                                                              y = current_Y,
                                                              alpha = 1,nfolds = 10,type.measure = "auc",family = "binomial")
}
## calculate training set and test set predictions
# predictions on training set
train.predictions.microb = list()
for (i in 1:5){
  train.predictions.microb[[i]] = 
    # predict(all_prediction_models[[paste0("pred.batch",i,".RData")]],
    #         newx = as.matrix(train_X.microb[[i]]),s="lambda.1se")[,1]
    predict(all_prediction_models[[paste0("pred.batch",i)]],
            newx = as.matrix(train_X[[i]]),s="lambda.min")[,1]
  
}
# predictions on test set
test.predictions.microb = list()
for (i in 1:5){
  test.predictions.microb[[i]] =
    
    predict(all_prediction_models[[paste0("pred.batch",i)]],
            newx = as.matrix(test_X[[i]]),s="lambda.min")[,1]
}
## calculate AUCs 
#microbial prediction
library(pROC)
auc.train.microb=list()
for (i in 1:5) {
  auc.train.microb[[i]]=
    auc(roc(train_Y[[i]],train.predictions.microb[[i]],direction = "<",levels = c(0,1)))
} 

auc.test.microb=list()
for (i in 1:5) {
  auc.test.microb[[i]]=
    auc(roc(test_Y[[i]],test.predictions.microb[[i]],direction = "<",levels = c(0,1)))
}


############1.6 Response prediction lasso model based on genus and genus ratio#################
genus_resp_pred <- data.frame(
  Oscillibacter=genus_abs_clr[c(PBCbase_resp_6m,PBCbase_noresp_6m),"k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Oscillospiraceae;g__Oscillibacter"],
  Bacteroides=genus_abs_clr[c(PBCbase_resp_6m,PBCbase_noresp_6m),'k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides']
)
genus_resp_pred$Veillonella=genus_abs_clr[c(PBCbase_resp_6m,PBCbase_noresp_6m),"k__Bacteria;p__Firmicutes;c__Negativicutes;o__Veillonellales;f__Veillonellaceae;g__Veillonella"]
row.names(genus_resp_pred)=c(PBCbase_resp_6m,PBCbase_noresp_6m)
genus_resp_pred$response=species_resp_pred[row.names(genus_resp_pred),"response"]
genus_resp_pred$BO_ratio=genus_resp_pred$Bacteroides-genus_resp_pred$Oscillibacter
genus_resp_pred$VO_ratio=genus_resp_pred$Veillonella-genus_resp_pred$Oscillibacter
genus_resp_pred=genus_resp_pred[,c(1,2,4:6,3)]

#determine cutoff and AUC
library("doSNOW")
library("doRNG")
library("tidyverse")
library(cutpointr)
cl <- makeCluster(2)
registerDoSNOW(cl)
set.seed(100)
opt_cut_b_BOratio=cutpointr(genus_resp_pred, BO_ratio, response, boot_runs = 1000,
                            silent = TRUE, allowParallel = TRUE)
stopCluster(cl)

summary(opt_cut_b_BOratio) #AUC=0.6026,optimal cutpoint=0.522

#lasso模型预测
library(caret)
library(glmnet)
set.seed(1234)
randomized_samples = sample(1:nrow(genus_resp_pred))
batches = split(randomized_samples,cut(seq_along(randomized_samples),5,labels = F))


train_X = lapply(batches, function(x) genus_resp_pred[-x,-6])
test_X = lapply(batches, function(x) genus_resp_pred[x,-6])

train_Y = lapply(batches, function(x) genus_resp_pred[-x,6])
test_Y = lapply(batches, function(x) genus_resp_pred[x,6])


## Build prediction models
all_prediction_models = list()
print(' >> training models')
for(i in 1:5) {
  print(paste0('X-val set ',i,' '))
  
  current_X = train_X[[i]]
  current_Y = train_Y[[i]]
  
  all_prediction_models[[paste0("pred.batch",i)]] = cv.glmnet(x = as.matrix(current_X),
                                                              y = current_Y,
                                                              alpha = 1,nfolds = 10,type.measure = "auc",family = "binomial")
}
## calculate training set and test set predictions
# predictions on training set
train.predictions.microb = list()
for (i in 1:5){
  train.predictions.microb[[i]] = 
    # predict(all_prediction_models[[paste0("pred.batch",i,".RData")]],
    #         newx = as.matrix(train_X.microb[[i]]),s="lambda.1se")[,1]
    predict(all_prediction_models[[paste0("pred.batch",i)]],
            newx = as.matrix(train_X[[i]]),s="lambda.min")[,1]
  
}
# predictions on test set
test.predictions.microb = list()
for (i in 1:5){
  test.predictions.microb[[i]] =
    
    predict(all_prediction_models[[paste0("pred.batch",i)]],
            newx = as.matrix(test_X[[i]]),s="lambda.min")[,1]
}
## calculate AUCs 
#microbial prediction
library(pROC)
auc.train.microb=list()
for (i in 1:5) {
  auc.train.microb[[i]]=
    auc(roc(train_Y[[i]],train.predictions.microb[[i]],direction = "<",levels = c(0,1)))
} 

auc.test.microb=list()
for (i in 1:5) {
  auc.test.microb[[i]]=
    auc(roc(test_Y[[i]],test.predictions.microb[[i]],direction = "<",levels = c(0,1)))
}





############################2. PBC HC classification############################
#1.species数据
#读入masslin分析的结果
masslin2_species_PBCHC_res=read.csv("masslin2_species_PBCHC_groupres.csv",header = T,row.names = 1)
row.names(masslin2_species_PBCHC_res)=masslin2_species_PBCHC_res$species
#生成预测的数据（masslin FDR<0.05,n=287）
species_PBCHC_classi=species_abs_clr_scale[c(PBCbase,HC),species_PBCHC_dif_maaslin]
species_PBCHC_classi=cbind(species_PBCHC_classi,metadata[c(PBCbase,HC),"group"])
colnames(species_PBCHC_classi)[288]="group"
#将group列factor转化为0，1因子数据
species_PBCHC_classi$group=as.character(species_PBCHC_classi$group)
species_PBCHC_classi$group[species_PBCHC_classi$group=="PBCbase"]=1
species_PBCHC_classi$group[species_PBCHC_classi$group=="HC"]=0
species_PBCHC_classi$group=as.factor(species_PBCHC_classi$group)

predictor=species_PBCHC_classi[,1:287]
outcome=species_PBCHC_classi[,288]%>% as.data.frame()
row.names(outcome)=row.names(species_PBCHC_classi)

#split PBC 和HC
train.x_PBCHC_classi=row.names(species_PBCHC_classi)[c(1:95,133:233)]
test.x_PBCHC_classi=row.names(species_PBCHC_classi)[c(96:132,234:282)]

#employ train dataset to choose best lambda value and predict in validation dataset
cv_lasso_train_all=cv.glmnet(as.matrix(predictor[train.x_PBCHC_classi,]),as.matrix(outcome[train.x_PBCHC_classi,]),family = "binomial",type.measure = "auc",nfolds = 5,alpha=1)#x,y is matrix not data frame
cv_lasso_train_all #see auc and corresponding lambda in train dataset
plot(cv_lasso_train_all) #绘制纳入变量数目与AUC的关系图
lambdamin=cv_lasso_train_all$lambda.min
lambdalse=cv_lasso_train_all$lambda.1se
cv_lasso_train_all #see auc in train dataset
lasso_train_all=glmnet(as.matrix(predictor[train.x_PBCHC_classi,]),as.matrix(outcome[train.x_PBCHC_classi,]),family = "binomial",type.measure = "auc",lambda = lambdamin)
lasso_train_all$df #看共有多少个feature纳入模型 df=31
predict.y=predict(lasso_train_all,newx = as.matrix(predictor[test.x_PBCHC_classi,]),s=lambdamin)
roc_all=roc(outcome[test.x_PBCHC_classi,],predict.y[,1])
plot(roc_all,print.auc=T,print.thres=F)

#2.血清代谢组数据
serum_PBCHC_classi=serum_0.7_log2UV[c(PBCbase_serum,HC_serum),serum_PBCHC_dif_lm]
serum_PBCHC_classi=cbind(serum_PBCHC_classi,metadata[c(PBCbase_serum,HC_serum),"group"])
colnames(serum_PBCHC_classi)[174]="group"
#因为三丙胺可以单独将PBCHC完全区分，将其删除
serum_PBCHC_classi=subset(serum_PBCHC_classi,select = -c(Tripropylamine_serum))

#将group列factor转化为0，1因子数据
serum_PBCHC_classi$group=as.character(serum_PBCHC_classi$group)
serum_PBCHC_classi$group[serum_PBCHC_classi$group=="PBCbase"]=1
serum_PBCHC_classi$group[serum_PBCHC_classi$group=="HC"]=0
serum_PBCHC_classi$group=as.factor(serum_PBCHC_classi$group)

predictor=serum_PBCHC_classi[,1:172]
outcome=serum_PBCHC_classi[,173]%>% as.data.frame()
row.names(outcome)=row.names(serum_PBCHC_classi)

#split PBC 和HC
train.x_PBCHC_classi=row.names(serum_PBCHC_classi)[c(1:75,112:202)]
test.x_PBCHC_classi=row.names(serum_PBCHC_classi)[c(76:111,203:241)]

#employ train dataset to choose best lambda value and predict in validation dataset
cv_lasso_train_all=cv.glmnet(as.matrix(predictor[train.x_PBCHC_classi,]),as.matrix(outcome[train.x_PBCHC_classi,]),family = "binomial",type.measure = "auc",nfolds = 5,alpha=1)#x,y is matrix not data frame
cv_lasso_train_all #see auc and corresponding lambda in train dataset
plot(cv_lasso_train_all) #绘制纳入变量数目与AUC的关系图
lambdamin=cv_lasso_train_all$lambda.min
lambdalse=cv_lasso_train_all$lambda.1se
cv_lasso_train_all #see auc in train dataset
lasso_train_all=glmnet(as.matrix(predictor[train.x_PBCHC_classi,]),as.matrix(outcome[train.x_PBCHC_classi,]),family = "binomial",type.measure = "auc",lambda = lambdamin)
lasso_train_all$df #看共有多少个feature纳入模型 df=12
lasso_train_all$beta
predict.y=predict(lasso_train_all,newx = as.matrix(predictor[test.x_PBCHC_classi,]),s=lambdamin)
roc_all=roc(outcome[test.x_PBCHC_classi,],predict.y[,1])
plot(roc_all,print.auc=T,print.thres=F)

#3.fecal metabolites
fecal_PBCHC_classi=fecal_0.7_log2UV[c(PBCbase_fecal,HC_fecal),fecal_PBCHC_dif_lm]
fecal_PBCHC_classi=cbind(fecal_PBCHC_classi,metadata[c(PBCbase_fecal,HC_fecal),"group"])
colnames(fecal_PBCHC_classi)[84]="group"
#将group列factor转化为0，1因子数据
fecal_PBCHC_classi$group=as.character(fecal_PBCHC_classi$group)
fecal_PBCHC_classi$group[fecal_PBCHC_classi$group=="PBCbase"]=1
fecal_PBCHC_classi$group[fecal_PBCHC_classi$group=="HC"]=0
fecal_PBCHC_classi$group=as.factor(fecal_PBCHC_classi$group)

predictor=fecal_PBCHC_classi[,1:83]
outcome=fecal_PBCHC_classi[,84]%>% as.data.frame()
row.names(outcome)=row.names(fecal_PBCHC_classi)

#split PBC 和HC
train.x_PBCHC_classi=row.names(fecal_PBCHC_classi)[c(1:90,128:226)]
test.x_PBCHC_classi=row.names(fecal_PBCHC_classi)[c(91:127,227:268)]

#employ train dataset to choose best lambda value and predict in validation dataset
cv_lasso_train_all=cv.glmnet(as.matrix(predictor[train.x_PBCHC_classi,]),as.matrix(outcome[train.x_PBCHC_classi,]),family = "binomial",type.measure = "auc",nfolds = 5,alpha=1)#x,y is matrix not data frame
cv_lasso_train_all #see auc and corresponding lambda in train dataset
plot(cv_lasso_train_all) #绘制纳入变量数目与AUC的关系图
lambdamin=cv_lasso_train_all$lambda.min
lambdalse=cv_lasso_train_all$lambda.1se
cv_lasso_train_all #see auc in train dataset
lasso_train_all=glmnet(as.matrix(predictor[train.x_PBCHC_classi,]),as.matrix(outcome[train.x_PBCHC_classi,]),family = "binomial",type.measure = "auc",lambda = lambdamin)
lasso_train_all$df #看共有多少个feature纳入模型 df=15
predict.y=predict(lasso_train_all,newx = as.matrix(predictor[test.x_PBCHC_classi,]),s=lambdamin)
roc_all=roc(outcome[test.x_PBCHC_classi,],predict.y[,1])
plot(roc_all,print.auc=T,print.thres=F)