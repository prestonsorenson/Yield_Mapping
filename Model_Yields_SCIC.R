library(ranger)
library(raster)
library(Metrics)
library(DescTools)
library(ggplot2)
library(brnn)
library(Cubist)
library(caret)
library(parallel) 
library(SpatialTools)
library(HiClimR)

# Calculate the number of cores
no_cores <- detectCores() - 1

#update as annual cross validation

library(doParallel)
# create the cluster for caret to use
cl <- makePSOCKcluster(no_cores)
registerDoParallel(cl)


setwd('/home/preston/OneDrive/Papers/Yield_Potential/Data')
#load('/home/preston/OneDrive/Papers/Yield_Potential/Analysis/Yield_Potential_2022_05_13.RData')

######functions#######
normalize <- function(x) {
  return ((x - min(na.omit(x))) / (max(na.omit(x)) - min(na.omit(x))))
}


##################set caret values#############
fitControl <- trainControl(## 5-fold CV
  method = "repeatedcv",
  number = 5,
  ## repeated 5 times
  repeats = 1,
  allowParallel = TRUE)

#model hyperparameters to test
cubist_grid=expand.grid(committees=seq(1,20, by=1), neighbors=c(0,1,3,5,7,9))
brnn_grid=expand.grid(neurons=c(1:10))
svm_grid=expand.grid(sigma=0.01,C=c(0.25,0.5,1,2,4,8,16,32,64))
gbm_grid=expand.grid(n.trees=100, interaction.depth=1, shrinkage=c(0.001, 0.01,0.05,0.1, 0.2,0.3,0.4,0.5), n.minobsinnode=1)

###############Process Data####################
#feature selection function
ranger_feature_selection <- function (x, y){
  #remove correlated features
  x.cor=fastCor(x[,-c(1:2)])
  x_sub=x[,-(findCorrelation(x.cor, cutoff=0.9, exact=TRUE))]
  
  model_temp=ranger(dependent.variable.name = y, data=x_sub, importance='impurity')
  
  #forward feature selection
  temp_var=names(sort(importance(model_temp), decreasing=TRUE))
  temp_var=c(y, temp_var)
  
  val=vector('list')
  features=vector('list')
  q=0
  for (i in 2:length(temp_var)){
    q=q+1
    temp=x[,temp_var]
    model=ranger(dependent.variable.name = y, data=temp, importance='impurity')
    val=c(val, model$prediction.error)
    features[[q]]=temp_var
    temp_var=c(y, names(sort(importance(model), decreasing=TRUE))[-length(sort(importance(model), decreasing=TRUE))])
    rm(model)
  }
  val=unlist(val)
  which.min(val)
  
  final_features=unlist(features[which.min(val)])[-1]
  return(final_features)
  }


#import data
yield_data=read.csv('scic_yield_1999_2019_yearly_avg.csv')
colnames(yield_data)[3]='year'

yield_data$id=paste(yield_data$LEGALLAND, yield_data$year, sep='_')

#monthly NDVI
ndvi_med=read.csv('/home/preston/OneDrive/Papers/Yield_Potential/Data/Satellite_Data/ndvi_med_1999_2020_monthly.csv')

#climate data
precip=read.csv('/home/preston/OneDrive/Papers/Yield_Potential/Data/Satellite_Data/precip_1999_2020_monthly.csv')
tmin=read.csv('/home/preston/OneDrive/Papers/Yield_Potential/Data/Satellite_Data/tmin_1999_2020_monthly.csv')
tmax=read.csv('/home/preston/OneDrive/Papers/Yield_Potential/Data/Satellite_Data/tmax_1999_2020_monthly.csv')

ndvi_med=ndvi_med[,-1]
precip=precip[,-1]
tmin=tmin[,-1]
tmax=tmax[,-1]

colnames(ndvi_med)[1]='LEGALLAND'
colnames(precip)[1]='LEGALLAND'
colnames(tmin)[1]='LEGALLAND'
colnames(tmax)[1]='LEGALLAND'

#other_covariates
covariates=read.csv('quarters_yields_predictors_2021_01_06.csv')
colnames(covariates)

covariates=covariates[,-c(1,3:18)]

#MODIS
modis=read.csv('Quarters_MODIS_NDVI_Data_2022_10_05.csv')


#remove non growing seasons for NDVI
ndvi_med=ndvi_med[,-c(2:7,22:25)]

#soil water content
soil_water=read.csv('Quarters_Soil_Moisture_Training_Data_2022_08_17.csv')
soil_water=soil_water[,c(12, 15:21)]

#create spatial grid
qrts=shapefile('/media/preston/My Book/Saskatchewan/SaskGrid_2015_QUARTERSECTION/SaskGrid_2015_QUARTERSECTION_CENTROIDS.shp')
qrts=data.frame(qrts)

x.min=min(qrts$coords.x1)
x.max=max(qrts$coords.x1)
y.min=min(qrts$coords.x2)
y.max=max(qrts$coords.x2)


#create prediction reference points
ref_point1=c(x.min, y.min)
ref_point2=c(x.min, y.max)
ref_point3=c(x.max, y.min)
ref_point4=c(x.max, y.max)
ref_point5=c((x.min+x.max)/2, (y.min+y.max)/2)


ref_points=t(data.frame(ref_point1,ref_point2,ref_point3,ref_point4,ref_point5))
colnames(ref_points)=c('coords.x1', 'coords.x2')

distances=dist2(as.matrix(qrts[,c('coords.x1', 'coords.x2')]), as.matrix(ref_points))
colnames(distances)=c('ref_point1', 'ref_point2', 'ref_point3', 'ref_point4', 'ref_point5')

qrts=data.frame(qrts, distances)

qrts=qrts[,c(10, 13, 14, 16:20)]


#########Hard Red Spring Wheat##################
yield_wheat=yield_data[yield_data$IU.Name=='Wheat - HRSW',]

#wheat model
yield_wheat=na.omit(yield_wheat)

#adjust to 2019 levels
#determine trend and adjustment for each year to normalize to 2019 yields for all years
yields_wheat_avg=aggregate(yield_wheat$kg_ha, by=list(yield_wheat$year), FUN=function (x) mean(na.omit(x)))

#determine annual trend
plot(yields_wheat_avg$x~yields_wheat_avg$Group.1)

colnames(yields_wheat_avg)=c("Year", 'Yield')

#determine trend and adjustment for each year to normalize to 2019 yields for all years
wheat_trend=abs(lm(yields_wheat_avg$Yield~yields_wheat_avg$Year)$fitted.values-max(lm(yields_wheat_avg$Yield~yields_wheat_avg$Year)$fitted.values))

names(wheat_trend)=yields_wheat_avg$Year

#adjust yield values
yield_wheat_adj=data.frame()
for (i in c(1999:2019)){
  temp=yield_wheat[yield_wheat$year==i,]
  temp$kg_ha_adj=temp$kg_ha+wheat_trend[names(wheat_trend)==i]
  yield_wheat_adj=rbind(yield_wheat_adj, temp)
  
}

#average yield per quarter section
yield_wheat_adj=yield_wheat_adj[,c(2,14)]
yield_wheat_adj=aggregate(yield_wheat_adj[,2], by=list(yield_wheat_adj$LEGALLAND), FUN=mean)

colnames(yield_wheat_adj)=c('LEGALLAND', 'kg_ha_adj')

#merge with covariate data
yield_wheat=merge(yield_wheat_adj, ndvi_med, by='LEGALLAND')
yield_wheat=merge(yield_wheat, precip, by='LEGALLAND')
yield_wheat=merge(yield_wheat, tmin, by='LEGALLAND')
yield_wheat=merge(yield_wheat, tmax, by='LEGALLAND')

#merge with covariates
yield_wheat=merge(yield_wheat, covariates, by='LEGALLAND')

#merge with soil water
yield_wheat=merge(yield_wheat, soil_water, by='LEGALLAND')

#merge with MODIS
yield_wheat=merge(yield_wheat, modis, by='LEGALLAND')

colnames(yield_wheat)
yield_wheat=yield_wheat[,-c(94:105)]

#train test splits
sample_size=round(nrow(yield_wheat)*0.5)
train = sample(seq_len(nrow(yield_wheat)),size = sample_size)

wheat_train=yield_wheat[train,]
wheat_test=yield_wheat[-train,]

#remove na and infinite values
wheat_train <- wheat_train[!is.infinite(rowSums(wheat_train[,-1])),]
wheat_train=wheat_train[complete.cases(wheat_train),]

wheat_test <- wheat_test[!is.infinite(rowSums(wheat_test[,-1])),]
wheat_test=wheat_test[complete.cases(wheat_test),]

#feature selection
wheat_features=ranger_feature_selection(wheat_train,'kg_ha_adj')

wheat_features=c('kg_ha_adj', wheat_features)

wheat_train=wheat_train[,wheat_features]

#determine optimal model and hyper parameters

#select optimal model type based on 5-fold CV and R2
cubistFit1 <- train(kg_ha_adj ~ ., data = wheat_train, 
                    method = "cubist", 
                    trControl = fitControl,
                    tuneGrid=cubist_grid,
                    verbose = TRUE)
cubistFit1
cubist_best=cubistFit1$results[which.max(cubistFit1$results$Rsquared),]
cubist_best$type='cubist'

brnnFit1 <-  train(kg_ha_adj ~ ., data = wheat_train, 
                   method = "brnn", 
                   trControl = fitControl,
                   tuneGrid=brnn_grid,
                   verbose = TRUE)

brnnFit1
brnn_best=brnnFit1$results[which.max(brnnFit1$results$Rsquared),]
brnn_best$type='brnn'

rangerFit1 <- train(kg_ha_adj ~ ., data = wheat_train, 
                    method = "ranger", 
                    trControl = fitControl,
                    verbose = TRUE)

rangerFit1
ranger_best=rangerFit1$results[which.max(rangerFit1$results$Rsquared),]
ranger_best$type='ranger'

svmFit1 <-  train(kg_ha_adj~ ., data = wheat_train, 
                  method = "svmRadial", 
                  trControl = fitControl,
                  tuneGrid=svm_grid,
                  verbose = TRUE)

svmFit1
svm_best=svmFit1$results[which.max(svmFit1$results$Rsquared),]
svm_best$type='svmRadial'


gbmFit1 <- train(kg_ha_adj ~ ., data = wheat_train, 
                 method = "gbm", 
                 trControl = fitControl,
                 tuneGrid=gbm_grid,
                 verbose = TRUE)

gbmFit1
gbm_best=gbmFit1$results[which.max(gbmFit1$results$Rsquared),]
gbm_best$type='gbm'


#build model
model_wheat=ranger(kg_ha_adj~., data=wheat_train, importance='impurity')
sort(importance(model_wheat), decreasing=TRUE)

wheat_test=wheat_test[!is.na(wheat_test$MODIS_NDVI_Aug),]
pred_wheat=predict(model_wheat, wheat_test)
pred_wheat=pred_wheat$predictions

summary(lm(pred_wheat~wheat_test$kg_ha_adj))
rmse(wheat_test$kg_ha_adj, pred_wheat)
CCC(wheat_test$kg_ha_adj, pred_wheat)$rho.c
bias(wheat_test$kg_ha_adj, pred_wheat)

#plot results
#create plot
plot_data_wheat=data.frame(wheat_test$kg_ha_adj, pred_wheat)
colnames(plot_data_wheat)=c("actual", 'predicted')
wheat_plot=ggplot(plot_data_wheat, (aes(x=actual, y=predicted))) + geom_point(size=1) + xlim(0, 8000) + ylim(0, 8000) + geom_abline() + xlab(expression(paste('Observed Wheat Yield (kg ha'^"-2",")"))) + ylab(expression(paste('Predicted Wheat Yield (kg ha'^"-2",")"))) + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), title=element_text(size=20),legend.key.size = unit(1, 'cm'), legend.text=element_text(size=20))+ggtitle("Hard Red Spring Wheat")
wheat_plot = wheat_plot + annotate("text", x=500, y = 7900, label=expression(paste("R"^"2 ", "= 0.71"))) + annotate("text", x=500, y=7700, label=expression(paste("RMSE = 539 kg ha"^"-2"))) + annotate("text", x=500, y=7500, label=expression(paste(rho['c'], "= 0.82"))) + annotate("text", x=500, y=7300, label='Bias = -0.02')
wheat_plot

ggsave('/home/preston/OneDrive/Papers/Yield_Potential/Figures/HRSW.png', plot=last_plot(), width=11, height=8.5)

write.csv(sort(importance(model_wheat), decreasing=TRUE), '/home/preston/OneDrive/Papers/Yield_Potential/Analysis/2022_05_12/wheat_feature_importance.csv')
write.csv(plot_data_wheat, '/home/preston/OneDrive/Papers/Yield_Potential/Analysis/2022_05_12/wheat_validation_results_rf.csv')


#################Barley################
yield_barley=yield_data[yield_data$IU.Name=='Barley',]

#merge with covariate data
yield_barley=merge(yield_barley, ndvi_med, by='id')
yield_barley=merge(yield_barley, precip, by='id')
yield_barley=merge(yield_barley, tmin, by='id')
yield_barley=merge(yield_barley, tmax, by='id')


#merge with covariates
yield_sat_barley=merge(yield_barley, covariates, by='LEGALLAND')

#barley model
yield_sat_barley=na.omit(yield_sat_barley)
yield_sat_barley=yield_sat_barley[!is.na(yield_sat_barley$b5),]

#export
write.csv(yield_sat_barley, '/home/preston/OneDrive/Papers/Yield_Potential/Data/scic_barley_covariate_data.csv')

#train-test split
sample_size=round(nrow(yield_sat_barley)*0.5)
train = sample(seq_len(nrow(yield_sat_barley)),size = sample_size)

barley_train=yield_sat_barley[train,]
barley_test=yield_sat_barley[-train,]

#subset barley_train
colnames(barley_train)
barley_train=barley_train[,c(15,18:80)]

#feature selection
barley_features=ranger_feature_selection(barley_train,'kg_ha')

barley_features=c('kg_ha', barley_features)

barley_train=barley_train[,barley_features]

#build model
model_barley=ranger(kg_ha~., data=barley_train, importance='impurity')
sort(importance(model_barley), decreasing=TRUE)

pred_barley=predict(model_barley, barley_test)
pred_barley=pred_barley$predictions

summary(lm(pred_barley~barley_test$kg_ha))
rmse(barley_test$kg_ha, pred_barley)
CCC(barley_test$kg_ha, pred_barley)$rho.c
bias(barley_test$kg_ha, pred_barley)

#plot results
#create plot
plot_data_barley=data.frame(barley_test$kg_ha, pred_barley)
colnames(plot_data_barley)=c("actual", 'predicted')
barley_plot=ggplot(plot_data_barley, (aes(x=actual, y=predicted))) + geom_point(size=1) + xlim(0, 10000) + ylim(0, 10000) + geom_abline() + xlab(expression(paste('Observed Barley Yield (kg ha'^"-2",")"))) + ylab(expression(paste('Predicted Barley Yield (kg ha'^"-2",")"))) + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), title=element_text(size=20),legend.key.size = unit(1, 'cm'), legend.text=element_text(size=20))+ggtitle("Barley")
barley_plot = barley_plot + annotate("text", x=1000, y =10000, label=expression(paste("R"^"2 ", "= 0.63"))) + annotate("text", x=1000, y=9600, label=expression(paste("RMSE = 740 kg ha"^"-2"))) + annotate("text", x=1000, y=9200, label=expression(paste(rho['c'], "= 0.76"))) + annotate("text", x=1000, y=8800, label='Bias = 3.59')
barley_plot

ggsave('/home/preston/OneDrive/Papers/Yield_Potential/Figures/barley.png', plot=last_plot(), width=11, height=8.5)


write.csv(sort(importance(model_barley), decreasing=TRUE), '/home/preston/OneDrive/Papers/Yield_Potential/Analysis/2022_05_12/barley_feature_importance.csv')
write.csv(plot_data_barley, '/home/preston/OneDrive/Papers/Yield_Potential/Analysis/2022_05_12/barley_validation_results_rf.csv')


################Canola#################
yield_canola=yield_data[yield_data$IU.Name=='Can/Rapeseed',]

#merge with covariate data
yield_canola=merge(yield_canola, ndvi_med, by='id')
yield_canola=merge(yield_canola, precip, by='id')
yield_canola=merge(yield_canola, tmin, by='id')
yield_canola=merge(yield_canola, tmax, by='id')

#merge with covariates
yield_sat_canola=merge(yield_canola, covariates, by='LEGALLAND')

#canola model
yield_sat_canola=na.omit(yield_sat_canola)
yield_sat_canola=yield_sat_canola[!is.na(yield_sat_canola$b5),]

#export
write.csv(yield_sat_canola, '/home/preston/OneDrive/Papers/Yield_Potential/Data/scic_canola_covariate_data.csv')

#train-test split
sample_size=round(nrow(yield_sat_canola)*0.5)
train = sample(seq_len(nrow(yield_sat_canola)),size = sample_size)

canola_train=yield_sat_canola[train,]
canola_test=yield_sat_canola[-train,]

#subset canola_train
colnames(canola_train)
canola_train=canola_train[,c(15,18:80)]

#feature selection
canola_features=ranger_feature_selection(canola_train,'kg_ha')

canola_features=c('kg_ha', canola_features)

canola_train=canola_train[,canola_features]


#build model
model_canola=ranger(kg_ha~., data=canola_train, importance='impurity')
sort(importance(model_canola), decreasing=TRUE)

pred_canola=predict(model_canola, canola_test)
pred_canola=pred_canola$predictions

summary(lm(pred_canola~canola_test$kg_ha))
rmse(canola_test$kg_ha, pred_canola)
CCC(canola_test$kg_ha, pred_canola)$rho.c
bias(canola_test$kg_ha, pred_canola)

#plot results
#create plot
plot_data_canola=data.frame(canola_test$kg_ha, pred_canola)
colnames(plot_data_canola)=c("actual", 'predicted')
canola_plot=ggplot(plot_data_canola, (aes(x=actual, y=predicted))) + geom_point(size=1) + xlim(0, 5500) + ylim(0, 5500) + geom_abline() + xlab(expression(paste('Observed Canola Yield (kg ha'^"-2",")"))) + ylab(expression(paste('Predicted Canola Yield (kg ha'^"-2",")"))) + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), title=element_text(size=20),legend.key.size = unit(1, 'cm'), legend.text=element_text(size=20))+ggtitle("Canola")
canola_plot = canola_plot + annotate("text", x=500, y =5500, label=expression(paste("R"^"2 ", "= 0.70"))) + annotate("text", x=500, y=5300, label=expression(paste("RMSE = 379 kg ha"^"-2"))) + annotate("text", x=500, y=5100, label=expression(paste(rho['c'], "= 0.81"))) + annotate("text", x=500, y=4900, label='Bias = 2.88')
canola_plot

ggsave('/home/preston/OneDrive/Papers/Yield_Potential/Figures/canola.png', plot=last_plot(), width=11, height=8.5)

write.csv(sort(importance(model_canola), decreasing=TRUE), '/home/preston/OneDrive/Papers/Yield_Potential/Analysis/2022_05_12/canola_feature_importance.csv')
write.csv(plot_data_canola, '/home/preston/OneDrive/Papers/Yield_Potential/Analysis/2022_05_12/canola_validation_results_rf.csv')


##################Field Peas############
yield_peas=yield_data[yield_data$IU.Name=='Field Peas',]

#merge with covariate data
yield_peas=merge(yield_peas, ndvi_med, by='id')
yield_peas=merge(yield_peas, precip, by='id')
yield_peas=merge(yield_peas, tmin, by='id')
yield_peas=merge(yield_peas, tmax, by='id')

#merge with covariates
yield_sat_peas=merge(yield_peas, covariates, by='LEGALLAND')

#peas model
yield_sat_peas=na.omit(yield_sat_peas)
yield_sat_peas=yield_sat_peas[!is.na(yield_sat_peas$b5),]

#export
write.csv(yield_sat_peas, '/home/preston/OneDrive/Papers/Yield_Potential/Data/scic_peas_covariate_data.csv')

#train-test split
sample_size=round(nrow(yield_sat_peas)*0.5)
train = sample(seq_len(nrow(yield_sat_peas)),size = sample_size)

peas_train=yield_sat_peas[train,]
peas_test=yield_sat_peas[-train,]

#subset peas_train
colnames(peas_train)
peas_train=peas_train[,c(15,18:80)]


#feature selection
peas_features=ranger_feature_selection(peas_train,'kg_ha')

peas_features=c('kg_ha', peas_features)

peas_train=peas_train[,peas_features]


#build model
model_peas=ranger(kg_ha~., data=peas_train, importance='impurity')
sort(importance(model_peas), decreasing=TRUE)

pred_peas=predict(model_peas, peas_test)
pred_peas=pred_peas$predictions

summary(lm(pred_peas~peas_test$kg_ha))
rmse(peas_test$kg_ha, pred_peas)
CCC(peas_test$kg_ha, pred_peas)$rho.c
bias(peas_test$kg_ha, pred_peas)

#plot results
#create plot
plot_data_peas=data.frame(peas_test$kg_ha, pred_peas)
colnames(plot_data_peas)=c("actual", 'predicted')
peas_plot=ggplot(plot_data_peas, (aes(x=actual, y=predicted))) + geom_point(size=1) + xlim(0, 7000) + ylim(0, 7000) + geom_abline() + xlab(expression(paste('Observed Field Peas Yield (kg ha'^"-2",")"))) + ylab(expression(paste('Predicted Field Peas Yield (kg ha'^"-2",")"))) + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), title=element_text(size=20),legend.key.size = unit(1, 'cm'), legend.text=element_text(size=20))+ggtitle("Field Peas")
peas_plot = peas_plot + annotate("text", x=500, y =7000, label=expression(paste("R"^"2 ", "= 0.62"))) + annotate("text", x=500, y=6700, label=expression(paste("RMSE = 570 kg ha"^"-2"))) + annotate("text", x=500, y=6400, label=expression(paste(rho['c'], "= 0.73"))) + annotate("text", x=500, y=6100, label='Bias = 4.18')
peas_plot

ggsave('/home/preston/OneDrive/Papers/Yield_Potential/Figures/peas.png', plot=last_plot(), width=11, height=8.5)

write.csv(sort(importance(model_peas), decreasing=TRUE), '/home/preston/OneDrive/Papers/Yield_Potential/Analysis/2022_05_12/peas_feature_importance.csv')
write.csv(plot_data_peas, '/home/preston/OneDrive/Papers/Yield_Potential/Analysis/2022_05_12/peas_validation_results_rf.csv')









