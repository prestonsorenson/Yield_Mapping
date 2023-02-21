library(raster)
library(prospectr)
library(ranger)
library(caret)
library(Metrics)
library(DescTools)
library(svMisc)

#normalize
normalize <- function(x) {
  return ((x - min(na.omit(x))) / (max(na.omit(x)) - min(na.omit(x))))
}

#feature selection script
psm_ranger_feature_selection <- function (x, y){
  #feature selection
  #select bare soil
  data_temp=y[,c(x,'b1', 'b2', 'b3', 'b4','b5','b7', 'precip', 'temperature')]
  model_temp_bare_soil=ranger(dependent.variable.name = x, data=data_temp, importance='impurity')
  model_temp_bare_soil
  sort(importance(model_temp_bare_soil), decreasing=TRUE)
  
  #forward feature selection
  bare_soil_var=names(sort(importance(model_temp_bare_soil), decreasing=TRUE))
  bare_soil_var=c(x, bare_soil_var)
  
  val=vector('list')
  bs_feature=vector('list')
  q=0
  for (i in 2:length(bare_soil_var)){
    q=q+1
    temp=data_temp[,bare_soil_var]
    model=ranger(dependent.variable.name = x, data=temp, importance='impurity')
    val=c(val, model$prediction.error)
    bs_feature[[q]]=bare_soil_var
    bare_soil_var=c(x, names(sort(importance(model), decreasing=TRUE))[-length(sort(importance(model), decreasing=TRUE))])
    rm(model)
  }
  val_bare_soil=unlist(val)
  which.min(val_bare_soil)
  
  bare_soil_var=unlist(bs_feature[which.min(val)])[-1]
  
  
  
  #band ratios
  data_temp=y[,c(x,'ari','ari_noBareSoil','CRSI','CRSI_noBareSoil','NDVI_July_Aug','NDVI_July_Aug_noBareSoil','NDVI_Sept_Oct','NDVI_Sept_Oct_noBareSoil','NDVI_SD','SAVI_Jul_Aug','SAVI_Jul_Aug_noBareSoil','precip','temperature')]
  model_temp_band_ratios=ranger(dependent.variable.name = x, data=data_temp, importance='impurity')
  model_temp_band_ratios
  sort(importance(model_temp_band_ratios), decreasing=TRUE)
  
  #forward feature selection
  band_ratios_var=names(sort(importance(model_temp_band_ratios), decreasing=TRUE))
  band_ratios_var=c(x, band_ratios_var)
  
  
  val=vector('list')
  band_ratios_feature=vector('list')
  q=0
  for (i in 2:length(band_ratios_var)){
    q=q+1
    temp=data_temp[,band_ratios_var]
    model=ranger(dependent.variable.name = x, data=temp, importance='impurity')
    val=c(val, model$prediction.error)
    band_ratios_feature[[q]]=band_ratios_var
    band_ratios_var=c(x, names(sort(importance(model), decreasing=TRUE))[-length(sort(importance(model), decreasing=TRUE))])
    rm(model)
  }
  val_band_ratios=unlist(val)
  which.min(val_band_ratios)
  
  band_ratios_var=unlist(band_ratios_feature[which.min(val)])[-1]
  
  
  #terrain attributes
  data_temp=y[,c(x,'precip', 'temperature','dem_3x3_sd3x3','dem_3x3_sd5x5','dem_3x3_sd9x9','dem_3x3_sd21x21','dem_9x9_sd21x21','dem_9x9_sd101x101','dem_3x3_tri','dem_5x5_tri','dem_9x9_tri','dem_9x9_tri_20', 'MidSlope_Pos_100m', 'Norm_Height_100m', 'Slope_Height_100m', 'Stand_Height_100m', 'SWI_100m', 'Valley_Depth_100m')]
  model_temp_terrain=ranger(dependent.variable.name = x, data=data_temp, importance='impurity')
  model_temp_terrain
  sort(importance(model_temp_terrain), decreasing=TRUE)
  
  #forward feature selection
  terrain_var=names(sort(importance(model_temp_terrain), decreasing=TRUE))
  terrain_var=c(x, terrain_var)
  
  val=vector('list')
  terrain_feature=vector('list')
  q=0
  for (i in 2:length(terrain_var)){
    q=q+1
    temp=data_temp[,terrain_var]
    model=ranger(dependent.variable.name = x, data=temp, importance='impurity')
    val=c(val, model$prediction.error)
    terrain_feature[[q]]=terrain_var
    terrain_var=c(x, names(sort(importance(model), decreasing=TRUE))[-length(sort(importance(model), decreasing=TRUE))])
    rm(model)
  }
  val_terrain=unlist(val)
  which.min(val_terrain)
  
  terrain_var=unlist(terrain_feature[which.min(val)])[-1]
  
  #final features
  temp_init_var=c(x, bare_soil_var, band_ratios_var, terrain_var)
  temp_init_var=temp_init_var[!duplicated(temp_init_var)]
  
  var_train_a_init=y[,temp_init_var]
  
  #remove correlated features
  var_cor=findCorrelation(cor(var_train_a_init[,-1]), cutoff=0.9)
  var_cor
  var_cor=var_cor+1
  
  var_train_a_init=var_train_a_init[,-var_cor]
  
  if(ncol(var_train_a_init)==0){
    var_train_a_init=y[,temp_init_var]
  }
  
  model_temp_final=ranger(dependent.variable.name = x, data=var_train_a_init, importance='impurity')
  model_temp_final
  sort(importance(model_temp_final), decreasing=TRUE)
  
  #forward feature selection
  final_var=names(sort(importance(model_temp_final), decreasing=TRUE))
  final_var=c(x, final_var)
  val=vector('list')
  final_feature=vector('list')
  q=0
  for (i in 2:length(final_var)){
    q=q+1
    temp=y[,final_var]
    model=ranger(dependent.variable.name = x, data=temp, importance='impurity')
    val=c(val, model$prediction.error)
    final_feature[[q]]=final_var
    final_var=c(x, names(sort(importance(model), decreasing=TRUE))[-length(sort(importance(model), decreasing=TRUE))])
    rm(model)
  }
  val_final=unlist(val)
  which.min(val_final)
  
  final_var=unlist(final_feature[which.min(val)])[-1]
  final_var=c(x, final_var)
  
  return(final_var)
}

#NPDB Data
soil_points=shapefile('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Data/Sask_NPDB_Chemical_Physical_Data.shp')
raster_values=read.csv('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Data/Historical_Soil_Properties_Training_Data_2021_10_01.csv')

#Calculate Weighted Averages for Soil Properties
soil_points$thickness=abs(as.numeric(soil_points$Lower_Dept)-as.numeric(soil_points$Upper_Dept))
soil_points=data.frame(soil_points)
soil_points_A=soil_points[soil_points$Master_Hor=='A',]
soil_points_B=soil_points[soil_points$Master_Hor=='B',]
soil_points_C=soil_points[soil_points$Master_Hor=='C',]

soil_points_solum=soil_points[soil_points$Master_Hor %in% c('A', 'B'),]

#soc
soc_a=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_A$PEDON_ID)){
  temp=soil_points_A[soil_points_A$PEDON_ID==i,]
  temp=temp[temp$Organic_Ca!="NA",]
  soc=round(sum(as.numeric(temp$Organic_Ca)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, soc))
  temp=t(data.frame(temp))
  soc_a=rbind(soc_a, temp)
}
colnames(soc_a)=c("PEDON_ID", 'soc')

#create training files
raster_values=raster_values[,c(3,54:91)]
raster_values=apply(raster_values, 2, FUN = function(x) as.numeric(x))
raster_values=data.frame(raster_values)

#merge with raster values
soc_a=merge(soc_a,raster_values, by='PEDON_ID')
soc_a=soc_a[,-1]

##########Sask OM Survey Data########
sask_om=read.csv('/home/preston/OneDrive/Papers/Saskatchewan_OM_Mapping/Data/Sask_OM_Training_Data_2022_09_13.csv')
sask_om=sask_om[,c(12,19:56)]

#convert OM to OC
sask_om$OM=sask_om$OM*0.58
colnames(sask_om)[1]='soc'

#combine datasets
soc_a=rbind(soc_a, sask_om)

#replace na with 0
soc_a[is.na(soc_a)] <- 0

###########create train test splits############
train_test_split=kenStone(soc_a[,-1], k=3792, pc=0.99, .center=TRUE, .scale=TRUE, metric='mahal')

soc_a=soc_a[soc_a$soc<7,]

train=soc_a[train_test_split$model,]
test=soc_a[train_test_split$test,]

train=train[train$soc<7,]
test=test[test$soc<7,]

#########Feature Selection##############
features=psm_ranger_feature_selection('soc', train)
features=c(features, "b5", "b7")

#########Build Model###########
train=train[,features]

model=ranger(soc~., data=train, importance='impurity')

pred=predict(model, data=test)
pred=pred$predictions

performace=data.frame(summary(lm(pred~test$soc))$r.squared,
           rmse(test$soc, pred),
           CCC(test$soc, pred)$rho.c,
           bias(test$soc, pred))

performace

setwd('/home/preston/OneDrive/Papers/Saskatchewan_OM_Mapping/Analysis/')
export_name=paste('soc_mapping_results_', Sys.Date(), '.csv', sep="")
write.csv(performace, export_name)

#export feature importance
imp_val=sort(importance(model), decreasing=TRUE)
imp_val_norm=sort(importance(model), decreasing=TRUE)/sum(sort(importance(model), decreasing=TRUE))
export_name=paste('soc_raw_importance_values_', Sys.Date(), '.csv', sep="")
write.csv(imp_val, export_name)

export_name=paste('soc_importance_values_', Sys.Date(), '.csv', sep="")
write.csv(imp_val_norm, export_name)

#create plot
plot_data=data.frame(test$soc, pred)
colnames(plot_data)=c("actual", 'predicted')
export_name=paste('soc_performance_', Sys.Date(), '.png', sep="")
png(file=export_name, width=800, height=600)
soc_plot=ggplot(plot_data, (aes(x=actual, y=predicted))) + geom_point(size=2) + xlim(0, 7) + ylim(0, 7) + geom_abline() + xlab("Observed Soil Organic Carbon (%)") + ylab("Predicted Soil Organic Carbon (%)") + theme(text = element_text(size = 20))
soc_plot
dev.off()


#############Make full province map################
bare_soil=stack('/media/preston/My Book/Saskatchewan/sk_bare_soil/sk_ls5_bare_soil_ndvi3_ndsi0_nbr1_focal10_filt/sk_ls5_bare_soil_ndsi_ndvi3_ndsi0_nbr1_focal10.tif')
b1=bare_soil[[1]]
b2=bare_soil[[2]]
b3=bare_soil[[3]]
b4=bare_soil[[4]]
b5=bare_soil[[5]]
b6=bare_soil[[6]]
b7=bare_soil[[7]]
ari=raster('/media/preston/My Book/Saskatchewan/sk_l5_ARI_median_focal10_100m/sk_l5_ARI_median_focal10_100m.tif')
ari_noBareSoil=raster('/media/preston/My Book/Saskatchewan/sk_l5_ARI__noBareSoil_median_focal10_100m/sk_l5_ARI_noBareSoil_median_focal10_100m.tif')
CRSI=raster('/media/preston/My Book/Saskatchewan/sk_l5_CSRI_median_focal10_100m/sk_l5_CSRI_median_focal10_100m.tif')
CRSI_noBareSoil=raster('/media/preston/My Book/Saskatchewan/sk_l5_CSRI_noBareS0il_median_focal10_100m/sk_l5_CSRI_noBareS0il_median_focal10_100m.tif')
NDVI_July_Aug=raster('/media/preston/My Book/Saskatchewan/sk_l5_ndvi_median_july_aug_focal10_100m/sk_l5_ndvi_median_july_aug_focal10_100m.tif')
NDVI_July_Aug_noBareSoil=raster('/media/preston/My Book/Saskatchewan/sk_l5_ndvi_noBareSoil_median_july_aug_focal10_100m/sk_l5_ndvi_noBareSoil_median_july_aug_focal10_100m.tif')
NDVI_Sept_Oct=raster('/media/preston/My Book/Saskatchewan/sk_l5_ndvi_median_sept_oct_focal10_100m/sk_l5_ndvi_median_sept_oct_focal10_100m.tif')
NDVI_Sept_Oct_noBareSoil=raster('/media/preston/My Book/Saskatchewan/sk_l5_ndvi_noBareSoil_median_sept_oct_focal10_100m/sk_l5_ndvi_noBareSoil_median_sept_oct_focal10_100m.tif')
NDVI_SD=raster('/media/preston/My Book/Saskatchewan/sk_l5_ndvi_sd_focal10_100m/sk_l5_ndvi_sd_focal10_100m.tif')
SAVI_Jul_Aug=raster('/media/preston/My Book/Saskatchewan/sk_l5_SAVI_median_july_aug_focal10_100m/sk_l5_SAVI_median_july_aug_focal10_100m.tif')
SAVI_Jul_Aug_noBareSoil=raster('/media/preston/My Book/Saskatchewan/sk_l5_SAVI_noBareSoil_median_july_aug_focal10_100m/sk_l5_SAVI_noBareSoil_median_july_aug_focal10_100m.tif')
SAVI_Sept_Oct=raster('/media/preston/My Book/Saskatchewan/sk_l5_SAVI_median_sept_oct_focal10_100m/sk_l5_SAVI_median_sept_oct_focal10_100m.tif')
SAVI_Sept_Oct_noBareSoil=raster('/media/preston/My Book/Saskatchewan/sk_l5_SAVI_noBareSoil_median_sept_oct_focal10_100m/sk_l5_SAVI_noBareSoil_median_sept_oct_focal10_100m.tif')
dem_3x3_sd3x3=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_3x3_sd3x3.tif')
dem_3x3_sd5x5=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_3x3_sd5x5.tif')
dem_3x3_sd9x9=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_3x3_sd9x9.tif')
dem_3x3_sd21x21=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_3x3_sd21x21.tif')
dem_9x9_sd21x21=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_9x9_sd21x21.tif')
dem_9x9_sd101x101=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_9x9_sd101x101.tif')
dem_3x3_tri=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_3x3_tri.tif') 
dem_5x5_tri=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_5x5_tri.tif')
dem_9x9_tri=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_9x9_tri.tif')
dem_9x9_tri_20=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_9x9_tri_20.tif')
precip=raster('/media/preston/My Book/Saskatchewan/Saskatchewan_Climate/sk_era5_tp_median_100m_bicubic.tif')
temperature=raster('/media/preston/My Book/Saskatchewan/Saskatchewan_Climate/sk_era5_2mt_median_100m_bicubic.tif') 
MidSlope_Pos_100m=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/Terrain_Derivatives/MidSlope_Pos_100m.tif')
Norm_Height_100m=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/Terrain_Derivatives/Norm_Height_100m.tif') 
Slope_Height_100m=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/Terrain_Derivatives/Slope_Height_100m.tif')
Stand_Height_100m=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/Terrain_Derivatives/Stand_Height_100m.tif') 
SWI_100m=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/Terrain_Derivatives/SWI_100m.tif')
Valley_Depth_100m=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/Terrain_Derivatives/Valley_Depth_100m.tif')

NDWI=readAll(raster('/media/preston/My Book/Saskatchewan/sk_sen2_ndwi_med_july/sk_sen2_ndwi_med_july_100m.tif'))
NDWI@data@values[NDWI@data@values>0] <- 1

raster_stack=stack(b1, b2, b3, b4, b5, b6, b7, ari, ari_noBareSoil, CRSI, CRSI_noBareSoil, NDVI_July_Aug, NDVI_July_Aug_noBareSoil, NDVI_Sept_Oct, NDVI_Sept_Oct_noBareSoil, NDVI_SD, SAVI_Jul_Aug, SAVI_Jul_Aug_noBareSoil, SAVI_Sept_Oct, SAVI_Sept_Oct_noBareSoil, dem_3x3_sd3x3, dem_3x3_sd5x5, dem_3x3_sd9x9, dem_3x3_sd21x21, dem_9x9_sd21x21, dem_9x9_sd101x101, dem_3x3_tri, dem_5x5_tri, dem_9x9_tri, dem_9x9_tri_20, precip, temperature, MidSlope_Pos_100m, Norm_Height_100m, Slope_Height_100m, Stand_Height_100m, SWI_100m, Valley_Depth_100m)

setwd('/home/preston/OneDrive/Papers/Saskatchewan_OM_Mapping/Analysis')

full_model_data=soc_a[,features]

model_full=ranger(soc~., data=full_model_data,importance='impurity', quantreg = TRUE)

tiles=shapefile("/home/preston/OneDrive/Papers/Historic_Soil_Properties/Analysis/analysis_tiles.shp")

setwd('/home/preston/OneDrive/Papers/Saskatchewan_OM_Mapping/Analysis/Maps')
for (i in 1:length(tiles)){
  progress(i, max.value=length(tiles))
  tile_sub=tiles[tiles$id==i,]
  raster_sub=crop(raster_stack, tile_sub)
  ndwi_sub=crop(NDWI, tile_sub)
  raster_sub=rasterToPoints(raster_sub)
  xy=raster_sub[,1:2]
  raster_sub=raster_sub[,-c(1:2)]
  raster_sub=data.frame(raster_sub)
  colnames(raster_sub)=c('b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'ari', 'ari_noBareSoil', 'CRSI', 'CRSI_noBareSoil', 'NDVI_July_Aug', 'NDVI_July_Aug_noBareSoil', 'NDVI_Sept_Oct', 'NDVI_Sept_Oct_noBareSoil', 'NDVI_SD', 'SAVI_Jul_Aug', 'SAVI_Jul_Aug_noBareSoil', 'SAVI_Sept_Oct', 'SAVI_Sept_Oct_noBareSoil', 'dem_3x3_sd3x3', 'dem_3x3_sd5x5', 'dem_3x3_sd9x9', 'dem_3x3_sd21x21', 'dem_9x9_sd21x21', 'dem_9x9_sd101x101', 'dem_3x3_tri', 'dem_5x5_tri', 'dem_9x9_tri','dem_9x9_tri_20', 'precip', 'temperature', 'MidSlope_Pos_100m', 'Norm_Height_100m', 'Slope_Height_100m', 'Stand_Height_100m', 'SWI_100m', 'Valley_Depth_100m')
  raster_sub[is.na(raster_sub)] <- 0
  raster_sub=do.call(data.frame, lapply(raster_sub, function(x) replace(x, is.infinite(x), 0)))
  pred_sub=predict(model_full, data=raster_sub, type='quantiles', quantiles=c(0.25, 0.5, 0.75))
  pred_sub=pred_sub$predictions
  pred_sub=data.frame(pred_sub)
  colnames(pred_sub)=c("pred_25", "pred_50", "pred_75")
  pred_sub$iqr=pred_sub$pred_75-pred_sub$pred_25
  pred_sub=data.frame(pred_sub, xy)
  
  pred_results_25=pred_sub[,c(1,5:6)]
  pred_results_50=pred_sub[,c(2,5:6)]
  pred_results_75=pred_sub[,c(3,5:6)]
  pred_results_iqr=pred_sub[,c(4,5:6)]
  pred_results_iqr_ratio=pred_sub[,c(4,5:6)]
  pred_results_iqr_ratio$iqr=normalize(pred_results_iqr_ratio$iqr)
  
  coordinates(pred_results_25)=~x+y
  coordinates(pred_results_50)=~x+y
  coordinates(pred_results_75)=~x+y
  coordinates(pred_results_iqr)=~x+y
  coordinates(pred_results_iqr_ratio)=~x+y
  
  pred_results_25=rasterFromXYZ(pred_results_25)
  pred_results_50=rasterFromXYZ(pred_results_50)
  pred_results_75=rasterFromXYZ(pred_results_75)
  pred_results_iqr=rasterFromXYZ(pred_results_iqr)
  pred_results_iqr_ratio=rasterFromXYZ(pred_results_iqr_ratio)
  
  pred_results_25=mask(pred_results_25, ndwi_sub, maskvalue=1)
  pred_results_50=mask(pred_results_50, ndwi_sub, maskvalue=1)
  pred_results_75=mask(pred_results_75, ndwi_sub, maskvalue=1)
  pred_results_iqr=mask(pred_results_iqr, ndwi_sub, maskvalue=1)
  pred_results_iqr_ratio=mask(pred_results_iqr_ratio, ndwi_sub, maskvalue=1)
  
  pred_stack=stack(pred_results_25, pred_results_50, pred_results_75, pred_results_iqr, pred_results_iqr_ratio)
  
  names(pred_stack)=c("pred_results_25","pred_results_50","pred_results_75","pred_results_iqr","pred_results_iqr_ratio")
  
  figure_name=paste("sk_predicted_soc_",i, sep='')
  crs(pred_stack)=crs(ari)
  writeRaster(pred_stack, figure_name, format="GTiff", overwrite=TRUE)
}
