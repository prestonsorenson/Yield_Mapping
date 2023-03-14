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
library(ggplot2)
library(effects)
library(lme4)
library(doParallel)
library(gridExtra)
cl=makeCluster(8)
registerDoParallel(cl)


#normalize function
normalize <- function(x) {
  return ((x - min(na.omit(x))) / (max(na.omit(x)) - min(na.omit(x))))
}


# Calculate the number of cores
no_cores <- detectCores() - 1

#update as annual cross validation

library(doParallel)
# create the cluster for caret to use
cl <- makePSOCKcluster(no_cores)
registerDoParallel(cl)


setwd('/home/preston/OneDrive/Papers/Yield_Potential/Data')
load('/home/preston/OneDrive/Papers/Yield_Potential/Analysis/Yield_Potential_2022_10_26.RData')

#import data
yield_data=read.csv('scic_yield_1999_2019_yearly_avg.csv')
colnames(yield_data)[3]='year'

yield_data$id=paste(yield_data$LEGALLAND, yield_data$year, sep='_')

#soil properties
soil=read.csv('/home/preston/OneDrive/Projects/Grassland_Carbon_Sequestration/Site_Selection/Data/quater_section_soil_land_use.csv')


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
yield_wheat_adj=yield_wheat_adj[,c(2,18)]
yield_wheat_adj=aggregate(yield_wheat_adj[,2], by=list(yield_wheat_adj$LEGALLAND), FUN=mean)

colnames(yield_wheat_adj)=c('LEGALLAND', 'kg_ha_adj')

#merge with soil properties
yield_wheat=merge(yield_wheat, soil, by='LEGALLAND')

#remove na
yield_wheat=yield_wheat[!is.na(yield_wheat$kg_ha_adj),]
yield_wheat=yield_wheat[!is.na(yield_wheat$soc),]
yield_wheat=yield_wheat[!is.na(yield_wheat$clay),]

#remove duplicates
yield_wheat=yield_wheat[!duplicated(yield_wheat$LEGALLAND),]

#########Barley##################
yield_barley=yield_data[yield_data$IU.Name=='Barley',]

#barley model
yield_barley=na.omit(yield_barley)

#adjust to 2019 levels
#determine trend and adjustment for each year to normalize to 2019 yields for all years
yields_barley_avg=aggregate(yield_barley$kg_ha, by=list(yield_barley$year), FUN=function (x) mean(na.omit(x)))

#determine annual trend
plot(yields_barley_avg$x~yields_barley_avg$Group.1)

colnames(yields_barley_avg)=c("Year", 'Yield')

#determine trend and adjustment for each year to normalize to 2019 yields for all years
barley_trend=abs(lm(yields_barley_avg$Yield~yields_barley_avg$Year)$fitted.values-max(lm(yields_barley_avg$Yield~yields_barley_avg$Year)$fitted.values))

names(barley_trend)=yields_barley_avg$Year

#adjust yield values
yield_barley_adj=data.frame()
for (i in c(1999:2019)){
  temp=yield_barley[yield_barley$year==i,]
  temp$kg_ha_adj=temp$kg_ha+barley_trend[names(barley_trend)==i]
  yield_barley_adj=rbind(yield_barley_adj, temp)
  
}

#average yield per quarter section
yield_barley_adj=yield_barley_adj[,c(2,18)]
yield_barley_adj=aggregate(yield_barley_adj[,2], by=list(yield_barley_adj$LEGALLAND), FUN=mean)

colnames(yield_barley_adj)=c('LEGALLAND', 'kg_ha_adj')

#merge with soil properties
yield_barley=merge(yield_barley, soil, by='LEGALLAND')

#remove na
yield_barley=yield_barley[!is.na(yield_barley$kg_ha_adj),]
yield_barley=yield_barley[!is.na(yield_barley$soc),]
yield_barley=yield_barley[!is.na(yield_barley$clay),]

yield_barley=yield_barley[!duplicated(yield_barley$LEGALLAND),]

#########canola##################
yield_canola=yield_data[yield_data$IU.Name=='Can/Rapeseed',]

#canola model
yield_canola=na.omit(yield_canola)

#adjust to 2019 levels
#determine trend and adjustment for each year to normalize to 2019 yields for all years
yields_canola_avg=aggregate(yield_canola$kg_ha, by=list(yield_canola$year), FUN=function (x) mean(na.omit(x)))

#determine annual trend
plot(yields_canola_avg$x~yields_canola_avg$Group.1)

colnames(yields_canola_avg)=c("Year", 'Yield')

#determine trend and adjustment for each year to normalize to 2019 yields for all years
canola_trend=abs(lm(yields_canola_avg$Yield~yields_canola_avg$Year)$fitted.values-max(lm(yields_canola_avg$Yield~yields_canola_avg$Year)$fitted.values))

names(canola_trend)=yields_canola_avg$Year

#adjust yield values
yield_canola_adj=data.frame()
for (i in c(1999:2019)){
  temp=yield_canola[yield_canola$year==i,]
  temp$kg_ha_adj=temp$kg_ha+canola_trend[names(canola_trend)==i]
  yield_canola_adj=rbind(yield_canola_adj, temp)
  
}

#average yield per quarter section
yield_canola_adj=yield_canola_adj[,c(2,18)]
yield_canola_adj=aggregate(yield_canola_adj[,2], by=list(yield_canola_adj$LEGALLAND), FUN=mean)

colnames(yield_canola_adj)=c('LEGALLAND', 'kg_ha_adj')


#merge with soil properties
yield_canola=merge(yield_canola, soil, by='LEGALLAND')

#remove na
yield_canola=yield_canola[!is.na(yield_canola$kg_ha_adj),]
yield_canola=yield_canola[!is.na(yield_canola$soc),]
yield_canola=yield_canola[!is.na(yield_canola$clay),]

yield_canola=yield_canola[!duplicated(yield_canola$LEGALLAND),]

#########peas##################
yield_peas=yield_data[yield_data$IU.Name=='Field Peas',]

#peas model
yield_peas=na.omit(yield_peas)

#adjust to 2019 levels
#determine trend and adjustment for each year to normalize to 2019 yields for all years
yields_peas_avg=aggregate(yield_peas$kg_ha, by=list(yield_peas$year), FUN=function (x) mean(na.omit(x)))

#determine annual trend
plot(yields_peas_avg$x~yields_peas_avg$Group.1)

colnames(yields_peas_avg)=c("Year", 'Yield')

#determine trend and adjustment for each year to normalize to 2019 yields for all years
peas_trend=abs(lm(yields_peas_avg$Yield~yields_peas_avg$Year)$fitted.values-max(lm(yields_peas_avg$Yield~yields_peas_avg$Year)$fitted.values))

names(peas_trend)=yields_peas_avg$Year

#adjust yield values
yield_peas_adj=data.frame()
for (i in c(1999:2019)){
  temp=yield_peas[yield_peas$year==i,]
  temp$kg_ha_adj=temp$kg_ha+peas_trend[names(peas_trend)==i]
  yield_peas_adj=rbind(yield_peas_adj, temp)
  
}

#average yield per quarter section
yield_peas_adj=yield_peas_adj[,c(2,18)]
yield_peas_adj=aggregate(yield_peas_adj[,2], by=list(yield_peas_adj$LEGALLAND), FUN=mean)

colnames(yield_peas_adj)=c('LEGALLAND', 'kg_ha_adj')

#merge with soil properties
yield_peas=merge(yield_peas, soil, by='LEGALLAND')

#remove na
yield_peas=yield_peas[!is.na(yield_peas$kg_ha_adj),]
yield_peas=yield_peas[!is.na(yield_peas$soc),]
yield_peas=yield_peas[!is.na(yield_peas$clay),]

yield_peas=yield_peas[!duplicated(yield_peas$LEGALLAND),]

#######################Normalize Yield and SOC within a RM to remove climate effects##################
#add in RM
rm_file=shapefile('/media/preston/My Book/Saskatchewan/SaskGrid_2015_QUARTERSECTION/Quarters_RM_Label/RM_SaskGrid_2015_QUARTERSECTION.shp')
rm_file=rm_file@data[,c(10, 16)]
colnames(rm_file)[1]='LEGALLAND'

rm_file=rm_file[!duplicated(rm_file$LEGALLAND),]

#wheat
yield_wheat=merge(rm_file, yield_wheat, by='LEGALLAND')
yield_wheat=yield_wheat[!duplicated(yield_wheat$LEGALLAND),]

#wheat
yield_wheat_norm=data.frame()
for (i in unique(yield_wheat$RMNO)){
  dat_sub=yield_wheat[yield_wheat$RMNO==i,]
  dat_sub$soc_norm=normalize(dat_sub$soc)
  dat_sub$kg_ha_adj_norm=normalize(dat_sub$kg_ha_adj)
  yield_wheat_norm=rbind(yield_wheat_norm, dat_sub)
}
yield_wheat_norm=yield_wheat_norm[!duplicated(yield_wheat_norm$LEGALLAND),]
yield_wheat_norm=yield_wheat_norm[!is.na(yield_wheat_norm$soc_norm),]

#lme with nlme
model_wheat_n=lmer(kg_ha_adj_norm~soc_norm + (soc_norm||RMNO),data=yield_wheat_norm, control = lmerControl(check.conv.grad= .makeCC("warning", tol = 1e-2, relTol = NULL)))
anova(model_wheat_n)
summary(model_wheat_n)

wheat_effects_n=effects::effect(term='soc_norm', mod=model_wheat_n)
wheat_effects_n=as.data.frame(wheat_effects_n)

#barley
yield_barley=merge(rm_file, yield_barley, by='LEGALLAND')
yield_barley=yield_barley[!duplicated(yield_barley$LEGALLAND),]

#barley
yield_barley_norm=data.frame()
for (i in unique(yield_barley$RMNO)){
  dat_sub=yield_barley[yield_barley$RMNO==i,]
  dat_sub$soc_norm=normalize(dat_sub$soc)
  dat_sub$kg_ha_adj_norm=normalize(dat_sub$kg_ha_adj)
  yield_barley_norm=rbind(yield_barley_norm, dat_sub)
}
yield_barley_norm=yield_barley_norm[!duplicated(yield_barley_norm$LEGALLAND),]
yield_barley_norm=yield_barley_norm[!is.na(yield_barley_norm$soc_norm),]

#lme with nlme
model_barley_n=lmer(kg_ha_adj_norm~soc_norm + (soc_norm||RMNO),data=yield_barley_norm, control = lmerControl(check.conv.grad= .makeCC("warning", tol = 1e-2, relTol = NULL)))
anova(model_barley_n)
summary(model_barley_n)

barley_effects_n=effects::effect(term='soc_norm', mod=model_barley_n)
barley_effects_n=as.data.frame(barley_effects_n)


#canola
yield_canola=merge(rm_file, yield_canola, by='LEGALLAND')
yield_canola=yield_canola[!duplicated(yield_canola$LEGALLAND),]

#canola
yield_canola_norm=data.frame()
for (i in unique(yield_canola$RMNO)){
  dat_sub=yield_canola[yield_canola$RMNO==i,]
  dat_sub$soc_norm=normalize(dat_sub$soc)
  dat_sub$kg_ha_adj_norm=normalize(dat_sub$kg_ha_adj)
  yield_canola_norm=rbind(yield_canola_norm, dat_sub)
}
yield_canola_norm=yield_canola_norm[!duplicated(yield_canola_norm$LEGALLAND),]
yield_canola_norm=yield_canola_norm[!is.na(yield_canola_norm$soc_norm),]

#lme with nlme
model_canola_n=lmer(kg_ha_adj_norm~soc_norm + (soc_norm||RMNO),data=yield_canola_norm, control = lmerControl(check.conv.grad= .makeCC("warning", tol = 1e-2, relTol = NULL)))
anova(model_canola_n)
summary(model_canola_n)

canola_effects_n=effects::effect(term='soc_norm', mod=model_canola_n)
canola_effects_n=as.data.frame(canola_effects_n)

#peas
yield_peas=merge(rm_file, yield_peas, by='LEGALLAND')
yield_peas=yield_peas[!duplicated(yield_peas$LEGALLAND),]

yield_peas_norm=data.frame()
for (i in unique(yield_peas$RMNO)){
  dat_sub=yield_peas[yield_peas$RMNO==i,]
  dat_sub$soc_norm=normalize(dat_sub$soc)
  dat_sub$kg_ha_adj_norm=normalize(dat_sub$kg_ha_adj)
  yield_peas_norm=rbind(yield_peas_norm, dat_sub)
}
yield_peas_norm=yield_peas_norm[!duplicated(yield_peas_norm$LEGALLAND),]
yield_peas_norm=yield_peas_norm[!is.na(yield_peas_norm$soc_norm),]

#lme with nlme
model_peas_n=lmer(kg_ha_adj_norm~soc_norm + (soc_norm||RMNO),data=yield_peas_norm, control = lmerControl(check.conv.grad= .makeCC("warning", tol = 1e-2, relTol = NULL)))
anova(model_peas_n)
summary(model_peas_n)

peas_effects_n=effects::effect(term='soc_norm', mod=model_peas_n)
peas_effects_n=as.data.frame(peas_effects_n)



################Plots###################
wheat_effects[1,1]=0.9

#absolute values
cols=c('Wheat'='red', 'Barley'='brown', 'Canola'='darkgoldenrod1', 'Peas'='blue')

fixed_effects_plot <- ggplot() + xlim(c(0.8,6)) + ylim(1500, 4000) + 
  geom_line(data=wheat_effects, aes(x=soc, y=fit, colour='Wheat'), linetype=4) +
  geom_ribbon(data= wheat_effects, aes(x=soc, ymin=lower, ymax=upper), alpha= 0.3, fill="red") +
  geom_line(data=barley_effects, aes(x=soc, y=fit, colour='Barley'), linetype=2) +
  geom_ribbon(data= barley_effects, aes(x=soc, ymin=lower, ymax=upper), alpha= 0.3, fill="brown") +
  geom_line(data=canola_effects, aes(x=soc, y=fit, colour='Canola'), linetype=1) +
  geom_ribbon(data= canola_effects, aes(x=soc, ymin=lower, ymax=upper), alpha= 0.3, fill="darkgoldenrod1") +
  geom_line(data=peas_effects, aes(x=soc, y=fit, colour='Peas'), linetype=3) +
  geom_ribbon(data= peas_effects, aes(x=soc, ymin=lower, ymax=upper), alpha= 0.3, fill="blue") +
  labs(x="Soil Organic Carbon (%)", y=expression(paste("Yield (kg ha"^"-1",")"))) +
  scale_linetype_manual(values = c("Wheat" = 1, "Barley" = 2, "Canola" = 3, 'Peas' = 4)) +
  scale_colour_manual(name="Crop Types",values=cols,guide = guide_legend(override.aes = list(linetype = c(4,2,1,3))))+
  theme(text = element_text(size = 20))

fixed_effects_plot

setwd('/home/preston/OneDrive/Papers/Yield_Potential/Figures/RM')
png(file=paste("Crop_Yield_Soil_Carbon", Sys.Date(), '.png', sep="_"),width=800, height=600)
fixed_effects_plot
dev.off()


#normalized values
cols=c('Wheat'='red', 'Barley'='brown', 'Canola'='darkgoldenrod1', 'Peas'='blue')

norm_fixed_effects_plot <- ggplot() + xlim(c(0,1)) + ylim(0.35,0.65) + 
  geom_line(data=wheat_effects_n, aes(x=soc_norm, y=fit, colour='Wheat'), linetype=4) +
  geom_ribbon(data=wheat_effects_n, aes(x=soc_norm, ymin=lower, ymax=upper), alpha= 0.3, fill="red") +
  geom_line(data=barley_effects_n, aes(x=soc_norm, y=fit, colour='Barley'), linetype=2) +
  geom_ribbon(data= barley_effects_n, aes(x=soc_norm, ymin=lower, ymax=upper), alpha= 0.3, fill="brown") +
  geom_line(data=canola_effects_n, aes(x=soc_norm, y=fit, colour='Canola'), linetype=1) +
  geom_ribbon(data= canola_effects_n, aes(x=soc_norm, ymin=lower, ymax=upper), alpha= 0.3, fill="darkgoldenrod1") +
  geom_line(data=peas_effects_n, aes(x=soc_norm, y=fit, colour='Peas'), linetype=3) +
  geom_ribbon(data= peas_effects_n, aes(x=soc_norm, ymin=lower, ymax=upper), alpha= 0.3, fill="blue") +
  labs(x="Normalized Soil Organic Carbon", y=("Normalized Yield")) +
  scale_linetype_manual(values = c("Wheat" = 1, "Barley" = 2, "Canola" = 3, 'Peas' = 4)) +
  scale_colour_manual(name="Crop Types",values=cols,guide = guide_legend(override.aes = list(linetype = c(4,2,1,3))))+
  theme(text = element_text(size = 20))

norm_fixed_effects_plot

setwd('/home/preston/OneDrive/Papers/Yield_Potential/Figures/RM')
png(file=paste("Crop_Yield_Soil_Carbon_Norm", Sys.Date(), '.png', sep="_"),width=800, height=600)
norm_fixed_effects_plot
dev.off()


##############check differences by township###########
wheat_min_vals=aggregate(yield_wheat[,c(107,3)], by=list(yield_wheat$RMNO), FUN=min)
wheat_max_vals=aggregate(yield_wheat[,c(107,3)], by=list(yield_wheat$RMNO), FUN=max)
wheat_max_min_yield=max_vals$kg_ha_adj-min_vals$kg_ha_adj
wheat_max_min_soc=max_vals$soc-min_vals$soc

barley_min_vals=aggregate(yield_barley[,c(107,3)], by=list(yield_barley$RMNO), FUN=min)
barley_max_vals=aggregate(yield_barley[,c(107,3)], by=list(yield_barley$RMNO), FUN=max)
barley_max_min_yield=max_vals$kg_ha_adj-min_vals$kg_ha_adj
barley_max_min_soc=max_vals$soc-min_vals$soc

canola_min_vals=aggregate(yield_canola[,c(107,3)], by=list(yield_canola$RMNO), FUN=min)
canola_max_vals=aggregate(yield_canola[,c(107,3)], by=list(yield_canola$RMNO), FUN=max)
canola_max_min_yield=max_vals$kg_ha_adj-min_vals$kg_ha_adj
canola_max_min_soc=max_vals$soc-min_vals$soc

peas_min_vals=aggregate(yield_peas[,c(107,3)], by=list(yield_peas$RMNO), FUN=min)
peas_max_vals=aggregate(yield_peas[,c(107,3)], by=list(yield_peas$RMNO), FUN=max)
peas_max_min_yield=max_vals$kg_ha_adj-min_vals$kg_ha_adj
peas_max_min_soc=max_vals$soc-min_vals$soc


hist(wheat_max_min_yield)
hist(wheat_max_min_soc)
median(wheat_max_min_yield)
median(wheat_max_min_soc)

median(yield_wheat$kg_ha_adj)


#############Back Calculate Average Increase in SOC relates to increase in Yield################
#output_ab = output_01 * (b - a) + a,
wheat_effects_n_rescaled=wheat_effects_n
wheat_effects_n_rescaled[,1]=wheat_effects_n[,1]*(mean(wheat_max_vals$soc)-mean(wheat_min_vals$soc))+mean(wheat_min_vals$soc)
wheat_effects_n_rescaled[,-1]=wheat_effects_n[,-1]*(mean(wheat_max_vals$kg_ha_adj)-mean(wheat_min_vals$kg_ha_adj))+mean(wheat_min_vals$kg_ha_adj)
wheat_gain=(max(wheat_effects_n_rescaled$fit)-min(wheat_effects_n_rescaled$fit))/(max(wheat_effects_n_rescaled$soc_norm)-min(wheat_effects_n_rescaled$soc_norm))

barley_effects_n_rescaled=barley_effects_n
barley_effects_n_rescaled[,1]=barley_effects_n[,1]*(mean(barley_max_vals$soc)-mean(barley_min_vals$soc))+mean(barley_min_vals$soc)
barley_effects_n_rescaled[,-1]=barley_effects_n[,-1]*(mean(barley_max_vals$kg_ha_adj)-mean(barley_min_vals$kg_ha_adj))+mean(barley_min_vals$kg_ha_adj)
barley_gain=(max(barley_effects_n_rescaled$fit)-min(barley_effects_n_rescaled$fit))/(max(barley_effects_n_rescaled$soc_norm)-min(barley_effects_n_rescaled$soc_norm))

canola_effects_n_rescaled=canola_effects_n
canola_effects_n_rescaled[,1]=canola_effects_n[,1]*(mean(canola_max_vals$soc)-mean(canola_min_vals$soc))+mean(canola_min_vals$soc)
canola_effects_n_rescaled[,-1]=canola_effects_n[,-1]*(mean(canola_max_vals$kg_ha_adj)-mean(canola_min_vals$kg_ha_adj))+mean(canola_min_vals$kg_ha_adj)
canola_gain=(max(canola_effects_n_rescaled$fit)-min(canola_effects_n_rescaled$fit))/(max(canola_effects_n_rescaled$soc_norm)-min(canola_effects_n_rescaled$soc_norm))

peas_effects_n_rescaled=peas_effects_n
peas_effects_n_rescaled[,1]=peas_effects_n[,1]*(mean(peas_max_vals$soc)-mean(peas_min_vals$soc))+mean(peas_min_vals$soc)
peas_effects_n_rescaled[,-1]=peas_effects_n[,-1]*(mean(peas_max_vals$kg_ha_adj)-mean(peas_min_vals$kg_ha_adj))+mean(peas_min_vals$kg_ha_adj)
peas_gain=(max(peas_effects_n_rescaled$fit)-min(peas_effects_n_rescaled$fit))/(max(peas_effects_n_rescaled$soc_norm)-min(peas_effects_n_rescaled$soc_norm))


#rescaled values plot
cols=c('Wheat'='red', 'Barley'='brown', 'Canola'='darkgoldenrod1', 'Peas'='blue')

rescaled_fixed_effects_plot <- ggplot() + xlim(c(1.6,3.5)) + ylim(1500,3200) + 
  geom_line(data=wheat_effects_n_rescaled, aes(x=soc_norm, y=fit, colour='Wheat'), linetype=4) +
  geom_ribbon(data= wheat_effects_n_rescaled, aes(x=soc_norm, ymin=lower, ymax=upper), alpha= 0.3, fill="red") +
  geom_line(data=barley_effects_n_rescaled, aes(x=soc_norm, y=fit, colour='Barley'), linetype=2) +
  geom_ribbon(data= barley_effects_n_rescaled, aes(x=soc_norm, ymin=lower, ymax=upper), alpha= 0.3, fill="brown") +
  geom_line(data=canola_effects_n_rescaled, aes(x=soc_norm, y=fit, colour='Canola'), linetype=1) +
  geom_ribbon(data= canola_effects_n_rescaled, aes(x=soc_norm, ymin=lower, ymax=upper), alpha= 0.3, fill="darkgoldenrod1") +
  geom_line(data=peas_effects_n_rescaled, aes(x=soc_norm, y=fit, colour='Peas'), linetype=3) +
  geom_ribbon(data= peas_effects_n_rescaled, aes(x=soc_norm, ymin=lower, ymax=upper), alpha= 0.3, fill="blue") +
  labs(x="Soil Organic Carbon (%)", y=expression(paste("Yield (kg ha"^"-1",")"))) +
  scale_linetype_manual(values = c("Wheat" = 1, "Barley" = 2, "Canola" = 3, 'Peas' = 4)) +
  scale_colour_manual(name="Crop Types",values=cols,guide = guide_legend(override.aes = list(linetype = c(4,2,1,3))))+
  theme(text = element_text(size = 20))

rescaled_fixed_effects_plot

setwd('/home/preston/OneDrive/Papers/Yield_Potential/Figures/RM')
png(file=paste("Crop_Yield_Soil_Carbon_rescaled", Sys.Date(), '.png', sep="_"),width=800, height=600)
rescaled_fixed_effects_plot
dev.off()

#export average gain per percent SOC
gains=c(wheat_gain, barley_gain, canola_gain, peas_gain)
names(gains)=c("Wheat", 'Barley', 'Canola', 'Peas')
write.csv(gains, paste("Average_Gains_Crop_Type_", Sys.Date(), '.csv', sep=""))

#prices
#https://tradingeconomics.com/commodity/canola
#HRSW = $638.87/tonne
#Barley = $371.27 / tonne
#Canola = $892.80 / tonne
#Peas  = $457 / tonne

wheat_gain*637.87/1000*65
barley_gain*371.27/1000*65
canola_gain*892.80/1000*65
peas_gain*457/1000*65


wheat_gain*637.87/1000*65/((median(yield_wheat$kg_ha_adj))*637.87/1000*65)
barley_gain*371.27/1000*65/((median(yield_barley$kg_ha_adj))*371.27/1000*65)
canola_gain*892.80/1000*65/((median(yield_canola$kg_ha_adj))*892.80/1000*65)
peas_gain*457/1000*65/((median(yield_peas$kg_ha_adj))*457/1000*65)



##################Performance by RM##############
wheat_data_counts=table(yield_wheat$RMNO)
barley_data_counts=table(yield_barley$RMNO)
canola_data_counts=table(yield_canola$RMNO)
peas_data_counts=table(yield_peas$RMNO)

quantile(wheat_data_counts)
quantile(barley_data_counts)
quantile(canola_data_counts)
quantile(peas_data_counts)

#drop those below 82 -> 1st percentile
wheat_twshp=names(wheat_data_counts[wheat_data_counts>500])
barley_twshp=names(barley_data_counts[barley_data_counts>500])
canola_twshp=names(canola_data_counts[canola_data_counts>500])
peas_twshp=names(peas_data_counts[peas_data_counts>500])


#rm shapefile
rm_shp=shapefile('/media/preston/My Book/Saskatchewan/RM_Shapefile/RM.shp')

i='126'
#Wheat SOC by RM
wheat_model_gain=foreach(i=unique(wheat_twshp), .combine=rbind) %dopar%
  {
    try({
    dat_sub=yield_wheat_norm[yield_wheat_norm$RMNO==i,]
    model_sub=lm(kg_ha_adj~soc ,data=dat_sub)
    model_results=c(i, model_sub$coefficients[2], median(dat_sub$soc), summary(model_sub)$coefficients[2,4])
    names(model_results)=c('RMNO', "yield_slope", 'median_soc', 'p-value')
    model_results
    })
    }

wheat_model_gain=apply(wheat_model_gain, 2, FUN=as.numeric)
wheat_model_gain=data.frame(wheat_model_gain)

barley_model_gain=foreach(i=unique(barley_twshp), .combine=rbind) %dopar%
  {
    try({
      dat_sub=yield_barley_norm[yield_barley_norm$RMNO==i,]
      model_sub=lm(kg_ha_adj~soc ,data=dat_sub)
      model_results=c(i, model_sub$coefficients[2], median(dat_sub$soc), summary(model_sub)$coefficients[2,4])
      names(model_results)=c('RMNO', "yield_slope", 'median_soc', 'p-value')
      model_results
    })
  }

barley_model_gain=apply(barley_model_gain, 2, FUN=as.numeric)
barley_model_gain=data.frame(barley_model_gain)


canola_model_gain=foreach(i=unique(canola_twshp), .combine=rbind) %dopar%
  {
    try({
      dat_sub=yield_canola_norm[yield_canola_norm$RMNO==i,]
      model_sub=lm(kg_ha_adj~soc ,data=dat_sub)
      model_results=c(i, model_sub$coefficients[2], median(dat_sub$soc), summary(model_sub)$coefficients[2,4])
      names(model_results)=c('RMNO', "yield_slope", 'median_soc', 'p-value')
      model_results
    })
  }

canola_model_gain=apply(canola_model_gain, 2, FUN=as.numeric)
canola_model_gain=data.frame(canola_model_gain)

peas_model_gain=foreach(i=unique(peas_twshp), .combine=rbind) %dopar%
  {
    try({
      dat_sub=yield_peas_norm[yield_peas_norm$RMNO==i,]
      model_sub=lm(kg_ha_adj~soc ,data=dat_sub)
      model_results=c(i, model_sub$coefficients[2], median(dat_sub$soc), summary(model_sub)$coefficients[2,4])
      names(model_results)=c('RMNO', "yield_slope", 'median_soc', 'p-value')
      model_results
    })
  }

peas_model_gain=apply(peas_model_gain, 2, FUN=as.numeric)
peas_model_gain=data.frame(peas_model_gain)


#ggplot
cols=c('Wheat'='red', 'Barley'='brown', 'Canola'='darkgoldenrod1', 'Peas'='blue')
soc_by_rm <- ggplot() +xlim(c(1.5, 3.5)) + 
  geom_smooth(data=wheat_model_gain, aes(x=median_soc, y=yield_slope, colour='Wheat', fill='Wheat'), linetype=4, span=2) +
  geom_smooth(data=barley_model_gain, aes(x=median_soc, y=yield_slope,  colour='Barley', fill='Barley'), linetype=2, span=2) +
  geom_smooth(data=canola_model_gain, aes(x=median_soc, y=yield_slope, colour='Canola', fill='Canola'), linetype=1, span=2) +
  geom_smooth(data=peas_model_gain, aes(x=median_soc, y=yield_slope,  colour='Peas', fill='Peas'), linetype=3, span=2) +
  labs(x="Soil Organic Carbon (%)", y=expression(paste("Yield Response (kg ha"^"-1",")"))) +
  guides(fill=FALSE)+
  scale_linetype_manual(values = c("Wheat" = 1, "Barley" = 2, "Canola" = 3, 'Peas' = 4)) +
  scale_colour_manual(name="Crop Types",values=cols, guide = guide_legend(override.aes = list(linetype = c(4,2,1,3))))+
  theme(text = element_text(size = 20))

soc_by_rm

setwd('/home/preston/OneDrive/Papers/Yield_Potential/Figures/RM')
png(file=paste("Crop_Yield_Response_By_RM", Sys.Date(), '.png', sep="_"),width=800, height=600)
soc_by_rm
dev.off()


# RMNO variables
#export where negative
wheat_model_gain$slope_type='positive'
wheat_model_gain$slope_type[wheat_model_gain$yield_slope<0] <- 'negative'

yield_testing=merge(yield_wheat_norm, wheat_model_gain, by='RMNO')

boxplot(yield_testing$precip~yield_testing$slope_type, ylim=c(400, 700))


##############temporal trend plot#########
plot(yields_wheat_avg$x~yields_wheat_avg$Group.1)
yields_wheat_avg
yields_barley_avg
yields_canola_avg
yields_peas_avg


#ggplot
cols=c('Wheat'='red', 'Barley'='brown', 'Canola'='darkgoldenrod1', 'Peas'='blue')
yields_trend <- ggplot() +
  geom_smooth(data=yields_wheat_avg, aes(x=Year, y=Yield, colour='Wheat', fill='Wheat'), linetype=4, se=FALSE, method=lm, formula=y~x) +
  geom_smooth(data=yields_barley_avg, aes(x=Year, y=Yield,  colour='Barley', fill='Barley'), linetype=2, , se=FALSE, method=lm, formula=y~x) +
  geom_smooth(data=yields_canola_avg, aes(x=Year, y=Yield, colour='Canola', fill='Canola'), linetype=1, , se=FALSE, method=lm, formula=y~x) +
  geom_smooth(data=yields_peas_avg, aes(x=Year, y=Yield, colour='Peas', fill='Peas'), linetype=1, , se=FALSE, method=lm, formula=y~x) +
  geom_point(data=yields_wheat_avg, aes(x=Year, y=Yield, colour='Wheat', fill='Wheat')) +
  geom_point(data=yields_barley_avg, aes(x=Year, y=Yield,  colour='Barley', fill='Barley')) +
  geom_point(data=yields_canola_avg, aes(x=Year, y=Yield, colour='Canola', fill='Canola')) +
  geom_point(data=yields_peas_avg, aes(x=Year, y=Yield, colour='Peas', fill='Peas')) +
  labs(x="Year", y=expression(paste("Yield (kg ha"^"-1",")"))) +
  guides(fill='none')+
  scale_linetype_manual(values = c("Wheat" = 1, "Barley" = 2, "Canola" = 3, 'Peas' = 4)) +
  scale_colour_manual(name="Crop Types",values=cols, guide = guide_legend(override.aes = list(linetype = c(4,2,1,3))))+
  theme(text = element_text(size = 20))

yields_trend

setwd('/home/preston/OneDrive/Papers/Yield_Potential/Figures/RM')
png(file=paste("Crop_Yield_Response_By_RM", Sys.Date(), '.png', sep="_"),width=800, height=600)
soc_by_rm
dev.off()



##########bayesian code##################
library(brms)

pr = prior(normal(0, 1), class = 'b')
model_barley=brm(kg_ha_adj_norm~soc_norm + (soc_norm||RMNO), data=yield_wheat_norm, prior=pr, cores=8)
summary(model_barley)
#hypothesis test
variables(model_barley)

hypothesis(model_barley, 'Intercept < Intercept + barley', class='b', alpha=0.1)

