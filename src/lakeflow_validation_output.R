# Produce validation stats and figures for rmd file to pull. 
# By: Ryan Riggs
# Date: 8/6/2024
################################################################################
inPath='C:/Users/rriggs/OneDrive - DOI/Research/LakeFlow_local/'
################################################################################
# Load Packages
################################################################################
library(foreign)
library(lubridate)
library(rstan)
library(data.table)
library(dplyr)
library(sf)
library(raster)
library(rstudioapi)
library(jsonlite)
library(httr)
library(BBmisc)
library(reticulate)
library(ggplot2)
library(USA.state.boundaries)
library(wesanderson)
library(ggpubr)
'%!in%' <- function(x,y)!('%in%'(x,y))

output_files = list.files(paste0(inPath, '/out/lf_results_dynamic/'), full.names=TRUE)
lf_outputs = rbindlist(lapply(output_files, fread),fill=TRUE)
lf_outputs$Date = as.Date(lf_outputs$Date)
lf_outputs$reach_id = as.character(lf_outputs$reach_id)
lf_outputs$lake_id=as.character(lf_outputs$lake_id)
lf_outputs = lf_outputs[exp(lf_outputs$width)>=50,]
flow = lf_outputs
lake_reach = unique(lf_outputs[,c('lake_id', 'reach_id')])
################################################################################
# Validate lakeflow observations with USGS data.  
################################################################################
usgs_gages = fread('C:\\Users\\rriggs\\OneDrive - DOI\\Research\\SWOT\\gage_sword_snapped.csv', colClasses = as.character('site_no'))
usgs_gages = usgs_gages[usgs_gages$reach_id%in%lf_outputs$reach_id,]

pull_gage = function(f){
  data = try(RivRetrieve::usa(site=f, variable='discharge'))
  if(is.error(data)){return(NA)}
  data$gage=f
  data=data[data$Date>=as.Date('2022-12-31'),]
  return(data)
}

gage_df_list = lapply(as.character(usgs_gages$site_no), pull_gage)
gage_df = rbindlist(gage_df_list[!is.na(gage_df_list)])
gage_df$reach_id = usgs_gages$reach_id[match(gage_df$gage, usgs_gages$site_no)]
gage_df$reach_id=as.character(gage_df$reach_id)
gage_df = gage_df[,meanQ:=mean(Q, na.rm=TRUE),by=gage]
gage_df = gage_df[meanQ>=10,]

lf_outputs$Date = as.Date(lf_outputs$date)
# Calculate the mean outputs for each reach - sometimes a reach can be found in multiple lakes. 
lf_outputs = lf_outputs[,-c('lake_id', 'type')][,lapply(.SD, mean), by=list(Date, reach_id)]
validation = left_join(gage_df, lf_outputs,by=c('Date', 'reach_id'))
full_gage_ts = validation
################################################################################
# Determine which locations are reservoirs and which are natural lakes. 
################################################################################
# Read in harmonized sword pld
sword_lake_db = 'C:/Users/rriggs/OneDrive - DOI/Research/LakeFlow_local/in/Intersected_SWORD_PLD_dataset.gdb'
layer_names= sf::st_layers(sword_lake_db)
sword_lake_db_files = layer_names$name

open_gdb = function(f){
  gdb = sf::st_read(sword_lake_db, layer=f, quiet=TRUE)
  gdb = gdb[gdb$U_reach_n>=0&gdb$D_reach_n>=0&as.character(gdb$lake_id)%in%lake_reach$lake_id,]
  if(nrow(gdb)==0){return(NA)}
  return(gdb)
}

gdb_list = lapply(sword_lake_db_files, FUN=open_gdb)
gdb = data.table::rbindlist(gdb_list[!is.na(gdb_list)], use.names=TRUE,fill=TRUE)
updated_pld = gdb
updated_pld$lake_id =  as.character(updated_pld$lake_id)
updated_pld = st_sf(updated_pld)
updated_pld = st_make_valid(updated_pld)
geodar = st_read('C:/Users/rriggs/OneDrive - DOI/Research/MISC/GeoDAR/GeoDAR_v10_v11/GeoDAR_v11_reservoirs.shp')
geodar = st_transform(geodar, st_crs(updated_pld))
geodar = st_make_valid(geodar)
jn = st_join(updated_pld, geodar, st_intersects)
jn$type = 'NA'
jn$type = ifelse(!is.na(jn$OBJECTID),'Reservoir','Natural')
lake_reach$lake_type = NA
lake_reach$lake_type = jn$type[match(lake_reach$lake_id, jn$lake_id)]
reach_type = unique(lake_reach[,c('reach_id','lake_type')])
validation$lake_type = reach_type$lake_type[match(validation$reach_id, reach_type$reach_id)]
################################################################################
# Error metrics. 
################################################################################
source('C:\\Users\\rriggs\\OneDrive - DOI\\Research\\Error_stats\\src\\error_functions.R')
validation = validation[!is.na(Q)&!is.na(q_lakeflow),]
validation[,OBS:=.N,by=reach_id]
validation[,meanQ:=mean(Q),by=reach_id]
validation = validation[OBS>=5&!is.na(meanQ),]
error_stats=data.table(validation)[,validate(q_lakeflow, Q),by=reach_id]
error_stats_model=data.table(validation)[,validate(q_model, Q),by=reach_id]

error_stats$type='LakeFlow'
error_stats_model$type='Geoglows'

error_combined = bind_rows(error_stats, error_stats_model)
error_combined$'|nBias|' = abs(error_combined$rBias)
error_tall = melt(error_combined, id.vars=c('type', 'reach_id'))
error_tall$meanQ =validation$meanQ[match(error_tall$reach_id, validation$reach_id)]
cl = wes_palette('AsteroidCity3')[3:4]
cl = rev(cl)

# cdf plots. 
a=ggplot(error_tall[variable%in%c('NRMSE'),])+
  stat_ecdf(aes(x=value,color=type),lwd=1.5)+
  facet_wrap(~variable, scales='free')+
  coord_cartesian(xlim=c(0,300))+
  theme_minimal()+
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=16),
        legend.position = 'top',
        legend.title = element_blank(),
        strip.text = element_text(size=18),
        legend.text = element_text(size=18))+
  xlab('Percent(%)')+
  ylab('CDF')+
  scale_color_manual(values=cl)
a
b=ggplot(error_tall[variable%in%c('KGE'),])+
  stat_ecdf(aes(x=value,color=type),lwd=1.5)+
  facet_wrap(~variable)+
  coord_cartesian(xlim=c(-1,1))+
  theme_minimal()+
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=16),
        legend.position = 'top',
        legend.title = element_blank(),
        strip.text = element_text(size=18),
        legend.text = element_text(size=18))+
  ylab('CDF')+
  scale_color_manual(values=cl)
b
c=ggplot(error_tall[variable%in%c('|nBias|')])+
  stat_ecdf(aes(x=abs(value),color=type),lwd=1.5)+
  facet_wrap(~variable)+
  coord_cartesian(xlim=c(0,100))+
  theme_minimal()+
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=16),
        legend.position = 'top',
        strip.text = element_text(size=18),
        legend.title = element_blank(),
        legend.text = element_text(size=18))+
  xlab('Percent(%)')+
  ylab('CDF')+
  scale_color_manual(values=cl)
c
pt = ggarrange(a,c,b, labels='auto',nrow=1,common.legend = TRUE)
################################################################################
# Write out validation data. 
error_tall$lake_type = lake_reach$lake_type[match(error_tall$reach_id, lake_reach$reach_id)]
error_tall$flow = flow$type[match(error_tall$reach_id, flow$reach_id)]
error_tall$gage = gage_df$gage[match(error_tall$reach_id, gage_df$reach_id)]
error_tall$obs = validation$OBS[match(error_tall$reach_id, validation$reach_id)]
error_tall$lake_id = lake_reach$lake_id[match(error_tall$reach_id, lake_reach$reach_id)]
error_tall$upstream = updated_pld$U_reach_n[match(error_tall$lake_id, updated_pld$lake_id)]
error_tall$downstream = updated_pld$D_reach_n[match(error_tall$lake_id, updated_pld$lake_id)]
error_tall$reaches = error_tall$upstream+error_tall$downstream

fwrite(error_tall, paste0(inPath, 'out/stats/', 'validation.csv'))
################################################################################
# Compare daily vs static prior Q results. 
################################################################################
output_files = list.files(paste0(inPath, '/out/lf_results_dynamic/'), full.names=TRUE)
lf_outputs = rbindlist(lapply(output_files, fread),fill=TRUE)
lf_outputs$Date = as.Date(lf_outputs$Date)
lf_outputs$reach_id = as.character(lf_outputs$reach_id)
lf_outputs$lake_id=as.character(lf_outputs$lake_id)
lf_outputs = lf_outputs[exp(lf_outputs$width)>=50,]
geoglows = lf_outputs
geoglows$model = 'dynamic'

output_files = list.files(paste0(inPath, '/out/lf_results_static/'), full.names=TRUE)
lf_outputs = rbindlist(lapply(output_files, fread),fill=TRUE)
lf_outputs$Date = as.Date(lf_outputs$Date)
lf_outputs$reach_id = as.character(lf_outputs$reach_id)
lf_outputs$lake_id=as.character(lf_outputs$lake_id)
lf_outputs = lf_outputs[exp(lf_outputs$width)>=50,]
grades = lf_outputs
grades$model = 'static'
comb = bind_rows(geoglows, grades)
comb$Date = as.Date(comb$date)
################################################################################
# Validate lakeflow observations with USGS data.  
################################################################################
usgs_gages = fread('C:\\Users\\rriggs\\OneDrive - DOI\\Research\\SWOT\\gage_sword_snapped.csv', colClasses = as.character('site_no'))
usgs_gages = usgs_gages[usgs_gages$reach_id%in%comb$reach_id,]

pull_gage = function(f){
  data = try(RivRetrieve::usa(site=f, variable='discharge'))
  if(is.error(data)){return(NA)}
  data$gage=f
  data=data[data$Date>=as.Date('2022-12-31'),]
  return(data)
}

gage_df_list = lapply(as.character(usgs_gages$site_no), pull_gage)
gage_df = rbindlist(gage_df_list[!is.na(gage_df_list)])
gage_df$reach_id = usgs_gages$reach_id[match(gage_df$gage, usgs_gages$site_no)]
gage_df$reach_id=as.character(gage_df$reach_id)
gage_df = gage_df[,meanQ:=mean(Q, na.rm=TRUE),by=gage]
gage_df = gage_df[meanQ>=10,]

# Calculate the mean outputs for each reach - sometimes a reach can be found in multiple lakes. 
comb = comb[,-c('lake_id', 'type')][,lapply(.SD, mean), by=list(Date, reach_id, model)]
validation = left_join(gage_df, comb,by=c('Date', 'reach_id'))
################################################################################
# Error metrics. 
################################################################################
source('C:\\Users\\rriggs\\OneDrive - DOI\\Research\\Error_stats\\src\\error_functions.R')
validation = validation[!is.na(Q)&!is.na(q_lakeflow),]
validation[,OBS:=.N,by=list(reach_id, model)]
validation[,meanQ:=mean(Q),by=list(reach_id, model)]
validation = validation[OBS>=5&!is.na(meanQ),]
error_stats=data.table(validation)[,validate(q_lakeflow, Q),by=list(reach_id,model)]
error_stats_model=data.table(validation)[,validate(q_model, Q),by=list(reach_id,model)]
################################################################################
# Summary error metrics plots. 
################################################################################
error_stats$type='LakeFlow'
error_stats_model$type='Model'

error_combined = bind_rows(error_stats, error_stats_model)
error_combined$'|nBias|' = abs(error_combined$rBias)
error_tall = melt(error_combined, id.vars=c('type', 'reach_id','model'))
# Write out comparison data. 
fwrite(error_tall, paste0(inPath, 'out/stats/', 'validation_comparison.csv'))
