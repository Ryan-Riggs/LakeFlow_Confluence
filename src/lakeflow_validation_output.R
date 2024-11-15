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
library(tidyhydat)
'%!in%' <- function(x,y)!('%in%'(x,y))

output_files = list.files(paste0(inPath, '/out/lf_results_na_3/'), full.names=TRUE)
lf_outputs = rbindlist(lapply(output_files, fread),fill=TRUE)
lf_outputs$Date = as.Date(lf_outputs$Date)
lf_outputs$reach_id = as.character(lf_outputs$reach_id)
lf_outputs$lake_id=as.character(lf_outputs$lake_id)
lf_outputs = lf_outputs[exp(lf_outputs$width)>=50&q_lakeflow>0&q_lakeflow<1e6,]
flow = lf_outputs
lake_reach = unique(lf_outputs[,c('lake_id', 'reach_id')])
################################################################################
# Validate lakeflow observations with USGS data.  
################################################################################
usgs_gages = fread('C:\\Users\\rriggs\\OneDrive - DOI\\Research\\SWOT\\gage_usgs_can_sword_snapped.csv', colClasses = as.character('site_no'))
usgs_gages = usgs_gages[usgs_gages$reach_id%in%lf_outputs$reach_id,]

pull_gage_usgs = function(f){
  data = try(RivRetrieve::usa(site=f, variable='discharge'))
  if(is.error(data)){return(NA)}
  data$gage=f
  data=data[data$Date>=as.Date('2022-12-31'),]
  return(data)
}

pull_gage_hydat = function(f){
  data = try(RivRetrieve::canada(site=f, variable='discharge'))
  if(is.error(data)){return(NA)}
  data$gage=f
  data=data[data$Date>=as.Date('2022-12-31'),]
  return(data)
}

pull_gage_hydat_nrt = function(f){
  can = realtime_ws(station_number = f, start_date = as.Date('2024-01-01'))
  can$Date = as.Date(can$Date)
  can_agg = data.table(can)[Unit=='m3/s',list(Q=mean(Value)),by=list(Date, STATION_NUMBER)]
  can_agg$gage = can_agg$STATION_NUMBER
  return(can_agg[,c('Date', 'Q', 'gage')])
}

pull_gage_quebec=function(f){
  website=paste0('https://www.cehq.gouv.qc.ca/depot/historique_donnees/fichier/',f,'_Q.txt')
  outpath=tempfile()
  downloading = try(download.file(website, outpath,quiet=TRUE),silent=TRUE)
  if(is.error(downloading)){return(NA)}
  data=fread(outpath, fill=TRUE)
  removeInfo=grep('Date', data$V2)
  data=data[(removeInfo+1):nrow(data),2:3]
  colnames(data)=c('Date','Q')
  data = data[!is.na(Date)&!is.na(Q),]
  if(nrow(data)<3){return(NA)}
  data$Date=as.Date(data$Date)
  data$Q=as.numeric(data$Q)
  data$gage=as.character(f)
  return(data)
}

gage_df_list_usgs = lapply(as.character(usgs_gages$site_no[usgs_gages$agency_cd=='USGS']), pull_gage_usgs)
gage_df_list_hydat = lapply(as.character(usgs_gages$site_no[usgs_gages$agency_cd=='Hydat']), pull_gage_hydat)
gage_df_list_hydat_nrt = list(pull_gage_hydat_nrt(usgs_gages$site_no[usgs_gages$agency_cd=='Hydat']))
gage_df_list_quebec = lapply(as.character(usgs_gages$site_no[usgs_gages$agency_cd=='Quebec']), pull_gage_quebec)

gage_df = rbindlist(c(gage_df_list_usgs[!is.na(gage_df_list_usgs)], gage_df_list_hydat[!is.na(gage_df_list_hydat)], gage_df_list_quebec[!is.na(gage_df_list_quebec)], gage_df_list_hydat_nrt[!is.na(gage_df_list_hydat_nrt)]))
gage_df$reach_id = usgs_gages$reach_id[match(gage_df$gage, usgs_gages$site_no)]
gage_df$reach_id=as.character(gage_df$reach_id)
gage_df = gage_df[,meanQ:=mean(Q, na.rm=TRUE),by=gage]
#gage_df = gage_df[meanQ>=10,]

lf_outputs$Date = as.Date(lf_outputs$date)
# Calculate the mean outputs for each reach - sometimes a reach can be found in multiple lakes. 
lf_outputs = lf_outputs[,-c('lake_id', 'type', 'prior_fit')][,lapply(.SD, mean), by=list(Date, reach_id)]
validation = left_join(gage_df, lf_outputs,by=c('Date', 'reach_id'))
full_gage_ts = validation
################################################################################
# Determine which locations are reservoirs and which are natural lakes. 
################################################################################
# Read in harmonized sword pld
sword_lake_db = 'C:/Users/rriggs/OneDrive - DOI/Research/LakeFlow_local/in/Intersected_SWORD_PLD_dataset_wo_ghost_rch.gdb'
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
validation = validation[!is.na(Q)&!is.na(q_lakeflow)&Q>1,]
validation = distinct(validation)
validation[,OBS:=.N,by=reach_id]
#validation[,meanQ:=mean(Q),by=reach_id]
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
a=ggplot(error_tall[variable%in%c('NSE'),])+
  stat_ecdf(aes(x=value,color=type),lwd=1.5)+
  facet_wrap(~variable, scales='free')+
  coord_cartesian(xlim=c(-5,1))+
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
error_tall$agency = usgs_gages$agency_cd[match(error_tall$gage, usgs_gages$site_no)]
error_tall = error_tall[meanQ>=10,] # Helps remove problem gages that are located nearby but not on mainstem. 

fwrite(error_tall, paste0(inPath, 'out/stats/', 'validation.csv'))
################################################################################
# Produce validation figures. 
################################################################################
library(ggpubr)
library(wesanderson)
library(see)
library(gghalves)
library(dplyr)
library(ggnewscale)
library(ggh4x)
comparison = fread('C:/Users/rriggs/OneDrive - DOI/Research/LakeFlow_local/out/stats/validation.csv')

comparison_2 = comparison#[type=='LakeFlow']
all = comparison_2
all$lake_type='Total'
comparison_2 = bind_rows(all, comparison_2)
comparison_2$class = comparison_2$lake_type
# comparison_2$value = ifelse(comparison_2$variable=='KGE'|comparison_2$variable=='Rvalue', comparison_2$value*100, comparison_2$value)
iqr = comparison_2[,list(first = quantile(value, .25), second=quantile(value,.75)),by=list(variable, class)]
iqr = iqr[,range:=second-first,by=list(variable, class)]
iqr = iqr[,list(first=first-(1.5*(range)),second=second+(1.5*(range))), by=list(variable, class)]
comparison_2 = left_join(comparison_2, iqr, by=c('variable', 'class'))
# comparison_2$variable = ifelse(comparison_2$variable=='KGE', 'KGE x 100', comparison_2$variable)
comparison_2$variable = ifelse(comparison_2$variable=='Rvalue', 'R', comparison_2$variable)
comparison_2 = comparison_2[variable%in%c('NRMSE','|nBias|','KGE', 'R'),]
comparison_2$variable = factor(comparison_2$variable, levels=c('R','KGE','|nBias|', 'NRMSE'))

static_meds = comparison[variable%in%c('NRMSE','|nBias|','KGE')&type=='LakeFlow',median(value),by=variable]

comparison_2$type = factor(comparison_2$type, levels=c('LakeFlow', 'Geoglows'))

nat_res=ggplot(comparison_2)+
  stat_ecdf(aes(color=lake_type,x=value, lty=type))+
  #coord_cartesian(xlim=c(-100, 100))+
  theme_minimal()+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10),
        strip.placement = "outside")+
  ylab('CDF')+
  xlab('')+
  scale_color_manual(values=c(wes_palette('Zissou1')[c(2,4)], wes_palette('AsteroidCity1')[4]))+
  facet_wrap(~variable, nrow=1, scales='free_x',strip.position="bottom")

cdf_perf=nat_res+facetted_pos_scales(x=list(
  scale_x_continuous(limits=c(-1,1)),
  scale_x_continuous(limits=c(-1,1)),
  scale_x_continuous(limits=c(0,100)),
  scale_x_continuous(limits=c(0,100))
), y=NULL)


high_perf = error_tall[variable=='Rvalue'&value>0.5&lake_type=='Reservoir'&flow=='outflow']
high_perf = error_tall[variable=='Rvalue'&value>0.5&lake_type=='Natural'&flow=='outflow']

nat = '08LF033'#'05247500'#'02YQ001'
res = '05331580'

plot_df=full_gage_ts[gage%in%c(nat,res)]#[gage%in%high_perf$gage,]
plot_df$type = NA
plot_df$type = ifelse(plot_df$gage==nat, 'Example Natural Lake Outflow', 'Example Reservoir Outflow')


hy_plt=ggplot(plot_df[Date>as.Date('2023-01-01'),])+
  geom_line(aes(x=Date, y=Q), col='grey50', lwd=0.25)+
  #geom_line(aes(x=Date,y=q_lakeflow,col='red'))+
  #geom_point(aes(x=Date,y=exp(bayes_q)), col='black')+
  #geom_hline(aes(yintercept=q_model,col=reach_id))+
  #geom_point(aes(x=Date,y=q_lakeflow),col='#006F41')+
  geom_point(aes(x=Date,y=q_lakeflow), col='red', size=0.75)+
  geom_point(aes(x=Date,y=q_model),col='lightpink', size=0.75)+
  #geom_line(aes(x=Date, y=q_model), col='grey50')+
  facet_wrap(~type,ncol=2, scales = 'free_x')+
  #scale_y_log10(n.breaks=4)+
  theme_minimal()+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10))+
  xlab('')+
  ylab('Discharge (cms)')
hy_plt

val_plot = ggarrange(hy_plt, cdf_perf, labels='auto', nrow=2)
val_plot


ggsave('C:\\Users\\rriggs\\OneDrive - DOI\\Research\\LakeFlow_local\\out\\figures\\manuscript\\validation_figure.pdf', val_plot, width=6.5, height=6.5, dpi=1000, units='in')






















