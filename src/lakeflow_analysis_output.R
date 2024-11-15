#**** Clean up code

# Produce output datasets for rmd file. 
# By: Ryan Riggs
# Date: 8/7/2024
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
library(ggpmisc)
'%!in%' <- function(x,y)!('%in%'(x,y))

output_files = list.files(paste0(inPath, '/out/lf_results_na_3/'), full.names=TRUE)
lf_outputs = rbindlist(lapply(output_files, fread),fill=TRUE)
lf_outputs$Date = as.Date(lf_outputs$Date)
lf_outputs$reach_id = as.character(lf_outputs$reach_id)
lf_outputs$lake_id=as.character(lf_outputs$lake_id)
lf_outputs = lf_outputs[exp(lf_outputs$width)>=50&q_lakeflow>0&q_lakeflow<1e6]
agg = lf_outputs[,list(q_lakeflow = sum(q_lakeflow, na.rm=TRUE)), by=list(date, type,lake_id)]
agg = agg[q_lakeflow>0,]
inflow = agg[type=='inflow',]
outflow = agg[type=='outflow',]
comb = merge(inflow, outflow, by=c('date', 'lake_id'))
comb = comb[,OBS:=.N,by=lake_id]
comb = comb[OBS>=5,]
comb$diff_q = comb$q_lakeflow.x - comb$q_lakeflow.y
agg = comb[,list(mean_dq = mean(diff_q, na.rm=TRUE), inflow = mean(q_lakeflow.x, na.rm=TRUE), outflow = mean(q_lakeflow.y, na.rm=TRUE), inflow_sd = sd(q_lakeflow.x, na.rm=TRUE), outflow_sd = sd(q_lakeflow.y, na.rm = TRUE)), by=lake_id]

# Read in harmonized sword pld
sword_lake_db = 'C:/Users/rriggs/OneDrive - DOI/Research/LakeFlow_local/in/Intersected_SWORD_PLD_dataset_wo_ghost_rch.gdb'
layer_names= sf::st_layers(sword_lake_db)
sword_lake_db_files = layer_names$name

open_gdb = function(f){
  gdb = sf::st_read(sword_lake_db, layer=f, quiet=TRUE)
  gdb = gdb[gdb$U_reach_n>=0&gdb$D_reach_n>=0&as.character(gdb$lake_id)%in%lf_outputs$lake_id,]
  if(nrow(gdb)==0){return(NA)}
  return(gdb)
}

gdb_list = lapply(sword_lake_db_files, FUN=open_gdb)
gdb = data.table::rbindlist(gdb_list[!is.na(gdb_list)], use.names=TRUE,fill=TRUE)
updated_pld = gdb
updated_pld$lake_id =  as.character(updated_pld$lake_id)
# Limit to North America to speed up code. 
updated_pld$continent = substring(updated_pld$lake_id, 1, 1)
#updated_pld = updated_pld[updated_pld$continent=='7',]
updated_pld = left_join(updated_pld, agg, by='lake_id')
updated_pld = st_sf(updated_pld)
state_boundaries_wgs84 = state_boundaries_wgs84[state_boundaries_wgs84$NAME%!in%c('Puerto Rico', 'Alaska', 'Hawaii','U.S. Virgin Islands'),]
updated_pld = st_make_valid(updated_pld)
updated_pld$ratio = updated_pld$inflow/updated_pld$outflow
updated_pld$ratioCap = ifelse(updated_pld$ratio>2, 2, updated_pld$ratio)
geodar = st_read('C:/Users/rriggs/OneDrive - DOI/Research/MISC/GeoDAR/GeoDAR_v10_v11/GeoDAR_v11_reservoirs.shp')
geodar = st_transform(geodar, st_crs(updated_pld))
geodar = st_make_valid(geodar)
jn = st_join(updated_pld, geodar, st_intersects)
jn$type = 'NA'
jn$type = ifelse(!is.na(jn$OBJECTID),'Reservoir','Natural')
updated_pld$type = jn$type[match(updated_pld$lake_id, jn$lake_id)]
#Add in countries/states for rate
world = rnaturalearth::ne_countries(scale = "medium", returnclass="sf")
world = st_join(updated_pld, world, st_within)
states = tigris::states()%>%st_transform(crs(updated_pld))
states = st_join(updated_pld, states, st_within)
#output data for rmd file. 
output = lf_outputs
output$lake_type = updated_pld$type[match(output$lake_id, updated_pld$lake_id)]
output$Country = world$adm0_a3[match(output$lake_id,world$lake_id)]
output$States = states$STUSPS[match(output$lake_id, states$lake_id)]
output$lake_type = updated_pld$type[match(output$lake_id, updated_pld$lake_id)]
fwrite(output, paste0(inPath, 'out/stats/', 'lf_outputs.csv'))
################################################################################
# Analysis figures. 
################################################################################
#Lakeflow applied across region
library(data.table)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(ggpmisc)
lf_outputs = fread('C:/Users/rriggs/OneDrive - DOI/Research/LakeFlow_local/out/stats/lf_outputs.csv')
agg = lf_outputs[,list(q_lakeflow = sum(q_lakeflow, na.rm=TRUE)), by=list(date, type,lake_id, lake_type)]
agg = agg[q_lakeflow>0,]
inflow = agg[type=='inflow',]
outflow = agg[type=='outflow',]
comb = merge(inflow, outflow, by=c('date', 'lake_id', 'lake_type'))
comb = comb[,OBS:=.N,by=lake_id]
comb = comb[OBS>=5,]
comb$diff_q = comb$q_lakeflow.x - comb$q_lakeflow.y

#theil sen
sen <- function(..., weights = NULL) {
  mblm::mblm(...)
}

sen_vals = comb[,data.table(mblm::mblm(q_lakeflow.y~q_lakeflow.x, data=.SD)$coefficients), by=lake_type]
sen_vals=data.table(slope=c(sen_vals$V1[2],sen_vals$V1[4]), int = c(sen_vals$V1[1],sen_vals$V1[3]), lake_type=c('Natural', 'Reservoir'))
vals = seq(min(lf_outputs$q_lakeflow, na.rm=TRUE), max(lf_outputs$q_lakeflow, na.rm=TRUE), 10)
sen_df = sen_vals[,list(x=vals, y=vals*slope+int), by=lake_type]

ts_plot=ggplot(comb[!is.na(q_lakeflow.x)&!is.na(q_lakeflow.y),],aes(x=q_lakeflow.x,y=q_lakeflow.y,fill=lake_type))+
  #stat_poly_line()+
  geom_abline(aes(slope=1, intercept=0), lty=2, col='grey75')+
  geom_point(shape=21,col='black',size=1, stroke=.1)+
  geom_line(data=sen_df, aes(x=x,y=y, group=lake_type), lwd=1, col=wes_palette('Zissou1')[5])+
  geom_text(data=sen_vals, aes(x=1e-1, y=1e5, label=paste0('y=',signif(slope,2), 'x+', signif(int,2))), col=wes_palette('Zissou1')[5], size=5)+
  # geom_text(data=comb[!is.na(q_lakeflow.x)&!is.na(q_lakeflow.y),mean(abs(q_lakeflow.y-q_lakeflow.x)),by=lake_type], aes(x=1e-1, y=1e3, label=paste0('MAE=',signif(V1,2))), col=wes_palette('Zissou1')[5])+
  #stat_poly_eq(method=sen, use_label(c('EQ')))+
  scale_x_log10(labels=scales::comma)+
  scale_y_log10(labels=scales::comma)+
  scale_fill_manual(values=wes_palette('Zissou1')[c(2,4)])+
  facet_wrap(~lake_type)+
  xlab('Same-day Inflow (cms)')+
  ylab('Same-day Outflow (cms)')+
  theme_minimal()+
  theme(legend.position = 'none',
        legend.title = element_blank(),
        strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10))

ts_plot
ggsave('C:\\Users\\rriggs\\OneDrive - DOI\\Research\\LakeFlow_local\\out\\figures\\manuscript\\ts_figure.pdf', ts_plot, width=6.5, height=2, dpi=1000, units='in')
################################################################################
#CV Discharge and seasonal Q cdfs.
################################################################################
library(data.table)
library(dplyr)
library(ggplot2)
library(wesanderson)
min_obs=5
lf_outputs = fread('C:/Users/rriggs/OneDrive - DOI/Research/LakeFlow_local/out/stats/lf_outputs.csv')
ryan = lf_outputs
ryan = ryan[!is.na(q_lakeflow),list(norm_q = ((q_lakeflow-min(q_lakeflow))/(max(q_lakeflow)-min(q_lakeflow))), date=date, q_lakeflow=q_lakeflow,obs=.N), by=list(lake_id, reach_id, type, lake_type)]
ryan = ryan[obs>=min_obs]
ryan$quarter = lubridate::quarter(ryan$date)
ryan$month = lubridate::month(ryan$date)
ryan$season = 'Winter'
ryan$season = ifelse(ryan$month>=3&ryan$month<=5, 'Spring', ryan$season)
ryan$season = ifelse(ryan$month>=6&ryan$month<=8, 'Summer', ryan$season)
ryan$season = ifelse(ryan$month>=9&ryan$month<=11, 'Fall', ryan$season)
ryan$season = factor(ryan$season, levels=c('Total','Fall', 'Winter', 'Spring', 'Summer'))

seasonal_q=ggplot(ryan)+
  stat_ecdf(aes(x=norm_q, color=lake_type, lty=type))+
  scale_x_continuous(labels=c(0,0.25,0.5,0.75,1))+
  #scale_x_log10()+
  #coord_cartesian(xlim=c(0.1,1e4))+
  facet_wrap(~season, nrow=1)+
  scale_color_manual(values=wes_palette('Zissou1')[c(2,4)])+
  xlab('Normalized Discharge')+
  ylab('CDF')+
  theme_minimal()+
  theme(legend.position = 'none',
        legend.title = element_blank(),
        strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10))
seasonal_q

seasonal_q_box=ggplot(ryan)+
  geom_boxplot(aes(x=lake_type,y=norm_q, fill=lake_type, lty=type))+
  #scale_x_continuous(labels=c(0,0.25,0.5,0.75,1))+
  #scale_x_log10()+
  #coord_cartesian(xlim=c(0.1,1e4))+
  facet_wrap(~season, nrow=1,strip.position = 'bottom')+
  scale_fill_manual(values=wes_palette('Zissou1')[c(2,4)])+
  ylab('Normalized Discharge')+
  xlab('')+
  theme_minimal()+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10),
        strip.placement = 'outside')
seasonal_q_box
ggsave('C:\\Users\\rriggs\\OneDrive - DOI\\Research\\LakeFlow_local\\out\\figures\\manuscript\\seasonal_nq_box.png',seasonal_q_box,width=6.5,height=4,units='in',dpi=1000)


#Lakeflow applied across region
library(data.table)
library(dplyr)
library(ggplot2)
library(wesanderson)
lf_outputs = fread('C:/Users/rriggs/OneDrive - DOI/Research/LakeFlow_local/out/stats/lf_outputs.csv')
agg = lf_outputs[!is.na(q_lakeflow),list(meanQ=mean(q_lakeflow), sdQ=sd(q_lakeflow),obs=.N, mannings_n=unique(n_lakeflow), 
                                         bathymetry=unique(a0_lakeflow),rangeQ=max(q_lakeflow)-min(q_lakeflow)),by=list(lake_id, lake_type, reach_id, type)]
agg = agg[obs>=min_obs,]

mn = agg
mn = mn[,signif(median(mannings_n),2),by=list(lake_type)]
#CV calculation
cv_dt = agg[obs>=min_obs,list(cv = (sdQ/meanQ)*100, sd=sdQ), by=list(lake_id, lake_type, reach_id, type)]
cv_median = cv_dt[,median(cv),by=list(lake_type,type)]

cdf_plt = ggplot(cv_dt)+
  stat_ecdf(aes(x=sd,color=lake_type, lty=type))+
  #facet_wrap(~type)+
  coord_cartesian(xlim=c(1,2000))+
  scale_x_log10(labels=scales::comma_format(accuracy=1))+
  xlab('Standard deviation of discharge (cms)')+
  ylab('CDF')+
  scale_color_manual(values=wes_palette('Zissou1')[c(2,4)])+
  theme_minimal()+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10))
cdf_plt

# SWOT obs sd
library(data.table)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(ggh4x)
lf_outputs = fread('C:/Users/rriggs/OneDrive - DOI/Research/LakeFlow_local/out/stats/lf_outputs.csv')
min_obs = 5
lf_outputs = lf_outputs[!is.na(q_lakeflow),obs:=.N,by=list(reach_id)]
lf_outputs = lf_outputs[obs>=min_obs,]

swot_obs = lf_outputs[,c('wse', 'da', 'slope2', 'width', 'lake_type','type','reach_id','lake_id')]
swot_obs$width = exp(swot_obs$width)
swot_obs$slope2 = exp(swot_obs$slope2)
swot_obs = melt(swot_obs,id.vars=c('lake_id','reach_id','type','lake_type'))
swot_obs$variable = as.character(swot_obs$variable)
swot_obs$variable = ifelse(swot_obs$variable=='wse', 'Water surface elevation (m)', swot_obs$variable)
swot_obs$variable = ifelse(swot_obs$variable=='da', 'Change in cross-sectional area (m2)', swot_obs$variable)
swot_obs$variable = ifelse(swot_obs$variable=='slope2', 'Slope (m/m)', swot_obs$variable)
swot_obs$variable = ifelse(swot_obs$variable=='width', 'Width (m)', swot_obs$variable)
swot_obs$variable = factor(swot_obs$variable, levels=c('Water surface elevation (m)',
                                                       'Slope (m/m)',
                                                       'Width (m)',
                                                       'Change in cross-sectional area (m2)'))


swot_var_plot = ggplot(swot_obs[,sd(value),by=list(lake_id,reach_id,type,lake_type,variable)])+
  stat_ecdf(aes(x=V1, color=lake_type,lty=type))+
  facet_wrap(~variable, scales = 'free_x', nrow=1, strip.position='bottom',labeller = label_wrap_gen())+
  xlab('Standard deviation of SWOT river observables')+
  ylab('CDF')+
  scale_color_manual(values=wes_palette('Zissou1')[c(2,4)])+
  theme_minimal()+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10),
        strip.placement = 'outside')

swot_var_plot = swot_var_plot+facetted_pos_scales(x=list(
  scale_x_log10(limits=c(1e-2,4), labels=scales::label_comma()),
  scale_x_log10(limits=c(1e-6,1e-1), labels=scales::label_scientific()),
  scale_x_log10(limits=c(1,750), labels=scales::label_scientific()),
  scale_x_log10(limits=c(1e1,1e3),labels=scales::label_comma())
), y=NULL)

##Add in lake observations for seperate panel. 
swot_obs_lk = lf_outputs[,c('dv','lake_type','lake_id','date')]
swot_obs_lk = distinct(swot_obs_lk)

pull_lake_data = function(feature_id){
  website = paste0('https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries?feature=PriorLake&feature_id=',feature_id, '&start_time=2023-01-01T00:00:00Z&end_time=2024-12-31T00:00:00Z&output=csv&fields=lake_id,time_str,wse,area_total,xovr_cal_q,partial_f,dark_frac,ice_clim_f')
  response = GET(website)
  pull = content(response, as='parsed')$results
  data = try(read.csv(textConnection(pull$csv), sep=','))
  if(is.error(data)){return(NA)}
  data$reach_id = feature_id
  return(data)
}

batch_download_SWOT_lakes <- function(obs_ids){
  plan(multisession, workers = 20)
  SWOT_data = future_lapply(unique(obs_ids),pull_lake_data)
  plan(sequential)
  return(SWOT_data)
}

start_tm = Sys.time()
files_filt = batch_download_SWOT_lakes(as.character(unique(swot_obs_lk$lake_id)))
combined = rbindlist(files_filt[!is.na(files_filt)])
end_tm = Sys.time()
end_tm-start_tm

#Testing to see if partial flags make a difference since we're only using wse for lakes at the moment. - partial wse seems great. 
#lakeData = combined[combined$ice_clim_f<2&combined$dark_frac<=0.5&combined$xovr_cal_q<2&combined$partial_f==0,]
lakeData = combined[combined$ice_clim_f<2&combined$dark_frac<=0.5&combined$xovr_cal_q<2&combined$time_str!='no_data'&combined$wse>(5000*-1),]
lakeData = lakeData%>%distinct(.keep_all=TRUE)
lakeData$lake_id = as.character(lakeData$lake_id)

#FIXME remove lakedata with multiple ids, what's causing that??
lakeData$lake_id_first = sub(";.*", "", lakeData$lake_id)
lakeData$lake_id=lakeData$lake_id_first
lakeData$date = as.Date(lakeData$time_str)
lakeData$lake_id = as.character(lakeData$lake_id)
lakeData = data.table(lakeData)[,list(wse=mean(wse)),by=list(lake_id,date)]
swot_obs_lk$lake_id = as.character(swot_obs_lk$lake_id)

swot_obs_lk_jn = merge(swot_obs_lk, lakeData, by=c('lake_id','date'))
swot_obs_lk_tall = melt(swot_obs_lk_jn, measurement.vars=c('dv','wse'))
swot_obs_lk_tall$variable = as.character(swot_obs_lk_tall$variable)
swot_obs_lk_tall = swot_obs_lk_tall[variable!='date',]
swot_obs_lk_tall$variable = ifelse(swot_obs_lk_tall$variable=='dv', 'Storage change (cms)', swot_obs_lk_tall$variable)
swot_obs_lk_tall$variable = ifelse(swot_obs_lk_tall$variable=='wse', 'Water surface elevation (m)', swot_obs_lk_tall$variable)
swot_obs_lk_tall$variable = factor(swot_obs_lk_tall$variable, levels=c('Water surface elevation (m)', 'Storage change (cms)'))

swot_lk_var_plot = ggplot(swot_obs_lk_tall[,sd(value),by=list(lake_id,lake_type,variable)])+
  stat_ecdf(aes(x=V1, color=lake_type))+
  facet_wrap(~variable, scales = 'free_x', nrow=1, strip.position='bottom',labeller = label_wrap_gen())+
  xlab('Standard deviation of SWOT lake observables')+
  ylab('CDF')+
  scale_color_manual(values=wes_palette('Zissou1')[c(2,4)])+
  theme_minimal()+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10),
        strip.placement = 'outside')

swot_lk_var_plot = swot_lk_var_plot+facetted_pos_scales(x=list(
  scale_x_log10(limits=c(1e-2,10), labels=scales::label_scientific()),
  scale_x_log10(limits=c(1e-3,1e3),labels=scales::label_scientific())
), y=NULL)

swot_raw_plots = ggpubr::ggarrange(plotlist = list(swot_lk_var_plot, swot_var_plot), labels='auto',
                         common.legend = FALSE, nrow=2)
swot_raw_plots
ggsave('C:\\Users\\rriggs\\OneDrive - DOI\\Research\\LakeFlow_local\\out\\figures\\manuscript\\swot_obs_cdf_figure.png',swot_raw_plots, width=6.5, height=6.5, dpi=1000, units='in',bg='white')
################################################################################
# Correlation map. 
################################################################################
library(data.table)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(ggpmisc)
lf_outputs = fread('C:/Users/rriggs/OneDrive - DOI/Research/LakeFlow_local/out/stats/lf_outputs.csv')
agg = lf_outputs[,list(q_lakeflow = sum(q_lakeflow, na.rm=TRUE)), by=list(date, type,lake_id, lake_type)]
agg = agg[q_lakeflow>0,]
inflow = agg[type=='inflow',]
outflow = agg[type=='outflow',]
comb = merge(inflow, outflow, by=c('date', 'lake_id', 'lake_type'))
comb = comb[,OBS:=.N,by=lake_id]
comb = comb[OBS>=5,]
comb$diff_q = comb$q_lakeflow.x - comb$q_lakeflow.y
cor_dt = comb[,cor(q_lakeflow.x,q_lakeflow.y,method='spearman'),by=list(lake_id,lake_type)]
sd_diff = comb[,list(sd_in = sd(q_lakeflow.x), sd_out = sd(q_lakeflow.y)),by=list(lake_id,lake_type)]
sd_diff$diff = ((sd_diff$sd_out - sd_diff$sd_in)/sd_diff$sd_in)*100
sd_diff$diff = sd_diff$sd_out - sd_diff$sd_in


lf_outputs = fread('C:/Users/rriggs/OneDrive - DOI/Research/LakeFlow_local/out/stats/lf_outputs.csv')
min_obs = 5
lf_outputs = lf_outputs[!is.na(q_lakeflow),obs:=.N,by=list(reach_id)]
#lf_outputs = lf_outputs[obs>=min_obs,]
lf_outputs$latitude=ifelse(lf_outputs$Country=='Canada','High','Mid')
earliest = as.character(format(min(lf_outputs$date),'%m/%d/%Y'))
latest = as.character(format(max(lf_outputs$date),'%m/%d/%Y'))
total_lakes = length(unique(lf_outputs$lake_id[!is.na(lf_outputs$q_lakeflow)]))
total_obs = length(na.omit(lf_outputs$q_lakeflow))
total_reaches = lf_outputs[!is.na(q_lakeflow), length(unique(reach_id)),by=type]
total_lakeflow_eligible_na_lakes = 2565

lake_types = lf_outputs[!is.na(q_lakeflow),length(unique(lake_id)),by=lake_type]
obs_per = signif(mean(lf_outputs[!is.na(q_lakeflow),.N,by=reach_id]$N), 2)#lf_outputs[,.N,by=list(reach_id, lake_type,type)][,mean(N),by=list(lake_type,type)]
# Read in harmonized sword pld
sword_lake_db = 'C:/Users/rriggs/OneDrive - DOI/Research/LakeFlow_local/in/Intersected_SWORD_PLD_dataset_wo_ghost_rch.gdb'
layer_names= sf::st_layers(sword_lake_db)
sword_lake_db_files = layer_names$name[grep('lakes', layer_names$name)]

open_gdb = function(f){
  gdb = sf::st_read(sword_lake_db, layer=f, quiet=TRUE)
  gdb$continent = substr(gdb$lake_id, 1,1)
  gdb = gdb[gdb$U_reach_n>=0&gdb$D_reach_n>=0&gdb$continent%in%c('7', '8'),]#&as.character(gdb$lake_id)%in%lf_outputs$lake_id,]
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
updated_pld$type = jn$type[match(updated_pld$lake_id, jn$lake_id)]
# ha_db = 'C:\\Users\\rriggs\\OneDrive - DOI\\Research\\MISC\\BasinATLAS\\BasinATLAS_v10.gdb'
# layer_names_ha= sf::st_layers(ha_db)
# ha_b2 = st_read(ha_db, layer='BasinATLAS_v10_lev03')
# ha_b2 = st_make_valid(ha_b2)
# jn = st_join(updated_pld, ha_b2, st_within)
# updated_pld = jn

library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
world = ne_countries(scale = "medium", returnclass="sf")
na = world[world$continent=='North America',]

n_per = lf_outputs[,max(obs, na.rm=TRUE),by=lake_id]
updated_pld$observations = n_per$V1[match(updated_pld$lake_id, n_per$lake_id)]
cv_out = lf_outputs[type=='outflow',list(cv = (sd(q_lakeflow)/mean(q_lakeflow)*100)),by=list(lake_id, reach_id, type, lake_type)]
updated_pld$cv_out = cv_out$cv[match(updated_pld$lake_id, cv_out$lake_id)]
updated_pld$observations = ifelse(is.na(updated_pld$observations), 0, updated_pld$observations)
pts = st_centroid(updated_pld[updated_pld$U_reach_n>0&updated_pld$D_reach_n>0,])
pts$observations = ifelse(is.na(pts$observations), 0, pts$observations)
pts$obs = cut(pts$observations, breaks=c(-Inf,0,4,10,15,Inf))
pts = pts%>%arrange(obs)
pts$cor = cor_dt$V1[match(pts$lake_id, cor_dt$lake_id)]
pts$cor_cat = cut(pts$cor,breaks=c(-1,0,0.59,1))
pts = pts%>%arrange(cor_cat)
pts$diff = sd_diff$diff[match(pts$lake_id, sd_diff$lake_id)]
pts$diff_cap = cut(pts$diff,breaks=c(-Inf,-5,-1,1,5,Inf))
#pts$diff_cap = cut(pts$diff,breaks=c(-Inf,-100,-50,-10,10,50,100,Inf))

mp = ggplot(NULL)+
  geom_sf(data=st_union(na), fill='white', color='black')+
  #geom_sf(data=state_boundaries_wgs84, fill='white', color='grey50')+
  #geom_sf(data=pts, aes(color=obs, shape=type))+
  geom_sf(data=pts[is.na(pts$diff),],aes(shape=type),col='grey90')+
  geom_sf(data=pts[!is.na(pts$diff),], aes(fill=diff_cap, shape=type))+
  scale_shape_manual(values=c(21,24))+
  #scale_color_manual(values=c('grey80', 'black','#8F009A', '#F99B00', '#FB0032'))+
  #geom_sf(data=st_centroid(updated_pld), aes(color=cv_out, shape=type))+
  #scale_color_gradient2(low='orange', high='green', mid='yellow', midpoint=1)+
  scale_fill_brewer(palette='Spectral')+
  # scale_color_manual(values=c(wes_palette('Darjeeling1')[1],
  #                             wes_palette('Darjeeling1')[4],
  #                             wes_palette('Darjeeling1')[2]),na.value = 'grey90')+
  #scale_color_viridis_d()+
  theme_classic()+
  theme(axis.title=element_blank(),axis.text=element_blank(),
        axis.ticks=element_blank(), axis.line = element_blank(),
        legend.position = 'top',
        legend.title =element_text(size=11),
        legend.text=element_text(size=9))+
  coord_sf(crs ="+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45",
           xlim = c(-5000000, 3000000))+
  guides(fill = guide_legend(override.aes = list(shape = 21, 24)))
mp  
ggsave('C:\\Users\\rriggs\\OneDrive - DOI\\Research\\LakeFlow_local\\out\\figures\\manuscript\\sd_diff_map.png',mp,width=6.5,height=6.5,units='in',dpi=1000)
################################################################################
##
library(data.table)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(ggpmisc)
lf_outputs = fread('C:/Users/rriggs/OneDrive - DOI/Research/LakeFlow_local/out/stats/lf_outputs.csv')
agg = lf_outputs[,list(q_lakeflow = sum(q_lakeflow, na.rm=TRUE)), by=list(date, type,lake_id, lake_type)]
agg = agg[q_lakeflow>0,]
inflow = agg[type=='inflow',]
outflow = agg[type=='outflow',]
comb = merge(inflow, outflow, by=c('date', 'lake_id', 'lake_type'))
comb = comb[,OBS:=.N,by=lake_id]
comb = comb[OBS>=5,]
comb$diff_q = comb$q_lakeflow.x - comb$q_lakeflow.y

cor_spearman_all = comb[!is.na(q_lakeflow.x)&!is.na(q_lakeflow.y),list(R_spearman = cor(q_lakeflow.x, q_lakeflow.y,method='spearman'),n=.N), by=list(lake_id)]
cor_spearman_all_med = signif(median(cor_spearman_all$R_spearman[cor_spearman_all$n>=min_obs]), 2)
cor_spearman = comb[!is.na(q_lakeflow.x)&!is.na(q_lakeflow.y),list(R_spearman = cor(q_lakeflow.x, q_lakeflow.y,method='spearman'),n=.N), by=list(lake_id, lake_type)]

spearman_box = ggplot(cor_spearman)+
  geom_boxplot(aes(x=lake_type,y=R_spearman, fill=lake_type), outlier.alpha = 0.1)+
  ylab('RS')+
  xlab('')+
  scale_fill_manual(values=wes_palette('Zissou1')[c(2,4)])+
  theme_minimal()+
  theme(legend.position = 'none',
        legend.title = element_blank(),
        strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10),
        strip.placement = 'outside',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
spearman_box

ggsave('C:\\Users\\rriggs\\OneDrive - DOI\\Research\\LakeFlow_local\\out\\figures\\manuscript\\spearman_box.pdf',spearman_box,width=1.75,height=2.5,units='in',dpi=1000)



################################################################################
# CV and catchment area. 
################################################################################
#Lakeflow applied across region
library(data.table)
library(dplyr)
library(ggplot2)
library(wesanderson)
lf_outputs = fread('C:/Users/rriggs/OneDrive - DOI/Research/LakeFlow_local/out/stats/lf_outputs.csv')
# CV by catchment area at outlet. 
cat = 'C:\\Users\\rriggs\\OneDrive - DOI\\Research\\SWORD_PLD_Harmonization\\Harmonized_SWORD_PLD_dataset.gdb'
cat_names = st_layers(cat)
cat_files = cat_names$name[grep('catchments', cat_names$name)]

open_gdb_cat = function(f){
  gdb = sf::st_read(cat, layer=f, quiet=TRUE)
  gdb = gdb[gdb$In_lake_id%in%lf_outputs$lake_id|gdb$Otlt_lk_id%in%lf_outputs$lake_id,]#&as.character(gdb$lake_id)%in%lf_outputs$lake_id,]
  if(nrow(gdb)==0){return(NA)}
  return(gdb)
}

cat_gdb_list = lapply(cat_files, FUN=open_gdb_cat)
cat_gdb = data.table::rbindlist(cat_gdb_list[!is.na(cat_gdb_list)], use.names=TRUE,fill=TRUE)
cat_shp = st_as_sf(cat_gdb)


cv_dt$cat_area = cat_shp$Cat_a_tot[match(cv_dt$lake_id,cat_shp$Otlt_lk_id)]

ggplot(cv_dt)+
  geom_boxplot(aes(x=cut(cat_area,breaks=c(-Inf,1e3,1e4,1e5,Inf)), y=sd, fill=lake_type,lty=type))+
  coord_cartesian(ylim=c(0,2000))



ggplot(seasonal_cv[season=='Total'])+
  geom_point(aes(x=cat_area,y=cv,color=lake_type))+
  scale_x_log10()+
  facet_wrap(~type)






################################################################################
#FLP Plots
################################################################################
library(data.table)
library(dplyr)
library(ggplot2)
library(wesanderson)
lf_outputs = fread('C:/Users/rriggs/OneDrive - DOI/Research/LakeFlow_local/out/stats/lf_outputs.csv')
min_obs = 5
lf_outputs = lf_outputs[!is.na(q_lakeflow),obs:=.N,by=list(reach_id)]
lf_outputs = lf_outputs[obs>=min_obs,]
lf_outputs$latitude=ifelse(lf_outputs$Country=='CAN'|lf_outputs$States=='AK','High','Mid')
agg = lf_outputs[!is.na(q_lakeflow),list(meanQ=mean(q_lakeflow),medQ=median(q_lakeflow), sdQ=sd(q_lakeflow),obs=.N, mannings_n=unique(n_lakeflow), 
                                         bathymetry=unique(a0_lakeflow)),by=list(lake_id, lake_type, reach_id, type)]

bath_fig=ggplot(agg)+
  geom_boxplot(aes(fill=lake_type,x=lake_type,y=bathymetry))+
  facet_wrap(~type)+
  coord_cartesian(ylim=c(0,300))+ 
  ylab(expression(paste("Bathymetry (", m^{2},')')))+
  xlab('')+
  scale_fill_manual(values=wes_palette('Zissou1')[c(2,4)])+
  theme_minimal()+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        strip.text = element_text(size=14),
        legend.text = element_text(size=14),
        axis.title = element_text(size=14),
        axis.text = element_text(size=10))


man_fig=ggplot(agg)+
  geom_boxplot(aes(fill=lake_type,x=lake_type,y=mannings_n))+
  facet_wrap(~type)+
  ylab(expression(paste("n (s/", m^{1/3},')')))+
  xlab('')+
  coord_cartesian(ylim=c(0.01,0.075))+
  scale_fill_manual(values=wes_palette('Zissou1')[c(2,4)])+
  theme_minimal()+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        strip.text = element_text(size=14),
        legend.text = element_text(size=14),
        axis.title = element_text(size=14),
        axis.text = element_text(size=10))

flp_plt=ggpubr::ggarrange(bath_fig,man_fig, labels='auto', common.legend = TRUE, nrow=2)
flp_plt





library(data.table)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(ggh4x)
lf_outputs = fread('C:/Users/rriggs/OneDrive - DOI/Research/LakeFlow_local/out/stats/lf_outputs.csv')
min_obs = 5
lf_outputs = lf_outputs[!is.na(q_lakeflow),obs:=.N,by=list(reach_id)]
lf_outputs = lf_outputs[obs>=min_obs,]

swot_obs = lf_outputs[,c('wse', 'da', 'slope2', 'width', 'lake_type','type','reach_id','lake_id')]
swot_obs$width = exp(swot_obs$width)
swot_obs$slope2 = exp(swot_obs$slope2)
swot_obs = melt(swot_obs,id.vars=c('lake_id','reach_id','type','lake_type'))

swot_var_plot = ggplot(swot_obs[,sd(value),by=list(lake_id,reach_id,type,lake_type,variable)])+
  stat_ecdf(aes(x=V1, color=lake_type,lty=type))+
  facet_wrap(~variable, scales = 'free_x', nrow=1)+
  xlab('Standard deviation of SWOT observables')+
  ylab('CDF')+
  scale_color_manual(values=wes_palette('Zissou1')[c(2,4)])+
  theme_minimal()+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10))

swot_var_plot+facetted_pos_scales(x=list(
  scale_x_log10(limits=c(1e-2,5)),
  scale_x_log10(limits=c(1e1, 1e3)),
  scale_x_log10(limits=c(1e-6,1e-1)),
  scale_x_log10(limits=c(1,750))
), y=NULL)

swot_obs_sd = swot_obs[,sd(value),by=list(lake_id,reach_id,type,lake_type,variable)] 

ks.test(swot_obs_sd$V1[swot_obs_sd$variable=='slope2'&swot_obs_sd$type=='outflow'&swot_obs_sd$lake_type=='Reservoir'],swot_obs_sd$V1[swot_obs_sd$variable=='slope2'&swot_obs_sd$type=='outflow'&swot_obs_sd$lake_type=='Natural'])

##All sig for res outflow besides width. 

################################################################################
#Gamma fit. 
################################################################################
library(fitdistrplus)
library(TidyDensity)
ryan = lf_outputs[!is.na(q_lakeflow),list(norm_q = ((q_lakeflow-min(q_lakeflow))/(max(q_lakeflow)-min(q_lakeflow))),z_score=((q_lakeflow-mean(q_lakeflow))/sd(q_lakeflow)), date=date, q_lakeflow=q_lakeflow,obs=.N), by=list(lake_id, reach_id, type, lake_type)]
ryan = ryan[!is.nan(norm_q),]

util_gamma_param_estimate(ryan$norm_q[ryan$lake_type=='Reservoir'&ryan$type=='outflow'])$parameter_tbl[1,c("shape","scale","shape_ratio")]
util_gamma_param_estimate(ryan$norm_q[ryan$lake_type=='Natural'&ryan$type=='outflow'])$parameter_tbl[1,c("shape","scale","shape_ratio")]

#generate 1,000 random values that follow gamma distribution
res <- dgamma(shape=1.56, rate=5.90)
nat = dgamma(shape=1.42, rate=5.28)

ggplot(NULL)+
  geom_density(data=NULL,aes(x=res), col='red')+
  geom_density(data=NULL, aes(x=nat),col='black')





ggplot(ryan)+
  geom_density(aes(norm_q,..scaled.., color=lake_type, lty=type))+
  scale_color_manual(values=wes_palette('Zissou1')[c(2,4)])+
  theme_minimal()+theme(legend.position = 'top')
  
ggplot(ryan)+
  geom_density(aes(z_score,..scaled.., color=lake_type, lty=type))+
  scale_color_manual(values=wes_palette('Zissou1')[c(2,4)])+
  theme_minimal()+theme(legend.position = 'top')

ggplot(ryan)+
  stat_ecdf(aes(x=z_score, lty=type,color=lake_type))


ks.test(ryan$norm_q[ryan$type=='outflow'&ryan$lake_type=='Natural'],ryan$norm_q[ryan$type=='outflow'&ryan$lake_type=='Reservoir'], alternative='greater')


