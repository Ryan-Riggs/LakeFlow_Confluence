# Create a variety of LakeFlow validation figures. 
# By: Ryan Riggs
# Date: 7/11/2024
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

output_files = list.files(paste0(inPath, '/out/lf_results_na/'), full.names=TRUE)
lf_outputs = rbindlist(lapply(output_files, fread),fill=TRUE)
lf_outputs$Date = as.Date(lf_outputs$Date)
lf_outputs$reach_id = as.character(lf_outputs$reach_id)
lf_outputs$lake_id=as.character(lf_outputs$lake_id)
lf_outputs = lf_outputs[exp(lf_outputs$width)>=50,]
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

pull_gage_can = function(f){
  data = try(RivRetrieve::canada(site=f, variable='discharge'))
  if(is.error(data)){return(NA)}
  data$gage=f
  data=data[data$Date>=as.Date('2022-12-31'),]
  return(data)
}



gage_df_list = lapply(as.character(usgs_gages$site_no[usgs_gages$agency_cd=='USGS']), pull_gage_usgs)
gage_df = rbindlist(gage_df_list[!is.na(gage_df_list)])
gage_df_list_can = lapply(as.character(usgs_gages$site_no[usgs_gages$agency_cd=='Hydat']), pull_gage_can)
gage_df_can = rbindlist(gage_df_list_can[!is.na(gage_df_list_can)])

gage_df = bind_rows(gage_df, gage_df_can)
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
################################################################################
# 1:1 plot
################################################################################
validation$flow = lf_outputs$type[match(validation$reach_id, lf_outputs$reach_id)]
s=ggplot(validation)+
  geom_point(aes(x=Q,y=q_lakeflow))+
  #geom_point(aes(x=Q,y=exp(bayes_q)),col='black',size=0.5)+
  #geom_point(aes(x=Q,y=q_model), col='grey50')+
  geom_abline(aes(slope=1,intercept=0))+
  facet_wrap(~lake_type, scales='free')+
  scale_x_log10()+
  scale_y_log10()+
  xlab('Gage Discharge (cms)')+
  ylab('LakeFlow Discharge (cms)')+
  coord_cartesian(xlim=c(1,1000),ylim=c(1,1000))+
  theme_minimal()+
  theme(axis.title=element_text(size=18),axis.text=element_text(size=16),
        legend.position = 'none')
s
#ggsave(paste0(inPath,'out/figures/scatterplot.png'),s,bg='white', width=8,height=8,dpi=1000,units='in')
cor(validation$Q, validation$q_lakeflow, method='spearman')

################################################################################
# Hydrograph plots. 
################################################################################

#rep = c(error_stats$reach_id[sample(length(error_stats), 3)])
rep = c(error_stats$reach_id[error_stats$KGE==max(error_stats$KGE)],error_stats$reach_id[error_stats$KGE==quantile(error_stats$KGE,0.75)],error_stats$reach_id[error_stats$KGE==median(error_stats$KGE)],error_stats$reach_id[error_stats$KGE==quantile(error_stats$KGE,0.25)],error_stats$reach_id[error_stats$KGE==min(error_stats$KGE)])
vals = c(error_stats$KGE[error_stats$KGE==max(error_stats$KGE)],error_stats$KGE[error_stats$KGE==quantile(error_stats$KGE,0.75)],error_stats$KGE[error_stats$KGE==median(error_stats$KGE)],error_stats$KGE[error_stats$KGE==quantile(error_stats$KGE,0.25)],error_stats$KGE[error_stats$KGE==min(error_stats$KGE)])
vals = signif(vals, 2)
labels = c(paste0('Max KGE: ', vals[1]),paste0('75th KGE: ', vals[2]),paste0('Median KGE: ', vals[3]),paste0('25th KGE: ', vals[4]),paste0('Min KGE: ', vals[5]))
plot_df=full_gage_ts[full_gage_ts$reach_id%in%rep,]
plot_df$reach_id = factor(plot_df$reach_id, levels=rep)
levels(plot_df$reach_id) = labels#c('Max KGE (0.9)','75th KGE(0.5)','Median KGE (0.2)','25th KGE (-0.3)', 'Min KGE (-60)')
g=ggplot(plot_df[Date>as.Date('2023-10-01'),])+
  geom_line(aes(x=Date, y=Q), col='grey50')+
  #geom_line(aes(x=Date,y=q_lakeflow,col='red'))+
  #geom_point(aes(x=Date,y=exp(bayes_q)), col='black')+
  #geom_hline(aes(yintercept=q_model,col=reach_id))+
  geom_point(aes(x=Date,y=q_lakeflow),col='red')+
  #geom_line(aes(x=Date, y=q_model), col='grey50')+
  facet_wrap(~reach_id, scales='free_y',ncol=1)+
  scale_y_log10(n.breaks=4)+
  theme_minimal()+
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=16),
              legend.position = 'none',
        strip.text = element_text(size=18))+
  xlab('')+
  ylab('Discharge (cms)')
g

#ggsave(paste0(inPath, 'out/figures/hydrograph.png'),g,bg='white', width=12,height=6,dpi=1000,units='in')

################################################################################
# Summary error metrics plots. 
################################################################################
error_stats$type='LakeFlow'
error_stats_model$type='Model'

error_combined = bind_rows(error_stats, error_stats_model)
error_combined$'|rBias|' = abs(error_combined$rBias)
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
c=ggplot(error_tall[variable%in%c('rBias')])+
  stat_ecdf(aes(x=value,color=type),lwd=1.5)+
  facet_wrap(~variable)+
  coord_cartesian(xlim=c(-100,100))+
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
pt

#ggsave(paste0(inPath, 'out/figures/error_cdfs.png'),pt,bg='white', width=12,height=6,dpi=1000,units='in')


# Boxplots. 
ggplot(error_tall[variable%in%c('rBias')])+
  geom_boxplot(aes(y=value,fill=type))+
  coord_cartesian(ylim=c(-200,200))


ggplot(error_tall[variable%in%c('RRMSE','NRMSE','|rBias|'),])+
  geom_boxplot(aes(y=value,fill=type))+
  facet_wrap(~variable, scales='free')+
  coord_cartesian(ylim=c(0,300))

ggplot(error_tall[variable%in%c('KGE','NSE', 'Rvalue'),])+
  geom_boxplot(aes(y=value,fill=type))+
  facet_wrap(~variable)+
  coord_cartesian(ylim=c(-2,1))
################################################################################
# Map of gage locations
################################################################################
shp = usgs_gages[usgs_gages$reach_id%in%validation$reach_id,]
shp = st_as_sf(shp, coords=c('dec_long_va','dec_lat_va'))
shp = left_join(shp, error_stats, by='reach_id')
data("state_boundaries_wgs84")
st_crs(state_boundaries_wgs84)
st_crs(shp) = crs(state_boundaries_wgs84)
state_boundaries_wgs84 = state_boundaries_wgs84[state_boundaries_wgs84$NAME%!in%c('Puerto Rico', 'Alaska', 'Hawaii','U.S. Virgin Islands'),]
shp$cap = cut(shp$KGE, breaks=c(-Inf, -.41, 0, .5, 1))
shp$capKGE = ifelse(shp$KGE<(1*-1), -1, shp$KGE)
#pal = rev(wes_palette('Zissou1', 100, 'continuous'))
pal = wes_palette('Zissou1')

mp = ggplot(NULL)+
  geom_sf(data=state_boundaries_wgs84, color='grey75', fill='grey50')+
  geom_sf(data=shp,col='black', lwd=1,size=2)+
  geom_sf(data=shp, aes(color=cap), lwd=1,size=1.75)+
  #scale_color_gradientn(colours=pal)+
  #scale_color_manual(values=pal)+
  #scale_color_gradient2(low='darkorange',mid='white',high='green',midpoint=-.41)+
  #scale_color_gradient2(low=pal[3],mid='white',high=pal[4],midpoint=-.41)+
  scale_color_manual(values=c(pal[5],'white',pal[2],pal[1]))+
  theme_classic()+
  theme(axis.title=element_blank(),axis.text=element_blank(),
        axis.ticks=element_blank(), axis.line = element_blank(),
        legend.title =element_text(size=11),
        legend.text=element_text(size=9))
mp
#ggsave(paste0(inPath, '/out/figures/map_gage.png'),mp,bg='white', width=12,height=8,dpi=1000,units='in')


library(tmap)
tmap_mode('view')

shp$meanQ = validation$meanQ[match(shp$site_no, validation$gage)]
tm_shape(shp[,c('site_no','capKGE','meanQ')])+
  tm_dots(col='capKGE')
################################################################################
# Map of lakes - currently rereading everything in . 
################################################################################
output_files = list.files(paste0(inPath, '/out/lf_results_na/'), full.names=TRUE)
lf_outputs = rbindlist(lapply(output_files, fread),fill=TRUE)
lf_outputs$Date = as.Date(lf_outputs$Date)
lf_outputs$reach_id = as.character(lf_outputs$reach_id)
lf_outputs$lake_id=as.character(lf_outputs$lake_id)
lf_outputs = lf_outputs[exp(lf_outputs$width)>=50,]

usgs_gages = fread('C:\\Users\\rriggs\\OneDrive - DOI\\Research\\SWOT\\gage_usgs_can_sword_snapped.csv', colClasses = as.character('site_no'))
shp = usgs_gages[usgs_gages$reach_id%in%lf_outputs$reach_id,]
shp = st_as_sf(shp, coords=c('dec_long_va','dec_lat_va'))
shp$flow = lf_outputs$type[match(shp$reach_id, lf_outputs$reach_id)]
data("state_boundaries_wgs84")
st_crs(state_boundaries_wgs84)
st_crs(shp) = crs(state_boundaries_wgs84)
state_boundaries_wgs84 = state_boundaries_wgs84[state_boundaries_wgs84$NAME%!in%c('Puerto Rico', 'Alaska', 'Hawaii','U.S. Virgin Islands'),]
#pal = rev(wes_palette('Zissou1', 100, 'continuous'))
shp$lake_id = lf_outputs$lake_id[match(shp$reach_id, lf_outputs$reach_id)]


# Read in harmonized sword pld
sword_lake_db = 'C:/Users/rriggs/OneDrive - DOI/Research/LakeFlow_local/in/Intersected_SWORD_PLD_dataset.gdb'
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
updated_pld = updated_pld[updated_pld$continent=='7',]
updated_pld = left_join(updated_pld, shp, by='lake_id')
updated_pld = st_sf(updated_pld)

# combine where there are both inflow and outflow gages.  
un_lks = unique(updated_pld$lake_id)
lks = list()
for(i in 1:length(un_lks)){
  sub = updated_pld[updated_pld$lake_id==un_lks[i],]
  flows = paste(unique(sub$flow),collapse=' and ')
  lks[[i]] = data.table(lake_id=un_lks[i],type=flows)
}
lks = rbindlist(lks)
updated_pld$flow = lks$type[match(updated_pld$lake_id, lks$lake_id)]
lake_shp = st_centroid(updated_pld)
lake_shp$flow = factor(lake_shp$flow, levels=c('inflow', 'outflow', 'inflow and outflow', 'NA'))

mp = ggplot(NULL)+
  geom_sf(data=state_boundaries_wgs84, fill='white', color='grey50')+
  geom_sf(data=lake_shp, aes(color=flow, size=Lake_area))+
  scale_color_manual(values=c('#6AC3F2','#FFBD00','#66A48B', '#A79D95'))+
  theme_classic()+
  theme(axis.title=element_blank(),axis.text=element_blank(),
        axis.ticks=element_blank(), axis.line = element_blank(),
        legend.position = 'top',
        legend.title =element_text(size=11),
        legend.text=element_text(size=9))
mp
#ggsave(paste0(inPath, '/out/figures/map_lakes.png'),mp,bg='white', width=12,height=8,dpi=1000,units='in')
################################################################################
# Performance vs River size. 
################################################################################
ggplot(error_tall[type=='LakeFlow'&variable%in%c('RRMSE', 'NRMSE', '|rBias|'),])+
  geom_point(aes(x=meanQ, y=value))+
  facet_wrap(~variable, scales='free')+
  coord_cartesian(ylim=c(0,200))


ggplot(error_tall[type=='LakeFlow'&variable%in%c('KGE', 'NSE', 'Rvalue'),])+
  geom_point(aes(x=meanQ, y=value))+
  facet_wrap(~variable, scales='free')+
  coord_cartesian(ylim=c(-2,1))+
  scale_x_log10()
  
################################################################################
# Performance vs inflow/outflow. 
################################################################################
output_files = list.files(paste0(inPath, '/out/lf_results_geoglows/'), full.names=TRUE)
lf_outputs_v2 = rbindlist(lapply(output_files, fread),fill=TRUE)
lf_outputs_v2$Date = as.Date(lf_outputs_v2$Date)
lf_outputs_v2$reach_id = as.character(lf_outputs_v2$reach_id)
lf_outputs_v2$lake_id=as.character(lf_outputs_v2$lake_id)


error_tall$flow = lf_outputs_v2$type[match(error_tall$reach_id, lf_outputs_v2$reach_id)]

ggplot(error_tall[type=='LakeFlow'&variable%in%c('RRMSE','NRMSE','|rBias|'),])+
  geom_boxplot(aes(x=flow,y=value,fill=type))+
  facet_wrap(~variable, scales='free')+
  coord_cartesian(ylim=c(0,300))

ggplot(error_tall[type=='LakeFlow'&variable%in%c('KGE','NSE', 'Rvalue'),])+
  geom_boxplot(aes(x=flow,y=value,fill=type))+
  facet_wrap(~variable)+
  coord_cartesian(ylim=c(-10,1))

################################################################################
# Calculate some stats. 
################################################################################

#Positive KGE
error_tall[variable=='KGE'&value>=(-0.41),(.N/length(unique(error_tall$reach_id)))*100,by=type]

#Positive NSE
error_tall[variable=='NSE'&value>0,(.N/length(unique(error_tall$reach_id)))*100,by=type]

#NRMSE
error_tall[variable=='NRMSE'&value<50,(.N/length(unique(error_tall$reach_id)))*100,by=type]


#Nbias
error_tall[variable=='rBias'&abs(value)<50,(.N/length(unique(error_tall$reach_id)))*100,by=type]

#
diff_df = merge(error_stats, error_stats_model, by='reach_id')

################################################################################
# AGU abstract stats. 
################################################################################
agu = merge(error_stats, error_stats_model, by='reach_id')
agu$width = lf_outputs$width[match(agu$reach_id, lf_outputs$reach_id)]
agu$KGE_diff = agu$KGE.x - agu$KGE.y
agu$KGE_pchange = (agu$KGE_diff/abs(agu$KGE.y))*100
agu$abs_rbias_diff = abs(agu$rBias.x)-abs(agu$rBias.y)
agu$abs_rbias_pchange = (agu$abs_rbias_diff/abs(agu$rBias.y))*100

median(agu$KGE.x)
median(abs(agu$rBias.x))

#Percent of KGE improvement over the model and how much it improves by
nrow(agu[agu$KGE_diff>0])/nrow(agu)
mean(agu$KGE_diff[agu$KGE_diff>0])
mean(agu$KGE_pchange[agu$KGE_diff>0])
median(agu$KGE_pchange)
median(agu$KGE_diff)


#Percent of abs nbias improvement over the model and how much it improves by
nrow(agu[agu$abs_rbias_diff<0])/nrow(agu)
mean(agu$abs_rbias_diff[agu$abs_rbias_diff<0])
mean(agu$abs_rbias_pchange[agu$abs_rbias_diff<0])
median(agu$abs_rbias_diff)

agu$width_grp = cut(exp(agu$width), breaks=c(0,100,Inf))
agu[KGE_diff>0,.N,by=width_grp]
agu[abs_rbias_diff<0, .N,by=width_grp]




################################################################################
# Break down by width groupings. 
################################################################################
error_tall$width = lf_outputs$width[match(error_tall$reach_id, lf_outputs$reach_id)]
error_tall$width = exp(error_tall$width)
error_tall$width_grp = cut(error_tall$width, breaks=c(0,100,Inf))
error_tall$value  = ifelse(error_tall$variable=='KGE', error_tall$value*100, error_tall$value)

ggplot(error_tall[variable%in%c('|rBias|', 'KGE', 'NRMSE')])+
  geom_boxplot(aes(x=width_grp, y=value, fill=type))+
  facet_wrap(~variable)+
  coord_cartesian(ylim=c(-100,250))+
  theme_minimal()



error_tall[]

################################################################################
# reservoir vs natural lake performance.
################################################################################
error_stats$lake_type = lake_reach$lake_type[match(error_stats$reach_id, lake_reach$reach_id)]


cl = wes_palette('AsteroidCity3')[3:4]
cl = rev(cl)

# cdf plots. 
a=ggplot(error_stats)+
  stat_ecdf(aes(x=NRMSE,color=lake_type),lwd=1.5)+
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
  scale_color_manual(values=wes_palette('Zissou1')[c(2,4)])
a
b=ggplot(error_stats)+
  stat_ecdf(aes(x=KGE,color=lake_type),lwd=1.5)+
  coord_cartesian(xlim=c(-1,1))+
  theme_minimal()+
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=16),
        legend.position = 'top',
        legend.title = element_blank(),
        strip.text = element_text(size=18),
        legend.text = element_text(size=18))+
  ylab('CDF')+
  scale_color_manual(values=wes_palette('Zissou1')[c(2,4)])
b
c=ggplot(error_stats)+
  stat_ecdf(aes(x=abs(rBias),color=lake_type),lwd=1.5)+
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
  scale_color_manual(values=wes_palette('Zissou1')[c(2,4)])
c
pt = ggarrange(a,c,b, labels='auto',nrow=1,common.legend = TRUE)
pt


error_stats[,median(KGE), by=lake_type]
error_stats[,median(NRMSE), by=lake_type]
error_stats[,median(Rvalue), by=lake_type]
error_stats[,median(abs(rBias)), by=lake_type]


error_stats$lake_id = lake_reach$lake_id[match(error_stats$reach_id, lake_reach$reach_id)]
error_stats$lake_area = updated_pld$Lake_area[match(error_stats$lake_id, updated_pld$lake_id)]

################################################################################
## Using LakeFlow estimates to create AHG - can LF extend beyond SWOT? 
################################################################################
ggplot(validation,aes(x=wse,y=Q, color=lake_type))+
  geom_point()+
  geom_smooth(method='lm',se=FALSE)+
  scale_x_log10()+
  scale_y_log10()+
  facet_wrap(~gage,scales='free')

gages = unique(validation$gage)
tab=list()
for(i in 1:length(gages)){
  sub = validation[gage==gages[i]]
  reg = lm(log(q_lakeflow)~log(wse),data=sub)
  pred = predict(reg, sub)
  tab[[i]] = data.table(lf = exp(pred), Q=sub$Q,raw_lf=sub$q_lakeflow, gage=gages[i],lake_type=sub$lake_type[1],date=sub$Date)
}
out = rbindlist(tab)
median(out[,hydroGOF::KGE(lf, Q),by=gage]$V1)

out[,hydroGOF::KGE(lf, Q),by=list(gage,lake_type)][,median(V1),by=lake_type]
out[,hydroGOF::KGE(raw_lf, Q),by=list(gage,lake_type)][,median(V1),by=lake_type]

ggplot(out)+
  geom_line(aes(x=date, y=Q))+
  geom_point(aes(x=date,y=raw_lf),col='blue')+
  geom_point(aes(x=date,y=lf),col='red')+
  facet_wrap(~gage, scales='free')



