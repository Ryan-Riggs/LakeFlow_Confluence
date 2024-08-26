#**** Clean up code

# Create a variety of LakeFlow analysis figures. 
# By: Ryan Riggs
# Date: 7/23/2024
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

output_files = list.files(paste0(inPath, '/out/lf_results_na/'), full.names=TRUE)
lf_outputs = rbindlist(lapply(output_files, fread),fill=TRUE)
lf_outputs$Date = as.Date(lf_outputs$Date)
lf_outputs$reach_id = as.character(lf_outputs$reach_id)
lf_outputs$lake_id=as.character(lf_outputs$lake_id)
lf_outputs = lf_outputs[exp(lf_outputs$width)>=50]
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

ggplot(NULL)+
  geom_sf(data=state_boundaries_wgs84, fill='white', color='grey50')+
  geom_sf(data=st_centroid(updated_pld[!is.na(updated_pld$ratio),]), aes(color=ratioCap, size=Lake_area))+
  #scale_color_gradient2(low='orange', high='green', mid='yellow', midpoint=1)+
  scale_color_viridis_c()+
  theme_classic()+
  theme(axis.title=element_blank(),axis.text=element_blank(),
        axis.ticks=element_blank(), axis.line = element_blank(),
        legend.position = 'top',
        legend.title =element_text(size=11),
        legend.text=element_text(size=9))

# percentage of capacity, etc. 



# Linear relationship at each lake bw inflow and outflow? 
comb$lake_type = updated_pld$type[match(comb$lake_id, updated_pld$lake_id)]

#theil sen
sen <- function(..., weights = NULL) {
  mblm::mblm(...)
}
ggplot(comb[!is.na(q_lakeflow.x)&!is.na(q_lakeflow.y),],aes(x=q_lakeflow.x,y=q_lakeflow.y,color=lake_type))+
  geom_point()+
  stat_poly_line()+
  geom_abline(aes(slope=1, intercept=0), lty=2)+
  stat_poly_eq(method=sen, use_label(c('EQ',"R2")))+
  scale_x_log10()+
  scale_y_log10()+
  scale_color_manual(values=wes_palette('Zissou1')[c(2,4)])+
  facet_wrap(~lake_type)+
  xlab('Inflow (cms)')+
  ylab('Outflow (cms)')+
  theme_minimal()+
  theme(legend.position = 'top')

reg = comb[!is.na(q_lakeflow.x)&!is.na(q_lakeflow.y),list(r2=summary(lm((q_lakeflow.x)~(q_lakeflow.y)))$r.squared, n=.N), by=list(lake_id, lake_type)]
reg = reg[n>5,]
updated_pld$r2 = reg$r2[match(updated_pld$lake_id, reg$lake_id)]
# sample a few individual plots.
rn = sample(reg$lake_id, 10)
ggplot(comb[!is.na(q_lakeflow.x)&!is.na(q_lakeflow.y)&lake_id%in%rn,],aes(x=(q_lakeflow.x),y=(q_lakeflow.y)))+
  geom_point(aes(color=lake_type))+
  stat_poly_line()+
  geom_abline(aes(slope=1, intercept=0), lty=2)+
  stat_poly_eq(use_label(c("R2")))+
  facet_wrap(~lake_id, scales='free',nrow=2)

ggplot(reg)+
  stat_ecdf(aes(x=r2, color=lake_type))+
  scale_color_manual(values=wes_palette('Zissou1')[c(2,4)])+
  theme_minimal()

ggplot(reg)+
  geom_boxplot(aes(x=lake_type,y=r2, fill=lake_type))+
  scale_fill_manual(values=wes_palette('Zissou1')[c(2,4)])+
  theme_minimal()


agg$lake_type = updated_pld$type[match(agg$lake_id, updated_pld$lake_id)]
ggplot(agg)+
  stat_ecdf(aes(x=outflow_sd, color=lake_type))


ggplot(NULL)+
  geom_sf(data=state_boundaries_wgs84, fill='white', color='grey50')+
  geom_sf(data=st_centroid(updated_pld[!is.na(updated_pld$r2),]), aes(color=r2, size=Lake_area))+
  #scale_color_gradient2(low='orange', high='green', mid='yellow', midpoint=1)+
  scale_color_viridis_c()+
  theme_classic()+
  theme(axis.title=element_blank(),axis.text=element_blank(),
        axis.ticks=element_blank(), axis.line = element_blank(),
        legend.position = 'top',
        legend.title =element_text(size=11),
        legend.text=element_text(size=9))

# Does time of year effect difference between inflow and outflow?
ryan = lf_outputs
ryan$lake_type = updated_pld$type[match(ryan$lake_id, updated_pld$lake_id)]
all_vals = ryan[,list(cv=sd(q_lakeflow,na.rm=TRUE)/mean(q_lakeflow,na.rm=TRUE)), by=list(lake_id, lake_type, type, reach_id)]
all_vals$season = 'Total'
ryan$quarter = lubridate::quarter(ryan$date)
ryan$month = lubridate::month(ryan$date)
ryan$season = 'Winter'
ryan$season = ifelse(ryan$month>=3&ryan$month<=5, 'Spring', ryan$season)
ryan$season = ifelse(ryan$month>=6&ryan$month<=8, 'Summer', ryan$season)
ryan$season = ifelse(ryan$month>=9&ryan$month<=11, 'Fall', ryan$season)
ryan = ryan[!is.na(q_lakeflow),list( sd=sd(q_lakeflow), cv = sd(q_lakeflow,na.rm = TRUE)/mean(q_lakeflow,na.rm=TRUE), obs=.N, meanQ=mean(q_lakeflow,na.rm=TRUE), max=max(q_lakeflow,na.rm = TRUE), min=min(q_lakeflow,na.rm = TRUE)),by=list(lake_type,season, lake_id, type,reach_id)]
ryan$quarter = paste0('Q',ryan$quarter)
ryan = bind_rows(all_vals, ryan)
ryan$season = factor(ryan$season, levels=c('Total','Fall', 'Winter', 'Spring', 'Summer'))
labels = ryan[,round(median(cv*100,na.rm = TRUE)),by=list(quarter, lake_type,type)]
#ryan = ryan[!is.na(q_lakeflow),list( sd=sd(q_lakeflow), cv = sd(q_lakeflow,na.rm = TRUE)/mean(q_lakeflow,na.rm=TRUE), obs=.N),by=list(lake_type,lake_id, type,reach_id)]
g=ggplot(ryan)+
  geom_boxplot(aes(x=season,y=cv*100,fill=lake_type))+
  facet_wrap(~type)+
  coord_cartesian(ylim=c(0,200))+
  ylab('Coefficient of variation (%)')+
  xlab('Annual quarter')+
  scale_fill_manual(values=wes_palette('Zissou1')[c(2,4)])+
  theme_minimal()+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        strip.text = element_text(size=14),
        legend.text = element_text(size=14),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))
g

all_medians=ryan[,median(cv*100,na.rm = TRUE),by=list(lake_type, as.factor(quarter), type)]
nat_medians = all_medians[lake_type=='Natural']
res_medians = all_medians[lake_type=='Reservoir']
comb_medians = merge(nat_medians, res_medians, by=c('as.factor','type'))
comb_medians$difference=comb_medians$V1.y-comb_medians$V1.x



#ggsave(paste0(inPath,'out/figures/cv_seasonal.png'),g,bg='white',width=8,height=4,units='in',dpi=1000)

library(gganimate)
# map of cv
ryan = lf_outputs
ryan$lake_type = updated_pld$type[match(ryan$lake_id, updated_pld$lake_id)]
all_vals = ryan[,list(cv=sd(q_lakeflow,na.rm=TRUE)/mean(q_lakeflow,na.rm=TRUE)), by=list(lake_id, lake_type, type, reach_id)]
all_vals$season = 'Total'
ryan$quarter = lubridate::quarter(ryan$date)
ryan$month = lubridate::month(ryan$date)
ryan$season = 'Winter'
ryan$season = ifelse(ryan$month>=3&ryan$month<=5, 'Spring', ryan$season)
ryan$season = ifelse(ryan$month>=6&ryan$month<=8, 'Summer', ryan$season)
ryan$season = ifelse(ryan$month>=9&ryan$month<=11, 'Fall', ryan$season)
ryan$season = ifelse(ryan$season=='Winter'|ryan$season=='Fall', 'Fall/Winter', 'Spring/Summer')
ryan = ryan[!is.na(q_lakeflow),list( sd=sd(q_lakeflow), cv = sd(q_lakeflow,na.rm = TRUE)/mean(q_lakeflow,na.rm=TRUE), obs=.N, meanQ=mean(q_lakeflow,na.rm=TRUE), max=max(q_lakeflow,na.rm = TRUE), min=min(q_lakeflow,na.rm = TRUE)),by=list(lake_type,season, lake_id, type,reach_id)]
ryan$quarter = paste0('Q',ryan$quarter)
ryan = bind_rows(all_vals, ryan)
all_in = data.table(ryan)[obs>=3,.N,by=list(reach_id, season)]
all_in = all_in[,sum(N),by=reach_id]
all_in = all_in[V1>=3]
#ryan$season = factor(ryan$season, levels=c('Total','Fall', 'Winter', 'Spring', 'Summer'))
shp_new = updated_pld
ryan$geometry = shp_new$Shape[match(ryan$lake_id, shp_new$lake_id)]
ryan = st_sf(ryan, sf_column_name='geometry')

library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
world = ne_countries(scale = "small", returnclass="sf")
na = world[world$continent=='North America',]



ggplot(NULL)+
  geom_sf(data=na, fill='white', color='grey50')+
  geom_sf(data=st_centroid(ryan[!is.na(ryan$cv)&ryan$type=='outflow'&ryan$season!='Total'&ryan$reach_id%in%all_in$reach_id,]), aes(color=cut(cv*100,breaks=c(-Inf, 25, 50, 75, Inf)), shape=lake_type))+
  #scale_color_gradient2(low='orange', high='green', mid='yellow', midpoint=1)+
  scale_color_viridis_d()+
  theme_classic()+
  theme(axis.title=element_blank(),axis.text=element_blank(),
        axis.ticks=element_blank(), axis.line = element_blank(),
        legend.position = 'top',
        legend.title =element_text(size=11),
        legend.text=element_text(size=9))+
  coord_sf(crs ="+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45")+
  facet_wrap(~season+lake_type, ncol=2)

ryan = lf_outputs
ryan$lake_type = updated_pld$type[match(ryan$lake_id, updated_pld$lake_id)]
ryan$quarter = lubridate::quarter(ryan$date)
ryan = ryan[!is.na(q_lakeflow)&order(date),list(diff = diff(q_lakeflow), dt_diff = diff(as.numeric(date)),meanQ=mean(q_lakeflow)),by=list(lake_type,quarter, lake_id, type,reach_id)]
ryan$quarter = paste0('Q',ryan$quarter)

ggplot(ryan)+
  geom_boxplot(aes(x=as.factor(quarter),y=(abs(diff))/meanQ,fill=lake_type))+
  facet_wrap(~type)+
  coord_cartesian(ylim=c(0,5))+
  ylab('|dQ|/meanQ')+
  xlab('Annual quarter')+
  scale_fill_manual(values=wes_palette('Zissou1')[c(2,4)])+
  theme_minimal()+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        strip.text = element_text(size=14),
        legend.text = element_text(size=14),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))


ryan = lf_outputs
ryan$lake_type = updated_pld$type[match(ryan$lake_id, updated_pld$lake_id)]
ryan = ryan[!is.na(q_lakeflow),list( sd=sd(q_lakeflow), cv = sd(q_lakeflow,na.rm = TRUE)/mean(q_lakeflow,na.rm=TRUE), obs=.N, meanQ=mean(q_lakeflow,na.rm=TRUE), max=max(q_lakeflow,na.rm = TRUE), min=min(q_lakeflow,na.rm = TRUE)),by=list(lake_type,lake_id, type,reach_id)]
#ryan = ryan[!is.na(q_lakeflow),list( sd=sd(q_lakeflow), cv = sd(q_lakeflow,na.rm = TRUE)/mean(q_lakeflow,na.rm=TRUE), obs=.N),by=list(lake_type,lake_id, type,reach_id)]
g=ggplot(ryan)+
  geom_boxplot(aes(x=as.factor(lake_type),y=cv*100,fill=lake_type))+
  facet_wrap(~type)+
  coord_cartesian(ylim=c(0,200))+
  ylab('Coefficient of variation (%)')+
  xlab('Annual quarter')+
  scale_fill_manual(values=wes_palette('Zissou1')[c(2,4)])+
  theme_minimal()+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        strip.text = element_text(size=14),
        legend.text = element_text(size=14),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))
g

################################################################################
# correlation with lag times. 
################################################################################
output_files = list.files(paste0(inPath, '/out/lf_results_na/'), full.names=TRUE)
lf_outputs = rbindlist(lapply(output_files, fread),fill=TRUE)
lf_outputs$Date = as.Date(lf_outputs$Date)
lf_outputs$reach_id = as.character(lf_outputs$reach_id)
lf_outputs$lake_id=as.character(lf_outputs$lake_id)
lf_outputs = lf_outputs[exp(lf_outputs$width)>=50]
agg = lf_outputs[,list(q_lakeflow = sum(q_lakeflow, na.rm=TRUE)), by=list(date, type,lake_id)]
agg$q_lakeflow = ifelse(agg$q_lakeflow==0,NA, agg$q_lakeflow)
inflow = agg[type=='inflow',]
outflow = agg[type=='outflow',]
comb = merge(inflow, outflow, by=c('date', 'lake_id'))
comb = comb[,OBS:=.N,by=lake_id]
comb = comb[OBS>=5,]
comb$lake_type =updated_pld$type[match(comb$lake_id, updated_pld$lake_id)]

compList = list()
corList = list()
ryan=as.data.table(comb)
# create a new column for each difference in days - only considers consective differences though
ryan = ryan[order(date),]
#ryan = ryan[order(date),dt_diff := c(NA, diff(date)),by=lake_id]
ryan = ryan[order(lake_id),]
lks = unique(ryan$lake_id)
# trying to calculate additional differences. e.g., dt[3] - dt[1]
output_list = list()
for(k in 1:length(lks)){
sub = ryan[lake_id==lks[k]]
tab=list()
for(j in 2:nrow(sub)){
  current_dt = sub[j]
  previous_dt = sub[sub$date<current_dt$date]
  dt_df = current_dt$date - previous_dt$date
  inflow = previous_dt$q_lakeflow.x
  outflow = current_dt$q_lakeflow.y
  tab[[j]] = data.table(inflow = inflow, outflow=outflow, lake_id=sub$lake_id[1],lake_type=sub$lake_type[1], dt_diff = dt_df)
}
output_list[[k]] = rbindlist(tab)
}
all_diff = rbindlist(output_list)

#lag_times = c(0,7,14,30,60,90,120)#unique(all_diff$dt_diff)
lag_times = c(0,30, 60)
for(i in 2:length(lag_times)){
  lag_dt = all_diff[dt_diff<=lag_times[i]&dt_diff>lag_times[i-1],]
  lag_dt = lag_dt[!is.na(inflow)&!is.na(outflow),]
  cor_dt = lag_dt[,list(cor=cor(inflow, outflow,method = "spearman"),lag=lag_times[i], obs=.N), by=list(lake_type, lake_id)]
  corList[[i]] = cor_dt
  compList[[i]] = data.table(inflow=lag_dt$inflow, outflow=lag_dt$outflow, lag=lag_times[i],lake_type = lag_dt$lake_type, lake_id=lag_dt$lake_id)
}
lag = rbindlist(corList)
noLag = comb[,list(cor=cor(q_lakeflow.x, q_lakeflow.y,method = "spearman"),lag=0, obs=.N),by=list(lake_type,lake_id)]
lag = bind_rows(lag, noLag)
lag = lag[obs>=5,]
labels = paste0("\u2264", lag_times[-1])
labels = c(paste0(lag_times,"\u003C", 'Days', labels))
labels = c('Days=0',labels)
lag$area = updated_pld$Lake_area[match(lag$lake_id, updated_pld$lake_id)]
#lag$quarter = lubridate::quarter(lag$end_date)
c=ggplot(lag)+
  geom_boxplot(aes(x=as.factor(lag), y=cor, fill=lake_type))+
  ylab('Correlation between inflow and outflow (R value)')+
  xlab('Lag time between inflow and outflow (Days)')+
  scale_fill_manual(values=wes_palette('Zissou1')[c(2,4)])+
  theme_minimal()+
  scale_x_discrete(labels = labels)+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        strip.text = element_text(size=14),
        legend.text = element_text(size=14),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=315))
c
#ggsave(paste0(inPath,'out/figures/inflow_outflow_lag_correlation.png'),c,bg='white',width=8,height=4,units='in',dpi=1000)


#Plot map by highest correlation lag time. 
highest_cor = lag[,.SD[which.max(cor)],by=lake_id]
highest_cor = highest_cor[cor>0,]
updated_pld$lag = highest_cor$lag[match(updated_pld$lake_id, highest_cor$lake_id)]
updated_pld$cor = highest_cor$cor[match(updated_pld$lake_id, highest_cor$lake_id)]

ggplot(NULL)+
  geom_sf(data=state_boundaries_wgs84, fill='white', color='grey50')+
  geom_sf(data=st_centroid(updated_pld[!is.na(updated_pld$lag),]), aes(color=as.factor(lag), size=cor,shape=type))+
  #scale_color_gradient2(low='orange', high='green', mid='yellow', midpoint=1)+
  #scale_color_viridis_d()+
  scale_color_brewer(palette='Spectral')+
  theme_classic()+
  theme(axis.title=element_blank(),axis.text=element_blank(),
        axis.ticks=element_blank(), axis.line = element_blank(),
        legend.position = 'top',
        legend.title =element_text(size=11),
        legend.text=element_text(size=9))

nat_lks = nrow(highest_cor[lake_type=='Natural'])
res = nrow(highest_cor) - nat_lks
highest_cor$total = nat_lks
highest_cor$total=ifelse(highest_cor$lake_type=='Reservoir', res, nat_lks)
bar = highest_cor[cor>0,(.N/total)*100,by=list(lake_type,lag, total)]

b=ggplot(bar)+
  geom_bar(aes(x=as.factor(lag), y=V1, fill=lake_type), stat='identity',
           position='dodge')+
  ylab('Highest correlation (%)')+
  xlab('Lag time between inflow and outflow (Days)')+
  scale_fill_manual(values=wes_palette('Zissou1')[c(2,4)])+
  theme_minimal()+
  scale_x_discrete(labels = labels)+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        strip.text = element_text(size=14),
        legend.text = element_text(size=14),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=315))
b

#ggsave(paste0(inPath,'out/figures/highest_inflow_outflow_lag_correlation.png'),b,bg='white',width=8,height=4,units='in',dpi=1000)


comparison = bind_rows(rbindlist(compList), data.table(inflow=comb$q_lakeflow.x, outflow=comb$q_lakeflow.y, lake_type=comb$lake_type, lake_id=comb$lake_id))
comparison = merge(comparison, highest_cor, by=c('lake_id', 'lake_type', 'lag'))

highest_cor$lake_area = updated_pld$Lake_area[match(highest_cor$lake_id, updated_pld$lake_id)]

ggplot(highest_cor)+
  geom_point(aes(x=lake_area, y=lag, color=lake_type))+
  scale_x_log10()+
  facet_wrap(~lake_type)


library(tmap)
tmap_mode('view')

tm_shape(updated_pld)+
  tm_polygons()


################################################################################
# Normalized flow = do lakes vs res have higher or lower flows?
################################################################################
ryan = lf_outputs
ryan$lake_type = updated_pld$type[match(ryan$lake_id, updated_pld$lake_id)]

ryan = ryan[!is.na(q_lakeflow),list(norm_q = ((q_lakeflow-min(q_lakeflow))/(max(q_lakeflow)-min(q_lakeflow))), date=date), by=list(lake_id, reach_id, type, lake_type)]
ryan$quarter = lubridate::quarter(ryan$date)
ryan$month = lubridate::month(ryan$date)
ryan$season = 'Winter'
ryan$season = ifelse(ryan$month>=3&ryan$month<=5, 'Spring', ryan$season)
ryan$season = ifelse(ryan$month>=6&ryan$month<=8, 'Summer', ryan$season)
ryan$season = ifelse(ryan$month>=9&ryan$month<=11, 'Fall', ryan$season)
ryan$season = factor(ryan$season, levels=c('Total','Fall', 'Winter', 'Spring', 'Summer'))

ggplot(ryan[type=='outflow'])+
  stat_ecdf(aes(x=norm_q, color=lake_type))+
  facet_wrap(~season)

ggplot(ryan[!is.na(norm_q),mean(norm_q),by=list(month, lake_type, type)])+
  geom_line(aes(x=month, y=V1, color=lake_type))+
  facet_wrap(~type)

#ryan = ryan[!is.na(norm_q),list(cv = sd(norm_q)/mean(norm_q)),by=list(type, lake_type, season)]

ggplot(ryan[!is.na(norm_q), list(cv=sd(norm_q)/mean(norm_q)),by=list(type, lake_type, season)])+
  geom_bar(aes(x=season, y=cv*100, fill=lake_type),stat = 'identity', position='dodge')+
  facet_wrap(~type)


slopes = ryan[,list(q_diff=norm_q-data.table::shift(norm_q),dt_diff = date-data.table::shift(date), meanQ=mean(q_lakeflow,na.rm=TRUE)),by=list(reach_id, lake_type, type, season)]

ggplot(slopes[!is.na(dt_diff), mean(abs(q_diff/dt_diff)),by=list(reach_id, type,lake_type, season)])+
  geom_boxplot(aes(x=season, y=V1, fill=lake_type))+
  coord_cartesian(ylim=c(0,.5))+
  facet_wrap(~type)

################################################################################
# Map out errors in lf flow law 
################################################################################
lf_outputs$lake_type = updated_pld$type[match(lf_outputs$lake_id, updated_pld$lake_id)]
fl_error = lf_outputs[!is.na(q_lakeflow),list(total=sum(q_lakeflow), dv= unique(dv), lateral = mean(tributary), et=mean(et)),by=list(date, lake_id, type,lake_type)]
fl_in = fl_error[type=='inflow']
fl_out = fl_error[type=='outflow']
fl_comb = merge(fl_in,fl_out, by=c('date', 'lake_id', 'lake_type'))


fl_comb[,error:=dv.x-(total.x-total.y+lateral.x-et.x),]
fl_mean = fl_comb[,list(median_error=median(abs(error)),mean_error=mean(abs(error))),by=list(lake_id,lake_type)]

updated_pld$fl_error = fl_mean$median_error[match(updated_pld$lake_id, fl_mean$lake_id)]

library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
world = ne_countries(scale = "small", returnclass="sf")
na = world[world$continent=='North America',]

ggplot(NULL)+
  geom_sf(data=na, fill='white', color='grey50')+
  geom_sf(data=st_centroid(updated_pld[!is.na(updated_pld$fl_error),]), aes(color=cut(fl_error,breaks=c(0,1,10,25,100,1000,Inf)), shape=type))+
  #scale_color_gradient2(low='orange', high='green', mid='yellow', midpoint=1)+
  #scale_color_viridis_c(transform='log10')+
  scale_color_brewer(palette='Spectral',direction=-1)+
  theme_classic()+
  theme(axis.title=element_blank(),axis.text=element_blank(),
        axis.ticks=element_blank(), axis.line = element_blank(),
        legend.position = 'top',
        legend.title =element_text(size=11),
        legend.text=element_text(size=9))






ggplot(fl_comb)+
  stat_ecdf(aes(x=error,color=lake_type))+
  coord_cartesian(xlim=c(-1e4,1e4))






##
