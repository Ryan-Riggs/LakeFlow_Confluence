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
#Add in hydrolakes for capacity information. 
hydrolakes = st_read('C:\\Users\\rriggs\\OneDrive - DOI\\Research\\MISC\\HydroLAKES\\HydroLAKES_polys_v10_shp\\HydroLAKES_polys_v10.shp')
hydrolakes = hydrolakes[hydrolakes$Continent=='North America',]
hl_valid = st_is_valid(hydrolakes)
hydrolakes = st_make_valid(hydrolakes)
hl_jn = st_join(updated_pld, hydrolakes, st_nearest_feature)
#output data for rmd file. 
output = lf_outputs
output$lake_type = updated_pld$type[match(output$lake_id, updated_pld$lake_id)]
output$capacity = hl_jn$Vol_total[match(output$lake_id, hl_jn$lake_id)]
output$Country = hl_jn$Country[match(output$lake_id, hl_jn$lake_id)]
fwrite(output, paste0(inPath, 'out/stats/', 'lf_outputs.csv'))
