##LakeFlow code for running locally
##By: Ryan Riggs and George Allen, June 2024.

################################################################################
# set Path to Lakeflow_local folder. 
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
library(geoBAMr)
library(future)
library(future.apply)
'%!in%' <- function(x,y)!('%in%'(x,y))
################################################################################
# Read in relevant files: Harmonized sword-pld, reservoirs of interest, swot lake data
################################################################################
updated_pld = fread(paste0(inPath,"in/SWORDv16_PLDv103_wo_ghost_rch.csv"))
updated_pld$lake_id =  as.character(updated_pld$lake_id)
updated_pld$continent = substr(updated_pld$lake_id, 1,1)
################################################################################
# Read in lake data via hydrocron. 
################################################################################
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

files_filt = batch_download_SWOT_lakes(updated_pld$lake_id[updated_pld$continent%in%c('7', '8')][1:100])
combined = rbindlist(files_filt[!is.na(files_filt)])
################################################################################
# Filter lake data. 
################################################################################
#Testing to see if partial flags make a difference since we're only using wse for lakes at the moment. - partial wse seems great. 
#lakeData = combined[combined$ice_clim_f<2&combined$dark_frac<=0.5&combined$xovr_cal_q<2&combined$partial_f==0,]
lakeData = combined[combined$ice_clim_f<2&combined$dark_frac<=0.5&combined$xovr_cal_q<2&combined$time_str!='no_data'&combined$wse>(5000*-1),]
lakeData = lakeData%>%distinct(.keep_all=TRUE)
lakeData$lake_id = as.character(lakeData$lake_id)

#Remove lakedata with multiple ids
lakeData$lake_id_first = sub(";.*", "", lakeData$lake_id)
lakeData$lake_id=lakeData$lake_id_first
lakeData = lakeData[lakeData$lake_id%in%updated_pld$lake_id,]
rm(combined)
################################################################################
# Function to pull SWOT reach data. 
################################################################################
pull_data = function(feature_id){
  website = paste0('https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries?feature=Reach&feature_id=',feature_id, '&start_time=2023-01-01T00:00:00Z&end_time=2024-12-31T00:00:00Z&output=csv&fields=reach_id,time_str,wse,width,slope,slope2,d_x_area,area_total,reach_q,p_width,xovr_cal_q,partial_f,dark_frac,ice_clim_f,wse_r_u,slope_r_u,reach_q_b')
  response = GET(website)
  pull = content(response, as='parsed')$results
  data = try(read.csv(textConnection(pull$csv), sep=','))
  if(is.error(data)){return(NA)}
  data$reach_id = feature_id
  return(data)
}
################################################################################
# Tukey test for removing outliers. 
################################################################################
tukey_test = function(ts){
  wseIQR = quantile(ts$wse, c(.25, .75))
  wseT_l = wseIQR[1] - (diff(wseIQR)*1.5)
  wseT_u = wseIQR[2] + (diff(wseIQR)*1.5)
  
  slopeIQR = quantile(ts$slope2, c(.25, .75))
  slopeT_l = slopeIQR[1] - (diff(slopeIQR)*1.5)
  slopeT_u = slopeIQR[2] + (diff(slopeIQR)*1.5)
  
  widthIQR = quantile(ts$width, c(.25, .75))
  widthT_l = widthIQR[1] - (diff(widthIQR)*1.5)
  widthT_u = widthIQR[2] + (diff(widthIQR)*1.5)
  
  ts_filt = ts[ts$width>=widthT_l&ts$width<=widthT_u&ts$wse>=wseT_l&ts$wse<=wseT_u&ts$slope2>=slopeT_l&ts$slope2<=slopeT_u,]
  return(ts_filt)
}

tukey_test_lake = function(ts){
  wseIQR = quantile(ts$wse, c(.25, .75))
  wseT_l = wseIQR[1] - (diff(wseIQR)*1.5)
  wseT_u = wseIQR[2] + (diff(wseIQR)*1.5)
  
  ts_filt = ts[ts$wse>=wseT_l&ts$wse<=wseT_u,]
  return(ts_filt)
}
################################################################################
# Function to filter SWOT reach data. 
################################################################################
filter_function = function(swot_ts){
  # Allowing partial obs to see if it degrades performances. 
  #dawg_filter = tukey_test(swot_ts[swot_ts$time_str!='no_data'&swot_ts$ice_clim_f<2&swot_ts$dark_frac<=0.5&swot_ts$xovr_cal_q<2&!is.na(swot_ts$slope2)&!is.infinite(swot_ts$slope2)&swot_ts$slope2>0&swot_ts$width>0,])
  dawg_filter = tukey_test(swot_ts[swot_ts$time_str!='no_data'&swot_ts$ice_clim_f<2&swot_ts$dark_frac<=0.5&swot_ts$xovr_cal_q<2&swot_ts$partial_f==0&!is.na(swot_ts$slope2)&!is.infinite(swot_ts$slope2)&swot_ts$slope2>0&swot_ts$width>0,])
  qual_filter = swot_ts[swot_ts$time_str!='no_data'&swot_ts$reach_q<=2,]
  ssf_filter = tukey_test(swot_ts[swot_ts$time_str!='no_data'&swot_ts$reach_q_b<=32768&swot_ts$dark_frac<=0.1&swot_ts$wse_r_u<=0.5&swot_ts$slope_r_u<=10e-5&swot_ts$ice_clim_f==0&swot_ts$xovr_cal_q<=1,])

  #return(qual_filter)
  return(dawg_filter)
  #return(ssf_filter)
}
################################################################################
# Pull in relevant SWOT river reach data and filter using above functions. 
################################################################################
batch_download_SWOT <- function(obs_ids){
  plan(multisession, workers = 20)
  SWOT_data = future_lapply(unique(obs_ids),pull_data)
  plan(sequential)
  return(SWOT_data)
}

# Filter to lakes with at least n observations. 
n = 5
lakes = unique(data.table(lakeData)[,.N,by=lake_id][N>=n]$lake_id)
lakeFilt = lakeData[lake_id%in%lakes,tukey_test_lake(.SD),by=lake_id]
lakes = unique(data.table(lakeFilt)[,.N,by=lake_id][N>=n]$lake_id)
lakeFilt = lakeFilt[lake_id%in%lakes,]
up_reaches = unlist(strsplit(updated_pld$U_reach_id[updated_pld$lake_id%in%lakes], ','))
dn_reaches = unlist(strsplit(updated_pld$D_reach_id[updated_pld$lake_id%in%lakes], ','))
reaches = c(up_reaches, dn_reaches)
swot_river_pull = batch_download_SWOT(reaches)
swot_river_filt = lapply(swot_river_pull[!is.na(swot_river_pull)], filter_function)
swot_river = rbindlist(swot_river_filt)
swot_river$time=as_datetime(swot_river$time_str)

# Work on adding in the date filtering approach. This step will help improve the functions speed. 
combining_lk_rv_obs = function(lake){
  #Pull in SWOT river data and subset predownloaded SWOT lake data. 
  upID = unlist(strsplit(updated_pld$U_reach_id[updated_pld$lake_id==lake], ','))
  dnID = unlist(strsplit(updated_pld$D_reach_id[updated_pld$lake_id==lake], ','))
  upObs_all = swot_river[swot_river$reach_id%in%upID,]
  dnObs_all = swot_river[swot_river$reach_id%in%dnID,]
  lakeObs_all = lakeFilt[lakeFilt$lake_id==lake,]
  lakeObs_all$time = as_datetime(lakeObs_all$time_str)
  
  # FIXME: changing lake areas to pld mean lake areas due to SWOT errors. 
  prior_area = updated_pld$Lake_area[updated_pld$lake_id==lake]
  lakeObs_all$area_total = prior_area
  
  if(nrow(lakeObs_all)<3){return(NA)}
  
  ################################################################################
  # Get dates in proper format and subset to matching dates. 
  ################################################################################
  lakeObs_all$date = as.Date(lakeObs_all$time)
  upObs_all$date = as.Date(upObs_all$time)
  dnObs_all$date = as.Date(dnObs_all$time)
  
  # FIXME: Aggregating lakes to mean values for multiple observations in one day. 
  lakeObs = data.table(lakeObs_all)[,c('wse', 'area_total', 'date')][,lapply(.SD, mean), by=date]
  upObs = data.table(upObs_all)[,c('wse', 'width', 'slope', 'slope2','reach_id', 'date')][,lapply(.SD, mean), by=list(date, reach_id)]
  dnObs = data.table(dnObs_all)[,c('wse', 'width', 'slope', 'slope2','reach_id', 'date')][,lapply(.SD, mean), by=list(date, reach_id)]
  
  lkDates = unique(lakeObs$date)
  upDts = upObs[,.N,by=date][N>=length(upID)] # limit to dates with obs for each upstream reach.
  dnDts = dnObs[,.N,by=date][N>=length(dnID)] # limit to dates with obs for each downstream reach. 
  
  #goodDates = lkDates[lkDates%in%upObs_all$date&lkDates%in%dnObs_all$date]
  goodDates = lkDates[lkDates%in%upDts$date&lkDates%in%dnDts$date]
  
  lakeObsGood = lakeObs[lakeObs$date%in%goodDates,]
  upObsGood = upObs[upObs$date%in%goodDates,]
  dnObsGood = dnObs[dnObs$date%in%goodDates,]
  
  lakeObs = lakeObsGood[order(lakeObsGood$date),]
  upObs = upObsGood[order(upObsGood$date),]
  dnObs = dnObsGood[order(dnObsGood$date),]
  
  upObs = upObs[order(upObs$reach_id),]
  dnObs = dnObs[order(dnObs$reach_id),]
  
  if(nrow(lakeObs)<4){return(NA)}
  output = list(lakeObs, upObs, dnObs)
  return(output)
}

# subset to lakes with enough lake and river SWOT obs to run LakeFlow. 
viable_data = lapply(lakes, combining_lk_rv_obs)
names(viable_data) = lakes
n_obs_lake = data.table(lake=lakes,obs=unlist(lapply(viable_data, length)))
viable_locations = n_obs_lake[obs>=3,] #Note this 3 is for ensuring an inflow, lake, and outlfow have viable data. 
################################################################################
# Function to pull tributary inflow Q estimates from Geoglows. 
################################################################################
start_date = '01-01-2023'

download_tributary = function(reaches, start_date=start_date){
  source_python(paste0(inPath, 'src/geoglows_aws_pull.py'))
  tributary_flow = pull_tributary(reach_id = reaches, start_date=start_date)
  tributary_flow$date = as.Date(row.names(tributary_flow))
  n_col = ncol(tributary_flow)
  tributary_aggregated = data.table(tributary_flow)[,tributary_total:=rowSums(.SD), .SDcols=-n_col]
  return(tributary_aggregated)
}
################################################################################
# Function to pull geoglows data for prior purposes
################################################################################
pull_geoglows = function(reaches, start_date){
  source_python(paste0(inPath, 'src/geoglows_aws_pull.py'))
  model_flow = pull_tributary(reach_id = reaches, start_date=start_date)
  model_flow$Date = as.Date(row.names(model_flow))
  return(data.table(model_flow))
}
################################################################################
# LakeFlow algorithm. 
################################################################################
lakeFlow = function(lake){
  
  # Use dynamic prior Q. False = SOS prior estimate from GRADES / MAF geoglows
  use_ts_prior=TRUE
  
  # Use modeled daily tributary flows. False = mean monthly grades tributaries / MAF geoglows
  use_ts_tributary=TRUE
  
  index=which(names(viable_data)==lake)
  relevant_data = viable_data[index][[1]]
  lakeObs = relevant_data[[1]]
  upObs = relevant_data[[2]]
  dnObs = relevant_data[[3]]
  
  upID = unlist(strsplit(updated_pld$U_reach_id[updated_pld$lake_id==lake], ','))
  dnID = unlist(strsplit(updated_pld$D_reach_id[updated_pld$lake_id==lake], ','))
  
  
  # The below chunk is now ran above to speed things up. 
  # #Pull in SWOT river data and subset predownloaded SWOT lake data. 
  # upID = unlist(strsplit(updated_pld$U_reach_id[updated_pld$lake_id==lake], ','))
  # dnID = unlist(strsplit(updated_pld$D_reach_id[updated_pld$lake_id==lake], ','))
  # upObs_all = try(lapply(upID, pull_data), silent=TRUE)
  # if(is.error(upObs_all)){return(NA)}
  # upObs_all = rbindlist(lapply(upObs_all, filter_function))
  # upObs_all$time=as_datetime(upObs_all$time_str)
  # dnObs_all= try(lapply(dnID, pull_data), silent=TRUE)
  # if(is.error(dnObs_all)){return(NA)}
  # dnObs_all = rbindlist(lapply(dnObs_all, filter_function))
  # dnObs_all$time = as_datetime(dnObs_all$time_str)
  # lakeObs_all = lakeData[lakeData$lake_id==lake,]
  # lakeObs_all$time = as_datetime(lakeObs_all$time_str)
  # lakeObs_all = tukey_test_lake(lakeObs_all)
  # 
  # # FIXME: changing SWOT widths to sword widths due to large errors. 
  # #upObs_all$width=upObs_all$p_width
  # #dnObs_all$width=dnObs_all$p_width
  # 
  # 
  # # FIXME: changing lake areas to pld areas due to some weird errors. 
  # prior_area = updated_pld$Lake_area[updated_pld$lake_id==lake]
  # lakeObs_all$area_total = prior_area
  # 
  # if(nrow(lakeObs_all)<3){return(NA)}
  # 
  # ################################################################################
  # # Get dates in proper format and subset to matching dates. 
  # ################################################################################
  # lakeObs_all$date = as.Date(lakeObs_all$time)
  # upObs_all$date = as.Date(upObs_all$time)
  # dnObs_all$date = as.Date(dnObs_all$time)
  # 
  # # FIXME: Aggregating lakes to mean values for multiple observations in one day. 
  # lakeObs = data.table(lakeObs_all)[,c('wse', 'area_total', 'date')][,lapply(.SD, mean), by=date]
  # upObs = data.table(upObs_all)[,c('wse', 'width', 'slope', 'slope2','reach_id', 'date')][,lapply(.SD, mean), by=list(date, reach_id)]
  # dnObs = data.table(dnObs_all)[,c('wse', 'width', 'slope', 'slope2','reach_id', 'date')][,lapply(.SD, mean), by=list(date, reach_id)]
  # 
  # lkDates = unique(lakeObs$date)
  # upDts = upObs[,.N,by=date][N>=length(upID)] # limit to dates with obs for each upstream reach.
  # dnDts = dnObs[,.N,by=date][N>=length(dnID)] # limit to dates with obs for each downstream reach. 
  # 
  # #goodDates = lkDates[lkDates%in%upObs_all$date&lkDates%in%dnObs_all$date]
  # goodDates = lkDates[lkDates%in%upDts$date&lkDates%in%dnDts$date]
  # 
  # lakeObsGood = lakeObs[lakeObs$date%in%goodDates,]
  # upObsGood = upObs[upObs$date%in%goodDates,]
  # dnObsGood = dnObs[dnObs$date%in%goodDates,]
  # 
  # lakeObs = lakeObsGood[order(lakeObsGood$date),]
  # upObs = upObsGood[order(upObsGood$date),]
  # dnObs = dnObsGood[order(dnObsGood$date),]
  # 
  # upObs = upObs[order(upObs$reach_id),]
  # dnObs = dnObs[order(dnObs$reach_id),]
  # 
  # if(nrow(lakeObs)<3){return(NA)}
  # 
  # FIXME: use upObs$d_x_area once it is released by JPL  
  d_x_area = function(wse, width){
    wse_from_min = wse - min(wse)
    x_area = width * wse_from_min
    area_dt = c(NA, diff(x_area))
    #return(area_dt)
    return(x_area)
  }
  
  # Apply d_x_area function from above. 
  upObs$d_x_area = upObs[,d_x_area(wse, width),by=reach_id]$V1
  dnObs$d_x_area = dnObs[,d_x_area(wse, width), by=reach_id]$V1
  
  # Calculate lake storage and storage change
   lake_wse_from_min = lakeObs$wse - min(lakeObs$wse)
   area_m2 = lakeObs$area_total * 1e6 # convert from km2 to m2
   lakeObs$storage = area_m2 * lake_wse_from_min
  # lakeObs$storage_dt = c(NA, diff(lakeObs$storage))
  ht_change = c(NA, diff(lakeObs$wse))
  area_m2 = lakeObs$area_total*1e6

  area_val = list()
  for(k in 2:length(area_m2)){
    current = area_m2[k]
    prior = area_m2[k-1]
    area_add = current+prior
    area_mult= current*prior
    area_sqrt = sqrt(area_mult)
    area_val[[k]] = area_add+area_sqrt
  }
  area_param = c(NA, unlist(area_val))
  lakeObs$storage_dt = (ht_change*area_param)/3
  
  
  
  # Put data into table matching LakeFlow code (synthetic dataset): 
  lakeObsOut = lakeObs
  upObsOut = upObs
  dnObsOut = dnObs
  names(lakeObsOut) = paste0(names(lakeObs), "_l")
  names(upObsOut) = paste0(names(upObs), "_u")
  names(dnObsOut) = paste0(names(dnObs), "_d")
  
  # Split dataframes by group into lists of dataframes.  
  up_df_list = split(upObsOut, by='reach_id_u')
  dn_df_list = split(dnObsOut, by='reach_id_d')
  
  # Used to assign ancillary data - not needed at the moment. 
  lakeObsOut$month = month(lakeObsOut$date_l)
  
  # add in tributary data:Either use geoglow (ts==TRUE or use GRADES-hydroDL mean monthly vals)
  if(use_ts_tributary==TRUE){
  tributary_locations = fread(paste0(inPath, 'in/ancillary/tributaries.csv'))
  tributary_locations = tributary_locations[tributary_locations$lake_id==(lake),]
  if(nrow(tributary_locations)==0){
    lakeObsOut$tributary_total=0
    }
  else{
  tributary_reaches = unique(tributary_locations$LINKNO[tributary_locations$lake_id==lake])
  tributary_reaches = as.list(tributary_reaches)
  tributary_data = download_tributary(tributary_reaches,'01-01-1940')
  mean_annual = mean(tributary_data[,year:=lubridate::year(date)][,mean(tributary_total),year]$V1)
  lakeObsOut$tributary_total = tributary_data$tributary_total[match(lakeObsOut$date_l, tributary_data$date)]
  }
  }else{
    # Removed GRADES monthly chunk and just using 
    # tributary_locations = fread(paste0(inPath, 'in/ancillary/tributaries_merit.csv'))
    # tributary_locations = tributary_locations[tributary_locations$lake_id==(lake),]
    # tributary_data = fread(paste0(inPath, 'in/ancillary/na_monthly_tributaries_grades.csv'))
    # tributary_data = tributary_data[tributary_data$lake_id==lake,]
    # if(nrow(tributary_data)==0){lakeObsOut$tributary_total=0}
    tributary_locations = fread(paste0(inPath, 'in/ancillary/tributaries.csv'))
    tributary_locations = tributary_locations[tributary_locations$lake_id==(lake),]
    if(nrow(tributary_locations)==0){
      lakeObsOut$tributary_total=0
    }else{
    tributary_reaches = unique(tributary_locations$LINKNO[tributary_locations$lake_id==lake])
    tributary_reaches = as.list(tributary_reaches)
    tributary_data = download_tributary(tributary_reaches,'01-01-1940')
    mean_annual = mean(tributary_data[,year:=lubridate::year(date)][,mean(tributary_total),year]$V1)
    lakeObsOut$tributary_total = mean_annual
  }
}
  # add in et data: 
  et = fread(paste0(inPath, 'in/ancillary/et.csv'))
  et_lake = et[et$lake_id==lake,]
  if(nrow(et_lake)==0){
    lakeObsOut$et = 0
  }
  et_lake$month = as.numeric(et_lake$month)
  lakeObsOut$et = et_lake$mean[match(lakeObsOut$month, et_lake$month)]
  
  # n days between observations. 
  lakeObsOut$n_days = c(NA, as.numeric(diff(date(lakeObsOut$date_l))))
  
  # Remove first day from all dataframes bc dv is missing. 
  lakeObsOut = lakeObsOut[-1,]
  up_df_list = lapply(up_df_list, function(f) f[-1,])
  dn_df_list = lapply(dn_df_list, function(f) f[-1,])
  
  # Transpose the data to get in proper format for stan. 
  up_df_wide = lapply(up_df_list, t)
  dn_df_wide = lapply(dn_df_list, t)
  
  # Combine the swot obs of multiple reaches back together.  
  up_df_out = do.call('rbind', up_df_wide)
  dn_df_out = do.call('rbind', dn_df_wide)
  
  # Create seperate matrices for each type of swot obs. Time across, reach down
  up_df_stan = lapply(names(upObsOut)[-1], function(f) as.matrix(up_df_out[grep(f,rownames(up_df_out)),]))
  names(up_df_stan) = names(upObsOut)[-1]
  up_df_stan = lapply(up_df_stan, apply, 2, as.numeric)
  
  dn_df_stan = lapply(names(dnObsOut)[-1], function(f) as.matrix(dn_df_out[grep(f,rownames(dn_df_out)),]))
  names(dn_df_stan) = names(dnObsOut)[-1]
  dn_df_stan = lapply(dn_df_stan, apply, 2, as.numeric)
  
  transposer <- function(df) {
    z<-t(df)
  }
  
  # If reach is only one value, transpose the data to match multiple inflow/outflow setups. 
  if(length(upID)==1){
    up_df_stan = lapply(up_df_stan, transposer)
  }
  
  if(length(dnID)==1){
    dn_df_stan = lapply(dn_df_stan, transposer)
  }
  
  # dA shift function. 
  da_shift_fun = function(x) median(x) - min(x)
  
  # dA scale to be greater than 0.
  da_scale = function(x) x - min(x)
  
  # Apply the above functions. 
  up_da_shift = array(apply(up_df_stan$d_x_area_u,1, da_shift_fun))
  d_x_area_scaled_u = t(as.matrix(apply(up_df_stan$d_x_area_u, 1, da_scale)))
  
  dn_da_shift = array(apply(dn_df_stan$d_x_area_d,1, da_shift_fun))
  d_x_area_scaled_d = t(as.matrix(apply(dn_df_stan$d_x_area_d, 1, da_scale)))
  
  # Function to extract priors from SOS. 
  sos_pull = function(reach_id){
    sos = paste0(inPath,"in/sos/constrained/na_sword_v15_SOS_priors.nc")
    sos_outflow = RNetCDF::open.nc(sos)
    reach_grp <- RNetCDF::grp.inq.nc(sos_outflow, "reaches")$self
    reach_ids <- RNetCDF::var.get.nc(reach_grp, "reach_id")
    index <- which(reach_ids==reach_id, arr.ind=TRUE)
    gbp_grp <- RNetCDF::grp.inq.nc(sos_outflow, "gbpriors/reach")$self
    geoA = RNetCDF::var.get.nc(gbp_grp, "logA0_hat")[index]
    geoN = RNetCDF::var.get.nc(gbp_grp, "logn_hat")[index]
    geoNlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_logn")[index]
    geoNupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_logn")[index]
    geoAlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_A0")[index]
    geoAupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_A0")[index]
    geoAsd = RNetCDF::var.get.nc(gbp_grp, "logA0_sd")[index]
    geoNsd = RNetCDF::var.get.nc(gbp_grp, "logn_sd")[index]
    model_grp <- RNetCDF::grp.inq.nc(sos_outflow, "model")$self
    qHat <- RNetCDF::var.get.nc(model_grp, "mean_q")[index]
    #Assign a prior of 1000 cms is q prior is NA. - probably makes sense to use some relationship between width and meanQ in the future.   
    qHat = ifelse(is.na(qHat), 1000, qHat)
    qUpper <- RNetCDF::var.get.nc(model_grp, "max_q")[index]
    qLower <- RNetCDF::var.get.nc(model_grp, "min_q")[index]
    qLower = ifelse(qLower<=0|is.na(qLower), 0.001, qLower)
    qUpper = ifelse(is.na(qUpper),500000, qUpper)
    sigma = RNetCDF::var.get.nc(gbp_grp, "sigma_man")[index]
    qSd=RNetCDF::var.get.nc(gbp_grp, "logQ_sd")[index]
    return(data.table(reach_id, geoA, geoN, geoNlower, geoNupper, geoAlower, geoAupper,
                      geoAsd, geoNsd, qHat, qUpper, qLower, sigma, qSd))
  }

  # Create our own geobam priors when SoS provides null priors - common issue bc lake influenced reaches often don't pass SoS QAQC.
  sos_fit = function(df){
    reach_id = unique(unlist(df[grep('reach_id', names(df))]))#[[1]]
    width_matrix = df[grep('width', names(df))][[1]]
    slope_matrix = df[grep('slope2', names(df))][[1]]
    geoA = estimate_logA0(width_matrix)
    geoN = estimate_logn(slope_matrix, width_matrix)
    geoNlower = estimate_lowerboundlogn(width_matrix)
    geoNupper = estimate_upperboundlogn(width_matrix)
    geoAlower = estimate_lowerboundA0(width_matrix)
    geoAupper = estimate_upperboundA0(width_matrix)
    geoAsd = estimate_A0SD(width_matrix)
    geoNsd = estimate_lognSD(width_matrix)
    #Talk to Craig about more informed estimates of the below. 
    qHat = rep(NA, length(reach_id))
    qLower = rep(NA, length(reach_id))
    qUpper = rep(NA, length(reach_id))
    qSd = rep(100, length(reach_id))
    sigma = rep(0.25, length(reach_id))
    #Assign a prior of 1000 cms is q prior is NA. - probably makes sense to use some relationship between width and meanQ in the future.   
    qHat = ifelse(is.na(qHat), 1000, qHat)
    qLower = ifelse(qLower<=0|is.na(qLower), 0.001, qLower)
    qUpper = ifelse(is.na(qUpper),500000, qUpper)
    return(list(data.table(reach_id, geoA, geoN, geoNlower, geoNupper, geoAlower, geoAupper,
                      geoAsd, geoNsd, qHat, qUpper, qLower, sigma, qSd)))
  }
  
  
  
  # Add u and d labels to each column of sos outputs - helps with organizing the data. 
  nms_paste = function(df, added_label){
    nms = colnames(df)
    new_nms = paste0(nms, '_', added_label)
    colnames(df) = new_nms
    return(df)
  }
  
  # Pull priors from sos
  up_sos = lapply(upID, sos_pull)
  sos_geobam = 'sos'
  if(any(unlist(lapply(up_sos, nrow))==0)){
    up_sos = sos_fit(up_df_stan)
    sos_geobam = 'geobam'
  }
  up_sos = lapply(up_sos, nms_paste, 'u')
  
  dn_sos = lapply(dnID, sos_pull)
  if(any(unlist(lapply(dn_sos, nrow))==0)){
    dn_sos = sos_fit(dn_df_stan)
    sos_geobam = 'geobam'
  }
  dn_sos = lapply(dn_sos, nms_paste, 'd')
  
  # Transpose the sos data to get in proper format for stan. 
  up_sos_wide = lapply(up_sos, t)
  dn_sos_wide = lapply(dn_sos, t)
  
  # Combine the sos priors of multiple reaches together.  
  up_sos_out = do.call('rbind', up_sos_wide)
  dn_sos_out = do.call('rbind', dn_sos_wide)
  
  # Create seperate matrices for each type of sos prior. reach down time across. 
  up_sos_stan = lapply(colnames(up_sos[[1]]), function(f) as.matrix(up_sos_out[grep(f,rownames(up_sos_out)),]))
  names(up_sos_stan) = colnames(up_sos[[1]])
  up_sos_stan = lapply(up_sos_stan, apply, 2, as.numeric)
  
  dn_sos_stan = lapply(colnames(dn_sos[[1]]), function(f) as.matrix(dn_sos_out[grep(f,rownames(dn_sos_out)),]))
  names(dn_sos_stan) = colnames(dn_sos[[1]])
  dn_sos_stan = lapply(dn_sos_stan, apply, 2, as.numeric)
  
  # Place in array to meet stan needs. 
  up_sos_stan = lapply(up_sos_stan, array)
  dn_sos_stan = lapply(dn_sos_stan, array)
  
  # Create matrix for select priors that need to be length of obs for stan. 
  convert_to_matrix = function(f){
    n = length(f)
    output = apply(f, 1, rep, nrow(lakeObsOut))
    return(t(output))
  }
  
  up_sos_stan$sigma_u=convert_to_matrix(up_sos_stan$sigma_u)
  dn_sos_stan$sigma_d=convert_to_matrix(dn_sos_stan$sigma_d)
  up_sos_stan$qSd_u = convert_to_matrix(up_sos_stan$qSd_u)
  dn_sos_stan$qSd_d = convert_to_matrix(dn_sos_stan$qSd_d)
  up_sos_stan$qHat_u = convert_to_matrix(up_sos_stan$qHat_u)
  dn_sos_stan$qHat_d = convert_to_matrix(dn_sos_stan$qHat_d)
  
  if(length(up_sos_stan$reach_id_u)==0){return(NA)}
  if(length(dn_sos_stan$reach_id_d)==0){return(NA)}
  
  
  # Pull in Modeled timeseries as prior Q rather than mean: #Updating alternative to using static Geoglows mean annual flow. 
  if(use_ts_prior==TRUE){
    sword_geoglows = fread(paste0(inPath, '/in/ancillary/sword_geoglows.csv'))
    sword_reaches = c(upID, dnID)
    sword_geoglows_filt = sword_geoglows[sword_geoglows$reach_id%in%sword_reaches,c('reach_id','LINKNO')]
    geoglows_reaches = unique(as.list(sword_geoglows$LINKNO[sword_geoglows$reach_id%in%sword_reaches]))
    # Pull in modeled geoglows data. 
    model_data = pull_geoglows(geoglows_reaches, '01-01-1940')
    # Convert bad Q to NA. 
    ind = ncol(model_data)-1
    model_data[,1:ind][model_data[,1:ind] < 0] <- 0.1
    # When model data is unavailable (NRT SWOT data), use mean model. 
    missing_dates = data.table(Date=seq.Date((max(model_data$Date)+1), Sys.Date(), by=1))
    model_data = bind_rows(model_data, missing_dates)
    model_data[] <- lapply(model_data, function(x) { 
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      x
    })
    model_data_all = model_data
    model_data = model_data[model_data$Date%in%lakeObsOut$date_l,]
    model_wide = melt(model_data, id.vars=c('Date'))
    model_wide$LINKNO = as.character(model_wide$variable)
    model_wide$reach_id = sword_geoglows_filt$reach_id[match(model_wide$LINKNO, sword_geoglows_filt$LINKNO)]
    model_list = split(model_wide, by='reach_id')
    
    up_ind = lapply(up_sos_stan$reach_id_u,function(x){which(as.character(x)==names(model_list))})
    dn_ind = lapply(dn_sos_stan$reach_id_d,function(x){which(as.character(x)==names(model_list))})
    
    if(length(unique(model_wide$LINKNO))==1){
      up_ind[[1]] = 1
      dn_ind[[1]] = 1
    }
    
    if(length(up_ind[[1]])==0){
      up_ind[[1]] = 2
    }
    
    up_qhat = lapply(unlist(up_ind), function(x){t(as.matrix(model_list[[x]]$value))})
    up_qhat = do.call(rbind, up_qhat)
    
    dn_qhat = lapply(unlist(dn_ind), function(x){t(as.matrix(model_list[[x]]$value))})
    dn_qhat = do.call(rbind, dn_qhat)
    
    up_sos_stan$qHat_u = up_qhat
    dn_sos_stan$qHat_d = dn_qhat
    
    # Replace upper and lower limit from SoS with Geoglows. Removed in case the geoglows estimate is wrong. LF can overcome bad prior but not bad low/up bounds.
    # model_upper_limit = apply(model_data_all, 2, max)
    # model_lower_limit = apply(model_data_all, 2, min)
    # model_limits = data.table(bind_rows(model_upper_limit, model_lower_limit))
    # model_limits$type=c('upper', 'lower')
    # model_limits = melt(model_limits,id.vars=c('type'))
    # model_limits = model_limits[model_limits$variable!='Date',]
    # model_limits$LINKNO = as.character(model_limits$variable)
    # model_limits$reach_id = sword_geoglows_filt$reach_id[match(model_limits$LINKNO, sword_geoglows_filt$LINKNO)]
    # model_limits$value = as.numeric(model_limits$value)
    # model_limits = split(model_limits, by='reach_id')
    # 
    # up_sos_stan$qUpper_u=as.array(unlist(lapply(unlist(up_ind), function(x){t(as.matrix(model_limits[[x]][type=='upper']$value))})))
    # up_sos_stan$qLower_u=as.array(unlist(lapply(unlist(up_ind), function(x){t(as.matrix(model_limits[[x]][type=='lower']$value))})))
    # dn_sos_stan$qUpper_d=as.array(unlist(lapply(unlist(dn_ind), function(x){t(as.matrix(model_limits[[x]][type=='upper']$value))})))
    # dn_sos_stan$qLower_d=as.array(unlist(lapply(unlist(dn_ind), function(x){t(as.matrix(model_limits[[x]][type=='lower']$value))})))
  }else{
    sword_geoglows = fread(paste0(inPath, '/in/ancillary/sword_geoglows.csv'))
    sword_reaches = c(upID, dnID)
    sword_geoglows_filt = sword_geoglows[sword_geoglows$reach_id%in%sword_reaches,c('reach_id','LINKNO')]
    geoglows_reaches = unique(as.list(sword_geoglows$LINKNO[sword_geoglows$reach_id%in%sword_reaches]))
    # Pull in modeled geoglows data. 
    model_data = pull_geoglows(geoglows_reaches, '01-01-1940')
    # Convert bad Q to NA. 
    ind = ncol(model_data)-1
    model_data[,1:ind][model_data[,1:ind] < 0] <- 0.1
    # When model data is unavailable (NRT SWOT data), use mean model. 
    missing_dates = data.table(Date=seq.Date((max(model_data$Date)+1), Sys.Date(), by=1))
    model_data = bind_rows(model_data, missing_dates)
    model_data[] <- lapply(model_data, function(x) { 
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      x
    })
    #model_data = model_data[model_data$Date%in%lakeObsOut$date_l,]
    model_annual = model_data[,year:=lubridate::year(Date),]
    model_wide = melt(model_data, id.vars=c('year'))
    model_wide = model_wide[,list(value=mean(value)),by=list(year, variable)][,list(value=mean(value)),by=variable][variable!='Date']
    model_wide = model_wide[rep(model_wide[,.I], nrow(lakeObsOut))]
    
    model_wide$LINKNO = as.character(model_wide$variable)
    model_wide$reach_id = sword_geoglows_filt$reach_id[match(model_wide$LINKNO, sword_geoglows_filt$LINKNO)]
    model_list = split(model_wide, by='reach_id')
    
    up_ind = lapply(up_sos_stan$reach_id_u,function(x){which(as.character(x)==names(model_list))})
    dn_ind = lapply(dn_sos_stan$reach_id_d,function(x){which(as.character(x)==names(model_list))})
    
    if(length(unique(model_wide$LINKNO))==1){
      up_ind[[1]] = 1
      dn_ind[[1]] = 1
    }
    
    up_qhat = lapply(unlist(up_ind), function(x){t(as.matrix(model_list[[x]]$value))})
    up_qhat = do.call(rbind, up_qhat)
    
    dn_qhat = lapply(unlist(dn_ind), function(x){t(as.matrix(model_list[[x]]$value))})
    dn_qhat = do.call(rbind, dn_qhat)
    
    up_sos_stan$qHat_u = up_qhat
    dn_sos_stan$qHat_d = dn_qhat
  }
  
  
  
  # Combine all data needed for stan. 
  stan_data = list(N = nrow(lakeObsOut),
                   n1=length(upID),
                   n2=length(dnID),
                   sigmaIn = up_sos_stan$sigma_u,
                   sigmaOut = dn_sos_stan$sigma_d,
                   qInSd = up_sos_stan$qSd_u,#rep(qInSd,nrow(data)),#qInSd used to multiply by 0.1?
                   qOutSd = dn_sos_stan$qSd_d,#rep(qOutSd,nrow(data)),#qOutSd used to multiply by 0.1?
                   q = log(up_sos_stan$qHat_u),#rep(log(qHatIn),nrow(data)),
                   sigma = up_sos_stan$sigma_u,#rep(0.25, nrow(data)),
                   da = d_x_area_scaled_u,#data$upArea_dt_u-min(data$upArea_dt_u, na.rm=T), # possible error from Ryan: dA is already change in A, but he's subtracting min dA 
                   w=log(up_df_stan$width_u),
                   s=log(up_df_stan$slope2_u), # using smoothed slope (slope2) to minimize negatives
                   da2=d_x_area_scaled_d,#da2=data$dnArea_dt_d-min(data$dnArea_dt_d, na.rm=T), # possible error from Ryan: dA is already change in A, but he's subtracting min dA 
                   w2=log(dn_df_stan$width_d),
                   s2=log(dn_df_stan$slope2_d), # = smoothed slope (slope2) produced similar results as raw slope 
                   q2=log(dn_sos_stan$qHat_d),#rep(log(qHatOut),nrow(data)),
                   #dv=data$storage_dt_l/86400, #seconds per day
                   et=lakeObsOut$et,#rep(0, nrow(lakeObsOut)), # filled with zero vals right now
                   lateral=lakeObsOut$tributary_total,#rep(0,nrow(lakeObsOut)), # filled with zero vals right now
                   dv_per = (lakeObsOut$storage_dt_l/86400)/lakeObsOut$n_days)
  
  
  # set NAs, NaNs, and Inf to 0 for stan (Ryan's code):
  stan_data$da = ifelse(is.na(stan_data$da)|is.infinite(stan_data$da), 0, stan_data$da)
  stan_data$da2 = ifelse(is.na(stan_data$da2)|is.infinite(stan_data$da2), 0, stan_data$da2)
  stan_data$s = ifelse(is.na(stan_data$s), 0, stan_data$s)
  stan_data$s2 = ifelse(is.na(stan_data$s2), 0, stan_data$s2)
  stan_data$lateral = ifelse(is.na(stan_data$lateral), 0, stan_data$lateral)
  stan_data$et = ifelse(is.na(stan_data$et), 0, stan_data$et)
  # stan_data$inc_none=none
  # stan_data$inc_et=et
  # stan_data$inc_lateral=q
  # stan_data$inc_et_lateral=et_q
  stan_data$nInlower = up_sos_stan$geoNlower_u #as.array(geoNInlower)
  stan_data$nInupper = up_sos_stan$geoNupper_u#as.array(geoNInupper)
  stan_data$aInlower= up_sos_stan$geoAlower_u#as.array(geoAInlower)
  stan_data$aInupper = up_sos_stan$geoAupper_u#as.array(geoAInupper)
  stan_data$nOutlower = dn_sos_stan$geoNlower_d#as.array(geoNOutlower)
  stan_data$nOutupper = dn_sos_stan$geoNupper_d#as.array(geoNOutupper)
  stan_data$aOutlower= dn_sos_stan$geoAlower_d#as.array(geoAOutlower)
  stan_data$aOutupper = dn_sos_stan$geoAupper_d#as.array(geoAOutupper)
  stan_data$daInShift = up_da_shift#apply(up_df_stan$d_x_area_u,1, da_shift_fun)#as.array(da_shift_fun(data$upArea_dt_u))
  stan_data$daOutShift = dn_da_shift#apply(dn_df_stan$d_x_area_d,1, da_shift_fun)#as.array(da_shift_fun(data$dnArea_dt_d))
  stan_data$nInHat = up_sos_stan$geoN_u#as.array(geoNin)
  stan_data$nInSd = up_sos_stan$geoNsd_u#as.array(0.25) # Ryan, what's this? 
  stan_data$aInHat = up_sos_stan$geoA_u#as.array(geoAin)
  stan_data$aInSd = up_sos_stan$geoAsd_u#as.array(geoAinSD)
  stan_data$qInupper=log(up_sos_stan$qUpper_u)#as.array(log(qInupper))
  stan_data$qInlower=log(up_sos_stan$qLower_u)#as.array(log(qInlower))
  stan_data$nOutHat = dn_sos_stan$geoN_d#as.array(geoNout) 
  stan_data$nOutSd = dn_sos_stan$geoNsd_d#as.array(0.25)
  stan_data$aOutHat = dn_sos_stan$geoA_d#as.array(geoAout)
  stan_data$aOutSd = dn_sos_stan$geoAsd_d#as.array(geoAoutSD)
  stan_data$qOutupper=log(dn_sos_stan$qUpper_d)#as.array(log(qOutupper))
  stan_data$qOutlower=log(dn_sos_stan$qLower_d)#as.array(log(qOutlower))
  
  # Apply stan code. 
  fit = try(stan(paste0(inPath, "src/lakeflow_stan_flexible.stan"),
                 data=stan_data,
                 chains=6,#4, #3
                 cores=6, #6
                 iter=4000,#4000, #iter=4000
                 control=list(stepsize=0.5,
                              adapt_delta=0.9)))
  if(is.error(fit)){next}
  
  # Manning's eqn:
  eqn1 = function(n, a, da, w, s){
    flow = (n^-1)*((a+da)^(5/3))*(w^(-2/3))*(s^0.5)
    return(flow)
  }
  
  # Pull outputs from stan. 
  mn = get_posterior_mean(fit)
  rw = row.names(mn)
  mn = data.table(mn)
  mn$type = rw
  if(nrow(mn)==0|ncol(mn)<6){return(NA)}
  mn = mn[, c('type','mean-all chains')]
  
  uncertainty = data.table(summary(fit)$summary)
  uncertainty$type = mn$type
  
  # Manning's n estimates from stan
  roughness_inflow = mn[startsWith(mn$type,'n[')]
  roughness_inflow_sd = uncertainty[startsWith(uncertainty$type,'n[')]
  roughness_outflow = mn[startsWith(mn$type,'nOut[')]
  roughness_outflow_sd = uncertainty[startsWith(uncertainty$type,'nOut[')]
  
  # A0 estimates from stan
  bath_inflow = mn[startsWith(mn$type,'a[')]
  bath_inflow_sd = uncertainty[startsWith(uncertainty$type,'a[')]
  bath_outflow = mn[startsWith(mn$type,'aOut[')]
  bath_outflow_sd = uncertainty[startsWith(uncertainty$type,'aOut[')]
  
  
  # Bayesian estimate of inflow and outflow - not used for our purposes. 
  bayes_inflow = mn[startsWith(mn$type,'logQ_in[')]
  bayes_inflow_sd = uncertainty[startsWith(uncertainty$type,'logQ_in[')]$sd
  bayes_outflow = mn[startsWith(mn$type,'logQ_out[')]
  bayes_outflow_sd = uncertainty[startsWith(uncertainty$type,'logQ_out[')]$sd
  
  # Quick way to organize the data in case of varying N of inflows/outflows - could be placed in a function at some point. 
  inflow_outputs = list()
  for(j in 1:length(upID)){
    q_estimate = eqn1(exp(roughness_inflow$`mean-all chains`[j]),bath_inflow$`mean-all chains`[j],
                      d_x_area_scaled_u[j,],up_df_stan$width_u[j,], up_df_stan$slope2_u[j,])
    q_bayes = bayes_inflow$`mean-all chains`[bayes_inflow$type]
    inflow_outputs[[j]] = data.table(q_lakeflow = q_estimate, reach_id=upID[j], n_lakeflow=exp(roughness_inflow$`mean-all chains`[j]),a0_lakeflow=bath_inflow$`mean-all chains`[j],date=lakeObsOut$date_l,lake_id=lake,q_model=up_sos_stan$qHat_u[j,],
                                     width = stan_data$w[j,], slope2 = stan_data$s[j,], da = d_x_area_scaled_u[j,], wse = up_df_stan$wse_u[j,],
                                     storage=lakeObsOut$storage_l,dv=stan_data$dv_per, tributary=stan_data$lateral, et=stan_data$et,type='inflow',n_lakeflow_sd = roughness_inflow_sd$sd[j],a0_lakeflow_sd =bath_inflow_sd$sd[j])
  }
  inflow_outputs = rbindlist(inflow_outputs)
  inflow_outputs$bayes_q = bayes_inflow$`mean-all chains`
  inflow_outputs$bayes_q_sd = bayes_inflow_sd
  
  outflow_outputs = list()
  for(j in 1:length(dnID)){
    q_estimate = eqn1(exp(roughness_outflow$`mean-all chains`[j]),bath_outflow$`mean-all chains`[j],
                      d_x_area_scaled_d[j,],dn_df_stan$width_d[j,], dn_df_stan$slope2_d[j,])
    outflow_outputs[[j]] = data.table(q_lakeflow = q_estimate, reach_id=dnID[j], n_lakeflow=exp(roughness_outflow$`mean-all chains`[j]),a0_lakeflow=bath_outflow$`mean-all chains`[j],date=lakeObsOut$date_l,lake_id=lake,q_model=dn_sos_stan$qHat_d[j,],
                                      width = stan_data$w2[j,], slope2 = stan_data$s2[j,], da = d_x_area_scaled_d[j,], wse = dn_df_stan$wse_d[j,],
                                      storage=lakeObsOut$storage_l,dv=stan_data$dv_per, tributary=stan_data$lateral, et=stan_data$et,type='outflow',n_lakeflow_sd = roughness_outflow_sd$sd[j],a0_lakeflow_sd =bath_outflow_sd$sd[j])
  }
  outflow_outputs = rbindlist(outflow_outputs)
  outflow_outputs$bayes_q = bayes_outflow$`mean-all chains`
  outflow_outputs$bayes_q_sd = bayes_outflow_sd
  
  
  
  output_df = bind_rows(inflow_outputs, outflow_outputs)
  output_df$prior_fit = sos_geobam
  fwrite(output_df, paste0(inPath, '/out/lf_results_na_3/', lake, '.csv'))
  return(output_df)
}

#Ignore: This is just for subsetting to lakes that haven't run yet if the code breaks during a run. 
# successful = list.files('C:\\Users\\rriggs\\OneDrive - DOI\\Research\\LakeFlow_local\\out\\lf_results_na_3\\')
# successful = gsub('.csv', '', successful)
# missing = viable_locations[viable_locations$lake%!in%successful,]
# viable_locations = missing

#Apply the LakeFlow code. - Could also be done using lapply but I like to see it print where I'm at.
output_list = list()
for(i in 1:nrow(viable_locations)){
  print(i)
  output_list[[i]] = lakeFlow(viable_locations$lake[i])
}
lf_outputs = rbindlist(output_list[!is.na(output_list)])
lf_outputs = lf_outputs[!is.na(lf_outputs$q_lakeflow),]