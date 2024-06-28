##LakeFlow code for running locally
##By: Ryan Riggs and George Allen, June 2024.

# *TODO: Test to see if all the DAWG flags help improve performance. 
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
################################################################################
# Read in relevant files: Harmonized sword-pld, reservoirs of interest, swot lake data
################################################################################
updated_pld = fread(paste0(inPath,"in/SWORDv16_PLDv103.csv"))
updated_pld$lake_id =  as.character(updated_pld$lake_id)
reservoirs = data.table::fread(paste0(inPath,'in/', 'filtered_res_61024_static_widths.csv'))
lakeData = fread(paste0(inPath, '/in/swot_lake_data.csv'))

# Below is for me so that I can easily add new swot observations to the folder. 
#swot_path = list.files(paste0(inPath, 'in/data_downloads/'), pattern='.dbf', full.names=TRUE)
swot_path = list.files('C:/Users/rriggs/OneDrive - DOI/Research/SWOT/src/data_downloads', pattern='.dbf', full.names=TRUE)


open_and_filter = function(f){
  file = foreign::read.dbf(f)
  file_sub = file#[file$lake_id%in%reservoirs$lake_id,]
  return(file_sub)
}

files_filt = lapply(swot_path, open_and_filter)
combined = data.table::rbindlist(files_filt)

##Filter to non flagged lake obs.
lakeData = combined[combined$quality_f=='0'&combined$partial_f=='0',]
lakeData = lakeData%>%distinct(.keep_all=TRUE)
lakeData$lake_id = as.character(lakeData$lake_id)

#FIXME remove lakedata with multiple ids, what's causing that??
lakeData$lake_id_first = sub(";.*", "", lakeData$lake_id)
lakeData$lake_id=lakeData$lake_id_first

################################################################################
# Function to pull reach data. 
################################################################################
pull_data = function(feature_id){
  website = paste0('https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries?feature=Reach&feature_id=',feature_id, '&start_time=2023-01-01T00:00:00Z&end_time=2024-12-31T00:00:00Z&output=csv&fields=reach_id,time_str,wse,width,slope,slope2,d_x_area,area_total,reach_q,p_width,xovr_cal_q,partial_f,dark_frac,ice_clim_f')
  response = GET(website)
  pull = content(response, as='parsed')$results
  data = read.csv(textConnection(pull$csv), sep=',')
  data$reach_id = feature_id
  return(data)
}################################################################################
# Function to filter SWOT reach data. 
################################################################################
filter_function = function(swot_ts){
  dawg_filter = swot_ts[swot_ts$time!='no_data'&swot_ts$ice_clim_f<2&swot_ts$dark_frac<=0.5&swot_ts$xovr_cal_q<2,]
  qual_filter = swot_ts[swot_ts$time!='no_data'&swot_ts$reach_q<=2,]
  return(qual_filter)
}
################################################################################
# Function to pull tributary inflow Q estimates. 
################################################################################
start_date = '20230101'

download_tributary = function(reaches){
  source_python(paste0(inPath, 'src/geoglows_aws_pull.py'))
  tributary_flow = pull_tributary(reach_id = reaches)
  tributary_flow$date = as.Date(row.names(tributary_flow))
  n_col = ncol(tributary_flow)
  tributary_aggregated = data.table(tributary_flow)[,tributary_total:=rowSums(.SD), .SDcols=-n_col]
  return(tributary_aggregated)
}
################################################################################
# LakeFlow algorithm. 
################################################################################
lakeFlow = function(lake){
  #Pull in SWOT river data and subset predownloaded SWOT lake data. 
  upID = unlist(strsplit(updated_pld$U_reach_id[updated_pld$lake_id==lake], ','))
  dnID = unlist(strsplit(updated_pld$D_reach_id[updated_pld$lake_id==lake], ','))
  upObs_all = rbindlist(lapply(upID, pull_data))
  upObs_all = filter_function(upObs_all)
  #upObs_all = upObs_all[upObs_all$time!='no_data'&upObs_all$reach_q<=2,]
  upObs_all$time=as_datetime(upObs_all$time_str)
  dnObs_all= rbindlist(lapply(dnID, pull_data))
  #dnObs_all = dnObs_all[dnObs_all$time!='no_data'&dnObs_all$reach_q<=2,]
  dnObs_all = filter_function(dnObs_all)
  dnObs_all$time = as_datetime(dnObs_all$time_str)
  lakeObs_all = lakeData[lakeData$lake_id==lake,]
  lakeObs_all$time = as_datetime(lakeObs_all$time_str)
  
  # FIXME: changing SWOT widths to sword widths due to large errors. 
  upObs_all$width=upObs_all$p_width
  dnObs_all$width=dnObs_all$p_width
  
  
  # FIXME: changing lake areas to pld areas due to some weird errors. 
  prior_area = updated_pld$Lake_area[updated_pld$lake_id==lake]
  lakeObs_all$area_total = prior_area
  
  if(nrow(lakeObs_all)<3){return(NA)}
  
  ################################################################################
  # Get dates in proper format and subset to matching dates. 
  ################################################################################
  lakeObs_all$date = as.Date(lakeObs_all$time)
  upObs_all$date = as.Date(upObs_all$time)
  dnObs_all$date = as.Date(dnObs_all$time)
  
  lkDates = unique(lakeObs_all$date)
  upDts = upObs_all[,.N,by=date][N>=length(upID)] # limit to dates with obs for each upstream reach.
  dnDts = dnObs_all[,.N,by=date][N>=length(dnID)] # limit to dates with obs for each downstream reach. 
  
  #goodDates = lkDates[lkDates%in%upObs_all$date&lkDates%in%dnObs_all$date]
  goodDates = lkDates[lkDates%in%upDts$date&lkDates%in%dnDts$date]
  
  lakeObsGood = lakeObs_all[lakeObs_all$date%in%goodDates,]
  upObsGood = upObs_all[upObs_all$date%in%goodDates,]
  dnObsGood = dnObs_all[dnObs_all$date%in%goodDates,]
  
  lakeObs = lakeObsGood[order(lakeObsGood$date),]
  upObs = upObsGood[order(upObsGood$date),]
  dnObs = dnObsGood[order(dnObsGood$date),]
  
  # FIXME: Aggregating lakes to mean values for multiple observations in one day. 
  lakeObs = data.table(lakeObs)[,c('wse', 'area_total', 'date')][,lapply(.SD, mean), by=date]
  upObs = data.table(upObs)[,c('wse', 'width', 'slope', 'slope2','reach_id', 'date')][,lapply(.SD, mean), by=list(date, reach_id)]
  dnObs = data.table(dnObs)[,c('wse', 'width', 'slope', 'slope2','reach_id', 'date')][,lapply(.SD, mean), by=list(date, reach_id)]
  
  upObs = upObs[order(upObs$reach_id),]
  dnObs = dnObs[order(dnObs$reach_id),]
  
  if(nrow(lakeObs)<3){return(NA)}
  
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
  lakeObs$storage_dt = c(NA, diff(lakeObs$storage))
  
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
  
  # add in tributary data:
  tributary_locations = fread(paste0(inPath, 'in/ancillary/tributaries.csv'))
  tributary_locations = tributary_locations[tributary_locations$lake_id==(lake),]
  if(nrow(tributary_locations)==0){
    lakeObsOut$tributary_total=0
    }
  else{
  tributary_reaches = unique(tributary_locations$LINKNO[tributary_locations$lake_id==lake])
  tributary_reaches = as.list(tributary_reaches)
  tributary_data = download_tributary(tributary_reaches)
  lakeObsOut$tributary_total = tributary_data$tributary_total[match(lakeObsOut$date_l, tributary_data$date)]
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
  
  # Add u and d labels to each column of sos outputs - helps with organizing the data. 
  nms_paste = function(df, added_label){
    nms = colnames(df)
    new_nms = paste0(nms, '_', added_label)
    colnames(df) = new_nms
    return(df)
  }
  
  # Pull priors from sos
  up_sos = lapply(upID, sos_pull)
  up_sos = lapply(up_sos, nms_paste, 'u')
  
  dn_sos = lapply(dnID, sos_pull)
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
  
  # Manning's n estimates from stan
  roughness_inflow = mn[startsWith(mn$type,'n[')]
  roughness_outflow = mn[startsWith(mn$type,'nOut[')]
  
  # A0 estimates from stan
  bath_inflow = mn[startsWith(mn$type,'a[')]
  bath_outflow = mn[startsWith(mn$type,'aOut[')]
  
  # Bayesian estimate of inflow and outflow - not used for our purposes. 
  bayes_inflow = mn[startsWith(mn$type,'logQ_in[')]
  bayes_outflow = mn[startsWith(mn$type,'logQ_out[')]
  
  # Quick way to organize the data in case of varying N of inflows/outflows - could be placed in a function at some point. 
  inflow_outputs = list()
  for(j in 1:length(upID)){
    q_estimate = eqn1(exp(roughness_inflow$`mean-all chains`[j]),bath_inflow$`mean-all chains`[j],
                      d_x_area_scaled_u[j,],up_df_stan$width_u[j,], up_df_stan$slope2_u[j,])
    q_bayes = bayes_inflow$`mean-all chains`[bayes_inflow$type]
    inflow_outputs[[j]] = data.table(q_lakeflow = q_estimate, reach_id=upID[j], n_lakeflow=exp(roughness_inflow$`mean-all chains`[j]),a0_lakeflow=bath_inflow$`mean-all chains`[j],date=lakeObsOut$date_l,lake_id=lake,q_model=up_sos_stan$qHat_u[j,],
                                     width = stan_data$w[j,], slope2 = stan_data$s[j,], da = d_x_area_scaled_u[j,], wse = up_df_stan$wse_u[j,],
                                     storage=lakeObsOut$storage_l,dv=stan_data$dv_per, tributary=stan_data$lateral, et=stan_data$et,type='inflow')
  }
  inflow_outputs = rbindlist(inflow_outputs)
  inflow_outputs$bayes_q = bayes_inflow$`mean-all chains`
  
  outflow_outputs = list()
  for(j in 1:length(dnID)){
    q_estimate = eqn1(exp(roughness_outflow$`mean-all chains`[j]),bath_outflow$`mean-all chains`[j],
                      d_x_area_scaled_d[j,],dn_df_stan$width_d[j,], dn_df_stan$slope2_d[j,])
    outflow_outputs[[j]] = data.table(q_lakeflow = q_estimate, reach_id=dnID[j], n_lakeflow=exp(roughness_outflow$`mean-all chains`[j]),a0_lakeflow=bath_outflow$`mean-all chains`[j],date=lakeObsOut$date_l,lake_id=lake,q_model=dn_sos_stan$qHat_d[j,],
                                      width = stan_data$w2[j,], slope2 = stan_data$s2[j,], da = d_x_area_scaled_d[j,], wse = dn_df_stan$wse_d[j,],
                                      storage=lakeObsOut$storage_l,dv=stan_data$dv_per, tributary=stan_data$lateral, et=stan_data$et,type='outflow')
  }
  outflow_outputs = rbindlist(outflow_outputs)
  outflow_outputs$bayes_q = bayes_outflow$`mean-all chains`
  
  output_df = bind_rows(inflow_outputs, outflow_outputs)
  fwrite(output_df, paste0(inPath, '/out/lf_results/', lake, '.csv'))
  return(output_df)
}


# Filter to lakes with at least n observations. 
n = 5
lakes = unique(data.table(lakeData)[,.N,by=lake_id][N>=n]$lake_id)
lakes = lakes#[lakes%in%reservoirs$lake_id]

# Apply LakeFlow at the first three lakes. 
#John Redmond '7420130653'
#Fontenelle '7720023433'
#Flaming Gorge '7720025003'
#Mohave '7720003433'
#Lake Chesdin '7310006793'
#Allatoona: '7320350863'

id_subset = '7420130653'
index=which(lakes==id_subset)
lf_results = lapply(lakes[index], lakeFlow)
lf_outputs = rbindlist(lf_results[!is.na(lf_results)])

# lf_results = lapply(lakes, lakeFlow)
# output_files = list.files(paste0(inPath, '/out/lf_results/'), full.names=TRUE)
# lf_outputs = rbindlist(lapply(output_files, fread))
# lf_outputs$Date = as.Date(lf_outputs$Date)
# lf_outputs$reach_id = as.character(lf_outputs$reach_id)
# lf_outputs$lake_id=as.character(lf_outputs$lake_id)

#fwrite(lf_outputs,paste0(inPath, '/out/lf_results.csv'))
################################################################################
# Validate lakeflow observations with USGS data.  
################################################################################
# Auto assign gages to reaches.  
library(dataRetrieval)
pull_geom = function(feature, feature_id){
  website = paste0('https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries?feature=', feature, '&feature_id=',feature_id, '&start_time=2023-01-01T00:00:00Z&end_time=2024-12-31T00:00:00Z&output=csv&fields=reach_id,geometry')
  response = GET(website)
  pull = content(response, as='parsed')$results
  data = read.csv(textConnection(pull$csv), sep=',')
  #data$time_str = as.Date(data$time_str)
  return(data)
}

geom_function = function(f){
  geom=pull_geom('Reach', f)[1,]
  line = st_as_sf(geom, wkt='geometry')
  return(line)
}
reaches = as.character(unique(lf_outputs$reach_id))
reach_geometries_list = lapply(reaches, geom_function)
reach_geometries = rbindlist(reach_geometries_list)
reach_geometries = st_as_sf(reach_geometries)
st_crs(reach_geometries) = 4326

# Gages in relevant states where there is a lakeflow reach. 
library(tigris)
states=states()
states_filt = st_join(reach_geometries%>%st_transform(crs(states)), states, st_within)
st = na.omit(unique(states_filt$STUSPS))

# Find gages in relevant states that have discharge data. 
state_list = list()
for(i in 1:length(st)){
  stations = whatNWISsites(stateCd=st[i],parameterCd='00060',hasDataTypeCd='iv')
  state_list[[i]] = stations
}
usgs = rbindlist(state_list)

# Link to gages. - distance is in meters. 
usgs = st_as_sf(usgs,coords=c('dec_long_va', 'dec_lat_va'))
st_crs(usgs) = 4326
reaches = unique(lf_outputs$reach_id)

nearest = st_nearest_feature(usgs,reach_geometries)
distance = st_distance(usgs, reach_geometries[nearest,],by_element = TRUE)
usgs$reach_id = reach_geometries$reach_id[nearest]
usgs$distance = as.numeric(distance)

# limit to usgs gages within 1km or whatever distance. 
usgs_filt = usgs[usgs$distance<=5000,]

# Remove gages with 'CREEK', etc. in the name. 
smallName = grep('CREEK|Creek|creek|STREAM|Stream|stream|INLET|Inlet|inlet',usgs_filt$station_nm)
usgs_filt=usgs_filt[-smallName,]
usgs_sites = unique(usgs_filt$site_no)

pull_gage = function(f){
  data = try(RivRetrieve::usa(site=f, variable='discharge'))
  if(is.error(data)){return(NA)}
  data$gage=f
  data=data[data$Date>=as.Date('2022-12-31'),]
  return(data)
}
gage_df_list = lapply(usgs_sites, pull_gage)
gage_df = rbindlist(gage_df_list[!is.na(gage_df_list)])
gage_df$reach_id = usgs_filt$reach_id[match(gage_df$gage, usgs_filt$site_no)]
gage_df$reach_id=as.character(gage_df$reach_id)

lf_outputs$Date = as.Date(lf_outputs$date)
validation = left_join(gage_df, lf_outputs,by=c('Date', 'reach_id'))
#fwrite(validation,paste0(inPath, 'out/gage_lf.csv'))


s=ggplot(validation)+
  geom_point(aes(x=Q,y=q_lakeflow,color=reach_id))+
  geom_point(aes(x=Q,y=exp(bayes_q)),col='black',size=0.5)+
  geom_abline(aes(slope=1,intercept=0))+
  facet_wrap(~gage)+
  scale_x_log10()+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = 'top')
s

g=ggplot(validation)+
  geom_line(aes(x=Date, y=Q))+
  geom_point(aes(x=Date,y=q_lakeflow,col=reach_id))+
  #geom_point(aes(x=Date,y=exp(bayes_q)), col='black')+
  geom_hline(aes(yintercept=q_model,col=reach_id))+
  facet_wrap(~gage, scales='free_y')+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = 'top')
g

#ggsave(paste0(inPath, 'figures/early_validation_v3.png'),g,bg='white', width=12,height=6,dpi=1000,units='in')
#ggsave(paste0(inPath,'figures/early_validation_scatter_v3.png'),s,bg='white', width=12,height=6,dpi=1000,units='in')

gages_of_interest = c('07182510')
voi = validation[validation$gage%in%gages_of_interest,]

ggplot(voi[!is.na(q_lakeflow),])+
  geom_line(aes(x=Date,y=Q))+
  geom_point(aes(x=Date,y=q_lakeflow),col='orange')+
  geom_point(aes(x=Date,y=exp(bayes_q)), col='green')+
  facet_wrap(~gage, scales='free')

reach_data = pull_data('74244700461')
reach_data = reach_data[reach_data$time!='no_data'&reach_data$reach_q<=2,]
reach_data$Date = as_datetime(reach_data$time_str)
reach_data$Date = as.Date(reach_data$Date)

stage_data = RivRetrieve::usa('07182510', 'stage')
stage_comp = merge(stage_data, reach_data, 'Date')
stage_comp$diff_h = (stage_comp$H-min(stage_comp$H))/(max(stage_comp$H)-min(stage_comp$H))
stage_comp$diff_swot = (stage_comp$wse-min(stage_comp$wse))/(max(stage_comp$wse)-min(stage_comp$wse))

ggplot(stage_comp)+
  geom_point(aes(x=diff_h, y= diff_swot))+
  geom_abline(aes(intercept=0,slope=1))
  


##John redmond lake, KS
lf_outputs$Date = as.Date(lf_outputs$date)
gage = '07182390' #inflow
gageData=RivRetrieve::usa(gage, 'discharge')
gageData$reach_id = '74244700461'
val_inflow = right_join(lf_outputs, gageData[gageData$Date>as.Date('2023-01-01'),], by=c('Date','reach_id'))
val_inflow$type='Inflow'

gage = '07182510' #outflow
gageData=RivRetrieve::usa(gage, 'discharge')
gageData$reach_id = '74244700431'
val_outflow = right_join(lf_outputs, gageData[gageData$Date>as.Date('2023-01-01'),], by=c('Date','reach_id'))
val_outflow$type='Outflow'
val = bind_rows(val_inflow, val_outflow)

in_out_plot = ggplot(val)+
  geom_line(aes(x=Date,y=Q))+
  geom_point(aes(x=Date, y=q_lakeflow), col='red')+
  geom_point(aes(x=Date, y= exp(bayes_q)), col='green')+
  geom_hline(aes(yintercept=q_model),col='red')+
  scale_y_log10()+#(limits=c(1,1000))+
  ylab('Discharge (cms)')+
  xlab('')+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5), 
        strip.text = element_text(size=14))+
  facet_wrap(~type)
in_out_plot

ggsave(paste0(inPath, 'out/figures/John_Redmond_Lake_static_widths.png'),in_out_plot,bg='white', width=8,height=6,dpi=1000,units='in')

##Error metrics. 
source('C:\\Users\\rriggs\\OneDrive - DOI\\Research\\Error_stats\\src\\error_functions.R')

error_stats=data.table(val)[type=='Outflow',validate(q_lakeflow, Q),by=type]

error_stats=data.table(val)[,validate(exp(bayes_q), Q),by=type]


################################################################################
##Fontenelle
lf_outputs$Date = as.Date(lf_outputs$date)
gage = '09209400' #inflow
gageData=RivRetrieve::usa(gage, 'discharge')
gageData$reach_id = '77299000171'
val_inflow = right_join(lf_outputs, gageData[gageData$Date>as.Date('2023-01-01'),], by=c('Date','reach_id'))
val_inflow$type='Inflow'

gage = '09211200' #outflow
gageData=RivRetrieve::usa(gage, 'discharge')
gageData$reach_id = '77299000131'
val_outflow = right_join(lf_outputs, gageData[gageData$Date>as.Date('2023-01-01'),], by=c('Date','reach_id'))
val_outflow$type='Outflow'
val = bind_rows(val_inflow, val_outflow)

ggplot(val)+
  geom_line(aes(x=Date,y=Q))+
  geom_point(aes(x=Date, y=q_lakeflow), col='red')+
  geom_point(aes(x=Date, y= exp(bayes_q)), col='green')+
  geom_hline(aes(yintercept=q_model),col='red')+
  scale_y_log10()+#(limits=c(1,1000))+
  ylab('Discharge (cms)')+
  xlab('')+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5), 
        strip.text = element_text(size=14))+
  facet_wrap(~type, scales='free')

ggplot(val)+
  geom_point(aes(x=Q, y=q_lakeflow), col='red')+
  geom_point(aes(x=Q, y= exp(bayes_q)), col='green')+
  ylab('Discharge (cms)')+
  xlab('')+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5), 
        strip.text = element_text(size=14))+
  facet_wrap(~type, scales='free')+
  geom_smooth(aes(x=Q,y=q_lakeflow), method=lm)

#ggsave(paste0(inPath, 'out/figures/Fontenelle_static_widths.png'),in_out_plot,bg='white', width=8,height=6,dpi=1000,units='in')
################################################################################
##Allatoona
lf_outputs$Date = as.Date(lf_outputs$date)
gage = '02392000' #inflow
gageData=RivRetrieve::usa(gage, 'discharge')
gageData$reach_id = '73289000241'
val_inflow = right_join(lf_outputs, gageData[gageData$Date>as.Date('2023-01-01'),], by=c('Date','reach_id'))
val_inflow$type='Inflow'

gage = '02394000' #outflow
gageData=RivRetrieve::usa(gage, 'discharge')
gageData$reach_id = '73289000194'
val_outflow = right_join(lf_outputs, gageData[gageData$Date>as.Date('2023-01-01'),], by=c('Date','reach_id'))
val_outflow$type='Outflow'
val = bind_rows(val_inflow, val_outflow)

ggplot(val)+
  geom_line(aes(x=Date,y=Q))+
  geom_point(aes(x=Date, y=q_lakeflow), col='red')+
  geom_point(aes(x=Date, y= exp(bayes_q)), col='green')+
  geom_hline(aes(yintercept=q_model),col='red')+
  scale_y_log10()+#(limits=c(1,1000))+
  ylab('Discharge (cms)')+
  xlab('')+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5), 
        strip.text = element_text(size=14))+
  facet_wrap(~type, scales='free')

ggplot(val)+
  geom_point(aes(x=Q, y=q_lakeflow), col='red')+
  geom_point(aes(x=Q, y= exp(bayes_q)), col='green')+
  ylab('Discharge (cms)')+
  xlab('')+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5), 
        strip.text = element_text(size=14))+
  facet_wrap(~type, scales='free')+
  geom_smooth(aes(x=Q,y=q_lakeflow), method=lm)

#ggsave(paste0(inPath, 'out/figures/Fontenelle_static_widths.png'),in_out_plot,bg='white', width=8,height=6,dpi=1000,units='in')


################################################################################
################################################################################
##Chesdin
lf_outputs$Date = as.Date(lf_outputs$date)
gage = '02040892' #inflow
gageData=RivRetrieve::usa(gage, 'discharge')
gageData$reach_id = '73190800051'
val_inflow = right_join(lf_outputs, gageData[gageData$Date>as.Date('2023-01-01'),], by=c('Date','reach_id'))
val_inflow$type='Inflow'

gage = '02041650' #outflow
gageData=RivRetrieve::usa(gage, 'discharge')
gageData$reach_id = '73190800031'
val_outflow = right_join(lf_outputs, gageData[gageData$Date>as.Date('2023-01-01'),], by=c('Date','reach_id'))
val_outflow$type='Outflow'
val = bind_rows(val_inflow, val_outflow)

ggplot(val)+
  geom_line(aes(x=Date,y=Q))+
  geom_point(aes(x=Date, y=q_lakeflow), col='red')+
  #geom_point(aes(x=Date, y= exp(bayes_q)), col='green')+
  geom_hline(aes(yintercept=q_model),col='red')+
  scale_y_log10()+#(limits=c(1,1000))+
  ylab('Discharge (cms)')+
  xlab('')+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5), 
        strip.text = element_text(size=14))+
  facet_wrap(~type, scales='free')

ggplot(val)+
  geom_point(aes(x=Q, y=q_lakeflow), col='red')+
  geom_point(aes(x=Q, y= exp(bayes_q)), col='green')+
  ylab('Discharge (cms)')+
  xlab('')+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5), 
        strip.text = element_text(size=14))+
  facet_wrap(~type, scales='free')+
  geom_smooth(aes(x=Q,y=q_lakeflow), method=lm)

#ggsave(paste0(inPath, 'out/figures/Fontenelle_static_widths.png'),in_out_plot,bg='white', width=8,height=6,dpi=1000,units='in')


################################################################################
################################################################################
##Barnett Reservoir
lf_outputs$Date = as.Date(lf_outputs$date)
gage = '02483500' #inflow
gageData=RivRetrieve::usa(gage, 'discharge')
gageData$reach_id = '74100600401'
val_inflow = right_join(lf_outputs, gageData[gageData$Date>as.Date('2024-01-01'),], by=c('Date','reach_id'))
val_inflow$type='Inflow'

gage = '02485735' #outflow
gageData=RivRetrieve::usa(gage, 'discharge')
gageData$reach_id = '74100600341'
val_outflow = right_join(lf_outputs, gageData[gageData$Date>as.Date('2024-01-01'),], by=c('Date','reach_id'))
val_outflow$type='Outflow'
val = bind_rows(val_inflow, val_outflow)

ggplot(val)+
  geom_line(aes(x=Date,y=Q))+
  geom_point(aes(x=Date, y=q_lakeflow), col='red')+
  #geom_point(aes(x=Date, y= exp(bayes_q)), col='green')+
  geom_hline(aes(yintercept=q_model),col='red')+
  scale_y_log10()+#(limits=c(1,1000))+
  ylab('Discharge (cms)')+
  xlab('')+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5), 
        strip.text = element_text(size=14))+
  facet_wrap(~type, scales='free')

ggplot(val)+
  geom_point(aes(x=Q, y=q_lakeflow), col='red')+
  geom_point(aes(x=Q, y= exp(bayes_q)), col='green')+
  ylab('Discharge (cms)')+
  xlab('')+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5), 
        strip.text = element_text(size=14))+
  facet_wrap(~type, scales='free')

#ggsave(paste0(inPath, 'out/figures/Fontenelle_static_widths.png'),in_out_plot,bg='white', width=8,height=6,dpi=1000,units='in')
