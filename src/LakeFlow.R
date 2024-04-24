##LakeFlow code for reading input data into Confluence. 
##By: Ryan Riggs and George Allen, February 2024.

args = commandArgs(trailingOnly = TRUE)

LakeFlow = function(inPath, outPath, lake_id){
  library(foreign)
  library(lubridate)
  library(rstan)
  library(data.table)
  library(dplyr)
  library(sf)
  library(raster)
  library(rstudioapi)
  #library(rgdal)
################################################################################
# Read in harmonized river-lake database: 
################################################################################
#pld = paste0(inPath, "in/SWORDv15_PLDv103_step2_harmonized.gdb")
#fc_list = ogrListLayers(pld)
# Read the feature class
#lakes = st_read(dsn=pld,layer="PLD_lakes_intersect")
lakes = fread(paste0(inPath, 'in/SWORDv15_PLDv103_step2_harmonized_table.csv'))
################################################################################
# Read in the SWOT data from processed csv files:
################################################################################ 

# read in processed LakeSP and RiverSP observations:

# Lakes: 
pldObs_inPath = paste0(inPath, "in/SWOT_FSP_SP_processed/SWOT_FSP_SP_processed_lakeObs.csv")
pldObs = read.csv(pldObs_inPath, header=T)
# Rivers: 
prdObs_inPath = paste0(inPath, "in/SWOT_FSP_SP_processed/SWOT_FSP_SP_processed_riverObs.csv")
prdObs = read.csv(prdObs_inPath, header=T)
# Ancillary - full of zeros right now: 
ancObs = fread(paste0(inPath, 'in/ancillary/ancillary_data.csv'))

################################################################################
# Process data for a given lake
################################################################################

# run code on Connecticut River Reservoir 
lakeID = lake_id
lakesFilt = lakes[lakes$lake_id==lakeID,]
upID = lakesFilt$us_rch_id 
dnID = lakesFilt$ds_rch_id 
ancData = ancObs[ancObs$lake_id==lakeID,]

# FIXME
# Manually assign nearby reaches due to issues with FSP data. 
#upID = 73120000271
#dnID = 73120000191 # note: 73120000241 is the reach directly downstream but it is missing from the SoS so I used the segment slightly lower

# process SWOT data: 

lakeObs_all = pldObs[pldObs$lake_id==lakeID,] 
upObs_all = prdObs[prdObs$reach_id==upID,]
dnObs_all = prdObs[prdObs$reach_id==dnID,]

# dates with observations of lake and all reaches:
lakeTime = as_datetime(lakeObs_all$time, origin="2000-01-01")
upTime = as_datetime(upObs_all$time, origin="2000-01-01")
dnTime = as_datetime(dnObs_all$time, origin="2000-01-01")
lakeHour = as.character(round_date(lakeTime, "hour"))
upHour = as.character(round_date(upTime, "hour"))
dnHour = as.character(round_date(dnTime, "hour"))

lakeUpMatch = upHour[match(lakeHour, upHour)]
datesWithGoodData = dnHour[match(lakeUpMatch, dnHour)]
datesWithGoodData = datesWithGoodData[!is.na(datesWithGoodData)]

lakeHour_goodInd = match(datesWithGoodData, lakeHour)
upHour_goodInd = match(datesWithGoodData, upHour)
dnHour_goodInd = match(datesWithGoodData, dnHour)

lakeObs = lakeObs_all[lakeHour_goodInd, ]
upObs = upObs_all[upHour_goodInd, ]
dnObs = dnObs_all[dnHour_goodInd, ]

# times should match within a few seconds: 
as_datetime(lakeObs$time, origin="2000-01-01") - as_datetime(upObs$time, origin="2000-01-01")
as_datetime(lakeObs$time, origin="2000-01-01") - as_datetime(dnObs$time, origin="2000-01-01")

# calculate change in lake volume: (Gao et al., 2012)
# FIXME: once the official SWOT parameter "p_storage" is released by JPL, 
# use lakeObs$p_storage once it instead of calculating it manually below 
lake_wse_from_min = lakeObs$wse - min(lakeObs$wse)
area_m2 = lakeObs$area_total * 1e6 # convert from km2 to m2
lakeObs$storage = area_m2 * lake_wse_from_min

# calculate change in river cross-sectional area(Durand et al.m 2014 JoH)
# FIXME: use upObs$d_x_area once it is released by JPL
up_wse_from_min = upObs$wse - min(upObs$wse)
upObs$x_area = upObs$width * up_wse_from_min

dn_wse_from_min = dnObs$wse - min(dnObs$wse)
dnObs$x_area = dnObs$width * dn_wse_from_min

# calculate dt change in volume & area between time steps:
lakeObs$storage_dt = c(NA, diff(lakeObs$storage))
upObs$upArea_dt = c(NA, diff(upObs$x_area))
dnObs$dnArea_dt = c(NA, diff(dnObs$x_area))

# put data into table matching LakeFlow code (synthetic dataset): 
lakeObsOut = lakeObs
upObsOut = upObs
dnObsOut = dnObs
names(lakeObsOut) = paste0(names(lakeObs), "_l")
names(upObsOut) = paste0(names(upObs), "_u")
names(dnObsOut) = paste0(names(dnObs), "_d")

data = cbind(lakeObsOut, upObsOut, dnObsOut)
data = cbind(lakeObsOut, upObsOut, dnObsOut)[-1,] # remove first row bc it has an NA

# add in the ancillary data based on month matchup:
data$month = month(data$time_str_d)
data = left_join(data, ancData, by='month')

# n days between observations. 
# FIXME double check when we have inconsistent samples. 
data$n_days = as.numeric(diff(date(datesWithGoodData)))

################################################################################
# Run LakeFlow snippet
################################################################################


sos = paste0(inPath,"in/sos/constrained/na_sword_v15_SOS_priors.nc")
sos_outflow = RNetCDF::open.nc(sos)
reach_grp <- RNetCDF::grp.inq.nc(sos_outflow, "reaches")$self
reach_ids <- RNetCDF::var.get.nc(reach_grp, "reach_id")
index <- which(reach_ids==upID, arr.ind=TRUE)
gbp_grp <- RNetCDF::grp.inq.nc(sos_outflow, "gbpriors/reach")$self
geoAin = RNetCDF::var.get.nc(gbp_grp, "logA0_hat")[index]
geoNin = RNetCDF::var.get.nc(gbp_grp, "logn_hat")[index]
geoNInlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_logn")[index]
geoNInupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_logn")[index]
geoAInlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_A0")[index]
geoAInupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_A0")[index]
geoAinSD = RNetCDF::var.get.nc(gbp_grp, "logA0_sd")[index]
geoNinSD = RNetCDF::var.get.nc(gbp_grp, "logn_sd")[index]
model_grp <- RNetCDF::grp.inq.nc(sos_outflow, "model")$self
qHatIn <- RNetCDF::var.get.nc(model_grp, "mean_q")[index]
sigmaIn = RNetCDF::var.get.nc(gbp_grp, "sigma_man")[index]
qInSd=RNetCDF::var.get.nc(gbp_grp, "logQ_sd")[index]
##outflow
index <- which(reach_ids==dnID, arr.ind=TRUE)
gbp_grp <- RNetCDF::grp.inq.nc(sos_outflow, "gbpriors/reach")$self
geoAout = RNetCDF::var.get.nc(gbp_grp, "logA0_hat")[index]
geoNout = RNetCDF::var.get.nc(gbp_grp, "logn_hat")[index]
geoNOutlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_logn")[index]
geoNOutupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_logn")[index]
geoAOutlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_A0")[index]
geoAOutupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_A0")[index]
geoAoutSD = RNetCDF::var.get.nc(gbp_grp, "logA0_sd")[index]
geoNoutSD = RNetCDF::var.get.nc(gbp_grp, "logn_sd")[index]
model_grp <- RNetCDF::grp.inq.nc(sos_outflow, "model")$self
qHatOut <- RNetCDF::var.get.nc(model_grp, "mean_q")[index]
qHatIn = ifelse(is.na(qHatIn), qHatOut, qHatIn)
qHatOut = ifelse(is.na(qHatOut), qHatIn, qHatOut)
sigmaOut = RNetCDF::var.get.nc(gbp_grp, "sigma_man")[index]
qOutSd=RNetCDF::var.get.nc(gbp_grp, "logQ_sd")[index]
# qHatIn = ifelse(is.na(qHatIn), mean(data$`in Q (m3/s)`), qHatIn)
# qHatOut = ifelse(is.na(qHatOut), mean(data$`out Q (m3/s)`), qHatOut)

da_shift_fun = function(x) median(x) - min(x)

stan_data = list(N = nrow(data),
                 sigmaIn = rep(sigmaIn, nrow(data)),
                 sigmaOut = rep(sigmaOut, nrow(data)),
                 qInSd = rep(log(qHatIn*.1),nrow(data)),#qInSd
                 qOutSd = rep(log(qHatOut*.1),nrow(data)),#qOutSd
                 q = rep(log(qHatIn),nrow(data)),
                 sigma = rep(0.25, nrow(data)),
                 da = data$upArea_dt_u,#-min(data$upArea_dt_u, na.rm=T), # possible error from Ryan: dA is already change in A, but he's subtracting min dA 
                 w=log(data$width_u),
                 s=log(data$slope2_u), # using smoothed slope (slope2) to minimize negatives
                 da2=data$dnArea_dt_d,#-min(data$dnArea_dt_d, na.rm=T), # possible error from Ryan: dA is already change in A, but he's subtracting min dA 
                 w2=log(data$width_d),
                 s2=log(data$slope2_d), # = smoothed slope (slope2) produced similar results as raw slope 
                 q2=rep(log(qHatOut),nrow(data)),
                 #dv=data$storage_dt_l/86400, #seconds per day
                 et=data$et, # filled with zero vals right now
                 lateral=data$lateral_q, # filled with zero vals right now
                 dv_per = (data$storage_dt_l/86400)/data$n_days)

# set NAs, NaNs, and Inf to 0 for stan (Ryan's code):
stan_data$da = ifelse(is.na(stan_data$da)|is.infinite(stan_data$da), 0, stan_data$da)
stan_data$da2 = ifelse(is.na(stan_data$da2)|is.infinite(stan_data$da), 0, stan_data$da2)
stan_data$s = ifelse(is.na(stan_data$s), 0, stan_data$s)
stan_data$s2 = ifelse(is.na(stan_data$s2), 0, stan_data$s2)
stan_data$lateral = ifelse(is.na(stan_data$lateral), 0, stan_data$lateral)
stan_data$et = ifelse(is.na(stan_data$et), 0, stan_data$et)
# stan_data$inc_none=none
# stan_data$inc_et=et
# stan_data$inc_lateral=q
# stan_data$inc_et_lateral=et_q
stan_data$nInlower = geoNInlower
stan_data$nInupper = geoNInupper
stan_data$aInlower= geoAInlower
stan_data$aInupper = geoAInupper
stan_data$nOutlower = geoNInlower
stan_data$nOutupper = geoNInupper
stan_data$aOutlower= geoAOutlower
stan_data$aOutupper = geoAOutupper
stan_data$daInShift = da_shift_fun(data$upArea_dt_u)
stan_data$daOutShift = da_shift_fun(data$dnArea_dt_d)
stan_data$nInHat = geoNin 
stan_data$nInSd = 0.25 # Ryan, what's this? 
stan_data$aInHat = geoAin
stan_data$aInSd = geoAinSD
stan_data$nOutHat = geoNout 
stan_data$nOutSd = 0.25
stan_data$aOutHat = geoAout
stan_data$aOutSd = geoAoutSD

fit = stan(paste0(inPath, "src/lakeflow_ifelse.stan"),
           data=stan_data,
           chains=3, #3
           cores=9, #6
           iter=4000, #iter=4000
           control=list(stepsize=0.5,
                        adapt_delta=0.9))

# Manning's eqn:
eqn1 = function(n, a, da, w, s){
  flow = (n^-1)*((a+da)^(5/3))*(w^(-2/3))*(s^0.5)
  return(flow)
}


mn = get_posterior_mean(fit)
rw = row.names(mn)
mn = data.table(mn)
mn$type = rw
mn = mn[, c('type','mean-all chains')]

model = eqn1(exp(mn$`mean-all chains`[mn$type=='n']), 
             (mn$`mean-all chains`[mn$type=='a']), 
             data$upArea_dt_u , 
             data$width_u, 
             data$slope2_u)
modelOut = eqn1(exp(mn$`mean-all chains`[mn$type=='nOut']), 
                mn$`mean-all chains`[mn$type=='aOut'], 
                data$dnArea_dt_d, 
                data$width_d, 
                data$slope2_d)  

# questions for Ryan: 
# 1) where is time variable included in stan? irregular observations 
# 2) is the SoS the correct version? Is there a shapefile? 
# 3) issue with subtracting dV from itself 

nObs=length(model)
n=exp(mn$`mean-all chains`[mn$type=='n'])
a=(mn$`mean-all chains`[mn$type=='a'])
no=exp(mn$`mean-all chains`[mn$type=='nOut'])
ao=(mn$`mean-all chains`[mn$type=='aOut'])

##output values
q=c(model,modelOut)
id=c(rep(upID, nObs),rep(dnID, nObs))
n=c(rep(n, nObs),rep(no, nObs))
a=c(rep(a, nObs),rep(ao, nObs))
date=rep(as_datetime(data$time_l, origin="2000-01-01"), 2)
output_df=data.table(q_lakeflow=q, reach_id=id, n_lakeflow=n, a0_lakeflow=a, date=date)
fwrite(output_df, paste0(outPath, 'test.csv'))
print(output_df)
}

LakeFlow(args[1], args[2], args[3])

