################################################################################################################################
#Run LakeFlow.  
#Ryan Riggs
#1/28/2022
#########################################################################################################
##Reading in, runnning LakeFlow, and exporting data. 
#########################################################################################################
source("Path/to/LakeFlow/input_data.R")
source("Path/to/LakeFlow/LakeFlow_function.R")
source("Path/to/LakeFlow/output_data.R")

########################################################################################################
##Inputs. 
########################################################################################################

##sos = "Path/to/unconstrained/0000/na_sword_v11_SOS.nc"
##data_dir = "Path/to/sample/netcdf/"
##lake_id = "7720003433"
##lake_reach_file = "Path/to/Lake_Database_sample.json"
##et = "Path/to/et.csv"
##lateralQ = "Path/to/lateralQ.csv"
##path = "output_path"

get_data = function(lake_id,lake_reach_file, data_dir,sos,et, lateralQ, output_dir){
  
##Get input data and priors from geobam. 
input = getInputdata(lake_id, lake_reach_file,data_dir,sos,et, lateralQ)
obs = input[[1]]
priors = input[[2]]
ids = input[[3]]

##Run Algorithm. 
data = LakeFlow(obs$win, obs$sin, obs$d_x_area_in, obs$wout, obs$sout,
                obs$d_x_area_out,priors$dA_lw_inflow, priors$dA_up_inflow,
                priors$dA_lw_outflow, priors$dA_up_outflow,priors$n_lw_inflow,
                priors$n_up_inflow, priors$n_lw_outflow,priors$n_up_outflow,
                obs$et,obs$sum,obs$delta_s_q,obs$Date) 
data$inflow = ids$reach_in
data$outflow = ids$reach_out

##Write output. 
write_nc(output_dir, as.numeric(ids$lake_id), as.numeric(data$inflow), as.numeric(data$outflow), data$n_inflow,
         data$A0_inflow, data$n_outflow, data$A0_outflow, data$n_inflow_sigma, data$A0_inflow_sigma,
         data$n_outflow_sigma, data$A0_outflow_sigma, data$Q_inflow, data$Q_outflow, as.numeric(data$Date))
}


