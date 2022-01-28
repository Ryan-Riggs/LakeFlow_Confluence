################################################################################################################################
#Lakeflow inputs.  
#Ryan Riggs
#1/28/2022
################################################################################################################################

getInputdata = function(lake_id, lake_reach_file,data_dir, sos, et, lateralQ){
  #read in lake dataset. 
  lake_file = rgdal::readOGR(lake_reach_file)@data
  lake = lake_file[lake_file$lake_id==lake_id]
  reach_in = lake$reach_id_in
  reach_out = lake$reach_id_out
  
  ##open .nc and extract obs. for lake. 
  lakeFile = paste0(data_dir, lake, "_SWOT.nc")
  lake_input = RNetCDF::open.nc(lakeFile)
  volume_change = RNetCDF::var.get.nc(lake_input, "delta_s_q")
  time_lake = as.Date(RNetCDF::var.get.nc(lake_input, "time_str"))
  RNetCDF::close.nc(lake_input)
  
  ##open .nc and extract obs. for inflow. 
  swot = paste0(data_dir, reach_in, "_SWOT.nc")
  swot_input = RNetCDF::open.nc(swot)
  nt_inflow = RNetCDF::var.get.nc(swot_input, "nt")
  reach_grp = RNetCDF::grp.inq.nc(swot_input, "reach")$self
  width_inflow = RNetCDF::var.get.nc(reach_grp, "width")
  dArea_inflow = RNetCDF::var.get.nc(reach_grp, "d_x_area")
  slope_inflow = RNetCDF::var.get.nc(reach_grp, "slope2")
  time_inflow = as.Date(RNetCDF::var.get.nc(reach_grp, "time_str"))
  sos_inflow = RNetCDF::open.nc(sos)
  reach_grp <- RNetCDF::grp.inq.nc(sos_inflow, "reaches")$self
  reach_ids <- RNetCDF::var.get.nc(reach_grp, "reach_id")
  index <- which(reach_ids==reach_in, arr.ind=TRUE)
  gbp_grp <- RNetCDF::grp.inq.nc(sos_inflow, "gbpriors/reach")$self
  db <- exp(RNetCDF::var.get.nc(gbp_grp, "logDb_hat")[index])
  dA_lw_inflow = RNetCDF::var.get.nc(gbp_grp, "lowerbound_A0")[index]
  dA_up_inflow = RNetCDF::var.get.nc(gbp_grp, "upperbound_A0")[index]
  n_lw_inflow = exp(RNetCDF::var.get.nc(gbp_grp, "lowerbound_logn")[index])
  n_up_inflow = exp(RNetCDF::var.get.nc(gbp_grp, "upperbound_logn")[index])
  
  #close .nc file. 
  RNetCDF::close.nc(swot_input)
  RNetCDF::close.nc(sos_inflow)
  
  ##open .nc and extract obs. for outflow. 
  swot = paste0(data_dir, reach_out, "_SWOT.nc")
  swot_input = RNetCDF::open.nc(swot)
  nt_outflow = RNetCDF::var.get.nc(swot_input, "nt")
  reach_grp = RNetCDF::grp.inq.nc(swot_input, "reach")$self
  width_outflow = RNetCDF::var.get.nc(reach_grp, "width")
  dArea_outflow = RNetCDF::var.get.nc(reach_grp, "d_x_area")
  slope_outflow = RNetCDF::var.get.nc(reach_grp, "slope2")
  time_outflow = as.Date(RNetCDF::var.get.nc(reach_grp, "time_str"))
  sos_outflow = RNetCDF::open.nc(sos)
  reach_grp <- RNetCDF::grp.inq.nc(sos_outflow, "reaches")$self
  reach_ids <- RNetCDF::var.get.nc(reach_grp, "reach_id")
  index <- which(reach_ids==reach_out, arr.ind=TRUE)
  gbp_grp <- RNetCDF::grp.inq.nc(sos_inflow, "gbpriors/reach")$self
  db <- exp(RNetCDF::var.get.nc(gbp_grp, "logDb_hat")[index])
  dA_lw_outflow = RNetCDF::var.get.nc(gbp_grp, "lowerbound_A0")[index]
  dA_up_outflow = RNetCDF::var.get.nc(gbp_grp, "upperbound_A0")[index]
  n_lw_outflow = exp(RNetCDF::var.get.nc(gbp_grp, "lowerbound_logn")[index])
  n_up_outflow = exp(RNetCDF::var.get.nc(gbp_grp, "upperbound_logn")[index])
  
  #close .nc file. 
  RNetCDF::close.nc(swot_input)
  RNetCDF::close.nc(sos_outflow)
  
  #apply exponents. 
  win = width_inflow^(-2/3)
  sin = slope_inflow^0.5
  wout = width_outflow^(-2/3)
  sout = slope_outflow^0.5
  
  #read in additional data. 
  et = data.table::fread(et)
  etDf = et[et$lake_id==lake_id,]
  etDf= reshape2::melt(etDf)
  lateralQ = data.table::fread(lateralQ)
  lateralQDf = lateralQ[lateralQ$lake_id==lake_id,]
  lateralQDf = reshape2::melt(etDf)
  
  ##Combine inflow, outflow, and lake and filter to same dates. 
  inflow = data.frame(win, sin, d_x_area_in = dArea_inflow, Date=time_inflow)
  outflow = data.frame(wout, sout, d_x_area_out=dArea_outflow, Date=time_outflow)
  lake = data.frame(delta_s_q = volume_change, Date=time_lake)
  all = merge(inflow, outflow, by="Date")
  all = merge(all, lake, by="Date")
  all$month = lubridate::month(all$Date)
  all$et = etDf$value[match(all$month, etDf$variable)]
  all$sum = lateralQDf$value[match(all$month, lateralQDf$variable)]
  
  data = all
  priors = data.frame(dA_lw_inflow,dA_lw_outflow,dA_up_inflow,dA_up_outflow,n_lw_inflow,n_lw_outflow,n_up_inflow,n_up_outflow)
  ids = data.frame(lake_id, reach_in, reach_out)
  export = list(data, priors, ids)
  return(export)
}






