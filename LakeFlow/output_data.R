################################################################################################################################
#output function.  
#Ryan Riggs
#1/28/2022
#########################################################################################################
##Write out as .nc file. 
#########################################################################################################
write_nc = function(output_dir,lake_id, reach_in, reach_out, n_inflow, A0_inflow, n_outflow, A0_outflow,
                    n_inflow_sigma, A0_inflow_sigma, n_outflow_sigma, A0_outflow_sigma, Q_inflow, Q_outflow, time){
  
id_def = ncdf4::ncdim_def("reach_id", "ID", c(reach_in,reach_out))
dimNi <- ncdf4::ncdim_def(name='n',units = "s/[m1/3]", vals=c(n_inflow,n_outflow))
dimAi <- ncdf4::ncdim_def(name='A0',units = "m2", vals=c(A0_inflow, A0_outflow))
dimNiSig <- ncdf4::ncdim_def(name='n_sigma',units = "s/[m1/3]", vals=c(n_inflow_sigma, n_outflow_sigma))
dimAiSig <- ncdf4::ncdim_def(name='A0_sigma',units = "m2", vals=c(A0_inflow_sigma, A0_outflow_sigma))
dimQin <- ncdf4::ncdim_def(name='Q',units = "m3/s", vals=c(Q_inflow, Q_outflow))
dimTime <- ncdf4::ncdim_def(name='Date',units = "Date", vals=c(time, time))


dimID <- ncdf4::ncdim_def(name='reach_id', units='ID', c(reach_in,reach_out))
varNi <- ncdf4::ncvar_def(name='n',units = "s/[m1/3]", dim=list(dimID), missval=-9999)
varAi <- ncdf4::ncvar_def(name='A0',units = "m2", dim=list(dimID), missval=-9999)
varNiSig <- ncdf4::ncvar_def(name='n_sigma',units = "s/[m1/3]", dim=list(dimID), missval=-9999)
varAiSig <- ncdf4::ncvar_def(name='A0_sigma',units = "m2", dim=list(dimID), missval=-9999)
varQin <- ncdf4::ncvar_def(name='Q',units = "m3/s", dim=list(dimID), missval=-9999)
varTime <- ncdf4::ncvar_def(name='Date',units = "Date", dim=list(dimID), missval=-9999)



vars <- list(varNi, varAi, varNiSig, varAiSig, varQin, varTime)
con <- ncdf4::nc_create(paste0(output_dir, lake_id,".nc"), vars)
ncdf4::ncatt_put(con, 'reach_id', 'axis', 'ID')
ncdf4::ncvar_put(con, varNi, c(n_inflow, n_outflow))
ncdf4::ncvar_put(con, varAi, c(A0_inflow, A0_outflow))
ncdf4::ncvar_put(con, varNiSig, c(n_inflow_sigma, n_outflow_sigma))
ncdf4::ncvar_put(con, varAiSig, c(A0_inflow_sigma, A0_outflow_sigma))
ncdf4::ncvar_put(con, varQin, c(Q_inflow, Q_outflow))
ncdf4::ncvar_put(con, varTime, c(time, time))

ncdf4::nc_close(con)
}






