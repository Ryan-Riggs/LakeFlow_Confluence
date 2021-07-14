

##################################################################################
#RODEO for Confluence. 
#Ryan Riggs
#7/14/2021
##################################################################################
##Set working directory. 
##################################################################################
library(ncdf4)
wd = "E:/RODEO/Confluence"
setwd(wd)
##################################################################################
##Read in example rating curves. 
##################################################################################
##widthDf: width component of rating curve lookup table. 
#Each width column corresponds to the percentile value for width at that node.
#Each row represents a different node. 
#Final Column: SWOT 'node_id'
widthDf = read.csv(paste0(wd, "/inputs/widthDf.csv"))
##dischargeDf: discharge component of rating curve lookup table. 
#Each discharge column corresponds to the percentile value for discharge at that node.
#Each row represents a different node. 
#Final Column: SWOT 'node_id'
dischargeDf = read.csv(paste0(wd, "/inputs/dischargeDf.csv"))
##################################################################################
##Read in synthetic SWOT data. 
##################################################################################
##width dataframe:
#row values: SWOT 'width'. 
#column names: SWOT 'time'. 
#Final column: SWOT 'node_id'.
width = read.csv(paste0(wd,"/inputs/ExampleSWOTdata.csv"))
##################################################################################
##Process data. 
##################################################################################
##Estimate discharge for a given width.
#For each 'node_id' and 'time', look up matching 'node_id' in the width/discharge lookup tables. 
#Use the lookup tables to approximate discharge for the observed SWOT 'width'
##################################################################################
nodes = unique(width$node_id)
estimateQ = function(id){
  Wrow = match(width$node_id[width$node_id==id], widthDf$node_id)
  Qrow = match(width$node_id[width$node_id==id], dischargeDf$node_id)
  rc = approxfun(widthDf[Wrow,], dischargeDf[Qrow,])
  estimatedQ = rc(width[width$node_id==id,1:ncol(width)])
  return(estimatedQ)
}
q = lapply(nodes, estimateQ)
##################################################################################
##Place data in the netcdf file and write it to outputs folder. 
##################################################################################
id_def = ncdim_def("node_id", "id", nodes)
time_def = ncdim_def("time", "time", 1:100)
dimCross <- ncdim_def(name='q', units='cms', longname='discharge', vals=unlist(q) )
dimTime <- ncdim_def('time', units='seconds', longname='measurement year', calendar="standard", vals=c(1:(ncol(width)-1)))
dimID <- ncdim_def(name='node_id', units='ID', nodes)
varQ <- ncvar_def(name='Q', units='cms', dim=list(dimID, dimTime), missval=-9999)
vars <- list(varQ)
con <- nc_create(paste0(wd, "/outputs/output.nc"), vars)
ncatt_put(con, 'time', 'standard_name', 'time')
ncatt_put(con, 'node_id', 'axis', 'ID')
ncvar_put(con, varQ, unlist(q))
nc_close(con)
nc = nc_open(paste0(wd, "/outputs/output.nc"))
Q = ncvar_get(nc, "Q")
