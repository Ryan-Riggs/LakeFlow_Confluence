##################################################################################
#RODEO for Confluence. 
#Ryan Riggs
#7/14/2021
##################################################################################
##Set working directory. 
##Within this directory, there should be an 'inputs' folder with the .csv input files. 
##There should be an empty 'outputs' folder, and a 'src' folder with this script. 
##################################################################################
wd = "E:/Lakeflow/Confluence"
setwd(wd)
lake_id = 77233000001#77233000000
###############################################################################################################################################################################
##Read in tables and synthetic datasets. 
###################################################################################
##Lookup table for reaches inflowing and outflowing of a reservoir. 
Lake_in_ou = read.csv(paste0(wd, "/inputs/Lake_in_ou.csv"))
Lake_in_ou = Lake_in_ou[Lake_in_ou$lake_id==lake_id,]
##Reach down, time across file of d_x_area for both inflowing reach (in) and outflowing reach (ou)
d_x_area = read.csv(paste0(wd, "/inputs/d_x_area.csv"))
##Reach down, time across file of slope for both inflowing reach (in) and outflowing reach (ou)
slope = read.csv(paste0(wd, "/inputs/slope.csv"))
##Reach down, time across file of width for both inflowing reach (in) and outflowing reach (ou)
width = read.csv(paste0(wd, "/inputs/width.csv"))
##Lake down, time across file of d_s_q (lake storage change - quadratic method)
delta_s_q = read.csv(paste0(wd, "/inputs/VolumeChange.csv"))
####################################################################################################
##Set up variables for the equation. 
####################################################################################################
##Mean storage change for reservoir of interest. 
mean_dv = mean(abs(as.numeric(delta_s_q[match(Lake_in_ou$lake_id, delta_s_q$lake_id),2:ncol(delta_s_q)])), na.rm = TRUE)
##Slope for inflowing reach. 
sin = sqrt(unlist(slope[match(Lake_in_ou$reach_id_in, slope$reach_id),2:ncol(slope)]))
##Slope for outflowing reach. 
sout = sqrt(unlist(slope[match(Lake_in_ou$reach_id_ou, slope$reach_id),2:ncol(slope)]))
##Width for inflowing reach. 
win = unlist(width[match(Lake_in_ou$reach_id_in, width$reach_id),2:ncol(width)])^(-2/3)
##Width for outflowing reach. 
wout = unlist(width[match(Lake_in_ou$reach_id_ou, width$reach_id),2:ncol(width)])^(-2/3)
##d_x_area for inflowing reach. 
d_x_area_in = unlist(d_x_area[match(Lake_in_ou$reach_id_in, width$reach_id),2:ncol(d_x_area)])
##d_x_area for outflowing reach. 
d_x_area_ou = unlist(d_x_area[match(Lake_in_ou$reach_id_ou, width$reach_id),2:ncol(d_x_area)])
##delta_s_q of lake. 
delta_s_q = unlist(delta_s_q[match(Lake_in_ou$lake_id, delta_s_q$lake_id),2:ncol(delta_s_q)])
################################################################################################################################
##Setting up outputs, t = number of valid estimates kept. 
##Provide a predetermined starting range of values for ni, no, ai, ao. 
################################################################################################################################
t = 200
total = t
rows = 0
outputs = as.data.frame(matrix(numeric(), nrow =total, ncol = 5))
colnames(outputs)= c("vals","ni", "ai","no","ao")
outputs$ni = c(.01,.1)
outputs$no = c(.01,.1)
outputs$ai = c(200,500)
outputs$ao = c(200,2000)
########################################################################################################
##Manning's equation. 
####################################################################################################################
eqn1 = function(n, a, da, w, s){
  flow = (n^-1)*((a+da)^(5/3))*w*s
  return(flow)
}
################################################################################################################
##Lakeflow equation. 
#Randomly estimate non SWOT observable parameters (ni, no, ai, ao) and keep a set of numbers (rw),
#that produce a mean d_s_q within a percentage (x) of the mean observed d_s_q.
############################################################################################################
lakeflow = function(x, rw){
  t = rw
  outputs = get("outputs", environment())
  outputs = outputs[1:t,]
  total = t
  rows = 0
  while(rows<total){
    ni = signif(runif(1, min(outputs$ni), max(outputs$ni)),2)
    no = signif(runif(1, min(outputs$no), max(outputs$no)),2)
    ai = round(runif(1, min(outputs$ai), max(outputs$ai)))
    ao = round(runif(1, min(outputs$ao), max(outputs$ao)))
    inflow = eqn1(ni, ai, d_x_area_in, win, sin)
    outflow = eqn1(no, ao, d_x_area_ou, wout, sout)
    diff = inflow - outflow
    vals = delta_s_q - diff
    vals = unlist(vals)
    if(mean(abs(vals))/mean_dv<x){
      print(rows)
      outputs$vals[rows+1] <- mean(abs(vals))
      outputs$ni[rows+1] = ni
      outputs$ai[rows+1] = ai
      outputs$no[rows+1] = no
      outputs$ao[rows+1] = ao
      rows = rows +1
    } else{next}
  }
return(assign("outputs", as.data.frame(outputs), envir = .GlobalEnv))
}

#####################################################################################################################
##Run lakeflow on a set of values where:
#values:x parameter. 
#rws: rw parameter. 
#As lakeflow runs, the random values of ni, no, ai, ao are created from the min/max of the previous iteration.  
#Processing time: ~0.5-1.5 minutes. 
#################################################################################################################
values = c(.35,.25,.15,.1,.075,.05,.025,.01,.0075,.005,.0025)
rws = c(200,200,200,200,200,100,100,100,100,100,10)
start = Sys.time()
for(i in 1:length(values)){
outputs = lakeflow(values[i], rws[i])
outputs = as.data.frame(outputs)
assign("outputs", as.data.frame(outputs), envir = .GlobalEnv)
}
end = Sys.time()
end - start
ni = median(outputs$ni)
no = median(outputs$no)
ai = median(outputs$ai)
ao = median(outputs$ao)
##################################################################################
##Place data in the netcdf file and write it to outputs folder. 
##################################################################################
id_def = ncdim_def("lake_id", "ID", lake_id)
dimNi <- ncdim_def(name='ni',units = "s/[m1/3]", vals=ni )
dimNo <- ncdim_def(name='no',units = "s/[m1/3]", vals=no )
dimAi <- ncdim_def(name='ai',units = "m2", vals=ai )
dimAo <- ncdim_def(name='ao',units = "m2", vals=ao )
dimID <- ncdim_def(name='lake_id', units='ID', lake_id)
varNi <- ncvar_def(name='ni',units = "s/[m1/3]", dim=list(dimID), missval=-9999)
varNo <- ncvar_def(name='no',units = "s/[m1/3]", dim=list(dimID), missval=-9999)
varAi <- ncvar_def(name='ai',units = "m2", dim=list(dimID), missval=-9999)
varAo <- ncvar_def(name='ao',units = "m2", dim=list(dimID), missval=-9999)
vars <- list(varNi, varNo, varAi, varAo)
con <- nc_create(paste0(wd, "/outputs/output.nc"), vars)
ncatt_put(con, 'lake_id', 'axis', 'ID')
ncvar_put(con, varNi, ni)
ncvar_put(con, varNo, no)
ncvar_put(con, varAi, ai)
ncvar_put(con, varAo, ao)
nc_close(con)




