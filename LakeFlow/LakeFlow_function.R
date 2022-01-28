################################################################################################################################
#Lakeflow function.  
#Ryan Riggs
#1/28/2022
#########################################################################################################
##Algorithm. 
#########################################################################################################
LakeFlow = function(win, sin, d_x_area_in, wout, sout, 
                    d_x_area_out,dA_lw_inflow, dA_up_inflow, 
                    dA_lw_outflow, dA_up_outflow, n_lw_inflow, 
                    n_up_inflow, n_lw_outflow, n_up_outflow,
                    et,sum,delta_s_q, time){
#Manning's equation. 
eqn1 = function(n, a, da, w, s){
  flow = (n^-1)*((a+da)^(5/3))*w*s
  return(flow)
}

##Objective function. 
fun = function(ni, no, ai, ao){
  
  inflow = (ni^-1)*((ai+d_x_area_in)^(5/3))*win*sin
  outflow = (no^-1)*((ao+d_x_area_out)^(5/3))*wout*sout
  
  inflow = inflow +sum
  outflow = outflow + et
  
  a = (delta_s_q-(inflow-outflow))
  a = sum(log(abs(a)))
  Score = -1*a
}
postfit <- function(object, ...)
{
  pop <- object@population
  # update info
  if(!exists(".pop", envir = globalenv()))
    assign(".pop", NULL, envir = globalenv())
  .pop <- get(".pop", envir = globalenv())
  assign(".pop", append(.pop, list(pop)), envir = globalenv()) 
  # output the input ga object (this is needed!!)
  object 
}

##Genetic Algorithm. 
t = GA::ga(type = "real-valued",
       fitness =function(x) fun(x[1], x[2], x[3], x[4]),
       lower = c(n_lw_inflow,n_lw_outflow,dA_lw_inflow,dA_lw_outflow),
       upper = c(n_up_inflow,n_up_outflow,dA_up_inflow,dA_up_outflow),
       popSize= 10000,
       maxiter =200,
       run = 50,
       parallel = FALSE,
       postFitness = postfit,
       elitism = round(length(win)*.20),
       seed = 357)
##Algorithm solution. 
outputs = summary(t)$solution

Q_inflow = (outputs[1]^-1)*((outputs[3]+d_x_area_in)^(5/3))*win*sin
Q_outflow =(outputs[2]^-1)*((outputs[4]+d_x_area_out)^(5/3))*wout*sout


##output. 
df = data.frame(Date = time, n_inflow = outputs[1], A0_inflow = outputs[3],
                n_outflow = outputs[2], A0_outflow = outputs[4], n_inflow_sigma = -9999, A0_inflow_sigma = -9999,
                n_outflow_sigma = -9999, A0_outflow_sigma = -9999, Q_inflow = Q_inflow, Q_outflow = Q_outflow)
return(df)
}                                                                                              











