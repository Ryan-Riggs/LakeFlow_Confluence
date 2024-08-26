//Make sure et/lateral are/not in log space. 
data {
  // Dimensions
  int N; //number of time steps
  int n1; //number of inflows
  int n2; //number of outflows
  
  // bounds on parameters - write as vectors to allow for multiple inflows/outflows. 
  vector[n1] nInlower;
  vector[n1] nInupper;
  vector[n2] nOutlower;
  vector[n2] nOutupper;
  vector[n1] aInlower;
  vector[n1] aInupper;
  vector[n2] aOutlower;
  vector[n2] aOutupper;
  vector[n1] daInShift;
  vector[n2] daOutShift;
  vector[n1] nInHat;
  vector[n1] nInSd;
  vector[n1] aInHat;
  vector[n1] aInSd;
  vector[n1] qInupper;
  vector[n1] qInlower;
  vector[n2] nOutHat;
  vector[n2] nOutSd;
  vector[n2] aOutHat;
  vector[n2] aOutSd;
  vector[n2] qOutupper;
  vector[n2] qOutlower;
  
  vector[N] sigmaIn[n1];
  vector[N] sigmaOut[n2];
  vector[N] qInSd[n1];
  vector[N] qOutSd[n2];
  vector[N] q[n1];
  
  vector[N] da[n1];
  vector[N] w[n1];
  vector[N] s[n1];
  vector[N] da2[n2];
  vector[N] w2[n2];
  vector[N] s2[n2];
  vector[N] q2[n2];
  vector[N] dv_per;
  vector[N] et;
  vector[N] lateral;
}

parameters {
  vector < lower=nInlower[n1], upper=nInupper[n1] >[n1] n; // mannning's n
  vector < lower=aInlower[n1], upper=aInupper[n1] >[n1] a; // Bathymetry
  vector < lower=nOutlower[n2], upper=nOutupper[n2] >[n2] nOut; // mannning's n
  vector < lower=aOutlower[n2], upper=aOutupper[n2] >[n2] aOut; // Bathymetry
  real < lower = 0 > sigma; // Error SD
  vector < lower = qInlower[n1], upper=qInupper[n1]>[N] logQ_in[n1]; //Inflow Q timeseries
  vector < lower = qOutlower[n2], upper=qOutupper[n2]>[N] logQ_out[n2]; //Outflow Q timeseries
  
  //vector[N] logQ_in[n1]; //Inflow Q timeseries
  //vector[N] logQ_out[n2]; //Outflow Q timeseries
}

transformed parameters {
  vector[N] lhsIn[n1]; // LHS for Manning likelihood
  vector[N] rhsIn[n1]; // log area for Manning's equation
  vector[N] lhsOut[n2]; // LHS for Manning likelihood
  vector[N] rhsOut[n2]; // log area for Manning's equation
  vector[N] lhsDV; // LHS for Manning likelihood
  vector<lower = sum(exp(qInlower)), upper = sum(exp(qInupper))>[N] sumIn;
  vector<lower = sum(exp(qOutlower)), upper = sum(exp(qOutupper))>[N] sumOut;
  //vector[N] sumIn;
  //vector[N] sumOut;
  vector[N] rhsDV; // log area for Manning's equation
  
  vector[N] a1[n1];
  vector[N] a2[n2];
  
  //For each inflow solve Manning's
  for(j in 1:n1){
    a1[j,] = log(da[j,]+a[j]);
    lhsIn[j,] = (4*w[j,])-(3*s[j,]);
    rhsIn[j,] = ((-6*n[j])+(10*a1[j,]))-(6*logQ_in[j,]);
  }

  //For each outflow solve Manning's
  for(k in 1:n2){
    a2[k,] = log(da2[k,]+aOut[k]);
    lhsOut[k,] = (4*w2[k,])-(3*s2[k,]);
    rhsOut[k,] = ((-6*nOut[k])+(10*a2[k,]))-(6*logQ_out[k,]);
  }

  
  lhsDV = (dv_per)-lateral+et;
  
  //For each timestep, sum inflows and outflows and then subtract them. 
  for(i in 1:N){
    sumIn[i] = sum(exp(logQ_in[,i]));
    sumOut[i] = sum(exp(logQ_out[,i]));
    rhsDV[i] = sumIn[i] - sumOut[i];
    //rhsDV[i] = sum(exp(logQ_in[,i])) - sum(exp(logQ_out[,i]));
  }

}


model {
  
  for(x in 1:n1){
    logQ_in[x,]~normal(q[x,],qInSd[x,]);
    n[x] ~ normal(nInHat[x],nInSd[x]);
    a[x]+daInShift[x] ~lognormal(aInHat[x], aInSd[x]);
    lhsIn[x,]~normal(rhsIn[x,], sigmaIn[x,]);
  }

  for(y in 1:n2){
    logQ_out[y,]~normal(q2[y,],qOutSd[y,]);
    nOut[y] ~ normal(nOutHat[y], nOutSd[y]);
    aOut[y] + daOutShift[y]~lognormal(aOutHat[y], aOutSd[y]);
    lhsOut[y,]~normal(rhsOut[y,], sigmaOut[y,]);
  }
  
  lhsDV~normal(rhsDV, sigma);
}

