data {
  int<lower=0> N;
  matrix[N,3]  y;
  matrix[N,6]  x;
  vector[N]   ones;       // technical support vector of 1's
  
}
parameters {
  vector[4]  tx;          //intercepts x
  vector[2]  ty;          //intercepts y
  vector[2]  muxi;        //means xi
  vector<lower=0>[4]  lx; //factor loadings x
  vector<lower=0>[2]  ly; //factor loadings y
  vector[4]           b1; //regression coefficients

  vector<lower=0>[6]      sigmax;   //residual sd x
  vector<lower=0>[3]      sigmay;   //residual sd y
  vector<lower=0>[2]      sigmaxi;  //sd xi
  real<lower=0>           sigmaeta; //residual sd eta
  cholesky_factor_corr[2] L1;       //choleski correlation matrix xi
  
  vector[N]  eta;              //dv   
  matrix[N,2] zi;              //iv choleski transformed
}
transformed parameters {
  
  
}

model {
  matrix[N,2] xi;            // transformed iv original metric
  matrix[N,6] mux;           // E[x|xi]
  matrix[N,3] muy;           // E[y|eta]
  vector[N]   mueta;         // E[eta|xi]
  matrix[N,2] muxi2;         // transformed mean in matrix
  
  muxi2[,1] = muxi[1]*ones;
  muxi2[,2] = muxi[2]*ones;
  
  // cholesky re-transformation
  xi = muxi2+ zi*diag_pre_multiply(sigmaxi,L1)';               
  
  ///////////////////////////////////
  // measurement models
  ///////////////////////////////////
  mux[,1] =             xi[,1];
  mux[,2] = tx[1]+lx[1]*xi[,1];
  mux[,3] = tx[2]+lx[2]*xi[,1];

  mux[,4] =             xi[,2];
  mux[,5] = tx[3]+lx[3]*xi[,2];
  mux[,6] = tx[4]+lx[4]*xi[,2];

  muy[,1] =             eta;
  muy[,2] = ty[1]+ly[1]*eta;
  muy[,3] = ty[2]+ly[2]*eta;
  
  ///////////////////////////////////
  // structural model
  ///////////////////////////////////
  mueta = b1[1]+b1[2]*xi[,1]+b1[3]*xi[,2]+b1[4]*(xi[,1].*xi[,2]);

  ///////////////////////////////////
  // distribution of observed variables  
  ///////////////////////////////////
  for(k in 1:6){x[,k] ~ normal(mux[,k],sigmax[k]);} 
  for(k in 1:3){y[,k] ~ normal(muy[,k],sigmay[k]);} 
  
  ///////////////////////////////////
  // distribution of latent variables
  ///////////////////////////////////
  eta ~ normal(mueta,sigmaeta);                      
  to_vector(zi) ~ normal(0,1);    // cholesky version

  ///////////////////////////////////
  // prior distributions
  ///////////////////////////////////
  // priors intercepts and means
  tx ~ normal(0,1);
  ty ~ normal(0,1);
  muxi ~ normal(0,1);
  
  // priors factor loadings and regression coefficients
  lx ~ normal(0,1);
  ly ~ normal(0,1);
  b1 ~ normal(0,1);
  
  // priors for (residual) SDs 
  sigmax ~ cauchy(0,2.5);   
  sigmay ~ cauchy(0,2.5);   
  sigmaxi ~ cauchy(0,2.5); 
  sigmaeta ~ cauchy(0,2.5); 
  
  // LKJ prior for correlation matrix
  L1 ~ lkj_corr_cholesky(1);

}

generated quantities{
  // covariance matrix xi
  matrix[2,2] phi;
  phi = diag_pre_multiply(sigmaxi,L1)*diag_pre_multiply(sigmaxi,L1)';
}




