functions{
    // ------------------------------------------------------
  //      WEIBULL SURVIVAL SUBMODEL
  // ------------------------------------------------------
    
  
    real integrand(real t,
             real xc,                
             array[] real theta,     
             array[] real x_r,                        
             array[] int x_i) {
      real rho_s = theta[1];
      real gamma_3 = theta[2];
      real U2i = theta[3];
      
      real h = t^rho_s*exp(gamma_3*U2i*t);
        
      return h;    
    }
  
  
  
    // Hazard function
    vector loghaz(vector t, matrix X, vector U1, vector U2, vector U3, vector gamma, real beta_21, real lambda, real rho_s){
         vector[num_elements(t)] out;
         real gamma_1 = gamma[1];
         real gamma_2 = gamma[2];
         real gamma_3 = gamma[3];
         for(i in 1:num_elements(t)){
            out[i] = log(lambda) + log(rho_s) + (rho_s-1)*log(t[i]) + beta_21*X[i,1] + gamma_1*U1[i] + gamma_2*U2[i] + gamma_3*(U1[i] + U2[i]*t[i]) + U3[i];
         }
         return out;
    }

    // Cumulative hazard function
    vector const_term_logcumhaz(vector t, matrix X, vector U1, vector U2, vector U3, vector gamma, real beta_21, real lambda, real rho_s){
         vector[num_elements(t)] out;
         real gamma_1 = gamma[1];
         real gamma_2 = gamma[2];
         real gamma_3 = gamma[3];
         for(i in 1:num_elements(t)){
            out[i] = log(lambda) + log(rho_s) + beta_21*X[i,1] + gamma_1*U1[i] + gamma_2*U2[i] + gamma_3*U1[i];
         }
         return out;
    }
  // ------------------------------------------------------ 
}


data{
  int n;
  int nobs;
  matrix[n,1] X;
  vector[n] times;
  int<lower=1,upper=n> indobs[nobs];
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters{
  real beta_21;
  vector[3] gamma;
  real<lower=0> lambda;
  real<lower=0> rho_s;
  real<lower=0> sigma_U[2];
  real<lower=-1, upper=1> rho;
  matrix[n,2] U;
  real<lower=0> sigma_U3;
  vector[n] U3;
}

transformed parameters{
  cov_matrix[2] Sigma;

  Sigma[1,1] = sigma_U[1];
  Sigma[2,2] = sigma_U[2];
  Sigma[1,2] = sqrt(sigma_U[1]*sigma_U[2])*rho;
  Sigma[2,1] = Sigma[1,2];
}


model {
  // ------------------------------------------------------
//          LOG-LIKELIHOOD FOR SURVIVAL SUBMODEL                
// ------------------------------------------------------
{
   vector[nobs] lhaz;
   vector[nobs] unc_times;
   
   unc_times = times[indobs];
   
   // Log-hazard function
   lhaz = loghaz(unc_times, X[indobs,], U[indobs,1], U[indobs,2], U3[indobs], gamma, beta_21, lambda, rho_s);
   // // Log-survival function
   // lsurv = -cumhaz(times, X, U[1:n,1], U[1:n,2], U3, gamma, beta_21, lambda, rho_s);
   
   vector[n] lsurv;
   
   lsurv = const_term_logcumhaz(times, X, U[1:n,1], U[1:n,2], U3, gamma, beta_21, lambda, rho_s);
   
   for (i in 1:n) {
     lsurv[i] = lsurv[i] + log(times[i]^rho_s*exp(gamma[3]*U[i,2]*times[i])-gamma[3]*U[i,2]*integrate_1d(integrand,
                                                    0.0,
                                                    times[i],
                                                    {rho_s, gamma[3], U[i,2]},
                                                    x_r,
                                                    x_i,
                                                    0.01
                                                    ))
                                                    -log(rho_s);
   }
   
   
   
   target += sum(lhaz) + sum(lsurv); 
   
   
   // PRIORS
   
      // Survival fixed effects
   target += normal_lpdf(beta_21 | 0, 100);

   // PGW scale parameter
   target += cauchy_lpdf(lambda | 0, 1);

   // PGW shape parameters
   target += cauchy_lpdf(rho_s | 0, 1);

   // Association parameters
   target += normal_lpdf(gamma | 0, 100);
   
   // Random-effects
   for(i in 1:n){ 
     target += multi_normal_lpdf(U[i,1:2] | rep_vector(0.0,2), Sigma); 
     target += normal_lpdf(U3[i] | 0, sigma_U3);
  }

   // Random-effects variance
   target += inv_gamma_lpdf(sigma_U | 0.01, 0.01);

   // Random-effects correlation
   target += beta_lpdf((rho+1)/2 | 0.5, 0.5);
}
}

