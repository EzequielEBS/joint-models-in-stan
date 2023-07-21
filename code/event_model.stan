functions{
  /**
  * Return a vector with each value out[i] corresponding to logarithm hazard 
  * function evaluated at point t[i].
  *
  * @param t Vector corresponding to observed time values  
  * @param x Matrix of covariates
  * @param u_3 Vector corresponding to random effects
  * @param beta_21 Value corresponding to fixed effects
  * @param lambda Scale of the Weibull baseline hazard
  * @param rho_s Shape of the Weibull baseline hazard
  * @param id Array corresponding to subjects to evaluate the function.
  * @return Vector corresponding to logarithm hazard function evaluated at 
  * observed times.
  */
  vector log_haz(vector t, 
                 matrix x,
                 vector u_3, 
                 real beta_21, 
                 real lambda, 
                 real rho_s,
                 array[] int id){
    int N = num_elements(id);
    vector[N] out;
    out = rep_vector(log(lambda), N) + rep_vector(log(rho_s), N) + (rho_s-1)*
          log(t[id]) + beta_21*x[id,1] + u_3[id];
    return out;
  }

  /**
  * Return a vector with each value out[i] corresponding to constant terms of
  * the cumulative hazard function evaluated at point t[i]
  *
  * @param t Vector corresponding to observed time values  
  * @param x Matrix of covariates
  * @param u_3 Vector corresponding to random effects
  * @param beta_21 Value corresponding to fixed effects
  * @param lambda Scale of the Weibull baseline hazard
  * @param rho_s Shape of the Weibull baseline hazard.
  * @return Vector corresponding to constant terms of the cumulative hazard 
  * function evaluated at observed times.
  */
  vector cum_haz(vector t, 
                 matrix x,  
                 vector u_3, 
                 real beta_21, 
                 real lambda, 
                 real rho_s){
    int N = num_elements(t);                         
    vector[N] lout;
    lout = rep_vector(log(lambda), N) + rho_s*log(t) + beta_21*x[1:N,1] + u_3;
    return exp(lout);
  }
}

data{
  int N;
  int n_unc_times;
  matrix[N,1] x;
  vector[N] times;
  array[n_unc_times] int<lower=1,upper=N> ind_unc_times;
}

parameters{
  real beta_21;
  real<lower=0> lambda;
  real<lower=0> rho_s;
  real<lower=0> var_u3;
  vector[N] u_3;
}

model {
  // ------------------------------------------------------
  //          LOG-LIKELIHOOD FOR SURVIVAL SUBMODEL                
  // ------------------------------------------------------
  
  vector[n_unc_times] lhaz;

  vector[N] lsurv;
  
  // Log-hazard function
  lhaz = log_haz(times,
                 x,
                 u_3,
                 beta_21,
                 lambda,
                 rho_s,
                 ind_unc_times);

  // Log-survival function
  lsurv = -cum_haz(times,
                   x,
                   u_3,
                   beta_21,
                   lambda,
                   rho_s);
   
  // Survival log-likelihood
  target += sum(lhaz) + sum(lsurv);
   
   
  // ------------------------------------------------------
  //          PRIOR DISTRIBUTIONS
  // ------------------------------------------------------

  // Survival fixed effects
  target += normal_lpdf(beta_21 | 1, 1);

  // Weibull scale parameter
  target += gamma_lpdf(lambda | 1, 1);

  // Weibull shape parameter
  target += gamma_lpdf(rho_s | 2, 1);

  // Random effects
  for(i in 1:N){
    target += normal_lpdf(u_3[i] | 0, sqrt(var_u3));
  }

  // Random effects variance
  target += inv_gamma_lpdf(var_u3 | 0.01, 0.01);
}
