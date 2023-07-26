functions{
  /**
  * Return a vector corresponding to linear predictor for longitudinal
  * sub-model.
  *
  * @param x Matrix of covariates
  * @param obs_times Vector of observed times
  * @param id Vector of integer values corrensponding to identifiers for each
           subject
  * @param beta_1 Vector corresponding to fixed effects
  * @param u Vector corresponding random effects.
  * @return Linear predictor vector.
  */
  vector linear_predictor(matrix x, 
                          vector obs_times, 
                          array[] int id, 
                          vector beta_1, 
                          matrix u){
    int N = num_elements(obs_times);
    vector[N] out;
    
    out = beta_1[1] + beta_1[2]*obs_times + beta_1[3]*x[id,1] + u[id,1] + 
          u[id,2] .* obs_times;
    
    return out;
  }
  
  /**
  * Return a real number corresponding to function to be integrated 
  * evaluated at point t.
  *
  * @param t Value to evaluate function  
  * @param xc A high precision version of the distance from x to the nearest 
  *        endpoint in a definite integral
  * @param theta Parameter values used to evaluate the integral
  * @param x_r Data values used to evaluate the integral
  * @param x_i Integer data used to evaluate the integral.
  * @return The value of the integrand evaluated at the point x.
  */
  real integrand_bp_stable(real t,
                           real xc,                
                           array[] real theta,     
                           array[] real x_r,                        
                           array[] int x_i) {
    real rho_s = theta[1];
    real gamma_3 = theta[2];
    real u_2i = theta[3];
    real lres = rho_s*log(t) + gamma_3*u_2i*t;
    return exp(lres);    
  }
  
  /**
  * Return a vector with each value out[i] corresponding to logarithm hazard 
  * function evaluated at point t[i].
  *
  * @param t Vector corresponding to observed time values  
  * @param x Matrix of covariates
  * @param u_1 Vector corresponding to random effects
  * @param u_2 Vector corresponding to random effects
  * @param u_3 Vector corresponding to random effects
  * @param gamma Vector corresponding to association parameters
  * @param beta_21 Value corresponding to fixed effects
  * @param lambda Scale of the Weibull baseline hazard
  * @param rho_s Shape of the Weibull baseline hazard
  * @param id Array corresponding to subjects to evaluate the function.
  * @return Vector corresponding to logarithm hazard function evaluated at 
  * observed times.
  */
  vector log_haz(vector t, 
                 matrix x, 
                 vector u_1, 
                 vector u_2, 
                 vector u_3, 
                 vector gamma, 
                 real beta_21, 
                 real lambda, 
                 real rho_s,
                 array[] int id){
    int N = num_elements(id);
    vector[N] out;
    real gamma_1 = gamma[1];
    real gamma_2 = gamma[2];
    real gamma_3 = gamma[3];
    out = rep_vector(log(lambda), N) + rep_vector(log(rho_s), N) + (rho_s-1)*
          log(t[id]) + beta_21*x[id,1] + gamma_1*u_1[id] + gamma_2*u_2[id] +
          gamma_3*(u_1[id] + u_2[id] .* t[id]) + u_3[id];
    return out;
  }

  /**
  * Return a vector with each value out[i] corresponding to constant terms 
  * of the cumulative hazard function evaluated at point t[i]
  *
  * @param t Vector corresponding to observed time values  
  * @param x Matrix of covariates
  * @param u_1 Vector corresponding to random effects
  * @param u_2 Vector corresponding to random effects
  * @param u_3 Vector corresponding to random effects
  * @param gamma Vector corresponding to association parameters
  * @param beta_21 Value corresponding to fixed effects
  * @param lambda Scale of the Weibull baseline hazard
  * @param rho_s Shape of the Weibull baseline hazard.
  * @return Vector corresponding to constant terms of the cumulative hazard 
  * function evaluated at observed times.
  */
  vector const_term_cum_haz(vector t, 
                            matrix x, 
                            vector u_1, 
                            vector u_2, 
                            vector u_3, 
                            vector gamma, 
                            real beta_21, 
                            real lambda, 
                            real rho_s){
    int N = num_elements(t);                         
    vector[N] lres;
    real gamma_1 = gamma[1];
    real gamma_2 = gamma[2];
    real gamma_3 = gamma[3];
    lres = rep_vector(log(lambda), N) + rep_vector(log(rho_s), N) + 
           beta_21*x[1:N,1] + gamma_1*u_1 + gamma_2*u_2 + gamma_3*u_1 + u_3;
    return exp(lres);
  }
}


data{
  // Total number of observations
  int n_obs;
  
  // Number of subjects
  int N;
  
  // Number of uncensored times
  int n_unc_times;
  
  // Longitudinal outcomes
  vector[n_obs] y;
  
  // Covariate
  matrix[N,1] x;
  
  // Subjects id
  array[n_obs] int<lower=1,upper=N> id;
  
  // Observed times for longitudinal outcomes
  vector[n_obs] obs_times;
  
  // Times to event
  vector[N] times;
  
  // Uncensored time indices
  array[n_unc_times] int<lower=1,upper=N> ind_unc_times;
}


transformed data {
  array[0] real x_r;
  array[0] int x_i;
}


parameters{
  // Longitudinal fixed effects
  vector[3] beta_1;
  
  // Survival fixed effects
  real beta_21;
  
  // Association parameters
  vector[3] gamma;
  
  // Weibull scale parameter
  real<lower=0> lambda;
  
  // Weibull shape parameter
  real<lower=0> rho_s;
  
  // Measurement errors variance
  real<lower=0> var_z;
  
  // Longitudinal random effects variance
  array[2] real<lower=0> var_u;
  
  // Longitudinal random effects correlation
  real<lower=-1, upper=1> rho;
  
  // Survival random effects variance
  real<lower=0> var_u3;
  
  // Longitudinal random effects
  matrix[N,2] u;
  
  // Survival random effects
  vector[N] u_3;
}

transformed parameters{
  // Covariance matrix of longitudinal random effects
  cov_matrix[2] sigma;

  sigma[1,1] = var_u[1];
  sigma[2,2] = var_u[2];
  sigma[1,2] = sqrt(var_u[1]*var_u[2])*rho;
  sigma[2,1] = sigma[1,2];
}


model{
  // ------------------------------------------------------
  //          LOG-LIKELIHOOD FOR SURVIVAL SUBMODEL                
  // ------------------------------------------------------
  
  // Log-hazard function evaluated at uncensored times
  vector[n_unc_times] lhaz;
  
  // Log-survival function evaluated at times
  vector[N] lsurv;
  
  // Cumulative hazard function evaluated at times
  vector[N] cum_haz;
  
  // Integral for cumulative hazard function
  vector[N] integral;
  
  // Log-hazard function
  lhaz = log_haz(times,
                 x,
                 u[1:N,1],
                 u[1:N,2],
                 u_3,
                 gamma,
                 beta_21,
                 lambda,
                 rho_s,
                 ind_unc_times);
                
  for (i in 1:N) {
    real lfirst_term_bp = rho_s*log(times[i]) + (gamma[3]*u[i,2]*times[i]);
    real integral_bp = integrate_1d(integrand_bp_stable,
                                    0.0,
                                    times[i],
                                    {rho_s, gamma[3], u[i,2]},
                                    x_r,
                                    x_i,
                                    1e-4);
    real second_term_bp = gamma[3] * u[i,2] * integral_bp;
    integral[i] = (1/rho_s)*(exp(lfirst_term_bp) - second_term_bp);
  } 

  cum_haz = const_term_cum_haz(times,
                               x,
                               u[1:N,1],
                               u[1:N,2],
                               u_3,
                               gamma,
                               beta_21,
                               lambda,
                               rho_s) .*
             integral;

  lsurv = -cum_haz;
  
  // Survival log-likelihood
  target += sum(lhaz) + sum(lsurv);
  
  // ------------------------------------------------------
  //        LOG-LIKELIHOOD FOR LONGITUDINAL SUBMODEL                
  // ------------------------------------------------------
  
  vector[n_obs] linpred; 

  // Linear predictor
  linpred = linear_predictor(x, obs_times, id, beta_1, u);

  // Longitudinal Normal log-likelihood
  target += normal_lpdf(y | linpred, sqrt(var_z));

  // ------------------------------------------------------
  //                       LOG-PRIORS                       
  // ------------------------------------------------------
    
  // Longitudinal fixed effects
  target += normal_lpdf(beta_1 | 0, 100);
  
  // Survival fixed effects
  target += normal_lpdf(beta_21 | 1, 5);

  // Weibull scale parameter
  target += gamma_lpdf(lambda | 5, 1);

  // Weibull shape parameter
  target += gamma_lpdf(rho_s | 2, 1);

  // Association parameters
  target += normal_lpdf(gamma | 0, 10);
  
  // Random effects
  for(i in 1:N){ 
    target += multi_normal_lpdf(u[i,1:2] | rep_vector(0.0,2), sigma);
    target += normal_lpdf(u_3[i] | 0, sqrt(var_u3));
  }

  // Random effects variance
  target += inv_gamma_lpdf(var_u | 0.01, 0.01);
  target += inv_gamma_lpdf(var_u3 | 0.01, 0.01);

  // Random effects correlation
  target += beta_lpdf((rho+1)/2 | 0.5, 0.5);
  
  // Residual error variance
  target += inv_gamma_lpdf(var_z | 0.01, 0.01);
}
