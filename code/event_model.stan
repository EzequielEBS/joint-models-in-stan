functions{
  /**
  * Return a real number corresponding to function function to be integrated 
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
  real integrand(real t, 
                 real xc, 
                 array[] real theta, 
                 array[] real x_r, 
                 array[] int x_i) {
    real rho_s = theta[1];
    real gamma_3 = theta[2];
    real u_2i = theta[3];
    
    real h = t^(rho_s-1)*exp(gamma_3*u_2i*t);
      
    return h;    
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
  real integrand_bp(real t,
                    real xc,                
                    array[] real theta,     
                    array[] real x_r,                        
                    array[] int x_i) {
    real rho_s = theta[1];
    real gamma_3 = theta[2];
    real u_2i = theta[3];
    
    real h = t^rho_s*exp(gamma_3*u_2i*t);
      
    return h;    
  }

  /**
  * Return a Laplace approximation of specific function.
  *
  * @param a Lower limit to evaluate integral  
  * @param b Upper limit to evaluate integral
  * @param rho_s Shape of the Weibull baseline hazard
  * @param gamma_3 3th association parameter
  * @param u_2i ith value of the second random effect.
  * @return The integral approximation.
  */
  real lap_app(real a, real b, real rho_s, real gamma_3, real u_2i){
    real t0 = (1-rho_s)/(gamma_3*u_2i);
    real h_t0 = (rho_s-1)*log(t0) + gamma_3*u_2i*t0;
    real ddh_t0 = (1-rho_s)/(t0^2);
    real lout = h_t0 + 0.5*(log(2) + log(pi()) - log(-ddh_t0)) +
                log(normal_cdf(b, t0, -ddh_t0) -
                normal_cdf(a, t0, -ddh_t0));
    // real out = exp(h_t0)*sqrt(2*pi()/(-ddh_t0))*
    //            (normal_cdf(b, t0, -ddh_t0) -
    //            normal_cdf(a, t0, -ddh_t0));
    return exp(lout);
  }

  /**
  * Return a lower incomplete gamma function aproximation.
  *
  * @param a Integrand parameter  
  * @param z Upper limit to evaluate integral
  * @param tol Aproximation tolerance
  * @return The approximation.
  */
  real low_inc_gamma(real a, real z, real tol) {
    real inc_gamma;
    
    if (z < -50) {
      int n = 0;
      real inc_gamma0 = (1-a)/(z^n);
      real inc_gamma1 = inc_gamma0 + 1;
      while (abs(inc_gamma1 - inc_gamma0) > tol){
        n += 1;
        inc_gamma1 = inc_gamma0;
        inc_gamma0 += (1-a)/(z^n);
      }
      inc_gamma = exp(z)/(z*tgamma(a))*inc_gamma0;
    }
    else {
      int n = 0;
      real inc_gamma0 = (-z)^n/(tgamma(n+1)*(a+n));
      real inc_gamma1 = inc_gamma0 + 1;
      while (abs(inc_gamma1 - inc_gamma0) > tol){
        n += 1;
        inc_gamma1 = inc_gamma0;
        inc_gamma0 += (-z)^n/(tgamma(n+1)*(a+n));
      }
      inc_gamma = 1/(tgamma(a))*inc_gamma0;
    }
    real out = z^a*tgamma(a)*inc_gamma;
    return out;
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
          gamma_3*(u_1[id] + rows_dot_product(u_2[id],t[id])) + u_3[id];
    return out;
  }

  /**
  * Return a vector with each value out[i] corresponding to constant terms of
  * the cumulative hazard function evaluated at point t[i]
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
    vector[N] out;
    real gamma_1 = gamma[1];
    real gamma_2 = gamma[2];
    real gamma_3 = gamma[3];
    out = lambda*rho_s*exp(beta_21*x[1:N,1] + gamma_1*u_1 + 
          gamma_2*u_2 + gamma_3*u_1 + u_3);
    return out;
  }
  
  /**
  * Return a vector with each value out[i] corresponding to log constant terms 
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
  * @return Vector corresponding to log constant terms of the cumulative hazard 
  * function evaluated at observed times.
  */
    
  vector lconst_term_cum_haz(vector t, 
                             matrix x, 
                             vector u_1, 
                             vector u_2, 
                             vector u_3, 
                             vector gamma, 
                             real beta_21, 
                             real lambda, 
                             real rho_s){
    int N = num_elements(t);
    vector[N] out;
    real gamma_1 = gamma[1];
    real gamma_2 = gamma[2];
    real gamma_3 = gamma[3];
    out = rep_vector(log(lambda), N) + rep_vector(log(rho_s), N) +
           beta_21*x[1:N,1] + gamma_1*u_1 + gamma_2*u_2 + gamma_3*u_1 + u_3;
    return out;
  }
}


data{
  int N;
  int n_unc_times;
  matrix[N,1] x;
  vector[N] times;
  array[n_unc_times] int<lower=1,upper=N> ind_unc_times;
}

transformed data {
  array[0] real x_r;
  array[0] int x_i;
}

parameters{
  real beta_21;
  vector[3] gamma;
  real<lower=0> lambda;
  real<lower=0> rho_s;
  array[2] real<lower=0> var_u;
  real<lower=-1, upper=1> rho;
  real<lower=0> var_u3;
  matrix[N,2] u;
  vector[N] u_3;
}

transformed parameters{
  cov_matrix[2] sigma;

  sigma[1,1] = var_u[1];
  sigma[2,2] = var_u[2];
  sigma[1,2] = sqrt(var_u[1]*var_u[2])*rho;
  sigma[2,1] = sigma[1,2];
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
                u[1:N,1],
                u[1:N,2],
                u_3,
                gamma,
                beta_21,
                lambda,
                rho_s,
                ind_unc_times);
                
  // Integration by parts

  // Log-survival function
  lsurv = -const_term_cum_haz(times,
                             x,
                             u[1:N,1],
                             u[1:N,2],
                             u_3,
                             gamma,
                             beta_21,
                             lambda,
                             rho_s);

  for (i in 1:N) {
    real integral = integrate_1d(integrand_bp,
                                 0.0,
                                 times[i],
                                 {rho_s, gamma[3], u[i,2]},
                                 x_r,
                                 x_i,
                                 1e-8);
    lsurv[i] *= (1/rho_s)*(times[i]^rho_s*exp(gamma[3]*u[i,2]*times[i]) -
                gamma[3]*u[i,2]);
  }

  // // Integration by parts with log
  // 
  // // Log-survival function
  // lsurv = lconst_term_cum_haz(times,
  //                             x,
  //                             u[1:N,1],
  //                             u[1:N,2],
  //                             u_3,
  //                             gamma,
  //                             beta_21,
  //                             lambda,
  //                             rho_s);
  // 
  // for (i in 1:N) {
  //   real integral = integrate_1d(integrand_bp,
  //                                0.0,
  //                                times[i],
  //                                {rho_s, gamma[3], u[i,2]},
  //                                x_r,
  //                                x_i,
  //                                0.01);
  //   lsurv[i] += -log(rho_s) +
  //               log((times[i]^rho_s*exp(gamma[3]*u[i,2]*times[i]) -
  //               gamma[3]*u[i,2]*integral));
  //   lsurv[i] = -exp(lsurv[i]);
  // }

  // Integration using incomplete gamma approximation
  // for (i in 1:N) {
  //   if (!is_nan(low_inc_gamma(rho_s, -gamma[3]*u[i,2]*times[i], 0.01))){
  //     lsurv[i] *= (-gamma[3]*u[i,2])^(-rho_s)*
  //                 low_inc_gamma(rho_s, -gamma[3]*u[i,2]*times[i], 0.01);
  //   }
  //   else {
  //     real integral = integrate_1d(integrand_bp,
  //                                  0.0,
  //                                  times[i],
  //                                  {rho_s, gamma[3], u[i,2]},
  //                                  x_r,
  //                                  x_i,
  //                                  0.01);
  //     lsurv[i] *= (1/rho_s)*
  //                 (times[i]^rho_s*exp(gamma[3]*u[i,2]*times[i]) -
  //                 gamma[3]*u[i,2]*integral);
  //   }
  // }

  // Direct integral
  // for (i in 1:N) {
  //   real integral = integrate_1d(integrand,
  //                                0.0,
  //                                times[i],
  //                                {rho_s, gamma[3], u[i,2]},
  //                                x_r,
  //                                x_i,
  //                                0.01);
  //   lsurv[i] = lsurv[i]*integral;
  // }
  
  // // Laplace Approximation
  // for (i in 1:N) {
  //   if (rho_s > 1){
  //     lsurv[i] *= lap_app(0, times[i], rho_s, gamma[3], u[i,2]);
  //   }
  //   else {
  //     real integral = integrate_1d(integrand_bp,
  //                                  0.0,
  //                                  times[i],
  //                                  {rho_s, gamma[3], u[i,2]},
  //                                  x_r,
  //                                  x_i,
  //                                  0.01);
  //     lsurv[i] *= (1/rho_s)*
  //                 (times[i]^rho_s*exp(gamma[3]*u[i,2]*times[i]) -
  //                 gamma[3]*u[i,2]*integral);
  //   }
  // }
   
  // Survival log-likelihood
  target += sum(lhaz) + sum(lsurv);
   
   
  // ------------------------------------------------------
  //          PRIOR DISTRIBUTIONS
  // ------------------------------------------------------

  // Survival fixed effects
  target += normal_lpdf(beta_21 | 1, 5);

  // Weibull scale parameter
  target += gamma_lpdf(lambda | 1, 1);

  // Weibull shape parameter
  target += gamma_lpdf(rho_s | 2, 1);

  // Association parameters
  target += normal_lpdf(gamma | 0, 10);

  // Random effects
  for(i in 1:N){
    target += multi_normal_lpdf(u[i,1:2] | rep_vector(0.0,2), sigma);
    target += normal_lpdf(u_3[i] | 0, var_u3);
  }

  // Random effects variance
  target += inv_gamma_lpdf(var_u | 0.01, 0.01);
  target += inv_gamma_lpdf(var_u3 | 0.01, 0.01);

  // Random effects correlation
  target += beta_lpdf((rho+1)/2 | 2, 1);
}
