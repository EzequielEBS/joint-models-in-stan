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
}


data{
  // Total number of observations
  int n_obs;
  
  // Number of subjects
  int N;
  
  // Longitudinal outcomes
  vector[n_obs] y;
  
  // Covariate
  matrix[N,1] x;
  
  // Subjects id
  array[n_obs] int<lower=1,upper=N> id;
  
  // Observed times for longitudinal outcomes
  vector[n_obs] obs_times;
}


parameters{
  // Longitudinal fixed effects
  vector[3] beta_1;
  
  // Measurement errors variance
  real<lower=0> var_z;
  
  // Longitudinal random effects variance
  array[2] real<lower=0> var_u;
  
  // Longitudinal random effects correlation
  real<lower=-1, upper=1> rho;
  
  // Longitudinal random effects
  matrix[N,2] u;
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
  
  // Random effects
  for(i in 1:N){ 
    target += multi_normal_lpdf(u[i,1:2] | rep_vector(0.0,2), sigma);
  }

  // Random effects variance
  target += inv_gamma_lpdf(var_u | 0.01, 0.01);

  // Random effects correlation
  target += beta_lpdf((rho+1)/2 | 0.5, 0.5);
  
  // Residual error variance
  target += inv_gamma_lpdf(var_z | 0.01, 0.01);
}
