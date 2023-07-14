functions{
  // ------------------------------------------------------
  //     LINEAR PREDICTOR FOR THE LONGITUDINAL SUBMODEL                
  // ------------------------------------------------------ 
    vector linear_predictor(matrix X, vector obs_times, int[] ID, vector beta_1, matrix U){
      int N = num_elements(obs_times);
      vector[N] out;
      
      out = beta_1[1] + beta_1[2]*obs_times + beta_1[3]*X[ID,1] + U[ID,1] + rows_dot_product(U[ID,2],obs_times);
      
      return out;
    } 
  // ------------------------------------------------------ 
}


data{
  int N;
  int n;
  vector[N] y;
  matrix[n,1] X1;
  int<lower=1,upper=n> ID[N];
  vector[N] obs_times;
}


parameters{
  vector[3] beta_1;
  real<lower=0> sigma_z;
  real<lower=0> sigma_U[2];
  real<lower=-1, upper=1> rho;
  matrix[n,2] U;
}

transformed parameters{
  cov_matrix[2] Sigma;

  Sigma[1,1] = sigma_U[1]^2;
  Sigma[2,2] = sigma_U[2]^2;
  Sigma[1,2] = sigma_U[1]*sigma_U[2]*rho;
  Sigma[2,1] = Sigma[1,2];
}


model{
  // ------------------------------------------------------
    //        LOG-LIKELIHOOD FOR LONGITUDINAL SUBMODEL                
  // ------------------------------------------------------
    {
      vector[N] linpred; 
      
      // Linear predictor
      linpred = linear_predictor(X1, obs_times, ID, beta_1, U);
      
      // Longitudinal Normal log-likelihood
      target += normal_lpdf(y | linpred, sigma_z);
    }  

  // ------------------------------------------------------
    //                       LOG-PRIORS                       
  // ------------------------------------------------------
    // Longitudinal fixed effects
  target += normal_lpdf(beta_1 | 0, 100);
  
  // Random-effects
   for(i in 1:n){ target += multi_normal_lpdf(U[i,1:2] | rep_vector(0.0,2), Sigma); }

   // Random-effects variance
   target += inv_gamma_lpdf(sigma_U | 0.01, 0.01);

   // Random-effects correlation
   target += beta_lpdf((rho+1)/2 | 0.5, 0.5);
  
  // Residual error variance
  target += inv_gamma_lpdf(sigma_z | 0.01, 0.01);   
}