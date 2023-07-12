functions{
// // ------------------------------------------------------
// //      WEIBULL SURVIVAL SUBMODEL                
// // ------------------------------------------------------
//     // Hazard function
//     vector loghaz(vector t, vector X, vector U1, vector U2, vector U3, vector gamma, real beta_2 real lambda, real rho_s){
//          vector[num_elements(t)] out;
//          real gamma_1 = gamma[1];
//          real gamma_2 = gamma[2];
//          real gamma_3 = gamma[3];
//          for(i in 1:num_elements(t)){
//             out[i] = log(lambda) + log(rho_s) + (rho_s-1)*log(t[i]) + beta_21*X[i] + gamma_1*U1[i] + gamma_2*U2[i] + gamma_3(U1[i] + U2[i]*t[i]) + U3[i];
//          }
//          return out;
//     }                                                                                     
// 
//     // Cumulative hazard function
//     vector cumhaz(vector t, real eta, real nu, real delta){
//          vector[num_elements(t)] out;  
//          for(i in 1:num_elements(t)){
//             out[i] = lambda*rho_s*exp(beta_21*X[i] + gamma_1*U1[i] + gamma_2*U2[i] + gamma_3*(U1[i]) + U3[i])*(-gamma_3*U2[i])^(-rho_s)*pracma::gammainc(-gamma_3*U2[i]*t[i], rho_s)[[1]]);
//          }
//          return out;
//     }
// ------------------------------------------------------ 

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

  Sigma[1,1] = sigma_U[1];
  Sigma[2,2] = sigma_U[2];
  Sigma[1,2] = sqrt(sigma_U[1]*sigma_U[2])*rho;
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
      target += normal_lpdf(y | linpred, sqrt(sigma_z));
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