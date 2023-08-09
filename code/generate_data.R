library("MASS") 
library("cmdstanr")
library("simsurv")

setwd(paste("C:/Users/ezequ/OneDrive - Fundacao Getulio Vargas - FGV/Mentoria/",
            "joint-models-in-stan", sep = ""))

#########################################################################
# Joint model simulation
#########################################################################
# N: sample size
# lambda > 0: scale for Weibull baseline hazard 
# rho_s > 0: shape for Weibull baseline hazard
# cens_time: censored time
# beta: vector of covariates
# gamma: vector of association coefficients
# var_u: vector with random effects variances
# var_z: measurement errors variance
# rho: correlation between u_1 and u_2
# n_rep_obs: number of repeated observations for each subject

sim_data <- function(N, 
                     lambda, 
                     rho_s, 
                     cens_time, 
                     beta, 
                     gamma, 
                     var_u, 
                     var_z, 
                     rho, 
                     n_rep_obs){
  
  id <- long_out <- x_total <- vector()
  
  beta_11 <- beta[1]
  beta_12 <- beta[2]
  beta_13 <- beta[3]
  beta_21 <- beta[4]
  
  gamma_1 <- gamma[1] 
  gamma_2 <- gamma[2]
  gamma_3 <- gamma[3]
  
  mu_u1 <- 0
  mu_u2 <- 0
  mu_u <- c(mu_u1,mu_u2)
  var_u1 <- var_u[1]
  var_u2 <- var_u[2]
  var_u3 <- var_u[3]
  sigma <- matrix(c(var_u1, sqrt(var_u1*var_u2)*rho, 
                    sqrt(var_u1*var_u2)*rho, var_u2), 2)
  bvn_u <- mvrnorm(N, mu = mu_u, Sigma = sigma)
  u_1 <- bvn_u[,1]
  u_2 <- bvn_u[,2]
  u_3 <- rnorm(N, 0, sqrt(var_u3))
  
  x <- rnorm(N, 0, 1)
  
 
  
  ###################
  # Survival process
  ###################
  
  # Simulating the times to event

  id_times <- c(1:N)
  
  covs <- data.frame(id = id_times, x = x, u_1 = u_1, u_2 = u_2, u_3 = u_3)
  
  event_proc <- simsurv(dist = "weibull", lambdas = lambda, gammas = rho_s, 
                        betas = c(x = beta_21, u_1 = gamma_1 + gamma_3, u_2 = 
                                  gamma_2, u_3 = 1),
                        x = covs, tde = c(u_2 = gamma_3), maxt = cens_time)

  times <- event_proc$eventtime
  status <- event_proc$status
  
  ##############################
  # Longitudinal process  
  ##############################
  
  obs_times_out <- vector()
  for(i in 1:N){
    # number of repeated observations for each individual
    obs_times <- seq(0, times[i], by = n_rep_obs) 
    
    
    
    x_t <- rep(x[i], length(obs_times))

    x_total <- c(x_total, x_t)
    z = rnorm(length(obs_times), 0, sqrt(var_z))
    y_t <- beta_11 + beta_12*obs_times + beta_13*rep(x[i], length(obs_times)) + 
          rep(u_1[i], length(obs_times)) + rep(u_2[i], length(obs_times))*
          obs_times + z

    long_out <- c(long_out, y_t)
    id <- c(id,rep(i, length(obs_times)))
    obs_times_out <- c(obs_times_out, obs_times)
  }
  
  #---------------------------------------------------------------------
  # Creating the longitudinal and survival processes object
  #---------------------------------------------------------------------
  N <- length(id_times)                # number of subjects
  n_obs <- length(long_out)            # total number of observations
  x <- as.matrix(x,1)                  # unique x
  obs_times <- obs_times_out           # visit times for repeated observations
  y <- long_out                        # longitudinal outcomes
  ind_unc_times <- which(status==1)    # uncensored times indicator
  n_unc_times <- length(ind_unc_times) # number of uncensored times
  
  long_data <- list(y=y,
                   N=N,
                   n_obs=n_obs,
                   x=x,
                   id=id,
                   obs_times=obs_times)
  
  event_data <- list(N=N,
                     x=x,
                     times=times,
                     ind_unc_times=ind_unc_times,
                     n_unc_times=n_unc_times)
  
  joint_data <- list(N=N,
                     n_obs=n_obs,
                     y=y,
                     id=id,
                     obs_times=obs_times,
                     x=x,
                     times=times,
                     ind_unc_times=ind_unc_times,
                     n_unc_times=n_unc_times)
  
  save(long_data, file = "data/long_data.RData")
  save(event_data, file = "data/event_data.RData")
  save(joint_data, file = "data/joint_data.RData")
}