library("MASS") 
library("cmdstanr")
library("pracma")

setwd(paste("C:/Users/ezequ/OneDrive - Fundacao Getulio Vargas - FGV/Mentoria/",
            "joint-models-in-stan", sep = ""))

set.seed(123)

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
  
  times <- id <- long_out <- x_total <- vector()
  
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
  u_3 <- rnorm(N, 0, var_u3)
  
  x <- rnorm(N, 0, 1)
  
  ###################
  # Survival process
  ###################
  
  # Simulating the times to event

  v <- runif(n=N)
  id_times <- c(1:N)
  
  for(i in 1:N){
    haz <- Vectorize(function(s) {
      lambda*rho_s*s^(rho_s-1)*exp(beta_21*x[i] + gamma_1*u_1[i] + 
      gamma_2*u_2[i] + gamma_3*(u_1[i] + u_2[i]*s) + u_3[i])
    })
    
    cum_haz <- Vectorize(function(t) integrate(haz, 0, t)$value)

    sv <- Vectorize(function(t) abs(exp(-cum_haz(t)) - v[i]))

    times[i] <- optim(1e-6, sv, lower = 0, upper = Inf, method = "L-BFGS-B")$par
  }
  
  status <- as.vector(times < cens_time)
  times <- as.vector(ifelse(status, times, cens_time))
  status <- as.numeric(status) # Censoring indicators (1=Observed, 0=Censored)
  
  ##############################
  # Longitudinal process  
  ##############################
  
  obs_times_out <- vector()
  for(i in 1:N){
    # number of repeated observations for each individual
    obs_times <- seq(0,times[i], by = n_rep_obs) 
    
    x_t <- rep(x[i], length(obs_times))

    x_total <- c(x_total,x_t)
    z = rnorm(length(obs_times), 0, sqrt(var_z))
    y_t <- beta_11 + beta_12*obs_times + beta_13*rep(x[i], length(obs_times)) + 
          rep(u_1[i], length(obs_times)) + rep(u_2[i], length(obs_times))*
          obs_times + z

    long_out <- c(long_out,y_t)
    id <- c(id,rep(i,length(obs_times)))
    obs_times_out <- c(obs_times_out,obs_times)
  }
  
  #---------------------------------------------------------------------
  # Creating the longitudinal and survival processes object
  #---------------------------------------------------------------------
  long_proc <- as.matrix(cbind(id, long_out, x_total, obs_times_out)) 
  surv_proc <- as.matrix(cbind(id_times, x, times, status)) 
  data <- list(long_proc, surv_proc)
  names(data) <- c("longitudinal","survival")
  
  return(data)
}

N<- 500
lambda <- 1
rho_s <- 2
cens_time <- 3
beta <- c(0,1,1,1)
gamma <- c(0,0,0)
var_u <- c(0.5,0.5,05)
var_z <- 0.25
rho <- 0.5
n_rep_obs <- 0.5

data <- sim_data(N, 
                 lambda, 
                 rho_s, 
                 cens_time, 
                 beta, 
                 gamma, 
                 var_u, 
                 var_z, 
                 rho, 
                 n_rep_obs)


# Required quantities for longitudinal model fitting

x <- as.matrix(data$survival[,2],1) # unique x
N <- size(x)[1]                     # total number of observations
y <- data$longitudinal[,2]          # longitudinal outcomes
n_obs <- length(y)                  # total number of longitudinal outcomes
id <- data$longitudinal[,1]         # patient IDs
obs_times <- data$longitudinal[,4]  # visit times for repeated observations

long_data = list(y=y,
                 N=N,
                 n_obs=n_obs,
                 x=x,
                 id=id,
                 obs_times=obs_times)
             
long_model <- cmdstan_model("code/long_model.stan")

long_mle <- long_model$optimize(data = long_data)
long_posterior_samples <- long_model$sample(data = long_data, chains = 1)

long_mle$summary()
long_posterior_samples$summary()

# Required quantities for event-time model fitting

x <- as.matrix(data$survival[,2],1)  # unique x
N <- size(x)[1]                      # total number of observations
status <- data$survival[,4]          # vital status (1 = dead, 0 = alive)
times <- data$survival[,3]           # times to event
ind_unc_times <- which(status==1)    # uncensored times indicator
n_unc_times <- length(ind_unc_times) # number of uncensored times

event_data <- list(N=N,
                   x=x,
                   times=times,
                   ind_unc_times=ind_unc_times,
                   n_unc_times=n_unc_times)


event_model <- cmdstan_model("code/event_model.stan")

event_mle <- event_model$optimize(data = event_data)
event_posterior_samples <- event_model$sample(data = event_data, chains = 1)

event_mle$summary()
event_posterior_samples$summary()
