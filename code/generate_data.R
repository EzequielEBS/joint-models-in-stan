library("MASS") 
library("cmdstanr")
library("pracma")

setwd("C:/Users/ezequ/OneDrive - Fundacao Getulio Vargas - FGV/Mentoria/")

set.seed(123)

#########################################################################
# Joint model simulation
#########################################################################
# m   : sample size
# lambda > 0: scale for Weibull baseline hazard 
# rho > 0: shape for Weibull baseline hazard
# cens_time: censored time
# beta: vector of covariates
# gamma: vector of association coefficients
# sigma_U: vector with variances of U1, U2 and U3
# sigma_z: standard deviation of measurement errors 
# rho: correlation between U1 and U2
# n_rep_obs: number of repeated observations for each individual



simDataJ <- function(m, lambda, rho_s, cens_time, beta, gamma, sigma_U, sigma_z, rho, n_rep_obs){
  
  times <- ID <- longit.out <- X_total <- vector()
  
  beta_11 <- beta[1]
  beta_12 <- beta[2]
  beta_13 <- beta[3]
  beta_21 <- beta[4]
  
  gamma_1 <- gamma[1] 
  gamma_2 <- gamma[2]
  gamma_3 <- gamma[3]
  
  mu_U1 <- 0
  mu_U2 <- 0
  mu_U <- c(mu_U1,mu_U2)
  sigma_U1 <- sigma_U[1]
  sigma_U2 <- sigma_U[2]
  sigma_U3 <- sigma_U[3]
  sigma_U <- matrix(c(sigma_U1^2, sigma_U1*sigma_U2*rho, sigma_U1*sigma_U2*rho, sigma_U2^2),
                    2)
  bvnU <- mvrnorm(m, mu = mu_U, Sigma = sigma_U)
  U1 <- bvnU[,1]
  U2 <- bvnU[,2]
  U3 <- rnorm(m, 0, sigma_U3)
  
  X <- rnorm(m, 0, 1)
  
  ###################
  # Survival process
  ###################
  
  # Simulating the times to event

  v <- runif(n=m)
  id_times <- c(1:m)
  
  for(i in 1:m){
    h_0_not_t <- lambda*rho_s
    h_not_t <- h_0_not_t*exp(beta_21*X[i] + gamma_1*U1[i] + gamma_2*U2[i] + gamma_3*(U1[i]) + U3[i])
    
    # add case gamma_3 == 0
    
    # result of integrate h(t). It works but it is not correct, because it needs -gamma_3*U2[i]*u >= 0
    H <- Vectorize(function(u) h_not_t*(-gamma_3*U2[i])^(-rho_s)*pracma::gammainc(-gamma_3*U2[i]*u, rho_s)[[1]])
    
    
    Sv <- Vectorize(function(u) abs(exp(-H(u)) - v[i]))
    
    # find t that minimizes absolute error: |S_i(t) - v[i]| (may not be a root)
    times[i] <- optim(0, Sv, lower = 0, upper = Inf)$par
    
  }
  status <- as.vector(times < cens_time)
  times <- as.vector(ifelse(status, times, cens_time))
  status <- as.numeric(status) # Censoring indicators (1=Observed, 0=Censored)
  
  ##############################
  # Longitudinal process  
  ##############################
  
  obs_times.out <- vector()
  for(i in 1:m){
    obs_times <- seq(0,times[i], by = n_rep_obs) # number of repeated observations for each individual
    
    X_t <- rep(X[i], length(obs_times))

    X_total <- c(X_total,X_t)
    Z = rnorm(length(obs_times), 0, sigma_z)
    yt <- beta_11 + beta_12*obs_times + beta_13*rep(X[i], length(obs_times)) + rep(U1[i], length(obs_times)) + rep(U2[i], length(obs_times))*obs_times + Z

    longit.out <- c(longit.out,yt)
    ID <- c(ID,rep(i,length(obs_times)))
    obs_times.out <- c(obs_times.out,obs_times)
  }
  
  #---------------------------------------------------------------------
  # Creating the longitudinal and survival processes object
  #---------------------------------------------------------------------
  long.proc <- as.matrix(cbind(ID, longit.out, X_total)) # Longitudinal process
  surv.proc <- as.matrix(cbind(id_times, times, status)) # Survival process
  obj <- list(long.proc,surv.proc,X_total,X, obs_times.out)
  names(obj) <- c("longitudinal","survival","X_total", "X", "obs_times")
  
  return(obj)
}


m<- 250
lambda <- 0.01
rho_s <- 1
cens_time <- 4
beta <- c(0,1,1,1)
gamma <- c(-1.5,0,2)
sigma_U <- c(0.5^0.5,1,0.25^0.5)
sigma_z <- 0.25^0.5
rho <- 0.3
n_rep_obs <- 0.5

obj <- simDataJ(m, lambda, rho_s, cens_time, beta, gamma, sigma_U, sigma_z, rho, n_rep_obs)


# Required quantities for model fitting
X <- obj$X                         # unique X    
X_total <- obj$longitudinal[,3]    # X with repeated observations
n <- nrow(X)                       # total number of observations
y <- obj$longitudinal[,2]          # longitudinal outcomes
ID <- obj$longitudinal[,1]         # patient IDs
nid <- length(unique(ID))          # number of patients
id_times <- obj$survival[,1]       # unique ids
status <- obj$survival[,3]         # vital status (1 = dead, 0 = alive)
times <- obj$survival[,2]          # times to event
obs_times <- obj$obs_times         # visit times for repeated observations
N <- length(y)                     # total number of longitudinal outcomes
nobs <- length(indobs)             # number of observed survival times

long_data <- data.frame(id = ID, y = y, obs_times = obs_times , X_total)


time_data <- data.frame(id=id_times,
                        times=times,
                        status=status,
                        X=X)



X1 <- as.matrix(X,ncol=1)
y <- long_data$y
n <- nrow(X1)
N <- length(y)
ID <- as.numeric(long_data$id)
obs_times <- long_data$obs_times


model <- cmdstan_model("code/long_model.stan")


mle <- model$optimize(data = list(y=y,N=N,n=n,X1=X1,ID=ID,obs_times=obs_times))
posterior_samples <- model$sample(data = list(y=y,N=N,n=n,X1=X1,ID=ID,obs_times=obs_times), chains = 1)

mle$summary()
posterior_samples$summary()