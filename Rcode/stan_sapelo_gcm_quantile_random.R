library(rstan)

cmdArgs <- commandArgs(trailingOnly = TRUE)
numCores <- as.integer(cmdArgs[1])
options(mc.cores=numCores)

filename1 <- cmdArgs[2]
filename2 <- cmdArgs[3]

tau1 <- as.numeric(cmdArgs[4])
tau2 <- as.numeric(cmdArgs[5])

load(file=filename1)

y<-sim_data

## prepare intitial values&fixed values
N<-length(y)/5
X<-rep(c(0:4),N)
X<-as.vector(X)

## use linear regression as intial value for MCMC 
y<-as.vector(t(y))
lm_est<-lm(y~X)

beta_0<-c(lm_est$coefficients[1],lm_est$coefficients[2])
Sigma_beta0<-matrix(100*c(1,0.0,0.0,0.1),ncol = 2)


Niter = 2000
max_treedepth <-  10
adapt_delta <- 0.99
n_chains = 4

dat<-list(N=N, T=5, y = y, X = X, tau1 = tau1, tau2 = tau2, 
           beta_0 = beta_0, Sigma_beta0 = Sigma_beta0, omega_prior = 3)
fit<-stan(file = "../stan/gcmQuant_random.stan", model_name = "GCM_Quantile", 
           pars = c('beta', 'delta1', 'delta2', 'sigma', 'Omega'),
           control = list(max_treedepth = max_treedepth,
                          adapt_delta = adapt_delta),
           data = dat, iter = Niter, chains = n_chains)

save(fit, file=filename2)




