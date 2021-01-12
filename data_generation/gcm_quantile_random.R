library(mvtnorm)

rm(list=ls())
nrep = 50

N = 500 # 100/200/300/500
T = 5
b_0 = c(6, 1.5)
sig_e2 = 0.5^2
rho = 0.3
Sig0 = diag(sqrt(c(0.5, 0.1)))%*%matrix(c(1,rho,rho,1),ncol=2)%*% diag(sqrt(c(0.5, 0.1)))
Lambda = cbind(rep(1,T), c(0:4))

for(i in 1:nrep){
  b_i = rmvnorm(N, mean = b_0, sigma = Sig0)
  e_i = rmvnorm(N, mean = rep(0,T), diag(sig_e2,T))
  sim_data = t(Lambda %*% t(b_i)) + e_i
  save(sim_data, b_0, sig_e2, rho, Sig0, Lambda, file=paste0("data/random/N", N, "/gcm_data_", i, ".rda"))
}

