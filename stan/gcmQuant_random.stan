data{
  int<lower=0> N;
  int<lower=0> T;
  vector[N*T] X;
  vector[N*T] y;
  real<lower = 0, upper = 1> tau1;
  real<lower = 0, upper = 1> tau2;
  real<lower = 0> omega_prior;
  vector[2] beta_0;
  cov_matrix[2] Sigma_beta0;
}

transformed data{
  vector[2] xi_vector;
  vector[2] sig_vector;
  matrix[2,2] L_tilde;
  int len;

  xi_vector[1] = (1-2*tau1)/(tau1*(1-tau1));
  xi_vector[2] = (1-2*tau2)/(tau2*(1-tau2));
    
  sig_vector[1] = sqrt(2/(tau1*(1-tau1)));
  sig_vector[2] = sqrt(2/(tau2*(1-tau2)));

  L_tilde = diag_matrix(sig_vector);
  
  len = N*T;
}

parameters{
  real<lower=0> sigma;
  real<lower=0> delta1;
  real<lower=0> delta2;
  vector<lower=0>[N] v;
  vector[2] beta;
  matrix[N,2] b;
  corr_matrix[2] Omega; 
  //cov_matrix[2] Sigma_beta;
}


transformed parameters{
  vector[len] mu;
  vector[2] mean_b[N];
  matrix[2,2] matrix_b[N];

  matrix[len,2] b_span;
  matrix[2,2] D;
  matrix[2,2] Sigma_b;
  
  vector[2] delta_vector; 
  real<lower=0> sig_y;
  
  for(i in 1:len){
    b_span[i,] = b[(i+T-1)/T,];
  }
  
  mu = (beta[1] + b_span[,1]) + (beta[2] + b_span[,2]) .* X;
  sig_y = sqrt(sigma);
  
  delta_vector[1] = delta1;
  delta_vector[2] = delta2;
    
  D = diag_matrix(delta_vector);
  Sigma_b = L_tilde*Omega*L_tilde;
  
  for(i in 1:N){
    mean_b[i] = D*xi_vector*v[i];
    matrix_b[i] = v[i]*D*Sigma_b*D';
  }
}

model{
  // model
  y ~ normal(mu, sig_y);
  
  // data augmentation
  v ~ exponential(1);
  
  // priors
  sigma ~ cauchy(0, 1);
  delta1 ~ cauchy(0, 1);
  delta2 ~ cauchy(0, 1);
  
  beta ~ multi_normal(beta_0, Sigma_beta0);
  
  for(i in 1:N){
    b[i,] ~ multi_normal(mean_b[i],matrix_b[i]);
  }
  
  Omega ~ lkj_corr(omega_prior); // LKJ prior on the correlation matrix 
}

generated quantities{
  vector[len] pp_y;
  for (i in 1:len) {
   pp_y[i] = normal_rng(mu[i], sig_y);
  }
}
