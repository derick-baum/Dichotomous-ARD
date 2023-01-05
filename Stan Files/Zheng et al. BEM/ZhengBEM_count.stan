data {             
  int N;
  int K;                    
  int y[N, K];                        
}

parameters {
  vector[N] alpha;                       
  vector[K] beta;                        
  vector<lower = 0 , upper = 1>[K] inv_omega;
  real mu_beta;
  real<lower=0> sigma_beta;                        
  real<lower=0> sigma_alpha;             
}
  
model {
  alpha ~ normal(0, sigma_alpha);  
  beta ~ normal(mu_beta, sigma_beta);                 
  for (k in 1:K) {
    real omega_k_m1;
    omega_k_m1 = inv(inv(inv_omega[k]) - 1);
    for (i in 1:N) {
      if (y[i,k] >= 0) {
        real xi_i_k;
        xi_i_k = omega_k_m1 * exp(alpha[i] + beta[k]);
        y[i,k] ~ neg_binomial(xi_i_k, omega_k_m1);
      }
    }
  }
}
