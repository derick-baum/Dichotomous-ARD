functions {
  real ll_a_lpdf(real xi, real inv_omega){
    return xi*log(inv_omega);
  }
  real ll_b_lpdf(real xi, real inv_omega){
    return log(xi) + xi * log(inv_omega) + log1m(inv_omega);
  }
  real ll_c_lpdf(real xi, real inv_omega){
    return log1m(exp(xi * log(inv_omega)) + exp(log(xi) + xi * log(inv_omega) + log1m(inv_omega)));
  }
}

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
      real xi_i_k;
      xi_i_k = omega_k_m1 * exp(alpha[i] + beta[k]);
      if (y[i, k] == 0) {
        target += ll_a_lpdf(xi_i_k|inv_omega[k]);
      } else if (y[i, k] == 1) {
        target += ll_b_lpdf(xi_i_k|inv_omega[k]);
      } else if (y[i, k] == 2) {
        target += ll_c_lpdf(xi_i_k|inv_omega[k]);
      }
    }
  }
}
