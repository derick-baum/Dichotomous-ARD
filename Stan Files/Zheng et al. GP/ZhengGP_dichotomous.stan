functions {
  real ll_a_lpdf(real alpha, real beta){
    return -exp(alpha + beta);
  }
  real ll_b_lpdf(real alpha, real beta){
    return log1m(exp(-exp(alpha + beta)));
  }
}

data {             
  int<lower=0> N;                        
  int<lower=0> K;                        
  int  y[N, K];
}
  
parameters {
  vector[N] alpha;                       
  vector[K] beta;                    
  real mu_beta;
  real<lower=0> sigma_beta;                        
  real<lower=0> sigma_alpha;              
}
  
model {
  alpha ~ normal(0, sigma_alpha);  
  for (i in 1:N) {
    for (k in 1:K) {
      if (y[i,k] == 0) {
        target += ll_a_lpdf(alpha[i]|beta[k]);
      } else if (y[i,k] == 1) {
        target += ll_b_lpdf(alpha[i]|beta[k]);
      }
    }
  }
}
