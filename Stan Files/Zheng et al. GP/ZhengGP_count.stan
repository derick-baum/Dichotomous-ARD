data {             
  int N;
  int K;                    
  int y[N, K];
}
  
parameters {
  vector[N] alpha;                       
  vector[K] beta;                    
  real mu_beta;
  real <lower=0> sigma_beta;                        
  real <lower=0> sigma_alpha;              
}
  
model {
  alpha ~ normal(0, sigma_alpha);  
  beta ~ normal(mu_beta, sigma_beta);  
  for (i in 1:N) {
    for (k in 1:K) {
      if (y[i, k] >= 0) {
        real rate = exp(alpha[i] + beta[k]);
        y[i, k] ~ poisson(rate);
      }
    }
  }
}
