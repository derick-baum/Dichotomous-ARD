data {             
  int N;
  int K;                        
  int y[N, K];
  vector[N] age;
  vector[N] male;
}
  
parameters {
  vector[N] alpha;                       
  vector[K] beta;                    
  real mu_beta;
  real<lower=0> sigma_beta;                        
  real<lower=0> sigma_alpha;
  vector[K] gamma_age; 
  vector[K] gamma_sex; 
  real<lower=0> sigma_gamma; 
}
  
model {
  alpha ~ normal(0, sigma_alpha);  
  beta ~ normal(mu_beta, sigma_beta); 
  gamma_age ~ normal(0, sigma_gamma);
  gamma_sex ~ normal(0, sigma_gamma);
  sigma_gamma ~ gamma(1, 0.01);
  for (i in 1:N) {
    for (k in 1:K) {
      if (y[i, k] >= 0) {
        real rate = exp(alpha[i] + beta[k] + gamma_age[k] * age[i] + gamma_sex[k] * male[i]);
        y[i, k] ~ poisson(rate);
      }
    }
  }
}
