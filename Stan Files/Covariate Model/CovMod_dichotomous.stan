functions {
  real ll_a_lpdf(real alpha, real beta, real gamma_age, real gamma_sex, real age, real male){
    return -exp(alpha + beta + gamma_age * age + gamma_sex * male);
  }
  real ll_b_lpdf(real alpha, real beta, real gamma_age, real gamma_sex, real age, real male){
    return log1m(exp(-exp(alpha + beta + gamma_age * age + gamma_sex * male)));
  }
}

data {             
  int N;
  int K;                                                
  matrix[N, K]  y;
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
// priors
  alpha ~ normal(0, sigma_alpha);  
  beta ~ normal(mu_beta, sigma_beta);
  gamma_age ~ normal(0, sigma_gamma);
  gamma_sex ~ normal(0, sigma_gamma);
  sigma_gamma ~ gamma(1, 0.01);
  for (i in 1:N) {
    for (k in 1:K) {
      if (y[i, k] == 0) {
        target += ll_a_lpdf(alpha[i]|beta[k], gamma_age[k], gamma_sex[k], age[i], male[i]);
      } else if (y[i,k] == 1) {
        target += ll_b_lpdf(alpha[i]|beta[k], gamma_age[k], gamma_sex[k], age[i], male[i]);
      }
    }
  }
}
