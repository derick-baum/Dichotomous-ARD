functions {
  real ll_a_lpdf(real m, real rho, real d){
    return lbeta(m * (1 / rho - 1), d + (1 - m) * (1 / rho - 1)) - lbeta(m * (1 / rho - 1), (1 - m) * (1 / rho - 1));
  }
  real ll_b_lpdf(real m, real rho, real d){
    return log1m(exp(lbeta(m * (1 / rho - 1), d + (1 - m) * (1 / rho - 1)) - lbeta(m * (1 / rho - 1), (1 - m) * (1 / rho - 1))));
  }
}

data {
  int K; 
  int N;
  matrix[N, K] y;
  vector <lower=0, upper=1>[K] m;
}

parameters {
  vector <lower=0>[N] d;
  real <lower=3, upper=8> mu;
  real <lower=1./4, upper=2> sigma;
  vector <lower=0, upper=1>[K] rho;
}

model {
  d ~ lognormal(mu, sigma);
  mu ~ uniform(3, 8);
  sigma ~ uniform(1./4, 2);
  rho ~ uniform(0, 1);
  for (i in 1:N){
    for (k in 1:K) {
      if (y[i, k] == 0) {
        target += ll_a_lpdf(m[k]|rho[k], d[i]);
      } else if (y[i, k] == 1) {
        target += ll_b_lpdf(m[k]|rho[k], d[i]);
      }
    }
  }
}
