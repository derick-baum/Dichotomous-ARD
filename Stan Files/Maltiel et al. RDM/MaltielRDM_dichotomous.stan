functions {
  real ll_a_lpdf(real m, real d){
    return d * log1m(m);
  }
  real ll_b_lpdf(real m, real d){
    return log1m(pow((1 - m), d));
  }
}

data {
  int N;
  int K; 
  matrix[N, K] y;
  vector<lower=0, upper=1>[K] m;
}

parameters {
  vector <lower=0>[N] d;
  real <lower=3,upper=8> mu;
  real <lower=1./4,upper=2> sigma;
}

model {
  d ~ lognormal(mu,sigma);
  mu ~ uniform(3,8);
  sigma ~ uniform(1./4,2);
  for (i in 1:N){
    for (k in 1:K) {
      if (y[i, k] == 0) {
        target += ll_a_lpdf(m[k]|d[i]);
      } else if (y[i, k] == 1) {
        target += ll_b_lpdf(m[k]|d[i]);
      }
    }
  }
}
