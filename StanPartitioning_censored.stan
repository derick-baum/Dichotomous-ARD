functions {
  real ll_a_lpdf(real m, real rho, real d){
    return lbeta(m*(1/rho-1),d+(1-m)*(1/rho-1))-lbeta(m*(1/rho-1),(1-m)*(1/rho-1));
  }
  real ll_b_lpdf(real m, real rho, real d){
    return log(d)+lbeta(m*(1/rho-1)+1,d+(1-m)*(1/rho-1)-1)-lbeta(m*(1/rho-1),(1-m)*(1/rho-1));
  }
  real ll_c_lpdf(real m, real rho, real d){
    return log1m(exp(lbeta(m*(1/rho-1),d+(1-m)*(1/rho-1))-lbeta(m*(1/rho-1),(1-m)*(1/rho-1)))+exp(log(d)+lbeta(m*(1/rho-1)+1,d+(1-m)*(1/rho-1)-1)-lbeta(m*(1/rho-1),(1-m)*(1/rho-1))));
  }
}


data {
  int K; 
  int N;
  int x[K,N];
  real <lower=0,upper=1> m[K];
  vector[N] L;
}
parameters {
  real <lower=1> d[N];
  real <lower=3,upper=8> mu;
  real <lower=1./4,upper=2> sigma;
  real <lower=0,upper=1> rho[K];
}


model {
  //Prior
  d ~ lognormal(mu,sigma);
  mu ~ uniform(3,8);
  sigma ~ uniform(1./4,2);
  rho ~ uniform(0,1);
  //Likelihood
  for (i in 1:N){
    for (k in 1:K) {
      if (x[k,i]==0) {
        target += ll_a_lpdf(m[k]|rho[k],d[i]);
      } else if (x[k,i]==1){
        target += ll_b_lpdf(m[k]|rho[k],d[i]);
      } else if (x[k,i]==2){
        target += ll_c_lpdf(m[k]|rho[k],d[i]);
      }
    }
  }
}

