//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> n1;
  int<lower=0> n2;
  vector [n2] wbe;
  int row [n2];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[n1] mu;
  real<lower=0> s1;
  real<lower=0> s2;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  int x [n2];
  mu ~ normal(0, 10);
  s1 ~ normal(0, 10);
  s2 ~ normal(0, 10);
  for(i in 3:n1){
  mu[i] ~ normal(2*mu[i-1] - mu[i-2], s1);
}
  for(k in 1:n2){
    x[k] = row[k];
  wbe[k] ~ normal(mu[x[k]], s2);
}
}






