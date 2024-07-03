
data {
  int<lower=0> n1;
  int incidence [n1];
  vector[n1] expect;
}

parameters {
  vector<lower=0.000001>[n1] Re;
}

model {
  Re ~ normal(0, 25);
  vector [n1] expect_2;
  for (m in 22:n1) {
  for (i in 0:21) {
      expect_2[m-i] = Re[m]*expect[m-i];
      incidence[m-i] ~ poisson(expect_2[m - i]);
    }
  }
}
