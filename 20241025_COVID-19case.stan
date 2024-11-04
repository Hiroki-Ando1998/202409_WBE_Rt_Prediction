
data {
  int<lower=0> n1;
  int incidence [n1];
}

parameters {
  vector<lower=0.01, upper =3>[n1] Re;
}

transformed parameters {
  // Generation interval (g2とは向きが違う)
  vector[10] g1;
  for (i in 1:10) {
      g1[i] = gamma_cdf(i, 4, 1) - gamma_cdf(i-1, 4, 1);  // 差分を計算
    }
  }
  

model {
  Re ~ normal(0, 3);
  vector [n1] expect;
  vector [n1] expect_2;
  
    // Newly infected people
  for (k in 11:n1) {
    // 条件に基づいてi1の値を準備
    expect[k] = g1[1]*incidence[k-1]+g1[2]*incidence[k-2]+g1[3]*incidence[k-3]+g1[4]*incidence[k-4]+g1[5]*incidence[k-5]+g1[6]*incidence[k-6]+
                g1[7]*incidence[k-7]+g1[8]*incidence[k-8]+g1[9]*incidence[k-9]+g1[10]*incidence[k-10];
      }

  for (m in 18:n1) {
  for (i in 0:6) {
      expect_2[m-i] = Re[m]*expect[m-i];
      incidence[m-i] ~ poisson(expect_2[m - i]);
    }
  }
}




