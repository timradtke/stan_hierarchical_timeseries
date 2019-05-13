data {
  int N;
  real yoy[N];
  real yoy_1[N];
  int N_ids;
  int id[N];
}
parameters {
  real mu;
  vector[N_ids] mu_i;
  vector<lower=0,upper=1>[N_ids] alpha_i;
  real<lower=0> sigma_mu_i;
  vector<lower=0>[N_ids] sigma_i;
  real<lower=0> sigma;
  real<lower=0> theta_1;
  real<lower=0> theta_2;
}
transformed parameters {
  vector[N] mu_t;
  for (i in 1:N) {
    mu_t[i] = (1 - alpha_i[id[i]]) * (mu + mu_i[id[i]]) + alpha_i[id[i]] * yoy_1[i];
  }
}
model {
  
  target += normal_lpdf( mu | 0, 0.05);
  target += normal_lpdf( mu_i | 0, sigma_mu_i);
  target += cauchy_lpdf( sigma_mu_i | 0, 0.05);
  target += beta_lpdf( alpha_i | theta_1, theta_2);
  target += exponential_lpdf( theta_1 | 1.0/16.0);
  target += exponential_lpdf( theta_2 | 1.0/32.0);
  target += cauchy_lpdf( sigma_i | 0, sigma);
  target += cauchy_lpdf( sigma | 0, 0.05);
  
  for (i in 1:N) {
    target += normal_lpdf( yoy[i] | mu_t[i], sigma_i[id[i]]);
  }
}
generated quantities {
  vector[N] yoy_hat;
  for (i in 1:N) {
    yoy_hat[i] = normal_rng(mu_t[i], sigma_i[id[i]]);
  }
}
