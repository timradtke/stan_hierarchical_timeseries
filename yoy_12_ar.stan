data {
  int N;
  real yoy[N];
  real yoy_1[N];
  real yoy_12[N];
  int yoy_12_available[N];
  int N_stations;
  int station[N];
}
parameters {
  real a;
  //real b;
  real c;
  real mu_b;
  real<lower=0> sigma;
  real<lower=0> sigma_station;
  real<lower=0> sigma_b;
  vector[N_stations] a_station;
  vector[N_stations] b;
}
model {
  vector[N] mu;
  for (i in 1:N) {
    if (yoy_12_available[i] == 0)
      mu[i] = a + a_station[station[i]] + b[station[i]] * yoy_1[i];
    else
      mu[i] = a + a_station[station[i]] + b[station[i]] * yoy_1[i] + c * yoy_12[i];
    }
  target += normal_lpdf( a | 0, 0.05);
  target += normal_lpdf( b | mu_b, sigma_b);
  target += normal_lpdf( c | 0, 1);
  target += normal_lpdf( a_station | 0, sigma_station);
  target += cauchy_lpdf( sigma | 0, 0.25);
  target += cauchy_lpdf( sigma_station | 0, 0.05);
  target += normal_lpdf( mu_b | 0, 0.5);
  target += cauchy_lpdf( sigma_b | 0, 0.25);
  target += normal_lpdf( yoy | mu, sigma);
}
generated quantities {
  vector[N] yoy_hat;
  {
    vector[N] mu;
    for (i in 1:N) {
      mu[i] = a + a_station[station[i]] + b[station[i]] * yoy_1[i] + c * yoy_12[i];
      yoy_hat[i] = normal_rng(mu[i], sigma);
    }
  }
}
