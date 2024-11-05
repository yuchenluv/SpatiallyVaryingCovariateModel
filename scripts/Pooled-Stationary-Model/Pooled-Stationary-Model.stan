// reference from R brms fun_gen_extreme_value.stan
// https://github.com/paul-buerkner/brms/blob/master/inst/chunks/fun_gen_extreme_value.stan
functions{
  real gev_lpdf(real y, real mu, real sigma, real xi){
    real x = (y - mu) / sigma;
     if (xi == 0) {
       return - log(sigma) - x - exp(-x);
     } 
     else {
       real t = 1 + xi * x;
       real inv_xi = 1 / xi;
       return - log(sigma) - (1 + inv_xi) * log(t) - pow(t, -inv_xi);
     }
  }
  real gev_lcdf(real y, real mu, real sigma, real xi) {
     real x = (y - mu) / sigma;
     if (xi == 0) {
       return - exp(-x);
     } else {
       return - pow(1 + xi * x, - 1 / xi);
     }
  }
  real gev_lccdf(real y, real mu, real sigma, real xi) {
     return log1m_exp(gev_lcdf(y | mu, sigma, xi));
  }
}

/* 
Referenced from https://mc-stan.org/docs/2_19/stan-users-guide/simulating-from-a-gaussian-process.html
*/

data {
  int<lower=1> n_station; //number of total stations
  int<lower=1> n_obs; //number of years
  vector[2] X[n_station]; //locations, longitude & latitude
  matrix[n_obs, n_station] y; //rainfall observations, indexed by [time, location]
}

parameters{
  real<lower=0> mu_rho; // kernel length parameter for mu
  real<lower=0> mu_alpha; // kernel std for the mu
  real<lower=0> logs_rho; // kernel length parameter for logs
  real<lower=0> logs_alpha; // kernel std for logs
  
  vector[n_station] mu;
  vector[n_station] logs;
  
  real<lower=0> xi; // GEV shape parameter, same for all locations
}

model{
  matrix[n_station, n_station] mu_K_xi; // kernel
  matrix[n_station, n_station] mu_K = gp_exponential_cov(X, mu_alpha, mu_rho);
  
  matrix[n_station, n_station] logs_K_xi; // kernel
  matrix[n_station, n_station] logs_K = gp_exponential_cov(X, logs_alpha, logs_rho);
  
  // diagonal elements
  for (n in 1:n_station)
    mu_K[n, n] = mu_K[n, n] + 10e-10;
  mu_K_xi = cholesky_decompose(mu_K);
  
  for (n in 1:n_station)
    logs_K[n, n] = logs_K[n, n] + 10e-10;
  logs_K_xi = cholesky_decompose(logs_K);
  
  xi ~ normal(0, 0.5);
  mu_rho ~ gamma(5, 1);
  mu_alpha ~ inv_gamma(5, 5);
  logs_rho ~ gamma(5, 1);
  logs_alpha ~ inv_gamma(5, 5);
  
  mu ~ multi_normal_cholesky(rep_vector(0, n_station), mu_K_xi);
  logs ~ multi_normal_cholesky(rep_vector(0, n_station), logs_K_xi);

  for (i in 1:n_obs){
    for (j in 1:n_station){
      if (y[i, j] > 0)
        y[i, j] ~ gev(mu[j], exp(logs[j]), xi);
    }
  }
}

