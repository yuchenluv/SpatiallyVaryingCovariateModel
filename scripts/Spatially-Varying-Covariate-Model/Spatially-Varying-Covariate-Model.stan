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
  int<lower=1> n_station; // number of stations
  int<lower=1> n_obs; // number of years
  vector[2] X[n_station]; // locations, longitude & latitude
  matrix[n_obs, n_station] y; // observations of rainfall, indexed by [time, location]
  vector[n_obs] x; // covariates, currently just using logCO2 anomalies
}

parameters{
  // kernel parameters
  real<lower=0> mu_rho; // kernel length scale
  real<lower=0> mu_alpha; // kernel std
  real<lower=0> logs_rho; // kernel length scale
  real<lower=0> logs_alpha; // kernel std
  real<lower=0> mu0_rho; // kernel length scale
  real<lower=0> mu0_alpha; // kernel std
  real<lower=0> logs0_rho; // kernel length scale
  real<lower=0> logs0_alpha; // kernel std
  
  real<lower=0> xi; // GEV shape parameter, same for all locations
  
  vector[n_station] mu_beta; // coefficients for mu
  vector[n_station] logs_beta; // coefficients for log(sigma)
  vector[n_station] mu0; // intercept for mu
  vector[n_station] logs0; // intercept for log(sigma)
}

model{
  // gp_exponential_cov from https://mc-stan.org/math/group__opencl_gac47f7d0f97a1035236ac74de3b965020.html
  matrix[n_station, n_station] mu_K_xi; // kernel for mu_beta
  matrix[n_station, n_station] mu_K = gp_exponential_cov(X, mu_alpha, mu_rho);
  
  matrix[n_station, n_station] logs_K_xi; // kernel for logs_beta
  matrix[n_station, n_station] logs_K = gp_exponential_cov(X, logs_alpha, logs_rho);
  
  matrix[n_station, n_station] mu0_K_xi; // kernel for mu0
  matrix[n_station, n_station] mu0_K = gp_exponential_cov(X, mu0_alpha, mu0_rho);
  
  matrix[n_station, n_station] logs0_K_xi; // kernel for logs0
  matrix[n_station, n_station] logs0_K = gp_exponential_cov(X, logs0_alpha, logs0_rho);
  
  matrix[n_obs, n_station] mu;
  matrix[n_obs, n_station] sigma;
  
  // diagonal elements
  for (n in 1:n_station)
    mu_K[n, n] = mu_K[n, n] + 10e-6;
  mu_K_xi = cholesky_decompose(mu_K);
  
  for (n in 1:n_station)
    logs_K[n, n] = logs_K[n, n] + 10e-6;
  logs_K_xi = cholesky_decompose(logs_K);
  
  for (n in 1:n_station)
    mu0_K[n, n] = mu0_K[n, n] + 10e-6;
  mu0_K_xi = cholesky_decompose(mu0_K);
  
  for (n in 1:n_station)
    logs0_K[n, n] = logs0_K[n, n] + 10e-6;
  logs0_K_xi = cholesky_decompose(logs0_K);
  
  mu_rho ~ gamma(5, 1);
  mu_alpha ~ inv_gamma(5, 5);
  logs_rho ~ gamma(5, 1);
  logs_alpha ~ inv_gamma(5, 5);
  mu0_rho ~ gamma(5, 1);
  mu0_alpha ~ inv_gamma(5, 5);
  logs0_rho ~ gamma(5, 1);
  logs0_alpha ~ inv_gamma(5, 5);
  
  xi ~ normal(0, 0.5);
  
  mu_beta ~ multi_normal_cholesky(rep_vector(0, n_station), mu_K_xi);
  mu0 ~ multi_normal_cholesky(rep_vector(0, n_station), mu0_K_xi);
  logs_beta ~ multi_normal_cholesky(rep_vector(0, n_station), logs_K_xi);
  logs0 ~ multi_normal_cholesky(rep_vector(0, n_station), logs0_K_xi);

  for (i in 1:n_obs){
    for (j in 1:n_station){
      mu[i, j] = mu0[j] + x[i] * mu_beta[j];
      sigma[i, j] = exp(logs0[j] + x[i] * logs_beta[j]);
      if (y[i, j] > 0){ // filter out missing values
        y[i, j] ~ gev(mu[i, j], sigma[i, j], xi);
      }
    }
  }
}
