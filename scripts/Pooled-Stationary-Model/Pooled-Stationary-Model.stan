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
  int<lower=1> n_station;
  int<lower=1> n_obs;
  vector[2] X[n_station]; //locations, longitude & latitude
  matrix[n_obs, n_station] y; //observations of rainfall, indexed by [time, location]
}

parameters{
  // mean for multivariate normal distribution
  // real<lower=0> mu_m;
  // real logs_m;
  
  // co-regionalization parameters, coefficients for the latent GP
  real mu_w;
  real logs_w;
  
  real<lower=0> rho; // kernel length parameter for the latent GP
  real<lower=0> alpha; // kernel variance for the latent GP
  
  vector[n_station] mu;
  vector[n_station] logs;
  
  real<lower=0> xi; // GEV shape parameter, same for all locations
}

model{
  matrix[n_station, n_station] K_xi; // kernel
  matrix[n_station, n_station] K = gp_exponential_cov(X, alpha, rho);
  
  // diagonal elements
  for (n in 1:n_station)
    K[n, n] = K[n, n] + 10e-10;
    
  K_xi = cholesky_decompose(K);
  
  mu_w ~ std_normal();
  logs_w ~ std_normal();
  // mu_m ~ normal(0, 5);
  // logs_m ~ std_normal();
  xi ~ normal(0, 0.5);
  rho ~ gamma(5, 1);
  alpha ~ inv_gamma(5, 5);
  
  mu ~ multi_normal_cholesky(rep_vector(0, n_station), mu_w * K_xi);
  logs ~ multi_normal_cholesky(rep_vector(0, n_station), logs_w * K_xi);

  for (i in 1:n_obs){
    for (j in 1:n_station){
      if (y[i, j] > 0)
        y[i, j] ~ gev(mu[j], exp(logs[j]), xi);
    }
  }
}

