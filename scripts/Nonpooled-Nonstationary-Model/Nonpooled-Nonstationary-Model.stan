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
  int<lower=1> n_obs;
  vector[n_obs] y; // observations of rainfall for one station
  vector[n_obs] x; // covariates, currently just logCO2 anomalies
}

parameters{
  real<lower=0> mu0;
  real logs0; // use logarithm to ensure positive scale parameter
  
  real<lower=0> xi; // GEV shape parameter
  
  real mu_beta;// coefficients for mu
  real logs_beta;// coefficients for sigma
}

model{
  vector[n_obs] mu;
  vector[n_obs] sigma;
  
  mu0 ~ normal(5, 5);
  logs0 ~ std_normal();
  
  mu_beta ~ std_normal();
  logs_beta ~ std_normal();
  
  xi ~ normal(0, 0.5);
  
  for (i in 1:n_obs){
    mu[i] = mu0 + x[i] * mu_beta;
    sigma[i] = exp(logs0 + x[i] * logs_beta);
    if (y[i] > 0){
      y[i] ~ gev(mu[i], sigma[i], xi);
    }
  }
}
