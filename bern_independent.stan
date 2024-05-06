// Simultaneous modeling of three independent response diseases,
// using Bernoulli distribution.

functions {
  // Function for computing Q of QR decomposition
  vector create_Q(int N) {
    vector[2 * N] Q;
    for (i in 1 : N) {
      Q[i] = -sqrt((N - i) / (N - i + 1.0));
      Q[i + N] = inv_sqrt((N - i) * (N - i + 1));
    }
    return Q;
  }
  // sum-to-zero constrained x from x_raw
  vector sum_to_zero(vector x_raw, vector Q, int N) {
    vector[N] x;
    real x_aux = 0;
    for (i in 1 : (N - 1)) {
      x[i] = x_aux + x_raw[i] * Q[i];
      x_aux = x_aux + x_raw[i] * Q[i + N];
    }
    x[N] = x_aux;
    return x;
  }
}
data {
  int<lower=0> T; // number of time points
  int<lower=0> N; // number of regions
  int<lower=0> n; // number of observations
  // response variable
  array[n] int<lower=0> y_p;
  array[n] int<lower=0> y_m;
  array[n] int<lower=0> y_s;
  // predictors
  int<lower=0> K_x; // number of predictors based on site's own characteristics
  int<lower=0> K_z; // number of predictors based on neighbors' characteristics
  matrix[n, K_x] x_p; // predictors related to focal site
  matrix[n, K_x] x_m; // predictors related to focal site
  matrix[n, K_x] x_s; // predictors related to focal site
  matrix[n, K_z] z_p; // predictors related to neighbors
  matrix[n, K_z] z_m; // predictors related to neighbors
  matrix[n, K_z] z_s; // predictors related to neighbors
  // indices for observed cells, same for all diseases in this case
  array[n] int time_id;
  array[n] int month_id;
  array[n] int region_id;
}
transformed data {
  // Scale parameter for the sum constrained prior
  real scaleN = inv(sqrt(1 - inv(N)));
  real scale12 = inv(sqrt(1 - inv(12)));
  // Create QR decomposition for sum-to-zero constraint
  vector[2 * N] QN = create_Q(N);
  vector[2 * 12] Q12 = create_Q(12);
}
parameters {
  // Population-level parameters
  vector[K_x] beta_p;
  vector[K_x] beta_m;
  vector[K_x] beta_s;
  vector[K_z] gamma_p;
  vector[K_z] gamma_m;
  vector[K_z] gamma_s;

  // Raw parameters for seasonal effects with sum-to-zero constraint
  // In case of omitting the seasonal term, comment out these three rows
  vector[11] seasonal_p_raw;
  vector[11] seasonal_m_raw;
  vector[11] seasonal_s_raw;

  // Raw parameters for factor loadings lambda
  vector[N - 1] lambda_p_raw;
  vector[N - 1] lambda_m_raw;
  vector[N - 1] lambda_s_raw;

  // Standard deviations of factor loadings lambda
  real<lower=0> sigma_lambda_p;
  real<lower=0> sigma_lambda_m;
  real<lower=0> sigma_lambda_s;

  // Latent factors tau for each three responses (p, m, s)
  array[T - 1] vector[3] tau;
  // Standard deviations of the latent factor noise terms
  vector<lower=0>[3] sigma_tau;

  // raw site-specific intercepts and coefficients
  vector[N] a_p_raw;
  vector[N] a_m_raw;
  vector[N] a_s_raw;
  vector[N] b_p_raw;
  vector[N] b_m_raw;
  vector[N] b_s_raw;
  vector[N] c_p_raw;
  vector[N] c_m_raw;
  vector[N] c_s_raw;

  // standard deviations of site-specific intercepts and coefficients
  real<lower=0> sigma_a_p;
  real<lower=0> sigma_a_m;
  real<lower=0> sigma_a_s;
  real<lower=0> sigma_b_p;
  real<lower=0> sigma_b_m;
  real<lower=0> sigma_b_s;
  real<lower=0> sigma_c_p;
  real<lower=0> sigma_c_m;
  real<lower=0> sigma_c_s;
}
transformed parameters {
  // set sum-to-zero constraints for seasonal effects
  vector[12] seasonal_p = sum_to_zero(seasonal_p_raw, Q12, 12);
  vector[12] seasonal_m = sum_to_zero(seasonal_m_raw, Q12, 12);
  vector[12] seasonal_s = sum_to_zero(seasonal_s_raw, Q12, 12);
  
  // fix means of lambda to 1
  vector[N] lambda_p = 1 + sigma_lambda_p * sum_to_zero(lambda_p_raw, QN, N);
  vector[N] lambda_m = 1 + sigma_lambda_m * sum_to_zero(lambda_m_raw, QN, N);
  vector[N] lambda_s = 1 + sigma_lambda_s * sum_to_zero(lambda_s_raw, QN, N);

  // transformations of a, b and c
  vector[N] a_p = sigma_a_p * a_p_raw; // implies a_p ~ normal(0, sigma_a_p)
  vector[N] a_m = sigma_a_m * a_m_raw;
  vector[N] a_s = sigma_a_s * a_s_raw;
  vector[N] b_p = 1 + sigma_b_p * b_p_raw; // implies b_p ~ normal(1, sigma_b_p)
  vector[N] b_s = 1 + sigma_b_s * b_s_raw;
  vector[N] b_m = 1 + sigma_b_m * b_m_raw;
  vector[N] c_p = 1 + sigma_c_p * c_p_raw; // implies c_p ~ normal(1, sigma_c_p)
  vector[N] c_s = 1 + sigma_c_s * c_s_raw;
  vector[N] c_m = 1 + sigma_c_m * c_m_raw;
}
model {
  // set the priors (mostly noninformative)
  
  // seasonal terms and lambdas need specific deviations due to the sum-to-zero constraints
  seasonal_p_raw ~ normal(0, scale12);
  seasonal_m_raw ~ normal(0, scale12);
  seasonal_s_raw ~ normal(0, scale12);
  lambda_p_raw ~ normal(0, scaleN);
  lambda_m_raw ~ normal(0, scaleN);
  lambda_s_raw ~ normal(0, scaleN);

  beta_p ~ normal(0, 2);
  beta_m ~ normal(0, 2);
  beta_s ~ normal(0, 2);
  gamma_p ~ normal(0, 2);
  gamma_m ~ normal(0, 2);
  gamma_s ~ normal(0, 2);

  sigma_lambda_p ~ gamma(2, 1);
  sigma_lambda_m ~ gamma(2, 1);
  sigma_lambda_s ~ gamma(2, 1);
  
  // non-correlated random walks for tau
  sigma_tau ~ gamma(2, 1);
  tau[1] ~ normal(-2, 2);
  for(t in 2:(T-1)) {
    tau[t] ~ normal(tau[t - 1], sigma_tau);
  }

  sigma_a_p ~ gamma(2, 1);
  sigma_a_m ~ gamma(2, 1);
  sigma_a_s ~ gamma(2, 1);

  sigma_b_p ~ gamma(2, 1);
  sigma_b_m ~ gamma(2, 1);
  sigma_b_s ~ gamma(2, 1);

  sigma_c_p ~ gamma(2, 1);
  sigma_c_m ~ gamma(2, 1);
  sigma_c_s ~ gamma(2, 1);

  // give standard normal priors here for the transformations
  a_p_raw ~ std_normal();
  a_m_raw ~ std_normal();
  a_s_raw ~ std_normal();
  b_p_raw ~ std_normal();
  b_m_raw ~ std_normal();
  b_s_raw ~ std_normal();
  c_p_raw ~ std_normal();
  c_m_raw ~ std_normal();
  c_s_raw ~ std_normal();
  
  // combine all the parts to form the complete model
  {
    vector[n] eta_p = lambda_p[region_id] .* to_vector(tau[time_id, 1]) + 
      seasonal_p[month_id] + 
      a_p[region_id] + 
      b_p[region_id] .* (x_p * beta_p) + 
      c_p[region_id] .* (z_p * gamma_p);
    vector[n] eta_m = lambda_m[region_id] .* to_vector(tau[time_id, 2]) + 
      seasonal_m[month_id] + 
      a_m[region_id] + 
      b_m[region_id] .* (x_m * beta_m) + 
      c_m[region_id] .* (z_m * gamma_m);
    vector[n] eta_s = lambda_s[region_id] .* to_vector(tau[time_id, 3]) + 
      seasonal_s[month_id] + 
      a_s[region_id] + 
      b_s[region_id] .* (x_s * beta_s) + 
      c_s[region_id] .* (z_s * gamma_s);
    y_p ~ bernoulli_logit(eta_p);
    y_m ~ bernoulli_logit(eta_m);
    y_s ~ bernoulli_logit(eta_s);
  }
}
//generated quantities {
  // get the log likelihood values
  // due to computational and memory issues we run separate generated
  // quantities block for this
//  vector[n] log_lik;
//  { 
//    vector[n] eta_p = lambda_p[region_id] .* to_vector(tau[time_id, 1]) + 
//      seasonal_p[month_id] + 
//      a_p[region_id] + 
//      b_p[region_id] .* (x_p * beta_p) + 
//      c_p[region_id] .* (z_p * gamma_p);
//    vector[n] eta_m = lambda_m[region_id] .* to_vector(tau[time_id, 2]) + 
//      seasonal_m[month_id] + 
//      a_m[region_id] + 
//      b_m[region_id] .* (x_m * beta_m) + 
//      c_m[region_id] .* (z_m * gamma_m);
//    vector[n] eta_s = lambda_s[region_id] .* to_vector(tau[time_id, 3]) + 
//      seasonal_s[month_id] + 
//      a_s[region_id] + 
//      b_s[region_id] .* (x_s * beta_s) + 
//      c_s[region_id] .* (z_s * gamma_s);
//      
//    for(i in 1:n) {
//      log_lik[i] = bernoulli_logit_lpmf(y_p[i] | eta_p[i]) +
//                   bernoulli_logit_lpmf(y_m[i] | eta_m[i]) +
//                   bernoulli_logit_lpmf(y_s[i] | eta_s[i]);
//    }
//  }
//}
