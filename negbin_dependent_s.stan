// Simultaneous modeling of three dependent response diseases without 
// seasonal effect using negative binomial distribution.

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
  // Variables below are here same for all responses
  // In general we could have different predictors and missing data patters
  int<lower=0> K_x; // number of predictors based on site's own characteristics
  int<lower=0> K_z; // number of predictors based on neighbors' characteristics
  matrix[n, K_x] x; // predictors related to focal site
  matrix[n, K_z] z; // predictors related to neighbors
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
  // Cholesky factor of the correlation matrix of latent factors
  cholesky_factor_corr[3] L;
  
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
  
  // dispersion parameters
  real alpha_phi_p;
  real alpha_phi_m;
  real alpha_phi_s;
  
  vector[N] phi_p_raw;
  vector[N] phi_m_raw;
  vector[N] phi_s_raw;

  real<lower=0> sigma_phi_p;
  real<lower=0> sigma_phi_m;
  real<lower=0> sigma_phi_s;
}
transformed parameters {
  // fix means of lambda to 1
  vector[N] lambda_p = 1 + sigma_lambda_p * sum_to_zero(lambda_p_raw, QN, N);
  vector[N] lambda_m = 1 + sigma_lambda_m * sum_to_zero(lambda_m_raw, QN, N);
  vector[N] lambda_s = 1 + sigma_lambda_s * sum_to_zero(lambda_s_raw, QN, N);
  
  // transformations of a, b and c
  vector[N] a_p = sigma_a_p * a_p_raw; // implies a_p ~ normal(0, sigma_a_p)
  vector[N] a_m = sigma_a_m * a_m_raw; // implies a_p ~ normal(0, sigma_a_p)
  vector[N] a_s = sigma_a_s * a_s_raw; // implies a_p ~ normal(0, sigma_a_p)
  vector[N] b_p = 1 + sigma_b_p * b_p_raw; // implies b_p ~ normal(1, sigma_b_p)
  vector[N] b_m = 1 + sigma_b_m * b_m_raw; // implies b_p ~ normal(1, sigma_b_p)
  vector[N] b_s = 1 + sigma_b_s * b_s_raw; // implies b_p ~ normal(1, sigma_b_p)
  vector[N] c_p = 1 + sigma_c_p * c_p_raw; // implies c_p ~ normal(1, sigma_c_p)
  vector[N] c_m = 1 + sigma_c_m * c_m_raw; // implies c_p ~ normal(1, sigma_c_p)
  vector[N] c_s = 1 + sigma_c_s * c_s_raw; // implies c_p ~ normal(1, sigma_c_p)
  
  // transformations phi
  vector[N] phi_p = sigma_phi_p * phi_p_raw;
  vector[N] phi_m = sigma_phi_m * phi_m_raw;
  vector[N] phi_s = sigma_phi_s * phi_s_raw;
}
model {
  // set the priors (mostly noninformative)
  
  // lambdas need specific deviations due to the sum-to-zero constraints
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
  
  // intercorrelated random walk for tau
  sigma_tau ~ gamma(2, 1);
  L ~ lkj_corr_cholesky(1);
  tau[1] ~ normal(-2, 2);
  tau[2:(T - 1)] ~ multi_normal_cholesky(tau[1:(T - 2)], diag_pre_multiply(sigma_tau, L));
  
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
  
  alpha_phi_p ~ std_normal();
  alpha_phi_m ~ std_normal();
  alpha_phi_s ~ std_normal();
  
  phi_p_raw ~ std_normal();
  phi_m_raw ~ std_normal();
  phi_s_raw ~ std_normal();
  
  sigma_phi_p ~ gamma(2, 1);
  sigma_phi_m ~ gamma(2, 1);
  sigma_phi_s ~ gamma(2, 1);
  // combine all the parts to form the complete model
  {
    vector[n] eta_p = lambda_p[region_id] .* to_vector(tau[time_id, 1]) + 
      a_p[region_id] + 
      b_p[region_id] .* (x * beta_p) + 
      c_p[region_id] .* (z * gamma_p);
    vector[n] eta_m = lambda_m[region_id] .* to_vector(tau[time_id, 2]) + 
      a_m[region_id] + 
      b_m[region_id] .* (x * beta_m) + 
      c_m[region_id] .* (z * gamma_m);
    vector[n] eta_s = lambda_s[region_id] .* to_vector(tau[time_id, 3]) + 
      a_s[region_id] + 
      b_s[region_id] .* (x * beta_s) + 
      c_s[region_id] .* (z * gamma_s);
    y_p ~ neg_binomial_2_log(eta_p, exp(alpha_phi_p + phi_p[region_id]));
    y_m ~ neg_binomial_2_log(eta_m, exp(alpha_phi_m + phi_m[region_id]));
    y_s ~ neg_binomial_2_log(eta_s, exp(alpha_phi_s + phi_s[region_id]));
  }
}
generated quantities {
  matrix[3, 3] corr_tau = multiply_lower_tri_self_transpose(L);
  
  // get the log likelihood values
  // due to computational and memory issues we run separate generated
  // quantities block for this and comment corr_tau out then
//  vector[n] log_lik;
//  { 
//    vector[n] eta_p = lambda_p[region_id] .* to_vector(tau[time_id, 1]) + 
//      a_p[region_id] + 
//      b_p[region_id] .* (x * beta_p) + 
//      c_p[region_id] .* (z * gamma_p);
//    vector[n] eta_m = lambda_m[region_id] .* to_vector(tau[time_id, 2]) + 
//      a_m[region_id] + 
//      b_m[region_id] .* (x * beta_m) + 
//      c_m[region_id] .* (z * gamma_m);
//    vector[n] eta_s = lambda_s[region_id] .* to_vector(tau[time_id, 3]) + 
//      a_s[region_id] + 
//      b_s[region_id] .* (x * beta_s) + 
//      c_s[region_id] .* (z * gamma_s);
//      
//    for(i in 1:n) {
//      log_lik[i] = bernoulli_logit_lpmf(y_p[i] | eta_p[i]) +
//                   bernoulli_logit_lpmf(y_m[i] | eta_m[i]) +
//                   bernoulli_logit_lpmf(y_s[i] | eta_s[i]);
//    }
//  }
}
