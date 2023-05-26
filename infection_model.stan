// Simultaneous modelling of three response diseases.

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
  // Response variables, first time point is fixed
  array[T - 1, N] int<lower=0, upper=1> pertussis;
  array[T - 1, N] int<lower=0, upper=1> measles;
  array[T - 1, N] int<lower=0, upper=1> smallpox;

  // Variables below are here same for all responses
  // In general we could have different predictors and missing data patters
  int<lower=0> K_x; // number of predictors based on site's own characteristics
  int<lower=0> K_z; // number of predictors based on neighbors' characteristics
  array[T - 1] matrix[N, K_x] x; // predictors
  array[T - 1] matrix[N, K_z] z; // predictors
  int<lower=0> n_obs; // number of non-missing observations
  array[n_obs] int<lower=0> ind; // indices of non-missing observations
}
transformed data {
  // Scale parameter for the sum constrained prior
  real scaleN = inv(sqrt(1 - inv(N)));
  real scale12 = inv(sqrt(1 - inv(12)));
  // Create QR decomposition for sum-to-zero constraint
  vector[2 * N] QN = create_Q(N);
  vector[2 * 12] Q12 = create_Q(12);
  // indexing variables for months
  array[T - 1] int month;
  for (t in 1:(T - 1)) {
    month[t] = 1 + t % 12;
  }
}
parameters {
  // Population-level parameters
  vector[K_x] beta_p;
  vector[K_x] beta_m;
  vector[K_x] beta_s;
  vector[K_z] gamma_p;
  vector[K_z] gamma_m;
  vector[K_z] gamma_s;

  // Standard deviations of factor loadings lambda
  real<lower=0> sigma_lambda_p;
  real<lower=0> sigma_lambda_m;
  real<lower=0> sigma_lambda_s;

  // Raw parameters for seasonal effects with sum-to-zero constraint
  vector[11] seasonal_p_raw;
  vector[11] seasonal_m_raw;
  vector[11] seasonal_s_raw;

  // Raw parameters for factor loadings lambda
  vector[N - 1] lambda_p_raw;
  vector[N - 1] lambda_m_raw;
  vector[N - 1] lambda_s_raw;

  // Latent factors tau for each three responses (p, m, s)
  array[T - 1] vector[3] tau;
  // Standard deviations of the latent factor noise terms
  vector<lower=0>[3] sigma_tau;
  // Cholesky factor of the correlation matrix of latent factors
  cholesky_factor_corr[3] L;

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
}
transformed parameters {
  vector[12] seasonal_p = sum_to_zero(seasonal_p_raw, Q12, 12);
  vector[12] seasonal_m = sum_to_zero(seasonal_m_raw, Q12, 12);
  vector[12] seasonal_s = sum_to_zero(seasonal_s_raw, Q12, 12);
  vector[N] lambda_p = 1 + sigma_lambda_p * sum_to_zero(lambda_p_raw, QN, N);
  vector[N] lambda_m = 1 + sigma_lambda_m * sum_to_zero(lambda_m_raw, QN, N);
  vector[N] lambda_s = 1 + sigma_lambda_s * sum_to_zero(lambda_s_raw, QN, N);

  vector[N] a_p = sigma_a_p * a_p_raw;
  vector[N] a_m = sigma_a_m * a_m_raw;
  vector[N] a_s = sigma_a_s * a_s_raw;
  vector[N] b_p = 1 + sigma_b_p * b_p_raw;
  vector[N] b_s = 1 + sigma_b_s * b_s_raw;
  vector[N] b_m = 1 + sigma_b_m * b_m_raw;
  vector[N] c_p = 1 + sigma_c_p * c_p_raw;
  vector[N] c_s = 1 + sigma_c_s * c_s_raw;
  vector[N] c_m = 1 + sigma_c_m * c_m_raw;
}
model {
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

  sigma_tau ~ gamma(2, 1);
  L ~ lkj_corr_cholesky(1);
  tau[1] ~ normal(-2, 2);
  tau[2:(T-1)] ~ multi_normal_cholesky(tau[1:(T - 2)], diag_pre_multiply(sigma_tau, L));

  sigma_a_p ~ gamma(2, 1);
  sigma_a_m ~ gamma(2, 1);
  sigma_a_s ~ gamma(2, 1);

  sigma_b_p ~ gamma(2, 1);
  sigma_b_m ~ gamma(2, 1);
  sigma_b_s ~ gamma(2, 1);

  sigma_c_p ~ gamma(2, 1);
  sigma_c_m ~ gamma(2, 1);
  sigma_c_s ~ gamma(2, 1);

  a_p_raw ~ std_normal();
  a_m_raw ~ std_normal();
  a_s_raw ~ std_normal();
  b_p_raw ~ std_normal();
  b_m_raw ~ std_normal();
  b_s_raw ~ std_normal();
  c_p_raw ~ std_normal();
  c_m_raw ~ std_normal();
  c_s_raw ~ std_normal();

  {
    vector[N * (T - 1)] logit_prob;
    for (t in 1:(T - 1)) {
      logit_prob[((t - 1) * N + 1) : (t * N)] =
         a_p + seasonal_p[month[t]] + lambda_p * tau[t, 1] +
         rows_dot_product(b_p, x[t] * beta_p) +
         rows_dot_product(c_p, z[t] * gamma_p);
    }
    to_array_1d(pertussis)[ind] ~ bernoulli_logit(logit_prob[ind]);

    for (t in 1:(T - 1)) {
      logit_prob[((t - 1) * N + 1) : (t * N)] =
        a_m + seasonal_m[month[t]] + lambda_m * tau[t, 2] +
        rows_dot_product(b_m, x[t] * beta_m) +
        rows_dot_product(c_m, z[t] * gamma_m);
    }
    to_array_1d(measles)[ind] ~ bernoulli_logit(logit_prob[ind]);

    for (t in 1:(T - 1)) {
      logit_prob[((t - 1) * N + 1) : (t * N)] =
        a_s + seasonal_s[month[t]] + lambda_s * tau[t, 3] +
        rows_dot_product(b_s, x[t] * beta_s) +
        rows_dot_product(c_s, z[t] * gamma_s);
    }
    to_array_1d(smallpox)[ind] ~ bernoulli_logit(logit_prob[ind]);
  }
}
generated quantities {
  matrix[3, 3] corr_tau = multiply_lower_tri_self_transpose(L);
}
