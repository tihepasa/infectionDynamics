## Bernoulli models
library(dplyr)
library(cmdstanr)
library(rstan)
library(loo)

set.seed(1)
# # read the data which includes:
# T, the number of weeks (time points)
# N, the number of towns (regions)
# n_pertussis, TxN matrix of the numbers of the observed deaths caused by pertussis
# n_measles, TxN matrix of the numbers of the observed deaths caused by measles
# n_smallpox, TxN matrix of the numbers of the observed deaths caused by smallpox
# is_pertussis, TxN matrix of the dichotomous observations of the deaths caused by pertussis
# is_measles, TxN matrix of the dichotomous observations of the deaths caused by measles
# is_smallpox, TxN matrix of the dichotomous observations of the deaths caused by smallpox
# missing, TxN matrix indicating the missingness of observations, 0 = observed, 1 = missing
# neighbors, NxM matrix of the row numbers of neighbors of each site, M = maximum number of neighbors
# n_neighbors, vector of length N indicating the number of neighbors of each site
# below, these are formatted to correspond the notation in the article
par <- readRDS("infectionpars_week.rds")

# dichotomous disease variables
obs_per <- as.matrix(par$is_pertussis)
obs_mea <- as.matrix(par$is_measles)
obs_sma <- as.matrix(par$is_smallpox)

# matrices for neighbors' values
p_per <- p_mea <- p_sma <- matrix(0, par$N, par$T)

for (t in 1:par$T) {
  for (i in 1:par$N) {
    p_per[i, t] <- mean(obs_per[par$neighbors[i, 1:par$n_neighbors[i]], t], na.rm = TRUE)
    p_mea[i, t] <- mean(obs_mea[par$neighbors[i, 1:par$n_neighbors[i]], t], na.rm = TRUE)
    p_sma[i, t] <- mean(obs_sma[par$neighbors[i, 1:par$n_neighbors[i]], t], na.rm = TRUE)
  }
}

# collect the responses and predictors, and set ids for sites, and 
# create time and month (= actually week in this case) indices
d <- data.frame(
  pertussis = unlist(par$is_pertussis),
  measles = unlist(par$is_measles),
  smallpox = unlist(par$is_smallpox),
  z_pertussis = c(p_per),
  z_measles = c(p_mea),
  z_smallpox = c(p_sma),
  id = factor(1:par$N),
  time = factor(rep(1:par$T, each = par$N)),
  #month = factor(rep(rep(1:12, length.out = par$T)[-1], each = par$N)),
  month = seq.Date(as.Date("1820-01-03"), by = "week", length.out = 1617) %>% 
    format("%m") %>% 
    as.numeric() %>% 
    rep(each = par$N) %>% 
    factor()
)

# mutate for the lags
d <- d %>% group_by(id) %>% 
  mutate(
    lag1_pertussis = dplyr::lag(pertussis),
    lag1_measles = dplyr::lag(measles),
    lag1_smallpox = dplyr::lag(smallpox),
    lag2_pertussis = dplyr::lag(pertussis, 2),
    lag2_measles = dplyr::lag(measles, 2),
    lag2_smallpox = dplyr::lag(smallpox, 2),
    lag3_pertussis = dplyr::lag(pertussis, 3),
    lag3_measles = dplyr::lag(measles, 3),
    lag3_smallpox = dplyr::lag(smallpox, 3),
    lag4_pertussis = dplyr::lag(pertussis, 4),
    lag4_measles = dplyr::lag(measles, 4),
    lag4_smallpox = dplyr::lag(smallpox, 4),
    z_lag1_pertussis = dplyr::lag(z_pertussis),
    z_lag1_measles = dplyr::lag(z_measles),
    z_lag1_smallpox = dplyr::lag(z_smallpox),
    z_lag2_pertussis = dplyr::lag(z_pertussis, 2),
    z_lag2_measles = dplyr::lag(z_measles, 2),
    z_lag2_smallpox = dplyr::lag(z_smallpox, 2), 
    z_lag3_pertussis = dplyr::lag(z_pertussis, 3),
    z_lag3_measles = dplyr::lag(z_measles, 3),
    z_lag3_smallpox = dplyr::lag(z_smallpox, 3),
    z_lag4_pertussis = dplyr::lag(z_pertussis, 4),
    z_lag4_measles = dplyr::lag(z_measles, 4),
    z_lag4_smallpox = dplyr::lag(z_smallpox, 4)) %>% 
  ungroup()


# filter complete cases
# in order to get the data for the model without the last two years, include
# the commented filtering part below and all the rest works as it is
d <- d %>% 
  select(!c(z_pertussis, z_measles, z_smallpox)) %>% 
  tidyr::drop_na() #%>%
#filter(as.numeric(time) < (par$T - 24))

# drop factor levels to get proper indices
d$time <- droplevels(d$time)

# select variables needed for fitting the model
# responses
y_p <- d[["pertussis"]]
y_m <- d[["measles"]]
y_s <- d[["smallpox"]]

# number of observations
n <- nrow(d)

# define how many lags to use
cond_x <- paste0("lag", 1:4)
cond_z <- paste0("z_lag", 1:4)

# focal predictors
x <- d %>% 
  select(starts_with(cond_x))
# neighbors' predictors
z <- d %>% 
  select(starts_with(cond_z))
# numbers of focal and neighboring covariates
K_x <- ncol(x)
K_z <- ncol(z)
# indices for observed months, overall times, and sites
time_id <- as.integer(d$time)
month_id <- as.integer(d$month)
region_id <- as.integer(d$id)

# collect the data to insert the stan model
standata <- tibble::lst(T = par$T, N = par$N, n,
                        y_p, y_m, y_s, K_x, K_z, x, z,
                        time_id, month_id, region_id)

# initial values for tau, based on data
tau_init <- log(
  0.01 + cbind(
    colMeans(par$n_pertussis[, -1], na.rm = TRUE),
    colMeans(par$n_measles[, -1], na.rm = TRUE),
    colMeans(par$n_smallpox[, -1], na.rm = TRUE)
  )
)

# set initial values based on earlier MCMC runs for faster adaptation, can be altered
# in case of omitting seasonal effect, comment out the rows with seasonal_
# in case of independent diseases, comment out the row with L = ...
inits <- replicate(
  4, 
  list(
    sigma_a_p = 0.7, sigma_a_m = 0.4, sigma_a_s = 0.3,
    sigma_b_p = 0.4, sigma_b_m = 0.3, sigma_b_s = 0.3,
    sigma_c_p = 0.6, sigma_c_m = 0.5, sigma_c_s = 0.4,
    sigma_lambda_p = 0.4, sigma_lambda_m = 0.2, sigma_lambda_s = 0.2,
    sigma_tau = c(0.2, 0.2, 0.2),
    seasonal_p_raw = rep(0, 11),
    seasonal_m_raw = rep(0, 11),
    seasonal_s_raw = rep(0, 11),
    lambda_p_raw = numeric(par$N - 1),
    lambda_m_raw = numeric(par$N - 1),
    lambda_s_raw = numeric(par$N - 1),
    tau = tau_init,
    L = matrix(c(1, 0.3, 0, 0, 0.95, 0.25, 0, 0, 0.95), 3, 3),
    a_p_raw = numeric(par$N),
    a_m_raw = numeric(par$N),
    a_s_raw = numeric(par$N),
    b_p_raw = numeric(par$N),
    b_m_raw = numeric(par$N),
    b_s_raw = numeric(par$N),
    c_p_raw = numeric(par$N),
    c_m_raw = numeric(par$N),
    c_s_raw = numeric(par$N),
    beta_p = rep(c(2, 0, 0), length.out = K_x),
    beta_m = rep(c(0, 2, 0), length.out = K_x),
    beta_s = rep(c(0, 0, 2), length.out = K_x),
    gamma_p = rep(c(2, 0, 0), length.out = K_z),
    gamma_m = rep(c(0, 2, 0), length.out = K_z),
    gamma_s = rep(c(0, 0, 2), length.out = K_z)
  ), simplify = FALSE)

# create a CmdStanModel object from the Stan program
model <- cmdstan_model("bern_dependent_s.stan", stanc_options = list("O1"))

# fit the model, 4 chains with 2500 warm-up iterations and 5000 posterior samples
fit <- model$sample(data = standata,
                    iter_sampling = 5000,
                    iter_warmup = 2500,
                    parallel_chains = 4,
                    chains = 4,
                    init = inits,
                    refresh = 100,
                    save_warmup = FALSE)

# save the model with a suitable name for future use as a CmdStan object or as a Stan object
# the Stan fits are needed at least for the predictions in ..._pred.R files,
# the CmdStan objects can be useful when calculating the leave-one-out cross-validations
fit$save_object(paste0("cmdstanfit_bern_week_lag4_noseason.rds"))
fit$diagnostic_summary()
#out <- rstan::read_stan_csv(fit$output_files())
#saveRDS(out, file = paste0("fit_bern_week_lag4_noseason.rds"))

