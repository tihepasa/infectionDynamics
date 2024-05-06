## Bernoulli models
library(dplyr)
library(cmdstanr)
library(rstan)
library(loo)

set.seed(1)
# read parameters
par <- readRDS("infectionpars.rds")

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
# create time and month indices
d <- data.frame(
  pertussis = unlist(par$n_pertussis[, -1]),
  measles = unlist(par$n_measles[, -1]),
  smallpox = unlist(par$n_smallpox[, -1]),
  lag_pertussis = unlist(par$is_pertussis[, -par$T]),
  lag_measles = unlist(par$is_measles[, -par$T]),
  lag_smallpox = unlist(par$is_smallpox[, -par$T]),
  z_pertussis = c(p_per[, -par$T]),
  z_measles = c(p_mea[, -par$T]),
  z_smallpox = c(p_sma[, -par$T]),
  id = factor(1:par$N),
  time = factor(rep(2:par$T, each = par$N)),
  month = factor(rep(rep(1:12, length.out = par$T)[-1], each = par$N))
)

# filter complete cases
# in order to get the data for the model without the last two years, include
# the commented filtering part below and all the rest works as it is
d <- d %>% 
  tidyr::drop_na() #%>%
#filter(as.numeric(time) < (par$T - 24))

# select variables needed for fitting the model
# responses
y_p <- d[["pertussis"]]
y_m <- d[["measles"]]
y_s <- d[["smallpox"]]

# number of observations
n <- nrow(d)

# focal predictors
x <- d %>% 
  select(starts_with("lag_"))
# neighbors' predictors
z <- d %>% 
  select(starts_with("z_"))
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
    alpha_phi_p = -0.5, alpha_phi_m = -1.1, alpha_phi_s = -0.8,
    sigma_phi_p = 0.7, sigma_phi_m = 0.5, sigma_phi_s = 0.5,
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
    phi_p_raw = numeric(par$N),
    phi_m_raw = numeric(par$N),
    phi_s_raw = numeric(par$N),
    beta_p = c(2, 0, 0),
    beta_m = c(0, 2, 0),
    beta_s = c(0, 0, 2),
    gamma_p = c(2, 0, 0),
    gamma_m = c(0, 2, 0),
    gamma_s = c(0, 0, 2)
  ), simplify = FALSE)

# create a CmdStanModel object from the Stan program
model <- cmdstan_model("bern_dependent.stan", stanc_options = list("O1"))
model <- cmdstan_model("bern_independent.stan", stanc_options = list("O1"))
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
#fit$save_object(paste0("cmdstanfit_negbin_dependent.rds"))
#fit$save_object(paste0("cmdstanfit_negbin_independent.rds"))
#fit$save_object(paste0("cmdstanfit_negbin_dependent_s.rds"))
#out <- rstan::read_stan_csv(fit$output_files())
#saveRDS(out, file = paste0("fit_negbin_dependent.rds"))
#saveRDS(out, file = paste0("fit_negbin_independent.rds"))
#saveRDS(out, file = paste0("fit_negbin_dependent_s.rds"))

## some results

# find dominant frequencies of the taus (= incidence factors) via spectral analysis
library(forecast)

# all the dates as a sequence
dates <- seq.Date(from = as.Date("1820-01-01"), to = as.Date("1850-12-01"), by = "month")

# extract taus from the fit
tau <- rstan::extract(fit, "tau")[[1]]

# take mean of the taus and form data frame of them and the dates (not necessary)
taut <- data.frame(mean = c(apply(tau, 2:3, mean)),
                   date = rep(dates[-1], 3))

# find the frequencies
forecast::findfrequency(taut$mean[1:(pars$T - 1)]) # pertussis
forecast::findfrequency(taut$mean[pars$T:(2 * (pars$T - 1))]) # measles
forecast::findfrequency(taut$mean[(2 * pars$T - 1):(3 * (pars$T - 1))]) # smallpox
