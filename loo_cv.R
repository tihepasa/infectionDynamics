# Fit the dependent and independent models (Bernoulli and negative binomial,
# complete data and the last two years) and compare them with leave-one-out cross-validation
# In case of memory issues, fit the models and/or their generated quantities separately, save them,
# and do the cross-validation afterwards for the saved objects.

library(dplyr)
library(cmdstanr)
library(rstan)
library(loo)

set.seed(1)
# read the data which includes:
# T, the number of months (time points)
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

# The code from here to the line where the loos are calculated (till row 213)
# has to be adapted and rerun according to the comments in order to gain the leave-one-out
# cross-validations for each individual model. The current version applies for dependent
# Bernoulli model.

# collect the responses and predictors, and set ids for sites, and
# create time and month indices
# for negative binomial models use pertussis, measles and smallpox on the rows now commented out
d <- data.frame(
  #pertussis = unlist(par$n_pertussis[, -1]),
  #measles = unlist(par$n_measles[, -1]),
  #smallpox = unlist(par$n_smallpox[, -1]),
  pertussis = unlist(par$is_pertussis[, -1]),
  measles = unlist(par$is_measles[, -1]),
  smallpox = unlist(par$is_smallpox[, -1]),
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
# for negative binomial model include the five rows commented out in the beginning of the list
inits <- replicate(
  4, 
  list(
    #alpha_phi_p = -0.5, alpha_phi_m = -1.1, alpha_phi_s = -0.8,
    #sigma_phi_p = 0.7, sigma_phi_m = 0.5, sigma_phi_s = 0.5,
    #phi_p_raw = numeric(par$N),
    #phi_m_raw = numeric(par$N),
    #phi_s_raw = numeric(par$N),
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
    beta_p = c(2, 0, 0),
    beta_m = c(0, 2, 0),
    beta_s = c(0, 0, 2),
    gamma_p = c(2, 0, 0),
    gamma_m = c(0, 2, 0),
    gamma_s = c(0, 0, 2)
  ), simplify = FALSE)

# create a CmdStanModel object from the Stan program
model <- cmdstan_model("bern_dependent.stan", stanc_options = list("O1"))
#model <- cmdstan_model("bern_independent.stan", stanc_options = list("O1"))
#model <- cmdstan_model("bern_dependent_s.stan", stanc_options = list("O1"))
#model <- cmdstan_model("negbin_dependent.stan", stanc_options = list("O1"))
#model <- cmdstan_model("negbin_independent.stan", stanc_options = list("O1"))
#model <- cmdstan_model("negbin_dependent_s.stan", stanc_options = list("O1"))

# fit the model, 4 chains with 2500 warm-up iterations and 5000 posterior samples
fit <- model$sample(data = standata,
                    iter_sampling = 5000,
                    iter_warmup = 2500,
                    parallel_chains = 4,
                    chains = 4,
                    init = inits,
                    refresh = 100,
                    save_warmup = FALSE)


## approximate leave-one-out cross-validation

# use separate generated quantities (to get the log likelihoods), use here a model that
# contains the desired calculations in the generated quantities block of the Stan file
# see also the comments on Stan files in the supplements
model_gq <- cmdstan_model("bern_dependent.stan", stanc_options = list("O1"))
#model_gq <- cmdstan_model("bern_independent.stan", stanc_options = list("O1"))
#model_gq <- cmdstan_model("bern_dependent_s.stan", stanc_options = list("O1"))
#model_gq <- cmdstan_model("negbin_dependent.stan", stanc_options = list("O1"))
#model_gq <- cmdstan_model("negbin_independent.stan", stanc_options = list("O1"))
#model_gq <- cmdstan_model("negbin_dependent_s.stan", stanc_options = list("O1"))

# use the complete data also for the predictive models
fit_gq <- model_gq$generate_quantities(fit, data = standata, parallel_chains = 4)
draws <- fit_gq$draws()

# from the above draws it is convenient to estimate the ELPDs
# relative effective sample sizes have to be added manually
r_eff <- relative_eff(exp(draws))

# make sure that above you have used the desired model and data and choose the
# corresponding row below
loo_bern_dep <- rstan::loo(draws, r_eff = r_eff, cores = 4)
loo_bern_indep <- rstan::loo(draws, r_eff = r_eff, cores = 4)
loo_bern_dep_s <- rstan::loo(draws, r_eff = r_eff, cores = 4)

loo_negbin_dep <- rstan::loo(draws, r_eff = r_eff, cores = 4)
loo_negbin_indep <- rstan::loo(draws, r_eff = r_eff, cores = 4)
loo_negbin_dep_s <- rstan::loo(draws, r_eff = r_eff, cores = 4)

# or in the case of the models omitting the last two years, select manually only the
# log likelihood values considering the last two years, thus the indices

### this part holds an error, the function to use should not be rstan::loo, but loo::elpd,
### the old version is here on comments and the corrected one below without comments
#r_eff <- relative_eff(exp(draws[, , 94344:101153]))
#loo_bern_dep <- rstan::loo(draws[, , 94344:101153], r_eff = r_eff, cores = 4)
#loo_bern_indep <- rstan::loo(draws[, , 94344:101153], r_eff = r_eff, cores = 4)
#loo_bern_dep_s <- rstan::loo(draws[, , 94344:101153], r_eff = r_eff, cores = 4)
#
#loo_negbin_dep <- rstan::loo(draws[, , 94344:101153], r_eff = r_eff, cores = 4)
#loo_negbin_indep <- rstan::loo(draws[, , 94344:101153], r_eff = r_eff, cores = 4)
#loo_negbin_dep_s <- rstan::loo(draws[, , 94344:101153], r_eff = r_eff, cores = 4)

# corrected code starts
loo_bern_dep <- loo::elpd(draws[, , 94344:101153])
loo_bern_indep <- loo::elpd(draws[, , 94344:101153])
loo_bern_dep_s <- loo::elpd(draws[, , 94344:101153])

loo_negbin_dep <- loo::elpd(draws[, , 94344:101153])
loo_negbin_indep <- loo::elpd(draws[, , 94344:101153])
loo_negbin_dep_s <- loo::elpd(draws[, , 94344:101153])
# corrected code ends

# results
loo_bern_dep
loo_bern_indep
loo_bern_dep_s

loo_negbin_dep
loo_negbin_indep
loo_negbin_dep_s

# compare the models
loo::loo_compare(loo_bern_dep, loo_bern_indep, loo_bern_dep_s,
                 loo_negbin_dep, loo_negbin_indep, loo_negbin_dep_s)

