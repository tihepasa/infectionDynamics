library(cmdstanr)

# read the data which includes:
# T, the number of months (time points)
# N, the number of towns (regions)
# pertussis, TxN matrix of the dichotomous observations of the deaths caused by pertussis
# measles, TxN matrix of the dichotomous observations of the deaths caused by measles
# smallpox, TxN matrix of the dichotomous observations of the deaths caused by smallpox
# K_x, the number of the covariates related to the town i
# K_z, the number of the covariates related to the neighborhood of the town i
# x, (T)x(N)x(K_x) array of the covariates considering the focal town i: observations of death occurrences
#   caused by pertussis, measles and smallpox in previous month
# z, (T)x(N)x(K_z) array of the covariates considering the neighborhood of town i: means of death occurrences
#   caused by pertussis, measles and smallpox in previous month over the neighborhood
# n_obs, the number of non-missing observations that are also preceded by non-missing observations
# ind, the indexes of the non-missing observations preceded by non-missing observations in the matrices
#   pertussis, measles and smallpox
pars <- readRDS("infectionpars.rds")

# compile the model
model <- cmdstan_model("infection_model.stan", stanc_options = list("O1"))

# set initial values based on earlier MCMC runs for faster adaptation, can be altered
inits <- replicate(4,
                   list(
                     betax_p = c(1.6, 0.1, 0),
                     betax_m = c(0.1, 1.9, 0.2),
                     betax_s = c(0.1, 0.2, 2.4),
                     betaz_p = c(1.3, 0.1, 0.1),
                     betaz_m = c(0.1, 2.2, 0.2),
                     betaz_s = c(0.2, 0.4, 2.6),
                     sigma_lambda_p = 0.3,
                     sigma_lambda_m = 0.2,
                     sigma_lambda_s = 0.2,
                     seasonal_p_raw = rep(0, 11),
                     seasonal_m_raw = rep(0, 11),
                     seasonal_s_raw = rep(0, 11),
                     lambda_p_raw = numeric(pars$N - 1),
                     lambda_m_raw = numeric(pars$N - 1),
                     lambda_s_raw = numeric(pars$N - 1),
                     tau = qlogis(
                       0.01 + cbind(
                         rowMeans(pars$pertussis),
                         rowMeans(pars$measles),
                         rowMeans(pars$smallpox)
                       )
                     ),
                     sigma_tau = c(0.15, 0.17, 0.2),
                     L = matrix(c(1, 0.32, 0.06, 0, 0.94, 0.23, 0, 0, 0.96), 3, 3),
                     sigma_a_p = 0.8,
                     sigma_a_m = 0.3,
                     sigma_a_s = 0.1,
                     sigma_bx_p = 0.4,
                     sigma_bx_m = 0.2,
                     sigma_bx_s = 0.2,
                     sigma_bz_p = 0.5,
                     sigma_bz_m = 0.4,
                     sigma_bz_s = 0.3,
                     a_p_raw = numeric(pars$N),
                     a_m_raw = numeric(pars$N),
                     a_s_raw = numeric(pars$N),
                     bx_p_raw = numeric(pars$N),
                     bx_m_raw = numeric(pars$N),
                     bx_s_raw = numeric(pars$N),
                     bz_p_raw = numeric(pars$N),
                     bz_m_raw = numeric(pars$N),
                     bz_s_raw = numeric(pars$N)
                   ),
                   simplify = FALSE)

set.seed(1)

# fit the model, 4 chains with 2500 warm-up iterations and 5000 posterior samples
fit <- model$sample(data = pars,
                    init = inits,
                    iter_sampling = 5000,
                    iter_warmup = 2500,
                    parallel_chains = 4,
                    refresh = 10,
                    save_warmup = FALSE)
