# Simulate the predictions by Bernoulli 2-month lag model
library(rstan)
library(cmdstanr)
library(dplyr)
library(tibble)
library(tidyr)

# read the parameter file
par <- readRDS("infectionpars.rds")

# read and extract the fit (cmdstanfit)
fit <- readRDS("cmdstanfit_bern_dependent_lag2.rds")

# include the next row when removing the last two years
#par$missing[, ((par$T - 23):par$T)] <- 1

# for the first time point we assume zeros for missing sites
obs_per <- par$is_pertussis %>% as.matrix()
obs_per[par$missing == 1] <- 0
obs_mea <- par$is_measles %>% as.matrix()
obs_mea[par$missing == 1] <- 0
obs_sma <- par$is_smallpox %>% as.matrix()
obs_sma[par$missing == 1] <- 0

# get indices for the observations which are used as a response and the observations 
# which are not missing but their covariates are so that they are not used as a response
# missingness pattern is the same for all diseases
obs <- obs_nocov <- matrix(0, nrow = par$N, ncol = (par$T - 1))
for (t in 2:par$T) {
  for (i in 1:par$N) {
    if (par$missing[i, t] == 0 && par$missing[i, t - 1] == 0 &&
        any(par$missing[par$neighbors[i, 1:par$n_neighbors[i]], t - 1] == 0)) {
      obs[i, t - 1] <- 1
    }
    
    if (par$missing[i, t] == 0 && (par$missing[i, t - 1] == 1 ||
                                   all(par$missing[par$neighbors[i, 1:par$n_neighbors[i]], t - 1] == 1))) {
      obs_nocov[i, t - 1] <- 1
    }
  }
}


# variables for the neighbor predictors
p_per <- p_mea <- p_sma <- matrix(0, ncol = par$T, nrow = par$N)

# calculate the means over neighbors' values
for (i in 1:par$N) {
  for (t in 1:par$T) {
    p_per[i, t] <- mean((1 * (obs_per[par$neighbors[i, 1:par$n_neighbors[i]], t])))
    p_mea[i, t] <- mean((1 * (obs_mea[par$neighbors[i, 1:par$n_neighbors[i]], t])))
    p_sma[i, t] <- mean((1 * (obs_sma[par$neighbors[i, 1:par$n_neighbors[i]], t])))
  }
}

# number of iterations during the MCMC estimation of the fit
Nsim <- nrow(fit$draws("a_p", format = "matrix"))

# matrix form of the neighbor predictors, needed for matrix multiplication
p_per_temp <- matrix(p_per[, 2], nrow = Nsim, ncol = par$N, byrow = TRUE)
p_mea_temp <- matrix(p_mea[, 2], nrow = Nsim, ncol = par$N, byrow = TRUE)
p_sma_temp <- matrix(p_sma[, 2], nrow = Nsim, ncol = par$N, byrow = TRUE)

p_per_temp_lag <- matrix(p_per[, 1], nrow = Nsim, ncol = par$N, byrow = TRUE)
p_mea_temp_lag <- matrix(p_mea[, 1], nrow = Nsim, ncol = par$N, byrow = TRUE)
p_sma_temp_lag <- matrix(p_sma[, 1], nrow = Nsim, ncol = par$N, byrow = TRUE)

# get the parameter estimates
b1 <- cbind(fit$draws("beta_p", format = "matrix")[, 1],
            fit$draws("gamma_p", format = "matrix")[, 1],
            fit$draws("beta_p", format = "matrix")[, 2],
            fit$draws("gamma_p", format = "matrix")[, 2],
            fit$draws("beta_p", format = "matrix")[, 3],
            fit$draws("gamma_p", format = "matrix")[, 3])

b1_lag <- cbind(fit$draws("beta_p", format = "matrix")[, 4],
                fit$draws("gamma_p", format = "matrix")[, 4],
                fit$draws("beta_p", format = "matrix")[, 5],
                fit$draws("gamma_p", format = "matrix")[, 5],
                fit$draws("beta_p", format = "matrix")[, 6],
                fit$draws("gamma_p", format = "matrix")[, 6])

b2 <- cbind(fit$draws("beta_m", format = "matrix")[, 1],
            fit$draws("gamma_m", format = "matrix")[, 1],
            fit$draws("beta_m", format = "matrix")[, 2],
            fit$draws("gamma_m", format = "matrix")[, 2],
            fit$draws("beta_m", format = "matrix")[, 3],
            fit$draws("gamma_m", format = "matrix")[, 3])

b2_lag <- cbind(fit$draws("beta_m", format = "matrix")[, 4],
                fit$draws("gamma_m", format = "matrix")[, 4],
                fit$draws("beta_m", format = "matrix")[, 5],
                fit$draws("gamma_m", format = "matrix")[, 5],
                fit$draws("beta_m", format = "matrix")[, 6],
                fit$draws("gamma_m", format = "matrix")[, 6])

b3 <- cbind(fit$draws("beta_s", format = "matrix")[, 1],
            fit$draws("gamma_s", format = "matrix")[, 1],
            fit$draws("beta_s", format = "matrix")[, 2],
            fit$draws("gamma_s", format = "matrix")[, 2],
            fit$draws("beta_s", format = "matrix")[, 3],
            fit$draws("gamma_s", format = "matrix")[, 3])

b3_lag <- cbind(fit$draws("beta_s", format = "matrix")[, 4],
                fit$draws("gamma_s", format = "matrix")[, 4],
                fit$draws("beta_s", format = "matrix")[, 5],
                fit$draws("gamma_s", format = "matrix")[, 5],
                fit$draws("beta_s", format = "matrix")[, 6],
                fit$draws("gamma_s", format = "matrix")[, 6])

a1 <- fit$draws("a_p", format = "matrix")
a2 <- fit$draws("a_m", format = "matrix")
a3 <- fit$draws("a_s", format = "matrix")

# seasonal effect, omit in case of seasonless model
s1 <- fit$draws("seasonal_p", format = "matrix")
s2 <- fit$draws("seasonal_m", format = "matrix")
s3 <- fit$draws("seasonal_s", format = "matrix")

tau1 <- fit$draws("tau", format = "matrix")[, 1:(par$T - 1)]
tau2 <- fit$draws("tau", format = "matrix")[, (par$T):(2 * (par$T - 1))]
tau3 <- fit$draws("tau", format = "matrix")[, (2 * (par$T - 1) + 1):(3 * (par$T - 1))]


lam1 <- fit$draws("lambda_p", format = "matrix")
lam2 <- fit$draws("lambda_m", format = "matrix")
lam3 <- fit$draws("lambda_s", format = "matrix")

# array form of regression coefficients for matrix multiplication
coefs1 <- array(c(fit$draws("b_p", format = "matrix"),
                  fit$draws("c_p", format = "matrix"),
                  fit$draws("b_p", format = "matrix"),
                  fit$draws("c_p", format = "matrix"),
                  fit$draws("b_p", format = "matrix"),
                  fit$draws("c_p", format = "matrix")),
                dim = c(Nsim, par$N, 6))

coefs2 <- array(c(fit$draws("b_m", format = "matrix"),
                  fit$draws("c_m", format = "matrix"),
                  fit$draws("b_m", format = "matrix"),
                  fit$draws("c_m", format = "matrix"),
                  fit$draws("b_m", format = "matrix"),
                  fit$draws("c_m", format = "matrix")),
                dim = c(Nsim, par$N, 6))

coefs3 <- array(c(fit$draws("b_s", format = "matrix"),
                  fit$draws("c_s", format = "matrix"),
                  fit$draws("b_s", format = "matrix"),
                  fit$draws("c_s", format = "matrix"),
                  fit$draws("b_s", format = "matrix"),
                  fit$draws("c_s", format = "matrix")),
                dim = c(Nsim, par$N, 6))

# matrix form of the focal site's predictors, needed for matrix multiplication
temp_obs_p1 <- matrix(obs_per[, 2], nrow = Nsim, ncol = par$N, byrow = TRUE)
temp_obs_m1 <- matrix(obs_mea[, 2], nrow = Nsim, ncol = par$N, byrow = TRUE)
temp_obs_s1 <- matrix(obs_sma[, 2], nrow = Nsim, ncol = par$N, byrow = TRUE)
temp_obs_p1_lag <- matrix(obs_per[, 1], nrow = Nsim, ncol = par$N, byrow = TRUE)
temp_obs_m1_lag <- matrix(obs_mea[, 1], nrow = Nsim, ncol = par$N, byrow = TRUE)
temp_obs_s1_lag <- matrix(obs_sma[, 1], nrow = Nsim, ncol = par$N, byrow = TRUE)

# average probabilities in time
# row = time, columns = mean, 2.5% and 97.5% quantile
pr_time_p <- pr_time_m <- pr_time_s <- matrix(NA, nrow = (par$T - 1), ncol = 3)

# temporary probabilites
pr_temp_p <- array(NA, dim = c(Nsim, par$N, par$T - 1))
pr_temp_m <- array(NA, dim = c(Nsim, par$N, par$T - 1))
pr_temp_s <- array(NA, dim = c(Nsim, par$N, par$T - 1))

# function to get mean and 95 percent (posterior) interval of a vector x
descript <- function(x) {c(mean = mean(x), quantile(x, probs = 0.025), quantile(x, probs = 0.975))}

for (t in 2:(par$T - 1)) {
  # number of missing observations in current time point
  n_missing <- sum(par$missing[, t + 1] == 1)
  
  ## pertussis
  # probabilities for all iterations and sites in current time point
  # in case of model without seasonal effect, remove s1/s2/s3
  pr_temp_p[,, t] <- plogis(a1 +
                              c(s1[, 1 + t %% 12]) +
                              lam1 * c(tau1[, t]) + 
                              c(b1[, 1]) * c(coefs1[, , 1]) * temp_obs_p1 + c(b1[, 2]) * c(coefs1[, , 2]) * p_per_temp +
                              c(b1[, 3]) * c(coefs1[, , 3]) * temp_obs_m1 + c(b1[, 4]) * c(coefs1[, , 4]) * p_mea_temp + 
                              c(b1[, 5]) * c(coefs1[, , 5]) * temp_obs_s1 + c(b1[, 6]) * c(coefs1[, , 6]) * p_sma_temp +
                              c(b1_lag[, 1]) * c(coefs1[, , 1]) * temp_obs_p1_lag + c(b1_lag[, 2]) * c(coefs1[, , 2]) * p_per_temp_lag +
                              c(b1_lag[, 3]) * c(coefs1[, , 3]) * temp_obs_m1_lag + c(b1_lag[, 4]) * c(coefs1[, , 4]) * p_mea_temp_lag + 
                              c(b1_lag[, 5]) * c(coefs1[, , 5]) * temp_obs_s1_lag + c(b1_lag[, 6]) * c(coefs1[, , 6]) * p_sma_temp_lag)
  
  # posterior interval and mean of the probabilities for current time point
  pr_time_p[t, ] <- descript(rowMeans(pr_temp_p[,, t]))
  
  # generate data to use at the next time point, use observed ones as they are and 
  # simulate missing ones, take transpose of the result later since now row = site and column = iteration
  temp_obs_p <- matrix(obs_per[, t + 1], ncol = Nsim, nrow = par$N, byrow = FALSE)
  temp_obs_p[par$missing[, t + 1] == 1, ] <- rbinom(Nsim * n_missing,
                                                    1,
                                                    (pr_temp_p[, par$missing[, t + 1] == 1, t])) > 0
  
  ## measles
  pr_temp_m[,, t] <- plogis(a2 +
                              c(s2[, 1 + t %% 12]) +
                              lam2 * c(tau2[, t]) + 
                              c(b2[, 1]) * c(coefs2[, , 1]) * temp_obs_p1 + c(b2[, 2]) * c(coefs2[, , 2]) * p_per_temp +
                              c(b2[, 3]) * c(coefs2[, , 3]) * temp_obs_m1 + c(b2[, 4]) * c(coefs2[, , 4]) * p_mea_temp + 
                              c(b2[, 5]) * c(coefs2[, , 5]) * temp_obs_s1 + c(b2[, 6]) * c(coefs2[, , 6]) * p_sma_temp +
                              c(b2_lag[, 1]) * c(coefs2[, , 1]) * temp_obs_p1_lag + c(b2_lag[, 2]) * c(coefs2[, , 2]) * p_per_temp_lag +
                              c(b2_lag[, 3]) * c(coefs2[, , 3]) * temp_obs_m1_lag + c(b2_lag[, 4]) * c(coefs2[, , 4]) * p_mea_temp_lag + 
                              c(b2_lag[, 5]) * c(coefs2[, , 5]) * temp_obs_s1_lag + c(b2_lag[, 6]) * c(coefs2[, , 6]) * p_sma_temp_lag)
  
  pr_time_m[t, ] <- descript(rowMeans(pr_temp_m[,, t]))
  
  temp_obs_m <- matrix(obs_mea[, t + 1], ncol = Nsim, nrow = par$N, byrow = FALSE)
  temp_obs_m[par$missing[, t + 1] == 1, ] <- rbinom(Nsim * n_missing,
                                                    1,
                                                    (pr_temp_m[, par$missing[, t + 1] == 1, t])) > 0
  
  ## smallpox
  pr_temp_s[,, t] <- plogis(a3 +
                              c(s3[, 1 + t %% 12]) +
                              lam3 * c(tau3[, t]) + 
                              c(b3[, 1]) * c(coefs3[, , 1]) * temp_obs_p1 + c(b3[, 2]) * c(coefs3[, , 2]) * p_per_temp +
                              c(b3[, 3]) * c(coefs3[, , 3]) * temp_obs_m1 + c(b3[, 4]) * c(coefs3[, , 4]) * p_mea_temp + 
                              c(b3[, 5]) * c(coefs3[, , 5]) * temp_obs_s1 + c(b3[, 6]) * c(coefs3[, , 6]) * p_sma_temp +
                              c(b3_lag[, 1]) * c(coefs3[, , 1]) * temp_obs_p1_lag + c(b3_lag[, 2]) * c(coefs3[, , 2]) * p_per_temp_lag +
                              c(b3_lag[, 3]) * c(coefs3[, , 3]) * temp_obs_m1_lag + c(b3_lag[, 4]) * c(coefs3[, , 4]) * p_mea_temp_lag + 
                              c(b3_lag[, 5]) * c(coefs3[, , 5]) * temp_obs_s1_lag + c(b3_lag[, 6]) * c(coefs3[, , 6]) * p_sma_temp_lag)
  
  pr_time_s[t, ] <- descript(rowMeans(pr_temp_s[,, t]))
  
  temp_obs_s <- matrix(obs_sma[, t + 1], ncol = Nsim, nrow = par$N, byrow = FALSE)
  temp_obs_s[par$missing[, t + 1] == 1, ] <- rbinom(Nsim * n_missing,
                                                    1,
                                                    (pr_temp_s[, par$missing[, t + 1] == 1, t])) > 0
  
  # focal site's past, take the transpose now
  temp_obs_p1_lag <- temp_obs_p1
  temp_obs_m1_lag <- temp_obs_m1
  temp_obs_s1_lag <- temp_obs_s1
  
  temp_obs_p1 <- t(temp_obs_p)
  temp_obs_m1 <- t(temp_obs_m)
  temp_obs_s1 <- t(temp_obs_s)
  
  # means over neighbors based on the current data, if clause is needed since rowMeans fails if n = 1
  p_per_temp_lag <- p_per_temp
  p_mea_temp_lag <- p_mea_temp
  p_sma_temp_lag <- p_sma_temp
  
  for (i in 1:par$N) {
    if (par$n_neighbors[i] == 1) {
      p_per_temp[, i] <- ((temp_obs_p1[, par$neighbors[i, 1:par$n_neighbors[i]]]))
      p_mea_temp[, i] <- ((temp_obs_m1[, par$neighbors[i, 1:par$n_neighbors[i]]]))
      p_sma_temp[, i] <- ((temp_obs_s1[, par$neighbors[i, 1:par$n_neighbors[i]]]))
    } else {
      p_per_temp[, i] <- rowMeans((temp_obs_p1[, par$neighbors[i, 1:par$n_neighbors[i]]]))
      p_mea_temp[, i] <- rowMeans((temp_obs_m1[, par$neighbors[i, 1:par$n_neighbors[i]]]))
      p_sma_temp[, i] <- rowMeans((temp_obs_s1[, par$neighbors[i, 1:par$n_neighbors[i]]]))
    }
  }
}

# take local averages over iterations and time points 
pr_site_p <- apply(pr_temp_p, c(1, 2), mean, na.rm = TRUE) %>% colMeans()
pr_site_m <- apply(pr_temp_m, c(1, 2), mean, na.rm = TRUE) %>% colMeans()
pr_site_s <- apply(pr_temp_s, c(1, 2), mean, na.rm = TRUE) %>% colMeans()

# now the results can be saved or further used, the variables worth inspecting are
# pr_time_ and pr_site_
saveRDS(list(pr_time_p, pr_time_m, pr_time_s, pr_site_p, pr_site_m, pr_site_s), "sims_dataa_lag2.rds")

