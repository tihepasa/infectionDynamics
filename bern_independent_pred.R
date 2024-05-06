# Simulate the predictions by Bernoulli model
library(rstan)
library(cmdstanr)
library(dplyr)
library(tibble)
library(tidyr)

# read the parameter file
par <- readRDS("infectionpars.rds")

# read and extract the fit (stanfit)
fit <- readRDS("fit_bern_independent.rds")

e <- rstan::extract(fit)

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
Nsim <- nrow(e$a_p)

# matrix form of the neighbor predictors, needed for matrix multiplication
p_per_temp <- matrix(p_per[, 1], nrow = Nsim, ncol = par$N, byrow = TRUE)
p_mea_temp <- matrix(p_mea[, 1], nrow = Nsim, ncol = par$N, byrow = TRUE)
p_sma_temp <- matrix(p_sma[, 1], nrow = Nsim, ncol = par$N, byrow = TRUE)

# get the parameter estimates
b1 <- cbind(e$beta_p,
            e$gamma_p)

b2 <- cbind(e$beta_m,
            e$gamma_m)

b3 <- cbind(e$beta_s, 
            e$gamma_s)

a1 <- e$a_p
a2 <- e$a_m
a3 <- e$a_s

s1 <- e$seasonal_p
s2 <- e$seasonal_m
s3 <- e$seasonal_s

tau1 <- e$tau[,, 1]
tau2 <- e$tau[,, 2]
tau3 <- e$tau[,, 3]

lam1 <- e$lambda_p
lam2 <- e$lambda_m
lam3 <- e$lambda_s

# array form of regression coefficients for matrix multiplication
coefs1 <- array(c(e$b_p,
                  e$c_p),
                dim = c(Nsim, par$N, 2))

coefs2 <- array(c(e$b_m,
                  e$c_m),
                dim = c(Nsim, par$N, 2))

coefs3 <- array(c(e$b_s,
                  e$c_s),
                dim = c(Nsim, par$N, 2))

# matrix form of the focal site's predictors, needed for matrix multiplication
temp_obs_p1 <- matrix(obs_per[, 1], nrow = Nsim, ncol = par$N, byrow = TRUE)
temp_obs_m1 <- matrix(obs_mea[, 1], nrow = Nsim, ncol = par$N, byrow = TRUE)
temp_obs_s1 <- matrix(obs_sma[, 1], nrow = Nsim, ncol = par$N, byrow = TRUE)

# average probabilities in time
# row = time, columns = mean, 2.5% and 97.5% quantile
pr_time_p <- pr_time_m <- pr_time_s <- matrix(NA, nrow = (par$T - 1), ncol = 3)

# temporary probabilites
pr_temp_p <- array(NA, dim = c(Nsim, par$N, par$T - 1))
pr_temp_m <- array(NA, dim = c(Nsim, par$N, par$T - 1))
pr_temp_s <- array(NA, dim = c(Nsim, par$N, par$T - 1))

# function to get mean and 95 percent (posterior) interval of a vector x
descript <- function(x) {c(mean = mean(x), quantile(x, probs = 0.025), quantile(x, probs = 0.975))}

for (t in 1:(par$T - 1)) {
  # number of missing observations in current time point
  n_missing <- sum(par$missing[, t + 1] == 1)
  
  ## pertussis
  # probabilities for all iterations and sites in current time point
  pr_temp_p[,, t] <- plogis(a1 +
                              s1[, 1 + t %% 12] +
                              lam1 * tau1[, t] + 
                              b1[, 1] * coefs1[, , 1] * temp_obs_p1 + b1[, 2] * coefs1[, , 2] * p_per_temp)
  
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
                              s2[, 1 + t %% 12] +
                              lam2 * tau2[, t] + 
                              b2[, 1] * coefs2[, , 1] * temp_obs_m1 + b2[, 2] * coefs2[, , 2] * p_mea_temp)
  
  pr_time_m[t, ] <- descript(rowMeans(pr_temp_m[,, t]))
  
  temp_obs_m <- matrix(obs_mea[, t + 1], ncol = Nsim, nrow = par$N, byrow = FALSE)
  temp_obs_m[par$missing[, t + 1] == 1, ] <- rbinom(Nsim * n_missing,
                                                              1,
                                                              (pr_temp_m[, par$missing[, t + 1] == 1, t])) > 0
  
  ## smallpox
  pr_temp_s[,, t] <- plogis(a3 +
                              s3[, 1 + t %% 12] +
                              lam3 * tau3[, t] + 
                              b3[, 1] * coefs3[, , 1] * temp_obs_s1 + b3[, 2] * coefs3[, , 2] * p_sma_temp)
  
  pr_time_s[t, ] <- descript(rowMeans(pr_temp_s[,, t]))
  
  temp_obs_s <- matrix(obs_sma[, t + 1], ncol = Nsim, nrow = par$N, byrow = FALSE)
  temp_obs_s[par$missing[, t + 1] == 1, ] <- rbinom(Nsim * n_missing,
                                                              1,
                                                              (pr_temp_s[, par$missing[, t + 1] == 1, t])) > 0
  
  # focal site's past, take the transpose now
  temp_obs_p1 <- t(temp_obs_p)
  temp_obs_m1 <- t(temp_obs_m)
  temp_obs_s1 <- t(temp_obs_s)
  
  # means over neighbors based on the current data, if clause is needed since rowMeans fails if n = 1
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
pr_site_p <- apply(pr_temp_p, c(1, 2), mean) %>% colMeans()
pr_site_m <- apply(pr_temp_m, c(1, 2), mean) %>% colMeans()
pr_site_s <- apply(pr_temp_s, c(1, 2), mean) %>% colMeans()

# save
saveRDS(list(pr_time_p, pr_time_m, pr_time_s,
             pr_site_p, pr_site_m, pr_site_s), "sims_dataa_bern_independent_short.rds")
