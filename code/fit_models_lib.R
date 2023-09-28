make_stan_data <- function(outcome, trt, unit, time, t_int, data) {
  
  # x is vector of times
  x <- data %>% distinct(!!time) %>% pull(!!time)
  
  out <- augsynth:::format_data(outcome, trt, unit, time,
                                t_int, data)
  
  pop <- augsynth:::format_data(quo(Population), trt, unit, time,
                                t_int, data)
  # rows of y are time periods, columns units
  y <- t(cbind(out$X, out$y))
  pop <- t(cbind(pop$X, pop$y))
  # get the heldout indices
  excl_idx <- matrix(1:length(y), nrow(y), ncol(y))[x >= t_int, out$trt == 1]
  cntrl_idx <- (1:length(y))[-excl_idx]
  treated_idx = which(out$trt == 1)
  return(list(x = x, y = y, control_idx = cntrl_idx, pop = pop, treated_idx=treated_idx))
}


get_post_pred_samples <- function(post, unit, time, data) {
  
  time <- enquo(time)
  unit <- enquo(unit)
  times <- data %>% distinct(!!time) %>% pull(!!time)
  units <- data %>% distinct(!!unit) %>% pull(!!unit)
  
  # post <- get_posterior(fit_mtgp)
  
  # get time-unit-sample pairs
  time <- as.character(time)[2]
  unit <- as.character(unit)[2]
  times_units_samples <- expand.grid(1:nrow(post), times, units)
  names(times_units_samples) <- c("sample", time, unit)
  times_units_samples$y <- c(post)
  return(times_units_samples)
}

#########################
# Normal fitting

fit_mtgp_normal <- function(outcome, trt, unit, time, t_int, data, n_k_f,
                            iter_warmup = 1000, iter_sampling = 1000,
                            chains = 4, parallel_chains = 4,
                            adapt_delta = 0.9, max_treedepth = 13) {
  # compile
  if(n_k_f > 0) {
    model <- cmdstan_model(stan_file = 'stan/normal.stan',
                           include_paths = ".")
  } else {
    model <- cmdstan_model(stan_file = 'stan/single_task_normal.stan',
                           include_paths = ".")
  }
  # format data
  out <- make_stan_data(enquo(outcome), enquo(trt), enquo(unit),
                        enquo(time), t_int, data)
  out$pop = out$pop / 1e5
  standata <- list(x = out$x,
                   y = out$y / out$pop,
                   inv_population = 1 / out$pop,
                   N = length(out$x),
                   D = ncol(out$y),
                   n_k_f = n_k_f,
                   control_idx = out$control_idx,
                   num_treated = length(out$y) - length(out$control_idx),
                   treatment_idx=out$treated_idx
  )
  # sample from model
  fit <- model$sample(data = standata,
                      iter_warmup = iter_warmup,
                      iter_sampling = iter_sampling,
                      chains = chains,
                      parallel_chains = parallel_chains,
                      adapt_delta = adapt_delta,
                      max_treedepth = max_treedepth)
  return(fit)
}

fit_mtgp_normal_cov <- function(outcome, trt, unit, time, t_int, data, n_k_f, cov,
                                iter_warmup = 1000, iter_sampling = 1000,
                                chains = 4, parallel_chains = 4,
                                adapt_delta = 0.9, max_treedepth = 13) {
  # compile
  model <- cmdstan_model(stan_file = 'stan/normal_covariate_adjustment.stan',
                         include_paths = ".")
  # format data
  out <- make_stan_data(enquo(outcome), enquo(trt), enquo(unit),
                        enquo(time), t_int, data)
  out$pop = out$pop / 1e5
  standata <- list(x = out$x,
                   y = out$y / out$pop,
                   inv_population = 1 / out$pop,
                   N = length(out$x),
                   D = ncol(out$y),
                   n_k_f = n_k_f,
                   control_idx = out$control_idx,
                   num_treated = length(out$y) - length(out$control_idx),
                   treatment_idx=out$treated_idx,
                   num_covariates=ncol(cov),
                   covariates=cov
  )
  # sample from model
  fit <- model$sample(data = standata,
                      iter_warmup = iter_warmup,
                      iter_sampling = iter_sampling,
                      chains = chains,
                      parallel_chains = parallel_chains,
                      adapt_delta = adapt_delta,
                      max_treedepth = max_treedepth)
  return(fit)
}

#########################
# Poisson fitting

fit_mtgp_pois <- function(outcome, trt, unit, time, t_int, data, n_k_f,
                          iter_warmup = 1000, iter_sampling = 1000,
                          chains = 4, parallel_chains = 4,
                          adapt_delta = 0.9, max_treedepth = 13) {
  # compile
  if(n_k_f > 0) {
    model <- cmdstan_model(stan_file = 'stan/poisson.stan',
                           include_paths = ".")
  }else{
    model <- cmdstan_model(stan_file = 'stan/single_task_poisson.stan',
                           include_paths = ".")
  }
  # format data
  out <- make_stan_data(enquo(outcome), enquo(trt), enquo(unit),
                        enquo(time), t_int, data)
  standata <- list(x = out$x,
                   y = c(out$y),
                   population = out$pop,
                   N = length(out$x),
                   D = ncol(out$y),
                   n_k_f = n_k_f,
                   control_idx = out$control_idx,
                   num_treated = length(out$y) - length(out$control_idx)
  )
  # sample from model
  fit <- model$sample(data = standata,
                      iter_warmup = iter_warmup,
                      iter_sampling = iter_sampling,
                      chains = chains,
                      parallel_chains = parallel_chains,
                      adapt_delta = adapt_delta,
                      max_treedepth = max_treedepth)
  return(fit)
}

fit_mtgp_pois_cov <- function(outcome, trt, unit, time, t_int, data, n_k_f, cov,
                              iter_warmup = 1000, iter_sampling = 1000,
                              chains = 4, parallel_chains = 4,
                              adapt_delta = 0.9, max_treedepth = 13) {
  # compile
  if(n_k_f > 0) {
    model <- cmdstan_model(stan_file = 'stan/poisson_covariate_adjustment.stan',
                           include_paths = ".")
  }else{
    model <- cmdstan_model(stan_file = 'stan/single_task_poisson_covariate_adjustment.stan',
                           include_paths = ".")
  }
  # format data
  out <- make_stan_data(enquo(outcome), enquo(trt), enquo(unit),
                        enquo(time), t_int, data)
  standata <- list(x = out$x,
                   y = c(out$y),
                   population = out$pop,
                   N = length(out$x),
                   D = ncol(out$y),
                   n_k_f = n_k_f,
                   control_idx = out$control_idx,
                   num_treated = length(out$y) - length(out$control_idx),
                   num_covariates=ncol(cov),
                   covariates=cov
  )
  # sample from model
  fit <- model$sample(data = standata,
                      iter_warmup = iter_warmup,
                      iter_sampling = iter_sampling,
                      chains = chains,
                      parallel_chains = parallel_chains,
                      adapt_delta = adapt_delta,
                      max_treedepth = max_treedepth)
  return(fit)
}


#################################
# Posterior functions

get_posterior_pois <- function(fit_mtgp) {
  draws <- as_draws_matrix(fit_mtgp$draws())
  means <- subset(draws, variable = "f")
  
  # sample from posterior with rate and intercepts
  posts <- matrix(sapply(c(means), rpois, n = 1), ncol = ncol(means))
  return(posts)
}

get_posterior_homoskedastic <- function(fit_mtgp, pops) {
  draws <- as_draws_matrix(fit_mtgp$draws())
  posts <- subset(draws, variable = "f_samples")
  return(posts)
}

