library(haven)
library(readxl)
library(tidyverse)
library(posterior)
library(cmdstanr)
library(augsynth)

setwd("/home/rgiordan/Documents/git_repos/replication/mtgp_panel/code")
source("fit_models_lib.R")

#Sys.setenv(STAN_NUM_THREADS=4)

start_year <- 1997
end_year <- 2018
t_int <- 2007

crimes <- load_crimes_data(start_year, end_year, t_int)


###################
# Poisson model


# ranks gets passed to the argument n_k_f (by omission)
# Recall that n_k_f is the rank of the time dependence kernel
# fit_pois <-
#   fit_mtgp_pois(
#     n_k_f=5,
#     outcome = homicide,
#     trt = treated,
#     unit = State, time = year, t_int = 2007,
#     data = crimes %>% filter(year >= start_year)
#   )



#model <- cmdstan_model(stan_file = 'stan/single_task_poisson.stan',
#                       include_paths = ".")
model <- cmdstan_model(stan_file = 'stan/poisson.stan',
                       include_paths = ".")

data <- crimes %>% filter(year >= start_year)


foo <- function(time, data) {
  # x is vector of times
  x <- data %>% distinct(!!time) %>% pull(!!time)
  #x <- data %>% select(!!time) %>% distinct() %>% pull(!!time)
  return(x)
}
#debug(foo)

# https://adv-r.hadley.nz/quasiquotation.html
#foo(time="year", data=data)
foo(time=quo(year), data=data)


time <- "year"
x <- data %>% select(!!time) %>% distinct() %>% pull(!!time)
out <- augsynth:::format_data(
  quo(homicide), 
  quo(treated),
  quo(State),
  quo(year),
  t_int, data)

out$X

data$State %>% unique() 

source("fit_models_lib.R")
#debug(make_stan_data)
out <- make_stan_data(
  outcome=quo(homicide),
  trt=quo(treated),
  unit=quo(State),
  time=quo(year), 
  t_int=2007, 
  data=data)


time <- "year"
data %>% pull(!!time)


n_k_f <- 7
samples <- 4000
chains <- 4

standata <- list(x = out$x,
                 y = c(out$y),
                 population = out$pop,
                 N = length(out$x),
                 D = ncol(out$y),
                 n_k_f = n_k_f,
                 control_idx = out$control_idx,
                 num_treated = length(out$y) - length(out$control_idx)
)


options(mc.cores=4)
# sample from model
fit_pois <- model$sample(
  data = standata,
  iter_warmup = samples,
  iter_sampling = samples,
  chains = chains,
  parallel_chains = chains)


# f_draws <-
#   as_draws_df(fit_pois$draws(variables="f")) %>%
#   pivot_longer(cols=c(-.chain, -.iteration, -.draw),
#                names_to="f")
# head(f_draws)


# Tidybayes is a another way to do the posterior and post predictive.
library(tidybayes)
library(GGally)


par_draws <-
  fit_pois %>%
  spread_draws(lengthscale_global, sigma_global, intercept, lp__) 

ggpairs(par_draws, columns=c("lengthscale_global", "sigma_global", "intercept", "lp__"))

par_draws_long <-
  fit_pois %>%
  gather_draws(lengthscale_global, sigma_global, intercept, lp__) 

# Why is lp__ discretized?  Numeric precision problems?
ggplot(par_draws_long) +
  geom_line(aes(x=.iteration, y=.value, color=.chain)) +
  facet_grid(.variable ~ ., scales="free")


fit_pois$draws(variables="lp__")

###########
# f draws

year_index <- data %>%
  distinct(year) %>%
  mutate(year_ind=1:n())

state_index <- data %>%
  distinct(State) %>%
  mutate(state_ind=1:n())

f_draws <-
  fit_pois %>%
  spread_draws(f[year_ind, state_ind]) %>%
  left_join(state_index, "state_ind") %>%
  left_join(year_index, "year_ind") %>%
  ungroup() %>%
  group_by(.chain, .iteration, .draw, State, year) %>%
  mutate(pois_draw=rpois(n=1, lambda=f)) %>%
  ungroup()
head(f_draws)

stopifnot(!any(is.na(f_draws$State)))
stopifnot(!any(is.na(f_draws$year)))

if (FALSE) {
  # Visually check mixing
  ggplot(f_draws %>% filter(State == "IL", year == 1997)) +
    geom_line(aes(x=.draw, y=f, color=.chain))
}


State <- "CA"
plot_data <-
  inner_join(
    filter(data, State == !!State) %>%
      select(year, State, homicide),
    filter(f_draws, State == !!State) %>%
      select(year, State, f, pois_draw),
    by=c("year", "State")
  )
ggplot(plot_data) +
  geom_line(aes(x=year, y=homicide)) +
  geom_point(aes(x=year, y=pois_draw), alpha=0.1) +
  geom_vline(aes(xintercept=2007), color="red")
