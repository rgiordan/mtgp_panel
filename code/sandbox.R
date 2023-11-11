library(haven)
library(readxl)
library(tidyverse)
library(posterior)
library(cmdstanr)
library(augsynth)

setwd("/home/rgiordan/Documents/git_repos/replication/mtgp_panel/code")
source("fit_models_lib.R")

start_year <- 1997
end_year <- 2018
t_int <- 2007


read_dta("../data/ucrthrough2018.dta") %>%
  mutate(across(c(violent_crime, homicide, rape_legacy, rape_revised, robbery,
                  aggravated_assault, property_crime, burglary, larceny,
                  motor_vehicle_theft),
                as.integer)) %>%
  rename(State = state_abbr, rape_rate = rape_legacy_rate, 
         murder_rate = homicide_rate, violent_rate = violent_crime_rate,
         assault_rate = aggravated_assault_rate,
         property_rate = property_crime_rate,
         mvt_rate = motor_vehicle_theft_rate, 
         Population = population) %>%
  filter(State != "DC") %>%
  # drop rape rate after 2015 since it switches to new definition
  mutate(rape_rate = ifelse(year > 2016, NA, rape_rate),
         treated = State == "CA", trt = treated * (year >= 2007)) -> crimes

if(start_year < 1995) {
  crimes <- crimes %>% filter(State != "MS")
}


###################
# Poisson model


# ranks gets passed to the argument n_k_f (by omission)
# Recall that n_k_f is the rank of the time dependence kernel
fit_pois <-
  fit_mtgp_pois(
    n_k_f=5,
    outcome = homicide,
    trt = treated,
    unit = State, time = year, t_int = 2007,
    data = crimes %>% filter(year >= start_year)
  )



model <- cmdstan_model(stan_file = 'stan/single_task_poisson.stan',
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
samples <- 500
chains <- 1

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
fit_pois <- model$sample(data = standata,
                    iter_warmup = samples,
                    iter_sampling = samples,
                    chains = chains)


f_draws <-
  as_draws_df(fit_pois$draws(variables="f")) %>%
  pivot_longer(cols=c(-.chain, -.iteration, -.draw),
               names_to="f")
head(f_draws)


# Tidybayes is a better way to do this.
library(tidybayes)

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
