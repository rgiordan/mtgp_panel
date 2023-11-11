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




pops <- make_stan_data(
  quo(homicide), quo(treated), quo(State), quo(year), 2007,
  crimes %>% filter(year >= start_year))$pop
posts_pois <- get_posterior_pois(fit_pois)
post_pred_pois <- get_post_pred_samples(
  posts_pois,
  unit=State, 
  time=year,
  data = crimes %>% filter(year >= start_year))

saveRDS(fit_ranks_pois, "sandbox_mtgp_fits_poisson.RDS")


