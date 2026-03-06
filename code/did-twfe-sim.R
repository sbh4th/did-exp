library(dplyr)
library(truncnorm)
library(progressr)
library(future.apply)
library(fastDummies)
library(fixest)

# --- Simulation parameters ---
n_units     <- 1500    # units per country
n_countries <- 38
n_waves     <- 4
n_sims      <- 50

# --- Base panel construction ---

# Full unit × wave panel
base <- expand.grid(
  unit    = 1:n_units, 
  country = 1:n_countries, 
  wave    = 1:n_waves)
  
# Cohort assignment: 24 never-treated, 
# 2 cohort 2, 4 cohort 3, 8 cohort 4
treat_taus <- data.frame(
  country = 1:n_countries,
  cohort_wave = sort(c(rep(1, 24), rep(2, 2), 
    rep(3, 4), rep(4, 8)))
)

# Full unit × wave panel
base <- expand.grid(
  unit    = 1:(n_units * n_countries),
  wave    = waves
) %>%
  mutate(country = (unit - 1) %/% n_units + 1) %>%  # assign country by blocks
  left_join(treat_taus, by = "country") %>%
  mutate(cohort_wave = ifelse(cohort_wave == 1, Inf, cohort_wave))

# Add dummy variables for cohort_wave & wave (fixed once)
base <- base %>%
  fastDummies::dummy_cols(select_columns = "cohort_wave",
                          remove_selected_columns = FALSE) %>%
  fastDummies::dummy_cols(select_columns = "wave",
                          remove_selected_columns = FALSE)

# --- Simulation function ---
sim_fun <- function(iter, p) {
  n_obs <- nrow(base)
  
  # unit FE: one draw per unit
  unit_fe <- rnorm(n_units * n_countries, 0, 1)
  unit_fe <- unit_fe[base$unit]
  
  # wave FE: one draw per wave
  wave_fe_lookup <- rnorm(length(waves), 0, 0.5)
  wave_fe <- wave_fe_lookup[base$wave]
  
  # idiosyncratic error
  error <- rnorm(n_obs, 0, 1)
  
  # treatment path
  d <- base %>%
    mutate(
      treat = ifelse((wave >= cohort_wave) & 
        (cohort_wave != Inf), 1, 0),
      tau = 0,
      tau_cum = ave(tau, unit, FUN = cumsum),
      y = round(rtruncnorm(
        n_obs, a = 0, b = 20,
        mean = 15 - cohort_wave + unit_fe + 
          wave_fe + tau_cum + error,
        sd = 4
      ))
    )
  
  twfe <- fixest::feols(y ~ treat | country + wave, 
    data = d, vcov = ~country)
  me <- tidy(twfe)
  est <- me$estimate[1]
  est_sig <- me$p.value[1] < 0.05
  
  p() # progress update
  return(list(estimate = est,
              significant = est_sig))
}

# --- Run simulations ---
plan(multisession, workers = 4)
handlers(global = TRUE)

with_progress({
  p <- progressor(steps = n_sims)
  results <- future_lapply(1:n_sims, sim_fun, p = p, future.seed = TRUE)
})



# build balanced panel
base <- expand.grid(unit = 1:(n_units*n_countries), wave = waves)
base$country <- rep(1:n_countries, each = n_units * length(waves))

# country fixed effects: one draw per country
country_fe <- rnorm(n_countries, mean = (1:n_countries)/10, sd = 1)
base$country_fe <- country_fe[base$country]

# wave fixed effects: one draw per wave
wave_fe <- rnorm(length(waves), mean = 0, sd = 0.5)
base$wave_fe <- wave_fe[base$wave]

