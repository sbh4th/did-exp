library(here)
library(tidyverse)
library(truncnorm)
library(etwfe)
library(marginaleffects)
library(modeldb)

set.seed(81164)
# unit fixed effects (unobserved heterogeneity)
unit <- tibble(
  unit = 1:57000,
  # generate clusters
  country = sample(1:38, 57000, replace = TRUE),
  unit_fe = rnorm(57000, country/10, 1),
  # generate instantaneous treatment effect
  #mu = rnorm(nobs, true_mu, 0.2)
  # mu = 2
)

# wave fixed effects (first part)
wave <- tibble(
  wave = 1:4,
  wave_fe = rnorm(length(wave), 0, 0.5)
)

# Put the clusters into treatment groups
treat_taus <- tibble(
  # sample the clusters randomly
  country = sample(1:38, 38, replace = FALSE),
  # place the randomly sampled countries into 1\{t \ge g \}G_g
  cohort_wave = sort(c(rep(1, 24), rep(2, 2), 
    rep(3, 4), rep(4, 8)))
)

# make main dataset
# full interaction of unit X year 
d <- expand_grid(unit = 1:57000, wave = 1:4) %>% 
  left_join(., unit, by="unit") %>% 
  left_join(., wave, by="wave") %>% 
  left_join(., treat_taus, by="country") %>% 
  # make error term and get treatment indicators and treatment effects
  # Also get cohort specific trends (modify time FE)
  mutate(error = rnorm(57000*4, 0, 1),
    treat = ifelse((wave >= 
      cohort_wave)*(cohort_wave != 1), 1, 0),
      # treatment effect = 1 if 2016, 2 if 2014, annually
      mu = ifelse(cohort_wave==2, 2, 
        ifelse(cohort_wave==3, 2, 2)),
      tau = ifelse(treat == 1, mu, 0),
      wave_fe = wave_fe - 0.5*(wave - cohort_wave)
  ) %>% 
  # calculate cumulative treatment effects
  group_by(unit) %>% 
  mutate(tau_cum = cumsum(tau)) %>% 
  ungroup() %>% 
  # calculate the dependent variable
  mutate(y = rtruncnorm(length(unit),
    a = 0, b = 20,
    mean = 15 - cohort_wave + 
     unit_fe + wave_fe + tau_cum + error,
    sd = 4),
    y = round(y)) %>%
  # Relabel 2018 cohort as never-treated
  mutate(cohort_wave = ifelse(
    cohort_wave == 1, Inf, cohort_wave)) %>%
  # treatment cohort dummies
  add_dummy_variables(cohort_wave, 
    values=c(Inf,2,3,4), remove_original = F) %>%
  # wave dummies
  add_dummy_variables(wave, 
    values=c(1,2,3,4), remove_original = F)
  
theme_set(theme_classic() + 
            theme(plot.background = element_blank()))

plot <- d %>% 
  ggplot(aes(x = wave, y = y, group = unit)) + 
  geom_line(alpha = 1/8, color = "grey") + 
  geom_line(data = d %>% 
    group_by(cohort_wave, wave) %>% 
    summarize(y = mean(y)),
   aes(x = wave, y = y, group = factor(cohort_wave),
    color = factor(cohort_wave)), linewidth = 2) + 
  labs(x = "", y = "Y",  color = "Treatment group   ") +
  scale_x_continuous(breaks = c(1, 2, 3, 4)) +
  geom_vline(xintercept = 2, 
    color = '#E41A1C', linewidth = 2) + 
  geom_vline(xintercept = 3, 
    color = '#377EB8', linewidth = 2) + 
  scale_color_brewer(palette = 'Set1') + 
  theme(legend.position = 'bottom',
        #legend.title = element_blank(), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  scale_color_manual(labels = 
    c("2", "3", "4", "Never-treated"),
  values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984ea3")) +
  ggtitle("Simulated data with heterogeneous treatment effect dynamics across cohorts \n and with a never-treated group")+
  theme(plot.title = element_text(hjust = 0.5, size=12))

plot 

# fit the model
etwfe_model = fixest::feols(
  y ~ treat:i(cohort_wave, i.wave, ref=Inf, ref2 = 1) | 
    cohort_wave + wave,
  data = d,
  vcov = ~country
)

etwfe <- lm(y ~
  treat:cohort_wave_2:wave_2 +
  treat:cohort_wave_2:wave_3 +
  treat:cohort_wave_2:wave_4 +
  treat:cohort_wave_3:wave_3 +
  treat:cohort_wave_3:wave_4 +
  treat:cohort_wave_4:wave_4 +
  wave_2 + wave_3 + wave_4 +
  cohort_wave_2 + cohort_wave_3 +
  cohort_wave_4,
  data = d)

# average marginal predictions
etwfe_mp <- avg_predictions(
  etwfe, 
  newdata = subset(d, treat==1), 
  var = "treat",
  vcov = ~country)

# average marginal effect
etwfe_me <- slopes(etwfe, 
  newdata = subset(d, treat==1), 
  var = "treat", 
  by = "treat",
  vcov = ~country)


# Simulate null just to check
coef_results <- c()
p_results <- c()
sig_results <- c()

for (i in 1:2000) {
# data
  unit <- tibble(
    unit = 1:57000,
    country = sample(1:38, 57000, replace = TRUE),
    unit_fe = rnorm(57000, country/10, 1))
  
  wave <- tibble(
    wave = 1:4,
    wave_fe = rnorm(length(wave), 0, 0.5))
  
  treat_taus <- tibble(
    country = sample(1:38, 38, replace = FALSE),
    cohort_wave = sort(c(rep(1, 24), rep(2, 2), 
                         rep(3, 4), rep(4, 8))))
  
  d <- expand_grid(unit = 1:57000, wave = 1:4) %>% 
    left_join(., unit, by="unit") %>% 
    left_join(., wave, by="wave") %>% 
    left_join(., treat_taus, by="country") %>% 
    mutate(
      error = rnorm(57000*4, 0, 1),
      treat = ifelse((wave >= 
        cohort_wave)*(cohort_wave != 1), 1, 0),
      mu = ifelse(cohort_wave==2, 0, 
        ifelse(cohort_wave==3, 0, 0)),
      tau = ifelse(treat == 1, mu, 0),
      wave_fe = wave_fe - 0.5*(wave - cohort_wave)) %>% 

    group_by(unit) %>% 
    mutate(tau_cum = cumsum(tau)) %>% 
    ungroup() %>% 
    
    mutate(
      y = rtruncnorm(length(unit),
        a = 0, b = 20,
        mean = 15 - cohort_wave + unit_fe + 
        wave_fe + tau_cum + error, sd = 4),
      y = round(y)) %>%

    mutate(cohort_wave = ifelse(
      cohort_wave == 1, Inf, cohort_wave)) %>%
    add_dummy_variables(cohort_wave, 
      values=c(Inf,2,3,4), remove_original = F) %>%
    add_dummy_variables(wave, 
      values=c(1,2,3,4), remove_original = F)

# model and estimates
etwfe <- lm(y ~ treat:cohort_wave_2:wave_2 +
  treat:cohort_wave_2:wave_3 +
  treat:cohort_wave_2:wave_4 +
  treat:cohort_wave_3:wave_3 +
  treat:cohort_wave_3:wave_4 +
  treat:cohort_wave_4:wave_4 +
  wave_2 + wave_3 + wave_4 +
  cohort_wave_2 + cohort_wave_3 + cohort_wave_4,
  data = d)

etwfe_me <- slopes(etwfe, 
  newdata = subset(d, treat==1), 
  var = "treat", 
  by = "treat",
  vcov = ~country)

# results
coef_results[i] <- etwfe_me$estimate
p_results[i] <- etwfe_me$p.value
sig_results[i] <- etwfe_me$p.value <= .05
  
}

# revised version
library(tidyverse)
library(truncnorm)
library(marginaleffects)
library(future.apply)

set.seed(123)

## Precompute static data
N <- 57000
n_countries <- 38
n_waves <- 4
n_sims <- 100

# Units: fixed across sims
units <- tibble(
  unit    = 1:N,
  country = sample(1:n_countries, N, replace = TRUE),
  unit_fe = rnorm(N, mean = (1:n_countries)[country] / 10, sd = 1)
)

# Waves: fixed across sims
waves <- tibble(
  wave    = 1:n_waves,
  wave_fe = rnorm(n_waves, 0, 0.5)
)

## Parallel plan
plan(multisession, workers = parallel::detectCores() - 1)

results <- future_lapply(1:n_sims, function(i) {
  # Random treatment timing assignment
  treat_taus <- tibble(
    country = sample(1:n_countries, n_countries, replace = FALSE),
    cohort_wave = sort(c(rep(1, 24), rep(2, 2), rep(3, 4), rep(4, 8)))
  )
  
  # Expand units × waves (faster than expand_grid + joins)
  d <- tibble(
    unit = rep(units$unit, each = n_waves),
    country = rep(units$country, each = n_waves),
    unit_fe = rep(units$unit_fe, each = n_waves),
    wave = rep(waves$wave, times = N),
    wave_fe = rep(waves$wave_fe, times = N)
  ) %>%
    left_join(treat_taus, by = "country") %>%
    mutate(
      error = rnorm(n(), 0, 1),
      treat = if_else((wave >= cohort_wave) & (cohort_wave != 1), 1, 0),
      mu = 0,
      tau = if_else(treat == 1, mu, 0),
      wave_fe = wave_fe - 0.5 * (wave - cohort_wave)
    ) %>%
    group_by(unit) %>%
    mutate(tau_cum = cumsum(tau)) %>%
    ungroup() %>%
    mutate(
      y = round(rtruncnorm(
        n(),
        a = 0, b = 20,
        mean = 15 - cohort_wave + unit_fe + wave_fe + tau_cum + error,
        sd = 4
      )),
      cohort_wave = if_else(cohort_wave == 1, Inf, cohort_wave)  # Inf = never treated
    )
  
  # Model: let lm() handle dummies with factor()
  etwfe_lm <- lm(
    y ~ treat * factor(cohort_wave) * factor(wave),
    data = d
  )
  
  # Estimate treatment effects
  etwfe_me <- slopes(
    etwfe,
    newdata = subset(d, treat == 1),
    var = "treat",
    by = "treat",
    vcov = ~country
  )
  
  list(
    coef = etwfe_me$estimate,
    pval = etwfe_me$p.value,
    sig  = etwfe_me$p.value <= 0.05
  )
})

## Collect results
coef_results <- sapply(results, `[[`, "coef")
p_results   <- sapply(results, `[[`, "pval")
sig_results <- sapply(results, `[[`, "sig")



### experimental
## Parameters
N <- 57000
n_countries <- 38
n_waves <- 4
n_sims <- 2000

# ---- Precompute static unit + wave information ----
units <- tibble(
  unit    = 1:N,
  country = sample(1:n_countries, N, replace = TRUE),
  unit_fe = rnorm(N, mean = (1:n_countries)[country] / 10, sd = 1)
)

waves <- tibble(
  wave    = 1:n_waves,
  wave_fe = rnorm(n_waves, 0, 0.5)
)

# Expand once: units × waves
base_d <- tibble(
  unit    = rep(units$unit, each = n_waves),
  country = rep(units$country, each = n_waves),
  unit_fe = rep(units$unit_fe, each = n_waves),
  wave    = rep(waves$wave, times = N),
  wave_fe = rep(waves$wave_fe, times = N)
)

# Add static wave dummies (never change)
base_d <- base_d %>%
  mutate(
    wave_2 = as.integer(wave == 2),
    wave_3 = as.integer(wave == 3),
    wave_4 = as.integer(wave == 4)
  )

# ---- Parallel plan ----
plan(multisession, workers = parallel::detectCores() - 1)

results <- future_lapply(1:n_sims, function(i) {
  
  # Treatment timing varies each sim
  treat_taus_d <- tibble(
    country = sample(1:n_countries, n_countries, replace = FALSE),
    cohort_wave = sort(c(rep(1, 24), rep(2, 2), rep(3, 4), rep(4, 8)))
  )
  
  # Merge into base_d
  dr <- base_d %>%
    left_join(treat_taus_d, by = "country") %>%
    mutate(
      error = rnorm(n(), 0, 1),
      treat = if_else((wave >= cohort_wave) & (cohort_wave != 1), 1L, 0L),
      mu = 0,
      tau = if_else(treat == 1, mu, 0),
      wave_fe = wave_fe - 0.5 * (wave - cohort_wave)
    ) %>%
    group_by(unit) %>%
    mutate(tau_cum = cumsum(tau)) %>%
    ungroup() %>%
    mutate(
      y = round(rtruncnorm(
        n(),
        a = 0, b = 20,
        mean = 15 - cohort_wave + unit_fe + wave_fe + tau_cum + error,
        sd = 4
      )),
      cohort_wave = if_else(cohort_wave == 1, Inf, cohort_wave),
      cohort_wave_2 = as.integer(cohort_wave == 2),
      cohort_wave_3 = as.integer(cohort_wave == 3),
      cohort_wave_4 = as.integer(cohort_wave == 4)
    )
  
  # ---- Exact model specification (your original version) ----
  etwfe_dr <- lm(
    y ~ treat:cohort_wave_2:wave_2 +
      treat:cohort_wave_2:wave_3 +
      treat:cohort_wave_2:wave_4 +
      treat:cohort_wave_3:wave_3 +
      treat:cohort_wave_3:wave_4 +
      treat:cohort_wave_4:wave_4 +
      wave_2 + wave_3 + wave_4 +
      cohort_wave_2 + cohort_wave_3 + cohort_wave_4,
    data = dr
  )
  
  
  
  library(tidyverse)
  library(truncnorm)
  library(marginaleffects)
  library(future.apply)
  library(progressr)
  
  set.seed(123)
  
  ## Parameters
  N <- 57000
  n_countries <- 38
  n_waves <- 4
  n_sims <- 100
  
  # ---- Precompute static unit + wave information ----
  units <- tibble(
    unit    = 1:N,
    country = sample(1:n_countries, N, replace = TRUE),
    unit_fe = rnorm(N, mean = (1:n_countries)[country] / 10, sd = 1)
  )
  
  waves <- tibble(
    wave    = 1:n_waves,
    wave_fe = rnorm(n_waves, 0, 0.5)
  )
  
  # Expand once: units × waves
  base_d <- tibble(
    unit    = rep(units$unit, each = n_waves),
    country = rep(units$country, each = n_waves),
    unit_fe = rep(units$unit_fe, each = n_waves),
    wave    = rep(waves$wave, times = N),
    wave_fe = rep(waves$wave_fe, times = N)
  )
  
  # Add static wave dummies (never change)
  base_d <- base_d %>%
    mutate(
      wave_2 = as.integer(wave == 2),
      wave_3 = as.integer(wave == 3),
      wave_4 = as.integer(wave == 4)
    )
  
  # ---- Parallel plan ----
  plan(multisession, workers = parallel::detectCores() - 1)
  
  # ---- Define simulation function ----
  sim_fun <- function(i, p) {
    # Treatment timing varies each sim
    treat_taus <- tibble(
      country = sample(1:n_countries, n_countries, replace = FALSE),
      cohort_wave = sort(c(rep(1, 24), rep(2, 2), rep(3, 4), rep(4, 8)))
    )
    
    dr <- base_d %>%
      left_join(treat_taus, by = "country") %>%
      mutate(
        error = rnorm(n(), 0, 1),
        treat = if_else((wave >= cohort_wave) & (cohort_wave != 1), 1L, 0L),
        mu = 0,
        tau = if_else(treat == 1, mu, 0),
        wave_fe = wave_fe - 0.5 * (wave - cohort_wave)
      ) %>%
      group_by(unit) %>%
      mutate(tau_cum = cumsum(tau)) %>%
      ungroup() %>%
      mutate(
        y = round(rtruncnorm(
          n(),
          a = 0, b = 20,
          mean = 15 - cohort_wave + unit_fe + wave_fe + tau_cum + error,
          sd = 4
        )),
        cohort_wave = if_else(cohort_wave == 1, Inf, cohort_wave),
        cohort_wave_2 = as.integer(cohort_wave == 2),
        cohort_wave_3 = as.integer(cohort_wave == 3),
        cohort_wave_4 = as.integer(cohort_wave == 4)
      )
    
    # Model: your exact specification
    etwfe <- lm(
      y ~ treat:cohort_wave_2:wave_2 +
        treat:cohort_wave_2:wave_3 +
        treat:cohort_wave_2:wave_4 +
        treat:cohort_wave_3:wave_3 +
        treat:cohort_wave_3:wave_4 +
        treat:cohort_wave_4:wave_4 +
        wave_2 + wave_3 + wave_4 +
        cohort_wave_2 + cohort_wave_3 + cohort_wave_4,
      data = dr
    )
    
    etwfe_me <- slopes(
      etwfe,
      newdata = subset(dr, treat == 1),
      var = "treat",
      by = "treat",
      vcov = ~country
    )
    
    p() # progress update
    
    list(
      coef = etwfe_me$estimate,
      pval = etwfe_me$p.value,
      sig  = etwfe_me$p.value <= 0.05
    )
  }
  
  # ---- Run with progress bar ----
  handlers(global = TRUE)   # enable progress bar handlers
  with_progress({
    p <- progressor(steps = n_sims)
    results <- future_lapply(1:n_sims, function(i) sim_fun(i, p),
      future.seed = TRUE)
  })
  
  # ---- Collect results ----
  coef_results <- sapply(results, `[[`, "coef")
  p_results   <- sapply(results, `[[`, "pval")
  sig_results <- sapply(results, `[[`, "sig")
  

  