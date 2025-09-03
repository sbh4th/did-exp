library(future.apply)
library(progressr)

simulate_once <- function(i, X, country, cohort_wave_map) {
  n_units <- length(unique(country))
  n_obs   <- nrow(X)
  
  # Random effects each iteration
  unit_fe <- rnorm(n_units, country/10, 1)
  wave_fe <- rnorm(4, 0, 0.5)   # varies each simulation
  error   <- rnorm(n_obs, 0, 1)
  
  # Map unit_fe + wave_fe to rows of X
  fe_vec <- unit_fe[rep(1:n_units, each=4)] + 
    wave_fe[rep(1:4, times=n_units)]
  
  # Treatment assignment from cohort map
  treat <- ifelse((X[,"wave"] >= X[,"cohort_wave"]) & 
                    (X[,"cohort_wave"] != 1), 1, 0)
  tau   <- ifelse(treat==1, 0, 0) # you can allow nonzero mu here
  
  tau_cum <- ave(tau, X[,"unit"], FUN=cumsum)
  
  # Outcome
  y <- truncnorm::rtruncnorm(
    n_obs, a=0, b=20,
    mean = 15 - X[,"cohort_wave"] + fe_vec + tau_cum + error,
    sd=4)
  y <- round(y)
  
  # Estimate
  fit <- lm.fit(X, y)
  coefs <- coef(fit)
  
  coefs
}

# ---------- SETUP (once) ----------
n_units <- 57000
waves   <- 1:4

# Cohort assignment (fixed)
treat_taus <- data.frame(
  country = 1:38,
  cohort_wave = sort(c(rep(1, 24), rep(2, 2), 
    rep(3, 4), rep(4, 8)))
)

# Base panel structure
base <- expand.grid(unit=1:n_units, wave=waves)

base <- merge(base, treat_taus, by="country", all.x=TRUE)

# Create model matrix once
X <- model.matrix(
  ~ treat:factor(cohort_wave):factor(wave) +
    factor(cohort_wave) + factor(wave),
  data = base
)

# Keep identifiers outside X
country <- base$country
unit    <- base$unit
cohort_wave_map <- base$cohort_wave

# ---------- SIMULATION ----------
plan(multisession, workers=8)
handlers("txtprogressbar")

n_sims <- 50
results <- with_progress({
  p <- progressor(along=1:n_sims)
  future_lapply(1:n_sims, function(i) {
    out <- simulate_once(i, X, country, cohort_wave_map)
    p()
    out
  }, future.seed=TRUE)
})







n_units <- 1500
n_countries <- 38
waves   <- 4

panel <- expand.grid(
  unit = 1:n_units, 
  country = 1:n_countries, 
  wave = 1:waves)

treat_taus <- data.frame(
    country = 1:n_countries,
    cohort_wave = sort(c(rep(1, 24), rep(2, 2), 
      rep(3, 4), rep(4, 8)))
  )

# unit → country mapping
unit_map <- data.frame(
  unit = rep(1:n_units, n_countries),
  country = sample(1:n_countries, n_units, replace = TRUE)
)

# base panel
base <- expand.grid(unit = 1:nrow(unit_map), 
                    wave = waves)

# add country
base <- merge(base, unit_map, by="unit")

# add cohort_wave
base <- merge(panel, treat_taus, by="country", all.x=TRUE)
