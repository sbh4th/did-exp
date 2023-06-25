library(cmdstanr)
library(brms)
library(tidybayes)
library(modeldb)
library(ggplot2)

options(mc.cores = 4,
        brms.backend = "cmdstanr")

db <- d %>%
  add_dummy_variables(cohort_year, 
    values=c(-Inf,2014,2016), remove_original = F) %>%
  add_dummy_variables(year, 
    values=c(2011:2020), remove_original = F)

## Model without random effects at state level
b1 <-
  brm(data = db, 
      family = gaussian(),
      y ~ 1 +  treat:cohort_year_2014:year_2014 + 
               treat:cohort_year_2014:year_2015 +
               treat:cohort_year_2014:year_2016 +
               treat:cohort_year_2014:year_2017 +
               treat:cohort_year_2014:year_2018 +
               treat:cohort_year_2014:year_2019 +
               treat:cohort_year_2014:year_2020 +
               treat:cohort_year_2016:year_2016 +
               treat:cohort_year_2016:year_2017 +
               treat:cohort_year_2016:year_2018 +
               treat:cohort_year_2016:year_2019 +
               treat:cohort_year_2016:year_2020 +
               cohort_year_2014 + cohort_year_2016 +
               year_2012 + year_2013 +
               year_2014 + year_2015 +
               year_2016 + year_2017 +
               year_2018 + year_2019 +
               year_2020,
      prior = c(prior(normal(0, 100), class = Intercept),
      prior(normal(0, 100), class = b)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes", seed = 2957,
      file = "code/fits/brms-fe")


## now model including random effects for state

b2 <-
  brm(data = db, 
      family = gaussian(),
      y ~ 1 + (1 | state) + treat:cohort_year_2014:year_2014 + 
               treat:cohort_year_2014:year_2015 +
               treat:cohort_year_2014:year_2016 +
               treat:cohort_year_2014:year_2017 +
               treat:cohort_year_2014:year_2018 +
               treat:cohort_year_2014:year_2019 +
               treat:cohort_year_2014:year_2020 +
               treat:cohort_year_2016:year_2016 +
               treat:cohort_year_2016:year_2017 +
               treat:cohort_year_2016:year_2018 +
               treat:cohort_year_2016:year_2019 +
               treat:cohort_year_2016:year_2020 +
               cohort_year_2014 + cohort_year_2016 +
               year_2012 + year_2013 +
               year_2014 + year_2015 +
               year_2016 + year_2017 +
               year_2018 + year_2019 +
               year_2020,
      prior = c(prior(normal(0, 100), class = Intercept),
        prior(normal(0, 100), class = b),
        prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes", seed = 3056,
      file = "code/fits/brms-re-flat")

# now RE model with more informative priors

b3 <-
  brm(data = db, 
      family = gaussian(),
      y ~ 1 + (1 | state) + treat:cohort_year_2014:year_2014 + 
               treat:cohort_year_2014:year_2015 +
               treat:cohort_year_2014:year_2016 +
               treat:cohort_year_2014:year_2017 +
               treat:cohort_year_2014:year_2018 +
               treat:cohort_year_2014:year_2019 +
               treat:cohort_year_2014:year_2020 +
               treat:cohort_year_2016:year_2016 +
               treat:cohort_year_2016:year_2017 +
               treat:cohort_year_2016:year_2018 +
               treat:cohort_year_2016:year_2019 +
               treat:cohort_year_2016:year_2020 +
               cohort_year_2014 + cohort_year_2016 +
               year_2012 + year_2013 +
               year_2014 + year_2015 +
               year_2016 + year_2017 +
               year_2018 + year_2019 +
               year_2020,
      prior = c(prior(normal(0, 10), class = Intercept),
        prior(normal(0, 5), class = b),
        prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes", seed = 49960,
      file = "code/fits/brms-re-inf")

# compare models using leave-one-out criterion
b1 <- add_criterion(b1, "loo")
b2 <- add_criterion(b2, "loo")
b3 <- add_criterion(b3, "loo")

l <- loo_compare(b1, b2, b3, criterion = "loo")
print(l, simplify = F)

# visualize the priors

b3_p <- prior_draws(b3) 

b3_p %>% ggplot(aes(x = b)) + 
  stat_halfeye(color= "black", fill =  '#1b9e77') +
  labs(x = "Parameter value", 
    y = "Density") +
  theme_classic()

# check the chains

b3 %>% 
  plot(variable = "b_Intercept", regex = T) +
  theme(axis.text.y = element_text(hjust = 0))

# posterior predictive checks
pp_check(b3, ndraws=1000)

# posterior distribution (needs revision)

pd <- as_draws_df(b3)
 
pd %>% select(ends_with("treat")) %>%
  pivot_longer(cols = everything(), 
    names_to = c("cohort", "year"), 
    names_pattern = "b_cohort_year_(\\d+):year_(\\d+):treat", 
    values_to = "estimate") %>%
  ggplot(aes(x = year, y = estimate)) +
    stat_halfeye(.width = c(.95, .5)) + facet_wrap(~cohort)





bme_c <- slopes(
  b5, 
  newdata   = subset(db, treat & cohort_year),
  variables = "treat", 
  by        = "cohort_year"
  )

bme_t <- slopes(
  b5, 
  newdata   = subset(db, treat & year>=2014),
  variables = "treat", 
  by        = "year"
  )

d <- d %>%
  mutate(event = year - cohort_year)

# Group effects for ETWFE (by time since treatment)
me_e <- slopes(
  etwfe, 
  newdata   = subset(d, treat & event>=0),
  variables = "treat", 
  by        = "event"
  )

# marginal effects 

bme_avg <- slopes(
  b3, 
  newdata   = subset(db, treat==1),
  variables = "treat", 
  by        = "treat"
  )

# wrangle aggregate ATTs for model summary table
bti <- data.frame(
  term = paste("ATT(", bme_avg$term, ")", sep = ""),
  estimate = bme_avg$estimate,
  conf.low = bme_avg$conf.low,
  conf.high = bme_avg$conf.high,
  std.error = abs(bme_avg$conf.high - 
                    bme_avg$conf.low) / (2 * 1.96)
)

gl <- data.frame()

betwfe_me_avg <- list(tidy = bti, glance = gl)
class(betwfe_me_avg) <- "modelsummary_list"

modelsummary(list("ETWFE Simple Average" = etwfe_me_avg, 
                  "Bayesian Simple Average" = betwfe_me_avg),
  shape = term ~ model + statistic, 
  statistic = c("({std.error})", "conf.int"),
  gof_omit ='._*')


bme_c <- slopes(
  b3, 
  newdata   = subset(db, treat & cohort_year),
  variables = "treat", 
  by        = "cohort_year"
  )

# wrangle aggregate ATTs for model summary table
bti_c <- data.frame(
  term = paste("ATT(g", bme_c$cohort_year, ")", sep=""), 
  estimate = bme_c$estimate,
  conf.low = bme_c$conf.low,
  conf.high = bme_c$conf.high,
  std.error = abs(bme_c$conf.high - 
                    bme_c$conf.low) / (2 * 1.96)
)

betwfe_me_c <- list(tidy = bti_c, glance = gl)
class(betwfe_me_c) <- "modelsummary_list"

modelsummary(list("ETWFE Cohort Average" = etwfe_me_c, 
                  "Bayesian Cohort Average" = betwfe_me_c),
  shape = term ~ model + statistic, 
  statistic = c("({std.error})", "conf.int"),
  gof_omit ='._*')


bme_t <- slopes(
  b3, 
  newdata   = subset(db, treat & year>=2014),
  variables = "treat", 
  by        = "year"
  )

# wrangle aggregate ATTs for model summary table
bti_t <- data.frame(
  term = paste("ATT(", bme_t$year, ")", sep=""), 
  estimate = bme_t$estimate,
  conf.low = bme_t$conf.low,
  conf.high = bme_t$conf.high,
  std.error = abs(bme_t$conf.high - 
                    bme_t$conf.low) / (2 * 1.96)
)

betwfe_me_t <- list(tidy = bti_t, glance = gl)
class(betwfe_me_t) <- "modelsummary_list"

modelsummary(list("ETWFE Cohort Average" = etwfe_me_t, 
                  "Bayesian Cohort Average" = betwfe_me_t),
  shape = term ~ model + statistic, 
  statistic = c("({std.error})", "conf.int"),
  gof_omit ='._*')


# rename coefficients from the models for table
cm5 <- c("b_cohort_year_2014:year_2014" = "ATT(2014,2014)",
      "b_treat:cohort_year_2014:year_2015" = "ATT(2014,2015)",
      "b_treat:cohort_year_2014:year_2016" = "ATT(2014,2016)",
      "b_treat:cohort_year_2014:year_2017" = "ATT(2014,2017)",
      "b_treat:cohort_year_2014:year_2018" = "ATT(2014,2018)",
      "b_treat:cohort_year_2014:year_2019" = "ATT(2014,2019)",
      "b_treat:cohort_year_2014:year_2020" = "ATT(2014,2020)",
      "b_treat:year_2016:cohort_year_2016" = "ATT(2016,2016)",
      "b_treat:year_2017:cohort_year_2016" = "ATT(2016,2017)",
      "b_treat:year_2018:cohort_year_2016" = "ATT(2016,2018)",
      "b_treat:year_2019:cohort_year_2016" = "ATT(2016,2019)",
      "b_treat:year_2020:cohort_year_2016" = "ATT(2016,2020)")

# build the table
modelsummary(b5, 
  coef_map = cm5, statistic="mad",
  gof_omit ='._*')


cm <- 
  c("b_cohort_year_2014:year_2014:treat" = "ATT(2014,2014)",
    "b_cohort_year_2014:year_2015:treat" = "ATT(2014,2015)",
    "b_cohort_year_2014:year_2016:treat" = "ATT(2014,2016)",
    "b_cohort_year_2014:year_2017:treat" = "ATT(2014,2017)",
    "b_cohort_year_2014:year_2018:treat" = "ATT(2014,2018)",
    "b_cohort_year_2014:year_2019:treat" = "ATT(2014,2019)",
    "b_cohort_year_2014:year_2020:treat" = "ATT(2014,2020)",
    "b_year_2016:cohort_year_2016:treat" = "ATT(2016,2016)",
    "b_year_2017:cohort_year_2016:treat" = "ATT(2016,2017)",
    "b_year_2018:cohort_year_2016:treat" = "ATT(2016,2018)",
    "b_year_2019:cohort_year_2016:treat" = "ATT(2016,2019)",
    "b_year_2020:cohort_year_2016:treat" = "ATT(2016,2020)")

# build the table
modelsummary(list("Fixed only" = b1, 
 "State RE (non-informative)" = b2, 
 "State RE (informative)" = b3), 
  coef_map = cm, statistic="mad",
  shape = term ~ model,
  gof_omit ='._*')



