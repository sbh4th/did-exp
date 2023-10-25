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
  ggplot(aes(x = estimate, y = year, 
    color = cohort, fill = cohort)) + 
    scale_color_brewer(palette = "Set1") +
    stat_slab(.width = .95) + theme_bw()


fe_model <- bf(y | subset(treat==0) ~ 1 + 
  year_2012 + year_2013 + year_2014 + year_2015 +
    year_2016 + year_2017 + year_2018 + year_2019 +
    year_2020 + cohort_year_2014 + cohort_year_2016)
w_model <- bf(w ~ 1 + e)

b14.6 <-
  brm(data = subdat_sim, 
      family = gaussian,
      e_model + w_model + set_rescor(TRUE),
      prior = c(# E model
                prior(normal(0, 0.2), class = Intercept, resp = e),
                prior(normal(0, 0.5), class = b, resp = e),
                prior(exponential(1), class = sigma, resp = e),
                
                # W model
                prior(normal(0, 0.2), class = Intercept, resp = w),
                prior(normal(0, 0.5), class = b, resp = w),
                prior(exponential(1), class = sigma, resp = w),
                
                # rho
                prior(lkj(2), class = rescor)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      seed = 14,
      file = "fits/b14.06")


b_first <-
  brm(data = subset(db, treat==0), 
      family = gaussian(),
      y ~ 1 + (1 | state) + 
        year_2012 + year_2013 + year_2014 + year_2015 +
        year_2016 + year_2017 + year_2018 + year_2019 +
        year_2020 + cohort_year_2014 + cohort_year_2016,
      prior = c(prior(normal(0, 10), class = Intercept),
        prior(normal(0, 5), class = b),
        prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes", seed = 3958,
      file = "code/fits/brms-re-fs")

res <- residuals(b_first)

newdata <- db

library(tidybayes)

db2 <- add_linpred_draws(object = b_first, 
  newdata = db, draws=1)

bpredict_values <- predict(b_first)

db2 <- NULL
